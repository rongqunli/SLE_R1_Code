import sys

sys.dont_write_bytecode = True
import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
from configuration import (
    CODE_ROOT,
    CODE_SUFFIXES,
    SOURCE_DEFAULTS,
    ORIGINAL_PARENT,
    resolved_workspace,
    source_paths,
    find_rscript,
    environment,
    expanded,
    values,
)
from catalog import catalog, select_stages, select_inputs


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1048576), b""):
            value.update(block)
    return value.hexdigest()


def repository_digest():
    value = hashlib.sha256()
    for path in sorted(CODE_ROOT.rglob("*")):
        if path.is_file() and path.suffix in CODE_SUFFIXES:
            value.update(path.relative_to(CODE_ROOT).as_posix().encode())
            value.update(path.read_bytes())
    return value.hexdigest()


def local_script(stage):
    rel = Path(stage["script"])
    if rel.is_absolute() or ".." in rel.parts:
        raise ValueError("Stage scripts must be safe relative paths")
    return CODE_ROOT / rel if rel.parts[0] == "steps" else CODE_ROOT / "src" / rel


def runtime_script(stage, root):
    rel = Path(stage["script"])
    return root / "_code" / rel if rel.parts[0] == "steps" else root / rel


def render_code(text, root):
    originals = [
        ORIGINAL_PARENT.replace("/", "\\\\"),
        ORIGINAL_PARENT.replace("/", "\\"),
        ORIGINAL_PARENT,
    ]
    pattern = (
        "(?P<prefix>"
        + "|".join((re.escape(value) for value in originals))
        + ")(?=[/\\\\\\\"'\\r\\n]|$)(?P<suffix>[^\\\"'\\r\\n]*)"
    )
    return re.sub(
        pattern, lambda match: root.as_posix() + re.sub("\\\\+", "/", match.group("suffix")), text
    )


def materialize_code(root):
    for source_dir, target_dir in [
        (CODE_ROOT / "src", root),
        (CODE_ROOT / "steps", root / "_code" / "steps"),
    ]:
        for src in source_dir.rglob("*"):
            if src.is_file() and src.suffix in CODE_SUFFIXES:
                target = target_dir / src.relative_to(source_dir)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(
                    render_code(src.read_text(encoding="utf-8-sig"), root),
                    encoding="utf-8",
                    newline="\n",
                )
    (root / "environment" / "R_library").mkdir(parents=True, exist_ok=True)


def validate_inputs(inputs, sources, root):
    results = []
    for item in inputs:
        if item["source"] not in sources:
            raise ValueError("Unknown source family: " + item["source"])
        src = (sources[item["source"]] / item["path"]).resolve()
        target = (root / item["target"]).resolve()
        if not target.is_relative_to(root):
            raise ValueError("Input target escapes the run workspace")
        present = src.is_dir() if item.get("kind") == "dir" else src.is_file()
        results.append(
            {
                "source": str(src),
                "target": str(target),
                "available": present,
                "optional": bool(item.get("optional")),
                "url": item.get("url"),
                "stages": item.get("stages", []),
            }
        )
    return results


def prepare_inputs(inputs, sources, root):
    statuses = validate_inputs(inputs, sources, root)
    missing = [
        item
        for item in statuses
        if not item["available"] and (not item["optional"]) and (not Path(item["target"]).exists())
    ]
    if missing:
        raise FileNotFoundError(
            "Missing external input resources:\n" + "\n".join((item["source"] for item in missing))
        )
    for item in statuses:
        src, dst = (Path(item["source"]), Path(item["target"]))
        if not item["available"]:
            item["sha256"] = digest(dst) if dst.is_file() else None
            continue
        if src.is_dir():
            raise ValueError(
                "Directory inputs must be expanded into explicit files before release: " + str(src)
            )
        dst.parent.mkdir(parents=True, exist_ok=True)
        source_hash = digest(src)
        if dst.exists():
            if (
                not dst.is_file()
                or src.stat().st_size != dst.stat().st_size
                or source_hash != digest(dst)
            ):
                raise ValueError(
                    "External input content differs from its existing runtime copy: "
                    + str(dst)
                    + "; use a new workspace for changed inputs"
                )
        else:
            shutil.copy2(src, dst)
            if source_hash != digest(dst):
                raise ValueError("External input changed while being copied: " + str(src))
        item["sha256"] = source_hash
    state = root / "_execution"
    state.mkdir(exist_ok=True)
    (state / "inputs.json").write_text(json.dumps(statuses, indent=2), encoding="utf-8")
    return statuses


def inputs_digest(statuses, stage_id):
    records = [
        {key: item.get(key) for key in ("source", "target", "sha256")}
        for item in statuses
        if stage_id in item.get("stages", [])
    ]
    records.sort(key=lambda item: item["target"])
    return hashlib.sha256(
        json.dumps(records, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def outputs_exist(stage, root):
    paths = stage.get("outputs", [])
    return bool(paths) and all(
        (
            (
                (root / value).exists()
                if not any((c in value for c in "*?["))
                else any(root.glob(value))
            )
            for value in paths
        )
    )


def execute(stages, inputs, args, root, sources):
    code_hash = repository_digest()
    inputs = select_inputs(inputs, stages)
    materialize_code(root)
    statuses = prepare_inputs(inputs, sources, root)
    directory = root / "_execution"
    completed = {}
    for stage in stages:
        dependencies = {}
        for parent in stage.get("depends", []):
            if parent not in completed:
                raise ValueError(
                    "Selected stages are not a complete dependency order: " + stage["id"]
                )
            dependencies[parent] = completed[parent]["run_id"]
        input_hash = inputs_digest(statuses, stage["id"])
        previous = directory / (stage["id"] + ".json")
        if args.resume and previous.exists():
            old = json.loads(previous.read_text(encoding="utf-8"))
            if (
                old.get("status") == "complete"
                and old.get("run_id")
                and (old.get("code_sha256") == code_hash)
                and (old.get("inputs_sha256") == input_hash)
                and (old.get("dependencies") == dependencies)
                and outputs_exist(stage, root)
            ):
                print("REUSE " + stage["id"], flush=True)
                completed[stage["id"]] = old
                continue
        script = runtime_script(stage, root)
        env = environment(root, sources, stage.get("env"))
        env["SLE_PYTHON"] = str(Path(args.python).resolve())
        env["SLE_RSCRIPT"] = find_rscript(args.rscript)
        libraries = getattr(args, "r_library", None) or []
        if libraries:
            paths = [str(Path(value).expanduser().resolve()) for value in libraries]
            if any((not Path(value).is_dir() for value in paths)):
                raise FileNotFoundError("An explicitly selected R package library does not exist")
            env["R_LIBS"] = os.pathsep.join(paths + ([env["R_LIBS"]] if env.get("R_LIBS") else []))
        replacements = values(root, sources)
        engine = stage["engine"].lower()
        if engine not in {"python", "r"}:
            raise ValueError("Unsupported stage engine: " + engine)
        executable = args.python if engine == "python" else find_rscript(args.rscript)
        command = (
            [executable]
            + (["--vanilla"] if engine == "r" else [])
            + [str(script)]
            + [expanded(value, replacements) for value in stage.get("args", [])]
        )
        cwd = Path(expanded(stage.get("cwd", str(script.parent)), replacements))
        if not cwd.is_absolute():
            cwd = root / cwd
        cwd.mkdir(parents=True, exist_ok=True)
        record = {
            "stage": stage["id"],
            "run_id": uuid.uuid4().hex,
            "status": "running",
            "command": command,
            "cwd": str(cwd),
            "code_sha256": code_hash,
            "inputs_sha256": input_hash,
            "dependencies": dependencies,
            "started": datetime.now(timezone.utc).isoformat(),
        }
        previous.write_text(json.dumps(record, indent=2), encoding="utf-8")
        print("RUN " + stage["id"], flush=True)
        with (directory / (stage["id"] + ".log")).open("wb") as stream:
            process = subprocess.run(
                command, cwd=cwd, env=env, stdout=stream, stderr=subprocess.STDOUT
            )
        record["returncode"] = process.returncode
        record["finished"] = datetime.now(timezone.utc).isoformat()
        record["status"] = (
            "complete" if process.returncode == 0 and outputs_exist(stage, root) else "failed"
        )
        previous.write_text(json.dumps(record, indent=2), encoding="utf-8")
        if record["status"] != "complete":
            raise RuntimeError(
                "Stage failed or did not produce its declared outputs: "
                + stage["id"]
                + "; see "
                + str(directory / (stage["id"] + ".log"))
            )
        completed[stage["id"]] = record


def main():
    parser = argparse.ArgumentParser(
        description="Run the manuscript analyses in a separate workspace. Use dependencies.py to check or install the frozen environment."
    )
    parser.add_argument(
        "command", choices=["check", "plan", "prepare", "run"], nargs="?", default="check"
    )
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--sle", default=str(SOURCE_DEFAULTS["SLE"]))
    parser.add_argument("--sle-v1", default=str(SOURCE_DEFAULTS["SLE_v1"]))
    parser.add_argument("--sle-r1", default=str(SOURCE_DEFAULTS["SLE_R1"]))
    parser.add_argument("--rscript")
    parser.add_argument(
        "--r-library",
        action="append",
        help="Additional existing R package library; repeat for multiple libraries",
    )
    parser.add_argument(
        "--python", help="Python executable; defaults to the workspace environment when installed"
    )
    parser.add_argument("--stage", action="append")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    root, sources = (resolved_workspace(args.workspace), source_paths(args))
    installed_python = (
        root
        / "environment"
        / "python"
        / ("Scripts/python.exe" if os.name == "nt" else "bin/python")
    )
    args.python = (
        args.python
        or os.environ.get("SLE_PYTHON")
        or (str(installed_python) if installed_python.is_file() else sys.executable)
    )
    if not Path(args.python).is_file():
        raise FileNotFoundError("Python executable was not found: " + args.python)
    for source in sources.values():
        if root == source or root.is_relative_to(source) or source.is_relative_to(root):
            raise ValueError("The run workspace overlaps an input source")
    stages, inputs = catalog()
    if not stages:
        raise ValueError("No analysis stages were registered")
    selected = select_stages(stages, args.stage)
    inputs = select_inputs(inputs, selected)
    for stage in stages:
        if not local_script(stage).is_file():
            raise FileNotFoundError(local_script(stage))
    if args.command == "plan":
        print(json.dumps(selected, ensure_ascii=False, indent=2))
    elif args.command == "check":
        status = validate_inputs(inputs, sources, root)
        missing = [
            item
            for item in status
            if not item["available"]
            and (not item["optional"])
            and (not Path(item["target"]).exists())
        ]
        print(
            json.dumps(
                {
                    "stages": len(stages),
                    "selected": [s["id"] for s in selected],
                    "external_inputs": len(inputs),
                    "missing": missing,
                    "rscript": find_rscript(args.rscript),
                    "workspace": str(root),
                },
                ensure_ascii=False,
                indent=2,
            )
        )
        if missing:
            raise SystemExit(2)
    elif args.command == "prepare":
        materialize_code(root)
        prepare_inputs(inputs, sources, root)
    else:
        execute(selected, inputs, args, root, sources)


if __name__ == "__main__":
    main()
