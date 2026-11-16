import sys

sys.dont_write_bytecode = True
import argparse
import csv
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from configuration import CODE_ROOT, find_rscript, resolved_workspace

PYTHON_VERSION = "3.12.10"
PYTHON_PACKAGES = {
    "matplotlib": "3.11.2",
    "numpy": "2.4.6",
    "openpyxl": "3.1.5",
    "pandas": "3.0.3",
    "Pillow": "13.0.0",
    "pypdf": "6.13.3",
    "pypdfium2": "5.13.0",
    "reportlab": "4.4.9",
    "scipy": "1.18.1",
    "statsmodels": "0.15.0",
}
PYTHON_IMPORTS = {"Pillow": "PIL"}


def python_binary(root):
    return root / ("Scripts/python.exe" if os.name == "nt" else "bin/python")


def execute(command, **kwargs):
    environment = os.environ.copy()
    environment.update({"PYTHONDONTWRITEBYTECODE": "1", "PYTHONUTF8": "1"})
    return subprocess.run(
        [str(x) for x in command],
        env=environment,
        text=True,
        encoding="utf-8",
        errors="replace",
        **kwargs
    )


def python_check(python):
    program = (
        "import sys,json,importlib,importlib.metadata as m; names="
        + repr({name: PYTHON_IMPORTS.get(name, name) for name in PYTHON_PACKAGES})
        + "; out={'python':sys.version.split()[0],'packages':{}}\nfor name,module in names.items():\n try:\n  version=m.version(name);importlib.import_module(module);out['packages'][name]={'version':version,'loadable':True}\n except Exception as error:out['packages'][name]={'version':None,'loadable':False,'error':str(error)}\nprint(json.dumps(out))"
    )
    result = execute([python, "-c", program], capture_output=True, timeout=120)
    if result.returncode:
        return {"python": None, "error": result.stderr, "packages": {}, "passed": False}
    report = json.loads(result.stdout)
    for name, expected in PYTHON_PACKAGES.items():
        record = report["packages"][name]
        record["expected"] = expected
        record["passed"] = record["loadable"] and record["version"] == expected
    report["expected_python"] = PYTHON_VERSION
    report["passed"] = report["python"] == PYTHON_VERSION and all(
        (row["passed"] for row in report["packages"].values())
    )
    return report


def tool_check():
    executable = shutil.which("pdftoppm")
    poppler = {"path": executable, "passed": False}
    if executable:
        result = execute([executable, "-v"], capture_output=True, timeout=15)
        poppler.update(
            {"version": (result.stdout + result.stderr).strip(), "passed": result.returncode == 0}
        )
    font_root = Path("C:/Windows/Fonts")
    fonts = [
        {"path": str(font_root / name), "passed": (font_root / name).is_file()}
        for name in ["arial.ttf", "arialbd.ttf"]
    ]
    return {
        "pdftoppm": poppler,
        "Arial": fonts,
        "passed": poppler["passed"] and all((row["passed"] for row in fonts)),
    }


def r_check(rscript, libraries, directory):
    report_file = directory / "r-package-check.tsv"
    command = [
        rscript,
        "--vanilla",
        CODE_ROOT / "dependencies.R",
        "--check",
        "--report=" + str(report_file),
    ]
    command.extend(("--library=" + str(path) for path in libraries))
    result = execute(command, capture_output=True, timeout=180)
    rows = []
    if report_file.exists():
        with report_file.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
    build_tools = next((row for row in rows if row["package"] == "R_build_tools"), None)
    return {
        "passed": result.returncode == 0,
        "packages": [row for row in rows if row["package"] != "R_build_tools"],
        "build_tools": build_tools,
        "stdout": result.stdout,
        "stderr": result.stderr,
    }


def install_python(python, workspace):
    destination = workspace / "environment/python"
    target = python_binary(destination)
    version = execute(
        [python, "-c", "import sys; print(sys.version.split()[0])"], capture_output=True, check=True
    ).stdout.strip()
    if version != PYTHON_VERSION:
        raise RuntimeError(
            "Use Python " + PYTHON_VERSION + " to create this environment; found " + version
        )
    if not target.exists():
        execute([python, "-m", "venv", destination], check=True)
    requirements = workspace / "environment/python-requirements.txt"
    requirements.write_text(
        "\n".join((name + "==" + version for name, version in PYTHON_PACKAGES.items())) + "\n",
        encoding="utf-8",
    )
    execute([target, "-m", "pip", "install", "--requirement", requirements], check=True)
    frozen = execute(
        [target, "-m", "pip", "freeze", "--all"], capture_output=True, check=True
    ).stdout
    (workspace / "environment/python-lock.txt").write_text(frozen, encoding="utf-8")
    return target


def install_r(rscript, workspace):
    library = workspace / "environment/R_library"
    library.mkdir(parents=True, exist_ok=True)
    execute(
        [
            rscript,
            "--vanilla",
            CODE_ROOT / "dependencies.R",
            "--install",
            "--library=" + str(library),
            "--report=" + str(workspace / "environment/r-install-check.tsv"),
            "--lock=" + str(workspace / "environment/r-lock.tsv"),
        ],
        check=True,
    )
    return library


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", nargs="?", default="check", choices=["check", "install"])
    parser.add_argument("--workspace")
    parser.add_argument("--report-dir")
    parser.add_argument("--python")
    parser.add_argument("--rscript")
    parser.add_argument("--r-library", action="append", default=[])
    parser.add_argument("--component", choices=["all", "python", "r"], default="all")
    args = parser.parse_args()
    workspace = resolved_workspace(args.workspace) if args.workspace else None
    if args.action == "install" and workspace is None:
        parser.error(
            "install requires --workspace outside the code and original project directories"
        )
    python = Path(args.python).resolve() if args.python else Path(sys.executable)
    if workspace and (not args.python) and python_binary(workspace / "environment/python").exists():
        python = python_binary(workspace / "environment/python")
    libraries = [Path(path).expanduser().resolve() for path in args.r_library]
    if workspace and (workspace / "environment/R_library").is_dir():
        libraries.insert(0, workspace / "environment/R_library")
    rscript = find_rscript(args.rscript) if args.component in ["all", "r"] else None
    if args.action == "install":
        (workspace / "environment").mkdir(parents=True, exist_ok=True)
        if args.component in ["all", "python"]:
            python = install_python(python, workspace)
        if args.component in ["all", "r"]:
            libraries = [install_r(rscript, workspace)]
    report = {
        "action": args.action,
        "python_executable": str(python),
        "rscript": rscript,
        "r_libraries": [str(path) for path in libraries],
        "tools": tool_check(),
    }
    with tempfile.TemporaryDirectory(prefix="sle_dependency_check_") as temporary:
        if args.component in ["all", "python"]:
            report["python"] = python_check(python)
        if args.component in ["all", "r"]:
            report["r"] = r_check(rscript, libraries, Path(temporary))
    required = [report[key] for key in ["python", "r"] if key in report]
    if args.component == "all":
        required.append(report["tools"])
    report["passed"] = all((section["passed"] for section in required))
    destination = (
        Path(args.report_dir).expanduser().resolve()
        if args.report_dir
        else workspace / "environment/reports" if workspace else None
    )
    if destination:
        if destination == CODE_ROOT or destination.is_relative_to(CODE_ROOT):
            raise ValueError("Dependency reports must be outside the code repository")
        destination.mkdir(parents=True, exist_ok=True)
        (destination / "dependency-check.json").write_text(
            json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
        )
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
