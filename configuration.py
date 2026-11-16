import sys

sys.dont_write_bytecode = True
import os
import shutil
from pathlib import Path

CODE_ROOT = Path(__file__).resolve().parent
SOURCE_DEFAULTS = {
    "SLE": Path("F:/program/Hei/SLE"),
    "SLE_v1": Path("F:/program/Hei/SLE_v1"),
    "SLE_R1": Path("F:/program/Hei/SLE_R1"),
}
ORIGINAL_PARENT = "F:/program/Hei"
CODE_SUFFIXES = {".py", ".R", ".r", ".ps1"}


def resolved_workspace(value):
    root = Path(value).expanduser().resolve()
    if any((character in str(root) for character in "\"'\r\n")):
        raise ValueError("Use a workspace path without quote or newline characters")
    protected = [CODE_ROOT, *SOURCE_DEFAULTS.values()]
    if any(
        (
            root == p.resolve()
            or root.is_relative_to(p.resolve())
            or p.resolve().is_relative_to(root)
            for p in protected
        )
    ):
        raise ValueError(
            "The run workspace must be separate from the code and original project directories"
        )
    return root


def find_rscript(value=None):
    candidates = [
        value,
        os.environ.get("RSCRIPT"),
        shutil.which("Rscript"),
        "D:/R/R-4.5.2/bin/Rscript.exe",
    ]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(Path(candidate).resolve())
    raise FileNotFoundError("Rscript was not found; use --rscript")


def source_paths(args):
    return {
        "SLE": Path(args.sle).expanduser().resolve(),
        "SLE_v1": Path(args.sle_v1).expanduser().resolve(),
        "SLE_R1": Path(args.sle_r1).expanduser().resolve(),
    }


def values(root, sources):
    result = {"runroot": root.as_posix(), "coderoot": CODE_ROOT.as_posix()}
    result.update({"source_" + key: path.as_posix() for key, path in sources.items()})
    return result


def expanded(value, replacements):
    result = str(value)
    for key, content in replacements.items():
        result = result.replace("{" + key + "}", content)
    return result


def environment(root, sources, extra=None):
    env = os.environ.copy()
    env.update(
        {
            "SLE_CODE_ROOT": str(CODE_ROOT),
            "SLE_RUN_ROOT": str(root),
            "SLE_WORK_ROOT": str(root),
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONUTF8": "1",
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
        }
    )
    for key, path in sources.items():
        env["SLE_SOURCE_" + key.upper()] = str(path)
    env["R_LIBS_USER"] = str(root / "environment" / "R_library")
    env["PYTHONPATH"] = str(CODE_ROOT) + os.pathsep + env.get("PYTHONPATH", "")
    for key, val in (extra or {}).items():
        env[key] = expanded(val, values(root, sources))
    return env
