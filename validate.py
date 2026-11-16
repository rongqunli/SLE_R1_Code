import sys

sys.dont_write_bytecode = True
import argparse
import ast
import io
import json
import subprocess
import tokenize
from pathlib import Path
from configuration import CODE_ROOT, CODE_SUFFIXES, find_rscript


def python_validation(path):
    source = path.read_text(encoding="utf-8-sig")
    tree = ast.parse(source, filename=str(path))
    compile(tree, str(path), "exec")
    comments = [
        item.start[0]
        for item in tokenize.generate_tokens(io.StringIO(source).readline)
        if item.type == tokenize.COMMENT
    ]
    docs = []
    for item in ast.walk(tree):
        if isinstance(
            item, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)
        ) and ast.get_docstring(item, clean=False):
            docs.append(getattr(item, "lineno", 1))
    return {
        "file": str(path.relative_to(CODE_ROOT)),
        "parsed": True,
        "comment_lines": comments,
        "docstring_lines": docs,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--rscript")
    parser.add_argument("--report")
    args = parser.parse_args()
    other = []
    python = []
    for path in sorted(CODE_ROOT.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix not in CODE_SUFFIXES:
            other.append(str(path.relative_to(CODE_ROOT)))
        elif path.suffix == ".py":
            python.append(python_validation(path))
    command = [
        find_rscript(args.rscript),
        "--vanilla",
        str(CODE_ROOT / "validate.R"),
        str(CODE_ROOT),
    ]
    process = subprocess.run(
        command, capture_output=True, text=True, encoding="utf-8", errors="replace"
    )
    result = {
        "python_files": python,
        "unexpected_files": other,
        "r_returncode": process.returncode,
        "r_stdout": process.stdout,
        "r_stderr": process.stderr,
    }
    if args.report:
        path = Path(args.report).resolve()
        if path.is_relative_to(CODE_ROOT):
            raise ValueError("Validation reports must be outside the code directory")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                "python_files": len(python),
                "unexpected_files": other,
                "python_comments": sum((len(x["comment_lines"]) for x in python)),
                "python_docstrings": sum((len(x["docstring_lines"]) for x in python)),
                "r_returncode": process.returncode,
                "r_summary": process.stdout,
            },
            ensure_ascii=False,
        )
    )
    if (
        other
        or process.returncode
        or any((x["comment_lines"] or x["docstring_lines"] for x in python))
    ):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
