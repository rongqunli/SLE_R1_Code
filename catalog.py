import sys

sys.dont_write_bytecode = True
import importlib.util
from configuration import CODE_ROOT


def catalog():
    stages = []
    inputs = []
    for path in sorted((CODE_ROOT / "steps").glob("*/plan.py")):
        spec = importlib.util.spec_from_file_location("sle_plan_" + path.parent.name, path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        stages.extend(module.STAGES)
        inputs.extend(getattr(module, "INPUTS", []))
    by_id = {}
    for stage in stages:
        if stage["id"] in by_id:
            raise ValueError("Duplicate stage: " + stage["id"])
        by_id[stage["id"]] = stage
        for item in stage.get("inputs", []):
            entry = dict(item)
            entry["stages"] = sorted(set(entry.get("stages", [])) | {stage["id"]})
            inputs.append(entry)
    for stage in stages:
        missing = set(stage.get("depends", [])) - set(by_id)
        if missing:
            raise ValueError(
                "Unknown dependencies for " + stage["id"] + ": " + ", ".join(sorted(missing))
            )
    ordered = []
    visiting = set()
    visited = set()

    def visit(key):
        if key in visiting:
            raise ValueError("Dependency cycle at " + key)
        if key in visited:
            return
        visiting.add(key)
        for parent in by_id[key].get("depends", []):
            visit(parent)
        visiting.remove(key)
        visited.add(key)
        ordered.append(by_id[key])

    for key in by_id:
        visit(key)
    unique = {}
    for item in inputs:
        owners = item.get("stages", [])
        if not isinstance(owners, list) or not owners:
            raise ValueError("External input has no consuming stage: " + item["target"])
        missing = set(owners) - set(by_id)
        if missing:
            raise ValueError(
                "Unknown input consumers for " + item["target"] + ": " + ", ".join(sorted(missing))
            )
        key = item["target"]
        entry = {name: value for name, value in item.items() if name != "stages"}
        if key in unique:
            previous = {name: value for name, value in unique[key].items() if name != "stages"}
            if previous != entry:
                raise ValueError("Conflicting input target: " + key)
            owners = list(set(unique[key]["stages"]) | set(owners))
        unique[key] = {**entry, "stages": sorted(set(owners))}
    return (ordered, list(unique.values()))


def select_stages(stages, requested):
    if not requested or requested == ["all"]:
        return stages
    by_id = {stage["id"]: stage for stage in stages}
    selected = set()

    def visit(key):
        if key not in by_id:
            raise ValueError("Unknown stage: " + key)
        if key in selected:
            return
        selected.add(key)
        for dependency in by_id[key].get("depends", []):
            visit(dependency)

    for key in requested:
        visit(key)
    return [stage for stage in stages if stage["id"] in selected]


def select_inputs(inputs, stages):
    selected = {stage["id"] for stage in stages}
    result = []
    for item in inputs:
        owners = item.get("stages", [])
        if not owners:
            raise ValueError("External input has no consuming stage: " + item["target"])
        if selected.intersection(owners):
            result.append(item)
    return result
