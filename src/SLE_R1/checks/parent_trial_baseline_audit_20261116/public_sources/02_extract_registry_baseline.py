import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import csv, json

ROOT = Path(__file__).resolve().parent
rows = []
overview = []
for nct in ["NCT01196091", "NCT01205438"]:
    data = json.loads((ROOT / f"{nct}.json").read_text(encoding="utf-8"))
    results = data["resultsSection"]
    b = results["baselineCharacteristicsModule"]
    groups = {x["id"]: x for x in b["groups"]}
    defaultden = {x["groupId"]: x["value"] for d in b.get("denoms", []) for x in d["counts"]}
    measureindex = []
    for mi, m in enumerate(b["measures"]):
        mden = {x["groupId"]: x["value"] for d in m.get("denoms", []) for x in d["counts"]}
        measureindex.append({k: v for k, v in m.items() if k != "classes"})
        for ci, c in enumerate(m.get("classes", [])):
            for ki, cat in enumerate(c.get("categories", [])):
                for ji, x in enumerate(cat.get("measurements", [])):
                    gid = x["groupId"]
                    rows.append(
                        dict(
                            nct=nct,
                            source_url=f"https://clinicaltrials.gov/study/{nct}?tab=results",
                            source_json=f"{nct}.json",
                            json_pointer=f"/resultsSection/baselineCharacteristicsModule/measures/{mi}/classes/{ci}/categories/{ki}/measurements/{ji}",
                            population_description=b.get("populationDescription", ""),
                            measure=m["title"],
                            description=m.get("description", ""),
                            param_type=m.get("paramType", ""),
                            dispersion_type=m.get("dispersionType", ""),
                            units=m.get("unitOfMeasure", ""),
                            measure_denom_json=json.dumps(m.get("denoms", [])),
                            class_title=c.get("title", ""),
                            category=cat.get("title", ""),
                            group_id=gid,
                            group_title=groups[gid]["title"],
                            default_baseline_n=defaultden.get(gid, ""),
                            explicit_measure_n=mden.get(gid, ""),
                            value=x.get("value", ""),
                            spread=x.get("spread", ""),
                            lower_limit=x.get("lowerLimit", ""),
                            upper_limit=x.get("upperLimit", ""),
                            comment=x.get("comment", ""),
                        )
                    )
    overview.append(
        dict(
            nct=nct,
            protocol_identification=data["protocolSection"]["identificationModule"],
            protocol_enrollment=data["protocolSection"]["designModule"].get("enrollmentInfo"),
            status=data["protocolSection"].get("statusModule"),
            baseline_population_description=b.get("populationDescription", ""),
            baseline_groups=b["groups"],
            baseline_denominators=b.get("denoms", []),
            baseline_measure_definitions=measureindex,
            participant_flow=results.get("participantFlowModule", {}),
            outcome_populations=[
                {
                    "title": m["title"],
                    "population": m.get("populationDescription", ""),
                    "denoms": m.get("denoms", []),
                }
                for m in results.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
            ],
        )
    )
with (ROOT / "registered_baseline_long.csv").open("w", encoding="utf-8-sig", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
(ROOT / "registry_population_and_definitions.json").write_text(
    json.dumps(overview, indent=2, ensure_ascii=False), encoding="utf-8"
)
with (ROOT / "registered_baseline_compact.txt").open("w", encoding="utf-8") as f:
    for trial in overview:
        nct = trial["nct"]
        f.write(
            f"\n{nct}\n{trial['baseline_population_description']}\n{trial['baseline_denominators']}\n"
        )
        for r in rows:
            if r["nct"] != nct or r["group_title"] != "Total":
                continue
            f.write(
                f"{r['measure']} | {r['class_title']} | {r['category']}: {r['value']} ({r['spread']}) [{r['param_type']};{r['dispersion_type']};{r['units']}; measure_n={r['explicit_measure_n']};default_n={r['default_baseline_n']}]\n"
            )
        f.write("\nPARTICIPANT FLOW\n" + json.dumps(trial["participant_flow"], indent=2))
print(f"Wrote {len(rows)} source rows without statistical transformation.")
