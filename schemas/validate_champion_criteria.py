#!/usr/bin/env python3
"""
Champion-criteria validator for the MSKKP RNA-seq gold-standard.

champion_criteria.yaml encodes aggregate/statistical checks over obs
(min/max/n_unique/n_rows, gated by `when` conditions) that LinkML itself
cannot express — this script is the "separate validator" that
mskkp_rnaseq_schema.yaml's header comment refers to.

It only looks at obs. Schema-level structural checks (required fields,
dtypes, enum membership) belong to linkml-validate against
mskkp_rnaseq_schema.yaml, not here.

Usage:
    python validate_champion_criteria.py --obs dataset.h5ad --tissue skeletal_muscle
    python validate_champion_criteria.py --obs obs.csv --tissue cartilage --min-cells-per-subtype 15

Exit code is 1 if any `required` rule fails, 0 otherwise. `recommended`
rules are reported but never affect the exit code.
"""

import argparse
import sys
from pathlib import Path

import yaml

try:
    import pandas as pd
except ImportError:
    sys.exit("This script requires pandas. Install with: pip install pandas")

SCHEMA_DIR = Path(__file__).resolve().parent

# obs.percent_mt is the LinkML slot name; real datasets (Seurat convention)
# often spell it with a dot instead of an underscore.
PERCENT_MT_ALIASES = ["percent_mt", "percent.mt"]


def load_obs(path):
    path = Path(path)
    if path.suffix == ".h5ad":
        try:
            import anndata
        except ImportError:
            sys.exit("Reading .h5ad requires anndata. Install with: pip install anndata")
        return anndata.read_h5ad(path).obs
    if path.suffix == ".csv":
        return pd.read_csv(path)
    sys.exit(f"Unsupported obs file type: {path.suffix} (use .h5ad or .csv)")


def resolve_column(obs, field):
    # field looks like "obs.nCount_RNA" — strip the leading "obs."
    name = field.split(".", 1)[1] if field.startswith("obs.") else field
    if name == "percent_mt":
        for alias in PERCENT_MT_ALIASES:
            if alias in obs.columns:
                return alias
        return None
    return name if name in obs.columns else None


def eval_when(obs, when):
    # only "obs.<field> == \"<value>\"" is supported — that's the only shape
    # champion_criteria.yaml uses.
    field, _, value = when.partition("==")
    field = field.strip()
    value = value.strip().strip('"').strip("'")
    col = resolve_column(obs, field)
    if col is None:
        return obs.iloc[0:0]
    return obs[obs[col] == value]


OPERATORS = {
    ">": lambda a, b: a > b,
    ">=": lambda a, b: a >= b,
    "<": lambda a, b: a < b,
    "<=": lambda a, b: a <= b,
}

AGGREGATES = {
    "min": lambda s: s.min(),
    "max": lambda s: s.max(),
    "n_unique": lambda s: s.nunique(dropna=True),
    "n_rows": lambda s: len(s),
}


class Result:
    def __init__(self, rule_id, severity, passed, detail):
        self.rule_id = rule_id
        self.severity = severity
        self.passed = passed
        self.detail = detail


def check_simple_rule(rule, obs):
    subset = eval_when(obs, rule["when"]) if "when" in rule else obs
    if subset.empty and "when" in rule:
        return Result(rule["id"], rule["severity"], True, "no rows match `when` — vacuously satisfied")

    aggregate = rule["aggregate"]
    if aggregate == "n_rows":
        value = len(subset)
    else:
        col = resolve_column(subset, rule["field"])
        if col is None:
            return Result(rule["id"], rule["severity"], False, f"column for {rule['field']} not found in obs")
        value = AGGREGATES[aggregate](subset[col])

    op = OPERATORS[rule["operator"]]
    passed = bool(op(value, rule["threshold"]))
    return Result(
        rule["id"], rule["severity"], passed,
        f"{aggregate}({rule['field']}) = {value} {rule['operator']} {rule['threshold']} -> {'PASS' if passed else 'FAIL'}",
    )


def check_annotation_depth(obs):
    has_celltype = "celltype" in obs.columns and obs["celltype"].notna().any()
    has_subtype = "cell_subtype" in obs.columns and obs["cell_subtype"].notna().any()
    levels = int(has_celltype) + int(has_subtype)
    passed = levels >= 2
    return Result("annotation_depth", "required", passed, f"levels present = {levels} (need >= 2)")


def check_clusters_well_defined(obs, min_cells_per_subtype):
    if "celltype" not in obs.columns or "cell_subtype" not in obs.columns:
        return Result("clusters_well_defined", "required", False, "celltype/cell_subtype columns missing")

    nests_ok = (obs.groupby("cell_subtype")["celltype"].nunique() <= 1).all()
    sizes = obs.groupby("cell_subtype").size()
    sizes_ok = bool((sizes >= min_cells_per_subtype).all())
    passed = bool(nests_ok and sizes_ok)
    smallest = sizes.min() if len(sizes) else 0
    detail = (
        f"every cell_subtype maps to <=1 celltype: {nests_ok}; "
        f"smallest cell_subtype has {smallest} cells (need >= {min_cells_per_subtype}): {sizes_ok}"
    )
    return Result("clusters_well_defined", "required", passed, detail)


def check_control_and_disease(obs):
    if "disease_status" not in obs.columns:
        return Result("control_and_disease", "recommended", False, "disease_status column missing")
    values = obs["disease_status"].dropna().unique().tolist()
    has_control = "Control" in values
    passed = has_control and len(values) >= 2
    return Result("control_and_disease", "recommended", passed, f"distinct values = {values}")


def run_all(obs, criteria, min_cells_per_subtype):
    results = []
    for rule in criteria:
        rid = rule["id"]
        if rid == "annotation_depth":
            results.append(check_annotation_depth(obs))
        elif rid == "clusters_well_defined":
            results.append(check_clusters_well_defined(obs, min_cells_per_subtype))
        elif rid == "control_and_disease":
            results.append(check_control_and_disease(obs))
        else:
            results.append(check_simple_rule(rule, obs))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--obs", required=True, help="Path to a .h5ad dataset or a .csv export of obs.")
    parser.add_argument("--tissue", required=True, help="tissue_id, e.g. skeletal_muscle — looked up in tissue_profiles.yaml for min_cells_per_subtype.")
    parser.add_argument("--criteria", default=str(SCHEMA_DIR / "champion_criteria.yaml"))
    parser.add_argument("--tissue-profiles", default=str(SCHEMA_DIR / "tissue_profiles.yaml"))
    parser.add_argument("--min-cells-per-subtype", type=int, default=None, help="Override the tissue_profile's min_cells_per_subtype.")
    args = parser.parse_args()

    obs = load_obs(args.obs)

    with open(args.criteria) as f:
        criteria = yaml.safe_load(f)["champion_criteria"]

    min_cells_per_subtype = args.min_cells_per_subtype
    if min_cells_per_subtype is None:
        with open(args.tissue_profiles) as f:
            profiles = yaml.safe_load(f)["tissue_profiles"]
        profile = profiles.get(args.tissue)
        if profile is None:
            sys.exit(f"Unknown tissue_id {args.tissue!r}. Known: {', '.join(profiles)}")
        min_cells_per_subtype = profile.get("min_cells_per_subtype", 20)

    results = run_all(obs, criteria, min_cells_per_subtype)

    print(f"Champion criteria for tissue={args.tissue} ({len(obs)} cells)\n")
    any_required_failed = False
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        if not r.passed and r.severity == "required":
            any_required_failed = True
        tag = "required" if r.severity == "required" else "recommended"
        print(f"  [{status}] {r.rule_id} ({tag})\n         {r.detail}")

    print()
    if any_required_failed:
        print("Result: NOT a champion dataset — one or more required rules failed.")
        sys.exit(1)
    print("Result: meets all required champion criteria.")
    sys.exit(0)


if __name__ == "__main__":
    main()
