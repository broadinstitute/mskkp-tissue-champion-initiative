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
    """Load a dataset's obs (cell metadata) table from disk.

    Args:
        path (str | pathlib.Path): Path to a ``.h5ad`` AnnData file or a
            ``.csv`` export of obs. The file extension determines how it's
            read; any other extension is rejected.

    Returns:
        pandas.DataFrame: The obs table, one row per cell.

    Raises:
        SystemExit: If the file extension is neither ``.h5ad`` nor ``.csv``,
            or if reading a ``.h5ad`` file is requested but the optional
            ``anndata`` dependency isn't installed.
    """
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
    """Map a champion_criteria.yaml field reference to an actual obs column name.

    Field references are written as e.g. ``"obs.nCount_RNA"``; this strips
    the leading ``"obs."`` and, for ``percent_mt``, also checks the
    alternate dotted spelling (``"percent.mt"``) that Seurat-derived
    datasets commonly use (see ``PERCENT_MT_ALIASES``).

    Args:
        obs (pandas.DataFrame): The obs table to look up the column in.
        field (str): A field reference from champion_criteria.yaml, e.g.
            ``"obs.nCount_RNA"`` or ``"obs.percent_mt"``. A bare column name
            without the ``"obs."`` prefix is also accepted.

    Returns:
        str | None: The matching column name in ``obs.columns``, or
        ``None`` if no matching column exists.
    """
    # field looks like "obs.nCount_RNA" — strip the leading "obs."
    name = field.split(".", 1)[1] if field.startswith("obs.") else field
    if name == "percent_mt":
        for alias in PERCENT_MT_ALIASES:
            if alias in obs.columns:
                return alias
        return None
    return name if name in obs.columns else None


def eval_when(obs, when):
    """Filter obs down to the rows matched by a rule's `when` condition.

    Only the ``obs.<field> == "<value>"`` shape is supported, since that's
    the only form champion_criteria.yaml uses (e.g. gating the mitochondrial
    thresholds on ``obs.assay == "snRNA-seq"``).

    Args:
        obs (pandas.DataFrame): The obs table to filter.
        when (str): A condition string of the form
            ``'obs.<field> == "<value>"'``.

    Returns:
        pandas.DataFrame: The subset of ``obs`` where the condition holds.
        Empty (zero rows) if the referenced field/column doesn't exist in
        ``obs``.
    """
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
    """The outcome of evaluating a single champion_criteria.yaml rule.

    Attributes:
        rule_id (str): The rule's ``id`` from champion_criteria.yaml (e.g.
            ``"min_counts"``).
        status (str): ``"required"`` or ``"recommended"``, copied from the
            rule definition (or hardcoded for the special-cased rules).
        passed (bool): Whether the rule's condition held for this dataset.
        detail (str): A human-readable explanation of the computed value(s)
            and why the rule passed or failed, printed in the report.
    """

    def __init__(self, rule_id, status, passed, detail):
        """Construct a Result.

        Args:
            rule_id (str): The rule's ``id`` from champion_criteria.yaml.
            status (str): ``"required"`` or ``"recommended"``.
            passed (bool): Whether the rule's condition held.
            detail (str): Human-readable explanation for the report.

        Returns:
            None
        """
        self.rule_id = rule_id
        self.status = status
        self.passed = passed
        self.detail = detail


def check_simple_rule(rule, obs):
    """Evaluate a field + aggregate + operator + threshold rule from champion_criteria.yaml.

    Handles every rule shape except ``annotation_depth``,
    ``celltype_well_defined``, and ``control_and_disease``, which have their
    own dedicated check functions and are dispatched separately by
    ``run_all``. If the rule has a ``when`` condition and no rows match it,
    the rule is treated as vacuously satisfied.

    Args:
        rule (dict): A single rule dict from champion_criteria.yaml, with
            keys ``id`` (str), ``status`` (str), ``aggregate`` (str, one
            of ``AGGREGATES``), ``field`` (str, e.g. ``"obs.nCount_RNA"``),
            ``operator`` (str, one of ``OPERATORS``), ``threshold``
            (int | float), and optionally ``when`` (str).
        obs (pandas.DataFrame): The obs table to evaluate the rule against.

    Returns:
        Result: The outcome of the rule for this dataset.
    """
    subset = eval_when(obs, rule["when"]) if "when" in rule else obs
    if subset.empty and "when" in rule:
        return Result(rule["id"], rule["status"], True, "no rows match `when` — vacuously satisfied")

    aggregate = rule["aggregate"]
    if aggregate == "n_rows":
        value = len(subset)
    else:
        col = resolve_column(subset, rule["field"])
        if col is None:
            return Result(rule["id"], rule["status"], False, f"column for {rule['field']} not found in obs")
        value = AGGREGATES[aggregate](subset[col])

    op = OPERATORS[rule["operator"]]
    passed = bool(op(value, rule["threshold"]))
    return Result(
        rule["id"], rule["status"], passed,
        f"{aggregate}({rule['field']}) = {value} {rule['operator']} {rule['threshold']} -> {'PASS' if passed else 'FAIL'}",
    )


def check_annotation_depth(obs):
    """Check the `annotation_depth` rule: both author_celltype and author_cell_subtype are present.

    Args:
        obs (pandas.DataFrame): The obs table to check. Must be checked for
            the presence of non-null values in the ``author_celltype`` and
            ``author_cell_subtype`` columns.

    Returns:
        Result: Always ``status="required"``. Passes if both the
        ``author_celltype`` and ``author_cell_subtype`` columns exist in
        ``obs`` and each has at least one non-null value.
    """
    has_celltype = "author_celltype" in obs.columns and obs["author_celltype"].notna().any()
    has_subtype = "author_cell_subtype" in obs.columns and obs["author_cell_subtype"].notna().any()
    levels = int(has_celltype) + int(has_subtype)
    passed = levels >= 2
    return Result("annotation_depth", "required", passed, f"levels present = {levels} (need >= 2)")


def check_celltype_well_defined(obs, min_cells_per_subtype):
    """Check the `celltype_well_defined` rule: clean nesting and minimum group size.

    Verifies that every ``author_cell_subtype`` value maps to at most one
    ``author_celltype`` value (clean Level-2-under-Level-1 nesting), and that
    every ``author_cell_subtype`` group has at least ``min_cells_per_subtype``
    cells.

    Args:
        obs (pandas.DataFrame): The obs table to check.
        min_cells_per_subtype (int): The minimum number of cells required
            per distinct ``author_cell_subtype`` value, typically the
            tissue's ``TissueProfile.min_cells_per_subtype`` (default 20).

    Returns:
        Result: Always ``status="required"``. Fails if the
        ``author_celltype`` or ``author_cell_subtype`` columns are missing,
        if any ``author_cell_subtype`` maps to more than one
        ``author_celltype``, or if any ``author_cell_subtype`` group has
        fewer than ``min_cells_per_subtype`` cells.
    """
    if "author_celltype" not in obs.columns or "author_cell_subtype" not in obs.columns:
        return Result("celltype_well_defined", "required", False, "author_celltype/author_cell_subtype columns missing")

    nests_ok = (obs.groupby("author_cell_subtype")["author_celltype"].nunique() <= 1).all()
    sizes = obs.groupby("author_cell_subtype").size()
    sizes_ok = bool((sizes >= min_cells_per_subtype).all())
    passed = bool(nests_ok and sizes_ok)
    smallest = sizes.min() if len(sizes) else 0
    detail = (
        f"every author_cell_subtype maps to <=1 author_celltype: {nests_ok}; "
        f"smallest author_cell_subtype has {smallest} cells (need >= {min_cells_per_subtype}): {sizes_ok}"
    )
    return Result("celltype_well_defined", "required", passed, detail)


def check_control_and_disease(obs):
    """Check the `control_and_disease` rule: Control plus at least one other group.

    Args:
        obs (pandas.DataFrame): The obs table to check.

    Returns:
        Result: Always ``status="recommended"`` (disease samples are only
        expected "if available"). Passes if the ``disease_status`` column
        exists and its non-null values include ``"Control"`` plus at least
        one other distinct value.
    """
    if "disease_status" not in obs.columns:
        return Result("control_and_disease", "recommended", False, "disease_status column missing")
    values = obs["disease_status"].dropna().unique().tolist()
    has_control = "Control" in values
    passed = has_control and len(values) >= 2
    return Result("control_and_disease", "recommended", passed, f"distinct values = {values}")


def run_all(obs, criteria, min_cells_per_subtype):
    """Evaluate every rule in champion_criteria.yaml against a dataset's obs.

    Dispatches each rule to its dedicated check function by ``id``
    (``annotation_depth``, ``celltype_well_defined``,
    ``control_and_disease``), falling back to ``check_simple_rule`` for
    every other rule.

    Args:
        obs (pandas.DataFrame): The obs table to evaluate all rules against.
        criteria (list[dict]): The ``champion_criteria`` list loaded from
            champion_criteria.yaml.
        min_cells_per_subtype (int): Passed through to
            ``check_celltype_well_defined``.

    Returns:
        list[Result]: One Result per rule in ``criteria``, in the same order.
    """
    results = []
    for rule in criteria:
        rid = rule["id"]
        if rid == "annotation_depth":
            results.append(check_annotation_depth(obs))
        elif rid == "celltype_well_defined":
            results.append(check_celltype_well_defined(obs, min_cells_per_subtype))
        elif rid == "control_and_disease":
            results.append(check_control_and_disease(obs))
        else:
            results.append(check_simple_rule(rule, obs))
    return results


def main():
    """CLI entry point: load obs + rule files, run all checks, and print a report.

    Parses command-line arguments (``--obs``, ``--tissue``, ``--criteria``,
    ``--tissue-profiles``, ``--min-cells-per-subtype``), loads the obs table
    and rule files, resolves ``min_cells_per_subtype`` (from the CLI flag or
    else the tissue's profile), runs every rule via ``run_all``, and prints
    a PASS/FAIL line per rule followed by an overall verdict.

    Args:
        None. Arguments are read from ``sys.argv`` via ``argparse``.

    Returns:
        None. This function never returns normally — it always terminates
        the process via ``sys.exit`` (see Raises).

    Raises:
        SystemExit: With status 0 if every ``required`` rule passed, or
            status 1 if any ``required`` rule failed. Also raised earlier,
            with a non-zero status and an error message, if ``--obs`` has
            an unsupported extension, if reading a ``.h5ad`` file requires
            the missing ``anndata`` package, or if ``--tissue`` doesn't
            match any key in the tissue_profiles.yaml file.
    """
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
        outcome = "PASS" if r.passed else "FAIL"
        if not r.passed and r.status == "required":
            any_required_failed = True
        print(f"  [{outcome}] {r.rule_id} ({r.status})\n         {r.detail}")

    print()
    if any_required_failed:
        print("Result: NOT a champion dataset — one or more required rules failed.")
        sys.exit(1)
    print("Result: meets all required champion criteria.")
    sys.exit(0)


if __name__ == "__main__":
    main()
