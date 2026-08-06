"""The Phase 3a gates.

Thresholds are frozen in protocol §3.3 and may not be edited. If a gate fails, the only
permitted remedy (§2.3) is extending signed-edge adjudication to more of the 248
un-adjudicated legacy records — literature curation performed without reference to any
perturbation outcome, completed and re-frozen before any effect value is read.

Prohibited remedies: relaxing the detectability rule; substituting inferred L1000 genes
for landmark genes; replacing curated cycles with database-enumerated cycles in the
confirmatory tier; lowering the context requirement.
"""

from __future__ import annotations

from types import MappingProxyType

GATES = MappingProxyType({
    "t_visible_testable_edges_pooled": 150,
    "completely_testable_cycles": 25,
    "contexts_with_50plus_testable_edges": 2,
    "nontransformed_contexts": 1,
    "cycles_testable_in_2plus_contexts": 10,
})

PERMITTED_REMEDY = (
    "Extend signed-edge adjudication to the 248 un-adjudicated legacy records, "
    "complete and re-freeze it before any effect value is read, then re-run Phase 3a."
)


def evaluate_gates(observed: dict) -> tuple[list[dict], str]:
    """Compare observed counts against the frozen thresholds.

    Returns (rows, ruling). Every gate in GATES must be present in `observed`; a missing
    gate is an error rather than a silent pass, because a silently skipped gate would
    read as a met gate in the decision file.
    """
    missing = [k for k in GATES if k not in observed]
    if missing:
        raise KeyError(f"gate(s) not evaluated: {missing}")

    rows = [{"gate": k, "threshold": v, "observed": int(observed[k]),
             "pass": int(observed[k]) >= v}
            for k, v in GATES.items()]
    ruling = "CONFIRMATORY_PERMITTED" if all(r["pass"] for r in rows) else "EXPLORATORY_ONLY"
    return rows, ruling


def upper_bound_ruling(upper: dict) -> tuple[list[dict], str]:
    """Rule on *upper bounds* rather than observed counts.

    A gate whose upper bound already sits below its threshold cannot be met by any
    resource coverage, so the ruling is decidable before the resources are downloaded.
    Gates whose upper bound clears the threshold are undecided, not passed.
    """
    rows = []
    for k, thr in GATES.items():
        if k not in upper:
            rows.append({"gate": k, "threshold": thr, "upper_bound": None,
                         "verdict": "NOT_BOUNDED_OFFLINE"})
            continue
        ub = int(upper[k])
        rows.append({"gate": k, "threshold": thr, "upper_bound": ub,
                     "verdict": "CANNOT_PASS" if ub < thr else "UNDECIDED"})
    ruling = ("EXPLORATORY_ONLY_ALREADY_DETERMINED"
              if any(r["verdict"] == "CANNOT_PASS" for r in rows)
              else "UNDETERMINED_PENDING_RESOURCES")
    return rows, ruling
