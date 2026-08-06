"""Testability of a curated edge, per protocol §3.1.

An edge A→B is testable in context c when all of:

  1. A is a perturbed gene in c's perturbation list;
  2. B is in c's measured gene space (landmark only for LINCS);
  3. B passes c's own published detection criterion;
  4. A and B are > 500 kb apart or on different chromosomes (CRISPRi cis exclusion);
  5. the curated edge has an assigned sign;
  6. the edge is mRNA-plausible under §3.2.

Attrition is attributed to the **most specific** cause available, so that the report
distinguishes "this edge can never be tested by any perturbation screen" from "this
particular screen does not carry the gene".
"""

from __future__ import annotations

from typing import Optional

from .graph import sym
from .nodes import STRUCTURAL_EXCLUSION_CLASSES

#: attrition reasons that hold in every context — no resource can rescue them
STRUCTURAL_REASONS = ("intercellular", "noncoding_rna", "complex_or_pathway",
                      "metabolite_or_process", "unmapped_endpoint", "cis_proximity")

CIS_WINDOW_BP = 500_000


def structural_reason(edge: dict) -> Optional[str]:
    """Resource-independent reason this edge can never be tested, or None.

    Does not consider any resource's gene space, so the answer is the same for every
    context and can be computed before a single byte of perturbation data exists.
    """
    classes = (edge.get("source_class"), edge.get("target_class"))
    if "cell_type_prefixed" in classes:
        return "intercellular"
    for cls in ("noncoding_rna", "complex_or_pathway", "metabolite_or_process"):
        if cls in classes:
            return cls
    if not bool(edge.get("mappable")) or sym(edge.get("source")) is None or sym(edge.get("target")) is None:
        return "unmapped_endpoint"
    if edge.get("cis_conflict") is True:
        return "cis_proximity"
    if edge.get("sign") not in (+1, -1):
        return "unsigned_edge"
    return None


def testable_in(ctx: dict, edge: dict) -> tuple[bool, str]:
    """Is this edge testable in context `ctx`? Returns (testable, reason)."""
    reason = structural_reason(edge)
    if reason is not None:
        return False, reason
    src, tgt = sym(edge["source"]), sym(edge["target"])
    if src not in ctx["perturbed"]:
        return False, "source_not_perturbed"
    if tgt not in ctx["measured"]:
        return False, "target_not_measured"
    return True, "testable"


def confirmatory_eligible(edge: dict) -> bool:
    """T-visible and structurally eligible — the pool the confirmatory set is drawn from.

    Protocol §3.2: only T-visible edges enter confirmatory inference. T-invisible edges
    are retained as the prespecified negative control (E8), not discarded.
    """
    return edge.get("visibility") == "T_visible" and structural_reason(edge) is None


__all__ = [
    "CIS_WINDOW_BP", "STRUCTURAL_REASONS", "STRUCTURAL_EXCLUSION_CLASSES",
    "structural_reason", "testable_in", "confirmatory_eligible",
]
