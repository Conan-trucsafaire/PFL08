"""Phase 3 — perturbation reciprocity and closure.

Governing protocol: docs/10_ANALYSIS_PROTOCOL_PHASE3_PERTURBATION_RECIPROCITY_FROZEN_2026-08-05.md
Charter:            00_PROJECT_CHARTER_FROZEN_2026-08-05.md (Drive, unchanged, still binding)

This package holds the *logic* of Phase 3a — circuit parsing, node resolution,
mRNA-visibility declaration, testability and the frozen gates — separated from the
Colab runner so that every rule the protocol states can be unit-tested without
touching a perturbation resource.

Nothing in this package reads an effect value. See pfl_phase3.blind.
"""

FREEZE_DATE = "2026-08-05"
PROTOCOL = "10_ANALYSIS_PROTOCOL_PHASE3_PERTURBATION_RECIPROCITY_FROZEN_2026-08-05.md"

__all__ = ["FREEZE_DATE", "PROTOCOL"]
