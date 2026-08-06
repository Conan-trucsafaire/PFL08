"""Node resolution and the prespecified mRNA-visibility declaration (protocol §3.2).

Two rules govern everything here, and both are deliberately conservative:

1.  **Visibility is a property of the edge's TARGET, not its source.** The readout for
    edge A→B is B's transcript, so it is B's mechanism descriptor that decides whether
    the edge can produce an mRNA signal at all. `MYC protein stability` as a target is
    T_invisible no matter who A is. Getting this backwards silently destroys the
    T-invisible negative control, which is the falsification arm of the whole design.

2.  **The classifier never promotes.** T_visible is asserted only when the descriptor
    names a transcript-level readout outright. Anything ambiguous is T_uncertain.
    Over-calling T_visible would both dilute the confirmatory set with edges that
    cannot produce a signal and weaken the negative control. Promotion to T_visible is
    a human act, recorded in the manual-review worklist — never a code act.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional

# --------------------------------------------------------------------------------------
# token vocabularies
# --------------------------------------------------------------------------------------

COMPLEX_TOKENS = [
    "mTORC1", "MAPK/ERK", "MAPK", "ERK", "NF-kappaB", "NF-κB", "JAK-STAT",
    "JAK2/STAT3", "TGF-beta pathway", "TGF-beta", "PI3K", "AKT", "Wnt", "IKK",
    "TLR2/TLR4", "IL1R", "CXCR2", "beta-catenin/TCF4", "TCF7L2/beta-catenin", "p38",
]

PROCESS_TOKENS = [
    "glycolysis", "glycolytic", "lactate", "lactylation", "NADPH",
    "pentose phosphate", "cholesterol synthesis", "m6A", "phase separation",
]

NONCODING_RE = re.compile(
    r"(^|[^A-Za-z0-9])(miR-|let-7|circ|lnc|LINC\d|MIR\d|SNHG|HCG\d|PVT1|CASC\d|LUCAT|RN7SK|"
    r"[A-Z0-9]+-AS\d?$|TPT1-AS|CSF3R-AS|MIR\d+HG|HEGBC)", re.I)

CELLTYPE_RE = re.compile(
    r"^(macrophage|tumou?r-cell|tumou?r|neuronal|mesothelial|CAF-like-cell|MDSC|stromal|immune)\b",
    re.I)

SYMBOL_RE = re.compile(r"\b([A-Z][A-Z0-9]{1,9}(?:[A-Z0-9-]{0,4}))\b")

# Amendment A1 (2026-08-06, made outcome-blind — no perturbation resource was consulted).
#
# The first-pass classifier matched complex/pathway tokens by bare prefix
# (`token.lower().startswith("wnt")`), which swallowed every real gene whose symbol
# begins with a pathway name: WNT7B -> "Wnt", AKT1/2/3 -> "AKT", MAPK1 -> "MAPK",
# PI3KCA -> "PI3K". Those nodes were emitted as `complex_or_pathway` with no gene
# assigned, so their edges were counted as permanently untestable and any cycle they
# closed was read as unclosed. That is a matching defect, not a curation choice.
#
# The fix requires a token boundary on both sides. The vocabulary itself is the
# author's curation and is left untouched.
COMPLEX_RE = [
    (c, re.compile(r"(?<![A-Za-z0-9])" + re.escape(c) + r"(?![A-Za-z0-9])", re.I))
    for c in COMPLEX_TOKENS
]


def _clean_symbol_hits(token: str) -> list[str]:
    """Symbol-shaped substrings of `token`, with trailing hyphens stripped.

    SYMBOL_RE's optional `[A-Z0-9-]{0,4}` tail is greedy across hyphens and `\\b` is
    satisfied straight after one, so `IL-1beta` yielded `IL-` and `DKK1-LRP6` yielded
    `DKK1-`. Neither is a gene symbol; both would silently occupy the symbol slot.
    """
    hits = []
    for h in SYMBOL_RE.findall(token):
        h = h.rstrip("-")
        if len(h) >= 2 and h not in hits:
            hits.append(h)
    return hits

# mechanism descriptor -> mRNA visibility of THIS node AS A READOUT TARGET
T_INVISIBLE_KEYS = [
    "protein stability", "protein abundance", "protein-interaction",
    "protein interaction", "phosphorylation", "palmitoylation", "lactylation",
    "methylation", "nuclear localization", "nuclear/transcriptional",
    "nuclear activity", "phase separation", "access to", "e3 activity",
    "translation", "secondary-structure", "receptor state",
    "methyltransferase activity", "promoter occupancy", "recruitment",
]

T_VISIBLE_KEYS = [
    "transcription", "mrna", "expression", "promoter activity", "biogenesis",
]

#: Node classes that make an edge untestable in a single-cell-line perturbation screen
#: regardless of which resource is used. These are resource-independent exclusions.
STRUCTURAL_EXCLUSION_CLASSES = (
    "cell_type_prefixed",
    "noncoding_rna",
    "complex_or_pathway",
    "metabolite_or_process",
    "unmapped",
)


# --------------------------------------------------------------------------------------
# visibility
# --------------------------------------------------------------------------------------

def classify_visibility(token: str) -> str:
    """Classify a node token as T_visible / T_invisible / T_uncertain (protocol §3.2)."""
    t = token.lower()
    inv = any(k in t for k in T_INVISIBLE_KEYS)
    vis = any(k in t for k in T_VISIBLE_KEYS)
    if vis and not inv:
        return "T_visible"
    if inv and not vis:
        return "T_invisible"
    return "T_uncertain"


# --------------------------------------------------------------------------------------
# node resolution
# --------------------------------------------------------------------------------------

def resolve_node(token: str, known: Optional[Iterable[str]] = None) -> dict:
    """Resolve one node token to a gene symbol and a node class.

    `known` is the gene universe (union of every resource's perturbed and measured
    spaces). Passing ``None`` selects **permissive mode**: any symbol-shaped token is
    accepted provisionally and the class is suffixed ``_unverified``. Permissive mode
    exists only for the resource-free audit, where the true universe is unavailable;
    it can only ever over-count mappable nodes, so counts derived from it are upper
    bounds and must be reported as such.
    """
    permissive = known is None
    universe = set() if permissive else set(known)

    t = token.strip()
    out = {
        "raw_token": t,
        "symbol": None,
        "node_class": None,
        "visibility": classify_visibility(t),
        "review_flag": False,
        "note": "",
    }

    if CELLTYPE_RE.match(t):
        out.update(node_class="cell_type_prefixed", review_flag=True,
                   note="intercellular/TME; cannot close inside one cell line")
        return out

    if any(p.lower() in t.lower() for p in PROCESS_TOKENS) and not SYMBOL_RE.search(t.split()[0]):
        out.update(node_class="metabolite_or_process", review_flag=True)
        return out

    if NONCODING_RE.search(t):
        out.update(node_class="noncoding_rna", review_flag=True,
                   note="non-coding; confirm membership in the screen's perturbation list")
        return out

    for c, rx in COMPLEX_RE:
        if rx.search(t):
            out.update(node_class="complex_or_pathway", review_flag=True,
                       note=f"complex/pathway token '{c}'; no gene assigned")
            return out

    candidates = _clean_symbol_hits(t)
    hits = candidates if permissive else [h for h in candidates if h in universe]
    suffix = "_unverified" if permissive else ""

    if len(hits) == 1:
        out.update(symbol=hits[0], node_class="single_symbol" + suffix,
                   review_flag=permissive,
                   note="provisional: no gene universe available" if permissive else "")
    elif len(hits) > 1:
        out.update(symbol=hits[0], node_class="multi_symbol" + suffix, review_flag=True,
                   note=f"multiple symbols {hits}; first taken provisionally")
    else:
        out.update(node_class="unmapped", review_flag=True,
                   note="no symbol-shaped token" if permissive
                        else "no symbol matched the gene universe")
    return out


def is_structurally_excluded(node_class: Optional[str]) -> bool:
    """True when this node class makes any incident edge untestable in every context."""
    return node_class in STRUCTURAL_EXCLUSION_CLASSES
