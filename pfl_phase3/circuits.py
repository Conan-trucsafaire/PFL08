"""Parsing of curated circuit strings into signed directed edges.

The frozen candidate table writes circuits as node/arrow chains, with four arrow
spellings in live use (`->`, `-\\>`, `→` for activation; `-|`, `┤` for inhibition) and
`;` separating several circuits reported by one publication.

Direction and sign are read from the arrows only. Protocol §2.1: a node list without
arrows is a component set — direction and closure are never inferred from token order.
"""

from __future__ import annotations

import re
from typing import Iterable, Sequence

# Order matters: `-\>` must be offered before `->` cannot match it, and `-|` before `┤`
# is irrelevant, but keeping the table explicit makes the sign map auditable.
ARROW_FORMS: dict[str, int] = {
    "->": +1,
    "-\\>": +1,
    "→": +1,
    "-|": -1,
    "┤": -1,
}

ARROW_RE = re.compile(r"\s*(" + "|".join(re.escape(a) for a in ARROW_FORMS) + r")\s*")

CIRCUIT_SEPARATOR = ";"


class CircuitParseError(ValueError):
    """Raised when a circuit string cannot be read as an arrow chain."""


def split_circuit(text: str) -> tuple[list[str], list[int]]:
    """Split one circuit string into its node tokens and its edge signs.

    Returns (nodes, signs) with len(signs) == len(nodes) - 1.
    """
    parts = ARROW_RE.split(text.strip())
    if len(parts) < 3:
        raise CircuitParseError(f"no arrow found in circuit: {text!r}")
    nodes = [parts[0].strip()]
    signs: list[int] = []
    for i in range(1, len(parts), 2):
        signs.append(ARROW_FORMS[parts[i]])
        nodes.append(parts[i + 1].strip())
    if any(not n for n in nodes):
        raise CircuitParseError(f"empty node token in circuit: {text!r}")
    return nodes, signs


def split_record(text: str) -> list[str]:
    """Split one table cell into the circuits it reports."""
    return [c.strip() for c in str(text).split(CIRCUIT_SEPARATOR) if c.strip()]


def sign_product(signs: Sequence[int]) -> int:
    """Product of edge signs. A curated positive loop must have product +1."""
    p = 1
    for s in signs:
        p *= s
    return p


def parse_circuit_table(records: Iterable[dict]) -> list[dict]:
    """Expand rows of the frozen candidate table into one dict per circuit.

    Each input record must carry at least `pmid` and `circuit`; `status`, `cancer`
    and `compartment` are carried through when present.
    """
    out: list[dict] = []
    for rec in records:
        for ci, circ in enumerate(split_record(rec.get("circuit", ""))):
            row = {
                "pmid": str(rec.get("pmid", "")).strip(),
                "status": rec.get("status", ""),
                "cancer": rec.get("cancer", ""),
                "compartment": rec.get("compartment", ""),
                "circuit_idx": ci,
                "raw": circ,
            }
            try:
                nodes, signs = split_circuit(circ)
            except CircuitParseError as exc:
                row.update(nodes=[], signs=[], parse_status="UNPARSED", parse_error=str(exc))
            else:
                row.update(
                    nodes=nodes,
                    signs=signs,
                    parse_status="OK",
                    parse_error="",
                    sign_product=sign_product(signs),
                )
            row["cycle_id"] = f"{row['pmid']}#{ci}"
            out.append(row)
    return out
