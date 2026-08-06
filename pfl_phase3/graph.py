"""Assembly of parsed circuits + resolved nodes into signed edges and cycles."""

from __future__ import annotations

import math
from typing import Any, Optional


def sym(x: Any) -> Optional[str]:
    """None-safe symbol read — the single place a symbol is turned into str-or-None.

    pandas turns ``None`` into ``NaN``, and ``NaN`` is both truthy and not-None. A naive
    ``x is not None`` test therefore marks every node that failed to map as mappable,
    which quietly pulls complexes, non-coding RNAs and TME nodes into the denominator.
    Every symbol read goes through here.
    """
    if x is None:
        return None
    if isinstance(x, float) and math.isnan(x):
        return None
    s = str(x).strip()
    return s if s and s.lower() != "nan" else None


def build_edges_and_cycles(circuits: list[dict], resolved: dict) -> tuple[list[dict], list[dict]]:
    """Turn parsed circuits into edge rows and cycle rows.

    `resolved` maps (cycle_id, position) -> the dict returned by nodes.resolve_node.

    Closure is asserted only when the first and last node resolve to the *same* symbol.
    An unresolved endpoint therefore never counts as closed.
    """
    edge_rows: list[dict] = []
    cycle_rows: list[dict] = []

    for c in circuits:
        if c.get("parse_status") != "OK":
            continue
        cid = c["cycle_id"]
        n = len(c["nodes"])
        res = [resolved[(cid, p)] for p in range(n)]
        syms = [sym(r["symbol"]) for r in res]
        closes = syms[0] is not None and syms[0] == syms[-1]

        for i in range(n - 1):
            s, t = res[i], res[i + 1]
            ssym, tsym = syms[i], syms[i + 1]
            edge_rows.append({
                "cycle_id": cid,
                "pmid": c["pmid"],
                "edge_idx": i,
                "source_token": s["raw_token"],
                "target_token": t["raw_token"],
                "source": ssym,
                "target": tsym,
                "sign": c["signs"][i],
                "source_class": s["node_class"],
                "target_class": t["node_class"],
                # protocol §3.2: visibility of an edge is the visibility of its TARGET
                "visibility": t["visibility"],
                "compartment": c["compartment"],
                "status": c["status"],
                "review_flag": bool(s["review_flag"] or t["review_flag"]),
                "mappable": ssym is not None and tsym is not None,
            })

        cycle_rows.append({
            "cycle_id": cid,
            "pmid": c["pmid"],
            "status": c["status"],
            "cancer": c["cancer"],
            "compartment": c["compartment"],
            "length": n - 1,
            "closure": "CLOSED" if closes else "NOT_CLOSED",
            "sign_product": c.get("sign_product"),
            "n_unmapped_nodes": int(sum(s is None for s in syms)),
            "raw": c["raw"],
        })

    return edge_rows, cycle_rows
