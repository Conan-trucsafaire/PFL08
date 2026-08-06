#!/usr/bin/env python3
"""Phase 3a — resource-free coverage audit (upper bounds on the frozen gates).

Governing protocol : docs/10_ANALYSIS_PROTOCOL_PHASE3_PERTURBATION_RECIPROCITY_FROZEN_2026-08-05.md
Runs on            : data/frozen/02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv
Requires           : nothing else. No perturbation resource, no network.

WHAT THIS DOES
    Phase 3a asks how many curated signed edges and cycles are measurable at all. Part
    of that answer does not depend on which screen is used: an edge whose endpoint is a
    cell-type-prefixed TME node, a pathway token, a metabolite, or a non-coding RNA is
    untestable in every single-cell-line perturbation screen, and an edge whose target
    is a protein-level readout is T-invisible in every mRNA screen. Those exclusions
    are computed here, giving a strict UPPER BOUND on each gate.

    An upper bound below a threshold settles that gate before 550 MB is downloaded.
    An upper bound above a threshold settles nothing — the gate stays UNDECIDED and
    needs the full run (scripts/run_phase3a.py).

WHY THE BOUND IS AN UPPER BOUND
    Node symbols are resolved in permissive mode (no gene universe available offline),
    so any symbol-shaped token is accepted provisionally. Real resolution against the
    resources' gene spaces can only ever reject more nodes, never fewer. Every count
    here is therefore >= the count the full run will produce.

ANTI-DRIFT (Charter §5)
    1. Question tested? None. Coverage only.
    2. Primary outcome? Upper bounds on the §3.3 gate quantities.
    3. What would refute? N/A — but an upper bound under a threshold binds Phase 3b to
       exploratory-only, exactly as a failed gate does.
    4. What does it exclude? That a later positive result rests on silently dropped edges.
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pfl_phase3 import FREEZE_DATE, PROTOCOL                      # noqa: E402
from pfl_phase3.circuits import parse_circuit_table               # noqa: E402
from pfl_phase3.gates import PERMITTED_REMEDY, upper_bound_ruling  # noqa: E402
from pfl_phase3.graph import build_edges_and_cycles               # noqa: E402
from pfl_phase3.nodes import resolve_node                         # noqa: E402
from pfl_phase3.testability import confirmatory_eligible, structural_reason  # noqa: E402

CIRCUIT_TSV = ROOT / "data" / "frozen" / "02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv"
OUT = ROOT / "results" / "phase3a_offline"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest().upper()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def counter_table(counts: Counter, label: str) -> str:
    if not counts:
        return f"  (none)"
    width = max(len(str(k)) for k in counts)
    return "\n".join(f"  {str(k):<{width}}  {v:>5d}" for k, v in counts.most_common())


def main() -> int:
    stamp = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"Phase 3a offline coverage audit · freeze {FREEZE_DATE} · run {stamp}")
    print(f"protocol: {PROTOCOL}")
    print("=" * 78)

    if not CIRCUIT_TSV.exists():
        print(f"STOP — frozen candidate table not found: {CIRCUIT_TSV}")
        return 1
    print(f"\ninput : {CIRCUIT_TSV.relative_to(ROOT)}")
    print(f"sha256: {sha256(CIRCUIT_TSV)}")

    # -- 1 · parse ---------------------------------------------------------------
    with open(CIRCUIT_TSV, encoding="utf-8") as f:
        records = list(csv.DictReader(f, delimiter="\t"))
    print(f"\nrecords: {len(records)}")
    print("status:")
    print(counter_table(Counter(r["status"] for r in records), "status"))

    circuits = parse_circuit_table(records)
    unparsed = [c for c in circuits if c["parse_status"] != "OK"]
    print(f"\ncircuits: {len(circuits)}   parsed OK: {len(circuits) - len(unparsed)}")
    for c in unparsed:
        print(f"  UNPARSED  {c['pmid']}  {c['raw']}  ({c['parse_error']})")

    ok = [c for c in circuits if c["parse_status"] == "OK"]
    bad_polarity = [c for c in ok if c["sign_product"] != +1]
    print(f"\nsign product +1 (positive loop): {sum(1 for c in ok if c['sign_product'] == 1)}"
          f" / {len(ok)}")
    for c in bad_polarity:
        print(f"  NEGATIVE-PRODUCT  {c['cycle_id']}  {c['raw']}")

    # -- 2 · resolve nodes (permissive) ------------------------------------------
    resolved = {}
    node_rows = []
    for c in ok:
        for pos, tok in enumerate(c["nodes"]):
            r = resolve_node(tok, known=None)   # permissive: upper-bound mode
            resolved[(c["cycle_id"], pos)] = r
            node_rows.append({**r, "cycle_id": c["cycle_id"], "pmid": c["pmid"],
                              "position": pos, "compartment": c["compartment"]})

    print(f"\nnode tokens: {len(node_rows)}   unique: "
          f"{len({r['raw_token'] for r in node_rows})}")
    print("node classes (permissive):")
    print(counter_table(Counter(r["node_class"] for r in node_rows), "class"))
    print("node visibility:")
    print(counter_table(Counter(r["visibility"] for r in node_rows), "visibility"))

    # -- 3 · edges and cycles ----------------------------------------------------
    edges, cycles = build_edges_and_cycles(ok, resolved)
    for e in edges:
        e["structural_reason"] = structural_reason(e) or "structurally_eligible"
        e["confirmatory_eligible"] = confirmatory_eligible(e)

    print(f"\nedges: {len(edges)}   cycles: {len(cycles)}")
    print("cycle closure (first node == last node):")
    print(counter_table(Counter(c["closure"] for c in cycles), "closure"))
    print("\nedge visibility (target-determined, protocol §3.2):")
    vis = Counter(e["visibility"] for e in edges)
    print(counter_table(vis, "visibility"))
    t_inv_share = vis["T_invisible"] / len(edges) if edges else 0
    print(f"  T_invisible share: {t_inv_share:.1%}"
          f"  (catalogue-level figure quoted in protocol §3.2: 44.4%)")

    print("\nRESOURCE-INDEPENDENT ATTRITION (holds in every context):")
    print(counter_table(Counter(e["structural_reason"] for e in edges), "reason"))

    # -- 4 · upper bounds --------------------------------------------------------
    eligible = [e for e in edges if e["structural_reason"] == "structurally_eligible"]
    conf_eligible = [e for e in edges if e["confirmatory_eligible"]]

    eligible_by_cycle = {}
    for e in edges:
        eligible_by_cycle.setdefault(e["cycle_id"], []).append(
            e["structural_reason"] == "structurally_eligible")
    closed_ids = {c["cycle_id"] for c in cycles if c["closure"] == "CLOSED"}
    complete_closed = [cid for cid, flags in eligible_by_cycle.items()
                       if cid in closed_ids and all(flags)]

    # Three nested ceilings for gate 1, so the verdict does not depend on how strict the
    # automated visibility classifier happens to be:
    #   B1  as classified          — T_visible AND structurally eligible
    #   B2  every T_uncertain promoted to T_visible by the human reviewer (§2.2)
    #   B3  absolute ceiling: every edge except the intercellular ones, which cannot
    #       close inside a single cell line under any resource
    b1 = len(conf_eligible)
    b2 = len([e for e in eligible if e["visibility"] != "T_invisible"])
    b3 = len([e for e in edges if e["structural_reason"] != "intercellular"])

    upper = {
        # B2 is used for the verdict: it already grants the human reviewer every
        # visibility promotion, so a CANNOT_PASS here survives the second review.
        "t_visible_testable_edges_pooled": b2,
        "completely_testable_cycles": len(complete_closed),
        "contexts_with_50plus_testable_edges": 2 if len(eligible) >= 50 else 0,
        "cycles_testable_in_2plus_contexts": len(complete_closed),
    }
    rows, ruling = upper_bound_ruling(upper)

    print("\nGATE 1 CEILINGS (T-visible testable edges pooled, threshold 150)")
    print(f"  B1  as classified (T_visible & eligible)                  {b1:>5}")
    print(f"  B2  + every T_uncertain promoted by the human reviewer    {b2:>5}   <- used for the verdict")
    print(f"  B3  absolute ceiling (all but intercellular edges)        {b3:>5}   "
          f"reachable only by relaxing the detectability rule, which §2.3 prohibits")

    print("\n" + "=" * 78)
    print("GATE UPPER BOUNDS (protocol §3.3 thresholds, not editable)")
    print("=" * 78)
    print(f"{'gate':<42}{'thr':>5}{'upper':>8}   verdict")
    for r in rows:
        ub = "n/a" if r["upper_bound"] is None else str(r["upper_bound"])
        print(f"{r['gate']:<42}{r['threshold']:>5}{ub:>8}   {r['verdict']}")
    print(f"\nRULING: {ruling}")
    if ruling.startswith("EXPLORATORY_ONLY"):
        print("\n  At least one gate cannot be met by any resource coverage.")
        print("  Permitted remedy (protocol §2.3):")
        print(f"    {PERMITTED_REMEDY}")
        print("  Prohibited: lowering a threshold, relaxing the detectability rule,")
        print("  substituting inferred L1000 genes, or swapping in enumerated cycles.")

    # -- 5 · deliverables --------------------------------------------------------
    write_tsv(OUT / "phase3a_offline_edge_inventory.tsv", edges, [
        "cycle_id", "pmid", "edge_idx", "source_token", "target_token", "source", "target",
        "sign", "source_class", "target_class", "visibility", "structural_reason",
        "confirmatory_eligible", "compartment", "status", "review_flag"])
    write_tsv(OUT / "phase3a_offline_cycle_inventory.tsv", cycles, [
        "cycle_id", "pmid", "status", "cancer", "compartment", "length", "closure",
        "sign_product", "n_unmapped_nodes", "raw"])
    write_tsv(OUT / "phase3a_offline_mrna_visibility_declaration.tsv", edges, [
        "cycle_id", "edge_idx", "source_token", "target_token", "visibility",
        "target_class", "review_flag"])
    # The §2.2 reviewer adjudicates *tokens*, not edge instances: one ruling on
    # "MAPK14/p38 activity" settles every edge that touches it. In permissive mode every
    # symbol is unverified, so an edge-level worklist would just be the whole edge list.
    by_token: dict[str, dict] = {}
    for r in node_rows:
        t = by_token.setdefault(r["raw_token"], {
            "raw_token": r["raw_token"], "occurrences": 0,
            "node_class": r["node_class"], "provisional_symbol": r["symbol"],
            "visibility": r["visibility"], "note": r["note"],
            "cycles": set(), "reviewer_symbol": "", "reviewer_visibility": "",
            "reviewer_decision": "",
        })
        t["occurrences"] += 1
        t["cycles"].add(r["cycle_id"])
    worklist = sorted(by_token.values(),
                      key=lambda t: (-t["occurrences"], t["raw_token"]))
    for t in worklist:
        t["cycles"] = ";".join(sorted(t["cycles"]))
    write_tsv(OUT / "phase3a_offline_NODE_REVIEW_WORKLIST.tsv", worklist, [
        "raw_token", "occurrences", "node_class", "provisional_symbol", "visibility",
        "note", "cycles", "reviewer_symbol", "reviewer_visibility", "reviewer_decision"])
    print(f"\nnode-level review worklist: {len(worklist)} unique tokens "
          f"(vs {len(node_rows)} token instances)")

    decision = {
        "phase": "3a-offline",
        "freeze_date": FREEZE_DATE,
        "run_stamp": stamp,
        "governing_protocol": PROTOCOL,
        "outcome_blind": True,
        "any_effect_value_read": False,
        "resources_consulted": [],
        "input": {"path": str(CIRCUIT_TSV.relative_to(ROOT)), "sha256": sha256(CIRCUIT_TSV),
                  "records": len(records)},
        "counts": {
            "circuits_parsed": len(ok),
            "circuits_unparsed": len(unparsed),
            "edges": len(edges),
            "cycles": len(cycles),
            "cycles_closed": len(closed_ids),
            "edges_structurally_eligible": len(eligible),
            "edges_T_visible_and_eligible": len(conf_eligible),
            "cycles_closed_and_fully_eligible": len(complete_closed),
            "edge_visibility": dict(vis),
            "structural_attrition": dict(Counter(e["structural_reason"] for e in edges)),
            "gate1_ceilings": {
                "B1_as_classified": b1,
                "B2_uncertain_promoted": b2,
                "B3_absolute_non_intercellular": b3,
            },
        },
        "gate_upper_bounds": rows,
        "ruling": ruling,
        "second_human_review_completed": False,
        "binding_notes": [
            "Symbols were resolved permissively (no gene universe offline); every count "
            "is an UPPER BOUND on what the full Phase 3a run can produce.",
            "Gate thresholds are frozen in protocol §3.3 and were not modified.",
            "A gate marked CANNOT_PASS is settled; a gate marked UNDECIDED still "
            "requires scripts/run_phase3a.py against R1/R2/R4/R6.",
            "Node resolution is automated and REQUIRES the human second review of "
            "protocol §2.2 before any confirmatory claim.",
        ],
    }
    (OUT / "phase3a_offline_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"\noutputs written to {OUT.relative_to(ROOT)}/")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name:<52}{f.stat().st_size:>9,d} B")
    print("\nNothing here may be reused to score edges — Phase 3b needs its own analysis lock.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
