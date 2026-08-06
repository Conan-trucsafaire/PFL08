#!/usr/bin/env python3
"""Phase 3a — outcome-blind feasibility gate (full run, needs the perturbation resources).

Governing protocol : docs/10_ANALYSIS_PROTOCOL_PHASE3_PERTURBATION_RECIPROCITY_FROZEN_2026-08-05.md
Charter            : 00_PROJECT_CHARTER_FROZEN_2026-08-05.md

Run in Colab with:

    from google.colab import drive; drive.mount('/content/drive')
    !git clone https://github.com/Conan-trucsafaire/PFL08 /content/PFL08
    %cd /content/PFL08
    !PFL_INPUTS=/content/drive/MyDrive/PFL08Claude版/00_inputs python scripts/run_phase3a.py

WHAT THIS DOES
    Counts how many curated signed edges and cycles are measurable at all in the frozen
    perturbation resources. It tests no hypothesis and reads no effect value.

BLINDNESS CONTRACT (enforced by code, not convention)
    pfl_phase3.blind.BlindH5 raises on any attempt to touch /X, /layers, /obsm, /varm,
    /obsp in an .h5ad, and the guard is self-tested at startup. DepMap is read with
    nrows=0 (header only). Every file access is logged to disk.
    Do NOT replace BlindH5 with anndata.read_h5ad — that loads X.

ANTI-DRIFT (Charter §5)
    1. Question tested? None. Feasibility only.
    2. Primary outcome? Counts of testable edges / cycles / contexts.
    3. What would refute? N/A — but a gate failure downgrades Phase 3b to exploratory,
       bindingly.
    4. What does it exclude? That a later positive result rests on silently dropped edges.
"""

from __future__ import annotations

import csv
import datetime
import hashlib
import json
import os
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pfl_phase3 import FREEZE_DATE, PROTOCOL                       # noqa: E402
from pfl_phase3.blind import AccessLog, BlindH5, self_test          # noqa: E402
from pfl_phase3.circuits import parse_circuit_table                 # noqa: E402
from pfl_phase3.gates import PERMITTED_REMEDY, evaluate_gates       # noqa: E402
from pfl_phase3.graph import build_edges_and_cycles, sym            # noqa: E402
from pfl_phase3.nodes import resolve_node                           # noqa: E402
from pfl_phase3.testability import CIS_WINDOW_BP, testable_in       # noqa: E402

RUN_STAMP = datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

IN = Path(os.environ.get("PFL_INPUTS", ROOT / "data" / "inputs"))
OUT = Path(os.environ.get("PFL_RESULTS", ROOT / "results" / "phase3a"))

CIRCUIT_TSV = ROOT / "data" / "frozen" / "02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv"
REPLOGLE = {
    "K562_gwps": IN / "K562_gwps_normalized_bulk_01.h5ad",
    "RPE1_essential": IN / "rpe1_normalized_bulk_01.h5ad",
    "K562_essential": IN / "K562_essential_normalized_bulk_01.h5ad",
}
DEPMAP_CRISPR = IN / "CRISPRGeneEffect.csv"
SOUTHARD_GENES = IN / "southard2025_crispra_perturbed_genes.txt"

ACCESS_LOG = AccessLog()
CLEAN_OUT = {"nan", "none", "non-targeting", "control", ""}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest().upper()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ======================================================================================
# 1 · input check — clean stop, not a traceback
# ======================================================================================

def check_inputs() -> bool:
    required = [
        (CIRCUIT_TSV, "in this repository (data/frozen/)"),
        (REPLOGLE["K562_gwps"], "Figshare+ DOI 10.25452/figshare.plus.20029387 (~375 MB)"),
        (REPLOGLE["RPE1_essential"], "Figshare+ DOI 10.25452/figshare.plus.20029387 (~95 MB)"),
    ]
    missing = [(p, s) for p, s in required if not p.exists()]
    if not missing:
        print("\n✓ required inputs present")
        return True

    print("\n" + "!" * 78)
    print("STOP — required inputs are missing. Put them in:")
    print(f"  {IN}")
    print("!" * 78)
    for p, s in missing:
        print(f"\n  [ ] {p.name}\n      from: {s}")
    print("\nOptional but recommended (absence is recorded, never silently skipped):")
    for p, s in [(REPLOGLE["K562_essential"], "same Figshare+ record (~80 MB)"),
                 (DEPMAP_CRISPR, "DepMap Public 26Q1 download page"),
                 (SOUTHARD_GENES, "Southard 2025, Zenodo 15213597 — one gene symbol per line")]:
        if not p.exists():
            print(f"  [ ] {p.name}\n      from: {s}")
    print("\nDownload ONLY the three *_normalized_bulk_01.h5ad files (~550 MB total).")
    print("Do NOT download the single-cell versions (65.8 GB) — nothing here uses them.")
    return False


# ======================================================================================
# 2 · resource gene spaces (blind)
# ======================================================================================

def gene_space(path: Path, label: str) -> dict:
    with BlindH5(path, ACCESS_LOG) as h:
        var_cols, obs_cols = h.columns("var"), h.columns("obs")
        print(f"\n[{label}] var cols: {var_cols}")
        print(f"[{label}] obs cols: {obs_cols}")
        measured = h.index_of("var")
        for c in ("gene_name", "gene_symbol", "symbol", "gene_names"):
            if c in var_cols:
                measured = h.column("var", c)
                print(f"[{label}] measured symbols <- var['{c}']")
                break
        pert = None
        for c in ("gene", "gene_symbol", "perturbation", "target_gene", "gene_id"):
            if c in obs_cols:
                pert = h.column("obs", c)
                print(f"[{label}] perturbed symbols <- obs['{c}']")
                break
        if pert is None:
            pert = h.index_of("obs")
            print(f"[{label}] ! no obs gene column found — using obs index as a fallback.")
            print(f"[{label}] ! CHECK the obs cols printed above before trusting this run.")

    def clean(seq):
        return {str(x).strip() for x in seq if x and str(x).strip().lower() not in CLEAN_OUT}

    return {"measured": clean(measured), "perturbed": clean(pert)}


def build_contexts() -> tuple[dict, set]:
    contexts = {}
    for label, path in REPLOGLE.items():
        if path.exists():
            contexts[label] = {
                "modality": "CRISPRi", "direction": "LoF",
                "cell_line": "K562" if label.startswith("K562") else "RPE1",
                "transformed": label.startswith("K562"),
                "tier": "primary" if label in ("K562_gwps", "RPE1_essential") else "sensitivity",
                **gene_space(path, label)}

    if SOUTHARD_GENES.exists():
        ACCESS_LOG.record(SOUTHARD_GENES, "gene list")
        contexts["RPE1_CRISPRa"] = {
            "modality": "CRISPRa", "direction": "GoF", "cell_line": "RPE1",
            "transformed": False, "tier": "primary",
            "perturbed": {l.strip() for l in open(SOUTHARD_GENES) if l.strip()},
            "measured": contexts.get("RPE1_essential", {}).get("measured", set()),
            "note": "measured space approximated by RPE1 CRISPRi var index — replace with the deposited var"}

    depmap_genes: set = set()
    if DEPMAP_CRISPR.exists():
        import pandas as pd
        ACCESS_LOG.record(DEPMAP_CRISPR, "header row only (nrows=0)")
        depmap_genes = {re.sub(r"\s*\(\d+\)$", "", c).strip()
                        for c in pd.read_csv(DEPMAP_CRISPR, nrows=0).columns[1:]}
    return contexts, depmap_genes


# ======================================================================================
# 3 · cis-proximity exclusion (protocol §3.1 rule 4)
# ======================================================================================

def ensembl_coords(symbols, batch: int = 250) -> tuple[dict, bool]:
    import requests
    coords, ok = {}, True
    syms = sorted({s for s in (sym(x) for x in symbols) if s})
    for i in range(0, len(syms), batch):
        try:
            r = requests.post("https://rest.ensembl.org/lookup/symbol/homo_sapiens",
                              headers={"Content-Type": "application/json",
                                       "Accept": "application/json"},
                              json={"symbols": syms[i:i + batch]}, timeout=60)
            r.raise_for_status()
            for s, d in r.json().items():
                if isinstance(d, dict) and "seq_region_name" in d:
                    coords[s] = (str(d["seq_region_name"]), int(d["start"]), int(d["end"]))
        except Exception as exc:
            print(f"  Ensembl lookup failed (batch {i // batch}): {exc}")
            ok = False
    return coords, ok


def cis_conflict(a, b, coords, ok):
    a, b = sym(a), sym(b)
    if not ok or a not in coords or b not in coords:
        return None
    ca, sa, ea = coords[a]
    cb, sb, eb = coords[b]
    if ca != cb:
        return False
    return (max(sa, sb) - min(ea, eb)) < CIS_WINDOW_BP


# ======================================================================================
# main
# ======================================================================================

def main() -> int:
    IN.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    print("=" * 78)
    print(f"Phase 3a · freeze {FREEZE_DATE} · run {RUN_STAMP}")
    print(f"protocol: {PROTOCOL}")
    print(f"inputs  : {IN}")
    print(f"outputs : {OUT}")
    print("=" * 78)

    if not check_inputs():
        return 2

    self_test(REPLOGLE["K562_gwps"], ACCESS_LOG)
    print("✓ blindness guard active (/X refused)")

    # -- resource manifest -------------------------------------------------------
    spec = [("R1", "Replogle K562 gwps", REPLOGLE["K562_gwps"], True),
            ("R2", "Replogle RPE1 essential", REPLOGLE["RPE1_essential"], True),
            ("R3", "Replogle K562 essential", REPLOGLE["K562_essential"], False),
            ("R4", "Southard CRISPRa gene list", SOUTHARD_GENES, False),
            ("R6", "DepMap 26Q1 header", DEPMAP_CRISPR, False),
            ("--", "Curated signed-edge circuits", CIRCUIT_TSV, True)]
    manifest = [{"resource_id": rid, "name": nm, "required": req, "path": str(p),
                 "present": p.exists(),
                 "bytes": p.stat().st_size if p.exists() else "",
                 "sha256": sha256(p) if p.exists() else "",
                 "coverage_failure_reason": "" if p.exists() else "not downloaded at run time"}
                for rid, nm, p, req in spec]
    write_tsv(OUT / "phase3a_resource_manifest.tsv", manifest,
              ["resource_id", "name", "required", "path", "present", "bytes", "sha256",
               "coverage_failure_reason"])
    for m in manifest:
        print(f"  {m['resource_id']:<4}{m['name']:<32}present={m['present']!s:<6}{m['bytes']}")

    contexts, depmap_genes = build_contexts()
    print()
    for k, v in contexts.items():
        print(f"  {k:16s} perturbed={len(v['perturbed']):6d} measured={len(v['measured']):6d} {v['tier']}")
    print(f"  DepMap gene space: {len(depmap_genes)}")

    # -- circuits ----------------------------------------------------------------
    with open(CIRCUIT_TSV, encoding="utf-8") as f:
        records = list(csv.DictReader(f, delimiter="\t"))
    ACCESS_LOG.record(CIRCUIT_TSV, "identifiers only")
    circuits = parse_circuit_table(records)
    ok_circuits = [c for c in circuits if c["parse_status"] == "OK"]
    print(f"\ncircuit table: {len(records)} records")
    print(f"circuits parsed OK: {len(ok_circuits)} / {len(circuits)}")
    for c in circuits:
        if c["parse_status"] != "OK":
            print(f"  UNPARSED {c['pmid']}  {c['raw']}")

    # -- node resolution ---------------------------------------------------------
    universe = set(depmap_genes)
    for v in contexts.values():
        universe |= v["measured"] | v["perturbed"]
    print(f"\nsymbol universe: {len(universe)}")
    if not universe:
        print("STOP — empty gene universe; every node would resolve as unmapped.")
        return 3

    resolved, node_rows = {}, []
    for c in ok_circuits:
        for pos, tok in enumerate(c["nodes"]):
            r = resolve_node(tok, known=universe)
            resolved[(c["cycle_id"], pos)] = r
            node_rows.append({**r, "cycle_id": c["cycle_id"], "pmid": c["pmid"],
                              "position": pos, "compartment": c["compartment"],
                              "status": c["status"]})
    write_tsv(OUT / "phase3a_node_resolution_REVIEW_REQUIRED.tsv", node_rows,
              ["cycle_id", "pmid", "position", "raw_token", "symbol", "node_class",
               "visibility", "review_flag", "note", "compartment", "status"])
    print("\nnode classes:")
    for k, v in Counter(r["node_class"] for r in node_rows).most_common():
        print(f"  {k:<26}{v:>5}")
    print("node visibility:")
    for k, v in Counter(r["visibility"] for r in node_rows).most_common():
        print(f"  {k:<26}{v:>5}")
    print(f"nodes flagged for human review: {sum(r['review_flag'] for r in node_rows)} / {len(node_rows)}")

    # -- edges, cycles, cis ------------------------------------------------------
    edges, cycles = build_edges_and_cycles(ok_circuits, resolved)
    coords, cis_ok = ensembl_coords({e["source"] for e in edges} | {e["target"] for e in edges})
    print(f"\ncoordinates resolved: {len(coords)}   lookup complete: {cis_ok}")
    for e in edges:
        e["cis_conflict"] = cis_conflict(e["source"], e["target"], coords, cis_ok)

    print(f"\nedges: {len(edges)}   cycles: {len(cycles)}")
    for k, v in Counter(c["closure"] for c in cycles).most_common():
        print(f"  {k:<14}{v:>5}")
    print("edge visibility:")
    for k, v in Counter(e["visibility"] for e in edges).most_common():
        print(f"  {k:<14}{v:>5}")

    # -- testability per context -------------------------------------------------
    tested = []
    for cname, ctx in contexts.items():
        for e in edges:
            ok_e, why = testable_in(ctx, e)
            tested.append({"context": cname, "cycle_id": e["cycle_id"], "edge_idx": e["edge_idx"],
                           "source": e["source"], "target": e["target"], "sign": e["sign"],
                           "visibility": e["visibility"], "testable": ok_e, "reason": why,
                           "tier": ctx["tier"], "direction": ctx["direction"],
                           "transformed": ctx["transformed"]})
    write_tsv(OUT / "phase3a_edge_testability_lock.tsv", tested,
              ["context", "cycle_id", "edge_idx", "source", "target", "sign", "visibility",
               "testable", "reason", "tier", "direction", "transformed"])

    print("\nATTRITION (all contexts pooled)")
    for k, v in Counter(t["reason"] for t in tested).most_common():
        print(f"  {k:<24}{v:>6}")

    # -- cycle completeness ------------------------------------------------------
    closure_of = {c["cycle_id"]: c["closure"] for c in cycles}
    per_cycle = {}
    for t in tested:
        per_cycle.setdefault((t["context"], t["cycle_id"]), []).append(t)
    cyc_rows = []
    for (cname, cid), g in per_cycle.items():
        cyc_rows.append({"context": cname, "cycle_id": cid, "n_edges": len(g),
                         "n_testable": sum(x["testable"] for x in g),
                         "complete": all(x["testable"] for x in g),
                         "n_T_visible": sum(x["visibility"] == "T_visible" for x in g),
                         "closure": closure_of.get(cid, "")})
    write_tsv(OUT / "phase3a_cycle_testability_lock.tsv", cyc_rows,
              ["context", "cycle_id", "closure", "n_edges", "n_testable", "complete", "n_T_visible"])
    closed_complete = [r for r in cyc_rows if r["complete"] and r["closure"] == "CLOSED"]

    # -- gates -------------------------------------------------------------------
    per_ctx = Counter(t["context"] for t in tested if t["testable"])
    vis_pairs = {(t["cycle_id"], t["edge_idx"]) for t in tested
                 if t["testable"] and t["visibility"] == "T_visible"}
    cycle_ctx = Counter(r["cycle_id"] for r in closed_complete)
    observed = {
        "t_visible_testable_edges_pooled": len(vis_pairs),
        "completely_testable_cycles": len({r["cycle_id"] for r in closed_complete}),
        "contexts_with_50plus_testable_edges": sum(1 for v in per_ctx.values() if v >= 50),
        "nontransformed_contexts": len({t["context"] for t in tested
                                        if t["testable"] and not t["transformed"]}),
        "cycles_testable_in_2plus_contexts": sum(1 for v in cycle_ctx.values() if v >= 2),
    }
    gate_rows, ruling = evaluate_gates(observed)
    if not cis_ok:
        ruling += "_PROVISIONAL_CIS_CHECK_NOT_PERFORMED"

    review = [e for e in edges if e["review_flag"]]
    write_tsv(OUT / "phase3a_mrna_visibility_declaration.tsv", edges,
              ["cycle_id", "edge_idx", "source_token", "target_token", "source", "target",
               "sign", "visibility", "target_class", "review_flag"])
    write_tsv(OUT / "phase3a_MANUAL_REVIEW_WORKLIST.tsv", review,
              ["cycle_id", "pmid", "edge_idx", "source_token", "target_token", "source",
               "target", "sign", "source_class", "target_class", "visibility"])
    write_tsv(OUT / "phase3a_access_log.tsv", list(ACCESS_LOG), ["ts", "path", "read"])

    decision = {
        "phase": "3a", "freeze_date": FREEZE_DATE, "run_stamp": RUN_STAMP,
        "governing_protocol": PROTOCOL,
        "outcome_blind": True, "any_effect_value_read": False, "cis_check_performed": cis_ok,
        "contexts": {k: {"tier": v["tier"], "direction": v["direction"],
                         "cell_line": v["cell_line"], "transformed": v["transformed"],
                         "n_perturbed": len(v["perturbed"]), "n_measured": len(v["measured"])}
                     for k, v in contexts.items()},
        "gates": gate_rows, "ruling": ruling, "second_human_review_completed": False,
        "binding_notes": [
            "Gate thresholds are frozen in protocol §3.3 and were not modified in this run.",
            f"If EXPLORATORY_ONLY: {PERMITTED_REMEDY}",
            "Node resolution is automated and REQUIRES the human second review of "
            "protocol §2.2 before any confirmatory claim.",
        ],
    }
    (OUT / "phase3a_gate_decision.json").write_text(
        json.dumps(decision, indent=2, ensure_ascii=False), encoding="utf-8")

    print("\n" + "=" * 78)
    print("GATE RESULT")
    print("=" * 78)
    print(f"{'gate':<42}{'thr':>5}{'obs':>8}   pass")
    for r in gate_rows:
        print(f"{r['gate']:<42}{r['threshold']:>5}{r['observed']:>8}   {r['pass']}")
    print(f"\nRULING: {ruling}")
    if ruling.startswith("EXPLORATORY"):
        print("\n  -> Phase 3b may run but produces NO confirmatory claim. This downgrade is binding.")
        print(f"  -> Permitted remedy: {PERMITTED_REMEDY}")
        print("  -> NOT: lowering a threshold here.")
    print(f"\nEdges needing human adjudication: {len(review)} / {len(edges)}")
    print(f"\nOutputs written to: {OUT}")
    for f in sorted(OUT.iterdir()):
        print(f"  {f.name:<52}{f.stat().st_size:>9,d} B")
    print("\nNEXT: read phase3a_gate_decision.json, then work phase3a_MANUAL_REVIEW_WORKLIST.tsv.")
    print("Nothing here may be reused to score edges — Phase 3b needs its own analysis lock.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
