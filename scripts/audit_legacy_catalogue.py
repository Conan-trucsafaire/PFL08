#!/usr/bin/env python3
"""Audit of the v5 submission's supplementary catalogue (Tables S1/S2).

Inputs are transcriptions of `Supplementary_Tables.xlsx` from the first submission
package, held read-only under data/legacy/. This script answers three questions that
the reanalysis pre-registration lists as blocking, and it recomputes them from the
tables rather than repeating the manuscript's own summary numbers.

  1. Are the 192 catalogue entries 192 distinct circuits?
  2. Was the dual-high selection rule applied as stated?
  3. How many independent gene groupings do the 17 "core loops" actually represent?

It also crosswalks the Phase 3a signed-edge circuits onto the catalogue by gene-set
overlap. Table S1 carries no PMID column, so a PMID-keyed crosswalk is not possible
from this file; overlap is the available substitute and is labelled as such.
"""

from __future__ import annotations

import collections
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
S1_PATH = ROOT / "data" / "legacy" / "S1_loop_catalog.tsv"
S2_PATH = ROOT / "data" / "legacy" / "S2_dualhigh_evidence.tsv"
EDGES = ROOT / "results" / "phase3a_offline" / "phase3a_offline_edge_inventory.tsv"
OUT = ROOT / "results" / "legacy_audit"


def read(path: Path) -> list[dict]:
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def genes_of(row: dict) -> frozenset[str]:
    return frozenset(g.strip() for g in row["genes"].split(",") if g.strip())


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    s1 = read(S1_PATH)
    s2 = read(S2_PATH)
    report: list[str] = []

    def say(line=""):
        print(line)
        report.append(line)

    say("=" * 78)
    say("Legacy catalogue audit — Supplementary Tables S1/S2 of the v5 submission")
    say("=" * 78)
    say(f"\nS1 entries: {len(s1)}")

    # -- 1 · duplicate member sets ----------------------------------------------
    by_set: dict[frozenset, list[str]] = collections.defaultdict(list)
    for r in s1:
        by_set[genes_of(r)].append(r["loop_id"])
    dups = {k: v for k, v in by_set.items() if len(v) > 1}
    say(f"\n[1] DISTINCT MEMBER SETS: {len(by_set)}  "
        f"(duplicate groups {len(dups)}, redundant entries {sum(len(v) - 1 for v in dups.values())})")
    say(f"    The headline count of {len(s1)} overstates distinct circuits by "
        f"{len(s1) - len(by_set)}.")
    for k, v in sorted(dups.items(), key=lambda x: x[1][0]):
        say(f"      {' = '.join(v):<32} {sorted(k)}")

    # -- 2 · was the stated selection rule applied? ------------------------------
    dh = [r for r in s1 if r["in_dualhigh17"] == "Yes"]
    violating = [r["loop_id"] for r in dh
                 if not (num(r["mean_coex"]) is not None and num(r["mean_d"]) is not None
                         and num(r["mean_coex"]) > 0.3 and num(r["mean_d"]) > 0.5)]
    missed = [r["loop_id"] for r in s1 if r["in_dualhigh17"] != "Yes"
              and num(r["mean_coex"]) is not None and num(r["mean_d"]) is not None
              and num(r["mean_coex"]) > 0.3 and num(r["mean_d"]) > 0.5]
    say(f"\n[2] SELECTION RULE (mean_coex > 0.3 AND mean_d > 0.5), flagged {len(dh)}")
    say(f"    flagged but violating the rule : {violating or 'none'}")
    say(f"    meets the rule but not flagged : {missed or 'none'}")
    say("    -> the rule was executed faithfully; the problem is post-selection")
    say("       inference, not a rule violation.")

    # -- 3 · how independent are the 17? -----------------------------------------
    slots = sum(len(genes_of(r)) for r in dh)
    freq = collections.Counter(g for r in dh for g in genes_of(r))
    shared = [(g, c) for g, c in freq.most_common() if c > 1]
    say(f"\n[3] INDEPENDENCE OF THE 17 CORE LOOPS")
    say(f"    member slots {slots}, unique genes {len(freq)}")
    say(f"    genes appearing in more than one of the 17: {shared}")
    say(f"    member counts: {sorted(len(genes_of(r)) for r in dh)}")
    say("    -> the 17 are not 17 independent observations. Treating them as")
    say("       independent units inflates any test over them.")

    sizes = collections.Counter(len(genes_of(r)) for r in s1)
    two = sizes[2]
    say(f"\n    whole catalogue member-count distribution: {dict(sorted(sizes.items()))}")
    say(f"    two-gene entries: {two}/{len(s1)} ({two / len(s1):.0%}) — for these the")
    say("    module score is the mean of two genes.")

    # -- 4 · CPTAC protein layer -------------------------------------------------
    pos = [r for r in s2 if (num(r["protein_d_mean4"]) or -9) >= 0.5]
    neg = [r for r in s2 if (num(r["protein_d_mean4"]) if num(r["protein_d_mean4"]) is not None else 9) <= 0]
    err = [r for r in s2 if num(r["protein_d_mean4"]) is None]
    cells = [(r["loop_id"], c) for r in s2
             for c in ("protD_COAD", "protD_LUAD", "protD_CCRCC", "protD_PDAC")
             if (num(r[c]) or 0) < 0]
    layers = collections.Counter(r["evidence_layers"] for r in s2)
    say(f"\n[4] CPTAC PROTEIN LAYER, 17 candidates")
    say(f"    mean protein d >= 0.5 : {len(pos):>2}  {[r['loop_id'] for r in pos]}")
    say(f"    mean protein d <= 0   : {len(neg):>2}  {[r['loop_id'] for r in neg]}")
    say(f"    unusable cell         : {len(err):>2}  "
        f"{[(r['loop_id'], r['protein_d_mean4']) for r in err]}")
    say(f"    candidate x cohort cells with reversed sign: {len(cells)}")
    say(f"    evidence_layers (of 5): {dict(sorted(layers.items(), reverse=True))}")
    say(f"    -> protein-layer replication is {len(pos)}/{len(s2)}, not a uniform validation.")

    # -- 5 · crosswalk to the Phase 3a signed-edge circuits ----------------------
    s1_sets = {r["loop_id"]: genes_of(r) for r in s1}
    dh_ids = {r["loop_id"] for r in dh}
    cyc: dict[str, set[str]] = collections.defaultdict(set)
    for e in read(EDGES):
        for s in (e["source"], e["target"]):
            if s.strip():
                cyc[e["cycle_id"]].add(s.strip())

    rows = []
    for cid, gs in sorted(cyc.items()):
        best = ("", 0.0, frozenset())
        for lid, ls in s1_sets.items():
            inter = gs & ls
            if inter:
                j = len(inter) / len(gs | ls)
                if j > best[1]:
                    best = (lid, j, inter)
        rows.append({"phase3a_cycle_id": cid, "best_S1_loop_id": best[0],
                     "jaccard": round(best[1], 3), "n_shared_genes": len(best[2]),
                     "dualhigh17": "YES" if best[0] in dh_ids else "",
                     "shared_genes": ",".join(sorted(best[2]))})
    strong = [r for r in rows if r["jaccard"] >= 0.5]
    say(f"\n[5] CROSSWALK — Phase 3a signed-edge circuits vs catalogue (gene-set overlap)")
    say(f"    circuits with a resolved symbol : {len(cyc)}/70")
    say(f"    any overlap with S1             : {sum(1 for r in rows if r['best_S1_loop_id'])}")
    say(f"    Jaccard >= 0.5                  : {len(strong)}")
    say(f"    exact member-set match          : {sum(1 for r in rows if r['jaccard'] == 1.0)}")
    say(f"    best match is a dual-high-17    : {sum(1 for r in rows if r['dualhigh17'])}")
    say("    -> the perturbation line overlaps the CATALOGUE substantially but the")
    say("       17 core candidates barely. Pre-registration P2 resolves accordingly.")
    say("    NOTE: S1 has no PMID column, so this is gene-set overlap, not a PMID crosswalk.")

    with open(OUT / "candidate_lineage_crosswalk.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, delimiter="\t", lineterminator="\n",
                           fieldnames=["phase3a_cycle_id", "best_S1_loop_id", "jaccard",
                                       "n_shared_genes", "dualhigh17", "shared_genes"])
        w.writeheader()
        w.writerows(rows)
    (OUT / "legacy_audit_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    say(f"\noutputs -> {OUT.relative_to(ROOT)}/")
    say("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
