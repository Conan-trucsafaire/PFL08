"""The frozen candidate table must stay parseable and unaltered.

data/frozen/ is an immutable input (protocol §12.4). These tests fail loudly if it is
edited, and they pin the structural facts the offline audit reports.
"""

import csv
from pathlib import Path

from pfl_phase3.circuits import parse_circuit_table
from pfl_phase3.graph import build_edges_and_cycles
from pfl_phase3.nodes import resolve_node

TSV = (Path(__file__).resolve().parents[1] / "data" / "frozen"
       / "02_SIGNED_EDGE_FIRST_PASS_CANDIDATES_FROZEN_2026-08-05.tsv")

EXPECTED_SHA256 = "3E9E9CFB1E2922A2D94C5AAC83E7619199E377FD0E4C4932334DDCD5165452D6"


def _records():
    with open(TSV, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def test_frozen_table_is_unmodified():
    import hashlib
    assert hashlib.sha256(TSV.read_bytes()).hexdigest().upper() == EXPECTED_SHA256


def test_record_count_and_columns():
    recs = _records()
    assert len(recs) == 53
    assert set(recs[0]) == {"pmid", "status", "second_review_pending", "cancer",
                            "compartment", "circuit"}


def test_every_record_is_still_pending_second_review():
    """Protocol §2.2 — manuscript-blocking until an independent reviewer has adjudicated."""
    assert all(r["second_review_pending"] == "TRUE" for r in _records())


def test_every_circuit_parses():
    circuits = parse_circuit_table(_records())
    unparsed = [c for c in circuits if c["parse_status"] != "OK"]
    assert unparsed == [], f"unparsed circuits: {[c['raw'] for c in unparsed]}"
    assert len(circuits) == 70


def test_exactly_one_curated_loop_has_a_negative_sign_product():
    """A curated *positive* loop should multiply to +1.

    42338478 reads `STRAP ┤ DKK1-LRP6 inhibition → β-catenin/TCF4 → STRAP`: the middle
    node is itself named as an inhibition, so the arrow signs multiply to −1 while the
    biology is a positive loop. This is a node-naming artefact, not a curation error,
    and it is exactly the kind of item the §2.2 second review must rule on. The test
    pins it so a second one cannot appear unnoticed.
    """
    circuits = parse_circuit_table(_records())
    negative = [c["cycle_id"] for c in circuits if c["sign_product"] != +1]
    assert negative == ["42338478#0"]


def test_edge_and_cycle_counts_are_stable():
    circuits = parse_circuit_table(_records())
    resolved = {(c["cycle_id"], p): resolve_node(tok, None)
                for c in circuits for p, tok in enumerate(c["nodes"])}
    edges, cycles = build_edges_and_cycles(circuits, resolved)
    assert len(cycles) == 70
    assert len(edges) == 229
