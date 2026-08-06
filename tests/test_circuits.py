"""Arrow parsing. All fixtures are verbatim strings from the frozen candidate table."""

import pytest

from pfl_phase3.circuits import (CircuitParseError, parse_circuit_table, sign_product,
                                 split_circuit, split_record)


@pytest.mark.parametrize("text,n_nodes,signs,product", [
    # ┤ and → in one string; two inhibitions -> positive loop
    ("HMGA1 → TRIP13 ┤ FBXW7 ┤ MYC → HMGA1", 5, [+1, -1, -1, +1], +1),
    # -| and -> mixed
    ("YTHDC1 activity -| SLC2A3/GLUT3 abundance -> RNF183 abundance -| YTHDC1 protein stability",
     4, [-1, +1, -1], +1),
    ("GATA4 → CDX2 → GATA4", 3, [+1, +1], +1),
    ("TTC13 -> JAK2/STAT3 activity -> TTC13", 3, [+1, +1], +1),
    ("SNHG17 ┤ miR-339-5p ┤ FOSL2 → SNHG17", 4, [-1, -1, +1], +1),
    # markdown-escaped arrow form, which appears in documents exported from Drive
    ("A -\\> B -\\> A", 3, [+1, +1], +1),
])
def test_split_circuit(text, n_nodes, signs, product):
    nodes, got = split_circuit(text)
    assert len(nodes) == n_nodes
    assert got == signs
    assert sign_product(got) == product


def test_nodes_are_stripped():
    nodes, _ = split_circuit("  A   ->   B  ")
    assert nodes == ["A", "B"]


def test_no_arrow_raises():
    with pytest.raises(CircuitParseError):
        split_circuit("MYC only")


def test_semicolon_splits_one_record_into_several_circuits():
    cell = ("MYC -> METTL3 -| MYC; "
            "MYC -> METTL3 -> LINC01006 -| miR-34a activity -| MYC")
    assert len(split_record(cell)) == 2


def test_parse_circuit_table_assigns_stable_cycle_ids():
    rows = parse_circuit_table([{"pmid": "123", "circuit": "A -> B -> A; C -> D -> C"}])
    assert [r["cycle_id"] for r in rows] == ["123#0", "123#1"]
    assert all(r["parse_status"] == "OK" for r in rows)


def test_unparsed_circuit_is_recorded_not_dropped():
    rows = parse_circuit_table([{"pmid": "9", "circuit": "no arrows here"}])
    assert len(rows) == 1
    assert rows[0]["parse_status"] == "UNPARSED"
    assert rows[0]["parse_error"]


def test_direction_is_never_inferred_from_token_order():
    """Protocol §2.1 — the same nodes with reversed arrows give reversed edges."""
    fwd, _ = split_circuit("A -> B")
    rev, _ = split_circuit("B -> A")
    assert fwd == ["A", "B"]
    assert rev == ["B", "A"]
