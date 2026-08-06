"""Node resolution and the mRNA-visibility declaration (protocol §3.2)."""

import pytest

from pfl_phase3.nodes import classify_visibility, resolve_node

UNIVERSE = {"MYC", "WNT7B", "AKT1", "MAPK14", "CTNNB1", "TCF7L2", "CSF3R", "DKK1",
            "LRP6", "STRAP", "RNF183", "SLC2A3", "YTHDC1", "CDCA3", "TRIM28"}


# ------------------------------------------------------------------ visibility

@pytest.mark.parametrize("token,expected", [
    ("RRAGB transcription/abundance", "T_visible"),
    ("CSF3R mRNA stability", "T_visible"),
    ("HEGBC transcription", "T_visible"),
    ("MYC protein stability", "T_invisible"),
    # a descriptor naming both layers stays uncertain — the code resolves neither way
    ("STAT3 Cys108 palmitoylation/transcriptional activity", "T_uncertain"),
    ("YAP1 nuclear/transcriptional activity", "T_uncertain"),
    # ambiguous descriptors are never promoted by code
    ("HIF1A activity", "T_uncertain"),
    ("RNF183 abundance", "T_uncertain"),
    ("MYC", "T_uncertain"),
    # both families present -> uncertain, not silently resolved either way
    ("KDM5B mRNA stability/protein abundance", "T_uncertain"),
])
def test_classify_visibility(token, expected):
    assert classify_visibility(token) == expected


def test_classifier_never_promotes_to_visible_on_ambiguity():
    """The human reviewer promotes; the code must not (protocol §3.2)."""
    for token in ["SOX2 abundance/activity", "EGFR activity", "GSK3B stability"]:
        assert classify_visibility(token) != "T_visible"


# ------------------------------------------------------------------ node classes

def test_tme_node_is_excluded():
    r = resolve_node("macrophage IL-1beta", UNIVERSE)
    assert r["node_class"] == "cell_type_prefixed"
    assert r["symbol"] is None


def test_noncoding_node_is_excluded():
    for token in ["miR-452-5p", "circPVT1", "LINC01605", "HCG18", "CSF3R-AS", "MIR17HG"]:
        assert resolve_node(token, UNIVERSE)["node_class"] == "noncoding_rna"


def test_pathway_token_gets_no_gene_assigned():
    for token in ["mTORC1 activity", "MAPK/ERK activity", "AKT activity"]:
        r = resolve_node(token, UNIVERSE)
        assert r["node_class"] == "complex_or_pathway"
        assert r["symbol"] is None


def test_process_node_is_excluded():
    assert resolve_node("glycolysis", UNIVERSE)["node_class"] == "metabolite_or_process"


def test_unmatched_symbol_is_unmapped_not_guessed():
    r = resolve_node("NOTAGENE1 activity", UNIVERSE)
    assert r["node_class"] == "unmapped"
    assert r["symbol"] is None


# ------------------------------------------------------------------ regressions

def test_gene_starting_with_a_pathway_name_is_not_a_pathway():
    """Amendment A1. Prefix matching swallowed WNT7B ('Wnt'), AKT1 ('AKT'), MAPK14 ('MAPK')."""
    for token, symbol in [("WNT7B", "WNT7B"), ("secreted WNT7B", "WNT7B"), ("AKT1", "AKT1")]:
        r = resolve_node(token, UNIVERSE)
        assert r["node_class"].startswith("single_symbol"), token
        assert r["symbol"] == symbol


def test_pathway_matching_still_fires_on_whole_tokens():
    assert resolve_node("AKT activity", UNIVERSE)["node_class"] == "complex_or_pathway"
    assert resolve_node("Wnt signaling", UNIVERSE)["node_class"] == "complex_or_pathway"


def test_trailing_hyphen_is_not_part_of_a_symbol():
    """Amendment A1. `IL-1beta` used to yield the symbol `IL-`, `DKK1-LRP6` yielded `DKK1-`."""
    r = resolve_node("DKK1-LRP6 inhibition", UNIVERSE)
    assert r["symbol"] == "DKK1"
    assert not r["symbol"].endswith("-")


# ------------------------------------------------------------------ permissive mode

def test_permissive_mode_is_an_upper_bound_on_strict_mode():
    """Permissive resolution may accept more nodes than strict, never fewer."""
    token = "NOTAGENE1 activity"
    assert resolve_node(token, UNIVERSE)["symbol"] is None
    assert resolve_node(token, None)["symbol"] == "NOTAGENE1"


def test_permissive_mode_flags_everything_it_provisionally_maps():
    r = resolve_node("NOTAGENE1 activity", None)
    assert r["review_flag"] is True
    assert r["node_class"].endswith("_unverified")


def test_permissive_mode_does_not_weaken_structural_exclusions():
    for token in ["miR-452-5p", "mTORC1 activity", "macrophage IL-1beta", "glycolysis"]:
        assert resolve_node(token, None)["symbol"] is None
