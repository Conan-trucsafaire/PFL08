"""Edge/cycle assembly, testability, blindness and the frozen gates."""

import pytest

from pfl_phase3.blind import BlindnessViolation, forbidden_top_level
from pfl_phase3.circuits import parse_circuit_table
from pfl_phase3.gates import GATES, evaluate_gates, upper_bound_ruling
from pfl_phase3.graph import build_edges_and_cycles, sym
from pfl_phase3.nodes import resolve_node
from pfl_phase3.testability import confirmatory_eligible, structural_reason
from pfl_phase3.testability import testable_in as is_testable_in

UNIVERSE = {"GA", "GB", "GC", "MYC", "CDCA3", "TRIM28"}


def _pipeline(circuit_text, universe=UNIVERSE):
    circuits = parse_circuit_table([{"pmid": "1", "circuit": circuit_text}])
    resolved = {(c["cycle_id"], p): resolve_node(tok, universe)
                for c in circuits for p, tok in enumerate(c["nodes"])}
    return build_edges_and_cycles(circuits, resolved)


# ------------------------------------------------------------------ sym()

def test_sym_treats_nan_as_missing():
    """pandas turns None into NaN, and NaN is truthy and not-None."""
    assert sym(float("nan")) is None
    assert sym(None) is None
    assert sym("nan") is None
    assert sym("  ") is None
    assert sym("  MYC ") == "MYC"


# ------------------------------------------------------------------ edges/cycles

def test_edge_visibility_comes_from_the_target_not_the_source():
    """Protocol §3.2 — the readout for A->B is B's transcript."""
    edges, _ = _pipeline("GA transcription -> GB protein stability")
    assert edges[0]["visibility"] == "T_invisible"
    edges, _ = _pipeline("GA protein stability -> GB transcription")
    assert edges[0]["visibility"] == "T_visible"


def test_cycle_closes_when_first_and_last_resolve_to_the_same_gene():
    _, cycles = _pipeline("CDCA3 -> TRIM28 protective activity toward MYC -> CDCA3")
    assert cycles[0]["closure"] == "CLOSED"


def test_cycle_with_unresolved_endpoint_is_never_called_closed():
    _, cycles = _pipeline("miR-1 -> GA -> miR-1")
    assert cycles[0]["closure"] == "NOT_CLOSED"


def test_edge_count_is_one_less_than_node_count():
    edges, cycles = _pipeline("GA -> GB -> GC -> GA")
    assert len(edges) == 3
    assert cycles[0]["length"] == 3


# ------------------------------------------------------------------ testability

def test_structural_reason_is_resource_independent():
    edges, _ = _pipeline("macrophage IL-1beta -> GA -> macrophage IL-1beta")
    assert structural_reason(edges[0]) == "intercellular"


def test_unmapped_endpoint_is_reported_specifically():
    edges, _ = _pipeline("NOTAGENE9 -> GA")
    assert structural_reason(edges[0]) == "unmapped_endpoint"


def test_cis_proximity_excludes_the_edge():
    edges, _ = _pipeline("GA -> GB")
    edges[0]["cis_conflict"] = True
    assert structural_reason(edges[0]) == "cis_proximity"


def test_testable_requires_source_perturbed_and_target_measured():
    edges, _ = _pipeline("GA -> GB")
    ctx = {"perturbed": {"GA"}, "measured": {"GB"}}
    assert is_testable_in(ctx, edges[0]) == (True, "testable")
    assert is_testable_in({"perturbed": set(), "measured": {"GB"}}, edges[0]) == (False, "source_not_perturbed")
    assert is_testable_in({"perturbed": {"GA"}, "measured": set()}, edges[0]) == (False, "target_not_measured")


def test_only_t_visible_edges_are_confirmatory_eligible():
    edges, _ = _pipeline("GA -> GB transcription")
    assert confirmatory_eligible(edges[0]) is True
    edges, _ = _pipeline("GA -> GB protein stability")
    assert confirmatory_eligible(edges[0]) is False


def test_t_invisible_edges_are_retained_as_the_negative_control():
    """Protocol §3.2/E8 — T-invisible edges are a control, not attrition."""
    edges, _ = _pipeline("GA -> GB protein stability")
    assert structural_reason(edges[0]) is None      # still structurally testable
    assert edges[0]["visibility"] == "T_invisible"  # just not confirmatory


# ------------------------------------------------------------------ gates

def test_gate_thresholds_are_the_frozen_ones():
    assert dict(GATES) == {
        "t_visible_testable_edges_pooled": 150,
        "completely_testable_cycles": 25,
        "contexts_with_50plus_testable_edges": 2,
        "nontransformed_contexts": 1,
        "cycles_testable_in_2plus_contexts": 10,
    }


def test_gates_cannot_be_mutated():
    with pytest.raises(TypeError):
        GATES["t_visible_testable_edges_pooled"] = 1


def test_a_missing_gate_raises_rather_than_silently_passing():
    with pytest.raises(KeyError):
        evaluate_gates({"t_visible_testable_edges_pooled": 999})


def test_all_gates_met_permits_confirmatory():
    _, ruling = evaluate_gates({k: v for k, v in GATES.items()})
    assert ruling == "CONFIRMATORY_PERMITTED"


def test_one_failed_gate_downgrades_the_whole_run():
    observed = dict(GATES)
    observed["completely_testable_cycles"] = 0
    _, ruling = evaluate_gates(observed)
    assert ruling == "EXPLORATORY_ONLY"


def test_upper_bound_below_threshold_settles_the_gate():
    rows, ruling = upper_bound_ruling({"t_visible_testable_edges_pooled": 95})
    verdicts = {r["gate"]: r["verdict"] for r in rows}
    assert verdicts["t_visible_testable_edges_pooled"] == "CANNOT_PASS"
    assert verdicts["completely_testable_cycles"] == "NOT_BOUNDED_OFFLINE"
    assert ruling == "EXPLORATORY_ONLY_ALREADY_DETERMINED"


def test_upper_bound_above_threshold_settles_nothing():
    rows, ruling = upper_bound_ruling({k: v + 1 for k, v in GATES.items()})
    assert all(r["verdict"] == "UNDECIDED" for r in rows)
    assert ruling == "UNDETERMINED_PENDING_RESOURCES"


# ------------------------------------------------------------------ blindness

def test_matrix_groups_are_forbidden():
    for key in ["X", "/X", "layers/counts", "obsm/X_pca", "varm/PCs", "obsp/dist", "raw/X"]:
        assert forbidden_top_level(key)


def test_annotation_groups_are_allowed():
    for key in ["obs", "var", "var/gene_name", "obs/gene", "uns"]:
        assert not forbidden_top_level(key)


def test_blindness_violation_is_an_error_not_a_warning():
    assert issubclass(BlindnessViolation, RuntimeError)
