"""E2b calibration — and the reason it cannot serve as a topology endpoint.

A module score over k genes is a function of exactly two things: the members' marginal
effects, and their correlation matrix. Its advantage over the best single member is an
averaging gain of roughly sqrt(k / (1 + (k-1)rho)). That gain is *largest* when members
are independent and shrinks to nothing as they become co-regulated.

So the increment runs backwards for the purpose it was proposed for: a genuinely
co-regulated set scores below an effect-matched random set, and a biologically
meaningless set of independent genes scores above it. These tests pin that behaviour so
the endpoint cannot be quietly reintroduced as evidence of coordination or closure.

What the statistic *can* do is answer whether a candidate's signal is reducible to one
member. That is a redundancy question, and it is the one the tests below endorse.
"""

import numpy as np
import pytest

from pfl_phase3.increment import (cohen_d, leave_one_out, increment_vs_null, loop_increment,
                                  member_effects, module_effect)

N_T, N_N = 120, 80
IS_T = np.r_[np.ones(N_T, bool), np.zeros(N_N, bool)]


def make_set(k, d, rho, seed):
    """k genes each with marginal effect d and exchangeable within-set correlation rho."""
    rng = np.random.default_rng(seed)
    shared = rng.normal(size=N_T + N_N)
    rows = []
    for _ in range(k):
        e = rng.normal(size=N_T + N_N)
        rows.append(np.sqrt(rho) * shared + np.sqrt(1 - rho) * e + d * IS_T)
    return np.array(rows)


def background(seed=7, n=400):
    rng = np.random.default_rng(seed)
    return np.array([rng.normal(size=N_T + N_N) + d * IS_T
                     for d in rng.uniform(-0.2, 1.0, n)])


# ----------------------------------------------------------------- basics

def test_cohen_d_sign_and_scale():
    rng = np.random.default_rng(0)
    a = rng.normal(size=100)
    assert cohen_d(a + 1, a) == pytest.approx(1.0, abs=0.25)
    assert cohen_d(a, a + 1) == pytest.approx(-1.0, abs=0.25)
    assert cohen_d(np.array([1.0]), np.array([0.0])) == 0.0, "n<2 is undefined, not a crash"
    assert cohen_d(np.ones(50), np.ones(50)) == 0.0, "zero pooled SD is undefined, not inf"


def test_member_and_module_effects_recover_the_signal():
    e = make_set(3, 0.8, 0.0, seed=1)
    assert member_effects(e, IS_T).mean() == pytest.approx(0.8, abs=0.25)
    assert module_effect(e, IS_T) > 0.8


def test_increment_needs_at_least_two_members():
    with pytest.raises(ValueError):
        loop_increment(np.random.default_rng(0).normal(size=(1, 200)), IS_T)


# ----------------------------------------------------------------- the flaw, pinned

@pytest.mark.parametrize("rho", [0.0, 0.3, 0.6, 0.9])
def test_module_advantage_shrinks_as_members_become_coregulated(rho):
    """Averaging gain falls with correlation. This is arithmetic, not biology."""
    r = loop_increment(make_set(3, 0.5, rho, seed=1), IS_T)
    ceiling = {0.0: (0.15, 0.45), 0.3: (0.02, 0.30), 0.6: (-0.05, 0.20), 0.9: (-0.15, 0.10)}[rho]
    assert ceiling[0] <= r.increment <= ceiling[1]


def test_increment_is_monotone_decreasing_in_correlation():
    inc = [loop_increment(make_set(3, 0.5, rho, seed=1), IS_T).increment
           for rho in (0.0, 0.3, 0.6, 0.9)]
    assert inc == sorted(inc, reverse=True)


def test_effect_matched_null_ranks_a_coregulated_set_BELOW_an_independent_one():
    """The endpoint runs backwards — the central reason E2b was withdrawn.

    An independent set carries no coordination at all and beats the null; a co-regulated
    set, which is what a real circuit would look like, falls to the bottom of it.
    """
    bg = background()
    indep = increment_vs_null(make_set(3, 0.5, 0.0, seed=2), bg, IS_T,
                           n_draws=400, rng=np.random.default_rng(3))
    coreg = increment_vs_null(make_set(3, 0.5, 0.6, seed=2), bg, IS_T,
                           n_draws=400, rng=np.random.default_rng(3))
    assert indep.percentile > 60.0
    assert coreg.percentile < 20.0
    assert coreg.percentile < indep.percentile


# ----------------------------------------------------------------- what survives

def test_single_carrier_set_is_detected_as_reducible():
    """One member carries the signal, the rest are noise: the set adds nothing."""
    rng = np.random.default_rng(11)
    carrier = rng.normal(size=N_T + N_N) + 1.5 * IS_T
    noise = [rng.normal(size=N_T + N_N) for _ in range(2)]
    r = loop_increment(np.array([carrier, *noise]), IS_T)
    assert r.d_best_member > 1.2
    assert r.increment < 0, "a set carried by one gene must not beat that gene"


def test_leave_one_out_identifies_the_carrier():
    rng = np.random.default_rng(12)
    carrier = rng.normal(size=N_T + N_N) + 1.5 * IS_T
    noise = [rng.normal(size=N_T + N_N) for _ in range(2)]
    loo = leave_one_out(np.array([carrier, *noise]), IS_T)
    assert len(loo) == 3
    assert loo[0] == min(loo), "dropping the carrier must cost the most"


def test_leave_one_out_undefined_for_two_member_sets():
    assert leave_one_out(make_set(2, 0.5, 0.0, seed=4), IS_T).size == 0


def test_null_pool_exhaustion_is_reported_not_hidden():
    tiny = background(seed=9, n=12)
    r = increment_vs_null(make_set(3, 0.5, 0.0, seed=5), tiny, IS_T,
                       n_draws=50, rng=np.random.default_rng(6))
    assert r.null_pool_exhausted is True
