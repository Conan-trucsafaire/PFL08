"""Decomposition of a curated member set into its members.

This module was written to implement a proposed endpoint — "does a member set beat its
own best member, by more than an effect-matched random set would?" — intended as
observational evidence that a circuit does something *as a circuit*. Implementing it
showed that it cannot serve that purpose, and the demonstration is pinned in
tests/test_increment.py.

**Why it fails as a coordination or topology endpoint.** A module score over k genes is
a function of exactly two quantities: the members' marginal effects and their correlation
matrix. Its advantage over the best single member is an averaging gain of roughly
sqrt(k / (1 + (k-1)rho)). That gain is *maximal* when the members are independent and
vanishes as they become co-regulated. So the statistic runs backwards: a co-regulated set
— what a real circuit should look like — scores near the bottom of the effect-matched
null, while a set of unrelated genes scores near the top. No choice of null repairs this,
because the ordering is arithmetic rather than statistical.

The same argument disposes of the wider hope: since a module score is determined by
marginals plus correlation, and both are properties of *membership*, no module statistic
over member expression can distinguish a closed cycle from an open chain on the same
genes. Topology is not recoverable from this observable at all.

**What the code is still for.** Asking whether a candidate's tumour signal is reducible
to a single member — a redundancy question, not a topology one. That is directly relevant
to candidates such as PSMD14+LDHA and TTK+BUB1B+LDHA, where one member may carry
everything (pre-registration amendment A1, finding F4). `loop_increment` and
`leave_one_out` answer it; `increment_vs_null` is retained only to reproduce the
calibration result above and must not be reported as evidence of coordination.

Nothing here reads a perturbation dataset. The unit of analysis is the deduplicated
member set (amendment A1, finding F2), not the catalogue entry.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def cohen_d(a: np.ndarray, b: np.ndarray) -> float:
    """Standardised mean difference with a pooled SD. Returns 0.0 if undefined."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return 0.0
    va, vb = a.var(ddof=1), b.var(ddof=1)
    pooled = ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
    if pooled <= 0:
        return 0.0
    return float((a.mean() - b.mean()) / np.sqrt(pooled))


def zscore_rows(x: np.ndarray) -> np.ndarray:
    """Z-score each gene across all samples, as the v5 module score did."""
    x = np.asarray(x, dtype=float)
    mu = x.mean(axis=1, keepdims=True)
    sd = x.std(axis=1, ddof=1, keepdims=True)
    sd = np.where(sd <= 0, np.nan, sd)
    return (x - mu) / sd


def member_effects(expr: np.ndarray, is_tumour: np.ndarray) -> np.ndarray:
    """Per-gene tumour-vs-normal Cohen's d, one value per row of `expr`."""
    is_tumour = np.asarray(is_tumour, dtype=bool)
    return np.array([cohen_d(row[is_tumour], row[~is_tumour]) for row in np.asarray(expr, float)])


def module_effect(expr: np.ndarray, is_tumour: np.ndarray) -> float:
    """Tumour-vs-normal Cohen's d of the mean z-scored member profile."""
    z = zscore_rows(expr)
    score = np.nanmean(z, axis=0)
    is_tumour = np.asarray(is_tumour, dtype=bool)
    return cohen_d(score[is_tumour], score[~is_tumour])


@dataclass
class IncrementResult:
    """Outcome of E2b for one candidate in one cohort."""

    n_members: int
    d_module: float
    d_members: np.ndarray
    d_best_member: float
    increment: float
    null_increments: np.ndarray = field(default_factory=lambda: np.array([]))
    percentile: float | None = None
    n_null_draws: int = 0
    null_pool_exhausted: bool = False

    @property
    def beats_best_member(self) -> bool:
        return self.increment > 0

    def exceeds_null(self, alpha: float = 0.05) -> bool:
        """Increment above the (1 - alpha) quantile of the effect-matched null."""
        if self.percentile is None:
            return False
        return self.percentile >= (1.0 - alpha) * 100.0


def loop_increment(expr: np.ndarray, is_tumour: np.ndarray) -> IncrementResult:
    """Increment of a member set over its own best member, before any null comparison."""
    expr = np.asarray(expr, dtype=float)
    if expr.ndim != 2 or expr.shape[0] < 2:
        raise ValueError("a member set needs at least two measured members")
    d_mem = member_effects(expr, is_tumour)
    d_mod = module_effect(expr, is_tumour)
    best = float(np.max(d_mem))
    return IncrementResult(n_members=expr.shape[0], d_module=d_mod, d_members=d_mem,
                           d_best_member=best, increment=float(d_mod - best))


def effect_matched_null(background: np.ndarray, is_tumour: np.ndarray,
                        target_effects: np.ndarray, n_draws: int = 1000,
                        tol: float = 0.1, rng: np.random.Generator | None = None,
                        exclude: set[int] | None = None) -> tuple[np.ndarray, bool]:
    """N8 — increments of random sets matched to `target_effects` gene by gene.

    For each member of the candidate, a background gene with a similar individual effect
    is drawn, and the resulting synthetic set is put through the identical statistic. The
    null therefore carries the same averaging gain and the same max-selection bias as the
    candidate, and differs from it only in whether the members are biologically coupled.

    Returns (increments, pool_exhausted). `pool_exhausted` is True when any member had
    fewer than two eligible partners, which is reported rather than silently tolerated.
    """
    rng = rng or np.random.default_rng(42)
    background = np.asarray(background, dtype=float)
    exclude = exclude or set()
    bg_d = member_effects(background, is_tumour)

    pools = []
    exhausted = False
    for t in np.asarray(target_effects, dtype=float):
        idx = np.flatnonzero(np.abs(bg_d - t) <= tol)
        idx = np.array([i for i in idx if i not in exclude], dtype=int)
        if len(idx) < 2:
            exhausted = True
            widened = np.argsort(np.abs(bg_d - t))
            idx = np.array([i for i in widened if i not in exclude][:50], dtype=int)
        pools.append(idx)

    out = np.empty(n_draws, dtype=float)
    for j in range(n_draws):
        pick = [int(rng.choice(p)) for p in pools]
        out[j] = loop_increment(background[pick], is_tumour).increment
    return out, exhausted


def increment_vs_null(expr: np.ndarray, background: np.ndarray, is_tumour: np.ndarray,
                   n_draws: int = 1000, tol: float = 0.1,
                   rng: np.random.Generator | None = None,
                   exclude: set[int] | None = None) -> IncrementResult:
    """Full E2b for one candidate: increment plus its effect-matched null percentile."""
    res = loop_increment(expr, is_tumour)
    null, exhausted = effect_matched_null(background, is_tumour, res.d_members,
                                          n_draws=n_draws, tol=tol, rng=rng, exclude=exclude)
    res.null_increments = null
    res.n_null_draws = n_draws
    res.null_pool_exhausted = exhausted
    res.percentile = float((null < res.increment).mean() * 100.0)
    return res


def leave_one_out(expr: np.ndarray, is_tumour: np.ndarray) -> np.ndarray:
    """Module effect with each member dropped in turn — exposes single-member carriage."""
    expr = np.asarray(expr, dtype=float)
    if expr.shape[0] < 3:
        return np.array([])
    return np.array([module_effect(np.delete(expr, i, axis=0), is_tumour)
                     for i in range(expr.shape[0])])
