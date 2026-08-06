# PFL08 — Phase 3: perturbation reciprocity and closure

Code and locked outputs for Phase 3 of the PFL project: testing whether literature-curated
pro-oncogenic positive feedback loops show **closed causal topology** in genome-scale
perturbation data, and whether validated closure — rather than loop membership — predicts
gene dependency.

The observable changed at this phase. Four controlled comparisons failed to show that bulk
expression of loop-member genes carries information beyond proliferation and composition,
and the root cause was not which layer was measured: **abundance is a state variable while a
loop is a causal object.** Co-expression of A and B is equally consistent with a cycle A⇄B,
a feedforward X→A / X→B, and a shared upstream program. Closure has one operational
definition — *perturb A, does B move? perturb B, does A move? do the signs match the curated
polarity?* — so Phase 3 measures perturbation response instead.

Governing documents and the clause → code map: [`docs/PROTOCOL_INDEX.md`](docs/PROTOCOL_INDEX.md).

## Current status

| | |
|---|---|
| Phase 3a offline coverage audit | **done** — [`docs/PHASE3A_OFFLINE_FINDING.md`](docs/PHASE3A_OFFLINE_FINDING.md) |
| Ruling | **`EXPLORATORY_ONLY_ALREADY_DETERMINED`** — gate 1 cannot be met by the current frozen edge set |
| Phase 3a full run (R1/R2/R4/R6) | not run — needs ~550 MB of resources, see below |
| §2.2 independent second review | **not done — manuscript-blocking** |
| Phase 3b | blocked: needs the second review, a calibration pass and its own analysis lock |

The short version: 53 adjudicated records give 229 edges, of which 110 survive the
resource-independent exclusions and at most **95** could ever count as T-visible testable
edges — against a frozen threshold of 150. Gate 1 is unreachable on this edge set no matter
which screen is used. The one permitted remedy (protocol §2.3) is extending signed-edge
adjudication to the 248 un-adjudicated legacy records, which is outcome-blind literature
work and carries no leakage risk.

## Layout

```
data/frozen/     immutable inputs; checksummed and pinned by test
pfl_phase3/      the rules — circuit parsing, node resolution, visibility, gates
scripts/         runnable entry points
results/         locked outputs, including negative rulings (protocol §12.5)
docs/            protocol pointer, amendments, findings
tests/           61 tests over the rules the protocol states
```

## Running it

```bash
pip install pandas numpy pytest      # h5py + requests additionally for the full run
python3 -m pytest tests -q
python3 scripts/run_phase3a_offline.py
```

`run_phase3a_offline.py` needs nothing but this repository — no network, no resources. It
computes the exclusions that hold in every context and turns them into upper bounds on the
frozen gates.

The full run needs the perturbation resources, which are blocked from this sandbox by
network policy and are downloaded in Colab:

```python
from google.colab import drive; drive.mount('/content/drive')
!git clone https://github.com/Conan-trucsafaire/PFL08 /content/PFL08
%cd /content/PFL08
!pip -q install h5py requests pandas numpy
!PFL_INPUTS=/content/drive/MyDrive/PFL08Claude版/00_inputs python scripts/run_phase3a.py
```

Required in `PFL_INPUTS` (download **only** the `*_normalized_bulk_*` files — the
single-cell versions are 65.8 GB and nothing here uses them):

| file | source | size |
|---|---|---|
| `K562_gwps_normalized_bulk_01.h5ad` | Figshare+ DOI `10.25452/figshare.plus.20029387` | ~375 MB |
| `rpe1_normalized_bulk_01.h5ad` | same record | ~95 MB |
| `K562_essential_normalized_bulk_01.h5ad` *(optional)* | same record | ~80 MB |
| `CRISPRGeneEffect.csv` *(optional)* | DepMap Public 26Q1 | — |
| `southard2025_crispra_perturbed_genes.txt` *(optional)* | Zenodo 15213597, one symbol per line | — |

## The disciplines this code enforces

Phase 3a is an **outcome-blind** step: it may read gene identifiers and resource gene lists,
and nothing else. That is enforced by code rather than convention. `pfl_phase3.blind.BlindH5`
raises `BlindnessViolation` on any attempt to reach `/X`, `/layers`, `/obsm`, `/varm`,
`/obsp` in an `.h5ad`; the guard is self-tested at startup by deliberately reaching for `/X`;
DepMap is read with `nrows=0`; every file access is logged to disk for later audit. Do not
substitute `anndata.read_h5ad` — it loads `X`, and the blindness claim in the decision file
becomes false.

Three further rules are load-bearing and easy to get backwards:

- **Edge visibility is a property of the target, not the source.** The readout for A→B is B's
  transcript, so B's mechanism descriptor decides whether the edge can produce an mRNA signal.
  Reversing this silently destroys the T-invisible negative control, which is the
  falsification arm of the design.
- **The classifier never promotes.** `T_visible` is asserted only when a descriptor names a
  transcript-level readout outright; anything ambiguous is `T_uncertain`. Promotion is a
  human act recorded in the review worklist.
- **Gate thresholds are frozen and immutable** (`gates.GATES` is a read-only mapping, and a
  gate absent from the observed counts raises rather than silently passing). When a gate
  fails, the only permitted remedy is extending curation — never lowering a threshold,
  relaxing detectability, substituting inferred L1000 genes, or swapping in
  database-enumerated cycles.

Negative and null results stay in `results/` (protocol §12.5). History is not overwritten:
corrections are recorded as numbered amendments stating whether they were made before or
after which outcomes were read — see [`docs/AMENDMENT_A1_node_resolution_defect.md`](docs/AMENDMENT_A1_node_resolution_defect.md).

## Standing prohibitions

Inherited verbatim from Charter §4 and extended to the abstract and title: no claim, under
any outcome, uses **necessary**, **sufficient**, **irreversible**, **prove**, **universal**,
"hysteresis", or "proves positive feedback causes cancer".
