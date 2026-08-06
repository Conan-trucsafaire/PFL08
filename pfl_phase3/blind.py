"""Outcome-blindness, enforced by code rather than by convention.

Phase 3a uses gene identifiers and resource gene lists only. No effect value, z-score
or p-value may be read (protocol §3). `BlindH5` opens an .h5ad and can reach /obs and
/var and nothing else; any attempt to touch a matrix raises.

Do NOT replace this with ``anndata.read_h5ad`` — that loads X, and the blindness claim
in the gate decision file becomes false.
"""

from __future__ import annotations

import datetime
from pathlib import Path

FORBIDDEN = ("X", "layers", "obsm", "varm", "obsp", "varp", "raw")


class BlindnessViolation(RuntimeError):
    """Raised when outcome-blind code reaches for an effect value."""


class AccessLog(list):
    """Append-only record of every file access, written to disk for later audit."""

    def record(self, path, what: str) -> None:
        self.append({
            "ts": datetime.datetime.now().isoformat(timespec="seconds"),
            "path": str(path),
            "read": what,
        })


def forbidden_top_level(key: str) -> bool:
    """Is this h5ad key inside a matrix group?"""
    return str(key).lstrip("/").split("/")[0] in FORBIDDEN


class BlindH5:
    """Read-only .h5ad accessor limited to /obs and /var."""

    def __init__(self, path, access_log: AccessLog | None = None):
        self.path = Path(path)
        self.access_log = access_log
        self._f = None

    # -- lifecycle ---------------------------------------------------------------
    def __enter__(self):
        import h5py
        self._f = h5py.File(self.path, "r")
        return self

    def __exit__(self, *exc):
        if self._f is not None:
            self._f.close()
            self._f = None

    # -- guarded access ----------------------------------------------------------
    def _check(self, key) -> None:
        if forbidden_top_level(key):
            top = str(key).lstrip("/").split("/")[0]
            raise BlindnessViolation(
                f"Phase 3a may not read '/{top}' in {self.path.name} "
                f"— outcome-blind step (protocol §3)")

    def get(self, key):
        self._check(key)
        if self.access_log is not None:
            self.access_log.record(self.path, key)
        return self._f[key]

    def top_level(self):
        return list(self._f.keys())

    # -- decoding ----------------------------------------------------------------
    @staticmethod
    def _decode(arr):
        import numpy as np
        import pandas as pd
        arr = np.asarray(arr)
        if arr.dtype.kind in ("S", "O"):
            return pd.Index([x.decode() if isinstance(x, bytes) else str(x) for x in arr])
        return pd.Index(arr)

    def index_of(self, group: str):
        g = self.get(group)
        name = g.attrs.get("_index", b"_index")
        if isinstance(name, bytes):
            name = name.decode()
        return self._decode(g[name][:])

    def column(self, group: str, col: str):
        import h5py
        import numpy as np
        import pandas as pd
        g = self.get(f"{group}/{col}")
        if isinstance(g, h5py.Group) and "categories" in g and "codes" in g:
            cats = self._decode(g["categories"][:])
            codes = np.asarray(g["codes"][:])
            return pd.Index([cats[c] if c >= 0 else None for c in codes])
        return self._decode(g[:])

    def columns(self, group: str):
        return [k for k in self.get(group).keys() if not k.startswith("_")]


def self_test(path, access_log: AccessLog | None = None) -> None:
    """Prove the guard fires before trusting any run that claims to be blind."""
    with BlindH5(path, access_log) as h:
        try:
            h.get("X")
        except BlindnessViolation:
            return
    raise AssertionError(f"GUARD FAILED — /X was reachable in {path}")
