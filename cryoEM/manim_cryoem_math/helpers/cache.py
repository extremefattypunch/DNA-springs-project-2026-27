"""Tiny .npz cache so scenes never recompute FFTs/projections at render time.

Cache directory: $MANIM_CRYOEM_CACHE (set by env/activate_manim.sh, on
netscratch), else ~/.cache/manim_cryoem_math. Delete it to force a rebuild;
everything is regenerable from fixed seeds.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np

CACHE_VERSION = 1


def cache_dir() -> Path:
    d = Path(os.environ.get("MANIM_CRYOEM_CACHE",
                            Path.home() / ".cache" / "manim_cryoem_math"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def cached(name: str, fn, **params) -> dict:
    """Return fn(**params) (a dict of arrays/scalars), cached on disk by name+params."""
    key = hashlib.sha1(json.dumps([CACHE_VERSION, params], sort_keys=True,
                                  default=str).encode()).hexdigest()[:12]
    path = cache_dir() / f"{name}-{key}.npz"
    if path.exists():
        with np.load(path, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}
    out = fn(**params)
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **{k: np.asarray(v) for k, v in out.items()})
    tmp.replace(path)
    return {k: np.asarray(v) for k, v in out.items()}
