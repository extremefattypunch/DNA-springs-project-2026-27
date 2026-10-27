"""Smoke-render the final frame of every clip at low quality (-ql -s).

Slow (about a minute per clip); deselect with:  pytest -m "not render"
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from main import SCENES  # noqa: E402


@pytest.mark.render
@pytest.mark.parametrize("module,cls,outname", SCENES, ids=[s[1] for s in SCENES])
def test_scene_renders_last_frame(module, cls, outname, tmp_path):
    cmd = [sys.executable, "-m", "manim", "render", "-ql", "-s", "--media_dir", str(tmp_path),
           str(ROOT / "scenes" / f"{module}.py"), cls]
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True,
                       env={**os.environ, "PYTHONPATH": str(ROOT)}, timeout=1800)
    assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-3000:]
    pngs = list(tmp_path.rglob(f"{cls}*.png"))
    assert pngs, "no still image written"
    assert pngs[0].stat().st_size > 10_000
