"""Scope checks (spec section 5, item 9): no biological or dataset references,
only elementary mathematical graphics, and the illustration label in every clip.

These are static checks of the source; the visual review of rendered frames
is done separately.
"""
import re
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from main import SCENES  # noqa: E402

SOURCES = sorted(list((ROOT / "scenes").glob("*.py")) + list((ROOT / "helpers").glob("*.py"))
                 + [ROOT / "main.py"])
SCENE_FILES = [ROOT / "scenes" / f"{m}.py" for m, _, _ in SCENES]

# words that would tie the clips to a specific molecule, sample, structure or dataset
BANNED = [
    r"protein", r"molecul", r"\batom", r"residue", r"amino", r"peptide", r"enzyme",
    r"ribbon", r"helix", r"helices", r"\bstrand", r"\bdomains?\b", r"subunit", r"chain",
    r"phosphat", r"PP2A", r"B56", r"ligand", r"virus", r"ribosom", r"nucleo", r"\bDNA\b",
    r"\bRNA\b", r"EMDB", r"EMPIAR", r"\bPDB\b", r"cryosparc", r"relion", r"angstrom",
    "Å", r"micrograph", r"specimen", r"biolog", r"\bcell\b", r"tissue",
]
# mathematical phrases that contain an otherwise banned word
ALLOWED = [r"(fourier|frequency|real|spatial|image)[ -]domain", r"markov chain"]
# Manim classes that would draw a continuous 3D body or import external assets
FORBIDDEN_CLASSES = [r"\bSurface\(", r"\bSphere\(", r"\bCylinder\(", r"\bTorus\(",
                     r"\bPolyhedron\(", r"\bSVGMobject\(", r"ImageMobject\(\s*['\"]",
                     r"OpenGLSurface", r"\bCone\("]


def _text(p):
    return p.read_text(encoding="utf-8")


@pytest.mark.parametrize("path", SOURCES, ids=lambda p: p.name)
def test_no_biological_or_dataset_vocabulary(path):
    txt = _text(path)
    for a in ALLOWED:
        txt = re.sub(a, "", txt, flags=re.IGNORECASE)
    hits = [w for w in BANNED if re.search(w, txt, flags=re.IGNORECASE)]
    assert not hits, f"{path.name}: {hits}"


@pytest.mark.parametrize("path", SOURCES, ids=lambda p: p.name)
def test_only_elementary_mathematical_graphics(path):
    txt = _text(path)
    hits = [c for c in FORBIDDEN_CLASSES if re.search(c, txt)]
    assert not hits, f"{path.name}: {hits}"


@pytest.mark.parametrize("path", SCENE_FILES, ids=lambda p: p.name)
def test_every_clip_carries_the_illustration_label(path):
    if not path.exists():
        pytest.skip("clip not written yet")
    txt = _text(path)
    assert "illustration_label(" in txt
    assert "title_card(" in txt and "takeaway_card(" in txt


def test_output_names_match_spec():
    assert [o for _, _, o in SCENES] == [
        "01_forward_model.mp4", "02_alignment_averaging.mp4", "03_orientation_fourier_slice.mp4",
        "04_iterative_refinement.mp4", "05_classification_pca.mp4", "06_fsc_resolution.mp4"]
