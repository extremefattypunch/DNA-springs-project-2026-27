"""Entry point for the six clips.

    manim -qh main.py Scene03OrientationFourierSlice   # one clip, 1080p30
    python main.py --list                             # module, class and output name per clip
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

# (module in scenes/, Scene class, output file name in output/)
SCENES = [
    ("s01_forward_model", "Scene01ForwardModel", "01_forward_model.mp4"),
    ("s02_alignment_averaging", "Scene02AlignmentAveraging", "02_alignment_averaging.mp4"),
    ("s03_orientation_fourier_slice", "Scene03OrientationFourierSlice",
     "03_orientation_fourier_slice.mp4"),
    ("s04_iterative_refinement", "Scene04IterativeRefinement", "04_iterative_refinement.mp4"),
    ("s05_classification_pca", "Scene05ClassificationPCA", "05_classification_pca.mp4"),
    ("s06_fsc_resolution", "Scene06FSCResolution", "06_fsc_resolution.mp4"),
]

if __name__ == "__main__" and "--list" in sys.argv:
    for m, c, o in SCENES:
        print(m, c, o)
else:
    import importlib

    for _m, _c, _o in SCENES:
        try:
            globals()[_c] = getattr(importlib.import_module(f"scenes.{_m}"), _c)
        except ModuleNotFoundError:     # clip not written yet
            pass
