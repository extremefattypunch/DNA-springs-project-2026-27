# Mathematics of single-particle cryo-EM: six Manim illustrations

Six silent 1080p30 clips built from abstract dots, arrows, planes, plots and
equations. The spec is `manim-gpt-animation-plan.md` (repo root). The full README
(commands, scene list, assumptions) is written once the clips are rendered.

```bash
source env/activate_manim.sh          # Manim env + TinyTeX (holylabs); caches on netscratch
python -m pytest -m "not render" -q   # numerical and scope checks
```
