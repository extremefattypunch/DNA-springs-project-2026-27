# The Mathematics of Single-Particle Cryo-EM — Manim animations

Work in progress: six silent 1080p30 clips that implement `manim-gpt-animation-plan.md`
(repo root). The full README (commands, scene list, assumptions) is written once the
clips are rendered.

```bash
source env/activate_manim.sh          # Manim env + TinyTeX (holylabs), caches on netscratch
python -m pytest -m "not render" -q   # numerical truth checks
```
