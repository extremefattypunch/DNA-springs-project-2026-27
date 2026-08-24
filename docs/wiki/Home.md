# DNA springs on sfGFP — all-atom MD

**[▶ Read the full report](https://extremefattypunch.github.io/DNA-springs-project-2026-27/)**
 · [no-setup mirror](https://htmlpreview.github.io/?https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/index.html)
 · [code and structures](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/sfGFP)
 · [response to review](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_doeke-response.md)
 · [how it was built](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_setup-tutorial.md)

Two amber codons put the tetrazine amino acid **Tet2-Et** at Asp133 and Asn149 on
adjacent strands of sfGFP's β-barrel; a doubly-sTCO-modified DNA strand clicks onto both,
closing a covalent loop. dsDNA is far too stiff to follow the barrel's curvature, so it
bends — and is supposed to push the two attachment points apart with a few piconewtons.

**4,050 ns over 81 replicates and 17 systems**, in two proteins with two different attachment
chemistries.

![the chimera](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/S3_spring27_overview_still.png)

## The headline

Attaching DNA to this protein pulls the two attachment sites **together** by about a third
of an ångström — and it makes **no measurable difference whether that DNA is a spring**.

10 of the 11 DNA-bearing constructs compress the sites; not one extends them,
which is what a spring pushing its anchors apart would do. And the size of the compression
tracks neither of the two obvious causes:

| against | correlation | p |
|---|---|---|
| the modelled spring force | r = **+0.08** (sign *positive*) | 0.81 |
| how much protein the DNA touches | r = **+0.32** | 0.34 |

Constructs with **no spring at all** — a slack single strand at provably zero tension, or
two arms that cannot circularise — compress by
0.28 Å. Those with a spring compress by
0.35 Å. The difference is not resolved
(p = 0.54).

This is the conclusion the companion [ubiquitin study](https://doi.org/10.64898/2026.07.05.736590)
from this lab reached from NMR — local protein–DNA interaction rather than spring action —
reproduced here independently, in a different protein, with a different chemistry, and
from a different observable.

| construct | protein | modelled force | Δ Cβ–Cβ vs its own clicked control (Å) | p |
|---|---|---|---|---|
| 40 bp nicked, closed loop | sfGFP | 2.49 | **-0.586** | 0.014 |
| 40 bp, closed loop | sfGFP | 3.14 | **-0.482** | 0.030 |
| 27 bp, closed loop | sfGFP | 6.64 | **-0.451** | 0.021 |
| 40 bp nicked, two-arm | sfGFP | 2.37 | **-0.418** | 0.024 |
| two 25 nt arms | ubiquitin | 0 — slack | **-0.418** | 0.072 |
| 50 nt ssDNA, closed loop | sfGFP | 0 — slack | **-0.322** | 0.063 |
| 50 bp dsDNA, closed loop | ubiquitin | 2.23 | **-0.315** | 0.056 |
| two 25 nt arms | sfGFP | 0 — slack | **-0.296** | 0.064 |
| 40 bp, two-arm | sfGFP | 3.04 | **-0.205** | 0.154 |
| 50 nt ssDNA, closed loop | ubiquitin | 0 — slack | **-0.094** | 0.439 |
| 27 bp, two-arm | sfGFP | 6.33 | **+0.021** | 0.878 |

## Why the protein barely moves

Clamping the same two Cβ with a *known* axial force gives the reason: sfGFP's stiffness
along that coordinate is **230–350 k<sub>B</sub>T nm⁻²**, two to three and a half times
Zocchi's generic ~100. At that stiffness a 4 pN spring moves the sites 0.03 Å — below what
50 ns of equilibrium fluctuation can resolve. Doeke Hekstra's physical objection was
right: a force couple this small averages out at equilibrium.

## What the springs deliver, by the model

| spring | span *x* | modelled force | stored energy | base pairs open |
|---|---|---|---|---|
| 27 bp, two-arm | 69.8 Å | **6.33 pN** | 5.3 k<sub>B</sub>T | 0.11 of 27 |
| 40 bp, two-arm | 65.9 Å | **3.04 pN** | 8.4 k<sub>B</sub>T | 0.07 of 40 |
| 40 bp nicked, two-arm | 62.9 Å | **2.37 pN** | 6.9 k<sub>B</sub>T | 0.03 of 40 |
| 27 bp, closed loop | 71.5 Å | **6.64 pN** | 5.1 k<sub>B</sub>T | 0.19 of 27 |
| 40 bp, closed loop | 72.3 Å | **3.14 pN** | 7.9 k<sub>B</sub>T | 0.47 of 40 |
| 40 bp nicked, closed loop | 72.8 Å | **2.49 pN** | 6.3 k<sub>B</sub>T | 0.36 of 40 |
| 50 bp dsDNA, closed loop | 53.3 Å | **2.23 pN** | 10.0 k<sub>B</sub>T | 0.83 of 50 |

These are Zocchi's *f(x)* evaluated at the end-to-end distance each simulation settles at
— a **model**, not a measurement.

## What they deliver, measured

A bare duplex, no protein, held at a series of separations by a
harmonic restraint between the two linker nitrogens. The mean restraint force
⟨k(r − x₀)⟩ *is* the force the duplex exerts — no model in between.

| duplex | separations sampled (Å) | measured force (pN) | Euler load π²B/L² (pN) | measured / Euler | p, slope vs separation |
|---|---|---|---|---|---|
| 27 bp | 33–63 | **25.7 ± 1.0** | 26.8 | 0.96 | 0.49 |
| 40 bp | 50–69 | **15.2 ± 2.3** | 11.9 | 1.27 | 0.29 |
| 40 bp, nicked | 50–69 | **14.0 ± 1.3** | 11.9 | 1.17 | 0.54 |

Two things to read.

**The force does not depend on how far you compress it.** No resolvable slope for any of
the three. A constant load independent of compression is the signature of a **buckled
rod**: past its critical load an elastic rod does not stiffen, it keeps folding at
roughly the same force.

**And the load it holds is the Euler load.** The 27 bp plateaus at
**25.7 ± 1.0 pN** against a classical π²B/L² of
26.8 pN — a ratio of 0.96, with nothing fitted. Textbook
elasticity out of an all-atom force field, with a persistence length nobody supplied.

Zocchi's kinked branch predicts **2.4–6.6 pN** for these same duplexes. We measure
25.7 pN — about **six times higher**, and on the *unyielded* branch, the one
his model says the duplex has already left. So on this timescale the duplex is not a
softened kinked spring; it is an elastic rod holding its buckling load, and the modelled
forces in the table above are if anything **underestimates**.

Each separation was run at two restraint stiffnesses, 33 and 200 kJ mol⁻¹ nm⁻². The soft one is barely a constraint — the duplex settles 2–7 Å past its target — so if the force were an artefact of how hard we pulled, the two would separate. Against the separation actually adopted they agree (p = 0.37, 0.40, 0.90). That is the control that makes this the duplex's force rather than the restraint's.

**What the nick does not do.** Zocchi's model treats the nick as the key softening control and predicts ~1.3× less force. Measured, there is no resolvable difference: 15.2 pN intact against 14.0 pN nicked, p = 0.66. Over 20–30 ns, removing one backbone bond from a rod held by 40 base pairs of stacking does not change its buckling load.

1 runs are **discarded**: the duplex sat close enough to its own periodic image that the restraint force is not its alone. The chimera's own separation (~70 Å for the 27 bp) is outside what these boxes can hold, so an earlier figure quoted at the chimera's geometry is **withdrawn**; repeating those windows in a larger box is queued.

![measured force against extension](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/fig6_force_extension.png)

## Does the duplex kink?

Not in the chimera. Our curvature is **4.8–5.2° per base-pair step**, and
in three 50 ns runs no interior base pair opened. The windows above reach
5.5–9.4°/bp, and there opening is a
probability rather than a threshold: it appears from **5.5°/bp**
and becomes certain — every run in the window — only above
**8.7°/bp**,
with the 8.6–9.3°/bp at which the one directly comparable all-atom calculation nucleates
a type II kink sitting where our probability is already 0.8.

The opening is not diffuse. **Base pair 17** accounts for 79% of the 14 openings seen in the 27 bp duplex, and its incoming step is **CpA** — one of the three pyrimidine–purine steps kinking is known to prefer.

So this is a measurement rather than an absence: the force field is entirely capable of
breaking base pairs, it does so around where the literature says it should, and our
construct sits just below the onset. That **excludes "the force field cannot kink"** as
the explanation — one of the three candidates put to us in review.

![base-pair opening against curvature](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/fig5_kink_threshold.png)

## Does closing the covalent loop matter?

The rebuild (below) made five matched pairs possible, differing only in whether a covalent
path runs through the DNA. All three **duplex** pairs say closure compresses more, by
0.31 Å on average
(paired p = 0.08); both
**single-stranded** pairs fail to reproduce it. A three-pair trend that reverses when two
more pairs are added is a hypothesis, not a finding. Six replicates per system would
settle it.

## What the modification itself does

Installing the two Tet2-Et residues shortens Cβ–Cβ by 0.62 Å before any DNA is attached —
fifteen to thirty times any force effect, which is exactly why using the *unclicked*
protein as the zero produced a spurious result. It is **not** a secondary-structure
change: DSSP is identical everywhere. Asp133 natively sits in the **left-handed α region**
(φ = +57°), which Asp, Asn and Gly favour and a phenylalanine-shaped side chain does not;
Tet2-Et is pushed 17.5° out and clicking gives it back.

![backbone conformation at the two sites](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/ramachandran.png)

The practical consequence: when choosing sites for this amino acid, avoid positive-φ
positions, and treat an aromatic neighbour as a liability — the ring stacks on Tyr200 in
some constructs and not others, and *differently before and after clicking*, so it will
masquerade as a mechanical effect.

> ### Corrections
>
> **The first build attached the DNA wrongly.** It put a linker on the 5′ end of *each of
> two different strands*. The experiment modifies **one** strand at both ends, which is
> what makes the chimera a covalent circle and what stops one chimera swapping arms with
> another to dimerise. As built, every tether bond was a bridge in the bond graph and no
> covalent cycle existed — the loop was closed by base pairing alone. Rebuilt, with a
> topology test the old build fails.
>
> **A claim about the spring's effect is retracted.** The report previously said the
> chimeras displaced the attachment sites *more* than an equivalent axial pull. That used
> the unclicked protein as its zero; against the matched clicked control the sign inverts.
>
> **Every long distance was measured with the wrong convention.** mdtraj applies the
> minimum image by default, which is only valid out to half the box width — 38–48 Å, while
> the springs' ends are 60–130 Å apart. Errors reached 42 Å. Fixed, and every force
> re-derived.
>
> **The duplex curvature was misreported** as 5.7–6.3°/bp; averaging base-pair centres
> over a full helical turn before differencing gives 4.8–5.2°/bp.
>
> The first two are due to Doeke Hekstra; the
> [full response](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_doeke-response.md) works through all
> nine of his comments.

## Still to run

Repaired force–extension windows for the 27 bp duplex inside a box large enough to hold
them, the 40 bp pair at the chimera's own separation, and three seeded type II kinks
(at base pairs 7, 15 and 21) to test whether a kink placed by hand survives. These are
waiting on the lab's GPU fairshare rather than on anything scientific.
