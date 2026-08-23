# DNA springs on sfGFP — all-atom MD

**[▶ Read the full report](https://extremefattypunch.github.io/DNA-springs-project-2026-27/)**
 · [no-setup mirror](https://htmlpreview.github.io/?https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/index.html)
 · [code and structures](https://github.com/extremefattypunch/DNA-springs-project-2026-27/tree/main/MD%20simulations/sfGFP)
 · [response to review](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_doeke-response.md)
 · [how it was built](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_setup-tutorial.md)

Two amber codons put the tetrazine amino acid **Tet2-Et** at Asp133 and Asn149 on
adjacent strands of sfGFP's β-barrel; a doubly-sTCO-modified DNA strand clicks onto both,
closing a covalent loop. dsDNA is far too stiff to follow the barrel's curvature, so it
bends — pushing the two attachment points apart with a few piconewtons.

![the chimera](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/S3_spring27_overview_still.png)

> ### Two corrections, August 2026
>
> **The first build attached the DNA wrongly.** It put a linker on the 5′ end of *each of
> two different strands*. The experiment modifies **one** strand at both ends, which is
> what makes the chimera a covalent circle and what stops one chimera swapping arms with
> another to form a dimer. As built, every tether bond was a bridge in the bond graph and
> no covalent cycle existed — the loop was closed by base pairing alone. Rebuilt, with a
> topology test the old build fails.
>
> **A claim about the spring's effect is retracted.** The report previously said the
> chimeras displaced the attachment sites more than an equivalent axial pull. That used
> the *unclicked* protein as its zero, and clicking alone moves the sites 0.65 Å — fifteen
> to thirty times any force effect. Against the matched clicked control the sign inverts.
> Both corrections are due to Doeke Hekstra; the [full
> response](https://github.com/extremefattypunch/DNA-springs-project-2026-27/blob/main/feedback/2026-8-22_doeke-response.md)
> works through all nine of his comments.

## What the springs deliver

These forces are Zocchi's *f(x)* evaluated at the end-to-end distance each simulation
settles at — a **model**, not a measurement. A direct measurement, holding a bare duplex
at a series of separations and reading off the restraint force, is running now.
2,550 ns across six systems.

| spring | span *x* | modelled force | stored energy | base pairs open |
|---|---|---|---|---|
| 27 bp | 52.4 Å | **4.82 ± 0.30 pN** | 7.6 k<sub>B</sub>T | 0.11 of 27 |
| 40 bp | 62.5 Å | **2.99 ± 0.06 pN** | 8.7 k<sub>B</sub>T | 0.07 of 40 |
| 40 bp, nicked | 60.8 Å | **2.34 ± 0.04 pN** | 7.0 k<sub>B</sub>T | 0.03 of 40 |

## What the protein does about it

Very little, and that is the result. Against a matched clicked-but-unloaded control, all
three chimeras sit **at or below** it — the springs compress the attachment sites or leave
them alone, and nothing survives a correction for the three comparisons.

Clamping the same two Cβ with a known axial force gives the reason: sfGFP's stiffness
along that coordinate is **230–350 k<sub>B</sub>T nm⁻²**, two to three and a half times
Zocchi's generic ~100. At that stiffness a 4 pN spring moves the sites 0.03 Å, which is
below what 50 ns of equilibrium fluctuation can resolve.

What *is* measurable is that the duplex touches the protein in **100 % of frames** — 20–26
protein heavy atoms within 4.5 Å. The spring is not acting at a distance. The most
compressed system is also the *nicked* one, which delivers the least force: a floppier
duplex lies against the barrel more readily. That ordering is what contact predicts and
the opposite of what force predicts, and it agrees with the companion
[ubiquitin study](https://doi.org/10.64898/2026.07.05.736590) from this lab, which reached
the same conclusion from NMR.

## What the modification itself does

Installing the two Tet2-Et residues shortens Cβ–Cβ by 0.62 Å before any DNA is attached.
It is **not** a secondary-structure change — DSSP is identical everywhere. Asp133 natively
sits in the **left-handed α region** (φ = +57°), which Asp, Asn and Gly like and a
phenylalanine-shaped side chain does not; Tet2-Et is pushed 17.5° out and clicking gives it
back.

![backbone conformation at the two sites](https://raw.githubusercontent.com/extremefattypunch/DNA-springs-project-2026-27/main/docs/figures/ramachandran.png)

The practical consequence: when choosing sites for this amino acid, avoid positive-φ
positions, and treat an aromatic neighbour as a liability — the ring stacks on it in some
constructs and not others, and it does so *differently before and after clicking*, so it
will masquerade as a mechanical effect.

## Does the duplex kink?

No, on this timescale — and that is the expected result rather than an anomaly. Our
curvature is 5.7–6.3° per base-pair step, against the 8.6–9.3°/bp at which the one
directly comparable all-atom calculation nucleates kinks. Zocchi's critical torque is the
most kink-friendly estimate in the field, placing the transition at a gentler bend than a
nucleosome where every independent estimate is 1.3–2.7× stiffer. Separately, base-pair
lifetimes of 1–40 ms predict 10⁻³ opening events in our 270 ns of sampling — so seeing
none bounds any stress-induced barrier reduction at about 5 k<sub>B</sub>T rather than
proving anything.

## Also running

The corrected closed-loop chimeras, a 50-nt ssDNA spring with a matched two-arm control
that cannot circularise, a direct force–extension measurement on the bare duplexes, and
the **ubiquitin** chimera (1UBQ F45W/T9C/K63C with the DBCO-maleimide crosslinker) in five
states, following the lab's preprint.
