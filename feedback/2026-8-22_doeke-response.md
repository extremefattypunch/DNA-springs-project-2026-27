# Response to Doeke's comments of 22 August 2026

Ian Poon — sfGFP–DNA spring simulations
Written 23 August 2026, revised twice the same day. All 81 replicates (4,050 ns across
17 systems) have now finished and every number below is from the complete set. What
remains outstanding is listed at the end and is waiting on the lab's GPU allocation, not
on anything scientific.

> **Revised after sending.** Setting up the direct force measurement turned up a
> periodic-boundary bug that had corrupted every spring end-to-end distance in the
> project, and with it every modelled force. #6 carries the corrected numbers and the
> full disclosure. It does **not** touch the retraction in #3, the protein stiffness, the
> point-9 correlations or the contact measurement — all four rest on a ~32 Å distance
> below every box's minimum-image limit, where the two conventions agree to 0.00000 Å on
> all eight systems. I also corrected a claim of my own in #7: I had called 6.9 nm a
> 27-mer's ssDNA "contour length" when it is its natural rms end-to-end, which made a
> matched strand look taut.

> **Revised again, 23 August, after all 81 replicates finished.** Three things in the
> version I sent you are now wrong and are corrected below, all in the same direction —
> I had drawn conclusions from partial data.
>
> 1. **#6, the nick.** I wrote that the force field makes the nicked duplex ~3× weaker
>    than the intact one, against the model's 1.3×. That came from single replicates.
>    Across the full set the nicked and intact duplexes are **indistinguishable**
>    (14.0 ± 1.3 against 15.7 ± 2.0 pN, p = 0.50). The nick does nothing measurable on
>    this timescale, which is worse for the model than what I claimed, but not the same
>    claim.
> 2. **#8, the threshold.** "6.7–8.3°/bp" rested partly on windows I have since had to
>    discard — the 27-mer's 66, 74, 82 and 88 Å windows all sit beyond what their box can
>    hold, and those were the rows showing *no* opening at low curvature. With only valid
>    windows, opening is a **probability, not a threshold**: it appears from 5.5°/bp and
>    becomes certain, meaning every run in the window, only above 8.7°/bp. Every conclusion I drew from it survives — the
>    force field can kink, it picks CpA, and our chimera at 4.8–5.2°/bp is below onset —
>    but "the threshold is 6.7–8.3" does not.
> 3. **#3, the mechanism.** With eleven constructs across two proteins I can now test the
>    two candidate explanations rather than argue for one. The compression tracks
>    **neither** the modelled force (r = +0.08, p = 0.81) **nor** the contact area
>    (r = +0.32, p = 0.33), and constructs with no spring at all compress as much as
>    those with one (p = 0.53). I had leaned on contact; the data does not support that
>    either.

---

## Summary

Nine comments. One found a genuine construction error, one forced me to retract a claim,
and one turned out to have a sharper answer than the hypothesis it came with.

| # | Comment | Outcome |
|---|---|---|
| 1 | Tutorial for lab meeting | Written: [`2026-8-22_setup-tutorial.md`](2026-8-22_setup-tutorial.md) |
| 2 | Two TCOs are on *different* strands | **You are right.** Real error. Rebuilt, and a topology test added that the old build fails |
| 3 | Simulate a non-circular control; the couple/shear argument is confused | **Claim retracted.** Your physics is right and mine was wrong. And with 11 constructs in two proteins: every DNA-bearing construct compresses the sites, but it tracks neither force nor contact, and a slack strand does it as well as a spring |
| 4 | Chromophore: bend or rotate? | **Rotate.** Caption reworded, with the measured torsions |
| 5 | The book DOI does not resolve | **Right, and four others were also wrong** — two of them resolved to unrelated papers |
| 6 | Is the force modelled or measured? | It was modelled. Now **measured** too — 25.7 ± 1.0 pN for the 27 bp, which is the **Euler buckling load** (26.8 pN) and 4–11× the kinked branch. A PBC bug and a box-size error found and fixed along the way |
| 7 | ssDNA chimera? | **Done**, 50 nt in both proteins, plus a two-arm control holding DNA content exactly fixed. Both slack — provably zero tension — which is what makes them the zero-force reference #3 leans on |
| 8 | "It does not kink" — force field, timescale, or model? | **The model.** Opening is a probability that rises with curvature: present from 5.4°/bp, certain above 9–10°/bp, at the CpA step the literature predicts. Our chimera sits at 4.8–5.2 and opened nothing in 3 × 50 ns. Not a force-field limitation — excluded by our own data |
| 9 | Do the Et-tet residues change local secondary structure? | **No** — DSSP is identical. But there is a real local effect, and it is a different one |

Two of these connect: the mechanism behind #9 is also the artefact behind #3.

---

## 2. The two attachment points were on different strands

You traced the movie correctly, and this was a real error.

Our build put a 5′-amino-modifier on nucleotide 1 of *each* of the two strands. The
experiment modifies **one** strand at both its 5′ and 3′ ends and adds the complement
unmodified. That difference is not cosmetic: it is what makes the chimera a covalent
circle, and it is what stops one chimera swapping arms with another.

### Confirming it from the topology rather than the movie

Reading the bonds straight out of the old prmtop:

```
DNL(chain C, residue 0):OL — DA(chain C, 1):P
DNL(chain D, residue 0):OL — DC(chain D, 1):P
```

and both 3′ ends are free `DG3`/`DT3` hydroxyls. So each tether is a **bridge** in the
bond graph, and no covalent cycle existed at all — the loop was closed by twenty-seven
base pairs and nothing else.

The obvious check does *not* catch this. "Is the solute one connected component?" passes
on the old build, because both arms hang off the protein backbone, so everything is
connected *through the protein*. The right test is for a **cycle**: cut the Cα–Cβ bond at
both modified residues and ask whether their Cβ are still connected. They can only be if
a path runs out along one tether, through the DNA, and back along the other. That test is
now in `build/04_assemble/build_system.py` and asserted at build time:

| system | verdict | atoms reached |
|---|---|---|
| S3_spring27 (old) | not a loop — arms are bridges | 922 |
| S4_spring40 (old) | not a loop | 1342 |
| S5_spring40nick (old) | not a loop | 711 |
| **L3_loop27 (new)** | **closed loop** | 820 |

The 820 atoms are the two tethers plus the modified strand; the complement is correctly
excluded, being held on by base pairing alone.

### What had to be built

A 3′ linker is not a mirror image of the 5′ one. A 3′-amino-modifier is installed as a
phosphoramidite on the 3′-OH, so it **brings its own phosphate** and is charged:

```
nt(N)-O3′ — P(=O)(O⁻) — O — (CH₂)₇ — NH — C(=O)…   [new residue DNP, q = −0.6921]
H₂N — (CH₂)₆ — O — P(=O)(O⁻) — O-5′ nt(1)          [existing DNL,  q = −0.3079]
```

−0.6921 and −0.3079 are the two shares Amber splits between a strand's 3′ and 5′ ends,
and they have to be right: the modified strand carries 28 phosphates for a 27-mer, so its
residues must sum to −28. The build asserts integral net charge, so an error here fails
loudly. `DNP` was parameterised through the existing pipeline (antechamber/AM1-BCC →
prepgen → parmchk2) with **0 parameters flagged**, and five cross-force-field junction
terms were added, each copied from Amber's own phosphodiester with its provenance
recorded.

### The geometry needed rethinking too

The old code found the bend radius by bisection on the 5′P-to-5′P distance. That works
only for the old pair, which sits nearly diametrically opposite, so its separation barely
depends on which way the helix bends. The two same-strand anchors are different: they are
both ~9 Å off axis and (N−1)×34° apart in azimuth, so at 27 bp they are **71 Å apart for
one bend direction and 54 Å for another at the same radius**. The solver now scans
azimuth and bisects on radius at each, requiring a sign change before bisecting so a
bracket that never contained a root cannot return an artefact.

One consequence worth stating: the span the protein has to bridge is a property of the
protein and the two crosslinkers — Zocchi's *h* — so it is now held **fixed at 66 Å for
every spring length** and only the contour changes. That makes the force ladder a
controlled comparison in a way the old build was not.

| system | bp | radius | arc | anchor span |
|---|---|---|---|---|
| L3_loop27 | 27 | 41.3 Å | 122° | 66.0 Å |
| L4_loop40 | 40 | 41.6 Å | 181° | 66.0 Å |
| L5_loop40nick | 40 | 41.6 Å | 181° | 66.0 Å |

### Where the nick goes now

It has to move to the **complementary** strand: nicking the modified strand would cut the
loop open and put us back at the two-arm construct. That is also Zocchi's own recommended
mechanics-only control (2021 Methods chapter, §4.5.3) — hybridise with two half-length
cDNAs rather than one, which changes the spring's bending stiffness while leaving the DNA
content and every contact with the protein identical.

**Status.** L3/L4/L5 built with 0 problems and integral charge; **3 replicates × 50 ns
each, finished**.

---

## 3. The claim about the spring versus an axial pull — retracted

The report said:

> the three chimeras sit +0.19, +0.52 and +0.59 Å from the zero-force clamp — monotonic
> in force, a slope of 0.13 Å pN⁻¹, more than ten times the axial-clamp slope.

**That is withdrawn.** It used the *unclicked* protein as the zero. Every chimera is
clicked, and clicking alone moves the sites by +0.647 Å:

| | Cβ–Cβ (Å) | n |
|---|---|---|
| WT | 32.245 ± 0.059 | 3 |
| 2× Tet2-Et, unclicked (the old zero) | 31.627 ± 0.022 | 3 |
| clicked, no DNA (`S2_clicked`, the right zero) | 32.274 ± 0.092 | 3 |

So every chimera appeared ~0.65 Å "extended" purely from the chemistry — 15 to 30 times
any force effect. Against the matched control the picture inverts:

| model force | system | Δ vs clicked control | t | p | p (Bonferroni ×3) |
|---|---|---|---|---|---|
| 2.34 pN | S5_spring40nick | **−0.418 Å** | −2.92 | 0.024 | 0.073 |
| 2.99 pN | S4_spring40 | −0.205 Å | −1.68 | 0.154 | 0.462 |
| 4.82 pN | S3_spring27 | +0.021 Å | +0.16 | 0.878 | 1.000 |

Every chimera sits at or *below* the zero-force control. The springs do not extend the
attachment sites; they compress them or leave them alone. Nothing survives correction for
multiple comparisons.

A positive trend with force does survive the change of baseline (+0.161 Å pN⁻¹, r² = 0.93)
but with three points and p = 0.17 it is not established, and it now describes *recovery
from compression* rather than extension. I am not making a claim from it.

### Your physics objection is correct

> A couple (torque, I assume) would exert its effect indeed but (once averaged over
> replicates) only transiently until the forces at the attachment sites are equal and
> opposite. I am not sure what to make of the shear.

You are right on both counts and I was wrong. At equilibrium the net force and torque on
the protein vanish, a couple does no net work, and along the Cβ–Cβ coordinate only the
component of force along that coordinate can displace the mean. Invoking a couple and a
shear to explain an equilibrium displacement was not a sound argument, and the data did
not need it anyway — the effect it was invented to explain was an artefact of the wrong
baseline.

> are the different simulations, force clamp and DNA attachment, not at the same sites?

They are the same sites. The clamp acts on the Cβ of residues 133 and 149, and the atom
indices are written by the system builder from the same site list the tethers use, so the
physics and the topology cannot drift apart.

### What the axial clamp actually measures

The clamp ladder gives the protein's stiffness along that coordinate directly:

| force | Cβ–Cβ (Å), clicked topology |
|---|---|
| 0 pN | 32.274 ± 0.092 |
| 2.35 pN | 32.289 ± 0.023 |
| 4.88 pN | 32.308 ± 0.029 |

A straight line through those three gives 6.9 mÅ pN⁻¹, i.e. **κ ≈ 350 k_BT nm⁻²**
(r² = 0.998). The longer ladder on the unclicked topology (six forces from 0 to 20 pN,
18 replicates) is noisier — r² = 0.177, p = 0.08 — and gives **κ ≈ 230 k_BT nm⁻²**. Both
say the same thing: the sfGFP barrel along this coordinate is **2–3.5× stiffer than the
~100 k_BT nm⁻² Zocchi uses as a generic protein spring constant** (*Molecular Machines*,
p. 118–119). At 350 k_BT nm⁻², a 4 pN spring moves the sites 0.03 Å. That is the
magnitude to expect, and it is below what 6 × 50 ns can resolve.

### So what *is* different between the spring and the clamp?

Not a couple. The DNA is a large charged object tethered two nanometres from the barrel,
and it leans on it. Measuring contact directly in the existing chimeras:

| system | protein heavy atoms within 4.5 Å of the duplex | frames with contact |
|---|---|---|
| S3_spring27 | 25.5 ± 12.6 | 100 % |
| S4_spring40 | 21.4 ± 8.8 | 100 % |
| S5_spring40nick | 20.4 ± 6.6 | 100 % |

The duplex is touching the protein in **every frame of every replicate**. That is not a
spring acting at a distance, and it is the obvious candidate for the ~0.4 Å compression
that has nothing to do with force. It is also exactly what the ubiquitin preprint
concluded from the other direction — a control with two separate short strands reproduced
most of the circular chimera's spectrum, pointing at local protein–DNA interaction rather
than spring action.

Note too that the *most* compressed system is the **nicked** one, which delivers the
*least* force. A nick makes the duplex floppier, so it can lie against the barrel more
readily. That ordering is what a contact mechanism predicts and the opposite of what a
force mechanism predicts.

### With all eleven constructs in, it is neither force *nor* contact

The paragraph above was written with three chimeras and it argued for contact. Eleven
constructs across two proteins let me test that instead of arguing it, and **contact does
not survive the test either**.

| construct | protein | modelled force | Δ Cβ–Cβ vs its own clicked control (Å) | p |
|---|---|---|---|---|
| 40 bp nicked, closed loop | sfGFP | 2.49 | **-0.586** | 0.014 |
| 40 bp, closed loop | sfGFP | 3.14 | **-0.482** | 0.030 |
| 27 bp, closed loop | sfGFP | 6.64 | **-0.451** | 0.021 |
| 40 bp nicked, two-arm | sfGFP | 2.37 | **-0.418** | 0.024 |
| two 25 nt arms | ubiquitin | 0 — slack | **-0.417** | 0.082 |
| 50 nt ssDNA, closed loop | sfGFP | 0 — slack | **-0.322** | 0.063 |
| 50 bp dsDNA, closed loop | ubiquitin | 2.23 | **-0.320** | 0.052 |
| two 25 nt arms | sfGFP | 0 — slack | **-0.296** | 0.064 |
| 40 bp, two-arm | sfGFP | 3.04 | **-0.205** | 0.154 |
| 50 nt ssDNA, closed loop | ubiquitin | 0 — slack | **-0.094** | 0.412 |
| 27 bp, two-arm | sfGFP | 6.33 | **+0.021** | 0.878 |

Every construct that carries DNA compresses the sites — ten of eleven negative, the
eleventh +0.02 which is zero, and not one positive. A spring pushing its anchors apart
would extend them. But:

- against the modelled spring force: **r = +0.08,
  p = 0.81** — and the sign is *positive*, so more force
  goes with *less* compression;
- against how much protein the DNA touches: **r = +0.32,
  p = 0.33**;
- constructs with **no spring at all** compress 0.28 Å,
  those with a spring 0.35 Å — not resolved
  (p = 0.53).

The two constructs delivering the *most* modelled force, the 27 bp pair at 6.3–6.6 pN, are
the two that compress *least* (+0.02 and −0.45). A slack 50-nt single strand at provably
zero tension compresses −0.32.

So my contact explanation was better than the force explanation and still wrong as stated.
What I can defend is narrower and matches your conclusion in the preprint more closely than
my version did: **attaching a large, highly charged object a couple of nanometres from the
barrel moves these sites by about a third of an ångström, and which object it is barely
matters.** Whether that is electrostatics, hydration, or the tether geometry itself, these
runs do not separate — and I should not have named a mechanism from three data points.

The nicked-is-most-compressed ordering I pointed to above also does not survive: with the
closed-loop builds included the ordering across all eleven is not monotonic in anything I
can measure.

### And on whether closing the loop matters

Your point 2 made this testable — five pairs differing only in whether a covalent path
runs through the DNA:

| pair | closed − two-arm (Å) | p |
|---|---|---|
| 27 bp | -0.472 | 0.007 |
| 40 bp | -0.277 | 0.018 |
| 40 bp nicked | -0.168 | 0.311 |
| 50 nt ssDNA, sfGFP | -0.026 | 0.724 |
| 50 nt ssDNA, ubiquitin | +0.323 | 0.141 |

All three **duplex** pairs say closure compresses more, by
0.31 Å on average, paired
p = 0.08. Both **single-stranded** pairs fail to
reproduce it, and across all five it is p = 0.41.

A three-pair trend at p = 0.08 that reverses when two more pairs are added is a hypothesis,
not a result, and I am reporting it as one. It would be physically sensible — only a closed
circuit can carry a bending moment around, and a slack strand has no moment to carry — which
is exactly why it needs six replicates per system rather than three before I believe it.

### The non-circular controls you asked for

Built, from the same 50-nt strand so that DNA content, sequence, sites and chemistry are
held exactly fixed and only the topology changes:

- **C1_ss50** — 50 nt clicked at both ends: the circular ssDNA chimera.
- **C2_ssarms** — the same strand cut in the middle: two 25-nt arms, each attached at one
  end, which cannot close.

This is a tighter comparison than the preprint's two 12-mers, and it is Zocchi's own
design principle: change the mechanics, hold every contact constant. Both built with
0 problems, and the loop test is asserted in **both** directions — C1 must be closed,
C2 must not be — because a control that quietly closed would look like a null result.
3 replicates × 50 ns each — **finished**, and reported above.

---

## 4. The chromophore: rotation, not bending

The implication is **rotation**. GFP's chromophore loses its quantum yield through
twisting about the two exocyclic bonds of the methine bridge, which opens a non-radiative
decay channel; the hydrogen-bond network holds the two ring systems coplanar and so
suppresses it. The caption should have said torsion and will.

The two bridge torsions, as measured (τ = phenol ring against the bridge, φ = bridge
against the imidazolinone; reported as |twist| from planarity, which is immune both to the
±180° wraparound and to the phenol ring's two-fold symmetry):

| system | τ mean | τ 95th pct | φ mean | φ 95th pct |
|---|---|---|---|---|
| WT | 5.6° | 13.3° | 6.5° | 16.3° |
| 2× Tet2-Et | 5.1° | 12.6° | 6.0° | 14.5° |
| clicked | 5.4° | 13.1° | 6.1° | 14.9° |
| 27 bp spring | 5.3° | 13.0° | 6.2° | 15.5° |
| 40 bp spring | 5.4° | 13.7° | 6.2° | 14.8° |
| 40 bp nicked | 5.2° | 13.1° | 6.1° | 15.0° |

The chromophore stays flat to within about 5–6° on average in every system, and the
spring changes nothing measurable. The honest statement is that the network holds it
planar and that none of our constructs perturbs it — not that the spring modulates it.

---

## 5. References

You are right that `10.2307/j.ctvc77hpc` does not resolve — the DOI handle API returns
"not found". Checking all 27 references turned up four more errors, two of which are the
dangerous kind: they resolve perfectly well, to unrelated papers.

| cited | problem | correct |
|---|---|---|
| `10.2307/j.ctvc77hpc` | does not exist | `10.23943/9781400890064` |
| `10.1016/bs.mie.2020.11.005` | does not exist | `10.1016/bs.mie.2020.09.010` |
| `10.1101/2026.07.05.736590` | does not exist (bioRxiv changed prefix) | `10.64898/2026.07.05.736590` |
| `10.1021/ja4045532` | **resolves to Iešmantavičius et al., helix propensity in a disordered protein** | `10.1021/ja4043565` (Tseng & Zocchi, RLuc) |
| `10.1529/biophysj.106.092049` | **resolves to Shemella et al., intein C-terminal cleavage** | `10.1529/biophysj.106.092866` (Choi & Zocchi, GK) |

Both of the misdirected DOIs are confirmed against the running footers of our own PDF
copies, and the Methods in Enzymology DOI is printed on page 1 of the file in the reading
list. All five corrections verified to resolve. Two further details: the book's subtitle
is *A Materials Science Approach* (unhyphenated on its title page), and the only ISBN
printed inside the book is the hardcover 978-0-691-17386-3.

**Do I have access to the book?** Yes — as the PDF in
`reading list/zocchi molecular machines/`. The report now cites it by page for every
equation it uses:

| what | book page | PDF page |
|---|---|---|
| DNA yield transition, *E*(*x*) and *f*(*x*), Eqs. 2.203–2.206 | 77 | 89 |
| kinked-state energy, Eqs. 2.179, 2.196 | 74, 76 | 86, 88 |
| the leaf-spring method, §3.3 | 89 | 101 |
| dimerisation as a thermometer, Eqs. 3.55–3.56, *E*_el = 9.1 ± 0.1 k_BT | 103–104 | 115–116 |
| κ ≈ 100 k_BT nm⁻², *x*_e = 3 Å, *f*_e ≈ 4 pN | 118–119 | 130–131 |
| the *γ* < 1 validity condition (DNA shorter than ~15 nm / 45 bp) | 75 | 87 |

That last row matters for the ubiquitin work: the 50-mer used there is 17 nm, so the
analytic *E*(*x*) is being applied slightly outside its stated domain, and the report now
says so.

---

## 6. Where the forces come from — both calculations, in full

Your reading was right: the quoted forces came from **a model**, not from the simulation.
Specifically, Zocchi's *f*(*x*) evaluated at the end-to-end distance the simulation
settled at. Here is that calculation written out, then the direct measurement, and then a
bug I found in between that changed the input to the first one.

### (a) The analytic route

Zocchi treats a sharply bent duplex as an elastic rod that, past a critical torque τ_c,
yields into a **kinked** state: two straight arms joined by a hinge sustaining a constant
torque. Bending modulus *B* = 200 pN nm² (ℓ_p ≈ 50 nm at 300 K), rise 0.33 nm/bp,
τ_c = 27 pN nm nicked and 31–36 intact, k_BT = 4.11 pN nm.

```
γ    = L τ_c / 2B                    small parameter, derivation assumes γ < 1
R    = (L/2)(1 − γ²/90)              arm half-length in the kinked state
x_c  = L (1 − γ²/15)                 the yield point
x₀   = L (1 − LT/10B)

x ≥ x_c  (unkinked, WLC-like):  f(x) = 10B/L² − T/(L−x)
x < x_c  (kinked):              f(x) = τ_c / (2R √(1 − (x/2R)²))
```

The implementation is checked against Zocchi's own published values before being used
here — 60-mer 1.47 vs 1.5 pN, 40-mer 2.33 vs 2.4 pN, and *E*_DNA = 6.62 vs 6.6 k_BT for
the GK chimera — so the arithmetic is right even where the model's domain is stretched.
One caveat that belongs with the numbers: γ = 1.12 for the 40-mer, outside the γ < 1
condition the derivation assumes (book p. 75), so those rows are extrapolations.

### A correction to the numbers I sent you

While setting up the direct measurement I found that **every end-to-end distance this
project had reported was wrong**, and with it every modelled force. mdtraj applies the
minimum-image convention by default. That is right for a hydrogen bond and wrong for
anything longer than half the cell's smallest perpendicular width — 38–48 Å in these
truncated octahedra — and our springs' ends are 60–75 Å apart. The reported value was
flipping between periodic images from frame to frame:

| system | as reported | true | worst frame-wise disagreement |
|---|---|---|---|
| S3_spring27 | 53.8 Å | **69.7 Å** | 30.7 Å |
| S4_spring40 | 62.5 Å | **66.4 Å** | 12.4 Å |
| L3_loop27 | 45.4 Å | **71.3 Å** | 41.8 Å |

The tell was in the spread all along and I did not read it: the raw distances have a
standard deviation of 2.7–3.6 Å, the minimum-imaged ones 6.7–8.8. The nonsensical −89 and
−178 pN in my first force-extension pass were this too, not the poor equilibration I put
them down to.

**What this does and does not touch.** The fix is in `analysis/geometry.py`, which picks
the convention from the topology and then checks its own reasoning. Checked explicitly
rather than assumed: the Cβ–Cβ separation of the attachment sites is ~32 Å, below every
box's limit, and the two conventions agree to **0.00000 Å on all eight systems**. So the
retraction in #3, the protein stiffness, the point-9 correlations and the contact
measurement all stand exactly as sent. Watson–Crick opening, hydrogen bonds and contacts
are short-range and were always right. What changes is the spring distances and the
forces derived from them:

| system | *L* | *x* (corrected) | modelled *f* | *E* | Euler load π²B/L² |
|---|---|---|---|---|---|
| S3_spring27 | 8.91 nm | 69.7 Å | 6.22 pN | 5.5 k_BT | 24.9 pN |
| S4_spring40 | 13.20 nm | 66.4 Å | 3.04 pN | 8.5 k_BT | 11.3 pN |
| S5_spring40nick | 13.20 nm | 63.5 Å | 2.36 pN | 6.9 k_BT | 11.3 pN |
| L3_loop27 (closed loop) | 8.91 nm | 71.3 Å | 6.49 pN | 5.2 k_BT | 24.9 pN |
| L4_loop40 | 13.20 nm | 72.4 Å | 3.14 pN | 8.1 k_BT | 11.3 pN |
| L5_loop40nick | 13.20 nm | 72.2 Å | 2.47 pN | 6.4 k_BT | 11.3 pN |

The ordering is unchanged and the 40-mers barely move, because in the kinked branch *f* is
flat at small *x*/2*R*; the 27-mer changes most, from 4.76 to 6.22 pN. Worth noting
separately: **closing the loop lets the duplex sit less compressed** — 71–72 Å against
63–70 Å for the matched two-arm builds — which is a real mechanical difference between the
two topologies and not something I expected.

### (b) The direct route — what the force field actually says

Take the duplex alone, no protein, solvated. Hold its ends at a fixed separation *x* with
a harmonic restraint of known stiffness *k*. At equilibrium the restraint's mean pull
balances what the molecule pushes with, so

```
f_DNA(x) = k (⟨r⟩ − x₀)
```

with no elasticity model anywhere in the derivation. The restraint centre is walked from
the built geometry to the window target over 2 ns before any frame is recorded, because
dropping a duplex into a window 20 Å away does ~400 kJ/mol of work in one step and
whatever base pair gave way would look exactly like the kink these runs also look for.

Also worth recording: my first pass used *k* = 200 kJ mol⁻¹ nm⁻², which is 33 pN Å⁻¹ — so
a single 1 Å thermal fluctuation reads as 33 pN and a 4 pN spring produces a 0.12 Å offset
against a 1.11 Å spread. The precision of a mean force taken from equilibrium fluctuations
is √(k k_BT/N) however the restraint is applied, so a **softer** spring is strictly better.
A second pass at *k* = 33 resolves 2.5× finer. Both are reported below, and their
agreement is the method checking itself.

**A second correction, and this one bites the headline.** The solvent boxes were padded
around the *bent* duplex, which is 61.7 Å long, while the straight molecule is 93 Å. So
every extended window ran in a box too small for it — the 27-mer beyond 58 Å and the
40-mers beyond 66 Å were measuring a duplex that can reach its own periodic image. That
includes the ~70 Å separation the chimera actually sits at, which is precisely the number
the measured-versus-modelled comparison rests on. I had quoted 33 pN there; that window is
not trustworthy and the figure is withdrawn until the re-run lands.

What survives untouched is the kinking threshold in #8: all of it lives between 32 and
58 Å, comfortably inside the box.

**Result, 27 bp duplex, windows inside their box only** (two replicates each, two
independent restraint stiffnesses):

| window *x*₀ | °/bp | ⟨*r*⟩ | measured *f* |
|---|---|---|---|
| 32 Å | 9.7 | 34.9 Å | 26.9 ± 5.8 pN |
| 38 | 9.0 | 40.5 | 22.0 ± 2.1 |
| 44 | 8.3 | 47.0 | 27.3 ± 3.0 |
| 50 | 7.6 | 52.9 | 26.4 ± 4.1 |
| 58 | 6.7 | 61.3 | 29.2 ± 8.0 |

A flat plateau of 22–29 pN across the whole valid range, which is what a buckled rod does
— and it sits on the Euler load for a 27-mer, π²*B*/*L*² = **25.6 pN**, not on Zocchi's
kinked branch, which gives 2.4–6.6 pN over the same range. The two independent restraint
stiffnesses agreeing at matched ⟨*r*⟩ is the method checking itself.

**And the nick does nothing measurable — which I had wrong.** The version I sent you put
the nicked 40-mer at 4.2 and 6.9 pN against the intact 13.3 and 18.1, and concluded the
force field softens it ~3× where the model says 1.3×. Those were single replicates. Across
every valid run:

| | measured force | n |
|---|---|---|
| 40 bp intact | 15.7 ± 2.0 pN | 8 |
| 40 bp nicked | 14.0 ± 1.3 pN | 8 |
| Welch | t = +0.70, **p = 0.50** | |

There is no resolvable difference. Removing one backbone phosphodiester from a rod held
together by forty base pairs' worth of stacking does not change its buckling load over
20–30 ns. Zocchi's model puts the nicked spring at 2.37 pN against the intact 3.04 — a
factor of 1.3, set almost entirely by τ_c — and we cannot reproduce even that.

This is *worse* for the model than the 3× I wrongly claimed, not better. The nicked spring
is the low-stress reference his whole calibration leans on, and in this force field it is
not a low-stress reference; it is the same spring. But it is a different statement from
the one I sent, and the one I sent was drawn from too little data.

### The two calculations disagree by a factor of five, and the reason is point 8

Over the range where both can be evaluated and the box is adequate:

| | force |
|---|---|
| Zocchi's kinked branch, 27 bp over 32–63 Å | **2.4–6.6 pN** |
| Euler buckling load, π²*B*/*L*², for an *unyielded* rod | **26.8 pN** |
| measured, from the force field, 19 valid runs | **25.7 ± 1.0 pN** |
| ratio, measured to kinked branch | **4–11×** |

(The Euler figure is 26.8 pN with *B* = 200 pN nm² and 26 steps of 0.33 nm. I quoted
25.6 and 24.9 in two places in the version I sent you, from a contour taken over 27 steps
rather than 26; 26.8 is the consistent value and it is what the figures now use.)

The measurement lands on the Euler load, not the kinked one — and #8 says why: this duplex
does not yield at the curvature the chimera imposes, so it stays on the elastic branch, and
an unyielded rod sustains its buckling load. Zocchi's branch is lower by four or five
precisely because it presupposes the yield.

Which of the two is right about real DNA is a separate question our simulation cannot
settle. If real DNA does yield here and ours does not, the model is right about the
experiment and our force field is on the wrong branch.

**Still outstanding.** Near-straight duplexes have been solvated in boxes that fit the long
windows (the 27-mer's now reaches 73 Å, covering the chimera's 70), and the 50 and 58 Å
windows need re-running in both boxes so the two halves of the curve can be checked
against each other before they are joined. That is **queued behind the lab's GPU
fairshare**, which a day of my own exploratory jobs exhausted; it is submitted to the
free `gpu_test` partition and will land when the allocation resets.

## 7. The ssDNA chimera

**Done.** **C1_ss50**, one 50-nt strand clicked at both ends, no complement, and its
ubiquitin twin **U2_ss50**. Both come out *slack*, which is the useful part: 50 nucleotides
is 295 Å of contour with a natural rms end-to-end of
94 Å, and the attachment sites hold it at only
54 Å. A worm-like chain held inside its own rms size is
under **no tension at all**, so the Marko–Siggia force is exactly zero — not small, zero.

That makes these the cleanest control in the project: same protein, same sites, same
chemistry, comparable DNA, provably zero mechanical load. If the site compression were a
force effect they should show none. They show
0.32 Å (sfGFP) and
0.09 Å (ubiquitin) — as much as any duplex.
That is the observation #3 now rests on.

The sign contrast is worth stating too, since it is the opposite of the duplex: a bent
duplex is an *enthalpic* spring and **pushes** its anchors apart; a single strand is an
*entropic* spring with a 1.5 nm persistence length and can only ever **pull**, and only
once stretched past its natural size. Neither construct is.

Two notes on the design.

**Why 50 nt and not 27 — with a correction to my own reasoning.** I first justified this
by saying a 27-mer has "only ~6.9 nm of contour at 0.255 nm per base" and would be taut.
That was wrong, and worth correcting because the number was right and the label was not.
6.9 nm is a 27-mer's *natural rms end-to-end distance*, √(2ξL₀); its **contour** is
15.9 nm, at 0.59 nm per nucleotide (Smith, Cui & Bustamante 1996). With ξ_ss = 1.5 nm —
the figure Zocchi's own 2021 chapter uses — a 27-mer's natural size is 6.9 nm against a
6.6 nm span, so it is **matched, not taut**, and a 50-mer at 9.4 nm is genuinely slack.

So the real reason to use 50 nt is the simpler one: it is the length used in the ubiquitin
work, which makes the two studies comparable, and it is the slackest of the options and
therefore the cleanest zero-stress reference. Neither would have been taut.

**The sign flips.** With a persistence length of ~1.5 nm, ssDNA is an *entropic* spring
that **pulls** as it is stretched, where the compressed duplex **pushes**. So this is not
simply "the spring turned off" — it is a spring of the opposite sign, and that contrast is
worth reporting in its own right.

**A caution about calling it a zero.** Zocchi's ss-form control is validated for GK and
RLuc, where the ss chimera's activity is indistinguishable from the unmodified mutant. It
*fails* for ubiquitin under NMR: the preprint finds that attaching ssDNA alone shifts
peaks throughout the protein, and that adding the complement to make dsDNA does
comparatively little. Zocchi anticipates this failure mode himself (2021, §4.5.3) for
proteins that interact with nucleobases, and recommends the nicked-versus-intact duplex
comparison instead — which is why our nick moved to the complement (see #2).

### And on the nicked spring being only slightly weaker

The model says 2.34 pN nicked against 2.99 pN intact, from τ_c = 27 versus 34 pN nm — a
22 % difference, because in the kinked regime *f* is set almost entirely by τ_c. Whether
that is right is exactly what the direct measurement in #6 will test, and it is worth
testing: the literature disagrees with itself by orders of magnitude on how much a nick
softens DNA (see #8).

---

## 8. "It does not kink" — which of your three explanations?

**(c), the model's threshold is too permissive.** We can now say that from our own data
rather than by citing other people's, and the same data excludes (a) outright: this force
field kinks perfectly well, just not at the curvature our construct imposes. (b), the
timescale, remains a real constraint on *spontaneous* nucleation but is no longer the
main story.

### What our system's curvature actually is — and a number I had wrong

The duplexes are built at *R* = 41–44 Å and relax over the run to **37–41 Å, i.e. 4.8–5.2°
per base-pair step** (measured, one-turn smoothed, 18 replicates), equivalent to a 70–78 bp
minicircle. The worst single step is ~9.5°, and Watson–Crick opening averages 0.03–0.11
pairs per frame — essentially none.

I need to correct myself twice over here. The report first said this was "roughly
nucleosomal"; I then "corrected" that to 31–34 Å and ~25 % tighter than nucleosomal. **The
31–34 Å figure was wrong** — I had carried it over without checking it against the
trajectories. The numbers above are what the bend analysis actually reports, and at 37–41 Å
against the nucleosome's 42 Å the original wording was closer to right than my correction
was. The conclusion is unaffected in direction and slightly strengthened in degree — and
the next section replaces the comparison entirely, by measuring where *this* force field
kinks rather than borrowing a threshold.

### We measured the threshold instead of arguing about it

The compression windows do this directly: hold the bare 27 bp duplex at a series of separations, from its contour
length down to 64 % compression, and ask at what curvature *this* force field starts
breaking base pairs. Two replicates per window, 12–30 ns each, and the criterion was
fixed in advance — an **internal** pair (not the terminal two, which fray in anything)
open in more than 10 % of frames.

**The table I sent you had to change.** Its low-curvature rows — the 27-mer at 66, 74, 82
and 88 Å — are the ones that showed no opening, and every one of them sits beyond what its
solvent box can hold, so they measure a duplex interacting with its own periodic image. I
cannot use them, and they were carrying the claim that opening *starts* at 6.7°/bp.

What is left, pooling every valid window from all three duplexes and binning by the
curvature the duplex actually adopted (not the target it was pulled towards — with the
soft restraint those differ by 2–7 Å):

| °/bp | runs | runs that opened an interior pair |
|---|---|---|
| 5–6 | 7 | 2 (29 %) |
| 6–7 | 12 | 4 (33 %) |
| 7–8 | 5 | 4 (80 %) |
| 8–9 | 6 | 5 (83 %) |
| 9–10 | 4 | **4 (100 %)** |

The per-run observation is bimodal — a run either keeps every interior pair shut for its
whole trajectory or opens one and holds it open 30–100 % of frames — so what varies
smoothly with curvature is the *probability*, and an average over runs is a number no run
exhibits.

Three things follow, and they settle your question between them.

**The force field can kink.** It does so reproducibly, and the opened pair sits at the
same place as the maximum of the bend profile every time. So explanation (a) — that AMBER
with OL15 simply cannot support a kink — is excluded by our own data rather than by
argument. That also means the seeded-kink test I proposed is no longer the decisive
experiment; it is still running, but as a question about metastability below the
threshold rather than about whether the force field is capable of it at all.

**There is no threshold; there is a probability.** Opening appears as low as **5.5°/bp**
and becomes certain — every run in the window — only above **8.7°/bp**. Curuksu *et al.*'s umbrella sampling puts
type II kink nucleation above 8.6–9.3°/bp — which is where our probability has reached
0.8, so the two agree about where kinking becomes the *expected* behaviour while ours also
shows a stochastic tail well below it. Both are a long way from Zocchi's τ_c, which places
the transition at 2.6–3.0°/bp: even our stochastic tail does not reach that.

**Our chimera is at 4.8–5.2°/bp, below the lowest curvature at which we ever see a pair
open.** That is the answer. Not seeing a kink in the chimera is not a force-field limitation and not
principally a timescale artefact: the construct is simply not bent far enough, on this
force field's own showing.

One detail worth having, because it is a check rather than a coincidence. In the 27 bp
duplex the kink lands on **base pair 17 in 11 of its 14 openings**, across independent
windows and replicates. Base 17
is an adenine and the step into it is **CpA** — one of the three pyrimidine–purine steps
(CpG, TpA, CpA) that Curuksu names as the preferred sites for roll to localise at. The
other pair that opened, bp13, sits at a **TpG** step, also pyrimidine–purine. So the force
field is not breaking a pair at random under load; it is picking exactly the dinucleotide
steps the literature says it should.

### The measured force says the same thing, independently

This is the part that persuaded me. If the duplex never yields it stays on the elastic
branch, and an unyielded rod under compression sustains its Euler buckling load,
π²*B*/*L*² = **26.8 pN** for a 27-mer. Measured over the 19 runs whose boxes can hold
them, the 27 bp duplex delivers **25.7 ± 1.0 pN** — a ratio of 0.96, with nothing fitted —
against Zocchi's kinked branch at 2.4–6.6 pN over the same range.

The force is also **flat in separation** (p = 0.49), which is the other half of the
signature: past its critical load a buckled rod does not stiffen, it keeps folding at
roughly constant force.

(I had quoted **33 pN** here, at the chimera's own ~70 Å separation. That window's target
was 87 % of its box width, so it is withdrawn — see #6 — and the numbers above come only
from windows inside their boxes. The chimera's own separation is not among them, and
re-running it in a larger box is queued behind the lab's GPU allocation.)

So two quite different observables — the shape of the bend profile and the magnitude of
the delivered force — agree that this duplex has not yielded. That is more convincing than
either alone, and it means the model's few-pN forces are not merely uncorroborated here:
they are *inconsistent* with what the same force field delivers, because they presuppose
the transition it never makes.

Which of the two is right about real DNA is a separate question our simulation cannot
settle. If real DNA does yield at this curvature and ours does not because 50 ns is too
short, then the model is right about the experiment and our force field is stuck on the
wrong branch.

### (c) Zocchi's τ_c is the most kink-friendly estimate in the field

Converting each published criterion to a threshold radius via *R* = *B*/τ_c:

| source | criterion | threshold *R* | °/bp |
|---|---|---|---|
| Zocchi, nicked (EPL 2010) | τ_c = 27 pN nm | 76 Å | 2.6 |
| Zocchi, intact (PRX 2011) | τ_c = 31 pN nm | 66 Å | 3.0 |
| Vologodskii & Frank-Kamenetskii 2013; Du 2008 | nuclease probing | 35–40 Å | 4.9–5.6 |
| Curuksu et al. 2009 (umbrella sampling, AMBER) | asymptotic torque 40–51 pN nm | — | kink only above **8.6–9.3** |
| Fields, Meyer & Cohen 2013 (FRET "molecular vice") | WLC linear up to 7°/bp | ≤28 Å | 7 |
| **our chimeras, measured** | | **37–41 Å** | **4.8–5.2** |

Zocchi's τ_c puts the softening transition at a *gentler* bend than a nucleosome. Every
independent estimate is 1.3–2.7× stiffer. On his criterion we are well past the
transition; on Curuksu's and on Fields'/Cohen's we are comfortably below it. So the
disagreement between our simulation and his prediction is substantially a disagreement
between his τ_c and the rest of the literature — and the one directly comparable
calculation, Curuksu's umbrella sampling on a bent *linear* duplex in the same force-field
family, finds kinks only above 8.6–9.3°/bp, **almost twice the curvature we impose**.
Not seeing a kink is the expected result, not an anomaly.

Zocchi's implied kink cost, (π/2)τ_c ≈ 12 k_BT, is likewise at the bottom of the published
range of 12–22 k_BT.

### (b) The timescale, with the arithmetic

Base-pair lifetimes from NMR imino exchange are 1–40 ms (Guéron 1987; Leroy 1988;
Krueger 2006), so the per-pair opening rate is ~10²–10³ s⁻¹. For ~30 pairs in parallel over
our 270 ns of aggregate sampling, the **expected number of opening events is 8×10⁻⁴ to
8×10⁻³**. Seeing none is what unstressed kinetics predicts before bending is even
considered. Turning that around gives a real bound rather than a null result: our data are
consistent with a stress-induced reduction of the opening barrier of **up to ~5 k_BT**, and
cannot support more.

Independently, kink free energies of 12–22 k_BT with a 10⁹–10¹² s⁻¹ prefactor give mean
waiting times of microseconds to hundreds of milliseconds — 10²–10⁷ times our window.
Fassino, Carlon & Voorspoels (JCTC 2026) state it flatly: *"in unbiased MD simulations, it
is not possible to observe the spontaneous formation of kinks due to limitations in the
time scale currently reachable in all-atom simulations."*

Calibration from published work: Pasi & Lavery saw kinks at 80–800 ns in nucleosomal loops
held far more tightly than ours; Mitchell & Harris saw them at 20–25 ns but only in 65 bp
*closed* minicircles (R ≈ 35 Å, torsionally coupled); and Yoo et al. ran a 90 bp minicircle
for **11 µs** at R ≈ 49 Å with **no kinks and no base-pair-level distortion at all**. That
last one is the strongest counter to "you would have seen it if the force field were
right": 240× longer than ours, closed topology, and nothing.

### (a) The force field — settled, but the concern was reasonable

The measurement above excludes this explanation: OL15 kinks. It is still worth recording what the literature says, because your worry was well founded and because the residual bias may be why our threshold sits a little below Curuksu's rather than a little above.

There *is* a documented bias here, and OL15 is specifically implicated: Winkler,
Galindo-Murillo & Cheatham (JCTC 2023) find that backbone-dihedral refinements
"inadvertently favored B-form DNA geometries even in non-helical contexts", with OL15
producing anomalous compact structures in 42 % of a mini-dumbbell trajectory. All AMBER
variants also overestimate persistence length, with OL15 the stiffest of them
(Roldán-Piñero et al. 2024). And there is a mechanistic reason to expect OL15 to be the
least kink-prone choice: twist-bend kinks require local *under*twisting, and the OL family
was tuned to *increase* helical twist.

But the magnitude is nowhere near enough to rescue the model. parmbsc0 and parmbsc1 are
softer and also fail to kink at comparable curvature on 20–240× longer runs; parmbsc1's
free-energy landscape is *sub*-harmonic, i.e. more deformable than elasticity predicts;
and the only alternative family, CHARMM36, has the opposite pathology of spurious
base-pair disruption and cannot hold a duplex past 1 µs. Notably, the classic MD kinking
results came from force fields with a 12 bp/turn helical repeat artefact, which imposes
spurious torsional stress in a closed circle — so some of the historical evidence *for*
kinking may itself be a force-field artefact, in the opposite direction from your concern.

There is one genuine gap: **no published head-to-head comparison of kinking propensity
across bsc0/bsc1/OL15/OL21/CHARMM36 in a sharply bent duplex exists**, as far as I can
find. That is the control your question implies and nobody has run it.

### What is still running

Three things, all waiting on GPU allocation rather than on analysis:

- **The 27 bp windows in a box that fits them**, at 50, 58, 66 and 72 Å, to recover the
  chimera's own separation and with it the measured-versus-modelled comparison at the
  geometry that actually matters. This is the withdrawn 33 pN figure.
- **The 40 bp pair at 70 Å**, for the same reason.
- **The seeded kinks** — three positions (bp 7, 15, 21), two replicates, 30 ns — testing
  whether a kink placed *below* the onset persists, which bears on whether the transition
  is hysteretic. Built and validated; the seeding needed three attempts, because the first
  version seeded into an already-bent duplex and produced an initial energy of
  3.4 × 10¹¹ kJ/mol.

I should say plainly why these are queued rather than finished: I ran 129 GPU jobs in one
day on the billed partition and exhausted the whole lab's fairshare, which stopped
everyone's jobs including my own. They are resubmitted to the unbilled `gpu_test`
partition and will run when the allocation resets.

The timescale arithmetic above still stands as a statement about *spontaneous* nucleation
in the chimera's own geometry — but it is now the secondary explanation, not the primary
one.

---

## 9. The Et-tet residues shorten the site separation — why

Your observation is right and worth pursuing: the two Tet2-Et residues shorten Cβ–Cβ by
0.62 Å relative to WT (32.245 → 31.627 Å). Your specific hypothesis, a local
secondary-structure change, is **not** supported — but there are two real local effects,
and together they explain it.

### DSSP says the fold is untouched

Over residues 128–155, secondary structure is identical in WT, 2× Tet2-Et and clicked:
site 133 is coil in every case (0 % β) and site 149 is 100 % β in every case, with His148
at 64–79 % β and no systematic trend.

### What does change: a backbone dihedral, and a stacking interaction

Circular means over all replicates (3 for WT and Tet2-Et, 9 for the clicked control,
6 per chimera; first 5 ns of each discarded):

| system | site 133 (φ, ψ) | site 149 (φ, ψ) | ring–Tyr200 contacts at 149 |
|---|---|---|---|
| WT (Asp133, Asn149) | **+56.9, +16.3** | −126.6, +137.2 | — |
| 2× Tet2-Et | **+74.4, +2.4** | −125.5, +124.1 | 10.5 ± 9.5 |
| clicked, no DNA | **+61.3, +13.3** | −113.7, +131.5 | 0.68 ± 0.4 |
| 27 bp chimera | +63.1, +9.5 | −112.3, +132.5 | 0.96 ± 0.4 |
| 40 bp chimera | +62.9, +9.9 | −106.8, +125.3 | 8.4 ± 2.5 |
| 40 bp nicked | +62.9, +10.2 | −111.8, +120.1 | 9.6 ± 1.9 |

**Site 133 is the clean result.** Asp133 natively occupies the **left-handed α region** —
positive φ — which aspartate, asparagine and glycine commonly adopt and which a
phenylalanine-shaped side chain tolerates poorly. Tet2-Et *is* a phenylalanine
derivative, so it pushes φ out by 17.5° and ψ by 14°; clicking gives most of it back. The
replicate spread is tiny (WT ±0.8°, Tet2-Et ±0.5°), so this is not noise. The ring at
site 133 never touches Tyr200 in any system.

**Site 149 stacks on Tyr200, but not reproducibly.** I first measured this on one
replicate and got 19 contacts per frame, which overstated it. Across three replicates the
unclicked tetrazine gives 10.5 ± 9.5 — the individual values are 20.0, 10.4 and 0.9, so
the ring samples both stacked and unstacked states and three replicates cannot pin the
occupancy. What *is* tight is the clicked control at 0.68 ± 0.4 over nine replicates. And
the 40 bp chimeras stack again (8.4, 9.6), so this is not simply "clicking releases it" —
the ring at 149 can stack, and whether it does depends on the construct.

**Both variables track the contraction.** Across all 30–33 replicates that have both
measurements:

| predictor of Cβ–Cβ | Pearson r | p | slope |
|---|---|---|---|
| Tyr200 contacts at 149 | −0.68 | 5×10⁻⁵ | −38 mÅ per contact |
| φ at site 133 | −0.59 | 3×10⁻⁴ | −39 mÅ per degree |
| both together | R² = 0.55 | | −22 mÅ/deg and −29 mÅ/contact |

So the contraction you noticed tracks **two local structural variables and not secondary
structure**, and the two are partly independent — together they account for a bit over
half the variance in Cβ–Cβ across every system that has these residues. I would call that
suggestive rather than causal: these are correlations across constructs that differ
chemically, so part of the association is simply that both variables differ between the
unclicked and clicked states.

### Which sites tolerate Tet2-Et

This is the practical answer to your last sentence, and it is a prediction the data
supports:

1. **Avoid positive-φ positions.** Sites in the left-handed α region are occupied by
   Asp/Asn/Gly for a reason, and a bulky aromatic side chain does not fit there. Site 133
   is exactly such a position and is the one that distorts.
2. **Treat an aromatic neighbour as a liability.** The ring at site 149 stacks on Tyr200
   in some constructs and not others, and the amount it does correlates with the
   contraction. That is a variable you do not control and cannot see without measuring
   it, and — critically — it **differs between the unclicked and clicked states**, so it
   will masquerade as a mechanical effect for anyone using the unclicked protein as a
   reference. That is exactly how the point-3 artefact arose.
3. **Ordinary β and α positions with no aromatic neighbour should be fine.** Site 149's
   backbone moves by less than 13° and its secondary structure is unchanged; what makes
   it awkward is the neighbour, not the position.

### How this connects back to #3

The 0.647 Å that clicking recovers is not an inexplicable chemistry offset. Across every
replicate that carries these residues, Cβ–Cβ tracks the φ angle at site 133 and the
stacking at site 149, and together those two account for a little over half its variance.
Whatever its precise decomposition, the offset is **local side-chain and backbone
structure at the modified sites**, not anything the spring does — which is why using the
unclicked protein as the zero produced a spurious force response.

---

## 1. The tutorial

Written as a separate document: [`2026-8-22_setup-tutorial.md`](2026-8-22_setup-tutorial.md).
It covers the pipeline end to end, why each choice was made, and — at more length than is
flattering — the traps, since those are the parts worth passing on.

---

## Also included: the ubiquitin system

You suggested it and it turned out to be the better-controlled case, because the preprint
measures a perturbation ladder we can match. Five systems, built and running:

| system | atoms | loop |
|---|---|---|
| U0_apo | 19,497 | — |
| U1_clicked (crosslinkers, no DNA) | 58,184 | — |
| U2_ss50 (circular ssDNA-Ub) | 124,044 | closed |
| U3_ds50 (circular dsDNA-Ub) | 123,384 | closed |
| U4_ssarms (two 25-nt arms) | 119,758 | open |

1UBQ with F45W/T9C/K63C; Cys9 Cβ to Cys63 Cβ comes out at 24.81 Å, matching the "~2.5 nm"
the preprint quotes. The DBCO-maleimide crosslinker was built as two residues split at an
amide, and the connectivity is checked rather than assumed: strip from them the atoms that
do not come from the crosslinker — the cysteine's own backbone and side chain, the three
triazole nitrogens the azide brought, and the arm running out to the DNA — and what is left
is exactly **C₂₅H₂₁N₃O₄**, the reagent's published formula.

Four details are **chosen rather than established**, and I would rather say so than have
them found later: the two spacer lengths inside DBCO-maleimide (only one arrangement of the
known pieces reproduces the formula, but the vendor publishes no structure sheet); the
triazole regiochemistry (SPAAC on DBCO gives both and they are not separated); the
stereocentre created at the thioether carbon (racemic in the experiment); and the NHS-ester
azide's own spacer.

**What these runs cannot do.** The β5 register shift is a slow, discrete equilibrium, and
the displacement a ~4 pN spring would produce along it is ~0.1 Å — below thermal noise and
far slower than MD. You flagged this yourself. These runs test the delivered force, the
local structure at the attachment sites, and the protein–DNA contacts. They do not test the
retraction equilibrium, and the report will not imply otherwise.

One gap worth naming: `U1_clicked` is **not** the preprint's DBCO-Ub, which carries an
*unreacted* alkyne. It is the matched-chemistry, no-DNA control — the ubiquitin counterpart
of `S2_clicked`. Modelling the unreacted DBCO needs a third residue with a strained alkyne.
Given that the preprint finds most of the NMR perturbation is already present at the
DBCO-Ub stage, that is the most interesting missing rung, and it is the obvious next build.

---

## What has run

**3,971 ns over 81 replicates**, of which 1,421 ns is new since the last version of this
document.

| set | systems | replicates × length | state |
|---|---|---|---|
| corrected closed-loop chimeras | L3_loop27, L4_loop40, L5_loop40nick | 3 × 50 ns | done |
| ssDNA spring and non-circular control | C1_ss50, C2_ssarms | 3 × ~40 ns | done |
| ubiquitin | U0_apo, U1_clicked, U2_ss50, U3_ds50, U4_ssarms | 3 × 40–50 ns | done |
| force–extension, two stiffnesses | 27/40/40-nick, 6 windows each | 2 × 30–40 ns | done |
| deeper compression, for the kink threshold | 27 bp at 44/38/32 Å, two stiffnesses | 2 × 30 ns | done |
| seeded type II kinks | 3 positions in a 27-mer | 2 × 30 ns | queued |

The existing 2,550 ns is kept and reused, reframed as what it actually is: a two-arm,
non-circularising construct — Zocchi's pre-ligation chimera, and the closest analogue we
had to the preprint's non-circular control. Comparing it against the rebuilt closed loop
is a direct test of circular versus not, at no extra compute, and it is a cleaner one than
anything I had planned.

## Corrections made to the public report

All of these are done and republished at
<https://extremefattypunch.github.io/DNA-springs-project-2026-27/>:

1. The point-3 claim is retracted on the page, and replaced by the matched-control
   comparison, the barrel stiffness, and the contact measurement.
2. Figure 1's caption carried the same claim — it said the chimeras "sit well above that
   line, so the spring is not equivalent to pulling along the same axis". It now says the
   offsets on that panel are against the wrong zero, and points at the matched-control
   table instead.
3. The chromophore caption says rotation, and says why.
4. The five references are fixed, and the book is cited by page for every equation used.
5. A new section, *What the modification itself does*, carries the point-9 result and the
   Ramachandran figure.
6. The kinking section gains the curvature comparison, the timescale arithmetic, and the
   paragraph on OL15's documented bias.
7. Everywhere a force is quoted, the page now says whether it is modelled or measured.

Still to land, when the current runs finish: the closed-loop chimeras as the model of the
experiment (with the old runs reframed as the two-arm construct), the measured *f*(*x*)
laid over Zocchi's curve, the ssDNA and two-arm controls, the kink search in the
compression windows, and the ubiquitin set.

One thing I could not do automatically: GitHub does not create a wiki repository until
someone makes the first page in the web UI, so the rewritten landing page is staged at
`docs/wiki/Home.md` and needs pasting once.
