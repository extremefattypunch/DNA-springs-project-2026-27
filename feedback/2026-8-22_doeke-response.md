# Response to Doeke's comments of 22 August 2026

Ian Poon — sfGFP–DNA spring simulations
Written 23 August 2026. Numbers marked **[running]** come from simulations that were
still finishing when this was written; everything else is final.

---

## Summary

Nine comments. One found a genuine construction error, one forced me to retract a claim,
and one turned out to have a sharper answer than the hypothesis it came with.

| # | Comment | Outcome |
|---|---|---|
| 1 | Tutorial for lab meeting | Written: [`2026-8-22_setup-tutorial.md`](2026-8-22_setup-tutorial.md) |
| 2 | Two TCOs are on *different* strands | **You are right.** Real error. Rebuilt, and a topology test added that the old build fails |
| 3 | Simulate a non-circular control; the couple/shear argument is confused | **Claim retracted.** Your physics is right and mine was wrong. Rebuilt against a matched control; the effect inverts sign |
| 4 | Chromophore: bend or rotate? | **Rotate.** Caption reworded, with the measured torsions |
| 5 | The book DOI does not resolve | **Right, and four others were also wrong** — two of them resolved to unrelated papers |
| 6 | Is the force modelled or measured? | It was modelled. Now measured directly from the force field as well; both calculations written out in full |
| 7 | ssDNA chimera? | Built and running, 50 nt, plus a two-arm control that holds DNA content exactly fixed |
| 8 | "It does not kink" — force field, timescale, or model? | Mostly **model and timescale**. Literature reviewed; our curvature is 40–60 % below the threshold at which all-atom MD nucleates kinks |
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

**Status.** L3/L4/L5 built with 0 problems and integral charge; 3 replicates × 50 ns each
**[running]**.

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
3 replicates × 50 ns each **[running]**.

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
settled at. Here is that calculation written out, and then the direct measurement that
replaces the need to trust it.

### (a) The analytic route

Zocchi treats a sharply bent duplex as an elastic rod that, past a critical torque τ_c,
yields into a **kinked** state: two straight arms joined by a hinge that sustains a
constant torque. Bending modulus *B* = 200 pN nm² (*ℓ*_p ≈ 50 nm at 300 K), rise
0.33 nm/bp, τ_c = 27 pN nm nicked and 31–36 pN nm intact, k_BT = 4.11 pN nm.

With *L* the contour length and *x* the end-to-end distance,

```
γ    = L τ_c / 2B                     small parameter, valid for γ < 1
R    = (L/2)(1 − γ²/90)                arm half-length in the kinked state
x_c  = L (1 − γ²/15)                   the yield point
x₀   = L (1 − LT/10B)

x ≥ x_c   (unkinked, WLC-like):  f(x) = 10B/L² − T/(L−x)
x < x_c   (kinked):              f(x) = τ_c / (2R √(1 − (x/2R)²))
```

Substituting the numbers each chimera actually settles at:

| system | *L* | τ_c | γ | *x* measured | regime | *f* |
|---|---|---|---|---|---|---|
| S3_spring27 | 8.91 nm | 34 | 0.76 | 52.4 ± 5.9 Å | kinked | **4.82 ± 0.30 pN** |
| S4_spring40 | 13.20 nm | 34 | 1.12 | 62.5 ± 4.2 Å | kinked | **2.99 ± 0.06 pN** |
| S5_spring40nick | 13.20 nm | 27 | 0.89 | 60.8 ± 3.9 Å | kinked | **2.34 ± 0.04 pN** |

Two caveats that belong with these numbers and were not in the report. First, γ = 1.12 for
the 40-mer is outside the γ < 1 condition the derivation assumes (book p. 75), so that
row is an extrapolation. Second, the implementation is verified against Zocchi's own
published values before it is used here — 60-mer 1.47 vs 1.5 pN, 40-mer 2.33 vs 2.4 pN,
and *E*_DNA = 6.62 vs 6.6 k_BT for the GK chimera — so the arithmetic is right even where
the model's domain is stretched.

### (b) The direct route — the force field's own answer

Take the duplex alone, no protein, solvated. Hold its two ends at a fixed separation *x*
with a harmonic restraint of known stiffness *k*. At equilibrium the restraint's mean pull
balances whatever the molecule is pushing with, so

```
f_DNA(x) = k (⟨r⟩ − x₀)
```

with no elasticity model anywhere in the derivation. Scanning *x*₀ traces out *f*(*x*)
directly, and it can be laid over Zocchi's curve.

Three constructs (27 bp, 40 bp, 40 bp nicked) × six windows from 50 Å to nearly straight.
The restraint centre is walked from the built geometry to the window target over 2 ns
before any frame is recorded — dropping a duplex into a window 20 Å away does ~400 kJ/mol
of work in one step, and whatever base pair gave way under that impulse would look exactly
like the kink these runs are also meant to detect. That ramp is working: the 27-mer went
from 74.26 Å to 50.64 Å before recording began.

One thing the first pass taught me, which is worth recording because it is a real design
error and not a detail. At *k* = 200 kJ mol⁻¹ nm⁻² the restraint is 33 pN Å⁻¹, so a single
1 Å thermal fluctuation reads as 33 pN and the few-pN forces we want are invisible — a
4 pN spring produces a 0.12 Å offset against a 1.11 Å spread. The precision of a mean
force taken from equilibrium fluctuations is √(k k_BT / N) however the restraint is
applied, so a **softer** spring is strictly better for the same sampling. A second pass at
*k* = 33 kJ mol⁻¹ nm⁻² (5.5 pN Å⁻¹) resolves 2.5× finer and — because it widens each
window from 1.1 Å to 2.7 Å — is also the only one of the two whose windows overlap enough
to reconstruct a free-energy profile. Both passes are running; the stiff one is kept
because pinning *x* tightly is what the kink search wants.

**Status: [running]** — 36 stiff windows at 30 ns and 36 soft ones at 40 ns. The measured
*f*(*x*) against Zocchi's curve will go into the web report. Whatever it shows will be
reported, including disagreement.

---

## 7. The ssDNA chimera

Built and running: **C1_ss50**, one 50-nt strand clicked at both ends, no complement.

Two notes on the design.

**Why 50 nt and not 27.** At 0.255 nm per base a 27-mer has only ~6.9 nm of contour
against the 6.6 nm the sites plus linkers impose. It would be taut, and would report
stretching rather than the slack, zero-stress reference the control is for. 50 nt is also
the length used in the ubiquitin work, which makes the two studies comparable.

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

Mostly **(c) the model's threshold is too permissive** and **(b) 45 ns is far too short**.
(a), a force-field limitation, is real but secondary and does not point the way that would
rescue the model. Here is the evidence.

### What our system's curvature actually is

Radius 31–34 Å, i.e. 5.7–6.3° per base-pair step, equivalent to a 57–63 bp minicircle. One
correction to how the report framed this: that is **~25 % tighter in radius than
nucleosomal DNA** (the nucleosome superhelix gives ~42 Å and ~4.6°/bp), not "roughly
nucleosomal". Our measured per-step bend in the chimeras is 4.8–5.2° with the worst single
step at ~9.5°, and Watson–Crick opening averages 0.03–0.11 pairs per frame — essentially
none.

### (c) Zocchi's τ_c is the most kink-friendly estimate in the field

Converting each published criterion to a threshold radius via *R* = *B*/τ_c:

| source | criterion | threshold *R* | °/bp |
|---|---|---|---|
| Zocchi, nicked (EPL 2010) | τ_c = 27 pN nm | 76 Å | 2.6 |
| Zocchi, intact (PRX 2011) | τ_c = 31 pN nm | 66 Å | 3.0 |
| Vologodskii & Frank-Kamenetskii 2013; Du 2008 | nuclease probing | 35–40 Å | 4.9–5.6 |
| Curuksu et al. 2009 (umbrella sampling, AMBER) | asymptotic torque 40–51 pN nm | — | kink only above **8.6–9.3** |
| Fields, Meyer & Cohen 2013 (FRET "molecular vice") | WLC linear up to 7°/bp | ≤28 Å | 7 |
| **our chimeras** | | **31–34 Å** | **5.7–6.3** |

Zocchi's τ_c puts the softening transition at a *gentler* bend than a nucleosome. Every
independent estimate is 1.3–2.7× stiffer. On his criterion we are well past the
transition; on Curuksu's and on Fields'/Cohen's we are comfortably below it. So the
disagreement between our simulation and his prediction is substantially a disagreement
between his τ_c and the rest of the literature — and the one directly comparable
calculation, Curuksu's umbrella sampling on a bent *linear* duplex in the same force-field
family, finds kinks only above 8.6–9.3°/bp, **40–60 % more curvature than we impose**.
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

### (a) The force field, honestly

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

### What we are doing about it

The compression windows from #6 do double duty: they squeeze the duplex far past where it
sits in the chimera (down to 50 Å end-to-end on a 92 Å contour), and `force_extension.py`
measures the per-step bend profile and the N1–N3 distance of every Watson–Crick pair in
each window. A type-II kink breaks exactly one central pair, so a per-pair measure is the
right one. What counts as a kink is stated in advance: a local step far above the smooth
background, or an opened pair. **[running]**

The cheapest decisive experiment, which I have not run and propose to: **seed a type-II
kink and test whether it survives.** Build the three-base-pair motif with the central pair
broken and propeller ~120°, place it at several positions in the bent duplex, run 45 ns × 3.
If it is metastable, the force field can support kinks and the answer is (b) + (c). If it
heals in a few ns, that is a force-field result worth reporting on its own. It costs no
more than one existing system.

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

### What does change: one backbone dihedral, and one stacking interaction

| system | site 133 (φ, ψ) | site 149 (φ, ψ) | Tyr200 contacts at 149 |
|---|---|---|---|
| WT (Asp133, Asn149) | **+57.7, +12.4** | −126.6, +139.9 | — |
| 2× Tet2-Et | **+74.3, +4.8** | −134.4, +119.7 | **19.0 / frame** |
| clicked (TDP) | **+61.0, +14.2** | −115.1, +130.4 | **0.5 / frame** |

*(one replicate each, first 5 ns discarded; circular statistics)*

**Site 133.** Asp133 natively occupies the **left-handed α region** — positive φ — which
is a conformation aspartate, asparagine and glycine commonly adopt and which a
phenylalanine-shaped side chain tolerates poorly. Tet2-Et *is* a phenylalanine derivative,
so it pushes φ out by 17° and ψ by 8°. Clicking hands it back to essentially WT. Site 149
is ordinary β and barely moves.

**Site 149, and this is the larger effect.** The **unclicked tetrazine at 149 is pinned
flat against Tyr200** — 19 heavy-atom contacts per frame — and clicking releases it to 0.5.
The ring at site 133 never touches Tyr200 at all. So the unclicked intermediate is held
against a neighbour that the clicked adduct cannot reach, and that pulls the site in.

### This is also the answer to #3

The 0.647 Å that clicking recovers is not a mysterious chemistry offset: it is mostly
Tyr200 letting go. That matters because the retracted claim in #3 was built on using the
unclicked protein as the zero. The artefact had a mechanism, and the mechanism is a
stacking interaction of the unclicked intermediate — nothing to do with the spring.

### Which sites tolerate Tet2-Et

This is the practical answer to your last sentence, and it is a prediction the data
supports:

1. **Avoid positive-φ positions.** Sites in the left-handed α region are occupied by
   Asp/Asn/Gly for a reason, and a bulky aromatic side chain does not fit there. Site 133
   is exactly such a position and is the one that distorts.
2. **Avoid positions with an aromatic neighbour.** The unclicked tetrazine will stack on
   it, and the resulting distortion **disappears on clicking** — which means it will look
   like a spring effect to anyone using the unclicked protein as a reference. That is
   precisely how the point-3 artefact arose.
3. **Ordinary β and α positions with no aromatic neighbour are fine.** Site 149's backbone
   moves by less than 10° and its secondary structure is unchanged.

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

## What is running

| set | systems | replicates × length |
|---|---|---|
| corrected closed-loop chimeras | L3_loop27, L4_loop40, L5_loop40nick | 3 × 50 ns |
| ssDNA and non-circular control | C1_ss50, C2_ssarms | 3 × 50 ns |
| force–extension, stiff | 27/40/40-nick × 6 windows | 2 × 30 ns |
| force–extension, soft | 27/40/40-nick × 6 windows | 2 × 40 ns |
| ubiquitin | U0–U4 | 3 × 50 ns |

The existing 2,550 ns is kept and reused, reframed as what it actually is: a two-arm,
non-circularising construct — Zocchi's pre-ligation chimera, and the closest analogue we
had to the preprint's non-circular control. Comparing it against the rebuilt closed loop
is a direct test of circular versus not, at no extra compute, and it is a cleaner one than
anything I had planned.

## Corrections that will be made to the public report

1. Retract the point-3 claim and replace it with the matched-control comparison and the
   contact measurement.
2. State that the old runs are the two-arm construct, and present the closed-loop rebuild
   as the model of the experiment.
3. Reword the chromophore caption to torsion, with the measured numbers.
4. Fix the five references; cite the book by page.
5. Add the direct force measurement alongside the analytic one, with both written out.
6. Add the kinking discussion, with the curvature comparison table and the timescale
   arithmetic.
7. Add the Ramachandran and Tyr200 result, and the guidance on which sites tolerate
   Tet2-Et.
8. Say plainly, wherever a force is quoted, whether it is modelled or measured.
