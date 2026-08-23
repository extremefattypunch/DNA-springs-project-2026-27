# Building a protein–DNA chimera for all-atom MD: how this was set up

Ian Poon — notes for lab meeting, written 23 August 2026 in response to Doeke's comment 1.

The goal here is that someone in the lab could rebuild this from scratch, or build a
different chimera, without rediscovering the same traps. So this is organised as the
pipeline in order, with the reasoning for each choice, and then a section of things that
went wrong — which is the part that is actually worth an hour of anyone's time.

Everything lives in `MD simulations/sfGFP/`. The whole build is
`build/build_all.sh` → `build/build_loop_systems.sh` → `build/build_ubiquitin.sh`, and
each is meant to be readable top to bottom.

---

## 0. The thing we are building

A protein with two engineered attachment points, joined by a piece of double-stranded DNA
short enough that it has to bend to span them. The bent duplex pushes the attachment
points apart — that is the "spring" — and the question is what that does to the protein.

For sfGFP: Tet2-Et (a tetrazine-bearing unnatural amino acid) at positions 133 and 149,
clicked to strained-alkene-modified DNA. For ubiquitin, following the lab's preprint:
cysteines at 9 and 63, DBCO-maleimide, and azide-modified DNA.

The single most important structural fact, and the one we got wrong first time: **one DNA
strand is modified at both its 5′ and 3′ ends**, and the complement is added unmodified.
That is what makes the chimera a covalent circle. Put a linker on one end of each of two
different strands and you get something that looks identical in a movie but is held
together by base pairing alone.

---

## 1. Environment

The lab's `~/.bashrc` still points at the retired `/n/hekstra_lab/...`, so there is a
pinned environment with its own activation script:

```bash
source ./env/activate.sh     # exports DNASPRING_PY, DNASPRING_ENV, DNASPRING_SCRATCH
```

Three storage rules on FASRC, learned the hard way:

- `/n/lab_storage/hekstra_lab` is **full**. Code and small results only — never
  trajectories.
- Trajectories go to `/n/netscratch/hekstra_lab/Lab/ian_poon/`, which is large and purged
  after 90 days.
- Conda environments go on `/n/holylabs`, not in your home directory.

---

## 2. Preparing the protein

`build/01_protein/prepare_protein.py` (sfGFP), `prepare_ubiquitin.py` (ubiquitin).

Three decisions get made here that people usually let a tool make silently:

**Histidine tautomers.** HID and HIE are different molecules with different hydrogen-bond
patterns, and the default is a guess. We assign from geometry: whichever ring nitrogen has
a hydrogen-bond partner close by is the acceptor, so the *other* one carries the proton.
For sfGFP's His148 this matters — it is one of the three hydrogen bonds holding the
chromophore planar. It comes out HID (ND1 2.98 Å from the phenolate acceptor, NE2 2.96 Å
from an Arg168 donor). When the two distances are within 0.3 Å the script says it cannot
tell and defaults to HIE, rather than pretending.

**Mutations.** Truncate the side chain to backbone + Cβ, rename the residue, and let tleap
build the rest from its own template. Safe for Thr→Cys (gains an S and an H). A real
approximation for Phe→Trp, which is buried and gains a fused indole where a phenyl was —
that is flagged rather than hidden.

**Caps.** Only heavy atoms are passed to tleap; a partial set of crystallographic
hydrogens confuses it more than it helps.

---

## 3. Parameters for the residues that do not exist yet

`build/02_params/`. This is most of the work, and it is where a build quietly goes wrong.

There is **one source of truth** — `residue_defs.py` — which lists, for each custom
residue, every heavy atom, every hydrogen and its parent, every bond with its order, the
net charge, and which atom is the head, the tail, and the bond out to the next residue.
Nothing downstream is allowed to hard-code an atom name.

```
TET   Tet2-Et, the unclicked tetrazine amino acid
TDP   the click product: dihydropyridazine fused to the sTCO cage
DNL   5′-amino-modifier C6 arm            q = −0.3079
DNP   3′-amino-modifier arm, with its own phosphate   q = −0.6921
DNH   the same arm ending in a free alcohol: the clicked, no-DNA control
CMS   cysteine + maleimide thioether  (ubiquitin)
CDB   the DBCO cage + triazole        (ubiquitin)
```

`validate()` checks every declared atom's valence against its element and the formula
against the mass, before anything else runs.

Then, per residue:

```
model_builder.py   build a capped model compound and write mol2
                   (mol2, not SDF or PDB: antechamber needs bond orders AND atom names)
build_residues.py  antechamber (AM1-BCC, gaff2) → prepgen → parmchk2
                   + harmonise_backbone, verify_tree, assert_elements
```

Three guards worth copying into any similar project:

- **`verify_tree`** — every parent link in prepgen's output tree must be a bond that
  `residue_defs.py` declares. Without this you can get a residue with the right atoms, the
  right charge, and the wrong topology, which nothing else will notice.
- **`assert_elements`** — antechamber infers elements from atom *names*. Map each assigned
  atom type back to an element and compare with the tables.
- **`harmonise_backbone`** — an amino-acid residue parameterised by AM1-BCC has GAFF
  backbone charges, which do not match the ff14SB residues on either side of it. Overwrite
  the six backbone atoms with ff14SB types and charges and spread the residual over the
  side chain, so the peptide bonds are described consistently.

The target is **zero parameters flagged by parmchk2**. Anything flagged means parmchk2
guessed, and a guess should be a decision you made, not one you inherited.

### Charges are not optional bookkeeping

Amber splits a DNA strand's charge so the terminal residues carry −0.3079 (5′) and
−0.6921 (3′), summing with the internal −1's to give one negative charge per phosphate.
A custom linker that replaces a terminus has to take over that exact share. For a 27-mer
modified at both ends: −0.3079 + 27×(−1) + (−0.6921) = −28, which is the 28 phosphates it
actually has. The build asserts integral net charge, so getting this wrong fails loudly
instead of producing a system with a fractional charge and a plausible-looking trajectory.

### Junctions between force fields

Where a GAFF2 residue meets the Amber DNA backbone, *neither* force field has the
parameters — gaff2.dat does not know about Amber's `P`, and parm10.dat does not know about
GAFF's `os`. parmchk2 cannot invent them because the model compound was capped with
GAFF atoms on both sides. These go in `junction.frcmod`, and every value is copied
verbatim from Amber's own phosphodiester with a comment saying where it came from:

```
OS-P     230.0  1.610   →  os-P and OS-p5
CT-OS-P  100.0  120.50  →  c3-os-P and C7-OS-p5
```

(OL15 types the sugar C3′ as `C7`, not `CT` — check, do not assume.)

---

## 4. Building the DNA

`build/03_dna/build_dna.py`. PyMOL's `fnab` makes a straight B-form duplex; the script
then normalises the naming, verifies the geometry, and bends it.

**Verify before you bend.** The script checks the rise, the twist, the C1′–C1′ distances,
that the numbering runs 5′→3′, and that the two 5′ ends are at opposite ends of the helix.
Any of these failing means the builder gave you something other than what you asked for.

**Bending is subtler than it looks.** Base-pair centres spiral around the helical axis by
about 1.9 Å. Rotate each base pair about *its own centre* and drop that centre onto an arc,
and you silently replace every pair's off-axis offset with the first pair's — a lateral
shift of up to 3.8 Å that tears the backbone open **even at zero curvature**. The pivot has
to be the point where the base pair's plane crosses the helical axis. The test to insist on:
the transform must be the exact identity as *R* → ∞. Ours gives O3′–P = 1.601 Å at
R = 566 Å.

**Choosing the radius.** The span the protein imposes is between the two linker anchors —
the 5′ P of nucleotide 1 and the 3′ O3′ of nucleotide N, **both on the modified strand**.
Those two atoms are each ~9 Å off axis and (N−1)×34° apart in azimuth, so the span depends
strongly on which way the helix bends: at 27 bp it is 71 Å for one bend direction and 54 Å
for another *at the same radius*. So the solver scans azimuth and bisects on radius at
each, and requires a sign change before bisecting — otherwise a bracket that never
contained a root returns a confident wrong answer.

Because that span is a property of the protein and the linkers (Zocchi's *h*), we hold it
**fixed across spring lengths** and vary only the contour. That is what makes a force
ladder a controlled comparison.

---

## 5. Assembly

`build/04_assemble/`.

**tleap cannot build polycyclic side chains.** Handed the sTCO cyclopropane it returned
sides of 5.28, 4.70 and 1.52 Å. So the whole tether — protein-side residue through to the
DNA — is built as one molecule in RDKit, MMFF-optimised, and superposed onto the protein
backbone. `validate_cage()` then asserts the ring geometry survived: cyclopropane sides,
and planarity of each aromatic ring, with a tolerance **per ring** (a 4,5-dihydropyridazine
has two sp³ carbons and is a half-chair — holding it to an aromatic flatness rejects a
perfectly good structure).

**Posing the tether** is coordinate descent over its rotatable torsions, with several
random restarts, scoring clash against the protein plus the length of the bond about to be
formed. Score the **bond length**, not the distance to a fixed target point: aiming at a
point fixes the bond's direction as well as its length, and the 17-torsion 3′ arm was
landing at 2.15 Å that way. Scoring the bond directly gives 1.61 Å.

**Ring detection matters** in `rotatable()`. The traversal that decides which atoms move
when a bond is turned must be allowed to *reach* the other end by another route — that is
the ring test. An earlier version skipped the far atom explicitly, which made the test
unfalsifiable: every aromatic and cyclopropane bond reported as rotatable with the whole
molecule downstream of it.

---

## 6. tleap, and how it will surprise you

`build/04_assemble/build_system.py` drives it in two passes (solvate and count ions, then
rebuild with the right ion numbers). Force fields:

```
ff14SB · DNA.OL15 · TIP3P · Joung–Cheatham ions · Li–Merz 12-6 Mg²⁺
gaff2.dat · leaprc.xFPchromophores (for the GFP chromophore CRO)
```

Mg²⁺ is 12-6 rather than 12-6-4 because OpenMM has no C4 term — using 12-6-4 parameters
without the C4 term would be worse than using the set that was fitted without it.

Things tleap does that you have to know about:

- **It cannot handle paths with spaces.** Our project directory is `MD simulations/`, so
  every tleap input uses a working directory it can reach. This will bite you.
- **It applies 5′/3′ terminal templates to chain-initial and chain-final residues.** If a
  linker sits before nucleotide 1, that nucleotide is no longer chain-initial and takes an
  *internal* template — which is exactly the charge bookkeeping above.
- **It gives a prep residue an implicit tail at its last main-chain atom**, whether or not
  you declared `TAIL_NAME`. So two custom residues written back to back get bonded to each
  other automatically. In one of our arms that duplicated a bond we had declared and tleap
  died with `1-4: cannot add bond`; in the other it would have made a *wrong* bond
  silently. A `TER` between them turns the automatic chaining off.
- **`MAIN_CHAIN` lists the atoms *between* head and tail.** A list that starts at the head
  makes prepgen emit a tree in which the head is its own parent.
- **A prmtop has no chain field, and no crystallographic numbering.** tleap renumbers every
  residue sequentially from 1. So `residue_map.json` is written at build time recording
  each residue's original chain, number and name, and every analysis and figure resolves
  through it. Skipping this produces figures whose labels point at unrelated residues, and
  DNA analyses that silently do nothing because mdtraj read the whole solute as one chain.

### The assertions worth making at build time

- integral, zero net charge;
- every declared bond present in the finished prmtop;
- every expected residue present (and not a default list — a leftover default of `CRO`
  quietly failed every ubiquitin system, since that protein has no chromophore);
- and the topology test: **cut the Cα–Cβ bond at both modified residues and check the two
  Cβ are still connected.** "Is the solute one connected component?" does *not*
  discriminate — the wrong build passes it, because both arms hang off the protein
  backbone. What you need is a covalent *cycle*.

---

## 7. Simulation

`runs/openmm_engine.py`, `slurm/md_chain.sbatch`.

```
minimise → NVT with restraints → NPT with the restraints ramped off → production
LangevinMiddle, 4 fs with hydrogen mass repartitioning on the solute only
PME 0.9 nm, Monte Carlo barostat, 300 K
```

HMR is applied with ParmEd's `HMassRepartition` with `dowater=False` — repartitioning
water masses as well changes its dynamics for no benefit.

Equilibration uses **two contexts** (NVT, then a new one with the barostat added) rather
than toggling a barostat on and off, which is not something OpenMM supports cleanly. The
System object is reused rather than rebuilt: rebuilding costs ~75 s of ParmEd work per
call for nothing.

Everything is **restart-safe**: each stage is skipped if its output exists, production
checkpoints every ten report intervals and records progress in `progress.json`, and the
sbatch script uses `--requeue` with append-mode logs. On a preemptible partition you lose
at most one checkpoint interval.

**Put every physical condition in the run directory name.** We run six force-clamp
magnitudes and eighteen restraint windows off a handful of topologies. The first version
wrote all six clamp forces to `data/runs/S6_clamp` — six different physical conditions
overwriting each other's trajectory. `TAG="$SYS${CLAMP_PN:+_${CLAMP_PN}pN}..."` fixes it,
and it is worth doing before it costs you a day.

**Enter a biased window gently.** When a restraint holds two atoms at a separation far
from where the structure starts, dropping straight in does hundreds of kJ/mol of work in a
single step, and whatever bond gives way under that impulse looks exactly like the physics
you were hoping to observe. Ramp the restraint centre over a nanosecond or two first, and
do not write frames while you do.

---

## 8. Analysis

`analysis/`, one module per question, each writing a per-frame CSV and a summary JSON, all
gathered by `collect.py` into one table.

The measurement details that turned out to matter:

- **Circular statistics for anything angular.** A dihedral near ±180° wraps, and an
  arithmetic mean and standard deviation are meaningless — the chromophore's phenol torsion
  reported a 175° standard deviation that way. Report the circular mean, the circular sd,
  and the resultant length *R* (which tells you whether the mean means anything).
- **Hydrogen bonds by H···A distance *and* D–H···A angle**, not by heavy-atom distance.
  sfGFP's His148 sits 3.2–4.9 Å from the phenolate by heavy-atom distance, which a 3.5 Å
  cutoff calls broken 90 % of the time, while its HD1 is 2.3–2.5 Å away with near-linear
  geometry — a perfectly good, long, mobile hydrogen bond.
- **Smooth the helical axis over one full turn before measuring bend.** Raw
  centre-to-centre vectors are dominated by the ~1.9 Å spiral: measured that way an
  ordinary B-form duplex reads 27° per step and several hundred degrees total, which says
  nothing about bending.
- **Derive topology from the bond graph, not from assumptions.** The spring's end-to-end
  distance is measured between the nucleotide atoms actually bonded to a linker, found by
  walking the bonds. That reads correctly for both the two-arm and the closed-loop builds
  and records which one it found.
- **Use neighbour searches, not full pair matrices.** The naive protein–DNA contact
  calculation builds an 1800 × 1700 distance array per frame and gets killed by the OOM
  reaper.

---

## 9. The traps, collected

If you remember nothing else:

1. tleap cannot handle a path containing a space.
2. tleap gives a prep residue an implicit tail at its last main-chain atom.
3. `MAIN_CHAIN` lists interior atoms, not the head.
4. A prmtop has no chains and no crystallographic numbering — keep a residue map.
5. Amber's terminal templates change a residue's charge; the linker must take over the
   right share.
6. antechamber infers elements from atom names. Never name a carbon `CL`-anything. (And
   the reverse: a checker that maps GAFF's `na` to sodium will reject a correct residue —
   in a GAFF prepin every lowercase `n*` is nitrogen and only `cl` and `br` are two-letter
   elements.)
7. prepgen's tree does not reproduce ring closures; verify it against your own bond table.
8. tleap cannot build polycyclic side chains. Build them elsewhere and superpose.
9. Bend a duplex about the helical axis, not about base-pair centres; test at *R* → ∞.
10. Put every physical condition in the output path.
11. OpenMM writes hybrid-36 atom serials above 99,999 atoms, which ParmEd cannot parse —
    read the serialised State instead.
12. ParmEd's `save(renumber=False)` can write every atom serial as −1, which collapses
    CONECT records; assign `atom.number` explicitly and assert uniqueness.
13. "One connected component" is not a test for a covalent loop.

And the ones that cost me the most time, added after the first draft — every one of these
produced a *plausible* wrong number rather than an error:

14. **mdtraj minimum-images distances by default.** `compute_distances(..., periodic=True)`
    is the default and is right for a hydrogen bond and wrong for anything longer than
    half the cell's smallest perpendicular width — 38–48 Å in a truncated octahedron of
    the size we use. Our springs' ends are 60–75 Å apart, so every reported end-to-end
    distance was flipping between periodic images from frame to frame, and every force
    derived from one was wrong. The give-away was the spread: raw distances had a
    standard deviation of 3 Å and the minimum-imaged ones 8. OpenMM wraps by whole
    *molecule* and never splits one, so within a molecule the raw coordinates are already
    contiguous and are the honest answer; across molecules, beyond the minimum image
    there is no answer to be had and the code should refuse rather than guess.
15. **A hardcoded site pair silently produces nothing.** `collect.py` passed residues
    133/149 to every system; ubiquitin is modified at 9/63. It does not fail — it just
    finds no such residues and writes an empty column, which reads exactly like a
    system that has not run yet. Derive per-system identifiers from what the builder
    wrote, never from a constant.
16. **Prefix globbing collides.** `D_40bp*` also matches `D_40bp_nick_*`, which paired
    one duplex's topology with another's trajectory — 77,397 atoms against 77,407, and
    ParmEd raised a reshape error thirty lines from the cause. Match the suffix pattern
    you actually generate, and check the atom count as well.
17. **`pkill -f <pattern>` matches the process issuing it.** The bracket trick
    (`[a]nalysis/collect.py`) does not save you if the same command line also contains
    the plain string somewhere else — mine did, in the relaunch that followed the kill,
    so it killed its own shell and the edit it was about to make never happened. Twice.
18. **A report file that only describes its last invocation.** `PARAMS_REPORT.json` holds
    the residues of the run that wrote it, so building a subset silently narrows what
    downstream text can quote — which is how a verification table came to say "all four
    custom residues" long after there were seven. Either always build the full set, or
    derive the count and let it move.
19. **Seeding a distortion into an already-distorted structure.** An 85° kink rolled into
    a duplex that was already bent 122° folded the helix into itself: sixty atom pairs
    closer than 1.5 Å, a guanine ring left stretched to 2.66 Å after minimisation, and
    NaN coordinates six minutes into production on all six runs. Minimisation relaxes
    *into* a tangle and cannot escape it. A clash check before writing the structure
    costs a second and would have saved the six jobs.
20. **Prefer a stretched bond to a threaded ring.** When a rigid-body seed has to break
    something, break it in a way minimisation can undo. A covalent bond stretched to 3 Å
    pulls straight back in; two rings passed through each other never separate.

---

## 10. What I would do differently

- **Write the topology assertion first.** The attachment error would have been caught on
  day one by a test that takes fifteen lines.
- **Decide the control before running the experiment.** The retracted claim in the report
  came from comparing chimeras against the *unclicked* protein when a matched clicked
  control existed. The chemistry difference was 15–30× the effect being claimed.
- **Be explicit about modelled versus measured.** Quoting a force from an analytic model
  is fine; not saying that is what it is, is not.
- **Say which parts of a structure are chosen rather than known.** For the ubiquitin
  crosslinker, four structural details are decisions I made because the vendor publishes a
  formula but not a structure sheet. Those belong in the write-up, not in a comment.
- **Check the units and conventions of every library call that returns a number you will
  publish.** The periodic-boundary bug above was one default argument, it changed every
  force in the project, and nothing anywhere failed. The pattern to watch for is a
  quantity whose *scatter* is larger than the physics should allow — that is usually the
  first symptom, and it is visible long before the mean looks wrong.
- **Prefer a check that discriminates to a check that passes.** Several assertions here
  were satisfied by the wrong build as happily as the right one: one connected component,
  integral net charge, zero flagged parameters. The test that mattered — cut two bonds and
  ask whether a cycle remains — took fifteen lines and would have caught the original
  error on day one.
