# Minimal cryo-EM mechanics feasibility test

## Goal

Perform a **small, low-compute cryo-EM experiment** to determine whether standard CryoSPARC analysis can recover a known force-associated conformational change in a protein.

This is a feasibility exercise only. Do **not** develop any new machine-learning method.

Use as little GPU time as reasonably possible so the experiment can run on a single low-priority/test GPU node such as Harvard's `gpu_test`.

The final output should answer:

1. Can a mixture of two mechanically different protein conformations be reconstructed from noisy cryo-EM projections?
2. Can CryoSPARC 3D Variability Analysis recover a coordinate that separates the two conformations?
3. Can CryoSPARC reconstruct intermediate/endpoint structures along that coordinate?
4. What information about the protein's mechanics can and cannot be obtained from this analysis alone?

---

# Biological system

Use F-actin structures from:

**Carl et al., Nature (2026), *Myosin forces remodel F-actin for mechanosensitive protein recognition*.**

Use these deposited maps:

* **EMD-46427** — no-ATP / low-force F-actin control
* **EMD-46426** — myosin-force-remodeled F-actin

These are experimentally reconstructed structures and therefore provide realistic mechanically distinct protein conformations.

Do **not** attempt to download or process tens of thousands of raw experimental movies.

For this test, generate a very small synthetic particle dataset directly from these two maps.

---

# Compute limits

Design everything around approximately:

* **1 GPU**
* preferably **<1 hour total GPU time**
* maximum roughly **2 hours**
* modest RAM
* no multi-GPU jobs
* no large movie processing
* no particle motion correction
* no high-resolution refinement
* no exhaustive parameter searches

Target only approximately **15–20 Å effective resolution**.

The purpose is to recover large-scale conformational mechanics, not atomic structure.

---

# Step 1 — Download and prepare two maps

Download:

```text
EMD-46427
EMD-46426
```

Call them:

```text
control.mrc
force.mrc
```

Prepare them identically.

Tasks:

1. Align the maps into the same reference frame.
2. Resample to the same voxel size.
3. Crop both tightly around a short F-actin segment.
4. Use the same box dimensions.
5. Normalize density amplitudes.
6. Apply the same soft mask.
7. Low-pass filter both to approximately:

```text
15 Å
```

There is no advantage to retaining high-frequency information for this test.

Prefer a small box such as approximately:

```text
96^3
```

or

```text
128^3
```

voxels.

Avoid boxes larger than necessary.

Generate:

```text
maps/
    control.mrc
    force.mrc
    difference.mrc
```

Also calculate:

* map-map correlation
* FSC between the two starting maps if convenient

The point is simply to verify that the two volumes are structurally different but well aligned.

---

# Step 2 — Generate a tiny synthetic cryo-EM dataset

Generate:

```text
2,500 control particles
2,500 force particles
```

for a total of:

```text
5,000 particles
```

This should be enough for the initial experiment.

If 5,000 turns out to be unnecessarily expensive, reduce to:

```text
1,500 + 1,500
```

before increasing GPU resources.

For every particle:

1. randomly choose a viewing orientation;
2. project the appropriate 3D volume;
3. apply a realistic but simple CTF;
4. add Gaussian noise;
5. optionally add very small random translations.

Use **exactly the same orientation, CTF and noise distributions for the two structural classes**.

The only systematic difference between them should be protein conformation.

Use a relatively easy SNR initially.

Do not make this a difficult benchmark.

Save:

```text
particles.mrcs
particles.star
ground_truth.csv
```

where `ground_truth.csv` contains:

```text
particle_id,true_state
1,control
2,force
...
```

The CryoSPARC workflow must **not** use `true_state`.

It is only for evaluating results afterward.

---

# Step 3 — Verify the simulation before using CryoSPARC

Before proceeding, make a small diagnostic panel showing approximately:

* 10 control projections
* 10 force projections

Confirm that the particles look like plausible noisy cryo-EM images.

Also generate class averages directly from the **known labels**, only as a simulation sanity check.

Do not spend more than a few minutes here.

---

# Step 4 — Import particles into CryoSPARC

Import the 5,000 shuffled particles.

Do not supply their ground-truth labels.

If orientations and CTF parameters can conveniently be supplied from the simulator, preserve them separately, but the preferred first test is to allow CryoSPARC to perform normal alignment.

Use the smallest practical particle box and Fourier crop further if appropriate.

---

# Step 5 — Produce a single consensus reconstruction

Run the minimum workflow needed to get an aligned particle set.

Suggested:

```text
Import Particles
        ↓
Ab-initio Reconstruction
        ↓
Homogeneous Refinement
```

Use:

```text
number of ab-initio classes = 1
```

Do not run:

* multiple ab-initio classes
* non-uniform refinement unless necessary
* local refinement
* CTF refinement
* per-particle motion correction
* high-resolution refinement

Terminate refinement once a sensible low-resolution F-actin volume is obtained.

The target is only approximately:

```text
15–20 Å
```

resolution.

Record:

```text
consensus.mrc
```

and the aligned particle set.

---

# Step 6 — Run a very small 3D Variability Analysis

Run **3D Variability Analysis** on the consensus-aligned particles.

Use:

```text
number of modes = 3
```

and approximately:

```text
filter resolution = 15–20 Å
```

Do not use a very high-resolution filter.

Use the same soft F-actin mask if required.

The main output needed is:

```text
per-particle coordinate along mode 1
per-particle coordinate along mode 2
per-particle coordinate along mode 3
```

---

# Step 7 — Compare 3DVA coordinates with the hidden true labels

After CryoSPARC finishes, join the per-particle 3DVA coordinates with:

```text
ground_truth.csv
```

For each variability component, produce one plot:

```text
3DVA coordinate
vs
true state
```

For example, violin plots or histograms of:

```text
control
force
```

along each coordinate.

Calculate for each component:

```text
ROC AUC
```

for discriminating `control` from `force`.

Also report:

```text
mean coordinate(control)
mean coordinate(force)
```

and optionally Cohen's $d$.

Identify:

```text
best_mode
```

as the mode with the highest force/control separation.

---

# Step 8 — Reconstruct states along the best mode

Using **3D Variability Display**, generate approximately:

```text
5
```

particle-supported reconstructions along `best_mode`.

For example:

```text
q1
q2
q3
q4
q5
```

Do not generate dozens of states.

Use `intermediates` or an equivalent mode that reconstructs from groups of actual particles.

Avoid interpreting purely interpolated density movies as experimental evidence.

Save:

```text
3dva_states/
    state_01.mrc
    state_02.mrc
    state_03.mrc
    state_04.mrc
    state_05.mrc
```

---

# Step 9 — Determine whether the structural trajectory corresponds to force

Compare each reconstructed state to both original maps:

```text
control.mrc
force.mrc
```

Calculate:

```text
corr(state_i, control)
corr(state_i, force)
```

or FSC if convenient.

Create a plot with variability coordinate on the $x$ axis and similarity on the $y$ axis.

The desired pattern is something like:

```text
low q → most similar to control
high q → most similar to force
```

or the reverse.

Also create a short movie of the five states if trivial to do.

---

# Step 10 — Measure one simple physical deformation

Do **not** attempt detailed atomic modelling.

Choose one very simple geometric observable that changes between the control and force maps.

Possible options:

* F-actin filament curvature
* filament radius
* inter-subunit axial spacing
* helical twist
* displacement of a recognizable density feature

Use whichever is easiest and most robust at approximately $15$ Å.

For each of the five reconstructed states measure:

$$
m(q)
$$

where $m$ is the chosen geometric quantity.

Plot:

$$
m
\quad\text{vs}\quad
q_{\mathrm{3DVA}}.
$$

The purpose is to determine whether the abstract CryoSPARC variability coordinate corresponds to a physically interpretable mechanical deformation.

Do not try to convert this to pN.

---

# Optional Step 11 — one continuous test

Only perform this if everything above works and there is still spare compute time.

Generate three interpolated structural states:

```text
control
middle
force
```

using map interpolation or experimentally available intermediate maps.

Generate:

```text
1,500 particles/state
```

for:

```text
4,500 particles total
```

Run the same workflow.

Assign:

```text
control = 0
middle = 0.5
force = 1
```

as hidden ground-truth coordinate:

$$
q_{\mathrm{true}}.
$$

Determine whether one 3DVA coordinate is monotonically correlated with this ordering.

Calculate:

$$
\rho_{\mathrm{Spearman}}
\left(
q_{\mathrm{3DVA}},
q_{\mathrm{true}}
\right).
$$

Do not do anything beyond this.

---

# Do not run these analyses

For this feasibility test, do **not** use:

* cryoDRGN
* RECOVAR
* 3DFlex
* neural-network training
* molecular dynamics
* atomistic refinement
* multiple GPU nodes
* raw movie processing
* CTF refinement
* Bayesian force inference
* thermodynamic free-energy estimation
* extensive hyperparameter sweeps

They are outside the scope of this experiment.

---

# Required final output

Create:

```text
minimal_force_cryoem_test/
├── README.md
├── maps/
├── particles/
├── cryosparc_outputs/
├── 3dva_states/
├── figures/
├── scripts/
└── results.csv
```

`results.csv` should contain at least:

```text
metric,value
map_correlation_control_force,...
AUC_mode1,...
AUC_mode2,...
AUC_mode3,...
best_mode,...
```

Include the continuous-test Spearman correlation only if that optional experiment was run.

---

# Required figures

Produce no more than **five figures**.

### Figure 1

Starting structures:

```text
control
force
difference
```

### Figure 2

Example noisy simulated particles.

### Figure 3

Distribution of the best 3DVA coordinate for:

```text
control
force
```

Include ROC AUC.

### Figure 4

Five reconstructed states along the best 3DVA coordinate.

### Figure 5

One mechanical structural observable plotted against the 3DVA coordinate.

---

# README questions

The final `README.md` should answer these questions directly.

## 1. Did standard CryoSPARC distinguish the mechanically different structures?

State:

* yes / partially / no
* best ROC AUC
* which 3DVA mode

## 2. Did the variability trajectory structurally connect the control and force maps?

Show whether map similarity changed monotonically from one endpoint toward the other.

## 3. Could the inferred coordinate be related to a physical structural deformation?

Describe the chosen geometric measurement and its relationship to the variability coordinate.

## 4. What does CryoSPARC actually provide?

Explain concretely what we can obtain from this experiment, for example:

```text
particle → latent variability coordinate
latent coordinate → reconstructed conformations
conformations → measurable structural deformation
```

## 5. What can it not provide from these data alone?

Explicitly discuss whether ordinary 3DVA provides:

* absolute force in pN;
* mechanical work;
* elastic energy;
* an intrinsic unbiased free-energy landscape;
* causal attribution of a structural mode to a known perturbation;
* prediction of the ensemble at an unseen force.

Do not speculate beyond what the experiment demonstrates.

---

# Success criterion

This test is successful if, using only approximately **5,000 particles and one small GPU**, CryoSPARC recovers a variability coordinate that substantially separates the known force/control conformations and whose reconstructed endpoints resemble the corresponding input structures.

As a rough diagnostic rather than a hard publication threshold:

$$
\mathrm{AUC}>0.8
$$

would indicate clear recovery of the mechanically associated structural variation.

For the optional three-state experiment,

$$
|\rho_{\mathrm{Spearman}}|>0.7
$$

would indicate that the inferred variability coordinate approximately recovers the known ordering of conformations.

If performance is poor, first check simulation/alignment quality rather than increasing computational resources.

---

# Priority

Complete the **binary control-versus-force experiment first**.

Do not spend time improving resolution once the mechanical conformational separation is visually and quantitatively clear.

The entire task should remain small enough to run comfortably as a test job on a single `gpu_test` GPU.
