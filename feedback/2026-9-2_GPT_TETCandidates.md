Yes. After combining **your Zocchi sources, your current sfGFP MD results, and the GFP structural literature**, I would **not choose 134/150 as the primary geometry if the explicit objective is to move the excitation and emission *wavelengths***.

The most promising region is the pair of β-strands immediately surrounding **His148 and Thr203**, because those residues directly control the chromophore. I would test **two spring geometries first: the N149/S202 pair and the S147/Q204 pair in canonical GFP numbering**. For your construct, the exact sequence numbers must be remapped before ordering because your repository already documents a numbering offset: your experimental D134/N150 corresponds to 2B3P D133/N149.

### My ranking

| Rank  | Spring anchors                 | Why                                                                                                             | Expected spectral effect                                           |
| ----- | ------------------------------ | --------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------ |
| **1** | **N149 ↔ S202** homologs       | Directly loads strands next to H148 and T203; lets you retain your existing N150-region handle                  | **Strong excitation-state shift; best practical first experiment** |
| **2** | **S147 ↔ Q204**                | Even more directly brackets H148/T203; proven spectroscopically sensitive geometry                              | **Potentially strongest excitation response**                      |
| **3** | **N149 ↔ Q204**                | Both surface-exposed; Q204 communicates with chromophore through structured water                               | Likely strong photophysical coupling, somewhat cleaner chemically  |
| **4** | D134 ↔ N149, your current pair | Mechanically active, but does not strongly perturb the most important chromophore restraints in your present MD | Better for intensity/proton-wire effects than large color shift    |

There is unusually strong experimental support for the first two geometries. **roGFP was engineered by putting a reversible disulfide across exactly these chromophore-adjacent strands.** S147C/Q204C gave about a **6.1-fold excitation-ratio dynamic range**, while N149C/S202C gave about **4.3-fold**. The structural reason is particularly relevant to you: residues 147 and 204 lie on adjacent β-strands next to H148 and T203. Constraining those strands alters chromophore protonation and redistributes the ~400- and ~490-nm excitation states. ([Springer][1])

That is almost a ready-made map of where to apply mechanical stress to GFP.

## Why I would move your 134 anchor to the ~202 region

Your present MD simulations already tell us something important.

For the current **D134/N150 spring**, the chromophore H-bond occupancies are approximately:

* H148 → chromophore: **72% WT**, spring-dependent but not systematically destroyed.
* T203 → chromophore: **~96–97% essentially throughout**.
* Arg96 → chromophore: **~100% throughout**.
* Ser205 → Glu222: **29% WT → 0% with the 40-bp spring**.

So your spring is clearly transmitting stress into the chromophore neighborhood, but the extraordinarily persistent **T203–chromophore and Arg96–chromophore interactions survive**.

That is almost exactly what I would expect for a geometry that changes **brightness/protonation/dynamics** but struggles to change the actual emitting electronic state.

Instead, consider:

$$
\boxed{\mathrm{N149 \longleftrightarrow S202}}
$$

These residues flank the pocket very close to

$$
\mathrm{H148 \longleftrightarrow CRO \longleftrightarrow T203}.
$$

The DNA spring then directly **pries apart/shears the two β-strands supporting the chromophore's two important phenolate-side interactions**.

The S147/Q204 design is even more aggressive:

$$
\boxed{\mathrm{S147 \longleftrightarrow Q204}}
$$

and has the strongest existing structural precedent: roGFP uses exactly S147C/Q204C, and closing the disulfide changes the 400/490-nm excitation equilibrium dramatically. ([PubMed Central (PMC)][2])

There is also independent evidence that **Q204 is a very attractive mechanical handle**. It is surface-exposed but sits close to the chromophore, and its backbone communicates with the chromophore through a conserved structured water molecule. Thus you can chemically attach something at the surface while mechanically coupling into the chromophore cavity. ([PubMed Central (PMC)][3])

### But note the direction of your DNA spring

This is important.

Although we casually say that the DNA spring "pulls on GFP," Zocchi's bent dsDNA spring actually acts as a **compressed leaf spring**. It wants its ends farther apart.

Therefore, approximately,

$$
\mathrm{site\ A}
\quad\longleftarrow\qquad\longrightarrow\quad
\mathrm{site\ B}.
$$

It **pushes the attachment sites apart**.

Zocchi describes exactly this: hybridization converts flexible ssDNA into a stiff dsDNA rod, producing a repulsive force between the two attachment points.  Your own implementation uses the same sign: the constant-force potential pushes the two attachment-site Cβ atoms apart.

So for S147/Q204 or N149/S202, you're essentially trying to **pry apart two adjacent chromophore-bearing β-strands**.

Interestingly, this is the *opposite sign* to the roGFP disulfide, which locks the corresponding strands together. That is not a problem. It actually makes a beautiful mechanistic experiment:

$$
\text{strand closure}
\leftrightarrow
\text{chromophore state A}
$$

versus

$$
\text{DNA-spring strand separation}
\leftrightarrow
\text{chromophore state B}.
$$

I would therefore predict a potentially large redistribution between the ~400 and ~490 nm excitation states, with the direction potentially opposite to oxidation of roGFP.

## How much force can you realistically apply?

Your sources give us a much firmer answer than the older approximate "10 pN" figure.

For Zocchi's experimentally calibrated geometries:

* **60-bp nicked spring:** ~**1.5 pN**
* **40-bp nicked spring:** ~**2.4 pN**
* intact springs can be somewhat stronger but are harder to calibrate exactly.

 

Zocchi's 2021 review explicitly describes the kink/softening transition as limiting the delivered force to **a few pN**. 

Your own current sfGFP simulations, after evaluating the Zocchi model at the actual relaxed spring geometry, span roughly

$$
\boxed{2.4\text{–}6.5\ {\rm pN}}
$$

across the spring designs, with the shorter 27-bp construction toward the high end. Your repository also correctly cautions against treating the **25–36 pN bare-duplex Euler-buckling force** as the expected experimental force—the latter corresponds to an unyielded rod rather than Zocchi's kinked/softened spring regime.

So I would design around roughly:

$$
\boxed{F_{\rm experimental}\sim2-7\ {\rm pN}}
$$

rather than 20–30 pN.

Your existing computational **0, 2, 4, 7 pN** clamp levels are therefore exactly the range I would use to screen alternative attachment geometries. The 12 and 20 pN simulations are useful stress tests but shouldn't be your assumed experimentally achievable spring regime.

One critical qualification: **a 27-bp spring will not necessarily give the new S147/Q204 pair the same force it gives D134/N150.** Force depends on the complete geometry:

$$
\text{protein site separation}
+
\text{linker lengths/orientations}
+
\text{DNA contour length}.
$$

Your own repository rejected 20 bp for the existing geometry because it became too short relative to the tether span.  So for each new pair, first determine the actual linker-to-linker end distance and then run the Zocchi spring calculation.

---

# There is a problem if you specifically want a **large emission-wavelength shift**

This is the most important caveat.

Changing H148 dramatically can make GFP switch from primarily

$$
\lambda_{\rm ex}\sim490\ {\rm nm}
$$

to

$$
\lambda_{\rm ex}\sim400\ {\rm nm},
$$

yet the protein can **still emit at ~510 nm**.

That has now been seen directly in sfGFP H148A/H148C. The mutations produce a dominant neutral ~400-nm ground-state population, but excitation at either ~400 or ~490 nm still leads to ~510-nm emission because excited-state proton transfer funnels the system into essentially the same emitting state. ([Nature][4])

So:

$$
\boxed{\text{big excitation shift}\not\Rightarrow\text{big emission shift}}
$$

This means even the perfect S147/Q204 DNA spring might give you a beautiful mechanical switch between 400 and 490 nm excitation **while emission remains ~510 nm**.

If your criterion is something like:

> "I want the spring to move emission by 10–20 nm, not just change CPS"

then I would change the protein design as well.

## My highest-risk/high-reward design: mechanically gate residue 203

Residue **203 is the strongest obvious lever for emission color**.

Changing T203 to Tyr creates YFP-like π-stacking between Tyr203 and the chromophore. Structural work shows that Tyr203 becomes nearly coplanar with the chromophore and even displaces it by roughly 1 Å. That electronic interaction shifts emission from green toward approximately **528 nm**. ([Wiley Online Library][5])

So instead of asking a WT sfGFP spring to somehow manufacture a new emitting state, I would build a mechanically sensitized construct:

$$
\boxed{
\text{T203Y sfGFP}
+
\text{spring across the 147/149 ↔ 202/204 region}
}
$$

The hypothesis would be:

**no/low stress**

$$
\text{Tyr203–chromophore π stacking}
\rightarrow
\lambda_{\rm em}\sim528\ {\rm nm}
$$

and then

**DNA spring pushes the surrounding β-strands apart**

$$
\rightarrow
\text{change Tyr203 position/orientation}
\rightarrow
\text{weaken/distort π stacking}
\rightarrow
\text{blue shift toward GFP-like emission}.
$$

If it worked strongly, that gives you the possibility of something qualitatively much more exciting:

$$
\boxed{\sim510\leftrightarrow528\ {\rm nm}}
$$

controlled by DNA hybridization.

That specific mechanical switching has **not** been demonstrated; it is my design inference. But unlike asking WT sfGFP to change emission color, it has a clear photophysical amplification mechanism: a few-pN protein deformation is being used to modulate a known chromophore π-stacking interaction whose spectral effect is already ~20 nm.

### What I would actually build

If I were prioritizing constructs for your project, I would make four:

1. **WT sfGFP + N149/S202-region spring** — best clean test of mechanical excitation-state switching.
2. **WT sfGFP + S147/Q204-region spring** — more aggressive chromophore coupling.
3. **T203Y sfGFP + N149/S202-region spring** — candidate mechanically tunable green↔yellow emission system.
4. **T203Y sfGFP + S147/Q204-region spring** — highest-risk/highest-potential geometry.

I would **not synthesize all four immediately**. Your existing MD framework can screen them first at **0, 2, 4 and 7 pN**.

For each candidate, rank the response by:

$$
d(\mathrm{H148-CRO}),
\quad
d/\theta(\mathrm{203-CRO}),
\quad
\mathrm{Ser205-Glu222},
$$

plus chromophore torsion/planarity, pocket water configuration, and ideally the **electric field projected across the chromophore**. For T203Y specifically, monitor the aromatic-ring/chromophore centroid distance and the interplanar angle.

Then do QM/MM or snapshot TD-DFT on the best trajectories. Classical MD can tell you **which geometry changes**, but it cannot reliably convert those changes into "the peak moves from 510 to 521 nm."

**So my single first-choice experimental redesign is: keep your existing N150-region attachment and move the 134 attachment to the S202-region homolog on the opposing chromophore-adjacent strand.** Then, if large *emission* shifts are the real objective rather than merely a 400↔490 excitation switch, repeat that architecture on a **T203Y sensitized sfGFP background**. This gives you a much stronger physical route from ~2–7 pN of DNA-spring force to a change in chromophore electronic structure.

[1]: https://link.springer.com/article/10.1007/s00109-013-1067-4?utm_source=chatgpt.com "Imaging ROS signaling in cells and animals | Journal of Molecular Medicine | Springer Nature Link"
[2]: https://pmc.ncbi.nlm.nih.gov/articles/PMC10611200/?utm_source=chatgpt.com "Genetically Encoded Biosensors for the Fluorescence Detection of O2 and Reactive O2 Species - PMC"
[3]: https://pmc.ncbi.nlm.nih.gov/articles/PMC7788595/?utm_source=chatgpt.com "Association of Fluorescent Protein Pairs and Its Significant Impact on Fluorescence and Energy Transfer - PMC"
[4]: https://www.nature.com/articles/s42004-025-01573-4?utm_source=chatgpt.com "Molecular dynamics guided identification of a brighter variant of superfolder Green Fluorescent Protein with increased photobleaching resistance | Communications Chemistry"
[5]: https://onlinelibrary.wiley.com/doi/10.1111/j.1365-2443.2009.01305.x?utm_source=chatgpt.com "Structural basis for red‐shifted emission of a GFP‐like protein from the marine copepod Chiridius poppei - Suto - 2009 - Genes to Cells - Wiley Online Library"
