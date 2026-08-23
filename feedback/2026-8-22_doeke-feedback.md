hi Ian, very interesting!. a few thoughts about your simulations:


1. When you present your first lab meeting, it would be great to include a short tutorial on how you set this up!.

2. If I am tracing the two strands of DNA correctly in your movie, the two TCO modifications are on two different strands of DNA. In our experiments one strand is modified on both ends while the complementary strand is unmodified. This prevents the chimera from forming dimers with other chimeras -- as Giovanni Zocchi has shown, such dimers are thermodynamically favorable.

3.

The spring moves the sites more than an equivalent axial pull
Ranked by the force their own spring delivers, the three chimeras sit +0.19, +0.52 and +0.59 Å from the zero-force clamp  monotonic in force, a slope of 0.13 Å pN{¹, more than ten times the axial-clamp slope. That is what Zocchi argues throughout his work8: forces applied at different places elicit different responses, not different magnitudes of one response. A spring anchored through two flexible tethers applies a couple and a shear as well as an extension, and the barrel evidently answers those more readily than a pure pull along C²C².
if you haven't yet, it would be helpful to also simulate the protein to two short ssDNA or dsDNA strands that do not form a circular chimera -- see https://www.biorxiv.org/content/10.64898/2026.07.05.736590v1 -- thinking about it, it may be interesting to simulate the ubiquitin system as well. More generally, I am a little confused about Claude's statement -- are the different simulations, force clamp and DNA attachment, not at the same sites? A couple (torque, I assume) would exert its effect indeed but (once averaged over replicates) only transiently until the forces at the attachment sites are equal and opposite. I am not sure what to make of the shear -- this may be dynamically important, but I don't see how it is important for the experienced forces when we are looking at equilibrium fluctuations of the system. (I could be wrong)

4. "Inside the barrel. The chromophore in green with the side chains lining its cavity, and dashed measures on the three hydrogen bonds that hold it planar" -- is the implication that the isolated chromophore would bend or that it would rotate?

5. Does Claude have access to Zocchi G. Molecular Machines: A Materials-Science Approach? The DOI it links to does not resolve.

6. "Zocchis f(x)8,9,10 evaluated at the end-to-end distance the simulation itself settles at" -- is the force calculated from a model of force exerted by DNA as a function of end-to-end distance or is this directly obtained from the force field calculations in the simulation, or by some other means? It would be helpful to understand the full calculation.

7. It is interesting to observe that the nicked DNA spring is weaker but not much weaker than the unnicked dsDNA spring. Have you run a simulation with an ssDNA chimera? -- here again it would be helpful to revisit the ubiquitin case (keeping in mind that the beta strand retraction equilibrium is too slow to cover using ordinary MD simulations).

8. "it does not kink" -- this is an important thing to look into. This can be a limitation of the MD force field, an MD timescale limitation (too short), or simply reflect errors in Giovanni's model of when exactly kinking will occur. It is worth scanning the literature more broadly to see how much is understood about DNA kinking dynamics.

9. simple observation: it is striking that inclusion of the two Et-tet residues reduces end-to-end distance -- perhaps it induces local changes in secondary structure at either/both sites? it may be good to see what the Ramachandran plots of these residues look like compared to canonical amino acids. Doing so may also tell us what kinds of mutations to Et-tet are the least invasive.