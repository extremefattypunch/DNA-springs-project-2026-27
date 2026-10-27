"""Clip 05 - Heterogeneity: consensus blur, illustrative PCA and soft 3D classes.

Storyboard beats (spec section 3, Scene 05) -> code sections:
    beat 1 (0-7 s)      title + question; two synthetic conformers (appendage at site A / B)
    beat 2 (7-15 s)     example particles -> one (idealised) consensus map (smeared appendage)
    beat 3 (15-24.5 s)  illustrative PC1/PC2 feature plot, raw -> low-pass features
    beat 4 (24.5-38 s)  EM soft assignment (ambiguous points) -> class maps vs consensus
    beat 5 (38-45 s)    what Figure S6 of Day et al. (2025) actually did
    freeze (45-50.5 s)  takeaway card under the consensus / class maps (4.5 s static)

Every number/curve/image comes from cached("classification", dm.classification)
or dm.cached_render(...). The class maps are rendered as the mixtures
w_c*rho_A + (1-w_c)*rho_B, with w_c the responsibility-weighted state-A fraction
of each class (classification()['class_wA']); the consensus uses
w = consensus_wA (the true state-A fraction of all particles).

Colours: state A / site A = green (TRUTH; a light tint marks the appendage on
a map), state B / site B = purple (ALT), estimated maps = orange (EST).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *

from helpers import demos as dm
from helpers import styling as st
from helpers import synthetic_phantom as sp
from helpers.cache import cached

# the view of the particle images in dm.classification(); both appendage sites
# are clearly separated from the scaffold in it
VIEW = dict(rot=-35.0, tilt=25.0, psi=10.0)
SITE_A = "#C4F5CF"          # light tint of TRUTH: appendage density at site A
SITE_B = st.ALT             # appendage density at site B


def site_densities(n: int, voxel: float, ws: list) -> dict:
    """Mean density inside the site-A and site-B appendage masks of the mixture
    w*rho_A + (1-w)*rho_B, for each w in ws (site A reproduces the helper's
    contrast_* numbers; site B is the same measurement at the other site)."""
    vA = sp.conformer("A", n, voxel, atoms=False)
    vB = sp.conformer("B", n, voxel, atoms=False)
    mA = sp.appendage_mask("A", n, voxel)
    mB = sp.appendage_mask("B", n, voxel)
    a, b = [], []
    for w in ws:
        v = w * vA + (1 - w) * vB
        a.append(float(v[mA].mean()))
        b.append(float(v[mB].mean()))
    return dict(siteA=np.array(a), siteB=np.array(b))


def render(state="A", wA=None):
    return dm.cached_render(VIEW["rot"], VIEW["tilt"], VIEW["psi"], state=state, wA=wA)


def _hex(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)], float)


class MapPanel:
    """Shaded isodensity render with the appendage pixels tinted by site and
    the transparent border cropped away (so `height` is the molecule's height).

    rA / rB: renders of the same density labelled with state A's / state B's
    blobs (blob 15 = the appendage at site A / site B); either may be None.
    footprint: optional {"A": mask, "B": mask} of the pure-state appendage
    pixels; tints are restricted to them so that a mixture's label map (whose
    argmax only knows one state's blobs) cannot tint the wrong site.
    """

    def __init__(self, rA, rB, color, height, crop, footprint=None):
        base_r = rA if rA is not None else rB
        shade, mask = base_r["shade"], base_r["mask"]
        rgb = np.broadcast_to(_hex(color), shade.shape + (3,)).copy()
        taken = np.zeros_like(mask)
        for key, r, col in (("A", rA, SITE_A), ("B", rB, SITE_B)):
            if r is not None:
                m = (r["label"] == 15) & mask & ~taken
                if footprint is not None:
                    m &= footprint[key]
                rgb[m] = _hex(col)
                taken |= m
        rgb = np.clip(rgb * shade[..., None], 0, 255)
        rgba = np.concatenate([rgb, 255 * mask[..., None].astype(float)], -1).astype(np.uint8)
        (i0, i1), (j0, j1) = crop
        self.crop = crop
        rgba = rgba[i0:i1, j0:j1]
        im = ImageMobject(np.stack([st.to_display(rgba[..., c]) for c in range(4)], -1))
        im.height = height
        self.image = im

    def site_point(self, mask_full):
        """Scene point of the centroid of a full-size pixel mask on this panel."""
        (i0, i1), (j0, j1) = self.crop
        ii, jj = np.nonzero(mask_full)
        fx = (ii.mean() - i0 + 0.5) / (i1 - i0) - 0.5
        fy = (jj.mean() - j0 + 0.5) / (j1 - j0) - 0.5
        return self.image.get_center() + np.array([fx * self.image.width,
                                                   fy * self.image.height, 0.0])


def dashed_ring(center, radius, color):
    return DashedVMobject(Circle(radius=radius, color=color, stroke_width=3),
                          num_dashes=18).move_to(center)


class Scene05Classification(Scene):
    def construct(self):
        d = cached("classification", dm.classification)
        states = d["states"].astype(int)
        wcons = float(d["consensus_wA"])
        cwA = d["class_wA"].astype(float)
        sites = cached("s05_heterogeneity_sites", site_densities, n=64,
                       voxel=float(d["voxel"]),
                       ws=[1.0, 0.0, wcons, float(cwA[0]), float(cwA[1])])
        rA, rB = render("A"), render("B")
        mixes = [wcons, float(cwA[0]), float(cwA[1])]            # consensus, class 1, class 2
        mix_r = [(render("A", w), render("B", w)) for w in mixes]
        # common crop: union of all molecule masks plus a small margin
        allm = np.any([r["mask"] for r in [rA, rB] + [x for p in mix_r for x in p]], 0)
        ii, jj = np.nonzero(allm)
        n = allm.shape[0]
        crop = ((max(ii.min() - 4, 0), min(ii.max() + 5, n)),
                (max(jj.min() - 4, 0), min(jj.max() + 5, n)))
        siteA_mask = rA["label"] == 15                                # true appendage sites
        siteB_mask = rB["label"] == 15
        col_state = [st.TRUTH, st.ALT]                                # known synthetic state

        badge = st.simulation_badge()

        # ------------------------------------------------ beat 1: two conformers
        title = st.title_card("Heterogeneity and 3D classification")
        question = st.T("Why does a single consensus 3D map become blurred\n"
                        "when the particles come from different structural states?",
                        st.BODY_SIZE - 2, st.MUTED)
        question.next_to(title, DOWN, aligned_edge=LEFT, buff=0.18)
        self.play(FadeIn(title), FadeIn(badge), run_time=1.0)
        self.play(FadeIn(question), run_time=1.0)

        H1 = 3.3
        pA = MapPanel(rA, None, st.TRUTH, H1, crop)
        pB = MapPanel(None, rB, st.TRUTH, H1, crop)
        pA.image.move_to([-3.2, -0.35, 0])
        pB.image.move_to([3.2, -0.35, 0])
        labA = st.T("conformer A", st.BODY_SIZE, st.TRUTH).next_to(pA.image, UP, buff=0.2)
        labB = st.T("conformer B", st.BODY_SIZE, st.ALT).next_to(pB.image, UP, buff=0.2)
        disc = st.phantom_disclaimer().to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(pA.image), FadeIn(pB.image), FadeIn(labA), FadeIn(labB), FadeIn(disc),
                  run_time=1.5)
        ringA = dashed_ring(pA.site_point(siteA_mask), 0.5, SITE_A)
        ringB = dashed_ring(pB.site_point(siteB_mask), 0.5, SITE_B)
        tA = st.T("appendage at site A", st.BODY_SIZE - 4, SITE_A).next_to(pA.image, DOWN, buff=0.2)
        tB = st.T("appendage at site B", st.BODY_SIZE - 4, SITE_B).next_to(pB.image, DOWN, buff=0.2)
        self.play(Create(ringA), Create(ringB), FadeIn(tA), FadeIn(tB), run_time=1.2)
        same = st.T("same 15-domain scaffold; only the appendage moves", st.BODY_SIZE - 2)
        same.move_to([0, -2.9, 0])
        self.play(FadeIn(same), run_time=0.8)
        self.wait(1.6)                                                   # ~9 s

        # ------------------------------------------------ beat 2: consensus blur
        self.play(FadeOut(VGroup(question, labA, labB, ringA, ringB, tA, tB, same, disc)),
                  FadeOut(pA.image), FadeOut(pB.image), run_time=0.8)
        tiles = Group()
        for k in range(6):
            im = st.array_to_image_mobject(d["examples"][k], height=1.4, nearest=True)
            fr = Rectangle(width=1.4, height=1.4, stroke_color=col_state[d["example_states"][k]],
                           stroke_width=4).move_to(im)
            tiles.add(Group(im, fr))
        tiles.arrange_in_grid(rows=2, cols=3, buff=0.18).move_to([-4.3, 0.35, 0])
        tl = st.T("example particles (one view shown)", st.BODY_SIZE - 4)
        tl.next_to(tiles, UP, aligned_edge=LEFT, buff=0.2)
        leg = st.legend([(st.TRUTH, "state A (known)"), (st.ALT, "state B (known)")],
                        size=st.BODY_SIZE - 4)
        leg.next_to(tiles, DOWN, buff=0.25)
        self.play(LaggedStart(*[FadeIn(t) for t in tiles], lag_ratio=0.15), FadeIn(tl),
                  FadeIn(leg), run_time=1.5)

        fp = {"A": siteA_mask, "B": siteB_mask}
        pC = MapPanel(*mix_r[0], st.EST, 3.2, crop, fp)
        pC.image.move_to([2.9, 0.6, 0])
        labC = st.T("consensus map (all particles)", st.BODY_SIZE - 2, st.EST)
        labC.next_to(pC.image, UP, buff=0.2)
        arr = Arrow(tiles.get_right() + 0.1 * RIGHT, pC.image.get_left() + 0.1 * LEFT,
                    color=st.MUTED, stroke_width=4, buff=0.05)
        avg = st.T("reconstruct\n(all views)", st.BODY_SIZE - 4, st.MUTED).next_to(arr, UP, buff=0.1)
        self.play(GrowArrow(arr), FadeIn(avg), run_time=0.8)
        self.play(FadeIn(pC.image), FadeIn(labC), run_time=1.2)
        rgA = dashed_ring(pC.site_point(siteA_mask), 0.42, SITE_A)
        rgB = dashed_ring(pC.site_point(siteB_mask), 0.42, SITE_B)
        eq_c = st.M(r"\bar\rho=%.3f\,\rho_A+%.3f\,\rho_B" % (wcons, 1 - wcons), size=42)
        eq_c.next_to(pC.image, DOWN, buff=0.25)
        idl = st.T("idealised: noise-free mixture\nat the true state-A fraction", st.SMALL_SIZE,
                   st.MUTED).next_to(eq_c, DOWN, buff=0.12)
        self.play(Create(rgA), Create(rgB), FadeIn(eq_c), FadeIn(idl), run_time=1.0)
        stmt = st.caption("A consensus reconstruction can average distinct states:\n"
                          "both appendage sites come out weak and shrunken.", st.TEXT)
        self.play(FadeIn(stmt), run_time=0.8)
        self.wait(1.7)                                                   # ~17 s

        # ------------------------------------------------ beat 3: PCA features
        self.play(FadeOut(Group(tiles, tl, leg, arr, avg, pC.image, labC, rgA, rgB, eq_c, idl, stmt)),
                  run_time=0.8)
        ax = Axes(x_range=[-4.5, 3.5, 1], y_range=[-3.5, 3.5, 1], x_length=5.2, y_length=4.8,
                  tips=False, axis_config=dict(color=st.MUTED, stroke_width=2,
                                               include_ticks=False))
        ax.move_to([-3.6, -0.65, 0])
        box = SurroundingRectangle(ax, color=st.MUTED, stroke_width=1, buff=0.0)
        xl = st.T("PC1", st.BODY_SIZE - 4, st.MUTED).next_to(ax.x_axis.get_right(), UL, buff=0.1)
        yl = st.T("PC2", st.BODY_SIZE - 4, st.MUTED).rotate(PI / 2)
        yl.next_to(box, LEFT, buff=0.08)
        head = st.T("PC1/PC2 - illustrative features", st.BODY_SIZE - 2)
        head.next_to(box, UP, buff=0.15)

        def dots_for(P):
            return VGroup(*[Dot(ax.c2p(*P[i]), radius=0.045, color=col_state[states[i]],
                                fill_opacity=0.8) for i in range(len(P))])

        dots = dots_for(d["pcs_raw"])
        dots_f = dots_for(d["pcs"])
        self.play(Create(ax), FadeIn(box), FadeIn(xl), FadeIn(yl), FadeIn(head), run_time=1.0)

        rx = -0.45
        leg2 = st.legend([(st.TRUTH, "known state A (synthetic)"),
                          (st.ALT, "known state B (synthetic)")], size=st.BODY_SIZE - 4)
        leg2.move_to([rx, 1.75, 0], aligned_edge=LEFT)
        f0 = st.T("features: raw pixels", st.BODY_SIZE).move_to([rx, 0.35, 0], aligned_edge=LEFT)
        f1 = st.T("features: low-pass-filtered pixels", st.BODY_SIZE)
        f1.move_to([rx, 0.35, 0], aligned_edge=LEFT)
        n0 = st.T("one dot = one particle image", st.BODY_SIZE - 4, st.MUTED)
        n0.next_to(f0, DOWN, aligned_edge=LEFT, buff=0.25)
        self.play(LaggedStart(*[FadeIn(dt) for dt in dots], lag_ratio=0.004), FadeIn(leg2),
                  FadeIn(f0), FadeIn(n0), run_time=1.6)
        ov = st.T("the two clouds overlap", st.BODY_SIZE, st.TEXT)
        ov.next_to(n0, DOWN, aligned_edge=LEFT, buff=0.4)
        self.play(FadeIn(ov), run_time=0.5)
        self.wait(1.0)
        sep = st.T("clouds become distinguishable\n(some overlap remains)", st.BODY_SIZE, st.TEXT)
        sep.next_to(n0, DOWN, aligned_edge=LEFT, buff=0.4)
        self.play(Transform(dots, dots_f), FadeTransform(f0, f1), FadeTransform(ov, sep),
                  run_time=2.4)
        g1 = st.caption("PCA: linear image coordinates, not a mechanism or an energy landscape")
        self.play(FadeIn(g1), run_time=0.6)
        self.wait(1.6)                                                   # ~28 s

        # ------------------------------------------------ beat 4: soft assignment
        self.play(FadeOut(VGroup(f1, n0, sep, g1)), run_time=0.5)
        eq_mix = st.M(r"p(I_i)=\sum_{c=1}^{K}\pi_c\,p(I_i\mid\rho_c)", size=40)
        eq_sum = st.M(r"\sum_c\pi_c=1", size=40)
        # the demo's responsibilities come from a 2-D Gaussian mixture on the PCA
        # features x_i, so gamma is written in that (feature-space) form
        eq_g = st.M(r"\gamma_{ic}=\frac{\pi_c\,\mathcal N(x_i\mid\mu_c,\Sigma_c)}"
                    r"{\sum_{c'}\pi_{c'}\,\mathcal N(x_i\mid\mu_{c'},\Sigma_{c'})}", size=40)
        eq_mix.move_to([rx, 0.6, 0], aligned_edge=LEFT)
        eq_sum.next_to(eq_mix, RIGHT, buff=0.4)
        eq_g.next_to(eq_mix, DOWN, aligned_edge=LEFT, buff=0.22)
        self.play(FadeIn(eq_mix), FadeIn(eq_sum), run_time=0.9)

        def ellipse(mu, cov, k=2.0):
            w, V = np.linalg.eigh(cov)
            L = V @ np.diag(np.sqrt(w)) * k
            return ParametricFunction(lambda t: ax.c2p(*(mu + L @ np.array([np.cos(t), np.sin(t)]))),
                                      t_range=[0, TAU], color=st.EST, stroke_width=4)
        mu = d["gmm_mu"]
        ells = VGroup(*[ellipse(mu[c], d["gmm_cov"][c]) for c in range(2)])
        cl1 = st.T("class 1", st.BODY_SIZE - 4, st.EST)
        cl1.next_to(ells[0].get_top(), UP, buff=0.08, aligned_edge=RIGHT).shift(0.1 * RIGHT)
        cl2 = st.T("class 2", st.BODY_SIZE - 4, st.EST)
        cl2.next_to(ells[1].get_top(), UP, buff=0.08, aligned_edge=LEFT).shift(0.1 * LEFT)
        self.play(Create(ells), FadeIn(cl1), FadeIn(cl2), FadeIn(eq_g), run_time=1.2)

        amb = d["ambiguous"].astype(int)
        rings = VGroup(*[Circle(radius=0.1, color=st.TEXT, stroke_width=2.5).move_to(
            ax.c2p(*d["pcs"][i])) for i in amb])
        notes = VGroup()
        for i, y in ((278, 3.0), (319, -3.05)):
            g = st.T(f"P(class 1) = {d['resp'][i]:.2f}", st.BODY_SIZE - 4, st.TEXT)
            g.add_background_rectangle(color=st.BG, opacity=0.85)
            g.move_to(ax.c2p(3.35, y), aligned_edge=RIGHT)
            p = ax.c2p(*d["pcs"][i])
            end = g.get_corner(DL if y > 0 else UL) + 0.2 * RIGHT
            notes.add(g, Line(p + 0.1 * normalize(end - p), end, color=st.TEXT, stroke_width=2))
        amb_t = st.T(f"{len(amb)} of {len(states)} stay ambiguous (0.2 < P < 0.8)",
                     st.BODY_SIZE - 4, st.TEXT)
        amb_t.next_to(eq_g, DOWN, aligned_edge=LEFT, buff=0.22)
        pi_t = st.T(f"fitted π = {d['gmm_pi'][0]:.2f}, {d['gmm_pi'][1]:.2f}: computational\n"
                    "assignments, not equilibrium populations", st.BODY_SIZE - 4, st.MUTED)
        pi_t.next_to(amb_t, DOWN, aligned_edge=LEFT, buff=0.15)
        self.play(Create(rings), FadeIn(notes), FadeIn(amb_t), run_time=1.0)
        self.play(FadeIn(pi_t), run_time=0.6)
        g2 = st.caption("demo: K = 2 Gaussians fitted by EM on these 2 features (soft, not k-means)")
        self.play(FadeIn(g2), run_time=0.5)
        self.wait(3.0)                                                   # ~33 s

        # class maps vs consensus
        self.play(FadeOut(VGroup(ax, box, xl, yl, head, dots, leg2, eq_mix, eq_sum, eq_g, ells,
                                 cl1, cl2, rings, notes, amb_t, pi_t, g2)), run_time=0.7)
        H3 = 2.5
        xs = [-4.45, 0.0, 4.45]
        panels = [MapPanel(*mix_r[k], st.EST, H3, crop, fp) for k in range(3)]
        maps = Group(*[p.image for p in panels])
        names = VGroup(*[st.T(nm, st.BODY_SIZE, st.EST) for nm in ("consensus", "class 1", "class 2")])
        pcts = VGroup(*[st.T((f"{100 * f:.0f}%" if f >= 0.01 else "<1%") + " state A",
                             st.BODY_SIZE - 4, st.EST) for f in mixes])
        top_y = title.get_bottom()[1] - 0.22
        for m, nm, pc, x in zip(maps, names, pcts, xs):
            nm.move_to([x, top_y - 0.3, 0], aligned_edge=DOWN)       # common baseline
            pc.move_to([x, top_y - 0.38, 0], aligned_edge=UP)
        y_map = min(pc.get_bottom()[1] for pc in pcts) - 0.12
        for m, x in zip(maps, xs):
            m.move_to([x, y_map, 0], aligned_edge=UP)
        site_rings = VGroup()
        for p in panels:
            site_rings.add(dashed_ring(p.site_point(siteA_mask), 0.3, SITE_A),
                           dashed_ring(p.site_point(siteB_mask), 0.3, SITE_B))
        # appendage density bars (mean density inside each site's mask)
        # one reference per row: the pure state that has its appendage at that site
        pure = {"siteA": float(sites["siteA"][0]), "siteB": float(sites["siteB"][1])}
        BL = 2.3
        bars = VGroup()
        for k, x in enumerate(xs):
            for row, (key, col, nm) in enumerate((("siteA", SITE_A, "site A"),
                                                  ("siteB", SITE_B, "site B"))):
                v = float(sites[key][2 + k])
                y = maps[0].get_bottom()[1] - 0.45 - 0.5 * row
                x0 = x - 0.85
                lab_ = st.T(nm, st.BODY_SIZE - 4, col).move_to([x0 - 0.12, y, 0], aligned_edge=RIGHT)
                track = Line([x0, y, 0], [x0 + BL, y, 0], color=st.PANEL, stroke_width=16)
                bar = Rectangle(width=max(BL * v / pure[key], 0.03), height=0.26, fill_color=col,
                                fill_opacity=0.95, stroke_width=0)
                bar.move_to([x0, y, 0], aligned_edge=LEFT)
                ref = DashedLine([x0 + BL, y - 0.19, 0], [x0 + BL, y + 0.19, 0], color=st.TEXT,
                                 stroke_width=2, dash_length=0.06)
                bars.add(VGroup(bar, VGroup(track, ref, lab_)))
        bar_t = st.T("appendage density at each site (dashed tick = a pure state)",
                     st.BODY_SIZE - 4, st.MUTED)
        bar_t.move_to([0, maps[0].get_bottom()[1] - 1.5, 0])
        foot = st.caption("maps idealised: w·ρA + (1-w)·ρB, w = state-A share\n"
                          "(consensus: all particles; classes: soft assignments)",
                          size=st.SMALL_SIZE)
        self.play(FadeIn(maps[0]), FadeIn(names[0]), FadeIn(pcts[0]), run_time=0.6)
        self.play(FadeIn(maps[1]), FadeIn(maps[2]), FadeIn(names[1:]), FadeIn(pcts[1:]),
                  FadeIn(site_rings), run_time=1.0)
        self.play(FadeIn(VGroup(*[b[1] for b in bars])),
                  LaggedStart(*[GrowFromEdge(b[0], LEFT) for b in bars], lag_ratio=0.08),
                  FadeIn(bar_t), FadeIn(foot), run_time=1.3)
        self.wait(2.5)                                                   # ~41 s

        # ------------------------------------------------ beat 5: Figure S6
        # shrink the three maps (labels stay at full size) to make room below
        H5 = 2.0
        anims = [FadeOut(VGroup(bars, bar_t, foot, pcts, site_rings))]
        for k, x in enumerate(xs):
            m = maps[k]
            grp = m
            ny = top_y - 0.3 - 0.15 - H5 / 2
            anims.append(grp.animate.scale(H5 / H3, about_point=m.get_center())
                         .shift(np.array([x, ny, 0]) - m.get_center()))
        self.play(*anims, run_time=0.9)
        s6_rows = [
            st.T("Day et al. (2025), Figure S6: real data, not this simulation",
                 st.BODY_SIZE - 2, st.BADGE, weight="BOLD"),
            st.T("•  3D classification with 10 initial classes", st.BODY_SIZE - 4),
            st.T("•  PCA initialisation (3 components); 6 Å filter resolution", st.BODY_SIZE - 4),
            st.T("•  Aβ: 2 classes with better B56γ density seeded a 2-class\n"
                 "    re-run on all particles", st.BODY_SIZE - 4),
            st.T("Not two validated conformations, not k-means on a 2D plot;\n"
                 "the 2-class demo above is deliberately simplified.", st.BODY_SIZE - 4,
                 st.MUTED, slant="ITALIC"),
        ]
        s6 = VGroup(*s6_rows).arrange(DOWN, aligned_edge=LEFT, buff=0.18)
        if s6.width > 12.4:
            s6.scale_to_fit_width(12.4)
        s6_box = RoundedRectangle(corner_radius=0.15, width=s6.width + 0.6, height=s6.height + 0.45,
                                  stroke_color=st.MUTED, stroke_width=1.5, fill_color=st.PANEL,
                                  fill_opacity=0.95)
        s6.move_to(s6_box)
        s6g = VGroup(s6_box, s6)
        s6g.move_to([0, -config.frame_height / 2 + st.MARGIN + s6g.height / 2, 0])
        self.play(FadeIn(s6g), run_time=0.8)
        self.wait(5.0)                                                   # ~45 s

        # ------------------------------------------------ takeaway freeze
        card = st.takeaway_card(
            r"p(I_i)=\sum\nolimits_{c}\pi_c\,p(I_i\mid\rho_c)",
            "Mixed states blur a consensus map; classes separate them.",
            "related question to multiconformer crystallography, but\n"
            "3D classes are inferred from single-particle projections.")
        card.move_to([0, -config.frame_height / 2 + st.MARGIN + card.height / 2, 0])
        self.play(FadeOut(s6g), FadeIn(card), run_time=0.9)
        self.wait(4.5)
