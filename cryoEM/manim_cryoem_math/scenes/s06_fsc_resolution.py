"""Clip 06 — FSC, global/local resolution and map-model comparisons.

Storyboard beats (spec section 3, Scene 06) -> code sections:
    B1  independent halves      observation cloud -> two disjoint halves -> two
                                independently estimated dot-position sets
    B2  Fourier shells and FSC  |F_A|, |F_B| on the k_z = 0 plane with a growing
                                shell |k| = s, FSC(s) traced from the 3D fields
    B3  resolution readout      FSC = 0.143 crossing, d = 1/s readout as s grows
    B4  global versus local     Region 1 (dots 1, 3, 4) vs Region 2 (dots 2, 5),
                                region-wise FSC curves
    B5  map-model               point model -> kernel density vs reconstructed
                                map; separate map-model curve, illustrative 0.5
    B6  summary + takeaway      three columns, then a static takeaway card

Every number on screen comes from cached("fsc_resolution", ne.fsc_resolution):
the FSC curves were computed on 64^3 numerical fields with independent noise in
the two halves (soft masks around the dots; empty shells return 0).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, noise_and_ctf as nc
from helpers.cache import cached

HALF_A = st.OBS        # cyan: half-map A (spec storyboard)
HALF_B = st.EST        # orange: half-map B
MODEL = st.TRUTH       # green: the ideal point model (the true dots)
CURVE = st.TEXT        # half-map and region-wise FSC curves


# ------------------------------------------------------------ local numbers
def layout_numbers(n_obs: int = 48, crop: int = 18) -> dict:
    """Positions of the schematic observation markers, their split into two
    disjoint halves, and the displayed map projections (cropped to the dots)."""
    g = nc.rng(600)
    ang = g.uniform(0, 2 * np.pi, n_obs)
    rad = np.sqrt(g.uniform(0.05, 1.0, n_obs))
    cloud = np.stack([2.3 * rad * np.cos(ang), 1.25 * rad * np.sin(ang)], 1)
    order = g.permutation(n_obs)
    d = cached("fsc_resolution", ne.fsc_resolution)
    c = d["proj_A"].shape[0] // 2
    sl = slice(c - crop, c + crop)
    full = 0.5 * (d["proj_A"] + d["proj_B"])
    return dict(cloud=cloud, half_of=(order >= n_obs // 2).astype(int),
                proj_map=full[sl, sl], proj_model=d["proj_model"][sl, sl])


def hull(points: np.ndarray) -> np.ndarray:
    """Convex hull (monotone order), counter-clockwise."""
    p = sorted(map(tuple, points))

    def turn(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    def half(seq):
        out = []
        for q in seq:
            while len(out) >= 2 and turn(out[-2], out[-1], q) <= 0:
                out.pop()
            out.append(q)
        return out
    lo, up = half(p), half(p[::-1])
    return np.array(lo[:-1] + up[:-1])


class Scene06FSCResolution(Scene):
    def construct(self):
        d = cached("fsc_resolution", ne.fsc_resolution)
        L = cached("s06_fsc_resolution_layout", layout_numbers)
        self.d = d
        s, nyq = d["s"], float(d["nyquist"])
        X = d["X"]

        self.tag = st.illustration_label()
        self.add(self.tag)
        tc = st.title_card("Three resolution measures",
                           "Which quantities are correlated, and why should\n"
                           "three resolution measurements not be conflated?")
        self.title, self.question = tc[0], tc[1]
        self.heading = None
        self.play(FadeIn(tc), run_time=1.0)

        # ============================================================ B1 halves
        n_obs = len(L["cloud"])
        cloud_c = np.array([0.0, -0.6, 0])
        # observations are images: drawn as tiny squares, not as dots (dots = source)
        obs = VGroup(*[Square(side_length=0.13, color=st.MUTED, fill_opacity=1,
                              stroke_width=0).move_to(cloud_c + [x, y, 0])
                       for x, y in L["cloud"]])
        obs_lab = st.T(f"{n_obs} observations (images)", 30, st.MUTED)
        obs_lab.next_to(obs, DOWN, buff=0.3)
        self.play(LaggedStart(*[FadeIn(o, scale=0.5) for o in obs], lag_ratio=0.02),
                  FadeIn(obs_lab), run_time=1.4)

        # split into two disjoint halves
        n_half = d["est_A"].shape[0]
        idxA = [i for i in range(n_obs) if L["half_of"][i] == 0][:n_half]
        idxB = [i for i in range(n_obs) if L["half_of"][i] == 1][:n_half]
        cA, cB = np.array([-4.0, -0.55, 0]), np.array([4.0, -0.55, 0])
        anims = [obs[i].animate.move_to(cA + 0.55 * (obs[i].get_center() - cloud_c))
                 .set_color(HALF_A) for i in idxA]
        anims += [obs[i].animate.move_to(cB + 0.55 * (obs[i].get_center() - cloud_c))
                  .set_color(HALF_B) for i in idxB]
        split_lab = st.T("two disjoint halves", 30, st.TEXT).move_to(cloud_c + 0.25 * UP)
        split_sub = st.T(f"{n_half} + {n_half}, independent noise", 26, st.MUTED)
        split_sub.next_to(split_lab, DOWN, buff=0.15)
        self.play(*anims, FadeOut(obs_lab), run_time=1.6)
        self.play(FadeIn(split_lab), FadeIn(split_sub), run_time=0.7)

        # each observation of a half -> one independently estimated set of 5 positions
        scatA = self.estimate_scatter(d["est_A"], cA + 0.15 * UP, 1.05, HALF_A)
        scatB = self.estimate_scatter(d["est_B"], cB + 0.15 * UP, 1.05, HALF_B)
        morph = []
        for idx, scat in ((idxA, scatA), (idxB, scatB)):
            for k, i in enumerate(idx):
                src = VGroup(*[obs[i].copy().scale(0.7) for _ in range(X.shape[0])])
                morph.append(ReplacementTransform(src, scat[k]))
        labA = st.T("Independent half-map A", 30, HALF_A).move_to([cA[0], -2.75, 0])
        labB = st.T("Independent half-map B", 30, HALF_B).move_to([cB[0], -2.75, 0])
        est_note = st.T("each half: its own estimated dot positions", 26, st.MUTED)
        est_note.to_edge(DOWN, buff=st.MARGIN)
        self.remove(*obs)
        self.play(*morph, run_time=2.0)
        self.play(FadeIn(labA), FadeIn(labB), FadeIn(est_note), run_time=0.8)
        self.wait(1.0)

        # ============================================================ B2 shells + FSC
        self.set_heading("Fourier shells and the FSC curve",
                         FadeOut(VGroup(scatA, scatB, labA, labB, est_note, split_lab,
                                        split_sub)), first=True)

        hm_h = 2.3
        ftA = st.heatmap(d["ft2_A"], height=hm_h).move_to([-5.15, 0.6, 0])
        ftB = st.heatmap(d["ft2_B"], height=hm_h).move_to([-2.45, 0.6, 0])
        frA, frB = st.framed(ftA, HALF_A), st.framed(ftB, HALF_B)
        fA_lab = st.M(r"|F_A|", size=36, color=HALF_A).next_to(ftA, UP, buff=0.12)
        fB_lab = st.M(r"|F_B|", size=36, color=HALF_B).next_to(ftB, UP, buff=0.12)
        ft_note = st.M(r"\log|F|\ \text{on the plane}\ k_z=0", size=36, color=st.MUTED)
        ft_note.next_to(Group(ftA, ftB), DOWN, buff=0.18)

        ax = self.fsc_axes()
        self.ax = ax
        s_t = ValueTracker(float(s[0]))
        npx = d["ft2_A"].shape[0]

        def ring(img):
            # zero frequency sits at pixel n/2 (half a pixel right of / above the centre)
            o = img.get_center() + 0.5 * img.height / npx * (RIGHT + UP)
            r = s_t.get_value() / nyq * img.height / 2
            c = Circle(radius=max(r, 0.01), color=st.TEXT, stroke_width=3).move_to(o)
            c.set_stroke(st.BG, width=7, background=True)
            return c
        ringA = always_redraw(lambda: ring(ftA))
        ringB = always_redraw(lambda: ring(ftB))
        shell_lab = st.M(r"\text{shell } S_s:\ |\mathbf k| = s", size=36)
        shell_lab.next_to(ft_note, DOWN, buff=0.25)

        half_curve = always_redraw(lambda: self.curve(d["fsc_half"], CURVE,
                                                      smax=s_t.get_value()))
        head = always_redraw(lambda: Dot(ax.c2p(s_t.get_value(),
                                                self.interp(d["fsc_half"], s_t.get_value())),
                                         radius=0.08, color=CURVE))
        s_read = VGroup(st.M(r"s =", size=36), DecimalNumber(s[0], 2, font_size=36),
                        st.T("1/units", 28, st.MUTED)).arrange(RIGHT, buff=0.15)
        s_read.move_to(ax.c2p(4.15, 0.85))
        s_read[1].add_updater(lambda m: m.set_value(s_t.get_value()))

        self.play(FadeIn(Group(ftA, frA, ftB, frB)), FadeIn(fA_lab), FadeIn(fB_lab),
                  FadeIn(ft_note), Create(ax), FadeIn(self.ax_labels), run_time=1.2)
        self.add(ringA, ringB, half_curve, head)
        self.play(FadeIn(shell_lab), FadeIn(s_read), run_time=0.6)
        self.play(s_t.animate.set_value(float(s[-1])), run_time=4.0, rate_func=linear)
        for m in (half_curve, head, s_read[1], ringA, ringB):
            m.clear_updaters()
        curve_lab = st.T("half-map\nA vs B", 26, CURVE)
        curve_lab.move_to(ax.c2p(0.15, 0.40), aligned_edge=LEFT)

        # the FSC definition replaces the shell label
        eq = st.M(r"\mathrm{FSC}(s)=\frac{\sum_{\mathbf k\in S_s} F_A(\mathbf k)\,"
                  r"F_B^*(\mathbf k)}{\sqrt{\sum_{S_s}|F_A|^2\,\sum_{S_s}|F_B|^2}}", size=36)
        for i in (12, 13, 30, 31):
            eq[0][i].set_color(HALF_A)
        for i in (17, 18, 19, 38, 39):
            eq[0][i].set_color(HALF_B)
        eq.move_to([-3.8, -2.4, 0])
        mask_note = st.T("FSC of 3D fields, same soft mask around the dots for A and B",
                         26, st.MUTED).to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeOut(VGroup(ringA, ringB, s_read, head)),
                  Succession(FadeOut(shell_lab, shift=0.1 * UP), FadeIn(eq, shift=0.1 * UP)),
                  FadeIn(curve_lab), FadeIn(mask_note), run_time=1.2)
        self.wait(1.2)

        # ============================================================ B3 readout
        thr = self.hline(0.143, CURVE)
        thr_lab = st.M(r"\mathrm{FSC}=0.143", size=36)
        thr_lab.next_to(ax.c2p(5.35, 0.143), UP, buff=0.1).align_to(ax.c2p(5.35, 0), RIGHT)
        s_half, d_half = float(d["s_half"]), float(d["d_half"])
        cross = Dot(ax.c2p(s_half, 0.143), radius=0.1, color=st.TEXT)
        cross_ring = Circle(radius=0.2, color=st.TEXT, stroke_width=3).move_to(cross)
        drop = DashedLine(ax.c2p(s_half, 0.143), ax.c2p(s_half, -0.2), color=st.MUTED,
                          stroke_width=2.5, dash_length=0.08)

        self.set_heading("Resolution readout",
                         FadeOut(Group(ftA, frA, ftB, frB, fA_lab, fB_lab, ft_note, eq,
                                       mask_note)))
        self.play(Create(thr), FadeIn(thr_lab), run_time=0.7)
        self.play(FadeIn(cross, scale=0.5), Create(cross_ring), Create(drop), run_time=0.6)

        d_eq = st.M(r"d = \frac{1}{s}", size=48).move_to([-3.8, 1.2, 0])
        r_t = ValueTracker(1.0)
        s_val = DecimalNumber(1.0, 2, font_size=40)
        d_val = DecimalNumber(1.0, 2, font_size=40)
        row_s = VGroup(st.M(r"s =", size=40), s_val, st.T("1/units", 30, st.MUTED))
        row_d = VGroup(st.M(r"d =", size=40), d_val, st.T("units", 30, st.MUTED))
        for r in (row_s, row_d):
            r.arrange(RIGHT, buff=0.18)
        rows = VGroup(row_s, row_d).arrange(DOWN, aligned_edge=LEFT, buff=0.3)
        rows.move_to([-3.8, -0.35, 0])
        s_val.add_updater(lambda m: m.set_value(r_t.get_value()))
        d_val.add_updater(lambda m: m.set_value(1.0 / r_t.get_value()))
        bar_scale = 3.0     # screen units per unit length
        bar_left = np.array([-5.3, -1.85, 0])

        def length_bar():
            e = bar_left + bar_scale / r_t.get_value() * RIGHT
            return VGroup(Line(bar_left, e, color=st.TEXT, stroke_width=6),
                          Line(bar_left + 0.12 * UP, bar_left + 0.12 * DOWN, color=st.TEXT,
                               stroke_width=3),
                          Line(e + 0.12 * UP, e + 0.12 * DOWN, color=st.TEXT, stroke_width=3))
        bar = always_redraw(length_bar)
        bar_lab = st.T("length d", 26, st.MUTED)
        bar_lab.next_to(bar_left, DOWN, buff=0.3, aligned_edge=LEFT)
        marker = always_redraw(lambda: Triangle(color=st.TEXT, fill_opacity=1, stroke_width=0)
                               .scale(0.11).rotate(PI).next_to(ax.c2p(r_t.get_value(), -0.2),
                                                               UP, buff=0.0))
        self.play(FadeIn(d_eq), FadeIn(rows), FadeIn(bar_lab), run_time=0.7)
        self.add(bar, marker)
        self.play(r_t.animate.set_value(s_half), run_time=3.0, rate_func=smooth)
        for m in (s_val, d_val, bar, marker):
            m.clear_updaters()
        crit = st.T("0.143: conventional criterion for two independent half-maps", 28,
                    st.MUTED).to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(crit), run_time=0.6)
        self.wait(1.0)

        # ============================================================ B4 global vs local
        self.set_heading("Global versus local",
                         FadeOut(VGroup(d_eq, rows, bar, bar_lab, crit, marker, cross,
                                        cross_ring, drop, curve_lab)),
                         half_curve.animate.set_stroke(st.GRID, width=3))

        cR = np.array([-3.7, 0.75, 0])
        sc4 = 1.0
        # interleave the two halves and keep the dots translucent so that the tight
        # Region 1 clusters show both colours (B drawn last would hide A completely)
        sA = self.estimate_scatter(d["est_A"], cR, sc4, HALF_A, radius=0.04)
        sB = self.estimate_scatter(d["est_B"], cR, sc4, HALF_B, radius=0.04)
        pool = VGroup(*[g for k in range(len(sA)) for g in (sA[k], sB[k])])
        pool.set_fill(opacity=0.6)
        reg1 = [int(j) for j in d["region1"]]
        reg2 = [int(j) for j in d["region2"]]
        box1 = self.region_outline(X[reg1, :2], cR, sc4, pad=0.18, dashed=False)
        box2 = self.region_outline(X[reg2, :2], cR, sc4, pad=0.42, dashed=True)
        r1_lab = st.T("Region 1", 28, st.TEXT).next_to(box1, UP, buff=0.1)
        r1_lab.align_to(box1, LEFT)
        r2_lab = st.T("Region 2", 28, st.TEXT).next_to(box2, RIGHT, buff=0.15)
        for lab in (r1_lab, r2_lab):
            lab.add_background_rectangle(color=st.BG, opacity=1, buff=0.05)
            lab.set_z_index(2)

        def d_row(sample, sub, val, dots):
            nums = "all dots" if len(dots) == X.shape[0] else \
                "dots " + ", ".join(str(j + 1) for j in dots)
            return VGroup(sample, st.M(rf"d_{{{sub}}} =", size=36),
                          DecimalNumber(val, 2, font_size=36), st.T("units", 28, st.MUTED),
                          st.T(f"({nums})", 28, st.MUTED)).arrange(RIGHT, buff=0.15)
        # legend row for the greyed global curve (masked over all five dots)
        row0 = d_row(Line(ORIGIN, 0.6 * RIGHT, color=st.GRID, stroke_width=4),
                     r"\mathrm{all}", float(d["d_half"]), list(range(X.shape[0])))
        row1 = d_row(Line(ORIGIN, 0.6 * RIGHT, color=CURVE, stroke_width=4), 1,
                     float(d["d_region1"]), reg1)
        row2 = d_row(DashedLine(ORIGIN, 0.6 * RIGHT, color=CURVE, stroke_width=4,
                                dash_length=0.1), 2, float(d["d_region2"]), reg2)
        drows = VGroup(row0, row1, row2).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
        drows.move_to([cR[0], -1.45, 0])
        concept = st.T("Estimate spread: a conceptual visual,\n"
                       "not a per-dot error bar.", 26, st.MUTED)
        concept.next_to(drows, DOWN, buff=0.18)

        self.play(FadeIn(pool), run_time=0.9)
        self.play(Create(box1), Create(box2), FadeIn(r1_lab), FadeIn(r2_lab), run_time=0.9)

        c1 = self.curve(d["fsc_region1"], CURVE)
        c2 = self.curve(d["fsc_region2"], CURVE, dashed=True)
        x1 = Dot(ax.c2p(float(d["s_region1"]), 0.143), radius=0.09, color=CURVE)
        x2 = Dot(ax.c2p(float(d["s_region2"]), 0.143), radius=0.09, color=CURVE)
        # the Region 2 curve rebounds above 0.143 later; d_2 uses the first crossing
        x2_lab = st.T("first\ncrossing", 26, st.MUTED).next_to(x2, DL, buff=0.06)
        x2_lab.add_background_rectangle(color=st.BG, opacity=1, buff=0.05)
        x2_lab.set_z_index(2)
        self.play(Create(c1), Create(c2), run_time=1.5)
        self.play(FadeIn(x1, scale=0.5), FadeIn(x2, scale=0.5), FadeIn(x2_lab),
                  FadeIn(drows), run_time=0.7)
        cap4 = st.T("A single global cutoff does not describe every location.", 30, st.TEXT)
        cap4.to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(cap4), FadeIn(concept), run_time=0.6)
        self.wait(1.4)

        # ============================================================ B5 map-model
        thr_lab5 = st.M(r"0.143", size=36, color=st.MUTED)
        thr_lab5.next_to(ax.c2p(5.35, 0.143), UP, buff=0.1).align_to(ax.c2p(5.35, 0), RIGHT)
        # the grey global curve keeps a label once the d_all legend row fades out
        gh_lab = st.T("global\nhalf-map", 26, st.MUTED)
        gh_lab.move_to(ax.c2p(0.15, 0.36), aligned_edge=LEFT)
        self.set_heading("Map–model agreement",
                         FadeOut(VGroup(pool, box1, box2, r1_lab, r2_lab, concept, c1, c2,
                                        x1, x2, x2_lab, drows, cap4)),
                         thr.animate.set_stroke(st.GRID),
                         ReplacementTransform(thr_lab, thr_lab5), FadeIn(gh_lab))

        hm2 = 1.7
        pt_c = np.array([-5.75, 0.95, 0])
        model_img = st.heatmap(L["proj_model"], height=hm2).move_to([-1.95, 0.95, 0])
        map_img = st.heatmap(L["proj_map"], height=hm2).move_to([-1.95, -1.6, 0])
        model_fr = st.framed(model_img, MODEL)
        map_fr = st.framed(map_img, st.MUTED)
        # point model: the true dots drawn on the same coordinates as the heatmap pixels
        npx2 = L["proj_model"].shape[0]
        h = 2 * 3.0 / d["proj_A"].shape[0]          # field spacing (extent 3, 64 samples)
        px = hm2 / npx2
        pts = VGroup(*[Dot(pt_c + (np.array([x, y, 0]) / h + 0.5 * (RIGHT + UP)) * px,
                           radius=0.08, color=MODEL) for x, y in X[:, :2]])
        pt_frame = Square(side_length=hm2, color=st.GRID, stroke_width=2).move_to(pt_c)
        pt_lab = st.T("point model", 28, MODEL).next_to(pt_frame, UP, buff=0.12)
        pt_lab.align_to(pt_frame, LEFT)
        md_lab = st.T("model density", 28, MODEL).next_to(model_fr, UP, buff=0.12)
        kern = Arrow(pt_frame.get_right(), model_fr.get_left(), buff=0.1, color=st.MUTED,
                     stroke_width=4, max_tip_length_to_length_ratio=0.25)
        kern_lab = st.T("Gaussian\nkernels", 26, st.MUTED).next_to(kern, DOWN, buff=0.15)
        map_lab = VGroup(st.T("reconstructed map", 28, st.TEXT),
                         st.M(r"(A+B)/2", size=36)).arrange(DOWN, buff=0.12)
        map_lab.next_to(map_fr, LEFT, buff=0.3).align_to(map_fr, UP)
        proj_note = st.T("projections shown;\nFSC uses masked\n3D fields", 26, st.MUTED)
        proj_note.next_to(map_lab, DOWN, buff=0.3).align_to(map_lab, RIGHT)
        cmp = DoubleArrow(model_fr.get_bottom(), map_fr.get_top(), buff=0.06,
                          color=MODEL, stroke_width=4, tip_length=0.18)
        cmp_lab = st.T("compare", 26, MODEL).next_to(cmp, RIGHT, buff=0.12)

        self.play(FadeIn(pts), Create(pt_frame), FadeIn(pt_lab), run_time=0.8)
        self.play(GrowArrow(kern), FadeIn(kern_lab), FadeIn(model_img), Create(model_fr),
                  FadeIn(md_lab), run_time=0.9)
        self.play(FadeIn(map_img), Create(map_fr), FadeIn(map_lab), GrowArrow(cmp),
                  FadeIn(cmp_lab), FadeIn(proj_note), run_time=0.9)

        mm_curve = self.curve(d["fsc_model"], MODEL)
        mm_lab = st.T("map vs model", 26, MODEL)
        mm_lab.move_to(ax.c2p(0.1, 0.74), aligned_edge=LEFT)
        thr5 = self.hline(0.5, MODEL)
        thr5_lab = VGroup(st.T("illustrative", 28, MODEL),
                          st.M(r"\mathrm{FSC}=0.5", size=36, color=MODEL)
                          ).arrange(DOWN, aligned_edge=RIGHT, buff=0.1)
        thr5_lab.next_to(ax.c2p(5.35, 0.5), UP, buff=0.12).align_to(ax.c2p(5.35, 0), RIGHT)
        s_model, d_model = float(d["s_model"]), float(d["d_model"])
        xm = Dot(ax.c2p(s_model, 0.5), radius=0.1, color=MODEL)
        dm = VGroup(st.M(r"d_{\mathrm{model}} =", size=44, color=MODEL),
                    DecimalNumber(d_model, 2, font_size=44, color=MODEL),
                    st.T("units", 28, st.MUTED)).arrange(RIGHT, buff=0.15)
        dm.next_to(ax.c2p(5.35, 1.0), UP, buff=0.08).align_to(ax.c2p(5.35, 0), RIGHT)
        self.play(Create(mm_curve), FadeIn(mm_lab), run_time=1.3)
        self.play(Create(thr5), FadeIn(thr5_lab), FadeIn(xm, scale=0.5), FadeIn(dm),
                  run_time=0.8)
        cap5 = st.T("A different statistic and cutoff convention from the half-map 0.143.",
                    28, st.TEXT).to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(cap5), run_time=0.6)
        self.wait(1.4)

        # ============================================================ B6 summary
        left_over = [m for m in self.mobjects if m is not self.title and m is not self.tag]
        self.play(*[FadeOut(m) for m in left_over], run_time=0.6)
        self.heading = None

        cols = VGroup()
        specs = [("Global", "half-map vs half-map", r"\mathrm{FSC}(A,B)\ \text{at}\ 0.143"),
                 ("Local", "spatially varying\ninformation", r"\text{region-wise FSC}"),
                 ("Model", "calculated model\ndensity vs map",
                  r"\mathrm{FSC}(\text{map},\,\text{model})")]
        for (hd, body, mth), col in zip(specs, [CURVE, CURVE, MODEL]):
            cols.add(VGroup(st.T(hd, 32, col, weight="BOLD"), st.T(body, 28, st.TEXT),
                            st.M(mth, size=36, color=st.MUTED)).arrange(DOWN, buff=0.16))
        cols.arrange(RIGHT, buff=1.0, aligned_edge=UP)
        card = st.takeaway_card(
            ["Half-map FSC compares two independent halves (0.143).",
             "Local resolution: reliable detail varies in space.",
             "Map–model FSC: model density vs map (illustrative 0.5).",
             "Neither cutoff is a positional error bar or physical limit."],
            math=r"d = 1/s")
        card.to_edge(DOWN, buff=st.MARGIN)
        cols.next_to(card, UP, buff=0.3)
        seps = VGroup(*[Line(UP, DOWN, color=st.GRID, stroke_width=2).scale(0.9).move_to(
            [(cols[i].get_right()[0] + cols[i + 1].get_left()[0]) / 2,
             cols.get_center()[1], 0]) for i in range(2)])
        self.play(LaggedStart(*[FadeIn(c, shift=0.1 * UP) for c in cols], lag_ratio=0.3),
                  FadeIn(seps), run_time=1.5)

        self.play(FadeIn(card), run_time=0.8)
        self.wait(4.0)

    # ------------------------------------------------------------ helpers
    def set_heading(self, text, *extra, first=False):
        """Swap the section heading; `extra` animations (fading the previous section's
        content) run in the same play so a heading never sits over old content."""
        new = st.T(text, 30, st.MUTED)
        new.next_to(self.title, DOWN, aligned_edge=LEFT, buff=0.18)
        old = self.question if first else self.heading
        self.play(Succession(FadeOut(old, shift=0.1 * UP), FadeIn(new, shift=0.1 * UP)),
                  *extra, run_time=0.8)
        self.heading = new

    def estimate_scatter(self, est, center, scale, color, radius=0.045):
        """One VGroup per estimate copy (its five dots), xy coordinates."""
        return VGroup(*[VGroup(*[Dot(center + scale * np.array([p[0], p[1], 0]),
                                     radius=radius, color=color, fill_opacity=0.85)
                                 for p in copy]) for copy in est])

    def region_outline(self, pts, center, scale, pad, dashed):
        """Straight-edged neighbourhood indicator: hull of the dots padded by an octagon."""
        ang = np.linspace(0, 2 * np.pi, 8, endpoint=False) + np.pi / 8
        ring = np.stack([np.cos(ang), np.sin(ang)], 1) * pad / np.cos(np.pi / 8)
        hv = hull((pts[:, None, :] + ring[None]).reshape(-1, 2))
        poly = Polygon(*[center + scale * np.array([x, y, 0]) for x, y in hv],
                       color=st.TEXT, stroke_width=2.5)
        return DashedVMobject(poly, num_dashes=40) if dashed else poly

    def fsc_axes(self):
        # Axes provides the coordinates; the s axis is drawn along the bottom (FSC = -0.2)
        # so that negative FSC values never run into the tick labels.
        ax = Axes(x_range=[0, 5.4, 1], y_range=[-0.2, 1.0, 0.2], x_length=5.5, y_length=3.9,
                  tips=False, axis_config=dict(color=st.GRID, stroke_width=2))
        ax.move_to([3.75, -0.15, 0])
        ax.x_axis.set_opacity(0)
        bottom = VGroup(Line(ax.c2p(0, -0.2), ax.c2p(5.4, -0.2), color=st.GRID,
                             stroke_width=2))
        for k in range(1, 6):
            bottom.add(Line(ax.c2p(k, -0.2), ax.c2p(k, -0.2) + 0.12 * UP, color=st.GRID,
                            stroke_width=2))
            bottom.add(st.T(str(k), 28, st.MUTED).next_to(ax.c2p(k, -0.2), DOWN, buff=0.12))
        for v, txt in ((0.0, "0"), (0.5, "0.5"), (1.0, "1")):
            bottom.add(st.T(txt, 28, st.MUTED).next_to(ax.c2p(0, v), LEFT, buff=0.15))
        xl = VGroup(st.M(r"s", size=36, color=st.MUTED), st.T("(1/units)", 28, st.MUTED))
        xl.arrange(RIGHT, buff=0.12).next_to(ax.c2p(2.7, -0.2), DOWN, buff=0.5)
        yl = st.M(r"\mathrm{FSC}", size=36, color=st.MUTED)
        yl.next_to(ax.c2p(0, 1.0), UP, buff=0.42)
        zero = DashedLine(ax.c2p(0, 0), ax.c2p(5.4, 0), color=st.GRID, stroke_width=1.5,
                          dash_length=0.06)
        self.ax_labels = VGroup(bottom, xl, yl, zero)
        return ax

    def interp(self, f, x):
        return float(np.interp(x, self.d["s"], f))

    def curve(self, f, color, smax=None, dashed=False):
        """Computed FSC values joined by straight segments (no smoothing)."""
        s = self.d["s"]
        if smax is None:
            smax = s[-1]
        keep = s <= smax
        xs = list(s[keep]) + ([smax] if smax > s[0] else [])
        ys = [self.interp(f, x) for x in xs]
        if len(xs) < 2:
            xs, ys = [s[0], s[0] + 1e-3], [f[0], f[0]]
        v = VMobject(color=color, stroke_width=4)
        v.set_points_as_corners([self.ax.c2p(x, y) for x, y in zip(xs, ys)])
        if dashed:
            return DashedVMobject(v, num_dashes=60, dashed_ratio=0.55)
        return v

    def hline(self, y, color):
        return DashedLine(self.ax.c2p(0, y), self.ax.c2p(5.4, y), color=color,
                          stroke_width=2.5, dash_length=0.12)
