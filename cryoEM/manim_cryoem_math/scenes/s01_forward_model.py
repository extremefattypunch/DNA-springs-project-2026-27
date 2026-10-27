"""Clip 01 — Projection, optical transfer and noise. Storyboard beats -> code sections (spec §3).

Beat 1 (projection)  : oblique 2D drawing of 3D axes, five dots, a viewing arrow and a
                       flat detector plane; vertical guide lines carry each dot to the plane.
Beat 2 (CTF)         : 1D illustration p(x) (y line integral of the projected image) with a
                       2D inset; oscillating CTF(k) with a zero and a negative lobe; three
                       frequency bands before/after the CTF (reversed, removed, kept).
Beat 3 (noise)       : s(x) + eps(x) = I(x) on three rows with one vertical scale; the 2D
                       inset gains noise from a precomputed series.
Beat 4 (statement)   : the Fourier-space observation equation, terms lit in turn, then the
                       inverse-problem takeaway.

Every number on screen is read from cached("forward_model", ne.forward_model).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne
from helpers.cache import cached

# ------------------------------------------------------------------ oblique drawing
# screen = O + SC * (x * EX + y * EY + z * EZ); z (the viewing axis) is drawn vertical,
# so projection lines along z are vertical guide lines on screen.
EX = np.array([1.0, 0.0, 0.0])
EY = np.array([0.52, 0.30, 0.0])
EZ = np.array([0.0, 1.0, 0.0])
SC = 0.92
O3 = np.array([-2.35, -0.05, 0.0])
Z_DET = -1.85          # detector plane height (detector-frame z)
HALF = 1.75            # detector half-width (units)


def obl(x, y, z):
    return O3 + SC * (x * EX + y * EY + z * EZ)


def curve(ax, x, y, color, width=4):
    """Polyline of sampled data on Axes ax (no smoothing: these are samples)."""
    pts = [ax.c2p(a, b) for a, b in zip(x, y)]
    return VMobject(stroke_color=color, stroke_width=width).set_points_as_corners(pts)


def clear_of_tag(mob, tag_left, gap=0.3):
    """Work-around: st.TITLE_MAX_W (9.6) lets a title run under the tag (left edge ~3.0)."""
    max_w = tag_left - gap - mob.get_left()[0]
    if mob.width > max_w:
        mob.scale(max_w / mob.width, about_edge=LEFT)
    return mob


def fmt(v, nd=2):
    s = f"{v:+.{nd}f}" if abs(v) >= 0.5 * 10 ** -nd else f"{0:.{nd}f}"
    return s.replace("-", "−")


class Scene01ForwardModel(Scene):
    def construct(self):
        d = cached("forward_model", ne.forward_model)
        tag = st.illustration_label()
        self.add(tag)
        self.tag_left = tag.get_left()[0]
        self.beat1_projection(d)
        self.beat2_ctf(d)
        self.beat3_noise(d)
        self.beat4_statement(d)

    # ================================================================ beat 1: 0–6.5 s
    def beat1_projection(self, d):
        title = st.title_card("Projection, optical transfer and noise",
                              "How do projection, the CTF and noise\nalter an observation?")
        clear_of_tag(title[0], self.tag_left)
        self.play(FadeIn(title, shift=0.1 * DOWN), run_time=0.8)

        Q = d["points3d"] @ d["R"].T            # dots in the detector frame (R x_j)
        P2 = d["proj2d"]                        # (R x_j)_xy
        assert np.allclose(Q[:, :2], P2)

        axes3 = VGroup()
        for end, lab in ((obl(1.55, 0, 0), "x"), (obl(0, 1.55, 0), "y"), (obl(0, 0, 1.45), "z")):
            a = Arrow(obl(0, 0, 0), end, buff=0, stroke_width=2.5, color=st.MUTED,
                      tip_length=0.16, max_tip_length_to_length_ratio=0.15)
            t = st.M(lab, size=36, color=st.MUTED).next_to(end, UR, buff=0.04)
            axes3.add(a, t)

        dots3 = VGroup(*[Dot(obl(*q), radius=0.1, color=st.TRUTH) for q in Q])
        corners = [obl(-HALF, -HALF, Z_DET), obl(HALF, -HALF, Z_DET),
                   obl(HALF, HALF, Z_DET), obl(-HALF, HALF, Z_DET)]
        plane = Polygon(*corners, stroke_color=st.MUTED, stroke_width=2,
                        fill_color=st.GRID, fill_opacity=0.45)
        plane_lab = st.T("detector plane", st.SMALL_SIZE + 2, st.MUTED)
        plane_lab.next_to(plane, DOWN, buff=0.12).align_to(plane, LEFT).shift(0.1 * RIGHT)

        vx, vy = -1.55, -1.25
        view = Arrow(obl(vx, vy, 2.0), obl(vx, vy, 0.55), buff=0, color=st.TEXT,
                     stroke_width=6, tip_length=0.25, max_tip_length_to_length_ratio=0.25)
        # no sign implied: the helper convention projects along +z (v = R^T e_z), and the
        # line integral is the same either way
        view_lab = st.T("projection\nalong z", st.SMALL_SIZE + 2, st.TEXT)
        view_lab.next_to(view, LEFT, buff=0.15)

        guides = VGroup(*[DashedLine(obl(*q), obl(q[0], q[1], Z_DET), dash_length=0.08,
                                     stroke_width=2, color=st.MUTED, stroke_opacity=0.7)
                          for q in Q])
        pdots = VGroup(*[Dot(obl(p[0], p[1], Z_DET), radius=0.09, color=st.OBS) for p in P2])

        eq = st.M(r"\mathbf p_j=(R\,\mathbf x_j)_{xy}", size=44)
        eq_note = st.T("z is integrated out:\none view cannot fix depth", st.SMALL_SIZE + 2,
                       st.MUTED)
        leg = st.legend([(st.TRUTH, "rotated source dots (3D)"), (st.OBS, "projected dots (2D)")],
                        size=st.SMALL_SIZE + 2)
        side = VGroup(eq, eq_note, leg).arrange(DOWN, aligned_edge=LEFT, buff=0.35)
        side.move_to([4.25, -0.35, 0])
        cap = st.caption("3D coordinates → 2D projection", color=st.TEXT, size=st.BODY_SIZE)

        self.play(FadeIn(axes3), LaggedStart(*[GrowFromCenter(p) for p in dots3], lag_ratio=0.15),
                  run_time=1.0)
        self.play(FadeIn(plane), FadeIn(plane_lab), GrowArrow(view), FadeIn(view_lab),
                  run_time=0.8)
        self.play(LaggedStart(*[Create(g) for g in guides], lag_ratio=0.12), run_time=1.0)
        self.play(LaggedStart(*[FadeIn(p, scale=0.5) for p in pdots], lag_ratio=0.12),
                  FadeIn(cap), FadeIn(side), run_time=0.9)
        self.wait(1.6)
        self.play(FadeOut(VGroup(title, axes3, dots3, plane, plane_lab, view, view_lab, guides,
                                 pdots, side, cap)), run_time=0.5)


    # ================================================================ beat 2: 6–16 s
    def beat2_ctf(self, d):
        x, p, s = d["x"], d["s"], d["s_ctf"]
        neg_col = st.DIVERGING[0][1]          # the signed-heatmap colour for negative values
        head = VGroup(st.T("Optical transfer:", st.BODY_SIZE, st.TEXT),
                      st.M(r"\hat s(k)=\mathrm{CTF}(k)\,\hat p(k)", size=40)
                      ).arrange(RIGHT, buff=0.25)
        head.to_corner(UL, buff=st.MARGIN)
        # the CTF is a toy curve: say so next to its equation (spec 0.6), with the units
        ctf_note = st.T("illustrative CTF, no envelope; x in units, k in 1/units",
                        st.SMALL_SIZE, st.MUTED).next_to(head, DOWN, aligned_edge=LEFT, buff=0.1)

        # ---- left: 2D inset and the 1D profile
        inset = st.heatmap(d["img"], height=1.7)
        inset.move_to([-5.55, 1.6, 0])
        inset_fr = st.framed(inset)
        inset_lab = st.T("2D image", st.SMALL_SIZE, st.MUTED).next_to(inset, DOWN, buff=0.08)
        one_d = st.T("1D illustration", st.BODY_SIZE - 2, st.TEXT, slant="ITALIC")
        one_d_eq = st.M(r"p(x)=\textstyle\int \mathrm{image}\,dy", size=38)
        one_d_grp = VGroup(one_d, one_d_eq).arrange(DOWN, aligned_edge=LEFT, buff=0.18)
        one_d_grp.next_to(inset, RIGHT, buff=0.4).align_to(inset, UP).shift(0.15 * DOWN)

        ymax = 1.1 * float(np.max(p))
        ax1 = Axes(x_range=[-2, 2, 1], y_range=[-ymax, ymax, ymax], x_length=5.9,
                   y_length=2.45, tips=False,
                   axis_config=dict(color=st.GRID, stroke_width=2, include_ticks=False))
        ax1.move_to([-3.55, -1.35, 0])
        ax1_lab = st.M("x", size=36, color=st.MUTED).next_to(ax1.x_axis, RIGHT, buff=0.1)
        p_curve = curve(ax1, x, p, st.TRUTH, 5)
        p_lab = st.M("p(x)", size=40, color=st.TRUTH).move_to(ax1.c2p(-1.5, 0.62 * ymax))

        # s(x) is much smaller than p(x) (CTF(0) is near zero); magnify for display
        mag = round(float(np.max(np.abs(p)) / np.max(np.abs(s))), 1)
        s_curve = curve(ax1, x, s * mag, st.TEXT, 5)
        # the label states the display factor: the white curve is mag * s(x)
        s_lab = st.M(rf"{mag:.1f}\times s(x)", size=40, color=st.TEXT)
        mag_lab = st.T("(magnified for display)", st.SMALL_SIZE, st.MUTED)
        # stacked just above the plot box, under the 1D-illustration block: the magnified
        # curve fills the box, so labels inside it would sit on the curve
        VGroup(s_lab, mag_lab).arrange(DOWN, aligned_edge=LEFT, buff=0.06).next_to(
            ax1, UP, buff=0.04).align_to(one_d_grp, LEFT)

        inset_ctf = st.heatmap(d["img_ctf"], height=1.7).move_to(inset)
        # display-rescaled (own min/max): the true image is ~14x weaker, see the 1D factor
        inset_ctf_lab = st.T("after CTF\n(rescaled)", st.SMALL_SIZE, st.TEXT).next_to(
            inset, DOWN, buff=0.08)

        # ---- right: CTF(k) and the three bands
        k, c = d["k"], d["ctf"]
        kmax = 6.0
        sel = k <= kmax
        ax2 = Axes(x_range=[0, kmax, 1], y_range=[-1.15, 1.15, 1], x_length=5.4, y_length=2.1,
                   tips=False, axis_config=dict(color=st.GRID, stroke_width=2,
                                                include_ticks=False))
        ax2.move_to([3.55, 1.3, 0])
        ctf_curve = curve(ax2, k[sel], c[sel], st.TEXT, 4)
        k_lab = st.M("k", size=36, color=st.MUTED).next_to(ax2.x_axis, RIGHT, buff=0.1)
        ones = VGroup(st.M("+1", size=30, color=st.MUTED).next_to(ax2.c2p(0, 1), LEFT, 0.1),
                      st.M("-1", size=30, color=st.MUTED).next_to(ax2.c2p(0, -1), LEFT, 0.1))
        # inside the plot, upper left: the curve is <= 0 there (first negative lobe)
        assert mag_lab.get_right()[0] < ax2.get_left()[0] - 0.6, mag_lab.get_right()
        assert s_lab.get_left()[0] > inset.get_right()[0] + 0.2
        assert s_lab.get_top()[1] < one_d_grp.get_bottom()[1] - 0.2
        ctf_lab = st.M(r"\mathrm{CTF}(k)", size=36)
        ctf_lab.next_to(ax2.c2p(0.22, 1.0), DR, buff=0).shift(0.05 * DOWN)
        # the dashed k_1 marker passes behind the label
        ctf_lab.add_background_rectangle(color=st.BG, opacity=1, buff=0.04)
        ctf_lab.set_z_index(2)

        z0 = float(d["ctf_zeros"][0])
        lobe = k <= z0
        lobe_poly = Polygon(*([ax2.c2p(a, b) for a, b in zip(k[lobe], c[lobe])]
                              + [ax2.c2p(z0, 0), ax2.c2p(0, 0)]),
                            stroke_width=0, fill_color=neg_col, fill_opacity=0.5)
        lobe_lab = st.T("negative lobe:\ncontrast reversed", st.SMALL_SIZE, st.TEXT)
        lobe_lab.next_to(ax2.c2p(1.4, -1.0), DOWN, buff=0.12).align_to(ax2.c2p(0.15, 0), LEFT)
        zero_dot = Dot(ax2.c2p(z0, 0), radius=0.1, color=st.TEXT).set_stroke(st.BG, 3,
                                                                              background=True)
        zero_lab = VGroup(st.T("zero at", st.SMALL_SIZE, st.TEXT),
                          st.M(rf"k={z0:.2f}", size=36)).arrange(RIGHT, buff=0.12)
        zero_lab.next_to(ax2.c2p(z0, 1.12), UP, buff=0.05)

        kb, cb = d["k_bands"], d["ctf_at_bands"]
        bmarks = VGroup(*[DashedLine(ax2.c2p(kk, -1.1), ax2.c2p(kk, 1.1), dash_length=0.07,
                                     stroke_width=2, color=st.MUTED) for kk in kb])
        bnums = VGroup(*[st.M(rf"k_{i + 1}", size=36, color=st.MUTED).next_to(
            ax2.c2p(kk, -1.1), DOWN, buff=0.06) for i, kk in enumerate(kb)])

        # the band has width 0.08 around the zero and the CTF changes sign across it, so a
        # small fraction survives: measure it on the displayed window and show it
        xw = (np.abs(x) <= 1.6)
        leak = [float(np.max(np.abs(d["bands_after"][i][xw]))
                      / np.max(np.abs(d["bands_before"][i][xw]))) for i in range(3)]
        verdict = ["reversed", "≈ removed", "kept"]
        vcol = [neg_col, st.TEXT, st.TRUTH]
        panels, before, after, plabs = VGroup(), VGroup(), VGroup(), VGroup()
        xs = (1.55, 3.6, 5.65)
        for i in range(3):
            bb, ba = d["bands_before"][i], d["bands_after"][i]
            nrm = float(np.max(np.abs(bb)))
            axb = Axes(x_range=[-1.6, 1.6, 1], y_range=[-1.15, 1.15, 1], x_length=1.8,
                       y_length=1.1, tips=False,
                       axis_config=dict(color=st.GRID, stroke_width=1.5, include_ticks=False))
            axb.move_to([xs[i], -1.5, 0])
            fr = Rectangle(width=1.95, height=1.3, stroke_color=st.GRID,
                           stroke_width=1.5).move_to(axb)
            panels.add(VGroup(fr, axb.x_axis))
            before.add(curve(axb, x[xw], bb[xw] / nrm, st.TRUTH, 3))
            after.add(curve(axb, x[xw], ba[xw] / nrm, st.TEXT, 3))
            top = st.M(rf"k_{i + 1}={kb[i]:.2f}", size=36).next_to(fr, UP, buff=0.1)
            bot = VGroup(st.M(fmt(float(cb[i])).replace("−", "-"), size=36),
                         st.T(verdict[i], st.SMALL_SIZE + 2, vcol[i], weight="BOLD"))
            if i == 1:
                bot.add(st.T(f"({100 * leak[1]:.0f}% left)", st.SMALL_SIZE, st.MUTED))
            bot.arrange(DOWN, buff=0.1).next_to(fr, DOWN, buff=0.1)
            plabs.add(VGroup(top, bot))
        ctf_row = st.T("CTF:", st.SMALL_SIZE, st.MUTED).next_to(plabs[0][1][0], LEFT, buff=0.2)
        # legend for the panels (and for green p(x) / white s(x) on the left), two lines at 26
        band_note = st.T("green = before CTF, white = after\n(each band normalised for display)",
                         st.SMALL_SIZE, st.MUTED)
        band_note.to_corner(DL, buff=st.MARGIN)

        # ---- animate
        self.play(FadeIn(head), FadeIn(Group(inset, inset_fr, inset_lab)),
                  FadeIn(one_d_grp), Create(ax1), FadeIn(ax1_lab), run_time=0.9)
        self.play(Create(p_curve), FadeIn(p_lab), run_time=0.9)
        self.play(Create(ax2), FadeIn(ctf_lab), FadeIn(k_lab), FadeIn(ones), FadeIn(ctf_note),
                  Create(ctf_curve, rate_func=linear), run_time=1.1)
        self.play(FadeIn(lobe_poly), FadeIn(lobe_lab), FadeIn(zero_dot, scale=0.5),
                  FadeIn(zero_lab), run_time=0.7)
        self.wait(0.6)
        # the callout sits where the k_1 marker labels appear: clear it first (no ghosting)
        self.play(FadeOut(lobe_lab), run_time=0.25)
        self.play(FadeIn(bmarks), FadeIn(bnums), FadeIn(panels),
                  FadeIn(VGroup(*[pl[0] for pl in plabs])),
                  *[Create(b) for b in before], FadeIn(band_note), run_time=0.85)
        self.wait(0.2)
        # after the CTF: white over the dimmed green band
        self.play(*[b.animate.set_stroke(opacity=0.4) for b in before],
                  *[TransformFromCopy(before[i], after[i]) for i in range(3)],
                  FadeIn(VGroup(*[pl[1] for pl in plabs])), FadeIn(ctf_row), run_time=1.2)
        self.wait(0.3)
        # the same multiplication applied to the whole profile and to the 2D inset
        self.play(Transform(p_curve, s_curve),
                  Succession(FadeOut(p_lab, shift=0.1 * UP, run_time=0.35),
                             FadeIn(s_lab, shift=0.1 * UP, run_time=0.35)),
                  FadeIn(mag_lab), FadeIn(inset_ctf),
                  Succession(FadeOut(inset_lab, run_time=0.35), FadeIn(inset_ctf_lab,
                                                                        run_time=0.35)),
                  run_time=1.0)
        self.remove(inset)
        self.wait(0.7)
        self.b2 = dict(head=head, inset=inset_ctf, inset_fr=inset_fr, inset_lab=inset_ctf_lab,
                       one_d=one_d, one_d_eq=one_d_eq, ax1=ax1, ax1_lab=ax1_lab, s_curve=p_curve,
                       s_lab=s_lab, mag_lab=mag_lab)
        self.play(FadeOut(VGroup(ax2, ctf_curve, ctf_lab, k_lab, ones, lobe_poly, zero_dot,
                                 zero_lab, bmarks, bnums, panels, before, after, plabs,
                                 band_note, ctf_row, ctf_note)), run_time=0.5)

    # ================================================================ beat 3: 16–23.5 s
    def beat3_noise(self, d):
        b = self.b2
        x, s, eps, I = d["x"], d["s_ctf"], d["eps"], d["I"]
        snr1 = float(np.std(s) / d["sigma_1d"])
        snr2 = float(np.std(d["img_ctf"]) / d["sigma_2d"])

        head = VGroup(st.T("Noise:", st.BODY_SIZE, st.TEXT),
                      st.M(r"I(x)=s(x)+\epsilon(x)", size=40)).arrange(RIGHT, buff=0.25)
        head.to_corner(UL, buff=st.MARGIN)

        yl = 1.08 * float(np.max(np.abs(I)))
        rows, labs = [], []
        ys = (1.6, -0.2, -2.0)
        for yy in ys:
            ax = Axes(x_range=[-2, 2, 1], y_range=[-yl, yl, yl], x_length=6.6, y_length=1.4,
                      tips=False, axis_config=dict(color=st.GRID, stroke_width=2,
                                                   include_ticks=False))
            ax.move_to([-1.45, yy, 0])
            rows.append(ax)
        cs = curve(rows[0], x, s, st.TEXT, 4)
        ce = curve(rows[1], x, eps, st.RESID, 3)
        ci = curve(rows[2], x, I, st.OBS, 3)
        names = [(r"s(x)", st.TEXT), (r"\epsilon(x)", st.RESID), (r"I(x)", st.OBS)]
        for ax, (tx, col) in zip(rows, names):
            labs.append(st.M(tx, size=40, color=col).next_to(ax, LEFT, buff=0.35))
        # operators on the curve rows, left of the labels, as in a written column sum
        ops = VGroup(st.M("+", size=44).next_to(labs[1], LEFT, buff=0.15),
                     st.M("=", size=44).next_to(labs[2], LEFT, buff=0.15))
        sub = [st.T("CTF-modulated signal, true scale", st.SMALL_SIZE, st.MUTED),
               st.T("independent, zero-mean Gaussian (idealised)", st.SMALL_SIZE, st.RESID),
               st.T(f"1D: amplitude SNR = std(s) / σ = {snr1:.2f}", st.SMALL_SIZE, st.OBS)]
        for ax, t in zip(rows, sub):
            t.next_to(ax, DOWN, buff=0.06).align_to(ax, LEFT)

        one_d = st.T("1D illustration (one vertical scale)", st.BODY_SIZE - 2, st.TEXT,
                     slant="ITALIC")
        one_d.next_to(rows[0], UP, buff=0.12).align_to(rows[0], RIGHT)

        # 2D inset on the right; the noise level steps through a precomputed series
        frames = d["noise_frames"]
        alphas = d["noise_alphas"]
        vmin, vmax = float(frames.min()), float(frames.max())
        H = 2.7
        ims = Group(*[st.heatmap(f, height=H, vmin=vmin, vmax=vmax).move_to([4.6, 0.25, 0])
                      for f in frames])
        fr = st.framed(ims[0])
        im_top = st.T("2D, after CTF (rescaled)", st.SMALL_SIZE, st.MUTED).next_to(fr, UP,
                                                                                buff=0.1)
        assert im_top.get_right()[0] <= 7.11 - st.MARGIN, im_top.get_right()
        tr = ValueTracker(0)
        n_fr = len(frames)

        def idx():
            return int(np.clip(np.floor(tr.get_value() + 0.5), 0, n_fr - 1))

        for j, im in enumerate(ims):
            im.set_opacity(1.0 if j == 0 else 0.0)
            im.add_updater(lambda m, j=j: m.set_opacity(1.0 if j == idx() else 0.0))
        a_txt = st.T("noise scale", st.SMALL_SIZE + 2, st.TEXT)
        a_sym = st.M(r"\alpha=", size=36)
        a_num = DecimalNumber(0.0, num_decimal_places=1, font_size=36, color=st.TEXT)
        a_row = VGroup(a_txt, a_sym, a_num).arrange(RIGHT, buff=0.15).next_to(fr, DOWN, buff=0.15)
        a_num.add_updater(lambda m: m.set_value(float(alphas[idx()])))
        # noise_frames[-1] has alpha = 1, i.e. noise sigma_2d: the SNR holds at alpha = 1
        assert np.isclose(alphas[-1], 1.0)
        snr2_lab = st.T(f"amplitude SNR = {snr2:.2f}\n(2D, at α = 1)", st.SMALL_SIZE + 2,
                        st.OBS)
        snr2_lab.next_to(a_row, DOWN, buff=0.12)

        # transition: keep s(x) and move it into the top row
        self.play(FadeOut(Group(b["head"], b["one_d"], b["one_d_eq"], b["mag_lab"],
                                b["inset_fr"], b["inset_lab"], b["inset"], b["ax1_lab"],
                                b["s_lab"])), run_time=0.4)
        self.play(FadeIn(head), Transform(b["ax1"], rows[0]), Transform(b["s_curve"], cs),
                  FadeIn(labs[0]), run_time=0.8)
        self.add(ims)
        self.play(FadeIn(sub[0]), FadeIn(one_d), FadeIn(fr), FadeIn(im_top),
                  FadeIn(a_row), run_time=0.6)
        self.play(Create(rows[1]), FadeIn(ops[0]), FadeIn(labs[1]), FadeIn(sub[1]),
                  run_time=0.5)
        # eps(x) is drawn while the 2D noise level steps up through the series
        self.play(Create(ce), tr.animate.set_value(n_fr - 1), run_time=1.4, rate_func=linear)
        for im in ims:
            im.clear_updaters()
        a_num.clear_updaters()
        self.play(FadeIn(snr2_lab), Create(rows[2]), FadeIn(ops[1]), FadeIn(labs[2]),
                  run_time=0.5)
        self.play(Create(ci), run_time=1.0)
        self.play(FadeIn(sub[2]), run_time=0.5)
        self.wait(1.1)
        self.play(FadeOut(Group(head, b["ax1"], b["s_curve"], labs[0], rows[1], rows[2], ce,
                                ci, labs[1], labs[2], ops, *sub, one_d, ims, fr,
                                im_top, a_row, snr2_lab)), run_time=0.5)

    # ================================================================ beat 4: 23.5–34 s
    def beat4_statement(self, d):
        head = st.T("Fourier-space observation model", st.BODY_SIZE + 4, st.TEXT, weight="BOLD")
        head.to_corner(UL, buff=st.MARGIN)
        clear_of_tag(head, self.tag_left)
        tex = (r"\widetilde I_i(\mathbf k)=\mathrm{CTF}_i(\mathbf k)\,\widehat{P_{R_i}\rho}"
               r"(\mathbf k)\,e^{-2\pi i\,\mathbf k\cdot\mathbf t_i}+\widetilde\epsilon_i(\mathbf k)")
        eq = st.M(tex, size=56)
        eq.move_to([0, 1.5, 0])
        g = eq[0]
        assert len(g) == 38, len(g)
        note = st.simplified_note("simplified, pedagogical model: no envelopes, no detector "
                                  "response")
        note.next_to(eq, DOWN, buff=0.3)

        # glyph groups (indices checked against the compiled formula)
        terms = [
            ("projection", "line integral along the viewing direction",
             [14, 15, 18, 19, 20, 21], st.TRUTH),
            ("orientation", "a rotation in SO(3), unknown per image", [16, 17], st.EST),
            ("CTF", "oscillating transfer with zeros and sign changes",
             list(range(7, 14)), st.TEXT),
            ("translation", "in-plane shift = phase ramp, unknown per image",
             list(range(22, 31)), st.EST),
            ("noise", "independent, zero-mean; real noise is more complex",
             list(range(32, 38)), st.RESID),
        ]
        rows = VGroup()
        for name, desc, _, col in terms:
            rows.add(VGroup(st.T(name, st.BODY_SIZE, col, weight="BOLD"),
                            st.T(desc, st.BODY_SIZE, st.TEXT)))
        name_w = max(r[0].width for r in rows)
        for r in rows:
            r[1].next_to(r[0], RIGHT, buff=0.3 + name_w - r[0].width)
        rows.arrange(DOWN, aligned_edge=LEFT, buff=0.22).move_to([0, -1.25, 0])

        g.set_opacity(0.35)
        self.play(FadeIn(head), FadeIn(eq), FadeIn(note), run_time=0.8)
        box, ctf_ul = None, VGroup()
        for (name, _, idx, col), r in zip(terms, rows):
            sub = VGroup(*[g[i] for i in idx])
            nb = SurroundingRectangle(sub, color=col, buff=0.07, stroke_width=3)
            anims = [sub.animate.set_opacity(1).set_color(col), FadeIn(r, shift=0.1 * RIGHT)]
            anims.append(Create(nb) if box is None else ReplacementTransform(box, nb))
            if name == "CTF":
                # white CTF glyphs look like the unlit ones once the box moves on: underline
                ctf_ul = Underline(sub, color=st.TEXT, buff=0.05, stroke_width=3)
                anims.append(Create(ctf_ul))
            self.play(*anims, run_time=0.45)
            box = nb
            self.wait(0.5)
        self.play(FadeOut(box), g.animate.set_opacity(1), run_time=0.4)
        self.wait(0.7)

        # takeaway freeze (3-5 s)
        lines = ["Projection integrates along the view: depth is lost.",
                 "The CTF scales, zeroes or sign-reverses each frequency.",
                 "Noise adds independent fluctuations; one image is not enough.",
                 "Inverse problem: recover a common source from many observations."]
        card = st.takeaway_card(lines, math=tex)
        card.move_to([0, -0.2, 0])
        self.play(FadeOut(VGroup(head, eq, note, rows, ctf_ul)), run_time=0.4)
        self.play(FadeIn(card), run_time=0.5)
        self.wait(3.2)
