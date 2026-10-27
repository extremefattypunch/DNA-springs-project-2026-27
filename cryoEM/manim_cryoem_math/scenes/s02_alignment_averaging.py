"""Clip 02 — 2D alignment, averaging and the square-root law.

Storyboard beats (spec section 3, Scene 02) -> code sections:
  beat 1  (0-8 s)   twelve noisy observations of the same 4-dot pattern; one is
                    enlarged with its rotation arc phi_i and shift arrow t_i
  beat 2  (8-15 s)  the unaligned mean smears the dots; the clean reference S
                    stays fixed beside it
  beat 3  (15-26 s) undo (phi_i, t_i) in every grid; the aligned mean is sharp;
                    red/green residual-versus-signal comparison
  beat 4  (26-37 s) aligned means for N = 1, 4, 16, 64, 256 and amplitude SNR
                    versus N: theory 0.5 sqrt(N) curve, measured markers
  beat 5  (37-42 s) freeze: unaligned vs aligned and the square-root law

All numbers come from cached("alignment", ne.alignment) (poses known exactly,
i.i.d. zero-mean Gaussian noise) or from local computations on those arrays.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, projections as pj
from helpers import point_constellation as pc
from helpers.cache import cached

N_PIX = 64
SEL = 2                      # the observation drawn enlarged with its arc and arrow
TILE_H = 1.2                 # small grid size (scene units)
TILE_GAP = 0.12


# ---------------------------------------------------------------- local numbers
def _clip02_numbers() -> dict:
    """Derived quantities shown on screen (all computed from the cached data)."""
    d = cached("alignment", ne.alignment)
    S, un, al = d["S"], d["unaligned_mean"], d["aligned_mean"]
    rms = lambda r: float(np.sqrt(np.mean(r ** 2)))
    # least-squares amplitude of S contained in a mean image: beta = <I, S> / <S, S>
    beta = lambda m: float(np.sum(m * S) / np.sum(S * S))
    return dict(
        rms_un=np.array([rms(m - S) for m in un]),
        rms_al=np.array([rms(m - S) for m in al]),
        beta_un=beta(un[-1]), beta_al=beta(al[-1]),
        peak_S=float(S.max()), peak_un=float(un[-1].max()), peak_al=float(al[-1].max()),
        std_S=float(np.std(S)),
    )


# ---------------------------------------------------------------- drawing helpers
def _img_to_scene(center, h, extent):
    """World (x, y) in image units -> scene point inside a heatmap of height h.

    Pixel i sits at (i - n/2) * dx, so the pixel block spans
    [-extent - dx/2, extent - dx/2] and its centre is at -dx/2.
    """
    dx = 2 * extent / N_PIX
    c = np.asarray(center, float)

    def f(xy):
        xy = np.asarray(xy, float)
        return c + np.array([(xy[0] + dx / 2) * h / (2 * extent),
                             (xy[1] + dx / 2) * h / (2 * extent), 0.0])
    return f


def _rgb_image(rgb01: np.ndarray, height: float) -> ImageMobject:
    """[x, y, 3] array in [0, 1] -> nearest-neighbour ImageMobject."""
    a = np.flipud(np.swapaxes(rgb01, 0, 1))
    im = ImageMobject((np.clip(a, 0, 1) * 255).astype(np.uint8))
    im.set_resampling_algorithm(RESAMPLING_ALGORITHMS["nearest"])
    im.height = height
    return im


def _residual_composite(mean, S, height):
    """Green = matched part beta*S of the mean, red = |mean - beta*S| (x3), scale of S.

    beta = <mean, S> / <S, S> is the least-squares amplitude of S in the mean.
    """
    vmax = S.max()
    beta = float(np.sum(mean * S) / np.sum(S * S))
    g = np.clip(beta * S / vmax, 0, 1)
    r = np.clip(np.abs(mean - beta * S) / vmax * 3.0, 0, 1)   # residual x3 so it is visible
    tg, tr = st._hex_to_rgb(st.TRUTH) / 255, st._hex_to_rgb(st.RESID) / 255
    rgb = g[..., None] * tg + r[..., None] * tr
    return _rgb_image(np.clip(rgb, 0, 1), height)


def _crosshair(center, h, opacity=0.45):
    c = np.asarray(center, float)
    return VGroup(
        Line(c + LEFT * h / 2, c + RIGHT * h / 2, stroke_width=1.2, color=st.MUTED,
             stroke_opacity=opacity),
        Line(c + DOWN * h / 2, c + UP * h / 2, stroke_width=1.2, color=st.MUTED,
             stroke_opacity=opacity))


class Scene02AlignmentAveraging(Scene):
    def construct(self):
        d = cached("alignment", ne.alignment)    # precomputed arrays; nothing per frame
        x = cached("s02_alignment_averaging_numbers_v2", _clip02_numbers)
        S, sigma, ext = d["S"], float(d["sigma"]), float(d["extent"])
        phi, t = d["phi"], d["t"]
        pattern = d["pattern"]
        ns = d["ns"].astype(int)
        snr_meas, snr1 = d["snr_measured"], float(d["snr1"])

        tag = st.illustration_label()
        self.add(tag)

        # ============================================================ beat 1
        title = st.title_card("2D alignment and averaging",
                              "Why align before averaging, and how fast does the noise fall?")
        self.play(FadeIn(title), run_time=0.8)

        lo, hi = -2 * sigma, float(S.max()) + sigma      # one fixed scale for all 12 grids
        tiles, frames, overlays, centers, xh = [], VGroup(), [], [], VGroup()
        x0, y0 = -6.66 + TILE_H / 2, 1.85 - TILE_H / 2
        for i in range(12):
            r, c = divmod(i, 4)
            ctr = np.array([x0 + c * (TILE_H + TILE_GAP), y0 - r * (TILE_H + TILE_GAP), 0])
            centers.append(ctr)
            im = st.heatmap(d["obs12"][i], TILE_H, vmin=lo, vmax=hi).move_to(ctr)
            tiles.append(im)
            frames.add(Rectangle(width=TILE_H, height=TILE_H, stroke_color=st.OBS,
                                 stroke_width=1.5).move_to(ctr))
            xh.add(_crosshair(ctr, TILE_H, 0.35))
            f = _img_to_scene(ctr, TILE_H, ext)
            P = pj.transform_2d(pattern, phi[i], t[i])
            dots = VGroup(*[Dot(f(p), radius=0.035, color=st.TRUTH) for p in P])
            for dt in dots:
                dt.set_stroke(st.BG, width=2, background=True)
            overlays.append(dots)
        tile_lab = st.T(f"12 of N = {len(d['phi'])} observations", st.SMALL_SIZE,
                        st.MUTED)
        tile_lab.next_to(frames, DOWN, buff=0.18).align_to(frames, LEFT)
        self.play(LaggedStart(*[FadeIn(Group(tiles[i], frames[i], xh[i])) for i in range(12)],
                              lag_ratio=0.08), FadeIn(tile_lab), run_time=1.6)
        self.play(LaggedStart(*[FadeIn(o) for o in overlays], lag_ratio=0.05), run_time=0.8)

        # enlarged copy of one observation with its rotation arc and shift arrow
        hl = Rectangle(width=TILE_H + 0.08, height=TILE_H + 0.08, stroke_color=st.EST,
                       stroke_width=4).move_to(centers[SEL])
        BIG = 3.4
        bctr = np.array([1.05, 0.65, 0])
        big = st.heatmap(d["obs12"][SEL], BIG, vmin=lo, vmax=hi).move_to(bctr)
        bframe = Rectangle(width=BIG, height=BIG, stroke_color=st.OBS,
                           stroke_width=2).move_to(bctr)
        fb = _img_to_scene(bctr, BIG, ext)
        o_pt, t_pt = fb((0, 0)), fb(t[SEL])
        bxh = _crosshair(o_pt, BIG, 0.55)
        Pb = pj.transform_2d(pattern, phi[SEL], t[SEL])
        bdots = VGroup(*[Dot(fb(p), radius=0.07, color=st.TRUTH) for p in Pb])
        for dt in bdots:
            dt.set_stroke(st.BG, width=3, background=True)
        ang = np.radians(phi[SEL])
        L = 1.05
        ref_ax = DashedLine(t_pt, t_pt + L * RIGHT, color=st.EST, stroke_width=3)
        body_ax = Arrow(t_pt, t_pt + L * np.array([np.cos(ang), np.sin(ang), 0]), buff=0,
                        color=st.EST, stroke_width=4, max_tip_length_to_length_ratio=0.22)
        arc = Arc(radius=0.45, start_angle=0, angle=ang, arc_center=t_pt, color=st.EST,
                  stroke_width=4)
        for m in (ref_ax, body_ax, arc):
            m.set_stroke(st.BG, width=9, background=True)
        # the shift is short (|t_i| ~ 0.3 scene units): thin outline, small head
        t_arrow = Arrow(o_pt, t_pt, buff=0, color=st.EST, stroke_width=6, tip_length=0.12,
                        max_tip_length_to_length_ratio=0.3, max_stroke_width_to_length_ratio=30)
        t_arrow.set_stroke(st.BG, width=3, background=True)
        o_dot = Dot(o_pt, radius=0.05, color=st.MUTED).set_stroke(st.BG, width=2,
                                                                 background=True)
        t_arrow.set_z_index(2)                  # keep the short shift above the body axis
        o_dot.set_z_index(2)
        geo = VGroup(ref_ax, body_ax, arc, t_arrow, o_dot)
        phi_lab = st.M(r"\phi_i", size=38, color=st.EST)
        phi_lab.move_to(t_pt + 1.0 * np.array([np.cos(np.radians(20)), np.sin(np.radians(20)), 0]))
        t_lab = st.M(r"\mathbf t_i", size=38, color=st.EST)
        t_lab.next_to(t_pt, DOWN, buff=0.14)       # below the tip: clear of the body axis
        for m in (phi_lab, t_lab):
            m.add_background_rectangle(color=st.BG, opacity=0.85, buff=0.04)

        ti = t[SEL]
        vals = VGroup(
            st.M(rf"\phi_i = {phi[SEL]:.0f}^\circ", size=38, color=st.EST),
            st.M(rf"\mathbf t_i = ({ti[0]:.2f},\,{ti[1]:.2f})", size=38, color=st.EST),
            st.M(rf"\sigma = {sigma:.2f}", size=38, color=st.RESID),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.28)
        vals.next_to(bframe, RIGHT, buff=0.4).align_to(bframe, UP).shift(DOWN * 0.2)
        vnote = VGroup(st.T("signal positions", st.SMALL_SIZE, st.TRUTH),
                       st.T("(known in this toy)", st.SMALL_SIZE, st.MUTED)
                       ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        vnote.next_to(vals, DOWN, buff=0.45, aligned_edge=LEFT)
        eq1 = st.M(r"I_i(\mathbf x)=S\big(R_{\phi_i}^{-1}(\mathbf x-\mathbf t_i)\big)"
                   r"+\epsilon_i(\mathbf x)", size=42)
        eq1.move_to([3.0, -1.65, 0])
        eq1_note = st.T("rotate by φ (ψ in 3D), shift, add i.i.d. noise",
                        st.SMALL_SIZE, st.MUTED).next_to(eq1, DOWN, buff=0.22)
        eq1_simp = st.simplified_note("toy 2D model: no CTF, one fixed view").next_to(
            eq1_note, DOWN, buff=0.12)

        self.play(Create(hl), run_time=0.4)
        self.play(FadeIn(Group(big, bframe, bxh), target_position=centers[SEL], scale=0.45),
                  FadeIn(bdots), run_time=1.0)
        self.play(FadeIn(o_dot), GrowArrow(t_arrow), FadeIn(t_lab), run_time=0.7)
        self.play(Create(ref_ax), Create(arc), GrowArrow(body_ax), FadeIn(phi_lab),
                  run_time=0.9)
        self.play(FadeIn(vals), FadeIn(vnote), run_time=0.6)
        self.play(Write(eq1), FadeIn(eq1_note), FadeIn(eq1_simp), FadeOut(title[1]),
                  run_time=1.0)
        self.wait(1.2)
        beat1_right = Group(big, bframe, bxh, bdots, geo, phi_lab, t_lab, vals, vnote, eq1,
                            eq1_note, eq1_simp, hl)

        # ============================================================ beat 2
        SH = 1.9
        slot_x = [0.15, 2.9, 5.65]
        slot_y = 0.6
        slots = [np.array([sx, slot_y, 0]) for sx in slot_x]

        def head(text, ctr, color=st.TEXT):
            return st.T(text, st.SMALL_SIZE, color).next_to(
                Rectangle(width=SH, height=SH).move_to(ctr), UP, buff=0.14)

        un_img = st.heatmap(d["unaligned_mean"][-1], SH).move_to(slots[0])
        un_fr = Rectangle(width=SH, height=SH, stroke_color=st.GRID, stroke_width=2
                          ).move_to(slots[0])
        un_head = head("unaligned mean", slots[0])
        ref_img = st.heatmap(S, SH).move_to(slots[2])
        ref_fr = Rectangle(width=SH, height=SH, stroke_color=st.TRUTH, stroke_width=2
                           ).move_to(slots[2])
        fr = _img_to_scene(slots[2], SH, ext)
        ref_dots = VGroup(*[Dot(fr(p), radius=0.045, color=st.TRUTH).set_stroke(
            st.BG, width=2, background=True) for p in pattern])
        ref_head = head("reference S", slots[2], st.TRUTH)
        ref_head.align_to(un_head, UP)          # same top line despite no descender
        ph = Rectangle(width=SH, height=SH, stroke_color=st.GRID, stroke_width=2
                       ).move_to(slots[1])
        ph_q = st.T("?", 40, st.MUTED).move_to(slots[1])
        n_lab = st.T(f"N = {ns[-1]}; each panel on its own grey scale", st.SMALL_SIZE,
                     st.MUTED)
        n_lab.move_to([slot_x[1], slot_y - SH / 2 - 0.3, 0])
        un_val = st.M(rf"\text{{unaligned peak }}{x['peak_un']:.2f}\ \text{{vs}}\ "
                      rf"{x['peak_S']:.2f}\text{{ in }}S", size=36, color=st.TEXT)
        un_val.move_to([slot_x[1], -1.3, 0])
        eq_naive = st.M(r"\overline I_N^{\,\mathrm{naive}}=\frac1N\sum_{i=1}^N I_i", size=42)
        eq_naive.move_to([slot_x[1], -2.3, 0])
        cap2 = st.caption("Averaging without alignment smears the dots into rings")

        self.play(FadeOut(beat1_right), run_time=0.6)
        # flyers travel underneath the tiles and text (z = -1) and shrink on the way
        flyers = [tiles[i].copy().set_z_index(-1) for i in range(12)]
        self.play(LaggedStart(*[f.animate.move_to(slots[0]).scale(0.5).set_opacity(0.0)
                                for f in flyers], lag_ratio=0.05),
                  FadeIn(Group(un_fr, un_head)), run_time=1.2)
        self.remove(*flyers)
        self.play(FadeIn(un_img), FadeIn(Group(ref_img, ref_fr, ref_dots, ref_head)),
                  FadeIn(Group(ph, ph_q)), run_time=0.8)
        self.play(Write(eq_naive), FadeIn(n_lab), FadeIn(cap2), run_time=0.9)
        self.play(FadeIn(un_val), run_time=0.5)
        self.wait(2.2)

        # ============================================================ beat 3
        cap3 = st.caption("Align: undo each rotation and shift (poses known exactly here)")
        self.play(FadeOut(un_val), FadeOut(eq_naive),
                  Succession(FadeOut(cap2, shift=0.1 * UP), FadeIn(cap3, shift=0.1 * UP)),
                  run_time=0.8)

        def undo_pose(mob, a):
            """Dots follow R_{(1-a) phi_i} x + (1-a) t_i: a = 0 observed, a = 1 aligned."""
            for i, o in enumerate(overlays):
                f = _img_to_scene(centers[i], TILE_H, ext)
                P = pj.transform_2d(pattern, (1 - a) * phi[i], (1 - a) * t[i])
                for dt, p in zip(o, P):
                    dt.move_to(f(p))

        al_tiles = [st.heatmap(d["aligned12"][i], TILE_H, vmin=lo, vmax=hi).move_to(centers[i])
                    for i in range(12)]
        dot_group = VGroup(*overlays)
        # draw order: images (0) < frames and crosshairs (1) < moving dots (3)
        for im in tiles + al_tiles:
            im.set_z_index(0)
        frames.set_z_index(1)
        xh.set_z_index(1)
        dot_group.set_z_index(3)
        # (a) the dots visibly rotate back by -phi_i and shift by -t_i over the observed grids
        self.play(UpdateFromAlphaFunc(dot_group, undo_pose), run_time=2.0, rate_func=smooth)
        # (b) then each grid's image is replaced by its aligned version (no cross-fade)
        #     as a sweep across the grid, so the whole grid never blanks at once
        self.play(LaggedStart(*[Succession(FadeOut(tiles[i], run_time=0.3),
                                           FadeIn(al_tiles[i], run_time=0.4))
                                for i in range(12)], lag_ratio=0.06), run_time=1.0)
        al_lab = st.T("12 grids, one common frame", st.SMALL_SIZE, st.MUTED)
        al_lab.move_to(tile_lab, aligned_edge=LEFT)
        self.play(Succession(FadeOut(tile_lab, shift=0.1 * UP), FadeIn(al_lab, shift=0.1 * UP)),
                  run_time=0.6)

        al_img = st.heatmap(d["aligned_mean"][-1], SH).move_to(slots[1])
        al_fr = Rectangle(width=SH, height=SH, stroke_color=st.GRID, stroke_width=2
                          ).move_to(slots[1])
        al_head = head("aligned mean", slots[1])
        al_head.align_to(un_head, UP)
        flyers = [al_tiles[i].copy().set_z_index(-1) for i in range(12)]
        self.play(FadeOut(ph_q), LaggedStart(*[f.animate.move_to(slots[1]).scale(0.5)
                                               .set_opacity(0.0)
                                               for f in flyers], lag_ratio=0.05),
                  run_time=1.1)
        self.remove(*flyers)
        self.play(FadeIn(al_img), FadeIn(al_head), FadeOut(ph), FadeIn(al_fr), run_time=0.6)
        eq_al = st.M(r"\overline I_N=\frac1N\sum_{i=1}^N I_i^{\mathrm{aligned}}", size=42)
        eq_al.move_to([slot_x[1], -1.75, 0])
        self.play(Write(eq_al), run_time=0.8)
        self.wait(0.6)

        # red/green: matched part beta*S of each mean (green) versus residual (red),
        # beta = <mean, S> / <S, S> (least-squares amplitude of S in that mean)
        un_cmp = _residual_composite(d["unaligned_mean"][-1], S, SH).move_to(slots[0])
        al_cmp = _residual_composite(d["aligned_mean"][-1], S, SH).move_to(slots[1])
        cmp_y = slot_y - SH / 2 - 0.32          # the row n_lab used (n_lab fades out)
        cmp_row = VGroup(
            st.M(rf"\beta={x['beta_un']:.2f}", size=36, color=st.TRUTH).move_to(
                [slot_x[0], cmp_y, 0]),
            st.M(rf"\beta={x['beta_al']:.2f}", size=36, color=st.TRUTH).move_to(
                [slot_x[1], cmp_y, 0]),
        )
        leg = st.legend([(st.TRUTH, "matched part βS"), (st.RESID, "|mean − βS| (×3)")],
                        size=st.SMALL_SIZE)
        b_note = st.T("β: least-squares amplitude of S", st.SMALL_SIZE, st.MUTED)
        leg_blk = VGroup(leg, b_note).arrange(DOWN, aligned_edge=LEFT, buff=0.16)
        leg_blk.next_to(cmp_row, DOWN, buff=0.3).align_to(un_fr, LEFT)
        cap3b = st.caption("Aligned: residual is only noise. Unaligned: signal leaks into it.")
        self.play(FadeTransform(un_img, un_cmp), FadeTransform(al_img, al_cmp), FadeOut(eq_al),
                  FadeOut(n_lab),
                  Succession(FadeOut(cap3, shift=0.1 * UP), FadeIn(cap3b, shift=0.1 * UP)),
                  run_time=1.0)
        self.play(FadeIn(cmp_row), FadeIn(leg_blk), run_time=0.7)
        self.wait(2.0)

        # ============================================================ beat 4
        beat3 = Group(*al_tiles, frames, xh, dot_group, al_lab, un_cmp, un_fr, un_head,
                      al_cmp, al_fr, al_head, ref_img, ref_fr, ref_dots, ref_head,
                      cmp_row, leg_blk, cap3b)
        self.play(FadeOut(beat3), run_time=0.7)

        BH = 3.0
        bc = np.array([-4.35, 0.1, 0])
        means = [st.heatmap(m, BH).move_to(bc) for m in d["aligned_mean"]]
        mfr = Rectangle(width=BH, height=BH, stroke_color=st.GRID, stroke_width=2).move_to(bc)
        mhead = VGroup(st.T("aligned mean", st.SMALL_SIZE + 2),
                       st.T("(own grey scale per panel)", st.SMALL_SIZE, st.MUTED)
                       ).arrange(DOWN, buff=0.08).next_to(mfr, UP, buff=0.15)

        def n_tex(N):
            return st.M(rf"N = {N}", size=40).next_to(mfr, DOWN, buff=0.2)

        def snr_tex(k):
            return st.M(rf"\text{{measured }}\mathrm{{SNR}}_{{\rm amp}} = {snr_meas[k]:.2f}",
                        size=38, color=st.OBS).next_to(n_tex(1), DOWN, buff=0.16)

        snr_def = st.T("measured: std(S) / std(mean − S)", st.SMALL_SIZE, st.MUTED)
        snr_def.next_to(snr_tex(0), DOWN, buff=0.12)
        # wider than the image column: left-align with the readout (inside the margin)
        snr_def.align_to(np.array([min(snr_tex(k).get_left()[0] for k in range(len(ns))),
                                   0, 0]), LEFT)

        # plot: SNR versus N on a base-4 log axis (u = log_4 N), linear SNR axis
        ax = Axes(x_range=[0, 4.3, 1], y_range=[0, 9, 2], x_length=6.6, y_length=4.0,
                  tips=False, axis_config=dict(color=st.GRID, stroke_width=2,
                                               include_ticks=True))
        ax.move_to([2.95, 0.2, 0])
        xt = VGroup(*[st.T(str(N), st.SMALL_SIZE, st.MUTED).next_to(ax.c2p(np.log(N) / np.log(4),
                                                                         0), DOWN, buff=0.15)
                      for N in ns])
        yt = VGroup(*[st.T(str(v), st.SMALL_SIZE, st.MUTED).next_to(ax.c2p(0, v), LEFT,
                                                                   buff=0.15)
                      for v in (0, 2, 4, 6, 8)])
        xl = st.T("N (log scale: ×4 per step)", st.SMALL_SIZE + 2, st.MUTED)
        xl.next_to(xt, DOWN, buff=0.16)
        yl = st.M(r"\mathrm{SNR}_{\rm amp}", size=36, color=st.MUTED)
        yl.next_to(ax.c2p(0, 9), UP, buff=0.12).align_to(ax, LEFT).shift(LEFT * 0.35)
        theory = ax.plot(lambda u: snr1 * np.sqrt(4.0 ** u), x_range=[0, 4.15],
                         color=st.TEXT, stroke_width=4)
        th_lab = st.M(rf"{snr1:.1f}\sqrt N", size=38).next_to(
            ax.c2p(3.0, snr1 * 2 ** 3.0), DOWN + RIGHT, buff=0.15)
        law = st.M(r"\mathrm{SNR}_{\rm amp}=\frac{\operatorname{std}(S)}"
                   r"{\sigma/\sqrt N}\propto\sqrt N", size=38)
        law.next_to(ax.c2p(0, 9), DOWN + RIGHT, buff=0.0).shift(RIGHT * 0.35 + DOWN * 0.05)
        pleg = st.legend([(st.TEXT, "theory"), (st.OBS, "measured")], size=st.SMALL_SIZE)
        pleg.next_to(ax.c2p(4.3, 0), UP + LEFT, buff=0.0).shift(UP * 0.25 + LEFT * 0.05)
        cap4 = st.caption("Assumes identical signal, exact poses (known here), i.i.d. "
                          "zero-mean Gaussian noise")
        cap4b = st.caption("Toy noise model: real detector (counting) noise is more complicated")

        u_tr = ValueTracker(0.0)
        # hollow ring in the theory colour that rides along the theory curve
        marker = always_redraw(lambda: Circle(radius=0.1, color=st.TEXT, stroke_width=3)
                               .set_stroke(st.BG, width=7, background=True)
                               .move_to(ax.c2p(u_tr.get_value(),
                                               snr1 * np.sqrt(4 ** u_tr.get_value())))
                               .set_z_index(3))           # ring stays above measured dots
        cur_img, cur_n, cur_s = means[0], n_tex(ns[0]), snr_tex(0)
        m0 = Dot(ax.c2p(0, snr_meas[0]), radius=0.08, color=st.OBS).set_z_index(2)
        self.play(FadeIn(Group(cur_img, mfr, mhead)), FadeIn(cur_n), FadeIn(cur_s),
                  FadeIn(snr_def), Create(ax), FadeIn(xt), FadeIn(yt), FadeIn(xl), FadeIn(yl),
                  FadeIn(cap4), run_time=1.0)
        self.play(Create(theory), FadeIn(th_lab), FadeIn(pleg), run_time=1.0)
        self.play(FadeIn(marker), FadeIn(m0), FadeIn(law), run_time=0.7)

        def swap(a, b):
            """Label swap: out 0.25 s, in 0.35 s, then hold (total 1.2 s)."""
            return Succession(FadeOut(a, shift=0.08 * UP, run_time=0.25),
                              FadeIn(b, shift=0.08 * UP, run_time=0.35),
                              Wait(run_time=0.6))

        for k in range(1, len(ns)):
            u = np.log(ns[k]) / np.log(4)
            nn, ss = n_tex(ns[k]), snr_tex(k)
            mk = Dot(ax.c2p(u, snr_meas[k]), radius=0.08, color=st.OBS).set_z_index(2)
            extra = [swap(cap4, cap4b)] if k == 2 else []
            self.play(u_tr.animate.set_value(u), FadeTransform(cur_img, means[k]),
                      swap(cur_n, nn), swap(cur_s, ss), *extra, run_time=1.2)
            self.play(FadeIn(mk, scale=1.2), run_time=0.4)
            cur_img, cur_n, cur_s = means[k], nn, ss
        self.wait(1.8)

        # ============================================================ beat 5
        marker.clear_updaters()
        self.play(*[FadeOut(m) for m in self.mobjects if m is not tag and m is not title[0]],
                  run_time=0.7)

        card = st.takeaway_card(
            ["Align first, then average: unaligned means smear the signal.",
             "Assumes same signal, exact poses, i.i.d. zero-mean noise;",
             "real detector (counting) noise is more complicated.",
             "Different viewing directions must not be averaged together:",
             "2D classification groups similar views while aligning them."],
            math=r"\operatorname{Std}(\overline\epsilon_N)=\frac{\sigma}{\sqrt N}"
                 r"\;\Rightarrow\;\mathrm{SNR}_{\rm amp}\propto\sqrt N", line_spacing=0.1)
        card.to_edge(DOWN, buff=st.MARGIN)
        gap = title[0].get_bottom()[1] - card.get_top()[1]
        FH = min(1.7, gap - 0.3)
        fy = (card.get_top()[1] + title[0].get_bottom()[1]) / 2
        fc_un, fc_al = np.array([-1.3, fy, 0]), np.array([1.3, fy, 0])
        # one shared grey scale, so the unaligned ring shows its true (low) amplitude
        lo_f = float(min(d["unaligned_mean"][-1].min(), d["aligned_mean"][-1].min()))
        hi_f = x["peak_al"]
        f_un = st.heatmap(d["unaligned_mean"][-1], FH, vmin=lo_f, vmax=hi_f).move_to(fc_un)
        f_al = st.heatmap(d["aligned_mean"][-1], FH, vmin=lo_f, vmax=hi_f).move_to(fc_al)
        f_fr = VGroup(Rectangle(width=FH, height=FH, stroke_color=st.GRID, stroke_width=2
                                ).move_to(fc_un),
                      Rectangle(width=FH, height=FH, stroke_color=st.GRID, stroke_width=2
                                ).move_to(fc_al))
        f_lab = VGroup(st.T(f"unaligned, N = {ns[-1]}", st.SMALL_SIZE + 2, st.MUTED
                            ).next_to(f_fr[0], LEFT, buff=0.25),
                       st.T(f"aligned, N = {ns[-1]}", st.SMALL_SIZE + 2, st.MUTED
                            ).next_to(f_fr[1], RIGHT, buff=0.25))
        f_note = st.T("one shared grey scale", st.SMALL_SIZE, st.MUTED)
        f_note.next_to(f_lab[0], DOWN, buff=0.12).align_to(f_lab[0], RIGHT)
        self.play(FadeIn(Group(f_un, f_al, f_fr, f_lab, f_note)), FadeIn(card), run_time=1.0)
        self.wait(4.5)
