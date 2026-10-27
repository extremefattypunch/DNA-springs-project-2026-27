"""Clip 02 — 2D alignment, averaging and the square-root-N law (~35 s).

Storyboard beats (spec §3, Scene 02) -> code sections:
    beat 1  0-7 s    12 particle windows: one 2D view S, random phi_i, t_i, i.i.d. noise
    beat 2  7-13 s   naive (unaligned) mean of the same 256 particles: a smear, not S
    beat 3  13-22 s  estimated phi_i, t_i removed: windows slide/rotate into one frame
    beat 4  22-29 s  aligned means for N = 1, 4, 16, 64, 256 + measured SNR(N) vs sqrt(N)
    beat 5  29-35 s  freeze: unaligned vs aligned (same N) + sqrt(N) law + takeaway

Every image, SNR point and correlation shown is computed in helpers/demos.alignment
(cached); nothing is hand-drawn. SNR convention (amplitude SNR): std(S) over the box /
noise sigma. Noise in the mean falls as sigma/sqrt(N) whether or not the copies are
aligned; alignment decides whether the signal they average to is S or a smear.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import styling as st, demos as dm
from helpers.cache import cached

# ---------------------------------------------------------------- layout constants
WIN_H = 1.30            # displayed height of one 96-px particle window
WIN_GAP = 1.47          # centre-to-centre spacing in the 3 x 4 grid
GRID_CENTER = np.array([-3.75, -0.62, 0.0])
SHIFT_EXAGGERATION = 3  # shift arrows are drawn 3x longer than the true offset
RAW_RANGE = (-35.0, 45.0)    # fixed grey scale for the raw windows (sigma = 13.7)
MEAN_RANGE = (-8.0, 34.0)    # fixed grey scale for every mean image (S max = 30.8)
UN_GAIN = 3                  # the unaligned mean (max ~10) is shown at 3x display contrast


def _circle_alpha(n: int, frac: float = 0.47) -> np.ndarray:
    """Soft-edged circular field-of-view mask (x, y)-indexed."""
    yy, xx = np.indices((n, n)) - n // 2
    r = np.sqrt(xx ** 2 + yy ** 2)
    return np.clip(frac * n - r + 0.5, 0, 1)


def _wrap(deg: float) -> float:
    return (deg + 180.0) % 360.0 - 180.0


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a.ravel(), b.ravel())[0, 1])


class Scene02AlignmentAveraging(Scene):
    def construct(self):
        d = cached("alignment", dm.alignment)
        S = d["signal"]
        n = S.shape[0]
        voxel = float(d["voxel"])
        ns = d["ns"].astype(int)
        snr_meas = d["snr_measured"]
        snr1 = float(snr_meas[0])
        alpha = _circle_alpha(n)
        a_per_unit = WIN_H / (n * voxel)          # scene units per angstrom

        def win_img(arr, height=WIN_H, rng=RAW_RANGE):
            return st.array_to_image_mobject(arr, height=height, vmin=rng[0], vmax=rng[1],
                                             nearest=True, alpha=alpha)

        def ring(center, height, color, width=2.5):
            return Circle(radius=0.47 * height + 0.02, color=color,
                          stroke_width=width).move_to(center)

        # ------------------------------------------------------------ opening
        title = st.title_card("2D alignment and averaging",
                              "Why does averaging help only after rotations and shifts are aligned?")
        badge = st.simulation_badge()
        self.add(badge)
        self.play(FadeIn(title, shift=0.15 * DOWN), run_time=1.0)

        # ------------------------------------------------------------ beat 1 (0-7 s)
        centers = []
        for r in range(3):
            for c in range(4):
                centers.append(GRID_CENTER + np.array([(c - 1.5) * WIN_GAP,
                                                       (1 - r) * WIN_GAP, 0]))
        raw = [win_img(d["images"][i]).move_to(centers[i]) for i in range(12)]
        rings = VGroup(*[ring(centers[i], WIN_H, st.OBS) for i in range(12)])
        grid_label = st.T("12 particle windows, one view", 28, st.OBS)
        grid_label.next_to(rings, UP, buff=0.18)

        model = st.M(r"I_i", r"=", r"T_{\mathbf t_i}R_{\phi_i}", r"S", r"+", r"\epsilon_i",
                     size=44)
        model[0].set_color(st.OBS)
        model[3].set_color(st.TRUTH)
        model[5].set_color(st.RESID)
        model.move_to([3.4, 1.6, 0])

        S_img = st.array_to_image_mobject(S, height=1.75, vmin=MEAN_RANGE[0],
                                          vmax=MEAN_RANGE[1], alpha=alpha)
        S_img.move_to([1.45, -0.3, 0])
        S_ring = ring(S_img.get_center(), 1.75, st.TRUTH)
        S_lab = VGroup(st.T("S: the same 2D view", 28, st.TRUTH),
                       st.T("in every window", 28, st.TRUTH)).arrange(DOWN, aligned_edge=LEFT,
                                                                       buff=0.1)
        S_lab.next_to(S_ring, RIGHT, buff=0.3).shift(0.35 * UP)
        noise_lab = st.M(r"\epsilon_i:\ \text{i.i.d. noise}", size=36, color=st.RESID)
        snr_lab = st.M(rf"\mathrm{{SNR}}={snr1:.2f}\ \text{{per image}}", size=36,
                       color=st.RESID)
        VGroup(noise_lab, snr_lab).arrange(DOWN, aligned_edge=LEFT, buff=0.15).next_to(
            S_lab, DOWN, aligned_edge=LEFT, buff=0.35)
        class_note = VGroup(st.T("In practice, 2D classification", 28, st.MUTED),
                            st.T("first groups particles by view.", 28, st.MUTED)).arrange(
            DOWN, aligned_edge=LEFT, buff=0.1)
        class_note.move_to([0, -2.2, 0]).align_to(S_ring, LEFT)
        disclaimer = st.phantom_disclaimer()
        disclaimer.to_edge(DOWN, buff=st.MARGIN)

        self.play(LaggedStart(*[FadeIn(Group(raw[i], rings[i])) for i in range(12)],
                              lag_ratio=0.08), FadeIn(grid_label), FadeIn(disclaimer),
                  run_time=1.3)
        self.play(Write(model), run_time=0.9)
        self.play(FadeIn(S_img), Create(S_ring), FadeIn(S_lab), run_time=1.0)
        self.play(FadeIn(noise_lab), FadeIn(snr_lab), FadeIn(class_note), run_time=0.7)
        self.wait(2.1)                                                     # t ~ 7.0

        # ------------------------------------------------------------ beat 2 (7-13 s)
        # faint misregistration arrows: true offsets t_i (unknown to the analyst)
        def shift_arrow(i, t_vec, color, opacity):
            tip = centers[i] + SHIFT_EXAGGERATION * a_per_unit * np.array([t_vec[0], t_vec[1], 0])
            a = Arrow(centers[i], tip, buff=0, color=color, stroke_width=4,
                      max_tip_length_to_length_ratio=0.35, stroke_opacity=opacity,
                      fill_opacity=opacity)
            return a.set_stroke(st.BG, width=9, background=True)   # dark halo over noise

        mis = VGroup(*[shift_arrow(i, d["t"][i], st.TEXT, 0.9) for i in range(12)])
        self.play(FadeOut(Group(S_img, S_ring, S_lab, noise_lab, snr_lab, class_note,
                                disclaimer)),
                  run_time=0.6)

        un_img = st.array_to_image_mobject(d["unaligned_mean"][4], height=2.3,
                                           vmin=MEAN_RANGE[0] / UN_GAIN,
                                           vmax=MEAN_RANGE[1] / UN_GAIN, alpha=alpha)
        un_img.move_to([3.55, -0.9, 0])
        un_ring = ring(un_img.get_center(), 2.3, st.EST)
        un_lab = st.T("Unaligned mean, N = 256", 30, st.EST)
        un_lab.next_to(un_ring, UP, buff=0.2)
        un_gain = st.T(f"(shown at {UN_GAIN}× contrast)", st.SMALL_SIZE, st.MUTED)
        un_gain.next_to(un_ring, DOWN, buff=0.15)
        un_eq = st.M(r"\overline I_N^{\,\mathrm{unal}}=\frac1N\sum_i I_i\ \longrightarrow\ "
                     r"\mathbb E_{\phi,\mathbf t}\!\left[T_{\mathbf t}R_{\phi}S\right]\neq S",
                     size=38)
        un_eq.move_to([3.4, 1.45, 0])
        mis_lab = VGroup(st.T("misregistered:", 28, st.MUTED),
                         st.M(r"\phi_i,\ \mathbf t_i\ \text{all differ}", size=36,
                              color=st.MUTED)).arrange(RIGHT, buff=0.2)
        mis_lab.next_to(rings, DOWN, buff=0.3)
        pipe = Arrow(rings.get_right() + 0.05 * RIGHT, un_ring.get_left() + 0.05 * LEFT,
                     color=st.MUTED, stroke_width=4, buff=0.1)

        self.play(FadeOut(model, shift=0.2 * UP), FadeIn(un_eq, shift=0.2 * UP), run_time=1.0)
        self.play(LaggedStart(*[GrowArrow(a) for a in mis], lag_ratio=0.05),
                  FadeIn(mis_lab), run_time=1.0)
        self.play(GrowArrow(pipe), FadeIn(un_img), Create(un_ring), FadeIn(un_lab),
                  FadeIn(un_gain), run_time=1.2)
        self.wait(1.6)                                                     # t ~ 13.0

        # ------------------------------------------------------------ beat 3 (13-22 s)
        self.play(FadeOut(Group(un_img, un_ring, un_lab, un_gain, pipe, mis_lab)),
                  run_time=0.6)
        al_eq = st.M(r"I_i^{\mathrm{aligned}}", r"=", r"R_{-\hat\phi_i}\,T_{-\hat{\mathbf t}_i}",
                     r"I_i", size=44)
        al_eq[0].set_color(st.OBS)
        al_eq[2].set_color(st.EST)
        al_eq[3].set_color(st.OBS)
        al_eq.move_to([3.4, 1.45, 0])
        est_lab = VGroup(st.T("estimated pose:", 30, st.EST),
                         st.M(r"\hat\phi_i\ \text{(angle)},\ \ \hat{\mathbf t}_i\ \text{(shift)}",
                              size=36, color=st.EST)).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
        est_lab.move_to([3.4, 0.05, 0])
        assume = VGroup(st.T("assumed accurate:", 28, st.MUTED),
                        st.T("errors ≈ 1.5°, 0.4 Å", 28, st.MUTED)).arrange(
            DOWN, aligned_edge=LEFT, buff=0.12)
        assume.next_to(est_lab, DOWN, aligned_edge=LEFT, buff=0.45)
        foot = st.T(f"shift arrows drawn {SHIFT_EXAGGERATION}× longer", st.SMALL_SIZE, st.MUTED)
        foot.next_to(rings, DOWN, buff=0.2)

        # orange estimates: angle tick + arc, shift arrow
        ticks, arcs, est_arrows, refs = VGroup(), VGroup(), VGroup(), VGroup()
        rr = 0.47 * WIN_H
        for i in range(12):
            ang = np.radians(_wrap(d["phi_est"][i]))
            ticks.add(Line(centers[i], centers[i] + rr * np.array([np.cos(ang), np.sin(ang), 0]),
                           color=st.EST, stroke_width=4).set_stroke(st.BG, width=8,
                                                                    background=True))
            arcs.add(Arc(radius=0.28, start_angle=0, angle=ang, arc_center=centers[i],
                         color=st.EST, stroke_width=3))
            est_arrows.add(shift_arrow(i, d["t_est"][i], st.EST, 1.0))
            refs.add(Line(centers[i] + (rr + 0.02) * RIGHT, centers[i] + (rr + 0.14) * RIGHT,
                          color=st.MUTED, stroke_width=3))

        # ticks + 0° marks stay above the aligned windows: they are the evidence of alignment
        ticks.set_z_index(2)
        refs.set_z_index(2)
        self.play(FadeOut(un_eq, shift=0.2 * UP), FadeIn(al_eq, shift=0.2 * UP), run_time=1.0)
        self.play(ReplacementTransform(mis, est_arrows), FadeIn(refs), Create(ticks),
                  Create(arcs), FadeIn(est_lab), FadeIn(foot), run_time=1.3)
        self.play(FadeIn(assume), run_time=0.6)
        self.wait(0.4)
        # undo the shift (move each window by -t_hat) ...
        self.play(*[raw[i].animate.shift(-a_per_unit * np.array([*d["t_est"][i], 0]))
                    for i in range(12)],
                  FadeOut(est_arrows), FadeOut(foot), run_time=0.9)
        # ... then undo the rotation about the box centre
        self.play(*[Rotate(raw[i], angle=-np.radians(_wrap(d["phi_est"][i])),
                           about_point=centers[i]) for i in range(12)],
                  *[Rotate(ticks[i], angle=-np.radians(_wrap(d["phi_est"][i])),
                           about_point=centers[i]) for i in range(12)],
                  FadeOut(arcs), run_time=2.0)
        aligned = [win_img(d["aligned_examples"][i]).move_to(centers[i]) for i in range(12)]
        common = st.T("one common reference frame", 28, st.EST)
        common.next_to(rings, UP, buff=0.18)
        self.play(*[FadeOut(raw[i]) for i in range(12)], *[FadeIn(aligned[i]) for i in range(12)],
                  FadeOut(grid_label), FadeIn(common), run_time=0.8)
        self.wait(0.9)                                                     # t ~ 22.0

        # ------------------------------------------------------------ beat 4 (22-29 s)
        self.play(FadeOut(Group(*aligned, rings, ticks, refs, common, est_lab, assume,
                                al_eq)), run_time=0.6)

        mean_c = np.array([-4.45, -0.55, 0])
        mean_h = 2.6
        means = [st.array_to_image_mobject(d["aligned_mean"][k], height=mean_h,
                                           vmin=MEAN_RANGE[0], vmax=MEAN_RANGE[1],
                                           alpha=alpha, nearest=True).move_to(mean_c)
                 for k in range(5)]
        mean_ring = ring(mean_c, mean_h, st.EST)
        mean_lab = st.M(r"\overline I_N=\frac1N\sum_{i=1}^{N}I_i^{\mathrm{aligned}}", size=38)
        mean_lab.next_to(mean_ring, UP, buff=0.2)

        def n_readout(k):
            g = VGroup(st.M(rf"N={ns[k]}", size=38, color=st.EST),
                       st.M(rf"\mathrm{{SNR}}={snr_meas[k]:.2f}", size=38, color=st.OBS))
            return g.arrange(RIGHT, buff=0.6).next_to(mean_ring, DOWN, buff=0.3)

        plot = self._snr_plot(d, center=(3.45, -0.45), width=5.3, height=2.5)
        law = st.M(r"\sigma_{\mathrm{mean}}=\frac{\sigma}{\sqrt N}", r",\qquad",
                   r"\mathrm{SNR}_{\mathrm{amp}}\propto\sqrt N", size=38)
        law.move_to([3.2, 1.85, 0])
        assumptions = VGroup(
            st.T("assumes i.i.d. zero-mean Gaussian noise (an approximation),",
                 st.SMALL_SIZE, st.MUTED),
            st.T("identical S, accurate alignment · amplitude SNR (power SNR ∝ N)",
                 st.SMALL_SIZE, st.MUTED)
        ).arrange(DOWN, buff=0.1)
        assumptions.to_edge(DOWN, buff=st.MARGIN)

        readout = n_readout(0)
        self.play(FadeIn(means[0]), Create(mean_ring), FadeIn(mean_lab), FadeIn(readout),
                  Create(plot["axes"]), FadeIn(plot["labels"]), FadeIn(law),
                  FadeIn(assumptions), run_time=1.1)
        self.play(Create(plot["curve"]), FadeIn(plot["legend"]),
                  FadeIn(plot["dots"][0], scale=0.5), run_time=0.9)
        for k in range(1, 5):
            new_read = n_readout(k)
            # new mean underneath, old one fades off the top: no dimmed midpoint
            self.add(means[k])
            self.bring_to_front(means[k - 1])
            self.play(FadeOut(means[k - 1]),
                      Succession(FadeOut(readout, shift=0.15 * UP),
                                 FadeIn(new_read, shift=0.15 * UP)),
                      FadeIn(plot["dots"][k], scale=0.5), run_time=0.9)
            readout = new_read
            self.wait(0.25)
        self.wait(0.6)                                                     # t ~ 29.0

        # ------------------------------------------------------------ beat 5 (29-35 s)
        fh = 2.5
        uc, ac = np.array([-5.25, 0.15, 0]), np.array([-2.35, 0.15, 0])
        un_f = st.array_to_image_mobject(d["unaligned_mean"][4], height=fh,
                                         vmin=MEAN_RANGE[0] / UN_GAIN,
                                         vmax=MEAN_RANGE[1] / UN_GAIN, alpha=alpha).move_to(uc)
        al_f = st.array_to_image_mobject(d["aligned_mean"][4], height=fh,
                                         vmin=MEAN_RANGE[0], vmax=MEAN_RANGE[1],
                                         alpha=alpha, nearest=True).move_to(ac)
        un_r, al_r = ring(uc, fh, st.EST), ring(ac, fh, st.EST)
        un_t = st.T("unaligned", 32, st.TEXT).next_to(un_r, UP, buff=0.18)
        al_t = st.T("aligned", 32, st.TEXT).next_to(al_r, UP, buff=0.18)
        mid_x = (uc[0] + ac[0]) / 2
        cc_t = st.T("correlation with the true S", 28, st.MUTED)
        cc_t.move_to([mid_x, un_r.get_bottom()[1] - 0.38, 0])
        cu = _corr(d["unaligned_mean"][4], S)
        ca = _corr(d["aligned_mean"][4], S)
        cc_u = st.M(rf"{cu:.2f}", size=44, color=st.TEXT).move_to([uc[0], cc_t.get_y() - 0.55, 0])
        cc_a = st.M(rf"{ca:.2f}", size=44, color=st.TEXT).move_to([ac[0], cc_t.get_y() - 0.55, 0])
        foot5 = VGroup(st.T(f"same {ns[-1]} particles in each mean;", st.SMALL_SIZE, st.MUTED),
                       st.T(f"unaligned shown at {UN_GAIN}× contrast", st.SMALL_SIZE,
                            st.MUTED)).arrange(DOWN, buff=0.1)
        foot5.move_to([mid_x, 0, 0]).to_edge(DOWN, buff=st.MARGIN)
        foot5.set_x(max(mid_x, -config.frame_width / 2 + st.MARGIN + foot5.width / 2))

        card = st.takeaway_card(
            r"\sigma_{\mathrm{mean}}=\frac{\sigma}{\sqrt N},\quad "
            r"\mathrm{SNR}_{\mathrm{amp}}\propto\sqrt N",
            "Noise falls as 1/√N either way;\nonly copies in register add up\nto S: align first, then average.",
            "\nlike raising redundancy, but images\nmust first be grouped by view.")
        card.to_edge(RIGHT, buff=st.MARGIN).set_y(-0.45)

        self.play(FadeOut(Group(means[4], mean_ring, mean_lab, readout, law, assumptions,
                                plot)), run_time=0.7)
        self.play(FadeIn(un_f), Create(un_r), FadeIn(un_t), FadeIn(al_f), Create(al_r),
                  FadeIn(al_t), run_time=0.9)
        self.play(FadeIn(cc_t), FadeIn(cc_u), FadeIn(cc_a), FadeIn(foot5),
                  FadeIn(card, shift=0.1 * UP), run_time=0.8)
        self.wait(4.0)                                                     # freeze

    # -------------------------------------------------------------------- plot
    def _snr_plot(self, d, center, width, height, tick=28):
        """Measured SNR(N) of the aligned means (cyan dots) against SNR_1 * sqrt(N)."""
        ns = d["ns"].astype(int)
        snr = d["snr_measured"]
        snr1 = float(snr[0])
        ax = Axes(x_range=[0, 8.4, 2], y_range=[0, 5.6, 1], x_length=width, y_length=height,
                  tips=False, axis_config=dict(color=st.MUTED, stroke_width=2,
                                               include_ticks=True, tick_size=0.06))
        ax.move_to([center[0], center[1], 0])
        labels = VGroup()
        xt = VGroup(*[st.T(str(N), tick, st.MUTED).next_to(ax.c2p(np.log2(N), 0), DOWN,
                                                          buff=0.14) for N in ns])
        yt = VGroup(*[st.T(str(v), tick, st.MUTED).next_to(ax.c2p(0, v), LEFT, buff=0.14)
                      for v in (0, 2, 4)])
        xl = st.T("N  (log scale)", tick, st.MUTED).next_to(xt, DOWN, buff=0.12)
        yl = st.M(r"\mathrm{SNR}=\mathrm{std}(S)\,/\,\sigma_{\mathrm{mean}}", size=36,
                  color=st.TEXT)
        yl.next_to(ax.c2p(0, 5.6), UP, buff=0.12).align_to(yt, LEFT)
        labels.add(xt, yt, xl, yl)
        curve = DashedVMobject(ax.plot(lambda x: snr1 * 2 ** (x / 2), x_range=[0, 8.2],
                                       color=st.TEXT, stroke_width=3), num_dashes=40)
        dots = VGroup(*[Dot(ax.c2p(np.log2(N), s), radius=0.085, color=st.OBS)
                        for N, s in zip(ns, snr)])
        leg = VGroup(VGroup(Dot(radius=0.085, color=st.OBS),
                            st.T("measured", tick, st.OBS)).arrange(
                                RIGHT, buff=0.2),
                     VGroup(DashedLine(ORIGIN, 0.45 * RIGHT, color=st.TEXT, stroke_width=3,
                                       dash_length=0.08),
                            st.M(rf"{snr1:.2f}\sqrt N\ \ \text{{prediction}}", size=36,
                                 color=st.TEXT)).arrange(RIGHT, buff=0.2)
                     ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        leg.next_to(ax.c2p(0.35, 5.5), DOWN + RIGHT, buff=0)
        return VDict({"axes": ax, "labels": labels, "curve": curve, "dots": dots,
                      "legend": leg})
