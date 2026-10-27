"""Clip 01 — Forward image formation: signal, CTF and noise (~33 s).

Storyboard beats (spec §3, Scene 01) -> code sections:
    beat 1  0-6 s    phantom + beam rays -> clean projection in the detector window
    beat 2  6-16 s   split screen: P_R rho | CTF(k) plot | CTF-modulated image;
                     three defocus variants, then one frequency band at a time:
                     CTF ~ -1 reverses contrast, CTF = 0 suppresses, CTF ~ +1 transfers
    beat 3  16-24 s  noise slider (sigma) from clear to low SNR, log-scale
                     signal-RMS / noise-sigma bars and a live SNR readout
    beat 4  24-33 s  Fourier-space forward model with highlighted terms
                     (structure / orientation / optics / shift / noise),
                     "inverse problem" label, takeaway freeze (~4 s)

Every image and curve is computed: images and CTF curves come from
dm.image_formation (cached); band-pass components from _bands() (cached
below); the moving CTF readout and SNR numbers from optics.ctf_1d and the
cached sigma / signal-RMS values.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *

from helpers import demos as dm
from helpers import optics as op
from helpers import styling as st
from helpers.cache import cached

# Optics (CTF) colour: not one of the reserved semantic colours; same blue family
# as the CTF-sign end of st.DIVERGING.
OPTICS = "#86AEF7"
WORK_DF = 1.0          # working defocus (um); data["noisy"] is built on this CTF image
SHOW_SIGNAL = True     # the observed panel is signal + noise; each term can be switched off
SHOW_NOISE = True


# ------------------------------------------------------------- clip-local data
def _bands(ks=(0.05, 0.07, 0.087), width=0.0045, defocus_um=WORK_DF) -> dict:
    """Narrow Gaussian band-pass components of the clean projection around each
    k0 (1/A), before and after multiplication by the CTF (optics.apply_ctf)."""
    d = cached("image_formation", dm.image_formation)
    clean, voxel = d["clean"], float(d["voxel"])
    n = clean.shape[0]
    K = np.hypot(*op.freq_grid(n, voxel))
    F = np.fft.fft2(clean)
    before, after = [], []
    for k0 in ks:
        b = np.real(np.fft.ifft2(F * np.exp(-0.5 * ((K - k0) / width) ** 2)))
        before.append(b)
        after.append(op.apply_ctf(b, voxel, defocus_um=defocus_um))
    return dict(ks=np.array(ks), before=np.array(before), after=np.array(after),
                ctf_at=op.ctf_1d(np.array(ks), defocus_um))


# ------------------------------------------------------------- small builders
def _img(arr, h, nearest=False, pct=None, **kw):
    """Array -> ImageMobject; pct = (lo, hi) percentile auto-contrast."""
    if pct is not None:
        lo, hi = np.percentile(arr, pct)
        kw.update(vmin=lo, vmax=hi)
    return st.array_to_image_mobject(arr, height=h, nearest=nearest, **kw)


def _frame(img, color, width=3):
    return Rectangle(width=img.width, height=img.height, stroke_color=color,
                     stroke_width=width).move_to(img)


def _header(text):
    h = st.T(text, 36, st.TEXT, weight="BOLD")
    h.to_corner(UL, buff=st.MARGIN)
    return h


def _tlab(text, tex, color):
    """Text label followed by a small math symbol, baseline-aligned."""
    g = VGroup(st.T(text, 28, color), st.M(tex, size=38, color=color)).arrange(RIGHT, buff=0.12)
    g[1].align_to(g[0], DOWN).shift(0.04 * DOWN)
    return g


def _cap_line(pre, tex, post, color):
    """Body-text line with one inline math symbol, e.g. '(small s_i)'."""
    g = VGroup(st.T(pre, 28, st.TEXT), st.M(tex, size=38, color=color), st.T(post, 28, st.TEXT))
    g.arrange(RIGHT, buff=0.1)
    g[2].shift(0.05 * LEFT)
    base = g[0][-1].get_bottom()[1]          # last letter of `pre` (l / e) sits on the baseline
    g[1].shift((base - g[1][0][0].get_bottom()[1]) * UP)   # main glyph on the same baseline
    g[2].align_to(g[0], DOWN)                # ')' matches '(' (both descend)
    return g


def _df_label(df):
    return st.T(f"Δf = {df:.1f} µm", 28, OPTICS)


class Scene01ImageFormation(Scene):
    def construct(self):
        d = cached("image_formation", dm.image_formation)
        k = d["k"]
        curves = d["ctf_curves"]
        defoci = [float(x) for x in d["defoci"]]
        j_work = defoci.index(WORK_DF)
        signal = d["ctf_images"][j_work]
        noise_part = d["noisy"] - signal[None]          # sigma_i * (one fixed white pattern)
        observed = (signal[None] if SHOW_SIGNAL else 0) + (noise_part if SHOW_NOISE else 0)

        # band frequencies from the working-defocus CTF: first trough (CTF ~ -1),
        # first zero, second peak (CTF ~ +1)
        c1 = curves[j_work]
        kz = float(d["first_zero"][j_work])
        k_rev = float(k[np.argmin(np.where(k < kz, c1, 9))])
        k_keep = float(k[np.argmax(np.where((k > kz) & (k < 0.1), c1, -9))])
        ks = tuple(round(x, 4) for x in (k_rev, kz, k_keep))
        bands = cached("s01_image_formation_bands", _bands, ks=ks, width=0.0045,
                       defocus_um=WORK_DF)

        badge = st.simulation_badge().set_z_index(20)
        self.add(badge)

        # ================================================ beat 1 (0-6 s): projection
        title = st.title_card("Forward image formation",
                              "Why are particle images noisy, and why must we model the microscope?")
        self.play(FadeIn(title, shift=0.2 * DOWN), run_time=0.8)

        # rendered at the same orientation R as the detector projection (dm.VIEW),
        # so the silhouette on the detector is recognisably the same object
        phantom = st.density_image(dm.cached_render(**dm.VIEW), st.TRUTH, height=2.9)
        phantom.move_to([-2.9, -0.25, 0]).set_z_index(2)
        disclaimer = st.phantom_disclaimer()
        disclaimer.move_to([-6.6 + disclaimer.width / 2, -2.0, 0])

        det_img = _img(d["clean"], 2.8).move_to([4.45, -0.25, 0])
        det_frame = _frame(det_img, st.TRUTH)
        det_lab = st.T("detector: one particle window", 28, st.TEXT)
        det_lab.next_to(det_frame, UP, buff=0.15).align_to(det_frame, RIGHT)

        x0, x1 = -6.55, det_img.get_left()[0] - 0.04
        rays = VGroup(*[Arrow([x0, y, 0], [x1, y, 0], buff=0, stroke_width=3,
                              color=st.MUTED, max_tip_length_to_length_ratio=0.03,
                              max_stroke_width_to_length_ratio=2)
                        for y in np.linspace(-1.3, 0.8, 5)])
        rays.set_opacity(0.75)
        beam_lab = st.T("electron beam (along z)", 28, st.MUTED)
        beam_lab.next_to(rays, UP, buff=0.15).align_to(rays, LEFT)
        proj_eq = st.M(r"P_R\rho(x,y)=\int \rho(R^{\top}\mathbf r)\,dz", size=38)
        proj_eq.next_to(det_frame, DOWN, buff=0.45).set_x(det_img.get_x())

        self.play(FadeIn(phantom), FadeIn(disclaimer), run_time=0.8)
        self.play(LaggedStart(*[GrowArrow(r) for r in rays], lag_ratio=0.08),
                  FadeIn(beam_lab), run_time=1.0)
        self.play(FadeIn(det_img), Create(det_frame), FadeIn(det_lab), Write(proj_eq),
                  run_time=1.0)
        head1 = st.T("The image is a projection, not the structure", 36, st.TEXT, weight="BOLD")
        head1.to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(head1, shift=0.15 * UP), run_time=0.6)
        self.wait(1.4)

        # ================================================ beat 2 (6-15 s): CTF
        hdr2 = _header("Optics: the CTF multiplies each frequency")
        eq2 = st.M(r"\widetilde I(\mathbf k)=", r"\mathrm{CTF}(k)", r"\,\widehat{P_R\rho}(\mathbf k)",
                   r",\qquad \mathrm{CTF}(k)=-\sin\chi(k)", size=40)
        eq2[1].set_color(OPTICS)
        eq2[2].set_color(st.TRUTH)
        eq2.move_to([0, 2.5, 0])

        img_y = -0.35
        left_pos = np.array([-4.85, img_y, 0])
        right_pos = np.array([4.85, img_y, 0])
        ph = 2.6

        self.play(*[FadeOut(m) for m in (title, phantom, disclaimer, rays, beam_lab, head1,
                                         proj_eq, det_lab)], run_time=0.6)
        self.play(Group(det_img, det_frame).animate.scale(ph / 2.8).move_to(left_pos),
                  FadeIn(hdr2), Write(eq2), run_time=0.9)
        clean_img, clean_frame = det_img, det_frame
        clean_lab = st.T("clean projection", 28, st.TRUTH).next_to(clean_frame, UP, buff=0.15)

        # CTF plot (centre)
        ax = Axes(x_range=[0, 0.2, 0.05], y_range=[-1, 1, 1], x_length=4.2, y_length=1.9,
                  tips=False, axis_config=dict(color=st.MUTED, stroke_width=2))
        ax.move_to([0.15, 0.95, 0])
        box = Rectangle(width=4.2, height=1.9, stroke_color=st.MUTED, stroke_width=1.2,
                        stroke_opacity=0.5).move_to(ax.c2p(0.1, 0))
        ylabs = VGroup(*[st.T(s, 26, st.MUTED).next_to(ax.c2p(0, v), LEFT, buff=0.12)
                         for s, v in (("+1", 1), ("0", 0), ("−1", -1))])
        ctf_lab = st.T("CTF", 28, OPTICS).rotate(PI / 2).next_to(ylabs, LEFT, buff=0.12)
        xlabs = VGroup(*[st.T(s, 26, st.MUTED).next_to(ax.c2p(kv, -1), DOWN, buff=0.1)
                         for s, kv in (("20 Å", 0.05), ("10 Å", 0.1), ("6.7 Å", 0.15),
                                       ("5 Å", 0.2))])
        xcap = st.T("spatial frequency k  (labels: d = 1/k)", 26, st.MUTED)
        xcap.next_to(xlabs, DOWN, buff=0.08).set_x(box.get_x())
        sel = k <= 0.2

        def ctf_curve(j):
            pts = [ax.c2p(a, b) for a, b in zip(k[sel], curves[j][sel])]
            return VMobject(stroke_color=OPTICS, stroke_width=3.5).set_points_as_corners(pts)

        plot = VGroup(box, ax, ylabs, ctf_lab, xlabs, xcap)
        curve = ctf_curve(len(defoci) - 1)

        ctf_imgs = [_img(d["ctf_images"][j], ph).move_to(right_pos) for j in range(len(defoci))]
        right_frame = _frame(ctf_imgs[0], OPTICS).set_z_index(10)
        right_lab = st.T("after the optics", 28, OPTICS)
        df_lab = _df_label(defoci[-1])
        VGroup(right_lab, df_lab).arrange(DOWN, buff=0.08).next_to(right_frame, UP, buff=0.15)
        foot2 = st.T("weak-phase model: 300 kV,\nCs 2.7 mm, 7% amp. contrast", 26, st.MUTED,
                     line_spacing=0.8)
        foot2.next_to(clean_frame, DOWN, buff=0.25).align_to(clean_frame, LEFT)

        cur_img = ctf_imgs[-1]
        self.play(FadeIn(clean_lab), FadeIn(plot), Create(curve), FadeIn(cur_img),
                  Create(right_frame), FadeIn(right_lab), FadeIn(df_lab), FadeIn(foot2),
                  run_time=1.0)
        # defocus variants: the zeros move and the image changes non-monotonically
        for j in range(len(defoci) - 2, j_work - 1, -1):
            new_img, new_df = ctf_imgs[j], _df_label(defoci[j]).move_to(df_lab)
            new_img.set_z_index(cur_img.z_index + 1)    # fade in over the old image: no dip
            self.play(Transform(curve, ctf_curve(j)), FadeTransform(df_lab, new_df),
                      FadeIn(new_img), run_time=0.65)
            self.remove(cur_img)
            cur_img, df_lab = new_img, new_df
        not_blur = st.T("oscillates, with zeros and sign flips: not a blur", 30, st.TEXT,
                        weight="BOLD")
        not_blur.move_to([0.15, -3.05, 0])        # below the weak-phase footnote
        self.play(FadeIn(not_blur), run_time=0.4)
        self.wait(0.6)

        # one frequency band at a time
        kt = ValueTracker(ks[0])
        ctf_w = lambda x: float(op.ctf_1d(np.array([x]), WORK_DF)[0])
        dot = always_redraw(lambda: Dot(ax.c2p(kt.get_value(), ctf_w(kt.get_value())),
                                        radius=0.08, color=st.TEXT))
        vline = always_redraw(lambda: DashedLine(ax.c2p(kt.get_value(), -1),
                                                 ax.c2p(kt.get_value(), 1),
                                                 stroke_width=2, color=st.TEXT,
                                                 dash_length=0.06))
        bh = 1.6
        by = -2.05
        bpos_l, bpos_r = np.array([-1.3, by, 0]), np.array([1.6, by, 0])
        crop = slice(24, 104)            # zoom on the particle (central 80 of 128 px)
        vmax = [np.abs(b).max() for b in bands["before"]]
        b_before = [_img(b[crop, crop], bh, symmetric=True, vmax=m).move_to(bpos_l)
                    for b, m in zip(bands["before"], vmax)]
        b_after = [_img(b[crop, crop], bh, symmetric=True, vmax=m).move_to(bpos_r)
                   for b, m in zip(bands["after"], vmax)]
        bfr = VGroup(_frame(b_before[0], st.TRUTH, 2.5), _frame(b_after[0], OPTICS, 2.5))
        bfr.set_z_index(10)
        barrow = Arrow(bpos_l + (bh / 2 + 0.06) * RIGHT, bpos_r + (bh / 2 + 0.06) * LEFT, buff=0,
                       color=OPTICS, stroke_width=4, max_tip_length_to_length_ratio=0.25)
        times = st.M(r"\times\mathrm{CTF}", size=36, color=OPTICS).next_to(barrow, UP, buff=0.1)
        band_cap = st.T("one frequency band:", 28, st.MUTED)
        band_cap.next_to(bfr[0], LEFT, buff=0.2)
        status_txt = ["CTF ≈ −1: contrast reversed", "CTF = 0 at this k: band suppressed",
                      "CTF ≈ +1: transferred"]
        status = [st.T(s, 30, st.TEXT, weight="BOLD").next_to(bfr, DOWN, buff=0.12)
                  for s in status_txt]

        self.play(FadeOut(not_blur), FadeOut(foot2), run_time=0.3)
        self.play(FadeIn(dot), FadeIn(vline), FadeIn(bfr), FadeIn(b_before[0]),
                  FadeIn(b_after[0]), GrowArrow(barrow), FadeIn(times), FadeIn(band_cap),
                  FadeIn(status[0]), run_time=0.7)
        self.wait(0.9)
        for i in (1, 2):
            b_before[i].set_z_index(1)
            b_after[i].set_z_index(1)
            self.play(kt.animate.set_value(ks[i]), FadeIn(b_before[i]), FadeIn(b_after[i]),
                      FadeOut(status[i - 1], run_time=0.3), run_time=0.8)
            self.remove(b_before[i - 1], b_after[i - 1])
            self.play(FadeIn(status[i]), run_time=0.3)
            self.wait(1.0)

        # ================================================ beat 3 (15-24 s): noise
        hdr3 = _header("Noise: an independent error in every pixel")
        eq3 = st.M(r"I_i(\mathbf x)=", r"s_i(\mathbf x)", "+", r"\epsilon_i(\mathbf x)",
                   r",\qquad \epsilon_i\sim\mathcal N(0,\sigma^2)\ \text{i.i.d.}", size=40)
        eq3[0].set_color(st.OBS)
        eq3[1].set_color(st.TRUTH)
        eq3[3].set_color(st.RESID)
        eq3.move_to([0, 2.5, 0])

        self.play(*[FadeOut(m) for m in (dot, vline, bfr, b_before[2], b_after[2], barrow, times,
                                         band_cap, status[2], plot, curve, df_lab, clean_img,
                                         clean_frame, clean_lab, right_lab, hdr2, eq2)],
                  run_time=0.6)
        self.remove(dot, vline)

        sig_pos = np.array([-4.85, -0.25, 0])
        obs_pos = np.array([-1.3, -0.25, 0])
        oh = 2.6
        sig_frame_new = _frame(cur_img, st.TRUTH).move_to(sig_pos).set_z_index(10)
        sig_lab = _tlab("signal", r"s_i", st.TRUTH).next_to(sig_frame_new, UP, buff=0.15)
        obs_imgs = [_img(x, oh, nearest=True, pct=(0.5, 99.5)).move_to(obs_pos) for x in observed]
        obs_frame = _frame(obs_imgs[0], st.OBS).set_z_index(10)
        obs_lab = _tlab("observed window", r"I_i", st.OBS).next_to(obs_frame, UP, buff=0.15)

        self.play(cur_img.animate.move_to(sig_pos), Transform(right_frame, sig_frame_new),
                  FadeIn(hdr3), Write(eq3), run_time=0.9)

        # right column: slider, SNR readout, log-scale amplitude bars
        sig = d["noise_sigma"]
        srms = float(d["signal_rms"])
        ls = np.log10(sig)
        lt = ValueTracker(ls[0])
        ty = 1.15
        track = Line([2.35, ty, 0], [5.75, ty, 0], color=st.MUTED, stroke_width=5)
        knob_x = lambda: track.get_left()[0] + (lt.get_value() - ls[0]) / (ls[-1] - ls[0]) * track.width
        knob = always_redraw(lambda: Circle(radius=0.14, color=st.RESID, fill_opacity=1,
                                            stroke_width=0).move_to([knob_x(), ty, 0]))
        sl_ends = VGroup(st.T("clear", 28, st.MUTED).next_to(track, DOWN, buff=0.14).align_to(track, LEFT),
                         st.T("noisy", 28, st.MUTED).next_to(track, DOWN, buff=0.14).align_to(track, RIGHT))

        def sig_readout():
            g = VGroup(st.T("noise", 28, st.RESID),
                       st.M(r"\sigma=", size=36, color=st.RESID),
                       DecimalNumber(10 ** lt.get_value(), num_decimal_places=2, font_size=36,
                                     color=st.RESID)).arrange(RIGHT, buff=0.12)
            g[0].align_to(g[1], DOWN)
            return g.next_to(track, UP, buff=0.2).align_to(track, LEFT)
        sig_val = always_redraw(sig_readout)

        def snr_readout():
            g = VGroup(st.M(r"\mathrm{SNR}_{\mathrm{amp}}=\frac{\mathrm{RMS}(s)}{\sigma}=", size=36),
                       DecimalNumber(srms / 10 ** lt.get_value(), num_decimal_places=2,
                                     font_size=40)).arrange(RIGHT, buff=0.1)
            return g.move_to([4.0, 0.05, 0])
        snr_val = always_redraw(snr_readout)

        # log axis 0.1 .. 100 (amplitude, image units)
        y0, dec = -2.55, 0.5
        ax_x = 2.85
        ylog = lambda v: y0 + (np.log10(v) + 1) * dec
        bax = Line([ax_x, y0, 0], [ax_x, ylog(100), 0], color=st.MUTED, stroke_width=2)
        bticks = VGroup()
        for v, s in ((0.1, "0.1"), (1, "1"), (10, "10"), (100, "100")):
            bticks.add(Line([ax_x - 0.08, ylog(v), 0], [ax_x, ylog(v), 0], color=st.MUTED,
                            stroke_width=2),
                       st.T(s, 26, st.MUTED).next_to([ax_x - 0.08, ylog(v), 0], LEFT, buff=0.08))
        base = Line([ax_x, y0, 0], [ax_x + 3.4, y0, 0], color=st.MUTED, stroke_width=2)
        bw = 0.7
        sx, nx = ax_x + 0.85, ax_x + 2.75
        s_bar = Rectangle(width=bw, height=ylog(srms) - y0, fill_color=st.TRUTH, fill_opacity=0.9,
                          stroke_width=0).move_to([sx, y0, 0], aligned_edge=DOWN)
        n_bar = always_redraw(lambda: Rectangle(
            width=bw, height=max(ylog(10 ** lt.get_value()) - y0, 1e-3), fill_color=st.RESID,
            fill_opacity=0.9, stroke_width=0).move_to([nx, y0, 0], aligned_edge=DOWN))
        bar_labs = VGroup(st.T("signal RMS", 28, st.TRUTH), st.T("noise σ", 28, st.RESID))
        bar_labs[0].next_to([sx, y0, 0], DOWN, buff=0.1)
        bar_labs[1].next_to([nx, y0, 0], DOWN, buff=0.1)
        bar_labs[1].align_to(bar_labs[0], UP)
        amp_lab = st.T("amplitude (log)", 28, st.MUTED).next_to(bax, UP, buff=0.08).align_to(bax, LEFT)

        self.play(FadeIn(sig_lab), FadeIn(obs_imgs[0]), Create(obs_frame), FadeIn(obs_lab),
                  Create(track), FadeIn(knob), FadeIn(sl_ends), FadeIn(sig_val), FadeIn(snr_val),
                  Create(bax), Create(base), FadeIn(bticks), FadeIn(s_bar), FadeIn(n_bar),
                  FadeIn(bar_labs), FadeIn(amp_lab), run_time=0.9)
        self.wait(0.3)
        for i in range(1, len(sig)):
            obs_imgs[i].set_z_index(i + 1)              # fade in over the old image: no dip
            self.play(lt.animate.set_value(ls[i]), FadeIn(obs_imgs[i]), run_time=0.6,
                      rate_func=smooth)
            self.remove(obs_imgs[i - 1])
        cap3 = VGroup(_cap_line("Low SNR = weak particle contrast (small", r"s_i", ")", st.TRUTH),
                      _cap_line("+ independent per-pixel noise (large", r"\sigma", ")", st.RESID)
                      ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        foot3 = st.T("Gaussian i.i.d. noise is an approximation;\nreal electron-counting noise is dose-dependent",
                     26, st.MUTED, slant="ITALIC", line_spacing=0.8)
        foot3.next_to(cap3, DOWN, buff=0.12).align_to(cap3, LEFT)
        VGroup(cap3, foot3).to_corner(DL, buff=st.MARGIN)
        self.play(FadeIn(cap3), FadeIn(foot3), run_time=0.6)
        self.wait(1.1)

        # ================================================ beat 4 (24-32 s): forward model
        self.play(*[FadeOut(m) for m in (cur_img, right_frame, sig_lab, obs_imgs[-1], obs_frame,
                                         obs_lab, track, knob, sl_ends, sig_val, bax, base,
                                         bticks, s_bar, n_bar, bar_labs, amp_lab, snr_val, cap3,
                                         foot3, hdr3, eq3)], run_time=0.6)
        self.remove(knob, sig_val, n_bar, snr_val)
        hdr4 = _header("The forward model in Fourier space")
        fm = st.M(r"\widetilde I_i(\mathbf k)", "=", r"\mathrm{CTF}_i(\mathbf k)", r"\,",
                  r"\widehat{P_{R_i}\rho}(\mathbf k)", r"\,", r"e^{-2\pi i\,\mathbf k\cdot\mathbf t_i}",
                  "+", r"\widetilde\epsilon_i(\mathbf k)", size=56)
        fm.move_to([0, 1.75, 0])
        proj = fm[4]          # glyphs: 0 hat, 1 P, 2 R, 3 i, 4 rho, 5-7 "(k)"
        R_i, rho = VGroup(proj[2], proj[3]), proj[4]
        self.play(FadeIn(hdr4), Write(fm), run_time=1.1)

        def below(mob, col, txt):
            br = Brace(mob, DOWN, buff=0.1, color=col)
            lab = st.T(txt, 28, col, line_spacing=0.75).next_to(br, DOWN, buff=0.08)
            return VGroup(br, lab)

        def above(mob, col, txt, side):
            top = proj[0].get_top()[1] + 0.08
            lead = Line([mob.get_x(), top, 0], [mob.get_x(), top + 0.3, 0], color=col,
                        stroke_width=3)
            lab = st.T(txt, 28, col).next_to(lead, UP, buff=0.06)
            lab.align_to(lead.get_center() - 0.04 * side, -side)   # extend away from the other
            return VGroup(lead, lab)

        steps = [
            (fm[0], st.OBS, below(fm[0], st.OBS, "measured")),
            (fm[2], OPTICS, below(fm[2], OPTICS, "optics (CTF)")),
            (rho, st.TRUTH, above(rho, st.TRUTH, "structure", RIGHT)),
            (R_i, st.EST, above(R_i, st.EST, "orientation", LEFT)),
            (fm[6], st.EST, below(fm[6], st.EST, "shift")),
            (fm[8], st.RESID, below(fm[8], st.RESID, "noise")),
        ]
        ref = steps[0][2][1]                      # one baseline row for the brace labels
        for _, _, g in (steps[1], steps[4], steps[5]):
            g[1].align_to(ref, UP)
        for mob, col, lab in steps:
            self.play(mob.animate.set_color(col), FadeIn(lab, shift=0.08 * UP), run_time=0.4)

        inv = VGroup(st.T("Inverse problem: recover", 34, st.TEXT, weight="BOLD"),
                     st.M(r"\rho", size=48, color=st.TRUTH),
                     st.T("from many imperfect", 34, st.TEXT, weight="BOLD"),
                     st.M(r"I_i", size=48, color=st.OBS)).arrange(RIGHT, buff=0.18)
        inv[1].align_to(inv[0], DOWN).shift(0.02 * DOWN)
        inv.move_to([0, 0.1, 0])
        self.play(FadeIn(inv, shift=0.1 * UP), run_time=0.6)

        card = st.takeaway_card(
            None,
            "Each image = CTF-filtered projection at an unknown pose + noise",
            "like a resolution-dependent scale on structure factors,\n"
            "but it oscillates through zero and flips sign")
        card.to_edge(DOWN, buff=st.MARGIN)
        foot4 = st.T("simplified weak-phase / parallel-projection model", 26, st.MUTED,
                     slant="ITALIC")
        foot4.next_to(inv, DOWN, buff=0.15)
        self.play(FadeIn(card, shift=0.15 * UP), FadeIn(foot4), run_time=0.6)
        self.wait(4.0)
