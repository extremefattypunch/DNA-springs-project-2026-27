"""Clip 06 — Gold-standard FSC, local resolution, map-model FSC (~44 s).

Storyboard beats (spec §3, Scene 06) -> code sections:
    0-6.6 s    beat 1  title + question; particle pool split -> half-maps A and B (+ disclaimer,
                       gold-standard caveat: shared mask / initial model / selection)
    6.6-16 s   beat 2  shells S_s on |F_A|, |F_B|; FSC curve drawn shell by shell; 0.143 line;
                       moving vertical readout to the crossing
    16-21 s    beat 3  s = 1/d (top d axis), coarse vs fine low-pass detail, d = 3.2 A (illustrative)
    21-27 s    beat 4  local resolution colouring by blob (local_res_A, LOCRES) + local FSCs
    27-33 s    beat 5  toy atomic model -> simulated rho_model; map-model FSC with 0.5 line
    33-39 s    beat 6  three-column summary + text-only "Paper reports" card
    39-44 s    takeaway freeze (4 s static; why 0.143 vs 0.5 differ)

Every curve and number comes from cached("fsc_resolution", dm.fsc_resolution);
the paper's Table S1 values appear only on the separate text-only card.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from scipy import ndimage

from helpers import styling as st, demos as dm, synthetic_phantom as sp
from helpers.cache import cached

M_SIZE = 38          # on-screen math (>= 36)
B_SIZE = 30          # body text
S_SIZE = 28          # secondary text
F_SIZE = 26          # footnotes / tick labels only

RENDER_EXTENT = 64.0     # dm.cached_render default: +-64 A, 220 px
RENDER_N = 220


# ---------------------------------------------------------------- clip-specific data
def lowpass_sections(d_list=(12.0, 3.2)) -> dict:
    """Ground-truth central section low-pass filtered at each d (hard 2D cutoff
    at s = 1/d); shows what 'finer detail at smaller d' means."""
    sec = cached("fsc_resolution", dm.fsc_resolution)["sec_truth"]
    n = sec.shape[0]
    F = np.fft.fftshift(np.fft.fft2(sec))
    k = (np.arange(n) - n // 2) / n                 # 1/A at 1 A/voxel
    KX, KY = np.meshgrid(k, k, indexing="ij")
    S = np.sqrt(KX ** 2 + KY ** 2)
    out = {}
    for d in d_list:
        out[f"d{d:g}"] = np.real(np.fft.ifft2(np.fft.ifftshift(F * (S <= 1.0 / d))))
    return out


def model_projection(n: int = 128, sigma_A: float = 1.6) -> dict:
    """Toy 'density simulated from the model': the toy model's pseudo-atoms
    projected along the display view (1 A pixels) and splatted with a Gaussian."""
    d = cached("fsc_resolution", dm.fsc_resolution)
    xy = (d["atoms"] @ dm.display_R().T)[:, :2]
    img = np.zeros((n, n))
    idx = np.round(xy + n // 2).astype(int)
    ok = np.all((idx >= 0) & (idx < n), axis=1)
    np.add.at(img, (idx[ok, 0], idx[ok, 1]), 1.0)
    return {"img": ndimage.gaussian_filter(img, sigma_A)}


def hexcol(rgb) -> str:
    return "#%02X%02X%02X" % tuple(int(c) for c in rgb)


def crop_window(render: dict, pad_A: float = 5.0):
    """Square window (centre x, centre y, side, all in A) around the rendered density."""
    ax = np.linspace(-RENDER_EXTENT, RENDER_EXTENT, RENDER_N)
    ii, jj = np.nonzero(render["mask"])
    x0, x1, y0, y1 = ax[ii.min()], ax[ii.max()], ax[jj.min()], ax[jj.max()]
    side = max(x1 - x0, y1 - y0) + 2 * pad_A
    return 0.5 * (x0 + x1), 0.5 * (y0 + y1), side


def crop_xy(arr: np.ndarray, coords: np.ndarray, win) -> np.ndarray:
    """Crop an (x, y)-indexed square array whose pixel centres sit at `coords` (A)."""
    cx, cy, side = win
    ix = np.nonzero(np.abs(coords - cx) <= side / 2)[0]
    iy = np.nonzero(np.abs(coords - cy) <= side / 2)[0]
    m = min(len(ix), len(iy))
    return arr[ix[0]:ix[0] + m, iy[0]:iy[0] + m]


class Scene06FSCResolution(Scene):
    def construct(self):
        data = cached("fsc_resolution", dm.fsc_resolution)
        lp = cached("s06_fsc_resolution_lowpass", lowpass_sections, d_list=[12.0, 3.2])
        mp = cached("s06_fsc_resolution_modelproj", model_projection, n=128, sigma_A=1.6)
        render_full = dm.cached_render()                 # default side-on display view
        win = crop_window(render_full)
        r_coords = np.linspace(-RENDER_EXTENT, RENDER_EXTENT, RENDER_N)
        render = {k: crop_xy(render_full[k], r_coords, win) for k in ("shade", "mask", "label")}
        r_side_A = render["shade"].shape[0] * (r_coords[1] - r_coords[0])
        r_center_A = np.array([crop_xy(np.tile(r_coords[:, None], (1, RENDER_N)), r_coords, win).mean(),
                               crop_xy(np.tile(r_coords[None, :], (RENDER_N, 1)), r_coords, win).mean()])

        s, fsc_half, fsc_model = data["s"], data["fsc_half"], data["fsc_model"]
        s_half = float(data["half_crossing_s"])
        s_model = float(data["model_crossing_s"])
        d_half, d_model = 1 / s_half, 1 / s_model
        loc = data["local_res_A"]
        pos_true, _ = sp.pseudo_atoms("A")
        rms_err = float(np.sqrt(((data["atoms"] - pos_true) ** 2).sum(1).mean()))

        def axes(w, h, y_rng, xnums, ynums, center):
            a = Axes(x_range=[0, 0.5, 0.1], y_range=y_rng, x_length=w, y_length=h, tips=False,
                     axis_config={"color": st.MUTED, "stroke_width": 2},
                     x_axis_config={"numbers_to_include": xnums, "font_size": F_SIZE,
                                    "decimal_number_config": {"num_decimal_places": 1}},
                     y_axis_config={"numbers_to_include": ynums, "font_size": F_SIZE,
                                    "decimal_number_config": {"num_decimal_places": 1}})
            return a.move_to(center)

        def s_label(a):
            lab = st.M(r"s\ (\text{\AA}^{-1})", size=36, color=st.MUTED)
            return lab.next_to(a.x_axis, DOWN, buff=0.15).align_to(a.x_axis, RIGHT)

        def curve_upto(xs, ys, a, t, color, width=5):
            m = xs <= t
            px, py = list(xs[m]), list(ys[m])
            if t > xs[0] and not m.all():
                j = int(m.sum())
                f = (t - xs[j - 1]) / (xs[j] - xs[j - 1])
                px.append(t)
                py.append(ys[j - 1] + f * (ys[j] - ys[j - 1]))
            pts = [a.c2p(u, v) for u, v in zip(px, py)]
            if len(pts) < 2:
                pts = [a.c2p(xs[0], ys[0])] * 2
            v = VMobject(stroke_color=color, stroke_width=width)
            v.set_points_as_corners(pts)
            return v

        # ------------------------------------------------ opening
        title = st.title_card("Resolution: global, local, model",
                              "What do global map resolution, local resolution, and\n"
                              "model resolution actually measure?")
        badge = st.simulation_badge()
        self.add(badge)
        self.play(FadeIn(title, shift=0.1 * DOWN), run_time=0.8)

        # ------------------------------------------------ beat 1: independent halves
        rng = np.random.default_rng(2027 + 6)
        parts = VGroup()
        for _ in range(24):
            parts.add(RoundedRectangle(corner_radius=0.04, width=0.34, height=0.34,
                                       stroke_color=st.OBS, stroke_width=2, fill_color=st.OBS,
                                       fill_opacity=0.15 + 0.25 * rng.random()))
        parts.arrange_in_grid(rows=4, cols=6, buff=0.12).move_to([-4.6, -0.3, 0])
        pool_lab = st.T("particle pool", S_SIZE, st.OBS).next_to(parts, UP, buff=0.25)

        img_vmin = min(data["sec_A"].min(), data["sec_B"].min())
        img_vmax = max(data["sec_A"].max(), data["sec_B"].max())
        tsec = data["sec_truth"]
        ii, jj = np.nonzero(tsec > 0.05 * tsec.max())
        half = max(ii.max() - ii.min(), jj.max() - jj.min()) // 2 + 10
        ci, cj = (ii.max() + ii.min()) // 2, (jj.max() + jj.min()) // 2
        sl = (slice(max(ci - half, 0), ci + half), slice(max(cj - half, 0), cj + half))
        hA = st.array_to_image_mobject(data["sec_A"][sl], height=2.1, vmin=img_vmin, vmax=img_vmax)
        hB = st.array_to_image_mobject(data["sec_B"][sl], height=2.1, vmin=img_vmin, vmax=img_vmax)
        rowA, rowB = 0.75, -1.65
        hA.move_to([-0.9, rowA, 0])
        hB.move_to([-0.9, rowB, 0])
        fA, fB = st.framed(hA, st.EST), st.framed(hB, st.EST)
        labA = st.T("Half-map A", B_SIZE, st.EST, weight="BOLD").next_to(hA, RIGHT, buff=0.35).shift(0.4 * UP)
        labB = st.T("Half-map B", B_SIZE, st.EST, weight="BOLD").next_to(hB, RIGHT, buff=0.35).shift(0.4 * UP)
        eqA = st.M(r"\rho_A=\rho+n_A", size=M_SIZE).next_to(labA, DOWN, aligned_edge=LEFT, buff=0.25)
        eqB = st.M(r"\rho_B=\rho+n_B", size=M_SIZE).next_to(labB, DOWN, aligned_edge=LEFT, buff=0.25)
        eqA[0][-2:].set_color(st.RESID)
        eqB[0][-2:].set_color(st.RESID)
        indep = VGroup(st.M(r"n_A,\ n_B", size=M_SIZE, color=st.RESID),
                       st.T("independent noise", S_SIZE, st.RESID)).arrange(DOWN, buff=0.12)
        indep.move_to([4.75, rowA - 0.05, 0])

        idxA = [k for k in range(24) if k % 2 == 0]
        idxB = [k for k in range(24) if k % 2 == 1]
        grpA = VGroup(*[parts[k].copy() for k in idxA]).arrange_in_grid(rows=3, cols=4, buff=0.1)
        grpB = VGroup(*[parts[k].copy() for k in idxB]).arrange_in_grid(rows=3, cols=4, buff=0.1)
        grpA.move_to([-5.4, rowA - 0.15, 0])
        grpB.move_to([-5.4, rowB - 0.15, 0])
        setA_lab = st.T("particle set A", S_SIZE, st.OBS).next_to(grpA, UP, buff=0.15)
        setB_lab = st.T("particle set B", S_SIZE, st.OBS).next_to(grpB, UP, buff=0.15)
        arA = Arrow(grpA.get_right(), hA.get_left(), buff=0.15, color=st.MUTED, stroke_width=3)
        arB = Arrow(grpB.get_right(), hB.get_left(), buff=0.15, color=st.MUTED, stroke_width=3)
        refine_lab = VGroup(*[st.T("refined\nseparately", S_SIZE, st.MUTED, line_spacing=0.8)
                              .next_to(ar, UP, buff=0.06).align_to(ar, RIGHT)
                              for ar in (arA, arB)])
        gold = st.T("Gold standard: separate refinement, so noise is not fitted in common",
                    S_SIZE, st.TEXT).move_to([0, -2.85, 0])
        gcave = st.T("caveat: shared mask,\ninitial model and\nparticle selection\n"
                     "still need validation", S_SIZE, st.MUTED, slant="ITALIC",
                     line_spacing=0.8).next_to(indep, DOWN, buff=0.3)
        disclaim = st.phantom_disclaimer().to_edge(DOWN, buff=st.MARGIN)

        self.play(LaggedStart(*[FadeIn(p, scale=0.6) for p in parts], lag_ratio=0.03),
                  FadeIn(pool_lab), run_time=1.1)
        moves = [Transform(parts[k], grpA[i]) for i, k in enumerate(idxA)]
        moves += [Transform(parts[k], grpB[i]) for i, k in enumerate(idxB)]
        self.play(*moves, FadeOut(pool_lab), FadeIn(setA_lab), FadeIn(setB_lab), run_time=1.2)
        self.play(GrowArrow(arA), GrowArrow(arB), FadeIn(refine_lab),
                  FadeIn(hA), FadeIn(hB), Create(fA), Create(fB),
                  FadeIn(labA), FadeIn(labB), FadeIn(disclaim), run_time=1.2)
        self.play(Write(eqA), Write(eqB), FadeIn(indep), run_time=1.0)
        self.play(FadeIn(gold), FadeIn(gcave), run_time=0.5)
        self.wait(0.8)

        # ------------------------------------------------ beat 2: shells and the FSC curve
        q = title[1]
        ftA = st.array_to_image_mobject(data["ft_A"], height=2.5, stops=st.MAGMA)
        ftB = st.array_to_image_mobject(data["ft_B"], height=2.5, stops=st.MAGMA)
        ftA.move_to([-5.05, -0.95, 0])
        ftB.move_to([-2.3, -0.95, 0])
        fFA, fFB = st.framed(ftA, st.EST), st.framed(ftB, st.EST)
        ftA_lab = st.M(r"|F_A|", size=M_SIZE, color=st.EST).next_to(ftA, UP, buff=0.15)
        ftB_lab = st.M(r"|F_B|", size=M_SIZE, color=st.EST).next_to(ftB, UP, buff=0.15)
        ft_note = VGroup(st.T("log |F| of the central slice (2D view)", S_SIZE, st.MUTED),
                         st.T("(the FSC uses 3D shells, masked)", S_SIZE, st.MUTED)
                         ).arrange(DOWN, buff=0.08)
        ft_note.next_to(Group(ftA, ftB), DOWN, buff=0.2)

        fsc_eq = st.M(r"\mathrm{FSC}(s)=\frac{\mathrm{Re}\sum_{\mathbf k\in S_s}F_A(\mathbf k)\,"
                      r"F_B(\mathbf k)^*}{\sqrt{\sum_{\mathbf k\in S_s}|F_A(\mathbf k)|^2"
                      r"\sum_{\mathbf k\in S_s}|F_B(\mathbf k)|^2}}", size=M_SIZE)
        fsc_eq.move_to([0, 1.95, 0])

        ax = axes(5.3, 3.5, [-0.1, 1.0, 0.2], [0.1, 0.2, 0.3, 0.4, 0.5],
                  [0, 0.2, 0.4, 0.6, 0.8, 1.0], [3.45, -0.95, 0])
        ax_xlab = s_label(ax)
        ax_ylab = st.M(r"\mathrm{FSC}", size=36, color=st.MUTED).rotate(PI / 2)
        ax_ylab.next_to(ax.y_axis, LEFT, buff=0.55)

        st_tr = ValueTracker(float(s[0]))
        half_curve = always_redraw(lambda: curve_upto(s, fsc_half, ax, st_tr.get_value(), st.EST))

        def shell_ring(img):
            def f():
                r = st_tr.get_value() / 0.5 * img.height / 2
                return Circle(radius=max(r, 0.01), stroke_color=st.TEXT, stroke_width=3).move_to(img)
            return always_redraw(f)

        faint = VGroup()
        for img in (ftA, ftB):
            for sv in np.arange(0.1, 0.5, 0.1):
                faint.add(Circle(radius=sv / 0.5 * img.height / 2, stroke_color=st.MUTED,
                                 stroke_width=1.2, stroke_opacity=0.6).move_to(img))
        ringA, ringB = shell_ring(ftA), shell_ring(ftB)
        shell_lab = st.M(r"S_s", size=M_SIZE, color=st.TEXT)
        shell_lab.next_to(ftB, RIGHT, buff=0.12).align_to(ftB, UP)

        self.play(FadeOut(Group(parts, setA_lab, setB_lab, arA, arB, refine_lab, eqA, eqB,
                                indep, gold, gcave, disclaim, q, labA, labB)), run_time=0.6)
        self.play(FadeTransform(hA, ftA), FadeTransform(hB, ftB),
                  ReplacementTransform(fA, fFA), ReplacementTransform(fB, fFB),
                  FadeIn(ftA_lab), FadeIn(ftB_lab), FadeIn(ft_note), run_time=1.0)
        self.play(Write(fsc_eq), Create(faint), run_time=1.2)
        self.play(Create(ax), FadeIn(ax_xlab), FadeIn(ax_ylab), run_time=0.7)
        self.add(half_curve, ringA, ringB, shell_lab)
        self.play(st_tr.animate.set_value(float(s[-1])), run_time=2.4, rate_func=linear)
        self.remove(ringA, ringB, shell_lab, half_curve)
        static_curve = curve_upto(s, fsc_half, ax, float(s[-1]) + 1, st.EST)
        self.add(static_curve)

        # threshold and the moving readout
        thr = DashedLine(ax.c2p(0, 0.143), ax.c2p(0.5, 0.143), color=st.TEXT, stroke_width=2.5,
                         dash_length=0.12)
        thr_lab = st.M(r"0.143", size=36).next_to(ax.c2p(0.5, 0.143), UP, buff=0.12).shift(0.45 * LEFT)
        ro = ValueTracker(0.12)

        def fsc_at(x):
            return float(np.interp(x, s, fsc_half))

        vline = always_redraw(lambda: DashedLine(ax.c2p(ro.get_value(), -0.1),
                                                 ax.c2p(ro.get_value(), 1.0),
                                                 color=st.EST, stroke_width=2.5, dash_length=0.1))
        vdot = always_redraw(lambda: Dot(ax.c2p(ro.get_value(), fsc_at(ro.get_value())),
                                         radius=0.08, color=st.EST))
        num = DecimalNumber(fsc_at(0.12), num_decimal_places=3, font_size=36, color=st.EST)
        num_lab = st.M(r"\mathrm{FSC}=", size=36, color=st.EST)
        readout = VGroup(num_lab, num)

        def place_readout(m):
            num.set_value(fsc_at(ro.get_value()))
            num.next_to(num_lab, RIGHT, buff=0.12)
            m.next_to(ax.c2p(ro.get_value(), 1.0), UP, buff=0.1)
        readout.add_updater(place_readout)
        self.play(Create(thr), FadeIn(thr_lab), run_time=0.5)
        self.add(vline, vdot, readout)
        self.play(ro.animate.set_value(s_half), run_time=1.4, rate_func=smooth)
        self.wait(0.8)

        # ------------------------------------------------ beat 3: s = 1/d and the crossing
        readout.clear_updaters()
        self.play(FadeOut(Group(ftA, ftB, fFA, fFB, ftA_lab, ftB_lab, ft_note, faint, fsc_eq,
                                readout)), run_time=0.5)
        # secondary d axis above the plot
        top_y = 1.13
        dax = VGroup(Line(ax.c2p(0, top_y), ax.c2p(0.5, top_y), color=st.MUTED, stroke_width=2))
        for dv in (10, 5, 4, 3):
            p = ax.c2p(1 / dv, top_y)
            dax.add(Line(p, p + 0.12 * UP, color=st.MUTED, stroke_width=2))
            dax.add(st.T(f"{dv:g}", F_SIZE, st.MUTED).next_to(p, UP, buff=0.17))
        d_lab = st.M(r"d\ (\text{\AA})", size=36, color=st.MUTED)
        d_lab.move_to(ax.c2p(0.45, top_y)).align_to(dax[-1], DOWN)
        sd_eq = st.M(r"s=\frac{1}{d}", size=48).move_to([-5.3, 2.0, 0])
        finer = VGroup(st.T("higher s  =  smaller d", B_SIZE),
                       st.T("=  finer detail", B_SIZE)).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        finer.next_to(sd_eq, RIGHT, buff=0.45)
        arrow_y = ax_xlab.get_y()
        dir_arrow = Arrow([ax.c2p(0.03, 0)[0], arrow_y, 0], [ax.c2p(0.33, 0)[0], arrow_y, 0],
                          buff=0, color=st.TEXT, stroke_width=3,
                          max_tip_length_to_length_ratio=0.08)

        lp_vmin = min(lp["d12"].min(), lp["d3.2"].min())
        lp_vmax = max(lp["d12"].max(), lp["d3.2"].max())
        im8 = st.array_to_image_mobject(lp["d12"], height=2.2, vmin=lp_vmin, vmax=lp_vmax)
        im32 = st.array_to_image_mobject(lp["d3.2"], height=2.2, vmin=lp_vmin, vmax=lp_vmax)
        im8.move_to([-5.2, -0.5, 0])
        im32.move_to([-2.4, -0.5, 0])
        fr8, fr32 = st.framed(im8, st.TRUTH), st.framed(im32, st.TRUTH)
        im8_lab = st.T("d = 12 Å", S_SIZE, st.TRUTH).next_to(im8, UP, buff=0.12)
        im32_lab = st.T("d = 3.2 Å", S_SIZE, st.TRUTH).next_to(im32, UP, buff=0.12)
        lp_note = st.T("ground truth, low-pass cut at d", S_SIZE, st.MUTED)
        lp_note.next_to(Group(im8, im32), DOWN, buff=0.15)

        cross_txt = VGroup(
            st.M(r"s=%.4f\ \text{\AA}^{-1}\ \Rightarrow\ d=1/s=%.1f\ \text{\AA}" % (s_half, d_half),
                 size=38, color=st.EST),
            st.T("illustrative: engineered on synthetic data", F_SIZE, st.MUTED, slant="ITALIC"),
        ).arrange(DOWN, buff=0.12)
        cross_box = VGroup(SurroundingRectangle(cross_txt, color=st.EST, buff=0.18,
                                                corner_radius=0.08, fill_color=st.PANEL,
                                                fill_opacity=0.95, stroke_width=2), cross_txt)
        cross_box.move_to([-3.25, -2.85, 0])

        self.play(Write(sd_eq), Create(dax), FadeIn(d_lab), run_time=0.9)
        self.play(FadeIn(finer), GrowArrow(dir_arrow), run_time=0.7)
        self.play(FadeIn(im8), Create(fr8), FadeIn(im8_lab),
                  FadeIn(im32), Create(fr32), FadeIn(im32_lab), FadeIn(lp_note), run_time=0.8)
        self.play(FadeIn(cross_box, shift=0.1 * UP), run_time=0.7)
        self.wait(1.5)

        # ------------------------------------------------ beat 4: local resolution
        vline.clear_updaters()
        vdot.clear_updaters()
        self.play(FadeOut(Group(sd_eq, finer, dir_arrow, im8, im32, fr8, fr32, im8_lab, im32_lab,
                                lp_note, cross_box, dax, d_lab, ax, ax_xlab, ax_ylab,
                                static_curve, thr, thr_lab, vline, vdot)), run_time=0.6)
        lo = np.floor(loc.min() * 10) / 10
        hi = np.ceil(loc.max() * 10) / 10
        lab_cols = {i: hexcol(st.colorize(np.array([(loc[i] - lo) / (hi - lo)]), st.LOCRES)[0])
                    for i in range(len(loc))}
        locimg = st.density_image(render, st.MUTED, height=5.0, label_colors=lab_cols)
        locimg.move_to([-4.2, -0.7, 0])
        bar_arr = np.tile(np.linspace(0, 1, 120)[None, :], (6, 1))
        bar = st.array_to_image_mobject(bar_arr, height=3.4, stops=st.LOCRES)
        bar.stretch_to_fit_width(0.3)
        bar.next_to(locimg, RIGHT, buff=0.3).set_y(-0.9)
        bar_fr = st.framed(bar, st.MUTED)
        bar_lo = st.T(f"{lo:.1f} Å", S_SIZE).next_to(bar, RIGHT, buff=0.12).align_to(bar, DOWN)
        bar_hi = st.T(f"{hi:.1f} Å", S_SIZE).next_to(bar, RIGHT, buff=0.12).align_to(bar, UP)
        bar_t = st.T("local res.", S_SIZE, st.MUTED).next_to(bar, UP, buff=0.15)
        loc_head = st.T("Different parts of the same map can be resolved differently",
                        B_SIZE, st.TEXT).move_to([0, 2.6, 0])
        disclaim2 = st.phantom_disclaimer().to_corner(DL, buff=st.MARGIN)

        core_i = int(np.argmin(loc))
        mid_i = 1
        app_i = 15
        lax = axes(4.5, 2.1, [-0.2, 1.0, 0.2], [0.1, 0.3, 0.5], [0, 0.5, 1.0], [4.05, 0.7, 0])
        lax_x = s_label(lax)
        lax_t = st.T("FSC inside a small window", S_SIZE, st.MUTED).next_to(lax, UP, buff=0.15)
        lthr = DashedLine(lax.c2p(0, 0.143), lax.c2p(0.5, 0.143), color=st.TEXT, stroke_width=2,
                          dash_length=0.1)
        lcurves = VGroup()
        legend_rows = []
        for i, name in ((core_i, "rigid core"), (mid_i, "partly flexible"), (app_i, "appendage")):
            lcurves.add(curve_upto(s, data["local_curves"][i], lax, 1.0, lab_cols[i], width=4))
            legend_rows.append((lab_cols[i], f"{name}: {loc[i]:.1f} Å"))
        leg = st.legend(legend_rows, size=S_SIZE)
        leg.next_to(lax, DOWN, buff=0.08).align_to(lax, LEFT)
        lthr_lab = st.M(r"0.143", size=36).next_to(lax.c2p(0.5, 0.143), UP, buff=0.08).shift(0.4 * LEFT)
        caveat = st.T("local res. is not a per-atom\nB factor or coordinate error",
                      S_SIZE, st.MUTED, slant="ITALIC", line_spacing=0.8)
        caveat.next_to(leg, DOWN, buff=0.1).to_edge(RIGHT, buff=st.MARGIN)

        # callouts: which coloured blob is the appendage and which is the rigid core
        b_xy = (data["blob_centers"] @ dm.display_R().T)[:, :2]
        l_unit = locimg.height / r_side_A
        l_cen = locimg.get_center()

        def blob_pt(i):
            return np.array([l_cen[0] + (b_xy[i, 0] - r_center_A[0]) * l_unit,
                             l_cen[1] + (b_xy[i, 1] - r_center_A[1]) * l_unit, 0.0])
        app_txt = st.T(f"appendage {loc[app_i]:.1f} Å", S_SIZE, st.MUTED).move_to([-5.15, 1.75, 0])
        core_txt = st.T(f"rigid core {loc[core_i]:.1f} Å", S_SIZE, st.MUTED).move_to([-4.55, -0.75, 0])
        app_ptr = Line(app_txt.get_corner(DR) + 0.05 * DOWN, blob_pt(app_i), stroke_color=st.MUTED,
                       stroke_width=2)
        core_ptr = Line(core_txt.get_corner(UR) + 0.05 * UP, blob_pt(core_i), stroke_color=st.MUTED,
                        stroke_width=2)
        callouts = VGroup(app_txt, core_txt, app_ptr, core_ptr,
                          Dot(blob_pt(app_i), radius=0.05, color=st.MUTED),
                          Dot(blob_pt(core_i), radius=0.05, color=st.MUTED))

        self.play(FadeIn(locimg), FadeIn(disclaim2), FadeIn(bar), Create(bar_fr), FadeIn(bar_lo),
                  FadeIn(bar_hi), FadeIn(bar_t), run_time=0.9)
        self.play(FadeIn(loc_head, shift=0.1 * DOWN), run_time=0.6)
        self.play(Create(lax), FadeIn(lax_x), FadeIn(lax_t), Create(lthr), FadeIn(lthr_lab),
                  run_time=0.6)
        self.play(LaggedStart(*[Create(c) for c in lcurves], lag_ratio=0.3), FadeIn(leg),
                  run_time=1.3)
        self.play(FadeIn(callouts), run_time=0.5)
        self.play(FadeIn(caveat), run_time=0.4)
        self.wait(1.3)

        # ------------------------------------------------ beat 5: map-model FSC
        self.play(FadeOut(Group(bar, bar_fr, bar_lo, bar_hi, bar_t, loc_head, lax, lax_x, lax_t,
                                lthr, lthr_lab, lcurves, leg, caveat, callouts, disclaim2,
                                locimg)), run_time=0.5)
        mapimg = st.density_image(render, st.MUTED, height=4.4)
        mapimg.move_to([-4.65, -0.35, 0])
        xy = (data["atoms"] @ dm.display_R().T)[:, :2]
        unit = mapimg.height / r_side_A
        cen = mapimg.get_center()
        P = np.column_stack([cen[0] + (xy[:, 0] - r_center_A[0]) * unit,
                             cen[1] + (xy[:, 1] - r_center_A[1]) * unit, np.zeros(len(xy))])
        owner = data["atom_owner"]
        pairs = set()
        for i in range(len(P)):          # rod from each atom to its nearest same-domain atom
            same = np.nonzero(owner == owner[i])[0]
            same = same[same != i]
            j = int(same[np.argmin(np.linalg.norm(P[same] - P[i], axis=1))])
            pairs.add((min(i, j), max(i, j)))
        rods = VGroup(*[Line(P[i], P[j], stroke_color=st.EST, stroke_width=1.4)
                        for i, j in sorted(pairs)])
        dots = VGroup(*[Dot(p, radius=0.022, color=st.EST) for p in P])
        model = VGroup(rods, dots)
        map_lab = st.T("map (smooth render)", S_SIZE, st.MUTED)
        model_lab = st.T("+ toy atomic model", S_SIZE, st.EST)
        map_lab_g = VGroup(map_lab, model_lab).arrange(DOWN, buff=0.08).next_to(mapimg, UP, buff=0.1)

        mp_coords = np.arange(128) - 64.0
        mimg = st.array_to_image_mobject(crop_xy(mp["img"], mp_coords, win), height=2.4,
                                         stops=st.GRAY)
        mimg.move_to([-1.15, 0.0, 0])
        mfr = st.framed(mimg, st.EST)
        m_lab = st.M(r"\rho_{\text{model}}", size=M_SIZE, color=st.EST).next_to(mimg, UP, buff=0.12)
        m_note = st.T("density simulated\nfrom the model", S_SIZE, st.MUTED, line_spacing=0.8)
        m_note.next_to(mimg, DOWN, buff=0.15)

        max_ = axes(5.0, 3.1, [-0.1, 1.0, 0.2], [0.1, 0.2, 0.3, 0.4, 0.5], [0, 0.5, 1.0],
                    [3.45, 0.15, 0])
        mx_lab = s_label(max_)
        my_lab = st.M(r"\mathrm{FSC}(\rho_{\text{model}},\ \tfrac{1}{2}(\rho_A+\rho_B))", size=36,
                      color=st.TEXT).next_to(max_, UP, buff=0.15)
        ghost = curve_upto(s, fsc_half, max_, 1.0, st.EST, width=3).set_stroke(opacity=0.4)
        mcurve = curve_upto(s, fsc_model, max_, 1.0, st.TEXT, width=5)
        mcurve_lab = st.T("model vs map", S_SIZE, st.TEXT)
        mcurve_lab.move_to(max_.c2p(0.295, 1.0), aligned_edge=UL)
        ghost_lab = VGroup(st.T("half-map", S_SIZE, st.EST),
                           st.T("(read at 0.143)", S_SIZE, st.EST)
                           ).arrange(DOWN, buff=0.06, aligned_edge=LEFT).set_opacity(0.75)
        ghost_lab.next_to(mcurve_lab, DOWN, buff=0.1, aligned_edge=LEFT)
        thr5 = DashedLine(max_.c2p(0, 0.5), max_.c2p(0.5, 0.5), color=st.TEXT, stroke_width=2.5,
                          dash_length=0.12)
        thr5_lab = st.M(r"0.5", size=36).next_to(max_.c2p(0.5, 0.5), DOWN, buff=0.1).shift(0.3 * LEFT)
        mv = DashedLine(max_.c2p(s_model, -0.1), max_.c2p(s_model, 1.0), color=st.TEXT,
                        stroke_width=2.5, dash_length=0.1)
        mdot = Dot(max_.c2p(s_model, 0.5), radius=0.08, color=st.TEXT)
        m_read = VGroup(st.M(r"d_{\text{model}}=%.1f\ \text{\AA}" % d_model, size=36),
                        st.T("(illustrative)", F_SIZE, st.MUTED, slant="ITALIC")
                        ).arrange(DOWN, buff=0.06)
        m_read.move_to(max_.c2p(s_model, 0.18), aligned_edge=RIGHT).shift(0.12 * LEFT)

        state1 = st.T("This tests model-map agreement; it is not the half-map FSC.", B_SIZE)
        state2 = st.T(f"0.5 cutoff ({d_model:.1f} Å) is not coordinate error: "
                      f"toy model RMS error = {rms_err:.1f} Å", S_SIZE, st.MUTED)
        states = VGroup(state1, state2).arrange(DOWN, buff=0.15).move_to([0, -3.05, 0])

        self.play(FadeIn(mapimg), FadeIn(map_lab), run_time=0.5)
        self.play(LaggedStart(Create(rods), FadeIn(dots), lag_ratio=0.2), FadeIn(model_lab),
                  run_time=0.9)
        self.play(FadeIn(mimg, shift=0.2 * RIGHT), Create(mfr), FadeIn(m_lab), FadeIn(m_note),
                  run_time=0.7)
        self.play(Create(max_), FadeIn(mx_lab), FadeIn(my_lab), FadeIn(ghost), FadeIn(ghost_lab),
                  run_time=0.6)
        self.play(Create(mcurve), FadeIn(mcurve_lab), Create(thr5), FadeIn(thr5_lab), run_time=0.9)
        self.play(Create(mv), FadeIn(mdot), FadeIn(m_read), FadeIn(state1), run_time=0.6)
        self.play(FadeIn(state2), run_time=0.4)
        self.wait(1.1)

        # ------------------------------------------------ beat 6: summary + paper card
        self.play(FadeOut(Group(map_lab_g, model, mimg, mfr, m_lab, m_note, max_,
                                mx_lab, my_lab, ghost, ghost_lab, mcurve, mcurve_lab, thr5,
                                thr5_lab, mv, mdot, m_read, states, mapimg)), run_time=0.5)

        def column(head, color, lines, x):
            h = st.T(head, B_SIZE, color, weight="BOLD")
            body = VGroup(*[st.T(t, S_SIZE) for t in lines]).arrange(DOWN, buff=0.1)
            g = VGroup(h, body).arrange(DOWN, buff=0.18)
            box = RoundedRectangle(corner_radius=0.12, width=4.15, height=g.height + 0.4,
                                   stroke_color=color, stroke_width=2.5, fill_color=st.PANEL,
                                   fill_opacity=0.9)
            g.move_to(box)
            return VGroup(box, g).move_to([x, 1.75, 0])

        col1 = column("Global map", st.EST, ["map vs map", "half-map A vs B", "FSC = 0.143"], -4.4)
        col2 = column("Local", lab_cols[core_i], ["region-varying detail", "FSC in small windows",
                                                 "a range, not one value"], 0.0)
        col3 = column("Model", st.TEXT, ["model vs map", "model map vs full map", "FSC = 0.5"], 4.4)
        cols = VGroup(col1, col2, col3)
        cc = st.T("half-map FSC ~ CC1/2 in spirit, not Rfree", S_SIZE, st.MUTED, slant="ITALIC")
        cc.next_to(cols, DOWN, buff=0.2)

        # text-only card with the paper's Table S1 values (published, not simulated)
        rows = [("", "Aα holoenzyme", "Aβ holoenzyme"),
                ("Global (half-map FSC = 0.143)", "3.4 Å", "3.2 Å"),
                ("Local resolution range", "3.04–5.87 Å", "2.82–5.22 Å"),
                ("Model (map-model FSC = 0.5)", "3.8 Å", "3.6 Å")]
        tab = VGroup()
        for r_i, r in enumerate(rows):
            col = st.MUTED if r_i == 0 else st.TEXT
            cells = VGroup()
            for c in r:          # trailing invisible '|' gives every cell the same line box
                t = st.T(c + "|", S_SIZE, col)
                t[-1].set_opacity(0)
                cells.add(t)
            tab.add(cells)
        widths = [max(row[c].width for row in tab) for c in range(3)]
        xs = [0, widths[0] + 0.7, widths[0] + widths[1] + 1.4]
        for r_i, row in enumerate(tab):
            for c in range(3):
                row[c].move_to([xs[c], -0.46 * r_i, 0], aligned_edge=LEFT)
        paper_head = VGroup(st.T("Paper reports (Day et al., 2025; Table S1)", S_SIZE, st.TEXT,
                                 weight="BOLD"),
                            st.T("published values - not from these simulations", F_SIZE,
                                 st.MUTED, slant="ITALIC")).arrange(DOWN, buff=0.08,
                                                                    aligned_edge=LEFT)
        paper_body = VGroup(paper_head, tab).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
        paper_box = RoundedRectangle(corner_radius=0.12, width=paper_body.width + 0.6,
                                     height=paper_body.height + 0.45, stroke_color=st.MUTED,
                                     stroke_width=1.5, fill_color=st.PANEL, fill_opacity=0.95)
        paper_body.move_to(paper_box)
        paper = VGroup(paper_box, paper_body).move_to([0, -1.85, 0])

        self.play(LaggedStart(FadeIn(col1, shift=0.1 * UP), FadeIn(col2, shift=0.1 * UP),
                              FadeIn(col3, shift=0.1 * UP), lag_ratio=0.25), run_time=1.1)
        self.play(FadeIn(cc), run_time=0.4)
        self.play(FadeIn(paper), run_time=0.6)
        self.wait(3.0)

        # ------------------------------------------------ takeaway freeze
        take = st.takeaway_card(
            r"\mathrm{FSC}_{A,B}=0.143\Rightarrow d_{\text{map}}\qquad"
            r"\mathrm{FSC}_{\text{model,map}}=0.5\Rightarrow d_{\text{model}}",
            "placeholder",
            "half-map FSC ~ CC1/2 in spirit, not Rfree")
        # rebuild the text row with real line spacing (Text's default is too tight for Å)
        t_body = take[1]
        t_lines = VGroup(*[st.T(t, S_SIZE) for t in (
            "Half-map FSC = 0.143 ~ full map vs a perfect reference at 0.5.",
            f"An imperfect model reaches 0.5 sooner ({d_model:.1f} Å vs {d_half:.1f} Å here),",
            f"and one global value hides worse local regions (up to {loc.max():.1f} Å here).")]
                         ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        t_body.submobjects[2] = t_lines
        t_body.arrange(DOWN, buff=0.2)
        take[0].stretch_to_fit_height(t_body.height + 0.45)
        take[0].stretch_to_fit_width(max(t_body.width + 0.6, 6))
        t_body.move_to(take[0])
        take.to_edge(DOWN, buff=st.MARGIN + 0.05)
        self.play(FadeOut(paper), FadeOut(cc), run_time=0.4)
        self.play(FadeIn(take, shift=0.1 * UP), run_time=0.6)
        self.wait(4.0)
