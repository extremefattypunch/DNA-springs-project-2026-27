"""Clip 04 — Iterative pose / source refinement and regularization.

Storyboard beats (spec section 3, Scene 04) -> code sections:
  1  0-9 s    initial mismatch: noisy observation vs prediction from the current
              estimate; orange estimated viewing direction vs green reference arrow;
              in-plane shift error arrow; schematic forward-model label
  2  9-19 s   pose update for the tracked observation: arrows converge, shift arrow
              shrinks, residual image and its norm fall to the noise level
  3  19-29 s  source update: all pose errors (bars) and the five estimated dots,
              alternating pose and source steps; computed loss vs iteration
  4  29-36 s  regularization: two 1D fits (no penalty vs penalty); smaller data
              misfit is not a smaller true error
  5  36-41 s  freeze: estimate poses -> update common source -> avoid overfitting

All numbers come from cached("refinement", ne.refinement): a toy alternating
minimisation of L = sum_i ||I_i - P_{R_i,t_i}[X]||^2 + lam * sum_j |x_j|^2 (CTF = 1).
The tracked pose and its shift (beats 1-2) are shown after ONE fixed superposition on the
reference (global_frame), so source steps leave the displayed pose unchanged. In beat 3 the
pose bars and the dots both use the per-record superposition of the data dict, so the two
panels share one reference frame at every step.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, projections as pj
from helpers.cache import cached


# ------------------------------------------------------------------ geometry
def cam(elev: float, azim: float) -> np.ndarray:
    """3x2 matrix: fixed orthographic view of 3D vectors (screen = v @ B)."""
    e, a = np.radians(elev), np.radians(azim)
    right = np.array([-np.sin(a), np.cos(a), 0.0])
    up = np.array([-np.sin(e) * np.cos(a), -np.sin(e) * np.sin(a), np.cos(e)])
    return np.stack([right, up], 1)


def to3(xy) -> np.ndarray:
    return np.array([xy[0], xy[1], 0.0])


def axes_3d(B: np.ndarray, origin: np.ndarray, scale: float, length: float = 1.25):
    """Three short grey coordinate axes in a fixed oblique view (lines + labels)."""
    g = VGroup()
    for v, name in zip(np.eye(3), "xyz"):
        end = origin + scale * length * to3(v @ B)
        ln = Line(origin, end, color=st.GRID, stroke_width=2.5)
        lab = st.M(name, size=36, color=st.MUTED)
        direction = end - origin
        lab.move_to(end + 0.28 * direction / max(np.linalg.norm(direction), 1e-6))
        g.add(ln, lab)
    return g


def axes_2d(x_range, y_range, w, h, x_ticks, y_ticks, x_fmt="{:g}", y_fmt="{:g}"):
    """Plain axes with readable (size 26) tick labels."""
    ax = Axes(x_range=x_range, y_range=y_range, x_length=w, y_length=h, tips=False,
              axis_config={"color": st.GRID, "stroke_width": 2, "include_ticks": False})
    labels = VGroup()
    for x in x_ticks:
        labels.add(st.T(x_fmt.format(x), st.SMALL_SIZE, st.MUTED)
                   .next_to(ax.c2p(x, y_range[0]), DOWN, buff=0.12))
    for y in y_ticks:
        labels.add(st.T(y_fmt.format(y), st.SMALL_SIZE, st.MUTED)
                   .next_to(ax.c2p(x_range[0], y), LEFT, buff=0.12))
    ticks = VGroup(*[Line(ax.c2p(x, y_range[0]), ax.c2p(x, y_range[0]) + 0.08 * DOWN,
                          color=st.GRID, stroke_width=2) for x in x_ticks],
                   *[Line(ax.c2p(x_range[0], y), ax.c2p(x_range[0], y) + 0.08 * LEFT,
                          color=st.GRID, stroke_width=2) for y in y_ticks])
    return VGroup(ax, ticks, labels)


def global_frame(d: dict) -> dict:
    """Re-express the tracked pose, its shift and the dots in ONE fixed frame.

    The loss is unchanged by a global rotation of the source compensated in every
    pose, so estimates are only defined up to that gauge. The data dict superposes
    every record on the reference separately (a fresh superposition per step), which
    makes poses appear to move during source steps. Here a single superposition (the
    one of the final source onto the reference) is applied to every record, so a
    source step leaves the displayed pose and shift exactly unchanged.
    """
    X_true = d["X_true"]
    Qf, mef, mtf = ne._kabsch(d["hist_X_raw"][-1], X_true)
    R_g, t_g, X_g = [], [], []
    for k in range(len(d["step"])):
        Qk, mek, mtk = ne._kabsch(d["hist_X_raw"][k], X_true)
        R_al = d["track_R"][k]
        R_raw = R_al @ Qk                                   # optimiser's own pose
        t_raw = d["track_t"][k] - (R_raw @ mek - R_al @ mtk)[:2]
        R_new = R_raw @ Qf.T
        R_g.append(R_new)
        t_g.append(t_raw + (R_raw @ mef - R_new @ mtf)[:2])
        X_g.append((d["hist_X_raw"][k] - mef) @ Qf.T + mtf)
    return dict(R=np.stack(R_g), t=np.stack(t_g), X=np.stack(X_g))


def label(text: str, tex: str | None = None, color=st.TEXT, size=28, tex_size=38):
    parts = [st.T(text, size, color)]
    if tex:
        parts.append(st.M(tex, size=tex_size, color=color))
    return VGroup(*parts).arrange(RIGHT, buff=0.15)


class Scene04IterativeRefinement(Scene):
    def construct(self):
        d = cached("refinement", ne.refinement)
        gf = global_frame(d)                      # one fixed superposition for all steps
        track = int(d["track"])
        steps = d["step"].astype(int)            # 0 init, 1 after pose step, 2 after source step
        n_rec = len(steps)
        n_iter = (n_rec - 1) // 2
        n_obs = d["imgs"].shape[0]
        n_pix = d["imgs"][0].size
        it = np.arange(n_rec) / 2.0               # iteration coordinate: pose steps at k+0.5

        tag = st.illustration_label()
        self.add(tag)
        title = st.title_card("Iterative refinement",
                              "What does iterative reconstruction optimize?")
        self.play(FadeIn(title), run_time=1.0)

        # ============================================================ beat 1
        H = 2.2
        obs = d["imgs"][track]
        lo, hi = float(obs.min()), float(obs.max())
        res_max = float(np.abs(d["track_res"]).max())

        def obs_img():
            return st.heatmap(obs, H, vmin=lo, vmax=hi)

        def pred_img(k):
            return st.heatmap(d["track_pred"][k], H, vmin=lo, vmax=hi)

        def res_img(k):
            # dark-centred signed map: zero residual = panel colour, not a bright square
            return st.heatmap(d["track_res"][k], H, vmax=res_max, symmetric=True,
                              stops=[(0.0, "#3B6FD8"), (0.5, st.PANEL), (1.0, st.RESID)])

        y_img = 0.75
        x_obs, x_pred, x_res = -5.45, -2.95, -0.45
        im_obs = obs_img().move_to([x_obs, y_img, 0])
        fr_obs = st.framed(im_obs, st.OBS)
        lab_obs = label("observation", "I_i", st.OBS).next_to(fr_obs, UP, buff=0.18)
        im_pred = pred_img(0).move_to([x_pred, y_img, 0])
        fr_pred = st.framed(im_pred, st.EST)
        lab_pred = label("prediction", r"\widehat I_i", st.EST).next_to(fr_pred, UP, buff=0.18)
        # common text baseline for the two headings (the hat makes the math taller);
        # align ascender tops, since only "prediction" has a descender
        lab_pred.shift((lab_obs[0].get_top()[1] - lab_pred[0].get_top()[1]) * UP)

        self.play(FadeIn(im_obs), Create(fr_obs), FadeIn(lab_obs), run_time=1.0)
        self.play(FadeIn(im_pred), Create(fr_pred), FadeIn(lab_pred), run_time=1.0)

        eq1 = st.M(r"\widehat I_i", "=", r"\mathrm{CTF}_i\,P_{R_i}\rho", "+", r"\text{shift}",
                   size=42)
        eq1.move_to([(x_obs + x_pred) / 2, -1.0, 0])
        note1 = st.simplified_note("schematic: CTF multiplies in Fourier space,\n"
                                   "shift = phase e^(-2πik·t); toy: CTF = 1")
        note1.next_to(eq1, DOWN, buff=0.2).align_to(fr_obs, LEFT)
        self.play(Write(eq1), FadeIn(note1), run_time=1.5)

        # viewing-direction inset (fixed oblique view of the unit vectors)
        Bv = cam(15, 25)
        d_true = pj.viewing_direction(d["R_true"][track])
        d_est = np.array([pj.viewing_direction(R) for R in gf["R"]])
        vd_err = np.degrees(np.arccos(np.clip(d_est @ d_true, -1, 1)))
        o_v = np.array([2.55, 1.1, 0])
        s_v = 1.55
        ax_v = axes_3d(Bv, o_v, 0.45, length=1.0)          # axes drawn at the arrow origin
        vt = ValueTracker(0.0)

        def vec_at(t):
            k0 = int(np.floor(t)); k1 = min(k0 + 1, n_rec - 1); f = t - k0
            v = (1 - f) * d_est[k0] + f * d_est[k1]
            return v / np.linalg.norm(v)

        ref_arrow = Arrow(o_v, o_v + s_v * to3(d_true @ Bv), buff=0, color=st.TRUTH,
                          stroke_width=11, max_tip_length_to_length_ratio=0.2)
        est_arrow = always_redraw(lambda: Arrow(
            o_v, o_v + s_v * to3(vec_at(vt.get_value()) @ Bv), buff=0, color=st.EST,
            stroke_width=5, max_tip_length_to_length_ratio=0.14))
        lab_v = st.T("viewing direction", 28).move_to([o_v[0], 2.2, 0])
        # colour key stacked in the free space left of the arrows (clear of the readout
        # and the residual-norm plot below)
        leg_v = VGroup(st.T("reference", 26, st.TRUTH), st.T("estimate", 26, st.EST)) \
            .arrange(DOWN, aligned_edge=LEFT, buff=0.08).move_to([1.6, 0.1, 0])

        def vd_text():
            # 2-DOF quantity: angle between the two viewing directions
            k = int(round(vt.get_value()))
            return label("direction error", rf"{vd_err[k]:.1f}^\circ", st.TEXT, 28, 36) \
                .move_to([o_v[0], -0.6, 0])
        vd_lab = always_redraw(vd_text)

        self.play(FadeIn(ax_v), FadeIn(lab_v), FadeIn(leg_v), run_time=0.8)
        self.play(GrowArrow(ref_arrow), run_time=0.8)
        self.add(est_arrow)
        self.play(FadeIn(est_arrow), FadeIn(vd_lab), run_time=0.8)

        # in-plane shift inset (zoomed): arrow from true shift to estimated shift
        t_true = d["t_true"][track]
        t_hist = gf["t"]
        sh_err = np.linalg.norm(t_hist - t_true, axis=1)
        o_s = np.array([5.35, 0.65, 0])
        zoom = 4.0                                         # screen units per data unit
        box_s = Square(1.9, color=st.GRID, stroke_width=2).move_to(o_s)
        cross = VGroup(Line(o_s + 0.95 * LEFT, o_s + 0.95 * RIGHT),
                       Line(o_s + 0.95 * DOWN, o_s + 0.95 * UP)).set_stroke(st.GRID, 1.5,
                                                                          opacity=0.6)
        dot_t = Dot(o_s, radius=0.11, color=st.TRUTH)

        def t_at(t):
            k0 = int(np.floor(t)); k1 = min(k0 + 1, n_rec - 1); f = t - k0
            return (1 - f) * t_hist[k0] + f * t_hist[k1]

        def shift_mob():
            p = o_s + zoom * to3(t_at(vt.get_value()) - t_true)
            g = VGroup()
            if np.linalg.norm(p - o_s) > 0.12:
                g.add(Arrow(o_s, p, buff=0, color=st.EST, stroke_width=5,
                            max_tip_length_to_length_ratio=0.3))
            g.add(Dot(p, radius=0.065, color=st.EST))
            return g
        shift_arrow = always_redraw(shift_mob)
        lab_s = st.T("in-plane shift", 28).move_to([o_s[0], 0, 0])
        lab_s.align_to(lab_v, DOWN)                         # same baseline as lab_v
        zoom_note = st.T(f"zoom x{zoom:g}", 26, st.MUTED).next_to(lab_s, DOWN, buff=0.08)

        def sh_text():
            k = int(round(vt.get_value()))
            return label("error", f"{sh_err[k]:.3f}", st.TEXT, 28, 36).next_to(box_s, DOWN,
                                                                              buff=0.18)
        sh_lab = always_redraw(sh_text)
        self.play(Create(box_s), FadeIn(cross), FadeIn(dot_t), FadeIn(lab_s),
                  FadeIn(zoom_note), run_time=0.8)
        self.add(shift_arrow)
        self.play(FadeIn(shift_arrow), FadeIn(sh_lab), run_time=0.6)
        self.wait(0.4)

        # ============================================================ beat 2
        im_res = res_img(0).move_to([x_res, y_img, 0])
        fr_res = st.framed(im_res, st.RESID)
        lab_res = st.M(r"r_i = I_i-\widehat I_i", size=38, color=st.RESID) \
            .next_to(fr_res, UP, buff=0.18)
        # residual-norm plot
        rn = d["track_res_norm"]
        noise_lvl = float(d["sigma"]) * np.sqrt(n_pix)
        y0, y1 = 2.4, 3.9
        rplot = axes_2d([0, n_iter, 1], [y0, y1, 0.5], 4.6, 1.7, range(0, n_iter + 1, 2),
                        [2.5, 3.0, 3.5], y_fmt="{:.1f}")
        rplot.move_to([4.0, -2.13, 0])
        rax = rplot[0]
        rtitle = label("residual norm", r"\|r_i\|", st.RESID)
        rtitle.next_to(rax.c2p(n_iter, y1), DOWN + LEFT, buff=0.12)
        rxlab = st.T("iteration", 26, st.MUTED).next_to(rplot, DOWN, buff=0.05) \
            .align_to(rplot, RIGHT)
        noise_line = DashedLine(rax.c2p(0, noise_lvl), rax.c2p(n_iter, noise_lvl),
                                color=st.MUTED, stroke_width=2.5)
        noise_lab = st.M(rf"\mathbb E\|\epsilon_i\|\approx\sigma\sqrt{{N_\mathrm{{pix}}}}"
                         rf"={noise_lvl:.2f}", size=36,
                         color=st.MUTED).next_to(noise_line, UP, buff=0.22) \
            .align_to(noise_line, RIGHT)          # right end: clear of the early descent
        r_pts = [rax.c2p(it[k], rn[k]) for k in range(n_rec)]
        r_dots = VGroup(*[Dot(p, radius=0.055, color=st.RESID) for p in r_pts])
        r_segs = VGroup(*[Line(r_pts[k - 1], r_pts[k], color=st.RESID, stroke_width=3)
                          for k in range(1, n_rec)])

        # bottom-left: what each step does
        def step_caption(k):
            if steps[k] == 0:
                s = "initial estimate"
            elif steps[k] == 1:
                s = f"iteration {(k + 1) // 2}: pose step (dots fixed)"
            else:
                s = f"iteration {k // 2}: source step (poses fixed)"
            return st.T(s, 28, st.MUTED).to_edge(LEFT, buff=st.MARGIN).set_y(-3.25)

        cap = step_caption(0)
        # gauge: the loss cannot fix a global rotation of the source (compensated in
        # every pose); the errors here use one fixed superposition on the reference
        gauge = st.T("errors: estimate superposed on reference\n"
                     "(the loss cannot fix a global rotation)", 26, st.MUTED)
        gauge.to_edge(LEFT, buff=st.MARGIN).set_y(-2.6)
        self.play(FadeIn(im_res), Create(fr_res), FadeIn(lab_res), run_time=1.0)
        self.play(FadeIn(rplot), FadeIn(rtitle), FadeIn(rxlab), Create(noise_line),
                  FadeIn(noise_lab), FadeIn(r_dots[0]), FadeIn(cap), FadeIn(gauge),
                  run_time=1.0)

        for k in range(1, n_rec):
            new_pred = pred_img(k).move_to(im_pred)
            new_res = res_img(k).move_to(im_res)
            self.remove(im_pred, im_res)
            self.add(new_pred, new_res)
            self.bring_to_front(fr_pred, fr_res)
            im_pred, im_res = new_pred, new_res
            # name the first two steps; the fast later steps share one caption
            if k <= 2:
                new_cap = step_caption(k)
            elif k == 3:
                new_cap = st.T(f"iterations 2–{n_iter}: pose and source steps alternate", 28,
                               st.MUTED).to_edge(LEFT, buff=st.MARGIN).set_y(-3.25)
            else:
                new_cap = cap
            if new_cap is not cap:
                self.remove(cap); self.add(new_cap); cap = new_cap
            rt = 1.3 if k <= 2 else 0.45
            self.play(vt.animate.set_value(k), Create(r_segs[k - 1]), FadeIn(r_dots[k]),
                      run_time=rt, rate_func=smooth)
            if k <= 2:
                self.wait(0.5)                  # let the named step caption be read
        self.wait(0.6)

        # ============================================================ beat 3
        est_arrow.clear_updaters(); shift_arrow.clear_updaters()
        vd_lab.clear_updaters(); sh_lab.clear_updaters()
        beat12 = Group(im_obs, fr_obs, lab_obs, im_pred, fr_pred, lab_pred, eq1, note1,
                       ax_v, ref_arrow, est_arrow, lab_v, leg_v, vd_lab, box_s, cross, dot_t,
                       shift_arrow, lab_s, zoom_note, sh_lab, im_res, fr_res, lab_res, rplot,
                       rtitle, rxlab, noise_line, noise_lab, r_dots, r_segs, cap, gauge)
        self.play(FadeOut(beat12), run_time=0.8)

        lam = float(d["lam"])
        eqL = st.M(r"\mathcal L(\rho,\{R_i,\mathbf t_i\})", "=",
                   r"\textstyle\sum_i\|I_i-\widehat I_i(\rho,R_i,\mathbf t_i)\|_2^2", "+",
                   r"\lambda\,\mathcal R(\rho)", size=42)
        eqL.move_to([0, 1.95, 0])
        noteL = st.simplified_note("illustrative regularized loss, not the exact objective "
                                   "of a specific program")
        noteL2 = st.M(r"\text{this toy: }\ \mathrm{CTF}=1,\quad \mathcal R(\rho)="
                      r"\textstyle\sum_j|\mathbf x_j|^2,\quad \lambda=" + f"{lam:g}",
                      size=36, color=st.MUTED)
        notes = VGroup(noteL, noteL2).arrange(DOWN, buff=0.12).move_to([0, -3.08, 0])
        self.play(Write(eqL), FadeIn(notes), run_time=1.5)

        PY = -0.45            # panel centre height
        PH = 3.6
        # panel A: full (3-DOF) rotation errors of all observations (bars)
        pe = d["pose_err"]
        pA = RoundedRectangle(corner_radius=0.1, width=3.9, height=PH, color=st.GRID,
                              stroke_width=2).move_to([-4.6, PY, 0])
        BW = 2.8
        bar_ax = axes_2d([0, n_obs, 6], [0, 40, 10], BW, 1.75, [], [0, 20, 40],
                         y_fmt="{:g}°")
        bar_ax.move_to(pA).shift(0.1 * RIGHT)
        bax = bar_ax[0]
        bar_ax.shift((PY - 1.1 - bax.c2p(0, 0)[1]) * UP)   # base 1.1 below panel centre
        bw = BW / n_obs * 0.7

        def bars(k):
            return VGroup(*[Rectangle(width=bw, height=max(bax.c2p(0, pe[k, i])[1]
                                                            - bax.c2p(0, 0)[1], 1e-3),
                                      fill_color=st.EST, fill_opacity=0.9, stroke_width=0)
                            .move_to(bax.c2p(i + 0.5, 0), aligned_edge=DOWN)
                            for i in range(n_obs)])
        bar_g = bars(0)
        tA = st.T(f"pose step: {n_obs} poses", 28)
        st.fit_width(tA, 3.6)
        tA.next_to(pA.get_top(), DOWN, buff=0.15)
        subA = st.T("rotation error (3 DOF)", 26, st.MUTED).next_to(tA, DOWN, buff=0.08)

        def meanA(k):
            return label("mean", rf"{pe[k].mean():.1f}^\circ", st.TEXT, 28, 36) \
                .next_to(pA.get_bottom(), UP, buff=0.15)
        mA = meanA(0)

        # panel B: shared source (five dots, fixed oblique view)
        # same per-record superposition as the pose bars (pose_err), so panels A and B
        # share one reference frame at every step
        hX = d["hist_X"]; Xt = d["X_true"]
        dot_err = np.linalg.norm(hX - Xt, axis=2).mean(1)
        pB = RoundedRectangle(corner_radius=0.1, width=3.9, height=PH, color=st.GRID,
                              stroke_width=2).move_to([-0.45, PY, 0])
        Bd = cam(22, -58)
        o_d = pB.get_center() + 0.6 * DOWN
        sd = 0.82
        ax_d = axes_3d(Bd, o_d, sd, length=0.8).set_opacity(0.6)
        true_dots = VGroup(*[Dot(o_d + sd * to3(x @ Bd), radius=0.12, color=st.TRUTH)
                             for x in Xt])
        est_dots = VGroup(*[Dot(o_d + sd * to3(x @ Bd), radius=0.07, color=st.EST)
                            for x in hX[0]])
        tB = st.T("source step: 5 dots", 28).next_to(pB.get_top(), DOWN, buff=0.15)
        legB = VGroup(st.T("true", 26, st.TRUTH), st.T("estimate (superposed)", 26, st.EST)) \
            .arrange(RIGHT, buff=0.25)
        if legB.width > 3.6:                       # keep size 26: stack instead of shrinking
            legB.arrange(DOWN, buff=0.06)
        legB.next_to(tB, DOWN, buff=0.08)

        def meanB(k):
            return label("mean error", f"{dot_err[k]:.3f}", st.TEXT, 28, 36) \
                .next_to(pB.get_bottom(), UP, buff=0.15)
        mB = meanB(0)

        # panel C: loss vs iteration
        L = d["loss"]; Ld = d["loss_data"]; nf = float(d["noise_floor"])
        lplot = axes_2d([0, n_iter, 1], [60, 420, 50], 4.0, 2.55, range(0, n_iter + 1, 2),
                        [100, 200, 300, 400])
        lplot.move_to([4.25, PY - 0.25, 0])
        lax = lplot[0]
        lxlab = st.T("iteration", 26, st.MUTED).next_to(lplot, DOWN, buff=0.02) \
            .align_to(lplot, RIGHT)
        # expected value of the red data term at the true solution (n_obs * N_pix * sigma^2)
        nf_line = DashedLine(lax.c2p(0, nf), lax.c2p(n_iter, nf), color=st.MUTED,
                             stroke_width=2.5)
        nf_lab = st.M(r"\mathbb E\textstyle\sum_i\|\epsilon_i\|^2", size=36, color=st.MUTED) \
            .next_to(nf_line, DOWN, buff=0.16).align_to(nf_line, LEFT).shift(0.1 * RIGHT)
        leg_L = VGroup(label("", r"\mathcal L", st.TEXT, 26, 36)[1],
                       label("", r"\textstyle\sum_i\|r_i\|^2", st.RESID, 26, 36)[1]) \
            .arrange(RIGHT, buff=0.35)
        leg_L.next_to(lax.c2p(n_iter, 420), DOWN + LEFT, buff=0.1)
        mk_legend = VGroup(
            VGroup(Circle(0.07, color=st.TEXT, stroke_width=2.5),
                   st.T("after pose step", 26, st.MUTED)).arrange(RIGHT, buff=0.1),
            VGroup(Dot(radius=0.07, color=st.TEXT), st.T("after source step", 26, st.MUTED))
            .arrange(RIGHT, buff=0.1)).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        mk_legend.next_to(leg_L, DOWN, buff=0.12).align_to(leg_L, RIGHT)
        # the initial record (iteration 0) has its own square symbol, labelled directly
        # next to the two squares (a third legend row would run into the curves)
        mk_init = st.T("initial", 26, st.MUTED) \
            .next_to(lax.c2p(0, (L[0] + Ld[0]) / 2), RIGHT, buff=0.22)
        mk_legend.add(mk_init)
        tC = st.T("loss (computed)", 28).next_to(lplot, UP, buff=0.12).align_to(lplot, LEFT)
        tC.shift(0.4 * RIGHT)

        def marker(k, y, col):
            p = lax.c2p(it[k], y)
            if steps[k] == 0:                       # initial record: its own symbol
                return Square(0.13, color=col, stroke_width=2.5).move_to(p).set_fill(col, 1)
            if steps[k] == 1:
                return Circle(0.07, color=col, stroke_width=2.5).move_to(p) \
                    .set_fill(st.BG, 1)
            return Dot(p, radius=0.07, color=col)
        L_pts = [lax.c2p(it[k], L[k]) for k in range(n_rec)]
        D_pts = [lax.c2p(it[k], Ld[k]) for k in range(n_rec)]
        L_segs = [Line(L_pts[k - 1], L_pts[k], color=st.TEXT, stroke_width=3)
                  for k in range(1, n_rec)]
        D_segs = [Line(D_pts[k - 1], D_pts[k], color=st.RESID, stroke_width=2.5)
                  for k in range(1, n_rec)]
        L_mk = [marker(k, L[k], st.TEXT) for k in range(n_rec)]
        D_mk = [marker(k, Ld[k], st.RESID) for k in range(n_rec)]

        self.play(Create(pA), Create(pB), FadeIn(bar_ax), FadeIn(bar_g), FadeIn(tA),
                  FadeIn(subA), FadeIn(mA), FadeIn(ax_d), FadeIn(true_dots),
                  FadeIn(est_dots), FadeIn(tB), FadeIn(legB), FadeIn(mB),
                  FadeIn(lplot), FadeIn(lxlab), Create(nf_line), FadeIn(nf_lab),
                  FadeIn(leg_L), FadeIn(mk_legend), FadeIn(tC),
                  FadeIn(L_mk[0]), FadeIn(D_mk[0]), run_time=1.2)
        self.wait(0.6)

        on, off = dict(color=st.TEXT, width=4), dict(color=st.GRID, width=2)
        for k in range(1, n_rec):
            anims = [Create(L_segs[k - 1]), Create(D_segs[k - 1]), FadeIn(L_mk[k]),
                     FadeIn(D_mk[k])]
            if steps[k] == 1:
                anims += [Transform(bar_g, bars(k)), pA.animate.set_stroke(**on),
                          pB.animate.set_stroke(**off)]
                new_m = meanA(k); self.remove(mA); self.add(new_m); mA = new_m
            else:
                anims += [est_dots[j].animate.move_to(o_d + sd * to3(hX[k][j] @ Bd))
                          for j in range(len(Xt))]
                anims += [pB.animate.set_stroke(**on), pA.animate.set_stroke(**off)]
                new_m = meanB(k); self.remove(mB); self.add(new_m); mB = new_m
            self.play(*anims, run_time=0.8 if k <= 2 else 0.5)
        # the one observation that stays far off: the local pose search is non-convex
        i7 = int(np.argmax(pe[-1]))
        lm = st.T("local\nminimum", 26, st.MUTED) \
            .next_to(bax.c2p(i7 + 0.5, pe[-1, i7]), RIGHT, buff=0.15)
        lim_l, lim_r = pA.get_left()[0] + 0.12, pA.get_right()[0] - 0.12
        if lm.get_right()[0] > lim_r:
            lm.shift((lim_r - lm.get_right()[0]) * RIGHT)
        if lm.get_left()[0] < lim_l:
            lm.shift((lim_l - lm.get_left()[0]) * RIGHT)
        self.play(pA.animate.set_stroke(**off), pB.animate.set_stroke(**off), FadeIn(lm),
                  run_time=0.4)
        self.wait(1.2)

        # ============================================================ beat 4
        beat3 = VGroup(pA, pB, bar_ax, bar_g, tA, subA, lm, mA, ax_d, true_dots, est_dots, tB,
                       legB, mB, lplot, lxlab, nf_line, nf_lab, leg_L, mk_legend, tC,
                       *L_segs, *D_segs, *L_mk, *D_mk, notes)
        self.play(FadeOut(beat3), run_time=0.7)
        pen_box = SurroundingRectangle(eqL[4], color=st.LABEL, buff=0.08, stroke_width=3)
        self.play(Create(pen_box), run_time=0.6)

        xs_d, f_d = d["reg_x"], d["reg_f"]
        xs, ys = d["reg_xs"], d["reg_y"]

        def fit_panel(fit, head, cx, data_res, true_err, pen):
            pl = axes_2d([-2, 2, 1], [-0.15, 0.45, 0.15], 4.0, 1.8, [-2, 0, 2], [0, 0.3])
            pl.move_to([cx, -0.25, 0])
            ax = pl[0]
            truth = ax.plot_line_graph(xs_d, f_d, add_vertex_dots=False,
                                       line_color=st.TRUTH, stroke_width=3)
            dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.05, color=st.OBS)
                            for x, y in zip(xs, ys)])
            fcurve = ax.plot_line_graph(xs_d, fit, add_vertex_dots=False, line_color=st.EST,
                                        stroke_width=4)
            hd = st.T(head, 28).next_to(pl, UP, buff=0.12)
            rows = VGroup(
                label("data misfit", rf"\textstyle\sum(y-\hat f)^2={data_res:.3f}",
                      st.TEXT, 28, 36),
                label("true error", rf"\mathrm{{rms}}(\hat f-f)={true_err:.3f}", st.TEXT,
                      28, 36),
                label("penalty", rf"\textstyle\int \hat f''^{{\,2}}dx={pen:.0f}", st.MUTED,
                      28, 36))
            for r in rows:
                r[0].set_color(st.MUTED)
            rows.arrange(DOWN, aligned_edge=LEFT, buff=0.06).next_to(pl, DOWN, buff=0.18)
            rows.set_x(pl[0].get_center()[0])
            lim = st.FRAME_W / 2 - st.MARGIN - 0.05          # keep the outer margin
            if rows.get_right()[0] > lim:
                rows.shift((lim - rows.get_right()[0]) * RIGHT)
            if rows.get_left()[0] < -lim:
                rows.shift((-lim - rows.get_left()[0]) * RIGHT)
            return pl, truth, dots, fcurve, hd, rows

        P1 = fit_panel(d["reg_fit_ls"], "no penalty: fits the noise", -4.15,
                       float(d["reg_data_res_ls"]), float(d["reg_true_err_ls"]),
                       float(d["reg_penalty_ls"]))
        P2 = fit_panel(d["reg_fit_rg"], "with smoothness penalty", 4.15,
                       float(d["reg_data_res_rg"]), float(d["reg_true_err_rg"]),
                       float(d["reg_penalty_rg"]))
        # the 1D panels use their own penalty (smoothness), not the boxed toy R(rho)
        lam_s = float(d["reg_lam_smooth"])
        e_s = int(np.floor(np.log10(lam_s)))
        lam_tex = rf"{lam_s / 10 ** e_s:g}\times10^{{{e_s}}}"
        legR = VGroup(st.T("1D analogue", 28, st.TEXT),
                      st.T(f"{len(xs)} noisy samples", 26, st.OBS),
                      st.T("true profile f", 26, st.TRUTH), st.T("fit", 26, st.EST),
                      st.M(r"\mathcal R(f)=\textstyle\int f''^2dx", size=36,
                           color=st.MUTED),
                      st.M(r"\lambda=" + lam_tex, size=36, color=st.MUTED)) \
            .arrange(DOWN, buff=0.14)
        # centre the legend in the gap between the two plots
        gap_c = (P1[0].get_right()[0] + P2[0].get_left()[0]) / 2
        legR.move_to([gap_c, -0.05, 0])
        self.play(*[FadeIn(m) for m in (P1[0], P1[2], P1[4], P2[0], P2[2], P2[4], legR)],
                  Create(P1[1]), Create(P2[1]), run_time=1.0)
        self.play(Create(P1[3]), Create(P2[3]), run_time=1.3)
        self.play(FadeIn(P1[5][0]), FadeIn(P2[5][0]), run_time=0.6)
        self.play(FadeIn(P1[5][1]), FadeIn(P2[5][1]), run_time=0.6)
        self.play(FadeIn(P1[5][2]), FadeIn(P2[5][2]), run_time=0.6)
        concl = st.T("no penalty: smaller data misfit but larger true error; "
                     "the penalty favours the stable fit", 28, st.TEXT)
        st.fit_width(concl, 13.3).move_to([0, -3.35, 0])   # 13.3 keeps size >= 26
        self.play(FadeIn(concl), run_time=0.6)
        self.wait(2.0)

        # ============================================================ beat 5
        everything = Group(*[m for m in self.mobjects if m is not tag])
        card = st.takeaway_card(
            ["estimate poses  \u2192  update common source  \u2192  avoid overfitting",
             "a smaller data residual is not automatically a better source",
             "sharpening is a separate post-refinement filter"],
            math=r"\mathcal L=\textstyle\sum_i\|I_i-\widehat I_i\|_2^2+\lambda\,\mathcal R(\rho)")
        card.move_to(ORIGIN + 0.2 * DOWN)
        self.play(FadeOut(everything), run_time=0.6)
        self.play(FadeIn(card), run_time=0.8)
        self.wait(3.6)
