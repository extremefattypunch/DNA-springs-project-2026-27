"""Clip 04 - Iterative refinement: poses, density, residuals, regularisation.

Storyboard beats (spec section 3, Scene 04) -> code sections:
  beat 1  ( 0-8 s)  _beat_mismatch     observed particle, crude map, wrong pose glyph
  beat 2  ( 8-17 s) _beat_pose_steps   pose/shift steps for one particle; residual shrinks
  beat 3  (17-26 s) _beat_loss         alternating pose/map steps, loss curve, before/after map
  beat 4  (26-31 s) _beat_regulariser  illustrative spatially varying regularisation
  beat 5  (31-35 s) _beat_takeaway     "align poses + update density + control overfitting"

All images, curves and numbers come from dm.refinement() (toy projection
matching) plus a scene-local replay of it (_refinement_track) that also records
the tracked particle's pose at every step and the map after every iteration.
The replay is checked against the shared cache (identical loss history).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from scipy import ndimage
from manim import *
from helpers import styling as st, demos as dm, projections as pj, fourier as fo, optics as op
from helpers import synthetic_phantom as sp
from helpers.cache import cached

TRACK = 0          # tracked particle (same as dm.refinement)


# ----------------------------------------------------------------- clip-local data
def _refinement_track(n: int = 48, voxel: float = 2.5, n_particles: int = 120, snr: float = 1.0,
                      init_err_deg: float = 20.0, iters: int = 8, smooth_sigma: float = 0.8,
                      lam: float = 2.0, steps=(8.0, 6.0, 4.0, 3.0, 2.0, 1.5, 1.0, 1.0)) -> dict:
    """Replay of dm.refinement (same seed, same steps) that additionally keeps the
    tracked particle's pose (R, t) at every record and the map after every
    iteration. Returns loss_total too so it can be checked against the cache."""
    rng = np.random.default_rng(dm.SEED + 4)
    R_true = pj.random_rotations(n_particles, rng)
    t_true = np.clip(rng.normal(0, 2.5, (n_particles, 2)), -6, 6)
    clean = np.stack([pj.project_analytic(R, "A", n, voxel, atoms=False, shift=t) / voxel
                      for R, t in zip(R_true, t_true)])
    sigma = op.sigma_for_snr(clean[0], snr)
    images = clean + sigma * rng.standard_normal(clean.shape)
    R_est = np.stack([R @ dm._random_axis_rotation(init_err_deg, rng) for R in R_true])
    t_est = np.zeros_like(t_true)

    def grad_penalty(v):
        return float(sum(np.sum(np.diff(v, axis=a) ** 2) for a in range(3)))

    def update_volume(R_set, t_set):
        unshifted = np.stack([pj.shift_image(im, -t, voxel) for im, t in zip(images, t_set)])
        v = fo.reconstruct_direct_fourier(unshifted, R_set)
        return ndimage.gaussian_filter(v, smooth_sigma)

    def loss_total(vol):
        F3 = fo.fft3c(fo.pad_volume(vol, 2))
        d = 0.0
        for im, R, t in zip(images, R_est, t_est):
            d += float(np.sum((im - pj.shift_image(fo.project_fourier(F3, R, n), t, voxel)) ** 2))
        return d + lam * grad_penalty(vol)

    vol = ndimage.gaussian_filter(update_volume(R_est, t_est), 2.5)
    vols = [vol.copy()]
    all_R, all_t, ltot = [], [], []

    def record():
        all_R.append(R_est.copy())
        all_t.append(t_est.copy())
        ltot.append(loss_total(vol))

    record()
    for it in range(iters):
        F3 = fo.fft3c(fo.pad_volume(vol, 2))
        cands = dm._small_rotations(steps[min(it, len(steps) - 1)])
        for i in range(n_particles):
            best = (np.inf, None, None)
            FI = np.fft.fft2(images[i])
            for dR in cands:
                Rc = R_est[i] @ dR
                pr = fo.project_fourier(F3, Rc, n)
                cc = np.real(np.fft.ifft2(FI * np.conj(np.fft.fft2(pr))))
                ix, iy = np.unravel_index(np.argmax(cc), cc.shape)
                sx = ((ix + n // 2) % n - n // 2) * voxel
                sy = ((iy + n // 2) % n - n // 2) * voxel
                if abs(sx) > 8 or abs(sy) > 8:
                    sx, sy = t_est[i]
                r = np.sum((images[i] - pj.shift_image(pr, (sx, sy), voxel)) ** 2)
                if r < best[0]:
                    best = (r, Rc, (sx, sy))
            R_est[i], t_est[i] = best[1], best[2]
        record()
        vol = update_volume(R_est, t_est)
        vols.append(vol.copy())
        record()
    all_R, all_t = np.array(all_R), np.array(all_t)
    err = np.array([[dm._angle_between(a, b) for a, b in zip(Rk, R_true)] for Rk in all_R])
    terr = np.linalg.norm(all_t - t_true[None], axis=-1)
    return dict(all_R=all_R, all_t=all_t, err=err, terr=terr, R_true=R_true, t_true=t_true,
                images=images, loss_total=np.array(ltot), vols=np.array(vols, dtype=np.float32))


def _renders(level_frac: float = 0.3, n: int = 160) -> dict:
    """Shaded isodensity views (display direction) of the map after each
    iteration and of the 48^3 truth, each contoured at level_frac * its max."""
    tr = cached("s04_refinement_track", _refinement_track)
    ref = cached("refinement", dm.refinement)
    R = dm.display_R()
    out = {}
    maps = list(tr["vols"]) + [ref["truth"]]
    for k, v in enumerate(maps):
        fn = pj.voxel_density_fn(np.asarray(v, float), float(ref["voxel"]))
        r = pj.render_density(R, state=None, n=n, nz=n, extent=60.0,
                              level=level_frac * float(np.max(v)), density_fn=fn)
        out[f"shade{k}"] = r["shade"]
        out[f"mask{k}"] = r["mask"]
        out[f"label{k}"] = r["label"]
    out["n_maps"] = len(maps)
    return out


def _example_particle(min_improvements: int = 3, min_glyph_deg: float = 15.0,
                      min_glyph_len: float = 0.7) -> dict:
    """Pick one representative particle for the close-up and compute its model
    projection and residual at every record (k = 0 start, odd = after pose step,
    even = after map step).

    Rule: among particles whose final pose error is within 1 deg of the median,
    whose shift error went down, and whose starting pose error is VISIBLE in the
    2D glyph (in-plane angle between the projected body axes of the estimated
    and true pose >= min_glyph_deg, projected true axis >= min_glyph_len of its
    full length), take the one with the most pose-step improvements (ties ->
    closest to the median). Particle 0 (dm's tracked particle) is recomputed
    the same way as a consistency check."""
    tr = cached("s04_refinement_track", _refinement_track)
    ref = cached("refinement", dm.refinement)
    n, voxel = tr["images"].shape[-1], float(ref["voxel"])
    e, terr = tr["err"], tr["terr"]
    fin, med = e[-1], np.median(e[-1])
    pose_rows = e[[0] + list(range(1, e.shape[0], 2))]
    n_impr = (np.diff(pose_rows, axis=0) < -0.3).sum(0)
    a = _body_axis()

    def inplane_deg(R, Rt):
        u, w = (R @ a)[:2], (Rt @ a)[:2]
        c = u @ w / (np.linalg.norm(u) * np.linalg.norm(w))
        return float(np.degrees(np.arccos(np.clip(c, -1, 1))))
    a0 = np.array([inplane_deg(tr["all_R"][0, i], tr["R_true"][i]) for i in range(e.shape[1])])
    L = np.array([np.linalg.norm((R @ a)[:2]) for R in tr["R_true"]])
    ok = ((np.abs(fin - med) < 1.0) & (terr[-1] < terr[0]) & (n_impr >= min_improvements)
          & (a0 >= min_glyph_deg) & (L >= min_glyph_len))
    score = np.where(ok, n_impr - 0.01 * np.abs(fin - med), -np.inf)
    p = int(np.argmax(score))
    F3s = [fo.fft3c(fo.pad_volume(np.asarray(v, float), 2)) for v in tr["vols"]]

    def model(i, k):
        return pj.shift_image(fo.project_fourier(F3s[k // 2], tr["all_R"][k, i], n),
                              tr["all_t"][k, i], voxel)

    K = e.shape[0]
    proj = np.stack([model(p, k) for k in range(K)])
    check0 = np.max(np.abs(np.stack([model(0, k) for k in range(K)]) - ref["track_proj"]))
    return dict(p=p, image=tr["images"][p], proj=proj, resid=tr["images"][p][None] - proj,
                err=e[:, p], terr=terr[:, p], R=tr["all_R"][:, p], t=tr["all_t"][:, p],
                R_true=tr["R_true"][p], t_true=tr["t_true"][p], check0=check0,
                median_final=med, glyph_deg0=a0[p])


def _reg_render(level_frac: float = 0.3, n: int = 200) -> dict:
    """Final map rendered at the display view with the illustrative
    regularisation strength reg_strength (dm.refinement) sampled on its surface."""
    ref = cached("refinement", dm.refinement)
    v = np.asarray(ref["vol_final"], float)
    voxel = float(ref["voxel"])
    r = pj.render_density(dm.display_R(), state=None, n=n, nz=n, extent=60.0,
                          level=level_frac * v.max(), density_fn=pj.voxel_density_fn(v, voxel))
    lam = np.full(r["mask"].shape, np.nan)
    lam[r["mask"]] = pj.voxel_density_fn(np.asarray(ref["reg_strength"], float), voxel)(r["hit"][r["mask"]])
    return dict(shade=r["shade"], mask=r["mask"], lam=lam)


def load_data():
    ref = cached("refinement", dm.refinement)
    tr = cached("s04_refinement_track", _refinement_track)
    assert np.allclose(tr["loss_total"], ref["loss_total"], rtol=1e-6), "replay != cached run"
    rd = cached("s04_refinement_renders", _renders)
    ex = cached("s04_refinement_example", _example_particle,     # explicit params -> cache key
                min_improvements=3, min_glyph_deg=15.0, min_glyph_len=0.7)
    assert float(ex["check0"]) < 1e-6 * np.abs(ref["track_proj"]).max(), "particle-0 check failed"
    rg = cached("s04_refinement_reg", _reg_render)
    return ref, tr, rd, ex, rg




# ----------------------------------------------------------------- scene
H_PANEL = 2.1                      # image panel height (scene units)
BOX_A = 48 * 2.5                   # image width in angstrom
FOG = "#6F8DFF"                    # regularisation strength (no reserved meaning elsewhere)
MID = "#E6E6E6"                    # weak end of the lambda colour scale (neutral)


def _body_axis() -> np.ndarray:
    """Pose glyph: unit vector from the scaffold centroid to the appendage."""
    dom = np.mean([b.center for b in sp.scaffold_domains()], axis=0)
    v = sp.appendage("A").center - dom
    return v / np.linalg.norm(v)


def _soft_shade(shade, mask, sigma: float = 1.5) -> np.ndarray:
    """Gaussian-smooth the isodensity shading inside the mask only (normalised
    convolution, so the silhouette edge does not darken). Removes the fine
    stripe texture of the coarse crude-map surface without changing the shape."""
    m = np.asarray(mask, float)
    num = ndimage.gaussian_filter(np.asarray(shade, float) * m, sigma)
    den = ndimage.gaussian_filter(m, sigma)
    return np.where(m > 0, num / np.maximum(den, 1e-6), 0.0)


def _rms(a) -> float:
    return float(np.sqrt(np.mean(np.asarray(a) ** 2)))


def _swap(old, new):
    """Cross-fade two mobjects (no glyph morphing)."""
    return [FadeOut(old), FadeIn(new)]


def _seq_swap(old, new, total: float):
    """Text change at the same spot: fade the old one out, then the new one in
    (never two captions half-visible on top of each other)."""
    return Succession(FadeOut(old, run_time=0.42 * total), FadeIn(new, run_time=0.58 * total))


class Scene04Refinement(Scene):
    def construct(self):
        self.ref, self.tr, self.rd, self.ex, self.rg = load_data()
        self.badge = st.simulation_badge()
        self.add(self.badge)
        self.title = st.title_card(
            "Iterative refinement",
            "How does the software improve an initial 3D map\nand the estimated particle orientations?")
        self.play(FadeIn(self.title, shift=0.1 * DOWN), run_time=1.0)
        self._beat_mismatch()
        self._beat_pose_steps()
        self._beat_loss()
        self._beat_regulariser()
        self._beat_takeaway()

    # ------------------------------------------------------------- builders
    def _render_img(self, k: int, color=st.EST, height=H_PANEL):
        rd = self.rd
        return st.density_image(dict(shade=_soft_shade(rd[f"shade{k}"], rd[f"mask{k}"]),
                                     mask=rd[f"mask{k}"], label=rd[f"label{k}"]),
                                color=color, height=height)

    def _particle_img(self, arr, height=H_PANEL):
        lo, hi = self.img_range
        return st.array_to_image_mobject(arr, height=height, vmin=lo, vmax=hi, nearest=True)

    def _resid_img(self, arr, height=H_PANEL):
        return st.array_to_image_mobject(arr, height=height, vmax=self.resid_max,
                                         symmetric=True, nearest=True)

    def _glyph(self, R, t, panel, color, dashed=False):
        """Projected body axis (60 A) drawn from the particle centre shifted by t."""
        s = panel.height / BOX_A
        v = (R @ self.axis)[:2] * 60.0
        start = panel.get_center() + np.array([t[0] * s, t[1] * s, 0])
        end = start + np.array([v[0] * s, v[1] * s, 0])
        if dashed:
            line = DashedLine(start, end, color=color, stroke_width=7, dash_length=0.07)
            return VGroup(line, Dot(end, radius=0.07, color=color), Dot(start, radius=0.06, color=color))
        return VGroup(Arrow(start, end, buff=0, color=color, stroke_width=4,
                            max_tip_length_to_length_ratio=0.25),
                      Dot(start, radius=0.06, color=color))

    def _angle_mark(self, R, t, panel):
        """Arc at the estimate's start between the estimated and the true
        projected axis (its in-plane angle; the 3D error is in the readout)."""
        s = panel.height / BOX_A
        start = panel.get_center() + np.array([t[0] * s, t[1] * s, 0])
        ve, vt = (R @ self.axis)[:2], (self.ex["R_true"] @ self.axis)[:2]
        a_e, a_t = np.arctan2(ve[1], ve[0]), np.arctan2(vt[1], vt[0])
        d = (a_e - a_t + np.pi) % (2 * np.pi) - np.pi
        d = d if abs(d) > 1e-3 else 1e-3
        # dark translucent wedge + white arc: reads on the grey particle image
        wedge = Sector(radius=0.75, start_angle=a_t, angle=d, arc_center=start,
                       fill_color=st.BG, fill_opacity=0.6, stroke_width=0)
        arc = Arc(radius=0.75, start_angle=a_t, angle=d, arc_center=start,
                  color=st.TEXT, stroke_width=5)
        return VGroup(wedge, arc).set_z_index(3)

    def _pose_readout(self, k):
        ex = self.ex
        g = VGroup(
            st.M(r"\Delta\theta_i^{\,\mathrm{3D}} = %.1f^\circ" % ex["err"][k], size=36, color=st.EST),
            st.M(r"|\Delta\mathbf t_i| = %.1f\,\text{\AA}" % ex["terr"][k], size=36, color=st.EST),
        ).arrange(DOWN, buff=0.12)
        return g.next_to(self.p_obs, DOWN, buff=0.2)

    def _resid_readout(self, k):
        g = VGroup(
            st.M(r"\mathrm{rms}(r_i) = %.2f" % _rms(self.ex["resid"][k]), size=36, color=st.RESID),
            st.M(r"\text{noise } \sigma = %.2f" % float(self.ref["sigma"]), size=36, color=st.MUTED),
        ).arrange(DOWN, buff=0.12)
        return g.next_to(self.p_res, DOWN, buff=0.2)

    @staticmethod
    def _start_caption():
        return st.T("start: crude map, every pose off by 20° (3D)", 28, st.TEXT)\
            .to_edge(DOWN, buff=st.MARGIN)

    @staticmethod
    def _steps_caption():
        """Static two-line caption for the whole of beat 2; the active step is
        highlighted by opacity only (no text swapping)."""
        g = VGroup(st.T("pose step: best rotation + shift per particle (map fixed)", 28, st.EST),
                   st.T("map step: rebuild ρ from all 120 particles + smooth (poses fixed)", 28, st.TEXT))\
            .arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        return g.to_edge(DOWN, buff=st.MARGIN)

    def _counter(self, it):
        return st.T(f"iteration {it} / 8", 32, st.TEXT).next_to(self.badge, DOWN, buff=0.25)\
            .align_to(self.badge, RIGHT)

    # ------------------------------------------------------------- beat 1
    def _beat_mismatch(self):
        """0-8 s: observed particle vs. projection of a crude map at a wrong pose."""
        ex = self.ex
        self.axis = _body_axis()
        img = ex["image"]
        self.img_range = (float(np.percentile(img, 0.5)), float(np.percentile(img, 99.5)))
        self.resid_max = float(np.percentile(np.abs(ex["resid"][0]), 99.5))
        y = 0.3
        xs = [-5.45, -1.95, 1.35, 4.65]

        def panel(mob, x, color, sw=3):
            mob.move_to([x, y, 0])
            return Rectangle(width=H_PANEL, height=H_PANEL, stroke_color=color,
                             stroke_width=sw).move_to(mob).set_z_index(2)   # frame over image
        self.map_img = self._render_img(0)
        self.p_map = panel(self.map_img, xs[0], st.EST, 2)
        self.mod_img = self._particle_img(ex["proj"][0])
        self.p_mod = panel(self.mod_img, xs[1], st.EST)
        self.obs_img = self._particle_img(img)
        self.p_obs = panel(self.obs_img, xs[2], st.OBS)
        self.res_img = self._resid_img(ex["resid"][0])
        self.p_res = panel(self.res_img, xs[3], st.RESID)

        def head(word, tex, color, pnl):
            g = VGroup(st.T(word, 28, color), st.M(tex, size=40, color=color)).arrange(RIGHT, buff=0.15)
            g[1].align_to(g[0], DOWN).shift(0.03 * DOWN)
            return g.next_to(pnl, UP, buff=0.16)
        self.h_map = head("current map", r"\rho", st.EST, self.p_map)
        self.h_mod = head("model", r"\hat I_i", st.EST, self.p_mod)
        self.h_obs = head("observed", r"I_i", st.OBS, self.p_obs)
        self.h_res = head("residual", r"r_i = I_i - \hat I_i", st.RESID, self.p_res)
        self.fwd_arrow = Arrow(self.p_map.get_right(), self.p_mod.get_left(), buff=0.1,
                               color=st.EST, stroke_width=5)
        self.fwd_lab = st.M(r"T_{\mathbf t_i}P_{R_i}", size=36, color=st.EST)\
            .next_to(self.fwd_arrow, UP, buff=0.08)
        self.vs = st.T("vs", 30, st.MUTED).move_to((self.p_mod.get_right() + self.p_obs.get_left()) / 2)

        disclaimer = st.phantom_disclaimer().to_edge(DOWN, buff=st.MARGIN)
        self.play(FadeIn(self.obs_img), Create(self.p_obs), FadeIn(self.h_obs),
                  FadeIn(disclaimer), run_time=1.0)
        self.wait(0.3)
        self.play(FadeIn(self.map_img), Create(self.p_map), FadeIn(self.h_map), run_time=0.8)
        self.play(GrowArrow(self.fwd_arrow), FadeIn(self.fwd_lab),
                  FadeIn(self.mod_img), Create(self.p_mod), FadeIn(self.h_mod),
                  FadeIn(self.vs), run_time=1.0)
        # pose glyphs on the observed image: truth (known only in simulation) vs estimate
        self.g_true = self._glyph(ex["R_true"], ex["t_true"], self.p_obs, st.TRUTH, dashed=True)
        self.g_true.set_z_index(5)         # green dashes drawn over the orange estimate
        self.g_est = self._glyph(ex["R"][0], ex["t"][0], self.p_obs, st.EST).set_z_index(4)
        self.g_ang = self._angle_mark(ex["R"][0], ex["t"][0], self.p_obs)
        self.leg = VGroup(
            VGroup(DashedLine(ORIGIN, 0.5 * RIGHT, color=st.TRUTH, stroke_width=6, dash_length=0.07),
                   st.T("true pose (simulation only)", 28, st.TRUTH)).arrange(RIGHT, buff=0.15),
            VGroup(Line(ORIGIN, 0.5 * RIGHT, color=st.EST, stroke_width=6),
                   VGroup(st.T("estimated pose", 28, st.EST),
                          st.M(r"R_i,\ \mathbf t_i", size=36, color=st.EST)).arrange(RIGHT, buff=0.15))
            .arrange(RIGHT, buff=0.15),
            st.T("glyph = projected body axis, from the\n"
                 "shifted centre; wedge = in-plane angle", 28, st.MUTED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        self.leg.to_edge(LEFT, buff=st.MARGIN).set_y(-1.75)
        self.pose_txt = self._pose_readout(0)
        start_cap = self._start_caption()
        self.play(Create(self.g_true), GrowFromPoint(self.g_est, self.g_est[1].get_center()),
                  Create(self.g_ang), FadeIn(self.leg), FadeIn(self.pose_txt),
                  _seq_swap(disclaimer, start_cap, 1.2), run_time=1.2)
        self.wait(1.6)
        self.cap = start_cap

    # ------------------------------------------------------------- beat 2
    def _beat_pose_steps(self):
        """8-17 s: pose + shift updates for this particle; the residual loses structure."""
        ex = self.ex
        self.res_txt = self._resid_readout(0)
        cap = self.cap
        self.play(FadeIn(self.res_img), Create(self.p_res), FadeIn(self.h_res),
                  FadeIn(self.res_txt), run_time=1.0)
        self.wait(0.2)
        K = ex["err"].shape[0]
        counter = None
        steps = self._steps_caption()
        DIM = 0.35
        for k in range(1, K):
            rt = 0.8 if k <= 3 else 0.3
            new_mod = self._particle_img(ex["proj"][k]).move_to(self.mod_img)
            new_res = self._resid_img(ex["resid"][k]).move_to(self.res_img)
            new_rt = self._resid_readout(k)
            anims = (_swap(self.mod_img, new_mod) + _swap(self.res_img, new_res)
                     + _swap(self.res_txt, new_rt))
            act, idle = (steps[0], steps[1]) if k % 2 == 1 else (steps[1], steps[0])
            if k == 1:         # one caption change: start -> static two-line step caption
                idle.set_opacity(DIM)
                counter = self._counter(1)
                anims += [_seq_swap(cap, steps, rt), FadeIn(counter)]
                cap = steps
            else:              # highlight the active step (opacity only, text stays put)
                anims += [act.animate.set_opacity(1.0), idle.animate.set_opacity(DIM)]
                if k % 2 == 1:     # a new iteration starts with its pose step
                    # same text except the digit -> Transform only morphs the digit
                    anims.append(Transform(counter, self._counter((k + 1) // 2)))
            if k % 2 == 1:     # pose step: glyph, angle arc and pose numbers move
                new_pt = self._pose_readout(k)
                anims += [Transform(self.g_est, self._glyph(ex["R"][k], ex["t"][k], self.p_obs, st.EST)),
                          Transform(self.g_ang, self._angle_mark(ex["R"][k], ex["t"][k], self.p_obs))]
                anims += _swap(self.pose_txt, new_pt)
                self.pose_txt = new_pt
            else:              # map step: the map panel changes
                new_map = self._render_img(k // 2).move_to(self.map_img)
                anims += _swap(self.map_img, new_map)
                self.map_img = new_map
            self.play(*anims, run_time=rt)
            self.mod_img, self.res_img, self.res_txt = new_mod, new_res, new_rt
            if k <= 3:
                self.wait(0.15)
        self.cap = VGroup(cap, counter)
        self.wait(0.8)

    # ------------------------------------------------------------- beat 3
    def _beat_loss(self):
        """17-26 s: alternate pose and map steps over all particles; loss curve."""
        ref = self.ref
        old = Group(self.map_img, self.p_map, self.h_map, self.fwd_arrow, self.fwd_lab, self.mod_img,
                    self.p_mod, self.h_mod, self.vs, self.obs_img, self.p_obs, self.h_obs, self.res_img,
                    self.p_res, self.h_res, self.g_true, self.g_est, self.g_ang, self.leg, self.pose_txt,
                    self.res_txt, self.cap)
        title_only = self.title[0]
        self.play(FadeOut(old), FadeOut(self.title[1]), run_time=0.6)

        loss = st.M(r"\mathcal L(\rho,\{R_i,\mathbf t_i\})", r"=", r"\sum_i",
                    r"\Big\|\,\tilde I_i", r"-", r"\mathrm{CTF}_i\,\widehat{P_{R_i}\rho}\;"
                    r"e^{-2\pi i\,\mathbf k\cdot\mathbf t_i}", r"\Big\|_2^2", r"+",
                    r"\lambda\,\mathcal R(\rho)", size=40)
        loss[3].set_color(st.OBS)
        loss[5].set_color(st.EST)
        loss.next_to(title_only, DOWN, buff=0.22).set_x(0)
        settings = VGroup(
            st.T("toy run: 120 particles,", 28, st.MUTED),
            st.M(r"N_{\mathrm{pix}} = 120\cdot 48^2,\ \ \mathrm{CTF}=1,\ \ "
                 r"\mathcal R(\rho)=\|\nabla\rho\|^2,\ \ \tilde I=\mathrm{FT}[I]",
                 size=36, color=st.MUTED),
        ).arrange(RIGHT, buff=0.15)
        settings[0].align_to(settings[1], DOWN).shift(0.04 * DOWN)
        note = VGroup(
            st.T("pedagogical objective — NOT the CryoSPARC loss", 28, st.TEXT),
            settings,
        ).arrange(DOWN, buff=0.1)
        note.next_to(loss, DOWN, buff=0.12)
        self.play(Write(loss), run_time=1.2)
        self.play(FadeIn(note), run_time=0.4)

        lt = ref["loss_total"] / 1e5
        # the shared cache predates the loss_noise_floor key -> same formula as dm.refinement
        floor = float(ref.get("loss_noise_floor",
                              self.tr["images"].size * float(ref["sigma"]) ** 2)) / 1e5
        ax = Axes(x_range=[0, 8, 1], y_range=[3, 13, 1], x_length=6.0, y_length=2.6,
                  tips=False, axis_config=dict(color=st.MUTED, stroke_width=2,
                                               include_numbers=True, font_size=28),
                  y_axis_config=dict(numbers_to_include=[4, 6, 8, 10, 12]))
        ax.shift([-6.65 + 0.75 - ax.c2p(0, 3)[0], -2.45 - ax.c2p(0, 3)[1], 0])
        xl = st.T("iteration", 28, st.MUTED).next_to(ax.x_axis, DOWN, buff=0.15)
        yl = st.M(r"\mathcal L/10^5", size=36, color=st.MUTED).next_to(ax.c2p(0, 13), UP, buff=0.12)
        fl = DashedLine(ax.c2p(0, floor), ax.c2p(8, floor), color=st.RESID, stroke_width=3)
        fl_row = VGroup(st.T("noise floor", 28, st.RESID),
                        st.M(r"N_{\mathrm{pix}}\,\sigma^2", size=36, color=st.RESID),
                        st.T("(perfect model)", 28, st.RESID)).arrange(RIGHT, buff=0.12)
        fl_row[1].align_to(fl_row[0], DOWN).shift(0.02 * DOWN)
        fl_lab = VGroup(fl_row, st.T("below it = fitting noise", 28, st.RESID))\
            .arrange(DOWN, aligned_edge=RIGHT, buff=0.08)\
            .next_to(fl, DOWN, buff=0.08).align_to(fl, RIGHT)
        leg = VGroup(
            VGroup(Line(ORIGIN, 0.45 * RIGHT, color=st.EST, stroke_width=6),
                   st.T("pose step (map fixed)", 28, st.EST)).arrange(RIGHT, buff=0.12),
            VGroup(Line(ORIGIN, 0.45 * RIGHT, color=st.TEXT, stroke_width=6),
                   st.T("map step (poses fixed)", 28, st.TEXT)).arrange(RIGHT, buff=0.12),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        leg.next_to(ax.c2p(8, 13), DL, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl), FadeIn(leg), run_time=0.6)

        # right: start map, current map, truth + live numbers
        hR = 1.6
        cx = [1.6, 3.85, 6.65 - hR / 2 - 0.05]
        ry = -0.55
        m0 = self._render_img(0, height=hR).move_to([cx[0], ry, 0])
        mk = self._render_img(0, height=hR).move_to([cx[1], ry, 0])
        mt = self._render_img(int(self.rd["n_maps"]) - 1, color=st.TRUTH, height=hR).move_to([cx[2], ry, 0])
        labs = VGroup(st.T("start", 28, st.EST), st.T("current", 28, st.EST), st.T("truth", 28, st.TRUTH))
        for lab, x in zip(labs, cx):
            lab.move_to([x, ry + hR / 2 + 0.25, 0])

        def cc_txt(k, x):
            return st.M(r"\mathrm{CC} = %.2f" % ref["map_corr"][k], size=36, color=st.EST)\
                .move_to([x, ry - hR / 2 - 0.25, 0])
        c0, ck = cc_txt(0, cx[0]), cc_txt(0, cx[1])
        cc_note = st.T("CC = correlation with truth", 26, st.MUTED).move_to([cx[1], -2.0, 0])

        def perr_txt(k):
            g = VGroup(st.T("mean pose error", 28, st.EST),
                       st.T("%.1f°" % ref["pose_err"][k], 32, st.EST)).arrange(RIGHT, buff=0.2)
            g[1].align_to(g[0], DOWN)
            return g.move_to([cx[1], -2.55, 0])
        pe = perr_txt(0)
        pe_note = st.T("(all 120 particles)", 26, st.MUTED).move_to([cx[1], -3.05, 0])
        dot = Dot(ax.c2p(0, lt[0]), radius=0.07, color=st.EST)
        self.play(FadeIn(Group(m0, mk, mt)), FadeIn(labs), FadeIn(c0), FadeIn(ck), FadeIn(cc_note),
                  FadeIn(pe), FadeIn(pe_note), FadeIn(dot), run_time=0.6)

        labels = ref["label"]
        segs = VGroup()
        for k in range(1, len(lt)):
            col = st.EST if labels[k] == 1 else st.TEXT
            seg = Line(ax.c2p((k - 1) / 2, lt[k - 1]), ax.c2p(k / 2, lt[k]), color=col, stroke_width=6)
            segs.add(seg)
            new_pe = perr_txt(k)
            anims = [Create(seg), dot.animate.move_to(ax.c2p(k / 2, lt[k])), *_swap(pe, new_pe)]
            pe = new_pe
            if labels[k] == 2:
                new = self._render_img(k // 2, height=hR).move_to(mk)
                new_ck = cc_txt(k, cx[1])
                anims += _swap(mk, new) + _swap(ck, new_ck)
                mk, ck = new, new_ck
            self.play(*anims, run_time=0.4 if k <= 4 else 0.28)
        self.play(Create(fl), FadeIn(fl_lab), run_time=0.6)
        self.wait(1.6)
        self.loss_eq, self.loss_note = loss, note
        self.shot3 = Group(ax, xl, yl, fl, fl_lab, leg, segs, dot, pe, pe_note, c0, ck, cc_note,
                           m0, mk, mt, labs)

    # ------------------------------------------------------------- beat 4
    def _beat_regulariser(self):
        """26-31 s: illustrative spatially varying regularisation (rationale of NU-refinement)."""
        rg = self.rg
        self.play(FadeOut(self.shot3), FadeOut(self.loss_note), FadeOut(self.loss_eq), run_time=0.6)
        # final map surface coloured by the illustrative lambda(r); solvent (lambda ~ 1) as fog
        lam = np.asarray(rg["lam"], float)
        m = np.asarray(rg["mask"], bool)
        lo, hi = np.nanpercentile(lam[m], [2, 98])
        u = np.clip((np.nan_to_num(lam, nan=hi) - lo) / (hi - lo), 0, 1)
        rgb = st.colorize(u, [(0.0, MID), (1.0, FOG)]).astype(float)
        rgb = np.clip(rgb * _soft_shade(rg["shade"], m)[..., None], 0, 255)
        rgba = np.concatenate([rgb, 255 * m[..., None].astype(float)], -1).astype(np.uint8)
        rgba = np.stack([st.to_display(rgba[..., ch]) for ch in range(4)], -1)
        hS = 3.2
        surf = ImageMobject(rgba)
        surf.height = hS
        surf.move_to([-4.6, -0.45, 0])
        fog = Rectangle(width=hS + 0.2, height=hS + 0.2, stroke_color=FOG, stroke_width=2,
                        fill_color=FOG, fill_opacity=0.28).move_to(surf)
        fog_lab = st.T("solvent: λ largest", 26, FOG)
        fog_lab.move_to(fog.get_corner(UL) + np.array([0.1 + fog_lab.width / 2, -0.1 - fog_lab.height / 2, 0]))
        head = st.M(r"\lambda\,\mathcal R(\rho)\ \longrightarrow\ \int \lambda(\mathbf r)\,"
                    r"|\nabla\rho(\mathbf r)|^2\,d\mathbf r", size=40)
        head.next_to(self.title[0], DOWN, buff=0.25).set_x(0.6)
        surf_plain = self._render_img(int(self.rd["n_maps"]) - 2, height=hS).move_to(surf)
        # colour bar (titled) and the footnote caption below it
        bar = Rectangle(width=1.0, height=0.26, stroke_width=0, fill_opacity=1)
        bar.set_fill(color=[FOG, MID])
        bar_lab = VGroup(st.M(r"\lambda(\mathbf r)", size=36), st.T("weak", 26, MID), bar,
                         st.T("strong", 26, FOG)).arrange(RIGHT, buff=0.15)
        bar_lab.next_to(fog, DOWN, buff=0.15)
        cap = st.T("final map, coloured by\nillustrative λ(r)", 26, st.MUTED)\
            .next_to(bar_lab, DOWN, buff=0.1)
        self.play(FadeIn(surf_plain), Write(head), FadeIn(cap), run_time=0.9)
        pts = VGroup(
            VGroup(Dot(color=FOG), st.T("strong smoothing where the signal is weak", 28, FOG))
            .arrange(RIGHT, buff=0.2),
            VGroup(Dot(color=MID), st.T("light touch in the reliable core", 28, MID))
            .arrange(RIGHT, buff=0.2),
            st.T("illustrates the rationale of non-uniform refinement,\n"
                 "not its objective; here λ(r) comes from the smoothed\n"
                 "map amplitude; the toy run used one uniform λ", 28, st.MUTED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2)
        pts.move_to([0, 0.75, 0]).align_to(fog, LEFT).shift((fog.width + 0.4) * RIGHT)
        fog.set_z_index(-1)                # fog sits behind the map
        self.play(FadeIn(fog), FadeIn(fog_lab), FadeIn(surf), FadeOut(surf_plain), FadeIn(bar_lab),
                  FadeIn(pts[0]), run_time=0.9)
        self.play(FadeIn(pts[1]), FadeIn(pts[2]), run_time=0.6)
        box_txt = VGroup(
            st.T("Regularization ≠ sharpening", 28, st.BADGE, weight="BOLD"),
            st.T("regularization: a prior used inside every\n"
                 "iteration; it changes which map is fitted", 28, st.TEXT),
            st.T("sharpening: post hoc boost of high-|k| amplitudes\n"
                 "of the finished map; poses and fit are unchanged", 28, st.TEXT),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.14)
        box = SurroundingRectangle(box_txt, buff=0.2, color=st.MUTED, stroke_width=1.5,
                                   fill_color=st.PANEL, fill_opacity=0.95, corner_radius=0.1)
        bx = VGroup(box, box_txt).next_to(pts, DOWN, buff=0.35).align_to(pts, LEFT)
        self.play(FadeIn(bx), run_time=0.6)
        self.wait(3.2)
        self.shot4 = Group(surf, fog, fog_lab, cap, head, bar_lab, pts, bx)

    # ------------------------------------------------------------- beat 5
    def _beat_takeaway(self):
        """31-35 s: freeze on the takeaway."""
        ref = self.ref
        self.play(FadeOut(self.shot4), run_time=0.6)
        card = st.takeaway_card(
            r"\min_{\rho,\,\{R_i,\mathbf t_i\}}\ \sum_i \|r_i\|^2 + \lambda\,\mathcal R(\rho)",
            "align poses + update density + control overfitting",
            "like iterative refinement against measured data,\n"
            "except the particle poses and density are being inferred")
        card.to_edge(DOWN, buff=st.MARGIN)
        hR = 1.75
        n_last = int(self.rd["n_maps"]) - 1
        maps = Group(self._render_img(0, height=hR), self._render_img(n_last - 1, height=hR),
                     self._render_img(n_last, color=st.TRUTH, height=hR))
        for mobj, x in zip(maps, [-4.6, -1.3, 2.0]):
            mobj.move_to([x, 0, 0])
        labs = VGroup(
            st.T("start, CC %.2f" % ref["map_corr"][0], 28, st.EST),
            st.T("8 iterations, CC %.2f" % ref["map_corr"][-1], 28, st.EST),
            st.T("truth", 28, st.TRUTH))
        for lab, mobj in zip(labs, maps):
            lab.next_to(mobj, DOWN, buff=0.08)
        arr = Arrow(maps[0].get_right(), maps[1].get_left(), buff=0.1, color=st.EST, stroke_width=5)
        pe = VGroup(st.T("mean pose error", 28, st.EST),
                    st.T("%.1f° → %.1f°" % (ref["pose_err"][0], ref["pose_err"][-1]),
                         32, st.EST)).arrange(DOWN, buff=0.1)
        top = Group(maps, labs, arr)
        top.shift((card.get_top()[1] + 0.3 - top.get_bottom()[1]) * UP)
        arr.put_start_and_end_on(maps[0].get_right() + 0.1 * RIGHT, maps[1].get_left() + 0.1 * LEFT)
        pe.move_to([5.0, maps[0].get_center()[1], 0])
        self.play(FadeIn(maps), FadeIn(labs), GrowArrow(arr), FadeIn(pe), FadeIn(card), run_time=1.0)
        self.wait(4.0)


if __name__ == "__main__":
    ref, tr, rd, ex, rg = load_data()
    print("particle", ex["p"], "err", np.round(ex["err"], 1))
    print("terr", np.round(ex["terr"], 2), "check0", ex["check0"])
    print("rms resid", [round(_rms(r), 3) for r in ex["resid"]], "sigma", ref["sigma"])
