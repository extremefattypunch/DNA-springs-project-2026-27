"""Clip 03 — Orientation, Fourier slice theorem and 3D reconstruction (~65 s).

Storyboard beats (spec §3, Scene 03) -> code sections:
    beat1  0-10 s   orientation: R in SO(3) = viewing arrow v = R^T e_z (2 DOF) + in-plane ring psi (1 DOF)
    beat2 10-19 s   projection P_R rho for 4 poses, thumbnails accumulate at the screen edge
    beat3 19-32 s   Fourier slice: central plane rotates with R; FT panel updates; numerical check
    beat4 32-45 s   coverage: many central planes, uniform vs preferred-orientation (biased) sets
    beat5 45-57 s   idealized direct-Fourier reconstruction with KNOWN poses, N = 2 ... 300, biased
    beat6 57-65 s   poses unknown -> inverse problem; takeaway freeze

Conventions (helpers/projections.py, verified in tests/test_math_helpers.py):
    rho_R(r) = rho(R^T r); beam along z; v = R^T e_z; FT2[P_R rho](k) = FT3[rho](R^T (kx, ky, 0)).
The 3D drawing lives in the OBJECT frame: the phantom is fixed and the pose
group (beam arrow, detector, in-plane ring / central Fourier plane) is the
canonical z-aligned drawing mapped by R^T.

Every displayed image, number and curve comes from cached("reconstruction",
dm.reconstruction) or the two small clip-local caches below.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *

from helpers import demos as dm
from helpers import drawing as dr
from helpers import fourier as fo
from helpers import projections as pj
from helpers import styling as st
from helpers.cache import cached

# camera and on-screen anchor of the 3D origin
PHI, THETA = 68 * DEGREES, -58 * DEGREES
ORIGIN_SCREEN = np.array([-3.3, -0.6, 0.0])
N_DRAWN = 60            # central planes drawn in the coverage beat (uniform and biased)


# ------------------------------------------------------------- clip-local data
def _coverage_stats(n_grid: int = 128, drawn: int = N_DRAWN, counts=(2, 6, 20, 60)) -> dict:
    """Fraction of Fourier-ball voxels |k| <= k_max touched by >= 1 (rounded) central plane.

    The number depends on the grid; n_grid = 128 is the 2x-padded grid that
    fo.reconstruct_direct_fourier actually fills for the 64^3 maps, and it is
    labelled as such on screen.
    """
    d = cached("reconstruction", dm.reconstruction)
    k = np.arange(n_grid) - n_grid // 2
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    ball = (KX ** 2 + KY ** 2 + KZ ** 2) <= (n_grid // 2 - 1) ** 2

    def frac(rots):
        return float((fo.sampling_coverage(rots, n_grid) & ball).sum() / ball.sum())

    return dict(counts=np.array(counts),
                uniform=np.array([frac(d["rotations"][:c]) for c in counts]),
                uniform_drawn=frac(d["rotations"][:drawn]),
                biased_drawn=frac(d["biased_rotations"][:drawn]))


def _slice_check(kmax: int = 40) -> dict:
    """Complex normalised correlation |<P2, S>| / (|P2| |S|) over |k| <= kmax at pose R2.

    P2 = FT2 of the ideal projection; S = central slice of FT3(truth) taken with the
    correct convention R^T (k, 0) and, for contrast, with the wrong one R (k, 0).
    """
    d = cached("reconstruction", dm.reconstruction)
    truth = d["truth"].astype(float)
    n, vox = truth.shape[0], float(d["voxel"])
    F3 = fo.fft3c(fo.pad_volume(truth, 2))
    R = d["show_R"][1]
    p = pj.project_analytic(R, "A", n, vox, atoms=False) / vox
    P2 = fo.fft2c(np.pad(p, n // 2))
    k = np.arange(P2.shape[0]) - P2.shape[0] // 2
    KX, KY = np.meshgrid(k, k, indexing="ij")
    m = KX ** 2 + KY ** 2 <= kmax ** 2

    def cc(M):          # fo.central_slice(F3, M) samples F3 at M^T (k, 0)
        S = fo.central_slice(F3, M)[m]
        return float(abs(np.vdot(P2[m], S)) / (np.linalg.norm(P2[m]) * np.linalg.norm(S)))

    return dict(corr_RT=cc(R), corr_R=cc(R.T))


def _recon_renders(n: int = 200) -> dict:
    """Shaded isosurfaces (display view) of truth, the 5 known-pose maps and the biased map.

    Every map is contoured at the same enclosed volume as the truth at level 0.35,
    so the comparison shows shape, not the overall scale of each map.
    """
    d = cached("reconstruction", dm.reconstruction)
    vox = float(d["voxel"])
    frac = float((d["truth"] > 0.35).mean())
    vols = [d["truth"]] + list(d["recons"]) + [d["biased_recon"]]
    R = dm.display_R()
    out = {}
    for i, v in enumerate(vols):
        v = v.astype(float)
        lev = float(np.quantile(v, 1 - frac))
        r = pj.render_density(R, state="A", n=n, extent=64.0, level=lev, nz=n,
                              density_fn=pj.voxel_density_fn(v, vox))
        out[f"shade{i}"] = r["shade"].astype(np.float32)
        out[f"mask{i}"] = r["mask"]
        out[f"label{i}"] = r["label"]
    return out


def _render_dict(rr, i):
    return dict(shade=rr[f"shade{i}"], mask=rr[f"mask{i}"], label=rr[f"label{i}"])


# ------------------------------------------------------------- geometry helpers
def disc_points(radius: float, n: int = 36) -> np.ndarray:
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.stack([radius * np.cos(a), radius * np.sin(a), np.zeros(n)], -1)


def central_disc(R: np.ndarray, radius: float, opacity: float, color=st.TRUTH,
                 stroke_opacity: float = 0.5) -> Polygon:
    """Central plane {R^T (kx, ky, 0)} drawn as a disc of radius k_max."""
    pts = disc_points(radius) @ R          # rows: (R^T p)^T = p^T R
    return Polygon(*pts, fill_color=color, fill_opacity=opacity, stroke_color=color,
                   stroke_width=0.8, stroke_opacity=stroke_opacity)


def box_edges(h: float, color=st.MUTED, opacity=0.45) -> VGroup:
    c = [np.array([sx, sy, sz]) * h for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
    g = VGroup()
    for i in range(8):
        for j in range(i + 1, 8):
            if np.sum(np.abs(c[i] - c[j]) > 1e-9) == 1:
                g.add(Line(c[i], c[j], stroke_color=color, stroke_width=1.5,
                           stroke_opacity=opacity))
    return g


class Scene03FourierSlice(ThreeDScene):
    # ------------------------------------------------------------- utilities
    def fx(self, *mobs):
        """Mark mobjects as fixed in frame (2D overlay) without adding them."""
        self.camera.add_fixed_in_frame_mobjects(*mobs)
        return mobs[0] if len(mobs) == 1 else mobs

    def world(self, sx: float, sy: float) -> np.ndarray:
        """World point (depth 0) that appears at screen offset (sx, sy) from the 3D origin."""
        return self.rotm.T @ np.array([sx, sy, 0.0])

    def swap(self, old, new, run_time=0.5):
        self.fx(new)
        return AnimationGroup(FadeOut(old), FadeIn(new), run_time=run_time)

    def panel(self, arr, height, frame_color, stops=st.GRAY, label=None, label_size=28,
              label_color=st.TEXT, nearest=False):
        img = st.array_to_image_mobject(arr, height=height, stops=stops, nearest=nearest)
        fr = st.framed(img, frame_color)
        g = Group(img, fr)
        if label is not None:
            lab = label if isinstance(label, Mobject) else st.T(label, label_size, label_color)
            lab.next_to(fr, DOWN, buff=0.15)
            g.add(lab)
        return g

    # ------------------------------------------------------------- scene
    def construct(self):
        d = cached("reconstruction", dm.reconstruction)
        cov = cached("s03_coverage128", _coverage_stats)
        sc = cached("s03_slice_check", _slice_check)
        rr = cached("s03_recon_renders", _recon_renders)

        self.set_camera_orientation(phi=PHI, theta=THETA)
        self.camera.reset_rotation_matrix()        # otherwise still the identity until the first frame
        self.rotm = self.camera.get_rotation_matrix()
        # the 3D origin is drawn at world point O3, which projects to ORIGIN_SCREEN at depth 0
        # (moving camera.frame_center would also shift the fixed-in-frame overlays)
        O3 = self.rotm.T @ ORIGIN_SCREEN
        self.O3 = O3

        # ------------------------------------------------ persistent overlays
        badge = self.fx(st.simulation_badge())
        title = st.title_card("Fourier slices and 3D reconstruction",
                              "How can 2D projections constrain a 3D density\n"
                              "without an indexed diffraction lattice?")
        self.fx(title)
        self.add(badge)
        self.play(FadeIn(title), run_time=1.0)

        # ================================================= beat 1: orientation (0-10 s)
        rot, tilt, psi = ValueTracker(0.0), ValueTracker(0.0), ValueTracker(0.0)

        def Rcur():
            return pj.rotation_from_euler(rot.get_value(), tilt.get_value(), psi.get_value())

        axes = ThreeDAxes(x_range=[-2, 2, 1], y_range=[-2, 2, 1], z_range=[-2, 2, 1],
                          x_length=4, y_length=4, z_length=4, tips=False,
                          axis_config=dict(stroke_color=st.MUTED, stroke_width=2,
                                           include_ticks=False)).shift(O3)
        ax_labels = VGroup(*[st.M(s, size=36, color=st.MUTED).move_to(O3 + p) for s, p in
                             [("x", 2.25 * RIGHT), ("y", 2.3 * UP + 0.35 * LEFT + 0.45 * OUT),
                              ("z", 2.25 * OUT)]])
        phantom = dr.phantom_3d("A", nsig=1.6, color=st.TRUTH, opacity=0.5,
                                scale=0.038, resolution=(12, 8)).shift(O3)
        disclaimer = self.fx(st.phantom_disclaimer().to_corner(DL, buff=st.MARGIN))

        # canonical pose drawing (identity pose: beam along +z), mapped by R^T
        canon = VGroup(
            Arrow3D(start=-2.35 * OUT, end=-1.45 * OUT, color=st.TEXT, thickness=0.03,
                    height=0.32, base_radius=0.1, resolution=8),
            Square(side_length=1.8, stroke_color=st.MUTED, stroke_width=2,
                   fill_color=st.PANEL, fill_opacity=0.35).shift(1.5 * OUT),
            Circle(radius=0.62, stroke_color=st.TEXT, stroke_width=3).shift(1.5 * OUT),
            Line(ORIGIN, 0.62 * RIGHT, stroke_color=st.TEXT, stroke_width=4).shift(1.5 * OUT),
            Dot(0.62 * RIGHT, radius=0.07, color=st.TEXT).shift(1.5 * OUT),
            # R^T turns the canonical drawing by -psi about +z: the arc points clockwise
            Arc(radius=0.8, start_angle=1.55, angle=-1.3, stroke_color=st.TEXT,
                stroke_width=3).add_tip(tip_length=0.16, tip_width=0.16).shift(1.5 * OUT),
        )
        pose = always_redraw(
            lambda: canon.copy().apply_matrix(Rcur().T, about_point=ORIGIN).shift(O3))

        eq_R = st.M(r"R\in SO(3)", r"\;=\;", r"R_z(\psi)\,R_y(\theta)\,R_z(\varphi)", size=40)
        dof1 = VGroup(st.M(r"\mathbf v = R^{\top}\mathbf e_z", size=38),
                      st.T("arrow: viewing direction (φ, θ) = 2 DOF", 28)
                      ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        dof2 = VGroup(st.M(r"\psi", size=40),
                      st.T("ring: in-plane rotation = 1 DOF", 28)
                      ).arrange(RIGHT, buff=0.2)
        cap1 = st.T("Orientation: 3 rotational\ndegrees of freedom", 32, st.TEXT, weight="BOLD")
        right1 = VGroup(eq_R, dof1, dof2, cap1).arrange(DOWN, aligned_edge=LEFT, buff=0.45)
        right1.move_to(np.array([3.5, -0.55, 0]))
        right1.shift(min(0, 6.6 - right1.get_right()[0]) * RIGHT)
        self.fx(right1)

        self.add_fixed_orientation_mobjects(*ax_labels)
        self.remove(*ax_labels)
        self.play(Create(axes), FadeIn(ax_labels), FadeIn(phantom), FadeIn(disclaimer),
                  run_time=1.5)
        self.play(FadeIn(pose), FadeIn(eq_R), run_time=0.8)
        self.play(FadeIn(dof1), tilt.animate.set_value(60), run_time=1.6)
        self.play(rot.animate.set_value(40), run_time=1.4)
        self.play(FadeIn(dof2), psi.animate.set_value(100), run_time=1.6)
        self.play(FadeIn(cap1), run_time=0.6)
        self.wait(1.2)
        # ---- t ~ 9.7 s

        # ================================================= beat 2: projection (10-19 s)
        q = title[1]
        eq_P = self.fx(st.M(r"(P_R\rho)(x,y)=\int \rho\!\left(R^{\top}(x,y,z)^{\top}\right)dz",
                            size=40).move_to(np.array([3.3, 1.9, 0])))
        self.play(FadeOut(right1), FadeOut(q), FadeOut(disclaimer), FadeIn(eq_P), run_time=0.9)

        poses = [tuple(a) for a in [(0, 0, 0), (40, 70, 20), (-60, 110, 75), (120, 35, -40)]]
        TS = 1.2          # thumbnail size; 2x2 grid at the right edge
        grid = [np.array([4.5 + 1.55 * (i % 2), 0.55 - 2.0 * (i // 2), 0]) for i in range(4)]
        # fill the right column first, so a moving image never passes over a placed thumbnail
        slots = [grid[j] for j in (1, 3, 0, 2)]
        thumbs = []
        for i, (a, b, c) in enumerate(poses):
            self.play(rot.animate.set_value(a), tilt.animate.set_value(b),
                      psi.animate.set_value(c), run_time=1.0)
            big = self.panel(d["show_proj"][i], 2.2, st.TRUTH).move_to(np.array([1.75, -0.35, 0]))
            lab = st.M(rf"P_{{R_{i + 1}}}\rho", size=36).next_to(big, DOWN, buff=0.12)
            self.fx(big, lab)
            self.play(FadeIn(big), FadeIn(lab), run_time=0.45)
            slab = st.M(rf"P_{{R_{i + 1}}}\rho", size=36).move_to(
                slots[i] + (TS / 2 + 0.28) * DOWN)
            self.fx(slab)
            self.play(big.animate.scale(TS / 2.2).move_to(slots[i]), ReplacementTransform(lab, slab),
                      run_time=0.55)
            thumbs.append(Group(big, slab))
        # ---- t ~ 18.6 s

        # ================================================= beat 3: Fourier slice (19-32 s)
        H = 1.55
        kbox = box_edges(H)
        kaxes = VGroup(*[Line(-1.95 * v, 1.95 * v, stroke_color=st.MUTED, stroke_width=1.5)
                         for v in (RIGHT, UP, OUT)]).shift(O3)
        kbox.shift(O3)
        k_labels = VGroup(*[st.M(s, size=36, color=st.MUTED).move_to(O3 + p) for s, p in
                            [("k_x", 2.25 * RIGHT), ("k_y", 2.7 * UP),
                             ("k_z", 2.3 * OUT + 0.3 * RIGHT)]])
        canonF = VGroup(
            Polygon(*disc_points(H), fill_color=st.TRUTH, fill_opacity=0.45,
                    stroke_color=st.TRUTH, stroke_width=2),
            Arrow3D(start=ORIGIN, end=1.6 * OUT, color=st.TEXT, thickness=0.025, height=0.28,
                    base_radius=0.09, resolution=8),
        )
        planeF = always_redraw(
            lambda: canonF.copy().apply_matrix(Rcur().T, about_point=ORIGIN).shift(O3))
        fs_head = self.fx(st.T("3D Fourier space  ρ̂(k)", 30, st.TEXT)
                          .move_to(np.array([-3.3, 2.05, 0])))
        self.add_fixed_orientation_mobjects(*k_labels)
        self.remove(*k_labels)
        self.play(FadeOut(phantom), FadeOut(pose), FadeOut(disclaimer), FadeOut(axes),
                  FadeOut(ax_labels), FadeOut(eq_P), *[FadeOut(t) for t in thumbs],
                  run_time=0.9)
        # start the Fourier view at pose R1
        rot.set_value(0), tilt.set_value(0), psi.set_value(0)
        self.play(Create(kbox), Create(kaxes), FadeIn(k_labels), FadeIn(planeF), FadeIn(fs_head),
                  run_time=1.2)

        eq_F = self.fx(st.M(r"\widehat{P_R\rho}(k_x,k_y)", r"=",
                            r"\hat\rho\!\left(R^{\top}(k_x,k_y,0)^{\top}\right)", size=40)
                       .move_to(np.array([3.35, 1.85, 0])))
        self.play(Write(eq_F), run_time=1.3)

        PX, FX, PY, PH = 2.1, 4.85, 0.05, 1.95

        def pair(i):
            p = self.panel(d["show_proj"][i], PH, st.TRUTH).move_to(np.array([PX, PY, 0]))
            f = self.panel(d["show_ft"][i], PH, st.TRUTH, stops=st.MAGMA).move_to(
                np.array([FX, PY, 0]))
            return Group(p, f)

        lab_p = st.M(r"P_R\rho", size=36).move_to(np.array([PX, PY - PH / 2 - 0.35, 0]))
        lab_f = st.M(r"\log\left|\widehat{P_R\rho}\right|", size=36).move_to(
            np.array([FX, PY - PH / 2 - 0.4, 0]))
        arrow_ft = Arrow(np.array([PX + PH / 2 + 0.05, PY, 0]), np.array([FX - PH / 2 - 0.05, PY, 0]),
                         buff=0.0, stroke_width=4, color=st.MUTED,
                         max_tip_length_to_length_ratio=0.35)
        ft2 = st.T("FT₂", 28, st.MUTED).next_to(arrow_ft, UP, buff=0.08)

        def pose_lab(i):
            return st.M(rf"R=R_{i + 1}", size=36).move_to(
                np.array([(PX + FX) / 2, PY - PH / 2 - 1.0, 0]))

        cur = self.fx(pair(0))
        plab = self.fx(pose_lab(0))
        self.fx(lab_p, lab_f, arrow_ft, ft2)
        self.play(FadeIn(cur), FadeIn(lab_p), FadeIn(lab_f), GrowArrow(arrow_ft), FadeIn(ft2),
                  FadeIn(plab), run_time=0.8)
        for i in (2, 3, 1):
            a, b, c = poses[i]
            self.play(rot.animate.set_value(a), tilt.animate.set_value(b),
                      psi.animate.set_value(c), run_time=1.3)
            new, nplab = pair(i), pose_lab(i)
            self.play(self.swap(cur, new, 0.45), self.swap(plab, nplab, 0.45))
            cur, plab = new, nplab
        # ---- t ~ 26.6 s : numerical check at pose R2
        c0, c1 = 32, 96
        CY = PY + 0.3          # lift the check panels so the captions and the numbers do not crowd
        chk_a = self.panel(d["slice_ft2"][c0:c1, c0:c1], PH, st.TRUTH, stops=st.MAGMA)
        chk_b = self.panel(d["slice_central"][c0:c1, c0:c1], PH, st.TRUTH, stops=st.MAGMA)
        chk_a.move_to(np.array([PX, CY, 0]))
        chk_b.move_to(np.array([FX, CY, 0]))
        chk_la = st.T("FT₂ of the\nprojection", 28).next_to(chk_a, DOWN, buff=0.15)
        chk_lb = st.T("central plane\nof the 3D FT", 28).next_to(chk_b, DOWN, buff=0.15)
        eqsign = st.M(r"\approx", size=48).move_to(np.array([(PX + FX) / 2, CY, 0]))
        chk_r = VGroup(
            st.T("complex corr. at R₂:", 28, st.TEXT),
            st.M(rf"R^{{\top}}\!:\ {sc['corr_RT']:.4f}", size=36),
            st.T(f"(R: {sc['corr_R']:.2f})", 28, st.MUTED),
        ).arrange(RIGHT, buff=0.18)
        chk_r.move_to(np.array([3.2, -2.05, 0]))
        chk_r.shift(min(0, 6.6 - chk_r.get_right()[0]) * RIGHT)
        note = st.T("The complex FT (amplitude and phase)\nlies on the plane; showing log|·| is\n"
                    "only a display convention.", 28, st.MUTED)
        note.to_edge(DOWN, buff=st.MARGIN).align_to(chk_r, LEFT)
        note.shift(min(0, 6.6 - note.get_right()[0]) * RIGHT)
        self.fx(chk_a, chk_b, chk_la, chk_lb, eqsign, chk_r, note)
        self.play(FadeOut(cur), FadeOut(lab_p), FadeOut(lab_f), FadeOut(arrow_ft), FadeOut(ft2),
                  FadeOut(plab), FadeIn(chk_a), FadeIn(chk_b), FadeIn(chk_la), FadeIn(chk_lb),
                  FadeIn(eqsign), run_time=0.8)
        self.play(FadeIn(chk_r), FadeIn(note), run_time=0.7)
        self.wait(3.0)
        # ---- t ~ 31.6 s

        # ================================================= beat 4: coverage (32-45 s)
        self.play(FadeOut(chk_a), FadeOut(chk_b), FadeOut(chk_la), FadeOut(chk_lb),
                  FadeOut(eqsign), FadeOut(chk_r), FadeOut(note), FadeOut(eq_F),
                  FadeOut(planeF), FadeOut(fs_head), run_time=0.8)
        ballU = VGroup(kbox, kaxes)
        discsU = [central_disc(R, H, 0.03).shift(O3) for R in d["rotations"][:N_DRAWN]]
        headU = self.fx(st.T("uniform orientations", 30, st.TEXT).move_to(np.array([-3.3, 2.4, 0])))

        def cov_text(i):
            n = int(cov["counts"][i])
            t = VGroup(st.T(f"{n} planes", 30, st.TEXT),
                       st.T(f"{100 * cov['uniform'][i]:.0f}% of Fourier-ball voxels touched",
                            28, st.TRUTH)
                       ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
            t.move_to(np.array([3.3, 0.2, 0]), aligned_edge=LEFT).shift(LEFT * 1.75)
            return t.shift(min(0, 6.6 - t.get_right()[0]) * RIGHT)

        # one line, so it sits below the bottom of the k-box (y ~ -3.03) before and after the slide
        cov_note = self.fx(st.T("% = Fourier-ball voxels hit by ≥ 1 plane, 128³ grid "
                                "(coarser grid → higher %)", st.SMALL_SIZE, st.MUTED))
        cov_note.to_edge(DOWN, buff=st.MARGIN).align_to(np.array([-6.65, 0, 0]), LEFT)
        cnt = self.fx(cov_text(0))
        self.play(FadeIn(headU), FadeIn(cnt), FadeIn(cov_note),
                  LaggedStart(*[FadeIn(x) for x in discsU[:2]], lag_ratio=0.5), run_time=1.0)
        bounds = [2, 6, 20, 60]
        for j in (1, 2, 3):
            new = cov_text(j)
            self.play(LaggedStart(*[FadeIn(x) for x in discsU[bounds[j - 1]:bounds[j]]],
                                  lag_ratio=0.15), self.swap(cnt, new, 0.6), run_time=1.6)
            cnt = new
        self.wait(0.4)
        # ---- t ~ 37.0 s : side-by-side with a preferred-orientation set
        groupU = VGroup(ballU, *discsU)
        dU = self.world(-0.6, -0.1)
        offB = np.array([4.35, -0.1, 0])
        dB = self.world(*offB[:2])
        ballB = VGroup(box_edges(H), VGroup(*[Line(-1.95 * v, 1.95 * v, stroke_color=st.MUTED,
                                                   stroke_width=1.5) for v in (RIGHT, UP, OUT)]))
        discsB = [central_disc(R, H, 0.02, stroke_opacity=0.3) for R in d["biased_rotations"][:N_DRAWN]]
        groupB = VGroup(ballB, *discsB).shift(O3 + dB)
        headB = self.fx(st.T("preferred (tilt ≤ 12°)", 30, st.TEXT)
                        .move_to(np.array([offB[0] - 3.3, 2.4, 0])))
        headU2 = st.T("uniform orientations", 30, st.TEXT).move_to(np.array([-3.9, 2.4, 0]))
        sumU = st.T(f"{N_DRAWN} planes: {100 * cov['uniform_drawn']:.0f}% of voxels", 28, st.TRUTH)
        sumB = st.T(f"{N_DRAWN} planes: {100 * cov['biased_drawn']:.0f}% of voxels", 28, st.TRUTH)
        sumU.move_to(np.array([-3.9, 1.95, 0]))
        sumB.move_to(np.array([offB[0] - 3.3, 1.95, 0]))
        # k_z stays labelled on both balls (rotation about k_z leaves the axis fixed)
        kzU = st.M("k_z", size=36, color=st.MUTED).move_to(O3 + dU + 2.3 * OUT + 0.3 * RIGHT)
        kzB = st.M("k_z", size=36, color=st.MUTED).move_to(O3 + dB + 2.3 * OUT + 0.3 * RIGHT)
        self.add_fixed_orientation_mobjects(kzU, kzB)
        self.remove(kzU, kzB)
        # the empty cone (half-angle 78 deg about k_z) lies above and below the slab of planes
        missing = VGroup(st.T("empty cone", 28, st.TEXT),
                         VGroup(st.T("around", 28, st.TEXT),
                                st.M("k_z", size=36, color=st.TEXT)).arrange(RIGHT, buff=0.12),
                         st.T("(never measured)", 28, st.TEXT)
                         ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        missing.move_to(np.array([3.75, 0.75, 0]), aligned_edge=LEFT)
        missing.shift(min(0, 6.6 - missing.get_right()[0]) * RIGHT)     # keep the 0.45 margin
        cB_scr = np.array([ORIGIN_SCREEN[0] + offB[0], ORIGIN_SCREEN[1] + offB[1], 0])
        up_kz = self.rotm @ (1.2 * OUT)
        kz_pt = cB_scr + np.array([up_kz[0], up_kz[1], 0]) + 0.08 * RIGHT   # k_z axis above the slab
        miss_arr = Arrow(missing.get_left() + 0.12 * LEFT, kz_pt, buff=0.05, stroke_width=3,
                         color=st.TEXT, max_tip_length_to_length_ratio=0.12)
        self.fx(headU2, sumU, sumB, missing, miss_arr)
        self.play(groupU.animate.shift(dU), FadeOut(cnt), FadeOut(k_labels),
                  Transform(headU, headU2), run_time=1.0)
        self.play(FadeIn(ballB), FadeIn(headB), FadeIn(kzU), FadeIn(kzB),
                  LaggedStart(*[FadeIn(x) for x in discsB], lag_ratio=0.03), run_time=1.6)
        self.play(FadeIn(sumU), FadeIn(sumB), run_time=0.6)
        cB = O3 + dB
        self.play(Rotate(groupU, 35 * DEGREES, axis=OUT, about_point=O3 + dU),
                  Rotate(groupB, 35 * DEGREES, axis=OUT, about_point=cB),
                  FadeIn(missing, run_time=1.0), GrowArrow(miss_arr, run_time=1.0),
                  run_time=3.2, rate_func=smooth)
        self.wait(0.8)
        # ---- t ~ 45.0 s

        # ================================================= beat 5: reconstruction (45-57 s)
        self.play(FadeOut(groupU), FadeOut(groupB), FadeOut(headU), FadeOut(headB),
                  FadeOut(sumU), FadeOut(sumB), FadeOut(missing), FadeOut(miss_arr),
                  FadeOut(kzU), FadeOut(kzB), FadeOut(cov_note), run_time=0.8)
        RH = 3.5
        truth_img = st.density_image(_render_dict(rr, 0), st.TRUTH, RH).move_to(
            np.array([-4.9, -0.3, 0]))
        truth_lab = st.T("ground truth (synthetic)", 28, st.TRUTH).move_to(
            np.array([-4.9, 1.65, 0])).align_to(np.array([-6.65, 0, 0]), LEFT)
        contour_note = st.T("each map contoured at the\ntruth's enclosed volume", st.SMALL_SIZE,
                            st.MUTED).move_to(np.array([-4.9, -2.45, 0]))
        contour_note.align_to(np.array([-6.65, 0, 0]), LEFT)
        counts = [int(c) for c in d["counts"]]
        rcorr = [float(c) for c in d["recon_corr"]]
        RX = -0.35

        def rec_img(i):
            return st.density_image(_render_dict(rr, i), st.EST, RH).move_to(
                np.array([RX, -0.3, 0]))

        def rec_lab(i):
            if i <= 5:
                t = f"{counts[i - 1]} views, known poses"
            else:
                t = "300 views, tilt ≤ 12°\n(missing cone)"
            lab_ = st.T(t, 28, st.EST).move_to(np.array([RX, 1.65, 0]))
            # the two-line biased label sits lower (the biased map's top is lower) to clear eq_rec
            return lab_.align_to(np.array([0, 1.5 if i <= 5 else 1.15, 0]), DOWN)

        eq_rec = st.T("insert each 2D FT on its plane, then inverse 3D FFT", 28, st.TEXT)
        eq_rec.move_to(np.array([0.0, 2.45, 0]))
        eq_rec.align_to(np.array([-6.65, 0, 0]), LEFT)
        ideal = st.T("idealized: assumes known orientations and no CTF / noise", 28, st.EST,
                     slant="ITALIC").to_edge(DOWN, buff=st.MARGIN)
        ideal.align_to(np.array([-6.65, 0, 0]), LEFT)

        # correlation plot (computed values)
        axp = Axes(x_range=[0, 2.6, 1], y_range=[0, 1.0, 0.5], x_length=3.6, y_length=2.6,
                   tips=False, axis_config=dict(stroke_color=st.MUTED, stroke_width=2,
                                                include_ticks=True, tick_size=0.05))
        axp.move_to(np.array([4.55, -0.25, 0]))
        xt = VGroup(*[st.T(s, 28, st.MUTED).next_to(axp.c2p(x, 0), DOWN, buff=0.12)
                      for s, x in [("1", 0), ("10", 1), ("100", 2)]])
        yt = VGroup(*[st.T(s, 28, st.MUTED).next_to(axp.c2p(0, y), LEFT, buff=0.12)
                      for s, y in [("0", 0), ("0.5", 0.5), ("1", 1)]])
        xl = st.T("views N (log scale)", 28, st.MUTED).next_to(axp, DOWN, buff=0.55)
        yl = st.T("corr. with truth", 28, st.MUTED).next_to(axp, UP, buff=0.2)
        pts = [axp.c2p(np.log10(n), c) for n, c in zip(counts, rcorr)]
        dots = [Dot(p, radius=0.07, color=st.EST) for p in pts]
        segs = [Line(pts[i], pts[i + 1], stroke_color=st.EST, stroke_width=3)
                for i in range(len(pts) - 1)]
        pB = axp.c2p(np.log10(300), float(d["biased_corr"]))
        dotB = Circle(radius=0.09, stroke_color=st.EST, stroke_width=3).move_to(pB)
        labB = st.T(f"biased: {float(d['biased_corr']):.2f}", 28, st.EST).next_to(dotB, LEFT, buff=0.15)
        plot = VGroup(axp, xt, yt, xl, yl)
        # object z maps to R e_z = (0.88, 0, 0.47) under dm.display_R(): mostly horizontal on screen
        smear = st.T("missing cone → smeared along z (≈ horizontal)", 28, st.EST).next_to(
            xl, DOWN, buff=0.22)
        smear.shift(min(0, 6.6 - smear.get_right()[0]) * RIGHT)
        self.fx(truth_img, truth_lab, contour_note, eq_rec, ideal, plot, *dots, *segs, dotB, labB,
                smear)

        img = rec_img(1)
        lab = rec_lab(1)
        val = st.T(f"corr = {rcorr[0]:.2f}", 28, st.EST).move_to(np.array([RX, -2.3, 0]))
        self.fx(img, lab, val)
        self.play(FadeIn(truth_img), FadeIn(truth_lab), FadeIn(eq_rec), FadeIn(ideal),
                  FadeIn(plot), FadeIn(contour_note), run_time=1.0)
        self.play(FadeIn(img), FadeIn(lab), FadeIn(val), FadeIn(dots[0]), run_time=0.7)
        for i in range(2, 6):
            nimg, nlab = rec_img(i), rec_lab(i)
            nval = st.T(f"corr = {rcorr[i - 1]:.2f}", 28, st.EST).move_to(val)
            self.fx(nimg, nlab, nval)
            self.play(FadeOut(img), FadeIn(nimg), FadeOut(lab), FadeIn(nlab), FadeOut(val),
                      FadeIn(nval), Create(segs[i - 2]), FadeIn(dots[i - 1]), run_time=1.3)
            img, lab, val = nimg, nlab, nval
        self.wait(0.6)
        nimg, nlab = rec_img(6), rec_lab(6)
        nval = st.T(f"corr = {float(d['biased_corr']):.2f}", 28, st.EST).move_to(val)
        self.fx(nimg, nlab, nval)
        self.play(FadeOut(img), FadeIn(nimg), FadeOut(lab), FadeIn(nlab), FadeOut(val),
                  FadeIn(nval), FadeIn(dotB), FadeIn(labB), FadeIn(smear), run_time=1.2)
        img, lab, val = nimg, nlab, nval
        self.wait(1.9)
        # ---- t ~ 57.4 s

        # ================================================= beat 6: unknown poses (57-65 s)
        self.play(*[FadeOut(m) for m in (truth_img, truth_lab, contour_note, eq_rec, ideal, plot,
                                         *dots, *segs, dotB, labB, smear, img, lab, val)],
                  run_time=0.7)
        TH = 1.1
        # observed images I_i: the same projections plus seeded i.i.d. Gaussian noise (SNR 1)
        rng = np.random.default_rng(2027 + 30)
        noisy = [pr + rng.normal(0, pr.std(), pr.shape) for pr in d["show_proj"]]
        row = Group()
        known, unknown, icons, qicons = VGroup(), VGroup(), VGroup(), VGroup()
        for i in range(4):
            p = self.panel(noisy[i], TH, st.OBS, nearest=True)
            k = st.M(rf"R_{i + 1}", size=36, color=st.TEXT)
            u = st.M(rf"R_{i + 1}=?", size=36, color=st.EST)
            v = pj.viewing_direction(d["show_R"][i])        # pose icon: v projected on screen
            dv = np.array([v[0], v[1], 0.0])
            dv = dv / np.linalg.norm(dv) if np.linalg.norm(dv) > 0.2 else UP
            ic = Arrow(-0.3 * dv, 0.3 * dv, buff=0, stroke_width=4, color=st.TEXT,
                       max_tip_length_to_length_ratio=0.4)
            row.add(Group(p, k, u, ic))
        for g in row:
            g[1].next_to(g[0], RIGHT, buff=0.18).shift(0.25 * UP)
            g[2].next_to(g[0], RIGHT, buff=0.18).shift(0.25 * UP)
            g[3].next_to(g[0], RIGHT, buff=0.25).shift(0.3 * DOWN)
        cells = Group(*[Group(g[0], g[2]) for g in row]).arrange(RIGHT, buff=0.55)
        cells.move_to(np.array([0, 2.2, 0]))
        ilabs = VGroup(*[st.M(rf"I_{i + 1}", size=36, color=st.OBS).next_to(g[0], DOWN, buff=0.1)
                         for i, g in enumerate(row)])
        for g in row:
            g[1].move_to(g[2], aligned_edge=LEFT)
            g[3].next_to(g[1], DOWN, buff=0.12).align_to(g[1], LEFT)
            known.add(g[1])
            unknown.add(g[2])
            icons.add(g[3])
            q = VGroup(g[3].copy().set_color(st.EST).rotate(1.9),
                       st.T("?", 30, st.EST, weight="BOLD"))
            q[1].next_to(q[0], RIGHT, buff=0.06)
            q.move_to(g[3], aligned_edge=LEFT)
            qicons.add(q)
        imgs = Group(*[g[0] for g in row])
        eq_inv = st.M(r"I_i = P_{R_i}\rho+\epsilon_i", r"\quad\Longrightarrow\quad",
                      r"\text{infer }\rho\text{ and all }R_i\text{ jointly}", size=38)
        eq_inv.move_to(np.array([0, 0.68, 0]))
        eq_inv[2].set_color(st.EST)
        sim_note = st.T("simulated: i.i.d. Gaussian noise, SNR 1; CTF and shifts tᵢ omitted here",
                        st.SMALL_SIZE, st.MUTED).to_edge(DOWN, buff=st.MARGIN)
        sim_note.align_to(np.array([-6.65, 0, 0]), LEFT)
        self.fx(imgs, ilabs, known, unknown, icons, qicons, eq_inv, sim_note)
        self.play(FadeIn(imgs), FadeIn(ilabs), FadeIn(known), FadeIn(icons), FadeIn(sim_note),
                  run_time=0.7)
        self.wait(0.8)
        self.play(FadeOut(known), FadeOut(icons), FadeIn(unknown), FadeIn(qicons), run_time=1.0)
        self.play(Write(eq_inv), run_time=1.0)
        self.wait(0.6)

        card = st.takeaway_card(
            r"\widehat{P_R\rho}(k_x,k_y)=\hat\rho\!\left(R^{\top}(k_x,k_y,0)^{\top}\right)",
            "Each ideal view samples one central plane of the 3D FT;\n"
            "real data: poses are unknown and must be inferred with ρ.",
            "a lattice indexes hkl; here each pose must be inferred.")
        card.move_to(np.array([0, -1.75, 0]))
        if card.get_top()[1] > 0.3:
            card.shift((0.3 - card.get_top()[1]) * UP)
        # keep the bottom edge outside the 0.45 frame margin
        card.shift(max(0, st.MARGIN + 0.05 - (card.get_bottom()[1] + 4.0)) * UP)
        self.fx(card)
        # the card covers the bottom footnote; the SIMULATION badge stays on
        self.play(FadeIn(card), FadeOut(sim_note), run_time=0.8)
        self.wait(4.0)
