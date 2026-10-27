"""Clip 03 — Rotation in SO(3), projections and the Fourier slice theorem.

Storyboard beats (spec section 3, Scene 03) -> code sections:
  B0  title, source axes and the five dots
  B1  rotation in 3D: viewing direction v (2 DOF) and a SEPARATE in-plane angle psi
  B2  projection onto a detector plane along v; depth loss; four views
  B3  Fourier slice theorem: central plane through the k-origin rotating with R
  B4  sampling coverage: well-distributed versus biased orientations
  B5  inversion with known poses: estimated dots approach the true dots
  B6  unknown poses: orientation labels become question marks; joint estimation
  B7  takeaway freeze

Convention (helpers/projections.py): the source stays fixed in its own frame
(axes x1, x2, x3; frequency axes k1, k2, k3). R maps it into the detector frame,
v = R^T e_z is the projection axis, the detector axes are R^T e_x and R^T e_y,
and the projected dots are p_j = (R x_j)_{x,y}. The camera never moves: only
arrows, planes and dots move.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, projections as pj
from helpers import point_constellation as pc, noise_and_ctf as nc
from helpers.cache import cached

PHI, THETA = 62 * DEGREES, -50 * DEGREES      # fixed camera
ORIGIN_SCREEN = (-3.5, -0.45)                 # where the 3D origin appears (beats 1-3, 5)
S = 1.25            # scene units per source unit (beats 1, 2, 3, 5)
D = 1.6             # detector distance from the origin along v (source units)
H_SMALL, H_DET = 0.7, 1.25     # detector half-size: orientation sketch / projection
EX, EY, EZ = np.eye(3)
XMAX, YMAX = st.FRAME_W / 2 - st.MARGIN, st.FRAME_H / 2 - st.MARGIN
TAG_LEFT = st.illustration_label().get_left()[0]


def frame_of(a, th, ps):
    """R from ZYZ Euler angles (deg) and its detector axes in the source frame."""
    R = pj.rotation_from_euler(a, th, ps)
    return R, R.T @ EX, R.T @ EY, R.T @ EZ


def disc_points(n, radius, k=48):
    """Vertices of a flat disc of given radius through the origin, normal n."""
    n = n / np.linalg.norm(n)
    u = np.cross(n, EX if abs(n[0]) < 0.9 else EY)
    u /= np.linalg.norm(u)
    w = np.cross(n, u)
    t = np.linspace(0, 2 * np.pi, k, endpoint=False)
    return radius * (np.outer(np.cos(t), u) + np.outer(np.sin(t), w))


def nested_inversion(ks=(1, 2, 3, 5, 10, 20), meas_noise=0.03, n_cov=60):
    """Least-squares dot estimates from the first k known views (beat 5).

    Same uniform rotations as the coverage data (nc.rng(30), first draw), but the
    views are nested: each larger k keeps the earlier views AND their noise, so the
    error curve shows what adding views does (the shared helper redraws the noise
    for every k). Small ridge 1e-3 toward the origin, as in the helper."""
    X = pc.points()
    Ru = pj.random_rotations(n_cov, nc.rng(30))
    noise = meas_noise * nc.rng(31).standard_normal((max(ks), len(X), 2))
    est = []
    for k in ks:
        A = np.concatenate([R[:2] for R in Ru[:k]])
        est.append([np.linalg.solve(A.T @ A + 1e-3 * np.eye(3), A.T @ np.concatenate(
            [(R @ x)[:2] + noise[i, j] for i, R in enumerate(Ru[:k])]))
            for j, x in enumerate(X)])
    est = np.array(est)
    return dict(ks=np.array(ks), est=est, err=np.linalg.norm(est - X, axis=-1).mean(-1),
                meas_noise=meas_noise)


def check_layout(m, name=""):
    """Warn if a static overlay leaves the safe area or runs into the tag."""
    l, r = m.get_left()[0], m.get_right()[0]
    b, t = m.get_bottom()[1], m.get_top()[1]
    if l < -XMAX - 1e-3 or r > XMAX + 1e-3 or b < -YMAX - 1e-3 or t > YMAX + 1e-3:
        print(f"LAYOUT WARNING {name}: x[{l:.2f},{r:.2f}] y[{b:.2f},{t:.2f}]")
    if r > TAG_LEFT - 0.1 and t > 2.95:
        print(f"LAYOUT WARNING {name}: touches the tag (right {r:.2f}, top {t:.2f})")


class Scene03OrientationFourierSlice(ThreeDScene):
    # ------------------------------------------------------------ utilities
    def fix(self, *mobs, name=""):
        """Register as fixed-in-frame without adding (so FadeIn can add them)."""
        self.add_fixed_in_frame_mobjects(*mobs)
        self.remove(*mobs)
        for m in mobs:
            check_layout(m, name)
        return mobs[0] if len(mobs) == 1 else mobs

    def world(self, sx, sy):
        """3D point that the fixed camera shows at frame position (sx, sy)."""
        return self.rot.T @ np.array([sx, sy, 0.0])

    def screen(self, p):
        q = self.renderer.camera.project_point(np.asarray(p, float))
        return np.array([q[0], q[1], 0.0])

    def arrow3(self, start_fn, end_fn, color, width=5, tip_len=0.24, tip_w=0.1):
        """Cheap 3D arrow: a projected line plus a flat tip that faces the camera,
        updated in place (rebuilding Arrow3D every frame is far too slow)."""
        shaft = Line(ORIGIN, RIGHT, stroke_color=color, stroke_width=width)
        tip = Polygon(ORIGIN, RIGHT, UP, fill_color=color, fill_opacity=1,
                      stroke_color=color, stroke_width=1)
        view = self.rot[2]

        def upd(g):
            s, e = np.asarray(start_fn(), float), np.asarray(end_fn(), float)
            u = (e - s) / max(np.linalg.norm(e - s), 1e-6)
            perp = np.cross(u, view)
            perp = perp / np.linalg.norm(perp) if np.linalg.norm(perp) > 1e-3 else self.rot[0]
            base = e - tip_len * u
            shaft.put_start_and_end_on(s, base)
            tip.set_points_as_corners([e, base + tip_w * perp, base - tip_w * perp, e])
        g = VGroup(shaft, tip)
        upd(g)
        g.add_updater(upd)
        return g

    def axes3d(self, origin, scale, half, names, color=st.GRID, pad=0.2):
        ax = ThreeDAxes(x_range=(-half, half, 1), y_range=(-half, half, 1),
                        z_range=(-half, half, 1), x_length=2 * half * scale,
                        y_length=2 * half * scale, z_length=2 * half * scale,
                        axis_config={"include_tip": False, "include_ticks": False,
                                     "stroke_color": color, "stroke_width": 2})
        ax.shift(origin - ax.c2p(0, 0, 0))
        labs = VGroup()
        for e, nm in zip(np.eye(3), names):
            lab = st.M(nm, size=36, color=st.MUTED)
            # x2 rises to the upper right with dot 5 just above it: label below the end
            off = 0.32 * DOWN + 0.05 * RIGHT if nm in ("x_2",) else ORIGIN
            lab.move_to(self.screen(origin + (half + pad) * scale * e) + off)
            lab.add_background_rectangle(color=st.BG, opacity=0.85, buff=0.04)
            labs.add(lab)
        self.fix(labs, name="axis labels")
        return ax, labs

    def header(self, text):
        h = st.fit_width(st.T(text, 36, weight="BOLD"), TAG_LEFT - 0.3 + st.FRAME_W / 2
                         - st.MARGIN)
        h.to_corner(UL, buff=st.MARGIN)
        return self.fix(h, name=text)

    def swap(self, old, new, rt=0.8):
        """Fade one fixed-in-frame label out, then the next one in (no ghosting)."""
        return Succession(FadeOut(old, shift=0.1 * UP, run_time=rt / 2),
                          FadeIn(new, shift=0.1 * UP, run_time=rt / 2))

    # ------------------------------------------------------------ the clip
    def construct(self):
        d = cached("orientation_fourier", ne.orientation_fourier)
        X, views, proj, imgs, ftm = d["X"], d["views"], d["proj"], d["imgs"], d["ft_mag"]
        for i, v3 in enumerate(views):     # the live geometry must reproduce the data
            assert np.allclose(pj.project_points(frame_of(*v3)[0], X), proj[i])

        self.set_camera_orientation(phi=PHI, theta=THETA)
        self.renderer.camera.reset_rotation_matrix()
        self.rot = self.renderer.camera.get_rotation_matrix()
        tag = st.illustration_label()
        self.add_fixed_in_frame_mobjects(tag)

        O = self.world(*ORIGIN_SCREEN)        # 3D origin, left half of the frame

        def P(x):
            return O + S * np.asarray(x)

        # ============================================================ B0 title
        title = st.title_card("Rotation in SO(3), projections\nand the Fourier slice theorem",
                              "How do unknown 3D orientations and\n"
                              "Fourier slices define reconstruction?")
        if title[0].get_right()[0] > TAG_LEFT - 0.3:          # keep clear of the tag
            title[0].scale((TAG_LEFT - 0.3 + XMAX) / title[0].width,
                           about_point=title[0].get_corner(UL))
            title[1].next_to(title[0], DOWN, aligned_edge=LEFT, buff=0.15)
        self.fix(title, name="title")
        ax, axlab = self.axes3d(O, S, 1.6, [r"x_1", r"x_2", r"x_3"])
        dz = ValueTracker(0.0)                # depth slide of dot 3 (beat 2)
        a, th, ps = ValueTracker(0.0), ValueTracker(0.0), ValueTracker(0.0)
        hs = ValueTracker(H_SMALL)

        def cur():
            return frame_of(a.get_value(), th.get_value(), ps.get_value())

        def v_now():
            return cur()[3]

        def src_pts():
            Y = X.copy()
            Y[2] = Y[2] + dz.get_value() * v_now()
            return Y

        dots = VGroup(*[Dot3D(P(x), radius=0.1, color=st.TRUTH) for x in X])
        for j, dt in enumerate(dots):
            dt.add_updater(lambda m, j=j: m.move_to(P(src_pts()[j])))

        self.play(FadeIn(title), Create(ax), run_time=1.0)
        self.play(LaggedStart(*[FadeIn(dt, scale=0.5) for dt in dots], lag_ratio=0.2),
                  run_time=1.2)
        self.add(dots)
        self.wait(0.8)

        # ============================================================ B1 rotation
        hdr = self.header("1  Orientation  R ∈ SO(3)")
        self.play(self.swap(title, hdr, 0.8), FadeIn(axlab, run_time=0.8))

        def detector():
            _, e1, e2, v = cur()
            h = hs.get_value()
            pts = [P(D * v + sx * h * e1 + sy * h * e2) for sx, sy in
                   ((-1, -1), (1, -1), (1, 1), (-1, 1))]
            return Polygon(*pts, stroke_color=st.MUTED, stroke_width=2,
                           fill_color=st.MUTED, fill_opacity=0.14)

        v_mob = self.arrow3(lambda: P(0 * v_now()), lambda: P((D + 0.55) * v_now()), st.EST)
        det_mob = always_redraw(detector)
        e1_mob = self.arrow3(lambda: P(D * v_now()), lambda: P(D * v_now() + 0.62 * cur()[1]),
                             st.EST, width=4, tip_len=0.17, tip_w=0.075)
        v_lab = self.fix(st.M(r"\mathbf v", size=40, color=st.EST))
        v_lab.add_updater(lambda m: m.move_to(self.screen(P((D + 0.85) * v_now()))))

        t_v = self.fix(VGroup(
            st.M(r"\mathbf v = R^{\top}\mathbf e_z", size=40, color=st.EST),
            st.M(r"R = R_z(\psi)\,R_y(\theta)\,R_z(\alpha)", size=36, color=st.MUTED),
            st.T("viewing direction: 2 DOF (α, θ)", 30),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([0.45, 2.3, 0], aligned_edge=UL),
            name="t_v")
        # the x3 label would sit inside the detector sketch: hide it while the detector shows
        self.play(FadeIn(v_mob), FadeIn(det_mob), FadeIn(e1_mob), FadeIn(v_lab),
                  FadeOut(axlab[2]),
                  FadeIn(t_v), run_time=0.8)
        self.play(a.animate.set_value(55), th.animate.set_value(60), run_time=2.4,
                  rate_func=smooth)

        # in-plane angle psi: v fixed, the detector frame spins about v
        def psi_arc():
            ang = ps.get_value()
            if ang < 1:
                return VMobject()
            pts = []
            for s in np.linspace(0, ang, 24):
                _, e1s, _, vv = frame_of(a.get_value(), th.get_value(), s)
                pts.append(P(D * vv + 0.5 * e1s))
            arc = VMobject(stroke_color=st.EST, stroke_width=4)
            arc.set_points_smoothly(pts)
            return arc

        arc_mob = always_redraw(psi_arc)
        psi_lab = st.M(r"\psi", size=40, color=st.EST)
        psi_lab.add_background_rectangle(color=st.BG, opacity=0.85, buff=0.04)
        self.fix(psi_lab)

        def psi_pos(m):
            _, e1s, _, vv = frame_of(a.get_value(), th.get_value(), ps.get_value() / 2)
            m.move_to(self.screen(P(D * vv + 1.2 * e1s)))
        psi_lab.add_updater(psi_pos)
        t_psi = self.fix(VGroup(
            st.T("ψ: in-plane rotation about v", 30),
            st.T("(1 DOF; v does not change)", 28, st.MUTED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12).next_to(t_v, DOWN, aligned_edge=LEFT,
                                                               buff=0.45), name="t_psi")
        self.add(arc_mob)
        self.play(FadeIn(t_psi), FadeIn(psi_lab), run_time=0.5)
        self.play(ps.animate.set_value(80), run_time=2.2, rate_func=smooth)

        dof_t = VGroup(st.T("Viewing direction", 30), st.T("+ in-plane rotation", 30),
                       st.T("= 3 orientation DOF", 30)).arrange(DOWN, aligned_edge=LEFT,
                                                                buff=0.1)
        dof_box = SurroundingRectangle(dof_t, buff=0.18, color=st.MUTED, stroke_width=1.5,
                                       corner_radius=0.08)
        dof = self.fix(VGroup(dof_box, dof_t).next_to(t_psi, DOWN, aligned_edge=LEFT,
                                                     buff=0.45), name="dof")
        so3 = self.fix(st.M(r"R^{\top}R = I,\;\; \det R = 1,\;\; R\in SO(3)", size=38)
                       .next_to(dof, DOWN, aligned_edge=LEFT, buff=0.4), name="so3")
        self.play(FadeIn(dof), run_time=0.7)
        self.play(FadeIn(so3), run_time=0.7)
        self.wait(1.6)

        # ============================================================ B2 projection
        hdr2 = self.header("2  Projection onto a detector")
        psi_lab.clear_updaters()
        self.play(self.swap(hdr, hdr2), FadeOut(VGroup(t_v, t_psi, dof, so3)),
                  FadeOut(psi_lab), FadeOut(e1_mob), run_time=0.8)
        self.remove(arc_mob)
        self.play(a.animate.set_value(views[0][0]), th.animate.set_value(views[0][1]),
                  ps.animate.set_value(views[0][2]), run_time=2.2,
                  rate_func=rate_functions.ease_in_out_sine)   # flatter peak speed than smooth

        def proj_now():
            return pj.project_points(cur()[0], src_pts())

        def rays():
            _, e1, e2, v = cur()
            g = VGroup()
            for x in src_pts():
                q = x + (D - x @ v) * v
                g.add(DashedLine(P(x), P(q), dash_length=0.08, stroke_width=1.8,
                                 stroke_color=st.MUTED, stroke_opacity=0.85))
            return g

        def det_pos(j):
            _, e1, e2, v = cur()
            p = proj_now()[j]
            return P(D * v + p[0] * e1 + p[1] * e2)

        ray_mob = always_redraw(rays)
        pdot_mob = VGroup(*[Dot3D(det_pos(j), radius=0.07, color=st.OBS, resolution=(6, 6))
                            for j in range(len(X))])
        for j, dt in enumerate(pdot_mob):
            dt.add_updater(lambda m, j=j: m.move_to(det_pos(j)))

        # 2D detector-coordinate panel (right)
        pan = Axes(x_range=(-1.5, 1.5, 0.5), y_range=(-1.5, 1.5, 0.5), x_length=2.5,
                   y_length=2.5, tips=False,
                   axis_config={"stroke_color": st.GRID, "include_ticks": False})
        pan_box = SurroundingRectangle(pan, buff=0.06, color=st.GRID, stroke_width=1.5)
        pan_g = VGroup(pan_box, pan).move_to([2.05, 0.95, 0])
        pdots2 = VGroup(*[Dot(pan.c2p(*p), radius=0.08, color=st.OBS) for p in proj[0]])
        for j, dt in enumerate(pdots2):
            dt.add_updater(lambda m, j=j: m.move_to(pan.c2p(*proj_now()[j])))
        pan_lab = st.T("detector coordinates (x, y)", 28, st.MUTED).next_to(pan_box, UP,
                                                                              buff=0.12)
        pan_lab.align_to(pan_box, LEFT)
        p_eq = st.M(r"\mathbf p_j = (R\,\mathbf x_j)_{x,y}", size=40, color=st.OBS)
        p_note = VGroup(st.T("depth", 28, st.MUTED),
                        st.M(r"(R\,\mathbf x_j)_z", size=36, color=st.MUTED)).arrange(RIGHT,
                                                                                     buff=0.12)
        p_note2 = st.T("is discarded", 28, st.MUTED)
        p_g = VGroup(p_eq, p_note, p_note2).arrange(DOWN, aligned_edge=LEFT, buff=0.15)
        p_g.next_to(pan_box, RIGHT, buff=0.3).align_to(pan_box, UP).shift(DOWN * 0.3)
        self.fix(pan_g, pdots2, pan_lab, p_g, name="panel")
        self.play(hs.animate.set_value(H_DET), FadeIn(pan_g), FadeIn(pan_lab), run_time=1.0)
        self.play(FadeIn(ray_mob), FadeIn(pdot_mob), FadeIn(pdots2), FadeIn(p_g), run_time=1.0)

        # depth loss: slide one dot along v; its projection does not move
        cap = self.fix(VGroup(st.T("Sliding a dot along v leaves", 28),
                              st.T("its image unchanged:", 28),
                              st.T("one projection loses depth", 28))
                       .arrange(DOWN, aligned_edge=LEFT, buff=0.1)
                       .move_to([0.45, -0.75, 0], aligned_edge=UL), name="cap")
        ring = always_redraw(lambda: Circle(radius=0.22, color=st.TRUTH, stroke_width=3)
                             .move_to(self.screen(P(src_pts()[2]))))
        self.fix(ring)
        self.play(FadeIn(cap), FadeIn(ring), run_time=0.5)
        self.play(dz.animate.set_value(1.25), run_time=1.4, rate_func=smooth)
        self.play(dz.animate.set_value(0.0), run_time=1.1, rate_func=smooth)
        self.play(FadeOut(ring), run_time=0.3)
        self.remove(ring)

        # four views -> four distinct 2D patterns (thumbnails = images P_R rho)
        thumbs = Group()
        for i in range(4):
            hm = st.heatmap(imgs[i], height=1.0)
            lab = st.M(rf"R_{i + 1}", size=36, color=st.EST)
            thumbs.add(Group(hm, st.framed(hm), lab))
        for i, g in enumerate(thumbs):
            g[0].move_to([1.75 + 1.3 * i, -2.25, 0])
            g[1].move_to(g[0])
            g[2].next_to(g[0], DOWN, buff=0.1)
        thumbs_lab = st.T("images", 28, st.MUTED).next_to(thumbs[0][0], LEFT, buff=0.2)
        self.fix(thumbs, thumbs_lab, name="thumbs")
        self.play(FadeOut(cap), FadeIn(thumbs[0]), FadeIn(thumbs_lab), run_time=0.6)
        for i in (1, 2, 3):
            rt = {1: 1.5, 2: 1.5, 3: 2.6}[i]      # R3 -> R4 is the largest turn: slower
            # ease_in_out_sine: peak slope pi/2 (smooth: 2.5), so no fast spin mid-move
            self.play(a.animate.set_value(views[i][0]), th.animate.set_value(views[i][1]),
                      ps.animate.set_value(views[i][2]), run_time=rt,
                      rate_func=rate_functions.ease_in_out_sine)
            self.play(FadeIn(thumbs[i]), run_time=0.4)
        self.wait(0.6)

        # ============================================================ B3 Fourier slice
        hdr3 = self.header("3  Fourier slice theorem")
        klab = VGroup()
        for e, nm in zip(np.eye(3), ["k_1", "k_2", "k_3"]):
            # plain axis-end positions (the x2 label's downward nudge is for the dots only)
            klab.add(st.M(nm, size=36, color=st.MUTED).move_to(self.screen(O + 2.0 * S * e))
                     .add_background_rectangle(color=st.BG, opacity=0.85, buff=0.04))
        self.fix(klab)
        # x3 (axlab[2]) is already hidden: fading the whole group would flash it back in
        self.play(self.swap(hdr2, hdr3), self.swap(VGroup(axlab[0], axlab[1]), klab, 0.9),
                  FadeOut(VGroup(dots, ray_mob, pdot_mob, det_mob, pdots2, pan_g, pan_lab, p_g,
                                 v_mob, v_lab)),
                  FadeOut(thumbs), FadeOut(thumbs_lab), run_time=0.9)
        for m in (dots, pdot_mob):
            for dt in m:
                dt.clear_updaters()
        # start from the pose B2 ends on (R4, oblique), then edge-on R3, then flat R1:
        # three plainly different planes and no hidden pose jump
        b3_views = (3, 2, 0)
        for tr, val in zip((a, th, ps), views[b3_views[0]]):
            tr.set_value(val)

        def central_plane():
            _, e1, e2, v = cur()
            r = 1.45
            pts = [P(r * (np.cos(t) * e1 + np.sin(t) * e2))
                   for t in np.linspace(0, 2 * np.pi, 48, endpoint=False)]
            disc = Polygon(*pts, stroke_color=st.OBS, stroke_width=2, fill_color=st.OBS,
                           fill_opacity=0.22)
            l1 = Line(P(-r * e1), P(r * e1), stroke_color=st.OBS, stroke_width=1.5,
                      stroke_opacity=0.7)
            l2 = Line(P(-r * e2), P(r * e2), stroke_color=st.OBS, stroke_width=1.5,
                      stroke_opacity=0.7)
            return VGroup(disc, l1, l2)

        plane_mob = always_redraw(central_plane)
        kv_mob = self.arrow3(lambda: P(0 * v_now()), lambda: P(1.3 * v_now()), st.EST)
        v_lab.clear_updaters()
        v_lab.add_background_rectangle(color=st.BG, opacity=0.85, buff=0.04)   # crosses the rim
        v_lab.add_updater(lambda m: m.move_to(self.screen(P(1.3 * v_now())) + 0.38 * RIGHT))

        plane_cap = self.fix(VGroup(
            st.T("Central plane ⟂ v of the 3D FT", 28),
            st.T("not a physical cut of the source", 28, st.MUTED),
            st.T("ideal projection: no CTF, shift or noise", 26, st.MUTED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to([0.45, 2.75, 0], aligned_edge=UL),
            name="plane_cap")

        # right: image -> FT2 -> log-magnitude display (central half of the frequency grid).
        # The column is stacked from plane_cap downwards with explicit gaps; the numerical
        # check goes to the free lower-left corner so the column can breathe.
        n = ftm.shape[1]
        c0, c1 = n // 4, 3 * n // 4
        H_PAIR = 1.5
        y_pair = plane_cap.get_bottom()[1] - 0.22 - H_PAIR / 2

        def pair(i):
            im = st.heatmap(imgs[i], height=H_PAIR).move_to([1.45, y_pair, 0])
            fm = st.heatmap(ftm[i][c0:c1, c0:c1], height=H_PAIR).move_to([4.95, y_pair, 0])
            return Group(im, st.framed(im), fm, st.framed(fm, st.OBS))

        pr = self.fix(pair(b3_views[0]), name="pair")
        ft_arrow = self.fix(VGroup(
            Arrow(LEFT * 0.6, RIGHT * 0.6, buff=0, color=st.TEXT, stroke_width=4),
            st.M(r"\mathcal F_2", size=36)).arrange(DOWN, buff=0.08).move_to([3.2, y_pair, 0]))
        lab_im = st.M(r"P_R\rho", size=36).next_to(pr[0], DOWN, buff=0.12)
        lab_ft = st.M(r"\log(1+|\widehat{P_R\rho}|)", size=36, color=st.OBS).next_to(
            pr[2], DOWN, buff=0.12)
        lab_bot = min(lab_im.get_bottom()[1], lab_ft.get_bottom()[1])
        mag_note = st.T("magnitude shown for display only; the\n"
                        "theorem holds for the complex values", 26, st.MUTED,
                        line_spacing=0.9).move_to([0.45, lab_bot - 0.2, 0], aligned_edge=UL)
        fs_eq = st.M(r"\widehat{P_R\rho}(k_x,k_y)=\widehat{\rho}\!\left(R^{\top}"
                     r"(k_x,k_y,0)^{\top}\right)", size=40)
        fs_eq.next_to(mag_note, DOWN, buff=0.25, aligned_edge=LEFT)
        self.fix(lab_im, lab_ft, mag_note, fs_eq, name="B3 right")

        self.play(FadeIn(plane_mob), FadeIn(kv_mob), FadeIn(v_lab), FadeIn(plane_cap),
                  run_time=0.9)
        self.play(FadeIn(pr), FadeIn(ft_arrow), FadeIn(lab_im), FadeIn(lab_ft), run_time=0.9)
        self.play(FadeIn(mag_note), run_time=0.5)
        self.play(FadeIn(fs_eq), run_time=0.8)
        self.wait(0.8)
        for i in b3_views[1:]:
            new = self.fix(pair(i))
            if i == 2:          # R4 -> R3: the 171 deg turn
                rt, fades = 2.8, (0.5, 1.6, 0.6)
            else:               # R3 -> R1 (alpha 120 -> 0, theta 50 -> 0)
                rt, fades = 2.2, (0.5, 1.1, 0.6)
            self.play(a.animate.set_value(views[i][0]), th.animate.set_value(views[i][1]),
                      ps.animate.set_value(views[i][2]),
                      Succession(FadeOut(pr, run_time=fades[0]), Wait(fades[1]),
                                 FadeIn(new, run_time=fades[2])),
                      run_time=rt, rate_func=rate_functions.ease_in_out_sine)
            pr = new
            self.wait(0.3)
        cRT, cR = d["corr_RT"], d["corr_R"]
        chk = self.fix(VGroup(
            # the data check covers views 2-4 (R1 = I, where R = R^T makes the test moot)
            st.T(f"Check, views 2-{len(cRT) + 1}: complex, linear interp.", 26, st.MUTED),
            # floor, not round: 0.99986 must read "≥ 0.9998", never "≥ 0.9999"
            st.T(f"correlation via Rᵀ: ≥ {np.floor(cRT.min() * 1e4) / 1e4:.4f}", 28),
            st.T(f"correlation via R: {cR.min():.2f}–{cR.max():.2f}", 28),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08).move_to([-6.5, -2.4, 0], aligned_edge=UL),
            name="chk")
        if chk.get_right()[0] > fs_eq.get_left()[0] - 0.2:      # shares a row with fs_eq
            print(f"LAYOUT WARNING chk runs into fs_eq: {chk.get_right()[0]:.2f}")
        self.play(FadeIn(chk), run_time=0.7)
        self.wait(1.8)

        # ============================================================ B4 coverage
        hdr4 = self.header("4  Coverage of Fourier space")
        v_lab.clear_updaters()
        self.play(self.swap(hdr3, hdr4),
                  FadeOut(VGroup(plane_mob, kv_mob, ax, klab, v_lab, plane_cap, ft_arrow,
                                 lab_im, lab_ft, mag_note, fs_eq, chk)),
                  FadeOut(pr), run_time=0.9)

        S2, half2, rad2 = 1.0, 1.55, 1.3
        OL, OR = self.world(-3.4, -0.2), self.world(3.15, -0.2)
        axL, labL = self.axes3d(OL, S2, half2, ["k_1", "k_2", "k_3"], pad=0.55)
        axR, labR = self.axes3d(OR, S2, half2, ["k_1", "k_2", "k_3"], pad=0.75)

        def discs(origin, normals):
            g = VGroup()
            for nv in normals:
                pts = [origin + S2 * p for p in disc_points(nv, rad2)]
                g.add(Polygon(*pts, stroke_color=st.OBS, stroke_width=1.5, stroke_opacity=0.75,
                              fill_color=st.OBS, fill_opacity=0.05))
            return g

        dL = discs(OL, d["normals_uniform"])
        dR = discs(OR, d["normals_biased"])
        tilt = float(d["max_theta_biased"])
        ttl_L = st.T("well-distributed views", 30).move_to([-3.4, -2.55, 0])
        ttl_R = VGroup(st.T(f"biased views (v within {tilt:.0f}° of", 30),
                       st.M("k_3", size=36), st.T(")", 30))
        # word space before k3, a tight gap before the closing parenthesis
        ttl_R[1].next_to(ttl_R[0], RIGHT, buff=0.15)
        ttl_R[2].next_to(ttl_R[1], RIGHT, buff=0.03)
        ttl_R[1].align_to(ttl_R[0], UP).shift(0.02 * DOWN)
        ttl_R[2].align_to(ttl_R[0], UP)
        ttl_R.move_to([3.15, -2.55, 0])
        ncov = int(d["n_cov"])
        cov_L = st.T(f"{ncov} planes sample {100 * float(d['cov_uniform']):.0f}% of the ball",
                     28, st.OBS).next_to(ttl_L, DOWN, buff=0.15)
        cov_R = st.T(f"{ncov} planes sample {100 * float(d['cov_biased']):.0f}% of the ball",
                     28, st.OBS).next_to(ttl_R, DOWN, buff=0.15)
        draw_note = st.T(f"(first {len(dL)} of {ncov} planes drawn; % = voxels of a 32³ "
                         "ball hit)", 26, st.MUTED)
        draw_note.move_to([-0.15, 2.45, 0])
        self.fix(ttl_L, ttl_R, cov_L, cov_R, draw_note, name="B4 text")
        self.play(Create(axL), Create(axR), FadeIn(labL), FadeIn(labR), FadeIn(ttl_L),
                  FadeIn(ttl_R), run_time=0.9)
        self.play(LaggedStart(*[FadeIn(p) for p in dL], lag_ratio=0.25),
                  LaggedStart(*[FadeIn(p) for p in dR], lag_ratio=0.25),
                  FadeIn(draw_note), run_time=3.0)
        self.play(FadeIn(cov_L), FadeIn(cov_R), run_time=0.7)

        # missing cone around k3 (biased side): boundary at (90 - tilt) deg from k3
        cone_half = np.radians(90 - tilt)
        sr = self.rot[0] - (self.rot[0] @ EZ) * EZ           # screen-right, horizontal
        sr /= np.linalg.norm(sr)
        cone = VGroup()
        for sgn_z in (1, -1):
            for sgn_r in (1, -1):
                dvec = np.cos(cone_half) * sgn_z * EZ + np.sin(cone_half) * sgn_r * sr
                cone.add(DashedLine(OR, OR + S2 * 1.9 * dvec, dash_length=0.08,
                                    stroke_color=st.TEXT, stroke_width=2.5))
        # the unsampled region itself: two translucent wedges between the rays and k3
        Lw = S2 * 1.9
        wedge = VGroup(*[Polygon(OR, OR + Lw * (np.cos(cone_half) * s * EZ + np.sin(cone_half) * sr),
                                 OR + Lw * s * EZ,
                                 OR + Lw * (np.cos(cone_half) * s * EZ - np.sin(cone_half) * sr),
                                 stroke_width=0, fill_color=st.MUTED, fill_opacity=0.18)
                         for s in (1, -1)])
        tip = self.screen(OR + S2 * 1.2 * EZ)
        miss = VGroup(st.T("missing cone", 28),
                      VGroup(st.T("around", 28), st.M("k_3", size=36)).arrange(RIGHT, buff=0.12),
                      st.T("no plane", 28), st.T("reaches it", 28)).arrange(DOWN, aligned_edge=LEFT,
                                                                           buff=0.12)
        miss.move_to([-0.1, tip[1] - 0.1, 0])
        pt_arrow = Arrow(miss.get_right() + RIGHT * 0.05, tip + LEFT * 0.1, buff=0.05,
                         color=st.TEXT, stroke_width=4, max_tip_length_to_length_ratio=0.3)
        self.fix(miss, pt_arrow, name="miss")
        self.play(FadeIn(wedge), Create(cone), FadeIn(miss), GrowArrow(pt_arrow), run_time=1.0)
        self.wait(2.2)

        # ============================================================ B5 inversion
        hdr5 = self.header("5  Inversion with known poses")
        self.play(self.swap(hdr4, hdr5),
                  FadeOut(VGroup(axL, axR, labL, labR, dL, dR, ttl_L, ttl_R, cov_L, cov_R,
                                 draw_note, cone, wedge, miss, pt_arrow)), run_time=0.9)
        ax5, axlab5 = self.axes3d(O, S, 1.6, [r"x_1", r"x_2", r"x_3"])
        # nested views: each k keeps the earlier views and their noise (see nested_inversion)
        inv = cached("s03_orientation_fourier_slice_inversion", nested_inversion)
        ks, est, err = inv["ks"], inv["est"], inv["err"]
        # true positions as green rings, estimates as orange dots
        rings = VGroup(*[Circle(radius=0.17, color=st.TRUTH, stroke_width=3.5)
                         .move_to(self.screen(P(x))) for x in X])
        est5 = VGroup(*[Dot3D(P(x), radius=0.09, color=st.EST) for x in est[0]])

        # error-vs-k plot (log-log), values from err_uniform
        pl = Axes(x_range=(0, 1.4, 0.2), y_range=(0, 2.4, 0.5), x_length=4.4, y_length=2.15,
                  tips=False, axis_config={"stroke_color": st.GRID, "include_ticks": False})
        pl.move_to([3.85, 1.35, 0])

        def lc(x, y):                     # (log10 k, log10 error) -> plot point
            return pl.c2p(x, y + 2.2)

        ticks = VGroup()
        for kk in (1, 2, 5, 10, 20):
            ticks.add(st.T(str(kk), 26, st.MUTED).next_to(lc(np.log10(kk), -2.2), DOWN,
                                                          buff=0.12),
                      Line(lc(np.log10(kk), -2.2), lc(np.log10(kk), -2.1),
                           stroke_color=st.GRID))
        for yv, ys in ((-2, "0.01"), (-1, "0.1"), (0, "1")):
            ticks.add(st.T(ys, 26, st.MUTED).next_to(lc(0, yv), LEFT, buff=0.12),
                      Line(lc(0, yv), lc(0.04, yv), stroke_color=st.GRID))
        xl = st.T("number of known views k", 28, st.MUTED).next_to(pl, DOWN, buff=0.5)
        yl = st.T("mean dot error (units)", 28, st.MUTED).next_to(pl, UP, buff=0.15)
        yl.align_to(pl, LEFT)
        pts5 = VGroup(*[Dot(lc(np.log10(k), np.log10(e)), radius=0.07, color=st.EST)
                        for k, e in zip(ks, err)])
        seg5 = VGroup(*[Line(pts5[i].get_center(), pts5[i + 1].get_center(),
                             stroke_color=st.EST, stroke_width=2.5) for i in range(len(ks) - 1)])

        def readout(i):
            word = "view" if ks[i] == 1 else "views"
            r = st.T(f"k = {ks[i]} {word}: mean error {err[i]:.3f}", 28, st.EST)
            return r.move_to([0.45, xl.get_bottom()[1] - 0.25, 0], aligned_edge=UL)

        rd = self.fix(readout(0), name="readout")

        def note(*lines):
            return VGroup(*[st.T(s, 26, st.MUTED) for s in lines]).arrange(
                DOWN, aligned_edge=LEFT, buff=0.08).next_to(rd, DOWN, buff=0.15,
                                                           aligned_edge=LEFT)
        # k = 1: the ridge (1e-3, toward the origin) puts the undetermined depth at 0
        k1_note = self.fix(note("k=1: depth along that view's v", "is undetermined (set to 0)"),
                           name="k1")
        # from k = 2: labelled dots are already fixed, unlike a density (beat 4)
        k2_note = self.fix(note("2 views already fix labelled dots;",
                                "a density needs many views (see 4)"), name="k2")
        meth = self.fix(VGroup(
            st.T("least squares; poses, dot labels known", 26, st.MUTED),
            st.T(f"cumulative views, noise s.d. {float(inv['meas_noise']):.2f} units", 26,
                 st.MUTED))
            .arrange(DOWN, aligned_edge=LEFT, buff=0.08).next_to(k1_note, DOWN, buff=0.2,
                                                                 aligned_edge=LEFT), name="meth")
        cap5 = self.fix(st.caption("Illustration of reconstruction; actual algorithms "
                                   "estimate a 3D density.", st.TEXT), name="cap5")
        leg = self.fix(VGroup(
            VGroup(Circle(radius=0.12, color=st.TRUTH, stroke_width=3.5),
                   st.T("true dots", 28)).arrange(RIGHT, buff=0.18),
            VGroup(Dot(radius=0.09, color=st.EST), st.T("estimates", 28)).arrange(RIGHT,
                                                                                 buff=0.21),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.15).move_to([-6.5, 2.55, 0], aligned_edge=UL),
            name="legend")
        self.fix(rings, pl, ticks, xl, yl, pts5, seg5, name="plot")
        self.play(Create(ax5), FadeIn(axlab5), FadeIn(rings), FadeIn(leg), FadeIn(pl),
                  FadeIn(ticks), FadeIn(xl), FadeIn(yl), FadeIn(cap5), run_time=0.9)
        self.play(FadeIn(est5), FadeIn(pts5[0]), FadeIn(rd), FadeIn(k1_note), FadeIn(meth),
                  run_time=0.7)
        self.wait(1.0)
        for i in range(1, len(ks)):
            new_rd = self.fix(readout(i))
            anims = [est5[j].animate.move_to(P(est[i][j])) for j in range(len(X))]
            out = ([Succession(FadeOut(k1_note, run_time=0.4), FadeIn(k2_note, run_time=0.4))]
                   if i == 1 else [])
            self.play(*anims, Create(seg5[i - 1]), FadeIn(pts5[i]), *out,
                      Succession(FadeOut(rd, run_time=0.3), FadeIn(new_rd, run_time=0.3)),
                      run_time=1.0)
            rd = new_rd
            self.wait(0.25)
        self.wait(1.0)

        # ============================================================ B6 unknown poses
        hdr6 = self.header("6  In practice the poses are unknown")
        self.play(self.swap(hdr5, hdr6),
                  FadeOut(VGroup(ax5, axlab5, rings, est5, leg, pl, ticks, xl, yl, pts5, seg5,
                                 rd, meth, cap5, k2_note)), run_time=0.9)
        row = Group()
        for i in range(4):
            hm = st.heatmap(imgs[i], height=1.6)
            row.add(Group(hm, st.framed(hm)))
        row.arrange(RIGHT, buff=0.6).move_to([0, 0.9, 0])
        lbls = VGroup(*[st.M(rf"R_{i + 1},\ \mathbf t_{i + 1}", size=36, color=st.EST)
                        .next_to(row[i], DOWN, buff=0.18) for i in range(4)])
        qs = VGroup(*[st.M(r"?\,,\ ?", size=40, color=st.EST).move_to(lbls[i])
                      for i in range(4)])
        img_lab = VGroup(*[st.M(rf"I_{i + 1}", size=36).next_to(row[i], UP, buff=0.15)
                           for i in range(4)])
        inv = st.M(r"\{I_i\}\ \longmapsto\ \{\rho,\ R_i,\ \mathbf t_i\}", size=48)
        inv.move_to([0, -1.6, 0])
        joint = st.T("Poses and density must be estimated jointly", 32).next_to(inv, DOWN,
                                                                                 buff=0.35)
        ctf_note = st.T("CTF and noise must also be modelled (not shown here)", 26,
                        st.MUTED).next_to(joint, DOWN, buff=0.2)
        clean_note = st.T("noise-free images for clarity", 26, st.MUTED).next_to(row, UP,
                                                                                  buff=0.55)
        self.fix(row, lbls, qs, img_lab, inv, joint, ctf_note, clean_note, name="B6")
        self.play(FadeIn(row), FadeIn(img_lab), FadeIn(lbls), FadeIn(clean_note), run_time=0.9)
        self.wait(0.5)
        self.play(LaggedStart(*[Succession(FadeOut(lbls[i], shift=0.1 * UP, run_time=0.35),
                                           FadeIn(qs[i], shift=0.1 * UP, run_time=0.35))
                                for i in range(4)], lag_ratio=0.25), run_time=1.4)
        self.play(FadeIn(inv), run_time=0.9)
        self.play(FadeIn(joint), FadeIn(ctf_note), run_time=0.7)
        self.wait(3.0)

        # ============================================================ B7 takeaway
        card = st.takeaway_card(
            ["Orientation: viewing direction (2 DOF) + in-plane angle (1)",
             "An ideal projection's 2D FT is a central plane of the 3D FT",
             "Biased views leave a cone of frequencies unsampled",
             "With unknown poses, density and poses are found jointly"],
            math=r"\widehat{P_R\rho}(k_x,k_y)=\widehat{\rho}\left(R^{\top}(k_x,k_y,0)^{\top}"
                 r"\right)")
        card.move_to([0, -0.2, 0])
        self.fix(card, name="card")
        self.play(FadeOut(Group(row, qs, img_lab, inv, joint, ctf_note, clean_note)),
                  FadeOut(hdr6), run_time=0.7)
        self.play(FadeIn(card), run_time=0.7)
        self.wait(4.0)
