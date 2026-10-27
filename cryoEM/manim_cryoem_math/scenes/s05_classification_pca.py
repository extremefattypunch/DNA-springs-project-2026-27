"""Clip 05 — Statistical classification and PCA. Storyboard beats -> code sections (spec §3).

Beat 1 (0-9 s)    two abstract states that differ only in dot 5
Beat 2 (9-17 s)   noisy observations at one matched view; the mixed mean is ambiguous
Beat 3 (17-29 s)  PCA scores (PC1, PC2): a linear basis of variation, not a trajectory
Beat 4 (29-39 s)  2-component Gaussian mixture: posterior bars for ambiguous points
Beat 5 (39-48 s)  freeze: mixed mean vs class-conditioned means, mixture equation

Every number on screen comes from cached("classification", ne.classification).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
from manim import *
from helpers import style as st, numeric_examples as ne, noise_and_ctf as nc
from helpers import projections as pj
from helpers.cache import cached

N_PIX = 48           # image edge used by ne.classification
EXTENT = 2.0         # images span [-EXTENT, EXTENT)
VIEW = (10.0, 20.0, 0.0)   # ZYZ Euler angles (deg) of the fixed view in ne.classification


def img_point(img, p, n=N_PIX, extent=EXTENT):
    """Scene position of 2D coordinate p on a heatmap of an [x, y] array
    (pixel i is centred at (i - n//2) * h, see pc.grid_coords)."""
    h = 2 * extent / n
    fx = (p[0] / h + n // 2 + 0.5) / n
    fy = (p[1] / h + n // 2 + 0.5) / n
    return np.array([img.get_left()[0] + fx * img.width,
                     img.get_bottom()[1] + fy * img.height, 0.0])


def ring(center, r=0.24, color=st.ALT, dashed=True, width=3):
    c = Circle(radius=r, color=color, stroke_width=width).move_to(center)
    return DashedVMobject(c, num_dashes=14) if dashed else c


def fmt(v, nd=2):
    return f"{v:.{nd}f}"


def tfmt(v, nd=2):
    """fmt for Text (typographic minus)."""
    return fmt(v, nd).replace("-", "\u2212")


class Scene05ClassificationPCA(Scene):
    def construct(self):
        d = cached("classification", ne.classification)
        tag = st.illustration_label().set_z_index(10)
        self.add(tag)

        labels = d["labels"].astype(int)
        n_obs = len(labels)
        pA = d["pA"]
        snr = nc.snr_amplitude(d["cleanA"], float(d["sigma"]))

        # ============================================================ beat 1: two states
        title = st.title_card("Classification and PCA",
                              "How do heterogeneity and classification interact?")
        self.play(FadeIn(title, shift=0.1 * DOWN), run_time=1.0)

        scale = 1.12
        panels, dots5, labs = [], [], []
        for k, (name, P) in enumerate((("State A", d["proj_A"]), ("State B", d["proj_B"]))):
            sq = Square(side_length=3.2, stroke_color=st.GRID, stroke_width=2)
            sq.move_to([(-3.3 if k == 0 else 3.3), 0.15, 0])
            c = sq.get_center()
            ds = VGroup(*[Dot(c + scale * np.array([x, y, 0]), radius=0.11,
                              color=(st.ALT if j == 4 else st.TRUTH))
                          for j, (x, y) in enumerate(P)])
            lab = st.T(name, st.BODY_SIZE, weight="BOLD").next_to(sq, UP, buff=0.15)
            panels.append(VGroup(sq, ds[:4]))
            dots5.append(ds[4])
            labs.append(lab)
        x5 = [d["points_A"][4], d["points_B"][4]]
        coords = VGroup(*[
            st.M(r"\mathbf x_5^{(" + c + r")}=(" + ",\\,".join(fmt(v) for v in x) + ")",
                 size=36, color=st.ALT).next_to(panels[k][0], DOWN, buff=0.18)
            for k, (c, x) in enumerate(zip("AB", x5))])
        # the panels draw the projected 2D positions; check the view angles against the data
        R5 = pj.rotation_from_euler(*VIEW)
        assert np.allclose(pj.project_points(R5, d["points_A"]), d["proj_A"])
        assert np.allclose(pj.project_points(R5, d["points_B"]), d["proj_B"])
        a_, t_, p_ = (f"{v:.0f}" for v in VIEW)
        view_note = VGroup(
            st.T("panels: 2D projection along", st.SMALL_SIZE, st.MUTED),
            st.M(r"\mathbf v=R^{\top}\mathbf e_z,\quad R=R_{ZYZ}(" + a_ + r"^\circ,"
                 + t_ + r"^\circ," + p_ + r"^\circ)", size=36, color=st.MUTED),
        ).arrange(RIGHT, buff=0.15).next_to(coords, DOWN, buff=0.12)
        mid = st.T("only dot 5\ndiffers", st.BODY_SIZE - 4, st.MUTED).move_to([0, 0.15, 0])
        # in state B: ghost of the state-A site and an arrow to the new site
        ghost_pos = panels[1][0].get_center() + scale * np.array([*d["proj_A"][4], 0])
        ghost = ring(ghost_pos, r=0.13, width=2)
        move = Arrow(ghost_pos, dots5[1].get_center(), buff=0.16, color=st.ALT,
                     stroke_width=4, max_tip_length_to_length_ratio=0.25)
        eq1 = st.M(r"X_c=\{\mathbf x_1,\dots,\mathbf x_5^{(c)}\},\quad c\in\{A,B\}", size=44)
        eq1.to_edge(DOWN, buff=st.MARGIN + 0.05)
        assert view_note.get_bottom()[1] > eq1.get_top()[1] + 0.1

        self.play(*[FadeIn(VGroup(labs[k], panels[k])) for k in range(2)], run_time=1.2)
        self.play(FadeIn(dots5[0], scale=1.5), FadeIn(dots5[1], scale=1.5),
                  FadeIn(coords), FadeIn(view_note), run_time=1.0)
        self.play(Create(ghost), GrowArrow(move), FadeIn(mid), run_time=1.2)
        self.play(Write(eq1), run_time=1.4)
        self.wait(3.2)
        beat1 = VGroup(*labs, *panels, *dots5, coords, view_note, mid, ghost, move, eq1)

        # ============================================================ beat 2: mixed mean
        self.play(FadeOut(beat1), FadeOut(title), run_time=0.8)
        imgs = d["imgs12"]
        idx_A = [i for i in range(len(imgs)) if labels[i] == 0][:3]
        idx_B = [i for i in range(len(imgs)) if labels[i] == 1][:3]
        order = [idx_A[0], idx_B[0], idx_A[1], idx_B[1], idx_A[2], idx_B[2]]
        thumbs = Group()
        for i in order:
            hm = st.heatmap(imgs[i], height=1.3)
            fr = st.framed(hm, color=(st.TRUTH if labels[i] == 0 else st.ALT))
            fr.set_stroke(width=4)
            thumbs.add(Group(hm, fr))
        thumbs.arrange_in_grid(rows=2, cols=3, buff=0.18).move_to([-4.0, 0.2, 0])
        obs_lab = st.T(f"6 of {n_obs} observations, one view", st.BODY_SIZE - 4)
        obs_lab.next_to(thumbs, UP, buff=0.2).align_to(thumbs, LEFT)
        obs_snr = st.T(f"amplitude SNR = {snr:.2f}:\npeaks hidden in each image",
                       st.SMALL_SIZE, st.MUTED)
        obs_snr.next_to(thumbs, DOWN, buff=0.2).align_to(thumbs, LEFT)
        leg = st.legend([(st.TRUTH, "frame: true state A"), (st.ALT, "frame: true state B")])
        leg.next_to(obs_snr, DOWN, buff=0.15).align_to(thumbs, LEFT)
        assert leg.get_bottom()[1] > -4 + st.MARGIN

        means = [d["mixed"], d["class_mean"][0], d["class_mean"][1]]
        vmin = min(float(m.min()) for m in means)
        vmax = max(float(m.max()) for m in means)
        mix = st.heatmap(d["mixed"], height=3.1, vmin=vmin, vmax=vmax).move_to([1.6, 0.05, 0])
        mix_fr = st.framed(mix, color=st.OBS)
        mix_lab = st.M(r"\bar I=\tfrac{1}{N}\textstyle\sum_{i=1}^{N} I_i", size=40)
        mix_lab.next_to(mix, UP, buff=0.2)
        arr = Arrow([thumbs.get_right()[0] + 0.1, mix.get_y(), 0], mix.get_left() + 0.1 * LEFT,
                    buff=0.1,
                    color=st.OBS, stroke_width=5)
        sA, sB = img_point(mix, d["site_A"]), img_point(mix, d["site_B"])
        rA, rB = ring(sA), ring(sB)
        pk = d["peaks"]
        tA = st.M(r"\mathbf x_5^{(A)}", size=36, color=st.ALT).next_to(mix, RIGHT, buff=0.35)
        tA.set_y(sA[1])
        tB = st.M(r"\mathbf x_5^{(B)}", size=36, color=st.ALT).next_to(mix, RIGHT, buff=0.35)
        tB.set_y(sB[1])
        lA = DashedLine(rA.get_right(), tA.get_left(), color=st.ALT, stroke_width=2)
        lB = DashedLine(rB.get_right(), tB.get_left(), color=st.ALT, stroke_width=2)
        vals = VGroup(
            st.T(f"value at site A: {fmt(pk[2, 0])}   site B: {fmt(pk[2, 1])}", st.SMALL_SIZE),
            st.T(f"(noise-free single state: {fmt(pk[0, 0])})", st.SMALL_SIZE, st.MUTED),
        ).arrange(DOWN, buff=0.1).next_to(mix, DOWN, buff=0.2)
        head2 = st.T("Mixed observations → ambiguous mean", st.BODY_SIZE, st.TEXT,
                     weight="BOLD").to_edge(DOWN, buff=st.MARGIN)

        self.play(LaggedStart(*[FadeIn(t) for t in thumbs], lag_ratio=0.15),
                  FadeIn(obs_lab), run_time=1.6)
        self.play(FadeIn(obs_snr), FadeIn(leg), run_time=0.6)
        self.play(GrowArrow(arr), FadeIn(mix), Create(mix_fr), FadeIn(mix_lab), run_time=1.3)
        self.play(Create(rA), Create(rB), Create(lA), Create(lB), FadeIn(tA), FadeIn(tB),
                  run_time=1.0)
        self.play(FadeIn(vals), FadeIn(head2, shift=0.1 * UP), run_time=0.8)
        self.wait(2.2)
        beat2 = Group(thumbs, obs_lab, obs_snr, leg, arr, mix, mix_fr, mix_lab,
                      rA, rB, tA, tB, lA, lB, vals, head2)

        # ============================================================ beat 3: PCA
        self.play(FadeOut(beat2), run_time=0.8)
        S = d["scores"]
        ax = Axes(x_range=[-7, 5, 2], y_range=[-6, 6, 2], x_length=6.0, y_length=5.0,
                  tips=False, axis_config=dict(color=st.GRID, stroke_width=2,
                                               include_ticks=True, tick_size=0.05))
        ax.move_to([-3.45, -0.5, 0])
        xl = st.T("PC1", st.BODY_SIZE - 4, st.MUTED).next_to(ax.x_axis, RIGHT, buff=0.15)
        # beside the top of the y-axis, so neither scatter legend (leg3, leg4) covers it
        yl = st.T("PC2", st.BODY_SIZE - 4, st.MUTED).next_to(ax.c2p(0, 6), LEFT, buff=0.12)
        yl.align_to(ax.c2p(0, 6), UP).shift(0.1 * DOWN)
        pts = VGroup(*[Dot(ax.c2p(*s), radius=0.045, color=st.MUTED, fill_opacity=0.85)
                       for s in S])
        head3 = VGroup(st.T("Each observation → a point", st.BODY_SIZE, weight="BOLD"),
                       st.M(r"(z_{i1},\,z_{i2})", size=40)).arrange(RIGHT, buff=0.2)
        head3.to_corner(UL, buff=st.MARGIN)
        # exact definition of the scores computed by ne.classification (orthonormal u_k)
        pca_eq = st.M(r"z_{ik}=\mathbf u_k^{\top}(I_i-\bar I),\quad k=1,2", size=40)
        pca_eq.move_to([3.85, 2.2, 0])
        # signed basis image: green = positive (site A), purple = negative (site B)
        u1 = st.heatmap(d["pc_images"][0], height=2.0, symmetric=True,
                        stops=[(0.0, st.ALT), (0.5, "#F2F2F2"), (1.0, st.TRUTH)])
        u1.move_to([1.85, 0.45, 0])
        # the SVD sign of u_1 is arbitrary: the labels below need u_1 > 0 at site A,
        # u_1 < 0 at site B and true-A scores at positive PC1 (nearest-pixel lookup)
        gx = (np.arange(N_PIX) - N_PIX // 2) * (2 * EXTENT / N_PIX)
        u1v = lambda p: d["pc_images"][0][np.argmin(abs(gx - p[0])), np.argmin(abs(gx - p[1]))]
        assert u1v(d["site_A"]) > 0 > u1v(d["site_B"]) and S[labels == 0, 0].mean() > 0
        u1_fr = st.framed(u1)
        u1_lab = VGroup(st.M(r"\mathbf u_1", size=40),
                        st.T("PC1 basis image", st.SMALL_SIZE, st.MUTED),
                        st.T("green +, purple −", st.SMALL_SIZE, st.MUTED)
                        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        u1_lab.next_to(u1, RIGHT, buff=0.3)
        ex = d["explained"]
        expl = st.T(f"explained variance: PC1 {100 * ex[0]:.2f} %;\n"
                    f"PC2–PC5 {100 * min(ex[1:5]):.2f}–{100 * max(ex[1:5]):.2f} % (noise level)",
                    st.SMALL_SIZE, st.MUTED)
        expl.next_to(u1, DOWN, buff=0.25).align_to(u1, LEFT)
        assert expl.get_right()[0] <= 7.11 - st.MARGIN - 0.05, expl.get_right()[0]
        caution3 = st.T("PCA: a linear basis of variation,\n"
                        "not a trajectory, not a state label", st.BODY_SIZE - 4, st.TEXT)
        caution3.next_to(expl, DOWN, buff=0.3).align_to(u1, LEFT)
        st.fit_width(caution3, 7.11 - st.MARGIN - caution3.get_left()[0])
        caution3.align_to(u1, LEFT)
        leg3 = st.legend([(st.TRUTH, "true state A"), (st.ALT, "true state B")])
        leg3.next_to(ax, UP, buff=0.05).align_to(ax, LEFT)
        leg3.add_background_rectangle(color=st.BG, opacity=1, buff=0.05).set_z_index(2)

        self.play(Create(ax), FadeIn(xl), FadeIn(yl), FadeIn(head3), run_time=1.2)
        self.play(LaggedStart(*[FadeIn(p) for p in pts], lag_ratio=0.01), run_time=1.8)
        self.play(Write(pca_eq), run_time=1.2)
        self.play(FadeIn(u1), Create(u1_fr), FadeIn(u1_lab), run_time=1.0)
        self.play(FadeIn(expl), run_time=0.7)
        self.play(*[p.animate.set_color(st.TRUTH if labels[i] == 0 else st.ALT)
                    for i, p in enumerate(pts)], FadeIn(leg3), run_time=1.3)
        self.play(FadeIn(caution3), run_time=0.8)
        self.wait(2.6)

        # ============================================================ beat 4: probabilities
        self.play(FadeOut(Group(pca_eq, u1, u1_fr, u1_lab, expl, caution3)),
                  FadeOut(head3, shift=0.1 * UP), run_time=0.7)
        head4 = st.T("Mixture of K = 2 Gaussians, fit by EM", st.BODY_SIZE,
                     weight="BOLD").to_corner(UL, buff=st.MARGIN)

        def ellipse(m, C, col, k=2.0):
            w, V = np.linalg.eigh(C)
            f = lambda t: ax.c2p(*(m + k * (V @ (np.sqrt(w) * np.array([np.cos(t), np.sin(t)])))))
            return ParametricFunction(f, t_range=[0, TAU], color=col, stroke_width=4)
        ells = VGroup(ellipse(d["em_means"][0], d["em_covs"][0], st.TRUTH),
                      ellipse(d["em_means"][1], d["em_covs"][1], st.ALT))
        ells.set_stroke(st.BG, width=9, background=True)
        cents = VGroup(*[Cross(scale_factor=0.12, stroke_color=c, stroke_width=5)
                         .move_to(ax.c2p(*d["em_means"][k]))
                         for k, c in enumerate((st.TRUTH, st.ALT))])
        pis = d["em_pis"]

        post = st.M(r"p(c\mid I_i)=\frac{\pi_c\,p(I_i\mid c)}{\sum_{c'}\pi_{c'}\,p(I_i\mid c')}",
                    size=38).move_to([3.75, 2.15, 0])
        # the mixture is fitted to the PCA scores z_i, so the likelihood is defined on them
        lik_def = st.M(r"p(I_i\mid c):=\mathcal N(\mathbf z_i;\,\boldsymbol\mu_c,\Sigma_c)",
                       size=36).next_to(post, DOWN, buff=0.15)
        # three ambiguous observations (0.25 < p(A|I) < 0.75), spread in p and on screen
        amb = [int(j) for j in d["ambiguous"]]
        pick = []
        for target in (0.7, 0.5, 0.33):
            ok = [j for j in amb if j not in pick and all(
                np.linalg.norm(pts[j].get_center() - pts[q].get_center()) > 0.75 for q in pick)]
            pick.append(min(ok, key=lambda j: abs(pA[j] - target)))
        rings = VGroup(*[Circle(radius=0.14, color=st.TEXT, stroke_width=3)
                         .move_to(pts[j].get_center()) for j in pick])
        rings.set_stroke(st.BG, width=7, background=True)
        # ring 1 sits where the two 2σ contours cross: its label goes right, inside the green
        # ellipse (clear of both contours); rings 2 and 3 are labelled on the left
        rnums = VGroup(*[st.T(str(k + 1), st.SMALL_SIZE)
                         .next_to(rings[k], RIGHT if k == 0 else LEFT, buff=0.06)
                         .set_stroke(st.BG, width=6, background=True)
                         .set_z_index(2) for k in range(3)])
        bar_w, x_lab, row_h = 2.0, 2.3, 0.42       # fixed columns so all bars line up
        rows = VGroup()
        y0 = lik_def.get_bottom()[1] - 0.35
        for k, j in enumerate(pick):
            y = y0 - k * (2 * row_h + 0.08)
            name = st.T(str(k + 1), st.BODY_SIZE - 4).move_to([1.3, y - row_h / 2, 0])
            grp = VGroup(name)
            for c, (col, pv) in enumerate(((st.TRUTH, pA[j]), (st.ALT, 1 - pA[j]))):
                yc = y - c * row_h
                lab = st.M(r"p(" + "AB"[c] + r"\mid I_i)", size=36, color=col)
                lab.move_to([x_lab, yc, 0])
                back = Rectangle(width=bar_w, height=0.24, stroke_color=st.GRID, stroke_width=1.5)
                back.move_to([3.15 + bar_w / 2, yc, 0])
                bar = Rectangle(width=max(bar_w * pv, 0.01), height=0.24, stroke_width=0,
                                fill_color=col, fill_opacity=0.9)
                bar.align_to(back, LEFT).set_y(yc)
                val = st.T(fmt(pv), st.SMALL_SIZE).next_to(back, RIGHT, buff=0.12)
                grp.add(lab, back, bar, val)
            rows.add(grp)
        acc = float(d["accuracy"])
        true_frac_A = float(np.mean(labels == 0))
        info = VGroup(
            st.M(r"\hat\pi_A=" + fmt(pis[0]) + r",\quad \hat\pi_B=" + fmt(pis[1]), size=36),
            st.T(f"true fraction A in this set: {fmt(true_frac_A)}", st.SMALL_SIZE, st.MUTED),
            st.T(f"assigned class = true state: {100 * acc:.0f} %", st.SMALL_SIZE, st.MUTED),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        info.next_to(rows, DOWN, buff=0.2).align_to(rows, LEFT)
        caution4 = st.caption("Schematic: a mixture fit on two PCA coordinates, "
                              "not any specific package's classifier")
        assert info.get_bottom()[1] > caution4.get_top()[1] + 0.1
        assert info.get_right()[0] <= 7.11 - st.MARGIN, info.get_right()[0]
        # scatter legend replacing the true-state legend: dot colour key and ellipse meaning
        key = Rectangle(width=1.6, height=0.18, stroke_width=0)
        key.set_fill(color=[st.ALT, st.TRUTH], opacity=1).set_sheen_direction(RIGHT)
        leg4 = VGroup(
            VGroup(st.T("dot colour:", st.SMALL_SIZE, st.MUTED),
                   st.M(r"p(A\mid I_i)", size=36, color=st.MUTED),
                   st.T("0", st.SMALL_SIZE, st.MUTED), key,
                   st.T("1", st.SMALL_SIZE, st.MUTED)).arrange(RIGHT, buff=0.12),
            VGroup(st.T("ellipses: 2σ contours of fitted", st.SMALL_SIZE, st.MUTED),
                   st.M(r"\mathcal N(\boldsymbol\mu_c,\Sigma_c)", size=36, color=st.MUTED)
                   ).arrange(RIGHT, buff=0.12),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        leg4.next_to(ax, UP, buff=0.05).align_to(ax, LEFT)
        leg4.add_background_rectangle(color=st.BG, opacity=1, buff=0.05).set_z_index(2)
        assert yl.get_top()[1] < leg4.get_bottom()[1] - 0.05
        assert yl.get_left()[0] > leg3.get_right()[0] + 0.05

        self.play(FadeIn(head4), Create(ells), FadeIn(cents), run_time=1.4)
        self.play(Write(post), FadeIn(lik_def), run_time=1.2)
        self.play(Create(rings), FadeIn(rnums), run_time=0.8)
        self.play(LaggedStart(*[FadeIn(r, shift=0.1 * LEFT) for r in rows], lag_ratio=0.3),
                  run_time=1.5)
        self.play(*[p.animate.set_color(interpolate_color(ManimColor(st.ALT),
                                                          ManimColor(st.TRUTH), float(pA[i])))
                    for i, p in enumerate(pts)],
                  Succession(FadeOut(leg3), FadeIn(leg4)), FadeIn(info[:2]), run_time=1.6)
        self.play(FadeIn(info[2]), FadeIn(caution4), run_time=0.8)
        self.wait(2.4)
        beat4 = VGroup(ax, xl, yl, pts, head4, ells, cents, post, lik_def, rings, rnums, rows,
                       info, leg4, caution4)

        # ============================================================ beat 5: freeze
        self.play(FadeOut(beat4), run_time=0.8)
        names = ["mixed mean", "class mean A", "class mean B"]
        cols = [st.OBS, st.TRUTH, st.ALT]
        tiles = Group()
        for k, (m, nm, col) in enumerate(zip(means, names, cols)):
            hm = st.heatmap(m, height=1.7, vmin=vmin, vmax=vmax)
            fr = st.framed(hm, color=col).set_stroke(width=3)
            rr = VGroup(ring(img_point(hm, d["site_A"]), r=0.15),
                        ring(img_point(hm, d["site_B"]), r=0.15))
            row = 2 + k
            cap = VGroup(st.T(nm, st.BODY_SIZE - 4, col, weight="BOLD"),
                         st.T(f"site A {tfmt(pk[row, 0])}   site B {tfmt(pk[row, 1])}",
                              st.SMALL_SIZE)).arrange(DOWN, buff=0.08)
            cap.next_to(hm, DOWN, buff=0.12)
            tiles.add(Group(hm, fr, rr, cap))
        tiles.arrange(RIGHT, buff=1.0, aligned_edge=UP)
        # captions (~3.7 units) are wider than the heatmaps: centre the row so it keeps the margins
        tiles.move_to([0, 0, 0]).align_to([0, 2.9, 0], UP)
        scale_note = st.T("common grey scale; dashed rings: dot-5 sites A (upper) and B (lower)",
                          st.SMALL_SIZE, st.MUTED)
        scale_note.next_to(tiles, DOWN, buff=0.1)
        assert tiles.get_left()[0] >= -7.11 + st.MARGIN and tiles.get_right()[0] <= 7.11 - st.MARGIN
        # class means are responsibility-weighted (soft) averages, as in ne.classification
        cm_eq = st.M(r"\bar I_c=\frac{\sum_i p(c\mid I_i)\,I_i}{\sum_i p(c\mid I_i)}", size=36)
        # shown below the two class-mean tiles, in the area the card takes later
        cm_eq.next_to(scale_note, DOWN, buff=0.4).set_x(0.5 * (tiles[1].get_x() + tiles[2].get_x()))
        assert cm_eq.get_right()[0] <= 7.11 - st.MARGIN
        card = st.takeaway_card(
            ["Mixing states smears dot 5; class means keep it.",
             "Fractions π of this selected set are not",
             "automatically equilibrium probabilities."],
            math=r"p(I_i)={\textstyle\sum_{c=1}^{K}}\pi_c\,p(I_i\mid c),\qquad "
                 r"{\textstyle\sum_{c=1}^{K}}\pi_c=1")
        card.to_edge(DOWN, buff=st.MARGIN)
        # the scale note stays through the freeze, so it must clear the card
        assert scale_note.get_bottom()[1] > card.get_top()[1] + 0.1

        self.play(FadeIn(tiles[0]), FadeIn(scale_note), run_time=0.8)
        self.play(FadeIn(tiles[1]), FadeIn(tiles[2]), FadeIn(cm_eq), run_time=1.0)
        self.wait(2.4)
        # one complex equation at a time: the soft-mean formula leaves before the card
        self.play(FadeOut(cm_eq), run_time=0.4)
        self.play(FadeIn(card, shift=0.1 * UP), run_time=1.0)
        self.wait(4.0)
