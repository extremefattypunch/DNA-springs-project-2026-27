"""Numerical truth checks (spec section 5, items 1-8). Run: pytest -m "not render" -q"""
import numpy as np
import pytest

from helpers import fourier as fo
from helpers import noise_and_ctf as nc
from helpers import numeric_examples as ne
from helpers import point_constellation as pc
from helpers import projections as pj

ANGLES = [(0, 0, 0), (40, 70, 20), (-60, 110, 75), (150, 30, -100)]


def corr(a, b):
    return np.corrcoef(np.ravel(a), np.ravel(b))[0, 1]


# 1. SO(3) -----------------------------------------------------------------------
@pytest.mark.parametrize("ang", ANGLES)
def test_rotation_is_proper(ang):
    R = pj.rotation_from_euler(*ang)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(R), 1.0)


def test_random_rotations_proper():
    Rs = pj.random_rotations(50, np.random.default_rng(0))
    assert np.allclose(np.einsum("nij,nik->njk", Rs, Rs), np.eye(3), atol=1e-10)
    assert np.allclose(np.linalg.det(Rs), 1.0)


def test_viewing_direction_differs_from_inplane_rotation():
    R0 = pj.rotation_from_euler(30, 50, 0)
    R_spin = pj.rotation_from_euler(30, 50, 70)        # psi only
    R_tilt = pj.rotation_from_euler(30, 80, 0)
    assert np.allclose(pj.viewing_direction(R0), pj.viewing_direction(R_spin))
    assert not np.allclose(pj.viewing_direction(R0), pj.viewing_direction(R_tilt))
    X = pc.points()
    # in-plane spin = 2D rotation of the projected dots; a tilt is not
    p0, ps, pt = (pj.project_points(R, X) for R in (R0, R_spin, R_tilt))
    assert np.allclose(pj.transform_2d(p0, 70), ps)
    dist = lambda P: np.sort(np.linalg.norm(P[:, None] - P[None], axis=-1).ravel())
    assert np.allclose(dist(p0), dist(ps)) and not np.allclose(dist(p0), dist(pt), atol=1e-3)


# 2. projection ------------------------------------------------------------------
@pytest.mark.parametrize("ang", ANGLES)
def test_projected_dots_match_direct_coordinates(ang):
    R = pj.rotation_from_euler(*ang)
    X = pc.points()
    direct = np.array([(R @ x)[:2] for x in X])
    assert np.allclose(pj.project_points(R, X), direct)


@pytest.mark.parametrize("ang", ANGLES)
def test_analytic_image_matches_numerical_line_integral(ang):
    R = pj.rotation_from_euler(*ang)
    n, ext = pc.FIELD_N, pc.FIELD_EXTENT
    vol = pc.field(n=n, extent=ext)
    assert corr(pj.project_image(R, n=n, extent=ext), pj.project_field(vol, R, ext)) > 0.995


def test_image_peaks_sit_on_projected_dots():
    R = pj.rotation_from_euler(*ANGLES[1])
    img = pj.project_image(R, n=128)
    ax = pc.grid_coords(128, pc.IMAGE_EXTENT)
    for p in pj.project_points(R, pc.points()):
        i, j = np.argmin(np.abs(ax - p[0])), np.argmin(np.abs(ax - p[1]))
        assert img[i, j] > 0.5 * img.max() * 0.3


def test_constellation_is_asymmetric():
    p = [pj.project_image(pj.rotation_from_euler(*a)) for a in ANGLES]
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            assert corr(p[i], p[j]) < 0.9
    assert corr(p[1], pj.rotate_image(p[1], 180)) < 0.8


# 3. Fourier slice theorem ----------------------------------------------------------
def test_fourier_slice_uses_R_transpose():
    d = ne.orientation_fourier()
    assert np.all(d["corr_RT"] > 0.995), d["corr_RT"]
    assert np.all(d["corr_R"] < 0.8), d["corr_R"]


def test_fourier_slice_is_complex_not_magnitude_only():
    """Shifting a projection changes the phases of its FT but not |FT|."""
    R = pj.rotation_from_euler(*ANGLES[1])
    a = fo.fft2c(pj.project_image(R))
    b = fo.fft2c(pj.project_image(R, t=(0.3, -0.2)))
    m = np.abs(a) > 1e-3 * np.abs(a).max()
    assert np.allclose(np.abs(a), np.abs(b), rtol=0, atol=0.03 * np.abs(a).max())
    assert np.abs(np.angle(a[m] * np.conj(b[m]))).mean() > 0.3


def test_reconstruction_from_known_poses_and_coverage():
    rng = np.random.default_rng(1)
    n, ext = 48, pc.FIELD_EXTENT
    h = pc.spacing(n, ext)
    vol = pc.field(n=n, extent=ext)
    Rs = pj.random_rotations(100, rng)
    imgs = np.stack([pj.project_image(R, n=n, extent=ext) / h for R in Rs])
    rec = fo.reconstruct_direct_fourier(imgs, Rs)
    assert corr(rec, vol) > 0.98
    d = ne.orientation_fourier()
    assert d["cov_biased"] < 0.5 * d["cov_uniform"]
    assert d["err_uniform"][-1] < 0.05 and d["err_uniform"][0] > 0.3       # 1 view: depth lost
    assert d["err_biased"][-1] > d["err_uniform"][-1]


# 4. CTF -----------------------------------------------------------------------------
def test_ctf_has_zero_and_sign_change():
    k = np.linspace(0, 8, 4000)
    c = nc.ctf_1d(k)
    assert np.sum(np.diff(np.sign(c)) != 0) >= 3
    assert c.max() > 0.95 and c.min() < -0.95
    z = nc.zeros(8.0)
    assert np.allclose(nc.ctf_1d(z), 0, atol=1e-6)


def test_ctf_is_not_a_blur():
    d = ne.forward_model()
    assert d["ctf_at_bands"][0] < -0.5 and d["ctf_at_bands"][2] > 0.5
    assert abs(d["ctf_at_bands"][1]) < 1e-6
    # the band in lobe 1 comes out sign-reversed, the band at the zero vanishes
    assert corr(d["bands_before"][0], d["bands_after"][0]) < -0.9
    assert np.std(d["bands_after"][1]) < 0.15 * np.std(d["bands_before"][1])


# 5. noise averaging -----------------------------------------------------------------------
def test_noise_std_scales_as_inverse_sqrt_n():
    rng = np.random.default_rng(3)
    for n_avg in (1, 4, 16, 64, 256):
        m = rng.normal(0, 2.0, (n_avg, 4000)).mean(0)
        assert np.isclose(m.std(), 2.0 / np.sqrt(n_avg), rtol=0.05)


def test_alignment_example_sqrt_n_and_misalignment():
    d = ne.alignment()
    slope = np.polyfit(np.log(d["ns"]), np.log(d["snr_measured"]), 1)[0]
    assert abs(slope - 0.5) < 0.03
    assert np.allclose(d["snr_measured"], d["snr_theory"], rtol=0.1)
    assert d["corr_aligned"][-1] > d["corr_unaligned"][-1] + 0.3


# 6. refinement -----------------------------------------------------------------------------
def test_refinement_values_are_computed_and_decrease():
    d = ne.refinement()
    assert np.all(np.diff(d["loss"]) <= 1e-6)
    assert np.allclose(d["loss"], d["loss_data"] + d["loss_reg"])
    # the stored data term really is sum ||I - P[X]||^2 for the stored estimates
    i = d["track"]
    assert np.isclose(np.sum(d["track_res"][-1] ** 2),
                      np.sum((d["imgs"][i] - d["track_pred"][-1]) ** 2))
    assert d["pose_err"][-1].mean() < 0.3 * d["pose_err"][0].mean()
    assert d["track_res_norm"][-1] < d["track_res_norm"][0]
    assert d["loss_data"][-1] > 0.95 * d["noise_floor"]          # does not fit the noise
    assert d["reg_data_res_ls"] < d["reg_data_res_rg"]
    assert d["reg_true_err_ls"] > d["reg_true_err_rg"]


# 7. classification -------------------------------------------------------------------------
def test_mixed_mean_is_ambiguous_and_class_means_are_not():
    d = ne.classification()
    cA, cB = d["cleanA"], d["cleanB"]
    mixed, (m0, m1) = d["mixed"], d["class_mean"]
    pk = d["peaks"]                                 # rows: A, B, mixed, class0, class1
    assert 0.25 < pk[2, 0] / pk[0, 0] < 0.8 and 0.25 < pk[2, 1] / pk[1, 1] < 0.8
    assert pk[3, 0] > 0.85 * pk[0, 0] and pk[4, 1] > 0.75 * pk[1, 1]
    assert corr(m0, cA) > corr(mixed, cA) and corr(m1, cB) > corr(mixed, cB)
    assert d["accuracy"] > 0.85 and len(d["ambiguous"]) > 0


# 8. FSC ------------------------------------------------------------------------------------
def test_fsc_identical_is_one_and_noise_is_zero():
    v = pc.field(n=48)
    s, f, cnt = fo.fsc(v, v)
    assert np.all(f[cnt > 0] > 0.999)
    rng = np.random.default_rng(4)
    s, f, cnt = fo.fsc(rng.standard_normal((48,) * 3), rng.standard_normal((48,) * 3))
    assert abs(f[4:].mean()) < 0.02 and np.all(np.abs(f[4:]) < 5 / np.sqrt(cnt[4:]))


def test_fsc_empty_shells_are_safe():
    a = np.zeros((32,) * 3)
    s, f, _ = fo.fsc(a, a)
    assert np.all(np.isfinite(f)) and np.all(f == 0)


def test_fsc_example_distinguishes_global_local_and_model():
    d = ne.fsc_resolution()
    assert 0.0 < d["s_half"] < d["nyquist"]
    assert d["d_region2"] > 1.3 * d["d_region1"]                 # local variation
    assert np.isclose(d["d_half"], 1 / d["s_half"])
    # independent half-map noise: FSC falls to ~0 at high frequency
    assert abs(d["fsc_half"][-8:].mean()) < 0.1
