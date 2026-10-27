"""Numerical truth checks (plan section 4, items 1-7). Run: pytest tests -q"""
import numpy as np
import pytest

from helpers import fourier as fo
from helpers import optics as op
from helpers import projections as pj
from helpers import synthetic_phantom as sp

N, V = 64, 2.0
ANGLES = [(0, 0, 0), (40, 70, 20), (-60, 110, 75), (150, 30, -100)]


def corr(a, b):
    return np.corrcoef(np.ravel(a), np.ravel(b))[0, 1]


@pytest.fixture(scope="module")
def vol():
    return sp.ground_truth(N, V, atoms=False)


# 1. rotations ---------------------------------------------------------------
@pytest.mark.parametrize("ang", ANGLES)
def test_rotation_is_proper(ang):
    R = pj.rotation_from_euler(*ang)
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(R), 1.0)


def test_random_rotations_proper():
    Rs = pj.random_rotations(50, np.random.default_rng(0))
    assert np.allclose(np.einsum("nij,nik->njk", Rs, Rs), np.eye(3), atol=1e-10)
    assert np.allclose(np.linalg.det(Rs), 1.0)


def test_inplane_spin_keeps_viewing_direction_but_tilt_changes_it():
    R0 = pj.rotation_from_euler(30, 50, 0)
    R_spin = pj.rotation_from_euler(30, 50, 70)     # psi only: in-plane rotation
    R_tilt = pj.rotation_from_euler(30, 80, 0)
    assert np.allclose(pj.viewing_direction(R0), pj.viewing_direction(R_spin))
    assert not np.allclose(pj.viewing_direction(R0), pj.viewing_direction(R_tilt))
    # an in-plane spin is a 2D rotation of the projection image
    p0 = pj.project_analytic(R0, "A", N, V, atoms=False)
    ps = pj.project_analytic(R_spin, "A", N, V, atoms=False)
    assert corr(pj.rotate_image(p0, 70), ps) > 0.995
    assert corr(p0, ps) < 0.9


# 2. projections ---------------------------------------------------------------
def test_single_gaussian_projection_matches_analytic_formula():
    b = sp.blob_params("A")[0]
    R = pj.rotation_from_euler(10, 40, 25)
    X, Y = pj.image_coords(N, V)
    mu, S = R @ b.center, R @ b.cov @ R.T
    S2 = S[:2, :2]
    mass = b.amplitude * (2 * np.pi) ** 1.5 * np.prod(b.sigmas)
    d = np.stack([X - mu[0], Y - mu[1]], -1)
    q = np.einsum("...i,ij,...j->...", d, np.linalg.inv(S2), d)
    expected = mass / (2 * np.pi * np.sqrt(np.linalg.det(S2))) * np.exp(-0.5 * q)
    # numerical line integral of the same single blob on a fine z grid
    zs = np.linspace(-80, 80, 801)
    num = np.zeros_like(X)
    for z in zs:
        pts = np.stack([X, Y, np.full_like(X, z)], -1) @ R      # R^T r
        num += sp._eval_blob(b, pts[..., 0], pts[..., 1], pts[..., 2])
    num *= zs[1] - zs[0]
    assert np.allclose(num, expected, rtol=1e-3, atol=1e-6 * expected.max())


@pytest.mark.parametrize("ang", ANGLES)
def test_analytic_projection_matches_voxel_projection(vol, ang):
    R = pj.rotation_from_euler(*ang)
    assert corr(pj.project_analytic(R, "A", N, V, atoms=False), pj.project_voxel(vol, R, V)) > 0.999


def test_phantom_is_not_symmetric():
    """Different views must look different (so rotations are visually distinguishable)."""
    p = [pj.project_analytic(pj.rotation_from_euler(*a), "A", N, V, atoms=False) for a in ANGLES]
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            assert corr(p[i], p[j]) < 0.9
    # and a 180 degree in-plane flip is not a symmetry
    assert corr(p[1], pj.rotate_image(p[1], 180)) < 0.8


# 3. Fourier slice theorem ------------------------------------------------------
@pytest.mark.parametrize("ang", ANGLES[1:])
def test_fourier_slice_uses_R_transpose(vol, ang):
    R = pj.rotation_from_euler(*ang)
    F3 = fo.fft3c(fo.pad_volume(vol, 2))
    p = pj.project_analytic(R, "A", N, V, atoms=False) / V
    P2 = fo.fft2c(np.pad(p, N // 2))
    S_rt = fo.central_slice(F3, R)          # samples F3 at R^T k
    S_r = fo.central_slice(F3, R.T)         # the other convention
    m = fo.shell_index(P2.shape) < N // 2   # compare where interpolation is accurate
    c_rt = np.abs(np.vdot(P2[m], S_rt[m])) / (np.linalg.norm(P2[m]) * np.linalg.norm(S_rt[m]))
    c_r = np.abs(np.vdot(P2[m], S_r[m])) / (np.linalg.norm(P2[m]) * np.linalg.norm(S_r[m]))
    assert c_rt > 0.995, c_rt
    assert c_r < 0.8, c_r
    # and the real-space round trip through the slice reproduces the projection
    assert corr(fo.project_fourier(F3, R, N), p) > 0.999


def test_direct_fourier_reconstruction_recovers_volume(vol):
    rng = np.random.default_rng(1)
    Rs = pj.random_rotations(120, rng)
    imgs = np.stack([pj.project_analytic(R, "A", N, V, atoms=False) / V for R in Rs])
    rec = fo.reconstruct_direct_fourier(imgs, Rs)
    assert corr(rec, vol) > 0.99
    few = fo.reconstruct_direct_fourier(imgs[:3], Rs[:3])
    assert corr(few, vol) < corr(rec, vol) - 0.1


def test_biased_orientations_cover_less_fourier_space():
    rng = np.random.default_rng(2)
    u = fo.sampling_coverage(pj.random_rotations(200, rng), 32).mean()
    b = fo.sampling_coverage(pj.biased_rotations(200, rng, 12), 32).mean()
    assert b < 0.6 * u


# 4. CTF -------------------------------------------------------------------------
@pytest.mark.parametrize("df", [0.6, 1.0, 1.6, 2.4])
def test_ctf_oscillates_crosses_zero_and_reverses_sign(df):
    k = np.linspace(0, 1 / 3, 4000)
    c = op.ctf_1d(k, df)
    assert c[0] < 0 and abs(c[0]) < 0.1                  # weak (amplitude) contrast at k=0
    flips = np.sum(np.diff(np.sign(c)) != 0)
    assert flips >= 4
    assert c.max() > 0.99 and c.min() < -0.99            # full contrast reversal
    z = op.first_zero(df)
    assert np.isclose(op.ctf_1d(np.array([z]), df)[0], 0, atol=1e-3)


def test_ctf_first_zero_moves_to_lower_frequency_with_more_defocus():
    zs = [op.first_zero(d) for d in (0.6, 1.0, 1.6, 2.4)]
    assert all(a > b for a, b in zip(zs, zs[1:]))


def test_wavelength_300kv():
    assert np.isclose(op.electron_wavelength(300), 0.01969, atol=2e-5)


def test_apply_ctf_is_not_a_blur():
    """A CTF-filtered image contains negated components (contrast reversal)."""
    p = pj.project_analytic(pj.rotation_from_euler(30, 55, 15), "A", 128, 1.5)
    c = op.apply_ctf(p, 1.5, defocus_um=1.0)
    F, G = np.fft.fft2(p), np.fft.fft2(c)
    ratio = np.real(G / np.where(np.abs(F) > 1e-6 * np.abs(F).max(), F, np.inf))
    assert (ratio > 0.5).any() and (ratio < -0.5).any()


# 5. noise averaging ----------------------------------------------------------------
def test_noise_std_of_mean_falls_as_one_over_sqrt_n():
    rng = np.random.default_rng(3)
    sigma = 2.0
    for n_avg in (1, 4, 16, 64, 256):
        means = rng.normal(0, sigma, (n_avg, 4000)).mean(0)
        assert np.isclose(means.std(), sigma / np.sqrt(n_avg), rtol=0.05)


def test_alignment_demo_follows_sqrt_n():
    from helpers.demos import alignment
    d = alignment(n=48, voxel=3.2, n_particles=64, ns=(1, 4, 16, 64))
    slope = np.polyfit(np.log(d["ns"]), np.log(d["snr_measured"]), 1)[0]
    assert abs(slope - 0.5) < 0.05
    # aligned averages beat unaligned ones against the true signal
    S = d["signal"]
    assert corr(d["aligned_mean"][-1], S) > corr(d["unaligned_mean"][-1], S) + 0.2


# 6. FSC ------------------------------------------------------------------------------
def test_fsc_identical_maps_is_one(vol):
    s, f, cnt = fo.fsc(vol, vol, V)
    assert np.all(f[cnt > 0] > 0.999)


def test_fsc_independent_noise_is_near_zero():
    rng = np.random.default_rng(4)
    s, f, cnt = fo.fsc(rng.standard_normal((48,) * 3), rng.standard_normal((48,) * 3))
    assert abs(f[4:].mean()) < 0.02
    assert np.all(np.abs(f[4:]) < 5 / np.sqrt(cnt[4:]))


def test_fsc_no_division_by_zero():
    a = np.zeros((32,) * 3)
    s, f, cnt = fo.fsc(a, a)
    assert np.all(np.isfinite(f))


def test_shaped_noise_hits_target_crossing():
    rng = np.random.default_rng(5)
    v = sp.ground_truth(64, 1.0)
    tf = lambda s: fo.target_fsc(s, 1 / 3.2)
    a = v + fo.shaped_noise(v, tf, rng, 1.0)
    b = v + fo.shaped_noise(v, tf, rng, 1.0)
    s, f, _ = fo.fsc(a, b, 1.0)
    assert abs(1 / fo.crossing(s, f, 0.143) - 3.2) < 0.35


# 7. heterogeneity -----------------------------------------------------------------------
def test_mixed_state_average_reduces_appendage_contrast():
    vA = sp.conformer("A", N, V, atoms=False)
    vB = sp.conformer("B", N, V, atoms=False)
    mask = sp.appendage_mask("A", N, V)
    mixed = 0.5 * (vA + vB)
    assert mixed[mask].mean() < 0.8 * vA[mask].mean()
    # scaffold far from both appendage positions is unchanged
    far = ~(sp.appendage_mask("A", N, V, 3.0) | sp.appendage_mask("B", N, V, 3.0))
    assert np.allclose(mixed[far], vA[far], rtol=0, atol=1e-2 * vA.max())   # Gaussian tails only


def test_classification_demo_separates_states_softly():
    from helpers.demos import classification
    d = classification(n=48, voxel=2.6, n_particles=200)
    assert d["accuracy"] > 0.8
    assert d["contrast_classes"][0] > d["contrast_consensus"]     # class A sharper
    assert len(d["ambiguous"]) > 0                                 # not all hard-assigned
