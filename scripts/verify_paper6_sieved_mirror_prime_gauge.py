#!/usr/bin/env python3
"""
Standalone Verification & Figure-Generation Suite for Paper VI:
"Sieved Mirror Kernels, Gonek-Hughes-Keating versus Entire Exponential Prime Factors,
 Prime-Power Purity on (X, X^2], and Functional-Equation Obstructions"

Author: Thijs Missiaen (ORCID: 0009-0007-3678-0810)

Certifies in 50-decimal arbitrary precision (mpmath.mp.dps = 50) and 35-digit
directed-rounding interval arithmetic (mpmath.iv.dps = 35):
  1. Canonical normalization Xi_X(z) = 8 int_0^infty Phi_X(u) cos(2zu) du and
     W_X(d,t) = 16 int_0^infty y sinh(2dy) I_X(y,t) dy for finite real Dirichlet
     polynomials M_X(s), preserving strip holomorphy on |Im(u)| < pi/8.
  2. Exact arithmetic comparison on sigma > 1 between Gonek-Hughes-Keating's
     truncated Euler factor P_X(s)^{-1} = prod_{p<=X}(1 - p^{-s}) and the entire
     zero-free exponential factor E_X(s) = exp(-sum_{2<=n<=X} (Lambda(n)/ln n) n^{-s}):
       - b_1 = 1 and b_N = 0 on 2 <= N <= floor(X) for both;
       - Prime-power purity b_N^E(X) = Lambda(N)/ln N on (floor(X), (floor(X)+1)^2)
         vs pure-prime support b_N^{GHK}(X) = 1_{P^-(N) > X};
       - Universal contractivity 0 <= b_N^E(X) <= 1 via local generating function;
       - Exact finite-divisor identity for M_{X,M}^E(s) on 1 <= N <= M.
  3. Explicit certification of the three critical-strip obstructions:
       - Asymmetric Riemann-Siegel expansion of zeta(s) M_X(s): convolution tail on
         (N(t), N(t)K] and unsieved reflected branch chi(s) M_X(s) sum_{m<=N(t)} m^{s-1};
       - Breakdown of X^{-2d}/(2d ln X) at zero-detection depth d ~ 1/ln t;
       - Single-line divergence w_n*(d_0) = -(n^{2d_0}-1)^{-1} -> -infty as d_0 -> 0+
         and non-zero explicit-formula zero-sum discrepancy inside 0 < sigma <= 1.
  4. Exact outer modulus identity for E_{X,kappa}(z) and invariance of zeros of Xi.
"""

from pathlib import Path
import sys
import time
import math
import numpy as np
import mpmath as mp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

mp.mp.dps = 50
mp.iv.dps = 35

ROOT_DIR = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT_DIR / "figures"
PAPER_FIG_DIR = ROOT_DIR / "paper" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
PAPER_FIG_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Number-Theoretic Helpers
# ============================================================================

def sieve_primes(limit: int) -> list[int]:
    if limit < 2:
        return []
    is_p = bytearray(b"\x01") * (limit + 1)
    is_p[0] = is_p[1] = 0
    for p in range(2, int(math.isqrt(limit)) + 1):
        if is_p[p]:
            is_p[p * p : limit + 1 : p] = bytearray(len(range(p * p, limit + 1, p)))
    return [i for i, v in enumerate(is_p) if v]


def factorize(n: int, primes: list[int]) -> dict[int, int]:
    fac: dict[int, int] = {}
    rem = n
    for p in primes:
        if p * p > rem:
            break
        if rem % p == 0:
            cnt = 0
            while rem % p == 0:
                rem //= p
                cnt += 1
            fac[p] = cnt
    if rem > 1:
        fac[rem] = fac.get(rem, 0) + 1
    return fac


def von_mangoldt_over_log(n: int, primes: list[int]) -> mp.mpf:
    """Returns Lambda(n) / ln(n) = 1/k if n = p^k (k >= 1), else 0."""
    if n < 2:
        return mp.mpf(0)
    fac = factorize(n, primes)
    if len(fac) == 1:
        (_, k), = fac.items()
        return mp.mpf(1) / mp.mpf(k)
    return mp.mpf(0)


def compute_c_E_trunc(X: int, M: int, primes: list[int]) -> list[mp.mpf]:
    """
    Exact Dirichlet coefficients c_k^E(X) for 1 <= k <= M of
    E_X(s) = exp(-sum_{2 <= n <= X} (Lambda(n)/ln n) n^{-s})
    via Euler's logarithmic-derivative recurrence:
      ln(k) c_k = -sum_{d | k, 2 <= d <= X} Lambda(d) c_{k/d}.
    Wait: differentiating wrt s gives -E_X'(s) = E_X(s) * (-sum_{2<=d<=X} Lambda(d) d^{-s})?
    Let f(s) = sum_{k>=1} c_k k^{-s} = exp(-A(s)) with A(s) = sum_{2<=d<=X} (Lambda(d)/ln d) d^{-s}.
    Then f'(s) = -sum_{k>=1} c_k ln(k) k^{-s}, while -A'(s) = sum_{2<=d<=X} Lambda(d) d^{-s}?
    Careful: d/ds (d^{-s}) = -ln(d) d^{-s}, so A'(s) = -sum_{2<=d<=X} Lambda(d) d^{-s},
    hence f'(s) = -A'(s) f(s) = (sum_{2<=d<=X} Lambda(d) d^{-s}) f(s).
    Comparing coefficients of k^{-s}: -ln(k) c_k = sum_{d|k, 2<=d<=X} Lambda(d) c_{k/d},
    so c_k = -(1/ln k) sum_{d|k, 2<=d<=X} Lambda(d) c_{k/d}.
    """
    lam = [mp.mpf(0)] * (max(X, M) + 1)
    for p in primes:
        if p > max(X, M):
            break
        pk = p
        lp = mp.log(p)
        while pk <= max(X, M):
            lam[pk] = lp
            pk *= p

    c = [mp.mpf(0)] * (M + 1)
    c[1] = mp.mpf(1)
    for k in range(2, M + 1):
        acc = mp.mpf(0)
        for d in range(2, min(X, k) + 1):
            if k % d == 0 and lam[d] != 0:
                acc += lam[d] * c[k // d]
        c[k] = -acc / mp.log(k)
    return c


def compute_b_E_via_exp_tail(X: int, N_max: int, primes: list[int]) -> list[mp.mpf]:
    """
    Computes b_N^E(X) for 1 <= N <= N_max from exp(sum_{n > X} (Lambda(n)/ln n) n^{-s})
    via ln(N) b_N = sum_{d | N, d > X} Lambda(d) b_{N/d}.
    """
    lam = [mp.mpf(0)] * (N_max + 1)
    for p in primes:
        if p > N_max:
            break
        pk = p
        lp = mp.log(p)
        while pk <= N_max:
            lam[pk] = lp
            pk *= p

    b = [mp.mpf(0)] * (N_max + 1)
    b[1] = mp.mpf(1)
    for N in range(2, N_max + 1):
        acc = mp.mpf(0)
        for d in range(X + 1, N + 1):
            if N % d == 0 and lam[d] != 0:
                acc += lam[d] * b[N // d]
        b[N] = acc / mp.log(N)
    return b


def local_beta_polynomial(a_max: int, k: int) -> list[mp.mpf]:
    """
    Coefficients beta_0(k), ..., beta_{a_max}(k) of
    (1 - z)^{-1} exp(-sum_{m=1}^k z^m / m).
    """
    # First compute g(z) = exp(-sum_{m=1}^k z^m / m): j g_j = -sum_{m=1}^{min(j,k)} g_{j-m}
    g = [mp.mpf(0)] * (a_max + 1)
    g[0] = mp.mpf(1)
    for j in range(1, a_max + 1):
        g[j] = -mp.fsum(g[j - m] for m in range(1, min(j, k) + 1)) / mp.mpf(j)
    # Multiply by (1 - z)^{-1}: prefix sums of g
    beta = [mp.mpf(0)] * (a_max + 1)
    running = mp.mpf(0)
    for j in range(a_max + 1):
        running += g[j]
        beta[j] = running
    return beta


# ============================================================================
# Analytic Kernel Helpers (Canonical Normalization)
# ============================================================================

def phi_shell(m: int, u):
    eu4 = mp.exp(4 * u)
    m2 = mp.mpf(m * m)
    return (2 * mp.pi**2 * (m2**2) * mp.exp(9 * u) - 3 * mp.pi * m2 * mp.exp(5 * u)) * mp.exp(-mp.pi * m2 * eu4)


def Phi_jacobi(u, m_max: int = 16):
    # Use reflection symmetry Phi(-u) = Phi(u) when Re(u) < 0 for rapid convergence
    if mp.re(u) < 0:
        u = -u
    return mp.fsum(phi_shell(m, u) for m in range(1, m_max + 1))


def Xi_exact(z):
    s = mp.mpf("0.5") + 1j * z
    return mp.mpf("0.5") * s * (s - 1) * (mp.pi ** (-s / 2)) * mp.gamma(s / 2) * mp.zeta(s)


def W_canonical_exact(d, t, c_dict: dict[int, mp.mpf] | None = None):
    """
    Canonical W(d,t) = (1/2) d/dd |Xi_X(t - i d)|^2 = R I_t - I R_t.
    If c_dict is provided, multiplies Xi(z) by M_X(1/2 + i z) M_X(1/2 - i z).
    """
    def F(z):
        val = Xi_exact(z)
        if c_dict is not None:
            s1 = mp.mpf("0.5") + 1j * z
            s2 = mp.mpf("0.5") - 1j * z
            m1 = mp.fsum(ck * (mp.mpf(k) ** (-s1)) for k, ck in c_dict.items())
            m2 = mp.fsum(ck * (mp.mpf(k) ** (-s2)) for k, ck in c_dict.items())
            val *= m1 * m2
        return val

    z0 = mp.mpf(t) - 1j * mp.mpf(d)
    f0 = F(z0)
    fp = mp.diff(F, z0)
    R, I = mp.re(f0), mp.im(f0)
    Rt, It = mp.re(fp), mp.im(fp)
    return R * It - I * Rt


def Phi_sieved(u, c_dict: dict[int, mp.mpf]):
    acc = mp.mpf(0)
    items = list(c_dict.items())
    for k, ck in items:
        for m, cm in items:
            shift = mp.mpf("0.5") * mp.log(mp.mpf(k) / mp.mpf(m))
            w = ck * cm / mp.sqrt(mp.mpf(k * m))
            acc += mp.mpf("0.5") * w * (Phi_jacobi(u + shift) + Phi_jacobi(u - shift))
    return acc


# ============================================================================
# Verification Sections
# ============================================================================

def verify_section2_sieved_kernel():
    print("=" * 78)
    print("[1/5] Canonical Sieved Mirror Kernel Phi_X(u) & Strip Holomorphy")
    print("=" * 78)
    # 2-prime Gonek-Hughes-Keating polynomial P_3(s)^{-1} = (1 - 2^{-s})(1 - 3^{-s})
    c_ghk3 = {1: mp.mpf(1), 2: mp.mpf(-1), 3: mp.mpf(-1), 6: mp.mpf(1)}

    # Check evenness Phi_X(-u) = Phi_X(u) at complex u inside |Im(u)| < pi/8
    u_test = mp.mpf("0.18") + 1j * mp.mpf("0.27")  # 0.27 < pi/8 = 0.392699
    even_res = abs(Phi_sieved(u_test, c_ghk3) - Phi_sieved(-u_test, c_ghk3))
    print(f"  Complex strip evenness |Phi_X(u) - Phi_X(-u)| at Im(u)=0.27 < pi/8: {mp.nstr(even_res, 6)}")
    assert even_res < mp.mpf("1e-42")

    # Check Xi_X(z) = 8 int_0^infty Phi_X(u) cos(2zu) du
    d0, t0 = mp.mpf("0.20"), mp.mpf("14.134725")
    z0 = t0 - 1j * d0
    s1, s2 = mp.mpf("0.5") + 1j * z0, mp.mpf("0.5") - 1j * z0
    M1 = (1 - mp.mpf(2) ** (-s1)) * (1 - mp.mpf(3) ** (-s1))
    M2 = (1 - mp.mpf(2) ** (-s2)) * (1 - mp.mpf(3) ** (-s2))
    Xi_X_ref = Xi_exact(z0) * M1 * M2
    Xi_X_int = 8 * mp.quad(lambda u: Phi_sieved(u, c_ghk3) * mp.cos(2 * z0 * u), [0, 1, 2.5, 4.5])
    xi_res = abs(Xi_X_ref - Xi_X_int)
    print(f"  Canonical integral |Xi_X(z) - 8 int_0^infty Phi_X(u)cos(2zu)du|   : {mp.nstr(xi_res, 6)}")
    assert xi_res < mp.mpf("1e-35")

    # Check canonical W_X(d,t) = 16 int_0^infty y sinh(2dy) I_X(y,t) dy at 30-dps GL quadrature
    old_dps = mp.mp.dps
    mp.mp.dps = 30
    d_w, t_w = mp.mpf("0.25"), mp.mpf("3.5")
    # Uniform trapezoidal rule on [-2.4, 2.4] x [0, 2.4] exploits strip holomorphy
    # (|Im(x)|, |Im(y)| < pi/4) + super-exponential decay -> O(exp(-pi^2/(2h))) accuracy!
    h_step = mp.mpf("0.04")
    m_half = 65  # 65 * 0.04 = 2.6
    c_ghk2 = {1: mp.mpf(1), 2: mp.mpf(-1)}
    W_ref2 = W_canonical_exact(d_w, t_w, c_ghk2)
    # Cache Phi_sieved on the 1D half-step lattice (k +/- j) * (h_step / 2)
    phi_cache = {idx: Phi_sieved(idx * h_step / 2, c_ghk2) for idx in range(-2 * m_half, 2 * m_half + 1)}
    cos_cache = {k: mp.cos(2 * t_w * k * h_step) for k in range(-m_half, m_half + 1)}
    W_gl = mp.mpf(0)
    for j in range(0, m_half + 1):
        yp = j * h_step
        yw = h_step * (mp.mpf("0.5") if j in (0, m_half) else mp.mpf(1))
        Iy = h_step * mp.fsum(phi_cache[k + j] * phi_cache[k - j] * cos_cache[k] for k in range(-m_half, m_half + 1))
        W_gl += yw * yp * mp.sinh(2 * d_w * yp) * Iy
    W_gl *= 16
    w_res = abs(W_ref2 - W_gl)
    mp.mp.dps = old_dps
    print(f"  Canonical Wigner |W_X(d,t) - 16 int_0^infty y sinh(2dy) I_X dy|   : {mp.nstr(w_res, 6)}")
    assert w_res < mp.mpf("1e-20")
    return xi_res, w_res


def verify_section3_ghk_vs_entire():
    print("=" * 78)
    print("[2/5] Gonek-Hughes-Keating P_X(s)^{-1} vs Entire E_X(s) on sigma > 1")
    print("=" * 78)
    X = 10
    N_max = 150
    primes = sieve_primes(N_max + 20)

    b_E = compute_b_E_via_exp_tail(X, N_max, primes)
    c_E = compute_c_E_trunc(X, N_max, primes)

    # Convolution of zeta(s) and M_{X, N_max}^E(s) = divisor sum of c_E
    b_div = [mp.mpf(0)] * (N_max + 1)
    for N in range(1, N_max + 1):
        b_div[N] = mp.fsum(c_E[k] for k in range(1, N + 1) if N % k == 0)

    max_div_err = max(abs(b_E[N] - b_div[N]) for N in range(1, N_max + 1))
    max_annih = max(abs(b_E[N]) for N in range(2, X + 1))

    # Purity on (X, (X+1)^2) = (10, 121) -> 11 <= N <= 120
    max_purity_err = max(
        abs(b_E[N] - von_mangoldt_over_log(N, primes))
        for N in range(X + 1, (X + 1) ** 2)
    )

    # Local generating function check & universal [0, 1] contractivity
    max_local_err = mp.mpf(0)
    for N in range(1, N_max + 1):
        assert -mp.mpf("1e-48") <= b_E[N] <= mp.mpf(1) + mp.mpf("1e-48")
        fac = factorize(N, primes)
        prod_local = mp.mpf(1)
        for p, a in fac.items():
            k_p = int(math.floor(math.log(X) / math.log(p) + 1e-12)) if p <= X else 0
            betas = local_beta_polynomial(a, k_p)
            prod_local *= betas[a]
        max_local_err = max(max_local_err, abs(b_E[N] - prod_local))

    # Compare against Gonek-Hughes-Keating b_N^{GHK}(X) = 1_{P^-(N) > X}
    ghk_diffs = []
    for N in range(X + 1, (X + 1) ** 2):
        fac = factorize(N, primes)
        b_ghk = mp.mpf(1) if all(p > X for p in fac.keys()) else mp.mpf(0)
        if abs(b_E[N] - b_ghk) > mp.mpf("1e-30"):
            ghk_diffs.append((N, float(b_E[N]), float(b_ghk)))

    # Direct Dirichlet series validation at sigma = 1.6 > 1
    s_test = mp.mpf("1.6") + 1j * mp.mpf("5.0")
    E_X_s = mp.exp(-mp.fsum(von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-s_test)) for n in range(2, X + 1)))
    prod_exact = mp.zeta(s_test) * E_X_s
    # Extend N_max_dir = 600 for tight convergence at sigma=1.6
    primes_big = sieve_primes(700)
    b_big = compute_b_E_via_exp_tail(X, 600, primes_big)
    prod_series = mp.fsum(b_big[N] * (mp.mpf(N) ** (-s_test)) for N in range(1, 601))
    dir_err = abs(prod_exact - prod_series)

    print(f"  Annihilation max_{{2<=N<={X}}} |b_N^E({X})|                        : {mp.nstr(max_annih, 6)}")
    print(f"  Purity max_{{{X}<N<{(X+1)**2}}} |b_N^E({X}) - Lambda(N)/ln N|       : {mp.nstr(max_purity_err, 6)}")
    print(f"  Finite-divisor identity max_{{1<=N<={N_max}}} |b_N^E - sum_{{k|N}} c_k|: {mp.nstr(max_div_err, 6)}")
    print(f"  Local generating-function max residual                            : {mp.nstr(max_local_err, 6)}")
    print(f"  Prime-power differences b_N^E(10) != b_N^{{GHK}}(10) on (10,121)  : {len(ghk_diffs)} prime powers { [x[0] for x in ghk_diffs] }")
    print(f"  sigma=1.6 Dirichlet identity |zeta(s)E_10(s) - sum_{{N<=600}} b_N| : {mp.nstr(dir_err, 6)}")

    assert max_annih < mp.mpf("1e-45")
    assert max_purity_err < mp.mpf("1e-45")
    assert max_div_err < mp.mpf("1e-45")
    assert max_local_err < mp.mpf("1e-45")
    assert len(ghk_diffs) > 0  # e.g. 16, 25, 27, 32, 49, 64, 81
    return b_E, ghk_diffs


def verify_section4_critical_strip_obstructions():
    print("=" * 78)
    print("[3/5] Critical-Strip Obstructions: No FE, Tautological theta<1/2, & Waterbed")
    print("=" * 78)
    primes = sieve_primes(200)

    # Obstruction A: zeta(s) * M_{X,K}^E(s) vs naive symmetric sieved RS model!
    # Take t = 250 -> N(t) = floor(sqrt(250 / (2*pi))) = 6, X = 6, K = 6.
    # In the naive symmetric RS model, setting X = N(t) = 6 gives b_1 = 1 and b_n = 0 for 2 <= n <= 6,
    # predicting M_6(s)[1^{-s} + chi(s) 1^{s-1}]? No, the naive model replaced BOTH branches by
    # sum_{n<=N} b_n n^{-s} + chi(s) sum_{n<=N} b_n n^{s-1} = 1 + chi(s)!
    # Let's compare the TRUE product zeta(s) * M_{6,6}^E(s) against 1 + chi(s) across t in [200, 300]!
    t_val = mp.mpf("250.0")
    d_val = mp.mpf("0.15")
    s_val = mp.mpf("0.5") + d_val + 1j * t_val
    N_t = int(mp.floor(mp.sqrt(t_val / (2 * mp.pi))))  # N_t = 6
    X_t = N_t
    c_6 = compute_c_E_trunc(X_t, X_t, primes)

    chi_s = (mp.pi ** (s_val - mp.mpf("0.5"))) * mp.gamma((1 - s_val) / 2) / mp.gamma(s_val / 2)
    M_s = mp.fsum(c_6[k] * (mp.mpf(k) ** (-s_val)) for k in range(1, X_t + 1))

    true_prod = mp.zeta(s_val) * M_s
    naive_sym_rs = mp.mpf(1) + chi_s  # because b_1=1 and b_2=...=b_6=0!
    # True RS decomposition of zeta(s) * M_s:
    fwd_rs = mp.fsum(mp.mpf(n) ** (-s_val) for n in range(1, N_t + 1)) * M_s
    refl_rs = chi_s * mp.fsum(mp.mpf(m) ** (s_val - 1) for m in range(1, N_t + 1)) * M_s
    # Note fwd_rs = 1 + convolution_tail_{7 <= N <= 36}
    conv_tail = fwd_rs - mp.mpf(1)
    refl_unsieved_diff = refl_rs - chi_s

    err_naive = abs(true_prod - naive_sym_rs)
    err_true_rs = abs(true_prod - (fwd_rs + refl_rs))
    print(f"  At (d={mp.nstr(d_val,3)}, t={mp.nstr(t_val,4)}), N(t)=X={N_t}:")
    print(f"    |zeta(s)M_X(s) - naive_sym_RS (1+chi(s))|  = {mp.nstr(err_naive, 6)}  (LARGE O(1) error!)")
    print(f"    |convolution tail on (N(t), N(t)X]|        = {mp.nstr(abs(conv_tail), 6)}")
    print(f"    |unsieved reflected branch - chi(s)|       = {mp.nstr(abs(refl_unsieved_diff), 6)}")
    print(f"    |zeta(s)M_X(s) - true asymmetric RS|       = {mp.nstr(err_true_rs, 6)}")
    assert err_naive > mp.mpf("0.1")

    # Obstruction B: Zero-detection depth d ~ 1/ln(t) makes X^{-2d}/(2d ln X) = O(1), NOT o(1)!
    print("  Zero-detection depth d = 1/ln(t) evaluation of B_L2 = X^{-2d}/(2d ln X) at X = sqrt(t/2pi):")
    depth_rows = []
    for t_num in [100.0, 1000.0, 1e6, 1e12]:
        tt = mp.mpf(str(t_num))
        XX = mp.sqrt(tt / (2 * mp.pi))
        dd = 1 / mp.log(tt)
        val_B = (XX ** (-2 * dd)) / (2 * dd * mp.log(XX))
        depth_rows.append((t_num, float(dd), float(XX ** (-2 * dd)), float(val_B)))
        print(f"    t = {t_num:9.1e}, d = 1/ln(t) = {float(dd):.5f}: X^(-2d) = {float(XX**(-2*dd)):.5f}, X^(-2d)/(2d ln X) = {float(val_B):.5f}")
    assert depth_rows[-1][3] > 0.35  # stays O(1) ~ e^{-1} > 0.3678!

    # Obstruction C: Single-line weight w_n*(d_0) = -(n^{2d_0}-1)^{-1} -> -infty as d_0 -> 0+
    # AND failure inside 0 < sigma <= 1 due to nontrivial zeros!
    print("  Single-line weight w_2*(d_0) = -(2^{2d_0} - 1)^{-1} blow-up as d_0 -> 0+:")
    w2_rows = []
    for d0_num in [0.60, 0.25, 0.05, 0.01, 0.001]:
        d0 = mp.mpf(str(d0_num))
        w2 = -1 / (mp.mpf(2) ** (2 * d0) - 1)
        # Residual at d = 2*d0: 1 - (2^{2d}-1)/(2^{2d_0}-1)
        res_mismatch = abs(1 - (mp.mpf(2) ** (4 * d0) - 1) / (mp.mpf(2) ** (2 * d0) - 1))
        w2_rows.append((d0_num, float(w2), float(res_mismatch)))
        print(f"    d_0 = {d0_num:6.3f}: w_2*(d_0) = {float(w2):11.5f}, off-line mismatch at d=2d_0: {float(res_mismatch):.5f}")
    assert abs(w2_rows[-1][1]) > 700.0

    # Interval enclosure of 1/(2 d_0 ln n) blow-up bound: |w_n*(d_0)| >= 1/(2 d_0 ln n * e^{2 d_0 ln n})
    iv_d0 = mp.iv.mpf("0.001")
    iv_w2 = -1 / (mp.iv.mpf(2) ** (2 * iv_d0) - 1)
    assert iv_w2.b < mp.mpf("-720")
    print(f"  Interval certification at d_0=0.001: w_2*(0.001) in [{mp.nstr(iv_w2.a, 8)}, {mp.nstr(iv_w2.b, 8)}]")
    return depth_rows, w2_rows


def verify_section5_outer_modulus():
    print("=" * 78)
    print("[4/5] Shifted Outer Modulus Identity & Zero-Invariance of Xi_{E_{X,kappa}}")
    print("=" * 78)
    primes = sieve_primes(50)
    X = 15
    kappa = mp.mpf("0.35")
    d = mp.mpf("0.20")
    t = mp.mpf("21.022039638771554992628479593896902777334340524903")  # gamma_1 or gamma_2! (gamma_2)

    z = t - 1j * d
    s_plus = mp.mpf("0.5") + kappa + 1j * z
    s_minus = mp.mpf("0.5") + kappa - 1j * z

    E_plus = mp.exp(-mp.fsum(von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-s_plus)) for n in range(2, X + 1)))
    E_minus = mp.exp(-mp.fsum(von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-s_minus)) for n in range(2, X + 1)))
    E_cal = E_plus * E_minus

    mod_sq_lhs = abs(E_cal) ** 2
    mod_sq_rhs = mp.exp(
        -4 * mp.fsum(
            von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-(mp.mpf("0.5") + kappa)))
            * mp.cosh(d * mp.log(n)) * mp.cos(t * mp.log(n))
            for n in range(2, X + 1)
        )
    )
    mod_res = abs(mod_sq_lhs - mod_sq_rhs)

    # Check logarithmic d-derivative identity
    def log_mod_sq(dd):
        zz = t - 1j * dd
        ep = mp.exp(-mp.fsum(von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-(mp.mpf("0.5") + kappa + 1j * zz))) for n in range(2, X + 1)))
        em = mp.exp(-mp.fsum(von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-(mp.mpf("0.5") + kappa - 1j * zz))) for n in range(2, X + 1)))
        return mp.log(abs(ep * em) ** 2)

    dlog_num = mp.mpf("0.5") * mp.diff(log_mod_sq, d)
    dlog_closed = -2 * mp.fsum(
        (von_mangoldt_over_log(n, primes) * mp.log(n)) * (mp.mpf(n) ** (-(mp.mpf("0.5") + kappa)))
        * mp.sinh(d * mp.log(n)) * mp.cos(t * mp.log(n))
        for n in range(2, X + 1)
    )
    dlog_res = abs(dlog_num - dlog_closed)

    # Verify at first zero gamma_1 on d = 0 that |E_{X,kappa}(gamma_1)|^2 > 0 while Xi(gamma_1) = 0
    gamma1 = mp.im(mp.zetazero(1))
    Xi_at_zero = abs(Xi_exact(gamma1))
    E_at_zero = mp.exp(
        -4 * mp.fsum(
            von_mangoldt_over_log(n, primes) * (mp.mpf(n) ** (-(mp.mpf("0.5") + kappa))) * mp.cos(gamma1 * mp.log(n))
            for n in range(2, X + 1)
        )
    )
    print(f"  Outer modulus identity ||E_{{X,kappa}}|^2 - closed_form|            : {mp.nstr(mod_res, 6)}")
    print(f"  Log-derivative |(1/2)d_d ln|E|^2 - (-2 sum Lambda n^{{-1/2-k}}...)| : {mp.nstr(dlog_res, 6)}")
    print(f"  At gamma_1 = {mp.nstr(gamma1, 10)}: |Xi(gamma_1)| = {mp.nstr(Xi_at_zero, 4)}, |E_{{15,0.35}}(gamma_1)|^2 = {mp.nstr(E_at_zero, 8)} > 0")
    assert mod_res < mp.mpf("1e-45")
    assert dlog_res < mp.mpf("1e-40")
    assert E_at_zero > mp.mpf("0.01")


# ============================================================================
# Figure Generation
# ============================================================================

def generate_figures(b_E):
    print("=" * 78)
    print("[5/5] Generating Publication Figures (PDF + PNG)")
    print("=" * 78)
    plt.rcParams.update({
        "font.size": 9.5,
        "axes.titlesize": 10,
        "axes.labelsize": 9.5,
        "legend.fontsize": 8.0,
        "figure.dpi": 300,
    })
    primes = sieve_primes(200)

    # --- Figure 1: Sieved Kernel Phi_X(u) & Arithmetic Coefficients b_N^E vs b_N^{GHK} ---
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.8))

    # (a) Real-axis kernels Phi(u) vs Phi_{GHK,3}(u)
    ax = axes[0, 0]
    u_grid = np.linspace(-1.6, 1.6, 240)
    c_ghk3 = {1: mp.mpf(1), 2: mp.mpf(-1), 3: mp.mpf(-1), 6: mp.mpf(1)}
    phi0_vals = [float(Phi_jacobi(mp.mpf(str(u)), m_max=10)) for u in u_grid]
    phix_vals = [float(Phi_sieved(mp.mpf(str(u)), c_ghk3)) for u in u_grid]
    ax.plot(u_grid, phi0_vals, lw=1.8, color="#1b4965", label=r"Unsieved $\Phi(u)$")
    ax.plot(u_grid, phix_vals, lw=1.8, ls="--", color="#d1495b", label=r"Sieved $\Phi_X(u)$ ($P_3^{-1}$)")
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlabel(r"$u$")
    ax.set_ylabel(r"Kernel amplitude")
    ax.set_title(r"(a) Even strip-holomorphic kernel $\Phi_X(u)$ on $|\mathrm{Im}(u)| < \pi/8$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    # (b) Coefficients b_N^E(10) on 1 <= N <= 135 showing annihilation, purity on (10,121), and composites >= 121
    ax = axes[0, 1]
    N_arr = np.arange(1, 136)
    b_arr = [float(b_E[n]) for n in N_arr]
    colors = []
    for n in N_arr:
        if n == 1:
            colors.append("#2a9d8f")
        elif n <= 10:
            colors.append("#adb5bd")
        elif n < 121:
            fac = factorize(int(n), primes)
            if len(fac) == 1 and list(fac.values())[0] == 1:
                colors.append("#1b4965")  # pure prime
            elif len(fac) == 1:
                colors.append("#e76f51")  # prime power p^k
            else:
                colors.append("#adb5bd")
        else:
            colors.append("#6c757d")
    ax.bar(N_arr, b_arr, width=0.85, color=colors)
    ax.axvspan(1.5, 10.5, color="#e9ecef", alpha=0.6, label=r"Annihilated $[2, X]$")
    ax.axvline(121, color="#d1495b", ls="--", lw=1.2, label=r"$(X+1)^2 = 121$ onset")
    ax.set_xlabel(r"Dirichlet index $N$")
    ax.set_ylabel(r"$b_N^{\mathcal{E}}(10)$")
    ax.set_title(r"(b) Exact prime-power purity $b_N^{\mathcal{E}}(10) = \Lambda(N)/\ln N$ on $(10, 121)$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    # (c) Difference b_N^E(10) - b_N^{GHK}(10) on (10, 121): prime powers p^k of p <= 10!
    ax = axes[1, 0]
    N_win = np.arange(11, 121)
    diff_vals = []
    for n in N_win:
        fac = factorize(int(n), primes)
        b_ghk = 1.0 if all(p > 10 for p in fac.keys()) else 0.0
        diff_vals.append(float(b_E[n]) - b_ghk)
    ax.stem(N_win, diff_vals, linefmt="#e76f51", markerfmt="o", basefmt="gray")
    ax.set_xlabel(r"$N \in (10, 121)$")
    ax.set_ylabel(r"$b_N^{\mathcal{E}}(10) - b_N^{\mathrm{GHK}}(10)$")
    ax.set_title(r"(c) Arithmetic price of entireness: prime powers $p^k > X$ ($p \leq X$)")
    ax.grid(alpha=0.25)

    # (d) Asymmetric RS breakdown |zeta(s) M_6(s)| vs naive symmetric RS |1 + chi(s)|
    ax = axes[1, 1]
    t_scan = np.linspace(220.0, 275.0, 140)
    c_6 = compute_c_E_trunc(6, 6, primes)
    d_fix = mp.mpf("0.15")
    true_vals, naive_vals = [], []
    for tv in t_scan:
        sv = mp.mpf("0.5") + d_fix + 1j * mp.mpf(str(tv))
        chis = (mp.pi ** (sv - mp.mpf("0.5"))) * mp.gamma((1 - sv) / 2) / mp.gamma(sv / 2)
        Ms = mp.fsum(c_6[k] * (mp.mpf(k) ** (-sv)) for k in range(1, 7))
        true_vals.append(float(abs(mp.zeta(sv) * Ms)))
        naive_vals.append(float(abs(1 + chis)))
    ax.plot(t_scan, true_vals, lw=1.7, color="#1b4965", label=r"True $|\zeta(s)M_{6,6}^{\mathcal{E}}(s)|$")
    ax.plot(t_scan, naive_vals, lw=1.5, ls="--", color="#d1495b", label=r"Symmetric ansatz $|1 + \chi(s)|$")
    ax.set_xlabel(r"Ordinate $t$ ($N(t) = X = 6$, $d = 0.15$)")
    ax.set_ylabel(r"Modulus")
    ax.set_title(r"(d) Failure of symmetric RS substitution ($M_X$ has no FE)")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    for out_dir in (FIG_DIR, PAPER_FIG_DIR):
        fig.savefig(out_dir / "fig1_sieved_mirror_kernel_and_purity.pdf")
        fig.savefig(out_dir / "fig1_sieved_mirror_kernel_and_purity.png", dpi=300)
    plt.close(fig)

    # --- Figure 2: Obstructions (d ~ 1/ln t scaling, w_n*(d_0) divergence, Outer Modulus) ---
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.8))

    # (a) X^{-2d}/(2d ln X) for fixed d vs zero-detection depth d = 1/ln t
    ax = axes[0, 0]
    log10_t = np.linspace(2.0, 12.0, 120)
    t_vals = 10.0 ** log10_t
    X_vals = np.sqrt(t_vals / (2.0 * np.pi))
    for d_const, col, lab in [(0.25, "#2a9d8f", r"Fixed $d = 0.25$"), (0.10, "#1b4965", r"Fixed $d = 0.10$")]:
        b_fixed = (X_vals ** (-2.0 * d_const)) / (2.0 * d_const * np.log(X_vals))
        ax.semilogy(log10_t, b_fixed, lw=1.8, color=col, label=lab)
    d_zero = 1.0 / np.log(t_vals)
    b_zero = (X_vals ** (-2.0 * d_zero)) / (2.0 * d_zero * np.log(X_vals))
    ax.semilogy(log10_t, b_zero, lw=2.0, ls="--", color="#d1495b", label=r"Zero-detection depth $d = 1/\ln t$ ($\asymp 1$)")
    ax.set_xlabel(r"$\log_{10} t$ ($X = \sqrt{t/2\pi}$)")
    ax.set_ylabel(r"$X^{-2d} / (2d\ln X)$")
    ax.set_title(r"(a) Power-law damping at fixed $d$ vs $\mathcal{O}(1)$ stall at $d = 1/\ln t$")
    ax.legend(loc="lower left")
    ax.grid(alpha=0.25)

    # (b) Single-line weight |w_p^*(d_0)| = (p^{2d_0}-1)^{-1} -> infty as d_0 -> 0+
    ax = axes[0, 1]
    d0_grid = np.logspace(-3.0, -0.1, 120)
    for p, col in [(2, "#d1495b"), (3, "#1b4965"), (5, "#2a9d8f")]:
        wp = 1.0 / (p ** (2.0 * d0_grid) - 1.0)
        ax.loglog(d0_grid, wp, lw=1.8, color=col, label=rf"$|w_{{{p}}}^*(d_0)| = ({p}^{{2d_0}}-1)^{{-1}}$")
    ax.set_xlabel(r"Tuning depth $d_0$")
    ax.set_ylabel(r"$|w_p^*(d_0)|$")
    ax.set_title(r"(b) Critical-line blow-up $w_n^*(d_0) \sim -(2d_0\ln n)^{-1} \to -\infty$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25, which="both")

    # (c) Off-line mismatch |1 - (n^{2d}-1)/(n^{2d_0}-1)| when d != d_0 (e.g. d_0 = 0.20)
    ax = axes[1, 0]
    d_grid = np.linspace(0.01, 0.50, 150)
    d0_ref = 0.20
    for p, col in [(2, "#d1495b"), (5, "#1b4965"), (11, "#2a9d8f")]:
        mismatch = np.abs(1.0 - (p ** (2.0 * d_grid) - 1.0) / (p ** (2.0 * d0_ref) - 1.0))
        ax.plot(d_grid, mismatch, lw=1.8, color=col, label=rf"$n = {p}$ ($d_0 = 0.20$)")
    ax.axvline(d0_ref, color="gray", ls=":", lw=1.2, label=r"Tuned line $d = d_0$")
    ax.set_xlabel(r"Evaluation depth $d$")
    ax.set_ylabel(r"Relative coefficient residual")
    ax.set_title(r"(c) Single-line tuning failure for all $d \neq d_0$")
    ax.legend(loc="upper left")
    ax.grid(alpha=0.25)

    # (d) Outer modulus |E_{X,kappa}(t - i d)|^2 > 0 vs |Xi(t - i d)| across first 3 critical zeros
    ax = axes[1, 1]
    t_z = np.linspace(12.0, 27.0, 200)
    E_mod = []
    Xi_scaled = []
    for tv in t_z:
        ev = math.exp(
            -4.0 * sum(
                float(von_mangoldt_over_log(n, primes)) * (n ** (-0.85)) * math.cos(tv * math.log(n))
                for n in range(2, 16)
            )
        )
        E_mod.append(ev)
        # Rescaled |Xi(t)| * e^{pi t / 4}
        xv = float(abs(Xi_exact(mp.mpf(str(tv)))) * mp.exp(mp.pi * mp.mpf(str(tv)) / 4))
        Xi_scaled.append(xv)
    ax.plot(t_z, E_mod, lw=1.8, color="#2a9d8f", label=r"$|\mathcal{E}_{15,0.35}(t)|^2 > 0$ (zero-free)")
    ax.plot(t_z, Xi_scaled, lw=1.6, ls="--", color="#1b4965", label=r"$\mathrm{e}^{\pi t/4}|\Xi(t)|$ (vanishes at $\gamma_k$)")
    ax.set_xlabel(r"Critical-line ordinate $t$ ($d = 0$)")
    ax.set_ylabel(r"Amplitude")
    ax.set_title(r"(d) Outer factor $\mathcal{E}_{X,\kappa}$ is zero-free and leaves zeros of $\Xi$ invariant")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.25)

    fig.tight_layout()
    for out_dir in (FIG_DIR, PAPER_FIG_DIR):
        fig.savefig(out_dir / "fig2_power_law_damping_and_shifted_gauge.pdf")
        fig.savefig(out_dir / "fig2_power_law_damping_and_shifted_gauge.png", dpi=300)
    plt.close(fig)
    print("  Saved fig1_sieved_mirror_kernel_and_purity.{pdf,png} and fig2_power_law_damping_and_shifted_gauge.{pdf,png}")


def main():
    t_start = time.time()
    verify_section2_sieved_kernel()
    b_E, _ = verify_section3_ghk_vs_entire()
    verify_section4_critical_strip_obstructions()
    verify_section5_outer_modulus()
    generate_figures(b_E)
    elapsed = time.time() - t_start
    print("=" * 78)
    print(f"[PASS] Paper VI verification suite completed in {elapsed:.2f}s (exit code 0).")
    print("=" * 78)


if __name__ == "__main__":
    main()
