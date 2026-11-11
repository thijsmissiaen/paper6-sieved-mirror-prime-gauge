---
name: sieved-mirror-prime-gauge
description: Scientific knowledge skill for Paper VI of the Transverse Wronskian series on Riemann's Xi-function. Contains the strip-holomorphic sieved mirror kernel Phi_X(u) under canonical Wronskian normalization, the exact sigma > 1 arithmetic comparison between Gonek-Hughes-Keating's truncated Euler product P_X(s)^{-1} and the entire zero-free exponential factor E_X(s) (exact annihilation on [2, X], prime-power purity b_N^E(X) = Lambda(N)/ln N on (X, (floor(X)+1)^2), universal [0,1] contractivity, and finite-divisor truncation), and the three critical-strip obstructions in 0 < sigma <= 1 (no functional equation / unsieved reflected RS branch, tautology of theta < 1/2 and stall at d ~ 1/ln t, and single-line divergence w_n*(d_0) -> -infty with Selberg's explicit-formula waterbed).
---

# Sieved Mirror Kernels, Gonek–Hughes–Keating vs Entire Exponential Prime Factors, and Critical-Strip Obstructions (Paper VI)

## 1. Canonical Sieved Kernel and Strip Holomorphy

For any finite real Dirichlet polynomial $M _ X(s) = \sum _ {k \in \mathcal{I} _ X} c _ k k^{-s}$ and $\mathcal{M} _ X(z) := M _ X(\frac{1}{2}+\mathrm{i} z)M _ X(\frac{1}{2}-\mathrm{i} z)$, the product $\Xi _ X(z) := \Xi(z)\mathcal{M} _ X(z)$ satisfies

$$
\Xi _ X(z) = 8\int _ 0^{\infty} \Phi _ X(u)\cos(2zu)\mathrm{d}u, \qquad W _ X(d,t) := \frac{1}{2}\frac{\partial}{\partial d}\left\lvert\Xi _ X(t - \mathrm{i} d)\right\rvert^2 = 16\int _ 0^{\infty} y\sinh(2dy) I _ X(y,t)\mathrm{d}y,
$$

where the even sieved kernel

$$
\Phi _ X(u) := \frac{1}{2}\sum _ {k,m \in \mathcal{I} _ X} \frac{c _ k c _ m}{\sqrt{km}}\left[\Phi\left(u + \frac{1}{2}\ln\frac{k}{m}\right) + \Phi\left(u - \frac{1}{2}\ln\frac{k}{m}\right)\right]
$$

remains holomorphic throughout the full strip $\lvert\mathrm{Im}(u)\rvert < \frac{\pi}{8}$.

## 2. Rigorous $\sigma > 1$ Arithmetic Purity: Gonek–Hughes–Keating vs $E _ X(s)$

On $\sigma = \mathrm{Re}(s) > 1$:
- **Gonek–Hughes–Keating (2007):** $P _ X(s)^{-1} := \prod _ {p \le X}(1 - p^{-s})$ has $\zeta(s)P _ X(s)^{-1} = \sum _ {N=1}^{\infty} b _ N^{\mathrm{GHK}}(X)N^{-s}$ with $b _ N^{\mathrm{GHK}}(X) = \mathbf{1} _ {P^-(N) > X} \in \lbrace 0,1\rbrace$ (pure primes on $(X, p _ {>X}^2)$, at the cost of zeros of $P _ X(s)^{-1}$ on $\sigma = 0$).
- **Entire Exponential Factor:** $E _ X(s) := \exp(-\sum _ {2 \le n \le X} \frac{\Lambda(n)}{\ln n}n^{-s})$ is entire and zero-free on $\mathbb{C}$, and on $\sigma > 1$ its product $\zeta(s)E _ X(s) = \sum _ {N=1}^{\infty} b _ N^{\mathcal{E}}(X)N^{-s}$ satisfies:
  1. $b _ 1^{\mathcal{E}}(X) = 1$ and $b _ N^{\mathcal{E}}(X) = 0$ for all $2 \le N \le \lfloor X\rfloor$;
  2. **Prime-Power Purity on $(\lfloor X\rfloor, (\lfloor X\rfloor+1)^2)$:** $b _ N^{\mathcal{E}}(X) = \frac{\Lambda(N)}{\ln N}$ ($1/k$ at prime powers $p^k$, $0$ whenever $\omega(N) \ge 2$), differing from $b _ N^{\mathrm{GHK}}(X)$ exclusively at prime powers $p^k > X$ of primes $p \le X$;
  3. **Universal Contractivity:** $0 \le b _ N^{\mathcal{E}}(X) = \prod _ {p^a \parallel N} \beta _ a(\lfloor\ln X/\ln p\rfloor) \le 1$ via $(1-z)^{-1}\exp(-\sum _ {m=1}^k z^m/m) = \exp(\sum _ {m=k+1}^{\infty} z^m/m)$;
  4. **Exact Finite-Divisor Identity:** Truncating $E _ X(s)$ at $k \le M$ preserves every coefficient $b _ N^{\mathcal{E}}(X)$ for $1 \le N \le M$.

## 3. Critical-Strip Obstructions ($0 < \sigma \le 1$)

1. **Divergence of $\mathcal{P} _ {>X}(s)$ and Lack of Functional Equation:** Inside $0 < \sigma \le 1$, $\mathcal{P} _ {>X}(s) = \sum _ {n > X} \frac{\Lambda(n)}{\ln n}n^{-s}$ diverges and $\zeta(s)E _ X(s)$ has the exact same zeros as $\zeta(s)$. Neither $E _ X(s)$ nor any finite $M _ X(s)$ satisfies a functional equation: multiplying $\zeta(s)$ by $M _ X(s)$ creates an incomplete convolution tail on $(N(t), N(t)K]$ and leaves the reflected Riemann–Siegel branch $\chi(s)M _ X(s)\sum _ {m \le N(t)} m^{s-1}$ unsieved with mixed exponents $k^{-s}m^{s-1}$.
2. **Tautology of $\theta < \frac{1}{2}$ and Zero-Detection Stall:** Setting $X = N(t) = \lfloor\sqrt{t/2\pi}\rfloor$ in a formal $N(t)$-term sum makes $b _ n^{\mathcal{E}}(N(t)) = 0$ ($2 \le n \le N(t){}$) tautological while ignoring the unsieved reflected branch; moreover, at zero-detection depth $d = c/\ln t$, $\frac{X^{-2d}}{2d\ln X} \to \frac{\mathrm{e}^{-c}}{c} = \mathcal{O}(1)$ does not decay.
3. **Single-Line Weight Blow-Up and Selberg's Explicit-Formula Waterbed:** Matching coefficients on $\sigma _ 0 = \frac{1}{2} + d _ 0 > 1$ forces $w _ n^{\ast}(d _ 0) = -(n^{2d _ 0}-1)^{-1}$, which fails for $d \ne d _ 0$, diverges as $-(2d _ 0\ln n)^{-1} \to -\infty$ as $d _ 0 \to 0^+$, and inside $0 < \sigma \le 1$ cannot cancel the meromorphic zero poles $\sum _ \rho \frac{\widetilde{\eta}(\rho - s)X^{\rho - s}}{\rho - s}$ ($\widetilde{\eta}(w) := -\int _ 0^{\infty} \eta'(u)u^w \mathrm{d}u$, $\widetilde{\eta}(0) = 1$) in Selberg's explicit formula.

## 4. Shifted Outer Modulus Identity

For $\mathcal{E} _ {X,\kappa}(z) := E _ X(\frac{1}{2}+\kappa+\mathrm{i} z)E _ X(\frac{1}{2}+\kappa-\mathrm{i} z)$ ($\kappa \ge 0$),

$$
\left\lvert\mathcal{E} _ {X,\kappa}(t - \mathrm{i} d)\right\rvert^2 = \exp\left(-4\sum _ {2 \le n \le X} \frac{\Lambda(n)}{n^{1/2+\kappa}\ln n}\cosh(d\ln n)\cos(t\ln n)\right) > 0,
$$

leaving the zero set of $\Xi(z)$ invariant.

## 5. Reproducible Verification

```bash
python scripts/verify_paper6_sieved_mirror_prime_gauge.py
```
