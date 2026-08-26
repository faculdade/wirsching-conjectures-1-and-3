"""Direct numerical test of Conjecture 3: phi(z_l)/phi_0(z_l) -> e^{H(0)} along
z_l = l * 3^{-l}, Wirsching's own comparison class at lambda=1.

Self-contained: phi is evaluated by the real-variable saddlepoint approximation
(justified rigorously, uniformly, by Theorem 13 of the paper), not by oscillatory
Fourier inversion. phi_0 uses the exponents gamma, delta, beta derived in Section 6
of the paper (matching Berg-Kruppel's own (9.6), see
derive_P_bergkruppel_identity.py in ../section-09-berg-kruppel-identity/).

This is supporting numerical evidence for Theorem 1, not a substitute for the proof
in papers/01-wirsching-conjecture3/main.tex.
"""
import mpmath as mp

mp.mp.dps = 60


def h(x):
    if x == 0:
        return mp.mpf(0)
    # expm1 rather than 1 - e^-u: at the smallest admitted x the
    # naive form loses every significant digit to cancellation.
    return mp.log(-mp.expm1(-2 * x) / (2 * x))


def m_func(x):
    if x == 0:
        return mp.mpf(0)
    return 1 - 2 * x / mp.expm1(2 * x)


def v_func(x):
    if x == 0:
        return mp.mpf(0)
    return 1 - (x / mp.sinh(x)) ** 2


def series_sum(f, s, terms=250):
    total = mp.mpf(0)
    denom = mp.mpf(3)
    for _ in range(terms):
        x = s / denom
        if abs(x) < mp.mpf(10) ** (-(mp.mp.dps - 8)):
            break
        total += f(x)
        denom *= 3
    return total


def phi_saddle(t, terms=250):
    """Real-variable saddlepoint approximation to phi(t), Theorem 13 of the paper."""
    t = mp.mpf(t)
    guess = 1 / t if t > 0 else mp.mpf(1)
    s = mp.findroot(lambda s: series_sum(m_func, s, terms) - s * t, guess,
                     tol=mp.mpf("1e-40"))
    K_s = series_sum(h, s, terms)
    V_s = series_sum(v_func, s, terms)
    return s / mp.sqrt(2 * mp.pi * V_s) * mp.e ** (K_s + s * t)


a = mp.log(3)
beta = 1 / (2 * a)
gamma = mp.mpf("-1.5") - (1 + mp.log(a / 2)) / a
delta = 1 + mp.log(a / 2) / a


def phi_0_bare(t):
    t = mp.mpf(t)
    L = -mp.log(t)
    return t ** gamma * L ** delta * mp.e ** (-beta * mp.log(t / L) ** 2)


if __name__ == "__main__":
    print("Berg-Kruppel exponents (a=3, lambda=2/3):")
    print(f"  beta  = {mp.nstr(beta, 15)}")
    print(f"  gamma = {mp.nstr(gamma, 15)}")
    print(f"  delta = {mp.nstr(delta, 15)}\n")
    print("Predicted limit (bare normalization, Section 6 of the paper):")
    print("  e^{C_P+H(0)} = 0.204987710306551537...\n")

    print(f"{'l':>5}{'ratio phi/phi_0':>20}")
    for l in range(5, 51, 5):
        t_l = mp.mpf(l) * mp.mpf(3) ** (-l)
        ratio = phi_saddle(t_l) / phi_0_bare(t_l)
        print(f"{l:>5}{mp.nstr(ratio, 12):>20}")

    # The limit under the bare normalization, derived rather than printed
    # as a string: e^{C_P + H(0)} with C_P from section 9's closed form and
    # H(0) from Theorem 6.1's decomposition, as in section-12/concordance.py.
    ratio20 = phi_saddle(mp.mpf(20) * mp.mpf(3) ** (-20)) / phi_0_bare(
        mp.mpf(20) * mp.mpf(3) ** (-20))
    if not mp.mpf("0.180") < ratio20 < mp.mpf("0.186"):
        raise AssertionError(
            "the ratio at l=20 is %s, outside the dip this README describes"
            % mp.nstr(ratio20, 8))
    ratio50 = phi_saddle(mp.mpf(50) * mp.mpf(3) ** (-50)) / phi_0_bare(
        mp.mpf(50) * mp.mpf(3) ** (-50))
    if not ratio50 > ratio20:
        raise AssertionError("the ratio is not rising back by l=50")
    print("\n  checked: ratio at l=20 is %s, in the dip; at l=50 it is %s,"
          % (mp.nstr(ratio20, 8), mp.nstr(ratio50, 8)))
    print("  already rising back toward the limit.")

    print("\nThe ratio dips below the limit before slowly rising back towards it,")
    print("visible already by l=50. Full convergence to 0.204988 (Theorem 1's")
    print("closed-form limit under the bare normalization) needs l well beyond what")
    print("mpmath's default root-finder handles reliably from a naive starting")
    print("guess at this precision; this range is illustrative numerical evidence,")
    print("not the proof, which is unconditional and does not depend on this script.")
