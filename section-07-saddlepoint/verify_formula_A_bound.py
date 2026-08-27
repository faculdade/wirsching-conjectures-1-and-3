"""Independent audit of the explicit constant chain in Theorem 13 (uniform
saddlepoint asymptotic), proof in Section 4 of
papers/01-wirsching-conjecture3/main.tex (mirroring
notes/H-006-formula-A-proof-2.md in the main project repository).

This script evaluates the analytic upper-bound expressions from the proof at
high precision and checks the numeric inequalities the proof states, exactly
as the proof states them (it is diagnostic, not a source of any inequality:
the written proof supplies the bounds, this script confirms the arithmetic).
It also directly evaluates kappa_2, kappa_3 (Lemma 9) against numerical
differentiation of the defining function, at real and complex arguments, and
evaluates the assembled error bound E(N) at a range of N to confirm it is
strictly decreasing, below 1 at N=19, and O(N^{-1/2}) with the stated limit.
"""

import decimal

import mpmath as mp

mp.mp.dps = 60


def ell(z):
    """g(z) = log((1-e^-z)/z), the log-density kernel."""
    return mp.log((1 - mp.exp(-z)) / z)


def kappa_direct(n, z):
    """kappa_n(z) = (-1)^n z^n g^{(n)}(z), by direct numerical differentiation."""
    return (-1) ** n * z ** n * mp.diff(ell, z, n)


def kappa_closed(z):
    """Closed forms for kappa_2, kappa_3 (Section 4 of the paper)."""
    x = z / 2
    return (
        1 - x ** 2 / mp.sinh(x) ** 2,
        2 - 2 * x ** 3 * mp.cosh(x) / mp.sinh(x) ** 3,
    )


def e_func(r):
    return r ** 2 / mp.sinh(r) ** 2


def f_func(r):
    return r ** 3 * mp.cosh(r) / mp.sinh(r) ** 3


def E_bound(N, alpha=mp.mpf(1) / 3):
    """The explicit bound E(N) assembled in the proof of Theorem 13."""
    A = mp.mpf("1.4269413069")
    Vlo = N - A
    Vup = N + mp.mpf("0.1283")
    B = (2 * N + mp.mpf("10.6")) / 6
    aN = mp.power(N, -alpha)
    eps = B * aN ** 3
    e1 = 4 * mp.e ** eps * B / (mp.sqrt(2 * mp.pi) * Vlo ** mp.mpf("1.5"))
    cc = aN * mp.sqrt(Vlo)
    e2 = (2 / mp.sqrt(2 * mp.pi)) * mp.e ** (-cc ** 2 / 2) / cc
    e3 = (
        # equation (14) carries e^{3/e} itself. 3.0152 is Lemma 7.3's
        # stated rounded upper bound for it, and substituting the rounded
        # value here is what made this script disagree with the paper's
        # printed E(18) and E(19) in the sixth significant digit.
        2 * mp.e ** (3 / mp.e) * mp.sqrt(Vup / (2 * mp.pi))
        * (1 + aN ** 2) ** (-(N - 2) / mp.mpf(2)) / (aN * (N - 2))
    )
    return e1 + e2 + e3


def main():
    print("Lemma 1: sector bounds M_2 <= 3, M_3 <= 7.657 (sup over |Im w| <= Re w)")
    M2, M3 = mp.mpf(0), mp.mpf(0)
    for b in [mp.mpf("0.1"), mp.mpf(1), mp.mpf(3), mp.mpf("4.076"), mp.mpf(10), mp.mpf(30)]:
        for u in [mp.mpf(0), mp.mpf("0.5"), mp.mpf(1)]:
            z = b * (1 + 1j * u)
            k2c, k3c = kappa_closed(z)
            k2d, k3d = kappa_direct(2, z), kappa_direct(3, z)
            assert abs(k2c - k2d) < mp.mpf("1e-40"), (b, u, "kappa_2 mismatch")
            assert abs(k3c - k3d) < mp.mpf("1e-40"), (b, u, "kappa_3 mismatch")
            M2, M3 = max(M2, abs(k2c)), max(M3, abs(k3c))
    print(f"  closed-form kappa_2, kappa_3 match numerical differentiation to 1e-40")
    print(f"  sampled sup |kappa_2| = {mp.nstr(M2, 8)} (bound 3)")
    print(f"  sampled sup |kappa_3| = {mp.nstr(M3, 8)} (bound 7.657)")
    assert M2 <= 3 and M3 <= mp.mpf("7.657")

    print("\nLemma 1': local bound at |w|<=2: |kappa_2|<=0.114|w|^2, |kappa_3|<=0.0119|w|^4")
    # The supremum of |kappa_2(w)/w^2| and |kappa_3(w)/w^4| over the closed
    # disc |w| <= 2 is attained on the boundary, by the maximum modulus
    # principle applied to each (both are analytic on Re w > 0 and the
    # relevant sector). Sampling interior points only, as an earlier version
    # of this script did, tests the bound where it is slack by 25 percent.
    print("  boundary |w| = 2, scanned by argument:")
    NGRID = 4000
    worst2 = worst3 = mp.mpf(0)
    arg2 = arg3 = None
    for i in range(NGRID + 1):
        theta = mp.pi * mp.mpf(i) / NGRID - mp.pi / 2      # Re w > 0 half
        theta = max(min(theta, mp.pi / 2 - mp.mpf("1e-6")), -mp.pi / 2 + mp.mpf("1e-6"))
        w = 2 * mp.e ** (1j * theta)
        k2b, k3b = kappa_closed(w)
        r2 = abs(k2b) / abs(w) ** 2
        r3 = abs(k3b) / abs(w) ** 4
        if r2 > worst2:
            worst2, arg2 = r2, theta
        if r3 > worst3:
            worst3, arg3 = r3, theta
    print("    max |kappa_2|/|w|^2 = %s at arg w = %s"
          % (mp.nstr(worst2, 9), mp.nstr(arg2, 6)))
    print("    max |kappa_3|/|w|^4 = %s at arg w = %s"
          % (mp.nstr(worst3, 9), mp.nstr(arg3, 6)))
    assert worst2 <= mp.mpf("0.114"), "the 0.114 of Lemma 7.1 fails on |w| = 2"
    assert worst3 <= mp.mpf("0.0119"), "the 0.0119 of Lemma 7.1 fails on |w| = 2"
    print("    both constants hold on the boundary, which is where the sup sits")
    print("  interior spot checks:")
    for b in [mp.mpf("0.01"), mp.mpf("0.5"), mp.mpf(2)]:
        z = b * (1 + 1j * mp.mpf("0.9"))
        if abs(z) > 2:
            continue
        k2c, k3c = kappa_closed(z)
        assert abs(k2c) <= mp.mpf("0.114") * abs(z) ** 2
        assert abs(k3c) <= mp.mpf("0.0119") * abs(z) ** 4
    print("  holds at sampled points")

    print("\nLemma 2: variance constant A = sum_{m>=0} e(3^m/2)")
    A = mp.nsum(lambda m: e_func(mp.mpf(3) ** m / 2), [0, mp.inf])
    print(f"  A = {mp.nstr(A, 15)}  (stated: 1.4269413069)")
    assert abs(A - mp.mpf("1.4269413069")) < mp.mpf("1e-9")

    print("\nLemma 4: tail constant F = sum_{m>=0} f(3^m/2)")
    F = mp.nsum(lambda m: f_func(mp.mpf(3) ** m / 2), [0, mp.inf])
    print(f"  F = {mp.nstr(F, 15)}")
    assert abs(F - mp.mpf("1.8635631489")) < mp.mpf("1e-9")

    # Lemma 7.4 no longer quotes F itself. It bounds Sigma by elementary
    # upper estimates of its first three terms plus an analytic tail, and
    # the chain below is exactly what the paper prints. Each link is
    # checked here rather than taken on the paper's word.
    print("\nLemma 7.4's chain to S < 2N + 10.6")
    heads = [(mp.mpf(1) / 2, "0.997"), (mp.mpf(3) / 2, "0.823"),
             (mp.mpf(9) / 2, "0.046")]
    for r, bound in heads:
        v = f_func(r)
        print(f"  f({mp.nstr(r,3)}) = {mp.nstr(v, 10)} < {bound}")
        assert v < mp.mpf(bound), f"f({r}) exceeds its stated bound {bound}"
    # the identity the paper derives, and the factor it bounds for r >= 5
    for r in [mp.mpf(1) / 2, mp.mpf(5), mp.mpf("13.5")]:
        lhs = f_func(r)
        rhs = 4 * r ** 3 * mp.e ** (-2 * r) * (1 + mp.e ** (-2 * r)) / (1 - mp.e ** (-2 * r)) ** 3
        assert abs(lhs - rhs) < mp.mpf("1e-30") * max(1, abs(lhs)), f"identity fails at r={r}"
    fac5 = (1 + mp.e ** (-10)) / (1 - mp.e ** (-10)) ** 3
    print(f"  (1+e^-2r)/(1-e^-2r)^3 at r=5: {mp.nstr(fac5, 10)} < 1.00025")
    assert fac5 < mp.mpf("1.00025")
    tail = mp.mpf("4.001") * mp.nsum(
        lambda m: (mp.mpf(3) ** m / 2) ** 3 * mp.e ** (-(mp.mpf(3) ** m)), [3, mp.inf])
    print(f"  tail from m=3: {mp.nstr(tail, 8)} < 2e-8")
    assert tail < mp.mpf("2e-8")
    # The purely decimal links are checked in exact decimal arithmetic:
    # 0.997 + 0.823 + 0.046 is EXACTLY 1.866, and rounding it through
    # binary floating point is how one convinces oneself otherwise.
    D = decimal.Decimal
    decimal.getcontext().prec = 50
    head_sum = D("0.997") + D("0.823") + D("0.046")
    sigma_bound = head_sum + D("2e-8")
    print(f"  head bounds sum to exactly {head_sum}, so Sigma < {sigma_bound}")
    assert head_sum == D("1.866"), "the three head bounds no longer sum to 1.866"
    assert sigma_bound == D("1.86600002")
    assert F < mp.mpf("1.86600002"), "the stated Sigma bound does not hold"
    assert D(2).sqrt() ** 5 * D("1.86600002") < D("10.5557")
    assert D("10.5557") + D("0.01704") < D("10.6")
    print("  2^(5/2)*1.86600002 < 10.5557, and 10.5557 + 0.01704 < 10.6: S < 2N + 10.6")

    print("\nE(N): the assembled error bound, Theorem 13")
    print(f"{'N':>8}{'E(N)':>16}{'sqrt(N)*E(N)':>16}")
    prev = None
    for N in [19, 20, 30, 50, 100, 1000, 10 ** 4, 10 ** 6]:
        e = E_bound(mp.mpf(N))
        print(f"{N:>8}{mp.nstr(e, 10):>16}{mp.nstr(mp.sqrt(N) * e, 10):>16}")
        if prev is not None:
            assert e < prev, f"E(N) not decreasing at N={N}"
        prev = e
    e18 = E_bound(mp.mpf(18))
    e19 = E_bound(mp.mpf(19))
    assert e18 > 1, "E(18) should be > 1 (N_0=19 is sharp)"
    assert e19 < 1, "E(19) should be < 1"
    limit = mp.sqrt(mp.mpf(10) ** 8) * E_bound(mp.mpf(10) ** 8)
    print(f"\nE(18) = {mp.nstr(e18, 6)} > 1, E(19) = {mp.nstr(e19, 6)} < 1: confirmed, N_0=19")
    print(f"sqrt(N)*E(N) at N=10^8: {mp.nstr(limit, 8)}  (claimed limit: 0.742358)")
    assert abs(limit - mp.mpf("0.742358")) < mp.mpf("0.01")

    # The paper NO LONGER claims either of the two properties below; the
    # sentence that did was cut, because nothing downstream used it and it
    # was a computational claim outside the three the paper declares. They
    # are kept here because they are true, cheap, and useful to anyone
    # reading E(N), and they are labelled as extra rather than as backing
    # any statement in the manuscript.
    print("\nExtra, backing no claim in the paper: E is strictly decreasing")
    print("and sqrt(N) E(N) <= 4.18 on every integer 19 <= N <= 5000.")
    prev = None
    worst = mp.mpf(0)
    worst_N = None
    for N in range(19, 5001):
        e = E_bound(mp.mpf(N))
        if prev is not None:
            assert e < prev, "E is not strictly decreasing at N = %d" % N
        prev = e
        sN = mp.sqrt(N) * e
        if sN > worst:
            worst, worst_N = sN, N
    print("  strictly decreasing across all 4982 integers: yes")
    print("  max sqrt(N) E(N) = %s at N = %d" % (mp.nstr(worst, 10), worst_N))
    assert worst <= mp.mpf("4.18"), "sqrt(N) E(N) exceeds 4.18 at N = %d" % worst_N

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
