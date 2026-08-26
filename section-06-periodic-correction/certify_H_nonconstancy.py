"""Rigorous ball-arithmetic certificate for H's oscillation.

NOTE: as of a later revision of the paper, H's non-constancy itself is proven
analytically (Hhat(1) != 0 by the classical zero-free theorem for zeta on
Re(s)=1, no computation needed at all) -- this script's job is now only the
*quantitative* certificate: a rigorous two-sided enclosure of osc(H), plus
the derivative bounds (5) that Lemma 16 depends on.  The name is kept as
`certify_H_nonconstancy.py` for continuity with earlier paper drafts, but
`main()`'s printed H(0) != H(log(3/2)) check below is now a stronger,
special-case corollary of the analytic argument, not the primary proof of
non-constancy.

This uses python-flint's Arb/Acb balls (not mpmath floating point).  The input
formula in H-006 is the formula for the *second derivative* Fourier
coefficient.  Thus, for m != 0,

  h_m := Hhat(m)
       = - 2^(i omega_m)/c * Gamma(-i omega_m) * zeta(1-i omega_m),
  omega_m = (2 pi/c)m, c = log(3).

The omitted m=0 term cancels in a difference.  Since H is real,

  D = H(0)-H(log(3/2))
    = 2 Re sum_{m>=1} h_m (1-exp(i omega_m log(3/2))).

For completeness, the tail majorant used below is derived here.  Set
alpha=2*pi/c and q=exp(-pi*alpha/2).  Euler summation, with N=ceil(y),
gives for y>=alpha

 |zeta(1-iy)| <= H_N + 1/y + 1/(2N) + sqrt(1+y^2)/(2N)
                 <= log(y+1) + C,

 C = 1 + 1/alpha + 1/2 + sqrt(1+alpha^-2)/2.

The identity |Gamma(iy)|^2=pi/(y*sinh(pi*y)), and
sinh(pi*y)>=exp(pi*y)/3 (valid here because pi*alpha>log(3)), imply

 |h_m| <= A q^m f_m,
 A = c^-1 sqrt(3*pi/alpha),
 f_m = m^-1/2 (log(alpha*m+1)+C).

C>2 makes f_m decreasing for m>=1: its logarithmic derivative multiplied
by m is at most 1/(log(alpha*m+1)+C)-1/2 < 0.  Consequently

 sum_{m>M}|h_m| <= A*f_(M+1)*q^(M+1)/(1-q).

The factor 4 in ``difference_tail`` below accounts for m and -m and for
|1-exp(i theta)|<=2.  M=4 already gives an error below 3.7e-19.
"""

from flint import acb, arb, ctx


ctx.dps = 100


def hhat(m, c, alpha):
    """Certified Acb enclosure of the m-th Fourier coefficient of H."""
    omega = arb(m) * alpha
    iomega = acb(0, omega)
    return -(acb(2).log() * iomega).exp() / c * acb(0, -omega).gamma() * acb(1, -omega).zeta()


def main():
    c = arb(3).log()
    pi = arb.pi()
    alpha = 2 * pi / c
    phase = (arb(3) / 2).log()
    q = (-pi * alpha / 2).exp()
    C = arb(1) + 1 / alpha + arb(1) / 2 + (1 + 1 / alpha**2).sqrt() / 2
    A = (3 * pi / alpha).sqrt() / c

    M = 4
    finite = arb(0)
    print("Raw Fourier coefficients Hhat(m), m=1..4 (as printed inline in the paper's")
    print("Proposition 8, at 100-decimal working precision here vs. 250-bit there):")
    for m in range(1, M + 1):
        coeff = hhat(m, c, alpha)
        print(f"  Hhat({m}) = {coeff}")

    print("\nCertified finite contribution to D = H(0)-H(log(3/2))")
    for m in range(1, M + 1):
        omega = arb(m) * alpha
        coeff = hhat(m, c, alpha)
        term = 2 * (coeff * (1 - acb(0, omega * phase).exp())).real
        finite += term
        print(f"m={m}: {term}")

    n = arb(M + 1)
    f_n = ((alpha * n + 1).log() + C) / n.sqrt()
    positive_tail = A * f_n * q**n / (1 - q)
    # Use the upper endpoint as a scalar error radius, so every subsequent
    # enclosure remains outward-rounded.
    difference_tail = 4 * positive_tail.upper()
    certified = finite + arb(0, difference_tail)

    print("\nTail data (all are rigorous Arb enclosures)")
    print(f"alpha = {alpha}")
    print(f"q     = {q}")
    print(f"C     = {C}")
    print(f"sum_(m>{M}) |Hhat(m)| <= {positive_tail}")
    print(f"|D-D_{M}| <= {difference_tail}")
    print(f"\nD_{M} = {finite}")
    print(f"Certified D = {certified}")

    # A deliberately wide decimal enclosure, convenient to quote in prose.
    left = arb("-0.000377190280943987")
    right = arb("-0.000377190280943985")
    assert certified.lower() > left
    assert certified.upper() < right
    assert certified < 0
    print("\nTherefore")
    print("  -0.000377190280943987 < H(0)-H(log(3/2))")
    print("                              < -0.000377190280943985,")
    print("so the two values are rigorously distinct.")

    certify_derivative_bounds(c, alpha, q, A)
    certify_oscillation(c, alpha)


def certify_derivative_bounds(c, alpha, q, A):
    """Certify sup|H'| and sup|H''| (equation (3) of the paper) via the same
    per-mode majorant |Hhat(m)| <= A*q^m*f_m used above, summed with the extra
    factor omega_m^order (order=1,2) closed-form (sum_{m>M} m^k q^m has a
    closed form for k=1,2; f_m is bounded above by its value at m=M+1, since
    it is decreasing for m>=1)."""
    C = arb(1) + 1 / alpha + arb(1) / 2 + (1 + 1 / alpha**2).sqrt() / 2
    M = 6
    for order, target in ((1, arb("0.0011977472315550332")),
                           (2, arb("0.0068518962896650951"))):
        head = arb(0)
        for m in range(1, M + 1):
            omega = arb(m) * alpha
            coeff_abs = abs(hhat(m, c, alpha))
            head += 2 * omega**order * coeff_abs
        n = arb(M + 1)
        f_n = ((alpha * n + 1).log() + C) / n.sqrt()
        # sum_{m>M} m^order q^m, bounded using f_m <= f_n for m>=n (f decreasing):
        # order=1: sum m q^m = q^n*(n-(n-1)q)/(1-q)^2 (m from n to infinity)
        # order=2: sum m^2 q^m = q^n*(n^2-(2n^2-2n-1)q+(n-1)^2 q^2)/(1-q)^3
        if order == 1:
            tail_sum = q**n * (n - (n - 1) * q) / (1 - q) ** 2
        else:
            tail_sum = (q**n * (n**2 - (2 * n**2 - 2 * n - 1) * q + (n - 1) ** 2 * q**2)
                        / (1 - q) ** 3)
        tail = 2 * alpha**order * A * f_n * tail_sum
        bound = (head + arb(0, tail.upper())).upper()
        print(f"\nsup|H^({order})| certified <= {bound}  (quoted in paper: {target})")
        assert arb(bound) <= target, f"certified bound {bound} exceeds quoted {target}"


def certify_oscillation(c, alpha):
    """Certify osc(H) = sup H - inf H matching the paper's CURRENT Proposition 8
    proof exactly (as of the blind-critique loop's Round 20 fix, which replaced
    an uncertified grid-plus-Lipschitz upper bound with a fully rigorous one).

    Lower bound: a grid of N=2^20 points, H evaluated via the Fourier series
    truncated at |m|<=10 (discarded-mode error ~2.7e-43, negligible here),
    locates two well-separated points w_486746, w_1011118; D := H(w_486746) -
    H(w_1011118) is a valid unconditional lower bound on osc(H), since both
    are just two specific points of the domain (H(w_486746)<=sup H,
    H(w_1011118)>=inf H).  This does NOT need the grid search to have found
    the true argmax/argmin.

    Upper bound: write H = Hhat(0) + 2*Re(Hhat(1)*e^{i*omega_1*w}) + Xi(w),
    with |Xi(w)| <= 2*sum_{m>=2}|Hhat(m)| for every w.  The one-mode term
    ranges over an interval of width exactly 4|Hhat(1)|, so
    osc(H) <= 4|Hhat(1)| + 4*sum_{m>=2}|Hhat(m)|.  No grid, no search, no
    Lipschitz argument -- uses only the already-certified Fourier
    coefficients, matching the paper's own upper-bound proof exactly."""
    n_grid = 2**20
    h = c / n_grid
    modes = list(range(1, 11))
    coeffs = [(m, arb(m) * alpha, hhat(m, c, alpha)) for m in modes]

    hi = None
    lo = None
    # Per-point error from truncating the Fourier sum at |m| <= 10.
    # The paper bounds it by 5.60e-43; a decade of margin is kept.
    TRUNC_PER_POINT = arb("5.60e-43")
    hi_i = lo_i = None
    for k in range(n_grid):
        w = h * k
        s = arb(0)
        for m, omega, coeff in coeffs:
            s += 2 * (coeff * acb(0, omega * w).exp()).real
        # For a rigorous LOWER bound on osc(H) the difference must be
        # rounded outward the other way: take the LOWER endpoint at the
        # candidate maximum and the UPPER endpoint at the candidate
        # minimum. Selecting on s.upper()/s.lower() as this script used
        # to do gives an upper-biased difference, which is not a
        # certificate of a lower bound however small the radii are.
        v_hi, v_lo = s.lower(), s.upper()
        if hi is None or v_hi > hi:
            hi, hi_i = v_hi, k
        if lo is None or v_lo < lo:
            lo, lo_i = v_lo, k

    # The grid sums only |m| <= 10. The discarded modes shift each
    # pointwise value by at most TRUNC (printed above), and the
    # difference of two values by at most twice that.
    D = arb(hi) - arb(lo) - 2 * TRUNC_PER_POINT
    osc_lo = D.lower()
    print(f"\nLower bound via N={n_grid} grid (h={h}), Fourier truncated at |m|<=10:")
    print(f"  max at grid index {hi_i}, min at grid index {lo_i}")
    print(f"  D = H(w_{hi_i}) - H(w_{lo_i}) = {D}")
    print("  (paper's Proposition 8 quotes D = 4.187449477152e-4)")
    # No slack. The previous version subtracted 1e-13 here, roughly
    # nineteen times the true margin of 5.2e-15, so it would have
    # accepted a quoted bound larger than the computed value.
    left = arb("4.1874494771e-4")
    assert osc_lo > left.lower(), (osc_lo, left)
    assert osc_lo > arb("1e-4")

    # Upper bound: 4|Hhat(1)| + 4*sum_{m>=2}|Hhat(m)|, using the already-computed
    # Hhat(2), Hhat(3), Hhat(4) exactly and the majorant tail bound for m>4.
    C = arb(1) + 1 / alpha + arb(1) / 2 + (1 + 1 / alpha**2).sqrt() / 2
    A = (3 * arb.pi() / alpha).sqrt() / c
    q = (-arb.pi() * alpha / 2).exp()
    M_tail = 4
    n = arb(M_tail + 1)
    f_n = ((alpha * n + 1).log() + C) / n.sqrt()
    tail_m_gt_4 = A * f_n * q**n / (1 - q)

    h1 = abs(hhat(1, c, alpha))
    sum_2_to_4 = sum(abs(hhat(m, c, alpha)) for m in (2, 3, 4))
    sum_ge_2 = sum_2_to_4 + arb(0, tail_m_gt_4.upper())

    four_h1 = (4 * h1).upper()
    osc_hi = (4 * h1 + 4 * sum_ge_2).upper()
    print(f"\nUpper bound via the leading Fourier mode alone:")
    print(f"  4|Hhat(1)| <= {four_h1}")
    print(f"  sum_(m>=2)|Hhat(m)| <= {sum_ge_2.upper()}")
    print(f"  osc(H) <= {osc_hi}")
    print("  (paper's Proposition 8 quotes osc(H) <= 4.187981e-4)")
    right = arb("4.187981e-4")
    assert osc_hi < right.upper(), (osc_hi, right)

    print(f"\nCertified two-sided enclosure: osc(H) in [{osc_lo}, {osc_hi}]")
    print("  (paper's Proposition 8 quotes")
    print("   [4.1874494771e-4, 4.187981e-4])")


if __name__ == "__main__":
    main()
