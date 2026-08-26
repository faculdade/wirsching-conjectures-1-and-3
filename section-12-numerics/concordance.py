"""Section 12: the extrapolated intercepts against the proved limit.

Nothing here is an input except the two fitted intercepts, which come
from experiment_conjecture3.py in this folder and are quoted with their
provenance below. In particular e^{H(0)} is DERIVED, not transcribed:
an earlier version of this file hardcoded it, and an audit showed the
script then passed with a value 12 percent wrong.

H(0) comes from Theorem 6.1's own decomposition,

    H(w) = L(w) - Q(w) + sum_{k>=0} log(1 - e^{-2 * 3^k * e^w}),
    L(w) = K(e^w),  K(s) = sum_{j>=1} g(2s/3^j),  g(b) = log((1-e^-b)/b),
    Q(w) = -w^2/(2c) + (1/2 - log2/c) w,   c = log 3,

at w = 0, where Q(0) = 0. This is a different route from the Fourier
series of Proposition 6.3 and from the Arb certificate of section 6,
so it is a genuine cross-check rather than a restatement; the script
confirms it against section 6's certified enclosure of
H(0) - H(log(3/2)).

Run:  python3 concordance.py
"""
import mpmath as mp

WORKING_DPS = 40

# From ../section-06-periodic-correction/certify_H_nonconstancy.py, which
# certifies this two-sided enclosure in Arb ball arithmetic.
CERT_D_LO = "-0.000377190280943987"
CERT_D_HI = "-0.000377190280943985"

# From experiment_conjecture3.py in this folder, run at its default
# --max-ell 500. Its printed fit block reads:
#   C/sqrt(ell)   L_inf=-0.618860  coeff=-0.7916  c=0.5386  max|resid|=2.14e-05
#   C/ln^2(ell)   L_inf=-0.599498  coeff=-2.1159  c=0.5491  max|resid|=5.51e-05
#   C/sqrt+1/ell  L_inf=-0.618438  coeff=-0.8067  c=0.5388  max|resid|=1.54e-05
FIT_SQRT = "-0.618860"
FIT_LN2 = "-0.599498"
FIT_SQRT_PLUS = "-0.618438"


class CheckFailed(Exception):
    """Raised instead of assert, so python -O cannot disable a check."""


def require(condition, message):
    if not condition:
        raise CheckFailed(message)


def g(b):
    return mp.log((1 - mp.e ** (-b)) / b)


def K(s):
    return mp.nsum(lambda j: g(2 * s / mp.mpf(3) ** j), [1, mp.inf])


def Q(w):
    c = mp.log(3)
    return -w ** 2 / (2 * c) + (mp.mpf(1) / 2 - mp.log(2) / c) * w


def H(w):
    s = mp.e ** w
    tail = mp.nsum(lambda k: mp.log(1 - mp.e ** (-2 * mp.mpf(3) ** k * s)), [0, mp.inf])
    return K(s) - Q(w) + tail


def main():
    mp.mp.dps = WORKING_DPS

    print("1. H(0), derived from Theorem 6.1's decomposition")
    h0 = H(mp.mpf(0))
    e_h0 = mp.e ** h0
    print("   H(0)     = %s" % mp.nstr(h0, 18))
    print("   e^{H(0)} = %s" % mp.nstr(e_h0, 18))
    print("   paper's Theorem 10.1 quotes 0.534122...")
    require(abs(e_h0 - mp.mpf("0.534122")) < mp.mpf("1e-6"),
            "the derived e^{H(0)} does not match the paper's 0.534122")

    print()
    print("2. cross-check against section 6's certified enclosure")
    d = h0 - H(mp.log(mp.mpf(3) / 2))
    lo, hi = mp.mpf(CERT_D_LO), mp.mpf(CERT_D_HI)
    print("   H(0) - H(log 3/2) here      = %s" % mp.nstr(d, 18))
    print("   section 6 certifies it in   [%s, %s]" % (CERT_D_LO, CERT_D_HI))
    require(lo <= d <= hi,
            "the decomposition route disagrees with the Arb certificate")
    print("   PASS: two independent routes to H agree inside the certified ball.")

    print()
    print("3. the fitted intercepts against the proved limit")
    log_true = h0
    fits = [("Lambda + b/sqrt(l) ", mp.mpf(FIT_SQRT)),
            ("Lambda + b'/ln^2(l)", mp.mpf(FIT_LN2))]
    print()
    print("   model                intercept    c = e^Lambda    overshoot")
    cs = []
    for name, lam in fits:
        c = mp.e ** lam
        cs.append(c)
        over = 100 * (c - e_h0) / e_h0
        print("   %s  %s   %s   %s %%"
              % (name, mp.nstr(lam, 8).ljust(10), mp.nstr(c, 8).ljust(11),
                 mp.nstr(over, 4)))
    lo_c, hi_c = min(cs), max(cs)
    print()
    print("   interval spanned: [%s, %s]" % (mp.nstr(lo_c, 7), mp.nstr(hi_c, 7)))
    require(not (lo_c <= e_h0 <= hi_c), "the interval contains the proved limit")
    print("   the proved limit lies OUTSIDE it")

    print()
    print("4. the gap against the largest within-model movement")
    movement = abs(mp.mpf(FIT_SQRT_PLUS) - mp.mpf(FIT_SQRT))
    gap = min(abs(lam - log_true) for _, lam in fits)
    ratio = gap / movement
    print("   adding a further 1/l term moved Lambda by  %s" % mp.nstr(movement, 6))
    print("   nearest gap to the true Lambda             %s" % mp.nstr(gap, 6))
    print("   ratio                                      %s" % mp.nstr(ratio, 6))
    require(ratio > 15, "the gap is not an order of magnitude beyond the movement")
    require(ratio < 20, "the ratio now exceeds 20; the paper says 19.6 and must "
                        "be updated if this changes")
    print("   PASS: %s times. The paper says 19.6." % mp.nstr(ratio, 3))
    print("   An earlier version of this script reported 'more than twenty")
    print("   times' by rounding the movement down to 4e-4 and the intercept")
    print("   to -0.6189. Both roundings favoured the claim, and it was false.")

    print()
    print("What this does NOT show: that either asymptotic expansion is wrong.")
    print("An approach Lambda_true + b l^-1/2 + o(l^-1/2) with a large")
    print("lower-order term produces exactly this. Only the fitted intercepts")
    print("over the range fitted are shown to be wrong.")
    print()
    print("all checks passed")


if __name__ == "__main__":
    main()
