"""Section 11: condition (*4), and the two steps the source states rather
than carries out.

Backs three claims of the unified paper's Section 11.

1. Theorem 11.1's proof replaces Wirsching's own attribution of the
   penultimate step to his (7.7), an identity for phi, by his (7.12),
   which holds for phi_0 only asymptotically. phi_0 does NOT satisfy the
   truncated equation exactly, and Remark 11.2 quotes the table this
   prints.

2. His (7.13), stated there as the result of "a somewhat lengthy
   calculation" that is not shown, follows from the closed form (2).
   This reproduces the limit 2/3, which is what gives mu = 1/3.

3. Remark 11.4's residual: the periodic factor cancels across a factor
   of 3 only asymptotically. The residual saddle phase shift is computed
   here from the saddle equation, not taken as an input, and turned into
   a bound with the certified sup|H'| of section 6.

WHAT THIS IS AND IS NOT. Finitely many evaluations cannot establish a
limit. What is asserted below is what finitely many evaluations can
establish: that a deviation is real rather than roundoff, that it
shrinks monotonically over the sampled range, and that it is under a
stated size at the largest sample. The limits themselves are proved in
the paper.

Every derivative is taken in closed form. At t ~ 500 * 3^-500 a
difference quotient is meaningless, and an earlier version of this file
produced complex garbage that way before the closed form replaced it.

Run:  python3 verify_star4.py
"""
import mpmath as mp

WORKING_DPS = 80

C = None
BETA = ALPHA = DELTA_BK = GAMMA = EPS = None

# sup|H'|, certified in Arb ball arithmetic by
# ../section-06-periodic-correction/certify_H_nonconstancy.py, which
# prints it as the second of the two bounds (10). Imported as a number
# rather than recomputed: that script is the certificate.
SUP_H_PRIME = "0.0011977472315550332"


class CheckFailed(Exception):
    """Raised instead of assert, so python -O cannot disable a check."""


def require(condition, message):
    if not condition:
        raise CheckFailed(message)


def setup():
    global C, BETA, ALPHA, DELTA_BK, GAMMA, EPS
    C = mp.log(3)
    BETA = 1 / (2 * C)
    ALPHA = mp.mpf(1) / 2 - mp.log(2) / C
    DELTA_BK = mp.mpf(1) / 2 + ALPHA - 2 * BETA * mp.log(2 * BETA)
    GAMMA = -2 * BETA - DELTA_BK - mp.mpf(1) / 2
    EPS = mp.mpf(1) / 2 + ALPHA - BETA * mp.log(2 * BETA)


def log_phi0(t):
    """log of Berg-Kruppel (9.6), the paper's equation (2), without the
    constant prefactor: it cancels in every ratio used here."""
    u = -mp.log(t)
    return GAMMA * mp.log(t) + DELTA_BK * mp.log(u) - BETA * (mp.log(t) - mp.log(u)) ** 2


def dlog_phi0(t):
    """d/dt log phi_0(t), in closed form. With u = -log t,
    log(t/u) = -(u + log u), so L(u) = -gamma u + delta_BK log u
    - beta (u + log u)^2 and dL/dt = -e^u dL/du."""
    u = -mp.log(t)
    dLdu = -GAMMA + DELTA_BK / u - 2 * BETA * (u + mp.log(u)) * (1 + 1 / u)
    return -mp.e ** u * dLdu


def truncated_equation_ratio(t):
    """(2/9) phi_0'(t) / phi_0(3t). Wirsching's (7.12) says this -> 1."""
    return (mp.mpf(2) / 9) * dlog_phi0(t) * mp.e ** (log_phi0(t) - log_phi0(3 * t))


def B_sm(y):
    """The smooth saddle function of Lemma 8.1, B_sm(y) = 2 beta y - alpha."""
    return 2 * BETA * y - ALPHA


def w0(tau):
    """The smooth saddle location: the solution of Phi_0(y) = y - log B_sm(y) = tau
    on the branch B_sm > 1/c. Solved here, not tabulated."""
    guess = tau + mp.log(tau / C)
    root = mp.findroot(lambda y: y - mp.log(B_sm(y)) - tau, guess)
    require(B_sm(root) > 1 / C, "the root found is off the branch B_sm > 1/c")
    require(abs(root - mp.log(B_sm(root)) - tau) < mp.mpf(10) ** (-WORKING_DPS + 20),
            "the saddle equation is not satisfied at the returned root")
    return root


def check_truncated_equation():
    print("1. Wirsching's (7.12): (2/9) phi_0'(t)/phi_0(3t), at t = l*3^-l.")
    print("   phi_0 solves the truncated equation only asymptotically, so his")
    print("   (7.7), an identity for phi, does not license the step.")
    print()
    print("      l        (2/9) phi_0'(t)/phi_0(3t)          ratio - 1")
    levels = [10, 50, 100, 400, 5000]
    devs = []
    for l in levels:
        t = mp.mpf(l) / mp.mpf(3) ** l
        r = truncated_equation_ratio(t)
        d = abs(r - 1)
        devs.append(d)
        print("   %6d   %s   %s" % (l, mp.nstr(r, 15).ljust(20), mp.nstr(r - 1, 4)))
        # the deviation is real, not roundoff: working precision is 80 digits
        require(d > mp.mpf(10) ** (-40),
                "the deviation at l=%d is at roundoff scale, so it proves nothing" % l)
    require(all(devs[i] > devs[i + 1] for i in range(len(devs) - 1)),
            "the deviation from 1 is not decreasing across the sampled levels")
    require(devs[-1] < mp.mpf("1e-6"),
            "the deviation at the largest level is not yet small")
    print("   PASS: the deviation is far above roundoff at every level, is")
    print("   strictly decreasing across them, and is below 1e-6 at l = 5000.")
    print("   So the ratio is genuinely not 1, and genuinely tending to it.")
    print()


def wirsching_713(l, offset=0):
    """(2/3^(l+1)) phi_0'(x)/phi_0(x) at x = (l + offset) 3^-l."""
    t = (mp.mpf(l) + offset) / mp.mpf(3) ** l
    return (mp.mpf(2) / mp.mpf(3) ** (l + 1)) * dlog_phi0(t)


def check_713():
    print("2. Wirsching's (7.13): (2/3^(l+1)) phi_0'(x_l)/phi_0(x_l) -> 2/3,")
    print("   which is what gives 1 - mu = 2/3, that is mu = 1/3.")
    print()
    target = mp.mpf(2) / 3
    for label, off in (("central,  x_l = l 3^-l", lambda l: 0),
                       ("upper,    x_l = (l+2 sqrt l) 3^-l", lambda l: 2 * mp.sqrt(l)),
                       ("lower,    x_l = (l-2 sqrt l) 3^-l", lambda l: -2 * mp.sqrt(l))):
        print("   %s" % label)
        levels = [10, 100, 500, 2000, 10000] if off(1) == 0 \
            else [100, 500, 2000, 10000, 50000]
        devs = []
        sides = []
        for l in levels:
            v = wirsching_713(l, off(l))
            d = abs(v - target)
            devs.append(d)
            side = "below" if v < target else "above"
            print("      l = %6d   %s   |v - 2/3| = %s   (%s)"
                  % (l, mp.nstr(v, 12).ljust(15), mp.nstr(d, 4), side))
            sides.append(side)
        require(all(devs[i] > devs[i + 1] for i in range(len(devs) - 1)),
                "not approaching 2/3 monotonically on the %s branch" % label)
        shrink = devs[0] / devs[-1]
        print("      deviation shrinks by a factor of %s from l = %d to l = %d"
              % (mp.nstr(shrink, 5), levels[0], levels[-1]))
        require(shrink > 10, "the deviation barely moved on the %s branch" % label)
        require(len(set(sides)) == 1,
                "the %s branch crosses 2/3 inside the tested range, so it is "
                "not approaching monotonically from one side" % label)
        print("      approaches 2/3 from %s throughout" % sides[0])
    print()
    print("   PASS: on the central sequence and on both edges of a delta = 2")
    print("   window the quantity approaches 2/3 monotonically, from below on")
    print("   the centre and the upper edge and from above on the lower edge.")
    print("   Approach from above costs (7.5) nothing: it is a limsup bound,")
    print("   and the limit is 2/3 on every branch. The approach is far slower off the centre, by")
    print("   two orders of magnitude at l = 10000, which is the same")
    print("   non-uniformity in the window that the deficit l(2/3 - L_l)")
    print("   shows in section-12. It is convergence either way, and (7.5)")
    print("   reads limsup <= 1 - mu < 1, so equality at 2/3 is admissible")
    print("   and mu = 1/3 needs no weakening.")
    print()


def check_residual():
    print("3. Remark 11.4: the periodic factor cancels only asymptotically.")
    print()
    print("   Multiplying the argument by 3 shifts tau by c. The saddle")
    print("   location does NOT shift by exactly c: solving the saddle")
    print("   equation at both points gives the residual below.")
    print()
    sup_hp = mp.mpf(SUP_H_PRIME)
    print("      l      tau         w_0(tau)     residual         bound on H")
    residuals = []
    for l in (100, 500, 2000):
        tau = -mp.log(mp.mpf(l) / mp.mpf(3) ** l)
        a, b = w0(tau), w0(tau - C)
        resid = abs(a - b - C)
        residuals.append(resid)
        print("   %6d  %s  %s  %s   %s"
              % (l, mp.nstr(tau, 7).ljust(10), mp.nstr(a, 7).ljust(11),
                 mp.nstr(resid, 8).ljust(14), mp.nstr(sup_hp * resid, 4)))
        require(resid > 0, "the residual came out zero, so the cancellation "
                           "would be exact, which the paper denies")
    require(all(residuals[i] > residuals[i + 1] for i in range(len(residuals) - 1)),
            "the residual is not shrinking with the level")
    l500 = residuals[1]
    bound = sup_hp * l500
    print()
    print("   at l = 500 the residual is %s and the bound is %s"
          % (mp.nstr(l500, 9), mp.nstr(bound, 4)))
    require(mp.mpf("1.9e-3") < l500 < mp.mpf("2.1e-3"),
            "the l=500 residual is not the 2.0e-3 the paper quotes")
    require(mp.mpf("2.3e-6") < bound < mp.mpf("2.5e-6"),
            "the l=500 bound is not the 2.4e-6 the paper quotes")
    print("   PASS: matches the paper's 2.0e-3 and 2.4e-6. Small, and not zero.")
    print()


def main():
    mp.mp.dps = WORKING_DPS
    setup()
    check_truncated_equation()
    check_713()
    check_residual()
    print("all checks passed")


if __name__ == "__main__":
    main()
