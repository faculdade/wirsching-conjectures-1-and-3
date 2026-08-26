"""Audit the exact P-versus-phi_0 constant calculation for H-006.

The proof is recorded in notes/H-006.md.  This script only evaluates the
closed formula at increasing tau, independently checking its numerical values.

Convention matters:
  * phi0_bare omits Berg--Krueppel (9.6)'s constant prefactor;
  * phi0_BK includes (2*beta)^epsilon/sqrt(2*pi).

For Q(w)=-w^2/(2h)+a*w, h=log(3), a=1/2-log(2)/h, the exact smooth
saddle is expressed through W_{-1}.  The derivation proves

  P(tau)-log(phi0_bare(exp(-tau))) -> log((2*beta)^epsilon/sqrt(2*pi)),
  P(tau)-log(phi0_BK(exp(-tau))) -> 0.
"""

import mpmath as mp


mp.mp.dps = 100
h = mp.log(3)
a = mp.mpf("0.5") - mp.log(2) / h
beta = 1 / (2 * h)
delta = mp.mpf("0.5") + a - 2 * beta * mp.log(2 * beta)
epsilon = mp.mpf("0.5") + a - beta * mp.log(2 * beta)
gamma = -2 * beta - delta - mp.mpf("0.5")
const_bare = epsilon * mp.log(2 * beta) - mp.log(2 * mp.pi) / 2


def smooth_saddle(tau):
    """Return w_0 and B=w_0/h-a, solving tau=w_0-log(B) exactly."""
    B = -mp.lambertw(-h * mp.exp(h * a - tau), -1) / h
    return h * (B + a), B


def P(tau):
    w, B = smooth_saddle(tau)
    return -w**2 / (2 * h) + a * w + B + w - mp.log(2 * mp.pi * (B - 1 / h)) / 2


def log_phi0_bare(tau):
    return -gamma * tau + delta * mp.log(tau) - beta * (tau + mp.log(tau)) ** 2


if __name__ == "__main__":
    print("Berg--Krueppel parameters for a=3, lambda=2/3")
    for name, value in (("beta", beta), ("gamma", gamma), ("delta", delta),
                        ("epsilon", epsilon), ("const_bare", const_bare)):
        print(f"{name:>11} = {mp.nstr(value, 40)}")
    print("\n tau       P-log(phi0_bare)       error from const_bare     P-log(phi0_BK)")
    for tau in (10, 100, 1_000, 10_000, 1_000_000, 100_000_000):
        diff = P(tau) - log_phi0_bare(tau)
        print(f"{tau:>10}  {mp.nstr(diff, 18):>22}  {mp.nstr(diff-const_bare, 18):>22}"
              f"  {mp.nstr(diff-const_bare, 18):>18}")
