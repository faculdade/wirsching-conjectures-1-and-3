#!/usr/bin/env python3
"""Finite checks for the generating-function proof of Wirsching Conjecture 1."""

from __future__ import annotations

from fractions import Fraction
from math import comb, exp, isqrt, log


log_ = log


def coin_values(ell: int) -> list[int]:
    return [1] + [2 * 3 ** (j - 1) for j in range(1, ell + 1)]


def coefficients_unbounded(coins: list[int], degree: int) -> list[int]:
    out = [0] * (degree + 1)
    out[0] = 1
    for coin in coins:
        for k in range(coin, degree + 1):
            out[k] += out[k - coin]
    return out


def coefficients_bounded(coins: list[int], degree: int) -> list[int]:
    out = [0] * (degree + 1)
    out[0] = 1
    for capacity in coins:
        nxt = [0] * (degree + 1)
        for total, count in enumerate(out):
            if count:
                for add in range(min(capacity - 1, degree - total) + 1):
                    nxt[total + add] += count
        out = nxt
    return out


def convolve_at(a: list[int], b: list[int], k: int) -> int:
    return sum(a[m] * b[k - m] for m in range(k + 1))


def main() -> None:
    print("ell max_degree identity tail_fraction_at_k=ell")
    for ell in range(2, 13):
        degree = 3 * ell
        coins = coin_values(ell)
        p = coefficients_unbounded(coins, degree)
        q = coefficients_bounded(coins, degree)
        for k in range(degree + 1):
            assert convolve_at(p, q, k) == comb(k + ell, ell)

        k = ell
        cutoff = max(1, isqrt(ell))
        denominator = comb(k + ell, ell)
        tail = sum(p[m] * q[k - m] for m in range(cutoff, k + 1))
        tail_fraction = Fraction(tail, denominator)
        print(ell, degree, "ok", f"{float(tail_fraction):.8f}")

    # A direct check of the elementary subexponential bound for p_infty(m).
    degree = 1000
    coins = [1]
    while coins[-1] <= degree:
        coins.append(2 if len(coins) == 1 else 3 * coins[-1])
    coins = [c for c in coins if c <= degree]
    p = coefficients_unbounded(coins, degree)
    # The constant must be fixed INDEPENDENTLY of the data. Setting it to
    # the maximum of the very ratios then tested makes the inequality true
    # by construction, for any sequence whatever: this script used to do
    # that, and it passed with p[m] = 2**m substituted in.
    #
    # equation (5)'s own elementary bound is p_l(m) <= (m+1)^(log_3(m/2)+2),
    # and 1.36 log^2(m+2) dominates its logarithm for every m >= 1
    # (see ../section-04-conjecture-1/verify_lemma_constants.py, check 1).
    C_FIXED = 1.36
    for m in range(2, degree + 1):
        elementary = (log_(m / 2) / log_(3) + 2) * log_(m + 1)
        if not log_(p[m]) <= elementary + 1e-12:
            raise AssertionError(
                "equation (5)'s elementary bound fails at m=%d" % m)
        if not p[m] <= exp(C_FIXED * log(m + 2) ** 2):
            raise AssertionError(
                "the subexponential form with C=%.2f fails at m=%d" % (C_FIXED, m))
    observed = max(log(max(1, p[m])) / (log(m + 2) ** 2)
                   for m in range(2, degree + 1))
    print(f"finite subexponential check through m={degree}: "
          f"fixed C={C_FIXED:.2f}, largest observed ratio={observed:.6f}")
    if not observed < C_FIXED:
        raise AssertionError("the observed ratio reached the fixed constant")


if __name__ == "__main__":
    main()
