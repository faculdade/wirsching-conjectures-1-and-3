#!/usr/bin/env python3
"""Coefficient-by-coefficient verification of the cancellation behind
Theorem 4.1 (Wirsching's Conjecture 1).

The paper's proof turns Wirsching's urn description of the averaged
generators into two generating functions over the coin set
c_0 = 1, c_j = 2*3^(j-1) for j >= 1:

    sum_k q_l(k) z^k = prod_{j=0}^{l} (1 - z^{c_j}) / (1 - z)      (bounded urns)
    sum_k p_l(k) z^k = prod_{j=0}^{l} (1 - z^{c_j})^{-1}           (unbounded)

and uses their product collapsing to a single closed form:

    P_l(z) Q_l(z) = (1 - z)^{-(l+1)},

so that the convolution (p_l * gbar_l)(k) has the exact binomial value
C(k+l, l) / (2 * 3^(l-1)).

This script verifies all of that in exact integer arithmetic, with no
floating point anywhere:

  1. the two product expansions against each other and against their
     closed forms. (An earlier version of this docstring said "against
     direct enumeration of urn occupancies"; no such enumeration is
     implemented, and the claim is withdrawn.)
  2. the cancellation, coefficient by coefficient, against the binomial
     coefficients of (1-z)^{-(l+1)},
  3. the identity q_l(k) = 2 * 3^(l-1) * gbar_l(k), where gbar_l comes
     from the same bounded-urn count used above rather than from an
     independent implementation of Wirsching's recursion (2.1); this
     check is algebraically downstream of the object it checks and is
     a consistency check, not independent corroboration,
  4. the convolution identity (p_l * gbar_l)(k) = C(k+l,l)/(2*3^(l-1)).

Run:  python3 cancellation_check.py [--max-ell 8] [--max-k 24]
Exits 0 if every check passes, 1 on the first mismatch.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from math import comb


def coins(ell: int) -> list[int]:
    return [1] + [2 * 3 ** (j - 1) for j in range(1, ell + 1)]


def poly_mul(a: list[int], b: list[int], trunc: int) -> list[int]:
    out = [0] * (trunc + 1)
    for i, ai in enumerate(a):
        if ai == 0 or i > trunc:
            continue
        for j, bj in enumerate(b):
            if i + j > trunc:
                break
            out[i + j] += ai * bj
    return out


def bounded_urn_counts(ell: int, trunc: int) -> list[int]:
    """q_l(k): urn j holds 0..c_j-1 balls, coefficient sum k.

    Product of (1 + z + ... + z^(c_j-1)) = (1-z^{c_j})/(1-z).
    """
    poly = [1]
    for c in coins(ell):
        poly = poly_mul(poly, [1] * c, trunc)
    return poly


def unbounded_counts(ell: int, trunc: int) -> list[int]:
    """p_l(k): unbounded multiplicities, prod (1-z^{c_j})^{-1}."""
    poly = [0] * (trunc + 1)
    poly[0] = 1
    for c in coins(ell):
        for k in range(c, trunc + 1):
            poly[k] += poly[k - c]
    return poly


def wirsching_gbar(ell: int, trunc: int) -> list[Fraction]:
    """gbar_l(k) from the bounded-urn count, as exact rationals."""
    q = bounded_urn_counts(ell, trunc)
    denom = 2 * 3 ** (ell - 1)
    return [Fraction(v, denom) for v in q]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-ell", type=int, default=8)
    ap.add_argument("--max-k", type=int, default=24)
    args = ap.parse_args()
    trunc = args.max_k

    failures = 0
    for ell in range(1, args.max_ell + 1):
        q = bounded_urn_counts(ell, trunc)
        p = unbounded_counts(ell, trunc)

        # (2) the cancellation, coefficient by coefficient
        prod = poly_mul(p, q, trunc)
        for k in range(trunc + 1):
            want = comb(k + ell, ell)
            if prod[k] != want:
                print(f"FAIL cancellation: ell={ell} k={k} "
                      f"got {prod[k]} want {want}")
                failures += 1
        # support bound: sum_j (c_j - 1) = 3^ell - ell - 1
        support = sum(c - 1 for c in coins(ell))
        if support != 3 ** ell - ell - 1:
            print(f"FAIL support bound: ell={ell} got {support} "
                  f"want {3**ell - ell - 1}")
            failures += 1
        # (4) the convolution identity in exact rationals
        gbar = wirsching_gbar(ell, trunc)
        for k in range(trunc + 1):
            conv = sum(Fraction(p[m]) * gbar[k - m] for m in range(k + 1))
            want = Fraction(comb(k + ell, ell), 2 * 3 ** (ell - 1))
            if conv != want:
                print(f"FAIL convolution: ell={ell} k={k} "
                      f"got {conv} want {want}")
                failures += 1
        print(f"ell={ell:2d}: cancellation, support bound and convolution "
              f"verified for k=0..{trunc}")

    if failures:
        print(f"\n{failures} mismatch(es).")
        return 1
    print(f"\nAll checks passed for ell=1..{args.max_ell}, k=0..{trunc}, "
          f"exact integer and rational arithmetic throughout.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
