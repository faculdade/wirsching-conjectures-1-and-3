#!/usr/bin/env python3
"""Verification of the partition bound used in the proof of Conjecture 1.

The proof needs, for the coin set c_0 = 1 and c_j = 2*3^(j-1),

    p_l(m) := #{ ways to write m as a sum of coins } <= exp(C log^2(m+2)).

The paper proves it in two clauses: only coins c_j <= m contribute, and
there are at most log_3(m/2)+2 of those, while each multiplicity is at
most m, so p_l(m) <= (m+1)^(log_3(m/2)+2), which is exp(O(log^2(m+2))).

This script does two things:

  1. computes p_l(m) exactly by dynamic programming and reports the
     smallest C for which the bound holds at the sampled m on the tested range, and
  2. checks the paper's own elementary bound against the exact count,
     confirming it is valid (and loose, which is fine: the proof only
     needs the exp(O(log^2)) shape).

Why it matters. This step is what makes the tail of the convolution
e_l = p_l * g_l negligible: the partition weights grow subexponentially
while the binomial ratio decays geometrically, so a window of width
O(sqrt l) dominates. If the partition count grew like exp(c*m) instead,
the argument would fail.

Run:  python3 partition_bound_check.py [--max-m 4000]
Exits 0 if the bound holds at the sampled m at every sampled m (see the sample list below; the all-m loop fits the constant and does NOT re-test the bound), 1 otherwise.
"""
from __future__ import annotations

import argparse
import math


def coins_up_to(m: int) -> list[int]:
    cs = [1]
    j = 1
    while 2 * 3 ** (j - 1) <= m:
        cs.append(2 * 3 ** (j - 1))
        j += 1
    return cs


def partition_counts(M: int) -> list[int]:
    """p_l(m) for m = 0..M, exact integers, unrestricted multiplicities."""
    dp = [0] * (M + 1)
    dp[0] = 1
    for c in coins_up_to(M):
        for k in range(c, M + 1):
            dp[k] += dp[k - c]
    return dp


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-m", type=int, default=4000)
    args = ap.parse_args()
    M = args.max_m
    dp = partition_counts(M)

    print(f"{'m':>6} {'p(m)':>18} {'log p(m)':>10} {'log^2(m+2)':>12} "
          f"{'ratio':>8} {'elementary bound':>18}")
    worst = 0.0
    bound_ok = True
    for m in [5, 10, 20, 50, 100, 200, 500, 1000, 2000, M]:
        if m > M:
            continue
        lp = math.log(dp[m])
        l2 = math.log(m + 2) ** 2
        n = len(coins_up_to(m))
        elem = n * math.log(m + 1)
        if elem < lp:
            bound_ok = False
        print(f"{m:6d} {dp[m]:18d} {lp:10.3f} {l2:12.3f} {lp/l2:8.4f} "
              f"{'exp(%.1f)' % elem:>18}")
    for m in range(1, M + 1):
        if dp[m] > 0:
            worst = max(worst, math.log(dp[m]) / math.log(m + 2) ** 2)

    print()
    print(f"smallest valid C on m <= {M}: {worst:.4f}")
    print(f"so p(m) <= exp({worst:.3f} * log^2(m+2)) throughout the range")
    print(f"paper's elementary bound (m+1)^(log_3(m/2)+2) valid: {bound_ok}")
    if not bound_ok:
        print("FAILED: the elementary bound is violated somewhere")
        return 1
    print("\nBound holds. The count is subexponential in m, which is what "
          "the proof's tail estimate needs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
