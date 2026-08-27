# Section 7: the uniform saddlepoint approximation

Backs Theorem 7.5 and the constant chain behind it.

| Script | Backs |
|---|---|
| `verify_formula_A_bound.py` | the explicit error budget `E(N)`, its constants, and the threshold `N >= 19` at which `E(N) < 1` |

Run: `python3 verify_formula_A_bound.py`

The theorem's content is the uniformity in `rho` over `[1,3)`, not the
size of the bound. This script audits the chain of constants entering
`E(N)` and reproduces `E(18) = 1.00207... >= 1 > E(19) = 0.957056...`,
which is where the theorem's
hypothesis comes from.

It also checks Lemma 7.4's chain to `S < 2N + 10.6` link by link: the
three head bounds `f(1/2) < 0.997`, `f(3/2) < 0.823`, `f(9/2) < 0.046`;
the identity `f(r) = 4 r^3 e^{-2r} (1+e^{-2r})/(1-e^{-2r})^3` and its
last factor below `1.00025` at `r = 5`; the tail from `m = 3` below
`2e-8`; and then the purely decimal links in exact decimal arithmetic,
because `0.997 + 0.823 + 0.046` is exactly `1.866` and rounding that
through binary floating point is how one talks oneself into
`Sigma < 1.866`, which does not follow. It follows that
`Sigma < 1.86600002`.
