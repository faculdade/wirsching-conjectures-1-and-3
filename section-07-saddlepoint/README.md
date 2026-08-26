# Section 7: the uniform saddlepoint approximation

Backs Theorem 7.5 and the constant chain behind it.

| Script | Backs |
|---|---|
| `verify_formula_A_bound.py` | the explicit error budget `E(N)`, its constants, and the threshold `N >= 19` at which `E(N) < 1` |

Run: `python3 verify_formula_A_bound.py`

The theorem's content is the uniformity in `rho` over `[1,3)`, not the
size of the bound. This script audits the chain of constants entering
`E(N)` and reproduces `E(18) >= 1 > E(19)`, which is where the theorem's
hypothesis comes from.
