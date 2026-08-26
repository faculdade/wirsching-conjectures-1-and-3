# Section 6: the periodic correction H

Backs Proposition 6.5 (non-constancy and certified oscillation) and the
derivative bounds of equation (10).

| Script | Backs |
|---|---|
| `certify_H_nonconstancy.py` | the rigorous interval-arithmetic enclosure of `osc(H)` and of `sup|H'|`, `sup|H''|` |

Run: `python3 certify_H_nonconstancy.py`

Non-constancy itself is proved analytically in the paper, from the
classical zero-free theorem for zeta on `Re s = 1`, with no computation:
`Hhat(1) != 0` follows directly. What this script certifies is the
*quantitative* part, the oscillation enclosure and the derivative
bounds, which the paper does use as computer-assisted results.

Arb ball arithmetic through python-flint, at a working precision stated
at the top of the script and higher than the 250-bit precision quoted
inline in Proposition 6.5's proof. The script prints its input balls and
its final enclosure, so the printed run output is the certificate.
