# Section 6: the periodic correction H

Backs Proposition 6.5 (non-constancy, certified oscillation, and the
enclosure of `H(0)`) and the derivative bounds of equation (10).

| Script | Backs |
|---|---|
| `certify_H_nonconstancy.py` | the rigorous interval-arithmetic enclosures of `osc(H)` and of `H(0)`, and of `sup|H'|`, `sup|H''|` |

The script has four certification stages, in the order the paper's
proof takes them: `main` certifies `H(0) - H(log(3/2))`, which is what
makes `H` provably non-constant quantitatively;
`certify_derivative_bounds` does equation (10);
`certify_oscillation` does the two-sided enclosure of `osc(H)`; and
`certify_H0` does `H(0)` itself, and hence `e^{H(0)}`, the constant
Theorem 10.1 converges to. That last one writes
`H(0) = Hhat(0) + 2 Re sum_{m>=1} Hhat(m)`: the zero mode comes from its
closed form with Arb's enclosures of the Euler-Mascheroni and first
Stieltjes constants, the first twelve nonzero modes from the same
`hhat` used throughout the file, and the rest from twice the majorant
the oscillation bound already uses, which at `M = 12` is below `1e-50`.

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
