# Section 10: Conjecture 3

Backs Theorem 10.1 numerically. The theorem does not depend on this.

| Script | Backs |
|---|---|
| `conjecture3_numerical_sweep.py` | `phi(z_l)/phi_0(z_l)` along Wirsching's own sequence `z_l = l * 3^-l`, `l = 5..50` |

Run: `python3 conjecture3_numerical_sweep.py`

Expect the ratio to dip to about `0.1827` near `l = 20` and then rise
slowly back toward the theorem's closed-form limit, `0.204988` in the
BARE normalization, that is without the leading constant
`(2 beta)^eps / sqrt(2 pi)`. That constant is `e^{C_P} = 0.3837844092...`,
computed in `../section-09-berg-kruppel-identity/`, and the limit in Berg
and Kruppel's own normalization, the one the paper's equation (2) uses,
is `e^{H(0)} = 0.5341220367...`, derived in
`../section-12-numerics/concordance.py`. Dividing the two ROUNDED values
`0.204988 / 0.3837844` gives `0.5341228`, whose last two digits are a
rounding artifact and not `e^{H(0)}`; an earlier version of this README
quoted that artifact.

Full numerical convergence needs `l` well beyond what a naive
root-finder handles reliably at this precision. This is illustrative
supporting evidence, not a proof.

`phi` here is evaluated by the real-variable saddlepoint approximation,
justified rigorously and uniformly by Theorem 7.5, not by oscillatory
Fourier inversion. For an evaluation of `phi` with zero truncation
error, see `../section-12-numerics/`.
