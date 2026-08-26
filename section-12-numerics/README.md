# Section 12: the evaluation of phi, and the concordance

Backs the two numbers Section 12 puts beside Theorem 10.1.

| Script | Backs |
|---|---|
| `experiment_conjecture3.py` | the evaluation of `phi` through depth `l = 500`, with zero truncation error at every point tested, and the two extrapolation fits |
| `concordance.py` | the comparison of those fitted intercepts against the proved limit `e^{H(0)}` |

## What is exact and what is not

The evaluation of `phi` is exact in the sense that matters: moments are
exact rationals from the self-similarity `X = (2U + X)/3`, combined with
an antiderivative reduction, so the truncation error is exactly zero at
every point reported, not merely small. Working precision is 100 digits,
so the printed decimals are a floating evaluation of an exact quantity.

The extrapolation is not certified and the paper says so. Two model
forms fit the tested range comparably:

```
Lambda + b/sqrt(l)   : Lambda = -0.618860, c = 0.53856, max residual 2.14e-5
Lambda + b'/ln^2(l)  : Lambda = -0.599498, c = 0.5491, max residual 5.51e-5
```

`concordance.py` derives `H(0)` from Theorem 6.1's own decomposition,
cross-checks it against section 6's certified enclosure of
`H(0) - H(log 3/2)`, and puts the fits against the resulting
`e^{H(0)} = 0.5341220367...`, in the same normalization of `phi_0`. It
hardcodes neither the limit nor the movement. The two overshoot
by `0.83%` and `2.80%`, the interval they span excludes the proved
value, and the gap on the `Lambda` scale is `19.6` times the largest
movement any within-model check produced, which is `4.22e-4` from
adding a further `1/l` term. An earlier version of this README and of
the paper said "more than twenty times"; that came from rounding the
movement down to `4e-4` and the intercept to `-0.6189`, both in the
claim's favour, and it was false.

That shows the fitted intercepts are wrong over the range fitted. It
does not show either asymptotic expansion is wrong, and the script says
so: an approach `Lambda_true + b l^-1/2 + o(l^-1/2)` with a large
lower-order term produces exactly this.

## Running

`experiment_conjecture3.py` is the long one. `concordance.py` is
instant and depends on nothing but the two intercepts printed above.
