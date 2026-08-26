# Section 8: the envelope lemma

| Script | Backs |
|---|---|
| `verify_envelope_lemma.py` | the true-versus-smooth saddle comparison of Lemma 8.1, and coarse rigorous bounds on `H'` and `H''` |

Run: `python3 verify_envelope_lemma.py`

## What it does and does not check

**It does** solve for the true and smooth saddles independently and
print `B_0 |w* - w_0|` at four values of `tau`, against Lemma 8.1's
`|w* - w_0| <= 0.00652 / B_0`. Newton's step here uses
`F' = (M - V) - e^{w-tau}`; an earlier version used `-V - e^{w-tau}`,
which is wrong because `M' = -L'' = M - V`, and it left about half the
printed digits incorrect while still converging.

**It does not** compute the five terms `T_1` to `T_5` of the envelope
estimate. An earlier version of this README said it did. It does not.

**The derivative bounds it certifies are the coarse ones**, `|H'| < 0.007`
and `|H''| < 0.04`, not equation (10)'s `0.0011977...` and `0.0068518...`.
That distinction is load bearing and the script now prints it: Lemma 8.1
needs `2 e eta_0 <= B_0 = 5` with `eta_0 = 1/c + sup|H''|`, and

```
certified 0.00685 : eta_0 = 0.91708, 2 e eta_0 = 4.986  <= 5   ok
coarse    0.04    : eta_0 = 0.95024, 2 e eta_0 = 5.166  >  5   fails
```

So the coarse bounds do not reproduce the lemma. Equation (10)'s
certified values do, and they are produced in
`../section-06-periodic-correction/`. The script asserts both halves of
that comparison, so if the certified bound ever loosened enough for the
lemma to fail, or if the coarse bound ever tightened enough to make the
comparison stale, it stops.
