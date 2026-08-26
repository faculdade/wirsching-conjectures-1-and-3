# Section 4: Wirsching's Conjecture 1

Backs Theorem 4.3, the implication `(*2) => (*1)`, its two lemmas, and
Corollary 4.4.

| Script | Backs |
|---|---|
| `cancellation_check.py` | the generating-function cancellation `P_l(z) Q_l(z) = (1-z)^{-(l+1)}` and the closed form for the averaged Elka functions, in exact integer and rational arithmetic |
| `check_generating_identity.py` | the same identity by an independent route, plus a finite subexponential check on the partition counts |
| `partition_bound_check.py` | the elementary bound `p_l(m) <= (m+1)^(log_3(m/2)+2)` of equation (5) |
| `verify_lemma_constants.py` | the explicit constants of Lemma 4.2 (`rho_0 = 3/5`, `rho = 7/10`, `A` finite), the tail bound against the tail actually summed, and Lemma 4.1's sequence construction |

Run each with `python3 <name>.py`. Each prints its own assertions and
exits non-zero if one fails.

## What is proof and what is instantiation

The proof of Theorem 4.3 is finite and self-contained; none of these
scripts is part of it. `cancellation_check.py` and
`check_generating_identity.py` exercise the algebraic identity the proof
turns on. `verify_lemma_constants.py` instantiates the displayed
inequalities of Lemma 4.2 at concrete `ell`, `delta` and `tau`, which is
what Rule 11c of the project's own working rules asks of any displayed
estimate.

One honest detail `verify_lemma_constants.py` prints rather than hides:
at `delta = 5` the binomial ratio is still above `3/5` at `ell = 50`. The
lemma says "for all large `ell`", and the script shows the ratio
decreasing monotonically in `ell` and passing under `3/5`, rather than
asserting a threshold it has not checked.
