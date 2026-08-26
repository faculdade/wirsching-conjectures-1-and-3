# Section 9: the Berg-Kruppel identity

Backs Proposition 9.1 and, through it, Corollary 10.3.

| Script | Backs |
|---|---|
| `derive_P_bergkruppel_identity.py` | the exact `P`-Berg-Kruppel identity and the closed-form constant `C_P` |

Run: `python3 derive_P_bergkruppel_identity.py`

Section 9 of the paper also records a sign discrepancy in the printed
formula for `f''` on p. 179 of Berg and Kruppel (1998), where direct
differentiation gives `-2*beta` and the printed line carries `+2*beta`.
Their own Proposition 9.1 is unaffected, the term being of lower order
there. The exact expression used in this paper is not, and the paper
carries the correction explicitly.
