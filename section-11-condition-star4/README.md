# Section 11: condition (*4)

Backs Theorem 11.1, Remark 11.2 and Remark 11.4.

| Script | Backs |
|---|---|
| `verify_star4.py` | the three numerical claims of Section 11 |

Run: `python3 verify_star4.py`

What it checks:

1. **Remark 11.2's table.** Wirsching attributes the penultimate step of
   his own `(*5) => (*4)` calculation to his (7.7), an identity for
   `phi`, on a line that carries `phi_0`. What licenses the step is his
   (7.12), which holds only asymptotically. The script evaluates
   `(2/9) phi_0'(t) / phi_0(3t)` at `t = l * 3^-l` in closed form and
   shows it converging to 1 without ever equalling it, which is why the
   citation matters.
2. **`mu = 1/3`.** His (7.13) is introduced with "a somewhat lengthy
   calculation shows" and the calculation is not shown. The script
   reproduces the limit `2/3`, on the central sequence and off it, which
   is what gives `1 - mu = 2/3`. Since (7.5) reads `limsup <= 1 - mu < 1`,
   equality at `2/3` is admissible and `mu = 1/3` needs no weakening.
3. **Remark 11.4's residual.** The periodic factor cancels across a
   factor of 3 only asymptotically: at `l = 500` the residual saddle
   phase shift is `2.0e-3`, and the certified `sup|H'|` from
   `../section-06-periodic-correction/` bounds its contribution by
   `2.4e-6`. Small, and not zero.

Every derivative is taken in closed form. At `t ~ 500 * 3^-500` a
difference quotient is meaningless, and an early draft of this check
produced complex garbage that way before the closed form replaced it.
