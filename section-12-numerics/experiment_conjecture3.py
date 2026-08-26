#!/usr/bin/env python3
# NOTE ON RECORDS. Hypothesis identifiers such as H-013, the hypothesis
# tracker and the literature index refer to the framework repository,
# which is not public. Everything a reader needs to reproduce the paper's
# computations is in this repository.
"""
H-001 - teste computacional certificado da Conjectura 3 de Wirsching
(2003, "On positive predecessor density in 3n+1 dynamics", DCDS 9(3)).

Objeto: a densidade invariante phi (ponto fixo unico do operador de
medias W_3 f(x) = (3/2)*integral_{3x-2}^{3x} f(t)dt sobre [0,1]) e' o
analogo base-3 da funcao de Fabius (densidade de X = soma_j 2*U_j*3^-j,
U_j uniformes iid em [0,1]. A identidade e' verificada abaixo:
X =d (2U+X)/3 reproduz exatamente a equacao de W_3).

Conjectura 3 (a mais concreta da cadeia de 3 conjecturas do Wirsching,
reduzindo densidade positiva de predecessores 3n+1): existe c>0 tal que

    lim_{l->inf} phi(z_l)/phi_0(z_l) = c > 0

uniformemente para sequencias (z_l) na classe A_delta (janela CLT
|l - k_l| <= delta*sqrt(l)), com phi_0 a assintotica fechada de
Berg-Kruppel (1998, Proposicao 9.1, pp.178-179, item 133 do INDEX):

    phi_0(t) ~ (2*beta)^eps/sqrt(2*pi) * t^gamma * (-ln t)^delta
               * exp(-beta*ln^2(t/(-ln t)))

com (a=3, lambda=2/3, nossos parametros):
    alpha = 1/2 - ln(2)/ln(3)
    beta  = 1/(2*ln(3))
    delta = 1/2 + alpha - 2*beta*ln(2*beta)
    gamma = -2*beta - delta - 1/2
    eps   = 1/2 + alpha - beta*ln(2*beta)

Metodo (SEM iterar W_3): os momentos M_i = E[X^i] de
phi sao RACIONAIS EXATOS via a autossimilaridade acima. Para calcular
phi(x) numa cauda extrema x~l*3^-l, reduz-se via primitivas iteradas:

    phi(x) = (3/2)^(m+1) * 3^(-m(m+1)/2) * F_{m+1}(3^(m+1)*x)

com m = maior inteiro tal que 3^m*x<=2/3, e F_j(y) = integral_0^min(y,1)
((y-t)^(j-1)/(j-1)!)*phi(t)dt. Para y em [1,2] (caso A), F_j(y) e' uma
soma binomial exata dos momentos - erro ZERO. Para y em (2/3,1) (caso
B), usa a simetria phi(t)=phi(1-t) (Berg-Kruppel Prop. 4.1) + reescala
recursiva, com cota certificada no truncamento.

Trabalha SEMPRE em log (phi decai como 3^(-l^2/2), exponenciar
destroi a precisao) - usa mpmath com dps alto.

Reproduzir: python3 experiment_conjecture3.py
"""
import argparse
import math
import sys
import time
from fractions import Fraction as Fr

from mpmath import mp, mpf, log as mlog, pi as mpi

sys.set_int_max_str_digits(0)  # denominadores crescem ~N^2 digitos (N~300 -> ~10k digitos)
mp.dps = 100

LOG3 = mlog(3)
LOG2 = mlog(2)

ALPHA = mpf('0.5') - LOG2 / LOG3
BETA = 1 / (2 * LOG3)
DELTA = mpf('0.5') + ALPHA - 2 * BETA * mlog(2 * BETA)
GAMMA = -2 * BETA - DELTA - mpf('0.5')
EPS = mpf('0.5') + ALPHA - BETA * mlog(2 * BETA)


def log_phi0(x_frac):
    """ln phi_0(t) - assintotica de Berg-Kruppel (9.6), t=x_frac (Fraction) pequeno."""
    t = mpf(x_frac.numerator) / mpf(x_frac.denominator)
    lnt = mlog(t)
    ln_neglnt = mlog(-lnt)
    ln_inner = lnt - ln_neglnt  # ln(t/(-ln t))
    return (EPS * mlog(2 * BETA) - mpf('0.5') * mlog(2 * mpi)
            + GAMMA * lnt + DELTA * ln_neglnt - BETA * ln_inner ** 2)


def moments(N):
    """Momentos racionais exatos M_i=E[X^i], X=d(2U+X)/3, U~Unif[0,1] iid."""
    M = [Fr(1)]
    for i in range(1, N + 1):
        s = sum(Fr(math.comb(i, k)) * Fr(2 ** k, k + 1) * M[i - k] for k in range(1, i + 1))
        M.append(s / (3 ** i - 1))
    return M


def E_pow(M, y, n):
    """E[(y-X)^n], y racional, exato."""
    return sum(Fr(math.comb(n, i)) * (-1) ** i * y ** (n - i) * M[i] for i in range(n + 1))


def F(M, j, y, depth=80):
    """F_j(y) = int_0^min(y,1) (y-t)^(j-1)/(j-1)! phi(t) dt.
    Retorna (valor_racional, cota_de_erro_racional). y em (2/3, 2]."""
    assert Fr(2, 3) < y <= 2, f"y fora do range: {y}"
    main_num = E_pow(M, y, j - 1) / math.factorial(j - 1)
    if y >= 1:
        return main_num, Fr(0)
    if depth == 0:
        w = 1 - y
        bound = w ** (j - 1) / math.factorial(j - 1)
        return main_num, bound
    w = 1 - y
    pref = Fr(1)
    jj, ww = j, w
    while ww <= Fr(2, 3):
        pref *= Fr(3, 2) * Fr(1, 3 ** jj)
        jj += 1
        ww *= 3
    sub, suberr = F(M, jj, ww, depth - 1)
    val = main_num + (-1) ** j * pref * sub
    return val, pref * suberr


def log_phi_reduced(M, x):
    """ln phi(x) exato (mpmath), + cota de erro relativo, + (m,y) usados."""
    m = 0
    while 3 ** (m + 1) * x <= Fr(2, 3):
        m += 1
    y = 3 ** (m + 1) * x
    val, err = F(M, m + 1, y)
    logval = mlog(mpf(val.numerator) / mpf(val.denominator))
    logphi = (m + 1) * mlog(mpf(3) / 2) - mpf(m * (m + 1)) / 2 * LOG3 + logval
    rel = float(err / val) if val != 0 else float('inf')
    return logphi, rel, m, y, val


def sample_points(ell, us=(-2, -1, -0.5, 0, 0.5, 1, 2), shifted=True):
    """pontos k_l = l + round(u*sqrt(l)) variando u (janela CLT, classe A_delta~2).

    `shifted` selects the evaluation point, and this is NOT cosmetic.

    Wirsching's (7.5) and this project's paper 05 both evaluate at
    x_l^+ := x_l + 3^(-l-1)  (main.tex, sec:conj3, following
    Wirsching 2003 section 7), not at the bare x_l = k/3^l. The two give
    different deficit coefficients for l*(2/3 - L_l): about 0.802 at the
    shifted point and about 0.580 at the bare one.

    This function returned the BARE point until 2026-08-17, while the
    paper quoted the shifted 0.802 and this experiment's own README quoted
    the bare 0.580, so the repository numerically contradicted the paper
    on the same named quantity. Found by critique round 1 (both critics
    independently). Default is now the paper's point; pass
    --bare-point to recover the old behaviour for comparison.
    """
    pts = []
    for u in us:
        k = ell + round(u * math.sqrt(ell))
        if k < 1:
            continue
        # x_l^+ = k/3^l + 1/3^(l+1) = (3k+1)/3^(l+1), exact
        x = Fr(3 * k + 1, 3 ** (ell + 1)) if shifted else Fr(k, 3 ** ell)
        if not (0 < x <= Fr(2, 3)):
            continue
        pts.append((u, k, x))
    return pts


def ratio_L(M, m, y):
    """L_l = 3^(1-l) * phi(3x+)/phi(x+) = 2*3^(m-l)*F_m(y)/F_{m+1}(y), mesmo y.
    Aqui como funcao auxiliar generica: retorna F_m(y) e F_{m+1}(y) exatos."""
    val_m, _ = F(M, m, y) if m >= 1 else (Fr(1), Fr(0))
    val_m1, _ = F(M, m + 1, y)
    return val_m, val_m1


def validate(M):
    print("=== Validacao ===")
    assert M[1] == Fr(1, 2) and M[2] == Fr(7, 24), "M1/M2 errados"
    for n in (3, 5, 7, 9):
        c = sum(Fr(math.comb(n, k)) * (Fr(-1, 2)) ** (n - k) * M[k] for k in range(n + 1))
        assert c == 0, f"momento central {n} != 0: {c}"
    print("  M1=1/2, M2=7/24 OK; momentos centrais impares (3,5,7,9) = 0 OK")
    # phi(1/2) deve ser exatamente 3/2 (Wirsching/BK)
    logphi_half, rel, m, y, val = log_phi_reduced(M, Fr(1, 2))
    phi_half = float(mp.e ** logphi_half) if logphi_half < 700 else float('inf')
    print(f"  phi(1/2) calculado = {phi_half:.10f} (esperado 1.5), rel_err_cota={rel:.1e}")
    assert abs(phi_half - 1.5) < 1e-8, "phi(1/2) deveria ser exatamente 3/2"
    print("  Validacao OK\n")


def fit_models(series):
    """Least squares for L_inf under two competing model forms.

    series: list of (ell, ln(phi/phi0)) along one u.

    Returns {name: (L_inf, coeff, max_abs_residual)}. Two forms are fitted
    because they are the paper's dominant systematic, not a refinement:
    E-001's own README records L = -0.619 +- 0.001(statistical)
    +- 0.015(functional form), and H-001 records that C/ln^2(l)
    is "equally stable ... indistinguishable from C/sqrt(l) in this range;
    it is the dominant source of systematic uncertainty". Critique round 1
    (2026-08-17) found the paper's abstract quoting only the 4e-4
    sub-range stability, which understates the real spread by more than an
    order of magnitude. Printing both fits puts the honest interval in the
    committed output.
    """
    import numpy as np
    out = {}
    for name, basis in (("C/sqrt(ell)", lambda l: 1.0 / math.sqrt(l)),
                        ("C/ln^2(ell)", lambda l: 1.0 / math.log(l) ** 2)):
        A = np.array([[1.0, basis(l)] for l, _ in series])
        b = np.array([v for _, v in series])
        (L_inf, coeff), *_ = np.linalg.lstsq(A, b, rcond=None)
        resid = float(np.max(np.abs(A @ np.array([L_inf, coeff]) - b)))
        out[name] = (float(L_inf), float(coeff), resid)
    # The paper also quotes the shift from adding a 1/l term to the
    # sqrt model. Critique round 3 found that figure with no run behind
    # it, so it is computed here rather than asserted.
    A = np.array([[1.0, 1.0 / math.sqrt(l), 1.0 / l] for l, _ in series])
    b = np.array([v for _, v in series])
    (L3, c1, c2), *_ = np.linalg.lstsq(A, b, rcond=None)
    resid3 = float(np.max(np.abs(A @ np.array([L3, c1, c2]) - b)))
    out["C/sqrt+1/ell"] = (float(L3), float(c1), resid3)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bare-point", action="store_true",
                    help="evaluate at the bare x_l instead of the paper's "
                         "x_l^+ = x_l + 3^(-l-1); for comparison only")
    ap.add_argument("--max-ell", type=int, default=500,
                    help="largest level to compute (default 500, the range "
                         "main.tex's Empirical Result quotes)")
    args = ap.parse_args()

    ELL_LIST = [e for e in
                [5, 10, 20, 30, 50, 75, 100, 150, 200, 250, 300,
                 350, 400, 450, 500]
                if e <= args.max_ell]
    # ELL_LIST stopped at 300 until 2026-08-17, so the ell=350..500 range
    # that main.tex's abstract and the numerical test of (*4) and (*5) depend on rested on
    # an inline run this file never reproduced (recorded in the README's own
    # Rule 9a block at the time). Critique round 1 raised it as major.
    N_MAX = max(ELL_LIST) + 10

    t0 = time.time()
    print(f"Calculando momentos exatos ate M_{N_MAX}...")
    M = moments(N_MAX)
    print(f"  Momentos prontos em {time.time()-t0:.1f}s\n")

    validate(M)

    print("=== Conjecture 3 test: ln(phi/phi0), and L_l which bears on (*4) only ===")
    print(f"{'ell':>5} {'u':>6} {'k_l':>8} {'m':>5} {'ln phi':>12} {'ln phi0':>12} "
          f"{'ln r':>10} {'L_l=phi(3x)/phi(x)*3^(1-l)':>28} {'rel_err':>10} {'tempo':>8}")

    shifted = not args.bare_point
    print(f"evaluation point: {'x_l^+ = x_l + 3^-(l+1)' if shifted else 'bare x_l'}"
          f"  ({'paper' if shifted else 'legacy, for comparison'})\n")
    lnr_by_u = {}
    for ell in ELL_LIST:
        for u, k, x in sample_points(ell, shifted=shifted):
            t1 = time.time()
            logphi, rel, m, y, val = log_phi_reduced(M, x)
            lp0 = log_phi0(x)
            lnr = float(logphi - lp0)

            # L_l: phi(3x)/phi(x) usa o MESMO y, so muda m->m-1
            # phi(x) = (3/2)^(m+1) 3^(-m(m+1)/2) F_{m+1}(y)
            # phi(3x)=(3/2)^m 3^(-(m-1)m/2) F_m(y)  [profundidade m-1 para 3x]
            if m >= 1:
                val_m, _ = F(M, m, y)
                # phi(3x)/phi(x) = [ (3/2)^m 3^{-m(m-1)/2} F_m(y) ] / [ (3/2)^{m+1} 3^{-m(m+1)/2} F_{m+1}(y) ]
                #               = (2/3) * 3^{m} * F_m(y)/F_{m+1}(y)
                log_ratio_33 = mlog(mpf(2) / 3) + m * LOG3 + mlog(mpf(val_m.numerator) / mpf(val_m.denominator)) \
                               - mlog(mpf(val.numerator) / mpf(val.denominator))
                log_Ll = log_ratio_33 + (1 - ell) * LOG3
                Ll = float(mp.e ** log_Ll) if abs(log_Ll) < 700 else float('inf') * (1 if log_Ll > 0 else -1)
            else:
                Ll = float('nan')

            lnr_by_u.setdefault(u, []).append((ell, lnr, Ll))
            dt = time.time() - t1
            print(f"{ell:5d} {u:6.1f} {k:8d} {m:5d} {float(logphi):12.3f} {float(lp0):12.3f} "
                  f"{lnr:10.4f} {Ll:28.6f} {rel:10.2e} {dt:8.2f}s", flush=True)

    # --- the quantities main.tex quotes, computed here rather than inline ---
    print("\n=== Deficit ell*(2/3 - L_l), central sequence u=0 ===")
    for ell, _, Ll in lnr_by_u.get(0, []):
        if Ll == Ll:  # not nan
            print(f"  ell={ell:4d}   ell*(2/3 - L_l) = {ell * (2/3 - Ll):.4f}")

    if any(e >= 200 for e, _, _ in lnr_by_u.get(0, [])):
        print("\n=== Spread of ln(phi/phi0) across u at the largest level ===")
        top = max(e for e, _, _ in lnr_by_u.get(0, []))
        vals = [lnr for u, s in lnr_by_u.items() for e, lnr, _ in s if e == top]
        print(f"  ell={top}: min={min(vals):.6f} max={max(vals):.6f} "
              f"spread={max(vals)-min(vals):.3e}")

        print("\n=== Deficit range across u at the largest level ===")
        drange = [(u, top * (2/3 - Ll)) for u, s in lnr_by_u.items()
                  for e, _, Ll in s if e == top and Ll == Ll]
        for u, d in sorted(drange):
            print(f"  u={u:+.1f}: {d:+.3f}")

        print("\n=== L_inf fits along u=0, two model forms (Rule 11: report both) ===")
        series = [(e, lnr) for e, lnr, _ in lnr_by_u[0] if e >= 200]
        if len(series) >= 3:
            fits = fit_models(series)
            for name, (L_inf, coeff, resid) in fits.items():
                print(f"  {name:12s}  L_inf={L_inf:+.6f}  coeff={coeff:+.4f}  "
                      f"c=exp(L_inf)={math.exp(L_inf):.4f}  max|resid|={resid:.2e}")
            Ls = [v[0] for v in fits.values()]
            print(f"  MODEL SPREAD on L_inf: {abs(Ls[0]-Ls[1]):.4f}  "
                  f"-> c in [{min(math.exp(L) for L in Ls):.4f}, "
                  f"{max(math.exp(L) for L in Ls):.4f}]")
            print("  This spread, not the sub-range stability, is the dominant")
            print("  systematic on c. Both forms fit this range comparably.")
            sub = [(e, v) for e, v in series if e >= 350]
            if len(sub) >= 3:
                f2 = fit_models(sub)["C/sqrt(ell)"][0]
                print(f"  sub-range check (ell>=350, C/sqrt): L_inf={f2:+.6f}  "
                      f"shift={abs(f2 - fits['C/sqrt(ell)'][0]):.2e}")

    print(f"\n=== Concluido em {time.time()-t0:.1f}s ===")
    print("Conjecture 3 asserts a LIMIT. Bounded oscillation is not a\n  weaker form of it: if ln r oscillates without settling, the limit\n  does not exist and Conjecture 3 is false. What is measured here is\n  the behaviour of ln r over the range computed, nothing beyond it.")
    print("L_l (prediction via phi_0: L_l -> 2/3) is evidence for (*4),\n  NOT for (*5): it cannot see the periodic alternative that would\n  falsify Conjecture 3. Only ln(phi/phi_0) bears on (*5).")


if __name__ == "__main__":
    main()
