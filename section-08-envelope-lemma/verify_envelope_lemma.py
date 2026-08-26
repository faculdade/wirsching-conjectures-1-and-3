"""Independent numerical audit of the H-006 envelope calculation.

This is deliberately standalone: it does not import the project's K, saddle, or
Fourier routines.  It evaluates the defining Laplace-product summands directly,
constructs H from the telescoping identity, and compares the two saddles.
"""
import mpmath as mp
from flint import arb, ctx

mp.mp.dps = 90
c = mp.log(3)
a = mp.mpf('0.5') - mp.log(2) / c
alpha = 2 * mp.pi / c


def qpoly(w):
    return -w * w / (2 * c) + a * w


def f(x):
    # log((1-exp(-2x))/(2x)); expm1 avoids cancellation at small x.
    return mp.log(-mp.expm1(-2 * x)) - mp.log(2 * x)


def mfun(x):
    # -x f'(x).  Terms below the cutoff have a geometric O(x) tail.
    return 1 - 2 * x / mp.expm1(2 * x)


def vfun(x):
    # x^2 f''(x), evaluated in the stable sinh form.
    return 1 - (x / mp.sinh(x)) ** 2


def laplace_log_and_m_v(w):
    """Return L(w), M=-L'(w), and V=L''(w)-L'(w), from the product."""
    s = mp.exp(w)
    ell = mp.mpf('0')
    mm = mp.mpf('0')
    vv = mp.mpf('0')
    x = s / 3
    for _ in range(700):
        if x < mp.mpf('1e-75'):
            break
        ell += f(x)
        mm += mfun(x)
        vv += vfun(x)
        x /= 3
    return ell, mm, vv


def tiny_d(w):
    return mp.log1p(-mp.exp(-2 * mp.exp(w)))


def H_from_telescoping(w):
    """H = (L-Q)+sum_{k>=0} d(w+kc), the exact periodic completion."""
    ell, _, _ = laplace_log_and_m_v(w)
    correction = mp.mpf('0')
    u = w
    for _ in range(20):
        term = tiny_d(u)
        correction += term
        if abs(term) < mp.mpf('1e-80'):
            break
        u += c
    return ell - qpoly(w) + correction


def hhat(n):
    om = alpha * n
    return -(mp.power(2, 1j * om) / c) * mp.gamma(-1j * om) * mp.zeta(1 - 1j * om)


def fourier_derivative_bound(order):
    # Bound (2) in H-006 plus log(alpha*m+1)<=m log(alpha+1), and
    # m^(order+1/2)<=m^3 for order 1,2.
    q = mp.exp(-mp.pi * alpha / 2)
    A = mp.sqrt(3 * mp.pi / alpha) / c
    C = 1 + 1 / alpha + mp.mpf('0.5') + mp.sqrt(1 + alpha ** -2) / 2
    D = mp.log(alpha + 1) + C
    sum_m3_qm = q * (1 + 4 * q + q * q) / (1 - q) ** 4
    return 2 * alpha ** order * A * D * sum_m3_qm


def certify_coarse_derivative_bounds():
    """Outward-rounded certification of the two coarse bounds used below."""
    ctx.dps = 90
    cc = arb(3).log()
    aa = 2 * arb.pi() / cc
    qq = (-arb.pi() * aa / 2).exp()
    AA = (3 * arb.pi() / aa).sqrt() / cc
    CC = arb(1) + 1 / aa + arb(1) / 2 + (1 + 1 / aa**2).sqrt() / 2
    DD = (aa + 1).log() + CC
    sm3 = qq * (1 + 4 * qq + qq**2) / (1 - qq)**4
    b1 = 2 * aa * AA * DD * sm3
    b2 = 2 * aa**2 * AA * DD * sm3
    # These are the coarse bounds. Lemma 8.1's own step needs equation
    # (10)'s certified values: with sup|H''| <= 0.04 the lemma's
    # 2 e eta_0 <= B_0 fails, so certifying only the coarse bounds does
    # not reproduce the lemma this script is said to back.
    assert b1.upper() < arb('0.007')
    assert b2.upper() < arb('0.04')
    # equation (10), certified in ../section-06-periodic-correction/.
    H1_CERT = mp.mpf('0.0011977472315550332')
    H2_CERT = mp.mpf('0.0068518962896650951')
    print("  coarse here: |H'| < 0.007, |H''| < 0.04")
    print("  equation (10), certified in section 6: |H'| <= %s, |H''| <= %s"
          % (mp.nstr(H1_CERT, 12), mp.nstr(H2_CERT, 12)))
    c_log3 = mp.log(3)
    eta_cert = 1 / c_log3 + H2_CERT
    eta_coarse = 1 / c_log3 + mp.mpf('0.04')
    print("  Lemma 8.1 needs 2 e eta_0 <= B_0 = 5, with eta_0 = 1/c + sup|H''|.")
    print("    with the certified 0.00685: eta_0 = %s, 2 e eta_0 = %s"
          % (mp.nstr(eta_cert, 8), mp.nstr(2 * mp.e * eta_cert, 8)))
    print("    with the coarse    0.04   : eta_0 = %s, 2 e eta_0 = %s"
          % (mp.nstr(eta_coarse, 8), mp.nstr(2 * mp.e * eta_coarse, 8)))
    if not 2 * mp.e * eta_cert <= 5:
        raise AssertionError("2 e eta_0 exceeds B_0 even with the certified bound")
    if 2 * mp.e * eta_coarse <= 5:
        raise AssertionError(
            "the coarse bound now suffices, so this comparison is stale")
    print("  So the coarse bounds this script certifies do NOT reproduce")
    print("  Lemma 8.1; equation (10)'s certified values do, and they come")
    print("  from section 6. This script backs the saddle comparison, not (10).")
    return b1, b2


print('Fourier/tail bounds obtained independently from the coefficient formula')
print('  coarse sup |H\'|  <=', mp.nstr(fourier_derivative_bound(1), 20))
print('  coarse sup |H\'\'| <=', mp.nstr(fourier_derivative_bound(2), 20))
b1_ball, b2_ball = certify_coarse_derivative_bounds()
print('  Arb certifies the convenient rigorous bounds |H\'|<0.007, |H\'\'|<0.04')

for theta in [mp.mpf('0'), mp.log(mp.mpf('1.5')), mp.mpf('0.7')]:
    # Evaluation at theta+12c makes the remaining telescoping correction invisible.
    print('H(', mp.nstr(theta, 8), ')=', mp.nstr(H_from_telescoping(theta + 12*c), 35))


def smooth_saddle(tau):
    # Solve the large (W_{-1}) branch B exp(-w)=exp(-tau), w=c(B+a).
    w = tau + mp.log(tau / c)
    for _ in range(30):
        B = w / c - a
        F = w - mp.log(B) - tau
        w -= F / (1 - 1 / (c * B))
    return w


def true_saddle(tau, guess):
    # Stationarity is M(exp(w))=exp(w-tau); M is summed from the defining product.
    w = guess
    for _ in range(30):
        _, M, V = laplace_log_and_m_v(w)
        F = M - mp.exp(w - tau)
        # F' = M' - exp(w-tau). With M = -L' and L'' = L' + V we get
        # M' = -L'' = M - V, NOT -V. The wrong sign here made Newton
        # contract at about 1/2 per step instead of quadratically, and
        # left roughly half the printed digits wrong.
        w -= F / ((M - V) - mp.exp(w - tau))
    return w


print('\nTrue versus smooth saddle (all quantities from this file only)')
print(' tau          B0*delta_w       B0*(g*-g0full)    B0*normalizer-diff')
for tau in [mp.mpf(20), mp.mpf(40), mp.mpf(80), mp.mpf(160)]:
    w0 = smooth_saddle(tau)
    ws = true_saddle(tau, w0)
    B0 = w0 / c - a
    L0, _, _ = laplace_log_and_m_v(w0)
    Ls, _, Vs = laplace_log_and_m_v(ws)
    g0full = L0 + mp.exp(w0 - tau)
    gs = Ls + mp.exp(ws - tau)
    normdiff = (ws - mp.log(2*mp.pi*Vs)/2) - (w0 - mp.log(2*mp.pi*(B0-1/c))/2)
    print(mp.nstr(tau, 4).rjust(4),
          mp.nstr(B0*(ws-w0), 16).rjust(22),
          mp.nstr(B0*(gs-g0full), 16).rjust(24),
          mp.nstr(B0*normdiff, 16).rjust(25))
