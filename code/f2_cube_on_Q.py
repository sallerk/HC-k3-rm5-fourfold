"""Does F restrict to a cube on the cubic curve E_Q = {Q = 0}, Q = (eps/2) x2^3 - b(x0, x1)?
(f1: Ann(v)_3 = <Q>.) Test 1: R(x0, x1) = Res_{x2}(F, Q), a binary form of degree 18; a cube if F|E_Q = h^3.
Test 2: solve F == c h^3 mod Q for a quadric h (exactly, symbolic parameters).
Usage: python -I f2_cube_on_Q.py
"""
import sympy as sp

x0, x1, x2 = XS = sp.symbols('x0 x1 x2')
P = sp.symbols('p0:7')
Bs = sp.symbols('b0:4')
eps = sp.Symbol('eps')


def family():
    D5 = lambda z, b: z ** 5 - 5 * b * z ** 3 + 5 * b ** 2 * z
    pb = sum(P[k] * x0 ** (6 - k) * x1 ** k for k in range(7))
    Bc = sum(Bs[k] * x0 ** (3 - k) * x1 ** k for k in range(4))
    num = sp.expand((pb * (x1 - eps * x2) + x1 ** 2 * D5(x2, Bc / x1)) * x1 ** 4)
    num = sp.expand(num.subs(P[0], 5 * Bs[0] ** 2 / eps) * eps)
    q, r = sp.div(sp.Poly(num, *XS), sp.Poly(x1 ** 5 * eps, *XS))
    assert r.is_zero
    return sp.expand(q.as_expr() + eps / 4 * x2 ** 6), pb, Bc


def main():
    F, pb, Bc = family()
    sub = {eps: 5, Bs[0]: 1, Bs[1]: 0, Bs[2]: 0, Bs[3]: 1, P[1]: 0, P[2]: 0, P[3]: 0, P[4]: 0, P[5]: 1, P[6]: 1}
    Fx = sp.expand(F.subs(sub))
    Q = sp.expand((eps / 2 * x2 ** 3 - Bc).subs(sub))
    R = sp.resultant(Fx, Q, x2)
    print('explicit member: Res_x2(F, Q) =', sp.factor(R))
    # symbolic: F mod Q as a polynomial of degree <= 2 in x2
    Qs = sp.expand(eps / 2 * x2 ** 3 - Bc)
    rem = sp.rem(sp.Poly(F, x2), sp.Poly(Qs, x2))
    print('F mod Q (coefficients in x2):')
    for k, c in enumerate(reversed(rem.all_coeffs())):
        print('  x2^%d:' % k, sp.factor(c))


if __name__ == '__main__':
    main()
