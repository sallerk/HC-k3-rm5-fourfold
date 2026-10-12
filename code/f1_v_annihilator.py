"""The first-order persistent vector v in (R_F)_5 (piece j = 3 of Y_F = {w^6 + F = 0}, the zeta_3
part, pulled back from Z_F = {u^3 + F = 0}) at an explicit member of the 8-dim family, and its annihilator in R_F.
P_5 = {v in (R_F)_5 : A v = 0 in (R_F)_11 for every tangent A (the 11 parameter derivatives)}.
For each persistent v: Ann(v)_d = kernel of (R_F)_d -> (R_F)_{5+d}, g -> g v, for d = 1..6 (dims and, in low degree,
the polynomials); dim T_HL(v) = dim {A in S_6 : A v = 0 in R_11}.
Members: 'ex' = the member of Lemma 2.1(4) of the K3 note (eps = 5, nu = x0^6 + x0 x1^5 + x1^6, b = x0^3 + x1^3); 'rnd' = random rational
member (seeded); control 'vgs' = eps = b0 = 0 member (hexaflex line x1 = 0 known).
Usage: python -I f1_v_annihilator.py prime seed
"""
import os
import random
import sys

import flint
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ivhs_ts as I
import b10_ivhs_new as B10

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
    return sp.expand(q.as_expr() + eps / 4 * x2 ** 6)


def vgs_family():
    D5 = lambda z, b: z ** 5 - 5 * b * z ** 3 + 5 * b ** 2 * z
    A = sp.symbols('a0:3')
    pb = sum(P[k] * x0 ** (6 - k) * x1 ** k for k in range(7))
    a = sum(A[k] * x0 ** (2 - k) * x1 ** k for k in range(3))
    return sp.expand(pb + x1 * D5(x2, a)), list(P) + list(A)


def mult(R, g, vec, a):
    """g (dict, degree d) times the element of (R_F)_a with std coordinates vec -> std coordinates in degree a + d."""
    p = R.p
    d = sum(next(iter(g)))
    acc = {}
    for c, j in zip(vec, R.std[a]):
        if c:
            m = R.mons[a][j]
            for e, cg in g.items():
                k = (e[0] + m[0], e[1] + m[1], e[2] + m[2])
                acc[k] = (acc.get(k, 0) + c * cg) % p
    return R.reduce(a + d, acc) if acc else [0] * R.dim(a + d)


def kernel(rows, ncols, p):
    if not rows:
        return [[1 if i == j else 0 for i in range(ncols)] for j in range(ncols)]
    M = flint.nmod_mat(rows, p)
    K, nk = M.nullspace()
    ent = [int(x) for x in K.entries()]
    return [[ent[i * ncols + k] for i in range(ncols)] for k in range(nk)]


def persistent(R, tans, a=5):
    p = R.p
    rows = []
    for A in tans:
        rows += R.mult_matrix(A, a)
    return kernel(rows, R.dim(a), p)


def poly_str(R, d, vec):
    terms = []
    for c, j in zip(vec, R.std[d]):
        if c:
            e = R.mons[d][j]
            c = c if c <= R.p // 2 else c - R.p
            terms.append('%d*x0^%d*x1^%d*x2^%d' % (c, *e))
    return ' + '.join(terms)


def analyse(name, R, tans, p):
    Ps = persistent(R, tans, 5)
    print('[%s] dim P_5 = %d' % (name, len(Ps)))
    for t, v in enumerate(Ps):
        dims = {}
        for d in range(1, 7):
            n = R.dim(d)
            # columns: images of the std basis of (R_F)_d times v
            cols = []
            for j in R.std[d]:
                cols.append(mult(R, {R.mons[d][j]: 1}, v, 5))
            rows = [[cols[k][i] for k in range(n)] for i in range(R.dim(5 + d))]
            ker = kernel(rows, n, p) if rows else [[1 if i == j else 0 for i in range(n)] for j in range(n)]
            dims[d] = len(ker)
            if d <= 3 and ker:
                for kv in ker[:4]:
                    print('   Ann(v%d)_%d element: %s' % (t, d, poly_str(R, d, kv)))
        # T_HL(v): A in S_6 (all monomials, not reduced) with A v = 0 in R_11
        mons6 = I.monomials(6)
        cols = [mult(R, {m: 1}, v, 5) for m in mons6]
        rows = [[cols[k][i] for k in range(len(mons6))] for i in range(R.dim(11))]
        rk = flint.nmod_mat(rows, p).rank()
        print('   v%d: dim Ann_d (d = 1..6) = %s | expected for generic v in R_5: dim R_d - dim R_{5+d} if > 0 |'
              ' dim T_HL(v) in S_6 = %d (codim %d)' % (t, dims, len(mons6) - rk, rk))


def main():
    p, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(48000 + seed)
    Fn = family()
    params = list(P[1:]) + list(Bs) + [eps]
    print('dims R_d (d = 0..12) at the explicit member:')
    sub_ex = {eps: 5, Bs[0]: 1, Bs[1]: 0, Bs[2]: 0, Bs[3]: 1, P[1]: 0, P[2]: 0, P[3]: 0, P[4]: 0, P[5]: 1, P[6]: 1}
    sub_rnd = {v: sp.Rational(rng.randrange(1, 10 ** 4), rng.randrange(1, 50)) for v in params}
    for name, sub in (('ex', sub_ex), ('rnd', sub_rnd)):
        F = B10.to_dict(sp.expand(4 * Fn.subs(sub)), XS, p)
        tans = [B10.to_dict(sp.expand(4 * sp.diff(Fn, v).subs(sub)), XS, p) for v in params]
        R = I.JacRing(F, p)
        if name == 'ex':
            print([R.dim(d) for d in range(13)])
        analyse(name, R, tans, p)
    Fv, pv = vgs_family()
    subv = {v: rng.randrange(1, 100) for v in pv}
    F = B10.to_dict(Fv.subs(subv), XS, p)
    tans = [B10.to_dict(sp.diff(Fv, v).subs(subv), XS, p) for v in pv]
    analyse('vgs', I.JacRing(F, p), tans, p)


if __name__ == '__main__':
    main()
