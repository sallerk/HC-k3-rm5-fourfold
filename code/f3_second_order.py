"""Second-order persistence of the first-order persistent vectors v in piece j = 3 of
Y_F = {w^6 + F = 0} (Hodge numbers (0, 18, 3); pairs with piece j = 1, Hodge numbers (3, 18, 0)).
A rational class lambda that stays Hodge along an arc F(s) has <lambda, omega_P(s)> == 0 for the holomorphic forms
omega_P(s) = Res(w P Omega / G(s)), P in S_1 (piece j = 1), G = w^6 + F. Since F^2 = 0 in piece j = 3 and
F^1 = everything in piece j = 1, all pairings are algebraic:
  order 1: <v, -P Fdot> (v in (R_F)_5, P Fdot in (R_F)_7, product in (R_F)_12);
  order 2: <v, (1/2) sum_i d_i a_i - P Fddot>, where 2 P Fdot^2 = sum_i a_i d_i F (degree 13; R_13 = 0), by Griffiths'
           reduction (sum B_i d_i G) Omega / G^{k+1} == (1/k) (sum d_i B_i) Omega / G^k.
For an arc p(s) = p0 + s pdot in the parameter space: Fdot = sum pdot_k dF/dp_k, Fddot = sum pdot_k pdot_l
d^2F/dp_k dp_l. The second-order persistent space is P^(2) = {v in P_5 : <v, c2(P, pdot)> = 0 for all P, pdot}
(linear in v); a Hodge class of the very general member lies in P^(2).
Controls: (i) GL_3-orbit arcs F(exp(s M) x): c2 == 0 in R_7 for every P (checks the factor 1/2 and the reduction);
(ii) vGS family: P_5 = 4, P^(2) must contain the class of the hexaflex lines (expect >= 1).
Usage: python -I f3_second_order.py prime seed
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
import f1_v_annihilator as F1

XS = F1.XS


def padd(Pd, Qd, p, c=1):
    out = dict(Pd)
    for e, v in Qd.items():
        out[e] = (out.get(e, 0) + c * v) % p
    return {e: v for e, v in out.items() if v}


def vec_to_poly(R, d, vec):
    return {R.mons[d][j]: c for c, j in zip(vec, R.std[d]) if c}


def solve_jac(Fd, A, p):
    """a_0, a_1, a_2 in S_{deg A - 5} with sum a_i d_i F = A (exists when R_{deg A} = 0)."""
    dA = sum(next(iter(A)))
    parts = [I.pderiv(Fd, i, p) for i in range(3)]
    mons_out = I.monomials(dA)
    idx = {e: k for k, e in enumerate(mons_out)}
    mons_in = I.monomials(dA - 5)
    cols = []
    for i in range(3):
        for m in mons_in:
            col = [0] * len(mons_out)
            for e, c in I.pmul(parts[i], {m: 1}, p).items():
                col[idx[e]] = c
            cols.append(col)
    nr, nc = len(mons_out), len(cols)
    rhs = [0] * nr
    for e, c in A.items():
        rhs[idx[e]] = c % p
    aug = flint.nmod_mat([[cols[k][r] for k in range(nc)] + [rhs[r]] for r in range(nr)], p)
    Rm, rk = aug.rref()
    ent = [int(x) for x in Rm.entries()]
    w = nc + 1
    sol = [0] * nc
    for r in range(rk):
        row = ent[r * w:(r + 1) * w]
        piv = next(j for j, v in enumerate(row) if v)
        assert piv < nc, 'not in the Jacobian ideal'
        sol[piv] = row[nc]
    a = []
    for i in range(3):
        a.append({m: sol[i * len(mons_in) + k] for k, m in enumerate(mons_in) if sol[i * len(mons_in) + k]})
    return a


def c2_class(R, Fd, P, Fdot, Fddot, p):
    """(1/2) sum_i d_i a_i - P Fddot in (R_F)_7 (std coordinates), 2 P Fdot^2 = sum a_i d_i F."""
    A = I.pmul(I.pmul(P, Fdot, p), Fdot, p)
    A = {e: 2 * c % p for e, c in A.items()}
    a = solve_jac(Fd, A, p)
    div = {}
    for i in range(3):
        div = padd(div, I.pderiv(a[i], i, p), p)
    half = pow(2, p - 2, p)
    tot = padd({e: c * half % p for e, c in div.items()}, I.pmul(P, Fddot, p), p, -1)
    return R.reduce(7, tot)


def pair(R, v, c):
    """<v, c> = coordinate of v * c in (R_F)_12 (one-dimensional)."""
    prod = I.pmul(vec_to_poly(R, 5, v), vec_to_poly(R, 7, c), R.p)
    r = R.reduce(12, prod)
    assert len(r) == 1
    return r[0]


def second_order_space(R, Fd, Ps, arcs, p):
    """P^(2): v in span(Ps) with <v, c2(P, arc)> = 0 for all P in S_1 and arcs (Fdot, Fddot)."""
    S1 = [{(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}]
    rows = []
    for Fdot, Fddot in arcs:
        for P in S1:
            c = c2_class(R, Fd, P, Fdot, Fddot, p)
            rows.append([pair(R, v, c) for v in Ps])
    ker = F1.kernel(rows, len(Ps), p)
    return len(ker), rows


def family_arcs(Fn, params, sub, p, rng, n):
    grads = [sp.diff(Fn, v) for v in params]
    hess = [[sp.diff(g, v) for v in params] for g in grads]
    G = [B10.to_dict(sp.expand(4 * g.subs(sub)), XS, p) for g in grads]
    H = [[B10.to_dict(sp.expand(4 * h.subs(sub)), XS, p) for h in row] for row in hess]
    arcs = []
    for _ in range(n):
        pd = [rng.randrange(p) for _ in params]
        Fdot, Fddot = {}, {}
        for k in range(len(params)):
            Fdot = padd(Fdot, G[k], p, pd[k])
            for l in range(len(params)):
                Fddot = padd(Fddot, H[k][l], p, pd[k] * pd[l] % p)
        arcs.append((Fdot, Fddot))
    return G, arcs


def gl3_arcs(Fn_point, p, rng, n):
    """Arcs F(exp(s M) x): Fdot = sum (Mx)_i d_i F, Fddot = sum (Mx)_i (Mx)_j d_i d_j F + sum (M^2 x)_i d_i F."""
    x0, x1, x2 = XS
    arcs = []
    for _ in range(n):
        M = sp.Matrix(3, 3, lambda i, j: rng.randrange(1, 50))
        X = sp.Matrix([x0, x1, x2])
        mx, m2x = M * X, M * M * X
        d1 = [sp.diff(Fn_point, v) for v in XS]
        Fdot = sum(mx[i] * d1[i] for i in range(3))
        Fddot = sum(mx[i] * mx[j] * sp.diff(d1[i], XS[j]) for i in range(3) for j in range(3)) \
            + sum(m2x[i] * d1[i] for i in range(3))
        arcs.append((B10.to_dict(sp.expand(Fdot), XS, p), B10.to_dict(sp.expand(Fddot), XS, p)))
    return arcs


def main():
    p, seed = int(sys.argv[1]), int(sys.argv[2])
    rng = random.Random(48300 + seed)
    Fn = F1.family()
    params = list(F1.P[1:]) + list(F1.Bs) + [F1.eps]
    sub_ex = {F1.eps: 5, F1.Bs[0]: 1, F1.Bs[1]: 0, F1.Bs[2]: 0, F1.Bs[3]: 1,
              F1.P[1]: 0, F1.P[2]: 0, F1.P[3]: 0, F1.P[4]: 0, F1.P[5]: 1, F1.P[6]: 1}
    sub_rnd = {v: sp.Rational(rng.randrange(1, 10 ** 4), rng.randrange(1, 50)) for v in params}
    for name, sub in (('ex', sub_ex), ('rnd', sub_rnd)):
        Fpt = sp.expand(4 * Fn.subs(sub))
        Fd = B10.to_dict(Fpt, XS, p)
        R = I.JacRing(Fd, p)
        G, arcs = family_arcs(Fn, params, sub, p, rng, 6)
        Ps = F1.persistent(R, G, 5)
        # control (i): GL_3 arcs, all of R_5 must pass
        allv = [[1 if i == j else 0 for i in range(R.dim(5))] for j in range(R.dim(5))]
        k_gl, _ = second_order_space(R, Fd, allv, gl3_arcs(Fpt, p, rng, 2), p)
        k2, rows = second_order_space(R, Fd, Ps, arcs, p)
        print('[%s] P_5 = %d | second-order persistent dim = %d | control GL3 arcs: %d of %d pass'
              % (name, len(Ps), k2, k_gl, R.dim(5)))
        print('   pairings <v, c2(P, arc)> (first rows):', [r[0] for r in rows[:9]])
    # control (ii): vGS
    Fv, pv = F1.vgs_family()
    subv = {v: rng.randrange(1, 100) for v in pv}
    Fpt = sp.expand(Fv.subs(subv))
    Fd = B10.to_dict(Fpt, XS, p)
    R = I.JacRing(Fd, p)
    grads = [sp.diff(Fv, v) for v in pv]
    G = [B10.to_dict(sp.expand(g.subs(subv)), XS, p) for g in grads]
    Ps = F1.persistent(R, G, 5)
    arcs = []
    hess = [[sp.diff(g, v) for v in pv] for g in grads]
    H = [[B10.to_dict(sp.expand(h.subs(subv)), XS, p) for h in row] for row in hess]
    for _ in range(6):
        pd = [rng.randrange(p) for _ in pv]
        Fdot, Fddot = {}, {}
        for k in range(len(pv)):
            Fdot = padd(Fdot, G[k], p, pd[k])
            for l in range(len(pv)):
                Fddot = padd(Fddot, H[k][l], p, pd[k] * pd[l] % p)
        arcs.append((Fdot, Fddot))
    k2, _ = second_order_space(R, Fd, Ps, arcs, p)
    print('[vgs] P_5 = %d | second-order persistent dim = %d (expect >= 1: hexaflex lines)' % (len(Ps), k2))


if __name__ == '__main__':
    main()
