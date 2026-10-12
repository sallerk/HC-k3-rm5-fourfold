"""Second-order bound for chain 1-7 of X_{F,F} = {F(x) + F(y) = 0} (and, by x <-> y, chain 5-11).
Chain 1-7 = H_2 (x) H_4 with H_2 = piece j = 1 of Y_F = {w^6 + F = 0} (Hodge (3,18,0); A = F^2 H_2 = H^{2,0} spanned
by omega_P = Res(w P Omega / G), P in S_1) and H_4 = piece j = 3 (Hodge (0,18,3)), H_4 = H_2^dual via cup product.
A class h in H_2 (x) H_4 is a map phi: H_2 -> H_2, and h in F^2 iff phi(A) in A. So a Hodge class of the very general
member gives a flat phi with phi(A_s) in A_s for all s. Write phi omega_P(s) = sum_Q c_PQ(s) omega_Q(s). Its Gr^2
image is (c, phibar), c = c(0) in End(S_1), phibar in End(H_2 / A) = End((R_F)_7). Differentiating along an arc:
  order 1: phibar(P Fdot) = (c P) Fdot in R_7, all P in S_1, all tangents Fdot (the first-order test of b10_ivhs_new.py; expect 22);
  order 2: phibar(c2(P)) = c2(c P) - 2 (cdot P) Fdot in R_7 for some cdot in End(S_1) (one unknown cdot per arc),
           c2(P) = class of nabla^2 omega_P mod A (f3.c2_class).
Bound = dim of the projection to (c, phibar) of the solution space. Target 1 (identity = the planes y = zeta x).
Controls: first order must give 22 (new family) and 76 at vGS (as b10_ivhs_new.py); vGS at second order must keep identity + the
hexaflex-line class (expect >= 2).
Usage: python -I f4_chain_second_order.py prime seed narcs [small]   (small: arc directions with entries in [-20, 20])
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
import f3_second_order as F3

XS = F1.XS
S1 = [{(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}]
SMALL = len(sys.argv) > 4 and sys.argv[4] == "small"  # arc directions in [-20, 20]


def tangents_and_arcs(Fn, params, sub, p, rng, narcs, scale):
    grads = [sp.diff(Fn, v) for v in params]
    G = [B10.to_dict(sp.expand(scale * g.subs(sub)), XS, p) for g in grads]
    H = [[B10.to_dict(sp.expand(scale * sp.diff(g, v).subs(sub)), XS, p) for v in params] for g in grads]
    arcs = []
    for _ in range(narcs):
        pd = [rng.randrange(-20, 21) % p if SMALL else rng.randrange(p) for _ in params]
        Fdot, Fddot = {}, {}
        for k in range(len(params)):
            Fdot = F3.padd(Fdot, G[k], p, pd[k])
            for l in range(len(params)):
                Fddot = F3.padd(Fddot, H[k][l], p, pd[k] * pd[l] % p)
        arcs.append((Fdot, Fddot))
    return G, arcs


def bound(R, Fd, tans, arcs, p):
    n7 = R.dim(7)
    assert R.dim(1) == 3
    nc, nphi, K = 9, n7 * n7, len(arcs)
    N = nc + nphi + 9 * K
    ic = lambda i, j: i * 3 + j                  # c(P_j) = sum_i c_ij P_i
    iphi = lambda a, b: nc + a * n7 + b          # phibar(x)_a = sum_b phi_ab x_b
    icd = lambda k, i, j: nc + nphi + 9 * k + i * 3 + j
    rows1 = []
    for Ft in tans:
        M = [R.reduce(7, I.pmul(P, Ft, p)) for P in S1]
        for j in range(3):
            for a in range(n7):
                row = [0] * N
                for b in range(n7):
                    row[iphi(a, b)] = M[j][b]
                for i in range(3):
                    row[ic(i, j)] = (row[ic(i, j)] - M[i][a]) % p
                rows1.append(row)
    rows2 = []
    for k, (Fdot, Fddot) in enumerate(arcs):
        m = [R.reduce(7, I.pmul(P, Fdot, p)) for P in S1]
        c2 = [F3.c2_class(R, Fd, P, Fdot, Fddot, p) for P in S1]
        for j in range(3):
            for a in range(n7):
                row = [0] * N
                for b in range(n7):
                    row[iphi(a, b)] = c2[j][b]
                for i in range(3):
                    row[ic(i, j)] = (row[ic(i, j)] - c2[i][a]) % p
                    row[icd(k, i, j)] = (row[icd(k, i, j)] + 2 * m[i][a]) % p
                rows2.append(row)

    def proj_dim(rows):
        ker = F1.kernel(rows, N, p)
        if not ker:
            return 0, 0
        return len(ker), flint.nmod_mat([v[:nc + nphi] for v in ker], p).rank()

    _, d1 = proj_dim(rows1)
    full2, d2 = proj_dim(rows1 + rows2)
    bound.full2 = full2  # dim of the full solution space (c, phibar, cdot); upper semicontinuous in the point
    return d1, d2


def main():
    p, seed, narcs = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(48400 + seed)
    Fn = F1.family()
    params = list(F1.P[1:]) + list(F1.Bs) + [F1.eps]
    sub_ex = {F1.eps: 5, F1.Bs[0]: 1, F1.Bs[1]: 0, F1.Bs[2]: 0, F1.Bs[3]: 1,
              F1.P[1]: 0, F1.P[2]: 0, F1.P[3]: 0, F1.P[4]: 0, F1.P[5]: 1, F1.P[6]: 1}
    sub_rnd = {v: sp.Rational(rng.randrange(1, 10 ** 4), rng.randrange(1, 50)) for v in params}
    for name, sub in (('ex', sub_ex), ('rnd', sub_rnd)):
        Fd = B10.to_dict(sp.expand(4 * Fn.subs(sub)), XS, p)
        R = I.JacRing(Fd, p)
        tans, arcs = tangents_and_arcs(Fn, params, sub, p, rng, narcs, 4)
        d1, d2 = bound(R, Fd, tans, arcs, p)
        print('[%s] chain 1-7: first-order bound %d (b10: 22) | second-order bound %d (%d arcs; full solution space '
              '%d)' % (name, d1, d2, narcs, bound.full2))
        if name == 'ex':
            rnd8 = [{e: rng.randrange(p) for e in I.monomials(6)} for _ in range(8)]
            c1, _ = bound(R, Fd, rnd8, [], p)
            print('   control: 8 random sextic directions at the same F, first order %d (b10: 1)' % c1)
    Fv, pv = F1.vgs_family()
    subv = {v: rng.randrange(1, 100) for v in pv}
    Fd = B10.to_dict(sp.expand(Fv.subs(subv)), XS, p)
    R = I.JacRing(Fd, p)
    tans, arcs = tangents_and_arcs(Fv, pv, subv, p, rng, narcs, 1)
    d1, d2 = bound(R, Fd, tans, arcs, p)
    print('[vgs] chain 1-7: first-order bound %d (b10: 76) | second-order bound %d (expect >= 2)' % (d1, d2))
    # sharp positive control: all sextics even in x2 (16 coefficients); g = (x2 -> -x2) commutes with mu_6, so the
    # planes y = zeta g(x) give a second Hodge class in the chain: bound >= 2 (id, g)
    mons = [e for e in I.monomials(6) if e[2] % 2 == 0]
    cs = sp.symbols('c0:%d' % len(mons))
    Fe = sum(c * x0 ** e[0] * x1 ** e[1] * x2 ** e[2] for c, e in zip(cs, mons) for x0, x1, x2 in [XS])
    sube = {c: rng.randrange(1, 100) for c in cs}
    Fd = B10.to_dict(sp.expand(Fe.subs(sube)), XS, p)
    R = I.JacRing(Fd, p)
    tans, arcs = tangents_and_arcs(Fe, list(cs), sube, p, rng, narcs, 1)
    d1, d2 = bound(R, Fd, tans, arcs, p)
    print('[even] chain 1-7: first-order bound %d | second-order bound %d (full %d; expect >= 2: id, g)'
          % (d1, d2, bound.full2))


if __name__ == '__main__':
    main()
