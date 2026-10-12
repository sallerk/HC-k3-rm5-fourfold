"""Structure of the 22-dimensional first-order solution space (C, Phibar) of chain 1-7 (f4 notation) at
the explicit member: dimension of the part with C = 0, and whether it equals {(0, r (x) v^dual) : r in R_7} for the
first-order persistent v in R_5 (Phibar(x) = <v, x> r, pairing R_5 x R_7 -> R_12).
Usage: python -I f6_split22.py prime
"""
import os
import sys

import flint
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ivhs_ts as I
import b10_ivhs_new as B10
import f1_v_annihilator as F1
import f3_second_order as F3

S1 = [{(1, 0, 0): 1}, {(0, 1, 0): 1}, {(0, 0, 1): 1}]


def main():
    p = int(sys.argv[1])
    Fn = F1.family()
    params = list(F1.P[1:]) + list(F1.Bs) + [F1.eps]
    sub = {F1.eps: 5, F1.Bs[0]: 1, F1.Bs[1]: 0, F1.Bs[2]: 0, F1.Bs[3]: 1,
           F1.P[1]: 0, F1.P[2]: 0, F1.P[3]: 0, F1.P[4]: 0, F1.P[5]: 1, F1.P[6]: 1}
    Fd = B10.to_dict(sp.expand(4 * Fn.subs(sub)), F1.XS, p)
    R = I.JacRing(Fd, p)
    tans = [B10.to_dict(sp.expand(4 * sp.diff(Fn, v).subs(sub)), F1.XS, p) for v in params]
    n7 = R.dim(7)
    N = 9 + n7 * n7
    rows = []
    for Ft in tans:
        M = [R.reduce(7, I.pmul(P, Ft, p)) for P in S1]
        for j in range(3):
            for a in range(n7):
                row = [0] * N
                for b in range(n7):
                    row[9 + a * n7 + b] = M[j][b]
                for i in range(3):
                    row[i * 3 + j] = (row[i * 3 + j] - M[i][a]) % p
                rows.append(row)
    full = F1.kernel(rows, N, p)
    c0 = F1.kernel(rows + [[1 if k == t else 0 for k in range(N)] for t in range(9)], N, p)
    cproj = flint.nmod_mat([v[:9] for v in full], p).rank()
    print('first-order solutions: %d; with C = 0: %d; rank of projection to C: %d' % (len(full), len(c0), cproj))
    # the r (x) v^dual solutions
    Ps = F1.persistent(R, tans, 5)
    assert len(Ps) == 1
    v = Ps[0]
    # functional x -> <v, x> on R_7 in std coordinates
    fun = [F3.pair(R, v, [1 if k == t else 0 for k in range(n7)]) for t in range(n7)]
    span = []
    for r in range(n7):
        vec = [0] * N
        for a in range(n7):
            if a == r:
                for b in range(n7):
                    vec[9 + a * n7 + b] = fun[b]
        span.append(vec)
    both = flint.nmod_mat(span + c0, p).rank()
    print('dim span{r (x) v^dual} = %d; dim (that + C=0 solutions) = %d' % (flint.nmod_mat(span, p).rank(), both))


if __name__ == '__main__':
    main()
