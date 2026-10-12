"""The system of Proposition 4.3 (chain 1-7) at c0 with the two arcs printed in Proposition 4.4 of the note,
d1 = (11, 2, 11, 1, -3, -17, -18, 2, 10, -18, -20), d2 = (2, -7, 18, 13, -12, 10, -19, -20, -6, -5, -11) in the
coordinates (nu1..nu6, b0..b3, eps), with the first implementation (f4_chain_second_order.bound) modulo p.
Prints the first-order bound, and the projection to (C, Phibar) and the full solution space with d1 and with d1, d2.
Expect 22, then 4 (d1) and 1 (d1, d2), full space 4 and 1.
Usage: python -I f7_printed_arcs.py prime
"""
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ivhs_ts as I
import b10_ivhs_new as B10
import f1_v_annihilator as F1
import f3_second_order as F3
import f4_chain_second_order as F4

D1 = [11, 2, 11, 1, -3, -17, -18, 2, 10, -18, -20]
D2 = [2, -7, 18, 13, -12, 10, -19, -20, -6, -5, -11]


def main():
    p = int(sys.argv[1])
    Fn = F1.family()
    params = list(F1.P[1:]) + list(F1.Bs) + [F1.eps]
    sub = {F1.eps: 5, F1.Bs[0]: 1, F1.Bs[1]: 0, F1.Bs[2]: 0, F1.Bs[3]: 1,
           F1.P[1]: 0, F1.P[2]: 0, F1.P[3]: 0, F1.P[4]: 0, F1.P[5]: 1, F1.P[6]: 1}
    Fd = B10.to_dict(sp.expand(4 * Fn.subs(sub)), F1.XS, p)
    R = I.JacRing(Fd, p)
    grads = [sp.diff(Fn, v) for v in params]
    G = [B10.to_dict(sp.expand(4 * g.subs(sub)), F1.XS, p) for g in grads]
    H = [[B10.to_dict(sp.expand(4 * sp.diff(g, v).subs(sub)), F1.XS, p) for v in params] for g in grads]
    arcs = []
    for pd in (D1, D2):
        Fdot, Fddot = {}, {}
        for k in range(len(params)):
            Fdot = F3.padd(Fdot, G[k], p, pd[k] % p)
            for l in range(len(params)):
                Fddot = F3.padd(Fddot, H[k][l], p, pd[k] * pd[l] % p)
        arcs.append((Fdot, Fddot))
    for n in (1, 2):
        d1, d2 = F4.bound(R, Fd, G, arcs[:n], p)
        print('p = %d, c0, printed arcs d1..d%d: first-order bound %d | second-order bound %d | full solution space %d'
              % (p, n, d1, d2, F4.bound.full2))


if __name__ == '__main__':
    main()
