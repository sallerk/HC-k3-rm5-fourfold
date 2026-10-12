"""First-order IVHS chain bounds for X_{F,F} at the explicit member (Lemma 2.1(4) of the K3 note: eps = 5,
nu = x0^6 + x0 x1^5 + x1^6, b = x0^3 + x1^3; we use 4F), with the 11 parameter derivatives as tangents
(ivhs_ts.kernel_dims, as in b10_ivhs_new.py, which uses a seeded random member). Expect {0-6-12: 122, 1-7: 22, 2-8: 1,
3-9: 1, 4-10: 1, 5-11: 22}.
Usage: python -I f5_chain_bounds_ex.py prime
"""
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ivhs_ts as I
import b10_ivhs_new as B10
import f1_v_annihilator as F1


def main():
    p = int(sys.argv[1])
    Fn = F1.family()
    params = list(F1.P[1:]) + list(F1.Bs) + [F1.eps]
    sub = {F1.eps: 5, F1.Bs[0]: 1, F1.Bs[1]: 0, F1.Bs[2]: 0, F1.Bs[3]: 1,
           F1.P[1]: 0, F1.P[2]: 0, F1.P[3]: 0, F1.P[4]: 0, F1.P[5]: 1, F1.P[6]: 1}
    Fd = B10.to_dict(sp.expand(4 * Fn.subs(sub)), F1.XS, p)
    R = I.JacRing(Fd, p)
    print('dim R_d, d = 0..13:', [R.dim(d) for d in range(14)])
    tans = [B10.to_dict(sp.expand(4 * sp.diff(Fn, v).subs(sub)), F1.XS, p) for v in params]
    kd = I.kernel_dims(R, R, [(t, t) for t in tans], p)
    print('p = %d, explicit member, chain kernels:' % p, {k: v['kernel'] for k, v in kd.items()})
    # comparison: a random integral sextic with all 28 monomial directions (the full family of plane sextics);
    # expect 1 per chain (six planes), so dim Hdg^4 = 7 for a very general sextic
    import random
    rng = random.Random(48500)
    coeffs = [rng.randrange(-9, 10) for _ in I.monomials(6)]
    Fg = {e: c % p for e, c in zip(I.monomials(6), coeffs)}
    print('random sextic: integer coefficients', coeffs, 'on the monomials x0^a x1^b x2^c, (a, b, c) =',
          I.monomials(6))
    Rg = I.JacRing(Fg, p)
    assert Rg.dim(13) == 0
    kd = I.kernel_dims(Rg, Rg, [({e: 1}, {e: 1}) for e in I.monomials(6)], p)
    print('p = %d, random sextic, all 28 directions:' % p, {k: v['kernel'] for k, v in kd.items()})


if __name__ == '__main__':
    main()
