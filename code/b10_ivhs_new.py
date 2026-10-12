"""IVHS (first-order Noether-Lefschetz) bound for X_{F,F} = {F(x) + F(y) = 0} in P^5 with F in the new 8-dim RM
family, per Thom-Sebastiani chain (uses ivhs_ts.py).
The kernel on a chain is an upper bound for the primitive Hodge classes of X_{F,F} in that chain at a very general
member of the family (mod-p kernel >= char-0 kernel at the reduction of a rational point).
Cases: new family (11 parameter derivatives), vGS (control: expect 76 / 17 / 2 on {1,7} / {2,8} / {3,9}),
random 8-dim subfamily at the same F (control: 1 on non-K3 chains, 1 + 11^2 on the K3 chain).
Usage: python -I b10_ivhs_new.py OUT.json prime seed
"""
import json
import os
import random
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ivhs_ts as I


def to_dict(expr, xs, p):
    d = {}
    for mon, c in sp.Poly(sp.expand(expr), *xs).as_dict().items():
        c = sp.Rational(c)
        v = int(c.p % p) * pow(int(c.q % p), p - 2, p) % p
        if v:
            d[tuple(int(u) for u in mon)] = v
    return d


def main():
    out, p, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    rng = random.Random(11000 + seed)
    x0, x1, x2 = xs = sp.symbols('x0 x1 x2')
    P = sp.symbols('p0:7')
    Bs = sp.symbols('b0:4')
    eps = sp.Symbol('eps')
    D5 = lambda z, b: z ** 5 - 5 * b * z ** 3 + 5 * b ** 2 * z
    pb = sum(P[k] * x0 ** (6 - k) * x1 ** k for k in range(7))
    Bc = sum(Bs[k] * x0 ** (3 - k) * x1 ** k for k in range(4))
    num = sp.expand((pb * (x1 - eps * x2) + x1 ** 2 * D5(x2, Bc / x1)) * x1 ** 4)
    num = sp.expand(num.subs(P[0], 5 * Bs[0] ** 2 / eps) * eps)
    q, r = sp.div(sp.Poly(num, *xs), sp.Poly(x1 ** 5 * eps, *xs))
    assert r.is_zero
    Fn = sp.expand(q.as_expr() + eps / 4 * x2 ** 6)
    params = list(P[1:]) + list(Bs) + [eps]
    sub = {v: sp.Rational(rng.randrange(1, 10 ** 5), rng.randrange(1, 10 ** 2)) for v in params}
    F = to_dict(Fn.subs(sub), xs, p)
    tans = [to_dict(sp.diff(Fn, v).subs(sub), xs, p) for v in params]
    R = I.JacRing(F, p)
    res = {"prime": p, "seed": seed, "smooth": I.smooth(R), "cases": []}
    assert res["smooth"], 'sextic not smooth'
    kd = I.kernel_dims(R, R, [(t, t) for t in tans], p)
    res["cases"].append({"case": "new8", "chains": kd})
    # control: random 8-dim subfamily at the same F
    allsext = I.monomials(6)
    rt = [{e: rng.randrange(p) for e in allsext} for _ in range(8)]
    kd = I.kernel_dims(R, R, [(t, t) for t in rt], p)
    res["cases"].append({"case": "random8_same_F", "chains": kd})
    # control: vGS
    F5, t5 = I.vgs(p, random.Random(12000 + seed))
    R5 = I.JacRing(F5, p)
    kd = I.kernel_dims(R5, R5, [(t, t) for t in t5], p)
    res["cases"].append({"case": "vgs", "chains": kd})
    for c in res["cases"]:
        print(c["case"], {k: v["kernel"] for k, v in c["chains"].items()})
    with open(out, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
