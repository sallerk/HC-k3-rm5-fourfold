"""IVHS (infinitesimal Noether-Lefschetz) bound for X = {F(x) + G(y) = 0} in P^5.

Jacobian ring of F(x) + G(y) is R_F (x) R_G (Thom-Sebastiani). H^{2,2}_prim(X) = R_12 = sum_a (R_F)_a (x) (R_G)_{12-a}.
A tangent direction of the family is a pair (A, B) of ternary sextics (A in x, B in y; either may be 0), acting on
R_12 -> R_18 by multiplication with A(x) + B(y). K = {lambda in R_12 : (A(x) + B(y)) lambda = 0 in R_18 for all
tangents}. Writing lambda = sum_a L_a with L_a in (R_F)_a (x) (R_G)_{12-a} (a matrix r_a x r_{12-a}), the
component of the product in bidegree (a+6, 12-a) is M_A^{(a)} L_a + L_{a+6} (M_B^{(6-a)})^T (a = 0..6); all other
components vanish for degree reasons (R_F, R_G are zero above degree 12). The bidegree chains
{0,6,12}, {1,7}, {2,8}, {3,9}, {4,10}, {5,11} decouple; we report the kernel dimension of each.
dim K is an upper bound for dim Hdg^2_prim at a very general point of the family.
All arithmetic is exact over F_p (python-flint nmod_mat).

Cases: two-block (control, expect 0), diagonal X_{F,F} (lower bound 6), vGS Sec. 5.9 diagonal (lower bound 8).
Usage: python -I ivhs_ts.py OUT.json [prime] [seeds]
"""
import itertools
import json
import random
import sys

import flint

TOP = 12  # socle degree of R_F for a ternary sextic


def monomials(n):
    return [e for e in itertools.product(range(n + 1), repeat=3) if sum(e) == n]


def padd(P, Q, p, c=1):
    R = dict(P)
    for e, v in Q.items():
        R[e] = (R.get(e, 0) + c * v) % p
    return {e: v for e, v in R.items() if v}


def pmul(P, Q, p):
    R = {}
    for e1, c1 in P.items():
        for e2, c2 in Q.items():
            e = (e1[0] + e2[0], e1[1] + e2[1], e1[2] + e2[2])
            R[e] = (R.get(e, 0) + c1 * c2) % p
    return {e: v for e, v in R.items() if v}


def pderiv(P, i, p):
    R = {}
    for e, c in P.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            f = tuple(f)
            R[f] = (R.get(f, 0) + c * e[i]) % p
    return {e: v for e, v in R.items() if v}


class JacRing:
    """R_F = F_p[x0,x1,x2]/(dF) in degrees 0..TOP, with a standard-monomial basis per degree."""

    def __init__(self, F, p):
        self.p = p
        parts = [pderiv(F, i, p) for i in range(3)]
        self.mons, self.idx, self.rows, self.piv, self.std, self.stdpos = {}, {}, {}, {}, {}, {}
        for n in range(TOP + 7):
            mons = monomials(n)
            idx = {e: k for k, e in enumerate(mons)}
            gens = []
            if n >= 5:
                for g in parts:
                    for m in monomials(n - 5):
                        row = [0] * len(mons)
                        for e, c in g.items():
                            row[idx[(e[0] + m[0], e[1] + m[1], e[2] + m[2])]] = c
                        gens.append(row)
            if gens:
                Rm, rk = flint.nmod_mat(gens, p).rref()
                ent = [int(x) for x in Rm.entries()]
                nc = len(mons)
                rows = [ent[i * nc:(i + 1) * nc] for i in range(rk)]
            else:
                rows = []
            piv = [next(j for j, v in enumerate(r) if v) for r in rows]
            pivset = set(piv)
            std = [j for j in range(len(mons)) if j not in pivset]
            self.mons[n], self.idx[n], self.rows[n], self.piv[n], self.std[n] = mons, idx, rows, piv, std
            self.stdpos[n] = {j: k for k, j in enumerate(std)}

    def dim(self, n):
        return len(self.std[n]) if 0 <= n <= TOP + 6 else 0

    def reduce(self, n, P):
        """Coordinates of the degree-n polynomial P in the standard basis of (R_F)_n."""
        p = self.p
        v = [0] * len(self.mons[n])
        for e, c in P.items():
            v[self.idx[n][e]] = (v[self.idx[n][e]] + c) % p
        for r, j in zip(self.rows[n], self.piv[n]):
            c = v[j]
            if c:
                v = [(a - c * b) % p for a, b in zip(v, r)]
        return [v[j] for j in self.std[n]]

    def mult_matrix(self, A, a):
        """Matrix (r_{a+6} x r_a) of multiplication by the sextic A: (R_F)_a -> (R_F)_{a+6}."""
        cols = []
        for j in self.std[a]:
            m = self.mons[a][j]
            cols.append(self.reduce(a + 6, pmul(A, {m: 1}, self.p)))
        ra6 = self.dim(a + 6)
        return [[cols[k][i] for k in range(len(cols))] for i in range(ra6)]


CHAINS = [(0, 6, 12), (1, 7), (2, 8), (3, 9), (4, 10), (5, 11)]


def kernel_dims(RF, RG, tangents, p):
    """tangents: list of (A, B) with A a sextic in x (dict) or {} and B a sextic in y (dict) or {}."""
    out = {}
    for chain in CHAINS:
        off, nv = {}, 0
        for a in chain:
            off[a] = nv
            nv += RF.dim(a) * RG.dim(TOP - a)
        rows = []
        for A, B in tangents:
            for a in chain:
                if a > 6:
                    continue
                b = a + 6
                ra, rb = RF.dim(a), RF.dim(b)
                sa, sb = RG.dim(TOP - a), RG.dim(TOP - b)  # sb = r^G_{6-a}
                MA = RF.mult_matrix(A, a) if A else None      # rb x ra
                MB = RG.mult_matrix(B, TOP - b) if B else None  # sa x sb
                for i in range(rb):
                    for j in range(sa):
                        row = [0] * nv
                        if MA:
                            for k in range(ra):
                                c = MA[i][k]
                                if c:
                                    row[off[a] + k * sa + j] = (row[off[a] + k * sa + j] + c) % p
                        if MB:
                            for l in range(sb):
                                c = MB[j][l]
                                if c:
                                    row[off[b] + i * sb + l] = (row[off[b] + i * sb + l] + c) % p
                        if any(row):
                            rows.append(row)
        rk = flint.nmod_mat(rows, p).rank() if rows else 0
        out["-".join(map(str, chain))] = {"unknowns": nv, "kernel": nv - rk}
    return out


def rand_sextic(p, rng):
    return {e: rng.randrange(p) for e in monomials(6)}


def vgs(p, rng):
    """vGS 2310.05196 Sec. 5.9: F = pb(x0,x1) + x1 * D5(x2, a), D5(x, a) = x^5 - 5 a x^3 + 5 a^2 x.
    Returns F and its 10 parameter derivatives (7 coefficients of pb, 3 of a)."""
    pb = {(i, 6 - i, 0): rng.randrange(p) for i in range(7)}
    amons = [(2, 0, 0), (1, 1, 0), (0, 2, 0)]
    a = {m: rng.randrange(p) for m in amons}
    x1 = {(0, 1, 0): 1}
    x2 = {(0, 0, 1): 1}
    x23 = {(0, 0, 3): 1}
    d5 = padd(padd({(0, 0, 5): 1}, pmul(a, x23, p), p, -5), pmul(pmul(a, a, p), x2, p), p, 5)
    F = padd(pb, pmul(x1, d5, p), p)
    tangents = [{(i, 6 - i, 0): 1} for i in range(7)]
    for m in amons:
        # d/da_m of D5 = -5 m x2^3 + 10 a m x2
        dd5 = padd({e: (-5 * c) % p for e, c in pmul({m: 1}, x23, p).items()},
                   {e: (10 * c) % p for e, c in pmul(pmul({m: 1}, a, p), x2, p).items()}, p)
        tangents.append(pmul(x1, dd5, p))
    return F, tangents


def tritangent(p, rng):
    """Positive control: F = c(x0,x2)^2 + x1 * H (the line x1 = 0 is tritangent to C_F), so Pic(S_F) contains
    H5-type classes a, b (a + b = h) and NS_prim (x) NS_prim gives one extra Hodge class on X_{F,F} (lower bound 7).
    Returns F and its parameter derivatives (2 c dc for dc in binary cubics, x1 dH for dH in ternary quintics)."""
    cmons = [(i, 0, 3 - i) for i in range(4)]
    c = {m: rng.randrange(p) for m in cmons}
    H = {e: rng.randrange(p) for e in monomials(5)}
    x1 = {(0, 1, 0): 1}
    F = padd(pmul(c, c, p), pmul(x1, H, p), p)
    tangents = [{e: (2 * v) % p for e, v in pmul(c, {m: 1}, p).items()} for m in cmons]
    tangents += [pmul(x1, {e: 1}, p) for e in monomials(5)]
    return F, tangents


def ab2_family(p, rng, parts):
    """F = A^2 + prod(parts) where A is a ternary cubic and parts are ternary forms of degrees summing to 6.
    Every factor curve {B = 0} meets C_F only where A = 0 (tangency), so its pullback to S_F splits.
    Returns F and all parameter derivatives (2 A dA, and dB_j * prod_{i != j} B_i)."""
    A = {e: rng.randrange(p) for e in monomials(3)}
    Bs = [{e: rng.randrange(p) for e in monomials(d)} for d in parts]

    def prod(lst):
        r = {(0, 0, 0): 1}
        for b in lst:
            r = pmul(r, b, p)
        return r

    F = padd(pmul(A, A, p), prod(Bs), p)
    tangents = [{e: (2 * v) % p for e, v in pmul(A, {m: 1}, p).items()} for m in monomials(3)]
    for j, d in enumerate(parts):
        rest = prod([b for i, b in enumerate(Bs) if i != j])
        tangents += [pmul(rest, {m: 1}, p) for m in monomials(d)]
    return F, tangents


def smooth(R):
    return [R.dim(n) for n in range(TOP + 2)] == [1, 3, 6, 10, 15, 18, 19, 18, 15, 10, 6, 3, 1, 0]


def main():
    outpath = sys.argv[1]
    p = int(sys.argv[2]) if len(sys.argv) > 2 else 32003
    seeds = int(sys.argv[3]) if len(sys.argv) > 3 else 3
    allsext = monomials(6)
    res = {"prime": p, "cases": []}
    for s in range(seeds):
        rng = random.Random(5000 + s)
        F, G = rand_sextic(p, rng), rand_sextic(p, rng)
        RF, RG = JacRing(F, p), JacRing(G, p)
        assert smooth(RF) and smooth(RG)
        two = [({m: 1}, {}) for m in allsext] + [({}, {m: 1}) for m in allsext]
        kd = kernel_dims(RF, RG, two, p)
        res["cases"].append({"case": "two_block", "seed": s, "chains": kd,
                             "K": sum(v["kernel"] for v in kd.values())})
        diag = [({m: 1}, {m: 1}) for m in allsext]
        kd = kernel_dims(RF, RF, diag, p)
        res["cases"].append({"case": "diagonal", "seed": s, "chains": kd,
                             "K": sum(v["kernel"] for v in kd.values())})
        F5, tans = vgs(p, random.Random(6000 + s))
        R5 = JacRing(F5, p)
        assert smooth(R5)
        kd = kernel_dims(R5, R5, [(t, t) for t in tans], p)
        res["cases"].append({"case": "vgs59_diagonal", "seed": s, "chains": kd,
                             "K": sum(v["kernel"] for v in kd.values())})
        Ft, tant = tritangent(p, random.Random(7000 + s))
        Rt = JacRing(Ft, p)
        assert smooth(Rt)
        kd = kernel_dims(Rt, Rt, [(t, t) for t in tant], p)
        res["cases"].append({"case": "tritangent_diagonal", "seed": s, "chains": kd,
                             "K": sum(v["kernel"] for v in kd.values())})
        # control: no tangent directions -> K = all of R_12 (1751)
        if s == 0:
            kd = kernel_dims(RF, RF, [], p)
            res["cases"].append({"case": "no_tangents_control", "seed": s, "chains": kd,
                                 "K": sum(v["kernel"] for v in kd.values())})
    for c in res["cases"]:
        print(c["case"], c["seed"], "K =", c["K"], {k: v["kernel"] for k, v in c["chains"].items()})
    with open(outpath, "w", newline="\n") as fh:
        fh.write(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
