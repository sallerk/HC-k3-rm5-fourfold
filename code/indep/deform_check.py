"""deform_check.py -- first- and second-order (C, Phi) systems for a family of plane sextics.

Independent implementation (see README.txt of the repository for the method).

Usage:
    python -I deform_check.py c0   [narcs] [seed]   # the 11-parameter family at the base point c0
    python -I deform_check.py even [narcs] [seed]   # control: all sextics even in x2 (16 parameters)
Defaults: narcs = 6, seed = 20261011.  Output is printed and saved in the current folder as
deform_check_<family>.txt.

Unknowns: C in M_3(Q) (C(P_j) = sum_i C_ij P_i on S_1), Phi in End(R_7), and per arc D in M_3(Q).
First order:  Phi([P_j T_k]) - sum_i C_ij [P_i T_k] = 0 in R_7     (all j, k).
Second order: Phi(c2(P_j)) - sum_i C_ij c2(P_i) + 2 sum_i D_ij [P_i Fdot] = 0 in R_7  (all j, per arc),
  c2(P) = [(1/2) div(a) - P Fddot],  2 P Fdot^2 = sum_i a_i dF/dx_i.
All ranks exact over Q (flint.fmpq_mat); ranks mod p = 1000003 are printed as a cross-check only.
"""
import os
import random
import sys
import time

import sympy as sp
from flint import fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jacring import (JacobianRing, PRIME, linform, monomials, nullspace_q, padd, pdeg,  # noqa: E402
                     pdiff, pmul, pscale, pstr, rank_mod_p, to_fmpq_mat, var)

HERE = os.path.dirname(os.path.abspath(__file__))
X = sp.symbols("x0 x1 x2")


class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")

    def __call__(self, *args):
        s = " ".join(str(a) for a in args)
        print(s, flush=True)
        self.f.write(s + "\n")
        self.f.flush()


def sym_to_poly(expr):
    P = sp.Poly(sp.expand(expr), *X, domain="QQ")
    out = {}
    for m, c in P.terms():
        c = sp.Rational(c)
        out[tuple(int(e) for e in m)] = fmpq(int(c.p), int(c.q))
    return out


# ---------------------------------------------------------------- families
def family_c0():
    x0, x1, x2 = X
    nu = sp.symbols("nu1:7")
    b = sp.symbols("b0:4")
    eps = sp.Symbol("eps")
    params = list(nu) + list(b) + [eps]
    nu0 = 5 * b[0] ** 2 / eps
    NU = nu0 * x0 ** 6 + sum(nu[k - 1] * x0 ** (6 - k) * x1 ** k for k in range(1, 7))
    B = sum(b[k] * x0 ** (3 - k) * x1 ** k for k in range(4))
    NUp = sum(nu[k - 1] * x0 ** (6 - k) * x1 ** (k - 1) for k in range(1, 7))
    Bp = sum(b[k] * x0 ** (3 - k) * x1 ** (k - 1) for k in range(1, 4))
    F = (NU - eps * NUp * x2 - 5 * B * x2 ** 3 + 10 * b[0] * x0 ** 3 * Bp * x2
         + 5 * x1 * Bp ** 2 * x2 + x1 * x2 ** 5 + sp.Rational(1, 4) * eps * x2 ** 6)
    c0 = {eps: 5, nu[0]: 0, nu[1]: 0, nu[2]: 0, nu[3]: 0, nu[4]: 1, nu[5]: 1,
          b[0]: 1, b[1]: 0, b[2]: 0, b[3]: 1}
    F0 = sym_to_poly(F.subs(c0))
    T = [sym_to_poly(sp.diff(F, p).subs(c0)) for p in params]
    H = [[sym_to_poly(sp.diff(F, p, q).subs(c0)) for q in params] for p in params]
    expected = sym_to_poly(4 * x0 ** 6 + 40 * x0 ** 3 * x1 ** 2 * x2 - 20 * x0 ** 3 * x2 ** 3
                           + 4 * x0 * x1 ** 5 - 20 * x0 * x1 ** 4 * x2 + 4 * x1 ** 6
                           - 20 * x1 ** 3 * x2 ** 3 + 4 * x1 * x2 ** 5 + 5 * x2 ** 6)
    names = [str(p) for p in params]
    return dict(F=F0, T=T, H=H, names=names, check4F=(pscale(F0, fmpq(4)) == expected))


def family_even(seed):
    mons = [m for m in monomials(6) if m[2] % 2 == 0]
    assert len(mons) == 16
    rng = random.Random(seed + 7)
    tries = 0
    while True:
        tries += 1
        coeffs = [rng.randint(-5, 5) for _ in mons]
        F0 = {m: fmpq(c) for m, c in zip(mons, coeffs) if c != 0}
        if pdeg(F0) == 6 and JacobianRing(F0).dim(13) == 0:
            break
    T = [{m: fmpq(1)} for m in mons]
    H = None  # linear family: Fddot = 0
    names = ["x0^%d x1^%d x2^%d" % m for m in mons]
    return dict(F=F0, T=T, H=H, names=names, coeffs=coeffs, tries=tries)


def arc_polys(fam, d):
    Fdot = {}
    for dk, Tk in zip(d, fam["T"]):
        Fdot = padd(Fdot, Tk, dk)
    Fddot = {}
    if fam["H"] is not None:
        for k, dk in enumerate(d):
            for l, dl in enumerate(d):
                if dk * dl:
                    Fddot = padd(Fddot, fam["H"][k][l], dk * dl)
    return Fdot, Fddot


# ---------------------------------------------------------------- systems
def build_rows(R, fam, arcs_data, nR):
    """Return (first-order rows, list of second-order row blocks), rows as {col: fmpq}."""
    P = [var(i) for i in range(3)]
    nT = len(fam["T"])
    v = [[R.reduce(pmul(P[j], fam["T"][k]), 7) for k in range(nT)] for j in range(3)]
    PHI = 9

    def phi_col(a, b):
        return PHI + a * nR + b

    first = []
    for j in range(3):
        for k in range(nT):
            for a in range(nR):
                row = {}
                for b in range(nR):
                    if v[j][k][b] != 0:
                        row[phi_col(a, b)] = v[j][k][b]
                for i in range(3):
                    if v[i][k][a] != 0:
                        row[i * 3 + j] = row.get(i * 3 + j, fmpq(0)) - v[i][k][a]
                first.append(row)
    second = []
    for t, (c2, w) in enumerate(arcs_data):
        base = PHI + nR * nR + 9 * t
        block = []
        for j in range(3):
            for a in range(nR):
                row = {}
                for b in range(nR):
                    if c2[j][b] != 0:
                        row[phi_col(a, b)] = c2[j][b]
                for i in range(3):
                    if c2[i][a] != 0:
                        row[i * 3 + j] = row.get(i * 3 + j, fmpq(0)) - c2[i][a]
                    if w[i][a] != 0:
                        row[base + i * 3 + j] = 2 * w[i][a]
                block.append(row)
        second.append(block)
    return first, second, v


def submatrix(rows, cols):
    pos = {c: k for k, c in enumerate(cols)}
    out = []
    for r in rows:
        out.append({pos[c]: val for c, val in r.items() if c in pos})
    return to_fmpq_mat(out, len(cols))


def analyse(rows, ncols, nR, narcs, log, label):
    t0 = time.perf_counter()
    allc = list(range(ncols))
    M = submatrix(rows, allc)
    r = M.rank()
    rp = rank_mod_p(M)
    null = ncols - r
    Dcols = list(range(9 + nR * nR, ncols))
    if Dcols:
        rD = submatrix(rows, Dcols).rank()
        nullD = len(Dcols) - rD
    else:
        nullD = 0
    noC = list(range(9, ncols))
    nullnoC = len(noC) - submatrix(rows, noC).rank()
    proj = null - nullD
    projC = null - nullnoC
    dt = time.perf_counter() - t0
    log(f"[{label}] arcs={narcs}: rows={len(rows)} cols={ncols} rank_Q={r} rank_mod_p={rp} "
        f"| dim solutions={null} | dim(D-only kernel)={nullD} | dim proj to (C,Phi)={proj} "
        f"| dim proj to C={projC} | {dt:.2f}s")
    return dict(rank=r, rank_p=rp, null=null, nullD=nullD, proj=proj, projC=projC)


def check_solution(rows, sol):
    """Max over rows of |row . sol| == 0 ?"""
    for r in rows:
        s = fmpq(0)
        for c, v in r.items():
            if c in sol:
                s += v * sol[c]
        if s != 0:
            return False
    return True


# ---------------------------------------------------------------- main
def main():
    family = sys.argv[1] if len(sys.argv) > 1 else "c0"
    narcs = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20261011
    suffix = "" if seed == 20261011 else f"_seed{seed}"
    log = Tee(f"deform_check_{family}{suffix}.txt")
    T0 = time.perf_counter()
    log(f"deform_check.py family={family} narcs={narcs} seed={seed}  (exact over Q; mod p={PRIME} ranks only as cross-check)")

    fam = family_c0() if family == "c0" else family_even(seed)
    log(f"family built in {time.perf_counter() - T0:.2f}s; parameters: {fam['names']}")
    if family == "c0":
        log(f"4F == 4x0^6+40x0^3x1^2x2-20x0^3x2^3+4x0x1^5-20x0x1^4x2+4x1^6-20x1^3x2^3+4x1x2^5+5x2^6 ? {fam['check4F']}")
    else:
        log(f"random even member (coefficients on the 16 even monomials, {fam['tries']} tries): {fam['coeffs']}")
    log(f"F = {pstr(fam['F'])}")

    R = JacobianRing(fam["F"])
    t = time.perf_counter()
    dims = [R.dim(d) for d in range(14)]
    log(f"(a) dim R_d, d=0..13: {dims}   [{time.perf_counter() - t:.2f}s]")
    expected = [1, 3, 6, 10, 15, 18, 19, 18, 15, 10, 6, 3, 1, 0]
    log(f"    equals complete-intersection (5,5,5) Hilbert function: {dims == expected}; R_13 = 0: {dims[13] == 0}")
    if dims[13] != 0:
        log("F is singular -- stop")
        return
    nR = dims[7]
    log(f"    standard monomials of R_7 (grevlex): {R.std_monomials(7)}")

    # ---- controls (i): GL_3 arcs
    rng_gl = random.Random(seed + 1)
    P = [var(i) for i in range(3)]
    for trial in range(2):
        Mx = [[rng_gl.randint(-5, 5) for _ in range(3)] for _ in range(3)]
        MX = [linform(Mx[i]) for i in range(3)]  # (Mx)_i
        M2 = [[sum(Mx[i][k] * Mx[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
        M2X = [linform(M2[i]) for i in range(3)]
        Fdot, Fddot = {}, {}
        for i in range(3):
            Fdot = padd(Fdot, pmul(MX[i], R.dF[i]))
            Fddot = padd(Fddot, pmul(M2X[i], R.dF[i]))
            for j in range(3):
                Fddot = padd(Fddot, pmul(pmul(MX[i], MX[j]), pdiff(R.dF[i], j)))
        rl = [rng_gl.randint(-9, 9) for _ in range(3)]
        tests = P + [linform(rl)]
        zero = [all(c == 0 for c in R.c2(Q, Fdot, Fddot)) for Q in tests]
        # sign/factor sanity: the same arc with Fddot scaled by 2 should NOT give zero
        bad = [all(c == 0 for c in R.c2(Q, Fdot, pscale(Fddot, fmpq(2)))) for Q in tests]
        log(f"(d)(i) GL3 arc M={Mx}: c2(P)=0 for P in x0,x1,x2,{rl}: {zero}; "
            f"(with Fddot doubled, zero? {bad} -- expected False)")

    # ---- arcs
    rng = random.Random(seed)
    nP = len(fam["T"])
    arcs = [[rng.randint(-20, 20) for _ in range(nP)] for _ in range(narcs)]
    arcs_data = []
    t = time.perf_counter()
    for d in arcs:
        Fdot, Fddot = arc_polys(fam, d)
        c2 = [R.c2(P[j], Fdot, Fddot) for j in range(3)]
        w = [R.reduce(pmul(P[i], Fdot), 7) for i in range(3)]
        arcs_data.append((c2, w))
        log(f"arc d = {d}: c2 nonzero? {[any(c != 0 for c in x) for x in c2]}")
    log(f"arcs computed in {time.perf_counter() - t:.2f}s")

    # well-definedness check: change the lift of 2 P Fdot^2 by random Koszul syzygies
    rk = random.Random(seed + 2)
    Fdot, Fddot = arc_polys(fam, arcs[0])
    same = []
    for j in range(3):
        a = R.lift(pscale(pmul(P[j], pmul(Fdot, Fdot)), fmpq(2)))
        for (p, q) in [(0, 1), (0, 2), (1, 2)]:
            h = {m: fmpq(rk.randint(-9, 9)) for m in monomials(3)}
            a[p] = padd(a[p], pmul(h, R.dF[q]))       # a + h*(F_q e_p - F_p e_q)
            a[q] = padd(a[q], pmul(h, R.dF[p]), -1)
        same.append(R.c2_from_lift(P[j], a, Fddot) == arcs_data[0][0][j])
    log(f"Koszul check (arc 1, lift changed by random Koszul syzygies): c2(P_j) unchanged: {same}")

    first, second, v = build_rows(R, fam, arcs_data, nR)
    span = to_fmpq_mat([vv for j in range(3) for vv in v[j]], nR).rank()
    log(f"rank of {{[P_j T_k]}} in R_7: {span} (of {nR})")
    n1 = 9 + nR * nR
    log(f"(b) FIRST ORDER: unknowns {n1} (C: 9, Phi: {nR}x{nR}), equations {len(first)}")
    res1 = analyse(first, n1, nR, 0, log, "first order")
    # explicit basis of the first-order C-projection (rref of the C-parts of a kernel basis)
    ker = nullspace_q(submatrix(first, list(range(n1))))
    Cs = to_fmpq_mat([k[:9] for k in ker], 9).rref()
    log(f"    first-order C-space (rref basis, C as row-major 3x3 [C_00,C_01,C_02,C_10,...]; dim {Cs[1]}):")
    for i in range(Cs[1]):
        log("      " + str([Cs[0][i, c] for c in range(9)]))

    # sanity (iii)
    def ident(narc):
        sol = {i * 3 + i: fmpq(1) for i in range(3)}
        for a in range(nR):
            sol[9 + a * nR + a] = fmpq(1)
        return sol
    ok = check_solution(first + [r for b in second for r in b], ident(narcs))
    log(f"(d)(iii) (C,Phi,D) = (I, I, 0) solves first + all second-order equations: {ok}")
    if family == "even":
        stdm = R.std_monomials(7)
        sol = {0: fmpq(1), 4: fmpq(1), 8: fmpq(-1)}
        for a, m in enumerate(stdm):
            sol[9 + a * nR + a] = fmpq((-1) ** m[2])
        ok2 = check_solution(first + [r for b in second for r in b], sol)
        log(f"(d)(ii) sigma: x2 -> -x2, (C,Phi,D) = (diag(1,1,-1), diag((-1)^e2 on std monomials), 0) "
            f"solves everything: {ok2}")

    log(f"(c) SECOND ORDER (first-order equations + second-order equations of the first n arcs):")
    results = {}
    for n in sorted(set([1, 2, 4, 6, narcs])):
        if n > narcs:
            continue
        rows = first + [r for b in second[:n] for r in b]
        results[n] = analyse(rows, n1 + 9 * n, nR, n, log, "second order")
    log(f"total time {time.perf_counter() - T0:.2f}s")


if __name__ == "__main__":
    main()
