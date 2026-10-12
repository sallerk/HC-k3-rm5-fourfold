"""jacring.py -- exact Jacobian-ring linear algebra for plane curves (independent implementation).

Polynomials in x0, x1, x2 are dicts {(e0, e1, e2): flint.fmpq}.  All arithmetic is exact over Q
(python-flint fmpq / fmpq_mat); a mod-p rank helper is provided only as a cross-check.

Normal forms in R = S/J, J = (dF/dx0, dF/dx1, dF/dx2):
  for each degree d the generators m*dF/dx_i (m a monomial of degree d-e, e = deg dF/dx_i)
  are written as rows over the monomials of S_d, listed in DESCENDING grevlex order
  (x0 > x1 > x2).  The reduced row echelon form of that matrix has pivots exactly at the
  grevlex leading monomials of J_d, so the non-pivot monomials are the grevlex standard
  monomials and "subtract the pivot rows" is the grevlex Groebner normal form, computed
  degree by degree with plain linear algebra.  Coordinates of a class in R_d are its
  coefficients on the standard monomials.

Lifting in the top degree: if R_D = 0, a form p of degree D is written p = sum a_i dF/dx_i
by choosing a maximal set of linearly independent generators m*dF/dx_i (pivot columns of the
rref of the transposed generator matrix) and inverting the resulting square matrix.

Usage: library only (imported by deform_check.py).
"""
import flint
from flint import fmpq, fmpq_mat, nmod_mat

PRIME = 1000003


# ---------------------------------------------------------------- polynomials
def monomials(d):
    """Exponent tuples of degree d in descending grevlex order (x0 > x1 > x2)."""
    mons = [(a, b, d - a - b) for a in range(d + 1) for b in range(d + 1 - a)]
    mons.sort(key=lambda e: (e[2], e[1]))
    return mons


def padd(p, q, c=1):
    """p + c*q."""
    r = dict(p)
    for m, v in q.items():
        r[m] = r.get(m, fmpq(0)) + c * v
    return {m: v for m, v in r.items() if v != 0}


def pscale(p, c):
    return {m: c * v for m, v in p.items() if c * v != 0}


def pmul(p, q):
    r = {}
    for m1, v1 in p.items():
        for m2, v2 in q.items():
            m = (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2])
            r[m] = r.get(m, fmpq(0)) + v1 * v2
    return {m: v for m, v in r.items() if v != 0}


def pdiff(p, i):
    r = {}
    for m, v in p.items():
        if m[i] > 0:
            mm = list(m)
            mm[i] -= 1
            r[tuple(mm)] = r.get(tuple(mm), fmpq(0)) + m[i] * v
    return {m: v for m, v in r.items() if v != 0}


def pdeg(p):
    degs = {sum(m) for m in p}
    if len(degs) > 1:
        raise ValueError("not homogeneous")
    return degs.pop() if degs else None


def var(i):
    e = [0, 0, 0]
    e[i] = 1
    return {tuple(e): fmpq(1)}


def linform(coeffs):
    """sum_i coeffs[i] * x_i."""
    return {tuple(1 if j == i else 0 for j in range(3)): fmpq(c) for i, c in enumerate(coeffs) if c != 0}


def pstr(p):
    names = ("x0", "x1", "x2")
    terms = []
    for m in monomials(pdeg(p)):
        if m in p:
            mon = "*".join(f"{names[i]}^{m[i]}" if m[i] > 1 else names[i] for i in range(3) if m[i])
            terms.append(f"({p[m]})*{mon}" if mon else f"({p[m]})")
    return " + ".join(terms) if terms else "0"


# ---------------------------------------------------------------- linear algebra helpers
def to_fmpq_mat(rows, ncols):
    """rows: list of dicts {col: value} or lists."""
    flat = []
    for r in rows:
        if isinstance(r, dict):
            line = [fmpq(0)] * ncols
            for c, v in r.items():
                line[c] = v
        else:
            line = list(r)
        flat.extend(line)
    return fmpq_mat(len(rows), ncols, flat)


def rank_mod_p(M, p=PRIME):
    """Rank of an fmpq_mat reduced mod p (denominators must be prime to p)."""
    m, n = M.nrows(), M.ncols()
    flat = []
    for i in range(m):
        for j in range(n):
            v = M[i, j]
            num, den = int(v.p), int(v.q)
            flat.append(num * pow(den, -1, p) % p)
    return nmod_mat(m, n, flat, p).rank()


def nullspace_q(M):
    """Basis (list of lists of fmpq) of the right kernel of an fmpq_mat, read off its rref."""
    A, r = M.rref()
    piv = pivots_of_rref(A, r)
    pset = set(piv)
    basis = []
    for f in range(M.ncols()):
        if f in pset:
            continue
        x = [fmpq(0)] * M.ncols()
        x[f] = fmpq(1)
        for i, pc in enumerate(piv):
            x[pc] = -A[i, f]
        basis.append(x)
    return basis


def pivots_of_rref(A, r):
    piv = []
    for i in range(r):
        for j in range(A.ncols()):
            if A[i, j] != 0:
                piv.append(j)
                break
    return piv


# ---------------------------------------------------------------- Jacobian ring
class JacobianRing:
    def __init__(self, F):
        self.F = F
        self.n = pdeg(F)
        self.e = self.n - 1
        self.dF = [pdiff(F, i) for i in range(3)]
        self._deg = {}
        self._lift = {}

    def gens(self, d):
        return [(m, i) for m in monomials(d - self.e) for i in range(3)] if d >= self.e else []

    def _data(self, d):
        if d in self._deg:
            return self._deg[d]
        mons = monomials(d)
        idx = {m: k for k, m in enumerate(mons)}
        G = self.gens(d)
        if G:
            rows = []
            for (m, i) in G:
                row = {}
                for mm, v in self.dF[i].items():
                    t = (m[0] + mm[0], m[1] + mm[1], m[2] + mm[2])
                    row[idx[t]] = v
                rows.append(row)
            A, r = to_fmpq_mat(rows, len(mons)).rref()
            piv = pivots_of_rref(A, r)
        else:
            A, r, piv = None, 0, []
        pset = set(piv)
        nonpiv = [k for k in range(len(mons)) if k not in pset]
        # N[i][s] = A[i, nonpiv[s]]
        N = [[A[i, s] for s in nonpiv] for i in range(r)] if r else []
        data = dict(mons=mons, idx=idx, rank=r, piv=piv, nonpiv=nonpiv, N=N)
        self._deg[d] = data
        return data

    def dim(self, d):
        return len(self._data(d)["nonpiv"])

    def std_monomials(self, d):
        D = self._data(d)
        return [D["mons"][k] for k in D["nonpiv"]]

    def reduce(self, p, d=None):
        """Coordinates (list of fmpq) of the class of p in R_d on the standard monomials."""
        if d is None:
            d = pdeg(p)
        D = self._data(d)
        idx = D["idx"]
        v = [fmpq(0)] * len(D["mons"])
        for m, c in p.items():
            if sum(m) != d:
                raise ValueError("degree mismatch")
            v[idx[m]] += c
        out = [v[k] for k in D["nonpiv"]]
        for i, pc in enumerate(D["piv"]):
            c = v[pc]
            if c != 0:
                Ni = D["N"][i]
                for s in range(len(out)):
                    out[s] -= c * Ni[s]
        return out

    def _lift_data(self, d):
        if d in self._lift:
            return self._lift[d]
        if self.dim(d) != 0:
            raise ValueError(f"R_{d} != 0, cannot lift")
        mons = monomials(d)
        idx = {m: k for k, m in enumerate(mons)}
        G = self.gens(d)
        # transpose: rows = monomials of S_d, cols = generators
        T = [[fmpq(0)] * len(G) for _ in mons]
        for g, (m, i) in enumerate(G):
            for mm, v in self.dF[i].items():
                t = (m[0] + mm[0], m[1] + mm[1], m[2] + mm[2])
                T[idx[t]][g] = v
        TM = to_fmpq_mat(T, len(G))
        A, r = TM.rref()
        assert r == len(mons)
        piv = pivots_of_rref(A, r)
        B = fmpq_mat(len(mons), len(piv), [TM[a, b] for a in range(len(mons)) for b in piv])
        Binv = B.inv()
        data = dict(mons=mons, idx=idx, gens=[G[b] for b in piv], Binv=Binv)
        self._lift[d] = data
        return data

    def lift(self, p):
        """Return [a0, a1, a2] with p = sum a_i dF/dx_i (requires R_deg(p) = 0)."""
        d = pdeg(p)
        L = self._lift_data(d)
        v = fmpq_mat(len(L["mons"]), 1, [p.get(m, fmpq(0)) for m in L["mons"]])
        sol = L["Binv"] * v
        a = [{}, {}, {}]
        for g, (m, i) in enumerate(L["gens"]):
            c = sol[g, 0]
            if c != 0:
                a[i][m] = a[i].get(m, fmpq(0)) + c
        # verify
        chk = {}
        for i in range(3):
            chk = padd(chk, pmul(a[i], self.dF[i]))
        assert padd(chk, p, -1) == {}, "lift failed"
        return a

    def c2(self, P, Fdot, Fddot):
        """Class in R of (1/2) div(a) - P*Fddot, where 2 P Fdot^2 = sum a_i dF/dx_i."""
        lhs = pscale(pmul(P, pmul(Fdot, Fdot)), fmpq(2))
        return self.c2_from_lift(P, self.lift(lhs), Fddot)

    def c2_from_lift(self, P, a, Fddot):
        """Class of (1/2) div(a) - P*Fddot for a given lift a = [a0, a1, a2]."""
        div = {}
        for i in range(3):
            div = padd(div, pdiff(a[i], i))
        val = padd(pscale(div, fmpq(1, 2)), pmul(P, Fddot), -1)
        return self.reduce(val, pdeg(P) + self.n)  # deg(div a) = deg(P) + 2n - (n-1) - 1 = deg(P) + n
