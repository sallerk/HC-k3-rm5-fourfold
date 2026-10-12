"""chain_kernels.py -- dimensions of K_chi in the fourfold Jacobian ring R (x) R at c0 (independent code).

Usage:
    python -I chain_kernels.py [seed]        (default seed 20261011; output -> chain_kernels.txt)

Setting: R = S/J for the plane sextic F (R_d = 0 for d >= 13), standard-monomial bases of R_d from
jacring.py (grevlex normal forms).  For chi in Z/6:
    L_chi = sum_{a+b=12, b = chi+3 mod 6} R_a (x) R_b,
    K_chi = { lam in L_chi : (T_k(x) + T_k(y)) lam = 0 in (R (x) R)_18, k = 1..m },
where T_k(x) acts on the first factor (R_a -> R_{a+6}), T_k(y) on the second (R_b -> R_{b+6}),
components of degree > 12 vanish.  dim K_chi = dim L_chi - rank of the stacked linear map.
Ranks are exact over Q (flint.fmpq_mat); ranks mod 1000003 are printed as a cross-check.

Runs: c0 with the 11 tangents T_k of deform_check.family_c0; control (i) no tangents (K = L);
control (ii) a random integral sextic (R_13 = 0 checked) with all 28 sextic monomials as tangents.
"""
import os
import random
import sys
import time

from flint import fmpq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from jacring import JacobianRing, monomials, pmul, rank_mod_p, to_fmpq_mat  # noqa: E402
from deform_check import Tee, family_c0  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))


def mult_matrices(R, T):
    """M[a][i'][i] = coefficient of std monomial i' of R_{a+6} in [T * (std monomial i of R_a)], a = 0..6."""
    out = {}
    for a in range(0, 7):
        std = R.std_monomials(a)
        cols = [R.reduce(pmul({m: fmpq(1)}, T), a + 6) for m in std]  # each: list over R_{a+6}
        out[a] = cols  # out[a][i][i']
    return out


def kernel_dims(R, tangents, log, label):
    dims = [R.dim(d) for d in range(13)]
    if tangents:
        rT = to_fmpq_mat([R.reduce(T, 6) for T in tangents], dims[6]).rank()
        log(f"[{label}] rank of the tangent classes [T_k] in R_6: {rT} (of {dims[6]})")
    mult = [mult_matrices(R, T) for T in tangents]
    res = []
    for chi in range(6):
        t0 = time.perf_counter()
        comps = [(12 - b, b) for b in range(13) if b % 6 == (chi + 3) % 6]
        offL, n = {}, 0
        for (a, b) in comps:
            offL[(a, b)] = n
            n += dims[a] * dims[b]
        # target blocks (a', b'), a' + b' = 18, a', b' <= 12
        tblocks = sorted({(a + 6, b) for (a, b) in comps if a + 6 <= 12}
                         | {(a, b + 6) for (a, b) in comps if b + 6 <= 12})
        offT, nt = {}, 0
        for (a, b) in tblocks:
            offT[(a, b)] = nt
            nt += dims[a] * dims[b]
        rows = []
        for Mk in mult:
            block = [dict() for _ in range(nt)]
            for (a, b) in comps:
                na, nb = dims[a], dims[b]
                for i in range(na):
                    for j in range(nb):
                        col = offL[(a, b)] + i * nb + j
                        if a + 6 <= 12:  # x-action: [T r_i] (x) s_j
                            base = offT[(a + 6, b)]
                            for ip, c in enumerate(Mk[a][i]):
                                if c != 0:
                                    r = base + ip * nb + j
                                    block[r][col] = block[r].get(col, fmpq(0)) + c
                        if b + 6 <= 12:  # y-action: r_i (x) [T s_j]
                            base = offT[(a, b + 6)]
                            nb6 = dims[b + 6]
                            for jp, c in enumerate(Mk[b][j]):
                                if c != 0:
                                    r = base + i * nb6 + jp
                                    block[r][col] = block[r].get(col, fmpq(0)) + c
            rows.extend(block)
        if rows and n:
            M = to_fmpq_mat(rows, n)
            rk = M.rank()
            rkp = rank_mod_p(M)
        else:
            rk = rkp = 0
        res.append(n - rk)
        log(f"[{label}] chi={chi}: components (a,b)={comps} dim L={n} | targets {tblocks} "
            f"rows={len(rows)} rank_Q={rk} rank_mod_p={rkp} | dim K={n - rk} | {time.perf_counter() - t0:.2f}s")
    log(f"[{label}] dim K_chi, chi=0..5: {res}   (symmetry K_chi = K_-chi: {res[1] == res[5] and res[2] == res[4]})")
    return res


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261011
    log = Tee("chain_kernels.txt" if seed == 20261011 else f"chain_kernels_seed{seed}.txt")
    T0 = time.perf_counter()
    log(f"chain_kernels.py seed={seed} (exact over Q; mod p=1000003 rank only as cross-check)")

    fam = family_c0()
    R = JacobianRing(fam["F"])
    dims = [R.dim(d) for d in range(14)]
    log(f"c0: 4F check {fam['check4F']}; dim R_d d=0..13: {dims}")
    log("control (i), no tangents:")
    kernel_dims(R, [], log, "c0, no tangents")
    log("c0 with the 11 tangents (nu1..nu6, b0..b3, eps):")
    kernel_dims(R, fam["T"], log, "c0, 11 tangents")

    rng = random.Random(seed)
    tries = 0
    while True:
        tries += 1
        coeffs = [rng.randint(-5, 5) for _ in monomials(6)]
        G = {m: fmpq(c) for m, c in zip(monomials(6), coeffs) if c != 0}
        RG = JacobianRing(G)
        if all(sum(m) == 6 for m in G) and RG.dim(13) == 0:
            break
    log(f"control (ii): random sextic ({tries} tries), coefficients on grevlex monomials of S_6: {coeffs}")
    log(f"  dim R_d d=0..13: {[RG.dim(d) for d in range(14)]}")
    kernel_dims(RG, [{m: fmpq(1)} for m in monomials(6)], log, "random sextic, 28 monomials")
    log(f"total time {time.perf_counter() - T0:.2f}s")


if __name__ == "__main__":
    main()
