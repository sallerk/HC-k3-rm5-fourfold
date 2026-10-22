Supporting code and data for the note
"The Hodge conjecture for the sextic fourfolds of an eight-dimensional family of K3 surfaces with real
multiplication"
==========================================================================================================

Status: accompanies the note (K. Saller, 2026); not peer reviewed.

The note is the sequel to "An explicit eight-dimensional family of K3 surfaces with real multiplication by
Q(sqrt5)" (the "K3 note", DOI 10.5281/zenodo.23309349), whose family it uses: F = F_c, c = (nu, b, eps) with
5 b0^2 = eps nu0, the eleven parameters nu1..nu6, b0..b3, eps (nu0 = 5 b0^2/eps eliminated), and the member c0
of Lemma 2.1(4) of the K3 note,
    eps = 5, nu = x0^6 + x0 x1^5 + x1^6, b = x0^3 + x1^3,
    4F = 4x0^6 + 40x0^3x1^2x2 - 20x0^3x2^3 + 4x0x1^5 - 20x0x1^4x2 + 4x1^6 - 20x1^3x2^3 + 4x1x2^5 + 5x2^6.

What the note proves (Theorem A): for very general c, the fourfold X = {F(x) + F(y) = 0} in P^5 has
dim Hdg^4(X) = 9, and the Hodge conjecture holds for X. Write Y = {w^6 + F = 0}, H_k = the zeta^k-eigenspace of
w -> zeta w on H^2(Y) (zeta = e^{2 pi i/6}); in the Jacobian ring R = C[x0,x1,x2]/J(F) the eigenspace H_{j+1}
corresponds to the coefficient of w^j, with Hodge pieces R_{2-j}, R_{8-j}, R_{14-j}. The chi-eigenspace of
y -> zeta y on H^4(X) is H_chi (x) H_{6-chi} for chi != 0, and the primitive (2,2)-part of H^4(X) is R_12 of the
Jacobian ring of F(x) + F(y), that is, the sum of R_a (x) R_b with a + b = 12 ("chains" a mod 6; chain 1-7 is
H_2 (x) H_4, chain 3-9 is chi = 0, chain 0-6-12 is H_3 (x) H_3 = H^2_prim(S) (x) H^2_prim(S)).

The proofs use the following computations (Section 7 of the note), all at c0:
(1) the first-order kernels K_chi(c0) of every chain: dimensions 1, 1, 22, 122, 22, 1 for chi = 0..5
    (Proposition 3.4) -- f5_chain_bounds_ex.py modulo 32003 and 31991; indep/chain_kernels.py exactly over Q;
(2) the solution spaces of the first- and second-order systems for chain 1-7 (Proposition 4.4): 22 at first
    order and 1 with two or more arcs -- f7_printed_arcs.py (the two arcs printed in the note) and
    f4_chain_second_order.py (random arcs) modulo 32003 and 31991; indep/deform_check.py exactly over Q.
The other scripts are further checks (not used in the proofs): see section 3.

Not checked by any script: the arguments of the note (Lemmas 2.1-2.2, 3.1-3.3, 4.1-4.2, Proposition 4.3, the
proof of Theorem A and the remarks) are proofs, not computations. In particular the step from the computed
kernel dimensions at c0 to the bounds for very general c (lower semicontinuity of the rank, reduction modulo p)
is an argument of the note.

Licence: the code in code/ is under the MIT licence; the outputs in outputs/ and this text are under CC BY 4.0
(see the LICENSE files of the repository). No third-party files are included.


1. Requirements
---------------
Python 3 with sympy and python-flint. Tested with Python 3.12.9, sympy 1.13.1 and python-flint 0.9.0. Every
script runs in a few seconds on a desktop PC.


2. How to run
-------------
The reference results are in outputs/. The scripts in code/ print their results to the screen (b10 also writes
a JSON file); the scripts in code/indep/ print their results and also save them in the current folder. From
the top folder of this repository:

    mkdir rerun
    cd rerun

Computations used in the proofs:

    python -I ../code/f5_chain_bounds_ex.py 32003 > f5_p32003.txt
    python -I ../code/f5_chain_bounds_ex.py 31991 > f5_p31991.txt
    python -I ../code/f4_chain_second_order.py 32003 0 4 small > f4_p32003_s0_4arcs.txt
    python -I ../code/f4_chain_second_order.py 31991 0 4 small > f4_p31991_s0_4arcs.txt
    python -I ../code/f4_chain_second_order.py 32003 0 2 small > f4_p32003_s0_2arcs.txt
    python -I ../code/f4_chain_second_order.py 31991 0 2 small > f4_p31991_s0_2arcs.txt
    python -I ../code/f4_chain_second_order.py 32003 2 2 small > f4_p32003_s2_2arcs.txt
    python -I ../code/f7_printed_arcs.py 32003 > f7_p32003.txt
    python -I ../code/f7_printed_arcs.py 31991 > f7_p31991.txt
    python -I ../code/indep/chain_kernels.py
    python -I ../code/indep/deform_check.py c0 6 20261011
    python -I ../code/indep/deform_check.py c0 6 1
    python -I ../code/indep/deform_check.py c0 6 2

(chain_kernels.py writes chain_kernels.txt; deform_check.py writes deform_check_c0.txt,
deform_check_c0_seed1.txt and deform_check_c0_seed2.txt.)

Further checks:

    python -I ../code/f1_v_annihilator.py 32003 0 > f1_p32003_s0.txt
    python -I ../code/f2_cube_on_Q.py > f2.txt
    python -I ../code/f3_second_order.py 32003 0 > f3_p32003_s0.txt
    python -I ../code/f6_split22.py 32003 > f6_p32003.txt
    python -I ../code/b10_ivhs_new.py b10_p32003_s0.json 32003 0 > b10_p32003_s0.txt
    python -I ../code/indep/deform_check.py even 6 20261011
    python -I ../code/indep/deform_check.py even 6 3
    python -I ../code/indep/chain_kernels.py 7

Then compare the files in rerun/ with those of the same name in outputs/ (for instance with "diff" or "fc").
They should be identical, except for the run times printed by the scripts in code/indep/.


3. What each script checks
--------------------------
Each script has a header with the details.

First implementation (code/; arithmetic modulo a prime p, python-flint nmod_mat):
ivhs_ts.py: library: polynomials as dictionaries, the Jacobian ring R (monomial bases of R_d by row
  reduction), and kernel_dims(R, R, tangents, p), the kernels of multiplication by T(x) + T(y) from R_12 of
  the Jacobian ring of F(x) + F(y) to R_18, chain by chain. (It can also be run on its own; usage in its
  header.)
b10_ivhs_new.py: the eleven-parameter family F_c; to_dict (reduction of a rational polynomial modulo p). As a
  script: the chain kernels at a random rational member of the family (expected 122, 22, 1, 1, 1, 22 for the
  chains 0-6-12, 1-7, ..., 5-11), and two controls: 8 random sextic directions at the same F (1 on every chain
  except 1 + 11^2 on 0-6-12), and the seven-dimensional family of van Geemen and Schuett (eps = b0 = 0).
f5_chain_bounds_ex.py: (1) at c0: dim R_d for d = 0..13 (R_13 = 0, so F is smooth modulo p), and the chain
  kernels K_chi(c0) with the eleven parameter derivatives as tangents. Also the first-order bounds for a random
  integral plane sextic with all 28 monomial directions (1 on every chain; Remark 6.1 of the note; its integer
  coefficients are printed).
f4_chain_second_order.py: (2), the system of Proposition 4.3 for chain 1-7: unknowns C (3 x 3, on S_1),
  Phibar (18 x 18, on R_7) and one D (3 x 3) per arc; equations (1) Phibar(P_j Fdot) = sum_i C_ij P_i Fdot for
  the eleven tangents Fdot, and (2) Phibar(c2(P_j)) = sum_i C_ij c2(P_i) - 2 sum_i D_ij P_i Fdot for each arc,
  c2(P) = (1/2) sum_i d a_i/d x_i - P Fddot with 2 P Fdot^2 = sum_i a_i dF/dx_i (Lemma 4.2). Prints the
  dimension of the projection of the solution space to (C, Phibar) at first order and with the arcs, and the
  dimension of the full solution space, for c0 ("ex") and a random rational member ("rnd"). With "small" the
  arc directions have entries in [-20, 20]. The arcs are drawn one after the other from the seed, so for c0 the
  run with 2 arcs uses the first two arcs of the run with 4 arcs (same prime and seed). Controls in the same run:
  8 random sextic directions at c0 (first order 1); the van Geemen-Schuett family (76 at first order, 20 with
  four arcs, 28 with two; Remark 6.3); the family of all sextics even in x2 (2 and 2: the identity and
  x2 -> -x2).
f7_printed_arcs.py: (2) with the two arcs d1, d2 of Proposition 4.4 (listed below under deform_check.py),
  using the function bound of f4: at c0 the first-order bound 22, then 4 with d1 and 1 with d1, d2 (also the
  dimension of the full solution space, including D).
f1_v_annihilator.py: the space P_5 of v in R_5 with Fdot v = 0 in R_11 for all eleven tangents (Remark 4.5):
  dimension 1 at c0 and at a random member, and Ann(v)_3 = <Q>, Q = (eps/2) x2^3 - b(x0, x1) (at c0 the
  printed element is x0^3 + x1^3 - (5/2) x2^3, as 15999 = -5/2 modulo 32003); 4 for the van Geemen-Schuett
  family.
f2_cube_on_Q.py: F restricted to the cubic curve Q = 0 is not a cube (resultant at c0; F modulo Q, exactly
  over Q(nu, b, eps)). Not used in the note; it records that v is not obviously the class of a curve.
f3_second_order.py: the pairing test of Remark 4.5: <v, c2(P)> for P in S_1 and six arcs; the space of v in
  P_5 that pass is 0 at c0 and at a random member, 1 for the van Geemen-Schuett family. Control: arcs
  F(e^{sM} x) (trivial deformations) give c2(P) = 0 in R_7 (18 of 18).
f6_split22.py: the 22 first-order solutions of f4 at c0: 18 of them have C = 0, and these are exactly the
  Phibar = r (x) v^dual (r in R_7, v as in f1); the projection to C has rank 4.

Second implementation (code/indep/; exact over Q, python-flint fmpq_mat; ranks modulo 1000003 are printed only
as a cross-check). It was written separately from the first, from the description of the systems in the note
(Sections 3 and 4), and imports nothing from code/:
jacring.py: library: normal forms in R_d with respect to the degree-reverse-lexicographic Groebner basis of J
  (x0 > x1 > x2), obtained degree by degree by row reduction; the lift sum a_i dF/dx_i = 2 P Fdot^2 in
  degree 13; exact kernels.
chain_kernels.py: (1) exactly over Q: K_chi(c0) for chi = 0..5 (1, 1, 22, 122, 22, 1); controls: no tangents
  (K_chi = all of L_chi: 200, 261, 333, 363, 333, 261) and a random integral sextic with all 28 monomial
  directions (1 for every chi; the coefficients are printed). The argument is the seed of that random sextic.
deform_check.py: (2) exactly over Q, with six arcs (directions with entries in [-20, 20], from the seed): the
  dimension of the solution space and of its projections to (C, Phi) and to C, at first order (22, of which 18
  with C = 0; projection to C of dimension 4) and with the first 1, 2, 4, 6 arcs (4, 1, 1, 1). The two arcs of
  Proposition 4.4 are the first two arcs for seed 20261011:
      d1 = [11, 2, 11, 1, -3, -17, -18, 2, 10, -18, -20],
      d2 = [2, -7, 18, 13, -12, 10, -19, -20, -6, -5, -11]
  in the coordinates (nu1, ..., nu6, b0, ..., b3, eps); the arc is c0 + s d. Controls in the same run: trivial
  deformations F(e^{sM} x) give c2 = 0 (and c2 != 0 if Fddot is doubled); c2(P) does not change when the a_i
  are changed by Koszul relations; (C, Phi, D) = (I, I, 0) solves all equations. With "even": the family of
  all sextics even in x2, which gives 2 at first and at second order.
