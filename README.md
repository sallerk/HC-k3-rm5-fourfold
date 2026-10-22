# HC-k3-rm5-fourfold

Code and data supporting the note *The Hodge conjecture for the sextic fourfolds of an eight-dimensional family of K3
surfaces with real multiplication* (K. Saller, 2026). It is the sequel to the note *An explicit eight-dimensional
family of K3 surfaces with real multiplication by Q(√5)* (the "K3 note", DOI
[10.5281/zenodo.23309349](https://doi.org/10.5281/zenodo.23309349); code in the repository
[HC-k3-rm5-maximal](https://github.com/sallerk/HC-k3-rm5-maximal)).

**Result.** Let F be a very general member of the eight-dimensional family of plane sextics of the K3 note, whose
double planes w^2 = F are K3 surfaces with real multiplication by Q(√5). The note proves the Hodge conjecture for the
sextic fourfold X = {F(x_0,x_1,x_2) + F(y_0,y_1,y_2) = 0} in P^5 and shows that its rational Hodge classes of degree
four form a space of dimension 9 (7 for a very general plane sextic). By the construction of Shioda and Katsura the
problem reduces to Hodge classes in tensor products H_k ⊗ H_{6-k} of the eigenspaces of the sextic surface
w^6 + F = 0. The first-order (Carlson–Griffiths) bounds settle all but two of these pieces; for H_2 ⊗ H_4 and
H_4 ⊗ H_2 the note uses a second-order condition (endomorphisms of H_2 that preserve the holomorphic forms).

The computations used in the proofs are kernel dimensions of explicit linear systems in the Jacobian ring of F at an
explicit member, computed modulo two primes (first implementation, `code/`) and exactly over Q (second, independently
written implementation, `code/indep/`). [`README.txt`](README.txt) gives the requirements, the rerun instructions and
what each script checks.

**Status:** not yet peer reviewed.

**Disclosure:** the arguments and computations were produced with the assistance of an AI system (Claude, by
Anthropic). They were checked with separately written code and in a separate internal review, which were also carried
out with AI assistance. What was not checked is listed in README.txt.

## Licences

- **Code** (`code/`): MIT License, see [`LICENSE`](LICENSE).
- **Everything else** (README files and computed outputs in `outputs/`): Creative Commons Attribution 4.0
  International (CC BY 4.0), see [`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0).
- No third-party files are included.
