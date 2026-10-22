# Quantum Anticodes

A Python library implementing parts of the mathematical framework in
[Quantum Anticodes].

Currently, it can:

- Build finite-field symplectic vectors and linear codes, and compute their
  symplectic pairings, dimensions, radicals, nondegenerate complements, and isorank.
- Construct coordinate anticodes and restrict, puncture, or shorten codes.
- Construct binary CSS stabilizer codes, recover their classical X/Z generators,
  and construct CSS encoder Choi-state generator matrices.
- Compute classical Schur products, puncturing, and shortening, and test the
  `trans_DJ` Schur-power containment condition.

General operations support finite fields through `galois`; CSS and Choi
functions currently require `GF(2)`. The library works with algebraic code
representations rather than quantum circuits.

**Status:** This project is unfinished and under active development. Code and
documentation may be incomplete or lack professional polish, and interfaces
may change.
