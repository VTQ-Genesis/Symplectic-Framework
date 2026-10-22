import numpy as np, galois
GF = galois.GF(2) # your field
G1 = GF([[1,1,0,0],[0,0,1,1],[1,0,1,0]]) # generator of C1 (rows)
G2 = GF([[1,1,1,1],[1,0,1,0],[1,1,0,0]]) # generator of C2 (rows)

print(G1.null_space())
print(G2.null_space())


H2 = G2.null_space()# basis of C2^perp
assert np.linalg.matrix_rank(np.vstack([G1, H2])) == np.linalg.matrix_rank(G1) # C2^perp ⊆ C1

# A1: rows of G1 extending span(H2) to C1
M, A1 = H2, []
for g in G1:
    N = np.vstack([M, g])
    if np.linalg.matrix_rank(N) > np.linalg.matrix_rank(M):
        M, A1 = N, A1 + [g]
A1 = np.vstack(A1); k = A1.shape[0]

Gt1 = np.vstack([np.hstack([H2, GF.Zeros((H2.shape[0], k))]),
    np.hstack([A1, GF.Identity(k)])])
Gt2 = Gt1.null_space() # = [[H1,0],[A2,-I]] up to row operations

print(Gt1); print(Gt2)