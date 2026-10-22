from quantum_anticodes import (
    GF,
    choi,
    puncture_classical,
    schur_product_codes,
    schur_product_vectors,
    shorten_classical,
    trans_DJ,
)


F2 = GF(2)

# Schur product of two vectors.
u = F2([1, 0, 1])
v = F2([0, 1, 1])
vector_product = schur_product_vectors(u, v)
print("Vector Schur product:", vector_product)

# Schur product of two codes, here the same code with itself.
G = F2([[1, 0, 1], [0, 1, 1]])
code_product = schur_product_codes(G, G)
print("\nCode Schur product:\n", code_product, sep="")

# J contains one-based labels of the coordinates to KEEP.
J = [1, 2]
punctured = puncture_classical(G, J)
shortened = shorten_classical(G, J)
print("\nClassical puncturing:\n", punctured, sep="")
print("\nClassical shortening:\n", shortened, sep="")

# Choice takes two binary generator matrices and returns two matrices.
Gx = F2([[1,1,1,1], [1,1,0,0]])
Gz = F2([[1,1,0,0], [0,0,1,1],[1,0,1,0]])
Lx, Lz = choi(Gx, Gz)
print("\nChoi Lx:\n", Lx, sep="")
print("\nChoi Lz:\n", Lz, sep="")

# trans_DJ function
C1 = F2([[1, 0, 1, 1], [0, 0, 1, 1]])
result1 = trans_DJ(C1, [1, 2], 1)
result2 = trans_DJ(C1, [1, 2], 5)
print("\ntrans_DJ with l=1:", result1)
print("trans_DJ with l=5:", result2)