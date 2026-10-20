from quantum_anticodes import (
    GF,
    choice,
    puncture_classical,
    schur_product_codes,
    schur_product_vectors,
    shorten_classical,
    trans_DJ,
)
from qldpc import codes

code = codes.quantum.C6Code()

S = code.symplectic_check_matrix() 

n = code.length()
X_matrix = S[:, :n]
Z_matrix = S[:, n:]

print(code.matrix_x)
print(code.matrix_z)
print("___________________")

Lx, Lz = choice(X_matrix, Z_matrix)
print("\nChoice Lx:\n", Lx, sep="")
print("\nChoice Lz:\n", Lz, sep="")

print(f"c6 logical ops:\n {code.get_logical_ops()}")