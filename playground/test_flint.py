from flint import nmod_mat

# Define matrix elements and the modulus (e.g., mod 23)
matrix_data = [[1, 2, 3, 4, 5],
               [2, 4, 6, 8, 10],
               [3, 6, 8, 10, 12]]
modulus = 5

# Create the modular matrix
A = nmod_mat(matrix_data, modulus)

# Compute the basis for the nullspace
# Returns: X (a matrix whose first 'nullity' columns form the basis) and nullity (dimension)
X, nullity = A.nullspace()

print(f"Nullity (Dimension of Null Space): {nullity}")
print("Basis vectors (columns of X):")
print(X)