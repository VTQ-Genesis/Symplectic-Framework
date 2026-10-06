# library for finte field computation
from itertools import product

from galois import GF

from quantum_anticodes import (
    LinearCode,
    SymplecticVector,
    VElement,
    linear_code,
    omega,
)


# define the field to work over
F = GF(2)

#define the basis vectors for each coord plane
e = VElement(1, 0, field=F)
f = VElement(0, 1, field=F)
#zero = VElement(0, 0, field=F)


# define arbitrary basis for code
b1 = SymplecticVector((e, 0,e+f), field=F)
b2 = SymplecticVector((f, 0,e+f), field=F)
b3 = SymplecticVector((0, e,f), field=F)

# print out the generators
print("Generators:")
print("b1 =", b1)
print("b2 =", b2)
print("b3 =", b3)
print("____________________________________________")

# compute pairings between each vector
print("\nPairings:")
print("omega(b1, b2) =", omega(b1, b2))
print("omega(b1, b3) =", omega(b1, b3))
print("omega(b2, b3) =", omega(b2, b3))
print("____________________________________________")

# compute the linear code, basis, give the dimension of the symplectic space and field
C = LinearCode([b1, b2, b3], field=F, n=3)

# print vector by vector as we dont compute it entirely
print("\nStored basis:")
for vector in C.basis:
    print(vector)
print("____________________________________________")


# print the field of the code and the symplectic dim
print("\nField:", C.field)
print("Length n:", C.n)

print("____________________________________________")



# dim over field, #symp pairs, isorank, #codewords
print("\ndimension_over_field =", C.dimension_over_field)
print("dimension            =", C.dimension)
print("isorank              =", C.isorank)
print("cardinality          =", C.cardinality)

print("____________________________________________")


# membership test
#inside = b1 + 2 * b2
#outside = SymplecticVector((0, f), field=F)

#print("\nMembership:")
#print(inside, "in C?", inside in C)
#print(outside, "in C?", outside in C)
#print("____________________________________________")



# compute the radical and symplectic split of the code
R, K = C.symplectic_decomposition()

print("\nRadical basis:")
for vector in R.basis:
    print(vector)
print("dim_F(radical) =", R.dimension_over_field)

print("____________________________________________")

print("\nSymplectic-part basis:")
for vector in K.basis:
    print(vector)
print("dim_F(symplectic part) =", K.dimension_over_field)

print("____________________________________________")

# print all codewords, we loop as we dont compute all
print("\nAll codewords:")
for word in C.codewords():
    print(word)

print("____________________________________________")








###############################################################
# The stuff below are just tests
class ClassicalCode:
    """The F_2-linear span of generators in F_2^n."""

    def __init__(self, generators, n):
        matrix = GF(2).Zeros((len(generators), n))
        for row, vector in enumerate(generators):
            matrix[row, :] = vector
        self.basis = matrix.row_space()
        self.field = GF(2)
        self.n = n

    @property
    def dimension(self):
        return self.basis.shape[0]

    @property
    def cardinality(self):
        return 2 ** self.dimension

    def codewords(self):
        for coefficients in product((0, 1), repeat=self.dimension):
            word = GF(2).Zeros(self.n)
            for coefficient, vector in zip(coefficients, self.basis):
                word = word + coefficient * vector
            yield word


def anticode(support: set[int], n: int, quantum: bool = True):
    if quantum:
        units = (VElement(1, 0, field=F), VElement(0, 1, field=F))
        generators = [
            SymplecticVector((unit if i == j else 0 for i in range(n)))
            for j in sorted(support)
            for unit in units
        ]
        return LinearCode(generators, field=F, n=n)
    generators = [
        [1 if i == j else 0 for i in range(n)]
        for j in sorted(support)
    ]
    return ClassicalCode(generators, n)

J = {0, 2}
A = anticode(J, n=3)
AR, AK = A.symplectic_decomposition()

B=anticode(J, n=3, quantum=False)

for vector in B.basis:
    print(vector)

print("____________________________________________")
