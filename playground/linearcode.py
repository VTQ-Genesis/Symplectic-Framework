from quantum_anticodes import GF, LinearCode, SymplecticVector, VElement, linear_code

# work over the finite field of 3
F = GF(3)

# define the e and f elements that act as the symplectic pairs
e = VElement(1, 0, field=F)
f = VElement(0, 1, field=F)

# define the 2 dimensional symplectic vector that will act as the geenrators for the code
u = SymplecticVector((e, 0))
v = SymplecticVector((0, f))

# create the code
code = linear_code([u, v, u + v], field=F, n=2)

print(code.dimension_over_field)  # actual basis of the code, ie number of rows
print(code.cardinality)           # num of codewords
print(u + v in code)              # checks inclusion

for word in code.codewords():
    print(word)
