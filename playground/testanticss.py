"""Run from the repository root:
PYTHONPATH=src .venv/bin/python practice/010-css-construction/test.py
"""

from quantum_anticodes import *


F2 = GF(2)
code = linear_code([
    SymplecticVector((e, 0, e)),
    SymplecticVector((0, f, 0)),
], F2, 3)
J = [1, 2] 

# instead of looping through the interator use * to print out all
print("Original basis:", *code.basis)

# Project every word onto J, keep the first and second component
punctured = code.puncture(J)
print("Punctured basis:", *punctured.basis)

# Keep words zero outside J, then project onto J
shortened = code.shorten(J)
print("Shortened basis:", *shortened.basis)

# Keep words zero outside J, retaining all 3 coordinate positions
restricted = code.restrict(J)
print("Restricted basis:", *restricted.basis)

# Every support and restriction
for support, restriction in code.restrictions():
    print("Restriction", support, "basis:", *restriction.basis)


# Construct all that are 0 outside J
ambient = anticode(J, F2, code.n)
print("Anticode basis:", *ambient.basis)
# gets all of the supports
for support, ambient_anticode in anticodes(F2, code.n):
    print("Anticode", support, "basis:", *ambient_anticode.basis)


# Classical generator matrices with G1 @ G2.T == 0 over F2.
G1 = F2([[1, 1, 0]])
G2 = F2([[1, 1, 0]])
S = css_construct(G1, G2)  # Stabilizer subspace e*C1 + f*C2.
print("CSS stabilizer basis:", *S.basis)

# Recover classical generator matrices.
H1, H2 = reverse_css(S)
print("Recovered G1:\n", H1)
print("Recovered G2:\n", H2)
