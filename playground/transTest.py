from quantum_anticodes import (
    GF,
    choi,
    puncture_classical,
    schur_product_codes,
    schur_product_vectors,
    shorten_classical,
    trans_DJ,
)

from qldpc import codes
F2 = GF(2)
code = codes.C6Code()


hx = code.matrix_x
hz = code.matrix_z

gx = hz.null_space()
gz = hx.null_space()
print("\nChoice gxx:\n", gx, sep="")
print("\nChoice gz:\n", gz, sep="")

Lx, Lz = choi(gx,gz)


print("\nChoice Lx:\n", Lx, sep="")
print("\nChoice Lz:\n", Lz, sep="")

print(Lz.null_space())
print(Lx.row_space())



# 2. Extract X-type logicals (Rows 0 and 2) and take the first 6 columns (the X half)
#logical_X = logical_ops[0::2, :6]

# 3. Extract Z-type logicals (Rows 1 and 3) and take the last 6 columns (the Z half)
#logical_Z = logical_ops[1::2, 6:]

#print(logical_X)
#print(logical_Z)




C1 = F2([[0, 0, 1, 1], [1,0,0,0]])
C2 = F2([[1,0,0,0], [0,1,0,0],[0,0,1,1]])
J = [1, 2]
l = 5


#Lx, Lz = choice(logical_X, logical_Z)
#print("\nChoice Lx:\n", Lx, sep="")
#print("\nChoice Lz:\n", Lz, sep="")

#print(Lx.row_space())
#print(Lz.null_space())
#print(code)
#print(logical_ops)
#result1 = trans_DJ(Lx, [1, 2], 1)


#result2 = trans_DJ(C1, [1, 2], 5)
#print(result1)
#print("trans_DJ with l=5:", result2)

