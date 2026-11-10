from quantum_anticodes import scan_css_code, choi, trans_DJ
from results.exact_css_scan import scan_css_code_exact
from qldpc import codes
import pprint


code = codes.SurfaceCode(2, 7, rotated=False)
print(code.matrix_x)
# print()
print(code.matrix_z)
print(code.dimension)
print(code.get_distance())
print(len(code))
print(code.get_logical_ops())

Lx, Lz = choi(code.matrix_z.null_space(),code.matrix_x.null_space())

#print(trans_DJ(Lx, [1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13], l=2))


# First grab the Hx and Hz stablizer check matrices. These are those that code words remain invariant
# under

# A Pauli operator keeps such a code word in the codespace if it commutes with the stablizers.
# gx only needs to commute with the z's as it commutes with itslef and vice versa. Then,
# gx and gz represesnet those generate valid paulis that keep words in the space

gx = code.matrix_z.null_space()
gz = code.matrix_x.null_space()

scan


# 2-> t, 1- > s

# can throw in other parameters as well
# result = scan_css_code(
#     gx, gz,
#     l=1,  
#     parameters={"constructor": "SurfaceCode(2, 7, rotated=True)", "arguments": {"l" : "1"}},
# )

# pprint.pprint(result)

