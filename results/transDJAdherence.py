from quantum_anticodes import scan_css_code
from qldpc import codes
import pprint

code = codes.C6Code()
# First grab the Hx and Hz stablizer check matrices. These are those that code words remain invariant
# under

# A Pauli operator keeps such a code word in the codespace if it commutes with the stablizers.
# gx only needs to commute with the z's as it commutes with itslef and vice versa. Then,
# gx and gz represesnet those generate valid paulis that keep words in the space

gx = code.matrix_z.null_space()
gz = code.matrix_x.null_space()

# can throw in other parameters as well
result = scan_css_code(
    gx, gz,
    l=2,  
    parameters={"constructor": "C6Code", "arguments": {"l" : "2"}},
)

pprint.pprint(result)