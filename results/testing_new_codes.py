from quantum_anticodes import scan_css_code
from results.exact_css_scan import scan_css_code_exact
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
result1 = scan_css_code(
    gx, gz,
    l=2,  
    parameters={"constructor": "C6Code", "arguments": {"l" : "2"}},
)
pprint.pprint(result1)


parameters = {
    "code_name": "C6Code",
    "Hx": code.matrix_x.tolist(),
    "Hz": code.matrix_z.tolist(),
    }

# quicker method
result2 = scan_css_code_exact(gx, gz, l=2, parameters=parameters)

pprint.pprint(result2)