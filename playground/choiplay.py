from quantum_anticodes import GF, choi, css_construct, reverse_css, trans_DJ

F2 = GF(2)

# Build the four-qubit code with stabilizers XXXX and ZZZZ.
Hx = F2([[1, 1, 1, 1]])
Hz = F2([[1, 1, 1, 1]])

S = css_construct(Hx, Hz)
print("Physical qubits:", S.n)
print("Independent stabilizers:", S.dimension_over_field)
print("Logical qubits:", S.n - S.dimension_over_field)

# Four physical qubits and two independent stabilizers leave two logical qubits.

# Recover the classical X and Z stabilizer generators.
recovered_Hx, recovered_Hz = reverse_css(S)
print("\nRecovered Hx:\n", recovered_Hx, sep="")
print("Recovered Hz:\n", recovered_Hz, sep="")
# Both matrices are [[1, 1, 1, 1]], recovering the original stabilizer spaces.

# Choi takes the opposite duals: these include the logical directions.
Gx = recovered_Hz.null_space()
Gz = recovered_Hx.null_space()
print("\nGx = Hz-perp:\n", Gx, sep="")
print("Gz = Hx-perp:\n", Gz, sep="")
# Both have three independent rows and generate the length-four even-weight code.

# Construct the X/Z generators of one Choi state, then its symplectic code.
Lx, Lz = choi(Gx, Gz)
choi_state = css_construct(Lx, Lz)
print("\nChoi Lx:\n", Lx, sep="")
print("Choi Lz:\n", Lz, sep="")
print("X/Z pairings:\n", Lx @ Lz.T, sep="")
print("Choi qubits:", choi_state.n)
print("Choi independent stabilizers:", choi_state.dimension_over_field)
# The pairings are zero. Six independent stabilizers describe one pure state
# on four physical qubits plus two reference qubits (positions 5 and 6).

# Test the X-sector; J lists one-based positions to KEEP.
J = [1, 2, 3, 4]
l = 2
result = trans_DJ(Lx, J, l)
print("\nFour-qubit example, physical J, l=2:", result)
# False: this valid Choi state does not satisfy the extra Schur-power containment.
# l=2 means two Schur factors; it is not ordinary entrywise matrix squaring.

# A second code has two physical qubits, stabilizer ZZ, and one logical qubit.
repetition = css_construct(F2.Zeros((0, 2)), F2([[1, 1]]))
rep_Hx, rep_Hz = reverse_css(repetition)
rep_Gx = rep_Hz.null_space()
rep_Gz = rep_Hx.null_space()
rep_Lx, rep_Lz = choi(rep_Gx, rep_Gz)
print("\nRepetition Choi Lx:\n", rep_Lx, sep="")
print("Repetition Choi Lz:\n", rep_Lz, sep="")
# Lx is [[1, 1, 1]]. Together these matrices describe the GHZ state
# (|000> + |111>)/sqrt(2), with reference qubit at position 3.

physical_result = trans_DJ(rep_Lx, [1, 2], 2)
full_result = trans_DJ(rep_Lx, [1, 2, 3], 2)
print("Repetition, physical J, l=2:", physical_result)
print("Repetition, physical + reference J, l=2:", full_result)
# True for the two physical positions, False when the reference is also kept.
# Changing J changes the containment being tested; it does not change the state.
