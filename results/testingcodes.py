from quantum_anticodes import choi
import json
from results.exact_css_scan import scan_css_code_exact


def build_codes():
    """Build the listed small codes and one example per general family."""
    from qldpc import abstract, codes
    from sympy.abc import x, y

    code_cases = [
        ("C4Code()", codes.C4Code()),
        ("C6Code()", codes.C6Code()),
        ("SteaneCode()", codes.SteaneCode()),
        ("TetrahedralCode()", codes.TetrahedralCode()),
        ("TetrahedralCode(algebraic=True)", codes.TetrahedralCode(algebraic=True)),
        ("QuantumHammingCode(3)", codes.QuantumHammingCode(3)),
        ("QuantumHammingCode(4)", codes.QuantumHammingCode(4)),
        ("ManyHypercubeCode(1)", codes.ManyHypercubeCode(1)),
        ("QuantumReedMullerCode(0, 2)", codes.QuantumReedMullerCode(0, 2)),
        ("QuantumReedMullerCode(0, 3)", codes.QuantumReedMullerCode(0, 3)),
    ]
    for n in (2, 4, 6, 8, 10, 12, 14):
        code_cases.append((f"IcebergCode({n})", codes.IcebergCode(n)))
    for n in range(1, 16):
        code_cases.append((f"TrivialCode({n})", codes.TrivialCode(n)))

    # Both rectangular orientations are included. The conditions on n give
    # exactly the side lengths in the earlier list, including unrotated r=1.
    for rows in range(1, 16):
        for cols in range(1, 16):
            if rows >= 2 and cols >= 2 and rows * cols <= 15:
                name = f"SurfaceCode({rows}, {cols}, rotated=True)"
                code_cases.append((name, codes.SurfaceCode(rows, cols, rotated=True)))
            if rows * cols + (rows - 1) * (cols - 1) <= 15:
                name = f"SurfaceCode({rows}, {cols}, rotated=False)"
                code_cases.append((name, codes.SurfaceCode(rows, cols, rotated=False)))
            if rows % 2 == 0 and cols % 2 == 0 and rows * cols <= 15:
                name = f"ToricCode({rows}, {cols}, rotated=True)"
                code_cases.append((name, codes.ToricCode(rows, cols, rotated=True)))
            if 2 * rows * cols <= 15:
                name = f"ToricCode({rows}, {cols}, rotated=False)"
                code_cases.append((name, codes.ToricCode(rows, cols, rotated=False)))

    for size, dim in ((2, 2), (3, 2), (2, 3), (1, 2)):
        name = f"GeneralizedSurfaceCode({size}, {dim})"
        code_cases.append((name, codes.GeneralizedSurfaceCode(size, dim)))
    # Nonperiodic size=1 gives the same one-qubit code for every dim>=2;
    # (1, 2) above represents that otherwise unbounded duplicate family.
    for size, dim in [(2, 2)] + [(1, dim) for dim in range(2, 16)]:
        name = f"GeneralizedSurfaceCode({size}, {dim}, periodic=True)"
        code_cases.append((name, codes.GeneralizedSurfaceCode(size, dim, periodic=True)))

    # One explicit lattice of determinant 2 gives a 12-qubit T4 example.
    lattice = [[2, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    code_cases.append(("T4Code(diag(2, 1, 1, 1))", codes.T4Code(lattice)))

    # One explicit binary example per general construction family.
    ring2 = abstract.GroupRing(abstract.CyclicGroup(2), field=2)
    ring3 = abstract.GroupRing(abstract.CyclicGroup(3), field=2)
    z = ring3.generators[0]
    quaternion_group = abstract.QuaternionGroup()
    i, j = quaternion_group.generators
    code_cases.extend([
        ("CSSCode([[1,1,1,1]], [[1,1,1,1]])",
         codes.CSSCode([[1, 1, 1, 1]], [[1, 1, 1, 1]], field=2)),
        ("TBCode([[1,1],[1,1]], [[1,1],[1,1]])",
         codes.TBCode([[1, 1], [1, 1]], [[1, 1], [1, 1]], field=2)),
        ("QCCode([3], 1+x, 1+x)", codes.QCCode([3], 1 + x, 1 + x, field=2)),
        ("BBCode((2,3), 1+x, 1+y)", codes.BBCode((2, 3), 1 + x, 1 + y, field=2)),
        ("HGPCode(RepetitionCode(3))",
         codes.HGPCode(codes.RepetitionCode(3, field=2), field=2)),
        ("CHGPCode((2,3), 1+x)", codes.CHGPCode((2, 3), 1 + x, field=2)),
        ("CRCode(2, 1+x)", codes.CRCode(2, 1 + x, field=2)),
        ("LPCode([[1,1]] over GF(2)[CyclicGroup(2)])",
         codes.LPCode(abstract.RingArray.build([[1, 1]], ring2))),
        ("GALACode([1,z], [z**2,1], 1) over GF(2)[CyclicGroup(3)]",
         codes.GALACode([ring3.one, z], [z**2, ring3.one], num_active_rows=1)),
        ("QTCode(QuaternionGroup, [i,~i], [j,~j], RepetitionCode(2), bipartite=False)",
         codes.QTCode([i, ~i], [j, ~j], codes.RepetitionCode(2, field=2),
                      field=2, bipartite=False)),
    ])

    # Keep the existing catalogue order. Aliases remain separate.
    code_cases.sort(key=lambda item: len(item[1]) + item[1].dimension)
    return code_cases


def main(code_cases=None, *, results_l1=None, results_l2=None):
    """Fill and return two result dictionaries; print/save nothing.

    Supplying dictionaries keeps completed results accessible if interrupted.
    """
    if code_cases is None:
        code_cases = build_codes()

    if results_l1 is None:
        results_l1 = {}
    if results_l2 is None:
        results_l2 = {}
    if results_l1 is results_l2:
        raise ValueError("Use separate result dictionaries for l=1 and l=2")
    for code_name, code in code_cases:
        if code_name in results_l1 or code_name in results_l2:
            raise ValueError(f"Duplicate code name: {code_name}")
        Hx = code.matrix_x
        Hz = code.matrix_z
        if not 1 <= Hx.shape[1] <= 15:
            raise ValueError(f"{code_name} must have between 1 and 15 physical qubits")
        Gx = Hz.null_space()
        Gz = Hx.null_space()
        Lx, Lz = choi(Gx, Gz)

        parameters = {
            "code_name": code_name,
            "Hx": Hx.tolist(),
            "Hz": Hz.tolist(),
        }
        results_l1[code_name] = scan_css_code_exact(Gx, Gz, l=1, parameters=parameters)
        results_l2[code_name] = scan_css_code_exact(Gx, Gz, l=2, parameters=parameters)

    with open("scan_results.json", "w") as file:
        json.dump(
            {"l1": results_l1, "l2": results_l2},
            file,
            indent=2,
        )
        

    passing_results_l1 = {}
    passing_results_l2 = {}

    for code_name, result in results_l1.items():
        if result["has_nonempty_pass"]:
            passing_results_l1[code_name] = result

    for code_name, result in results_l2.items():
        if result["has_nonempty_pass"]:
            passing_results_l2[code_name] = result

    with open("scan_results_nonempty_l1.json", "w") as file:
        json.dump(passing_results_l1, file, indent=2)

    with open("scan_results_nonempty_l2.json", "w") as file:
        json.dump(passing_results_l2, file, indent=2)




    return results_l1, results_l2


if __name__ == "__main__":
    # These remain available in interactive Python even after Ctrl-C.
    results_l1 = {}
    results_l2 = {}
    main(results_l1=results_l1, results_l2=results_l2)

