import itertools
from sympy.abc import x, y
from qldpc import codes

def get_canonical_shift(poly_exponents, ell, m):
    """
    Returns the lexicographically smallest cyclic shift of a polynomial.
    """
    all_shifts = []
    for origin_i, origin_j in poly_exponents:
        shifted_poly = []
        for i, j in poly_exponents:
            new_i = (i - origin_i) % ell
            new_j = (j - origin_j) % m
            shifted_poly.append((new_i, new_j))
        
        shifted_poly.sort()
        all_shifts.append(tuple(shifted_poly))
        
    all_shifts.sort()
    return all_shifts[0]

def generate_canonical_polynomials(ell, m, weight=3):
    """
    Generates the canonical shift representations of weight-3 polynomials.
    """
    grid = [(i, j) for i in range(ell) for j in range(m) if (i, j) != (0, 0)]
    canonical_polys = set()
    
    for terms in itertools.combinations(grid, weight - 1):
        poly = ((0, 0),) + terms
        canonical = get_canonical_shift(poly, ell, m)
        canonical_polys.add(canonical)
        
    # Sort to guarantee deterministic generation order
    return sorted(list(canonical_polys))

def to_sympy(poly_tuple):
    """
    Converts a tuple of coordinates back to a sympy expression for qldpc.
    """
    return sum((x**i) * (y**j) for i, j in poly_tuple)

def generate_small_bb_codes(max_n=20, weight=3):
    """
    Generator that yields (code, n, k, ell, m, A_expr, B_expr) 
    for all unique Bivariate Bicycle codes up to max_n physical qubits.
    """
    max_lm = max_n // 2
    
    # We restrict to ell <= m because an (ell x m) grid generates the exact 
    # same codes as an (m x ell) grid, just with the x and y axes swapped.
    for ell in range(2, max_lm + 1):
        for m in range(ell, max_lm // ell + 1):
            n = 2 * ell * m
            
            # 1. Get the reduced polynomial set
            canonical_polys = generate_canonical_polynomials(ell, m, weight=weight)
            
            # 2. Convert to sympy ONCE per polynomial to avoid overhead in the inner loop
            sympy_polys = [to_sympy(p) for p in canonical_polys]
            
            # 3. combinations_with_replacement safely avoids testing (B, A) if (A, B) was tested
            for A_expr, B_expr in itertools.combinations_with_replacement(sympy_polys, 2):
                
                # Instantiate the qLDPC code
                code = codes.BBCode({x: ell, y: m}, A_expr, B_expr)
                
                # Calculate dimension k (fast integer matrix rank check)
                k = code.dimension
                
                # Yield the object for external testing pipelines
                yield code, n, k, ell, m, A_expr, B_expr

if __name__ == "__main__":
    print("Searching for unique Bivariate Bicycle codes with n <= 20...")
    
    # Example Testing Pipeline
    for code, n, k, ell, m, A, B in generate_small_bb_codes(max_n=20):
        
        # Test 1: Only pass codes that successfully encode logical qubits
        if k > 0:
            
            # Test 2: Evaluate the exact distance using Integer Linear Programming (ILP).
            # For codes this small, ILP solves in milliseconds. 
            d = code.get_distance() 
            
            # Test 3: Filter out trivial distance-1 codes
            if d > 1:
                print(f"[[{n}, {k}, {d}]] on {ell}x{m} | A = {A}, B = {B}")

