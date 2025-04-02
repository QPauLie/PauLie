from paulie.common.pauli_string_factory import get_pauli_string
from paulie.application.pauli_compiler import pauli_compiler, nested_adjoint 


def test_pauli_compiler():
    """Test the PauliCompiler with examples"""
    # Define parameters
    N = 4  # Total number of qubits
    k = 3  # Number of qubits in first subsystem (V acts on k qubits)
    
    # Create the example sets
    # A is a set of Pauli strings on k qubits
    A_set = [
        get_pauli_string("XII"),  # X₁
        get_pauli_string("IXI"),  # X₂
        get_pauli_string("IIX"),  # X₃
        get_pauli_string("ZII"),  # Z₁
        get_pauli_string("IZI"),  # Z₂
        get_pauli_string("IIZ"),  # Z₃
        get_pauli_string("YII"),  # Y₁
        get_pauli_string("IYI"),  # Y₂
        get_pauli_string("IIY")   # Y₃
    ]
    
    # B is a set of Pauli strings on (N-k) qubits
    B_set = [
        get_pauli_string("X"),  # X₄
        get_pauli_string("Z"),  # Z₄
        get_pauli_string("Y")   # Y₄
    ]
    
    # A_prime is A⊗I^{⊗(N-k)}
    A_prime_set = [
        get_pauli_string("XIII"),  # X₁⊗I₄
        get_pauli_string("IXII"),  # X₂⊗I₄
        get_pauli_string("IIXI"),  # X₃⊗I₄
        get_pauli_string("ZIII"),  # Z₁⊗I₄
        get_pauli_string("IZII"),  # Z₂⊗I₄
        get_pauli_string("IIZI"),  # Z₃⊗I₄
        get_pauli_string("YIII"),  # Y₁⊗I₄
        get_pauli_string("IYII"),  # Y₂⊗I₄
        get_pauli_string("IIYI")   # Y₃⊗I₄
    ]
    
    # B_prime is U_B⊗B
    B_prime_set = [
        get_pauli_string("XIIX"),  # X₁⊗X₄
        get_pauli_string("XIIZ"),  # X₁⊗Z₄
        get_pauli_string("XIIY"),  # X₁⊗Y₄
        get_pauli_string("ZIIX"),  # Z₁⊗X₄
        get_pauli_string("ZIIZ"),  # Z₁⊗Z₄
        get_pauli_string("ZIIY")   # Z₁⊗Y₄
    ]
    
    print("=" * 50)
    print(f"Testing PauliCompiler with {N} qubits, k={k}")
    print("=" * 50)
    
    # Example 1: XYZI
    print("\n" + "=" * 50)
    print("EXAMPLE 1: Compiling XYZI")
    print("=" * 50)
    
    target1 = get_pauli_string("XYZI")
    sequence1 = pauli_compiler(target1, A_set, A_prime_set, B_set, B_prime_set, k, N)
    
    print("\nCompiled sequence:")
    for i, pauli in enumerate(sequence1):
        print(f"  G_{i+1}: {pauli}")
    
    result1 = nested_adjoint(sequence1[:-1], sequence1[-1])
    if result1:
        print(f"\nResult: {result1}")
        print(f"Matches target: {result1 == target1}")
    else:
        print("\nThe sequence produced a zero commutator")
    # Example 2: XIII (a case where W=I)
    print("\n" + "=" * 50)
    print("EXAMPLE 2: Compiling XIII (W=I case)")
    print("=" * 50)
    
    target2 = get_pauli_string("XIII")
    sequence2 = pauli_compiler(target2, A_set, A_prime_set, B_set, B_prime_set, k, N)
    
    print("\nCompiled sequence:")
    for i, pauli in enumerate(sequence2):
        print(f"  G_{i+1}: {pauli}")
    
    result2 = nested_adjoint(sequence2[:-1], sequence2[-1])
    if result2:
        print(f"\nResult: {result2}")
        print(f"Matches target: {result2 == target2}")
    else:
        print("\nThe sequence produced a zero commutator")
    
    # Example 3: IIIZ (a case where V=I)
    print("\n" + "=" * 50)
    print("EXAMPLE 3: Compiling IIIZ (V=I case)")
    print("=" * 50)
    
    target3 = get_pauli_string("IIIZ")
    sequence3 = pauli_compiler(target3, A_set, A_prime_set, B_set, B_prime_set, k, N)
    
    print("\nCompiled sequence:")
    for i, pauli in enumerate(sequence3):
        print(f"  G_{i+1}: {pauli}")
    
    result3 = nested_adjoint(sequence3[:-1], sequence3[-1])
    if result3:
        print(f"\nResult: {result3}")
        print(f"Matches target: {result3 == target3}")
    else:
        print("\nThe sequence produced a zero commutator")
    
    # More examples to test robustness
    print("\n" + "=" * 50)
    print("EXAMPLE 4: Compiling XZYI (general case)")
    print("=" * 50)
    
    target4 = get_pauli_string("XZYI")
    sequence4 = pauli_compiler(target4, A_set, A_prime_set, B_set, B_prime_set, k, N)
    
    print("\nCompiled sequence:")
    for i, pauli in enumerate(sequence4):
        print(f"  G_{i+1}: {pauli}")
    
    result4 = nested_adjoint(sequence4[:-1], sequence4[-1])
    if result4:
        print(f"\nResult: {result4}")
        print(f"Matches target: {result4 == target4}")
    else:
        print("\nThe sequence produced a zero commutator")


if __name__ == "__main__":
    test_pauli_compiler()
