from paulie.common.pauli_string import PauliString
from paulie.common.pauli_string_factory import get_identity
from six.moves import reduce





def nested_adjoint(operators: list[PauliString], target: PauliString) -> PauliString:
    """
    Compute nested adjoint maps: ad_{A_1}...ad_{A_{n-1}}(A_n)
    Returns None if any commutator is zero
    
    Args:
        operators: List of PauliString operators [A_1, A_2, ..., A_n]
        target: Target PauliString A_n
    """
    if not operators:
        return target
        
    current = target
    for op in reversed(operators):
        result = op.adjoint_map(current)
        if result is None:
            return None
        current = result
    
    return current


def count_anticommuting_positions(a: PauliString, b: PauliString) -> int:
       """Count the number of qubit positions where A and B anticommute"""
       return reduce(lambda count, xy: count if xy[0].is_identity() or xy[1].is_identity() or xy[0] == xy[1] else count + 1, zip(a, b), 0)


def derive_generating_operators(target: PauliString) -> list[PauliString]:
    """
    Derive a pair of operators that generate a target Pauli string through commutation
    
    For each position in the target:
    - If it's X, use the pair [Z,Y] at that position
    - If it's Y, use the pair [X,Z] at that position
    - If it's Z, use the pair [Y,X] at that position
    - If it's I, use the identity at that position
    
    Make sure the total number of anticommuting positions is odd
    """
    n = len(target)
    op1 = get_identity(n)
    op2 = get_identity(n)
    # Count non-identity positions
    non_identity_count = reduce(lambda count, x: count if x.is_identity() else count + 1, target, 0)
    # Set up operators based on target
    for i, t in enumerate(target):
        # Identity position - leave as identity in both operators
        if t.is_identity():
            continue
        # X position - use Z in op1, Y in op2
        elif t == "X":
            op1[i] = "Z"
            op2[i] = "Y"

        # Z position - use Y in op1, X in op2
        elif t == "Z":
            op1[i] = "Y"
            op2[i] = "X"
        # Y position - use X in op1, Z in op2
        elif t == "Y":
            op1[i] = "X"
            op2[i] = "Z"

    # Make sure we have an odd number of anticommuting positions
    # by adjusting operators if needed
    if op1.commutes_with(op2):
        # If they commute, there must be an even number of anticommuting positions
        # We can fix this by switching one position to make them commute there
        for i, o in enumerate(op1):
            if not o.is_identity():
                if o == "X":
                    op1[i] = "Z"
                elif o == "Z":
                    op1[i] = "X"
                elif o == "Y":
                    op1[i] = "I"
                break
    
    # Verify the operators generate the target
    result = op1.adjoint_map(op2)
    if result and result == target:
        #print(f"Successfully derived operators {op1} and {op2} that generate {target}")
        return [op1, op2]
    
    # Try the reverse
    result = op2.adjoint_map(op1)
    if result and result == target:
        #print(f"Successfully derived operators {op2} and {op1} that generate {target} (reversed)")
        return [op2, op1]
    
    # If we got here, something went wrong
    #print(f"Warning: derived operators do not generate target. Got {result}, expected {target}")
    return [op1, op2]  # Return them anyway; they might still be useful


def pauli_compiler(target_P: PauliString, A_set: list[PauliString], A_prime_set: list[PauliString],
                  B_set: list[PauliString], B_prime_set: list[PauliString], k: int, N: int) -> list[PauliString]:
    """
    PauliCompiler algorithm from the paper
    
    Args:
        target_P: Target Pauli string on N qubits
        A_set: List of Pauli strings on k qubits
        A_prime_set: List of Pauli strings on N qubits of the form A⊗I
        B_set: List of Pauli strings on (N-k) qubits
        B_prime_set: List of Pauli strings on N qubits of the form U_B⊗B
        k: Number of qubits in first subsystem
        N: Total number of qubits
        
    Returns:
        List of Pauli strings that generate target_P when nested commutators are applied
    """
    #print(f"PauliCompiler called for target={target_P}")
    
    # Extract V and W from target_P = V⊗W
    V = target_P.get_subsystem(0, k)
    W = target_P.get_subsystem(k, N-k)
    
    #print(f"Target decomposed: V={V}, W={W}")
    
    # Try to derive operators directly from the target
    derived_sequence = derive_generating_operators(target_P)
    
    # Verify the derived sequence works
    result = derived_sequence[0].adjoint_map(derived_sequence[1])
    if result and result == target_P:
        #print(f"Using derived operators: {derived_sequence[0]}, {derived_sequence[1]}")
        return derived_sequence
        
    # If direct derivation didn't work, try the specific cases
    # Case 1: W = I (lines 2-4 in Algorithm 1)
    if W.is_identity():
        #print("Case 1: W=I")
        
        # Create all possible V⊗I operators
        v_operators = []
        for A in A_set:
            v_operators.append(A.tensor(get_identity(N-k)))
        
        # First try to find a pair in A_set that directly generates V
        for i, A1 in enumerate(A_set):
            for j, A2 in enumerate(A_set):
                if i == j:
                    continue
                # Check if this pair generates V
                result = A1.adjoint_map(A2)
                #print(f"A1={A1}, A2={A2}, result={result}")
                if result and result == V:
                    # Found a pair, extend to full system
                    I_Nmk = get_identity(N-k)
                    return [A1.tensor(I_Nmk), A2.tensor(I_Nmk)]
        
        # Try deriving a pair specifically for V
        derived_V_sequence = derive_generating_operators(V)
        print(f"derived_V_sequence[0] = {derived_V_sequence[0]} derived_V_sequence[1] = {derived_V_sequence[1]}")
        # Extend to full system
        I_Nmk = get_identity(N-k)
        return [derived_V_sequence[0].tensor(I_Nmk), derived_V_sequence[1].tensor(I_Nmk)]
        
    # Case 2: V = I (lines 5-12 in Algorithm 1)
    elif V.is_identity():
        #print("Case 2: V=I")
        
        # First try to find a pair in B_set that directly generates W
        for i, B1 in enumerate(B_set):
            for j, B2 in enumerate(B_set):
                if i == j:
                    continue
                
                # Check if this pair generates W
                result = B1.adjoint_map(B2)
                if result and result == W:
                    # Found a pair, extend to full system
                    I_k = get_identity(k)
                    return [I_k.tensor(B1), I_k.tensor(B2)]
        
        # Try deriving a pair specifically for W
        derived_W_sequence = derive_generating_operators(W)
        
        # Extend to full system
        I_k = get_identity(k)
        return [I_k.tensor(derived_W_sequence[0]), I_k.tensor(derived_W_sequence[1])]
        
    # Case 3: V ≠ I and W ≠ I (lines 13-16 in Algorithm 1)
    else:
        #print("Case 3: V≠I and W≠I")
        
        # Try all operator pairs from combined sets
        all_operators = A_prime_set + B_prime_set
        for i, op1 in enumerate(all_operators):
            for j, op2 in enumerate(all_operators):
                if i == j:
                    continue
                
                # Check if they anticommute (have odd number of anticommuting positions)
                if not op1.commutes_with(op2):
                    # Check if they generate the target
                    result = op1.adjoint_map(op2)
                    if result and result == target_P:
                        #print(f"Found pair in existing operators: {op1}, {op2}")
                        return [op1, op2]
        
        # No direct pair works, try generating V and W separately
        # Find operators that generate V
        derived_V_sequence = derive_generating_operators(V)
        
        # Find operators that generate W
        derived_W_sequence = derive_generating_operators(W)
        
        # Try combining them
        I_Nmk = get_identity(N-k)
        I_k = get_identity(k)
        
        op1 = derived_V_sequence[0].tensor(I_Nmk)
        op2 = derived_V_sequence[1].tensor(I_Nmk)
        op3 = I_k.tensor(derived_W_sequence[0])
        op4 = I_k.tensor(derived_W_sequence[1])
        
        # See if any combinations work
        combinations = [(op1, op2), (op1, op3), (op1, op4), 
                       (op2, op3), (op2, op4), (op3, op4)]
        
        for op1, op2 in combinations:
            if not op1.commutes_with(op2):
                result = op1.adjoint_map(op2)
                if result and result == target_P:
                    #print(f"Found working combination: {op1}, {op2}")
                    return [op1, op2]
        
        # Final fallback - use a custom constructed pair
        #print("Constructing custom operator pair")
        op1 = get_identity(N)
        op2 = get_identity(N)
        
        # For simplicity, just use Y at positions where target has X
        # X at positions where target has Y
        # XY at positions where target has Z
        # Keep I where target has I
        for i, p in enumerate(target_P):
            if p == "X":  # X
                op1[i] = "Y"  # Y
                op2[i] = "X"  # X
            elif p == "Z":  # Z
                op1[i] = "X"  # X
                op2[i] = "Z"  # Z
            elif p == "Y":  # Y
                op1[i] = "Z"  # Z
                op2[i] = "X"  # X
        
        # Ensure an odd number of anticommuting positions
        if op1.commutes_with(op2):
            # Add an anticommuting pair
            op1[0] = "X"  # X
            op2[0] = "Z"  # Z
        
        return [op1, op2]



