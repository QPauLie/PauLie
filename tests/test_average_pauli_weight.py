import numpy as np
from paulie.application.average_pauli_weight import quantum_fourier_entropy, avg_pauli_weights_from_strings, average_pauli_weight, get_pauli_weights
from matplotlib import pyplot as plt
from itertools import product
from tqdm import tqdm
import concurrent.futures


# --- Testing the Conjecture ---

def generate_hermitian_operator(num_qubits: int) -> np.ndarray:
    """
    Generates a random operator O such that O is Hermitian and O^2 = I.
    """
    dim = 2**num_qubits
    # 1. Create a random unitary matrix U (eigenvectors)
    rand_mat = np.random.rand(dim, dim) + 1j * np.random.rand(dim, dim)
    q, _ = np.linalg.qr(rand_mat)
    U = q

    # 2. Create a diagonal matrix D with eigenvalues +-1
    eigenvalues = np.random.choice([-1, 1], size=dim)
    D = np.diag(eigenvalues)

    # 3. Construct O = U D U_dagger
    O = U @ D @ U.conj().T
    return O

def compute_c_for_weights(args):
    n_qubits, weights = args
    O = generate_hermitian_operator(n_qubits)
    I = average_pauli_weight(O, weights=weights)
    H = quantum_fourier_entropy(O)
    if I < 1e-12:
        return None  # Mark as zero
    return H / I



def compute_c_for_operator(args):
    n_qubits, pauli_strings = args
    O = generate_hermitian_operator(n_qubits)
    I = avg_pauli_weights_from_strings(O, pauli_strings=pauli_strings)
    H = quantum_fourier_entropy(O)
    if I < 1e-12:
        return None  # Mark as zero
    return H / I


if __name__ == "__main__":
    print("--- Testing the Conjecture: H(O) <= c * I(O) ---")

    trial_count = 10000

    c_list_nqbits = {}
    max_ratio_nqbits = {}
    min_ratio_nqbits = {}
    number_of_zeros_nqbits = {}
    frequent_c_values = {}

    for n_qubits in range(1, 7):
        print(f"\nNumber of qubits: {n_qubits}")
        #pauli_strings = [''.join(p) for p in product('XYZI', repeat=n_qubits)]
        #args = [(n_qubits, pauli_strings)] * trial_count
        weights = get_pauli_weights(n_qubits, identity_pos=0)
        args = [(n_qubits, weights)] * trial_count

        c_list = []
        number_of_zeros = 0
        max_ratio = 0
        min_ratio = float('inf')

        with concurrent.futures.ProcessPoolExecutor() as executor:
            #results = list(tqdm(executor.map(compute_c_for_operator, args), total=trial_count))
            results = list(tqdm(executor.map(compute_c_for_weights, args), total=trial_count))

        for c in results:
            if c is None:
                number_of_zeros += 1
                continue
            c_list.append(c)
            if c > max_ratio:
                max_ratio = c
            if c < min_ratio:
                min_ratio = c
        c_list_nqbits[n_qubits] = c_list
        max_ratio_nqbits[n_qubits] = max_ratio
        number_of_zeros_nqbits[n_qubits] = number_of_zeros
        min_ratio_nqbits[n_qubits] = min_ratio
    # plot histogram of c values
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs = axs.flatten()
    for n_qubits in range(1, 7):
        res=axs[n_qubits-1].hist(c_list_nqbits[n_qubits], bins=100, alpha=0.7, color='blue')
        axs[n_qubits-1].set_title(f'Histogram of c values for {n_qubits} qubits')
        axs[n_qubits-1].set_xlabel('c = H/I')
        axs[n_qubits-1].set_ylabel('Frequency')
        # get x values of highest bin
        max_bin_index = np.argmax(res[0])
        max_bin_value = res[1][max_bin_index]
        frequent_c_values[n_qubits]=max_bin_value
        axs[n_qubits-1].axvline(x=max_bin_value, color='red', linestyle='--', label='Frequent Bin')
        axs[n_qubits-1].legend()  
    


    plt.tight_layout()
    plt.savefig("c_value_histograms.png", dpi=300)
    plt.show()
    for n_qubits in range(1, 7):
        max_ratio = max_ratio_nqbits[n_qubits]
        number_of_zeros = number_of_zeros_nqbits[n_qubits]
        min_ratio = min_ratio_nqbits[n_qubits]
        print(f"\nNumber of qubits: {n_qubits}")
        print(f"Most frequent c value: {frequent_c_values[n_qubits]:.4f}")
        print(f"Max ratio H/I: {max_ratio:.4f}")
        print(f"Min ratio H/I: {min_ratio:.4f}")
        print(f"Number of zeros I: {number_of_zeros}")

