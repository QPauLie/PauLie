import re
import numpy as np


def get_so_basis(N: int) -> np.ndarray:
    dim = (N * (N - 1)) // 2
    basis = np.zeros((dim, N, N), dtype=float)
    idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            basis[idx, i, j] = 1.0
            basis[idx, j, i] = -1.0
            idx += 1
    return basis


def get_su_basis(N: int) -> np.ndarray:
    dim = N * N - 1
    basis = np.zeros((dim, N, N), dtype=complex)
    idx = 0
    for i in range(N):
        for j in range(i + 1, N):
            basis[idx, i, j] = 1.0j
            basis[idx, j, i] = 1.0j
            idx += 1
    for i in range(N):
        for j in range(i + 1, N):
            basis[idx, i, j] = -1.0
            basis[idx, j, i] = 1.0
            idx += 1
    for k in range(1, N):
        mat = np.zeros((N, N), dtype=complex)
        for i in range(k):
            mat[i, i] = 1.0j
        mat[k, k] = -k * 1.0j
        norm = np.sqrt(2.0 / (k * (k + 1)))
        basis[idx] = mat * norm
        idx += 1
    return basis


def get_sp_basis(N2: int) -> np.ndarray:
    N = N2 // 2
    dim = N * (2 * N + 1)
    basis = np.zeros((dim, N2, N2), dtype=float)
    idx = 0
    for i in range(N):
        for j in range(N):
            mat = np.zeros((N2, N2), dtype=float)
            mat[i, j] = 1.0
            mat[N + j, N + i] = -1.0
            basis[idx] = mat
            idx += 1
    for i in range(N):
        for j in range(i, N):
            mat = np.zeros((N2, N2), dtype=float)
            if i == j:
                mat[i, N + j] = 1.0
            else:
                mat[i, N + j] = 1.0
                mat[j, N + i] = 1.0
            basis[idx] = mat
            idx += 1
    for i in range(N):
        for j in range(i, N):
            mat = np.zeros((N2, N2), dtype=float)
            if i == j:
                mat[N + i, j] = 1.0
            else:
                mat[N + i, j] = 1.0
                mat[N + j, i] = 1.0
            basis[idx] = mat
            idx += 1
    return basis


def get_u_basis(N: int) -> np.ndarray:
    if N == 1:
        return np.array([[[1.0j]]], dtype=complex)
    su_b = get_su_basis(N)
    eye_mat = np.eye(N, dtype=complex) * 1.0j
    eye_mat = eye_mat / np.sqrt(N)
    return np.vstack([su_b, [eye_mat]])


def get_algebra_basis_from_label(
    label: str,
) -> list[np.ndarray]:
    n1 = 0
    clean_label = label.strip()
    if "*" in clean_label:
        parts = clean_label.split("*")
        if parts[0] == "2":
            n1 = 1
        clean_label = parts[-1]
    elif "^" in clean_label:
        match_n1 = re.match(r"2\^(\d+)", clean_label)
        if match_n1:
            n1 = int(match_n1.group(1))
        clean_label = clean_label.split()[-1]
    match_fam = re.match(r"(so|su|sp|u)\((\d+)\)", clean_label)
    if not match_fam:
        raise ValueError(
            f"Nieznany format algebry Liego: {label}"
        )
    family = match_fam.group(1)
    N = int(match_fam.group(2))
    num_summands = 2**n1
    if family == "so":
        single_basis = get_so_basis(N)
    elif family == "su":
        single_basis = get_su_basis(N)
    elif family == "sp":
        single_basis = get_sp_basis(N)
    elif family == "u":
        single_basis = get_u_basis(N)
    else:
        raise ValueError(
            f"Nieobslugiwana rodzina: {family}"
        )
    return [single_basis.copy() for _ in range(num_summands)]
