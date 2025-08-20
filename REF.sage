import cProfile
sys.path.insert(0, "./gauss_gf2")
import gauss_gf2
import numpy as np
import time

def pack_matrix(M):
    nrows, ncols = M.nrows(), M.ncols()
    words_per_row = (ncols + 63) // 64  # ceil division
    packed = np.zeros(nrows * words_per_row, dtype=np.uint64)
    
    for i in range(nrows):
        for j in range(ncols):
            if M[i,j] != 0:
                word_idx = j // 64
                bit_idx = j % 64
                packed[i*words_per_row + word_idx] |= np.uint64(1) << np.uint64(bit_idx)
    
    return packed, words_per_row

def unpack_matrix(packed, nrows, ncols):
    words_per_row = (ncols + 63) // 64
    M_new = Matrix(GF(2), nrows, ncols)
    
    for i in range(nrows):
        for j in range(ncols):
            word_idx = j // 64
            bit_idx = j % 64
            # Conversion explicite en int pour éviter les problèmes de types
            word_value = int(packed[i*words_per_row + word_idx])
            if (word_value >> bit_idx) & 1:
                M_new[i,j] = 1
    
    return M_new

def REF(M):
    R_piv, C_piv, R_nonpiv, C_nonpiv = gauss_gf2.find_pivots_cython(M)

    A = M[R_piv, C_piv]
    B = M[R_piv, C_nonpiv]
    C_block = M[R_nonpiv, C_piv]
    D = M[R_nonpiv, C_nonpiv]

    # Transformations
    B = A.inverse() * B
    D = D - C_block * B

    t1 = time.time()
    packed, words_per_row = gauss_gf2.pack_matrix(D)
    gauss_gf2.gaussian_elim(packed, D.nrows(), words_per_row)
    D = gauss_gf2.unpack_matrix(packed, D.nrows(), D.ncols())
    t2 = time.time()
    print(f"Temps de calcul: {t2 - t1}s")

    # Reconstruct the matrix
    M_reconstructed = Matrix(GF(2), M.nrows(), M.ncols(), 0)
    A_id = identity_matrix(GF(2), len(R_piv))

    for i, r in enumerate(R_piv):
        for j, c in enumerate(C_piv):
            M_reconstructed[r, c] = A_id[i, j]

    for i, r in enumerate(R_piv):
        for j, c in enumerate(C_nonpiv):
            M_reconstructed[r, c] = B[i, j]

    for i, r in enumerate(R_nonpiv):
        for j, c in enumerate(C_nonpiv):
            M_reconstructed[r, c] = D[i, j]

    return M_reconstructed

M = load("matrix_test_3.sobj")
cProfile.run("REF(M)", sort='cumtime')