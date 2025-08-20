# cython: language_level=3

import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint64_t

# Import SageMath nécessaires (en mode Python, pas cimport)
from sage.matrix.constructor import Matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.matrix_mod2_dense cimport Matrix_mod2_dense
from sage.libs.m4ri cimport mzd_t

ctypedef np.uint64_t UINT64_t
ctypedef np.uint8_t BOOL_t

cdef extern from "gauss_core.h":
    void gaussian_elim_packed_cache_optimized(UINT64_t *matrix, int nrows, int words_per_row)

@cython.boundscheck(False)
@cython.wraparound(False)
def gaussian_elim(np.ndarray[UINT64_t, ndim=1, mode="c"] mat, nrows, words_per_row):
    if not mat.flags['C_CONTIGUOUS']:
        mat = np.ascontiguousarray(mat)

    cdef np.uint64_t[:] mat_view = mat  # typed memoryview for fast access

    # Get pointer to first element of contiguous matrix data
    cdef UINT64_t *ptr = &mat_view[0]

    gaussian_elim_packed_cache_optimized(ptr, nrows, words_per_row)

@cython.boundscheck(False)
@cython.wraparound(False)
def pack_matrix_sage(M):
    """Version optimisée Cython pour pack_matrix"""
    cdef int nrows = M.nrows()
    cdef int ncols = M.ncols()
    cdef int words_per_row = (ncols + 63) // 64
    
    cdef np.ndarray[UINT64_t, ndim=1] packed = np.zeros(nrows * words_per_row, dtype=np.uint64)
    
    cdef int i, j, word_idx, bit_idx
    cdef uint64_t bit_mask
    cdef uint64_t *packed_ptr = <uint64_t*>np.PyArray_DATA(packed)
    
    for i in range(nrows):
        for j in range(ncols):
            # Utiliser l'accès standard SageMath
            if M[i, j] != 0:
                word_idx = j >> 6  # Division par 64 optimisée
                bit_idx = j & 63   # Modulo 64 optimisé
                bit_mask = <uint64_t>1 << bit_idx
                packed_ptr[i * words_per_row + word_idx] |= bit_mask
    
    return packed, words_per_row

@cython.boundscheck(False)
@cython.wraparound(False)
def unpack_matrix_sage(np.ndarray[UINT64_t, ndim=1] packed, int nrows, int ncols):
    """Version optimisée Cython pour unpack_matrix"""
    cdef int words_per_row = (ncols + 63) // 64
    
    # Créer la matrice SageMath
    M_new = Matrix(GF(2), nrows, ncols)
    
    cdef int i, j, word_idx, bit_idx
    cdef uint64_t word_value, bit_mask
    cdef uint64_t *packed_ptr = <uint64_t*>np.PyArray_DATA(packed)
    
    for i in range(nrows):
        for j in range(ncols):
            word_idx = j >> 6  # Division par 64 optimisée
            bit_idx = j & 63   # Modulo 64 optimisé
            word_value = packed_ptr[i * words_per_row + word_idx]
            bit_mask = <uint64_t>1 << bit_idx
            
            if word_value & bit_mask:
                M_new[i, j] = 1  # Utilisation de l'accès standard
    
    return M_new

@cython.boundscheck(False)
@cython.wraparound(False)
def pack_matrix(np.ndarray[BOOL_t, ndim=2] M):
    """
    Pack une matrice booléenne NumPy (nrows x ncols) 
    en un tableau de uint64 (bitset par ligne).
    """
    cdef int nrows = M.shape[0]
    cdef int ncols = M.shape[1]
    cdef int words_per_row = (ncols + 63) // 64
    
    cdef np.ndarray[UINT64_t, ndim=1] packed = np.zeros(nrows * words_per_row, dtype=np.uint64)
    cdef uint64_t *packed_ptr = <uint64_t*>np.PyArray_DATA(packed)

    cdef int i, j, word_idx, bit_idx
    cdef uint64_t bit_mask
    
    for i in range(nrows):
        for j in range(ncols):
            if M[i, j] != 0:
                word_idx = j >> 6
                bit_idx = j & 63
                bit_mask = <uint64_t>1 << bit_idx
                packed_ptr[i * words_per_row + word_idx] |= bit_mask
    
    return packed, words_per_row

@cython.boundscheck(False)
@cython.wraparound(False)
def unpack_matrix(np.ndarray[UINT64_t, ndim=1] packed, int nrows, int ncols):
    """
    Décompresse un tableau packé en une matrice booléenne NumPy.
    """
    cdef int words_per_row = (ncols + 63) // 64
    cdef np.ndarray[BOOL_t, ndim=2] M_new = np.zeros((nrows, ncols), dtype=np.bool_)
    
    cdef uint64_t *packed_ptr = <uint64_t*>np.PyArray_DATA(packed)
    cdef int i, j, word_idx, bit_idx
    cdef uint64_t word_value, bit_mask
    
    for i in range(nrows):
        for j in range(ncols):
            word_idx = j >> 6
            bit_idx = j & 63
            word_value = packed_ptr[i * words_per_row + word_idx]
            bit_mask = <uint64_t>1 << bit_idx
            if word_value & bit_mask:
                M_new[i, j] = 1  # bool = True
    
    return M_new

@cython.boundscheck(False)
@cython.wraparound(False)
def find_pivots_cython(M):
    """Version Cython ultra-rapide pour find_pivots"""
    cdef int nrows = M.nrows()
    cdef int ncols = M.ncols()
    
    # Listes pour les résultats
    cdef list C_piv = []
    cdef list R_piv = []
    cdef list R_nonpiv = []
    cdef list C_nonpiv = []
    
    # Array pour marquer les colonnes utilisées
    cdef np.ndarray[np.uint8_t, ndim=1] col_used = np.zeros(ncols, dtype=np.uint8)
    cdef np.ndarray[np.uint8_t, ndim=1] row_used = np.zeros(nrows, dtype=np.uint8)
    
    cdef int i, j, lead_col
    
    # Première passe : trouver les pivots
    for i in range(nrows):
        lead_col = -1
        
        # Recherche du premier élément non-nul
        for j in range(ncols):
            if M[i, j] != 0:
                lead_col = j
                break
        
        if lead_col >= 0 and col_used[lead_col] == 0:
            C_piv.append(lead_col)
            R_piv.append(i)
            col_used[lead_col] = 1
            row_used[i] = 1
    
    # Construire les listes des non-pivots
    for i in range(nrows):
        if row_used[i] == 0:
            R_nonpiv.append(i)
    
    for j in range(ncols):
        if col_used[j] == 0:
            C_nonpiv.append(j)
    
    return R_piv, C_piv, R_nonpiv, C_nonpiv