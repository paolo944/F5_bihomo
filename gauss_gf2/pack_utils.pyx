# pack_utils.pyx - Cython implementation
import numpy as np
cimport numpy as cnp
cimport cython
from libc.stdint cimport uint64_t, uint8_t

ctypedef cnp.uint64_t DTYPE_t
ctypedef cnp.uint8_t BOOL_t

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline uint64_t pack_64_bits(BOOL_t[:] vec, int start_idx, int max_idx) nogil:
    """Pack up to 64 bits from a boolean array into a uint64"""
    cdef uint64_t word = 0
    cdef int i, bit_pos
    cdef int end_idx
    
    if start_idx + 64 < max_idx:
        end_idx = start_idx + 64
    else:
        end_idx = max_idx
    
    for i in range(start_idx, end_idx):
        if vec[i]:
            bit_pos = i - start_idx
            word |= <uint64_t>1 << bit_pos
    
    return word

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline uint64_t pack_64_bits_2d(BOOL_t[:, :] vec, int row, int start_idx, int max_idx) nogil:
    """Pack up to 64 bits from a 2D boolean array into a uint64"""
    cdef uint64_t word = 0
    cdef int i, bit_pos
    cdef int end_idx
    cdef int max_bits = 64
    
    if start_idx + 64 < max_idx:
        end_idx = start_idx + 64
    else:
        end_idx = max_idx
        max_bits = max_idx - start_idx
    
    for i in range(max_bits):
        if start_idx + i < max_idx and vec[row, start_idx + i]:
            word |= <uint64_t>1 << i
    
    return word

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pack_vector_to_row(cnp.ndarray[BOOL_t, ndim=1] vec, 
                           cnp.ndarray[DTYPE_t, ndim=1] packed_matrix,
                           int row_idx, 
                           int words_per_row):
    """
    Optimized Cython version of pack_vector_to_row
    """
    cdef int ncols = vec.shape[0]
    cdef int w, start_idx
    cdef uint64_t word
    cdef int base_offset = row_idx * words_per_row
    
    cdef BOOL_t[:] vec_view = vec
    cdef DTYPE_t[:] packed_view = packed_matrix
    
    with nogil:
        for w in range(words_per_row):
            start_idx = w * 64
            if start_idx < ncols:
                word = pack_64_bits(vec_view, start_idx, ncols)
                packed_view[base_offset + w] = word
            else:
                packed_view[base_offset + w] = 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def pack_multiple_vectors_fast(cnp.ndarray[BOOL_t, ndim=2] vectors,
                              cnp.ndarray[DTYPE_t, ndim=1] packed_matrix,
                              int start_row,
                              int words_per_row):
    """
    Pack multiple vectors at once for even better performance
    """
    cdef int nrows = vectors.shape[0]
    cdef int ncols = vectors.shape[1]
    cdef int row, w, start_idx
    cdef uint64_t word
    cdef int base_offset
    
    cdef BOOL_t[:, :] vec_view = vectors
    cdef DTYPE_t[:] packed_view = packed_matrix
    
    with nogil:
        for row in range(nrows):
            base_offset = (start_row + row) * words_per_row
            for w in range(words_per_row):
                start_idx = w * 64
                if start_idx < ncols:
                    word = pack_64_bits_2d(vec_view, row, start_idx, ncols)
                    packed_view[base_offset + w] = word
                else:
                    packed_view[base_offset + w] = 0

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def bi_degree(monomial, int nx, int ny):
    """
    Calcule le bidegré d’un monôme (optimisé Cython).
    Renvoie (d1, d2) où :
      - d1 = somme des exposants sur les nx premières variables
      - d2 = somme des exposants sur les ny suivantes
    """
    cdef object exponents
    cdef int i, d1 = 0, d2 = 0

    # monomial.exponents() renvoie un tuple de tuples, on prend le premier
    exponents = monomial.exponents()[0]

    for i in range(nx):
        d1 += exponents[i]

    for i in range(nx, nx + ny):
        d2 += exponents[i]

    return (d1, d2)
