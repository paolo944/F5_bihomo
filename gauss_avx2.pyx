# gauss_avx2.pyx
cimport cython
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free

cdef extern from "avx2_core.h":
    void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords)
    int ctzll(uint64_t x)

@cython.boundscheck(False)
@cython.wraparound(False)
def gauss_elim_avx2(uint64_t[:, ::1] matrix not None, int nrows, int ncols_words):
    cdef int i, j, w, bit_pos
    cdef int pivot_word, pivot_bit
    cdef uint64_t pivot_mask, row_word
    cdef list zeros = []
    
    print("Pivots for opti Gauss\n [", end="")    
    for i in range(nrows):
        pivot_word = -1
        pivot_bit = -1

        # Find first set bit in the row
        for w in range(ncols_words):
            row_word = matrix[i, w]
            if row_word != 0:
                bit_pos = ctzll(row_word)
                pivot_word = w
                pivot_bit = bit_pos
                break

        if pivot_word == -1:
            zeros.append(i)
            continue  # Entire row is zero

        # Compute global bit index (column index)
        print(pivot_word * 64 + pivot_bit, end=", ")
        pivot_mask = 1 << pivot_bit

        for j in range(i + 1, nrows):
            if (matrix[j, pivot_word] & pivot_mask) != 0:
                xor_rows_avx2(&matrix[j, 0], &matrix[i, 0], ncols_words)
    
    print("]")
    print("Zero rows:", zeros)
    return zeros  # Optional return if you want to use it in Python

