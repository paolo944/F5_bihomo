# gauss_avx2.pyx
cimport cython
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free

cdef extern from "avx2_core.h":
    void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords)

@cython.boundscheck(False)
@cython.wraparound(False)
def gauss_elim_avx2(uint64_t[:, ::1] matrix not None, int nrows, int ncols_words):
    cdef int i, j, w, pos
    cdef int pivot_word, pivot_bit
    cdef uint64_t pivot_mask, row_word

    print("Pivots for opti Gauss")    
    for i in range(nrows):
        pivot_word = -1
        pivot_bit = -1
        for w in range(ncols_words):
            row_word = matrix[i, w]
            if row_word != 0:
                pos = 0
                while ((row_word >> pos) & 1) == 0:
                    pos += 1
                pivot_word = w
                pivot_bit = pos
                break

        print(pivot_word * 64 + pivot_bit, end=" ")
        if pivot_word == -1:
            continue  # zero row

        pivot_mask = 1 << pivot_bit

        for j in range(i + 1, nrows):
            if (matrix[j, pivot_word] & pivot_mask) != 0:
                # XOR row i into row j
                xor_rows_avx2(&matrix[j, 0], &matrix[i, 0], ncols_words)