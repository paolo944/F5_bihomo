# gauss_avx2.pyx

cimport cython
from libc.stdint cimport uint64_t
from libc.stdlib cimport malloc, free

cdef extern from "avx2_core.h":
    void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords)
    int first_pivot(uint64_t x)

@cython.boundscheck(False)
@cython.wraparound(False)
def gauss_elim_avx2(uint64_t[:, ::1] matrix not None, int nrows, int ncols_words):
    cdef int i, j, w, bit_pos
    cdef int pivot_word, pivot_bit
    cdef uint64_t pivot_mask, row_word
    cdef list zeros = []

    # Convert the matrix to raw pointer
    cdef uint64_t* mat_ptr = &matrix[0, 0]
    cdef int stride = matrix.strides[0] // sizeof(uint64_t)  # row stride in words
    cdef uint64_t* row_j
    cdef uint64_t* row_i

    for i in range(nrows):
        pivot_word = -1
        pivot_bit = -1

        # Get row pointer for current row
        row_i = mat_ptr + i * stride

        for w in range(ncols_words):
            row_word = row_i[w]
            if row_word != 0:
                pivot_word = w
                pivot_bit = first_pivot(row_word)
                break

        if pivot_word == -1 or pivot_bit == -1:
            continue

        pivot_mask = 1 << pivot_bit

        for j in range(i + 1, nrows):
            row_j = mat_ptr + j * stride

            if (row_j[pivot_word] & pivot_mask) != 0:
                xor_rows_avx2(row_j, row_i, ncols_words)

    return zeros
