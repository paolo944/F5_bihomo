# cython: language_level=3

import numpy as np
cimport numpy as np
ctypedef np.uint64_t UINT64_t

cdef extern from "gauss_core.h":
    void gaussian_elim_packed(UINT64_t *matrix, int nrows, int words_per_row)
    void gaussian_elim_packed_minimal_optimized(UINT64_t *matrix, int nrows, int words_per_row)
    void gaussian_elim_packed_cache_optimized(UINT64_t *matrix, int nrows, int words_per_row)
    void gaussian_elim_packed_prefetch_optimized(UINT64_t *matrix, int nrows, int words_per_row)

def gaussian_elim(np.ndarray[UINT64_t, ndim=1, mode="c"] mat, nrows, words_per_row):
    if not mat.flags['C_CONTIGUOUS']:
        mat = np.ascontiguousarray(mat)

    cdef np.uint64_t[:] mat_view = mat  # typed memoryview for fast access

    # Get pointer to first element of contiguous matrix data
    cdef UINT64_t *ptr = &mat_view[0]

    gaussian_elim_packed_cache_optimized(ptr, nrows, words_per_row)
