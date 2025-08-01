#ifndef GAUSS_CORE_H
#define GAUSS_CORE_H

#include <stdint.h>

void gaussian_elim_packed(uint64_t *matrix, int nrows, int words_per_row);

void gaussian_elim_packed_minimal_optimized(uint64_t *matrix, int nrows, int words_per_row);

void gaussian_elim_packed_cache_optimized(uint64_t *matrix, int nrows, int words_per_row);

void gaussian_elim_packed_prefetch_optimized(uint64_t *matrix, int nrows, int words_per_row);

#endif