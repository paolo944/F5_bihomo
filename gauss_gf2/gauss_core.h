#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>


void gaussian_elim_packed_cache_optimized(uint64_t *matrix, int nrows, int words_per_row);

#ifdef __cplusplus
}
#endif