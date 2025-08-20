#include <immintrin.h>
#include "gauss_core.h"

static inline void xor_rows_avx2(uint64_t *restrict a, uint64_t *restrict b, int words) {
    int i = 0;
    
    // Process 8 words (512 bits) at a time using two AVX2 operations
    for (; i + 8 <= words; i += 8) {
        // Prefetch next cache lines
        _mm_prefetch((char*)(a + i + 8), _MM_HINT_T0);
        _mm_prefetch((char*)(b + i + 8), _MM_HINT_T0);
        
        __m256i va1 = _mm256_loadu_si256((__m256i *)(a + i));
        __m256i vb1 = _mm256_loadu_si256((__m256i *)(b + i));
        __m256i va2 = _mm256_loadu_si256((__m256i *)(a + i + 4));
        __m256i vb2 = _mm256_loadu_si256((__m256i *)(b + i + 4));
        
        __m256i res1 = _mm256_xor_si256(va1, vb1);
        __m256i res2 = _mm256_xor_si256(va2, vb2);
        
        _mm256_storeu_si256((__m256i *)(b + i), res1);
        _mm256_storeu_si256((__m256i *)(b + i + 4), res2);
    }
    
    // Handle remaining 4-word chunks
    for (; i + 4 <= words; i += 4) {
        __m256i va = _mm256_loadu_si256((__m256i *)(a + i));
        __m256i vb = _mm256_loadu_si256((__m256i *)(b + i));
        __m256i res = _mm256_xor_si256(va, vb);
        _mm256_storeu_si256((__m256i *)(b + i), res);
    }
    
    // Handle remaining words
    for (; i < words; i++) {
        b[i] ^= a[i];
    }
}

static inline int find_pivot_bit_fast(uint64_t *row, int words_per_row) {
    for (int w = 0; w < words_per_row; ++w) {
        if (row[w]) {
            return __builtin_ctzll(row[w]) + 64 * w;
        }
    }
    return -1;
}

// Cache-optimized version that maintains identical algorithm behavior
void gaussian_elim_packed_cache_optimized(uint64_t *matrix, int nrows, int words_per_row) {
    for (int i = 0; i < nrows - 1; ++i) {
        uint64_t * __restrict pivot = matrix + i * words_per_row;
        int pivot_bit = -1;
        
        // EXACT same pivot finding logic as original
        for (int w = 0; w < words_per_row; ++w) {
            if (pivot[w]) {
                pivot_bit = __builtin_ctzl(pivot[w]) + 64 * w;
                break;
            }
        }
        
        if (pivot_bit == -1) continue;
        
        int pivot_word = pivot_bit / 64;
        uint64_t pivot_mask = 1ULL << (pivot_bit % 64);
        
        // Prefetch the pivot row to improve cache hit rate
        for (int w = 0; w < words_per_row; w += 8) {
            _mm_prefetch((char*)(pivot + w), _MM_HINT_T0);
        }
        
        // Process rows in blocks to improve cache locality while maintaining order
        const int CACHE_BLOCK_SIZE = 16;
        
        for (int block_start = i + 1; block_start < nrows; block_start += CACHE_BLOCK_SIZE) {
            int block_end = (block_start + CACHE_BLOCK_SIZE < nrows) ? 
                           block_start + CACHE_BLOCK_SIZE : nrows;
            
            // Process each row in the block sequentially (maintains identical behavior)
            for (int j = block_start; j < block_end; ++j) {
                uint64_t * __restrict row = matrix + j * words_per_row;
                
                // Prefetch next row while processing current one
                if (j + 1 < block_end) {
                    uint64_t *next_row = matrix + (j + 1) * words_per_row;
                    _mm_prefetch((char*)next_row, _MM_HINT_T0);
                }
                
                if (row[pivot_word] & pivot_mask) {
                    xor_rows_avx2(pivot, row, words_per_row);
                }
            }
        }
    }
}

