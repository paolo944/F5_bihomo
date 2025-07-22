// gauss_avx2.c
#include <immintrin.h>
#include <stdint.h>

void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords) {
    int i;
    int limit = (nwords / 4) * 4;

    for (i = 0; i < limit; i += 4) {
        __m256i d = _mm256_loadu_si256((__m256i*)(dest + i));
        __m256i s = _mm256_loadu_si256((__m256i*)(src + i));
        d = _mm256_xor_si256(d, s);
        _mm256_storeu_si256((__m256i*)(dest + i), d);
    }

    for (; i < nwords; i++) {
        dest[i] ^= src[i];
    }
}