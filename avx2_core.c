#include <stdint.h>
#include <immintrin.h>

void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords) {
    int i = 0;
    int nw = nwords;

    for (; i + 4 <= nw; i += 4) {
        __m256i d = _mm256_loadu_si256((__m256i*)(dest + i));
        __m256i s = _mm256_loadu_si256((__m256i*)(src + i));
        d = _mm256_xor_si256(d, s);
        _mm256_storeu_si256((__m256i*)(dest + i), d);
    }

    for (; i < nw; ++i) {
        dest[i] ^= src[i];
    }
}

int first_pivot(uint64_t x) {
    uint64_t mask;
    for (int i = 0; i < 64; i++) {
        mask = 1 << i;
        if ((x & mask) != 0)
            return i;
    }
    return -1;
}