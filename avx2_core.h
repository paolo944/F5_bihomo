// gauss_avx2.h
#ifndef GAUSS_AVX2_H
#define GAUSS_AVX2_H

#include <stdint.h>

void xor_rows_avx2(uint64_t* dest, const uint64_t* src, int nwords);
int first_pivot(uint64_t x);

#endif