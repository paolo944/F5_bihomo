import numpy as np
import time
from statistics import mean
sys.path.insert(0, "./gauss_gf2")
import gauss_gf2

def pack_matrix(matrix, nrows, ncols):
    words_per_row = (ncols + 63) // 64
    packed = np.zeros(nrows * words_per_row, dtype=np.uint64)
    for i in range(nrows):
        row = matrix[i * ncols : (i + 1) * ncols]
        for w in range(words_per_row):
            chunk = row[w * 64 : (w + 1) * 64]
            word = 0
            for bit in range(len(chunk)):
                word |= int(chunk[bit]) << bit
            packed[i * words_per_row + w : (i + 1) * words_per_row + w] = word
    return (packed, words_per_row)
    
def unpack_matrix(packed, nrows, ncols, words_per_row):
    matrix = np.zeros(nrows * ncols, dtype=bool)
    for i in range(nrows):
        for w in range(words_per_row):
            word = int(packed[i * words_per_row + w])
            for bit in range(64):
                idx = w * 64 + bit
                if idx < ncols:
                    matrix[i * ncols + idx] = (word >> bit) & 1
    return matrix

def xor_gauss_elim(matrix, nrows, ncols):
    for i in range(nrows - 1):
        for k in range(ncols):
            if matrix[i * ncols + k]:
                break
        for j in range(i + 1, nrows):
            if matrix[j * ncols + k]:
                start_j = j * ncols
                start_i = i * ncols
                matrix[start_j : start_j + ncols] = np.logical_xor(matrix[start_i : start_i + ncols], matrix[start_j : start_j + ncols])


def time_run(func, *args):
    start = time.time()
    func(*args)
    return time.time() - start


def run_benchmarks(sizes, runs_per_size=5, seed=42):
    print(f"{'Rows':>6} {'Cols':>6} | {'XOR Mean':>10}  {'AVX Mean':>10}  | {'XOR Best':>10}  {'AVX Best':>10}  | {'XOR Worst':>10}  {'AVX Worst':>10} | {'Rapport':>10}")
    print("-" * 90)
    for nrows, ncols in sizes:
        xor_times = []
        avx_times = []

        for _ in range(runs_per_size):
            np.random.seed(seed)
            flat = np.random.randint(0, 2, size=nrows * ncols, dtype=np.uint8)

            # NumPy XOR
            m1 = flat.copy()
            xor_times.append(time_run(xor_gauss_elim, m1, nrows, ncols))

            # AVX version
            m2 = flat.copy()
            start = time.time()
            m2 = gauss_gf2.gaussian_elim(m2, nrows, ncols)
            avx_times.append(time.time() - start)

            # Optional: check correctness
            assert np.array_equal(m1, m2), "Mismatch between XOR and AVX"

        print(f"{nrows:6d} {ncols:6d} | "
              f"{mean(xor_times):10.6f}  {mean(avx_times):10.6f}  | "
              f"{min(xor_times):10.6f}  {min(avx_times):10.6f}  | "
              f"{max(xor_times):10.6f}  {max(avx_times):10.6f} | {mean(xor_times) / mean(avx_times) :10.6f}")


if __name__ == "__main__":
    # Customize sizes and repetitions here
    matrix_sizes = [
        (128, 512),
        (256, 1024),
        (512, 2048),
        (1024, 4096),
        (2058, 4096),
        (4096, 4096),
        (8192, 4096),
    ]
    run_benchmarks(matrix_sizes, runs_per_size=5)


"""
First version
  Rows   Cols |   XOR Mean    AVX Mean  |   XOR Best    AVX Best  |  XOR Worst   AVX Worst |    Rapport
------------------------------------------------------------------------------------------
   128    512 |   0.020190    0.000054  |   0.019893    0.000052  |   0.020519    0.000059 | 371.094654
   256   1024 |   0.085047    0.000245  |   0.083915    0.000221  |   0.086931    0.000319 | 346.592888
   512   2048 |   0.365872    0.001149  |   0.361377    0.001113  |   0.368110    0.001180 | 318.350635
  1024   4096 |   1.577385    0.007052  |   1.561237    0.006802  |   1.612986    0.007231 | 223.668704
  2058   4096 |   6.684518    0.030219  |   6.400812    0.028579  |   6.838416    0.031475 | 221.200865
  4096   4096 |  29.730632    0.163382  |  25.932879    0.124252  |  34.164297    0.210109 | 181.970372
"""