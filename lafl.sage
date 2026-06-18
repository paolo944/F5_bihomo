import time
import os
import psutil

def find_pivots(M):
    used_cols = set()
    C_piv = []
    R_piv = []

    current_row = -1
    for i, j in M.nonzero_positions():
        if i == current_row:
            continue

        if j not in used_cols:
            C_piv.append(j)
            R_piv.append(i)
            used_cols.add(j)
            current_row = i
        else:
            current_row = i

    sorted_pairs = sorted(zip(C_piv, R_piv))
    C_piv = [c for c, r in sorted_pairs]
    R_piv = [r for c, r in sorted_pairs]

    return C_piv, R_piv

def construct_blocks(M, C_piv, R_piv, C_npiv, R_npiv):
    A = M[R_piv, C_piv]
    B = M[R_piv, C_npiv]
    C = M[R_npiv, C_piv]
    D = M[R_npiv, C_npiv]
    return A, B, C, D

def Trsm(A, B, density_threshold=0.2):
    b_density = B.density()
    
    if b_density >= density_threshold:
        B_prepared = B.dense_matrix()
    else:
        B_prepared = B.sparse_matrix()

    return A.solve_right(B)

def Schur(B, C, D):
    return D - C*B

M = load("Matrices/test_matrix_7.sobj")
print(f"Matrice creuse : {M.is_sparse()}")

t_start = time.time()

C_piv, R_piv = find_pivots(M)

all_rows = set(range(M.nrows()))
all_cols = set(range(M.ncols()))

R_npiv = sorted(list(all_rows - set(R_piv)))
C_npiv = sorted(list(all_cols - set(C_piv)))

print(f"Nombre de pivots (colonnes) : {len(C_piv)}")
print(f"Nombre de pivots (lignes) : {len(R_piv)}")

A, B, C, D = construct_blocks(M, C_piv, R_piv, C_npiv, R_npiv)

B = Trsm(A, B)
print(f"densité A: {float(A.density())}, densité B: {float(B.density())}, densité C: {float(C.density())}, densité D: {float(D.density())}") 
D = Schur(B, C, D)
print(f"Nouvelle densité D: {float(D.density())}")
print(f"Rank of M = {len(R_piv) + D.rank()}")

t_end = time.time()

print(f"Temps: {t_end - t_start:.4f}")

process = psutil.Process(os.getpid())
mem_bytes = process.memory_info().rss
mem_mb = mem_bytes / (1024 ** 2)
mem_gb = mem_bytes / (1024 ** 3)
if mem_gb >= 1.0:
    print(f"Mémoire RAM max utilisée : {float(mem_gb):.2f} Go")
else:
    print(f"Mémoire RAM max utilisée : {float(mem_mb):.2f} Mo")