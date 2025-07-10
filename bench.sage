# Test de performance
from time import time

# Matrice de base
M = random_matrix(QQ, 50000, 500)
new_row = [1]*500

### Méthode 1 : liste + reconstruction
start = time()
rows = M.rows()
rows.append(new_row)
M1 = Matrix(QQ, rows)
print("Liste + reconstruction :", time() - start)

### Méthode 2 : block_matrix
start = time()
new_row_matrix = Matrix(QQ, [new_row])
M2 = block_matrix([[M], [new_row_matrix]])
print("block_matrix :", time() - start)

### Méthode 3 : transpose + augment
start = time()
M3 = M.transpose().augment(Matrix(QQ, [[x] for x in new_row])).transpose()
print("transpose + augment :", time() - start)