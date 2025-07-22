import numpy as np

def compact_matrix(mat_bool):
    """
    Convertit une matrice booléenne en version compacte avec des bits.
    
    Args:
        mat_bool: matrice numpy booléenne
    
    Returns:
        matrice numpy uint64 compacte où chaque bit représente un booléen
    """
    nrows, ncols = mat_bool.shape
    ncols_words = (ncols + 63) // 64
    compact = np.zeros((nrows, ncols_words), dtype=np.uint64)
    
    for i in range(nrows):
        for j in range(ncols):
            if mat_bool[i, j]:
                word = j // 64
                bit = j % 64
                compact[i, word] |= np.uint64(1) << bit
                
    return compact

def compact_matrix_optimized(mat_bool):
    """
    Version optimisée utilisant la vectorisation numpy.
    """
    nrows, ncols = mat_bool.shape
    ncols_words = (ncols + 63) // 64
    compact = np.zeros((nrows, ncols_words), dtype=np.uint64)
    
    # Traiter par blocs de 64 colonnes
    for word_idx in range(ncols_words):
        start_col = word_idx * 64
        end_col = min(start_col + 64, ncols)
        
        # Extraire le bloc de colonnes
        block = mat_bool[:, start_col:end_col]
        
        # Convertir chaque colonne en bits
        for bit_pos in range(end_col - start_col):
            mask = block[:, bit_pos]
            compact[:, word_idx] |= (mask.astype(np.uint64) << bit_pos)
    
    return compact

def decomp_matrix(compact, original_shape):
    """
    Décompacte une matrice compacte vers sa forme booléenne originale.
    """
    nrows, ncols = original_shape
    mat_bool = np.zeros((nrows, ncols), dtype=bool)
    
    for i in range(nrows):
        for j in range(ncols):
            word = j // 64
            bit = j % 64
            if compact[i, word] & (np.uint64(1) << bit):
                mat_bool[i, j] = True
                
    return mat_bool

# Test de la fonction
if __name__ == "__main__":
    # Créer une matrice test
    test_matrix = np.array([
        [True, False, True, False, True],
        [False, True, False, True, False],
        [True, True, False, False, True]
    ], dtype=bool)
    
    print("Matrice originale:")
    print(test_matrix)
    
    # Tester la version corrigée
    compact1 = compact_matrix(test_matrix)
    print(f"\nMatrice compacte (version 1): shape={compact1.shape}")
    print(compact1)
    
    # Tester la version optimisée
    compact2 = compact_matrix_optimized(test_matrix)
    print(f"\nMatrice compacte (version 2): shape={compact2.shape}")
    print(compact2)
    
    # Vérifier que les deux versions donnent le même résultat
    print(f"\nLes deux versions sont identiques: {np.array_equal(compact1, compact2)}")
    
    # Décompacter pour vérifier
    decompacted = decomp_matrix(compact1, test_matrix.shape)
    print(f"\nDécompactage correct: {np.array_equal(test_matrix, decompacted)}")
    print("Matrice décompactée:")
    print(decompacted)