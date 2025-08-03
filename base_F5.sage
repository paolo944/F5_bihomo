import copy
import time
import numpy as np
sys.path.insert(0, "./gauss_gf2")
import gauss_gf2
import pack_utils
from collections import defaultdict

def special_identity_matrix(nrows, ncols):
    m = []
    for i in range(min(nrows, ncols)):
        row = [0] * ncols
        row[i] = 1
        m.append(row)
    return m

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

class BaseMac:
    def __init__(self, F, monomials_str, estimated_total_size=None, estimated_nrows=None):
        self.matrix = None  # The Macaulay matrix
        self.sig = []  # Index of the polynomial used in each row

        self.crit = defaultdict(set)

        # For GF(2), we'll work directly with packed matrices
        if self.poly_ring.characteristic() == 2:
            self._ncols = None  # Will be set when we know the number of columns
            self._words_per_row = None
            self._current_row = 0
            self._estimated_nrows = estimated_nrows or 1000  # Default estimate
            
            if estimated_nrows:
                # Pre-allocate packed matrix buffer
                self._packed_buffer = None  # Will be allocated when ncols is known
                self._use_packed_buffer = True
            else:
                self._packed_rows = []  # List to store packed rows
                self._use_packed_buffer = False
        else:
            # For other fields, keep the original approach
            if estimated_total_size:
                self._buffer = np.empty(estimated_total_size, dtype=np.float64)
                self._current_size = 0
                self._use_buffer = True
            else:
                self._row_list = []
                self._use_buffer = False

    def _init_packed_matrix(self, ncols):
        """Initialize packed matrix structure when we know the number of columns"""
        if self.poly_ring.characteristic() == 2 and self._ncols is None:
            self._ncols = ncols
            self._words_per_row = (ncols + 63) // 64
            
            if self._use_packed_buffer:
                self._packed_buffer = np.zeros(
                    self._estimated_nrows * self._words_per_row, 
                    dtype=np.uint64
                )

    def polynomial_to_vector(self, f):
        """
        Convert polynomial f to vector according to monomial ordering
        Fonction testé -> Correcte
        """
        n = len(self.monomial_hash_list)
        if self.poly_ring.characteristic() == 2:
            vec = np.zeros(n, dtype=bool)
            for monomial, coeff in zip(f.monomials(), f.coefficients()):
                try:
                    index = self.monomial_hash_list[monomial]
                    vec[index] = coeff == 1
                except KeyError:
                    print(f"Erreur: monôme {monomial} non trouvé dans hash_list")
                    continue
            return vec
        else:
            vec = [0] * n
            for monomial, coeff in zip(f.monomials(), f.coefficients()):
                try:
                    index = self.monomial_hash_list[monomial]
                    vec[index] = coeff
                except KeyError:
                    print(f"Erreur: monôme {monomial} non trouvé dans hash_list")
                    continue
            return vec

    def row_lm(self, i):
        """
        Returns the leading monomial of row i
        Fonction testé -> Correcte
        """
        if self.poly_ring.characteristic() == 2:
            # For packed matrix, we need to check bits directly
            if hasattr(self, '_words_per_row') and self._words_per_row is not None:
                for w in range(self._words_per_row):
                    word = int(self.matrix[i * self._words_per_row + w])
                    if word != 0:
                        # Find first set bit in this word
                        for bit in range(64):
                            if (word & (1 << bit)) != 0:
                                j = w * 64 + bit
                                if j < self._ncols:
                                    return self.monomial_inverse_search[j]
                return 0
            else:
                # Fallback for unpacked matrix
                ncols = len(self.monomial_hash_list)
                for j in range(ncols):
                    if self.matrix[i * ncols + j]:
                        return self.monomial_inverse_search[j]
                return 0
        else:
            positions = self.matrix.nonzero_positions_in_row(i)
            if positions:
                return self.monomial_inverse_search[positions[0]]
            return 0
    
    def add_row(self, vec):
        """
        Add a row to the matrix - optimized for GF(2) to work directly with packed format
        Fonction testé -> Correcte
        """
        if self.poly_ring.characteristic() == 2:
            # Initialize packed matrix if not done yet
            if self._ncols is None:
                self._init_packed_matrix(len(vec))
            
            if self._use_packed_buffer:
                if self._current_row < self._estimated_nrows:
                    # Pack vector directly into the buffer
                    pack_utils.pack_vector_to_row(vec, self._packed_buffer, self._current_row, self._words_per_row)
                    self._current_row += 1
                else:
                    # Buffer overflow - fallback to list approach
                    packed_row = np.zeros(self._words_per_row, dtype=np.uint64)
                    pack_utils.pack_vector_to_row(vec, packed_row, 0, self._words_per_row)
                    
                    if not hasattr(self, '_packed_rows'):
                        # Move buffer contents to list
                        self._packed_rows = []
                        for row_idx in range(self._current_row):
                            row_start = row_idx * self._words_per_row
                            row_end = (row_idx + 1) * self._words_per_row
                            self._packed_rows.append(self._packed_buffer[row_start:row_end].copy())
                    
                    self._packed_rows.append(packed_row)
                    self._use_packed_buffer = False
            else:
                # Pack vector into a single row
                packed_row = np.zeros(self._words_per_row, dtype=np.uint64)
                pack_utils.pack_vector_to_row(vec, packed_row, 0, self._words_per_row)
                self._packed_rows.append(packed_row)
                self._current_row += 1
        else:
            # Original logic for non-GF(2) fields
            if self.d < 4:
                new_row_matrix = matrix(self.poly_ring.base_ring(), [vec])
            else:
                new_row_matrix = matrix(self.poly_ring.base_ring(), [vec], sparse=True)

            if self.matrix is None:
                self.matrix = new_row_matrix
            else:
                self.matrix = self.matrix.stack(new_row_matrix)
        return

    def finalize_matrix(self):
        """Finalize the matrix construction"""
        if self.poly_ring.characteristic() == 2:
            if self._use_packed_buffer:
                # Trim buffer to actual size
                actual_size = self._current_row * self._words_per_row
                self.matrix = self._packed_buffer[:actual_size]
            elif hasattr(self, '_packed_rows'):
                # Concatenate all packed rows
                self.matrix = np.concatenate(self._packed_rows)
            
            # Store matrix dimensions for later use
            self._nrows = self._current_row
            
            # Cleanup
            if hasattr(self, '_packed_buffer'):
                delattr(self, '_packed_buffer')
            if hasattr(self, '_packed_rows'):
                delattr(self, '_packed_rows')
        else:
            # Original logic for non-GF(2)
            if self._use_buffer:
                self.matrix = self._buffer[:self._current_size]
            elif hasattr(self, '_row_list'):
                self.matrix = np.concatenate(self._row_list)
            
            # Cleanup
            if hasattr(self, '_buffer'):
                delattr(self, '_buffer')
            if hasattr(self, '_row_list'):
                delattr(self, '_row_list')

    def vector_to_polynomial(self, i):
        """
        Convert row i of matrix back to polynomial
        Fonction testé -> Correcte
        """
        poly = 0
        for j in self.matrix.nonzero_positions_in_row(i):
            if j < len(self.monomial_inverse_search):
                poly += self.matrix[i, j] * self.monomial_inverse_search[j]
            else:
                print(f"Index {j} hors limites pour monomial_inverse_search")
        return poly

    def F5_criterion(self, M_prev):
        """
        F5 criterion: returns True if row with signature (u, f_i) 
        will reduce to 0
        """
        if M_prev is None:
            return

        if self.poly_ring.characteristic() == 2:
            nrows = M_prev._nrows if hasattr(M_prev, '_nrows') else len(M_prev.matrix) // M_prev._words_per_row
        else:
            nrows = M_prev.matrix.nrows()
        
        for j in range(nrows):
            _, sig_i = M_prev.sig[j]
            self.crit[sig_i].add(M_prev.row_lm(j))
        return

    def verify_reductions_zero(self):
        """
        Count zero rows after reduction - optimized for packed matrices
        Fonction testé -> Correcte
        """
        counter = 0
        lignes_0 = []
        
        if self.poly_ring.characteristic() == 2:
            if hasattr(self, '_words_per_row') and self._words_per_row is not None:
                # Working with packed matrix
                nrows = self._nrows if hasattr(self, '_nrows') else len(self.matrix) // self._words_per_row
                
                for i in range(nrows):
                    # Check if entire row is zero by checking all words
                    row_is_zero = True
                    for w in range(self._words_per_row):
                        word = int(self.matrix[i * self._words_per_row + w])
                        if word != 0:
                            row_is_zero = False
                            break
                    
                    if row_is_zero:
                        counter += 1
                        lignes_0.append(i)
            else:
                # Fallback for unpacked matrix
                ncols = len(self.monomial_hash_list)
                nrows = len(self.matrix) // ncols
                for i in range(nrows):
                    if not np.any(self.matrix[i * ncols : (i + 1) * ncols]):
                        counter += 1
                        lignes_0.append(i)
        else:
            for i in range(self.matrix.nrows()):
                if self.matrix.nonzero_positions_in_row(i) == []:
                    counter += 1
                    lignes_0.append(i)
        return counter, lignes_0

    def gauss(self, debug=True):
        """
        Simple Gauss sans pivot et sans backtracking avec l'élimination
        Now works directly with packed matrices for GF(2)
        """
        t0 = time.time()
        
        if self.poly_ring.characteristic() == 2:
            nrows = self._nrows if hasattr(self, '_nrows') else self._current_row
            ncols = self._ncols

            gauss_gf2.gaussian_elim(self.matrix, nrows, self._words_per_row)
        else:
            # Original logic for non-GF(2)
            self.matrix.echelonize()
 
        t1 = time.time()
        nrows = self._nrows if hasattr(self, '_nrows') else (self.matrix.nrows() if hasattr(self.matrix, 'nrows') else self._current_row)
        ncols = self._ncols if hasattr(self, '_ncols') else (self.matrix.ncols() if hasattr(self.matrix, 'ncols') else len(self.monomial_hash_list))
        print(f"[TIMER] Temps pour Gauss (matrice {nrows}x{ncols}) : {t1 - t0:.4f} s")

    def corank(self):
        """
        Compute corank of the matrix
        Fonction testé -> Correcte
        """
        nnz_columns = 0
        nnz_rows = 0

        ncols = len(self.monomial_hash_list)
        nrows = len(self.matrix) // ncols
        
        #for i in range(self.matrix.ncols()):
        #    if self.matrix.nonzero_positions_in_column(i) != []:
        #        nnz_columns += 1

        empty_row = True
        for i in range(nrows):
            for k in range(ncols):
                if self.matrix[i * ncols + k]:
                    empty_row = False
                    break
            if empty_row:
                nnz_rows += 1

        return ncols - nnz_rows

def update_gb(gb, Md, Mtilde):
    if Md.matrix.nrows() != Mtilde.matrix.nrows():
        print("Euuuh erreur pas normal, Mtilde.nrows() == Md.nrows() normalement")
    for i in range(Md.matrix.nrows()):
        Md_lm = Mtilde.row_lm(i)
        if Md.row_lm(i) != Md_lm:
            poly = Mtilde.vector_to_polynomial(i)
            already_in_gb = False
            for p in gb:
                try:
                    already_in_gb = already_in_gb or p.lm() == Md_lm
                except AttributeError:
                    if p == 1:
                        already_in_gb = already_in_gb or 1 == Md_lm
                    elif p == 0:
                        continue
                    else:
                        print("Oulalala erreur")
            if not already_in_gb:
                gb.append(poly)
    return

def generating_bardet_series(system):
    """
    Returns the generating series of bardet of a supposed 
    semi-regular system without taking into account the 
    field equations
    """
    series_ring.<z> = PowerSeriesRing(ZZ)
    n = system[0].parent().ngens()
    term1 = 1
    for i in system:
        term1 *= (1-z**i.degree())

    term2 = (1-z)**n
    return term1 / term2

def homogenized_ideal(system):
    """
    Returns the homogenized system that is supposed
    quadratic
    """
    system2 = []
    for i in system:
        try:
            system2.append(i.homogeneous_components()[2])
        except KeyError:
            system2.append(i.homogeneous_components()[1])

    return system2

def doit_bilinear(n_x, n_y, m):
    """
    Génère m polynômes bilinéaires homogènes de bi-degré (1,1)
    en n_x variables x_i et n_y variables y_j sur GF(2),
    avec une solution plantée qui annule tous les polynômes.
    """
    K = GF(2)

    x_vars = ['x{}'.format(i) for i in range(1, n_x + 1)]
    y_vars = ['y{}'.format(j) for j in range(1, n_y + 1)]
    R = PolynomialRing(K, x_vars + y_vars, order='degrevlex')

    # Solution plantée
    x_sol = vector(K, [K.random_element() for _ in range(n_x)])
    y_sol = vector(K, [K.random_element() for _ in range(n_y)])

    x_gens = vector(R, R.gens()[:n_x])
    y_gens = vector(R, R.gens()[n_x:])

    polynomials = []

    while len(polynomials) < m:
        A = random_matrix(K, n_x, n_y)
        f = x_gens * A * y_gens

        if f(*x_sol.list(), *y_sol.list()) == 0:
            polynomials.append(f)

    return polynomials

def doit(n, m):
    """
    Generate random system of n variables and m
    polynomials on GF(2)
    Stolen from hpXbred :)
    """
    # planted solution
    V = GF(2)**n 
    x = V.random_element() 
    I = []

    monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
    #monomials_str = ['x'+str(i) for i in range(1, n+1)]
    R = PolynomialRing(GF(2), monomials_str, order='degrevlex')
    
    def random_quad_poly(R):
        K = R.base_ring()
        v = vector(R.gens())
        n = len(v) 
        Mq = matrix.random(K, n, n)
        Ml = matrix.random(K, 1, n)
        f = v * Mq * v + (Ml*v)[0] + K.random_element()
        return f
    
    # m random polynomials
    for _ in range(m): 
        f = random_quad_poly(R) 
        f += f(*x) 
        I.append(f)

    return I
