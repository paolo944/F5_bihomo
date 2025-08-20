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
    def __init__(self, F, monomials_str, h):
        self.sig = []  # Index of the polynomial used in each row

        self.crit = defaultdict(set)

        self.n = len(self.monomial_hash_list)
        self.is_char2 = (self.poly_ring.characteristic() == 2)
        self.ncols = len(self.monomial_inverse_search)
        self.nrows = self.ncols - h

        self.matrix = np.zeros((self.nrows, self.ncols), dtype=bool) # The Macaulay matrix
        self.matrix_index = 0

        if self.is_char2:
            self.words_per_row = None

    def polynomial_to_vector(self, f):
        """
        Convert polynomial f to vector according to monomial ordering
        Fonction testé -> Correcte
        """
        n = len(self.monomial_hash_list)
        monomial_hash = self.monomial_hash_list
        if self.is_char2:
            vec = np.zeros(n, dtype=bool)
            for monomial, coeff in zip(f.monomials(), f.coefficients()):
                index = monomial_hash.get(monomial)
                if index is not None:
                    vec[index] = (coeff == 1)
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
        if self.is_char2:
            position = np.argmax(self.matrix[i])
        else:
            position = np.argmax(self.matrix[i] != 0)
        if self.matrix[i, position] or self.matrix[i, position] != 0:
            return self.monomial_inverse_search[position]
        return 0
    
    def add_row(self, vec):
        """
        Add a row to the matrix - optimized for GF(2) to work directly with packed format
        Fonction testé -> Correcte
        """
        if self.matrix_index < self.nrows:
            self.matrix[self.matrix_index, :] = vec
        else:
            self.matrix = np.vstack([self.matrix, vec])
            #print("erreur, pas la bonne taille de matrice")    
        self.matrix_index += 1
        return

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
        
        for j in range(M_prev.nrows):
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
        for i in range(self.nrows):
            if np.all(self.matrix[i] == 0):
                counter += 1
                lignes_0.append(i)
        return counter, lignes_0

    def gauss(self, debug=True):
        """
        Simple Gauss sans pivot et sans backtracking avec l'élimination
        Now works directly with packed matrices for GF(2)
        """
        if self.is_char2:
            packed, words_per_row = gauss_gf2.pack_matrix(self.matrix)
            gauss_gf2.gaussian_elim(packed, self.nrows, words_per_row)
            self.matrix = gauss_gf2.unpack_matrix(packed, self.nrows, self.ncols)
        else:
            # Original logic for non-GF(2)
            self.matrix.echelonize()

    def corank(self):
        """
        Compute corank of the matrix
        Fonction testé -> Correcte
        """
        nnz_columns = 0
        nnz_rows = 0
        
        #for i in range(self.matrix.ncols()):
        #    if self.matrix.nonzero_positions_in_column(i) != []:
        #        nnz_columns += 1

        empty_row = True
        for i in range(self.nrows):
            for k in range(self.ncols):
                if self.matrix[i * self.ncols + k]:
                    empty_row = False
                    break
            if empty_row:
                nnz_rows += 1

        return self.ncols - nnz_rows

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
