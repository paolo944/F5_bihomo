import copy
import time

class BaseMac:
    def __init__(self, F, monomials_str):
        self.matrix = None  # The Macaulay matrix
        self.sig = []  # Index of the polynomial used in each row
        self.monomial_hash_list = {}  # Monomial -> column index
        self.monomial_inverse_search = []  # Column index -> monomial

        self.poly_ring = PolynomialRing(F, monomials_str, order='degrevlex')
        variables = [self.poly_ring(m) for m in monomials_str]
        self.quotient_ring = self.poly_ring  # No quotienting by field equations
        self.variables = [self.quotient_ring(v) for v in variables]

    def polynomial_to_vector(self, f):
        """
        Convert polynomial f to vector according to monomial ordering
        Fonction testé -> Correcte
        """
        n = len(self.monomial_hash_list)
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
        row = self.matrix.row(i)
        positions = self.matrix.nonzero_positions_in_row(i)
        if positions:
            return self.monomial_inverse_search[positions[0]]
        return None
    
    def add_row(self, vec):
        """
        Add a row to the matrix
        Fonction testé -> Correcte
        """
        if self.poly_ring.characteristic() == 2:
            if self.matrix == None:
                self.matrix = [vec]
                return
            else:
                self.matrix.append(vec)
                return
        else:
            if self.d < 4:
                new_row_matrix = matrix(self.poly_ring.base_ring(), [vec])
            else:
                new_row_matrix = matrix(self.poly_ring.base_ring(), [vec], sparse=True)

            if self.matrix is None:
                self.matrix = new_row_matrix
            else:
                self.matrix = self.matrix.stack(new_row_matrix)
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

    def F5_criterion(self, u, f_i, i, M_prev):
        """
        F5 criterion: returns True if row with signature (u, f_i) 
        will reduce to 0
        """
        if M_prev is None:
            return False
        
        # Selon Bardet: si u est un terme de tête dans M̃_{d-d_i,i-1}
        for j in range(M_prev.matrix.nrows()):
            sig_u, sig_i = M_prev.sig[j]
            if sig_i < i and M_prev.row_lm(j) == u:
                return True
            elif sig_i >= i:
                return False
        return False

    def verify_reductions_zero(self):
        """
        Count zero rows after reduction
        Fonction testé -> Correcte
        """
        counter = 0
        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) == []:
                counter += 1
        return counter

    def gauss(self):
        """
        Simple Gauss sans pivot et sans backtracking avec l'élimination
        Fonction testé -> Correcte
        """
        t0 = time.time()
        if self.poly_ring.characteristic() == 2:
            self.matrix.echelonize(algorithm="m4ri", reduced=False)

        else:
            self.matrix.echelonize()
 
        t1 = time.time()        
        print(f"[TIMER] Temps pour Gauss (matrice {self.matrix.nrows()}x{self.matrix.ncols()}) : {t1 - t0:.4f} s")

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

        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) != []:
                nnz_rows += 1

        return self.matrix.ncols() - nnz_rows

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