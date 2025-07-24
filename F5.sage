import copy
import time
import numpy as np

class Mac:
    def __init__(self, n, d, F):
        self.matrix = None #The Macaulay matrix
        self.sig = [] #Array containing to know the index of the polynomial used in the line of this array's index
        self.monomial_hash_list = {} #To know which column index correspond to this monomial
        self.monomial_inverse_search = [] #To know which monomial correspond the column i

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        #monomials_str = ['x'+str(i) for i in range(1, n + 1)]
        self.poly_ring = PolynomialRing(F, monomials_str, order='degrevlex')
        variables = [self.poly_ring(monom) for monom in monomials_str]
        field_eq = [mon**2 + mon for mon in variables]
        #self.quotient_ring = self.poly_ring.quotient(field_eq)
        self.quotient_ring = self.poly_ring 
        self.variables = [self.quotient_ring(monom) for monom in variables]
        self.d = d
        return

    """
    def monomial_ordered_list_deg_d(self, d):
        
        Generate all the monomials of degree d in
        GF(2).<x1,..,xn> and add them to the hash
        list of the existant monomials

        R = self.quotient_ring
        n = R.ngens()
        subsets = Subsets(range(n), d)

        monomials = []
        for subset in subsets:
            exponents = [0]*n
            for i in subset:
                exponents[i] = 1
            monomials.append(R(self.poly_ring.monomial(*exponents)))

        monomials.sort()

        hash_size = len(self.monomial_hash_list)

        self.monomial_hash_list = {self.quotient_ring(m): i + hash_size for i, m in enumerate(monomials)} | self.monomial_hash_list
        self.monomial_inverse_search += [self.quotient_ring(m) for m in monomials]
        return
    """

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all monomials of degree d and add them to hash list
        Fonction testé -> Correcte
        """
        R = self.quotient_ring
        poly_ring = self.poly_ring

        monomials = poly_ring.monomials_of_degree(d)
        monomials_in_R = sorted([R(m) for m in monomials], reverse=True)

        hash_size = len(self.monomial_hash_list)
        self.monomial_hash_list = {m: i for i, m in enumerate(monomials_in_R)} | self.monomial_hash_list
        self.monomial_inverse_search += monomials_in_R
    
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
            if sig_i >= i:
                return False

            lm_j = M_prev.row_lm(j)

            try:
                if u % lm_j == 0:
                    return True
            except (TypeError, ZeroDivisionError):
                if lm_j == u:
                    return True
        return False

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        Fonction testé -> Correcte
        """
        vec = self.polynomial_to_vector(self.quotient_ring(u*f))
        self.add_row(vec)
        self.sig.append((u, i))
        return

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

    def add_lines(self, f_i, i, Mac_d_1, Mac_d_2):
        """
        Adds all the lines of signature (u, f_i)
        s.t. deg(uf_i) == d
        """
        #print(f"add lines for i:{i} f_i:{f_i}")
        for row_i, (e, f_ii) in enumerate(Mac_d_1.sig):
            if f_ii < i:
                continue
            elif f_ii > i:
                return
            if e == 1:
                x_lambda = 1
            else:
                x_lambda = e.monomials()[0].variables()[0] #biggest variable in e
                #print(f"{e.monomials()}")
            #print(f"e: {e}, x_lambda:{x_lambda}")
            for x_i in self.variables:
                if x_i >= x_lambda:
                    if f_i.total_degree() == 1:
                        #if not self.F5_criterion(x_i*e, f_i, i, Mac_d_1):
                        self.add_line(f_i, i, x_i*e)
                    elif f_i.total_degree() == 2:
                        #if not self.F5_criterion(x_i*e, f_i, i, Mac_d_2):
                        self.add_line(f_i, i, x_i*e)
                            #print(f"added ({x_i*e}, {f_i}) = {f_i*x_i*e}")
            #print()
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
    """
    Takes a list of polynomials that composes the Gröbner basis
    and updates it with Mac and Mac_tilde
    Fonction testé -> Correcte
    """
    if Md.matrix.nrows() != Mtilde.matrix.nrows():
        print("Euuuh erreur pas normal, Mtilde.nrows() == Md.nrows() normalement")
    #print(f"gb before: {gb}")
    for i in range(Md.matrix.nrows()):
        Md_lm = Mtilde.row_lm(i)
        if Md.row_lm(i) != Md_lm:
            poly = Mtilde.vector_to_polynomial(i)
            """
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
                print(f"add poly in gb: {poly}")
                gb.append(poly)
            """
            #print(f"add poly in gb: {poly}")
            gb.append(poly)
    return

def F5Matrix(F, dmax):
    """
    F is homogeneous polynomials ordered such that deg(f_i) < deg(f_j) forall i, j
    such that i < j
    dmax is the maximal degree to obtain a d_max-grobner basis of F
    F is supposed to be a quadratic system, so deg(f_i) <= 2 forall 0 <= i < m
    We follow F5Matrix by Bardet in her PhD thesis p.21
    """
    m = len(F)
    R = F[0].base_ring()
    n = F[0].parent().ngens()
    Mac_d = None
    Mac_d_1 = None
    Mac_d_2 = None
    gb = copy.copy(F)

    #print("Polynomials in F:")
    #for i in F:
    #    print(i)

    print(f"F5 for d={F[0].total_degree()}...{dmax}")

    for d in range(F[0].total_degree(), dmax+1):
        print(f"\n{'-' * 20} d={d} {'-' * 20}\n")
        t_deg_start = time.time()
        Mac_d = Mac(n, d, R)
        Mac_d.monomial_ordered_list_deg_d(d)
        t0 = time.time()
        for i in range(0, m):
            f_i = F[i]
            if F[i].total_degree() == d:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, Mac_d_1, Mac_d_2)
        t1 = time.time()
        print(f"[TIMER] Temps pour add_lines : {t1 - t0:.4f} s")
        #tmp_Mac = copy.deepcopy(Mac_d)
        Mac_d.matrix = matrix(GF(2), Mac_d.matrix)
        print("done conversion")
        Mac_d.gauss()
        #update_gb(gb, tmp_Mac, Mac_d)
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree {d} with normal Gauss: {reductions_to_zero} / {Mac_d.matrix.nrows()}")
        print(f"Corank of degree {d}: {Mac_d.corank()}")
        Mac_d_2 = Mac_d_1
        Mac_d_1 = Mac_d

        t_deg_end = time.time()
        print(f"[TIMER] Temps total pour degré {d} : {t_deg_end - t_deg_start:.2f} s")

        """
        #Verifier que les conversions sont bonnes:
        for (i, r) in enumerate(Mac_d.matrix.rows()):
            p = Mac_d.vector_to_polynomial(i)
            v = Mac_d.polynomial_to_vector(p)
            print(f"test ligne {i}: {v == r}")
        """

    return gb

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

if __name__ == '__main__':
    F = homogenized_ideal(doit(14, 15))
    """
    R.<x1, x2, x3> = PolynomialRing(GF(5), order='degrevlex')
    F = [x2^2 + 4*x2*x3,
    2*x1^2 + 3*x1*x2 + 4*x2^2 + 3*x3^2,
    3*x1^2 + 4*x1*x2 + 2*x2^2]
    """
    #F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_12_13.sobj"))
    #series_ring.<z> = PowerSeriesRing(ZZ)
    #hilbert_series = series_ring(Ideal(F).hilbert_series())
    #print(f"Hilbert Series: {hilbert_series}")
    D = Ideal(F).degree_of_semi_regularity()
    print(generating_bardet_series(F))
    #for i in F:
    #    print(i.total_degree())
    #print(f"degree of semi-regularity of F: {D}")

    gb = F5Matrix(F, D)

    #for i in F:
    #    print(i in gb)

    #gb = [lift(p) for p in gb]
    print(len(gb))
    #print(gb)

    gb2 = Ideal(gb).groebner_basis()

    print(Ideal(gb).basis_is_groebner())

    gb = Ideal(F).groebner_basis()
    gb3 = Ideal(F).groebner_basis()

    #print(gb)
    print(len(gb))
    print(len(gb3))
    print(len(gb2))

    """
    M = Mac(4)

    M.monomial_ordered_list_deg_d(1)
    M.monomial_ordered_list_deg_d(2)
    M.monomial_ordered_list_deg_d(3)
    M.monomial_ordered_list_deg_d(4)

    print(M.monomial_hash_list)
    print(M.monomial_inverse_search)
    print(M.quotient_ring)

    #def F5_criterion(sig, )

    #def F5_bihom(F, d1, d2):
    """