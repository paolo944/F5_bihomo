import copy

class Mac:
    def __init__(self, n, d1, d2):
        self.matrix = None #The Macaulay matrix
        self.sig = [] #Array containing to know the index of the polynomial used in the line of this array's index
        self.monomial_hash_list = {} #To know which column index correspond to this monomial
        self.monomial_inverse_search = [] #To know which monomial correspond the column i

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        self.poly_ring = PolynomialRing(GF(2), monomials_str, order='degrevlex')
        variables = [self.poly_ring(monom) for monom in monomials_str]
        #field_eq = [mon**2 + mon for mon in variables]
        #self.quotient_ring = self.poly_ring.quotient(field_eq)
        #self.quotient_ring = self.poly_ring 
        self.variables = [self.poly_ring(monom) for monom in variables]
        self.d = (d1, d2)
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

    def monomial_ordered_list_deg_d(self, d1, d2):
        """
        Generate all the monomials of degree d in
        GF(2).<x1,..,xn> using Sage's monomials_of_degree(d),
        and add them to the hash list of the existing monomials.
        """
        R = self.poly_ring
        n = R.ngens()
        monomials = []

        for exponents_x in IntegerVectors(d1, n//2):
            for exponents_y in IntegerVectors(d2, n//2):
                monomials.append(R.monomial(*(list(exponents_x) + list(exponents_y))))

        monomials.sort()
        hash_size = len(self.monomial_hash_list)
        self.monomial_hash_list = {m: i + hash_size for i, m in enumerate(monomials)} | self.monomial_hash_list
        self.monomial_inverse_search += [m for m in monomials]
        return
    

    def polynomial_to_vector(self, f):
        """
        Convert polynomial f to vector according to monomial ordering
        """
        n = len(self.monomial_hash_list)
        vec = [0] * n
        for monomial in f.monomials():
            try:
                index = self.monomial_hash_list[monomial]
                vec[index] = 1
            except KeyError:
                print(f"Erreur: monôme {monomial} non trouvé dans hash_list")
                continue
        return vector(GF(2), vec)

    def row_lm(self, i):
        """
        Returns the leading monomial of row i
        """
        row = self.matrix.row(i)
        positions = self.matrix.nonzero_positions_in_row(i)
        if positions:
            return self.monomial_inverse_search[positions[0]]
        return None

    def add_row(self, vec):
        """
        Add a row to the matrix
        """
        d1, d2 = self.d
        if d1 + d2 < 4:
            new_row_matrix = matrix(GF(2), [vec])
        else:
            new_row_matrix = matrix(GF(2), [vec], sparse=True)
        
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
            if sig_i < i and M_prev.row_lm(j) == u:
                return True
        return False

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        """
        vec = self.polynomial_to_vector(u*f)
        self.add_row(vec)
        self.sig.append((u, i))
        return

    def gauss(self):
        """
        Gaussian elimination following Bardet's description
        Only operations allowed: row_i ← c * row_i + c' * row_{i-j} with j > 0
        """
        nrows = self.matrix.nrows()
        for i in range(nrows):
            try:
                lead_i = self.matrix.nonzero_positions_in_row(i)[0]
            except IndexError:
                continue  # ligne nulle
            
            # Réduire les lignes suivantes par la ligne i
            for j in range(i+1, nrows):
                try:
                    positions_j = self.matrix.nonzero_positions_in_row(j)
                    if positions_j and positions_j[0] == lead_i:
                        self.matrix.add_multiple_of_row(j, i, 1)  # row_j += row_i
                except IndexError:
                    continue

    def verify_reductions_zero(self):
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
        for row_i, (e, f_ii) in enumerate(Mac_d_1.sig):
            if f_ii < i:
                continue
            elif f_ii > i:
                return
            if e == 1:
                x_lambda = 1
            else:
                x_lambda = e.monomials()[0].variables()[0] #biggest variable in e
            for x_i in self.variables:
                if x_i >= x_lambda:
                    if f_i.total_degree() == 1:
                        if not self.F5_criterion(x_i*e, f_i, i, Mac_d_1):
                            self.add_line(f_i, i, x_i*e)
                    elif f_i.total_degree() == 2:
                        if not self.F5_criterion(x_i*e, f_i, i, Mac_d_2):
                            self.add_line(f_i, i, x_i*e)
        return
    
    def vector_to_polynomial(self, i):
        """
        Convert the vector of row i of the Macaulay matrix
        into the polynomial correspondong to it
        """        
        n = len(self.monomial_inverse_search)
        poly = 0
        for j in self.matrix.nonzero_positions_in_row(i):
            try:
                poly += self.monomial_inverse_search[len(self.monomial_inverse_search) - j - 1]
            except:
                print(f"{len(self.monomial_inverse_search) - j} / {len(self.monomial_inverse_search)}")
                return
        return poly

    def corank(self):
        nnz_columns = 0
        nnz_rows = 0
        for i in range(self.matrix.ncols()):
            if self.matrix.nonzero_positions_in_column(i) != []:
                nnz_columns += 1

        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) != []:
                nnz_rows += 1

        return nnz_columns - nnz_rows

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

def integer_vectors_at_most(max_sum):
    """Génère tous les vecteurs d'entiers de longueur donnée avec somme ≤ max_sum"""
    result = []
    for s in range(max_sum + 1):
        for v in IntegerVectors(s, 2, min_part=1):
            result.append(v)
    return result

def F5Matrix(F, dmax):
    """
    F is homogeneous polynomials ordered such that deg(f_i) < deg(f_j) forall i, j
    such that i < j
    dmax is the maximal degree to obtain a d_max-grobner basis of F
    F is supposed to be a quadratic system, so deg(f_i) <= 2 forall 0 <= i < m
    We follow F5Matrix by Bardet in her PhD thesis p.21
    """
    m = len(F)
    n = F[0].parent().ngens()
    Mac_d = None
    Mac_d_1 = None
    Mac_d_2 = None
    gb = []

    degs = integer_vectors_at_most(dmax)

    for i in degs:
        print(i)

    print(f"F5 bihomogeneous for d={F[0].total_degree()}...{dmax}")

    for (d1, d2) in degs:
        print(f"d=({d1}, {d2})")
        Mac_d = Mac(n, d1, d2)
        Mac_d.monomial_ordered_list_deg_d(d1, d2)
        #Mac_d.monomial_ordered_list_deg_d(d-1)
        #Mac_d.monomial_ordered_list_deg_d(d-2)
        for i in range(0, m):
            f_i = F[i]
            if d1+d2 == 2:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, Mac_d_1, Mac_d_2)
        tmp_Mac = copy.deepcopy(Mac_d)
        Mac_d.gauss()
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree ({d1}, {d2}): {reductions_to_zero} / {Mac_d.matrix.nrows()}")
        print(f"Corank of degree ({d1}, {d2}): {Mac_d.corank()}")
        update_gb(gb, tmp_Mac, Mac_d)
        Mac_d_2 = Mac_d_1
        Mac_d_1 = Mac_d
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
    """
    test = integer_vectors_at_most(5)
    for i in test:
        print(i)
    """
    #F = homogenized_ideal(doit(8, 9))
    F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_8_9.sobj"))
    D = Ideal(F).degree_of_semi_regularity()
    print(generating_bardet_series(F))
    for i in F:
        print(i.total_degree())
    print(f"degree of semi-regularity of F: {D}")

    gb = F5Matrix(F, 3)

    """
    print(len(gb))

    print(Ideal(gb).basis_is_groebner())

    gb = Ideal(F).groebner_basis('msolve')

    print(len(gb))

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