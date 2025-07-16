import copy

class Mac:
    def __init__(self, n, d):
        self.matrix = None #The Macaulay matrix
        self.sig = [] #Array containing to know the index of the polynomial used in the line of this array's index
        self.monomial_hash_list = {} #To know which column index correspond to this monomial
        self.monomial_inverse_search = [] #To know which monomial correspond the column i

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        self.poly_ring = PolynomialRing(GF(2), monomials_str, order='degrevlex')
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
        Generate all the monomials of degree d in
        GF(2).<x1,..,xn> using Sage's monomials_of_degree(d),
        and add them to the hash list of the existing monomials.
        """
        R = self.quotient_ring
        poly_ring = self.poly_ring

        # Generate monomials of exact degree d using Sage built-in method
        monomials = poly_ring.monomials_of_degree(d)

        # Convert each monomial to the quotient ring
        monomials_in_R = [R(m) for m in monomials]

        # Sort them (if order matters for hashing)
        monomials_in_R.sort()

        # Add to hash list
        hash_size = len(self.monomial_hash_list)
        self.monomial_hash_list = {m: i + hash_size for i, m in enumerate(monomials_in_R)} | self.monomial_hash_list
        self.monomial_inverse_search += monomials_in_R
    

    def polynomial_to_vector(self, f):
        """
        Convert the polynomial f into a vector
        with columns corresponding to the monomial
        list ordering
        """
        n = len(self.monomial_hash_list)
        vec = [0] * n
        for monomial in f.monomials():
            try:
                index = self.monomial_hash_list[monomial]
                vec[n - 1 - index] = 1
            except Exception as error:
                print(f"Erreur polynomial_to_vector {error}")
                sys.exit()
        return vector(GF(2), vec)

    def row_lm(self, i):
        """
        Returns the lead monomial of the row i in the Macaulay
        matrix
        """
        row = self.matrix.row(i)
        for j, n in enumerate(row):
            if n != 0:
                break
        try:
            return self.monomial_inverse_search[len(self.monomial_inverse_search) - j - 1]
        except:
            print(f"{len(self.monomial_inverse_search) - j} / {len(self.monomial_inverse_search)}")

    def add_row(self, vec):
        if self.d < 4:
            new_row_matrix = matrix(GF(2), [vec])
        else:
            new_row_matrix = matrix(GF(2), [vec], sparse=True)
        self.matrix = self.matrix.stack(new_row_matrix)
        return

    def F5_frobenius_criterion(self, u, f, i, M):
        """
        Returns True if the row of signature (u, f_i)
        will reduce to 0 in the Macaulay martix
        """
        if M == None:
            return False
        for j, (_, f_index) in enumerate(M.sig):
            if f_index <= i:
                return M.row_lm(f_index) == u
            else:
                return False
        return False

    def add_line(self, f, i, u):
        """
        add the line (u, f_i) to the Macaulay
        matrix
        """
        vec = self.polynomial_to_vector(self.quotient_ring(u*f))
        if self.matrix == None:
            self.matrix = matrix(GF(2), vec)
        else:
            self.add_row(vec)
        self.sig.append((u, i))
        return

    def gauss(self):
        """
        Simple Gauss without pivoting
        """
        for i in range(self.matrix.nrows()):
            try:
                k = self.matrix.nonzero_positions_in_row(i)[0]
            except:
                continue

            for j in range(i+1, self.matrix.nrows()):
                try:
                    kp = self.matrix.nonzero_positions_in_row(j)[0]
                    if kp == k:
                        self.matrix.add_multiple_of_row(j, i, 1)
                except:
                    continue
        return

    def verify_reductions_zero(self):
        counter = 0
        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) == []:
                counter += 1
        return counter

    def add_lines(self, f_i, i, d, Mac_d_1, Mac_d_2):
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
                if x_i > x_lambda:
                    if f_i.total_degree() == 1:
                        if not self.F5_frobenius_criterion(x_i*e, f_i, i, Mac_d_1):
                            self.add_line(f_i, i, x_i*e)
                    elif f_i.total_degree() == 2:
                        if not self.F5_frobenius_criterion(x_i*e, f_i, i, Mac_d_2):
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
            already_in_gb = any(p.lm() == Md_lm for p in gb)
            if not already_in_gb:
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
    n = F[0].parent().ngens()
    Mac_d = None
    Mac_d_1 = None
    Mac_d_2 = None
    gb = []

    print(f"F5 for d={F[0].total_degree()}...{dmax}")

    for d in range(F[0].total_degree(), dmax+1):
        print(f"d={d}")
        Mac_d = Mac(n, d)
        Mac_d.monomial_ordered_list_deg_d(d)
        Mac_d.monomial_ordered_list_deg_d(d-1)
        Mac_d.monomial_ordered_list_deg_d(d-2)
        for i in range(0, m):
            f_i = F[i]
            if F[i].total_degree() == d:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, d, Mac_d_1, Mac_d_2)
        tmp_Mac = copy.deepcopy(Mac_d)
        Mac_d.gauss()
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree {d}: {reductions_to_zero} / {Mac_d.matrix.nrows()}")
        print(f"Corank of degree {d}: {Mac_d.corank()}")
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
    F = homogenized_ideal(doit(8, 9))
    D = Ideal(F).degree_of_semi_regularity()
    print(generating_bardet_series(F))
    print(f"degree of semi-regularity of F: {D}")

    gb = F5Matrix(F, D)

    gb = [lift(p) for p in gb]
    print(len(gb))

    print(Ideal(gb).basis_is_groebner())

    gb = Ideal(F).groebner_basis('msolve')

    print(len(gb))

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