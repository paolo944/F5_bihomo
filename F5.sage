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
        self.quotient_ring = self.poly_ring.quotient(field_eq)
        self.variables = [self.quotient_ring(monom) for monom in variables]
        self.d = d
        return

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all the monomials of degree d in
        GF(2).<x1,..,xn> and add them to the hash
        list of the existant monomials
        """
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
        for i, n in enumerate(row):
            if n != 0:
                break
        return self.monomial_inverse_search[len(self.monomial_inverse_search) - i]

    def add_row(self, vec):
        if self.d < 4:
            new_row_matrix = matrix(GF(2), [vec])
        else:
            new_row_matrix = matrix(GF(2), [vec], sparse=True)
        self.matrix = block_matrix([ [self.matrix], [new_row_matrix] ])
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
        print(self.matrix)
        for i in range(self.matrix.nrows()):
            try:
                k = self.matrix.nonzero_positions_in_row(i)[0]
            except:
                print("Erreur, la ligne est nulle")
                #sys.exit()

            for j in range(i+1, self.matrix.nrows()):
                try:
                    kp = self.matrix.nonzero_positions_in_row(j)[0]
                    if kp == k:
                        self.matrix.add_multiple_of_row(j, i, 1)
                except:
                    print("Erreur, la ligne est nulle")
                    #sys.exit()
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
            if e == 1:
                x_lambda = 1
            else:
                print(self.poly_ring(e))
                x_lambda = self.quotient_ring(list(self.poly_ring(e).variables()).sort()[0]) #biggest variable in e
            for x_i in self.variables:
                if x_i > x_lambda:
                    if f_i.total_degree() == 1:
                        if not self.F5_frobenius_criterion(x_i*e, f_i, i, Mac_d_1):
                            self.add_line(f_i, i, x_i*e)
                    elif f_i.total_degree() == 2:
                        if not self.F5_frobenius_criterion(x_i*e, f_i, i, Mac_d_2):
                            self.add_line(f_i, i, x_i*e)
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

    print(f"F5 for d={F[0].total_degree()}...{dmax}")

    for d in range(F[0].total_degree(), dmax+1):
        print(f"d={d}")
        Mac_d = Mac(n, d)
        Mac_d.monomial_ordered_list_deg_d(d)
        Mac_d.monomial_ordered_list_deg_d(d-1)
        Mac_d.monomial_ordered_list_deg_d(d-2)
        print(Mac_d)
        print(Mac_d_1)
        print(Mac_d_2)
        for i in range(0, m):
            f_i = F[i]
            if F[i].total_degree() == d:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, d, Mac_d_1, Mac_d_2)
        Mac_d.gauss()
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree {d}: {reductions_to_zero}")
        Mac_d_2 = Mac_d_1
        Mac_d_1 = Mac_d
    return

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

if __name__ == '__main__':
    F = homogenized_ideal(doit(4, 5))

    F5Matrix(F, 5)

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