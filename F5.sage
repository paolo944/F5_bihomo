class Mac:
    def __init__(self, n):
        self.matrix = None #The Macaulay matrix
        self.sig = [] #Array containing to know the index of the polynomial used in the line of this array's index
        self.monomial_hash_list = {} #To know which column index correspond to this monomial
        self.monomial_inverse_search = [] #To know which monomial correspond the column i

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        self.poly_ring = PolynomialRing(GF(2), monomials_str, order='degrevlex')
        self.variables = [self.poly_ring(monom) for monom in monomials_str]
        field_eq = [mon**2 + mon for mon in self.variables]
        self.quotient_ring = self.poly_ring.quotient(field_eq)
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

        self.monomial_hash_list = {m: i + hash_size for i, m in enumerate(monomials)} | self.monomial_hash_list
        self.monomial_inverse_search += [m for m in monomials]
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
            except Exception as e:
                print(e)
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
        new_row_matrix = matrix(GF(2), [vec])
        self.matrix = block_matrix([ [self.matrix], [new_row_matrix] ])
        return

    def F5_frobenius_criterion(self, u, f, i, M):
        """
        Returns True if the row of signature (u, f_i)
        will reduce to 0 in the Macaulay martix
        """
        for j, (_, f_index) in enumerate(self.sig):
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
                print("Erreur, la ligne est nulle")
                sys.exit()

            for j in range(i+1, self.matrix.nrows()):
                try:
                    kp = self.matrix.nonzero_positions_in_row(j)[0]
                    if kp == k:
                        self.matrix.add_multiple_of_row(j, i, 1)
                except:
                    print("Erreur, la ligne est nulle")
                    sys.exit()
        return

    def verify_reductions_zero(self):
        counter = 0
        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) == []:
                counter += 1
        return counter

    def add_lines(self, f_i, i, d, Mac_d_1):
        """
        Adds all the lines of signature (u, f_i)
        s.t. deg(uf_i) == d
        """
        #for row_i, (e, f_ii) in enumerate(Mac_d_1.sig):
        return

def F5Matrix(F, dmax):
    """
    F is homogeneous polynomials ordered such that deg(f_i) < deg(f_j) forall i, j
    such that i < j
    dmax is the maximal degree to obtain a d_max-grobner basis of F
    F is supposed to be a quadratic system, so deg(f_i) <= 2 forall 0 <= i < m
    """
    m = len(F)
    n = F[0].parent().ngens()
    Mac_d = Mac(n)
    Mac_d_1 = None
    Mac_d_2 = None

    for d in range(F[0].total_degree(), dmax):
        for i in range(0, m):
            f_i = F[i]
            if f_i.total_degree() == d:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, d, Mac_d_1)
        Mac_d.gauss()
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree {d}: {reductions_to_zero}")
        Mac_d_1 = Mac_d
        Mac_d_2 = Mac_d_1

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

