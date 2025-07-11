class Mac:
    def __init__(self, n):
        self.matrix = matrix([])
        self.sig = []
        self.monomial_hash_list = {}
        self.monomial_inverse_search = []

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        self.poly_ring = PolynomialRing(GF(2), monomials_str, order='degrevlex')
        self.variables = [self.poly_ring(monom) for monom in monomials_str]
        return

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all the monomials of degree d in
        GF(2).<x1,..,xn> and add them to the hash
        list of the existant monomials
        """
        R = self.poly_ring
        n = R.ngens()
        subsets = Subsets(range(n), d)

        monomials = []
        for subset in subsets:
            exponents = [0]*n
            for i in subset:
                exponents[i] = 1
            monomials.append(R.monomial(*exponents))

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
        return vector(vec)

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
        new_row_matrix = Matrix(QQ, [vec])
        self.matrix = block_matrix([[self.matrix], [new_row_matrix]])
        return

    def F5_frobenius_criterion(self, u, f, i):
        for sig in self.sig:
            if 

        return

    def add_line(self, f, i, u):
        """
        add the line (u, f_i) to the Macaulay
        matrix
        """
        vec = self.polynomial_to_vector(u*f)


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

M = Mac(4)

M.monomial_ordered_list_deg_d(1)
M.monomial_ordered_list_deg_d(2)
M.monomial_ordered_list_deg_d(3)
M.monomial_ordered_list_deg_d(4)

print(M.monomial_hash_list)
print(M.monomial_inverse_search)

#def F5_criterion(sig, )

#def F5_bihom(F, d1, d2):

