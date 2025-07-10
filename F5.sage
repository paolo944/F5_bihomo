class Mac:
    def __init__(self, n):
        self.matrix = matrix([])
        self.sig = []
        self.monomial_hash_list = {}

        monomials_str = ['x'+str(i) for i in range(1, n//2 + 1)] + ['y'+str(i) for i in range(1, n//2 + 1)]
        self.poly_ring = PolynomialRing(GF(2), monomials_str, order='degrevlex')
        self.variables = [self.poly_ring(monom) for monom in monomials_str]

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

    def add_line(self, f, i, u):
        """
        add the line (u, f_i) to the Macaulay
        matrix
        """
        return

M = Mac(4)

M.monomial_ordered_list_deg_d(1)
M.monomial_ordered_list_deg_d(2)
M.monomial_ordered_list_deg_d(3)
M.monomial_ordered_list_deg_d(4)

print(M.monomial_hash_list)

#def F5_criterion(sig, )

#def F5_bihom(F, d1, d2):

