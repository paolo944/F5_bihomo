import copy
import time

def bi_degree(monomial, nx, ny):
    exponents = monomial.exponents()[0]
    d1 = sum(exponents[:nx])
    d2 = sum(exponents[nx:nx+ny])
    return (d1, d2)

class Mac:
    def __init__(self, d1, d2, F, nx, ny):
        self.matrix = None #The Macaulay matrix
        self.sig = [] #Array containing to know the index of the polynomial used in the line of this array's index
        self.monomial_hash_list = {} #To know which column index correspond to this monomial
        self.monomial_inverse_search = [] #To know which monomial correspond the column i
        self.nx = nx
        self.ny = ny

        monomials_str = ['x'+str(i) for i in range(1, nx + 1)] + ['y'+str(i) for i in range(1, ny + 1)]
        #monomials_str = ['x'+str(i) for i in range(1, n + 1)]
        self.poly_ring = PolynomialRing(F, monomials_str, order='degrevlex')
        variables = [self.poly_ring(monom) for monom in monomials_str]
        #field_eq = [mon**2 + mon for mon in variables]
        #self.quotient_ring = self.poly_ring.quotient(field_eq)
        self.quotient_ring = self.poly_ring 
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

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all monomials of degree d and add them to hash list
        Fonction testé -> Correcte
        """
        R = self.quotient_ring
        poly_ring = self.poly_ring

        monomials = poly_ring.monomials_of_degree(d)
        monomials_in_R = []
        for monomial in monomials:
            if bi_degree(monomial, self.nx, self.ny) == self.d:
                monomials_in_R.append(R(monomial))
        
        monomials_in_R = sorted(monomials_in_R, reverse=True)

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
        return vector(f.base_ring(), vec)

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
        Fonction testé -> Correcte
        """
        if self.poly_ring.characteristic() == 2:
            #print("char 2")
            if self.matrix == None:
                self.matrix = matrix(GF(2), 1, len(vec), [vec], sparse=False)
                return
            else:
                self.matrix = self.matrix.transpose().augment(vec).transpose()
                return
        else:
            #print("not char 2")
            if self.d[0] + self.d[1] < 4:
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
            if sig_i < i and M_prev.row_lm(j) == u:
                return True
            elif sig_i >= i:
                return False
        return False

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        """
        poly = u*f
        if bi_degree(poly, self.nx, self.ny) != self.d:
            #print("not bi-degree")
            return
        vec = self.polynomial_to_vector(u*f)
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
            nrows = self.matrix.nrows()
            for i in range(nrows):
                try:
                    lead_i = self.matrix.nonzero_positions_in_row(i)[0]
                except IndexError:
                    continue  # ligne nulle
                for j in range(i+1, nrows):
                    try:
                        positions_j = self.matrix.nonzero_positions_in_row(j)
                        if lead_i in positions_j:
                            factor = self.matrix[j, positions_j.index(lead_i)] / self.matrix[i, lead_i]
                            self.matrix.add_multiple_of_row(j, i, factor)
                    except IndexError:
                        continue
 
        t1 = time.time()        
        print(f"[TIMER] Temps pour Gauss (matrice {self.matrix.nrows()}x{self.matrix.ncols()}) : {t1 - t0:.4f} s")

    def verify_reductions_zero(self):
        counter = 0
        for i in range(self.matrix.nrows()):
            if self.matrix.nonzero_positions_in_row(i) == []:
                counter += 1
        return counter
    
    def biggest_x(self, monomial):
        if monomial == 1 or monomial == 0:
            return monomial
        exp = monomial.exponents()[0][:self.nx]
        for i, k in enumerate(exp):
            if k == 1:
                return self.variables[i]

    def biggest_y(self, monomial):
        if monomial == 1 or monomial == 0:
            return monomial
        exp = monomial.exponents()[0][self.ny:]
        for i, k in enumerate(exp):
            if k == 1:
                return self.variables[self.ny + i]

    def add_lines(self, f_i, i, Mac_d_1):
        """
        Adds all the lines of signature (u, f_i)
        s.t. deg(uf_i) == d
        """
        if Mac_d_1.d == (self.d[0] - 1, self.d[1]):
            degree_to_add = 0
            variables = self.variables[:self.poly_ring.ngens() // 2]
        elif Mac_d_1.d == (self.d[0], self.d[1] - 1):
            degree_to_add = 1
            variables = self.variables[self.poly_ring.ngens() // 2:]

        for row_i, (e, f_ii) in enumerate(Mac_d_1.sig):
            if f_ii < i:
                continue
            elif f_ii > i:
                return
            if degree_to_add == 0:
                x_lambda = self.biggest_x(e)
            elif degree_to_add == 1:
                x_lambda = self.biggest_y(e)
            else:
                print("Oula !! grosse erreur")

            for x_i in variables:
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

        #return nnz_columns - nnz_rows
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

def integer_vectors_at_most(max_sum):
    """Génère tous les vecteurs d'entiers de longueur donnée avec somme ≤ max_sum"""
    result = []
    for s in range(max_sum + 1):
        for v in IntegerVectors(s, 2):
            result.append(v)
    return result

def F5Matrix(F, dmax, nx, ny):
    """
    F is homogeneous polynomials ordered such that deg(f_i) < deg(f_j) forall i, j
    such that i < j
    dmax is the maximal degree to obtain a d_max-grobner basis of F
    F is supposed to be a quadratic system, so deg(f_i) <= 2 forall 0 <= i < m
    We follow F5Matrix by Bardet in her PhD thesis p.21
    """
    from sage.rings.polynomial.polydict import ETuple
    m = len(F)
    R = F[0].base_ring()
    n = F[0].parent().ngens()
    Mac_d = None
    Mac_d_old = []
    gb = []
    bi_hilbert = hilbert_biseries(nx, ny, len(F)).monomial_coefficients()
    
    degs = integer_vectors_at_most(dmax)
    degs.pop(0)

    print(f"F5 bihomogeneous for d={F[0].total_degree()}...{dmax}")

    for (d1, d2) in degs:
        print(f"\n{'-' * 20} d=({d1}, {d2}) {'-' * 20}\n")
        t_deg_start = time.time()
        Mac_d = Mac(d1, d2, R, nx, ny)
        Mac_d.monomial_ordered_list_deg_d(d1 + d2)
        Mac_d_1 = None
        t0 = time.time()
        if (d1 == 0 and d2 != 0) or (d2 == 0 and d1 != 0):
            try:
                print(f"Corank of degree ({d1}, {d2}): {len(Mac_d.monomial_inverse_search)} ----- Value in Hilbert Biseries: {bi_hilbert[ETuple([d1, d2])]}")
            except KeyError:
                continue
            continue
        elif d1+d2 != 2:
            for M in Mac_d_old:
                if M.d == (d1 - 1, d2) or M.d == (d1, d2 - 1):
                    Mac_d_1 = M
                    print(f"Matrix used for Mac_d_1 is of degree: {Mac_d_1.d}")
                    break
        for i in range(0, m):
            f_i = F[i]
            if d1+d2 == 2:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, Mac_d_1)
        t1 = time.time()
        print(f"[TIMER] Temps pour add_lines : {t1 - t0:.4f} s")
        #tmp_Mac = copy.deepcopy(Mac_d)
        Mac_d.gauss()
        reductions_to_zero = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree ({d1}, {d2}): {reductions_to_zero} / {Mac_d.matrix.nrows()}")
        try:
            print(f"Corank of degree ({d1}, {d2}): {Mac_d.corank()} ----- Value in Hilbert Biseries: {bi_hilbert[ETuple([d1, d2])]}")
        except KeyError:
            print(f"Corank of degree ({d1}, {d2}): {Mac_d.corank()} ----- Value in Hilbert Biseries: 0")
        #update_gb(gb, tmp_Mac, Mac_d)
        Mac_d_old.append(Mac_d)

        t_deg_end = time.time()
        print(f"[TIMER] Temps total pour bi-degré ({d1}, {d2}) : {t_deg_end - t_deg_start:.2f} s")
    return gb

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

def Nm(n, m, t1, t2):
    """
    internal function for the hilbert_biseries
    """
    sum_l = 0
    for l in range(1, m - n + 1):
        sum_k = 0
        for k in range(1, n + 1):
            sum_k += (t1^(n - k))*binomial(l + n - k - 1, n - k)
        bracket = 1 - (1 - t1)**l * sum_k
        term = ((1 - t1*t2)^(m - n - l)) * t1*t2*((1 - t2)^(n))
        sum_l += term * bracket
    return sum_l

def hilbert_biseries(nx, ny, m):
    """
    Returns the hilbert bi-series for a quadratic
    system of nx + ny variables and m polynomials
    result from https://arxiv.org/abs/1001.4004
    """
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) +2)
    denom = ((1 - tx)^(nx)) * ((1 - ty)^(ny))
    num = (1 - tx*ty)^m + Nm(ny, m, tx, ty) + Nm(nx, m, ty, tx)
    return num / denom

if __name__ == '__main__':
    """
    test = integer_vectors_at_most(5)
    for i in test:
        print(i)
    """
    F = homogenized_ideal(doit_bilinear(7, 5, 12))
    """
    R.<x1, x2, x3, y1, y2, y3, y4> = PolynomialRing(GF(7), order='degrevlex')
    F = [x1*y1 + 5*x2*y1 + 4*x3*y1 + 5*x1*y2 + 3*x2*y2 + x1*y3 + 4*x2*y3 + 5*x3*y3 + 5*x1*y4 + x2*y4 + 2*x3*y4,
    2*x1*y1 + 4*x2*y1 + 6*x3*y1 + 2*x1*y2 + 5*x2*y2 + 6*x1*y3 + 4*x3*y3 + 3*x1*y4 + 2*x2*y4 + 4*x3*y4,
    5*x1*y1 + 5*x2*y1 + 2*x3*y1 + 4*x1*y2 + 6*x2*y2 + 4*x3*y2 + 6*x2*y3 + 4*x3*y3 + x1*y4 + x2*y4 + 5*x3*y4,
    6*x1*y1 + 5*x3*y1 + 4*x1*y2 + 5*x2*y2 + x3*y2 + x1*y3 + x2*y3 + 6*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4,
    6*x1*y1 + 3*x2*y1 + 6*x3*y1 + 3*x1*y2 + 5*x3*y2 + 2*x1*y3 + 4*x2*y3 + 5*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4
    ]
    """
    #F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_14_15.sobj"))
    #D = Ideal(F).degree_of_semi_regularity()
    print(f"---------------Generating Serie Bardet: {generating_bardet_series(F)}\n\n")
    series_ring.<z> = PowerSeriesRing(ZZ)
    hilbert_series = series_ring(Ideal(F).hilbert_series())
    print(f"---------------Hilbert Series: {hilbert_series}\n\n")

    nx = 7
    ny = 5
    print(f"---------------degree of semi-regularity of F: {D}\n\n")
    print(f"---------------Série génératrice bilinéaire: {hilbert_biseries(nx, ny, len(F))}\n\n")

    gb = F5Matrix(F, min(nx, ny) + 2, nx, ny)

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