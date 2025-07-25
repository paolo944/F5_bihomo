load("base_F5.sage")
import copy
import time

def bi_degree(monomial, nx, ny):
    exponents = monomial.exponents()[0]
    d1 = sum(exponents[:nx])
    d2 = sum(exponents[nx:nx+ny])
    return (d1, d2)

def number_of_monomials_bidegree(nx, ny, d1, d2):
    y = binomial(ny + d2 - 1, d2)
    x = binomial(nx + d1 - 1, d1)
    return x*y

class MacBiHom(BaseMac):
    def __init__(self, d1, d2, F, nx, ny):
        self.d = (d1, d2)
        self.nx = nx
        self.ny = ny
        monomials_str = ['x' + str(i) for i in range(1, nx + 1)] + ['y' + str(i) for i in range(1, ny + 1)]
        super().__init__(F, monomials_str)

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

        self.monomial_hash_list = {m: i for i, m in enumerate(monomials_in_R)} | self.monomial_hash_list
        self.monomial_inverse_search += monomials_in_R

        assert len(self.monomial_inverse_search) == len(self.monomial_hash_list) and len(self.monomial_hash_list) == number_of_monomials_bidegree(self.nx, self.ny, self.d[0], self.d[1]), "Erreur dans la génération de monomes, le compte n'est pas bon" 

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
        Mac_d = MacBiHom(d1, d2, R, nx, ny)
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
        Mac_d.matrix = matrix(GF(2), Mac_d.matrix, sparse=False)
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

def hilbert_biseries(nx, ny, m):
    """
    Returns the hilbert bi-series for a quadratic
    system of nx + ny variables and m polynomials
    result from https://arxiv.org/abs/1001.4004
    """
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
    #F = homogenized_ideal(doit_bilinear(7, 5, 12))
    """
    R.<x1, x2, x3, y1, y2, y3, y4> = PolynomialRing(GF(7), order='degrevlex')
    F = [x1*y1 + 5*x2*y1 + 4*x3*y1 + 5*x1*y2 + 3*x2*y2 + x1*y3 + 4*x2*y3 + 5*x3*y3 + 5*x1*y4 + x2*y4 + 2*x3*y4,
    2*x1*y1 + 4*x2*y1 + 6*x3*y1 + 2*x1*y2 + 5*x2*y2 + 6*x1*y3 + 4*x3*y3 + 3*x1*y4 + 2*x2*y4 + 4*x3*y4,
    5*x1*y1 + 5*x2*y1 + 2*x3*y1 + 4*x1*y2 + 6*x2*y2 + 4*x3*y2 + 6*x2*y3 + 4*x3*y3 + x1*y4 + x2*y4 + 5*x3*y4,
    6*x1*y1 + 5*x3*y1 + 4*x1*y2 + 5*x2*y2 + x3*y2 + x1*y3 + x2*y3 + 6*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4,
    6*x1*y1 + 3*x2*y1 + 6*x3*y1 + 3*x1*y2 + 5*x3*y2 + 2*x1*y3 + 4*x2*y3 + 5*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4
    ]
    """
    F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_14_15.sobj"))
    #D = Ideal(F).degree_of_semi_regularity()
    print(f"---------------Generating Serie Bardet: {generating_bardet_series(F)}\n\n")
    series_ring.<z> = PowerSeriesRing(ZZ)
    hilbert_series = series_ring(Ideal(F).hilbert_series())
    print(f"---------------Hilbert Series: {hilbert_series}\n\n")

    nx = 7
    ny = 7
    print(f"---------------degree of semi-regularity of F: {D}\n\n")
    print(f"---------------Série génératrice bilinéaire: {hilbert_biseries(nx, ny, len(F))}\n\n")

    gb = F5Matrix(F, min(nx, ny) + 2, nx, ny)