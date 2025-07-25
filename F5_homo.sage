load("base_F5.sage")
import copy
import time

class MacHom(BaseMac):
    def __init__(self, n, d, F):
        self.n = n
        self.d = d
        monomials_str = ['x' + str(i) for i in range(1, n//2 + 1)] + ['y' + str(i) for i in range(1, n//2 + 1)]
        super().__init__(F, monomials_str)

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

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        Fonction testé -> Correcte
        """
        vec = self.polynomial_to_vector(self.quotient_ring(u*f))
        self.add_row(vec)
        self.sig.append((u, i))
        return

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

    hilbert_series = generating_bardet_series(F).list()

    #print("Polynomials in F:")
    #for i in F:
    #    print(i)

    print(f"F5 for d={F[0].total_degree()}...{dmax}")

    for d in range(F[0].total_degree(), dmax+1):
        print(f"\n{'-' * 20} d={d} {'-' * 20}\n")
        t_deg_start = time.time()
        Mac_d = MacHom(n, d, R)
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
        print(f"Corank of degree {d}: {Mac_d.corank()} ----- Value in Hilbert Series: {hilbert_series[d]}")
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