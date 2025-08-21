load("base_F5.sage")
#load("REF.sage")
import copy
import time
import cProfile

class MacHom(BaseMac):
    def __init__(self, n, d, F, h):
        self.n = n
        self.d = d
        #monomials_str = ['x' + str(i) + 'bar' for i in range(1, n//2 + 1)] + ['y' + str(i) + 'bar' for i in range(1, n//2 + 1)]
        monomials_str = ['x' + str(i) for i in range(1, n//2 + 1)] + ['y' + str(i) for i in range(1, n//2 + 1)]
        self.monomial_hash_list = {}  # Monomial -> column index
        self.monomial_inverse_search = []  # Column index -> monomial
        self.poly_ring = PolynomialRing(F, monomials_str, order='degrevlex')
        self.variables = [self.poly_ring(m) for m in monomials_str]
        self.monomial_ordered_list_deg_d(d)
        super().__init__(F, monomials_str, h)

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all monomials of degree d and add them to hash list
        Fonction testé -> Correcte
        """
        poly_ring = self.poly_ring

        monomials = poly_ring.monomials_of_degree(d)
        monomials_in_R = sorted([poly_ring(m) for m in monomials], reverse=True)

        hash_size = len(self.monomial_hash_list)
        self.monomial_hash_list = {m: i for i, m in enumerate(monomials_in_R)} | self.monomial_hash_list
        self.monomial_inverse_search += monomials_in_R

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        Fonction testé -> Correcte
        """
        vec = self.polynomial_to_vector(self.poly_ring(u*f))
        self.add_row(vec)
        self.sig.append((u, i))
        return

    def add_lines(self, f_i, i, Mac_d_1):
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
            
            first = False
            for x_i in self.variables:
                if x_i >= x_lambda:
                    u = x_i * e
                    if not any(u in self.crit[ii] for ii in range(i)):
                        if first:
                            self.matrix[self.matrix_index] = np.roll(self.matrix[self.matrix_index - 1], 1)
                        else:
                            self.add_line(f_i, i, u)
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
        Mac_d = MacHom(n, d, R, hilbert_series[d])
        Mac_d.monomial_ordered_list_deg_d(d)
        t0 = time.time()
        Mac_d.F5_criterion(Mac_d_2)
        for i in range(0, m):
            f_i = F[i]
            if f_i.total_degree() == d:
                Mac_d.add_line(f_i, i, 1)
            elif f_i.total_degree() > d:
                continue
            else:
                Mac_d.add_lines(f_i, i, Mac_d_1)
        t1 = time.time()
        print(f"[TIMER] Temps pour add_lines : {t1 - t0:.4f} s")
        #tmp_Mac = copy.deepcopy(Mac_d)
        t0 = time.time()
        Mac_d.gauss(debug=True)
        t1 = time.time()
        print(f"[TIMER] Temps pour Gauss (matrice {Mac_d.nrows}x{Mac_d.ncols}) : {t1 - t0:.4f} s")
        #update_gb(gb, tmp_Mac, Mac_d)
        reductions_to_zero, lignes_0 = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree {d} with normal Gauss: {reductions_to_zero} / {Mac_d.nrows}")
        print(f"Corank of degree {d}: {Mac_d.ncols- Mac_d.nrows + reductions_to_zero} ----- Value in Hilbert Series: {hilbert_series[d]}")

        if reductions_to_zero > 0:
            for i in lignes_0:
                print(Mac_d.sig[i])
                #print(depend_matrix[i])
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
        
        # Vérification row_lm() cohérence avec le vecteur original
        for idx in range(Mac_d.matrix.nrows()):
            lm_vec = Mac_d.row_lm(idx)
            poly = Mac_d.vector_to_polynomial(idx)
            try:
                lm_poly = poly.lm()
            except AttributeError:
                lm_poly = 0
            if lm_vec != lm_poly:
                print(f"[ERREUR row_lm] ligne {idx} : row_lm()={lm_vec}, poly.lm()={lm_poly}")
            #else:
            #    print(f"[OK] row_lm ligne {idx} correcte : {lm_vec}")
        """

    return gb

if __name__ == '__main__':
    """
    R.<x1, x2, x3> = PolynomialRing(GF(5), order='degrevlex')
    F = [x2^2 + 4*x2*x3,
    2*x1^2 + 3*x1*x2 + 4*x2^2 + 3*x3^2,
    3*x1^2 + 4*x1*x2 + 2*x2^2]
    """
    
    F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_16_17.sobj"))
    #F = homogenized_ideal(doit_bilinear(4, 4, 5))
    #F = homogenized_ideal(doit(10, 11))
    #series_ring.<z> = PowerSeriesRing(ZZ)
    #hilbert_series = series_ring(Ideal(F).hilbert_series())
    #print(f"Hilbert Series: {hilbert_series}")
    #D = Ideal(F).degree_of_semi_regularity()
    #print(generating_bardet_series(F))
    #print(f"degree of semi-regularity of F: {D}")

    F = sorted(F, key=lambda f: f.degree())
    F5Matrix(F, 5)
    #cProfile.run("F5Matrix(F, D)", sort='cumtime')

    G = Ideal(F).groebner_basis(algorithm="msolve")
    deg_max = 0
    for p in G:
        d = p.total_degree()
        if d > deg_max:
            deg_max = d
    print(deg_max)

    """
    print(len(gb))
    print(Ideal(gb).basis_is_groebner())


    #Test linéarisation avec Frobenius
    F = load("../MPCitH_SBC/system/sage/system_bilin_36_37_test.sobj")
    m = len(F)
    R = F[0].base_ring()
    n = F[0].parent().ngens()
    Mac = MacHom(n, 2, R, m)
    Mac.monomial_ordered_list_deg_d(2)
    Mac.monomial_ordered_list_deg_d(1)
    Mac.monomial_ordered_list_deg_d(0)
    for i, f in enumerate(F):
        Mac.add_line(f, i, 1)
    Mac.matrix = matrix(GF(2), Mac.matrix)
    print(Mac.matrix.nrows())
    print(Mac.matrix.ncols())
    print(Mac.matrix.rank())
    #print(F)
    """