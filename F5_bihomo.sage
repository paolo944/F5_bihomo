load("base_F5.sage")
import copy
import time
import cProfile
import pack_utils

def number_of_monomials_bidegree(nx, ny, d1, d2):
    y = binomial(ny + d2 - 1, d2)
    x = binomial(nx + d1 - 1, d1)
    return x*y

def complexity_bi(nx, ny, d, omega):
    degs = integer_vectors_at_most(d)
    degs.pop(0)
    c = 0
    for (d1, d2) in degs:
        c += number_of_monomials_bidegree(nx, ny, d1, d2) ^ omega
    return c

def complexity_homo(nx, ny, d, omega):
    c = binomial(nx + ny + d - 1, d) ^ omega
    return c

def integer_vectors_at_most(max_sum):
    """Génère tous les vecteurs d'entiers de longueur donnée avec somme ≤ max_sum"""
    result = []
    for s in range(max_sum + 1):
        for v in IntegerVectors(s, 2):
            result.append(v)
    return result

def verif_antisym(A):
    m = A.nrows()
    n = A.ncols()
    min_dim = min(m, n)

    for i in range(min_dim):
        for j in range(min_dim):
            if A[i,j] != -A[j,i]:
                print(f"A[{i},{j}] = {A[i,j]} n'est pas l'opposé de A[{j},{i}] = {A[j,i]}")
                #return False
    return True

class MacBiHom(BaseMac):
    def __init__(self, d1, d2, F, nx, ny, h):
        self.d = (d1, d2)
        self.nx = nx
        self.ny = ny
        monomials_str = ['x' + str(i) for i in range(1, nx + 1)] + ['y' + str(i) for i in range(1, ny + 1)]
        self.monomial_hash_list = {}  # Monomial -> column index
        self.monomial_inverse_search = []  # Column index -> monomial
        self.poly_ring = PolynomialRing(F, monomials_str, order='degrevlex')
        self.variables = [self.poly_ring(mon) for mon in monomials_str]
        self.monomial_ordered_list_deg_d((d1, d2))
        super().__init__(F, monomials_str, h)

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all monomials of degree d and add them to hash list
        Fonction testé -> Correcte
        """
        d1, d2 = d
        poly_ring = self.poly_ring
        nx, ny = self.nx, self.ny

        # Génération directe des exposants de bidegré (d1, d2)
        Lx = IntegerVectors(d1, nx)
        Ly = IntegerVectors(d2, ny)

        monomials_in_R = []
        for expx in Lx:
            for expy in Ly:
                exponents = list(expx) + list(expy)
                monomial = poly_ring.monomial(*exponents)
                monomials_in_R.append(monomial)

        # Tri décroissant
        monomials_in_R.sort(reverse=True)

        # Mise à jour des structures de hachage
        self.monomial_hash_list = {m: i for i, m in enumerate(monomials_in_R)} | self.monomial_hash_list
        self.monomial_inverse_search += monomials_in_R

        #assert len(self.monomial_inverse_search) == len(self.monomial_hash_list) and len(self.monomial_hash_list) == number_of_monomials_bidegree(self.nx, self.ny, self.d[0], self.d[1]), "Erreur dans la génération de monomes, le compte n'est pas bon" 

    def add_line(self, f, i, u):
        """
        Add line (u, f_i) to Macaulay matrix
        """
        poly = u*f
        if pack_utils.bi_degree(poly, self.nx, self.ny) != self.d:
            return
        vec = self.polynomial_to_vector(u*f)
        self.add_row(vec)
        self.sig.append((u, i))
        return
    
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

    def add_lines(self, f_i, i, Mac_d_1, Mac_d_2):
        """
        Adds all the lines of signature (u, f_i)
        s.t. deg(uf_i) == d
        """
        if Mac_d_1.d == (self.d[0] - 1, self.d[1]):
            degree_to_add = 0
            variables = self.variables[:self.nx]
        elif Mac_d_1.d == (self.d[0], self.d[1] - 1):
            degree_to_add = 1
            variables = self.variables[self.nx:self.nx+self.ny]

        for row_i, (e, f_ii) in enumerate(Mac_d_1.sig):
            if f_ii < i:
                continue
            elif f_ii > i:
                return
            
            x_lambda = (
            self.biggest_x(e) if degree_to_add == 0 else self.biggest_y(e)
            )

            for x_i in variables:
                if x_i >= x_lambda:
                    u = x_i * e
                    if not any(u in self.crit[ii] for ii in range(i)):
                        self.add_line(f_i, i, u)
        return

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
    Mac_d_2 = []
    gb = []
    bi_hilbert = hilbert_biseries(nx, ny, len(F)).monomial_coefficients()
    
    degs = integer_vectors_at_most(dmax)
    degs.pop(0)

    print(f"F5 bihomogeneous for d={F[0].total_degree()}...{dmax}")

    for (d1, d2) in degs:
        try: 
            h = bi_hilbert[ETuple([d1, d2])]
        except KeyError:
            h = 0
        print(f"\n{'-' * 20} d=({d1}, {d2}) {'-' * 20}\n")
        t_deg_start = time.time()
        Mac_d = MacBiHom(d1, d2, R, nx, ny, h)
        #Mac_d.monomial_ordered_list_deg_d(d1 + d2)
        Mac_d_1 = None
        t0 = time.time()
        if (d1 == 0 and d2 != 0) or (d2 == 0 and d1 != 0):
            print(f"Corank of degree ({d1}, {d2}): {len(Mac_d.monomial_inverse_search)} ----- Value in Hilbert Biseries: {h}")
            continue
        elif d1+d2 != 2:
            for M in Mac_d_old:
                if M.d == (d1 - 1, d2) or M.d == (d1, d2 - 1):
                    Mac_d_1 = M
                    break
        if d1 + d2 > 3:
            if Mac_d_old[-1].d[0] + Mac_d_old[-1].d[1] == d1 + d2:
                print("Je réutilise crit")
                Mac_d.crit = Mac_d_old[-1].crit
            else:
                for M in Mac_d_old:
                    if M.d[0] + M.d[1] == d1 + d2 - 2:
                        print(f"Added Mac_{M.d[0]}_{M.d[1]}")
                        Mac_d.F5_criterion(M)
        for i in range(0, m):
            f_i = F[i]
            if d1 + d2 == 2:
                Mac_d.add_line(f_i, i, 1)
            else:
                Mac_d.add_lines(f_i, i, Mac_d_1, Mac_d_old)
        t1 = time.time()
        print(f"[TIMER] Temps pour add_lines : {t1 - t0:.4f} s")
        #tmp_Mac = copy.deepcopy(Mac_d)
        #Mac_d.matrix = matrix(GF(2), Mac_d.matrix, sparse=False)
        print(f"Final matrix before gaussian elimination: {Mac_d.nrows}x{Mac_d.ncols}")
        t0 = time.time()
        Mac_d.gauss()
        t1 = time.time()
        print(f"[TIMER] Temps pour Gauss (matrice {Mac_d.nrows}x{Mac_d.ncols}) : {t1 - t0:.4f} s")
        reductions_to_zero, lignes_a_0 = Mac_d.verify_reductions_zero()
        print(f"number of reductions to 0 in degree ({d1}, {d2}): {reductions_to_zero} / {Mac_d.nrows}")
        print(f"Total non zero lines in: {Mac_d.nrows - reductions_to_zero} expected: {h + Mac_d.ncols}")
        print(f"Corank of degree ({d1}, {d2}): {Mac_d.ncols - Mac_d.nrows + reductions_to_zero} ----- Value in Hilbert Biseries: {h}")
        #update_gb(gb, tmp_Mac, Mac_d)
        #if reductions_to_zero > 0:
        #    for i in lignes_a_0:
        #        print(Mac_d.sig[i])

        Mac_d_old.append(Mac_d)

        t_deg_end = time.time()
        print(f"[TIMER] Temps total pour bi-degré ({d1}, {d2}) : {t_deg_end - t_deg_start:.2f} s")
    return gb

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
    
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) + 5)
    denom = ((1 - tx)^(nx)) * ((1 - ty)^(ny))
    num = (1 - tx*ty)^m + Nm(ny, m, tx, ty) + Nm(nx, m, ty, tx)
    return num / denom

if __name__ == '__main__':
    """
    test = integer_vectors_at_most(5)
    for i in test:
        print(i)
    """
    """
    R.<x1, x2, x3, y1, y2, y3, y4> = PolynomialRing(GF(7), order='degrevlex')
    F = [x1*y1 + 5*x2*y1 + 4*x3*y1 + 5*x1*y2 + 3*x2*y2 + x1*y3 + 4*x2*y3 + 5*x3*y3 + 5*x1*y4 + x2*y4 + 2*x3*y4,
    2*x1*y1 + 4*x2*y1 + 6*x3*y1 + 2*x1*y2 + 5*x2*y2 + 6*x1*y3 + 4*x3*y3 + 3*x1*y4 + 2*x2*y4 + 4*x3*y4,
    5*x1*y1 + 5*x2*y1 + 2*x3*y1 + 4*x1*y2 + 6*x2*y2 + 4*x3*y2 + 6*x2*y3 + 4*x3*y3 + x1*y4 + x2*y4 + 5*x3*y4,
    6*x1*y1 + 5*x3*y1 + 4*x1*y2 + 5*x2*y2 + x3*y2 + x1*y3 + x2*y3 + 6*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4,
    6*x1*y1 + 3*x2*y1 + 6*x3*y1 + 3*x1*y2 + 5*x3*y2 + 2*x1*y3 + 4*x2*y3 + 5*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4
    ]
    """
    nx = 6
    ny = 6
    m = 13

    F = doit_bilinear(nx, ny, m)
    print(F)
    F = homogenized_ideal(F)

    #F = homogenized_ideal(load(f"../MPCitH_SBC/system/sage/system_bilin_{nx + ny}_{2*(nx + ny) + 1}.sobj"))
    #D = Ideal(F).degree_of_semi_regularity()
    #print(f"---------------Generating Serie Bardet: {generating_bardet_series(F)}\n\n")
    series_ring.<z> = PowerSeriesRing(ZZ)
    hilbert_series = series_ring(Ideal(F).hilbert_series())
    print(f"---------------Hilbert Series: {hilbert_series}\n\n")#

    
    #print(f"---------------degree of semi-regularity of F: {D}\n\n")
    print(f"---------------Série génératrice bilinéaire: {hilbert_biseries(nx, ny, len(F))}\n\n")

    F5Matrix(F, min(nx, ny) + 2, nx, ny)
    #cProfile.run("F5Matrix(F, min(nx, ny) + 2, nx, ny)", sort='cumtime')

    """
    for nx in range(12, 16):
        F = homogenized_ideal(load(f"../MPCitH_SBC/system/sage/system_bilin_{nx * 2}_{nx * 2 + 1}.sobj"))
        ny = nx
        print(f"---------------Série génératrice bilinéaire: {hilbert_biseries(nx, ny, len(F))}\n\n")
        F5Matrix(F, nx - 5, nx, ny)
    """