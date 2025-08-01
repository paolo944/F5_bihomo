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
    def __init__(self, d1, d2, F, nx, ny, m):
        self.d = (d1, d2)
        self.nx = nx
        self.ny = ny
        monomials_str = ['x' + str(i) for i in range(1, nx + 1)] + ['y' + str(i) for i in range(1, ny + 1)]
        super().__init__(F, monomials_str, m)

    def monomial_ordered_list_deg_d(self, d):
        """
        Generate all monomials of degree d and add them to hash list
        Fonction testé -> Correcte
        """
        poly_ring = self.poly_ring

        monomials = poly_ring.monomials_of_degree(d)
        monomials_in_R = []
        for monomial in monomials:
            if bi_degree(monomial, self.nx, self.ny) == self.d:
                monomials_in_R.append(poly_ring(monomial))
        
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
                        u = x_i * e
                        crit = False
                        for (ii, lm) in self.crit:
                            if ii < i and lm == u:
                                crit = True
                        if not crit:
                            self.add_line(f_i, i, u)
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
    Mac_d_2 = []
    gb = []
    bi_hilbert = hilbert_biseries(nx, ny, len(F)).monomial_coefficients()
    
    degs = integer_vectors_at_most(dmax)
    degs.pop(0)

    print(f"F5 bihomogeneous for d={F[0].total_degree()}...{dmax}")

    for (d1, d2) in degs:
        print(f"\n{'-' * 20} d=({d1}, {d2}) {'-' * 20}\n")
        t_deg_start = time.time()
        Mac_d = MacBiHom(d1, d2, R, nx, ny, m)
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
                    break
        if d1 + d2 > 3:
            for M in Mac_d_old:
                if M.d[0] + M.d[1] == d1 + d2 - 2:
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
        Mac_d.gauss()
        reductions_to_zero, lignes_a_0 = Mac_d.verify_reductions_zero()
        ncols = len(Mac_d.monomial_hash_list)
        nrows = len(Mac_d.matrix) // ncols
        print(f"number of reductions to 0 in degree ({d1}, {d2}): {reductions_to_zero} / {nrows}")
        try:
            print(f"Corank of degree ({d1}, {d2}): {ncols - nrows + reductions_to_zero} ----- Value in Hilbert Biseries: {bi_hilbert[ETuple([d1, d2])]}")
        except KeyError:
            print(f"Corank of degree ({d1}, {d2}): {ncols - nrows + reductions_to_zero} ----- Value in Hilbert Biseries: 0")
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
    #F = homogenized_ideal(doit_bilinear(8, 8, 17))
    """
    R.<x1, x2, x3, y1, y2, y3, y4> = PolynomialRing(GF(7), order='degrevlex')
    F = [x1*y1 + 5*x2*y1 + 4*x3*y1 + 5*x1*y2 + 3*x2*y2 + x1*y3 + 4*x2*y3 + 5*x3*y3 + 5*x1*y4 + x2*y4 + 2*x3*y4,
    2*x1*y1 + 4*x2*y1 + 6*x3*y1 + 2*x1*y2 + 5*x2*y2 + 6*x1*y3 + 4*x3*y3 + 3*x1*y4 + 2*x2*y4 + 4*x3*y4,
    5*x1*y1 + 5*x2*y1 + 2*x3*y1 + 4*x1*y2 + 6*x2*y2 + 4*x3*y2 + 6*x2*y3 + 4*x3*y3 + x1*y4 + x2*y4 + 5*x3*y4,
    6*x1*y1 + 5*x3*y1 + 4*x1*y2 + 5*x2*y2 + x3*y2 + x1*y3 + x2*y3 + 6*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4,
    6*x1*y1 + 3*x2*y1 + 6*x3*y1 + 3*x1*y2 + 5*x3*y2 + 2*x1*y3 + 4*x2*y3 + 5*x3*y3 + 2*x1*y4 + 4*x2*y4 + 5*x3*y4
    ]
    """
    F = homogenized_ideal(load("../MPCitH_SBC/system/sage/system_bilin_20_21.sobj"))
    #D = Ideal(F).degree_of_semi_regularity()
    print(f"---------------Generating Serie Bardet: {generating_bardet_series(F)}\n\n")
    #series_ring.<z> = PowerSeriesRing(ZZ)
    #hilbert_series = series_ring(Ideal(F).hilbert_series())
    #print(f"---------------Hilbert Series: {hilbert_series}\n\n")

    nx = 10
    ny = 10
    #print(f"---------------degree of semi-regularity of F: {D}\n\n")
    print(f"---------------Série génératrice bilinéaire: {hilbert_biseries(nx, ny, len(F))}\n\n")

    gb = F5Matrix(F, min(nx, ny) + 2, nx, ny)

    """
    ---------------Generating Serie Bardet: 1 + 20*z + 189*z^2 + 1120*z^3 + 4655*z^4 + 14364*z^5 + 33915*z^6 + 62016*z^7 + 87210*z^8 + 90440*z^9 + 58786*z^10 - 58786*z^12 - 90440*z^13 - 87210*z^14 - 62016*z^15 - 33915*z^16 - 14364*z^17 - 4655*z^18 - 1120*z^19 + O(z^20)


---------------Série génératrice bilinéaire: 1 + 10*tx + 10*ty + 55*tx^2 + 79*tx*ty + 55*ty^2 + 220*tx^3 + 340*tx^2*ty + 340*tx*ty^2 + 220*ty^3 + 715*tx^4 + 1045*tx^3*ty + 1135*tx^2*ty^2 + 1045*tx*ty^3 + 715*ty^4 + 2002*tx^5 + 2530*tx^4*ty + 2650*tx^3*ty^2 + 2650*tx^2*ty^3 + 2530*tx*ty^4 + 2002*ty^5 + 5005*tx^6 + 5005*tx^5*ty + 4675*tx^4*ty^2 + 4545*tx^3*ty^3 + 4675*tx^2*ty^4 + 5005*tx*ty^5 + 5005*ty^6 + 11440*tx^7 + 8008*tx^6*ty + 6160*tx^5*ty^2 + 5400*tx^4*ty^3 + 5400*tx^3*ty^4 + 6160*tx^2*ty^5 + 8008*tx*ty^6 + 11440*ty^7 + 24310*tx^8 + 9295*tx^7*ty + 5005*tx^6*ty^2 + 3465*tx^5*ty^3 + 3060*tx^4*ty^4 + 3465*tx^3*ty^5 + 5005*tx^2*ty^6 + 9295*tx*ty^7 + 24310*ty^8 + 48620*tx^9 + 2860*tx^8*ty - 1430*tx^7*ty^2 - 2310*tx^6*ty^3 - 2520*tx^5*ty^4 - 2520*tx^4*ty^5 - 2310*tx^3*ty^6 - 1430*tx^2*ty^7 + 2860*tx*ty^8 + 48620*ty^9 + 92378*tx^10 - 24310*tx^9*ty - 14300*tx^8*ty^2 - 10725*tx^7*ty^3 - 9240*tx^6*ty^4 - 8820*tx^5*ty^5 - 9240*tx^4*ty^6 - 10725*tx^3*ty^7 - 14300*tx^2*ty^8 - 24310*tx*ty^9 + 92378*ty^10 + 167960*tx^11 - 97240*tx^10*ty - 28600*tx^9*ty^2 - 17160*tx^8*ty^3 - 13200*tx^7*ty^4 - 11760*tx^6*ty^5 - 11760*tx^5*ty^6 - 13200*tx^4*ty^7 - 17160*tx^3*ty^8 - 28600*tx^2*ty^9 - 97240*tx*ty^10 + 167960*ty^11 + O(tx, ty)^12


F5 bihomogeneous for d=2...12

-------------------- d=(1, 0) --------------------

Corank of degree (1, 0): 10 ----- Value in Hilbert Biseries: 10

-------------------- d=(0, 1) --------------------

Corank of degree (0, 1): 10 ----- Value in Hilbert Biseries: 10

-------------------- d=(2, 0) --------------------

Corank of degree (2, 0): 55 ----- Value in Hilbert Biseries: 55

-------------------- d=(1, 1) --------------------

[TIMER] Temps pour add_lines : 0.0070 s
[TIMER] Temps pour Gauss (matrice 21x100) : 0.0007 s
number of reductions to 0 in degree (1, 1): 0 / 21
Corank of degree (1, 1): 79 ----- Value in Hilbert Biseries: 79
[TIMER] Temps total pour bi-degré (1, 1) : 0.01 s

-------------------- d=(0, 2) --------------------

Corank of degree (0, 2): 55 ----- Value in Hilbert Biseries: 55

-------------------- d=(3, 0) --------------------

Corank of degree (3, 0): 220 ----- Value in Hilbert Biseries: 220

-------------------- d=(2, 1) --------------------

[TIMER] Temps pour add_lines : 0.0457 s
[TIMER] Temps pour Gauss (matrice 210x550) : 0.0300 s
number of reductions to 0 in degree (2, 1): 0 / 210
Corank of degree (2, 1): 340 ----- Value in Hilbert Biseries: 340
[TIMER] Temps total pour bi-degré (2, 1) : 0.11 s

-------------------- d=(1, 2) --------------------

[TIMER] Temps pour add_lines : 0.0434 s
[TIMER] Temps pour Gauss (matrice 210x550) : 0.0192 s
number of reductions to 0 in degree (1, 2): 0 / 210
Corank of degree (1, 2): 340 ----- Value in Hilbert Biseries: 340
[TIMER] Temps total pour bi-degré (1, 2) : 0.10 s

-------------------- d=(0, 3) --------------------

Corank of degree (0, 3): 220 ----- Value in Hilbert Biseries: 220

-------------------- d=(4, 0) --------------------

Corank of degree (4, 0): 715 ----- Value in Hilbert Biseries: 715

-------------------- d=(3, 1) --------------------

[TIMER] Temps pour add_lines : 0.6132 s
[TIMER] Temps pour Gauss (matrice 1155x2200) : 0.7996 s
number of reductions to 0 in degree (3, 1): 0 / 1155
Corank of degree (3, 1): 1045 ----- Value in Hilbert Biseries: 1045
[TIMER] Temps total pour bi-degré (3, 1) : 1.62 s

-------------------- d=(2, 2) --------------------

[TIMER] Temps pour add_lines : 1.9418 s
[TIMER] Temps pour Gauss (matrice 1890x3025) : 1.1065 s
number of reductions to 0 in degree (2, 2): 0 / 1890
Corank of degree (2, 2): 1135 ----- Value in Hilbert Biseries: 1135
[TIMER] Temps total pour bi-degré (2, 2) : 3.26 s

-------------------- d=(1, 3) --------------------

[TIMER] Temps pour add_lines : 0.3688 s
[TIMER] Temps pour Gauss (matrice 1155x2200) : 0.3840 s
number of reductions to 0 in degree (1, 3): 0 / 1155
Corank of degree (1, 3): 1045 ----- Value in Hilbert Biseries: 1045
[TIMER] Temps total pour bi-degré (1, 3) : 0.97 s

-------------------- d=(0, 4) --------------------

Corank of degree (0, 4): 715 ----- Value in Hilbert Biseries: 715

-------------------- d=(5, 0) --------------------

Corank of degree (5, 0): 2002 ----- Value in Hilbert Biseries: 2002

-------------------- d=(4, 1) --------------------

[TIMER] Temps pour add_lines : 19.5812 s
[TIMER] Temps pour Gauss (matrice 5565x7150) : 20.5796 s
number of reductions to 0 in degree (4, 1): 945 / 5565
Corank of degree (4, 1): 2530 ----- Value in Hilbert Biseries: 2530
[TIMER] Temps total pour bi-degré (4, 1) : 41.12 s

-------------------- d=(3, 2) --------------------

[TIMER] Temps pour add_lines : 75.8808 s
[TIMER] Temps pour Gauss (matrice 9450x12100) : 26.9756 s
number of reductions to 0 in degree (3, 2): 0 / 9450
Corank of degree (3, 2): 2650 ----- Value in Hilbert Biseries: 2650
[TIMER] Temps total pour bi-degré (3, 2) : 103.85 s

-------------------- d=(2, 3) --------------------

[TIMER] Temps pour add_lines : 75.7943 s
[TIMER] Temps pour Gauss (matrice 9450x12100) : 21.3546 s
number of reductions to 0 in degree (2, 3): 0 / 9450
Corank of degree (2, 3): 2650 ----- Value in Hilbert Biseries: 2650
[TIMER] Temps total pour bi-degré (2, 3) : 98.18 s

-------------------- d=(1, 4) --------------------

[TIMER] Temps pour add_lines : 12.6168 s
[TIMER] Temps pour Gauss (matrice 5565x7150) : 7.7504 s
number of reductions to 0 in degree (1, 4): 945 / 5565
Corank of degree (1, 4): 2530 ----- Value in Hilbert Biseries: 2530
[TIMER] Temps total pour bi-degré (1, 4) : 21.31 s

-------------------- d=(0, 5) --------------------

Corank of degree (0, 5): 2002 ----- Value in Hilbert Biseries: 2002

-------------------- d=(6, 0) --------------------

Corank of degree (6, 0): 5005 ----- Value in Hilbert Biseries: 5005

-------------------- d=(5, 1) --------------------

[TIMER] Temps pour add_lines : 944.1129 s
[TIMER] Temps pour Gauss (matrice 26355x20020) : 399.0586 s
number of reductions to 0 in degree (5, 1): 11342 / 26355
Corank of degree (5, 1): 5007 ----- Value in Hilbert Biseries: 5005
[TIMER] Temps total pour bi-degré (5, 1) : 1347.16 s

-------------------- d=(4, 2) --------------------

[TIMER] Temps pour add_lines : 4563.1160 s
[TIMER] Temps pour Gauss (matrice 41763x39325) : 576.5253 s
number of reductions to 0 in degree (4, 2): 7445 / 41763
Corank of degree (4, 2): 5007 ----- Value in Hilbert Biseries: 4675
[TIMER] Temps total pour bi-degré (4, 2) : 5143.68 s
"""