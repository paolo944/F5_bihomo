from sage.all import *
import time

###############################################################################
# 1. GÉNÉRATION DU SYSTÈME (Anciennement base_F5.sage)
###############################################################################

def generate_bilinear_system(n_x, n_y, m):
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

    return R, polynomials

def count_monomials_bidegree(n_vars, degree):
    """ Calcule (n+d-1) choose d """
    if degree < 0: return 0
    return binomial(n_vars + degree - 1, degree)

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

def get_matrix_dimensions(F, nx, ny, target_d):
    """
    Calcule à l'avance nrows et ncols pour le bidegré target_d=(d1, d2)
    """
    from sage.rings.polynomial.polydict import ETuple
    d1, d2 = target_d
    # Colonnes : tous les monômes de bidegré (d1, d2)
    ncols = count_monomials_bidegree(nx, d1) * count_monomials_bidegree(ny, d2)
    
    bi_hilbert = hilbert_biseries(nx, ny, len(F)).monomial_coefficients()
    
    h = bi_hilbert[ETuple([d1, d2])]
    nrows = int(ncols) - h
            
    return nrows, int(ncols)

###############################################################################
# 3. ALGORITHME DE CALCUL (Structure F5 simplifiée)
###############################################################################

def solve_bihomo_fast(R, F, dmax):
    nx = len([v for v in R.gens() if str(v).startswith('x')])
    ny = len([v for v in R.gens() if str(v).startswith('y')])

    # 1. Calcul des dimensions et pré-allocation
    nrows, ncols = get_matrix_dimensions(F, nx, ny, (dmax, 1))
    if nrows == 0: return None
    
    print(f"Allocation Matrix [{nrows}x{ncols}] pour bidegré ({dmax},1)...")
    # On crée la matrice directement en mémoire GF(2)
    M = matrix(GF(2), nrows, ncols)

    # 2. Génération des monômes de colonnes (Triés pour l'ordre monoïdal)
    Lx = [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(dmax, nx)]
    Ly = [R.monomial(*([0]*nx + list(v))) for v in IntegerVectors(1, ny)]
    cols = sorted([mx * my for mx in Lx for my in Ly], reverse=True)
    monom_to_idx = {m: i for i, m in enumerate(cols)}

    # 3. Remplissage direct de la matrice
    current_row = 0
    for f in F:
        du_x, du_y = dmax - 1, 0
        
        if du_x < 0 or du_y < 0: continue
        
        # Multiplicateurs u
        Lx_u = [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(du_x, nx)]
        Ly_u = [R.monomial(*([0]*nx + list(v))) for v in IntegerVectors(du_y, ny)]
        multipliers = [mx * my for mx in Lx_u for my in Ly_u]
        
        for u in multipliers:
            poly = u * f
            # Optimisation : on remplit la ligne i de la matrice M
            for coeff, monom in poly:
                if monom in monom_to_idx:
                    # M[row, col] = val est efficace sur GF(2) dans Sage
                    M[current_row, monom_to_idx[monom]] = coeff
            current_row += 1

    # 4. Élimination de Gauss (M4RI automatique)
    print(f"Lancement de M4RI...")
    t_start = time.time()
    M.echelonize()
    rank = M.rank()
    print(f"Rang: {rank} (Pleine: {rank == min(nrows, ncols)}) nblignes = {nrows} ncols = {ncols}")
    t_end = time.time()
    
    print(f"Terminé en {t_end - t_start:.4f}s")
    return M

###############################################################################
# EXÉCUTION
###############################################################################

if __name__ == "__main__":
    from sage.rings.polynomial.polydict import ETuple    
    
    NX, NY, M = 9, 9, 19
    R, F = generate_bilinear_system(NX, NY, M)
    
    bi_hilbert = hilbert_biseries(NX, NY, M).monomial_coefficients()
    dmax = 0
    for i in range(NX + 2):
        if bi_hilbert[ETuple([i, 1])] <= 0:
            dmax = i
            break

    print(f"dmax = {dmax}")

    print(f"Système généré: {M} équations, {NX} variables x, {NY} variables y")
    
    # On teste sur quelques bidegrés
    solve_bihomo_fast(R, F, dmax)