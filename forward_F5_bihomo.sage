import time
import os
import psutil
import cProfile

p = 2

def generate_bilinear_system(n_x, n_y, m):
    """
    Génère m polynômes bilinéaires homogènes de bi-degré (1,1)
    en n_x variables x_i et n_y variables y_j sur GF(2),
    avec une solution plantée qui annule tous les polynômes.
    """
    print(f"generating random {m} polynomials on GF({p}) with {n_x} variables in\
    x and {n_y} variables in y")

    K = GF(p)

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
    
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) + 10)
    denom = ((1 - tx)^(nx)) * ((1 - ty)^(ny))
    num = (1 - tx*ty)^m + Nm(ny, m, tx, ty) + Nm(nx, m, ty, tx)
    return num / denom

def series_briaud(nx, ny, m):
    R.<tx,ty> = PowerSeriesRing(ZZ, default_prec=max(nx, ny) + 10)
    Ad = ((1 - tx*ty)**(m)) / ((1 - tx)^(nx) * (1 - ty)^(ny))
    # first_term = 1/((1 - tx)^(nx) * (1 - ty)^(ny))
    second_num = tx * (1 - ty)^(nx) - ty * (1 - tx)^(ny)
    second_denom = (tx - ty) * (1 - tx)^(nx) * (1 - ty)^(ny)
    Bd = second_num/second_denom
    return (Ad, Bd)

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
    print(f"h: {h}, nrows: {nrows}")
            
    return nrows, int(ncols)

def get_matrix_stats(M):
    nb_lignes_nz = M.rank()

    cols_actives = set()
    for i in range(nb_lignes_nz):
        cols_actives.update(M[i].nonzero_positions())
    
    nb_cols_nz = len(cols_actives)

    difference = nb_cols_nz - nb_lignes_nz

    return {
        "lignes_non_nulles": nb_lignes_nz,
        "colonnes_non_nulles": nb_cols_nz,
        "difference": difference
    }

def timing_rank(M):
    return M.rank('m4ri')

def construct_matrix(M, multipliers, monom_to_idx):
    current_row = 0
    for u in multipliers:
        for f in F:
            poly = u * f
            for coeff, monom in poly:
                try:
                    M[current_row, monom_to_idx[monom]] = coeff
                except IndexError:
                    print(current_row)
                    return
            current_row += 1

def generate_monomials(R, nx, ny, dx, dy):
    Lx = [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(dx, nx)]
    Ly = [R.monomial(*([0]*nx + list(v))) for v in IntegerVectors(dy, ny)]
    cols = sorted([mx * my for mx in Lx for my in Ly], reverse=True)
    monom_to_idx = {m: i for i, m in enumerate(cols)}
    return monom_to_idx

def get_multipliers(R, dmax, nx, ny):
    return [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(dmax - 1, nx)]

def solve_bihomo_fast(R, F, dmax):
    nx = len([v for v in R.gens() if str(v).startswith('x')])
    ny = len([v for v in R.gens() if str(v).startswith('y')])

    # 1. Calcul des dimensions et pré-allocation
    nrows, ncols = get_matrix_dimensions(F, nx, ny, (dmax, 1))
    if nrows == 0: return None
    
    print(f"Allocation Matrix [{nrows}x{ncols}] pour bidegré ({dmax},1)...")
    # On crée la matrice directement en mémoire GF(2)
    M = matrix(R.base_ring(), nrows, ncols)

    monom_to_idx = generate_monomials(R, nx, ny, dmax, 1)
    
    multipliers = get_multipliers(R, dmax, nx, ny)

    construct_matrix(M, multipliers, monom_to_idx)

    #####Ecriture format pour SPASM#####
    # A = matrix(M, sparse=True)
    # A.save(f"test_matrix_{dmax}.sobj")
    # out = open(f"Matrices/A_{dmax}.sms", "w")
    # out.write("{0} {1} M\n".format(nrows, ncols))
    # for (i,j) in A.nonzero_positions():
    #     out.write("{0} {1} {2}\n".format(i + 1, j + 1, A[i,j]))
    # out.write("0 0 0\n")
    # out.close()

    print(f"Lancement de l'algèbre linéaire...")
    t_start = time.time()
    # print(f"density before rref {M.density(approx=True)}")
    rank = timing_rank(M)
    print(f"Rang: {rank} (Pleine: {rank == min(nrows, ncols)}) nblignes = {nrows} ncols = {ncols}")
    t_end = time.time()
    
    print(f"Terminé en {t_end - t_start:.4f}s")

    # p = M.visualize_structure(maxsize=1024)
    # p.save('/tmp/matrix_echelon.png')

    return M, rank == min(nrows, ncols)


def solve_bihomo_sparse(R, F, dmax):
    nx = len([v for v in R.gens() if str(v).startswith('x')])
    ny = len([v for v in R.gens() if str(v).startswith('y')])

    nrows, ncols = get_matrix_dimensions(F, nx, ny, (dmax, 1))
    if nrows == 0: 
        return None, False
    
    Lx = [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(dmax, nx)]
    Ly = [R.monomial(*([0]*nx + list(v))) for v in IntegerVectors(1, ny)]
    cols = sorted([mx * my for mx in Lx for my in Ly], reverse=True)
    monom_to_idx = {m: i for i, m in enumerate(cols)}

    entries = {}
    current_row = 0
    
    Lx_u = [R.monomial(*(list(v) + [0]*ny)) for v in IntegerVectors(dmax - 1, nx)]
    
    for f in F:
        for u in Lx_u:
            poly = u * f
            for coeff, monom in poly:
                if monom in monom_to_idx:
                    entries[(current_row, monom_to_idx[monom])] = 1
            current_row += 1

    M = matrix(GF(2), nrows, ncols, entries, sparse=True)

    t_start = time.time()
    
    density = M.density()
    print(f"Densité avant calcul du rang : {float(density):.6f}")
    
    rank = M.rank() 
    
    t_end = time.time()
    print(f"Terminé en {t_end - t_start:.4f}s")
    
    return M, rank == min(nrows, ncols)

if __name__ == "__main__":
    from sage.rings.polynomial.polydict import ETuple    
    
    NX = 8
    NY = NX
    M = NY + NX + 1
    R, F = generate_bilinear_system(NX, NY, M)
    # F = homogenized_ideal(load(f"../MPCitH_SBC/system/sage/system_bilin_{NX + NY}_{M}.sobj"))
    # print(f"SBC Nx = Ny = {NX} M = {M}")
    R = F[0].parent()
    
    F = homogenized_ideal(F)
    bi_hilbert = hilbert_biseries(NX, NY, M)
    bi_hilbert = bi_hilbert.monomial_coefficients()
    Ad, Bd = series_briaud(NX, NY, M)
    Ad = Ad.monomial_coefficients()
    Bd = Bd.monomial_coefficients()
    
    dmax = 0
    for i in range(NX + NY):
        if bi_hilbert[ETuple([i, 1])] <= 0:
            dmax = i
            break

    dmax2 = 0
    for i in range(1, NX + NY):
        if Ad[ETuple([i, 1])] <= 0 or Bd[ETuple([i, 1])] <= 0:
            dmax2 = i
            break
    
    for i in range(1, max(dmax2, dmax)):
        print(f"Ad[{i}, 1] = {Ad[ETuple([i, 1])]} Bd[{i}, 1] = {Bd[ETuple([i, 1])]} HBS[{i}, 1] = {bi_hilbert[ETuple([i, 1])]}")

    estimation = ceil((NX*NY - NY)/(M - NY))

    print(f"dmax = {dmax} estimation = {estimation} briaud/romaric predicted = {dmax2}")

    # print(f"Système généré: {M} équations, {NX} variables x, {NY} variables y")
    
    # On teste sur quelques bidegrés

    F = sorted(F, key=lambda f: f.lm())

    # for i in range(2, 10):
    #     print(f"d = {i}")
    #     M, rank = solve_bihomo_fast(R, F, i)
    #     if not rank:
    #         print(f"dreg = {i+1}")
    #         break

    cProfile.run("solve_bihomo_fast(R, F, 7)", sort='cumtime')

    M.save("test_matrix.sobj")

process = psutil.Process(os.getpid())
mem_bytes = process.memory_info().rss
mem_mb = mem_bytes / (1024 ** 2)
mem_gb = mem_bytes / (1024 ** 3)
if mem_gb >= 1.0:
    print(f"Mémoire RAM max utilisée : {float(mem_gb):.2f} Go")
else:
    print(f"Mémoire RAM max utilisée : {float(mem_mb):.2f} Mo")