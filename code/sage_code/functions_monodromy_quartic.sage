def find_curves(surface, dmax, NS=None, g=0):
    IP = surface.intersection_product
    hyperplane_class = surface.hyperplane_class

    if NS == None:
        NS = IntegerRelations(surface.holomorphic_period_matrix.transpose()).basis
    IP = NS*IP*NS.transpose()
    hyperplane_class = NS.solve_left(hyperplane_class)
    
    lat = IntegralLattice(IP)
    lat.signature_pair()
    r = IP.nrows()
    h2 = hyperplane_class*IP*hyperplane_class
    DK = (surface.degree-4)
    a,b,c = h2^2*(2-2*g), h2^2*DK, h2

    Pic0 = lat.orthogonal_complement([hyperplane_class]).basis()
    Pic4plusH = matrix(ZZ, [h2*v for v in identity_matrix(r).rows()]+[hyperplane_class]).image().basis()
    Lambda = lat.sublattice(Pic0).intersection(lat.sublattice(Pic4plusH)).basis()

    QF = QuadraticForm(-matrix(Lambda) * IP * matrix(Lambda).transpose())
    S = QF.short_vector_list_up_to_length((a+b*dmax+c*dmax^2)/2+1)

    M = [[]] + [S[(a+b*d+c*d^2)/2] for d in range(1,dmax+1) if (a+b*d+c*d^2)/2 in ZZ]
    for d in range(1, dmax+1):
        classes = [(v*matrix(Lambda)+d*hyperplane_class)/h2 for v in M[d]]
        classes = [v for v in classes if v in identity_matrix(ZZ,r).image()]
        M[d] = classes
    N = [[]]
    for d in range(1, dmax+1):
        add = M[d].copy()
        for d2 in range(1, d):
            for v2 in N[d2]:
                add = [v for v in add if v*IP*v2>=0]
        N += [add]
    N = [[v*NS for v in add] for add in N]
    return N
    
def get_ini_conds(w, fam, basepoint):
    GM, denom = fam.gaussmanin()
    L = fam.picard_fuchs_equation(w)
    Dt = L.parent().gen(0)
        
    n = L.order()
    derivatives = [w]
    for k in range(n-1):
        derivatives += [denom*derivatives[-1].derivative(t) +  derivatives[-1]*GM]
    derivatives_at_basepoint2 = [w(basepoint) for w in derivatives]
    integration_correction = diagonal_matrix([1/ZZ(factorial(k)) for k in range(n)])
    
    matders = []
    for i in range(n):
        matders += [[c(basepoint) for c in ((denom * Dt)^i).coefficients(sparse=False)]+[0]*(n-1-i)]
    matders = matrix(matders)
    
    derivatives_at_basepoint = matders.inverse() * matrix(derivatives_at_basepoint2).change_ring(QQ)
    initial_conditions = integration_correction * derivatives_at_basepoint

    return initial_conditions

def pick_indices(Ls, inicondss):
    indices = [i for i in range(len(Ls))]
    indices.sort(key=lambda i: (Ls[i].order(), Ls[i].degree()))
    kept_indices = []
    for i in indices:
        CB = block_matrix([[inicond.transpose() for i, inicond in enumerate(inicondss) if i in kept_indices]])
        if any([w not in CB.transpose().image() for w in inicondss[i].rows()]):
            kept_indices += [i]
    for j in range(len(kept_indices)-1,-1,-1): # we remove low orders covered by higher
        i = kept_indices[j]
        CB = block_matrix([[inicond.transpose() for k, inicond in enumerate(inicondss) if k in kept_indices and k!=i]])
        if all([w in CB.transpose().image() for w in inicondss[i].rows()]):
            kept_indices.remove(i)
    return kept_indices

def map_from_gen(monomials, g):
    vs = monomials[0].parent().gens()
    coef_mat = matrix([[matrix_action(g, P).coefficient(m) for m in monomials] for P in monomials]).inverse().change_ring(ZZ)
    return coef_mat
def matrix_action(g, P):
    vs = vector(P.parent().gens())
    return P(list(g*vs))
def invariant_polynomials(group, monomials):
    Lambda = identity_matrix(len(monomials)).image()
    for g in group:
        Lambda = Lambda.intersection((map_from_gen(monomials, g)-1).kernel())
    return list(Lambda.basis_matrix() * vector(monomials))
def S5_to_GL4(p):
    vecs = identity_matrix(4).rows() + [vector([-1,-1,-1,-1])]
    return matrix([vecs[p[i]] for i in range(4)]).transpose()

def matrix_action_on_cohomology(g, fibre):
    coho = Cohomology(fibre.P)
    vs = vector(fibre.P.parent().gens())
    M = matrix([coho.coordinates(matrix_action(g,w)) for w in fibre.cohomology_internal])
    e = g.det()
    M = block_diagonal_matrix([e*M, identity_matrix(1)])
    return (fibre.period_matrix^-1 * M * fibre.period_matrix).change_ring(ZZ)