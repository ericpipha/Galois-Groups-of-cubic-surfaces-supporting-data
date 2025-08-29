def map_from_gen(monomials, g):
    ''' Given a matrix `g` of size `n`, compute the induced action of `g` of a polynomial `P` in `n` variables.
    '''
    vs = monomials[0].parent().gens()
    coef_mat = matrix([[matrix_action(g, P).coefficient(m) for m in monomials] for P in monomials]).inverse().change_ring(ZZ)
    return coef_mat

def matrix_action(g, P):
    ''' Given a matrix `g` of size `n`, compute the induced action of `g` of a polynomial `P` in `n` variables.
    '''
    vs = vector(P.parent().gens())
    return P(list(g*vs))

def invariant_polynomials(group, monomials):
    ''' Given a matrix group `group` and a list of monomials `monomials`, computes the set of polynomials invariant under the action of the group.
    '''
    Lambda = identity_matrix(len(monomials)).image()
    for g in group:
        Lambda = Lambda.intersection((map_from_gen(monomials, g)-1).kernel())
    return list(Lambda.basis_matrix() * vector(monomials))

def S5_to_GL4(p):
    ''' Given a list `p` or distinct elements of `[0,1,2,3,4]`, returns the corresponding linear automorphism of QQ^4
    '''
    vecs = identity_matrix(4).rows() + [vector([-1,-1,-1,-1])]
    return matrix([vecs[p[i]] for i in range(4)]).transpose()

def find_curves(surface, dmax, g=0):
    ''' Given a cubic surface `surface` and a degree `dmax`, and a genus `g`, 
    returns a list of size `dmax+1` where the d-th entry is the list of homology classes
    corresponding to smooth curves of genus g and degree d. Only guaranteed to work for d=1, g=0.
    '''
    IP = surface.intersection_product
    hyperplane_class = surface.hyperplane_class
    
    lat = IntegralLattice(IP)
    lat.signature_pair()
    r = IP.nrows()
    h2 = hyperplane_class*IP*hyperplane_class
    DK = (surface.degree-4)
    a,b,c = h2^2*(2-2*g), h2^2*DK, h2

    Pic0 = lat.orthogonal_complement([hyperplane_class]).basis()
    Pic4plusH = matrix([h2*v for v in identity_matrix(r).rows()]+[hyperplane_class]).image().basis()
    Lambda = lat.sublattice(Pic0).intersection(lat.sublattice(Pic4plusH)).basis()
    QF = QuadraticForm(-matrix(Lambda) * IP * matrix(Lambda).transpose())
    S = QF.short_vector_list_up_to_length((a+b*dmax+c*dmax^2)+1)

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
    return N

def get_permutation(vs, Ms):
    """ Given a list of vectors `vs` and a list of matrices `Ms` acting on `vs` by permutation, returns the corresponding permutations.
    """
    permutations = []
    for M in Ms:
        permutations += [Permutation([vs.index(M*v)+1 for v in vs])]
    return permutations

def matrix_action_on_cohomology(g, fibre):
    """ Given a variety `fibre` defined by an equation f and a linear transformation `g` leaving f invariant, computes the induced action on the homology of `fibre`
    """
    coho = Cohomology(fibre.P)
    vs = vector(fibre.P.parent().gens())
    M = matrix([coho.coordinates(matrix_action(g,w)) for w in fibre.cohomology])
    e = g.det()
    M = block_diagonal_matrix([e*M, identity_matrix(1)])
    return (fibre.period_matrix^-1 * M * fibre.period_matrix).change_ring(ZZ)