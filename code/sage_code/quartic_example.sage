os.environ["SAGE_NUM_THREADS"] = '10'
from lefschetz_family.fibration import Fibration
from lefschetz_family.fibration import Hypersurface
from lefschetz_family.numperiods.integerRelations import IntegerRelations
from lefschetz_family.util import Util
from lefschetz_family.numperiods.family import Family
import os.path

import logging
logging.basicConfig()
logging.getLogger('lefschetz_family.hypersurface').setLevel(logging.INFO)
logging.getLogger('lefschetz_family.fibration').setLevel(logging.INFO)

load("functions_monodromy_quartic.sage")

R.<x,y,z,w> = QQ[]
S.<t> = R[]
monomials = R.monomials_of_degree(4)

invariant_pols = {}
i=1
for key, G in crystallographic_groups.items():
    print(i,"/",len(crystallographic_groups),":",key, " "*10, end="\r")
    i+=1
    invariant_pols[key] = invariant_polynomials(G.gens(), monomials)
save(invariant_pols, "data_crystallographic/invariant_pols")

crystallographic_groups = load("crystallographic_groups")
invariant_pols = load("invariant_pols")

generically_singular_families = []
k=1
for key, span in invariant_pols.items():
    print(k,"/",len(crystallographic_groups),":",key, end="\r")
    k+=1
    smooth = False
    for i in range(100):
        if smooth:
            break
        P = vector([randint(-100,100) for i in range(len(span))]) * vector(span)
        if P.jacobian_ideal().dimension()==0:
            smooth = True
    if not smooth:
        generically_singular_families += [key]
len(generically_singular_families)

possible_keys = [key for key in list(invariant_pols.keys()) if key not in generically_singular_families]

seen = []
pols = []
for key, ip in invariant_pols.items():
    if ip in pols:
        seen[pols.index(ip)] += [key]
    else:
        seen += [[key]]
        pols += [ip]
possible_keys = [key for key in possible_keys if key in [l[0] for l in seen]]

key = "32.9"
symmetrygroup = crystallographic_groups[key]
span = invariant_pols[key]

vals = 50, -89, 12, -95, 77, 45, 96, 99, 81, 79, 99, 44, 22, -20, -62, -75, 55, -65, -95, 14, -63, -76, -21, -9, 95, 48, 64, -56, -24, -40, 15, -31, 7, -62, -76
vals2 = 5,3,7,-9,13,-2,4, -6, -1, -9, -11, -3, -14, -13, 4, -14, 10, 2, -1, -10, -10, -9, -13, -14, -7, -8, 10, 7, -1, -9, -12, -12, -11

basepoint = 0
P0 = vector(vals[:len(span)])*vector(span)
Pinf = vector(vals2[:len(span)])*vector(span)
Pt = P0 + t*Pinf

# this can take quite long (expected to take one hour) -- we recommend saving the result!
couldnotopen = False
try:
    fibre = load("fibre_"+key)
except:
    couldnotopen = True
if couldnotopen or fibre.P != Pt(basepoint):
    fibration = [vector(ZZ, [10, -7, 6, 1]), vector(ZZ, [9, -6, -10, 6]), vector(ZZ, [-10, -2, -7, -7])]
    fibre = Hypersurface(Pt(basepoint), fibration=fibration, nbits=1000)
_ =  fibre.period_matrix

group_action = [matrix_action_on_cohomology(matrix(g), fibre) for g in symmetrygroup]
Picsublattice = matrix(ZZ, fibre.hyperplane_class).image()
for a in group_action:
    if ZZ((fibre.period_matrix * a*fibre.period_matrix.inverse())[0,0]) == 1:
        Picsublattice += ( (a-1).right_kernel_matrix() * fibre.intersection_product ).right_kernel() 
    else:
        Picsublattice += (a-1).right_kernel()
Picsublattice = Picsublattice.saturation()
generic_Picard_lattice = Picsublattice.basis_matrix()

conics = find_curves(fibre, 2, NS=generic_Picard_lattice)[2]

fam = Family(Pt, basepoint=basepoint)
GM, denom = fam.gaussmanin()
critical_values = []
for P, _ in denom.factor():
    critical_values += [P.roots(QQbar, multiplicities=False)]
critical_values = flatten(critical_values)

# We generate a list of candidate cyclic forms -- we are looking for forms with small Picard-Fuchs equations that generate the cohomology

group_action = [matrix_action_on_cohomology(matrix(g), fibre) for g in symmetrygroup.gens()]
NS = IntegerRelations(fibre.holomorphic_period_matrix.transpose()).basis
TrX = (NS * fibre.intersection_product).right_kernel_matrix()

candidates = [identity_matrix(22).row(0)]
for c in conics:
    ws = IntegerRelations(matrix(fibre.period_matrix * TrX.stack(c).transpose())).basis.rows()
    for w in ws:
        if w not in candidates and all([abs(ZZ(v))<=2 for v in w]) and w!=identity_matrix(22).row(21):
            candidates+=[w]
more_candidates = []
for c in candidates:
    restr_permat = (matrix([c] + [identity_matrix(22).row(0)]) * fibre.period_matrix)
    orthcomp = (IntegerRelations(restr_permat.transpose()).basis * fibre.intersection_product).right_kernel_matrix()
    ws = IntegerRelations(fibre.period_matrix * orthcomp.transpose()).basis[1:]
    for w in ws:
        if w not in candidates+more_candidates and all([abs(ZZ(v))<=2 for v in w]) and w!=identity_matrix(22).row(21):
            more_candidates+=[w]
candidates = candidates + more_candidates

# We pick out our cyclic forms among the candidates
iniconds = [get_ini_conds(w[:21].change_ring(fam.upolring), fam, basepoint) for w in candidates]
picard_fuchs_equations = [fam.picard_fuchs_equation(w[:21]) for w in candidates]
cyclic_forms = [candidates[i][:21].change_ring(fam.upolring) for i in pick_indices(picard_fuchs_equations, iniconds)]

fib = Fibration(Pt, family=fam, fibre=fibre, basepoint=basepoint, cyclic_forms=cyclic_forms, nbits=400)

%time fib.monodromy_matrices