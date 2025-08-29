os.environ["SAGE_NUM_THREADS"] = '10'

from lefschetz_family.fibration import Fibration
from lefschetz_family.numperiods.cohomology import Cohomology
from lefschetz_family import Hypersurface
import time
load("functions_monodromy_cubic.sage")

d=3
R.<x,y,z,w> = QQ[]
monomials = R.monomials_of_degree(3)
symmetric_subgroups = {
    'S4':[[2,3,4,1], [2,1,3,4]],
    'S3':[[2,3,1,4], [2,1,3,4]],
    'A4':[[2,3,1,4], [2,1,4,3]],
    'D4':[[3,4,2,1], [2,1,3,4]],
    'K4':[[2,1,4,3], [3,4,1,2]],
    'nK4':[[2,1,3,4], [1,2,4,3]],
    'Z4':[[3,4,2,1]],
    'Z3':[[2,3,1,4]],
    'Z2':[[2,1,3,4]],
    'DT':[[2,1,4,3]],
    'Z1':[[1,2,3,4]],
}
for key, l in symmetric_subgroups.items():
    symmetric_subgroups[key] = MatrixGroup([Permutation(perm).to_matrix() for perm in l])

subgroups_S5 = {
    "S5": [[3,0,1,4,2], [1,0,2,3,4]],
    "A5": [[3,0,1,4,2],[1,2,0,3,4]],
    "Z5": [[3,0,1,4,2]],
    "D5": [[3,0,1,4,2], [1,0,3,2,4]],
    "tS3": [[1,4,2,3,0], [1,0,3,2,4]],
    "F5": [[1,2,3,4,0],[0,2,4,1,3]],
    "D6": [[1,2,0,3,4], [1,0,2,3,4], [0,1,2,4,3]],
    "Z6": [[1,2,0,3,4], [0,1,2,4,3]]
}
for key, gens in subgroups_S5.items():
    symmetric_subgroups[key] = MatrixGroup([S5_to_GL4(g).transpose() for g in gens])

invariant_pols = {}
for key, G in symmetric_subgroups.items():
    invariant_pols[key] = invariant_polynomials(G.gens(), monomials)

invariant_pols_distinct = []
invariant_pols_indices = []
for key, l in invariant_pols.items():
    if len(l)<2: # either the list is empty or it's a point
        continue
    if l not in invariant_pols_distinct:
        invariant_pols_distinct += [l]
        invariant_pols_indices += [[key]]
    else:
        invariant_pols_indices[invariant_pols_distinct.index(l)] += [key]


# we compute the periods of the Clebsch cubic surface
print("Computing periods of Clebsch surface")
Pbase = invariant_pols["S5"][0]
fibration = [vector(ZZ, [4, -3, 8, 1]), vector(ZZ, [-1, 3, 4, -10]), vector(ZZ, [8, 4, -1, -8])]
fibre = Hypersurface(Pbase, nbits=1000, fibration=fibration)
_ = fibre.period_matrix

# we compute the monodromy matrices along a random line in L_{S_4}
print("Computing monodromy along line in S4-invariant cubics")
random_vector = vector([-41, -2, -95, 11, -58, -9, -94, -81, -93, 38, -86, -13, 41, 36, 62, -13, 82, -80, 31, 25, 55, -74, -51, -9, 27])
span = invariant_pols["S4"]
basepoint = 0
S.<t> = R[]
Pt = Pbase + (t-basepoint) * random_vector[:len(span)]*vector(span)
fib = Fibration(Pt, fibre=fibre, basepoint=basepoint, nbits=200)
fib.monodromy_matrices

# we recover the 27 lines
lines = find_curves(fibre, 1)[1]
assert len(lines)==27, "did not find 27 lines"

# and the permutation of the lines induced by monodromy
print("Monodromy on lines of S4-symmetric cubic surfaces:")
print([Permutation([lines.index(M*L) +1 for L in lines]) for M in fib.monodromy_matrices])

# we recover the Klein 4-group
print("Group id of Galois group of lines of S4-symmetric cubic surfaces:")
print(PermutationGroup([Permutation([lines.index(M*L) +1 for L in lines]) for M in fib.monodromy_matrices]).group_id())



# Now we compute all the Galois groups of symmetric families

print("\n\n Computing monodromy of all G-invariant families for G a subgroup of S5")

# this cell sorts the lines
action_matrices = [matrix_action_on_cohomology(matrix(g), fibre) for g in symmetric_subgroups["S4"].gens()]
S4_action = get_permutation(lines, action_matrices)
lattices = []
for o in PermutationGroup(S4_action).orbits():
    L = [lines[i-1] for i in o]
    fixed = [L[0]]
    while len(fixed) < len(L):
        L.sort(key = lambda l: matrix(fixed) * fibre.intersection_product * l)
        fixed += [L[len(fixed)]]
    lattices += [fixed]
for L in lattices:
    if len(L) ==3:
        triple = L
    elif matrix(L[:6]) * fibre.intersection_product * matrix(L[:6]).transpose() == -1:
        o1 = L
    else:
        o2 = L
CB = o1[:6]
CB = [((fibre.hyperplane_class+sum(CB))/3).change_ring(ZZ)] + CB
CB = matrix(CB).transpose()
o2.sort(key = lambda l: CB^-1 * l)
lines = o1 + o2 + triple
o1[4], o1[5] = o1[5], o1[4]
CB = o1[:6]
CB = [((fibre.hyperplane_class+sum(CB))/3).change_ring(ZZ)] + CB
CB = matrix(CB).transpose()
o11 = o1[:6]
o12 = o1[6:]
CB = o11
CB = [((fibre.hyperplane_class+sum(CB))/3).change_ring(ZZ)] + CB
CB = matrix(CB).transpose()
o12.sort(key=lambda l: -CB^-1*l)
o1 = o11+o12
lines = o1 + o2 + triple
lines[24], lines[25], lines[26] = lines[26], lines[24], lines[25]


monodromy_action = [None for i in range(len(invariant_pols_distinct))]
parametrised_lines = [None for i in range(len(invariant_pols_distinct))]
group_action = {}


for key, group in symmetric_subgroups.items():
    action_matrices = [matrix_action_on_cohomology(matrix(g), fibre) for g in group.gens()]
    group_action[key] = get_permutation(lines, action_matrices)

S.<t> = R[]
for i, span in enumerate(invariant_pols_distinct):
    if monodromy_action[i]!=None:
        continue
    try:
        print(len([l for l in monodromy_action if l!=None])+1,"/",len(invariant_pols_distinct),":",invariant_pols_indices[i])
        print(span)
        basepoint = 0
        # Pt = pick_random_line(S, span, Pbase=Pbase)
        Pt = Pbase + (t-basepoint) * random_vector[:len(span)]*vector(span)

        # Compute the group action on the lines
        for key in invariant_pols_indices[i]:
            action_matrices = [matrix_action_on_cohomology(matrix(g), fibre) for g in symmetric_subgroups[key].gens()]
            group_action[key] = get_permutation(lines, action_matrices)
            
        fib = Fibration(Pt, fibre=fibre, basepoint=basepoint, nbits=400)
        begin = time.time()
        monodromy_matrices = fib.monodromy_matrices
        end = time.time()
        duration_str = time.strftime("%H:%M:%S",time.gmtime(end-begin))
        print("computed monodromy in", duration_str)
        
        monodromy_action[i] = get_permutation(lines, monodromy_matrices)
        parametrised_lines[i] = Pt
        print("\n")
        # break
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except Exception as e:
        print("failed " + str(i))
        print(e, "\n")
        continue

for key in symmetric_subgroups.keys():
    if key not in group_action.keys():
        action_matrices = [matrix_action_on_cohomology(matrix(g), fibre) for g in symmetric_subgroups[key].gens()]
        group_action[key] = get_permutation(lines, action_matrices)