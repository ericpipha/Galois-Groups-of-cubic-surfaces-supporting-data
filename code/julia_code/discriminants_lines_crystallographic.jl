using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.instantiate()

using UnirationalDiscriminants
using Oscar

PATHOUTPUT = "./"
G = ARGS[1]

D = get_dict_T8_lines(;K = QQ)

io = open(PATHOUTPUT*"discriminants_$(G)_T8_lines_QQ.txt","w")
io2 = open(PATHOUTPUT*"discriminants_$(G)_T8_lines_short_QQ.txt","w")
pols = D[G]
println(G)
println(PATHOUTPUT*"discriminants_T8_lines_QQ.txt")
println(PATHOUTPUT*"discriminants_T8_lines_short_QQ.txt")

write(io,G*"\n")
write(io2,G*"\n")
if length(pols) == 1
    println("0-diml linear system")
else
    res = @timed disc = naive_elimination_no_saturation(pols)
    MP = minimal_primes(ideal(disc))
    println("computing disc took $(res.time) seconds")
    write(io,"computing disc took $(res.time) seconds \n")
    write(io2,"computing disc took $(res.time) seconds \n")
    println("its components have degree $(degree.(MP))")
    write(io,"its components have degree $(degree.(MP))\n")
    write(io2,"its components have degree $(degree.(MP))\n")
    for p in MP
        if dim(p) == length(pols)-1
            println("disc = "*string(gens(p)[1]))
            write(io,string(gens(p)[1])*"\n")
            write(io,"----\n")
        else
            println("dim = "*string(dim(p))*" instead of "*string(length(pols)-1))
        end
    end
end
println("-------------------")
write(io,"------------------------------------------------------ \n")
write(io2,"------------------------------------------------------ \n")
