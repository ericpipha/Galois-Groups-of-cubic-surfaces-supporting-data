using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.instantiate()

using UnirationalDiscriminants
using Oscar

PATHOUTPUT = "/usr/people/pichonphara/discriminant_symmetric/outputs/overQQ/"
G = ARGS[1]

D = get_dict_T4(;K =QQ)

io = open(PATHOUTPUT*"discriminants_T4_$(G).txt","w")
io2 = open(PATHOUTPUT*"discriminants_T4_$(G)_short.txt","w")
pols = D[G]
println(G)
write(io,G*"\n")
write(io2,G*"\n")
if length(pols) == 1
    println("0-diml linear system")
elseif G == "Z1"
    println("trivial group")
else
    res = @timed prY1 = get_prY1(pols)
    println("computing prY1 took $(res.time) seconds")
    write(io,"computing prY1 took $(res.time) seconds \n")
    write(io2,"computing prY1 took $(res.time) seconds \n")
    println("its components have degree $(degree.(prY1))")
    write(io,"its components have degree $(degree.(prY1))\n")
    write(io2,"its components have degree $(degree.(prY1))\n")
    for p in prY1
        write(io,string(p)*"\n")
        write(io,"----\n")
    end
    multigrading = get_multigrading_T4(G)
    #multideg = get_multideg_disc_T1(G)
    res = @timed prY0 = get_prY0(pols)
    if prY0 != -1
        println("computing prY0 took $(res.time) seconds")
        write(io,"computing prY0 took $(res.time) seconds \n")
        write(io2,"computing prY0 took $(res.time) seconds \n")
        println("it has degree $(degree(prY0))")
        write(io,"it has degree $(degree(prY0)) \n")
        write(io2,"it has degree $(degree(prY0)) \n")
        println("certified: $(certify_prY0(prY0,pols))")
        write(io,"certified: $(certify_prY0(prY0,pols)) \n")
        write(io2,"certified: $(certify_prY0(prY0,pols)) \n")
        write(io,string(prY0)*"\n")
    else
        write(io,"prY0 is not a hypersurface \n")
        write(io2,"prY0 is not a hypersurface \n")
        println("it has dimension $(dim_prY0(pols))")
        write(io,"it has dimension $(dim_prY0(pols)) \n")
        write(io2,"it has dimension $(dim_prY0(pols)) \n")
    end
end
println("-------------------")
write(io,"------------------------------------------------------ \n")
write(io2,"------------------------------------------------------ \n")
close(io)
close(io2)
