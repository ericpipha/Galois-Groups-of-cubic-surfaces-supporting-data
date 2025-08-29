using UnirationalDiscriminants
using Oscar

# We start by illustrating some functionalities of the code.

# Here is a linear system I made up. 
R, vrs = polynomial_ring(QQ,["x";"y";"z";"w"])
x,y,z,w = vrs
h = [x^3 + y^3 + z^3 + w^3; x^2*y + z^2*w; x*y*z]

# We compute the discriminant in 3 ways.

# Use naive elimination from the critical equations, with saturation:
@time naive_elimination(h)
# Avoid saturation by doing the computation once on each chart of projective space:
@time naive_elimination_no_saturation(h)
# Now split the discriminant up into its two contributions, like in the text:
res1 = @timed prY1 = get_prY1(h)
res2 = @timed prY0 = get_prY0(h)
[res1.value;res2.value]

# Access all linear systems from Table 1
T1 = get_dict_T1()
# Here is how to ask for a basis for a specific group: 
h = T1["S3"]
# The following command returns a list of all groups in this table: 
all_groupnames("T1")
# You can do this over any field 
T1 = get_dict_T1(;K=GF(1993))
h = T1["S3"]
# Discriminants are computed over the ground field you specified 
get_prY1(h)

# There are similar commands for T3, T4. 
T4 = get_dict_T4()
h = T4["S5"]
D = naive_elimination_no_saturation(h)
factor(D)

# There are also tables for the restrictions to generic lines. 
T4_lines = get_dict_T4_lines(; K = GF(11)) # strangely, this may take a minute 
Gs = all_groupnames("T4_lines")
for G in Gs
    println(G)
    h = T4_lines[G]
    res = @timed disc = naive_elimination_no_saturation(h)
    MP = minimal_primes(ideal(disc))
    println("degrees: $(degree.(MP))")
end

#############################################

# The code below was used to generate the txt files. 

io = open("discriminants_T1.txt","w")
io2 = open("discriminants_T1_short.txt","w")
for G in all_groupnames("T1")
    pols = get_basis(G,"T1"; K=QQ)
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
        multigrading = get_multigrading_T1(G)
        multideg = get_multideg_disc_T1(G)
        res = @timed prY0 = get_prY0(pols; multigrading=multigrading,multideg=multideg)
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
end
close(io)
close(io2)

io = open("discriminants_T1_lines.txt","w")
io2 = open("discriminants_T1_lines_short.txt","w")
for G in all_groupnames("T1_lines")
    pols = get_basis(G,"T1_lines"; K=QQ)
    println(G)
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
            write(io,string(gens(p)[1])*"\n")
            write(io,"----\n")
        end
    end
    println("-------------------")
    write(io,"------------------------------------------------------ \n")
    write(io2,"------------------------------------------------------ \n")
end
close(io)
close(io2)


io = open("discriminants_T3.txt","w")
io2 = open("discriminants_T3_short.txt","w")
for G in all_groupnames("T3")
    pols = get_basis(G,"T3"; K=QQ)
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
        multigrading = get_multigrading_T3(G)
        #multideg = get_multideg_disc_T1(G)
        res = @timed prY0 = get_prY0(pols; multigrading=multigrading)
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
end
close(io)
close(io2)


io = open("discriminants_T3_lines.txt","w")
io2 = open("discriminants_T3_lines_short.txt","w")
for G in all_groupnames("T3_lines")
    pols = get_basis(G,"T3_lines"; K=QQ)
    println(G)
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
end
close(io)
close(io2)

io = open("discriminants_T4.txt","w")
io2 = open("discriminants_T4_short.txt","w")
D = get_dict_T4(;K = GF(1993))
for G in all_groupnames("T4")
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
        res = @timed prY0 = get_prY0(pols; multigrading=multigrading)
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
end
close(io)
close(io2)

io = open("discriminants_T4_lines_GF(105929).txt","w")
io2 = open("discriminants_T4_lines_short_GF(105929).txt","w")
D = get_dict_T4_lines(;K = GF(105929))
for G in all_groupnames("T4_lines")
    pols = D[G]
    println(G)
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
end
close(io)
close(io2)
