libgap.LoadPackage("CrystCat")
groups = {}
nrsystems = int(str(gap("NrCrystalSystems( 4 )")))
for system in range(1,nrsystems+1):
    nrqclass = int(str(gap("NrQClassesCrystalSystem(4,"+str(system)+")")))
    for qclass in range(1,nrqclass+1):
        nrzclass = int(str(gap("NrZClassesQClass(4,"+str(system)+","+str(qclass)+")")))
        nrzclass = 1;
        for zclass in range(1, nrzclass+1):
            res = gap("GeneratorsOfGroup(MatGroupZClass( 4,"+str(system)+","+str(qclass)+","+str(zclass)+ "))")
            gens = [matrix(M) for M in eval(str(res))]
            key = str(system)+"."+str(qclass)
            groups[key] = MatrixGroup(gens if gens!=[] else [identity_matrix(4)])
save(groups, "crystallographic_groups")