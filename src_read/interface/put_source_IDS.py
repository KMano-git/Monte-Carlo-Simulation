# impor imas mod

try:
    import imas
except Exception as e:
    print("Required IMAS support library not available on this system.")

# import 
from dmyIDS import impMome # implity momentum data

# ======================================================================
# store source data for IDS
# 1st argument : IDS database
# 2nd argument : source datas
# 3rd argument : data size[X, Y]
# 4th argument : element data
# ======================================================================
def es2IDS(imas_obj, data, ds, elm):
    [wssn, wssp, wswi] = data[:3]
    [twci1, twci2, twrd1, twrd2] = data[3:7]
    [tionZ01, tionZ1, tionZ02, tionZ2, trecZ01, trecZ1, trecZ02, trecZ2] = data[7:15] # x[0]:0,x[1]:1~nzmx
    [tfrz01, tfrz1, tfrz02, tfrz2, tthz01, tthz1, tthz02, tthz2, tdnz01, tdnz1, tdnz02, tdnz2] = data[15:] # x[0]:0,x[1]:1~nzmx
    
    [nx, ny] = ds
    dsize = nx * ny

    [ename, nzmx] = [[elm[0][0], elm[1][0]], [elm[0][1], elm[1][1]]]

    tdnz = [tdnz01, tdnz1, tdnz02, tdnz2]
    tfrz = [tfrz01, tfrz1, tfrz02, tfrz2]
    tthz = [tthz01, tthz1, tthz02, tthz2]
    
    [[imome01, imome02], [imome1, imome2]] = impMome([tdnz, tfrz, tthz], nzmx)

    ntl_size = 3 + 3 + 3 # ntl(3), Ar(3), C(3)
    ion_size = (nzmx[0]+nzmx[1])*3 + 2 # (Ar + C)*(tionZ, irecZ, imome) + twci

    ntlname = ["wssn", "wswi", "wssp"]
    for i in range(len(ename)):
        ntlname.append("tionZ_" + ename[i])
        ntlname.append("trecZ_" + ename[i])
        ntlname.append("imome_" + ename[i])

    impname = []
    for i in range(len(ename)):
        impname.append([])
        impname[i].append("tionZ_" + ename[i] + "_")
        impname[i].append("trecZ_" + ename[i] + "_")
        impname[i].append("imome_" + ename[i] + "_")
        impname[i].append("twci_"  + ename[i])

    #--------------------------------------
    # edge_sources
    imas_obj.edge_sources.grid_ggd.resize(1)
    imas_obj.edge_sources.time.resize(1)
    imas_obj.edge_sources.time[0] = 1
    imas_obj.edge_sources.ids_properties.homogeneous_time = 1
    imas_obj.edge_sources.source.resize(1)
    imas_obj.edge_sources.source[0].ggd.resize(1)
    imas_obj.edge_sources.source[0].ggd[0].neutral.resize(ntl_size) # n1~3:ntl, n4~7:im1, n8~11:im2
    imas_obj.edge_sources.source[0].ggd[0].ion.resize(ion_size) # i1~55:im1, i56~74:im2
 
    nind = 0
    # neutral particles(wssn)
    esnp0 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnp0.particles.resize(1)
    esnp0.grid_subset_index = 1
    esnp0.particles[0].values.resize(dsize)
    esnp0.particles[0].coefficients.resize([dsize, 2])
    esnp0.label = ntlname[nind]
    nind += 1

    # neutral energy(wswi)
    esne0 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esne0.energy.resize(1)
    esne0.grid_subset_indes = 1
    esne0.energy[0].values.resize(dsize)
    esne0.energy[0].coefficients.resize([dsize, 2])
    esne0.label = ntlname[nind]
    nind += 1

    # neutral momentum(wssp)
    esnm0 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnm0.momentum.resize(1)
    esnm0.grid_subset_indes = 1
    esnm0.momentum[0].radial.resize(dsize)
    esnm0.momentum[0].radial_coefficients.resize([dsize, 2])
    esnm0.label = ntlname[nind]
    nind += 1

    # neutral particles(tionZ01)
    esnp011 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnp011.particles.resize(1)
    esnp011.grid_subset_index = 1
    esnp011.particles[0].values.resize(dsize)
    esnp011.particles[0].coefficients.resize([dsize, 2])
    esnp011.label = ntlname[nind]
    nind += 1

    # neutral particles(trecZ01)
    esnp012 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnp012.particles.resize(1)
    esnp012.grid_subset_index = 1
    esnp012.particles[0].values.resize(dsize)
    esnp012.particles[0].coefficients.resize([dsize, 2])
    esnp012.label = ntlname[nind]
    nind += 1

    # ion0 momentum(imome0[0])
    esnm1 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnm1.momentum.resize(1)
    esnm1.grid_subset_indes = 1
    esnm1.momentum[0].radial.resize(dsize)
    esnm1.momentum[0].radial_coefficients.resize([dsize, 2])
    esnm1.label = ntlname[nind]
    nind += 1

    # neutral particles(tionZ02)
    esnp021 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnp021.particles.resize(1)
    esnp021.grid_subset_index = 1
    esnp021.particles[0].values.resize(dsize)
    esnp021.particles[0].coefficients.resize([dsize, 2])
    esnp021.label = ntlname[nind]
    nind += 1

    # neutral particles(trecZ02)
    esnp022 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnp022.particles.resize(1)
    esnp022.grid_subset_index = 1
    esnp022.particles[0].values.resize(dsize)
    esnp022.particles[0].coefficients.resize([dsize, 2])
    esnp022.label = ntlname[nind]
    nind += 1

    # ion0 momentum(imome[1])
    esnm2 = imas_obj.edge_sources.source[0].ggd[0].neutral[nind]
    esnm2.momentum.resize(1)
    esnm2.grid_subset_indes = 1
    esnm2.momentum[0].radial.resize(dsize)
    esnm2.momentum[0].radial_coefficients.resize([dsize, 2])
    esnm2.label = ntlname[nind]
    nind += 1

    for iind in range(nzmx[0]):
        # ion particles(tionZ1)
        esip11 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 0 + iind]
        esip11.particles.resize(1)
        esip11.grid_subset_index = 1
        esip11.particles[0].values.resize(dsize)
        esip11.particles[0].coefficients.resize([dsize, 2])
        esip11.label = impname[0][0] + str("{:>02d}".format(iind + 1))

        # ion particles(trecZ1)
        esip12 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 1 + iind]
        esip12.particles.resize(1)
        esip12.grid_subset_index = 1
        esip12.particles[0].values.resize(dsize)
        esip12.particles[0].coefficients.resize([dsize, 2])
        esip12.label = impname[0][1] + str("{:>02d}".format(iind + 1))

        # ion momentum(imome[0])
        esim1 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 2 + iind]
        esim1.momentum.resize(1)
        esim1.grid_subset_indes = 1
        esim1.momentum[0].radial.resize(dsize)
        esim1.momentum[0].radial_coefficients.resize([dsize, 2])
        esim1.label = impname[0][2] + str("{:>02d}".format(iind + 1))

    # ion energy(twci1)
    esie1 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 3]
    esie1.energy.resize(1)
    esie1.grid_subset_index = 1
    esie1.energy[0].values.resize(dsize)
    esie1.energy[0].coefficients.resize([dsize, 2])
    esie1.label = impname[0][3]

    for iind in range(nzmx[1]):
        # ion particles(tionZ1)
        esip21 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 3 + 1 + nzmx[1] * 0 + iind]
        esip21.particles.resize(1)
        esip21.grid_subset_index = 1
        esip21.particles[0].values.resize(dsize)
        esip21.particles[0].coefficients.resize([dsize, 2])
        esip21.label = impname[1][0] + str("{:>02d}".format(iind + 1))

        # ion particles(trecZ1)
        esip22 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 3 + 1 + nzmx[1] * 1 + iind]
        esip22.particles.resize(1)
        esip22.grid_subset_index = 1
        esip22.particles[0].values.resize(dsize)
        esip22.particles[0].coefficients.resize([dsize, 2])
        esip22.label = impname[1][1] + str("{:>02d}".format(iind + 1))

        # ion momentum
        esim2 = imas_obj.edge_sources.source[0].ggd[0].ion[nzmx[0] * 3 + 1 + nzmx[1] * 2 + iind]
        esim2.momentum.resize(1)
        esim2.grid_subset_indes = 1
        esim2.momentum[0].radial.resize(dsize)
        esim2.momentum[0].radial_coefficients.resize([dsize, 2])
        esim2.label = impname[1][2] + str("{:>02d}".format(iind + 1))

    # ion energy(twci2)
    esie2 = imas_obj.edge_sources.source[0].ggd[0].ion[(nzmx[0] + nzmx[1]) * 3 + 1]
    esie2.energy.resize(1)
    esie2.grid_subset_index = 1
    esie2.energy[0].values.resize(dsize)
    esie2.energy[0].coefficients.resize([dsize, 2])
    esie2.label = impname[1][3]

    # data insert
    for i in range(dsize):
        # neutrals
        esnp0.particles[0].values[i] = wssn[i]
        esnm0.momentum[0].radial[i] = wssp[i]
        esne0.energy[0].values[i] = wswi[i]
        esnp011.particles[0].values[i] = tionZ01[i]
        esnp012.particles[0].values[i] = trecZ01[i]
        esnp021.particles[0].values[i] = tionZ02[i]
        esnp022.particles[0].values[i] = trecZ02[i]
        esnm1.momentum[0].radial[i] = imome01[i]
        esnm2.momentum[0].radial[i] = imome02[i]

        # ions
        ipath = imas_obj.edge_sources.source[0].ggd[0]
        for j in range(nzmx[0]):
            ipath.ion[nzmx[0] * 0 + j].particles[0].values[i] = tionZ1[j][i]
            ipath.ion[nzmx[0] * 1 + j].particles[0].values[i] = trecZ1[j][i]
            ipath.ion[nzmx[0] * 2 + j].momentum[0].radial[i]  = imome1[j][i]
        ipath.ion[nzmx[0] * 3].energy[0].values[i] = twci1[i]
        for j in range(nzmx[1]):
            ipath.ion[nzmx[0] * 3 + 1 + nzmx[1] * 0 + j].particles[0].values[i] = tionZ2[j][i]
            ipath.ion[nzmx[0] * 3 + 1 + nzmx[1] * 1 + j].particles[0].values[i] = trecZ2[j][i]
            ipath.ion[nzmx[0] * 3 + 1 + nzmx[1] * 2 + j].momentum[0].radial[i]  = imome2[j][i]
        ipath.ion[(nzmx[0] + nzmx[1]) * 3 + 1].energy[0].values[i] = twci2[i]

    return imas_obj