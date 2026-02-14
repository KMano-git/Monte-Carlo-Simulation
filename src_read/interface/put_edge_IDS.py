#! /usr/bin/env python
#Python 3.5

# -----------------------------------------------------------------------------
# Author: Dejan Penko
#
# DESCRIPTION
# This Python script is used to read geometry from b2fgmtry file together with
# electron density, electron temperature and ion temperature scalars from
# b2fstati file. The same data is then written to IDS together by
# creating "Cells" and "Nodes" grid subsets.
#
# Basic environment settings (terminal commands on hpc iter.org)
# $ module load imas/3.7.4/ual/3.4.0
# $ imasdb solps-iter
# -----------------------------------------------------------------------------

# impor imas mod
try:
    import imas
except Exception as e:
    print("Required IMAS support library not available on this system.")

# import 
import getopt
import sys
import os # add 2021/05/17

# dumy data 
from dmyIDS import dmyIDS

# put IDS
from put_source_IDS import es2IDS
from put_radiation_IDS import rd2IDS

# --------------------------------------------------
# read mesh data for b2fgmtry
# 1st argument: path, string
# --------------------------------------------------
def readB2fgmtry(file_path):
    # Read geometry from file b2fgmtry (x and y coordinates of nodes)
    # and insert them into array for later use

    crx = []
    cry = []

    found_nxyx = 0
    found_crx = 0
    found_cry = 0
    with open(file_path + "/b2fgmtry", 'r') as infile:
        for line in infile:
            lineSplit = line.split() # line split in form
                          # ['*cf:', 'int', '2', 'nx,ny']
            if len(lineSplit) >= 4:
                if lineSplit[3] == "nx,ny":
                    found_nxyx = 1
                    continue
            if (found_nxyx == 1):
                nx = int(lineSplit[0])
                ny = int(lineSplit[1])
                found_nxyx = 0

            if (found_crx == 1 and found_cry == 0):
            # Write coordinates between 'crx' and next 'cf*:' to crx array
                for j in range(len(lineSplit)):
                    if lineSplit[j] == "*cf:":
                        found_crx = 0
                        found_cry = 0
                        break
                    else:
                        crx.append(float(lineSplit[j].split('E')[0]) * \
                                   pow(10, int(lineSplit[j].split('E')[1])))
            if (found_crx == 0 and found_cry == 1):
            # Write coorindates between 'cry' and next 'cf*:' to cry array
                for j in range(len(lineSplit)):
                    if lineSplit[j] == "*cf:":
                        found_crx = 0
                        found_cry = 0
                        break
                    else:
                        cry.append(float(lineSplit[j].split('E')[0]) * \
                                   pow(10, int(lineSplit[j].split('E')[1])))
            for i in range(len(lineSplit)):
                if lineSplit[i] == "crx":
                    found_crx = 1
                    found_cry = 0
                    break
                if lineSplit[i] == "cry":
                    found_crx = 0
                    found_cry = 1
                    break
    print("nx: \%d | ny: \%d" % (nx,ny))

    return crx, cry, nx, ny

# --------------------------------------------------
# read plasma data for b2fstate
# 1st argument: path, string
# --------------------------------------------------
def readB2fstate(file_path):
    # Read values from file b2fstate
    # (electron density (ne), electron temperature(te), ion temperature(ti))
    # and insert them into array for later use
    [vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea] = [[], [], [], [], [], [], [], [], [], [], []]

    # find flg
    # 0x001:vna, 0x002:vne, 0x004:vni, 0x008:vnezef
    # 0x010:vzf, 0x020:vva, 0x040:vve, 0x080:vti
    # 0x100:vte, 0x200:vcs, 0x400:vea
    fdf = 0x000
    
    with open(file_path + "/b2fstate", 'r') as infile:
        for line in infile:
            lineSplit = line.split() # line split in form
                          # ['*cf:', 'int', '2', 'nx,ny']
            
            for j in range(len(lineSplit)):
                if lineSplit[j] == "*cf:":
                    fdf = 0x000
                    break
                elif fdf != 0x000:
                    ld = float(lineSplit[j].split('E')[0]) * pow(10, int(lineSplit[j].split('E')[1]))
                    if fdf == 0x001:
                        vna.append(ld)
                    if fdf == 0x002:
                        vne.append(ld)
                    if fdf == 0x004:
                        vni.append(ld)
                    if fdf == 0x008:
                        vnezef.append(ld)
                    if fdf == 0x010:
                        vzf.append(ld)
                    if fdf == 0x020:
                        vva.append(ld)
                    if fdf == 0x040:
                        vve.append(ld)
                    if fdf == 0x080:
                        vti.append(ld)
                    if fdf == 0x100:
                        vte.append(ld)
                    if fdf == 0x200:
                        vcs.append(ld)
                    if fdf == 0x400:
                        vea.append(ld)
            
            for i in range(len(lineSplit)):
                if lineSplit[i] == "vna":
                    fdf = 0x001
                    break
                if lineSplit[i] == "vne":
                    fdf = 0x002
                    break
                if lineSplit[i] == "vni":
                    fdf = 0x004
                    break
                if lineSplit[i] == "vnezef":
                    fdf = 0x008
                    break
                if lineSplit[i] == "vzf":
                    fdf = 0x010
                    break
                if lineSplit[i] == "vva":
                    fdf = 0x020
                    break
                if lineSplit[i] == "vve":
                    fdf = 0x040
                    break
                if lineSplit[i] == "vti":
                    fdf = 0x080
                    break
                if lineSplit[i] == "vte":
                    fdf = 0x100
                    break
                if lineSplit[i] == "vcs":
                    fdf = 0x200
                    break
                if lineSplit[i] == "vea":
                    fdf = 0x400
                    break
    pls = [vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea]     
    return pls

# --------------------------------------------------
# read neutral data
# 1st argument: path, string
# --------------------------------------------------
def readNeutral(file_path):
    [wssn, wssp, wswe, wswi, wden, tdnz, tionZ, trecZ, tfrz, tthz] = [[], [], [], [], [], [], [], [], [], []]
    fdf = 0x000 # 0x001:wssn,0x002:wssp,0x004:wswe,0x008:wswi,0x010:wden
                # 0x020:tdnz,0x040:tionZ,0x080:trecZ,0x100:tfrz,0x200:tthz

    with open(file_path + "/neutral", "r") as nt:
        for line in nt:
            lsplit = line.split()
            
            for j in range(len(lsplit)):
                if lsplit[j] == "*cf:": 
                    fdf = 0x000
                    break
                elif fdf != 0x000:
                    ld = float(lsplit[j].split('E')[0]) * pow(10, int(lsplit[j].split('E')[1]))
                    if fdf == 0x001:
                        wssn.append(ld)
                    if fdf == 0x002:
                        wssp.append(ld)
                    if fdf == 0x004:
                        wswe.append(ld)
                    if fdf == 0x008:
                        wswi.append(ld)
                    if fdf == 0x010:
                        wden.append(ld)
                    if fdf == 0x020:
                        tdnz[len(tdnz)-1].append(ld)
                    if fdf == 0x040:
                        tionZ[len(tionZ)-1].append(ld)
                    if fdf == 0x080:
                        trecZ[len(trecZ)-1].append(ld)
                    if fdf == 0x100:
                        tfrz[len(tfrz)-1].append(ld)
                    if fdf == 0x200:
                        tthz[len(tthz)-1].append(ld)

            for i in range(len(lsplit)):
                if lsplit[i] == "wssn":
                    fdf = 0x001
                    break
                if lsplit[i] == "wssp":
                    fdf = 0x002
                    break
                if lsplit[i] == "wswe":
                    fdf = 0x004
                    break
                if lsplit[i] == "wswi":
                    fdf = 0x008
                    break
                if lsplit[i] == "wden":
                    fdf = 0x010
                    break
                if lsplit[i].find("tdnz")>=0:
                    fdf = 0x020
                    tdnz.append([])
                    break
                if lsplit[i].find("tionZ")>=0:
                    fdf = 0x040
                    tionZ.append([])
                    break
                if lsplit[i].find("trecZ")>=0:
                    fdf = 0x080
                    trecZ.append([])
                    break
                if lsplit[i].find("tfrz")>=0:
                    fdf = 0x100
                    tfrz.append([])
                    break
                if lsplit[i].find("tthz")>=0:
                    fdf = 0x200
                    tthz.append([])
                    break

    vm  = [wssn, wssp, wswe, wswi, wden]
    vm += [tdnz, tionZ, trecZ, tfrz, tthz]
    return vm

# --------------------------------------------------
# read imprity data
# 1st argument: path, string
# --------------------------------------------------
def readImprity(file_path):
    [twrd, tdnz, tfrz, tthz, tionZ, trecZ, twci] = [[], [], [], [], [], [], []]
    fdf = 0x00 # 0x01:twrd,0x02:tdnz,0x04:tionZ,0x08:trecZ,0x10:twci,0x20:tfrz,0x40:tthz
    with open(file_path + "/impurity1", "r") as im:
        for line in im:
            lsplit = line.split()

            for j in range(len(lsplit)):
                if lsplit[j] == "*cf:": 
                    fdf = 0x00
                    break
                elif fdf != 0x00:
                    ld = float(lsplit[j].split('E')[0])*pow(10,int(lsplit[j].split('E')[1]))
                    if fdf == 0x01:
                        twrd.append(ld)
                    if fdf == 0x02:
                        tdnz[len(tdnz)-1].append(ld)
                    if fdf == 0x04:
                        tionZ[len(tionZ)-1].append(ld)
                    if fdf == 0x08:
                        trecZ[len(trecZ)-1].append(ld)
                    if fdf == 0x10:
                        twci.append(ld)
                    if fdf == 0x20:
                        tfrz[len(tfrz)-1].append(ld)
                    if fdf == 0x40:
                        tthz[len(tthz)-1].append(ld)

            for i in range(len(lsplit)):
                if lsplit[i].find("twrd")>=0:
                    fdf = 0x01
                    break
                if lsplit[i].find("tdnz")>=0:
                    fdf = 0x02
                    tdnz.append([])
                    break
                if lsplit[i].find("tionZ")>=0:
                    fdf = 0x04
                    tionZ.append([])
                    break
                if lsplit[i].find("trecZ")>=0:
                    fdf = 0x08
                    trecZ.append([])
                    break
                if lsplit[i].find("twci")>=0:
                    fdf = 0x10
                    break
                if lsplit[i].find("tfrz")>=0:
                    fdf = 0x20
                    tfrz.append([])
                    break
                if lsplit[i].find("tthz")>=0:
                    fdf = 0x40
                    tthz.append([])
                    break

    im1 = [twrd, tdnz, tfrz, tthz, tionZ, trecZ, twci]
    [twrd, tdnz, tfrz, tthz, tionZ, trecZ, twci] = [[], [], [], [], [], [], []]
    fdf = 0x00 # 0x01:twrd,0x02:tdnz,0x04:tionZ,0x08:trecZ,0x10:twci,0x20:tfrz,0x40:tthz
    with open(file_path + "/impurity2", "r") as im:
        for line in im:
            lsplit = line.split()

            for j in range(len(lsplit)):
                if lsplit[j] == "*cf:": 
                    fdf = 0x00
                    break
                elif fdf != 0x00:
                    ld = float(lsplit[j].split('E')[0])*pow(10,int(lsplit[j].split('E')[1]))
                    if fdf == 0x01:
                        twrd.append(ld)
                    if fdf == 0x02:
                        tdnz[len(tdnz)-1].append(ld)
                    if fdf == 0x04:
                        tionZ[len(tionZ)-1].append(ld)
                    if fdf == 0x08:
                        trecZ[len(trecZ)-1].append(ld)
                    if fdf == 0x10:
                        twci.append(ld)
                    if fdf == 0x20:
                        tfrz[len(tfrz)-1].append(ld)
                    if fdf == 0x40:
                        tthz[len(tthz)-1].append(ld)

            for i in range(len(lsplit)):
                if lsplit[i].find("twrd")>=0:
                    fdf = 0x01
                    break
                if lsplit[i].find("tdnz")>=0:
                    fdf = 0x02
                    tdnz.append([])
                    break
                if lsplit[i].find("tionZ")>=0:
                    fdf = 0x04
                    tionZ.append([])
                    break
                if lsplit[i].find("trecZ")>=0:
                    fdf = 0x08
                    trecZ.append([])
                    break
                if lsplit[i].find("twci")>=0:
                    fdf = 0x10
                    break
                if lsplit[i].find("tfrz")>=0:
                    fdf = 0x20
                    tfrz.append([])
                    break
                if lsplit[i].find("tthz")>=0:
                    fdf = 0x40
                    tthz.append([])
                    break

    im2 = [twrd, tdnz, tfrz, tthz, tionZ, trecZ, twci]

    IM  = [im1[0], im1[1], im1[2], im1[3], im1[4], im1[5], im1[6]]
    IM += [im2[0], im2[1], im2[2], im2[3], im2[4], im2[5], im2[6]]
    
    return IM

# --------------------------------------------------
# write data to IMAS
# 1st argument:input(shot, run, user, device, version), list(int, int, str, str, str)
# 2nd argument:geometry data, list(float 1D, float 1D, int, int)
# 3rd argument:plasma data(ne, te, ti), list(float 1D, float 1D, float 1D)
# 4th argument:neutral data(wssn, wssp, wswe, wswi, wden, tdnz[0], tionZ[0], trecZ[0]), list(float 1D's)
# 5th argument:implity data*2(twrd, tdnz, tionZ, trecZ), list(float 1D's)
# 6th argument:elements(element name, element number), list(str, int)
# --------------------------------------------------
def B2toIDS(inp, ncxy, pls, ntl, imp, elm): 
    # set data
    [shot, run, user, device, version] = inp
    [xc, yc, nx, ny] = ncxy
    [vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea] = pls
    [wssn, wssp, wswe, wswi, wden] = ntl[:5]
    [tdnz0, tionZ0, trecZ0, tfrz0, tthz0] = ntl[5:]
    [tdnz01, tdnz02] = tdnz0[:]
    [tionZ01, tionZ02] = tionZ0[:]
    [trecZ01, trecZ02] = trecZ0[:]
    [tfrz01, tfrz02] = tfrz0[:]
    [tthz01, tthz02] = tthz0[:]
    [twrd1, tdnz1, tfrz1, tthz1, tionZ1, trecZ1, twci1, twrd2, tdnz2, tfrz2, tthz2, tionZ2, trecZ2, twci2] = imp
    
    # Write previously found data in b2fgmtry and b2fstati to IDS database
    print('Writing IDS: ')
    time = 1
    interp = 1

    # --Create IDS database--
    imas_obj = imas.ids(shot, run, shot, run)

    # imas_obj.create() # Create the data entry
    imas_obj.create_env(user, device, version)

    if imas_obj.isConnected():
        print('Creation of data entry OK!')
    else:
        print('Creation of data entry FAILED!')
        sys.exit()

    # --Basic IDS space allocation--
    imas_obj.edge_profiles.profiles_1d.resize(1)
    imas_obj.edge_profiles.ggd.resize(1)
    imas_obj.edge_profiles.grid_ggd.resize(1)
    #imas_obj.edge_profiles.putNonTimed()
    imas_obj.edge_profiles.time.resize(1)
    imas_obj.edge_profiles.time[0] = 1
    imas_obj.edge_profiles.ids_properties.homogeneous_time = 1

    # dumy
    dirpath = os.getcwd()
    
    # Set IDS grid description
    grid_description = "This is IDS" + \
      " shot=" + str(shot) + " run=" + str(run) + " user=" + str(user) + \
      " device=" + str(device) + " version=" + str(version) + \
      " written by put_edge_ids using b2fgmtry and b2fstati files found " +\
      "in directory" + dirpath + "."
    # Put IDS grid description
    #imas_obj.edge_profiles.ggd[0].grid.identifier.description = grid_description
    #imas_obj.edge_profiles.ggd[0].grid.space.resize(1)
    
    imas_obj.edge_profiles.grid_ggd[0].identifier.description = grid_description
    imas_obj.edge_profiles.grid_ggd[0].space.resize(1)
    
    # Set (IDS substructure shortcut variable) space0
    #space0 = imas_obj.edge_profiles.ggd[0].grid.space[0]
    space0 = imas_obj.edge_profiles.grid_ggd[0].space[0]
    space0.objects_per_dimension.resize(3)

    num_obj_0D = len(xc) # Number of nodes (len(xc) == len(yc))
    # # have 2D coordinates, P(x,y)
    num_coord = len(xc) + len(yc) # Number of all available coordinates
    num_gridSubsets = 2 # Number of grid subsets to write (Cells and Nodes)

    space0.coordinates_type.resize(2)
    space0.coordinates_type[0] = 1 # X
    space0.coordinates_type[1] = 2 # Y

    #imas_obj.edge_profiles.ggd[0].grid.grid_subset.resize(num_gridSubsets)
    
    imas_obj.edge_profiles.grid_ggd[0].grid_subset.resize(num_gridSubsets)
   
    # -- Put DATA FOR GRID SUBSET "Nodes" --
    # (grid subset index: 2, objects forming the grid subset: nodes, 0D)
    # Note: All indices must be put in Fortran index notation
    # (starting with 1), not in C++/python index notation(starts with 0)!
    # So in our case: Python_Index == Fortran_Index -1 !

    gridSubset_index = 2 # grid subset index of grid subset Nodes
                # (Indexing of grid subsets follows the IDS
                # examples :
                # shot: 1, run:1, # device: iter; and
                # shot: 16151, run: 1000; device: aug
    gridSubset_name = "Nodes"
    gridSubset_dim_index = 1 # Grid subset Nodes consists of
                  # points -> 0D objects -> dimension index = 1
                  # (edges -> 1D objects -> dimension index = 2
                  # cells -> 2D objects -> dimension index = 3)


    # Write all available 0D objects (all of them form the grid subset Nodes
    # Set (IDS substructure shortcut variable) dim0
    dim0 = space0.objects_per_dimension[gridSubset_dim_index - 1]
    dim0.object.resize(num_obj_0D)
    for i in range(num_obj_0D):
        dim0.object[i].nodes.resize(1)
        dim0.object[i].nodes[0] = i
        dim0.object[i].geometry.resize(2)
        dim0.object[i].geometry[0] = xc[i]
        dim0.object[i].geometry[1] = yc[i]

    # Set(IDS substructure shortcut variable) gridSubsetBaseData
    gridSubsetBaseData = \
        imas_obj.edge_profiles.grid_ggd[0].grid_subset[gridSubset_index-1]
    # Put base grid subset data/parameters (name, index)
    gridSubsetBaseData.identifier.name = gridSubset_name
    gridSubsetBaseData.identifier.index = gridSubset_index
    # Put grid subset element and element object data
    gridSubsetBaseData.element.resize(num_obj_0D)
    for i in range(num_obj_0D):
        gridSubsetBaseData.element[i].object.resize(1)
        gridSubsetBaseData.element[i].object[0].space = 0 + 1
        gridSubsetBaseData.element[i].object[0].dimension = gridSubset_dim_index
        gridSubsetBaseData.element[i].object[0].index = i + 1

    # -- Put DATA FOR GRUD SUBSET "Cells"
    # (grid subset index: 1, objects forming the grid subset: cells, 2D)

    # !!!
    # numCellsX = nx + 2
    # numCellsY = ny + 2
    numCellsX = nx
    numCellsY = ny
    num_cells = numCellsX * numCellsY
    num_obj_2D = num_cells
    gridSubset_index = 1
    gridSubset_name = "Cells"
    gridSubset_dim_index = 3

    # Write all available 2D objects (all of them form the grid subset Cells
    dim2 = space0.objects_per_dimension[gridSubset_dim_index - 1]
    dim2.object.resize(num_obj_2D)
    for i in range(num_obj_2D):
        dim2.object[i].nodes.resize(4)
    cellId = 1
    for j in range(numCellsY):
        for i in range(numCellsX):
            dim2.object[cellId - 1].nodes[0] = cellId+0*numCellsX*numCellsY
            dim2.object[cellId - 1].nodes[1] = cellId+1*numCellsX*numCellsY
            dim2.object[cellId - 1].nodes[2] = cellId+3*numCellsX*numCellsY
            dim2.object[cellId - 1].nodes[3] = cellId+2*numCellsX*numCellsY
            cellId += 1

    # Set (IDS substructure shortcut variable) subgridDaseData for
    # Cells grid subset
    gridSubsetBaseData = \
        imas_obj.edge_profiles.grid_ggd[0].grid_subset[gridSubset_index - 1]
    # Put base grid subset data/parameters (name, index)
    gridSubsetBaseData.identifier.name = gridSubset_name
    gridSubsetBaseData.identifier.index = gridSubset_index
    # Put grid subset element and element object data
    gridSubsetBaseData.element.resize(num_obj_2D)
    for i in range(num_obj_2D):
        gridSubsetBaseData.element[i].object.resize(1)
        gridSubsetBaseData.element[i].object[0].space = 0 + 1
        gridSubsetBaseData.element[i].object[0].dimension = gridSubset_dim_index
        gridSubsetBaseData.element[i].object[0].index = i + 1

    # PUT VALUES (ne, te, ti) for "Cells" grid subset
    # Put ne (electron density)
    num_ne_gridSubset = 1
    num_ne_values = len(vne)
    imas_obj.edge_profiles.ggd[0].electrons.density.resize(num_ne_gridSubset)
    nePath = \
        imas_obj.edge_profiles.ggd[0].electrons.density[num_ne_gridSubset - 1]
    nePath.grid_subset_index = gridSubset_index
    nePath.values.resize(num_ne_values)
    for n in range(num_ne_values):
        nePath.values[n] = vne[n] 
    
    # Put te (electron temperature)
    num_te_gridSubset = 1
    num_te_values = len(vte)
    imas_obj.edge_profiles.ggd[0].electrons.temperature.resize(num_te_gridSubset)
    tePath = \
        imas_obj.edge_profiles.ggd[0].electrons.temperature[num_te_gridSubset - 1]
    tePath.grid_subset_index = gridSubset_index
    tePath.values.resize(num_te_values)
    for n in range(num_te_values):
    # convert to eV (1 J = 6.242e18 eV)
        tePath.values[n] = vte[n] * (6.242e18)
    
    # Put ti (ion temperature)
    num_ti_gridSubset = 1
    num_ti_values = len(vti)
    num_ti_species = 1 + elm[0][1] + elm[1][1] # Number of ion species, as in number of different ion charges.
    ion_specie = 1
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].ion.resize(num_ti_species)
    imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].temperature.resize(num_ti_gridSubset)
    imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].label = "ION"
    tiPath = imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].temperature[num_ti_gridSubset - 1]
    tiPath.grid_subset_index = gridSubset_index
    tiPath.values.resize(num_ti_values)
    for n in range(num_ti_values):
    # convert to eV (1 J = 6.242e18 eV)
        tiPath.values[n] = vti[n] * (6.242e18)
    
    # Put ve
    num_ve_gridSubset = 1
    num_ve_values = len(vve)
    imas_obj.edge_profiles.ggd[0].electrons.velocity.resize(num_ve_gridSubset)
    vePath = imas_obj.edge_profiles.ggd[0].electrons.velocity[num_ve_gridSubset - 1]
    vePath.grid_subset_index = gridSubset_index
    vePath.parallel.resize(num_ve_values)
    for n in range(num_ve_values):
        vePath.parallel[n] = vve[n]

    # Put va
    num_va_gridSubset = 1
    num_va_values = len(vva)
    ion_specie = 1
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].velocity.resize(num_va_gridSubset)
    vaPath = imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].velocity[num_va_gridSubset - 1]
    vaPath.grid_subset_index = gridSubset_index
    vaPath.parallel.resize(num_va_values)
    for n in range(num_va_values):
        vaPath.parallel[n] = vva[n]

    # Put ni
    num_ni_gridSubset = 1
    num_ni_values = len(vni)
    ion_specie = 1
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].density.resize(num_ni_gridSubset)
    niPath = imas_obj.edge_profiles.ggd[0].ion[ion_specie - 1].density[num_ni_gridSubset - 1]
    niPath.grid_subset_index = gridSubset_index
    niPath.values.resize(num_ni_values)
    for n in range(num_ni_values):
        niPath.values[n] = vni[n]
   
    # PUT VALUES (tdnz[1~]) for "Cells" grid subset
    # Put tdnz[1~]
    num_tdnz1_gridSubset = elm[0][1]
    ion_specie = 2

    tdnz = []
    for i in range(num_tdnz1_gridSubset):
        tdnz.append(tdnz1[i])

    for i in range(ion_specie - 1, num_tdnz1_gridSubset + 1):
        imas_obj.edge_profiles.ggd[0].ion[i].density.resize(1)

        num_tdnz_values = len(tdnz[i - 1])

        imas_obj.edge_profiles.ggd[0].ion[i].label = "TDNZ_" + elm[0][0] + "_" + str("{:>02d}".format(i))
        tdnzPath = imas_obj.edge_profiles.ggd[0].ion[i].density[0]
        tdnzPath.grid_subset_index = gridSubset_index
        tdnzPath.values.resize(num_tdnz_values)
        for n in range(num_tdnz_values):
            tdnzPath.values[n] = tdnz[i - 1][n]

    # Put tdnz[1~]
    num_tdnz2_gridSubset = elm[1][1]
    ion_specie = 1 + num_tdnz1_gridSubset
    tdnz = []
    for i in range(num_tdnz2_gridSubset):
        tdnz.append(tdnz2[i])
        
    for i in range(ion_specie, ion_specie + num_tdnz2_gridSubset):
        imas_obj.edge_profiles.ggd[0].ion[i].density.resize(1)

        num_tdnz_values = len(tdnz[i - ion_specie])

        imas_obj.edge_profiles.ggd[0].ion[i].label = "TDNZ_" + elm[1][0] + "_" + str("{:>02d}".format(i - ion_specie + 1))
        tdnz2Path = imas_obj.edge_profiles.ggd[0].ion[i].density[0]
        tdnz2Path.grid_subset_index = gridSubset_index
        tdnz2Path.values.resize(num_tdnz_values)
        for n in range(num_tdnz_values):
            tdnz2Path.values[n] = tdnz[i - ion_specie][n]

   # PUT VALUES (wden) for "Cells" grid subset
    # Put wden
    num_wden_gridSubset = 1
    num_wden_values = len(wden)
    num_wden_species = 4 # Number of ion species, as in
              # number of different ion charges.
    nutlal_specie = 2
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].neutral.resize(num_wden_species)
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density.resize(num_wden_gridSubset)
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].label = "WDEN"
    wdenPath = imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density[num_wden_gridSubset - 1]
    wdenPath.grid_subset_index = gridSubset_index
    wdenPath.values.resize(num_wden_values)
    for n in range(num_wden_values):
        wdenPath.values[n] = wden[n]

    # PUT VALUES (tdnz[0]) for "Cells" grid subset
    # Put tdnz[0](Ar)
    num_tdnz_gridSubset = 1
    num_tden_values = len(tdnz01)
    num_wden_species = 1 # Number of ion species, as in
              # number of different ion charges.
    nutlal_specie = 3
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density.resize(num_tdnz_gridSubset)
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].label = "TDNZ_" + elm[0][0]
    tdnzPath = imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density[num_tdnz_gridSubset - 1]
    tdnzPath.grid_subset_index = gridSubset_index
    tdnzPath.values.resize(num_tden_values)
    for n in range(num_tden_values):
        tdnzPath.values[n] = tdnz01[n]

    # PUT VALUES (tdnz[0]) for "Cells" grid subset
    # Put tdnz[0](Ar)
    num_tdnz_gridSubset = 1
    num_tden_values = len(tdnz02)
    num_wden_species = 1 # Number of ion species, as in
              # number of different ion charges.
    nutlal_specie = 4
    # Ion specie is linked with the ion density of each ion charge,
    # as ion temperature is taken as the same for all ion charges.
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density.resize(num_tdnz_gridSubset)
    imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].label = "TDNZ_" + elm[1][0]
    tdnzPath = imas_obj.edge_profiles.ggd[0].neutral[nutlal_specie - 1].density[num_tdnz_gridSubset - 1]
    tdnzPath.grid_subset_index = gridSubset_index
    tdnzPath.values.resize(num_tden_values)
    for n in range(num_tden_values):
        tdnzPath.values[n] = tdnz02[n]

    # edge_source
    data = [wssn, wssp, wswi, twci1, twci2, twrd1, twrd2, tionZ01, tionZ1, tionZ02, tionZ2, \
           trecZ01, trecZ1, trecZ02, trecZ2, tfrz01, tfrz1, tfrz02, tfrz2, tthz01, tthz1, tthz02, tthz2, \
           tdnz01, tdnz1, tdnz02, tdnz2]
    imas_obj = es2IDS(imas_obj, data, [nx, ny], elm)

    # radiation
    imas_obj = rd2IDS(imas_obj, [wswe, twrd1, twrd2], [nx, ny], elm)

    # dumy data
    imas_obj = dmyIDS(imas_obj, [nx, ny], [xc, yc])

    # Write all put data do IDS
    imas_obj.edge_profiles.putSlice()
    imas_obj.edge_sources.putSlice()
    imas_obj.edge_transport.putSlice()
    imas_obj.radiation.putSlice()

    # Close IDS
    imas_obj.close()
    print("IDS write finished")
    print("IDS closed")

if __name__ == "__main__":
    
    # For launching python script directly from treminal with python command
    try:
        opts, args = getopt.getopt(sys.argv[1:], "srutvh", ["dirpath=",
                                  "shot=", "run=",
                                  "user=", "device=",
                                  "version=", "help"])

        for opt, arg in opts:
            #print opt, arg
            if opt in ("-fp", "--dirpath"):
                dirpath = arg
            elif opt in ("-s", "--shot"):
                shot = int(arg)
            elif opt in ("-r", "--run"):
                run = int(arg)
            elif opt in ("-u", "--user"):
                user = arg
            elif opt in ("-t", "--device"):
                device = arg
            elif opt in ("-v", "--version"):
                version = arg

            if opt in ("-h", "--help"):
                print("In order to run b2read file path, shot, run, user,"
                  "device and version variables must be defined."
                  "Example (terminal): "
                  "python3.5 put_edge_ids.py "
                  "--dirpath=/home/ITER/penkod/solps-iter/runs/AUG_16151_D/"
                  "baserun "
                  "--shot=1000 --run=1 --user=penkod --device=solps-iter "
                  "--version=3")
                sys.exit()

        dirpath, shot, run, user, device, version
    except getopt.GetoptError:
        print ('Supplied option not recognized!')
        print ('For help: b2read -h / --help')
        sys.exit(2)

    # few paths to example files for testing
    # /home/ITER/tomsicp/solps-iter/runs/AUG_16151_D/baserun
    # /home/ITER/tomsicp/solps-iter-devel/runs/ITER_535_D+He+Ar/baserun

    cnxy = readB2fgmtry(dirpath)
    pls = readB2fstate(dirpath)
    ntl = readNeutral(dirpath)
    imp = readImprity(dirpath)
    inp = [shot, run, user, device, version]
    elem = [["Ar", 18], ["C", 6]]

    B2toIDS(inp, cnxy, pls, ntl, imp, elem)
