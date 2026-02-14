# ======================================================================
# -*- name  : sonic2imas.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import os

from getData import getYAML, getEnv, DTPLV_MSH

# new import 2021
from getName import getName, getElem
from checkInp import checkInp
from getData2 import DTIMP_MSH, DTNTL_MSH, IMPSY_MSH
# from getData2 import DTNTL_MSH, IMPSY_MSH
from readMSH2 import readMSH2

# import IDS editing library
from put_edge_IDS import readB2fgmtry, readB2fstate, B2toIDS, readImprity, readNeutral

# ======================================================================
# main program
# store sonic outputs in IMAS
# ======================================================================
def sonic2imas():

    try:
        config_data = getYAML()
        env = config_data
        config_paths = config_data[0]['Paths']
        #Much better to explicitly point data file name in the config file.
        fileList[0] = config_paths['bin_pls'] #plasma prof
        fileList[1] = config_paths['bin_mtrc'] # mesh data
        fileList[2] = config_paths['asc_wall'] # wall data
        fileList[3] = config_paths['bin_ntl'] # neutrals prof
        fileList[4] = config_paths['bin_imp1'] # impurity prof 1st spc
        fileList[5] = config_paths['bin_imp2'] # impurity prof 2nd spc
        fileList[6] = config_paths['bin_ims1'] # additional impurity prof 1st spc
        fileList[7] = config_paths['bin_ims2'] # additional impurity prof 2nd spc

    except:
        # getting data
        env = getEnv()

        # setting paths
        # DTPLV, dat.mtrc, wall_gui.dat, DTNTL, DTIMP1, DTIMP2, IMPSY1, IMPSY2, dat.mont
        fileList = env[0]

    # file path is NULL
    for i in range(len(fileList)-1):
        if fileList[i]=='':
            print('Paths undefined in path.txt.')
            print('sys exit.')
            sys.exit()
    
    print()
    
    # get pls
    [X, Y, Pls, ind] = DTPLV_MSH(fileList[0], fileList[1])

    plName = ["vna", "vne", "vni", "vnezef", "vzf", "vva", "vve", "vti", "vte", "vcs", "vea"]

    pl = []
    for i in range(len(Pls)):
        [p, lx, ly] = create_data(Pls[i])
        if plName[i].find("vt") >= 0:
            pl.append(eV2J(p))
        else:
            pl.append(p)
    print("plasma data created.\n")
    
    # get mont mesh data
    [cntgrd, cntwal, cntgat, com_size] = readMSH2(fileList[8])

    ncmax = cntgrd[32] # 5806
    iaxs = cntgrd[44] # 62
    ind.append(iaxs)

    # get impurity data
    [Imp1, nzmx1] = DTIMP_MSH(fileList[4], ind, ncmax) # twci
    [Imp2, nzmx2] = DTIMP_MSH(fileList[5], ind, ncmax) # twci
    [Ims1, nzmx1] = IMPSY_MSH(fileList[6], ind, ncmax) # twrd, tdnz, tfrz, tthz, tionZ, trecZ
    [Ims2, nzmx2] = IMPSY_MSH(fileList[7], ind, ncmax) # twrd, tdnz, tfrz, tthz, tionZ, trecZ
    im = []
    Ims1.append(Imp1[len(Imp1)-1])
    Ims2.append(Imp2[len(Imp2)-1])
    IM = [Ims1, Ims2]
    for i in range(len(IM)):
        im.append([])
        for j in range(len(IM[i])):
            im[i].append(create_data(IM[i][j])[0])
            # im = [[twrd1, tdnz1, tfrz1, tthz1, tionZ1, trecZ1, twci1]
            #       [twrd2, tdnz2, tfrz2, tthz2, tionZ2, trecZ2, twci2]]
        print("impurity data", i+1, "created.\n")
    
    # create element list
    elem =[]
    elem.append(getElem(nzmx1))
    elem.append(getElem(nzmx2))

    # create imprity name list # 20220303 fix
    [imsName1,imsName2] = [[],[]]
    nzmx = [nzmx1, nzmx2]
    for i in range(len(nzmx)):
        inl = []
        for j in range(1+(nzmx[i]+1)*5+1):
            inl.append("")
        inl[0] = "twrd_" + elem[i]
        for j in range(nzmx[i]+1):
            inl[1+(nzmx[i]+1)*0+j] = "tdnz_"  + elem[i] + "_" + str(j)
            inl[1+(nzmx[i]+1)*1+j] = "tfrz_"  + elem[i] + "_" + str(j)
            inl[1+(nzmx[i]+1)*2+j] = "tthz_"  + elem[i] + "_" + str(j)
            inl[1+(nzmx[i]+1)*3+j] = "tionZ_" + elem[i] + "_" + str(j)
            inl[1+(nzmx[i]+1)*4+j] = "trecZ_" + elem[i] + "_" + str(j)
        inl[len(inl)-1] = "twci_" + elem[i]
        if i == 0: imsName1 = inl
        if i == 1: imsName2 = inl

    imName = [imsName1, imsName2]

    # get nutral data
    # make ntl_NameList
    ntName  = ["wssn", "wssp", "wswe", "wswi", "wsbr", "wden"]
    for i in range(len(IM)):
        ntName.append("tdnz_"  + elem[i])
        ntName.append("tfrz_" + elem[i])
        ntName.append("tthz_" + elem[i])
        ntName.append("tionZ_" + elem[i])
        ntName.append("trecZ_" + elem[i])

    Nt = DTNTL_MSH(fileList[3], ind, ncmax)
    nt = []
    for i in range(len(Nt)):
        nt.append(create_data(Nt[i])[0])
    #nzmx = [nzmx1, nzmx2]
    for i in range(len(im)): # len(im) = 2
        for j in range(5):
            nt.append(im[i][1+nzmx[i]*j+j])
            #del im[i][1+nzmx[i]*j]
        for j in range(4, -1, -1): # 20220303 fix
            del im[i][1+nzmx[i]*j+j]
            del imName[i][1+nzmx[i]*j+j]
    print("nutral data created.\n")

    """
    # create imprity name list
    [imsName1,imsName2] = [[],[]]

    # imp1 Ar
    imsName1.append("twrd_" + elem[0])
    for i in range(1, nzmx1+1):
        imsName1.append("tdnz_" + elem[0] + "_" + str(i))
        imsName1.append("tfrz_" + elem[0] + "_" + str(i))
        imsName1.append("tthz_" + elem[0] + "_" + str(i))
        imsName1.append("tionZ_" + elem[0] + "_" + str(i))
        imsName1.append("trecZ_" + elem[0] + "_" + str(i))
    imsName1.append("twci_" + elem[0])

    # imp2 C
    imsName2.append("twrd_" + elem[1])
    for i in range(1, nzmx2+1):
        imsName2.append("tdnz_" + elem[1] + "_" + str(i))
        imsName2.append("tfrz_" + elem[1] + "_" + str(i))
        imsName2.append("tthz_" + elem[1] + "_" + str(i))
        imsName2.append("tionZ_" + elem[1] + "_" + str(i))
        imsName2.append("trecZ_" + elem[1] + "_" + str(i))
    imsName2.append("twci_" + elem[1])
    
    imName = [imsName1, imsName2]
    """
    # get x,y
    x = create_geometry(X)
    y = create_geometry(Y)
    print("geometry data created.\n")
    
    # change sonic -> solps
    S2SP(pl, plName, lx, ly)
    S2SN(nt, ntName, lx, ly)
    S2SI(im, imName, lx, ly)
    S2SG(x, y, lx, ly)
    

    # get device, version
    [ddev, dver] = getName(env[1])

    # set loop flg
    flg = 1

    while flg == 1:
        # user
        user = input("user >>>")
        if user=="":
            user = env[1] # default
            print("WARNING : The default value ( login name[", user, "] ) for \"user\" has been set.")

        # device
        device = input("device >>>")
        if device=="":
            device = ddev[0] # default
            print("WARNING : The default value (", device, ") for \"device\" has been set.")

        # version
        version = input("version >>>")
        if version=="":
            version = dver[0][0][0] # default
            print("WARNING : The default value (", version, ") for \"version\" has been set.")

        # shot
        try:
            shot = int(input("shot[0 to 214748] >>>"))
        except ValueError:
            shot = -1
        if shot<0 or 214748<shot:
            shot = 1 # default
            print("WARNING : The default value (", shot, ") for \"shot\" has been set.")

        # run
        try:
            run = int(input("run[0 to 99999] >>>"))
        except ValueError:
            run = -1
        if run<0 or 99999<run:
            run = 1 # default
            print("WARNING : The default value (", run, ") for \"run\" has been set.")

        # check input
        flg = checkInp(user, device, version, run, shot)
        # print(flg)
        if flg == 1:
            print("WARNING : shot =", shot, "and run =", run, "are used!")
            while 1:
                ans = input("\nDo you want to overwrite?[y / n] >>>")
                if ans == 'y' or ans == 'Y' or ans == 'n' or ans == 'N':
                    break
            if ans == "y" or ans == "Y":
                flg = 0
                print("overwrite with shot =", shot, "and run =", run)
            else:
                print("Please re-enter.")

    # get data from file
    dirpath = os.getcwd() # get current directory(full path)
    xc, yc, nx, ny = readB2fgmtry(dirpath)
    pls = readB2fstate(dirpath)
    imp = readImprity(dirpath)
    ntl = readNeutral(dirpath)

    # set elements
    element = [[elem[0], nzmx1], [elem[1], nzmx2]]

    # write IDS
    inp  = [shot, run, user, device, version]
    cnxy = [xc, yc, nx, ny]
    B2toIDS(inp, cnxy, pls, ntl, imp, element)

#-------------------------------------------------------------------------------
# formatting plasma data
# 1st argument:1D array(electron density), float
# 2nd argument:1D array(electron temperature), float
# 4th argumene:X length, int
# 5th argument:Y length, int
#-------------------------------------------------------------------------------
def S2SP(data, name, lx, ly):
    # setting data
    nx = lx
    ny = ly - 1
    ns = 1

    # deta format
    fmtd = "{:>22,.13E}"

    # write data
    with open("b2fstate","w") as spo:
        # length and size
        spo.write("*cf:    int                3    nx,ny,ns              \n")
        spo.write("        " + str(nx) + "          " + str(ny) + "          " + str(ns) + "\n")

        for i in range(len(data)):
            cnt = 0
            spo.write("*cf:    real               " + str(len(data[i])) + "    " + name[i] + "\n")
            for j in range(len(data[i])):
                if cnt == 6:
                    cnt = 0
                    spo.write("\n")
                spo.write(fmtd.format(data[i][j]))
                cnt += 1
            spo.write("\n")

#-------------------------------------------------------------------------------
# formatting geometry data
# 1st argument:1D array(X geometry), float
# 2nd argument:1D array(Y geometry), float
# 3rd argument:X length, int
# 4th argumene:Y length, int
#-------------------------------------------------------------------------------
def S2SG(x, y, lx, ly):
    # setting data
    nx = lx
    ny = ly - 1

    # data format
    fmtd = "{:>22,.13E}"

    # write data
    with open("b2fgmtry","w") as sgo:
        # length
        sgo.write("*cf:    int                2    nx,ny              \n")
        sgo.write("        " + str(nx) + "          " + str(ny)+"\n")

        # crx
        data = (len(x))
        sgo.write("*cf:    real               " + str(data) + \
                  "    crx             \n")

        lp = 0
        for i in range(len(x)):
            sgo.write(fmtd.format(x[i]))
            lp += 1
            # 6 data/line
            if lp == 6:
                lp =0
                sgo.write("\n")

        # cry
        data = len(y)
        sgo.write("\n*cf:    real               " + str(data) + \
                  "    cry             \n")

        lp = 0
        for i in range(len(y)):
            sgo.write(fmtd.format(y[i]))
            lp += 1
            # 6 data/line
            if lp == 6:
                lp =0
                sgo.write("\n")

        sgo.write("\n")

#-------------------------------------------------------------------------------
# formatting neutral data
# 1st argument:neutral data list
# 2nd argument:data namelist
# 3rd argumene:X length, int
# 4th argument:Y length, int
#-------------------------------------------------------------------------------
def S2SN(data, name, lx, ly):
    # setting data
    nx = lx
    ny = ly - 1
    ns = 1
    
    # data format
    fmt = "{:>22,.13E}"

    # write data
    with open("neutral","w") as ntl:
        # length
        ntl.write("*cf:    int                3    nx,ny,ns              \n")
        ntl.write("        " + str(nx) + "          " + str(ny) + "          " + str(ns)+"\n")

        # write neutral data
        for i in range(len(data)):
            cnt = 0
            ntl.write("*cf:    real               " + str(len(data[i])) + "    " + name[i] + "\n")
            for j in range(len(data[i])):
                if cnt == 6:
                    cnt = 0
                    ntl.write("\n")
                ntl.write(fmt.format(data[i][j]))
                cnt += 1
            ntl.write("\n")

#-------------------------------------------------------------------------------
# formatting impurity data
# 1st argument:impurity data list
# 2nd argument:data namelist
# 3rd argumene:X length, int
# 4th argument:Y length, int
#-------------------------------------------------------------------------------
def S2SI(data, name, lx, ly):
    # setting data
    nx = lx
    ny = ly - 1
    ns = 1

    # data format
    fmt = "{:>22,.13E}"

    # open files "w"mod
    file1 = open("impurity1", "w")
    file2 = open("impurity2", "w")

    # length
    file1.write("*cf:    int                3    nx,ny,ns              \n")
    file1.write("        " + str(nx) + "          " + str(ny) + "          " + str(ns)+"\n")
    file2.write("*cf:    int                3    nx,ny,ns              \n")
    file2.write("        " + str(nx) + "          " + str(ny) + "          " + str(ns)+"\n")
    
    # write imprity data1
    for i in range(len(data[0])):
        cnt = 0
        file1.write("*cf:    real               " + str(len(data[0][i])) + "    " + name[0][i] + "\n")
        for j in range(len(data[0][i])):
            if cnt == 6:
                cnt = 0
                file1.write("\n")
            file1.write(fmt.format(data[0][i][j]))
            cnt += 1
        file1.write("\n")

    # write imprity data2
    for i in range(len(data[1])):
        cnt = 0
        file2.write("*cf:    real               " + str(len(data[1][i])) + "    " + name[1][i] + "\n")
        for j in range(len(data[1][i])):
            if cnt == 6:
                cnt = 0
                file2.write("\n")
            file2.write(fmt.format(data[1][i][j]))
            cnt += 1
        file2.write("\n")
    
    # close files
    file1.close()
    file2.close()       

#-------------------------------------------------------------------------------
# switching data
# 1st argument:1D array
#-------------------------------------------------------------------------------
def switch_array(data):
    # data len
    ld = len(data)
    for i in range(int(ld/2)):
        dm = data[i]
        data[i] = data[ld-1-i]
        data[ld-1-i] = dm

    # return switched data
    return data

#-------------------------------------------------------------------------------
# nD array -> 1D array
# 1st argument:nD array
#-------------------------------------------------------------------------------
def create_data(data):
    cd = []
    nd = []

    # put together 3 area's(nD array -> 2D array(37*jcmax))
    for y in range(len(data[0])):    # 37
        cd.append([])
        for i in range(len(data)):   # 3
            for x in range(len(data[i][y])):
                cd[y].append(data[i][y][x])
                
    # 1D array
    for y in range(1, len(cd), 1):
        for x in range(len(cd[y])):
            # add data
            nd.append(cd[y][x])
    
    lx = len(cd[0]) 
    ly = len(cd) # 37

    # switch data
    nd = switch_array(nd)
    
    # return 1D array, X length, Y length
    return nd, lx, ly

#-------------------------------------------------------------------------------
# create geometry data
# 1st argument:nD array(geometry data)
#-------------------------------------------------------------------------------
def create_geometry(data):
    cd = []
    nd = []

    # put together 3 area's(nD array -> 2D array(37*jcmax))
    for y in range(len(data[0])): # 37
        cd.append([])
        for i in range(len(data)): # 3
            for x in range(len(data[i][y])):
                cd[y].append(data[i][y][x])
            # last line
            if i == len(data)-1:
                cd[y].append(data[i][y][x])
            
    # length
    ldy = len(cd) # 37
    ldx = len(cd[0])
    
    # point 1
    for y in range(ldy - 1, 0, -1):
        for x in range(ldx - 1, 0, -1):
            nd.append(cd[y][x])

    # point 2
    for y in range(ldy - 1, 0, -1):
        for x in range(ldx - 2, -1, -1):
            nd.append(cd[y][x])
    # point 4
    for y in range(ldy - 2, -1, -1):
        for x in range(ldx - 1, 0, -1):
            nd.append(cd[y][x])

    # point 3
    for y in range(ldy - 2, -1, -1):
        for x in range(ldx -2, -1, -1):
            nd.append(cd[y][x])

    # return 1D array
    return nd
#-------------------------------------------------------------------------------
# change eV -> J
# 1st argument:1D array
#-------------------------------------------------------------------------------
def eV2J(data):
    for i in range(len(data)):
        data[i] = data[i] / 6.242e18

    return data

# ======================================================================
# running the program
# ======================================================================
if __name__ == "__main__":
    print("start program \"sonic2imas.py\"\n")
    sonic2imas()
    print("\nend program \"sonic2imas.py\"")
