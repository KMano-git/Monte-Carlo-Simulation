# ======================================================================
# -*- name  : getData2.py
# -*- coding: utf-8 -*-
# ======================================================================

# import binary data reading library
#from os import write
from readDTNTL import readDTNTL
from readDTIMP import readDTIMP, readIMPSY
from makeDumy import makeDumy

import numpy as np

# --------------------------------------------------
# neutral and mesh
# 1st argument: path, string
# 2nd argument:radiation data name, stirng
# 3rd argument:data name, string
# 4th argument:mesh index data, integer
# 5th argument:parametar name, string
# 6th argument:index of nDarray, integer
# --------------------------------------------------
def DTNTL_MSH(fileDTNTL, ind, ncmax):

    [jcxp1, jcxp2, jcmax, icwl1, icwl2, nszx, nszy, iaxs] = ind # 30 91 120 1 37 160 80, 56or62
    # get neutral data
    # vmwork
    w0 = readDTNTL(fileDTNTL)[2]
    # creating plot data
    [X, Y, W, WW] = [[], [], [], []]

    [vmworks, wflx] = dep_vmwork(w0, ncmax, ind)

    for k in range(len(vmworks)):
        W = []
        # area1
        W.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addW=[]
            # loop
            for j in range(1, jcxp1-1+1): # jcxp1=30
                addW.append(vmworks[k][j-1][i-1])
            # last line
            addW.append(vmworks[k][jcxp1-1-1][i-1])
            addW.append(vmworks[k][jcxp1-1-1][i-1])
        
            # add the created column to the plot data
            W[0].append(addW)

        # --------------------
        # area2
        W.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addW=[]
            # loop
            for j in range(jcxp1-1+1, jcxp2-1-1+1): # jcxp2=91
                addW.append(vmworks[k][j-1][i-1])
            # last line
            addW.append(vmworks[k][jcxp2-1-2][i-1])
        
            # add the created column to the plot data
            W[1].append(addW)
    
        # --------------------
        # area3
        W.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addW=[]
            # loop
            for j in range(jcxp2-1, jcmax-2+1):
                addW.append(vmworks[k][j-1][i-1])
            # last line
            addW.append(vmworks[k][jcmax-2-1][i-1])
            addW.append(vmworks[k][jcmax-2-1][i-1])
        
            # add the created column to the plot data
            W[2].append(addW)
        WW.append(W)
        
    # return data list(wssn, wssp, wswe, wswi, wsbr, wden)
    return WW

# --------------------------------------------------
# deployment vmwork
# 1st argument: 
# 2nd argument:radiation data name, stirng
# --------------------------------------------------
def dep_vmwork(data, ncmax, ind):

    [jcxp1, jcxp2, jcmax, icwl1, icwl2, ndx, ndy, iaxs] = ind # 30 91 120 1 37 160 80, 56or62
    wflx = data[0]
    [wssn, wssp, wswe, wswi, wsbr, wden] = [[], [], [], [], [], []]
    dl = 6500 + 1
    # deployment data
    for i in range(dl):
        wssn.append(data[1+0*dl+i])
        wssp.append(data[1+1*dl+i])
        wswe.append(data[1+2*dl+i])
        wswi.append(data[1+3*dl+i])
        wsbr.append(data[1+4*dl+i])
        wden.append(data[1+5*dl+i])
    vm = [wssn, wssp, wswe, wswi, wsbr, wden]

    # sum data
    Vm = []
    for i in range(len(vm)):
        Vm.append([])
        for j in range(len(vm[i])):
            sum = 0
            for k in range(len(vm[i][j])):
                sum += vm[i][j][k]
            Vm[i].append(sum)

    nd = []
    VM = []
    # dumy array
    for i in range(0,ndx):
        nd.append([])            
        for j in range(0,ndy):
            nd[i].append(0)

    for k in range(len(Vm)):
        i = 0
        j = 0
        # 1D -> 2D
        e1m = icwl2*(jcxp1-1)
        e2m = e1m+(jcxp2-jcxp1-1)*(iaxs-1)
        for cnt in range(ncmax):
            if cnt >0:
                # area1
                if (cnt<e1m) and (cnt%icwl2 == 0):
                    i+=1
                    j=0
                # area2 
                elif (cnt>=e1m and cnt<e2m) and ((cnt-e1m)%(iaxs-1) == 0):
                    i+=1
                    j=0
                # area3
                elif (cnt>=e2m and cnt<ncmax) and ((cnt-e2m)%icwl2 == 0):
                    i+=1
                    j=0
            nd[i][j]=Vm[k][cnt]
            j+=1
        VM.append(nd)
        
    # return neutral datalist, wflx
    return VM, wflx

# --------------------------------------------------
# impurity and mesh
# 1st argument:impurity datafile path, string
# 2nd argument:mesh datafile path, string
# 3rd argument:data name, string
# 4th argument:mesh index data, integer
# 5th argument:index of nDarray, integer
# --------------------------------------------------
def DTIMP_MSH(fileIMP, ind, ncmax):

    [jcxp1, jcxp2, jcmax, icwl1, icwl2, nszx, nszy, iaxs] = ind # 30 91 120 1 37 160 80, 56or62
    # get imprity data
    # nsizp nsizs nsizc nsput nzmx ncmx fsput twrd tdnz csput twci
    t0 = readDTIMP(fileIMP)

    # twci not out
    if len(t0) == 10:
        # set dumy
        twci = makeDumy(6500, 0, 0)
    else:
        # set read data
        twci = t0[10]

    nzmx = t0[4][0]    # nzmx = 6 or 18

    # creating plot data
    [X, Y, T, TT] = [[], [], [], []]

    data = []
    data.append(trance_data(t0[7], ncmax, ind)) # twrd
    for i in range(nzmx+1):
        data.append(trance_data(t0[8][i], ncmax, ind))  # tdnz[0~ncmax]
    data.append(trance_data(twci, ncmax, ind)) # twci

    for k in range(len(data)):
        # area1
        T = []
        T.append([])
    
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            for j in range(1, jcxp1-1+1): # jcxp1=30
                addT.append(data[k][j-1][i-1])
            # last line
            addT.append(data[k][jcxp1-1-1][i-1])
            addT.append(data[k][jcxp1-1-1][i-1]) # dumy insert
        
            # add the created column to the plot data
            T[0].append(addT)

        # --------------------
        # area2
        T.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            for j in range(jcxp1-1+1, jcxp2-1-1+1): # jcxp2=91
                addT.append(data[k][j-1][i-1])
            # last line
            addT.append(data[k][jcxp2-1-2][i-1])
    
            # add the created column to the plot data
            T[1].append(addT)
    
        # --------------------
        # area3
        T.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            for j in range(jcxp2-1, jcmax-2+1):
                addT.append(data[k][j-1][i-1])
            # last line
            addT.append(data[k][jcmax-2-1][i-1])
            addT.append(data[k][jcmax-2-1][i-1])    # dumy insert
    
            # add the created column to the plot data
            T[2].append(addT)
        
        TT.append(T)
        
    # return data list(TT[0]:twrd,TT[1~nzmx]:tdnz[0~nzmx-1])
    return TT, nzmx

# --------------------------------------------------
# impurity and mesh
# 1st argument:impurity datafile path, string
# 2nd argument:mesh datafile path, string
# 3rd argument:data name, string
# 4th argument:mesh index data, integer
# 5th argument:index of nDarray, integer
# 6th argument:switch function flag(0->DTIMP, 1->IMPSY), integer
# --------------------------------------------------
def IMPSY_MSH(fileIMP, ind, ncmax):

    [jcxp1, jcxp2, jcmax, icwl1, icwl2, nszx, nszy, iaxs] = ind # 30 91 120 1 37 160 80, 56or62
    # get imprity data
    # nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput, tfrz, tthz, tvlz, tionZ, trecZ, traidZ, tradliZ, tradrZ
    t0 = readIMPSY(fileIMP)
    nzmx = t0[4][0]

    # creating plot data
    [X, Y, T, TT] = [[], [], [], []]

    data = []
    # 20220303 fix
    data.append(trance_data(t0[7], ncmax, ind)) # twrd
    for i in range(nzmx+1):
        data.append(trance_data(t0[8][i] , ncmax, ind)) # tdnz
    for i in range(nzmx+1):
        data.append(trance_data(t0[10][i], ncmax, ind)) # tfrz
    for i in range(nzmx+1):
        data.append(trance_data(t0[11][i], ncmax, ind)) # tthz
    for i in range(nzmx+1):
        data.append(trance_data(t0[13][i], ncmax, ind)) # tionZ 
    for i in range(nzmx+1):
        data.append(trance_data(t0[14][i], ncmax, ind)) # trecZ

    for k in range(len(data)):
        T = []
        # area1
        T.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            #addT.append(data[0][i-1])
            for j in range(1, jcxp1+1): # jcxp1=30
                addT.append(data[k][j-1][i-1])
            # last line
            addT.append(data[k][jcxp1][i-1]) #2022.03.08 SY omit faulty plot above Xp 
        
           # add the created column to the plot data
            T[0].append(addT)

        # --------------------
        # area2
        T.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            for j in range(jcxp1+1, jcxp2-1+1): # jcxp2=91
                addT.append(data[k][j-1][i-1])
            # last line SY 
            addT.append(data[k][jcxp1][i-1])
        
            # add the created column to the plot data
            T[1].append(addT)
    
        # --------------------
        # area3
        T.append([])
        for i in range(icwl1, icwl2+1): # icwl1=1, icwl2=37
            addT=[]
            # loop
            #addT.append(data[k][90][i-1])
            for j in range(jcxp2, jcmax+1):
                addT.append(data[k][j-1][i-1])
            # last line
            addT.append(data[k][jcmax-1][i-1])
        
            # add the created column to the plot data
            T[2].append(addT)
        
        TT.append(T)
        
    # return data list(TT[0]:twrd,TT[1~nzmx]:tdnz[0~nzmx-1],TT[nzmx+1~2*nzmx]:tionZ[0~nzmx-1],TT[2*nzmx+1~3*nzmx]:trecZ[0~nzmx-1])
    return TT, nzmx

# --------------------------------------------------
# 1D -> 2D
# 1st argument:1D array, double
# 2nd argument:max len, integer
# --------------------------------------------------
def trance_data(data, ncmax, ind):
    nd=[]
    cnt=0
    [jcxp1, jcxp2, jcmax, icwl1, icwl2, ndx, ndy, iaxs] = ind # 30 91 120 1 37 160 80, 56or62
    # create dumy array
    for i in range(ndx):
        nd.append([])
        for j in range(ndy):
            nd[i].append(0)       
    
    for j in range(icwl2):
        nd[0][j] = 0.0
        nd[jcmax-1][j] = 0.0
    i = 1
    j = 0
    # 1D -> 2D
    e1m = icwl2*(jcxp1-1)
    e2m = e1m+(jcxp2-jcxp1-1)*(iaxs-1)
    # loop ncmax
    for cnt in range(ncmax):
        if cnt >0:
            # area1
            if (cnt<e1m) and (cnt%icwl2 == 0):
                i+=1
                j=0
            # area2
            elif (cnt>=e1m and cnt<e2m) and ((cnt-e1m)%(iaxs-1) == 0):
                i+=1
                j=0
            # area3
            elif (cnt>=e2m and cnt<ncmax) and ((cnt-e2m)%icwl2 == 0):
                i+=1
                j=0
        nd[i][j]=data[cnt]
        j+=1
        
    #  return implty data(2D array)
    return nd

# ======================================================================
# running the program
# ======================================================================
if __name__ == "__main__":

    # use function getData
    from getData import getEnv
    
    # get data
    env = getEnv()
    
    # setting data path
    fileDTNTL  = env[0][3]  # DTNTL
    fileDTIMP1 = env[0][4]  # DTIMP1
    fileDTIMP2 = env[0][5]  # DTIMP2
    fileIMPSY1 = env[0][6]  # IMPSY1
    fileIMPSY2 = env[0][7]  # IMPSY2
    
    # program description
    print('getData2.py   2021/08/24')
    print('--------------------')
    #print(' > function [DTNTL] get neutralData')
    #print('   , using function [readDTNTL.py readDTNTL]')
    #print(' > function [DTIMP] get impurityData')
    #print('   , using function [readDTIMP.py readDTIMP]')
    #print(' > function [IMPSY] get impurityData')
    #print('   , using function [readDTIMP.py readIMPSY]')
    print()
    
    # faile path is NULL
    if (fileDTNTL=='') | (fileDTIMP1=='') | (fileIMPSY1=='') | (fileDTIMP2=='') | (fileIMPSY2==''):
        # | (fileDTWRD==''):
        print('Paths undefined in path.txt.')
        print('sys exit.')
        sys.exit()

    w=DTNTL_MSH(fileDTNTL, (30, 91, 120, 1, 37), 0, 5806) # vmwork[0]
    #t=DTIMP_MSH(fileDTIMP1, (30, 91, 120, 1, 37), 5806)       # DTIMP1
    #t=IMPSY_MSH(fileIMPSY1, (30, 91, 120, 1, 37), 5806)       # IMPSY1
    #print(len(t))
    #print(len(t[0]))
    print('======================================================================')
