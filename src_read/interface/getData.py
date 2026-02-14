# ======================================================================
# -*- name  : getData.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import os
import yaml
import numpy as np
from math import log

# import binary data reading library
from readDTPLV import readDTPLV
from readMSH import readMSH

# import setting file reading library
import inpfig as inpfig

# --------------------------------------------------
# plasma data
# --------------------------------------------------
def DTPLV(fileDTPLV, vName):

    logFlg=False
    if 'log' in vName:
        logFlg=True
        vName=vName.replace('log','')
        vName=vName.replace('(','')
        vName=vName.replace(')','')
    
    # reading plasma data
    [nszx, nszy, nszs, vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea]=readDTPLV(fileDTPLV)

    # returns the variable of vName
    vList = [vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea]
    nameList = ['vna', 'vne', 'vni', 'vnezef', 'vzf', 'vva', 'vve', 'vti', 'vte', 'vcs', 'vea']
    
    v=[]
    for i in range(0, len(vList)):
        # remove whitespace
        vName = vName.replace(' ','')
        nameList[i] = nameList[i].replace(' ','')
        
        # argument and name match
        if (vName == nameList[i]):
            v=vList[i]
            break
    
    return nszx, nszy,nszs, v, logFlg

# --------------------------------------------------
# mesh data
# --------------------------------------------------
def MSH(fileMSH, nszx, nszy):
    
    # reading mesh data
    [cver, nndq, nndp, cgdcom, nndx, nndy, cplmet]=readMSH(fileMSH, nszx, nszy)
    
    # grid data
    # data[3]
    [grdx, grdy] = cgdcom[17:19] #17~18

    # grid data
    # data[6]
    [hgdx, hgdy] = cplmet[39:41] #39~40
    [jcxp1, jcxp2, jcmax] = cplmet[5:8] #5~7
    icwl1= cplmet[8] #8
    icwl2= cplmet[10] #10
    [kgdx, kgdy] = cplmet[25:27] #25~26

    # print(jcxp1,jcxp2,jcmax,icwl1,icwl2) # 30 91 120 1 37

    return hgdx, hgdy, grdx, grdy, jcxp1, jcxp2, jcmax, icwl1, icwl2, kgdx, kgdy

# --------------------------------------------------
# extraction of non-zero elements
# --------------------------------------------------
# 2021-02-22 not used
def NonZero(hgdx, hgdy, dtplv, logFlg):
    
    [X,  Y,  V ] = [[], [], []]# record a group of multiple lines in multiple parts
    [xs, ys, vs] = [[], [], []]# a group of multiple lines
    [x,  y , v ] = [[], [], []]# coordinate sequence of each line
    
    beforeLen=0
    for i in range(0, hgdx.shape[0]):
        for j in range(0, hgdx.shape[1]):
            if ((hgdx[i][j]!=0) & (hgdy[i][j]!=0) & (dtplv[i][j]!=0)):
                x.append(float(hgdx[i][j]))
                y.append(float(hgdy[i][j]))
                
                # check logFlg
                if logFlg:
                    v.append(float(log(dtplv[i][j],10)))
                else:
                    v.append(float(dtplv[i][j]))
        
            
        # if the number of coordinate sequences of each line is different from the previous one, separate it there
        if((i!=0) & (beforeLen!=len(x))):
            X.append(xs)
            Y.append(ys)
            V.append(vs)
            [xs, ys, vs] = [[], [], []]# array initialization
        beforeLen=len(x)
        
        # insert x,y,v into xs,ys,vs
        if((len(x)!=0) & (len(y)!=0) & (len(v)!=0)):
            xs.append(x)
            ys.append(y)
            vs.append(v)
            [x,  y , v ] = [[], [], []]# array initialization
    
    return X, Y, V

# ======================================================================
# function to create plot data
# ======================================================================
# --------------------------------------------------
# plasma & mesh data
# --------------------------------------------------
#def DTPLV_MSH(fileDTPLV, fileMSH, vName, vIndex):
def DTPLV_MSH(fileDTPLV, fileMSH):

    # gets the specified variable plasma data
    pls = readDTPLV(fileDTPLV)
    [nszx, nszy, nszs, vna0, vne, vni, vnezef, vzf, vva0, vve, vti, vte, vcs, vea0] = pls

    [vna, vva, vea] = [[], [], []]
    for i in range(nszx):
        vna.append([])
        vva.append([])
        vea.append([])
        for j in range(nszy):
            vna[i].append(vna0[i][j][0])
            vva[i].append(vva0[i][j][0])
            vea[i].append(vea0[i][j][0])
    
    # get mesh data
    msh = readMSH(fileMSH, nszx, nszy)
    # grid data
    # cgdcom
    [grdx, grdy] = msh[3][17:19] #17~18

    # grid data
    # cplmet
    [hgdx, hgdy] = msh[6][39:41] #39~40
    [jcxp1, jcxp2, jcmax] = msh[6][5:8] #5~7
    icwl1= msh[6][8] #8
    icwl2= msh[6][10] #10
    [kgdx, kgdy] = msh[6][25:27] #25~26

    # add 2021.11.29
    ind = []
    ind.append(jcxp1) # 30 div range
    ind.append(jcxp2) # 91  div range + core range
    ind.append(jcmax) # 120 div + core + div
    ind.append(icwl1) # 1 
    ind.append(icwl2) # 37
    ind.append(nszx)  # 160
    ind.append(nszy)  # 80

    # creating plot data
    # --------------------------------------------------
    [X, Y, V, VV] = [[], [], [], []]

    data = [vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea]

    # --------------------
    for k in range(len(data)):
        [X, Y, V] = [[], [], []]
        # area1
        X.append([])
        Y.append([])
        V.append([])
        for i in range(icwl1, icwl2+1):
            addX=[]
            addY=[]
            addV=[]
            # loop
            for j in range(1, jcxp1+1):
                addX.append(grdx[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addY.append(grdy[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addV.append(data[k][j-1][i-1])
           # TODO: above definition of grid position is not the way that we have intended
           # last line
            addX.append(grdx[kgdx[jcxp1-1,i-1,2]-1,kgdy[jcxp1-1,i-1,2]-1])
            addY.append(grdy[kgdx[jcxp1-1,i-1,2]-1,kgdy[jcxp1-1,i-1,2]-1])
            addV.append(data[k][jcxp1][i-1])
        
           # add the created column to the plot data
            X[0].append(addX)
            Y[0].append(addY)
            V[0].append(addV)
    
        # --------------------
        # area2
        X.append([])
        Y.append([])
        V.append([])
        for i in range(icwl1, icwl2+1):
            addX=[]
            addY=[]
            addV=[]
            # loop
            for j in range(jcxp1+1, jcxp2-1+1):
                addX.append(grdx[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addY.append(grdy[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addV.append(data[k][j-1][i-1])
            # last line
            addX.append(grdx[kgdx[jcxp2-2,i-1,2]-1,kgdy[jcxp2-2,i-1,2]-1])
            addY.append(grdy[kgdx[jcxp2-2,i-1,2]-1,kgdy[jcxp2-2,i-1,2]-1])
            addV.append(data[k][jcxp1][i-1])
        
           # add the created column to the plot data
            X[1].append(addX)
            Y[1].append(addY)
            V[1].append(addV)
    
        # --------------------
        # area3
        X.append([])
        Y.append([])
        V.append([])
        for i in range(icwl1, icwl2+1):
            addX=[]
            addY=[]
            addV=[]
            # loop
            for j in range(jcxp2, jcmax+1):
                addX.append(grdx[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addY.append(grdy[kgdx[j-1,i-1,3]-1,kgdy[j-1,i-1,3]-1])
                addV.append(data[k][j-1][i-1])
            # last line
            addX.append(grdx[kgdx[jcmax-1,i-1,2]-1,kgdy[jcmax-1,i-1,2]-1])
            addY.append(grdy[kgdx[jcmax-1,i-1,2]-1,kgdy[jcmax-1,i-1,2]-1])
            addV.append(data[k][jcmax-1][i-1])
        
            # add the created column to the plot data
            X[2].append(addX)
            Y[2].append(addY)
            V[2].append(addV)
        
        VV.append(V)
        
    return X, Y, VV, ind

# --------------------------------------------------
# wall data
# --------------------------------------------------
def WALL(fileName):
    f=open(fileName)
    data = f.read()
    f.close()

    [xs, ys, vs] = [[], [], []]# a group of multiple lines
    [x,  y , v ] = [[], [], []]# coordinate sequence of each line
    
    lines = data.split('\n')
    for i in range(0,len(lines)):
        line = lines[i].split()
        if(len(line) == 2):
            x.append(float(line[0]))
            y.append(float(line[1]))
        elif((len(x)!=0) & (len(y)!=0)):
            xs.append(x)
            ys.append(y)
            [x,  y , v ] = [[], [], []]# array initialization
    
    # It cannot be an array because the number of columns in the first and second rows is different
    # list as return
    
    return xs, ys

# --------------------------------------------------
# setting
# --------------------------------------------------
def INPFIG(cNo, upcntv):
    # graph display area
    csizText = inpfig.upfgsz.csiz       # text to display in the selector
    memo_xy = csizText[cNo].split(':')  # divide the string into notes and coordinates
    xy = memo_xy[1].split(',')          # divide the coordinate string into left,right,bottom,top
    csiz = [float(xy[0]),float(xy[1]),float(xy[2]),float(xy[3])]# extract as a numerical value
    
    # Physical quantity : Log display corrected
    # Choice is log(), but inpfig is defined by log_ notation
    if 'log' in upcntv:
        upcntv = upcntv.replace('log(','')
        upcntv = upcntv.replace(')','')
        upcntv = upcntv.replace(' ','')
        upcntv = 'log_'+upcntv
    
    # Color bar delimiter [v0,v1,v2, ... ,vN]
    bounds = []# initial value
    for key, value in inpfig.upcntv.__dict__.items():
        if (upcntv == key):
            bounds = value
    
    # return value
    return csiz, csizText, bounds

# ======================================================================
# get the information you need
# ======================================================================
def getYAML():
    
    # --------------------------------------------------
    # get login name
    loginUser=''
    try:
        # Windows
        loginUser = str(os.environ['USER'])
    except:
        # Linux
        loginUser = str(os.environ['USERNAME'])
    else:
        pass
    
    # get execution directory (absolute path of this program)
    workDir = os.path.dirname(os.path.abspath(__file__))
    workDir = workDir.replace('\\\\','/')
    workDir = workDir.replace('\\','/')
    
    # --------------------------------------------------
    # setting path
    # --------------------------------------------------
    # set the path of each data

    # common record (reference)
    pathTxt = workDir+'/config.yaml'
        
    # personal record (recording destination) *one level higher
    pathTxt_user=os.path.dirname(workDir)+'/config_'+loginUser+'.yaml'
        
    # if you have personal records, see personal records
    if os.path.exists(pathTxt_user):
        pathTxt=pathTxt_user
        
    f=open(pathTxt, 'r')
    data = f.read()
    f.close()
    yamlcon = yaml.safe_load(data)
    
    # --------------------------------------------------
    # *If you want to add a physical quantity option, add it to inpfig.py's upcntv class variable as well

    # Physical quantity choice
    nameList = ['vni', 'vne', 'vti', 'vte', 'vva', 'vna', 'vnezef', 'vzf', 'vve',  'vcs', 'vea']
    
    return yamlcon, loginUser


# ======================================================================
# get the information you need
# ======================================================================
def getEnv():
    
    # --------------------------------------------------
    # get login name
    loginUser=''
    try:
        # Windows
        loginUser = str(os.environ['USER'])
    except:
        # Linux
        loginUser = str(os.environ['USERNAME'])
    else:
        pass
    
    # get execution directory (absolute path of this program)
    workDir = os.path.dirname(os.path.abspath(__file__))
    workDir = workDir.replace('\\\\','/')
    workDir = workDir.replace('\\','/')
    
    # --------------------------------------------------
    # setting path
    try:
        # --------------------------------------------------
        # set the path of each data

        # common record (reference)
        pathTxt = workDir+'/path.txt'
        
        # personal record (recording destination) *one level higher
        pathTxt_user=os.path.dirname(workDir)+'/path_'+loginUser+'.txt'
        
        # if you have personal records, see personal records
        if os.path.exists(pathTxt_user):
            pathTxt=pathTxt_user
        
        f=open(pathTxt)
        data = f.read()
        f.close()
        
        lines = data.split('\n')
        
        # unification of path
        for i in range(0, len(lines)):
            # \\ -> \
            lines[i]=lines[i].replace('\\\\','\\')
            # \ -> /
            lines[i]=lines[i].replace('\\','/')
        
        # add 2021.10.18
        for i in range(len(lines)-1):
            if (os.path.exists(lines[i]))==False:
                print('Paths undefined in path.txt.')
                print('sys exit.')
                sys.exit()
        
    except:
        # no record
        lines = ''
    else:
        pass
    
    # --------------------------------------------------
    # *If you want to add a physical quantity option, add it to inpfig.py's upcntv class variable as well

    # Physical quantity choice
    nameList = ['vni', 'vne', 'vti', 'vte', 'vva', 'vna', 'vnezef', 'vzf', 'vve',  'vcs', 'vea']
    
    return lines, loginUser

# ======================================================================
# output plotdata
# ======================================================================
def outputXYV(X, Y, V, vNamee, outFile):
    # hgdx, hgdy, output plasma data
    for h in range(0, len(X)):
        x = np.array(X[h])
        y = np.array(Y[h])
        v = np.array(V[h])
        s=''
        # change open() and close() to with open() 2021/03/29
        with open(outFile,"w") as f:
            for i in range(0, x.shape[0]):
                for j in range(0, x.shape[1]):
                    s+='  '+'{:.10E}'.format(x[i][j]).rjust(18)+'{:.10E}'.\
                            format(y[i][j]).rjust(18)+'{:.10E}'.format(v[i][j]).rjust(18)+'\n'
                    s+=' \n'
            f.write(s)
    # successful
    print(outFile + " was created.")


# ======================================================================
# running the program
# ======================================================================
if __name__ == "__main__":
    
    # get data
    env = getEnv()
    
    # setting data path
    fileDTPLV = env[0][0]
    fileMSH   = env[0][1]
    fileWALL  = env[0][2]
    
    # program description
    print('getData.py   2021/05/10')
    print('--------------------')
    print(' > function [DTPLV_MSH] get plotData (X,Y,V) , is Not0')
    print('   , using function [DTPLV] and [MSH] and [NonZero]')
    print('       function [DTPLV], using function [readDTPLV.py readDTPLV]')
    print('       function [MSH], , using function [readMSH.py readMSH]')
    print(' > function [WALL] get plotData (xs,ys) ')
    print(' > function [INPFIG] get plotInfo from inpfig.py ')
    print(' > function [getEnv] get Info')
    print()
    
    if ((fileDTPLV=='') | (fileMSH=='') | (fileWALL=='')):
        print('Paths undefined in path.txt.')
        print('sys exit.')
        sys.exit()
    
    # read wall data
    [wx, wy] = WALL(fileWALL)
    print('wallData xy -- OK ')
    print()
    
    # read setting file
    plotsize=input('plotsize >>>')
    vName=input('upcntv >>>')
    vIndex=int(input('vIndex >>>'))
    [csiz, csizText, bounds] = INPFIG(int(plotsize),vName)
    print()
    
    # get data
    [X, Y, V, ind] = DTPLV_MSH(fileDTPLV, fileMSH)[:4]
    print(ind)
    print(fileDTPLV)
    
    # V == [], plasma data is Invalid
    if(V == []):
        print('input [V] error')
    else:
        # output text data
        print('output TextData ? ')
        outFlg=int(input('0:No or 1:Yes >>>'))
        print()
        if outFlg==1:
            out = input('outDir or exit >>>')
            if out=='exit':
                pass
            else:
                outFile = out+'/fort_'+str(vName)+'.txt'
                outputXYV(X, Y, V, vName, outFile)
        
    print('======================================================================')
