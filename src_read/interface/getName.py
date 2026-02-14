# ======================================================================
# -*- name  : getName.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import os

#-------------------------------------------------------------------------------
# get device, version
# 1st argument:loggin user name, string
#-------------------------------------------------------------------------------
def getName(user):
    
    dbpath = '/home/' + user + '/public/imasdb'
    #dbpath = '../../' + user + '/public/imasdb/' # debug

    # check path
    if os.path.exists(dbpath) == False:
        print("ERROR : DBpath \"", dbpath, "\" is not exist!")
        sys.exit()

    # DBname list
    dev = os.listdir(dbpath)

    devlist=[]
    for d in dev:
        devlist.append(os.path.join(dbpath,d))

    ver=[]
    for i in range(len(devlist)):
        ver.append([])
        ver[i] = os.listdir(devlist[i])

    return dev, ver

#-------------------------------------------------------------------------------
# get element name
# 1st argment:element No, integer
#-------------------------------------------------------------------------------
def getElem(nzmax):

    if nzmax < 1 or 74 < nzmax:
        print("ERROR: Z number", nzmax, "is not handled!")
        sys.exit()

    # create element list
    name  = ["H","He","Li","Be","B","C","N","O","F","Ne"]
    name += ["Na","Mg","Al","Si","P","S","Cl","Ar"]
    name += ["K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"]
    name += ["Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pb","Ag","Cb","In","Sn","Sb","Te","I","Xe"]
    name += ["Cs","Ba"]
    name += ["La","Ce","Pr","Nb","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"]
    name += ["Hf","Ta","W"]

    return name[nzmax-1]

# ======================================================================
# running the program
# ======================================================================
if __name__ == "__main__":
    # inport
    from getData import getEnv
    
    print("start program \"getName.py\"\n")
    user = getEnv()[1]
    data = getName(user)
    print(data[0])
    print(data[1])

    nzmax = 18
    print(getElem(nzmax))
    nzmax = 6
    print(getElem(nzmax))
    print("\nend program \"getName.py\"")
