# ======================================================================
# -*- name  : readDTIMP.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import numpy as np

# use function readBinary
from readBinary import getBinary, readBinary, resize

#-------------------------------------------------------------------------------
# reading data for DTIMP's
# 1st argument:file path
#-------------------------------------------------------------------------------
def readDTIMP(fileName):
    # Binary data partial extraction
    # And extract the original data part
    with open(fileName,"rb") as f:
        [data, headerList, footerList] = getBinary(f, len(f.read()))
        # print(headerList) # [12, 3952072, 52000]

    # -----------------------------------------------
    # Variable definition
    index = 0

    # Let f be the data extracted from the binary file
    f=data

    #
    [[nsizp, nsizs, nsizc],index] = readBinary(f, index, headerList[0], '>i')
    # print(index) # 12

    [ndmp, ndis, ndmc] = [nsizp, nsizs, nsizc]

    # nsput
    [b,index] = readBinary(f, index, 4, '>i')
    nsput = b
    # print(index) # 16

    # nzmx
    [b,index] = readBinary(f, index, 4, '>i')
    nzmx = b
    # print(index) # 20

    # ncmx
    [b,index] = readBinary(f, index, 4, '>i')
    ncmx = b
    # print(index) # 24

    # fsput
    [b,index] = readBinary(f, index, 8*5, '>d')
    fsput = resize(b,1*5)
    # print(index) # 64

    # twrd
    [b,index] = readBinary(f, index, 8*ndmc, '>d')
    twrd = resize(b,1*ndmc)
    # print(index) # 52064

    # tdnz
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tdnz = resize(b,((ndis+1), ndmc))
    # print(index) # 3952064

    # csput
    [b,index] = readBinary(f, index, 4*5, '>c')
    csput = np.reshape(b,(5,4))
    # print(index) # 3592084

    # twci
    try:
        [b,index] = readBinary(f, index, headerList[2], '>d')
        twci = b
        return nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput, twci
    
    except:
        print("twci couldn't read.")

        return nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput

#-------------------------------------------------------------------------------
# reading data for IMPSY's
# 1st argument:file path
#-------------------------------------------------------------------------------
def readIMPSY(fileName):
    # Binary data partial extraction
    # And extract the original data part
    with open(fileName,"rb") as f:
        [data, headerList, footerList] = getBinary(f, len(f.read()))
        # print(headerList) # [12, 3952072, 11700000, 7800000, 11700000]
        
    # -----------------------------------------------
    # Variable definition
    index = 0

    # Let f be the data extracted from the binary file
    f=data

    # ndmp, ndis, ndmc
    [[nsizp, nsizs, nsizc],index] = readBinary(f, index, headerList[0], '>i')
    # print(index) # 12

    [ndmp, ndis, ndmc] = [nsizp, nsizs, nsizc]

    # nsput
    [b,index] = readBinary(f, index, 4, '>i')
    nsput = b

    # nzmx
    [b,index] = readBinary(f, index, 4, '>i')
    nzmx = b

    # ncmx
    [b,index] = readBinary(f, index, 4, '>i')
    ncmx = b

    # fsput
    [b,index] = readBinary(f, index, 8*5, '>d')
    fsput = resize(b,1*5)

    # twrd
    [b,index] = readBinary(f, index, 8*ndmc, '>d')
    twrd = resize(b,1*ndmc)

    # tdnz
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tdnz = resize(b,((ndis+1), ndmc))

    # csput
    [b,index] = readBinary(f, index, 4*5, '>c')
    csput = np.reshape(b,(5,4))
    # print(index) # 3952084

    # tfrz
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tfrz = resize(b,((ndis+1), ndmc))
    # print(index) # 7852084

    # tthz
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tthz = resize(b,((ndis+1), ndmc))
    # print(index) # 11752084

    # tvlz
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tvlz = resize(b,((ndis+1), ndmc))
    # print(index) # 15652084

    # tionZ
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tionZ = resize(b,((ndis+1), ndmc))
    # print(index) # 19552084

    # trecZ
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    trecZ = resize(b,((ndis+1), ndmc))
    # print(index) # 23452084

    # tradiZ
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tradiZ = resize(b,((ndis+1), ndmc))
    # print(index) # 27352084

    # tradliZ
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tradliZ = resize(b,(75, ndmc))
    # print(index) # 31252084

    # tradrZ
    [b,index] = readBinary(f, index, 8*(ndis+1)*ndmc, '>d')
    tradrZ = resize(b,((ndis+1), ndmc))
    # print(index) # 3515208

    return nsizp, nsizs, nsizc, nsput, nzmx, ncmx, fsput, twrd, tdnz, csput, tfrz, tthz, tvlz, tionZ, trecZ, tradiZ, tradliZ, tradrZ

# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':
    
    # use function getData
    from getData import getEnv

    # get data
    env = getEnv()
    
    # setting path
    #fileDTIMP = env[0][4] # DTIMP1
    fileDTIMP = env[0][5] # DTIMP2
    #fileIMPSY = env[0][6] # IMPSY1
    fileIMPSY = env[0][7] # IMPSY2
    
    # program description
    print('readDTIMP.py')
    print('--------------------')
    print(' > function [readDTIMP] and [readIMPSY] read DTIMP\'s or IMPSY\'s from datafile')
    print('   , using function [readBinary.py getBinary & readBinary & resize]')
    print('--------------------')
    print('[test] read this file ...')
    print(fileDTIMP.format('01'))
    print(fileIMPSY.format('01'))
    print()
    
    # read data test
    #data=readDTIMP(fileDTIMP)
    data=readIMPSY(fileIMPSY)

    #with open("dw_DTIMP1","w") as f:
    #with open("dw_DTIMP2","w") as f:
    #with open("dw_IMPSY1","w") as f:
    with open("dw_IMPSY2","w") as f:
        
    # print of read data
        for i in range(0, len(data)):
            if (data[i].size==1):
                f.write('No.{}: Value={}'.format(i+1, data[i])+"\n")
            else:
                f.write('No.{}: Array{}'.format(i+1, data[i].shape)+"\n")

    print('======================================================================')
