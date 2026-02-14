# ======================================================================
# -*- name  : readDTNTL.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import numpy as np

# use function readBinary
from readBinary import getBinary, readBinary, resize

#-------------------------------------------------------------------------------
# reading data for DTNTL
# 1st argument:file path
#-------------------------------------------------------------------------------
def readDTNTL(fileName):
    # Binary data partial extraction
    # And extract the original data part
    with open(fileName,"rb") as f:
        [data, headerList, footerList] = getBinary(f, len(f.read()))
        # print(headerList) # [4224365] or [4288685]

    # -----------------------------------------------
    # Variable definition
    index = 0

    # Let f be the data extracted from the binary file
    ld=data

    # vtime
    [b, index] = readBinary(ld, index, 8, '>d')
    vtime = b
    
    # vflux
    [b, index] = readBinary(ld, index, 64, '>d')
    vflux = resize(b, (1, 8))
    
    # vmwork
    [b, index] = readBinary(ld, index, 8*89342*6, '>d')
    vmwork = resize(b, (89342, 6))

    # vitim
    [b, index] = readBinary(ld, index, 4, '>i')
    vitim = b
    
    # vitnt
    [b, index] = readBinary(ld, index, 4, '>i')
    vitnt = b
    
    # vnsrc
    [b, index] = readBinary(ld, index, 4, '>i')
    vnsrc = b
    
    # vsmty
    [b, index] = readBinary(ld, index, 4, '>i')
    vsmty = b
    
    # vsmno
    [b, index] = readBinary(ld, index, 4, '>i')
    vsmno = b
    
    # visrc
    [b, index] = readBinary(ld, index, 4*8, '>i')
    visrc = b
    
    # vkflx
    [b, index] = readBinary(ld, index, 4*8, '>i')
    vkflx = resize(b, (1, 8))
    
    # vksty
    [b, index] = readBinary(ld, index, 4*8, '>i')
    vksty = resize(b, (1, 8))
    
    # vnsmp
    [b, index] = readBinary(ld, index, 4*8, '>i')
    vnsmp = resize(b, (1, 8))
    
    # vcsrc
    [b, index] = readBinary(ld, index, 8*6, '>c')
    # vcsrc = resize(b, (8, 6))
    vcsrc = np.reshape(b, (8, 6))
    
    # cntmnt_emrk
    [b, index] = readBinary(ld, index, 1*1, '>c')
    cntmnt_emrk = b

    return vtime, vflux, vmwork, vitim, vitnt, vnsrc, vsmty, vsmno, visrc, vkflx, vksty, vnsmp, vcsrc, cntmnt_emrk

# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':
    

    # setting path
    # fileDTNTL = "./inp/DTNTL"
    fileDTNTL = "./sonicV4/src/interface/input/DTNTL"
    
    # program description
    print('readDTNTL.py')
    print('--------------------')
    print(' > function [readDTNTL] read DTNTL from datafile')
    print('   , using function [readBinary.py getBinary & readBinary & resize]')
    print('--------------------')
    print('[test] read this file ...')
    print(fileDTNTL.format('01'))
    print()
    
    # read data test
    data=readDTNTL(fileDTNTL)

    with open("dw_DTNTL","w") as f:
        
    # print of read data
        for i in range(0, len(data)):
            if (data[i].size==1):
                f.write('No.{}: Value={}'.format(i+1, data[i])+"\n")
            else:
                f.write('No.{}: Array{}'.format(i+1, data[i].shape)+"\n")

    print('======================================================================')
