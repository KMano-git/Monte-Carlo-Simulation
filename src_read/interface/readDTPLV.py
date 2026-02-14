# ======================================================================
# -*- name  : readDTPLV.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import numpy as np

# use function readBinary
from readBinary import getBinary, readBinary, resize

# ======================================================================
# readDTPLV(fileName):
# ----------------------------------
# read prasma data
# ======================================================================
def readDTPLV(fileName):
    # ======================================================================
    # Binary data partial extraction
    # ======================================================================
    f = open(fileName, "rb")
    # Extract the original data part
    [data, headerList, footerList]=getBinary(f, len(f.read()))
    f.close()

    # -----------------------------------------------
    # Variable definition
    index = 0

    # Let f be the data extracted from the binary file
    f=data
    
    # -----------------------------------------------
    [[nszx, nszy, nszs],index] = readBinary(f, index, headerList[0], '>i')
    
    [ndx, ndy, ndsp]=[nszx, nszy, nszs]
    
    # vna   ( ndx, ndy, ndsp )
    [b,index] = readBinary(f, index, 8*ndx*ndy*ndsp, '>d')
    vna = resize(b,(ndx, ndy, ndsp))
    
    # vne   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vne = resize(b,(ndx, ndy))

    # vni   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vni = resize(b,(ndx, ndy))

    # vnezef( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vnezef = resize(b,(ndx, ndy))

    # vzf   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vzf = resize(b,(ndx, ndy))

    # vva   ( ndx, ndy, ndsp )
    [b,index] = readBinary(f, index, 8*ndx*ndy*ndsp, '>d')
    vva = resize(b,(ndx, ndy, ndsp))

    # vve   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vve = resize(b,(ndx, ndy))

    # vti   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vti = resize(b,(ndx, ndy))

    # vte   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vte = resize(b,(ndx, ndy))

    # vcs   ( ndx, ndy       )
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    vcs = resize(b,(ndx, ndy))

    # vea   ( ndx, ndy, ndsp )
    [b,index] = readBinary(f, index, 8*ndx*ndy*ndsp, '>d')
    vea = resize(b,(ndx, ndy, ndsp))
    
    # ---------------------------------------------
    # return all of read data
    
    return nszx, nszy, nszs, vna, vne, vni, vnezef, vzf, vva, vve, vti, vte, vcs, vea

# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':
    
    from getData import getEnv
    # get information
    env = getEnv()

    # setting path
    fileDTPLV = env[0][0]
    
    # program description
    print('readDTPLV.py')
    print('--------------------')
    print(' > function [readDTPLV] read DTPLV from datafile')
    print('   , using function [readBinary.py getBinary & readBinary & resize]')
    print('--------------------')
    print('[test] read this file ...')
    print(fileDTPLV.format('01'))
    print()
    
    # read data test
    data=readDTPLV(fileDTPLV)

    # print of read data
    for i in range(0, len(data)):
        
        if (data[i].size==1):
            print('No.{}: Value={}'.format(i+1, data[i]))
        else:
            print('No.{}: Array{}'.format(i+1, data[i].shape))

    print('======================================================================')

