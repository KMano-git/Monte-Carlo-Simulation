# ======================================================================
# -*- name  : readMSH.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import sys
import numpy as np

# use function readBinary
from readBinary import getBinary, readBinary, resize

# ======================================================================
# readMSH(fileName, ndx, ndy):
# ----------------------------------
# read mesh data
# ======================================================================
def readMSH(fileName, ndx, ndy):
    # ======================================================================
    # Binary data partial extraction
    # ======================================================================
    f = open(fileName, "rb")
    # Extract the original data part
    [data, headerList, footerList]=getBinary(f, len(f.read()))
    f.close()
    # print(headerList)
    # 2021: 80 8 513988 8 3229540
    # 2020: 80 8 513988 8 3229540

    # -----------------------------------------------
    # Variable definition
    # [ndq, ndp]=[ndx, ndy]
    index = 0
    
    # Let f be the data extracted from the binary file
    f=data
    
    # -----------------------------------------------
    # character : cver
    [cverS,index] = readBinary(f, index, headerList[0], '>c')

    # Recreate the string
    cver=''
    for s in cverS:
        cver+=str(s.decode())
    # -----------------------------------------------
    
    # integer : nndq, nndp
    # [[nndq, nndp],index] = readBinary(f, index, 4*2, '>i')
    [[ndq, ndp],index] = readBinary(f, index, 4*2, '>i')
    # print(ndq,ndp)
    # 160 80
    # 160 80

    # Defined in cgdcom integer:
    #     >  mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx
    #     >  ,mpw1, mpsp, mpw2, mpax, npmx
    #     >  ,mpsol, mpprv, mpman
    # Defined in cgdcom real*8:
    #     >  ,grdx(ndq,ndp), grdy(ndq,ndp), hbr(ndq,ndp), hbz(ndq,ndp), hbt(ndq,ndp)
    #     >  ,pssol(ndp), psprv(ndp), psman(ndp)
    [[mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx,mpw1, mpsp, mpw2, mpax, npmx, mpsol, mpprv, mpman], index] = readBinary(f, index, 4*17, '>i')
    # print(mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx,mpw1, mpsp, mpw2, mpax, npmx, mpsol, mpprv, mpman)
    # 1 30 31 91 53 70 92 121 121 1 31 38 62 62 30 7 32
    # 1 31 32 92 62 92 93 123 123 1 20 38 56 56 0 0 0

    [b,index] = readBinary(f, index, 8*ndq*ndp, '>d')
    grdx = resize(b,(ndq,ndp))

    [b,index] = readBinary(f, index, 8*ndq*ndp, '>d')
    grdy = resize(b,(ndq,ndp))

    [b,index] = readBinary(f, index, 8*ndq*ndp, '>d')
    hbr = resize(b,(ndq,ndp))

    [b,index] = readBinary(f, index, 8*ndq*ndp, '>d')
    hbz = resize(b,(ndq,ndp))

    [b,index] = readBinary(f, index, 8*ndq*ndp, '>d')
    hbt = resize(b,(ndq,ndp))

    [pssol,index] = readBinary(f, index, 8*ndp, '>d')
    [psprv,index] = readBinary(f, index, 8*ndp, '>d')
    [psman,index] = readBinary(f, index, 8*ndp, '>d')

    # integer : nndx, nndy
    [[nndx,nndy],index] = readBinary(f, index, 4*2, '>i')
    # print(nndx,nndy)

    # Defined in cplmet integer:
    #     >   nosol, noprv, nompl
    #     >  ,jcdp1, jcdp2, jcxp1, jcxp2, jcmax
    #     >  ,icwl1, icspx, icwl2, icmps, icmpe, icaxs
    #     >  ,itsls, itsle, itpvs, itpve, itmps, itmpe, itmax
    #     >  ,icmin(ndx), icmax(ndx), jtmin(ndy), jtmax(ndy)
    #     >  ,kgdx(ndx,ndy,4), kgdy(ndx,ndy,4), kreg(ndx,ndy)
    #     >  ,jcel(ndx,ndy), icel(ndx,ndy), jnxp(ndx,ndy), jnxm(ndx,ndy)
    #     >  ,kce, kcw, kcn, kcs
    # Defined in cplmet  real*8:
    #     >  ,vlmn(ndx), romn(ndy)
    #     >  ,hvol(ndx,ndy), hgdx(ndx,ndy), hgdy(ndx,ndy)
    #     >  ,hdxm(ndx,ndy), hdxp(ndx,ndy)
    #     >  ,hdsp(ndx,ndy,4), hdsv(ndx,ndy,4), hvsb(ndx,ndy)
    #     >  ,hare(ndx,ndy), hpit(ndx,ndy)
    #     >  ,hwtm(ndx,ndy), hwtp(ndx,ndy)
    #     >  ,gdsv(ndx,ndy,4), gare(ndx,ndy)
    #     >  ,gwtm(ndx,ndy), gwtp(ndx,ndy)

    [[nosol, noprv, nompl],index] = readBinary(f, index, 4*3, '>i')
    [[jcdp1, jcdp2, jcxp1, jcxp2, jcmax],index] = readBinary(f, index, 4*5, '>i')
    [[icwl1, icspx, icwl2, icmps, icmpe, icaxs],index] = readBinary(f, index, 4*6, '>i')
    [[itsls, itsle, itpvs, itpve, itmps, itmpe, itmax],index] = readBinary(f, index, 4*7, '>i')
    
    [icmin,index] = readBinary(f, index, 4*ndx, '>i')
    [icmax,index] = readBinary(f, index, 4*ndx, '>i')
    [jtmin,index] = readBinary(f, index, 4*ndy, '>i')
    [jtmax,index] = readBinary(f, index, 4*ndy, '>i')

    [b,index] = readBinary(f, index, 4*ndx*ndy*4, '>i')
    kgdx = resize(b,(ndx,ndy,4))
    
    [b,index] = readBinary(f, index, 4*ndx*ndy*4, '>i')
    kgdy = resize(b,(ndx,ndy,4))
    
    [b,index] = readBinary(f, index, 4*ndx*ndy, '>i')
    kreg = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 4*ndx*ndy, '>i')
    jcel = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 4*ndx*ndy, '>i')
    icel = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 4*ndx*ndy, '>i')
    jnxp = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 4*ndx*ndy, '>i')
    jnxm = resize(b,(ndx,ndy))

    [[kce, kcw, kcn, kcs],index] = readBinary(f, index, 4*4, '>i')

    [vlmn,index] = readBinary(f, index, 8*ndx, '>d')
    [romn,index] = readBinary(f, index, 8*ndy, '>d')

    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hvol = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hgdx = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hgdy = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hdxp = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hdxm = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 8*ndx*ndy*4, '>d')
    hdsp = resize(b,(ndx,ndy,4))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy*4, '>d')
    hdsv = resize(b,(ndx,ndy,4))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hvsb = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hare = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hpit = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hwtp = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    hwtm = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 8*ndx*ndy*4, '>d')
    gdsv = resize(b,(ndx,ndy,4))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    gare = resize(b,(ndx,ndy))

    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    gwtp = resize(b,(ndx,ndy))
    
    [b,index] = readBinary(f, index, 8*ndx*ndy, '>d')
    gwtm = resize(b,(ndx,ndy))
    
    # -----------------------------------------------
    # Places with a lot of data are grouped by record。
    cgdcom  = [mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx] #0~8
    cgdcom += [mpw1, mpsp, mpw2, mpax, npmx] #9~13
    cgdcom += [mpsol, mpprv, mpman] #14~16
    cgdcom += [grdx, grdy, hbr, hbz, hbt] #17~21
    cgdcom += [pssol, psprv, psman] #22~24

    cplmet  = [nosol, noprv, nompl] #0~2
    cplmet += [jcdp1, jcdp2, jcxp1, jcxp2, jcmax] #3~7
    cplmet += [icwl1, icspx, icwl2, icmps, icmpe, icaxs] #8~13
    cplmet += [itsls, itsle, itpvs, itpve, itmps, itmpe, itmax] #14~20
    cplmet += [icmin, icmax, jtmin, jtmax] #21~24
    cplmet += [kgdx, kgdy, kreg] #25~27
    cplmet += [jcel, icel, jnxp, jnxm] #28~31
    cplmet += [kce, kcw, kcn, kcs] #32~35
    cplmet += [vlmn, romn] #36~37
    cplmet += [hvol, hgdx, hgdy] #38~40
    cplmet += [hdxm, hdxp] #41~42
    cplmet += [hdsp, hdsv, hvsb] #43~45
    cplmet += [hare, hpit] #46~47
    cplmet += [hwtm, hwtp] #48~49
    cplmet += [gdsv, gare] #50~51
    cplmet += [gwtm, gwtp] #52~53
   
    return cver, ndq, ndp, cgdcom, nndx, nndy, cplmet
    
# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':

    from getData import getEnv
    # get information
    env = getEnv()

    # setting path
    fileMSH = env[0][1]
    
    # program discription
    print('readMSH.py')
    print('--------------------')
    print(' > function [readMSH] read MSH_***/dat.mtrc from datafile')
    print('   , using function [readBinary.py getBinary & readBinary & resize]')
    print('--------------------')
    print('[test] read this file ...')
    print(fileMSH)
    print()
    
    # read data test
    [ndx, ndy]=[160,80]
    [cver, nndq, nndp, cgdcom, nndx, nndy, cplmet]=readMSH(fileMSH, ndx, ndy)

    # print of read data
    print('cver', cver)
    print('nndq, nndp', nndq, nndp)

    # cgdcom
    [mqd1, mqx1, mqs1, mqs2, mqh1, mqh2, mqx2, mqd2, nqmx] = cgdcom[0:9]#0~8
    [mpw1, mpsp, mpw2, mpax, npmx] = cgdcom[9:14] #9~13
    [mpsol, mpprv, mpman] = cgdcom[14:17] #14~16
    [grdx, grdy, hbr, hbz, hbt] = cgdcom[17:22] #17~21
    [pssol, psprv, psman] = cgdcom[22:] #22~24

    # cplmet
    [nosol, noprv, nompl] = cplmet[0:3]#0~2
    [jcdp1, jcdp2, jcxp1, jcxp2, jcmax] = cplmet[3:8]#3~7
    [icwl1, icspx, icwl2, icmps, icmpe, icaxs] = cplmet[8:14]#8~13
    [itsls, itsle, itpvs, itpve, itmps, itmpe, itmax] = cplmet[14:21]#14~20
    [icmin, icmax, jtmin, jtmax] = cplmet[21:25]#21~24
    [kgdx, kgdy, kreg] = cplmet[25:28]#25~27
    [jcel, icel, jnxp, jnxm] = cplmet[28:32]#28~31
    [kce, kcw, kcn, kcs] = cplmet[32:36]#32~35
    [vlmn, romn] = cplmet[36:38]#36~37
    [hvol, hgdx, hgdy] = cplmet[38:41]#38~40
    [hdxm, hdxp] = cplmet[41:43]#41~42
    [hdsp, hdsv, hvsb] = cplmet[43:46]#43~45
    [hare, hpit] = cplmet[46:48]#46~47
    [hwtm, hwtp] = cplmet[48:50]#48~49
    [gdsv, gare] = cplmet[50:52]#50~51
    [gwtm, gwtp] = cplmet[52:]#52~53

    print('hgdx.shape, hgdy.shape', hgdx.shape, hgdy.shape) # 160,80 160,80
    print('======================================================================')

