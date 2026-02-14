# ======================================================================
# -*- name  : readMSH2.py
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
def readMSH2(fileName):
    with open(fileName, "rb") as f:
        [data, header, footer] = getBinary(f, len(f.read()))

    index = 0

    f = data

    # -----------------------------------------------
    # character : cver
    [cverS,index] = readBinary(f, index, header[0], '>c')
    cver=''
    for s in cverS:
        cver+=str(s.decode())

    [[ndmg, ndmc, ndms, ndx, ndy, nvxd, nvyd],index] = readBinary(f, index, header[1], '>i')

    [xpnt,index] = readBinary(f, index, 8*ndmg, '>d')
    [ypnt,index] = readBinary(f, index, 8*ndmg, '>d')

    [volm,index] = readBinary(f, index, 8*ndmc, '>d')
    
    [bvx,index] = readBinary(f, index, 8*(1+ndmc), '>d')
    [bvy,index] = readBinary(f, index, 8*(1+ndmc), '>d')
    [bvz,index] = readBinary(f, index, 8*(1+ndmc), '>d')

    [[grax, gzax, grxp, gzxp],index] = readBinary(f, index, 8*4, '>d')

    [mclx,index] = readBinary(f, index, 4*ndx, '>i')
    [mcly,index] = readBinary(f, index, 4*ndy, '>i')
    [mcly2,index] = readBinary(f, index, 4*ndy, '>i')
    [mxjw,index] = readBinary(f, index, 4*ndx, '>i')
    [mxjw2,index] = readBinary(f, index, 4*ndx, '>i')
    [myit,index] = readBinary(f, index, 4*ndy, '>i')

    [b,index] = readBinary(f, index, 4*(1+ndx)*(1+ndy), '>i')
    nogd = resize(b, (1+ndx,1+ndy))
    [b,index] = readBinary(f, index, 4*(1+ndx)*(1+ndy), '>i')
    nocl = resize(b, (1+ndx,1+ndy))
    [b,index] = readBinary(f, index, 4*(1+ndx)*(1+ndy), '>i')
    mcel = resize(b, (1+ndx,1+ndy))
    [b,index] = readBinary(f, index, 4*(1+nvxd)*(1+nvyd), '>i')
    nogdv = resize(b, (1+nvxd,1+nvyd))
    [b,index] = readBinary(f, index, 4*(1+nvxd)*(1+nvyd), '>i')
    noclv = resize(b, (1+nvxd,1+nvyd))

    [mseg,index] = readBinary(f, index, 4*ndmc, '>i')
    [b,index] = readBinary(f, index, 4*ndmc*(ndms+1), '>i')
    mgrd = resize(b, (ndmc,(ndms+1)))
    [b,index] = readBinary(f, index, 4*ndmc*ndms, '>i')
    mknd = resize(b, (ndmc,ndms))

    [mrgn,index] = readBinary(f, index, 4*ndmc, '>i')
    [migx,index] = readBinary(f, index, 4*ndmc, '>i')
    [migy,index] = readBinary(f, index, 4*ndmc, '>i')
    [b,index] = readBinary(f, index, 4*ndmc*ndms, '>i')
    Next = resize(b, (ndmc,ndms))
    [iplx,index] = readBinary(f, index, 4*(1+ndmc), '>i')
    [iply,index] = readBinary(f, index, 4*(1+ndmc), '>i')

    [[ngmax, ngmax2, ncmax, ncmax2],index] = readBinary(f, index, 4*4, '>i')

    [[jdp1, jxp1, jsl1, jsl2, jxp2, jdp2],index] = readBinary(f, index, 4*6, '>i')

    [[iwl1, ispx, iwl2, impl, iaxs, ivs1, ivs2],index] = readBinary(f, index, 4*7, '>i')

    [[jvb1,jvb2, ivb1, ivb2],index] = readBinary(f, index, 4*4, '>i')

    [[ndwl, ndwh, ndwp],index] = readBinary(f, index, header[3], '>i')

    [rcywl,index] = readBinary(f, index, 8*ndwh, '>d')
    [temwl,index] = readBinary(f, index, 8*ndwh, '>d')

    [rcwl,index] = readBinary(f, index, 8*ndwp, '>d')
    [e0wl,index] = readBinary(f, index, 8*ndwp, '>d')
    [cswl,index] = readBinary(f, index, 8*ndwp, '>d')
    [snwl,index] = readBinary(f, index, 8*ndwp, '>d')

    [cwalS,index] = readBinary(f, index, ndwl*4, '>c')
    cwal=''
    for s in cwalS:
        cwal+=str(s.decode())

    [chwlS,index] = readBinary(f, index, ndwh*4, '>c')
    chwl=''
    for s in chwlS:
        chwl+=str(s.decode())

    [tywlS,index] = readBinary(f, index, ndwp*4, '>c')
    tywl=''
    for s in tywlS:
        tywl+=str(s.decode())

    [npsw,index] = readBinary(f, index, 4*ndwl, '>i')
    [npew,index] = readBinary(f, index, 4*ndwl, '>i')

    [inwl,index] = readBinary(f, index, 4*ndwp, '>i')
    [ikwl,index] = readBinary(f, index, 4*ndwp, '>i')
    [ipwl,index] = readBinary(f, index, 4*ndwp, '>i')
    [icwl,index] = readBinary(f, index, 4*ndwp, '>i')
    [isxw,index] = readBinary(f, index, 4*ndwp, '>i')

    [isyw,index] = readBinary(f, index, 4*ndwp, '>i')
    [ixyw,index] = readBinary(f, index, 4*ndwp, '>i')
    [igtw,index] = readBinary(f, index, 4*ndwp, '>i')
    [ihwl,index] = readBinary(f, index, 4*ndwp, '>i')
    [ipmp,index] = readBinary(f, index, 4*ndwp, '>i')

    [[nowl, npwl],index] = readBinary(f, index, 4*2, '>i')

    [[npmp, nowl2, npwl2, nhwl],index] = readBinary(f, index, 4*4, '>i')

    [[ndgt, ndgtp, ndad, ndpm],index] = readBinary(f, index, header[5], '>i')

    [b,index] = readBinary(f, index, 8*ndgt*2, '>d')
    xsgt = resize(b, (ndgt,2))
    [b,index] = readBinary(f, index, 8*ndgt*2, '>d')
    xegt = resize(b, (ndgt,2))
    [b,index] = readBinary(f, index, 8*ndgt*2, '>d')
    ysgt = resize(b, (ndgt,2))
    [b,index] = readBinary(f, index, 8*ndgt*2, '>d')
    yegt = resize(b, (ndgt,2))

    [b,index] = readBinary(f, index, 4*ndgt*2, '>i')
    nsgt = resize(b, (ndgt,2))
    [b,index] = readBinary(f, index, 4*ndgt*2, '>i')
    negt = resize(b, (ndgt,2))

    [iwgt,index] = readBinary(f, index, 4*ndgtp, '>i')
    [ipgt,index] = readBinary(f, index, 4*ndgtp, '>i')
    [icgt,index] = readBinary(f, index, 4*ndgtp, '>i')

    [[nogt],index] = readBinary(f, index, 4, '>i')
    
    cntgrd  = [xpnt,ypnt,volm,bvx,bvy,bvz,grax,gzax,grxp,gzxp] # 1~10
    cntgrd += [mclx,mcly,mcly2,mxjw,mxjw2,myit] # 11~16
    cntgrd += [nogd,nocl,mcel,nogdv,noclv,mseg,mgrd,mknd] # 17~24
    cntgrd += [mrgn,migx,migy,Next,iplx,iply] # 25~30
    cntgrd += [ngmax,ngmax2,ncmax,ncmax2] # 31~34
    cntgrd += [jdp1,jxp1,jsl1,jsl2,jxp2,jdp2] # 35~40
    cntgrd += [iwl1,ispx,iwl2,impl,iaxs,ivs1,ivs2] # 41~47
    cntgrd += [jvb1, jvb2, ivb1, ivb2] # 48~51

    cntwal  = [rcywl, temwl, rcwl, e0wl, cswl, snwl] # 1~6
    cntwal += [cwal, chwl, tywl] # 7~9
    cntwal += [npsw, npew] # 10~11
    cntwal += [inwl, ikwl, ipwl] # 12~14
    cntwal += [icwl, isxw, isyw, ixyw, igtw] # 15~19
    cntwal += [ihwl, ipmp] # 20~21
    cntwal += [nowl, npwl, npmp, nowl2, npwl2, nhwl] # 22~27

    cntgat  = [xsgt,xegt,ysgt,yegt] # 1~4
    cntgat += [nsgt,negt,iwgt,ipgt,icgt,nogt] # 5~10

    com_size  = [ndmg, ndmc, ndms, ndx, ndy, nvxd, nvyd] # 1~7
    com_size += [ndwl, ndwh, ndwp] # 8~10
    com_size += [ndgt, ndgtp, ndad, ndpm] # 11~14

    return cntgrd, cntwal, cntgat, com_size

# ======================================================================
# running program
# ======================================================================
if __name__ == '__main__':

    from getData import getEnv
    # get information
    env =  getEnv()

    fileMSH = env[0][8]

    [cntgrd, cntwal, cntgat, com_size]=readMSH2(fileMSH)