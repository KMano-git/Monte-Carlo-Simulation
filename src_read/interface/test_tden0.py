from readDTNTL import readDTNTL
from readMSH2 import readMSH2
from getData  import getEnv

def tden0():
    flist = getEnv()[0]
    fileDTNTL = flist[3]
    fileMSH2 = flist[8]
    data = readDTNTL(fileDTNTL)
    vmwork = data[2]
    visrc = data[8]
    [cntgrd, cntwal] = readMSH2(fileMSH2)[:2]
    npwl2 = cntwal[24]
    volm = cntgrd[2]
    Mrgn = cntgrd[24]
    ncmax = cntgrd[30]
    ncmax2 = cntgrd[31]

    ngas = 1
    mfmax = 8
    nwkmp_nt = len(vmwork)

    tden0 = []

    wmwork = []
    for m in range(mfmax):
        den0 = []
        tden0.append([])
        isrc = visrc[m]
        if 0<isrc & isrc<7:
            for i in range(nwkmp_nt):
                if isrc > 0:
                    wmwork.append(vmwork[i][isrc-1])
            wwal = []
            for j in range(480):
                try:
                    wwal.append(vmwork[84994 + j][isrc-1])
                except:
                    pass
            wssn = []
            wden = []
            dl = 6500 + 1
            for j in range(dl):
                try:
                    wssn.append(vmwork[1 + 0 * dl + j][isrc-1])
                    wden.append(vmwork[1 + 5 * dl + j][isrc-1])
                except:
                    pass

            swabs = 0
            for iw in range(len(wwal)):
                swabs += wwal[iw]

            swreg = [0,0,0,0,0,0,0,0,0,0,0]
            for ic in range(ncmax2):
                ir = Mrgn[ic]
                if ir < 6 & 0 < ir:
                    swreg[ir] = swreg[ir] + wssn[ic]

            swion = swreg[1] + swreg[2] + swreg[3] + swreg[4] + swreg[5] + swreg[6] + swreg[7]
            
            swnrm = swion + swabs
            for i in range(ncmax2):
                den0.append(wden[i]*wmwork[0]/swnrm/volm[i])
        tden0[m] = den0
    return tden0

if __name__ == "__main__":
    a = tden0()
    for k in range(len(a)):
        with open("tden0_"+str(k)+".csv","w") as f:
            for i in range(len(a[k])):
                f.write(str(a[k][i]) + ",")
   