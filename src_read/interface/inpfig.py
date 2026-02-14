# ======================================================================
# -*- name  : inpfig.py
# -*- coding: utf-8 -*-
# ======================================================================
# Rewrite namelist for Fortran as a class variable for Python
# Use the setting of inpfig_SA
# ------------------------------------------------------------------
# *Notation notes
# In Fortran, it is written as csiz (1), but in Python, it is obtained by csiz [1].
# Arrays are counted from 0 in Python
# In Python, you don't need "," to separate variables
#   When writing multiple variables on one line, add ";" at the end of each variable.
# In Python, np.arange (0.0, 2.0, 0.25) can be used to create an array in which "0.0 to 2.0"
#    is separated by "0.25".
# In Fortran, it is written as "-3.0d4", but in Python, double-precision
#    real numbers are written as e, so it is written as "-3.0e4".
# ======================================================================

# import
import sys
import numpy as np

class upfgsz:
    csiz = []
    csiz.append("Main  :  0.5,   5.0,   -3.2,   3.2") # ( 1) => csize[0]
    csiz.append("Div   :  1.6,   3.2,   -2.9,  -1.9") # ( 2)
    csiz.append("MTop :   1.6,   4.6,    0.0,   3.0") # ( 3)
    csiz.append("I-div : 1.85,  2.15,   -2.5,  -2.2") # ( 4)
    csiz.append("O-div : 2.35,  2.65,  -2.85, -2.55") # ( 5)
   
    cleg = []
    cleg.append("5.6,  1.0") # ( 1)
    cleg.append("2.8, -2.4") # ( 2)
    cleg.append("3.6,  2.5") # ( 3)
    cleg.append("2.2, -2.3") # ( 4)
    cleg.append("2.7, -2.6") # ( 4)
   
class upfgmx:
    fpnemx_sol =  99.0
    fpnemx_prv =  99.0
    fpnemx_man =  99.0
   
    fptemx_sol =  99.0
    fptemx_prv =  99.0
    fptemx_man =  99.0
   
    fpvpmx_sol =   5.0
    fpvpmx_prv =   5.0
    fpvpmx_man =   1.0
   
    fvnemx_sol =  12.0
    fvnemx_prv = 150.0
    fvnemx_man =  99.0
   
    fvtemx_sol = 4.0e3
    fvtemx_prv = 120.0
    fvtemx_man =  99.0
   
    fvvpmx_sol =   5.0
    fvvpmx_prv =   5.0
    fvvpmx_man =   1.0
   
    fhnemx_sol =  99.0
    fhnemx_prv =  99.0
    fhnemx_man =   5.0
   
    fhtemx_sol =  99.0
    fhtemx_prv =   5.0
    fhtemx_man =  99.0
   
    fhvpmx_sol =   1.0
    fhvpmx_prv =  99.0
    fhvpmx_man =   4.0
   
    fhflmn = -0.3e22
    fhflmx = 1.5e22
   
    fhitmn = 0.0
    fhitmx = 99.0

class upcntv:
    Ni     = [0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0]
    log_Ni = np.arange(0.0, 2.0, 0.25)
    Ti     = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
    log_Ti = np.arange(0.0, 3.2, 0.4)
    Te     = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
    log_Te = np.arange(0.0, 3.2, 0.4)
    Vf     = np.arange(-3.0e4, 3.0e4, 0.5e4)

    vni     = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0])*1e19
    log_vni = np.arange(0.0, 2.0, 0.25)*19
    vne     = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 10.0, 20.0, 40.0])*1e19
    vti     = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
    log_vti = np.arange(0.0, 3.2, 0.4)
    vte     = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0]
    log_vte = np.arange(0.0, 3.2, 0.4)
    vva     = np.arange(-3.0e4, 3.0e4, 0.5e4)

    vwrt     = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    log_vwrt = np.arange(0.0, 1.6, 0.2)
    vwrm     = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    log_vwrm = np.arange(0.0, 1.6, 0.2)
    vwrc     = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    log_vwrc = np.arange(0.0, 1.6, 0.2)

    N0     = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
    log_N0 = np.arange(-3.0, 2.0, 1.0)
    E0     = [2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1.0e3]
    log_E0 = np.arange(0.0, 3.2, 0.4)
    Ng     = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
    log_Ng = np.arange(-3.0, 2.0, 1.0)
    Eg     = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    log_Eg = np.arange(0.0, 3.2, 0.4)

    Wrd     = np.arange(2.0, 20.0, 2.0)
    log_Wrd = np.arange(0.0, 1.6, 0.2)
    Nz0     = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
    log_Nz0 = np.arange(-3.0, 2.0, 1.0)
    Nzi     = [0.001, 0.01, 0.1, 1.0, 10.0, 100.0]
    log_Nzi = np.arange(-3.0, 2.0, 1.0)

class uptub:
    tbtub = []
    tbtub.append("30")  # ( 1)
    tbtub.append("29")  # ( 2)
    tbtub.append("28")  # ( 3)
    tbtub.append("25")  # ( 4)
    tbtub.append("20")  # ( 5)
    tbtub.append("15")  # ( 6)
    tbtub.append("10")  # ( 7)
    tbtub.append(" 5")  # ( 8)
    tbtub.append("38")  # ( 9)
    tbtub.append("39")  # (10)
    tbtub.append("40")  # (11)
    tbtub.append("  ")  # (12)

# ======================================================================
# running program
# ======================================================================
if __name__ == "__main__":
    # csiz
    print('calss csiz:')
    print('------------------------------')
    csiz = upfgsz.csiz
    for i in range(0, len(csiz)):
        print(i, csiz[i])
    print('------------------------------')
    print()
   
    # upcntv
    print('calss upcntv:')
    print('------------------------------')
    for key, value in upcntv.__dict__.items():
        print(key, ':', value)
    print('------------------------------')
    print()
