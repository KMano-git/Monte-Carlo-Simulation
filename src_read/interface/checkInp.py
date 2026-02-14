# ======================================================================
# -*- name  : checkInp.py
# -*- coding: utf-8 -*-
# ======================================================================

# import
import os
import sys

#-------------------------------------------------------------------------------
# check input and file search
# 1st argument:login user name
# 2nd argument:DB name(device)
# 3rd argument:IMAS major version No.
# 4th argument:run No.(0 to 99999)
# 5th argument:shot No.(0 to 214748)
#-------------------------------------------------------------------------------
def checkInp(user, device, version, run, shot):
    # check directry
    dirname = "/home/" + user + "/public/imasdb/" + device
    #dirname = "../../" + user + "/public/imasdb/" + device # debug
    if os.path.exists(dirname) == False:
        print("ERROR : device \"" + device + "\" is an invalid input!")
        sys.exit()

    dirname = dirname + "/" + version
    if os.path.exists(dirname) == False:
        print("ERROR : version \"" + version + "\" is an invalid input!")
        sys.exit()

    # directry name
    dirname = dirname + "/" + str(int(run / 10000))
    
    # file name
    r = run % 10000

    if shot == 0 :
        if r < 1000:
            file = "ids_" + str("{:>03d}".format(r))
        else:
            file = "ids_" + str("{:>04d}".format(r))
    else:
        file = "ids_" + str(shot) + str("{:>04d}".format(r))
    flist = []
    flist.append(file + ".characteristics")
    flist.append(file + ".datafile")
    flist.append(file + ".tree")

    flg = 0
    # check exist
    for file in flist:
        fname = dirname + "/" + file
        if os.path.exists(fname):
            # return find flag(found)
            flg = 1
            break

    # return find flag(not found)
    return flg

# ======================================================================
# running the program
# ======================================================================
if __name__ == "__main__":
    # import 
    from getData import getEnv
    from getName import getName

    print("start program \"checkInp.py\"\n")
    user = getEnv()[1]
    data = getName(user)
    dev = input("device >>>")
    ver = input("version >>>")
    run = int(input("shot[0 to 214748] >>>"))
    shot = int(input("run[0 to 99999] >>>"))
    flg = checkInp(user, dev, ver, run, shot)
    # flg=0, not found file
    # flg=1, found file
    print(flg)
    print("\nend program \"checkInp.py\"")
