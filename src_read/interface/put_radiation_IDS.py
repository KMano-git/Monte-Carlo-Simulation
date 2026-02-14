# impor imas mod

try:
    import imas
except Exception as e:
    print("Required IMAS support library not available on this system.")

# ======================================================================
# store radiation data for IDS
# 1st argument : IDS database
# 2nd argument : radiation datas
# 3rd argument : data size[X, Y]
# 4th argument : element data
# ======================================================================
def rd2IDS(imas_obj, data, ds, elm):
    [wswe, twrd1, twrd2] = data
    dsize = ds[0] * ds[1]
    radname = ["wswe"]
    for i in range(len(elm)):
        radname.append("twrd_" + elm[i][0])
    #--------------------------------------
    # radiation
    imas_obj.radiation.grid_ggd.resize(1)
    imas_obj.radiation.time.resize(1)
    imas_obj.radiation.time[0] = 1
    imas_obj.radiation.ids_properties.homogeneous_time = 1
    imas_obj.radiation.process.resize(1)
    imas_obj.radiation.process[0].ggd.resize(1)
    imas_obj.radiation.process[0].ggd[0].neutral.resize(1) 
    imas_obj.radiation.process[0].ggd[0].ion.resize(2) 

    # neutral emissivity(wswe)
    ntlse = imas_obj.radiation.process[0].ggd[0].neutral[0]
    ntlse.state.resize(1)
    ntlse.state[0].label = radname[0]
    ntlse.state[0].emissivity.resize(1)
    ntlse = ntlse.state[0].emissivity[0]
    ntlse.values.resize(dsize)

    # ion emissivity(twrd1)
    inse1 = imas_obj.radiation.process[0].ggd[0].ion[0]
    inse1.state.resize(1)
    inse1.state[0].label = radname[1]
    inse1.state[0].emissivity.resize(1)
    inse1 = inse1.state[0].emissivity[0]
    inse1.values.resize(dsize)

    # ion emissivity(twrd2)
    inse2 = imas_obj.radiation.process[0].ggd[0].ion[1]
    inse2.state.resize(1)
    inse2.state[0].label = radname[2]
    inse2.state[0].emissivity.resize(1)
    inse2 = inse2.state[0].emissivity[0]
    inse2.values.resize(dsize)

    # data insert
    for i in range(dsize):
        # neutral
        ntlse.values[i] = wswe[i]

        # ion
        inse1.values[i] = twrd1[i]
        inse2.values[i] = twrd2[i]

    return imas_obj
