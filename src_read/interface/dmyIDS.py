# impor imas mod
try:
    import imas
except Exception as e:
    print("Required IMAS support library not available on this system.")

# import 
from makeDumy import makeDumy
# ======================================================================
# store dumy data for IDS
# 1st argument : IDS database
# 2nd argument : data size[X, Y]
# 3rd argument : grid data[X, Y]
# ======================================================================
def dmyIDS(imas_obj, ds, grd):
    # grid data
    xc = grd[0]
    yc = grd[1]
    nx = ds[0]
    ny = ds[1]
    dsize = nx*ny

    # 
    shot = 1
    run = 1
    user = "nais5001_713"
    device = "iter"
    version = "3"
    dirpath = "."

    # make dumy data
    dmy = makeDumy(dsize, 0, 0)

    #--------------------------------------
    # edge_profiles(electrons)
    pePath = imas_obj.edge_profiles.ggd[0].electrons

    # temperature
    pePath.temperature[0].coefficients.resize([dsize, 2])

    # density
    pePath.density[0].coefficients.resize([dsize, 2])

    # density_fast
    pePath.density_fast.resize(1)
    pePath.density_fast[0].grid_subset_index = 1
    pePath.density_fast[0].values.resize(dsize)
    pePath.density_fast[0].coefficients.resize([dsize, 2])

    # pressure
    pePath.pressure.resize(1)
    pePath.pressure[0].grid_subset_index = 1
    pePath.pressure[0].values.resize(dsize)
    pePath.pressure[0].coefficients.resize([dsize, 2])

    # pressure_fast_perpendicular
    pePath.pressure_fast_perpendicular.resize(1)
    pePath.pressure_fast_perpendicular[0].grid_subset_index = 1
    pePath.pressure_fast_perpendicular[0].values.resize(dsize)
    pePath.pressure_fast_perpendicular[0].coefficients.resize([dsize, 2])

    # pressure_fast_parallel
    pePath.pressure_fast_parallel.resize(1)
    pePath.pressure_fast_parallel[0].grid_subset_index = 1
    pePath.pressure_fast_parallel[0].values.resize(dsize) 
    pePath.pressure_fast_parallel[0].coefficients.resize([dsize, 2])

    # velocity
    pePath.velocity[0].radial.resize(dsize)
    pePath.velocity[0].radial_coefficients.resize([dsize, 2])
    pePath.velocity[0].diamagnetic.resize(dsize)
    pePath.velocity[0].diamagnetic_coefficients.resize([dsize, 2])
    pePath.velocity[0].parallel_coefficients.resize([dsize, 2])
    pePath.velocity[0].poloidal.resize(dsize)
    pePath.velocity[0].poloidal_coefficients.resize([dsize, 2])
    pePath.velocity[0].toroidal.resize(dsize)

    #--------------------------------------
    # edge_profiles(ion)
    piPath = imas_obj.edge_profiles.ggd[0].ion[0]

    # temperature
    piPath.temperature[0].coefficients.resize([dsize, 2])

    # density
    piPath.density[0].coefficients.resize([dsize, 2])

    # density_fast
    piPath.density_fast.resize(1)
    piPath.density_fast[0].grid_subset_index = 1
    piPath.density_fast[0].values.resize(dsize)
    piPath.density_fast[0].coefficients.resize([dsize, 2])

    # pressure
    piPath.pressure.resize(1)
    piPath.pressure[0].grid_subset_index = 1
    piPath.pressure[0].values.resize(dsize)
    piPath.pressure[0].coefficients.resize([dsize, 2])

    # pressure_fast_perpendicular
    piPath.pressure_fast_perpendicular.resize(1)
    piPath.pressure_fast_perpendicular[0].grid_subset_index = 1
    piPath.pressure_fast_perpendicular[0].values.resize(dsize)
    piPath.pressure_fast_perpendicular[0].coefficients.resize([dsize, 2])

    # pressure_fast_parallel
    piPath.pressure_fast_parallel.resize(1)
    piPath.pressure_fast_parallel[0].grid_subset_index = 1
    piPath.pressure_fast_parallel[0].values.resize(dsize) 
    piPath.pressure_fast_parallel[0].coefficients.resize([dsize, 2])

    # velocity
    piPath.velocity[0].radial.resize(dsize)
    piPath.velocity[0].radial_coefficients.resize([dsize, 2])
    piPath.velocity[0].diamagnetic.resize(dsize)
    piPath.velocity[0].diamagnetic_coefficients.resize([dsize, 2])
    piPath.velocity[0].parallel_coefficients.resize([dsize, 2])
    piPath.velocity[0].poloidal.resize(dsize)
    piPath.velocity[0].poloidal_coefficients.resize([dsize, 2])
    piPath.velocity[0].toroidal.resize(dsize)

    # energy_density_kinetic
    piPath.energy_density_kinetic.resize(1)
    piPath.energy_density_kinetic[0].grid_subset_index = 1
    piPath.energy_density_kinetic[0].values.resize(dsize) 
    piPath.energy_density_kinetic[0].coefficients.resize([dsize, 2])

    # multiple_states_flag
    piPath.multiple_states_flag = 0

    #--------------------------------------
    # edge_profiles(neutral)
    pnPath = imas_obj.edge_profiles.ggd[0].neutral[0]

    # label
    pnPath.label = "NEUTRAL"

    # temperature
    pnPath.temperature.resize(1)
    pnPath.temperature[0].grid_subset_index = 1
    pnPath.temperature[0].values.resize(dsize)
    pnPath.temperature[0].coefficients.resize([dsize, 2])

    # density
    pnPath.density.resize(1)
    pnPath.density[0].grid_subset_index = 1
    pnPath.density[0].values.resize(dsize)
    pnPath.density[0].coefficients.resize([dsize, 2])

    # density_fast
    pnPath.density_fast.resize(1)
    pnPath.density_fast[0].grid_subset_index = 1
    pnPath.density_fast[0].values.resize(dsize)
    pnPath.density_fast[0].coefficients.resize([dsize, 2])

    # pressure
    pnPath.pressure.resize(1)
    pnPath.pressure[0].grid_subset_index = 1
    pnPath.pressure[0].values.resize(dsize)
    pnPath.pressure[0].coefficients.resize([dsize, 2])

    # pressure_fast_perpendicular
    pnPath.pressure_fast_perpendicular.resize(1)
    pnPath.pressure_fast_perpendicular[0].grid_subset_index = 1
    pnPath.pressure_fast_perpendicular[0].values.resize(dsize)
    pnPath.pressure_fast_perpendicular[0].coefficients.resize([dsize, 2])

    # pressure_fast_parallel
    pnPath.pressure_fast_parallel.resize(1)
    pnPath.pressure_fast_parallel[0].grid_subset_index = 1
    pnPath.pressure_fast_parallel[0].values.resize(dsize) 
    pnPath.pressure_fast_parallel[0].coefficients.resize([dsize, 2])

    # velocity
    pnPath.velocity.resize(1)
    pnPath.velocity[0].radial.resize(dsize)
    pnPath.velocity[0].radial_coefficients.resize([dsize, 2])
    pnPath.velocity[0].diamagnetic.resize(dsize)
    pnPath.velocity[0].diamagnetic_coefficients.resize([dsize, 2])
    pnPath.velocity[0].parallel.resize(dsize)
    pnPath.velocity[0].parallel_coefficients.resize([dsize, 2])
    pnPath.velocity[0].poloidal.resize(dsize)
    pnPath.velocity[0].poloidal_coefficients.resize([dsize, 2])
    pnPath.velocity[0].toroidal.resize(dsize)

    # energy_density_kinetic
    pnPath.energy_density_kinetic.resize(1)
    pnPath.energy_density_kinetic[0].grid_subset_index = 1
    pnPath.energy_density_kinetic[0].values.resize(dsize) 
    pnPath.energy_density_kinetic[0].coefficients.resize([dsize, 2])

    # multiple_states_flag
    pnPath.multiple_states_flag = 0

    #--------------------------------------
    # edge_sources
    
    # edge_sources(electrons)
    sePath = imas_obj.edge_sources.source[0].ggd[0].electrons

    # particles
    sePath.particles.resize(1)
    sePath.particles[0].values.resize(dsize)
    sePath.particles[0].grid_subset_index = 1
    sePath.particles[0].coefficients.resize([dsize, 2])

    # energy
    sePath.energy.resize(1)
    sePath.energy[0].grid_subset_index = 1
    sePath.energy[0].values.resize(dsize)
    sePath.energy[0].coefficients.resize([dsize, 2])

    #--------------------------------------
    # edge_transport
    imas_obj.edge_transport.grid_ggd.resize(1)
    imas_obj.edge_transport.time.resize(1)
    imas_obj.edge_transport.time[0] = 1
    imas_obj.edge_transport.ids_properties.homogeneous_time = 1

    # grid_ggd

    # edge_transport(electrons)
    imas_obj.edge_transport.model.resize(1)
    imas_obj.edge_transport.model[0].ggd.resize(1)
    tePath = imas_obj.edge_transport.model[0].ggd[0].electrons

    # particles
    tePath.particles.d.resize(1)
    tePath.particles.d[0].values.resize(dsize)
    tePath.particles.d[0].coefficients.resize([dsize, 2])
    tePath.particles.v.resize(1)
    tePath.particles.v[0].values.resize(dsize)
    tePath.particles.v[0].coefficients.resize([dsize, 2])
    tePath.particles.flux.resize(1)
    tePath.particles.flux[0].values.resize(dsize)
    tePath.particles.flux[0].coefficients.resize([dsize, 2])
    tePath.particles.flux_limiter.resize(1)
    tePath.particles.flux_limiter[0].values.resize(dsize)
    tePath.particles.flux_limiter[0].coefficients.resize([dsize, 2])

    # energy
    tePath.energy.d.resize(1)
    tePath.energy.d[0].values.resize(dsize)
    tePath.energy.d[0].coefficients.resize([dsize, 2])
    tePath.energy.v.resize(1)
    tePath.energy.v[0].values.resize(dsize)
    tePath.energy.v[0].coefficients.resize([dsize, 2])
    tePath.energy.flux.resize(1)
    tePath.energy.flux[0].values.resize(dsize)
    tePath.energy.flux[0].coefficients.resize([dsize, 2])
    tePath.energy.flux_limiter.resize(1)
    tePath.energy.flux_limiter[0].values.resize(dsize)
    tePath.energy.flux_limiter[0].coefficients.resize([dsize, 2])

    #--------------------------------------
    # edge_transport(ion)
    imas_obj.edge_transport.model[0].ggd[0].ion.resize(1)
    tiPath = imas_obj.edge_transport.model[0].ggd[0].ion[0]

    tiPath.particles.d.resize(1)
    tiPath.particles.d[0].values.resize(dsize)
    tiPath.particles.d[0].coefficients.resize([dsize, 2])
    tiPath.particles.v.resize(1)
    tiPath.particles.v[0].values.resize(dsize)
    tiPath.particles.v[0].coefficients.resize([dsize, 2])
    tiPath.particles.flux.resize(1)
    tiPath.particles.flux[0].values.resize(dsize)
    tiPath.particles.flux[0].coefficients.resize([dsize, 2])
    tiPath.particles.flux_limiter.resize(1)
    tiPath.particles.flux_limiter[0].values.resize(dsize)
    tiPath.particles.flux_limiter[0].coefficients.resize([dsize, 2])

    # energy
    tiPath.energy.d.resize(1)
    tiPath.energy.d[0].values.resize(dsize)
    tiPath.energy.d[0].coefficients.resize([dsize, 2])
    tiPath.energy.v.resize(1)
    tiPath.energy.v[0].values.resize(dsize)
    tiPath.energy.v[0].coefficients.resize([dsize, 2])
    tiPath.energy.flux.resize(1)
    tiPath.energy.flux[0].values.resize(dsize)
    tiPath.energy.flux[0].coefficients.resize([dsize, 2])
    tiPath.energy.flux_limiter.resize(1)
    tiPath.energy.flux_limiter[0].values.resize(dsize)
    tiPath.energy.flux_limiter[0].coefficients.resize([dsize, 2])

    # insert
    for m in range(dsize):
        #--------------------------------------
        # edge_profiles
        # electrons
        pePath.temperature[0].coefficients[m][0] = dmy[m] * (6.242e18)
        pePath.temperature[0].coefficients[m][1] = dmy[m]
        pePath.density[0].coefficients[m][0] = dmy[m]
        pePath.density[0].coefficients[m][1] = dmy[m]
        pePath.density_fast[0].values[m] = dmy[m]
        pePath.density_fast[0].coefficients[m][0] = dmy[m]
        pePath.density_fast[0].coefficients[m][1] = dmy[m]
        pePath.pressure[0].values[m] = dmy[m]
        pePath.pressure[0].coefficients[m][0] = dmy[m]
        pePath.pressure[0].coefficients[m][1] = dmy[m]
        pePath.pressure_fast_perpendicular[0].values[m] = dmy[m]
        pePath.pressure_fast_perpendicular[0].coefficients[m][0] = dmy[m]
        pePath.pressure_fast_perpendicular[0].coefficients[m][1] = dmy[m]
        pePath.pressure_fast_parallel[0].values[m] = dmy[m]
        pePath.pressure_fast_parallel[0].coefficients[m][0] = dmy[m]
        pePath.pressure_fast_parallel[0].coefficients[m][1] = dmy[m]
        pePath.velocity[0].radial[m] = dmy[m]
        pePath.velocity[0].radial_coefficients[m][0] = dmy[m]
        pePath.velocity[0].radial_coefficients[m][1] = dmy[m]
        pePath.velocity[0].diamagnetic[m] = dmy[m]
        pePath.velocity[0].diamagnetic_coefficients[m][0] = dmy[m]
        pePath.velocity[0].diamagnetic_coefficients[m][1] = dmy[m]
        pePath.velocity[0].parallel_coefficients[m][0] = dmy[m]
        pePath.velocity[0].parallel_coefficients[m][1] = dmy[m]
        pePath.velocity[0].poloidal[m] = dmy[m]
        pePath.velocity[0].poloidal_coefficients[m][0] = dmy[m]
        pePath.velocity[0].poloidal_coefficients[m][1] = dmy[m]
        pePath.velocity[0].toroidal[m] = dmy[m]

        # ion
        piPath.temperature[0].coefficients[m][0] = dmy[m] * (6.242e18)
        piPath.temperature[0].coefficients[m][1] = dmy[m]
        piPath.density[0].coefficients[m][0] = dmy[m]
        piPath.density[0].coefficients[m][1] = dmy[m]
        piPath.density_fast[0].values[m] = dmy[m]
        piPath.density_fast[0].coefficients[m][0] = dmy[m]
        piPath.density_fast[0].coefficients[m][1] = dmy[m]
        piPath.pressure[0].values[m] = dmy[m]
        piPath.pressure[0].coefficients[m][0] = dmy[m]
        piPath.pressure[0].coefficients[m][1] = dmy[m]
        piPath.pressure_fast_perpendicular[0].values[m] = dmy[m]
        piPath.pressure_fast_perpendicular[0].coefficients[m][0] = dmy[m]
        piPath.pressure_fast_perpendicular[0].coefficients[m][1] = dmy[m]
        piPath.pressure_fast_parallel[0].values[m] = dmy[m]
        piPath.pressure_fast_parallel[0].coefficients[m][0] = dmy[m]
        piPath.pressure_fast_parallel[0].coefficients[m][1] = dmy[m]
        piPath.velocity[0].radial[m] = dmy[m]
        piPath.velocity[0].radial_coefficients[m][0] = dmy[m]
        piPath.velocity[0].radial_coefficients[m][1] = dmy[m]
        piPath.velocity[0].diamagnetic[m] = dmy[m]
        piPath.velocity[0].diamagnetic_coefficients[m][0] = dmy[m]
        piPath.velocity[0].diamagnetic_coefficients[m][1] = dmy[m]
        piPath.velocity[0].parallel_coefficients[m][0] = dmy[m]
        piPath.velocity[0].parallel_coefficients[m][1] = dmy[m]
        piPath.velocity[0].poloidal[m] = dmy[m]
        piPath.velocity[0].poloidal_coefficients[m][0] = dmy[m]
        piPath.velocity[0].poloidal_coefficients[m][1] = dmy[m]
        piPath.velocity[0].toroidal[m] = dmy[m]
        piPath.energy_density_kinetic[0].values[m] = dmy[m]
        piPath.energy_density_kinetic[0].coefficients[m][0] = dmy[m]
        piPath.energy_density_kinetic[0].coefficients[m][1] = dmy[m]

        # neutral
        pnPath.temperature[0].values[m] = dmy[m] * (6.242e18)
        pnPath.temperature[0].coefficients[m][0] = dmy[m] * (6.242e18)
        pnPath.temperature[0].coefficients[m][1] = dmy[m]
        pnPath.temperature[0].values[m] = dmy[m]
        pnPath.density[0].values[m] = dmy[m]
        pnPath.density[0].coefficients[m][0] = dmy[m]
        pnPath.density[0].coefficients[m][1] = dmy[m]
        pnPath.density_fast[0].values[m] = dmy[m]
        pnPath.density_fast[0].coefficients[m][0] = dmy[m]
        pnPath.density_fast[0].coefficients[m][1] = dmy[m]
        pnPath.pressure[0].values[m] = dmy[m]
        pnPath.pressure[0].coefficients[m][0] = dmy[m]
        pnPath.pressure[0].coefficients[m][1] = dmy[m]
        pnPath.pressure_fast_perpendicular[0].values[m] = dmy[m]
        pnPath.pressure_fast_perpendicular[0].coefficients[m][0] = dmy[m]
        pnPath.pressure_fast_perpendicular[0].coefficients[m][1] = dmy[m]
        pnPath.pressure_fast_parallel[0].values[m] = dmy[m]
        pnPath.pressure_fast_parallel[0].coefficients[m][0] = dmy[m]
        pnPath.pressure_fast_parallel[0].coefficients[m][1] = dmy[m]
        pnPath.velocity[0].radial[m] = dmy[m]
        pnPath.velocity[0].radial_coefficients[m][0] = dmy[m]
        pnPath.velocity[0].radial_coefficients[m][1] = dmy[m]
        pnPath.velocity[0].diamagnetic[m] = dmy[m]
        pnPath.velocity[0].diamagnetic_coefficients[m][0] = dmy[m]
        pnPath.velocity[0].diamagnetic_coefficients[m][1] = dmy[m]
        pnPath.velocity[0].parallel_coefficients[m][0] = dmy[m]
        pnPath.velocity[0].parallel_coefficients[m][1] = dmy[m]
        pnPath.velocity[0].poloidal[m] = dmy[m]
        pnPath.velocity[0].poloidal_coefficients[m][0] = dmy[m]
        pnPath.velocity[0].poloidal_coefficients[m][1] = dmy[m]
        pnPath.velocity[0].toroidal[m] = dmy[m]
        pnPath.energy_density_kinetic[0].values[m] = dmy[m]
        pnPath.energy_density_kinetic[0].coefficients[m][0] = dmy[m]
        pnPath.energy_density_kinetic[0].coefficients[m][1] = dmy[m]    

        #--------------------------------------
        # edge_sources
        # electrons
        sePath.particles[0].values[m] = dmy[m] 
        sePath.particles[0].coefficients[m][0] = dmy[m] 
        sePath.energy[0].values[m] = dmy[m] 
        sePath.energy[0].coefficients[m][0] = dmy[m] 

        #--------------------------------------
        # edge_transport
        # electrons
        tePath.particles.d[0].values[m] = dmy[m] + 0.3
        tePath.particles.d[0].coefficients[m][0] = dmy[m] 
        tePath.particles.v[0].values[m] = dmy[m] 
        tePath.particles.v[0].coefficients[m][0] = dmy[m] 
        tePath.particles.flux[0].values[m] = dmy[m] 
        tePath.particles.flux[0].coefficients[m][0] = dmy[m] 
        tePath.particles.flux_limiter[0].values[m] = dmy[m] 
        tePath.particles.flux_limiter[0].coefficients[m][0] = dmy[m] 
        tePath.energy.d[0].values[m] = dmy[m] + 1
        tePath.energy.d[0].coefficients[m][0] = dmy[m] 
        tePath.energy.v[0].values[m] = dmy[m] 
        tePath.energy.v[0].coefficients[m][0] = dmy[m] 
        tePath.energy.flux[0].values[m] = dmy[m] 
        tePath.energy.flux[0].coefficients[m][0] = dmy[m] 
        tePath.energy.flux_limiter[0].values[m] = dmy[m] 
        tePath.energy.flux_limiter[0].coefficients[m][0] = dmy[m] 

        #--------------------------------------
        # ion
        tiPath.particles.d[0].values[m] = dmy[m] + 0.3
        tiPath.particles.d[0].coefficients[m][0] = dmy[m] 
        tiPath.particles.v[0].values[m] = dmy[m] 
        tiPath.particles.v[0].coefficients[m][0] = dmy[m] 
        tiPath.particles.flux[0].values[m] = dmy[m] 
        tiPath.particles.flux[0].coefficients[m][0] = dmy[m] 
        tiPath.particles.flux_limiter[0].values[m] = dmy[m] 
        tiPath.particles.flux_limiter[0].coefficients[m][0] = dmy[m] 
        tiPath.energy.d[0].values[m] = dmy[m] + 1
        tiPath.energy.d[0].coefficients[m][0] = dmy[m] 
        tiPath.energy.v[0].values[m] = dmy[m] 
        tiPath.energy.v[0].coefficients[m][0] = dmy[m] 
        tiPath.energy.flux[0].values[m] = dmy[m] 
        tiPath.energy.flux[0].coefficients[m][0] = dmy[m] 
        tiPath.energy.flux_limiter[0].values[m] = dmy[m] 
        tiPath.energy.flux_limiter[0].coefficients[m][0] = dmy[m] 

    # imas_obj = testDMY(imas_obj)
    return imas_obj

def testDMY(imas_obj):
    imas_obj.edge_profiles.ggd[0].t_i_average.resize(1)
    imas_obj.edge_profiles.ggd[0].n_i_total_over_n_e.resize(1)
    imas_obj.edge_profiles.ggd[0].zeff.resize(1)
    imas_obj.edge_profiles.ggd[0].pressure_thermal.resize(1)
    imas_obj.edge_profiles.ggd[0].pressure_perpendicular.resize(1)
    imas_obj.edge_profiles.ggd[0].pressure_parallel.resize(1)
    imas_obj.edge_profiles.ggd[0].j_inertial.resize(1)
    imas_obj.edge_profiles.ggd[0].j_anomalous.resize(1)
    imas_obj.edge_profiles.ggd[0].j_ion_neutral_friction.resize(1)
    imas_obj.edge_profiles.ggd[0].j_parallel_viscosity.resize(1)
    imas_obj.edge_profiles.ggd[0].j_perpendicular_viscosity.resize(1)
    imas_obj.edge_profiles.ggd[0].j_heat_viscosity.resize(1)
    imas_obj.edge_profiles.ggd[0].j_pfirsch_schlueter.resize(1)
    imas_obj.edge_profiles.ggd[0].j_diamagnetic.resize(1)
    imas_obj.edge_profiles.ggd[0].e_field.resize(1)
    imas_obj.edge_profiles.ggd[0].phi_potential.resize(1)

    imas_obj.edge_profiles.ggd_fast.resize(1)
    return imas_obj

# ======================================================================
# store dumy data for IDS
# 1st argument : imprity data
# 2nd argument : nzmax[Ar,C]
# ======================================================================
def impMome(data, nzmax):
    [tdnz, tfrz, tthz] = data
    [tdnz01, tdnz1, tdnz02, tdnz2] = tdnz
    [tfrz01, tfrz1, tfrz02, tfrz2] = tfrz
    [tthz01, tthz1, tthz02, tthz2] = tthz

    [imome01, imome02,  imome1, imome2] = [[], [], [], []]

    for i in range(len(tdnz01)):
        imome01.append(tdnz01[i]*(tfrz01[i]+tthz01[i]))
        imome02.append(tdnz02[i]*(tfrz02[i]+tthz02[i]))

    for i in range(nzmax[0]):
        imome1.append([])
        for j in range(len(tdnz1[i])):
            imome1[i].append(tdnz1[i][j]*(tfrz1[i][j]+tthz1[i][j]))

    for i in range(nzmax[1]):
        imome2.append([])
        for j in range(len(tdnz2[i])):
            imome2[i].append(tdnz2[i][j]*(tfrz2[i][j]+tthz2[i][j]))
    
    return [imome01, imome02], [imome1, imome2]
