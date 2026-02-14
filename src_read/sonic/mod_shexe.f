!.. The variables in this module were stored in shexe.
!.. These variables were moved in inppls.
!***********************************************************************
      module mod_shexe
!***********************************************************************
      use mod_mpicomm, only : port_inppls
      use cntwedf,     only : nene_mx, nsmx, nwmx, nene_mx_line
      use csize,       only : ndbr, ndgs, ndmc, ndmg, ndwp, ndx, ndxy
     >    , ndy, nvyd, nwkmp_nt, ndgtp
      use mod_externalgrid, only: use_exdata, grid_size, vac_grid_size
     > , pri_grid_size, cell_size, vac_ele_size, pri_ele_size
     > , is_compareEX
      implicit none

      logical:: use_namelist

      integer:: OMP_NUM_THREADS

      character(20):: BTYPE, INP_PLS, INP_NTL
      character(20):: INP_IMP, ERROR
      character:: WHAT_DAY*11, WHAT_TIM*8, WHAT_LOD*80
      character:: WHAT_DIR*80, WHAT_JOB*10

      character(50):: REFIT, REQIL, MESH_MTR, MESH_NTL, MESH_SIZE
     > , MESH_DIR, DEGAS_REFL, DEGAS_ATMHY
     > , DEGAS_MOLHY, DEGAS_ELATM, DEGAS_ELMOL
      data REFIT, REQIL, MESH_MTR, MESH_NTL
     > , MESH_SIZE /"","","","",""/
      integer :: npe_log = 1

! impmc_model : IMPMC model ( = 0 : steady-state, = 1 : time dependent )
      integer :: impmc_model = 0

      namelist /impmc/ impmc_model

      namelist /comsize/ ndmc, ndmg, ndx, ndy, nvyd, ndgtp

      namelist /shexe/
     > REFIT, REQIL, MESH_MTR, MESH_NTL, MESH_SIZE, MESH_DIR,
     > use_exdata,
     > DEGAS_REFL, DEGAS_ATMHY, DEGAS_MOLHY,
     > DEGAS_ELATM, DEGAS_ELMOL, OMP_NUM_THREADS,
     > BTYPE, INP_PLS, INP_NTL, INP_IMP, ERROR,
     > WHAT_DAY, WHAT_TIM, WHAT_LOD,
     > WHAT_DIR, WHAT_JOB
     > , ndmc, ndmg, ndx, ndy, nvyd, ndgtp ! dynamic allocation of arrays
     > , npe_log
     > , grid_size,vac_grid_size,pri_grid_size
     > , cell_size,vac_ele_size,pri_ele_size
     > , is_compareEX

!**********************************************************************
      contains
!***********************************************************************
      subroutine set_shexe_variables
!***********************************************************************
      integer :: error_int
      logical :: file_exist

        open(unit=port_inppls,file='inppls')
        read(unit=port_inppls, nml=shexe, iostat=error_int)
        if(error_int == 0) then
           use_namelist = .TRUE.
           !::set mesh data path
           if(trim(REQIL)=="") then
             REQIL = trim(MESH_DIR)//"/dat.eqil"
             inquire(file=REQIL,EXIST=file_exist)
             if(.not. file_exist) REQIL = ""
           endif
           !::
           if(trim(REFIT)=="") REFIT = trim(MESH_DIR)//"/dat.efit"
           if(trim(MESH_MTR)=="") 
     >       MESH_MTR = trim(MESH_DIR)//"/dat.mtrc"
           if(trim(MESH_NTL)=="") 
     >       MESH_NTL = trim(MESH_DIR)//"/dat.mont"
           if(trim(MESH_SIZE)=="") 
     >       MESH_SIZE = trim(MESH_DIR)//"/comsize"
        else
           use_namelist = .FALSE.
           MESH_DIR = ""
           REFIT = ""
           REQIL = ""
           MESH_MTR = ""
           MESH_NTL = ""
           MESH_SIZE = ""
           use_exdata = .False.
           DEGAS_REFL = ""
           DEGAS_ATMHY = ""
           DEGAS_MOLHY = ""
           DEGAS_ELATM = ""
           DEGAS_ELMOL = ""
           OMP_NUM_THREADS = 1
           BTYPE = ""
           INP_PLS = ""
           INP_NTL = ""
           INP_IMP = ""
           ERROR = ""
           WHAT_DAY = ""
           WHAT_TIM = ""
           WHAT_LOD = ""
           WHAT_DIR = ""
           WHAT_JOB = ""
           npe_log = 1
        end if
! read values from inppls
        rewind port_inppls
        read(unit=port_inppls, nml=impmc, iostat=error_int)
        if( impmc_model /= 0 ) impmc_model = 1
        close(port_inppls)
! read ndmc from comsize
        open(unit=port_inppls,file=trim(MESH_SIZE),iostat=error_int)
        read(unit=port_inppls, nml=comsize, iostat=error_int)
        close(port_inppls)
      end subroutine set_shexe_variables

!***********************************************************************
      subroutine set_csize
!***********************************************************************
      ndbr = ndmc
      ndwp = 2 * ( ndx + ndy )
      if( ndmc == 7600 .and. ndmg == 7900
     >    .and. ndx == 160 .and. ndy == 80 ) ndwp = 600
      ndxy = max( ndx, ndy )
      nwkmp_nt = 1                   ! wflx
     >   + (ndmc+1)*(ndgs*2+2+ndgs)  ! wssn(,ndgs),wssp(,ndgs),wswe,wswi,wsbr(,ndgs)
     >   + (ndmc+1)*(ndgs*3)         ! wden(,ndgs),weng(,ndgs),wvlp(,ndgs)
     >   + (ndmc+1)*(ndgs*3)         ! wnfl0x(,ndgs),wnfl0y(,ndgs),wnfl0z(,ndgs) 
     >   + (ndmc+1)*(ndgs+2*2)       ! wtion(,ndgs),wgden(,2),wgeng(,2)
     >   + (ndmc+1)*(2*3)            ! wnflgx(,2),wnflgy(,2),wnflgz(,2)
     >   + ndwp*2+10+1+30*3*2        ! wemt(),wwal(),wend(10),wtot,whta(30,3),whtm(30,3)
     >   + ndwp*4                    ! wfhta(),wfhtm(),wehta(),wehtm()
     >   + 11+6                      ! wreg(0:10),wsum,wnrm,wion,wabs,wpmp,werr
     >   + nene_mx*nwmx*nsmx         ! wedf, see Cntwedf
     >   + ndwp*3                    ! weema(), weemm(), wees()
     >   + nene_mx_line*nwmx*nsmx    ! wedf_line, see Cntwedf

      end subroutine set_csize

      end module mod_shexe
