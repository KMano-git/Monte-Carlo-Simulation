! bcast variables in mod_externalgrid.f
!**********************************************************************
      subroutine exgrid_lsr
!**********************************************************************
      use mod_externalgrid
      use csize, only : ndmc, ndmg, ndms, ndwp, ndx, ndy, ndgtp,ndwl
      use cunit, only : lmspe, mywld
      use mpi
      implicit none

! local variables
      integer :: merr,  ns

! com_ntgrd
      call mpi_bcast( j_save,   ndmg,   MPI_INTEGER,  lmspe
     >                , mywld, merr )
      call mpi_bcast( i_save,   ndmg,   MPI_INTEGER,     lmspe
     >                , mywld, merr )
      call mpi_bcast( ipgt2,   ndgtp,   MPI_INTEGER,     lmspe
     >                , mywld, merr )
      call mpi_bcast( ipwl2,   ndwp,    MPI_INTEGER,     lmspe
     >                , mywld, merr )

! com_vcdat
      call mpi_bcast( vgx_EX,  grid_size,   MPI_REAL8,   lmspe
     >                , mywld, merr )
      call mpi_bcast( vgy_EX,  grid_size,   MPI_REAL8,   lmspe
     >                , mywld, merr )
      call mpi_bcast( mseg_subdiv, cell_size, MPI_INTEGER, lmspe
     >                , mywld, merr )
      ns = cell_size*ndms
      call mpi_bcast( subdiv_cell, ns, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( nvm, 1, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( ksvm, ndem, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( kevm, ndem, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( ivem, ndeg, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( kvem, ndeg, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( cvnm, ndem*2, MPI_CHARACTER, lmspe
     >                , mywld, merr )
      call mpi_bcast( pvem, ndeg*ndem, MPI_INTEGER, lmspe
     >                , mywld, merr )

! externalgrid
      call mpi_bcast( vac_grid_x, vac_grid_size, MPI_REAL8, lmspe
     >                , mywld, merr )
      call mpi_bcast( vac_grid_y, vac_grid_size, MPI_REAL8, lmspe
     >                , mywld, merr )
      call mpi_bcast( pri_grid_x, pri_grid_size, MPI_REAL8, lmspe
     >                , mywld, merr )
      call mpi_bcast( pri_grid_y, pri_grid_size, MPI_REAL8, lmspe
     >                , mywld, merr )
      call mpi_bcast( mseg_vacume, vac_ele_size, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( mseg_pri, pri_ele_size, MPI_INTEGER, lmspe
     >                , mywld, merr )
      ns = vac_ele_size*ndms
      call mpi_bcast( vac_element, ns, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( boundary_vac, ns*2, MPI_CHARACTER, lmspe
     >                , mywld, merr )
      ns = pri_ele_size*ndms
      call mpi_bcast( pri_element, ns, MPI_INTEGER, lmspe
     >                , mywld, merr )
      call mpi_bcast( boundary_pri, ns*2, MPI_CHARACTER, lmspe
     >                , mywld, merr )
      return
      end
