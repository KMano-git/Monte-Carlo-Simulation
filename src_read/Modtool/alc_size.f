! added dynamic allocation of arrays by kamata 2022/05/29
! set the size required for allocate
      subroutine alc_size
      use cimcom,      only : nwkmp_im1, nwkmp_im1y, nwkmps_im1
     >    , nwkmps_im1y
      use csize,       only : mdsp, ndbr, ndmc, ndmg, ndmis, ndsp, ndwp
     >    , ndx, ndxy, ndy, nvyd, nwkmp_nt, nzsmx
      use mod_mpicomm, only : grpprg, m6, mygrp, ngrp
      use mod_shexe,   only : set_csize
      implicit none
! local variables
      integer    i

! set from the value read by &shexe
!- default size is csize_68
      call set_csize

! set nzsmx for IMPMC
      nzsmx = 0
      do i = 1, ngrp
        if( grpprg(i)(1:5) == 'IMPMC' ) nzsmx = nzsmx + 1
      enddo
      nzsmx = max( nzsmx, 1 ) ! for SOLDOR

! set mdsp in csize
      mdsp  = ndmis * nzsmx + ndsp

! set cimcom
      nwkmp_im1   = 12 + (ndmis+1) * (ndmc+1) * 3 + (ndmc+1)
      nwkmp_im1y  = (ndmis+1) * (ndmc+1) * 4
      nwkmps_im1  = 10 + (ndmis+1) * (ndmc+1) * 2 + (ndmc+1)
      nwkmps_im1y = (ndmis+1) * (ndmc+1) * 5

! output array size
      if( mygrp == 1 ) then ! SOLDOR
        write(m6,900) 'ndbr  = ', ndbr, 'ndmc  = ', ndmc
     >    , 'ndmg  = ', ndmg
        write(m6,900) 'ndwp  = ', ndwp, 'nzsmx = ', nzsmx
     >    , 'mdsp  = ', mdsp
        write(m6,900) 'ndx   = ', ndx,  'ndy   = ', ndy
     >    , 'ndxy  = ', ndxy
        write(m6,900) 'nvyd  = ', nvyd, 'nwkmp_nt = ', nwkmp_nt
      endif

      return
  900 format( 2( a, i10, 2x ), a, i10 )
      end
