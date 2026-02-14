!**********************************************************************
      subroutine dtypset( dtyp, mtyp )
!**********************************************************************
      implicit none

      character(*), intent(in) :: dtyp
      integer, intent(in) :: mtyp ! 1:pre-Bcst,2:post-bcast,3:Send,4:Recv

! added 29 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
! local variables
      integer         mj, nc
      character(2)    np
! mj : number of characters
! nc : code number
! np : processing number

! check dtyp
      if( mtyp /= 4 ) return
      mj = len_trim( dtyp )
      if( mj < 5 ) return

! get nc and np from 'IMP(nc)_(np)' and execution
      if( dtyp(1:3) == 'IMP' ) then
! get processing number and code number
        call getncnp( dtyp, nc, np )

! execution
        if( nc > 0 .and. np /= ' ' ) then
          select case( np )
          case( '1 ' )
            call set_plimp_1( nc )
          case( '2 ' )
            call set_plimp_2( nc )
          case default
            call wexit( 'dtypset', 'dtyp is wrong '//trim( dtyp ) )
          end select
        endif
      endif

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   character(lnnam) :: ctyp

! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   ctyp = trim(dtyp)

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   if(ctyp == "IMP_1"  .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   if(dtyp == "IMP_1"  .and. mtyp == 4) then
!ik      call set_plimp_1(1)
! deleted 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik      call set_atcom_A(1)
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   elseif(ctyp == "IMP_2"  .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   elseif(dtyp == "IMP_2"  .and. mtyp == 4) then
!ik      call set_plimp_2(1)
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   elseif(ctyp == "IMP2_1" .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   elseif(dtyp == "IMP2_1" .and. mtyp == 4) then
!ik      call set_plimp_1(2)
! deleted 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik      call set_atcom_A(2)
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   elseif(ctyp == "IMP2_2" .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   elseif(dtyp == "IMP2_2" .and. mtyp == 4) then
!ik      call set_plimp_2(2)

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   elseif(ctyp == "IMP3_1" .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   elseif(dtyp == "IMP3_1" .and. mtyp == 4) then
!ik      call set_plimp_1(3)
! deleted 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik      call set_atcom_A(3)
!
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   elseif(ctyp == "IMP3_2" .and. mtyp == 4) then
! deleted 2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   elseif(dtyp == "IMP3_2" .and. mtyp == 4) then
!ik      call set_plimp_2(3)

! deleted 1 line treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   endif
!   2018/10/18 Y.Homma to write out plasma parameters.
!     Overwritten at each iteration, and the latest communicated
!     data to IMPMC remain finally.
!     USE FOR FINAL SMOOTHING PHASE
!      if(ctyp == "PLS_2" .and. mtyp == 3) then
!         call yh_sld_plasma_parameters_1_plot
!         call heat_flux_along_flux_tube(24) ! Separatrix
!         call heat_flux_radial_output(1, 0, 1) ! Odiv
!         call heat_flux_radial_output(1, 0, 40) ! OX
!         call heat_flux_radial_output(1, 0, 50) ! OX-Omid
!         call heat_flux_radial_output(1, 0, 65) ! Omid
!         call heat_flux_radial_output(1, 0, 84) ! Imid
!         call heat_flux_radial_output(1, 0, 100) ! IX-Imid
!         call heat_flux_radial_output(1, 0, 109) ! IX
!         call heat_flux_radial_output(1, 0, 148) ! Idiv
!         call heat_flux_radial_output(0, 1, 1) ! Write_all, all cells
!      endif
!      if(mtyp == 1)then     ! pre-Bcast
!      if(mtyp == 2)then     ! post-Bcast
!      if(mtyp == 3)then     ! Send
!      if(mtyp == 4)then     ! Recv
!      endif

      return
      end
