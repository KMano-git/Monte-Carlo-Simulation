!**********************************************************************
      subroutine dtypdef(ctyp,ierr)
!**********************************************************************
      use mod_dtypdef, only : nddt, typ_dnam, typ_prog
      use mod_sizedef, only : lnnam
      implicit none

! modified 2/2 lines organize local variables and include files by kamata 2021/05/31
!ik   character(lnnam) :: ctyp
!ik   integer :: ierr
      character(lnnam), intent(in)  :: ctyp
      integer,          intent(out) :: ierr

! local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: kerr, i, ityp
      integer :: kerr
! added 4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
      integer      nc
      character    code*10, np*2
! code : code name
! nc   : code number
! np   : processing number

      select case(trim(ctyp))
      case("MST_1")
        call def_MST_1( kerr )  ! MST_1  lstep  Mnkall_ini.f
      case("MST_2")
        call def_MST_2( kerr )  ! MST_2  time   Mnkall_tim.f

      case("PLS_1")
        call def_PLS_1( kerr )  ! PLS_1  aion   Mnkntl_ini.f
      case("PLS_1B")
        call def_PLS_1B( kerr )  ! PLS_1B  mdl_wrd  Mnkimp_ini.f
      case("PLS_2")
        call def_PLS_2( kerr )  ! PLS_2   vna      Mnkimp_pls.f
      case("PLS_2B")
        call def_PLS_2B( kerr )  ! PLS_2B  vna      Mnkntl_pls.f
      case("PLS_3")
        call def_PLS_3( kerr )       ! PLS_3   vna      Mnkntl_pls.f

      case("NTL_1")
        call def_NTL_1( kerr )  ! NTL_1  xpnt   Mnkntl_iniSnd.f
      case("NTL_2")
        call def_NTL_2( kerr )     ! NTL_2  vtime  Mnkntl_cal.f

! deleted 6 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   case("IMP_1")
!ik     call def_lnkimp_iniSnd( kerr )  ! IMP_1  xdaty  Mnkimp_iniSnd.f
!ik   case("IMP_2")
!ik     call def_lnkimp_cal( kerr )     ! IMP_2  nsput  Mnkimp_cal.f
!ik   case("IMP_2B")
!ik     call def_IMP_2B( kerr )         ! IMP_2B zw_catmz

! deleted 4 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   case("IMP2_1")
!ik     call def_IMP2_1( kerr )         ! IMP_1  xdaty  Mnkimp_iniSnd.f
!ik   case("IMP2_2")
!ik     call def_IMP2_2( kerr )          ! IMP_2  nsput  Mnkimp_cal.f

! deleted 6 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik   case("IMP3_1")
!ik     call def_IMP3_1( kerr )         ! IMP_1  xdaty  Mnkimp_iniSnd.f
!ik   case("IMP3_2")
!ik     call def_IMP3_2( kerr )          ! IMP_2  nsput  Mnkimp_cal.f
!ik   case("IMP3_2B")
!ik     call def_IMP3_2B( kerr )         ! IMP_2B zw_catmz2

      case default
! added 29 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
        if(  ctyp(1:3) == 'IMP' ) then
! get code number and processing number
          call getncnp( ctyp, nc, np )

! execution
          if( nc > 0 .and. np /= ' ' ) then
            if( nc == 1 ) then
              code = 'IMPMC'
            else
              write(code,'(i5)') nc
              code = adjustl( code )
              code = 'IMPMC' // trim( code )
            endif
            select case( np )
            case( '1 ' )
              call def_IMP_1( ctyp, code, kerr )  ! IMP_1  xdaty  Mnkimp_iniSnd.f
            case( '2 ' )
              call def_IMP_2( ctyp, code, kerr )     ! IMP_2  nsput  Mnkimp_cal.f
            case( '2B' )
              nddt = nddt + 1
              typ_dnam(nddt) = ctyp
              typ_prog(nddt) = code
              kerr = 0
            case default
              ierr = 9
              call wf_exit("dtypdef","no found dtypdef "//trim(ctyp))
            end select
          endif
        else
          ierr = 9
          call wf_exit("dtypdef","no found dtypdef "//trim(ctyp))
! added 1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
        endif
      end select

!xx      ntyp = ityp

!::debug write
!xx      write(n6,'(2x,"ntyp =",i5,"  ndtyp =",i5)') ntyp, ndtyp
!xx      do i = 1, ntyp
!xx        write(n6,'(2x,i3,2x,a,2x,a,i3,i15,i2,i15)') i,typ_dnam(i),
!xx     >  typ_prog(i),typ_igrp(i),typ_itag(i),typ_ierr(i),typ_adr0(i)
!xx      enddo

      ierr = 0
      if( kerr /= 0 ) ierr = 1

      return
      end


!**********************************************************************
      subroutine entry_dtyp(ityp)
!**********************************************************************
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: ityp
      integer, intent(inout) :: ityp

      ityp = ityp + 1

! deleted 5 lines generalization of derived data type by kamata 2020/10/31
!ik   typ_dnam(ityp) = tydnam
!ik   typ_prog(ityp) = typrog
!ik   typ_itag(ityp) = tyitag
!ik   typ_ierr(ityp) = tyierr
!ik   typ_adr0(ityp) = tyadr0

      return
      end
