!**********************************************************************
      subroutine dtyplnk(dtyp,sgrp,rgrp,kcon)
!**********************************************************************
!
!::Error  adr1 2383290 => variables  Ghd_div
!     call MPI_Send( adr1,  1, itag, rpe, mark, nwld_cmm, merr )
!
!     Do NOT use : adr0, typ_adr0
!
!---------------------------------------------------------------------
      use cimctl,      only : tmimp
      use cimden,      only : nsput
      use cntcom,      only : grax
      use cntmnt,      only : vtime
      use cplcom,      only : aion, mdl_wrd, tfps
!  PLS_1:aion   PLS_1B:mdl_wrd  PLS_2:vna  PLS_2B:vna
      use cplwrd,      only : wfac
      use csonic,      only : lstep, time
      use cunit,       only : mygrp
      use mod_dtflw,   only : ntyp, typ_dnam, typ_itag
      use mod_mpicomm, only : nmst_grp, nrnk_cmm, nwld_cmm, nwld_grp
      use mod_sizedef, only : lnnam
      use mpi!,         only : mpi_bcast, mpi_recv, mpi_send
!    >    , mpi_status_size
      implicit none

      character(*), intent(in)  :: dtyp
      integer,      intent(in)  :: sgrp, rgrp
      integer,      intent(out) :: kcon ! = 0 only

! local variables
      integer :: ind, itag
      integer :: ityp
      integer :: spe, rpe, merr, mark, istat(MPI_STATUS_SIZE)
      integer         nc
      character(2)    np
! nc  : code number
! np  : processing number

      kcon = 0
      call tbfind( ntyp, typ_dnam, ityp, dtyp )
      if( ityp <= 0 ) then
      call wf_exit("dtyplnk","no found dtyp in typ_dnam "//trim(dtyp))
      endif

      itag = typ_itag(ityp)
      spe = sgrp
      rpe = rgrp
      mark = 300 + ityp

!::MPI_Bcast
      if( spe == rpe .and. mygrp == sgrp ) then
      call dtypset(dtyp,1)
      spe = nmst_grp
      select case(trim(dtyp))
      case("MST_1")
                  call MPI_Bcast(lstep, 1, itag,spe,nwld_grp,merr)
      case("MST_2")
                  call MPI_Bcast(time,  1, itag,spe,nwld_grp,merr)
      case("PLS_1")
                  call MPI_Bcast(aion,  1, itag,spe,nwld_grp,merr)
                  call cgdcom_sr( 1, spe )
                  call cplcom_sr( dtyp, 1, spe )
                  call cplmet_sr( 1, spe )
                  call cpmpls_sr( 1, spe )
      case("PLS_1B")
                  call MPI_Bcast(mdl_wrd, 1, itag,spe,nwld_grp,merr)
      case("PLS_2")
                  call MPI_Bcast(tfps,  1, itag,spe,nwld_grp,merr)
                  call cplcom_sr( dtyp, 1, spe )
                  call cntwfl_sr( 1, spe )
      case("PLS_2B")
                  call MPI_Bcast(tfps,  1, itag,spe,nwld_grp,merr)
                  call cplcom_sr( dtyp, 1, spe )
      case("PLS_3")
                  call MPI_Bcast(wfac,   1, itag,spe,nwld_grp,merr)
      case("NTL_1")
                  call MPI_Bcast(grax,  1, itag,spe,nwld_grp,merr)
                  call cntcom_sr( 1, spe )
      case("NTL_2")
                  call MPI_Bcast(vtime, 1, itag,spe,nwld_grp,merr)
                  call cntmnt_sr( 1, spe )
      case default
        if(  dtyp(1:3) == 'IMP' ) then
          call getncnp( dtyp, nc, np )
! execution
          if( nc > 0 .and. np /= ' ' ) then
            select case( np )
            case( '1 ' )
              call MPI_Bcast(tmimp, 1, itag,spe,nwld_grp,merr)
              call catcom_sr( 1, spe, 1, 1 )
              call cimcom_sr( 1, spe )
            case( '2 ' )
              call MPI_Bcast(nsput, 1, itag,spe,nwld_grp,merr)
              call cimden_sr( 1, spe )
            case( '2B' )
              call czwflx_sr( 1, spe, nc, nc )
            case default
              call wf_exit("dtyplnk","no found MPI_Bcast("//trim(dtyp))
            end select
          endif
        else
           call wf_exit("dtyplnk","no found MPI_Bcast("//trim(dtyp))
        endif
      end select
      call dtypset(dtyp,2)

!::MPI_Send
      elseif( nrnk_cmm == spe ) then
      call dtypset(dtyp,3)
      select case(trim(dtyp))
      case("MST_1")
                  call MPI_Send(lstep, 1,itag,rpe,mark,nwld_cmm,merr)
      case("MST_2")
                  call MPI_Send(time,  1,itag,rpe,mark,nwld_cmm,merr)
      case("PLS_1")
                  call MPI_Send(aion, 1,itag,rpe,mark,nwld_cmm,merr)
                  call cgdcom_sr( 2, rpe )
                  call cplcom_sr( dtyp, 2, rpe )
                  call cplmet_sr( 2, rpe )
                  call cpmpls_sr( 2, rpe )
      case("PLS_1B")
                  call MPI_Send(mdl_wrd, 1,itag,rpe,mark,nwld_cmm,merr)
      case("PLS_2")
                  call MPI_Send(tfps, 1,itag,rpe,mark,nwld_cmm,merr)
                  call cplcom_sr( dtyp, 2, rpe )
                  call cntwfl_sr( 2, rpe )
      case("PLS_2B")
                  call MPI_Send(tfps, 1,itag,rpe,mark,nwld_cmm,merr)
                  call cplcom_sr( dtyp, 2, rpe )
      case("PLS_3")
                  call MPI_Send(wfac, 1,itag,rpe,mark,nwld_cmm,merr)
      case("NTL_1")
                  call MPI_Send(grax,1,itag,rpe,mark,nwld_cmm,merr)
                  call cntcom_sr( 2, rpe )
      case("NTL_2")
                  call MPI_Send(vtime,1,itag,rpe,mark,nwld_cmm,merr)
                  call cntmnt_sr( 2, rpe )
      case default
        if(  dtyp(1:3) == 'IMP' ) then
! get code number and processing number
          call getncnp( dtyp, nc, np )
! execution
          if( nc > 0 .and. np /= ' ' ) then
            select case( np )
            case( '1 ' )
              call MPI_Send(tmimp,1,itag,rpe,mark,nwld_cmm,merr)
              ind = nc
              if( spe > 2 ) ind = 1
              call catcom_sr( 2, rpe, ind, 1 )
              call cimcom_sr( 2, rpe )
            case( '2 ' )
              call MPI_Send(nsput,1,itag,rpe,mark,nwld_cmm,merr)
              call cimden_sr( 2, rpe )
            case( '2B' )
              call czwflx_sr( 2, rpe, nc, 1 )
            case default
              call wf_exit("dtyplnk","no found MPI_Send("//trim(dtyp))
            end select
          endif
        else
          call wf_exit("dtyplnk","no found MPI_Send("//trim(dtyp))
        endif
      end select

!::MPI_Recv
      elseif( nrnk_cmm == rpe ) then
      select case(trim(dtyp))
      case("MST_1")
             call MPI_Recv(lstep, 1,itag,spe,mark,nwld_cmm,istat,merr)
      case("MST_2")
             call MPI_Recv(time, 1,itag,spe,mark,nwld_cmm,istat,merr)
      case("PLS_1")
             call MPI_Recv(aion, 1,itag,spe,mark,nwld_cmm,istat,merr)
             call cgdcom_sr( 3, spe )
             call cplcom_sr( dtyp, 3, spe )
             call cplmet_sr( 3, spe )
             call cpmpls_sr( 3, spe )
      case("PLS_1B")
             call MPI_Recv(mdl_wrd, 1,itag,spe,mark,nwld_cmm,istat,merr)
      case("PLS_2")
             call MPI_Recv(tfps, 1,itag,spe,mark,nwld_cmm,istat,merr)
             call cplcom_sr( dtyp, 3, spe )
             call cntwfl_sr( 3, spe )
      case("PLS_2B")
             call MPI_Recv(tfps, 1,itag,spe,mark,nwld_cmm,istat,merr)
             call cplcom_sr( dtyp, 3, spe )
      case("PLS_3")
             call MPI_Recv(wfac, 1,itag,spe,mark,nwld_cmm,istat,merr)
      case("NTL_1")
             call MPI_Recv(grax,1,itag,spe,mark,nwld_cmm,istat,merr)
             call cntcom_sr( 3, spe )
      case("NTL_2")
             call MPI_Recv(vtime,1,itag,spe,mark,nwld_cmm,istat,merr)
             call cntmnt_sr( 3, spe )
      case default
        if(  dtyp(1:3) == 'IMP' ) then
          call getncnp( dtyp, nc, np )

! execution
          if( nc > 0 .and. np /= ' ' ) then
            select case( np )
            case( '1 ' )
              call MPI_Recv(tmimp,1,itag,spe,mark,nwld_cmm,istat,merr)
              ind = nc
              if( rpe > 2 ) ind = 1
              call catcom_sr( 3, spe, 1, ind )
              call cimcom_sr( 3, spe )
            case( '2 ' )
              call MPI_Recv(nsput,1,itag,spe,mark,nwld_cmm,istat,merr)
              call cimden_sr( 3, spe )
            case( '2B' )
              call czwflx_sr( 3, spe, 1, nc )
            case default
              call wf_exit("dtyplnk","no found MPI_Recv("//trim(dtyp))
            end select
          endif
        else
          call wf_exit("dtyplnk","no found MPI_Recv("//trim(dtyp))
        endif
      end select
      call dtypset(dtyp,4)

      endif

      kcon = 0
      return
      end subroutine dtyplnk
