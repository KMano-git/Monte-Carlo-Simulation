!***********************************************************************
      subroutine set_adas
!***********************************************************************
!
!     iclass  dsn           meaning
!  ----------------------------------------------------
!      (1)  acd96_c.data : recombination rate
!      (2)  scd96_c.data : ionization rate
!      (3)  ccd96_c.data : CX-transfer
!      (4)  prb96_c.data : Power recombination
!      (5)  prc96_c.data : Power continuum
!      (6)  qcd96_c.data : ??
!      (7)  xcd96_c.data : ??
!      (8)  plt96_c.data : Power ionization
!      (9)  pls96_c.data : Power charge transfer
!  ----------------------------------------------------
!
! Header of data
!    6   24   30    1    6     /CARBON             /GCR PROJECT
!
! Header of data  acd96r_he.dat
!    2   24   30    1    2     /HELIUM             /GCR PROJECT
!    2    1    1          meta-stable
!-----------------------------------------------------------------------
      use catcom, only : ched, dratm, ip_cxr, ip_ion, ip_plt, ip_prb
     >    , ip_prc, ip_rec, nfatm, rtyp, tbdsn
      use cimctl, only : nprg
      use cunit,  only : cdgr, lmspe, lmype, lnope, mygrp, n6
      implicit none
!
!::mpi variables
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  nls, nle, nbf, ierr
      integer  ierr
!
!::local variables
! modified 2/1 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ipp, ityp, iz, i
!ik   character(80) :: cdir
      integer  ipp, i
      character(5), dimension(7) :: comt
! modified 1/1 lines multiple IMPMC execution bug by kamata 2022/04/15
!ik   character(10) :: dtyp
      character(10) :: cid, dtyp
!
      write(n6,'(/2x,"*** set_adas ***")')
      write(n6,'(2x,"ionization & recombination cross section",
     >  " calculated by ADAS-code  2001/4/20 12/01/06")')
!
!::tbdsn
      if( lmype.eq.lmspe ) then
! added 2 lines multiple IMPMC execution bug by kamata 2022/04/15
        write(cid,'(i9)') nprg
        cid = adjustl( cid )
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call getls(dratm,"f",nfatm,tbdsn,10)
      call getls(dratm(1),"f",nfatm,tbdsn,10,'_'//cid)
!
!::table
!
      comt=(/"RT-i ","RT-r ","RT-cx","PW-i ","NOuse","PW-r ","PW-cx"/)
      rtyp=(/"ION", "REC", "CXR",  "PLT", "RAD",  "PRB", "PRC"/)
      ched=(/"scd", "acd", "ccd",  "plt", "plt",  "prb", "prc"/)
!     ipNO=   ion    rec    cxr     plt    xxx     prb    prc
!
      write(n6,'(2x,"comt = ",10(a5,2x))') (comt(i),i=1,7)
      write(n6,'(2x,"rtyp = ",10(a5,2x))') (rtyp(i),i=1,7)
      write(n6,'(2x,"ched = ",10(a5,2x))') (ched(i),i=1,7)
!
      write(n6,'(2x)')
      write(n6,'(2x,"adas-file = ",10(a,2x))')
     >   (trim(tbdsn(i)),i=1,nfatm)
!
!::example : how to get dsn-name
!xx  call atm_fnam("ATOM_ION",dsn)
      endif
!
!::clear
! modified 6/6 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ip_ion = 0
!ik   ip_rec = 0
!ik   ip_cxr = 0
!ik   ip_plt = 0
!ik   ip_prb = 0
!ik   ip_prc = 0
      ip_ion(1) = 0
      ip_rec(1) = 0
      ip_cxr(1) = 0
      ip_plt(1) = 0
      ip_prb(1) = 0
      ip_prc(1) = 0
!
!::read file of adas-data
      if( lmype.eq.lmspe ) then
      call atm_copy("clear",ipp)
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_copy("ATOM_ION",ip_ion)  ! rate_ion
!ik   call atm_copy("ATOM_REC",ip_rec)  ! rate_rec
!ik   call atm_copy("ATOM_CXR",ip_cxr)  ! rate_cxr
      call atm_copy("ATOM_ION",ip_ion(1))  ! rate_ion
      call atm_copy("ATOM_REC",ip_rec(1))  ! rate_rec
      call atm_copy("ATOM_CXR",ip_cxr(1))  ! rate_cxr
!xx   call atm_copy("ATOM_RAD",ip_rad)  ! <== No use
! modified 3/3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_copy("ATOM_PLT",ip_plt)  ! rad_line
!ik   call atm_copy("ATOM_PRB",ip_prb)  ! rad_rec
!ik   call atm_copy("ATOM_PRC",ip_prc)  ! rad_cxr  <== No use 12/01/06
      call atm_copy("ATOM_PLT",ip_plt(1))  ! rad_line
      call atm_copy("ATOM_PRB",ip_prb(1))  ! rad_rec
      call atm_copy("ATOM_PRC",ip_prc(1))  ! rad_cxr  <== No use 12/01/06
      write(n6,'(2x,"------ atm_copy in set_adas")')
      endif
!KH      write(n6,'(2x,"------ atm_copy in set_adas")')
!KH      call flush(n6)
!
!::[MPI_Bcast in set_adas]  catcom (xdaty,emrk)  10/04/21
!xx      if( lnope.gt.1 ) then
!xx      nls = loc( xdaty(1) )
!xx      nle = loc( catcom_emrk ) + 1
!xx      nbf = nle - nls
!xx      call MPI_Bcast( xdaty, nbf, MPI_BYTE, lmspe, mywld, ierr )
!xx      write(n6,'(2x,"------ MPI_Bcast in set_adas")')
!xx      endif

      if( lnope  > 1 ) then
!      call dtyplnk("IMP_1",cdgr(3),cdgr(3),ierr)
! modified 3/6 lines treat 4 or more impurities with IMPMC by kamata 2022/04/21
!ik     if(mygrp.eq.cdgr(3)) dtyp = "IMP_1"
!ik     if(mygrp.eq.cdgr(4)) dtyp = "IMP2_1"
!ik     if(mygrp.eq.cdgr(5)) dtyp = "IMP3_1"
        cid = ' '
        if( nprg > 1 ) then
          write(cid,'(i5)') nprg
          cid = adjustl( cid )
        endif
        dtyp = 'IMP' // trim( cid ) // '_1'
        call dtyplnk(dtyp,mygrp,mygrp,ierr)
      endif
!
!::debug write
      write(n6,'(2x,"  ip_ion =",i2,"  ip_rec =",i2,"  ip_cxr =",i2,
     > "  ip_plt =",i2,"  ip_prb =",i2,"  ip_prc =",i2)')
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik  >  ip_ion, ip_rec, ip_cxr, ip_plt, ip_prb, ip_prc
     >  ip_ion(1), ip_rec(1), ip_cxr(1), ip_plt(1), ip_prb(1), ip_prc(1)
      call flush(n6)
!
!::ip_prc (Power CXR) No use
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   if( ip_ion.eq.0 .or. ip_rec.eq.0 .or.
!ik  >    ip_plt.eq.0 .or. ip_prb.eq.0 ) goto 910
      if( ip_ion(1) == 0 .or. ip_rec(1) == 0 .or.
     >    ip_plt(1) == 0 .or. ip_prb(1) == 0 ) goto 910
!
!::set interp-variables
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_pset
      call atm_pset( 1 )
      write(n6,'(2x,"------ atm_pset in set_adas")')
!
      return
!
 910  continue
      call wexit("set_adas","ip_ion/rec/plt/prb/ = 0")
      end
