!**********************************************************************
      subroutine imwcon(kk)
!**********************************************************************
!
!       kk = 0  loop set back   value for 1PE
!       kk = 1  last sumup      value for total PE
!
!
!     imhist : delete   2012/11/17  K. Shimizu
!
!      imhist  :: wght => hcnpt, hcnwt  /cimcom_12b/  wwork3
!      imcndv  :: wght => wvacn         /cimcom_12c/  wwork4
!
!      imwcon  :: wtemt, whit           /cimcom_12/   wwork2
!      imwcon  :: wght => hcnpt, hcnwt  /comcom_12b/  wwork3
!
!
!     do k = 1, nsput
!
!     do lt = 1, ltmax
!       call imrqout(lt,iprf,ihst)
!       call imstep
!        if( ihst.eq.1 ) then
!          call imcndv
!     !xx  call imhist(lt,kend)
!        endif
!        if( iprf.eq.1 ) then
!          call improf   ! scoreing denZ temZ for each source
!        endif
!     enddo
!
!      call imcndv
!      call imhist(lt,kend)  ! Note put back wwork3 for lmspe
!      call imwcon(1)        ! Note after imhist
!
!      call imdens
!      call imsave(sptyc,"w")
!      enddo   ! source
!
!
!c::pumped case
!c      ln = il(ip)
!c      iw = -ln
!c      ih = ihwl(iw)
!c         if( chwl(ih)(1:1).eq."a" .or.     !  a1
!c   >         chwl(ih)(3:3).eq."a" ) then   !  V1a1 <== no meaning
!
!----------------------------------------------------------------------
!::/cimcom_12/ wtemt, wpemt, whit, wexh
!      nwkmp_im2    wwork2,  vwork2      at imwcon (last cal.)
!      these data is important.     Normalization of Monte-result
!
! These scored variables cimcom_12 should be summed after calculation.
!
!
!::/cimcom_12b/ hmype, htout, hstp, ... hrgpa, hrgpc, hcnpt, hcnwt
!      nwkmp_im3    wwork3,  vwork3      at  imhist (loop)
!                                        and imwcon (last)
!
!::/cimcom_12c/ wvacn, wvaci
!      nwkmp_im4    wwork4,  vwork4      at imcndv (loop)
!
!      set the values before call imwcon
!
!        flemt : impurity flux
!        wtemt : total weight of sample particles
!
!----------------------------------------------------------------------
      use cimcom, only : hcnpt, hcnwt, hrgpa, hrgpc, idcal, idemt, iderr
     >    , idetc, idexh, idstk, ien, ir, ismax, itimz, lhst, lpcm, lprf
     >    , npmax, ns2 => nwkmp_im2, nst => nwkmp_im3, nss => nwkmps_im3
     >    , sptyc, timez, wexh, wght, wpemt, wtemt, mrgn_max
      use cimctl, only : lstdy
      use cntcom, only : mrgn, nhwl
      use cunit,  only : lmspe, lmype, lnope, mype, mywld, n6
      use mod_shexe, only : impmc_model
      use mpi!,    only : mpi_real8, mpi_reduce, mpi_sum
      implicit none
!
!::argument
      integer, intent(in) :: kk
!
!::mpi variables
      real(8), dimension(ns2) :: wwork2, vwork2, twork2
      real(8), dimension(nst) :: wwork3, vwork3, twork3 ! nst > nss
!
!::local variables
      integer :: irg, ip, ic, iz, ih, icnd
      integer :: ierr
      real(8) :: zwt, zsma, zsmc, ztexh
      real(8) :: wgt0
      integer    ns3
      integer :: ih60 = 0
      save ::  ih60
!
      write(n6,'(2x,"*** imwcon ***  lpcm =",i10,"  lprf =",i2,
     >  "  lhst =",i2)') lpcm, lprf, lhst
!
!----------------------------------------------------------------------
!::classify test particles  for conservation
!----------------------------------------------------------------------
!   emt, cal, exh/pmp, err/etc, stk   same process in imhist
!
!::wght when emitted
      if( impmc_model == 0 ) then
        wgt0 = 1.0d0    ! see imemit
!
!::hrgpa, hrgpc  rg:region, pa:allpartcile  pc:cal-partcile
!::hcnpt, hcnwt  cn:condition, pt:partcile, wt:weight
        hrgpa(0:10) = 0.0d0   ! region for all particle
        hrgpc(0:10) = 0.0d0   ! region for cal-particle  ien = 0
        hcnpt(0:10) = 0.0d0   ! condition of all particles
        hcnwt(0:10) = 0.0d0   ! condition of all particles
!
        hcnpt(idemt) = wpemt
        hcnwt(idemt) = wtemt
!
        do ip = 1, npmax
          icnd = ien(ip)
          zwt  = wght(ip)
          hcnwt(idexh) = hcnwt(idexh) + wgt0-zwt
          if( icnd.eq.0 ) then
            hcnpt(idcal) = hcnpt(idcal) + 1
            hcnwt(idcal) = hcnwt(idcal) + zwt
          elseif( icnd.eq.6 .or. icnd.eq.7 ) then
            hcnpt(idstk) = hcnpt(idstk) + 1
            hcnwt(idstk) = hcnwt(idstk) + zwt
            hcnwt(idexh) = hcnwt(idexh) + zwt       ! No trace
          elseif( icnd.eq.9 ) then
            hcnpt(iderr) = hcnpt(iderr) + 1
            hcnwt(iderr) = hcnwt(iderr) + zwt
          else
            hcnpt(idetc) = hcnpt(idetc) + 1
            hcnwt(idetc) = hcnwt(idetc) + zwt
          endif
          ic  = ir(ip)
          irg = mrgn(ic)
          hrgpa(irg) = hrgpa(irg) + zwt
          if( icnd.eq.0 ) hrgpc(irg) = hrgpc(irg) + zwt
        enddo
!
!--------------------------------------
!::debug write   (imwcon for each PE)
!-------------------------------------
        if( lprf.eq.1 ) then
          write(n6,'(2x,"1PE_P: emt =",1pe12.4,"  cal =",1pe12.4,
     >    "  exh/stk =",1p2e12.4,"  err/etc =",1p2e12.4)')
     >     hcnpt(idemt), hcnpt(idcal), hcnpt(idexh), hcnpt(idstk),
     >     hcnpt(iderr), hcnpt(idetc)
        endif
      endif
!
!::save  work2 work3
      call setvwork2( 2, wwork2 )
      if( impmc_model == 0 ) then
        ns3 = nss
        call setvwork3( 0, 2, wwork3, ns3 )
        twork2(1:ns2) = wwork2(1:ns2)
        twork3(1:ns3) = wwork3(1:ns3)
      else
        ns3 = nst
        call setvwork3( 1, 2, wwork3, ns3 )
      endif

!----------------------------------------------------------------------
!::sumup
!----------------------------------------------------------------------
!::[MPI_Reduce in imwcon]  cimcom_12   (wtemt,wexh)  10/05/26
      call MPI_Reduce( wwork2, vwork2, ns2, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
!
!::[MPI_Reduce in imwcon]  cimcom_12b  (hmype,hcnwt) 10/05/26
      call MPI_Reduce( wwork3, vwork3, ns3, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
!
      if( lmype.ne.lmspe ) goto 100
!
!::Do not sumup lmype, tout
      vwork3(1) = wwork3(1)  ! lmype
      vwork3(2) = wwork3(2)  ! tout
      vwork3(3) = vwork3(3)/dfloat(lnope)  ! averaged step number
!
      call setvwork2( 1, vwork2 )
      call setvwork3( impmc_model, 1, vwork3, ns3 )
!
!----------------------------------------------------------------------
!::conservation of wght
!----------------------------------------------------------------------
!
!::total
      if( impmc_model == 0 ) then
        zsma = 0.0d0
        zsmc = 0.0d0
        do irg = 0, mrgn_max
          zsma = zsma + hrgpa(irg)
          zsmc = zsmc + hrgpc(irg)
        enddo
!
!::total exhaust
        ztexh = 0.0d0
        do ih = 1, nhwl+1
          do iz = 0, ismax
            ztexh = ztexh + wexh(iz,ih)
          enddo
        enddo
!
!--------------------------------------
!::debug write   (imwcon for all PE)
!-------------------------------------
        if( lprf.eq.1 ) then
          write(n6,'(2x,"Tot_P: emt =",1pe12.4,"  cal =",1pe12.4,
     >    "  exh/stk =",1p2e12.4,"  err/etc =",1p2e12.4)')
     >     hcnpt(idemt), hcnpt(idcal), hcnpt(idexh), hcnpt(idstk),
     >     hcnpt(iderr), hcnpt(idetc)
!
          write(n6,'(2x,"Tot_W: emt =",1pe12.4,"  cal =",1pe12.4,
     >      "  exh/stk =",1p2e12.4,"  err/etc =",1p2e12.4,
     >      " :  wex/nrm =",1p2e12.4)')
     >       hcnwt(idemt), hcnwt(idcal), hcnwt(idexh), hcnwt(idstk),
     >       hcnwt(iderr), hcnwt(idetc), ztexh, hcnwt(idcal)+ztexh
!
          write(n6,'(2x,"TRg_a: tot/err =",1p2e12.4,"  odv/sol/idv =",
     >    1p3e12.4,"  opv/ipv =",1p2e12.4,"  edg/cor/vac =",1p3e12.4)')
     >    zsma, hrgpa(0), (hrgpa(irg),irg=1,8)
!
          write(n6,'(2x,"TRg_c: tot/err =",1p2e12.4,"  odv/sol/idv =",
     >    1p3e12.4,"  opv/ipv =",1p2e12.4,"  edg/cor/vac =",1p3e12.4)')
     >    zsmc, hrgpc(0), (hrgpc(irg),irg=1,8)
        endif
      endif
!
 100  continue
      if( impmc_model .ne. 0 ) return
      if( lstdy.ge.1 .and. lmype.eq.lmspe .and. lhst.eq.1 ) then
        if( ih60.eq.0 ) then
          ih60 = 380000 + 1000 + mype
          open(unit=ih60,file="Hist_imwcon.txt")
          write(ih60,'(2x,"sptyc",5x,"lpcm",3x,"itimZ",2x,"timeZ",7x,
     >        "TPemt",7x,"TPcal",7x,"TPexh",7x,"TPstk",7x,
     >        "TPerr",7x,"TPetc",7x,
     >        "TWemt",7x,"TWcal",7x,"TWexh",7x,"TWstk",7x,
     >        "TWerr",7x,"TWetc",7x,"TWwex",7x,"TWnrm",7x,
     >        "RAall",7x,"RAerr",7x,"RAodv",7x,"RAsol",7x,"RAidv",7x
     >        "RAopv",7x,"RAipv",7x,"RAedg",7x,"RAcor",7x,"RAvac",7x,
     >        "RCall",7x,"RCerr",7x,"RCodv",7x,"RCsol",7x,"RCidv",7x
     >        "RCopv",7x,"RCipv",7x,"RCedg",7x,"RCcor",7x,"RCvac")')
        endif
!
        write(ih60,'(2x,a4,i10,i7,1pe12.4,  1p6e12.4, 1p8e12.4,
     >      1p10e12.4, 1p10e12.4)')
     >     sptyc, lpcm, itimZ, timeZ
     >    ,hcnpt(idemt), hcnpt(idcal), hcnpt(idexh), hcnpt(idstk),
     >      hcnpt(iderr), hcnpt(idetc)
     >    ,hcnwt(idemt), hcnwt(idcal), hcnwt(idexh), hcnwt(idstk),
     >      hcnwt(iderr), hcnwt(idetc), ztexh, hcnwt(idcal)+ztexh
     >    ,zsma, hrgpa(0), (hrgpa(irg),irg=1,8)
     >    ,zsmc, hrgpc(0), (hrgpc(irg),irg=1,8)
      endif

!::set back  wwork2, wwork3
      if( kk.eq.0 ) then
        call setvwork2( 1, twork2 )
        call setvwork3( 0, 1, twork3, ns3 )
      endif

      return
      end
