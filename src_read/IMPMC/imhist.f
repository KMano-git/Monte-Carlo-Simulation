!***********************************************************************
      subroutine imhist(lptm,kend)
!***********************************************************************

!::Teese data are produced based on the particle information
!::sumup in imhist     wght(ip) => xcnpt, xcnwt

!!    /cimcom_12b/
!!      real*8  hmype, htout, hstep, hcpum, hpmax
!!      real*8  hrgpa(0:10), hrgpc(0:10) ! region of partciles
!!      real*8  hcnpt(0:10), hcnwt(0:10) ! condition of particle
!!      parameter (nwkmp_im3 = 5 + 11*4 )

!::Note
!     Variables in cimcom_12b (hcnpt,hcnwt) should be
!       put back after sumup. because sumup is carried out in imwcon.
!
!     But  cimcom_12b ==> com_imhist (local common) in imhist
!        no need to put back   ???   2012/11/17   Shimizu
!
!     istpm : maximum step number ==> averaged step number
!     donot sum up   w3(1) = lmype, w3(2) = tout
!-----------------------------------------------------------------------
      use cimcom, only : fwcal, fwstp, hcpum, hgrpn, hgrpp, hgrpw, hmype
     >    , hpmax, hregn, hregp, hregw, hstep, htout, idcal, idemt
     >    , iderr, idetc, idstk, igcal, igemt, igerr, igexh, igmin
     >    , ignul, igpmp, igrun, igsce, igsgt, lhst, nst => nwkmp_im3
     >    , nss => nwkmps_im3, sptyc, twabs, twcal, twemt, twerr
      use csonic, only : itim
      use cunit,  only : lmspe, lmype, lnope, mywld, n6
      use clocal, only : xcnpt, xcnwt, xcpum, xmype, xpmax, xrgpa, xrgpc
     >    , xstep, xtout
      use mod_shexe, only : impmc_model
      use mpi!,    only : mpi_bcast, mpi_real8, mpi_recv, mpi_reduce
!    >    , mpi_send, mpi_status_size, mpi_sum
      implicit none

!::argument
      integer, intent(in)  :: lptm ! dummy
      integer, intent(out) :: kend

!::local variables
      integer :: i60
      integer :: imdp, idbg, k, i, nsiz
      integer :: ipcal, ipwab, iperr, ipetc, ipemt, iptlm, imype
      real(8) :: istpm, xtmp(4)
      integer :: imypeb

      real(8), dimension(nst) :: xwork3, vwork3, twork3 ! nst > nss

!::mpi variables
      integer  nbf, spe, rpe, ierr
      integer  istat(MPI_STATUS_SIZE)
      integer lhist

      if( impmc_model == 0 ) then
        nbf = nss
        lhist = 1
      else
        nbf = nst
        lhist = 0 ! cannot be changed
      endif

      i60 = n6

      imdp = lnope/5
      if( imdp.eq.0 ) imdp = 1

!::classify test particles
      if( impmc_model == 0 ) then
        call imwclas(xrgpa,xrgpc,xcnpt,xcnwt)

!::add information of imstep
        xmype = hmype
        xtout = htout
        xstep = hstep
        xcpum = hcpum
        xpmax = hpmax
      else
        call imp_wclas(hgrpp,hgrpw,hgrpn,hregp,hregw,hregn)
      endif
      call setvwork3( impmc_model*2-1, 2, xwork3, nbf )
! first argument = -1 if impmc_model = 0, = 1 if impmc_model = 1

!::clear V and save W
      vwork3(1:nbf) = 0.0d0
      twork3(1:nbf) = xwork3(1:nbf)

      if(lhist.eq.0)then

        imypeb = -1
        do k = 1, lnope
          spe = k - 1
          rpe = lmspe
          if( k.eq.1 ) goto 100

!::[MPI_Send in imhist]  com_imhist (xmype,xcnwt)        12/08/09
          if( lmype.eq.spe ) then
            call MPI_Send( xwork3, nbf, MPI_REAL8, rpe
     >       , 3002,mywld,ierr )
          endif

!::[MPI_Recv in imhist]  com_imhist (xmype,xcnwt)        12/08/09
          if( lmype.eq.rpe ) then
            call MPI_Recv( xwork3, nbf, MPI_REAL8, spe, 3002,mywld,
     >        istat,ierr )
            call setvwork3( impmc_model*2-1, 1, xwork3, nbf )
! first argument = -1 if impmc_model = 0, = 1 if impmc_model = 1
          endif

 100      continue
!::ERR 2010/07/02  mod(hmyoe,imdp) ==>
          idbg = 0
          if( k.eq.1 .or. mod(k,imdp).eq.0 ) idbg = 1

!::output
          if( idbg.eq.1 ) then
!::see define these variables at last imstep
!xx   xmype, xtout, xstep, xcpum, xpmax
            if( impmc_model == 0 ) then
              ipcal = int(xcnpt(idcal))
              ipwab = int(xcnpt(idstk))  ! last condition
              iperr = int(xcnpt(iderr))
              ipetc = int(xcnpt(idetc))
              ipemt = int(xcnpt(idemt))
              istpm = xstep
              iptlm = int(xpmax + 1.0d-8)
              imype = int(xmype + 1.0d-8)
            else
              istpm = hstep
              imype = int(hmype + 1.0d-8)
            endif

            if( imype.ne.imypeb .and. lhst.eq.1 ) then
              if( impmc_model == 0 ) then
                write(i60,'(3x,a,i6,2x,a,i7,1pe12.3,2x,a,3i7,2x,a,3i7,
     >            2x,a,0pf8.1,0pf9.3)')
     >            "imhistH  "//sptyc, imype, "tim", itim, xtout,
     >            "np/cl/st", iptlm, ipcal, ipwab,  ! npmax, cal, stk
     >            "er/et/em", iperr, ipetc, ipemt,  ! error, etc, emit
     >            "avst/cpu",  istpm, xcpum
              else
                write(i60,'(1x,a,i6,2x,a,i7,1pe12.3, 2x,a,2i6, 2x,a,
     >            2i6,2x,a,4i6,2x,a,2i6,2x,a,0pf8.1,0pf10.3)')
     >            "imhist  "//sptyc, imype, "tim", itim, htout,
     >            "em/cl", int(hgrpp(igemt)), int(hgrpp(igcal)),
     >            "gt/ce", int(hgrpp(igsgt)), int(hgrpp(igsce)),
     >            "mn/er/ex/pm", int(hgrpp(igmin)), int(hgrpp(igerr)),
     >             int(hgrpp(igexh)), int(hgrpp(igpmp)),
     >            "nl/rn",  int(hgrpp(ignul)), int(hgrpp(igrun)),
     >            "avst/cpu", istpm, hcpum
             endif
            endif ! lhst
          endif ! idbg

          if( lmype.eq.lmspe ) then
            vwork3(1) = xwork3(1)                  ! lmype
            vwork3(2) = xwork3(2)                  ! tout
            do i = 3, nbf
              vwork3(i) = vwork3(i) + xwork3(i)
            enddo
          endif

          imypeb = imype
        enddo  ! loop(k)

      else   ! lhist = 1
        call MPI_Reduce( xwork3, vwork3, nbf, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
      endif  ! lhist

      if( lmype.ne.lmspe ) goto 200

!::common variables of lmype
      call setvwork3( impmc_model*2-1, 1, vwork3, nbf )
! first argument = -1 if impmc_model = 0, = 1 if impmc_model = 1

!::total
      if( impmc_model == 0 ) then
        xcpum = xcpum/dfloat(lnope)
        ipcal = int(xcnpt(idcal))
        ipwab = int(xcnpt(idstk))
        iperr = int(xcnpt(iderr))
        ipetc = int(xcnpt(idetc))
        ipemt = int(xcnpt(idemt))
        istpm = xstep/dfloat(lnope)
        iptlm = int(xpmax + 1.0d-8)

        fwcal = xcnwt(idcal)/xcnwt(idemt)

        if( lhst.eq.1 .or. lhist.eq.1) then
          write(i60,'(3x,a,a6,2x,a,i7,1pe12.3,2x,a,3i7,2x,a,3i7,
     >    2x,a,0pf8.1,0pf9.3)')
     >    "imhistL  "//sptyc, "   T", "tim", itim, xtout,
     >    "np/cl/st", iptlm, ipcal, ipwab,  ! npmax, cal, stk
     >    "er/et/em", iperr, ipetc, ipemt,  ! error, etc, emit
     >    "avst/cpu",  istpm, xcpum
        endif

!::set back
        call setvwork3( -1, 1, twork3, nbf )
      else
!::total  xcpum => hcpum  xnpt => hgrpp
        hcpum = hcpum/dfloat(lnope)
        istpm = hstep/dfloat(lnope)
        twemt = hgrpw(igemt)
        twcal = hgrpw(igcal)
        twabs = hgrpw(igsgt)+hgrpw(igsce)+hgrpw(igmin)+hgrpw(igexh)
        twerr = hgrpw(igerr)

        xtmp(1) = twemt
        xtmp(2) = twcal
        xtmp(3) = twabs
        xtmp(4) = twerr

        if( lhst.eq.1 ) then
          write(n6,'(1x,a,a6,2x,a,i7,1pe12.3, 2x,a,2i6, 2x,a,2i6,
     >    2x,a,4i6,2x,a,2i6,2x,a,0pf8.1,0pf10.3)')
     >    "imhist  "//sptyc, "   T", "tim", itim, htout,
     >    "em/cl", int(hgrpp(igemt)), int(hgrpp(igcal)),
     >    "gt/ce", int(hgrpp(igsgt)), int(hgrpp(igsce)),
     >    "mn/er/ex/pm", int(hgrpp(igmin)), int(hgrpp(igerr)),
     >       int(hgrpp(igexh)), int(hgrpp(igpmp)),
     >    "nl/rn",  int(hgrpp(ignul)), int(hgrpp(igrun)),
     >    "avst/cpu", istpm, hcpum
        endif
      endif

 200  continue
      if( impmc_model == 0 ) then
        nsiz = 1
        call MPI_Bcast( fwcal, nsiz, MPI_REAL8, lmspe, mywld, ierr )
      else
        nsiz = 4
        call MPI_Bcast( xtmp, nsiz, MPI_REAL8, lmspe, mywld, ierr )
  
        twemt = xtmp(1)
        twcal = xtmp(2)
        twabs = xtmp(3)
        twerr = xtmp(4)
        if( twemt .ne. 0.0d0) then
          fwcal = twcal/twemt
        else
          fwcal = 0.0d0
        endif
        write(n6,*) 'fwcal', fwcal, fwstp
      endif

      kend = 0
      if( fwcal.le.fwstp ) kend = 1

      if( impmc_model == 0 ) then
        if( lhst.eq.1 .or. kend.eq.1 .or. lhist.eq.1) then
          write(i60,'(3x,a,2x,"fwcal/fwstp =",0p2f8.4,"  kend =",i2)')
     >    "imhist", fwcal, fwstp, kend
        endif
      else
        if( lhst.eq.1 .or. kend.eq.1 ) then
          write(i60,'(1x,a,a6,2x,a,i7,
     >      2x,a,0p2f7.4,2x,a,i2,2x,a,1p8e11.4)')
     >      "imhist  "//sptyc, " CON", "tim", itim,
     >      "fwcal =", fwcal, fwstp, "kend =", kend,
     >      "hgrpw = ",hgrpw(1:8)
        endif

!::set back
!   because  sumup weit and partcile number at every time
        if( lmype == lmspe ) then
          call setvwork3( 1, 1, twork3, nbf )
        endif
      endif

      return
      end
