!***********************************************************************
      subroutine imbkflx
!***********************************************************************
      use cimcom, only : bkfrac, flximpamp, fwcal, fwstp, htout
     >    , lgstk, ltmax, lwstk
      use cimctl, only : lstdy
      use cimden, only : csput, npsput, nsput
      use cimpuf, only : bk_mag, bk_npt, bk_nty, lbkstw, pf_mag, pf_nty
      use cntcom, only : chwl, nhwl, nogt
      use csize,  only : ndwh
      use csonic, only : itim
      use cunit,  only : n6
      use mpi!,    only : mpi_wtime
      implicit none
!
!::local variables
      integer :: sv_nsput
      integer, dimension(5) :: sv_npsput
!
      integer :: sv_lgstk,  sv_ltmax
      real(8) :: sv_fwstp
      integer, dimension(ndwh) :: sv_lwstk
!
      real(8) :: totpf,  cptm0, cptm1
      integer :: ih, i
      character(4) :: cmsg
!
      cptm0 = MPI_WTIME()
!
!::variables for output
      bk_nty = nogt
      totpf = 0.0d0
      if(csput(1).eq."Arpf") then
        do i = 1, pf_nty
          totpf = totpf + pf_mag(i)
        enddo
      elseif(csput(1).eq."Cedg")then
        totpf = flxImpAmp
      endif
!
!::input data
      if( bk_npt.le.0 ) then
        cmsg = "inpt"
        write(n6,'(2x,"imbkflx  itim =",i7,2x,a,"  totpf =",1pe12.4,
     >  "  bk_mag =",1p10e12.4)') itim, cmsg, totpf, bk_mag(1:bk_nty)
      return
      endif
!
!::input
      sv_nsput = nsput
      sv_npsput(1:5) = npsput(1:5)
      sv_lgstk  = lgstk
      sv_fwstp  = fwstp
      sv_ltmax  = ltmax
      sv_lwstk(1:ndwh) = lwstk(1:ndwh)
!
      nsput = 1
      npsput(1) = bk_npt
      lgstk  = 0
      lbkstw = 0
      fwstp  = 0.20   !  stop condition
!
      lwstk(1:ndwh) = 0
      if( lgstk.eq.1 ) then
        do ih = 1, nhwl
          if( chwl(ih)(3:3).eq."g" ) lwstk(ih) = 1
          if( chwl(ih)(1:1).eq."g" ) lwstk(ih) = 1
        enddo
      endif
!
      write(n6,'(/2x,"*** imbkflx ***  cal.  lstdy =",i2)') lstdy
      write(n6,'(2x,5(1x,a4,1x,i7,2x))') (csput(i),npsput(i),i=1,nsput)
      write(n6,'(2x,"lgstk  = ",i2,"  fwstp = ",f8.4,"  ltmax =",i8)')
     >   lgstk, fwstp, ltmax
      write(n6,'(2x,"chwl = ",20(2x,a,1x))') (chwl(ih),ih=1,nhwl)
      write(n6,'(2x,"lwstk= ",20(2x,i2.2,3x))') (lwstk(ih),ih=1,nhwl)
!
!::calculation
      call immont
!
!::back flow
      if( bk_npt.gt.0 ) then
        do i = 1, bk_nty
          bk_mag(i) =  totpf*bkfrac(i)
        enddo
      endif
!
      cptm1 = MPI_WTIME()
      cmsg  = "cal."
      write(n6,'(2x,"imbkflx  itim =",i7,2x,a,"  htout =",f8.3,
     >  "  cpu =",f12.3,"  fwcal/fwstp =",2f8.3,"  totpf =",1pe12.4,
     >  "  bk_mag =",1p10e12.4)') itim, cmsg, htout, cptm1-cptm0,
     >  fwcal, fwstp, totpf, bk_mag(1:bk_nty)
!
!::input  (recovery)
      nsput = sv_nsput
      npsput(1:5) = sv_npsput(1:5)
      lgstk  = sv_lgstk
      fwstp  = sv_fwstp
      ltmax  = sv_ltmax
      lwstk(1:ndwh) = sv_lwstk(1:ndwh)
!
      write(n6,'(/2x,"*** imbkflx ***  input data of  imp_cal")')
      write(n6,'(2x,5(1x,a4,1x,i7,2x))') (csput(i),npsput(i),i=1,nsput)
      write(n6,'(2x,"lgstk  = ",i2,"  fwstp = ",f8.4,"  ltmax =",i8)')
     >   lgstk, fwstp, ltmax
      write(n6,'(2x,"chwl = ",15(2x,a,1x))') (chwl(ih),ih=1,nhwl)
      write(n6,'(2x,"lwstk= ",15(2x,i2.2,3x))') (lwstk(ih),ih=1,nhwl)
!
      return
      end
