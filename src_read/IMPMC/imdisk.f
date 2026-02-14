!***********************************************************************
      subroutine imdisk(cdsn,cact)
!***********************************************************************
      use cimcom, only : ismax, ndis, ndmp, sflux, swtot
      use cimden, only : csput, fsput, ncmx, nsput, nzmx, tdnz, tfrz
     >    , tionz, tradiz, tradliz, tradrz, trecz, tthz, tvlz, twrd
     >    , twci, tengz, tprz
      use cntcom, only : ncmax, ncmax2
      use cplmet, only : itmpe, itmps, itpve, itpvs, itsle, itsls
      use csize,  only : ndmc, ndy, ndmis
      use cunit,  only : lmspe, lmype, n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      character, intent(in) :: cdsn*(*), cact*(*)
!
!::local variables
      integer  nft
      integer  nsizp, nsizs, nsizc
      integer  it, i, ic, iz, ierr, i_col, i_row
      real*8   twr, twr_rg(10), twr_it(ndy)
      real(8) :: zsum

      if( lmype.ne.lmspe ) then
        call wexit("imdisk","call lmype .ne. lmspe")
      endif
!
      nft = 21
!
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) then
      open(unit=nft, file=cdsn, form="unformatted", action="write")
      write(n6,'(/2x,"*** imdisk ***  ",a,2x,a)') cact, cdsn
!
!     tprz:impurity pressure,  tengz:impurity energy
!     tprz is calculated by ideal gas equation
!     tprz[N/m^2] = e[J/eV] * tdnz[1/m^3] * tengz[eV]
      do i_col = 1, ndmc
        do i_row = 0, ndmis
          tprz(i_row,i_col) = 1.602176634d-19 * 
     >      tdnz(i_row,i_col) * tengz(i_row,i_col)
        enddo
      enddo

      write(nft)  ndmp, ndis, ndmc
      write(nft)  nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
     > ,twci ! add twci for plot tdnz and twrd 22/11/4 yamamoto 
     > ,tengz, tprz
      if( impmc_model == 1 ) write(nft)  tvlz
      close(nft)
!
      if( impmc_model == 1 ) then
        zsum = 0.0d0
        do ic = 1, ncmax2
          do iz = 0, ismax
            zsum = zsum + tvlz(iz,ic)
          enddo
        enddo
        write(n6,'(2x,"### dbg_denz",2x,
     >    "sflux, swtot, tot_tvlz =",1p3e12.3)')
     >    sflux, swtot, zsum
      endif
      return
      endif
!-----------------------------------------------------------------------
!::Yamoto
!-----------------------------------------------------------------------
      if( cact(1:1).eq."y" ) then
      open(unit=nft, file=cdsn, form="unformatted", action="write")
      write(n6,'(/2x,"*** imdisk ***  ",a,2x,a)') cact, cdsn
!
      write(nft)  ndmp, ndis, ndmc
      write(nft)  nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
      if( impmc_model == 0 ) then
        write(nft)  tfrz, tthz, tvlz
        write(nft)  tionZ, trecZ
        write(nft)  tradiZ, tradliZ, tradrZ
      else
        write(nft)  tfrz
        write(nft)  tthz
        write(nft)  tionZ, trecZ
      endif
      close(nft)
!
      return
      endif
!
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      if( cact(1:1).eq."r" ) then
      open(unit=nft, file=cdsn, form="unformatted", action="read")
      write(n6,'(/2x,"*** imdisk ***  ",a,2x,a)') cact, cdsn
!
      read(nft) nsizp, nsizs, nsizc
      write(n6,'(2x,"ndmp =",2i7,"  ndis =",2i4,"  ndmc =",2i7)')
     >   ndmp, nsizp, ndis, nsizs, ndmc, nsizc
      if( ndmp.ne.nsizp .or. ndis.ne.nsizs .or. ndmc.ne.nsizc ) then
      call wexit("imdisk","dimension size error")
      endif
!
      read(nft,iostat=ierr)  nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
     > ,twci, tengz, tprz
      if( impmc_model == 1 ) read(nft,iostat=ierr)  tvlz
      close(nft)
      endif

!-----------------------------------------------------------------------
!::read : Yamoto
!-----------------------------------------------------------------------

      if( cact(1:1).eq."z" ) then
      open(unit=nft, file=cdsn, form="unformatted", action="read")
      write(n6,'(/2x,"*** imdisk ***  ",a,2x,a)') cact, cdsn
!
!
      read(nft)  nsizp, nsizs, nsizc
      write(n6,'(2x,"ndmp =",2i7,"  ndis =",2i4,"  ndmc =",2i7)')
     >   ndmp, nsizp, ndis, nsizs, ndmc, nsizc
      if( ndmp.ne.nsizp .or. ndis.ne.nsizs .or. ndmc.ne.nsizc ) then
      call wexit("imdisk","dimension size error")
      endif
      read(nft)  nsput, nzmx, ncmx, fsput, twrd, tdnz, csput
      if( impmc_model == 0 ) then
        read(nft)  tfrz, tthz, tvlz
        read(nft)  tionZ, trecZ
        read(nft)  tradiZ, tradliZ, tradrZ
      else
        read(nft)  tfrz
        read(nft)  tthz
        read(nft)  tionZ, trecZ
      endif
      close(nft)
!
      return
      endif

!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      call totrgn(twrd, twr, twr_rg, twr_it)
      write(n6,'(2x,"nsput =",i2,"  nzmx =",i2,"  ncmx =",i5)')
     >   nsput, nzmx, ncmx
      write(n6,'(2x,5(a," =",1pe12.3,2x))')
     >  (csput(i), fsput(i), i = 1, nsput)
      write(n6,'(2x,"totwr =",1pe12.3)')   twr
      write(n6,'(2x,"wr_rg =",1p10e12.3)') twr_rg
      write(n6,'(2x,"wr_sl =",1p10e12.3)') (twr_it(it),it=itsls,itsle)
      write(n6,'(2x,"wr_pv =",1p10e12.3)') (twr_it(it),it=itpvs,itpve)
      write(n6,'(2x,"wr_ed =",1p10e12.3)') (twr_it(it),it=itmps,itmpe)
      write(n6,'(2x,"twrd  =",1p10e12.3)') (twrd(ic),ic=1,ncmax,100)

      if( impmc_model == 1 ) then
        zsum = 0.0_8
        do ic = 1, ncmax2
          do iz = 0, ismax
            zsum = zsum + tvlz(iz,ic)
          enddo
        enddo
        write(n6,'(2x,"### dbg_denz",2x,
     >    "sflux, swtot, tot_tvlz =",1p3e12.3)')
     >    sflux, swtot, zsum
      endif

      return
      end
