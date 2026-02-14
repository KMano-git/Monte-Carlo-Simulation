!**********************************************************************
      subroutine imp_exec
!**********************************************************************
!)
!)   sample particles are followed dring dtscr (5 msec)
!)
!)   input of cutbk
!)       tmcut :   time start     0.3 sec
!)       dtcut :   every dtcut    0.1 sec
!)       nctst :   npmax > nctst  1200   ct:cut,  st:start
!)
!---------------------------------------------------------------------
      use cimcom, only : dtscr, hgrpn, hgrpp, hgrpw, hregn, hregp, hregw
     >    , itmmc, nctst, ndmp, npmax, scpum, sitmz, stimz, tmmcz
      use cimctl, only : icalz
      use cimden, only : nsput, nsrc_spt
      use cimfls, only : lfls
      use cimntl, only : ndmp2
      use cimpuf, only : bk_npt, bk_nty
      use csonic, only : itim, limp
      use cunit,  only : n6
      use mpi!,    only : mpi_wtime
      implicit none

!::local variables
      real(8) :: tnow
      integer :: itm
      real(8) :: cpu0, cpu1
      integer :: sv_bk_npt

      cpu0=0.d0

      if( limp == 0 ) return

!::dimension check
      if( ndmp /= ndmp2 ) then
        write(n6,'(/2x,"dimension error in imp_prof  ndmp /= ndmp2",
     >    2i7)') ndmp, ndmp2
        call wexit("imp_prof","ndmp /~ ndmp2")
      endif

!::loop start
      tnow = tmmcz
      itm = itmmc

!:: Backflow
      if(bk_nty.gt.0 .and. nsrc_spt.eq.1 .and. mod(icalZ,10).eq.2) then
        sv_bk_npt = bk_npt
!-------
        write(n6,'(2x,"SS-impcal imbkflx ",3i6)') itim, icalZ, bk_npt
!-------
        call imstpbk(tnow,dtscr)
        bk_npt = sv_bk_npt
      else
        write(n6,'(2x,"Skip SS-impcal imbkflx ",2i6)') itim, icalZ
      endif

!::trace => denZ
      call imp_step(tnow,dtscr)

!::local var. in sub. imp_exec
!::incremented in the last nsput loop
      if(nsrc_spt.eq.nsput)then
        itm  = itm + 1
        tnow = tnow + dtscr
!::common in IMPMC
        tmmcz = tnow
        itmmc = itm
        stimz = tnow
        sitmz = itm

        cpu1  = MPI_WTIME()
        scpum = (cpu1-cpu0)/60.0d0
        cpu0 = cpu1
      endif

!::particle conservation
      call imp_wclas(hgrpp,hgrpw,hgrpn,hregp,hregw,hregn)
      call imp_print("pcon")  ! 1PE
      call imp_mpsum("pcon")
      call imp_print("pcnT")  ! TPE

!::denZ(MC scoreing) => dnz, wrd
      call imp_mpsum("denz")
      call imp_mpsum("simp")
      call imp_mpsum("flos")
      call imwcon( -1 ) !SY 14072020 ! added dummy argument

!::averaged quantities during dtsiz
      call imp_phyv ! imdens relevant
      call imp_print("wrad")
      call imwflx(nsrc_spt)
      if(lfls.eq.1) call lstflos("flux")
      call imp_pack
! Test particle # control
      if( npmax > nctst ) then
        call imp_cutbk
      endif

!::loop end

!::preserve output file
      call flush(n6)

      return
      end
