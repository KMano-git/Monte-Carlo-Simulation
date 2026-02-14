!***********************************************************************
      subroutine imp_prof
!***********************************************************************
!)
!)   initial loading of IMPMC  time evolution version
!)
!)   input of cutbk
!)       tmcut :   time start     0.3 sec
!)       dtcut :   every dtcut    0.1 sec
!)       nctst :   npmax > nctst  1200   ct:cut,  st:start
!)
!---------------------------------------------------------------------
      use cimcom,      only : dtscr, dttgt, hgrpn, hgrpp, hgrpw, hregn
     >    , hregp, hregw, i6_scat, itmmc, nctst, ndmp, npmax, scpum
     >    , sitmz, stimz, tmmcz, tmref, tmsim
      use cimctl,      only : icalz
      use cimden,      only : nsput, nsrc_spt
      use cimfls,      only : lfls
      use cimntl,      only : ndmp2
      use cimpuf,      only : bk_npt, bk_nty
      use csonic,      only : itim, limp
      use cunit,       only : n6
      use mod_loc_tim, only : timem
      use mpi!,         only : mpi_wtime
      implicit none
!
!::local variables
      real(8) :: tnow
      integer :: itm, i6
      real(8) :: cpu0, cpu1
      integer :: sv_bk_npt
      integer :: bk_call
!
      if( limp == 0 ) return

!::dimension check
      if( ndmp /= ndmp2 ) then
        write(n6,'(/2x,"dimension error in imp_prof  ndmp /= ndmp2",
     >    2i7)') ndmp, ndmp2
        call wexit("imp_prof","ndmp /~ ndmp2")
      endif

      call imp_file("open")


!::partcile number
      bk_call = 0
      tnow = 0.0d0
      itm = 0
      cpu0 = MPI_WTIME()
      write(n6,*) '=== imp_prof: START ===', tnow, nsrc_spt

!:: SY load at first loop of nsput
      if(nsrc_spt.eq.1)then
        call imp_rstfl(0)
        call imp_rstfl(2)
      endif

!:: loop time
      do
        if( dttgt.eq.0.d0) goto 380
        if( tnow+1.0d-12 > dttgt) exit

!::trace => denZ
!::cal. bk_mag
        if( bk_nty.gt.0 .and. nsrc_spt.eq.1 .and. bk_call.eq.0) then
          sv_bk_npt = bk_npt
!-------
          write(n6,'(2x,"SS-impcal imbkflx ",3i6)') itim, icalZ, bk_npt
          call imstpbk(tnow,dtscr)
          bk_npt = sv_bk_npt
          bk_call = 1
        endif

        call imp_step(tnow,dtscr)

        itm = itm + 1
        tnow = tnow + dtscr
        stimz = tnow
        sitmz = itm
        cpu1  = MPI_WTIME()
        scpum = (cpu1-cpu0)/60.0d0
        cpu0 = cpu1
 380    continue
        call lst_mscat(tnow)

!::particle conservation
        call imp_wclas(hgrpp,hgrpw,hgrpn,hregp,hregw,hregn)
        call imp_print("pcon")  ! 1PE
        call imp_mpsum("pcon")
        call imp_print("pcnT")  ! TPE

!::denZ(MC scoreing) => dnz, wrd
        call imp_mpsum("denz")
        call imp_mpsum("flos")

!::averaged quantities during dtsiz
        call imp_phyv !SY imdens relevant
        call imp_print("wrad")
        call imwflx(nsrc_spt)
        if(lfls.eq.1) call lstflos("flux")
!SY test particle number control
!SY imp_pack MUST be called in advance
!SY otherwise wexit
        call imp_pack ! SY kill dead particles (ien(ip).ge.6 for instance)
        if( npmax > nctst ) then
          call imp_cutbk
!        tmcut = tmcut + dtcut !SY this does never work with multiple source
        endif
        if( dttgt.eq.0.d0) exit
      enddo  !  loop till tim > dttgt

!::preserve output file
      call flush(n6)

!::local time of IMPMC => common time of soldor
!SY this routine also requires imp_pack.
!SY this routine must be called after finishing whole the nsput loop
      if(nsrc_spt.eq.nsput) call imp_reset

!::check
      i6 = i6_scat
      if( i6 > 0 ) then
        call lst_scat(tnow,1,i6)  ! called in imp_prof for HK
        close(i6)
      endif

!::time control  in imp_exec step
      tmsim = timeM
      tmref = tmsim
      tmmcz = 0.0d0
      itmmc = 0

      write(n6,'(2x,"imp_prof  tnow =",1pe12.3,"  tmref/tmsim =",
     >  1pe15.8,1pe12.4,"  tmmcz =",1pe12.4,"  itmmc =",i6)')
     >   tnow, tmref, tmsim, tmmcz, itmmc

      return
      end
