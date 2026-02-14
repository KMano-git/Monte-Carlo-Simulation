!**********************************************************************
      subroutine plrecm(tfsrc)
!**********************************************************************
!
!        cstyp        vli       vol
!    recombination  H+ + e  ==>  H0
!
!    core boundary
!                                    itmax
!            |___|___|___|___|___|___|_O_|
!                                   B.C.
!       region of O  it=itmax
!              kreg(j,i)=6, mrgn(ic)=7, mrgnp(ic)=7
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, iplx, iply, mrgn, ncmax, temin_rec
      use cntwcn, only : swnrm, wabs, wden, wpmp, wssn, wssp, wswe, wswi
     >    , wtot
      use cphcns, only : cev
      use cplcom, only : ama, flvl, nion, nlp, nlpmn, vne, vte, vti, vva
      use cplmet, only : kreg
      use cplqcn, only : mrgnp
      use csize,  only : ndgs, ndmc
      use csonic, only : itim, time
      use cunit,  only : mype
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: tfsrc
      real(8), intent(in) :: tfsrc

! modified 3/3 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ic, ia, j, i, ig
!ik   real(8) :: ane, ate, zsgv, zels, zesr, zwi, zwe, zsn, zsn1
!ik   real(8) :: svhrcm, elhrcm
      integer ic, ia, j, i
      real(8) :: ane, ate, zels, zesr, zwi, zwe, zsn
      real(8) :: elhrcm
      integer :: ift2, ii, ir0, ir1, ir2

!::vli  iflx = 8  trace of D0 in vol
      wden(0:ndmc,1:ndgs) = 0.0d0

!::zero clear  (2006/02/17)
      wssn(0:ndmc,1:ndgs) = 0.0d0
      wssp(0:ndmc,1:ndgs) = 0.0d0
      wswi(0:ndmc) = 0.0d0
      wswe(0:ndmc) = 0.0d0

      do ic = 1, ncmax
        j = iplx(ic)
        i = iply(ic)
        ir2 = mrgn(ic)
        if( j.le.0 .or. i.le.0 .or. ir2 >= 7 ) cycle
        if( vne(j,i).le.0.0d0  ) cycle
        ane  = vne(j,i)
        ate  = vte(j,i)
        ate  = dmax1(ate,temin_rec)  ! svhrcm
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik     zsgv = svhrcm(ate,ane)
        zels = elhrcm(ate,ane)
        zesr = 13.595

        zwi = 0.0d0
        zwe = 0.0d0
        do ia = 1, nion
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik       ig = ia
!xx       wden(ic,ig) = 0.0d0  ! see top of sub. plrecm

!::Note  zsn = -Ne*Ni*<sigv>

!xx       zsn = -wsbr(ic,ig)
          zsn = -flvl(j,i,ia)/(tfsrc/wtot)
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik       zsn1 = zsn
!xx       if( vte(j,i).le.0.6d0 ) zsn1 = 0.0d0     ! 08/06/20

          wssn(ic,ia) = zsn
          wssp(ic,ia) = ama(ia)*vva(j,i,ia)*zsn
          zwi=zwi+(1.5d0*vti(j,i)*cev+0.5d0*ama(ia)*vva(j,i,ia)**2)*zsn
          zwe=zwe+(zels-zesr)*cev*zsn
        enddo  ! loop (ia)

        wswi(ic) = zwi
        wswe(ic) = zwe
      enddo  ! loop(ic)

!::pumping
      wabs = 0.0d0
      wpmp = 0.0d0
      swnrm = wtot


!::debug write
      if( itim == 100 .and. nlp == nlpmn ) then
      ift2 = 22000 + mype
      write(ift2,'(/2x,"*** plsorcA ***  itim =",i7,"  time =",
     >  1pe14.6,"  nlp =",i3,2x,a)') itim, time, nlp, cstyp
      write(ift2,'(4x,a)')
     >  "ic  ir      wden(ic)  wssn(ic)   wssp(ic)   wswi   wswe"
      ii = 0
      do ic = 1, ncmax
        if( wswi(ic) == 0.0d0 ) cycle
        j = iplx(ic)
        i = iply(ic)
        ir0 = kreg(j,i)
        ir1 = mrgn(ic)
        ir2 = mrgnp(ic)

        if( j <= 0 .or. i <= 0 ) cycle
        ii = ii + 1
        if( ii > 30 ) cycle
          write(ift2,'(i6,3i3,1p5e12.4)') ic, ir0, ir1, ir2,
     >      wden(ic,1), wssn(ic,1), wssp(ic,1), wswi(ic), wswe(ic)
      enddo
      endif

!::debug write
      if( mype == 1 .and. itim == 100 .and. nlp == nlpmn ) then
      ift2 = 22000 + mype
      write(ift2,'(/2x,"*** plsorcA ***  itim =",i7,"  time =",
     >  1pe14.6,"  nlp =",i3,2x,a)') itim, time, nlp, cstyp
      write(ift2,'(4x,a)')
     >  "ic  ir      wden(ic)  wssn(ic)   wssp(ic)   wswi   wswe"
      ii = 0
      do ic = 1, ncmax
        if( wswi(ic) == 0.0d0 ) cycle
        j = iplx(ic)
        i = iply(ic)
        ir0 = kreg(j,i)
        ir1 = mrgn(ic)
        ir2 = mrgnp(ic)

        if( j <= 0 .or. i <= 0 ) cycle
        ii = ii + 1
        if( ii > 30 ) cycle
          write(ift2,'(i6,3i3,1p5e12.4)') ic, ir0, ir1, ir2,
     >      wden(ic,1), wssn(ic,1), wssp(ic,1), wswi(ic), wswe(ic)
      enddo
      endif
!==========================

      return
      end
