!***********************************************************************
      subroutine pldens(cmsg)
!***********************************************************************
!
!  mid-plane (out)  xmid =  4.0726  ymid = -0.0693  jc,ic =   52   30
!  mcel( 3,30), mcel(52,30), mcel(52,35), mcel( 3,32)
!  tden0, teng0 (in NTL) => tN0, tT0 (in PLS)  tN0=sum dotN*wden(ic)
!  teng0 = 2/3*1/2*M*V^2/[eV] = T0
!
!----------------------------------------------------------------------
      use cntcom, only : cstyp, den0, deng, eng0, engg, ncmax2, vlp0
      use cntmnt, only : mcflx, mfmax, visrc
      use cplcom, only : nion, snflx, tN0, tNg, tT0, tTg, tV0
      use csize,  only : ndgs, ndmc
      implicit none

      character, intent(in) :: cmsg*(*) ! dummy

      real(8) :: tfsrc
      integer :: ic, ig, iflx, isrc, i

!::clear common variables (N0, E0, V0)
      tN0(0:ndmc,1:ndgs) = 0.0d0
      tT0(0:ndmc,1:ndgs) = 0.0d0
      tV0(0:ndmc,1:ndgs) = 0.0d0

      tNg(0:ndmc,1:2) = 0.0d0
      tTg(0:ndmc,1:2) = 0.0d0

!::sum up enutral density
      do iflx = 1, mfmax
        isrc  = visrc(iflx)
        cstyp = mcflx(iflx)
        tfsrc = snflx(iflx)
        if( tfsrc == 0.0d0 ) cycle
        if( isrc <= 0 ) cycle
        call setvmwork( 1, isrc )
        call ntwcon

        if( iflx == 8 ) then
          call ntdens(0.0d0)   ! vli  H+ + e => H0
        else
          call ntdens(tfsrc)
        endif

        do ic = 1, ncmax2  ! include vac. region
          do ig = 1, nion
            tN0(ic,ig) = tN0(ic,ig) + den0(ic,ig)
            tT0(ic,ig) = tT0(ic,ig) + den0(ic,ig)*eng0(ic,ig)
            tV0(ic,ig) = tV0(ic,ig) + den0(ic,ig)*vlp0(ic,ig)
          enddo
          tNg(ic,1) = tNg(ic,1) + deng(ic,1)
          tTg(ic,1) = tTg(ic,1) + deng(ic,1)*engg(ic,1)
          tNg(ic,2) = tNg(ic,2) + deng(ic,2)
          tTg(ic,2) = tTg(ic,2) + deng(ic,2)*engg(ic,2)
        enddo   ! loop(ic)
      enddo     ! loop(iflx)

!::temperature
      do ic = 1, ncmax2
        do ig = 1, nion
          if(tN0(ic,ig) > 0.0d0) then
            tT0(ic,ig) = tT0(ic,ig)/tN0(ic,ig)
            tV0(ic,ig) = tV0(ic,ig)/tN0(ic,ig)
          endif
        enddo
        if(tNg(ic,1) > 0.0d0) tTg(ic,1) = tTg(ic,1)/tNg(ic,1)
        if(tNg(ic,2) > 0.0d0) tTg(ic,2) = tTg(ic,2)/tNg(ic,2)
      enddo
!
      return
      end
