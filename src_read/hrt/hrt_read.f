!-------------------------------------------------------------
      subroutine hrt_read()
!-------------------------------------------------------------
!
!     read H Radiation Trapping data from Sawada-sensei
!   
!               K.Hoshino 2015/8/12

      use chrt, only : cdhrt, facn0, facne, facng, fact0, facte, lfn0
     >    , lfne, lfng, lft0, lfte, mxn0, mxne, mxng, mxt0, mxte, sgalf
     >    , sgion, sgrec, sn0, sn0max_hrt, sn0min_hrt, sne, snemax_hrt
     >    , snemin_hrt, sng, sngmax_hrt, sngmin_hrt, st0, st0max_hrt
     >    , st0min_hrt, ste, stemax_hrt, stemin_hrt
      implicit none
!
      integer kne,kte, kn0, kt0, kng
      integer i1,i2,i3,i4,i5
      real*8  tmp1,tmp2, tmp3

      lfte = 9
      lfne = 5
      lft0 = 5
      lfn0 = 5
      lfng = 5

      mxte = 3 * lfte
      mxne = 3 * lfne
      mxt0 = 2 * lft0
      mxn0 = 3 * lfn0
      mxng = 3 * lfng

      facte = exp(log(10.0d0)/dble(lfte-1))
      facne = exp(log(10.0d0)/dble(lfne-1))
      fact0 = exp(log(10.0d0)/dble(lft0-1))
      facn0 = exp(log(10.0d0)/dble(lfn0-1))
      facng = exp(log(10.0d0)/dble(lfng-1))

      stemin_hrt = 0.1d0  
      snemin_hrt = 1.0d20 
      st0min_hrt = 0.1d0  
      sn0min_hrt = 1.0d18
      sngmin_hrt = 1.0d18

      stemax_hrt = stemin_hrt * facte**dble(mxte-1)
      snemax_hrt = snemin_hrt * facne**dble(mxne-1)
      st0max_hrt = st0min_hrt * fact0**dble(mxt0-1)
      sn0max_hrt = sn0min_hrt * facn0**dble(mxn0-1)
      sngmax_hrt = sngmin_hrt * facng**dble(mxng-1)

      do i1=1,mxte
        ste(i1) = dlog(stemin_hrt * facte**dble(i1-1))
      enddo
      do i1=1,mxne
        sne(i1) = dlog(snemin_hrt * facne**dble(i1-1))
      enddo
      do i1=1,mxt0
        st0(i1) = dlog(st0min_hrt * fact0**dble(i1-1))
      enddo
      do i1=1,mxn0
        sn0(i1) = dlog(sn0min_hrt * facn0**dble(i1-1))
      enddo
      do i1=1,mxng
        sng(i1) = dlog(sngmin_hrt * facng**dble(i1-1))
      enddo

      open(unit=100, file=cdhrt)

      do i1=1,mxte
      do i2=1,mxne
      do i3=1,mxt0
      do i4=1,mxn0
      do i5=1,mxng
        read(100,*)kte,kne,kt0,kn0,kng,tmp1, tmp2, tmp3
        sgrec(kte,kne,kt0,kn0,kng) = dmax1(tmp1,1.0d-40)
        sgion(kte,kne,kt0,kn0,kng) = dmax1(tmp2,1.0d-40)
        sgalf(kte,kne,kt0,kn0,kng) = dmax1(tmp3,1.0d-40)
      enddo
      enddo
      enddo
      enddo
      enddo

      close(100)

      return

      end subroutine
