!------------------------------------------------------------
      real*8 function hrt_eval(ctyp, zte, zne, zt0, zn0, zng)
!------------------------------------------------------------
!
!     liner interpolat of  H Radiation Trapping data
!
!               K.Hoshino 2015/8/12
!
      use chrt,  only : mxn0, mxne, mxng, mxt0, mxte, sgalf, sgion
     >    , sgrec, sn0, sn0max_hrt, sn0min_hrt, sne, snemax_hrt
     >    , snemin_hrt, sng, sngmax_hrt, sngmin_hrt, st0, st0max_hrt
     >    , ste, stemax_hrt, stemin_hrt
      use cunit, only : n6
      implicit none
!
      character*3, intent(in) :: ctyp
      real(8),     intent(in) :: zte, zne, zt0, zn0, zng

!
      real(8) zzte, zzne, zzt0, zzn0, zzng
      integer kte, kne, kt0, kn0, kng
      real*8 ate(0:2), ane(0:2), at0(0:2), an0(0:2), ang(0:2)
      real*8 ss(2,2,2,2,2)

      real*8 zlte, zlne, zlt0, zln0, zlng
!      real*8 zite, zine, zit0, zin0, zing
      integer ite, ine, it0, in0, ing, i
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik   real*8 tmp
!
! finction
      real(8)    hrt_fit
!
      zzte = zte
      zzne = zne
      zzt0 = zt0
      zzn0 = zn0
      zzng = zng
!:: index te
      zzte = dmax1( zzte, stemin_hrt )
      zzte = dmin1( zzte, stemax_hrt )
      zlte = dlog(zzte)
!      zite = (zlte-stemin)/stedlt + 1.0d0
!      kte  = ding(zite)
      do i = 1, mxte-1
        if(zlte .le. ste(i+1))then
          kte = i
          ate(0) = zlte
          ate(1) = ste(i)
          ate(2) = ste(i+1)
          exit
        endif
      enddo

!:: index ne
      zzne = dmax1( zzne, snemin_hrt )
      zzne = dmin1( zzne, snemax_hrt )
      zlne = dlog(zzne)
!      zine = (zlne-snemin)/snedlt + 1.0d0
!      kne  = ding(zine)
      do i = 1, mxne-1
        if(zlne .le. sne(i+1))then
          kne = i
          ane(0) = zlne
          ane(1) = sne(i)
          ane(2) = sne(i+1)
          exit
        endif
      enddo

!:: index t0
!      zt0 = dmax1( zt0, st0min )
!      if( zt0 .le. 0.5d0 )then
!        write(n6,*)"HRT: zt0=",zn0
!      endif
      zzt0 = dmax1( zzt0, 0.5d0 )
      zzt0 = dmin1( zzt0, st0max_hrt )
      zlt0 = dlog(zzt0)
!      zit0 = (zlte-st0min)/st0dlt + 1.0d0
!      kt0  = ding(zit0)
      do i = 1, mxt0-1
        if(zlt0 .le. st0(i+1) )then
          kt0 = i
          at0(0) = zlt0
          at0(1) = st0(i)
          at0(2) = st0(i+1)
          exit
        endif
      enddo

!:: index n0
      zzn0 = dmax1( zzn0, sn0min_hrt )
!      if( zn0 .ge. sn0max )then
!        write(n6,*)"HRT: zn0=",zn0
!      endif
      zzn0 = dmin1( zzn0, sn0max_hrt )
      zln0 = dlog(zzn0)
!      zin0 = (zlte-sn0min)/sn0dlt + 1.0d0
!      kn0  (zin0)
      do i = 1, mxn0-1
        if(zln0 .le. sn0(i+1) )then
          kn0 = i
          an0(0) = zln0
          an0(1) = sn0(i)
          an0(2) = sn0(i+1)
          exit
        endif
      enddo

!:: index ng
      zzng = dmax1( zzng, sngmin_hrt )
      zzng = dmin1( zzng, sngmax_hrt )
      zlng = dlog(zzng)
!      zing = (zlng-sngmin)/sngdlt + 1.0d0
!      kng  = ding(zing)
      do i = 1, mxng-1
        if(zlng .le. sng(i+1) )then
          kng = i
          ang(0) = zlng
          ang(1) = sng(i)
          ang(2) = sng(i+1)
          exit
        endif
      enddo


      if(ctyp(1:3).eq."ion")then
        do ing = 1,2
        do in0 = 1,2
        do it0 = 1,2
        do ine = 1,2
        do ite = 1,2
          ss(ite, ine, it0, in0, ing) =
     >        dlog(sgion(kte+ite-1, kne+ine-1,
     >                   kt0+it0-1, kn0+in0-1, kng+ing-1))
        enddo
        enddo
        enddo
        enddo
        enddo
      elseif(ctyp(1:3).eq."alp")then
        do ing = 1,2
        do in0 = 1,2
        do it0 = 1,2
        do ine = 1,2
        do ite = 1,2
          ss(ite, ine, it0, in0, ing) =
     >        dlog(sgalf(kte+ite-1, kne+ine-1,
     >                   kt0+it0-1, kn0+in0-1, kng+ing-1))
        enddo
        enddo
        enddo
        enddo
        enddo
      elseif(ctyp(1:3).eq."rec")then
        do ing = 1,2
        do in0 = 1,2
        do it0 = 1,2
        do ine = 1,2
        do ite = 1,2
          ss(ite, ine, it0, in0, ing) =
     >        dlog(sgrec(kte+ite-1, kne+ine-1,
     >                   kt0+it0-1, kn0+in0-1, kng+ing-1))
        enddo
        enddo
        enddo
        enddo
        enddo
      else
        write(n6,*) "HRT: type error: ",ctyp
        return
      endif

      hrt_eval = exp(hrt_fit(ss,ate,ane,at0,an0,ang))

      end function
!
!
!
!------------------------------------------------------------
      real*8 function hrt_fit(ss,ate,ane,at0,an0,ang)
!------------------------------------------------------------
!
!  interporation of H radiation trapping data
!
      implicit none
!
! modified 2/3 lines organize local variables and include files by kamata 2021/06/16
!ik   real*8 ss(2,2,2,2,2)
!ik   real*8 ate(0:2), ane(0:2), at0(0:2), an0(0:2), ang(0:2)
      real(8), intent(in) :: ss(2,2,2,2,2)
      real(8), intent(in) :: ate(0:2), ane(0:2), at0(0:2), an0(0:2)
     >                     , ang(0:2)

! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer kte, kne, kt0, kn0, kng
      integer kne, kt0, kn0, kng
      real*8 zz(3,3,3,3)

      zz=0.0d0

      do kng = 1,2
      do kn0 = 1,2
      do kt0 = 1,2
      do kne = 1,2
        zz(kne,kt0,kn0,kng) =
     >     (ss(1,kne,kt0,kn0,kng)-ss(2,kne,kt0,kn0,kng))
     >     /(ate(1) - ate(2)) * (ate(0)-ate(1)) + ss(1,kne,kt0,kn0,kng)
      enddo !kne
      enddo !kt0
      enddo !kn0
      enddo !kng

      do kng = 1,2
      do kn0 = 1,2
      do kt0 = 1,2
        zz(3,kt0,kn0,kng) =
     >     (zz(1,kt0,kn0,kng)-zz(2,kt0,kn0,kng))
     >     /(ane(1) - ane(2)) * (ane(0)-ane(1)) + zz(1,kt0,kn0,kng)
      enddo !kt0
      enddo !kn0
      enddo !kng

      do kng = 1,2
      do kn0 = 1,2
        zz(3,3,kn0,kng) =
     >     (zz(3,1,kn0,kng)-zz(3,2,kn0,kng))
     >     /(at0(1) - at0(2)) * (at0(0)-at0(1)) + zz(3,1,kn0,kng)
      enddo !kn0
      enddo !kng

      do kng = 1,2
        zz(3,3,3,kng) =
     >     (zz(3,3,1,kng)-zz(3,3,2,kng))
     >     /(an0(1) - an0(2)) * (an0(0)-an0(1)) + zz(3,3,1,kng)
      enddo !kng

      zz(3,3,3,3) =
     >     (zz(3,3,3,1)-zz(3,3,3,2))
     >     /(ang(1) - ang(2)) * (ang(0)-ang(1)) + zz(3,3,3,1)

      hrt_fit = zz(3,3,3,3)

      return
      end
