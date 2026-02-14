!***********************************************************************
      subroutine plqcon_cl
!***********************************************************************
      use cplcom, only : nion
      use cplqcn, only : qfx_cd, qfx_cv, qfx_df, qfx_vh, qfy_cd, qfy_df
     >    , qfy_vh, qvl_al, qvl_cl, qvl_dt, qvl_pe, qvl_pi, qvl_sc
      use cplvpn, only : qfy_vp
      use csize,  only : ndx, ndy
      implicit none
!
      integer  k1, k2, k3, nequ
!
      nequ = 2*nion + 2
      do k3=1,nequ
      do k2=1,ndy
      do k1=1,ndx
      qfx_cv(k1,k2,k3) = 0.0d0
      qfx_df(k1,k2,k3) = 0.0d0
      qfy_df(k1,k2,k3) = 0.0d0
      qfx_cd(k1,k2,k3) = 0.0d0
      qfy_cd(k1,k2,k3) = 0.0d0
      qvl_pe(k1,k2,k3) = 0.0d0
      qvl_cl(k1,k2,k3) = 0.0d0
      qvl_sc(k1,k2,k3) = 0.0d0
      qvl_al(k1,k2,k3) = 0.0d0
      qvl_dt(k1,k2,k3) = 0.0d0
!
      qfy_vp(k1,k2,k3) = 0.0d0
      enddo
      enddo
      enddo
!
      do k3=1,nion
      do k2=1,ndy
      do k1=1,ndx
      qvl_pi(k1,k2,k3) = 0.0d0
      enddo
      enddo
      enddo
!
      do k2=1,ndy
      do k1=1,ndx
      qfx_vh(k1,k2) = 0.0d0
      qfy_vh(k1,k2) = 0.0d0
      enddo
      enddo
!
      return
      end
!
!***********************************************************************
      subroutine plqcon_sx(jst,jen,ist,ien,iout)
!***********************************************************************
!
!   integ flux against x-direction from jst to jen
!
!    main region   call plqcon_sx(jcxp1+1,jcxp2-1,icspx,icmpe,iout)
!    odiv region   call plqcon_sx( 2,     jcxp1,  icwl1,icwl2,iout)
!    sol  region   call plqcon_sx( 2,     jcmax-1,icwl1,icspx,iout)
!    oprv region   call plqcon_sx( 2,     jcxp1,  icspx+1, icwl2, iout)
!
!
!   cell name in metric calculation
!
!                       +icaxs--------------------+
!                       !                         !
!                       !icmpe~~~~~~~icmax(jc)~~~~!
!      icwl2 /+^^^^^^^^^+                         +^^^^^^^^^+/
!            /!  (4)    !icmps  Main plasma (6)   !    (5)  !/
!            /!---------!-------------------------!---------!/
!      icspx /!  OUT    !                         !  IN     !/
!            /!   DIV   !      Scrape-Off         !   DIV   !/
!            /!  (1)    !                  (2)    !    (3)  !/
!   A  icwl1 /j+^^^^^^^j^^^^^^^^^^^^^icmin(jc)^^^^^j^^^^^^^^j j
! i !         c        c                           c        c c
! c !         d        x                           x        d m
!   !-->    1=p        p                           p        p=a
!     jc      1        1                           2        2 x
!
!                do jc = 1, jcmax                ir = kreg(jc,ic)
!                do ic = icmin(jc), icmax(jc)
!
!----------------------------------------------------------------------
      use cplcom, only : ama, nion
      use cplqcn, only : qsx_cd, qsx_df, qsx_vh, qfy_cd, qfy_df, qfy_vh
      use cplvpn, only : qfy_vp, qsx_vp
      use cunit,  only : n6
      use csize,  only : ndeq
      implicit none
!
!::argument
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer  jst,jen,ist,ien,iout
      integer, intent(in) :: jst, jen, ist, ien, iout
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer ic, jc, i, nsiz
      integer ic, jc
      integer ia, m1a, m2a, m3, m4
!
      m3 = 2*nion + 1
      m4 = 2*nion + 2
!
      do ic = ist, ien
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   i = ic
!
! modified 4/3 lines organize local variables and include files by kamata 2021/05/31
!ik   nsiz = nqfm
!ik   call setd( qsx_df, nsiz, 0.0d0 )
!ik   call setd( qsx_cd, nsiz, 0.0d0 )
!ik   call setd( qsx_vp, nsiz, 0.0d0 )
      call setd( qsx_df, ndeq, 0.0d0 )
      call setd( qsx_cd, ndeq, 0.0d0 )
      call setd( qsx_vp, ndeq, 0.0d0 )
      qsx_vh = 0.0d0
!
      do jc = jst, jen
!
      do ia = 1, nion
      m1a = 2*ia - 1
      m2a = 2*ia
      qsx_df(m1a) = qsx_df(m1a) + qfy_df(jc,ic,m1a)/ama(ia)
      qsx_vp(m1a) = qsx_vp(m1a) + qfy_vp(jc,ic,m1a)/ama(ia)
      qsx_vp(m2a) = qsx_vp(m2a) + qfy_vp(jc,ic,m1a)
      enddo
      qsx_df(m3)  = qsx_df(m3)  + qfy_df(jc,ic,m3)
      qsx_df(m4)  = qsx_df(m4)  + qfy_df(jc,ic,m4)
      qsx_vp(m3)  = qsx_vp(m3)  + qfy_vp(jc,ic,m3)
      qsx_vp(m4)  = qsx_vp(m4)  + qfy_vp(jc,ic,m4)
      qsx_cd(m3)  = qsx_cd(m3)  + qfy_cd(jc,ic,m3)
      qsx_cd(m4)  = qsx_cd(m4)  + qfy_cd(jc,ic,m4)
      qsx_vh      = qsx_vh      + qfy_vh(jc,ic)
!
      enddo  ! loop(jc)
!
!::debug write
      if( iout.eq.1 ) then
      ia = 1
      m1a = 2*ia - 1
      write(n6,'(2x,"[plqfly]  i =",i3,2x,i3," =>",i3,1pe12.3,
     >  2x,1p4e12.3,2x,1p4e12.3,2x,1p3e12.3)')
     >  ic, jst, jen, qsx_df(m1a),
     >  qsx_df(m3),qsx_cd(m3),qsx_vh, qsx_df(m3)+qsx_cd(m3)+qsx_vh,
     >  qsx_df(m4),qsx_cd(m4),0.0d0,  qsx_df(m4)+qsx_cd(m4),
     >  qsx_vp(m1a), qsx_vp(m3), qsx_vp(m4)
      endif
!
      enddo  !  loop(ia)
!
      return
      end
