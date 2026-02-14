!***********************************************************************
      subroutine plqcon_edge
!***********************************************************************
!
!     dq/dt + div(F) + div(G) + Hp = S
!
!      V*dq/dt = - F^psi(i+1/2) + F^psi(i-1/2) - V*Hp + V*S
!
!             j-integrate in the core edge
!
!       qvl_dt(j,i,m1a) = (q1a(j,i,ia)-wq1a(j,i,ia))/dtvl1
!       qvl_dt(j,i,m2a) = (q2a(j,i,ia)-wq2a(j,i,ia))/dtvl2
!       qvl_dt(j,i,m3)  = (q3(j,i)-wq3(j,i))/dtvl3
!       qvl_dt(j,i,m4)  = (q4(j,i)-wq4(j,i))/dtvl4
!
!-----------------------------------------------------------------------
      use cplcom, only : ama, nion, q1a, q2a, q3, q4, wq1a, wq2a, wq3
     >    , wq4
      use cplmet, only : hvol, icmpe, itmpe, itmps, icspx, jcxp1, jcxp2
      use cplqcn, only : qfy_cd, qfy_df, qvl_cl, qvl_dt, qvl_pe, qvl_sc
      use csonic, only : dtim, itim, time
      implicit none
!
!::local variables
      integer  jst, jen, k, ift, i, j, it
      integer  ia, m1a, m2a, m3, m4, mm
      real*8   zvy_df, zvy_cd, zvl_pe, zvl_cl, zvl_sc, zvl_dt
      real*8   zvym_df, zvym_cd, zsm_q, zsm_wq, zsm_dt, w
!
      return   ! Not passed plqcon_edge
!xx   if( mod(itim,100).ne.10 ) return
!
!::index
      jst = jcxp1+1
      jen = jcxp2-1
      ia = 1
      m1a = 2*ia - 1
      m2a = 2*ia
      m3  = 2*nion + 1
      m4  = 2*nion + 2
!
!::q1a, q2a, q3, q4
      do k = 1, 4
      if( k.eq.1 ) mm = m1a
      if( k.eq.2 ) mm = m2a
      if( k.eq.3 ) mm = m3
      if( k.eq.4 ) mm = m4
      ift = 87000 + k
      write(ift,'(/4x,"itim",3x,"time",10x,"dtim",8x,2x,"i",2x,
     >   "vym_df",6x,"vym_cd",6x,"vy_df",7x,"vy_cd",7x,"vl_pe",7x,
     >   "vl_cl",7x,"vl_sc",7x,"vl_dt",7x,"sm_q",8x,"dsm_q")')
!
      zvy_df = 0.0d0
      zvy_cd = 0.0d0
      zvl_pe = 0.0d0
      zvl_cl = 0.0d0
      zvl_sc = 0.0d0
      zvl_dt = 0.0d0
!
      do i = icspx, icmpe-1
!
      zvym_df = zvy_df
      zvym_cd = zvy_cd
!
      zvy_df = 0.0d0
      zvy_cd = 0.0d0
      zvl_pe = 0.0d0
      zvl_cl = 0.0d0
      zvl_sc = 0.0d0
      zvl_dt = 0.0d0
!
      zsm_q  = 0.0d0
      zsm_wq = 0.0d0
!
      do j = jst, jen
      zvy_df = zvy_df + qfy_df(j,i,mm)
      zvy_cd = zvy_cd + qfy_cd(j,i,mm)
      zvl_pe = zvl_pe + qvl_pe(j,i,mm)
      zvl_cl = zvl_cl + qvl_cl(j,i,mm)
      zvl_sc = zvl_sc + qvl_sc(j,i,mm)
      zvl_dt = zvl_dt + qvl_dt(j,i,mm)
!
      if( k.eq.1 ) then
        zsm_q  = zsm_q  +  q1a(j,i,ia)*hvol(j,i)
        zsm_wq = zsm_wq + wq1a(j,i,ia)*hvol(j,i)
      elseif( k.eq.2 ) then
        zsm_q  = zsm_q  +  q2a(j,i,ia)*hvol(j,i)
        zsm_wq = zsm_wq + wq2a(j,i,ia)*hvol(j,i)
      elseif( k.eq.3 ) then
        zsm_q  = zsm_q  +  q3(j,i)*hvol(j,i)
        zsm_wq = zsm_wq + wq3(j,i)*hvol(j,i)
      elseif( k.eq.4 ) then
        zsm_q  = zsm_q  +  q4(j,i)*hvol(j,i)
        zsm_wq = zsm_wq + wq4(j,i)*hvol(j,i)
      endif
      enddo
!
      zsm_dt = - (zvy_df+zvy_cd) + (zvym_df+zvym_cd)
     >         - zvl_pe - zvl_cl + zvl_sc
!
      w = 1.0d0
      if( k.eq.1 ) w = 1.0d0/ama(ia)
!
      write(ift,'(2x,i7,1pe14.6,1pe12.3,2x,i3,1p11e12.3)')
     >  itim, time, dtim, i, zvym_df*w, zvym_cd*w, zvy_df*w,
     >  zvy_cd*w, zvl_pe*w, zvl_cl*w, zvl_sc*w,
     >  zvl_dt*w, zsm_dt*w,
     >  zsm_q*w, (zsm_q-zsm_wq)/dtim*w
      enddo
      enddo
!
!::plqcon
      do it = itmps, itmpe
      call plqcon_flx(it)
      enddo
!
      return
      end
