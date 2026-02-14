!***********************************************************************
      subroutine imgpuff
!***********************************************************************
      use cimcom, only : arbz, csbz, flbz, gfbz, icbz, iebz, ikbz, isbz
     >    , iwbz, nbz, ndbz, prbz, sflux, snbz, tfbz
      use cimpuf, only : pf_are, pf_deg, pf_flx, pf_iiw, pf_imx, pf_ity
     >    , pf_mag, pf_nty, use_gaspufftime_IMP
      use cntcom, only : csps, icps, ikps, ipps, snps, xpnt, ypnt
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::local variables
      integer  ii, i, m, iw, ic
      real*8   zsmn, zfn, zsum, zmax
      real*8   x0, y0, x1, x2, y1, y2
      real*8   dgbz(ndbz)    !  degree
! function
      integer  imox, imoy
!
      write(n6,'(/2x,"*** imgpuff ***  pf_nty =",i2)') pf_nty
!
      if( pf_nty.le.0 ) then
        call wexit("imgpuff","pf_nty = 0")
      endif

!:: time depending impurity gaspuff rate
      if(use_gaspufftime_IMP) then
        call gaspuff_update_IMP
      endif
!
      ii = 0
      zmax = -1.0d20
      zsmn = 0.0d0
      do i = 1, pf_imx
        m    = pf_ity(i)
        iw   = pf_iiw(i)
        zfn  = pf_flx(i)*pf_mag(m)
        zsmn = zsmn + zfn
        if( zfn.le.0.0d0 ) cycle
!
        zmax = dmax1( zmax, zfn )
        ii = ii + 1
        if( ii.gt.ndbz ) call wexit("imgphys","ii.gt.ndbz")
        prbz(ii) = zsmn
        flbz(ii) = zfn
        arbz(ii) = pf_are(i)
        icbz(ii) = icps(iw)
        ikbz(ii) = ikps(iw)
        isbz(ii) = ipps(iw)
        iebz(ii) = ipps(iw+1)
        iwbz(ii) = iw
        csbz(ii) = csps(iw)
        snbz(ii) = snps(iw)
        dgbz(ii) = pf_deg(i)
      enddo
!
!::return
      nbz = ii
      tfbz = zsmn
      gfbz = zmax
      if(nbz .eq. 0) then
        sflux = 0.0d0
        prbz  = 0.0d0
        return
      else
!
!::normalization
        do i = 1, nbz
          prbz(i) = prbz(i)/zsmn
        enddo
        prbz(nbz) = 1.0001d0
!
!::output
        write(n6,'(4x,"i",4x,"iw",4x,"x0",7x,"y0",7x,"th",7x,"prb",9x,
     >    "spt",9x,"are",9x,"ix",3x,"iy")')
        zsum = 0.0d0
        do i = 1, nbz
          iw = iwbz(i)
          ic = icbz(i)
          if(use_exdata) then
            call get_position_byic(ic,isbz(i),.true.,x1)
            call get_position_byic(ic,iebz(i),.true.,x2)
            call get_position_byic(ic,isbz(i),.false.,y1)
            call get_position_byic(ic,iebz(i),.false.,y2)
          else
            x1 = xpnt(isbz(i))
            x2 = xpnt(iebz(i))
            y1 = ypnt(isbz(i))
            y2 = ypnt(iebz(i))
          endif
          x0 = 0.5d0*(x1+x2)
          y0 = 0.5d0*(y1+y2)
          zsum = zsum + flbz(i)
          write(n6,'(2x,2i5,2f9.4,f9.3,1p3e12.3,2i5)') i,iw,x0,y0
     >     , dgbz(i),prbz(i), flbz(i)/arbz(i)
     >     , arbz(i), imox(ic), imoy(ic)
        enddo
!
        sflux = zsum
        write(n6,'(2x,10x,18x,2x,"total",2x,1p2e12.3)') tfbz, sflux
      endif
!
      return
!
!**********************************************************************
      contains
      subroutine gaspuff_update_IMP
!**********************************************************************
      use cimpuf, only: ndpfz, pf_mag
      implicit none
!::local variables
      integer m
      real*8 interpolate_value
      do m = 1, ndpfz
        call timevalue_inter_IMP(interpolate_value,m)
        pf_mag(m) = interpolate_value
      enddo

      end subroutine gaspuff_update_IMP

!**********************************************************************
      subroutine timevalue_inter_IMP(output_val,target_m)
!**********************************************************************
      use csonic, only: time
      use cimpuf, only: pufftime_size_IMP, pufftime_IMP, pf_mag_time
      use cimctl, only: nprg
      use cplcom, only: qtim_ini
      implicit none
!
!::arguments
      integer, intent(in) :: target_m
      real*8, intent(out) :: output_val
!
!::local variables
      integer i_time
      real*8 t_next,t_preb,val_next,val_preb,time_accu

      if(time.ne.0.0d0) then
        time_accu = time
      else
        time_accu = qtim_ini ! first step for continue calculation
      endif

      ! time less than time-profile
      if(time_accu .le. pufftime_IMP(1)) then
        output_val = pf_mag_time(target_m,1,nprg)
        return
      endif
      ! time over than time-profile
      if(time_accu .ge. pufftime_IMP(pufftime_size_IMP)) then
        output_val = pf_mag_time(target_m,pufftime_size_IMP,nprg)
        return
      endif

      do i_time = 2, pufftime_size_IMP
        t_next = pufftime_IMP(i_time)
        if(time_accu .le. t_next) then
          t_preb = pufftime_IMP(i_time-1)
          val_next = pf_mag_time(target_m,i_time,nprg)
          val_preb = pf_mag_time(target_m,i_time-1,nprg)
          ! liner interpolate of value
          output_val = val_preb +
     >      (val_next-val_preb)/(t_next-t_preb) * (time_accu-t_preb)
          return
        endif
      enddo
      end subroutine timevalue_inter_IMP
!::
      end subroutine imgpuff
