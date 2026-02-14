!***********************************************************************
      subroutine plstep
!***********************************************************************
      use cntctl,     only : itmnt, nclnt
      use cplcom,     only : dlpq, dtlmt, dtmax, dtmin, dtnxt, edtmx
     >    , edtq, elpmx, itfix, nion, q1a, q2a, q3, q4, wq1a, wq2a, wq3
     >    , wq4
      use csize,      only : ndx, ndy
      use csonic,     only : dtim, itim, limp, lstop, time
      use cunit,      only : n6
      use topics_mod, only : dtcal, dtduc, dtim_s, lgtpc
      implicit none

!::local variables
      integer ia, j, i
      integer kout, kcnv
      real*8  dtstp, tfac, sfac, dtim0, fclp, fcdt, ftmp
!
      real*8  fcdtmx
!
!::input data
      fcdtmx = 1.8d0
!
!-----------------------------------------------------------------------
!::store q(N)
!-----------------------------------------------------------------------
      do ia = 1, nion
      do j  = 1, ndx
      do i  = 1, ndy
      wq1a(j,i,ia) = q1a(j,i,ia)
      wq2a(j,i,ia) = q2a(j,i,ia)
      enddo
      enddo
      enddo
!
      do j = 1, ndx
      do i = 1, ndy
      wq3(j,i) = q3(j,i)
      wq4(j,i) = q4(j,i)
      enddo
      enddo
!
!-----------------------------------------------------------------------
!::solve q(N+1,l+1)
!-----------------------------------------------------------------------
      do ! loop
        kout = 2
        dtstp = dtim
        call pwstep(kout,dtstp,kcnv)
!
!::converge
        if( kcnv.eq.0 ) exit
!
!::fixed dtim
        if( kcnv.eq.2 .and. itim.gt.itfix ) then
          dtim0 = dtim
          tfac = 1.0d0
          fclp = elpmx/dlpq
          fcdt = edtmx/edtq
          write(n6,'(a,i6,i4,2x,i8,"  dtim =",1pe10.3,
     >    "  _nx =",1pe11.3,"  tfac =",0pf6.3,
     >    "  fclp,fcdt =",1p2e10.2)')
     >    "<dtcntl-pl> dt itfix ",itim,nclnt,itmnt,dtim0,dtim,
     >    tfac,fclp,fcdt
          exit
        endif
!
!-----------------------------------------------------------------------
!::recalculation with new dtim due to negative value or large variation
!-----------------------------------------------------------------------
!::new time-step
        if( kcnv.eq.1 ) then        ! iteration
          tfac = 0.80d0
        elseif( kcnv.eq.2 ) then    ! dt-step
          tfac = edtq/edtmx*1.2d0
          tfac = 1.0d0/tfac
          tfac = dmax1(tfac,0.70d0)
        elseif( kcnv.eq.3 ) then    ! negative value
          lstop = 0
          tfac = 0.60d0
        endif
!----
        dtim0 = dtim
        dtim_s = dtim * tfac
        if( lgtpc == 1 .and. dtim_s > dtduc ) then
          if( dtduc > 0.0_8 ) then
            dtim   = dtduc
          else
            dtim   = min( dtim_s, dtcal )
          endif
        else
          dtim  = dtim_s
        endif

        fclp = elpmx/dlpq
        fcdt = edtmx/edtq
        write(n6,'(a,i6,i4,2x,i8,"  dtim =",1pe10.3,"  _nx =",1pe11.3,
     >  "  tfac =",0pf6.3,"  fclp,fcdt =",1p2e10.2)')
     >  "<dtcntl-pl> small dt ",itim,nclnt,itmnt,dtim0,dtim,
     >  tfac,fclp,fcdt
!----
!
!::too small dtim
        if( dtim.lt.dtmin ) then
          write(n6,'(/2x,"*** serious error at sub, plstep ***"
     >     /5x,"dtim is too small !  ",1p2e12.3)') dtim,dtmin
          lstop = 1
          goto 900
        endif
!
!::q(N+1,l=0) = q(N)
        do ia = 1, nion
        do j  = 1, ndx
        do i  = 1, ndy
        q1a(j,i,ia) = wq1a(j,i,ia)
        q2a(j,i,ia) = wq2a(j,i,ia)
        enddo
        enddo
        enddo
!
        do j = 1, ndx
        do i = 1, ndy
        q3(j,i) = wq3(j,i)
        q4(j,i) = wq4(j,i)
        enddo
        enddo
!
!::aux-variables
        call plauxv
!
      enddo
!
!-----------------------------------------------------------------------
!::dtim for next step (dtnxt)
!-----------------------------------------------------------------------
      call plpcon_out     !  hisotry (zflx)
      call plout_reg      !  history (zreg)
!
!::KH170202 particle flux to wall for C-genaration in IMPMC
      call ntwflx
      call plwflx
!
!!!old-version of aasonic
      time  = time + dtim
      dtim0 = dtim
      dtnxt = dtim
!
      fclp = elpmx/dlpq
      fcdt = edtmx/edtq
!
!::time step control (OLD type)
      if( fcdt.ge.fcdtmx ) then
      tfac = 1.2d0+0.2d0*dmax1(fcdt-3.0d0,0.0d0)
      tfac = dmin1( tfac, 2.0d0 )
      dtnxt = dtnxt*tfac
      dtnxt = dmin1( dtnxt, dtmax )
      endif

      tfac = dtnxt/dtim0
      write(n6,'(a,i6,i4,2x,i8,"  dtim =",1pe10.3,"  _nx =",1pe11.3,
     >  "  tfac =",0pf6.3,"  fclp,fcdt =",1p2e10.2)')
     >  "<dtcntl-P3> large dt ",itim,nclnt,itmnt,dtim0,dtnxt,
     >  tfac,fclp,fcdt
!
      return
!
!-----------------------------------------------------------------------
!::serious error
!-----------------------------------------------------------------------
 900  continue
      write(n6,'(/2x,"*** plstep ***  set old value due to lstop =",
     >   i2)')  lstop
!
!::q(N)
      do ia = 1, nion
      do j  = 1, ndx
      do i  = 1, ndy
      q1a(j,i,ia) = wq1a(j,i,ia)
      q2a(j,i,ia) = wq2a(j,i,ia)
      enddo
      enddo
      enddo
!
      do j = 1, ndx
      do i = 1, ndy
      q3(j,i) = wq3(j,i)
      q4(j,i) = wq4(j,i)
      enddo
      enddo
!
      call plauxv
!
      return
      end
