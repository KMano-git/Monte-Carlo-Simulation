!***********************************************************************
      module mod_soldorTimeSeries
      implicit none
!***********************************************************************
      integer timeNum ! = 1(defoult). Number of time steps for time series
      integer interNum !  = 1(defoult). Time step interval to time series sampling
      integer,parameter :: dataNum = 10 ! Amount of save time series 
      integer,parameter :: icNum = 100 ! icNum > itpve. should be allocatable
      integer :: seriesNow = 0 ! Current amount of data in timeSeries_i,o
      integer :: interNow = 0 ! Time steps since the last time the results were saved

      real*8,dimension(dataNum,icNum) :: latestVal_o ! for Qdpl_o.txt
      real*8,dimension(dataNum,icNum) :: latestVal_i ! for Qdpl_i.txt

!    dimension(timeNum,dataNum,icNum). The smaller the first argument, the older result.
      real*8,allocatable,dimension(:,:,:)
     >  :: timeSeries_o 
      real*8,allocatable,dimension(:,:,:)
     >  :: timeSeries_i 
!***********************************************************************
      contains
      subroutine calcLatestVal()
!***********************************************************************
!     simmilar to pst_qdpl.f
      use cgdcom, only : grdx, grdy
      use cphcns, only : cev, cpi
      use cplcom, only : ama, nion, vcs, vne, vte, vti, vva, vni
      use cplmet, only : hvsb, icel, itpve, itsls, jcel, jtmax
      use cplpst, only : dwrad, farei, fareo, fbdpi, fbdpo, fldpi, fldpo
     >    , flsmi, flsmo, tflsmi, tflsmo
      use cplqcn, only : qfx_cd, qfx_cv, qfx_vh 
      use cntwfl, only : dwntl
      use csize,  only : ndy
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer  nt1, nt2, m, it, nt, i, jmax, jw, j, jp, jm, jsf, kbc
      integer  ia, m1a, m3, m4 
      integer  k
      real*8   sigv, zcs
!
      integer  mx1, mx2, mx3, mx4, my1, my2, my3, my4
      real*8   r1, r2, r3, r4, z1, z2, z3, z4, dl
      real*8   zsv, zsp
      real*8   fsg, zflx, erec, hcsw

      m3   = 2*nion + 1
      m4   = 2*nion + 2
      nt1 = itsls + 1
      nt2 = itpve - 1

!      ::clear
      do m = 1, m4
            do it = 1, ndy
                  fbdpi(it,m) = 0.0d0
                  fbdpo(it,m) = 0.0d0
                  fldpi(it,m) = 0.0d0
                  fldpo(it,m) = 0.0d0
            enddo
      enddo      
      do i = 1, 10
            tflsmi(i) = 0.0d0
            tflsmo(i) = 0.0d0
            do it = 1, ndy
                  flsmi(it,i) = 0.0d0
                  flsmo(it,i) = 0.0d0
            enddo
      enddo

! ::loop (it)
      do nt = nt1, nt2
            it   = nt
            jmax = jtmax(it)
!
!-----------------------------------------------------------------------
!::boundary condition j = 1
!-----------------------------------------------------------------------
            jw   = 1
            i    = icel(jw,it)
            j    = jcel(jw,it)
            jp   = jcel(jw+1,it)
            jsf  = j
            sigv = -1.0d0
            zcs  = sigv*vcs(j,i)   !  for all species
            kbc  = 1
!
!::are
            call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
            r2 = grdx(mx2,my2)
            r3 = grdx(mx3,my3)
            z2 = grdy(mx2,my2)
            z3 = grdy(mx3,my3)
            dl = dsqrt((r3-r2)**2+(z3-z2)**2)
            fareo(it) = 2.0d0*cpi*0.5d0*(r2+r3)*dl
  
!::New type
            fsg = -1.0d0
            flsmo(it,1) = fsg*qfx_cv(jsf,i,m4)
            flsmo(it,2) = fsg*qfx_cd(jsf,i,m4)
            flsmo(it,3) = fsg*qfx_cv(jsf,i,m3)
            flsmo(it,4) = fsg*(qfx_cd(jsf,i,m3)+qfx_vh(jsf,i))
      
            ia  = 1
            m1a = 2*ia - 1
            zflx = qfx_cv(jsf,i,m1a)/ama(ia)
            erec = 13.595
            flsmo(it,5) = fsg*zflx*erec*cev
            flsmo(it,6) = flsmo(it,5)+dwrad(i,4)
            flsmo(it,7) = flsmo(it,6)+dwntl(i,4)

!-----------------------------------------------------------------------
!::boundary condition j = jmax
!-----------------------------------------------------------------------
            jw   = jmax
            i    = icel(jw,it)
            j    = jcel(jw,it)
            jm   = jcel(jw-1,it)
            jsf  = jcel(jw-1,it)
            sigv = +1.0d0
            zcs  = sigv*vcs(j,i)   !  for all species
            kbc  = 2
!::are
            call mcpnt(j,i,mx1,mx2,mx3,mx4,my1,my2,my3,my4)
            r4 = grdx(mx4,my4)
            r1 = grdx(mx1,my1)
            z4 = grdy(mx4,my4)
            z1 = grdy(mx1,my1)
            dl = dsqrt((r4-r1)**2+(z4-z1)**2)
            farei(it) = 2.0d0*cpi*0.5d0*(r4+r1)*dl
!
!::New type
            fsg = +1.0d0
            flsmi(it,1) = fsg*qfx_cv(jsf,i,m4)
            flsmi(it,2) = fsg*qfx_cd(jsf,i,m4)
            flsmi(it,3) = fsg*qfx_cv(jsf,i,m3)
            flsmi(it,4) = fsg*(qfx_cd(jsf,i,m3)+qfx_vh(jsf,i))

            ia  = 1
            m1a = 2*ia - 1
            zflx = qfx_cv(jsf,i,m1a)/ama(ia)
            erec = 13.595

            flsmi(it,5) = fsg*zflx*erec*cev
            flsmi(it,6) = flsmi(it,5) + dwrad(i,2)
            flsmi(it,7) = flsmi(it,6) + dwntl(i,2)
      enddo  !  loop (it)
     
!     For "Qdpl_o.txt"
      do it = 2, itpve-1
            j   = jcel(1,it)
            i   = icel(1,it)
            zsv = fareo(it)    ! S^sita 
            zsp = hvsb(j,i)  ! S^sita*h*cosw
            hcsw = zsp/zsv

            latestVal_o(1,it-1) = flsmo(it,1)/zsp*hcsw ! qt_ev
            latestVal_o(2,it-1) = flsmo(it,2)/zsp*hcsw ! qt_ed
            latestVal_o(3,it-1) = flsmo(it,3)/zsp*hcsw ! qt_iv
            latestVal_o(4,it-1) = flsmo(it,4)/zsp*hcsw ! qt_id
            latestVal_o(5,it-1) = flsmo(it,5)/zsp*hcsw ! qt_rec
            latestVal_o(6,it-1) = vne(j,i) ! Ned
            latestVal_o(7,it-1) = vte(j,i) ! Ted
            latestVal_o(8,it-1) = vti(j,i) ! Tid
            latestVal_o(9,it-1) = vva(j,i,ia) ! Vid
            latestVal_o(10,it-1) = vni(j,i) !vni
      enddo

!     For "Qdpl_i.txt"
      do it = 2, itpve-1
            jmax = jtmax(it)
            j   = jcel(jmax,it)
            i   = icel(jmax,it)
            jsf = jcel(jmax-1,it)
            zsv = farei(it)    ! S^sita 
            zsp = hvsb(jsf,i)  ! S^sita*h*cosw
            hcsw = zsp/zsv

            latestVal_i(1,it-1) = flsmi(it,1)/zsp*hcsw ! qt_ev
            latestVal_i(2,it-1) = flsmi(it,2)/zsp*hcsw ! qt_ed
            latestVal_i(3,it-1) = flsmi(it,3)/zsp*hcsw ! qt_iv
            latestVal_i(4,it-1) = flsmi(it,4)/zsp*hcsw ! qt_id
            latestVal_i(5,it-1) = flsmi(it,5)/zsp*hcsw ! qt_rec
            latestVal_i(6,it-1) = vne(j,i) ! Ned
            latestVal_i(7,it-1) = vte(j,i) ! Ted
            latestVal_i(8,it-1) = vti(j,i) ! Tid
            latestVal_i(9,it-1) = vva(j,i,ia) ! Vid
            latestVal_i(10,it-1) = vni(j,i) !vni
      enddo

      end subroutine calcLatestVal

!***********************************************************************
      subroutine saveLatest
!***********************************************************************
!     FIFO:Store the result of latestVal_o,i in timeSeries_o,i
      implicit none
      integer i_time

      if(seriesNow>0) then
         do i_time = 1, timeNum-1
            timeSeries_o(i_time,:,:) = timeSeries_o(i_time+1,:,:)
            timeSeries_i(i_time,:,:) = timeSeries_i(i_time+1,:,:)
         enddo
      endif
      timeSeries_o(timeNum,:,:) = latestVal_o
      timeSeries_i(timeNum,:,:) = latestVal_i
      if(seriesNow < timeNum) then
            seriesNow = seriesNow + 1
      endif

      end subroutine saveLatest

      end module mod_soldorTimeSeries