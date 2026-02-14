!***********************************************************************
      subroutine set_ratez
!***********************************************************************
!
!     direct ionization rate for any atom    : cdi (/sec)
!     radiative & dielectronic recombination : cdr (/sec)
!
!    Note
!       change ionization & recombination data
!                       from kato-data to ADAS-data
!                                                 2001/4/20
!
!::KSFUJI   add ndis in arguments
!     call atm_eval2(ip_ion,zne,zte,acof,isno,ndis)
!
!-----------------------------------------------------------------------
      use catcom, only : ip_cxr, ip_ion, ip_plt, ip_prb, ip_rec, xdsn
      use cimcom, only : amz, cdi, cdr, cdlz_cx, cdlz_i, cdlz_r, irmax
     >    , ismax, ndis
      use cntpls, only : dene, teme, temi
      use cphcns, only : cev
      use cplcom, only : tn0, tt0
      use csize,  only : ndmc
      use csonic, only : itim
      use cunit,  only : n6
      implicit none
!
      integer  ic, iz, isno, iz1, lpout, mm, ie, ig
      integer  nemx, izmx
      real*8   zero, zne, zte
      real*8   emin, emax, edlt, erel, zvl2,  zsgv
      real*8   svmx(0:ndis), sgv(0:ndis), svcx(0:ndis)
      real*8   acof(ndis), acof2(ndis)
      integer, parameter :: icN0 = 67   ! see dbg_snflx.f
! function
      integer    imox, imoy
!
      write(n6,'(2x)')
      write(n6,'(2x,"*** set_ratez ***")')
      write(n6,'(2x,"ismax =",i3,"  ndis =",i3,"  irmax =",i6,
     >"  icN0 =",i5)') ismax, ndis, irmax, icN0
      if( ismax.gt.ndis ) call wexit("set_ratez","ismax.ft.ndis")
!
!-----------------------------------------------------------------------
!::clear
!-----------------------------------------------------------------------
      zero = 1.0d-30
      do ic = 0, ndmc
        do iz = 0, ndis
          cdi(iz,ic) = zero
          cdr(iz,ic) = zero
          cdLz_i(iz,ic) = 0.0d0
          cdLz_r(iz,ic) = 0.0d0
          cdLz_cx(iz,ic) = 0.0d0
        enddo
      enddo
!
      acof (1:ndis) = 0.0d0
      acof2(1:ndis) = 0.0d0
!
!-----------------------------------------------------------------------
!::ionization rate  <sv>  [m3/s]
!-----------------------------------------------------------------------
      if( ip_ion(1) > 0 ) then
        write(n6,'(2x,"-- ionization rate --  ",a)')
     >    trim(xdsn(ip_ion(1),1))
        do ic = 1, irmax
          zne = dene(ic)
          zte = teme(ic)
          if( zte.gt.0.0d0 ) then
            call atm_eval2( ip_ion(1), zne, zte, acof, isno, ndis, 1 )
            do iz = 1, ismax
              iz1 = iz-1
              cdi(iz1,ic) = dmax1(zne*acof(iz),zero)
            enddo
          endif
        enddo
      endif
!
!-----------------------------------------------------------------------
!::recombination rate  <sv>  [m3/s]
!-----------------------------------------------------------------------
      if( ip_rec(1) > 0 ) then
        write(n6,'(2x,"-- recombination rate --  ",a)')
     >    trim(xdsn(ip_rec(1),1))
        do ic = 1, irmax
          zne = dene(ic)
          zte = teme(ic)
          if( zte.gt.0.0d0 ) then
            call atm_eval2( ip_rec(1), zne, zte, acof, isno, ndis, 1 )
            do iz = 1, ismax
              cdr(iz,ic) = dmax1(zne*acof(iz),zero)
            enddo
          endif
        enddo
      endif
!
!-----------------------------------------------------------------------
!::CXR rate  <sv>(Er)  [m3/s]    Er = 1/2*mH*vrel^2
!-----------------------------------------------------------------------
      if( ip_cxr(1).gt.0 ) then
        write(n6,'(2x,"-- CXR rate --  ",a)') trim(xdsn(ip_cxr(1),1))
!
        call sgvcxr_set
!
        izmx = min0( ismax, 10 )
        nemx = 21
        emin = 0.1d0
        emax = 2.0d3
        edlt = dlog(emax/emin)/dfloat(nemx-1)
        write(n6,'(4x,"ie",3x,"erel",8x,10("sigv_",i2.2,5x,:))')
     >    (iz,iz=1,izmx)
        do ie = 1, nemx, 4
          erel = edlt*dfloat(ie-1)
          erel = emin*exp(erel)
          do iz = 1, ismax
            call sgvcxr(0,iz,erel,sgv(iz))
          enddo
          write(n6,'(2x,i4,1p12e12.3)')
     >      ie, erel, (sgv(iz),iz=1,izmx)
        enddo
      endif
!
!-----------------------------------------------------------------------
!::Lz(Te)
!-----------------------------------------------------------------------
      do ic = 1, irmax
        zne = dene(ic)
        zte = teme(ic)
        if( zte.gt.0.0d0 ) then
          call atm_eval2( ip_plt(1), zne, zte, acof, isno, ndis, 1 )
          call atm_eval2( ip_prb(1), zne, zte, acof2, isno, ndis, 1 )
          do iz = 1, ismax
            iz1 = iz-1
            cdLz_i(iz1,ic) = acof(iz)
            cdLz_r(iz, ic) = acof2(iz)
          enddo
        endif
      enddo
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      write(n6,'(2x)')
      lpout = 5
      izmx = min0( ismax, 10 )
      mm = irmax/lpout
      if( mm.le.0 ) mm = 1
!
      write(n6,'(2x,"ionization rate  itim =",i8)') itim
      do iz = 0, ismax
        svmx(iz) = 0.0d0
      enddo
      do ic = 1, irmax
        do iz = 0, ismax
          svmx(iz) = dmax1(svmx(iz),cdi(iz,ic))
        enddo
        if(teme(ic).le.0.0) cycle
        if(mod(ic,mm).eq.0 .or. ic.eq.icN0) then
          write(n6,'(2x,i6,2i5,1p15e11.2)') ic, imox(ic), imoy(ic),
     >     dene(ic), teme(ic),(cdi(iz,ic),iz=0,izmx)
        endif
      enddo
      write(n6,'(2x,32x,a,1p15e11.2)') "max : ",(svmx(iz),iz=0,izmx)
!
      write(n6,'(2x,"recombination rate  itim =",i8)') itim
      do iz = 0, ismax
        svmx(iz) = 0.0d0
      enddo
      do ic = 1, irmax
        do iz = 0, ismax
          svmx(iz) = dmax1(svmx(iz),cdr(iz,ic))
        enddo
        if(teme(ic).le.0.0) cycle
        if(mod(ic,mm).eq.0 .or. ic.eq.icN0) then
          write(n6,'(2x,i6,2i5,1p15e11.2)') ic, imox(ic), imoy(ic),
     >     dene(ic), teme(ic),(cdr(iz,ic),iz=0,izmx)
        endif
      enddo
      write(n6,'(2x,32x,a,1p15e11.2)') "max : ",(svmx(iz),iz=0,izmx)
!
      if( ip_cxr(1) > 0 ) then
        write(n6,
     >   '(2x,"charge exchange recombination rate  itim =",i8)') itim
        ig = 1
        do iz = 0, ismax
          svmx(iz) = 0.0d0
        enddo
        do ic = 1, irmax
          zvl2 = 2.0d0*temi(ic)*cev/amz
          do iz = 0, ismax
            zsgv = zero
            if( iz.ne.0 ) then
              call sgvcxr(ic,iz,zvl2,zsgv)
            endif
            svcx(iz) = dmax1( zsgv, zero )
            svmx(iz) = dmax1( svmx(iz), svcx(iz) )
          enddo
          if(teme(ic).le.0.0) cycle
          if(mod(ic,mm).eq.0 .or. ic.eq.icN0) then
            write(n6,'(2x,i6,2i5,1p15e11.2)') ic, imox(ic), imoy(ic),
     >      tN0(ic,ig), tT0(ic,ig), temi(ic), (svcx(iz),iz=0,izmx)
          endif
        enddo
        write(n6,'(2x,32x,a,1p15e11.2)') "max : ",(svmx(iz),iz=0,izmx)
      endif
!
      write(n6,'(2x,"ion-radiation rate Lz_i(Te)  itim =",i8)') itim
      do ic = 1, irmax
        if(teme(ic).le.0.0) cycle
        if(mod(ic,mm).eq.0 .or. ic.eq.icN0) then
          write(n6,'(2x,i6,2i5,1p15e11.2)') ic, imox(ic), imoy(ic),
     >     dene(ic), teme(ic), (cdLz_i(iz,ic),iz=0,izmx)
        endif
      enddo
!
      write(n6,'(2x,"rec-radiation rate Lz_r(Te)  itim =",i8)') itim
      do ic = 1, irmax
        if(teme(ic).le.0.0) cycle
        if(mod(ic,mm).eq.0 .or. ic.eq.icN0) then
          write(n6,'(2x,i6,2i5,1p15e11.2)') ic, imox(ic), imoy(ic),
     >     dene(ic), teme(ic), (cdLz_r(iz,ic),iz=0,izmx)
        endif
      enddo
!
      return
      end
