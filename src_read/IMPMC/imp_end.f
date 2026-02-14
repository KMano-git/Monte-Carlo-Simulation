!***********************************************************************
      subroutine imp_end
!***********************************************************************
      use cimcom, only : ismax, ndis, sflux, sptyc, swtot
      use cimctl, only : cdirz, ftav, icalz, lstdy
      use cimden, only : bkspt, csput, fsput, nsput, tdnz, tfrz, tionz
     >    , tradiz, tradliz, tradrz, trecz, tthz, tvlz, twci, twrd
     >    , wtsput, tengz
      use cimpuf, only : bk_mag
      use cntcom, only : nogt
      use csize,  only : ndmc
      use csonic, only : is_rsta, itim, limp, time
      use cunit,  only : n6
      implicit none
!
!::local variables
      integer :: it, nv, k
      integer :: lwrt
      real(8) :: wrd(ndmc), wci(ndmc), dnz(0:ndis,ndmc)
     > ,engz(0:ndis,ndmc)
      real(8) :: frz(0:ndis,ndmc), thz(0:ndis,ndmc)
      real(8) :: vlvlz(0:ndis,ndmc)
      real(8) :: ldenZ(0:ndis,ndmc), tdenZ(0:ndis,ndmc)
      real(8) :: lionZ(0:ndis,ndmc)
      real(8) :: lrecZ(0:ndis,ndmc)
      real(8) :: radi(0:ndis,ndmc), radli(0:ndis,ndmc)
      real(8) :: radr(0:ndis,ndmc)
      real(8) :: zsum
!
!KH 130404 for time smoothing
      real(8) ::   wrd0(ndmc), wci0(ndmc), dnz0(0:ndis,ndmc), ftav0
      real(8) ::   frz0(0:ndis,ndmc), thz0(0:ndis,ndmc)
      real(8) ::   vlvlz0(0:ndis,ndmc)
      real(8) ::   ionZ0(0:ndis,ndmc), recZ0(0:ndis,ndmc)
      real(8) :: radi0(0:ndis,ndmc), radli0(0:ndis,ndmc)
      real(8) :: radr0(0:ndis,ndmc)

      if( limp.ne.3 ) return
!
      write(n6,'(/30("*"),"  imp_end(start)   itim =",i6,
     >  "  time =",1pe12.4)') itim, time
!
!KH 130404 for time smoothing
      ftav0 = ftav
      if((itim.eq.0).and.(.not.is_rsta)) ftav0 = 0.0d0

!::previous results
!)  clear twrd, twci and tdnz in sub. imp_ini
      wrd0(1:ndmc) = twrd(1:ndmc)
      wci0(1:ndmc) = twci(1:ndmc)
      dnz0(0:ndis,1:ndmc) = tdnz(0:ndis,1:ndmc)
      frz0(0:ndis,1:ndmc) = tfrz(0:ndis,1:ndmc)
      thz0(0:ndis,1:ndmc) = tthz(0:ndis,1:ndmc)
      vlvlz0(0:ndis,1:ndmc) = tvlz(0:ndis,1:ndmc)
      ionZ0(0:ndis,1:ndmc) = tionZ(0:ndis,1:ndmc)
      recZ0(0:ndis,1:ndmc) = trecZ(0:ndis,1:ndmc)
      radi0(0:ndis,1:ndmc) = tradiZ(0:ndis,1:ndmc)
      radli0(0:ndis,1:ndmc) = tradliZ(0:ndis,1:ndmc)
      radr0(0:ndis,1:ndmc) = tradrZ(0:ndis,1:ndmc)

!::clear
      twrd(1:ndmc) = 0.0d0
      twci(1:ndmc) = 0.0d0
      tdnz(0:ndis,1:ndmc) = 0.0d0
      tengz(0:ndis,1:ndmc)= 0.0d0
      tfrz(0:ndis,1:ndmc) = 0.0d0
      tthz(0:ndis,1:ndmc) = 0.0d0
      tvlz(0:ndis,1:ndmc) = 0.0d0
      tdenZ(0:ndis,1:ndmc)= 0.0d0
      tionZ(0:ndis,1:ndmc)= 0.0d0
      tradiZ(0:ndis,1:ndmc)= 0.0d0
      tradliZ(0:ndis,1:ndmc)= 0.0d0
      tradrZ(0:ndis,1:ndmc)= 0.0d0
!
!::sum up for generations
      do k = 1, nsput
      sptyc = csput(k)
      sflux = fsput(k)
      swtot = wtsput(k)
!
      write(n6,'(/2x,i3,2x,a,2x,"sflux =",1pe12.4,"  swtot =",
     >   1pe12.4)') k, sptyc, sflux, swtot
!
      call imsave(sptyc,"r",k)
      call imsave(sptyc,"z",k)
      call calwrd("ionrec",nv,wrd,wci,dnz,engz)
      call calflw("ionrec",nv,frz,thz,vlvlz,ldenZ,lionZ,lrecZ,dnz,
     > radi,radli,radr)

      call dbg_vzpz("imp_end")

      twrd(1:nv) = twrd(1:nv) + wrd(1:nv)
      twci(1:nv) = twci(1:nv) + wci(1:nv)
      tdnz(0:ismax,1:nv) = tdnz(0:ismax,1:nv) + dnz(0:ismax,1:nv)
      tengz(0:ismax,1:nv)= tengz(0:ismax,1:nv)+ engz(0:ismax,1:nv)
      tdenZ(0:ismax,1:nv)= tdenZ(0:ismax,1:nv)+ ldenZ(0:ismax,1:nv)
      tfrz(0:ismax,1:nv) = tfrz(0:ismax,1:nv) + frz(0:ismax,1:nv)
      tthz(0:ismax,1:nv) = tthz(0:ismax,1:nv) + thz(0:ismax,1:nv)
      tvlz(0:ismax,1:nv) = tvlz(0:ismax,1:nv) + vlvlz(0:ismax,1:nv)
      tionZ(0:ismax,1:nv) = tionZ(0:ismax,1:nv) + lionZ(0:ismax,1:nv)
      trecZ(0:ismax,1:nv) = trecZ(0:ismax,1:nv) + lrecZ(0:ismax,1:nv)
      tradiZ(0:ismax,1:nv)= tradiZ(0:ismax,1:nv)+ radi(0:ismax,1:nv)
      tradliZ(0:ismax,1:nv)= tradliZ(0:ismax,1:nv)+ radli(0:ismax,1:nv)
      tradrZ(0:ismax,1:nv)= tradrZ(0:ismax,1:nv)+ radr(0:ismax,1:nv)
      enddo   !  loop(sputtering)

      ! Averaged over denZ Yamoto.
      ! Source averaging is also included.
      do it=1,nv
         do k=0,ismax
            if(tdenZ(k,it).gt.0.d0)then
               tfrz(k,it) = tfrz(k,it)/tdenZ(k,it)
               tthz(k,it) = tthz(k,it)/tdenZ(k,it)
               tvlz(k,it) = tvlz(k,it)/tdenZ(k,it)
            else
               tfrz(k,it) = 0.0d0
               tthz(k,it) = 0.0d0
               tvlz(k,it) = 0.0d0
            end if
         end do
      end do

      zsum = 0.d0
      do it=1,nv
         do k=0,ismax
            zsum = zsum + trecZ(k,it)
         enddo
      enddo
            write(n6,'(2x,"### dbg_trecZ",2x,a,2x,a,2x,
     >  "sflux, swtot, tot_trecZ =",1p3e12.3)')
     >  trim(sptyc), "averaged", sflux, swtot, zsum


!::smoothing new and old reults (default ftav0 = 0.0)
      twrd(1:nv) = (1.0d0-ftav0)*twrd(1:nv) + ftav0*wrd0(1:nv)
      twci(1:nv) = (1.0d0-ftav0)*twci(1:nv) + ftav0*wci0(1:nv)
      tdnz(0:ismax,1:nv) = (1.0d0-ftav0)*tdnz(0:ismax,1:nv)
     >       + ftav0*dnz0(0:ismax,1:nv)
      tfrz(0:ismax,1:nv) = (1.0d0-ftav0)*tfrz(0:ismax,1:nv)
     >       + ftav0*frz0(0:ismax,1:nv)
      tthz(0:ismax,1:nv) = (1.0d0-ftav0)*tthz(0:ismax,1:nv)
     >       + ftav0*thz0(0:ismax,1:nv)
      tvlz(0:ismax,1:nv) = (1.0d0-ftav0)*tvlz(0:ismax,1:nv)
     >       + ftav0*vlvlz0(0:ismax,1:nv)
      tionZ(0:ismax,1:nv) = (1.0d0-ftav0)*tionZ(0:ismax,1:nv)
     >       + ftav0*ionZ0(0:ismax,1:nv)
      trecZ(0:ismax,1:nv) = (1.0d0-ftav0)*trecZ(0:ismax,1:nv)
     >       + ftav0*recZ0(0:ismax,1:nv)
      tradiZ(0:ismax,1:nv) = (1.0d0-ftav0)*tradiZ(0:ismax,1:nv)
     >       + ftav0*radi0(0:ismax,1:nv)
      tradliZ(0:ismax,1:nv) = (1.0d0-ftav0)*tradliZ(0:ismax,1:nv)
     >       + ftav0*radli0(0:ismax,1:nv)
      tradrZ(0:ismax,1:nv) = (1.0d0-ftav0)*tradrZ(0:ismax,1:nv)
     >       + ftav0*radr0(0:ismax,1:nv)

      zsum = 0.d0
      do it=1,nv
         do k=0,ismax
            zsum = zsum + trecZ(k,it)
         enddo
      enddo
            write(n6,'(2x,"### dbg_trecZ",2x,a,2x,a,2x,
     >  "sflux, swtot, tot_trecZ =",1p3e12.3)')
     >  trim(sptyc), "ftav", sflux, swtot, zsum

!
!::new variables to soldor
      bkspt(1:nogt) = bk_mag(1:nogt)
!
!::select time
      if( lstdy.eq.1 ) lwrt = 1
!
!::mark
      write(n6,'(/30("="),"  IMPMC-END     itim =",i6,"  time =",
     >   1pe14.6,"  icalZ =",i4,"  cdirZ = ",a)') itim, time,
     >   icalZ, cdirZ
!
      return
      end
