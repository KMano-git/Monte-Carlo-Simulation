      module ntpfctl

      implicit none

      contains
!**********************************************************************
      subroutine pfctl
!**********************************************************************
      use cntcom,   only : mcel
      use cntpfctl, only : fac_pfctl, fpfctl, lpfctl, lpfctl_rst
     >    , nsep_pfctl, rst_pfctl
      use cntpls,   only : dene
      use cplmet,   only : icel, itsle, jcel
      use cpmpls,   only : jmd1
      use csonic,   only : dtim, itim
      use cunit,    only : n6
      implicit none

      integer i,j,ic, stat
      real(8),save :: zfac = 1.0d0
      real(8) :: frac, tmp
      integer :: nf = 120
      logical, save :: lini=.true.

      fpfctl = 1.0d0

!      if(mygrp /= cdgr(2)) return

      if(lini)then
        tmp = 1.0d0
        if(lpfctl_rst == 1)then
          open(nf,file="rst_pfctl.txt",status='old',iostat=stat)
          read(nf,*,iostat=stat) tmp
          close(nf)
          if(stat==0) then
            rst_pfctl = tmp
            write(n6,'(/2x,"*** ntpfctl *** restart factor from file",
     >                 ": factor = ",1pe13.5)') tmp
          else
            write(n6,'(/2x,"*** ntpfctl cannot read rst_pfctl.txt.
     >                    rst_pfctl in inppls is used.")')
          endif
        endif

        if(rst_pfctl > 0.0d0) then
          zfac = rst_pfctl
          rst_pfctl = 0.0d0
        endif
        lini=.false.
      endif

      if(itim<1) then
        fpfctl(lpfctl) = zfac
        return
      endif

      j  = jcel(jmd1,itsle)
      i  = icel(jmd1,itsle)
      ic = mcel(j,i)

      frac = 1.0d0 - dtim*fac_pfctl*(dene(ic)-nsep_pfctl)/nsep_pfctl
      zfac = zfac*frac
      zfac = max(zfac, 1.0d-3)
      zfac = min(zfac, 1.0d+2)
      fpfctl(lpfctl) = zfac

!     write(n6,'(2x,"*pfctl itim,fac,fac_l,nse,nsep_pfctl: ",
!    >     i7,i4,4(1pe13.5),",")')
!    >     itim, lpfctl,zfac,frac,dene(ic),nsep_pfctl

      return
      end
!
!
!**********************************************************************
      subroutine out_pfctl
!**********************************************************************
      use cntcom,   only : mcel, pfmag
      use cntpfctl, only : fpfctl, lpfctl, nsep_pfctl
      use cntpls,   only : dene
      use cplmet,   only : icel, itsle, jcel
      use cpmpls,   only : jmd1
      use csonic,   only : itim
      use cunit,    only : n6, wh_dir, lmype, lmspe
      implicit none

! modified 2/2 lines organize local variables and include files by kamata 2021/11/18
!ik   integer nf
!ik   integer i,j,ic, stat
      integer :: nf = 120
      integer i,j,ic
      real(8) puff

! deleted 1 line organize local variables and include files by kamata 2021/11/18
!ik   nf = 120

! added 1 line for bug by kamata 2021/11/18
      if(lpfctl==0) return

      j  = jcel(jmd1,itsle)
      i  = icel(jmd1,itsle)
      ic = mcel(j,i)

      puff = pfmag(lpfctl)*fpfctl(lpfctl)

! deleted 1 line for bug by kamata 2021/11/18
!ik   if(lpfctl==0) return

      write(n6,'(/2x,"*** ntpfctl ***", "  nsep_pfctl =",1pe11.4,
     >      ", fpfctl(",i1,") =",1pe12.5,", gas puff rate = ",1pe12.5)')
     >     dene(ic), lpfctl,fpfctl(lpfctl),puff

      if(lmype==lmspe)then
        open(nf,file="../rst_pfctl.txt")
        write(nf,'(1pe16.8)')fpfctl(lpfctl)
        write(nf,*)
        write(nf,'("wroking dir: ",a,"  itim = ",i7)')trim(wh_dir),itim
        write(nf,'(3x,"lpfctl = ",i2,", nsep_pfctl =",1pe11.3,
     >      ", rst_pfctl =", 1pe16.8,",")')
     >     lpfctl, nsep_pfctl, fpfctl(lpfctl)
        write(nf,'(3x,"gas puff rate = ",1pe13.5)') puff

        close(nf)
      endif

      return
      end subroutine

      end module ntpfctl
