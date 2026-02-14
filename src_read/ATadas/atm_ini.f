!*********************************************************************
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   subroutine atm_ini
      subroutine atm_ini( kis )
!*********************************************************************
!
!      original lnkimp_iniSnd.f => Mnkimp_iniSnd  (MPI_Send/Recv)
!                                  atm_ini (test)
!---------------------------------------------------------------------
      use catcom, only : ip_plt, ip_prb, ip_prc
      use cimcom, only : ndis
      use cplimp, only : ismaxl
      use cunit,  only : lmspe, lmype, n6
      implicit none

! added 3 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
! arguments
      integer, intent(in) :: kis
! kis : impurity number

!::local variables
      integer :: ntmax, i, k, iz, izmx, ip, it, isno
      real(8) :: anefix, atemin, atemax, atedlt, zte, zne
      real(8) :: tbte(201), tbne(201)
      real(8) :: acof(ndis)

      if( lmype /= lmspe ) goto 100

      call trmark("atm_ini","atm_pset")
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_pset
      call atm_pset( kis )
      write(n6,'(2x,"------ atm_pset in atm_ini")')
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   write(n6,'(2x,"ismax =",i3)') ismax
      write(n6,'(2x,"ismax =",i3)') ismaxl(kis)
!
      ntmax  = 11     ! <== 21
      anefix = 1.0d20
      atemin = 0.1d0
      atemax = 1.0d3
!
      atedlt = dlog10(atemax/atemin)/dfloat(ntmax-1)
      do i = 1, ntmax
      zte = atedlt*dfloat(i-1)
      zte = atemin*10.0d0**zte
      tbte(i) = zte
      tbne(i) = anefix
      enddo
!
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   izmx = min0( ismax, 10 )
      izmx = min0( ismaxl(kis), 10 )
      do k  = 1, 3
      if( k.eq.1 ) then
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ip = ip_plt
!ik   write(n6,'(2x,"Lz(Te) of ionization  ip_plt =",i2)') ip_plt
      ip = ip_plt(kis)
      write(n6,'(2x,"Lz(Te) of ionization  ip_plt =",i2)') ip_plt(kis)
      write(n6,'(2x,1x,"it",3x,"Ne",10x,"Te",10x,12(i2,10x))')
     >   (iz,iz=0,izmx)
      endif
      if( k.eq.2 ) then
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ip = ip_prb
!ik   write(n6,'(2x,"Lz(Te) of recomb.  ip_prb =",i2)') ip_prb
      ip = ip_prb(kis)
      write(n6,'(2x,"Lz(Te) of recomb.  ip_prb =",i2)') ip_prb(kis)
      write(n6,'(2x,1x,"it",3x,"Ne",10x,"Te",10x,12(i2,10x))')
     >   (iz,iz=1,izmx)
      endif
      if( k.eq.3 ) then
! modified 2/2 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   ip = ip_prc
!ik   write(n6,'(2x,"Lz(Te) of CXR  ip_prc =",i2)') ip_prc
      ip = ip_prc(kis)
      write(n6,'(2x,"Lz(Te) of CXR  ip_prc =",i2)') ip_prc(kis)
      cycle
      write(n6,'(2x,1x,"it",3x,"Ne",10x,"Te",10x,12(i2,10x))')
     >   (iz,iz=1,izmx)
      endif
!
      do it = 1, ntmax
      zne = tbne(it)
      zte = tbte(it)
      if( zte.le.0.0d0 ) cycle
! modified 1/1 lines treat 4 or more impurities with IMPMC by kamata 2022/04/15
!ik   call atm_eval2(ip,zne,zte,acof,isno,ndis)
      call atm_eval2(ip,zne,zte,acof,isno,ndis,kis)
      write(n6,'(2x,i3,1p15e12.3)')
     >   it, zne, zte, (acof(iz),iz=1,izmx)
      enddo
      enddo
      call trmark("atm_ini","return")
!
 100  continue
!
      return
      end
