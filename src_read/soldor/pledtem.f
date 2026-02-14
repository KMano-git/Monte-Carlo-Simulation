!**********************************************************************
      subroutine pledtem
!**********************************************************************
      use cplcom, only : vte, vti
      use cplmet, only : icel, itpve, itsle, itsls, jcel, jtmax, jtmin
      use cunit,  only : n6
      implicit none
!
!::local variables
! modified 1/1 lines organize local variables and include files by kamata 2021/05/31
!ik   integer :: it, jt, j, i, ia, jc, ic
      integer :: it, jt, j, i
      real(8) :: hte0, hti0
! deleted 1 line organize local variables and include files by kamata 2021/05/31
!ik   real(8) :: zne, zq3
!
      hte0 = 2.0d0
      hti0 = 2.0d0
      write(n6,'(/2x,"*** pledtem ***  temin/timin =",2f8.3)')
     >   hte0, hti0
!
      write(n6,'(2x,"Before edit")')
      write(n6,'(2x,"itsle =",i3)') itsle
      write(n6,'(2x,"Odiv Te =",0p15f8.3)')
     > (vte(jcel(jtmin(it),it),icel(jtmin(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Odiv Ti =",0p15f8.3)')
     > (vti(jcel(jtmin(it),it),icel(jtmin(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Idiv Te =",0p15f8.3)')
     > (vte(jcel(jtmax(it),it),icel(jtmax(it),it)),it=itsls,itpve)
      write(n6,'(2x,"Idiv Ti =",0p15f8.3)')
     > (vti(jcel(jtmax(it),it),icel(jtmax(it),it)),it=itsls,itpve)
!
      do it = itsls, itpve
      do jt = jtmin(it), jtmax(it)
        j = jcel(jt,it)
        i = icel(jt,it)
        vte(j,i) = dmax1( vte(j,i), hte0 )
        vti(j,i) = dmax1( vti(j,i), hti0 )
      enddo
      enddo
!
      return
      end
