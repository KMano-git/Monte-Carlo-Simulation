!***********************************************************************
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   subroutine imlimt(ipnm,icnd)
      subroutine imlimt(ip,icnd)
!***********************************************************************
      use cimcom, only : amz, ipcm, ir, is, lpcm, r1lim, r2lim, rr, vv
     >    , vz, z1lim, z2lim, zz
      use cphcns, only : cev
      use cunit,  only : n6
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer  ipnm, icnd
      integer, intent(in)  :: ip
      integer, intent(out) :: icnd
!
!::local variables
! modified 1/2 lines organize local variables and include files by kamata 2021/06/28
!ik   integer i, ip, imox, imoy
! function
      integer    imox, imoy
!
!
      icnd = 0
!
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   ip = ipnm
      if( (rr(ip)-r1lim)*(rr(ip)-r2lim).gt.0.0d0 .or.
     >    (zz(ip)-z1lim)*(zz(ip)-z2lim).gt.0.0d0 ) then
!-----
! modified 2/2 lines with TOPICS by kamata 2021/12/22
!ik   write(n6,'(2x,"error imlimt  out of system  ",
!ik  >  "lpcm =",i6,"  ipcm =",i6,"  is =",i2,"  ic =",
      write(n6,'(2x,"error imlimt  out of system  ",
     >  "lpcm =",i10,"  ipcm =",i6,"  is =",i2,"  ic =",
     >  i7,2i5,"  rr,zz =",1p2e12.3,"  vv,vz,Eng =",1p3e12.3)')
     >  lpcm, ipcm, is(ip), ir(ip),imox(ir(ip)),imoy(ir(ip)),
     >   rr(ip), zz(ip), vv(ip), vz(ip), 0.5d0*amz*vv(ip)/cev
!-----
      icnd = 1
      endif
!
      return
      end
!
!***********************************************************************
      subroutine set_limt
!***********************************************************************
!     z >  zg(nz-5)   in SA-calculation    2009/12/19
!-----------------------------------------------------------------------
      use cimcom,    only : r1lim, r2lim, z1lim, z2lim
      use com_eqdat, only : nr, nz, rg, zg
      use cunit,     only : n6
      implicit none
!

      write(n6,'(/2x,"*** imlimt ***")')
      r1lim = rg(1)
      r2lim = rg(nr)
      z1lim = zg(1)
      z2lim = zg(nz)
      write(n6,'(2x,"r1lim,r2lim,z1lim,z2lim =",4f8.3,"  eq_region =",
     >  4f8.3)') r1lim,r2lim,z1lim,z2lim, rg(1),rg(nr),zg(1),zg(nz)
!
      return
      end
