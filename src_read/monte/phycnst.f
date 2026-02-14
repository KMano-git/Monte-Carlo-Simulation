!**********************************************************************
      subroutine phycnst
!**********************************************************************
      use com_phycnst, only : cev, cme, cmp, cpi, cspi, csq2
      use cunit,       only : n6
      implicit none
!
!::physical constant [MKS]
      cpi = 4.0d0*datan(1.0d0)
! modified 3/3 lines minor bug by kamata 2022/08/08
!ik   cev = 1.60210e-19
!ik   cmp = 1.67252e-27
!ik   cme = 9.10908e-31
      cev = 1.60210d-19
      cmp = 1.67252d-27
      cme = 9.10908d-31
      csq2 = sqrt(2.0d0)
      cspi = sqrt(cpi)
!
      write(n6,'(/2x,"*** phycnst ***")')
      write(n6,'(2x,"cpi =",1pe12.4,"  cev =",1pe12.4)') cpi, cev
!
      return
      end
