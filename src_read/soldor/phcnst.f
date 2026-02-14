!**********************************************************************
      subroutine phcnst
!**********************************************************************
      use cphcns, only : c_eps_zero, c_ten_pi_sqrt_pi, cev, cme, cmp
     >    , cmpev, cpi
      use cunit,  only : n6
      implicit none
!
      cpi = 4.0d0*datan(1.0d0)
      cev = 1.60210d-19
      cmp = 1.67252d-27
      cme = 9.10908d-31
      cmpev = cmp/cev
!
      c_ten_pi_sqrt_pi = 10.d0 * cpi * sqrt(cpi)
      c_eps_zero = 8.854187817d-12 
!
      write(n6,'(/2x,"*** phcnst ***")')
      write(n6,'(5x,"cpi =",1pe25.16)') cpi
      write(n6,'(5x,"physical  con. : cev,cmp,cme,cmpev =",1p4e12.3)')
     >   cev,cmp,cme,cmpev
      write(n6,
     >'(5x,"physical con. : c_ten_pi_sqrt_pi, c_eps_zero =",1p3e12.3)')
     >   c_ten_pi_sqrt_pi, c_eps_zero
      write(n6,*)
!
      return
      end
