!***********************************************************************
      subroutine yh_heat_flux_plot
!***********************************************************************
      use csize
      use cimcom
      use cntcom
      use cplcom
      use cplmet
      use cplqcn
      use cplvpn
      use cpmpls
      use cgdcom
      use csonic
      implicit none

  ! itmax
                        ! hdxm, hdxp, hpit
  ! vdda, vdet, vdxi, vdxe
  ! nqfx, nqfy
  ! vdvp
  ! adamp, aximp, axemp

 ! grdx, grdy
!
 ! mcel, ipls, iply
!
 ! heat_flux_para_at_ic_center
                       ! flmxi_at_impmc
!
!::local variables
      integer  it, jt, jts, jte, irg, ia, iy
      integer  j, i 
      integer  jchsl, icosl, icisl, i6
      real*8   fcenh, fcslw

      integer ic

c$$$      real*8 permittivity_eps_zero
c$$$      real*8 tau_ii, kappa_para, q_para_helander
c$$$      real*8 coulomb_log(ndx,ndy)
!
!::clear
!
c$$$      permittivity_eps_zero = 8.854187817d-12

!::heat flux q_para pm3d data
      open(3010,file='imc_heat_flux_para_pm3d.dat',status='replace')
      
      do it = 1, itmax
         jts = jtmin(it)
         jte = jtmax(it)
         if( it.ge.itmps .and. it.le.itmpe ) then
            jts = jts + 1
            jte = jte - 1
         endif
!     
         do jt = jts, jte
            j  = jcel(jt,it)
            i  = icel(jt,it)
            ic = mcel(j,i)
!            irg = kreg(j,i)
!     
!       q// helander evaluate (at cell center) for check 
!                                             2018/05/22 Y.Homma
c$$$            tau =
c$$$            kappa = 
c$$$            q_para_helander = 
!

            write(3010, '(10(E16.4e2))') 
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)), 
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)), 
     >           heat_flux_para_at_ic_center(ic) 
            write(3010, '(10(E16.4e2))') 
     >           grdx(kgdx(j,i,2), kgdy(j,i,2)), 
     >           grdy(kgdx(j,i,2), kgdy(j,i,2)),
     >           heat_flux_para_at_ic_center(ic) 
            write(3010, '(10(E16.4e2))') 
     >           grdx(kgdx(j,i,3), kgdy(j,i,3)), 
     >           grdy(kgdx(j,i,3), kgdy(j,i,3)),
     >           heat_flux_para_at_ic_center(ic) 
            write(3010, '(10(E16.4e2))') 
     >           grdx(kgdx(j,i,4), kgdy(j,i,4)), 
     >           grdy(kgdx(j,i,4), kgdy(j,i,4)),
     >           heat_flux_para_at_ic_center(ic) 
            write(3010, '(10(E16.4e2))')
     >           grdx(kgdx(j,i,1), kgdy(j,i,1)), 
     >           grdy(kgdx(j,i,1), kgdy(j,i,1)),
     >           heat_flux_para_at_ic_center(ic) 

            write(3010, '(2(E16.4e2))') 
         enddo

         write(3010, '(2(E16.4e2))') 
      enddo


      close(3010)


!
!::heat flux q_para poloidal plot date

      call tube_plot_immesh(42, heat_flux_para_at_ic_center,
     >     'imc_heat_flux_para_pol_edge_42.dat')

      call tube_plot_immesh(25, heat_flux_para_at_ic_center,
     >     'imc_heat_flux_para_pol_sepx_25.dat')

      call tube_plot_immesh(02, heat_flux_para_at_ic_center,
     >     'imc_heat_flux_para_pol_outermost_02.dat')

      call tube_plot_immesh(24, flmxi_at_impmc,
     >     'imc_flmxi_pol_sepx_24.dat')
!
!::heat flux q// elec send test  2018/10/10
      call tube_plot_immesh(24, heat_flux_para_at_ic_center,
     >     'impmc_qi_pol_24.dat')
      call tube_plot_immesh(24, heat_flux_para_e_ic_center,
     >     'impmc_qe_pol_24.dat')


!      stop

      return
      end subroutine 
!
!
!
!
!***********************************************************************
      subroutine tube_plot_immesh(it,array_immesh,name)
!***********************************************************************
      use csize
      use cimcom
      use cntcom
      use cplcom
      use cplmet
      use cplqcn
      use cplvpn
      use cpmpls
      use cgdcom
      use csonic
      implicit none

  ! itmax
                        ! hdxm, hdxp, hpit
 ! mcel, ipls, iply
!
      integer, intent(in) :: it
      real*8, intent(in) :: array_immesh(0:ndmc)
      character(*),intent(in) :: name
!::local variables
      integer  jt, jts, jte, irg, ia, iy
      integer  j, i 
      integer ic

      real*8 dlnx(ndx)
      

      ! Mesh connecting length
      call gdlen(it,dlnx)

      ! open file
      open(3020,file=name,status='replace')
      write(3020,*) "### it = ", it
      
      jts = jtmin(it)
      jte = jtmax(it)
      if( it.ge.itmps .and. it.le.itmpe ) then
         jts = jts + 1
         jte = jte - 1
      endif
!     
      do jt = jts, jte
         j  = jcel(jt,it)
         i  = icel(jt,it)
         ic = mcel(j,i)

         write(3020, '(2i,3(E16.4e2))')
     >        j, i, dlnx(jt), sum(dlnx(1:jt)), array_immesh(ic)

      enddo                     !     do jt = jts, jte

      ! close file
      close(3020)

      return
      end subroutine
