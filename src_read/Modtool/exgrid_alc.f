! allocate variables in exgrid_alc
      subroutine exgrid_alc( kk )
      use mod_externalgrid
      use csize,  only : ndms, ndwp, ndgtp, ndmg
      implicit none
! arguments
      integer, intent(in) :: kk
! kk : flag = 1: allocate,    = 2: deallocate

! local variables
      integer    istat
      character  cmsg*80

      select case( kk )
      case( 1 )
          if( .not. allocated( vgx_EX ) ) then
              allocate(
     >          vgx_EX(grid_size),vgy_EX(grid_size)
     >        , vac_grid_x(vac_grid_size),vac_grid_y(vac_grid_size)
     >        , pri_grid_x(pri_grid_size),pri_grid_y(pri_grid_size)
     >        , mseg_vacume(vac_ele_size),mseg_pri(pri_ele_size)
     >        , vac_element(vac_ele_size,ndms)
     >        , pri_element(pri_ele_size,ndms)
     >        , mseg_subdiv(cell_size)
     >        , subdiv_cell(cell_size,ndms)
     >        , ipgt2(ndgtp),ipwl2(ndwp)
     >        , boundary_vac(vac_ele_size,ndms)
     >        , boundary_pri(pri_ele_size,ndms)
     >        , j_save(ndmg), i_save(ndmg)
     >        , stat = istat )
              if( istat /= 0 ) then
                 write(cmsg,'(a,i6)')
     >            'vgx_EX allocate error in exgrid_alc, istat = ', istat
                 call wexit( 'exgrid_alc', trim( cmsg ) )
              endif

! initial set
              vgx_EX(1:grid_size)= 0.0_8
              vgy_EX(1:grid_size)= 0.0_8
              vac_grid_x(1:vac_grid_size) = 0.0_8
              vac_grid_y(1:vac_grid_size) = 0.0_8
              pri_grid_x(1:pri_grid_size)= 0.0_8
              pri_grid_y(1:pri_grid_size)= 0.0_8
              mseg_vacume(1:vac_ele_size)=0
              mseg_pri(1:pri_ele_size)=0
              vac_element(1:vac_ele_size,1:ndms)=0
              pri_element(1:pri_ele_size,1:ndms)=0
              mseg_subdiv(1:cell_size)=0
              subdiv_cell(1:cell_size,1:ndms)=0
              ipgt2(1:ndgtp)=0
              ipwl2(1:ndwp) = 0
              boundary_vac(1:vac_ele_size,1:ndms)="--"
              boundary_pri(1:pri_ele_size,1:ndms)="--"
              j_save(1:ndmg)=0
              i_save(1:ndmg)=0
          endif
      case( 2 )
! deallocate
          if( allocated( vgx_EX ) ) then
              deallocate(vgx_EX,vgy_EX,vac_grid_x,vac_grid_y
     >        ,pri_grid_x,pri_grid_y
     >        , mseg_vacume,mseg_pri,vac_element
     >        , pri_element, mseg_subdiv
     >        , subdiv_cell, boundary_vac, boundary_pri
     >        , ipgt2, ipwl2, j_save, i_save
     >        , stat = istat )
              if( istat /= 0 ) then
                  write(cmsg,'(a,i6)')
     >             'vgx_EX deallocate error in exgrid, istat = ', istat
                  call wexit( 'exgrid_alc', trim( cmsg ) )
              endif
          endif
      end select

      return
      end
