!***********************************************************************
      module mod_externalgrid
!***********************************************************************
      use csize,only:ndwl
      implicit none
!----------------------------------------------------------------------
      logical :: use_exdata = .False.

      character(len=14), parameter :: 
     > vac_exdata_filename = "vac_mesh.txt"
     > ,pri_exdata_filename = "pri_mesh.txt"
      integer :: vac_grid_size = 0, pri_grid_size = 0
     > ,vac_ele_size = 0, pri_ele_size = 0
!:: grid
      real(8), allocatable :: vac_grid_x(:), vac_grid_y(:)
     > ,pri_grid_x(:), pri_grid_y(:)
!:: cell elements(cell index,grid index)
      integer, allocatable :: vac_element(:,:),pri_element(:,:)
!:: grid number in eath cell, may 3(=triangle) or 4(=quad)
      integer, allocatable :: mseg_vacume(:),mseg_pri(:)
!:: indices of adjacent cells
      integer, allocatable :: adjacent_vac(:,:),adjacent_pri(:,:)
!:: Boundary infomation of cells
      character(len=2), allocatable :: boundary_pri(:,:)
      character(len=2), allocatable :: boundary_vac(:,:)
      character(len=2), parameter :: label_none="--"
      character(len=2), parameter :: label_subdivwall="WD"
      character(len=2), parameter :: label_gate1="G1"
      character(len=2), parameter :: label_gate2="G2"
      character(len=2), parameter :: label_indiv="id"
      character(len=2), parameter :: label_outdiv="od"
      character(len=2), parameter :: label_priwall_other="Wo"
      character(len=2), parameter :: label_vac_vessel="VV"
      character(len=2), parameter :: label_SOLbound="SO"
      integer, parameter :: kindOfLabel = 8
      character(len=2), parameter :: label_list(kindOfLabel) = 
     >  (/ label_subdivwall,label_gate1,label_gate2,label_indiv
     > ,label_outdiv,label_priwall_other,label_vac_vessel
     > ,label_SOLbound /)
!
      logical :: is_compareEX = .False.
!
!*************************************************
!:: com_vcdat
!*************************************************
      integer, parameter :: ndem  =  20
      integer, parameter :: ndeg  = 300
      integer, parameter :: ndwpl = 400
      integer, parameter :: ndwv  = 200

!::subdivertor grid
      integer grid_size
      real*8, allocatable :: vgx_EX(:), vgy_EX(:)
!::subdivertor cell (indices of the grids comprising the cell)
      integer cell_size
      integer, allocatable :: subdiv_cell(:,:)
      integer, allocatable :: mseg_subdiv(:)
!::Indices of adjacent subdivertor cells
      integer, allocatable :: adjacent_subdiv(:,:)
!
!::division
      integer :: ivxs(10)=0, ivxe(10)=0, nvv=0
!
!::vacume object
      integer :: ksvm(ndem)=0, kevm(ndem)=0, nvm=0
      integer :: kvem(ndeg)=0, ivem(ndeg)=0, pvem(ndeg,ndem)=0 ! ivem should be allocatable
      integer :: nov=0, nog=0, nop=0, noc=0, noa=0, nos=0
      character :: cvnm(ndem)*2 = ''
      character :: cvem(ndeg)*2  = ''

!::plasma wall
      integer :: nwpl = 0
      real*8 :: wpx(ndwpl)=0.0_8, wpy(ndwpl)=0.0_8
      character :: wpc(ndwpl)*4=''

!::vacume wall
      integer :: nwv=0
      integer :: iwvs(20)=0, iwve(20)=0
      integer :: wvb(ndwv)=0

!----------------------------------------------------------------------
!::grid information for neut2d
!----------------------------------------------------------------------
      integer, allocatable :: j_save(:),i_save(:),ipgt2(:)
     > ,ipwl2(:)

!***********************************************************************
      contains
!***********************************************************************
      subroutine get_position_byic(ic,ig,is_x,xy_out)
!***********************************************************************
      use cntcom, only: xpnt, ypnt, ncmax
      implicit none
!:: arguments
      integer, intent(in) :: ic,ig
      logical, intent(in) :: is_x
      real(8), intent(out) :: xy_out

      if(ic.le.ncmax) then ! SOL region
        if(is_x) then
          xy_out = xpnt(ig)
        else
          xy_out = ypnt(ig)
        endif
      elseif(ic.le.ncmax+vac_ele_size) then ! vacume region
        if(is_x) then
          xy_out = vac_grid_x(ig)
        else
          xy_out = vac_grid_y(ig)
        endif
      elseif(ic.le.ncmax+vac_ele_size+pri_ele_size) then ! private region
        if(is_x) then
          xy_out = pri_grid_x(ig)
        else
          xy_out = pri_grid_y(ig)
        endif
      else ! sub-diverter
        if(is_x) then
          xy_out = vgx_EX(ig)
        else
          xy_out = vgy_EX(ig)
        endif
      endif
      end subroutine get_position_byic
!::

      end module mod_externalgrid