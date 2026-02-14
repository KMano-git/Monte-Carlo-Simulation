!**********************************************************************
      subroutine ntfolw(zincx,zint1,zint2,ztim1,ztim2)
!**********************************************************************
!
!      folow sample particle for atom and molecular
!      is_atom_or_mole = 0:for atom, = 1:for molecular
!
!         istp : event of a particle
!         ievt : event in a cell                 !  2003/01/20
!
!::need many steps to hit wall in vacume region
!::need more steps to include EL-collision
!
!----------------------------------------------------------------------
      use cntcom, only : icpo, ievt
     >    , ilon, istp, ixpo, iypo, lsmax
     >    , ltrc, migx, migy
     >    , tmas, trct, vel
     >    , velx, vely, velz, weit, xpos, ypos, ncmax, mseg
     >    , is_atom_or_mole, ntmont_flg
      use cntwcn, only : wend
      use cunit,  only : n6
      use mod_externalgrid
      implicit none
!
!::arguments
      real*8,intent(inout) :: zincx(0:1), zint1(0:1), zint2(0:1)
     > , ztim1(0:1), ztim2(0:1)
!
!::local variables
      real*8  wtim, wpox, wpoy, zdlt
      integer icpn, ker, iend, i, ntmont_flg_bk, ic_tmp
      integer, parameter :: loop_max = 100000
      character(len=4) ntrcChar
!
!:: for atom
      character cend(6)*3
      data cend/"ion","abs","hwl","stp","err","end"/
!:: for molecular
      integer mtyp
      character cend_mole(5)*3
      data    cend_mole/"atm","stp","abs","2H+","err"/

!----------------------------------------------------------------------
!:: initialization
!----------------------------------------------------------------------
      if(is_atom_or_mole==0) then ! for atom
        ntmont_flg_bk = ntmont_flg
        ntmont_flg = 0
        if(ntmont_flg_bk==0) ztim1(0) = 0.0d0
        if( weit.le.0.0d0 ) then
          wpox = xpos + velx*ztim1(0)
          wpoy = ypos + vely*ztim1(0)
          call ntrc680( istp,wpox,wpoy,weit,icpo,"folw_"//cend(6) )
          return
        endif
        if(ntmont_flg_bk == 0) ievt = 0
        ntrcChar = "folw"
      else ! for molecular
      !::initial position of H2
        mtyp = 1
        ntrcChar = "mole"
      endif
      iend = 0
!
!----------------------------------------------------------------------
!::start after lounch or charge exchange event
!----------------------------------------------------------------------
      call initialize_rand(is_atom_or_mole,zincx,zint1,zint2,ztim1
     >  ,ztim2)
!
!----------------------------------------------------------------------
!::loop of track
!----------------------------------------------------------------------
      do i = 1, loop_max
        istp = istp + 1
        ievt = ievt + 1
!
!--debug
        wtim = ztim1(is_atom_or_mole)
        wpox = xpos + velx*wtim
        wpoy = ypos + vely*wtim
        ixpo = migx(icpo)
        iypo = migy(icpo)
!
        if(is_atom_or_mole==0) then ! for atom
        !--collision time
          call ntrc680( istp,wpox,wpoy,weit,icpo,ntrcChar )
          call ntcrs_AH0(icpo,vel)
        else
        !--cross section data ("D2")
          call ntcrs_MH0(icpo,vel)
          call ntrc680( istp,wpox,wpoy,weit,icpo,ntrcChar )
        endif
!
!--spend time in cell
        icpn = icpo
!
        if(.not. use_exdata) then ! ordinary mesh
          call sptim(xpos,ypos,velx,vely,icpn,ilon,wtim,wpox,wpoy,ker)
        else ! the mesh made by Gmesh and external tool
          if(icpn.le.ncmax) then ! core and SOL
            call sptim_ex(xpos,ypos,velx,vely
     >            ,icpn,ilon,wtim,wpox,wpoy,ker,mseg(icpn))
          elseif(icpn .le. ncmax+vac_ele_size) then ! vacume
            ic_tmp = icpn-ncmax
            call sptim_grid(xpos,ypos,velx,vely
     >            ,icpn,ilon,wtim,wpox,wpoy,ker
     >            ,vac_element,vac_ele_size
     >            ,mseg_vacume(ic_tmp)
     >            ,vac_grid_size,vac_grid_x,vac_grid_y
     >            ,ic_tmp)
          elseif(icpn .le. ncmax+vac_ele_size+pri_ele_size) then ! private
            ic_tmp = icpn-(ncmax+vac_ele_size)
            call sptim_grid(xpos,ypos,velx,vely
     >            ,icpn,ilon,wtim,wpox,wpoy,ker
     >            ,pri_element,pri_ele_size
     >            ,mseg_pri(ic_tmp)
     >            ,pri_grid_size,pri_grid_x,pri_grid_y
     >            ,ic_tmp)
          else !subdiverter
            ic_tmp = icpn-(ncmax+vac_ele_size+pri_ele_size)
            call sptim_grid(xpos,ypos,velx,vely
     >            ,icpn,ilon,wtim,wpox,wpoy,ker
     >            ,subdiv_cell,cell_size
     >            ,mseg_subdiv(ic_tmp)
     >            ,grid_size,vgx_EX,vgy_EX
     >            ,ic_tmp)
          endif
        endif
!
        if( ker.ne.0 ) then
          call terminate_ntfolw(is_atom_or_mole,iend,5,5)
          exit
        endif
!
        ztim1(is_atom_or_mole) = wtim
        zdlt  = ztim1(is_atom_or_mole) - ztim2(is_atom_or_mole)
        zint1(is_atom_or_mole) = zint1(is_atom_or_mole) + zdlt*trct
        if( zint1(is_atom_or_mole).ge.zincx(is_atom_or_mole) ) then
          !::collision occur
          call collision_ntfolw(is_atom_or_mole,iend
     >       ,zincx,zint1,zint2,ztim1,ztim2,mtyp)
          !::loop end or continue
          if(iend>0)then
            exit
          else
            cycle
          endif
        endif

!--scoreing density & temperature
        call scoreing_value(is_atom_or_mole,zdlt,mtyp)
!
!--track estimator  ! before H2 move to next cell
        call ntscor_tr(zdlt)
!
!--next cell
        ztim2(is_atom_or_mole) = ztim1(is_atom_or_mole)
        zint2(is_atom_or_mole) = zint1(is_atom_or_mole)
        icpo  = icpn
!
!--hit the wall
        if( ilon.lt.0 ) then
          call hitwall_ntfolw(is_atom_or_mole,iend,zincx,zint1,zint2
     >       ,ztim1,ztim2,wpox,wpoy,mtyp)
          if(ntmont_flg==2) then
            return ! turn atom to mole in ntrefl
          endif
          if(iend > 0) then
            exit
          else
            cycle
          endif
        endif
!
!--too many loop
        if((ievt.gt.lsmax).or.(istp*ltrc.gt.lsmax)) then
          call terminate_ntfolw(is_atom_or_mole,iend, 4, 2)
          exit
        endif
!
        ievt = 0
      enddo
!:: end of main loop

!::too many loop
      if(i .ge. loop_max) then
        if(is_atom_or_mole==0) then
          wend(4) = wend(4) + weit
          iend = 4
        else
          wend(4) = wend(4) + weit*2.0d0
          weit  = 0.0d0
          iend = 2
        endif
        write(n6,'(2x,"too many step  istp,ievt,loop_max ",3i6,
     >      "  ixpo,iypo =",2i6)') istp,ievt,loop_max,ixpo,iypo
      endif
!
!----------------------------------------------------------------------
!::return
!----------------------------------------------------------------------
      if(is_atom_or_mole==0) then
        wpox = xpos + velx*ztim1(0)
        wpoy = ypos + vely*ztim1(0)
        call ntrc680(istp,wpox,wpoy,weit,icpo,"folw_"//cend(iend))
      else
        call ntrc680(istp,xpos,ypos,weit,icpo,"mole_"//cend_mole(iend))
      endif
      return
      end subroutine ntfolw

!**********************************************************************
      subroutine terminate_ntfolw(is_atom,iend,in_atom,in_mole)
!**********************************************************************
      implicit none
! argument
      integer, intent(in) :: is_atom ! 0:atom, 1:mole
      integer, intent(inout) :: iend
      integer, intent(in) :: in_atom,in_mole
!
      if(is_atom == 0) then
        call terminate_atom(iend,in_atom)
      else 
        call terminate_mole(iend,in_mole)
      endif
      end subroutine terminate_ntfolw

!**********************************************************************
      subroutine collision_ntfolw(is_atom,iend,zincx,zint1,zint2
     >      ,ztim1,ztim2,mtyp)
!**********************************************************************
      implicit none
! argument
      integer, intent(in) :: is_atom ! 0:atom, 1:mole
      integer, intent(inout) :: iend,mtyp
      real*8, intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     >  ,ztim1(0:1),ztim2
!
      if(is_atom == 0) then
        call collision_atom(iend,zincx,zint1,zint2,ztim1,ztim2)
      else 
        call collision_mole(iend,zincx,zint1,zint2,ztim1,ztim2,mtyp)
      endif
      end subroutine collision_ntfolw

!**********************************************************************
      subroutine hitwall_ntfolw(is_atom,iend,zincx,zint1,zint2,ztim1
     >   ,ztim2,wpox,wpoy,mtyp)
!**********************************************************************
      implicit none
! argument
      integer, intent(in) :: is_atom ! 0:atom, 1:mole
      integer, intent(inout) :: iend,mtyp
      real*8, intent(inout) :: zincx(0:1),zint1(0:1),zint2(0:1)
     >  ,ztim1(0:1),ztim2(0:1)
      real*8, intent(in) :: wpox,wpoy
!
      if(is_atom == 0) then
        call hitwall_atom(iend,zincx,zint1,zint2,ztim1,ztim2
     >       ,wpox,wpoy)
      else 
        call hitwall_mole(iend,zincx,zint1,zint2,ztim1,ztim2
     >       ,wpox,wpoy,mtyp)
      endif
      end subroutine hitwall_ntfolw

!**********************************************************************
      subroutine scoreing_value(is_atom,zdlt,mtyp)
!**********************************************************************
      use cntcom, only : bvx, bvy, bvz
     >    , icpo, igas
     >    , vel2, velp
     >    , velx, vely, velz, weit
      use cntwcn, only : wden, weng, wvlp, wnfl0x, wnfl0y, wnfl0z
     >    , wgden, wgeng, wnflgx, wnflgy, wnflgz
      implicit none

! arguments
      integer, intent(in) :: is_atom ! 0:atom, 1:mole
      integer, intent(in) :: mtyp
      real*8, intent(in) :: zdlt
! local
      real*8 zdn, zde, zdv

      velp = velx*bvx(icpo)+vely*bvy(icpo)+velz*bvz(icpo)
      zdn = zdlt*weit
      zde = zdn*vel2
      if(is_atom == 0) then
        zdv = zdn*velp
      !--scoreing density & temperature
        wden(icpo,igas) = wden(icpo,igas) + zdn
        weng(icpo,igas) = weng(icpo,igas) + zde
        wvlp(icpo,igas) = wvlp(icpo,igas) + zdv
      ! count neutral flow, 2016.06.23 toku
        wnfl0x(icpo,igas)=wnfl0x(icpo,igas)+zdn*velx
        wnfl0y(icpo,igas)=wnfl0y(icpo,igas)+zdn*vely
        wnfl0z(icpo,igas)=wnfl0z(icpo,igas)+zdn*velz
      else
      !--scoreing density & temperature
        wgden(icpo,mtyp) = wgden(icpo,mtyp) + zdn
        wgeng(icpo,mtyp) = wgeng(icpo,mtyp) + zde
      ! count neutral flow, 2016.06.23 toku
        wnflgx(icpo,igas)=wnflgx(icpo,igas)+zdn*velx
        wnflgy(icpo,igas)=wnflgy(icpo,igas)+zdn*vely
        wnflgz(icpo,igas)=wnflgz(icpo,igas)+zdn*velz
      endif
      end subroutine scoreing_value

!**********************************************************************
      subroutine initialize_rand(i_atom,zincx,zint1,zint2,ztim1,ztim2)
!**********************************************************************
!----------------------------------------------------------------------
!::start after lounch or charge exchange event
!----------------------------------------------------------------------
      implicit none
! arguments
      integer, intent(in) :: i_atom
      real*8, intent(out) :: zincx(0:1), zint1(0:1), zint2(0:1)
     > , ztim1(0:1), ztim2(0:1)
! function
      real(8)    random
!::
      zincx(i_atom) = -dlog(random(0))
      zint1(i_atom) = 0.0d0
      zint2(i_atom) = 0.0d0
      ztim1(i_atom) = 0.0d0
      ztim2(i_atom) = 0.0d0
      end subroutine initialize_rand