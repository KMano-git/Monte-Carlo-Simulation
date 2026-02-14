!**********************************************************************
      subroutine ntwinp(nft)
!**********************************************************************
!
!      set rcywl(ndwl), temwl(ndwl) from input data
!
!       lrf_vs, lrf_gt(4)  = 1   Atom ==> Mole reflection    06/02/01
!       lrf_sh  No use ==> use                               07/02/07
!
!----------------------------------------------------------------------
      use cntcom, only : chwl, irfw, nhwl, rcywl, temwl
      use csize,  only : ndwh
      use cunit,  only : n6
      implicit none
!
!::argument
      integer, intent(in) :: nft
!
!::local variables
      real*8 :: rcy_wl=0.0_8, rcy_vs=0.0_8, rcy_di=0.0_8
     >      , rcy_do=0.0_8, rcy_ws=0.0_8, tem_wl=0.0_8
     >      , tem_vs=0.0_8, tem_di=0.0_8, tem_do=0.0_8
     >      , rcy_gt(4)=0.0_8, rcy_pa(4)=0.0_8
     >      , rcy_sh(4)=0.0_8, tem_gt(4)=0.0_8
     >      , tem_al(4)=0.0_8, tem_pa(4)=0.0_8, tem_sh(4)=0.0_8
      real*8  :: rcy_al(4) = 0.98d0
      integer :: lrf_gt(4) = 0, lrf_vs = 0
      integer :: lrf_sh(4) = 0, lrf_wl = 0
!
      character cty*4
      integer iset(ndwh), i, kn, lenx, no
!
      namelist/uwlinp/
     >        rcy_wl, rcy_vs, rcy_di, rcy_do, rcy_ws
     >       ,tem_wl, tem_vs, tem_di, tem_do
     >       ,rcy_gt, rcy_al, rcy_pa, rcy_sh
     >       ,tem_gt, tem_al, tem_pa, tem_sh
     >       ,lrf_gt, lrf_vs
     >       ,lrf_sh, lrf_wl
!
      write(n6,'(/2x,"*** ntwinp ***   11/06/29")')
!
!::clear
      call setd( rcywl, ndwh, 0.0d0 )
      call setd( temwl, ndwh, 0.0d0 )
      call seti( irfw,  ndwh, 0 )
      call seti( iset,  ndwh, 0 )
!
!::input data
      read(nft,uwlinp)
!
      write(n6,'(2x,"wall type in mesh data   nhwl =",i3)') nhwl
      write(n6,'(10(2x,a))') (chwl(i),i=1,nhwl)
!
      do i = 1, nhwl
      cty = chwl(i)
!
!::plasma wall
      if( cty(1:1).eq."W" ) then
        read( cty(2:2),* ) kn
        if( lenx(cty).eq.2 ) then
          iset(i)  = 1
          rcywl(i) = rcy_wl
          temwl(i) = tem_wl
          irfw(i)  = lrf_wl
!---------
          if( cty(1:2).eq."W1" ) then
            rcywl(i) = rcy_ws
          endif
!---------
        elseif( cty(3:3).eq."d" .and. kn.eq.2 ) then
          iset(i)  = 1
          rcywl(i) = rcy_di
          temwl(i) = tem_di
          irfw(i)  = lrf_wl
        elseif( cty(3:3).eq."d" .and. kn.eq.4 ) then
          iset(i)  = 1
          rcywl(i) = rcy_do
          temwl(i) = tem_do
          irfw(i)  = lrf_wl
        elseif( cty(3:3).eq."g" ) then
          read( cty(4:4),* ) no
          if( no.le.0 .or. no.gt.4 ) goto 910
          iset(i)  = 1
          rcywl(i) = rcy_gt(no)
          temwl(i) = tem_gt(no)
          irfw(i)  = lrf_gt(no)
        else
          write(n6,'(2x,"invalid wall of plasma wall  ",a)') cty
        endif
!
!::vacume wall
      else
      if( cty(1:1).eq."V" ) then
        if( cty(3:3).ne." " ) cty = cty(3:4)
      endif
        if( cty(1:1).eq."V" ) then
           iset(i) = 1
           rcywl(i) = rcy_vs
           temwl(i) = tem_vs
           irfw(i)  = lrf_vs
        elseif( cty(1:1).eq."p" ) then
           read( cty(2:2),* ) no
           if( no.le.0 .or. no.gt.4 ) goto 910
           iset(i) = 1
           rcywl(i) = rcy_pa(no)
           temwl(i) = tem_pa(no)
           irfw(i)  = lrf_vs
        elseif( cty(1:1).eq."g" ) then
           read( cty(2:2),* ) no
           if( no.le.0 .or. no.gt.4 ) goto 910
           iset(i)  = 1
           rcywl(i) = rcy_gt(no)
           temwl(i) = tem_gt(no)
           irfw(i)  = lrf_gt(no)
        elseif( cty(1:1).eq."a" ) then
           read( cty(2:2),* ) no
           if( no.le.0 .or. no.gt.4 ) goto 910
           iset(i)  = 1
           rcywl(i) = rcy_al(no)
           temwl(i) = tem_al(no)
           irfw(i)  = lrf_vs
         elseif( cty(1:1).eq."s" ) then
           read( cty(2:2),* ) no
           iset(i) = 1
           rcywl(i) = rcy_sh(no)
           temwl(i) = tem_sh(no)
           irfw(i)  = lrf_sh(no)  ! lrf_vs ==> lrf_sh
         elseif( cty(3:3).eq."c" ) then
           call wexit("ntwinp","case of crio-pump under construction")
         else
           call wexit("ntwinp","Cannot understand")
         endif
      endif
      enddo
!
      return
!
 910  continue
      call wexit("ntwinp","no < 1 .or. no > 4")
      end
