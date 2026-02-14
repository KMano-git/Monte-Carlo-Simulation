!**********************************************************************
      subroutine ntgdsk(nft,cact)
!**********************************************************************
      use cntcom,      only : bvx, bvy, bvz, chwl, cswl, cwal, e0wl
     >    , grax, grxp, gzax, gzxp, iaxs, icgt, icwl, igtw, ihwl, ikwl
     >    , impl, inwl, ipgt, iplx, iply, ipmp, ipwl, ispx, isxw, isyw
     >    , ivb1, ivb2, ivs1, ivs2, iwgt, iwl1, iwl2, ixyw, jdp1, jdp2
     >    , jsl1, jsl2, jvb1, jvb2, jxp1, jxp2, mcel, mclx, mcly, mcly2
     >    , mgrd, migx, migy, mknd, mrgn, mseg, mxjw, mxjw2, myit, ncmax
     >    , ncmax2, negt, next, ngmax, ngmax2, nhwl, nocl, nogd
     >    , nogdv, nogt, nowl, nowl2, npew, npmp, npsw, npwl, npwl2
     >    , nsgt, rcwl, rcywl, snwl, temwl, tywl, volm, xegt, xpnt, xsgt
     >    , yegt, ypnt, ysgt
      use csize,       only : ndad, ndgt, ndgtp, ndmc, ndmg, ndms, ndpm
     >    , ndwh, ndwl, ndwp, ndx, ndy, nvxd, nvyd
      use cunit,       only : lmspe, lmype, lnope, mywld, n6
      use mod_dtypdef, only : nddt, typ_dnam, typ_itag
      use mpi!,         only : mpi_bcast
      use mod_externalgrid
      implicit none
!
!::argument
      integer,   intent(in) :: nft
      character, intent(in) :: cact*(*)
!
!::local variables
      integer   nndmg, nndmc, nndms, nndx, nndy, nnvxd, nnvyd
      integer   nndwl, nndwh, nndwp
      integer   nndgt, nndgtp, nndad, nndpm
      integer   grid_size_tmp, cell_size_tmp, vac_ele_size_tmp
      integer   pri_ele_size_tmp,vac_grid_size_tmp,pri_grid_size_tmp
      integer   iw, ic, mmigx, mmigy, ks, kg, ngs, nge, ng
      integer   ip, icx, icy, i
      real*8    bdir, x_out, y_out
      character cver*80

! functions
      integer   lenx
!
!::MPI
      integer  ierr, itag, ityp
!
!::dummy
      integer noclv(0:nvxd,0:nvyd)
!
      noclv = 0
!
!::version
      cver = "/home/shimizu/Gmesh/soc/ntgdsk.f   04/04/30"
!
      write(n6,'(/2x,"*** ntgdsk ***  ",a)') cact(1:lenx(cact))
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))
!
      if( lmype.ne.lmspe ) goto 1000
!
!::master
      if( nft.le.0 ) then
      write(n6,'(2x,"nft =",i3)') nft
      call wexit("ntgdsk","nft.le.0")
      endif
!
!-----------------------------------------------------------------------
!::write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) then
      write(nft) cver
!-----
      write(nft) ndmg, ndmc, ndms, ndx, ndy, nvxd, nvyd
      write(nft) ! noclv,mxjw,mxjw2,mcly2 is not used in SONIC
     >    xpnt,ypnt,volm,bvx,bvy,bvz,grax,gzax,grxp,gzxp
     >   ,mclx,mcly,mcly2,mxjw,mxjw2,myit
     >   ,nogd,nocl,mcel,nogdv,noclv,mseg,mgrd,mknd
     >   ,mrgn,migx,migy,next,iplx,iply
     >   ,ngmax,ngmax2,ncmax,ncmax2
     >   ,jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
     >   ,iwl1,ispx,iwl2,impl,iaxs,ivs1,ivs2
     >   ,jvb1, jvb2, ivb1, ivb2
!----
      write(nft) ndwl, ndwh, ndwp
      write(nft) ! ixyw is not used in SONIC
     >   rcywl, temwl, rcwl, e0wl, cswl, snwl
     >  ,cwal, chwl, tywl
     >  ,npsw, npew, inwl, ikwl, ipwl
     >  ,icwl, isxw, isyw, ixyw, igtw
     >  ,ihwl, ipmp
     >  ,nowl, npwl, npmp, nowl2, npwl2, nhwl
!-----
      write(nft) ndgt, ndgtp, ndad, ndpm
      write(nft) ! xsgt,xegt,ysgt,yegt is not used in SONIC
     >   xsgt,xegt,ysgt,yegt
     >  ,nsgt,negt,iwgt,ipgt,icgt,nogt
      endif
!
!-----------------------------------------------------------------------
!::read
!-----------------------------------------------------------------------
      if( cact(1:1).eq."r" ) then
      read(nft) cver
!
      write(n6,'(2x,"version : ",a)') cver(1:lenx(cver))
!-----
      read(nft) nndmg, nndmc, nndms, nndx, nndy, nnvxd, nnvyd
!
      write(n6,'(2x,"dimension size (grid)")')
      write(n6,'(2x,"ndmg/dsk =",2i6,"  ndmc/dsk =",2i6,
     > "  ndms/dsk =",2i6)') ndmg, nndmg, ndmc, nndmc, ndms, nndms
      write(n6,'(2x,"ndx/dsk  =",2i6,"  ndy/dsk  =",2i6,
     > "  nvxd/dsk =",2i6,"  nvyd/dsk =",2i6)')
     > ndx,nndx, ndy,nndy, nvxd,nnvxd,nvyd,nnvyd
      if( ndmg.ne.nndmg .or. ndmc.ne.nndmc .or. ndms.ne.nndms .or.
     >    ndx.ne.nndx   .or. ndy.ne.nndy   .or.
     >    nvxd.ne.nnvxd .or. nvyd.ne.nnvyd ) goto 910
!
!::KSFUJI  confirm    mcel(j,i) = 0 for dummy cell
      read(nft)  ! noclv,mxjw,mxjw2,mcly2 is not used in SONIC
     >    xpnt,ypnt,volm,bvx,bvy,bvz,grax,gzax,grxp,gzxp
     >   ,mclx,mcly,mcly2,mxjw,mxjw2,myit
     >   ,nogd,nocl,mcel,nogdv,noclv,mseg,mgrd,mknd
     >   ,mrgn,migx,migy,next,iplx,iply
     >   ,ngmax,ngmax2,ncmax,ncmax2
     >   ,jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
     >   ,iwl1,ispx,iwl2,impl,iaxs,ivs1,ivs2
     >   ,jvb1, jvb2, ivb1, ivb2
!----
      read(nft) nndwl, nndwh, nndwp
!
      write(n6,'(2x,"dimension size (wall)")')
      write(n6,'(2x,"ndwl =",2i5,"  ndwh =",2i5,"  ndwp =",2i5)')
     >   ndwl,nndwl,ndwh,nndwh,ndwp,nndwp
      if( ndwl.ne.nndwl .or. ndwh.ne.nndwh .or. ndwp.ne.nndwp ) goto 910
!
      read(nft) !ixyw,npwl,nowl is not used in SONIC
     >   rcywl, temwl, rcwl, e0wl, cswl, snwl
     >  ,cwal, chwl, tywl
     >  ,npsw, npew, inwl, ikwl, ipwl
     >  ,icwl, isxw, isyw, ixyw, igtw
     >  ,ihwl, ipmp
     >  ,nowl, npwl, npmp, nowl2, npwl2, nhwl
!-----
      read(nft) nndgt, nndgtp, nndad, nndpm
!
      write(n6,'(2x,"dimension size (gate)")')
      write(n6,'(2x,"ndgt =",2i5,"  ndgtp =",2i5,"  ndad =",2i5,
     > "  ndpm =",2i5)') ndgt,nndgt,ndgtp,nndgtp,ndad,nndad,ndpm,nndpm
      if( ndgt.ne.nndgt .or. ndgtp.ne.nndgtp .or. ndad.ne.nndad .or.
     >    ndpm.ne.nndpm ) goto 910

      read(nft) ! xsgt,xegt,ysgt,yegt is not used in SONIC
     >   xsgt,xegt,ysgt,yegt
     >  ,nsgt,negt,iwgt,ipgt,icgt,nogt
!:: use mesh data of Gmesh with external tool
      if(use_exdata) then
          !:: external grid data
          read(nft) grid_size_tmp, cell_size_tmp, vac_ele_size_tmp
     >        ,pri_ele_size_tmp
          read(nft) vac_grid_size_tmp, pri_grid_size_tmp
          write(n6,'(2x,"dimension size (exgrid)")')
          write(n6,'(2x,"grid_size =",2i5,"  cell_size =",2i5
     >         ,"  vac_ele_size =",2i5,"  pri_ele_size =",2i5
     >         ,"  vac_grid_size =",2i5,"  pri_grid_size =",2i5)')
     >          grid_size,grid_size_tmp
     >         ,cell_size,cell_size_tmp
     >         ,vac_ele_size,vac_ele_size_tmp
     >         ,pri_ele_size,pri_ele_size_tmp
     >         ,vac_grid_size,vac_grid_size_tmp
     >         ,pri_grid_size,pri_grid_size_tmp
          if( grid_size.ne.grid_size_tmp
     >      .or. cell_size.ne.cell_size_tmp 
     >      .or. vac_ele_size.ne.vac_ele_size_tmp
     >      .or. pri_ele_size.ne.pri_ele_size_tmp 
     >      .or. vac_grid_size.ne.vac_grid_size_tmp
     >      .or. pri_grid_size.ne.pri_grid_size_tmp) then
              call wexit("ntgdsk","exdata dimension size error")
          endif
          read(nft)
     >       vgx_EX,vgy_EX,vac_grid_x,vac_grid_y,pri_grid_x,pri_grid_y
     >       ,mseg_vacume,mseg_pri,vac_element,pri_element
     >       ,mseg_subdiv,subdiv_cell
     >       ,ipgt2,ipwl2
     >       ,boundary_vac,boundary_pri
     >       ,j_save,i_save
     >       ,nvm,ksvm,kevm,ivem,kvem,cvnm,pvem
          endif ! use exdata or not
      endif ! read or write
!
!-----------------------------------------------------------------------
!::correction of b-direction
!-----------------------------------------------------------------------
      write(n6,'(2x,"[direction of b-field]")')
!-----
      call ntbdir(bdir)
!-----
      write(n6,'(2x)')
      if( bdir.lt.0.0 ) then
        write(n6,'(2x,"[change b-direction so as to b*e_sita > 0")')
        do i = 0, ndmc
          bvx(i) = -bvx(i)
          bvy(i) = -bvy(i)
        enddo
      endif
!
!::KSFUJI correction of gate
      if(.not.use_exdata) then
        call edtgate
        call lstgate
      endif
!
!----------------------------------------------------------------------
!::send ntgdsk
!----------------------------------------------------------------------
 1000 continue
      if( lnope.gt.1 ) then
        write(n6,'(5x,"passed MPI_Bcast")')
!
!:[:MPI_Bcast in ntgdsk]   cntcom (xpnt,emrk)  10/04/21
        call tbfind( nddt, typ_dnam, ityp, 'CNTGRD' )
        itag = typ_itag(ityp)
        call MPI_Bcast( grax, 1, itag, lmspe, mywld, ierr )
        call cntcom_lsr( 1 )
!
!::[MPI_Bcast in ntgdsk]  cntcom (rcywl,emrk)  10/04/21
        call tbfind( nddt, typ_dnam, ityp, 'CNTWAL' )
        itag = typ_itag(ityp)
        call MPI_Bcast( rcywl, 1, itag, lmspe, mywld, ierr )
        call cntcom_lsr( 2 )
!
!::[MPI_Bcast in ntgdsk]   cntcom (nsgt,emrk)  10/04/21
        call tbfind( nddt, typ_dnam, ityp, 'CNTGAT' )
        itag = typ_itag(ityp)
        call MPI_Bcast( nsgt, 1, itag, lmspe, mywld, ierr )
!
!::[MPI_Bcast in ntgdsk]   con_ntgrd,mod_exterlnalgrid,con_vcdat
        if(use_exdata) then
          call exgrid_lsr
        endif
      endif
!
!-----------------------------------------------------------------------
!::debug write
!-----------------------------------------------------------------------
      if( cact(1:1).eq."w" ) return
!
!--grid
      write(n6,'(5x,"jdp1 =",i3,"  jxp1 =",i3,"  jsl1 =",i3,
     >  "  jsl2 =",i3,"  jxp2 =",i3,"  jdp2 =",i3)')
     >   jdp1,jxp1,jsl1,jsl2,jxp2,jdp2
      write(n6,'(5x,"ivs1 =",i3,"  iwl1 =",i3,"  ispx =",i3,
     >  "  iwl2 =",i3,"  ivs2 =",i3,"  impl =",i3,"  iaxs =",i3)')
     >   ivs1,iwl1,ispx,iwl2,ivs2,impl,iaxs
      write(n6,'(5x,"jvb1 =",i3,"  jvb2 =",i3,"  ivb1 =",i3,
     >  "  ivb2 =",i3)') jvb1,jvb2,ivb1,ivb2
      write(n6,'(5x,"axis =",2f8.4,"  xpnt =",2f8.4)')
     >   grax,gzax,grxp,gzxp
      write(n6,'(5x,"nowl =",3i6,"  npwl =",3i6,"  npmp =",i5,
     > "  nhwl =",i3)') nowl,nowl2,ndwl,npwl,npwl2,ndwp,npmp,nhwl
!--wall
      if(.not. use_exdata) then
        write(n6,'(6x,"iw",4x,"knd",3x,"ip",4x,"ic",4x,"igx",3x,
     >     "igy",3x,"ipx",3x,"ipy",3x,"xp",6x,"yp",6x,"jdp",3x,"idp",
     >     3x,"rcy",3x,"ivid",2x,"igtw")')
        do iw = 1, npwl2, 50
          ic = icwl(iw)
          mmigx = 0
          mmigy = 0
          if( ic.ne.0 ) then
            mmigx = migx(ic)
            mmigy = migy(ic)
          endif
          write(n6,'(2x,8i6,2f8.3,2i6,f8.3,i6,2x,a,2x,a)')
     >     iw,ikwl(iw),ipwl(iw),icwl(iw),mmigx,mmigy,
     >     iplx(icwl(iw)),iply(icwl(iw)),xpnt(ipwl(iw)),ypnt(ipwl(iw)),
     >     isxw(iw),isyw(iw),rcwl(iw),igtw(iw),tywl(iw)
     >    ,chwl(ihwl(iw))
        enddo
      endif
!
!--gate
      write(n6,'(2x,"nogt =",i3)') nogt
      do ks = 1, 2
        if( ks.eq.1 ) write(n6,'(2x,a)') "gate in plasma wall"
        if( ks.eq.2 ) write(n6,'(2x,a)') "gate in vacume wall"
        do kg = 1, nogt
          ngs = nsgt(kg,ks)
          nge = negt(kg,ks)
          write(n6,'(2x,i2,2x,2i4,2x,2f9.4,2x,2f9.4)') kg
     >    ,ngs,nge,xsgt(kg,ks),ysgt(kg,ks),xegt(kg,ks),yegt(kg,ks)
          if( ngs.eq.0 ) cycle
          do ng = ngs, nge
            iw = iwgt(ng)
            ip = ipgt(ng)
            ic = icgt(ng)
            if( ic.gt.0 .and. .not.use_exdata) then
              icx = migx(ic)
              icy = migy(ic)
            else
              icx = 0
              icy = 0
            endif
            !
            if(.not. use_exdata) then
              x_out = xpnt(ip)
              y_out = ypnt(ip)
            elseif(ic.le.ncmax+vac_ele_size+pri_ele_size) then
              x_out = pri_grid_x(ip)
              y_out = pri_grid_y(ip)
            else
              x_out = vgx_EX(ip)
              y_out = vgy_EX(ip)    
            endif
            write(n6,'(2x,3i6,2f10.4,3i6)')
     >         ikwl(iw),iw,ip,x_out,y_out,ic,icx,icy
          enddo
          !:: add last grid point
          if(use_exdata) then
            iw = iwgt(nge)
            ip = ipgt2(nge)
            ic = icgt(nge)
            if(ic.le.ncmax+vac_ele_size+pri_ele_size) then
              x_out = pri_grid_x(ip)
              y_out = pri_grid_y(ip)
            else
              x_out = vgx_EX(ip)
              y_out = vgy_EX(ip)    
            endif
            write(n6,'(2x,3i6,2f10.4,3i6)')
     >         ikwl(iw),iw,ip,x_out,y_out,ic,0,0
          endif
        enddo
      enddo
!
      return
!
!::error
 910  continue
      call wexit("ntgdsk","dimension size error")
      end
