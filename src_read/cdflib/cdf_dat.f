!**********************************************************************
      subroutine cdf_dat(cvar,dnam,devl,drnk,daty
     >         ,wvar,wnum,wmin,wmax,wmlt,wunt,wspc, ndrw,ndty)
!**********************************************************************
!
!      cvar    : (i) char  variable name
!
!      dnam    : (o) char  data set name
!      devl    : (o) char  evaluation type of y-data
!      drnk    : (o) int   rank
!      daty    : (o) r8    y-data

!      wvar(4) : (o) char  name of    y & x-data
!      wnum(4) : (o) int   number of  y & x-data
!      wmin(4) : (o) r8    min of     y & x-data
!      wmax(4) : (o) r8    max of     y & x-data
!      wmlt(4) : (o) r8    mult of    y & x-data
!      wunt(4) : (o) char  unit of    y & x-data
!      wspc(4) : (o) char  spacing of y & x-data
!
!      ndrw    : (i) int   rank of y & x-data   ndrw = 4
!      ndty    : (i) int   dimension size of y-data
!
!----------------------------------------------------------------------
      use cdfcom, only : nfcd
      implicit none
!
!::argument
! modified 6/8 lines organize local variables and include files by kamata 2021/06/16
!ik   integer    ndrw,ndty
!ik   character  cvar*(*), dnam*(*), devl*(*)
!ik   character  wvar(ndrw)*(*), wunt(ndrw)*(*), wspc(ndrw)*(*)
!ik   integer    drnk, wnum(ndrw)
!ik   real*8     wmin(ndrw), wmax(ndrw), wmlt(ndrw)
!ik   real*8     daty(ndty)
      integer,   intent(in)  :: ndrw, ndty
      character, intent(in)  :: cvar*(*)
      character, intent(out) :: dnam*(*), devl*(*)
      character, intent(out) :: wvar(ndrw)*(*), wunt(ndrw)*(*)
     >                        , wspc(ndrw)*(*)
      integer,   intent(out) :: drnk, wnum(ndrw)
      real(8),   intent(out) :: wmin(ndrw), wmax(ndrw), wmlt(ndrw)
      real(8),   intent(out) :: daty(ndty)

!::local variables
      character ckey*40
      integer   idim; parameter (idim=100)
! modified 2/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   ivar, lenx, num, i0, i, itmp(idim)
!ik   character ctmp(idim)*40; integer mtmp(idim)
      integer   ivar, i
!
      integer   inum, ybas, yinc, ynum, ynum2, nw
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   integer   i1, i2, i3, i4
      integer   i1, i4
      character cmsg*80
!
      integer ndvd; parameter (ndvd=4)
      character cvdt(ndvd)*40
      integer nvdt,mvdt(ndvd)
! added 2 lines organize local variables and include files by kamata 2021/06/16
! function
      integer    lenx
!
!x      write(6,'(/2x,"*** cdf_dat ***   cvar = ",a)') cvar(1:lenx(cvar))
!
!::check
      write(6,'(2x,a,2x,a," [",a,"]",2x,a,i3)') "[DBG] cdf_dat",
     >   "cvar =",trim(cvar), "ndrw =",ndrw
!
!::ivar
      call cdf_inqv(cvar,ivar)
      if( ivar.le.0 ) then
        call wexit("cdf_dat","no found cvar  "//cvar(1:lenx(cvar)))
      endif
!
      rewind nfcd
 100  continue
!
!::rank
      ckey = "xs_rank"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      read(cvdt(1),*) drnk
      nw = drnk + 1
!x      write(6,'(4x,"drnk =",i2)') drnk
!
!::default
      wnum(1:ndrw) = 1
!
!::tab_index
      ckey = "xs_tab_index"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      ynum2 = 1
      do i = 1, drnk
      read(cvdt(i),*) inum
      ynum2 = ynum2*inum
      wnum(i+1) = inum
      enddo
      wnum(1) = ynum2
!x    write(6,'(4x,"wnum =",10i6)') (wnum(i),i=1,nw)
      write(6,'(2x,a,2x,a,10i7)') "[DBG] cdf_dat", "wnum =",wnum
!
!::data_base & data_inc
      ckey = "xs_data_base"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      read(cvdt(1),*) ybas
      ckey = "xs_data_inc"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      read(cvdt(1),*) yinc
!x      write(6,'(4x,"ybas =",i6,"  yinc =",i6)') ybas,yinc
!
!::name
      ckey = "xs_name"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      dnam = cvdt(1)
!x      write(6,'(4x,"dnam = [",a,"]")') dnam(1:lenx(dnam))
!
!::spacing
      ckey = "xs_spacing"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      do i = 1, nw
      wspc(i) = cvdt(i)
      enddo
!x      write(6,'(4x,"wspc =",5(2x,a))')
!x     >   (wspc(i)(1:lenx(wspc(i))),i=1,nw)
!
!::var
      ckey = "xs_var"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      do i = 1, nw
      wvar(i) = cvdt(i)
      enddo
!x      write(6,'(4x,"wvar =",5(2x,a))')
!x    >   (wvar(i)(1:lenx(wvar(i))),i=1,nw)
!
!::unit
      ckey = "xs_units"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      do i = 1, nw
      wunt(i) = cvdt(i)
      enddo
!x      write(6,'(4x,"wunt =",5(2x,a))')
!x     >   (wunt(i)(1:lenx(wunt(i))),i=1,nw)
!
!::eval
      ckey = "xs_eval_name"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      devl = cvdt(1)
!x      write(6,'(4x,"devl =",2x,a)') devl
!
!::min
      ckey = "xs_min"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      wmin(1) = 0.0d0
      if( drnk.ge.1 ) then
      do i = 1, drnk
      read(cvdt(i),*) wmin(i+1)
      enddo
      endif
!x      write(6,'(4x,"wmin =",1p10e11.3)') (wmin(i),i=1,nw)
!
!::max
      ckey = "xs_max"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      wmax(1) = 0.0d0
      if( drnk.ge.1 ) then
      do i = 1, drnk
      read(cvdt(i),*) wmax(i+1)
      enddo
      endif
!x      write(6,'(4x,"wmax =",1p10e11.3)') (wmax(i),i=1,nw)
!
!::mult
      ckey = "xs_mult"
      call cdf_getky(ckey)
      call cdf_getvr(ivar,nvdt,cvdt,mvdt,ndvd)
      do i = 1, nw
      read(cvdt(i),*) wmlt(i)
      enddo
!x      write(6,'(4x,"wmlt =",1p10e11.3)') (wmlt(i),i=1,nw)
!
!::dta_tab
      ckey = "xs_data_tab"
      call cdf_getdt(ckey,ybas,yinc,ynum,daty,ndty)
! modified 2/2 lines organize local variables and include files by kamata 2021/06/16
!ik   i1 = 1; i2 = min0(i1+5,ynum)
!ik   i3 = max0( 1, ynum-5 ); i4 = ynum
      i1 = 1
      i4 = ynum
!
!x      write(6,'(4x,"ynum =",2i7)') ynum,ynum2
!x      write(6,'(4x,2i7,2x,1p6e12.3)') i1,i2,(daty(i),i=i1,i2)
!x      if( ynum.gt.6 )
!x     >write(6,'(4x,2i7,2x,1p6e12.3)') i3,i4,(daty(i),i=i3,i4)
!
!::debug write
      return
!
      write(6,'(/2x,"*** cdf_dat ***  [",a,"]  cvar =",a
     > ,2x,"devl =",a12,"  drnk =",i2,"  ndty =",i7)')
     >  dnam(1:lenx(dnam)), cvar(1:lenx(cvar)),devl,drnk,ndty
      do i = 1, nw
      write(6,'(4x,"wvar =",a18,"  wnum =",i6,"  wmin =",1pe10.3,
     >  "  wmax =",1pe10.3,"  wmlt =",1pe10.3,"  wunt =",a8,
     >  "  wspc =",a8)')
     >  wvar(i),wnum(i),wmin(i),wmax(i),wmlt(i),wunt(i),wspc(i)
      enddo
      if( ynum.gt.1 ) then
      write(6,'(4x,"ybas =",i6,"  yinc =",i6,"  ynum =",2i6,
     >     "  daty =",i6,1x,1pe10.3,"  daty =",i6,1x,1pe10.3)')
     >   ybas,yinc,ynum,ynum2,i1,daty(i1),i4,daty(i4)
      else
      write(6,'(4x,"ybas =",i6,"  yinc =",i6,"  ynum =",2i6,
     >     "  daty =",i6,1x,1pe10.3)')
     >   ybas,yinc,ynum,ynum2,i1,daty(i1)
      endif
!
!::error
      if( ynum.ne.ynum2 ) then
        write(cmsg,'("number of data ynum ",2i7)') ynum,ynum
        call wexit("cdf_dat",cmsg)
      endif
!
      return
      end
