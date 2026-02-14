!***********************************************************************
      subroutine gdpl2dN(cnam,nv,var,creg)
!***********************************************************************
!
!       cnam : "n0"
!       nv   : size 
!       var  : varaibles
!       creg : "pla", "vac", "pla;vac","  "
!
!-----------------------------------------------------------------------
      use csize
      use cntcom
      use csonic
      use com_gmsiz
      use cunit
      implicit none

!::argument
      character ::  cnam*(*)
      character ::  creg*(*)
      integer   ::  nv
      real(4)   ::  var(nv)

!::local variables
      character :: ched*80
      integer :: kreg, ist, ien
      integer :: ncl, kclr(31)
      integer :: i, ic, kmax, ip, ik, iclr, nd, nsw
      real(4) :: hcl(30), cxp, cyp
      real(4) :: zx(10), zy(10), wx(1000), wy(1000)


      kreg = 0
      if( index(creg,"pla") > 0 ) kreg = kreg + 1
      if( index(creg,"vac") > 0 ) kreg = kreg + 2

      if( kreg == 0 ) then
         ist = 1
         ien = nv
      elseif( kreg == 1 ) then
         ist = 1
         ien = ncmax
      elseif( kreg == 2 ) then
         ist = ncmax + 1
         ien = ncmax2
      elseif( kreg == 3 ) then
         ist = 1
         ien = ncmax2
      endif

      write(6,'(//2x,"***  gdpl2dN ***  cnam = [",a,"]","  creg = [",a,
     >  "]")') trim(cnam), trim(creg)
      write(6,'(2x,a,2x,"kreg =",i2,"  ist,ien =",2i6)') 
     >  trim(creg), kreg, ist, ien


!
!::size
      call gmsiz(rmi,rmx,zmi,zmx,xln,yln)
!
!::colour value
      call gclval(cnam,nv,var,ncl,hcl)
!
!::colour
      call gdfps_new( ncl, kclr, 31 )
!
      write(n6,'(2x,"kclr =",10i4)') (kclr(i),i=1,ncl+1)
      write(n6,'(2x,"hcl  =",1p10e11.2)') (hcl(i),i=1,ncl)
!
!::loop
      do ic = 1, nv
!
!::void
      if( var(ic).lt.-0.1d30 ) cycle
!
      kmax = mseg(ic)
      do i = 1, kmax
      ip = mgrd(ic,i)
      zx(i) = xpnt(ip)
      zy(i) = ypnt(ip)
      enddo
!
      do i = 1, ncl
      ik = i
      if( var(ic).lt.hcl(i) ) goto 110
      enddo
      ik = ncl + 1
 110  continue
      iclr = kclr(ik)
!
!::debug
!xx      if( zx(1).ge.2.6d0  .and. zx(1).le.2.9 .and.
!xx     >    zy(1).ge.-2.4d0 .and. zy(1).le.-2.0d0 ) then
!xx        write(6,'(2x,"CHK gdpl2dN ",i7,2i5,i4,1pe12.3,1x,i3,1pe12.3)')
!xx     >    ic, imox(ic), imoy(ic), mrgn(ic), var(ic), iclr, hcl(iclr)
!xx      endif
!
      if( iclr.le.0 ) cycle
      call gdfpl( 1, kmax, zx, zy, iclr )
      enddo
!
!::legend (physical) ==> (paper)
!xx      call gdput("indicate the position of legend")
!xx      call gdmaus2(-1,ckey,cxp,cyp)
!-----
      call gdpcnv(1,pslgx,pslgy,cxp,cyp)
!-----
      call gdfleg_new(cxp,cyp,hcl,kclr,ncl,1)
      call gdfpe(1)
      call newpen(1)
!
!::wall
      nd = 1000
      call wslwid(2)  ! pen width
      call getwal("swal",nsw,wx,wy,nd)
      call gdplt(1,nsw,wx,wy,1,"  ")
!xx   call getwal("vwal",nsw,wx,wy,nd)
!xx   call gdplt(1,nsw,wx,wy,1,"  ")
      call getwal("mspx",nsw,wx,wy,nd)
      call gdplt(1,nsw,wx,wy,1,"  ")
      call getwal("psol",nsw,wx,wy,nd)
      call gdplt(1,nsw,wx,wy,1,"  ")
      call getwal("pprv",nsw,wx,wy,nd)
      call gdplt(1,nsw,wx,wy,1,"  ")
      call wslwid(1)
!
!::header & page
      write(ched,'(a,f9.3," msec")') trim(cnam),time*1.0d3
      call gdtbpo(4.0,18.2)
      call gdtbv1(trim(ched),99.0,' ')
      call gdpagh(1)
!
      return
      end

      
