!***********************************************************************
      subroutine gclval(cnam,nv,var,ncl,hcl)
!***********************************************************************
      use csize
      use cntcom
      use com_gmsiz
      use com_gpsta
      implicit none
!
!::argument
      character(*) :: cnam
      integer ::  nv, ncl
      real(4) ::  var(nv), hcl(30)
!
!::local variables
      character(80) :: clin, cinp, cpvl, ctmp
      integer ic, lenx, mji, jp1, jp2, i, ii, j, ip
      real*4  zmin, zmax, va, vb, dv, zhsv
      integer mjs(30), mje(30), kmx
!
!::min, max
 100  continue
      zmin =  1.0d30
      zmax = -1.0d30
      ii = 0
      do ic = 1, nv
      ip = mgrd(ic,1)
      if( (xpnt(ip)-rmi)*(xpnt(ip)-rmx).gt.0.0d0 ) cycle
      if( (ypnt(ip)-zmi)*(ypnt(ip)-zmx).gt.0.0d0 ) cycle
      if( var(ic).le.-0.1e30 ) cycle
      if( var(ic).ne.0.0d0 ) then
      zmin = min(zmin,var(ic))
      zmax = max(zmax,var(ic))
      ii = ii + 1
      endif
      enddo
!
!::contour values
      call linsep( cnam, " ,", kmx, mjs, mje, 30 )
      ctmp = cnam(mjs(1):mje(1))
      cpvl = " "
      do i = 1, nvar
      if( trim(ctmp).eq.trim(ctnm(i)) ) then
        cpvl = ctcv(i)
        exit
      endif
      enddo
!
!::message
      write(clin,'(a,"  ",i6,1p2e10.2)') 
     >   cnam(1:lenx(cnam)),ii,zmin, zmax
      write(6,'(2x,a)') trim(clin)
      call gdput("contour of "//clin(1:lenx(clin)))
!
      call gdput("  default :"//trim(cpvl))
      call gdget(">> enter val ==> ",cinp,mji)
      if( mji.eq.0 ) cinp = cpvl
      mji = lenx(cinp)
      write(6,'(2x," mji =",i2,"  cinp = [",a,"]")')
     >   mji, cinp(1:lenx(cinp))
      jp1 = index(cinp,"(")
      jp2 = index(cinp,")")
!
!::input "2.0,5.0,10.0,20.0,50.0"
      if( jp1.eq.0 .and. jp2.eq.0 ) then
      ii = 1
      do j = 1, mji
      if( cinp(j:j).eq."," ) ii = ii + 1
      enddo
      ncl = ii
      if( ii.gt.30 ) ncl = 30
      read(cinp(1:mji),*) (hcl(i),i=1,ncl)
!
      do j=1,ncl-1
        do i=j+1,ncl
          if( hcl(j) .gt. hcl(i) ) then
            zhsv   = hcl(i)
            hcl(i) = hcl(j)
            hcl(j) = zhsv
          endif
        enddo
      enddo
!
      goto 200
      endif
!
!::input "1.0,8.0,(1.0)"
      read(cinp(1:jp1-1),*) va, vb
      read(cinp(jp1+1:jp2-1),*) dv
      ii = 0
 50   continue
      ii = ii + 1
      hcl(ii) = va + dv*float(ii-1)
      if( hcl(ii).gt.vb-dv*0.9 ) goto 60
      if( ii.le.30 ) goto 50
 60   continue
      ncl = ii
!
!::debug write
 200  continue
      write(6,'(/2x,"*** gclval ***  [",a,"]")') cinp(1:mji)
      write(6,'(4x,"hcl =",1p10e10.2)') (hcl(i),i=1,ncl)
!
      return
      end
