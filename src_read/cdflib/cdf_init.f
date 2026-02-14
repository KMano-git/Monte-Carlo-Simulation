!**********************************************************************
      subroutine cdf_init( nft, dsn )
!**********************************************************************
      use cdfcom, only : hvar, hvrk, hvty, mhvr, mxvr, ndat_max, ndxvr
     >    , nfcd, nhvr, nrnk_max, nvar_max, nxvr, xvar
      implicit none
!
!::argument
! modified 1/2 lines organize local variables and include files by kamata 2021/06/16
!ik   integer nft; character dsn*(*)
      integer,   intent(in) :: nft
      character, intent(in) :: dsn*(*)

!::local variables
      character clin*150,ckey*40
      integer  i,mj,i6,lenx,ktyp,ndim,irk
      character cdnm*40; integer inum
      character comt(5)*120;   integer ncmt
!
      character csep*7
      data  csep/" =;,()/"/
!
      integer ndsp;  parameter (ndsp=201)
      integer mjs(ndsp),mje(ndsp),nk
!
!::unit number
      nfcd = nft
      i6   = 6
!
      ktyp = 0
      ncmt = 0; nhvr = 0;
!
 100  continue
      read( nfcd,'(a)',end=910 ) clin
      call linsep(clin,csep,nk,mjs,mje,ndsp)
      if( nk.eq.0 ) goto 100
!
      ckey = clin(mjs(nk):mje(nk))
      mj   = lenx(ckey)
!
      if( ckey(1:mj).eq."dimensions:" ) then
        ktyp = 1
      elseif( ckey(1:mj).eq."variables:") then
        ktyp = 2
      elseif( ckey(1:mj).eq."attributes:" ) then
        ktyp = 3
      elseif( ckey(1:mj).eq."data:" ) then
        ktyp = 4
      else
      endif
!
      if( ktyp.eq.3 ) call linsep(clin,'=";/',nk,mjs,mje,ndsp)
!
 200  continue
!
!::debug write
!x      if( nk.eq.1 ) write(i6,'( )')
!x      write(i6,'(2x,i1,2x,10(i2,":","/"a,"/",2x))')
!x     >  ktyp,(i,clin(mjs(i):mje(i)),i=1,nk)
!
      if( ktyp.eq.4 ) goto 300
      if( nk.eq.1 )   goto 100
!
!::dimensions
      if( ktyp.eq.1 ) then
      cdnm = clin(mjs(1):mje(1))
      read( clin(mjs(2):mje(2)),* ) inum
      if( index(cdnm,"rank_ind0").gt.0 )   nrnk_max = inum
      if( index(cdnm,"dep_var_ind").gt.0 ) nvar_max = inum
      if( index(cdnm,"xs_data_ind").gt.0 ) ndat_max = inum
!
!::variables
      elseif( ktyp.eq.2 ) then
      irk = 0
      if( nk.lt.3 ) then
      elseif( nk.eq.3 ) then
       if( clin(mjs(3):mje(3)).eq."dep_var_ind" ) irk = 1
      else
       if( clin(mjs(3):mje(3)).eq."dep_var_ind" ) irk = 1
       if( clin(mjs(4):mje(4)).eq."rank_ind"  )   irk = nrnk_max-1
       if( clin(mjs(4):mje(4)).eq."rank_ind0" )   irk = nrnk_max
      endif
!
      nhvr = nhvr + 1
      hvty(nhvr) = clin(mjs(1):mje(1))
      hvar(nhvr) = clin(mjs(2):mje(2))
      hvrk(nhvr) = irk
      mhvr(nhvr) = lenx(hvar(nhvr))
!
!::attribute
      elseif( ktyp.eq.3 ) then
      ncmt = ncmt + 1
      comt(ncmt) = clin(mjs(2):mje(2))
      endif

      goto 100
!
 300  continue
!
!::variables
      ndim = ndxvr
      call cdf_var(nxvr,xvar,ndim)
      do i = 1, nxvr
      mxvr(i) = lenx(xvar(i))
      enddo
!
!::summary
      write(i6,'(/2x,"*** cdf_ini ***")')
      write(i6,'(2x,"file name =",a)') dsn(1:lenx(dsn))
      write(i6,'(2x,"comment")')
      do i = 1, ncmt
      write(i6,'(5x,a)') comt(i)(1:lenx(comt(i)))
      enddo
      write(i6,'(2x,"dimensions")')
      write(i6,'(5x,"nrnk_max =",i2,"  nvar_max =",i3,
     >  "  ndat_max =",i6)') nrnk_max, nvar_max, ndat_max
      write(i6,'(2x,"variables   nxvr =",i3)') nxvr
      write(i6,'(5(4x,a))') (xvar(i)(1:mxvr(i)),i=1,nxvr)
      write(i6,'(2x,"xs_variables  nhvr =",i3)') nhvr
      write(i6,'(4(4x,a15,2x,a,2x,i1))')
     >   (hvar(i)(1:mhvr(i)),hvty(i)(1:1),hvrk(i),i=1,nhvr)
!
!::file
      rewind nfcd
      return
!
 910  continue
      call wexit("cdf_init","wrong header")

      end
