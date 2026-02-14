!***********************************************************************
      subroutine cdf_getky(ckey)
!***********************************************************************
      use cdfcom, only : ctab, mtab, ndtb, nfcd, ntab, nvar_max, tbkey
     >    , tbrnk, tbtyp
      implicit none
!
!::arguments
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   character ckey*(*)
      character, intent(in) :: ckey*(*)

!::local varaibales
! modified 3/5 lines organize local variables and include files by kamata 2021/06/16
!ik   integer  ii, kk, ml, ns, n, lenx, me, i, nk, nlst
!ik   integer  mjs(110), mje(110), nrnk, numd
!ik   character  clin*150, cmsg*120, ctyp
      integer  ii, kk, ml, ns, n, me, i, nk, nlst
      integer  mjs(110), mje(110), numd
      character  clin*150, cmsg*120
! function
      integer    lenx
!
!::common table
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   call cdf_attr(ckey,ctyp,nrnk)
      call cdf_attr(ckey,tbtyp,tbrnk)
      tbkey = ckey
! deleted 2 lines organize local variables and include files by kamata 2021/06/16
!ik   tbtyp = ctyp
!ik   tbrnk = nrnk
!
      ml = lenx(ckey)
      ii = 0
      kk = 0
!
 100  continue
      read(nfcd,'(a)',end=910) clin
      call linsep(clin,'=",;',nk,mjs,mje,110)
      if( nk.le.0 ) goto 100
      if( clin(mjs(1):mje(1)).eq.ckey(1:ml) ) kk = 1
      if( kk.gt.0 ) then
      ns = 1; if(kk.eq.1) ns = 2
!--
      do n = ns, nk
      ii = ii + 1
      if( ii.gt.ndtb ) goto 920
      ctab(ii) = clin(mjs(n):mje(n))
      mtab(ii) = lenx(ctab(ii))
      enddo
!--
      kk = kk + 1
      me = lenx(clin)
      if( clin(me:me).eq.";" ) goto 110
      endif
      goto 100
!
 110  continue
      ntab = ii
!
!::check
      numd = tbrnk*nvar_max
      if( tbrnk.eq.0 ) numd = 1
      if( ntab.ne.numd ) then
      write(6,'(2x,"*** cdf_getky ***  ",a18,"  typ =",a,"  rnk =",i2,
     >  "  ntab =",2i6)') tbkey(1:lenx(tbkey)),tbtyp,tbrnk,ntab,numd
      call wexit("cdf_getky","error data size")
      endif
!
      return
!
      write(6,'(2x,"*** cdf_getky ***  ",a18,"  typ =",a,"  rnk =",i2,
     >  "  ntab =",2i6)') tbkey(1:lenx(tbkey)),tbtyp,tbrnk,ntab,numd
!
      write(6,'(2x,"*** cdf_getky ***  ",a,2x,i4)')
     >  tbkey(1:lenx(tbkey)), ntab
      write(6,'(4x,"typ =",a,"  rnk =",i2)') tbtyp,tbrnk
      nlst = ntab
      write(6,'(4(4x,a20))') (ctab(i)(1:mtab(i)),i=1,nlst)
!
      return
!
!::error
 910  continue
! modified 1/1 lines organize local variables and include files by kamata 2021/06/16
!ik   call wexit("cdf_getky","no found "//ckey(1:lenx(ckey)))
      call wexit("cdf_getky","no found "//ckey(1:ml))
 920  continue
      write(cmsg,'("too many data  ii > ndtb  ",2i7)') ii,ndtb
      call wexit("cdf_getky",cmsg)
      end
