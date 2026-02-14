!**********************************************************************
      subroutine ntrc680(ist,xps,yps,wet,ic,cnd)
!**********************************************************************
      use BGProf, only : den0BG, dengBG, eng0BG, enggBG
      use cntcom, only : evel, ievt, iptl, irct, lnnel, ltrc, migx, migy
     >    , mrct, mrgn, nrct, srct, tmas, velx, vely
      use cntpls, only : dene, teme, temi
      use cphcns, only : cmp
      use csonic, only : itim
      use cunit,  only : cdrout, cmype, lmype, mjpe, n6
      implicit none
!
!::argument
      integer, intent(in)   :: ist, ic
      real(8), intent(in)   :: xps, yps, wet
      character, intent(in) :: cnd*(*)
!
!::local variables
      character cdsn*120,cmt*80
      real*8    zvlp
      integer   ix, iy, ir, i
!
      character crgn(0:9)*3
      data  crgn/"???","Odv","Sol","Idv","Opv",
     >           "Ipv","Edg","Man","Vac","???"/
      save  crgn
!
      integer nf
      data nf/ 87 /; save nf
      integer:: ncnt=0 
      save ncnt

      character ciptl*5
! function
      integer    lenx
!
      if( ltrc.eq.0 ) return
      if( lnnel==1 .and. itim<1)return
      if( lmype > 0) return
!
!::file allocation
      if( index(cnd,"ntrc_").gt.0 ) then
        if(iptl==1) ncnt=0
        ncnt=ncnt+1
        if( ncnt > 20) return
        write(ciptl,'(I0)')ncnt
        cdsn = cdrout(1:lenx(cdrout))//
     >         cnd(1:lenx(cnd))//"_"//cmype(mjpe:8)//"_"//trim(ciptl)
        write(n6,'("Open file for ntrc680  ",a)') cdsn(1:lenx(cdsn))
        open( unit=nf, file=cdsn )
        !ncnt = 0
        return
      endif
!
!      if( ncnt > ncntmx) return
!
!::start
      if( index(cnd,"ntsta").gt.0 ) then
!        ncnt = ncnt +1
!        if( ncnt > ncntmx) return
!
        write(nf,'(/"*** ",a," ***  iptl =",i6)') cnd(1:lenx(cnd)),iptl
        write(nf,'(a7,2x,a4,3a10,2x,a4,2x,a2,2a5,2x,a10,2x,a6,
     >   8a10,1x,a4,1x,a4,2x,a)')
     >    "ist","ievt","xps","yps","wet","crgn","ir","ix","iy"
     >    ,"cnd","M","evel","dene","teme","temi"
     >    ,"den0","eng0","denG","engG","mrct","irct","mfp"
        return
      endif

! default
      zvlp = dsqrt(velx*velx+vely*vely)
!
      cmt = cnd(1:lenx(cnd))
      if( ievt.gt.0 .and. mod(ievt,500).eq.0 )
     >  cmt = cmt(1:lenx(cmt))//"#"
!
      if( ic.le.0 ) return
!
      ix = migx(ic)
      iy = migy(ic)
      ir = mrgn(ic)
!
      write(nf,680) ist,ievt,xps,yps,wet,crgn(ir),ir,ix,iy
     >  ,cmt(1:lenx(cmt)),tmas/cmp,evel
     >  ,dene(ic),teme(ic),temi(ic)
     >  ,den0BG(ic,1),eng0BG(ic,1),dengBG(ic,1),enggBG(ic,1)
     >  ,mrct,(irct(i),zvlp/srct(i),i=1,nrct)
!
 680  format(i7,2x,i4,3f10.5,2x,a4,2x,i2,2i5
     > ,2x,a10,2x,f6.2,1pe10.2
     > ,1p3e10.2
     > ,1p4e10.2
     > ,2x,i3,2x,8(i4,1pe10.2))
!
      i=index(cmt,"folw_end")
      i=i+index(cmt,"mole_end")
      i=i+index(cmt,"folw_ion")
      if(i>0) close(nf)
      return
      end
