!***********************************************************************
      subroutine iminit(emttyp)
!***********************************************************************
      use cimcom, only : cdi, cdi0, dtemt, emflx, emptl, emwgt
     >    , iml, irmax, itg, ndmp, npemt, npmax, sflux, sptyc, sptyi
     >    , swtot, tfbz
      use csize,  only : ndmc
      use cunit,  only : lmype, lnope, n6
      use mod_shexe, only : impmc_model
      implicit none
!
!::argument
      character, intent(in) :: emttyp*(*)
!
!::local variables
      real*8     zero
      integer    ip, ic
      integer    iemtyp, nc
      character(2) nc_char
! function
      real(8)    rnlast
!
      write(n6,'(/2x,"*** iminit ***")')
      write(n6,'(2x,"first  random in iminit  ran =",0pf12.9)')
     >   rnlast(0)
!
!::tag-name
      do ip = 1, npmax
        itg(ip) = lmype*npmax + ip
      enddo
!
!::impurity generation model
      if( emttyp.eq."Cchp" ) then
         iemtyp = 11
         call imgchem("idp,prv,odp")

      elseif( emttyp.eq."Cchw" ) then
         iemtyp = 12
         call imgchem("sol,puf,vol")

      elseif( emttyp.eq."Cchd" ) then
         iemtyp = 13
         call imgchem("prv")

      elseif( emttyp.eq."Cchm" ) then
         iemtyp = 14
         call imgchem("idp,prv,odp,sol,puf,vol")

      elseif( emttyp.eq."Cphy" ) then
         iemtyp = 15
         call imgphys("phy")

      elseif( (emttyp(1:2).eq."Zm") ) then
         iemtyp = 16
         read(emttyp(3:len(emttyp)),*) nc
         write(nc_char,"(i0)") nc
         call imgzspt( "zsp"//trim(nc_char), nc )

      !****************!
      ! for old-version
      !****************!
      elseif( emttyp.eq."Zphy" .or. emttyp.eq."Zspt") then
         iemtyp = 16
         call imgzspt( "zsp1", 1 )
      !****************!
      elseif( emttyp.eq."Zph2" .or. emttyp.eq."Zsp2") then
         iemtyp = 16
         call imgzspt( "zsp2", 3 )
      !****************!

      elseif( emttyp.eq."Arpf" ) then
         iemtyp = 21
         call imgpuff

      elseif( emttyp.eq."Arbk" .or. emttyp.eq."Bflw") then
         iemtyp = 22
         call imgbkgt

      elseif( emttyp.eq."Ipnt" ) then
         iemtyp = 31
         call imgpont("---")

      elseif( emttyp.eq."Cedg" ) then
         iemtyp = 41
         call imgcedg("---")

      else
         call wexit("iminit","No found sputtering "//emttyp)
      endif
!
!::sputtering type for disk data /cimcom_11b/
      sptyi = iemtyp
      sptyc = emttyp
!
!::no sputtering
      if( tfbz.le.0.0d0 ) return
!
!::cross section of ionization of neutral
      zero = 1.0d-40
      cdi0(0:ndmc) = zero
      do ic = 1, irmax
        cdi0(ic) = cdi(0,ic)
      enddo
!
!::CD4 ==> C+ (simple model)  flag iml(ip)
      iml(1:ndmp) = 0
      if( iemtyp.ge.11 .and. iemtyp.le.14 ) then
        call set_sgiz(1)
        iml(1:ndmp) = 1
      endif
!
!::initial loading of imp ions (imsput) ==> imgene
!     see sub. imemit

      if( impmc_model == 1 ) then
!::normalization factor  for 1PE (Not lnope) ??? (IMPV5)
!::Note   pemt(ip) = emptl in imemt_wall.f
        emflx = sflux
        emwgt = npemt*lnope
        emptl = emflx*dtemt/emwgt
!DBG Yamoto
        swtot = emwgt

        write(n6,'(2x,"## def.", a ,2x,
     >    "pemt(ip) = emptl ",2x,"emflx =",1pe12.4,
     >    "  emwgt =",1pe12.4,"  emptl =",1pe12.4)')
     >     emttyp, emflx, emwgt, emptl
      endif

      return
      end