!**********************************************************************
      subroutine imp_file(typ)
!**********************************************************************
!)
!)   "pcon"   : partcile conservation for 1PE
!)   "scat"   : position of test partciles
!)   "trak"   : trace of a partcile
!)   "pcnT"   : particle conservation for TPE
!)   "wrad"   : radiation in each section
!)   "hflxb"  : history of flxb
!)
!)   How to edit the code
!)    cimcom      real/common  i6_wrad
!)    imp_file    input data, case("wrad")
!)    lst_open    header
!)    imp_print   print data
!)    imp_file    close
!)
!)--------------------------------------------------------------------
      use cimcom, only : i6_flxb, i6_hflxb, i6_mscat, i6_pcnT, i6_pcon
     >    , i6_scat, i6_trak, i6_wrad, ip_trak
      use cunit,  only : n5, n6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   character :: typ*(*)
      character, intent(in) :: typ*(*)

      integer :: ierr

      namelist /uimlst/ ip_trak, i6_trak, i6_scat,
     >   i6_pcon, i6_pcnT, i6_flxb, i6_wrad, i6_hflxb, i6_mscat

      select case(typ)

      case("setunt")
         ip_trak = 0
         i6_trak = 0
         i6_scat = 0
         i6_pcon = 0
         i6_pcnT = 0
         i6_flxb = 0
         i6_wrad = 0
         i6_hflxb = 0
         i6_mscat = 0
        read(n5,uimlst,iostat=ierr)
        if(ierr==0)then
        if( ip_trak == 0 ) i6_trak = 0
        write(n6,'(/2x,"*** imp_file ***  setunt")')
        write(n6,'(4x,"ip_trak =",i6,"  i6_trak =",i5,
     >   "  i6_scat =",i5,"  i6_pcon ~",i5,"  i6_pcnT =",i5,
     >   "  i6_flxb =",i5,"  i6_wrad =",i5)')
     >    ip_trak, i6_trak, i6_scat, i6_pcon, i6_pcnT,
     >    i6_flxb, i6_wrad
        write(n6,'(4x,"i6_hflxb=",i6,"  i6_mscat=",i6)')
     >    i6_hflxb, i6_mscat
        endif ! ierr

      case("open")
        write(n6,'(2x,"*** imp_file *** : open")')
        if( i6_pcon  > 0 ) call lst_open("pcon")
        if( i6_scat  > 0 ) call lst_open("scat")
        if( i6_trak  > 0 ) call lst_open("trak")
        if( i6_pcnT  > 0 ) call lst_open("pcnT")
!        if( i6_flxb  > 0 ) call lst_open("flxb")
        if( i6_wrad  > 0 ) call lst_open("wrad")
!        if( i6_hflxb > 0 ) call lst_open("hflxb")
!xx     if( i6_mscat > 0 ) call lst_open("mscat")  ! Not open/close

      case("close")
        write(n6,'(2x,"*** imp_file *** : close")')
        if( i6_pcon  > 0 ) close(i6_pcon)
        if( i6_scat  > 0 ) close(i6_scat)
        if( i6_trak  > 0 ) close(i6_trak)
        if( i6_pcnT  > 0 ) close(i6_pcnT)
!        if( i6_flxb  > 0 ) close(i6_flxb)
        if( i6_wrad  > 0 ) close(i6_wrad)
!        if( i6_hflxb > 0 ) close(i6_hflxb)
!xx     if( i6_mscat > 0 ) close(i6_mscat)

      case default
      end select

      return
      end
