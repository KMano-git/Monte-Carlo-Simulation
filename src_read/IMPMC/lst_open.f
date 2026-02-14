!***********************************************************************
      subroutine lst_open(typ)
!***********************************************************************
!)
!) "pcon"    : particle conservation  cpu, fpuff, fexh
!) "scat"    : positions of test particles
!) "trak"    : trace of some partciles for debug
!)
!) "pcnT"    : MPI_reduce  cpu, fpuff, fexh   T(total)
!) "flxb"    : MPI_reduce  time evolution of fbc_flxI/O
!) "wrad"    : MPI_reduce  time evolution of radiation in each section
!)
!)----------------------------------------------------------------------
      use cimcom, only : i6_mscat, i6_pcnt, i6_pcon, i6_scat, i6_trak
     >    , i6_wrad
      use cunit,  only : lmspe, lmype, n6, ufl6
      implicit none

! modified 1/1 lines organize local variables and include files by kamata 2021/07/04
!ik   character :: typ*(*)
      character, intent(in) :: typ*(*)

      character :: dsn*20

!::dsn
      write(dsn,'("lst_",a,"_",i5.3)') trim(typ), lmype

      select case(typ)
!-----
      case("trak")    ! Note not "trac"
        i6_trak = 0
        if( ufl6 /= "/dev/null" ) then
          i6_trak = 161
          open(unit=i6_trak,file=dsn)
          write(i6_trak,'(1x, 2x,a, 2x,a, 3x,a, 2x,a, 5x,a, 4x,a,
     >      6x,a, 8x,a, 10x,a, 6x,a,
     >      4x,a, 2x,a, 2x,a, 1x,a, 4x,a, 1x,a, 2x,a, 2x,a,
     >      8x,a, 8x,a, 6x,a, 4x,a, 9x,a, 8x,a)')
     >      "comt", "ksp", "lp", "ich", "ist", "ip",
     >      "ptmx", "ptim", "tt", "pdtm",
     >      "ic", "ix", "iy", "reg","ln", "ien", "ko", "is",
     >      "rr", "zz", "ro", "wght", "vz", "Eng"
        endif
      write(n6,602) "lst_open", "i6_trak  =",i6_trak
 602    format(2x,a,2x,a,i5)

!-----
      case("pcon")
        i6_pcon = 0
        if( ufl6 /= "/dev/null" ) then
          i6_pcon = 162
          open(unit=i6_pcon,file=dsn)
          write(i6_pcon,'(3x,a,3x,a,7x,a,5x,a,6x,15(a,5x))')
     >     "itm", "stimz", "npmax", "cpuM",
     >     "pc_tst", "pc_ten", "pc_tbl", "pc_emt", "pc_exh", "pc_pmp",
     >     "tot_rw", "rw_odv", "rw_sol", "rw_idv", "rw_opv", "rw_ipv",
     >     "rw_edg", "rw_man", "rw_vac"
!xx      write(i6,'(1x,i6,1pe12.4,1p2e10.2,1p6e11.3,1p9e11.3)')
        endif
      write(n6,602) "lst_open", "i6_pcon  =",i6_pcon

!-----
      case("scat")
      i6_scat = 0
      if( ufl6 /= "/dev/null" ) then
        i6_scat = 163
        open(unit=i6_scat,file=dsn)
!::move these statement into sub. lst_scat  (print more than once)
!xx      write(i6_scat,'(/2x,"particle information  tnow =",1pe12.4,
!xx     >  "  npmax =",i6)') tnow, npmax
!xx      write(i6_scat,'(4x,"ip",2x,"tt",11x,"ic",4x,"ix",3x,"iy",2x,"ien",3x,
!xx     > "il",3x,"is",3x,"rr",8x,"zz",7x,"wgt",6x,"Vz",9x,"Vl",9x,"Ez")')
      endif
      write(n6,602) "lst_open", "i6_scat  =",i6_scat

!-----
      case("pcnT")
        i6_pcnT = 0
        if( lmype == lmspe ) then
          i6_pcnT = 164
          open(unit=i6_pcnT,file=dsn)
          write(i6_pcnT,'(3x,a,3x,a,7x,a,5x,a,6x,15(a,5x))')
     >     "itm", "stimz", "npmax", "cpuM",
     >     "pc_tst", "pc_ten", "pc_tbl", "pc_emt", "pc_exh", "pc_pmp",
     >     "tot_rw", "rw_odv", "rw_sol", "rw_idv", "rw_opv", "rw_ipv",
     >     "rw_edg", "rw_man", "rw_vac"
!xx      write(i6,'(1x,i6,1pe12.4,1p2e10.2,1p6e11.3,1p9e11.3)')
        endif
      write(n6,602) "lst_open", "i6_pcnT  =",i6_pcnT

!-----
!$$$      case("flxb")
!$$$        i6_flxb = 0
!$$$        if( lmype == lmspe ) then
!$$$          i6_flxb = 165
!$$$          open(unit=i6_flxb,file=dsn)
!$$$        endif
!$$$      write(n6,602) "lst_open", "i6_flxb  =",i6_flxb

!-----
      case("wrad")
        i6_wrad = 0
        if( lmype == lmspe ) then
          i6_wrad = 166
          open(unit=i6_wrad,file=dsn)
          write(i6_wrad,'(4x,a,2x,a,8x,a,6x,9(a,6x))')
     >    "itm", "stimz", "npmax", "wr_tot", "wr_odv", "wr_sol",
     >    "wr_idv", "wr_opv", "wr_ipv", "wr_edg", "wr_cor", "wr_vac"
        endif
      write(n6,602) "lst_open", "i6_wrad  =",i6_wrad

!-----
!$$$      case("hflxb")
!$$$        i6_hflxb = 0
!$$$        if( lmype == lmspe ) then
!$$$          i6_hflxb = 167
!$$$          open(unit=i6_hflxb,file=dsn)
!$$$          write(i6_hflxb,'(4x,a,2x,a,7x,a,9x,8(a,7x))')
!$$$     >      "itm","stimz","rho","fI_15","fO_15",
!$$$     >      "fI_16","fO_16","fI_17","fO_17"
!$$$!)SC  "nI_16", "nI_18"
!$$$        endif
!$$$      write(n6,602) "lst_open", "i6_hflxb =",i6_hflxb

!-----
      case("mscat")  ! dsn = "lst_mscat_tm0050"
      i6_mscat = 0
      if( lmype == lmspe ) then
        i6_mscat = 168
      endif
!xx   not open/close
      write(n6,602) "lst_open", "i6_mscat  =",i6_mscat

!-----
      case default
      end select

      return
      end
