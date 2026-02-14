!**********************************************************************
      subroutine imp_rstfl(kact)
!**********************************************************************
!)
!)    kact = 0 : set job_number  rst_drW = "../wxdr_24"
!)         = 1 : write restart file for IMPMC
!)         = 2 : read  restart file for IMPMC
!)
!)    data of ion & neutral particle data  not dtyp data
!)
!)::ion
!)  integer npmax, npemt, jpmax
!)  integer ir(ndmp), is(ndmp), il(ndmp), ivs(ndmp), ivd(ndmp)
!)  integer itg(ndmp), igc(ndmp), igs(ndmp), igi(ndmp)
!)  integer iml(ndmp), ien(ndmp), igtno(ndmp), igtem(ndmp)
!)  integer ntags(ndmp)
!)  real(8) pemt(ndmp)
!)  real(8) tt(ndmp), rr(ndmp), zz(ndmp), fi(ndmp)
!)  real(8) vr(ndmp), vz(ndmp), vv(ndmp), wght(ndmp), wstp(ndmp)
!)  real(8) vvr(ndmp), vvz(ndmp), v(ndmp)
!)  character  cimcom2_emrk*1
!)
!)::neutral
!)  integer ndmp2; parameter (ndmp2=24000)
!)  real*8  sx0(ndmp2), sy0(ndmp2), sz0(ndmp2)
!)  real*8  svx(ndmp2), svy(ndmp2), svz(ndmp2)
!)  real*8  srn(ndmp2), sint(ndmp2), stm(ndmp2), stb(ndmp2)
!)
!)::rstfile(0)  rst_drW(wxdr_25)  rst_drR(wxdr_10)
!)  GL know   GM  donot know  rst_drW rst_drR ==> ALL PE know
!)
!)  Note. inquire
!)   [../IMPV5/imfile.f]    lex = T
!)   [  ../IMPV5/imfile.f]  lex = T
!)    ls "  ../IMPV5_imfile.f"  No such file or directory
!)
!)--------------------------------------------------------------------
      use cimcom,      only : cftrz, cftwz, fi, ien, igi, igs
     >    , igtem, igtno, il, iml, ir, is, itg, ivd, ivs
     >    , npmax, ntags, pemt, rr, tt, v, vr, vv, vvr, vvz, vz, wght
     >    , wstp, zz, ndmp
      use cimntl,      only : sint, srn, stb, stm, svx, svy, svz, sx0
     >    , sy0, sz0
      use cunit,       only : lmype, n6
      use mod_mpicomm, only : nope_grp, nrnk_grp, nwld_grp
      use mpi!,         only : mpi_barrier
      implicit none

!::argument
      integer, intent(in) :: kact

!::local variables
      integer :: nft, mx, ip, i, ipe, ierr, nrnk_grpX, ip0
      character(80) :: dsn
      logical :: lex
!::dummy
      integer :: igc(ndmp) = 0

      select case(kact)

!::restart file (see iminpt.f)
      case(0)
        write(n6,'(/a,i1,2x,"cftwZ = [",a,"]",2x,"cftrZ = [",a,"]")')
     >   "RSTMC", kact, trim(cftwZ), trim(cftrZ)
        return

!::write
      case(1)
        nft = 21
        dsn = cftwZ
        write(n6,'(a,i1,2x,"wrt-dsn = ",a)') "RSTMC", kact, trim(dsn)
        
        do i = 1, nope_grp
          ipe = i-1

          if( nrnk_grp == ipe ) then  ! err 69
            open(unit=nft,file=dsn,action="write",form="unformatted",
     >        position="append")

            write(nft) nrnk_grp, npmax
            mx = npmax
            ntags(1:mx) = 1 ! force ntags=1 for BK_ON -> OFF capability
            write(nft)
     >       ir(1:mx), is(1:mx), il(1:mx), ivs(1:mx), ivd(1:mx),
     >       itg(1:mx), igc(1:mx), igs(1:mx), igi(1:mx),
     >       iml(1:mx), ien(1:mx), igtno(1:mx), igtem(1:mx),
     >       ntags(1:mx), pemt(1:mx),
     >       tt(1:mx), rr(1:mx), zz(1:mx), fi(1:mx),
     >       vr(1:mx), vz(1:mx), vv(1:mx), wght(1:mx), wstp(1:mx),
     >       vvr(1:mx), vvz(1:mx), v(1:mx),
     >       sx0(1:mx), sy0(1:mx), sz0(1:mx),
     >       svx(1:mx), svy(1:mx), svz(1:mx),
     >       srn(1:mx), sint(1:mx), stm(1:mx), stb(1:mx)
            close(nft)
          endif
          call MPI_Barrier(nwld_grp,ierr)  ! important
        enddo

!::read
      case(2)
        nft = 21
        dsn = cftrZ
        write(n6,'(a,i1,2x,"red-dsn = ",a)')
     >    "RSTMC", kact, trim(dsn)

        inquire(file=trim(dsn), exist=lex)
        if( dsn(1:1) == " " ) lex = .FALSE.
        if( lex ) then
          write(n6,'(a,i1,2x,"execute RESTART cal. from dsn ",a)')
     >    "RSTMC", kact, "["//trim(dsn)//"]"
        else
          write(n6,'(a,i1,2x,"execute INITIAL cal. nonexistent dsn "
     >     ,a)') "RSTMC", kact, "["//trim(dsn)//"]"
          return
        endif
        call flush(n6)

        open(unit=nft,file=dsn,action="read",form="unformatted")

        do i = 1, nope_grp
          read(nft) nrnk_grpX, npmax
          mx = npmax

          read(nft)
     >      ir(1:mx), is(1:mx), il(1:mx), ivs(1:mx), ivd(1:mx),
     >      itg(1:mx), igc(1:mx), igs(1:mx), igi(1:mx),
     >      iml(1:mx), ien(1:mx), igtno(1:mx), igtem(1:mx),
     >      ntags(1:mx), pemt(1:mx),
     >      tt(1:mx), rr(1:mx), zz(1:mx), fi(1:mx),
     >      vr(1:mx), vz(1:mx), vv(1:mx), wght(1:mx), wstp(1:mx),
     >      vvr(1:mx), vvz(1:mx), v(1:mx),
     >      sx0(1:mx), sy0(1:mx), sz0(1:mx),
     >      svx(1:mx), svy(1:mx), svz(1:mx),
     >      srn(1:mx), sint(1:mx), stm(1:mx), stb(1:mx)
          if( nrnk_grp == nrnk_grpX ) exit
        enddo
        close(nft)

      case default
        write(n6,'(2x,"invalid kact =",i3,"  in sub. imp_rstfil")')
     >     kact
        call wexit("imp_rstfl","invalid value of kact")
      end select

!::debug write
      ip0 = 0  !  3   sample particle number
      if( ip0 > 0 ) then
        write(n6,'(a,i1,2x,i6,i7,i4,i4,1p4e12.3)') "RSTMC", kact,
     >  lmype,npmax,ip0,is(ip0),tt(ip0),rr(ip0),zz(ip0),wght(ip0)
      endif

      write(n6,'(2x,"kact =",i2,"  nrnk_grp =",i5,"  npmax =",i6)')
     >  kact, nrnk_grp, npmax
      write(n6,'(6x,a)') "ip    ic  is  ien  ntags   pemt        tt"//
     >  "          rr          zz          vz          vv"//
     >  "          wght        wstp        stb         stm"
      do ip = 1, npmax
        if( ip < 4 .or. ip > npmax-3 ) then
          write(n6,'(2x,i6,i6,i4,i4,i4,1p10e12.3)')
     >      ip,ir(ip),is(ip),ien(ip),ntags(ip),
     >      pemt(ip),tt(ip),rr(ip),zz(ip),
     >      vz(ip),vv(ip),wght(ip),wstp(ip), stb(ip), stm(ip)
        endif
      enddo

!::set tt = 0.0
      if( kact == 2 ) then
        write(n6,'(a,i1,2x,"imp_reset  set tt,sx0,sy0,sz0")')
     >    "RSTMC", kact
        call imp_pack ! SY
        call imp_reset
        write(n6,'(a,i1,2x,a)') "RSTMC"
     >    ,kact,"finish to initial setting"
      endif

      return
      end
