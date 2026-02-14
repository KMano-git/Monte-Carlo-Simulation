!**********************************************************************
      subroutine imp_mpsum(typ)
!**********************************************************************
!)
!)## "debz"
!)  MPI_Reduce  denZ, temZ in /cimcoom_11/
!)
!)   parameter ( nwkmp_im1 = 10+(ndmis+1)*(ndmc+1)*2+(ndmc+1) )
!)   real(8) :: sflux, swtot, sptcl, stimp, sitmp, stimz, sitmz, sdtmz
!)   real(8) :: scpum, sdumy
!)   real(8) :: denZ(0:ndmis,0:ndmc), temZ(0:ndmis,0:ndmc), wsct(0:ndmc)
!)
!)   Not sumup : sflux stimp(time_soldor) sitmp(itim_soldor)
!)      stimz(time_IMPMC) sitmz(time_IMPMC) sdtmz(dtimZ)
!)
!)## "pcon"
!)  MPI_Reduce  hpcon in /cimcom_12/
!)
!)   parameter ( nwkmp_im3 = 5 + 10*5 )
!)   real(8) :: hmype, htout, hstep, hcpum, hpmax
!)   real(8) :: hgrpp(10), hgrpw(10) ! conditon
!)   real(8) :: hregp(10), hregw(10) ! region
!)   real(8) :: hpcon(10)            ! conservation
!)
!)## "flxb"
!)  MPI_Reduce  fbc_flxI/O in /cimfbc_2/ of inc file(cimfbc)
!)
!)   parameter ( nwkmp_im5 = (ndmis+1)*5*3 )
!)   real(8) :: fbc_flxI(0:ndmis,5), fbc_flxO(0:ndmis,5)
!)   real(8) :: fbc_denc(0:ndmis,5)
!)
!)---------------------------------------------------------------------
      use cimcom, only : ns => nwkmp_im1, nsy => nwkmp_im1y
     >    , ns3 => nwkmp_im3, scpum, sdtmz, sflux, sitmp, sitmz, stimp
     >    , stimz
      use cimfls, only : fls_flxi, fls_flxo, lfls
      use cunit,  only : lmspe, lmype, lnope, mywld, n6
      use mpi!,    only : mpi_real8, mpi_reduce, mpi_sum
      use csize,  only : ndmis
      implicit none

      character, intent(in) :: typ*(*)

!::mpi variables "denz"
      real(8), dimension(ns) :: wwork, vwork
!::mpi variables "simp"
      real(8), dimension(nsy) :: wwork2, vwork2
!::mpi variables "pcon"
      real(8), dimension(ns3) :: wwork3, vwork3

!::mpi variables "flxb"
      integer, parameter :: nwk_fls = (ndmis+1)*5*3 !cimedg_3 fls_flxI, fls_flxO, fls_dnz
      real(8), dimension(nwk_fls)  :: wfls, vfls
      integer ierr
      real(8)  svtmp(10)
      integer :: i, iz

!::No need MPI_Reduce
      if( lnope == 1 ) then
        write(n6,'(2x,"  Not passed MPI_Reduce [",a,"] due to ",
     >     "lnope =",i2)') trim(typ), lnope
        return
      endif

      vwork(1:ns)     = 0.0_8
      vwork2(1:nsy)   = 0.0_8
      vwork3(1:ns3)   = 0.0_8
      vfls(1:nwk_fls) = 0.0_8

!::variables
      select case(typ)

!----------------------------------------------------------------------
!::MPI_Reduce in imp_mpsum("denz")  cimcom_11 size (sflux,wsct)
!----------------------------------------------------------------------
      case("denz")
        call setvwork( 1, 2, wwork, ns )
        svtmp(1:10) = wwork(1:10)
        write(n6,'(2x,"1PE wwork(1:1) =",1p10e12.4)') wwork(1:10)
        call MPI_Reduce( wwork,  vwork,  ns, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
        call setvwork( 1, 1, vwork, ns )

!::uncountable variabes
        if( lmype == lmspe ) then
          sflux = svtmp(1)
          stimp = svtmp(4)
          sitmp = svtmp(5)
          stimz = svtmp(6)
          sitmz = svtmp(7)
          sdtmz = svtmp(8)
          scpum = svtmp(9)/dfloat(lnope)
          write(n6,'(2x,"TPE wwork(1:1) =",1p10e12.4)') vwork(1:9)
        endif

!----------------------------------------------------------------------
!::MPI_Reduce in imp_mpsum("simp")  cimcom_11y
!----------------------------------------------------------------------
      case("simp")
        call setvworky( 1, 2, wwork2, nsy )
        call MPI_Reduce( wwork2, vwork2, nsy, MPI_REAL8, MPI_SUM,
     >                lmspe, mywld, ierr )
        call setvworky( 1, 1, vwork2, nsy )
        write(n6,*) typ

!----------------------------------------------------------------------
!::MPI_Reduce in imp_mpsum("pcon")  cimcom_11 size (sflux,wsct)
!----------------------------------------------------------------------
      case("pcon")
        call setvwork3( 1, 2, wwork3, ns3 )
        svtmp(1:5) = wwork3(1:5)
        write(n6,'(2x,"1PE wwork3(16:25) =",1p10e12.4)') wwork3(16:25)
        call MPI_Reduce( wwork3, vwork3, ns3, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
        call setvwork3( 1, 1, vwork3, ns3 )
        write(n6,'(2x,"TPE wwork3(16:25) =",1p10e12.4)') vwork3(16:25)

!----------------------------------------------------------------------
!::MPI_Reduce in imp_mpsum("flxb")  cimcom_11 size (sflux,wsct)
!----------------------------------------------------------------------
      case("flos")
        if( lfls.eq.1 ) then
          iz = 10
          write(n6,'(2x,a,2x,i3,2x,1p10e11.3)')
     >      "fls_1PE",iz, fls_flxI(iz,1:5), fls_flxO(iz,1:5)
          call setvfls( 2, wfls, nwk_fls )
          call MPI_Reduce( wfls, vfls, nwk_fls, MPI_REAL8, MPI_SUM,
     >                   lmspe, mywld, ierr )
          call setvfls( 1, vfls, nwk_fls )
        endif
!-----
      case default
        call wexit("imp_sumup","invalid  typ = "//typ)

      end select
      return
      end
