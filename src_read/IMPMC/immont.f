!***********************************************************************
      subroutine immont
!***********************************************************************
!
!   sputtering type   nwspt*4
!      "Cchp" : chemical sputtering from target plate  iemtyp = 11
!      "Cchw" : chemical sputtering from wall          iemtyp = 12
!      "Cchd" : physical sputtering from target plate  iemtyp = 13
!      "Cphy" : self sputtering from target plate      iemtyp = 14
!      "ZmXX" : physical sputtering from target plate  iemtyp = 16
!      "Arpf" : Ar puff                                iemtyp = 21
!      "Arbk" : Ar back flow from gate                 iemtyp = 22
!      "Ipnt" : point source for check                 iemtyp = 31
!
!      typspt = "Cchp, -Cchw, Cphy" => "Cchp, Cphy"
!       iemtyp => sptyi (i:integer)
!                                      see  sub. iminit & iminpt
!-----------------------------------------------------------------------
      use cimcom, only : amz, azmas, denflxztowall, dtimz, e0bkz, e0emtz
     >    , e0inz, e0tmz, eneflxztowall, fwcal, fwstp, irmax, is, ismax
     >    , itimz, ldist, lhst, lpcm, lprf, ltmax, npmax
     >    , ns => nwkmps_im1, nsy => nwkmps_im1y, sflux, sptyc, swtot
     >    , tfbz, timez, v0emtz, weflx, wpflx
      use cimctl, only : icalz, nprg
      use cimden, only : csput, eipot, fsput, ncmx, npsput, nsput
     >    , nwspt, nzmx, wtsput
      use cimfls, only : lfls
      use cimi60, only : i60_imcndv, i60_imhist
      use cphcns, only : cev
      use cplwrd, only : wfac
      use csize,  only : ndmc, ndmis
      use csonic, only : itim
      use cunit,  only : lmspe, lnope, mywld, n6
      use mpi!,    only : mpi_real8, mpi_reduce, mpi_sum
      implicit none
!
!::local variables
      real*8  wsflx
      integer lt, i, k, ii, ip, kend
      real(8), dimension(ns) :: wwork, vwork
      real(8), dimension(nsy) :: wworky, vworky

      integer, parameter :: nwk_fls = (ndmis+1)*5*3 !cimedg_2 fls_flxI, fls_flxO, fls_dnz
      real(8), dimension(nwk_fls)  :: wfls, vfls
      integer    ierr
!
      call ran_init
!
!
      if(wfac.le.0.0d0) wfac = 1.0d0
      if(nprg.ne.1) wfac = 1.0d0
!
!::clear  sputtered flux   KSFUJI   referr in imdisk
      fsput(1:5) = 0.0d0
!:: clear
      wpflx(:,:,:) = 0.0d0
      weflx(:,:,:) = 0.0d0
      vwork(1:ns)     = 0.0_8
      vworky(1:nsy)   = 0.0_8
      vfls(1:nwk_fls) = 0.0_8
!
!----------------------------------------------------------------------
!::generation mechanism
!----------------------------------------------------------------------
      ii = 0
      do k = 1, nsput
        nwspt  = csput(k)
        npmax  = npsput(k)
!
        e0emtZ = e0inZ
        if( nwspt.eq."Arpf" ) e0emtZ = e0inZ
        if( nwspt.eq."Arbk" ) e0emtZ = e0bkZ
        if( nwspt.eq."Cphy" ) e0emtZ = e0tmZ
        v0emtZ = sqrt(2.0d0*e0emtZ*cev/amz)

        write(n6,'(/30("="),"  immont",2x,a,"  itim =",i6,
     >   "  icalZ =",i4,
     >   "  npmax =",i8,"  e0emtZ =",1pe12.3,"  v0emtZ =",1pe12.3)')
     >   nwspt, itim, icalZ, npmax, e0emtZ, v0emtZ
!
!----------------------------------------------------------------------
!::ionization process of sputtered neutral impurity
!----------------------------------------------------------------------
        call imclear("particle")
        call imclear("score")
        call imclear("floss")
! toku: initialize
        denFlxZtoWall = 0.0d0
        eneFlxZtoWall = 0.0d0
!
        call iminit(nwspt)
        ii = ii + 1
        sptyc = csput(ii)
        fsput(ii) = sflux
!
        write(n6,'(2x,"impurity generation ",i2,2x,a,2x,1pe12.3,
     >    "  azmas =",0pf7.3,"  ltmax =",i6,"  npmax =",i6)')
     >   ii, csput(ii), fsput(ii), azmas, ltmax, npmax
        write(n6,'(2x,"eipot =",10f10.3)') (eipot(i),i=0,ismax)
!
!----------------------------------------------------------------------
!::trace of impurity ion & neutral    dtimZ*dfloat(lout)
!----------------------------------------------------------------------
        if( tfbz.gt.0.0d0 ) then ! tfbz<0 => no sputtering
!::No define sputtered particle
          do ip = 1, npmax
            is(ip) = -1
          enddo
!
          timeZ = 0.0d0
          itimZ = 0
          lpcm  = 0
!
          do lt = 1, ltmax
            write(n6, '(/2x,"==== immont lt loop ==== lt/ltmax = "
     >       , i5,i5)') lt, ltmax
            lpcm = lt
            call imrqout(lt,lprf,lhst)
            call imstep
!-----
            if( lprf.eq.1 ) then
              if( ldist.eq.1 ) call imdist(lt,dtimZ)
            endif
            if( lhst.eq.1 ) then
              call imcndv
              call imwcon(0)
            endif
!-----
            call imhist(lpcm,kend)
            if( kend.eq.1 ) exit

          enddo !lt
!
!::reserve sflux  before mpi_reduce
          wsflx = sflux
          call dbg_denz("1PE")
!
!----------------------------------------------------------------------
!::MPIsumup
!----------------------------------------------------------------------
!::[MPI_Reduce in immont]  cimcom_11 size (sflux,wsct)   10/05/27
          if( lnope.gt.1 ) then
            call setvwork( 0, 2, wwork, ns )
            call MPI_Reduce( wwork, vwork, ns, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
            call setvwork( 0, 1, vwork, ns )
            ! Collecting additional impurity information: Yamoto
            call setvworky( 0, 2, wworky, nsy )
            call MPI_Reduce( wworky, vworky, nsy, MPI_REAL8, MPI_SUM,
     >                 lmspe, mywld, ierr )
            call setvworky( 0, 1, vworky, nsy )
          else
            write(n6,'(2x,"  Not passed MPI_Reduce & MPI_Bcast",
     >         " due to lnope =",i2)') lnope
          endif

          if(lfls.eq.1)then
            if( lnope.ge.1 ) then
              call setvfls( 2, wfls, nwk_fls )
              call MPI_Reduce( wfls, vfls, nwk_fls, MPI_REAL8
     >                   ,MPI_SUM,lmspe, mywld, ierr )
              call setvfls( 1, vfls, nwk_fls )
            endif
          endif
!
          sflux = wsflx
          call dbg_denz("TPE")
!
!
!::KSFUJI (use dnz  Not use tdnz)
          i60_imcndv = n6   ! output imcndv & imhist
          i60_imhist = n6
          call imcndv
          call dbg_denz("imcndv")
          call imhist(lpcm,kend)
          call dbg_denz("imhist")
          call imwcon(1)    ! 1:keep sum-up data!
          call dbg_denz("imwcon")

!::impurity density
          call imdens
          call dbg_denz("imdens")
          if(lfls.eq.1) call lstflos("flux")
          wtsput(ii) = swtot
!
!::Flux to Wall
          call imwflx(ii)
!
          write(n6,'(1x,a,6x,2x,a,2x,i7,2x,a,2x,i10,i5,
     >      2x,a,2x,1pe12.4,i7,2x,a,2x,0p2f8.4,i3)')
     >      "imhistC  "//sptyc, "itim",itim, "lpcm/icalZ", lpcm, 
     >       icalZ,"sflux/npmax", sflux, npmax,"  fwcal/fwstp/kend",
     >       fwcal, fwstp, kend
!
        endif ! tfbz > 0
!
!----------------------------------------------------------------------
!::prorog
!----------------------------------------------------------------------
        call imsumry(2)
        call imsave(sptyc,"w",k)
        call imsave(sptyc,"y",k)
        call dbg_denz("imsave")
        call dbg_vzpz("imsave")
        write(n6,'(2x,"### dbg_denz")')
      enddo ! k: sptyp loop
!
      nsput = ii
      ncmx  = irmax
      nzmx  = ismax
!
      write(n6,'(/2x,"generation mechanism   nsput =",i2)') nsput
      write(n6,'(2x," i",2x,"csput",2x,"wtsput",2x,"fsput")')
      do i = 1, nsput
        write(n6,'(2x,i2,2x,a,2x,1p2e12.4)') 
     >   i,csput(i),wtsput(i),fsput(i)
      enddo
      write(n6,'(30("="),"  immont  end")')

!
      return
      end