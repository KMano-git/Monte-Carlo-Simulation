!**********************************************************************
      subroutine imwflx(ii)
!**********************************************************************
      use cimcom, only : denflxztowall, eneflxztowall, ismax, ndis
     >    , sflux, weflx, wpflx, wtemt
      use cimden, only : csput
      use cntcom, only : npew, npsw
      use cntwfl, only : xare
      use csize,  only : ndwp
      use cunit,  only : lmspe, lmype, mype, mywld, n6
      use mod_shexe, only : impmc_model
      use mpi!,    only : mpi_real8, mpi_reduce, mpi_sum
      use mod_externalgrid, only : use_exdata
      implicit none
!
      integer, intent(in) :: ii
!
!::mpi variables
      integer :: ierr
!
!::local variables
      integer :: iw

! Toku, Impurity Flux onto Wall
      integer :: nArrayLength
      real*8  denFlxZreceiver(0:ndis,ndwp),eneFlxZreceiver(0:ndis,ndwp)
!
      integer :: iws, iwe, nw
      real(8) :: fac
      real(8) :: totf(0:ndis),tote(0:ndis)
!
      denFlxZreceiver=0.0d0
      eneFlxZreceiver=0.0d0
      totf=0.0d0
      tote=0.0d0
!----------------------------------------------------------------------
!::sumup
!----------------------------------------------------------------------
!
![TOKU] Broadcast impurity flux data from local master PE to all (as well as in the ntmont)
!      call MPI_Bcast( vwork2, ns2, MPI_REAL8, lmspe, mywld, ierr )
      nArrayLength = ( ndis + 1)*ndwp
      call MPI_Reduce( denFlxZtoWall, denFlxZreceiver, nArrayLength,
     >                 MPI_REAL8, MPI_SUM, lmspe, mywld, ierr )
      call MPI_Reduce( eneFlxZtoWall, eneFlxZreceiver, nArrayLength,
     >                 MPI_REAL8, MPI_SUM, lmspe, mywld, ierr )
!
      if( lmype.ne.lmspe ) return
!
!----------------------------------------------------------------------
      write(n6,*)
      if( impmc_model == 0 ) then
        write(n6,"(2x,'*** imwflx  ***',i2,2x,a5,2x,'mype=',i5,
     &           '  sflux=',e16.8,'  wtemt=',e16.8)")
     &      ii,csput(ii),mype,sflux,wtemt
      else
        write(n6,"('Flux to Wall: ',i2,2x,a5,2x,'mype=',i5,
     &           '  sflux=',e16.8,'  wtemt=',e16.8)")
     &      ii,csput(ii),mype,sflux,wtemt
      endif

      if( impmc_model == 0 ) fac = sflux / wtemt
      do nw = 1, 4
        iws = npsw(nw)
        iwe = npew(nw)
        if( impmc_model == 0 ) then
          do iw = iws, iwe
            if(.not.use_exdata .and. iw==iwe) cycle
            wpflx(0:ismax,iw,ii) = 
     >       denFlxZreceiver(0:ismax,iw)*fac/xare(iw)
            weflx(0:ismax,iw,ii) = 
     >       eneFlxZreceiver(0:ismax,iw)*fac/xare(iw) !* cev
            totf(0:ismax) = totf(0:ismax)+wpflx(0:ismax,iw,ii)*xare(iw)
            tote(0:ismax) = tote(0:ismax)+weflx(0:ismax,iw,ii)*xare(iw)
          enddo
        else
!SY 2020.03.11    fac = sflux / wtemt: Already considered in imwrfi/n
          do iw = iws, iwe
            if(.not.use_exdata .and. iw==iwe) cycle
            wpflx(1:ismax,iw,ii) = denFlxZreceiver(1:ismax,iw)/xare(iw)
            weflx(1:ismax,iw,ii) = eneFlxZreceiver(1:ismax,iw)/xare(iw) !* cev
          enddo
        endif
      enddo

      if( impmc_model == 0 ) then
        write(n6,'(4x,"total partle flux (iz=1-5) = ",75e12.4)')
     >    totf(0:5)
        write(n6,'(4x,"total energy flux (iz=1-5) = ",75e12.4)')
     >    tote(0:5)
        write(n6,*)
      endif
!
      return
      end
