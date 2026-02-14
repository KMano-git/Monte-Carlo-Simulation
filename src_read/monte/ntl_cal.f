!**********************************************************************
      subroutine ntl_cal
!**********************************************************************
      use cntcom,  only : lntmd
      use cntmnt,  only : vcsrc, vnsrc
      use cputim,  only : cntcal
      use csonic,  only : itim
      use cunit,   only : n6
      use mpi!,     only : mpi_wtime
      use ntpfctl, only : out_pfctl
      implicit none
!
!::local varaiables
      integer i
      real(8) :: cptm0, cptm1
!
!::cpu
      cptm0 = MPI_WTIME()
!
      if( lntmd.eq.0 ) goto 900
!
      call ntmont
!
      if(mod(itim,100).eq.0 ) call out_pfctl
!
!::output
      write(n6,'(2x,"<cpu-time>  ",i2,10(2x,a6,1pe11.2))')
     >  vnsrc,(vcsrc(i),cntcal(i),i=1,vnsrc)
      write(n6,'(30("="),"  END  ntl_cal  itim =",i6)') itim
!
 900  continue
      cptm1 = MPI_WTIME()
      cntcal(0) = cptm1-cptm0
!
      return
      end
