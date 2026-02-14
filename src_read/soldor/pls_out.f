!***********************************************************************
      subroutine pls_out
!***********************************************************************
      use cplcom, only : cftw, nftw
      use csonic, only : itim, kdsk, khst, lntl
      use cunit,  only : n6
      implicit none
!
      if( kdsk.eq.1 ) then
        write(n6,'(2x)')
        if( lntl.eq.3 ) call ntoutp(n6,'conv')
        call ntoutp(n6,'tots')
      endif
!
      if( itim.le.5 .or. mod(itim,100).eq.0 )  call plmflx
!
      call pmonit       ! moniter
!
      if( kdsk.eq.1 ) then
         call pldisk(nftw,1,cftw)
         call pst_rprf
      endif
!
      if( khst.eq.1 ) call plhist(2)
!
      if(mod(itim,1000).eq.0)call ploutp(3,25) !KH for FEC25
!
!xxx  call plrscu       ! rescue
!
!::selected time
!xx      call plevol(itim)
!
      return
      end
