!**********************************************************************
      subroutine ntl_end
!**********************************************************************
      use cntcom, only : cftnw, nftnw
      use cunit,  only : n6
      implicit none
!
      integer itemp
!
      write(n6,'(/2x,"%%% ntl_end %%%")')
!
      itemp = nftnw
      nftnw = 39
      if( cftnw(1:1).eq." " ) nftnw = 0
!
      write(n6,'(4x,"nftnw =",2i5,"  cftnw =",a)') ! KCSFUJI
     >     itemp,nftnw,trim(cftnw)
!
      if( nftnw.gt.0 ) then
      call ntdisk(nftnw,1,cftnw)
      endif
!
!xxx      call ntpump(4)
      call ntpump(3)   ! KH150213
!
!KH move to ntl_ini      call ntwdeg
      call ntwflx
      call plwflx    ! when limp = 0  KSFUJI
      call lswflx("Wnflx.txt")
      call lst_wedf
!
      return
      end
