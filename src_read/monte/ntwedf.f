!**********************************************************************
      subroutine ntwedf(iw, zweit, ene, ns)
!**********************************************************************
!
!    neutral energy distribution function on the wall
!
!      wedf(nene, nw, ns)
!        ns: 1:atom 2:mol
!        nw: No of masured position. iw is specied by iw_edf(nw)
!        nene: No of energy slot
!
!      wedf is included vmwork. see  cntwcn
!      MPI_reduce and smoothing will be done as vmwork
!
!      weight will be multiplied in ntwflx
!
!----------------------------------------------------------------------
      use cntcom,  only : weit
      use cntwedf, only : d_ene, enemax, enemin, iw_edf, nwmx, wedf
     > ,d_ene_line, enemax_line, enemin_line, wedf_line
      implicit none

      integer, intent(in) :: iw, ns
      real(8), intent(in) :: zweit,ene

      integer nw, lsv, nene
      real(8) zene

      lsv = 0
      do nw = 1,nwmx
        if(iw_edf(nw).eq.0) exit
        if(iw.eq.iw_edf(nw))then
          lsv = 1
          exit
        endif
      enddo
      if(lsv.ne.1) return

! log version
      zene = ene
      zene = dmax1(zene,enemin)
      zene = dmin1(zene,enemax)
      zene = dlog10(zene)
      nene = int((zene-dlog10(enemin))/d_ene) + 1
      wedf(nene,nw,ns) = wedf(nene,nw,ns) + weit

! liner version
      zene = ene
      zene = dmax1(zene,enemin_line)
      zene = dmin1(zene,enemax_line)
      nene = int((zene-enemin_line)/d_ene_line) + 1
      wedf_line(nene,nw,ns) = wedf_line(nene,nw,ns) + weit  
      return
      end

!**********************************************************************
      subroutine lst_wedf
!**********************************************************************
!
!    neutral energy distribution function on the wall
!
!----------------------------------------------------------------------
      use cntwedf, only : d_ene, enemin, iw_edf, nene_mx, nsmx, nwmx
     > , xwedf, deg_iw
     > , d_ene_line, enemin_line, nene_mx_line, xwedf_line
      use cntwfl,  only : xare, xdeg
      use cunit,   only : lmspe, lmype, n6
      use csize,   only : ndwp
      use cntwfl,  only : iywl
      use cntcom,  only : npew
      implicit none

      real(8) wk1, wk2, wk3
      integer nene, nw, ns, iw
      integer nf

      character(:),allocatable :: filename
      integer i
      integer, allocatable :: indexListInner(:)
      integer, allocatable :: indexListOuter(:)
      integer :: inner_upperbound = 0, inner_lowerbound = 0
      integer :: outer_upperbound = 0, outer_lowerbound = 0
      logical is_log

      if(lmype.ne.lmspe ) return
      write(n6,'(/2x,"*** lst_wedf *** ")')

      do ns = 1, nsmx
       do nw = 1, nwmx
        iw = iw_edf(nw)
        if(iw .le. 0) exit 
        wk1 = 1.0d0/xare(iw)
        ! for log version
        do nene = 1, nene_mx
         xwedf(nene,nw,ns) = xwedf(nene,nw,ns) * wk1
        enddo
        ! for liner version
        do nene = 1, nene_mx_line
         xwedf_line(nene,nw,ns) = xwedf_line(nene,nw,ns) * wk1
        enddo
       enddo
      enddo

      nf = 21
      
!     output from deg_iw (defined by inc/cntwedf)
      open(nf,file="Wntedf.txt")
      write(nf,'(a1,3x,3(5x,a3),10(2x,a3,i3,f6.1))') ! this assume deg_iw = 5
     >     "#","=<E","E<","Eav",
     > ("atm",iw_edf(nw),xdeg(iw_edf(nw)),nw=1,size(deg_iw)),
     > ("mol",iw_edf(nw),xdeg(iw_edf(nw)),nw=1,size(deg_iw))

      do nene = 1, nene_mx
        wk1 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)-1.0d0))
        wk2 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)      ))
        wk3 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)-0.5d0))
        write(nf,'(i4, 3f8.1, 14e14.6)') 
     >   nene, wk1,wk2,wk3, 
     >  (xwedf(nene,nw,1),nw=1,size(deg_iw)),
     >  (xwedf(nene,nw,2),nw=1,size(deg_iw))
      enddo
      close(nf)

!     The position of diverter is found by "iywn" of exe/XXX/XXX/wxdr_XX/Wntl.txt.
      do iw = 2, ndwp
       if(iywl(iw-1) == 0 .and. .not. iywl(iw) == 0
     >   .and. inner_upperbound == 0) then
        inner_upperbound = iw
       elseif(iywl(iw) == 0 .and. .not. iywl(iw-1) == 0
     >   .and. inner_lowerbound == 0 
     >   .and. .not. inner_upperbound == 0) then
        inner_lowerbound = iw-1
       endif

       if(iywl(iw-1) == 0 .and. .not. iywl(iw) == 0
     >   .and. .not. inner_lowerbound == 0
     >   .and. outer_upperbound == 0
     >   .and. iw > npew(3)) then
         outer_upperbound = iw
       elseif(iywl(iw) == 0 .and. .not. iywl(iw-1) == 0
     >   .and. outer_lowerbound == 0
     >   .and. .not. outer_upperbound == 0
     >   .and. iw > npew(3)) then
         outer_lowerbound = iw-1
       endif
      enddo

      do i = 1, nwmx
        if(iw_edf(i) .le. 0) exit
        if(inner_upperbound <= iw_edf(i) .and. 
     >    iw_edf(i) <= inner_lowerbound) then
            if(allocated(indexListInner)) then
              indexListInner = [indexListInner,iw_edf(i)]
            else
              allocate(indexListInner(1))
              indexListInner(1) = iw_edf(i)
            endif
        elseif(outer_upperbound <= iw_edf(i) .and.
     >    iw_edf(i) <= outer_lowerbound) then
          if(allocated(indexListOuter)) then
            indexListOuter = [indexListOuter,iw_edf(i)]
          else
            allocate(indexListOuter(1))
            indexListOuter(1) = iw_edf(i)
          endif
        endif
      enddo
!
!:: log version
      is_log = .true.
      filename = "Wntedf_dia.txt"
      call WriteWntedf(filename,"atm",indexListInner,is_log)
      filename = "Wntedf_dim.txt"
      call WriteWntedf(filename,"mol",indexListInner,is_log)
      filename = "Wntedf_doa.txt"
      call WriteWntedf(filename,"atm",indexListOuter,is_log)
      filename = "Wntedf_dom.txt"
      call WriteWntedf(filename,"mol",indexListOuter,is_log)
!
!:: liner version
      is_log = .false.
      filename = "Wntedf_dia_deng.txt"
      call WriteWntedf(filename,"atm",indexListInner,is_log)
      filename = "Wntedf_dim_deng.txt"
      call WriteWntedf(filename,"mol",indexListInner,is_log)
      filename = "Wntedf_doa_deng.txt"
      call WriteWntedf(filename,"atm",indexListOuter,is_log)
      filename = "Wntedf_dom_deng.txt"
      call WriteWntedf(filename,"mol",indexListOuter,is_log)
!
!**********************************************************************
      contains
      subroutine WriteWntedf(filename,isAtm,indexList,is_log)
!**********************************************************************
      implicit none
!
!:: arguments
      character(:), intent(in), allocatable :: filename
      character(3), intent(in) :: isAtm
      integer, intent(in) :: indexList(:)
      logical, intent(in) :: is_log
!
!:: local variables
      character(:), allocatable :: myWriteFormatTitle
      character(:), allocatable :: myWriteFormatResult
      character(20)::nwmxCharacter
      integer atm_mole_index
      
      open(nf,file=filename)

      write(nwmxCharacter,"(I0)") size(indexList)
      allocate(myWriteFormatTitle,
     > source = "(a1,3x,7x,a4,6x,a5,8x,a3,"
     > //trim(nwmxCharacter)
     > //"(2x,a4,i3.3,a1,4x))")
!     22/9/15 yamamoto, uncomment if you output degree (1/4)
!     > //"(2x,a3,a1,i3.3,a1,f7.3))"      

      if(isAtm=="atm") then
       write(nf,myWriteFormatTitle)
     >     "#","Elow","Ehigh","Eav",
     >    ("atm(",indexList(nw),")",
     >     nw=1,size(indexList))
!     22/9/15 yamamoto, uncomment if you output degree (2/4)
!     >     xdeg(indexList(nw))),nw=1,size(indexList))
      else
       write(nf,myWriteFormatTitle)
     >     "#","Elow","Ehigh","Eav",
     >    ("mol(",indexList(nw),")",   
     >     nw=1,size(indexList))
!     22/9/15 yamamoto, uncomment if you output degree (3/4)
!     >     xdeg(indexList(nw))),nw=1,size(indexList))
      endif

      write(nwmxCharacter,"(I0)") size(indexList)+4
      allocate(myWriteFormatResult,
     > source = "(i4,3f11.3,"
     > //trim(nwmxCharacter)
     > //"e14.6)")
!     22/9/15 yamamoto, uncomment if you output degree (4/4)
!      myWriteFormatResult = myWriteFormatResult//"e17.6)"
!
      if(isAtm=="atm") then
        atm_mole_index = 1 ! for atom
      else
        atm_mole_index = 2 ! for molecular
      endif
!
      if(is_log) then
        do nene = 1, nene_mx
          wk1 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)-1.0d0))
          wk2 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)      ))
          wk3 = 10.0d0**(dlog10(enemin) + d_ene*(dble(nene)-0.5d0))
          write(nf,myWriteFormatResult)
     >    nene, wk1,wk2,wk3
     >    ,(xwedf(nene,getIndex(indexList(nw)),atm_mole_index)
     >    ,nw=1,size(indexList))
        enddo
!
      else
        do nene = 1, nene_mx_line
          wk1 = enemin_line + d_ene_line*(dble(nene)-1.0d0)
          wk2 = enemin_line + d_ene_line*(dble(nene))
          wk3 = enemin_line + d_ene_line*(dble(nene)-0.5d0)
          write(nf,myWriteFormatResult)
     >    nene, wk1,wk2,wk3
     >    ,(xwedf_line(nene,getIndex(indexList(nw)),atm_mole_index)
     >    ,nw=1,size(indexList))
        enddo
      endif
!
      close(nf)
      end subroutine

!**********************************************************************
      integer function getIndex(numIndexList)
!**********************************************************************
      integer numIndexList,j
      getIndex = 0
      do j = 1, nwmx
        if(iw_edf(j)==numIndexList) then
            getIndex = j
        endif
      enddo 
      end function
      end subroutine


!**********************************************************************
      subroutine ntwedf_init
!**********************************************************************i
!     neutral energy distribution function on the wall
!---
      use cntwedf, only : iw_edf, nwmx, deg_iw
      use cunit,   only : n6
      use cntcom,  only : npsw, npew
      use cntwfl,  only : xdeg
      use csize,   only : ndwp
      use cntwfl,  only : iywl
      use mod_externalgrid, only : use_exdata
      implicit none

      real*8 wrk, wrk0
      integer iw, iw0, i
      integer iws, iwe

      integer :: inner_upperbound = 0, inner_lowerbound = 0
      integer :: outer_upperbound = 0, outer_lowerbound = 0
      integer inner_iw_quant, outer_iw_quant

      write(n6,'(/2x,"*** ntwedf_init ***")')

!     get index consist of deg_iw. deg_iw is defined by inc/cntwedf
      iws = npsw(1)
      iwe = npew(1)
      do i = 1, size(deg_iw)
        if(deg_iw(i).lt.0.0d0)exit
        wrk  = 999.0d0
        wrk0 = 999.0d0
        iw0 = 0
        do iw = iws, iwe
          if(.not.use_exdata .and. iw==iwe) cycle
          wrk = abs(xdeg(iw)-deg_iw(i))
          if(wrk.gt.180.0d0) wrk = abs(xdeg(iw)-360.0d0-deg_iw(i))
          if(wrk.lt.wrk0)then
            wrk0=wrk
            iw0=iw
          endif
        enddo ! iw
        iw_edf(i) = iw0
      enddo
      write(n6,'(2x,a2,a6,3x,a3,a7)') "i","deg_iw","iw_edf","degree"
      do i = 1, size(deg_iw)
        write(n6,'(2x,i2,x,f5.1,a3,i3,f7.2)') 
     >      i, deg_iw(i), "=>", iw_edf(i), xdeg(iw_edf(i))
      enddo
      
!     The position of diverter is found by "iywn" of exe/XXX/XXX/wxdr_XX/Wntl.txt.
      do iw = 2, ndwp
       if(iywl(iw-1) == 0 .and. .not. iywl(iw) == 0
     >   .and. inner_upperbound == 0) then
        inner_upperbound = iw
       elseif(iywl(iw) == 0 .and. .not. iywl(iw-1) == 0
     >   .and. inner_lowerbound == 0 
     >   .and. .not. inner_upperbound == 0) then
        inner_lowerbound = iw-1
       endif

       if(iywl(iw-1) == 0 .and. .not. iywl(iw) == 0
     >   .and. .not. inner_lowerbound == 0
     >   .and. outer_upperbound == 0
     >   .and. iw > npew(3)) then
         outer_upperbound = iw
       elseif(iywl(iw) == 0 .and. .not. iywl(iw-1) == 0
     >   .and. outer_lowerbound == 0
     >   .and. .not. outer_upperbound == 0
     >   .and. iw > npew(3)) then
         outer_lowerbound = iw-1
       endif
      enddo

      inner_iw_quant = inner_lowerbound - inner_upperbound + 1
      outer_iw_quant = outer_lowerbound - outer_upperbound + 1
      if(inner_iw_quant + outer_iw_quant + size(deg_iw)> nwmx) then
        call wexit("ntwedf_init","nwmx less than divertor element")
      endif

      do i = 1, inner_iw_quant
        iw_edf(i+size(deg_iw)) = i + inner_upperbound - 1
      enddo    
      do i = 1, outer_iw_quant 
        iw_edf(i+size(deg_iw)+inner_iw_quant) = i + outer_upperbound - 1
      enddo

      write(n6,'(2x,a2,2x,a6)') "i","iw_edf"
      do i = 1+size(deg_iw), nwmx
        write(n6,'(2x,i4,x,i5)') i, iw_edf(i)
      enddo

      return
      end
