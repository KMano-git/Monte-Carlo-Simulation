!> \brief sight line view for variaous properties.
!>  accumulate values property on the mesh along sight lines
!>  and draw the accumulated values on the display.
!>
!> \param[in] prop_name   property name
!> \param[in] prop_size   size of data array
!> \param[in] prop        data array
!>                        kcn=0: normal end, kcn &gt; 0 otherwise 
!> \param[in] scal        scaling factor
!> \param[out] kcn        end status, kcn=0; normal end, kcn &gt; 0 otherwize
!**********************************************************************
      subroutine chord_view(prop_name, prop_size, prop)
!**********************************************************************
      
      implicit none
      
      character(len=*), intent(in)  :: prop_name
      integer(kind=4),  intent(in)  :: prop_size
      real(kind=8),     intent(in)  :: prop(prop_size)
!
      real(kind=8)  :: scal = 1.0d0
      integer       :: kcn

!>  \brief capacity of sight lines
      integer, parameter :: ndslmx = 100
  
!>  \brief the start and end points of sight lines.
!>         the i-th sight line starts at the point (xst(i), yst(i)) and
!>                             ends   at the point (xen(i), yen(i))
      real*8 :: xst(ndslmx), yst(ndslmx), xen(ndslmx), yen(ndslmx)
!> \brief channel index of sight lines
      real(kind=4) :: channel(ndslmx)
  
!> \brief the number of valid sight lines
      integer :: nsl
  
!> \brief accumulated value on sight lines
      real(kind=4), allocatable :: accv(:)
  
  
! \brief the capacity of the intermedium points on a sight line
      integer, parameter :: ndim=400
! \brief a list of the coordinates of the intermedium points on a sight line
      real(kind=8) :: posx(ndim), posy(ndim)
! \brief a list of cell indecies
      integer :: ipos(ndim)
! \brief valid number of array elements
      integer :: ncol
  
!> the unit No. for a output file
      integer :: nft
  
!> \brief string buffer for building a file name
      character(len=80) :: dsn
  
!> \brief user input string
      character(len=80) :: cinp
!> \brief the number of characters of the user input string
      integer :: mji
!> \brief flag that indicate whether drawing or not
      logical :: draw_flg
  
      integer :: i
      integer :: ipen
      real*8 :: integ_v
      real*4 :: xmx
      character(len=3) :: l_no

      logical lex

      interface
          subroutine read_chord(unitno, capacity,
     >                                xst, yst, xen, yen,
     >                                nlines, error)
              integer, intent(in)  :: unitno
              integer, intent(in)  :: capacity
              real*8, intent(out)  :: xst(capacity)
              real*8, intent(out)  :: yst(capacity)
              real*8, intent(out)  :: xen(capacity)
              real*8, intent(out)  :: yen(capacity)
              integer, intent(out) :: nlines
              integer, intent(out), optional :: error
          end subroutine read_chord
      end interface

      kcn=0

      ! wall
      call set_pwall
  
      ! preparation
      !!call set_HAL(nsl, xha1, yha1, xha2, yha2, ndslmx)
      nft=1
      open(unit=nft, file="inpfig")
      call read_chord(nft, ndslmx, xst, yst, xen, yen, nsl)
      close(unit=nft)
  
      allocate(accv(1:nsl))
  
      ! calcurate the accumulated values and output them to file
      nft=21
      open(unit=nft, file=trim(adjustl(prop_name))//'_prof.txt')
      ! write title line
      write(unit=nft,
     >    fmt='(2x,"i",2x,"ncl",3x,"xp1",9x,"yp1",9x,"xp2",9x,"yp2",
     >    9x,a3,9x)') trim(adjustl(prop_name))
      do i=1, nsl
          ! trace sight line on the mesh
          call pth_line(xst(i), yst(i), xen(i), yen(i), ndim,
     >                  ncol, posx, posy, ipos)
!D         write(n6,*) 'traced a sight line'
          ! accumulate a property value along a sight line
          call pth_intg(ncol, posx, posy, ipos, prop, integ_v)
!          accv(i)=integ_v*scal
          accv(i)=integ_v
!D         write(n6,*) 'integrated along a sight line'
  
          ! write title line
          write(unit=nft, fmt='(2x, i3, i5, 1p5e12.3)')
     >        i, ncol, posx(1), posy(1), posx(ncol), posy(ncol), integ_v
      end do
      close(unit=nft)
!D     write(n6,*) '**** write integrated data ****'
  
      forall(i=1:nsl) channel(i)=real(i)
      call gdsiz(3.0, 2.0, 20.0, 15.0)   ! set a drawing region
      xmx = (int(nsl/10.0d0) + 0.999d0)*10.0d0
      call gdaxs(1, 0.0, xmx, "channel")
      call gdays(1, 0.0, 99.0, trim(adjustl(prop_name))//" intensity")
      call gdplt(1, nsl, channel, accv, 1, trim(adjustl(prop_name)))
      call gdpag(1)
!     write(n6,*) '**** draw integrated data ****'

      deallocate(accv)
  
      ! path line
      nft=21
      do i=1, nsl
          write(dsn, '(a,"_pline_",i2.2)') trim(adjustl(prop_name)), i
          if (nft > 0) then
              open(unit=nft, file=trim(adjustl(dsn)))
          endif
  
          call pth_line(xst(i), yst(i), xen(i), yen(i), ndim,
     >                  ncol, posx, posy, ipos)
          call pth_list(ncol, posx, posy, ipos, trim(adjustl(prop_name))
     >                , prop, nft)
  
          if (nft > 0) then
              close(nft)
          end if
      end do
  
  
      ! plot sight lines
      call gdget(">> need plot ? (y/<rtn>:n ==> ", cinp, mji)
      draw_flg=.false.
      if (mji > 0 .and. cinp(1:1) == 'y') draw_flg=.true.
      if (draw_flg) then

          nft = 21
          call nopen(nft,"inpfig","text",dsn,lex)
          call gpsta(nft)
          close(nft)
          call gtfgmx

          call msplot(0, 1)
          do i=1, nsl
              ipen=2
              l_no=' '
              if (mod(i,10) == 0) ipen=3
              if (mod(i,10) == 1) write(l_no, '(i3)') i
              call pth_line(xst(i), yst(i), xen(i), yen(i), ndim,
     >                      ncol, posx, posy, ipos)
              call gdpltd(1, ncol, posx, posy, ipen,trim(adjustl(l_no)))
          end do
          call gdpag(1)
      end if
  
      return
      end subroutine chord_view 
