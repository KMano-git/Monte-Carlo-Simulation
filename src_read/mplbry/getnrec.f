! added for use with ncopy.f by kamata 2022/05/29
! get number of records in file
      subroutine getnrec( fnm, fnw, nr )
      implicit none
! arguments
      character, intent(in)  :: fnm*(*), fnw*(*)
      integer,   intent(out) :: nr
! fnm : target file name
! fnw : work   file name
! nr  : number of records in fnm

! local variables
      integer :: ios, nft= 21
      real(4) :: wtm = 0.1_4
! ios : i/o status
! nft : unit number for result file
! wtm : sleep time

! executing wc command
!     print *,'system start ',trim( fnw )
      call system( 'wc ' // trim( fnm ) // ' > ' // trim( fnw ) )
!     print *,'system end ',trim( fnw )

!- read results
      do
!       call sleep( wtm )
        open( nft, file=fnw, form='formatted', status='old'
     >    , iostat=ios )
!       print *,ios,trim( fnw )
        if( ios /= 0 ) cycle
        read(nft,*,iostat=ios) nr
!       print *,'nr = ',nr,trim( fnw )
        if( ios /= 0 ) then
          close( nft )
          cycle
        endif
        close( nft, status='delete'  )
        exit
      enddo

      return
      end
