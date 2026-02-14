! remade TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
! added integrated calculation with SONIC by kamata 2020/11/16
! sent / received data to file
      subroutine dbg_top_sr( kk, ktyp, gtim, gmion, n6 )
      use topics_mod, only : gcnm, gnro, gro, groh, gvr, ndmk
      implicit none
! arguments
      integer, intent(in) :: gmion, kk, ktyp, n6
      real(8), intent(in) :: gtim
! gmion : number of ions
! gtim  : data time
! ktyp  : data type
!         = 1 : TOPICS profile data, = 2 : SONIC edge data
! kk    : send / recieve flag
!         = 1 : recieve, = 2 : send
! n6    : unit number for log

! local variables
      character :: fnam(2,2)*10 =
     >    reshape( (/ 'Dprfrcv', 'Dedgrcv', 'Dprfsnd', 'Dedgsnd' /)
     >           , (/ 2, 2 /) ) ! fnam(ktyp,kk)
      integer   :: io = 1, is, k, n
      real(8)   :: ztot(ndmk) = 0.0_8

! file open
      open( io, file=fnam(ktyp,kk), form='formatted'
     >    , position='append' )

      write(io,'(/a,es15.7)') 'time = ', gtim
      do is = 0, gmion
        write(io,'(a,i2)') 'is = ', is
        write(io,'(3x,"n",3x,"ro",6x,"roh",5x,20(a5,6x))')
     >      (gcnm(k),k=1,ndmk)

        ztot(9:ndmk) = 0.0_8
        do n = 1, gnro
          write(io,'(1x,i3,2f8.4,20es11.3)')
     >        n, gro(n), groh(n), (gvr(n,is,k),k=1,ndmk)
          do k = 9, ndmk-1
            ztot(k) = ztot(k) + gvr(n,is,k)
          enddo  ! loop(k)
        enddo  ! loop(n)
        write(io,'(1x,3x,"Tot",13x,20es11.3)') (ztot(k),k=1,ndmk)
      enddo  ! loop(is)

! file close
      close( io )

      return
      end
