! modified generalization of derived data type by kamata 2020/10/31
! 1) ierr to lerr in dtyp_addv 
! 2) ntyp to nddt ( number of derived data type )
! original : sonicV4/MPMD5/WFL/mod_dtypdef.f

!  In: integer, Ch: character, R4: real(4), R8:real(8)
!  0d: scalar,  1d: 1 dimension, 2d: 2 dimension, 3d, 4d, 5d
!  dtyp_addv_0dIn, 1dIn, 2dIn, 3dIn
!  dtyp_addv_0dR8, 1dR8, 2dR8, 3dR8
!  dtyp_addv_0dCh,
!  dtyp_addv_0dR4,
!
!   How to use dtyp_addv routine
!     dtyp_addv(vne)  Array
!     dtyp_addv(te0)  Scalar
!
!   Note 
!   integer(kind=MPI_Address_Kind) :: addr ! for call MPI_GET_Address
!----------------------------------------------------------------------
      module mod_dtypdef
!#ifdef MPI
      use mpi
      use mod_sizedef, only : lnnam, ndtyp, ndvio
      use mod_mpicomm, only : m6
      implicit none

! modified 5/4 lines generalization of derived data type by kamata 2020/10/31
!ik   character(lnnam) :: tydnam
!ik   character(lnnam) :: typrog
!ik   integer  :: tyitag
!ik   integer  :: tyierr
!ik   integer(kind=MPI_Address_Kind) :: tyadr0
      integer :: nddt = 0 ! number of derived data type
      character(lnnam) :: typ_dnam(ndtyp)
      character(lnnam) :: typ_prog(ndtyp)
      integer          :: typ_itag(ndtyp)
      integer(kind=MPI_Address_Kind) :: displace(ndvio)

      integer :: blocklen(ndvio)
! deleted 1 line generalization of derived data type by kamata 2020/10/31
!ik   integer :: dispint4(ndvio)
      integer :: typelist(ndvio)
      integer, save :: nvio
! added 2 lines generalization of derived data type by kamata 2020/10/31
      integer(kind=MPI_Address_Kind) :: addr
      integer :: lerr

      interface dtyp_addv
        module procedure 
     >    dtyp_0dIn, dtyp_1dIn, dtyp_2dIn, dtyp_3dIn,
     >    dtyp_0dR8, dtyp_1dR8, dtyp_2dR8, dtyp_3dR8,
     >    dtyp_0dCh, dtyp_1dCh, dtyp_2dCh, dtyp_3dCh,
     >    dtyp_0dR4, dtyp_1dR4, dtyp_2dR4, dtyp_3dR4
      end interface

      contains
!---------------------------------------------------------------------
      subroutine dtyp_init( dnam, cprg )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      character(*), intent(in) :: dnam, cprg

      nvio  = 0
      blocklen(1:ndvio) = 0
      displace(1:ndvio) = 0
      typelist(1:ndvio) = MPI_REAL8

! modified 2/3 lines generalization of derived data type by kamata 2020/10/31
!ik   tydnam = dnam
!ik   typrog = cprg
      nddt = nddt + 1
! added 6 lines check nddt by kamata 2022/08/08, yamamoto 22/10/07
      if( nddt > ndtyp ) then
       write(m6,
     >    '(2x,"** dtyp_init **  nddt (",i5," ) > ndtyp (",i5," )")')
     >  nddt, ndtyp
        call wf_exit( 'dtyp_init', 'not created new type struct' )
      endif
      typ_dnam(nddt) = dnam
      typ_prog(nddt) = cprg

      return
      end subroutine dtyp_init

!---------------------------------------------------------------------
      subroutine dtyp_0dIn( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      integer, intent(in) :: var

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = 1
      typelist(nvio) = MPI_INTEGER
      displace(nvio) = addr
      return
      end subroutine dtyp_0dIn

!---------------------------------------------------------------------
      subroutine dtyp_1dIn( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      integer, intent(in) :: var(:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blockLen(nvio) = size( var )
      typelist(nvio) = MPI_INTEGER
      displace(nvio) = addr
      return
      end subroutine dtyp_1dIn

!---------------------------------------------------------------------
      subroutine dtyp_2dIn( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      integer, intent(in) :: var(:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_INTEGER
      displace(nvio) = addr
      return
      end subroutine dtyp_2dIn

!---------------------------------------------------------------------
      subroutine dtyp_3dIn( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      integer, intent(in) :: var(:,:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_INTEGER
      displace(nvio) = addr
      return
      end subroutine dtyp_3dIn

!---------------------------------------------------------------------
      subroutine dtyp_0dR8( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(8), intent(in) :: var

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = 1
      typelist(nvio) = MPI_REAL8
      displace(nvio) = addr
      return
      end subroutine dtyp_0dR8

!---------------------------------------------------------------------
      subroutine dtyp_1dR8( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(8), intent(in) :: var(:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL8
      displace(nvio) = addr
      return
      end subroutine dtyp_1dR8

!---------------------------------------------------------------------
      subroutine dtyp_2dR8( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(8), intent(in) :: var(:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL8
      displace(nvio) = addr
      return
      end subroutine dtyp_2dR8

!---------------------------------------------------------------------
      subroutine dtyp_3dR8( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(8), intent(in) :: var(:,:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL8
      displace(nvio) = addr
      return
      end subroutine dtyp_3dR8

!---------------------------------------------------------------------
      subroutine dtyp_0dCh( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      character(*), intent(in) :: var

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = len( var )
      typelist(nvio) = MPI_CHARACTER
      displace(nvio) = addr
      return
      end subroutine dtyp_0dCh

!---------------------------------------------------------------------
      subroutine dtyp_1dCh( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      character(*), intent(in) :: var(:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = len( var ) * size( var )
      typelist(nvio) = MPI_CHARACTER
      displace(nvio) = addr
      return
      end subroutine dtyp_1dCh

!---------------------------------------------------------------------
      subroutine dtyp_2dCh( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      character(*), intent(in) :: var(:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = len( var ) * size( var )
      typelist(nvio) = MPI_CHARACTER
      displace(nvio) = addr
      return
      end subroutine dtyp_2dCh

!---------------------------------------------------------------------
      subroutine dtyp_3dCh( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      character(*), intent(in) :: var(:,:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = len( var ) * size( var )
      typelist(nvio) = MPI_CHARACTER
      displace(nvio) = addr
      return
      end subroutine dtyp_3dCh

!---------------------------------------------------------------------
      subroutine dtyp_0dR4( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(4), intent(in) :: var

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = 1
      typelist(nvio) = MPI_REAL4
      displace(nvio) = addr
      return
      end subroutine dtyp_0dR4

!---------------------------------------------------------------------
      subroutine dtyp_1dR4( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(4), intent(in) :: var(:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL4
      displace(nvio) = addr
      return
      end subroutine dtyp_1dR4

!---------------------------------------------------------------------
      subroutine dtyp_2dR4( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(4), intent(in) :: var(:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL4
      displace(nvio) = addr
      return
      end subroutine dtyp_2dR4

!---------------------------------------------------------------------
      subroutine dtyp_3dR4( var )
!---------------------------------------------------------------------
!     use mpi
      implicit none
! arguments
      real(4),  intent(in) :: var(:,:,:)

      nvio = nvio + 1
      if( nvio > ndvio ) call dtyp_error

      call MPI_GET_Address( var, addr, lerr )
      blocklen(nvio) = size( var )
      typelist(nvio) = MPI_REAL4
      displace(nvio) = addr
      return
      end subroutine dtyp_3dR4

!---------------------------------------------------------------------
! modified 1/1 lines generalization of derived data type by kamata 2020/10/31
!ik   subroutine dtyp_term( ierr )
      subroutine dtyp_term( kerr )
!---------------------------------------------------------------------
      use mod_mpicomm, only : m6
      implicit none
! arguments
! modified 1/1 lines generalization of derived data type by kamata 2020/10/31
!ik   integer, intent(OUT) :: ierr
      integer, intent(out) :: kerr
! local variables
      integer(kind=MPI_Address_Kind) :: adrs0
      integer :: ktype
! added 1 line generalization of derived data type by kamata 2020/10/31
      integer    dispint4(ndvio)
! added 2 lines get the code number from the run-time parameters by kamata 2021/08/18
      integer    i
      integer(kind=MPI_Address_Kind) :: adrsr, i4max = 2**31

      adrs0 = displace(1)
! modified 2/1 lines generalization of derived data type by kamata 2020/10/31
!ik   displace(1:nvio) = displace(1:nvio) - adrs0
!ik   dispint4(1:nvio) = displace(1:nvio)
      dispint4(1:nvio) = displace(1:nvio) - adrs0
!      do i = 1, nvio
!        adrsr = displace(i) - adrs0
!        if( i4max > adrsr .and. adrsr >=  -i4max ) then
!          dispint4(i) = adrsr
!        else
!          write(m6,
!     >      '(2x,"** dtyp_term **  [",a,"]  displace =",i20,
!     >        "  relative address =",i20)')
!     >      trim( typ_dnam(nddt) ), displace(i), adrsr
!          call wf_exit( 'dtyp_term', 'not created type struct' )
!        endif
!      enddo

! modified 2/3 lines generalization of derived data type by kamata 2020/10/31
!ik   call MPI_Type_Struct( nvio,blockLen,dispint4,typelist,ktype,ierr )
!ik   call MPI_Type_Commit( ktype,ierr )
      call MPI_Type_Struct( nvio, blocklen, dispint4, typelist, ktype
     >                    , kerr )
      call MPI_Type_Commit( ktype, kerr )

! modified 3/6 lines generalization of derived data type by kamata 2020/10/31
!ik   if( ierr /= 0 )  write(m6,
!ik  >  '(2x,"*** dtyp_term ***  [",a,"]  ktype =",i12,"  ierr =",i4)')
!ik  > trim(tydnam), ktype, ierr
      if( kerr /= 0 )  then
        write(m6,
     >    '(2x,"** dtyp_term **  [",a,"]  ktype =",i12,"  kerr =",i4)')
     >    trim( typ_dnam(nddt) ), ktype, kerr
        call wf_exit( 'dtyp_term', 'not created type struct' )
      endif

! modified 3/1 lines generalization of derived data type by kamata 2020/10/31
!ik   tyitag = ktype
!ik   tyierr = ierr
!ik   tyadr0 = adrs0
      typ_itag(nddt) = ktype

      return
      end subroutine dtyp_term
!
!
      subroutine dtyp_error
      implicit none

      write(m6,'(/2x,"Dimension Error at sub. dtyp_addv  ",a,2x,2i5)')
! modified 3/3 lines generalization of derived data type by kamata 2020/10/31
!ik  > trim(tydnam), nvio, ndvio
!ik   call wf_exit("dtyp_addv",
!ik  >  "Dimension Error nvio > ndvio  "//trim(tydnam) )
     > trim( typ_dnam(nddt) ), nvio, ndvio
      call wf_exit( 'dtyp_addv', 'Dimension Error nvio > ndvio  ' //
     >    trim( typ_dnam(nddt) ) )
      end subroutine dtyp_error
!#endif

      end module mod_dtypdef
