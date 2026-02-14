! added integrated calculation with TOPICS by kamata 2020/11/16
      module topics_mod
      implicit none

! control
      integer :: ktoprnk = -1, lgtpc = 0
! ktoprnk : rank of TOPICS in MPI_COMM_WORLD
! lgtpc   : TOPICS and SONIC integrated calculation flags
!           = 0 : not use TOPICS, = 1 : integrated calculation with TOPICS

! time control
      real(8) :: dtduc = 0.0_8, dtim_s = 0.0_8, tduc = 0.0_8
! dtduc  : time to tduc [s] ( = tduc - time )
! dtim_s : backup of dtim set with SOLDOR [s]
! tduc   : time of data send / recieve cycle of TOPICS and SONIC [s]
!          ( data update cycle )

! mesh
      integer :: iybf
!        icmps                      icmpe
!         V                          V
!        |31|32|33|34|35|36|37|38|39|40|41|42|43|....   |78|79|80|
!
!        nroblk-1                   nrbc (iybc)
!         V                          V
!        |50|49|48|47|46|45|44|43|42|41|40|39|38|....   |03|02|01|
!                                   A  A
!                                  42  41 
!                                 (iybf)

! buffer for sending and receiving data
! modified 1/1 lines addition and reconsider of passed data by kamata 2022/02/23
!ik   integer, parameter :: ndma = 2, ndmk = 17, ndmr = 101
      integer, parameter :: ndmk = 14
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
!ik   integer     :: gmion = 0, gnion = 0, kcon = 0, lexdt = 0
      integer     :: gmion = 0, gnion = 0, gnro = 0, kcon = 0, lexdt = 0
     >             , lhist = 0, nro = 0, nroblk = 0
      real(8)     :: dtcal = 0.0_8, gtim = 0.0_8
! modified 1/1 lines addition and reconsider of passed data by kamata 2022/02/23
!ik  >             , gro(ndmr) = 0.0_8, gvr(ndmr,0:ndma,ndmk) = 0.0_8
! modified 1/2 lines density of SONIC impurities to TOPICS by kamata 2022/03/12
!ik   real(8), allocatable :: gro(:), gvr(:,:,:)
! modified 1/1 lines TOPICS calculation using arbitrary SONIC mesh data by kamata 2022/08/14
! returnd 1/1 lines changed base of sent/received data from ρ_vol to ρ_psi by kamata 2022/11/24
      real(8), allocatable :: gdn(:,:), gro(:), gvr(:,:,:)
!ik   real(8), allocatable :: gdn(:,:), gro(:), gvl(:), gvr(:,:,:)
! added 1 line bug for SONIC + IMPACT by kamata 2022/09/05
      real(8), allocatable :: groh(:)
      integer        nddnz(3)
      character   :: gcnm(ndmk)*5 =
     >    (/ 'den  ', 'tem  ', 'Pflxh', 'Hflxh'
     >     , 'Dh   ', 'Xh   ', 'Vinh ', 'Vhph '
! modified 3/2 lines addition and reconsider of passed data by kamata 2022/02/23
!ik  >     , 'Snth ', 'Wnth ', 'Sexh ', 'Wexh '
!ik  >     , 'Wrdh ', 'Wjlh ', 'Sexh ', 'Wexh '
!ik  >     , 'dnnbi' /)
     >     , 'Snth ', 'Wnth ', 'Wrdh ', 'Sexh '
     >     , 'Wexh ', 'dnnbi' /)

! dtcal  : time to exchange data between TOPICS and SONIC
! gcnm   : data name
! gdn    : impurity density for each valence near the core edge
! gmion  : number of ion type ( bulk ions + impurities )
! gnion  : number of bulk ion type
! gnro   : SONIC on mesh size for TOPICS ( = icaxs - icmps + 1 + 1 )
! gro    : volume rho on mesh [-]
! groh   : volume rho on half mesh [-]
! gtim   : time
! gvl    : volume on mesh [m^3]
! gvr    : buffer area for data to be sent and received
! kcon   : SONIC execution status flag
!          = 0 : normal, = 1 : abnormal, = -1 : time limit
! lexdt  : data transfer flag between TOPICS and SONIC
!          = 1 : ON, = 0 : OFF, = -1 : end of integrated calculation
! lhist  : flag to output the value of sent / received data to a file
!          = 1 : output, = 0 : no output
! nddnz  : array size
!          (1) first array size of gdn ( rho mesh )
!          (2) second array size of gdn ( alloacte valence )
!          (3) total number of valences
! ndma   : dimension size for ion type
! ndmk   : number of variables to send and receive
! ndmr   : number of rho mesh to send and receive
! nro    : number of rho mesh
! nroblk : number of rho mesh in the bulk area

!:: name
!      p:  partcile, h: heat
!      sn: neutral  sb: beam  sr: radiation   sj: joule
!      sc: source except neutral
      integer, parameter :: kden =  1, ktem =  2, kpfl =  3, khfl =  4
      integer, parameter :: kpdf =  5, khdf =  6, kpvp =  7, khvp =  8
! modified 3/2 lines addition and reconsider of passed data by kamata 2022/02/23
!ik   integer, parameter :: kpsn =  9, khsn = 10, kpsb = 11, khsb = 12
!ik   integer, parameter :: khsr = 13, khsj = 14, kpsc = 15, khsc = 16
!ik   integer, parameter :: kdnb = 17
      integer, parameter :: kpsn =  9, khsn = 10, khsr = 11, kpsc = 12
      integer, parameter :: khsc = 13, kdnb = 14

      end module topics_mod
