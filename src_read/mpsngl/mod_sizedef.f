! modified for TOPICS by kamata 2019/08/05
! original : sonicV4/MPMD5/WFL/mod_sizedef.f
      module mod_sizedef

      integer, parameter :: ndnpe = 40000 ! PE number
      integer, parameter :: ndgrp =  6    ! prg. number

      integer, parameter :: lnnam =  10   ! length of prg-name or data name
      integer, parameter :: lnexd = 200   ! length of spec(exchange data)

      integer, parameter :: ldfil =  20
      integer, parameter :: ldlin = 120

!ik modified 1/1 lines combining TOPICS and IMPACT with MPMD model by kamata 2019/08/05
!ik   integer, parameter :: ndvio = 200   ! displace, typelist
      integer, parameter :: ndvio = 300   ! displace, typelist

! modified 1/1 lines replace all include files with module files by kamata 2021/08/18
!ik   integer, parameter :: ndtyp = 20    ! TOK, DIV, NTL, ...
      integer, parameter :: ndtyp = 50    ! TOK, DIV, NTL, ...
      integer, parameter :: ndprc = 100   ! group * (step + ~2 )
      integer, parameter :: ndexd = 50    ! TOK 1 => 2 exec

      end module  mod_sizedef
