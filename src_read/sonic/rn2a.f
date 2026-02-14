! real number to ascii
      subroutine rn2a( rnum, ndec, lntxt, atxt )
      use cunit, only : n6
      implicit none

! arguments
      real(8),   intent(in)  :: rnum         ! real number
      integer,   intent(in)  :: ndec         ! number of decimals
      integer,   intent(in)  :: lntxt        ! length of atxt
      character, intent(out) :: atxt*(lntxt) ! ascii text

! local variables
      character  fmt*20
      integer    ip, jw

! check decimals
      if( lntxt < abs(ndec)+3 ) then
        write(n6,*) 'lntxt < ndec+3, lntxt = ', lntxt, ' ndec = ', ndec
     >    , ' at rn2a.'
        call wexit("rn2a","lntxt < abs(ndec)+3")
      endif

! real number to ascii
      jw   = max( ndec, 0 )
      write(fmt,'(a,i2.2,a,i1,a)') '(f',lntxt,'.', jw, ')'
      write(atxt,fmt) rnum

! 'XXX.' --> 'XXX'
      ip = len_trim( atxt )
      if( atxt(ip:ip) == '.' ) atxt(ip:ip) = ' '

! left justify
      atxt = adjustl( atxt )

      return
      end
