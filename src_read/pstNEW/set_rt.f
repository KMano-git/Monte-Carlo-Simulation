!***********************************************************************
      subroutine set_rt(cnam, nty, var, nv, cvnm)
!***********************************************************************
      use csize
      use cntcom
      use cntpls
      use cntsrc
      use cplimp
      use cunit
      use catcom
      implicit none

!::argument
      character(20), intent(in)  :: cnam
      integer, intent(in) :: nty
      real*8, intent(out) :: var(ndmc)
      integer, intent(out) :: nv
      character(30), intent(out) :: cvnm

      integer :: nrt, is0, is1, iz, ic, ia, ig
      real*8 :: sum1, sum2
      character(30) :: ctmp

      nrt=0
      if(cnam(3:3).eq."0")nrt=0 ! (Z0-ismax)/Ni
      if(cnam(3:3).eq."1")nrt=1 ! (Z1-ismax)/Ni
      if(cnam(3:3).eq."a")nrt=2 ! (Z0-ismax)/(Ni+N0+Ng)
      if(cnam(3:3).eq."n")nrt=3 ! Z0/(N0+Ng)

      is0 = 0
      if(nrt==1)is0=1 
      is1 = ismaxL(nty)
      if(nrt==3) is1 = 0

      nv = ncmax
      if(nrt>1)nv=ncmax2

      ia = 1
      ig = 1
      do ic = 1, nv
        sum1 = 0.0d0
        sum2 = 0.0d0
        do iz = is0, is1
          sum1 = sum1 + tdnzL(iz,ic,nty)
        enddo

        if(nrt==2)then
          sum2 = deni(ic,ia)+tden0(ic,ig)+tdeng(ic,ig)
        elseif(nrt==3)then
          sum2 = tden0(ic,ig)+tdeng(ic,ig)
        else
          sum2 = deni(ic,ia)
        endif

        if( sum2.gt.void_ni ) then
          var(ic) = sum1/sum2
        else
          var(ic) = 0.0d0
        endif
      enddo

      if(nrt==0)then
        write(ctmp,'(a2,"(0-",i2,")/Ni")')
!     >                 catmzL(nty),ismaxL(nty)
     >                 catmz(nty),ismaxL(nty)
      elseif(nrt==1) then
        write(ctmp,'(a2,"(1-",i2,")/Ni")')
!     >                 catmzL(nty),ismaxL(nty)
     >                 catmz(nty),ismaxL(nty)
      elseif(nrt==2) then
        write(ctmp,'(a2,"(0-",i2,")/N(i+0+g)")')
!     >                 catmzL(nty),ismaxL(nty)
     >                 catmz(nty),ismaxL(nty)
      elseif(nrt==3) then
        write(ctmp,'(a2,"(0)/N(0+g)")')
!     >                 catmzL(nty)
     >                 catmz(nty)
      endif

      cvnm = trim(cnam)//" "//trim(ctmp)

      write(n6,*)"ctmp: ",ctmp
      write(n6,*)"cnam: ",cnam
      write(n6,*)"cvnm: ",cvnm

      end subroutine

