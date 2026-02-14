c########################################################################
c
c     Atomic Hydrogen dat by CR model 
c
c     Data from amujuel Apr 28 2018    
c     Equation and parameter for rate coefficient: 
c             sum(log(Te)**j * sum(B_ij * log(E)**i),i=0..8),j=0..8
c
c     Yamoto 2021/11/12
c########################################################################

c      program amjuel
c        call sigv_test
c        stop
c      end

!########################################################################
      subroutine asigv_test
!########################################################################

      implicit none
      
      integer k,m,n
      real*8 te,ne
      real*8 sigvv(1:2)
      real*8 sv_atomhcr

      call set_atomhcr

c      m=2
c      do i=0,8
c       write(*,'(i4,9e17.9)')i,(bij(j,i,m),j=0,8)
c      enddo
c      write(*,*)

      ne=1.0d20 

      do m=1,2
        write(10+m,'(9a18)')
     &     "#te","1e8","1e10","1e12","1e14","1e16","1e18"
        do k=0,40
          te = 10.0d0**(dble(k-10)*0.1d0)

          do n=1,2
            ne=10.0d0**dble(12+n*2)
            
c            write(*,'(2e18.9)')te,ne

             sigvv(n) = sv_atomhcr(te,ne,m)
c            write(*,*)te,sigv
           enddo
        write(10+m,'(8e18.9)')te,sigvv
        enddo
      enddo

      return
      end

!########################################################################
      real*8 function sv_atomhcr(te,ne,nreac)
!########################################################################
!     Equation and parameter for rate coefficient: 
!             sum(log(Te)**j * sum(B_ij * log(E)**i),i=0..8),j=0..8      
!########################################################################
      use AM_data, only : a_bij
      implicit none
      real*8 sigv
      real*8 te, ne, zte, zne
      integer n

      integer i,j,nreac
      real*8 logte,logne,cte,cne,fac
      
      zne=ne
      zne = max(1.0d14, zne)
      zne = min(1.0d22, zne)
      zte=te
      zte = max(0.1d0,  zte)
      zte = min(2.0d4,  zte)
c       TEMIN = 0.10000D 00 EV
c       TEMAX = 2.00000D 04 EV
c       NEMIN = 1.00000D 08 1/CM3 1D14/m3
c       NEMAX = 1.00000D 16 1/CM3 1D22/m3

      logne = dlog(ne * 1.0d-8 * 1.0d-6) ! ne 10^8 cm^-3
      logte = dlog(te)

      sigv = 0.0d0
      cte = 1.0d0
      do j=0,8
        cne = 1.0d0
        fac = 0.0d0
        do i=0,8
          fac = fac+a_bij(i,j,nreac)*cne
          cne = cne*logne
        enddo
        sigv = sigv+cte*fac
        cte = cte*logte
      enddo
      sigv=dexp(sigv)*1.0d-6
      sigv=max(sigv,1.0d-30)

!      if(nreac==1) write(700,*) zte, sigv
!      if(nreac==2) write(701,*) zte, sigv

      sv_atomhcr = sigv
      
      end function

!########################################################################
      subroutine set_atomhcr
!########################################################################
      use cunit, only : n6
      use AM_data, only : a_bij
      implicit none
      
      integer m
c
      write(n6,'(2x, "  atomhcr is used")')
      
ccc EI: H + e -> H+ + 2e
ccc AMJUEL H4 2.1.5
ccc    Effective hydrogenic ionization rate    Data: K. Sawada/T. Fujimoto (redone: 2016,
ccc   extend Te range of fit validity from 0.1 -- 1e3 to 0.1 -- 2e4 eV)

      m=1
      
      a_bij(0:8,0,m)  = (/
     &    -3.248025D+01, -5.440669D-02,  9.048888D-02, 
     &    -4.054079D-02,  8.976514D-03, -1.060334D-03,  
     &     6.846238D-05, -2.242955D-06,  2.890438D-08 /)
c
      a_bij(0:8,1,m)  = (/
     &     1.425332D+01, -3.594347D-02, -2.014729D-02,
     &     1.039774D-02, -1.771792D-03,  1.237467D-04,
     &    -3.130184D-06, -3.051995D-08,  1.888148D-09 /)
c
      a_bij(0:8,2,m)  = (/
     &    -6.632235D+00,  9.255558D-02, -5.580210D-03,
     &    -5.902219D-03,  1.295610D-03, -1.056722D-04,
     &     4.646310D-06, -1.479612D-07,  2.852251D-09 /)
c
      a_bij(0:8,3,m)  = (/
     &     2.059544D+00, -7.562462D-02,  1.519596D-02, 
     &     5.803498D-04, -3.527285D-04,  3.201534D-05,  
     &    -1.835197D-06,  9.474014D-08, -2.342506D-09 /)
c
      a_bij(0:8,4,m)  = (/
     &    -4.425370D-01,  2.882634D-02, -7.285771D-03, 
     &     4.643390D-04,  1.145701D-06,  8.493663D-07,  
     &    -1.001033D-08, -1.476839D-08,  6.047700D-10 /)
c
      a_bij(0:8,5,m)  = (/
     &     6.309382D-02, -5.788687D-03,  1.507383D-03,  
     &    -1.201551D-04,  6.574488D-06, -9.678783D-07, 
     &     5.176266D-08,  1.291552D-09, -9.685157D-11 /)
c
      a_bij(0:8,6,m)  = (/
     &    -5.620092D-03,  6.329106D-04, -1.527778D-04,  
     &     8.270125D-06,  3.224102D-08,  4.377403D-08,  
     &    -2.622922D-09, -2.259663D-10,  1.161439D-11 /)
c
      a_bij(0:8,7,m)  = (/
     &     2.812017D-04, -3.564133D-05,  7.222727D-06, 
     &     1.433019D-07, -1.097431D-07,  7.789032D-09,  
     &    -4.197729D-10,  3.032260D-11, -8.911077D-13 /)
c
      a_bij(0:8,8,m)  = (/
     &    -6.011143D-06,  8.089651D-07, -1.186213D-07,  
     &    -2.381081D-08,  6.271174D-09, -5.483010D-10, 
     &     3.064612D-11, -1.355903D-12,  2.935080D-14 /)
     
c       TEMIN = 0.10000D 00 EV
c       TEMAX = 2.00000D 04 EV
c       NEMIN = 1.00000D 08 1/CM3 1D14
c       NEMAX = 1.00000D 16 1/CM3 1D22


ccc EI: H + e -> H+ + 2e
ccc AMJUEL H10 2.1.5
ccc Electron energy loss weighted rate coefficient. Data: Sawada/Fujimoto, \cite{kn:Sawada}
ccc (redone May 2016: extend Te range of fit validity from 0.1 -- 1e3 now to 0.1 -- 2e4 eV)
      m=2
  
      a_bij(0:8,0,m)  = (/
     &    -2.497580D+01,  1.081654D-03, -7.358936D-04, 
     &     4.122399D-04, -1.408153D-04,  2.469731D-05,  
     &    -2.212824D-06,  9.648140D-08, -1.611904D-09 /)
c
      a_bij(0:8,1,m)  = (/
     &     1.004449D+01, -3.189475D-03,  2.510128D-03, 
     &    -7.707041D-04,  1.031310D-04, -3.716939D-06,  
     &    -4.249705D-07,  4.164961D-08, -9.893424D-10 /)
c
      a_bij(0:8,2,m)  = (/
     &    -4.867953D+00, -5.852268D-03,  2.867459D-03,  
     &    -8.328668D-04,  2.056134D-04, -3.301571D-05, 
     &     2.831740D-06, -1.164969D-07,  1.785440D-09 /)
c
      a_bij(0:8,3,m)  = (/
     &     1.689422D+00,  7.744372D-03, -3.087364D-03, 
     &     4.707676D-04, -5.508612D-05,  7.305868D-06,  
     &    -6.000116D-07,  2.045212D-08, -1.790313D-10 /)
c
      a_bij(0:8,4,m)  = (/
     &    -4.103532D-01, -3.622291D-03,  1.327415D-03, 
     &    -1.424079D-04,  3.307340D-06,  5.256680D-09,  
     &     7.597020D-10,  1.799505D-09, -9.280890D-11 /)
c
      a_bij(0:8,5,m)  = (/
     &     6.469718D-02,  8.268568D-04, -2.830940D-04,  
     &     2.411848D-05,  5.707985D-07, -1.016946D-07, 
     &     3.517155D-09, -4.453196D-10,  2.002478D-11 /)
c
      a_bij(0:8,6,m)  = (/
     &    -6.215861D-03, -9.836596D-05,  3.017297D-05,  
     &    -1.474254D-06, -2.397869D-07,  1.518743D-08,  
     &     4.149085D-10, -6.803200D-12, -1.151856D-12 /)
c
      a_bij(0:8,7,m)  = (/
     &     3.289810D-04,  5.845698D-06, -1.479324D-06, 
     &    -4.633029D-08,  3.337390D-08, -1.770252D-09,  
     &    -5.289806D-11,  3.864395D-12, -8.694979D-15 /)
c
      a_bij(0:8,8,m)  = (/
     &    -7.335808D-06, -1.367574D-07,  2.423236D-08,  
     &     5.733871D-09, -1.512778D-09,  8.733801D-11, 
     &     7.196799D-13, -1.441034D-13,  1.734769D-15 /)

      return
      end subroutine

