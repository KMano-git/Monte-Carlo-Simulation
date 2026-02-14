      MODULE Phelps
      IMPLICIT NONE
!
      integer,parameter :: nn=41
      real(8),save ::El0(1:nn),Qm0(1:nn), log10El0(1:nn), log10Qm0(1:nn)
      real(8),save :: log10EcmMin =-2.d0, log10EcmMax =2.d0,
     &                log10ElabMin=10.d0, log10ElabMax=0.d0

! Database from Phelps ( Elab[eV] & Qm0[E-20 m^2] and their log10 )
c$$$  DATA El0      /  0.10000000E+00,  0.13340000E+00,  0.17760000E+00,  0.23710000E+00,
c$$$                   0.31620000E+00,  0.42170000E+00,  0.56230000E+00,  0.74990000E+00,
c$$$                   0.10000000E+01,  0.13340000E+01,  0.17760000E+01,  0.23710000E+01,
c$$$                   0.31620000E+01,  0.42170000E+01,  0.56230000E+01,  0.74990000E+01,
c$$$                   0.10000000E+02,  0.13340000E+02,  0.17760000E+02,  0.23710000E+02,
c$$$                   0.31620000E+02,  0.42170000E+02,  0.56230000E+02,  0.74990000E+02,
c$$$                   0.10000000E+03,  0.13340000E+03,  0.17760000E+03,  0.23710000E+03,
c$$$                   0.31620000E+03,  0.42170000E+03,  0.56230000E+03,  0.74990000E+03,
c$$$                   0.10000000E+04,  0.13340000E+04,  0.17760000E+04,  0.23710000E+04,
c$$$                   0.31620000E+04,  0.42170000E+04,  0.56230000E+04,  0.74990000E+04,
c$$$                   0.10000000E+05 /
c$$$
c$$$  DATA Qm0      /  0.19800000E+02,  0.18940000E+02,  0.18200000E+02,  0.17400000E+02,
c$$$                   0.16550000E+02,  0.15640000E+02,  0.14700000E+02,  0.13600000E+02,
c$$$                   0.12500000E+02,  0.11400000E+02,  0.10350000E+02,  0.93000000E+01,
c$$$                   0.82500000E+01,  0.71500000E+01,  0.61500000E+01,  0.51500000E+01,
c$$$                   0.42500000E+01,  0.34300000E+01,  0.26800000E+01,  0.21000000E+01,
c$$$                   0.15900000E+01,  0.11900000E+01,  0.87000000E+00,  0.63000000E+00,
c$$$                   0.44500000E+00,  0.31500000E+00,  0.21200000E+00,  0.13600000E+00,
c$$$                   0.88000000E-01,  0.55000000E-01,  0.34000000E-01,  0.20200000E-01,
c$$$                   0.11400000E-01,  0.64500000E-02,  0.34400000E-02,  0.19000000E-02,
c$$$                   0.10300000E-02,  0.56000000E-03,  0.30000000E-03,  0.16000000E-03,
c$$$                   0.85000000E-04 /

      DATA log10El0 /
     &-1.00000000E+00,-8.75000000E-01,-7.50000000E-01, -6.25000000E-01,
     &-5.00000000E-01,-3.75000000E-01,-2.50000000E-01, -1.25000000E-01,
     & 0.00000000E+00, 1.25000000E-01, 2.50000000E-01,  3.75000000E-01,
     & 5.00000000E-01, 6.25000000E-01, 7.50000000E-01,  8.75000000E-01,
     & 1.00000000E+00, 1.12500000E+00, 1.25000000E+00,  1.37500000E+00,
     & 1.50000000E+00, 1.62500000E+00, 1.75000000E+00,  1.87500000E+00,
     & 2.00000000E+00, 2.12500000E+00, 2.25000000E+00,  2.37500000E+00,
     & 2.50000000E+00, 2.62500000E+00, 2.75000000E+00,  2.87500000E+00,
     & 3.00000000E+00, 3.12500000E+00, 3.25000000E+00,  3.37500000E+00,
     & 3.50000000E+00, 3.62500000E+00, 3.75000000E+00,  3.87500000E+00,
     & 4.00000000E+00 /

      DATA log10Qm0 /
     & 0.12966652E+01, 0.12773800E+01, 0.12600714E+01,  0.12405492E+01,
     & 0.12187980E+01, 0.11942367E+01, 0.11673173E+01,  0.11335389E+01,
     & 0.10969100E+01, 0.10569049E+01, 0.10149403E+01,  0.96848295E+00,
     & 0.91645395E+00, 0.85430604E+00, 0.78887512E+00,  0.71180723E+00,
     & 0.62838893E+00, 0.53529412E+00, 0.42813479E+00,  0.32221929E+00,
     & 0.20139712E+00, 0.75546961E-01,-0.60480747E-01, -0.20065945E+00,
     &-0.35163999E+00,-0.50168945E+00,-0.67366414E+00, -0.86646109E+00,
     &-0.10555173E+01,-0.12596373E+01,-0.14685211E+01, -0.16946486E+01,
     &-0.19430951E+01,-0.21904403E+01,-0.24634416E+01, -0.27212464E+01,
     &-0.29871628E+01,-0.32518120E+01,-0.35228787E+01, -0.37958800E+01,
     &-0.40705811E+01 /
!
      public :: sigma_H2xH2_t
!
      CONTAINS
! ----------------------------------------
!      Subroutine sigma_H2xH2_m(lnE_cm,m_alfa,m_beta,sig)
      Subroutine sigma_H2xH2_t(lnE_cm,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE_cm !, m_alfa, m_beta
      REAL(8),INTENT(OUT) :: sig
      real(8) :: log10El, log10Qm ! E=E_Lab=(m_a+m_b)/m_b*E_cm
      real(8) :: X, Y, Y1, Y0, X1, X0, flag
      integer :: i,j,ii
      real(8) :: fac=2.d0
!
c$$$      log10El=log10( (m_alfa+m_beta)/m_beta*exp(lnE_cm) ) ! convert ln(E_CM) -> log10(E_lab)
c$$$      if(log10ElabMin==10.d0)then
c$$$         log10ElabMin=log10( (m_alfa+m_beta)/m_beta*(10.d0**log10EcmMin) ) ! convert E_CM -> E_lab
c$$$         log10ElabMax=log10( (m_alfa+m_beta)/m_beta*(10.d0**log10EcmMax) ) ! convert E_CM -> E_lab
c$$$         write(6,*)log10ElabMin, log10El0(1), log10El0(nn), log10ElabMax
c$$$      endif
      log10El=log10( fac*exp(lnE_cm) ) ! convert ln(E_CM) -> log10(E_lab)
      if(log10ElabMin==10.d0)then
         log10ElabMin=log10( fac*(10.d0**log10EcmMin) ) ! convert E_CM -> E_lab
         log10ElabMax=log10( fac*(10.d0**log10EcmMax) ) ! convert E_CM -> E_lab
!!         write(6,*)log10ElabMin, log10El0(1), log10El0(nn), log10ElabMax
      endif
!
!---  Linear Interpolation
      X=log10El
      Lp_domain: do ii=1,1
      IF(X<=log10ElabMin .AND. X <=log10El0(1)) THEN
         j=1
         X=MIN(log10ElabMin,log10El0(1))
         X0=log10El0(1); X1=log10El0(2)
         Y0=log10Qm0(1); Y1=log10Qm0(2)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         exit Lp_domain
      ELSEIF(X>log10ElabMin .AND. X<=log10El0(1) )then
         j=2
         X0=log10El0(1); X1=log10El0(2)
         Y0=log10Qm0(1); Y1=log10Qm0(2)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         exit Lp_domain
      ELSEIF(X>log10El0(1) .AND. X<log10El0(nn))then
         j=3
         do i=1,nn-1
            X0=log10El0(i); X1=log10El0(i+1)
            Y0=log10Qm0(i); Y1=log10Qm0(i+1)
            flag=(X0-X)*(X1-X)
!     if(flag==0.d0) Y=Y0
            if(flag <=0.d0) then
!     Y=(Y1-Y0)*(X1-X0)*(X-X0)+Y0
               Y=Y1*(X-X0)/(X1-X0)-Y0*(X-X1)/(X1-X0)
               exit Lp_domain
            endif
         enddo
      ELSEIF(X>=log10El0(nn) .AND. X<log10ElabMax)then
         j=4
         X0=log10El0(nn-1); X1=log10El0(nn)
         Y0=log10Qm0(nn-1); Y1=log10Qm0(nn)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         exit Lp_domain
      ELSEIF(X>=log10ElabMax .AND. X>=log10El0(nn))then
         j=5
         X=MAX(log10ElabMax,log10El0(nn))
         X0=log10El0(nn-1); X1=log10El0(nn)
         Y0=log10Qm0(nn-1); Y1=log10Qm0(nn)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         exit Lp_domain
      ELSE
         Y=0.d0
         exit Lp_domain
      ENDIF
      ENDDO Lp_domain
      log10Qm=Y
!
      sig=10.d0**(log10Qm)*1.d-20 ![m^2]

      Return
      END Subroutine sigma_H2xH2_t
!
!     ----------------------------------------------------------------
      END MODULE Phelps
