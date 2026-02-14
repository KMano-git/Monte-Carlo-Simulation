      MODULE KrsticSchultz3
      IMPLICIT NONE
      real(8) :: extbase=1.d-2
      real(8),parameter,private :: uni=2.80028d-21 ! [m^2]
!
      public :: sigma_DxD_t,   sigma_DxD_m,
     &          sigma_DxD2_t,  sigma_DxD2_m,
     &          sigma_DixD_t,  sigma_DixD_m,
     &          sigma_pxH_t_k, sigma_pxH_m_k
      real(8),parameter :: Emin=1.d-1, Emax=1.d2 ! [eV]
      real(8) :: lnEmin0=log(Emin), lnEmax1=log(Emax), lnEmax0, lnEmin1
!
      CONTAINS
! ----------------------------------------
      Subroutine sigma_DxD_t(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
! D+D ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.187295d+03
      a(1) =  0.720980d+02
!
      b(1) =  0.496047d+00
      b(2) =  0.610483d-01
      b(3) =  0.210781d-02
      b(4) = -0.108431d-02
      b(5) =  0.492770d-03
! ----------------------
      sig=0.d0

      IF(lnE < lnEmin0) then
! constant
!!$       X0=lnEmin0
!!$       atmp=0.d0
!!$       btmp=0.d0
!!$       do j=0,jj
!!$          atmp=atmp+a(j)*X0**j
!!$       enddo
!!$       do j=1,jj
!!$          btmp=btmp+b(j)*X0**j
!!$       enddo
!!$       sig=atmp/(1.d0+btmp)*uni
! linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
! linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
!
      Return
      END Subroutine sigma_DxD_t
!
! ----------------------------------------
      Subroutine sigma_DxD_m(lnE,sig) ! Identical with sigma_DxD_t
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
! D+D ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.187295d+03
      a(1) =  0.720980d+02
!
      b(1) =  0.496047d+00
      b(2) =  0.610483d-01
      b(3) =  0.210781d-02
      b(4) = -0.108431d-02
      b(5) =  0.492770d-03
! ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
! constant
!!$       X0=lnEmin0
!!$       atmp=0.d0
!!$       btmp=0.d0
!!$       do j=0,jj
!!$          atmp=atmp+a(j)*X0**j
!!$       enddo
!!$       do j=1,jj
!!$          btmp=btmp+b(j)*X0**j
!!$       enddo
!!$       sig=atmp/(1.d0+btmp)*uni
! linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
! linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
      Return
      END Subroutine sigma_DxD_m
!
! ----------------------------------------
      Subroutine sigma_DxD2_t(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     D+D2  ----------------
      a=0.d0;  b=0.d0
      a(0) =  0.188725d+03
      a(1) = -0.139781d+02
      a(2) = -0.811171d+00
!     ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,2
            atmp=atmp+a(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
!----------------------------------------------------------------
      atmp=0.d0
      btmp=0.d0
      do j=0,2
         atmp=atmp+a(j)*lnE**j
      enddo
      sig=atmp/(1.d0+btmp)*uni
!----------------------------------------------------------------
      Return
      END Subroutine sigma_DxD2_t
!
!     ----------------------------------------
      Subroutine sigma_DxD2_m(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     D+D2 Momentum Transfer
      a=0.d0;  b=0.d0
      a(0) =  0.197634d+02
      a(1) = -0.186056d+02
      a(2) =  0.899988d+01
      a(3) = -0.156354d+01
      a(4) = -0.195343d+00
      a(5) =  0.915858d-01
      a(6) = -0.766847d-02
!
      b(1) =  0.224238d+00
      b(2) =  0.295438d+00
      b(3) =  0.461722d-01
!     ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,6
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,3
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
      Return
      END Subroutine sigma_DxD2_m
!
!     ----------------------------------------
      Subroutine sigma_pxH_t_k(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     p+H ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.591039d+03
      a(1) = -0.877354d+02
      a(2) =  0.256830d+01
!     ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
      Return
      END Subroutine sigma_pxH_t_k
!
!     ----------------------------------------
      Subroutine sigma_pxH_m_k(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     p+H ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.324832d+03
      a(1) = -0.392017d+02
      a(2) =  0.124924d+01
!     ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
      Return
      END Subroutine sigma_pxH_m_k
!
!     ----------------------------------------
      Subroutine sigma_DixD_t(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     Di+D ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.637254d+03
      a(1) = -0.125310d+03
      a(2) =  0.167078d+02
      a(3) = -0.990209d+00
!
!     ----------------------
      sig=0.d0

      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
!
      Return
      END Subroutine sigma_DixD_t
!
!     ----------------------------------------
      Subroutine sigma_DixD_m(lnE,sig)
      IMPLICIT NONE
      REAL(8),INTENT( IN) :: lnE
      REAL(8),INTENT(OUT) :: sig
      integer,parameter :: jj=6
      real(8) :: a(0:jj),b(0:jj)
      real(8) :: atmp, btmp
      integer :: j
      real(8) :: X, Y, X0, X1, Y0, Y1
!
!     Di+D ------------------
      a=0.d0;  b=0.d0
      a(0) =  0.347935d+03
      a(1) =  0.419075d+03
      a(2) =  0.117426d+03
      a(3) = -0.324821d+01
      a(4) = -0.105819d+01
!
      b(1) =  0.131735d+01
      b(2) =  0.482389d+00
      b(3) =  0.401043d-01
!     ----------------------
      sig=0.d0
      IF(lnE < lnEmin0) then
!     constant
!     !$       X0=lnEmin0
!     !$       atmp=0.d0
!     !$       btmp=0.d0
!     !$       do j=0,jj
!     !$          atmp=atmp+a(j)*X0**j
!     !$       enddo
!     !$       do j=1,jj
!     !$          btmp=btmp+b(j)*X0**j
!     !$       enddo
!     !$       sig=atmp/(1.d0+btmp)*uni
!     linear extrapolation in log-space
         lnEmin1=(1.d0+extbase)*lnEmin0
         X0=lnEmin0; X1=lnEmin1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X0)+Y0
         sig=exp(Y)
!
      ELSEIF(lnE > lnEmax1) then
!     linear extrapolation in log-space
         lnEmax0=(1.d0-extbase)*lnEmax1
         X0=lnEmax0; X1=lnEmax1; X=lnE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X0**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X0**j
         enddo
         Y0=log(atmp/(1.d0+btmp)*uni)
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*X1**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*X1**j
         enddo
         Y1=log(atmp/(1.d0+btmp)*uni)
         Y=(Y1-Y0)/(X1-X0)*(X-X1)+Y1
         sig=exp(Y)
      ELSE
         atmp=0.d0
         btmp=0.d0
         do j=0,jj
            atmp=atmp+a(j)*lnE**j
         enddo
         do j=1,jj
            btmp=btmp+b(j)*lnE**j
         enddo
         sig=atmp/(1.d0+btmp)*uni
      ENDIF
      Return
      END Subroutine sigma_DixD_m
!     ----------------------------------------------------------------
      END MODULE KrsticSchultz3
