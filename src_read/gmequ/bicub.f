C
C***********************************************************************
        SUBROUTINE BICUBE(C,X0,Y0,S,SX,SY,SXY,SXX,SYY)
C***********************************************************************
! modified 2/4 lines organize local variables and include files by kamata 2021/06/16
!ik     implicit real*8 (a-h,o-z)
!ik     REAL*8 C(4,4),X0,Y0,S,SX,SY,SXY,SXX,SYY
      use clocal, only : xim1, yjm1
      implicit none
! arguments
      real(8), intent(in)  :: c(4,4), x0, y0
      real(8), intent(out) :: s, sx, sxx, sxy, sy, syy
C************
C
C   SUBROUTINE BICUBE
C
C
C       THIS SUBROUTINE EVALUATES A BICUBIC POLYNOMIAL S, AND ITS
C       DERIVATIVES THROUGH ORDER 2. SPECIFICALLY, GIVEN VALUES
C       X=X0 AND Y=Y0 IN THE RECTANGLE:
C               R(I,J): X(I-1).LE.X.LE.X(I)
C                       Y(J-1).LE.Y.LE.Y(J)
C       AND GIVEN THE 16 COEFFICIENTS OF THE BICUBIC INTERPOLATING
C       POLYNOMIAL ASSOCIATED WITH R(I,J), THIS SUBROUTINE COMPUTES
C       THE VALUE:
C                   S = S(X0,Y0)
C       AND THE DERIVATIVES:
C                  SX = DS/DX
C                  SY = DS/DY
C                 SXY = DS/DXDY
C                 SXX = DDS/DXDX
C                 SYY = DDS/DYDY
C       EVALUATED AT THE POINT X0,Y0.
C
C   INPUT
C
C       C   IS A 4 BY 4 ARRAY WHICH CONTAINS THE COEFFICIENTS OF
C           THE BICUBIC POLYNOMIAL S. C IS NOT ALTERED BY THE
C           SUBROUTINE.
C       X0  IS THE VALUE OF THE X-ARGUMENT AT WHICH S IS TO BE
C           EVALUATED. X0 IS NOT ALTERED BY THE SUBROUTINE.
C       Y0  IS THE VALUE OF THE Y-ARGUMENT AT WHICH S IS TO BE
C           EVALUATED. Y0 IS NOT ALTERED BY THE SUBROUTINE.
C  OUTPUT
C
C       S,SX,SY,SXY,SXX,SYY   AS DEFINED ABOVE.

C
C
C   P.W.GAFFNEY   6TH. MARCH 1978
C
C************
C
C       THE FOLLOWING ARRAYS ARE USED FOR WORKSPACE
! modified 1/5 lines organize local variables and include files by kamata 2021/06/16
!ik     REAL*8 X(1,4),Y(4,1),DX(1,4),DY(4,1),D2X(1,4),D2Y(4,1)
! local variables
      real(8) :: cij, d2x(1,4), d2y(4,1), dx(1,4), dy(4,1)
     >         , one = 1.0_8, two = 2.0_8, x(1,4), x1i, xbar
     >         , y(4,1), ybar, yj1, zero = 0.0_8
      integer    i, j

C
C       IN THE FOLLOWING COMMON BLOCK THE QUANTITIES
C       XIM1=X(I-1) AND YJM1=Y(J-1) ARE SET IN SUBROUTINE COEFF.
! added 1 line organize local variables and include files by kamata 2021/06/16
! deleted 2 lines replace all include files with module files by kamata 2021/08/18
!ik     real(8)    xim1, yjm1
!ik     COMMON/COBI/XIM1,YJM1
C
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik     DATA ZERO,ONE,TWO/0.D0,1.D0,2.D0/
C
        XBAR=X0-XIM1
        YBAR=Y0-YJM1
C
        X(1,1)=1.0D0
        X(1,2)=XBAR
        X(1,3)=XBAR*XBAR
        X(1,4)=XBAR**3
C
        Y(1,1)=1.0D0
        Y(2,1)=YBAR
        Y(3,1)=YBAR*YBAR
        Y(4,1)=YBAR**3
C

        DX(1,1)=ZERO
        DX(1,2)=ONE
        DX(1,3)=2.D0*XBAR
        DX(1,4)=3.D0*XBAR*XBAR
        DY(1,1)=ZERO
        DY(2,1)=ONE
        DY(3,1)=2.D0*YBAR
        DY(4,1)=3.D0*YBAR*YBAR
C
        D2X(1,1)=ZERO
        D2X(1,2)=ZERO
        D2X(1,3)=TWO
        D2X(1,4)=6.D0*XBAR
        D2Y(1,1)=ZERO
        D2Y(2,1)=ZERO
        D2Y(3,1)=TWO
        D2Y(4,1)=6.D0*YBAR
C
        S=ZERO
        SX=ZERO
        SY=ZERO
        SXY=ZERO
        SXX=ZERO
        SYY=ZERO
        DO 2 J=1,4
        YJ1=Y(J,1)
        DO 1 I=1,4
        X1I=X(1,I)
        CIJ=C(I,J)
        S=S+X1I*CIJ*YJ1
        SX=SX+DX(1,I)*CIJ*YJ1
        SY=SY+X1I*CIJ*DY(J,1)
        SXY=SXY+DX(1,I)*CIJ*DY(J,1)
        SXX=SXX+D2X(1,I)*CIJ*YJ1
    1   SYY=SYY+X1I*CIJ*D2Y(J,1)
    2   CONTINUE
        RETURN
        END
C
C***********************************************************************
        SUBROUTINE COEFF(N,M,X,Y,LDU,U,UX,UY,UXY,I,J,C)
C***********************************************************************
! modified 4/9 lines organize local variables and include files by kamata 2021/06/16
!ik     implicit real*8 (a-h,o-z)
!ik     REAL*8 X(N),Y(M),U(LDU,M),UX(LDU,M),UY(LDU,M),UXY(LDU,M),
!ik  *           C(4,4)
!ik     INTEGER N,M,I,J
      use clocal, only : xim1, yjm1
      implicit none
! arguments
      integer, intent(in)  :: i, j, ldu, m, n
      real(8), intent(in)  :: u(ldu,m), ux(ldu,m), uxy(ldu,m), uy(ldu,m)
     >    , x(n), y(m)
      real(8), intent(out) :: c(4,4)
! local variables
      real(8)    hx, hy, wk(8)
      integer    i1, ii, irow, j1, jcol, jj
C************
C
C   SUBROUTINE COEFF
C
C
C       THIS SUBROUTINE CALCULATES THE 16 COEFFICIENTS OF
C       THE INTERPOLATING BICUBIC POLYNOMIAL FOR THE
C       RECTANGLE R(I,J):
C                          X(I-1).LE.X.LE.X(I)
C                          Y(J-1).LE.Y.LE.Y(J)
C
C   INPUT
C
C       N   IS THE NUMBER OF DATA POINTS IN THE X-DIRECTION. N IS
C           NOT ALTERED BY THE SUBROUTINE.
C       M   IS THE NUMBER OF DATA POINTS IN THE Y-DIRECTION. M IS
C           NOT ALTERED BY THE SUBROUTINE.
C       X   IS A LINEAR ARRAY WHICH CONTAINS THE N GRID POINTS
C           X(1),...,X(N) ARRANGED IN ASCENDING ORDER. X IS NOT
C           ALTERED BY THE SUBROUTINE.
C       Y   IS A LINEAR ARRAY WHICH CONTAINS THE M GRID POINTS
C           Y(1),...,Y(M) ARRANGED IN ASCENDING ORDER. Y IS NOT
C           ALTERED BY THE SUBROUTINE.
C       LDU IS THE LEADING DIMENSION OF THE ARRAY U. LDU IS NOT
C           ALTERED BY THE SUBROUTINE.
C       U   IS A TWO DIMENSIONAL ARRAY WHICH CONTAINS THE GIVEN
C           FUNCTION VALUES U(X(IR),Y(JR)), IR=1 TO N, JR=1 TO M.
C           U IS NOT ALTERED BY THE SUBROUTINE.
C       UX  IS A TWO DIMENSIONAL ARRAY WHICH CONTAINS THE
C           DERIVATIVE OF U WITH RESPECT TO X EVALUATED AT
C           THE GRID POINTS X(IR),Y(JR). UX IS NOT ALTERED BY
C           THE SUBROUTINE.
C       UY  IS A TWO DIMENSIONAL ARRAY WHICH CONTAINS THE
C           DERIVATIVE OF U WITH RESPECT TO Y EVALUATED AT
C           THE GRID POINTS X(IR),Y(JR). UY IS NOT ALTERED BY
C           THE SUBROUTINE.
C       UXY IS A TWO DIMENSIONAL ARRAY WHICH CONTAINS THE
C           DERIVATIVE OF U WITH RESPECT TO X AND Y EVALUATED
C           AT THE GRID POINTS X(IR),Y(JR).  UXY IS NOT ALTERED BY
C           THE SUBROUTINE.
C       I,J ARE THE INTEGERS WHICH SPECIFY THE RECTANGLE R(I,J).
C           I AND J ARE NOT ALTERED BY THE SUBROUTINE.
C
C   OUTPUT
C
C       C   IS A 4 BY 4 ARRAY WHICH CONTAINS THE COEFFICIENTS OF THE
C           BICUBIC POLYNOMIAL:
C                           S(XBAR,YBAR) = XBAR*C*YBAR
C
C           WHERE XBAR IS A ROW VECTOR WITH COMPONENTS:
C
C                           XBAR(L) = (X-X(I-1))**(L-1) L=1,2,3,4
C
C           AND YBAR IS A COLUMN VECTOR WITH COMPONENTS:
C
C                           YBAR(L) = (Y-Y(J-1))**(L-1) L=1,2,3,4.
C
C
C   P.W.GAFFNEY  6TH. MARCH 1978
C
C************
C
C
! deleted 1 line organize local variables and include files by kamata 2021/06/16
!ik     REAL*8 WK(8)
C
C       IN THE FOLLOWING COMMON BLOCK THE VARIABLES XIM1 AND YJM1 ARE
C       DEFINED IN THIS SUBROUTINE AND USED IN SUBROUTINE BICUBE.
! added 1 line organize local variables and include files by kamata 2021/06/16
! deleted 2 lines replace all include files with module files by kamata 2021/08/18
!ik     real(8)    xim1, yjm1
!ik     COMMON /COBI/XIM1,YJM1
        XIM1=X(I-1)
        YJM1=Y(J-1)
C       WRITE(6,611) I,J,XIM1,YJM1
C       WRITE(6,612) X(I),Y(J)
 611    FORMAT(1H //,'  I,J,XIM1,YJM1',2I5,1P2E14.4)
 612    FORMAT( 1P2E14.4)
C
        DO 2 IROW=1,3,2
                I1=I-1
        IF(IROW.EQ.3)I1=I
        DO 1 JCOL=1,3,2
                J1=J-1
        IF(JCOL.EQ.3)J1=J
                C(IROW,JCOL)=U(I1,J1)
                C(IROW,JCOL+1)=UY(I1,J1)
                C(IROW+1,JCOL)=UX(I1,J1)
                C(IROW+1,JCOL+1)=UXY(I1,J1)
    1   CONTINUE
    2   CONTINUE
C
        HX=X(I)-XIM1
        DO 3 JJ=1,4
        WK(JJ)=(-3.D0*C(1,JJ)-2.D0*HX*C(2,JJ)+3.D0*C(3,JJ)
     &  -HX*C(4,JJ))/(HX*HX)
        WK(JJ+4)=(2.D0*C(1,JJ)+HX*C(2,JJ)-2.D0*C(3,JJ)
     &  +HX*C(4,JJ))/(HX**3)
    3   CONTINUE
C
        DO 4 JJ=1,4
        C(3,JJ)=WK(JJ)
        C(4,JJ)=WK(JJ+4)
    4   CONTINUE
C
        HY=Y(J)-YJM1
        DO 5 II=1,4
        WK(II)=(-3.D0*C(II,1)-2.D0*HY*C(II,2)+3.D0*C(II,3)
     &  -HY*C(II,4))/(HY*HY)
        WK(II+4)=(2.D0*C(II,1)+HY*C(II,2)-2.D0*C(II,3)
     &  +HY*C(II,4))/(HY**3)
    5   CONTINUE
C
        DO 6 II=1,4
        C(II,3)=WK(II)
        C(II,4)=WK(II+4)
    6   CONTINUE
        RETURN
        END
C
C***********************************************************************
        SUBROUTINE DERIVS(N,M,X,Y,LDU,U,UX,UY,UXY)
C***********************************************************************
!
!     working are wk(iwk) is deleted    2010/05/19  K. Shimizu
!
!-----------------------------------------------------------------------
! modified 4/8 lines organize local variables and include files by kamata 2021/06/16
!ik     implicit real*8 (a-h,o-z)
!ik     REAL*8 X(N),Y(M),U(LDU,M),
!ik  *        UX(LDU,M),UY(LDU,M),UXY(LDU,M)
!ik     INTEGER N,M
      implicit none
! arguments
      integer, intent(in)    :: ldu, m, n
      real(8), intent(in)    :: u(ldu,m), x(n), y(m)
! output variables
      real(8), intent(out) :: ux(ldu,m), uxy(ldu,m), uy(ldu,m)
! local variables
      integer    n1
C************
C
C   SUBROUTINE DERIVS
C
C       THE PURPOSE OF THIS SUBROUTINE IS TO COMPUTE APPROXIMATIONS
C       TO THE DERIVATIVES OF A FUNCTION U(X,Y) GIVEN VALUES OF THE
C       FUNCTION AT THE GRID POINTS OF A RECTANGULAR MESH R:
C
C                         X(1).LE.X.LE.X(N)
C                         Y(1).LE.Y.LE.Y(M).
C
C       THE APPROXIMATIONS ARE COMPUTED USING CUBIC SPLINE INTERPOLATION
C       IN THE X AND Y DIRECTIONS.
C
C       THE SUBROUTINE MAY BE USED FOR COMPUTING A BICUBIC SPLINE
C       FUNCTION WHICH INTERPOLATES THE GIVEN DATA. IN THIS
C       CASE IT SHOULD ONLY BE CALLED ONCE FOR EACH MESH R.
C
C
C       THIS IS A DRIVER SUBROUTINE WHICH CALLS DERIVD IN ORDER
C       TO COMPUTE THE NECESSARY DERIVATIVES AT THE MESH POINTS,
C       OF THE RECTANGULAR GRID.
C
C   INPUT
C
C       N   IS THE NUMBER OF DATA POINTS IN THE X-DIRECTION. N IS
C           NOT ALTERED BY THE SUBROUTINE.
C       M   IS THE NUMBER OF DATA POINTS IN THE Y-DIRECTION. M IS
C           NOT ALTERED BY THE SUBROUTINE.
C       X   IS A LINEAR ARRAY WHICH CONTAINS THE N GRID POINTS
C           X(1),...,X(N) ARRANGED IN ASCENDING ORDER. X IS NOT
C           ALTERED BY THE SUBROUTINE.
C       Y   IS A LINEAR ARRAY WHICH CONTAINS THE M GRID POINTS
C           Y(1),...,Y(M) ARRANGED IN ASCENDING ORDER.  Y IS NOT
C           ALTERED BY THE SUBROUTINE.
C       LDU IS THE LEADING DIMENSION OF THE ARRAY U. LDU IS NOT
C           ALTERED BY THE SUBROUTINE.
C       U   IS A TWO DIMENSIONAL ARRAY WHICH CONTAINS THE GIVEN
C           FUNCTION VALUES  U(X(I),Y(J)) I=1 TO N, J=1 TO M.
C           U IS NOT ALTERED BY THE SUBROUTINE.
C       IWK IS THE LENGTH OF THE WORK ARRAY WK. THE VALUE OF IWK
C           SHOULD BE AT LEAST 6 * MAX(N,M). IWK IS NOT ALTERED
C           BY THE SUBROUTINE.
C       WK  IS A WORK ARRAY OF LENGTH IWK.
C
C   OUTPUT
C
C       UX  CONTAINS DU/DX AT THE MESH POINTS
C       UY  CONTAINS DU/DY AT THE MESH POINTS
C       UXY CONTAINS DU/DXDY AT THE MESH POINTS
C
C
C   P.W.GAFFNEY    6TH.MARCH 1978
C
!::local variables
      integer  ndm
      parameter (ndm=513)
      real*8, dimension(ndm) :: h, th, t, alpha, betas, beta
!
C************
C
C       PARTITION WORKSPACE ARRAY
C       N1=MAX(N,M)
C       WK(I)=H(I)
C       WK(N1+I)=TH(I)
C       WK(2*N1+I)=T(I)
C       WK(3*N1+I)=ALPHA(I)
C       WK(4*N1+I)=BETAS(I)
C       WK(5*N1+I)=BETA(I)
C
!x        N1=MAX0(N,M)
!x        N2=2*N1
!x        N3=3*N1
!x        N4=4*N1
!x        N5=5*N1
!
      n1 = max0(n,m)
!KH   if( n1.ne.ndm ) then  ! KH111108
      if( n1.gt.ndm ) then
        write(6,'(2x,"stop at sub. derivs in bicub.f")')
        write(6,'(2x,"please change parameter size ndm =",i6)') ndm
        write(6,'(2x," n, m =",2i6,"  ndm => max(n,m)")') n, m
        call wexit("derivs","dimension size")
      endif
!
!KH   CALL DERIVD(N,M,X,Y,LDU,U,N1,h,th,t,alpha,beta,betas,
!KH  $              UX,UY,UXY)
!
      CALL DERIVD(N,M,X,Y,LDU,U,N1,h(1:n1),th(1:n1),t(1:n1),
     $   alpha(1:n1),beta(1:n1),betas(1:n1),UX,UY,UXY)
!
        RETURN
        END
C
C***********************************************************************
        SUBROUTINE DERIVD(N,M,X,Y,LDU,U,N1,H,TH,T,ALPHA,BETA,BETAS,
     $                    P,Q,S)
C***********************************************************************
! modified 4/11 lines organize local variables and include files by kamata 2021/06/16
!ik     implicit real*8 (a-h,o-z)
!ik     REAL*8 X(N),Y(M),U(LDU,M),H(N1),TH(N1),T(N1),ALPHA(N1),
!ik  $       BETA(N1),BETAS(N1),P(LDU,M),Q(LDU,M),S(LDU,M)
!ik     DATA E11,E18,E9,E2,E6/11.D0,18.D0,9.D0,2.D0,6.D0/
      implicit none
! arguments
      integer, intent(in)  :: ldu, m, n, n1
      real(8), intent(in)  :: u(ldu,m), x(n), y(m)
      real(8), intent(out) :: alpha(n1), beta(n1), betas(n1), h(n1)
     >    , p(ldu,m), q(ldu,m), s(ldu,m), t(n1), th(n1)
! local variables
      real(8) :: e11 = 11.0_8, e18 = 18.0_8, e2 = 2.0_8
     >         , e6 = 6.0_8, e6hx1, e6hxn, e6hy1, e6hym
     >         , e9 = 9.0_8
      integer    i, ib, im1, ip1, j, jb, mm1, mm3, nm1, nm3
C
C*******COMPUTE BOUNDARY CONDITIONS, BY APPROXIMATING DERIVATIVES
C       WITH 3RD-ORDER FORWARD DIFFERENCE FORMULA. IT IS ASSUMED
C       THAT THE FIRST THREE MESH POINTS AWAY FROM THE BOUNDARIES
C       ARE EQUALLY SPACED.
C
C
        E6HX1=E6*(X(2)-X(1))
        E6HXN=E6*(X(N)-X(N-1))
C
C       COMPUTE THE APPROXIMATION TO DU(X,Y)/DX  AT
C       X=X(1), Y=Y(J),J=1,...,M AND
C       X=X(N), Y=Y(J),J=1,...,M.
C
        DO 100 J=1,M
        P(1,J)=(-E11*U(1,J)+E18*U(2,J)-E9*U(3,J)+E2*U(4,J))/E6HX1
 100    P(N,J)=(E11*U(N,J)-E18*U(N-1,J)+E9*U(N-2,J)-E2*U(N-3,J))/E6HXN
        E6HY1=E6*(Y(2)-Y(1))
        E6HYM=E6*(Y(M)-Y(M-1))
C
C       COMPUTE THE APPROXIMATION TO  DU(X,Y)/DY  AT
C       X=X(I),I=1,...,N,  Y=Y(1)
C       X=X(I),I=1,...,N  Y=Y(M)
C
        DO 200 I=1,N
        Q(I,1)=(-E11*U(I,1)+E18*U(I,2)-E9*U(I,3)+E2*U(I,4))/E6HY1
 200    Q(I,M)=(E11*U(I,M)-E18*U(I,M-1)+E9*U(I,M-2)-E2*U(I,M-3))/E6HYM
C
C       COMPUTE THE APPROXIMATION TO  D2U(X,Y)/DXDY AT THE FOUR CORNERS
C       OF THE RECTANGLE.
C
        S(1,1)=(-E11*Q(1,1)+E18*Q(2,1)-E9*Q(3,1)+E2*Q(4,1))/E6HX1
        S(N,1)=(E11*Q(N,1)-E18*Q(N-1,1)+E9*Q(N-2,1)-E2*Q(N-3,1))
     $          /E6HXN
        S(1,M)=(E11*P(1,M)-E18*P(1,M-1)+E9*P(1,M-2)-E2*P(1,M-3))/E6HYM
        S(N,M)=(E11*P(N,M)-E18*P(N,M-1)+E9*P(N,M-2)-E2*P(N,M-3))/E6HYM
C*******END OF COMPUTING BOUNDARY DERIVATIVES.
C
C*******COMPUTE APPROXIMATIONS TO THE DERIVATIVES OF U, AT THE INTERIOR
C       GRID POINTS OF THE MESH, USING CUBIC SPLINE INTERPOLATION.
C
C       SET UP THE BIDIAGONAL MATRIX, WHICH IS THE RESULT OF PERFORMING
C       GAUSSIAN ELIMINATION ON THE TRIDIAGONAL MATRIX WHOSE ELEMENTS
C       ARE FUNCTIONS OF X.
C       ALPHA=DIAGONAL
C       H=OFF-DIAGONAL
C       AT THE SAME TIME STORE THE QUANTITIES H(I-1)/H(I),I=2,...,N-1
C       WHICH ARE REQUIRED FOR SETTING UP THE RIGHT HAND SIDE, AND ALSO
C       STORE THE QUANTITIES H(I)/ALPHA(I-1),I=3,...,N-1, WHICH
C       APPEAR IN THE FORWARD AND BACKWARD SUBSTITUTION.
C
        H(1)=X(2)-X(1)
        H(2)=X(3)-X(2)
        TH(2)=H(1)/H(2)
        ALPHA(2)=2.D0*(X(3)-X(1))
C
        NM1=N-1
        DO 1I=3,NM1
        IP1=I+1
        IM1=I-1
        H(I)=X(IP1)-X(I)
        TH(I)=H(IM1)/H(I)
        T(I)=H(I)/ALPHA(IM1)
        ALPHA(I)=2.D0*(X(IP1)-X(IM1))-H(I-2)*T(I)
 1      CONTINUE
C
        DO 11J=1,M
C
C       SET UP THE TRANSFORMED RIGHT HAND SIDE IN PREPARATION
C       FOR COMPUTING P
C
C       BETA=TRANSFORMED RIGHT HAND SIDE
C
        BETA(2)=3.D0*(TH(2)*(U(3,J)-U(2,J))+(U(2,J)-U(1,J))/TH(2))
     $          -H(2)*P(1,J)
        DO 2 I=3,NM1
        BETA(I)=3.D0*(TH(I)*(U(I+1,J)-U(I,J))+(U(I,J)-U(I-1,J))/
     $          TH(I))-BETA(I-1)*T(I)
 2      CONTINUE
        BETA(NM1)=BETA(NM1)-H(N-2)*P(N,J)
C
C       BEGIN BACK SUBSTITUTION FOR P
C
        P(NM1,J)=BETA(NM1)/ALPHA(NM1)
        NM3=N-3
        DO 3 IB=1,NM3
        I=NM1-IB
        P(I,J)=(BETA(I)-H(I-1)*P(I+1,J))/ALPHA(I)
 3      CONTINUE
C       END OF COMPUTING P
 11     CONTINUE
C
C       SET UP THE TRANSFORMED RIGHT HAND SIDE IN PREPARATION
C       FOR COMPUTING S ALONG THE LINES Y=Y(1) AND Y=Y(M)
C
        MM1=M-1
        DO 12J=1,M,MM1
        BETA(2)=3.D0*(TH(2)*(Q(3,J)-Q(2,J))+(Q(2,J)-Q(1,J))/TH(2))
     $          -H(2)*S(1,J)
        DO 4I=3,NM1
        BETA(I)=3.D0*(TH(I)*(Q(I+1,J)-Q(I,J))+(Q(I,J)-Q(I-1,J))/
     $          TH(I))-BETA(I-1)*T(I)
 4      CONTINUE
        BETA(NM1)=BETA(NM1)-H(N-2)*S(N,J)
C
C       BEGIN BACK SUBSTITUTION FOR S
C
        S(NM1,J)=BETA(NM1)/ALPHA(NM1)
        NM3=N-3
        DO 5 IB=1,NM3
        I=NM1-IB
        S(I,J)=(BETA(I)-H(I-1)*S(I+1,J))/ALPHA(I)
 5      CONTINUE
C       END OF COMPUTING S ALONG Y=Y(1) AND Y=Y(M)
 12     CONTINUE
C
C       SIMILARLY SET UP THE BIDIAGONAL MATRIX, WHICH IS THE RESULT OF
C       PERFORMING GAUSSIAN ELIMINATION ON THE TRIDIAGONAL MATRIX
C       WHOSE ELEMENTS ARE FUNCTIONS OF Y.
        H(1)=Y(2)-Y(1)
        H(2)=Y(3)-Y(2)
        TH(2)=H(1)/H(2)
        ALPHA(2)=2.D0*(Y(3)-Y(1))
        DO 13 I=3,MM1
        IP1=I+1
        IM1=I-1
        H(I)=Y(IP1)-Y(I)
        TH(I)=H(IM1)/H(I)
        T(I)=H(I)/ALPHA(IM1)
        ALPHA(I)=2.D0*(Y(IP1)-Y(IM1))-H(I-2)*T(I)
 13     CONTINUE
C
        DO 14 I=1,N
C
C       SET UP THE TRANSFORMED RIGHT HAND SIDE IN PREPARATION FOR
C       COMPUTING Q
C
        BETA(2)=3.D0*(TH(2)*(U(I,3)-U(I,2))+(U(I,2)-U(I,1))/TH(2))
     $          -H(2)*Q(I,1)
        BETAS(2)=3.D0*(TH(2)*(P(I,3)-P(I,2))+(P(I,2)-P(I,1))/TH(2))
     $           -H(2)*S(I,1)
        DO 15 J=3,MM1
        BETA(J)=3.D0*(TH(J)*(U(I,J+1)-U(I,J))+(U(I,J)-U(I,J-1))/TH(J))
     $          -BETA(J-1)*T(J)
        BETAS(J)=3.D0*(TH(J)*(P(I,J+1)-P(I,J))+(P(I,J)-P(I,J-1))/TH(J))
     $           -BETAS(J-1)*T(J)
 15     CONTINUE
        BETA(MM1)=BETA(MM1)-H(M-2)*Q(I,M)
        BETAS(MM1)=BETAS(MM1)-H(M-2)*S(I,M)
C
C       BEGIN BACK SUBSTITUTION FOR Q AND S
C
        Q(I,MM1)=BETA(MM1)/ALPHA(MM1)
        S(I,MM1)=BETAS(MM1)/ALPHA(MM1)
        MM3=M-3
        DO 16 JB=1,MM3
        J=MM1-JB
        Q(I,J)=(BETA(J)-H(J-1)*Q(I,J+1))/ALPHA(J)
        S(I,J)=(BETAS(J)-H(J-1)*S(I,J+1))/ALPHA(J)
 16     CONTINUE
C       END OF COMPUTING Q AND S
 14     CONTINUE
C*******END OF COMPUTING DERIVATIVES.
        RETURN
        END
