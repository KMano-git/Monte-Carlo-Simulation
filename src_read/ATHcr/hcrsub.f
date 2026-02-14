C**********************************************************************
      SUBROUTINE EINSTN(F,A,LIM)
C      CALCULATION OF OSCILLATOR STRENGTH AND EINSTEIN COEFFICIENT
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 F(40,40)
C
      REAL*8 A(40,40)
      D3=3.0
      DO 101 I=1,LIM-1
      DO 102 J=I+1,LIM
      AI=I
      AJ=J
      X=1.-(AI/AJ)**2
      G=(.9935+.2328/AI-.1296/AI**2)
     *-(.6282-.5598/AI+.5299/AI**2)/(AI*X)
     *+(.3387-1.181/AI+1.470/AI**2)/(AI*X)**2
      IF(I.NE.1.AND.I.NE.2) GO TO 300
      G=1.0785-.2319/X+.02947/X**2
  200 IF(I.NE.1) GO TO 300
      G=1.1330-.4059/X+.07014/X**2
  300 F(I,J)=2.**6/(3.*DSQRT(D3)*3.1416)*(AI/AJ)**3/(2.*AI**2)*G/X**3
      A(J,I)=8.03E9*AI**2/AJ**2*(AI**(-2)-AJ**(-2))**2*F(I,J)
  102 CONTINUE
  101 CONTINUE
      RETURN
      END
C
C#################################################################
C
      SUBROUTINE CLSAHA(TEMP,SAHA)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 SAHA(40)
C
      TE=TEMP*1.1605E4
C
      DO 101 I=1,40
      P=I
      UION=13.595/TEMP/P**2
  101 SAHA(I)=P**2*DEXP(UION)/2.414E15/DSQRT(TE**3)
C
      RETURN
      END
C
C#################################################################
C
      SUBROUTINE RATCOF(TEMP,OSC,SAHA,C,F,S,ALPHA,BETA,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 OSC(40,40),C(40,40),F(40,40),U(40,40)
      REAL*8 SAHA(40),S(40),ALPHA(40),BETA(40),UION(40)
C
      DO 1 I=1,40
      S(I)=0.0
      ALPHA(I)=0.0
      BETA(I)=0.0
      UION(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
      F(I,J)=0.0
      U(I,J)=0.0
    1 CONTINUE
C
      TE=TEMP*1.1605E4
      UH=13.595
C
      DO 101 I=1,40
      P=I
  101 UION(I)=13.595/TEMP/P**2
      DO 102 I=1,40
      DO 102 J=1,40
  102 U(I,J)=UION(I)-UION(J)
C
C*****                                       *********************
C     EXCITATION RATE COFFICIENT FROM JHONSON
C                                 AND VRIENTS & SMEETS
C                     FIT TO PATHAK ET.AL
C*****                                       *********************

      CALL EXCOFF(U,OSC,TEMP,C,F,S,K)
C     IF(ICON.NE.0) GO TO 1000
C
      DO 105 I=1,40
  105 ALPHA(I)=S(I)*SAHA(I)
C
      DO 602 I=1,40
      P=I
      XP=UH/TEMP/P**2
      CALL CLBETA(XP,P,XS)
      BETA(I)=5.197E-14*(UH/TEMP)**.5/P*XS
C
CCC      WRITE(6,*) BETA(I)
C
  602 CONTINUE
C
      RETURN
      END
C
C  #################################################################
C
C     BEGIN OF EXCOFF ROUTINE
C
      SUBROUTINE EXCOFF(U,OSC,TEMP,C,F,S,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 U(40,40),OSC(40,40),C(40,40),F(40,40)
      REAL*8 S(40),ALPHA(40),G(0:2,40)
      DO 1 I=1,40
      S(I)=0.0
      DO 1 J=1,40
      C(I,J)=0.0
    1 F(I,J)=0.0
C
      TE=TEMP*1.1605E4
C
C*********  1 -> J
      I=1
      P=I
      DO 100 J=2,40
      Q=J
      CALL COF1N(U(I,J),OSC(I,J),TE,C1,I,J,K)
      C(I,J)=C1
  100 F(J,I)=P**2/Q**2*DEXP(U(I,J))*C(I,J)
C
C*********  2-10 -> J
      DO 110 I=2,10
      P=I
      DO 110 J=I+1,40
      Q=J
      CALL COFJO(U(I,J),OSC(I,J),TE,CJ,I,J)
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      GG=((P-2.)/8.)**0.25
      C(I,J)=(1.-GG)*CJ+GG*CV
110   F(J,I)=P**2/Q**2*DEXP(U(I,J))*C(I,J)
C
C*********  2 -> 3
C
C*********  I(>11) -> J
C
      DO 120 I=11,39
      P=I
      DO 120 J=I+1,40
      Q=J
      CALL COFVR(U(I,J),OSC(I,J),TEMP,CV,I,J)
      C(I,J)=CV
  120 F(J,I)=P**2/Q**2*DEXP(U(I,J))*C(I,J)
C
C*********  S  1 ->
      I=1
      CALL COFJS(TE,S1,I)
      S(I)=S1
C
C*********  S  2-10 ->
      DO 210 I=2,10
      CALL COFJS(TE,SJ,I)
      CALL COFVS(TEMP,SV,I)
      GGG=((P-2.)/8.)**0.25
  210 S(I)=(1-GGG)*SJ+GGG*SV
C
C*********  S  I(>11) ->
      DO 220 I=11,40
      CALL COFVS(TEMP,SV,I)
  220 S(I)=SV
C
      RETURN
      END
C
C***** SUBROUTINS ****************************************
C
C *******    C(1,J)  FROM   * JOHNSON *    *******
C
      SUBROUTINE COF1N(U,OSC,TE,C,I,J,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y1,Z1
C
      P=I
      Q=J
      X=1-P**2/Q**2
      Y=U
C**
      R=0.45*X
C     IF (K.EQ.1) THEN
C     RX=-0.95
      IF(J.EQ.2) THEN
      RX=-0.95
      ELSE IF(J.EQ.3) THEN
      RX=-0.95
      ELSE IF(J.EQ.4) THEN
      RX=-0.95
      ELSE IF(J.EQ.5) THEN
      RX=-0.95
      ELSE IF(J.EQ.6) THEN
      RX=-0.94
      ELSE IF(J.EQ.7) THEN
      RX=-0.93
      ELSE IF(J.EQ.8) THEN
      RX=-0.92
      ELSE IF(J.EQ.9) THEN
      RX=-0.91
      ELSE
      RX=-0.9
      END IF
C**
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X-0.603/X**2)
      C1=2*P**2/X
      B=B-A*DLOG(C1)
      Y1=-Y
      Z1=-Z
      CALL DEXPI(Y1,E1Y,ICON)
      CALL DEXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
C
      E2Y=DEXP(-Y)-Y*E1Y
      E2Z=DEXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y+RX*(1/Z+0.5)*E1Z
      E2=E2Y/Y+RX*E2Z/Z
C
      C=1.093E-10*DSQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)
C
      RETURN
      END
C
C *******    C(I,J)  FROM  JOHNSON    *******
C
      SUBROUTINE COFJO(U,OSC,TE,C,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 Y1,Z1
C
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
      Q=J
      X=1-P**2/Q**2
      Y=U
      R=1.94/P**1.57*X
      Z=R+Y
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      C1=2*P**2/X
      B=B-A*DLOG(C1)
      Y1=-Y
      Z1=-Z
      CALL DEXPI(Y1,E1Y,ICON)
      CALL DEXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
C
      E2Y=DEXP(-Y)-Y*E1Y
      E2Z=DEXP(-Z)-Z*E1Z
      E1=(1/Y+0.5)*E1Y-(1/Z+0.5)*E1Z
      E2=E2Y/Y-E2Z/Z
C
      C=1.093E-10*DSQRT(TE)*P**2/X*Y**2*(A*E1+B*E2)
C
      RETURN
      END

C *******    C(I,J)  FROM  VRIENS & SMEETS   *******
C
      SUBROUTINE COFVR(U,OSC,TEMP,C,I,J)
      IMPLICIT REAL*8 (A-H,O-Z)
      UH=13.595
      P=I
      BN=1.4/P*DLOG(P)-0.7/P-0.51/P**2+1.16/P**3-0.55/P**4
      Q=J
      S1=Q-P
      X=1-P**2/Q**2
      A=2*P**2*OSC/X
      B=4*P**4/Q**3/X**2*(1+1.3333/X+BN/X**2)
      DELTA=DEXP(-B/A)+0.06*S1**2/P**2/Q
      G1=1+TEMP/UH*P**3
      G2=3+11*S1**2/P**2
      G3=6+1.6*S1*Q+0.3/S1**2+0.8*Q**1.5/S1**0.5*(S1-0.6)
      GAMMA=UH*DLOG(G1)*G2/G3
      C1=1.6E-7*TEMP**0.5/(TEMP+GAMMA)*DEXP(-U)
      C2=0.3*TEMP/UH+DELTA
C
      C=C1*(A*DLOG(C2)+B)
C
      RETURN
      END
C
C *******    S(I)  FROM  JOHNSON   *******
C
      SUBROUTINE COFJS(TE,S,I)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 G(0:2,40)
      REAL*8 Y1,Z1
C
      G(0,1)= 1.1330
      G(1,1)=-0.4059
      G(2,1)= 0.07014
      G(0,2)= 1.0785
      G(1,2)=-0.2319
      G(2,2)= 0.02947
      DO 350 N=3,40
      G(0,N)=0.9935+0.2328/N-0.1296/N**2
      G(1,N)=-0.6282/N+0.5598/N**2-0.5299/N**3
  350 G(2,N)=0.3887/N**2-1.181/N**3+1.470/N**4
C
      IF (I.EQ.1) THEN
C
      P=1.0
      Y=1.57770E5/TE
      R=0.45
      Z=R+Y
C     RX=-0.74
      RX=-0.59
C
      A=0.0
      DO 223 K=0,2
  223 A=A+G(K,1)/(K+3)
      A=A*1.9603*P
C
      B=0.66667*P**2*(5-0.603)
      C1=2*P**2
      B=B-A*DLOG(C1)
      Y1=-Y
      Z1=-Z
      CALL DEXPI(Y1,E1Y,ICON)
      CALL DEXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
C
      E2Y=DEXP(-Y)-Y*E1Y
      E2Z=DEXP(-Z)-Z*E1Z
      E1=1/Y*E1Y+RX*1/Z*E1Z
      E0Y=DEXP(-Y)/Y
      E0Z=DEXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY+RX*EGZ
C
      S=1.093E-10*DSQRT(TE)*P**2*Y**2*(A*E1+B*E2)
C
      ELSE
C
      P=I
      BN=(4.0-18.63/P+36.24/P**2-28.09/P**3)/P
C     X=1.0
      Y=1.57770E5/TE/P**2
      R=1.94/P**1.57
      Z=R+Y
C *                                 AN
      A=0.0
      DO 222 K=0,2
  222 A=A+G(K,I)/(K+3)
      A=A*1.9603*P
C *
      B=0.66667*P**2*(5+BN)
      C1=2*P**2
      B=B-A*DLOG(C1)
      Y1=-Y
      Z1=-Z
      CALL DEXPI(Y1,E1Y,ICON)
      CALL DEXPI(Z1,E1Z,ICON)
      E1Y=-E1Y
      E1Z=-E1Z
C
      E2Y=DEXP(-Y)-Y*E1Y
      E2Z=DEXP(-Z)-Z*E1Z
      E1=1/Y*E1Y-1/Z*E1Z
      E0Y=DEXP(-Y)/Y
      E0Z=DEXP(-Z)/Z
      EGY=E0Y-2*E1Y+E2Y
      EGZ=E0Z-2*E1Z+E2Z
      E2=EGY-EGZ
C
      S=1.093E-10*DSQRT(TE)*P**2*Y**2*(A*E1+B*E2)
C
      END IF
      RETURN
      END
C
C
C *******    S(I)  FROM  VRIENS & SMEETS   *******
C
      SUBROUTINE COFVS(TEMP,S,I)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      P=I
      UI=13.595/TEMP/P**2
      UIZ=UI**2.33+4.38*UI**1.72+1.32*UI
      S=9.56E-6/TEMP**1.5*DEXP(-UI)/UIZ
C
      RETURN
      END
C
C~~~~~END~OF~EXCOFF~ROUTINE~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
C ######################################################################
C
      SUBROUTINE CLBETA(XP,P,S)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PP,XPP,A,B,EPSA,EPSR,V,GAUN
      REAL*8 ST,EN,STP
      EXTERNAL GAUN
      COMMON PP,XPP
      S=0.0
      II=INT(P)
      PP=P
      XPP=XP
      ST=0.0
      EN=20.0
      ND=100
      STP=(EN-ST)/DBLE(ND)
      DO 1 I=1,ND
      A=STP*(I-1)
      B=A+STP
      N=16
      CALL DGAUSP(GAUN,A,B,N,V)
      S=S+DBLE(V)
    1 CONTINUE
CCC      WRITE (6,*) II,S
      RETURN
      END
C
C##################################################################
C
      SUBROUTINE MOLCOF(TEMP,H2)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 H2(40),S(40),C(40),A(40),Z2(40),Z4(40),ZZ(40)
      REAL*8 M2(40),M,JL,RR(40),B(40),CN(40)
      REAL*8 X1,X2
C
      JL=6.242E18
      TE=TEMP*1.60210E-19
      PAI=3.14159
      A0=5.29167E-11
      R=13.6
      R0=0.01265
      M=9.1091E-31
C
C     ***********  M**2  ************
C
      M2(2)=(0.0911+0.0547)*0.69*0.89
      M2(3)=4.07E-3*1.74
      M2(4)=3.09E-4*2.39
      M2(5)=4.67E-5*2.73
C     M2(6)=2.23E-5*2.83
      DO 50 I=6,40
      P=I
      M2(I)=11.67/P**6.961
   50 CONTINUE
C
C     ***********   S    ************
C
      S(2)=1.2347
      S(3)=1.3469
      S(4)=1.4845
      S(5)=1.5116
C     S(6)=1.49239
      S(6)=1.5432
      DO 49 I=6,8
      S(I)=1.5388+0.00289*(I-6)
   49 CONTINUE
      DO 51 I=9,12
      S(I)=1.5422
   51 CONTINUE
      DO 151 I=13,40
      S(I)=1.5408-0.0016975*(I-13)
  151 CONTINUE
C
C     ***********   R    ************
C
      DO 52 I=2,40
      RR(I)=1.70E-2
   52 CONTINUE
C
C     **********    C    ************
C
      C(2)=-1.2031
      C(3)=1.9316
      C(4)=13.1813
      C(5)=27.9954
C     C(6)=20.0966
      DO 3 I=6,40
      P=I
      C(I)=66.3884*DLOG(0.3049*P)
    3 CONTINUE
C
C     **********    CN   ************
C
      CN(2)=0.57*0.69*1.48*0.89E-21
      CN(3)=0.5975E-22*1.74
      CN(4)=0.06775E-22*2.39
      CN(5)=0.015225E-22*2.73
      DO 13 I=6,40
      P=I
      CN(I)=6.693/P**5.99*1.0E-20
   13 CONTINUE
      DO 54 I=2,40
      P=I
      IF (I.EQ.2) THEN
CC    ES=17.0*1.6021E-19
CC    E0=17.0*1.6021E-19
      ES=16.565*1.6021E-19
      E0=16.565*1.6021E-19
      ELSE
C     ES=17.0*1.6021E-19
      ES=(13.6*(1.0-(1.0/(P+1.0)**2))+4.476)*1.6021E-19
      E0=28.0*1.6021E-19
      END IF
      R0=RR(I)
      X1=-E0/TE
      X2=-(R0*JL+1.0/TE)*E0
      CALL DEXPI(X1,E1,ICON)
      CALL DEXPI(X2,E2,ICON)
      E1=-E1
      E2=-E2
C
      Z1=DEXP(-E0/TE)*DLOG(JL*E0)*TE+E1*TE
      Z2(I)=S(I)*DEXP(-(R0*JL+1.0/TE)*E0)*DLOG(JL*E0)
     *      *(1.0/(R0*JL+1.0/TE))
     *      +S(I)*E2*(1.0/(R0*JL+1.0/TE))
      Z3=C(I)*TE*DEXP(-E0/TE)
      Z4(I)=S(I)*C(I)*(1.0/(R0*JL+1.0/TE))*DEXP(-(R0*JL+1.0/TE)*E0)
      ZZ(I)=(Z1-Z2(I)+Z3-Z4(I))/JL
   54 CONTINUE
      DO 55 I=2,40
      P=I
      A(I)=4*PAI*A0**2*R*M2(I)*4*PAI*2.0*DSQRT((M/(2.0*PAI*TE))**3)/M**2
      B(I)=4*PAI*2.0*DSQRT((M/(2.0*PAI*TE))**3)/M**2
     *     *(TE*DEXP(-ES/TE)*(ES+TE)-TE*DEXP(-E0/TE)*(E0+TE))
      H2(I)=(A(I)*ZZ(I)+B(I)*CN(I))*1.0E6
CCC   WRITE(6,*) H2(I)
   55 CONTINUE
C
C
      RETURN
      END
C
C######################################################################
C
      SUBROUTINE POPCOF(DENSEL,SAHA,C,F,H2,S,A,ALPHA,BETA,
     *                 LUP,LIM,R0,R1,R2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 C(40,40),F(40,40),A(40,40),W(40,40),H2(40)
      REAL*8 S(40),ALPHA(40),BETA(40),SAHA(40),HH2(40)
      DIMENSION IP(40)
      REAL*8 BLAX(40),VW(40),WA(40,40),WK(40)
      REAL*8 R0(40),R1(40),R2(40),NA,NI,NMOL
      REAL*8 EPSZ
C
      DO 201 K=2,LUP-1
      DO 202 L=2,K
  202 W(K,L)=C(L,K)*DENSEL
      SUMF=0.
      DO 301 I=1,K-1
  301 SUMF=SUMF+F(K,I)
      SUMC=0.
      DO 302 I=K+1,LIM
  302 SUMC=SUMC+C(K,I)
      SUMA=0.
      DO 303 I=1,K-1
  303 SUMA=SUMA+A(K,I)
      W(K,K)=-(DENSEL*(SUMF+SUMC+S(K))+SUMA)
      DSF = DENSEL*(SUMF+SUMC+S(K))
      DO 203 L=K+1,LUP
  203 W(K,L)=DENSEL*F(L,K)+A(L,K)
  201 CONTINUE
      DO 211 L=2,LUP-1
  211 W(LUP,L)=C(L,LUP)*DENSEL
      SUMF=0.
      DO 311 I=1,LUP-1
  311 SUMF=SUMF+F(LUP,I)
      SUMC=0.0
      DO 313 I=LUP+1,LIM
  313 SUMC=SUMC+C(LUP,I)
      SUMA=0.
      DO 312 I=1,LUP-1
  312 SUMA=SUMA+A(LUP,I)
      W(LUP,LUP)=-(DENSEL*(SUMF+SUMC+S(LUP))+SUMA)
      DO 550 K=2,LUP
      SUMF=0.0
      DO 500 I=LUP+1,LIM
  500 SUMF=SUMF+SAHA(I)*F(I,K)
      SUMAS=0.0
      DO 501 I=LUP+1,LIM
  501 SUMAS=SUMAS+SAHA(I)*A(I,K)
C
      W(K,LUP+1)=-((DENSEL*SUMF+SUMAS)+(DENSEL*ALPHA(K)
     *            +BETA(K)))
C
      W(K,LUP+2)=-C(1,K)
      W(K,LUP+3)=-H2(K)
  550 CONTINUE
      DO 3001 II=LUP,LUP+2
      DO 402 I=1,LUP-1
      DO 402 J=1,LUP+2
  402 WA(I,J)=W(I+1,J+1)
      DO 3000 J=1,LUP-1
      BLAX(J)=WA(J,II)
 3000 CONTINUE
      EPSZ=0.0
C
      CALL DLAX(WA,40,LUP-1,BLAX,EPSZ,1,IS,VW,IP,ICON)
C
      DO 3010 J=1,LUP-1
      IF(II.EQ.LUP) THEN
      R0(J+1)=BLAX(J)
      ELSE IF (II.EQ.LUP+1) THEN
      R1(J+1)=BLAX(J)
      ELSE
      R2(J+1)=BLAX(J)
      END IF
 3010 CONTINUE
 3001 CONTINUE
C
      RETURN
      END
C##################################################################
C
      SUBROUTINE CLPOP(C,F,S,A,LUP,R1,DENSEL,SCR)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R0(40),R1(40),R2(40),C(40,40),F(40,40),S(40),A(40,40)
      SCR=S(1)*DENSEL
      DO 100 I=2,LUP
      SCR=SCR+C(1,I)*DENSEL-(F(I,1)*DENSEL+A(I,1))*R1(I)*DENSEL
  100 CONTINUE
      SCR=SCR/DENSEL
      RETURN
      END
C
C##################################################################
C
      SUBROUTINE CLPOP1(LUP,R0,R1,R2,PD,DENSEL)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 PD(40),SAHA(40)
C     REAL*8 PG(40)
      REAL*8 R0(40),R1(40),R2(40)
      REC=0.
      DO 100 I=2,LUP
      P=I
      PD(I)=(R0(I)+R1(I)+R2(I))*DENSEL
C     PG(I)=PD(I)/2./P**2
  100 CONTINUE

      RETURN
      END
C***********************************************************************
      SUBROUTINE CLPOP2(R0,R1,R2,DENSEL,C,F,S,A,SS,LA,LB,LG)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 R0(40),R1(40),R2(40)
      REAL*8 C(40,40),F(40,40),S(40),A(40,40),LA,LB,LG
C
      SUMS=0.0
      SUMF=0.0
      SUMA=0.0
      SUML=0.0
C
      DO 100 I=2,35
      DO 101 II=36,40
      SUMS=SUMS+R1(I)*DENSEL*C(I,II)
  101 CONTINUE
      SUMS=SUMS+R1(I)*DENSEL*S(I)
      SUMA=SUMA+R1(I)*A(I,1)
      SUMF=SUMF+R1(I)*DENSEL*F(I,1)
  100 CONTINUE
C
      SS=SUMS
      LA=R1(2)*A(2,1)
      LB=R1(3)*A(3,1)
      LG=R1(4)*A(4,1)
C
      RETURN
      END
C
C##################################################################
      SUBROUTINE DGAUSP ( FUNC, A, B, N, V)
*****************************************************
*        Legendre-Gauss quadratuer formula          *
*                    over (A,B)                     *
*      COPYRIGHT : M.Mori  JUNE 30 1989  V.1        *
*   ---- input parameters ----                      *
*     FUNC = name of the function subprogram for    *
*            the integrand                          *
*     A = lower bound of integration                *
*     B = upper bound of integration                *
*     N = number of points of the formula           *
*         not larger than 16                        *
*   ---- output parameter ----                      *
*     V = result of integration                     *
*****************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
*
      REAL*8 FUNC
      EXTERNAL FUNC
      REAL*8 C(0:16,16),W(0:16,16)
*
      DATA (C(I, 1),I=0, 0) /
     $   0.0000000000000000D+00 /
      DATA (W(I, 1),I=0, 0) /
     $   0.2000000000000000D+01 /
      DATA (C(I, 2),I=1, 1) /
     $   0.5773502691896257D+00 /
      DATA (W(I, 2),I=1, 1) /
     $   0.9999999999999998D+00 /
      DATA (C(I, 3),I=0, 1) /
     $   0.0000000000000000D+00,
     $   0.7745966692414832D+00 /
      DATA (W(I, 3),I=0, 1) /
     $   0.8888888888888888D+00,
     $   0.5555555555555558D+00 /
      DATA (C(I, 4),I=1, 2) /
     $   0.3399810435848563D+00,
     $   0.8611363115940525D+00 /
      DATA (W(I, 4),I=1, 2) /
     $   0.6521451548625459D+00,
     $   0.3478548451374539D+00 /
      DATA (C(I, 5),I=0, 2) /
     $   0.0000000000000000D+00,
     $   0.5384693101056831D+00,
     $   0.9061798459386639D+00 /
      DATA (W(I, 5),I=0, 2) /
     $   0.5688888888888888D+00,
     $   0.4786286704993666D+00,
     $   0.2369268850561892D+00 /
      DATA (C(I, 6),I=1, 3) /
     $   0.2386191860831969D+00,
     $   0.6612093864662644D+00,
     $   0.9324695142031520D+00 /
      DATA (W(I, 6),I=1, 3) /
     $   0.4679139345726909D+00,
     $   0.3607615730481386D+00,
     $   0.1713244923791703D+00 /
      DATA (C(I, 7),I=0, 3) /
     $   0.0000000000000000D+00,
     $   0.4058451513773972D+00,
     $   0.7415311855993944D+00,
     $   0.9491079123427584D+00 /
      DATA (W(I, 7),I=0, 3) /
     $   0.4179591836734694D+00,
     $   0.3818300505051189D+00,
     $   0.2797053914892767D+00,
     $   0.1294849661688699D+00 /
      DATA (C(I, 8),I=1, 4) /
     $   0.1834346424956498D+00,
     $   0.5255324099163289D+00,
     $   0.7966664774136267D+00,
     $   0.9602898564975361D+00 /
      DATA (W(I, 8),I=1, 4) /
     $   0.3626837833783622D+00,
     $   0.3137066458778873D+00,
     $   0.2223810344533746D+00,
     $   0.1012285362903763D+00 /
      DATA (C(I, 9),I=0, 4) /
     $   0.0000000000000000D+00,
     $   0.3242534234038089D+00,
     $   0.6133714327005903D+00,
     $   0.8360311073266357D+00,
     $   0.9681602395076260D+00 /
      DATA (W(I, 9),I=0, 4) /
     $   0.3302393550012598D+00,
     $   0.3123470770400028D+00,
     $   0.2606106964029356D+00,
     $   0.1806481606948573D+00,
     $   0.8127438836157435D-01 /
      DATA (C(I,10),I=1, 5) /
     $   0.1488743389816312D+00,
     $   0.4333953941292473D+00,
     $   0.6794095682990244D+00,
     $   0.8650633666889845D+00,
     $   0.9739065285171717D+00 /
      DATA (W(I,10),I=1, 5) /
     $   0.2955242247147527D+00,
     $   0.2692667193099963D+00,
     $   0.2190863625159819D+00,
     $   0.1494513491505806D+00,
     $   0.6667134430868799D-01 /
      DATA (C(I,11),I=0, 5) /
     $   0.0000000000000000D+00,
     $   0.2695431559523450D+00,
     $   0.5190961292068118D+00,
     $   0.7301520055740493D+00,
     $   0.8870625997680953D+00,
     $   0.9782286581460569D+00 /
      DATA (W(I,11),I=0, 5) /
     $   0.2729250867779006D+00,
     $   0.2628045445102467D+00,
     $   0.2331937645919905D+00,
     $   0.1862902109277342D+00,
     $   0.1255803694649046D+00,
     $   0.5566856711617373D-01 /
      DATA (C(I,12),I=1, 6) /
     $   0.1252334085114689D+00,
     $   0.3678314989981802D+00,
     $   0.5873179542866174D+00,
     $   0.7699026741943046D+00,
     $   0.9041172563704747D+00,
     $   0.9815606342467192D+00 /
      DATA (W(I,12),I=1, 6) /
     $   0.2491470458134029D+00,
     $   0.2334925365383548D+00,
     $   0.2031674267230659D+00,
     $   0.1600783285433463D+00,
     $   0.1069393259953185D+00,
     $   0.4717533638651187D-01 /
      DATA (C(I,13),I=0, 6) /
     $   0.0000000000000000D+00,
     $   0.2304583159551348D+00,
     $   0.4484927510364469D+00,
     $   0.6423493394403402D+00,
     $   0.8015780907333098D+00,
     $   0.9175983992229779D+00,
     $   0.9841830547185881D+00 /
      DATA (W(I,13),I=0, 6) /
     $   0.2325515532308739D+00,
     $   0.2262831802628972D+00,
     $   0.2078160475368886D+00,
     $   0.1781459807619456D+00,
     $   0.1388735102197873D+00,
     $   0.9212149983772848D-01,
     $   0.4048400476531587D-01 /
      DATA (C(I,14),I=1, 7) /
     $   0.1080549487073436D+00,
     $   0.3191123689278898D+00,
     $   0.5152486363581540D+00,
     $   0.6872929048116854D+00,
     $   0.8272013150697650D+00,
     $   0.9284348836635735D+00,
     $   0.9862838086968123D+00 /
      DATA (W(I,14),I=1, 7) /
     $   0.2152638534631578D+00,
     $   0.2051984637212956D+00,
     $   0.1855383974779379D+00,
     $   0.1572031671581936D+00,
     $   0.1215185706879031D+00,
     $   0.8015808715976016D-01,
     $   0.3511946033175195D-01 /
      DATA (C(I,15),I=0, 7) /
     $   0.0000000000000000D+00,
     $   0.2011940939974345D+00,
     $   0.3941513470775634D+00,
     $   0.5709721726085388D+00,
     $   0.7244177313601700D+00,
     $   0.8482065834104272D+00,
     $   0.9372733924007058D+00,
     $   0.9879925180204854D+00 /
      DATA (W(I,15),I=0, 7) /
     $   0.2025782419255613D+00,
     $   0.1984314853271116D+00,
     $   0.1861610000155623D+00,
     $   0.1662692058169940D+00,
     $   0.1395706779261542D+00,
     $   0.1071592204671719D+00,
     $   0.7036604748810814D-01,
     $   0.3075324199611710D-01 /
      DATA (C(I,16),I=1, 8) /
     $   0.9501250983763744D-01,
     $   0.2816035507792589D+00,
     $   0.4580167776572274D+00,
     $   0.6178762444026437D+00,
     $   0.7554044083550029D+00,
     $   0.8656312023878317D+00,
     $   0.9445750230732326D+00,
     $   0.9894009349916499D+00 /
      DATA (W(I,16),I=1, 8) /
     $   0.1894506104550686D+00,
     $   0.1826034150449236D+00,
     $   0.1691565193950025D+00,
     $   0.1495959888165767D+00,
     $   0.1246289712555339D+00,
     $   0.9515851168249285D-01,
     $   0.6225352393864789D-01,
     $   0.2715245941175410D-01 /
*
      IF ((N. LT. 1) .OR. (N .GT. 16)) GO TO 901
*
      C1 = (B + A) / 2
      C2 = (B - A) / 2
*
      IF (MOD(N,2) .EQ. 0) THEN
        NH = N / 2
        V = 0
      ELSE
        NH = (N - 1) / 2
        V = W(0,N) * FUNC(C1)
      END IF
*
      DO 10 I = 1, NH
        V = V + W(I,N) * (FUNC(C1 + C2 * C(I,N))
     $                  + FUNC(C1 - C2 * C(I,N)))
   10 CONTINUE
*
      V = C2 * V
*
      RETURN
*
  901 CONTINUE
      WRITE (6,2001) N
 2001 FORMAT (' (SUBR.DGAUSP) invalid argument',
     $        ' N =',I4)
      RETURN
*
      END




C*    *** LAX    ***   *   *********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
C     *                                                                *
C     *   A22-11-0101    LAX,DLAX,QLAX            VERSION-1            *
C     *                                                                *
C     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
C     *             K.ITO               1976                           *
C     *                                                                *
C     *   'LAX' SOLVES THE LINEAR EQUATIONS BY CROUT'S METHOD.         *
C     *   AX=B,WHERE A IS REGULAR MATRIX WITH REAL*8 COEFFICIENT.        *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL LAX(A,K,N,B,EPSZ,ISW,IS,VW,IP,ICON)                *
C     *          A   ....GIVEN REGULAR COEFFICIENT MATRIX.             *
C     *              ....RESULTANT LU-DECOMPOSED ENTRYS.               *
C     *                    2 REAL*8AL ARRAY AS A(K,N),SIZE IS N.    *
C     *          K   ....GIVEN ADJUSTABLE REAL*8 FOR ARRAY A.       *
C     *          N   ....GIVEN ORDER OF MATRIX A.                      *
C     *              ....RESULTANT SOLUTION VECTOR.                    *
C     *          B   ....GIVEN CONSTANT VECTOR.                        *
C     *              ....RESULTANT SOLUTION VECTOR.                    *
C     *                    1 REAL*8AL ARRAY,SIZE IS N.              *
C     *          EPSZ....GIVEN VALUE WHICH IS REFERED WHEN  RELATIVE   *
C     *                    ZERO IS DETECTED,POSITIVE.                  *
C     *                    DEFAULT VALUE IS SETTED IF ZERO ASSIGNED    *
C     *          ISW ....GIVEN CONTROL INFORMATION.                    *
C     *                    IF 1,SOLVE EQUATIONS ENTIRLY.               *
C     *                    IF 2,SOLVE EQUATIONS WITH LAST              *
C     *                    LU-DECOMPOSED ENTRYS.                       *
C     *          IS  ....RESULTANT SIGN OF DETERMINANT A.              *
C     *          IP  ....AUXILIARY 1 REAL*8ED ARRAY,SIZE IS N.      *
C     *                    TRANSPOSITION VECTOR WHICH REPRESENTS       *
C     *                    ROW-EXCHANGING BY PARTIAL PIVOTING.         *
C     *          ICON....RESULTANT CONDITION CODE.                     *
C     *                                                                *
C     *   METHOD                                                       *
C     *        CROUT'S METHOD                                          *
C     *                                                                *
C     *   SLAVE SUBROUTINE                                             *
C     *        ALU,LUX,AMACH                                           *
C     *                                                                *
C     *   REFERENCE CODE                                               *
C     *                                                                *
C     *   REMARK                                                       *
C     *        WHEN PARTIAL PIVOTING,DO ROW EQUILIBRATED PARTIAL       *
C     *        PIVOTING                                                *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     ******************************************************************
      SUBROUTINE DLAX(A,K,N,B,EPSZ,ISW,IS,VW,IP,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(K,N),B(N),VW(N)
      DIMENSION IP(N)
!      character mcode(6)*2                       !  DIMENSION MCODE(6)
!      DATA MCODE/2HA2,2H2-,2H11,2H-0,2H10,2H1 /  ! <==  Warning 963
      character(2) :: mcode(6)
      data mcode /"A2","2-","11","-0","10","1"/
      IF(ISW.EQ.1) GO TO 1000
      IF(ISW.EQ.2) GO TO 1100
      ICON=30000
      GO TO 8000
 1000 CALL DALU(A,K,N,EPSZ,IP,IS,VW,ICON)
      IF(ICON.NE.0) GO TO 8000
 1100 CALL DLUX(B,A,K,N,1,IP,ICON)
 8000 CALL MGSSL(ICON,MCODE)
      RETURN
      END
C
C
C
C*    *** ALU    ***   *   *********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1982,1983        *
C     *                                                                *
C     *   A22-11-0202    ALU,DALU,QALU               VERSION-1         *
C     *                                                                *
C     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
C     *             K.ITO               1976                           *
C     *                                                                *
C     *   'ALU' DECOMPOSES THE NON SINGULAR REAL*8 MATRIX BY CROUT'S     *
C     *   METHOD(LU-DECOMPOSITION).                                    *
C     *                                                                *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL ALU(A,K,N,EPSZ,IP,IS,VW,ICON)                      *
C     *          A   ....GIVEN NON SINGULAR REAL*8 MATRIX.               *
C     *              ....RESULTANT LU-DECOMPOSED MATRIX L AND U.       *
C     *                    ARRAY A INCLUDES LOWER TRIANGULAR MATRIX L  *
C     *                    ON LOWER TRIANGULAR AREA AND  UNIT UPPER    *
C     *                    TRIANGULAR MATRIX U ON UPPER TRIANGULAR     *
C     *                    AREA WITH NO DIAGONAL ELEMENTS OF U.        *
C     *                    2 REAL*8AL ARRAY AS A(K,N).              *
C     *          K   ....GIVEN ADJUSTABLE REAL*8 FOR ARRAY A.       *
C     *          N   ....GIVEN ORDER OF MATRIX A.                      *
C     *          EPSZ....GIVEN MAXIMUM RELATIVE ERROR.                 *
C     *                    IF ZERO ASSIGNED,BE ASSUMED AS IF DEFAULT   *
C     *                    VALUE WAS ASSIGNED(UNIT ROUND OFF*16).      *
C     *          IP  ....RESULTANT TRANSPOSITION VECTOR WHICH          *
C     *                    REPRESENTS ROW-EXCHANGING BY PARTIAL        *
C     *                    PIVOTING.                                   *
C     *                    1 REAL*8AL ARRAY,SIZE IS N.              *
C     *          IS  ....RESULTANT SIGN OF DETERMINANT A.              *
C     *          VW  ....AUXILIARY 1 REAL*8AL ARRAY,SIZE IS N.      *
C     *          ICON....RESULTANT CONDITION CODE.                     *
C     *                                                                *
C     *   METHOD                                                       *
C     *        CROUT'S METHOD                                          *
C     *                                                                *
C     *   SLAVE SUBROUTINE                                             *
C     *        AMACH                                                   *
C     *                                                                *
C     *   REMARK                                                       *
C     *        IF ASSIGNED 10**(-S) TO EPSZ, IT MEANS THAT IF          *
C     *        CANCELLATION MORE THAN S OCCURED,SO THAT S DIGITS WERE  *
C     *        LOSSED,THEN THE VALUE OF PIVOT IS ASSUMED TO BE ZERO.   *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     ******************************************************************
      SUBROUTINE DALU(A,K,N,EPSZ,IP,IS,VW,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 A(K,N),VW(N)
      DIMENSION IP(N)
!      character mcode(6)*2                       !  DIMENSION MCODE(6)
!      DATA MCODE/2HA2,2H2-,2H11,2H-0,2H20,2H2 /  ! <== Warning 1034
      character(2) :: mcode(6)
      data mcode / "A2","2-","11","-0","20","2" /
      REAL*8    SUM
      REAL*8    DMACH
      EXTERNAL DMACH
C    ------------------------------------------------------------------
C     ENTRY POINT.
C     SET CONDITION CODE TO DEFAULT VALUE.
C    ------------------------------------------------------------------
      ICON=0
      IP(1)=1
C    ------------------------------------------------------------------
C     CHECK ARGUMENTS
C    ------------------------------------------------------------------
      IF(K.LT.N .OR.
     *   EPSZ.LT.0.0D0) GO TO 10
      IF(N.GT.1) GO TO 1000
      IF(N.EQ.1) GO TO 15
C    ------------------------------------------------------------------
C     ERROR  ICON=30000
C    ------------------------------------------------------------------
   10 ICON=30000
      GO TO 8000
   15 IF(A(1,1)) 8000,20,8000
C    ------------------------------------------------------------------
C     ERROR  ICON=2000
C    ------------------------------------------------------------------
   20 ICON=20000
      GO TO 8000
C    ------------------------------------------------------------------
C     DETERMINE SCALING FACTOR FOR EACH ROW.
C    ------------------------------------------------------------------
 1000 DO 25 I=1,N
      VW(I)=DABS(A(I,1))
   25 CONTINUE
      DO 30 J=2,N
      DO 27 I=1,N
      IF(VW(I).LT.DABS(A(I,J))) VW(I)=DABS(A(I,J))
   27 CONTINUE
   30 CONTINUE
C    ------------------------------------------------------------------
C     CHECK IF SCALING FACTOR IS ZERO.
C    ------------------------------------------------------------------
      OGMAX=0.0
      DO 40 I=1,N
      TEMP1=VW(I)
      IF(TEMP1.EQ.0.0) GO TO 20
      IF(TEMP1.GT.OGMAX) OGMAX=TEMP1
      VW(I)=1.0/TEMP1
   40 CONTINUE
C    ------------------------------------------------------------------
C     SET EPS TO EPSZ. IF EPSZ IS EQUAL TO ZERO,SET EPS TO DEFAULT
C     VALUE.
C    ------------------------------------------------------------------
      EPS=EPSZ
      IF(EPSZ.EQ.0.0D0) EPS=DMACH(EPSZ)*16.0D0
      IS=1
C.   ------------------------------------------------------------------
C.    DETERMINE THE PIVOT OF 1ST COLUMN.
C.   ------------------------------------------------------------------
      OG=A(1,1)
      S=ABS(OG*VW(1))
      DO 80 I=2,N
      IF(S.GE.ABS(A(I,1)*VW(I))) GO TO 80
      OG=A(I,1)
      S=ABS(OG*VW(I))
      IP(1)=I
   80 CONTINUE
C.   ------------------------------------------------------------------
C.    BEGIN TRIANGULAR DECOMPOSITION  WITH PARTIAL PIVOTING.
C.   ------------------------------------------------------------------
 1100 DO 180 J=2,N
C.   ------------------------------------------------------------------
C.    EXCHANGE ROWS FOR PARTIAL PIVOTING.
C.   ------------------------------------------------------------------
      LOCTM1=J-1
      IF(IP(J-1).EQ.(J-1)) GO TO 91
      IS=-IS
      M=IP(J-1)
      VW(M)=VW(J-1)
      DO 90 I=1,LOCTM1
      L=J-I
      T=A(LOCTM1,L)
      A(LOCTM1,L)=A(M,L)
      A(M,L)=T
   90 CONTINUE
   91 DO 100 I=1,LOCTM1
      IF(IP(I).EQ.I) GO TO 100
      M=IP(I)
      T=A(I,J)
      A(I,J)=A(M,J)
      A(M,J)=T
  100 CONTINUE
C.   ------------------------------------------------------------------
C.    CHECK IF THE PIVOT IS ZERO.
C.   ------------------------------------------------------------------
      T=A(J-1,J-1)
      IF(T.EQ.0.0  .OR.
     1   ABS(T).LT.ABS(OG*EPS)  .OR.
     2   (OG.EQ.0.0  .AND.
     3    ABS(T).LE.(OGMAX*EPS))) GO TO 20
      VW(J-1)=1.0/T
C.   ------------------------------------------------------------------
C.    LU-DECOMPOSITION BY COLUMNWISE
C.    AT FIRST STAGE,CALCULATE MATRIX U ENTRYS(:U(I,J))
C.   ------------------------------------------------------------------
 1200 A(1,J)=A(1,J)*VW(1)
      IK=1
      IF(J.EQ.2) GO TO 1300
      DO 130 I=2,LOCTM1
      SUM=0.0D0
      M=I-1
      IF(IK.EQ.1) M=1
      I900=I-1
      DO 120 II=1,I900
      SUM=A(I,M )*A(M ,J)+SUM
  120 M=M+IK
      IK=-IK
      A(I,J)=(A(I,J)-SUM)*VW(I)
  130 CONTINUE
C.   ------------------------------------------------------------------
C.    AT STAGE AFTER CALCULATING J-TH COLUMN U ELEMENTS,
C.    CALCULATE J-TH COLUMN L ELEMENTS(:L(I,J)).
C.   ------------------------------------------------------------------
 1300 S=-1.0
      DO 170 I=J,N
      SUM=0.0D0
      U=A(I,J)
      M=J-1
      IF(IK.EQ.1) M=1
      I901=J-1
      DO 140 II=1,I901
      SUM=A(I,M )*A(M ,J)+SUM
  140 M=M+IK
      A(I,J)=U-SUM
      T=ABS(A(I,J)*VW(I))
      IF(S.GE.T) GO TO 160
      S=T
      OG=U
      IP(J)=I
  160 IK=-IK
  170 CONTINUE
  180 CONTINUE
C.   ------------------------------------------------------------------
C.    CHECK IF DIAGONAL ELEMENT A(N,N) IS RELATIVE ZERO OR
C.    ABSOLUTE ZERO.
C.   ------------------------------------------------------------------
 1400 T=A(N,N)
      IF(T.EQ.0.0  .OR.
     1   ABS(T).LT.ABS(OG*EPS)  .OR.
     2  (OG.EQ.0.0.AND.ABS(T).LE.(OGMAX*EPS))) GO TO 20
 8000 CALL MGSSL(ICON,MCODE)
      RETURN
      END
C
C
C
C
C*    *** LUX    ***   *   *********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1982,1983        *
C     *                                                                *
C     *   A22-11-0302    LUX,DLUX,QLUX                VERSION-1        *
C     *                                                                *
C     *   AUTHOR....J.MIKAMI,H.MURAOKA,                                *
C     *             K.ITO               1976                           *
C     *                                                                *
C     *   'LUX' SOLVES THE LINEAR EQUATIONS LUX=PB.  WHERE L AND U ARE *
C     *   LOWER AND UNIT UPPER TRIANGULAR REAL*8 MATRICES,X IS           *
C     *   SOLUTION VECTOR,P IS PERMUTATION MATRIX AND B IS REAL*8        *
C     *   CONSTANT VECTOR.                                             *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL LUX(B,FA,K,N,ISW,IP,ICON)                          *
C     *          B   ....GIVEN CONSTANT VECTOR.                        *
C     *              ....RESULTANT SOLUTION VECTOR.                    *
C     *                    1 REAL*8AL ARRAY,SIZE IS N.              *
C     *          FA  ....GIVEN LU-DECOMPOSED(FACTORIZED) MATRIX WHICH  *
C     *                    CONSIST OF MATRIX L AND U.                  *
C     *                    2 REAL*8AL ARRAY AS FA(K,N).             *
C     *          K   ....GIVEN ADJUSTABLE REAL*8 FOR ARRAY FA.      *
C     *          N   ....GIVEN ORDER OF THE MARIX A.                   *
C     *          ISW ....GIVEN CONTROL INFORMATION.                    *
C     *                    IF 1,SOLVE LUX=PB.                          *
C     *                    IF 2,SOLVE LY=PB.                           *
C     *                    IF 3,SOLVE UZ=B.                            *
C     *          IP  ....GIVEN TRANSPOSITION VECTOR WHICH              *
C     *                    REPRESENTS ROW-EXCHANGING BY PARTIAL        *
C     *                    PIVOTING.                                   *
C     *                    1 REAL*8AL ARRAY,SIZE IS N.              *
C     *          ICON....RESULTANT CONDITION CODE.                     *
C     *                                                                *
C     *   METHOD                                                       *
C     *        BACK SUBSTITUTION,FORWARD SUBSTITUTION.                 *
C     *                                                                *
C     *   SLAVE SUBROUTINE                                             *
C     *        NONE.                                                   *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     ******************************************************************
C    ------------------------------------------------------------------
C     FOR FURTHER SIMPLICITY,WE SUBSTITUTE ARRAY FA BY ARRAY A
C     IN BELOW ALGORITHM
C    ------------------------------------------------------------------
      SUBROUTINE DLUX(B,A,K,N,ISW,IP,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8  A(K,N),B(N)
      DIMENSION IP(N)
!      character mcode(6)*2                       !  DIMENSION MCODE(6)
!      DATA MCODE/2HA2,2H2-,2H11,2H-0,2H30,2H2 /
      character(2) :: mcode(6)
      data mcode /"A2","2-","11","-0","30","2"/
      REAL*8     SUM
C    ------------------------------------------------------------------
C     ENTRY POINT
C     SET CONDITION CODE TO DEFAULT VALUE.
C    ------------------------------------------------------------------
      ICON=0
C    ------------------------------------------------------------------
C     CHECK ARGUMENTS
C    ------------------------------------------------------------------
      IF(N.LE.0   .OR.
     1   K.LT.N   .OR.
     2   ISW.LE.0 .OR.
     3   ISW.GE.4) GOTO 9000
      DO 10 I=1,N
      IF(IP(I).LT.I .OR.
     1   IP(I).GT.N) GOTO 9000
   10 CONTINUE
C.   ------------------------------------------------------------------
C.    IF ISW ^= 3 THEN COMPUTE LY=PB BY FORWARD SUBSTITUTION.
C.   ------------------------------------------------------------------
 1000 IF(ISW.EQ.3) GO TO 1100
      IF(A(1,1).NE.0.0) GO TO 15
   14 ICON=20000
      GO TO 8000
   15 IF(IP(1).EQ.1) GO TO 16
      M=IP(1)
      T=B(1)
      B(1)=B(M)
      B(M)=T
   16 B(1)=B(1)/A(1,1)
      IF(N.EQ.1) GO TO 1100
      M=1
      IK=1
      DO 40 I=2,N
      IF(A(I,I).EQ.0.0) ICON=20000
      IF(ICON.NE.0) GO TO 40
      W=A(I,I)
      IF(IP(I).EQ.I) GO TO 20
      M=IP(I)
      T=B(I)
      B(I)=B(M)
      B(M)=T
   20 SUM=0.0D0
      M=I-1
      IF(IK.EQ.1) M=1
      I900=I-1
      DO 30 IJ=1,I900
      SUM=A(I,M )*B(M )+SUM
   30 M=M+IK
      IK=-IK
      B(I)=(B(I)-SUM)/W
   40 CONTINUE
C.   ------------------------------------------------------------------
C.    CHECK IF ANY DIAGONAL ELEMENTS ARE ZERO.
C.   ------------------------------------------------------------------
      IF(ICON.NE.0) GO TO 8000
C    ------------------------------------------------------------------
C     IF ISW^=2 THEN COMPUTE UZ=B BY BACK SUBSTITUTION.
C    ------------------------------------------------------------------
 1100 IF(ISW.EQ.2 .OR.
     1   N.EQ.1) GO TO 8000
      IK=-1
      DO 60 J=2,N
      I=N-J+1
      SUM=0.0D0
      M=N
      IF(IK.EQ.1) M=I+1
      I901=I+1
      DO 50 IJ=I901,N
      SUM=A(I,M )*B(M )+SUM
   50 M=M+IK
      IK=-IK
      B(I)=B(I)-SUM
   60 CONTINUE
C    ------------------------------------------------------------------
C
C    ------------------------------------------------------------------
 8000 CALL MGSSL(ICON,MCODE)
      RETURN
 9000 ICON=30000
      GOTO 8000
      END


C ######################################################################
C
      REAL*8 FUNCTION GAUN(X)
      REAL*8 PP,XPP,U,B,X
      COMMON PP,XPP
      U=X/XPP
      B=PP
      GAUN=(1./(U+1.)+0.1728*(U-1.)/B**(2./3.)/(U+1.)**(5./3.)-0.0496*
     *(U**2+4./3.*U+1.)/B**(4./3.)/(U+1.)**(7./3.))*DEXP(-X)
      RETURN
      END
C
C
C
C*    *** DEXPI   ***   *   ********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
C     *                                                                *
C     *   I11-31-0101    DEXPI                          VERSION-1      *
C     *                                                                *
C     *   AUTHOR....K.GYODA             1976.05                        *
C     *                                                                *
C     *   'DEXPI' COMPUTES THE DEXPONENTIAL INTEGRAL.                  *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL DEXPI(X,EI,ICON)                                   *
C     *          X   ....GIVEN ARGUMENT OF DEXPONENTIAL INTEGRAL.      *
C     *          EI  ....RESULTANT OF DEXPONENTIAL INTEGRAL.           *
C     *          ICON....RESULTANT CONDITION CODE.                     *
C     *                        0...NORMAL END.                         *
C     *                    20000...ARGUMENT VALUE IS UNREASONABLE (X). *
C     *                    30000...PARAMETER ERROR (X).                *
C     *   METHOD                                                       *
C     *        CHEBYSHEV APPROXIMATION.                                *
C     *   SLAVE SUBROUTINE                                             *
C     *        NONE.                                                   *
C     *   REMARKS                                                      *
C     *        THE RANGE OF X MUST BE AS FOLLOWS;                      *
C     *        X.GE.-174.673.AND.X.LE.174.673.AND.X.NE.0               *
C     ******************************************************************
C
      SUBROUTINE DEXPI(XX,EI,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 VV0,VV1,VV2,VV3,VV4,VV5,WW0,WW1,WW2,WW3,WW4,XX
      REAL*8 DMACH
!      character mcode(6)*2                       !  DIMENSION MCODE(6)
!      DATA MCODE/2HI1,2H1-,2H31,2H-0,2H10,2H1 /
      character(2) :: mcode(6)
      data mcode /"I1","1-","31","-0","10","1"/
      EXTERNAL DMACH
C    ------------------------------------------------------------------
C         P,Q,R,S,T,U TAKEN FROM MATHEMATICS OF COMPUTATION
C         VOL22,JULY 1968 P641-649
C    ------------------------------------------------------------------
      DATA P0 /-9.9999 9540  D-01/,
     1     P1 /-8.1325 5727  D 00/,
     2     P2 /-1.2451 4812  D 01/,
     3     P3 /-6.4250 5778  D-01/
      DATA Q1 /+1.0132 4792  D 01/,
     1     Q2 /+2.6721 0710  D 01/,
     2     Q3 /+1.7146 3138  D 01/
      DATA R0 /+1.7750 5462  D-05/,
     1     R1 /+9.9948 5008  D-01/,
     2     R2 /+4.2022 4422  D 00/,
     3     R3 /+3.6644 3110  D 00/,
     4     R4 /+5.3937 0234  D-01/
      DATA S1 /+5.1948 3844  D 00/,
     1     S2 /+6.9325 1459  D 00/,
     2     S3 /+2.4923 0348  D 00/,
     3     S4 /+1.5227 2584  D-01/
      DATA T0 /-1.6433 0274  D 02/,
     1     T1 /+2.2125 5491  D 02/,
     2     T2 /+2.9173 0565  D 01/,
     3     T3 /+4.3231 8059  D 00/
      DATA U0 /+2.8469 4759  D 02/,
     1     U1 /+1.0990 5660  D 02/,
     2     U2 /+1.6560 3618  D 01/
      DATA V0 /-4.3409 3034  D 01/,
     1     V1 /+2.8101 6387  D 01/,
     2     V2 /-3.1073 2436  D 01/,
     3     V3 /+2.2932 1941  D 00/,
     4     V4 /-2.7185 9075  D 00/,
     5     V5 /-4.2358 8307  D-01/
      DATA W0 /-3.9412 5513  D 01/,
     1     W1 /+8.3963 0455  D 01/,
     2     W2 /-7.5849 1696  D 01/,
     3     W3 /+3.6356 5454  D 01/,
     4     W4 /-9.2520 0680  D 00/
      DATA VV0/+7.6269 41797 55560 01160  D 01/,
     1     VV1/+1.7013 41431 36980 00575  D 02/,
     2     VV2/+1.2986 93223 08004 01519  D 02/,
     3     VV3/+3.8935 48293 38035 59125  D 01/,
     4     VV4/+3.7302 55272 46073 18563  D 00/,
     5     VV5/+2.7873 59112 66279 36858  D-02/
      DATA WW0/+7.6269 41797 55560 01160  D 01/,
     1     WW1/+2.0826 88521 24758 00633  D 02/,
     2     WW2/+2.0858 06090 45197 68375  D 02/,
     3     WW3/+9.2870 19124 20386 74802  D 01/,
     4     WW4/+1.7451 81064 78266 41644  D 01/
C    ------------------------------------------------------------------
C         A,B,C,D,E,F,X0 TAKEN FROM MATHEMATICS OF COMPUTATION
C         VOL23 APRIL 1969 P289-303
C    ------------------------------------------------------------------
      DATA A0 /+1.0021 5312  D 00/,
     1     A1 /-5.4799 3904  D 00/,
     2     A2 /+3.3879 3318  D 00/,
     3     A3 /-1.3967 8752  D 01/,
     4     A4 /+2.7114 0236  D 00/
      DATA B0 /+8.6616 6401  D-01/,
     1     B1 /+3.9509 9333  D 01/,
     2     B2 /+1.6364 6404  D 01/,
     3     B3 /+9.4202 7521  D 01/
      DATA C0 /+9.9961 7704  D-01/,
     1     C1 /+2.3561 1092  D 00/,
     2     C2 /+4.8766 0509  D 01/,
     3     C3 /-1.2148 0216  D 00/,
     4     C4 /-1.5070 4426  D 01/
      DATA D0 /+1.0581 5220  D 00/,
     1     D1 /-2.4698 5543  D 02/,
     2     D2 /+1.1868 9737  D 01/,
     3     D3 /+7.9086 5047  D 01/
      DATA E0 /+1.0000 0095  D 00/,
     1     E1 /-3.0103 0272  D 00/,
     2     E2 /-6.7666 5604  D 00/
      DATA F0 /+1.9997 4970  D 00/,
     1     F1 /-2.6592 3656  D 00/
      DATA X0 /+3.7250 7410  D-01/
C    ------------------------------------------------------------------
C         RATIONAL CHEBYSHEV APPROXIMATIONS FOR
C         THE DEXPONENTIAL INTEGRAL
C    ------------------------------------------------------------------
      X=DBLE(XX)
C
      ICON=0
      IF(X.GE.0.0E00) GO TO 30
      Y=-X
      IF(X.LT.-1.0E00) GO TO 10
C
C    ---- X .GE. -1 .AND. X .LT. 0 ------------------------------------
C
      EI=DLOG(Y)-(((T3*Y+T2)*Y+T1)*Y+T0)/(((Y+U2)*Y+U1)*Y+U0)
      RETURN
   10 IF(X.LT.-4.0E00) GO TO 20
C
C    ---- X .GE. -4 .AND. X .LT. -1 -----------------------------------
C
      EI=-DEXP(X)*((((R0*Y+R1)*Y+R2)*Y+R3)*Y+R4)/((((Y+S1)*Y+S2)*Y+S3)
     1   *Y+S4)
      RETURN
   20 IF(X.LT.-174.673E00) GO TO 9000
C
C    ---- X .GE. -174.673 .AND. X .LT. -4 -----------------------------
C
      EI=DEXP(X)*(1.0E00+(((P0*Y+P1)*Y+P2)*Y+P3)/((((Y+Q1)*Y+Q2)
     1   *Y+Q3)*Y))/X
      RETURN
   30 IF(X.EQ.0.0E00) GO TO 9001
C    ------------------------------------------------------------------
C         CHEBYSHEV APPROXIMATIONS FOR THE DEXPONENTIAL INTEGRAL
C    ------------------------------------------------------------------
      IF(X.GT.6.0E00) GO TO 50
C
C    ---- X .GT. 0 .AND. X .LE. 6 -------------------------------------
C
      W=X/6.0E00
      R=(((((V5*W+V4)*W+V3)*W+V2)*W+V1)*W+V0)
     1  /(((((W+W4)*W+W3)*W+W2)*W+W1)*W+W0)
      XX0=X-0.37250741078136663446 D00
      IF(ABS(XX0).GE.0.037E00) GO TO 40
      W=XX0/X0
      EI=((((((VV5*W+VV4)*W+VV3)*W+VV2)*W+VV1)*W+VV0)
     1   /((((((W+WW4)*W+WW3)*W+WW2)*W+WW1)*W+WW0)*X0)+R)*XX0
      RETURN
   40 EI=DLOG(X/X0)+XX0*R
      RETURN
   50 IF(X.GT.12.0E00) GO TO 60
C
C    ---- X .GT. 6 .AND. X .LE. 12 ------------------------------------
C
      EI=DEXP(X)*(A0+B0/(A1+X+(B1/(A2+X+(B2/(A3+X+(B3/(A4+X))))))))/X
      RETURN
   60 IF(X.GT.24.0E00) GO TO 70
C
C    ---- X .GT. 12 .AND. X .LE. 24 -----------------------------------
C
      EI=DEXP(X)*(C0+D0/(C1+X+(D1/(C2+X+(D2/(C3+X+(D3/(C4+X))))))))/X
      RETURN
   70 IF(X.GT.174.673E00) GO TO 9002
C
C    ---- X .GT. 24 .AND. X .LE. 174.673 ------------------------------
C
      EI=DEXP(X)/X*(1.0E00+(E0+F0/(E1+X+(F1/(E2+X))))/X)
      RETURN
 9000 ICON=20000
      EI=0.0E00
      GO TO 8000
 9001 ICON=30000
      EI=0.0E00
      GO TO 8000
 9002 ICON=20000
      WRITE(6,*) ICON
      call wexit("hcrsub","X.GT.174.673E00")
C     EI=7.2E75
 8000 CALL MGSSL(ICON,MCODE)
      RETURN
      END
C
C*    *** AQC8   ***   *   *********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
C     *                                                                *
C     *   G23-11-0301    AQC8,DAQC8,QAQC8              VERSION-5       *
C     *                                                                *
C     *   AUTHOR....TAKEMITSU HASEGAWA,TATSUO TORII                    *
C     *                                 1978.6                         *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL AQC8(A,B,FUN,EPSA,EPSR,NMIN,NMAX,S,ERR,N,ICON)     *
C     *          A   ....GIVEN LOWER LIMIT OF INTEGRATION.             *
C     *          B   ....GIVEN UPPER LIMIT OF INTEGRATION.             *
C     *          FUN ....GIVEN INTEGRAND.                              *
C     *          EPSA....GIVEN REQUIRED TOLERANCE(ABSOLUTE ERROR).     *
C     *          EPSR....GIVEN REQUIRED TOLERANCE(RELATIVE ERROR).     *
C     *          NMIN....GIVEN MINIMUM NUMBER OF FUNCTION EVALUATIONS. *
C     *          NMAX....GIVEN MAXIMUM NUMBER OF FUNCTION EVALUATIONS. *
C     *          S   ....RESULTANT APPROXIMATION TO THE INTEGRAL.      *
C     *          ERR ....RESULTANT ESTIMATION OF ABSOLUTE ERROR BOUND. *
C     *          N   ....RESULTANT NUMBER OF FUNCTION EVALUATIONS.     *
C     *          ICON....RESULTANT CONDITION CODE.                     *
C     *                 0....NORMAL END.                               *
C     *             10000....THE REQUIRED ACCURACY OF S IS NOT         *
C     *                      ATTAINED. THE ERROR OF RESULTANT S IS AS  *
C     *                      SMALL AS THE ACCUMULATED ROUNDOFF ERROR.  *
C     *             20000....THE REQUIRED ACCURACY OF S HAS NOT YET    *
C     *                      ATTAINED WHEN 511 SAMPLE POINTS HAVE BEEN *
C     *                      EXHAUSTED.                                *
C     *             30000....PARAMETER ERROR.                          *
C     *                                                                *
C     *   METHOD                                                       *
C     *          AUTOMATIC QUADRATURE WITH OPEN FORMULA BY A           *
C     *          GENERALIZED CHEBYSHEV INTERPOLATION PROCEDURE         *
C     *          INCREASING 8 SAMPLE POINTS AT A TIME.                 *
C     *                                                                *
C     *   SLAVE SUBROUTINE                                             *
C     *          AMACH                                                 *
C     *                                                                *
C     *   REMARK                                                       *
C     *          THIS SUBROUTINE WORKS EFFECTIVELY FOR SUFFICIENTLY    *
C     *          SMOOTH FUNCTIONS ON THE CLOSED INTERVAL OR FOR        *
C     *          OSCILLATING FUNCTIONS.                                *
C     *                                                                *
C     *                                                                *
C     *                                                                *
C     ******************************************************************
C     ******************************************************************
C     ******************************************************************
      SUBROUTINE DAQC8(A,B,EPSA,EPSR,NMIN,NMAX,S,ERR,N,ICON)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 C(256),W(257),V(256),FF(8),Z(4),U(4)
!      character mcode(6)*2                       !  DIMENSION MCODE(6)
!      DATA MCODE/2HG2,2H3-,2H11,2H-0,2H30,2H1 /
      character(2) :: mcode(6)
      data mcode /"G2","3-","11","-0","30","1"/
      REAL*8 FUN
      EXTERNAL FUN
      DATA ISW,L,M/0,4,8/,HALF/0.5D0/
      DATA FAC1,FAC2/20.D0,100.D0/
      ICON=30000
      IF(EPSA.LT.0.0D0.OR.EPSR.LT.0.0D0.OR.NMIN.LT.0.OR.NMAX.LT.NMIN)
     *     GO TO 8000
      NNMAX=MAX0(NMAX,15)
      IF(NMAX.GT.511) NNMAX=511
      NNMAX=(NNMAX-7)/8*8+7
      NNMIN=MIN0(NMIN,511)
      IF(ISW.EQ.0) GO TO 320
   10 COE=(B-A)*0.5D0
      CON=(B+A)*0.5D0
      ICON=0
      IF(A.EQ.B) GO TO 310
      EPSN=DABS(EPSA/COE)*8.0D0
C     ******************************************************************
C     *   INITIALIZATION.                                              *
C     *   SET CONSTANT.                                                *
C     ******************************************************************
      LL=-4
      L=-1
      L1=0
      J0=0
      JJ=1
      TEM=0.0D0
      SIN2A=C(1)
      I2=2
      I2D=1
      IP=-1
C     ******************************************************************
C     *   CALCULATION OF THE 4 CHEBYSHEV COEFFICIENTS C1...C4 FOR      *
C     *   FIRST 7 SAMPLE POINTS.                                       *
C     ******************************************************************
      X=COE*C(2)
      FF(1)=FUN(X+CON)
      FF(2)=FUN(CON-X)
      F1=FF(1)+FF(2)
      XX=COE*C(1)
      FF(3)=FUN(CON+XX)
      FF(4)=FUN(CON-XX)
      F2=FF(3)+FF(4)
      XX=COE*C(3)
      FF(5)=FUN(CON+XX)
      FF(6)=FUN(CON-XX)
      FF(7)=FUN(CON)
      F3=FF(5)+FF(6)
      F11=(F1-F3)*C(1)
      F4=FF(7)+FF(7)
      F3=F1+F3-F11
      F22=F2+F4
      F2=F2-F4
      C1=F3+F22
      C2=F11+F2
      C3=F11-F2
      C4=F3-F22
      AINTEG=C1*W(1)+C2*W(2)+C3*W(3)+C4*W(4)
      AMAXF=DMAX1(DABS(FF(1)),DABS(FF(2)),DABS(FF(3)),DABS(FF(4)),DABS(F
     *F(5)),DABS(FF(6)),DABS(FF(7)))
      AMINF=DMIN1(DABS(FF(1)),DABS(FF(2)),DABS(FF(3)),DABS(FF(4)),DABS(F
     *F(5)),DABS(FF(6)),DABS(FF(7)))
      ERROR1=(DABS(C3)+DABS(C4))*DABS(W(5))
      LN=7
      EPS1=AMAXF*AZERO
      EPS=DMAX1(DABS(EPSR*AINTEG),EPSN,EPS1)
      SS1=AINTEG
      IS=0
      IF(AMAXF.GT.AMINF*FAC1) GO TO 30
      IF(LN.LT.NNMIN) GO TO 30
      IF(ERROR1.LE.EPS) GO TO 280
      GO TO 30
C     ******************************************************************
C     *   PREPARATION FOR THE NEXT STEP WITH THE NUMBER OF THE SAMPLE  *
C     *   POINTS GREATER THAN OR EQUAL TO 15                           *
C     ******************************************************************
   20 CL0=CL1
      IP=IP+1
      IF(IP.LT.I2) GO TO 30
      IP=0
      I2D=I2
      I2=I2+I2
   30 L=L1
      IF(LN.GE.NNMAX) GO TO 300
      L1=L1+1
      L2=L1+L1
      L4=L2+L2
      LL=LL+4
      SINA=-C(L2+1)
      COSA=C(L2)
      COS2A=C(L1)
      IF(L1.EQ.1) GO TO 60
      IF(J0-JJ) 40,50,50
   40 J0=JJ
      SIN2A=-C(L1+1)
      TEM=C(J0)
      GO TO 60
   50 JJ=J0+1
      SIN2A=C(L)
      TEM=-C(J0)
   60 CL1=TEM+TEM
C     ******************************************************************
C     *   EVALUATION OF THE FUNCTIONAL VALUES AT 8 SAMPLE POINTS.      *
C     ******************************************************************
      XX=COE*C(L4)
      FF(1)=FUN(XX+CON)
      FF(2)=FUN(CON-XX)
      F0=FF(1)+FF(2)
      XA=COE*C(L4+1)
      FF(3)=FUN(CON-XA)
      FF(4)=FUN(CON+XA)
      F2=FF(3)+FF(4)
      X=COE*C(L4+2)
      FF(5)=FUN(X+CON)
      FF(6)=FUN(CON-X)
      F1=FF(5)+FF(6)
      XB=C(L4+3)*COE
      FF(7)=FUN(CON-XB)
      FF(8)=FUN(CON+XB)
      LN=LN+8
      AMAXF=DMAX1(DABS(FF(1)),DABS(FF(2)),DABS(FF(3)),DABS(FF(4)),DABS(F
     *F(5)),DABS(FF(6)),DABS(FF(7)),DABS(FF(8)),AMAXF)
      F3=FF(7)+FF(8)
      F02=F0+F2
      F0=F0-F2
      F13=F1+F3
      F1=F1-F3
      F00=F13+F02
      FC=(F02-F13)*COS2A
      F01=F0*COSA-F1*SINA
      COS4=COS2A+COS2A
      F0=(F0*COSA+F1*SINA)*COS4
      SN=COS4*SIN2A
      SNSN=1.0D0/SN
      IF(L.NE.0) GO TO 80
C     ******************************************************************
C     *   CALCULATION OF THE 4 PSEUDO CHEBYSHEV COEFFICIENTS V(1)...   *
C     *   V(4)                                                         *
C     ******************************************************************
      V(1)=(F0-F01-C4)*SNSN
      V(2)=(FC+F01-F0-C3)*SNSN
      V(3)=(F01-FC-C2)*SNSN
      V(4)=(F00-F01-C1)*SNSN
      RR=V(1)*W(5)+V(2)*W(6)+V(3)*W(7)+V(4)*W(8)
      AINTEG=AINTEG+RR
      PQ=DABS(RR)
      IF(ERROR1.LT.AZERO) ERROR1=AZERO
      PPQ=PQ/ERROR1
      PQ2=PPQ
      IF(PQ.GT.ERROR1) IS=1
      ERROR1=(DABS(V(4))+DABS(V(3)))*DABS(W(9))
      ERRORN=ERROR1
      SS1=AINTEG
      AMINF=DMIN1(DABS(FF(1)),DABS(FF(2)),DABS(FF(3)),DABS(FF(4)),DABS(F
     *F(5)),DABS(FF(6)),DABS(FF(7)),DABS(FF(8)),AMINF)
      IF(AMAXF.GT.AMINF*FAC2) GO TO 20
      GO TO 240
C     ******************************************************************
C     *   PREPARATION FOR THE CALCULATION OF V(LL+1)...V(LL+4).        *
C     ******************************************************************
   80 Z(1)=((C1-F00)*TEM-C4-F01+F0)*SNSN
      Z(2)=(TEM*C2-C3-F0+FC+F01)*SNSN
      Z(3)=(TEM*C3-C2-FC+F01)*SNSN
      Z(4)=(TEM*C4-C1-F01+F00)*SNSN
      IF(I2.GT.2) GO TO 120
      IF(IP.EQ.1) GO TO 100
C     ******************************************************************
C     *   CALCULATION OF 4 PSEUDO CHEBYSHEV COEFFICIENTS V(5)...V(8)   *
C     ******************************************************************
      DO 90 I=1,4
   90 V(I+4)=Z(I)-V(I)*SN
      GO TO 220
C     ******************************************************************
C     *   CALCULATION OF 4 PSEUDO CHEBYSHEV COEFFICIENTS V(9)...V(12)  *
C     ******************************************************************
  100 CO=1.0D0/(CL1-CL0)
      DO 110 I=1,4
  110 V(I+8)=(Z(I)-V(I)*SN-V(I+4))*CO
      GO TO 220
C     ******************************************************************
C     *   START FOR THE CALCULATION OF V(LL+1)...V(LL+4)               *
C     ******************************************************************
  120 I24=I2*4-8
      IA=I2D-1
      CO=CL1-C(IA)-C(IA)
      DO 130 I=1,4
      II=I24+I
  130 U(I)=V(II)*CO+V(II-4)
      IT=I24-8
      IF(I2D.LT.3) GO TO 150
      DO 140 J=3,I2D
      IA=IA-1
      CO=CL1+C(IA)+C(IA)
      COO=CL1-C(IA)-C(IA)
      U(1)=(U(1)*CO+V(IT+1))*COO+V(IT-3)
      U(2)=(U(2)*CO+V(IT+2))*COO+V(IT-2)
      U(3)=(U(3)*CO+V(IT+3))*COO+V(IT-1)
      U(4)=(U(4)*CO+V(IT+4))*COO+V(IT)
  140 IT=IT-8
  150 DO 160 I=1,4
  160 U(I)=Z(I)-(U(I)*CL1+V(I))*SN
      IF(IP.EQ.0) GO TO 200
      I=I24+5
      J01=J0-1
      IF(J01.LT.I2D) GO TO 180
      DO 170 J=I2D,J01
      CO=0.5D0/(TEM-C(J))
      COO=0.5D0/(TEM+C(J))
      U(1)=((U(1)-V(I))*CO-V(I+4))*COO
      U(2)=((U(2)-V(I+1))*CO-V(I+5))*COO
      U(3)=((U(3)-V(I+2))*CO-V(I+6))*COO
      U(4)=((U(4)-V(I+3))*CO-V(I+7))*COO
  170 I=I+8
  180 IF(J0.EQ.JJ) GO TO 200
      IM1=I-1
      CO=1.0D0/(CL1-CL0)
      DO 190 J=1,4
      IJ=IM1+J
      LLJ=LL+J
  190 V(LLJ)=(U(J)-V(IJ))*CO
      GO TO 220
  200 DO 210 I=1,4
      LLI=LL+I
  210 V(LLI)=U(I)
  220 CONTINUE
C     ******************************************************************
C     *   V(LL+1)...V(LL+4) HAVE BEEN OBTAINED.                        *
C     ******************************************************************
      QQ=V(LL+1)*W(LL+5)+V(LL+2)*W(LL+6)+V(LL+3)*W(LL+7)+V(LL+4)*W(LL+8)
C     ******************************************************************
C     *   MULTIPLE OF 'AINTEG' GIVES THE ESTIMATE FOR THE INTEGRAL.    *
C     ******************************************************************
      AINTEG=AINTEG+QQ
      ERROR1=(DABS(V(LL+4))+DABS(V(LL+3)))*DABS(W(LL+9))
      IF(L1.EQ.3.OR.L1.EQ.7) GO TO 250
      IF(L1.EQ.15.OR.L1.EQ.31.OR.L1.EQ.63) GO TO 250
  240 IF(IS.EQ.1) ERROR1=ERROR1*PPQ
      GO TO 270
  250 PQ=DABS(SS1-AINTEG)
      SS1=AINTEG
      IF(ERRORN.LT.AZERO) ERRORN=AZERO
      ERRRN=ERRORN*PPQ
      PQ1=PQ2
      PQ2=PQ/ERRORN
      PPQ=PQ2
      IF(PQ2.GT.PQ1) PPQ=PPQ/PQ1*PQ2
      IF(IS.EQ.1.AND.ERRRN.LT.PQ) PPQ=PPQ*PQ/ERRRN
      IS=0
      ERRORN=ERROR1
      IF(PPQ.LE.1.0D0) GO TO 270
      IS=1
      GO TO 240
  270 CONTINUE
      IF(LN.LT.NNMIN) GO TO 20
      EPS1=DFLOAT(L1)*AMAXF*AZERO
      EPS=DMAX1(DABS(EPSR*AINTEG),EPSN,EPS1)
C     ******************************************************************
C     *   CONVERGENCE TEST.                                            *
C     ******************************************************************
      IF(LN.EQ.15.AND.ERROR1.LE.EPS) GO TO 280
      IF(ERROR1.GT.EPS) GO TO 20
  280 CONTINUE
      IF(EPS.EQ.EPS1) ICON=10000
C     ******************************************************************
C     *   CONVERGENCE IS ATTAINED.                                     *
C     ******************************************************************
  290 S=AINTEG*COE*0.125D0
      IF(ERROR1.LT.EPS1) ERROR1=EPS1
      ERR=ERROR1*COE*0.125D0
      N=LN
      GO TO 8000
  300 ICON=20000
      GO TO 290
  310 ERR=0.0D0
      S=0.0D0
      N=0
      GO TO 8000
C     ******************************************************************
C     *   TABULATION OF TRIGONOMETRIC FUNCTION.                        *
C     ******************************************************************
  320 C(1)=DSQRT(HALF)
      C(2)=DSQRT((C(1)+1.0D0)*0.5D0)
      C(3)=-C(1)*0.5D0/C(2)
      C(4)=DSQRT((C(2)+1.0D0)*0.5D0)
      C(5)=C(3)*0.5D0/C(4)
      C(6)=(C(5)+C(4))*C(1)
      C(7)=(C(5)-C(4))*C(1)
      N0=2
      DO 340 K=4,M
      N1=N0+N0
      N2=N1+N1
      N4=N2+N2
      N3=N4+N2
      CA=DSQRT((C(N1)+1.0D0)*0.5D0)
      C(N2)=CA
      CP=0.5D0/CA
      SA=C(N1+1)*CP
      C(N2+1)=SA
      C(N4-1)=-(C(N1)+C(N0))*CP
      C(N4-2)=-(C(N1+1)+C(N0+1))*CP
      LP=N2+2
      MP=N2+N1-2
      DO 330 J=LP,MP,2
      J1=J-N1
      J2=N3-J
      P=C(J1)*CA
      Q=C(J1+1)*SA
      R=C(J1+1)*CA
      S=C(J1)*SA
      C(J)=P+Q
      C(J+1)=R-S
      C(J2-1)=Q-P
      C(J2-2)=-R-S
  330 CONTINUE
      N0=N1
  340 CONTINUE
C     ******************************************************************
C     *   TABULATION OF WEIGHT FUNCTION.                               *
C     ******************************************************************
      NN=L+L
      FNN=DFLOAT(NN)
      FFF=FNN+FNN
      FN2=FNN*FNN
      DO 350 I=1,L
      P=DFLOAT(I+I-1)
      W(I)=2.0D0/P
      IL=I+L
  350 W(IL)=FFF/(FN2-P*P)
      NNM1=NN-1
      IF(NNM1.LT.3) GO TO 410
      I2=NN
      I1=L
      DO 400 I=3,NNM1
      I2P=I2+I2
      P=DFLOAT(I2P+I2P)
      DO 360 K=1,I2
      KK=K+K-1
      IK=I2+K
      PP=DFLOAT(I2P-KK)*DFLOAT(I2P+KK)
  360 W(IK)=P/PP
      IP1=1
      IIPD=I1
      DO 390 IP=3,I
      LS=I2
      IP21=IP1+IP1-1
      DO 380 IS=IP1,IP21
      LS1=LS+IIPD
      CD=C(IS)+C(IS)
      LSP1=LS+1
      IILS=LS1+LSP1
      DO 370 K=LSP1,LS1
      LK=IIPD+K
      LK1=IILS-K
  370 W(LK)=W(LK)+W(LK1)-CD*W(K)
      LS=LS1+IIPD
  380 CONTINUE
      IP1=IP1+IP1
      IIPD=IIPD/2
  390 CONTINUE
      I1=I2
      I2=I2P
  400 CONTINUE
  410 ISW=1
      W(257)=1024.D0/(511.D0*513.D0)
      AZERO=DMACH(AZERO)*8.0D0
      GO TO 10
 8000 CALL MGSSL(ICON,MCODE)
      RETURN
      END
C
C ######################################################################
C
      REAL*8 FUNCTION FUN(X)
      REAL*8 PP,XPP,U,B,X
      COMMON PP,XPP
      U=X/XPP
      B=PP
      FUN=(1./(U+1.)+0.1728*(U-1.)/B**(2./3.)/(U+1.)**(5./3.)-0.0496*
     *(U**2+4./3.*U+1.)/B**(4./3.)/(U+1.)**(7./3.))*DEXP(-X)
      RETURN
      END
C
C*    *** MGSSL  ***   *   *********************************************
C     *                                                                *
C     *   SSL2 COPYRIGHT FUJITSU LTD. 1979                             *
C     *                                                                *
C     *   SLAVE          MGSSL                         VERSION-1       *
C     *                                                                *
C     *   AUTHOR....J.MIKAMI            1976.7                         *
C     *   CODER.....J.MIKAMI            1976.7                         *
C     *   IMPROVER..J.MIKAMI            1981.12        LEVEL-003       *
C     *                                                                *
C     *   'MGSSL' PRINTS THE CONDITION MESSAGE OF SSL2 SUBROUTINE .    *
C     *   THIS ROUTINE HAS SUBENTRIES NAMED 'MGSET' AND 'SSL2VL' .     *
C     *   'MGSET' SETS THE CONTROL INFORMATION FOR 'MGSSL' .           *
C     *   'SSL2VL' PRINTS THE CURRENT VERSION/LEVEL .                  *
C     *                                                                *
C     *   USAGE                                                        *
C     *        CALL MGSSL(ICON,MCODE)                                  *
C     *          ICON....GIVEN CONDITION CODE .                        *
C     *          MCODE...GIVEN IDENTIFICATION CODE OF SSL2 SUBROUTINE .*
C     *        CALL MGSET(ISET,IFLE)                                   *
C     *          ISET....GIVEN CONTROL INFORMATION OF PRINT LEVEL .    *
C     *          IFLE....GIVEN IDENTIFICATION NUMBER OF DATA SET FOR   *
C     *                    PRINT OUT .                                 *
C     *        CALL SSL2VL                                             *
C     *          (THIS ENTRY HAS NO PARAMETER)                         *
C     *   REMARK                                                       *
C     *        NORMALLY 'MGSSL' PRINTS NOTHING .                       *
C     *                                                                *
C     ******************************************************************
      SUBROUTINE    MGSSL(ICON,MCODE)
!      IMPLICIT REAL*8 (A-H,O-Z)
      implicit none
!      DIMENSION     MCODE(6)
      integer, intent(in) :: icon
      character(2), intent(out) :: mcode(6)
      integer, parameter :: iflg = -1, ifil = 6
!      DATA          IFLG/-1/,IFIL/6/
      IF(IFLG.LT.0)    RETURN
      IF(ICON.LT.IFLG) RETURN
      WRITE(IFIL,300)  MCODE,ICON
      RETURN
CC    ENTRY         MGSET(ISET,IFLE)
CC    IFLG  = ISET * 10000
CC    IF(ISET.LT.0)    RETURN
CC    IFIL  = IFLE
CC    RETURN
CC    ENTRY         SSL2VL
CC    WRITE(IFIL,200)
CC    RETURN
  100 FORMAT(1H0,   42H**** SSL2 COPYRIGHT FUJITSU LTD. 1979 ****)
  200 FORMAT(1H0/1X,
     *  50H**** LIBRARY STATUS (OSIV/F4 MSP SSL2 V10L20) ****/)
  300 FORMAT(1H ,   12H**** SSL2 ( ,6A2,12H) CONDITION ,I5,5H ****)
      END
C
C
C
C*    *** AMACH  ***   *   *********************************************
C     *                                                                *
C     *   SSL2,SSL2/VP COPYRIGHT FUJITSU LIMITED 1979,1983             *
C     *                                                                *
C     *   SLAVE          AMACH                         VERSION-1       *
C     *                                                                *
C     *   ---------- BASE = 16 , DIGIT =  6 ----------                 *
C     *                                                                *
C     *   AUTHOR....J.MIKAMI            1976.7                         *
C     *                                                                *
C     *   'AMACH' SETS THE UNIT ROUND OFF FOR SINGLE PRECISION         *
C     *           ARITHMETIC.                                          *
C     *                                                                *
C     *   USAGE(FUNCTION SUBPROGRAM)                                   *
C     *        AMACH(EPS)                                              *
C     *             EPS...GIVEN, DUMMY ARGUMENT.                       *
C     *                                                                *
C     *   REMARK                                                       *
C     *        AMACH IS MACHIN DEPENDENT ROUTINE.                      *
C     *                                                                *
C     ******************************************************************
      REAL*8 FUNCTION      DMACH(EPS)
      REAL*8 EPS
      DATA          O/1.0D-15/
      DMACH   = O
      RETURN
      END
