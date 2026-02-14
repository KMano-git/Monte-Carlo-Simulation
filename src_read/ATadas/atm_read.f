!***********************************************************************
      subroutine atm_read( dsninc, iunt12, iz0,   iclass, lpart, ifail,
     >  isdimd, izdimd, itdimd, iddimd,
     >  ismaxd, izmaxd, itmaxd, idmaxd,
     >  zdata,  dtevd,  ddensd,  xdaty,   ndty )
!***********************************************************************
!
!   copy titan /data/home2/shimizu/adas405/adas4xx/adaslib/dxrdnm
!   /home/shimizu/OLDWS/adas405/soc/adas4xx/adaslib/dxrdnm.for
!     add dimension size of denisty  iddimd
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  ****************** FORTRAN77 SUBROUTINE: DXRDNM *********************
!
! PURPOSE : TO EXTRACT COLLISIONAL DIELECTRONIC DATA  FROM
!           EITHER PARTIAL (METASTABLE/PARENT RESOLVED) OR STANDARD
!           (UNRESOLVED) ISONUCLEAR MASTER FILES
!
! NOTE    : THE SOURCE DATA IS CONTAINED AS SEQUENTIAL DATASETS
!           WITH THE FOLLOWING NAMING CONVENTIONS:
!
!                   (1) JETSHP.ACD<YR>#<EL).<CODE>DATA
!                   (2) JETSHP.SCD<YR>#<EL>.<CODE>DATA
!                   (3) JETSHP.CCD<YR>#<EL>.<CODE>DATA
!                   (4) JETSHP.PRB<YR>#<EL>.<FILT>.<CODE>DATA
!                   (5) JETSHP.PRC<YR>#<EL>.<FILT>.<CODE>DATA
!                   (6) JETSHP.QCD<YR>#<EL>.<CODE>DATA
!                   (7) JETSHP.XCD<YR>#<EL>.<CODE>DATA
!                   (8) JETSHP.PLT<YR>#<EL>.<CODE>DATA
!                   (9) JETSHP.PLS<YR>#<EL>.<CODE>DATA
!
!       WHERE, <YR>   = TWO DIGIT YEAR NUMBER
!              <EL>   = ONE OR TWO CHARACTER ELEMENT SYMBOL
!              <CODE> = R       => PARTIAL DATA
!                       U       => PARTIAL DATA
!                       OMITTED => STANDARD DATA
!              <FILT> = SIX CHARACTER POWER FILTER CODE
!
!       AND DATA OF CLASSES 6 AND 7 DO NOT EXIST FOR THE PARTIAL CASE.
!
!
! INPUT  : (C*120) DSNINC   = ISONUCLEAR MASTER FILE NAME - VERIFIED
!                           AND READY FOR DYNAMIC ALLOCATION.
! INPUT  : (L*4)  LPART     = .TRUE.  => PARTIAL (RESOLVED) MASTER DATA
!                            . FALSE. => UNSRESOLVED MASTER DATA
! INPUT  : (I*4)  IZ0       = NUCLEAR CHARGE
! INPUT  : (I*4)  NPART()   = METASTABLE PARTITION.  I.E. NUMBER OF
!                             METASTABLES FROM CHARGE STATE IZ1MIN-1 TO
!                            IZ1MAX ON INPUT
! INPUT  : (I*4)  IPRTD     = REQUIRED PARENT INDEX
! INPUT  : (I*4)  IGRDD     = REQUIRED GROUND INDEX
! INPUT  : (I*4)  ICLASS    = CLASS OF DATA (1 - 9 )
! INPUT  : (I*4)  IZ1       = REQUIRED ION CHARGE + 1
! INPUT  : (I*4)  ITMAX     = NUMBER OF ( DTEV() , DDENS() ) PAIRS
! INPUT  : (I*4)  ISDIMD    = MAXIMUM NUMBER OF (CHARGE, PARENT, GROUND)
!                             BLOCKS IN ISONUCLEAR MASTER FILES
! INPUT  : (I*4)  IZDIMD    = MAXIMUM NUMBER OF CHARGE STATES
!                             IN ISONUCLEAR MASTER FILES
! INPUT  : (I*4)  ITDIMD    = MAXIMUM NUMBER OF TEMP OR DENS VALUES IN
!                             ISOELECTRONIC MASTER FILES
! INPUT  : (I*4)  IDDIMD    = Maximum Number of DENS values ! KS100315
!                             Isoelectronic master files    ! KS100315
! INPUT  : (R*8)  DTEV()    = DLOG10(ELECTRON TEMPERATURES (EV))
! INPUT  : (R*8)  DDENS()   = DLOG10(ELECTRON DENSITIES (CM-3))
!
! OUTPUT : (I*4)  IFAIL     = 0    IF ROUTINE SUCCESSFUL - DATA FOR THE
!                                  REQUESTED YEAR USED.
!                           = 1    IF ROUTINE OPEN STATEMENT FAILED
!                           = 2    IF FILE EXISTS BUT REQUIRED DATA
!                                  BLOCK DOES NOT
! OUTPUT : (I*4)  ISMAXD    = NUMBER OF (CHARGE, PARENT, METASTABLE)
!                             BLOCKS IN SELECTED MASTER FILE
! OUTPUT : (I*4)  IZMAXD    = NUMBER OF ZDATA() VALUES IN SELECTED
!                             MASTER FILE
! OUTPUT : (I*4)  ITMAXD    = NUMBER OF DTEVD() VALUES IN SELECTED
!                             MASTER FILE
! OUTPUT : (I*4)  IDMAXD    = NUMBER OF DDENSD() VALUES IN SELECTED
!                             MASTER FILE
! OUTPUT : (I*4)  NPARTR()  = METASTABLE PARTITION.  I.E. NUMBER OF
!                             METASTABLES FROM CHARGE STATE IZ1MIN-1 TO
!                             IZ1MAX FOUND IN MASTER FILE
! OUTPUT : (R*8)  DTEVD()   = DLOG10(DATA ELECTRON TEMPERATURES (EV))
!                             IN SELECTED MASTER FILE
! OUTPUT : (R*8)  DDENSD()  = DLOG10(DATA ELECTRON DENSITIES (CM-3))
!                             IN SELECTED MASTER FILE
!-----------------------------------------------------------------------
! OUTPUT : (R*8)  DRCOFD(,,)= DLOG10(DATA RATE COEFFICIENTS (CM-3/S))
!                             IN SELECTED MASTER FILE
!                             1ST DIM: (CHARGE,META,GRD) BLOCK INDEX
!                             2ND DIM: TEMPERATURE INDEX
!                             3RD DIM: DENSITY INDEX
!-----------------------------------------------------------------------
! OUTPUT : (r*8)  xdaty()   = dlog10(data rate coefficients (cm-3/s))
!                 xdaty(jj) = drcofd(is,it,id)
!                    jj = id + idmaxd*(it-1) + idmaxd*itmaxd*(is-1)
! INPUT  : (i*4)  ndty      = maximum number of xdaty
!-----------------------------------------------------------------------
! OUTPUT : (R*8)  ZDATA()   = CHARGE + 1 FOR IONS IN SELECTED MASTER
!                             FILE
!                             1ST DIM: (CHARGE,META,GRD) BLOCK INDEX
! OUTPUT : (R*8)  DRCOFI()  = INTERPOLATION OF DRCOFD(,,) FOR
!                             DTEV() & DDENS()
!
! PROGRAM: (C*80) DSNOLD    = FILE NAME USED IN PREVIOUS CALL
!          (C*80) CLINE     = GENERAL CHARACTER VARIABLE
!          (C*80) CTERM     = TERMINATOR LINE - '-' FILLED VARIABLE
!          (C*4)) CPATRN()  = PATTERN USED TO DETECT DATA CLASS
!          (I*4)  IZ0D      = NUCLEAR CHARGE READ FROM MASTER FILE
!          (I*4)  IZ1MIN    = MINIMUM CHARGE+1 READ FROM MASTER FILE
!          (I*4)  IZ1MAX    = MAXIMUM CHARGE+1 READ FROM MASTER FILE
!          (I*4)  IABT      = ABORT CODE
!          (I*4)  INDSEL    = LOCATION OF (CHARGE,PRNT,GRND)
!                             DATA BLOCK IN FILE
!          (I*4)  IZDAT     = CURRENT DATA BLOCK ION CHARGE +1
!          (I*4)  ISEL      = GENERAL INDEX
!          (I*4)  I         = GENERAL INDEX
!          (I*4)  IT        = GENERAL INDEX
!          (I*4)  ID        = GENERAL INDEX
!          (I*4)  IZCHK     = INDEX TO VERIFY DATA Z1 SET COMPLETE
!          (I*4)  IPRTR()   = PARENT INDICES IN DATA SET
!          (I*4)  IGRDR()   = GROUND INDICES IN DATA SET
!          (I*4)  LCK       = MUST BE GREATER THAN 'ITMAXD' & 'IDMAXD'
!                             & 'ITMAX' - ARRAY SIZE FOR SPLINE CALCS.
!          (R*8)  A()       = GENERAL ARRAY
!          (R*8)  DRCOF0(,) = INTERPOLATION OF DRCOFD(,,) W.R.T DTEV()
!          (L*8)  LEXIST    = TRUE --- FILE TO OPEN EXISTS ELSE NOT
!          (I*4)  L1      = PARAMETER = 1
!          (I*4)  IOPT    = DEFINES THE BOUNDARY DERIVATIVES FOR THE
!                             SPLINE ROUTINE 'XXSPLN', SEE 'XXSPLN'.
!          (L*4)  LSETX   = .TRUE.  => SET UP SPLINE PARAMETERS RELATING
!                                      TO X-AXIS.
!                           .FALSE. => DO NOT SET UP SPLINE PARAMETERS
!                                      RELATING TO X-AXIS.
!                                      (I.E. THEY WERE SET IN A PREVIOUS
!                                            CALL )
!                           (VALUE SET TO .FALSE. BY 'XXSPLN')
!          (R*8)  DY()    = SPLINE INTERPOLATED DERIVATIVES
!
!
! ROUTINES:
!          ROUTINE    SOURCE    BRIEF DESCRIPTION
!          ------------------------------------------------------------
!          I4UNIT     ADAS      FETCH UNIT NUMBER FOR OUTPUT OF MESSAGES
!          I4FCTN     ADAS      CONVERT STRING TO INTEGER FORM
!
!          (R*8 ADAS FUNCTION - 'R8FUN1' ( X -> X) )
!
! AUTHOR : H. P. SUMMERS, JET
!          K1/1/57
!          JET  EXT. 4941
!
! DATE   : 24/04/94
!
! UPDATE : 21/07/94 - HPS - BYPASS CHECK ON CHARGE STATE COMPLETENESS
!                           FOR XCD AND QCD FILES
!
! UNIX-IDL PORT:
!
! VERSION: 1.1                          DATE: 08-11-95
! MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
!               - FIRST RELEASE
!
! VERSION: 1.2                          DATE: 22-11-95
! MODIFIED: TIM HAMMOND (TESSELLA SUPPORT SERVICES PLC)
!               - CHANGED TEST FOR ION LIMITS SLIGHTLY FROM
!                 IZ1MIN.GE.IZ1MAX TO IZ1MIN.GT.IZ1MAX TO ALLOW
!                 RUNS FOR HYDROGEN TO PROCEED.
!
! VERSION: 1.3				DATE: 13-10-99
! MODIFIED: Martin O'Mullane
!		- PRB definition has been changed and they are now
!                 summed over the parents. This necessitates accessing
!                 the data more like PLT/PRC/PLS than the others.
!               - DSNOLD made same size as DSNINC
!
!
!-----------------------------------------------------------------------
      character, intent(in)  :: dsninc*120
      integer,   intent(in)  :: iclass, iddimd, isdimd, itdimd, iunt12
     >                        , iz0, izdimd, ndty
      logical,   intent(in)  :: lpart
      integer,   intent(out) :: ifail, idmaxd, ismaxd, izmaxd, itmaxd
      real(8),   intent(out) :: ddensd(iddimd), dtevd(itdimd)
     >                        , xdaty(ndty), zdata(isdimd)
       integer, parameter :: LCK = 100, L1 = 1
!-----------------------------------------------------------------------
       INTEGER   IZ0D    , IZ1MIN  , IZ1MAX  , IABT    , IZDAT
       INTEGER   I       , IT      , ID
       INTEGER   IZCHK
       integer   jst, jj
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
       INTEGER   NPART(IZDIMD)     , NPARTR(IZDIMD)
       INTEGER   IGRDR(LCK) , IPRTR(LCK)
       CHARACTER CLINE*80  , CTERM*80   ,  CHINDI*4
!-----------------------------------------------------------------------
       CHARACTER CPATRN(9)*4
       LOGICAL   LEXIST
!-----------------------------------------------------------------------
!xx       EXTERNAL  R8FUN1
!-----------------------------------------------------------------------
       DATA      CTERM /'-----------------------------------------------
     &---------------------------------'/
       DATA      CPATRN/'IPRT','IPRT','IPRT','IPRT','IPRT',
     &                  'IGRD','IPRT','IGRD','IGRD'/
!-----------------------------------------------------------------------
       SAVE      IZ1MIN
!-----------------------------------------------------------------------
      integer  n6
! function
      integer   i4unit, i4fctn, lenx
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! ON FIRST ENTRY: MAKE SURE SET ARRAY BOUNDS ARE VALID
!-----------------------------------------------------------------------
!
            IF (LCK.LT.ITDIMD) call wexit("atom_read",
     &                   ' DXRDNM ERROR: ARRAY DIMENSION LCK < ITDIMD')
            IF (LCK.LT.IDDIMD) call wexit("atom_read",      ! KS100315
     &       ' DXRDNM ERROR: ARRAY DIMENSION LCK < IDDIMD') ! KS100315
            IF (LCK.LT.ISDIMD) call wexit("atom_read",
     &                   ' DXRDNM ERROR: ARRAY DIMENSION LCK < ISDIMD')
!
      n6 = 10
      write(n6,'(/2x,"*** atm_read ***")')
      write(n6,'(2x,"itdimd =",i3,"  iddimd =",i3,"  isdimd =",
     >  i3,"  <","  lck =",i3)') itdimd,iddimd,isdimd,lck  ! KS100315
!
!-----------------------------------------------------------------------
!      BYPASS MASTER FILE READ IF PREVIOUSLY DONE
!-----------------------------------------------------------------------
!
       IFAIL  = 0
!
!-----------------------------------------------------------------------
!      CONFIRMATORY CHECK THAT ISONUCLEAR MASTER FILE FILE
!-----------------------------------------------------------------------
!
       INQUIRE(FILE=DSNINC,EXIST=LEXIST)
!
       IF( .NOT.LEXIST ) GOTO 9999
       write(n6,'(2x,"dsn =",a)') dsninc(1:lenx(dsninc))
!
!-----------------------------------------------------------------------
!      READ FILE.  REVERIFY THAT TYPE IS CORRECLY IDENTIFIED
!-----------------------------------------------------------------------
!
!
       READ(IUNT12,'(A80)')CLINE
       IZ0D = I4FCTN(CLINE(1:5),IABT)
       IF(IABT.GT.0.OR.IZ0D.NE.IZ0) THEN
            WRITE(I4UNIT(-1),2001)'INCORRECT NUCLEAR CHARGE'
            WRITE(I4UNIT(-1),2002)
            call wexit("atom_read","incorrect nuclear charge")
       ENDIF
       IDMAXD = I4FCTN(CLINE(6:10),IABT)
       IF(IABT.GT.0.OR.IDMAXD.LE.0.OR.IDMAXD.GT.IDDIMD) THEN ! KS100315
              WRITE(I4UNIT(-1),2001)'INVALID NUMBER OF DENSITIES'
              WRITE(I4UNIT(-1),2002)
            call wexit("atom_read","invalid number of densities")
       ENDIF
       ITMAXD = I4FCTN(CLINE(11:15),IABT)
       IF(IABT.GT.0.OR.ITMAXD.LE.0.OR.ITMAXD.GT.ITDIMD) THEN
            WRITE(I4UNIT(-1),2001)'INVALID NUMBER OF TEMPERATURES'
            WRITE(I4UNIT(-1),2002)
            call wexit("atom_read","invalid number of temperatures")
       ENDIF
       IZ1MIN = I4FCTN(CLINE(16:20),IABT)
       IZ1MAX = I4FCTN(CLINE(21:25),IABT)
       IF(IABT.GT.0.OR.IZ1MIN.GT.IZ1MAX.OR.IZ1MIN.LT.1
     &     .OR.IZ1MAX.GT.IZ0) THEN
           WRITE(I4UNIT(-1),2001)'INCORRECT ION LIMITS'
           WRITE(I4UNIT(-1),2002)
           call wexit("atom_read","incorrect ion limits")
       ENDIF
!
       READ(IUNT12,'(A80)') CLINE
       IF(LPART) THEN
           READ(IUNT12,'(A80)') CLINE
           IF (INDEX(CLINE,'.').GT.0)THEN
              WRITE(I4UNIT(-1),2001)'INCORRECT STRUCTURE'
              WRITE(I4UNIT(-1),2002)
              call wexit("atom_read","incorrect structure")
           ELSE
               READ(CLINE,'(16I5)') (NPARTR(I),I=1,
     &                               MIN0(IZ1MAX-IZ1MIN+2,16))
           ENDIF
           IF(IZ1MAX+IZ1MAX-2.GT.16) THEN
               READ(IUNT12,'(16I5)') (NPARTR(I),I=17,IZ1MAX-IZ1MIN+2)
           ENDIF
!::debug
       write(6,'(2x,"iz1min, iz1max, imax =",3i5)')
     >    iz1min, iz1max, iz1max-iz1min+2
       write(6,'(2x,"npartr =",10i5)') (NPARTR(I),I=1,IZ1MAX-IZ1MIN+2)
!------
           DO 5 I=1,IZ1MAX-IZ1MIN+2
            npart(i) = npartr(i)       !  <== add 1-line  2003/01/07
            IF(NPART(I).NE.NPARTR(I))THEN
                WRITE(I4UNIT(-1),2001)'PARTITION MISMATCH'
                WRITE(I4UNIT(-1),2002)
                call wexit("atom_read","partition mismatch")
            ENDIF
    5      CONTINUE
           READ(IUNT12,'(A80)') CLINE
           IF(CLINE(1:5).NE.'-----') THEN
              WRITE(I4UNIT(-1),2001)'INCORRECT STRUCTURE'
              WRITE(I4UNIT(-1),2002)
              call wexit("atom_read","incorrect structure")
           ENDIF
       ENDIF
!
       READ(IUNT12,1040) ( DDENSD(ID) , ID = 1 , IDMAXD )
       READ(IUNT12,1040) ( DTEVD(IT) , IT = 1 , ITMAXD )
!
       CHINDI = CPATRN(ICLASS)
!
       jst = 0
       ISMAXD = 0
   10  READ(IUNT12,'(A80)',END=17)CLINE
       if( cline(1:2).eq."  " ) goto 10
       IF(CLINE(2:80).NE.CTERM(2:80)) THEN
           IF(LPART) THEN
               IF( CLINE(24:27) .EQ. CHINDI) THEN
                  ISMAXD = ISMAXD + 1
                  READ(CLINE,1070) IPRTR(ISMAXD), IGRDR(ISMAXD), IZDAT
                  ZDATA(ISMAXD) = DFLOAT( IZDAT )
               ELSE
                  WRITE(I4UNIT(-1),2001)'INCORRECT CLASS IPRT CODE',
     &                                  ICLASS,CHINDI
                  WRITE(I4UNIT(-1),2002)
                  call wexit("atom_read","incorrect class iprt code")
              ENDIF
           ELSE
              ISMAXD = ISMAXD + 1
              READ(CLINE,1071) IZDAT
              ZDATA(ISMAXD) = DFLOAT( IZDAT )
           ENDIF
!------
           do it = 1, itmaxd
             read(iunt12,1040) (xdaty(jj),jj=jst+1, jst+idmaxd)
             jst = jst + idmaxd
           enddo
!------
           GO TO 10
       ENDIF
!
   17  continue
       CLOSE(IUNT12)
!
!-----------------------------------------------------------------------
!      VERIFY Z1 SET IN MASTER FILE CONSISTENT WITH IZ1MIN AND IZ1MAX
!      EXCEPT FOR XCD AND QCD CASES
!-----------------------------------------------------------------------
!       xxxxxxx
!
       IF(ICLASS.EQ.6.OR.ICLASS.EQ.7) THEN
           IZMAXD = IZ1MAX-IZ1MIN+1
       ELSE
           IZCHK = IZ1MIN
           DO 18 I=1,ISMAXD
            IZDAT = INT(ZDATA(I))
            IF(IZDAT.NE.IZCHK) THEN
                IZCHK = IZCHK+1
            ENDIF
 18               CONTINUE
           IF(IZCHK.EQ.IZ1MAX) THEN
                IZMAXD = IZ1MAX-IZ1MIN+1
           ELSE
                WRITE(I4UNIT(-1),2001)'INCONSISTENT Z1 SET IN FILE'
                WRITE(I4UNIT(-1),2002)
                call wexit("atm_read","IZCHK.nq.IZ1MAX")
           ENDIF
       ENDIF
!
      write(n6,'(2x,"izmaxd =",i5)') izmaxd
      goto 500
!
!-----------------------------------------------------------------------
!      SELECT CORRECT DATA BLOCK REQUESTED
!-----------------------------------------------------------------------
!       xxxxxxxx
!
!-----------------------------------------------------------------------
!      INTERPOLATE USING SPLINES
!-----------------------------------------------------------------------
!       xxxxxx
 500   continue
       write(n6,'(2x,30("-"),"  END : atm_read")')
!
       RETURN
!
!-----------------------------------------------------------------------
! DATA SET OPENING/EXISTENCE ERROR HANDLING
!-----------------------------------------------------------------------
!
 9999  IFAIL  = 1
       RETURN
!
!-----------------------------------------------------------------------
!
 1000  FORMAT('FILE = ',1A30)
 1010  FORMAT(5I5)
 1040  FORMAT(8F10.5)
 1050  FORMAT(I6)
 1060  FORMAT(1X,'NOTE: REQUESTED DATASET - ',A30,' DOES NOT EXIST.'/
     1        7X,      'USING DEFAULT YEAR (',A2,') DATASET INSTEAD'/)
 1070  FORMAT(28X,I2,9X,I2,16X,I2)
 1071  FORMAT(57X,I2)
!
 2001 FORMAT(1X,32('*'),' DXRDNM ERROR ',32('*')//
     &       1X,'FAULT IN MASTER DATA FILE: ',A,I3,A)
 2002 FORMAT(/1X,28('*'),'  PROGRAM TERMINATED   ',27('*'))
!
 2010  FORMAT(1H ,'INDSEL=',I2,3X,'IZDAT =',I2,3X,
     &            'IPRTD =',I2,3X,'IGRDD =',I2,/
     &        /7F10.5/7F10.4/7F10.5)
!-----------------------------------------------------------------------
!
       END
