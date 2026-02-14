!X ULTRIX PORT - SCCS info: Module @(#)i4fctn.for	1.2 Date 04/24/95
!X
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   FUNCTION I4FCTN( STR , IABT )
      integer function i4fctn( str , iabt )
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
!  *************** FORTRAN77 INTEGER*4 FUNCTION: I4FCTN ****************
!
!  FUNCTION: TO CONVERT AN INTEGER NUMBER STORED IN THE STRING 'STR'
!            INTO A INTEGER*4 VARIABLE, USING INTERNAL READ.
!            INTIALLY THE PROGRAM CHECKS TO SEE IF THE NUMBER IS OF A
!            VALID FORM.
!
!  CALLING PROGRAM: GENERAL USE
!
!  FUNCTION:
!
!          (I*4)   I4FCTN  = FUNCTION NAME
!          (C*(*)) STR     = STRING CONTAINING SINGLE FLOATING POINT NO.
!          (I*4)   IABT    = RETURN CODE:
!                               0 => NO ERROR
!                               1 => ERROR (A VALUE 'I4FCTN=0' WILL BE
!                                           RETURNED).
!
!          (C*1)   CH0      = PARAMETER = '0'
!          (C*1)   CH9      = PARAMETER = '9'
!          (C*1)   BLANK    = PARAMETER = ' '
!          (C*1)   CPLUS    = PARAMETER = '+'
!          (C*1)   CMINUS   = PARAMETER = '-'
!
!          (I*4)   ILEN     = LENGTH OF 'STR' STRING IN BYTES
!          (I*4)   ILAST    = POSITION OF LAST BYTE OF IDENTIFIED NUMBER
!          (I*4)   I1       = STARTING BYTE IN 'STR' OF NUMBER
!                             INCLUDING SIGN IF PRESENT
!          (I*4)   IS       = 0 => NUMBER HAS NO SIGN
!                             1 => NUMBER HAS A SIGN
!          (I*4)   ICH0     = ICHAR('0')
!          (I*4)   ICH9     = ICHAR('9')
!          (I*4)   ISTR     = ICHAR(CURRENT BYTE POSITION IN 'STR')
!          (I*4)   I        = GENERAL USE
!
!          (L*4)   LFOUND   = .TRUE.  => ALL OF THE INPUT NUMBER BYTES
!                                        HAVE BEEN ASSESSED.
!                             .FALSE. => INPUT NUMBER BYTES STILL BEING
!                                        ASSESSED.
!          (L*4)   LSTART   = .TRUE.  => THE FIRST DIGIT HAS BEEN FOUND
!                             .FALSE. => THE FIRST DIGIT HAS NOT YET
!                                        BEEN REACHED.
!
!          (C*5)   CFORM5   = FORMAT FOR INTERNAL READING OF INTEGER
!
!
! NOTE:     AN ERROR WILL OCCUR (IABT=1) IF THERE IS MORE THAN ONE
!           NUMBER OCCURING IN THE STRING 'STR()'
!
!
! AUTHOR:   PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
!           K1/0/37
!           JET EXT. 2520
!
! DATE:     11/07/90
!
! UPDATE:   11/02/92 - PE BRIDEN: BLANKS NOW ALLOWED BETWEEN SIGN AND
!                                 FIRST DIGIT. LSTART VARIABLE ADDED.
!                                 VARIABLE I2 REMOVED.
!                                 + SOME MINOR RECODING - (IF STRING
!                                 ENTERED IS BLANK IABT IS NOW SET TO 1)
!
! UPDATE:   16/08/93 - PE BRIDEN: CORRECTED BUG TO ALLOW BLANKS BETWEEN
!                                 SIGN AND FIRST DIGIT (SEE ABOVE).
!                                 1) ILAST VARIABLE ADDED.
!                                 2) FORMATTED READ USED INSTEAD OF *
!                                    WHEN CONVERTING IDENTIFIED INTEGER
!                                    USING THE INTERNAL READ. (THIS
!                                    RESTRICTS IDENTIFIED NUMBER TO BE
!                                    < 100 BYTES IN LENGTH!)
!                                 3) EXCLUDE TRAILING BLANKS IN THE
!                                    INTERNAL READING OF THE INTEGER
!                                    I.E. STR(I1:ILAST) INSTEAD OF
!                                         STR(I1:ILEN)
!
! UPDATE:   07/03/95 - PE BRIDEN: INSTEAD OF USING FORMAT SPECIFIER I99
!                                 WHEN INTERNALLY READING THE INTEGER
!                                 CREATE THE APPROPRIATE SPECIFIER
!                                 WITHIN CFORM5 AND USE THIS.
!
!-----------------------------------------------------------------------
!
! added 2 lines organize local variables and include files by kamata 2021/06/28
      character, intent(in)  :: str*(*)
      integer,   intent(out) :: iabt
!-----------------------------------------------------------------------
      CHARACTER  CH0*1   , CH9*1   , BLANK*1   , CPLUS*1   , CMINUS*1
!-----------------------------------------------------------------------
      PARAMETER( CH0='0', CH9='9', BLANK=' ', CPLUS='+', CMINUS='-' )
!-----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   CHARACTER  STR*(*)
!-----------------------------------------------------------------------
! deleted 1 line organize local variables and include files by kamata 2021/06/28
!ik   INTEGER    I4FCTN  , IABT
      INTEGER    I1      , IS      , ILEN      , ILAST     ,
     &           ICH0    , ICH9    , ISTR      , I
!-----------------------------------------------------------------------
      LOGICAL    LSTART  , LFOUND
!-----------------------------------------------------------------------
      CHARACTER  CFORM5*5
!-----------------------------------------------------------------------
      DATA       CFORM5 / '(I??)' /
!-----------------------------------------------------------------------
      SAVE       CFORM5
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! INITIALIZE VALUES
!-----------------------------------------------------------------------
!
      I4FCTN = 0
      IABT   = 0
      I1     = 0
      IS     = 0
      LSTART = .FALSE.
      LFOUND = .FALSE.
      ICH0   = ICHAR(CH0)
      ICH9   = ICHAR(CH9)
      ILEN   = LEN(STR)
      ILAST  = ILEN
!
!-----------------------------------------------------------------------
! FIND STARTING BYTE OF NUMBER
!-----------------------------------------------------------------------
!
         DO 1 I=1,ILEN
               IF ( STR(I:I).NE.BLANK ) THEN
                  I1 = I
                  GOTO 2
               ENDIF
    1    CONTINUE
!
!-----------------------------------------------------------------------
! IDENTIFY IF NUMBER HAS A SIGN
!-----------------------------------------------------------------------
!
    2    IF (I1.EQ.0) THEN
            IABT = 1
            RETURN
         ENDIF
!
      IF   ( ( STR(I1:I1).EQ.CPLUS  )
     &                  .OR.
     &       ( STR(I1:I1).EQ.CMINUS ) ) IS=1
!
!-----------------------------------------------------------------------
! IDENTIFY IF NUMBER IS OF A VALID FORM
!-----------------------------------------------------------------------
!
         DO 3 I=I1+IS,ILEN

               IF (LFOUND) THEN
!
!-----------------------------------------------------------------------
! INPUT NO. COMPLETELY DEFINED: IDENTIFY IF EXTRA NON-BLANK BYTES EXIST
!-----------------------------------------------------------------------
!
                  IF (STR(I:I).NE.BLANK) IABT=1
!-----------------------------------------------------------------------
               ELSEIF (STR(I:I).EQ.BLANK) THEN
                  LFOUND = LSTART
!-----------------------------------------------------------------------
               ELSE
                  LSTART = .TRUE.
                  ILAST  = I
                  ISTR   = ICHAR(STR(I:I))
                  IF ( (ISTR.LT.ICH0) .OR. (ISTR.GT.ICH9) ) IABT=1
               ENDIF
!
!-----------------------------------------------------------------------
! RETURN ERROR CODE IF ERROR FOUND
!-----------------------------------------------------------------------
!
            IF (IABT.NE.0) RETURN
    3    CONTINUE
!
!-----------------------------------------------------------------------
! IDENTIFY IF VALID NUMBER FOUND (RECODED: PEB 11/02/92)
!                                (RECODED: PEB 07/03/95 - ADDED CFORM5)
! YES => USE INTERNAL READ TO OBTAIN THE INTEGER NUMBER
! NO  => RETURN ERROR CODE IF ERROR FOUND
!-----------------------------------------------------------------------
!
         IF (LSTART) THEN
            I      = 1 + ILAST - I1
            I      = MIN0(I,99)
            WRITE(CFORM5(3:4),'(I2.2)') I
            READ(STR(I1:ILAST),CFORM5) I4FCTN
         ELSE
            IABT=1
         ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
