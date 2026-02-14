!X UNIX PORT - SCCS info: Module @(#)i4unit.for	1.1 Date 04/10/95
!X
! modified 1/1 lines organize local variables and include files by kamata 2021/06/28
!ik   FUNCTION I4UNIT( IUNIT )
      integer function i4unit( iunit )
      IMPLICIT NONE
!-----------------------------------------------------------------------
!
!  ************** FORTRAN77 INTEGER*4 FUNCTION: I4UNIT *****************
!
!  PURPOSE: TO RESET OR RETURN A STORED INTEGER*4 VALUE GREATER THAN OR
!           EQUAL TO ZERO.
!           THIS IS USED WITHIN ADAS TO STORE THE STREAM/UNIT NUMBER
!           FOR THE OUTPUT OF ERROR MESSAGES (TO THE SCREEN).
!
!           BY DEFAULT THE STORED VALUE WILL BE 6, AND WILL BE RETURNED
!           BY THE FUNCTION IF IUNIT ON INPUT < 0.
!
!           TO RESET THE STORED VALUE THEN SET IUNIT TO THE REQUIRED
!           POSITIVE INTEGER (INC. ZERO). THIS VALUE WILL ALSO BE
!           RETURNED BY THE FUNCTION.
!
!                 IUNIT VALUE               RETURNED FUNCTION VALUE
!                 -----------               -----------------------
!                 IUNIT <  0            = CURRENT STORED INTEGER VALUE
!                                         (6 BY DEFAULT).
!                 IUNIT >= 0            = IUNIT , AND RESETS THE STORED
!                                                 VALUE TO IUNIT.
!
!
!  CALLING PROGRAM: GENERAL USE
!
!  SUBROUTINE:
!
!  O     : (I*4)  I4UNIT   = FUNCTION NAME - (SEE ABOVE)
!
!  I     : (I*4)  IUNIT    = FUNCTION ARGUMENT - (SEE ABOVE)
!
!          (I*4)  IDEFLT   = PARAMETER = DEFAULT STORED INTEGER VALUE
!
!          (I*4)  ICURNT   = CURRENT STORED INTEGER VALUE
!
!
! ROUTINES:
!          ROUTINE    SOURCE    BRIEF DESCRIPTION
!          ------------------------------------------------------------
!
!
! AUTHOR:  PAUL E. BRIDEN (TESSELLA SUPPORT SERVICES PLC)
!          K1/0/37
!          JET EXT. 5023
!
! DATE:    23/04/93
!
! UPDATE:  24/05/93 - PE BRIDEN - ALLOWED 0 TO BE A VALID STORED NUMBER
!
!-----------------------------------------------------------------------
!ik   INTEGER    I4UNIT     , IDEFLT
!-----------------------------------------------------------------------
!xx   PARAMETER( IDEFLT = 6 )
!      PARAMETER( IDEFLT = 31 )
! -----------------------------------------------------------------------
!      INTEGER    IUNIT      , ICURNT
!-----------------------------------------------------------------------
!      SAVE       ICURNT
      integer, intent(in) :: iunit
      integer, parameter :: ideflt = 31
      integer, save :: ICURNT
!-----------------------------------------------------------------------
      DATA       ICURNT / IDEFLT /
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
! RETRIEVE OR RESET STORED INTEGER VALUE ACCORDINGLY
!-----------------------------------------------------------------------
!
      IF (IUNIT.LT.0) THEN
        I4UNIT = ICURNT
      ELSE
        ICURNT = IUNIT
        I4UNIT = IUNIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      RETURN
      END
