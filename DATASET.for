C  ################### SUBROUTINE ##################################
C
C  ###########################################################
C  #                                                         #
C  #       PREPARE EARTH MODEL DATA  FOR INTEGRATION         #
C  #     AT 'I1234=?', THE ? STEP OF RUNGE-KUTTA METHOD      #
C  #                                                         #
C  ###########################################################
C
      SUBROUTINE DATASET(I1234,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THIS SUBROUTINE PREPARES THE DATA SET FROM THE GIVEN EARTH MODEL,
C  USED FOR RUNGE-KUTTA'S INTEGRATION.
C
C    ARGUMENTS: (ALL ARE INPUT)
C               I1234 -- ORDER NUMBER OF THE FOUR RUNGE-KUTTA'S
C                        COEFFICIENTS: Z = Z + (K1+2K2+2K3+K4)/6
C               MSTEP -- MAXIMUM LAYER OF ADOPTED EARTH MODEL
C               GRAVITY, RHO, LAMDA, MU, RADIUS -- EARTH MODEL
C    OUTPUT IN: COMMON /MDI/
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8    LAMDA,MU,LAMI,MUI
      DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
      COMMON    /INT1/ II
      COMMON    /INT2/ JJ
      COMMON    /INT3/ IJ

C        WRITE(*,'(2I6)') II,JJ

      IF (I1234.EQ.1) THEN

         RHOI = RHO(2*II-1)    + (RHO(2*II)    -RHO(2*II-1))
     &     /IJ*(JJ-1.0)
         LAMI = LAMDA(2*II-1)  + (LAMDA(2*II)  -LAMDA(2*II-1))
     &     /IJ*(JJ-1.0)
         MUI  = MU(2*II-1)     + (MU(2*II)     -MU(2*II-1))
     &     /IJ*(JJ-1.0)
         GRAVI= GRAVITY(2*II-1)+ (GRAVITY(2*II)-GRAVITY(2*II-1))
     &     /IJ*(JJ-1.0)
         RADI = RADIUS(2*II-1) + (RADIUS(2*II) -RADIUS(2*II-1))
     &     /IJ*(JJ-1.0)
         STEP =              (RADIUS(2*II)-RADIUS(2*II-1))
     &     /IJ

      ELSEIF(I1234.EQ.2) THEN

         RHOI = RHO(2*II-1)    + (RHO(2*II)    -RHO(2*II-1))
     &     /IJ*(JJ-0.5)
         LAMI = LAMDA(2*II-1)  + (LAMDA(2*II)  -LAMDA(2*II-1))
     &     /IJ*(JJ-0.5)
         MUI  = MU(2*II-1)     + (MU(2*II)     -MU(2*II-1))
     &     /IJ*(JJ-0.5)
         GRAVI= GRAVITY(2*II-1)+ (GRAVITY(2*II)-GRAVITY(2*II-1))
     &     /IJ*(JJ-0.5)
         RADI = RADIUS(2*II-1) + (RADIUS(2*II) -RADIUS(2*II-1))
     &     /IJ*(JJ-0.5)
         STEP =              (RADIUS(2*II)-RADIUS(2*II-1))
     &     /IJ

      ELSEIF(I1234.EQ.4) THEN

         RHOI = RHO(2*II-1)    + (RHO(2*II)    -RHO(2*II-1))
     &     /IJ*JJ
         LAMI = LAMDA(2*II-1)  + (LAMDA(2*II)  -LAMDA(2*II-1))
     &     /IJ*JJ
         MUI  = MU(2*II-1)     + (MU(2*II)     -MU(2*II-1))
     &     /IJ*JJ
         GRAVI= GRAVITY(2*II-1)+ (GRAVITY(2*II)-GRAVITY(2*II-1))
     &     /IJ*JJ
         RADI = RADIUS(2*II-1) + (RADIUS(2*II) -RADIUS(2*II-1))
     &     /IJ*JJ
         STEP =              (RADIUS(2*II)-RADIUS(2*II-1))
     &     /IJ

      ENDIF

       RADI=RADI/6371.0;STEP=STEP/6371.0;LAMI=LAMI/13080.0
       MUI=MUI/13080.0;RHOI=RHOI/13.08848;GRAVI=GRAVI/981.56

      RETURN
      END
