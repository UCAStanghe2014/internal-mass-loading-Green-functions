      SUBROUTINE LOVE0(DEP,DEFORM,DLA,DMU)
        IMPLICIT REAL*8(A-H,O-Z)
        IMPLICIT INTEGER*8(I-N)
        REAL*8 LAMDAC,LAMI,MUI,INITS,INITST,INITL,MU,LAMDA
        REAL*8 minstep
        PARAMETER (MSTEP = 10000)
        DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
        COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
        COMMON    /MDI/   RADI, STEP, RHOI, LAMI, MUI, GRAVI
        COMMON    /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
        COMMON    /INT1/ II
        COMMON    /INT2/ JJ
        COMMON    /INT3/ IJ
        common    /info_source/ grs

        DIMENSION SUMS(2,2),SUML(2,2),SUML1(2,2)
        DIMENSION SUMS1(2,2),SUMS2(2,2),SUMS3(2,2),SUMS4(2,2),SUMS5(2,2)
        DIMENSION INITS(2,1),INITS1(2,1),INITL(2,1)
        DIMENSION SVT(2,1),SHT(2,1)
        DIMENSION SURFA1(2,1),SURFA2(2,1),SURFA3(2,2)
        DIMENSION AAA(2,2),BBB(2),SURFA4(2,1),SURFA5(2,1)
        DIMENSION SURINITS(2,1),SOURCE(2,1),YLOVE(6)
C       minstep = 0.01d0  ! default value
C       minstep = 0.15
C       minstep =5d-2

        call int_step(DEP,minstep)

        !WRITE(*,*) 'minstep',minstep

        als = (6371d3)**2 * 1308d9


        N=0
        PRINT 111, N
 111    FORMAT ('   DEGREE OF INTEGRAL = ', 1I8)

C+++++++++++++++++++INNER CORE+++++++++++++++++++++++++++++++++++
        OPEN(101,FILE='model1_incore.dat',ACTION='READ')
        READ(101,*) M
        DO I=1,2*M
          READ(101,*) RADIUS(I),GRAVITY(I),RHO(I),LAMDA(I),MU(I)
        END DO
        CLOSE(101)
        DO II = M, M
          IJ=ceiling((RADIUS(2*II)-RADIUS(2*II-1))/minstep)
          DO JJ= 1,IJ
            CALL DATASET(1,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            IF (II.EQ.M .AND. JJ.EQ.1)  CALL INITI0 (INITS,SUMS,N)
            CALL SOLID0(SUMS,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
          END DO
        END DO
        CALL MATRIX(2,2,1,SUMS,INITS,INITL)
        CALL UNIT(SUML, 2)

C+++++++++++++++++++OUTER CORE+++++++++++++++++++++++++++++++++++
        OPEN(102,FILE='model2_outcore.dat',ACTION='READ')
        READ(102,*) M
        DO I=1,2*M
          READ(102,*) RADIUS(I),GRAVITY(I),RHO(I),LAMDA(I),MU(I)
        END DO
        CLOSE(102)
        DO II = 1, M
          IJ=ceiling((RADIUS(2*II)-RADIUS(2*II-1))/minstep)
          DO JJ= 1,IJ
            CALL DATASET(1,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            CALL LIQUID(N,SUML,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
          END DO
        END DO
        CALL MATRIX(2,2,1,SUML,INITL,INITS1)

C+++++++++From the core-mantle boundary to the surface++++++++++++
        OPEN(103,FILE='model3_mantle.dat',ACTION='READ')
        READ(103,*) M
        DO I=1,2*M
          READ(103,*) RADIUS(I),GRAVITY(I),RHO(I),LAMDA(I),MU(I)
        END DO
        CLOSE(103)
        DO II = 1, M
          IJ=ceiling((RADIUS(2*II)-RADIUS(2*II-1))/minstep)
          DO JJ= 1,IJ
            CALL DATASET(1,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            IF (II.EQ.1 .AND. JJ.EQ.1) CALL INITI0 (INITS,SUMS,N)
            CALL SOLID0(SUMS,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            IF(ABS(6371.0-RADI*6371.0-DEFORM).LT.0.01) THEN
              CALL MATRIX1(2,2,SUMS,SUMS4)
              DLA=LAMI*13080.0;DMU=MUI*13080.0
            END IF
            IF(ABS(6371.0-RADI*6371.0-DEP).LT.0.01) THEN
              CALL MATRIX1(2,2,SUMS,SUMS2)
              CALL DISLO0(3,SVT,SUMS)
              CALL DISLO0(4,SHT,SUMS)
            END IF
          END DO
        END DO

C+++++++++From the surface to the core-mantle boundary++++++++++++++
        OPEN(104,FILE='model3_mantle.dat',ACTION='READ')
        READ(104,*) M
        DO I=1,2*M
          READ(104,*) RADIUS(2*M+1-I),GRAVITY(2*M+1-I)
     &      ,RHO(2*M+1-I),LAMDA(2*M+1-I),MU(2*M+1-I)
        END DO
        CLOSE(104)
        DO II = 1, M
          IJ=ceiling((RADIUS(2*II-1)-RADIUS(2*II))/minstep)
          DO JJ= 1,IJ
            CALL DATASET(1,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            IF (II.EQ.1 .AND. JJ.EQ.1) CALL UNIT(SUMS,2)
            CALL SOLID0(SUMS,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
            IF(ABS(6371.0-RADI*6371.0-DEFORM).LT.0.01) THEN
              CALL MATRIX1(2,2,SUMS,SUMS5)
            END IF
            IF(DEFORM.LT.0.01) CALL UNIT(SUMS5,2)
            IF(ABS(6371.0-RADI*6371.0-DEP).LT.0.01) THEN
              CALL MATRIX1(2,2,SUMS,SUMS3)
              grs = GRAVI*9.8156d0
            END IF
          END DO
        END DO

C+++++++++++++Solve the constant coefficient+++++++++++++++++++++
        CALL MATRIX(2,2,1,SUMS2,INITS1,SURFA1)
C       CALL MATRIX(2,2,1,SUMS2,INITS,SURFA1)
        CALL ZERO(SURINITS, 2,1)
        SURINITS(1,1)=1.0
        CALL MATRIX(2,2,1,SUMS3,SURINITS,SURFA2)

        DO I=1,2
          SURFA3(I,1)=SURFA1(I,1)
          SURFA3(I,2)=SURFA2(I,1)
        END DO

        DO IJS=3,4
          IF(IJS.EQ.3) CALL MATRIX1(2,1,SVT,SOURCE)
          IF(IJS.EQ.4) CALL MATRIX1(2,1,SHT,SOURCE)

          DO I=1,2
            DO J=1,2
              AAA(I,J)=SURFA3(I,J)
            END DO
            BBB(I)=SOURCE(I,1)
          END DO
          CALL GS4(2,AAA,BBB,.1D-30,KWJI)
          COEFS1=-1.0*BBB(1);COEFS2=BBB(2)

          IF(DEP.LT.DEFORM) THEN
            CALL MATRIX(2,2,1,SUMS4,INITS1,SURFA4)
C           CALL MATRIX(2,2,1,SUMS4,INITS,SURFA4)
            YLOVE1=SURFA4(1,1)*COEFS1
            YLOVE2=SURFA4(2,1)*COEFS1
          END IF
          IF(DEP.GT.DEFORM) THEN
            CALL MATRIX(2,2,1,SUMS5,SURINITS,SURFA5)
            YLOVE1=SURFA5(1,1)*COEFS2
            YLOVE2=SURFA5(2,1)*COEFS2
          END IF

C        !// z-numbers to y-numbers of unit mass loading by a factor
          factor = 6.67D-11/((6371d3**2)*9.8156d0)
          YLOVE1 = (YLOVE1*factor) * 6371d3
          YLOVE2 = (YLOVE2*factor) * 1308d9

          if (grs<1d-13) then
            write(*,*) ' Please chenge the integration step at line 25 of file "LOVE0.for". '
          end if
          !WRITE(*,*) 'grs = ',grs

          !WRITE(*,*) YLOVE1,YLOVE2
          !READ(*,*)

          !// write  y-numbers to file
          !// Without any regularization, dimension is consistent with the physical meaning of y.
          WRITE(700+IJS,'(1I10,10(2x,1E20.10E5))') N,
     &      YLOVE1,YLOVE2,
     &      0.0,0.0,6.67D-11/6371D3,0.0,0.0,0.0
        END DO

      END
