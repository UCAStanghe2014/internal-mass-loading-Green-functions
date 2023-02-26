      SUBROUTINE LOVE1(DEP,DEFORM)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8 LAMDAC,LAMI,MUI,INITS,INITST,INITL,MU,LAMDA
          REAL*8 minstep
      PARAMETER (MSTEP = 10000)
      DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &,GRAVITY(MSTEP)
      COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON    /MDI/ RADI, STEP, RHOI, LAMI, MUI, GRAVI
      COMMON    /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON    /INT1/ II
      COMMON    /INT2/ JJ
      COMMON    /INT3/ IJ
      common /info_source/ grs

      DIMENSION SUMS(6,6),SUMST(2,2)
      DIMENSION SUML(2,2),SUML1(2,2)
      DIMENSION SUMS1(6,6),SUMST1(2,2)
      DIMENSION SUMS2(6,6),SUMST2(2,2)
      DIMENSION SUMS3(6,6),SUMST3(2,2)
      DIMENSION SUMS4(6,6),SUMST4(2,2)
      DIMENSION SUMS5(6,6),SUMST5(2,2)
      DIMENSION INITS(6,3),INITS1(6,3),INITST(2,1),INITL(2,1)
      DIMENSION SSS(6,1),SVD(6,1),SVT(6,1),SHT(6,1),ZS(6,1)
      DIMENSION SSST(2,1),SVDT(2,1),SVTT(2,1),SHTT(2,1),ZST(2,1)
      DIMENSION SURFA1(6,3),SURFA2(6,3),SURFA3(6,3)
      DIMENSION AAA(3,3),BBB(3),COEFS(3)
      DIMENSION SOURCE(6,1),SOURCE1(6,1),YLOVE(6),YLOVE1(6)
        ! minstep = 0.01d0  ! default value
C       minstep =1d-2
        call int_step(DEP,minstep)
        !WRITE(*,*) 'minstep',minstep

      als = (6371d3)**2 * 1308d9

      N=1
      PRINT 111, N
 111  FORMAT ('   DEGREE OF INTEGRAL = ', 1I8)

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
          IF (II.EQ.M .AND. JJ.EQ.1)
     &    CALL INITIAL (INITS,INITST,SUMS,SUMST,N)
          CALL SOLIDST(N,SUMS,SUMST,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
        END DO
      END DO

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
          IF (II.EQ.1 .AND. JJ.EQ.1)
     &    CALL BNDRYSL(SUML,INITL,INITS,SUMS)
          CALL LIQUID(N,SUML,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
        END DO
      END DO
      CALL BNDRYLS(SUML,INITL,INITS1,SUMS)

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
          IF (II.EQ.1 .AND. JJ.EQ.1)
     &    CALL INITIAL (INITS,INITST,SUMS,SUMST,N)

          CALL SOLIDST(N,SUMS,SUMST,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
          IF(ABS(6371.0-RADI*6371.0-DEFORM) .LT. 0.01) THEN
            CALL MATRIX1(6,6,SUMS,SUMS4)
            CALL MATRIX1(2,2,SUMST,SUMST4)
          END IF

          IF(ABS(6371.0-RADI*6371.0-DEP).LT.0.01) THEN
            grs = GRAVI*9.8156d0
            s_gravi  = GRAVI
            s_source = RADI
            CALL MATRIX1(6,6,SUMS,SUMS2)
            CALL MATRIX1(2,2,SUMST,SUMST2)
            CALL DISLO(2,N,SVD,SVDT,SUMS,SUMST)
            CALL DISLO(3,N,SVT,SVTT,SUMS,SUMST)
            CALL DISLO(4,N,SHT,SHTT,SUMS,SUMST)
          END IF

        END DO
      END DO

      CALL MATRIX1(6,6,SUMS,SUMS5)
      CALL MATRIX1(2,2,SUMST,SUMST5)


C+++++++++++++Solve the constant coefficient+++++++++++++++++++++
      CALL MATRIX(6,6,3,SUMS2,INITS1,SURFA1)
C     CALL MATRIX(6,6,3,SUMS2,INITS,SURFA1)
      CALL MATRIX(6,6,3,SUMS5,SURFA1,SURFA2)

      DO IJS=2,4
        IF(IJS.EQ.2) CALL MATRIX(6,6,1,SUMS5,SVD,SOURCE)
        IF(IJS.EQ.3) CALL MATRIX(6,6,1,SUMS5,SVT,SOURCE)
        IF(IJS.EQ.4) CALL MATRIX(6,6,1,SUMS5,SHT,SOURCE)

        ! 两组基本积分解
        DO J=1,2
          AAA(1,J) = SURFA2(2,J)
          AAA(2,J) = SURFA2(4,J)
          AAA(3,J) = SURFA2(5,J)
        END DO
        ! 第三组解为刚性平移解
        AAA(1,3) = 0.d0
        AAA(2,3) = 0.d0
        AAA(3,3) = 1.d0

        ! z2(a)=z4(a)=0, and z5(a)=s_source/g1(s_source)/(4*pi*constantA)
        BBB(1) = -1.0*SOURCE(2,1) + 0d0
        BBB(2) = -1.0*SOURCE(4,1) + 0d0
        !BBB(3) = -1.0*SOURCE(5,1) + ((RADI/s_gravi)/(CONSTA*4D0*PAI)) ! old
        BBB(3) = -1.0*SOURCE(5,1) + RADI
        CALL GS4(3,AAA,BBB,0.1D-30,KWJI)

        ! 前2组解的系数
        DO J=1,2
          COEFS(J) = BBB(J)
        END DO
        COEFS(3)  =  0.d0

        ! 3个组合系数
        YLOVE1(1) = BBB(3)
        YLOVE1(3) = BBB(3)
        YLOVE1(5) = BBB(3)

        YLOVE1(2) = 0.d0
        YLOVE1(4) = 0.d0
        YLOVE1(6) = 0.d0

        VALUES=((6371.0-DEFORM)/(6371.0-DEP))**N

        IF(DEP .LT. DEFORM) THEN
          CALL MATRIX(6,6,3,SUMS4,SURFA1,SURFA3)
          CALL MATRIX(6,3,1,SURFA3, COEFS, YLOVE)
          DO I=1,6
            YLOVE(I)=(YLOVE(I)+YLOVE1(I))*VALUES
          END DO
        END IF

        IF(DEP .GT. DEFORM) THEN
          CALL MATRIX(6,6,3,SUMS4,SURFA1,SURFA3)
          CALL MATRIX(6,3,1,SURFA3, COEFS, YLOVE)
          IF(IJS.EQ.2) CALL MATRIX(6,6,1,SUMS4,SVD,SOURCE1)
          IF(IJS.EQ.3) CALL MATRIX(6,6,1,SUMS4,SVT,SOURCE1)
          IF(IJS.EQ.4) CALL MATRIX(6,6,1,SUMS4,SHT,SOURCE1)
          DO I=1,6
            YLOVE(I)=(YLOVE(I)+YLOVE1(I)+SOURCE1(I,1))/VALUES
          END DO
        END IF

        YLOVET1=0.0
        YLOVET2=0.0

        IF(IJS .EQ. 2) THEN
          SB   = 3438.3D0 / 6371.D0
          DPTH = (6371.d0-DEP) / 6371.d0
          TOR  = - 3.D0 / 16.D0 / 3.1415926D0 / DPTH**3 * (DPTH**5
     &    -SB**5) / (1.D0 - SB**5)
          C1=TOR/6371.0
          C2=(1-DPTH**5)/(SB**5-DPTH**5)*C1
          IF(DEP.GT.DEFORM) YLOVET1=C1*(6371.d0-DEFORM)
          IF(DEP.LT.DEFORM) YLOVET1=C2*(6371.d0-DEFORM)
          YLOVET2=0.0
        END IF
        ! 环型解=0
        YLOVET1 = 0D0
        YLOVET2 = 0D0

        !// z-numbers to y-numbers of unit mass loading by a factor
        call z2y_unit_mass_loading(YLOVE,YLOVET1,YLOVET2)

        YLOVE(5)=6.67D-11/6371D3*(1D0-DEP/6371D0)

        !// write  y-numbers to file
        !// Without any regularization, dimension is consistent with the physical meaning of y.
        WRITE(700+IJS,'(1I10,8(2X,1E20.10E5))') N,
     &  (YLOVE(I),I=1,6),YLOVET1,YLOVET2
      END DO

      END



