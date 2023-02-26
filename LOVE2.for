      SUBROUTINE LOVE2(ND,DEP,DEFORM)
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
        DIMENSION SURFA1(6,3),SURFAT1(2,1)
        DIMENSION SURFA2(6,3),SURFAT2(2,1)
        DIMENSION SURFA3(6,6),SURFAT3(2,2)
        DIMENSION SURFA4(6,3),SURFAT4(2,1)
        DIMENSION SURFA5(6,3),SURFAT5(2,1)
        DIMENSION AAA(6,6),BBB(6),AAAT(2,2),BBBT(2),COEFS1(3),COEFS2(3)
        DIMENSION SURINITS(6,3),SURINITST(2,1)
        DIMENSION SOURCE(6,1),SOURCET(2,1),YLOVE(6)
        ! minstep = 0.01d0  ! default value
C       minstep =5d-2
        call int_step(DEP,minstep)
        !WRITE(*,*) 'minstep',minstep


          als = (6371d3)**2 * 1308d9

          N=ND
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
            IF (II.EQ.M .AND. JJ.EQ.1)
     &      CALL INITIAL (INITS,INITST,SUMS,SUMST,N)
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
     &      CALL BNDRYSL(SUML,INITL,INITS,SUMS)
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
     &      CALL INITIAL (INITS,INITST,SUMS,SUMST,N)
         CALL SOLIDST(N,SUMS,SUMST,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
           IF(ABS(6371.0-RADI*6371.0-DEFORM).LT.0.01) THEN
             CALL MATRIX1(6,6,SUMS,SUMS4)
             CALL MATRIX1(2,2,SUMST,SUMST4)
           END IF
           IF(ABS(6371.0-RADI*6371.0-DEP).LT.0.01) THEN
             CALL MATRIX1(6,6,SUMS,SUMS2)
             CALL MATRIX1(2,2,SUMST,SUMST2)
             CALL DISLO(1,N,SSS,SSST,SUMS,SUMST)
             CALL DISLO(2,N,SVD,SVDT,SUMS,SUMST)
             CALL DISLO(3,N,SVT,SVTT,SUMS,SUMST)
             CALL DISLO(4,N,SHT,SHTT,SUMS,SUMST)
           END IF
         END DO
        END DO

C+++++++++From the surface to the core-mantle boundary++++++++++++++
      OPEN(104,FILE='model3_mantle.dat',ACTION='READ')
           READ(104,*) M
         DO I=1,2*M
           READ(104,*) RADIUS(2*M+1-I),GRAVITY(2*M+1-I)
     &       ,RHO(2*M+1-I),LAMDA(2*M+1-I),MU(2*M+1-I)
           END DO
      CLOSE(104)
         DO II = 1, M
            IJ=ceiling((RADIUS(2*II-1)-RADIUS(2*II))/minstep)
           DO JJ= 1,IJ
          CALL DATASET(1,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
           IF (II.EQ.1 .AND. JJ.EQ.1) CALL UNIT(SUMS,6)
           IF (II.EQ.1 .AND. JJ.EQ.1) CALL UNIT(SUMST,2)
          CALL SOLIDST1(N,SUMS,SUMST,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
           IF(ABS(6371.0-RADI*6371.0-DEFORM).LT.0.01) THEN
             CALL MATRIX1(6,6,SUMS,SUMS5)
             CALL MATRIX1(2,2,SUMST,SUMST5)
           END IF
           IF(DEFORM.LT.0.01) CALL UNIT(SUMS5,6)
           IF(DEFORM.LT.0.01) CALL UNIT(SUMST5,2)
           IF(ABS(6371.0-RADI*6371.0-DEP).LT.0.01) THEN
           grs = GRAVI*9.8156d0
           CALL MATRIX1(6,6,SUMS,SUMS3)
               CALL MATRIX1(2,2,SUMST,SUMST3)
         END IF
         END DO
         END DO

C+++++++++++++Solve the constant coefficient+++++++++++++++++++++
        IF(N.LE.1000) CALL MATRIX(6,6,3,SUMS2,INITS1,SURFA1)
        IF(N.GT.1000) CALL MATRIX(6,6,3,SUMS2,INITS,SURFA1)
C         CALL MATRIX(6,6,3,SUMS2,INITS,SURFA1)
         CALL MATRIX(2,2,1,SUMST2,INITST,SURFAT1)
         CALL ZERO(SURINITS, 6,3)
         CALL ZERO(SURINITST, 2,1)
            SURINITS(1,1)=1.0;SURINITS(3,2)=1.0
          SURINITS(5,3)=1.0;SURINITST(1,1)=1.0
         CALL MATRIX(6,6,3,SUMS3,SURINITS,SURFA2)
         CALL MATRIX(2,2,1,SUMST3,SURINITST,SURFAT2)

          DO J=1,3
           DO I=1,6
            SURFA3(I,J)=SURFA1(I,J)
            SURFA3(I,J+3)=SURFA2(I,J)
           END DO
          END DO
          DO I=1,2
            SURFAT3(I,1)=SURFAT1(I,1)
            SURFAT3(I,2)=SURFAT2(I,1)
          END DO

         DO IJS=1,4
            IF(IJS.EQ.1) CALL MATRIX1(6,1,SSS,SOURCE)
            IF(IJS.EQ.1) CALL MATRIX1(2,1,SSST,SOURCET)
            IF(IJS.EQ.2) CALL MATRIX1(6,1,SVD,SOURCE)
            IF(IJS.EQ.2) CALL MATRIX1(2,1,SVDT,SOURCET)
            IF(IJS.EQ.3) CALL MATRIX1(6,1,SVT,SOURCE)
            IF(IJS.EQ.3) CALL MATRIX1(6,1,SVTT,SOURCET)
            IF(IJS.EQ.4) CALL MATRIX1(6,1,SHT,SOURCE)
            IF(IJS.EQ.4) CALL MATRIX1(6,1,SHTT,SOURCET)

              DO I=1,6
               DO J=1,6
                 AAA(I,J)=SURFA3(I,J)
             END DO
               BBB(I)=SOURCE(I,1)
            END DO
              CALL GS4(6,AAA,BBB,.1D-30,KWJI)
              DO J=1,3
              COEFS1(J)=-1.0*BBB(J)
              COEFS2(J)=BBB(J+3)
            END DO

              DO I=1,2
               DO J=1,2
                 AAAT(I,J)=SURFAT3(I,J)
             END DO
               BBBT(I)=SOURCET(I,1)
            END DO
              CALL GS4(2,AAAT,BBBT,.1D-30,KWJI)
              COEFST1=-1.0*BBBT(1);COEFST2=BBBT(2)

             VALUES=((6371.0-DEFORM)/(6371.0-DEP))**N
           IF(DEP.LT.DEFORM) THEN
        IF(N.LE.1000) CALL MATRIX(6,6,3,SUMS4,INITS1,SURFA4)
        IF(N.GT.1000) CALL MATRIX(6,6,3,SUMS4,INITS,SURFA4)
C             CALL MATRIX(6,6,3,SUMS4,INITS,SURFA4)
             CALL MATRIX(2,2,1,SUMST4,INITST,SURFAT4)
             CALL MATRIX(6,3,1,SURFA4, COEFS1, YLOVE)
             DO I=1,6
             YLOVE(I)=YLOVE(I)*VALUES
             END DO
             YLOVET1=SURFAT4(1,1)*COEFST1*VALUES
             YLOVET2=SURFAT4(2,1)*COEFST1*VALUES
           END IF
           IF(DEP.GT.DEFORM) THEN
             CALL MATRIX(6,6,3,SUMS5,SURINITS,SURFA5)
             CALL MATRIX(2,2,1,SUMST5,SURINITST,SURFAT5)
             CALL MATRIX(6,3,1,SURFA5, COEFS2, YLOVE)
             DO I=1,6
             YLOVE(I)=YLOVE(I)/VALUES
             END DO
             YLOVET1=SURFAT5(1,1)*COEFST2/VALUES
             YLOVET2=SURFAT5(2,1)*COEFST2/VALUES
           END IF
            !// z-numbers to y-numbers of unit mass loading by a factor
            call z2y_unit_mass_loading(YLOVE,YLOVET1,YLOVET2)

          !// write  y-numbers to file
          !// Without any regularization, dimension is consistent with the physical meaning of y.
            WRITE(700+IJS,'(1I10,8(2X,1E20.10E5))') N,
     &          (YLOVE(I),I=1,6),YLOVET1,YLOVET2
         END DO

      END

