      SUBROUTINE LOVE(DEP,DEFORM,M)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      CHARACTER*30 INPUTFILE(5),OUTPUTFILE(5)
      DIMENSION MDEGREE(1000),NDEGREE(100000)
      DIMENSION DEGREEM(1000),DEGREEN(100000)
      DIMENSION YVT1(100),YVT2(100),YHT1(100),YHT2(100)
      DIMENSION YVD(8),YVT(6),YHT(6)
      DIMENSION Y1(1000),Y2(1000),YY1(100000),YY2(100000)
      DIMENSION Y3(1000),Y4(1000),YY3(100000),YY4(100000)
      DIMENSION Y5(1000),Y6(1000),YY5(100000),YY6(100000)
      DIMENSION YT1(1000),YT2(1000),YYT1(100000),YYT2(100000)

      INPUTFILE(1)='LOVE12.DAT' ; OUTPUTFILE(1)='DSLOV12.DAT'
      INPUTFILE(2)='LOVE32.DAT' ; OUTPUTFILE(2)='DSLOV32.DAT'
      INPUTFILE(3)='LOVE220.DAT'; OUTPUTFILE(3)='DSLOV220.DAT'
      INPUTFILE(4)='LOVE33.DAT' ; OUTPUTFILE(4)='DSLOV33.DAT'

      DO  ITYPE = 1,4
        OPEN(100+ITYPE,FILE=INPUTFILE(ITYPE),ACTION='READ')
        OPEN(200+ITYPE,FILE=OUTPUTFILE(ITYPE),ACTION='READ')
CCCCCCCCCCCCCC  DEGREE=0 OR 1  CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (ITYPE.EQ.2) THEN
          READ(100+ITYPE,*) NVD1,(YVD(J),J=1,8)
          WRITE(200+ITYPE,753) NVD1,(YVD(J),J=1,8)
        END IF
        IF (ITYPE.EQ.3) THEN
          READ(100+ITYPE,*) NVT1,(YVT(J),J=1,6)
          WRITE(200+ITYPE,753) NVT1,(YVT(J),J=1,6)
          READ(100+ITYPE,*) NVT1,(YVT(J),J=1,6)
          WRITE(200+ITYPE,753) NVT1,(YVT(J),J=1,6)
        END IF
        IF (ITYPE.EQ.4) THEN
          READ(100+ITYPE,*) NHT1,(YHT(J),J=1,6)
          WRITE(200+ITYPE,753) NHT1,(YHT(J),J=1,6)
          READ(100+ITYPE,*) NHT1,(YHT(J),J=1,6)
          WRITE(200+ITYPE,753) NHT1,(YHT(J),J=1,6)
        END IF

CCCCCCCCCCCCCC  DEGREE>=2 FOR SS12 CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IF (ITYPE.EQ.1) THEN
        DO I=1,M-1
          READ(100+ITYPE,*) MDEGREE(I),Y1(I),Y2(I)
     &    ,Y3(I),Y4(I),Y5(I),Y6(I),YT1(I),YT2(I)
          DEGREEM(I)=MDEGREE(I)*1.0D0
          FACTOR1=(6371.0D0/(6371.0D0-ABS(DEP-DEFORM)))**MDEGREE(I)
          Y1(I) =Y1(I)*FACTOR1
          Y2(I) =Y2(I)/MDEGREE(I)*FACTOR1
          Y3(I) =Y3(I)*MDEGREE(I)*FACTOR1
          Y4(I) =Y4(I)*FACTOR1
          Y5(I) =Y5(I)*MDEGREE(I)*FACTOR1
          Y6(I) =Y6(I)*FACTOR1
          YT1(I)=YT1(I)*MDEGREE(I)**2*FACTOR1
          YT2(I)=YT2(I)*MDEGREE(I)*FACTOR1
        END DO
        N=MDEGREE(M-1)
        DO I=1,N-1
          NDEGREE(I)=I+1
          DEGREEN(I)=(I+1)*1.0D0
        END DO
        CALL SPLINE3(M-1,DEGREEM,Y1, N-1,DEGREEN,YY1)
        CALL SPLINE3(M-1,DEGREEM,Y2, N-1,DEGREEN,YY2)
        CALL SPLINE3(M-1,DEGREEM,Y3, N-1,DEGREEN,YY3)
        CALL SPLINE3(M-1,DEGREEM,Y4, N-1,DEGREEN,YY4)
        CALL SPLINE3(M-1,DEGREEM,Y5, N-1,DEGREEN,YY5)
        CALL SPLINE3(M-1,DEGREEM,Y6, N-1,DEGREEN,YY6)
        CALL SPLINE3(M-1,DEGREEM,YT1,N-1,DEGREEN,YYT1)
        CALL SPLINE3(M-1,DEGREEM,YT2,N-1,DEGREEN,YYT2)
        DO I=1,N-1
C          WRITE(301,753) NDEGREE(I),YY1(I),YY2(I)
C      &   ,YY3(I),YY4(I),YY5(I),YY6(I),YYT1(I),YYT2(I)
              FACTOR2=((6371.0D0-ABS(DEP-DEFORM))/6371.0D0)**NDEGREE(I)
              YY1(I) =YY1(I)*FACTOR2
              YY2(I) =YY2(I)*NDEGREE(I)*FACTOR2
              YY3(I) =YY3(I)/NDEGREE(I)*FACTOR2
              YY4(I) =YY4(I)*FACTOR2
              YY5(I) =YY5(I)/NDEGREE(I)*FACTOR2
              YY6(I) =YY6(I)*FACTOR2
              YYT1(I)=YYT1(I)/NDEGREE(I)**2*FACTOR2
              YYT2(I)=YYT2(I)/NDEGREE(I)*FACTOR2
              WRITE(200+ITYPE,753) NDEGREE(I),YY1(I),YY2(I)
     &        ,YY3(I),YY4(I),YY5(I),YY6(I),YYT1(I),YYT2(I)
            END DO
          END IF


CCCCCCCCCCCCCC  DEGREE>=2 FOR VD32 CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IF (ITYPE.EQ.2) THEN
          DO I=1,M-1
            READ(100+ITYPE,*) MDEGREE(I),Y1(I),Y2(I)
     &      ,Y3(I),Y4(I),Y5(I),Y6(I),YT1(I),YT2(I)
            DEGREEM(I)=MDEGREE(I)*1.0D0
            FACTOR1=(6371.0D0/(6371.0D0-ABS(DEP-DEFORM)))**MDEGREE(I)
            Y1(I) =Y1(I)/MDEGREE(I)*FACTOR1
            Y2(I) =Y2(I)/MDEGREE(I)**2*FACTOR1
            Y3(I) =Y3(I)*FACTOR1
            Y4(I) =Y4(I)/MDEGREE(I)*FACTOR1
            Y5(I) =Y5(I)*FACTOR1
            Y6(I) =Y6(I)/MDEGREE(I)*FACTOR1
            YT1(I)=YT1(I)*MDEGREE(I)*FACTOR1
            YT2(I)=YT2(I)*FACTOR1
          END DO
          DO I=1,N-1
            NDEGREE(I)=I+1
            DEGREEN(I)=(I+1)*1.0D0
          END DO
          CALL SPLINE3(M-1,DEGREEM,Y1, N-1,DEGREEN,YY1)
          CALL SPLINE3(M-1,DEGREEM,Y2, N-1,DEGREEN,YY2)
          CALL SPLINE3(M-1,DEGREEM,Y3, N-1,DEGREEN,YY3)
          CALL SPLINE3(M-1,DEGREEM,Y4, N-1,DEGREEN,YY4)
          CALL SPLINE3(M-1,DEGREEM,Y5, N-1,DEGREEN,YY5)
          CALL SPLINE3(M-1,DEGREEM,Y6, N-1,DEGREEN,YY6)
          CALL SPLINE3(M-1,DEGREEM,YT1,N-1,DEGREEN,YYT1)
          CALL SPLINE3(M-1,DEGREEM,YT2,N-1,DEGREEN,YYT2)
          DO I=1,N-1
C          WRITE(302,753) NDEGREE(I),YY1(I),YY2(I)
C      &   ,YY3(I),YY4(I),YY5(I),YY6(I),YYT1(I),YYT2(I)
                FACTOR2=((6371.0D0-ABS(DEP-DEFORM))/6371.0D0)**NDEGREE(I)
                YY1(I) =YY1(I)*NDEGREE(I)*FACTOR2
                YY2(I) =YY2(I)*NDEGREE(I)**2*FACTOR2
                YY3(I) =YY3(I)*FACTOR2
                YY4(I) =YY4(I)*NDEGREE(I)*FACTOR2
                YY5(I) =YY5(I)*FACTOR2
                YY6(I) =YY6(I)*NDEGREE(I)*FACTOR2
                YYT1(I)=YYT1(I)/NDEGREE(I)*FACTOR2
                YYT2(I)=YYT2(I)*FACTOR2
                WRITE(200+ITYPE,753) NDEGREE(I),YY1(I),YY2(I)
     &          ,YY3(I),YY4(I),YY5(I),YY6(I),YYT1(I),YYT2(I)
              END DO
            END IF


CCCCCCCCCCCCCC  DEGREE>=2 FOR VT220 CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IF (ITYPE.EQ.3) THEN
            DO I=1,M-1
              READ(100+ITYPE,*) MDEGREE(I),Y1(I),Y2(I)
     &        ,Y3(I),Y4(I),Y5(I),Y6(I)
              DEGREEM(I)=MDEGREE(I)*1.0D0
              FACTOR1=(6371.0D0/(6371.0D0-ABS(DEP-DEFORM)))**MDEGREE(I)
              Y1(I) =Y1(I)/MDEGREE(I)**2*FACTOR1
              Y2(I) =Y2(I)/MDEGREE(I)**3*FACTOR1
              Y3(I) =Y3(I)/MDEGREE(I)*FACTOR1
              Y4(I) =Y4(I)/MDEGREE(I)**2*FACTOR1
              Y5(I) =Y5(I)/MDEGREE(I)*FACTOR1
              Y6(I) =Y6(I)/MDEGREE(I)**2*FACTOR1
            END DO
            DO I=1,N-1
              NDEGREE(I)=I+1
              DEGREEN(I)=(I+1)*1.0D0
            END DO
            CALL SPLINE3(M-1,DEGREEM,Y1, N-1,DEGREEN,YY1)
            CALL SPLINE3(M-1,DEGREEM,Y2, N-1,DEGREEN,YY2)
            CALL SPLINE3(M-1,DEGREEM,Y3, N-1,DEGREEN,YY3)
            CALL SPLINE3(M-1,DEGREEM,Y4, N-1,DEGREEN,YY4)
            CALL SPLINE3(M-1,DEGREEM,Y5, N-1,DEGREEN,YY5)
            CALL SPLINE3(M-1,DEGREEM,Y6, N-1,DEGREEN,YY6)
            DO I=1,N-1
C          WRITE(303,753) NDEGREE(I),YY1(I),YY2(I)
C      &   ,YY3(I),YY4(I),YY5(I),YY6(I)
                  FACTOR2=((6371.0D0-ABS(DEP-DEFORM))/6371.0D0)**NDEGREE(I)
                  YY1(I) =YY1(I)*NDEGREE(I)**2*FACTOR2
                  YY2(I) =YY2(I)*NDEGREE(I)**3*FACTOR2
                  YY3(I) =YY3(I)*NDEGREE(I)*FACTOR2
                  YY4(I) =YY4(I)*NDEGREE(I)**2*FACTOR2
                  YY5(I) =YY5(I)*NDEGREE(I)*FACTOR2
                  YY6(I) =YY6(I)*NDEGREE(I)**2*FACTOR2
                  WRITE(200+ITYPE,753) NDEGREE(I),YY1(I),YY2(I)
     &            ,YY3(I),YY4(I),YY5(I),YY6(I)
                END DO
              END IF


CCCCCCCCCCCCCC  DEGREE>=2 FOR HT33 CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            IF (ITYPE.EQ.4) THEN
              DO I=1,M-1
                READ(100+ITYPE,*) MDEGREE(I),Y1(I),Y2(I)
     &          ,Y3(I),Y4(I),Y5(I),Y6(I)
                DEGREEM(I)=MDEGREE(I)*1.0D0
                FACTOR1=(6371.0D0/(6371.0D0-ABS(DEP-DEFORM)))**MDEGREE(I)
                Y1(I) =Y1(I)/MDEGREE(I)**2*FACTOR1
                Y2(I) =Y2(I)/MDEGREE(I)**3*FACTOR1
                Y3(I) =Y3(I)/MDEGREE(I)*FACTOR1
                Y4(I) =Y4(I)/MDEGREE(I)**2*FACTOR1
                Y5(I) =Y5(I)/MDEGREE(I)*FACTOR1
                Y6(I) =Y6(I)/MDEGREE(I)**2*FACTOR1
              END DO
              DO I=1,N-1
                NDEGREE(I)=I+1
                DEGREEN(I)=(I+1)*1.0D0
              END DO
              CALL SPLINE3(M-1,DEGREEM,Y1, N-1,DEGREEN,YY1)
              CALL SPLINE3(M-1,DEGREEM,Y2, N-1,DEGREEN,YY2)
              CALL SPLINE3(M-1,DEGREEM,Y3, N-1,DEGREEN,YY3)
              CALL SPLINE3(M-1,DEGREEM,Y4, N-1,DEGREEN,YY4)
              CALL SPLINE3(M-1,DEGREEM,Y5, N-1,DEGREEN,YY5)
              CALL SPLINE3(M-1,DEGREEM,Y6, N-1,DEGREEN,YY6)
              DO I=1,N-1
C          WRITE(304,753) NDEGREE(I),YY1(I),YY2(I)
C      &   ,YY3(I),YY4(I),YY5(I),YY6(I)
                    FACTOR2=((6371.0D0-ABS(DEP-DEFORM))/6371.0D0)**NDEGREE(I)
                    YY1(I) =YY1(I)*NDEGREE(I)**2*FACTOR2
                    YY2(I) =YY2(I)*NDEGREE(I)**3*FACTOR2
                    YY3(I) =YY3(I)*NDEGREE(I)*FACTOR2
                    YY4(I) =YY4(I)*NDEGREE(I)**2*FACTOR2
                    YY5(I) =YY5(I)*NDEGREE(I)*FACTOR2
                    YY6(I) =YY6(I)*NDEGREE(I)**2*FACTOR2
                    WRITE(200+ITYPE,753) NDEGREE(I),YY1(I),YY2(I)
     &              ,YY3(I),YY4(I),YY5(I),YY6(I)
                  END DO
                END IF

                CLOSE(100+ITYPE)
                CLOSE(200+ITYPE)
              END DO

753           FORMAT(1I10,7(2X,1E25.15e5))
            END



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C    ***************************************************************************
C    *   THE CUBIC SPLINE INTERPOLATION BASED ON THE THREE MOMENT FUNCTION     *
C    *                     FREE BOUNDARY CONDITION                             *
C    ***************************************************************************
C
            SUBROUTINE SPLINE3(N,X,Y,M,T,U)
            IMPLICIT REAL*8(A-H,O-Z)
            IMPLICIT INTEGER*8(I-N)
            DIMENSION X(N),Y(N),T(M),U(M),A(N),B(N),C(N),D(N),UU(M),V(M)

            A(1)= 0.0d0
            A(N)= 1.0d0
            D(1)= 0.0d0
            D(N)= 0.0d0
            C(1)=-1.0d0
            C(N)= 0.0d0
            B(1)= 1.0d0
            B(N)=-1.0d0
            L   =N-1

            DO 5 I=2,L
              A(I)=(X(I)-X(I-1))/6.0d0
              C(I)=(X(I+1)-X(I))/6.0d0
              B(I)=2.0d0*(A(I)+C(I))
5           D(I)=(Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))/(X(I)-X(I-1))
            C(1)=C(1)/B(1)
            D(1)=D(1)/B(1)

            DO 10 I=2,N
              C(I)=C(I)/(B(I)-A(I)*C(I-1))
10          D(I)=(D(I)-A(I)*D(I-1))/(B(I)-A(I)*C(I-1))
            A(N)=D(N)

            DO 15 I=1,L
              J=N-I
15          A(J)=D(J)-C(J)*A(J+1)

            DO  30 J1=1,M
              F=T(J1)
              DO 20 J2=1,N-1
                IF(X(J2).LE.F.AND.F.LE.X(J2+1)) GOTO 25
20            CONTINUE
              GOTO 30
25            E=X(J2+1)-X(J2)
              RR=(A(J2)*(X(J2+1)-F)**3+A(J2+1)*(F-X(J2))**3)/6.0d0/E
              SS=(X(J2+1)-F)*(Y(J2)/E-A(J2)*E/6.0d0)
              TT=(F-X(J2))*(Y(J2+1)/E-A(J2+1)*E/6.0d0)
              U(J1)=RR+SS+TT

30          CONTINUE

            RETURN
            END
