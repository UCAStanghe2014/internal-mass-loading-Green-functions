C  ################### SUBROUTINE ##################################
C
C      *************************************************
C      *                                               *
C      *    INITIAL VALUES FOR INTEGRATING OUT CORE    *
C      *   CALCULATED FROM THE RESULTS OF INNER CORE   *
C      *       BY IN-OUT CORE BOUNDARY CONDITION       *
C      *                                               *
C      *************************************************
C
      SUBROUTINE BNDRYSL(SUML,INITL,INITS,SUMS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    THIS SUBROUTINE PREPARES INITIAL VALUES USED FOR INTEGRATION
C    OF THE LIQUID CORE, BY TRANSFORMING THE INTEGRATION RESULTS
C    OF THE INNER CORE, THROUGH THE IN-OUT CORE BOUNDARY CONDITIONS.
C
C    THE 4 ARGUMENTS ARE INPUTS AND OUTPUTS:
C
C    INPUTS:
C       INITS(6,3) - INITIAL VALUES AT THE EARTH CENTER
C       SUMS(6,6)  - INTEGRATION RESULTS OF THE INNER CORE
C
C    OUTPUTS:
C       INITL(2,1) - INITIAL VALUES FOR LIQUID CORE
C       SUML(2,2)  - UNIT TENSOR TO BE INSTORE INTEGRATION RESULTS
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8    INITS,INITL,LAMI,MUI
      DIMENSION SUML(2,2),INITL(2,1),INITS(6,3),SUMS(6,6),WORK(6,3)
      COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      CALL MATRIX(6,6,3,SUMS,INITS,WORK)

      A = WORK(2,1)-CONSTC*RHOI*(GRAVI*WORK(1,1)-WORK(5,1))
      B = WORK(2,2)-CONSTC*RHOI*(GRAVI*WORK(1,2)-WORK(5,2))
      C = WORK(2,3)-CONSTC*RHOI*(GRAVI*WORK(1,3)-WORK(5,3))
      E = WORK(4,1)/WORK(4,3)
      F = WORK(4,2)/WORK(4,3)
      D = (A/C-E)/(B/C-F)
      G = -A/C+B/C*D

      INITL(1,1) = WORK(5,1)-D*WORK(5,2)+G*WORK(5,3)
      INITL(2,1) = WORK(6,1)-D*WORK(6,2)+G*WORK(6,3)
     &           +(WORK(2,1)-D*WORK(2,2)+G*WORK(2,3))/CONSTA/GRAVI

      CALL UNIT(SUML,2)
      END

C  ################### SUBROUTINE ##################################
C
C        **********************************************
C        *                                            *
C        *    INITIAL VALUES FOR INTEGRATING MANTLE   *
C        *   CALCULATED FROM THE RESULTS OF OUT CORE  *
C        *    BY BOUNDARY CONDITION OF CORE-MANTAL    *
C        *                                            *
C        **********************************************
C
              SUBROUTINE BNDRYLS(SUML,INITL,INITS,SUMS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C    THIS SUBROUTINE PREPARES INITIAL VALUES USED FOR INTEGRATION
C    OF THE MANTAL, BY TRANSFORMING THE INTEGRATION RESULTS
C    OF THE LIQUID CORE, THROUGH THE CORE-MANTAL BOUNDARY CONDITIONS.
C
C    THE 4 ARGUMENTS ARE INPUTS AND OUTPUTS:
C
C    INPUTS:
C       INITL(2,1) - INITIAL VALUES AT THE BOTOM OF THE LIQUID CORE
C       SUML(2,2)  - INTEGRATION RESULTS OF THE LIQUID CORE
C
C    OUTPUTS:
C       INITS(6,3) - INITIAL VALUES AT THE BOTTOM OF THE MANTAL
C       SUMS(6,6)  - UNIT TENSOR TO INSTORE THE INTEGRATION RESULTS
C                    OF THE MANTAL
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8    INITS,INITL,LAMI,MUI
      DIMENSION SUML(2,2),INITL(2,1),INITS(6,3),SUMS(6,6),WORK(2,1)
      COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      CALL MATRIX(2,2,1,SUML,INITL,WORK)

      CALL ZERO(INITS,6,3)

      INITS(2,1) = -CONSTC*RHOI*WORK(1,1)
      INITS(5,1) =  WORK(1,1)
      INITS(6,1) =  WORK(2,1)+CONSTB*RHOI/GRAVI*WORK(1,1)
      INITS(1,2) =  1.D0
      INITS(2,2) =  CONSTC*RHOI*GRAVI
      INITS(6,2) = -CONSTB*RHOI
      INITS(3,3) =  1.D0

      CALL UNIT(SUMS,6)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C        ********************************************
C        *                                          *
C        *     CALCULATION OF THE COEFFICENTS OF    *
C        *     THE EQUATIONS OF MOTION IN LIQUID    *
C        *          RESULT IN 'EQUIL(2,2)'          *
C        *                                          *
C        ********************************************
C
      SUBROUTINE COEFIL(EQUIL,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT   REAL*8(A-H,O-Z)
      REAL*8     LAMI,LAMDAC,MUI
      DIMENSION  EQUIL(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      EQUIL(1,1) = (CONSTB*RHOI/GRAVI-(N+1)/RADI)-N/RADI
      EQUIL(1,2) = 1.D0
      EQUIL(2,1) = 2.D0*(N-1)/RADI*CONSTB*RHOI/GRAVI
      EQUIL(2,2) = (N-1)/RADI-CONSTB*RHOI/GRAVI-N/RADI

      IF(N.EQ.0) THEN
        EQUIL(1,1) = -2.D0/RADI
        EQUIL(1,2) = 1.D0/LAMI
        EQUIL(2,1) = -4.D0*CONSTC*RHOI*GRAVI/RADI
        EQUIL(2,2) = 0.D0
      ENDIF

      DO 30 I = 1, 2
         DO 30 J = 1, 2
 30         EQUIL(I,J) = STEP*EQUIL(I,J)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C        ********************************************
C        *                                          *
C        *     CALCULATION OF THE COEFFICENTS OF    *
C        *     THE EQUATIONS OF MOTION IN SOLID     *
C        *          RESULT IN 'EQUIS(6,6)'          *
C        *                                          *
C        ********************************************
C
      SUBROUTINE COEFIS(EQUIS,EQUIST,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8     LAMI,MUI,LAMDAC
      DIMENSION  EQUIS(6,6), EQUIST(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
C
      CALL ZERO(EQUIS, 6,6)
      CALL ZERO(EQUIST,2,2)

      EQUIS(1,1) = -2.D0*LAMI/(LAMI+2.D0*MUI)/RADI-N/RADI
      EQUIS(1,2) = 1.D0/(LAMI+2.D0*MUI)
      EQUIS(1,3) = N*(N+1.D0)*LAMI/(LAMI+2.D0*MUI)/RADI
      EQUIS(2,1) = (-FREQ**2*RHOI*CONSTC*RADIS*RADI*RADI/GRAVIS
     &             -4.D0*CONSTC*RHOI*GRAVI*RADI+4.D0*MUI*(3.D0*LAMI
     &             +2.D0*MUI)/(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(2,2) = -4.D0*MUI/(LAMI+2.D0*MUI)/RADI-N/RADI
      EQUIS(2,3) = N*(N+1.D0)*(CONSTC*RHOI*GRAVI*RADI-2.D0*MUI
     &             *(3.D0*LAMI+2.D0*MUI)/(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(2,4) = N*(N+1.D0)/RADI
      EQUIS(2,5) = CONSTC*RHOI*(N+1.D0)/RADI
      EQUIS(2,6) = -CONSTC*RHOI
      EQUIS(3,1) = -1.D0/RADI
      EQUIS(3,3) = 1.D0/RADI-N/RADI
      EQUIS(3,4) = 1.D0/MUI
      EQUIS(4,1) = (CONSTC*RHOI*GRAVI*RADI-2.D0*MUI*(3.D0*LAMI+2.D0*MUI)
     &             /(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(4,2) = -LAMI/(LAMI+2.D0*MUI)/RADI
      EQUIS(4,3) = -FREQ**2*RHOI*CONSTC*RADIS/GRAVIS+2.D0*MUI
     &             /(LAMI+2.D0*MUI)
     &             /RADI/RADI*(LAMI*(2.D0*N*N+2.D0*N-1.D0)+2.D0*MUI
     &             *(N*N+N-1.D0))
      EQUIS(4,4) = -3.D0/RADI-N/RADI
      EQUIS(4,5) = -CONSTC*RHOI/RADI
      EQUIS(5,1) = CONSTB*RHOI
      EQUIS(5,5) = -(N+1.D0)/RADI-N/RADI
      EQUIS(5,6) = 1.D0
      EQUIS(6,1) = CONSTB*(N+1.D0)*RHOI/RADI
      EQUIS(6,3) = -N*(N+1.D0)*CONSTB*RHOI/RADI
      EQUIS(6,6) = (N-1.D0)/RADI-N/RADI

      EQUIST(1,1) = 1.D0/RADI-N/RADI
      EQUIST(1,2) = 1.D0/MUI
      EQUIST(2,1) = (N+2.D0)*(N-1.D0)*MUI/RADI/RADI
      EQUIST(2,2) = -3.D0/RADI-N/RADI

      DO 20 I = 1, 6
         DO 20 J= 1, 6
 20         EQUIS(I,J) = STEP*EQUIS(I,J)
      DO 21 I = 1, 2
         DO 21 J= 1, 2
 21         EQUIST(I,J) = STEP*EQUIST(I,J)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C        ********************************************
C        *                                          *
C        *     CALCULATION OF THE COEFFICENTS OF    *
C        *     THE EQUATIONS OF MOTION IN SOLID     *
C        *          RESULT IN 'EQUIS(6,6)'          *
C        *                                          *
C        ********************************************
C
      SUBROUTINE COEFIS1(EQUIS,EQUIST,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8     LAMI,MUI,LAMDAC
      DIMENSION  EQUIS(6,6), EQUIST(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
C
      CALL ZERO(EQUIS, 6,6)
      CALL ZERO(EQUIST,2,2)

      EQUIS(1,1) = -2.D0*LAMI/(LAMI+2.D0*MUI)/RADI+N/RADI
      EQUIS(1,2) = 1.D0/(LAMI+2.D0*MUI)
      EQUIS(1,3) = N*(N+1.D0)*LAMI/(LAMI+2.D0*MUI)/RADI
      EQUIS(2,1) = (-FREQ**2*RHOI*CONSTC*RADIS*RADI*RADI/GRAVIS
     &             -4.D0*CONSTC*RHOI*GRAVI*RADI+4.D0*MUI*(3.D0*LAMI
     &             +2.D0*MUI)/(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(2,2) = -4.D0*MUI/(LAMI+2.D0*MUI)/RADI+N/RADI
      EQUIS(2,3) = N*(N+1.D0)*(CONSTC*RHOI*GRAVI*RADI-2.D0*MUI
     &             *(3.D0*LAMI+2.D0*MUI)/(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(2,4) = N*(N+1.D0)/RADI
      EQUIS(2,5) = CONSTC*RHOI*(N+1.D0)/RADI
      EQUIS(2,6) = -CONSTC*RHOI
      EQUIS(3,1) = -1.D0/RADI
      EQUIS(3,3) = 1.D0/RADI+N/RADI
      EQUIS(3,4) = 1.D0/MUI
      EQUIS(4,1) = (CONSTC*RHOI*GRAVI*RADI-2.D0*MUI*(3.D0*LAMI+2.D0*MUI)
     &             /(LAMI+2.D0*MUI))/RADI/RADI
      EQUIS(4,2) = -LAMI/(LAMI+2.D0*MUI)/RADI
      EQUIS(4,3) = -FREQ**2*RHOI*CONSTC*RADIS/GRAVIS+2.D0*MUI
     &             /(LAMI+2.D0*MUI)
     &             /RADI/RADI*(LAMI*(2.D0*N*N+2.D0*N-1.D0)+2.D0*MUI
     &             *(N*N+N-1.D0))
      EQUIS(4,4) = -3.D0/RADI+N/RADI
      EQUIS(4,5) = -CONSTC*RHOI/RADI
      EQUIS(5,1) = CONSTB*RHOI
      EQUIS(5,5) = -(N+1.D0)/RADI+N/RADI
      EQUIS(5,6) = 1.D0
      EQUIS(6,1) = CONSTB*(N+1.D0)*RHOI/RADI
      EQUIS(6,3) = -N*(N+1.D0)*CONSTB*RHOI/RADI
      EQUIS(6,6) = (N-1.D0)/RADI+N/RADI

      EQUIST(1,1) = 1.D0/RADI+N/RADI
      EQUIST(1,2) = 1.D0/MUI
      EQUIST(2,1) = (N+2.D0)*(N-1.D0)*MUI/RADI/RADI
      EQUIST(2,2) = -3.D0/RADI+N/RADI

      DO 20 I = 1, 6
         DO 20 J= 1, 6
 20         EQUIS(I,J) = STEP*EQUIS(I,J)
      DO 21 I = 1, 2
         DO 21 J= 1, 2
 21         EQUIST(I,J) = STEP*EQUIST(I,J)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C        ********************************************
C        *                                          *
C        *     CALCULATION OF COEFFICENTS OF        *
C        *      EQUATIONS OF MOTION IN SOLID        *
C        *      RESULT IN 'EQUIS(2,2)' (N=0)        *
C        *                                          *
C        ********************************************
C
      SUBROUTINE COEFS0(EQUIS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT   REAL*8(A-H,O-Z)
      REAL*8     LAMI,MUI
      DIMENSION  EQUIS(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      CALL ZERO(EQUIS,2,2)

        EQUIS(1,1) = -2.D0*LAMI/(LAMI+2.D0*MUI)/RADI
        EQUIS(1,2) = 1.D0/(LAMI+2.D0*MUI)
        EQUIS(2,1) = (-4.D0*CONSTC*RHOI*GRAVI*RADI+4.D0*MUI*(3.D0*LAMI
     &               +2.D0*MUI)/(LAMI+2.D0*MUI))/RADI/RADI
        EQUIS(2,2) = -4.D0*MUI/(LAMI+2.D0*MUI)/RADI

        DO 20 I=1,2
           DO 20 J=1,2
 20           EQUIS(I,J) = STEP*EQUIS(I,J)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C       **************************************************
C       *                                                *
C       *      SOURCE FUNCTIONS OF DISLOCATIONS          *
C       *      AS INITIAL VALUES OF INTEGRATION          *
C       *             AT THE SOURCE POINT                *
C       *                                                *
C       **************************************************
C
      SUBROUTINE DISLO(ITYPE,N,INITS,INITST,SUMS,SUMST)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      ARGUMENTS: N     -- N DEGREE
C                 ITYPE -- FOUR TYPE DISLOCATIONS
C                 SUMS  -- A UNIT TENSOR
C                 INITS -- OUTPUT
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8    MUI,LAMI,INITS,INITST
      DIMENSION INITS(6,1),SUMS(6,6),INITST(2,1),SUMST(2,2)
      COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      CALL ZERO(INITST,2,1)
      CALL ZERO(INITS,6,1)

C  ITYPE1: X2=-1,Y1=1; ITYPE2: X2=-1,Y3=1; ITYPE3: X2=Y2=1; ITYP4: X3=Y3=1

      IF(ITYPE.EQ.1) THEN
         INITS(4,1) =-(2.D0*N+1.D0)/4.D0/PAI/N/(N+1.D0)/2.D0/RADI**3*MUI
         INITST(2,1)=-(2.D0*N+1.D0)/4.D0/PAI/N/(N+1.D0)/2.D0/RADI**3*MUI
      ENDIF

      IF(ITYPE.EQ.2) THEN
         INITS(3,1) = (2.D0*N+1.D0)/4.D0/PAI/N/(N+1.D0)/2.D0/RADI**2
         INITST(1,1)=-(2.D0*N+1.D0)/4.D0/PAI/N/(N+1.D0)/2.D0/RADI**2
      ENDIF

      IF(ITYPE.EQ.3) THEN
         INITS(1,1) = (2.D0*N+1.D0)/4.D0/PAI/RADI**2*
     &                LAMI/(LAMI+2.D0*MUI)
         INITS(2,1) = -(2.D0*N+1.D0)/2.D0/PAI/RADI**3
     &                *MUI*(3*LAMI+2*MUI)/(LAMI+2.D0*MUI)
         INITS(4,1) = (2.D0*N+1.D0)/4.D0/PAI/RADI**3
     &                *MUI*(3*LAMI+2*MUI)/(LAMI+2.D0*MUI)
      ENDIF

      IF(ITYPE.EQ.4) THEN
C       INITS(1,1) =  (2.D0*N+1.D0)/4.D0/PAI/RADI**2

C       ! unit mass internal loading new, 2023/01/31
        INITS(1,1) =  0.D0
          INITS(2,1) =  (2.D0*N+1.D0)/RADI**2*(CONSTA*GRAVI)
        INITS(6,1) = -(2.D0*N+1.D0)/RADI**2

      ENDIF

      CALL UNIT(SUMS, 6)
      CALL UNIT(SUMST,2)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C       **************************************************
C       *                                                *
C       *      SOURCE FUNCTIONS OF DISLOCATIONS          *
C       *      AS INITIAL VALUES OF INTEGRATION          *
C       *         AT THE SOURCE POINT (N=0)              *
C       *                                                *
C       **************************************************
C
      SUBROUTINE DISLO0(ITYPE,INITS,SUMS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT   REAL*8(A-H,O-Z)
      REAL*8     MUI,LAMI,INITS
      DIMENSION  INITS(2,1),SUMS(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      IF(ITYPE.EQ.3) THEN
         INITS(1,1) = 1.D0/4.D0/PAI/RADI**2*LAMI/(LAMI+2.D0*MUI)
         INITS(2,1) = -1.D0/2.D0/PAI/RADI**3
     &                *MUI*(3*LAMI+2*MUI)/(LAMI+2.D0*MUI)
      ENDIF
      IF(ITYPE.EQ.4) THEN
C         INITS(1,1) = 1.D0/4.D0/PAI/RADI**2
C         INITS(2,1) = 0.D0


C        ! unit mass internal loading new
         INITS(1,1) = 0.D0
         INITS(2,1) = 1.D0/RADI**2*(CONSTA*GRAVI)




      ENDIF

      CALL UNIT(SUMS,2)

      RETURN
      END

C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THIS SUBROUTINE IS USED FOR SOLVE LINEAR EQUATIONS, A GAUSS METHOD.
C
C  ARGUMENTS:
C     N      - THE SIZE OF THE DIMENSION OF THE EQUATIONS.
C     A(N,N) - THE COEFFICIENT MATRIX OF THE EQUATION
C     B(N)   - THE RIGHT SIDE TERM OF THE LINEAR EQUATIONS
C              IT IS AN INPUT FIRST, THEN IS USED TO STORE OUTPUT FINALY.
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE GS4(N,A,B,EP,KWJI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(N,N),B(N),M(100)

      DO 10 I=1,N
 10      M(I) = I
      DO 20 K=1,N
         P=0.D0
         DO 30 I=K,N
            DO 30 J=K,N
               IF(DABS(A(I,J)).LE.DABS(P)) GOTO 30
               P = A(I,J)
               IO = I
               JO = J
 30      CONTINUE
         IF(DABS(P)-EP) 200,200,300
 200     KWJI = 1
         RETURN
 300     IF(JO.EQ.K) GOTO 400
         DO 40 I = 1,N
            T = A(I,JO)
            A(I,JO) = A(I,K)
 40         A(I,K) = T
         J = M(K)
         M(K) =M(JO)
         M(JO) = J
 400     IF(IO.EQ.K) GOTO 500
         DO 50 J = K,N
            T = A(IO,J)
            A(IO,J) = A(K,J)
 50         A(K,J) = T
         T = B(IO)
         B(IO) = B(K)
         B(K) = T
 500     P = 1.D0/P
         IN = N-1
         IF(K.EQ.N) GOTO 600
         DO 60 J = K,IN
 60         A(K,J+1) = A(K,J+1)*P
 600     B(K) = B(K)*P
         IF(K.EQ.N) GOTO 20
         DO 70 I = K,IN
            DO 80 J = K,IN
 80            A(I+1,J+1)=A(I+1,J+1)-A(I+1,K)*A(K,J+1)
 70         B(I+1)=B(I+1)-A(I+1,K)*B(K)
 20   CONTINUE
      DO 90 I1 = 2,N
         I=N+1-I1
         DO 90 J=I,IN
 90         B(I)=B(I)-A(I,J+1)*B(J+1)
      DO 1 K = 1,N
         I = M(K)
 1       A(1,I) = B(K)
      DO 2 K =1,N
 2       B(K) = A(1,K)
      KWJI = 0
      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C  ********************************************************
C  *                                                      *
C  *             INITIAL VALUES FOR INTEGRATION           *
C  *               IN THE CENTER OF THE EARTH             *
C  *            FORMULA FROM 'SEISMIC SURFACE WAVES'      *
C  *              BY TAKEUCHI & SAITO (PP.243-245)        *
C  *                                                      *
C  *            ALF,BAT   -- WAVE VELOCITIES              *
C  *           INITS(6,3) -- INITIAL VALUES               *
C  *                 ****(N=0)****                        *
C  *                                                      *
C  ********************************************************
C
      SUBROUTINE INITI0(INITS,SUMS,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT   REAL*8(A-H,O-Z)
      REAL*8     MUI,LAMI,K,INITS,LAMDAC
      DIMENSION  INITS(2,1),SUMS(2,2)
      COMMON     /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON     /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON     /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
C
      ALF = DSQRT((LAMI+2.D0*MUI)/RHOI*GRAVIS*RADIS/CONSTC)
      BAT = DSQRT(MUI/RHOI*GRAVIS*RADIS/CONSTC)
      GAM = CONSTB*GRAVIS*RHOI/3.D0/RADIS
      K   = (((FREQ**2+4.D0*GAM)/ALF**2+(FREQ/BAT)**2
     &      +DSQRT(((FREQ/BAT)**2-(FREQ**2+4.D0*GAM)/ALF**2)**2
     &      +4.D0*N*(N+1.D0)*GAM**2/(ALF*BAT)**2))/2.D0)
      F   = BAT**2/GAM*(K-(FREQ/BAT)**2)
      H   = F-(N+1)
      X   = K*(RADIS*RADI)**2

      FIN = 1.D0
      DO 1 IJK = 1,100
         HKK = 1.D0
         HJJ = 1.D0
         DO 2 IIJ = 1,IJK
 2          HKK = HKK*IIJ
         DO 3 IIJ = 2*N+1+2*IJK, 2*N+1+1, -2
 3          HJJ = HJJ*IIJ
         FINS = FIN
         FIN  = FIN+(-X/2.D0)**IJK/HKK/HJJ
       IF (DABS(FIN-FINS).LT.1D-15) GOTO 4
 1     CONTINUE
 4     CONTINUE

       DAF = 2.D0*(2.D0*N+3)/X*(1.D0-FIN)

       NIJ  = N+1
       FIN1 = 1.D0
       DO 91 IJK = 1, 100
         HKK = 1.D0
         HJJ = 1.D0
         DO 92 IIJ = 1, IJK
 92          HKK = HKK*IIJ
         DO 93 IIJ= 2*NIJ+1+2*IJK, 2*NIJ+1+1, -2
 93          HJJ = HJJ*IIJ
         FINS = FIN1
         FIN1 = FIN1+(-X/2.D0)**IJK/HKK/HJJ
         IF (DABS(FIN1-FINS).LT.1D-15) GOTO 94
 91    CONTINUE
 94    CONTINUE

      INITS(1,1) = -1.D0/(2.D0*N+3.D0)*(N*H*DAF/2.D0+F*FIN1)
      INITS(2,1) = -(LAMI+2.D0*MUI)/RADI*F*FIN+MUI/RADI/(2.D0*N+3.D0)
     &             *(-N*(N-1.D0)*H*DAF+2.D0*(2.D0*F+N*(N+1.D0))*FIN1)

      CALL UNIT(SUMS,2)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C   *******************************************************
C   *                                                     *
C   *            INITIAL VALUES AT THE CENTER EARTH       *
C   *               FOR INTEGRATING INNER CORE            *
C   *           FORMULA FROM 'SEISMIC SURFACE WAVES'      *
C   *             BY TAKEUCHI & SAITO (PP.243-245)        *
C   *                                                     *
C   *            ALF,BAT    -- WAVE VELOCITIES            *
C   *            INITS(6,3) -- INITIAL VALUES             *
C   *                                                     *
C   *******************************************************
C
      SUBROUTINE INITIAL(INITS,INITST,SUMS,SUMST,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  THIS SUBROUTIN CALCULATE THE INITIAL SOLUTIONS (VALUES) AT THE CENTER
C  OF THE EARTH. ACTUALLY IT IS A VERY SMALL HOMOGENOUS SPHERE, SO THAT
C  THERE EXIST THREE INDEPENDENT ANALYTICAL SOLUTIONS, REF. THE FORMULAS
C  FROM 'SEISMIC SURFACE WAVES' BY TAKEUCHI & SAITO (PP.243-245).
C  THE HOMOGENOUS SPHERE CAN BE TAKEN SMALL ENOUGH SO THAT THE INTEGRATION
C  OVER THE REALISTIC EARTH MODEL CAN NOT BE AFFECTED AT ALL.
C
C  ARGUMENTS:
C     INITS(6,3) - OUTPUT, TO BE PUT THE THREE SETS OF INITIAL SOLUTIONS
C     SUMS(6,6)  - A UNIT TERSOR USED FOR INNER CORE INTEGRATION
C     N          - THE HARMONIC DEGREE 'N'.
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      IMPLICIT  INTEGER*8(I-N)
      REAL*8    MUI,LAMI,K,INITS,INITST,LAMDAC
      DIMENSION INITS(6,3),SUMS(6,6),INITST(2,1),SUMST(2,2)
      COMMON    /ABC/ CONSTA,CONSTB,CONSTC,FREQ,PAI
      COMMON    /NUM/ RHOC,LAMDAC,GRAVIS,RADIS,GNEWTN,RHOA,SURFMU,SM,SL
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      XN = N

      ALF = DSQRT((LAMI+2.D0*MUI)/RHOI*GRAVIS*RADIS/CONSTC)
      BAT = DSQRT(MUI/RHOI*GRAVIS*RADIS/CONSTC)
      GAM = CONSTB*GRAVIS*RHOI/3.D0/RADIS
      K   = (((FREQ**2+4.D0*GAM)/ALF**2+(FREQ/BAT)**2
     &      +DSQRT(((FREQ/BAT)**2-(FREQ**2+4.D0*GAM)/ALF**2)**2
     &      +4.D0*XN*(XN+1.D0)*GAM**2/(ALF*BAT)**2))/2.D0)
      HK  = K*2.D0
      F   = BAT**2/GAM*(K-(FREQ/BAT)**2)
      H   = F-(XN+1.D0)
      X   = K*(RADIS*RADI)**2

      FIN = 1.D0
      DO 1 IJK = 1,100
         HKK = 1.D0
         HJJ = 1.D0
         DO 2 IIJ = 1,IJK/2
 2          HKK = (-X/2.D0)*HKK/IIJ
         HJJ =HKK
         DO 3 IIJ = 2*N+1+2*IJK, 2*N+1+1, -2
 3          HJJ = HJJ/IIJ
         HKK =HJJ
         DO 72 IIJ = IJK/2+1,IJK
 72          HKK = (-X/2.D0)*HKK/IIJ
         FINS = FIN
         FIN  = FIN+HKK
       IF (DABS(FIN-FINS).LT.1.D-30) GOTO 4
 1     CONTINUE
 4     CONTINUE

       DAF = 2.D0*(2.D0*XN+3.D0)/X*(1.D0-FIN)

       NIJ  = N+1
       FIN1 = 1.D0
       DO 91 IJK = 1, 100
         HKK = 1.D0
         HJJ = 1.D0
         DO 92 IIJ = 1, IJK/2
 92          HKK = (-X/2.D0)*HKK/IIJ
             HJJ = HKK
         DO 93 IIJ= 2*NIJ+1+2*IJK, 2*NIJ+1+1, -2
 93          HJJ = HJJ/IIJ
             HKK = HJJ
         DO 792 IIJ = IJK/2+1, IJK
 792          HKK = (-X/2.D0)*HKK/IIJ
         FINS = FIN1
         FIN1 = FIN1+HKK
         IF (DABS(FIN1-FINS).LT.1.D-30) GOTO 94
 91    CONTINUE
 94    CONTINUE

      INITS(1,1) = XN
      INITS(2,1) = 2.D0*MUI*XN*(XN-1.D0)/RADI
      INITS(3,1) = 1.D0
      INITS(4,1) = 2.D0*MUI*(XN-1.D0)/RADI
      INITS(5,1) = (XN*GAM-FREQ**2)/GRAVIS*RADIS*RADI
      INITS(6,1) = (2.D0*XN+1.D0)*INITS(5,1)/RADI
     &             -3.D0*XN*GAM*RADIS/GRAVIS

      INITS(1,2) = -1.D0/(2.D0*XN+3.D0)*(XN*H*DAF/2.D0+F*FIN1)
      INITS(2,2) = -(LAMI+2.D0*MUI)/RADI*F*FIN+MUI/RADI/(2.D0*XN+3.D0)
     &             *(-XN*(XN-1.D0)*H*DAF+2.D0*(2.D0*F
     &             +XN*(XN+1.D0))*FIN1)
      INITS(3,2) = -1.D0/(2.D0*XN+3.D0)*(H*DAF/2.D0-FIN1)
      INITS(4,2) = MUI/RADI*(FIN-((XN-1.D0)*H*DAF+2.D0*(F+1.D0)*FIN1)
     &             /(2.D0*XN+3.D0))
      INITS(5,2) = RADIS*RADI/GRAVIS*((ALF**2*F-(XN+1.D0)*BAT**2)
     &             /(RADIS*RADI)**2-3.D0*GAM*F*DAF/2.D0/(2.D0*XN+3.D0))
      INITS(6,2) = (2.D0*XN+1.D0)/RADI*INITS(5,2)+3.D0*XN*GAM*H*RADIS
     &             /2.D0/(2.D0*XN+3.D0)/GRAVIS*DAF

      K = (((FREQ**2+4.D0*GAM)/ALF**2+(FREQ/BAT)**2)**2-(((FREQ/BAT)**2
     &    -(FREQ**2+4.D0*GAM)/ALF**2)**2+4.D0*XN*(XN+1.D0)*GAM**2
     &    /(ALF*BAT)**2))/2.D0/HK
      F = BAT**2/GAM*(K-(FREQ/BAT)**2)
      H = F-(XN+1.D0)
      X = K*RADIS*RADI*RADIS*RADI


      FIN = 1.D0
      DO 981 IJK = 1, 100
         HKK = 1.D0
         HJJ = 1.D0
         DO 982 IIJ = 1, IJK/2
 982          HKK = (-X/2.D0)*HKK/IIJ
              HJJ = HKK
         DO 983 IIJ = 2*N+1+2*IJK, 2*N+1+1, -2
 983          HJJ = HJJ/IIJ
              HKK = HJJ
         DO 987 IIJ = IJK/2+1, IJK
 987          HKK = (-X/2.D0)*HKK/IIJ
         FINS = FIN
         FIN  = FIN+HKK
         IF (DABS(FIN-FINS).LT.1.D-30) GOTO 984
 981     CONTINUE
 984  CONTINUE

       DAF = 2.D0*(2.D0*XN+3.D0)/X*(1.D0-FIN)

      NIJ = N+1
      FIN1 = 1.D0
      DO 991 IJK = 1, 100
         HKK = 1.D0
         HJJ = 1.D0
         DO 992 IIJ = 1, IJK/2
 992          HKK = (-X/2.D0)*HKK/IIJ
              HJJ = HKK
         DO 993 IIJ = 2*NIJ+1+2*IJK, 2*NIJ+1+1, -2
 993          HJJ = HJJ/IIJ
              HKK = HJJ
         DO 997 IIJ = IJK/2+1, IJK
 997          HKK = (-X/2.D0)*HKK/IIJ
         FINS = FIN1
         FIN1 = FIN1+HKK
         IF(DABS(FIN1-FINS).LT.1.D-30) GOTO 994
 991     CONTINUE
 994  CONTINUE


      INITS(1,3) = -1.D0/(2.D0*XN+3.D0)*(XN*H*DAF/2.D0+F*FIN1)
      INITS(2,3) = -(LAMI+2.D0*MUI)/RADI*F*FIN+MUI/RADI/(2.D0*XN+3.D0)
     &             *(-XN*(XN-1.D0)*H*DAF+2.D0*(2.D0*F+XN*(XN+1.D0))
     &             *FIN1)
      INITS(3,3) = -1.D0/(2.D0*XN+3.D0)*(H*DAF/2.D0-FIN1)
      INITS(4,3) = MUI/RADI*(FIN-((XN-1.D0)*H*DAF+2.D0*(F+1.D0)*FIN1)
     &             /(2.D0*XN+3.D0))
      INITS(5,3) = RADIS*RADI/GRAVIS*((ALF**2*F-(XN+1)*BAT**2)
     &             /(RADIS*RADI)**2-3.D0*GAM*F*DAF/2.D0/(2.D0*XN+3.D0))
      INITS(6,3) = (2.D0*XN+1.D0)/RADI*INITS(5,3)+3.D0*XN*GAM*H*RADIS
     &             /2.D0/(2.D0*XN+3.D0)/GRAVIS*DAF

      INITST(1,1) = 1.D0
C      INITST(2,1) = MUI * (XN - 1.D0) / RADI
      INITST(2,1) = 0.D0

      CALL UNIT(SUMS, 6)
      CALL UNIT(SUMST,2)

C 321  FORMAT(I5,5D13.6)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C    ****************************************************
C    *  MULTIPLICATION OF MATRIX IN RUNGE-KUTTA METHOD  *
C    ****************************************************
C
      SUBROUTINE KZ1(N,WORK7,EQUIL,WORK5,SUML)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  WORK7(N,N),EQUIL(N,N),WORK5(N,N),SUML(N,N)

      DO 10 I = 1, N
         DO 10 J = 1, N
            WORK7(I,J) = 0.D0
            DO 30 IJ = 1, N
               WORK7(I,J) = WORK7(I,J) + EQUIL(I,IJ) * WORK5(IJ,J)/2.D0
 30         CONTINUE
 10         WORK7(I,J) = WORK7(I,J) + SUML(I,J)

      DO 40 I = 1, N
         DO 40 J = 1, N
 40         WORK5(I,J) = WORK7(I,J)

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C    *********************************************
C    *     CALCULATING  A(N,N)=2*B(N,N)*C(N,N)   *
C    *********************************************
C
      SUBROUTINE KZ3(N,RONGKL,EQUIL,WORK7)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  RONGKL(N,N),EQUIL(N,N),WORK7(N,N)

      DO 10 I = 1, N
         DO 10 J = 1, N
            DO 10 IJ = 1, N
 10            RONGKL(I,J) = RONGKL(I,J) + 2.D0*EQUIL(I,IJ)*WORK7(IJ,J)
      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C    ******************************************************
C    *                                                    *
C    *    MULTIPLICATION OF MATRIX: A(L,N)=B(L,M)*C(M,N)  *
C    *                                                    *
C    ******************************************************
C
      SUBROUTINE MATRIX(L,M,N,B,C,A)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  B(L,M), C(M,N), A(L,N)

      DO 10 I = 1, L
         DO 10 J = 1, N
            A(I,J) = 0.D0
            DO 10 IJ = 1, M
 10            A(I,J) = A(I,J) + B(I,IJ) * C(IJ,J)
      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C    ******************************************************
C    *                                                    *
C    *    MULTIPLICATION OF MATRIX: A(L,N)=B(L,M)*C(M,N)  *
C    *                                                    *
C    ******************************************************
C
      SUBROUTINE MATRIX1(N,M,B,A)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  B(N,M), A(N,M)

      DO 10 I = 1, N
         DO 10 J = 1, M
 10            A(I,J) = B(I,J)
      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C  *****************************************************************
C  *                                                               *
C  *  INTEGRATING SOLID PARTS (INNER CORE OR MENTLE)               *
C  *  RUNGE-KUTTA METHOD  :    Z = Z + (K1+2K2+2K3+K4)/6           *
C  *  SUBROUTINES:  COEFS0 -- COEFFICIENT OF EQUATIONS IN SOLID    *
C  *                KZ1    -- MULTIPLYING MATRIX IN THE METHOD     *
C  *                MATRIX -- A(L,N)=B(L,M)*C(M,N)                 *
C  *                KZ3    -- A(N,N)=2*B(N,N)*C(N,N)               *
C  *                DATASET-- PREPARY DATA (I1234=1,2,4)           *
C  *     ****(N=0)****                                             *
C  *                                                               *
C  *****************************************************************
C
      SUBROUTINE SOLID0(SUMS,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8    LAMDA,MU,LAMI,MUI
      DIMENSION SUMS(2,2),WORK3(2,2),EQUIS(2,2),RONGKS(2,2),WORK6(2,2)
     &          ,RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      DO 90 I1234=1,4

      IF (I1234.NE.1) GOTO 20
         DO 70 I=1,2
            DO 70 J=1,2
 70            WORK3(I,J) = SUMS(I,J)
      CALL COEFS0(EQUIS)
      CALL KZ1(2,WORK6,EQUIS,WORK3,SUMS)
      CALL MATRIX(2,2,2,EQUIS,SUMS,RONGKS)
      GOTO 90

 20   IF (I1234.NE.2) GOTO 100
      CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
      CALL COEFS0(EQUIS)
      CALL KZ3(2,RONGKS,EQUIS,WORK6)
      CALL KZ1(2,WORK6,EQUIS,WORK3,SUMS)
      GOTO 90

 100  IF (I1234.NE.3) GOTO 110
      CALL KZ3(2,RONGKS,EQUIS,WORK6)
      DO 120 I=1,2
         DO 120 J=1,2
            WORK6(I,J)=0.D0
            DO 130 KI=1,2
 130           WORK6(I,J)=WORK6(I,J)+EQUIS(I,KI)*WORK3(KI,J)
 120           WORK6(I,J)=WORK6(I,J)+SUMS(I,J)
      GOTO 90

110   CONTINUE
      CALL DATASET(4,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
      CALL COEFS0(EQUIS)
      DO 180 I=1,2
         DO 180 J=1,2
            DO 180 KI=1,2
 180        RONGKS(I,J) = RONGKS(I,J)+EQUIS(I,KI)*WORK6(KI,J)
 90   CONTINUE

      DO 190 I=1,2
         DO 190 J=1,2
 190        SUMS(I,J) = SUMS(I,J)+RONGKS(I,J)/6.D0

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C  ****************************************************************
C  *                                                              *
C  *  INTEGRATING SOLID PARTS (INNER CORE OR MENTLE)              *
C  *  RUNGE-KUTTA METHOD  :    Z = Z + (K1+2K2+2K3+K4)/6          *
C  *  SUBROUTINES: COEFIS -- COEFFICIENT OF EQUATIONS IN SOLID    *
C  *             KZ1    -- MULTIPLYING MATRIX IN THE METHOD       *
C  *             MATRIX -- A(L,N)=B(L,M)*C(M,N)                   *
C  *             KZ3    -- A(N,N)=2*B(N,N)*C(N,N)                 *
C  *             DATASET-- PREPARY DATA (I1234=1,2,4)             *
C  *                                                              *
C  ****************************************************************
C
      SUBROUTINE SOLIDST(N,SUMS,SUMST
     &       ,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTEGRATION IN SOLID PART. SPHEROIDAL AND TOROIDAL PARTS TOGETHER.
C  WITHIN THE SUBROUTINE, THE EARTH MODEL VALUES AT THE CALCULATING
C  POINT ARE CHANGED (RESIGNED THROUGH
C      CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C  WHILE THE VALUES ARE PASSED ON GLOBALY BY THE COMMON CENTANSE
C      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
C  THEREFORE, TO CALL THIS SUBROUTINE, SPHEROIDAL AND TOROIDAL PARTS
C  SHOULD COME TOGETHER, TO AVOID MIXING UP.
C
C      ARGUMENTS: N     -- N DEGREE
C                 MSTEP -- MAXIMUM LAYER OF ADOPTED EARTH MODEL
C                 GRAVITY, RHO, LAMDA, MU, RADIUS -- EARTH MODEL
C                 EQUIS -- COEFFICIENTS OF EQUATION BY SUBROUTINE 'COEFIS'
C                 WORK3,WORK6,RONGKS -- FOR WORK
C
C                 EQUIST,WORK3T,WORK6T,RONGKST -- FOR TOROIDAL CASE
C
C      OUTPUT IN: SUMS  -- SPHEROIDAL RESULTS
C                 SUMST -- TOROIDAL RESULTS
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8    LAMDA,MU,LAMI,MUI
      DIMENSION SUMS(6,6),WORK3(6,6),EQUIS(6,6),RONGKS(6,6),WORK6(6,6)
      DIMENSION SUMST(2,2),WORK3T(2,2),EQUIST(2,2),RONGKST(2,2)
     &          ,WORK6T(2,2)
      DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      DO 90 I1234=1,4

      IF (I1234.NE.1) GOTO 20

         DO 70 I=1,6
            DO 70 J=1,6
 70            WORK3(I,J) = SUMS(I,J)
         DO 71 I=1,2
            DO 71 J=1,2
 71            WORK3T(I,J) = SUMST(I,J)
      CALL COEFIS(EQUIS,EQUIST,N)
      CALL KZ1(6,WORK6,EQUIS,WORK3,SUMS)
      CALL MATRIX(6,6,6,EQUIS,SUMS,RONGKS)
      CALL KZ1(2,WORK6T,EQUIST,WORK3T,SUMST)
      CALL MATRIX(2,2,2,EQUIST,SUMST,RONGKST)

      GOTO 90

 20   IF (I1234.NE.2) GOTO 100

      CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)

      CALL COEFIS(EQUIS,EQUIST,N)
      CALL KZ3(6,RONGKS,EQUIS,WORK6)
      CALL KZ1(6,WORK6,EQUIS,WORK3,SUMS)
      CALL KZ3(2,RONGKST,EQUIST,WORK6T)
      CALL KZ1(2,WORK6T,EQUIST,WORK3T,SUMST)

      GOTO 90

 100  IF (I1234.NE.3) GOTO 110

      CALL KZ3(6,RONGKS,EQUIS,WORK6)
      DO 120 I=1,6
         DO 120 J=1,6
            WORK6(I,J)=0.D0
            DO 130 KI=1,6
 130           WORK6(I,J)=WORK6(I,J)+EQUIS(I,KI)*WORK3(KI,J)
 120           WORK6(I,J)=WORK6(I,J)+SUMS(I,J)

      CALL KZ3(2,RONGKST,EQUIST,WORK6T)
      DO 121 I=1,2
         DO 121 J=1,2
            WORK6T(I,J)=0.D0
            DO 131 KI=1,2
 131           WORK6T(I,J)=WORK6T(I,J)+EQUIST(I,KI)*WORK3T(KI,J)
 121           WORK6T(I,J)=WORK6T(I,J)+SUMST(I,J)

      GOTO 90

110   CONTINUE

      CALL DATASET(4,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)

      CALL COEFIS(EQUIS,EQUIST,N)
      DO 180 I=1,6
         DO 180 J=1,6
            DO 180 KI=1,6
 180        RONGKS(I,J) = RONGKS(I,J)+EQUIS(I,KI)*WORK6(KI,J)
      DO 181 I=1,2
         DO 181 J=1,2
            DO 181 KI=1,2
 181        RONGKST(I,J) = RONGKST(I,J)+EQUIST(I,KI)*WORK6T(KI,J)

 90   CONTINUE

      DO 190 I=1,6
         DO 190 J=1,6
 190        SUMS(I,J) = SUMS(I,J)+RONGKS(I,J)/6.D0

      DO 191 I=1,2
         DO 191 J=1,2
 191        SUMST(I,J) = SUMST(I,J)+RONGKST(I,J)/6.D0

       RETURN
      END

C  ################### SUBROUTINE ##################################
C
C  ****************************************************************
C  *                                                              *
C  *  INTEGRATING SOLID PARTS (INNER CORE OR MENTLE)              *
C  *  RUNGE-KUTTA METHOD  :    Z = Z + (K1+2K2+2K3+K4)/6          *
C  *  SUBROUTINES: COEFIS -- COEFFICIENT OF EQUATIONS IN SOLID    *
C  *             KZ1    -- MULTIPLYING MATRIX IN THE METHOD       *
C  *             MATRIX -- A(L,N)=B(L,M)*C(M,N)                   *
C  *             KZ3    -- A(N,N)=2*B(N,N)*C(N,N)                 *
C  *             DATASET-- PREPARY DATA (I1234=1,2,4)             *
C  *                                                              *
C  ****************************************************************
C
      SUBROUTINE SOLIDST1(N,SUMS,SUMST
     &       ,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INTEGRATION IN SOLID PART. SPHEROIDAL AND TOROIDAL PARTS TOGETHER.
C  WITHIN THE SUBROUTINE, THE EARTH MODEL VALUES AT THE CALCULATING
C  POINT ARE CHANGED (RESIGNED THROUGH
C      CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C  WHILE THE VALUES ARE PASSED ON GLOBALY BY THE COMMON CENTANSE
C      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI
C  THEREFORE, TO CALL THIS SUBROUTINE, SPHEROIDAL AND TOROIDAL PARTS
C  SHOULD COME TOGETHER, TO AVOID MIXING UP.
C
C      ARGUMENTS: N     -- N DEGREE
C                 MSTEP -- MAXIMUM LAYER OF ADOPTED EARTH MODEL
C                 GRAVITY, RHO, LAMDA, MU, RADIUS -- EARTH MODEL
C                 EQUIS -- COEFFICIENTS OF EQUATION BY SUBROUTINE 'COEFIS'
C                 WORK3,WORK6,RONGKS -- FOR WORK
C
C                 EQUIST,WORK3T,WORK6T,RONGKST -- FOR TOROIDAL CASE
C
C      OUTPUT IN: SUMS  -- SPHEROIDAL RESULTS
C                 SUMST -- TOROIDAL RESULTS
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8    LAMDA,MU,LAMI,MUI
      DIMENSION SUMS(6,6),WORK3(6,6),EQUIS(6,6),RONGKS(6,6),WORK6(6,6)
      DIMENSION SUMST(2,2),WORK3T(2,2),EQUIST(2,2),RONGKST(2,2)
     &          ,WORK6T(2,2)
      DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      DO 90 I1234=1,4

      IF (I1234.NE.1) GOTO 20

         DO 70 I=1,6
            DO 70 J=1,6
 70            WORK3(I,J) = SUMS(I,J)
         DO 71 I=1,2
            DO 71 J=1,2
 71            WORK3T(I,J) = SUMST(I,J)
      CALL COEFIS1(EQUIS,EQUIST,N)
      CALL KZ1(6,WORK6,EQUIS,WORK3,SUMS)
      CALL MATRIX(6,6,6,EQUIS,SUMS,RONGKS)
      CALL KZ1(2,WORK6T,EQUIST,WORK3T,SUMST)
      CALL MATRIX(2,2,2,EQUIST,SUMST,RONGKST)

      GOTO 90

 20   IF (I1234.NE.2) GOTO 100

      CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)

      CALL COEFIS1(EQUIS,EQUIST,N)
      CALL KZ3(6,RONGKS,EQUIS,WORK6)
      CALL KZ1(6,WORK6,EQUIS,WORK3,SUMS)
      CALL KZ3(2,RONGKST,EQUIST,WORK6T)
      CALL KZ1(2,WORK6T,EQUIST,WORK3T,SUMST)

      GOTO 90

 100  IF (I1234.NE.3) GOTO 110

      CALL KZ3(6,RONGKS,EQUIS,WORK6)
      DO 120 I=1,6
         DO 120 J=1,6
            WORK6(I,J)=0.D0
            DO 130 KI=1,6
 130           WORK6(I,J)=WORK6(I,J)+EQUIS(I,KI)*WORK3(KI,J)
 120           WORK6(I,J)=WORK6(I,J)+SUMS(I,J)

      CALL KZ3(2,RONGKST,EQUIST,WORK6T)
      DO 121 I=1,2
         DO 121 J=1,2
            WORK6T(I,J)=0.D0
            DO 131 KI=1,2
 131           WORK6T(I,J)=WORK6T(I,J)+EQUIST(I,KI)*WORK3T(KI,J)
 121           WORK6T(I,J)=WORK6T(I,J)+SUMST(I,J)

      GOTO 90

110   CONTINUE

      CALL DATASET(4,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)

      CALL COEFIS1(EQUIS,EQUIST,N)
      DO 180 I=1,6
         DO 180 J=1,6
            DO 180 KI=1,6
 180        RONGKS(I,J) = RONGKS(I,J)+EQUIS(I,KI)*WORK6(KI,J)
      DO 181 I=1,2
         DO 181 J=1,2
            DO 181 KI=1,2
 181        RONGKST(I,J) = RONGKST(I,J)+EQUIST(I,KI)*WORK6T(KI,J)

 90   CONTINUE

      DO 190 I=1,6
         DO 190 J=1,6
 190        SUMS(I,J) = SUMS(I,J)+RONGKS(I,J)/6.D0

      DO 191 I=1,2
         DO 191 J=1,2
 191        SUMST(I,J) = SUMST(I,J)+RONGKST(I,J)/6.D0

       RETURN
      END

C  ################### SUBROUTINE ##################################
C
C  ****************************************************************
C  *                                                              *
C  * INTEGRATING  LIQUID OUT CORE                                 *
C  * RUNGE-KUTTA METHOD  : Z = Z + (K1+2K2+2K3+K4)/6              *
C  * SUBROUTINES: COEFIL -- COEFFICIENT OF EQUATIONS IN LIQUID    *
C  *              KZ1    -- MULTIPLYING MATRIX IN THE METHOD      *
C  *              MATRIX -- A(L,N)=B(L,M)*C(M,N)                  *
C  *              KZ3    -- A(N,N)=2*B(N,N)*C(N,N)                *
C  *              DATASET-- PREPARE DATA (I1234=1,2,4)            *
C  *                                                              *
C  ****************************************************************
C
      SUBROUTINE LIQUID(N,SUML
     &                  ,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C      ARGUMENTS: N     -- N DEGREE
C                 MSTEP -- MAXIMUM LAYER OF ADOPTED EARTH MODEL
C                 GRAVITY, RHO, LAMDA, MU, RADIUS -- EARTH MODEL
C                 EQUIL -- COEFFICIENTS OF EQUATION IN LIQUID PART
C                          BY SUBROUTINE 'COEFIS'
C                 WORK5,WORK7,RONGKL -- FOR WORK
C      OUTPUT IN: SUML
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      REAL*8    LAMI,MUI,LAMDA,MU
      DIMENSION SUML(2,2),WORK5(2,2),EQUIL(2,2),WORK7(2,2),RONGKL(2,2)
      DIMENSION RADIUS(MSTEP),RHO(MSTEP),MU(MSTEP),LAMDA(MSTEP)
     &          ,GRAVITY(MSTEP)
      COMMON    /MDI/ RADI,STEP,RHOI,LAMI,MUI,GRAVI

      DO 120 I1234 = 1, 4

         IF (I1234.NE.1) GOTO 90
           DO 100 I=1,2
              DO 100 J=1,2
 100             WORK5(I,J) = SUML(I,J)
           CALL COEFIL(EQUIL,N)
           CALL KZ1(2,WORK7,EQUIL,WORK5,SUML)
           CALL MATRIX(2,2,2,EQUIL,SUML,RONGKL)
           GOTO 120

 90      IF (I1234.NE.2) GOTO 130
           CALL DATASET(2,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
           CALL COEFIL(EQUIL,N)
           CALL KZ3(2,RONGKL,EQUIL,WORK7)
           CALL KZ1(2,WORK7,EQUIL,WORK5,SUML)
           GOTO 120

 130     IF (I1234.NE.3) GOTO 140
           CALL KZ3(2,RONGKL,EQUIL,WORK7)
           DO 150 I=1,2
              DO 150 J=1,2
                 WORK7(I,J) = 0.D0
                 DO 160 IJ=1,2
 160                WORK7(I,J) = WORK7(I,J)+EQUIL(I,IJ)*WORK5(IJ,J)
 150             WORK7(I,J) = WORK7(I,J)+SUML(I,J)
           GOTO 120

 140     CONTINUE
           CALL DATASET(4,MSTEP,GRAVITY,RHO,LAMDA,MU,RADIUS)
           CALL COEFIL(EQUIL,N)
           DO 170 I=1,2
              DO 170 J=1,2
                 DO 170 IJ=1,2
 170                RONGKL(I,J) = RONGKL(I,J)+EQUIL(I,IJ)*WORK7(IJ,J)
 120     CONTINUE

      DO 180 I=1,2
         DO 180 J=1,2
 180        SUML(I,J) = SUML(I,J)+RONGKL(I,J)/6.D0

      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C     ***********************************************
C     *                                             *
C     *       SET  AN UNIT TENSOR TO A(N,N)         *
C     *                                             *
C     ***********************************************
C
      SUBROUTINE UNIT(A,N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(N,N)
      DO 1 I = 1, N
         DO 1 J = 1, N
            IF(I.EQ.J) A(I,J) = 1.D0
 1          IF(I.NE.J) A(I,J) = 0.D0
      RETURN
      END

C  ################### SUBROUTINE ##################################
C
C     ***********************************************
C     *                                             *
C     *       SET   0   TO  DIMENSION  A(N,M)       *
C     *                                             *
C     ***********************************************
C
      SUBROUTINE ZERO(A,N,M)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION  A(N,M)
      DO 1 I = 1, N
         DO 1 J = 1, M
 1          A(I,J) = 0.D0
      RETURN
      END
