      SUBROUTINE GREEN0KM(DEP,DEFORM,NUP,SL,SM)
!     DEP  = depth of source
!     DEFORM = depth of obs

      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      REAL*8    LOVEH,LOVEL,LOVEK,LOVET,LOVEH1,LOVEL1,LOVEK1,LOVET1
      PARAMETER ( MAXN=100005 )
      PARAMETER ( MLOVE=100005 )
      PARAMETER ( MANGLE=151 )
      PARAMETER(PAI=3.14159265358979d0)
      DIMENSION NI(MLOVE),LOVEH(MLOVE),LOVEL(MLOVE),LOVEK(MLOVE)
     &,LOVET(MLOVE),LOVEH1(MLOVE),LOVEL1(MLOVE)
     &,LOVEK1(MLOVE),LOVET1(MLOVE),ANGLE(MANGLE)
     &,YY1(MAXN),YY2(MAXN),YY3(MAXN),YY6(MAXN),YY7(MAXN)
     &,YY9(MAXN),YY10(MAXN),YY11(MAXN),YY12(MAXN)
     &,HL(4,2),disp2(2)
           common /info_source/ grs

      SL  =  SL / 13080.0
      SM  =  SM / 13080.0
      RADIS = 6371.0
          tt = 1d0-DEP/RADIS

      OPEN (41, FILE = 'DSLOV12.DAT',ACTION='READ')
      OPEN (42, FILE = 'DSLOV32.DAT',ACTION='READ')
      OPEN (43, FILE = 'DSLOV220.DAT',ACTION='READ')
      OPEN (44, FILE = 'DSLOV33.DAT',ACTION='READ')
      OPEN (51, FILE = 'GRNFN12.DAT')
      OPEN (52, FILE = 'GRNFN32.DAT')
      OPEN (53, FILE = 'GRNFN220.DAT')
      OPEN (54, FILE = 'GRNFN33.DAT')

C!      DATA ANGLE/
C!     &0.001,    0.002,    0.003,    0.004,    0.005,
C!     &0.006,    0.007,    0.008,    0.009,    0.010,
C!     &0.020,    0.030,    0.040,    0.050,    0.060,
C!     &0.070,    0.080,    0.090,    0.100,    0.110,
C!     &0.120,    0.130,    0.140,    0.150,    0.160,
C!     &0.170,    0.180,    0.190,    0.200,    0.210,
C!     &0.220,    0.230,    0.240,    0.250,    0.260,
C!     &0.270,    0.280,    0.290,    0.300,    0.310,
C!     &0.320,    0.330,    0.340,    0.350,    0.360,
C!     &0.370,    0.380,    0.390,    0.400,    0.410,
C!     &0.420,    0.430,    0.440,    0.450,    0.460,
C!     &0.470,    0.480,    0.490,    0.500,    0.550,
C!     &0.600,    0.650,    0.700,    0.750,    0.800,
C!     &0.850,    0.900,    0.950,    1.000,    1.100,
C!     &1.200,    1.300,    1.400,    1.500,    1.600,
C!     &2.000,    2.500,    3.000,    4.000,    5.000,
C!     &6.000,    7.000,    8.000,    9.000,   10.000,
C!     &12.000,   16.000,   20.000,   25.000,   30.000,
C!     &40.000,   50.000,   60.000,   70.000,   80.000,
C!     &90.000,  100.000,  110.000,  120.000,  130.000,
C!     &140.000,  150.000,  160.000,  170.000,  179.900 /

       DATA ANGLE/
     & 1d-5, 2d-5, 3d-5, 4d-5, 5d-5, 6d-5, 8d-5, 9d-5,
     & 0.0001d0,  0.0002d0,  0.0003d0,  0.0004d0, 0.0005d0,
     & 0.0006d0,  0.0007d0,  0.0008d0,  0.0009d0, 0.00095d0,
     & 0.001d0,   0.0015d0,  0.0018d0,  0.002d0,  0.003d0, 0.005d0,
     & 0.006d0,   0.007d0,   0.008d0,   0.009d0,  0.010d0, 0.015d0,
     & 0.020d0,   0.025d0,   0.030d0,   0.040d0,  0.050d0, 0.060d0,
     & 0.070d0,   0.080d0,   0.090d0,   0.100d0,  0.110d0,
     & 0.120d0,   0.130d0,   0.140d0,   0.150d0,  0.160d0,
     & 0.170d0,   0.180d0,   0.190d0,   0.200d0,  0.210d0,
     & 0.220d0,   0.230d0,   0.240d0,   0.250d0,  0.260d0,
     & 0.270d0,   0.280d0,   0.290d0,   0.300d0,  0.310d0,
     & 0.320d0,   0.330d0,   0.340d0,   0.350d0,  0.360d0,
     & 0.370d0,   0.380d0,   0.390d0,   0.400d0,  0.410d0,
     & 0.420d0,   0.430d0,   0.440d0,   0.450d0,  0.460d0,
     & 0.470d0,   0.480d0,   0.490d0,   0.500d0,  0.550d0,
     & 0.600d0,   0.650d0,   0.700d0,   0.750d0,  0.800d0,
     & 0.850d0,   0.900d0,   0.950d0,   1.000d0,  1.100d0,
     & 1.200d0,   1.300d0,   1.400d0,   1.500d0,  1.600d0,
     & 2.000d0,   2.500d0,   3.000d0,   4.000d0,  5.000d0,
     & 6.000d0,   7.000d0,   8.000d0,   9.000d0, 10.000d0,
     &  12.000d0,  16.000d0,  20.000d0,  25.000d0, 30.000d0,
     &  40.000d0,  50.000d0,  60.000d0,  70.000d0, 80.000d0,
     &  90.000d0, 100.000d0, 110.000d0, 120.000d0,  130.000d0,
     & 140.000d0, 150.000d0, 160.000d0, 170.000d0,  171.000d0,
     & 172.000d0, 173.000d0, 174.000d0, 175.000d0,  176.000d0,
     & 177.000d0, 178.000d0, 179.000d0, 179.100d0,  179.200d0,
     & 179.300d0, 179.400d0, 179.500d0, 179.600d0,  179.700d0,
     & 179.800d0, 179.900d0, 179.950d0, 179.980d0,  179.990d0,
     & 179.992d0, 179.994d0, 179.996d0, 179.998d0,  179.999d0/



      DO 399 ITYPE=4,4
        IF (ITYPE .EQ. 1) I12 = 2
        IF (ITYPE .EQ. 2) I12 = 1
        IF (ITYPE .EQ. 3) I12 = 1
        IF (ITYPE .EQ. 4) I12 = 1
        IF (ITYPE .EQ. 3) READ(43,*) NL,HN0,HN01,SL3,SL4,SL5,SL6
        IF (ITYPE .EQ. 4) READ(44,*) NL,HN0,HN01,SL3,SL4,SL5,SL6

        DO 543 I = I12, NUP
          IF (ITYPE.EQ.1)
     &    READ(41,*) NI(I), LOVEH(I),LOVEH1(I),LOVEL(I),LOVEL1(I)
     &    ,LOVEK(I),LOVEK1(I),LOVET(I),LOVET1(I)
          IF (ITYPE.EQ.2)
     &    READ(42,*) NI(I), LOVEH(I),LOVEH1(I),LOVEL(I),LOVEL1(I)
     &    ,LOVEK(I),LOVEK1(I),LOVET(I),LOVET1(I)
          IF (ITYPE.EQ.3)
     &    READ(43,*) NI(I), LOVEH(I),LOVEH1(I),LOVEL(I),LOVEL1(I)
     &    ,LOVEK(I),LOVEK1(I)
          IF (ITYPE.EQ.4)
     &    READ(44,*) NI(I), LOVEH(I),LOVEH1(I),LOVEL(I),LOVEL1(I)
     &    ,LOVEK(I),LOVEK1(I)
 543    CONTINUE

C ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C ! Asymptotic solution of y-numbers

C        IF (ITYPE .EQ. 4) THEN
C          call asymptotic_love(
C     &    SL*13080*1d8,SM*13080*1d8,(6371-DEP)*1D3,6371D3,HL)
C
C           ! write(*,'(4(ES16.6))') HL(:,1)
C         ! write(*,'(4(ES16.6))') HL(:,2)
C
C          ASY_H = 0d0
C         ASY_L = 0d0
C
C          ! 从2阶项开始, 减去高阶渐近级数
C         !write(*,*)'NUP = ',NUP
C         !write(*,*) '0.1^5000 = ', (1d0-DEP/6371d0)**NUP
C         !write(*,*) 'delta = rs/ra',(1d0-DEP/6371d0)
C
C          DO I = 2, NUP
C           ! write(*,*) 'delta = (rs/ra)^n',I, (1d0-DEP/6371d0)**NI(I)
C            ! write(*,*)I, NI(I)
C
C            ! 单位力高阶yi
C            ASY_H = HL(1,1)*NI(I)+HL(2,1)+HL(3,1)/
C     &              dble(NI(I))+HL(4,1)/(dble(NI(I))*dble(NI(I))-1d0)
C
C            ASY_L = HL(1,2)*NI(I)+HL(2,2)+HL(3,2)/
C     &              dble(NI(I))+HL(4,2)/(dble(NI(I))*dble(NI(I))-1d0)
C
C            ! 单位质量的yi
C            ASY_H = ASY_H * ((1d0-DEP/6371d0)**NI(I)) * grs
C            ASY_L = ASY_L * ((1d0-DEP/6371d0)**NI(I)) * grs
C
C            LOVEH(I) = LOVEH(I) - ASY_H
C            LOVEL(I) = LOVEL(I) - ASY_L
C
C            ! WRITE(123,'(2(2X,I10),4(2X,ES16.5))') I,NI(I),LOVEH(I),ASY_H,LOVEL(I),ASY_L
C
C          END DO
C        END IF
C

C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  CALCULATING GREEN'S FUNCTIONS FOR EACH ANGLE:
C
        DO 90 IANGLE = 1, MANGLE
          ANGLEI = ANGLE(IANGLE)

          U    = 0.D0
          VCT  = 0.D0
          VLM  = 0.D0
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  INITIAL VALUES FOR LEGENDRE FUNCTIONS:
C
          CT   = ANGLEI*PAI/180.D0
          D    = DCOS(CT/30.D0)
          X    = DCOS(CT)
          Y    = DSIN(CT)
          C2   = DCOS(CT/2.D0)
          S2   = DSIN(CT/2.D0)
          CO2  = DCOS(CT*2.D0)

          P  = 0.D0
          PA = 1.D0
          PB = X
          P1 = 1.D0
          P2 = D
          PJ = 0.D0
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
          DO 31 I = 1, MAXN
            YY1(I) = 0.D0
            YY2(I) = 0.D0
            YY3(I) = 0.D0
            YY6(I) = 0.D0
            YY7(I) = 0.D0
            YY9(I) = 0.D0
            YY10(I)= 0.D0
            YY11(I)= 0.D0
 31       YY12(I)= 0.D0
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
          IY1 = 1
          IY2 = 1
          IY3 = 1
          IY6 = 1
          IY7 = 1
          IY9 = 1
          IY10= 1
          IY11= 1
          IY12= 1

          IF (ITYPE.EQ.1) KKK = 2
          IF (ITYPE.NE.1) KKK = 1

          IHLK = NI(NUP)
C        IF(ANGLEI.GT.2.0 .AND. IHLK.GT.15000) IHLK=15000
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  SUMMING UP FOR EVERY DEGREE 'N':
C
 40       DO 90 N = KKK, IHLK

C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  SINCE LOVE NUMBERS ARE NOMALIZED, THEY SHOULD BE UNNOMALIZED. THE
C  FOLLOWING ARE UNNOMALIZING FACTORS, TO BE USED BELOW:
C
            XN = N
C
            H0   = LOVEH(N)
            HL0  = LOVEL(N)
            HK0  = LOVEK(N)
            HLT  = LOVET(N)
            H01  = LOVEH1(N)
            HL01 = LOVEL1(N)
            HK01 = LOVEK1(N)
            HLT1 = LOVET1(N)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  SPECIAL TREATMENT FOR LEGENDRE FUNCTIONS BEGINING WITH N=2:
C
            IF (KKK.EQ.2.AND.N.EQ.2) THEN
              P  = 1.D0
              PA = X
              PB = 1.D0/2.D0*(3.D0*X**2-1.D0)
              P1 = D
              P2 = 1.D0/2.D0*(3.D0*D**2-1.D0)
              PI = 1.D0
              PJ =-(2.D0*(1.D0+CO2)/(1.D0-CO2)+2.D0)*PA
     &        +4.D0*X/(1.D0-CO2)*P
            ENDIF
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  CALCULATING EVERY LEGENDRE FUNCTION:
C
      IF (ANGLEI.GT.0.0035) THEN
C        FOR M=0:
          PP   = P
          P    = PA
          PA   = PB
          PB   = ((2.D0*XN+1.D0)*X*PA-XN*P)/(XN+1.D0)
          PD   = XN/Y*(X*PA-P)
          PDD  = XN/Y/Y*(X*P-(XN*Y*Y+1.D0)*PA)

C        FOR M=1:
          PN1  = - PD
          PN1D = 1.D0/Y**2*(-XN**2*X**2*PA
     &    +(2.D0*XN*XN-1.D0)*X*P-(XN*XN-1.D0)*PP)
          PN1DD= XN/Y**3*(X*(-2+XN*XN*Y*Y)*PA
     &    + (1-XN-XN*XN+(1+XN+XN*XN)*X*X)*P)

C        FOR M=2:
          PI   = PJ
          PJ   =-(2.D0*XN*(1.D0+CO2)/(1.D0-CO2)+XN*(XN+1.D0))*PA
     &    +4.D0*XN*X/(1.D0-CO2)*P
          PJD  = 1.D0/Y*(-(XN+2.D0)*PI+XN*X*PJ)
          PJDD = XN/Y**4*((-4*(XN+1)+XN*(XN+1)**2*Y**2+4*(XN-2)*X**2
     &    -XN*XN*(XN-1)*X**2*Y**2)*PA
     &    +(2.D0*(1+3*XN-XN*XN)+(8-9*XN-XN*XN)*Y**2
     &    +2*(XN*XN-3*XN+5)*X**2)*X*P)
        ENDIF

        IF (ANGLEI.LT.0.0035) THEN
C        FOR M=0:
            PZ1=1-(XN*(1+XN))/4.0*CT**2
     &      +XN*(-2+XN*(1+3*XN*(2+XN)))/192.0*CT**4
            PZ2=-(XN*(1+XN))/2.0*CT
     &      +XN*(-2+XN*(1+3*XN*(2+XN)))/48.0*CT**3
            PZ3=-XN*(1+XN)/2.0+XN*(1+XN)*(2+XN+XN**2)/16.0*CT**2
     &      -(XN*(1+XN)*(-8+XN*(1+XN)*(6+XN+XN**2)))/384.0*CT**4
            PZ4=-XN*(1+XN)/2.0+XN*(-2+XN*(1+3*XN*(2+XN)))/16.0*CT**2
     &      -(XN*(1+XN)*(8+5*(-1+XN)*XN*(1+XN)*(2+XN)))/384.0*CT**4

C        FOR M=1:
            PO1=XN*(1+XN)/2.0*CT-(XN*(-2+XN*(1+3*XN*(2+XN))))/48.0*CT**3
            PO2=XN*(1+XN)/2.0-(XN*(-2+XN*(1+3*XN*(2+XN))))/16.0*CT**2
     &      +XN*(1+XN)*(8+5*(-1+XN)*XN*(1+XN)*(2+XN))/384.0*CT**4
            PO3=-(XN*(-2+XN*(1+3*XN*(2+XN))))/8.0*CT
     &      +XN*(1+XN)*(8+5*(-1+XN)*XN*(1+XN)*(2+XN))/96.0*CT**3
            PO4=(XN*(1+XN))/(2.0*CT)
     &      -(XN*(1+XN)*(-2+3*XN)*(5+3*XN))/48.0*CT
     &      +XN*(1+XN)*(296+15*XN*(1+XN)*(-22+5*XN*(1+XN)))/5760.0*CT**3
            PO5=(XN*(1+XN))/(2.0*CT)
     &      -(XN*(1+XN)*(1+3*XN)*(2+3*XN))/48.0*CT
     &      +XN*(1+XN)*(-184+15*XN*(1+XN)*(14+5*XN*(1+XN)))/5760.0*CT**3
            PO6=XN*(1+XN)/2.0-((-1+XN)*XN*(1+XN)*(2+XN))/16.0*CT**2
     &      +(-1+XN)*XN*(1+XN)*(2+XN)*(-4+XN+XN**2)/384.0*CT**4
            PO7=XN*(1+XN)/2.0-(XN*(1+XN)*(2+XN+XN**2))/16.0*CT**2
     &      +XN*(1+XN)*(-8+XN*(1+XN)*(6+XN+XN**2))/384.0*CT**4
            PO8=(XN*(1+XN))/(2.0*CT)
     &      -(XN*(1+XN)*(-10+3*XN*(1+XN)))/48.0*CT
     &      +XN*(1+XN)*(296+15*XN*(1+XN)*(-10+XN+XN**2))/5760.0*CT**3
            PO9=(XN*(1+XN))/(2.0*CT)
     &      -(XN*(1+XN)*(2+3*XN*(1+XN)))/48.0*CT
     &      +XN*(1+XN)*(-184+15*XN*(1+XN)*(2+XN+XN**2))/5760.0*CT**3
            PO10=XN*(1+XN)/4.0*CT-(XN*(1+XN)*(-8+9*XN*(1+XN)))/96.0*CT**3
            PO11=XN*(1+XN)/4.0*CT**2+XN*(4+XN-6*XN**2-3*XN**3)/96.0*CT**4
            PO12=XN*(1+XN)/4.0*CT-(XN*(1+XN)*(-8+3*XN*(1+XN)))/96.0*CT**3

C        FOR M=2:
            PT1=(-1+XN)*XN*(1+XN)*(2+XN)/8.0*CT**2
     &      -(XN*(1+XN)*(-2+XN+XN**2)**2)/96.0*CT**4
            PT2=(-1+XN)*XN*(1+XN)*(2+XN)/4.0*CT
     &      -(XN*(1+XN)*(-2+XN+XN**2)**2)/24.0*CT**3
            PT3=(-1+XN)*XN*(1+XN)*(2+XN)/4.0
     &      -(XN*(1+XN)*(-2+XN+XN**2)**2)/8.0*CT**2
     &      +(-1+XN)*XN*(1+XN)*(2+XN)*(136+5*XN*(1+XN)*(-14+3*XN*(1+XN)))
     &      /1536.0*CT**4
            PT4=(-1+XN)*XN*(1+XN)*(2+XN)/4.0
     &      -((-1+XN)*XN*(1+XN)*(2+XN)*(-3+XN+XN**2))/24.0*CT**2
     &      +(-1+XN)*XN*(1+XN)*(2+XN)*(168+XN*(1+XN)*(-74+9*XN*(1+XN)))
     &      /4608.0*CT**4
            PT5=(-1+XN)*XN*(1+XN)*(2+XN)/4.0
     &      -(XN**2*(1+XN)**2*(-2+XN+XN**2))/24.0*CT**2
     &      +(-1+XN)*XN*(1+XN)*(2+XN)*(-72+XN*(1+XN)
     &      *(22+9*XN*(1+XN)))/4608.0*CT**4
            PT6=(-1+XN)*XN*(1+XN)*(2+XN)/8.0*CT
     &      -((-1+XN)*XN*(1+XN)*(2+XN)*(-4+XN+XN**2))/96.0*CT**3
            PT7=(-1+XN)*XN*(1+XN)*(2+XN)/8.0*CT
     &      -((-1+XN)*XN*(1+XN)*(2+XN)*(2+XN+XN**2))/96.0*CT**3
            PT8=(-1+XN)*XN*(1+XN)*(2+XN)/8.0
     &      -((-2+XN)*(-1+XN)*XN*(1+XN)*(2+XN)*(3+XN))/96.0*CT**2
     &      +(-2+XN)*(-1+XN)*XN*(1+XN)*(2+XN)*(3+XN)*(-28+3*XN*(1+XN))
     &      /9216.0*CT**4
            PT9=(-1+XN)*XN*(1+XN)*(2+XN)/8.0
     &      -(XN**2*(1+XN)**2*(-2+XN+XN**2))/96.0*CT**2
     &      +(-1+XN)*XN*(1+XN)*(2+XN)*(-72+XN*(1+XN)*(2+3*XN*(1+XN)))
     &      /9216.0*CT**4
            PT10=(-1+XN)*XN*(1+XN)*(2+XN)/8.0*CT**2
     &      -((-1+XN)*XN*(1+XN)*(2+XN)*(-5+2*XN*(1+XN)))/96.0*CT**4
            PT11=(-1+XN)*XN*(1+XN)*(2+XN)/16.0*CT**3
            PT12=(-1+XN)*XN*(1+XN)*(2+XN)/16.0*CT**2
     &      -((-1+XN)*XN*(1+XN)*(2+XN)*(-5+XN+XN**2))/192.0*CT**4
          ENDIF

C        FOR DISK FACTOR:
          P0   = P1
          P1   = P2
          P2   = ((2.D0*XN+1.D0)*D*P1-XN*P0)/(XN+1.D0)
          IF  (N.GT.0) XX = (P2-D*P1)/XN/(D-1.D0)
      If (int(XN)==1) HK0 = 6.67d-11/6371d3*tt

C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF (ANGLEI.LT.0.0035) THEN
            IF (ITYPE.EQ.1) THEN
              U    = -2.D0 * H0 *PT1                                 *XX
              VCT  = -2.D0 * (HL0 *PT2 + 2.D0*HLT*PT6)               *XX
              VLM  = -2.D0 * (2.D0*HL0 *PT6 + HLT*PT2)               *XX
              ERR  =  2.D0*SL/(SL+2.D0*SM)*(H0*2.D0-HL0*XN*(XN+1.D0))*PT1*XX
              ECC  =  2.D0*(-HL0*PT3-H0*PT1-2.D0*HLT*(PT4-PT9))      *XX
              EFF  =  2.D0*(HL0*(4*PT8-PT5)-H0*PT1+2*HLT*(PT4-PT9))  *XX
              ECF  =  2.D0*(HL0*4*(-PT4+PT9)+HLT*(PT5-4*PT8-PT3))    *XX
              FI   = -2.D0 * HK0 *PT1                 *XX
              GRAV = -2.D0 * HK0*(XN+1) *PT1          *XX
            ELSEIF (ITYPE.EQ.2) THEN
              U    = -2.D0 * H0 *PO1                                 *XX
              VCT  = -2.D0 * (HL0 *PO2  - HLT*PO6)                   *XX
              VLM  = -2.D0 * (HL0 *PO6 - HLT*PO2)                    *XX
              ERR  = 2.D0*SL/(SL+2.D0*SM)*(H0*2.D0-HL0*XN*(XN+1.D0))*PO1*XX
              ECC  = 2.D0*(-HL0*PO3-H0*PO1+HLT*(PO4-PO9))            *XX
              EFF  = 2.D0*(HL0*(PO8-X*PO5)-H0*PO1-HLT*(PO4-PO9))     *XX
              ECF  =  2.D0*(HL0*2*(-PO4+PO9)-HLT*(PO5-PO8-PO3))      *XX
              FI   = -2.D0 * HK0 *PO1                 *XX
              GRAV = -2.D0 * HK0*(XN+1) *PO1          *XX
            ELSEIF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
              U    = H0  *PZ1                                        *XX
              VCT  = HL0 *PZ2                                        *XX
              VLM  = 0.D0
              ERR  =  SL/(SL+2.D0*SM)*(HL0*XN*(XN+1.D0)-2.D0*H0) * PZ1*XX
              ECC  =  (HL0 * PZ4 +H0*PZ1)                             *XX
              EFF  =  (HL0*PZ3+H0*PZ1)                                *XX
              ECF  =   0.D0
              FI   =  HK0 *PZ1                 *XX
C             GRAV =  HK0*(XN+1) *PZ1          *XX
C           GRAV =  (2*H0-(XN+1)*HK0) *PZ1   *XX
C             ! gravity change due to internal loading
              GRAV = (2*H0*9.8156D0/6371D3-(XN+1)/6371D3
     &               *(HK0-tt**(int(XN)*(6.67D-11)/6371D3)))*PZ1   *XX

            ENDIF
          ENDIF
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (ANGLEI.GT.0.0035) THEN
          IF (ITYPE.EQ.1) THEN
            U    = -2.D0 * H0 *PJ                                   *XX
            VCT  = -2.D0 * (HL0 *PJD + 2.D0*HLT*PJ/Y)               *XX
            VLM  = -2.D0 * (2.D0*HL0 *PJ/Y + HLT*PJD)               *XX
            ERR  =  2.D0*SL/(SL+2.D0*SM)*(H0*2.D0-HL0*XN*(XN+1.D0))*PJ*XX
            ECC  =  2.D0*(-HL0*PJDD-H0*PJ-2.D0*HLT*(PJD/Y-X/Y/Y*PJ))*XX
            EFF  =  2.D0*(HL0/Y*(4*PJ/Y-X*PJD)-H0*PJ
     &      +2*HLT/Y*(PJD-X/Y*PJ))                      *XX
            ECF  =  2.D0*(HL0/Y*4*(-PJD+X/Y*PJ)
     &      +HLT*(X/Y*PJD-4*PJ/Y/Y-PJDD))*XX
            FI   =  -2.D0*HK0 *PJ                           *XX
            GRAV =  -2.D0*(XN+1)*HK0 *PJ                    *XX
          ELSEIF (ITYPE.EQ.2) THEN
            U    = -2.D0 * H0 *PN1                                  *XX
            VCT  = -2.D0 * (HL0 *PN1D  - HLT*PN1/Y)                 *XX
            VLM  = -2.D0 * (HL0 *PN1/Y - HLT*PN1D)                  *XX
            ERR  = 2.D0*SL/(SL+2.D0*SM)*(H0*2.D0-HL0*XN*(XN+1.D0))*PN1*XX
            ECC  = 2.D0*(-HL0*PN1DD-H0*PN1+HLT*(PN1D/Y-X/Y/Y*PN1))  *XX
            EFF  = 2.D0*(HL0/Y*(PN1/Y-X*PN1D)-H0*PN1
     &      -HLT/Y*(PN1D-X/Y*PN1))                      *XX
            ECF  =  2.D0*(HL0/Y*2*(-PN1D+X/Y*PN1)
     &      -HLT*(X/Y*PN1D-PN1/Y/Y-PN1DD))              *XX
            FI   =  -2.D0*HK0 *PN1                           *XX
            GRAV =  -2.D0*(XN+1)*HK0 *PN1                    *XX
          ELSEIF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
            U    = H0  *PA                                          *XX
            VCT  = HL0 *PD                                          *XX
            VLM  = 0.D0
            ERR  =  SL/(SL+2.D0*SM)*(HL0*XN*(XN+1.D0)-2.D0*H0) * PA*XX
            ECC  =  (HL0 * PDD +H0*PA)                              *XX
            EFF  =  (X/Y*HL0*PD+H0*PA)                              *XX
            ECF  =   0.D0
            FI   =  HK0 *PA                           *XX
C           GRAV =  ((XN+1)*HK0) *PA             *XX
C         GRAV =  (2*H0-(XN+1)*HK0) *PA             *XX

C           ! gravity change due to internal loading
            GRAV = (2*H0*9.8156D0/6371D3-(XN+1)/6371D3
     &               *(HK0-(tt**(int(XN))*(6.67D-11)/6371D3)))*PA   *XX


          ENDIF
        ENDIF
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF((IY1+2).GE.MAXN.OR.(IY2+2).GE.MAXN.OR.(IY3+2).GE.MAXN) GOTO 777
      IF(U   .GE.0.D0.AND.DSIGN(1.D0,TU)   .NE. DSIGN(1.D0,U))
     &IY1 = IY1+2
      IF(VCT .GE.0.D0.AND.DSIGN(1.D0,TVCT) .NE. DSIGN(1.D0,VCT))
     &IY2 = IY2+2
      IF(VLM .GE.0.D0.AND.DSIGN(1.D0,TVLM) .NE. DSIGN(1.D0,VLM))
     &IY3 = IY3+2

      IF(U.GE.0.D0)    YY1(IY1)   = YY1(IY1)  +U
      IF(U.LT.0.D0)    YY1(IY1+1) = YY1(IY1+1)+U
      IF(VCT.GE.0.D0)  YY2(IY2)   = YY2(IY2)  +VCT
      IF(VCT.LT.0.D0)  YY2(IY2+1) = YY2(IY2+1)+VCT
      IF(VLM.GE.0.D0)  YY3(IY3)   = YY3(IY3)  +VLM
      IF(VLM.LT.0.D0)  YY3(IY3+1) = YY3(IY3+1)+VLM

      TU    = U
      TVCT  = VCT
      TVLM  = VLM

 777  CONTINUE

C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      IF((IY6+2).GE.MAXN.OR.(IY7+2).GE.MAXN.OR.(IY9+2).GE.MAXN.OR.
     &(IY10+2).GE.MAXN.OR.(IY11+2).GE.MAXN.OR.
     &(IY12+2).GE.MAXN) GOTO 888
      IF(ERR.GE.0.D0.AND.DSIGN(1.D0,TERR).NE.DSIGN(1.D0,ERR))
     &IY6 = IY6+2
      IF(ECC.GE.0.D0.AND.DSIGN(1.D0,TECC).NE.DSIGN(1.D0,ECC))
     &IY7 = IY7+2
      IF(EFF.GE.0.D0.AND.DSIGN(1.D0,TEFF).NE.DSIGN(1.D0,EFF))
     &IY9 = IY9+2
      IF(FI.GE.0.D0.AND.DSIGN(1.D0,TFI).NE.DSIGN(1.D0,FI))
     &IY10 = IY10+2
      IF(GRAV.GE.0.D0.AND.DSIGN(1.D0,TGRAV).NE.DSIGN(1.D0,GRAV))
     &IY11 = IY11+2
      IF(ECF.GE.0.D0.AND.DSIGN(1.D0,TECF).NE.DSIGN(1.D0,ECF))
     &IY12 = IY12+2

      IF(ERR.GE.0.D0)   YY6(IY6)     = YY6(IY6)    +ERR
      IF(ERR.LT.0.D0)   YY6(IY6+1)   = YY6(IY6+1)  +ERR
      IF(ECC.GE.0.D0)   YY7(IY7)     = YY7(IY7)    +ECC
      IF(ECC.LT.0.D0)   YY7(IY7+1)   = YY7(IY7+1)  +ECC
      IF(EFF.GE.0.D0)   YY9(IY9)     = YY9(IY9)    +EFF
      IF(EFF.LT.0.D0)   YY9(IY9+1)   = YY9(IY9+1)  +EFF
      IF(FI.GE.0.D0)    YY10(IY10)   = YY10(IY10)  +FI
      IF(FI.LT.0.D0)    YY10(IY10+1) = YY10(IY10+1)+FI
      IF(GRAV.GE.0.D0)  YY11(IY11)   = YY11(IY11)  +GRAV
      IF(GRAV.LT.0.D0)  YY11(IY11+1) = YY11(IY11+1)+GRAV
      IF(ECF.GE.0.D0)   YY12(IY12)   = YY12(IY12)  +ECF
      IF(ECF.LT.0.D0)   YY12(IY12+1) = YY12(IY12+1)+ECF

      TERR  = ERR
      TECC  = ECC
      TEFF  = EFF
      TFI  = FI
      TGRAV  = GRAV
      TECF  = ECF

 888  CONTINUE

C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
      IF (N.EQ.IHLK) THEN
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  EUTER'S TRANSFORM (TO SPEED UP CONVERGENCE):
C
                CALL  WJNG1(YY1, IY1+1, 1.D-56,U,   IEPSI1)
                CALL  WJNG1(YY2, IY2+1, 1.D-56,VCT, IEPSI2)
                CALL  WJNG1(YY3, IY3+1, 1.D-56,VLM, IEPSI3)
                CALL  WJNG1(YY6, IY6+1, 1.D-56,ERR, IEPSI6)
                CALL  WJNG1(YY7, IY7+1, 1.D-56,ECC, IEPSI7)
                CALL  WJNG1(YY9, IY9+1, 1.D-56,EFF, IEPSI9)
                CALL  WJNG1(YY10,IY10+1,1.D-56,FI,  IEPSI10)
                CALL  WJNG1(YY11,IY11+1,1.D-56,GRAV,IEPSI11)
                CALL  WJNG1(YY12,IY12+1,1.D-56,ECF, IEPSI12)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  NOTE THAT: 0 DEGREE OF LOVE NUMBER (HN0) FOR TYPES 3,4 IS CONSIDERED
C  IN THE FOLLOWING: (I.E. U = U + HN0).
C
                ERR2  =  ERR
                ECC2  =  ECC
                EFF2  =  EFF
                EFI2  =  FI
                EGRAV2=  GRAV
                ECF2  =  ECF
                IF(ITYPE.EQ.3.OR.ITYPE.EQ.4) THEN
                  U = U + HN0
                  ERR2   =  ERR + SL/(SL+2.D0*SM)*(-2*HN0)
                  ECC2   =  ECC + HN0
                  EFF2   =  EFF + HN0
                  EFI2   =  FI
C                 EGRAV2 =  GRAV
                  EGRAV2 =  GRAV+(6.67D-11-SL5)/(6371D3**2)+2d0*9.8156D0*HN0/6371D3
                  ECF2   =  ECF
                ENDIF
                ERRSUM   = ERR2   / (RADIS-DEFORM)
                ECCSUM   = ECC2   / (RADIS-DEFORM)
                EFFSUM   = EFF2   / (RADIS-DEFORM)
                EFISUM   = EFI2
C               EGRAVSUM = EGRAV2 * 9.8156D0/(6371d0**3)
                EGRAVSUM = EGRAV2
                ECFSUM   = ECF2   / (RADIS-DEFORM)



C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C       IF (ITYPE .EQ. 4) THEN
C         call asymptotic_love(
C    &         SL*13080*1d8,SM*13080*1d8,(6371.0-DEP)*1D3,6371D3,HL)
C
C         delta = (6371d0-DEP)/6371d0
C         call asymptotic_green(HL,delta,grs,ANGLEI,disp2)
C         ! write(126,'(F10.4,2(2X,ES16.5))') ANGLEI,disp2
C
C         U   = U   + disp2(1)
C         VCT = VCT + disp2(2)
C
C         ! write(*,'(A,2(2X,ES16.5))') 'disp2 = ',disp2
C
C       END IF
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C  OUTPUT RESULTS INTO FILES:
C
                OUTANGLE = ANGLEI

C      IF (IANGLE .EQ. MANGLE) OUTANGLE = 180.0


                IF(ITYPE.EQ.1) WRITE(51,100) OUTANGLE,U,VCT,VLM
     &          ,ERRSUM,ECCSUM,EFFSUM,ECFSUM,EFISUM
C     &               ,EGRAVSUM
     &          ,EGRAVSUM
                IF(ITYPE.EQ.2) WRITE(52,100) OUTANGLE,U,VCT,VLM
     &          ,ERRSUM,ECCSUM,EFFSUM,ECFSUM,EFISUM
C     &               ,EGRAVSUM
     &          ,EGRAVSUM
                IF(ITYPE.EQ.3) WRITE(53,100) OUTANGLE,U,VCT,VLM
     &          ,ERRSUM,ECCSUM,EFFSUM,ECFSUM,EFISUM
C     &               ,EGRAVSUM
     &          ,EGRAVSUM-2.0d0*U
                IF(ITYPE.EQ.4) WRITE(54,100) OUTANGLE,U,VCT,VLM
     &          ,ERRSUM,ECCSUM,EFFSUM,ECFSUM,EFISUM
C     &               ,EGRAVSUM
     &          ,EGRAVSUM

              ENDIF

 90       CONTINUE
399     CONTINUE

 333    FORMAT(10I6)
 123    FORMAT(4(1E18.8,2X))
 100    FORMAT(1X,1F11.6,1X,9(2x,1E16.6))

        END
C  ################### SUBROUTINE ##################################
C
C  *************************************
C  *                                   *
C  *         EULER'S TRANSFORM (2)     *
C  *                                   *
C  *************************************
C
      SUBROUTINE WJNG1(A,N,CRT,SUM,ICODE)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     CALCULATES  SUM OF ALTERNATING SERIES ;
C
C     SUM = A(1)+A(2)+A(3)+.....
C
C     ICODE =-1 :          TERMS LESS THAN 20, DIRECT SUM
C                          INSTEAD OF EULER TRANSFORMATION.
C     ICODE = 0 :          CONVERGED (ALTERNATING SERIES)
C     ICODE = 1 : NOT YET  CONVERGED (ALTERNATING SERIES)
C
C  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT  REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      PARAMETER (MXWORK=30005)
      DIMENSION A(N),WORK(MXWORK)

      IF(N.LT.20) GO TO 100
      RATIO = 1.D7
      M     = 0
      L     = 0
      SUM   = 0.D0

      DO 1 I=1,MXWORK
    1 WORK(I)=0.D0
      WORK(1)= A(1)
      SUM    = 0.5D0*A(1)

   10 MM=M+2
      JJ=L+MM
      IF(JJ.GT.N.OR.MM.GT.MXWORK) GO TO 99
      FNEXT=A(JJ)
      DO 20 I=1,MM
        TT     = 0.5D0*(WORK(I)+FNEXT)
        WORK(I)= FNEXT
        FNEXT  = TT
   20 CONTINUE

      SOLD=SUM
      IF(DABS(WORK(MM)).GT.DABS(WORK(MM-1)) ) GO TO 30
      M   = M+1
      SUM = SUM+0.5D0*WORK(MM)
      GO TO 40
   30 L=L+1
      SUM=SUM+WORK(MM)
   40 CONTINUE

      IF(SOLD.EQ.0.D0) GO TO 50
      RATIO=SUM/SOLD-1.D0
   50 CONTINUE

      IF(DABS(RATIO).GE.CRT) GO TO 10
      ICODE = 0
      RETURN

   99 ICODE = 1
      RETURN
   60 FORMAT(1H ,T95,' CONV.  L=',I6,' M=',I6)
   61 FORMAT(//1H ,T5,'NOT YET CONV. L=',I6,' M=',I6/)

  100 CONTINUE
      ICODE = -1
      SUM=0.D0
      DO 120 I=1,N
  120 SUM=SUM+A(I)
   80 FORMAT(1H ,T95,'ALL SUMMED',I6)
      END
