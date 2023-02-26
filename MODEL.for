      SUBROUTINE MODEL(FILEMODEL)
      IMPLICIT REAL*8(A-H,O-Z)
      IMPLICIT INTEGER*8(I-N)
      CHARACTER*30  FILEMODEL
      DIMENSION DR1(100),DG1(100),DL1(100),DP1(100)
      DIMENSION DR2(100),DG2(100),DL2(100),DP2(100)
      DIMENSION DMM1(100),DMK1(100),DEM1(100),DEK1(100)
      DIMENSION DMM2(100),DMK2(100),DEM2(100),DEK2(100)

          GNEWTN = 6.67D-8
          PAI    = 3.14159265358979D0

      OPEN(100,FILE=FILEMODEL,ACTION='READ')
      OPEN(103,FILE='model1_incore.dat')
      OPEN(104,FILE='model2_outcore.dat')
      OPEN(105,FILE='model3_mantle.dat')
           READ(100,*)
           READ(100,*) N
           READ(100,*)
           READ(100,*)
           READ(100,*)
          DO I=1,N
            J=N+1-I
           READ(100,*) UNM,DR2(J),DP2(J)
     &       ,DL2(J),DMM2(J),DMK2(J),DEM2(J),DEK2(J)
           READ(100,*) NUM,DR1(J),DP1(J)
     &       ,DL1(J),DMM1(J),DMK1(J),DEM1(J),DEK1(J)
           DR1(J)=6371.0-DR1(J);DR2(J)=6371.0-DR2(J)
           DL1(J)=1.0D1*DL1(J);DL2(J)=1.0D1*DL2(J)
           DMM1(J)=1.0D1*DMM1(J);DMM2(J)=1.0D1*DMM2(J)
           DMK1(J)=1.0D1*DMK1(J);DMK2(J)=1.0D1*DMK2(J)
          END DO

          DO J=1,N-1
             IF(DMM1(J+1).LT.0.001 .AND.DMM2(J).GT.0.001)  M1=J
             IF(DMM1(J+1).GT.0.001 .AND.DMM2(J).LT.0.001)  M2=J
          END DO
           WRITE(*,'(2I5)') M1,M2

          DO J=1,N
            DR=DR1(J);DP=DP1(J);UR=DR2(J);UP=DP2(J)
          IF(J.EQ.1) DM0=4.D0/3.D0*PAI*DP*DR**3
          DG0=DM0*GNEWTN/(DR**2)*1.0D5
          DG1(J)=DG0
            DPK=(UP+DP)/2
          DM0=DM0+4.0*PAI*DPK*UR**3/3.0-4.0*PAI*DPK*DR**3/3.0
            DG0=DM0*GNEWTN/(UR**2)*1.0D5
          DG2(J)=DG0
          END DO

             WRITE(103,'(1I8)') M1
          DO J=1,M1
            WRITE(103,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR1(J),DG1(J),DP1(J),DL1(J),DMM1(J),DMK1(J),DEM1(J),DEK1(J)
            WRITE(103,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR2(J),DG2(J),DP2(J),DL2(J),DMM2(J),DMK2(J),DEM2(J),DEK2(J)
          END DO

             WRITE(104,'(1I8)') M2-M1
          DO J=M1+1,M2
            WRITE(104,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR1(J),DG1(J),DP1(J),DL1(J),DMM1(J),DMK1(J),DEM1(J),DEK1(J)
            WRITE(104,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR2(J),DG2(J),DP2(J),DL2(J),DMM2(J),DMK2(J),DEM2(J),DEK2(J)
          END DO

             WRITE(105,'(1I8)') N-M2
          DO J=M2+1,N
            WRITE(105,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR1(J),DG1(J),DP1(J),DL1(J),DMM1(J),DMK1(J),DEM1(J),DEK1(J)
            WRITE(105,'(1F12.4,1F14.6,2E16.6,2F14.6,2E10.2)')
     &       DR2(J),DG2(J),DP2(J),DL2(J),DMM2(J),DMK2(J),DEM2(J),DEK2(J)
          END DO

       CLOSE(105)
       CLOSE(104)
       CLOSE(103)
       CLOSE(100)

      END
