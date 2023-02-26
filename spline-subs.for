c     http://astronomy.nmsu.edu/cwc/Teaching/Codes/spline-subs.f
c     this module contains four Fortran subroutines useful for 1D and
c     2D cubic spline interpolation; the subroutines are

c     1D interpolation subroutines
c     SPLINE(X,Y,N,YP1,YPN,Y2)
c     SPLINT(XA,YA,Y2A,N,X,Y)

c     2D interpolation subroutines
c     SPLIE2(X1A,X2A,YA,M,N,Y2A)
c     SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

c     how these work: consider 1D interpolation, for which you use the
c     SPLINE and SPLINT.  The interpolation is performed in a two step
c     process.  Step 1 is that you compute the 2nd derivatives of the Y
c     array data at each X array data point.  This is accomplished using
c     a single call to subroutine SPLINE.  Step 2 is the interpolation
c     itself, which is accomplished by calling SPLINT.  Multiple
c     successive calls to SPLINT can be made to interpolate the Y array
c     data to obtain an interpolated y value at various x values.

c     further explanations for their use are given in comments within
c     the subroutines

c..............................................................................

      SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)

c     SPLINE use: given an 1D array of X data and an array of Y data,
c     both of length N, this routine computes the 2nd derivatives, Y2 at
c     each X data point.  The user needs to specify the values of YP1
c     and YP2, which flags how the Y2 are computed at the edges.  For
c     natural spline fitting (recommended), set YP1 and YPN to numbers
c     greater than 1.0E+30.

c     this routine called once, prior to using routine SPLINT, as a set
c     up for using routine SPLINT, which performs the actual
c     interpolation

c     IMPORTANT NOTE: the X data values in array X must be in ascending
c     order or the interpolation will fail

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NMAX=600)
      DIMENSION X(N),Y(N),Y2(N),U(NMAX)

c     if YP1>1.0E+30 use natural spline, otherwise estimate Y2 at the
c     first point

      IF (YP1.GT..99D30) THEN
        Y2(1)=0.
        U(1)=0.
      ELSE
        Y2(1)=-0.5
        U(1)=(3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF

c     store intermediate values of terms in the expansion series

      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+2.
        Y2(I)=(SIG-1.)/P
        U(I)=(6.*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))
     *      /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE

c     if YPN>1.0E+30 use natural spline, otherwise estimate Y2 at the
c     last point point

      IF (YPN.GT..99D30) THEN
        QN=0.
        UN=0.
      ELSE
        QN=0.5
        UN=(3./(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+1.)

c     compute the Y2 from the 2nd order expansion series

      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE

      RETURN
      END

c..............................................................................

      SUBROUTINE SPLINT(XA,YA,Y2A,N,X,Y)

c     SPLINT use: given an 1D array of XA data, an array of YA data, and
c     an array of the 2nd derivatives Y2A, all of length N, this routine
c     performs cubic spline interpolation, returning the interpolated
c     value Y at the user input value X.  The Y2A are computed in
c     routine SPLINE, which is called once before calling SPLINT.

c     IMPORTANT NOTE: the X data values in array X must be in ascending
c     order or the interpolation will fail

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XA(N),YA(N),Y2A(N)

      KLO=1
      KHI=N

c     determine the indices of array XA that bracket the input X value

1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF

c     determine the finite difference along the X dimension

      H=XA(KHI)-XA(KLO)
      IF (H.EQ.0.) STOP 'Bad XA input in routine SPLINE.'

c     interpolate

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     *      ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/6.

      RETURN
      END



c..............................................................................

      SUBROUTINE SPLIE2(X1A,X2A,YA,M,N,Y2A)

c     SPLIE2 use: given an array of X1A data of length M, and an array
c     of X2A data of length N, this routine computes the 2nd
c     derivatives, Y2A, at each X1A,X2A data point.  Thus Y2A has
c     dimensions Y2A(M,N).  Natural spline fitting is assumed.

c     this routine called once, prior to using routine SPLIN2, as a set
c     up for using routine SPLIN2, which performs the actual
c     interpolation.

c     Uses routines: SPLINE

c     IMPORTANT NOTE: the X1A and X2A data values must both be in
c     ascending order or the interpolation will fail

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=600)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN)

      DO 13 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
11      CONTINUE
        CALL SPLINE(X2A,YTMP,N,1.D30,1.D30,Y2TMP)
        DO 12 K=1,N
          Y2A(J,K)=Y2TMP(K)
12      CONTINUE
13    CONTINUE
      RETURN
      END

c..............................................................................

      SUBROUTINE SPLIN2(X1A,X2A,YA,Y2A,M,N,X1,X2,Y)

c     SPLIN2 use: given an array of X1A data of length M, an array of
c     X2A data of length N, and an array of 2nd derivatives, Y2A at each
c     X1A,X2A data point, dimensioned Y2A(M,N), this routine performs 2D
c     interpolation, returning the interpolated value Y at user input
c     values X1 and X2.  Natural spline fitting is assumed.

c     Uses routines: SPLINT, SPLINE

c     IMPORTANT NOTE: the X1A and X2A data values must both be in
c     ascending order or the interpolation will fail


      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NN=600)
      DIMENSION X1A(M),X2A(N),YA(M,N),Y2A(M,N),YTMP(NN),Y2TMP(NN),YYTMP(
     *NN)
      DO 12 J=1,M
        DO 11 K=1,N
          YTMP(K)=YA(J,K)
          Y2TMP(K)=Y2A(J,K)
11      CONTINUE
        CALL SPLINT(X2A,YTMP,Y2TMP,N,X2,YYTMP(J))
12    CONTINUE
      CALL SPLINE(X1A,YYTMP,M,1.D30,1.D30,Y2TMP)
      CALL SPLINT(X1A,YYTMP,Y2TMP,M,X1,Y)
      RETURN
      END

c     eof
