C 
C**HP**C ATTENTION: Normalization has been changed by HP
C
C SUB PDFZBETA
C
C ===================================================
      SUBROUTINE PDFZBETA(ZMEAN,ZVAR,Z,PDFZ,NJ)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C      EXTERNAL DBETA
      DOUBLE PRECISION   DZ(5000)
      DOUBLE PRECISION   log,exp
      DIMENSION PDFZ(NJ),Z(NJ)
      IF(NJ.GT.5000) THEN
        WRITE(6,*)'******* DIMENSIONS IN PDFZBETA TOO SMALL: ', NJ
        STOP
      END IF
      ZM=ZMEAN
      ZV=ZVAR
C      IF (ZV.LT.1.0D-5) ZV=0.0D0
      DZ(1)=0.5D0*(Z(2)-Z(1))
      DZ(NJ)=0.5D0*(Z(NJ)-Z(NJ-1))
      DO 12 J=2,NJ-1
      DZ(J)=0.5D0*(Z(J+1)-Z(J-1))
12    CONTINUE
C ---  SET LIMIT TO ZM AND ZV
      SMDELT=1.D-10
      SMZV=ZM*SMDELT
      SMZV1=(1.D0-ZM)*SMDELT
      SMZV=DMIN1(SMZV,SMZV1)
C -------- FOR THE SPECIAL POINT ZM=0 OR 1
17    IF(ZM.LT.SMDELT) THEN
        DO 14 J=2,NJ
        PDFZ(J)=0.D0
14      CONTINUE
        PDFZ(1)=1.D0
      ELSE IF(ZM.GT.(1.D0-SMDELT)) THEN
        DO 16 J=1,NJ-1
        PDFZ(J)=0.D0
16      CONTINUE
        PDFZ(NJ)=1.D0
      ELSE 
C   --- FOR ZV=0 BUT ZM=(0,1)
        IF(ZV.LT.SMZV) THEN
          DO 18 J=2,NJ
          IF(ZM.GT.Z(J)) GOTO 18
          J1=J
          J1M=J-1
          GOTO 20
18        CONTINUE
20        CONTINUE
          DO 22 J=1,NJ
          PDFZ(J)=0.D0
22        CONTINUE
          PDFZ(J1)=(ZM-Z(J1M))/(Z(J1)-Z(J1M))
          PDFZ(J1M)=1.D0-PDFZ(J1)
        ELSE
C    --- FOR THE OTHER CASES IN WHICH ZM,ZV NOT EQUAL TO SPECIAL POINT.
          GAMA=ZM*(1.D0-ZM)/ZV-1.D0
          IF(GAMA.LT.SMDELT) THEN
            DO 24 J=2,NJ-1
            PDFZ(J)=0.D0
24          CONTINUE
            PDFZ(1)=1.0D0-ZM
            PDFZ(NJ)=ZM
          ELSE
            IF (1.ge.0) THEN
              A=ZM*GAMA
              B=GAMA-A
              A1=A-1.0D0
              B1=B-1.0D0
C              VBETA=DBETA(A,B)
CCCHP            VBETA=f77beta(A,B)
              VBETA=f77gammaln(a)+f77gammaln(b)-f77gammaln(a+b)
              IF(VBETA.EQ.0.D0) THEN
                DO 1118 J=2,NJ
                IF(ZM.GT.Z(J)) GOTO 1118
                J1=J
                J1M=J-1
                GOTO 1200
1118             CONTINUE
1200             CONTINUE
                DO 1222 J=1,NJ
                 PDFZ(J)=0.D0
1222             CONTINUE
                PDFZ(J1)=(ZM-Z(J1M))/(Z(J1)-Z(J1M))
                PDFZ(J1M)=1.D0-PDFZ(J1)
                GOTO 223
              ENDIF
              PDFZ(1)=0.D0
              PDFZ(NJ)=0.D0
              PDFNJ2=0.d0
              DO 130 J=2,NJ-1
C23456789012345678901234567890123456789012345678901234567890123456789012
              PDFZ(J)=exp(A1*log(Z(J)) + B1*log(1.D0-Z(J)) - VBETA) 
     &                  * DZ(J)
              PDFNJ2=PDFNJ2+PDFZ(J)
130            CONTINUE
C    ---    FOR THE SPECIAL POINT Z(1) AND Z(NJ)
              IF(A1.GE.0.D0.AND.B1.GE.0.D0) THEN
                PDFZ(1)=DMAX1(Z(1),0.0D0)**A1*(1.D0-Z(1))**B1/VBETA
     &                  * DZ(1)
                PDFZ(NJ)=Z(NJ)**A1*DMAX1(1.D0-Z(NJ),0.0D0)**B1
     &                   /VBETA*DZ(NJ)
              ELSE 
                IF(A1.LT.0.D0) PDFZ(1)=1.D0-PDFNJ2
                IF(B1.LT.0.D0) PDFZ(NJ)=1.D0-PDFNJ2
                IF(A1.LT.0.D0.AND.B1.LT.0.D0) THEN
                  PDFZ(1)=(1.D0-ZM)*(1.D0-PDFNJ2)
                  PDFZ(NJ)=ZM*(1.D0-PDFNJ2)
                END IF
              END IF
            ELSE
              A=ZM*GAMA
              B=GAMA-A
              A1=A-1.0D0
              B1=B-1.0D0
C              VBETA=DBETA(A,B)
              VBETA=f77beta(A,B)
              IF(VBETA.EQ.0.D0) THEN
                DO 118 J=2,NJ
                IF(ZM.GT.Z(J)) GOTO 118
                J1=J
                J1M=J-1
                GOTO 200
118             CONTINUE
200             CONTINUE
                DO 222 J=1,NJ
                 PDFZ(J)=0.D0
222             CONTINUE
                PDFZ(J1)=(ZM-Z(J1M))/(Z(J1)-Z(J1M))
                PDFZ(J1M)=1.D0-PDFZ(J1)
                GOTO 223
              ENDIF
              PDFZ(1)=0.D0
              PDFZ(NJ)=0.D0
              PDFNJ2=0.d0
              DO 30 J=2,NJ-1
              PDFZ(J)=Z(J)**A1*(1.D0-Z(J))**B1/VBETA*DZ(J)
              PDFNJ2=PDFNJ2+PDFZ(J)
30            CONTINUE
C    ---    FOR THE SPECIAL POINT Z(1) AND Z(NJ)
              IF(A1.GE.0.D0.AND.B1.GE.0.D0) THEN
                PDFZ(1)=DMAX1(Z(1),0.0D0)**A1*(1.D0-Z(1))**B1/VBETA
     &                *DZ(1)
                PDFZ(NJ)=Z(NJ)**A1*DMAX1(1.D0-Z(NJ),0.0D0)**B1
     &                   /VBETA*DZ(NJ)
              ELSE 
                IF(A1.LT.0.D0) PDFZ(1)=1.D0-PDFNJ2
                IF(B1.LT.0.D0) PDFZ(NJ)=1.D0-PDFNJ2
                IF(A1.LT.0.D0.AND.B1.LT.0.D0) THEN
                  PDFZ(1)=(1.D0-ZM)*(1.D0-PDFNJ2)
                  PDFZ(NJ)=ZM*(1.D0-PDFNJ2)
                END IF
              END IF
            END IF
          END IF
        END IF 
      END IF
C  ------  CHECK THE SUMMER OF PDFZ(J) WHICH MUST BE EQUAL TO 1.0
223   PDFSUM=0.D0
      DO 32 J=1,NJ
      IF(PDFZ(J).LT.-1.D-30) THEN
        WRITE(6,*)'******* PDFZ().LT.0',PDFZ(J),J,ZM,ZV
C        STOP
      END IF
      PDFSUM=PDFSUM+PDFZ(J)
32    CONTINUE
      IF(PDFSUM.LE.0.D0) THEN
      	ZV = 0.0D0
      	GOTO 17
C        WRITE(6,*)'******* PDFSUM=',PDFSUM
C        STOP
      END IF
      IF(DABS(PDFSUM-1.D0).GT.SMDELT) THEN
       DO 34 J=2,NJ
         PDFZ(J) = PDFZ(J) / PDFSUM
34     CONTINUE
C**HP**Cc ----  find the interval which ZM lies
C**HP**C       DO 34 J=2,NJ
C**HP**C          IF(ZM.GT.Z(J)) GOTO 34
C**HP**C          J1=J
C**HP**C          J1M=J-1
C**HP**C          GOTO 35
C**HP**C34     CONTINUE
C**HP**C35     CONTINUE
C**HP**C       PSUM=1.D0-PDFSUM
C**HP**C       ALF1=(ZM-Z(J1M))/(Z(J1)-Z(J1M))
C**HP**C       ALF2=(Z(J1)-ZM)/(Z(J1)-Z(J1M))
C**HP**C       PDFZ(J1)=ALF1*PSUM+PDFZ(J1)
C**HP**C       PDFZ(J1M)=ALF2*PSUM+PDFZ(J1M)
      END IF
      RETURN
      END
C ====================== END OF PDFZBETA

C new function from RIF Library
      SUBROUTINE f77betapdf(pdfz,zgrid,zmean,zvari,nz,nzmax)
c
c ======================================================================
c
c      compute beta pdf
c
c      f77betapdf is CALLed by:  species
c
c      f77betapdf CALLs the following subroutines and functions:  
c      none
c
c      author: Hardo Barths, ITM, RWTH Aachen, august 1995
c      modified by: Christian Hasse,ITM,RWTH Aaachen, April 1999
c
c ======================================================================
      IMPLICIT none
c----+------------------------------------------------------------------
      INTEGER nz,nzmax
      DOUBLE PRECISION    pdfz(nzmax),zgrid(nzmax)
      DOUBLE PRECISION    zmean,zvari
c----+------------------------------------------------------------------
      INTEGER n,n1,n1m
      DOUBLE PRECISION    dz
      DOUBLE PRECISION    a,am1,alf1,alf2,b,bm1,gamma,pdfnn2,pdfsum,
     &        psum,smdelt,smzv,smzv1,vbeta,zm,zv
      DOUBLE PRECISION    f77beta,f77gammaln
      DOUBLE PRECISION    min,max,abs
c----+------------------------------------------------------------------
      zm=zmean
      zv=zvari
c
c ---  set limit to zm and zv
c
      smdelt=15.0D-08
      smzv=zm*smdelt
      smzv1=(1.0D+00-zm)*smdelt
      smzv=min(smzv,smzv1)
c
c -------- for the special point zm=0 or 1
c
      IF(zm.LT.smdelt) THEN
        DO 20 n=2,nz
        pdfz(n)=0.0D+00
   20   CONTINUE
        pdfz(1)=1.0D+00
      ELSE IF(zm.GT.(1.0D+00-smdelt)) THEN
        DO 30 n=1,nz-1
        pdfz(n)=0.0D+00
   30   CONTINUE
        pdfz(nz)=1.0D+00
      ELSE 
c
c   --- for zv=0 but zm=(0,1)
c
        IF(zv.LT.smzv) THEN
          DO 40 n=2,nz
          IF(zm.GT.zgrid(n)) GOTO 40
          n1=n
          n1m=n-1
          GOTO 50
  40      CONTINUE
  50      CONTINUE
          DO 60 n=1,nz
          pdfz(n)=0.0D+00
  60      CONTINUE
          pdfz(n1)=(zm-zgrid(n1m))/(zgrid(n1)-zgrid(n1m))
          pdfz(n1m)=1.0D+00-pdfz(n1)
        ELSE
c
c    --- for the other cases in which zm,zv not equal to special point.
c
          gamma=zm*(1.0D+00-zm)/zv-1.0D+00
          IF(gamma.LT.smdelt) THEN
            DO 70 n=2,nz-1
            pdfz(n)=0.0D+00
  70        CONTINUE
            pdfz(1)=1.0D+000-zm
            pdfz(nz)=zm
          ELSE
            a=zm*gamma
            b=gamma-a
            am1=a-1.0D+00
            bm1=b-1.0D+00
c            vbeta=f77beta(a,b)
            vbeta=exp(f77gammaln(a)+f77gammaln(b)-f77gammaln(a+b))
            IF(vbeta.LT.1.D-30) vbeta=1.D-30
            pdfz(1)=0.0D+00
            pdfz(nz)=0.0D+00
            pdfnn2=0.0D+00
            DO 80 n=2,nz-1
              dz=0.5D+00*(zgrid(n+1)-zgrid(n-1))
              pdfz(n)=exp(am1*log(zgrid(n))+bm1
     &                      *log(1.0D+00-zgrid(n)))/vbeta*dz
              pdfnn2=pdfnn2+pdfz(n)
  80        CONTINUE
c
c    ---    for the special point zgrid(1) and zgrid(nz)
c
            IF(am1.GE.0.0D+00.AND.bm1.GE.0.0D+00) THEN
              dz=0.5D+00*(zgrid(2)-zgrid(1))
              pdfz(1)=zgrid(1)**am1 * exp( bm1
     &                *log(1.0D+00-zgrid(1)))/vbeta*dz
              dz=0.5D+00*(zgrid(nz)-zgrid(nz-1))
              pdfz(nz)=exp(am1*log(zgrid(nz)))*
     &                max(0.0D0,1.0D0-zgrid(nz))**bm1*dz/vbeta
            ELSE 
              IF(am1.LT.0.0D+00) pdfz(1)=1.0D+00-pdfnn2
              IF(bm1.LT.0.0D+00) pdfz(nz)=1.0D+00-pdfnn2
              IF(am1.LT.0.0D+00.AND.bm1.LT.0.0D+00) THEN
                pdfz(1)=(1.0D+00-zm)*(1.0D+00-pdfnn2)
                pdfz(nz)=zm*(1.0D+00-pdfnn2)
              ENDIF
            ENDIF
          ENDIF
        ENDIF 
      ENDIF
c
c  ------  check the sum of pdfz(n) which must be equal to 1.0
c
      pdfsum=0.0D+00
      DO 90 n=1,nz
      pdfsum=pdfsum+pdfz(n)
  90  CONTINUE
      IF(abs(pdfsum-1.0D+00).GT.smdelt) THEN
c
c ----  find the interval which zm lies in
c
       DO 100 n=2,nz
          IF(zm.GT.zgrid(n)) GOTO 100
          n1=n
          n1m=n-1
          GOTO 110
 100   CONTINUE
 110   CONTINUE
       psum=1.0D+00-pdfsum
       alf1=(zm-zgrid(n1m))/(zgrid(n1)-zgrid(n1m))
       alf2=(zgrid(n1)-zm)/(zgrid(n1)-zgrid(n1m))
       pdfz(n1)=alf1*psum+pdfz(n1)
       pdfz(n1m)=alf2*psum+pdfz(n1m)
      ENDIF
      RETURN
      END

      DOUBLE PRECISION FUNCTION f77beta(a,b)
c
c ======================================================================
c
c      compute beta function
c
c      f77beta is CALLed by:  f77betapdf
c
c      f77beta CALLs the following subroutines and functions:  
c      gammaln
c
c      author: Hardo Barths, ITM, RWTH Aachen, august 1995
c
c ======================================================================
      IMPLICIT none
c----+------------------------------------------------------------------
      DOUBLE PRECISION   a,b
      DOUBLE PRECISION   f77gammaln,exp
c----+------------------------------------------------------------------
      f77beta=exp(f77gammaln(a)+f77gammaln(b)-f77gammaln(a+b))
      RETURN
      END
      DOUBLE PRECISION FUNCTION f77gammaln(xx)
c
c ======================================================================
c
c      compute gamma function
c
c      f77gammaln is CALLed by:  f77beta
c
c      f77gammaln CALLs the following subroutines and functions:  
c      none
c
c      author: Hardo Barths, ITM, RWTH Aachen, august 1995
c
c ======================================================================
      IMPLICIT none
c----+------------------------------------------------------------------
      DOUBLE PRECISION    xx
c----+------------------------------------------------------------------
      INTEGER j,nx
      PARAMETER (nx=6)
      DOUBLE PRECISION    ser,tmp,x
      DOUBLE PRECISION    cof(nx)
      DOUBLE PRECISION    log
c      DATA    cof /76.18009173,-86.50532033,24.01409822,-1.231739516,
c     &             0.120858003e-02,-0.536382e-05 /
      DATA    cof /76.180091D+00,-86.50532D+00,24.014098D+00,
     &             -1.2317395D+00,0.120858D-02,-0.536382D-05 /
c----+------------------------------------------------------------------
      x=xx-1.0D+00
      tmp=x+5.5D+00
      tmp=tmp-(x+0.5D+00)*log(tmp)
      ser=1.0D+00
      DO 10 j=1,nx
        x=x+1.0D+00
        ser=ser+cof(j)/x
   10 CONTINUE
c      f77gammaln=log(2.50662827465*ser)-tmp
      f77gammaln=log(2.5066282D+00*ser)-tmp
      RETURN
      END
