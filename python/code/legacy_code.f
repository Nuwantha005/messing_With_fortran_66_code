COMMON SRW
      
      DIMENSION AL(21,21), BETA(21,21), CAL(21,21), CBETA(21,21),
     1CURV(21,21), DN(21,21), PRS(21,21), R(21,21), Z(21,21), SM(21,21),
     2SA(21,21), SB(21,21), SC(21,21), SD(21,21), SAL(21,21), 
     3SBETA(21,21),TN(21,21), TT(21,21), WA(21,21), WTR(21,21)

      DIMENSION AB(21), AC(21), AD(21), BA(21), DELBTA(21), DRDM(21),
     1DTDR(21), DTDZ(21), DWMDM(21), DWTDM(21), RH(21), RS(21), ZH(21),
     2ZS(21),THTA(21), WTFL(21), XR(21), XT(21), XZ(21)
      
      INTEGER RUNO, TYPE, BCDP, SRW
      
      RUNO=0
      READ(5,1010) MX, KMX, MR, MZ, W, WT, XN, GAM, AR
      ITNO = 1
      
      RUNO=RUNO+1
      
      WRITE(6,1020) RUNO
      
      WRITE(6,1010) MX, KMX, MR, MZ, W, WT, XN, GAM, AR
      
      READ(5,1010) TYPE, BCDP, SRW, NULL, TEMP, ALM, RHO,
     1TOLER, PLOSS, WTOLER
      WRITE(6,1010) TYPE, BCDP, SRW, NULL, TEMP, ALM, RHO,
     1TOLER, PLOSS, WTOLER
      READ(5,1010) MTHTA, NPRT, ITER, NULL, SFACT, ZSPLIT, BETIN, RB,
     1CORFAC
      WRITE(6,1010) MTHTA, NPRT, ITER, NULL, SFACT, ZSPLIT, BETIN, RB,
     1CORFAC
      READ(5,1030) (ZS(I), I=1, MX)
      
      WRITE(6,1030) (ZS(I), I=1, MX)
      
      READ(5,1030) (ZH(I), I=1, MX)
      
      WRITE(6,1030) (ZH(I), I=1, MX)
      
      READ(5,1030) (RS(I), I=1, MX)
      
      WRITE(6,1030) (RS(I), I=1, MX)
      
      READ(5,1030) (RH(I), I=1, MX)
      
      WRITE(6,1030) (RH(I), I=1, MX)
      
      DO 20 I=1, MX
        ZS(I)=ZS(I)/12.0
        ZH(I)=ZH(I)/12.0
        RS(I)=RS(I)/12.0
        RH(I)=RH(I)/12.0
20    CONTINUE
      
      IF(TYPE.NE.0) GO TO 40
      
      WA(1,1) = WT/RHO/(ZS(1)-ZH(1))/3.14/(RS(1)+RH(1))
      
      DO 30 I=1, MX
        DN(I,KMX)=SQRT((ZS(I)-ZH(I))**2+(RS(I)-RH(I))**2)
        
        DO 30 K=1, KMX
          DN(I,K)=FLOAT(K-1)/FLOAT(KMX-1)*DN(I,KMX)
          WA(I,K)=WA(1,1)
          Z(I,K)=DN(I,K)/DN(I,KMX)*(ZS(I)-ZH(I))+ZH(I)
          R(I,K)=DN(I,K)/DN(I,KMX)*(RS(I)-RH(I))+RH(I)
30    CONTINUE
      
      GO TO 50
      
40    IF(TYPE.NE.1) GO TO 145
      
      CALL BCREAD(DN(1,1), DN(21,21))
      CALL BCREAD(WA(1,1), WA(21,21))
      CALL BCREAD(Z(1,1), Z(21,21))
      CALL BCREAD(R(1,1), R(21,21))
      
50    WRITE(6,1040)
      
      READ(5,1030) (THTA(I), I=1, MTHTA)
      
      WRITE(6,1030) (THTA(I), I=1, MTHTA)
      
      READ(5,1030) (XT(I), I=1, MTHTA)
      
      WRITE(6,1030) (XT(I), I=1, MTHTA)
      
      DO 60 K=1, MR
        READ(5,1030) (TN(I,K), I=1, MZ)
        WRITE(6,1030) (TN(I,K), I=1, MZ)
60    CONTINUE
      
      READ(5,1030) (XZ(I), I=1, MZ)
      WRITE(6,1030) (XZ(I), I=1, MZ)
      READ(5,1030) (XR(I), I=1, MR)
      WRITE(6,1030) (XR(I), I=1, MR)
      
C     END OF INPUT STATEMENTS
      
C     SCALING-CHANGE INCHES TO FEET AND PSI TO LB/SQ FT.
      
      DO 90 K=1, MR
        DO 80 I=1, MZ
          TN(I,K) = TN(I,K)/12.0
80      CONTINUE
        XR(K) = XR(K)/12.0
90    CONTINUE
      
      DO 100 I=1, MZ
        XZ(I) = XZ(I)/12.0
100   CONTINUE
      
      DO 110 I=1, MX
        DO 110 K=1, KMX
          SM(I,K) = 0.0
110   CONTINUE
      
      BA(1) = 0.0
      
      DO 120 K=2, KMX
        BA(K) = FLOAT(K-1)*WT/FLOAT(KMX-1)
120   CONTINUE
      
      DO 130 I=1, MX
        DN(I,1) = 0.0
130   CONTINUE
      
      DO 140 I=1, MTHTA
        XT(I) = XT(I)/12.0
140   CONTINUE
      
      ROOT = SQRT(2.0)
      
145   CONTINUE
      
      TOLER = TOLER/12.0
      RB = RB/12.0
      ZSPLIT = ZSPLIT/12.0
      PLOSS = PLOSS*144.0
      
      CI = SQRT(GAM*AR*TEMP)
      WRITE(6,1050) CI
      
      KMXM1 = KMX-1
      CP = AR*GAM/(GAM-1.0)
      
      EXPON = 1.0/(GAM-1.0)
      
      BETIN = -BETIN/57.29577
      RINLET = (RS(1)+RH(1))/2.0
      CEF = SIN(BETIN)/COS(BETIN)/RINLET/(RINLET-RB)**2
      ERROR = 100000.0
      
C     BEGINNING OF LOOP FOR ITERATIONS
      
150   IF(ITER.EQ.0) WRITE(6,1060) ITNO
      IF(ITER.EQ.0) WRITE(6,1070)
      ERRORI = ERROR
      
      ERROR = 0.0
      
C     START CALCULATION OF PARAMETERS
      
      DO 230 K=1, KMX
C       INITIALIZE.
        DO 160 I=1, MX
          AB(I) = (Z(I,K)-R(I,K))/ROOT
          AC(I) = (Z(I,K)+R(I,K))/ROOT
160     CONTINUE
        
        CALL SPLINE(AB, AC, MX, AL(1,K), CURV(1,K))
        
        DO 170 I=1, MX
          CURV(I,K) = CURV(I,K)/(1.0+AL(I,K)**2)**1.5
          AL(I,K) = ATAN(AL(I,K))-0.785398
          CAL(I,K) = COS(AL(I,K))
          SAL(I,K) = SIN(AL(I,K))
170     CONTINUE
        
        DO 180 I=2, MX
          SM(I,K) = SM(I-1,K)+SQRT((Z(I,K)-Z(I-1,K))**2+
     1    (R(I,K)-R(I-1,K))**2)
180     CONTINUE
        
        CALL SPLDER(XT(1), THTA(1), MTHTA, Z(1,K), MX, DTDZ(1))
        
        DO 220 I=1, MX
          CALL LININT(Z(I,K), R(I,K), XZ, XR, TN, 21, 21, T)
          
          IF(R(I,K).LE.RB) GO TO 200
          DTDR(I) = CEF*(R(I,K)-RB)**2
          GO TO 210
200       DTDR(I) = 0.0
          
210       TQ = R(I,K)*DTDR(I)
          TP = R(I,K)*DTDZ(I)
          TT(I,K) = T*SQRT(1.0+TP*TP)
          BETA(I,K) = ATAN(TP*CAL(I,K)+TQ*SAL(I,K))
          SBETA(I,K) = SIN(BETA(I,K))
          CBETA(I,K) = COS(BETA(I,K))
          
          SA(I,K) = CBETA(I,K)**2*CAL(I,K)*CURV(I,K)-
     1    SBETA(I,K)**2/R(I,K)+SAL(I,K)*CBETA(I,K)*SBETA(I,K)*DTDR(I)
          
          SC(I,K) = -SAL(I,K)*CBETA(I,K)**2*CURV(I,K)
     1    +SAL(I,K)*CBETA(I,K)*SBETA(I,K)*DTDZ(I)
220     CONTINUE
        
        AB(I) = WA(I,K)*CBETA(I,K)
        AC(I) = WA(I,K)*SBETA(I,K)
        
        CALL SPLINE(SM(1,K), AB, MX, DWMDM, AD)
        CALL SPLINE(SM(1,K), AC, MX, DWTDM, AD)
        
        IF((ITER.LE.0).AND.(MOD(K-1,NPRT).EQ.0)) WRITE(6,1080) K
        
        DO 230 I=1, MX
          SB(I,K) = SAL(I,K)*CBETA(I,K)*DWMDM(I)-2.0*W*SBETA(I,K)
     1    +DTDR(I)*R(I,K)*CBETA(I,K)*(DWTDM(I)+2.0*W*SAL(I,K))
          
          SD(I,K) = CAL(I,K)*CBETA(I,K)*DWMDM(I)+DTDZ(I)*
     1    R(I,K)*CBETA(I,K)*(DWTDM(I)+2.0*W*SAL(I,K))
          
          IF((ITER.GT.0).OR.(MOD(K-1,NPRT).NE.0)) GO TO 230
          
          A = AL(I,K)*57.29577
          B = SM(I,K)*12.0
          E = TT(I,K)*12.0
          G = BETA(I,K)*57.29577
          
          WRITE(6,1090) A, CURV(I,K), B, G, E, SA(I,K), SB(I,K), SC(I,K), SD(I,K)
230   CONTINUE
      
C     END OF LOOP - PARAMETER CALCULATION
      
C     CALCULATE BLADE SURFACE VELOCITIES (AFTER CONVERGENCE)
      
      IF(ITER.NE.0) GO TO 260
      
      DO 250 K=1, KMX
        CALL SPLINE(SM(1,K), TT(1,K), MX, DELBTA, AC)
        A = XN
        
        DO 240 I=1, MX
          AB(I) = (R(I,K)*W+WA(I,K)*SBETA(I,K))*
     1    (6.283186*R(I,K)/A-TT(I,K))
240     CONTINUE
        
        CALL SPLINE(SM(1,K), AB, MX, DRDM, AC)
        
        IF(SFACT.LE.1.0) GO TO 245
        A = SFACT*XN
        
        DO 244 I=1, MX
          AB(I) = (R(I,K)*W+WA(I,K)*SBETA(I,K))*
     1    (6.283186*R(I,K)/A-TT(I,K))
244     CONTINUE
        
        CALL SPLINE(SM(1,K), AB, MX, AD, AC)
        
245     DO 250 I=1, MX
          BETAD = BETA(I,K)-DELBTA(I)/2.0
          BETAT = BETAD+DELBTA(I)
          COSBD = COS(BETAD)
          COSBT = COS(BETAT)
          
          IF(Z(I,K).LT.ZSPLIT) DRDM(I) = AD(I)
          
          WTR(I,K) = COSBD*COSBT/(COSBD+COSBT)*(2.0*WA(I,K)/COSBD+R(I,K)
     1    *W*(BETAD-BETAT)/CBETA(I,K)**2+DRDM(I))
250   CONTINUE
      
C     END OF BLADE SURFACE VELOCITY CALCULATIONS
      
C     START CALCULATION OF WEIGHT FLOW VS. DISTANCE FROM HUB
      
260   DO 370 I=1, MX
        IND = 1
        DO 270 K=1, KMX
          AC(K) = DN(I,K)
270     CONTINUE
        
        GO TO 290
        
280     WA(I,1) = 0.5*WA(I,1)
        
290     DO 300 K=2, KMX
          J = K-1
          HR = R(I,K)-R(I,J)
          HZ = Z(I,K)-Z(I,J)
          WAS = WA(I,J)*(1.0+SA(I,J)*HR+SC(I,J)*HZ)
          WASS = WA(I,J)+WAS*(SA(I,K)*HR+SC(I,K)*HZ)
     1    +SB(I,K)*HR+SD(I,K)*HZ
          WA(I,K) = (WAS+WASS)/2.0
300     CONTINUE
        
310     DO 340 K=1, KMX
          TIP = 1.0-(WA(I,K)**2+2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
          
          IF(TIP.LT.0.0) GO TO 280
          
          TPPIP = 1.0-(2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
          
          DENSTY = TIP**EXPON*RHO-(TIP/TPPIP)**EXPON*
     1    PLOSS/AR/TPPIP/TEMP*32.17*SM(I,K)/SM(MX,K)
          
          PRS(I,K) = DENSTY*AR*TIP*TEMP/32.17/144.0
          
          IF(ZS(I).LE.ZH(I)) GO TO 320
          PSI = ATAN((RS(I)-RH(I))/(ZS(I)-ZH(I)))-1.5708
          GO TO 330
          
320       PSI = ATAN((ZH(I)-ZS(I))/(RS(I)-RH(I)))
          
330       WTHRU = WA(I,K)*CBETA(I,K)*COS(PSI-AL(I,K))
          A = XN
          
          IF(Z(I,K).LT.ZSPLIT) A = SFACT*XN
          
          C = 6.283186*R(I,K)-A*TT(I,K)
          AD(K) = DENSTY*WTHRU*C
340     CONTINUE
        
        CALL INTGRL(AC(1), AD(1), KMX, WTFL(1))
        
        IF(ABS(WT-WTFL(KMX)).LE.WTOLER) GO TO 350
        
        CALL CONTIN(WA(I,1), WTFL(KMX), IND, I, WT)
        
        IF(IND.NE.0) GO TO 290
        
350     CALL SPLINT(WTFL, AC, KMX, BA, KMX, AB)
        
        DO 360 K=1, KMX
          DELTA = ABS(AB(K)-DN(I,K))
          DN(I,K) = (1.0-CORFAC)*DN(I,K)+CORFAC*AB(K)
          IF(DELTA.GT.ERROR) ERROR = DELTA
360     CONTINUE
        
370   CONTINUE
      
C     END OF LOOP - WEIGHT FLOW CALCULATION
      
C     CALCULATE STREAMLINE COORDINATES FOR NEXT ITERATION
      
      DO 380 K=2, KMXM1
        DO 380 I=1, MX
          Z(I,K) = DN(I,K)/DN(I,KMX)*(ZS(I)-ZH(I))+ZH(I)
          R(I,K) = DN(I,K)/DN(I,KMX)*(RS(I)-RH(I))+RH(I)
380     CONTINUE
390   CONTINUE
      
      IF((ERROR.GE.ERRORI).OR.(ERROR.LE.TOLER)) ITER = ITER-1
      
      IF(ITER.GT.0) GO TO 410
      
      WRITE(6,1100)
      
      DO 400 K=1, KMX, NPRT
        WRITE(6,1080) K
        
        DO 390 I=1, MX
          AB(I) = (Z(I,K)-R(I,K))/ROOT
          AC(I) = (Z(I,K)+R(I,K))/ROOT
          
          CALL SPLINE(AB, AC, MX, AL(I,K), CURV(I,K))
          
          DO 400 I=1, MX
            CURV(I,K) = CURV(I,K)/(1.0+AL(I,K)**2)**1.5
            
            A = DN(I,K)*12.0
            B = Z(I,K)*12.0
            D = R(I,K)*12.0
            
            WRITE(6,1110) A, B, D, WA(I,K), PRS(I,K),
     1      WTR(I,K), CURV(I,K)
400   CONTINUE
      
      WRITE(6,1130)
      A = ERROR*12.0
      WRITE(6,1120) ITNO, A
      ITNO = ITNO+1
      
410   IF(ITER.GE.0) GO TO 150
      
      IF(BCDP.NE.1) GO TO 10
      
      CALL BCDUMP(DN(1,1), DN(21,21))
      CALL BCDUMP(WA(1,1), WA(21,21))
      CALL BCDUMP(Z(1,1), Z(21,21))
      CALL BCDUMP(R(1,1), R(21,21))
      
      GO TO 10
      
1010  FORMAT(4I5, 5F10.4)
1020  FORMAT(8H1RUN NO.I3, 19X, 25HINPUT DATA CARD LISTING     )
1030  FORMAT(7F10.4)
1040  FORMAT(10X, 24HBCD CARDS FOR DN,WA,Z,R- )
1050  FORMAT(36HK STAGE SPEED OF SOUND AT INLET ,F9.2)
1060  FORMAT(///5X, 13HITERATION NO.I3)
1070  FORMAT(1H, 6X5HAL, 9X5HRC, 9X5HSM, 9X5HBETA, 9X5HTT, 9X5HSA,
     1 1X5HSB, 9X5HSC, 9X5HSD     )
1080  FORMAT(2X, 10HSTREAMLINEI3)
1090  FORMAT(9F14.6)
1100  FORMAT(1H1, 9X5HDN, 15X5HZ, 15X5HR, 15X5HWA, 15X5HPRS, 14X3HW
     1TR, 14X3HRC     )
1110  FORMAT(6F19.6, F18.6)
1120  FORMAT(18H ITERATION NO. I3, 10X, 24HMAX. STREAMLINE CHANGE = ,
     1F10.6)
1130  FORMAT(1H1)
      END

      SUBROUTINE SPLINE(X, Y, N, SLOPE, EM)
      DIMENSION X(50), Y(50), S(50), A(50), B(50), C(50), F(50), W(50), 
     1           SB(50),G(50), EM(50), SLOPE(50)
      COMMON Q
      INTEGER Q
      DO 10 I=2,N
   10 S(I) = X(I) - X(I-1)
      NO = N-1
      DO 20 I=2,NO
      A(I) = S(I)/6.
      B(I) = (S(I) + S(I+1))/3.
      C(I) = S(I+1)/6.
   20 F(I) = (Y(I+1) - Y(I))/S(I+1) - (Y(I) - Y(I-1))/S(I)
      A(N) = -.5
      B(1) = 1.
      B(N) = 1.
      C(1) = -.5
      F(1) = 0.
      F(N) = 0.
      W(1) = B(1)
      SB(1) = C(1)/W(1)
      G(1) = 0.
      DO 30 I=2,N
      W(I) = B(I) - A(I)*SB(I-1)
      SB(I) = C(I)/W(I)
   30 G(I) = (F(I) - A(I)*G(I-1))/W(I)
      EM(N) = G(N)
      DO 40 I=2,N
      K = N+1-I
   40 EM(K) = G(K) - SB(K)*EM(K+1)
      SLOPE(1) = -S(2)/6.*(2.*EM(1) + EM(2)) + (Y(2) - Y(1))/S(2)
      DO 50 I=2,N
   50 SLOPE(I) = S(I)/6.*(2.*EM(I) + EM(I-1)) + (Y(I) - Y(I-1))/S(I)
      IF(Q.EQ.13) WRITE(6,100) N, (X(I), Y(I), SLOPE(I), EM(I), I=1,N)
  100 FORMAT(2X, 15HNO. OF POINTS , I3/10X, 5HX    , 15X, 5HY    ,
     1       15X, 5HSLOPE, 15X, 5HEM   /(4F20.8))
      RETURN
      END

      SUBROUTINE SPLINT(X, Y, N, Z, MAX, YINT)
      DIMENSION X(50), Y(50), S(50), A(50), B(50), C(50), F(50), W(50), 
     1           SB(50),G(50), EM(50), Z(50), YINT(50)
      COMMON Q
      INTEGER Q
      DO 10 I=2,N
   10 S(I) = X(I) - X(I-1)
      NO = N-1
      DO 20 I=2,NO
      A(I) = S(I)/6.0
      B(I) = (S(I) + S(I+1))/3.0
      C(I) = S(I+1)/6.0
   20 F(I) = (Y(I+1) - Y(I))/S(I+1) - (Y(I) - Y(I-1))/S(I)
      A(N) = -.5
      B(1) = 1.0
      B(N) = 1.0
      C(1) = -.5
      F(1) = 0.0
      F(N) = 0.0
      W(1) = B(1)
      SB(1) = C(1)/W(1)
      G(1) = 0.0
      DO 30 I=2,N
      W(I) = B(I) - A(I)*SB(I-1)
      SB(I) = C(I)/W(I)
   30 G(I) = (F(I) - A(I)*G(I-1))/W(I)
      EM(N) = G(N)
      DO 40 I=2,N
      K = N+1-I
   40 EM(K) = G(K) - SB(K)*EM(K+1)
      DO 90 I=1,MAX
      K = 2
      IF(Z(I) - X(1)) 60, 50, 70
   50 YINT(I) = Y(1)
      GO TO 90
   60 IF(Z(I).LT.(1.1*X(1) - 0.1*X(2))) WRITE(6, 1000) Z(I)
      GO TO 85
 1000 FORMAT(17H OUT OF RANGE Z =, F10.6)
   65 IF(Z(I).GT.(1.1*X(N) - 0.1*X(N-1))) WRITE(6, 1000) Z(I)
      K = N
      GO TO 85
   70 IF(Z(I) - X(K)) 85, 75, 80
   75 YINT(I) = Y(K)
      GO TO 90
   80 K = K+1
      IF(K-N) 70, 70, 65
   85 YINT(I) = EM(K-1)*(X(K) - Z(I))**3/6./S(K) + EM(K)*(Z(I) - X(K-1))
     1          **3/6./S(K) + (Y(K)/S(K) - EM(K)*S(K)/6.)*(Z(I) - 
     2          X(K-1)) + (Y(K-1)/S(K) - EM(K-1)*S(K)/6.)*(X(K) - Z(I))
   90 CONTINUE
      IF(Q.EQ.16) WRITE(6, 1010) N,MAX,(X(I),Y(I),Z(I),YINT(I), I=1,N)
 1010 FORMAT(2X, 21HNO. OF POINTS GIVEN =, I3, 30H, NO. OF INTERPOLATED 
     1       POINTS =, I3, /10X, 5HX    , 15X, 5HY    , 12X, 
     2       11HX-INTERPOL.,9X, 11HY-INTERPOL./(4E20.8))
      RETURN
      END

      SUBROUTINE SPLDER(X, Y, N, Z, MAX, DYDX)
      DIMENSION X(50), Y(50), S(50), A(50), B(50), C(50), F(50), W(50),
     1           SB(50),G(50), EM(50), Z(50), DYDX(50)
      DO 10 I=2,N
   10 S(I) = X(I) - X(I-1)
      NO = N-1
      DO 20 I=2,NO
      A(I) = S(I)/6.0
      B(I) = (S(I) + S(I+1))/3.0
      C(I) = S(I+1)/6.0
   20 F(I) = (Y(I+1) - Y(I))/S(I+1) - (Y(I) - Y(I-1))/S(I)
      A(N) = -.5
      B(1) = 1.0
      B(N) = 1.0
      C(1) = -.5
      F(1) = 0.0
      F(N) = 0.0
      W(1) = B(1)
      SB(1) = C(1)/W(1)
      G(1) = 0.0
      DO 30 I=2,N
      W(I) = B(I) - A(I)*SB(I-1)
      SB(I) = C(I)/W(I)
   30 G(I) = (F(I) - A(I)*G(I-1))/W(I)
      EM(N) = G(N)
      DO 40 I=2,N
      K = N+1-I
   40 EM(K) = G(K) - SB(K)*EM(K+1)
      DO 90 I=1,MAX
      K = 2
      IF(Z(I) - X(1)) 60, 70, 70
   60 WRITE(6, 1000) Z(I)
 1000 FORMAT(17H OUT OF RANGE Z =, F10.6)
      GO TO 85
   65 WRITE(6, 1000) Z(I)
      K = N
      GO TO 85
   70 IF(Z(I) - X(K)) 85, 85, 80
   80 K = K+1
      IF(K-N) 70, 70, 65
   85 DYDX(I) = -EM(K-1)*(X(K) - Z(I))**2/2.0/S(K) + EM(K)*(Z(I) -
     1           X(K-1))**2/2.0/S(K) + (Y(K) - Y(K-1))/S(K) - 
     2           (EM(K) - EM(K-1))*S(K)/6. 
   90 CONTINUE
      RETURN
      END

      SUBROUTINE INTGRL(X, Y, N, SUM)
      DIMENSION X(50), Y(50), S(50), A(50), B(50), C(50), F(50), W(50),
     1           SB(50),G(50), EM(50), SUM(50)
      DO 10 I=2,N
   10 S(I) = X(I) - X(I-1)
      NO = N-1
      DO 20 I=2,NO
      A(I) = S(I)/6.0
      B(I) = (S(I) + S(I+1))/3.0
      C(I) = S(I+1)/6.0
   20 F(I) = (Y(I+1) - Y(I))/S(I+1) - (Y(I) - Y(I-1))/S(I)
      A(N) = -.5
      B(1) = 1.0
      B(N) = 1.0
      C(1) = -.5
      F(1) = 0.0
      F(N) = 0.0
      W(1) = B(1)
      SB(1) = C(1)/W(1)
      G(1) = 0.0
      DO 30 I=2,N
      W(I) = B(I) - A(I)*SB(I-1)
      SB(I) = C(I)/W(I)
   30 G(I) = (F(I) - A(I)*G(I-1))/W(I)
      EM(N) = G(N)
      DO 40 I=2,N
      K = N+1-I
   40 EM(K) = G(K) - SB(K)*EM(K+1)
      SUM(1) = 0.0
      DO 50 K=2,N
   50 SUM(K) = SUM(K-1) + S(K)*(Y(K) + Y(K-1))/2.0 - S(K)**3*(EM(K) 
     1+ EM(K-1))/24.0 
      RETURN
      END
C     -----------------------------------------------------------------------
      SUBROUTINE LININT(X1, Y1, X, Y, TN, MX, MY, F)
      COMMON K
  DIMENSION X(MX), Y(MY), TN(MX, MY)
  
  DO 10 J3=1, MX
10 IF(X1.LE.X(J3)) GO TO 20
  J3=MX
20 DO 30 J4=1, MY
30 IF(Y1.LE.Y(J4)) GO TO 40
  J4=MY
40 J1=J3-1
  J2=J4-1
  EPS1=(X1-X(J1))/(X(J3)-X(J1))
  EPS2=(Y1-Y(J2))/(Y(J4)-Y(J2))
  EPS3=1.-EPS1
  EPS4=1.-EPS2
  F=TN(J1,J2)*EPS3*EPS4+TN(J3,J2)*EPS1*EPS4+TN(J1,J4)*EPS2*EPS3+
     1TN(J3,J4)*EPS1*EPS2
  IF(K.EQ.14) WRITE(6,1) X1, Y1, F, J1, J2, EPS1, EPS2
1 FORMAT(8H LININT,3F10.5,2I3,2F10.5)
  K=0
  RETURN
END


SUBROUTINE CONTIN(WA, WTFL, IND, I, WT)
  DIMENSION SPEED(3), WEIGHT(3)
135 GO TO (140,150,210,270,370), IND
140 SPEED(1) = WA
  WEIGHT(1) = WTFL
  WA = WT/WTFL*WA
  IND = 2
  RETURN
150 IF((WTFL-WEIGHT(1))/(WA-SPEED(1))) 180,180,160
160 SPEED(2) = WA
  WA = (WT-WTFL)/(WTFL-WEIGHT(1))*(WA-SPEED(1))+WA
  IF(ABS(WA-SPEED(2))-100.0) 166,166,161
161 IF(WA-SPEED(2)) 163,163,162
162 WA = SPEED(2)+100.0
  GO TO 166
163 WA = SPEED(2)-100.0
166 SPEED(1) = SPEED(2)
  WEIGHT(1) = WTFL
  RETURN
170 WRITE(6,1000) I, WTFL
  IND = 6
  RETURN
180 IND = 3
  IF(WTFL.GE.WT) GO TO 140
  IF(SPEED(1)-WA) 190,200,200
190 SPEED(2) = SPEED(1)
  SPEED(1) = 2.0*SPEED(1)-WA
  SPEED(3) = WA
  WEIGHT(2) = WEIGHT(1)
  WEIGHT(3) = WTFL
  WA = SPEED(1)
  RETURN
200 SPEED(2) = WA
  SPEED(3) = SPEED(1)
  SPEED(1) = 2.0*WA-SPEED(1)
  WEIGHT(2) = WTFL
  WEIGHT(3) = WEIGHT(1)
  WA = SPEED(1)
  RETURN
210 WEIGHT(1) = WTFL
  IF(WTFL.GE.WT) GO TO 140
  IF(WEIGHT(1)-WEIGHT(2)) 230,380,220
220 WEIGHT(3) = WEIGHT(2)
  WEIGHT(2) = WEIGHT(1)
  SPEED(3) = SPEED(2)
  SPEED(2) = SPEED(1)
  SPEED(1) = 2.0*SPEED(2)-SPEED(3)
  WA = SPEED(1)
  RETURN
230 IF(SPEED(3)-SPEED(1)-10.0) 170,170,240
240 IND = 4
245 IF(WEIGHT(3)-WEIGHT(1)) 260,260,250
250 WA = (SPEED(1)+SPEED(2))/2.0
  RETURN
260 WA = (SPEED(3)+SPEED(2))/2.0
  RETURN
270 IF(SPEED(3)-SPEED(1)-10.0) 170,170,280
280 IF(WTFL-WEIGHT(2)) 320,350,290
290 IF(WA-SPEED(2)) 310,300,300
300 SPEED(1) = SPEED(2)
  SPEED(2) = WA
  WEIGHT(1) = WEIGHT(2)
  WEIGHT(2) = WTFL
  GO TO 245
310 SPEED(3) = SPEED(2)
  SPEED(2) = WA
  WEIGHT(3) = WEIGHT(2)
  WEIGHT(2) = WTFL
  GO TO 245
320 IF(WA-SPEED(2)) 340,330,330
330 WEIGHT(3) = WTFL
  SPEED(3) = WA
  GO TO 245
340 WEIGHT(1) = WTFL
  SPEED(1) = WA
  GO TO 245
350 IND = 5
  IF(WA-SPEED(2)) 380,360,360
360 SPEED(1) = SPEED(2)
  WEIGHT(1) = WEIGHT(2)
  SPEED(2) = (SPEED(1)+SPEED(3))/2.0
  WA = SPEED(2)
  RETURN
370 IND = 4
  WEIGHT(2) = WTFL
  WA = (SPEED(1)+SPEED(2))/2.0
  RETURN
380 IND = 5
  WEIGHT(3) = WEIGHT(2)
  SPEED(3) = SPEED(2)
  SPEED(2) = (SPEED(1)+SPEED(3))/2.0
  WA = SPEED(2)
  RETURN
1000 FORMAT(/12H FIXED LINE ,I2,12H , MAX WT = ,F10.6)
END