COMMON SRW
      
      DIMENSION AL(21,21), BETA(21,21), CAL(21,21), CBETA(21,21),
     1 CURV(21,21), DN(21,21), PRS(21,21), R(21,21), Z(21,21), SM(21,21),
     2 SA(21,21), SB(21,21), SC(21,21), SD(21,21), SAL(21,21), SBETA(21,21),
     3 TN(21,21), TT(21,21), WA(21,21), WTR(21,21)
      
      DIMENSION AB(21), AC(21), AD(21), BA(21), DELBTA(21), DRDM(21),
     1 DTDR(21), DTDZ(21), DWMDM(21), DWTDM(21), RH(21), RS(21), ZH(21), ZS(21),
     2 THTA(21), WTFL(21), XR(21), XT(21), XZ(21)
      
      INTEGER RUNO, TYPE, BCDP, SRW
      
      RUNO=0
      READ(5,1010) MX, KMX, MR, MZ, W, WT, XN, GAM, AR
      ITNO = 1
      
      RUNO=RUNO+1
      
      WRITE(6,1020) RUNO
      
      WRITE(6,1010) MX, KMX, MR, MZ, W, WT, XN, GAM, AR
      
      READ(5,1010) TYPE, BCDP, SRW, NULL, TEMP, ALM, RHO, TOLER, PLOSS, WTOLER
      WRITE(6,1010) TYPE, BCDP, SRW, NULL, TEMP, ALM, RHO, TOLER, PLOSS, WTOLER
      READ(5,1010) MTHTA, NPRT, ITER, NULL, SFACT, ZSPLIT, BETIN, RB, CORFAC
      WRITE(6,1010) MTHTA, NPRT, ITER, NULL, SFACT, ZSPLIT, BETIN, RB, CORFAC
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
          SM(I,K) = SM(I-1,K)+SQRT((Z(I,K)-Z(I-1,K))**2+(R(I,K)-R(I-1,K))**2)
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
          
          SA(I,K) = CBETA(I,K)**2*CAL(I,K)*CURV(I,K)-SBETA(I,K)**2/R(I,K)+
     1    SAL(I,K)*CBETA(I,K)*SBETA(I,K)*DTDR(I)
          
          SC(I,K) = -SAL(I,K)*CBETA(I,K)**2*CURV(I,K)+SAL(I,K)*CBETA(I,K)*
     1    SBETA(I,K)*DTDZ(I)
220     CONTINUE
        
        AB(I) = WA(I,K)*CBETA(I,K)
        AC(I) = WA(I,K)*SBETA(I,K)
        
        CALL SPLINE(SM(1,K), AB, MX, DWMDM, AD)
        CALL SPLINE(SM(1,K), AC, MX, DWTDM, AD)
        
        IF((ITER.LE.0).AND.(MOD(K-1,NPRT).EQ.0)) WRITE(6,1080) K
        
        DO 230 I=1, MX
          SB(I,K) = SAL(I,K)*CBETA(I,K)*DWMDM(I)-2.0*W*SBETA(I,K)+DTDR(I)*
     1    R(I,K)*CBETA(I,K)*(DWTDM(I)+2.0*W*SAL(I,K))
          
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
          AB(I) = (R(I,K)*W+WA(I,K)*SBETA(I,K))*(6.283186*R(I,K)/A-TT(I,K))
240     CONTINUE
        
        CALL SPLINE(SM(1,K), AB, MX, DRDM, AC)
        
        IF(SFACT.LE.1.0) GO TO 245
        A = SFACT*XN
        
        DO 244 I=1, MX
          AB(I) = (R(I,K)*W+WA(I,K)*SBETA(I,K))*(6.283186*R(I,K)/A-TT(I,K))
244     CONTINUE
        
        CALL SPLINE(SM(1,K), AB, MX, AD, AC)
        
245     DO 250 I=1, MX
          BETAD = BETA(I,K)-DELBTA(I)/2.0
          BETAT = BETAD+DELBTA(I)
          COSBD = COS(BETAD)
          COSBT = COS(BETAT)
          
          IF(Z(I,K).LT.ZSPLIT) DRDM(I) = AD(I)
          
          WTR(I,K) = COSBD*COSBT/(COSBD+COSBT)*(2.0*WA(I,K)/COSBD+R(I,K)*W*
     1    (BETAD-BETAT)/CBETA(I,K)**2+DRDM(I))
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
          WASS = WA(I,J)+WAS*(SA(I,K)*HR+SC(I,K)*HZ)+SB(I,K)*HR+SD(I,K)*HZ
          WA(I,K) = (WAS+WASS)/2.0
300     CONTINUE
        
310     DO 340 K=1, KMX
          TIP = 1.0-(WA(I,K)**2+2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
          
          IF(TIP.LT.0.0) GO TO 280
          
          TPPIP = 1.0-(2.0*W*ALM-(W*R(I,K))**2)/2.0/CP/TEMP
          
          DENSTY = TIP**EXPON*RHO-(TIP/TPPIP)**EXPON*PLOSS/AR/TPPIP/TEMP*
     1    32.17*SM(I,K)/SM(MX,K)
          
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
      
      DO 390 K=2, KMXM1
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
            
            WRITE(6,1110) A, B, D, WA(I,K), PRS(I,K), WTR(I,K), CURV(I,K)
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
1020  FORMAT(8H1RUN NO.I3, 19X, 25HINPUT DATA CARD LISTING )
1030  FORMAT(7F10.4)
1040  FORMAT(10X, 24HBCD CARDS FOR DN,WA,Z,R- )
1050  FORMAT(36HK STAGE SPEED OF SOUND AT INLET ,F9.2)
1060  FORMAT(///5X, 13HITERATION NO.I3)
1070  FORMAT(1H, 6X5HAL, 9X5HRC, 9X5HSM, 9X5HBETA, 9X5HTT, 9X5HSA,
     1 1X5HSB, 9X5HSC, 9X5HSD)
1080  FORMAT(2X, 10HSTREAMLINEI3)
1090  FORMAT(9F14.6)
1100  FORMAT(1H1, 9X5HDN, 15X5HZ, 15X5HR, 15X5HWA, 15X5HPRS, 14X3HW
     1TR, 14X3HRC)
1110  FORMAT(6F19.6, F18.6)
1120  FORMAT(18H ITERATION NO. I3, 10X, 24HMAX. STREAMLINE CHANGE = ,
     1F10.6)
1130  FORMAT(1H1)
      END

