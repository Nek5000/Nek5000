      SUBROUTINE HLEGN(H,ZH,ZLEG,NZZ)
      REAL H(1)
      REAL ZLEG(1)
      parameter (nxm=64)
      REAL    ALFAN,PNZZ(NXM)
      SAVE    ALFAN,PNZZ
      INTEGER NZ0
      SAVE    NZ0
      DATA    NZ0 /0/
C
      IF (NZZ.NE.NZ0) THEN
         NZ0 = NZZ
         ALFAN = FLOAT(NZZ*(NZZ-1))
         DO 10 I=1,NZZ
            PNZZ(I) = ALFAN*PNZ(Zleg(I),NZZ)
   10    CONTINUE
      ENDIF
C
      IF (NZZ.EQ.1) THEN
         H(1)=1.0
         RETURN
      ENDIF
C
      EPS = 1.E-5
      PD  =  - (1.-ZH*ZH)*PND(ZH,NZZ)
C
      DO 20 I=1,NZZ
         DZ = ZH - ZPTS(I)
         IF (ABS(DZ) .LT. EPS) THEN
            H(I) = 1.0
         ELSE
            H(I) = PD/(PNZZ(I)*DZ)
         ENDIF
   20 CONTINUE
C
      RETURN
      END
