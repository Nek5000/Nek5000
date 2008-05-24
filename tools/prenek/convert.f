      SUBROUTINE VCNVERT(A,N)
      REAL*4 A(1)
      DO 100 I=1,N
         CALL CONVERT(A(I))
  100 CONTINUE
      RETURN
      END
      SUBROUTINE CONVERT(T)
C
      CHARACTER*12 VC
      CHARACTER*1  V1(12)
      EQUIVALENCE (V1,VC)
      CHARACTER*1  CSIGN
C
      REAL*4 W
      CHARACTER*4  WC
      CHARACTER*1  W1(4)
      EQUIVALENCE (W1,W)
      EQUIVALENCE (W1,WC)
C
      CHARACTER*1 ALPH64(0:63)
      SAVE        ALPH64
      DATA        ALPH64 
     $  /'1','2','3','4','5','6','7','8','9','0'
     $   ,'a','b','c','d','e','f','g','h','i','j'
     $   ,'k','l','m','n','o','p','q','r','s','t'
     $   ,'u','v','w','x','y','z'
     $   ,'A','B','C','D','E','F','G','H','I','J'
     $   ,'K','L','M','N','O','P','Q','R','S','T'
     $   ,'U','V','W','X','Y','Z','+','-'/
C
C     Find out the usual decimal format for T
C
      WRITE(VC,10) T
   10 FORMAT(E12.5)
C
C     Begin converting the mantissa to base 64
C
      READ(VC,11) MANTIS
   11 FORMAT(3X,I5)
C
C     Sign?
C
      READ(VC,12) CSIGN
   12 FORMAT(A1)
      IF (CSIGN.EQ.'-') MANTIS=-MANTIS
      MANTIS=MANTIS+131072
C
C     ONES,TENS, HUNDREDS
C
      IONE=MOD(MANTIS,64)
      ITMP=MANTIS/64
      ITEN=MOD(ITMP,64)
      ITMP=ITMP/64
      IHUN=MOD(ITMP,64)
C
C     Exponent
C
      READ(VC,21) IEXP
   21 FORMAT(9X,I3)
C     We assume that the exponent is bounded by 31.
      IEXP=IEXP+31
C
C     Compute alpha equivalent
C
      W1(1)=ALPH64(IHUN)
      W1(2)=ALPH64(ITEN)
      W1(3)=ALPH64(IONE)
      W1(4)=ALPH64(IEXP)
C
C     Convert the input value
      T=W
      RETURN
      END
      SUBROUTINE VRNVERT(A,N)
      REAL*4 A(1)
      DO 100 I=1,N
         A(I)=REVERT(A(I))
  100 CONTINUE
      RETURN
      END
      FUNCTION REVERT(T)
C
      REAL*4 REVERT
      CHARACTER*4 T
      CHARACTER*12 VC
      CHARACTER*1  V1(12)
      EQUIVALENCE (V1,VC)
C
      INTEGER IC
      CHARACTER*1  C(4)
      EQUIVALENCE (IC,C)
C
      CHARACTER*4  WC
      CHARACTER*1  W1(4)
      EQUIVALENCE (W1,WC)
C
      CHARACTER*1 ALPSRT(4,0:63)
      CHARACTER*4 ALPSR4(0:63)
      SAVE        ALPH64,ALPSRT
      INTEGER INDEX(0:63),INTALP(0:63),LOG64
      SAVE    INDEX,INTALP,LOG64
      EQUIVALENCE (INTALP,ALPSRT)
      EQUIVALENCE (INTALP,ALPSR4)
C
      CHARACTER*1 ALPH64(0:63)
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
      DATA        ALPH64 
     $  /'1','2','3','4','5','6','7','8','9','0'
     $   ,'a','b','c','d','e','f','g','h','i','j'
     $   ,'k','l','m','n','o','p','q','r','s','t'
     $   ,'u','v','w','x','y','z'
     $   ,'A','B','C','D','E','F','G','H','I','J'
     $   ,'K','L','M','N','O','P','Q','R','S','T'
     $   ,'U','V','W','X','Y','Z','+','-'/
C
C     On first call, set up sorted list of Alpha values so
C     that we can use binary chop to evaluate the integer
C     values.
C
      IF (ICALLD.EQ.0) THEN
         ICALLD=1
         CALL BLANK(ALPSRT,256)
         DO 10 I=0,63
            ALPSRT(4,I)=ALPH64(I)
   10    CONTINUE
         CALL ISORT(INTALP,INDEX,64)
         DO 20 I=0,63
   20    INDEX(I)=INDEX(I)-1
         LOG64=7
      ENDIF
C
C     Copy T to W
C
      WC=T
      CALL BLANK(C,4)
C
C     Begin divide and conquer search
C
      ILO=0
      IHI=63
      C(4)=W1(1)
      DO 100 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IHUN=INDEX(I)
            GOTO 101
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  100 CONTINUE
  101 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(2)
      DO 200 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            ITEN=INDEX(I)
            GOTO 201
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  200 CONTINUE
  201 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(3)
      DO 300 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IONE=INDEX(I)
            GOTO 301
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  300 CONTINUE
  301 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(4)
      DO 400 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IEXP=INDEX(I)
            GOTO 401
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  400 CONTINUE
  401 CONTINUE
C
c     MANTIS=4096*IHUN+64*ITEN+IONE
      MANTIS=4096*IHUN+64*ITEN+IONE-131072
      IEXP  =IEXP-31
C
C     kludge (5-12-91)
c     IF (ABS(MANTIS).GT.99999.OR.ABS(IEXP).GT.25) THEN
c        MANTIS=0
c        IEXP=0
c     ENDIF
      IF (MANTIS.LT.0) THEN
         MANTIS=ABS(MANTIS)
         WRITE(VC,1001) MANTIS,IEXP
      ELSE
         WRITE(VC,1002) MANTIS,IEXP
      ENDIF
 1001 FORMAT('-0.',I5.5,'E',I3.2)
 1002 FORMAT(' 0.',I5.5,'E',I3.2)
C
      READ(VC,2001,ERR=3000) W
 2001 FORMAT(E12.5)
 2002 FORMAT(2X,A12)
C
      REVERT=W
 2003 FORMAT(2X,A4,F30.12)
      RETURN
C
 3000 CONTINUE
      IF (MANTIS.LT.0) THEN
         MANTIS=ABS(MANTIS)
         WRITE(6,1001) MANTIS,IEXP
      ELSE
         WRITE(6,1002) MANTIS,IEXP
      ENDIF
      WRITE(6,3001) T,MANTIS,IEXP,Ihun,Iten,Ione
 3001 FORMAT(' Problem in REVERT: ',A4,5I10)
      stop
      END
