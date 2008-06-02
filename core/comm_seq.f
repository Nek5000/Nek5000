c----------------------------------------------------------------------
      SUBROUTINE INIPROC
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      COMMON /CTMP0/ IPG(LP),IPG1(LP),ITMP(LP)
      INTEGER CD1
C
C------------------------------------------------------
C     Initialize cube process parameters for iPSC
C------------------------------------------------------
C
      call init_nek_comm (nid,np,wdsize)
      flag_gs_init = 0

      CD  = LOG2(NP)
      PID = 0
      NODE= NID+1
C
C     Check cube dimension 
C
      IF (NP.GT.LP) THEN
         WRITE(6,*) 
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         WRITE(6,*) 
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         WRITE(6,*) 
     $   'Aborting in routine INIPROC.'
      ENDIF
C
c     These flags added for native 64-bit machines  (pff 10/1/98)
c
      ifdblas = .false.
      if (wdsize.eq.8) ifdblas = .true.
c
c     This is to determine the message length of integers
c
      isize = 4
c
c     On the Cray:   ifdblas = .false., since sblas = 64 bit
c                    isize  is 8
c
c     ifdblas = .false.
c     isize   = 8
c
c
      MANAGER=0
      ALLNODES=-1
      NULLPID=0
      NODE0=0
C------------------------------------------------------
C     GOP types
C------------------------------------------------------
C
      TYPE  = 99929917
      TYPER = 98917722
C
C------------------------------------------------------
C     Ring pass data
C------------------------------------------------------
C
C     Initialize parallel processor parameters (for serial code)
C
      CALL GRAY(IPG,IPG1,ITMP,NP)
      IGNODE=IPG1(NODE)-1
      DO 100 IP=1,NP
         IGNODE=IGNODE+1
         IF (IGNODE.GT.NP) IGNODE=IGNODE-NP
         IPRING(IP)=IPG(IGNODE)-1
  100 CONTINUE
      LFTNBR=IPRING(NP)
      RGTNBR=IPRING(2)
C
C------------------------------------------------------
C     All done
C------------------------------------------------------
C
      RETURN
      END

      subroutine init_nek_comm(nido,npo,wdsize)
C------------------------------------------------------
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer wdsize
C
      nekcomm = 1
      nekgroup = 1
      nid  = mynode()
      np   = numnodes()
      nido = nid
      npo  = np
c
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) wdsize=8
      nekreal = mpi_real
      if (wdsize.eq.8) nekreal = mpi_double_precision
c
      return
      end
c
c-----------------------------------------------------------------------
      SUBROUTINE RRING(WORK,X,N,IRG)
      DIMENSION X(N),WORK(N)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
C
      IF (IRG.EQ.1) THEN
         CALL COPY(WORK,X,N)
      ELSE
         LEN   = WDSIZE*N
         LEN2  = 2*LEN
         ITYPE = 6000+NID
         CALL CSEND(ITYPE,WORK,LEN,LFTNBR,NULLPID)
         JTYPE = 6000+RGTNBR
         CALL CRECV(JTYPE,WORK,LEN2)
      ENDIF
      RETURN
      END
      SUBROUTINE CRING(WORK,X,N,IRG)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      DIMENSION X(N),WORK(N)
      CHARACTER*1 X,WORK
C
      PARAMETER (LTOT2=LX1*LY1*LZ1*LELT*2)
      PARAMETER (LTOT8=LTOT2*4)
      COMMON /CTMP0/ CW1(LTOT8)
      CHARACTER*1    CW1
      DIMENSION      RW1(LTOT2)
      EQUIVALENCE   (RW1,CW1)
C
      IF (IRG.EQ.1) THEN
         CALL CHCOPY(WORK,X,N)
      ELSE
C
C        enough work space?
         IF (N.GT.LTOT8) THEN
            WRITE(6,100) NID,N,LTOT8
  100       FORMAT(2X,I5,
     $     'WARNING: In routine CRING, not enough work space.'
     $     ,/,2X,'Required # of words:',I7,' Supplied:',I7)
            CALL EXITT
         ENDIF
C
         CALL CHCOPY(CW1,WORK,N)
         LEN   = N
         LEN2  = 2*LEN
         ITYPE = 6000+NID
         CALL CSEND(ITYPE,RW1,LEN,LFTNBR,NULLPID)
         JTYPE = 6000+RGTNBR
         CALL CRECV(JTYPE,RW1,LEN2)
         ICNT = INFOCOUNT()
         CALL CHCOPY(WORK,CW1,ICNT)
      ENDIF
      RETURN
      END
      SUBROUTINE IRING(WORK,X,N,IRG)
      INTEGER X(N),WORK(N)
C
C     Pass data X around a ring network, and receive corresponding
C     data into WORK.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
C
      IF (IRG.EQ.1) THEN
         CALL ICOPY(WORK,X,N)
      ELSE
         LEN   = ISIZE*N
         LEN2  = 2*LEN
         ITYPE = 6000+NID
         CALL CSEND(ITYPE,WORK,LEN,LFTNBR,NULLPID)
         JTYPE = 6000+RGTNBR
         CALL CRECV(JTYPE,WORK,LEN2)
      ENDIF
      RETURN
      END
      SUBROUTINE GOP( X, WORK, OP, N)
c
c  Global vector commutative operation using spanning tree.
c
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'INPUT'
      INCLUDE 'CTIMER'
C
      REAL X(N), WORK(N)
      CHARACTER*3 OP
      DIMENSION ORGNL(100)
C
      INTEGER BIT, BYTES, CNT, DIFF, SPSIZE, I, 
     *   PARENT, TROOT, XOR, ROOT
      LOGICAL IFGOT
C
c  All participating processes must have the same process id (PID).
c
c  Input..
c
c    X         the input vector to be used in the operation.  
c    N         the length of the vector.
c    OP        '+'  sum
c              '*'  product
c              'M'  maximum
c              'm'  minimum
c
c  Output..
c
c    X         for all processes, X contains the desired result.
c
c  Workspace
c
c    WORK      used to receive other contributions.  
c
c  Calls:  MYNODE, RECVW, SENDW, XOR
c
c
c     Find temporary root (either the real root, or the lowest
c     numbered node in the active subcube--found by zeroing the
c     CD lowest bits in NID).
c
C
      if (icalld.eq.0) then
        tgop =0.0
        ngop =0
        icalld=1
      endif
      ngop = ngop + 1
      etime1=dclock()
C
c     IF (IFGPRNT) THEN
c        A1=0.0
c        TIME1=SECOND(A1)
c        N100=MIN(9,N)
c        CALL COPY(ORGNL,X,N100)
c     ENDIF
      TYPE  = TYPE+100
      IF (TYPE.GT.99924000) TYPE=TYPE-92000
      TYPER = TYPE-1
      IF (NP.GT.1) THEN
C
C        perform global operation...
C
      ROOT  = NODE0
      TROOT = MAX0((NID/NP)*NP, ROOT)
c
      DIFF = XOR(NID,TROOT)
      IF (DIFF .GE. NP) THEN
         WRITE(6,*) NID,'GOP: CALLED BY NON PARTICIPANT'
         RETURN
      ENDIF
c
c     Accumulate contributions from children, if any
c
      BYTES = WDSIZE*N
      if (ifgprnt.and.nid.eq.302) 
     $   write(6,302) bytes,type,x(1),op,ngop
  302 format(' gop:',2i12,g17.8,1x,a3,i12)
C
      level2=1
    5 continue
      level=level2
      level2=level+level
      IF (mod(nid,level2).ne.0) GO TO 20
C
         CALL CRECV(  TYPE,WORK,BYTES               )
C
C
         CNT = INFOCOUNT()
         IF (CNT .GT. BYTES) 
     $      WRITE(6,*) NID,CNT,'GOP: LONG MESSAGE'
         IF (CNT .LT. BYTES) THEN
            WRITE(6,8) TYPE,BYTES,CNT,level,NID
    8       FORMAT('GOP: SHORT',5I6)
         ENDIF
C
C
         IF (OP.EQ.'+  '.AND.N.GE.20) THEN
            CALL ADD2(X,WORK,N)
         ELSEIF (OP.EQ.'SUM'.AND.N.GE.20) THEN
            CALL ADD2(X,WORK,N)
         ELSEIF (OP.EQ.'MUL'.AND.N.GE.20) THEN
            CALL COL2(X,WORK,N)
         ELSEIF (OP.EQ.'*  '.AND.N.GE.20) THEN
            CALL COL2(X,WORK,N)
         ELSEIF ( (OP.EQ.'m  '.OR.OP.EQ.'M  ') .OR.
     $            (OP.EQ.'MAX'.OR.OP.EQ.'MIN') .OR.
     $            (OP.EQ.'MXA'.OR.OP.EQ.'MNA') .OR.
     $            (OP.EQ.'SUM'.OR.OP.EQ.'MUL') .OR.
     $            (OP.EQ.'+  '.OR.OP.EQ.'*  ')    ) THEN
            DO 10 I = 1, N
               IF (OP .EQ. 'MXA'.and.ABS(WORK(I)).GT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MNA'.and.ABS(WORK(I)).LT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MUL')
     $            X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'SUM')
     $            X(I) = X(I) + WORK(I)
               IF (OP .EQ. '+  ') X(I) = X(I) + WORK(I)
               IF (OP .EQ. '*  ') X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'M  ' .OR. OP .EQ. 'MAX')
     $            X(I) = MAX(X(I),WORK(I))
               IF (OP .EQ. 'm  ' .OR. OP .EQ. 'MIN')
     $            X(I) = MIN(X(I),WORK(I))
   10       CONTINUE
         ELSEIF (OP .EQ. 'ISm') THEN
C           ISAMIN
C           Isolate first occurance of MIN(WORK(I))
            DO 11 I=1,N/2,2
c              IF (WORK(I+1).LT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).LE.X(I+1)) THEN
                  IF (WORK(I+1).LT.X(I+1)) THEN
                     X(I)  =WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
c                 write(6,88) nid,work(i+1),x(i+1),work(i),x(i),i
   88             format(' ISm:',i4,4E12.4,I5)
               ENDIF
   11       CONTINUE
         ELSEIF (OP .EQ. 'ISM') THEN
C           ISAMAX
            DO 12 I=1,N/2,2
c              IF (WORK(I+1).GT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).GE.X(I+1)) THEN
                  IF (WORK(I+1).GT.X(I+1)) THEN
                     X(I)=WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
               ENDIF
   12       CONTINUE
         ENDIF
      IF (LEVEL2.LT.NP) GOTO 5
c
c     Pass result back to parent
c
   20 CONTINUE
C
      IF (nid .NE. 0) THEN
         PARENT = nid-level
         CALL CSEND(TYPE,X,BYTES,PARENT,NULLPID)
      ENDIF
C
C---------------------------------------------------------------
C     AWAIT FINAL ANSWER FROM NODE 0
C---------------------------------------------------------------
C
      IF (PARAM(183).EQ.0) THEN
C
C         We do this via log_2 fan out
C
          LEVEL=NP/2
          IFGOT=.FALSE.
          IF (NID.EQ.ROOT) IFGOT=.TRUE.
C
          DO 1000 I=1,CD
           IF (IFGOT) THEN
              JNID=NID+LEVEL
              CALL CSEND(TYPER,X,BYTES,JNID,NULLPID)
           ELSEIF (MOD(NID,LEVEL).EQ.0) THEN
              CALL CRECV(TYPER,X,BYTES)
              IFGOT=.TRUE.
           ENDIF
           LEVEL=LEVEL/2
 1000     CONTINUE
C
      ELSE
C
C         Use global broadcast
C
          IF (NID.EQ.ROOT) THEN
             CALL CSEND ( TYPE,X,BYTES,ALLNODES,NULLPID)
          ELSE
             CALL CRECV ( TYPE,X,BYTES)
          ENDIF
      ENDIF
C
C
C     End of parallel section....
      ENDIF
C
C     Diagnostics?
C
c     IF (IFGPRNT) THEN
c        TIME2=SECOND(A1)
c        GTIME=TIME2-TIME1
c        ETIME=TIME2-TIME0
c        DO 100 I=1,N100
c           WRITE(6,101) NID,OP,TYPE,I,N,X(I),ORGNL(I),ETIME,GTIME
c 100    CONTINUE
c 101    FORMAT(I3,' GOP ',A3,I6,2I2,4G12.4)
c     ENDIF
      tgop =tgop +(dclock()-etime1)
      if (ifgprnt.and.nid.eq.302) 
     $   write(6,303) bytes,typer,x(1),op,ngop
  303 format(' gop2',2i12,g17.8,1x,a3,i12)
      RETURN
      END
      SUBROUTINE OGOP( X, WORK, OP, N)
c
c  Global vector commutative operation using spanning tree.
c
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'INPUT'
C
      REAL X(N), WORK(N)
      CHARACTER*3 OP
      DIMENSION ORGNL(100)
C
      INTEGER BIT, BYTES, CNT, DIFF, SPSIZE, I, 
     *   PARENT, TROOT, XOR, ROOT
      LOGICAL IFGOT
C
c  All participating processes must have the same process id (PID).
c
c  Input..
c
c    X         the input vector to be used in the operation.  
c    N         the length of the vector.
c    OP        '+'  sum
c              '*'  product
c              'M'  maximum
c              'm'  minimum
c
c  Output..
c
c    X         for all processes, X contains the desired result.
c
c  Workspace
c
c    WORK      used to receive other contributions.  
c
c  Calls:  MYNODE, RECVW, SENDW, XOR
c
c
c     Find temporary root (either the real root, or the lowest
c     numbered node in the active subcube--found by zeroing the
c     CD lowest bits in NID).
c
      IF (IFGPRNT) THEN
         A1=0.0
c         TIME1=SECOND(A1)
         TIME1=dclock()
         N100=MIN(9,N)
         CALL COPY(ORGNL,X,N100)
      ENDIF
      TYPE  = TYPE+100
      IF (TYPE.GT.99924000) TYPE=TYPE-92000
      TYPER = TYPE-1
C
      IF (NP.GT.1) THEN
C
C        perform global operation...
C
      ROOT  = NODE0
      TROOT = MAX0((NID/NP)*NP, ROOT)
c
      DIFF = XOR(NID,TROOT)
      IF (DIFF .GE. NP) THEN
         WRITE(6,*) NID,'GOP: CALLED BY NON PARTICIPANT'
         RETURN
      ENDIF
c
c     Accumulate contributions from children, if any
c
      BIT = NP/2
      BYTES = WDSIZE*N
    5 IF (BIT .LE. DIFF) GO TO 20
         CALL CRECV(   TYPE,WORK,BYTES               )
         CNT = INFOCOUNT()
         IF (CNT .GT. BYTES) 
     $      WRITE(6,*) NID,CNT,'GOP: LONG MESSAGE'
         IF (CNT .LT. BYTES) THEN
            WRITE(6,8) TYPE,BYTES,CNT,BIT,NID
    8       FORMAT('GOP: SHORT',5I6)
         ENDIF
         IF (OP.EQ.'+  '.AND.N.GE.20) THEN
            CALL ADD2(X,WORK,N)
         ELSEIF (OP.EQ.'SUM'.AND.N.GE.20) THEN
            CALL ADD2(X,WORK,N)
         ELSEIF (OP.EQ.'MUL'.AND.N.GE.20) THEN
            CALL COL2(X,WORK,N)
         ELSEIF (OP.EQ.'*  '.AND.N.GE.20) THEN
            CALL COL2(X,WORK,N)
         ELSEIF ( (OP.EQ.'m  '.OR.OP.EQ.'M  ') .OR.
     $            (OP.EQ.'MAX'.OR.OP.EQ.'MIN') .OR.
     $            (OP.EQ.'MXA'.OR.OP.EQ.'MNA') .OR.
     $            (OP.EQ.'SUM'.OR.OP.EQ.'MUL') .OR.
     $            (OP.EQ.'+  '.OR.OP.EQ.'*  ')    ) THEN
            DO 10 I = 1, N
               IF (OP .EQ. 'MXA'.and.ABS(WORK(I)).GT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MNA'.and.ABS(WORK(I)).LT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MUL')
     $            X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'SUM')
     $            X(I) = X(I) + WORK(I)
               IF (OP .EQ. '+  ') X(I) = X(I) + WORK(I)
               IF (OP .EQ. '*  ') X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'M  ' .OR. OP .EQ. 'MAX')
     $            X(I) = MAX(X(I),WORK(I))
               IF (OP .EQ. 'm  ' .OR. OP .EQ. 'MIN')
     $            X(I) = MIN(X(I),WORK(I))
   10       CONTINUE
         ELSEIF (OP .EQ. 'ISm') THEN
C           ISAMIN
C           Isolate first occurance of MIN(WORK(I))
            DO 11 I=1,N/2,2
c              IF (WORK(I+1).LT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).LE.X(I+1)) THEN
                  IF (WORK(I+1).LT.X(I+1)) THEN
                     X(I)  =WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
c                 write(6,88) nid,work(i+1),x(i+1),work(i),x(i),i
   88             format(' ISm:',i4,4E12.4,I5)
               ENDIF
   11       CONTINUE
         ELSEIF (OP .EQ. 'ISM') THEN
C           ISAMAX
            DO 12 I=1,N/2,2
c              IF (WORK(I+1).GT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).GE.X(I+1)) THEN
                  IF (WORK(I+1).GT.X(I+1)) THEN
                     X(I)=WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
               ENDIF
   12       CONTINUE
         ENDIF
         BIT = BIT/2
      GOTO 5
c
c     Pass result back to parent
c
   20 CONTINUE
C
      IF (BIT .NE. 0) THEN
         PARENT = XOR(NID, BIT)
         CALL CSEND(TYPE,X,BYTES,PARENT,NULLPID)
      ELSE
         IF (ROOT.LT.0) CALL CSEND(TYPE,X,BYTES,MANAGER,NULLPID)
      ENDIF
C
C     AWAIT FINAL ANSWER FROM NODE 0
C
      if (param(183).eq.0) then
C
        LEVEL=NP/2
        IFGOT=.FALSE.
        IF (NID.EQ.ROOT) IFGOT=.TRUE.
C
        DO 1000 I=1,CD
         IF (IFGOT) THEN
            JNID=NID+LEVEL
c           write(6,*) nid,'sending to :',jnid,i,level
            CALL CSEND(TYPER,X,BYTES,JNID,NULLPID)
         ELSEIF (MOD(NID,LEVEL).EQ.0) THEN
            CALL CRECV(TYPER,X,BYTES)
            IFGOT=.TRUE.
         ENDIF
         LEVEL=LEVEL/2
 1000   CONTINUE
C
      ELSE
        IF (NID.EQ.ROOT) THEN
           CALL CSEND ( TYPE,X,BYTES,ALLNODES,NULLPID)
        ELSE
           CALL CRECV ( TYPE,X,BYTES)
        ENDIF
      ENDIF
C
C
C     End of parallel section....
      ENDIF
C
C     Diagnostics?
C
      IF (IFGPRNT) THEN
c         TIME2=SECOND(A1)
         TIME2=dclock()
         GTIME=TIME2-TIME1
         ETIME=TIME2-TIME0
         DO 100 I=1,N100
            WRITE(6,101) NID,OP,TYPE,I,N,X(I),ORGNL(I),ETIME,GTIME
  100    CONTINUE
  101    FORMAT(I3,' GOP ',A3,I6,2I2,4G12.4)
      ENDIF
      RETURN
      END
      SUBROUTINE DGOP( X, WORK, OP, N)
c
c  Global vector commutative operation using spanning tree.
c
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INCLUDE 'INPUT'
C
      REAL*8 X(N), WORK(N)
      CHARACTER*3 OP
C
      INTEGER BIT, BYTES, CNT, DIFF, SPSIZE, I, 
     *   PARENT, TROOT, XOR, ROOT
      LOGICAL IFGOT
C
c  All participating processes must have the same process id (PID).
c
c  Input..
c
c    X         the input vector to be used in the operation.  
c    N         the length of the vector.
c    OP        '+'  sum
c              '*'  product
c              'M'  maximum
c              'm'  minimum
c
c  Output..
c
c    X         for all processes, X contains the desired result.
c
c  Workspace
c
c    WORK      used to receive other contributions.  
c
c  Calls:  MYNODE, RECVW, SENDW, XOR
c
c
      TYPE  = TYPE+100
      IF (TYPE.GT.99924000) TYPE=TYPE-92000
      TYPER = TYPE-1
      IF (NP.EQ.1) RETURN
C
C        perform global operation...
C
      ROOT  = NODE0
      TROOT = MAX0((NID/NP)*NP, ROOT)
c
      DIFF = XOR(NID,TROOT)
      IF (DIFF .GE. NP) THEN
         WRITE(6,*) NID,'DGOP: CALLED BY NON PARTICIPANT'
         RETURN
      ENDIF
c
c     Accumulate contributions from children, if any
c
      BYTES = 8*N
C
      level2=1
    5 continue
      level=level2
      level2=level+level
      IF (mod(nid,level2).ne.0) GO TO 20
C
         CALL CRECV(   TYPE,WORK,BYTES               )
C
C
         CNT = INFOCOUNT()
         IF (CNT .GT. BYTES) 
     $      WRITE(6,*) NID,CNT,'DGOP: LONG MESSAGE'
         IF (CNT .LT. BYTES) THEN
            WRITE(6,8) TYPE,BYTES,CNT,level,NID
    8       FORMAT('DGOP: SHORT',5I6)
         ENDIF
C
C
         IF     ( (OP.EQ.'m  '.OR.OP.EQ.'M  ') .OR.
     $            (OP.EQ.'MAX'.OR.OP.EQ.'MIN') .OR.
     $            (OP.EQ.'MXA'.OR.OP.EQ.'MNA') .OR.
     $            (OP.EQ.'SUM'.OR.OP.EQ.'MUL') .OR.
     $            (OP.EQ.'+  '.OR.OP.EQ.'*  ')    ) THEN
            DO 10 I = 1, N
               IF (OP .EQ. 'MXA'.and.ABS(WORK(I)).GT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MNA'.and.ABS(WORK(I)).LT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MUL')
     $            X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'SUM')
     $            X(I) = X(I) + WORK(I)
               IF (OP .EQ. '+  ') X(I) = X(I) + WORK(I)
               IF (OP .EQ. '*  ') X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'M  ' .OR. OP .EQ. 'MAX')
     $            X(I) = MAX(X(I),WORK(I))
               IF (OP .EQ. 'm  ' .OR. OP .EQ. 'MIN')
     $            X(I) = MIN(X(I),WORK(I))
   10       CONTINUE
         ELSEIF (OP .EQ. 'ISm') THEN
C           ISAMIN
C           Isolate first occurance of MIN(WORK(I))
            DO 11 I=1,N/2,2
c              IF (WORK(I+1).LT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).LE.X(I+1)) THEN
                  IF (WORK(I+1).LT.X(I+1)) THEN
                     X(I)  =WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
               ENDIF
   11       CONTINUE
         ELSEIF (OP .EQ. 'ISM') THEN
C           ISAMAX
            DO 12 I=1,N/2,2
c              IF (WORK(I+1).GT.X(I+1)) X(I)=WORK(I)
               IF (WORK(I+1).GE.X(I+1)) THEN
                  IF (WORK(I+1).GT.X(I+1)) THEN
                     X(I)=WORK(I)
                     X(I+1)=WORK(I+1)
                  ELSE
                     X(I)=MIN(X(I),WORK(I))
                  ENDIF
               ENDIF
   12       CONTINUE
         ENDIF
      IF (LEVEL2.LT.NP) GOTO 5
c
c     Pass result back to parent
c
   20 CONTINUE
C
      IF (nid .NE. 0) THEN
         PARENT = nid-level
         CALL CSEND(TYPE,X,BYTES,PARENT,NULLPID)
      ENDIF
C
C---------------------------------------------------------------
C     AWAIT FINAL ANSWER FROM NODE 0
C---------------------------------------------------------------
C
      IF (PARAM(183).EQ.0) THEN
C
C         We do this via log_2 fan out
C
          LEVEL=NP/2
          IFGOT=.FALSE.
          IF (NID.EQ.ROOT) IFGOT=.TRUE.
C
          DO 1000 I=1,CD
           IF (IFGOT) THEN
              JNID=NID+LEVEL
              CALL CSEND(TYPER,X,BYTES,JNID,NULLPID)
           ELSEIF (MOD(NID,LEVEL).EQ.0) THEN
              CALL CRECV(TYPER,X,BYTES)
              IFGOT=.TRUE.
           ENDIF
           LEVEL=LEVEL/2
 1000     CONTINUE
C
      ELSE
C
C         Use global broadcast
C
          IF (NID.EQ.ROOT) THEN
             CALL CSEND ( TYPE,X,BYTES,ALLNODES,NULLPID)
          ELSE
             CALL CRECV ( TYPE,X,BYTES)
          ENDIF
      ENDIF
C
C     End of parallel section....
C
      RETURN
      END
      SUBROUTINE IGOP( X, WORK, OP, N)
c
c  Global vector commutative operation using spanning tree.
c
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      INTEGER ORGNL(100)
C
      INTEGER X(N), WORK(N)
      CHARACTER*3 OP
C
c  All participating processes must have the same process id (PID).
c
c  Input..
c
c    X         the input vector to be used in the operation.  
c    N         the length of the vector.
c    OP        '+'  sum
c              '*'  product
c              'M'  maximum
c              'm'  minimum
c
c  Output..
c
c    X         for all processes, X contains the desired result.
c
c  Workspace
c
c    WORK      used to receive other contributions.  
c
c  Calls:  MYNODE, RECVW, SENDW, XOR
c
      INTEGER BIT, BYTES, CNT, DIFF, SPSIZE, I, 
     *   PARENT, TROOT, XOR, ROOT
c
c     Find temporary root (either the real root, or the lowest
c     numbered node in the active subcube--found by zeroing the
c     CD lowest bits in NID).
c
      IF (IFGPRNT) THEN
         A1=0.0
c         TIME1=SECOND(A1)
         TIME1=dclock()
         N100=MIN(9,N)
         CALL ICOPY(ORGNL,X,N100)
      ENDIF
      TYPE = TYPE+100
      IF (TYPE.GT.99924000) TYPE=TYPE-92000
      TYPER = TYPE-1
C
      IF (NP.GT.1) THEN
C
C        perform global operation...
C
      ROOT  = NODE0
      TROOT = MAX0((NID/NP)*NP, ROOT)
c
      DIFF = XOR(NID,TROOT)
      IF (DIFF .GE. NP) THEN
         WRITE(6,*) NID,'GOP: CALLED BY NON PARTICIPANT'
         RETURN
      ENDIF
c
c     Accumulate contributions from children, if any
c
      BIT = NP/2
      BYTES = ISIZE*N
    5 IF (BIT .LE. DIFF) GO TO 20
         CALL CRECV(   TYPE,WORK,BYTES               )
         CNT = INFOCOUNT()
         IF (CNT .GT. BYTES) 
     $      WRITE(6,*) NID,CNT,'GOP: LONG MESSAGE'
         IF (CNT .LT. BYTES) THEN
            WRITE(6,8) TYPE,BYTES,CNT,BIT,NID,op
    8       FORMAT('IGOP: SHORT',5I6,1x,a3)
         ENDIF
         IF ( (OP.EQ.'m  '.OR.OP.EQ.'M  ') .OR.
     $        (OP.EQ.'MXA'.OR.OP.EQ.'MNA') .OR.
     $        (OP.EQ.'SUM'.OR.OP.EQ.'MUL') .OR.
     $        (OP.EQ.'+  '.OR.OP.EQ.'*  ') ) THEN
            DO 10 I = 1, N
               IF (OP .EQ. 'MXA'.and.ABS(WORK(I)).GT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MNA'.and.ABS(WORK(I)).LT.ABS(X(I)))
     $            X(I) = WORK(I)
               IF (OP .EQ. 'MUL')
     $            X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'SUM')
     $            X(I) = X(I) + WORK(I)
               IF (OP .EQ. '+  ') X(I) = X(I) + WORK(I)
               IF (OP .EQ. '*  ') X(I) = X(I) * WORK(I)
               IF (OP .EQ. 'M  ') X(I) = MAX(X(I),WORK(I))
               IF (OP .EQ. 'm  ') X(I) = MIN(X(I),WORK(I))
   10       CONTINUE
         ELSE
            IF (OP .EQ. 'COM') CALL COMBIN2(X,WORK,N)
         ENDIF
         BIT = BIT/2
      GO TO 5
c
c     Pass result back to parent
c
   20 CONTINUE
C
      IF (BIT .NE. 0) THEN
         PARENT = XOR(NID, BIT)
         CALL CSEND(TYPE,X,BYTES,PARENT,NULLPID)
      ELSE
         IF (ROOT.LT.0) CALL CSEND(TYPE,X,BYTES,MANAGER,NULLPID)
      ENDIF
C
C     AWAIT FINAL ANSWER FROM NODE 0
C
      IF (NID.EQ.ROOT) THEN
         CALL CSEND ( TYPE,X,BYTES,ALLNODES,NULLPID)
      ELSE
         CALL CRECV ( TYPE,X,BYTES)
      ENDIF
C
C     End of global operation
C
      ENDIF
C
C     Diagnostics?
C
      IF (IFGPRNT) THEN
c         TIME2=SECOND(A1)
         TIME2=dclock()
         GTIME=TIME2-TIME1
         ETIME=TIME2-TIME0
         DO 100 I=1,N100
            WRITE(6,101) NID,OP,TYPE,I,N,X(I),ORGNL(I),ETIME,GTIME
  100    CONTINUE
  101    FORMAT(I3,' IGOP ',A3,I6,2I2,2I12,2G12.4)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION XOR(M,N)
c
c  If NOT running on a parallel processor, it is sufficient to 
c  have this routine return a value of XOR=1.
c
c  Pick one of the following:
c
c  UNIX 4.2, f77:
      XOR = OR(M,N)-AND(M,N)
c
c  Intel FTN286:
c     XOR = M.NEQV.N
c
c  Ryan-McFarland Fortran
C      XOR = IEOR(M,N)
c
c     XOR = 0
c     IF(M.EQ.1 .OR.  N.EQ.1) XOR=1
c     IF(M.EQ.0 .AND. N.EQ.0) XOR=0
c     IF(M.EQ.1 .AND. N.EQ.1) XOR=0
c     IF(M.GT.1 .OR.N.GT.1 .OR.M.LT.0.OR.N.LT.0) THEN
c        PRINT*,'ERROR IN XOR'
c        STOP
c     ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
c      FUNCTION SECOND(A)
c      REAL*8 DCLOCK
c      SECOND=DCLOCK(A)
c      SECOND=DCLOCK()
c      RETURN
c      END
      SUBROUTINE LBCAST(IFIF)
C
C  Broadcast logical variable to all processors.
C
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
      LOGICAL IFIF
C
      INTEGER MTYPE
      SAVE    MTYPE
      DATA    MTYPE /400/
C
      IF (NP.EQ.1) RETURN
C
      MTYPE=MTYPE-399
      MTYPE=MOD(MTYPE,6000)+400
      LEN=4
      ITEM=0
      IF (IFIF) ITEM=1
      IF (NID.EQ.NODE0) THEN
         CALL CSEND(MTYPE, ITEM, LEN, ALLNODES ,NULLPID)
      ELSE
         CALL CRECV(MTYPE, ITEM, LEN            )
      ENDIF
      IFIF=.FALSE.
      IF (ITEM.EQ.1) IFIF=.TRUE.
      RETURN
      END
C
c-----------------------------------------------------------------------
      FUNCTION CPUTIME(DUMMY)
      REAL*8 cpu2,dclock
C
C     this function returns the cpu time in seconds
C
      cpu2=dclock()
      CPUTIME = cpu2
      RETURN
      END
C========================================================================
C     Below are the machine specific cube simulation routines
C     which used to be in csim.f
C========================================================================
C
C     When running on the cube, the following MUST be commented out.
C     When running in parallel on the dec's with the Reactive Kernal, 
C          you must use the csimrk.f file in place of the below subroutines.
C
c-----------------------------------------------------------------------
      SUBROUTINE CSEND
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CRECV3(   JTYPE,BUF ,LENR , LENM  )
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CRECV(   JTYPE,BUF ,LEN         )
      INCLUDE 'SIZE'
      COMMON /RECDAT/  LACT
      DIMENSION BUF(1)
      PARAMETER (LCTM0=10)
      COMMON /CREC1/ CDATA(LCTM0)
      DIMENSION CCDATA(LCTM0)
      DIMENSION LCDATA(LCTM0)
      COMMON /CCSIM/ TEXT1
      CHARACTER*10 TEXT,TEXT1
      CHARACTER*1 CCDATA
      LOGICAL     LCDATA
      EQUIVALENCE (CCDATA,CDATA)
      EQUIVALENCE (LCDATA,CDATA)
      CHARACTER*40 FMAT
c
         IIN=10
         READ(IIN,100,END=1500) NDE,TEXT,IC,ITYPE,LACT,NODE,IPID,FMAT
         LACT4=LACT/4
c
         IF (ABS(ITYPE).NE.ABS(JTYPE)) THEN  
            WRITE(6,*) 'Warning: error in receiving data.'
            write(6,*) 'nde,ic,itype,LACT,node,ipid,TEXT',nde,
     $                  ic,itype,LACT,node,ipid,TEXT
            stop
         endif
c
         IF (ITYPE.GT.0) READ(IIN,FMAT,END=2000) (CDATA(I),I=1,LACT4)
         IF (ITYPE.LT.0.and.ic.ne.917) 
     $       READ(IIN,FMAT,END=2000) (CCDATA(I),I=1,LACT)
         IF (ITYPE.LT.0.and.ic.eq.917) 
     $       READ(IIN,FMAT,END=2000) (LCDATA(I),I=1,LACT)
C
      CALL COPY(BUF,CDATA,LACT4)
 1500 CONTINUE
 2000 CONTINUE
  100 FORMAT(///,18X,I4,/,13X,A10,/,12X,I6,/,12X,I6,/,12X,I6,/
     $          ,12X,I6,/,12X,I6,/,14X,A40,/)
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION INFOCOUNT()
      COMMON /RECDAT/  LACT
      INFOCOUNT=LACT
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION IPROBE()
      IPROBE = 0
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION NUMNODES()
      NUMNODES = 1
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION MYHOST()
      MYHOST = 0
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION MYNODE()
      MYNODE = 0
      RETURN
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION MYPID()
      MYPID = 0
      RETURN
      END
c-----------------------------------------------------------------------
      REAL*8 FUNCTION DCLOCK()
      real*4 A(2),etime,s1
      A(1)=0.0
      A(2)=0.0
      S1=ETIME(A)
c     S1=0.
      DCLOCK = S1
c     DCLOCK = 0.0
      RETURN
      END
c-----------------------------------------------------------------------
      function irecv(iegi,a,leni)
      irecv = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine GSSUM(w,nn,w2)
      return
      end
c-----------------------------------------------------------------------
      integer function gihigh(x,n)
      integer x(1),xmax
      xmax = x(1)
      do i=1,n
         xmax = max(xmax,x(i))
      enddo
      gihigh = xmax
      return
      end
c-----------------------------------------------------------------------
      integer function gisum(x,n)
      integer x(1),xsum
      xsum = 0.
      do i=1,n
         xsum = xsum + x(i)
      enddo
      gisum = xsum
      return
      end
c-----------------------------------------------------------------------
      function isend()
      isend=0
      return
      end
c-----------------------------------------------------------------------
      subroutine gsync
      return
      end
c-----------------------------------------------------------------------
      subroutine bcast
      return
      end
c-----------------------------------------------------------------------
      subroutine exitt
      include 'SIZE'
      include 'TOTAL'
      write(6,*) 'quittin'

c     write(6,*) 'quittin2',z,b   ! comment out these 2 lines
c     call exit                   ! to trap the exitt() call
 
      z = -nx1
      z = sqrt(z)
      y = 1./(nx1-lx1)
      y = 0.*y
      a = 1./y
      b = 1./y
      write(6,*) 'quittin3',z,b
 
      call exit

      return
      end
c-----------------------------------------------------------------------
