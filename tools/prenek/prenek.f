C------------------------------------------------------------------------------
C
C                          NEKTON 2.6  2/8/90
C
C			Copyright (C) 1990, by the 
C
C		Massachusetts Institute of Technology  and Nektonics, Inc.
C
C All Rights Reserved
C
C This program is a licenced product of MIT and Nektonics, Inc.,  and it is 
C not to be disclosed to others, copied, distributed, or displayed 
C without prior authorization.
C
C------------------------------------------------------------------------------
C
      SUBROUTINE FPREP
C23456789012345678901234567890123456789012345678901234567890123456789012
C       PreProcessor for Spectral element code.  Inputs geometry and flow
C       data from keyboard or graphics tablet interactively or from instruction
C       file.  - Ed Bullister
C
C   !!?? MININUM NX=3;NY=3; BUT NZ=4!! FOR NEKTON2
      INCLUDE 'basics.inc'
      LOGICAL IFALTR
      PARAMETER(NCBC=NELM*MAXFLD*6)
      CHARACTER FILE*17
      common /cfilold/ filold
      CHARACTER FILOLD*17
      CHARACTER*40 s40
      COMMON/INOUT/  IEXT
      DATA IFALTR /.FALSE./
C     Replaced Data stmts with assignment statements
C
      CALL DATA
      CALL INIT
      CALL GNABLE
      CALL SCROLL(5)
C
1     CONTINUE
      IFCEN=.TRUE.
      nchoic =  2
      ITEM(1)='TYPE IN  NEW  PARAMETERS'
      ITEM(2)='READ PREVIOUS PARAMETERS'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'READ PARAMETER')
      IF(CHOICE.EQ.'TYPE IN  NEW  PARAMETERS')THEN
C        Interactive session
         IN=5
         IFREAD=.FALSE.
         CALL PRS(' $')
         CALL PRS(' Assemble form of equation to be solved $')
         CALL PRS
     $   (' using mouse to toggle logical switches on and off.$')
      ELSE IF(CHOICE.EQ.'READ PREVIOUS PARAMETERS')THEN
C        Read from file
         CALL PRS(' Enter name of previous session$')
         CALL RES(LINE,70)
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         IFREAD=.TRUE.
         FILOLD=LINE
         file=sesion
C        Check How many letters in read name
         lastch=ltrunc(FILOLD,17)
         do 9 i=1,lastch
            int=ichar(FILOLD(i:i))
            if(int.ge.65 .and. int.le.90) int=int+32
            FILOLD(i:i)=char(int)
9        continue
         m=lastch+1
         n=m+3
         FILOLD(m:n) ='.rea'
         IF(FILE.EQ.arg(2)) FILOLD(n+1:n+4)=OLDVER(1:4)
         CALL PRS(' Will Read Parameters From  '//FILOLD//'$')
         CALL OPENF(9,FILOLD,'OLD',3,IERR)
         IF(IERR.NE.0)THEN
            CALL PRS(' Can''t open file '//FILOLD//'   Try again.$')
            IN=5
            GO TO 1
         ENDIF
      ENDIF
C
C     Here is where you jump to to begin setting stuff ??(But what if READing?)
C     Reset stuff back to compiled state
2     CONTINUE
C
      IF(.NOT.IFREAD.OR.IFALTR) THEN
C       Not reading in old session (or ALTERing read parameters)
C       First, set equation type
300     CALL SETEQT
C       Next, Set Required Parameters for this eqtype, nktonv
C       FIX IT IF THEY TRY TO USE INTYPE -1 FOR VERSION 1
        IF(NKTONV.EQ.1 .AND.PARAM(19).EQ.-1)PARAM(19)=1
        CALL SETPAR
      ELSE
C       read in stuff from previous session
C
C       Do this read twice to get parameter values and comments
C
        READ(9,'(a1)',ERR=59)ans
        READ(9,*,ERR=59) VNEKOLD
        VNEKTON=2.61
        NKTONV=VNEKTON
        READ(9,*,ERR=59) NDIM
        IF(NDIM.EQ.2)IF3D=.FALSE.
        IF(NDIM.EQ.3)IF3D=.TRUE.
        READ(9,*,ERR=59) NPARMO
C       Read in Parameters
        DO 1051 IP=1,NPARMO
          READ(9,*,ERR=59) PARAM(IP)
          IF(CPARAM(IP).EQ.'TORDER' .AND.VNEKOLD .LE.2.5)PARAM(IP)=2.0
          IF(CPARAM(IP).EQ.'IOCOMM' .AND.VNEKOLD .LE.2.5)PARAM(IP)=20
          IF(CPARAM(IP).EQ.'COURANT'.AND.VNEKOLD .LE.2.5)PARAM(IP)=0.25
          IF(CPARAM(IP).EQ.'VISCOS '.AND.VNEKOLD .LE.2.5)
     $    PARAM(IP)=PARAM(IP)*PARAM(1)
1051    CONTINUE
C
        rewind(9)
C
        READ(9,'(a1)',ERR=59) cdum
        READ(9,*,ERR=59) dum
        READ(9,*,ERR=59) dum
        READ(9,*,ERR=59) NPARMO
C       Read in Parameters
        DO 1052 IP=1,NPARMO
          call blank(s40,40)
          READ(9,1053,ERR=59) s40
          if (nindx1(cparam(ip),' ',1).eq.0) 
     $       call chcopy(cparam(ip),s40,40)
 1052   CONTINUE
 1053   format(19x,a40)
C
        NPARAM=MAX(NPARAM,NPARMO)
        NPSCAL=PARAM(23)
        READ(9,*,ERR=59)
        READ(9,*,ERR=59) (PCOND (I),I=3,11)
        READ(9,*,ERR=59) (PRHOCP(I),I=3,11)
C       IFFLOW,IFHEAT,IFTRAN,IFNAV
        READ(9,*,ERR=59)NLOGIC
        READ(9,*,ERR=59)IFFLOW
        READ(9,*,ERR=59)IFHEAT
        READ(9,*,ERR=59)IFTRAN
C       IFADVC(1)=IFNAV
        READ(9,*,ERR=59)(IFADVC(I),I=1,NPSCAL+2)
        IF(VNEKOLD.LT.2.6)READ(9,*,ERR=59)(IFTMSH(I),I=1,NPSCAL+2)
        IF(VNEKOLD.GE.2.6)READ(9,*,ERR=59)(IFTMSH(I),I=0,NPSCAL+2)
        IF(NLOGIC.GE.6)READ(9,*,ERR=59)IFAXIS
        IF(NLOGIC.GE.7)READ(9,*,ERR=59)IFSTRS
        IF(NLOGIC.GE.8)READ(9,*,ERR=59)IFSPLIT
        IF(NLOGIC.GE.9)READ(9,*,ERR=59)IFMGRID
        IF(NLOGIC.GE.10)READ(9,*,ERR=59)IFMODEL
        IF(NLOGIC.GE.11)READ(9,*,ERR=59)IFKEPS
        IF(NLOGIC.GE.12)READ(9,*,ERR=59)IFMVBD
        IF(NLOGIC.GE.13)READ(9,*,ERR=59)IFCHAR
C        IF(VNEKOLD.GE.2.4)THEN
C           READ(9,*,ERR=59)NTEXTSW
C           DO 143 I=1,NTEXTSW
C              READ(9,'(2A40)',err=59)TEXTSW(I,1),TEXTSW(I,2)
C143        CONTINUE
C        ENDIF
        NFLDS=1
        IF(IFHEAT)NFLDS=NFLDS+1
        NFLDS=NFLDS+NPSCAL
C       Read Lots more stuff !!??
        CLOSE(UNIT=9)
C       Check that all req'd parameters are set to legal values
C       If something missing, then go to 300
        CALL CKPAR('CHECK',IERR)
      endif
C     Do this interactively independently of whether we read in parameters
      NZ=PARAM(20)
      CALL LEGEND(ZPTS,WGHT,NZ)
310   NCHOIC =  6
      ITEM(1)='ALTER PARAMETERS'
      ITEM(2)='SHOW PARAMETERS'
      ITEM(3)='BUILD INTERACTIVELY'
      ITEM(4)='BUILD FROM FILE'
      ITEM(5)='IMPORT UNIVERSAL FILE'
      ITEM(6)='ABORT (SAVING PARAMETERS)'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'CENTRAL')
      IF(CHOICE.EQ.'SHOW PARAMETERS') THEN
C        Print Out Parameters
         CALL PRS('****** PARAMETERS *****$')
         IF(AXIS.EQ.0)THEN
            WRITE(S,'('' NEKTON VERSION '',G14.1,
     $      '' ;    '',I3,'' DIMENSIONAL '')')VNEKTON,NDIM
            CALL PRS(S//'$')
         ELSE
            CALL PRSRS(' NEKTON VERSION $',VNEKTON,
     $      ' ; AXISYMMETRIC$')
            CALL PRS(S//'$')
         ENDIF
         CALL SCROLL(12)
         CALL CKPAR('SHOW ',IERR)
         GO TO 310
      ELSE IF(CHOICE.EQ.'ALTER PARAMETERS') THEN
         IFALTR=.TRUE.
         GO TO 300
      ENDIF
C
C     Make logicals rational
      IF(.NOT. IFFLOW)THEN
         IFNAV =.FALSE.
         IFSTRS=.FALSE.
         IFSPLIT=.FALSE.
      ENDIF
      IF(.NOT. IFHEAT)THEN
         IFADVC(2)=.FALSE.
      ENDIF
C
C     All parameters set; Ready to build
C     How does build know if we are interactive or reading?
C     On end-of-file on read file, will we always jump to the right place?
C     Do this interactively independently of whether we read in parameters
      IF(CHOICE.EQ.'ABORT (SAVING PARAMETERS)') THEN
         CALL WRTPAR('ABORT     ')
         CLOSE(UNIT=10)
         CALL DEVEX
         CALL EXITT
      ENDIF
      IF(NDIM.EQ.2)THEN
          NSIDES=4
          NEDGES=4
      ELSE IF(NDIM.EQ.3)THEN
          NSIDES=6
          NEDGES=8
      ENDIF
      CALL SCROLL(5)
      IF(CHOICE.EQ.'BUILD INTERACTIVELY') THEN
         IFREAD=.FALSE.
      ELSE IF(CHOICE.EQ.'IMPORT UNIVERSAL FILE') THEN
         IFUNIV = .TRUE.
1070     CALL PRS(
     $   ' Enter full name (including extension) of universal file:$')
         CALL RES(LINE,70)
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         IFREAD=.TRUE.
         FILENM=LINE
         CALL OPENF(8,FILENM,'OLD',3,IERR)
         IF(IERR.NE.0)THEN
            CALL PRS(' Can''t open file '//FILENM//'   Try again.$')
            IN=5
            GOTO 310
         ENDIF
         CALL REAUNI
      ELSE IF(CHOICE.EQ.'BUILD FROM FILE') THEN
C        Read from file
1080     CALL PRS(' Enter name of previous session$')
         CALL PRS('(Default= '//FILOLD//')$')
         CALL RES(LINE,70)
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         IFREAD=.TRUE.
         IF (LINE.NE.' ') THEN
            FILENM=LINE
C           Check How many letters in read name
            lastch=ltrunc(filenm,17)
            do 19 i=1,lastch
               int=ichar(FILENM(i:i))
               if(int.ge.65 .and. int.le.90) int=int+32
               FILENM(I:I)=CHAR(INT)
19           continue
            m=lastch+1
            n=m+3
            FILENM(M:N) ='.rea'
         ELSE
C           Use default file
            FILENM=FILOLD
         ENDIF
         CALL PRS(' Will Read Mesh and B.C. data from  '//FILENM//'$')
         CALL OPENF(9,FILENM,'OLD',3,IERR)
         IF(IERR.NE.0)THEN
            CALL PRS(' Can''t open file $'//FILENM//'   Try again.$')
            IN=5
c           GO TO 1080
            GOTO 310
         ENDIF
      ENDIF
      IF(IFNOSEG)CALL DRCOVR(13)
      CALL BUILD
C     {If interactive} Ask with Menu whether you want to:
C        Specify Outputs
C        Specify History Points
C        Put in Driving Force
C        Exit Prenek [default]
C     Make output specification flag consistent with equation type
      IF(.NOT. IFHEAT) IFTO = .FALSE.
      IF(.NOT. IFHEAT) IFTGO= .FALSE.
      IF(.NOT. IFFLOW) IFVO = .FALSE.
      IF(.NOT. IFFLOW) IFPO = .FALSE.
      IF(NPSCAL.EQ.0 ) IPSCO= 0
C
      NLINE=NLINR+NLINP
      IF (NLINE.GT.0)THEN
         CALL BEEP
         CALL BEEP
         CALL PRS('** CURRENT RESTART PRESOLVE OPTIONS **$')
         CALL PRS('   USE INITIAL COND TO CHANGE/UPDATE.$')
      ENDIF
      DO 105 I=1,NLINE
         CALL PUTS(INITP(I),80)
  105 CONTINUE
      DO 106 I=1,NFLDS
         LINE=INITC(I)
         IF(I.GT.1)LINE=INITC(I+3)
         IF(LINE(1:9).NE.'C Default')THEN
            CALL PRSIS
     $      ('** NOTE: ** CURRENT INITIAL CONDITIONS FOR FIELD$'
     $,     i,'ARE$')
            CALL PRS(LINE//'$')
            CALL PRS(' USE INITIAL COND TO CHANGE/UPDATE.$')
         ENDIF
106   CONTINUE
108   NCHOIC =  3
      ITEM(1)='EXIT'
      ITEM(2)='OUTPUT'
      ITEM(3)='DRIVE FORCE'
      IF(IFTRAN)THEN
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='HISTORY'
      ENDIF
C     The initial conditions are allowed to give a first guess at steady
C     temperature, or, to give a flow field for steady forced convection.
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='INITIAL COND'
      IF(IFTRAN .AND. IFHEAT .AND. IFADVC(2) .AND. (.NOT.IFFLOW))THEN
C        This kludge for velocity field for forced convection calls
C        the initial condition menu
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='FORCING VELOCITY'
      ENDIF
C
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='INTEGRAL QUANTITY'
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='OBJECT'
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='ZOOM'

c     NCHOIC=NCHOIC+1
c     ITEM(NCHOIC)='RSB'

C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'OPTIONS')
      IF(CHOICE.EQ.'HISTORY'.OR.CHOICE.EQ.'OUTPUT')THEN
C        Find out which fields to save
         CALL SETFLD
      ENDIF
      IF(CHOICE.EQ.'EXIT')THEN
C        Exit Prenek
         call prexit
c        call rsb_xxt_set(2,nel)
         call session_exit
      ELSEIF(CHOICE.EQ.'RSB')THEN
         call rsb_xxt_set(1,nel)
      ELSE IF(CHOICE.EQ.'HISTORY')THEN
         CALL HISTRY
      ELSE IF(CHOICE.EQ.'INTEGRAL QUANTITY')THEN
         CALL INTEGQ
      ELSE IF(CHOICE.EQ.'OUTPUT')THEN
         CALL OUTPUT
      ELSE IF(CHOICE.EQ.'OBJECT')THEN
         CALL OBJEC
      ELSE IF(CHOICE.EQ.'DRIVE FORCE')THEN
         CALL DRIVEF
      ELSE IF(CHOICE.EQ.'ZOOM')THEN
         CALL SETZOOM
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         DO 175 IEL=1,NEL
            CALL DRAWEL(IEL)
  175    CONTINUE
      ELSE IF(CHOICE.EQ.'INITIAL COND')THEN
         CALL INTCND
      ELSE IF(CHOICE.EQ.'FORCING VELOCITY')THEN
         CALL PRS
     $   ('Forcing velocity is put in via the initial conditions$')
         CALL PRS
     $   ('menu.  This velocity will remain constant with time.$')
         CALL INTCND
      ELSE
      ENDIF
      GO TO 108
59    CALL PRS(' ERROR READING PARAMETERS FROM FILE.$')
      STOP
60    CALL PRS(' ERROR WRITING; CHECK DISK QUOTA$')
      STOP
      END
C
      SUBROUTINE REAUNI
      INCLUDE 'basics.inc'
      LOGICAL IFSTRT,IFEND
      CHARACTER*80 LINE80
      INTEGER IICORN(8)
C
C     First read in table of nodes
      READ(8,'(A80)',ERR=13,END=13)
      DO 1 I=1,100000
         READ(8,'(A80)',ERR=13,END=13)LINE80
         IF(LINE80(1:6) .EQ. '    -1')IFSTRT = .TRUE.
         IF(LINE80(1:6) .NE. '    -1')IFSTRT = .FALSE.
         IF(IFSTRT)THEN
            READ(8,*,ERR=13,END=13)
            READ(8,*,ERR=13,END=13)IDSET
            IF(IDSET .EQ. 15) THEN
C              NODES
               DO 2 II = 1,NELM*NXM*NYM
                  READ(8,'(A80)',ERR=13,END=13)LINE80
                  IF(LINE80(1:6) .EQ. '    -1')IFEND = .TRUE.
                  IF(LINE80(1:6) .NE. '    -1')IFEND = .FALSE.
                  IF(IFEND) BACKSPACE(8)
                  IF(IFEND) BACKSPACE(8)
                  IF(IFEND) GO TO 3
C                 USE XPTS as scratch array for nodes
                  READ(LINE80,'(4I10,1P3E13.5)',ERR=14,END=14)
     $            IDUM1,IDUM2,IDUM3,IDUM4,
     $            XPTS(II,1,1),YPTS(II,1,1),ZPTN(II,1,1)
2              CONTINUE
               CALL PRIS
     $         (II,' Nodes is too many for this configuration$')
               CALL PRS(' Contact your NEKTON Support Rep.$')
3              CONTINUE
               NNODES = II-1
               CALL PRSIS('Read in $',NNODES,' Nodes$')
            ENDIF
            IF(IDSET .EQ. 71) THEN
C              ELEMENTS
               DO 42 IEL = 1,NELM
                  READ(8,'(A80)',ERR=13,END=13)LINE80
                  IF(LINE80(1:6) .EQ. '    -1')IFEND = .TRUE.
                  IF(LINE80(1:6) .NE. '    -1')IFEND = .FALSE.
                  IF(IFEND) GO TO 43
C                 Fill Up element arrays
                  READ(LINE80,'(7I10)',ERR=14,END=14)
     $            IDUM1,IDUM2,IDUM3,IDUM4,IDUM5,IDUM6,NCORNS
                  IF(IF3D .AND. NCORNS .NE. 8)THEN
                     CALL PRSIS
     $               ('ERROR: ELEMENT HAS $',NCORNS,' CORNERS$')
                     CALL PRS('3-D ELEMENTS MUST HAVE 8 CORNERS$')
                     STOP
                  ELSE IF((.NOT. IF3D) .AND. NCORNS .NE. 4)THEN
                     CALL PRSIS
     $               ('ERROR: ELEMENT HAS $',NCORNS,' CORNERS$')
                     CALL PRS('2-D ELEMENTS MUST HAVE 4 CORNERS$')
                     STOP
                  ENDIF
                  READ(8,'(A80)',ERR=13,END=13)LINE80
                  READ(LINE80,'(8I10)',ERR=14,END=14)
     $            (IICORN(IC),IC=1,NCORNS)
C
                  DO 50 IC=1,NCORNS
                     X(IEL,IC) = XPTS(IICORN(IC),1,1)
                     Y(IEL,IC) = YPTS(IICORN(IC),1,1)
                     Z(IEL,IC) = ZPTN(IICORN(IC),1,1)
50                CONTINUE
                  NUMAPT(IEL) = 1
                  ILETAP = ILETAP + 1
                  IF(ILETAP.LE.122)LETAPT(IEL) = CHAR(ILETAP)
                  IF(ILETAP.GT.122)LETAPT(IEL) = ' '
42             CONTINUE
               CALL PRSIS
     $         ('$',IEL,' Elements is too many for this configuration$')
               CALL PRS(' Contact your NEKTON Support Rep.$')
43             CONTINUE
               NEL = IEL-1
               CALL PRSIS('Read in $',NEL,' Elements$')
               GO TO 15
            ENDIF
         ENDIF
1     CONTINUE
13    CONTINUE
      CALL PRS('Error Reading Data$')
C     FILL UP Stuff READAT would normally have done
15    NDIM = 2
      IF(IF3D)NDIM=3
      XFAC = ABS( RMAXX(X,NELM,NEL,NCORNS)-RMINX(X,NELM,NEL,NCORNS) )*2.
      YFAC = ABS( RMAXX(Y,NELM,NEL,NCORNS)-RMINX(Y,NELM,NEL,NCORNS) )*2.
      IF(XFAC/YFAC .LT. 1.5 .OR. XFAC/YFAC .GT. 0.66) THEN
C        MAKE THEM THE SAME
         YFAC = XFAC
      ENDIF
      XZERO = RBARX(X,NELM,NEL,NCORNS) - XFAC/2.
      YZERO = RBARX(Y,NELM,NEL,NCORNS) - YFAC/2.
      RETURN
14    CONTINUE
      CALL PRS('Error Reading Internal Data$')
      STOP
      END
      FUNCTION RMAXX(R,NELM,NEL,NC)
      REAL R(NELM,8)
      RM=R(1,1)
      DO 1 I=1,NEL
      DO 1 J=1,NC
         IF(R(I,J) .GT. RM) RM=R(I,J)
1     CONTINUE
      RMAXX = RM
      RETURN
      END
      FUNCTION RMINX(R,NELM,NEL,NC)
      REAL R(NELM,8)
      RM=R(1,1)
      DO 1 I=1,NEL
      DO 1 J=1,NC
         IF(R(I,J) .LT. RM) RM=R(I,J)
1     CONTINUE
      RMINX = RM
      RETURN
      END
      FUNCTION RBARX(R,NELM,NEL,NC)
      REAL R(NELM,8)
      RB =0.0
      DO 1 I=1,NEL
      DO 1 J=1,NC
         RB = RB + R(I,J)
1     CONTINUE
      RBARX = RB / (NEL*NC)
      RETURN
      END
      FUNCTION RMAX(R,N)
      REAL R(N)
      RM=R(1)
      DO 1 I=1,N
         IF(R(I) .GT. RM) RM=R(I)
1     CONTINUE
      RMAX = RM
      RETURN
      END
      FUNCTION RMIN(R,N)
      REAL R(N)
      RM=R(1)
      DO 1 I=1,N
         IF(R(I) .LT. RM) RM=R(I)
1     CONTINUE
      RMIN = RM
      RETURN
      END
      FUNCTION RBAR(R,N)
      REAL R(N)
      RB =0.0
      DO 1 I=1,N
         RB = RB + R(I)
1     CONTINUE
      RBAR = RB / N
      RETURN
      END
      SUBROUTINE WRTPAR(CFLAG)
      INCLUDE 'basics.inc'
      PARAMETER(NCBC=NELM*MAXFLD*6)
      CHARACTER FILE*10,CFLAG*10
      CHARACTER*1 s401(40)
      COMMON/INOUT/  IEXT
C
C     Write out parameter stuff
C     First check B.C.'s to set logical switches
      IFMVBD    = .FALSE.
      IFTMSH(0) = .FALSE.
      IF(CFLAG.EQ.'FULL DUMP ')THEN
         IF(.NOT.IF3D) NSIDES=4
         IF(     IF3D) NSIDES=6
         DO 1 IEL=1,NEL
            DO 1 ISIDE=1,NSIDES
               CBC3 = CBC(ISIDE,IEL,1)
               IF( CBC3(1:1).EQ.'M' .OR. CBC3(1:1).EQ.'m')IFMVBD=.TRUE.
C              Test for moving mesh in solid
               CBC3 = CBC(ISIDE,IEL,2)
               IF(CBC3(1:1).EQ.'M'.OR.CBC3(1:1).EQ.'m')IFTMSH(0)=.TRUE.
1        CONTINUE
      ENDIF
C
      M=IEXT
      n=m+3
      filenm = sesion
C
      FILENM(M:N) ='.rea'
      CALL OPENF(10,FILENM,'NEW',1,IERR)
      CALL PRS('Writing Parameters to file$')
      write(10,*,err=60)'****** PARAMETERS *****'
      WRITE(10,*,err=60)VNEKTON,' NEKTON VERSION '
      WRITE(10,*,err=60)NDIM,   ' DIMENSIONAL RUN'
      WRITE(10,*,err=60)NPARAM, ' PARAMETERS FOLLOW'
C     Print Out Parameters
      DO 1050 IP=1,NPARAM
         call chcopy(s401,cparam(ip),40)
         l=ltrunc(s401,40)
         WRITE(10,'(G14.6,5X,40A1)',ERR=60) PARAM(IP),(s401(j),j=1,l)
1050  CONTINUE
      WRITE(10,*)'     4  Lines of passive scalar data follows',
     $'2 CONDUCT; 2RHOCP'
      write(10,'(5g14.6)',ERR=60)(PCOND (I),I=3,11)
      write(10,'(5g14.6)',ERR=60)(PRHOCP(I),I=3,11)
      NLOGIC=13
      WRITE(10,*)NLOGIC,'  LOGICAL SWITCHES FOLLOW'
      WRITE(10,*)IFFLOW,'     IFFLOW'
      WRITE(10,*)IFHEAT,'     IFHEAT'
      WRITE(10,*)IFTRAN,'     IFTRAN'
C      WRITE(10,*)IFNAV ,'     IFNAV '
      WRITE(10,*)IFADVC,
     $' IFNAV & IFADVC (convection in P.S. fields)'
      WRITE(10,*)IFTMSH,
     $' IFTMSH (IF mesh for this field is T mesh)'
      WRITE(10,*)IFAXIS,'     IFAXIS'
      WRITE(10,*)IFSTRS,'     IFSTRS'
      WRITE(10,*)IFSPLIT,'     IFSPLIT'
      WRITE(10,*)IFMGRID,'     IFMGRID'
      WRITE(10,*)IFMODEL,'     IFMODEL'
      WRITE(10,*)IFKEPS ,'     IFKEPS'
      WRITE(10,*)IFMVBD ,'     IFMVBD'
      WRITE(10,*)IFCHAR ,'     IFCHAR'
C
      RETURN
60    CALL PRS(' ERROR WRITING; CHECK DISK QUOTA$')
      STOP
      END
      SUBROUTINE SETEQT
      INCLUDE 'basics.inc'
      CHARACTER*26 ITEMON
      CHARACTER*120 LSET
C     Sets equation type based on logical parameters
C     Toggle parameters on and off
C     Check if combinations of switches is legal before exiting menu
C
1     CONTINUE
C     Write current settings
      LSET=' '
      ISET=1
      IF(IFSTRS)THEN
         LSET(ISET:ISET+9)=' STRESSFUL'
         ISET=ISET+10
      ENDIF
      IF(IFTRAN)THEN
         LSET(ISET:ISET+8)=' UNSTEADY'
         ISET=ISET+9
      ELSE
         LSET(ISET:ISET+6)=' STEADY'
         ISET=ISET+7
      ENDIF
      IF(IFAXIS)THEN
         LSET(ISET:ISET+12)=' AXISYMMETRIC'
         ISET=ISET+13
      ELSE IF(IF3D)THEN
         LSET(ISET:ISET+17)=' THREE DIMENSIONAL'
         ISET=ISET+18
      ELSE
         LSET(ISET:ISET+15)=' TWO DIMENSIONAL'
         ISET=ISET+16
      ENDIF
      IF(IFFLOW)THEN
         IF(IFNAV)THEN
            LSET(ISET:ISET+13)=' NAVIER-STOKES'
            ISET=ISET+14
         ELSE
            LSET(ISET:ISET+6)=' STOKES'
            ISET=ISET+7
         ENDIF
      ENDIF
      IF(IFHEAT)THEN
         IF(IFFLOW)THEN
            LSET(ISET:ISET+17)=' AND HEAT TRANSFER'
            ISET=ISET+18
         ELSE
            LSET(ISET:ISET+13)=' HEAT TRANSFER'
            ISET=ISET+14
         ENDIF
      ENDIF
      IF(NPSCAL.GT.0)THEN
         LSET(ISET:ISET+22)=' WITH   PASSIVE SCALARS'
         WRITE(LSET(ISET+6:ISET+6),'(I1)')NPSCAL
         ISET=ISET+23
      ENDIF
      CALL PRS(' CURRENT SETTING:'//LSET(1:60)//'$')
      IF(ISET.GT.60)CALL PRS(LSET(61:120)//'$')
      ITEM(1)='ACCEPT CURRENT SWITCHES'
C      ITEM(2)='NEKTON VERSION'
      ITEM(2)='              '
      ITEM(3)='DIMENSION     '
      IF (.NOT.IFTRAN) ITEM(4)='STEADY      '
      IF (     IFTRAN) ITEM(4)='UNSTEADY      '
      ITEM(5)='FLUID FLOW    '
      IF(IFFLOW)THEN
         ITEM(6)='ADVECTION         '
         ITEM(7)='STRESS FORMULATION'
         ITEM(8)='SPLIT  FORMULATION'
         ITEM(9)='TURBULENCE MODEL  '
      ELSE
         ITEM(6)='  '
         ITEM(7)='  '
         ITEM(8)=' '
         ITEM(9)=' '
      ENDIF
      ITEM(10)='HEAT TRANSFER '
      ITEM(11)='ADDL PASSIVE SCALARS'
      ITEM(12)='MULTIGRID SOLUTION'
      NCHOIC=12
C
C      ITEMON=ITEM(2)
C      IF(NKTONV.EQ.2)ITEMON(25:25)='2'
C      IF(NKTONV.EQ.1)ITEMON(25:25)='1'
C      ITEM(2)=ITEMON
      DO 20 I=3,12
         ITEMON=ITEM(I)
         IF(I.EQ.3)THEN
            IF(     IF3D)ITEMON(25:25)='3'
            IF(.NOT.IF3D)ITEMON(25:25)='2'
            IF(   IFAXIS)ITEMON(25:25)='A'
c        ELSE IF(I.EQ.4)THEN
c           IF(     IFTRAN)ITEMON(25:25)='Y'
c           IF(.NOT.IFTRAN)ITEMON(25:25)='N'
         ELSE IF(I.EQ.5)THEN
            IF(     IFFLOW)ITEMON(25:25)='Y'
            IF(.NOT.IFFLOW)ITEMON(25:25)='N'
         ELSE IF(I.EQ.6 .AND. IFFLOW)THEN
            IF(     IFNAV )ITEMON(25:25)='Y'
            IF(.NOT.IFNAV )ITEMON(25:25)='N'
         ELSE IF(I.EQ.7 .AND. IFFLOW)THEN
            IF(     IFSTRS)ITEMON(25:25)='Y'
            IF(.NOT.IFSTRS)ITEMON(25:25)='N'
         ELSE IF(I.EQ.8 .AND. IFFLOW)THEN
            IF(     IFSPLIT)ITEMON(25:25)='Y'
            IF(.NOT.IFSPLIT)ITEMON(25:25)='N'
         ELSE IF(I.EQ.9 .AND. IFFLOW)THEN
            ITEMON(24:25)=TURBMOD(1:2)
            IF(TURBMOD(1:2).EQ.'NO')ITEMON(24:25)=' N'
         ELSE IF(I.EQ.10)THEN
            IF(     IFHEAT)ITEMON(25:25)='Y'
            IF(.NOT.IFHEAT)ITEMON(25:25)='N'
         ELSE IF(I.EQ.11)THEN
            WRITE(ITEMON(25:25),'(I1)')NPSCAL
C           MPSCAL=9 IN PARAMETER STMT
         ELSE IF(I.EQ.12)THEN
            IF(     IFMGRID)ITEMON(25:25)='Y'
            IF(.NOT.IFMGRID)ITEMON(25:25)='N'
         ENDIF
         ITEM(I)=ITEMON
20    CONTINUE
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOECHO')
C     Toggle logical flags based on choices
C     Do the opposite of what flag currently is, so it gets changed.
C     Round robin for dimensions
      IF(CHOICE.EQ.'DIMENSION               A')IF3D  =.TRUE.
      IF(CHOICE.EQ.'DIMENSION               A')IFAXIS=.FALSE.
      IF(CHOICE.EQ.'DIMENSION               2')IF3D  =.FALSE.
      IF(CHOICE.EQ.'DIMENSION               2')IFAXIS=.TRUE.
      IF(CHOICE.EQ.'DIMENSION               3')IF3D  =.FALSE.
      IF(CHOICE.EQ.'DIMENSION               3')IFAXIS=.FALSE.
      IF(CHOICE.EQ.'UNSTEADY                 ')IFTRAN=.FALSE.
      IF(CHOICE.EQ.'STEADY                   ')IFTRAN=.TRUE.
      IF(CHOICE.EQ.'HEAT TRANSFER           Y')IFHEAT=.FALSE.
      IF(CHOICE.EQ.'HEAT TRANSFER           N')IFHEAT=.TRUE.
      IF(CHOICE.EQ.'FLUID FLOW              Y')IFFLOW=.FALSE.
      IF(CHOICE.EQ.'FLUID FLOW              N')IFFLOW=.TRUE.
      IF(CHOICE.EQ.'ADVECTION               Y')IFNAV =.FALSE.
      IF(CHOICE.EQ.'ADVECTION               N')IFNAV =.TRUE.
      IF(CHOICE.EQ.'STRESS FORMULATION      Y')IFSTRS=.FALSE.
      IF(CHOICE.EQ.'STRESS FORMULATION      N')THEN
         IFSTRS = .TRUE.
         IFSPLIT= .FALSE.
      ENDIF
      IF(CHOICE.EQ.'SPLIT  FORMULATION      Y')THEN
         IFSPLIT=.FALSE.
         CALL PRS('Time stepping scheme set to Second Order$')
         PARAM(27) = 2
         REQD (27) = 'P'
      ENDIF
      IF(CHOICE.EQ.'SPLIT  FORMULATION      N')THEN
         IFSPLIT=.TRUE.
         IFSTRS =.FALSE.
C        Set TORDER
         CALL PRS('Time stepping scheme set to First Order$')
         PARAM(27) = 1
         REQD (27) = ' '
      ENDIF
      IF(CHOICE(1:10).EQ.'TURBULENCE')THEN
C         CALL PRS
C     $   ('**** WARNING **** Turbulence Models under development.$')
C         CALL PRS
C     $   ('Use at your own risk. By hitting a carriage return,$')
C         CALL PRS('you acknowledge having read this warning.$')
C         CALL RES(LINE,0)
          CALL PRS('Turbulence Models not yet implemented$')
      ENDIF
      IF(.FALSE.)THEN
       IF(CHOICE.EQ.'TURBULENCE MODEL        N')TURBMOD='KEPSRNG'
       IF(CHOICE.EQ.'TURBULENCE MODEL       KE')TURBMOD='MIXLRNG'
       IF(CHOICE.EQ.'TURBULENCE MODEL       MI')TURBMOD='NONE'
C
       IF(TURBMOD.EQ.'KEPSRNG'.OR.TURBMOD.EQ.'MIXLRNG') IFMODEL=.TRUE.
       IF(TURBMOD.EQ.'KEPSRNG'                        ) IFKEPS =.TRUE.
      ENDIF
C
C      IF(CHOICE.EQ.'MULTIGRID SOLUTION      Y')IFMGRID=.FALSE.
C      IF(CHOICE.EQ.'MULTIGRID SOLUTION      N')IFMGRID=.TRUE.
      IF(CHOICE.EQ.'MULTIGRID SOLUTION      N')THEN
         CALL PRS('Multigrid not implemented in this version$')
      ENDIF
      IF(CHOICE(1:24).EQ.'ADDL PASSIVE SCALARS    ')NPSCAL=NPSCAL+1
      IF(NPSCAL.GT.MPSCAL)NPSCAL=0
      PARAM(23)=NPSCAL
C
      IF(CHOICE.EQ.'HEAT TRANSFER           Y' .OR.
     $   CHOICE.EQ.'HEAT TRANSFER           N' .OR.
     $   CHOICE.EQ.'FLUID FLOW              Y' .OR.
     $   CHOICE.EQ.'FLUID FLOW              N')THEN
C        Reset the IFADVC switches to the default value for the new eq type
         DO 402 I=1,11
            IF(.NOT.IFFLOW)IFADVC(I)=.FALSE.
            IF(     IFFLOW.AND.I.NE.1)IFADVC(I)=.TRUE.
402      CONTINUE
      ENDIF
C
      IF(CHOICE(1:23).EQ.'ACCEPT CURRENT SWITCHES')THEN
C         Set appropriate stuff (eqtype=; check for reasonableness)
          IF(     IFAXIS)AXIS=1
          IF(.NOT.IFAXIS)AXIS=0
          IF(     IF3D  )NDIM=3
          IF(.NOT.IF3D  )NDIM=2
          NFLDS=1
          IF(     IFHEAT)NFLDS=NFLDS+1
          NFLDS=NFLDS+NPSCAL
          IF(NKTONV.EQ.2)THEN
C            Nekton 2.0
C            Check if requested eqtype is available as of yet
             IF(IFNAV.and..NOT.IFTRAN)THEN
                 CALL PRS(' ****** ERROR ********$')
                 CALL PRS(' Navier-Stokes solver calculates for $')
                 CALL PRS
     $           (' steady flows by solving the time dependent$')
                 CALL PRS(' equations until the flow reaches steady$')
                 CALL PRS
     $           (' state.  You must turn on the unsteady option$')
                 CALL PRS(' if the advection term is turned on.$')
                 CALL BEEP
                 GO TO 1
             ENDIF
             IF(IFMGRID .AND. .NOT. IFFLOW)THEN
                 CALL PRS(' ****** ERROR ********$')
                 CALL PRS(' Multigrid only works with fluid flow$')
                 CALL PRS(' in which Nu is constant.$')
                 CALL BEEP
                 GO TO 1
             ENDIF
             IF(IFFLOW.AND.IFSPLIT.and..NOT.IFTRAN)THEN
                 CALL PRS(' ****** ERROR ********$')
                 CALL PRS(' Split Formulation only works for unsteady$')
                 CALL PRS(' Flows$')
                 CALL BEEP
                 GO TO 1
             ENDIF
           ENDIF
C          Check for errors
           CALL CHKEQT(IERR)
           IF(IERR.NE.0) GO TO 1
C          Now decide whether there is convection in temperature and passive
C          scalar fields
C          No steady-state forced convection yet.
           IFNEEDC=.FALSE.
           IF(IFFLOW .AND. IFADVC(1))IFNEEDC=.TRUE.
           IF(IFHEAT.AND.IFTRAN)THEN
C             IFADVC(2) is for heat
              DO 50 I=0,NPSCAL
                 IFLD=I+2
49               IF(IFLD.EQ.2)THEN
                    CALL PRS('Convection in Temperature Field$')
                    IF(.NOT.IFFLOW)THEN
                     CALL PRS('Forced convection is turned on here if$')
                     CALL PRS('you want a forcing term from a steady$')
                     CALL PRS('velocity field.  You will define that$')
                     CALL PRS('field later as a fortran function or$')
                     CALL PRS('to be read in from a field file.$')
                    ENDIF
                 ENDIF
                 IF(IFLD.GT.2)CALL PRSIS(' PASSIVE SCALAR$',I,'$')
                 NCHOIC=2
C                Passive scalars currently Must use temperature mesh
C                 NCHOIC=3
                 ITEM(1)='ACCEPT CURRENT SWITCHES'
                 IF(     IFADVC(IFLD))ITEM(2)=
     $           'WITH CONVECTION         Y'
                 IF(.NOT.IFADVC(IFLD))ITEM(2)=
     $           'WITH CONVECTION         N'
                 IF(.NOT.IFTMSH(IFLD)) ITEM(3)='FLUID       MESH'
                 IF(     IFTMSH(IFLD)) ITEM(3)='TEMPERATURE MESH'
                 CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOECHO')
                 IF(CHOICE.EQ.'WITH CONVECTION         Y')THEN
                    IFADVC(IFLD)=.FALSE.
                    GO TO 49
                 ELSE IF(CHOICE.EQ.'WITH CONVECTION         N')THEN
                    IFADVC(IFLD)=.TRUE.
                    GO TO 49
                 ELSE IF(CHOICE.EQ.'FLUID       MESH')THEN
                    IFTMSH(IFLD)=.TRUE.
                    GO TO 49
                 ELSE IF(CHOICE.EQ.'TEMPERATURE MESH')THEN
                    IFTMSH(IFLD)=.FALSE.
                    GO TO 49
                 ELSE IF(CHOICE.EQ.'ACCEPT CURRENT SWITCHES')THEN
                    CONTINUE
                    IF(IFADVC(IFLD))IFNEEDC=.TRUE.
                 ENDIF
50            CONTINUE
           ENDIF
C          Characteristics
           IF(IFNEEDC)THEN
149           CONTINUE
              ITEM(1)='ACCEPT CURRENT SWITCHES'
              IF(     IFCHAR)ITEM(2)='CHARACTERISTICS         Y'
              IF(.NOT.IFCHAR)ITEM(2)='CHARACTERISTICS         N'
              NCHOIC=2
              CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOECHO')
              IF(CHOICE.EQ.'CHARACTERISTICS         Y')THEN
                 IFCHAR=.FALSE.
                 PARAM(26) = 0.25
                 CALL PRS
     $           ('*** Default COURANT Number reset to 0.25 ****$')
                 CALL PRS('Modify in Set Parameter Menu if desired$')
                 GO TO 149
              ELSE IF(CHOICE.EQ.'CHARACTERISTICS         N')THEN
                 IFCHAR=.TRUE.
                 PARAM(26) = 1.0
                 CALL PRS
     $           ('*** Default COURANT Number set to 1.0 ****$')
                 CALL PRS('Modify in Set Parameter Menu if desired$')
                 GO TO 149
              ELSE IF(CHOICE.EQ.'ACCEPT CURRENT SWITCHES')THEN
                 CONTINUE
              ENDIF
           ELSE
              IFCHAR=.FALSE.
           ENDIF
C
           IF(IFNOSEG) CALL DRCOVR(13)
           RETURN
      ENDIF
      GO TO 1
      END
      SUBROUTINE CHKPAR(IERR)
      include 'basics.inc'
      IERR=0
      IF(IFMGRID)THEN
         CALL PRS('Automatically setting multigrid parameters$')
         CALL GENNMG(PARAM)
      ENDIF
      RETURN
      END
      SUBROUTINE CHKEQT(IERR)
      include 'basics.inc'
C
C     Checks for illegal equation settings.  IERR=1 if any were found.
      IERR=0
      IF(NPSCAL.GT.0 .AND. .NOT. IFHEAT) THEN
         CALL PRS('**ERROR** Cannot have extra passive scalars$')
         CALL PRS('without heat transfer.  Try again.$')
         CALL BEEP
         IERR=1
      ENDIF
      IF(.NOT. IFFLOW .AND. .NOT. IFHEAT) THEN
         CALL PRS('**ERROR** Not specifying any variables makes for a'//
     $   ' poor simulation.$')
         CALL PRS('Specify heat and/or fluid.$')
         CALL BEEP
         IERR=1
      ENDIF
      RETURN
      END
      SUBROUTINE SETFLD
      include 'basics.inc'
C     Default OUTPUTs
C     Sets up which quantities to be dumped
C     Boxes with XX, etc inside change color when specified.
      CALL PRS(' Use mouse to toggle Output fields off and on$')
C     Draw Boxes
CC      DO 1 I=1,5
C         YB=0.05
C         YT=0.15
C         XL=1.0+(I-1)/20.0
C         XR=1.0+I    /20.0
C         CALL BEGINB(XPHY(XL),YPHY(YB))
C         CALL DRAWSC(XR,YB)
C         CALL DRAWSC(XR,YT)
C         CALL DRAWSC(XL,YT)
C         CALL ENDP
CC         CALL GWRITE()
C1     CONTINUE
C      IF(CHOICE.EQ.'OUTPUT')THEN
CC        Put in output array
C      ELSE
CC        Put in History array
C      ENDIF
      RETURN
      END
      SUBROUTINE SETREQ
      INCLUDE 'basics.inc'
C
C     REQD(IPARAM)  =I Irrelevant
C                           =O Optional
C                           =R Required NonZero
C                           =P Required NonZero and positive
      DO 10 I=1,40
         REQD(I)='I'
10    CONTINUE
C
C     DENSITY
      IF(IFFLOW)              REQD(1)='O'
C     Absolute VISCOS
      IF(IFFLOW)              REQD(2)='R'
C     BETAG
      IF(IFFLOW .AND. IFHEAT) REQD(3)='O'
C     GTHETA
      IF(IFFLOW .AND. IFHEAT) REQD(4)='O'
C     PGRADX
      IF(IFFLOW .AND. NKTONV.EQ.1) REQD(5)='O'
C     FLOWRATE
      IF(IFFLOW .AND. NKTONV.EQ.1) REQD(6)='O'
C     RHOCP
      IF(IFHEAT .AND. IFTRAN) REQD(7)='P'
C     CONDUCT
      IF(IFHEAT             ) REQD(8)='P'
C     QVOL
      IF(IFHEAT) REQD(9)='O'
C     FINTIME
      IF(IFTRAN .AND. NKTONV.GT.1)REQD(10)='O'
C     NSTEPS
      IF(IFTRAN             )REQD(11)='O'
C     DT
      IF(IFTRAN             )REQD(12)='O'
C     PCOURANT
      IF(IFTRAN .and. NKTONV.EQ.1)REQD(13)='O'
C     IOTIME
      IF(IFTRAN.AND.NKTONV.GT.1)REQD(14)='O'
C     IOSTEP
      IF(IFTRAN             )REQD(15)='O'
C     EQTYPE
                             REQD(16)='I'
C     AXIS
                             REQD(17)='I'
C     GRID
                             REQD(18)='O'
C     INTYPE
C      IF(IFTRAN)             REQD(19)='O'
C     NORDER
                             REQD(20)='P'
C     DIVERGENCE
      IF(NKTONV.NE.1 .AND. IFFLOW .AND. .NOT.IFTRAN)REQD(21)='O'
C     HELMHOLTZ
      IF(NKTONV.NE.1 .AND. .NOT.(IFFLOW .AND. IFTRAN))REQD(22)='O'
      IF(NKTONV.NE.1 .AND. IFHEAT                    )REQD(22)='O'
C     TOLREL
      IF(NKTONV.NE.1 .AND.IFFLOW .AND. IFTRAN) REQD(24)='O'
C     TOLABS
      IF(NKTONV.NE.1 .AND.IFFLOW .AND. IFTRAN) REQD(25)='O'
C     COURANT
      IF(IFNEEDC) REQD(26)='O'
C     TORDER
      IF(IFTRAN             )REQD(27)='P'
C
C     XMAGNET
C      IF(IFFLOW) REQD(36)='O'
C     MultiGrid
      IF(IFMGRID)REQD(37)='M'
      IF(IFMGRID)REQD(38)='M'
      IF(IFMGRID)REQD(39)='M'
      IF(IFMGRID)REQD(40)='M'
      IF(IFMGRID)REQD(41)='M'
      RETURN
      END
      SUBROUTINE SETPAR
C     Based on EQTYPE and  and IF3D, Prompt for Parameters to be set
      include 'basics.inc'
      CHARACTER*19 BOXLAB
      CALL GINDIS
C     4107 Tends to send <cr> when you disable gin.  This reads it.
      IF(IF4107)CALL RES(LINE,70)
      IF(IF4107)THEN
      IF(IFLEARN)WRITE(3,'(A70)') LINE
      CALL PRS
     $('SET PARAMETER menu.  TYPE in with KEYBOARD (drop mouse)$')
      CALL PRS
     $('Type in new value or <cr> to keep current (default) value$')
      ELSE IF(IFSUN) THEN
         CALL PRS('SET PARAMETER menu.  KEY in new value with KEYPAD$')
         CALL PRS('KEY "R" (Return) with no number '//
     $   'to keep current (default) value$')
      ENDIF
C     First Set Required Parameters
C!!?? Make option to add custom parameters.  Make nparam=nparam+1... Back door?
      CALL SETREQ
1     CONTINUE
      NRLOOPS=2
C      IF(IFMGRID)NRLOOPS=NRLOOPS+1
      DO 40 IREQ=1,nrloops
       IF(IREQ.EQ.3)CALL SETREQ
       IF(IREQ.EQ.1)CALL PRS(' REQUIRED PARAMETERS$')
       IF(IREQ.EQ.2)CALL PRS(' OPTIONAL PARAMETERS$')
C       IF(IREQ.EQ.3)CALL PRS(' MULTIGRID PARAMETERS$')
C       IF(IREQ.EQ.3)THEN
c          IF(PARAM(37).EQ.0.0)THEN
C            Mulitgrid Parameters Haven't been set; Make defaults
C             IF(NORDER.LT.5)CALL PRS('ERROR; you need NORDER at least'//
C     $       '5 to use mulitgrid$')
C             IF(NORDER.LT.7)NGRIDS=2
C             IF(NORDER.GE.7)NGRIDS=3
C             param(37)=NGRIDS
C             IF(NGRIDS.eq.2)then
C                PARAM(38)=NORDER/2+1
C             ELSE if(ngrids.eq.3)THEN
C                PARAM(38)=NORDER/2+2
C                PARAM(39)=NORDER/2
C             ENDIF
C          ENDIF
C       ENDIF
       DO 40 LOOP=1,3
         IBOX=0
         NLOOPS=NPARAM+NPSCAL*2
         DO 5 I=1,NLOOPS
            IEQ=EQTYPE
            IF((IREQ.EQ.1.AND.(REQD(I).EQ.'R'.OR.REQD(I).EQ.'P'))
     $     .OR.(IREQ.EQ.2.AND. REQD(I).EQ.'O')
     $     .OR.(IREQ.EQ.3.AND. REQD(I).EQ.'M')
     $     .OR.(IREQ.EQ.1.AND.I.GT.NPARAM) ) THEN
C              It's Required
               IBOX =IBOX+1
C              Pretend it's a menu: use drmenu to Draw Boxes
               XCH=XLMEN+.175
               YCH=(YBS(IBOX)+YTS(IBOX))/2 - 0.006
               IF(I.GT.NPARAM)THEN
C                 Setting Extra Passive scalars
                  IPSCAL=I-NPARAM
                  IF(I.GT.NPARAM+NPSCAL)IPSCAL=I-(NPARAM+NPSCAL)
                  IF(LOOP.EQ.1) THEN
                     IF(I.GT.NPARAM+NPSCAL) THEN
C                       RHOCP
                        LINE='RHOCP'
                     ELSE
                        LINE='CONDUCT'
                     ENDIF
                     WRITE(LINE(9:9),'(I1)')IPSCAL
                     item(iBOX)=LINE
                  ELSE IF(LOOP.EQ.2) THEN
C                    WRITE VALUES
                     IF(I.GT.NPARAM+NPSCAL) THEN
C                       RHOCP
                        write(BOXLAB,'(G14.6)')PRHOCP(IPSCAL+2)
                     ELSE
                        write(BOXLAB,'(G14.6)')PCOND(IPSCAL+2)
                     ENDIF
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  ENDIF
               ELSE
                  IF(LOOP.EQ.1) item(iBOX)=cparam(i)
                  IF(LOOP.EQ.2) THEN
C                    WRITE VALUES
                     IF(I.LE.NPARAM)THEN
                        write(BOXLAB,'(G14.6)')PARAM(I)
                     ELSE IF(I.GT.NPARAM+NPSCAL) THEN
                        write(BOXLAB,'(G14.6)')PRHOCP(IPSCAL+2)
                        CALL PRSIR('IPSCAL,PRHOCP$'
     $                  ,IPSCAL,PRHOCP(IPSCAL))
                     ELSE IF(I.LE.NPARAM+NPSCAL) THEN
                        write(BOXLAB,'(G14.6)')PCOND (IPSCAL+2)
                        CALL PRSIR('IPSCAL,PCOND $'
     $                  ,IPSCAL,PCOND(IPSCAL+2))
                     ENDIF
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  ENDIF
               ENDIF
               IF(LOOP.EQ.3) THEN
C                 Blink Current box
                  CALL COLOR(5)
                  IF(IFBWGKS)CALL COLOR(0)
                  CALL MOVESC(XLMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YBS(IBOX))
C                 Change Values
301               continue
                  IF(IF4107.OR.IFXWIN) CALL RES(LINE,70)
                  IF(IFLEARN)WRITE(3,'(A70)') LINE
                  IF(IFSUN) CALL KEYPAD(VAL)
                  IF(IFGKS) CALL KEYPAD(VAL)
                  IF(LINE.EQ.'      ') THEN
C                    Do nothing
                  ELSE
                     IF(IF4107.OR.IFXWIN)CALL READER(VAL,IERR)
                     IF(IFSUN) IERR=0
                     IF(IERR.NE.0)THEN
                        CALL PRS(' The line you typed in:$')
                        CALL PRS(LINE//'$')
                        CALL PRS(' Was not understood.  Try again.$')
                        GO TO 301
                     ENDIF
C                    Successfully got new value; cover old one &Write it in box
                     CALL COLOR(2)
                     CALL FILLP(-4)
                     CALL BEGINB (XPHY(XCH),YPHY(YBS(IBOX)))
                     CALL DRAWSC (XRMEN ,YBS(IBOX))
                     CALL DRAWSC (XRMEN ,YTS(IBOX))
                     CALL DRAWSC (XCH,YTS(IBOX))
                     CALL DRAWSC (XCH,YBS(IBOX))
                     CALL ENDP
                     IF(I.GT.NPARAM+NPSCAL) THEN
                        PRHOCP(IPSCAL+2)=val
                     ELSE IF(I.GT.NPARAM) THEN
                        PCOND (IPSCAL+2)=val
                     ELSE
                        PARAM(I)=VAL
                     ENDIF
                     write(BOXLAB,'(G14.6)')VAL
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  ENDIF
C                 Check if new value is reasonable
C                 ??!! Make more extensive tests
                  IF(REQD(I).EQ.'R' .AND.PARAM(I).EQ.0.0) THEN
                     CALL PRS
     $               (' **ERROR** '//CPARAM(I)//'Must be NONZERO$')
                     CALL PRS(S//'$')
                     CALL PRS(' Try again: type in a nonzero value$')
                     CALL BEEP
                     GO TO 301
                  ELSE IF(REQD(I).EQ.'P' .AND.PARAM(I).LE.0.0) THEN
                     CALL PRS
     $               (' **ERROR** '//CPARAM(I)//'Must be POSITIVE$')
                    CALL PRS(S//'$')
                    CALL PRS(' Try again: type in a POSITIVE value$')
                    CALL BEEP
                    GO TO 301
                  ENDIF
                  IF(CPARAM(I).EQ.'GRID'.AND.GRID.GE.0.5) then
                     CALL PRS(' Grid should be less than 0.5$')
                     CALL PRS
     $               (' Grid > 0.5 gives only 4 possible points.$')
                     GO TO 301
                  ENDIF
                  NORDER=PARAM(20)
                  IF (NORDER.GT.NXM) THEN
                     WRITE(S,104) NX,NXM
  104                FORMAT(' Warning, current NORDER exceeds NXM'
     $                     ,' resetting from',I3,' to',I3,'.$')
                     CALL PRS(S)
                     NORDER=NXM
                     PARAM(20)=NORDER
                  ENDIF
                  IF(CPARAM(I).EQ.'NORDER'.AND.NORDER.GE.NXM)THEN
                     CALL PRSI(' ERROR: Maximum NORDER set to $',NXM)
                     CALL PRS(' Please contact your NEKTONICS '//
     $               'representative to increase it.$')
                     CALL BEEP
                     GO TO 301
                  ENDIF
C                 UN-Blink Current box
                  CALL COLOR(1)
                  CALL MOVESC(XLMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YBS(IBOX))
               ENDIF
C              Check that Value is Legal
               IF(IREQ.EQ.2 .AND.VAL.EQ.0.0) THEN
                  CALL PRS(' $')
C                 Now Re-do box with new (or default) Value
               endif
            ENDIF
5        CONTINUE
         NCHOIC=IBOX
         IF(LOOP.EQ.1)CALL DRMENU('SET PARAMETER')
         IF(LOOP.EQ.3.AND.IFNOSEG)CALL DRCOVR(13)
40    CONTINUE
      CALL CHKPAR(IERR)
      IF(IERR.NE.0) GO TO 1
      CALL GNABLE
C     Now get rid (permanently) of the stuff in segment 10 (the gwrited stuff)
      CALL CLSSEG
      CALL DELSEG(10)
      CALL OPNSEG(10)
      RETURN
      END
      SUBROUTINE CKPAR(MODE,IERR)
      INCLUDE 'basics.inc'
C     This was stolen from above setpar routine
      CHARACTER*26 BOXLAB
      CHARACTER*5 MODE
      DIMENSION ISHOW(20)
      CALL SETREQ
      IERR=0
      IEQTYPE=PARAM(16)
      IF(MODE.EQ.'SHOW ')THEN
         DO 4 IREQ=1,2
            IF(IREQ.EQ.1)CALL PRS(' *** REQUIRED PARAMETERS ***$')
            IF(IREQ.EQ.2)CALL PRS(' *** OPTIONAL PARAMETERS ***$')
            NSHOW=0
            DO 3 I=1,NPARAM
               IF( (REQD(I).EQ.'R'.OR.REQD(I).EQ.'P')
     $         .AND.IREQ.EQ.1) THEN
                  NSHOW=NSHOW+1
                  ISHOW(NSHOW)=I
               ELSE IF(REQD(I).EQ.'O' .AND.IREQ.EQ.2)THEN
                  NSHOW=NSHOW+1
                  ISHOW(NSHOW)=I
               ENDIF
3           CONTINUE
            IF(NSHOW.GT.0)THEN
               NLOOPS=((NSHOW-1)/3)+1
               DO 1 ILOOP=1,NLOOPS
                  II1=1+(ILOOP-1)*3
                  II2=II1+2
                  IF(II2.GT.NSHOW)II2=NSHOW
                  WRITE(S,'(3(1X,A10,1X,G14.6))')
     $            (CPARAM(ISHOW(II)),PARAM(ISHOW(II)),II=II1,II2)
               CALL PRS(S//'$')
 1             CONTINUE
            ENDIF
4        CONTINUE
c     ELSE IF(MODE.EQ.'CHECK')THEN
c        DO 5 I=1,NPARAM
c           IF(REQD(I).EQ.'R' .AND.PARAM(I).EQ.0.0) THEN
c              CALL PRS(' **ERROR** '//CPARAM(I)//'Must be NONZERO$')
c              CALL PRS(' Returning to beginning of parameter menu$')
c              CALL PRS(' Go back and enter a nonzero value$')
c              CALL BEEP
c              IERR=1
c           ELSE IF(REQD(I).EQ.'P' .AND.PARAM(I).LE.0.0) THEN
c pff          CALL PRS(' **ERROR** '//CPARAM(I)//'Must be POSITIVE$')
c annoyed      CALL PRS(' Returning to beginning of parameter menu$')
c 3/12/98      CALL PRS(' Go back and enter a POSITIVE value$')
c              CALL BEEP
c              IERR=1
c           ENDIF
5        CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE READLN(LOCATE)
      INCLUDE 'basics.inc'
C
C       Reads a line of text.  Jumps to appropriate place on help.
C
C       Writes out line to read file; don't worry, if you don't like it
C       you can backspace the read file and write a cleaner line
C
C
      CHARACTER HFILE*20,HLINE*80
C
C1     READ(IN,'(A70)',END=50,ERR=50) LINE
1     CALL RES(LINE,70)
      IF(IFLEARN)WRITE(3,'(A70)') LINE
      IF(IN.EQ.9) CALL PRS(LINE//'$')
      IF(LOCATE.NE.9999)CALL CAPIT(LINE,1)
2     CONTINUE
C
C     Don't Write out line in readfile.  Can backspace
      WRITE (X13,'(A70)') LINE
c      WRITE(10,'(A70)',ERR=60)  LINE
      RETURN
50    IN=5
      GO TO 1
60    CALL EREXIT
      END
c-----------------------------------------------------------------------
      SUBROUTINE PREXIT
      INCLUDE 'basics.inc'
      CHARACTER CTEMP*80,CHAR1*1,CHTEMP*3
c      LOGICAL IFMVBD
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      COMMON/INOUT/  IEXT
C
      CALL CLEARA
      CALL WRTPAR('FULL DUMP ')
      sesion(11:14) ='   '
      WRITE(10,'(4G14.6,'' XFAC,YFAC,XZERO,YZERO'')')
     $XFACO,YFACO,XZEROO,YZEROO
      IF(.NOT.IF3D)WRITE(10,*)
     $'**MESH DATA** 1st line is X of corner 1,2,3,4. 2nd line is Y.'
      IF(     IF3D)WRITE(10,*)
     $'**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8'
      NELV=NEL-NCOND
C
C  formatted/unformatted rea file  8/1/93 pff
C
c     CALL PRS('Formatted / Unformatted .rea file? (1/0)$')
c     CALL REI(Ioopt)
      Ioopt = 1
                      IFFMTIN=.TRUE.
      IF (Ioopt.eq.0) IFFMTIN=.FALSE.
C
      M=IEXT
      n=m+3
      filenm = sesion
C
      IF (.not.IFFMTIN) THEN
         FILENM(M:N) ='.re2'
         write(6,*) 'opening unformatted file ',filenm
         OPEN(UNIT=11,FILE=FILENM,STATUS='NEW',FORM='UNFORMATTED')
      ENDIF
      NELsgn=NEL
      IF (.not.IFFMTIN) NELsgn = -NEL
      WRITE(10,'(3I9,'' NEL,NDIM,NELV'')')NELsgn,NDIM,NELV
C
      DO 98 IEL=1,NEL
         IF (IEL.GT.52) LETAPT(IEL) = 'A'
         IF (IFFMTIN) THEN
C           formatted
            WRITE(10,'(A15,I9,A2,I5,A1,A10,I6)')
     $   '      ELEMENT  ',IEL,' [',NUMAPT(IEL),LETAPT(IEL),
     $   ']    GROUP',IGROUP(IEL)
            IF(NDIM.EQ.2)THEN
               WRITE(10,'(4G15.7)',ERR=60)(X(IEL,IC),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IEL,IC),IC=1,4)
            ELSE IF(NDIM.EQ.3)THEN
               WRITE(10,'(4G15.7)',ERR=60)(X(IEL,IC),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IEL,IC),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Z(IEL,IC),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(X(IEL,IC),IC=5,8)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IEL,IC),IC=5,8)
               WRITE(10,'(4G15.7)',ERR=60)(Z(IEL,IC),IC=5,8)
            ENDIF
         ELSE
C           Unformatted
            WRITE(11) IGROUP(IEL)
            IF(NDIM.EQ.2)THEN
               WRITE(11,ERR=60)(X(IEL,IC),IC=1,4)
               WRITE(11,ERR=60)(Y(IEL,IC),IC=1,4)
            ELSE IF(NDIM.EQ.3)THEN
               WRITE(11,ERR=60)(X(IEL,IC),IC=1,4)
               WRITE(11,ERR=60)(Y(IEL,IC),IC=1,4)
               WRITE(11,ERR=60)(Z(IEL,IC),IC=1,4)
               WRITE(11,ERR=60)(X(IEL,IC),IC=5,8)
               WRITE(11,ERR=60)(Y(IEL,IC),IC=5,8)
               WRITE(11,ERR=60)(Z(IEL,IC),IC=5,8)
            ENDIF
         ENDIF
C
98    CONTINUE
      IF (IFFMTIN) THEN
         WRITE(10,*)' ***** CURVED SIDE DATA *****'
         WRITE(10,'(I6,A20,A33)')NCURVE,' Curved sides follow',
     $   ' IEDGE,IEL,CURVE(I),I=1,5, CCURVE'
         IF(NCURVE.GT.0)THEN
            DO 19 IEL=1,NEL
            DO 19 IEDGE=1,8
               IF(CCURVE(IEDGE,IEL).NE.' ')THEN
                  IF (NEL.LT.1000) THEN
                     WRITE(10,'(I3,I3,5G14.6,1X,A1)')IEDGE,IEL,
     $               (CURVE(I,IEDGE,IEL),I=1,5),CCURVE(IEDGE,IEL)
                  ELSE
                     WRITE(10,'(I2,I6,5G14.6,1X,A1)')IEDGE,IEL,
     $               (CURVE(I,IEDGE,IEL),I=1,5),CCURVE(IEDGE,IEL)
                  ENDIF
               ENDIF
19          CONTINUE
         ENDIF
C
      ELSE
C
C Unformatted read
C
         WRITE(11) NCURVE
         IF(NCURVE.GT.0)THEN
            DO 119 IEL=1,NEL
            DO 119 IEDGE=1,8
               IF(CCURVE(IEDGE,IEL).NE.' ')THEN
                  WRITE(11) IEDGE,IEL,
     $            (CURVE(I,IEDGE,IEL),I=1,5),CCURVE(IEDGE,IEL)
               ENDIF
119         CONTINUE
         ENDIF
      ENDIF
C
      NSIDES=4
      IF(IF3D)NSIDES=6
      IF (IFFMTIN) WRITE(10,*)' ***** BOUNDARY CONDITIONS *****'
      NN=2
      IF(NFLDS.GT.2)NN=NFLDS
      DO 86 IFLD=1,NN
         IP=IFLD-2
         IF( (IFLD.EQ.1.AND.IFFLOW).OR.(IFLD.EQ.2.AND.IFHEAT).OR.
     $   IFLD.GT.2 )THEN
            IF(IFLD.EQ.1 .and. IFFMTIN )WRITE(10,*)
     $      ' ***** FLUID   BOUNDARY CONDITIONS *****'
            IF(IFLD.EQ.2 .and. IFFMTIN )write(10,*)
     $      ' ***** THERMAL BOUNDARY CONDITIONS *****'
            IF(IFLD.GT.2 .and. IFFMTIN )write(10,*)
     $      ' ***** PASSIVE SCALAR',ip,' BOUNDARY CONDITIONS *****'
C           Fluid and/or thermal
C           !!?? NELF DIFFERENT FROM NEL??
            DO 85 IEL=1,NEL
               DO 85 ISIDE=1,NSIDES
                  CHTEMP='   '
                  IF(IFLD.EQ.1 .OR. (IFLD.EQ.2 .AND. .NOT. IFFLOW))
     $            CHTEMP = CBC(ISIDE,IEL,0)
                  IF (NEL.LT.1000.and.IFFMTIN) THEN
                     WRITE(10,'(A1,A3,2I3,5G14.6)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  ELSEIF (IFFMTIN) THEN
                     WRITE(10,'(A1,A3,I5,I1,5G14.6)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  ELSE
                     WRITE(11,ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  ENDIF
                  ICBC=ICHAR(CBC(ISIDE,IEL,IFLD))
                  IF(ICBC.GE.97 .AND. ICBC.LE.122 )THEN
C                    Small letter for fortran b.c.'s
C                    Extra line(s) for Inflow.  BC(1) has # of lines; BC(2) addr
                     CBC3=CBC(ISIDE,IEL,IFLD)
                     IF(CBC3(3:3) .EQ. 'i') THEN
C                       Special storage of pointers in Internal boundaries
                        NLINES=BC(4,ISIDE,IEL,IFLD)
                        LINE1 =BC(5,ISIDE,IEL,IFLD)
                     ELSE
                        NLINES=BC(1,ISIDE,IEL,IFLD)
                        LINE1 =BC(2,ISIDE,IEL,IFLD)
                     ENDIF
                     IF (IFFMTIN) THEN
                        DO 82 I=1,NLINES
                           WRITE(10,'(A70)',ERR=60)INBC(LINE1+I-1)
82                      CONTINUE
                     ENDIF
                  ENDIF
85          CONTINUE
         ELSE
C           NO B.C.'s for this field
            IF(IFLD.EQ.1 .and. IFFMTIN )WRITE(10,*)
     $      ' ***** NO FLUID   BOUNDARY CONDITIONS *****'
            IF(IFLD.EQ.2 .and. IFFMTIN )write(10,*)
     $      ' ***** NO THERMAL BOUNDARY CONDITIONS *****'
         ENDIF
86    CONTINUE
C
      NSKIP=NLINP+NLINR
      WRITE(10,*)NSKIP,' PRESOLVE/RESTART OPTIONS  *****'
      DO 234 I=1,NSKIP
         WRITE(10,'(A80)')INITP(I)
234   CONTINUE
      NSKIP=NFLDS+3
      IF(NSKIP.LT.7)NSKIP=7
      WRITE(10,*)NSKIP,'         INITIAL CONDITIONS *****'
      DO 235 I=1,NSKIP
         WRITE(10,'(A80)')INITC(I)
235   CONTINUE
      WRITE(10,*)
     $' ***** DRIVE FORCE DATA ***** BODY FORCE, FLOW, Q'
      NSKIP=(4+NPSCAL)
      WRITE(10,*)NSKIP,
     $'                 Lines of Drive force data follow'
      WRITE(10,'(A80)')(DRIVC(I),I=1,4+NPSCAL)
C     First determine which fields have which mesh
      WRITE(10,*)' ***** Variable Property Data *****',
     $' Overrrides Parameter data.'
      NFLDSC=0
      DO 117 IFLD=2,NFLDS
         IF(IFTMSH(IFLD))NFLDSC=NFLDSC+1
117   CONTINUE
      IF1=1
      IF(.NOT.IFFLOW)IF1=2
      MFLDS=NFLDS
      IF(.NOT.IFFLOW) MFLDS=NFLDS-1
      NPACKS=0
      DO 182 IF=IF1,NFLDS
         DO 182 I=-5,10
            IF(MATYPE(I,IF).NE.0)THEN
               NPACKS=NPACKS+1
            ENDIF
182   CONTINUE
      NSKIP=NPACKS*4 +1
      WRITE(10,*)NSKIP,' Lines follow.'
      WRITE(10,*)NPACKS,' PACKETS OF DATA FOLLOW'
         DO 18 IGRP=-5,10
            DO 118 IF=IF1,NFLDS
               IF(MATYPE(IGRP,IF).NE.0)THEN
                  ITYPE=2
                  CTEMP=VPROP(IGRP,IF,1)
                  IF(CTEMP(1:1).EQ.'C')ITYPE=1
                  WRITE(10,'(3I5,A44)') IGRP,IF,MATYPE(IGRP,IF),
     $            '  GROUP, FIELD, TYPE OF DATA (1=CONST;2=FTN)'
                  DO 115 IPROP=1,3
                     IF(MATYPE(IGRP,IF).EQ.1)
     $               WRITE(10,'(G14.6)')CPROP(IGRP,IF,IPROP)
                     IF(MATYPE(IGRP,IF).EQ.2)
     $               WRITE(10,'(A80)'  )VPROP(IGRP,IF,IPROP)
115               CONTINUE
               ENDIF
118         CONTINUE
18       CONTINUE
C
      WRITE(10,*)' ***** HISTORY AND INTEGRAL DATA ***** '
      WRITE(10,*)NHIS,'   POINTS.  Hcode, I,J,H,IEL'
C     Sort so that integrals are last
      IF(NHIS.GT.0)THEN
         DO 50 ITYPE=1,2
         DO 50 I=1,NHIS
            IF(LOCHIS(3,I).EQ.0) LOCHIS(3,I)=1
            IF(LOCHIS(1,I).NE.0.AND.ITYPE.EQ.1 .OR.
     $         LOCHIS(1,I).EQ.0.AND.ITYPE.EQ.2 )
     $      WRITE(10,'(1X,11A1,1X,4I5)')
     $      (HCODE(II,I),II=1,11),(LOCHIS(I2,I),I2=1,4)
50       CONTINUE
      ENDIF
      WRITE(10,*)' ***** OUTPUT FIELD SPECIFICATION *****'
      NOUTS=6+NPSCAL
      WRITE(10,*)NOUTS,' SPECIFICATIONS FOLLOW'
      WRITE(10,*)IFXYO,'      COORDINATES'
      WRITE(10,*)IFVO, '      VELOCITY'
      WRITE(10,*)IFPO, '      PRESSURE'
      WRITE(10,*)IFTO, '      TEMPERATURE'
      WRITE(10,*)IFTGO,'      TEMPERATURE GRADIENT'
      WRITE(10,*)NPSCAL,'      PASSIVE SCALARS'
      IF(NPSCAL.GT.0)THEN
         DO 61 I=1,NPSCAL
            WRITE(10,'(1X,L1,1X,A5)')IFPSCO(I),PSNAME(I)
61       CONTINUE
      ENDIF
      WRITE(10,*)' ***** OBJECT SPECIFICATION *****'
      WRITE(10,'(2X,I6,A17)')NSOBJS, ' Surface Objects '
      IF(NSOBJS .GT. 0)THEN
          DO 62 IOBJ=1,NSOBJS
             WRITE(10,'(4X,I4,A35,A20,A1)')NFACE(IOBJ),
     $       ' element faces for surface object "',SOBJ(IOBJ),'"'
             DO 63 IFACE=1,NFACE(IOBJ)
                WRITE(10,*)ILSURF(1,IOBJ,IFACE),ILSURF(2,IOBJ,IFACE)
63           CONTINUE
62        CONTINUE
      ENDIF
      WRITE(10,'(2X,I6,A17)')NVOBJS, ' Volume  Objects '
      WRITE(10,'(2X,I6,A17)')NEOBJS, ' Edge    Objects '
      WRITE(10,'(2X,I6,A17)')NPOBJS, ' Point   Objects '
C
      CLOSE(UNIT=10)
c
      return
c
60    CALL EREXIT
105   FORMAT(I6,3G15.7,' Con Ele#, cdivity,RhoCp,QVOL')
112   FORMAT(G15.7,5X,A10)
      END
c-----------------------------------------------------------------------
      subroutine session_exit
c
C     Exit session
C
      CALL PRS('Deleting tmp.* files.$')
      CALL DELTMP
C     write where it will stay
      PRINT*,' Exiting session ',sesion
      CALL DEVEX
      CALL EXITT
c
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE EREXIT
      CALL BEEP
      CALL PRS(' **** ERROR DURING WRITE *****$')
      CALL PRS(' Incomplete writing to files$')
      CALL PRS(' Check your personal Quota and $')
      CALL PRS(' and if the disk is full$')
      CALL DEVEX
      CALL EXITT
      END
      SUBROUTINE INIT
      INCLUDE 'basics.inc'
      CHARACTER HLBNM*40,OPSYS*10,TERMIN*10
      CHARACTER*10 ARGS(5)
      CHARACTER*10 s10
      COMMON /file_pref/ file_prefix
      character*80 file_prefix
      COMMON/INOUT/  IEXT
      COMMON /PFFLG/  ILGRNG,ISCND,IENTRP,INEWTX
C
      XPHY0=0.0
      YPHY0=0.0
      ILGRNG=0
      ISCND =0
      IENTRP=0
c
      if_unstructure = .false.
C
      IFGRID = .TRUE.
      IFGRDC = .TRUE.
      IFGRDH = .FALSE.
      IFGRDP = .FALSE.
      IFOBJS = .FALSE.
      DO 1 I=1,20
         IFOBJG(I)=.FALSE.
    1 CONTINUE
      NOBJS = 0
C
      DO 2 IEL=1,NELM
         DO 2 IFIELD=1,MAXFLD
              MASKEL(IEL,IFIELD) = 1
2     CONTINUE
C
      XFAC=1.0
      YFAC=1.0
      CALL DEVINI('GENERIC')
      CALL SETENV
      CALL NOHARD
      Call CLRMAP('MESH')
C
C     Draw big red box
c
      CALL CLEAR
      CALL CLSSEG
      CALL OPNSEG(9)
      XL=1.0
      XR=1.3
      YT=0.99
      YB=0.01
      CALL COLOR(1)
      CALL FILLP(-2)
      CALL BEGINB (XPHY(XL),YPHY(YB))
      CALL DRAWSC (XR,YB)
      CALL DRAWSC (XR,YT)
      CALL DRAWSC (XL,YT)
      CALL DRAWSC (XL,YB)
      CALL ENDP
C
C     Open segment for Plots
c
      CALL GINDIS
C
5     call blank(sesion,14)
      CALL PRS(' Choose a Name for This Session:$')
      CALL RES(LINE,70)
      IF(line.eq.'!')THEN
         CALL PRS('In DEMO Mode$')
         IFDEMO=.true.
         GO TO 5
      ENDIF
      IF(line.eq.'$')THEN
         CALL PRS('In Learning Mode$')
         IFLEARN=.true.
         OPEN(UNIT=3,FILE='DEMO.COM',STATUS='NEW')
         WRITE(3,*)'$R NEKTON:PRENEK'
         WRITE(3,*)'!'
         WRITE(3,*)'DEMO'
         GO TO 5
      ENDIF
      IF(LINE.NE.'          ') THEN
         sesion=line
         OPEN(UNIT=4,FILE='session.name',STATUS='NEW',ERR=6)
         write(4,'(a14)')sesion
      ELSE
C        NO DEFAULT ON PRENEK
         GO TO 5
      ENDIF
6     CONTINUE
C
      CALL CLSSEG
      CALL OPNSEG(10)
C     WRITE ON SCREEN
C     Draw Keypad
      CALL DRKEY
C     KEYPAD OFF; COVER ON
      CALL SGVIS(12,0)
      CALL SGVIS(11,1)
C
      call blank(file_prefix,80)
      call chcopy(file_prefix,sesion,14)
c
      WRITE(S,'(A20,A10)') ' Beginning Session ',sesion
      CALL PRS(S//'$')
C     Check How many letters in sesion name
      do 7 i=10,1,-1
                if(sesion(i:i).ne.' ')then
                   lastch=i
                   go to 8
                endif
7     continue
8     continue
      do 9 i=1,lastch
          int=ichar(sesion(i:i))
          if(int.ge.65 .and. int.le.90) int=int+32
          FILENM(I:I)=CHAR(INT)
9     continue
      m=lastch+1
      IEXT=M
      n=m+3
C     Check for previous version of read file
      call findfl(sesion,oldver)
      sesion(11:14)='    '
c
C     Save of sesion to REAd later and open DRAwing
      FILENM(M:N) ='.dra'
      CALL OPENF(45,FILENM,'NEW',3,IERR)
C
      IF(IFLASE )CALL INITQMS
      IF(IFPOSTS)CALL INITPS
c      sesion(11:14)='$'
      call hard
      CALL GWRITE(1.03,0.05,1.3,sesion//'$')
      CALL GWRITE(XPHY(0.51),YPHY(0.95),1.5,'NEKTON 2.6$')
      CALL GWRITE(XPHY(0.71),YPHY(.95),1.5,DATE(5:16)//DATE(24:28)//'$')
      call nohard
C                                      Parameters
      DO 917  IJKLMN=1,500
         CALL BLANK(CPARAM(IJKLMN),40)
  917 CONTINUE
C
      CPARAM( 1)='DENSITY'
      CPARAM( 2)='VISCOS'
      CPARAM( 3)='BETAG'
      CPARAM( 4)='GTHETA'
      CPARAM( 5)='PGRADX'
      CPARAM( 6)='FLOWRATE'
      CPARAM( 7)='RHOCP'
      CPARAM( 8)='CONDUCT'
      CPARAM( 9)='QVOL'
      CPARAM(10)='FINTIME'
      CPARAM(11)='NSTEPS'
      CPARAM(12)='DT'
      CPARAM(13)='IOCOMM'
      CPARAM(14)='IOTIME'
      CPARAM(15)='IOSTEP'
C! INTERVAL BETWEEN STEPS
      CPARAM(16)='EQTYPE'
C! 1=N.S; 2=NS NO HEAT 3=STOKES NO HEAT 4=FORCED
C                             CONV W STEADY VEL; 5=TRANSIENT COND; 6=SSCOND
      CPARAM(17)='AXIS'
C! 0=CART; 1=AXISYMMETRIC   Y=R X=X
      CPARAM(18)='GRID'
      CPARAM(19)='INTYPE'
      CPARAM(20)='NORDER'
      CPARAM(21)='DIVERGENCE'
      CPARAM(22)='HELMHOLTZ'
      CPARAM(23)='NPSCAL'
C
      CPARAM(24)='TOLREL'
      CPARAM(25)='TOLABS'
      CPARAM(26)='COURANT'
      CPARAM(27)='TORDER'
C
      CPARAM(36)='XMAGNET'
      CPARAM(37)='NGRIDS'
      CPARAM(38)='NORDER2'
      CPARAM(39)='NORDER3'
      CPARAM(40)='NORDER4'
      CPARAM(41)='NORDER5'
      CPARAM(52)='HISTEP'
      CPARAM(53)='HISTIME'
C
      NPARAM=53
      DO 10 I=1,NPARAM
C          DEFAULTS
           PARAM(I)=0.0
C
           call chcopy(s10,cparam(i),10)
           IF( s10 .EQ.'RHOCP     ')PARAM(I) = 1
           IF( s10 .EQ.'CONDUCT   ')PARAM(I) = 1
           IF( s10 .EQ.'VISCOS    ')PARAM(I) = 1
           IF( s10 .EQ.'NSTEPS    ')PARAM(I) = 1
           IF( s10 .EQ.'TORDER    ')PARAM(I) = 2
           IF( s10 .EQ.'COURANT   ')PARAM(I) = 0.25
           IF( s10 .EQ.'DENSITY   ')PARAM(I) = 1
           IF( s10 .EQ.'IOCOMM    ')PARAM(I) = 20
           IF( s10 .EQ.'INTYPE    '.AND.NKTONV.EQ.1)PARAM(I) = 1
           IF( s10 .EQ.'INTYPE    '.AND.NKTONV.EQ.2)PARAM(I) =-1
           IF( s10 .EQ.'NORDER    ')PARAM(I) = 5
           IF( s10 .EQ.'GRID      ')PARAM(I) = 0.05
           IF( s10 .EQ.'DIVERGENCE')PARAM(I) = 1.0E-2
           IF( s10 .EQ.'HELMHOLTZ ')PARAM(I) = 1.0E-4
           IF( s10 .EQ.'DT        ')PARAM(I) = 1.0
           IF( s10 .EQ.'TOLREL    ')PARAM(I) = 0.01
           IF( s10 .EQ.'TOLABS    ')PARAM(I) = 0.01
10    CONTINUE
      DO 20 I=1,11
         PCOND (I)=1.0
         PRHOCP(I)=1.0
20    CONTINUE
      CWRITE = 1
      CALL SCROLL(32)
C
C                                       Commands IS THIS NECESSARY?
      CMD( 1)='  '
      CMD( 2)='SET'
      CMD( 3)='SHOW'
      CMD( 4)='EXIT'
      CMD( 5)='QUIT'
      CMD( 6)='READ'
      CMD( 7)='HELP'
      CMD( 8)='BUILD'
      CMD( 9)='HISTORY'
      CMD(10)='FAMILY'
      CMD(11)='ABORT'
      CMD(12)='INITCOND'
      CMD(13)='DRIVFORC'
      CMD(14)='OUTPUT'
      NCMD=14
C
C     Initialize transformation matrix
      THETAr=THETA*PI/180.0
      PHIr  =PHI  *PI/180.0
C
C     Direction toward eye (normal vector)
      VHOBS(1)=COS(THETAR)*COS(PHIR)
      VHOBS(2)=COS(THETAR)*SIN(PHIR)
      VHOBS(3)=SIN(THETAR)
C       Perpendicular on x-y plane
      XHOBS(1)=-1.0*SIN(PHIR)
      XHOBS(2)=COS(PHIR)
      XHOBS(3)=0.0
C       Perpendicular to above two
      YHOBS(1)=-1.0*SIN(THETAR)*COS(PHIR)
      YHOBS(2)=-1.0*SIN(THETAR)*SIN(PHIR)
      YHOBS(3)=COS(THETAR)
C
      RETURN
      END
      SUBROUTINE OUTPUT
      INCLUDE 'basics.inc'
      CHARACTER KEY,STRING*5,NEWNAM*5
      CHARACTER*26 ITEMON
C! DEFAULT
C       Sets the parameters set in the field dumps
C
      IF(NEL.EQ.0)THEN
              CALL PRS(' Do Outputs AFTER you Build mesh.$')
              RETURN
      ENDIF
      CALL PRS(' Choose menu item to toggle storage on and off: $')
C     MMUST DECIDE HOW TO SET DEFAULTS BASED ON OCODES
      ITEM(1)=               'MAIN MENU'
      IF(     IFXYO) ITEM(2)='COORDINATES  YES'
      IF(.NOT.IFXYO) ITEM(2)='COORDINATES  NO '
1     NCHOIC=2
      IF(IFFLOW)THEN
         ITEMON='VELOCITY'
         IF(      IFVO) ITEMON(14:16)='YES'
         IF(.NOT. IFVO) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
C
         ITEMON='PRESSURE'
         IF(      IFPO) ITEMON(14:16)='YES'
         IF(.NOT. IFPO) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
      IF(IFHEAT)THEN
         ITEMON='TEMPERATURE'
         IF(      IFTO) ITEMON(14:16)='YES'
         IF(.NOT. IFTO) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
C
c         ITEMON='TGRADIENT'
c         IF(      IFTGO) ITEMON(14:16)='YES'
c         IF(.NOT. IFTGO) ITEMON(14:16)='NO '
c         NCHOIC=NCHOIC+1
c         ITEM(NCHOIC)=ITEMON
      ENDIF
      NNOPS=NCHOIC
      IF(NPSCAL.GT.0)THEN
         DO 5 IPSCAL=1,NPSCAL
C            ITEMON='PASSIVE SCALAR  '
            ITEMON=PSNAME(IPSCAL)
            IF(      IFPSCO(IPSCAL)) ITEMON(14:16)='YES'
            IF(.NOT. IFPSCO(IPSCAL)) ITEMON(14:16)='NO '
            NCHOIC=NCHOIC+1
            ITEM(NCHOIC)=ITEMON
5        CONTINUE
      ENDIF
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'OUTPUT')
      IF(CHOICE.NE.'MAIN MENU') THEN
C        Toggle
         IPSCAL=ICHOIC-nnops
         ITEMON=ITEM(ICHOIC)
         IF(ITEMON(14:16).EQ.'YES')then
            ITEMON(14:16)='NO '
            IF(IPSCAL .GT.0)IFPSCO(IPSCAL) = .FALSE.
         ELSE IF(ITEMON(14:16).EQ.'NO ')THEN
            ITEMON(14:16)='YES'
            IF(IPSCAL .GT.0)IFPSCO(IPSCAL) = .TRUE.
            IF(IPSCAL .GT.0)THEN
C              Passive Scalar
               CALL PRSIS
     $         ('Current label for passive scalar $',IPSCAL,' :$')
               CALL PRS(' '//PSNAME(IPSCAL)//'$')
               CALL PRS
     $         ('<CR> to keep or type in new name <= 5 Characters$')
               CALL RES(NEWNAM,5)
               IF(IFLEARN)WRITE(3,'(A5)')newnam
               IF(NEWNAM.NE.'     ')PSNAME(IPSCAL)=NEWNAM
               ITEMON(1:5)=PSNAME(IPSCAL)
            ENDIF
         ENDIF
         ITEM(ICHOIC)=ITEMON
         GO TO 1
      ENDIF
C     Before returning to main menu,
      IF(NPSCAL.GT.0)THEN
         IPSAVE=0
         DO 101 I=1,NPSCAL
           IF(IFPSCO(I))IPSAVE=IPSAVE+1
101      CONTINUE
         IF(IPSAVE.GT.3)THEN
            CALL PRS(
     $      ' *** ERROR *** MAXIMUM OF 3 PASSIVE SCALARS CAN BE SAVED$')
            CALL PRS('Please turn off some of the output switches$')
            GO TO 1
         ENDIF
      ENDIF
      DO 10 I=2,NCHOIC
         ITEMON=ITEM(I)
         IF(ITEMON(1:13).EQ.'COORDINATES')THEN
            IF(NKTONV.EQ.1)THEN
               CALL PRS('ERROR: Coordinates always saved in version 1$')
            ELSE
               IF(ITEMON(14:16).EQ.'YES')IFXYO=.TRUE.
               IF(ITEMON(14:16).EQ.'NO ')IFXYO=.FALSE.
            ENDIF
         ELSE IF(ITEMON(1:13).EQ.'VELOCITY')THEN
            IF(ITEMON(14:16).EQ.'YES')IFVO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFVO=.FALSE.
         ELSE IF(ITEMON(1:13).EQ.'PRESSURE')THEN
            IF(ITEMON(14:16).EQ.'YES')IFPO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPO=.FALSE.
         ELSE IF(ITEMON(1:13).EQ.'TEMPERATURE')THEN
            IF(ITEMON(14:16).EQ.'YES')IFTO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTO=.FALSE.
         ELSE IF(ITEMON(1:13).EQ.'TGRADIENT')THEN
            IF(ITEMON(14:16).EQ.'YES')IFTGO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTGO=.FALSE.
         ENDIF
9     CONTINUE
10    CONTINUE
      RETURN
      END
      SUBROUTINE HISTRY
      INCLUDE 'basics.inc'
      CHARACTER*26 ITEMON
C     Local logicals
      LOGICAL IFPART,IFPNT,IFVH,IFTH,IFPH,IFPS1H,IFPS2H
     $,IFPS3H,IFPS4H,IFXH
      CHARACTER KEY,STRING*5
      CALL PRS('                  *** HISTORIES ***$')
C
      IFPART=.FALSE.
      IFPNT =.TRUE.
C     Default history switches to output switches
      IFVH = IFVO
      IFXH = .NOT. IFVH
      IFTH = IFTO
      IFPH = IFPO
      IFXYH= IFXYO
      IPSCH= IPSCO
C     Hilite historical points (If any are already input)
      IF(NHIS.GT.0)THEN
         CALL PENW(7)
         CALL COLOR(3)
         DO 5 I=1,NHIS
            IH  = LOCHIS(1,I)
            JH  = LOCHIS(2,I)
            IELH= LOCHIS(4,I)
            IF(HCODE(10,I).EQ.'H'.OR.HCODE(10,I).EQ.'P')
     $      CALL XDOT(XPTS(IH,JH,IELH),YPTS(IH,JH,IELH))
5        CONTINUE
         CALL COLOR(1)
         CALL PENW(3)
      ENDIF
C
1     CONTINUE
C
      CALL PRS(' Choose menu item to toggle storage on and off: $')
C     MMUST DECIDE HOW TO SET DEFAULTS BASED ON OCODES
      ITEM(1)='MAIN MENU'
      ITEM(2)='ENTER POINT'
      NCHOIC=2
      IF(NLEVEL.GT.1)THEN
         nchoic = nchoic+1
         ITEM(nchoic)='UP   LEVEL'
         nchoic = nchoic+1
         ITEM(nchoic)='DOWN LEVEL'
      ENDIF
      IF(IFFLOW)THEN
         ITEMON='VELOCITY'
         IF(      IFVH) ITEMON(14:16)='YES'
         IF(.NOT. IFVH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
C
C
         ITEMON='PRESSURE'
         IF(      IFPH) ITEMON(14:16)='YES'
         IF(.NOT. IFPH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
      IF(IFMVBD)THEN
         ITEMON='MESH LOCATION'
         IF(      IFXH) ITEMON(14:16)='YES'
         IF(.NOT. IFXH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
      IF(IFHEAT)THEN
         ITEMON='TEMPERATURE'
         IF(      IFTH) ITEMON(14:16)='YES'
         IF(.NOT. IFTH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
      IF(NPSCAL.GT.0)THEN
         ITEMON='PASSIVE SCALAR  '
         WRITE(ITEMON(16:16),'(I1)')IPSCH
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
C
      IF(NKTONV.EQ.2 .AND. IFFLOW)THEN
         ITEMON='STILL POINT'
         IF(      IFPNT) ITEMON(14:16)='YES'
         IF(.NOT. IFPNT) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
C
         ITEMON='TRACK PARTCL'
         IF(      IFPART) ITEMON(14:16)='YES'
         IF(.NOT. IFPART) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      ENDIF
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'HISTORY')
      IF(CHOICE.EQ.'UP   LEVEL' .OR.CHOICE.EQ.'DOWN LEVEL')THEN
C        Erase old mesh& draw new
         IF(CHOICE.EQ. 'UP   LEVEL'.AND.ILEVEL.EQ.NLEVEL)THEN
            CALL PRS(' ERROR: AT TOP LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         ENDIF
         IF(CHOICE.EQ. 'DOWN LEVEL'.AND.ILEVEL.EQ.1)THEN
            CALL PRS(' ERROR: AT FIRST LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         ENDIF
         IF(CHOICE.EQ.'UP   LEVEL')
     $   CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
         IF(CHOICE.EQ.'DOWN LEVEL')
     $   CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
         DO 500 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
500      CONTINUE
         IF(CHOICE.EQ. 'UP   LEVEL')ILEVEL=ILEVEL+1
         IF(CHOICE.EQ. 'DOWN LEVEL')ILEVEL=ILEVEL-1
         DO 510 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL)THEN
            CALL DRAWEL( I)
            ENDIF
510      CONTINUE
      else IF(CHOICE.NE.'MAIN MENU'.AND.CHOICE.NE.'ENTER POINT') THEN
C        Toggle
         ITEMON=ITEM(ICHOIC)
         IF(ITEMON(14:16).EQ.'YES')then
            ITEMON(14:16)='NO '
         ELSE IF(ITEMON(14:16).EQ.'NO ')THEN
            ITEMON(14:16)='YES'
         ENDIF
         ITEM(ICHOIC)=ITEMON
C
         DO 10 I=3,NCHOIC
          ITEMON=ITEM(I)
          IF(ITEMON(1:13).EQ.'COORDINATES')THEN
            IF(NKTONV.EQ.1)THEN
               CALL PRS('ERROR: Coordinates always saved in version 1$')
            ELSE
               IF(ITEMON(14:16).EQ.'YES')IFXYH=.TRUE.
               IF(ITEMON(14:16).EQ.'NO ')IFXYH=.FALSE.
            ENDIF
          ELSE IF(ITEMON(1:13).EQ.'VELOCITY')THEN
            IF(ITEMON(14:16).EQ.'YES')IFVH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFVH=.FALSE.
            IFXH = .NOT. IFVH
          ELSE IF(ITEMON(1:13).EQ.'MESH LOCATION')THEN
            IF(ITEMON(14:16).EQ.'YES')IFXH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFXH=.FALSE.
            IFVH = .NOT. IFXH
          ELSE IF(ITEMON(1:13).EQ.'PRESSURE')THEN
            IF(ITEMON(14:16).EQ.'YES')IFPH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPH=.FALSE.
          ELSE IF(ITEMON(1:13).EQ.'TEMPERATURE')THEN
            IF(ITEMON(14:16).EQ.'YES')IFTH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTH=.FALSE.
          ELSE IF(ITEMON(1:13).EQ.'STILL POINT')THEN
            IF(ITEMON(14:16).EQ.'YES')IFPNT=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPNT=.FALSE.
          ELSE IF(ITEMON(1:13).EQ.'TRACK PARTCL')THEN
            IF(ITEMON(14:16).EQ.'YES')IFPART=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPART=.FALSE.
          ELSE IF(ITEMON(1:14).EQ.'PASSIVE SCALAR')THEN
            IPSCH=IPSCH+1
            IF(IPSCH .GT. NPSCAL) IPSCH=0
            WRITE(ITEMON(16:16),'(I1)')IPSCH
          ENDIF
10       CONTINUE
         GO TO 1
      ELSE IF(  CHOICE.EQ.'MAIN MENU') THEN
         RETURN
      ELSE IF(CHOICE.EQ.'ENTER POINT')THEN
C       How do you enter 3-d points from preprocessor?
        CALL PRS(' Enter x,y coordinates with mouse:$')
        CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
        IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C          He apparently is trying to use the keypad
           Call Prs(' Enter X-coordinate with keypad:$')
           CALL KEYPAD(XMOUSE)
           Call Prs(' Now enter Y-coordinate with keypad:$')
           CALL KEYPAD(YMOUSE)
        ENDIF
        IF(IF3D)THEN
           CALL PRS(' Enter z coordinate.  Use keypad.$')
           CALL KEYPAD(ZMOUSE)
        ENDIF
C       Now that we have everything, put it in arrays.
        NHIS=NHIS+1
C       Character flags
        DO 15 I=1,11
           HCODE(I,NHIS)=' '
15      CONTINUE
        IF(IFVH)THEN
          HCODE(1,NHIS)='U'
          HCODE(2,NHIS)='V'
          IF(IF3D)HCODE(3,NHIS)='W'
        ENDIF
        IF(IFXH)THEN
          HCODE(1,NHIS)='X'
          HCODE(2,NHIS)='Y'
          IF(IF3D)HCODE(3,NHIS)='Z'
        ENDIF
        IF(IFPH) HCODE(4,NHIS)='P'
        IF(IFTH) HCODE(5,NHIS)='T'
        WRITE(HCODE(6,NHIS),'(I1)')IPSCH
        IF(IPSCH.EQ.0) HCODE(6,NHIS)=' '
C       Currently only historical points supported
        IF(IFPNT ) HCODE(10,NHIS)='H'
        IF(IFPART) HCODE(10,NHIS)='P'
C       Use closest collocation point
        RMIN=1.0E6
        DO 250 IEL=1,NEL
          DO 250 I=1,NX
            DO 250 J=1,NY
              DO 250 K=1,NZ
C               WHAT'S XPTS??
                IF(.NOT. IF3D) THEN
                   R=SQRT((XPTS(I,J,IEL)-XMOUSE)**2+
     $                    (YPTS(I,J,IEL)-YMOUSE)**2)
                ELSE
                   R=SQRT((XPTS(I,J,IEL)-XMOUSE)**2+
     $                    (YPTS(I,J,IEL)-YMOUSE)**2+
     $             (ZMOUSE- ( (Z(IEL,5)+Z(IEL,1))/2+
     $          (Z(IEL,5)-Z(IEL,1))/2*ZPTS(K) ) )**2)
                ENDIF
                IF(R.LT.RMIN) THEN
                    RMIN=R
                    IH=I
                    JH=J
                    KH=K
                    IELH=IEL
                ENDIF
250     CONTINUE
        LOCHIS(1,NHIS)=IH
        LOCHIS(2,NHIS)=JH
        LOCHIS(3,NHIS)=KH
        LOCHIS(4,NHIS)=IELH
C       Hilite historical point
        CALL COLOR(3)
        CALL PENW(7)
        CALL XDOT(XPTS(IH,JH,IELH),YPTS(IH,JH,IELH))
        CALL COLOR(1)
        CALL PENW(3)
        GO TO 1
      ENDIF
      GO TO 1
100   FORMAT(3I5,4X,4a1,
     $  '  I,J,IEL, and CODES of Point for History(3I5,4X,A4)')
      END
      SUBROUTINE INTEGQ
      INCLUDE 'basics.inc'
1     CONTINUE
      DO 18 ISOBJ=1,NSOBJS
         ITEM(ISOBJ)=SOBJ(ISOBJ)
         DO 17 IHIS=1,NHIS
            IF(LOCHIS(1,IHIS) .EQ.ISOBJ) ITEM(ISOBJ) = 'I '//SOBJ(ISOBJ)
17       CONTINUE
18    CONTINUE
      NCHOIC=NSOBJS
      NCHOIC=NCHOIC+1
C
      ITEM(NCHOIC)='MAIN MENU'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'INTEGRAL')
      IF(CHOICE.EQ.'MAIN MENU')RETURN
C
      ISOBJ = ICHOIC
      DO 5 IHIS=1,NHIS
         IF(LOCHIS(1,IHIS) .EQ.ISOBJ)THEN
            CALL BEEP
            CALL PRS('Object '//SOBJ(ISOBJ)//' Already Chosen $')
            GO TO 1
         ENDIF
5     CONTINUE
      CALL PRS('Adding '//SOBJ(ISOBJ)//'$')
      NHIS=NHIS+1
C     Character flags
      DO 15 I=1,11
           HCODE(I,NHIS)=' '
15    CONTINUE
      IF(IFFLOW)THEN
         HCODE(1,NHIS)='F'
         HCODE(2,NHIS)='F'
         IF(IF3D)HCODE(3,NHIS)='F'
      ENDIF
      IF(IFHEAT)HCODE(5,NHIS)='Q'
      DO 20 IPSCAL=1,NPSCAL
         IF(IPSCAL.LE.4)WRITE(HCODE(IPSCAL+5,NHIS),'(I1)')IPSCAL
20    CONTINUE
      HCODE(10,NHIS)='I'
      LOCHIS(1,NHIS)=ISOBJ
      LOCHIS(2,NHIS)=0
      LOCHIS(3,NHIS)=0
      LOCHIS(4,NHIS)=0
      GO TO 1
      END
      SUBROUTINE OBJEC
      INCLUDE 'basics.inc'
      CHARACTER*26 ITEMON
      CHARACTER*80 LINE80,ITEMD
      CHARACTER CHFACE*1
C     Local logicals
      LOGICAL IFPART,IFPNT,IFVH,IFTH,IFPH,IFPS1H,IFPS2H
     $,IFPS3H,IFPS4H
      CHARACTER KEY,STRING*5
      CALL PRS('                  *** SURFACE OBJECTS ***$')
      CALL PRS(
     $' Define a surface object $')
C
C     Hilite Current surface objects
      IF(NSOBJS .EQ. 0)CALL PRS(' No Surface Objects Yet Defined$')
      IF(NSOBJS .GT. 0)CALL PRS(' Current Surface Objects Defined:$')
C
      DO 5 IOBJ=1,NSOBJS
         CALL PRSI(' '//SOBJ(IOBJ)//' $',IOBJ)
         CALL PRS(S//'$')
         DO 6 IFACE=1,NFACE(IOBJ)
C           Annotate Object
            IEL  = ILSURF(1,IOBJ,IFACE)
            ISIDE= ILSURF(2,IOBJ,IFACE)
            WRITE(CHFACE,'(I1)',ERR=1)IOBJ
            XSTAR = SIDES(IEL,ISIDE,1)
            YSTAR = SIDES(IEL,ISIDE,2)
            CALL GWRITE(XSTAR,YSTAR,1.0,'o'//CHFACE//'$')
6        CONTINUE
5     CONTINUE
1     CONTINUE
C
      ITEM(1)='MAIN MENU'
      ITEM(2)='DELETE OBJECT'
      ITEM(3)='MODIFY OBJECT'
      ITEM(4)='ADD    OBJECT'
      NCHOIC=4
C
      IF(NLEVEL.GT.1)THEN
         NCHOIC = NCHOIC+1
         ITEM(NCHOIC)='UP   LEVEL'
         NCHOIC = NCHOIC+1
         ITEM(NCHOIC)='DOWN LEVEL'
      ENDIF
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'OBJECT')
      IF(CHOICE.EQ.'UP   LEVEL' .OR.CHOICE.EQ.'DOWN LEVEL')THEN
C        Erase old mesh& draw new
         IF(CHOICE.EQ. 'UP   LEVEL'.AND.ILEVEL.EQ.NLEVEL)THEN
            CALL PRS(' ERROR: AT TOP LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         ENDIF
         IF(CHOICE.EQ. 'DOWN LEVEL'.AND.ILEVEL.EQ.1)THEN
            CALL PRS(' ERROR: AT FIRST LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         ENDIF
         IF(CHOICE.EQ.'UP   LEVEL')
     $   CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
         IF(CHOICE.EQ.'DOWN LEVEL')
     $   CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
         DO 500 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
500      CONTINUE
         IF(CHOICE.EQ. 'UP   LEVEL')ILEVEL=ILEVEL+1
         IF(CHOICE.EQ. 'DOWN LEVEL')ILEVEL=ILEVEL-1
         DO 510 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL)THEN
            CALL DRAWEL( I)
            ENDIF
510      CONTINUE
      ELSE IF(  CHOICE.EQ.'MAIN MENU') THEN
         RETURN
      ELSE IF(CHOICE.EQ.'MODIFY OBJECT'
     $.OR.    CHOICE.EQ.'ADD    OBJECT')THEN
C        NSOBJS:        Number of objects
C        ISOBJ :        Current Object
C        NFACE(ISOBJ): Number of faces in current object
C
         IF(CHOICE.EQ.'ADD    OBJECT')THEN
            CALL PRS('Enter Object Name (20 characters or less):$')
            NSOBJS=NSOBJS+1
            ISOBJ =NSOBJS
            NFACE(ISOBJ) = 0
            CALL RES(SOBJ(ISOBJ),20)
         ENDIF
         IF(CHOICE.EQ.'MODIFY OBJECT')THEN
            DO 11 IOBJ=1,NSOBJS
              ITEM(IOBJ)=SOBJ(IOBJ)
11          CONTINUE
            NCHOIC=NSOBJS
            CALL MENU(XMOUSE,YMOUSE,BUTTON,'MODIFY OBJECT')
            ISOBJ = ICHOIC
         ENDIF
C
2        CONTINUE
         ITEM(1)='END OBJECT'
         ITEM(2)='ADD    FACE'
         ITEM(3)='DELETE FACE'
         NCHOIC=3
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'ADD    OBJECT')
         IF(CHOICE.EQ.'END OBJECT')THEN
            GO TO 1
         ELSE IF(CHOICE.EQ.'ADD    FACE')THEN
            CALL GETFAC(IFACE,IEL,'ADD   ',ISOBJ)
            NFACE(ISOBJ) = NFACE(ISOBJ) + 1
            ILSURF(1,ISOBJ,NFACE(ISOBJ)) = IEL
            ILSURF(2,ISOBJ,NFACE(ISOBJ)) = IFACE
         ELSE IF(CHOICE.EQ.'DELETE FACE')THEN
            CALL GETFAC(IFACE,IEL,'DELETE',ISOBJ)
            DO 12 I=1,NFACE(ISOBJ)
               IF(ILSURF(1,ISOBJ,I) .EQ. IEL   .AND.
     $            ILSURF(2,ISOBJ,I) .EQ. IFACE )THEN
C                 Copy highest face into this storage location and reduce number
C                 of faces by one
                  ILSURF(1,ISOBJ,I)=ILSURF(1,ISOBJ,NFACE(ISOBJ))
                  ILSURF(2,ISOBJ,I)=ILSURF(2,ISOBJ,NFACE(ISOBJ))
                  NFACE(ISOBJ)= NFACE(ISOBJ)-1
                  GO TO 13
               ENDIF
12          CONTINUE
            CALL PRS('Error: Selected face does not match object face$')
13          CONTINUE
         ENDIF
         GO TO 2
      ELSE IF(CHOICE.EQ.'DELETE OBJECT')THEN
         DO 18 IOBJ=1,NSOBJS
           ITEM(IOBJ)=SOBJ(IOBJ)
18       CONTINUE
         NCHOIC=NSOBJS
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='DO NOT DELETE'
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'DELETE OBJECT')
         IF(CHOICE.NE.'DO NOT DELETE')THEN
            ISOBJ = ICHOIC
            CALL PRS('Deleting '//SOBJ(ISOBJ)//'$')
            CALL DELOBJ(ISOBJ)
         ENDIF
      ENDIF
      GO TO 1
      END
      SUBROUTINE DELOBJ(ISOBJT)
      INCLUDE 'basics.inc'
      CHARACTER CHFACE
C
C     Label new number of former highest; blank location of former ISOBJ
      DO 1 IFACE = 1,NFACE(NSOBJS)
         IEL   = ILSURF(1,NSOBJS,IFACE)
         ISIDE = ILSURF(2,NSOBJS,IFACE)
         XSTAR = SIDES(IEL,ISIDE,1)
         YSTAR = SIDES(IEL,ISIDE,2)
         WRITE(CHFACE,'(I1)',ERR=11)ISOBJ
11       CALL GWRITE(XSTAR,YSTAR,1.0,'o'//CHFACE//'$')
1     CONTINUE
      DO 2 IFACE = 1,NFACE(ISOBJ)
         IEL   = ILSURF(1,ISOBJ,IFACE)
         ISIDE = ILSURF(2,ISOBJ,IFACE)
         XSTAR = SIDES(IEL,ISIDE,1)
         YSTAR = SIDES(IEL,ISIDE,2)
         CALL GWRITE(XSTAR,YSTAR,1.0,'  $')
2     CONTINUE
C     Copy Highest object into ISOBJ; reduce NSOBJS by 1
      SOBJ (ISOBJ) = SOBJ (NSOBJS)
      NFACE(ISOBJ) = NFACE(NSOBJS)
      DO 3 IFACE = 1,NFACE(ISOBJ)
         ILSURF(1,ISOBJ,IFACE) = ILSURF(1,NSOBJS,IFACE)
         ILSURF(2,ISOBJ,IFACE) = ILSURF(2,NSOBJS,IFACE)
3     CONTINUE
C     Change pointer to highest object to ISOBJ; Warn if any had been set to
C     ISOBJ??
      DO 4 IHIS=1,NHIS
         IF(LOCHIS(1,IHIS).EQ.NSOBJS) LOCHIS(1,IHIS)=ISOBJ
4     CONTINUE
      NSOBJS = NSOBJS - 1
      RETURN
      END
      SUBROUTINE GETFAC(IFACE,IEL,CFLAG,ISOBJT)
      INCLUDE 'basics.inc'
      CHARACTER*6 CFLAG
      CHARACTER CHFACE
      CALL PRS(' Enter element side with mouse:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C          He apparently is trying to use the keypad
           Call Prs(' Enter X-coordinate with keypad:$')
           CALL KEYPAD(XMOUSE)
           Call Prs(' Now enter Y-coordinate with keypad:$')
           CALL KEYPAD(YMOUSE)
      ENDIF
C       Use closest side
        RMIN=1.0E6
        DO 250 IEL=1,NEL
           IF(NUMAPT(IEL).EQ.ILEVEL)THEN
              DO 240 I=1,NSIDES
C               WHAT'S XPTS??
                R=SQRT((SIDES(IEL,I,1)-XMOUSE)**2+
     $                 (SIDES(IEL,I,2)-YMOUSE)**2)
                IF(R.LT.RMIN) THEN
                    RMIN=R
                    IQSIDE= I
                    IQEL  = IEL
                    IF(IQSIDE.EQ.5.OR.IQSIDE.EQ.6)THEN
C                      CHECK IF IT'S ABOVE OR BELOW CENTER
                       IF(YMOUSE.GT.SIDES(IQEL,IQSIDE,2))THEN
                          IQSIDE=6
                       ELSE
                          IQSIDE=5
                       ENDIF
                    ENDIF
                ENDIF
240           CONTINUE
           ENDIF
250     CONTINUE
        IFACE=IQSIDE
        IEL  =IQEL
C       Hilite side
        CALL COLOR(4)
        IF(IQSIDE.EQ.5.OR.IQSIDE.EQ.6)THEN
           IF(IQSIDE.EQ.5)THEN
              IED1=1
              IED2=4
           ELSE IF(IQSIDE.EQ.6)THEN
              IED1=1
              IED2=4
           ENDIF
           CALL MOVE(X(IQEL,IED1),Y(IQEL,IED1))
           DO 230 IEDGE=IED1,IED2
              CALL DRAWED(IQEL,IEDGE,1)
230        CONTINUE
        ELSE
           CALL MOVE(X(IQEL,IQSIDE),Y(IQEL,IQSIDE))
           CALL DRAWED(IQEL,IQSIDE,1)
           WRITE(CHFACE,'(I1)',ERR=111)ISOBJ
111        CONTINUE
           XSTAR = SIDES(IQEL,IQSIDE,1)
           YSTAR = SIDES(IQEL,IQSIDE,2)
           IF(CFLAG.EQ.'ADD   ')
     $     CALL GWRITE(XSTAR,YSTAR,1.0,'o'//CHFACE//'$')
           IF(CFLAG.EQ.'DELETE')
     $     CALL GWRITE(XSTAR,YSTAR,1.0,'  $')
        ENDIF
        CALL COLOR(1)
        CALL PENW(3)
      RETURN
      END
      SUBROUTINE DRIVEF
      INCLUDE 'basics.inc'
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      IDRIVF=1
      CALL PRS('          *** DRIVE FORCE MENU ***$')
      IF(IFFLOW)THEN
C               You can drive the fluid
              IPSCAL=0
100           CONTINUE
              CALL PRS('Possible Fluid Driving Functions:$')
              CALL PRS('B to specify Body Force'//
     $        ' (Forcing function)$')
              CALL PRS('N for No Forcing Function$')
              CALL PRS(' Fluid Driver >$')
              CALL READLN(142)
              READ(X13,'(A70)') LINE
              CALL CAPIT(ANS,1)
              IF(ANS.EQ.'B'.OR.ANS.EQ.'N')THEN
                CALL DRIVFC(IPSCAL)
              ELSE
                CALL PRS('TRY AGAIN; B OR N$')
                GO TO 100
              ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE DRIVFC(IPSCAL)
      INCLUDE 'basics.inc'
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      character*80 CH
      DATA CH /' '/
C
      IF(ANS.EQ.'N')RETURN
C         ZERO OUT DRIVING FORCES
          DRIVC(1)=' '
          DRIVC(2)=' '
          DRIVC(3)=' '
      IF(ANS.EQ.'P'.OR.ANS.EQ.'B')THEN
          IPFLAG=1
          if(nktonv.eq.1)then
          else if(nktonv.eq.2)then
              call prs
     $        ('Enter 3 fortran statements for FFX,FFY,FFZ.  Each$')
              call prs('Can be a function of x,y,z, and time.$')
              CALL PRS('Example:$')
              CALL PRS('FFX = SIN(TIME*10)*Y*Z$')
              CALL PRS('FFY = SIN(TIME*10)*X*Z$')
              IF(IF3D)CALL PRS('FFZ = SIN(TIME*10)*X*Y$')
              CALL PRS('FFX:$')
         ENDIF
      ENDIF
      IF(ANS.EQ.'F')THEN
          IFFLAG=1
          CALL PRS(
     $    'FLOW is total flow rate through the geometry and can be $')
          CALL PRS(
     $    'a function of TIME. It is input in the style of a fortran $')
          CALL PRS('statement.  Start in column 1.  Use one line.$')
          CALL PRS('Example:$')
          CALL PRS('FLOW = 1.0 + 1.0/(TIME + 1.0)$')
          CALL PRS('Enter statement For FLOW:$')
      ENDIF
      IF(ANS.EQ.'Q')THEN
          IQFLAG=1
          IF(IPSCAL.GT.0)THEN
             CALL PRSI('PASSIVE SCALAR$',IPSCAL)
             CALL PRS(S//'$')
          ENDIF
          CALL PRS(
     $    'QVOL internal heat source/volume (or area).  It can be a$')
          CALL PRS(
     $    'function of X, Y,[Z], TIME, TEMP, and PASS. It is input as$')
          CALL PRS(
     $    'a fortran statement.  Start in column 1.  Use one line.$')
          CALL PRS('Example:$')
          CALL PRS('QVOL = (1.0 -X) * Y + EXP(-TIME) + TEMP$')
          CALL PRS('Enter statement For QVOL:$')
      ENDIF
C
C       Read in fortran line; Write to fortran subroutine
      CALL RES(CH(8:80),73)
      IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
      IF(ANS.EQ.'P'.OR.ANS.EQ.'B')DRIVC(1)=CH
      IF((ANS.EQ.'P'.OR.ANS.EQ.'B').AND.nktonv.eq.2)then
C        Read 2 more lines for forcing function
         CALL PRS('FFY:$')
         CALL RES(CH(8:80),73)
         IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
         DRIVC(2)=CH
         IF(IF3D)THEN
            CALL PRS('FFZ:$')
            CALL RES(CH(8:80),73)
            IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
            DRIVC(3)=CH
         ENDIF
      ENDIF
      IF(ANS.EQ.'F')DRIVC(4)=CH
      IF(ANS.EQ.'Q')DRIVC(4+IPSCAL)=CH
      RETURN
      END
C     End of PRENEK ****************************************
      SUBROUTINE GENNMG(PARAM)
      REAL PARAM(100)
      integer n(5)
      CHARACTER*80 SLOC
c
       JMAX = 5
       DO 10 J=1,JMAX
          N(J) = 0
 10    CONTINUE
       NORDER=PARAM(20)
       N(1) = NORDER
       DO 100 J=2,JMAX
          NPLOLD = N(J-1) - 1
          NPLNEW = NPLOLD/2
          IF (NPLNEW*2 .NE. NPLOLD) NPLNEW = NPLNEW+1
          IF (NPLNEW .EQ. 1) THEN
             NGRIDS = J-1
             GO TO 200
          ENDIF
          N(J) = NPLNEW+1
 100   CONTINUE
 200   CONTINUE
c
       PARAM(37)=NGRIDS
       DO 300 J=2,NGRIDS
          NPOLY = N(J)-1
          PARAM(36+J)=N(J)
          WRITE(SLOC,'(3i3)') J,NPOLY,N(J)
          CALL PRS(SLOC//'$')
 300   CONTINUE
      RETURN
      END
      SUBROUTINE INTCND
      INCLUDE 'basics.inc'
C
C     Initial condions menu
C
      CHARACTER*26 INITEM(5)
      SAVE        INITEM
      DATA        INITEM / 'MAIN MENU                 '
     $                    ,'DEFAULT=0                 '
     $                    ,'FUNCTION                  '
     $                    ,'PRESOLVE                  '
     $                    ,'RESTART                   '/
C
C
   10 CONTINUE
      NCHOIC=5
      CALL CHCOPY(ITEM,INITEM,200)
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'INIT. CONDS.')
C
      IF (CHOICE.EQ.INITEM(1)) RETURN
      IF (CHOICE.EQ.INITEM(2)) CALL INITDEF
      IF (CHOICE.EQ.INITEM(3)) CALL INITFUN
      IF (CHOICE.EQ.INITEM(4)) CALL INITPRE
      IF (CHOICE.EQ.INITEM(5)) CALL INITRES
      GOTO 10
C
      END
      SUBROUTINE INITFUN
      INCLUDE 'basics.inc'
      character*80 CH
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
C
      IF (NLINF.GT.0)
     $   CALL PRS(' Current FORTRAN function initial conditions:$')
      DO 100 IFLD=1,NLINF
            CALL PUTS(INITC(IFLD),80)
  100 CONTINUE
C
      CALL PRS
     $   (' Enter initial conditions for each field as a function$')
      CALL PRS
     $   (' of X and Y (Z) in the style of a Fortran statement.$')
C
      DO 1000 IFLD=1,NFLDS
C
         CALL BLANK(CH,80)
         CALL PRS
     $   (' Starting in column 1, use one line for $')
         LINE=' '
         IF(IFLD.EQ.1)THEN
              IF(IF3D)
     $        CALL PRS(' UX one line for UY and one for UZ.  Example:$')
              IF(.NOT.IF3D)
     $        CALL PRS(' UX and one for UY.  Example:$')
              CALL PRS(' UX = 0.667 * (1.0-Y**2)$')
              CALL PRS(' UY = 0.0$')
              IF(IF3D)CALL PRS('UZ = 0.0$')
              CALL PRS('  Note: Check FORTRAN syntax carefully.$')
         ELSE IF(IFLD.EQ.2)THEN
              CALL PRS(' TEMP.  Example:$')
              CALL PRS(' TEMP= 0.667 * (1.0-Y**2)$')
         ELSE
              IDUMMY=IFLD-2
              CALL PRiS(IDUMMY,'th PS.  Example:$')
              CALL PRS(' PS = 0.667 * (1.0-Y**2)$')
         ENDIF
C
         IF(IFLD.EQ.1)THEN
              CALL PRS(' Enter statement for UX:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(1)=CH
              CALL PRS(' Now enter statement for UY:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(2)=CH
              IF(NDIM.EQ.3)THEN
                 CALL PRS(' Now enter statement for UZ:$')
                 CALL RES(CH(8:80),73)
                 IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
                 INITC(3)=CH
              ENDIF
         ELSEIF (IFLD.EQ.2)THEN
              CALL PRS(' Enter statement for TEMP:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(IFLD+3)=CH
         else
              CALL PRS(' Enter statement for PS:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(IFLD+3)=CH
         ENDIF
 1000 CONTINUE
      NLINF=NFLDS
      RETURN
      END
      SUBROUTINE INITPRE
      INCLUDE 'basics.inc'
      CHARACTER*1 OPTS(80),INITP1(80)
      EQUIVALENCE (INITP1,INITP)
C
      NLIN=NLINP+NLINR
      IF (NLIN.GT.0) CALL PRS(' Current restart/presolve status:$')
      DO 100 I=1,NLIN
         CALL PUTS(INITP(I),80)
  100 CONTINUE
C
      IF (NLINP.EQ.0) THEN
         NLINP=1
         NLINR=MIN(NLINR,14)
         DO 200 IC=15,2,-1
            INITP(IC)=INITP(IC-1)
  200    CONTINUE
         INITP(1)='PRESOLVE'
C
         CALL PRS(' The default is to presolve all fields.$')
         IF(IFFLOW.AND.IFSPLIT) CALL PRS
     $ ('** Note:  For Split Scheme, velocity PRESOLV not applicable$')
  300    CALL PRS(' Enter <cr> to accept, field specification, or ?.$')
         CALL RES(OPTS,80)
         CALL LJUST(OPTS)
         CALL CAPIT(OPTS,80)
C
         IF (LTRUNC(OPTS,80).EQ.0) RETURN
         IF (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITPRE',1)
            GOTO 300
         ELSE
            CALL CHCOPY(INITP1(11),OPTS,70)
         ENDIF
      ENDIF
C
  400 CALL PRS(' Current presolve status:$')
      CALL PUTS(INITP1,80)
  500 CALL PRS
     $(' Enter <cr> to accept, M to modify, R to remove, or ?.$')
      CALL RES(OPTS,80)
      CALL LJUST(OPTS)
      CALL CAPIT(OPTS,80)
C
      IF (LTRUNC(OPTS,80).EQ.0) RETURN
C
      IF (OPTS(1).EQ.'R') THEN
         CALL PRS('Removing presolve option.$')
         NLINP=0
         DO 550 IC=1,14
            I1=IC+1
            CALL CHCOPY(INITP(IC),INITP(I1),80)
  550    CONTINUE
         CALL BLANK(INITP(15),80)
         RETURN
      ENDIF
C
      IF (OPTS(1).EQ.'M') THEN
  600    CALL PRS(' Enter field specification or ?.$')
         CALL RES(OPTS,80)
         IF (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITPRE',1)
            GOTO 600
         ELSE
            CALL CHCOPY(INITP1(11),OPTS,70)
         ENDIF
         GOTO 400
      ENDIF
C
      IF (OPTS(1).EQ.'?') THEN
         CALL HELPER('INITPRE',2)
         GOTO 400
      ENDIF
C
C     Else, invalid entry, re-prompt:
      GOTO 500
C
      END
      SUBROUTINE INITRES
      INCLUDE 'basics.inc'
      CHARACTER*80 FNAME
      CHARACTER*1  OPTS(80),INITP1(80,15)
      EQUIVALENCE (INITP1,INITP)
      LOGICAL IFADD,IFREMV
C
      IFREMV=.FALSE.
   10 CONTINUE
      NLIN=NLINP+NLINR
      IF (NLIN.GT.0) CALL PRS(' Current restart/presolve status:$')
      DO 100 I=1,NLIN
         CALL PUTS(INITP(I),80)
  100 CONTINUE
C
C     Add/review?
C
      IFADD=.FALSE.
      IF (NLINR.GT.0.OR.IFREMV) THEN
         CALL PRS(' Input A to add a file, M to modify/remove, $')
         CALL PRS(' or <cr> to return.$')
         CALL RES(OPTS,80)
         CALL PRS(' $')
         CALL LJUST(OPTS)
         CALL CAPIT(OPTS,80)
C
         IF (LTRUNC(OPTS,80).EQ.0) RETURN
         IF (OPTS(1).EQ.'A') IFADD=.TRUE.
      ENDIF
C
C     Add a restart file
C
      IF (IFADD.OR.NLINR.EQ.0) THEN
         CALL PRS(' Enter restart session or file name:$')
         CALL RES(FNAME,80)
         IF (LTRUNC(FNAME,80).EQ.0) RETURN
         NLINR=NLINR+1
         NLIN =NLINR+NLINP
         CALL LJUST(FNAME)
         INITP(NLIN)=FNAME
C
         CALL PRS(' The default is to restart all fields.$')
  300    CALL PRS(' Enter <cr> to accept, field specification, or ?.$')
         CALL RES(OPTS,80)
         CALL LJUST(OPTS)
         CALL CAPIT(OPTS,80)
C
         IF (LTRUNC(OPTS,80).EQ.0) RETURN
C
         IF (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITRES',1)
            GOTO 300
         ELSE
C           specifying restart options
            IOPT=LTRUNC(INITP(NLIN),80)+2
            LOPT=80-IOPT+1
            CALL CHCOPY(INITP1(IOPT,NLIN),OPTS,LOPT)
            RETURN
         ENDIF
      ENDIF
C
C     Alter/Remove restart requests.
C
      ILIN=1+NLINP
      IF (NLINR.GT.1) THEN
         CALL PRS(' Enter restart file # to be altered:$')
         CALL REI(ILIN)
      ENDIF
C
C     Alter/accept current restart request
C
  400 CALL PRS(' Current restart status:$')
      CALL PUTS(INITP(ILIN),80)
  500 CALL PRS
     $(' Enter <cr> to accept, M to modify, R to remove, or ?.$')
      CALL RES(OPTS,80)
      CALL LJUST(OPTS)
      CALL CAPIT(OPTS,80)
C
      IF (LTRUNC(OPTS,80).EQ.0) RETURN
C
      IF (OPTS(1).EQ.'?') THEN
         CALL HELPER('INITRES',2)
         GOTO 500
      ENDIF
C
      IF (OPTS(1).EQ.'M') THEN
  600    CALL PRS(' Enter field specification or ?.$')
         CALL RES(OPTS,80)
         IF (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITRES',1)
            GOTO 600
         ELSE
            CALL LJUST(INITP(ILIN))
            IOPT=INDX1(INITP(ILIN),' ',1)+1
            LOPT=80-IOPT+1
            CALL CHCOPY(INITP1(IOPT,ILIN),OPTS,LOPT)
         ENDIF
         GOTO 400
      ENDIF
C
      IF (OPTS(1).EQ.'R') THEN
         DO 700 I=ILIN,14
            I1=I+1
            CALL CHCOPY(INITP(I),INITP(I1),80)
  700    CONTINUE
         CALL BLANK(INITP(15),80)
         NLINR=NLINR-1
c        kludge city
         IFREMV=.TRUE.
         GOTO 10
      ENDIF
C
      RETURN
      END
      SUBROUTINE HELPER(ROUTINE,IENTRY)
      INCLUDE 'basics.inc'
      CHARACTER ROUTINE*(*)
C
      IF (ROUTINE.EQ.'INITPRE'.AND.IENTRY.EQ.1) THEN
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn.$')
        CALL PRS
     $  (' If presolve is requested without specification, all$')
        CALL PRS
     $  (' fields will start with a steady state solution as an$')
        CALL PRS
     $  (' initial condition.$')
        RETURN
      ENDIF
C
      IF (ROUTINE.EQ.'INITPRE'.AND.IENTRY.EQ.2) THEN
        CALL PRS
     $  (' M - allows you to modify the field specifications.$')
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn.$')
        CALL PRS(' $')
        CALL PRS
     $  (' R - remove all presolve requests.$')
        RETURN
      ENDIF
C
      IF (ROUTINE.EQ.'INITRES'.AND.IENTRY.EQ.1) THEN
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn,$')
        CALL PRS
     $  (' and/or TIME=ttt.tt, and/or d#.$')
        CALL PRS
     $  (' If no options are specified, the default is that all$')
        CALL PRS
     $  (' fields will start from the last full dump in the specified$')
        CALL PRS
     $  (' file. Initial value of TIME will be that of the last dump.$')
        CALL PRS(' $')
        CALL PRS
     $  (' Example:  session1        U$')
        CALL PRS
     $  ('           file1.old    T   TIME=10.0 $')
        CALL PRS
     $  (' ...would take the velocity from session1.fld and the$')
        CALL PRS
     $  ('    temperature from file1.old.  Initial time would be 10.$')
        RETURN
      ENDIF
C
      IF (ROUTINE.EQ.'INITRES'.AND.IENTRY.EQ.2) THEN
        CALL PRS
     $  (' M - allows you to modify the field specifications.$')
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn, d#,$')
        CALL PRS
     $  (' and TIME=ttt.tt (d# is the dump number, ttt.ttt is time)$')
        CALL PRS(' $')
        CALL PRS
     $  (' R - remove this restart request.$')
        RETURN
      ENDIF
C
      RETURN
      END
      SUBROUTINE INITDEF
      INCLUDE 'basics.inc'
      CHARACTER*1 OPTS(80),INITP1(80)
      EQUIVALENCE (INITP1,INITP)
C
      CALL PRS(' Set all initial conditions to zero? (y/n)$')
      CALL RES(OPTS,80)
      CALL LJUST(OPTS)
      CALL CAPIT(OPTS,80)
      IF (OPTS(1).NE.'Y') RETURN
C
      NLINR=0
      NLINP=0
      CALL BLANK(INITP,1200)
      CALL BLANK(INITC,1200)
      DO 10 I=1,15
         INITC(I)='C Default'
   10 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine exitt
      write(6,*) 'stopping in exitt'
      stop
      end
c-----------------------------------------------------------------------
      subroutine exit
      write(6,*) 'stopping in exit'
      stop
      end
c-----------------------------------------------------------------------
