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
c-----------------------------------------------------------------------
      subroutine fprep

c     PreProcessor for Spectral element code.  Inputs geometry and 
c     flow data from keyboard or graphics tablet interactively or 
c     from instruction file.  - Ed Bullister

      include 'basics.inc'

      character file*17
      common /cfilold/ filold
      character filold*17
      common/inout/  iext

      call data     !  Replaced Data stmts with assignment statements
      call init
      call gnable

      file=sesion

      nchoic = 3
      item(1)='READ PREVIOUS PARAMETERS'
      item(2)='TYPE IN  NEW  PARAMETERS'
      item(3)='CONJ. HEAT TRANSFER MERGE'

      ifcen=.true.
      call menu(xmouse,ymouse,button,'READ PARAMETER')

    1 continue

      if (choice.eq.'TYPE IN  NEW  PARAMETERS') then ! Interactive query
         in=5
         CALL PRS(' $')
         CALL PRS(' Assemble form of equation to be solved $')
         CALL PRS
     $   (' using mouse to toggle logical switches on and off.$')

         call seteqt !  First, set equation type
         call setpar !  Next, Set Required Parameters for this eqtype

      else

         if (choice.eq.'READ PREVIOUS PARAMETERS') then

            call prs(' Enter name of previous session$')

         elseif (choice.eq.'CONJ. HEAT TRANSFER MERGE') then

            ifmerge      = .true.
            ifconj_merge = .true.

            call prs(' Enter name of fluid session$')

         endif

         call blank(line,80)
         call res  (line,80)
         if (indx1(line,'.rea',4).eq.0) then ! Append .rea, if not present
            len = ltrunc(line,80)
            call chcopy (lines(len+1),'.rea',4)
         endif

         ifread=.true.
         filold=line
         call openf(9,filold,'OLD',3,ierr)

         if (ierr.ne.0) then
            call prs(' Can''t open file '//FILOLD//'   Try again.$')
            call prs(' $')
            in=5
            goto 1
         endif

         call read_params

      endif

      call build_options
      call build
      call after_build

      stop
      end
c-----------------------------------------------------------------------
      subroutine read_params
      include 'basics.inc'
      character*40 s40

C       Do this read twice to get parameter values and comments

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

        rewind(9)

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
C        IF(VNEKOLD.GE.2.4) then
C           READ(9,*,ERR=59)NTEXTSW
C           DO 143 I=1,NTEXTSW
C              READ(9,'(2A40)',err=59)TEXTSW(I,1),TEXTSW(I,2)
C143        CONTINUE
C        endif
        NFLDS=1
        IF(IFHEAT)NFLDS=NFLDS+1
        NFLDS=NFLDS+NPSCAL
C       Read Lots more stuff !!??
        CLOSE(UNIT=9)
C       Check that all req'd parameters are set to legal values
C       If something missing, then go to 300
        CALL CKPAR('CHECK',IERR)

      return

   59 call prs(' Error reading parameters from file.$')
      stop

      end
c-----------------------------------------------------------------------
      subroutine build_options
      include 'basics.inc'

      common /cfilold/ filold
      character filold*17

      if (ifmerge) then        ! No queries - just use existing file
         choice='BUILD FROM FILE'   ! For build routine
         filenm=filold
         call prs(' Will Read Mesh and B.C. data from  '//FILENM//'$')
         call openf(9,filenm,'old',3,ierr)
         return
      endif

310   NCHOIC =  2
      ITEM(1)='BUILD FROM FILE'
      ITEM(2)='BUILD INTERACTIVELY'
c     ITEM(2)='ALTER PARAMETERS'
c     ITEM(3)='SHOW PARAMETERS'
c     ITEM(4)='BUILD INTERACTIVELY'
c     ITEM(5)='IMPORT UNIVERSAL FILE'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'CENTRAL')


      IF (choice.eq.'SHOW PARAMETERS') THEN
C        Print Out Parameters
         CALL PRS('****** PARAMETERS *****$')
         IF(AXIS.EQ.0) then
            WRITE(S,'('' NEKTON VERSION '',G14.1,
     $      '' ;    '',I3,'' DIMENSIONAL '')')VNEKTON,NDIM
            CALL PRS(S//'$')
         else
            CALL PRSRS(' NEKTON VERSION $',VNEKTON,
     $      ' ; AXISYMMETRIC$')
            CALL PRS(S//'$')
         endif
         CALL CKPAR('SHOW ',IERR)
         GO TO 310
      elseif (choice.eq.'ALTER PARAMETERS') THEN

         call seteqt !  First, set equation type
         call setpar !  Next, Set Required Parameters for this eqtype

      endif
C
C     Make logicals rational
      IF(.NOT. IFFLOW) then
         IFNAV =.FALSE.
         IFSTRS=.FALSE.
         IFSPLIT=.FALSE.
      endif
      IF(.NOT. IFHEAT) then
         IFADVC(2)=.FALSE.
      endif
C
C     All parameters set; Ready to build
C     How does build know if we are interactive or reading?
C     On end-of-file on read file, will we always jump to the right place?
C     Do this interactively independently of whether we read in parameters

      IF(NDIM.EQ.2) then
          NSIDES=4
          NEDGES=4
      elseif (NDIM.EQ.3) then
          NSIDES=6
          NEDGES=8
      endif
      IF(choice.eq.'BUILD INTERACTIVELY') THEN
         IFREAD=.FALSE.
      elseif (choice.eq.'IMPORT UNIVERSAL FILE') THEN
         IFUNIV = .TRUE.
1070     CALL PRS(
     $   ' Enter full name (including extension) of universal file:$')
         CALL RES(LINE,70)
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         IFREAD=.TRUE.
         FILENM=LINE
         CALL OPENF(8,FILENM,'OLD',3,IERR)
         IF(IERR.NE.0) then
            CALL PRS(' Can''t open file '//FILENM//'   Try again.$')
            IN=5
            GOTO 310
         endif
         CALL REAUNI
      elseif (choice.eq.'BUILD FROM FILE') THEN
C        Read from file
1080     CALL PRS(' Enter name of previous session$')
         CALL PRS('(Default= '//FILOLD//')$')
         CALL RES(LINE,70)
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         IFREAD=.TRUE.
         if (LINE.NE.' ') THEN
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
         else
C           Use default file
            FILENM=FILOLD
         endif
         CALL PRS(' Will Read Mesh and B.C. data from  '//FILENM//'$')
         CALL OPENF(9,FILENM,'OLD',3,IERR)
         IF(IERR.NE.0) then
            CALL PRS(' Can''t open file $'//FILENM//'   Try again.$')
            IN=5
c           GO TO 1080
            GOTO 310
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine after_build
      include 'basics.inc'


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
      if (NLINE.GT.0) then
         CALL BEEP
         CALL BEEP
         CALL PRS('** CURRENT RESTART PRESOLVE OPTIONS **$')
         CALL PRS('   USE INITIAL COND TO CHANGE/UPDATE.$')
      endif
      DO 105 I=1,NLINE
         CALL PUTS(INITP(I),80)
  105 CONTINUE
      DO 106 I=1,NFLDS
         LINE=INITC(I)
         IF(I.GT.1)LINE=INITC(I+3)
         IF(LINE(1:9).NE.'C Default') then
            CALL PRSIS
     $      ('** NOTE: ** CURRENT INITIAL CONDITIONS FOR FIELD$'
     $,     i,'ARE$')
            CALL PRS(LINE//'$')
            CALL PRS(' USE INITIAL COND TO CHANGE/UPDATE.$')
         endif
106   CONTINUE
108   NCHOIC =  3
      ITEM(1)='EXIT'
      ITEM(2)='OUTPUT'
      ITEM(3)='DRIVE FORCE'


c     IF(IFTRAN) then             ! Commented out on 3/19/2010 pff
c        NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='HISTORY'
c     endif



C     The initial conditions are allowed to give a first guess at steady
C     temperature, or, to give a flow field for steady forced convection.
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='INITIAL COND'
      IF(IFTRAN .AND. IFHEAT .AND. IFADVC(2) .AND. (.NOT.IFFLOW)) then
C        This kludge for velocity field for forced convection calls
C        the initial condition menu
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='FORCING VELOCITY'
      endif

      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='INTEGRAL QUANTITY'
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='OBJECT'
      NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='ZOOM'

c     NCHOIC=NCHOIC+1
c     ITEM(NCHOIC)='RSB'


      if (ifconj_merge) then
         choice = item(1)
      else
         call menu(xmouse,ymouse,button,'OPTIONS')
      endif

      IF(choice.eq.'EXIT') then
C        Exit Prenek
         call prexit
c        call rsb_xxt_set(2,nel)
         call session_exit
c     elseif(choice.eq.'RSB') then
c        call rsb_xxt_set(1,nel)
      elseif (choice.eq.'HISTORY') then
         CALL HISTRY
      elseif (choice.eq.'INTEGRAL QUANTITY') then
         CALL INTEGQ
      elseif (choice.eq.'OUTPUT') then
         CALL OUTPUT
      elseif (choice.eq.'OBJECT') then
         CALL OBJEC
      elseif (choice.eq.'DRIVE FORCE') then
         CALL DRIVEF
      elseif (choice.eq.'ZOOM') then
         CALL SETZOOM
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         DO 175 IEL=1,NEL
            CALL DRAWEL(IEL)
  175    CONTINUE
      elseif (choice.eq.'INITIAL COND') then
         CALL INTCND
      elseif (choice.eq.'FORCING VELOCITY') then
         CALL PRS
     $   ('Forcing velocity is put in via the initial conditions$')
         CALL PRS
     $   ('menu.  This velocity will remain constant with time.$')
         CALL INTCND
      endif
      GO TO 108
59    CALL PRS(' ERROR READING PARAMETERS FROM FILE.$')
      stop
60    CALL PRS(' ERROR WRITING; CHECK DISK QUOTA$')
      stop
      end
c-----------------------------------------------------------------------
      subroutine reauni
      include 'basics.inc'
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
         IF(IFSTRT) then
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
            endif
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
                  IF(IF3D .AND. NCORNS .NE. 8) then
                     CALL PRSIS
     $               ('ERROR: ELEMENT HAS $',NCORNS,' CORNERS$')
                     CALL PRS('3-D ELEMENTS MUST HAVE 8 CORNERS$')
                     STOP
                  elseif ((.NOT. IF3D) .AND. NCORNS .NE. 4) then
                     CALL PRSIS
     $               ('ERROR: ELEMENT HAS $',NCORNS,' CORNERS$')
                     CALL PRS('2-D ELEMENTS MUST HAVE 4 CORNERS$')
                     STOP
                  endif
                  READ(8,'(A80)',ERR=13,END=13)LINE80
                  READ(LINE80,'(8I10)',ERR=14,END=14)
     $            (IICORN(IC),IC=1,NCORNS)
C
                  DO 50 IC=1,NCORNS
                     x(ic,iel) = xpts(iicorn(ic),1,1)
                     y(ic,iel) = ypts(iicorn(ic),1,1)
                     z(ic,iel) = zptn(iicorn(ic),1,1)
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
            endif
         endif
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
      endif
      XZERO = RBARX(X,NELM,NEL,NCORNS) - XFAC/2.
      YZERO = RBARX(Y,NELM,NEL,NCORNS) - YFAC/2.
      return
14    CONTINUE
      CALL PRS('Error Reading Internal Data$')
      stop
      end
c-----------------------------------------------------------------------
      function rmaxx(r,nelm,nel,nc)
      real r(8,nelm)
      rm=r(1,1)
      do 1 i=1,nel
      do 1 j=1,nc
         if (r(j,i) .gt. rm) rm=r(j,i)
1     continue
      rmaxx = rm
      return
      end
c-----------------------------------------------------------------------
      function rminx(r,nelm,nel,nc)
      real r(8,nelm)
      rm=r(1,1)
      do 1 i=1,nel
      do 1 j=1,nc
         if (r(j,i) .lt. rm) rm=r(j,i)
1     continue
      rminx = rm
      return
      end
c-----------------------------------------------------------------------
      function rbarx(r,nelm,nel,nc)
      real r(8,nelm)
      rb=0.
      do 1 i=1,nel
      do 1 j=1,nc
         rb = rb+r(i,j)
1     continue
      rbarx = rb
      return
      end
c-----------------------------------------------------------------------
      function rmax(r,n)
      REAL R(N)
      RM=R(1)
      DO 1 I=1,N
         IF(R(I) .GT. RM) RM=R(I)
1     CONTINUE
      RMAX = RM
      return
      end
c-----------------------------------------------------------------------
      function rmin(r,n)
      REAL R(N)
      RM=R(1)
      DO 1 I=1,N
         IF(R(I) .LT. RM) RM=R(I)
1     CONTINUE
      RMIN = RM
      return
      end
c-----------------------------------------------------------------------
      function rbar(r,n)
      REAL R(N)
      RB =0.0
      DO 1 I=1,N
         RB = RB + R(I)
1     CONTINUE
      RBAR = RB / N
      return
      end
c-----------------------------------------------------------------------
      subroutine wrtpar(cflag)
      include 'basics.inc'

      CHARACTER FILE*10,CFLAG*10
      CHARACTER*1 s401(40)
      COMMON/INOUT/  IEXT
C
C     Write out parameter stuff
C     First check B.C.'s to set logical switches
      IFMVBD    = .FALSE.
      IFTMSH(0) = .FALSE.
      IF(CFLAG.EQ.'FULL DUMP ') then
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
      endif
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
         write(10,1049) param(ip),ip,(s401(j),j=1,l)
 1049    format(g14.6,1X,'p',i3.3,' ',40a1)
 1050 CONTINUE
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
      return
60    CALL PRS(' ERROR WRITING; CHECK DISK QUOTA$')
      STOP
      end
c-----------------------------------------------------------------------
      subroutine seteqt
      include 'basics.inc'
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
      IF(IFSTRS) then
         LSET(ISET:ISET+9)=' STRESSFUL'
         ISET=ISET+10
      endif
      IF(IFTRAN) then
         LSET(ISET:ISET+8)=' UNSTEADY'
         ISET=ISET+9
      else
         LSET(ISET:ISET+6)=' STEADY'
         ISET=ISET+7
      endif
      IF(IFAXIS) then
         LSET(ISET:ISET+12)=' AXISYMMETRIC'
         ISET=ISET+13
      elseif (IF3D) then
         LSET(ISET:ISET+17)=' THREE DIMENSIONAL'
         ISET=ISET+18
      else
         LSET(ISET:ISET+15)=' TWO DIMENSIONAL'
         ISET=ISET+16
      endif
      IF(IFFLOW) then
         IF(IFNAV) then
            LSET(ISET:ISET+13)=' NAVIER-STOKES'
            ISET=ISET+14
         else
            LSET(ISET:ISET+6)=' STOKES'
            ISET=ISET+7
         endif
      endif
      IF(IFHEAT) then
         IF(IFFLOW) then
            LSET(ISET:ISET+17)=' AND HEAT TRANSFER'
            ISET=ISET+18
         else
            LSET(ISET:ISET+13)=' HEAT TRANSFER'
            ISET=ISET+14
         endif
      endif
      IF(NPSCAL.GT.0) then
         LSET(ISET:ISET+22)=' WITH   PASSIVE SCALARS'
         WRITE(LSET(ISET+6:ISET+6),'(I1)')NPSCAL
         ISET=ISET+23
      endif
      CALL PRS(' CURRENT SETTING:'//LSET(1:60)//'$')
      IF(ISET.GT.60)CALL PRS(LSET(61:120)//'$')
      ITEM(1)='ACCEPT CURRENT SWITCHES'
C      ITEM(2)='NEKTON VERSION'
      ITEM(2)='              '
      ITEM(3)='DIMENSION     '
      if (.NOT.IFTRAN) ITEM(4)='STEADY      '
      if (     IFTRAN) ITEM(4)='UNSTEADY      '
      ITEM(5)='FLUID FLOW    '
      IF(IFFLOW) then
         ITEM(6)='ADVECTION         '
         ITEM(7)='STRESS FORMULATION'
         ITEM(8)='SPLIT  FORMULATION'
         ITEM(9)='TURBULENCE MODEL  '
      else
         ITEM(6)='  '
         ITEM(7)='  '
         ITEM(8)=' '
         ITEM(9)=' '
      endif
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
         IF(I.EQ.3) then
            IF(     IF3D)ITEMON(25:25)='3'
            IF(.NOT.IF3D)ITEMON(25:25)='2'
            IF(   IFAXIS)ITEMON(25:25)='A'
c        elseif (I.EQ.4) then
c           IF(     IFTRAN)ITEMON(25:25)='Y'
c           IF(.NOT.IFTRAN)ITEMON(25:25)='N'
         elseif (I.EQ.5) then
            IF(     IFFLOW)ITEMON(25:25)='Y'
            IF(.NOT.IFFLOW)ITEMON(25:25)='N'
         elseif (I.EQ.6 .AND. IFFLOW) then
            IF(     IFNAV )ITEMON(25:25)='Y'
            IF(.NOT.IFNAV )ITEMON(25:25)='N'
         elseif (I.EQ.7 .AND. IFFLOW) then
            IF(     IFSTRS)ITEMON(25:25)='Y'
            IF(.NOT.IFSTRS)ITEMON(25:25)='N'
         elseif (I.EQ.8 .AND. IFFLOW) then
            IF(     IFSPLIT)ITEMON(25:25)='Y'
            IF(.NOT.IFSPLIT)ITEMON(25:25)='N'
         elseif (I.EQ.9 .AND. IFFLOW) then
            ITEMON(24:25)=TURBMOD(1:2)
            IF(TURBMOD(1:2).EQ.'NO')ITEMON(24:25)=' N'
         elseif (I.EQ.10) then
            IF(     IFHEAT)ITEMON(25:25)='Y'
            IF(.NOT.IFHEAT)ITEMON(25:25)='N'
         elseif (I.EQ.11) then
            WRITE(ITEMON(25:25),'(I1)')NPSCAL
C           MPSCAL=9 IN PARAMETER STMT
         elseif (I.EQ.12) then
            IF(     IFMGRID)ITEMON(25:25)='Y'
            IF(.NOT.IFMGRID)ITEMON(25:25)='N'
         endif
         ITEM(I)=ITEMON
20    CONTINUE
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOECHO')
C     Toggle logical flags based on choices
C     Do the opposite of what flag currently is, so it gets changed.
C     Round robin for dimensions
      IF(choice.eq.'DIMENSION               A')IF3D  =.TRUE.
      IF(choice.eq.'DIMENSION               A')IFAXIS=.FALSE.
      IF(choice.eq.'DIMENSION               2')IF3D  =.FALSE.
      IF(choice.eq.'DIMENSION               2')IFAXIS=.TRUE.
      IF(choice.eq.'DIMENSION               3')IF3D  =.FALSE.
      IF(choice.eq.'DIMENSION               3')IFAXIS=.FALSE.
      IF(choice.eq.'UNSTEADY                 ')IFTRAN=.FALSE.
      IF(choice.eq.'STEADY                   ')IFTRAN=.TRUE.
      IF(choice.eq.'HEAT TRANSFER           Y')IFHEAT=.FALSE.
      IF(choice.eq.'HEAT TRANSFER           N')IFHEAT=.TRUE.
      IF(choice.eq.'FLUID FLOW              Y')IFFLOW=.FALSE.
      IF(choice.eq.'FLUID FLOW              N')IFFLOW=.TRUE.
      IF(choice.eq.'ADVECTION               Y')IFNAV =.FALSE.
      IF(choice.eq.'ADVECTION               N')IFNAV =.TRUE.
      IF(choice.eq.'STRESS FORMULATION      Y')IFSTRS=.FALSE.
      IF(choice.eq.'STRESS FORMULATION      N') then
         IFSTRS = .TRUE.
         IFSPLIT= .FALSE.
      endif
      IF(choice.eq.'SPLIT  FORMULATION      Y') then
         IFSPLIT=.FALSE.
         CALL PRS('Time stepping scheme set to Second Order$')
         PARAM(27) = 2
         REQD (27) = 'P'
      endif
      IF(choice.eq.'SPLIT  FORMULATION      N') then
         IFSPLIT=.TRUE.
         IFSTRS =.FALSE.
C        Set TORDER
         CALL PRS('Time stepping scheme set to First Order$')
         PARAM(27) = 1
         REQD (27) = ' '
      endif
      IF(CHOICE(1:10).EQ.'TURBULENCE') then
C         CALL PRS
C     $   ('**** WARNING **** Turbulence Models under development.$')
C         CALL PRS
C     $   ('Use at your own risk. By hitting a carriage return,$')
C         CALL PRS('you acknowledge having read this warning.$')
C         CALL RES(LINE,0)
          CALL PRS('Turbulence Models not yet implemented$')
      endif
      IF(.FALSE.) then
       IF(choice.eq.'TURBULENCE MODEL        N')TURBMOD='KEPSRNG'
       IF(choice.eq.'TURBULENCE MODEL       KE')TURBMOD='MIXLRNG'
       IF(choice.eq.'TURBULENCE MODEL       MI')TURBMOD='NONE'
C
       IF(TURBMOD.EQ.'KEPSRNG'.OR.TURBMOD.EQ.'MIXLRNG') IFMODEL=.TRUE.
       IF(TURBMOD.EQ.'KEPSRNG'                        ) IFKEPS =.TRUE.
      endif
C
C      IF(choice.eq.'MULTIGRID SOLUTION      Y')IFMGRID=.FALSE.
C      IF(choice.eq.'MULTIGRID SOLUTION      N')IFMGRID=.TRUE.
      IF(choice.eq.'MULTIGRID SOLUTION      N') then
         CALL PRS('Multigrid not implemented in this version$')
      endif
      IF(CHOICE(1:24).EQ.'ADDL PASSIVE SCALARS    ')NPSCAL=NPSCAL+1
      IF(NPSCAL.GT.MPSCAL)NPSCAL=0
      PARAM(23)=NPSCAL
C
      IF(choice.eq.'HEAT TRANSFER           Y' .OR.
     $   choice.eq.'HEAT TRANSFER           N' .OR.
     $   choice.eq.'FLUID FLOW              Y' .OR.
     $   choice.eq.'FLUID FLOW              N') then
C        Reset the IFADVC switches to the default value for the new eq type
         DO 402 I=1,11
            IF(.NOT.IFFLOW)IFADVC(I)=.FALSE.
            IF(     IFFLOW.AND.I.NE.1)IFADVC(I)=.TRUE.
402      CONTINUE
      endif
C
      IF(CHOICE(1:23).EQ.'ACCEPT CURRENT SWITCHES') then
C         Set appropriate stuff (eqtype=; check for reasonableness)
          IF(     IFAXIS)AXIS=1
          IF(.NOT.IFAXIS)AXIS=0
          IF(     IF3D  )NDIM=3
          IF(.NOT.IF3D  )NDIM=2
          NFLDS=1
          IF(     IFHEAT)NFLDS=NFLDS+1
          NFLDS=NFLDS+NPSCAL
          IF(NKTONV.EQ.2) then
C            Nekton 2.0
C            Check if requested eqtype is available as of yet
             IF(IFNAV.and..NOT.IFTRAN) then
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
             endif
             IF(IFMGRID .AND. .NOT. IFFLOW) then
                 CALL PRS(' ****** ERROR ********$')
                 CALL PRS(' Multigrid only works with fluid flow$')
                 CALL PRS(' in which Nu is constant.$')
                 CALL BEEP
                 GO TO 1
             endif
             IF(IFFLOW.AND.IFSPLIT.and..NOT.IFTRAN) then
                 CALL PRS(' ****** ERROR ********$')
                 CALL PRS(' Split Formulation only works for unsteady$')
                 CALL PRS(' Flows$')
                 CALL BEEP
                 GO TO 1
             endif
           endif
C          Check for errors
           CALL CHKEQT(IERR)
           IF(IERR.NE.0) GO TO 1
C          Now decide whether there is convection in temperature and passive
C          scalar fields
C          No steady-state forced convection yet.
           IFNEEDC=.FALSE.
           IF(IFFLOW .AND. IFADVC(1))IFNEEDC=.TRUE.
           IF(IFHEAT.AND.IFTRAN) then
C             IFADVC(2) is for heat
              DO 50 I=0,NPSCAL
                 IFLD=I+2
49               IF(IFLD.EQ.2) then
                    CALL PRS('Convection in Temperature Field$')
                    IF(.NOT.IFFLOW) then
                     CALL PRS('Forced convection is turned on here if$')
                     CALL PRS('you want a forcing term from a steady$')
                     CALL PRS('velocity field.  You will define that$')
                     CALL PRS('field later as a fortran function or$')
                     CALL PRS('to be read in from a field file.$')
                    endif
                 endif
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
                 IF(choice.eq.'WITH CONVECTION         Y') then
                    IFADVC(IFLD)=.FALSE.
                    GO TO 49
                 elseif (choice.eq.'WITH CONVECTION         N') then
                    IFADVC(IFLD)=.TRUE.
                    GO TO 49
                 elseif (choice.eq.'FLUID       MESH') then
                    IFTMSH(IFLD)=.TRUE.
                    GO TO 49
                 elseif (choice.eq.'TEMPERATURE MESH') then
                    IFTMSH(IFLD)=.FALSE.
                    GO TO 49
                 elseif (choice.eq.'ACCEPT CURRENT SWITCHES') then
                    CONTINUE
                    IF(IFADVC(IFLD))IFNEEDC=.TRUE.
                 endif
50            CONTINUE
           endif
C          Characteristics
           IF(IFNEEDC) then
149           CONTINUE
              ITEM(1)='ACCEPT CURRENT SWITCHES'
              IF(     IFCHAR)ITEM(2)='CHARACTERISTICS         Y'
              IF(.NOT.IFCHAR)ITEM(2)='CHARACTERISTICS         N'
              NCHOIC=2
              CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOECHO')
              IF(choice.eq.'CHARACTERISTICS         Y') then
                 IFCHAR=.FALSE.
                 PARAM(26) = 0.25
                 CALL PRS
     $           ('*** Default COURANT Number reset to 0.25 ****$')
                 CALL PRS('Modify in Set Parameter Menu if desired$')
                 GO TO 149
              elseif (choice.eq.'CHARACTERISTICS         N') then
                 IFCHAR=.TRUE.
                 PARAM(26) = 1.0
                 CALL PRS
     $           ('*** Default COURANT Number set to 1.0 ****$')
                 CALL PRS('Modify in Set Parameter Menu if desired$')
                 GO TO 149
              elseif (choice.eq.'ACCEPT CURRENT SWITCHES') then
                 CONTINUE
              endif
           else
              IFCHAR=.FALSE.
           endif
C
           IF(IFNOSEG) CALL DRCOVR(13)
           return
      endif
      GO TO 1
      end
c-----------------------------------------------------------------------
      subroutine chkeqt(ierr)
      include 'basics.inc'
C
C     Checks for illegal equation settings.  IERR=1 if any were found.
      IERR=0
      IF(NPSCAL.GT.0 .AND. .NOT. IFHEAT) THEN
         CALL PRS('**ERROR** Cannot have extra passive scalars$')
         CALL PRS('without heat transfer.  Try again.$')
         CALL BEEP
         IERR=1
      endif
      IF(.NOT. IFFLOW .AND. .NOT. IFHEAT) THEN
         CALL PRS('**ERROR** Not specifying any variables makes for a'//
     $   ' poor simulation.$')
         CALL PRS('Specify heat and/or fluid.$')
         CALL BEEP
         IERR=1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine setreq
      include 'basics.inc'
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
      return
      end
c-----------------------------------------------------------------------
      subroutine setpar
C     Based on EQTYPE and  and IF3D, Prompt for Parameters to be set
      include 'basics.inc'
      CHARACTER*19 BOXLAB

      return   !  edit file to set params


      CALL GINDIS
C     4107 Tends to send <cr> when you disable gin.  This reads it.
      IF(IF4107)CALL RES(LINE,70)
      IF(IF4107) then
         IF(IFLEARN)WRITE(3,'(A70)') LINE
         CALL PRS
     $   ('SET PARAMETER menu.  TYPE in with KEYBOARD (drop mouse)$')
         CALL PRS
     $   ('Type in new value or <cr> to keep current (default) value$')
      endif
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
C       IF(IREQ.EQ.3) then
c          IF(PARAM(37).EQ.0.0) then
C            Mulitgrid Parameters Haven't been set; Make defaults
C             IF(NORDER.LT.5)CALL PRS('ERROR; you need NORDER at least'//
C     $       '5 to use mulitgrid$')
C             IF(NORDER.LT.7)NGRIDS=2
C             IF(NORDER.GE.7)NGRIDS=3
C             param(37)=NGRIDS
C             IF(NGRIDS.eq.2)then
C                PARAM(38)=NORDER/2+1
C             else if(ngrids.eq.3) then
C                PARAM(38)=NORDER/2+2
C                PARAM(39)=NORDER/2
C             endif
C          endif
C       endif
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
               IF(I.GT.NPARAM) then
C                 Setting Extra Passive scalars
                  IPSCAL=I-NPARAM
                  IF(I.GT.NPARAM+NPSCAL)IPSCAL=I-(NPARAM+NPSCAL)
                  IF(LOOP.EQ.1) THEN
                     IF(I.GT.NPARAM+NPSCAL) THEN
C                       RHOCP
                        LINE='RHOCP'
                     else
                        LINE='CONDUCT'
                     endif
                     WRITE(LINE(9:9),'(I1)')IPSCAL
                     item(iBOX)=LINE
                  elseif (LOOP.EQ.2) THEN
C                    WRITE VALUES
                     IF(I.GT.NPARAM+NPSCAL) THEN
C                       RHOCP
                        write(BOXLAB,'(G14.6)')PRHOCP(IPSCAL+2)
                     else
                        write(BOXLAB,'(G14.6)')PCOND(IPSCAL+2)
                     endif
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  endif
               else
                  IF(LOOP.EQ.1) item(iBOX)=cparam(i)
                  IF(LOOP.EQ.2) THEN
C                    WRITE VALUES
                     IF(I.LE.NPARAM) then
                        write(BOXLAB,'(G14.6)')PARAM(I)
                     elseif (I.GT.NPARAM+NPSCAL) THEN
                        write(BOXLAB,'(G14.6)')PRHOCP(IPSCAL+2)
                        CALL PRSIR('IPSCAL,PRHOCP$'
     $                  ,IPSCAL,PRHOCP(IPSCAL))
                     elseif (I.LE.NPARAM+NPSCAL) THEN
                        write(BOXLAB,'(G14.6)')PCOND (IPSCAL+2)
                        CALL PRSIR('IPSCAL,PCOND $'
     $                  ,IPSCAL,PCOND(IPSCAL+2))
                     endif
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  endif
               endif
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
                  IF(IFSUN) CALL RER(VAL)
                  IF(IFGKS) CALL RER(VAL)
                  IF(LINE.EQ.'      ') THEN
C                    Do nothing
                  else
                     IF(IF4107.OR.IFXWIN)CALL READER(VAL,IERR)
                     IF(IFSUN) IERR=0
                     IF(IERR.NE.0) then
                        CALL PRS(' The line you typed in:$')
                        CALL PRS(LINE//'$')
                        CALL PRS(' Was not understood.  Try again.$')
                        GO TO 301
                     endif
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
                     elseif (I.GT.NPARAM) THEN
                        PCOND (IPSCAL+2)=val
                     else
                        PARAM(I)=VAL
                     endif
                     write(BOXLAB,'(G14.6)')VAL
                     BOXLAB(19:19)='$'
                     CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
                  endif
C                 Check if new value is reasonable
C                 ??!! Make more extensive tests
                  IF(REQD(I).EQ.'R' .AND.PARAM(I).EQ.0.0) THEN
                     CALL PRS
     $               (' **ERROR** '//CPARAM(I)//'Must be NONZERO$')
                     CALL PRS(S//'$')
                     CALL PRS(' Try again: type in a nonzero value$')
                     CALL BEEP
                     GO TO 301
                  elseif (REQD(I).EQ.'P' .AND.PARAM(I).LE.0.0) THEN
                     CALL PRS
     $               (' **ERROR** '//CPARAM(I)//'Must be POSITIVE$')
                    CALL PRS(S//'$')
                    CALL PRS(' Try again: type in a POSITIVE value$')
                    CALL BEEP
                    GO TO 301
                  endif
                  IF(CPARAM(I).EQ.'GRID'.AND.GRID.GE.0.5) then
                     CALL PRS(' Grid should be less than 0.5$')
                     CALL PRS
     $               (' Grid > 0.5 gives only 4 possible points.$')
                     GO TO 301
                  endif
                  NORDER=PARAM(20)
                  if (NORDER.GT.NXM) THEN
                     WRITE(S,104) NX,NXM
  104                FORMAT(' Warning, current NORDER exceeds NXM'
     $                     ,' resetting from',I3,' to',I3,'.$')
                     CALL PRS(S)
                     NORDER=NXM
                     PARAM(20)=NORDER
                  endif
                  IF(CPARAM(I).EQ.'NORDER'.AND.NORDER.GE.NXM) then
                     CALL PRSI(' ERROR: Maximum NORDER set to $',NXM)
                     CALL PRS(' Please contact your NEKTONICS '//
     $               'representative to increase it.$')
                     CALL BEEP
                     GO TO 301
                  endif
C                 UN-Blink Current box
                  CALL COLOR(1)
                  CALL MOVESC(XLMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YBS(IBOX))
                  CALL DRAWSC(XRMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YTS(IBOX))
                  CALL DRAWSC(XLMEN,YBS(IBOX))
               endif
C              Check that Value is Legal
               IF(IREQ.EQ.2 .AND.VAL.EQ.0.0) THEN
                  CALL PRS(' $')
C                 Now Re-do box with new (or default) Value
               endif
            endif
5        CONTINUE
         NCHOIC=IBOX
         IF(LOOP.EQ.1)CALL DRMENU('SET PARAMETER')
         IF(LOOP.EQ.3.AND.IFNOSEG)CALL DRCOVR(13)
40    CONTINUE
      IF(IERR.NE.0) GO TO 1
      CALL GNABLE
C     Now get rid (permanently) of the stuff in segment 10 (the gwrited stuff)
      CALL CLSSEG
      CALL DelseG(10)
      CALL OPNSEG(10)
      return
      end
c-----------------------------------------------------------------------
      subroutine ckpar(mode,ierr)
      include 'basics.inc'
C     This was stolen from above setpar routine
      CHARACTER*26 BOXLAB
      CHARACTER*5 MODE
      DIMENSION ISHOW(20)
      CALL SETREQ
      IERR=0
      IEQTYPE=PARAM(16)
      IF(MODE.EQ.'SHOW ') then
         DO 4 IREQ=1,2
            IF(IREQ.EQ.1)CALL PRS(' *** REQUIRED PARAMETERS ***$')
            IF(IREQ.EQ.2)CALL PRS(' *** OPTIONAL PARAMETERS ***$')
            NSHOW=0
            DO 3 I=1,NPARAM
               IF( (REQD(I).EQ.'R'.OR.REQD(I).EQ.'P')
     $         .AND.IREQ.EQ.1) THEN
                  NSHOW=NSHOW+1
                  ISHOW(NSHOW)=I
               elseif (REQD(I).EQ.'O' .AND.IREQ.EQ.2) then
                  NSHOW=NSHOW+1
                  ISHOW(NSHOW)=I
               endif
3           CONTINUE
            IF(NSHOW.GT.0) then
               NLOOPS=((NSHOW-1)/3)+1
               DO 1 ILOOP=1,NLOOPS
                  II1=1+(ILOOP-1)*3
                  II2=II1+2
                  IF(II2.GT.NSHOW)II2=NSHOW
                  WRITE(S,'(3(1X,A10,1X,G14.6))')
     $            (CPARAM(ISHOW(II)),PARAM(ISHOW(II)),II=II1,II2)
               CALL PRS(S//'$')
 1             CONTINUE
            endif
4        CONTINUE
c     elseif (MODE.EQ.'CHECK') then
c        DO 5 I=1,NPARAM
c           IF(REQD(I).EQ.'R' .AND.PARAM(I).EQ.0.0) THEN
c              CALL PRS(' **ERROR** '//CPARAM(I)//'Must be NONZERO$')
c              CALL PRS(' Returning to beginning of parameter menu$')
c              CALL PRS(' Go back and enter a nonzero value$')
c              CALL BEEP
c              IERR=1
c           elseif (REQD(I).EQ.'P' .AND.PARAM(I).LE.0.0) THEN
c pff          CALL PRS(' **ERROR** '//CPARAM(I)//'Must be POSITIVE$')
c annoyed      CALL PRS(' Returning to beginning of parameter menu$')
c 3/12/98      CALL PRS(' Go back and enter a POSITIVE value$')
c              CALL BEEP
c              IERR=1
c           endif
5        CONTINUE
      endif
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine readln(locate)
      include 'basics.inc'
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
      return
50    IN=5
      GO TO 1
60    CALL EREXIT
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine prexit
      include 'basics.inc'
      CHARACTER CTEMP*80,CHAR1*1,CHTEMP*3
c      LOGICAL IFMVBD
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      COMMON/INOUT/  IEXT


      write(6,*) nel,' this is nel in prexit'
      call curcnt  ! Recount number of curved sides


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
      if (Ioopt.eq.0) IFFMTIN=.FALSE.
C
      M=IEXT
      n=m+3
      filenm = sesion
C
      if (.not.IFFMTIN) THEN
         FILENM(M:N) ='.re2'
         write(6,*) 'opening unformatted file ',filenm
         OPEN(UNIT=11,FILE=FILENM,STATUS='NEW',FORM='UNFORMATTED')
      endif
      NELsgn=NEL
      if (.not.IFFMTIN) NELsgn = -NEL
      WRITE(10,'(3I10,'' NEL,NDIM,NELV'')')NELsgn,NDIM,NELV
C
      DO 98 IEL=1,NEL
         if (IEL.GT.52) LETAPT(IEL) = 'A'
         if (IFFMTIN) THEN
C           formatted
            WRITE(10,'(A15,I9,A2,I5,A1,A10,I6)')
     $   '      ELEMENT  ',IEL,' [',NUMAPT(IEL),LETAPT(IEL),
     $   ']    GROUP',IGROUP(IEL)
            IF(NDIM.EQ.2) then
               WRITE(10,'(4G15.7)',ERR=60)(X(IC,IEL),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IC,IEL),IC=1,4)
            elseif (NDIM.EQ.3) then
               WRITE(10,'(4G15.7)',ERR=60)(X(IC,IEL),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IC,IEL),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(Z(IC,IEL),IC=1,4)
               WRITE(10,'(4G15.7)',ERR=60)(X(IC,IEL),IC=5,8)
               WRITE(10,'(4G15.7)',ERR=60)(Y(IC,IEL),IC=5,8)
               WRITE(10,'(4G15.7)',ERR=60)(Z(IC,IEL),IC=5,8)
            endif
         else
C           Unformatted
            WRITE(11) IGROUP(IEL)
            IF(NDIM.EQ.2) then
               WRITE(11,ERR=60)(X(IC,IEL),IC=1,4)
               WRITE(11,ERR=60)(Y(IC,IEL),IC=1,4)
            elseif (NDIM.EQ.3) then
               WRITE(11,ERR=60)(X(IC,IEL),IC=1,4)
               WRITE(11,ERR=60)(Y(IC,IEL),IC=1,4)
               WRITE(11,ERR=60)(Z(IC,IEL),IC=1,4)
               WRITE(11,ERR=60)(X(IC,IEL),IC=5,8)
               WRITE(11,ERR=60)(Y(IC,IEL),IC=5,8)
               WRITE(11,ERR=60)(Z(IC,IEL),IC=5,8)
            endif
         endif
C
98    CONTINUE
      if (IFFMTIN) THEN
         WRITE(10,*)' ***** CURVED SIDE DATA *****'
         WRITE(10,'(I10,A20,A33)')NCURVE,' Curved sides follow',
     $   ' IEDGE,IEL,CURVE(I),I=1,5, CCURVE'
         if(ncurve.gt.0)then
            do 19 iel=1,nel
            do 19 iedge=1,12
               if(ccurve(iedge,iel).ne.' ')then
                  if (nel.lt.1000) then
                     write(10,'(i3,i3,5g14.6,1x,a1)')iedge,iel,
     $               (curve(i,iedge,iel),i=1,5),ccurve(iedge,iel)
                  elseif (nel.lt.1000000) then
                     write(10,'(i2,i6,5g14.6,1x,a1)')iedge,iel,
     $               (curve(i,iedge,iel),i=1,5),ccurve(iedge,iel)
                  else
                     write(10,'(i2,i12,5g14.6,1x,a1)')iedge,iel, ! Feb, 2013
     $               (curve(i,iedge,iel),i=1,5),ccurve(iedge,iel)
                  endif
               endif
19          continue
         endif
C
      else
C
C Unformatted read
C
         WRITE(11) NCURVE
         IF(NCURVE.GT.0) then
            DO 119 IEL=1,NEL
            DO 119 IEDGE=1,12
               IF(CCURVE(IEDGE,IEL).NE.' ') then
                  WRITE(11) IEDGE,IEL,
     $            (CURVE(I,IEDGE,IEL),I=1,5),CCURVE(IEDGE,IEL)
               endif
119         CONTINUE
         endif
      endif
C
      NSIDES=4
      IF(IF3D)NSIDES=6
      if (IFFMTIN) WRITE(10,*)' ***** BOUNDARY CONDITIONS *****'
      NN=2
      IF(NFLDS.GT.2)NN=NFLDS
      DO 86 IFLD=1,NN
         IP=IFLD-2
         IF( (IFLD.EQ.1.AND.IFFLOW).OR.(IFLD.EQ.2.AND.IFHEAT).OR.
     $   IFLD.GT.2 ) then
            IF(IFLD.EQ.1 .and. IFFMTIN )WRITE(10,*)
     $      ' ***** FLUID   BOUNDARY CONDITIONS *****'
            IF(IFLD.EQ.2 .and. IFFMTIN )write(10,*)
     $      ' ***** THERMAL BOUNDARY CONDITIONS *****'
            IF(IFLD.GT.2 .and. IFFMTIN )write(10,*)
     $      ' ***** PASSIVE SCALAR',ip,' BOUNDARY CONDITIONS *****'
C           Fluid and/or thermal
C           !!?? NELF DIFFERENT FROM NEL??
            nelb = nelv
            if (ifld.eq.2) nelb = nel
            do 85 iel=1,nelb
               DO 85 ISIDE=1,NSIDES
                  CHTEMP='   '
                  IF(IFLD.EQ.1 .OR. (IFLD.EQ.2 .AND. .NOT. IFFLOW))
     $            CHTEMP = CBC(ISIDE,IEL,0)
                  if (NEL.LT.1000.and.IFFMTIN) THEN
                     WRITE(10,'(A1,A3,2I3,5G14.6)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  elseif (nel.lt.100000.and.iffmtin) then
                     WRITE(10,'(A1,A3,I5,I1,5G14.6)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  elseif (nel.lt.1 000 000.and.iffmtin) then
                     WRITE(10,'(A1,A3,I6,5G14.6)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  elseif (iffmtin) then
                     WRITE(10,'(A1,A3,I12,5G18.11)',ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  else
                     WRITE(11,ERR=60)
     $               CHTEMP,
     $               CBC(ISIDE,IEL,IFLD),IEL,ISIDE,
     $               (BC(II,ISIDE,IEL,IFLD),II=1,5)
                  endif
c                 ICBC=ICHAR(CBC(ISIDE,IEL,IFLD))
c                 IF(ICBC.GE.97 .AND. ICBC.LE.122 ) then
C                    Small letter for fortran b.c.'s
C                    Extra line(s) for Inflow.  BC(1) has # of lines; BC(2) addr
c                    CBC3=CBC(ISIDE,IEL,IFLD)
c                    IF(CBC3(3:3) .EQ. 'i') THEN
C                       Special storage of pointers in Internal boundaries
c                       NLINES=BC(4,ISIDE,IEL,IFLD)
c                       LINE1 =BC(5,ISIDE,IEL,IFLD)
c                    else
c                       NLINES=BC(1,ISIDE,IEL,IFLD)
c                       LINE1 =BC(2,ISIDE,IEL,IFLD)
c                    endif
c                    if (IFFMTIN) THEN
c                       DO 82 I=1,NLINES
c                          WRITE(10,'(A70)',ERR=60)INBC(LINE1+I-1)
c 82                    CONTINUE
c                    endif
c                 endif
85          CONTINUE
         else
C           NO B.C.'s for this field
            IF(IFLD.EQ.1 .and. IFFMTIN )WRITE(10,*)
     $      ' ***** NO FLUID   BOUNDARY CONDITIONS *****'
            IF(IFLD.EQ.2 .and. IFFMTIN )write(10,*)
     $      ' ***** NO THERMAL BOUNDARY CONDITIONS *****'
         endif
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
            IF(MATYPE(I,IF).NE.0) then
               NPACKS=NPACKS+1
            endif
182   CONTINUE
      NSKIP=NPACKS*4 +1
      WRITE(10,*)NSKIP,' Lines follow.'
      WRITE(10,*)NPACKS,' PACKETS OF DATA FOLLOW'
         DO 18 IGRP=-5,10
            DO 118 IF=IF1,NFLDS
               IF(MATYPE(IGRP,IF).NE.0) then
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
               endif
118         CONTINUE
18       CONTINUE
C
      WRITE(10,*)' ***** HISTORY AND INTEGRAL DATA ***** '
      WRITE(10,*)NHIS,'   POINTS.  Hcode, I,J,H,IEL'
C     Sort so that integrals are last
      IF(NHIS.GT.0) then
         DO 50 ITYPE=1,2
         DO 50 I=1,NHIS
            IF(LOCHIS(3,I).EQ.0) LOCHIS(3,I)=1
            IF(LOCHIS(1,I).NE.0.AND.ITYPE.EQ.1 .OR.
     $         LOCHIS(1,I).EQ.0.AND.ITYPE.EQ.2 )
     $      WRITE(10,'(1X,11A1,1X,4I5)')
     $      (HCODE(II,I),II=1,11),(LOCHIS(I2,I),I2=1,4)
50       CONTINUE
      endif
      WRITE(10,*)' ***** OUTPUT FIELD SPECIFICATION *****'
      NOUTS=6+NPSCAL
      WRITE(10,*)NOUTS,' SPECIFICATIONS FOLLOW'
      WRITE(10,*)IFXYO,'      COORDINATES'
      WRITE(10,*)IFVO, '      VELOCITY'
      WRITE(10,*)IFPO, '      PRESSURE'
      WRITE(10,*)IFTO, '      TEMPERATURE'
      WRITE(10,*)IFTGO,'      TEMPERATURE GRADIENT'
      WRITE(10,*)NPSCAL,'      PASSIVE SCALARS'
      IF(NPSCAL.GT.0) then
         DO 61 I=1,NPSCAL
            WRITE(10,'(1X,L1,1X,A5)')IFPSCO(I),PSNAME(I)
61       CONTINUE
      endif
      WRITE(10,*)' ***** OBJECT SPECIFICATION *****'
      WRITE(10,'(2X,I6,A17)')NSOBJS, ' Surface Objects '
      IF(NSOBJS .GT. 0) then
          DO 62 IOBJ=1,NSOBJS
             WRITE(10,'(4X,I4,A35,A20,A1)')NFACE(IOBJ),
     $       ' element faces for surface object "',SOBJ(IOBJ),'"'
             DO 63 IFACE=1,NFACE(IOBJ)
                WRITE(10,*)ILSURF(1,IOBJ,IFACE),ILSURF(2,IOBJ,IFACE)
63           CONTINUE
62        CONTINUE
      endif
      WRITE(10,'(2X,I6,A17)')NVOBJS, ' Volume  Objects '
      WRITE(10,'(2X,I6,A17)')NEOBJS, ' Edge    Objects '
      WRITE(10,'(2X,I6,A17)')NPOBJS, ' Point   Objects '
C
      CLOSE(UNIT=10)

      call session_exit

      return
c
60    CALL EREXIT
105   FORMAT(I6,3G15.7,' Con Ele#, cdivity,RhoCp,QVOL')
112   FORMAT(G15.7,5X,A10)
      end
c-----------------------------------------------------------------------
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
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine erexit
      CALL BEEP
      CALL PRS(' **** ERROR DURING WRITE *****$')
      CALL PRS(' Incomplete writing to files$')
      CALL PRS(' Check your personal Quota and $')
      CALL PRS(' and if the disk is full$')
      CALL DEVEX
      CALL EXITT
      end
c-----------------------------------------------------------------------
      subroutine init
      include 'basics.inc'
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

      ifmerge        = .false.
      ifconj_merge   = .false.
      if_unstructure = .false.

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
      IF(line.eq.'!') then
         CALL PRS('In DEMO Mode$')
         IFDEMO=.true.
         GO TO 5
      endif
      IF(line.eq.'$') then
         CALL PRS('In Learning Mode$')
         IFLEARN=.true.
         OPEN(UNIT=3,FILE='DEMO.COM',STATUS='NEW')
         WRITE(3,*)'$R NEKTON:PRENEK'
         WRITE(3,*)'!'
         WRITE(3,*)'DEMO'
         GO TO 5
      endif
      IF(LINE.NE.'          ') THEN
         sesion=line
         OPEN(UNIT=4,FILE='session.name',STATUS='NEW',ERR=6)
         write(4,'(a14)')sesion
      else
C        NO DEFAULT ON PRENEK
         GO TO 5
      endif
6     CONTINUE
C
      CALL CLSSEG
      CALL OPNSEG(10)
C     WRITE ON SCREEN
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

      call init_param ! Default parameters

      CWRITE = 1

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

      do i=1,nelm
         letapt(i)='a'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine output
      include 'basics.inc'
      CHARACTER KEY,STRING*5,NEWNAM*5
      CHARACTER*26 ITEMON
C! DEFAULT
C       Sets the parameters set in the field dumps
C
      IF(NEL.EQ.0) then
              CALL PRS(' Do Outputs AFTER you Build mesh.$')
              return
      endif
      CALL PRS(' Choose menu item to toggle storage on and off: $')
C     MMUST DECIDE HOW TO SET DEFAULTS BASED ON OCODES
      ITEM(1)=               'MAIN MENU'
      IF(     IFXYO) ITEM(2)='COORDINATES  YES'
      IF(.NOT.IFXYO) ITEM(2)='COORDINATES  NO '
1     NCHOIC=2
      IF(IFFLOW) then
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
      endif
      IF(IFHEAT) then
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
      endif
      NNOPS=NCHOIC
      IF(NPSCAL.GT.0) then
         DO 5 IPSCAL=1,NPSCAL
C            ITEMON='PASSIVE SCALAR  '
            ITEMON=PSNAME(IPSCAL)
            IF(      IFPSCO(IPSCAL)) ITEMON(14:16)='YES'
            IF(.NOT. IFPSCO(IPSCAL)) ITEMON(14:16)='NO '
            NCHOIC=NCHOIC+1
            ITEM(NCHOIC)=ITEMON
5        CONTINUE
      endif
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'OUTPUT')
      IF(CHOICE.NE.'MAIN MENU') THEN
C        Toggle
         IPSCAL=ICHOIC-nnops
         ITEMON=ITEM(ICHOIC)
         IF(ITEMON(14:16).EQ.'YES')then
            ITEMON(14:16)='NO '
            IF(IPSCAL .GT.0)IFPSCO(IPSCAL) = .FALSE.
         elseif (ITEMON(14:16).EQ.'NO ') then
            ITEMON(14:16)='YES'
            IF(IPSCAL .GT.0)IFPSCO(IPSCAL) = .TRUE.
            IF(IPSCAL .GT.0) then
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
            endif
         endif
         ITEM(ICHOIC)=ITEMON
         GO TO 1
      endif
C     Before returning to main menu,
      IF(NPSCAL.GT.0) then
         IPSAVE=0
         DO 101 I=1,NPSCAL
           IF(IFPSCO(I))IPSAVE=IPSAVE+1
101      CONTINUE
         IF(IPSAVE.GT.3) then
            CALL PRS(
     $      ' *** ERROR *** MAXIMUM OF 3 PASSIVE SCALARS CAN BE SAVED$')
            CALL PRS('Please turn off some of the output switches$')
            GO TO 1
         endif
      endif
      DO 10 I=2,NCHOIC
         ITEMON=ITEM(I)
         IF(ITEMON(1:13).EQ.'COORDINATES') then
            IF(NKTONV.EQ.1) then
               CALL PRS('ERROR: Coordinates always saved in version 1$')
            else
               IF(ITEMON(14:16).EQ.'YES')IFXYO=.TRUE.
               IF(ITEMON(14:16).EQ.'NO ')IFXYO=.FALSE.
            endif
         elseif (ITEMON(1:13).EQ.'VELOCITY') then
            IF(ITEMON(14:16).EQ.'YES')IFVO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFVO=.FALSE.
         elseif (ITEMON(1:13).EQ.'PRESSURE') then
            IF(ITEMON(14:16).EQ.'YES')IFPO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPO=.FALSE.
         elseif (ITEMON(1:13).EQ.'TEMPERATURE') then
            IF(ITEMON(14:16).EQ.'YES')IFTO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTO=.FALSE.
         elseif (ITEMON(1:13).EQ.'TGRADIENT') then
            IF(ITEMON(14:16).EQ.'YES')IFTGO=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTGO=.FALSE.
         endif
9     CONTINUE
10    CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine histry
      include 'basics.inc'
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
      IF(NHIS.GT.0) then
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
      endif
C
1     CONTINUE
C
      CALL PRS(' Choose menu item to toggle storage on and off: $')
C     MMUST DECIDE HOW TO SET DEFAULTS BASED ON OCODES
      ITEM(1)='MAIN MENU'
      ITEM(2)='ENTER POINT'
      NCHOIC=2
      IF(NLEVEL.GT.1) then
         nchoic = nchoic+1
         ITEM(nchoic)='UP   LEVEL'
         nchoic = nchoic+1
         ITEM(nchoic)='DOWN LEVEL'
      endif
      IF(IFFLOW) then
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
      endif
      IF(IFMVBD) then
         ITEMON='MESH LOCATION'
         IF(      IFXH) ITEMON(14:16)='YES'
         IF(.NOT. IFXH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      endif
      IF(IFHEAT) then
         ITEMON='TEMPERATURE'
         IF(      IFTH) ITEMON(14:16)='YES'
         IF(.NOT. IFTH) ITEMON(14:16)='NO '
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      endif
      IF(NPSCAL.GT.0) then
         ITEMON='PASSIVE SCALAR  '
         WRITE(ITEMON(16:16),'(I1)')IPSCH
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)=ITEMON
      endif
C
      IF(NKTONV.EQ.2 .AND. IFFLOW) then
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
      endif
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'HISTORY')
      IF(choice.eq.'UP   LEVEL' .OR.CHOICE.EQ.'DOWN LEVEL') then
C        Erase old mesh& draw new
         IF(choice.eq. 'UP   LEVEL'.AND.ILEVEL.EQ.NLEVEL) then
            CALL PRS(' ERROR: AT TOP LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         endif
         IF(choice.eq. 'DOWN LEVEL'.AND.ILEVEL.EQ.1) then
            CALL PRS(' ERROR: AT FIRST LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         endif
         IF(choice.eq.'UP   LEVEL')
     $   CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
         IF(choice.eq.'DOWN LEVEL')
     $   CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
         DO 500 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
500      CONTINUE
         IF(choice.eq. 'UP   LEVEL')ILEVEL=ILEVEL+1
         IF(choice.eq. 'DOWN LEVEL')ILEVEL=ILEVEL-1
         DO 510 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) then
            CALL DRAWEL( I)
            endif
510      CONTINUE
      else IF(CHOICE.NE.'MAIN MENU'.AND.CHOICE.NE.'ENTER POINT') THEN
C        Toggle
         ITEMON=ITEM(ICHOIC)
         IF(ITEMON(14:16).EQ.'YES')then
            ITEMON(14:16)='NO '
         elseif (ITEMON(14:16).EQ.'NO ') then
            ITEMON(14:16)='YES'
         endif
         ITEM(ICHOIC)=ITEMON
C
         DO 10 I=3,NCHOIC
          ITEMON=ITEM(I)
          IF(ITEMON(1:13).EQ.'COORDINATES') then
            IF(NKTONV.EQ.1) then
               CALL PRS('ERROR: Coordinates always saved in version 1$')
            else
               IF(ITEMON(14:16).EQ.'YES')IFXYH=.TRUE.
               IF(ITEMON(14:16).EQ.'NO ')IFXYH=.FALSE.
            endif
          elseif (ITEMON(1:13).EQ.'VELOCITY') then
            IF(ITEMON(14:16).EQ.'YES')IFVH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFVH=.FALSE.
            IFXH = .NOT. IFVH
          elseif (ITEMON(1:13).EQ.'MESH LOCATION') then
            IF(ITEMON(14:16).EQ.'YES')IFXH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFXH=.FALSE.
            IFVH = .NOT. IFXH
          elseif (ITEMON(1:13).EQ.'PRESSURE') then
            IF(ITEMON(14:16).EQ.'YES')IFPH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPH=.FALSE.
          elseif (ITEMON(1:13).EQ.'TEMPERATURE') then
            IF(ITEMON(14:16).EQ.'YES')IFTH=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFTH=.FALSE.
          elseif (ITEMON(1:13).EQ.'STILL POINT') then
            IF(ITEMON(14:16).EQ.'YES')IFPNT=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPNT=.FALSE.
          elseif (ITEMON(1:13).EQ.'TRACK PARTCL') then
            IF(ITEMON(14:16).EQ.'YES')IFPART=.TRUE.
            IF(ITEMON(14:16).EQ.'NO ')IFPART=.FALSE.
          elseif (ITEMON(1:14).EQ.'PASSIVE SCALAR') then
            IPSCH=IPSCH+1
            IF(IPSCH .GT. NPSCAL) IPSCH=0
            WRITE(ITEMON(16:16),'(I1)')IPSCH
          endif
10       CONTINUE
         GO TO 1
      elseif (  choice.eq.'MAIN MENU') THEN
         return
      elseif (choice.eq.'ENTER POINT') then
C       How do you enter 3-d points from preprocessor?
        CALL PRS(' Enter x,y coordinates with mouse:$')
        CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
        IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C          apparently trying to enter
           call prs(' Type in X and Y-coordinates:$')
           call rerr(xmouse,ymouse)
        endif
        IF(IF3D) then
           CALL PRS(' Type Z coordinate:$')
           CALL RER(ZMOUSE)
        endif
C       Now that we have everything, put it in arrays.
        NHIS=NHIS+1
C       Character flags
        DO 15 I=1,11
           HCODE(I,NHIS)=' '
15      CONTINUE
        IF(IFVH) then
          HCODE(1,NHIS)='U'
          HCODE(2,NHIS)='V'
          IF(IF3D)HCODE(3,NHIS)='W'
        endif
        IF(IFXH) then
          HCODE(1,NHIS)='X'
          HCODE(2,NHIS)='Y'
          IF(IF3D)HCODE(3,NHIS)='Z'
        endif
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
                else
                   R=SQRT((XPTS(I,J,IEL)-XMOUSE)**2+
     $                    (YPTS(I,J,IEL)-YMOUSE)**2+
     $             (ZMOUSE- ( (Z(5,IEL)+Z(1,IEL))/2+
     $          (Z(5,IEL)-Z(1,IEL))/2*ZPTS(K) ) )**2)
                endif
                IF(R.LT.RMIN) THEN
                    RMIN=R
                    IH=I
                    JH=J
                    KH=K
                    IELH=IEL
                endif
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
      endif
      GO TO 1
100   FORMAT(3I5,4X,4a1,
     $  '  I,J,IEL, and CODES of Point for History(3I5,4X,A4)')
      end
c-----------------------------------------------------------------------
      subroutine integq
      include 'basics.inc'
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
      IF(choice.eq.'MAIN MENU')return
C
      ISOBJ = ICHOIC
      DO 5 IHIS=1,NHIS
         IF(LOCHIS(1,IHIS) .EQ.ISOBJ) then
            CALL BEEP
            CALL PRS('Object '//SOBJ(ISOBJ)//' Already Chosen $')
            GO TO 1
         endif
5     CONTINUE
      CALL PRS('Adding '//SOBJ(ISOBJ)//'$')
      NHIS=NHIS+1
C     Character flags
      DO 15 I=1,11
           HCODE(I,NHIS)=' '
15    CONTINUE
      IF(IFFLOW) then
         HCODE(1,NHIS)='F'
         HCODE(2,NHIS)='F'
         IF(IF3D)HCODE(3,NHIS)='F'
      endif
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
      end
c-----------------------------------------------------------------------
      subroutine objec
      include 'basics.inc'
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
      IF(NLEVEL.GT.1) then
         NCHOIC = NCHOIC+1
         ITEM(NCHOIC)='UP   LEVEL'
         NCHOIC = NCHOIC+1
         ITEM(NCHOIC)='DOWN LEVEL'
      endif
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'OBJECT')
      IF(choice.eq.'UP   LEVEL' .OR.CHOICE.EQ.'DOWN LEVEL') then
C        Erase old mesh& draw new
         IF(choice.eq. 'UP   LEVEL'.AND.ILEVEL.EQ.NLEVEL) then
            CALL PRS(' ERROR: AT TOP LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         endif
         IF(choice.eq. 'DOWN LEVEL'.AND.ILEVEL.EQ.1) then
            CALL PRS(' ERROR: AT FIRST LEVEL ALREADY$')
            CALL BEEP
            GO TO 1
         endif
         IF(choice.eq.'UP   LEVEL')
     $   CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
         IF(choice.eq.'DOWN LEVEL')
     $   CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
         DO 500 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
500      CONTINUE
         IF(choice.eq. 'UP   LEVEL')ILEVEL=ILEVEL+1
         IF(choice.eq. 'DOWN LEVEL')ILEVEL=ILEVEL-1
         DO 510 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) then
            CALL DRAWEL( I)
            endif
510      CONTINUE
      elseif (  choice.eq.'MAIN MENU') THEN
         return
      elseif (choice.eq.'MODIFY OBJECT'
     $.OR.    choice.eq.'ADD    OBJECT') then
C        NSOBJS:        Number of objects
C        ISOBJ :        Current Object
C        NFACE(ISOBJ): Number of faces in current object
C
         IF(choice.eq.'ADD    OBJECT') then
            CALL PRS('Enter Object Name (20 characters or less):$')
            NSOBJS=NSOBJS+1
            ISOBJ =NSOBJS
            NFACE(ISOBJ) = 0
            CALL RES(SOBJ(ISOBJ),20)
         endif
         IF(choice.eq.'MODIFY OBJECT') then
            DO 11 IOBJ=1,NSOBJS
              ITEM(IOBJ)=SOBJ(IOBJ)
11          CONTINUE
            NCHOIC=NSOBJS
            CALL MENU(XMOUSE,YMOUSE,BUTTON,'MODIFY OBJECT')
            ISOBJ = ICHOIC
         endif
C
2        CONTINUE
         ITEM(1)='END OBJECT'
         ITEM(2)='ADD    FACE'
         ITEM(3)='DELETE FACE'
         NCHOIC=3
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'ADD    OBJECT')
         IF(choice.eq.'END OBJECT') then
            GO TO 1
         elseif (choice.eq.'ADD    FACE') then
            CALL GETFAC(IFACE,IEL,'ADD   ',ISOBJ)
            NFACE(ISOBJ) = NFACE(ISOBJ) + 1
            ILSURF(1,ISOBJ,NFACE(ISOBJ)) = IEL
            ILSURF(2,ISOBJ,NFACE(ISOBJ)) = IFACE
         elseif (choice.eq.'DELETE FACE') then
            CALL GETFAC(IFACE,IEL,'DELETE',ISOBJ)
            DO 12 I=1,NFACE(ISOBJ)
               IF(ILSURF(1,ISOBJ,I) .EQ. IEL   .AND.
     $            ILSURF(2,ISOBJ,I) .EQ. IFACE ) then
C                 Copy highest face into this storage location and reduce number
C                 of faces by one
                  ILSURF(1,ISOBJ,I)=ILSURF(1,ISOBJ,NFACE(ISOBJ))
                  ILSURF(2,ISOBJ,I)=ILSURF(2,ISOBJ,NFACE(ISOBJ))
                  NFACE(ISOBJ)= NFACE(ISOBJ)-1
                  GO TO 13
               endif
12          CONTINUE
            CALL PRS('Error: Selected face does not match object face$')
13          CONTINUE
         endif
         GO TO 2
      elseif (choice.eq.'DELETE OBJECT') then
         DO 18 IOBJ=1,NSOBJS
           ITEM(IOBJ)=SOBJ(IOBJ)
18       CONTINUE
         NCHOIC=NSOBJS
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='DO NOT DELETE'
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'DELETE OBJECT')
         IF(CHOICE.NE.'DO NOT DELETE') then
            ISOBJ = ICHOIC
            CALL PRS('Deleting '//SOBJ(ISOBJ)//'$')
            CALL DELOBJ(ISOBJ)
         endif
      endif
      GO TO 1
      end
c-----------------------------------------------------------------------
      subroutine delobj(isobjt)
      include 'basics.inc'
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
      return
      end
c-----------------------------------------------------------------------
      subroutine getfac(iface,iel,cflag,isobjt)
      include 'basics.inc'
      CHARACTER*6 CFLAG
      CHARACTER CHFACE
      CALL PRS(' Enter element side with mouse:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C          apparently is trying to type in
           call prs(' Type X and Y coordinates:$')
           call rerr(xmouse,ymouse)
      endif
C       Use closest side
        RMIN=1.0E6
        DO 250 IEL=1,NEL
           IF(NUMAPT(IEL).EQ.ILEVEL) then
              DO 240 I=1,NSIDES
C               WHAT'S XPTS??
                R=SQRT((SIDES(IEL,I,1)-XMOUSE)**2+
     $                 (SIDES(IEL,I,2)-YMOUSE)**2)
                IF(R.LT.RMIN) THEN
                    RMIN=R
                    IQSIDE= I
                    IQEL  = IEL
                    IF(IQSIDE.EQ.5.OR.IQSIDE.EQ.6) then
C                      CHECK if IT'S ABOVE OR BELOW CENTER
                       IF(YMOUSE.GT.SIDES(IQEL,IQSIDE,2)) then
                          IQSIDE=6
                       else
                          IQSIDE=5
                       endif
                    endif
                endif
240           CONTINUE
           endif
250     CONTINUE
        IFACE=IQSIDE
        IEL  =IQEL
C       Hilite side
        CALL COLOR(4)
        IF(IQSIDE.EQ.5.OR.IQSIDE.EQ.6) then
           IF(IQSIDE.EQ.5) then
              IED1=1
              IED2=4
           elseif (IQSIDE.EQ.6) then
              IED1=1
              IED2=4
           endif
           CALL MOVE(X(IQEL,IED1),Y(IQEL,IED1))
           DO 230 IEDGE=IED1,IED2
              CALL DRAWED(IQEL,IEDGE,1)
230        CONTINUE
        else
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
        endif
        CALL COLOR(1)
        CALL PENW(3)
      return
      end
c-----------------------------------------------------------------------
      subroutine drivef
      include 'basics.inc'
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      IDRIVF=1
      CALL PRS('          *** DRIVE FORCE MENU ***$')
      IF(IFFLOW) then
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
              IF(ANS.EQ.'B'.OR.ANS.EQ.'N') then
                CALL DRIVFC(IPSCAL)
              else
                CALL PRS('TRY AGAIN; B OR N$')
                GO TO 100
              endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine drivfc(ipscal)
      include 'basics.inc'
      COMMON/FORTRN/ IDRIVF,INITCS,IPFLAG,IFFLAG,IQFLAG
      character*80 CH
      DATA CH /' '/
C
      IF(ANS.EQ.'N')return
C         ZERO OUT DRIVING FORCES
          DRIVC(1)=' '
          DRIVC(2)=' '
          DRIVC(3)=' '
      IF(ANS.EQ.'P'.OR.ANS.EQ.'B') then
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
         endif
      endif
      IF(ANS.EQ.'F') then
          IFFLAG=1
          CALL PRS(
     $    'FLOW is total flow rate through the geometry and can be $')
          CALL PRS(
     $    'a function of TIME. It is input in the style of a fortran $')
          CALL PRS('statement.  Start in column 1.  Use one line.$')
          CALL PRS('Example:$')
          CALL PRS('FLOW = 1.0 + 1.0/(TIME + 1.0)$')
          CALL PRS('Enter statement For FLOW:$')
      endif
      IF(ANS.EQ.'Q') then
          IQFLAG=1
          IF(IPSCAL.GT.0) then
             CALL PRSI('PASSIVE SCALAR$',IPSCAL)
             CALL PRS(S//'$')
          endif
          CALL PRS(
     $    'QVOL internal heat source/volume (or area).  It can be a$')
          CALL PRS(
     $    'function of X, Y,[Z], TIME, TEMP, and PASS. It is input as$')
          CALL PRS(
     $    'a fortran statement.  Start in column 1.  Use one line.$')
          CALL PRS('Example:$')
          CALL PRS('QVOL = (1.0 -X) * Y + EXP(-TIME) + TEMP$')
          CALL PRS('Enter statement For QVOL:$')
      endif
C
c-----------------------------------------------------------------------
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
         IF(IF3D) then
            CALL PRS('FFZ:$')
            CALL RES(CH(8:80),73)
            IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
            DRIVC(3)=CH
         endif
      endif
      IF(ANS.EQ.'F')DRIVC(4)=CH
      IF(ANS.EQ.'Q')DRIVC(4+IPSCAL)=CH
      return
      end
C     End of PRENEK ****************************************
c-----------------------------------------------------------------------
      subroutine intcnd
      include 'basics.inc'
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
      if (choice.eq.INITEM(1)) return
      if (choice.eq.INITEM(2)) CALL INITDEF
      if (choice.eq.INITEM(3)) CALL INITFUN
      if (choice.eq.INITEM(4)) CALL INITPRE
      if (choice.eq.INITEM(5)) CALL INITRES
      GOTO 10
C
      end
c-----------------------------------------------------------------------
      subroutine initfun
      include 'basics.inc'
      character*80 CH
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
C
      if (NLINF.GT.0)
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
         IF(IFLD.EQ.1) then
              IF(IF3D)
     $        CALL PRS(' UX one line for UY and one for UZ.  Example:$')
              IF(.NOT.IF3D)
     $        CALL PRS(' UX and one for UY.  Example:$')
              CALL PRS(' UX = 0.667 * (1.0-Y**2)$')
              CALL PRS(' UY = 0.0$')
              IF(IF3D)CALL PRS('UZ = 0.0$')
              CALL PRS('  Note: Check FORTRAN syntax carefully.$')
         elseif (IFLD.EQ.2) then
              CALL PRS(' TEMP.  Example:$')
              CALL PRS(' TEMP= 0.667 * (1.0-Y**2)$')
         else
              IDUMMY=IFLD-2
              CALL PRiS(IDUMMY,'th PS.  Example:$')
              CALL PRS(' PS = 0.667 * (1.0-Y**2)$')
         endif
C
         IF(IFLD.EQ.1) then
              CALL PRS(' Enter statement for UX:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(1)=CH
              CALL PRS(' Now enter statement for UY:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(2)=CH
              IF(NDIM.EQ.3) then
                 CALL PRS(' Now enter statement for UZ:$')
                 CALL RES(CH(8:80),73)
                 IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
                 INITC(3)=CH
              endif
         elseif (IFLD.EQ.2) then
              CALL PRS(' Enter statement for TEMP:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(IFLD+3)=CH
         else
              CALL PRS(' Enter statement for PS:$')
              CALL RES(CH(8:80),73)
              IF(IFLEARN)WRITE(3,'(A73)')CH(8:80)
              INITC(IFLD+3)=CH
         endif
 1000 CONTINUE
      NLINF=NFLDS
      return
      end
c-----------------------------------------------------------------------
      subroutine initpre
      include 'basics.inc'
      CHARACTER*1 OPTS(80),INITP1(80)
      EQUIVALENCE (INITP1,INITP)
C
      NLIN=NLINP+NLINR
      if (NLIN.GT.0) CALL PRS(' Current restart/presolve status:$')
      DO 100 I=1,NLIN
         CALL PUTS(INITP(I),80)
  100 CONTINUE
C
      if (NLINP.EQ.0) THEN
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
         if (LTRUNC(OPTS,80).EQ.0) return
         if (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITPRE',1)
            GOTO 300
         else
            CALL CHCOPY(INITP1(11),OPTS,70)
         endif
      endif
C
  400 CALL PRS(' Current presolve status:$')
      CALL PUTS(INITP1,80)
  500 CALL PRS
     $(' Enter <cr> to accept, M to modify, R to remove, or ?.$')
      CALL RES(OPTS,80)
      CALL LJUST(OPTS)
      CALL CAPIT(OPTS,80)
C
      if (LTRUNC(OPTS,80).EQ.0) return
C
      if (OPTS(1).EQ.'R') THEN
         CALL PRS('Removing presolve option.$')
         NLINP=0
         DO 550 IC=1,14
            I1=IC+1
            CALL CHCOPY(INITP(IC),INITP(I1),80)
  550    CONTINUE
         CALL BLANK(INITP(15),80)
         return
      endif
C
      if (OPTS(1).EQ.'M') THEN
  600    CALL PRS(' Enter field specification or ?.$')
         CALL RES(OPTS,80)
         if (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITPRE',1)
            GOTO 600
         else
            CALL CHCOPY(INITP1(11),OPTS,70)
         endif
         GOTO 400
      endif
C
      if (OPTS(1).EQ.'?') THEN
         CALL HELPER('INITPRE',2)
         GOTO 400
      endif
C
C     Else, invalid entry, re-prompt:
      GOTO 500
C
      end
c-----------------------------------------------------------------------
      subroutine initres
      include 'basics.inc'
      CHARACTER*80 FNAME
      CHARACTER*1  OPTS(80),INITP1(80,15)
      EQUIVALENCE (INITP1,INITP)
      LOGICAL IFADD,IFREMV
C
      IFREMV=.FALSE.
   10 CONTINUE
      NLIN=NLINP+NLINR
      if (NLIN.GT.0) CALL PRS(' Current restart/presolve status:$')
      DO 100 I=1,NLIN
         CALL PUTS(INITP(I),80)
  100 CONTINUE
C
C     Add/review?
C
      IFADD=.FALSE.
      if (NLINR.GT.0.OR.IFREMV) THEN
         CALL PRS(' Input A to add a file, M to modify/remove, $')
         CALL PRS(' or <cr> to return.$')
         CALL RES(OPTS,80)
         CALL PRS(' $')
         CALL LJUST(OPTS)
         CALL CAPIT(OPTS,80)
C
         if (LTRUNC(OPTS,80).EQ.0) return
         if (OPTS(1).EQ.'A') IFADD=.TRUE.
      endif
C
C     Add a restart file
C
      if (IFADD.OR.NLINR.EQ.0) THEN
         CALL PRS(' Enter restart session or file name:$')
         CALL RES(FNAME,80)
         if (LTRUNC(FNAME,80).EQ.0) return
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
         if (LTRUNC(OPTS,80).EQ.0) return
C
         if (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITRES',1)
            GOTO 300
         else
C           specifying restart options
            IOPT=LTRUNC(INITP(NLIN),80)+2
            LOPT=80-IOPT+1
            CALL CHCOPY(INITP1(IOPT,NLIN),OPTS,LOPT)
            return
         endif
      endif
C
C     Alter/Remove restart requests.
C
      ILIN=1+NLINP
      if (NLINR.GT.1) THEN
         CALL PRS(' Enter restart file # to be altered:$')
         CALL REI(ILIN)
      endif
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
      if (LTRUNC(OPTS,80).EQ.0) return
C
      if (OPTS(1).EQ.'?') THEN
         CALL HELPER('INITRES',2)
         GOTO 500
      endif
C
      if (OPTS(1).EQ.'M') THEN
  600    CALL PRS(' Enter field specification or ?.$')
         CALL RES(OPTS,80)
         if (OPTS(1).EQ.'?') THEN
            CALL HELPER('INITRES',1)
            GOTO 600
         else
            CALL LJUST(INITP(ILIN))
            IOPT=INDX1(INITP(ILIN),' ',1)+1
            LOPT=80-IOPT+1
            CALL CHCOPY(INITP1(IOPT,ILIN),OPTS,LOPT)
         endif
         GOTO 400
      endif
C
      if (OPTS(1).EQ.'R') THEN
         DO 700 I=ILIN,14
            I1=I+1
            CALL CHCOPY(INITP(I),INITP(I1),80)
  700    CONTINUE
         CALL BLANK(INITP(15),80)
         NLINR=NLINR-1
c        kludge city
         IFREMV=.TRUE.
         GOTO 10
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine helper(routine,ientry)
      include 'basics.inc'
      CHARACTER ROUTINE*(*)
C
      if (ROUTINE.EQ.'INITPRE'.AND.IENTRY.EQ.1) THEN
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn.$')
        CALL PRS
     $  (' If presolve is requested without specification, all$')
        CALL PRS
     $  (' fields will start with a steady state solution as an$')
        CALL PRS
     $  (' initial condition.$')
        return
      endif
C
      if (ROUTINE.EQ.'INITPRE'.AND.IENTRY.EQ.2) THEN
        CALL PRS
     $  (' M - allows you to modify the field specifications.$')
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn.$')
        CALL PRS(' $')
        CALL PRS
     $  (' R - remove all presolve requests.$')
        return
      endif
C
      if (ROUTINE.EQ.'INITRES'.AND.IENTRY.EQ.1) THEN
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
        return
      endif
C
      if (ROUTINE.EQ.'INITRES'.AND.IENTRY.EQ.2) THEN
        CALL PRS
     $  (' M - allows you to modify the field specifications.$')
        CALL PRS
     $  (' Possible field specifications are U T P1 P2 ... Pn, d#,$')
        CALL PRS
     $  (' and TIME=ttt.tt (d# is the dump number, ttt.ttt is time)$')
        CALL PRS(' $')
        CALL PRS
     $  (' R - remove this restart request.$')
        return
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine initdef
      include 'basics.inc'
      CHARACTER*1 OPTS(80),INITP1(80)
      EQUIVALENCE (INITP1,INITP)
C
      CALL PRS(' Set all initial conditions to zero? (y/n)$')
      CALL RES(OPTS,80)
      CALL LJUST(OPTS)
      CALL CAPIT(OPTS,80)
      if (OPTS(1).NE.'Y') return
C
      NLINR=0
      NLINP=0
      CALL BLANK(INITP,1200)
      CALL BLANK(INITC,1200)
      DO 10 I=1,15
         INITC(I)='C Default'
   10 CONTINUE

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt
      include 'basics.inc'
      write(6,*) 'stopping in exitt'
      nz = 1/(nx-ny)
      zzz = -nel
      zzz = sqrt(zzz)
      write(6,*) zzz
      
      stop
      end
c-----------------------------------------------------------------------
      subroutine exit
      write(6,*) 'stopping in exit'
      stop
      end
c-----------------------------------------------------------------------
      subroutine init_param ! Initialize parameters
      include 'basics.inc'

      call init_cparam ! Initialize parameter names
      call init_vparam ! Initialize parameter values

      return
      end
c-----------------------------------------------------------------------
      subroutine init_cparam ! Initialize parameters
      include 'basics.inc'

      do i=1,500
         call blank(cparam(i),40)
      enddo

      cparam( 1)='DENSITY'
      cparam( 2)='VISCOS'
c     cparam( 3)='BETAG'
c     cparam( 4)='GTHETA'     ! Unused params commented out, 11/19/2010
c     cparam( 5)='PGRADX'
c     cparam( 6)='FLOWRATE'
      cparam( 7)='RHOCP'
      cparam( 8)='CONDUCT'
c     cparam( 9)='QVOL'
      cparam(10)='FINTIME'
      cparam(11)='NSTEPS'
      cparam(12)='DT'
      cparam(13)='IOCOMM'
      cparam(14)='IOTIME'
      cparam(15)='IOSTEP'
      cparam(16)='PSSOLVER: 0=default'
c     cparam(17)='AXIS'
      cparam(18)='GRID < 0 --> # cells on screen'
      cparam(19)='INTYPE'
      cparam(20)='NORDER'
      cparam(21)='DIVERGENCE'
      cparam(22)='HELMHOLTZ'
      cparam(23)='NPSCAL'

      cparam(24)='TOLREL'
      cparam(25)='TOLABS'
      cparam(26)='COURANT/NTAU'
      cparam(27)='TORDER'
      cparam(28)='TORDER: mesh velocity (0: p28=p27)'

      cparam(29)='= magnetic visc if > 0, = -1/Rm if < 0'
      cparam(30)='> 0 ==> properties set in uservp()'
      cparam(31)='NPERT: #perturbation modes'
      cparam(32)='#BCs in re2 file, if > 0'


      cparam(41)='1-->multiplicative SEMG'
      cparam(42)='0=gmres/1=pcg'
      cparam(43)='0=semg/1=schwarz'
      cparam(44)='0=E-based/1=A-based prec.'
      cparam(45)='Relaxation factor for DTFS'
      cparam(46)='reserved'
      cparam(47)='vnu: mesh matieral prop.'
c     cparam(49)='reserved'
      cparam(52)='IOHIS'

      cparam(54)='|p54|=1,2,3-->fixed flow rate dir=x,y,z'
      cparam(54)='fixed flow rate dir: |p54|=1,2,3=x,y,z'
      cparam(55)='vol.flow rate (p54>0) or Ubar (p54<0)'

      cparam(59)='!=0 --> full Jac. eval. for each el.'
      cparam(60)='!=0 --> init. velocity to small nonzero'

      cparam(62)='>0 --> force byte_swap for output'
      cparam(63)='=8 --> force 8-byte output'
      cparam(64)='=1 --> perturbation restart'

      cparam(65)='#iofiles (eg, 0 or 64); <0 --> sep. dirs'
      cparam(66)='output : <0=ascii, else binary'
      cparam(67)='restart: <0=ascii, else binary'
      cparam(68)='iastep: freq for avg_all'
      cparam(68)='iastep: freq for avg_all (0=iostep)'

      cparam(74)='verbose Helmholtz'


      cparam(84)='!=0 --> sets initial timestep if p12>0'
      cparam(85)='dt ratio if p84 !=0, for timesteps>0'
      cparam(86)='reserved' !=0 --> use skew-symm form (convection)'

      cparam(93)='Number of previous pressure solns saved'
      cparam(94)='start projecting velocity after p94 step'
      cparam(95)='start projecting pressure after p95 step'

      cparam(99)='dealiasing: <0--> off/3--> old/4--> new'
      cparam(101)='Number of additional modes to filter'
      cparam(102)='Dump out divergence at each time step'
      cparam(103)='weight of stabilizing filter (.01)'

      cparam(107)='!=0 --> add to h2 array in hlmhotz eqn'

      cparam(116)='!=0: x elements for fast tensor product'
      cparam(117)='!=0: y elements for fast tensor product'
      cparam(118)='!=0: z elements for fast tensor product'

c
c     Check these following ones
c
c     125,126,128 for perturbation code
c
c
c     150 - maxiter
c     151 - maxiter
c
c
c

      nparam=118

      return
      end
c-----------------------------------------------------------------------
      subroutine init_vparam ! Initialize parameters
      include 'basics.inc'

      call rzero(param,500)

      param( 1)= 1.0  ! 'DENSITY'
      param( 2)= .01  ! 'VISCOS'
      param( 3)=  0   ! 'BETAG'
c     param( 4)=  0   ! 'GTHETA'     ! Unused params commented out, 11/19/2010
c     param( 5)=  0   ! 'PGRADX'
c     param( 6)=  0   ! 'FLOWRATE'
      param( 7)=  1   ! 'RHOCP'
      param( 8)=  .01 ! 'CONDUCT'
c     param( 9)=  0   ! 'QVOL'
      param(10)=  0   ! 'FINTIME'
      param(11)=  1   ! 'NSTEPS'
      param(12)=  .01 ! 'DT'
      param(13)=  0   ! 'IOCOMM'
      param(14)=  0   ! 'IOTIME'
      param(15)=  100 ! 'IOSTEP'
      param(16)=  0   ! 'PSSOLVER: 0=default'
c     param(17)=  0   ! 'AXIS'
      param(18)=  -20 ! 'GRID'
      param(19)=  -1  ! 'INTYPE'
      param(20)=  4   ! 'NORDER'
      param(21)=  1.e-6 ! 'DIVERGENCE'
      param(22)=  1.e-8 ! 'HELMHOLTZ'
      param(23)=  0     ! 'NPSCAL'

      param(24)=  .01   ! 'TOLREL'
      param(25)=  .01   ! 'TOLABS'
      param(26)=  0.5   ! 'COURANT/NTAU'
      param(27)=  2     ! 'TORDER'
      param(28)=  1     ! 'TORDER: mesh velocity (0: p28=p27)'

      param(29)=  0     ! '= magnetic visc if > 0, = -1/Rm if < 0'
      param(30)=  0     ! '> 0 ==> properties set in uservp()'
      param(31)=  0     ! 'NPERT: #perturbation modes'
      param(32)=  0     ! '#BCs in re2 file, if > 0'


      param(41)=  0     ! '1-->multiplicative SEMG'
      param(42)=  0     ! '0=gmres/1=pcg'
      param(43)=  0     ! '0=semg/1=schwarz'
      param(44)=  0     ! '0=E-based/1=A-based prec.'
      param(45)=  0     ! 'Relaxation factor for DTFS'
      param(46)=  0     ! 'reserved'
      param(47)=  0     ! 'vnu: mesh matieral prop.'
      param(49)=  0     ! 'reserved'
      param(52)=  0     ! 'IOHIS'

      param(54)=  0     ! '|p54|=1,2,3-->fixed flow rate dir=x,y,z'
      param(54)=  0     ! 'fixed flow rate dir: |p54|=1,2,3=x,y,z'
      param(55)=  0     ! 'vol.flow rate (p54>0) or Ubar (p54<0)'

      param(59)=  0     ! '!=0 --> full Jac. eval. for each el.'
      param(60)=  0     ! '!=0 --> init. velocity to small nonzero'

      param(62)=  0     ! '>0 --> force byte_swap for output'
      param(63)=  0     ! '=8 --> force 8-byte output'
      param(64)=  0     ! '=1 --> perturbation restart'

      param(65)=  0     ! '#iofiles (eg, 0 or 64); <0 --> sep. dirs'
      param(66)=  0     ! 'output : <0=ascii, else binary'
      param(67)=  0     ! 'restart: <0=ascii, else binary'
      param(68)=  0     ! 'iastep: freq for avg_all'
      param(68)=  0     ! 'iastep: freq for avg_all (0=iostep)'

      param(74)=  0     ! 'verbose Helmholtz'


      param(84)=  0     ! '!=0 --> sets initial timestep if p12>0'
      param(85)=  0     ! 'updates dt if p84 !=0, for timesteps>0'
      param(86)=  0     ! '!=0 --> use skew-symm form (convection)'

c     param(87)=  0     ! '=1 --> use old convection operator'
c     param(88)=  0     ! '=1 --> skips pressure precond solver'
c     param(89)=  0     ! '!=0 --> set alpha for pressure precond'
c     param(91)=  0     ! '!=0 --> set tol for preconditioner'

      param(93)=  20    ! 'Number of previous pressure solns saved'
      param(94)=  0     ! 'start projecting velocity after p94 step'
      param(95)=  5     ! 'start projecting pressure after p95 step'

      param(99) =  3     ! 'dealiasing: <0--> off/3--> old/4--> new'
      param(100)=  0     ! 'CG pressure precond:0=Jacobi/>0=Schwarz'
      param(101)=  0     ! 'Number of additional modes to filter'
      param(102)=  1     ! 'Dump out divergence at each time step'
      param(103)=  .01   ! 'weight of stabilizing filter'
      param(104)=  0     ! 'debugging option in Jacobi precond'
      param(105)=  0     ! '=2 --> an extra orthogonalization iter'
      param(106)=  0     ! '!=0 --> dump Unfiltered-Filtered vars'
      param(107)=  0     ! '!=0 --> add to h2 array in hlmhotz eqn'

      param(116)=  0     ! '!=0: x elements for fast tensor product'
      param(117)=  0     ! '!=0: y elements for fast tensor product'
      param(118)=  0     ! '!=0: z elements for fast tensor product'

      nparam=118

      call rzero(pcond ,11)
      call rzero(prhocp,11)

      return
      end
c-----------------------------------------------------------------------
