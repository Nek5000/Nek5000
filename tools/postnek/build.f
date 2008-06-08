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
      SUBROUTINE BUILD
C
C     Menu-based module that prompts the user to input corners.
C     It then establishes connectivity, sets global node numbers, and
C     constructs elemental mesh. Calls routines that set boundary
C     conditions.
C234567890123456789012345678901234567890123456789012345678901234567890123
C
      INCLUDE 'basics.inc'
      DIMENSION ICRVS(4)
C     ELEMENT,CORNER,COOORDINATE#,LEVEL (OLD OR CURRENT)
      CHARACTER KEY,STRING*6,LETO,CHAR1
      LOGICAL IFTMP
      COMMON /SPLITT/ ENEW(NELM),IND(NELM)
C
      PI=4.0*ATAN(1.0)
      NX=PARAM(20)
      IF (NX.GT.NXM) THEN
         WRITE(S,104) NX,NXM
  104    FORMAT(' Warning, current N exceeds NXM, resetting from',I3
     $         ,' to',I3,'.$')
         CALL PRS(S)
         NX=NXM
      ENDIF
      NY=NX
      NZ=1
      IF(IF3D)NZ=NX
      IF(NKTONV.EQ.1) THEN
         DO 10 I=1,NX
            XXPTS(I) = -COS(PI*FLOAT(I-1)/FLOAT(NX-1))
            YYPTS(I) = -COS(PI*FLOAT(I-1)/FLOAT(NX-1))
   10    CONTINUE
      ELSE
C        Need legendre points for mesh drawn on screen
         CALL LEGEND(XXPTS,WGHT,NX)
         CALL LEGEND(YYPTS,WGHT,NY)
      ENDIF
      NSIDES=4
      IFCEN=.FALSE.
      IF(IF3D)NSIDES=6
      IF(CHOICE.EQ.'BUILD FROM FILE'.OR.
     $   CHOICE.EQ.'IMPORT UNIVERSAL FILE')THEN
         IF(CHOICE.EQ.'BUILD FROM FILE')CALL READAT
C        Set grid stuff based on old data
         GRIDDX=XFAC*GRID
         GRIDDY=YFAC*GRID
C        Set clipping window based on scale factors
         WT = YPHY(.995)
         WB = YPHY(.005)
         WL = XPHY(.005)
         WR = XPHY(.995)
C        Save scale factors for UnZoom function
         XFACO=XFAC
         YFACO=YFAC
         XZEROO=XZERO
         YZEROO=YZERO
C        Now set stuff based on data from readat
         DO 20 IEL=1,NEL
            HEIGHT(NUMAPT(IEL))=Z(IEL,5)-Z(IEL,1)
            XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.0
            YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.0
   20    CONTINUE
C        Display Mesh (or at least 1st floor)
         CALL DRGRID
         NLEVEL=0
C        Sort out the elements so that they will be drawn correctly
         CALL SORTEL
         DO 30 IEL=1,NEL
            CALL DRAWIS(ISRT(IEL))
            IF(NUMAPT(IEL).EQ.1)CALL DRAWEL(IEL)
            IF(NUMAPT(IEL).GT.NLEVEL)NLEVEL=NUMAPT(IEL)
   30    CONTINUE
         IF(NDIM.EQ.2)NLEVEL=1
C        Display Elevator 1st floor hilighted
         IF(NLEVEL.GT.1)THEN
            DO 40 I=NLEVEL,2,-1
               CALL DRELEV(I-1,I,'     ')
   40       CONTINUE
         ENDIF
         CALL SCROLL(5)
         CALL GNABLE
         NCHOIC =  2
         ITEM(1)='ACCEPT MESH'
         ITEM(2)='REVIEW/MODIFY'
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'ACCEPT/REVIEW')
         IF(IFNOSEG)CALL DRCOVR(13)
         IF(CHOICE.EQ.'ACCEPT MESH')THEN
            GOTO 330
         ELSE IF(CHOICE.EQ.'REVIEW/MODIFY')THEN
C           PREPARE TO MODIFY FLOOR BY FLOOR
         ENDIF
      ELSE
C        Interactive Input
         CALL SCROLL(5)
         CALL GNABLE
         CALL SETSCL
      ENDIF
C     Just in case it didn't get set in setscl
      IFGRID=.TRUE.
      IF(NLEVEL.EQ.0)NLEVEL=1
      CALL PRS('                    *** BUILD MENU ***$')
C
C
      IFCEN=.FALSE.
      IF(.NOT.IFREAD)CALL DRGRID
      IF(.NOT.IFREAD)NEL=0
C     Start at First level
      ILEVEL=1
      IF(IF3D) THEN
         IF(IFREAD)THEN
C           Put in bogus height of first level
         ELSE
C           For interactive sesion, Demand Height of First Level
            CALL PRSI('Please Enter Height of Level$',ILEVEL)
   50       CALL KEYPAD(HEIGHT(ILEVEL))
            IF(HEIGHT(ILEVEL).LE.0.0)THEN
               CALL PRS('HEIGHT must be a positive number!$')
               GOTO 50
               ENDIF
         ENDIF
C        2nd arg: 0 for no fade; negative for bottom only (for floor modify)
         CALL DRELEV(ILEVEL,0,'     ')
      ENDIF
C     Draw FIRST Level
      IF(NEL.GT.0) THEN
         DO 60 IEL=1,NEL
C            CALL DRAWIS(ISRT(IEL))
C            IF(NUMAPT(IEL).EQ.ILEVEL)CALL DRAWEL(IEL)
            IF (NUMAPT(IEL).EQ.ILEVEL)
     $         MAXLET=MAX(MAXLET,ICHAR(LETAPT(IEL)))
   60    CONTINUE
      ENDIF
      ILETAP=MAXLET+96-32
C     ! ??!!
C
 1000 CONTINUE
C***  BIG DECISION POINT  ****
C     Draw menu
      nchoic = 0
      IF(IF3D)THEN
         nchoic = nchoic+1
         ITEM(nchoic)        =             'END    ELEMENTS'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'MODIFY ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'GLOBAL REFINE'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'UP   LEVEL'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'DOWN LEVEL'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'CURVE SIDES'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'DELETE ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'CEILING'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'REDRAW ISOMETRIC'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'REDRAW MESH'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'SET GRID'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'ZOOM'
c        nchoic = nchoic+1
c        ITEM(nchoic)        =             'DEFINE OBJECT'
         IF(IFCEIL)THEN
C           Restrict choices for person on CEILING
            ITEM(1)='OFF CEILING'
            ITEM(2)='MODIFY ELEMENT'
            ITEM(3)='CURVE SIDES'
            NCHOIC=3
         ENDIF
      ELSE
C        2-D
         nchoic = nchoic+1
         ITEM(nchoic)       =             'END    ELEMENTS'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'MODIFY ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'GLOBAL REFINE'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'CURVE SIDES'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'DELETE ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'ZOOM'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'SET GRID'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'DEFINE OBJECT'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'REDRAW MESH'
      ENDIF
C
C     Menu's all set, prompt for user input:
C
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
C
      IF(CHOICE.EQ.'ADD    ELEMENT') THEN
         IF(NEL.GE.NELM-3)THEN
            CALL PRS('TOO MANY ELEMENTS$')
            CALL PRS('SEE YOUR NEKTONICS REPRESENTATIVE TO ADJUST$')
            GOTO 1000
         ENDIF
         NEL=NEL+1
         NUMAPT(NEL)=ILEVEL
         IF(ILETAP.LE.122)LETAPT(NEL)=CHAR(ILETAP)
         IF(ILETAP.GT.122)LETAPT(NEL)=' '
         CALL PRS('Enter element corners. Use mouse or$')
         CALL PRS('use keypad for x then y coordinate$')
C        Turn on Keypad
         CALL NEWEL(NEL,XMOUSE,YMOUSE,BUTTON,IERR)
         IF (IERR.EQ.1) THEN
            CALL DRAWEL(-NEL)
            NEL=NEL-1
         ENDIF
      ELSE IF(CHOICE.EQ.'CEILING') THEN
C        Make sure everythingh is drawn on ceiling.  And that
C        MODEL and CURVE know about it, too
         DO 80 I=1,NEL
             IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
   80    CONTINUE
         IFCEIL=.TRUE.
         DO 90 I=1,NEL
             IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
   90    CONTINUE
         CALL DRELEV(-(ILEVEL+1),ILEVEL,'     ')
         GOTO 1000
      ELSE IF(CHOICE.EQ.'OFF CEILING') THEN
C        Make sure everythingh is drawn on ceiling.  And that
C        MODEL and CURVE know about it, too
         DO 100 I=1,NEL
             IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
  100    CONTINUE
         IFCEIL=.FALSE.
         DO 110 I=1,NEL
             IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
  110    CONTINUE
         CALL DRELEV(ILEVEL,-(ILEVEL+1),'     ')
         GOTO 1000
      ELSE IF(CHOICE.EQ.'CURVE SIDES')THEN
         CALL CURVES
         GOTO 1000
      ELSE IF(CHOICE.EQ.'MODIFY ELEMENT')THEN
C        Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
            GOTO 1000
         ELSE
            CALL MODEL(NEL)
         ENDIF
      ELSE IF(CHOICE.EQ.'GLOBAL REFINE')THEN
C        Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
            GOTO 1000
         ELSE
C           Only floor of elevator hilighted during modify
            CALL DRELEV(-ILEVEL,ILEVEL,'     ')
            CALL GLOMOD
            CALL DRELEV(ILEVEL,0,'     ')
         ENDIF
      ELSE IF(CHOICE.EQ.'DELETE ELEMENT')THEN
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to delete$')
            GOTO 1000
         ELSE
C           Find out which element to delete
            CALL PRS(
     $      'Enter (with mouse) 2 points in element to be deleted,$')
            CALL PRS(
     $      'or, 2 points framing a box containing elements.$')
            CALL PRS(
     $      'Enter in menu area to abort DELETE ELEMENT operation.$')
            IFTMP =IFGRID
            IFGRID=.FALSE.
  120       CONTINUE
            CALL PRS('Enter 1st point:$')
            CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
            IF (XMOUSE.GT.XPHY(1.0)) THEN
C              look for a keypad input
               CALL PRS(
     $         'Select new menu item, or DELETE ELEMENT to try again.$')
               CALL BEEP
               IFGRID=IFTMP 
               GOTO 1000
            ELSE
               CALL PRS('Enter 2nd point:$')
               CALL MOUSE(XMOUS2,YMOUS2,BUTTON)
               IF (XMOUS2.GT.XPHY(1.0)) THEN
                CALL PRS(
     $         'Select new menu item, or DELETE ELEMENT to try again.$')
                CALL BEEP
                IFGRID=IFTMP 
                GOTO 1000
               ENDIF
            ENDIF
C
C           We successfully input 2 points in the build area
C              1) count number of centroids in the bounding box
C              2) if ncntrd=0, find nearest element and anihilate it.
C
C           Box
            XMAX=MAX(XMOUSE,XMOUS2)
            XMIN=MIN(XMOUSE,XMOUS2)
            YMAX=MAX(YMOUSE,YMOUS2)
            YMIN=MIN(YMOUSE,YMOUS2)
C
C           Check box to see if it contains any elements, and delete them.
C
C           This is a gross N^2 algorithm, but it beats entering them all
C           by hand....
C
            NUMDEL=0
  124       CONTINUE
            NELT=NEL
            DO 125 JEL=1,NELT
               IF (JEL.LE.NEL        .AND.
     $             XMIN.LE.XCEN(JEL) .AND.
     $             XCEN(JEL).LE.XMAX .AND.
     $             YMIN.LE.YCEN(JEL) .AND.
     $             YCEN(JEL).LE.YMAX ) THEN
                      IF (NUMDEL.EQ.0) 
     $                CALL DRWBOX(XMIN,YMIN,XMAX,YMAX,1)
                      CALL DELEL(JEL)
                      NUMDEL=NUMDEL+1
                      GOTO 124
               ENDIF
  125       CONTINUE
            IF (NUMDEL.GT.0) THEN
C              redraw the mesh
               CALL REFRESH
               CALL DRMENU('NOCOVER')
               CALL DRGRID
               DO 128 IEL=1,NEL
                  CALL DRAWEL(IEL)
  128          CONTINUE
               IFGRID=IFTMP 
               GOTO 1000
            ELSE
C              Look for closest element (standard delete element option)
               XMOUSE=(XMIN+XMAX)/2.0
               YMOUSE=(YMIN+YMAX)/2.0
               RMIN=1.0E20
               DO 130 JEL=1,NEL
                  IF (.NOT.IF3D .OR. NUMAPT(JEL).EQ.ILEVEL) THEN
                     R2= (XCEN(JEL)-XMOUSE)**2 + (YCEN(JEL)-YMOUSE)**2
                     IF(R2.LT.RMIN) THEN
                        RMIN  = R2
                        IELMIN= JEL
                     ENDIF
                  ENDIF
  130          CONTINUE
            ENDIF
C           Check if it$')s OK to delete
            IF(IF3D)THEN
               IEQ=0
               IGT=0
               DO 140 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL  )IEQ=IEQ+1
                  IF(NUMAPT(I).EQ.ILEVEL+1)IGT=IGT+1
  140          CONTINUE
               IF(IEQ.EQ.1)THEN
C                 Last element on floor
                  IF(IGT.GT.0)THEN
                     CALL PRS('**ERROR** You cannot delete all the$')
                     CALL PRS('elements on a level when there are $')
                     CALL PRS('elements above$')
                     IFGRID=IFTMP 
                     GOTO 1000
                  ELSE IF(ILEVEL.GT.1)THEN
                     CALL PRS(' ** EMPTY LEVEL ** Please ADD ELEMENT$')
                     CALL PRS('Go DOWN LEVEL to continue$')
                     nlevel=nlevel-1
                  ENDIF
               ENDIF
            ENDIF
            CALL PRSI('Deleting Element $',IELMIN)
C           Now, black out old element
            CALL DRAWEL(-IELMIN)
            CALL DRAWIS(-IELMIN)
C           Now draw any elements with higher priority
C           ISRT contains the numbers of the elements in order of
C           plotting priority.  Element (ISRT(NEL)) is most visible.
            CALL DELEL(IELMIN)
C           Sort out the elements so that they will be drawn correctly
            CALL SORTEL
C           Now redraw all the isometric elements.
            DO 150 I=1,NEL
C               CALL DRAWIS(ISRT(I))
  150       CONTINUE
            IF(.NOT.IF3D .AND.IELMIN.LE.NEL)CALL DRAWEL(IELMIN)
            IFGRID=IFTMP 
         ENDIF
      ELSE IF(CHOICE.EQ.'REDRAW ISOMETRIC')THEN
            CALL SORTEL
C           Now redraw all the isometric elements.
            DO 160 I=1,NEL
               CALL DRAWIS(ISRT(I))
  160       CONTINUE
      ELSE IF(CHOICE.EQ.'REDRAW MESH')THEN
           CALL REFRESH
           CALL DRMENU('NOCOVER')
           CALL DRGRID
           DO 170 IEL=1,NEL
              CALL DRAWEL(IEL)
  170      CONTINUE
           GOTO 1000
      ELSE IF(CHOICE.EQ.'ZOOM')THEN
           CALL SETZOOM
           CALL REFRESH
           CALL DRMENU('NOCOVER')
           CALL DRGRID
           DO 175 IEL=1,NEL
              CALL DRAWEL(IEL)
  175      CONTINUE
           GOTO 1000
      ELSE IF(CHOICE.EQ.'SET GRID')THEN
           CALL SETGRD
           GOTO 1000
      ELSE IF(CHOICE.EQ.'DEFINE OBJECT')THEN
           CALL SETOBJ
           GOTO 1000
      ELSE IF(CHOICE.EQ.'END    ELEMENTS')THEN
C        WHAT ELSE TO DO WHEN 2-D PROBLEM?
         IF(NEL.EQ.0) THEN
            CALL PRS('A mesh without elements makes for a poor $')
            CALL PRS('Simulation$')
            CALL PRS('ERROR: Can''t "END ELEMENTS" when no elements$')
            CALL PRS('are input. ^C if you want to give up$')
         ELSE
            NLEVEL=0
            DO 180 I=1,NEL
               NLEVEL=MAX(NLEVEL,NUMAPT(I))
  180       CONTINUE
            GOTO 320
         ENDIF
      ELSE IF(CHOICE.EQ.'UP   LEVEL') THEN
C        Make Checks
         NTHISL=0
         NNEXTL=0
C        !!?? maybe 96-32??  also, sometimes up gives wrong ele
         MAXLET=0
         DO 190 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL  )NTHISL=NTHISL+1
            IF(NUMAPT(I).EQ.ILEVEL+1)NNEXTL=NNEXTL+1
            IF(NUMAPT(I).EQ.ILEVEL  )MAXLET=
     $         MAX(MAXLET,ICHAR(LETAPT(I)))
  190    CONTINUE
         IF(NTHISL.EQ.0)THEN
            CALL PRSI('*** ERROR *** No elements on level$',ILEVEL)
            CALL PRS('Can''t make level above empty one$')
            GOTO 1000
         ENDIF
         IF(NNEXTL.EQ.0)THEN
            CALL PRS('** New Level **  Enter Height:$')
            CALL PRS('Negative to Abort start of new level$')
            CALL KEYPAD(HEIGHT(ILEVEL+1))
C           !!?? HOW TO MODIFY HEIGHT ONCE IT'S IN?
            IF(HEIGHT(ILEVEL+1) .LE. 0.0) THEN
C              Abort level change
               CALL PRSI('Aborting level change. Still on $',ilevel)
               HEIGHT(ILEVEL+1)=0
               GOTO 1000
            ELSE IF(HEIGHT(ILEVEL+1) .GT. 0.0) THEN
C              Really go up one
               CALL PRSIS('Using Ceiling of Level$',ILEVEL,
     $         '  Mesh as default$')
C              Erase old mesh
               DO 200 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
  200          CONTINUE
               ILEVEL=ILEVEL+1
               CALL DRELEV(ILEVEL,ILEVEL-1,'     ')
               ILETAP=MAXLET
               NLEVEL=NLEVEL+1
            ENDIF
            IELNOW=NEL
            DO 250 I=1,IELNOW
               IF(NUMAPT(I).EQ.ILEVEL-1)THEN
                   NEL=NEL+1
                   CALL COPYEL(I,NEL)
                   NUMAPT(NEL)=ILEVEL
C                  Correct Stuff affected if walls are not vertical
                   XCEN(NEL)=0.0
                   YCEN(NEL)=0.0
C                  Remove periodic b.c.'s
                   DO 210 IF=1,NFLDS
                      DO 210 ISIDE=1,NSIDES
                         IF(  CBC( ISIDE, NEL,IF).EQ.'P')THEN
                            CBC(   ISIDE, NEL,IF)=' '
                            BC (1, ISIDE, NEL,IF)= 0
                            BC (2, ISIDE, NEL,IF)= 0
                            BC (3, ISIDE, NEL,IF)= 0
                         ENDIF
  210              CONTINUE
                   DO 230 IC=1,4
                      X(NEL,IC)=X(I,IC+4)
                      Y(NEL,IC)=Y(I,IC+4)
                      Z(NEL,IC)=Z(NEL,IC)+HEIGHT(ILEVEL-1)
                      SIDES(NEL,IC,3)=SIDES(NEL,IC,3)+
     $                (HEIGHT(ILEVEL-1)+HEIGHT(ILEVEL))/2.
                      XCEN(NEL)=XCEN(NEL)+X(NEL,IC)/4.0
                      YCEN(NEL)=YCEN(NEL)+Y(NEL,IC)/4.0
C                     Curved side stuff
                      CCURVE(IC,NEL)=CCURVE(IC+4,I)
                      DO 220 II=1,6
                         CURVE(II,IC,NEL)=CURVE(II,IC+4,I)
  220                 CONTINUE
  230              CONTINUE
C                  RAISE Z OF CEILING OF NEW ELEMENTS
                   DO 240 IC=5,8
                      Z(NEL,IC)=Z(NEL,IC)+HEIGHT(ILEVEL)
  240              CONTINUE
                   SIDES(NEL,5,3)=SIDES(NEL,5,3)+HEIGHT(ILEVEL-1)
                   SIDES(NEL,6,3)=SIDES(NEL,6,3)+HEIGHT(ILEVEL)
                   CALL DRAWEL(NEL)
               ENDIF
  250       CONTINUE
C           Sort out the elements so that they will be drawn correctly
            CALL SORTEL
            DO 260 I=1,NEL
               IF(NUMAPT(ISRT(I)).EQ.ILEVEL)CALL DRAWIS(ISRT(I))
  260       CONTINUE
            GOTO 1000
         ELSE
C           Already have elements on next level. Erase old mesh& draw new
            CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
            DO 270 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
  270       CONTINUE
            ILEVEL=ILEVEL+1
            DO 280 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
  280       CONTINUE
            GOTO 1000
         ENDIF
      ELSE IF(CHOICE.EQ.'DOWN LEVEL')THEN
C        Make Checks
         IF(ILEVEL.EQ.1)THEN
            CALL PRS('You already are on Level 1.  You cannot$')
            CALL PRS('go down any further.  Choose something else.$')
            GOTO 1000
         ENDIF
         NTHISL=0
         NDOWNL=0
         MAXLET=0
         DO 290 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL  )NTHISL=NTHISL+1
            IF(NUMAPT(I).EQ.ILEVEL-1)NDOWNL=NDOWNL+1
            IF(NUMAPT(I).EQ.ILEVEL-1)MAXLET=
     $         MAX(MAXLET,ICHAR(LETAPT(I)))
  290    CONTINUE
C        Go down one level.  Erase old mesh& draw new
         CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
         DO 300 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
  300    CONTINUE
         ILEVEL=ILEVEL-1
         DO 310 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
  310    CONTINUE
         GOTO 1000
      ELSE
          CALL PRS('CHOICE:'//CHOICE//'NOT IN MENU$')
      ENDIF
      GOTO 1000
  320 CONTINUE
C     Now, cover menu
      CALL SGVIS(13,0)
      CALL SGVIS(13,1)
      NELF= NEL
C
C     Finished with gridlatch
      IFGRID=.FALSE.
  330 CONTINUE
C     Variable properties put in here
      CALL VPROPS
C     JUMP IN HERE IF YOU DON'T WANT TO MODIFY MESH YOU JUST READ IN
C     Find Sides' Midpoints
      CALL MKSIDE
C     Find Sides' Overlaps (used for finding boundary sides)
C!! Check here for Irrational B.c.'s: unrequited connectivities;
C     Warn if physical boundaries internally.  This check is not exhaustive!!??
C     ZERO OUT OLD ELEMENTAL CONNECTIVITY
      DO 335 Ifld=1,NFLDS
      DO 335 IEL=1,NEL
      DO 335 ISIDE=1,NSIDES
         IF (CBC( ISIDE, IEL,Ifld).EQ.'E') THEN
             CBC(   ISIDE, IEL,Ifld)=' '
             BC (1, ISIDE, IEL,Ifld)= 0
             BC (2, ISIDE, IEL,Ifld)= 0
             BC (3, ISIDE, IEL,Ifld)= 0
         ENDIF
  335 CONTINUE
C
C     Make internal side comparisons
C
      CALL GENCEN
      DO 370 Ifld=1,NFLDS
       DO 365 IEL=1,NEL
        if (nel.gt.100.and.mod(iel,100).eq.0)
     $  write(6,*) 'checking el:',iel
C
        IF (MASKEL(IEL,Ifld).EQ.0) GOTO 365
        DO 360 JEL=1,IEL-1
          IF (MASKEL(JEL,Ifld).EQ.0) GOTO 360
C
          D2 = SQRT( ( xcen(iel)-xcen(jel) )**2 +
     $               ( ycen(iel)-ycen(jel) )**2 +
     $               ( zcen(iel)-zcen(jel) )**2 )
          IF (D2.GT.(RCEN(IEL)+RCEN(JEL))) GOTO 360
C
          DELTA = .001 * MIN(rcen(iel),rcen(jel))
C
          DO 350 ISIDE=1,NSIDES
             IF (CBC(Iside,Iel,Ifld).eq.'E') GOTO 350
             INTERN=0
C
             DO 340 JSIDE=1,NSIDES
                IF (CBC(Jside,Jel,Ifld).eq.'E') GOTO 340
                DELTAX = ABS(SIDES(IEL,ISIDE,1)-SIDES(JEL,JSIDE,1))
                DELTAY = ABS(SIDES(IEL,ISIDE,2)-SIDES(JEL,JSIDE,2))
                DELTAZ = ABS(SIDES(IEL,ISIDE,3)-SIDES(JEL,JSIDE,3))
                IF  (DELTAX .LT. DELTA .AND.
     $               DELTAY .LT. DELTA .AND.
     $               DELTAZ .LT. DELTA ) THEN
C                  BC Array used to define neighboring elements
C                  For want of better notation, 'E' means 
C                                     elemental (internal) bdry
C                  1st two reals are element & side of neighbor; 
C                                     3rd is orientation
C
C                  Re-do internal connectivity
C                  "Normal" internal side.
                   CBC3=CBC( ISIDE, IEL, Ifld)
                   IF(CBC3(3:3) .NE. 'I' .AND. CBC3(3:3) .NE. 'i')THEN
C                     Overlapping edges not Internal B.C.'s are elemental b.c.'s
                      CBC(ISIDE,IEL,Ifld)='E'
                      CBC(JSIDE,JEL,Ifld)='E'
                   ENDIF
                   IORIEN = 1
                   BC (1, ISIDE, IEL, Ifld) = JEL
                   BC (2, ISIDE, IEL, Ifld) = JSIDE
C                  BC (3, ISIDE, IEL, Ifld) = IORIEN
                   BC (1, JSIDE, JEL, Ifld) = IEL
                   BC (2, JSIDE, JEL, Ifld) = ISIDE
C                  BC (3, JSIDE, JEL, Ifld) = IORIEN
                   GOTO 350
                ENDIF
                INTERN=1
  340        CONTINUE
C            If the side is not internal but has an internal B.C., zero it out
             CBC3 = CBC(   ISIDE, IEL,Ifld)
             IF(INTERN.EQ.0 .AND.
     $       (CBC3.EQ.'E' .OR. CBC3(3:3).EQ.'I'.OR. CBC3(3:3).EQ.'i')
     $       )THEN
               CBC(   ISIDE, IEL,Ifld)=' '
               BC (1, ISIDE, IEL,Ifld)= 0
               BC (2, ISIDE, IEL,Ifld)= 0
C              BC (3, ISIDE, IEL,Ifld)= 0
            ENDIF
  350     CONTINUE
  360   CONTINUE
  365  CONTINUE
  370 CONTINUE
      CALL GINGRD(0.0)
C     TURN OFF GRIDDING
      CALL SCROLL(7)
C  ??!! MUST PUT P IN OTHER PERIODIC B.C. IF NEK2.  MUST ACCOUNT FOR
C  CONDUCTION ELEMENTS-- THEY DON'T NEED B.C'S FOR THEIR FLUID
      NSIDES=4
      IF(IF3D)NSIDES=6
      DO 380 Ifld=1,NFLDS
      DO 380 IEL=1,NEL
      DO 380 ISIDE=1,NSIDES
         IF (CBC(ISIDE,IEL,Ifld).EQ.' ') then
           NEEDBC=1
           write(6,*) 'bc:',iel,iside,enew(iel)
         ENDIF
  380 CONTINUE
      CALL BOUND
C
      IFMVBD=.FALSE.
      DO 381 Ifld=1,NFLDS
      DO 381 IEL=1,NEL
      DO 381 ISIDE=1,NSIDES
         CHAR1 = CBC(ISIDE,IEL,Ifld)
         IF(CHAR1.EQ.'M'.OR.CHAR1.EQ.'m') IFMVBD=.TRUE.
  381 CONTINUE
C     Force a dump of mesh location, unless the user specifiaclly turn it
C     off in the output menu
      IF(IFMVBD) IFXYO = .TRUE.
C
      CALL MESGEN
C
      RETURN
      END
      SUBROUTINE READAT
      INCLUDE 'basics.inc'
C                                  Paul's stuff
C     REAL PARAM(100)
      logical iffold,ifhold
      CHARACTER CHTEMP*3
C     NPARAM=20
      NLINF=0
      NLINP=0
      NLINR=0
C
C     Read in all the build,bound,curve,history,output,initcond data
C     Read Dummy Parameters
      READ(9,*,ERR=33)
      READ(9,*,ERR=33)VNEKOLD
      nktold=VNEKOLD
      READ(9,*,ERR=33)
      READ(9,*,ERR=33)NDUM
      DO 188 IDUM=1,NDUM
         read(9,*,err=33)xxx
         if(idum.eq.23)npsold=xxx
188    CONTINUE
C     Read Passive scalar data
      read(9,*,err=33)NSKIP
      IF(NSKIP.GT.0) THEN
         DO 177 I=1,NSKIP
            read(9,*,err=33)
177      CONTINUE
      endif
C
C     READ DUMMY LOGICALS
      READ(9,*,ERR=33)NLOGIC
      DO 189 I=1,NLOGIC
         if(i.eq.1)read(9,*,err=33)iffold
         if(i.eq.2)read(9,*,err=33)ifhold
         if(I.ne.1.and.i.ne.2)READ(9,*,ERR=33)
189   CONTINUE
      NoLDS=1
      IF(ifhold)NoLDS=2+NPSold
C
C     Read Elemental Mesh data
      READ(9,*,ERR=33,end=34)XFAC,YFAC,XZERO,YZERO
      READ(9,*,ERR=33)
      READ(9,*,ERR=33,END=33)NEL,NDIM
c
                    iffmtin = .true.
      IF (NEL.lt.0) iffmtin = .false.
      IF (NEL.lt.0) NEL = -NEL
c
      IF(IF3D.AND.NDIM.EQ.2)THEN
         CALL PRS
     $   ('WARNING: PARAMETER SET TO 2-D BUT MESH DATA IS FOR 3-D$')
         CALL PRS('PLEASE ENTER Z-HEIGHT: $')
         CALL RER( ZHGHT)
      ENDIF
      IF(.NOT.IF3D.AND.NDIM.EQ.3)
     $CALL PRS('ERROR: PARAMETER SET TO 3-D BUT MESH DATA IS FOR 2-D$')
      NCOND=0
      nel50=nel/50
      if (nel.lt.100) nel50=500
      DO 98 IEL=1,NEL
         READ(9,'(20X,I4,4X,I3,A1,11x,i5)',ERR=33,END=33)
     $   IDUM,NUMAPT(IEL),LETAPT(IEL),IGROUP(IEL)
         IF(IGROUP(IEL).GT.0) NCOND=NCOND+1
         IF (NEL.LT.100.or.mod(iel,nel50).eq.0) THEN
            WRITE(S,'(A20,I4,A4,I3,A1,A1)',ERR=33)
     $    ' Reading Element    ',IDUM,'   [',NUMAPT(IEL),LETAPT(IEL),']'
            CALL PRS(S//'$')
         ENDIF
         IF(NDIM.EQ.2)THEN
            READ(9,*,ERR=33,END=33)(X(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(Y(IEL,IC),IC=1,4)
            IF (IF3D) THEN
               DO 197 IC=1,4
                  I4 = IC+4
                  Z(IEL,IC) = 0.0
                  X(IEL,I4) = X(IEL,IC)
                  Y(IEL,I4) = Y(IEL,IC)
                  Z(IEL,I4) = ZHGHT
  197          CONTINUE
            ENDIF
         ELSE IF(NDIM.EQ.3)THEN
            READ(9,*,ERR=33,END=33)(X(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(Y(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(Z(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(X(IEL,IC),IC=5,8)
            READ(9,*,ERR=33,END=33)(Y(IEL,IC),IC=5,8)
            READ(9,*,ERR=33,END=33)(Z(IEL,IC),IC=5,8)
C           For Flat floors only.  On write will flatten floors.
         ENDIF
98    CONTINUE
      NELF=NEL-NCOND
C     Read curved side data
      READ(9,*,ERR=57,END=57)
      READ(9,*,ERR=57,END=57)NCURVE
      IF(NCURVE.GT.0)THEN
         DO 19 I=1,NCURVE
            if (nel.lt.1000) then
               READ(9,'(I3,I3,5G14.6,1X,A1)',ERR=57,END=57)
     $         IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            else
               READ(9,'(I2,I6,5G14.6,1X,A1)',ERR=57,END=57)
     $         IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            endif
            CALL DDUMMY(IEDGE,IEL)
            CURVE(1,IEDGE,IEL)=R1
            CURVE(2,IEDGE,IEL)=R2
            CURVE(3,IEDGE,IEL)=R3
            CURVE(4,IEDGE,IEL)=R4
            CURVE(5,IEDGE,IEL)=R5
            CCURVE( IEDGE,IEL)=ANS
19       CONTINUE
      ENDIF
C     Check here if the old data has the same fields as the new stuff.
C     if not, skip the b.c.'s, etc to avoid confusion
c...skip, pff        IF( (IFFLOW.AND..NOT.IFFOLD).OR.(.NOT.IFFLOW.AND.IFFOLD)
c...skip, pff       $.OR.(IFHEAT.AND..NOT.IFHOLD).OR.(.NOT.IFHEAT.AND.IFHOLD)
c...skip, pff       $.OR. NPSOLD .NE. NPSCAL .or. NKTOLD .ne. NKTONV) THEN
c...skip, pff           CALL PRS
c...skip, pff       $ ('Old B.C.''s and options were for different equations.$')
c...skip, pff           CALL PRS
c...skip, pff       $ ('Using mesh data only; ignoring old B.C''s and options$')
c...skip, pff           CALL PRS('Difference, New, Old:$')
c...skip, pff           IF(NPSOLD.NE.NPSCAL) THEN
c...skip, pff               CALL PRS('Passive Scalars $')
c...skip, pff               CALL PRII(npsold,NPSCAL)
c...skip, pff           ENDIF
c...skip, pff           IF(NKTOLD.NE.NKTONV) THEN
c...skip, pff              CALL PRS('Nekton Version $')
c...skip, pff              CALL PRII(nktold,nktonv)
c...skip, pff           ENDIF
c...skip, pff           RETURN
c...skip, pff        ENDIF
C     Read Boundary Conditions (and connectivity data)
      if (iffmtin) READ(9,*,ERR=44,END=44)
      DO 90 IFLD=1,NoLDS
C        Fluid and/or thermal
C        !!?? NELF DIFFERENT FROM NEL??
         if (iffmtin) READ(9,*,ERR=44,END=44)
         IF( (IFLD.EQ.1.AND.IFFOLD) .OR. IFLD.GT.1)THEN
C           FIX UP FOR WHICH OF FIELDS TO BE USED
            DO 88 IEL=1,NEL
            DO 88 ISIDE=1,NSIDES
C              !Fix to a4,i2 when you make cbc character*4
               IF(VNEKOLD .LE. 2.5) NBCREA = 3
               IF(VNEKOLD .GE. 2.6) NBCREA = 5
               IF (NEL.LT.1000.and.iffmtin) THEN
                  READ(9,'(1X,A3,1x,I2,I3,5G14.6)',ERR=44,END=44)
     $            CBC(ISIDE,IEL,IFLD),ID,ID,
     $            (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
               ELSEIF (iffmtin) then
                  READ(9,'(1X,A3,I5,I1,5G14.6)',ERR=44,END=44)
     $            CBC(ISIDE,IEL,IFLD),ID,ID,
     $            (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
               ELSE
                  READ(8,ERR=44,END=44) chtmp3,
     $            CBC(ISIDE,IEL,IFLD),ID,ID,
     $            (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
               ENDIF
               CBC1=CBC(ISIDE,IEL,IFLD)
               IF(CBC1.EQ.'M'.OR.CBC1.EQ.'m')THEN
                 IFMOVB=.TRUE.
                 IFGNGM=.TRUE.
               ENDIF
               ICBC1=ICHAR(CBC1)
               IF(ICBC1.GE.97 .AND. ICBC1.LE.122)THEN
                  CBC3 = CBC(ISIDE,IEL,IFLD)
                  IF(CBC3(3:3) .EQ. 'i')THEN
C                    Special storage for internal boundaries
                     NLINES=BC(4,ISIDE,IEL,IFLD)
                     BC(5,ISIDE,IEL,IFLD)=LOCLIN
                  ELSE
                     NLINES=BC(1,ISIDE,IEL,IFLD)
                     BC(2,ISIDE,IEL,IFLD)=LOCLIN
                  ENDIF
                  DO 86 I=1,NLINES
                    READ(9,'(A70)',ERR=44,END=44)INBC(LOCLIN)
                    LOCLIN=LOCLIN+1
86                CONTINUE
               ENDIF
88          CONTINUE
            DO 89 IEL=1,NEL
              IF(IFMOVB)ICRV(IEL)=1+4+9+16
              DO 89 IEDGE=1,8
                 IF(IFMOVB)CCURVE(IEDGE,IEL)='M'
89          CONTINUE
         ENDIF
90    CONTINUE
      IF(NFLDS.EQ.1.and.iffmtin)READ(9,*,err=45,end=45)
c... newxC     Read Boundary Conditions (and connectivity data)
c... newx      READ(9,*,ERR=35,END=35)
c... newx      DO 90 IFLD=1,NoLDS
c... newxC        Fluid and/or thermal
c... newxC        !!?? NELF DIFFERENT FROM NEL??
c... newx         if ((iffold.and.ifld.eq.1).or.ifld.gt.1) then
c... newx         READ(9,*,ERR=35,END=35)
c... newxC        No Fluid
c... newx         IF(IFLD.EQ.1 .AND.  .NOT. IFFLOW)GO TO 89
c... newxC        FIX UP FOR WHICH OF FIELDS TO BE USED
c... newx         DO 88 IEL=1,NEL
c... newx         DO 88 ISIDE=1,NSIDES
c... newx            IF(VNEKOLD .LE. 2.5) NBCREA = 3
c... newx            IF(VNEKOLD .GE. 2.6) NBCREA = 5
c... newx            if (nel.lt.1000) then
c... newx               READ(9,'(A1,A3,1X,I2,I3,5G14.6)',ERR=35,END=35)
c... newx     $         CHTEMP,
c... newx     $         CBC(ISIDE,IEL,IFLD),ID,ID,
c... newx     $         (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
c... newx            else
c... newx               READ(9,'(A1,A3,I5,I1,5G14.6)',ERR=35,END=35)
c... newx     $         CHTEMP,
c... newx     $         CBC(ISIDE,IEL,IFLD),ID,ID,
c... newx     $         (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
c... newx            endif
c... newxC
c... newx            IF(IFLD.EQ.1 .OR. (IFLD.EQ.2 .AND. .NOT. IFFLOW))
c... newx     $      CBC(ISIDE,IEL,0)= CHTEMP
c... newxC
c... newx            CBC1=CBC(ISIDE,IEL,IFLD)
c... newx            CBC3=CBC(ISIDE,IEL,IFLD)
c... newx            ICBC1=ICHAR(CBC1)
c... newx            IF(ICBC1.GE.97 .AND. ICBC1.LE.122)THEN
c... newxC              Read in # of lines and pointer to first. THen read in lines.
c... newx               IF(CBC3(3:3) .EQ. 'i')THEN
c... newxC                 Internal B.C.'s, unfortunately, must be handled specially
c... newxC                 BC locations 1 and 2 are used for direct stiffness pointers
c... newxC                 4 and 5 are the only locations available for line skips.
c... newx                  NLINES=BC(4,ISIDE,IEL,IFLD)
c... newxC                 Re-store the LOCLIN pointer (dont believe the read file)
c... newx                  BC(5,ISIDE,IEL,IFLD)=LOCLIN
c... newx               ELSE
c... newxC                 "Normal" storage locations for line skip pointers
c... newx                  NLINES=BC(1,ISIDE,IEL,IFLD)
c... newxC                 Re-store the LOCLIN pointer (dont believe the read file)
c... newx                  BC(2,ISIDE,IEL,IFLD)=LOCLIN
c... newx               ENDIF
c... newx               DO 86 I=1,NLINES
c... newx                 READ(9,'(A70)',ERR=35,END=35)INBC(LOCLIN)
c... newx                 LOCLIN=LOCLIN+1
c... newx                 IF(LOCLIN.GE.MAXIBC)
c... newx     $           CALL PRS('***WARNING**** RAN OUT OF'//
c... newx     $           'SPACE FOR INFLOW B.C.s.  NO MORE ALLOWED!$')
c... newx86             CONTINUE
c... newx            ENDIF
c... newx88        CONTINUE
c... newx89        CONTINUE
c... newx90    CONTINUE
c... newxC     No Thermal B.C.'s
c... newx      IF(NFLDS.EQ.1) READ(9,*,ERR=35,END=35)
C     Read Initial Conditions
      IF(VNEKOLD.LT.2.5)THEN
         CALL PRS('***** OLD VERSION RESTART DATA ****$')
         CALL PRS(
     $   'Please Re-enter Initial conditions from the options menu$')
          READ(9,*,ERR=50,END=50)
          DO 21 I=1,7
21          READ(9,'(A80)',ERR=50,END=50)
      ELSE
C         Check for 2.6 style read file - separate restart/fortran sections.
          READ(9,'(A70)',ERR=50,END=50) LINE
          IF (INDX1(LINE,'RESTART',7).NE.0) THEN
             REWIND(13)
             WRITE(13,'(A70)')LINE
             REWIND(13)
             READ (13,*,ERR=50,END=50) NLINR
             REWIND(13)
             NSKIP=NLINR
             DO 201 I=1,NSKIP
                READ(9,'(A80)')INITP(I)
                LINE=INITP(I)
                CALL CAPIT(LINE,70)
                IF (INDX1(LINE,'PRESOLV',7).NE.0) THEN
                   NLINP=1
                   NLINR=NLINR-1
                ENDIF
201          CONTINUE
             READ(9,'(A70)',ERR=50,END=50) LINE
          ENDIF
C         Read fortran initial condition data.
          REWIND(13)
          WRITE(13,'(A70)')LINE
          REWIND(13)
          READ(13,*,ERR=50,END=50)NLINF
          REWIND(13)
          DO 22 I=1,NLINF
             READ(9,'(A80)')INITC(I)
22        CONTINUE
      ENDIF
C     Read drive force data
      READ(9,*,ERR=50,END=50)
      READ(9,*,ERR=50,END=50)nskip
      do 23 i=1,nskip
23          READ(9,'(A80)',ERR=50,END=50)DRIVC(I)
C     READ CONDUCTION ELEMENT DATA
      READ(9,*,ERR=60,END=60)
      NFLDSC=0
      DO 117 IFLD=1,NFLDS
         IF(IFTMSH(IFLD))NFLDSC=NFLDSC+1
117   CONTINUE
      READ(9,*,ERR=60,END=60)NSKIP
      READ(9,*,ERR=60,END=60)NPACKS
      IF(NPACKS.GT.0)THEN
      DO 218 IIG=1,NPACKS
            READ(9,*)IGRP,IF,ITYPE
            MATYPE(IGRP,IF)=ITYPE
            DO 218 IPROP=1,3
               IF(ITYPE.EQ.1) READ(9,*      ) CPROP(IGRP,IF,IPROP)
               IF(ITYPE.EQ.2) READ(9,'(A80)') VPROP(IGRP,IF,IPROP)
218   CONTINUE
      ENDIF
C
C     Read history data
      READ(9,*,ERR=50,END=50)
      READ(9,*,ERR=50,END=50)NHIS
C     HCODE(11) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
      IF(NHIS.GT.0)THEN
         DO 51 I=1,NHIS
            READ(9,'(1X,11A1,1X,4I5)',ERR=50,END=50)
     $      (HCODE(II,I),II=1,11),(LOCHIS(II,I),II=1,4)
51       CONTINUE
      ENDIF
C     Read output specs
      READ(9,*,ERR=50,END=50)
      READ(9,*,ERR=50,END=50)NOUTS
      READ(9,*,ERR=50,END=50)IFXYO
      READ(9,*,ERR=50,END=50)IFVO
      READ(9,*,ERR=50,END=50)IFPO
      READ(9,*,ERR=50,END=50)IFTO
      READ(9,*,ERR=50,END=50)IFTGO
      READ(9,*,ERR=50,END=50)IPSCO
      IF(IPSCO .GT.0)THEN
         DO 1221 I=1,IPSCO
            READ(9,'(1X,L1,1X,A5)',ERR=50,END=50) IFPSCO(I),PSNAME(I)
1221     CONTINUE
      ENDIF
C     OBJECT SPECIFICATION DATA
      READ(9,*,ERR=50,END=50)
      READ(9,*,ERR=50,END=50)NSOBJS
      IF(NSOBJS .GT. 0)THEN
          DO 62 IOBJ=1,NSOBJS
             READ(9,'(4X,I4,35X,A20)',ERR=50,END=50)
     $       NFACE(IOBJ),SOBJ(IOBJ)
             DO 63 IFACE=1,NFACE(IOBJ)
                READ(9,*,ERR=50,END=50)
     $          ILSURF(1,IOBJ,IFACE),ILSURF(2,IOBJ,IFACE)
63           CONTINUE
62        CONTINUE
      ENDIF
      READ(9,*,ERR=50,END=50)NVOBJS
      READ(9,*,ERR=50,END=50)NEOBJS
      READ(9,*,ERR=50,END=50)NPOBJS
C
      RETURN
33    CALL PRS
     $('ERROR: MESH DATA MISSING OR CORRUPTED.  CHECK READ FILE.$')
      CALL BEEP
      CALL PRS('<CR> TO CONTINUE$')
      CALL RES(LINE,0)
      CALL DEVEX
      CALL EXIT
34    CALL PRS(
     $'ERROR: FILE CONTAINS ONLY PARAMETER DATA. MESH DATA MISSING$')
      CALL BEEP
      CALL PRS('<CR> TO CONTINUE$')
      CALL RES(LINE,0)
      CALL DEVEX
      CALL EXIT
35    CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $'  CONTINUEING WITH MESH ONLY.$')
      RETURN
44    WRITE(S,'(1X,A36,I4,A6,I3)')
     $'Error reading B.C. data for element ',IEL,' side ',ISIDE
      CALL PRS(S//'$')
      RETURN
45    CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $'  CONTINUEING WITH MESH ONLY.$')
      RETURN
40    CALL PRS('CURVED SIDES NOT READ CORRECTLY.  PLEASE CHECK$')
      CALL PRS('SIDES IN CURVED SIDE MENU$')
50    CALL PRS
     $('OPTIONS MAY NOT HAVE BEEN READ CORRECTLY.  PLEASE CHECK$')
      CALL PRS('THEM IN OPTIONS MENU$')
      RETURN
57    CALL PRS
     $('CURVED SIDES HAVE NOT BEEN READ CORRECTLY.  PLEASE REVIEW$')
      CALL PRS('THEM.$')
      RETURN
60    CALL PRS
     $('CONDUCTION ELEMENTS MAY NOT HAVE BEEN READ CORRECTLY.  $')
      CALL PRS('PLEASE CHECKTHEM IN OPTIONS MENU$')
70    CALL PRS('NO OBJECTS$')
      RETURN
      END
      SUBROUTINE SETSCL
C     Sets scale factors
C
      INCLUDE 'basics.inc'
      LOGICAL IFMENU
      REAL XSC(4),YSC(4)
      GRID=PARAM(18)
3     CONTINUE
      CALL CLEAR
      CALL REFRESH
      CALL DRGRID
c     CALL DRMENU('NOCOVER')
      CALL PRS(' SET SCALE FACTORS $')
c     NCHOIC =  3
      ITEM(1)='MOUSE'
      ITEM(2)='TYPE'
      ITEM(3)='USE NORMALIZED'
      CALL PRS(
     $'Scale Factors: Input thru Mouse, type in, or use normalized?$')
      IFCEN=.TRUE.
C
C     10/7/92 pff .... always use MOUSE
C
c     CALL MENU(XMOUSE,YMOUSE,BUTTON,'SCALE FACTOR')
c     IF(IFNOSEG)CALL DRCOVR(13)
      CHOICE=ITEM(1)
c
c
c
      IFCEN=.FALSE.
C     Begin using grid latch filter
      IFGRID=.TRUE.
      IF(CHOICE.EQ.'TYPE') THEN
          CALL PRS('Enter Scale Factors Xfac,Yfac:$')
          CALL PRS(
     $    'They are the lengths that correspond to 1 screen height$')
          CALL RERR(XFAC,YFAC)
          IF(IFLEARN)WRITE(3,*)XFAC,YFAC
          XFAC=ABS(XFAC)
          YFAC=ABS(YFAC)
C
          CALL PRS
     $('Enter the coordinates of Left Bottom corner of screen$')
          CALL RERR(XZERO,YZERO)
          IF(IFLEARN)WRITE(3,*)XZERO,YZERO
C         Mark Zero
          CALL COLOR(7)
      ELSE IF(CHOICE.EQ.'MOUSE') THEN
   17     CALL PRS('Push yellow button with mouse at 2 different'//
     $    ' x (horizontal)$')
          CALL PRS('locations$')
          DO 18 I=1,4
              IF(I.EQ.2)CALL PRS('2nd x location$')
              IF(I.EQ.4)CALL PRS('2nd y location$')
              IF(I.EQ.3)THEN
                 CALL PRS('Enter corresponding x length$')
                 CALL GINGRD(0.0)
                 CALL KEYPAD(XLEN)
                 CALL GINGRD(GRID)
                 CALL PRS('Push yellow button with mouse at 2'//
     $           ' different y (vertical) locations$')
                 CALL PRS('or hit keypad for equal X-Y scaling.$')
              ENDIF
              CALL MOUSE(XSC(I),YSC(I),BUTTON)
C
C             Check for keypad entry?
C
              IF (IFMENU(XSC(I),YSC(I)).AND.I.LE.2) THEN
                 CALL PRS('Please enter two points in build area.$')
                 GOTO 17
              ENDIF
              IF (IFMENU(XSC(I),YSC(I)).AND.I.GT.2) THEN
                 CALL PRS('Assuming equal scaling in X and Y.$')
                 YSC(4)=XSC(2)
                 YSC(3)=XSC(1)
                 YLEN=XLEN
                 ICOUNT=I-1
                 GOTO 19
              ENDIF
C
              CALL COLOR(7)
              CALL MOVE(XSC(I)-.02,YSC(I))
              CALL DRAW(XSC(I)+.02,YSC(I))
              CALL MOVE(XSC(I),YSC(I)-.02)
              CALL DRAW(XSC(I),YSC(I)+.02)
              ICOUNT=I
18        CONTINUE
          CALL PRS('Enter in corresponding y length$')
          CALL GINGRD(0.0)
          CALL KEYPAD(YLEN)
          CALL GINGRD(GRID)
19        CONTINUE
          IF(XSC(2).EQ.XSC(1) .OR. YSC(4).EQ.YSC(3))THEN
             CALL PRS('ERROR: MUST BE NONZERO GAP.  START OVER$')
             GO TO 3
          ENDIF
C
C         Compute and check scale factors
C
          XFACT=ABS(XLEN/(XSC(2)-XSC(1)))
          YFACT=ABS(YLEN/(YSC(4)-YSC(3)))
          IF(XFACT.EQ.0.0 .OR. YFACT.EQ.0.0)THEN
              CALL PRS('ERROR: MUST BE NONZERO SCALE FACTORS '//
     $                 'START OVER$')
              GO TO 3
          ENDIF
C
          DO 20 I=1,ICOUNT
C             Erase markers
              CALL COLOR(0)
              CALL MOVEC(XSC(I)-.02,YSC(I))
              CALL DRAWC(XSC(I)+.02,YSC(I))
              CALL MOVEC(XSC(I),YSC(I)-.02)
              CALL DRAWC(XSC(I),YSC(I)+.02)
20        CONTINUE
C
C===============================================================
C         Set Origin
C===============================================================
C
          CALL PRS (' $')
          CALL PRS
     $    ('Push yellow button with mouse at known point in x,y$')
          CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
C
          CALL PRS('Enter X coordinate of point:$')
          CALL GINGRD(0.0)
          CALL KEYPAD(XPHYS)
          CALL PRS('Enter Y coordinate of point:$')
          CALL KEYPAD(YPHYS)
          CALL GINGRD(GRID)
C
          XFAC  = ABS(XLEN/(XSC(2)-XSC(1)))
          YFAC  = ABS(YLEN/(YSC(4)-YSC(3)))
          XZERO = XPHYS - XFAC * XMOUSE
          YZERO = YPHYS - YFAC * YMOUSE
C===============================================================
      ELSE IF(CHOICE.EQ.'USE NORMALIZED') THEN
          XFAC=  2.0
          YFAC=  2.0
          XZERO =-1.
          YZERO =-1.
      ELSE
         CALL PRS('Error Reading Input$')
         GO TO 3
      ENDIF
C
C===============================================================
C     Set grid stuff
C===============================================================
      GRIDDX=XFAC*GRID
      GRIDDY=YFAC*GRID
C     Set clipping window based on scale factors
      WT = YPHY(.995)
      WB = YPHY(.005)
      WL = XPHY(.005)
      WR = XPHY(.995)
C     Save scale factors for UnZoom function
      XFACO=XFAC
      YFACO=YFAC
      XZEROO=XZERO
      YZEROO=YZERO
C
      CALL REFRESH
      CALL DRGRID
      RETURN
      END
