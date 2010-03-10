c-----------------------------------------------------------------------
C     
C     NEKTON 2.6  2/8/90
C     
C     Copyright (C) 1990, by the 
C     
C     Massachusetts Institute of Technology  and Nektonics, Inc.
C     
C     All Rights Reserved
C     
C     This program is a licenced product of MIT and Nektonics, Inc.,  and it is 
C     not to be disclosed to others, copied, distributed, or displayed 
C     without prior authorization.
C     
C     
c-----------------------------------------------------------------------
      subroutine build
C     
C     Menu-based module that prompts the user to input corners.
C     It then establishes connectivity, sets global node numbers, and
C     constructs elemental mesh. Calls routines that set boundary
C     conditions.
C     234567890123456789012345678901234567890123456789012345678901234567890123
C     
      include 'basics.inc'
      DIMENSION ICRVS(4)
C     ELEMENT,CORNER,COOORDINATE#,LEVEL (OLD OR CURRENT)
      CHARACTER KEY,STRING*6,LETO,CHAR1
      LOGICAL IFTMP
      COMMON /SPLITT/ ENEW(NELM),IND(NELM)
C     
      PI=4.0*ATAN(1.0)
c     NX=PARAM(20)
      nx=nxm
      IF (NX.GT.NXM) THEN
         WRITE(S,104) NX,NXM
 104     FORMAT(' Warning, current N exceeds NXM, resetting from',I3
     $        ,' to',I3,'.$')
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
 10      CONTINUE
      ELSE
C     Need legendre points for mesh drawn on screen
         CALL LEGEND(XXPTS,WGHT,NX)
         CALL LEGEND(YYPTS,WGHT,NY)
      ENDIF
      NSIDES=4
      IFCEN=.FALSE.
      IF(IF3D)NSIDES=6
      IF(CHOICE.EQ.'BUILD FROM FILE'.OR.
     $     CHOICE.EQ.'IMPORT UNIVERSAL FILE')THEN
         IF(CHOICE.EQ.'BUILD FROM FILE') then
           call readat
           close(unit=9)
           if (.not.if3d) call chk_right_hand(nel)
         endif
C     Set grid stuff based on old data
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
C     Now set stuff based on data from readat
         DO 20 IEL=1,NEL
            HEIGHT(NUMAPT(IEL))=Z(IEL,5)-Z(IEL,1)
            XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.0
            YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.0
 20      CONTINUE
C     Display Mesh (or at least 1st floor)
         CALL DRGRID
         NLEVEL=0
C     Sort out the elements so that they will be drawn correctly
         CALL SORTEL
         DO 30 IEL=1,NEL
            CALL DRAWIS(ISRT(IEL))
            IF(NUMAPT(IEL).EQ.1)CALL DRAWEL(IEL)
            IF(NUMAPT(IEL).GT.NLEVEL)NLEVEL=NUMAPT(IEL)
 30      CONTINUE
C     
         IF(NDIM.EQ.2)NLEVEL=1
C     Display Elevator 1st floor hilighted
         IF(NLEVEL.GT.1)THEN
            DO 40 I=NLEVEL,2,-1
               CALL DRELEV(I-1,I,'     ')
 40         CONTINUE
         ENDIF
         CALL SCROLL(5)
         CALL GNABLE
         ITEM(1)='ACCEPT MESH'
         ITEM(2)='REVIEW/MODIFY'

         goto 331  ! quick hack


c        ITEM(3)='Edit Mesh'
         NCHOIC =  2
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'ACCEPT/REVIEW')
         IF(IFNOSEG)CALL DRCOVR(13)
         IF(CHOICE.EQ.'ACCEPT MESH')THEN
            GOTO 330
         ELSE IF(CHOICE.EQ.'Edit Mesh')THEN
            call mesh_edit
         ELSE IF(CHOICE.EQ.'REVIEW/MODIFY')THEN
c           prepare to modify floor by floor
         ENDIF

  331    continue   ! quick hack


      ELSE
C     Interactive Input
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
C     Put in bogus height of first level
         ELSE
C     For interactive sesion, Demand Height of First Level
            CALL PRSI('Please Enter Height of Level$',ILEVEL)
 50         call rer(HEIGHT(ILEVEL))
            IF(HEIGHT(ILEVEL).LE.0.0)THEN
               CALL PRS('HEIGHT must be a positive number!$')
               GOTO 50
            ENDIF
         ENDIF
C     2nd arg: 0 for no fade; negative for bottom only (for floor modify)
         CALL DRELEV(ILEVEL,0,'     ')
      ENDIF
C     Draw FIRST Level
      IF(NEL.GT.0) THEN
         DO 60 IEL=1,NEL
C     CALL DRAWIS(ISRT(IEL))
C     IF(NUMAPT(IEL).EQ.ILEVEL)CALL DRAWEL(IEL)
            IF (NUMAPT(IEL).EQ.ILEVEL)
     $           MAXLET=MAX(MAXLET,ICHAR(LETAPT(IEL)))
 60      CONTINUE
      ENDIF
      ILETAP=MAXLET+96-32
C     ! ??!!
C     
 1000 CONTINUE
      call gencen
C     
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
c        nchoic = nchoic+1
c        ITEM(nchoic)        =             'REDRAW ISOMETRIC'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'REDRAW MESH'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'SET GRID'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'ZOOM'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'IMPORT MESH'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'IMPORT vtx MESH'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'IMPORT vtk MESH'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'REFLECT MESH '
c        nchoic = nchoic+1
c        ITEM(nchoic)        =             'RSB'
         nchoic = nchoic+1
         ITEM(nchoic)        =             'DEFINE OBJECT'
         IF(IFCEIL)THEN
C     Restrict choices for person on CEILING
            ITEM(1)='OFF CEILING'
            ITEM(2)='MODIFY ELEMENT'
            ITEM(3)='CURVE SIDES'
            NCHOIC=3
         ENDIF
      ELSE
C     2-D
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
c        ITEM(nchoic)       =             'Edit Mesh'
c        nchoic = nchoic+1
         ITEM(nchoic)       =             'REDRAW MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'IMPORT MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'IMPORT vtk MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'IMPORT vtx MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =             'REFLECT MESH '
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
         CALL PRS('Enter element corners. Use mouse,$')
         CALL PRS('or click menu area to type x-y coordinate$')
C     Turn on Keypad
         CALL NEWEL(NEL,XMOUSE,YMOUSE,BUTTON,IERR)
         IF (IERR.EQ.1) THEN
            CALL DRAWEL(-NEL)
            NEL=NEL-1
         ENDIF
      ELSE IF(CHOICE.EQ.'CEILING') THEN
C     Make sure everythingh is drawn on ceiling.  And that
C     MODEL and CURVE know about it, too
         DO 80 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 80      CONTINUE
         IFCEIL=.TRUE.
         DO 90 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
 90      CONTINUE
         CALL DRELEV(-(ILEVEL+1),ILEVEL,'     ')
         GOTO 1000
      ELSE IF(CHOICE.EQ.'OFF CEILING') THEN
C     Make sure everythingh is drawn on ceiling.  And that
C     MODEL and CURVE know about it, too
         DO 100 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 100     CONTINUE
         IFCEIL=.FALSE.
         DO 110 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
 110     CONTINUE
         CALL DRELEV(ILEVEL,-(ILEVEL+1),'     ')
         GOTO 1000
      ELSE IF(CHOICE.EQ.'IMPORT MESH')THEN
         call imp_mesh(.true.)
         call redraw_mesh_small
         GOTO 1000
      ELSE IF(CHOICE.EQ.'REFLECT MESH')THEN
         CALL REFLECT_MESH
         GOTO 1000
      ELSE IF(CHOICE.EQ.'IMPORT vtk MESH')THEN
         call imp_mesh_vtk
         GOTO 1000
      ELSE IF(CHOICE.EQ.'IMPORT vtx MESH')THEN
         call imp_mesh_vtx
         GOTO 1000
      ELSE IF(CHOICE.EQ.'CURVE SIDES')THEN
         CALL CURVES
         GOTO 1000
      ELSE IF(CHOICE.EQ.'MODIFY ELEMENT')THEN
C     Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
            GOTO 1000
         ELSE
            CALL MODEL(NEL)
         ENDIF
      ELSE IF(CHOICE.EQ.'GLOBAL REFINE')THEN
C     Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
            GOTO 1000
         ELSE
C     Only floor of elevator hilighted during modify
            CALL DRELEV(-ILEVEL,ILEVEL,'     ')
            CALL GLOMOD
            CALL DRELEV(ILEVEL,0,'     ')
         ENDIF
      ELSE IF(CHOICE.EQ.'DELETE ELEMENT')THEN
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to delete$')
            GOTO 1000
         ELSE
C     Find out which element to delete
            CALL PRS(
     $           'Enter (w/mouse) 2 points in element to be deleted,$')
            CALL PRS(
     $           'or, 2 points framing a box containing elements.$')
            CALL PRS(
     $           'Enter in menu area to abort DELETE ELEMENT op.$')
            IFTMP =IFGRID
            IFGRID=.FALSE.
 120        CONTINUE
            CALL PRS('Enter 1st point:$')
            CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
            IF (XMOUSE.GT.XPHY(1.0)) THEN
C     look for a keypad input
               ifgrid=iftmp 
               call beep
               call prs(
     $              'Abort DELETE? (y/s=special):$')
c    $              'Select item, or DELETE ELEMENT to try again.$')
               call res(ans,1)
               if (ans.eq.'y'.or.ans.eq.'Y') goto 1000
               if (ans.eq.'s'.or.ans.eq.'S') call special_delete
               goto 1000
            ELSE
               CALL PRS('Enter 2nd point:$')
               CALL MOUSE(XMOUS2,YMOUS2,BUTTON)
               IF (XMOUS2.GT.XPHY(1.0)) THEN
                  ifgrid=iftmp 
                  call beep
                  call prs(
     $                 'Abort DELETE? (y/s=special):$')
c    $                 'Select item, or DELETE ELEMENT to try again.$')
                  call res(ans,1)
                  if (ans.eq.'y'.or.ans.eq.'Y') goto 1000
                  if (ans.eq.'s'.or.ans.eq.'S') call special_delete
                  goto 1000
               ENDIF
            ENDIF
C     
C     We successfully input 2 points in the build area
C     1) count number of centroids in the bounding box
C     2) if ncntrd=0, find nearest element and anihilate it.
C     
C     Box
            XMAX=MAX(XMOUSE,XMOUS2)
            XMIN=MIN(XMOUSE,XMOUS2)
            YMAX=MAX(YMOUSE,YMOUS2)
            YMIN=MIN(YMOUSE,YMOUS2)
C     
C     Check box to see if it contains any elements, and delete them.
C     
C     This is a gross N^2 algorithm, but it beats entering them all
C     by hand....
C     
            NUMDEL=0
 124        CONTINUE
            NELT=NEL
            DO 125 JEL=1,NELT
               IF (JEL.LE.NEL        .AND.
     $              XMIN.LE.XCEN(JEL) .AND.
     $              XCEN(JEL).LE.XMAX .AND.
     $              YMIN.LE.YCEN(JEL) .AND.
     $              YCEN(JEL).LE.YMAX ) THEN
                  IF (NUMDEL.EQ.0) 
     $                 CALL DRWBOX(XMIN,YMIN,XMAX,YMAX,1)
                  CALL DELEL(JEL)
                  NUMDEL=NUMDEL+1
                  GOTO 124
               ENDIF
 125        CONTINUE
            IF (NUMDEL.GT.0) THEN
C     redraw the mesh
               CALL REFRESH
               CALL DRMENU('NOCOVER')
               CALL DRGRID
               DO 128 IEL=1,NEL
                  CALL DRAWEL(IEL)
 128           CONTINUE
               IFGRID=IFTMP 
               GOTO 1000
            ELSE
C     Look for closest element (standard delete element option)
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
 130           CONTINUE
            ENDIF
C     Check if it$')s OK to delete
            IF(IF3D)THEN
               IEQ=0
               IGT=0
               DO 140 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL  )IEQ=IEQ+1
                  IF(NUMAPT(I).EQ.ILEVEL+1)IGT=IGT+1
 140           CONTINUE
               IF(IEQ.EQ.1)THEN
C     Last element on floor
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
C     Now, black out old element
            CALL DRAWEL(-IELMIN)
            CALL DRAWIS(-IELMIN)
C     Now draw any elements with higher priority
C     ISRT contains the numbers of the elements in order of
C     plotting priority.  Element (ISRT(NEL)) is most visible.
            CALL DELEL(IELMIN)
C     Sort out the elements so that they will be drawn correctly
            CALL SORTEL
C     Now redraw all the isometric elements.
            DO 150 I=1,NEL
C     CALL DRAWIS(ISRT(I))
 150        CONTINUE
            IF(.NOT.IF3D .AND.IELMIN.LE.NEL)CALL DRAWEL(IELMIN)
            IFGRID=IFTMP 
         ENDIF
      ELSE IF(CHOICE.EQ.'REDRAW ISOMETRIC')THEN
         CALL SORTEL
C     Now redraw all the isometric elements.
         DO 160 I=1,NEL
            CALL DRAWIS(ISRT(I))
 160     CONTINUE
      ELSE IF(CHOICE.EQ.'REDRAW MESH')THEN
c     
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         DO 170 IEL=1,NEL
            CALL DRAWEL(IEL)
 170     CONTINUE
c     
C     Now redraw all the isometric elements.
         call sortel
         do i=1,nel
            call drawis(isrt(i))
         enddo
c     
         GOTO 1000
c     
      ELSE IF(CHOICE.EQ.'ZOOM')THEN
         CALL SETZOOM
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         DO 175 IEL=1,NEL
            CALL DRAWEL(IEL)
 175     CONTINUE
         GOTO 1000
      ELSE IF(CHOICE.EQ.'SET GRID')THEN
         CALL SETGRD
         GOTO 1000
      ELSE IF(CHOICE.EQ.'Edit Mesh')THEN
         call mesh_edit
         GOTO 1000
      ELSE IF(CHOICE.EQ.'DEFINE OBJECT')THEN
         CALL SETOBJ
         GOTO 1000
      ELSE IF(CHOICE.EQ.'END    ELEMENTS')THEN
C     WHAT ELSE TO DO WHEN 2-D PROBLEM?
         IF(NEL.EQ.0) THEN
            CALL PRS('ERROR: Can''t "END ELEMENTS" when no elements$')
            CALL PRS('are input. ^C if you want to give up$')
         ELSE
            NLEVEL=0
            DO 180 I=1,NEL
               NLEVEL=MAX(NLEVEL,NUMAPT(I))
 180        CONTINUE
            GOTO 320
         ENDIF
      ELSE IF(CHOICE.EQ.'UP   LEVEL') THEN
C     Make Checks
         NTHISL=0
         NNEXTL=0
C     !!?? maybe 96-32??  also, sometimes up gives wrong ele
         MAXLET=0
         DO 190 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL  )NTHISL=NTHISL+1
            IF(NUMAPT(I).EQ.ILEVEL+1)NNEXTL=NNEXTL+1
            IF(NUMAPT(I).EQ.ILEVEL  )MAXLET=
     $           MAX(MAXLET,ICHAR(LETAPT(I)))
 190     CONTINUE
         IF(NTHISL.EQ.0)THEN
            CALL PRSI('*** ERROR *** No elements on level$',ILEVEL)
            CALL PRS('Can''t make level above empty one$')
            GOTO 1000
         ENDIF
         IF(NNEXTL.EQ.0)THEN
            CALL PRS('** New Level **  Enter Height:$')
            CALL PRS('Negative to Abort start of new level$')
            call rer(HEIGHT(ILEVEL+1))
C     !!?? HOW TO MODIFY HEIGHT ONCE IT'S IN?
            IF(HEIGHT(ILEVEL+1) .LE. 0.0) THEN
C     Abort level change
               CALL PRSI('Aborting level change. Still on $',ilevel)
               HEIGHT(ILEVEL+1)=0
               GOTO 1000
            ELSE IF(HEIGHT(ILEVEL+1) .GT. 0.0) THEN
C     Really go up one
               CALL PRSIS('Using Ceiling of Level$',ILEVEL,
     $              '  Mesh as default$')
C     Erase old mesh
               DO 200 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 200           CONTINUE
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
C     Correct Stuff affected if walls are not vertical
                  XCEN(NEL)=0.0
                  YCEN(NEL)=0.0
C     Remove periodic b.c.'s
                  DO 210 IF=1,NFLDS
                     DO 210 ISIDE=1,NSIDES
                        IF(  CBC( ISIDE, NEL,IF).EQ.'P  ')THEN
                           CBC(   ISIDE, NEL,IF)=' '
                           BC (1, ISIDE, NEL,IF)= 0
                           BC (2, ISIDE, NEL,IF)= 0
                           BC (3, ISIDE, NEL,IF)= 0
                        ENDIF
 210                 CONTINUE
                     DO 230 IC=1,4
                        X(NEL,IC)=X(I,IC+4)
                        Y(NEL,IC)=Y(I,IC+4)
                        Z(NEL,IC)=Z(I,IC+4)
                        SIDES(NEL,IC,3)=SIDES(NEL,IC,3)+
     $                       (HEIGHT(ILEVEL-1)+HEIGHT(ILEVEL))/2.
                        XCEN(NEL)=XCEN(NEL)+X(NEL,IC)/4.0
                        YCEN(NEL)=YCEN(NEL)+Y(NEL,IC)/4.0
C     Curved side stuff
                        CCURVE(IC,NEL)=CCURVE(IC+4,I)
                        DO 220 II=1,6
                           CURVE(II,IC,NEL)=CURVE(II,IC+4,I)
 220                    CONTINUE
 230                 CONTINUE
C     RAISE Z OF CEILING OF NEW ELEMENTS
                     DO 240 IC=5,8
                        Z(NEL,IC)=Z(NEL,IC)+HEIGHT(ILEVEL)
 240                 CONTINUE
                     SIDES(NEL,5,3)=SIDES(NEL,5,3)+HEIGHT(ILEVEL-1)
                     SIDES(NEL,6,3)=SIDES(NEL,6,3)+HEIGHT(ILEVEL)
                     CALL DRAWEL(NEL)
                  ENDIF
 250           CONTINUE
C     Sort out the elements so that they will be drawn correctly
               CALL SORTEL
               DO 260 I=1,NEL
                  IF(NUMAPT(ISRT(I)).EQ.ILEVEL)CALL DRAWIS(ISRT(I))
 260           CONTINUE
               GOTO 1000
            ELSE
C     Already have elements on next level. Erase old mesh& draw new
               CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
               DO 270 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 270           CONTINUE
               ILEVEL=ILEVEL+1
               DO 280 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
 280           CONTINUE
               GOTO 1000
            ENDIF
         ELSE IF(CHOICE.EQ.'DOWN LEVEL')THEN
C     Make Checks
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
     $              MAX(MAXLET,ICHAR(LETAPT(I)))
 290        CONTINUE
C     Go down one level.  Erase old mesh& draw new
            CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
            DO 300 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 300        CONTINUE
            ILEVEL=ILEVEL-1
            DO 310 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
 310        CONTINUE
            GOTO 1000
         ELSE
            CALL PRS('CHOICE:'//CHOICE//'NOT IN MENU$')
         ENDIF
         GOTO 1000
 320     CONTINUE
C     Now, cover menu
         CALL SGVIS(13,0)
         CALL SGVIS(13,1)
         NELF= NEL
C     
C     Finished with gridlatch
         IFGRID=.FALSE.
 330     CONTINUE
C     Variable properties put in here
         CALL VPROPS
C     JUMP IN HERE IF YOU DON'T WANT TO MODIFY MESH YOU JUST READ IN
C     Find Sides' Midpoints
         CALL MKSIDE
C     Find Sides' Overlaps (used for finding boundary sides)
C!!   Check here for Irrational B.c.'s: unrequited connectivities;
C     Warn if physical boundaries internally.  This check is not exhaustive!!??
C     ZERO OUT OLD ELEMENTAL CONNECTIVITY
         DO 335 Ifld=1,NFLDS
         DO 335 IEL=1,NEL
         DO 335 ISIDE=1,NSIDES
            IF (CBC( ISIDE, IEL,Ifld).EQ.'E  ') then
c           IF (CBC( ISIDE, IEL,Ifld).EQ.'E  ' .or.
c    $          CBC( ISIDE, IEL,Ifld).EQ.'P  ' .or.
c    $          CBC( ISIDE, IEL,Ifld).EQ.'SP ' .or.
c    $          CBC( ISIDE, IEL,Ifld).EQ.'J  ') then
               CBC(   ISIDE, IEL,Ifld)=' '
               BC (1, ISIDE, IEL,Ifld)= 0
               BC (2, ISIDE, IEL,Ifld)= 0
               BC (3, ISIDE, IEL,Ifld)= 0
            ENDIF
 335     CONTINUE
C     
C     Make internal side comparisons
C     


         icount=1
         do jcount=1,8
            if (icount*100.gt.nel) goto 345
            icount = icount*10
         enddo
  345    continue

         call gencen
         DO 370 Ifld=1,NFLDS
            DO 365 IEL=1,NEL
               if (nel.gt.100.and.mod(iel,icount).eq.0)
     $            write(6,*) 'checking el:',iel
C     
                  IF (MASKEL(IEL,Ifld).EQ.0) GOTO 365
                  DO 360 JEL=1,IEL-1
                     IF (MASKEL(JEL,Ifld).EQ.0) GOTO 360
C     
                     D2 = SQRT( ( xcen(iel)-xcen(jel) )**2 +
     $                  ( ycen(iel)-ycen(jel) )**2 +
     $                  ( zcen(iel)-zcen(jel) )**2 )
                     IF (D2.GT.(RCEN(IEL)+RCEN(JEL))) GOTO 360
C     
                     DELTA = .001 * MIN(rcen(iel),rcen(jel))
C     
                     DO 350 ISIDE=1,NSIDES
                        IF (CBC(Iside,Iel,Ifld).eq.'E') GOTO 350
                        IF (CBC(Iside,Iel,Ifld).eq.'SP') GOTO 350
                        IF (CBC(Iside,Iel,Ifld).eq.'J') GOTO 350
                        INTERN=0
C     
                        DO 340 JSIDE=1,NSIDES
                           IF (CBC(Jside,Jel,Ifld).eq.'E') GOTO 340
                           DELTAX = ABS(SIDES(IEL,ISIDE,1)-
     $                             SIDES(JEL,JSIDE,1))
                           DELTAY = ABS(SIDES(IEL,ISIDE,2)-
     $                             SIDES(JEL,JSIDE,2))
                           DELTAZ = ABS(SIDES(IEL,ISIDE,3)-
     $                             SIDES(JEL,JSIDE,3))
                           IF  (DELTAX .LT. DELTA .AND.
     $                          DELTAY .LT. DELTA .AND.
     $                          DELTAZ .LT. DELTA ) THEN
C     BC Array used to define neighboring elements
C     For want of better notation, 'E' means 
C     elemental (internal) bdry
C     1st two reals are element & side of neighbor; 
C     3rd is orientation
C     
C     Re-do internal connectivity
C     "Normal" internal side.
                              CBC3=CBC( ISIDE, IEL, Ifld)
                              IF(CBC3(3:3) .NE. 'I' .AND. CBC3(3:3) 
     $                                .NE. 'i')THEN
C     Overlapping edges not Internal B.C.'s are elemental b.c.'s
                                 CBC(ISIDE,IEL,Ifld)='E'
                                 CBC(JSIDE,JEL,Ifld)='E'
                              ENDIF
                              IORIEN = 1
                              BC (1, ISIDE, IEL, Ifld) = JEL
                              BC (2, ISIDE, IEL, Ifld) = JSIDE
C                             BC (3, ISIDE, IEL, Ifld) = IORIEN
                              BC (1, JSIDE, JEL, Ifld) = IEL
                              BC (2, JSIDE, JEL, Ifld) = ISIDE
C                             BC (3, JSIDE, JEL, Ifld) = IORIEN
                              GOTO 350
                           ENDIF
                           INTERN=1
 340                    CONTINUE
C     If the side is not internal but has an internal B.C., zero it out
                        CBC3 = CBC(   ISIDE, IEL,Ifld)
                        IF (INTERN.EQ.0 .AND.
     $                     (CBC3.EQ.'E' .OR. CBC3(3:3).EQ.'I'.OR. 
     $                      CBC3(3:3).EQ.'i')
     $                          )THEN
                            CBC(   ISIDE, IEL,Ifld)=' '
                            BC (1, ISIDE, IEL,Ifld)= 0
                            BC (2, ISIDE, IEL,Ifld)= 0
C                           BC (3, ISIDE, IEL,Ifld)= 0
                        ENDIF
 350                 CONTINUE
 360              CONTINUE
 365           CONTINUE
 370  CONTINUE
      CALL GINGRD(0.0)
C     TURN OFF GRIDDING
      CALL SCROLL(7)
C     ??!! MUST PUT P IN OTHER PERIODIC B.C. IF NEK2.  MUST ACCOUNT FOR
C     CONDUCTION ELEMENTS-- THEY DON'T NEED B.C'S FOR THEIR FLUID
      NSIDES=4
      IF(IF3D)NSIDES=6
      DO 380 Ifld=1,NFLDS
      DO 380 IEL=1,NEL
      DO 380 ISIDE=1,NSIDES
         IF (CBC(ISIDE,IEL,Ifld).EQ.' ') then
            NEEDBC=1
            write(6,*) 'bc:',iel,iside,enew(iel)
         ENDIF
 380  CONTINUE
c
c     All 'E-E' boundaries are set.  Check for nonconforming (SP-J)
c
      call check_spj
c
      CALL BOUND
c
c     All internal boundaries are set.  Query remaining bcs
c
      CALL BOUND
C     
      IFMVBD=.FALSE.
      DO 381 Ifld=1,NFLDS
      DO 381 IEL=1,NEL
      DO 381 ISIDE=1,NSIDES
                              CHAR1 = CBC(ISIDE,IEL,Ifld)
                              IF(CHAR1.EQ.'M'.OR.CHAR1.EQ.'m') 
     $                             IFMVBD=.TRUE.
 381  CONTINUE
C     Force a dump of mesh location, unless the user specifiaclly turn it
C     off in the output menu
      IF(IFMVBD) IFXYO = .TRUE.
C     
      CALL MESGEN
C     
      return
      end
c-----------------------------------------------------------------------
      subroutine readat
      include 'basics.inc'
C     Paul's stuff
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
 188  CONTINUE
C     Read Passive scalar data
      read(9,*,err=33)NSKIP
      IF(NSKIP.GT.0) THEN
         DO 177 I=1,NSKIP
            read(9,*,err=33)
 177     CONTINUE
      endif
C     
C     READ DUMMY LOGICALS
      READ(9,*,ERR=33)NLOGIC
      DO 189 I=1,NLOGIC
         if(i.eq.1)read(9,*,err=33)iffold
         if(i.eq.2)read(9,*,err=33)ifhold
         if(I.ne.1.and.i.ne.2)READ(9,*,ERR=33)
 189  CONTINUE
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
     $        ('WARNING: PARAMETER SET TO 2D BUT DATA IS FOR 3D$')
         CALL PRS('PLEASE ENTER Z-HEIGHT: $')
         CALL RER( ZHGHT)
      ENDIF
      IF(.NOT.IF3D.AND.NDIM.EQ.3) THEN
         CALL PRS(
     $        'ERROR: PARAMETER SET TO 3-D BUT DATA IS FOR 2-D$')
      ENDIF
      NCOND=0
      nel50=nel/50
      if (nel.lt.100) nel50=500
      DO 98 IEL=1,NEL
         READ(9,'(20X,I4,4X,I3,A1,11x,i5)',ERR=33,END=33)
     $        IDUM,NUMAPT(IEL),LETAPT(IEL),IGROUP(IEL)
         IF(IGROUP(IEL).GT.0) NCOND=NCOND+1
         IF (NEL.LT.100.or.mod(iel,nel50).eq.0) THEN
            WRITE(S,'(A20,I4,A4,I3,A1,A1)',ERR=33)
     $           ' Reading Element',
     $           IDUM,' [',NUMAPT(IEL),LETAPT(IEL),']'
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
 197           CONTINUE
            ENDIF
         ELSE IF(NDIM.EQ.3)THEN
            READ(9,*,ERR=33,END=33)(X(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(Y(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(Z(IEL,IC),IC=1,4)
            READ(9,*,ERR=33,END=33)(X(IEL,IC),IC=5,8)
            READ(9,*,ERR=33,END=33)(Y(IEL,IC),IC=5,8)
            READ(9,*,ERR=33,END=33)(Z(IEL,IC),IC=5,8)
C     For Flat floors only.  On write will flatten floors.
         ENDIF
 98   CONTINUE
      NELF=NEL-NCOND
C     Read curved side data
      READ(9,*,ERR=57,END=57)
      READ(9,*,ERR=57,END=57)NCURVE
      IF(NCURVE.GT.0)THEN
         DO 19 I=1,NCURVE
            if (nel.lt.1000) then
               READ(9,'(I3,I3,5G14.6,1X,A1)',ERR=57,END=57)
     $              IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            else
               READ(9,'(I2,I6,5G14.6,1X,A1)',ERR=57,END=57)
     $              IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            endif
            CALL DDUMMY(IEDGE,IEL)
            CURVE(1,IEDGE,IEL)=R1
            CURVE(2,IEDGE,IEL)=R2
            CURVE(3,IEDGE,IEL)=R3
            CURVE(4,IEDGE,IEL)=R4
            CURVE(5,IEDGE,IEL)=R5
            CCURVE( IEDGE,IEL)=ANS
 19      CONTINUE
      ENDIF
C     Check here if the old data has the same fields as the new stuff.
C     if not, skip the b.c.'s, etc to avoid confusion
c...  skip, pff        IF( (IFFLOW.AND..NOT.IFFOLD).OR.(.NOT.IFFLOW.AND.IFFOLD)
c...  skip, pff       $.OR.(IFHEAT.AND..NOT.IFHOLD).OR.(.NOT.IFHEAT.AND.IFHOLD)
c...  skip, pff       $.OR. NPSOLD .NE. NPSCAL .or. NKTOLD .ne. NKTONV) THEN
c...  skip, pff           CALL PRS
c...  skip, pff       $ ('Old B.C.''s and options were for different equations.$')
c...  skip, pff           CALL PRS
c...  skip, pff       $ ('Using mesh data only; ignoring old B.C''s and options$')
c...  skip, pff           CALL PRS('Difference, New, Old:$')
c...  skip, pff           IF(NPSOLD.NE.NPSCAL) THEN
c...  skip, pff               CALL PRS('Passive Scalars $')
c...  skip, pff               CALL PRII(npsold,NPSCAL)
c...  skip, pff           ENDIF
c...  skip, pff           IF(NKTOLD.NE.NKTONV) THEN
c...  skip, pff              CALL PRS('Nekton Version $')
c...  skip, pff              CALL PRII(nktold,nktonv)
c...  skip, pff           ENDIF
c...  skip, pff           return
c...  skip, pff        ENDIF
C     Read Boundary Conditions (and connectivity data)
      if (iffmtin) READ(9,*,ERR=44,END=44)
      DO 90 IFLD=1,NoLDS
C     Fluid and/or thermal
C     !!?? NELF DIFFERENT FROM NEL??
         if (iffmtin) READ(9,*,ERR=44,END=44)
         IF( (IFLD.EQ.1.AND.IFFOLD) .OR. IFLD.GT.1)THEN
C     FIX UP FOR WHICH OF FIELDS TO BE USED
            DO 88 IEL=1,NEL
               DO 88 ISIDE=1,NSIDES
C     !Fix to a4,i2 when you make cbc character*4
                  IF(VNEKOLD .LE. 2.5) NBCREA = 3
                  IF(VNEKOLD .GE. 2.6) NBCREA = 5
                  IF (NEL.LT.1000.and.iffmtin) THEN
                     READ(9,'(1X,A3,1x,I2,I3,5G14.6)',ERR=44,END=44)
     $                    CBC(ISIDE,IEL,IFLD),ID,ID,
     $                    (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
                  ELSEIF (iffmtin) then
                     READ(9,'(1X,A3,I5,I1,5G14.6)',ERR=44,END=44)
     $                    CBC(ISIDE,IEL,IFLD),ID,ID,
     $                    (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
                  ELSE
                     READ(8,ERR=44,END=44) chtmp3,
     $                    CBC(ISIDE,IEL,IFLD),ID,ID,
     $                    (BC(II,ISIDE,IEL,IFLD),II=1,NBCREA)
                  ENDIF
                  CBC1=CBC(ISIDE,IEL,IFLD)
                  IF(CBC1.EQ.'M'.OR.CBC1.EQ.'m')THEN
                     IFMOVB=.TRUE.
                     IFGNGM=.TRUE.
                  ENDIF
                  ICBC1=ICHAR(CBC1)
                  IF(ICBC1.GE.97 .AND. ICBC1.LE.122)THEN
                     CBC3 = CBC(ISIDE,IEL,IFLD)
c
c                    pff  no nlines!!
                     NLINES=0
                     do k=1,5
                        BC(k,ISIDE,IEL,IFLD) = 0
                     enddo
ccc  c
ccc                       IF(CBC3(3:3) .EQ. 'i')THEN
ccc  C     Special storage for internal boundaries
ccc                          NLINES=BC(4,ISIDE,IEL,IFLD)
ccc                          BC(5,ISIDE,IEL,IFLD)=LOCLIN
ccc                       ELSE
ccc                          NLINES=BC(1,ISIDE,IEL,IFLD)
ccc                          BC(2,ISIDE,IEL,IFLD)=LOCLIN
ccc                       ENDIF
ccc                       DO 86 I=1,NLINES
ccc                          READ(9,'(A70)',ERR=44,END=44)INBC(LOCLIN)
ccc                          LOCLIN=LOCLIN+1
ccc   86                  CONTINUE
                  ENDIF
 88            CONTINUE
c              DO 89 IEL=1,NEL                       ! NO MORE !
c                 IF(IFMOVB)ICRV(IEL)=1+4+9+16       ! pff, Aug.15,2009
c                 DO 89 IEDGE=1,8
c                    IF(IFMOVB)CCURVE(IEDGE,IEL)='M'
c89               CONTINUE
               ENDIF
 90         CONTINUE
            IF(NFLDS.EQ.1.and.iffmtin)READ(9,*,err=45,end=45)
c
C           Read Initial Conditions
            IF(VNEKOLD.LT.2.5)THEN
               CALL PRS('***** OLD VERSION RESTART DATA ****$')
               CALL PRS(
     $              'Please Re-enter IC''s from the options menu$')
               READ(9,*,ERR=50,END=50)
               DO 21 I=1,7
 21               READ(9,'(A80)',ERR=50,END=50)
               ELSE
C     Check for 2.6 style read file - separate restart/fortran sections.
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
 201                 CONTINUE
                     READ(9,'(A70)',ERR=50,END=50) LINE
                  ENDIF
C     Read fortran initial condition data.
                  REWIND(13)
                  WRITE(13,'(A70)')LINE
                  REWIND(13)
                  READ(13,*,ERR=50,END=50)NLINF
                  REWIND(13)
                  DO 22 I=1,NLINF
                     READ(9,'(A80)')INITC(I)
 22               CONTINUE
               ENDIF
C     Read drive force data
               READ(9,*,ERR=50,END=50)
               READ(9,*,ERR=50,END=50)nskip
               do 23 i=1,nskip
 23               READ(9,'(A80)',ERR=50,END=50)DRIVC(I)
C     READ CONDUCTION ELEMENT DATA
                  READ(9,*,ERR=60,END=60)
                  NFLDSC=0
                  DO 117 IFLD=1,NFLDS
                     IF(IFTMSH(IFLD))NFLDSC=NFLDSC+1
 117              CONTINUE
                  READ(9,*,ERR=60,END=60)NSKIP
                  READ(9,*,ERR=60,END=60)NPACKS
                  IF(NPACKS.GT.0)THEN
                     DO 218 IIG=1,NPACKS
                        READ(9,*)IGRP,IF,ITYPE
                        MATYPE(IGRP,IF)=ITYPE
                        DO 218 IPROP=1,3
                           IF(ITYPE.EQ.1) READ(9,*      ) 
     $                          CPROP(IGRP,IF,IPROP)
                           IF(ITYPE.EQ.2) READ(9,'(A80)') 
     $                          VPROP(IGRP,IF,IPROP)
 218                    CONTINUE
                     ENDIF
C     
C     Read history data
                     READ(9,*,ERR=50,END=50)
                     READ(9,*,ERR=50,END=50)NHIS
C     HCODE(11) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
                     IF(NHIS.GT.0)THEN
                        DO 51 I=1,NHIS
                           READ(9,'(1X,11A1,1X,4I5)',ERR=50,END=50)
     $                          (HCODE(II,I),II=1,11),
     $                          (LOCHIS(II,I),II=1,4)
 51                     CONTINUE
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
                           READ(9,'(1X,L1,1X,A5)',ERR=50,END=50) 
     $                          IFPSCO(I),PSNAME(I)
 1221                   CONTINUE
                     ENDIF
C     OBJECT SPECIFICATION DATA
                     READ(9,*,ERR=50,END=50)
                     READ(9,*,ERR=50,END=50)NSOBJS
                     IF(NSOBJS .GT. 0)THEN
                        DO 62 IOBJ=1,NSOBJS
                           READ(9,'(4X,I4,35X,A20)',ERR=50,END=50)
     $                          NFACE(IOBJ),SOBJ(IOBJ)
                           DO 63 IFACE=1,NFACE(IOBJ)
                              READ(9,*,ERR=50,END=50)
     $                             ILSURF(1,IOBJ,IFACE),
     $                             ILSURF(2,IOBJ,IFACE)
 63                        CONTINUE
 62                     CONTINUE
                     ENDIF
                     READ(9,*,ERR=50,END=50)NVOBJS
                     READ(9,*,ERR=50,END=50)NEOBJS
                     READ(9,*,ERR=50,END=50)NPOBJS
C     
                     return
 33                  CALL PRS
     $                    ('ERR: MESH DATA MISSING OR CORRUPTED.$')
                     CALL BEEP
                     CALL PRS('<CR> TO CONTINUE$')
                     CALL RES(LINE,0)
                     CALL DEVEX
                     CALL EXITT
 34                  CALL PRS('ERR: FILE CONTAINS ONLY PARAM DATA.$')
                     CALL PRS('ERR: MESH DATA MISSING$')
                     CALL BEEP
                     CALL PRS('<CR> TO CONTINUE$')
                     CALL RES(LINE,0)
                     CALL DEVEX
                     CALL EXITT
 35                  CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $                    '  CONTINUEING WITH MESH ONLY.$')
                     return
 44                  WRITE(S,'(1X,A36,I4,A6,I3)')
     $                    'Err reading B.C. data for elm',IEL,
     $                    ' side ',ISIDE
                     CALL PRS(S//'$')
                     return
 45                  CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $                    '  CONTINUEING WITH MESH ONLY.$')
                     return
 40                  CALL PRS(
     $                    'CURVED SIDES NOT READ CORRECTLY. PLS CHECK$')
                     CALL PRS('SIDES IN CURVED SIDE MENU$')
 50                  CALL PRS
     $                    ('OPTIONS MAY NOT HAVE BEEN READ CORRECTLY.$')
                     CALL PRS('PLS CHECK THEM IN OPTIONS MENU$')
                     return
 57                  CALL PRS
     $                    ('CURVED SIDES HAVE NOT BEEN READ CORRECTLY$')
                     CALL PRS('PLS CHECK THEM.$')
                     return
 60                  CALL PRS
     $                    ('CONDUCTION ELEMENTS MAY NOT HAVE BEEN$')
                     CALL PRS('READ CORRECTLY$')
                     CALL PRS('PLEASE CHECK THEM IN OPTIONS MENU$')
 70                  CALL PRS('NO OBJECTS$')
      return
      END
c-----------------------------------------------------------------------
      subroutine setscl
C     Sets scale factors
C     
      include 'basics.inc'
      LOGICAL IFMENU
      REAL XSC(4),YSC(4)
      GRID=PARAM(18)
 3    CONTINUE
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
     $     'Scale Factors: Input thru Mouse, type, or use normalized?$')
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
     $        'They are the lengths that correspond to 1 scr height$')
         CALL RERR(XFAC,YFAC)
         IF(IFLEARN)WRITE(3,*)XFAC,YFAC
         XFAC=ABS(XFAC)
         YFAC=ABS(YFAC)
C     
         CALL PRS
     $        ('Enter the coordinates of Left Bottom corner of screen$')
         CALL RERR(XZERO,YZERO)
         IF(IFLEARN)WRITE(3,*)XZERO,YZERO
C     Mark Zero
         CALL COLOR(7)
      ELSE IF(CHOICE.EQ.'MOUSE') THEN
 17      CALL PRS('Push left button with mouse at 2 different'//
     $        ' x (horizontal)$')
         CALL PRS('locations$')
         DO 18 I=1,3
            IF(I.EQ.2)CALL PRS('2nd x location$')
            IF(I.EQ.4)CALL PRS('2nd y location$')
            IF(I.EQ.3)THEN
               CALL PRS('Enter corresponding x length$')
               CALL GINGRD(0.0)
               call rer(XLEN)
               CALL GINGRD(GRID)
c              CALL PRS('Push left button with mouse at 2'//
c    $              ' different y (vertical) locations$')
c              CALL PRS('or click in menu area for equal X-Y scaling.$')
            ENDIF
            if (i.lt.3) CALL MOUSE(XSC(I),YSC(I),BUTTON)
C     
C     Check for keypad entry?
C     
c           IF (IFMENU(XSC(I),YSC(I)).AND.I.LE.2) THEN
c              CALL PRS('Please enter two points in build area.$')
c              GOTO 17
c           ENDIF
c           IF (IFMENU(XSC(I),YSC(I)).AND.I.GT.2) THEN
            IF (I.GT.2) THEN
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
 18      CONTINUE
         CALL PRS('Enter in corresponding y length$')
         CALL GINGRD(0.0)
         call rer(YLEN)
         CALL GINGRD(GRID)
 19      CONTINUE
         IF(XSC(2).EQ.XSC(1) .OR. YSC(4).EQ.YSC(3))THEN
            CALL PRS('ERROR: MUST BE NONZERO GAP.  START OVER$')
            GO TO 3
         ENDIF
C     
C     Compute and check scale factors
C     
         XFACT=ABS(XLEN/(XSC(2)-XSC(1)))
         YFACT=ABS(YLEN/(YSC(4)-YSC(3)))
         IF(XFACT.EQ.0.0 .OR. YFACT.EQ.0.0)THEN
            CALL PRS('ERROR: MUST BE NONZERO SCALE FACTORS '//
     $           'START OVER$')
            GO TO 3
         ENDIF
C     
         DO 20 I=1,ICOUNT
C     Erase markers
            CALL COLOR(0)
            CALL MOVEC(XSC(I)-.02,YSC(I))
            CALL DRAWC(XSC(I)+.02,YSC(I))
            CALL MOVEC(XSC(I),YSC(I)-.02)
            CALL DRAWC(XSC(I),YSC(I)+.02)
 20      CONTINUE
C     
C===============================================================
C     Set Origin
C===============================================================
C     
         CALL PRS (' $')
         CALL PRS
     $        ('Push left button with mouse at known point in x,y$')
         CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
C     
         CALL PRS('Enter X coordinate of point:$')
         CALL GINGRD(0.0)
         call rer(XPHYS)
         CALL PRS('Enter Y coordinate of point:$')
         call rer(YPHYS)
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
      return
      end

c-----------------------------------------------------------------------
      subroutine imp_mesh(ifquery_displace)

      include 'basics.inc'
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)

      logical ifquery_displace,ifdisplace

      real    xyzbox(6)

      call prs('input name of new .rea file$')
      call blank(fname,70)
      call res  (fname,70)

      if (indx1(fname,'base',4).eq.1) then
         call imp_mesh_special
         return
      endif

      ifdisplace = .false.
      if (ifquery_displace) then
       call prs('Would you like to displace existing elements in box?$')
       call res  (ans,1)
       if (ans.eq.'y'.or.ans.eq.'Y') ifdisplace = .true.
      endif

      if (indx1(fname,'.rea',4).eq.0) then !  Append .rea, if not present
         len = ltrunc(fname,70)
         call chcopy (fnam1(len+1),'.rea',4)
      endif
      open(unit=47,file=fname,status='old',err=1000)

c     Successfully opened file, scan until 'MESH DATA' is found
      call readscan('MESH DATA',9,47)
      read (47,*) nelin,ndimn
      neln = nelin + nel

      if (ndimn.ne.ndim) then
         call prs('Dimension must match dimension of current session.$')
         call prs('Returning.$')
         close(47)
         return
      endif

      if (neln.ge.nelm) then
         call prs('Sorry, number of elements would exceed nelm.$')
         call prs('Returning.$')
         call prii(neln,nelm)
         close(47)
         return
      endif

c     Read geometry

      nelo = nel
      nels = nel+1

      ierr=imp_geom(x(nels,1),y(nels,1),z(nels,1),nelm
     $             ,numapt(nels),letapt(nels),igroup(nels)
     $             ,ndim,nelin,47)
      if (ierr.ne.0) then
         call prs('Error reading geometry... returning.$')
         close(47)
         return
      endif
c
      ierr=imp_curv(nc,ccurve,curve,ndim,nelin,nel,47)
      if (ierr.ne.0) then
         call prs('Error reading curve side info... returning.$')
         close(47)
         return
      endif
c
c
c     Read BC info
c
c
      read(47,3) d
    3 format(a3)
c
      ifld0 = 1
      if (.not.ifflow) ifld0=2
      write(6,*) 'IFLD:',ifld0,nflds,ifflow,ifheat
c
      do ifld=ifld0,nflds
         if (.not.ifflow) read(47,*) ans  ! dummy read
         ierr
     $   =imp_bc(cbc(1,nels,ifld),bc(1,1,nels,ifld),ndim,nelin,nel,47)
         if (ierr.ne.0) then
            call prsii('nelin,ifld:$',nelin,ifld)
            call prs('Error reading boundary conditions. Returning.$')
            close(47)
            return
         endif
      enddo
c
c     ALL DONE
c
      close(47)
c
      write(6,*) 'This is nel,ncurve old:',nel,ncurve
      nel = neln
      ncurve = ncurve + nc
      write(6,*) 'This is nel,ncurve new:',nel,ncurve
c
      call gencen
c
c
      if (ifdisplace) then
c
         ncrnr = 2**ndim
c
         xyzbox(1) = 9.e21
         xyzbox(2) =-9.e21
         xyzbox(3) = 9.e21
         xyzbox(4) =-9.e21
         xyzbox(5) = 9.e21
         xyzbox(6) =-9.e21
c
         do ie=nels,nel   ! Find box containing imported mesh
         do i=1,ncrnr
            xyzbox(1) = min(xyzbox(1),x(ie,i))
            xyzbox(2) = max(xyzbox(2),x(ie,i))
            xyzbox(3) = min(xyzbox(3),y(ie,i))
            xyzbox(4) = max(xyzbox(4),y(ie,i))
            xyzbox(5) = min(xyzbox(5),z(ie,i))
            xyzbox(6) = max(xyzbox(6),z(ie,i))
         enddo
         enddo
c
         write(s,110) 'X',xyzbox(1),xyzbox(2)
         call prs(s)
         write(s,110) 'Y',xyzbox(3),xyzbox(4)
         call prs(s)
         if (if3d) write(s,110) 'Z',xyzbox(5),xyzbox(6)
         if (if3d) call prs(s)
  110    format(' ',a1,' min/max of new elements:',2g15.7,'$')
c
         call substitute_el(xyzbox,nelo)
      endif
c
      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      close(47)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine readscan(key,len_key,io)
c
c     Read a file until "key" is found or eof is found.
c
c
      character*80 key
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
   80 format(a80)
   81 format(80a1)
c
      do i=1,90000
         call blank(string,80)
         read (io,80,end=100,err=100) string
         len = ltrunc(string,80)
         if (indx1(string,key,len_key).ne.0) return
      enddo
  100 continue
      return
      end
c-----------------------------------------------------------------------
      function imp_geom(x,y,z,ld,numapt,letapt,igroup,ndim,nel,io)
      real        x(ld,8),y(ld,8),z(ld,8)
      integer     numapt(1),igroup(1)
      character*1 letapt(1)
c
      character*1 c
c
c
c
      if (ndim.eq.3) then
         do ie=1,nel
c
c           read(io,1) c
c   1       format(a1)
c
c           read(io,'(20x,i4,4x,i3,a1,11x,i5)',err=9,end=9)
c    $      idum,numapt(ie),letapt(ie),igroup(ie)
            read(io,*)
            numapt(ie) = 1
            letapt(ie) = 'A'
            igroup(ie) = 0
c
            read(io,*,end=9,err=9) (x(ie,k),k=1,4)
            read(io,*,end=9,err=9) (y(ie,k),k=1,4)
            read(io,*,end=9,err=9) (z(ie,k),k=1,4)
c
            read(io,*,end=9,err=9) (x(ie,k),k=5,8)
            read(io,*,end=9,err=9) (y(ie,k),k=5,8)
            read(io,*,end=9,err=9) (z(ie,k),k=5,8)
c
         enddo
      else
         do ie=1,nel
c
c           read(io,1) c
c
c           read(io,'(20x,i4,4x,i3,a1,11x,i5)',err=9,end=9)
c    $      idum,numapt(ie),letapt(ie),igroup(ie)
            read(io,*)
            numapt(ie) = 1
            letapt(ie) = 'A'
            igroup(ie) = 0
c
            read(io,*,end=9,err=9) (x(ie,k),k=1,4)
            read(io,*,end=9,err=9) (y(ie,k),k=1,4)
         enddo
      endif
c
      imp_geom = 0
      return
c
    9 continue
      imp_geom = 1
      return
      end
c-----------------------------------------------------------------------
      function imp_curv(nc,cc,c,ndim,nel,nelo,io)
      character*1 cc(12,1)
      real        c(6,12,1)
c
      character*1 d
      logical iffmtin
c
c
c
      iffmtin = .true.
c
      read(io,1) d
    1 format(a1)
c
      read(io,*) nc
      write(6,*) 'Found',nc,' curve sides.'
c
      do i=1,nc
         if (nel.lt.1000.AND.iffmtin) then
            read(io,'(2i3,5g14.6,1x,a1)',err=9,end=9)
     $                    iedge,ie,r1,r2,r3,r4,r5,d
         elseif (iffmtin) THEN
            read(io,'(i2,i6,5g14.6,1x,a1)',err=9,end=9)
     $                    iedge,ie,r1,r2,r3,r4,r5,d
         else
            read(io,err=9,end=9) iedge,ie,r1,r2,r3,r4,r5,d
         endif
         ie = nelo + ie
         c(1,iedge,ie)=r1
         c(2,iedge,ie)=r2
         c(3,iedge,ie)=r3
         c(4,iedge,ie)=r4
         c(5,iedge,ie)=r5
         c(6,iedge,ie)=0.
         cc( iedge,ie)=d
      enddo
c
      imp_curv = 0
      return
c
    9 continue
      imp_curv = 1
      return
      end
c-----------------------------------------------------------------------
      function imp_bc(cbc,bc,ndim,neln,nel,io)
      character*3  cbc   (6,1)
      real          bc (5,6,1)
c
      character*80 a80
      character*1  a81(1)
      equivalence(a81,a80)
c
      imp_bc = 0
c
      call blank (a80,80)
      read(io,80) a80
   80 format(a80)
c
      call capit(a80,80)
      len = ltrunc(a80,80)
      a81(len+1) = '$'
      call prs(a80)
      write(6,80) a80
c
      i1 = indx1(a80,'NO',2)
      write(6,*) i1,' indx1 '
      if (i1.ne.0) return
c     if (indx1(a80,'NO',2).ne.0) return
c
      nsides = 2*ndim
      nbcrea = 5
c
      do ie=1,neln
      do iside = 1,nsides
         if (neln.lt.1000) then
            read(io,'(1x,a3,1X,i2,i3,5g14.6)',err=9,end=9)
     $      cbc(iside,ie),id,jd,
     $      (bc(ii,iside,ie),ii=1,nbcrea)
         else
            read(io,'(1x,a3,i5,i1,5g14.6)',err=9,end=9)
     $      cbc(iside,ie),id,jd,
     $      (bc(ii,iside,ie),ii=1,nbcrea)
         endif
c        Adjust periodic boundary conditions
         if (cbc(iside,ie).eq.'P  ') bc(1,iside,ie) = bc(1,iside,ie)+nel
      enddo
      enddo
c
      imp_bc = 0
      return
c
    9 continue
      write(6,*) 'made it to:',ie,iside
      imp_bc = 1
      return
      end
c
c-----------------------------------------------------------------------
      subroutine new_out(fname,lname,gv,ncrnr,xyz,ind)
c
c     Read a file until "key" is found or eof is found.
c
      include 'basics.inc'

      parameter (maxv = 8*nelm)
      real      xp(maxv),yp(maxv),zp(maxv)
      common /c_xyz/ xp,yp,zp
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      integer gv(ncrnr,1),ind(1)
      real    xyz(1)
c
      integer ecrnr(8)
      save    ecrnr
      data    ecrnr  / 1,2,4,3,5,6,8,7 /
c
   80 format(a80)
   81 format(80a1)
c
      call blank(string,80)
      call chcopy(string,fname,lname)
      call chcopy(string1(lname+1),'.map',4)
c
c     Get current vertex map info
c
      open(unit=10,file=string,err=999)
      read(10,*) neln
      do ie=1,neln
         read(10,*,end=998,err=998) idum,(gv(k,ie),k=1,ncrnr)
      enddo
      close(unit=10)
c
c
c     Sort data and gridpoints by global vertex number
c
      l = 1
      do ie = 1,neln
         do k=1,ncrnr
            j     = ecrnr(k)
            xp(l) = x (ie,j)
            yp(l) = y (ie,j)
            zp(l) = z (ie,j)
            l=l+1
         enddo
      enddo
c
      nv = ncrnr*neln
      call irank(gv,ind,nv)
      call swap (xp,xyz,ind,nv)
      call swap (yp,xyz,ind,nv)
      call swap (zp,xyz,ind,nv)
c
c     Compress lists
c
      nnv=1
      do i=2,nv
         if (gv(i,1).ne.gv(nnv,1)) then
            nnv       = nnv+1
            gv(nnv,1) = gv(i,1)
            xp(nnv)   = xp(i)
            yp(nnv)   = yp(i)
            zp(nnv)   = zp(i)
         endif
      enddo
      nv = nnv
c
c     Output vertex info
c
      l=1
      if (ndim.eq.3) then
         do i=1,nv
            xyz(l  ) = xp(i)
            xyz(l+1) = yp(i)
            xyz(l+2) = zp(i)
            l=l+3
         enddo
      else
         do i=1,nv
            xyz(l  ) = xp(i)
            xyz(l+1) = yp(i)
            l=l+2
         enddo
      endif
      npts_o = nv*ndim
c
      call chcopy(string1(lname+1),'.vtx',4)
      open(unit=11,file=string,err=999)
      write(11,*) npts_o
      if (ndim.eq.3) then
         write(11,33) (xyz(i),i=1,npts_o)
      else
         write(11,32) (xyz(i),i=1,npts_o)
      endif
   32 format(1p2e15.7)
   33 format(1p3e15.7)
      return
c
  998 continue
      call prs('end of .map file reached')
      return
c
  999 continue
      call prs('Could not open .map file.')
      return
      end
c-----------------------------------------------------------------------
      function i_finds(list,ig,n)
c     Find if there is a match of ig in *sorted* list
      integer i_finds,ig,list(1),n
      integer hi,lo,m
c
      lo = 1
      hi = n
      icount  = 0
      i_finds = 0
c
    1 continue 
      m = (hi+lo)/2
      if (ig.eq.list(m)) then
         i_finds = m
         return
      elseif (ig.lt.list(m)) then
         hi = m
      else
         lo = max(m,lo+1)
      endif
      icount = icount+1
      if (icount.gt.n.or.lo.gt.hi) then
         write(6,*) 'ERROR IN i_finds:',ig,n,hi,lo,m,' ABORT'
         j = i_findu(list,ig,n)
         i_finds = j
         write(6,*) 'result of i_findu:',j,list(j)
         call exitt
         return
      else
         goto 1
      endif
      return
      end
c-----------------------------------------------------------------------
      function i_findu(list,ig,n)
c     Find if there is a match of ig in *unsorted* list
      integer i_findu,ig,list(1),n
      integer hi,lo,m
c
      i_findu = 1
c
      do m=1,n
         if (ig.eq.list(m)) then
            i_findu = m
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_vtx
c
c     Read vtk-like unstructured mesh format
c
      include 'basics.inc'
c
      parameter (maxv = 8*nelm)
      integer        vv(8,nelm),vnum(maxv)
      common /c_vtx/ vv,vnum

      real      xp(maxv),yp(maxv),zp(maxv)
      common /c_xyz/ xp,yp,zp

      integer      kcell(8),indx(8)
      save         indx
      data         indx  / 1,2,4,3,5,6,8,7 /
C
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)
c
      call prsi('Current number of elements:$',nel)
      call prs('input name of new vertex/cell file$')
      call blank(fname,70)
      call res  (fname,70)
      ikill = indx1(fname,' ',1)
      nkill = 70-ikill+1
      call blank(fnam1(ikill),nkill)
c
      open(unit=47,file=fname,status='old',err=1000)
c
c     Read vertex info
c
      read(47,*) nvtx 
c
      if (nvtx.gt.maxv) then
         nvtx = nvtx+1
         call prsi('Too many pts. Inc. maxv to:$',nvtx)
         close(47)
         return
      endif
c
      if (if3d) then
        do i=1,nvtx
c        read(47,*) vnum(i),xp(i),yp(i),zp(i)
         read(47,*) xp(i),yp(i),zp(i)
         vnum(i) = i
        enddo
      else
        do i=1,nvtx
c        read(47,*) vnum(i),xp(i),yp(i)
         read(47,*) xp(i),yp(i)
         vnum(i) = i
        enddo
      endif
      call prsis('Found$',nvtx,' points.$')
      xmin = glmin(xp,nvtx)
      xmax = glmax(xp,nvtx)
      ymin = glmin(yp,nvtx)
      ymax = glmax(yp,nvtx)
      call prsrr('xmin xmax:',xmin,xmax)
      call prsrr('ymin ymax:',ymin,ymax)
c
c
c     Read cell data
c
      read(47,*) ncell
c
      if (ncell+nel.gt.nelm) then
         write(6,*) ncell,nel,nelm,' ncell,nel,nelm'
         ncell = ncell+nel
         call prsi('Too many elements. Increase nelm to:$',ncell)
         close(47)
         return
      endif
c
      npt=2**ndim
      do ie=1,ncell
c        read(47,*) ii1,ii2,ii3,(kcell(k),k=1,npt)
         read(47,*) (kcell(k),k=1,npt)
         do k=1,npt
c           kcell(k) = kcell(k) + 1   ! no need to add 1
            vv(k,ie) = i_findu(vnum,kcell(k),nvtx)
c           write(6,*) ie,k,vv(k,ie),kcell(k),nvtx
c           vv(k,ie) = i_finds(vnum,kcell(k),nvtx)
         enddo
         NUMAPT(IE)=1
         LETAPT(IE)='A'
      enddo
      call prs('continue?$')
      call res(ans,1)
c
c     At this point, we have the cell data in the "new" format.
c
c     Now we chuck it in favor of the worthless nekton format.
c
c     4/17/01 -- first, try to figure out the local vertex ordering
c
c
      call prs('Flip for Gambit? (Yes=Gambit, No=Seung)$')
      call res(ans,1)
c
c
      if (ans.eq.'n'.or.ans.eq.'N') then
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               x(nel,k)=xp(l)
               y(nel,k)=yp(l)
               z(nel,k)=zp(l)
            enddo
         enddo
      else                     ! flip
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               j = indx(k)
               x(nel,j)=xp(l)
               y(nel,j)=yp(l)
               z(nel,j)=zp(l)
            enddo
         enddo
      endif
c        do ie=1,ncell
c           do k=1,npt
c              write(6,*) ie,k,x(ie,k),y(ie,k)
c           enddo
c        enddo
      if_unstructure=.true.
      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      call prs(fname)
      close(47)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_vtk
c
c     Read vtk-like unstructured mesh format
c
      include 'basics.inc'
c
      parameter (maxv = 8*nelm)
      integer        vv(8,nelm),vnum(maxv)
      common /c_vtx/ vv,vnum

      real      xp(maxv),yp(maxv),zp(maxv)
      common /c_xyz/ xp,yp,zp

      integer      kcell(8),indx(8)
      save         indx
      data         indx  / 1,2,4,3,5,6,8,7 /
C
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)
c
      character*1  s81(80)
      character*80 s80
      equivalence (s80,s81)
c
      call prsi('Current number of elements:$',nel)
      call prs('input name of new vertex/cell file$')
      call blank(fname,70)
      call res  (fname,70)
      ikill = indx1(fname,' ',1)
      nkill = 70-ikill+1
      call blank(fnam1(ikill),nkill)
c
      open(unit=47,file=fname,status='old',err=1000)
c
c     Successfully opened file, scan until 'MESH DATA' is found
c
c     Read vertex info
c
c     read(47,*) 
c     read(47,*) 
c     read(47,*) 
c
      call blank(s80,80)
      read(47,80) s80
      ib = indx1(s80,' ',1)
      open (unit=48,file='my_vtk.tmp')
      write(48,81) (s81(k),k=ib+1,80)
      rewind (unit=48)
      read (48,*) nvtx
      close (unit=48)
c
      if (nvtx.gt.maxv) then
         nvtx = nvtx+1
         call prsi('Too many pts. Inc. maxv to:$',nvtx)
         close(47)
         return
      endif
c
      if (if3d) then
        do i=1,nvtx
c        read(47,*) xp(i),yp(i),zp(i)
c        vnum(i) = i-1
         read(47,*) vnum(i),xp(i),yp(i),zp(i)
         write(6,*) vnum(i),xp(i),yp(i),zp(i),i
        enddo
      else
        do i=1,nvtx
         read(47,*) xp(i),yp(i)
         vnum(i) = i-1
        enddo
      endif
c
c
c     Read cell data
c
      call blank(s80,80)
      read(47,80) s80
      write(6,*) s80
      ib = indx1(s80,' ',1)
      open (unit=48,file='my_vtk.tmp')
      write(48,81) (s81(k),k=ib+1,80)
      rewind (unit=48)
      read (48,*) ncell
      close (unit=48)
c
c
      if (ncell+nel.gt.nelm) then
         ncell = ncell+nel
         call prsi('Too many elements. Increase nelm to:$',ncell)
         close(47)
         return
      endif
c
      nelo = 0
      npt=2**ndim
      do ie=1,ncell
c        read(47,*) nv_per_cell,(kcell(k),k=1,npt)
         write(6,*) ie,npt,ncell,'  KCELL'
         read(47,*) (kcell(k),k=1,npt)
         write (6,*) (kcell(k),k=1,npt)
         do k=1,npt
            vv(k,ie) = i_findu(vnum,kcell(k),nvtx)
c           vv(k,ie) = i_finds(vnum,kcell(k),nvtx)
         enddo
      enddo
c
c     At this point, we have the cell data in the "new" format.
c
c     Now we chuck it in favor of the worthless nekton format.
c
c     4/17/01 -- first, try to figure out the local vertex ordering
c
c
c     call prs('Flip for Gambit? (Yes=Gambit, No=Seung)$')
c     call res(ans,1)
      ans = 'y'
c
c
      if (ans.eq.'n'.or.ans.eq.'N') then
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               x(nel,k)=xp(l)
               y(nel,k)=yp(l)
               z(nel,k)=zp(l)
            enddo
         enddo
      else                     ! flip
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               j = indx(k)
               x(nel,j)=xp(l)
               y(nel,j)=yp(l)
               z(nel,j)=zp(l)
            enddo
         enddo
      endif
      if_unstructure=.true.

      call set_std_bcs(nelo)


      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      call prs(fname)
      close(47)
c
   80 format(a80)
   81 format(80a1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine check_spj
c
c     Check for split (parent) - join (child) nonconforming interfaces.
c
      include 'basics.inc'
      logical ifclose,ifok
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /
c
c     Currently designed only for 2D, NO CURVE nonconf.
c
      if (ndim.eq.3) return
c
      nsides = 2*ndim
      do kf=1,nflds
      do ie=1,nel
       if (maskel(ie,kf).gt.0) then
        do is=1,nsides
c         write(6,*) 'CBC 1:',cbc(is,ie,kf),ie,is
          if (cbc(is,ie,kf).ne.'E  ') then   ! n^2 loop time!
            do je=1,nel
             if (maskel(je,kf).gt.0) then
              do js=1,nsides
               if (cbc(js,je,kf).ne.'E  '.and.ie.ne.je) then
                  if (ifclose(ie,is,je,js)) then
c                    Currently designed only for 2D, NO CURVE nonconf.
                     if (if3d) then
                       call chk_segs(ifok,icond,ie,is,je,js)
                       if (ifok) then
                          cbc(is,ie,kf)='SP '
                          bc (5,is,ie,kf)=bc (5,is,ie,kf)+1
                          cbc(js,je,kf)='J  '
                          bc (1,js,je,kf)=ie
                          bc (2,js,je,kf)=is
c                         write(6,*) 'Match:',ie,is,je,js,icond
c                       write(6,*) 'Match:',cbc(is,ie,kf),cbc(js,je,kf)
                       else
                       endif
                     else
                       call chk_seg(g1,e1,g2,e2,ie,is,je,js)
                       if (abs(e1).le.0.001          .and.
     $                     abs(e2).le.0.001          .and.
     $                     -.01.le. g1.and.g1.le.1.01 .and.
     $                     -.01.le. g2.and.g2.le.1.01) then  ! segment
                          cbc(is,ie,kf)='SP '
                          bc (5,is,ie,kf)=bc (5,is,ie,kf)+1
                          cbc(js,je,kf)='J  '
                          bc (1,js,je,kf)=ie
                          bc (2,js,je,kf)=is
                          bc (3,js,je,kf)=g1
                          bc (4,js,je,kf)=g2
                       endif
                     endif
                  endif
               endif
              enddo
             endif
            enddo
         endif
c        write(6,*) 'CBC 2:',cbc(is,ie,kf),ie,is
        enddo
       endif
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      logical function ifclose(ie,is,je,js)
c
c
      include 'basics.inc'
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /

c
      ifclose = .true.
c
      nv=2**(ndim-1)
      dz=0.
      do jv=1,nv
         dx = sides(ie,is,1)-x(je,jvs(jv,js))
         dy = sides(ie,is,2)-y(je,jvs(jv,js))
         if (if3d) dz = sides(ie,is,3)-z(je,jvs(jv,js))
         ds = dx*dx+dy*dy+dz*dz
         if (ds.le.rcen(ie)*rcen(ie)) return
      enddo
      ifclose = .false.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_seg(g1,e1,g2,e2,ie,is,je,js)
c
c     Upon return,  (g1,e1) and (g2,e2) are the transformed
c                   coordinates of the two vertices associated with
c                   side (je,js) in the coordinate system defined
c                   by the two vertices associated with side (ie,is).
c
c     For example, if side (je,js) is *coincident* with (ie,is), then either
c                  (g1,e1) = (0.,0.) and  (g2,e2) = (1.,0.),  or
c                  (g1,e1) = (1.,0.) and  (g2,e2) = (0.,0.).
c
c
      include 'basics.inc'
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /
c
c
c     find the (x,y) coordinates of the two endpoints of side (ie,is)
      x1 = x(ie,jvs(1,is))
      y1 = y(ie,jvs(1,is))
      x2 = x(ie,jvs(2,is))
      y2 = y(ie,jvs(2,is))
c
c     find the (x,y) coordinates of the two endpoints of side (je,js)
      x3 = x(je,jvs(1,js))
      y3 = y(je,jvs(1,js))
      x4 = x(je,jvs(2,js))
      y4 = y(je,jvs(2,js))
c
c     use these endpoints to set up the transformation.
c     the 3x3 matrix is not explicitly formed, just the closed 
c     form result.
c
      x2_x1 = x2 - x1
      y2_y1 = y2 - y1
      hyp2 = ( x2_x1*x2_x1 + y2_y1*y2_y1 )
c
c
      scale  = 1./hyp2
      g1 = scale*( x3*x2_x1 + y3*y2_y1 - x1*x2_x1 - y1*y2_y1 )
      e1 = scale*(-x3*y2_y1 + y3*x2_x1 + x1*y2_y1 - y1*x2_x1 )
      g2 = scale*( x4*x2_x1 + y4*y2_y1 - x1*x2_x1 - y1*y2_y1 )
      e2 = scale*(-x4*y2_y1 + y4*x2_x1 + x1*y2_y1 - y1*x2_x1 )
c
c     write(8,1) ie,is,je,js,g1,e1,g2,e2
c   1 format(4i4,1p4e12.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_segs(ifok,icond,ie,is,je,js)
c
c     3D routine for determining if face(je,js) is a subset of
c     face(ie,is)
c
c     Currently predicated on assumption that faces are PLANAR.
c
c
      include 'basics.inc'
c
      logical ifok,ifxoeq0
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,6,5 , 2,3,7,6 , 3,4,8,7
     $            , 4,1,5,8 , 4,3,2,1 , 5,6,7,8 /  ! EB cclockws ordering
c
      real n(3),nhat(3),xx(3,4),xj(3),xo(3),d1(3),d2(3),norm
      real vnt(3),vv(3),ee(3),op(3),xjp(3)
c
c     First, get canonical coordinates of plane containing face i
c
c     if (ie.eq.1.and.je.gt.1.and.is.eq.6.and.js.eq.5)
c    $    write(6,*) 'ie,is,je,js:',ie,is,je,js
c     write(6,*) 'Continue? (any key)'
c     read (5,*) ans
c
      do i=1,4
         xx(1,i) = x(ie,jvs(i,is))
         xx(2,i) = y(ie,jvs(i,is))
         xx(3,i) = z(ie,jvs(i,is))
      enddo
c     write(6,1) (i,(xx(k,i),k=1,3),' xx',i=1,4)
    1 format(/,4(i4,3f10.4,a3,/))
c
      do k=1,3
        xo(k) = 0.
        do i=1,4
           xo(k) = xo(k)+xx(k,i)
        enddo
        xo(k) = 0.25*xo(k)
      enddo
c     write(6,1) i,(xo(k),k=1,3),' xo'
c
      call sub3(d1,xx(1,3),xx(1,1),3) ! vectors connecting diagonally
      call sub3(d2,xx(1,2),xx(1,4),3) ! opposed vertices of face i
      call vcross_normal(nhat,d2,d1)
c     write(6,1) i,(nhat(k),k=1,3),' nh'
      d1n = dotprod(d1,d1)
      d2n = dotprod(d2,d2)
      clen = max(d1n,d2n)
      clen = sqrt(clen)
c
c     Generate basis vectors for plane i
c
      call norm3d(d1)
      call vcross_normal(d2,nhat,d1)
c     write(6,1) i,(d1(k),k=1,3),' d1'
c     write(6,1) i,(d2(k),k=1,3),' d2'
c
c     Now for each vertex of face j, check if it lies on the plane of face i
c     and if it's on the edge of or interior to face i.
c
      ifok  = .false.
      icond = 0
      do j=1,4
         xj(1) = x(je,jvs(j,js))
         xj(2) = y(je,jvs(j,js))
         xj(3) = z(je,jvs(j,js))
c        write(6,1) j,(xj(k),k=1,3),' xj'
c
c        map to canonical plane of face i
c
         call sub3(xjp,xj,xo,3)
         xp1 = dotprod(d1  ,xjp)
         xp2 = dotprod(d2  ,xjp)
         xp3 = dotprod(nhat,xjp)
         xjp(1) = xp1
         xjp(2) = xp2
         xjp(3) = xp3
c        write(6,1) j,(xjp(k),k=1,3),'xjp'
c
c        Now check for any disqualifications 
c
         icond = 2
         epsr  = 1.e-4
         rel_err = abs(xjp(3))/clen
c        write(6,*) j,rel_err,xjp(3),clen,'rel_err'
         if (rel_err.gt.epsr) return
c
c        Vertex is co-planar with face i, check if interior to face i
c
         epsn = 1.e-4*clen
         sign = 0.
         do i=1,4
            i1 = i+1
            i1 = mod1(i1,4)
            call sub3(ee,xx(1,i1),xx(1,i),3)
            call sub3(vv,xj(1)   ,xx(1,i),3)
c
            vnorm = dotprod(vv,vv)
            if (vnorm.gt.0) vnorm = sqrt(vnorm)
c
            if (vnorm.gt.epsn) then
               call norm3d(vv)
               call norm3d(ee)
               scal = 1.0 - abs(dotprod(vv,ee))  ! check for co-linearity
               if (scal.gt.epsr) then
                  call cross(vnt,ee,vv)
                  dotp = dotprod(vnt,nhat)
                  if (sign.eq.0) then
                     sign = dotp
                  elseif (sign*dotp.lt.0) then
                     icond=3+i
                     return
                  endif
               endif
            else
               goto 11
            endif
         enddo
   11    continue
      enddo
      ifok  = .true.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_right_hand(nl)
      include 'basics.inc'
c
      real xyz(2,4)
c
      do iel=1,nl
         DO 210 I=1,4
            XYZ(1,I)= x(iel,i)
            XYZ(2,I)= y(iel,i)
  210    CONTINUE
C
C        CRSS2D(A,B,O) = (A-O) X (B-O)
C                      (note that the notation for corner number here
C                       differs slightly from that used in NEKTON the code.)
         C1=CRSS2D(XYZ(1,2),XYZ(1,4),XYZ(1,1))
         C2=CRSS2D(XYZ(1,3),XYZ(1,1),XYZ(1,2))
         C3=CRSS2D(XYZ(1,4),XYZ(1,2),XYZ(1,3))
         C4=CRSS2D(XYZ(1,1),XYZ(1,3),XYZ(1,4))
C
         IERR=0
         IF (C1.LE.0.0.AND.C2.LE.0.0.AND.
     $       C3.LE.0.0.AND.C4.LE.0.0 ) THEN
C            cyclic permutation (counter clock-wise):  reverse
             call prsi('Reversing Element: $',iel)
             x(iel,2) = xyz(1,4)
             y(iel,2) = xyz(2,4)
             x(iel,4) = xyz(1,2)
             y(iel,4) = xyz(2,2)
             IF (IF3D) THEN
                DO 400 I=1,4
                   x(iel,i+4)=x(iel,i)
                   y(iel,i+4)=y(iel,i)
  400           CONTINUE
             ENDIF
         ELSEIF (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $           C3.LE.0.0.OR.C4.LE.0.0 ) THEN
             CALL PRSi('ERROR in entering element.  Re-enter.$',iel)
             IERR=1
         ENDIF
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine special_delete
      include 'basics.inc'
c
      integer dflag(nelm),e,emin,ecount
c
      open(unit=49,file='cyls.dat',status='old',err=999)
c
      read(49,*) ncyl
c
c     Delete all elements within radius
c
c
      call izero(dflag,nel)
c
      do k=1,ncyl
         read(49,*) rad,x0,y0
         r1 = .99*rad
         r1 = r1*r1
         do e=1,nel
         do i=1,4
            dx = x0-x(e,i)
            dy = y0-y(e,i)
            r2 = dx*dx + dy*dy
            if (r2.lt.r1) dflag(e) = 1
         enddo
         enddo
      enddo
      close(49)
c
      emin = nel+1
      do e=1,nel
         if (dflag(e).ne.0) then
            emin=e
            goto 10
         endif
      enddo
   10 continue
c
      ecount = emin-1
      do e=emin+1,nel
         if (dflag(e).eq.0) then
            ecount = ecount + 1
            call copyel(e,ecount)
         endif
      enddo
      ndel = nel - ecount
      call prsis('Special delete of$',ndel,' elements.$')
      nel  = ecount
c
  999 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_trans(iname,xtr,ytr,ztr)
c
      include 'basics.inc'
c
      character*70 iname
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)
      integer e,f
c
      neli = nel
c
      call chcopy(fname,iname,70)
      if (indx1(fname,'.rea',4).eq.0) then  ! append .rea, if not present
         len = ltrunc(fname,70)
         call chcopy (fnam1(len+1),'.rea',4)
      endif
      open(unit=47,file=fname,status='old',err=1000)
c
c
c     Scan until 'MESH DATA' is found:
c
      call readscan('MESH DATA',9,47)
      read (47,*) nelin,ndimn
      neln = nelin + nel
c
      if (ndimn.ne.ndim) then
         call prs('Dimension must match dimension of current session.$')
         call prs('Returning.$')
         close(47)
         return
      endif
c
      if (neln.ge.nelm) then
         call prs('Sorry, number of elements would exceed nelm.$')
         call prs('Returning.$')
         call prii(neln,nelm)
         close(47)
         return
      endif
c
c     Read geometry
c
      nelo = nel
      nels = nel+1
c
      ierr=imp_geom(x(nels,1),y(nels,1),z(nels,1),nelm
     $             ,numapt(nels),letapt(nels),igroup(nels)
     $             ,ndim,nelin,47)
      if (ierr.ne.0) then
         call prs('Error reading geometry... returning.$')
         close(47)
         return
      endif
c
      ierr=imp_curv(nc,ccurve,curve,ndim,nelin,nel,47)
      if (ierr.ne.0) then
         call prs('Error reading curve side info... returning.$')
         close(47)
         return
      endif
c
c     Translate mesh
c
      write(6,*) neli,neln,xtr,ytr,ztr,'  XTR'
      do e=neli+1,neln
         do i=1,8
            x(e,i) = x(e,i) + xtr
            y(e,i) = y(e,i) + ytr
            z(e,i) = z(e,i) + ztr
         enddo
         do f=1,6
            if (ccurve(f,e).eq.'s') then
               curve(1,f,e) = curve(1,f,e) + xtr
               curve(2,f,e) = curve(2,f,e) + ytr
               curve(3,f,e) = curve(3,f,e) + ztr
            endif
         enddo
      enddo
c
c
c     Read BC info
c
c
      read(47,3) d
    3 format(a3)
c
      ifld0 = 1
      if (.not.ifflow) ifld0=2
      write(6,*) 'IFLD:',ifld0,nflds,ifflow,ifheat
c
      do ifld=ifld0,nflds
         if (.not.ifflow) read(47,*) ans  ! dummy read
         ierr
     $   =imp_bc(cbc(1,nels,ifld),bc(1,1,nels,ifld),ndim,nelin,nel,47)
         if (ierr.ne.0) then
            call prsii('nelin,ifld:$',nelin,ifld)
            call prs('Error reading boundary conditions. Returning.$')
            close(47)
            return
         endif
      enddo
c
c     ALL DONE
c
      close(47)
c
      nel = neln
      ncurve = ncurve + nc
c
      call gencen
c     call redraw_mesh
c
      return
c
 1000 continue
      write(6,*) iname
      write(6,*) fname
      call prs('Unable to open file.  Returning.$')
      close(47)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_special
c
      include 'basics.inc'
      character*70 iname
c
      open(unit=84,file='t.dat',status='old',err=999)
      read(84,*) nsph 
c
      call blank(iname,70)
      write(iname,10) 'base'
   10 format(a4)
c
      ztr = 0.
      do i=1,nsph
         read(84,*) xtr,ytr
         write(6,*) xtr,ytr,' xtr,ytr'
         call imp_mesh_trans(iname,xtr,ytr,ztr)
      enddo
c
  999 return
      end
c-----------------------------------------------------------------------
      subroutine set_std_bcs(nelo)
      include 'basics.inc'
      integer e

      do e=nelo+1,nel
         cbc(1,e,1) = 'W  '
         cbc(2,e,1) = 'W  '
         cbc(3,e,1) = 'W  '
         cbc(4,e,1) = 'W  '
         cbc(5,e,1) = 'v  '
         cbc(6,e,1) = 'O  '
         numapt(e) = 1
         letapt(e) = 'A'
      enddo

      return
      end
c-----------------------------------------------------------------------
