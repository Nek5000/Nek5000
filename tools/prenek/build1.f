      SUBROUTINE NEWEL(IEL,XMOUSE,YMOUSE,BUTTON,IERR)
C     Adds new element to mesh
      INCLUDE 'basics.inc'
      DIMENSION XYZ(2,4)
      LOGICAL IFCURV
      CHARACTER KEY,STRING*6
      DIMENSION IOBJS(8)
C
C
      DO 500 ICORN=1,4
          GO TO 100
   90     CALL PRS('Error reading input.  To input corner, use $')
          CALL PRS('mouse$')
  100     IF(ICORN.EQ.1)CALL PRS(' Corner I  >$')
          IF(ICORN.EQ.2)CALL PRS(' Corner II >$')
          IF(ICORN.EQ.3)CALL PRS(' Corner III>$')
          IF(ICORN.EQ.4)CALL PRS(' Corner IV >$')
          IF(XMOUSE.LE.XPHY(1.0) .AND. ICORN.EQ.1) THEN
C            Special Kludge for "ADD ELEMENT"
C            First element corner already input
          ELSE
             XOLD=XMOUSE
             YOLD=YMOUSE
             CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
             IF (XMOUSE.EQ.XOLD.AND.YMOUSE.EQ.YOLD) GOTO 90
          ENDIF
C !!?? Do we want a cover on keypad now?
          IF (XMOUSE.GT.XPHY(1.0)) THEN
C               look for a keypad input
C               BUG IN KEYOLD IF OLD KEY WAS OFF KEYPAD: IT WONT LET YOU REPAIR!
                CALL PRS('SWITCHING TO KEYPAD INPUTS$')
                CALL SGVIS(11,0)
                CALL SGVIS(12,1)
                CALL PRS('ENTER X:$')
                CALL KEYPAD(XXX)
                CALL PRS('NOW ENTER Y:$')
                CALL KEYPAD(YYY)
                XMOUSE=XXX
                YMOUSE=YYY
                BUTTON='LEFT'
C               TURN OFF KEYPAD HILITE
                CALL SGVIS(12,0)
                CALL SGVIS(11,1)
          ENDIF
          IF(BUTTON.EQ.'LEFT') THEN
c           new unique corner
C
C           Round New Unique Corners - 7-25-90 pff
C
            XMOUSE=ROUND(XMOUSE) 
            YMOUSE=ROUND(YMOUSE)
            X(IEL,ICORN)=XMOUSE
            Y(IEL,ICORN)=YMOUSE
            IOBJS(ICORN)=IOBJCT
          ELSE IF(BUTTON.EQ.'RIGHT')THEN
C           Latch to closest point
            XC=XMOUSE
            YC=YMOUSE
C           Find Closest Corner (NOT IN SAME ELEMENT)
            RMIN = 1.0E8
            DO 150 IIEL=1,IEL-1
               IF(IF3D .AND. NUMAPT(IIEL).NE.NUMAPT(IEL)) GO TO 140
               DO 130 IICORN=1,4
                  XT=X(IIEL,IICORN)
                  YT=Y(IIEL,IICORN)
                  R=(XT-XC)**2+(YT-YC)**2
                  IF(R.LT.RMIN) THEN
                     RMIN  = R
                     IELMIN= IIEL
                     ICMIN = IICORN
                  ENDIF
  130           CONTINUE
  140         CONTINUE
  150       CONTINUE
            X(IEL,ICORN)=X(IELMIN,ICMIN)
            Y(IEL,ICORN)=Y(IELMIN,ICMIN)
C...if its attached to an object, define it as so.
            IEDG1 = ICMIN-1
            IEDG2 = ICMIN
            IF (IEDG1.EQ.0) IEDG1=4
C  NOTE: a bug in the logic here, one pt. CAN be a member of two objects..!
C        a SIDE can be a member of only ONE object, must make distinction.
            IOBJS(ICORN)=0
            IF (CCURVE(IEDG1,IELMIN).EQ.'O')
     $         IOBJS(ICORN)=INT(CURVE(1,IEDG1,IELMIN))
            IF (CCURVE(IEDG2,IELMIN).EQ.'O')
     $         IOBJS(ICORN)=INT(CURVE(1,IEDG2,IELMIN))
          ELSE
            CALL PRS('Error reading mouse input$')
            CALL BEEP
            GOTO 90
          ENDIF
          IF(IF3D)THEN
             X(IEL,ICORN+4)=X(IEL,ICORN)
             Y(IEL,ICORN+4)=Y(IEL,ICORN)
          ENDIF
          CALL color(3)
          CALL fillp(-15)
          XS=X(IEL,ICORN)
          YS=Y(IEL,ICORN)
C
          IF(ICORN.EQ.1) THEN
            CALL BEGINB(XS,YS)
            XS1=XS
            YS1=YS
            IF(IN.EQ.5)CALL RUBBER
          ELSE
            CALL DRAWC(XS,YS)
          ENDIF
          IF(ICORN.EQ.4) THEN
C           Turn off rubberbanding
            CALL OFFRUB
            CALL DRAWC(XS1,YS1)
            CALL ENDP
C           Check Here for points not entered counter-clockwise
            DO 210 I=1,4
               XYZ(1,I)= x(iel,i)
               XYZ(2,I)= y(iel,i)
  210       CONTINUE
C
C           CRSS2D(A,B,O) = (A-O) X (B-O)
C                         (note that the notation for corner number here
C                          differs slightly from that used in NEKTON the code.)
            C1=CRSS2D(XYZ(1,2),XYZ(1,4),XYZ(1,1))
            C2=CRSS2D(XYZ(1,3),XYZ(1,1),XYZ(1,2))
            C3=CRSS2D(XYZ(1,4),XYZ(1,2),XYZ(1,3))
            C4=CRSS2D(XYZ(1,1),XYZ(1,3),XYZ(1,4))
C
            IERR=0
            IF (C1.LE.0.0.AND.C2.LE.0.0.AND.
     $          C3.LE.0.0.AND.C4.LE.0.0 ) THEN
C               cyclic permutation (counter clock-wise):  reverse
                x(iel,2) = xyz(1,4)
                y(iel,2) = xyz(2,4)
                x(iel,4) = xyz(1,2)
                y(iel,4) = xyz(2,2)
                IF (IF3D) THEN
                   DO 400 I=1,4
                      x(iel,i+4)=x(iel,i)
                      y(iel,i+4)=y(iel,i)
  400              CONTINUE
                ENDIF
            ELSEIF (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $              C3.LE.0.0.OR.C4.LE.0.0 ) THEN
                CALL PRSI('ERROR in entering element.  Re-enter.$',iel)
                IERR=1
                RETURN
            ENDIF
          ENDIF
  500 CONTINUE
c     Label Corners
      XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
      YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
      ILETAP=ILETAP+1
      IF(ILETAP.LE.122)LETAPT(IEL)=CHAR(ILETAP)
      IF(ILETAP.GT.122)LETAPT(IEL)=' '
C     New elements have their Z-coordinate put in here
      IF(NDIM.EQ.3)THEN
         ZCEIL=0.0
         DO 620 I=1,NUMAPT(IEL)
            ZCEIL=ZCEIL+HEIGHT(I)
  620    CONTINUE
         ZFLOOR=ZCEIL-HEIGHT(NUMAPT(IEL))
         DO 630 IC=1,8
           IF(IC.LE.4)Z(IEL,IC)=ZFLOOR
           IF(IC.GT.4)Z(IEL,IC)=ZCEIL
  630    CONTINUE
      ENDIF
C
C     Check for curve sides here if using polar coordinates or objects
C
      IFCURV=.FALSE.
      DO 700 IC=1,4
         IC1 = IC+1
         IC1 = MOD1(IC1,4)
         IF (IOBJS(IC).EQ.IOBJS(IC1).AND.IOBJS(IC).NE.0) IFCURV=.TRUE.
  700 CONTINUE
      IF (IFCURV) THEN
         CALL DRAWEL(-IEL)
         DO 710 IC=1,4
            IC1 = IC+1
            IC1 = MOD1(IC1,4)
            IF (IOBJS(IC).EQ.IOBJS(IC1).AND.IOBJS(IC).NE.0) THEN
               iobjc=iobjs(ic)
               CCURVE(IC,IEL)=ccobjs(iobjc)
               call copy(CURVE(1,IC,IEL),cobjs(1,iobjc),6)
               IF (IF3D) THEN
                  IC4 = IC+4
                  CCURVE(IC4,IEL)=CCURVE(IC,IEL)
                  call copy(CURVE(1,IC4,IEL),CURVE(1,IC,IEL),6)
               ENDIF
            ENDIF
  710    CONTINUE
c     end of object element check
      ENDIF
      IF (IFGRDP) THEN
c     beginning of polar grid check
         RADUSC = 0.0
         DO 800 IC=1,4
            DX = X(IEL,IC)-GRIDXP
            DY = Y(IEL,IC)-GRIDYP
            RADUSC = RADUSC + SQRT(DX**2+DY**2)
  800    CONTINUE
         RADUSC = RADUSC/4.0
C
C        check for any curve sides, if yes, delete straight sided image.
C        (two steps required to delete proper shaped element)
         IFCURV=.FALSE.
         DO 801 IC=1,4
            IC1 = IC+1
            IF (MOD(IC,4).EQ.0) IC1 = IC-3
            DX0=X(IEL,IC )-GRIDXP
            DY0=Y(IEL,IC )-GRIDYP
            DX1=X(IEL,IC1)-GRIDXP
            DY1=Y(IEL,IC1)-GRIDYP
            RADUS0 = SQRT(DX0**2+DY0**2)
            RADUS1 = SQRT(DX1**2+DY1**2)
            DRAD = (RADUS1-RADUS0)/(RADUS1+RADUS0)
            IF (ABS(DRAD).LT.0.0001) IFCURV=.TRUE.
  801    CONTINUE
c
c        Check for element being adjacent to another which already
c        has a curved side
c
         do 819 ic = 1,4
            jc = ic+1
            jc = mod1(jc,4)
            xm = 0.5*(x(iel,ic)+x(iel,jc))
            ym = 0.5*(y(iel,ic)+y(iel,jc))
            delta = (x(iel,ic)-x(iel,jc))**2
     $            + (y(iel,ic)-y(iel,jc))**2
            call getside(jel,jic,xm,ym)
            jjc = jic+1
            jjc = mod1(jjc,4)
            xj = 0.5*(x(jel,jic)+x(jel,jjc))
            yj = 0.5*(y(jel,jic)+y(jel,jjc))
            delta2 = ((xm-xj)**2 + (ym-yj)**2)/delta
            if (delta2.lt.0.002) then
c              matched sides, give curve side attributes to new
c              element.... pff 3/27/94:  only "C" curve sides for now.
c
               if (CCURVE(jic,jel).eq.'C') then
                  CCURVE(ic,iel)='C'
                  CURVE(1,ic,iel) =  -CURVE(1,jjc,jel)
                  ifcurv = .true.
               endif
            endif
  819    continue
C
         IF (IFCURV) THEN
            CALL DRAWEL(-IEL)
C           Now curve side and update element image.
            DO 802 IC=1,4
               IC1 = IC+1
               IF (MOD(IC,4).EQ.0) IC1 = IC-3
               DX0=X(IEL,IC )-GRIDXP
               DY0=Y(IEL,IC )-GRIDYP
               DX1=X(IEL,IC1)-GRIDXP
               DY1=Y(IEL,IC1)-GRIDYP
               RADUS0 = SQRT(DX0**2+DY0**2)
               RADUS1 = SQRT(DX1**2+DY1**2)
               DRAD = (RADUS1-RADUS0)/(RADUS1+RADUS0)
               IF (ABS(DRAD).LT.0.0001) THEN
C                 curve side
                  IF (RADUS0.LT.RADUSC) RADUS0 = -RADUS0
                  CURVE(1,IC,IEL)=RADUS0
                  CCURVE(IC,IEL)='C'
                  IF (IF3D) THEN
                     IC4 = IC+4
                     CURVE(1,IC4,IEL)=RADUS0
                     CCURVE(IC4,IEL)='C'
                  ENDIF
               ENDIF
  802       CONTINUE
         ENDIF
C     End of polar grid check
      ENDIF
C
C     Element all set , draw it.
C
 1000 CONTINUE
      CALL DRAWEL(IEL)
      CALL SORTEL
C     Figure out which to draw
      DO 1040 I=1,NEL
         IF(ISRT(I).EQ.NEL) IBEGIN=I
 1040 CONTINUE
      DO 1050 I=IBEGIN,NEL
         CALL DRAWIS(ISRT(I))
 1050 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine delelq(idel)
C     Deletes element from mesh.  But really, it copies a null element to
C     the one to be deleted.
      INCLUDE 'basics.inc'
C
c     CALL BLANK(S,80)
c     WRITE(S,10,ERR=20) IDEL,NEL
c  10 FORMAT(' Deleting element',I5,' of',I5,'.$')
c     CALL PRS(S)
   20 CONTINUE
C
      if (idel.ne.nel) CALL COPYEL(NEL,IDEL)
      CALL COPYEL(NELM,NEL)
C     If the last element was a conduction element, reduce NCOND
      IF(MASKEL(NEL,1).EQ.0) NCOND=NCOND-1
      NEL=NEL-1
C     ILETAP=ILETAP-1
C
C     Recount the number of curved sides
C
      NCURVE=0
      DO 100 IE=1,NEL
      DO 100 IEDGE=1,8
         IF(CCURVE(IEDGE,IE).NE.' ') NCURVE=NCURVE+1
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine delel(idel)
C     Deletes element from mesh.  But really, it copies a null element to
C     the one to be deleted.
      INCLUDE 'basics.inc'
C
      CALL BLANK(S,80)
      WRITE(S,10,ERR=20) IDEL,NEL
   10 FORMAT(' Deleting element',I5,' of',I5,'.$')
      CALL PRS(S)
   20 CONTINUE
C
      if (idel.ne.nel) CALL COPYEL(NEL,IDEL)
      CALL COPYEL(NELM,NEL)
C     If the last element was a conduction element, reduce NCOND
      IF(MASKEL(NEL,1).EQ.0) NCOND=NCOND-1
      NEL=NEL-1
C     ILETAP=ILETAP-1
C
C     Recount the number of curved sides
C
      NCURVE=0
      DO 100 IE=1,NEL
      DO 100 IEDGE=1,8
         IF(CCURVE(IEDGE,IE).NE.' ') NCURVE=NCURVE+1
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine copyel(isrc,ides)
      INCLUDE 'basics.inc'
      CHARACTER KEY,STRING*6
C
C     Copies everything vaguely related to element
      IF(IDES.EQ.NELM-2)CALL PRS
     $('***WARNING*** NO ADDITIONAL ELEMENTS ALLOWED: ARRAYS FULL$')
      IF(IDES.EQ.NELM  )THEN
        CALL PRS
     $('***ERROR*** TOO MANY ELEMENTS!  SEE YOUR NEKTONICS REP.$')
        write(6,*) 'ides,nelm:',IDES,NELM
        RETURN
      ENDIF
C
      DO 100 IFLD=0,MAXFLD
      DO 100 I=1,6
         CBC(I,IDES,IFLD) = CBC(I,ISRC,IFLD)
         DO 100 J=1,5
            BC(J,I,IDES,IFLD) = BC(J,I,ISRC,IFLD)
  100 CONTINUE
C
      DO 200 I =1,6
      DO 200 II=1,4
            SIDES(IDES,I,II)=SIDES(ISRC,I,II)
  200 CONTINUE
C
      DO 300 IC=1,8
         X(IDES,IC)=X(ISRC,IC)
         Y(IDES,IC)=Y(ISRC,IC)
         Z(IDES,IC)=Z(ISRC,IC)
  300 CONTINUE
      IGROUP(IDES)=IGROUP(ISRC)
C
C     Correct NCOND
C
      IF(MASKEL(IDES,1).EQ.0 .AND. MASKEL(ISRC,1).EQ.1)NCOND=NCOND-1
      IF(MASKEL(IDES,1).EQ.1 .AND. MASKEL(ISRC,1).EQ.0)NCOND=NCOND+1
      DO 500 IF=1,MAXFLD
         MASKEL(IDES,IF)=MASKEL(ISRC,IF)
  500 CONTINUE
      ISRT  (IDES) = ISRT (ISRC)
      ICRV  (IDES) = ICRV (ISRC)
      XCEN  (IDES) = XCEN (ISRC)
      YCEN  (IDES) = YCEN (ISRC)
      ZCEN  (IDES) = ZCEN (ISRC)
      RCEN  (IDES) = RCEN (ISRC)
      NUMAPT(IDES) =NUMAPT(ISRC)
      LETAPT(IDES) =LETAPT(ISRC)
      DO 600 IEDGE=1,8
         CCURVE(IEDGE,IDES) = CCURVE(IEDGE,ISRC)
         DO 550 IDIR=1,3
            EDGES(IEDGE,IDIR,IDES)=EDGES(IEDGE,IDIR,ISRC)
  550    CONTINUE
  600 CONTINUE
      DO 610 IEDGE=1,8
      DO 610 I=1,6
         CURVE(I,IEDGE,IDES) = CURVE(I,IEDGE,ISRC)
  610 CONTINUE
  611 FORMAT(' curve',4I5,2E14.5)
      IF(NHIS.GT.0)THEN
C        Delete any integral sides
         IDEL=0
         DO 700 I=1,NHIS
            IF(HCODE(10,I).EQ.'I')IDEL=1
            IF(HCODE(10,I).EQ.'I')HCODE(10,I)=' '
  700    CONTINUE
         IF(IDEL.NE.0)THEN
            CALL PRS('***RESETTING INTEGRAL QUANTITY***$')
            CALL PRS(
     $   ' YOU MUST REDEFINE SIDES TO GET INTEGRAL FLUX, LIFT,OR DRAG$')
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE MODEL(IEL)
C     Modifies element by moving point in mesh
C     Modified neighbor points iff they had been latched to point in questio
      INCLUDE 'basics.inc'
      DIMENSION IELMOV(NELM)
      REAL XCHECK(8),YCHECK(8)
      CHARACTER KEY,STRING*6
      LOGICAL IFTMP
C
C
C     FIND CLOSEST ELEMENT AND CORNER
C     Erase whole isometric surface
c
      CALL PRS
     $('Push left button close to the corner you want changed.$')
      IFTMP =IFGRID
      IFGRID=.FALSE.
1     CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0) THEN
C        He apparently is trying to use the keypad
         CALL PRS('Enter X-coordinate with keypad:$')
         CALL KEYPAD(XMOUSE)
         CALL PRS('Now enter Y-coordinate with keypad:$')
         CALL KEYPAD(YMOUSE)
      ELSE IF(BUTTON.EQ.'RIGHT')THEN
C        Latch to closest element vertex
         CALL BLATCH(XMOUSE,YMOUSE)
      ENDIF
      XM=XMOUSE
      YM=YMOUSE
      IFGRID=IFTMP
C     If the areas from the triangles formed by each side and the mouse point
C     are all positive (ccw) we must be inside the element!
      IELMV=0
      DO 2 I=1,IEL
C        Skip apartments on other floors
         IF(.NOT. IF3D .OR. NUMAPT(I).EQ.ILEVEL)THEN
            IF(IFCEIL)THEN
               A1 = (x(i,5)-XM)*(y(i,6)-YM)-(x(i,6)-XM)*(y(i,5)-YM)
               A2 = (x(i,6)-XM)*(y(i,7)-YM)-(x(i,7)-XM)*(y(i,6)-YM)
               A3 = (x(i,7)-XM)*(y(i,8)-YM)-(x(i,8)-XM)*(y(i,7)-YM)
               A4 = (x(i,8)-XM)*(y(i,5)-YM)-(x(i,5)-XM)*(y(i,8)-YM)
            ELSE
               A1 = (x(i,1)-XM)*(y(i,2)-YM)-(x(i,2)-XM)*(y(i,1)-YM)
               A2 = (x(i,2)-XM)*(y(i,3)-YM)-(x(i,3)-XM)*(y(i,2)-YM)
               A3 = (x(i,3)-XM)*(y(i,4)-YM)-(x(i,4)-XM)*(y(i,3)-YM)
               A4 = (x(i,4)-XM)*(y(i,1)-YM)-(x(i,1)-XM)*(y(i,4)-YM)
            ENDIF
4           CONTINUE
C           We're inside the element
            IF(A1.GE.0.0.AND.A2.GE.0.0.AND.A3.GE.0.0.AND.A4.GE.0.0)
     $      IELMV=I
         ENDIF
2     CONTINUE
      IF(IELMV.EQ.0)THEN
C        If it still = 0 then the point entered must have been outside all
C        elements.
         CALL PRS
     $   ('Enter a point INSIDE one of the elements in order to$')
         CALL PRS('move one of the element''s corners.  Try again.$')
         GO TO 1
      ENDIF
C
      IF(IFCEIL)THEN
         IC1=5
         IC2=8
      ELSE
         IC1=1
         IC2=4
      ENDIF
C     FIND CLOSEST CORNER
      RMIN=1.0E8
      DO 3 I=IC1,IC2
         R=SQRT((X(IELMV,I)-XMOUSE)**2+(Y(IELMV,I)-YMOUSE)**2)
         IF(R.LT.RMIN) THEN
            RMIN=R
            ICOMV=I
         ENDIF
3     CONTINUE
      XPICKED=X(IELMV,ICOMV)
      YPICKED=Y(IELMV,ICOMV)
c
      CALL PRS
     $('Enter new point to which element corner is to be moved:$')
      CALL MOUSE(XMOVED,YMOVED,BUTTON)
      IF(XSCR(XMoved).GT.1.0) THEN
C        He apparently is trying to use the keypad
         CALL PRS('Enter X-coordinate with keypad:$')
         CALL KEYPAD(XMoved)
         CALL PRS('Now enter Y-coordinate with keypad:$')
         CALL KEYPAD(YMoved)
      ELSEIF (BUTTON.EQ.'RIGHT') THEN
C        Latch to closest element vertex
         CALL BLATCH(XMOVED,YMOVED)
      ENDIF
C
C     HAVE IELMV,ICOMV  , NOW MOVE THE APPROPRIATE POINTS AND REDRAW ELEMENTS
C
C     First loop checks legality, second moves and erases, third redraws
      CALL DRAWIS(-IELMV)
      DO 10 ILOOP=1,3
       DO 10 IIEL=1,IEL
         IF(ILOOP.EQ.1) IELMOV(IIEL)=0
         DO 10 IICORN=IC1,IC2
C           Only move corners on same floor
            IF(NUMAPT(IIEL).EQ.NUMAPT(IELMV)) THEN
C              Check floor for same global #'s  NO: Floors automatically
C              Modify overlapping corners; ceilings move 1 at a time.
C              Skip this check on 3rd loop; the redraw will know which to do
               IF(ILOOP.NE.3)THEN
C                    Only move corners close to the one picked
                     IF (ABS(X(IIEL,IICORN)-XPICKED).GT.XFAC/100.
     $               .OR.ABS(Y(IIEL,IICORN)-YPICKED).GT.YFAC/100.)
     $               GO TO 9
c                     CALL PRSII('MOVING$',IIEL,IICORN)
               ENDIF
               IF(ILOOP.EQ.1)THEN
C
C                 First check if move is legal (has 4 acute <'s for corners
                  XTEST=X(IIEL,IICORN)
                  YTEST=Y(IIEL,IICORN)
                  DO 8 I=IC1,IC2
                     XCHECK(I)=X(IIEL,I)
                     YCHECK(I)=Y(IIEL,I)
8                 CONTINUE
                  XCHECK(IICORN)=XMOVED
                  YCHECK(IICORN)=YMOVED
C                 Check that the 4 cross products (areas) of each adjacent
C                 side pair is positive (angles between 0 and 180 degrees)
                  DO 109 ICC=1,4
                     I1=MOD(ICC-1,4)+1
                     I2=MOD(ICC  ,4)+1
                     I3=MOD(ICC+1,4)+1
                     IF(IFCEIL)THEN
                        I1=I1+4
                        I2=I2+4
                        I3=I3+4
                     ENDIF
                     X1=XCHECK(I2)-XCHECK(I1)
                     X2=XCHECK(I3)-XCHECK(I2)
                     Y1=YCHECK(I2)-YCHECK(I1)
                     Y2=YCHECK(I3)-YCHECK(I2)
                     AREA=X1*Y2-X2*Y1
                     IF(AREA.LE.0.0) THEN
                        CALL PRSIS(
     $                  '**ERROR** Element$',IIEL,'Cannot be modified$')
                        CALL PRSIS(
     $                  'Angle at corner $',I2,' Would have been $')
                        CALL PRS('illegal (>180 degrees)$')
                        CALL PRS('No Elements were modified.$' )
                        DO 52 I=1,NEL
c                           CALL DRAWIS(ISRT (I))
52                      CONTINUE
                        RETURN
                     ENDIF
109                 CONTINUE
C
               ELSE IF(ILOOP.EQ.2)THEN
C                 Now, erase old element lines
                  CALL DRAWEL(-IIEL)
                  CALL DRAWIS(-IIEL)
C                 Move corners
                  IELMOV(IIEL)=1
                  IF(IFCEIL)THEN
C                  First move points on floor of upper level
                    DO 18 I=1,IEL
                     MOVEDE=0
                     DO 16 IC=1,4
                       IF(NUMAPT(I).EQ.ILEVEL+1)THEN
                         IF  (ABS(X(I,IC)-X(IIEL,IICORN)).LT.XFAC/100.
     $                   .AND.ABS(Y(I,IC)-Y(IIEL,IICORN)).LT.YFAC/100.)
     $                   THEN
c                            CALL PRS('moving upper$')
                            MOVEDE=1
                            X(I,IC)=XMOVED
                            Y(I,IC)=YMOVED
                         ENDIF
                       ENDIF
16                    CONTINUE
                      IF(MOVEDE.EQ.1)CALL DRAWIS(-IEL)
18                    CONTINUE
                  ELSE IF(.NOT.IFCEIL)THEN
C                    First move points on ceiling of lower level
                     IF(ILEVEL.GT.1)THEN
                      DO 14 I=1,IEL
                       MOVEDE=0
                       DO 12 IC=5,8
                         IF(NUMAPT(I).EQ.ILEVEL-1)THEN
                           IF(ABS(X(I,IC)-X(IIEL,IICORN)).LT.XFAC/100.
     $                   .AND.ABS(Y(I,IC)-Y(IIEL,IICORN)).LT.YFAC/100.)
     $                     THEN
                              X(I,IC)=XMOVED
                              Y(I,IC)=YMOVED
                              MOVEDE=1
                           ENDIF
                         ENDIF
12                      CONTINUE
                        IF(MOVEDE.EQ.1)CALL DRAWIS(-IEL)
14                     CONTINUE
                     ENDIF
C                    Here we also move the points on the ceiling of the element,
C                    But only if there are no elements on any levels above
C                    the current one. (to avoid complications)
                     DO 11 I=1,IEL
                        IF(NUMAPT(I).GT.ILEVEL) THEN
                           CALL PRS('Presence of elements on higher '//
     $                     'floors inhibits modifying ceiling also$')
                           GO TO 15
                        ENDIF
11                   CONTINUE
C                    Since we got here, there must be no elements on upper floo
                     X(IIEL,IICORN+4)=XMOVED
                     Y(IIEL,IICORN+4)=YMOVED
                  ENDIF
                  IF(IFCEIL.AND.IICORN.EQ.5)THEN
C                    Manually erase old ceiling
                     CALL color(0)
                     CALL MOVEC(X(IIEL,8),Y(IIEL,8))
                     DO 114 IC=5,8
                        CALL DRAWC(X(IIEL,IC),Y(IIEL,IC))
114                  CONTINUE
                     CALL color(10)
                  ENDIF
C                 Now move points on this level
15                X(IIEL,IICORN)=XMOVED
                  Y(IIEL,IICORN)=YMOVED
C                  CALL PRS('IIEL,IICORN,XMOVED$')
C                  CALL PRS(IIEL,IICORN,XMOVED
C                 Move Center
                  IF(.NOT.IFCEIL)THEN
                  XCEN(IIEL)=(X(IIEL,1)+X(IIEL,2)+X(IIEL,3)+X(IIEL,4))/4
                  YCEN(IIEL)=(Y(IIEL,1)+Y(IIEL,2)+Y(IIEL,3)+Y(IIEL,4))/4
                  ENDIF
               ELSE IF(ILOOP.EQ.3)THEN
C                 Now, draw new elements for elements that were modified
                  IF(IICORN.NE.IC2)GO TO 9
                  IF(IELMOV(IIEL).NE.0)CALL DRAWEL(IIEL)
               ENDIF
            endif
9           CONTINUE
10    CONTINUE
C     Draw ISOMETRICALLY ONLY ELEMENT MOVED
      CALL DRAWIS(IELMV)
      RETURN
      END
      SUBROUTINE DRISED(IEL,IEDGE,IFLIP)
C     DRaws ISometric EDge.  IFLIP CAuses to draw from end to beginning
C     POSITIVE MEANS CCW; ON VERTICAL STRUTS 9-12, MEANS UPWARD
C     Draws only.  No move or fill.
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      INCLUDE 'basics.inc'
C

      IF(IEDGE.GT.8)THEN
C        Vertical strut
         IC=IEDGE-4
         XISOM(IC)=XPHY(XSCR(X(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         YISOM(IC)=YPHY(YSCR(Y(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         IC=IEDGE-8
         XISOM(IC)=XPHY(XSCR(X(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         YISOM(IC)=YPHY(YSCR(Y(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         IF(IFLIP.EQ.1)THEN
            CALL MOVEC(XISOM(IEDGE-8),YISOM(IEDGE-8))
            CALL DRAWC(XISOM(IEDGE-4),YISOM(IEDGE-4))
         ELSE
            CALL MOVEC(XISOM(IEDGE-4),YISOM(IEDGE-4))
            CALL DRAWC(XISOM(IEDGE-8),YISOM(IEDGE-8))
         ENDIF
      ELSE IF(CCURVE(IEDGE,IEL).EQ.' ')THEN
C        Draw straight side
         IF(IFLIP.EQ.1)THEN
            IC=IEDGE+1
            IF(IC.EQ.5)IC=1
            IF(IC.EQ.9)IC=5
         ELSE
            IC=IEDGE
         ENDIF
         XISOM(IC)=XPHY(XSCR(X(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         YISOM(IC)=YPHY(YSCR(Y(IEL,IC))/5.0 + 0.8)-Z(IEL,IC)/20.
         CALL DRAWC(XISOM(IC),YISOM(IC))
      ELSE
C        Draw curved side
         NPOINT=10
         DO 18 I=1,NPOINT
            CSPACE(I)=(I-1.0)/(NPOINT-1.0)
18       CONTINUE
         CALL GETPTS(NPOINT,CSPACE,IEL,IEDGE,XCRVED,YCRVED)
         IF(IFLIP.EQ.1)THEN
            IBEGIN=1
            IEND  =NPOINT
         ELSE
            IBEGIN=NPOINT
            IEND  =1
         ENDIF
         DO 118 I=IBEGIN,IEND,IFLIP
            XI=XPHY(XSCR(XCRVED(I))/5.0 + 0.8)-Z(IEL,IEDGE)/20.
            YI=YPHY(YSCR(YCRVED(I))/5.0 + 0.8)-Z(IEL,IEDGE)/20.
            CALL DRAWC(XI,YI)
118      CONTINUE
      ENDIF
      RETURN
      END
      SUBROUTINE DRAWIS(IEL)
      INCLUDE 'basics.inc'
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      IIEL=IABS(IEL)
      IF(.NOT.IF3D)RETURN
      IF(nel.gt.1000)RETURN
C        Now draw isometric view  (??! RESCALE??)
         IF(IEL.GT.0)call color(10)
         IF(IEL.le.0)call color(0)
         DO 7 IC=1,8
         XISOM(IC)=XPHY(XSCR(X(IIEL,IC))/5.0 + 0.8)-Z(IIEL,IC)/20.
         YISOM(IC)=YPHY(YSCR(Y(IIEL,IC))/5.0 + 0.8)-Z(IIEL,IC)/20.
7        CONTINUE
         CALL BEGINB(XISOM(1),YISOM(1))
         DO 8 IEDGE=1,4
            CALL DRISED(IIEL,IEDGE,1)
8        CONTINUE
         CALL ENDP
C        Now draw side panels
         IF(IEL.GT.0)CALL fillp(-14)
         IF(IEL.LE.0)CALL fillp(0)
         CALL BEGINB(XISOM(2),YISOM(2))
         CALL DRISED(IIEL, 2, 1)
         CALL DRISED(IIEL,11, 1)
         CALL DRISED(IIEL, 6,-1)
         CALL DRISED(IIEL,10,-1)
         CALL ENDP
C
         CALL BEGINB(XISOM(4),YISOM(4))
         CALL DRISED(IIEL,12, 1)
         CALL DRISED(IIEL, 7,-1)
         CALL DRISED(IIEL,11,-1)
         CALL DRISED(IIEL, 3, 1)
         CALL ENDP
C
C        Draw Ceiling panel
         IF(IEL.GT.0)CALL fillp(-15)
         IF(IEL.LE.0)CALL fillp(0)
         CALL BEGINB(XISOM(5),YISOM(5))
         DO 9 I=5,8
           CALL DRISED(IIEL,I, 1)
9        CONTINUE
         CALL ENDP
      RETURN
      END
      SUBROUTINE DRAWEL(IEL)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
      INCLUDE 'basics.inc'
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      IIEL=IABS(IEL)
C

C HMT TRACE
C      write (6,*) 'DRAWEL :: iel = ',IEL

      IFSPHR=.FALSE.
      DO 7 IEDGE=1,8
         IF (CCURVE(IEDGE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         XYZCTR(1,1)=X(IE,1)+X(IE,4)+X(IE,5)+X(IE,8)
         XYZCTR(1,2)=X(IE,2)+X(IE,3)+X(IE,6)+X(IE,7)
         XYZCTR(1,3)=X(IE,1)+X(IE,2)+X(IE,5)+X(IE,6)
         XYZCTR(1,4)=X(IE,3)+X(IE,4)+X(IE,7)+X(IE,8)
         XYZCTR(1,5)=X(IE,1)+X(IE,2)+X(IE,3)+X(IE,4)
         XYZCTR(1,6)=X(IE,5)+X(IE,6)+X(IE,7)+X(IE,8)
         XYZCTR(2,1)=Y(IE,1)+Y(IE,4)+Y(IE,5)+Y(IE,8)
         XYZCTR(2,2)=Y(IE,2)+Y(IE,3)+Y(IE,6)+Y(IE,7)
         XYZCTR(2,3)=Y(IE,1)+Y(IE,2)+Y(IE,5)+Y(IE,6)
         XYZCTR(2,4)=Y(IE,3)+Y(IE,4)+Y(IE,7)+Y(IE,8)
         XYZCTR(2,5)=Y(IE,1)+Y(IE,2)+Y(IE,3)+Y(IE,4)
         XYZCTR(2,6)=Y(IE,5)+Y(IE,6)+Y(IE,7)+Y(IE,8)
         XYZCTR(3,1)=Z(IE,1)+Z(IE,4)+Z(IE,5)+Z(IE,8)
         XYZCTR(3,2)=Z(IE,2)+Z(IE,3)+Z(IE,6)+Z(IE,7)
         XYZCTR(3,3)=Z(IE,1)+Z(IE,2)+Z(IE,5)+Z(IE,6)
         XYZCTR(3,4)=Z(IE,3)+Z(IE,4)+Z(IE,7)+Z(IE,8)
         XYZCTR(3,5)=Z(IE,1)+Z(IE,2)+Z(IE,3)+Z(IE,4)
         XYZCTR(3,6)=Z(IE,5)+Z(IE,6)+Z(IE,7)+Z(IE,8)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            IF (J.EQ.1) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.2) THEN
               CALL BEGINB(XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL ENDP
            ENDIF
            IF (J.EQ.3) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL ENDP
            ENDIF
            IF (J.EQ.4) THEN
               CALL BEGINB(XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.5) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL ENDP
            ENDIF
            IF (J.EQ.6) THEN
               CALL BEGINB(XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
   72    CONTINUE
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(-14)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
               call fillp(-15)
c              i15 = mod(iel,15)+1
c              call fillp(i15)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(-2)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp( -(IGROUP(IEL)+2) )
            ELSE
               CALL fillp( -(IGROUP(IEL)+3) )
            ENDIF
            IF(IFBWGKS)call fillp(0)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(0)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8
      IF(IEL.GT.10000)THEN
         IIEL=IIEL-10000
         CALL MOVEC(x(Iiel,IC),y(Iiel,IC))
      ELSE
         CALL BEGINB(x(Iiel,IC),y(Iiel,IC))
      ENDIF
      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          YCENTER=YCENTER+Y(IIEL,IC)/4.
          XCENTER=XCENTER+X(IIEL,IC)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          IF(CCURVE(IEDGE,IIEL).EQ.' ')THEN
             CALL DRAWC(X(IIEL,IC),Y(IIEL,IC))
          ELSE
C            Draw curved side
             CALL GETPTS(NPOINT,CSPACE,IIEL,IEDGE,XCRVED,YCRVED)
             DO 118 I=1,NPOINT
                CALL DRAWC(XCRVED(I),YCRVED(I))
118          CONTINUE
          ENDIF
6     CONTINUE
      IF(IEL.LT.10000)CALL ENDP
C
      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
         WRITE(STRING,'(I3)')IEL
         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
C HMT color TRACE - seems to be the outline
      call color(1)
      RETURN
      END
      SUBROUTINE MKSIDE
      INCLUDE 'basics.inc'
C
C     Find Sides' Midpoints
C
      DO 25 IEL=1,NEL
      DO 25 ISIDE=1,NSIDES
         IC1=ISIDE
         IC2=ISIDE+1
         IF(ISIDE.EQ.4)IC2=1
C        This stuff only relevant for 3d
         IC3=IC1+4
         IC4=IC2+4
         IF (ISIDE.EQ.5) THEN
            IC1=1
            IC2=2
            IC3=3
            IC4=4
         ELSEIF (ISIDE.EQ.6) THEN
            IC1=1+4
            IC2=2+4
            IC3=3+4
            IC4=4+4
         ENDIF
         IF (IF3D) THEN
            XS =( X(IEL,IC1)+X(IEL,IC2)
     $         +  X(IEL,IC3)+X(IEL,IC4) )/4.
            YS =( Y(IEL,IC1)+Y(IEL,IC2)
     $         +  Y(IEL,IC3)+Y(IEL,IC4) )/4.
            ZS =( Z(IEL,IC1)+Z(IEL,IC2)
     $         +  Z(IEL,IC3)+Z(IEL,IC4) )/4.
         ELSE
            XS =( X(IEL,IC1)+X(IEL,IC2) )/2.
            YS =( Y(IEL,IC1)+Y(IEL,IC2) )/2.
            ZS = 0.0
         ENDIF
         SIDES (IEL,ISIDE,1)=XS
         SIDES (IEL,ISIDE,2)=YS
         SIDES (IEL,ISIDE,3)=ZS
25    CONTINUE
      RETURN
      END
      FUNCTION CRSS2D(XY1,XY2,XY0)
      REAL XY1(2),XY2(2),XY0(2)
C
         V1X=XY1(1)-XY0(1)
         V2X=XY2(1)-XY0(1)
         V1Y=XY1(2)-XY0(2)
         V2Y=XY2(2)-XY0(2)
         CRSS2D = V1X*V2Y - V1Y*V2X
C
      RETURN
      END
      FUNCTION ROUND(X)
C
C     Try to Round X to fractional integer - eg .05 .1 .15 - if it's within 10-6
C
      EPS=1.0E-4
      XTMP=20.0*X
      XTMP2=XTMP+0.5
      ITMP=INT(XTMP2)
      XTMP2=FLOAT(ITMP)
      IF (ABS(XTMP-XTMP2).LT.EPS) THEN
         ROUND=XTMP2/20.0
      ELSE
         ROUND=X
      ENDIF
      RETURN
      END
      SUBROUTINE BLATCH(XMOUSE,YMOUSE)
C     Latch to closest element vertex
      INCLUDE 'basics.inc'
      XC=XMOUSE
      YC=YMOUSE
C     Find Closest Corner (NOT IN SAME ELEMENT)
      RMIN = 1.0E8
      DO 150 IEL=1,NEL
         DO 130 IICORN=1,4
            ICORN=IICORN
            IF (IFCEIL) ICORN=ICORN+4
            XT=X(IEL,ICORN)
            YT=Y(IEL,ICORN)
            R=(XT-XC)**2+(YT-YC)**2
            IF (R.LT.RMIN) THEN
               RMIN  = R
               IELMIN= IEL
               ICMIN = ICORN
            ENDIF
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
      XMOUSE=X(IELMIN,ICMIN)
      YMOUSE=Y(IELMIN,ICMIN)
      RETURN
      END
C     End of BUILD subroutines, preprocessor *****
      SUBROUTINE SORTEL
      INCLUDE 'basics.inc'
C     Sorts elements according to their visibility, i.e., ISRT (1) is behind
C     all the others and gets drawn first; element ISRT (NEL) is in front
C     and is most visible.
C     Elements with lowest zplane have lowqest sorted index.
C     Given the same zplane, elements with lowest xcen+ycen are lowest.
C
      INTEGER IND(NELM)
      EQUIVALENCE (IND,ZDEPTH(1,3))
C
      DO 1 IE=1,NEL
         ZDEPTH(IE,1) = 100.0*ZCEN(IE) + XCEN(IE)+YCEN(IE)
         ISRT(IE) = IE
1     CONTINUE
C
      CALL SORT(ZDEPTH,IND,NEL)
      CALL ISWAP(ISRT,ZDEPTH,IND,NEL)
C
      RETURN
      END
      SUBROUTINE FILLPF(Icolor)
C     SELECTS FILL panel
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
c
      CALL PRS('Input color:$')
      CALL REI(ic)
      CALL XSETFILL(IC)
c
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE DRAWEL_color(IEL,icolor)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
      INCLUDE 'basics.inc'
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      IIEL=IABS(IEL)
C
      IFSPHR=.FALSE.
      DO 7 IEDGE=1,8
         IF (CCURVE(IEDGE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         XYZCTR(1,1)=X(IE,1)+X(IE,4)+X(IE,5)+X(IE,8)
         XYZCTR(1,2)=X(IE,2)+X(IE,3)+X(IE,6)+X(IE,7)
         XYZCTR(1,3)=X(IE,1)+X(IE,2)+X(IE,5)+X(IE,6)
         XYZCTR(1,4)=X(IE,3)+X(IE,4)+X(IE,7)+X(IE,8)
         XYZCTR(1,5)=X(IE,1)+X(IE,2)+X(IE,3)+X(IE,4)
         XYZCTR(1,6)=X(IE,5)+X(IE,6)+X(IE,7)+X(IE,8)
         XYZCTR(2,1)=Y(IE,1)+Y(IE,4)+Y(IE,5)+Y(IE,8)
         XYZCTR(2,2)=Y(IE,2)+Y(IE,3)+Y(IE,6)+Y(IE,7)
         XYZCTR(2,3)=Y(IE,1)+Y(IE,2)+Y(IE,5)+Y(IE,6)
         XYZCTR(2,4)=Y(IE,3)+Y(IE,4)+Y(IE,7)+Y(IE,8)
         XYZCTR(2,5)=Y(IE,1)+Y(IE,2)+Y(IE,3)+Y(IE,4)
         XYZCTR(2,6)=Y(IE,5)+Y(IE,6)+Y(IE,7)+Y(IE,8)
         XYZCTR(3,1)=Z(IE,1)+Z(IE,4)+Z(IE,5)+Z(IE,8)
         XYZCTR(3,2)=Z(IE,2)+Z(IE,3)+Z(IE,6)+Z(IE,7)
         XYZCTR(3,3)=Z(IE,1)+Z(IE,2)+Z(IE,5)+Z(IE,6)
         XYZCTR(3,4)=Z(IE,3)+Z(IE,4)+Z(IE,7)+Z(IE,8)
         XYZCTR(3,5)=Z(IE,1)+Z(IE,2)+Z(IE,3)+Z(IE,4)
         XYZCTR(3,6)=Z(IE,5)+Z(IE,6)+Z(IE,7)+Z(IE,8)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            IF (J.EQ.1) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.2) THEN
               CALL BEGINB(XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL ENDP
            ENDIF
            IF (J.EQ.3) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL ENDP
            ENDIF
            IF (J.EQ.4) THEN
               CALL BEGINB(XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.5) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL ENDP
            ENDIF
            IF (J.EQ.6) THEN
               CALL BEGINB(XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
   72    CONTINUE
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
C              call fillp(icolor)
               i15 = mod(iel,15)+1
               call fillp(icolor)
C               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8
      IF(IEL.GT.10000)THEN
         IIEL=IIEL-10000
         CALL MOVEC(x(Iiel,IC),y(Iiel,IC))
      ELSE
         CALL BEGINB(x(Iiel,IC),y(Iiel,IC))
      ENDIF
      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          YCENTER=YCENTER+Y(IIEL,IC)/4.
          XCENTER=XCENTER+X(IIEL,IC)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          IF(CCURVE(IEDGE,IIEL).EQ.' ')THEN
             CALL DRAWC(X(IIEL,IC),Y(IIEL,IC))
          ELSE
C            Draw curved side
             CALL GETPTS(NPOINT,CSPACE,IIEL,IEDGE,XCRVED,YCRVED)
             DO 118 I=1,NPOINT
                CALL DRAWC(XCRVED(I),YCRVED(I))
118          CONTINUE
          ENDIF
6     CONTINUE
      IF(IEL.LT.10000)CALL ENDP
C
      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
         WRITE(STRING,'(I3)')IEL
         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
      call color(1)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE HMT_DRAWEL_color(IEL,icolor)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
      INCLUDE 'basics.inc'
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      IIEL=IABS(IEL)
C
      IFSPHR=.FALSE.
      DO 7 IEDGE=1,8
         IF (CCURVE(IEDGE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         XYZCTR(1,1)=X(IE,1)+X(IE,4)+X(IE,5)+X(IE,8)
         XYZCTR(1,2)=X(IE,2)+X(IE,3)+X(IE,6)+X(IE,7)
         XYZCTR(1,3)=X(IE,1)+X(IE,2)+X(IE,5)+X(IE,6)
         XYZCTR(1,4)=X(IE,3)+X(IE,4)+X(IE,7)+X(IE,8)
         XYZCTR(1,5)=X(IE,1)+X(IE,2)+X(IE,3)+X(IE,4)
         XYZCTR(1,6)=X(IE,5)+X(IE,6)+X(IE,7)+X(IE,8)
         XYZCTR(2,1)=Y(IE,1)+Y(IE,4)+Y(IE,5)+Y(IE,8)
         XYZCTR(2,2)=Y(IE,2)+Y(IE,3)+Y(IE,6)+Y(IE,7)
         XYZCTR(2,3)=Y(IE,1)+Y(IE,2)+Y(IE,5)+Y(IE,6)
         XYZCTR(2,4)=Y(IE,3)+Y(IE,4)+Y(IE,7)+Y(IE,8)
         XYZCTR(2,5)=Y(IE,1)+Y(IE,2)+Y(IE,3)+Y(IE,4)
         XYZCTR(2,6)=Y(IE,5)+Y(IE,6)+Y(IE,7)+Y(IE,8)
         XYZCTR(3,1)=Z(IE,1)+Z(IE,4)+Z(IE,5)+Z(IE,8)
         XYZCTR(3,2)=Z(IE,2)+Z(IE,3)+Z(IE,6)+Z(IE,7)
         XYZCTR(3,3)=Z(IE,1)+Z(IE,2)+Z(IE,5)+Z(IE,6)
         XYZCTR(3,4)=Z(IE,3)+Z(IE,4)+Z(IE,7)+Z(IE,8)
         XYZCTR(3,5)=Z(IE,1)+Z(IE,2)+Z(IE,3)+Z(IE,4)
         XYZCTR(3,6)=Z(IE,5)+Z(IE,6)+Z(IE,7)+Z(IE,8)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            IF (J.EQ.1) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.2) THEN
               CALL BEGINB(XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL ENDP
            ENDIF
            IF (J.EQ.3) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL ENDP
            ENDIF
            IF (J.EQ.4) THEN
               CALL BEGINB(XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.5) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL ENDP
            ENDIF
            IF (J.EQ.6) THEN
               CALL BEGINB(XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
   72    CONTINUE
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
C              call fillp(icolor)
               i15 = mod(iel,15)+1
C HMT               call fillp(icolor)
               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8
      IF(IEL.GT.10000)THEN
         IIEL=IIEL-10000
         CALL MOVEC(x(Iiel,IC),y(Iiel,IC))
      ELSE
         CALL BEGINB(x(Iiel,IC),y(Iiel,IC))
      ENDIF
      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          YCENTER=YCENTER+Y(IIEL,IC)/4.
          XCENTER=XCENTER+X(IIEL,IC)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          IF(CCURVE(IEDGE,IIEL).EQ.' ')THEN
             CALL DRAWC(X(IIEL,IC),Y(IIEL,IC))
          ELSE
C            Draw curved side
             CALL GETPTS(NPOINT,CSPACE,IIEL,IEDGE,XCRVED,YCRVED)
             DO 118 I=1,NPOINT
                CALL DRAWC(XCRVED(I),YCRVED(I))
118          CONTINUE
          ENDIF
6     CONTINUE
      IF(IEL.LT.10000)CALL ENDP
C
C HMT
C      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
C HMT
C         WRITE(STRING,'(I3)')IEL
C         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
C         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
C         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
C      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
      call color(1)
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE HMT_DRAWEL_color_nln(IEL,icolor)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
      integer icolor
      INCLUDE 'basics.inc'
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      IIEL=IABS(IEL)
C
      IFSPHR=.FALSE.
      DO 7 IEDGE=1,8
         IF (CCURVE(IEDGE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         XYZCTR(1,1)=X(IE,1)+X(IE,4)+X(IE,5)+X(IE,8)
         XYZCTR(1,2)=X(IE,2)+X(IE,3)+X(IE,6)+X(IE,7)
         XYZCTR(1,3)=X(IE,1)+X(IE,2)+X(IE,5)+X(IE,6)
         XYZCTR(1,4)=X(IE,3)+X(IE,4)+X(IE,7)+X(IE,8)
         XYZCTR(1,5)=X(IE,1)+X(IE,2)+X(IE,3)+X(IE,4)
         XYZCTR(1,6)=X(IE,5)+X(IE,6)+X(IE,7)+X(IE,8)
         XYZCTR(2,1)=Y(IE,1)+Y(IE,4)+Y(IE,5)+Y(IE,8)
         XYZCTR(2,2)=Y(IE,2)+Y(IE,3)+Y(IE,6)+Y(IE,7)
         XYZCTR(2,3)=Y(IE,1)+Y(IE,2)+Y(IE,5)+Y(IE,6)
         XYZCTR(2,4)=Y(IE,3)+Y(IE,4)+Y(IE,7)+Y(IE,8)
         XYZCTR(2,5)=Y(IE,1)+Y(IE,2)+Y(IE,3)+Y(IE,4)
         XYZCTR(2,6)=Y(IE,5)+Y(IE,6)+Y(IE,7)+Y(IE,8)
         XYZCTR(3,1)=Z(IE,1)+Z(IE,4)+Z(IE,5)+Z(IE,8)
         XYZCTR(3,2)=Z(IE,2)+Z(IE,3)+Z(IE,6)+Z(IE,7)
         XYZCTR(3,3)=Z(IE,1)+Z(IE,2)+Z(IE,5)+Z(IE,6)
         XYZCTR(3,4)=Z(IE,3)+Z(IE,4)+Z(IE,7)+Z(IE,8)
         XYZCTR(3,5)=Z(IE,1)+Z(IE,2)+Z(IE,3)+Z(IE,4)
         XYZCTR(3,6)=Z(IE,5)+Z(IE,6)+Z(IE,7)+Z(IE,8)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            IF (J.EQ.1) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.2) THEN
               CALL BEGINB(XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL ENDP
            ENDIF
            IF (J.EQ.3) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL ENDP
            ENDIF
            IF (J.EQ.4) THEN
               CALL BEGINB(XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
            IF (J.EQ.5) THEN
               CALL BEGINB(XISO(X(IE,1),Y(IE,1),Z(IE,1)) 
     $                    ,YISO(X(IE,1),Y(IE,1),Z(IE,1)))
               CALL DRAWC (XISO(X(IE,2),Y(IE,2),Z(IE,2)) 
     $                    ,YISO(X(IE,2),Y(IE,2),Z(IE,2)))
               CALL DRAWC (XISO(X(IE,3),Y(IE,3),Z(IE,3)) 
     $                    ,YISO(X(IE,3),Y(IE,3),Z(IE,3)))
               CALL DRAWC (XISO(X(IE,4),Y(IE,4),Z(IE,4)) 
     $                    ,YISO(X(IE,4),Y(IE,4),Z(IE,4)))
               CALL ENDP
            ENDIF
            IF (J.EQ.6) THEN
               CALL BEGINB(XISO(X(IE,5),Y(IE,5),Z(IE,5)) 
     $                    ,YISO(X(IE,5),Y(IE,5),Z(IE,5)))
               CALL DRAWC (XISO(X(IE,6),Y(IE,6),Z(IE,6)) 
     $                    ,YISO(X(IE,6),Y(IE,6),Z(IE,6)))
               CALL DRAWC (XISO(X(IE,7),Y(IE,7),Z(IE,7)) 
     $                    ,YISO(X(IE,7),Y(IE,7),Z(IE,7)))
               CALL DRAWC (XISO(X(IE,8),Y(IE,8),Z(IE,8)) 
     $                    ,YISO(X(IE,8),Y(IE,8),Z(IE,8)))
               CALL ENDP
            ENDIF
   72    CONTINUE
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
C HMT Looking for outline
C         call color(icolor)
         call hmt_color(icolor)
         IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
c              call fillp(icolor)
               i15 = mod(iel,15)+1
C HMT HACK - found where he fills elms
C               write (6,*) 'icolor = ',icolor
C               write (6,*) 'i15 = ',i15
               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8
      IF(IEL.GT.10000)THEN
         IIEL=IIEL-10000
         CALL MOVEC(x(Iiel,IC),y(Iiel,IC))
      ELSE
         CALL BEGINB(x(Iiel,IC),y(Iiel,IC))
      ENDIF
      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          YCENTER=YCENTER+Y(IIEL,IC)/4.
          XCENTER=XCENTER+X(IIEL,IC)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          IF(CCURVE(IEDGE,IIEL).EQ.' ')THEN
             CALL DRAWC(X(IIEL,IC),Y(IIEL,IC))
          ELSE
C            Draw curved side
             CALL GETPTS(NPOINT,CSPACE,IIEL,IEDGE,XCRVED,YCRVED)
             DO 118 I=1,NPOINT
                CALL DRAWC(XCRVED(I),YCRVED(I))
118          CONTINUE
          ENDIF
6     CONTINUE
      IF(IEL.LT.10000)CALL ENDP
C
      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
C HMT Turn numbers off
C         WRITE(STRING,'(I3)')IEL
C         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
C         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
C         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
C      call color(1)
      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine flipel(ieg,fplane)
c
c     This routine flips about x, y, or z plane
c
      INCLUDE 'basics.inc'
      character key,string*6
      character*1 fplane
      character*3 cbt
c
      integer efac(6)
      save    efac
      data    efac  / 4,2,1,3,5,6 /
c
      logical ifswap(8)
c
      i1 = efac(3)
      i2 = efac(4)
      if (if3d) then
         i1 = efac(5)
         i2 = efac(6)
      endif
c
      DO IFLD=0,MAXFLD
         cbt              = CBC(i1,ieg,IFLD)
         CBC(i1,ieg,IFLD) = CBC(i2,ieg,IFLD)
         CBC(i2,ieg,IFLD) = cbt
         if (CBC(i1,ieg,IFLD).eq.'P  ') CBC(i1,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'P  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'E  ') CBC(i1,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'E  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'SP ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'J  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'SP ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'J  ') CBC(i2,ieg,IFLD) = '   '
         do j=1,5
            Bt                = BC(J,i1,ieg,IFLD)
            BC(J,i1,ieg,IFLD) = BC(J,i2,ieg,IFLD)
            BC(J,i2,ieg,IFLD) = bt
         enddo
      enddo
C
      DO II=1,4
         st              =sides(ieg,i1,II)
         sides(ieg,i1,II)=sides(ieg,i2,II)
         sides(ieg,i2,II)=st
      enddo
C
      if (if3d) then
         do ic=1,4
            i4 = ic+4
            xt       =x(ieg,ic)
            x(ieg,ic)=x(ieg,I4)
            x(ieg,I4)=xt
            yt       =y(ieg,ic)
            y(ieg,ic)=y(ieg,I4)
            y(ieg,I4)=yt
            zt       =z(ieg,ic)
            z(ieg,ic)=z(ieg,I4)
            z(ieg,I4)=zt
         enddo
      else
         xt      =x(ieg,1)
         x(ieg,1)=x(ieg,4)
         x(ieg,4)=xt
         xt      =x(ieg,2)
         x(ieg,2)=x(ieg,3)
         x(ieg,3)=xt
         yt      =y(ieg,1)
         y(ieg,1)=y(ieg,4)
         y(ieg,4)=yt
         yt      =y(ieg,2)
         y(ieg,2)=y(ieg,3)
         y(ieg,3)=yt
      endif
c
      if (fplane.eq.'x') then
         idir = 1
         xcen(ieg) = -xcen(ieg)
         do ic=1,8
            x(ieg,ic) = -x(ieg,ic)
         enddo
      elseif (fplane.eq.'y') then
         idir = 2
         ycen(ieg) = -ycen(ieg)
         do ic=1,8
            y(ieg,ic) = -y(ieg,ic)
         enddo
      else
         idir = 3
         zcen(ieg) = -zcen(ieg)
         do ic=1,8
            z(ieg,ic) = -z(ieg,ic)
         enddo
      endif
C
C     Correct NCOND
C
      IF(MASKEL(ieg,1).EQ.0 .AND. MASKEL(ieg,1).EQ.1)NCOND=NCOND-1
      IF(MASKEL(ieg,1).EQ.1 .AND. MASKEL(ieg,1).EQ.0)NCOND=NCOND+1
      DO 500 IF=1,MAXFLD
         MASKEL(ieg,IF)=MASKEL(ieg,IF)
  500 CONTINUE
c     ISRT  (ieg) = ISRT (ieg)
c     ICRV  (ieg) = ICRV (ieg)
c     NUMAPT(ieg) =NUMAPT(ieg)
c     LETAPT(ieg) =LETAPT(ieg)
      if (if3d) then
         do iface=1,6
            if (ccurve(iface,ieg).eq.'s') 
     $      curve(idir,iface,ieg)=-curve(idir,iface,ieg)
         enddo
         if (ccurve(5,ieg).eq.'s'.or.ccurve(6,ieg).eq.'s') then
            cbt           = ccurve(6,ieg)
            ccurve(6,ieg) = ccurve(5,ieg)
            ccurve(5,ieg) = cbt
            do i=1,5
               ct             = curve(i,6,ieg)
               curve(i,6,ieg) = curve(i,5,ieg)
               curve(i,5,ieg) = ct
            enddo
         endif
c
         do iedge=1,8
            ifswap(iedge)=.false.
         enddo
c
         do iedge=1,4
            iedg4=iedge+4
            if (ccurve(iedge,ieg).eq.'C'.and. .not.ifswap(iedge)) then
c
               ifswap(iedge) = .true.
               ifswap(iedg4) = .true.
c
               cbt           = ccurve(6,ieg)
               ccurve(6,ieg) = ccurve(5,ieg)
               ccurve(5,ieg) = cbt
               do i=1,5
                  ct                 = curve(i,iedg4,ieg)
                  curve(i,iedg4,ieg) = curve(i,iedge,ieg)
                  curve(i,iedge,ieg) = ct
               enddo
            endif
         enddo
c
         do iedge=5,8
            iedg4=iedge-4
            if (ccurve(iedge,ieg).eq.'C'.and. .not.ifswap(iedge)) then
c
               ifswap(iedge) = .true.
               ifswap(iedg4) = .true.
c
               cbt           = ccurve(6,ieg)
               ccurve(6,ieg) = ccurve(5,ieg)
               ccurve(5,ieg) = cbt
               do i=1,5
                  ct                 = curve(i,iedg4,ieg)
                  curve(i,iedg4,ieg) = curve(i,iedge,ieg)
                  curve(i,iedge,ieg) = ct
               enddo
            endif
         enddo
      else   ! 2D:  swap 1 & 3
         cbt           = ccurve(1,ieg)
         ccurve(1,ieg) = ccurve(3,ieg)
         ccurve(3,ieg) = cbt
c
         do i=1,5
            ct             = curve(i,1,ieg)
            curve(i,1,ieg) = curve(i,3,ieg)
            curve(i,3,ieg) = ct
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine reflect_mesh
c
c     Reflect about x-, y-, or z-plane
c
      INCLUDE 'basics.inc'
      character*1 fplane
c
1     CONTINUE
      ITEM(1)='BUILD MENU'
      ITEM(2)='REFLECT ABOUT X=0'
      ITEM(3)='REFLECT ABOUT Y=0'
      ITEM(4)='REFLECT ABOUT Z=0'
      NCHOIC=3
      if (if3d) NCHOIC=4
c
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'Reflect Mesh')
c
3013  CONTINUE
      IF (CHOICE.EQ.'BUILD MENU') return
      IF (CHOICE.EQ.'REFLECT ABOUT X=0') fplane = 'x'
      IF (CHOICE.EQ.'REFLECT ABOUT Y=0') fplane = 'y'
      IF (CHOICE.EQ.'REFLECT ABOUT Z=0') fplane = 'z'
c
      do ieg=1,nel
c        write(6,*) fplane,'  flip: ',ieg
         call flipel(ieg,fplane)
      enddo
c
      goto 1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine substitute_el(xyzbox,nelold)
C     Delete elements in xyzbox, provied e < nelold + 1
      INCLUDE 'basics.inc'
      real xyzbox(6)
c
      common /cdell/ ifkeep(nelm)
      logical ifkeep
      integer e
C
      do e=1,nel
         ifkeep(e) = .true.
      enddo
c
      do e=1,nelold
         if (if3d) then
            if (xyzbox(1).le.xcen(e).and.xcen(e).le.xyzbox(2) .and.
     $          xyzbox(3).le.ycen(e).and.ycen(e).le.xyzbox(4) .and.
     $          xyzbox(5).le.zcen(e).and.zcen(e).le.xyzbox(6)) then
                ifkeep(e)=.false.
            endif
         else
            if (xyzbox(1).le.xcen(e).and.xcen(e).le.xyzbox(2)  .and.
     $          xyzbox(3).le.ycen(e).and.ycen(e).le.xyzbox(4)) then
                ifkeep(e)=.false.
            endif
         endif
c        if (ifkeep(e)) write(6,*) 'KEEP:',e
         if (.not.ifkeep(e)) write(6,*) 'NO KEEP:',e,if3d
     $                                   ,xcen(e),ycen(e),zcen(e)
      enddo
c
c     Now, shift all elements down, to compress out deleted elements
c
      jskip = 0
      je    = 0
      do e=1,nel
         if (ifkeep(e)) then
            je = je+1
            if (je.lt.e) then
               call copyel(e,je)
C              If the last element was a conduction element, reduce NCOND
               if(maskel(je,1).eq.0) ncond=ncond-1
            endif
         else
            write(s,10) e,nel
   10       format(' Deleting element',I9,' of',I9,'.$')
            call prs(s)
         endif
      enddo
      nel = je
C
C     Recount the number of curved sides
C
      ncurve=0
      do 100 ie=1,nel
      do 100 iedge=1,8
         if(ccurve(iedge,ie).ne.' ') ncurve=ncurve+1
  100 continue
C
      return
      end
c-----------------------------------------------------------------------
