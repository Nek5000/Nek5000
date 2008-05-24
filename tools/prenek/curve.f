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
      SUBROUTINE CURVES
      INCLUDE 'basics.inc'
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
      LOGICAL IFTMP
C
      IFTMP=IFGRID
      IFGRID=.FALSE.
      DO 30 IEL=1,NEL
         DO 30 IEDGE=1,NEDGES
            IP1 = IEDGE+1
            IF(IP1.EQ.5) IP1=1
            IF(IP1.EQ.9) IP1=5
            EDGES(IEDGE,1,IEL) = (X(IEL,IEDGE)+X(IEL,IP1))/2.0
            EDGES(IEDGE,2,IEL) = (Y(IEL,IEDGE)+Y(IEL,IP1))/2.0
            EDGES(IEDGE,3,IEL) = (Z(IEL,IEDGE)+Z(IEL,IP1))/2.0
30    CONTINUE
1     CONTINUE
      ITEM(1)='BUILD MENU'
      ITEM(2)='STRAIGHTEN CURVE'
      ITEM(3)='MAKE CIRCLE'
      ITEM(4)='J''ADOUBE'
      ITEM(5)='Tile with hexagons'
c     ITEM(6)='TRANSITION HEXAGONS'
      ITEM(6)='Refine hexagons'
c     ITEM(4)='MAKE SPLINE'
c     ITEM(5)='MAKE SINE WAVE'
C     ITEM(6)='FORTRAN FUNCTION'
      NCHOIC=6
      if (if3d) then
         nchoic=nchoic+1
         ITEM(nchoic)='SPHERICAL MESH'
         nchoic=nchoic+1
         ITEM(nchoic)='STRETCH SPHERE'
         nchoic=nchoic+1
c        ITEM(nchoic)='GEN TEST MESH'
         ITEM(nchoic)='NEAR-WALL SPHERE'
c     else
c        nchoic=nchoic+1
c        item(nchoic)='CIRC MESH'
      endif
      nchoic=nchoic+1
      ITEM(nchoic)='RENUMBER ELEMENTS'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'CURVE SIDES')
3013  CONTINUE
      IF(CHOICE.EQ.'MAKE SPLINE')THEN
          CALL GETEDG(ISID,IELS)
C         GET 2 Auxiliary control points
          CALL COLOR(1)
          CALL xdot(X(IELS,ISID),Y(IELS,ISID))
          CALL PRS(
     $    'Enter First Control point (adjacent to hilighted corner):$')
          CALL MOUSE(XMOUS1,YMOUS1,BUTTON)
          CALL COLOR(8)
          CALL xdot(XMOUS1,YMOUS1)
          CALL PRS(
     $    'Enter Second Control point (adjacent to other corner):$')
          CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
          CALL COLOR(8)
          CALL xdot(XMOUSE,YMOUSE)
          CALL COLOR(1)
          CALL DRAWEL(-IELS)
          CCURVE (ISID,IELS)='S'
          CURVE(1,ISID,IELS)=XMOUS1
          CURVE(2,ISID,IELS)=YMOUS1
          CURVE(3,ISID,IELS)=XMOUSE
          CURVE(4,ISID,IELS)=YMOUSE
C         Export copies curve ICURVE to all overlapping sides
          CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'RENUMBER ELEMENTS') THEN
           CALL RENUM
           GOTO 1
      ELSE IF(CHOICE.EQ.'GEN TEST MESH') THEN
           CALL GENMESH
           GOTO 1
      ELSE IF(CHOICE.EQ.'J''ADOUBE')THEN
           CALL VERTADJ
           GOTO 1
      ELSE IF(CHOICE.EQ.'Refine hexagons') THEN
           call hexagon_refine
           GOTO 1
      ELSE IF(CHOICE.EQ.'TRANSITION HEXAGONS') THEN
           call hex_transition_3d
           GOTO 1 
      ELSE IF(CHOICE.EQ.'Tile with hexagons') THEN
           call hexagon_tile
           GOTO 1
      ELSE IF(CHOICE.EQ.'STRETCH SPHERE') THEN
           call stretch_sphere
           GOTO 1
c     ELSE IF(CHOICE.EQ.'CIRC MESH') THEN
c          call circ_mesh 
c          GOTO 1
      ELSE IF(CHOICE.EQ.'SPHERICAL MESH') THEN
           CALL SPHMESH
           GOTO 1
      ELSE IF(CHOICE.EQ.'NEAR-WALL SPHERE') THEN
           call nw_sphmesh
           goto 1
      ELSE IF(CHOICE.EQ.'STRAIGHTEN CURVE')THEN
          CALL GETEDG(ISID,IELS)
          CALL DRAWEL(-IELS)
          DO 46 I=1,6
             CURVE(I,ISID,IELS)=0.0
46        CONTINUE
          CCURVE (ISID,IELS)=' '
          WRITE(s,'(1X,A14,I3,I4)')' STRAIGHTENING',ISID,IELS
          CALL PRS(S//'$')
C         Export copies curve ICURVE to all overlapping sides
          CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'MAKE CIRCLE')THEN
            CALL GETEDG(ISID,IELS)
            XCENT=(X(IELS,1)+X(IELS,2)+X(IELS,3)+X(IELS,4))/4.0
            YCENT=(Y(IELS,1)+Y(IELS,2)+Y(IELS,3)+Y(IELS,4))/4.0
            CALL GWRITE(XCENT,YCENT,1.0,'*$')
            CALL PRS(' Enter Circle Radius for side of element *$')
            CALL PRS(' >0 for convex element, <0 for concave$')
            CALL KEYPADm(RADIUS,button,xmouse,ymouse)
            if (button.eq.'RIGHT') then
c              clicking near an element 
               call getside(jel,jside,xmouse,ymouse)
               RADIUS =  -CURVE(1,jside,jel)
            endif
            CALL DRAWEL(-IELS)
            CURVE(1,ISID,IELS)=RADIUS
            CCURVE(ISID,IELS)='C'
C
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'FORTRAN FUNCTION')THEN
            CALL GETEDG(ISID,IELS)
            CALL DRAWEL(-IELS)
            CCURVE(ISID,IELS)='F'
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'MAKE SINE WAVE')THEN
            CALL GETEDG(ISID,IELS)
            CALL PRS(' **WARNING** sine wave cannot be rotated.$')
            CALL PRS(' **WARNING** sine wave cannot be refined.$')
            CALL PRS('It must be of the form:   Y=Yo + A*sin(alpha x)$')
            CALL PRS
     $      ('First Enter the amplitude  A: (Negative to abort)$')
            CALL KEYPAD(CURVE(1,ISID,IELS))
            IF(CURVE(1,ISID,IELS).LE.0.0)GO TO 1
            CALL PRS('Enter the wave number  alpha:$')
            CALL KEYPAD(CURVE(2,ISID,IELS))
            CALL PRS('Enter the Y-offset Yo:$')
            CALL KEYPAD(CURVE(3,ISID,IELS))
            CALL DRAWEL(-IELS)
            CCURVE(ISID,IELS)='W'
C
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'BUILD MENU')THEN
            NCURVE=0
            DO 114 IE=1,NEL
            DO 114 IEDGE=1,8
               IF (CCURVE(IEDGE,IE).NE.' ') THEN
                  NCURVE=NCURVE+1
               ENDIF
114         CONTINUE
            IFGRID=IFTMP
            RETURN
      ENDIF
C     After modifiacation of floor, check if there are any elements in floors
C     above.  If not, modify ceiling also.
      IF(.NOT.IFCEIL .AND.IF3D)THEN
         DO 11 I=1,IEL
            IF(NUMAPT(I).GT.ILEVEL) THEN
               CALL PRS(' Presence of elements on higher '//
     $         'floors inhibits modifying ceiling also$')
               RETURN
            ENDIF
11       CONTINUE
C        Since we got here, there must be no elements on upper floors
C        Copy modifications.
         DO 12 I=1,6
             CURVE(I,ISID+4,IELS)=CURVE(I,ISID,IELS)
12        CONTINUE
          CCURVE(ISID+4,IELS)=CCURVE(ISID,IELS)
C         Export copies curve ICURVE to all overlapping sides
          CALL EXPORT(ISID+4,IELS)
      ENDIF
      GO TO 1
      END
      SUBROUTINE MESGEN
      INCLUDE 'basics.inc'
      DIMENSION DXCHEB(20),DYCHEB(20),XCHEB(20),YCHEB(20)
      REAL CSPACE(100),XCRVED(100),YCRVED(100),XFRAC(20)
      REAL XPTSEL(NXM,NYM),YPTSEL(NXM,NYM)
C     Make mesh and draw on screen
      CALL PRSIS(' Generating Mesh Points for $',NEL,' Elements$')
      IF(NEL .GT. NELM-2 ) CALL PRS(' ***** ERROR ******* $')
      IF(NEL .GT. NELM-2 ) CALL PRSIS('EXCEEDED MAXIMUM OF$',NELM,
     $' ELEMENTS$')
      call color(8)
      IF (NX.GT.NXM) THEN
         WRITE(S,104) NX,NXM
  104    FORMAT(' Warning, current N exceeds NXM, resetting from',I3
     $         ,' to',I3,'.$')
         CALL PRS(S)
         NX=NXM
         NY=NYM
         NZ=NZM
      ENDIF
      DO 105 I=1,NX
         XFRAC(I)=XXPTS(I)/2.0 + 0.5
105   CONTINUE
      DO 108 IEL = 1,NEL
         DO 106 I = 1,NX
           DO 106 J = 1,NY
C           Calculate original xpts based on straight sides
            XPTS(I,J,IEL) =  1./4.* (
     $      X(IEL,1)*(XXPTS(I)-1.0)*(YYPTS(J)-1.0)-
     $      X(IEL,2)*(XXPTS(I)+1.0)*(YYPTS(J)-1.0)+
     $      X(IEL,3)*(XXPTS(I)+1.0)*(YYPTS(J)+1.0)-
     $      X(IEL,4)*(XXPTS(I)-1.0)*(YYPTS(J)+1.0))
            YPTS(I,J,IEL) =  1./4.* (
     $      Y(IEL,1)*(XXPTS(I)-1.0)*(YYPTS(J)-1.0)-
     $      Y(IEL,2)*(XXPTS(I)+1.0)*(YYPTS(J)-1.0)+
     $      Y(IEL,3)*(XXPTS(I)+1.0)*(YYPTS(J)+1.0)-
     $      Y(IEL,4)*(XXPTS(I)-1.0)*(YYPTS(J)+1.0))
C           Initial values for temporary array
            XPTSEL(I,J) = XPTS(I,J,IEL)
            YPTSEL(I,J) = YPTS(I,J,IEL)
106       CONTINUE
           ICRV(IEL)=0
C          For 2-d only
           DO 165 IE=1,4
             IF(CCURVE(IE,IEL).NE.' ')ICRV(IEL)=ICRV(IEL)+IE**2
165        CONTINUE
           IF(ICRV(IEL).NE.0)THEN
C             Shift internal points to account for curviness
              DO 107 ISID=1,4
                IF(CCURVE(ISID,IEL).NE.' ')THEN
                  NPOINT=NX
                  CALL GETPTS(NPOINT,XFRAC,IEL,ISID,XCRVED,YCRVED)
                  DO 117 I = 1,NX
                    IF(ISID.EQ.1)THEN
                      II=I
                      JJ=1
                    ELSE IF(ISID.EQ.2)THEN
                      II=NX
                      JJ=I
                    ELSE IF(ISID.EQ.3)THEN
                      II=NX-I+1
                      JJ=NX
                    ELSE IF(ISID.EQ.4)THEN
                      II=1
                      JJ=NX-I+1
                     ENDIF
                    DXCHEB(I)=XCRVED(I)-XPTS(II,JJ,IEL)
                    DYCHEB(I)=YCRVED(I)-YPTS(II,JJ,IEL)
117               CONTINUE
                  DO 120 I=1,NX
                    DO 120 J=1,NY
                      IF(ISID.EQ.2)III=J
                      IF(ISID.EQ.4)III=NX-J+1
                      IF(ISID.EQ.1)III=I
                      IF(ISID.EQ.3)III=NX-I+1
                      IF(ISID.EQ.1) FAC = XFRAC(NX-J+1)
                      IF(ISID.EQ.2) FAC = XFRAC(I)
                      IF(ISID.EQ.3) FAC = XFRAC(J)
                      IF(ISID.EQ.4) FAC = XFRAC(NX-I+1)
C
                      XPTSEL(I,J)=XPTSEL(I,J)+DXCHEB(III)*FAC
                      YPTSEL(I,J)=YPTSEL(I,J)+DYCHEB(III)*FAC
120               CONTINUE
                ENDIF
107           CONTINUE
C             Now that the points in the temporary array are moved, copy to xpts
              DO 109 I=1,NX
                DO 109 J=1,NY
                  XPTS(I,J,IEL)=XPTSEL(I,J)
                  YPTS(I,J,IEL)=YPTSEL(I,J)
109           CONTINUE
            endif
108       CONTINUE
C
        DXMIN=1.0E10
        DYMIN=1.0E10
        DO 119 IEL=1,NEL
              DX=SQRT((XPTS(2,1,IEL)-XPTS(1,1,IEL))**2
     $               +(YPTS(2,1,IEL)-YPTS(1,1,IEL))**2)
              DY=SQRT((XPTS(1,2,IEL)-XPTS(1,1,IEL))**2
     $               +(YPTS(1,2,IEL)-YPTS(1,1,IEL))**2)
              DXMIN=AMIN1(DXMIN,DX)
              DYMIN=AMIN1(DYMIN,DY)
119      CONTINUE
C
        CALL HARD
        IFHARD=.TRUE.
        call color(8)
        DO 112 IEL = 1,NEL
C         Draw only those on this floor
          IF(NUMAPT(IEL).EQ.ILEVEL)THEN
            IF(ICRV(IEL).EQ.0)THEN
               DO 209 I = 1,NX
                 CALL MOVEC(XPTS(I,1 ,IEL),YPTS(I,1 ,IEL))
                 CALL DRAWC(XPTS(I,NY,IEL),YPTS(I,NY,IEL))
209            CONTINUE
               DO 211 J = 1,NY
                 CALL MOVEC(XPTS(1 ,J,IEL),YPTS(1 ,J,IEL))
                 CALL DRAWC(XPTS(NX,J,IEL),YPTS(NX,J,IEL))
211            CONTINUE
            ELSE
C              Draw it the hard, slow, curvy way
               DO 1109 J=1,NY-1
                   DO 1109 I = 1,NX
                     IF(CCURVE(2,IEL).NE.' '.AND. I.EQ.NX) GO TO 1108
                     IF(CCURVE(4,IEL).NE.' '.AND. I.EQ. 1) GO TO 1108
                     CALL MOVEC(XPTS(I,J  ,IEL),YPTS(I,J  ,IEL))
                     CALL DRAWC(XPTS(I,J+1,IEL),YPTS(I,J+1,IEL))
1108                CONTINUE
1109            CONTINUE
               DO 111 I=1,NX-1
                   DO 111 J = 1,NY
                     IF(CCURVE(1,IEL).NE.' '.AND. J.EQ. 1) GO TO 110
                     IF(CCURVE(3,IEL).NE.' '.AND. J.EQ.NY) GO TO 110
                     CALL MOVEC(XPTS(I  ,J,IEL),YPTS(I  ,J,IEL))
                     CALL DRAWC(XPTS(I+1,J,IEL),YPTS(I+1,J,IEL))
110                  CONTINUE
111            CONTINUE
C              Draw curved sides for hardcopy
               DO 130 IEDGE=1,4
                  IF(CCURVE(IEDGE,IEL).NE.' ')THEN
                     CALL MOVEC(X(IEL,IEDGE),Y(IEL,IEDGE))
                     CALL DRAWED(IEL,IEDGE,1)
                  ENDIF
130            CONTINUE
            ENDIF
          ENDIF
112      continue
         CALL NOHARD
        IFHARD=.FALSE.
         NCURVE=0
         DO 114 IEL=1,NEL
            DO 114 IEDGE=1,8
               IF(CCURVE(IEDGE,IEL).NE.' ') NCURVE=NCURVE+1
114      CONTINUE
         RETURN
      END
      SUBROUTINE GETEDG(ISID,IELS)
      INCLUDE 'basics.inc'
C
1     CALL PRS(' Enter element side.$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      XH=XMOUSE
      YH=YMOUSE
      IF(XSCR(XMOUSE) .GT.1.0)THEN
         CALL PRS
     $   (' **ERROR**: Expecting you to define side.  Define side$')
         CALL PRS(' by hitting mouse button in mesh area.  Try again.$')
         go to 1
      ENDIF
C     Here IF(IFCEIL) we are on the ceiling;else on the bottom of a given floor.
      ZH=0.0
      IF(ILEVEL.GT.1)THEN
         DO 10 I=1,ILEVEL-1
            ZH=ZH+HEIGHT(I)
10       CONTINUE
      ENDIF
      IF(IFCEIL)ZH=ZH+HEIGHT(ILEVEL)
      RMIN=1.0E10
      DO 50 IEL=1,NEL
         IF(IFCEIL)THEN
            II1=5
            II2=8
         ELSE
            II1=1
            II2=4
         ENDIF
         DO 40 IEDGE=II1,II2
C         DO 40 IEDGE=1,nedges
            DXH=XH-EDGES(IEDGE,1,IEL)
            DYH=YH-EDGES(IEDGE,2,IEL)
            DZH=ZH-EDGES(IEDGE,3,IEL)
            R=SQRT(DXH**2+DYH**2+DZH**2)
            IF(R.LT.RMIN) THEN
               RMIN=R
               IELS=IEL
               ISID=IEDGE
            ENDIF
40       CONTINUE
50    CONTINUE
      RETURN
      END
      SUBROUTINE EXPORT(IEDGE,IEL)
      INCLUDE 'basics.inc'
      DIMENSION IOVER(2,20)
C     Takes characteristics of curve ICURVE and copies it to any overlapping
C     edges
C
C     Use apartment letters to find if there are any downstairs ceilings
C     What was wrong with going up to roof level?
C     First Check for overlapping side
C
C     Find all elements and sides that get changed
C     Put original side in IOVER 1
      NOVER = 1
      IOVER(1,NOVER)=IEDGE
      IOVER(2,NOVER)=IEL
      EPSI=(EDGES(1,1,IEL)-EDGES(3,1,IEL))**2
     $    +(EDGES(1,2,IEL)-EDGES(3,2,IEL))**2
     $    +(EDGES(1,3,IEL)-EDGES(3,3,IEL))**2
      DO 10 IELS=1,NEL
         EPSJ=(EDGES(1,1,IELS)-EDGES(3,1,IELS))**2
     $       +(EDGES(1,2,IELS)-EDGES(3,2,IELS))**2
     $       +(EDGES(1,3,IELS)-EDGES(3,3,IELS))**2
         DO 10 IEDG=1,NEDGES
            DXH=EDGES(IEDGE,1,IEL)-EDGES(IEDG,1,IELS)
            DYH=EDGES(IEDGE,2,IEL)-EDGES(IEDG,2,IELS)
            DZH=EDGES(IEDGE,3,IEL)-EDGES(IEDG,3,IELS)
            R=DXH**2+DYH**2+DZH**2
            EPS=0.001*MIN(EPSI,EPSJ)
            IF(R.LT.EPS) THEN
               IF(IEDG.NE.IEDGE .OR. IEL.NE.IELS)THEN
                  NOVER = NOVER+1
                  IOVER(1,NOVER)=IEDG
                  IOVER(2,NOVER)=IELS
      write(s,11) nover,iedg,iels,iedge,iel,r,eps,epsi,epsj
11    format('jiedgel,r,epsij:',5i3,4e10.3,'$')
      call prs(s)
               ENDIF
            ENDIF
10    CONTINUE
      IF(NOVER.GE.2)THEN
C        Do this only for overlapping elements on this level, not original
         DO 20 I=2,NOVER
            IELS=IOVER(2,I)
            IF(NUMAPT(IELS).EQ.NUMAPT(IEL)) CALL DRAWEL(-IELS)
20       CONTINUE
         DO 30 I=2,NOVER
            IEDG=IOVER(1,I)
            IELS=IOVER(2,I)
            DO 25 II=1,6
               CURVE(II,IEDG,IELS)=CURVE(II,IEDGE,IEL)
25          CONTINUE
            CCURVE (IEDG,IELS)=CCURVE( IEDGE,IEL)
C           Convex/concave for circle (except on different levels)
C            IF(NUMAPT(IEL).EQ.NUMAPT(IELS))THEN
            IF(LETAPT(IEL).NE.LETAPT(IELS))THEN
               IF(CCURVE(IEDG,IELS).EQ.'C')
     $         CURVE(1,IEDG,IELS)=-CURVE(1,IEDGE,IEL)
               IF(CCURVE(IEDG,IELS).EQ.'S')THEN
C              Spline on same level switch auxiliary points
                  CURVE(1,IEDG,IELS)=CURVE(3,IEDGE,IEL)
                  CURVE(2,IEDG,IELS)=CURVE(4,IEDGE,IEL)
                  CURVE(3,IEDG,IELS)=CURVE(1,IEDGE,IEL)
                  CURVE(4,IEDG,IELS)=CURVE(2,IEDGE,IEL)
               ENDIF
            ENDIF
30       CONTINUE
      ENDIF
      write(s,11) nover,iedg,iels,iedge,iel,r,eps,epsi,epsj
      call prs(s)
      DO 40 I=1,NOVER
         IELS=IOVER(2,I)
         IF(NUMAPT(IELS).EQ.NUMAPT(IEL)) CALL DRAWEL(IELS)
40    CONTINUE
      RETURN
      END
      subroutine getside(jel,jside,xp,yp)
      INCLUDE 'basics.inc'
C     Find closest element, side that we want to duplicate
      RMIN=1.e22
      DO 133 IIEL=1,NEL
C        Only try those on same floor
         IF(NUMAPT(IIEL).EQ.ILEVEL)THEN
            DO 131 IISIDE=1,NSIDES
               R=    ((XP-SIDES(IIEL,IISIDE,1))**2
     $         +      (YP-SIDES(IIEL,IISIDE,2))**2)
               IF (R.LT.RMIN) THEN
                  RMIN=R
                  jside = IISIDE
                  jel   = IIEL
               ENDIF
131         CONTINUE
         ENDIF
133   CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine stretch_sphere
c
c     Stretch all points in spherical shell
c
      INCLUDE 'basics.inc'
C
      call prs(' Input sphere center (x,y,z):$')
      call rerrr(x0,y0,z0)
c
      call prs(' Input inner cut-off radius, r0:$')
      call rer(r0)
c
      call prs(' Input outer cut-off radius, r1:$')
      call rer(r1)
c
      call prs(' Input expansion factor ( =< 0 to abort):$')
      call rer(sfact)
      if (sfact.le.0) return
C
C     Stretch for hemisphere problem
C
C     1)  For all elements with z < z0, just stretch in X-Y
C
C     2)  For all elements with z > z0, stretch in X-Y and Z
C
C
C     Take care of pts first
C
      nvts=8
c     sfact = 1.25
c     z0 = 0.0346
c     rr = 0.5*0.5 + .0001
c     re = 0.5 + .0001
c
      z00= -9.e9
      z0 = 0.
      rr = 0.
      re = 0.
C
      r02 = r0*r0
      r12 = r1*r1
      do ie=1,nel
      do i=1,nvts
         xs = x(ie,i)-x0
         ys = y(ie,i)-y0
         zs = z(ie,i)-z0
         r2 = xs*xs + ys*ys + zs*zs
         if (r02.lt.r2 .and. r2 .lt. r12) then
            xs = sfact*xs
            ys = sfact*ys
            zs = sfact*zs
c
            x(ie,i) = x0 + xs
            y(ie,i) = y0 + ys
            z(ie,i) = z0 + zs
         endif
      enddo
      enddo
C
C     Take care of curved sides
C
      eps = .01*r1
      do ie = 1,nel
      do is = 1,8
c
         if (ccurve(is,ie).eq.'s') then
c           spherical side, x0  = curve(1,is,ie)
c           spherical side, y0  = curve(2,is,ie)
c           spherical side, z0  = curve(3,is,ie)
c           spherical side, rad = curve(4,is,ie)
c
            xs  = curve(1,is,ie)-x0
            ys  = curve(2,is,ie)-y0
            zs  = curve(3,is,ie)-z0
            re2 = xs*xs + ys*ys + zs*zs
            rad = curve(4,is,ie)
            if (re2.lt.eps .and. r0 .lt.rad .and. rad.le. r1) 
     $         curve(4,is,ie) = sfact*curve(4,is,ie)
         endif
c
c        if (ccurve(is,ie).eq.'C') then
c           cylindrical side, rad = curve(1,is,ie)
c           if (abs(curve(1,is,ie)).gt.re)
c    $         curve(1,is,ie) = sfact*curve(1,is,ie)
c        endif
c
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine find_long_edges ! Find long edge of a 2D element
c
      include 'basics.inc'
      integer e,en,et
c
      do e=1,nel
         long_edge(e) = long_edge_e(e)
      enddo
      return
      end
c-----------------------------------------------------------------------
      function long_edge_e(e) ! Find long edge of a 2D element
c
      include 'basics.inc'
      integer e,en,et
c
c
      dx = x(e,1)-x(e,4)
      dy = y(e,1)-y(e,4)
      dm = dx*dx + dy*dy
      long_edge_e = 4
c
      do i=1,3
         dx = x(e,i+1)-x(e,i)
         dy = y(e,i+1)-y(e,i)
         d2 = dx*dx + dy*dy
         if (d2.gt.dm) then
            long_edge_e = i
            dm = d2
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_tile
c
c     Tile a 2D region w/ hexagons (2D only, pff 10/10/05)
c
      include 'basics.inc'
      real hexscale(8,2),hexbase(8,2)
      save hexscale,hexbase
c
      integer e
      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) then
         icalld = icalld + 1
         one = 1.
         pi  = 4.*atan(one)
         pi3 = pi/3.
         call rzero(hexbase,16)
         hexbase(1,1) = -2.
c
         do i=0,5
            angle = pi + i*pi3
            k     = 2+i
            hexbase(k,1) = cos(angle)
            hexbase(k,2) = sin(angle)
         enddo
         hexbase(8,1) = hexbase(2,1) ! close loop
         hexbase(8,2) = hexbase(2,2)
      endif
      call outmat(hexbase(2,1),1,8,'hxb x',1)
      call outmat(hexbase(2,2),1,8,'hxb y',1)
c
C
      CALL PRS('Enter (w/mouse) 2 points for framing box.$')

c     IFTMP =IFGRID
c     IFGRID=.FALSE.
 
      call prs('Enter 1st point:$')
         call mouse(xmouse,ymouse,button)
         xmouse = min(xmouse,xphy(1.))
         xmouse = max(xmouse,xphy(0.))
         ymouse = min(ymouse,yphy(1.))
         ymouse = max(ymouse,yphy(0.))

      call prs('Enter 2nd point:$')
         call mouse(xmous2,ymous2,button)
         xmous2 = min(xmous2,xphy(1.))
         xmous2 = max(xmous2,xphy(0.))
         ymous2 = min(ymous2,yphy(1.))
         ymous2 = max(ymous2,yphy(0.))
C     
C     We successfully input 2 points in the build area
C     1) count number of centroids in the bounding box
C     2) if ncntrd=0, find nearest element and anihilate it.
C     
C     Box
      xmax=max(xmouse,xmous2)
      xmin=min(xmouse,xmous2)
      ymax=max(ymouse,ymous2)
      ymin=min(ymouse,ymous2)
      call drwbox(xmin,ymin,xmax,ymax,1)
c
      dxhex = 3*griddx
      dyhex = 3.
      dyhex = sqrt(dyhex)*griddx
      imin = xmin/dxhex
      imax = xmax/dxhex
      jmin = ymin/dyhex
      jmax = ymax/dyhex
      e = nel
      write(6,*) imin,imax,jmin,jmax,' imnx',dxhex
      write(6,*) xmin,xmax,ymin,ymax,' xmnx',dyhex
      do j=jmin,jmax
      do i=imin,imax
         xshift = 0
         yshift = 0
         do k=1,2
            xi = i*dxhex + xshift
            yj = j*dyhex + yshift
            call cmult2(hexscale,hexbase,griddx,16)
            call cadd  (hexscale(1,1),xi,8)
            call cadd  (hexscale(1,2),yj,8)
c
            e = e + 1
            do iv=1,4
               x(e,iv) = hexscale(1+iv,1) ! kth hex, lower half
               y(e,iv) = hexscale(1+iv,2)
            enddo
            long_edge(e) = 3
            call drawel(e)
c
            e = e + 1
            do iv=1,4
               x(e,iv) = hexscale(4+iv,1) ! kth hex, upper half
               y(e,iv) = hexscale(4+iv,2)
            enddo
            long_edge(e) = 1
            call drawel(e)
c
            xshift = dxhex/2
            yshift = dyhex/2
         enddo
      enddo
      enddo
      nel = e
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine hex_transition_3d_e(el,en)  
c
c     Refine 3D hex prisms in transition zone, i.e. the zone between  
c     coarse and fine 3D hex meshes
c
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)

c     Element el is partitioned into 5 elements: a, b, c, d, e. Element a
c     replaces original element el. 
c
      parameter (m3=-3,m2=-2,n2=2,n3=3)
      real    rsvx(-3:3) 
      integer rptrx(8,5), rptry(8,5), rptrz(8,5) ! 8 pts, 5 elmts, 3 dir.
      save rptr
c
c     Coordinates for each of a, b, c, d and e (in x, y and z) 
      data rptrx / m2,n2,n3,m3,m2,n2,n3,m3
     $            ,-1,m2,m3,-1,-1,m2,m3,-1
     $            ,m3,n3, 1,-1,m3,n3, 1,-1
     $            ,n2, 1, 1,n3,n2, 1, 1,n3
     $            ,-1, 1, 1,-1,m2,n2,n3,m3  /
c      
      data rptry / -1,-1, 0, 0,-1,-1, 0, 0
     $            ,-1,-1, 0, 1,-1,-1, 0, 1
     $            , 0, 0, 1, 1, 0, 0, 1, 1
     $            ,-1,-1, 1, 0,-1,-1, 1, 0
     $            ,-1,-1, 1, 1,-1,-1, 0, 0  /
c
      data rptrz /  0, 0,n3,n3, 1, 1, 1, 1 
     $            ,-1, 0,n3,-1, 1, 1, 1, 1
     $            ,n3,n3,-1,-1, 1, 1, 1, 1
     $            , 0,-1,-1,n3, 1, 1, 1, 1
     $            ,-1,-1,-1,-1, 0, 0,n3,n3  /

      integer el,en,et
c
      rsvx( 0) = 0
      do i=1,3
        rsvx( i) =  1./i
        rsvx(-i) = -1./i
      end do
c
      do k=0,3
        call copyel(el,en+k)
      enddo
c
      nx = nxm
      ny = nym
      nz = nzm
      call genxyz (xp,yp,zp,el,1,nx,ny,nz)
c
      rrl(3) = 0 
      do k=1,5    ! 5 new elements
        et = en + (k-2)
        if (k.eq.1) et = el
        do i=1,8  ! 8 points
          rrl(1) = rsvx(rptrx(i,k))
          rrl(2) = rsvx(rptry(i,k))
          rrl(3) = rsvx(rptrz(i,k))
          call evalsc (x(et,i),xp,rrl,1)
          call evalsc (y(et,i),yp,rrl,0)
          call evalsc (z(et,i),zp,rrl,0)
        enddo

        if (k.eq.3) then
           ccurve(1,et) = ' '
           ccurve(5,et) = ' '
        endif

        call drawel(et)
      enddo
c
      nel = nel + 4
c
      return
      end
c-----------------------------------------------------------------------              
      subroutine hex_transition_3d ! Refine 3D hexes in transition zone
c
      include 'basics.inc'
      integer e,en,et
c
      integer icalld
      save icalld
      data icalld  /0/
c
      call find_long_edges
      call gencen
c
      nelo = nel
      do e=1,nelo
         en = nel+1
         call rotate_element(e,long_edge(e))
         call hex_transition_3d_e(e,en) ! refine 3D hexes in a zone
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_refine_e(e,en) ! refine 2D hexes in a zone
c
      include 'basics.inc'
      common /chexb/ hexscale(8,2),hexbase(8,2)
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
c
      parameter (m3=-3,m2=-2,n2=2,n3=3)
      real    rsvx(-3:3)
      integer rptr(4,2,4,4)
      save rptr
      data rptr / m2,n2,n3,m3,-1,-1, 0, 0
     $          , -1,m2,m3,-1,-1,-1, 0, 1
     $          , n2, 1, 1,n3,-1,-1, 1, 0
     $          , m3,n3, 1,-1, 0, 0, 1, 1
c
     $          ,  0, 1, 1, 0,m3,m2,n2,n3
     $          , -1, 1, 1, 0,-1,-1,m2,m3
     $          ,  0, 1, 1,-1,n3,n2, 1, 1
     $          , -1, 0, 0,-1,-1,m3,n3, 1
c
     $          , m3,n3,n2,m2, 0, 0, 1, 1
     $          , n3, 1, 1,n2, 0,-1, 1, 1
     $          , -1,m3,m2,-1,-1, 0, 1, 1
     $          , -1, 1,n3,m3,-1,-1, 0, 0
c
     $          , -1, 0, 0,-1,m2,m3,n3,n2
     $          , -1, 0, 1,-1,n2,n3, 1, 1
     $          , -1, 1, 0,-1,-1,-1,m3,m2
     $          ,  0, 1, 1, 0,m3,-1, 1,n3   /

      integer e,en,et
c
      rsvx( 0) = 0
      do i=1,3
         rsvx( i) =  1./i
         rsvx(-i) = -1./i
      enddo
c
      do k=0,2
         call copyel(e,en+k)
      enddo
c
      if (long_edge(e).eq.1) then
         long_edge(en  ) = 4
         long_edge(en+1) = 2
         long_edge(en+2) = 3
      elseif (long_edge(e).eq.2) then
         long_edge(en  ) = 1
         long_edge(en+1) = 3
         long_edge(en+2) = 4
      elseif (long_edge(e).eq.3) then
         long_edge(en  ) = 2
         long_edge(en+1) = 4
         long_edge(en+2) = 1
      elseif (long_edge(e).eq.4) then
         long_edge(en  ) = 3
         long_edge(en+1) = 1
         long_edge(en+2) = 2
      endif
c
      do j=0,maxfld
         cbc(long_edge(e),en+2,j) = ' '
      enddo
c
      nx = nxm
      ny = nym
      nz = 1
      call genxyz (xp,yp,zp,e,1,nx,ny,nz)
c
      rrl(3) = 0.
      do k=1,4     ! 4 new elements
         et = en + (k-2)
         if (k.eq.1) et = e
         do i=1,4  ! 4 points
            rrl(1) = rsvx(rptr(i,1,k,long_edge(e)))
            rrl(2) = rsvx(rptr(i,2,k,long_edge(e)))
            call evalsc (x(et,i),xp,rrl,1)
            call evalsc (y(et,i),yp,rrl,0)
                        z(et,i) = 0
         enddo
         call drawel(et)
      enddo
      nel = nel + 3
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_refine ! refine 2D hexes in a zone
c
      include 'basics.inc'
      integer e,en,et
c
      integer elist(nelm),marked(nelm)
c
      integer icalld
      save    icalld
      data    icalld  /0/
c
      call get_el2(nlist,elist)
c
c     This is ugly (potentially ON^2), but effective
c
      k  = 0
      kk = 0
   10 continue
c
      call find_long_edges
      call gencen
      call gen_neigh
c
      k  = k+1
      kk = kk+1
      nnew = 0
      do e=1,nel
         if (elist(e).ne.0) then ! pending request
            l  = long_edge(e)
            en = neighb(l,e)
            write(6,*) nnew,en,e,l,' Enew'
            if (en.gt.0) then
               ln = long_edge(en)
               et = neighb(ln,en)
               if (et.eq.e) then ! reciprocal pair, ok to refine
                  elist(e ) = -1
                  elist(en) = -1
               elseif (elist(en).eq.0) then
                  nnew = nnew+1
                  elist(en) = 1
               endif
            else
               elist(e) = -1
            endif
         endif
      enddo
      write(6,*) 'KPASS:',k,nnew,kk
      if (nnew.gt.0) goto 10
      k  = 0
c
      nelo = nel
      do e=1,nelo
         if (elist(e).lt.0) then
            en = nel+1
            call hexagon_refine_e(e,en) ! refine 2D hexes in a zone
            elist(e) = 0
         endif
      enddo
c
      icount = 0
      do e=1,nelo
         if (elist(e).ne.0) icount = icount+1
      enddo
      if (icount.gt.0) goto 10
c
      call redraw_mesh
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_element(e,kface)
c
c     Rotate element e so that face kface is mapped to face 1
c
c     We assume that kface  is in [1,4]
c
      include 'basics.inc'
      integer e,et
c
      if (kface.eq.1) return
      if (kface.lt.1.or.kface.gt.4) then
         call prsii('Cannot rotate element. kface out of range.$'
     $              ,e,kface)
         return
      endif
c
      et = nelm-1
      call copyel(e,et)
c
      do it=1,4
c
         i = it-(kface-1)
         if (i.lt.1) i = i + 4
         x(e,i  ) = x(et,it)
         y(e,i  ) = y(et,it)
         z(e,i  ) = z(et,it)
         x(e,i+4) = x(et,it+4)
         y(e,i+4) = y(et,it+4)
         z(e,i+4) = z(et,it+4)
c
         do ifld=0,maxfld
            cbc(i,e,ifld) = cbc(it,et,ifld)
            call copy(bc(1,i,e,ifld),bc(1,it,et,ifld),5)
         enddo
         do ii=1,4
            sides(e,i,ii)=sides(et,it,ii)
         enddo

         ccurve(i  ,e) = ccurve(it  ,et)
         call copy(curve(1,i,e),curve(1,it,et),5)
 
         ccurve(i+4,e) = ccurve(it+4,et)
         call copy(curve(1,i+4,e),curve(1,it+4,et),5)

         edges (i,1,e) = edges (it,1,et)
         edges (i,2,e) = edges (it,2,et)
         edges (i,3,e) = edges (it,3,et)
      enddo

      do k=5,6
         if (ccurve(k,et).eq.'s') then
            ccurve(k,e) = ccurve(k,et)
            call copy(curve(1,k,e),curve(1,k,et),5)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_el2(nlist,elist)
c
c     Grab a set of elements
c
      include 'basics.inc'
      integer elist(1),e,tlist(nelm)
c
      call get_els(nlist,tlist)
c
c     convert list to flags
c
      call izero  (elist,nelm)
      do i=1,nlist
         e = tlist(i)
         elist(e) = 1
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_els(nlist,elist)
c
c     Grab a set of elements
c
      include 'basics.inc'
      integer elist(1),e
      logical iftmp
C
      CALL PRS('Enter (w/mouse) 2 points for framing box.$')

      IFTMP =IFGRID
      IFGRID=.FALSE.
 
      call prs('Enter 1st point:$')
      call mouse(xmouse,ymouse,button)
      if (xmouse.lt.xphy(1.)) then
         xmouse = min(xmouse,xphy(1.))
         xmouse = max(xmouse,xphy(0.))
         ymouse = min(ymouse,yphy(1.))
         ymouse = max(ymouse,yphy(0.))

         call prs('Enter 2nd point:$')
            call mouse(xmous2,ymous2,button)
            xmous2 = min(xmous2,xphy(1.))
            xmous2 = max(xmous2,xphy(0.))
            ymous2 = min(ymous2,yphy(1.))
            ymous2 = max(ymous2,yphy(0.))
C     
C        We successfully input 2 points in the build area
C        1) count number of centroids in the bounding box
C        2) if ncntrd=0, find nearest element and anihilate it.
C     
C        Box
         xmax=max(xmouse,xmous2)
         xmin=min(xmouse,xmous2)
         ymax=max(ymouse,ymous2)
         ymin=min(ymouse,ymous2)
         call drwbox(xmin,ymin,xmax,ymax,1)
      else
         call prs('Enter point (drop mouse!) (99,99 for file input):$')
         call rerr(xmin,ymin)
         xmax = xmin
         ymax = ymin
         if (ymin.eq.99.and.ymax.eq.99) then
            call get_els_bc(nlist,elist)
            return
         endif
      endif
c
      nlist = 0
      do e=1,nel
         if ( xmin.le.xcen(e) .and. xcen(e).le.xmax .and.
     $        ymin.le.ycen(e) .and. ycen(e).le.ymax ) then
            nlist = nlist + 1
            elist(nlist) = e
         endif
      enddo
      if (nlist.eq.0) then ! find closest element
         xm = .5*(xmax+xmin)
         ym = .5*(ymax+ymin)
         dx = xm-xcen(1)
         dy = ym-ycen(1)
         dm = dx*dx + dy*dy
         im = 1
         do e=1,nel
            dx = xm-xcen(e)
            dy = ym-ycen(e)
            d2 = dx*dx + dy*dy
            if (d2.lt.dm) then
               dm = d2
               im = e
            endif
         enddo
         nlist = 1
         elist(1) = im
      endif
c
      IFGRID=IFTMP
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_els_bc(nlist,elist)
c
c     Grab a set of elements
c
      include 'basics.inc'
      integer elist(1),e
c
      ifld = 2
      if (ifflow) ifld = 1
c
      nlist = 0
      do e=1,nel
         do k=1,4
            if (cbc(k,e,ifld).eq.'v'.or.cbc(k,e,ifld).eq.'t') then
               nlist = nlist+1
               elist(nlist) = e
               goto 10
            endif
         enddo
   10    continue
      enddo
c
      return
      end
c-----------------------------------------------------------------------
