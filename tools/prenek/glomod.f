C------------------------------------------------------------------------------
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
c-----------------------------------------------------------------------
      subroutine glomod
      include 'basics.inc'
      include 'basicsp.inc'
      INTEGER ISPLIT(NELM)
      LOGICAL IFTMP
C
      IFTMP=IFGRID
      IFGRID=.FALSE.
1     NCHOIC=1
      ITEM(NCHOIC)='END GLOBAL REFINE'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='ZIPPER'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='Non-conf SPLIT'
c
c     if (.not.if3d) then          ! pff 8/10/05 (make room in menu)
c             NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='SPIDER WEB'
c             NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='PICTURE FRAME'
c             NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='CORNER FRAME'
c           NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='SMOOTH'
c     endif
c
      IF(IF3D)THEN
c        NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='SPLIT FLOOR'
         NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='NEW SPLIT'
c        NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='3D CORNER'
      ENDIF
c
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='REPLICATE/ROTATE'
           NCHOIC=NCHOIC+1
      if (ndim.eq.2) then
         ITEM(NCHOIC)='Refine Hexagons'
              NCHOIC=NCHOIC+1
      else
         ITEM(NCHOIC)='Hex transition'
              NCHOIC=NCHOIC+1
      endif
      ITEM(NCHOIC)='OCT/Multi-SPLIT'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='Mesh Edit'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='SHIFT'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='STRETCH'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='CLIP DOMAIN'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='Clean up vertices'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'GLOBAL REFINE')
C     Update sides array
      CALL MKSIDE
C     UPDATE MAXLET (which could have been changed in the refine operation)
      MAXLET=0
      IF(NEL.GT.0) THEN
         DO 950 IEL=1,NEL
            MAXLET=AMAX0(MAXLET,ICHAR(LETAPT(IEL)))
950      CONTINUE
      ENDIF
C
      IF(CHOICE.EQ.'END GLOBAL REFINE')THEN
         call redraw_mesh
         ifgrid=iftmp
         return
      ELSE IF(CHOICE.EQ.'SPIDER WEB')THEN
         CALL PARKER(ISPLIT)
      ELSE IF(CHOICE.EQ.'Clean up vertices')THEN
         CALL VERTADJ
      ELSE IF(CHOICE.EQ.'CLIP DOMAIN')THEN
         CALL CLPDOM
      ELSE IF(CHOICE.EQ.'STRETCH')THEN
         CALL STRETCH
      ELSE IF(CHOICE.EQ.'Mesh Edit')THEN
         call mesh_edit
      ELSE IF(CHOICE.EQ.'SHIFT')THEN
         CALL SHIFT
      ELSE IF(CHOICE.EQ.'SMOOTH')THEN
         call smooth
      ELSE IF(CHOICE.EQ.'NEW SPLIT')THEN
         CALL ZIP2
      ELSE IF(CHOICE.EQ.'REPLICATE/ROTATE')THEN
         CALL rep_rot
      ELSE IF(CHOICE.EQ.'Hex transition')THEN
         call hex_transition_3d 
      ELSE IF(CHOICE.EQ.'Refine Hexagons')THEN
         call hexagon_refine
      ELSE IF(CHOICE.EQ.'OCT/Multi-SPLIT')THEN
         CALL OCTSPL
      ELSE IF(CHOICE.EQ.'3D CORNER')THEN
         CALL WIND3D
      ELSE IF(CHOICE.EQ.'SPLIT FLOOR')THEN
         CALL SPLITF
      ELSE IF(CHOICE.EQ.'PICTURE FRAME')THEN
         CALL FRAME(ISPLIT,2)
      ELSE IF(CHOICE.EQ.'CORNER FRAME')THEN
         CALL FRAME(ISPLIT,1)
      ELSE IF(CHOICE.EQ.'ZIPPER')THEN
C        Initiate crack in 1st element  
C        Element IELCRK a fraction SFRAC from IFAcrk
         CALL CRACK(IELCRK,IFACRK,SFRAC)
         IF(SFRAC.GT.0.AND.SFRAC.LT.1.0)THEN
C           Mark where CRACK propagates in vector ISPLIT
            CALL MARK (IELCRK,IFACRK,SFRAC,ISPLIT)
C           Split elements marked in vector ISPLIT
            CALL SPLIT(SFRAC,ISPLIT)
         ELSE
            CALL PRS('Aborting zipper operation.$')
         ENDIF
      ELSE IF(CHOICE.EQ.'Non-conf SPLIT')THEN
C        Initiate crack in 1st element  
C        Element IELCRK a fraction SFRAC from IFAcrk
         CALL CRACK(IELCRK,IFACRK,SFRAC)
         IF(SFRAC.GT.0.AND.SFRAC.LT.1.0)THEN
C           Mark where CRACK propagates in vector ISPLIT
            CALL MARK2 (IELCRK,IFACRK,SFRAC,ISPLIT)
C           Split elements marked in vector ISPLIT
            CALL SPLIT(SFRAC,ISPLIT)
         ELSE
            CALL PRS('Aborting zipper operation.$')
         ENDIF
      ENDIF
c
      GO TO 1
      END
c-----------------------------------------------------------------------
      subroutine frame(isplit,nframe)
C     Refines corner elements
      include 'basics.inc'
      DIMENSION IOVER(2),IOSIDE(2),IFEL(2),ICIN(2),ISPLIT(NELM)
      IF(NFRAME.EQ.1)THEN
         CALL PRS('Enter Element to be divided:$')
         CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
         RMIN=1.0E10
         DO 8 IEL=1,NEL
C           Only for 2-d
            RAD=SQRT( (XMOUSE-XCEN(IEL))**2 + (YMOUSE-YCEN(IEL))**2 )
            IF(RAD.LT.RMIN)THEN
               RMIN=RAD
               IFEL(1) = IEL
            ENDIF
8       CONTINUE
         CALL PRS
     $   ('Enter Corner where smaller element is to be created:$')
         CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
         RMIN=1.0E10
         DO 50 IC=1,4
C           Only for 2-d
            RAD=SQRT( (XMOUSE-X(IFEL(1),IC))**2
     $      +         (YMOUSE-Y(IFEL(1),IC))**2 )
            IF(RAD.LT.RMIN)THEN
               RMIN=RAD
               ICIN(1) = IC
            ENDIF
50       CONTINUE
      ELSE IF(NFRAME.EQ.2)THEN
         CALL PRS
     $   ('For a Picture Frame corner mesh refinement, there must$')
         CALL PRS
     $   ('Be two elements which share a side.  Each element must$')
         CALL PRS(
     $   'have a boundary side adjacent as illustrated in the manual$')
         CALL PRS('Enter inside corner of picture frame:$')
         CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
         DO 200 IIEL=1,NFRAME
            RMIN=1.0E10
            DO 150 IC=1,4
              DO 100 IEL=1,NEL
C                Only for 2-d
                 RAD=SQRT((XMOUSE-X(IEL,IC))**2 + (YMOUSE-Y(IEL,ic))**2)
                 IF(RAD.LT.RMIN)THEN
                    IF(IIEL.EQ.1)THEN
C                      Pick any old element without further checking
                       RMIN=RAD
                       IFEL(IIEL) =IEL
                       ICIN(IIEL) =IC
                    ELSE IF(IIEL.EQ.2)THEN
                       IF(IEL.EQ.IFEL(1)) THEN
C                         Skip over element we already found
                       ELSE
C                         Make sure that this candidate for 2nd element
C                         is adjacent to 1st element, then pick it
                          LAP=0
                          DO 90 IS1=1,4
                            DO 90 IS2=1,4
                             ROVER = SQRT(
     $                       (SIDES(IFEL(1),IS1,1)-SIDES(IEL,IS2,1))**2+
     $                       (SIDES(IFEL(1),IS1,2)-SIDES(IEL,IS2,2))**2+
     $                       (SIDES(IFEL(1),IS1,3)-SIDES(IEL,IS2,3))**2)
                             IF(ROVER.LT.XFAC/1000.) THEN
C                               Set overlaps, we found a common side.
                                LAP=1
                                IOVER(1)=IS1
                                IOVER(2)=IS2
                             ENDIF
90                        CONTINUE
                          IF(LAP.EQ.1)THEN
                             RMIN=RAD
                             IFEL(IIEL) =IEL
                             ICIN(IIEL) =IC
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDIF
100           CONTINUE
150         CONTINUE
200      CONTINUE
      ENDIF
      CALL PRS
     $('Enter the Shrink factor for the new elements.  (0<s<1).$')
      CALL PRS
     $('I recommend that these new elements be shrunk to about .3$')
      CALL PRS('<0 to abort. $')
      call rer(fac)
      IF(FAC.LE.0.0)THEN
         CALL PRS('Aborting Picture Frame mesh refinement.$')
         RETURN
      ELSE IF(FAC.GT.1.0)THEN
        CALL PRS('Error; You cannot make new refined elements larger $')
        CALL PRS('than the old ones.  ABORTING.$')
         RETURN
      ENDIF
C     Mark elements on other levels for spider web
      DO 10 I=1,NELM
         ISPLIT(I)=0
10    CONTINUE
      ISPLIT(IFEL(1))=1
      IF(NFRAME.EQ.2)ISPLIT(IFEL(2))=2
      IF(IF3D)THEN
C        Look above and below for emelents to be marked for cracking
         DO 40 Isweep=1,NEL
            DO 30 IEL=1,NEL
               IF(ISPLIT(IEL).NE.0)THEN
C                 Check neighbors of this marked element for splits
                  DO 20 IIEL=1,NEL
                     RADB=SQRT(
     $               (SIDES(IEL,6,1)-SIDES(IIEL,5,1))**2 +
     $               (SIDES(IEL,6,2)-SIDES(IIEL,5,2))**2 +
     $               (SIDES(IEL,6,3)-SIDES(IIEL,5,3))**2 )
                     RADT=SQRT(
     $               (SIDES(IEL,5,1)-SIDES(IIEL,6,1))**2 +
     $               (SIDES(IEL,5,2)-SIDES(IIEL,6,2))**2 +
     $               (SIDES(IEL,5,3)-SIDES(IIEL,6,3))**2 )
                     IF(RADB.LT.XFAC/1000.)ISPLIT(IIEL)=ISPLIT(IEL)
                     IF(RADT.LT.XFAC/1000.)ISPLIT(IIEL)=ISPLIT(IEL)
20                CONTINUE
               ENDIF
30          CONTINUE
40       CONTINUE
      ENDIF
C     Do a few optional checks
C     If an element on top of another is rotated, you may have problems!
      DO 313 IIEL=1,NFRAME
            IF(ICIN(IIEL).EQ.IOVER(IIEL))THEN
C              FOR "TOP" element (inside corner# = overlapping side #)
C              "FREE" side should be overlap-1
               IFREE=IOVER(IIEL) - 1
               IF(IFREE.LT.1)IFREE=IFREE+4
            ELSE
C              FOR "SIDE" element (inside corner# NOT= overlapping side #)
C              "FREE" side should be overlap+1
               IFREE=IOVER(IIEL) + 1
               IF(IFREE.GT.4)IFREE=IFREE-4
            ENDIF
C           Now check if free side is truly free
            DO 310 IEL=1,NEL
              DO 310 IS2=1,4
C               Single element has 2 sides that need checking; Each double has 1
                IF(NFRAME.EQ.1) NCHECK=2
                IF(NFRAME.EQ.2) NCHECK=1
                DO 310 ICHECK=1,NCHECK
                  IF(NFRAME.EQ.1 .AND. ICHECK.EQ.1) IFREE=ICIN(1)
                  IF(NFRAME.EQ.1 .AND. ICHECK.EQ.2) IFREE=ICIN(1)-1
                  IF(IFREE.LT.1)IFREE=IFREE+4
                  ROVER = SQRT(
     $            (SIDES(IFEL(IIEL),IFREE,1)-SIDES(IEL,IS2,1))**2+
     $            (SIDES(IFEL(IIEL),IFREE,2)-SIDES(IEL,IS2,2))**2+
     $            (SIDES(IFEL(IIEL),IFREE,3)-SIDES(IEL,IS2,3))**2)
                  IF(ROVER.LT.XFAC/1000. .AND.IEL.NE.IFEL(IIEL)) THEN
C                   Problem.  "Free" side overlaps another element.
                    CALL PRSIS('***WARNING*** Element$',IFEL(IIEL),
     $              ' is being split $')
                    CALL PRSI ('where it is adjacent to element$',IEL)
                    CALL PRSIS('You must split$',IEL,' to match.$')
                  ENDIF
310         continue
313   CONTINUE
C      DO 400 IIEL=1,NFRAME
       DO 400 IIIEL=1,NEL
         IIEL=ISPLIT(IIIEL)
         IF(IIEL.GT.0)THEN
           IFEL(IIEL)=IIIEL
C          IFEL=ELEMENT#, ICIN=INSIDE CORNER #, IOVER=OVERLAPPING SIDE #
           IC=ICIN(IIEL)
           ICP1=IC+1
           ICP2=IC+2
           ICP3=IC+3
           IF(ICP1.GT.4) ICP1=ICP1-4
           IF(ICP2.GT.4) ICP2=ICP2-4
           IF(ICP3.GT.4) ICP3=ICP3-4
C          Straighten out element to be split
           DO 160 IEDGE=1,8
              CCURVE(IEDGE,IFEL(IIEL))=' '
160        CONTINUE
           CALL COPYEL(IFEL(IIEL),NEL+1)
           CALL COPYEL(IFEL(IIEL),NEL+2)
c
           ILETAP=MAXLET
           IF(IIEL.EQ.2)ILETAP=MAXLET+2
           DO 320 I=1,2
                 ILETAP=ILETAP+1
                 LETAPT(NEL+I)=CHAR(MIN0(ILETAP,122))
320        CONTINUE
         DO 390 ICEILF=1,2
           IF(ICEILF.EQ.2)THEN
               IC  =IC  +4
               ICP1=ICP1+4
               ICP2=ICP2+4
               ICP3=ICP3+4
           ENDIF
           IF(ICIN(IIEL).EQ.IOVER(IIEL))THEN
C            FOR "TOP" element (inside corner# = overlapping side #)
C
             X(NEL+1,IC  )=X(NEL+1,IC)+(X(NEL+1,ICP1)-X(NEL+1,IC))*FAC
             Y(NEL+1,IC  )=Y(NEL+1,IC)+(Y(NEL+1,ICP1)-Y(NEL+1,IC))*FAC
C            Calculate Middle point
             XDIAG=X(NEL+1,IC)
             YDIAG=Y(NEL+1,IC)
             XMID=XDIAG+( X(IFEL(IIEL),ICP2)-X(IFEL(IIEL),ICP1) ) * FAC
             YMID=YDIAG+( Y(IFEL(IIEL),ICP2)-Y(IFEL(IIEL),ICP1) ) * FAC
C
            X(NEL+1,ICP3)=XMID
            Y(NEL+1,ICP3)=YMID
C
            X(NEL+2,IC  )=X(NEL+2,IC)+(X(NEL+2,ICP3)-X(NEL+2,IC))*FAC
            Y(NEL+2,IC  )=Y(NEL+2,IC)+(Y(NEL+2,ICP3)-Y(NEL+2,IC))*FAC
            X(NEL+2,ICP1)=XMID
            Y(NEL+2,ICP1)=YMID
C
            X(IFEL(IIEL),ICP1)=X(IFEL(IIEL),IC)+
     $     (X(IFEL(IIEL),ICP1)-X(IFEL(IIEL),IC))*FAC
            Y(IFEL(IIEL),ICP1)=Y(IFEL(IIEL),IC)+
     $     (Y(IFEL(IIEL),ICP1)-Y(IFEL(IIEL),IC))*FAC
            X(IFEL(IIEL),ICP2)=XMID
            Y(IFEL(IIEL),ICP2)=YMID
            X(IFEL(IIEL),ICP3)=X(IFEL(IIEL),IC)+
     $     (X(IFEL(IIEL),ICP3)-X(IFEL(IIEL),IC))*FAC
            Y(IFEL(IIEL),ICP3)=Y(IFEL(IIEL),IC)+
     $     (Y(IFEL(IIEL),ICP3)-Y(IFEL(IIEL),IC))*FAC
         ELSE
C           FOR "other" element (inside corner# NOT= overlapping side #)
            X(NEL+1,IC  )=X(NEL+1,IC)+(X(NEL+1,ICP3)-X(NEL+1,IC))*FAC
            Y(NEL+1,IC  )=Y(NEL+1,IC)+(Y(NEL+1,ICP3)-Y(NEL+1,IC))*FAC
C
C           Calculate Middle point
            XDIAG=X(NEL+1,IC)
            YDIAG=Y(NEL+1,IC)
            XMID=XDIAG+( X(IFEL(IIEL),ICP2)-X(IFEL(IIEL),ICP3) ) * FAC
            YMID=YDIAG+( Y(IFEL(IIEL),ICP2)-Y(IFEL(IIEL),ICP3) ) * FAC
C
            X(NEL+1,ICP1)=XMID
            Y(NEL+1,ICP1)=YMID
C
            X(NEL+2,IC  )=X(NEL+2,IC)+(X(NEL+2,ICP1)-X(NEL+2,IC))*FAC
            Y(NEL+2,IC  )=Y(NEL+2,IC)+(Y(NEL+2,ICP1)-Y(NEL+2,IC))*FAC
            X(NEL+2,ICP3)=XMID
            Y(NEL+2,ICP3)=YMID
C
            X(IFEL(IIEL),ICP1)=X(IFEL(IIEL),IC)+
     $     (X(IFEL(IIEL),ICP1)-X(IFEL(IIEL),IC))*FAC
            Y(IFEL(IIEL),ICP1)=Y(IFEL(IIEL),IC)+
     $     (Y(IFEL(IIEL),ICP1)-Y(IFEL(IIEL),IC))*FAC
            X(IFEL(IIEL),ICP2)=XMID
            Y(IFEL(IIEL),ICP2)=YMID
            X(IFEL(IIEL),ICP3)=X(IFEL(IIEL),IC)+
     $     (X(IFEL(IIEL),ICP3)-X(IFEL(IIEL),IC))*FAC
            Y(IFEL(IIEL),ICP3)=Y(IFEL(IIEL),IC)+
     $     (Y(IFEL(IIEL),ICP3)-Y(IFEL(IIEL),IC))*FAC
         ENDIF
390      CONTINUE
C
         NEL=NEL+2
       ENDIF
400   CONTINUE
C     Recalculate centers
      IF(NFRAME.EQ.1)NBACK=1
      IF(NFRAME.EQ.2)NBACK=3
      DO 110 IEL=1,nel
        XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
        YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
        CALL DRAWEL(IEL)
110   CONTINUE
      DO 120 IIEL=1,NFRAME
        IEL=IFEL(IIEL)
        XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
        YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
        CALL DRAWEL(IEL)
120   CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine splitf
C     Splits floor in mesh
      include 'basics.inc'
      CALL PRS('Which floor do you wish to split?$')
      call rer(floor)
      ISPLIT=FLOOR
      IF(ISPLIT.GT.NLEVEL) THEN
         CALL PRS
     $   ('ERROR- ONLY ',NLEVEL,'Floors Exist.  Aborting Split.$')
         RETURN
      ENDIF
      CALL PRS('At what height in this floor do you want the crack?$')
      CALL PRS('0<h<1; 0 is at floor, 1 is at ceiling.'//
     $'  negative to abort.$')
      call rer(fac)
      IF(FAC.GE.1.0)THEN
         CALL PRS('ERROR- Split cant be above ceiling$')
         RETURN
      ELSE IF (FAC.LE.0.0)then
         CALL PRS('Aborting Floor Split$')
         RETURN
      ENDIF
C
      NELOLD=NEL
      NLEVEL=NLEVEL+1
      ILEVEL=ILEVEL+1
      CALL DRELEV(ILEVEL,ILEVEL-1,'     ')
      ILETAPT=MAXLET
      DO 10 IEL=1,NELOLD
         IF(NUMAPT(IEL).EQ.ISPLIT)THEN
           NEL=NEL+1
           CALL COPYEL(IEL,NEL)
           DO 5 IC=1,4
              X(IEL,IC+4)=X(IEL,IC  ) + (X(IEL,IC+4)-X(IEL,IC  )) * FAC
              X(NEL,IC  )=X(IEL,IC+4)
              Y(IEL,IC+4)=Y(IEL,IC  ) + (Y(IEL,IC+4)-Y(IEL,IC  )) * FAC
              Y(NEL,IC  )=Y(IEL,IC+4)
              Z(IEL,IC+4)=Z(IEL,IC  ) + (Z(IEL,IC+4)-Z(IEL,IC  )) * FAC
              Z(NEL,IC  )=Z(IEL,IC+4)
C             Now interpolate curved sides
              DO 4 II=1,6
                  CURVE(II,IC+4,IEL)=CURVE(II,IC  ,IEL) +
     $           (CURVE(II,IC+4,IEL)-CURVE(II,IC  ,IEL))* FAC
                  CURVE(II,IC  ,NEL)=CURVE(II,IC+4,IEL)
4             CONTINUE
C             Make new line straight unless both old ones were circles
C             or both old ones were splines
              IF( (CCURVE(IC,NEL).EQ.'C'.AND.CCURVE(IC+4,NEL).EQ.'C')
     $        .OR.(CCURVE(IC,NEL).EQ.'S'.AND.CCURVE(IC+4,NEL).EQ.'S') )
     $        THEN
C                Leave them alone
              ELSE
C                Straighten them
                 CCURVE(IC+4,IEL)=' '
              ENDIF
              CCURVE(IC  ,NEL)=CCURVE(IC+4,IEL)
C
5          continue
           NUMAPT(NEL)=ISPLIT+1
           XCEN(NEL)=(X(NEL,1)+X(NEL,2)+X(NEL,3)+X(NEL,4))/4.
           YCEN(NEL)=(Y(NEL,1)+Y(NEL,2)+Y(NEL,3)+Y(NEL,4))/4.
         ELSE IF(NUMAPT(IEL).GT.ISPLIT)THEN
              NUMAPT(IEL)=NUMAPT(IEL)+1
         ENDIF
10    CONTINUE
      CALL SORTEL
C     Draw whole new isometric surface
      DO 51 IEL=1,NEL
c        CALL DRAWIS(ISRT(IEL))
         IF(.NOT.IF3D  .OR. NUMAPT(IEL).EQ.ILEVEL) CALL DRAWEL(IEL)
51    CONTINUE
C     Fix Height Vector
      DO 80 I=NLEVEL,1,-1
         IF(I.GT.ISPLIT+1)HEIGHT(I)=HEIGHT(I-1)
         IF(I.EQ.ISPLIT+1)HEIGHT(I)=HEIGHT(ISPLIT)*(1.0-FAC)
         IF(I.EQ.ISPLIT  )HEIGHT(I)=HEIGHT(ISPLIT)*(    FAC)
80    continue
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine parker(isplit)
C     Makes spider web
      include 'basics.inc'
      INTEGER ISPLIT(NELM)
C
      CALL PRS('Push left button to enter element to be split:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C        apparently is trying to type
         CALL PRS('Type X and Y coordinates:$')
         call rerr(XMOUSE,YMOUSE)
      ENDIF
      RMIN=1.0E10
      DO 100 IEL=1,NEL
         RAD=SQRT( (XMOUSE-XCEN(IEL))**2 + (YMOUSE-YCEN(IEL))**2 )
         IF(RAD.LT.RMIN .AND. NUMAPT(IEL).EQ.ILEVEL)THEN
            RMIN=RAD
            IELCRK=IEL
         ENDIF
100   CONTINUE
C
C     Mark elements on other levels for spider web
      DO 10 I=1,NELM
         ISPLIT(I)=0
10    CONTINUE
      ISPLIT(IELCRK)=1
      IF(IF3D)THEN
C        Look above and below for emelents to be marked for cracking
         DO 40 Isweep=1,NEL
            DO 30 IEL=1,NEL
               IF(ISPLIT(IEL).NE.0)THEN
C                 Check neighbors of this marked element for splits
                  DO 20 IIEL=1,NEL
                     RADB=SQRT(
     $               (SIDES(IEL,6,1)-SIDES(IIEL,5,1))**2 +
     $               (SIDES(IEL,6,2)-SIDES(IIEL,5,2))**2 +
     $               (SIDES(IEL,6,3)-SIDES(IIEL,5,3))**2 )
                     RADT=SQRT(
     $               (SIDES(IEL,5,1)-SIDES(IIEL,6,1))**2 +
     $               (SIDES(IEL,5,2)-SIDES(IIEL,6,2))**2 +
     $               (SIDES(IEL,5,3)-SIDES(IIEL,6,3))**2 )
                     IF(RADB.LT.XFAC/1000.)ISPLIT(IIEL)=1
                     IF(RADT.LT.XFAC/1000.)ISPLIT(IIEL)=1
20                CONTINUE
               ENDIF
30          CONTINUE
40       CONTINUE
      ENDIF
      NELOLD=NEL
      DO 500 ICHECK=1,NELOLD
         IF(ISPLIT(ICHECK).NE.0)THEN
            IELCRK=ICHECK
            CALL COPYEL(IELCRK,NEL+1)
            CALL COPYEL(IELCRK,NEL+2)
            CALL COPYEL(IELCRK,NEL+3)
            CALL COPYEL(IELCRK,NEL+4)
            ILETAP=MAXLET
            DO 320 I=1,4
                 ILETAP=ILETAP+1
                 LETAPT(NEL+I)=CHAR(MIN0(ILETAP,122))
320         CONTINUE
            IF(IF3D)THEN
              XCENT=(X(IELCRK,5)+X(IELCRK,6)+X(IELCRK,7)+X(IELCRK,8))/4
              YCENT=(Y(IELCRK,5)+Y(IELCRK,6)+Y(IELCRK,7)+Y(IELCRK,8))/4
            ENDIF
            DO 200 IP=1,4
C              Straighten internal curved sides
               DO 160 IEDGE=1,8
                  CCURVE(IEDGE,IELCRK)=' '
                  IF(IEDGE.NE.IP .AND. IEDGE.NE.IP+4)
     $            CCURVE(IEDGE,NEL+IP)=' '
                  iiiel=nel+ip
160            CONTINUE
               DO 200 INEW=1,2
C              Element adjacent to side IP
               IC=IP+INEW+1
               IF(IC.GT.4)IC=IC-4
               IF(IC.LE.0)IC=IC+4
               IF(INEW.EQ.1)IOUT=IC+3
               IF(INEW.EQ.2)IOUT=IC+1
               IF(IOUT.GT.4)IOUT=IOUT-4
               X    (NEL+IP,IC)=X(NEL+IP,IOUT)+
     $         (XCEN(NEL+IP) -  X(NEL+IP,IOUT))*.65
               Y    (NEL+IP,IC)=Y(NEL+IP,IOUT)+
     $         (YCEN(NEL+IP) -  Y(NEL+IP,IOUT))*.65
               IF(IF3D)THEN
                  X(NEL+IP,IC+4)=X(NEL+IP,IOUT+4)+
     $            (XCENT    -    X(NEL+IP,IOUT+4))*.65
                  Y(NEL+IP,IC+4)=Y(NEL+IP,IOUT+4)+
     $            (YCENT    -    Y(NEL+IP,IOUT+4))*.65
               ENDIF
               IF(INEW.EQ.2)THEN
C                 Move Corners of center element
                  X(IELCRK,IP)=X(NEL+IP,IC)
                  Y(IELCRK,IP)=Y(NEL+IP,IC)
                  IF(IF3D)THEN
                     X(IELCRK,IP+4)=X(NEL+IP,IC+4)
                     Y(IELCRK,IP+4)=Y(NEL+IP,IC+4)
                  ENDIF
               ENDIF
200         CONTINUE
C           Recalculate centers
            DO 110 IEL=NEL+1,nel+4
              XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
              YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
110         CONTINUE
            IEL=IELCRK
            XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
            YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
C
            NEL=NEL+4
         ENDIF
500   CONTINUE
C     Draw whole new isometric surface
      CALL SORTEL
      DO 51 IEL=1,NEL
c        CALL DRAWIS(ISRT(IEL))
         IF(.NOT.IF3D  .OR. NUMAPT(IEL).EQ.ILEVEL) CALL DRAWEL(IEL)
51    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine crack(ielcrk,ifacrk,sfrac)
C     Initiates crack in 1st element of split.  Output is Element #, Side#, and
C     Fraction away from side.
      include 'basics.inc'
      DIMENSION IELMOV(NELM),IELAB(2),ICORAB(2),ICOMON(2,2,16)
      REAL XCHECK(8),YCHECK(8),XAB(2),YAB(2)
      CHARACTER KEY,STRING*6,MOVPAR*9
      logicAL CHGLOB
C
C
C     FIND CLOSEST ELEMENT AND CORNER
      SFRAC=0.5
      CALL PRS('Push left button to enter element to be split:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C        apparently trying to enter
         CALL PRS('Type X-Y coordinates:$')
         call rerr(XMOUSE,YMOUSE)
      ENDIF
      RMIN=1.0E10
      DO 100 IEL=1,NEL
         RAD=SQRT( (XMOUSE-XCEN(IEL))**2 + (YMOUSE-YCEN(IEL))**2 )
C        Only look at elements on this floor
         IF(RAD.LT.RMIN.AND.(.NOT.IF3D .OR. NUMAPT(IEL).EQ.ILEVEL))THEN
            RMIN=RAD
            IELCRK=IEL
         ENDIF
100   CONTINUE
      CALL PRS('Enter a side of this element that is to be split$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C        apparently trying to enter
         CALL PRS('Type X-Y coordinates:$')
         call rerr(XMOUSE,YMOUSE)
      ENDIF
      RMIN=1.0E10
      DO 2 IS=1,4
C        We're inside the element
         RAD=SQRT( (XMOUSE-SIDES(IELCRK,IS,1))**2
     $    +        (YMOUSE-SIDES(IELCRK,IS,2))**2)
         IF(RAD.LT.RMIN)THEN
            RMIN=RAD
            ISS=IS
         ENDIF
2     CONTINUE
9999  CONTINUE
C     Find adjacent sides
      IFACRK=ISS+1
      IF(IFACRK.EQ.5)IFACRK=1
      IFBCRK=ISS-1
      IF(IFBCRK.EQ.0)IFBCRK=4
      CALL COLOR(1)
      CALL GWRITE(SIDES(IELCRK,IFACRK,1),SIDES(IELCRK,IFACRK,2),1.,'a$')
      CALL GWRITE(SIDES(IELCRK,IFBCRK,1),SIDES(IELCRK,IFBCRK,2),1.,'b$')
c
      CALL PRS('Enter S, the location of split.  0<S<1 .  S < 0.5 is $')
      CALL PRS(
     $'Closer to a; S > 0.5 is closer to b. <0 to abort.$')
      call rer(sfrac)
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine mark (ielcrk,ifacrk,sfrac,isplit)
C     Mark where CRACK propagates in vector ISPLIT
      include 'basics.inc'
      INTEGER ISPLIT(NELM)
C
      DO 10 I=1,NELM
         ISPLIT(I)=0
10    CONTINUE
      ISPLIT(IELCRK)=IFACRK
C     Now Check each element in mesh to see if it is adjacent to a side with
C     a crack on it.  If it is, enter it into the ISPLIT array.  If it is
C     already split, make sure the new split wouldn't conflict with old.
C     After going through all elements, repeat until no new splits occur
C     after checking all elements.
      DO 300 IELCHK=1,NEL
         NEWSPL=0
         DO 200 MARKED=1,NEL
            IF(ISPLIT(MARKED).NE.0)THEN
               DO 100 IEL=1,NEL
C                 Run through faces already cracked
                  DO 90 IFMAR=1,NSIDES
C                    Skip marked face and its opposite
                     ISFACE=ISPLIT(MARKED)
                     IF(ISFACE.EQ.1) IOPP  = 3
                     IF(ISFACE.EQ.2) IOPP  = 4
                     IF(ISFACE.EQ.3) IOPP  = 1
                     IF(ISFACE.EQ.4) IOPP  = 2
                     IF(ISFACE.EQ.5) IOPP  = 6
                     IF(ISFACE.EQ.6) IOPP  = 5
                     IF(IFMAR.NE.ISFACE .AND. IFMAR .NE. IOPP) THEN
C                     Check against all new faces
                      DO 80 IFNEW=1,NSIDES
C                      Check if any new faces overlap any split faces
                       IF(ABS(SIDES(IEL,IFNEW,1)-SIDES(MARKED,IFMAR,1))
     $                 .LT. XFAC/1000.0    .AND.
     $                 ABS(SIDES(IEL,IFNEW,2)-SIDES(MARKED,IFMAR,2))
     $                 .LT. YFAC/1000.0    .AND.
     $                 ABS(SIDES(IEL,IFNEW,3)-SIDES(MARKED,IFMAR,3))
     $                 .LT. XFAC/1000.0  .AND. MARKED .NE.IEL  ) THEN
C                          We found an overlap
C                          Calculate what new value of ISPLIT would be
                           IF(IFMAR .EQ.5 .OR. IFMAR.EQ.6) THEN
C                            Different level
C                            Find rotation between elements
C                            IMAR=5:Check floor of marked with ceiling of new
C                            IMAR=6:Check ceiling of marked with floor of new
                             RMIN=1.0E10
                             DO 120 IC=1,4
                                RAD5=SQRT((X(MARKED,1)-X(IEL,IC+4))**2+
     $                                    (Y(MARKED,1)-Y(IEL,IC+4))**2)
                                RAD6=SQRT((X(MARKED,5)-X(IEL,IC  ))**2+
     $                                    (Y(MARKED,5)-Y(IEL,IC  ))**2)
                                IF(IFMAR.EQ.5)RAD=RAD5
                                IF(IFMAR.EQ.6)RAD=RAD6
                                IF(RAD .LT. RMIN) THEN
C                                  Overlapping point
                                   RMIN=RAD
                                   IROTAT=IC-1
                                ENDIF
120                          CONTINUE
                             NSPLIT=ISPLIT(MARKED) + IROTAT
                             IF(NSPLIT.GT.4)NSPLIT=NSPLIT-4
                           ELSE
C                            Touching on one of 4 sides
                             IF( ISFACE-IFMAR .EQ. 1
     $                       .OR.ISFACE-IFMAR .EQ.-3)THEN
C                               Isface IS CCW FROM overlap
                                NSPLIT=IFNEW-1
                                IF(IFNEW-1.EQ.0)NSPLIT=4
                             ELSE IF(ISFACE-IFMAR .EQ.-1
     $                       .OR.    ISFACE-IFMAR .EQ. 3)THEN
C                               Isface IS CCW FROM overlap
                                NSPLIT=IFNEW+1
                                IF(IFNEW+1.EQ.5)NSPLIT=1
                             ELSE
                                CALL PRS('ERROR: CANT ORIENT SPLIT$' )
                             ENDIF
                          ENDIF
C
                           IF(ISPLIT(IEL).EQ.0)THEN
C                            We'll split a new, unsplit element
                             NEWSPL=NEWSPL+1
                             ISPLIT(IEL)=NSPLIT
                           ELSE
C                            Make sure old split element doesn't conflict
                             IF(ISPLIT(IEL).NE.NSPLIT)THEN
                               CALL PRSIS('ERROR: ELEMENT$',IEL,
     $                         'WOULD NEED TO BE SPLIT TWICE.  '//
     $                         'ABORTING SPLIT.$')
                               call prsii('isplit,nsplit$'
     $                             ,isplit(iel),nsplit)

                               call prs('Continue splitting?$')
                               call res(ans,1)

                               if (ans.eq.'n'.or.ans.eq.'N') then
                                 do i=1,nelm
                                   ISPLIT(I)=0
                                 enddo
                                 goto 301
                               endif

                             ENDIF
                          ENDIF
                        ENDIF
80                   CONTINUE
                    ENDIF
90               CONTINUE
100            CONTINUE
            ENDIF
200      CONTINUE
C        We're done if no new elements got split
         IF(NEWSPL.EQ.0) GO TO 301
300   CONTINUE
      CALL PRS
     $('ERROR; SOME ELEMENTS MAY BE LEFT UNSPLIT ! PLEASE CHECK.$')
301   CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine mark2(ielcrk,ifacrk,sfrac,isplit)
C     Mark where CRACK propagates in vector ISPLIT
      include 'basics.inc'
      INTEGER ISPLIT(NELM)
C
      DO 10 I=1,NELM
         ISPLIT(I)=0
10    CONTINUE
      ISPLIT(IELCRK)=IFACRK
C
      return
      end
c-----------------------------------------------------------------------
      subroutine split(sfrac,isplit)
C     Split elements marked in vector ISPLIT
      include 'basics.inc'
      INTEGER ISPLIT(1000),IIPOINT(2,3,2)
      CHARACTER LETNEW(1000)
C     IPOINT: 1st subscript-which end of a; 2nd subscript- a or b; 3rd-top/bot
C     First we modify element.  Move oposite face toward IFACRK.  Insert
C     new element in gap.
c
C     Set up stuff for letter of daughter element
      MAXNW=MAXLET
      DO 20 I=1,MAXLET
         LETNEW(I)=' '
20    CONTINUE
      IA=1
      IB=2
      NELOLD=NEL
      DO 100 IEL=1,NELOLD
         IF(ISPLIT(IEL).NE.0)THEN
C           First make copy
            WRITE(S,'(1X,A18,I5,A6,I5,I5)')
     $      'SPLITTING ELEMENT ',IEL,'  SIDE',ISPLIT(I),NEL
            CALL PRS(S//'$')
            NEL=NEL+1
            CALL COPYEL(IEL,NEL)
C           Give Daughter element a letter
            IF(LETNEW(ICHAR(LETAPT(IEL))) .EQ. ' ')THEN
C              A parent with this letapt hasn't split yet; assign new pointer
               MAXNW=MAXNW+1
               LETNEW(ICHAR(LETAPT(IEL))) = CHAR(MIN0(122,MAXNW))
            ENDIF
            LETAPT(NEL) = LETNEW(ICHAR(LETAPT(IEL)))
C           That was, Each old letter points to a new letter for the daughter
C           Establish pointers to the element corners
            IIPOINT(1,IA,1)=ISPLIT(IEL)
            IIPOINT(2,IA,1)=ISPLIT(IEL)+1
            IF(IIPOINT(2,IA,1).EQ.5) IIPOINT(2,IA,1)=1
C
            IIPOINT   (2,IB,1)     = IIPOINT(2,IA,1)+1
            IF(IIPOINT(2,IB,1).EQ.5) IIPOINT(2,IB,1)=1
            IIPOINT   (1,IB,1)     = IIPOINT(1,IA,1)-1
            IF(IIPOINT(1,IB,1).EQ.0) IIPOINT(1,IB,1)=4
C
            DO 50 I12=1,2
               DO 50 IABC=1,3
                  IIPOINT(I12,IABC,2) = IIPOINT(I12,IABC,1) + 4
50          CONTINUE
C           Now move corners
            DO 60 I12=1,2
            DO 60 ITB=1,2
              ICA=IIPOINT(I12,IA,ITB)
              ICB=IIPOINT(I12,IB,ITB)
              IF(I12.EQ.1) IEDGE=ICB
              IF(I12.EQ.2) IEDGE=ICA
C             Move corners of new element
              IF(CCURVE(IEDGE,IEL).EQ.' ')THEN
C                Straight line
                 X(NEL,ICA)=X(NEL,ICA) + (X(NEL,ICB)-X(NEL,ICA))*SFRAC
                 Y(NEL,ICA)=Y(NEL,ICA) + (Y(NEL,ICB)-Y(NEL,ICA))*SFRAC
                 Z(NEL,ICA)=Z(NEL,ICA) + (Z(NEL,ICB)-Z(NEL,ICA))*SFRAC
              ELSE
                IF(I12.EQ.2)
     $          CALL GETPTS(1,    SFRAC,IEL,IEDGE,X(NEL,ICA),Y(NEL,ICA))
                IF(I12.EQ.1)
     $          CALL GETPTS(1,1.0-SFRAC,IEL,IEDGE,X(NEL,ICA),Y(NEL,ICA))
              ENDIF
              Z(NEL,ICA)=Z(NEL,ICA) + (Z(NEL,ICB)-Z(NEL,ICA))*SFRAC
C             Move corners of old element
              X(IEL,ICB)=X(NEL,ICA)
              Y(IEL,ICB)=Y(NEL,ICA)
              Z(IEL,ICB)=Z(NEL,ICA)
              IF(CCURVE(IEDGE,IEL).EQ.'S')THEN
C                Modify Control points
                 IF(I12.EQ.2)THEN
                    CURVE(3,IEDGE,IEL)=X(NEL,ICB)
                    CURVE(4,IEDGE,IEL)=Y(NEL,ICB)
                    CURVE(1,IEDGE,NEL)=X(IEL,ICA)
                    CURVE(2,IEDGE,NEL)=Y(IEL,ICA)
                 ELSE IF(I12.EQ.1)THEN
                    CURVE(3,IEDGE,NEL)=X(IEL,ICA)
                    CURVE(4,IEDGE,NEL)=Y(IEL,ICA)
                    CURVE(1,IEDGE,IEL)=X(NEL,ICB)
                    CURVE(2,IEDGE,IEL)=Y(NEL,ICB)
                 ENDIF
              ENDIF
60          CONTINUE
C           Make reasonable b.c.'s along split
C           Any splines present kill curve, TWO circles get average;
C           One circle?? currently yields straight.
            DO 65 ITB=1,2
C              Edges
               NELCUT=IIPOINT(1,IA,ITB)
               IELCUT=IIPOINT(2,IB,ITB)
               IF(CCURVE(NELCUT,NEL).EQ.'C'.AND.
     $            CCURVE(IELCUT,IEL).EQ.'C')THEN
C                 Get new radius--Tricky.
                  CURVE(1,IELCUT,IEL)=-CURVE(1,NELCUT,IEL)+
     $           (CURVE(1,IELCUT,IEL)+ CURVE(1,NELCUT,IEL))*SFRAC
                  CURVE(1,NELCUT,NEL)=-CURVE(1,IELCUT,IEL)
               ELSE IF(CCURVE(NELCUT,NEL).EQ.'S'.OR.
     $                 CCURVE(IELCUT,IEL).EQ.'S')THEN
                  CCURVE(IELCUT,IEL)=' '
                  CCURVE(NELCUT,NEL)=' '
               ELSE IF( (CCURVE(NELCUT,NEL).EQ.'C'.AND.
     $                   CCURVE(IELCUT,IEL).EQ.' ')
     $          .OR.    (CCURVE(NELCUT,NEL).EQ.' '.AND.
     $                   CCURVE(IELCUT,IEL).EQ.'C') )THEN
                  CURVE(1,IELCUT,IEL)=CURVE(1,NELCUT,IEL)+
     $           (CURVE(1,IELCUT,IEL)-CURVE(1,NELCUT,IEL))*SFRAC
                  CCURVE(IELCUT,IEL)=' '
                  CCURVE(NELCUT,NEL)=' '
               ELSE IF(  CCURVE(NELCUT,NEL).EQ.'O' .OR.
     $                   CCURVE(IELCUT,IEL).EQ.'O'  ) THEN
                  CCURVE(IELCUT,IEL)=' '
                  CCURVE(NELCUT,NEL)=' '
               ENDIF
65          CONTINUE
C           Delete periodic b.c.'s to avoid confusion
            DO 70 IS=1,NSIDES
               DO 70 IF=1,2
                  IF(CBC (IS,NEL,IF).EQ.'P  ')THEN
                       CBC (IS,NEL,IF)=' '
                       BC(1,IS,NEL,IF)=0.0
                       BC(2,IS,NEL,IF)=0.0
                       BC(3,IS,NEL,IF)=0.0
                  ENDIF
70          CONTINUE
C
         ENDIF
100   CONTINUE
C     Recalculate centers
      DO 110 IEL=1,NEL
        XCEN(IEL)=(X(IEL,1)+X(IEL,2)+X(IEL,3)+X(IEL,4))/4.
        YCEN(IEL)=(Y(IEL,1)+Y(IEL,2)+Y(IEL,3)+Y(IEL,4))/4.
110   CONTINUE
C
      CALL SORTEL
C     Draw whole new isometric surface
      DO 51 IEL=1,NEL
c        CALL DRAWIS(ISRT(IEL))
         IF(.NOT.IF3D  .OR. NUMAPT(IEL).EQ.ILEVEL) CALL DRAWEL(IEL)
51    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine vertadj
C
      include 'basics.inc'
      integer nelold
      save    nelold
      data    nelold /0/
C
      integer icalld
      save    icalld
      data    icalld/0/
c
      icalld=icalld+1
      io = icalld+75
C
      CALL GENCEN
      RMIN=GLMIN(RCEN,NEL)
C
      WRITE(S,10) RMIN
   10 FORMAT(2X,'RMIN:',g13.6,'$')
      CALL PRS(S)
C
c     CALL PRS('Input tolerance, relative to el. size (e.g., 0.5:)$')
c     CALL RER(EPS)
      eps = 0.05
C
      WRITE(S,20) EPS,RMIN
   20 FORMAT(2X,'EPS,RMIN:',2E13.6,'$')
      CALL PRS(S)
      EPS2=EPS**2
      epsr=eps*rmin
C
      ncrnr=2**ndim
      ichk=10
      if (nel.gt.200) ichk=50
      iadj=0
C
      if (nelold.ne.nel) then
         nelold=nel
         DO 500 IE=1,NEL
         DO 500 IC=1,NCRNR
            X(IE,IC)=ROUND(X(IE,IC))
            Y(IE,IC)=ROUND(Y(IE,IC))
            Z(IE,IC)=ROUND(Z(IE,IC))
  500    CONTINUE
      endif
C
      epsr=.0001*rmin
      DO 1000 IE=1,NEL
C
         if (mod(ie,ichk).eq.0.and.iadj.eq.0) then
            WRITE(6,1001) IE
         else
            iadj=0
         endif
C
         DO 100 IC=1,NCRNR
C
            if (abs(x(ie,ic)).lt.epsr.and.x(ie,ic).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing x:',ie,ic,x(ie,ic)
               CALL PRS(S//'$')
               x(ie,ic)=0.0
               iadj=1
            endif
C
            if (abs(y(ie,ic)).lt.epsr.and.y(ie,ic).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing y:',ie,ic,y(ie,ic)
               CALL PRS(S//'$')
               y(ie,ic)=0.0
               iadj=1
            endif
C
            if (abs(z(ie,ic)).lt.epsr.and.z(ie,ic).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing z:',ie,ic,z(ie,ic)
               CALL PRS(S//'$')
               z(ie,ic)=0.0
               iadj=1
            endif
  100    CONTINUE
C
         DO 800 JE=IE+1,NEL
            DIST1=(XCEN(IE)-XCEN(Je))**2+(YCEN(IE)-YCEN(Je))**2
     $           +(ZCEN(IE)-ZCEN(Je))**2
c           DIST2=       RCEN(IE)**2 + RCEN(JE)**2
            DIST2= 1.01*( RCEN(IE) + RCEN(JE) )**2
c           if (ie.eq.514.and.je.eq.3377) then
c              write(6,*) ie,xcen(ie),ycen(ie),zcen(ie),rcen(ie)
c              write(6,*) je,xcen(je),ycen(je),zcen(je),rcen(je)
c              write(6,*) dist1,dist2
c              write(6,*) 'go?'
c              call res(s,1)
c           endif
            IF (DIST1.LE.DIST2) THEN
c              epsrr = .01*min(rcen(ie),rcen(je))
               epsrr = eps2*min(rcen(ie),rcen(je))
               epsrr = epsrr**2
               DO 700 IC=1,NCRNR
               DO 700 JC=1,NCRNR
                  DIST1=(X(IE,IC)-X(JE,JC))**2+(Y(IE,IC)-Y(JE,JC))**2
     $                 +(Z(IE,IC)-Z(JE,JC))**2
                  IF (DIST1.LT.epsrr.and.dist1.ne.0.0) THEN
                     if (nel.lt.2000) then
                        write(s,'(1x,a10,2(i6,i3))')
     $                  'Adjusting:',je,jc,ie,ic
                        CALL PRS(S//'$')
c
c                   diag.
c                      if (ie.eq.514.and.je.eq.3377) then
c                         write(6,*) ic,x(ie,ic),y(ie,ic),z(ie,ic),dist1
c                         write(6,*) jc,x(ie,jc),y(ie,jc),z(ie,jc),eps2
c                         write(6,*) 'go?'
c                         call res(s,1)
c                      endif
c                      write(io,75) ie,ic,x(ie,ic),y(ie,ic),z(ie,ic),je,jc
c                      write(io,76) je,jc,x(je,jc),y(je,jc),z(je,jc),ie,ic
c                   diag.
c
                        call prrrr(x(ie,ic),y(ie,ic),z(ie,ic))
                        call prrrr(x(je,jc),y(je,jc),z(je,jc))
                     endif
c
                     X(JE,JC)=X(IE,IC)
                     Y(JE,JC)=Y(IE,IC)
                     Z(JE,JC)=Z(IE,IC)
C
c                diag.
c                  write(io,77) je,jc,x(je,jc),y(je,jc),z(je,jc),ie,ic
c                  call prrrr(x(je,jc),y(je,jc),z(je,jc))
C
   75 format('I',i6,i2,3g18.9,I8,i2)
   76 format('J',i6,i2,3g18.9,I8,i2)
   77 format('j',i6,i2,3g18.9,I8,i2)
                     iadj=1
                  ENDIF
  700          CONTINUE
            ENDIF
  800    CONTINUE
 1000 CONTINUE
 1001 FORMAT('  Checking',I5)
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine clpdom
      include 'basics.inc'
      integer  del_list(nelm)
      real     xyznorm(3),xyzclip(3)
C
    1 CONTINUE
      nchoic = 1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Clip X'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Clip Y'
      IF (IF3D) THEN
         nchoic = nchoic+1
         ITEM(nchoic)    =             'Clip Z'
      ENDIF
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Clip N'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Input DEL list'
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
C
      IF (CHOICE.EQ.'UP MENU') RETURN
      IF (CHOICE.EQ.'Clip X') THEN
         CALL PRS(
     $   'Input location of X-clipping plane.$')
         CALL RER(Xclip)
c        CALL DRAWLINE(Xclip,Ymax,Xclip,Ymin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') RETURN
         CALL Clipper(Xclip,ANS,X)
      ELSEIF (CHOICE.EQ.'Clip Y') THEN
         CALL PRS(
     $   'Input location of Y-clipping plane.$')
         CALL RER(Yclip)
c        CALL DRAWLINE(Yclip,Xmax,Yclip,Xmin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') RETURN
         CALL Clipper(Yclip,ANS,Y)
      ELSEIF (CHOICE.EQ.'Clip Z') THEN
         CALL PRS(
     $   'Input location of Z-clipping plane.$')
         CALL RER(Zclip)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') RETURN
         CALL Clipper(Zclip,ANS,Z)
      ELSEIF (CHOICE.EQ.'Clip N') THEN
         CALL PRS('Allows to clip below/above arbitrary plane.$')
         if (if3d) then
            CALL PRS('Input xyz location of clipping plane.$')
            call rerrr(xyzclip(1),xyzclip(2),xyzclip(3))
            CALL PRS('Input xyz normal of clipping plane.$')
            call rerrr(xyznorm(1),xyznorm(2),xyznorm(3))
         else
            CALL PRS('Input xy location of clipping plane.$')
            call rerr(xyzclip(1),xyzclip(2))
            CALL PRS('Input xy normal of clipping plane.$')
            call rerr(xyznorm(1),xyznorm(2))
            xyzclip(3)= 0.
            xyznorm(3)= 0.
         endif
c
c     normalize
         xyzl = xyznorm(1)*xyznorm(1)
     $        + xyznorm(2)*xyznorm(2)
     $        + xyznorm(3)*xyznorm(3)
         if (xyzl.le.0) return
         xyzl = sqrt(xyzl)
         xyznorm(1) = xyznorm(1)/xyzl
         xyznorm(2) = xyznorm(2)/xyzl
         xyznorm(3) = xyznorm(3)/xyzl
c
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') RETURN
         CALL ClipperN(xyzclip,xyznorm,ANS,x,y,z,del_list)
      ELSEIF (CHOICE.EQ.'Input DEL list') THEN
         call mult_del(del_list)
      ENDIF
      GOTO 1
      END
c-----------------------------------------------------------------------
      subroutine clipper(Clip,DIR,pts)
      include 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,YESNO
C
      Nvts = 4
      IF (IF3D) Nvts=8
C
      IF (DIR.eq.'>') THEN
C
C        Count number
         NUMDEL=0
         DO 101 JE=1,NEL
            DO 100  i=1,Nvts
               if (je.le.nel .and. pts(je,i).gt.clip) then
                   NUMDEL=NUMDEL+1
                   GOTO 101
               endif
  100       CONTINUE
  101    CONTINUE
         write(s,102) numdel,nel
  102    format(' You will be eliminating',i5,' of ',i5,
     $   ' elements ABOVE clipping plane.$')
         call prs(s)
         call prs(' OK? (Y/N)$')
         call res(yesno,1)
C
         if (yesno.eq.'y'.or.yesno.eq.'Y') THEN
C
C           If any part of an element is above the clipping plane, delte it.
C           This is a gross N^2 algorithm, but it beats entering them all
C           by hand....
C
            NUMDEL=0
  110       CONTINUE
            NELT=NEL
            DO 120 JE=1,NELT
            DO 120  i=1,Nvts
               if (je.le.nel .and. pts(je,i).gt.clip) then
                   call delelq(je)
                   NUMDEL=NUMDEL+1
                   GOTO 110
               ENDIF
  120       CONTINUE
            IF (NUMDEL.GT.0.and.nel.le.100) THEN
C              redraw the mesh
               CALL REFRESH
               CALL DRMENU('NOCOVER')
               CALL DRGRID
               write(s,500) numdel,nel
               call prs(s)
               DO 130 IEL=1,NEL
                  CALL DRAWEL(IEL)
  130          CONTINUE
            ENDIF
         ENDIF
C
      ELSEIF (DIR.eq.'<') THEN
C
C        If any part of an element is below the clipping plane, delte it.
C        Count number
         NUMDEL=0
         DO 201 JE=1,NEL
            DO 200  i=1,Nvts
               if (je.le.nel.and.pts(je,i).lt.clip) then
                   NUMDEL=NUMDEL+1
                   GOTO 201
               endif
  200       CONTINUE
  201    CONTINUE
         write(s,202) numdel,nel
  202    format(' You will be eliminating',i5,' of ',i5,
     $   ' elements BELOW clipping plane.$')
         call prs(s)
         call prs(' OK? (Y/N)$')
         call res(yesno,1)
C
         if (yesno.eq.'y'.or.yesno.eq.'Y') THEN
            NUMDEL=0
  210       CONTINUE
            NELT=NEL
            DO 220 JE=1,NELT
            DO 220  i=1,Nvts
               if (je.le.nel .and. pts(je,i).lt.clip) then
                   call delelq(je)
                   NUMDEL=NUMDEL+1
                   GOTO 210
               ENDIF
  220       CONTINUE
            IF (NUMDEL.GT.0.and.nel.le.100) THEN
C              redraw the mesh
               CALL REFRESH
               CALL DRMENU('NOCOVER')
               CALL DRGRID
               write(s,500) numdel,nel
               call prs(s)
               DO 230 IEL=1,NEL
                  CALL DRAWEL(IEL)
  230          CONTINUE
            ENDIF
         ENDIF
      ENDIF
  500 format(i7,' elements deleted.',i7,' elements remaining.$')
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine stretch
      call shift2
      return
      end
c-----------------------------------------------------------------------
      subroutine smoother(cell,nvc,ia,ja,ww,b,g,x0,x1,y0,y1,z0,z1)
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      real ww(0:1)
      integer ia(1),ja(1)
c
      integer cell(nvc,1),b(1),g(1)
c
c
      integer edge(2,12)
      save    edge
      data    edge  / 1,2 , 2,3 , 3,4 , 4,1
     $              , 5,6 , 6,7 , 7,8 , 8,5
     $              , 1,5 , 2,6 , 3,7 , 4,8 /
c
      integer iebad(nelm)
c
c
      icount= 0
      nebad = 0
      call izero(iebad,nelm)
 1000 continue
      icount = icount+1
c
c     use connectivity graph to smooth mesh
c
      n1 = nvc*nel
      call get_gxyz(nv,cell,nvc,ncell,ww(n1),ww,ierr)
c
c     Identify boundary and curve side vertices
      nvc=2**ndim
      nfaces=2*ndim
      call izero(b,nv)
      call find_bc_crv(b,nv,cell,nvc,ncell)
c
c     Identify dof's inside smoothing box
      call find_sm_box(b,nv,x0,x1,y0,y1,z0,z1)
c
c     Identify dof's connected to elements yielding bad Jacobians
      call find_bd_jac(b,nv,cell,nvc,nel,iebad,nebad)
c
c
c
c     Build adjaceny graph
c
c
c     Build adjaceny graph
      ne = 2*ndim*(ndim-1)
      n1 = 14*ne*ncell
      n2 = 14*ne*ncell + n1
      call build_adj  (ja,ia,n,cell,nvc,ncell,edge,ne,ww,ww(n1),ww(n2))
c     call out_matlab (ja,ia,n,6)
c
c     Construct a connectivity list telling the topological distance from
c     boundary points
c
      call boundary_graph(g,b,ja,ia,n)
c
c
c     Smooth the mesh
c
      eps = 0.10
      call prs
     $    ('input relaxation params (s,eps,gg) (typ 0.9,.01,1-999):$')
      call rerrr(ss,eps,gg)
c
      nnz = ia(n+1)-ia(1)
      do iter=1,80
         dxi=gs_smooth(xp,yp,zp,b,ww,ww(nnz),ja,ia,n,g,gg,ss,eps,iter)
         write(6,*) iter,' smooth:',dxi
      enddo
      call chk_jacob(nebad,iebad,b,xp,yp,zp,cell,nvc)
      if (nebad.gt.0.and.icount.lt.3) then
         write(6,*) 'Bad jacobian',nebad,icount
c        do je=1,nebad
c        write(6,*) 
c        do j =1,nvc
c           ie = iebad(je)
c           write(6,3) ie,xp(cell(j,ie)),yp(cell(j,ie)),zp(cell(j,ie))
c        enddo
c        enddo
    3    format(i9,1p3e14.5)
         write(6,*) 'Bad jacobian',nebad,icount
         goto 1000
      elseif (nebad.gt.0) then
         write(6,*) 'Still Bad jacobian, no smoothing',nebad,icount
      endif
c
      do ie=1,nel
      do i=1,nvc
        xo=x(ie,i)
        if (b(cell(i,ie)).le.0) then
           x(ie,i) = xp(cell(i,ie))
           y(ie,i) = yp(cell(i,ie))
           z(ie,i) = zp(cell(i,ie))
        endif
        jc=cell(i,ie)
        write(66,*) b(jc),' xo:',xo,x(ie,i),i,ie,jc,g(jc)
      enddo
      enddo
    1 format('xc:',2i4,i9,4f12.5)
    2 format(a3,2i4,i9,2f12.5)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine find_bc_crv(b,nv,cell,nvc,ncell)
c
      include 'basics.inc'
      integer cell(nvc,ncell),b(nv)
c
      integer efc(4,6)
      save    efc
      data    efc  / 1,2,5,6
     $             , 2,3,6,7
     $             , 3,4,7,8
     $             , 4,1,8,5
     $             , 1,2,3,4
     $             , 5,6,7,8 /
c
c
c     Identify any elements on the boundary or having curved sides
c
      call izero(b,nv)
c
      nfaces = 2*ndim
      ncrnf  = 2**(ndim-1)
      ifld = 2
      if (ifflow) ifld = 1
c
      do ie=1,nel
      do is=1,nfaces
c
         if (cbc(is,ie,ifld).ne.'E  ') then
            do iv=1,ncrnf
               b(cell(efc(iv,is),ie)) =  1
            enddo
         endif
c
         if (cbc(is,ie,ifld).eq.'SYM') then
c           Later this will be fixed to handle each spatial dimension
            do iv=1,ncrnf
               b(cell(efc(iv,is),ie)) = -2
            enddo
         endif
c
c        Check for curve sides (spherical in particular)
         if (ccurve(is,ie).eq.'s') then
            do iv=1,ncrnf
               b(cell(efc(iv,is),ie)) =  1
            enddo
         endif
       enddo
       enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_adj(ja,ia,n,cell,nvc,ncell,edge,ne,w,wk,ind)
c
c     Build an adjacency graph from cell-based data, given
c     edge array
c
c
      integer ja(1),ia(1),cell(nvc,ncell),w(2,1),ind(1),wk(1)
      integer edge(2,ne)
      integer key (2)
      integer w1,w2
c
      k=0
      do icell=1,ncell
      do ie   =1,ne
         w1 = cell(edge(1,ie),icell)
         w2 = cell(edge(2,ie),icell)
         k=k+1
         w(1,k) = min(w1,w2)
         w(2,k) = max(w1,w2)
      enddo
      enddo
6     format(a3,8i4)
      nv = k
c
c
      key(1) = 1
      key(2) = 2
      call ituple_s_merge(nw,w,2,nv,key,2,ind,wk)
      write(6,*) 'this is nw:',nw,nv
c
c
c     Get symmetric other-half
c
      w1 = 0
      w2 = 0
      do i=1,nw
         w(1,i+nw) = w(2,i)
         w(2,i+nw) = w(1,i)
         w1 = max(w1,w(1,i))
         w2 = max(w2,w(2,i))
      enddo
      nv = 2*nw
      write(6,*) 'this is nw:',nw,nv,w1,w2
      call ituple_s_merge(nw,w,2,nv,key,2,ind,wk)
      write(6,*) 'this is nw3',nw,nv,w1,w2
c
c
c     Take merged list of vertex pairs and build csr format adjancy array
c
      n     = 1
      nnz   = 1
      ia(1) = 1
      ja(1) = w(1,1)
      ilast = w(1,1)
      nnz   = nnz+1
      ja(2) = w(2,1)
      do i=2,nv
         if (w(1,i).eq.ilast) then
            nnz = nnz+1
            ja(nnz) = w(2,i)
         else
            nnz   = nnz+1
            n     = n+1
            ia(n) = nnz
            ja(nnz) = w(1,i)
            ilast   = w(1,i)
            nnz = nnz+1
            ja(nnz) = w(2,i)
         endif
      enddo
      ia(n+1) = nnz+1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ituple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_ialtb
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
c           aa  = a  (l)
            call icopy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call icopy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call icopy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call icopy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
c              if ( a(j).lt.a(j+1) ) j=j+1
               if (iftuple_ialtb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
c           if (aa.lt.a(j)) then
            if (iftuple_ialtb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call icopy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call icopy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine tuple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_altb
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
c           aa  = a  (l)
            call copy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call copy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call copy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call copy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
c              if ( a(j).lt.a(j+1) ) j=j+1
               if (iftuple_altb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
c           if (aa.lt.a(j)) then
            if (iftuple_altb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call copy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call copy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      logical function iftuple_ialtb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_ialtb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_ialtb = .false.
            return
         endif
      enddo
      iftuple_ialtb = .false.
      return
      end
c-----------------------------------------------------------------------
      logical function iftuple_altb(a,b,key,nkey)
      real a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_altb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_altb = .false.
            return
         endif
      enddo
      iftuple_altb = .false.
      return
      end
c-----------------------------------------------------------------------
      subroutine ituple_s_merge(na,a,lda,n,key,nkey,ind,aa)
c
      integer a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_ialtb
c
c     Sort a(tuple,i) and merge by removing repeated entries
c
      call ituple_sort(a,lda,n,key,nkey,ind,aa)
c
c     Now merge
c
      na = 1
      do j=2,n
         if (iftuple_ialtb(a(1,j-1),a(1,j),key,nkey)) then
c           different, bump pointer
            na = na+1
            do l=1,lda
               a(l,na) = a(l,j)
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      function gs_smooth(x,y,z,b,c,a,ja,ia,n,g,gg,s,eps,icall)
      real     c(1),a(1),x(1),y(1),z(1)
      integer  ja(1),ia(1),b(1),g(1)
c
c
c     Input parameters
c
c
c     b(i)  -- boundary (and box) information, b(i)>0 ==> point doesn't move
c     icall -- 1 ==> compute initial distance,  a(ij)
c     s     -- stretch.  s=.85 implies that a(ij)/d_final = 1./0.85.
c     eps   -- scaling factor for amount of movement on this pass
c     g(i)  -- graph distance from node i to boundary
c     gg    -- gain adjustment factor, such that s_ij = s*gg/(gi+gj+gg-1)
c              gg = 1--10 typical.  gg=1e9 ==> ~no boundary gain
c
c
      if (icall.eq.1) then
c
c        Build stiffness matrix based on original xj-xi distances
c
         do i=1,n
            j0 = ia(i)
            j1 = ia(i+1)-1
            do ij=j0,j1
               j = ja(ij)
               if (j.ne.i) then
                  dx  = (x(i)-x(j))
                  dy  = (y(i)-y(j))
                  dz  = (z(i)-z(j))
                  dx2 = dx*dx + dy*dy + dz*dz
                  ds  = sqrt(dx2)
                  a(ij) = ds
                  c(ij) = 1.
c                 Increase the stiffness for links connected to the boundary?
c                 if (b(i).eq.1.or.b(j).eq.1) c(ij) = 2.
                  if (gg.gt.1 .or. g(i).ne.0 .or. g(j).ne.0) then
                     c(ij) = gg/(g(i)+g(j)+gg-1.)
                  else
                     c(ij) = 0.
                  endif
               else
                  a(ij) = 0.
               endif
c              write(6,1) i,j,g(i),g(j),gg,c(ij),a(ij)
c  1           format(4i6,1p3e12.5)
            enddo
         enddo
      endif
c
      dmx = 0.
      do i=1,n
         if (b(i).le.0) then
            j0 = ia(i)
            j1 = ia(i+1)-1
            dxa=0.
            dya=0.
            dza=0.
            dsm=1.e22
            do ij=j0,j1
               j = ja(ij)
               if (j.ne.i) then
                  dx  = -(x(i)-x(j))
                  dy  = -(y(i)-y(j))
                  dz  = -(z(i)-z(j))
                  dx2 = dx*dx + dy*dy + dz*dz
                  ds  = sqrt(dx2)
                  dsm = min(ds,dsm)
c
c                 Strength of link is proportional to 
c
c                 c(ij) := 1./(g(i)+g(j)),  g(i) := dist | v(i) - boundary |
c
                  scale = c(ij)*(1.-s*a(ij)/ds)
                  dxa = dxa + scale*dx
                  dya = dya + scale*dy
                  dza = dza + scale*dz
               endif
            enddo
            xo = x(i)
            yo = y(i)
            zo = z(i)
c
c           This will ultimately allow planar constraints
c
            if (b(i).eq.-1) dxa = 0
            if (b(i).eq.-2) dya = 0
            if (b(i).eq.-3) dza = 0
            if (b(i).gt.0) then
               dxa = 0
               dya = 0
               dza = 0
            endif
c
            x(i) = x(i) + eps*dxa
            y(i) = y(i) + eps*dya
            z(i) = z(i) + eps*dza
            dxn = dxa*dxa + dya*dya + dza*dza
            dxn = eps*sqrt(dxn)
            dmx = max(dmx,dxn)
         endif
      enddo
      gs_smooth = dmx
      return
      end
c-----------------------------------------------------------------------
      subroutine get_gxyz(nv,cell,nvc,ncell,ww,ind,ierr)
c
c     Read a file until "key" is found or eof is found.
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      common /cfilold/ filold
      character filold*17
c
      integer cell(nvc,1),ind(1),ww(0:1)
c
      common /file_pref/ file_prefix
      character*80 file_prefix
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
c
      integer ecrnr(8),kcell(8)
      save    ecrnr
      data    ecrnr  / 1,2,4,3,5,6,8,7 /
c
   80 format(a80)
   81 format(80a1)
c
c
c     Get map file from old input
c
      length = indx1(filold,'.',1)-1
      if (length.le.0) length = indx1(filold,' ',1)-1
      call blank (string,80)
      call chcopy(string,filold,length)
      call chcopy(string1(length+1),'.map',4)
c
c     Get current vertex map info
c
      open(unit=10,file=string,err=999)
      read(10,*) ncell
      do ie=1,ncell
c        read(10,*,end=998,err=998) idum,(cell(k,ie),k=1,nvc)
c
c        HMT's data (.map) is in the good h-cube ordering
         read(10,*,end=998,err=998) idum,(kcell(k),k=1,nvc)
         do k=1,nvc
            j=ecrnr(k)
            cell(k,ie) = kcell(j)
         enddo
      enddo
      close(unit=10)
c
c
c     Sort data and gridpoints by global vertex number
c
      lmax  = 0
      do ie = 1,ncell
         do k=1,nvc
c           j     = ecrnr(k)
            l     = cell(k,ie)
            xp(l) = x (ie,k)
            yp(l) = y (ie,k)
            zp(l) = z (ie,k)
            lmax  = max(l,lmax)
c           write(6,6) 'x:',ie,k,l,xp(l),yp(l)
         enddo
      enddo
6     format(a2,3i9,2f12.4)
      nv = lmax
c
      ierr = 0
      return
c
  998 continue
      ierr = 1
      call prs('end of .map file reached')
      return
c
  999 continue
      ierr = 1
      call prs('Could not open .map file.')
      return
      end
c-----------------------------------------------------------------------
      subroutine out_matlab (ja,ia,n,io)
      integer ia(1),ja(1),n
c
      nnz = ia(n+1) - ia(1)
      write(6,*) 'out matl:',ia(1),n,nnz
      do i=1,n
         j0 = ia(i)
         j1 = ia(i+1)
         do j=j0,j1-1
            write(io,10) i,ja(j),j
         enddo
      enddo
   10 format(3i12)
      return
      end
c-----------------------------------------------------------------------
      subroutine outc(cell,nvc,ncell,k)
c
      integer cell(nvc,ncell)
      write(6,9) k,nvc,ncell
      do i=1,ncell
         write(6,8) k,i,(cell(j,i),j=1,nvc)
8        format(10i8)
      enddo
9     format('CELL OUT:',3i8)
      return
      end
c-----------------------------------------------------------------------
      subroutine smooth
      include 'basics.inc'
      include 'basicsp.inc'
c
      x0 = -9.e19
      y0 = -9.e19
      z0 = -9.e19
      x1 =  9.e19
      y1 =  9.e19
      z1 =  9.e19
C
    1 CONTINUE
c
      nchoic = 0
      nchoic = nchoic+1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Redraw mesh'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Set smoothing box'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Smooth'
      nchoic = nchoic+1
c
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
c
C
      IF (CHOICE.EQ.'UP MENU') RETURN
      IF (CHOICE.EQ.'Redraw mesh') then
         call redraw_mesh
      ELSEIF (CHOICE.EQ.'Set smoothing box') THEN
         CALL PRS(
     $   'Input X-locations defining smoothing box.$')
         CALL RERR(x0,x1)
         CALL PRS(
     $   'Input Y-locations defining smoothing box.$')
         CALL RERR(y0,y1)
         if (if3d) CALL PRS(
     $   'Input Z-locations defining smoothing box.$')
         if (if3d) CALL RERR(z0,z1)
c
      ELSEIF (CHOICE.EQ.'Smooth') THEN
       nvc = 2**ndim
c      call smoother(cell,nvc,ia ,ja  ,ww        ,b    )
c     subroutine smoother(cell,nvc,ia,ja,ww,b,g,x0,x1,y0,y1,z0,z1)
       call smoother
     $      (rxm1,nvc,wv1,wkv1,sdump(1,4),sdump(1,5),jacm1
     $                                         ,x0,x1,y0,y1,z0,z1)
      ENDIF
c
      GOTO 1
      END
c-----------------------------------------------------------------------
      subroutine find_sm_box(b,nv,x0,x1,y0,y1,z0,z1)
c
      include 'basics.inc'
      include 'basicsp.inc'
c
      integer b(1)
c
      call qb(b,nv,xp,x0,x1)
      call qb(b,nv,yp,y0,y1)
      if (if3d) call qb(b,nv,zp,z0,z1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine qb(b,nv,x,x0,x1)
c
      integer b(1)
      real    x(1)
c
      xn=min(x0,x1)
      xx=max(x0,x1)
      do i=1,nv
         if (x(i).lt.xn) b(i) = 10
         if (x(i).gt.xx) b(i) = 10
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine find_bd_jac(b,nv,cell,nvc,nel,iebad,nebad)
c
      integer b(nv),cell(nvc,nel),iebad(nebad)
c
      do je=1,nebad
         ie=iebad(je)
         do j=1,nvc
            i=cell(j,ie)
            b(i) = 1
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_jacob(nebad,iebad,b,xp,yp,zp,cell,nvc)
      include 'basics.inc'
c
      integer iebad(1),cell(nvc,1),b(1)
      real xp(1),yp(1),zp(1)
c
      parameter (lx3=2)
      common  /cjac/ dxm3(lx3,lx3),dxm3t(lx3,lx3),z3(lx3),w3(lx3)
     $  , xm3 (lx3,lx3,lx3), ym3 (lx3,lx3,lx3), zm3 (lx3,lx3,lx3)
     $  , xrm3(lx3,lx3,lx3), yrm3(lx3,lx3,lx3), zrm3(lx3,lx3,lx3)
     $  , xsm3(lx3,lx3,lx3), ysm3(lx3,lx3,lx3), zsm3(lx3,lx3,lx3)
     $  , xtm3(lx3,lx3,lx3), ytm3(lx3,lx3,lx3), ztm3(lx3,lx3,lx3)
     $  , jacm3(lx3,lx3,lx3)
      real jacm3
c
      INTEGER INV (8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
C
      nx3 = 2
      ny3 = 2
      nz3 = 2
      if (ndim.eq.2) nz3=1
c
      call legend(z3,w3,nx3)
      call dgll  (dxm3,dxm3t,z3,nx3,nx3)
c
      NXY3  = NX3*NY3
      NYZ3  = NY3*NZ3
      NXYZ3 = NX3*NY3*NZ3
      NTOT3 = NXYZ3
      nebad = 0
C
C
C     Compute isoparametric partials.
C
      IF (NDIM.EQ.2) THEN
C
C        Two-dimensional case
C
         do ie=1,nel
            do j=1,4
               i=inv(j)
               xm3(i,1,1) = x(ie,j)
               ym3(i,1,1) = y(ie,j)
            enddo
c           Overload with new mesh data
            do j=1,4
               if (b(cell(j,ie)).le.0) then
                  i=inv(j)
                  xm3(i,1,1) = xp(cell(j,ie))
                  ym3(i,1,1) = yp(cell(j,ie))
               endif
            enddo
c
c
C           Use the appropriate derivative- and interpolation operator in
C           the y-direction (= radial direction if axisymmetric).
            CALL MXM(DXM3,NX3,XM3,NX3,XRM3,NY3)
            CALL MXM(DXM3,NX3,YM3,NX3,YRM3,NY3)
            CALL MXM(XM3,NX3,DXM3t,NY3,XSM3,NY3)
            CALL MXM(YM3,NX3,DXM3t,NY3,YSM3,NY3)
            CALL RZERO   (JACM3,NTOT3)
            CALL ADDCOL3 (JACM3,XRM3,YSM3,NTOT3)
            CALL SUBCOL3 (JACM3,XSM3,YRM3,NTOT3)
            CALL CHKJAC  (ierr,JACM3,NXYZ3)
            if (ierr.gt.0) then
               nebad=nebad+1
               iebad(nebad) = ie
            endif
         enddo
      ELSE
C
C        Three-dimensional case
C
         do ie=1,nel
c
            do j=1,8
               i=inv(j)
               xm3(i,1,1) = x(ie,j)
               ym3(i,1,1) = y(ie,j)
               zm3(i,1,1) = z(ie,j)
            enddo
c           Overload with new mesh data
            do j=1,8
               i=inv(j)
               if (b(cell(j,ie)).le.0) then
                  xm3(i,1,1) = xp(cell(j,ie))
                  ym3(i,1,1) = yp(cell(j,ie))
                  zm3(i,1,1) = zp(cell(j,ie))
               endif
            enddo
c
c
            CALL MXM(DXM3,NX3,XM3,NX3,XRM3,NYZ3)
            CALL MXM(DXM3,NX3,YM3,NX3,YRM3,NYZ3)
            CALL MXM(DXM3,NX3,ZM3,NX3,ZRM3,NYZ3)
            do iz=1,nz3
               CALL MXM(XM3(1,1,IZ),NX3,DXM3t,NY3,XSM3(1,1,IZ),NY3)
               CALL MXM(YM3(1,1,IZ),NX3,DXM3t,NY3,YSM3(1,1,IZ),NY3)
               CALL MXM(ZM3(1,1,IZ),NX3,DXM3t,NY3,ZSM3(1,1,IZ),NY3)
            enddo
            CALL MXM(XM3,NXY3,DZTM3,NZ3,XTM3,NZ3)
            CALL MXM(YM3,NXY3,DZTM3,NZ3,YTM3,NZ3)
            CALL MXM(ZM3,NXY3,DZTM3,NZ3,ZTM3,NZ3)
c
            CALL RZERO   (JACM3,NTOT3)
            CALL ADDCOL4 (JACM3,XRM3,YSM3,ZTM3,NTOT3)
            CALL ADDCOL4 (JACM3,XTM3,YRM3,ZSM3,NTOT3)
            CALL ADDCOL4 (JACM3,XSM3,YTM3,ZRM3,NTOT3)
            CALL SUBCOL4 (JACM3,XRM3,YTM3,ZSM3,NTOT3)
            CALL SUBCOL4 (JACM3,XSM3,YRM3,ZTM3,NTOT3)
            CALL SUBCOL4 (JACM3,XTM3,YSM3,ZRM3,NTOT3)
            CALL CHKJAC  (ierr,JACM3,NXYZ3)
            if (ierr.gt.0) then
               nebad=nebad+1
               iebad(nebad) = ie
               do i=1,8
               write(6,4) ie,xm3(i,1,1),xp(cell(inv(i),ie)),jacm3(i,1,1)
               write(6,4) ie,ym3(i,1,1),yp(cell(inv(i),ie)),jacm3(i,1,1)
               write(6,4) ie,zm3(i,1,1),zp(cell(inv(i),ie)),jacm3(i,1,1)
                  write(6,4) 
               enddo
4              format(i5,1p4e14.6)
               write(6,*) ie,nebad
               call prs('continue?$')
               call res(string,1)
            endif
         enddo
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine chkjac(IERR,jac,n)
C
C     Check the array JAC for a change in sign.
C
      REAL JAC(N)
      ierr=1
      SIGN = JAC(1)
      DO 100 I=2,N
         IF (SIGN*JAC(I).LE.0.0) RETURN
  100 CONTINUE
      ierr=0
      return
      end
c-----------------------------------------------------------------------
      subroutine boundary_graph(g,b,ja,ia,n)
      integer  g(1),ja(1),ia(1),b(1)
      logical  ifdone
c
c
c     Build stiffness matrix based on original xj-xi distances
c
      call ifill(g,n,n)
      do i=1,n
         if (b(i).eq.1) g(i)=0
      enddo
      do ipass = 1,20
         ifdone=.true.
         do i=1,n
            j0 = ia(i)
            j1 = ia(i+1)-1
            do ij=j0,j1
               j = ja(ij)
c              a_ij is a link
               if (abs(g(i)-g(j)).gt.1) ifdone = .false.
               if (g(i).gt.g(j)) g(i)=g(j)+1
               if (g(j).gt.g(i)) g(j)=g(i)+1
c              write(6,1) i,j,ij,g(i),g(j)
c   1          format(5i7,'graph')
            enddo
         enddo
         if (ifdone) goto 10
      enddo
   10 continue
      maxgraph = iglmax(g,n)
      call prsii('npasses, graph depth:$',ipass,maxgraph)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine mult_del(del_list)
      include 'basics.inc'
      character*40 list_file
      integer del_list(nelm)
      character*1 yesno
c
      call prs('Input file containing list to be deleted:$')
      call res(list_file,40)
      open (unit=77,file=list_file,err=999)
C
C
C     Count number
      numdel=0
      do je=1,nel
         read(77,*,end=101) del_list(je)
         numdel = numdel+1
      enddo
  101 continue
c
      write(s,102) numdel,nel
  102 format('You will be eliminating',i6,' of ',i6,
     $' elements. OK? (Y/N)$')
      call prs(s)
      call res(yesno,1)
      if (yesno.ne.'y'.and.yesno.ne.'Y') return
      call prs('Deleting elements.$')
c
c
c     Sort into descending element order for "delel" routine
c
      call iusort('d',numdel,del_list,info)
      if (info.ne.0) then
         call prsii('Problem w/ sort:$',info,numdel)
         call prs('No elements deleted.$')
         return
      endif
c
      do i=1,numdel
         je = del_list(i)
         call delel(je)
      enddo
c
      IF (NUMDEL.GT.0) THEN
C        redraw the mesh
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         write(s,500) numdel,nel
         call prs(s)
         DO 130 IEL=1,NEL
            CALL DRAWEL(IEL)
  130    CONTINUE
      ENDIF
  500 format(i7,' elements deleted.',i7,' elements remaining.$')
      RETURN
  999 continue
      call prs('Could not open file.  Returning.$')
      return
      END
c-----------------------------------------------------------------------
      subroutine ipsort( ID, N, D, INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      integer   D(1)
*     ..
*
*  Purpose
*  =======
*
*  Sort the numbers in D in increasing order (if ID = 'I') or
*  in decreasing order (if ID = 'D' ).
*
*  Use Quick Sort, reverting to Insertion sort on arrays of
*  size <= 20. Dimension of STACK limits N to about 2**32.
*
*  Arguments
*  =========
*
*  ID      (input) CHARACTER*1
*          = 'I': sort D in increasing order;
*          = 'D': sort D in decreasing order.
*
*  N       (input) INTEGER
*          The length of the array D.
*
*  D       (input/output) INTEGER array, dimension (N)
*          On entry, the array to be sorted.
*          On exit, D has been sorted into increasing order
*          (D(1) <= ... <= D(N) ) or into decreasing order
*          (D(1) >= ... >= D(N) ), depending on ID.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER   DIR, ENDD, I, J, START, STKPNT
      INTEGER   D1, D2, D3, DMNMX, TMP
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input paramters.
*
      INFO = 0
      DIR = -1
      IF (ID.eq.'D' .or. ID.eq.'d') THEN
         DIR = 0
      ELSEIF (ID.eq.'I' .or. ID.eq.'i') THEN
         DIR = 1
      ENDIF
      IF (DIR.EQ.-1 ) THEN
         INFO = -1
      ELSEIF (N.LT.0 ) THEN
         INFO = -2
      ENDIF
      IF (INFO.NE.0 ) THEN
         RETURN
      ENDIF
*
*     Quick return if possible
*
      IF (N.LE.1) RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF (ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on D( START:ENDD )
*
         IF (DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF (D( J ).GT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 30
                  ENDIF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF (D( J ).LT.D( J-1 ) ) THEN
                     DMNMX = D( J )
                     D( J ) = D( J-1 )
                     D( J-1 ) = DMNMX
                  ELSE
                     GO TO 50
                  ENDIF
   40          CONTINUE
   50       CONTINUE
*
         ENDIF
*
      ELSEIF (ENDD-START.GT.SELECT ) THEN
*
*        Partition D( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         D1 = D( START )
         D2 = D( ENDD )
         I = ( START+ENDD ) / 2
         D3 = D( I )
         IF (D1.LT.D2 ) THEN
            IF (D3.LT.D1 ) THEN
               DMNMX = D1
            ELSEIF (D3.LT.D2 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D2
            ENDIF
         ELSE
            IF (D3.LT.D2 ) THEN
               DMNMX = D2
            ELSEIF (D3.LT.D1 ) THEN
               DMNMX = D3
            ELSE
               DMNMX = D1
            ENDIF
         ENDIF
*
         IF (DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF (D( J ).LT.DMNMX )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF (D( I ).GT.DMNMX )
     $         GO TO 80
            IF (I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 60
            ENDIF
            IF (J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            ENDIF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF (D( J ).GT.DMNMX )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF (D( I ).LT.DMNMX )
     $         GO TO 110
            IF (I.LT.J ) THEN
               TMP = D( I )
               D( I ) = D( J )
               D( J ) = TMP
               GO TO 90
            ENDIF
            IF (J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            ENDIF
         ENDIF
      ENDIF
      IF (STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRT
*
      END
c-----------------------------------------------------------------------
      subroutine clipperN(xyzc,xyzn,DIR,xcl,ycl,zcl,del_list)
      include 'basics.inc'
      real xyzc(3),xyzn(3)
      real xcl(nelm,8),ycl(nelm,8),zcl(nelm,8)
      integer del_list(1)
      real dxyz(3)
      CHARACTER*1 DIR,YESNO
C
      Nvts = 4
      IF (IF3D) Nvts=8
c
C
C
C     Count number
      NUMDEL=0
      DO 101 JE=1,NEL
         DO 100  i=1,Nvts
            dxyz(1) = xcl(je,i)-xyzc(1)
            dxyz(2) = ycl(je,i)-xyzc(2)
            dxyz(3) = zcl(je,i)-xyzc(3)
            dn    = dxyz(1)*xyzn(1)+dxyz(2)*xyzn(2)+dxyz(3)*xyzn(3)
            if (dir.eq.'>') then
               if (je.le.nel .and. dn .gt. 0) then
                   numdel=numdel+1
                   del_list(numdel) = je
                   goto 101
               endif
            else
               if (je.le.nel .and. dn .lt. 0) then
                   numdel=numdel+1
                   del_list(numdel) = je
                   goto 101
               endif
            endif
  100    CONTINUE
  101 CONTINUE
c
      if (dir.eq.'>') then
         write(s,102) numdel,nel
  102    format(' You will be eliminating',i5,' of ',i5,
     $   ' elements ABOVE clipping plane.$')
      else
         write(s,103) numdel,nel
  103    format(' You will be eliminating',i5,' of ',i5,
     $   ' elements BELOW clipping plane.$')
      endif
      call prs(s)
      call prs(' OK? (Y/N)$')
      call res(yesno,1)
      if (yesno.ne.'y'.and.yesno.ne.'Y') return
      call prs('Deleting elements.$')
c
c
c     Sort into descending element order for "delel" routine
c
      call iusort('d',numdel,del_list,info)
      if (info.ne.0) then
         call prsii('Problem w/ sort:$',info,numdel)
         call prs('No elements deleted.$')
         return
      endif
c
      do i=1,numdel
         je = del_list(i)
         call delel(je)
      enddo
c
      IF (NUMDEL.GT.0) THEN
C        redraw the mesh
         CALL REFRESH
         CALL DRMENU('NOCOVER')
         CALL DRGRID
         write(s,500) numdel,nel
         call prs(s)
         DO 130 IEL=1,NEL
            CALL DRAWEL(IEL)
  130    CONTINUE
      ENDIF
  500 format(i7,' elements deleted.',i7,' elements remaining.$')
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine iusort( ID, N, list, INFO )
c
c     Return unique list of input length N, output length N=N'
c
      character*1 id
      integer list(1)
c
      if (n.le.1) return
      call ipsort (id,n,list,info)
c
      if (info.eq.0) then
c        compress list
         nold = n
         n = 1
         do i=2,nold
            if (list(i).ne.list(n)) then
               n=n+1
               list(n) = list(i)
            endif
         enddo
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine stretch_rad
      include 'basics.inc'
c
      logical if_sph_str
C
      call prs(' Input expansion factor ( =< 0 to abort):$')
      call rer(sfact)
      if (sfact.le.0) return
      call prs(' Spherical or Cylindrical (x-y) stretch? (s,c)$')
      call res(ans,1)
      if_sph_str = .true.
      if (ans.eq.'c' .or. ans.eq.'C') if_sph_str = .false.

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
      do 100 ie=1,nel
      do 100 i=1,nvts
         if (z(ie,i).gt.z00) then
c           zt = z(ie,i)-z0
            zt = z(ie,i)
            rt = zt*zt + x(ie,i)**2 + y(ie,i)**2
            if (rt.gt.rr) then
               if (if_sph_str) z(ie,i) = z0 + sfact*zt
               x(ie,i) = sfact*x(ie,i)
               y(ie,i) = sfact*y(ie,i)
               if (if_sph_str) z(ie,i) = sfact*z(ie,i)
            endif
         else
            rt = x(ie,i)**2 + y(ie,i)**2
            if (rt.gt.rr) then
               x(ie,i) = sfact*x(ie,i)
               y(ie,i) = sfact*y(ie,i)
            endif
         endif
  100 continue
C
C     Take care of curved sides
C
      do 200 ie = 1,nel
      do 200 is = 1,8
         if (ccurve(is,ie).eq.'s') then
C           spherical side, rad = curve(4,is,ie)
c           if (curve(4,is,ie).gt.re)
c    $         curve(4,is,ie) = sfact*curve(4,is,ie)
            do i=1,4
               curve(i,is,ie) = sfact*curve(i,is,ie)
            enddo
         endif
         if (ccurve(is,ie).eq.'C') then
C           cylindrical side, rad = curve(1,is,ie)
            if (abs(curve(1,is,ie)).gt.re)
     $         curve(1,is,ie) = sfact*curve(1,is,ie)
         endif
  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine clean_spheres
      include 'basics.inc'
      integer e,f
c
c     see routine crn3d ...
c
c     do e=1,nel
c     do f=1,6
c        if (ccurve(f,e).eq.'s') then
c           curve(1,f,e) = 
c           curve(2,f,e) = 
c           curve(3,f,e) = 
c        endif
c     enddo
c     enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rep_rot  !  replication (translation) rotation option
      include 'basics.inc'
c
      call prs('Rep. (1), Rot. (2), Rep/Rot (3) ? (0=abort)$')
      call res(ans,1)
c
      if (ans.eq.'1') call replicate_mesh
      if (ans.eq.'2'.and.ndim.eq.2) call rotate_mesh_2d
      if (ans.eq.'2'.and.ndim.eq.3) call rotate_mesh_3d
      if (ans.eq.'3') call rep_rotate_mesh
c
      return
      end
c-----------------------------------------------------------------------
      subroutine replicate_mesh
      include 'basics.inc'
c
      if (if3d) then
         call prs('Input translation vector: x,y,z$')
         call rerrr(xt0,yt0,zt0)
      else
         call prs('Input translation vector: x,y$')
         call rerr (xt0,yt0)
      endif
c
      call prs('Input number of reps (e.g., 1 --> double mesh size)$')
      call rei(nrep)
c
      do i=1,nrep
         ie0 = 1 + i*nel
         ie1 = ie0 + nel - 1
         call copy_sub_mesh(1,nel,ie0)
         xt = i*xt0
         yt = i*yt0
         zt = i*zt0
         call translate_sub_mesh(ie0,ie1,xt,yt,zt)
      enddo
      nel = ie1
      return
      end
c-----------------------------------------------------------------------
      subroutine rep_rotate_mesh
      include 'basics.inc'
c
      call prs('Input rotation angle (deg):$')
      call rer(angle_deg)
c
      call prs('Input number of reps (e.g., 1 --> double mesh size)$')
      call rei(nrep)
c
      ie0 = 1
      ie1 = nel
      ie2 = ie1+1
      ie3 = ie2 + nel-1
      do i=1,nrep

         ie2 = ie1+1
         ie3 = ie2 + nel-1

c        call copy_sub_mesh(ie0,ie1,ie2)
c        call rotate_submesh_2d(ie2,ie3,angle_deg)

         call copy_sub_mesh(1,nel,ie2)
         angle_deg_i = i*angle_deg
         call rotate_submesh_2d(ie2,ie3,angle_deg_i)

         ie0 = ie0 + nel
         ie1 = ie0 + nel-1

      enddo

      nel = ie3

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_sub_mesh(ie0,ie1,inew)
      include 'basics.inc'
c
      je = inew
      do ie=ie0,ie1
         call copyel(ie,je)
         je = je+1
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine translate_sub_mesh(e0,e1,xt,yt,zt)
      include 'basics.inc'
      integer e,e0,e1,f
c
      do e=e0,e1
c
         do i =1,8
            x(e,i) = x(e,i) + xt
            y(e,i) = y(e,i) + yt
            z(e,i) = z(e,i) + zt
         enddo
c
         do f=1,6
            if (ccurve(f,e).eq.'s') then
               curve(1,f,e) = curve(1,f,e) + xt
               curve(2,f,e) = curve(2,f,e) + yt
               curve(3,f,e) = curve(3,f,e) + zt
            endif
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_mesh_3d
      include 'basics.inc'
      real a(3,3)
      integer e,f
c
      call rzero(a,9)
      call prs('Opening file rot.mat$')
      open  (unit=10,file='rot.mat',status='old',err=999)
      read  (10,*) ((a(i,j),j=1,ndim),i=1,ndim)
      close (10)
c
      do e=1,nel
      do i=1,2**ndim
         xt = x(e,i)
         yt = y(e,i)
         zt = z(e,i)
         x(e,i) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
         y(e,i) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
         z(e,i) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
      enddo
      enddo
c
      do e=1,nel
      do f=1,6
         if (ccurve(f,e).eq.'s') then
            xt = curve(1,f,e)
            yt = curve(2,f,e)
            zt = curve(3,f,e)
            curve(1,f,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
            curve(2,f,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
            curve(3,f,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
         endif
      enddo
      enddo
c
      return
c
  999 call prs('File "rot.mat" not found...$')
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_mesh_2d
      include 'basics.inc'
      real a(3,3)
      integer e,f
c
      call prs('Input angle to rotate (in degrees):$')
      call rer(angle)
      angle = angle*pi/180.
      ca    = cos(angle)
      sa    = sin(angle)
c
      do e=1,nel
      do i=1,2**ndim
         xt = x(e,i)
         yt = y(e,i)
         x(e,i) = ca*xt-sa*yt
         y(e,i) = sa*xt+ca*yt
      enddo
      enddo

      call rzero(a,9)
      a(1,1) =  ca
      a(1,2) = -sa
      a(2,1) =  sa
      a(2,2) =  ca
      a(3,3) =  1.

      do e=1,nel
      do f=1,6
         if (ccurve(f,e).eq.'s') then
            xt = curve(1,f,e)
            yt = curve(2,f,e)
            zt = curve(3,f,e)
            curve(1,f,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
            curve(2,f,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
            curve(3,f,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
         endif
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_submesh_2d(e0,e1,angle_deg)
      include 'basics.inc'
      integer e,f,e0,e1
      real a(3,3)

      angle = pi*angle_deg/180.
      ca    = cos(angle)
      sa    = sin(angle)

      call rzero(a,9)
      a(1,1) =  ca
      a(1,2) = -sa
      a(2,1) =  sa
      a(2,2) =  ca
      a(3,3) =  1.

      do e=e0,e1
      do i=1,2**ndim
         xt = x(e,i)
         yt = y(e,i)
         x(e,i) = ca*xt-sa*yt
         y(e,i) = sa*xt+ca*yt
      enddo
      enddo

      do e=e0,e1
      do f=1,6
         if (ccurve(f,e).eq.'s') then
            xt = curve(1,f,e)
            yt = curve(2,f,e)
            zt = curve(3,f,e)
            curve(1,f,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
            curve(2,f,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
            curve(3,f,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
