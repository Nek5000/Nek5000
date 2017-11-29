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
#     include "basics.inc"
      include 'basicsp.inc'

      common /cisplit/ isplit(nelm)
      logical iftmp

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
              NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='SPIDER WEB'
              NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='PICTURE FRAME'
              NCHOIC=NCHOIC+1
         ITEM(NCHOIC)='CORNER FRAME'
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
c     if (ndim.eq.2) then
c             NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='Refine Hexagons'
c     else
c             NCHOIC=NCHOIC+1
c        ITEM(NCHOIC)='Hex transition'
c     endif

           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='OCT/Multi-SPLIT'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='Mesh Edit'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='SHIFT'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='STRETCH'
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='CLIP DOMAIN'
      if (ndim.eq.2) then
           NCHOIC=NCHOIC+1
        ITEM(NCHOIC)='SMOOTH'
      endif
           NCHOIC=NCHOIC+1
      ITEM(NCHOIC)='Clean up vertices'

      if (ndim.eq.2) call redraw_mesh
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
         call rep_rot
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
      end
c-----------------------------------------------------------------------
      subroutine frame(isplit,nframe)
C     Refines corner elements
#     include "basics.inc"
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
            rad=sqrt( (xmouse-x(ic,ifel(1)))**2
     $      +         (ymouse-y(ic,ifel(1)))**2 )
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
                 rad=sqrt((xmouse-x(ic,iel))**2 + (ymouse-y(ic,iel))**2)
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
         return
      ELSE IF(FAC.GT.1.0)THEN
        CALL PRS('Error; You cannot make new refined elements larger $')
        CALL PRS('than the old ones.  ABORTING.$')
         return
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
           DO 160 IEDGE=1,12
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
           if (icin(iiel).eq.iover(iiel)) then
             ifl = ifel(iiel)
             ne1 = nel+1
             ne2 = nel+2

C            FOR "TOP" element (inside corner# = overlapping side #)

             x(ic,ne1)=x(ic,ne1)+(x(icp1,ne1)-x(ic,ne1))*fac
             y(ic,ne1)=y(ic,ne1)+(y(icp1,ne1)-y(ic,ne1))*fac

C            Calculate Middle point
             xdiag=x(ic,ne1)
             ydiag=y(ic,ne1)
             xmid=xdiag+( x(icp2,ifl)-x(icp1,ifl) ) * fac
             ymid=ydiag+( y(icp2,ifl)-y(icp1,ifl) ) * fac

             x(icp3,ne1)=xmid
             y(icp3,ne1)=ymid

             x(ic,ne1)=x(ic,ne1)+(x(icp3,ne1)-x(ic,ne1))*fac
             y(ic,ne1)=y(ic,ne1)+(y(icp3,ne1)-y(ic,ne1))*fac
             x(icp1,ne2)=xmid
             y(icp1,ne2)=ymid

             x(icp1,ifl)=x(ic,ifl)+(x(icp1,ifl)-x(ic,ifl))*fac
             y(icp1,ifl)=y(ic,ifl)+(y(icp1,ifl)-y(ic,ifl))*fac
             x(icp2,ifl)=xmid
             y(icp2,ifl)=ymid
             x(icp3,ifl)=x(ic,ifl)+(x(icp3,ifl)-x(ic,ifl))*fac
             y(icp3,ifl)=y(ic,ifl)+(y(icp3,ifl)-y(ic,ifl))*fac
         else
C            FOR "other" element (inside corner# NOT= overlapping side #)
             x(ic,ne1)=x(ic,ne1)+(x(icp3,ne1)-x(ic,ne1))*fac
             y(ic,ne1)=y(ic,ne1)+(y(icp3,ne1)-y(ic,ne1))*fac

C            Calculate Middle point
             xdiag=x(ic,ne1)
             ydiag=y(ic,ne1)
             xmid=xdiag+( x(icp2,ifl)-x(icp3,ifl) ) * fac
             ymid=ydiag+( y(icp2,ifl)-y(icp3,ifl) ) * fac

             x(icp1,ne1)=xmid
             y(icp1,ne1)=ymid

             x(ic  ,ne2)=x(ic,ne2)+(x(icp1,ne2)-x(ic,ne2))*fac
             y(ic  ,ne2)=y(ic,ne2)+(y(icp1,ne2)-y(ic,ne2))*fac
             x(icp3,ne2)=xmid
             y(icp3,ne2)=ymid

             x(icp1,ifl)=x(ic,ifl)+(x(icp1,ifl)-x(ic,ifl))*fac
             y(icp1,ifl)=y(ic,ifl)+(y(icp1,ifl)-y(ic,ifl))*fac
             x(icp2,ifl)=xmid
             y(icp2,ifl)=ymid
             x(icp3,ifl)=x(ic,ifl)+(x(icp3,ifl)-x(ic,ifl))*fac
             y(icp3,ifl)=y(ic,ifl)+(y(icp3,ifl)-y(ic,ifl))*fac
         endif
390      continue

         nel=nel+2
       endif
400   continue
C     Recalculate centers
      IF(NFRAME.EQ.1)NBACK=1
      IF(NFRAME.EQ.2)NBACK=3
      do 110 iel=1,NEL
        xcen(iel)=(x(1,iel)+x(2,iel)+x(3,iel)+x(4,iel))/4.
        ycen(iel)=(y(1,iel)+y(2,iel)+y(3,iel)+y(4,iel))/4.
        call drawel(iel)
110   continue
      do 120 iiel=1,nframe
        iel=ifel(iiel)
        xcen(iel)=(x(1,iel)+x(2,iel)+x(3,iel)+x(4,iel))/4.
        ycen(iel)=(y(1,iel)+y(2,iel)+y(3,iel)+y(4,iel))/4.
        call drawel(iel)
120   continue
C
      return
      end
c-----------------------------------------------------------------------
      subroutine splitf
C     Splits floor in mesh
#     include "basics.inc"
      CALL PRS('Which floor do you wish to split?$')
      call rer(floor)
      ISPLIT=FLOOR
      IF(ISPLIT.GT.NLEVEL) THEN
         CALL PRS
     $   ('ERROR- ONLY ',NLEVEL,'Floors Exist.  Aborting Split.$')
         return
      ENDIF
      CALL PRS('At what height in this floor do you want the crack?$')
      CALL PRS('0<h<1; 0 is at floor, 1 is at ceiling.'//
     $'  negative to abort.$')
      call rer(fac)
      IF(FAC.GE.1.0)THEN
         CALL PRS('ERROR- Split cant be above ceiling$')
         return
      ELSE IF (FAC.LE.0.0)then
         CALL PRS('Aborting Floor Split$')
         return
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
              x(ic+4,iel)=x(ic  ,iel) + (x(ic+4,iel)-x(ic  ,iel)) * fac
              x(ic  ,nel)=x(ic+4,iel)
              y(ic+4,iel)=y(ic  ,iel) + (y(ic+4,iel)-y(ic  ,iel)) * fac
              y(ic  ,nel)=y(ic+4,iel)
              z(ic+4,iel)=z(ic  ,iel) + (z(ic+4,iel)-z(ic  ,iel)) * fac
              z(ic  ,nel)=z(ic+4,iel)
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
           xcen(nel)=vlsum(x(1,nel),4)/4.
           ycen(nel)=vlsum(y(1,nel),4)/4.
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
      return
      end
c-----------------------------------------------------------------------
      subroutine parker(isplit)
C     Makes spider web
#     include "basics.inc"
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
              xcent=vlsum(x(5,ielcrk),4)/4.
              ycent=vlsum(y(5,ielcrk),4)/4.
            ENDIF
            DO 200 IP=1,4
C              Straighten internal curved sides
               DO 160 IEDGE=1,12
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
               x    (ic,nel+ip)=x(iout,nel+ip)+
     $         (xcen(nel+ip) -  x(iout,nel+ip))*.65
               y    (ic,nel+ip)=y(iout,nel+ip)+
     $         (ycen(nel+ip) -  y(iout,nel+ip))*.65
               if (if3d) then
                  x(ic+4,nel+ip)=x(iout+4,nel+ip)+
     $            (xcent    -    x(iout+4,nel+ip))*.65
                  y(ic+4,nel+ip)=y(iout+4,nel+ip)+
     $            (ycent    -    y(iout+4,nel+ip))*.65
               endif
               if (inew.eq.2) then
C                 Move Corners of center element
                  x(ip,ielcrk)=x(ic,nel+ip)
                  y(ip,ielcrk)=y(ic,nel+ip)
                  if(if3d)then
                     x(ip+4,ielcrk)=x(ic+4,nel+ip)
                     y(ip+4,ielcrk)=y(ic+4,nel+ip)
                  endif
               endif
200         continue
C           Recalculate centers
            do iel=nel+1,NEL+4
              xcen(iel)=.25*vlsum(x(1,iel),4)
              ycen(iel)=.25*vlsum(y(1,iel),4)
            enddo
            iel=ielcrk
            xcen(iel)=.25*vlsum(x(1,iel),4)
            ycen(iel)=.25*vlsum(y(1,iel),4)
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
      return
      end
c-----------------------------------------------------------------------
      subroutine crack(ielcrk,ifacrk,sfrac)
C     Initiates crack in 1st element of split.  Output is Element #, Side#, and
C     Fraction away from side.
#     include "basics.inc"
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
      return
      end
c-----------------------------------------------------------------------
      subroutine mark (ielcrk,ifacrk,sfrac,isplit)
C     Mark where CRACK propagates in vector ISPLIT
#     include "basics.inc"
      INTEGER ISPLIT(NELM)

c     write(6,*) ielcrk,ifacrk,sfrac,nel,' Cracking mark'
      call gen_neigh

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
                                rad5=sqrt((x(marked,1)-x(ic+4,iel))**2+
     $                                    (y(marked,1)-y(ic+4,iel))**2)
                                rad6=sqrt((x(marked,5)-x(ic  ,iel))**2+
     $                                    (y(marked,5)-y(ic  ,iel))**2)
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
c                              call prsii('isplit,nsplit$'
c    $                             ,isplit(iel),nsplit)

                               call prs('Continue splitting?$')
                               call res(ans,1)

                               if (.not.ifgraf) goto 301
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
      return
      end
c-----------------------------------------------------------------------
      subroutine mark2(ielcrk,ifacrk,sfrac,isplit)
C     Mark where CRACK propagates in vector ISPLIT
#     include "basics.inc"
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
#     include "basics.inc"
      integer isplit(1000),iipoint(2,3,2),aface
      character letnew(1000)
      logical if_special_curve

c     IPOINT: 1st subscript-which end of a; 2nd subscript- a or b; 3rd-top/bot
c     First we modify element.  Move oposite face toward IFACRK.  Insert
c     new element in gap.

      maxnw=maxlet  ! Set up stuff for letter of daughter element
      do 20 i=1,maxlet
         letnew(i)=' '
20    continue

c     write(6,222) (isplit(k),k=1,nel)
c 222 format(' isplit: ',20i4)

      ia=1
      ib=2
      nelold=nel
      do 100 iel=1,nelold
         if (isplit(iel).ne.0) then 
           if_special_curve = .false.
           nedge = 4 + 8*(ndim-2)
           do ic=1,nedge
              if (ccurve(ic,iel).eq.'s') if_special_curve = .true.
              if (ccurve(ic,iel).eq.'m') if_special_curve = .true.
c             write(6,*) iel,ic,ccurve(ie,iel),if_special_curve
           enddo 
           if (if_special_curve) then
              aface = isplit(iel)
              ratio = (1-sfrac)/sfrac
              ratii = 1./ratio
c             write(6,*) iel,aface,sfrac,ratio,nel,' MSPLIT$'
c             write(s,*) iel,aface,sfrac,ratio,nel,' MSPLIT$'
c             call prs(s)
              if (aface.eq.1) call msplite(iel,1,2,1,1.,ratio,1.)
              if (aface.eq.2) call msplite(iel,2,1,1,ratii,1.,1.)
              if (aface.eq.3) call msplite(iel,1,2,1,1.,ratii,1.)
              if (aface.eq.4) call msplite(iel,2,1,1,ratio,1.,1.)
           else
            write(s,'(1x,a18,I9,a6,I9,I9)')
     $      'SPLITTING ELEMENT ',IEL,'  SIDE',ISPLIT(I),NEL
            call prs(s//'$')
            nel=nel+1
            call copyel(iel,nel)
c           Give Daughter element a letter
            if(letnew(ichar(letapt(iel))) .eq. ' ')then
c              A parent with this letapt hasn't split yet; assign new pointer
               maxnw=maxnw+1
               letnew(ichar(letapt(iel))) = char(min0(122,maxnw))
            endif
            letapt(nel) = letnew(ichar(letapt(iel)))
c           That was, Each old letter points to a new letter for the daughter
c           Establish pointers to the element corners
            iipoint(1,ia,1)=isplit(iel)
            iipoint(2,ia,1)=isplit(iel)+1
            if(iipoint(2,ia,1).eq.5) iipoint(2,ia,1)=1

            iipoint   (2,ib,1)     = iipoint(2,ia,1)+1
            if(iipoint(2,ib,1).eq.5) iipoint(2,ib,1)=1
            iipoint   (1,ib,1)     = iipoint(1,ia,1)-1
            if(iipoint(1,ib,1).eq.0) iipoint(1,ib,1)=4

            do 50 i12=1,2
            do 50 iabc=1,3
               iipoint(i12,iabc,2) = iipoint(i12,iabc,1) + 4
50          continue
C           Now move corners
            do 60 i12=1,2
            do 60 itb=1,2
              ica=iipoint(i12,ia,itb)
              icb=iipoint(i12,ib,itb)
              if (i12.eq.1) iedge=icb
              if (i12.eq.2) iedge=ica
C             Move corners of new element
              if (ccurve(iedge,iel).eq.' ') then  ! straight line
                 x(ica,nel)=x(ica,nel) + (x(icb,nel)-x(ica,nel))*sfrac
                 y(ica,nel)=y(ica,nel) + (y(icb,nel)-y(ica,nel))*sfrac
                 z(ica,nel)=z(ica,nel) + (z(icb,nel)-z(ica,nel))*sfrac
              else
                if(i12.eq.2)
     $          call getpts(1,    sfrac,iel,iedge,x(ica,nel),y(ica,nel))
                if(i12.eq.1)
     $          call getpts(1,1.0-sfrac,iel,iedge,x(ica,nel),y(ica,nel))
              endif
              z(ica,nel)=z(ica,nel) + (z(icb,nel)-z(ica,nel))*sfrac
C             Move corners of old element
              x(icb,iel)=x(ica,nel)
              y(icb,iel)=y(ica,nel)
              z(icb,iel)=z(ica,nel)
              if(ccurve(iedge,iel).eq.'S')then ! Modify Control points
                 if (i12.eq.2)then
                    curve(3,iedge,iel)=x(icb,nel)
                    curve(4,iedge,iel)=y(icb,nel)
                    curve(1,iedge,nel)=x(ica,iel)
                    curve(2,iedge,nel)=y(ica,iel)
                 elseif (i12.eq.1)then
                    curve(3,iedge,nel)=x(ica,iel)
                    curve(4,iedge,nel)=y(ica,iel)
                    curve(1,iedge,iel)=x(icb,nel)
                    curve(2,iedge,iel)=y(icb,nel)
                 endif
              endif
60          continue
cuuuu
cuuuu
cuuuu  Here, we need to add support for midside nodes
cuuuu
cuuuu  Easiest way is to check for any 'm' curves, then
c      generate x27, etc. for that element only (or, "x125" as it were)
c      or... I guess it depends on "s" ... and edge.... so better think
cuuuu
cuuuu
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
                       ibc(is,nel,if)=0
                  ENDIF
70          CONTINUE

         ENDIF
         ENDIF
100   CONTINUE
C     Recalculate centers
      DO 110 IEL=1,NEL
        xcen(iel)=(x(1,iel)+x(2,iel)+x(3,iel)+x(4,iel))/4.
        ycen(iel)=(y(1,iel)+y(2,iel)+y(3,iel)+y(4,iel))/4.
110   CONTINUE
C
      CALL SORTEL
C     Draw whole new isometric surface
      DO 51 IEL=1,NEL
c        CALL DRAWIS(ISRT(IEL))
         IF(.NOT.IF3D  .OR. NUMAPT(IEL).EQ.ILEVEL) CALL DRAWEL(IEL)
51    CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine vertadj

#     include "basics.inc"

      integer icalld,nelold,ff,e
      save    icalld,nelold
      data    icalld,nelold /2*0/

      icalld=icalld+1
      io = icalld+75

      call gencen
      rmin=glmin(rcen,nel)

      call cell_cell_connectivity

      WRITE(S,10) RMIN
   10 FORMAT(2X,'RMIN:',g13.6,'$')
      CALL PRS(S)

c     CALL PRS('Input tolerance, relative to el. size (e.g., 0.5:)$')
c     CALL RER(EPS)
      eps = 0.05


      WRITE(S,20) EPS,RMIN
   20 FORMAT(2X,'EPS,RMIN:',2E13.6,'$')
      CALL PRS(S)
      EPS2=EPS**2
      epsr=eps*rmin


      ncrnr=2**ndim
      ichk=10
      if (nel.gt.200) ichk=50
      if (nel.gt.2000) ichk=500
      if (nel.gt.20000) ichk=5000
      iadj=0

      if (nelold.ne.nel) then
         nelold=nel
         do 500 ie=1,nel
         do 500 ic=1,ncrnr
            x(ic,ie)=round(x(ic,ie))
            y(ic,ie)=round(y(ic,ie))
            z(ic,ie)=round(z(ic,ie))
  500    continue
      endif

      epsr=.0001*rmin
      do 1000 ie=1,nel

         if (mod(ie,ichk).eq.0.and.iadj.eq.0) then
            write(6,1001) ie
         else
            iadj=0
         endif

         do 100 ic=1,ncrnr

            if (abs(x(ic,ie)).lt.epsr.and.x(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing x:',ic,ie,x(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               x(ic,ie)=0.0
               iadj=1
            endif
C
            if (abs(y(ic,ie)).lt.epsr.and.y(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing y:',ic,ie,y(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               y(ic,ie)=0.0
               iadj=1
            endif
C
            if (abs(z(ic,ie)).lt.epsr.and.z(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing z:',ic,ie,z(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               z(ic,ie)=0.0
               iadj=1
            endif
  100    continue

         nfaces = 2*ndim
         do ff=1,nfaces
           if (cbc(ff,ie,1).eq.'E  ') then
            je = ibc(ff,ie,1)
            dist1=(xcen(ie)-xcen(jE))**2+(ycen(ie)-ycen(jE))**2
     $           +(zcen(ie)-zcen(jE))**2
            dist2= 1.01*( rcen(ie) + rcen(je) )**2
            if (dist1.le.dist2) then
               epsrr = eps2*min(rcen(ie),rcen(je))
               epsrr = epsrr**2
               do 700 ic=1,ncrnr
               do 700 jc=1,ncrnr
                  dist1=(x(ic,ie)-x(jc,je))**2+(y(ic,ie)-y(jc,je))**2
     $                 +(z(ic,ie)-z(jc,je))**2
                  if (dist1.lt.epsrr.and.dist1.ne.0.0) then
                     if (nel.lt.2000) then
                        write(s,'(1x,a10,2(i6,i3))')
     $                  'Adjusting:',je,jc,ie,ic
                        CALL PRS(S//'$')
                        call prrrr(x(ic,ie),y(ic,ie),z(ic,ie))
                        call prrrr(x(jc,je),y(jc,je),z(jc,je))
                     endif
c
                     x(jc,je)=x(ic,ie)
                     y(jc,je)=y(ic,ie)
                     z(jc,je)=z(ic,ie)
                     iadj=1
                  endif
  700          continue
            endif
           endif
         enddo
 1000 continue
 1001 format('  Checking',i9)
      return
      end
c-----------------------------------------------------------------------
      subroutine vertadj_old
C
#     include "basics.inc"
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
      if (nel.gt.2000) ichk=500
      if (nel.gt.20000) ichk=5000
      iadj=0
C
      if (nelold.ne.nel) then
         nelold=nel
         do 500 ie=1,nel
         do 500 ic=1,ncrnr
            x(ic,ie)=round(x(ic,ie))
            y(ic,ie)=round(y(ic,ie))
            z(ic,ie)=round(z(ic,ie))
  500    continue
      endif
C
      epsr=.0001*rmin
      do 1000 ie=1,nel
C
         if (mod(ie,ichk).eq.0.and.iadj.eq.0) then
            write(6,1001) ie
         else
            iadj=0
         endif
C
         do 100 ic=1,ncrnr
C
            if (abs(x(ic,ie)).lt.epsr.and.x(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing x:',ic,ie,x(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               x(ic,ie)=0.0
               iadj=1
            endif
C
            if (abs(y(ic,ie)).lt.epsr.and.y(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing y:',ic,ie,y(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               y(ic,ie)=0.0
               iadj=1
            endif
C
            if (abs(z(ic,ie)).lt.epsr.and.z(ic,ie).ne.0.0) then
               write(s,'(1x,a10,i6,i4,g16.8)') 
     $         'zeroing z:',ic,ie,z(ic,ie)
               if (nel.lt.10000.or.mod(ie,1000).eq.0) CALL PRS(S//'$')
               z(ic,ie)=0.0
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
               do 700 ic=1,ncrnr
               do 700 jc=1,ncrnr
                  dist1=(x(ic,ie)-x(jc,je))**2+(y(ic,ie)-y(jc,je))**2
     $                 +(z(ic,ie)-z(jc,je))**2
                  IF (DIST1.LT.epsrr.and.dist1.ne.0.0) THEN
                     if (nel.lt.2000) then
                        write(s,'(1x,a10,2(i6,i3))')
     $                  'Adjusting:',je,jc,ie,ic
                        CALL PRS(S//'$')
c
c                   diag.
c                      if (ie.eq.514.and.je.eq.3377) then
c                         write(6,*) ic,x(ic,ie),y(ic,ie),z(ic,ie),dist1
c                         write(6,*) jc,x(jc,ie),y(jc,ie),z(jc,ie),eps2
c                         write(6,*) 'go?'
c                         call res(s,1)
c                      endif
c                      write(io,75) ic,ie,x(ic,ie),y(ic,ie),z(ic,ie),je,jc
c                      write(io,76) je,jc,x(jc,je),y(jc,je),z(jc,je),ic,ie
c                   diag.
c
                        call prrrr(x(ic,ie),y(ic,ie),z(ic,ie))
                        call prrrr(x(jc,je),y(jc,je),z(jc,je))
                     endif
c
                     x(jc,je)=x(ic,ie)
                     y(jc,je)=y(ic,ie)
                     z(jc,je)=z(ic,ie)
C
c                diag.
c                  write(io,77) jc,je,x(jc,je),y(jc,je),z(jc,je),ic,ie
c                  call prrrr(x(jc,je),y(jc,je),z(jc,je))
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
 1001 FORMAT('  Checking',i9)
      return
      end
c-----------------------------------------------------------------------
      subroutine clpdom
#     include "basics.inc"
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
      ITEM(nchoic)       =             'Clip R'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Input DEL list'
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
C
      IF (CHOICE.EQ.'UP MENU') return
      IF (CHOICE.EQ.'Clip X') THEN
         CALL PRS(
     $   'Input location of X-clipping plane.$')
         CALL RER(Xclip)
c        CALL DRAWLINE(Xclip,Ymax,Xclip,Ymin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
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
         IF (ANS.eq.'=') return
         CALL Clipper(Yclip,ANS,Y)
      ELSEIF (CHOICE.EQ.'Clip Z') THEN
         CALL PRS(
     $   'Input location of Z-clipping plane.$')
         CALL RER(Zclip)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL Clipper(Zclip,ANS,Z)
      ELSEIF (CHOICE.EQ.'Clip R') THEN
         CALL PRS
     $   ('Type in center and radius of cylinder cut (e.g., 0 0 .5):$')
         call rerrr(xcyl,ycyl,rcyl)
         call prs(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         call clipper_r(xcyl,ycyl,rcyl,ans)
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

         call normalize(xyznorm,xyzl,3)

         CALL PRS(
     $   'Input "<" or ">" to indicate desired clip section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL ClipperN(xyzclip,xyznorm,ANS,x,y,z,del_list)
      ELSEIF (CHOICE.EQ.'Input DEL list') THEN
         call mult_del(del_list)
      ENDIF
      GOTO 1
      end
c-----------------------------------------------------------------------
      subroutine clipper(Clip,DIR,pts)
#     include "basics.inc"
      real pts(8,nelm)
      character*1 dir,yesno
      common /ctmp0/ idel(nelm),edel(nelm)
      integer edel,slot,e

      call izero(idel,nel)

      nvts = 4
      if (if3d) nvts=8

      numdel=0
      do e=1,nel
         do i=1,nvts
            if ((dir.eq.'>'.and.pts(i,e).gt.clip)  .or.
     $          (dir.eq.'<'.and.pts(i,e).lt.clip)) then
                numdel=numdel+1
                idel(e)=1
                edel(numdel)=e
                goto 100
            endif
         enddo
  100    continue
      enddo

      if (dir.eq.'>') write(s,101) numdel,nel
  101 format(' You will be eliminating',i8,' of ',i8,
     $' elements ABOVE clipping plane.$')

      if (dir.eq.'<') write(s,102) numdel,nel
  102 format(' You will be eliminating',i8,' of ',i8,
     $' elements BELOW clipping plane.$')

      call prs(s)
      call prs(' OK? (Y/N)$')
      call res(yesno,1)
      if (yesno.eq.'n'.or.yesno.eq.'N') return
      if (numdel.eq.0) return

      slot = edel(1)
      do e=edel(1)+1,nel
         if (idel(e).eq.0) then
            call copyel(e,slot)  ! e-->slot
            slot = slot+1
         endif
      enddo
      nel = slot-1

      call copyel(nel,nelm)

C     Recount the number of curved sides
      ncurve=0
      do e=1,nel
      do iedge=1,12
         if(ccurve(iedge,e).ne.' ') ncurve=ncurve+1
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine clipper_r(xcyl,ycyl,rcyl,dir)
#     include "basics.inc"
      character*1 dir,yesno
      common /ctmp0/ idel(nelm),edel(nelm)
      integer edel,slot,e

      call izero(idel,nel)

      nvts = 4
      if (if3d) nvts=8

      rcyl2 = rcyl**2

      numdel=0
      do e=1,nel
         rr=0
         do i=1,nvts
            rr = rr+(x(i,e)-xcyl)**2+(y(i,e)-ycyl)**2
         enddo
         rr = rr/nvts
         if ((dir.eq.'>'.and.rr.gt.rcyl2)  .or.
     $       (dir.eq.'<'.and.rr.lt.rcyl2) ) then
                numdel=numdel+1
                idel(e)=1
                edel(numdel)=e
                goto 100
         endif
  100    continue
      enddo

      if (dir.eq.'>') write(s,101) numdel,nel
  101 format(' You will be eliminating',i8,' of ',i8,
     $' elements OUTSIDE the clipping radius.$')

      if (dir.eq.'<') write(s,102) numdel,nel
  102 format(' You will be eliminating',i8,' of ',i8,
     $' elements INSIDE the clipping radius.$')

      call prs(s)
      call prs(' OK? (Y/N)$')
      call res(yesno,1)
      if (yesno.eq.'n'.or.yesno.eq.'N') return
      if (numdel.eq.0) return

      slot = edel(1)
      do e=edel(1)+1,nel
         if (idel(e).eq.0) then
            call copyel(e,slot)  ! e-->slot
            slot = slot+1
         endif
      enddo
      nel = slot-1

      call copyel(nel,nelm)

C     Recount the number of curved sides
      ncurve=0
      do e=1,nel
      do iedge=1,12
         if(ccurve(iedge,e).ne.' ') ncurve=ncurve+1
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine stretch
      call shift2
      return
      end
c-----------------------------------------------------------------------
      subroutine smoother2(ncell)
c
#     include "basics.inc"
      include 'basicsp.inc'
      real ww(0:1)
      integer eln(4,ncell),i,ncell
      real x4(4,ncell),y4(4,ncell)
      real x4b(4,ncell),y4b(4,ncell)
      real bndr(ncell*4),shrf,nodc(ncell*4)
      real bndrl(ncell*4),nodcl(ncell*4)
      real d(4,ncell)
      integer nv,nouter,nlap,nopt,nnpts
      parameter(lpts=8*nelm)
      common /arrayi/ i_n(lpts)   , j_n(4*lpts)
     $              , j_o(4*lpts)
     $              , ee(lpts)   , cell(lpts)
      integer cell,ee

      common /arrayr/ wk (lpts),dx(4*lpts)

      integer clist_i(lpts),clist_j(4*lpts)
      equivalence (clist_i,i_n)
      equivalence (clist_j,j_n)
 
      nvc = 2**ndim
      call set_d2a(dx,x,y,z,nvc,ncell,ndim)
      call makecell  (eln,nvc,ncell,dx,ndim,wk,i_n,j_n,j_o)
      call chker2    (eln,nvc,ncell)

      n=0
      nv = 0
      do i=1,nel
      do j=1,4
        nv = max(nv,eln(j,i))
      enddo
      enddo
      write(6,*) 'Element connectivity formed'

      nfaces=2**ndim
      call rone(bndr,nv)
      call find_bc_crv2(bndr,nv,eln,nvc,ncell)

      write(6,*) 'Smoother - boundary points identified'

c     eln is preproc notation right now
c     b is based on global node numbering from map file.. the number of 
c     values in b is equivalent to nv

c     so far the boundary points and elnod matrix have been constructed
c     need to be able to construct the Q matrix

c     Start generation Q
c     Q has nel*4 rows and nv columns

      call gennodc(nodc,eln,nv,ncell)      
      nvv = nv 
      write(6,*) 'Starting smoother'

      call prs ('Input the number of iterations outer,lap,opt:$')
      call prs ('suggested values are 20,10,10:$')
      call reiii(nouter,nlap,nopt)

      call xtox4(x4,y4,ncell,0)
      call copy(x4b,x4,ncell*4)
      call copy(y4b,y4,ncell*4)
      nnpts = 4
      call globtoloc(eln,nvv,ncell,nodcl,nodc,nnpts)
      call globtoloc(eln,nvv,ncell,bndrl,bndr,nnpts)

      call cheap_dist(d,1,'W  ',2,eln,nvv,ncell,nnpts)
      call disfun(d,0.1,nnpts,ncell)

      do i=1,nouter
       call runlapsm2d(x4,y4,eln,Q,nvv,bndrl,nodcl,d,ncell,nlap)
       call runoptsm2d(x4,y4,eln,Q,nvv,bndrl,nodcl,d,ncell,1,nopt)
      enddo

      call restbndrlay(x4,y4,x4b,y4b,d,ncell)
      call xtox4(x4,y4,ncell,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine 
     $       runlapsm2d(x4,y4,eln,Q,nvv,bndrl,nodcl,d,ncell,itmax)
#     include "basics.inc"
      integer nvv,ncell
      integer eln(4,ncell)
      real x4(4,ncell),y4(4,ncell)
      real bndrl(ncell*4),nodcl(ncell*4)
      real xa(4,ncell),ya(4,ncell)
      real dxl(ncell*4),xpm,ypm,dyl(ncell*4),d(4,ncell)
      real dxg(nvv),shrf
      integer n1,n2,i,j,k,e,itmax,nnpts

      nnpts = 4
      shrf = 0.99

      n1 = ncell*4 
      call rzero(dxl,n1)
      call rzero(dyl,n1)
      call rzero(dxg,nvv)

      do k=1,itmax
      call rzero(xa,n1)
      call rzero(ya,n1)
      n1 = 0
      do e=1,ncell
        xpm = vlsum(x4(1,e),4)/4.
        ypm = vlsum(y4(1,e),4)/4.
        do i=1,4
          n1 = n1+1
          dxl(n1) = xpm+shrf*(x4(i,e)-xpm)-x4(i,e) 
          dyl(n1) = ypm+shrf*(y4(i,e)-ypm)-y4(i,e)
        enddo
      enddo

      n1 = ncell*4
       call col2(dxl,bndrl,n1)
       call col2(dxl,d,n1)
       call qqtavg(eln,nvv,ncell,dxl,nodcl,nnpts)
       call add2(x4,dxl,n1)

       call col2(dyl,bndrl,n1)
       call col2(dyl,d,n1)
       call qqtavg(eln,nvv,ncell,dyl,nodcl,nnpts)
       call add2(y4,dyl,n1)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine runoptsm2d
     $      (x4,y4,eln,Q,nvv,bndrl,nodcl,d,ncell,opt,itmax)
#     include "basics.inc"
      integer nvv,ncell,opt
      integer eln(4,ncell),d(4,ncell)
      real bndrl(ncell*4),x4(4,ncell),y4(4,ncell)
      real nodcl(ncell*4),xa(4,ncell),ya(4,ncell)
      real scalek,hval
      real dxl(ncell*4),dyl(ncell*4),dval(4,ncell)
      integer n1,n2,i,j,k,e,itmax,iter,nnpts
      real f2,dfdx(4,ncell,2)

      n1 = 2**2
      n2 = ncell*n1

      nnpts = 4
 

      call get_nodscale(scalek,x4,y4,ncell)
      hval = scalek
      siz = 4 !same for jacobian and len^2 for 2D :) 

      do iter=1,itmax
        call gradf2d(f2,dfdx,x4,y4,siz,opt,hval,ncell,nvv,eln)
        dumc1 =  glamax(dfdx(1,1,1),2*n2)

        do e=1,ncell
        do j=1,nnpts
           dval(j,e) = max(abs(dfdx(j,e,1)),abs(dfdx(j,e,2)))
           if (dval(j,e).eq.0) dval(j,e) = (1.e+6)*(dumc1)
        enddo
        enddo

        call invcol1(dval,n2)
        call cmult(dval,scalek,n2)
        one = -1.
        call cmult(dval,one,n2)

        call col2(dfdx(1,1,1),bndrl,n2)
        call col2(dfdx(1,1,2),bndrl,n2)
  
        call col2(dfdx(1,1,1),d,n2)
        call col2(dfdx(1,1,2),d,n2)

        call add2col2(x4,dfdx(1,1,1),dval,n2)
        call add2col2(y4,dfdx(1,1,2),dval,n2)

      enddo
        write(6,*) f2,'globsum'

      return
      end
c-----------------------------------------------------------------------
      subroutine gradf2d(f1,dfdx,x4,y4,siz,opt,hval,ncell,nvv,eln)
#     include "basics.inc"
      integer n1,n2,i,j,k,e,opt,siz,nvv
      integer eln(4,ncell)
      real par(siz)
      real f1,f2,dfdx(4,ncell,2),fl
      real x4(4,ncell),y4(4,ncell)
      real xt(4),yt(4)
      integer nnpts

      nnpts = 2**2

      f1 = 0
      do e=1,ncell
       if (opt.eq.1) call get_jac2d(fl,x4(1,e),y4(1,e),siz)
       f1 = f1+fl
       call copy(xt,x4(1,e),4)
       call copy(yt,y4(1,e),4)
       do j=1,4
         xt(j) = x4(j,e)+hval
         if (opt.eq.1) call get_jac2d(f2,xt,yt,siz)
         dfdx(j,e,1) = (f2-fl)/hval
         xt(j) = x4(j,e)
         
         yt(j) = y4(j,e)+hval
         if (opt.eq.1) call get_jac2d(f2,xt,yt,siz)
         dfdx(j,e,2) = (f2-fl)/hval
         yt(j) = y4(j,e)
       enddo
      enddo

      call qqt(eln,nvv,ncell,dfdx(1,1,1),nnpts) 
      call qqt(eln,nvv,ncell,dfdx(1,1,2),nnpts) 


      return
      end
c-----------------------------------------------------------------------
      subroutine get_jac2d(val,xx,yy,siz)
#     include "basics.inc"
      real fl,xx(4),yy(4)
      real fr1,fr2,jac(4),jm(2,2),jin(2,2)
      real frn(4),sum1,par(4)
      real val
      integer siz,i,j,ind1,ind2
      integer bzindx(24),czindx(24)
      SAVE bzindx
      DATA bzindx / 2,3,5, 1,4,6, 4,1,7, 3,2,8,
     $             6,7,1, 5,8,2, 8,5,3, 7,6,4 /
c     bzindx tells which node is node connected to in r,s,t direction
c     example: node 1 is connected to 2,3,5; 2 to 1,4,6 and so on
      SAVE czindx
      DATA czindx / 1,1,1,  -1,1,1, 1,-1,1, -1,-1,1,
     $              1,1,-1, -1,1,-1, 1,-1,-1, -1,-1,-1 /


      do i=1,4
        ind1=(i-1)*3
        do j=1,2
         ind2 = ind1+j
         jm(1,j) = 0.5*czindx(ind2)*(xx(bzindx(ind2))-xx(i))
         jm(2,j) = 0.5*czindx(ind2)*(yy(bzindx(ind2))-yy(i))
        enddo
        jac(i)=jm(1,1)*jm(2,2)-jm(2,1)*jm(1,2)
        jin(1,1) = jm(2,2)
        jin(2,1) = -jm(2,1)

        jin(1,2) = -jm(1,2)
        jin(2,2) = jm(1,1)
 
        dumc = 1./jac(i)
        call cmult(jin,dumc,4)

        call rzero(frn,4)
        call col3(frn,jm,jm,4) !square the entries
        sum1 = vlsum(frn,4)
        fr1 = SQRT(sum1)           !squareroot

        call rzero(frn,4)
        call col3(frn,jin,jin,4)
        sum1 = vlsum(frn,4)
        fr2 = SQRT(sum1)

        par(i) = fr1*fr2
        par(i) = (par(i)/2)**2

      enddo

      val = vlsum(par,4)
      val = val/4

      return
      end
c-----------------------------------------------------------------------
      subroutine get_nodscale(scalek,x4,y4,ncell)
#     include "basics.inc"
      integer ncell
      real x4(4,ncell),y4(4,ncell)
      real scalek
      integer i,e,nod1,nod2
      real dlmin,lv,dx,dy
      integer efc2(2,4)
      data    efc2 / 1,2
     $             , 2,3
     $             , 3,4
     $             , 4,1 /

      dlmin = 1.e+7
      do e=1,ncell
      do i=1,4
         nod1 = efc2(1,i)
         nod2 = efc2(2,i)
         dx = x4(nod1,e)-x4(nod2,e)
         dy = y4(nod1,e)-y4(nod2,e)
         lv = (dx**2+dy**2)**0.5
         if (lv.lt.dlmin) dlmin=lv 
      enddo
      enddo

      scalek = 0.01*dlmin
      

      return
      end
c-----------------------------------------------------------------------
      subroutine xtox4(x4,y4,ncell,flag)
#     include "basics.inc"
      integer ncell,flag
      real x4(4,ncell),y4(4,ncell)
      integer i,e
c  if flag is 0, x to x4, otherwise x4 to x

      if (flag.eq.0) then
       do e=1,ncell
       do i=1,4
            x4(i,e) = x(i,e)
            y4(i,e) = y(i,e)
       enddo
       enddo
      else
       do e=1,ncell
       do i=1,4
            x(i,e) = x4(i,e)
            y(i,e) = y4(i,e)
       enddo
       enddo
      endif
       
      return
      end
c-----------------------------------------------------------------------
      subroutine gennodc(nodc,eln,nv,ncell)
#     include "basics.inc"
      integer eln(4,ncell)
      integer nv,ncell,n1,n2,n3,i,k
      real nodc(nv)

      n1 = ncell*(4)
      call izero(nodc,nv)
      n2 = 0
      do i=1,ncell
      do k=1,2**ndim
         n2 = n2+1
         j1 = eln(k,i)
         n3 = (j1-1)*n1      !this puts at the end of previous column 
         nodc(j1) = nodc(j1)+1.
      enddo
      enddo 
      write(6,*) 'Q matrix generated'
       

      return
      end
c-----------------------------------------------------------------------
      subroutine find_bc_crv2(b,nv,eln,nvc,ncell)
c
#     include "basics.inc"
      integer eln(nvc,ncell)
      real b(nv)
c
      integer efc2(2,4)
      save    efc2
      data    efc2 / 1,2
     $             , 2,3
     $             , 3,4
     $             , 4,1 /
c
c
c     Identify any elements on the boundary or having curved sides
c
      call rone(b,nv)
c
      nfaces = 2*ndim
      ifld = 2
      if (ifflow) ifld = 1
c
      do ie=1,ncell
      do is=1,nfaces
c         write(6,*) is,ie,'k10faceel'
         if (cbc(is,ie,ifld).ne.'E  ') then
            do iv=1,2
               b(eln(efc2(iv,is),ie)) =  0*b(eln(efc2(iv,is),ie))
            enddo
         endif
       enddo
       enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine smoother(cell,nvc,ia,ja,ww,b,g,x0,x1,y0,y1,z0,z1)
c
#     include "basics.inc"
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
        xo=x(i,ie)
        if (b(cell(i,ie)).le.0) then
           x(i,ie) = xp(cell(i,ie))
           y(i,ie) = yp(cell(i,ie))
           z(i,ie) = zp(cell(i,ie))
        endif
        jc=cell(i,ie)
        write(66,*) b(jc),' xo:',xo,x(i,ie),i,ie,jc,g(jc)
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
#     include "basics.inc"
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
c         if (cbc(is,ie,ifld).eq.'SYM') then
c           Later this will be fixed to handle each spatial dimension
c            do iv=1,ncrnf
c               b(cell(efc(iv,is),ie)) = -2
c            enddo
c         endif
c
c        Check for curve sides (spherical in particular)
c         if (ccurve(is,ie).eq.'s') then
c            do iv=1,ncrnf
c               b(cell(efc(iv,is),ie)) =  1
c            enddo
c         endif
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
#     include "basics.inc"
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
            xp(l) = x(k,ie)
            yp(l) = y(k,ie)
            zp(l) = z(k,ie)
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
#     include "basics.inc"
      include 'basicsp.inc'
      real xback(4,nelm),yback(4,nelm)
      integer ncell
      integer smflag
      SAVE smflag
      data smflag /0/
c
      x0 = -9.e19
      y0 = -9.e19
      z0 = -9.e19
      x1 =  9.e19
      y1 =  9.e19
      z1 =  9.e19
C
    1 CONTINUE
c      call xtox4(xback,yback,nel,0)  !creates a backup of the mesh
c
      nchoic = 0
      nchoic = nchoic+1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Redraw mesh'
c      nchoic = nchoic+1
c      ITEM(nchoic)       =             'Set smoothing box'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Smooth'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Undo'
c
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
c
C
      IF (CHOICE.EQ.'UP MENU') return
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
       smflag = 1
        call xtox4(xback,yback,nel,0)
        ncell = nel
        call convert_m_to_c_allr(10000000.)
        call smoother2(ncell)
        call redraw_mesh
      ELSEIF (CHOICE.EQ.'Undo') THEN
       if (smflag.eq.0) then
         write(6,*) 'Nothing to undo here'
       else
         write(6,*) 'Undoing mesh smoothing'
         call xtox4(xback,yback,nel,1)
       endif
       call redraw_mesh
      ENDIF
c
      GOTO 1
      end
c-----------------------------------------------------------------------
      subroutine find_sm_box(b,nv,x0,x1,y0,y1,z0,z1)
c
#     include "basics.inc"
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
#     include "basics.inc"
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
               xm3(i,1,1) = x(j,ie)
               ym3(i,1,1) = y(j,ie)
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
               xm3(i,1,1) = x(j,ie)
               ym3(i,1,1) = y(j,ie)
               zm3(i,1,1) = z(j,ie)
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
         IF (SIGN*JAC(I).LE.0.0) return
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
      call jfill(g,n,n)
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
#     include "basics.inc"
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
      return
  999 continue
      call prs('Could not open file.  Returning.$')
      return
      end
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
         return
      ENDIF
*
*     Quick return if possible
*
      IF (N.LE.1) return
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
      return
*
*     End of DLASRT
*
      end
c-----------------------------------------------------------------------
      subroutine clipperN(xyzc,xyzn,DIR,xcl,ycl,zcl,del_list)
#     include "basics.inc"
      real xyzc(3),xyzn(3)
      real xcl(8,1),ycl(8,1),zcl(8,1)
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
            dxyz(1) = xcl(i,je)-xyzc(1)
            dxyz(2) = ycl(i,je)-xyzc(2)
            dxyz(3) = zcl(i,je)-xyzc(3)
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
      return
      end
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
#     include "basics.inc"
c
      logical if_sph_str
C
      call prs(' Input expansion factor ( =< 0 to abort):$')
      call rer(sfact)
      if (sfact.le.0) return

      if_sph_str = .false.
      if (if3d) then
         if_sph_str = .true.
         call prs(' Spherical or Cylindrical (x-y) stretch? (s,c)$')
         call res(ans,1)
         if (ans.eq.'c' .or. ans.eq.'C') if_sph_str = .false.
      endif

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
      rr = 0.
      re = 0.

      do 100 ie=1,nel
      do 100 i=1,nvts
         if (if_sph_str) then
            rt = x(i,ie)**2 + y(i,ie)**2 + z(i,ie)**2
            if (rt.gt.rr) then
               x(i,ie) = sfact*x(i,ie)
               y(i,ie) = sfact*y(i,ie)
               z(i,ie) = sfact*z(i,ie)
            endif
         else
            rt = x(i,ie)**2 + y(i,ie)**2
            if (rt.gt.rr) then
               x(i,ie) = sfact*x(i,ie)
               y(i,ie) = sfact*y(i,ie)
            endif
         endif
  100 continue

C     Take care of curved sides
C
      do 200 ie = 1,nel
      do 200 is = 1,12
         if (ccurve(is,ie).eq.'s') then
            do i=1,4
               curve(i,is,ie) = sfact*curve(i,is,ie)
            enddo
         endif
         if (ccurve(is,ie).eq.'C') then
C           cylindrical side, rad = curve(1,is,ie)
            if (abs(curve(1,is,ie)).gt.re)
     $         curve(1,is,ie) = sfact*curve(1,is,ie)
         endif
         if (ccurve(is,ie).eq.'m') then
            curve(1,is,ie) = sfact*curve(1,is,ie)
            curve(2,is,ie) = sfact*curve(2,is,ie)
            if (if_sph_str) curve(3,is,ie) = sfact*curve(3,is,ie)
         endif

  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine stretch_theta
#     include "basics.inc"
c
      logical if_sph_str
C
      call prs(' Input expansion factor ( =< 0 to abort):$')
      call rer(sfact)
      if (sfact.le.0) return

C
C     Take care of pts first
C
      nvts=8

      do 100 ie=1,nel
      do 100 i=1,nvts

         rad_o   = y(i,ie)**2 + x(i,ie)**2
         if (rad_o.gt.0) rad_o = sqrt(rad_o)
         theta_o = atan2(y(i,ie),x(i,ie))
         theta_n = sfact*theta_o
         x(i,ie) = rad_o*cos(theta_n)
         y(i,ie) = rad_o*sin(theta_n)

  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine clean_spheres
#     include "basics.inc"
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
#     include "basics.inc"
c
      call prs
     $ ('Rep. (1), Rot. (2), Rep/Rot (3) Template (4)? (0=abort)$')
      call res(ans,1)

      if (ans.eq.'1') call replicate_mesh
      if (ans.eq.'2'.and.ndim.eq.2) call rotate_mesh_2d
      if (ans.eq.'2'.and.ndim.eq.3) call rotate_mesh_3d
      if (ans.eq.'3') call rep_rotate_mesh
      if (ans.eq.'4') call template_mesh
c      if (ans.eq.'4') call template_mesh2

      return
      end
c-----------------------------------------------------------------------
      subroutine replicate_mesh
#     include "basics.inc"
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
#     include "basics.inc"
      character*1 axisr

      call prs('Input rotation angle (deg):$')
      call rer(angle_deg)

      call prs('Input number of reps (e.g., 1 --> double mesh size)$')
      call rei(nrep)

      axisr = 'Z'
      if (if3d) then
         call prs('Input axis of rotation (x,y,z):$')
         call res(axisr,1)
         call capit(axisr,1)
      endif

      ie0 = 1
      ie1 = nel
      ie2 = ie1+1
      ie3 = ie2 + nel-1
      do i=1,nrep

         ie2 = ie1+1
         ie3 = ie2 + nel-1

         call copy_sub_mesh(1,nel,ie2)
         angle_deg_i = i*angle_deg
c        call rotate_submesh_2d(ie2,ie3,angle_deg_i)
         call rotate_submesh_3d(ie2,ie3,angle_deg_i,axisr)

         ie0 = ie0 + nel
         ie1 = ie0 + nel-1

      enddo

      nel = ie3

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_sub_mesh(ie0,ie1,inew)
#     include "basics.inc"

c     Map (ie0:ie1) to (inew:inew+nie)

      je = inew
      do ie=ie0,ie1
         call copyel_p(ie,je) ! ie --> je; preserve 'P  ' bc w/ shift
         je = je+1
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine translate_sub_mesh(e0,e1,xt,yt,zt)
#     include "basics.inc"
      integer e,e0,e1,f
c
      do e=e0,e1
c
         do i =1,8
            x(i,e) = x(i,e) + xt
            y(i,e) = y(i,e) + yt
            z(i,e) = z(i,e) + zt
         enddo
c
         do f=1,12
            if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
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
      subroutine get_t1(t,n)
      real t(3),n(3)

      call copy(t,n,3)

      tmin = glmin(t,3)
      tmax = glmax(t,3)
      tnrm = glamax(t,3)

      eps = 1.e-4

      if (tnrm.lt.eps) then
         t(1) = 1.
         call orthogonalize(t,n,3)
         call normalize (t,alpha,3)
c        call outmat(n,1,3,'nveca',1)
c        call outmat(t,1,3,'tveca',1)
         return
      endif

      do i=2,3
         if (abs(t(i)-t(1))/tnrm.gt.eps) then
            t1=t(1)
            t(1) = t(i)
            t(i) = t1
            call orthogonalize(t,n,3)
            call normalize (t,alpha,3)
c           call outmat(n,1,3,'nvecb',1)
c           call outmat(t,1,3,'tvecb',1)
            return
         endif
      enddo

c     If we get here, all 3 components are same and nonzero

      t(1) = 0
      call orthogonalize(t,n,3)
      call normalize (t,alpha,3)
c     call outmat(n,1,3,'nvecc',1)
c     call outmat(t,1,3,'tvecc',1)

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_rotate_mat_3d(a,normal,angle)

      real a(3,3),normal(3)

      real r(3,3),o(3,3),rt(3,3)

c     call outmat(normal,1,3,'norml',1)

      call copy      (r,normal,3)  ! Column 1 of R is normal
      call normalize (r,alpha,3)

c     call outmat    (r,1,3,'rnrml',1)

      call get_t1       (r(1,2),r)    ! arbritrary t1 from normal
      call vcross_normal(r(1,3),sine,r,r(1,2))!t2=nxt1-->n=t1xt2
      call transpose_r (rt,3,r,3)

      one   = 1.
      pi    = 4.*atan(one)
      theta = pi*angle/180.
      c     = cos(theta)
      s     = sin(theta)

      call rzero(o,9)
      o(1,1) = 1
      o(2,2) = c
      o(3,3) = c
      o(2,3) = s
      o(3,2) = -s

c     call outmat(r ,3,3,'r rot',1)
c     call outmat(rt,3,3,'rtrot',2)
c     call outmat(o ,3,3,'omega',3)

      call mxm(o,3,rt,3,a,3)
      call copy(o,a,9)
      call mxm(r,3,o,3,a,3)

c     call outmat(a ,3,3,'a  rot',1)

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_el_vec(e,normal,angle)

c     Rotate element e about axis normal emanating from origin

#     include "basics.inc"
      real normal(3)
      integer e
      real a(3,3)

      call gen_rotate_mat_3d(a,normal,angle)
      call rotate_el_3d(e,a)

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_el_3d(e,a)
#     include "basics.inc"
      real a(3,3)
      integer e,f

      do i=1,2**ndim
         xt = x(i,e)
         yt = y(i,e)
         zt = z(i,e)
         x(i,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
         y(i,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
         z(i,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
      enddo

      do i=1,27
         xt = x27(i,e)
         yt = y27(i,e)
         zt = z27(i,e)
         x27(i,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
         y27(i,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
         z27(i,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
      enddo

      do f=1,12
         if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
            xt = curve(1,f,e)
            yt = curve(2,f,e)
            zt = curve(3,f,e)
            curve(1,f,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
            curve(2,f,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
            curve(3,f,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_mesh_3d_mat
#     include "basics.inc"
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
           xt = x(i,e)
           yt = y(i,e)
           zt = z(i,e)
           x(i,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
           y(i,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
           z(i,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
        enddo

        do i=1,27
           xt = x27(i,e)
           yt = y27(i,e)
           zt = z27(i,e)
           x27(i,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
           y27(i,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
           z27(i,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
        enddo

      enddo
c
      do e=1,nel
      do f=1,12
         if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
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
      subroutine rotate_mesh_3d
#     include "basics.inc"
      real normal(3)
      integer e,f

      call prs('Enter vector, (e.g., 1,1,1):$')
      call rerrr(normal(1),normal(2),normal(3))

      call prs('Enter angle:$')
      call rer(angle)

      do e=1,nel
         call rotate_el_vec(e,normal,angle)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_mesh_2d
#     include "basics.inc"
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
         xt = x(i,e)
         yt = y(i,e)
         x(i,e) = ca*xt-sa*yt
         y(i,e) = sa*xt+ca*yt
      enddo
      enddo

      call rzero(a,9)
      a(1,1) =  ca
      a(1,2) = -sa
      a(2,1) =  sa
      a(2,2) =  ca
      a(3,3) =  1.

      do e=1,nel
      do f=1,12
         if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
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
#     include "basics.inc"
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
         xt = x(i,e)
         yt = y(i,e)
         x(i,e) = ca*xt-sa*yt
         y(i,e) = sa*xt+ca*yt
      enddo
      enddo

      do e=e0,e1
      do f=1,12
         if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
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
      subroutine rotate_submesh_3d(e0,e1,angle_deg,axisr)
#     include "basics.inc"
      character*1 axisr
      integer e,f,e0,e1
      real a(3,3)

      angle = pi*angle_deg/180.
      ca    = cos(angle)
      sa    = sin(angle)

      call rzero(a,9)
      if (axisr.eq.'Z') then
         a(1,1) =  ca
         a(1,2) = -sa
         a(2,1) =  sa
         a(2,2) =  ca
         a(3,3) =  1.
      elseif (axisr.eq.'X') then
         a(1,1) =  1.
         a(3,2) = -sa
         a(2,3) =  sa
         a(2,2) =  ca
         a(3,3) =  ca
      else
         a(1,1) =  ca
         a(1,3) = -sa
         a(3,1) =  sa
         a(3,3) =  ca
         a(2,2) =  1.
      endif

      do e=e0,e1
      do i=1,2**ndim
         xt = x(i,e)
         yt = y(i,e)
         zt = z(i,e)
         x(i,e) = a(1,1)*xt+a(1,2)*yt+a(1,3)*zt
         y(i,e) = a(2,1)*xt+a(2,2)*yt+a(2,3)*zt
         z(i,e) = a(3,1)*xt+a(3,2)*yt+a(3,3)*zt
      enddo
      enddo

      do e=e0,e1
      do f=1,12
         if (ccurve(f,e).eq.'s'.or.ccurve(f,e).eq.'m') then
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
      subroutine template_to_forms

#     include "basics.inc"
      integer e,ef,et,e0

c     Form each element within Template into each element within a 
c     2D or 3D block comprising elements in the subset "Form"
c
c       . Template is stored in etemplate0:etemplate1
c       . Form is stored in     eform0:eform1
c
c     Typically, these are stored at the end of the allocated arrays, 
c     so the upper bound on number of elements is
c
c     nel < nelm - (eform1-eform0) - (etemplate1-etemplate0)
c
c
c     Process:  nel_tmpl = number of elements in template
c               nel_form = number of elements in form
c
c     Result:   nel_form x nel_tmpl elements

      ifmid=.true.

      e0=nel
      do et=etemplate0,etemplate1
         e0 = e0+1
         e  = e0 ! mold a copy of et into ef counterpart, put in e
         do ef=eform0,eform1 
               call template_shape(e,et,ef) 
               e=e+nel_templ
         enddo
      enddo

      nel = nel + nel_templ*nel_form

      return
      end
c-----------------------------------------------------------------------
      subroutine template_shape(e,et,ef) ! Transform et to ef; Store in e
#     include "basics.inc"
      integer e,ef,et

      real jr(27*3),js(27*3),jt(27*3)
      save jr,js,jt
      integer etl
      save etl
      data etl /0/

      call copyel_p(et,e) ! preserve periodic bcs, with (e-et) shift

c     We assume that standard curved boundaries are lost through
c     the transformation and that only midside node definitions are
c     preserved.
c
c     The template and form, however, may be defined by their original
c     curvature information - but only a 27 x 27 transformation matrix
c     is used.
c
      call blank    (ccurve(1,e) ,12)
      call rzero    (curve(1,1,e),72)
c
c     Because The transformation matrix is expensive to generate, we
c     precompute it on the first call.
      if (et.ne.etl) call gen_jrst(jr,js,jt,et)
      etl = et

      if (if3d) then
         call apply_form_3d(x27(1,e),x27(1,ef),jr,js,jt) ! generate 27-noded bricks
         call apply_form_3d(y27(1,e),y27(1,ef),jr,js,jt)
         call apply_form_3d(z27(1,e),z27(1,ef),jr,js,jt)
      else
         call apply_form_2d(x27(1,e),x27(1,ef),jr,js)
         call apply_form_2d(y27(1,e),y27(1,ef),jr,js)
      endif

      call q_to_neklin  (x(1,e),1,x27(1,e),if3d)
      call q_to_neklin  (y(1,e),1,y27(1,e),if3d)
      call q_to_neklin  (z(1,e),1,z27(1,e),if3d)

      call fix_m_curve(e)

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_form_3d(xo,xi,jr,js,jt) ! generate 27-noded bricks

      real xo(27),xi(3,3,3),jr(27,3),js(27,3),jt(27,3)

      do l=1,27      ! this could probably be slightly improved...but..
         xo(l) = 0
         do k=1,3
         do j=1,3
         do i=1,3
           xo(l) = xo(l) + jr(l,i)*js(l,j)*jt(l,k)*xi(i,j,k)
         enddo
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine apply_form_2d(xo,xi,jr,js) ! generate 9-noded bricks

      real xo(9),xi(3,3),jr(9,3),js(9,3)

      do l=1,9      ! this could probably be slightly improved...but..
         xo(l) = 0
         do j=1,3
         do i=1,3
           xo(l) = xo(l) + jr(l,i)*js(l,j)*xi(i,j)
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_jrst(jr,js,jt,et)
#     include "basics.inc"

      real jr(27*3),js(27*3),jt(27*3)
      integer et

      common /ctmp0/ rt(27),st(27),tt(27),rr(3),wk(27*3)

      if (if3d) then
         do i=1,27
            rt(i) = 2*(x27(i,et)-xtmpl_min)/(xtmpl_max-xtmpl_min) - 1.
            st(i) = 2*(y27(i,et)-ytmpl_min)/(ytmpl_max-ytmpl_min) - 1.
            tt(i) = 2*(z27(i,et)-ztmpl_min)/(ztmpl_max-ztmpl_min) - 1.
         enddo
      else
         do i=1,9
            rt(i) = 2*(x27(i,et)-xtmpl_min)/(xtmpl_max-xtmpl_min) - 1.
            st(i) = 2*(y27(i,et)-ytmpl_min)/(ytmpl_max-ytmpl_min) - 1.
         enddo
      endif

      rr(1) = -1
      rr(2) =  0
      rr(3) =  1

      n = 3**ndim
      call gen_int_gz(jr,wk,rt,n,rr,3)
      call gen_int_gz(js,wk,st,n,rr,3)
      call gen_int_gz(jt,wk,tt,n,rr,3)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_template
#     include "basics.inc"

      nelo = nel
      call imp_mesh(.false.)     ! Get template as imported mesh

      nel_templ  = nel-nelo
      etemplate1 = nelm-3
      etemplate0 = etemplate1  - (nel_templ-1)
      call copy_sub_mesh(nelo+1,nel,etemplate0)

c     write(6 ,*) etemplate0,etemplate1,nel_templ,nelo,nel,' etm'

      nel = nelo
      
      xtmpl_min = glmin(x27(1,etemplate0),27*nel_templ)
      xtmpl_max = glmax(x27(1,etemplate0),27*nel_templ)
      ytmpl_min = glmin(y27(1,etemplate0),27*nel_templ)
      ytmpl_max = glmax(y27(1,etemplate0),27*nel_templ)
      ztmpl_min = glmin(z27(1,etemplate0),27*nel_templ)
      ztmpl_max = glmax(z27(1,etemplate0),27*nel_templ)

      return
      end
c-----------------------------------------------------------------------
      subroutine template_mesh
#     include "basics.inc"
c
c     Form each element within Template into each element within a 2D/3D 
c     block comprising elements in the subset "Form"
c
c       . Template is stored in etemplate0:etemplate1
c       . Form is stored in     eform0:eform1
c
c     Typically, these are stored at the end of the allocated arrays, 
c     so the upper bound on number of elements is
c
c     nel < nelm - (eform1-eform0) - (etemplate1-etemplate0)
c
c
c     Process:  nel_tmpl = number of elements in template
c               nel_form = number of elements in form
c
c     Result:   nel_form x nel_tmpl elements


      call get_template


c     Move forms (for now, assumed to be whole existant mesh)

      iform0 = 1
      iform1 = nel

      nel_form = 1 + iform1-iform0 
      eform1   = etemplate0-1
      eform0   = eform1 - (nel_form-1)
      call copy_sub_mesh(iform0,iform1,eform0)

      nel = 0                 ! Nothing left after form-shift

      call template_to_forms  ! Map template to forms

      return
      end
c-----------------------------------------------------------------------
      function blend_circ_in_box2(x,y,r0,x0i,x1i,y0i,y1i)
      real p(2),v(2,5)
      real o(2)
      save o
      data o / 0. , 0. /

      blend_circ_in_box2 = 1.

      r=x*x+y*y
      if (r.le.r0*r0) return

      ey = 1.e-6*(y1i-y0i)   ! Put a slight tolerance on box size
      ex = 1.e-6*(x1i-x0i)
      x0 = x0i-ex
      x1 = x1i+ex
      y0 = y0i-ey
      y1 = y1i+ey

      r = sqrt(r)
      t = atan2(y,x)

      p(1)=x
      p(2)=y

      do i=1,5
         v(1,i)=x0
         v(2,i)=y0
      enddo
      v(1,2)=x1
      v(1,3)=x1
      v(2,3)=y1
      v(2,4)=y1
      call cmult(v,2.0,10) ! Make domain bigger for "in_triangle" check

      ks=0
      do k=1,4
         in = in_triangle2(p,o,v(1,k),v(1,k+1)) ! triangle: [ o vk vk1 ]
         if (in.gt.0) then
           if (k.eq.1) then       !  Lower y boundary
              yc = r0*sin(t)
              b  = (y-y0)/(yc-y0)
           elseif (k.eq.2) then   ! Right x boundary
              xc = r0*cos(t)
              b  = (x-x1)/(xc-x1)
           elseif (k.eq.3) then   ! Upper y boundary
              yc = r0*sin(t)
              b  = (y-y1)/(yc-y1)
           else                   !  Left x boundary
              xc = r0*cos(t) 
              b  = (x-x0)/(xc-x0) 
           endif
           blend_circ_in_box2 = max(b,0.)
           return
         endif
      enddo
      close(78)
      close(83)
      t = sqrt(-t)
      t = sqrt(-t)

      return
      end
c-----------------------------------------------------------------------
      function in_triangle2(p,a,b,c)
      real p(2),a(2),b(2),c(2)
      real v0(2),v1(2),v2(2),invdenom

      call sub3(v0,c,a,2) ! v0 = c-a
      call sub3(v1,b,a,2) ! v1 = b-a
      call sub3(v2,p,a,2) ! v2 = p-a

      d00 = vlsc2(v0,v0,2)
      d01 = vlsc2(v0,v1,2)
      d02 = vlsc2(v0,v2,2)
      d11 = vlsc2(v1,v1,2)
      d12 = vlsc2(v1,v2,2)

c     Compte barycentric coordinates

      invdenom = 1./(d00*d11-d01*d01)
      u=(d11*d02-d01*d12)*invdenom
      v=(d00*d12-d01*d02)*invdenom

      in_triangle2 = 0
      if (u.ge.0.and.v.ge.0.and.(u+v).lt.1) in_triangle2 = 1
      write(83,1) in_triangle2,p(1),p(2),a(1),a(2),b(1),b(2),c(1),c(2)
    1 format(i3,1p8e11.3)

      return
      end
c-----------------------------------------------------------------------
      subroutine stretch_outside_circ
#     include "basics.inc"

      call gencen  ! Generate x27

      n = 27*nel
      xmn=glmin(x27,n)
      xmx=glmax(x27,n)
      ymn=glmin(y27,n)
      ymx=glmax(y27,n)
      zmn=glmin(z27,n)
      zmx=glmax(z27,n)

      call prs('Input protected radius:$')
      call rer(r0)

      call prsrr('Input new xmin/xmax (0,0 to scale): $',xmn,xmx)
      call rerr(x0,x1)
      if (x0.eq.0..and.x1.eq.0.) then
         call prs('Input x-scale factor:$')
         call rer(scale)
         x0 = scale*xmn
         x1 = scale*xmx
      endif


      call prsrr('Input new ymin/ymax (0,0 to scale): $',ymn,ymx)
      call rerr(y0,y1)
      if (y0.eq.0..and.y1.eq.0.) then
         call prs('Input y-scale factor:$')
         call rer(scale)
         y0 = scale*ymn
         y1 = scale*ymx
      endif

      call stretch_outside_circ2(r0,x0,x1,y0,y1,xmn,xmx,ymn,ymx)

      return
      end
c-----------------------------------------------------------------------
      subroutine stretch_outside_circ2(r0,x0,x1,y0,y1,xmn,xmx,ymn,ymx)
#     include "basics.inc"

      n = 27*nel

      scalex = (x1-x0)/(xmx-xmn)
      scaley = (y1-y0)/(ymx-ymn)

      xbr  = (xmn+xmx)/2  ! Midpoint of current data
      ybr  = (ymn+ymx)/2
      xbrt = (x0+x1)/2    ! Target midpoint
      ybrt = (y0+y1)/2

      do i=1,n
         xx = x27(i,1)
         yy = y27(i,1)

         xnt= xbrt + scalex*(xx-xbr) ! Temporary new point
         ynt= ybrt + scaley*(yy-ybr)

         bb = blend_circ_in_box(xx,yy,r0,xmn,xmx,ymn,ymx)

         dx = (xnt-xx)*(1-bb)
         dy = (ynt-yy)*(1-bb)

         xn = xx+dx   ! Now in a proper rectangle to combine
         yn = yy+dy   ! two arc-segment transformations

         x27(i,1) = xn
         y27(i,1) = yn

      enddo

      call x27_to_e  ! Converts all edges to 'm'
      call flatten_edges_outside_circle(r0)

      call redraw_mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine flatten_edges_outside_circle(r0)
#     include "basics.inc"
      integer e

      nedge = 4 + 8*(ndim-2)

      re = 1.01*r0  ! Allow for a bit of wobble
      r2 = re*re

      do e=1,nel
      do i=1,nedge
         if (ccurve(i,e).eq.'m') then
            xx=curve(1,i,e)
            yy=curve(2,i,e)
            rr = xx*xx+yy*yy
            if (rr.gt.r2) ccurve(i,e) = ' '
            if (rr.gt.r2) call rzero(curve(1,i,e),6)
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gencen2(nel1,nel2)
#     include "basics.inc"
c     added by k10 here to do gencen for specific element range

      integer e,nel1,nel2

      common /ctmp2/ xp(nxm,nxm,nxm),yp(nxm,nxm,nxm),zp(nxm,nxm,nxm)

      if (mod(nxm,2).eq.0) then
         write(6,*) 'ERROR: Recompile with nxm odd in basics.inc'
         call prexit(0)
      endif

      nxh = (nxm+1)/2
      nh1 = nxh-1

      write(6,*) 'inside gencen ',nel,if3d,nxh

C     Generate the element centers
      do e=nel1,nel2
         call genxyz_e (xp,yp,zp,e,nxm,nxm,nxm)

         if (if3d) then
            xcen(e)=xp(nxh,nxh,nxh)
            ycen(e)=yp(nxh,nxh,nxh)
            zcen(e)=zp(nxh,nxh,nxh)
            l=0
            do k=1,nxm,nh1
            do j=1,nxm,nh1
            do i=1,nxm,nh1
               l=l+1
               x27(l,e) = xp(i,j,k)
               y27(l,e) = yp(i,j,k)
               z27(l,e) = zp(i,j,k)
            enddo
            enddo
            enddo
         else
            xcen(e)=xp(nxh,nxh,1)
            ycen(e)=yp(nxh,nxh,1)
            zcen(e)=0
            l=0
            do j=1,nxm,nh1
            do i=1,nxm,nh1
               l=l+1
               x27(l,e) = xp(i,j,1)
               y27(l,e) = yp(i,j,1)
               z27(l,e) = 0
            enddo
            enddo
         endif

      enddo

C     Compute the maximum radius from the center

      call rzero(rcen,nel2)
      if (if3d) then
         do 300 e=nel1,nel2
         do 300 j=1,8
            rad=(x(j,e)-xcen(e))**2 + (y(j,e)-ycen(e))**2
     $         +(z(j,e)-zcen(e))**2
            rcen(e)=max(rcen(e),rad)
  300    continue
      else
         do 400 e=nel1,nel2
         do 400 j=1,4
            rad=(x(j,e)-xcen(e))**2 + (y(j,e)-ycen(e))**2
            rcen(e)=max(rcen(e),rad)
  400    CONTINUE
      ENDIF
      call vsqrt(rcen,nel)

      write(6,*) 'done gencen ',nel2-nel1+1

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_internal_bcs
#     include "basics.inc"
      integer n,e,i,f,j,ef,et
      integer n1,e1,e2,f1,f2

      write(6,*) nel,'neigh being found for these many elements'
      call gen_neigh
      n1 = nel 
      n=0
      write(6,*) 'Trying to fix internal boundary conditions'
      do j=1,1 
      do e=1,nel
      do f=1,6
       if(cbc (f,e,j).ne.'E  ') then
         call get_neigh(e,f,e2,f2)
         if (e2.ne.0) cbc(f,e,j) = 'E  '
       endif
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2col2(a,b,c,n)
      real a(1),b(1),c(1)
c
      do i=1,n
         a(i) = a(i) + b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE INVCOL1 (A,N)
      REAL A(1)
      do I=1,N
         A(I)=1./A(I)
      enddo
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine cheap_dist(d,ifld,b,dims,eln,nvv,ncell,nnpts)
c     added by k10

c     Finds a pseudo-distance function.

c     INPUT:  ifld - field type for which distance function is to be found.
c             ifld = 1 for velocity
c             ifld = 2 for temperature, etc.

c     OUTPUT: d = "path" distance to nearest wall

c     This approach has a significant advantage that it works for
c     periodict boundary conditions, whereas most other approaches
c     will not.

#     include "basics.inc"

      integer e,eg,f,nx1,ny1,nz1,i,j,nnods,nnpts,nchange,nnface,nvv
      integer dims,ncell,ifld,eln(nnpts,ncell),ipass
      real d(2*2*(dims-1),ncell)
      character*3 b  ! Boundary condition of interest
      integer efc3(4,6)
      data    efc3 / 1,2,6,5
     $             , 2,3,7,6
     $             , 3,4,8,7
     $             , 4,1,5,8
     $             , 1,2,3,4
     $             , 5,6,7,8 /

      nx1 = 2
      ny1 = 2
      nz1 = dims-1
      n = nx1*ny1*nz1*ncell
      call domain_sizek(xmin,xmax,ymin,ymax,zmin,zmax,dims,ncell)

      xmn = min(xmin,ymin)
      xmx = max(xmax,ymax)
      if (dims.eq.3) xmn = min(xmn ,zmin)
      if (dims.eq.3) xmx = max(xmx ,zmax)

      big = 10*(xmx-xmn)
      call rone(d,n) 
      call cmult(d,big,n)

      nnface = 2*dims
      nnods = 2+(dims-2)*2 !number of nodes per edge/face

       do e=1,ncell     ! Set d=0 on walls
       do f=1,nnface
          if (cbc(f,e,ifld).eq.b) then
            do j=1,nnods
              d(efc3(j,f),e) = 0.*d(efc3(j,f),e)
            enddo
          endif
       enddo
       enddo

      do ipass=1,10000
         dmax    = 0
         nchange = 0
         do e=1,nel
          do i=1,nnpts
           do j=1,nnpts
             if (dims.eq.3) then
         dtmp = d(j,e)+dist3d(x(j,e),y(j,e),z(j,e),x(i,e),y(i,e),z(i,e))
             else
         dtmp = d(j,e)+dist2d(x(j,e),y(j,e),x(i,e),y(i,e))
             endif
             if (dtmp.lt.d(i,e)) then
               d(i,e) = dtmp
               nchange = nchange+1
               dmax = max(dmax,d(i,e))
              endif
           enddo
          enddo
        enddo


         call qqtmin(eln,nvv,ncell,d,nnpts)
         dmax = glmax(dmax,1)
         write(6,1) ipass,nchange,dmax,b
    1    format(i9,i12,1pe12.4,' max distance b: ',a3)
         if (nchange.eq.0) goto 1000
      enddo
 1000 return
      end
c-----------------------------------------------------------------------
      subroutine domain_sizek(xmin,xmax,ymin,ymax,zmin,zmax,dims,ncell)
#     include "basics.inc"
      integer n,dims,e,i
      real xvv,yvv,xmin,xmax,ymin,ymax,zmin,zmax 

      n = 2*2*(dims-1)*ncell

      if (dims.eq.3) then 
        xmin = glmin(x,n)
        xmax = glmax(x,n)
        ymin = glmin(y,n)
        ymax = glmax(y,n)
        zmin = glmin(zm1,n)
        zmax = glmax(zm1,n)
      else
         xmin = 1.e+7
         xmax = -1.e+7
         ymin = 1.e+7
         ymax = -1.e+7
         do e=1,ncell
         do i=1,4
           xvv = x(i,e)
           yvv = y(i,e)
           if (xvv.gt.xmax) xmax = xvv
           if (yvv.gt.ymax) ymax = yvv
           if (xvv.lt.xmin) xmin = xvv
           if (yvv.lt.ymin) ymin = yvv
         enddo
         enddo
         zmin = 0.
         zmax = 0.
      endif

      return
      end
c-------------------------------------------------------------------------
      function dist3d(a,b,c,x,y,z)

      d = (a-x)**2 + (b-y)**2 + (c-z)**2

      dist3d = 0.
      if (d.gt.0) dist3d = sqrt(d)

      return
      end
c-----------------------------------------------------------------------
      function dist2d(a,b,x,y)

      d = (a-x)**2 + (b-y)**2

      dist2d = 0.
      if (d.gt.0) dist2d = sqrt(d)

      return
      end
c-----------------------------------------------------------------------
      subroutine qqt(eln,nv,ncell,qloc,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)

      call loctoglob(eln,nv,ncell,qloc,qglob,nnpts)
      call globtoloc(eln,nv,ncell,qloc,qglob,nnpts)

      return
      end
c-----------------------------------------------------------------------
      subroutine qqtavg(eln,nv,ncell,qloc,nodcl,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell),nodcl(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)

      call loctoglob(eln,nv,ncell,qloc,qglob,nnpts)
      call globtoloc(eln,nv,ncell,qloc,qglob,nnpts)
      call invcol2(qloc,nodcl,ncell*nnpts)

      return
      end
c-----------------------------------------------------------------------
      subroutine loctoglob(eln,nv,ncell,qloc,qglob,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)
      call rzero(qglob,nv)
      do e=1,ncell
      do i=1,nnpts
            igl = eln(i,e)
            qglob(igl) = qglob(igl)+qloc(i,e)
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine globtoloc(eln,nv,ncell,qloc,qglob,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)

      call rzero(qloc,nnpts*ncell)

      do e=1,ncell
      do i=1,nnpts
            igl = eln(i,e)
            qloc(i,e) = qglob(igl)
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine loctoglobmin(eln,nv,ncell,qloc,qglob,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv),maxv
      maxv = 1.e+11
      call rone(qglob,nv)
      call cmult(qglob,maxv,nv)
      do e=1,ncell
      do i=1,nnpts
            igl = eln(i,e)
            qglob(igl) = min(qglob(igl),qloc(i,e))
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine loctoglobmax(eln,nv,ncell,qloc,qglob,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv),minv
      minv = -1.e+11
      call rone(qglob,nv)
      call cmult(qglob,minv,nv)

      do e=1,ncell
      do i=1,nnpts
            igl = eln(i,e)
            qglob(igl) = max(qglob(igl),qloc(i,e))
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine qqtmin(eln,nv,ncell,qloc,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)

      call loctoglobmin(eln,nv,ncell,qloc,qglob,nnpts)
      call globtoloc(eln,nv,ncell,qloc,qglob,nnpts)

      return
      end
c-----------------------------------------------------------------------
      subroutine qqtmax(eln,nv,ncell,qloc,nnpts)
#     include "basics.inc"
      integer nv,ncell,i,e,nnpts
      integer eln(nnpts,ncell)
      real qloc(nnpts,ncell)
      real qglob(nv)

      call loctoglobmax(eln,nv,ncell,qloc,qglob,nnpts)
      call globtoloc(eln,nv,ncell,qloc,qglob,nnpts)

      return
      end
c-----------------------------------------------------------------------
      subroutine disfun(d,delta,nnpts,ncell)
#     include "basics.inc"
      integer nnpts,ncell
      real d(ncell*nnpts),dd2(ncell*nnpts)
      real dis(ncell*nnpts),delta
      integer i,n
      real dscale,dmax

      write(6,*) 'distance function being calculated'

      n = nnpts*ncell
      dmax   = glamax(d,n)
      dscale = 1./dmax
      call cmult(d,dscale,n)
 
      do i=1,ncell*nnpts
        dd2(i) = (1-EXP(-d(i)/delta))
c        dd2(i) = 0.5*(tanh(10.*(dis(i)-0.3))+1)
      enddo
      call copy(d,dd2,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine restbndrlay(x4,y4,x4b,y4b,d,ncell)
#     include "basics.inc"
      integer nnpts,ncell
      real x4(4,ncell),y4(4,ncell),x4b(4,ncell),y4b(4,ncell)
      real d(4,ncell)
      integer i,j,n
! x4b is the original mesh
! x4 is the smooth mesh
! d is the weightin function

      n = 4*ncell
      call sub2(x4,x4b,n) !x4 = x4-x4b
      call sub2(y4,y4b,n)

      call col2(x4,d,n)  !x4 = x4*d = d(x4-x4b)
      call col2(y4,d,n)

      call add2(x4,x4b,n) !x4 = d(x4-x4b)+x4b
      call add2(y4,y4b,n) !x4 = d(x4) + (1-d)x4b

      return
      end
c-----------------------------------------------------------------------
      subroutine get_neigh(e1,f1,e2,f2)
#     include "basics.inc"
      integer j,e,e1,f1,e2,f2

      f2 = 0
      e2 = neighb(f1,e1)
      do j=1,6
         if (neighb(j,e2).eq.e1) f2=j
      enddo
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------

