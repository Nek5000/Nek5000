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
c-----------------------------------------------------------------------
      subroutine BOUND
C     Sets boundary conditions
      include 'basics.inc'
      CHARACTER heatbc,velbc,KEY,STRING*5,ELE(4),BCLAB(2),REPLY,QUES(2)
     $,REPLY2,IPAPT,BCCHOICE*26,MODE*15,CMENU*15
      INTEGER FCORNS (4,6),ICALLD
      SAVE    FCORNS,ICALLD
      LOGICAL IFTMP
      DATA CMENU /' '/
      DATA FCORNS / 1,2,6,5,
     $              2,3,7,6,
     $              3,4,8,7,
     $              4,1,5,8,
     $              1,2,3,4,
     $              5,6,7,8 /
      DATA ICALLD / 0 /
C     Turn of grid latch for improved mousy resolution.
      IFGRID=.FALSE.
C       PROMPTS USER FOR B.C.'S
C
C!!?? AFTER READING IN MESH, DISPLAY B.C.'S.  CHECK FOR SUSPECTED ONES.
C     Boundary condition storage
C     CBC array of characters. (6,NELM,MAXFLD) Side (or face); Element; Field
C     BC  array of reals.    (3,6,NELM,MAXFLD) 1st Subscript: 1-temp;2-flux;3-?
C     BC  array of reals.    For Internal Bndry:neighbor side; ele; orientation
1     CONTINUE
      BCLAB(2)='$'
      QUES(1)='?'
      QUES(2)='$'
      IF(.NOT.IF3D) NSIDES=4
      IF(     IF3D) NSIDES=6
      IF1=1
      IF(.NOT. IFFLOW) IF1=2
C! No Fluid B.C.'s
C     Check if b.c.'s are necessary
C     ??!! Also check if old b.c.'s are reasonable?? Put in ability to edit
C     the .rea file to make overlapping sides have branch-cut b.c.'s??
      ELEC= 0.95
      ELEW= 0.02
      ELEB= 0.1
      ELEDY=0.1
      NEEDBC=0
      MAXLEV=0
C     Make sure B.C.'s for nonexistent elements are blanked out.
      ifld0=if1
      ifld1=nflds
      if (ifconj_merge) then
         ifld0=2
         ifld1=2
      endif
      do 28 if=ifld0,ifld1
        IF(NELF.NE.NEL)THEN
           DO 27 IEL=NELF+1,NEL
              DO 27 ISIDE=1,NSIDES
C                 IF(.NOT.IFTMSH(IF))
C     $           PRINT*,'iel,iside,if',iel,iside,if
                 IF(.NOT.IFTMSH(IF))CBC(ISIDE,IEL,IF)=' '
27         CONTINUE
        ENDIF
28    continue

      do 29 if=ifld0,ifld1
        IF(     IFTMSH(IF))NNEL=NEL
        IF(.NOT.IFTMSH(IF))NNEL=NELF
        DO 29 IEL=1,NNEL
          IF(NUMAPT(IEL).GT.MAXLEV)MAXLEV=NUMAPT(IEL)
          DO 29 ISIDE=1,NSIDES
          IF(CBC(ISIDE,IEL,IF).EQ.' ')NEEDBC=1
29    CONTINUE
C     Display B.C.'s already here (for 2-d case)
C
      if(nlevel.eq.1) then
        do 30 if=ifld0,ifld1
        do 30 iel=1,nel
        do 30 iside=1,nsides
           IF(CBC(ISIDE,IEL,IF).NE.' ')
     $         CALL LETBC(ISIDE,IEL,IF,CBC(ISIDE,IEL,IF))
30      continue
      endif
C
c
c
c     Dump element centroids, w/ element numbers for positional sorting
c                                   pff 4/8/99
      call out_cent
c
      IF(NEEDBC.EQ.1)THEN
C         Ya gotta have b.c.'s
         CALL PRS('              *** MIDWAY BREAK MENU ***$')
         CALL PRS('ABORT to write build data to file and exit$')
         CALL PRS('SET B.C''S to continue$')
         ITEM(1)='SET BCs'
         ITEM(2)='ABORT'
         nchoic =  2
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'MIDWAY BREAK')
         IF(CHOICE.EQ.'ABORT')THEN
             CALL MESGEN
             CALL PREXIT
         ENDIF
      ELSE
         ITEM(1)='ACCEPT B.C.''s'
         ITEM(2)='REVIEW/MODIFY'
         nchoic = 2
         if (ifconj_merge) then
            choice = item(1)
         else
            CALL MENU(XMOUSE,YMOUSE,BUTTON,'ACCEPT/REVIEW')
         endif
         IF(CHOICE.EQ.'ACCEPT B.C.''s')THEN
C           Where do you jump if he/she wants to accept??
            CALL CHKBCS
c
            DO IF=IF1,NFLDS
               call period_check(if)
            ENDDO
c
            RETURN
         ENDIF
      ENDIF
      IF(CHOICE.EQ.'SET BCs GLOBALLY'.OR.CHOICE.EQ.'GLOBAL MODIFY')
     $THEN
         ITEM(1)='ALL BOUNDARIES'
         ITEM(2)='ALL FLOORS'
         ITEM(3)='ALL CEILINGS'
         ITEM(4)='ONE LEVEL'
         NCHOIC=4
         CALL MENU(XMOUSE,YMOUSE,BUTTON,'TITLE$')
         IF(CHOICE.EQ.'ALL BOUNDARIES')THEN


            CALL PRS('ENTER B.C.TYPE$')
            CALL RES(ANS,1)
            IF(IFLEARN)WRITE(3,'(A1)') ANS
            IF(ANS.EQ.'W'.OR.ANS.EQ.'w')THEN
               DO 39 IF=IF1,NFLDS
                 DO 39 IEL=1,NEL
                   DO 39 ISIDE=1,NSIDES
                     IF(CBC(ISIDE,IEL,IF).NE.'E')
     $               CBC(ISIDE,IEL,IF)='W'
39             CONTINUE
            ELSE
               CALL PRS('NOT IMPLEMENTED YET. USE w for wall.$')
            ENDIF
         ELSE
            CALL PRS('NOT IMPLEMENTED YET. USE ALL BOUNDARIES.$')
         ENDIF
         GO TO 1
      ENDIF
C
      CALL PRS('                   *** BOUNDARY CONDITION MENU ***$')
C     Set up gray elevator
      IF(NDIM.EQ.3)THEN
         IF(MAXLEV.GE.2)THEN
            DO 50 I=MAXLEV,2,-1
               CALL DRELEV(I-1,I,'     ')
50          CONTINUE
         ELSE
            CALL DRELEV(I,0,'     ')
         ENDIF
      ENDIF
C
C     On each level First demand all unset B.C.'s; then give option to modify
      DO 2000 IF=IF1,NFLDS
C
         IF(IF.EQ.1) THEN
             CALL PRS('Begin Inputting Fluid Boundary Conditions$')
         ELSE IF(IF.EQ.2) THEN
             CALL PRS('Begin Inputting Thermal Boundary Conditions$')
         ELSE
             CALL PRS
     $       ('Begin Inputting PASSIVE SCALAR Boundary Conditions$')
             NN=IF-2
             CALL PRSI('PASSIVE SCALAR #$',NN)
         ENDIF
C
         DO 1000 ILEVEL=1,NLEVEL
1010        CONTINUE
            IF(NLEVEL.GT.1) THEN
C              Black out whole level
               CALL COLOR(0)
               CALL FILLP(0)
               CALL BEGINP(XPHY(0.0),YPHY(0.0))
               CALL DRAW  (XPHY(.90),YPHY(0.0))
               CALL DRAW  (XPHY(.90),YPHY(.99))
               CALL DRAW  (XPHY(0.0),YPHY(.99))
               CALL DRAW  (XPHY(0.0),YPHY(0.0))
               CALL ENDP
C              Make Isometric drawing
               IF (ICALLD.EQ.0) THEN
                  ICALLD=1
                  DO 122 IEL=1,NEL
                     CALL DRAWIS(ISRT(IEL))
122               CONTINUE
               ENDIF
C              Draw mesh OR WHOLE ELEMENT?  ??Delete previous mesh(+bclabels?)
               CALL COLOR(1)
               NTHISL=0
               DO 300 IEL=1,NEL
                   IF(NUMAPT(IEL).EQ.ILEVEL)THEN
                      NTHISL=NTHISL+1
                      CALL DRAWEL(IEL)
                      DO 200 ISIDE=1,NSIDES
                         IF(CBC(ISIDE,IEL,IF).NE.' ')
     $                   CALL LETBC(ISIDE,IEL,IF,CBC(ISIDE,IEL,IF))
200                   CONTINUE
                   ENDIF
300            CONTINUE
               IF(NTHISL.EQ.0) GO TO 1999
               IF(NDIM.EQ.3)CALL DRELEV(ILEVEL,ILEVEL-1,'     ')
            ENDIF
C           Go thru normal stuff
800         CONTINUE
            MODE='NORMAL'
            DO 900 ISIDE=1,NSIDES
               DO 900 IEL=1,NEL
                  IF(NUMAPT(IEL) .NE. ILEVEL )GO TO 900
                  IF(CBC(ISIDE,IEL,IF).NE.' ')GO TO 900
                  IF(MASKEL(IEL,IF) .EQ. 0   )GO TO 900
C                 Side which is on Boundary and Not periodic with Previous el
C
                  IISIDE=ISIDE+1
                  IF(IISIDE.EQ.5) IISIDE = 1
C                 Flash Elevator
                  IFLASH=ILEVEL
                  IF(ISIDE.EQ.5)IFLASH= -ILEVEL
                  IF(ISIDE.EQ.6)IFLASH=-(ILEVEL+1)
                  IFADE=IFLASH
                  CALL DRELEV(IFLASH,IFADE,'BLINK')
C                 Flash element side
C                 Hilight side with Blinking lights
                  CALL COLOR(2)
                  IF(IFBWGKS)CALL COLOR(0)
                  CALL HILITES(IEL,ISIDE)
                  CALL LETBC(ISIDE,IEL,IF,QUES)
c                 
                  IF(IF3D)THEN
C                    FLASH SIDE ON ISOMETRIC
                     DO 223 ICORN=0,4
                        IF(ICORN.EQ.0) IC=FCORNS(4    ,ISIDE)
                        IF(ICORN.NE.0) IC=FCORNS(ICORN,ISIDE)
                        XISOM=XPHY(XSCR(X(IEL,IC))/5.0 + 0.8)
     $                  -Z(IEL,IC)/20.
                        YISOM=YPHY(YSCR(Y(IEL,IC))/5.0 + 0.8)
     $                  -Z(IEL,IC)/20.
                        IF(ICORN.EQ.0) CALL MOVEC(XISOM,YISOM)
                        IF(ICORN.NE.0) CALL DRAWC(XISOM,YISOM)
223                  CONTINUE
                  ENDIF
C
113               CONTINUE
C                 Error Reading Input
                  FLUX=0.0
                  TEMP=0.0
                  VX  =0.0
                  VY  =0.0
                  IF(MODE.NE.'SKIP MENU')CALL PRSii(' B.C.>$',iel,iside)
                  IDUP=0
C                 SET UP MENU CHOICES
                  IF(IF.EQ.1) THEN
C                    FLUID B.C.'s
                     nchoic=1
                     ITEM(NCHOIC)='PERIODIC-AUTO'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='PERIODIC'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='WALL'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='VELOCITY'
                     nchoic=nchoic+1
                     ITEM(nchoic)='OUTFLOW'
                     nchoic=nchoic+1
                     ITEM(nchoic)='OUTFLOW/N'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='VELOCITY (LOCAL)'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='SYMMETRY'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='MOVING WALL'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='MELTING FRONT'
c                    nchoic=nchoic+1
c                    ITEM(NCHOIC)='QUASI-SYMMETRY'
                     IF(IFAXIS)THEN
                        nchoic=nchoic+1
                        ITEM(NCHOIC)='AXIS'
                     ENDIF
                     IF(IFSTRS)THEN
                        nchoic=nchoic+1
                        ITEM(NCHOIC)='STRESS B.C.s'
                     ENDIF
                  ELSE IF(IF.GE.2)THEN
C                    THERMAL B.C.'s  USE FOR PASSICVE SCALARS
                     nchoic=1
                     ITEM(NCHOIC)='PERIODIC-AUTO'
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='PERIODIC'
                     nchoic=nchoic+1
                     ITEM(nchoic)='INSULATED'
                     nchoic=nchoic+1
                     ITEM(nchoic)='FLUX'
                     nchoic=nchoic+1
                     ITEM(nchoic)='TEMP'
C                    Convection IS now ready for v2
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='CONVECTION'
                     IF(IFFLOW.OR.IFADVC(IF))THEN
                        nchoic=nchoic+1
                        ITEM(NCHOIC)='OUTFLOW'
                     ENDIF
                     IF(IFAXIS)THEN
                        nchoic=nchoic+1
                        ITEM(NCHOIC)='AXIS'
                     ENDIF
                     nchoic=nchoic+1
                     ITEM(NCHOIC)='RADIATION'
                  ENDIF
                  nchoic=nchoic+1
                  ITEM(NCHOIC)='SET ENTIRE LEVEL'
                  nchoic=nchoic+1
                  ITEM(nchoic)='ZOOM'
C
                  IF(CMENU.EQ.'INTERNAL')THEN
C                    Start over with short menu
                     IF(IF.EQ.1)THEN
                        nchoic=1
                        ITEM(NCHOIC)='FLUID LAYERS'
                     ELSE IF(IF.EQ.2)THEN
C                        Fill up solid & liquid front when melting front is
C                        chosen
C                        NCHOIC       =1
C                        ITEM(NCHOIC)='SOLID FRONT'
C                        nchoic=NCHOIC+1
C                        ITEM(NCHOIC)='LIQUID FRONT'
C
c                        nchoic=NCHOIC+1   No internal radiation yet
c                        ITEM(NCHOIC)='RADIATION'
                        nchoic=1
                        ITEM(NCHOIC)='MELTING FRONT'
C                        CALL PRS(
c     $                  'No internal boundaries for temperature$')
c                        GO TO 313
                     ELSE IF(IF.GT.2)THEN
                        CALL PRS(
     $                  'No internal boundaries for passive scalars$')
                        GO TO 313
                     ENDIF
                  ENDIF
C
C                 We already got info for first value; skip menu now.
                  IF(MODE.EQ.'GET FIRST')MODE='SKIP MENU'
63                IF(MODE.NE.'SKIP MENU')THEN
                     IF(IFAXIS .AND.
     $               ABS(SIDES(IEL,ISIDE,2)).LT.YFAC/1000.)THEN
                        CALL PRS
     $                  ('**Axis B.C. Strongly Recommended here**$')
                        CALL BEEP
                     ENDIF
                     CALL MENU (XMOUSE,YMOUSE,BUTTON,'NOCOVER')
                  ENDIF
C
                  IF(CHOICE.eq.'STRESS B.C.s')THEN
C
                        NCHOIC=       1
                        ITEM(NCHOIC)='MAIN B.C. MENU'
                        NCHOIC=NCHOIC+1
                        ITEM(NCHOIC)='SHEAR'
                        NCHOIC=NCHOIC+1
                        ITEM(NCHOIC)='SHEAR    (LOCAL)'
                        NCHOIC=NCHOIC+1
                        ITEM(NCHOIC)='FREE SURFACE'
C                    New B.C.'s From Lee Ho
                        NCHOIC=NCHOIC+1
                        ITEM(NCHOIC)='TRACTION'
                        NCHOIC=NCHOIC+1
                        ITEM(NCHOIC)='TRACTION (LOCAL)'
C                        Catch at internal b.c. menu
c                        NCHOIC=NCHOIC+1
c                        ITEM(NCHOIC)='FLUID LAYERS'
                        CALL MENU (XMOUSE,YMOUSE,BUTTON,'NOCOVER')
                        IF(CHOICE.EQ.'MAIN B.C. MENU') GO TO 113
                  ENDIF
C                 End of patch
C
                 IF(MODE.EQ.'GET FIRST'.AND.CHOICE.EQ.'PERIODIC') then
c                IF((MODE.EQ.'GET FIRST'.AND.CHOICE.EQ.'PERIODIC').OR.
c    $              (MODE.EQ.'GET FIRST'.AND.
c    $               CHOICE.EQ.'PERIODIC-AUTO'))THEN
                     CALL PRS
     $               ('*** ERROR ***  YOU CANNOT SET WHOLE LEVEL$')
                    CALL PRS('TO BE PERIODIC; CHOOSE A DIFFERENT B.C.$')
                     GO TO 63
                  elseif (MODE.EQ.'GET FIRST'.AND.
     $                    CHOICE.EQ.'PERIODIC-AUTO') THEN
                     call autoperiod
                  ENDIF
                  if(choice.eq.'ZOOM')then
                     CALL SETZOOM
                     CALL REFRESH
                     CALL DRMENU('NOCOVER')
                     CALL DRGRID
                     DO 170 IE=1,NEL
                        CALL DRAWEL(IE)
  170                CONTINUE
C
C                    Re-highlight the current element boundaries
C
                     CALL COLOR(2)
                     IF(ISIDE.EQ.5.OR.ISIDE.EQ.6) THEN
                       IF(ISIDE.EQ.5)ILEV=ILEVEL
                       IF(ISIDE.EQ.6)ILEV=ILEVEL+1
                       CALL MOVESC(ELEC-ELEW,ELEB+ELEDY*(ILEV-1))
                       CALL DRAWSC(ELEC+ELEW,ELEB+ELEDY*(ILEV-1))
C                      Flash mesh sides
                       IF(ISIDE.EQ.5)THEN
                          IED1=1
                          IED2=4
                       ELSE IF(ISIDE.EQ.6)THEN
                          IED1=1
                          IED2=4
                       ENDIF
                       CALL MOVEC(X(IEL, IED1),Y(IEL, IED1))
                       DO 172 IEDGE=IED1,IED2
                          CALL DRAWED(IEL,IEDGE,1)
  172                  CONTINUE
                     ELSE
C                       SIDES 1-4  (SAME AS EDGES IN THIS CASE)
                        CALL MOVEC(X(IEL, ISIDE),Y(IEL, ISIDE))
                        CALL DRAWED(IEL,ISIDE,1)
                     ENDIF
                     GO TO 63
                  ENDIF
                  if(choice.eq.'SET ENTIRE LEVEL')then
C                    First get b.c., then copy for rest of loops
                     MODE='GET FIRST'
                     CHOICE=' '
                     IELDUP=IEL
                     ISIDUP=ISIDE
                     CALL PRS
     $               ('Choose B.C.''s for all remaining sides '//
     $               'in entire level$')
                     GO TO 63
                  ENDIF
                  IF(MODE.EQ.'SKIP MENU')BUTTON='RIGHT'
                  IF(BUTTON.EQ.'RIGHT') THEN
                     BCCHOICE='DUPLICATE'
                     IF(XSCR(XMOUSE).GE.1.0.AND.MODE.NE.'SKIP MENU')THEN
                        CALL PRS
     $                  ('Right button not used in red menu area.$')
                        CALL PRS('Try again.$')
                        IF(MODE.NE.'NORMAL')CALL PRS(
     $                  'Global mode reset.  Start over.$')
                        MODE='NORMAL'
                        GO TO 113
                     ENDIF
                     IF(MODE.NE.'SKIP MENU')THEN
                       XP=XMOUSE
                       YP=YMOUSE
C                      Find closest element, side that we want to duplicate
                       RMIN=10000.0
                       DO 133 IIEL=1,NEL
C                         Only try those on same floor
                          IF(NUMAPT(IIEL).EQ.ILEVEL)THEN
                             DO 131 IISIDE=1,NSIDES
                                R=    ((XP-SIDES(IIEL,IISIDE,1))**2
     $                          +      (YP-SIDES(IIEL,IISIDE,2))**2)
                                IF(R.LT.RMIN) THEN
                                   RMIN=R
                                   ISIDUP = IISIDE
                                   IELDUP = IIEL
                                ENDIF
131                          CONTINUE
                          ENDIF
133                    CONTINUE
                       IF(ISIDUP.EQ.5.OR.ISIDUP.EQ.6)THEN
C                         CHECK IF IT'S ABOVE OR BELOW CENTER
                          IF(YP.GT.SIDES(IELDUP,ISIDUP,2))THEN
                              ISIDUP=6
                          ELSE
                              ISIDUP=5
                          ENDIF
                       ENDIF
                     ENDIF
                     IF (CBC(ISIDUP,IELDUP,IF).EQ.' ') THEN
                        call prsi (
     $                      'B.C. not set for isidup$',isidup,
     $                      'of element$',ieldup)
                        go to 113
                     ENDIF
                     CBC3=CBC(ISIDUP,IELDUP,IF)
                     IF (CBC(ISIDUP,IELDUP,IF).EQ.'E') THEN
                        WRITE(S,'(A5,I4,A12,I4,A1)')'Side ',isidup,
     $                  ' of element ',IELDUP,'$'
                        CALL PRS(S)
                        CALL PRS
     $                  ('is internal side.  I won''t duplicate.$')
                        go to 113
                     ELSE IF (CBC(ISIDUP,IELDUP,IF).EQ.'P  ') THEN
                        WRITE(S,'(A5,I4,A12,I4,A1)')'Side ',isidup,
     $                  ' of element ',IELDUP,'$'
                        CALL PRS(S)
                        CALL PRS
     $                  ('is Periodic side.  I won''t duplicate.$')
                        go to 113
                     ELSE IF(CBC3(3:3).EQ.'i'.OR.CBC3(3:3).EQ.'I')THEN
                        WRITE(S,'(A5,I4,A12,I4,A1)')'Side ',isidup,
     $                  ' of element ',IELDUP,'$'
                        CALL PRS(S)
                        CALL PRS
     $                  ('is Internal side.  I won''t duplicate.$')
                        go to 113
                     ENDIF
C                    Need flag to know if b.c. is fortran function.
C                    Special handling puts new fortran in next available lines
                     CBC3 = CBC(ISIDUP,IELDUP,IF)
                     ICBC = ICHAR(CBC3(1:1))
                     IF(ICBC.GE.97.AND.ICBC.LE.122)THEN
C                       Small letter signifies fortran function
                        IF(CBC3(3:3).eq.'i')THEN
C                          Special storage locations
                           NLINES=BC(4,ISIDUP,IELDUP,IF)
                           LINE1 =BC(5,ISIDUP,IELDUP,IF)
                           BC(4,ISIDE,IEL,IF)=NLINES
                           BC(5,ISIDE,IEL,IF)=LOCLIN
                        ELSE
                           NLINES=BC(1,ISIDUP,IELDUP,IF)
                           LINE1 =BC(2,ISIDUP,IELDUP,IF)
                           BC(1,ISIDE,IEL,IF)=NLINES
                           BC(2,ISIDE,IEL,IF)=LOCLIN
                        ENDIF
                        CBC (ISIDE,IEL,IF)=CBC( ISIDUP,IELDUP,IF)
                        DO 777 I=1,NLINES
                           INBC(LOCLIN)=INBC(LINE1+I-1)
                           LOCLIN=LOCLIN+1
777                     CONTINUE
                     ELSE
                        BC(1,ISIDE,IEL,IF)=BC(1,ISIDUP,IELDUP,IF)
                        BC(2,ISIDE,IEL,IF)=BC(2,ISIDUP,IELDUP,IF)
                        BC(3,ISIDE,IEL,IF)=BC(3,ISIDUP,IELDUP,IF)
                        BC(4,ISIDE,IEL,IF)=BC(4,ISIDUP,IELDUP,IF)
                        BC(5,ISIDE,IEL,IF)=BC(5,ISIDUP,IELDUP,IF)
                        CBC (ISIDE,IEL,IF)=CBC( ISIDUP,IELDUP,IF)
                     ENDIF
                     BCLAB(1) = CBC(ISIDUP,IELDUP,IF)
                  ELSE IF(BUTTON.EQ.'LEFT') THEN
                     BCCHOICE=CHOICE
                     BCLAB(1)          = BCCHOICE(1:1)
                     CBC(ISIDE,IEL,IF) = BCCHOICE(1:1)
           IF(CHOICE.EQ.'VELOCITY'        )CBC(ISIDE,IEL,IF)='V'
           IF(CHOICE.EQ.'VELOCITY (LOCAL)')CBC(ISIDE,IEL,IF)='VL'
           IF(CHOICE.EQ.'MOVING WALL     ')CBC(ISIDE,IEL,IF)='mv'
           IF(CHOICE.EQ.'STRESS'          )CBC(ISIDE,IEL,IF)='S'
           IF(CHOICE.EQ.'STRESS   (LOCAL)')CBC(ISIDE,IEL,IF)='SL'
           IF(CHOICE.EQ.'SHEAR'           )CBC(ISIDE,IEL,IF)='SH'
           IF(CHOICE.EQ.'SHEAR    (LOCAL)')CBC(ISIDE,IEL,IF)='SHL'
           IF(CHOICE.EQ.'SYMMETRY'        )CBC(ISIDE,IEL,IF)='SYM'
           IF(CHOICE.EQ.'QUASI-SYMMETRY'  )CBC(ISIDE,IEL,IF)='QSM'
           IF(CHOICE.EQ.'FREE SURFACE'    )CBC(ISIDE,IEL,IF)='MS'
           IF(CHOICE.EQ.'OUTFLOW/N'       )CBC(ISIDE,IEL,IF)='ON'
           IF(CHOICE.EQ.'TRACTION'        )CBC(ISIDE,IEL,IF)='S'
           IF(CHOICE.EQ.'TRACTION (LOCAL)')CBC(ISIDE,IEL,IF)='SL'
           IF(CHOICE.EQ.'FLUID LAYERS'    )CBC(ISIDE,IEL,IF)='MSI'
           IF(CHOICE.EQ.'MELTING FRONT'   )CBC(ISIDE,IEL,IF)='MF'
           IF(CHOICE.EQ.'SOLID FRONT'     )CBC(ISIDE,IEL,IF)='MCI'
           IF(CHOICE.EQ.'LIQUID FRONT'    )CBC(ISIDE,IEL,IF)='MLI'
C
                 CBC3=CBC(ISIDE,IEL,IF)
                 BCLAB(1)= CBC3(1:1)
c
                     IF(CHOICE.EQ.'TEMP'         .OR.
     $               CHOICE.EQ.'FLUX'            .OR.
     $               CHOICE.EQ.'VELOCITY'        .OR.
     $               CHOICE.EQ.'VELOCITY (LOCAL)'.OR.
     $               CHOICE.EQ.'STRESS'          .OR.
     $               CHOICE.EQ.'STRESS   (LOCAL)'.OR.
     $               CHOICE.EQ.'SHEAR'           .OR.
     $               CHOICE.EQ.'SHEAR    (LOCAL)'.OR.
     $               CHOICE.EQ.'MOVING WALL     '.OR.
     $               CHOICE.EQ.'CONVECTION'      .OR.
     $               CHOICE.EQ.'OUTFLOW'         .OR.
     $               CHOICE.EQ.'OUTFLOW/N'       .OR.
     $               CHOICE.EQ.'TRACTION'        .OR.
     $               CHOICE.EQ.'TRACTION (LOCAL)'.OR.
     $               CHOICE.EQ.'RADIATION'       .OR.
     $               CHOICE.EQ.'FREE SURFACE'    .OR.
     $               CHOICE.EQ.'FLUID LAYERS'    .OR.
     $               CHOICE.EQ.'MELTING FRONT'   .OR.
     $               CHOICE.EQ.'SOLID FRONT'     .OR.
     $               CHOICE.EQ.'LIQUID FRONT'    )THEN
                        IF(CHOICE.EQ.'MOVING WALL') THEN
C                          Can only have fortran function
                           CHOICE='FORTRAN FUNCTION'
                        ELSE IF(CHOICE.EQ.'OUTFLOW')THEN
C                          Can only have constant
                           CHOICE='CONSTANT'
                        ELSE IF(CHOICE.EQ.'OUTFLOW/N')THEN
C                          Can only have constant
                           CHOICE='CONSTANT'
                        ELSE IF(CHOICE.EQ.'MELTING FRONT') THEN
C                          Can only have constant
                           CHOICE='CONSTANT'
                        ELSE
                           CALL PRS(
     $                     ' Choose Format for Relevant Parameters$')
                           NCHOIC=2
                           ITEM(1)='CONSTANT'
                           ITEM(2)='FORTRAN FUNCTION'
                           CALL MENU (XMOUSE,YMOUSE,BUTTON,
     $                     'FORTRAN FUNCTION')
                        ENDIF
                        IF(CHOICE.EQ.'FORTRAN FUNCTION')THEN
                         DO 111 I=1,3
                            III=ICHAR(CBC3(I:I))+32
                            IF(III.GE.97.AND.III.LE.122)THEN
                               CBC3(I:I)=CHAR(ICHAR(CBC3(I:I))+32)
                            ENDIF
111                      CONTINUE
                         CBC(ISIDE,IEL,IF)=CBC3
                         BCLAB(1)          = CBC3(1:1)
                         CALL INFLOW(IEL,ISIDE,IF,CBC(ISIDE,IEL,IF))
                         call rzero(bc(1,iside,iel,if),5)
                        ELSE IF(CHOICE.EQ.'CONSTANT')THEN
                           IF(BCCHOICE.EQ.'SHEAR    (LOCAL)')THEN
                             IF(IF3D)THEN
                               CALL PRS
     $                         (' Type #1-comp of Shear:$')
                               call rer(BC(2,ISIDE,IEL,IF))
                               CALL PRS
     $                         (' Type #2-comp of Shear:$')
                               call rer(BC(3,ISIDE,IEL,IF))
                             ELSE
                               CALL PRS
     $                         (' Type Shear:$')
                               call rer(BC(2,ISIDE,IEL,IF))
                             ENDIF
                           ELSE IF(BCCHOICE(10:16).EQ.'(LOCAL)')THEN
                              CALL PRS(' Type:$')
                              CALL PRS(' NORMAL-'//BCCHOICE//'$')
                              call rer(BC(1,ISIDE,IEL,IF))
                              IF(IF3D)THEN
                                 CALL PRS
     $                           (' #1 TANGENTIAL '//BCCHOICE//'$')
                                 call rer(BC(2,ISIDE,IEL,IF))
                                 CALL PRS
     $                           (' #2 TANGENTIAL '//BCCHOICE//'$')
                                 call rer(BC(3,ISIDE,IEL,IF))
                              ELSE
                                 CALL PRS
     $                           (' TANGENTIAL '//BCCHOICE//'$')
                                 call rer(BC(2,ISIDE,IEL,IF))
                              ENDIF
                           ELSE IF(BCCHOICE.EQ.'TEMP')THEN
                              CALL PRS('Type TEMP:$')
                              call rer(BC(1,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'FLUX')THEN
                              CALL PRS('Type FLUX:$')
                              call rer(BC(1,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'OUTFLOW/N' .OR.
     $                             BCCHOICE.EQ.'OUTFLOW')THEN
                              IF(.NOT.IFSTRS)THEN
                                 CALL PRS(' Exit pressure set to zero$')
                                 BC(1,ISIDE,IEL,IF) = 0.0
                              ELSE
                                 CALL PRS(
     $                           ' Type Exit Pressure:$')
                                 call rer(BC(1,ISIDE,IEL,IF))
                              ENDIF
                           ELSE IF(BCCHOICE.EQ.'FLUID LAYERS')THEN
                              WRITE(S,
     $                      '('' Type Surface Tension:'')')
                            CALL PRS(S//'$')
                            call rer(BC(4,ISIDE,IEL,IF))
                            CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
                            CBC(  ISIDEO,IELO,IF) = 'MPI'
                            BC(4,ISIDEO,IELO,IF)=BC(4,ISIDE,IEL,IF)
                           ELSE IF(BCCHOICE.EQ.'SHEAR')THEN
                             IF(IF3D)THEN
                               CALL PRS
     $                         (' Type X-comp of Shear:$')
                               call rer(BC(1,ISIDE,IEL,IF))
                               CALL PRS
     $                         (' Type Y-comp of Shear:$')
                               call rer(BC(2,ISIDE,IEL,IF))
                               CALL PRS
     $                         (' Type Z-comp of Shear:$')
                               call rer(BC(3,ISIDE,IEL,IF))
                             ELSE
                               CALL PRS
     $                         (' Type X-comp of Shear:$')
                               call rer(BC(1,ISIDE,IEL,IF))
                               CALL PRS
     $                         (' Type Y-comp of Shear:$')
                               call rer(BC(2,ISIDE,IEL,IF))
                             ENDIF
                           ELSE IF(BCCHOICE.EQ.'CONVECTION') THEN
                             CALL PRS
     $                       (' Heat Transfer Coefficient (h)>$')
C                            2nd storage Contains h
                             call rer(BC(2,ISIDE,IEL,IF))
                             CALL PRS
     $                       (' Temperature at Infinity (Tinf)>$')
                             call rer(BC(1,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'FREE SURFACE') THEN
                             CALL PRS
     $                       (' Pressure of Gas (Normal Traction)>$')
                             call rer(BC(1,ISIDE,IEL,IF))
                             CALL PRS
     $                       (' Shear (Tangential Traction)>$')
                             call rer(BC(2,ISIDE,IEL,IF))
                             IF(IF3D)THEN
                             CALL PRS
     $                       (' #2 Component of Shear '//
     $                       '(Tangential Traction)>$')
                             call rer(BC(3,ISIDE,IEL,IF))
                             ENDIF
                             CALL PRS(' Surface tension Coefficient>$')
                             call rer(BC(4,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'RADIATION') THEN
                             CALL PRS
     $                       (' Temperature at Infinity (Tinf)>$')
                             call rer(BC(1,ISIDE,IEL,IF))
                             CALL PRS(' Type Product of Emissivity$')
                             CALL PRS(' and Boltzmanns Constant >$')
                             call rer(BC(2,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'SOLID FRONT'
     $                     .OR.    BCCHOICE.EQ.'LIQUID FRONT')THEN
                             CALL PRS(' Freezing Temperature>$')
                             call rer(BC(4,ISIDE,IEL,IF))
                             CALL PRS(' Rho * Latent Heat>$')
                             call rer(BC(5,ISIDE,IEL,IF))
                           ELSE IF(BCCHOICE.EQ.'MELTING FRONT')THEN
C                            Fill up both sides of TEMPERATURE B.C.
                             CALL PRS(' Freezing Temperature>$')
                             call rer(BC(4,ISIDE,IEL,2))
                             CALL PRS(' Rho * Latent Heat>$')
                             call rer(BC(5,ISIDE,IEL,2))
C                            Export to overlapping edge
                             CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
                             BC (4,ISIDEO,IELO,2)=BC (4,ISIDE,IEL,2)
                             BC (5,ISIDEO,IELO,2)=BC (5,ISIDE,IEL,2)
                             IF(IIEL.GT.IEL)THEN
                                CBC(ISIDE ,IEL ,2)='MLI'
                                CBC(ISIDEO,IELO,2)='MCI'
                             ELSE
                                CBC(ISIDE ,IEL ,2)='MCI'
                                CBC(ISIDEO,IELO,2)='MLI'
                             ENDIF
                           ELSE
C                             B.C. requiring 3 x,y,z components
                              WRITE(S,'('' Type:'')')
                              CALL PRS(S//'$')
                              WRITE(S,'('' X-'',A26)')BCCHOICE
                              CALL PRS(S//'$')
                              call rer(BC(1,ISIDE,IEL,IF))
                              WRITE(S,'('' Y-'',A26)')BCCHOICE
                              CALL PRS(S//'$')
                              call rer(BC(2,ISIDE,IEL,IF))
                              IF(IF3D)THEN
                                 WRITE(S,'('' Z-'',A26)')BCCHOICE
                                 CALL PRS(S//'$')
                                 call rer(BC(3,ISIDE,IEL,IF))
                              ENDIF
                           ENDIF
                        ENDIF
                     ELSE IF(BCCHOICE.EQ.'PERIODIC-AUTO') THEN
C                    Automatic selection of Periodic bc
                        CALL FNDSIDA(JSIDE,JEL,ISIDE,IEL,IF)
                        IF (JSIDE.EQ.0) THEN
                         CALL PRS(
     $  'Cannot find a corresponding side, choose a new BC.$')
                         GOTO 63
                        ELSE
                          CALL LETBC(JSIDE,JEL,IF,BCLAB)
                          CALL LETBC(ISIDE,IEL,IF,BCLAB)
                          CBC( ISIDE,IEL,IF) = 'P  '
                          BC(1,ISIDE,IEL,IF) = JEL
                          BC(2,ISIDE,IEL,IF) = JSIDE
                          CBC( JSIDE,JEL,IF) = 'P  '
                          BC(1,JSIDE,JEL,IF) = IEL
                          BC(2,JSIDE,JEL,IF) = ISIDE
                        ENDIF
                     ELSE IF(BCCHOICE.EQ.'PERIODIC') THEN
C                    Periodic B.C.
C                    !!?? In general, check if old b.c. was a periodic or
C                    !!?? B.c. before writing over it!
                     CALL PRS
     $               ('Choose corresponding periodic with mouse$')
                     CALL PRS
     $               ('(use red menu area if side on other level)$')
C
                     go to 124
c23                  CALL PRS('ERROR: To input side, Use Mouse $')
124                  WRITE(S,'('' Periodic Side>'')')
                     CALL PRS(S//'$')
C                    Look for mouse input for periodic side
116                  CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
                     XP=XMOUSE
                     YP=YMOUSE
                     IF(BUTTON.NE.'LEFT') THEN
                        CALL PRS
     $                  ('Use left button to specify periodic side.$')
                        GO TO 124
                     ENDIF
                     IOTHER = 0
                     IF(XSCR(XMOUSE).GT.1.0)THEN
                        IOTHER = 1
                        CALL PRS('Type Level:$')
                        CALL REI(IPLEV)
                        IF(IFLEARN)WRITE(3,*) IPLEV
                        CALL PRS(
     $                  'Type Apartment letter (with correct case):$')
                        CALL RES(IPAPT,1)
                        IF(IFLEARN)WRITE(3,*) IPAPT
                        CALL PRS('Type Side number:$')
                        CALL REI(IPSIDE)
                        IF(IFLEARN)WRITE(3,*) IPSIDE
C                       Find Closest Side
C                       Use data from other level
                        ISIDEP = Ipside
                        DO 103 IIEL=1,NEL
                           IF(NUMAPT(IIEL).EQ.IPLEV .AND.
     $                        LETAPT(IIEL).EQ.IPAPT     ) IELP=IIEL
103                     CONTINUE
                        CALL PRSI('Displaying BC label From level$'
     $                  ,iplev)
                        CALL PRS('On this level.  Do not be alarmed.$')
                     ELSE
C                       Look on this level
                        RMIN=10000.0
                        DO 105 IIEL=1,NEL
                           DO 105 IISIDE=1,NSIDES
                              IF(NUMAPT(IIEL).NE.ILEVEL) GO TO 105
                              R=    ((XP-SIDES(IIEL,IISIDE,1))**2
     $                        +      (YP-SIDES(IIEL,IISIDE,2))**2)
                              IF(R.LT.RMIN) THEN
                                 RMIN=R
                                 ISIDEP = IISIDE
                                 IELP   = IIEL
                              ENDIF
105                     CONTINUE
                        IF(ISIDEP.EQ.5.OR.ISIDEP.EQ.6)THEN
C                          CHECK IF IT'S ABOVE OR BELOW CENTER
                           IF(YP.GT.SIDES(IELP,ISIDEP,2))THEN
                               ISIDEP=6
                           ELSE
                               ISIDEP=5
                           ENDIF
                        ENDIF
                        CALL PRSI('ISIDEP=$',ISIDEP)
                     ENDIF
                     IF(IELP.EQ.IEL .AND. ISIDEP.EQ.ISIDE) THEN
                        CALL PRS
     $                  ('** ERROR **  A side cannot be periodic$')
                        CALL PRS
     $                  ('with itself !  Choose the (different) $')
                        CALL PRS('side that the flashing side is $')
                        CALL PRS('periodic with.$')
                        go to 124
                     ENDIF
                     WRITE(S,110) IELP,ISIDEP,IEL,ISIDE
                     CALL PRS(S//'$')
110                  FORMAT(' ',2I7,'  IS PERIODIC WITH ',2I7)
C                    Now Enter Information into Periodicity Array
C                    !?? WHATS THE ORIENTATION??
                     IORIEN=1
                     CALL LETBC(ISIDEP,IELP,IF,BCLAB)
                     CBC( ISIDEP,IELP,IF) = 'P  '
                     BC(1,ISIDE ,IEL ,IF) = IELP
                     BC(2,ISIDE ,IEL ,IF) = ISIDEP
C                     BC(3,ISIDE ,IEL ,IF) = IORIEN
                     BC(1,ISIDEP,IELP,IF) = IEL
                     BC(2,ISIDEP,IELP,IF) = ISIDE
C                     BC(3,ISIDEP,IELP,IF) = IORIEN
                  ENDIF
C                 Post-pick warnings
                  IF(BCCHOICE.EQ.'AXIS')THEN
C                    Axisymmetric B.C.
                     IF(ABS(Y(IEL,ISIDE)/YFAC).GT.0.01)
     $               CALL PRS('** WARNING ** AXIS B.C. FOR Y=0 ONLY!$')
                     IF(ISIDE .NE. 1) THEN
                        CALL PRS('*** Error!!! Gauss-Jacobi Mesh$')
                        CALL PRS('allows only side 1 to be on axis.$')
                        CALL PRSI('Please delete element$',IEL)
                        CALL PRS('and Re-enter with axial side first.$')
                        CALL PRS('Hit <cr> to continue$')
                        CALL RES(LINE,0)
                     ENDIF
                  ENDIF
115               FORMAT(A1,'   <-- B.C. FOR ELE ',I3,' SIDE',I3)
                 ENDIF
                 CALL COLOR(3)
C                Turn Off Der Blinkenlights
                 IISIDE=ISIDE+1
                 IF(IISIDE.EQ.5) IISIDE = 1
                 CALL DRELEV(IFLASH,IFADE,'     ')
                 IF(ISIDE.EQ.5.OR.ISIDE.EQ.6) THEN
C                   Un-Flash Elevator at element floor or ceiling
                    IF(ISIDE.EQ.5)ILEV=ILEVEL
                    IF(ISIDE.EQ.6)ILEV=ILEVEL+1
                    CALL MOVESC(ELEC-ELEW,ELEB+ELEDY*(ILEV-1))
                    CALL DRAWSC(ELEC+ELEW,ELEB+ELEDY*(ILEV-1))
C                   Un-Flash mesh sides
                    IF(ISIDE.EQ.5)THEN
                       IED1=1
                       IED2=4
                    ELSE IF(ISIDE.EQ.6)THEN
                       IED1=1
                       IED2=4
                    ENDIF
                    CALL MOVEC(X(IEL, IED1),Y(IEL, IED1))
                    DO 230 IEDGE=IED1,IED2
                       CALL DRAWED(IEL,IEDGE,1)
230                 CONTINUE
                 ELSE
                    IF(IF3D)THEN
C                      Un-Flash at elevator at walls of current level
                       CALL MOVESC(ELEC     ,ELEB+ELEDY*(ILEVEL-1))
                       CALL DRAWSC(ELEC     ,ELEB+ELEDY*(ILEVEL  ))
                    ENDIF
C                   Un-Flash mesh side
                    CALL MOVEC(X(IEL, ISIDE),Y(IEL, ISIDE))
C                    CALL DRAW(X(IEL,IISIDE),Y(IEL,IISIDE))
                    CALL DRAWED(IEL,ISIDE,1)
                 ENDIF
                 IF(IF3D)THEN
C                   UN-FLASH SIDE ON ISOMETRIC
                    DO 233 ICORN=0,4
                        IF(ICORN.EQ.0) IC=FCORNS(4    ,ISIDE)
                        IF(ICORN.NE.0) IC=FCORNS(ICORN,ISIDE)
                        XISOM=XPHY(XSCR(X(IEL,IC))/5.0 + 0.8)
     $                  -Z(IEL,IC)/20.
                        YISOM=YPHY(YSCR(Y(IEL,IC))/5.0 + 0.8)
     $                  -Z(IEL,IC)/20.
                        IF(ICORN.EQ.0) CALL MOVEC(XISOM,YISOM)
                        IF(ICORN.NE.0) CALL DRAWC(XISOM,YISOM)
233                 CONTINUE
                 ENDIF
C                Label Boundary Side
                 CALL LETBC(ISIDE,IEL,IF,BCLAB)
900         CONTINUE
C
C==================================================================
C           Entering MODIFY BC MENU area (apparently)
C==================================================================
C
313         CONTINUE
            ITEM(1)='END  LEVEL'
            ITEM(2)='MODIFY B.C.'
            ITEM(3)='SHOW   B.C.'
            ITEM(4)='USER   B.C.'
            IF(IF.eq.1 .OR. IF.EQ.2 .AND. .NOT. IFFLOW)THEN
               ITEM(5)='INTERNAL B.C.'
               NCHOIC=5
            ELSE
               NCHOIC=4
            ENDIF
            NCHOIC=NCHOIC+1
            ITEM(NCHOIC)='ZOOM'
            CALL MENU (XMOUSE,YMOUSE,BUTTON,'MODIFY B.C.')
C
            IF (CHOICE.EQ.'USER   B.C.') then
               call user_bc
               return
            endif
c
            IF(CHOICE.EQ.'INTERNAL B.C.') THEN
               CMENU = 'INTERNAL'
            ELSE
               CMENU = ' '
            ENDIF
            IF(CHOICE.EQ.'MODIFY B.C.' .OR. CHOICE.EQ.'SHOW   B.C.'
     $     .OR.CHOICE.EQ.'INTERNAL B.C.')THEN
                CALL PRS('Enter side with mouse$')
                CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
C               Find closest BOUNDARY side
                XP=XMOUSE
                YP=YMOUSE
                RMIN=10000.0
                DO 106 IIEL=1,NEL
                   DO 106 IISIDE=1,NSIDES
                      IF(NUMAPT(IIEL).EQ.ILEVEL .AND.
     $                (CBC(IISIDE,IIEL,IF).NE.'E'.OR.
     $                 CHOICE.EQ.'INTERNAL B.C.'))THEN
                         R=    ((XP-SIDES(IIEL,IISIDE,1))**2
     $                   +      (YP-SIDES(IIEL,IISIDE,2))**2)
                         IF(R.LT.RMIN) THEN
C                         One last check: above or below center- 5 or 6
                          IF(IISIDE.EQ.5.AND.YP.GT.SIDES(IIEL,IISIDE,2))
     $                    GO TO 1113
                          IF(IISIDE.EQ.6.AND.YP.LT.SIDES(IIEL,IISIDE,2))
     $                    GO TO 1113
                            RMIN=R
                            ISIMOD = IISIDE
                            IELMOD = IIEL
1113                        CONTINUE
                         ENDIF
                      ENDIF
106             CONTINUE
c
                IF(CHOICE.EQ.'SHOW   B.C.')THEN
C                  Show   B.C.
                   WRITE(S,'(A7,I4,A6,I3,A1)')'Element',IELMOD,
     $             '  Side',ISIMOD,'$'
                   CALL PRS(S)
                   WRITE(S,'(a5,5g12.5,a1)')
     $             CBC(isimod,IELMOD,IF),
     $             BC(1,ISIMOD,IELMOD,IF),
     $             BC(2,ISIMOD,IELMOD,IF),
     $             BC(3,ISIMOD,IELMOD,IF),
     $             BC(4,ISIMOD,IELMOD,IF),
     $             BC(5,ISIMOD,IELMOD,IF),'$'
                   CALL PRS(S)
                   GO TO 313
                ELSE IF(CHOICE.EQ.'MODIFY B.C.'
     $           .OR.   CHOICE.EQ.'INTERNAL B.C.')THEN
C                  Reset to zero and send back into loop.  Loop will demand bc
                   IF(CBC(isimod,ielmod,IF).EQ.'P  ')THEN
C                     We need to zero out the periodic side also
                      IELPER=BC(1,ISIMOD,IELMOD,IF)
                      ISIPER=BC(2,ISIMOD,IELMOD,IF)
                      CALL LETBC(ISIPER,IELPER,IF,' ')
                      CBC( ISIPER,IELPER,IF)=' '
                      BC(1,ISIPER,IELPER,IF)=0
                      BC(2,ISIPER,IELPER,IF)=0
                      BC(3,ISIPER,IELPER,IF)=0
                      BC(4,ISIPER,IELPER,IF)=0
                      BC(5,ISIPER,IELPER,IF)=0
                   ENDIF
                   CALL LETBC(ISIMOD,IELMOD,IF,' ')
                   CBC( isimod,ielmod,IF)=' '
                   IF(CHOICE.NE.'INTERNAL B.C.')THEN
C                     Retain old elemental connectivity for internal B.C.'s
                      BC(1,ISIMOD,IELMOD,IF)=0
                      BC(2,ISIMOD,IELMOD,IF)=0
                      BC(3,ISIMOD,IELMOD,IF)=0
                      BC(4,ISIMOD,IELMOD,IF)=0
                      BC(5,ISIMOD,IELMOD,IF)=0
                   ENDIF
                   GO TO 800
                ENDIF
             ELSE IF(CHOICE.EQ. 'END  LEVEL') THEN
C               First check to see that all elements on this level are
C               Consistent in coordinates; a given element cannot have both
C               local and global definition of velocity and stress.
C
                DO 1500 IEL=1,NEL
                   IF(NUMAPT(IEL).EQ.ILEVEL)THEN
                      ITYPE=0
                      DO 1200 ISIDE=1,NSIDES
                         IF( CBC(ISIDE,IEL,IF).EQ.'VL'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'vl'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'SL'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'sl')THEN
C                           There is at least one BC in local coordinates
                            IF(ITYPE.EQ.1)THEN
C                              There is a conflict
                               CALL PRS
     $                         ('Error: You cannot mix local and'
     $                         //' cartesian coordinates in B.C.''s in '
     $                         //'the same element.$')
                               WRITE(S,'(1X,A7,I4,A2,I3,A2,A2)')
     $                         'Element',IEL,' [',NUMAPT(IEL),
     $                         LETAPT(IEL),']$'
                               CALL PRS(S)
                               CALL PRS(' has such a mixture.$')
                               CALL PRS(
     $'Please MODIFY velocity and stress B.C.''s to consistent '
     $//'coordinates.$')
                               GO TO 1010
                            ELSE
                               ITYPE=2
                            ENDIF
                         ENDIF
                         IF( CBC(ISIDE,IEL,IF).EQ.'V'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'v'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'S'
     $                   .OR.CBC(ISIDE,IEL,IF).EQ.'s')THEN
C                           There is at least one BC in global coordinates
                            IF(ITYPE.EQ.2)THEN
C                              There is a conflict
                               CALL PRS
     $                         ('Error: You cannot mix local and'
     $                         //' cartesian coordinates in B.C.''s in '
     $                         //'the same element.$')
                               WRITE(S,'(1X,A7,I4,A2,I3,A2,A2)')
     $                         'Element',IEL,' [',NUMAPT(IEL),
     $                         LETAPT(IEL),']$'
                               CALL PRS(S)
                               CALL PRS('has such a mixture.$')
                               CALL PRS('Please MODIFY velocity '
     $                         //'and stress B.C.''s to consistent '
     $                         //'coordinates.$')
                               GO TO 1010
                            ELSE
                               ITYPE=1
                            ENDIF
                         ENDIF
1200                  CONTINUE
                   ENDIF
1500            CONTINUE
             ELSE IF(CHOICE.EQ. 'ZOOM') THEN
                CALL SETZOOM
                CALL REFRESH
                CALL DRMENU('NOCOVER')
                CALL DRGRID
                DO 1600 IEL=1,NEL
                   CALL DRAWEL(IEL)
 1600           CONTINUE
                GOTO 313
             ENDIF
C            Just keep on going in loop (to 1000)
1000     CONTINUE
1999  CONTINUE
      ILEVEL=NLEVEL
2000  CONTINUE
c
c     Check and reset stray period bcs  , pff  4/7/99
c
      DO IF=IF1,NFLDS
         call period_check(if)
      ENDDO
c
      RETURN
c
c     This is the auto-bc setting loop  pff
c
      END
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine overlap(iel,iside,ielo,isideo)
C
      include 'basics.inc'

C     Find closest element, side we want to duplicate, If there is no overlap,
C     returns original element and side
C
      RMIN=10000.0 * XFAC
      RMAX=    0.0
C
      IELO   = IEL
      ISIDEO = ISIDE
      DO 3 IIEL=1,NEL
         DO 1 IISIDE=1,NSIDES
            R=    (
     $       (SIDES(IEL,ISIDE,1)-SIDES(IIEL,IISIDE,1))**2
     $      +(SIDES(IEL,ISIDE,2)-SIDES(IIEL,IISIDE,2))**2
     $      +(SIDES(IEL,ISIDE,3)-SIDES(IIEL,IISIDE,3))**2
     $      )
            IF(IEL.NE.IIEL.OR.IISIDE.NE.ISIDE)THEN
               IF(R.GT.RMAX) RMAX = R
               IF(R.LT.RMIN) THEN
                  RMIN   = R
                  ISIDEO = IISIDE
                  IELO   = IIEL
               ENDIF
            ENDIF
1        CONTINUE
3     CONTINUE
C
      IF(RMIN .EQ. 10000.0 * XFAC)THEN
         CALL PRS('ERROR FINDING OVERLAPPING SIDE FOR B.C.$')
      ELSE IF (RMIN .GE. RMAX/100.)THEN
         CALL PRS('ERROR FINDING OVERLAPPING SIDE FOR B.C.$')
         CALL PRS('NO ADJACENT ELEMENT$')
         IELO   = IEL
         ISIDEO = ISIDE
      ENDIF
C
      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine INFLOW(IEL,ISIDE,IF,CBCI)
C! REAL KLUDGE USED NOT TO CORRUPT LINE  !!??DO THE POINTERS GET READ IN RIGHT?
      include 'basics.inc'
      CHARACTER*3 CBCI
C
C     Don't read inline stuff anymore... pff 8/30/93.
C
      if (nel.gt.0) return
C
      CBC3 = CBCI
      IF(CBCI .NE. 'on ')THEN
         CALL PRS(' Enter boundary conditions for this element $')
         CALL PRS
     $   (' as a function of X,Y [Z]& TIME in the style of a fortran$')
         CALL PRS(' statement.  Start in column 1.  Use one line for $')
      ENDIF
      IF(CBCI.EQ.'s  ')THEN
C        stress (Traction)
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRSIS(' each of the $',NLINES,
     $   ' TRACTION components$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TRY=SIN(X)$')
         IF(IF3d)CALL PRS(' TRZ=1.0$')
      ELSE IF(CBCI.EQ.'sl ')THEN
C        stress (Traction)
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRSIS(' each of the $',NLINES,
     $   ' TRACTION components$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRN=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
         IF(IF3D)CALL PRS(' TR2=SIN(X)$')
      ELSE IF(CBCI.EQ.'sh ')THEN
C        shear
         IF(.NOT.IF3D)THEN
            NLINES=2
            CALL PRS(' each of the 2 shear (TRaction) components $')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TRY=    Y  +  3 * TIME)$')
         ELSE IF(IF3D)THEN
            NLINES=3
            CALL PRS(' each of the 3 shear (TRaction) componnts$')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TRY=    Y  +  3 * TIME)$')
            CALL PRS(' TRZ=    Z  +  3 * TIME)$')
         ENDIF
      ELSE IF(CBCI.EQ.'shl')THEN
C        shear local
         IF(.NOT.IF3D)THEN
            NLINES=1
            CALL PRS(' shear (TRaction 1) $')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
         ELSE IF(IF3D)THEN
            NLINES=2
            CALL PRS(' each of the 2 shear (TRaction) componnts$')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TR2=X      + EXP( TIME)$')
         ENDIF
      ELSE IF(CBCI.EQ.'sl ' .OR. CBCI .EQ. 'ms ')THEN
C        stress
         IF(.NOT.IF3D)NLINES=3
         IF(     IF3D)NLINES=4
         CALL PRSIS(' each of the $',NLINES,
     $   ' stress (TRACTION) components $')
         CALL PRS(' And Surface Tension$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRN=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TR1=SIN(X)$')
         IF(IF3D)CALL PRS(' TR2=1.0$')
         CALL PRS(' SIGMA=1.2 * TEMP$')
      ELSE IF(CBCI .EQ. 'msi')THEN
C        Fluid Layers
         NLINES = 1
         CALL PRS(' the surface tension coefficient$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' SIGMA = 1.2 * TEMP$')
      ELSE IF(CBCI .EQ. 'on ')THEN
C        Outflow
         NLINES = 1
         CALL PRS(' Enter exit pressure (TRactioN)function$')
         CALL PRS(' EXAMPLE:  TRN = 2. * TIME $')
      ELSE IF(CBCI .EQ. 'mf ')THEN
C        Melting Front
         NLINES = 2
         CALL PRS(' the freezing temp and the product Rho*Latent heat$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TEMP??? What is it Lee? = 1.2$')
         CALL PRS(' RHOL =What is it LEE??? 3 + $')
      ELSE IF(CBCI(1:1).EQ.'f')THEN
C        flux
         NLINES=1
         CALL PRS(' the flux.  EXAMPLE: FLUX=1.0$')
      ELSE IF(CBCI(1:1).EQ.'t')THEN
C        temperature
         NLINES=1
         CALL PRS(' TEMP.  Example:$')
         CALL PRS(' TEMP = 0.667 * (1.0-Y**2) +EXP(TIME)$')
         CALL PRS(' Enter statement For TEMP:$')
      ELSE IF(CBCI(1:1).EQ.'v'.OR.CBCI(1:1).EQ.'m')THEN
C        velocity
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRS(' each of the $',NLINES,' velocity components as $')
         CALL PRS(' functions of space and time.  Example:$')
         IF(CBCI .EQ. 'vl') THEN
            CALL PRS(' UN=SIN(X) + EXP(-TIME)$')
            CALL PRS(' U1=SIN(Y) + EXP(-TIME)$')
            IF(IF3D)
     $      CALL PRS(' U2 = 1 - TIME$')
         ELSE
            CALL PRS(' UX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' UY=1$')
            IF(IF3D)
     $      CALL PRS(' UZ = 0.0$')
         ENDIF
      ELSE IF(CBCI(1:1).EQ.'c')THEN
C        convection
         NLINES=2
         CALL PRS(' Tinf & one line for HC.  Example:$')
         CALL PRS(' TINF=SIN(X) + EXP(TIME)$')
         CALL PRS(' HC=5.0$')
      ELSE IF(CBCI(1:1).EQ.'r')THEN
         CALL PRS(
     $   ' The HRAD (Sigma*F) and TINF (the absolute Temperature).$')
         CALL PRS('  EXAMPLE: HRAD = 1.0e-4$')
         NLINES = 2
      ELSE
         CALL PRS('New BC Type.$')
         CALL PRS('Enter BC in the Following three lines $')
         CALL PRS('Ignoring Prompts for components$')
         NLINES=3
      ENDIF
      IF(CBC3(3:3) .EQ. 'i')THEN
         BC(4,ISIDE,IEL,IF)=NLINES
         BC(5,ISIDE,IEL,IF)=LOCLIN
      ELSE
         BC(1,ISIDE,IEL,IF)=NLINES
         BC(2,ISIDE,IEL,IF)=LOCLIN
      ENDIF
      DO 100 I=1,NLINES
         IF(NLINES.GT.1 .AND. CBCI .eq. 'ms ')THEN
            IF(I.EQ.1)CALL PRS(' Normal Component:$')
            IF(IF3d)THEN
               IF(I.EQ.2)CALL PRS(' #1 Tangential Component:$')
               IF(I.EQ.3)CALL PRS(' #2 Tangential Component:$')
               IF(I.EQ.4)CALL PRS(' Surface Tension SIGMA:$')
            ELSE
               IF(I.EQ.2)CALL PRS(' #1 Tangential Component:$')
               IF(I.EQ.3)CALL PRS(' Surface Tension SIGMA:$')
            ENDIF
         ELSE IF(NLINES.GT.1 .AND. CBCI .eq. 'shl')THEN
            IF(I.EQ.1)CALL PRS(' 1st component:$')
            IF(I.EQ.2)CALL PRS(' 2nd component:$')
         ELSE IF(NLINES.NE.1 .AND.(CBCI(2:2).EQ.'l'
     $   .OR.CBCI(3:3).eq.'l'))THEN
            IF(I.EQ.1)CALL PRS(' Enter Normal component:$')
            IF(I.EQ.2.AND.NLINES.EQ.2)
     $      CALL PRS(' Enter Tangential component:$')
            IF(I.EQ.2.AND.NLINES.EQ.3)
     $      CALL PRS(' Enter Tangential component #1:$')
            IF(I.EQ.3)CALL PRS(' Enter Tangential component #2:$')
         ELSE IF(NLINES.NE.1 .AND.CBCI(1:1).EQ.'r'.AND.I.EQ.2)THEN
            CALL PRS(' Now TINF.   EXAMPLE: TINF = 1.0e-4$')
         ELSE IF(NLINES.NE.1 .AND.CBCI(1:1).NE.'c')THEN
            IF(I.EQ.1)CALL PRS(' Enter X-component:$')
            IF(I.EQ.2)CALL PRS(' Enter Y-component:$')
            IF(I.EQ.3)CALL PRS(' Enter Z-component:$')
         ENDIF
 99      CALL RES(INBC(LOCLIN),60)
         IF(IFLEARN)WRITE(3,'(A60)')INBC(loclin)
         IF(INBC(LOCLIN) .EQ. ' ')THEN
            CALL PRS('Please enter a non-blank fortran line:$')
            GO TO 99
         ENDIF
         LOCLIN=LOCLIN+1
100   CONTINUE
C
      IF(CBC3(3:3) .EQ. 'i')THEN
C        EXPORT B.C. FOR INTERNAL Boundary
         CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
         CBC(  ISIDEO,IELO,IF) =CBC(  ISIDE,IEL,IF)
         IF(CBC3 .EQ. 'msi') CBC(  ISIDEO,IELO,IF) = 'mpi'
         BC (4,ISIDEO,IELO,IF) = BC(4,ISIDE,IEL,IF)
         DO 121 I=LOCLIN-1,LOCLIN-1+NLINES
            INBC(I+NLINES)=INBC(I)
121      CONTINUE
          LOCLIN=LOCLIN+NLINES
          BC(5,ISIDEO,IELO,IF) = LOCLIN
      ENDIF

      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine LETBC(ISIDE,IEL,IF,BCLAB)
      include 'basics.inc'
      CHARACTER BCLAB,BCLAB2*2
      BCLAB2(1:1)=BCLAB
      BCLAB2(2:2)='$'
      YSHIFT= 0.0
      IF(ISIDE.EQ.5) YSHIFT=-YFAC*0.05
      IF(ISIDE.EQ.6) YSHIFT= YFAC*0.02
      XSHIFT= XFAC*(-3.0/200.) +(IF-1)*XFAC/50.0
      XLAB=SIDES(IEL,ISIDE,1)+XSHIFT
      YLAB=SIDES(IEL,ISIDE,2)+YSHIFT
c     CALL PRSii(' B.C.?$',iel,iside)
      IF(BCLAB.EQ.' ')THEN
C        Draw empty box
         DX=XFAC/50.
         DY=YFAC/40.
         CALL FILLP(-15)
         CALL BEGINP(XLAB,YLAB)
         CALL DRAW(XLAB+DX,YLAB   )
         CALL DRAW(XLAB+DX,YLAB+DY)
         CALL DRAW(XLAB   ,YLAB+DY)
         CALL DRAW(XLAB   ,YLAB   )
         CALL ENDP
      ELSE
C        Normal Write
         IF(CWRITE.NE.0. .AND.BCLAB2 .NE.'E$')
     $   CALL GWRITE(XLAB+XFAC/300.,YLAB+YFAC/300.,1.0,BCLAB2)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine CHKBCS
      include 'basics.inc'
      CHARACTER*3 CBCTMP
      CHARACTER*1 YESNO
      LOGICAL IFFAIL
C
      IF=1
      DO 100 IEL=1,NELF
         DO 100 ISIDE=1,NSIDES
            IF( CBC(ISIDE,IEL,IF).EQ.'MSI'
     $      .OR.CBC(ISIDE,IEL,IF).EQ.'msi')THEN
               CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
               IF(CBC (ISIDEO,IELO,IF).NE. 'MPI'
     $        .and.CBC(ISIDEO,IELO,IF).NE. 'mpi')THEN
                 CALL BEEP
                 CALL PRS('Inconsistency in Fluid Layer B.C.s$')
                 CALL PRS('Attempting Correction....$')
                 IF(CBC(ISIDE,IEL,IF).EQ.'MSI')
     $              CBC(ISIDEO,IELO,IF)='MPI'
                 IF(CBC(ISIDE,IEL,IF).EQ.'msi')
     $              CBC(ISIDEO,IELO,IF)='mpi'
               ENDIF
               IF(IGROUP(IELO) .EQ. IGROUP(IEL))THEN
                  CALL BEEP
                  CALL PRS('Inconsistency in Fluid Layer B.C.s$')
                  CALL PRII(IEL,IELO)
                  CALL PRS
     $            ('Are 2 elements of same fluid that straddle a$')
                  CALL PRS
     $            ('Fluid Layer.  Different fluids (groups) are reqd$')
               ENDIF
               IF(IGROUP(IELO) .GT. IGROUP(IEL))THEN
C                 Swap
                  CBCTMP              = CBC(ISIDEO,IELO,IF)
                  CBC(ISIDEO,IELO,IF) = CBC(ISIDE ,IEL ,IF)
                  CBC(ISIDE ,IEL ,IF) = CBCTMP
                  CALL PRS('Normal Fluid layer B.C. swap$')
               ENDIF
            ENDIF
  100 CONTINUE
C
C     Check for curve side consistency
C
  101 CONTINUE
      IF (IFFLOW) Ifld=1
      IF (IFHEAT) Ifld=2
C
      IFfail = .false.
      DO 400 Ie=1,Nel
      DO 400 ISIDE=1,NSIDES
         IF (CBC(iside,ie,Ifld).eq.'E'     .or.
     $       CBC(iside,ie,Ifld).eq.'P  ' )   THEN
C
             Je    = INT ( BC(1,ISIDE,IE,Ifld) )
             Jside = INT ( BC(2,ISIDE,IE,Ifld) )
C
C           Check edges
            IF (ISIDE.le.4) THEN
             Do 250 Ilev=0,4,4
               Ied = Iside + Ilev
               IF (CCURVE(Ied,Ie).eq.'C') THEN
                Jed = Jside+Ilev
                IF (CCURVE(Jed,Je).ne.CCURVE(Ied,Ie)    .or.
     $             -CURVE(1,Jed,Je).ne.CURVE(1,Ied,Ie)) THEN
                 call qchk(ie,iside,je,jside,IF)
                 IFfail = .true.
                 CALL PRS('ERROR: Curve side consistency failure.$')
  200            CALL PRS(
     $          'Enter element number to be changed (other=abort).$')
                 WRITE(S,210) Ie,Ied,CCURVE(Ied,Ie),Curve(1,Ied,Ie)
                 CALL PRS(S)
                 WRITE(S,210) Je,Jed,CCURVE(Jed,Je),Curve(1,Jed,Je)
                 CALL PRS(S)
  210            FORMAT(
     $           ' El:',I5,3x,'Edge:',I2,3x,'C: ',A3,' Rad:',G14.6,'$')
cccc             CALL REI(Ke)
                 IF (Ke.eq.Ie) THEN
                    CCURVE(Ied,Ie)   = CCURVE(Jed,je)
                    CURVE (1,Ied,Ie) = -CURVE (1,Jed,je)
                 ELSE IF(Ke.eq.Je) THEN
                    CCURVE(Jed,Je)   = CCURVE(Ied,Ie)
                    CURVE (1,Jed,Je) = -CURVE (1,Ied,Ie)
                 ELSE
                  CALL PRS(
     $            'Are you sure you want to keep it this way?$')
cccc              CALL RES(YESNO,1)
cccc              IF (YESNO.ne.'y'.and.YESNO.ne.'Y') GOTO 200
                 ENDIF
                ENDIF
               ENDIF
  250        CONTINUE
            ENDIF
C
            IF (CCURVE(Iside,Ie).eq.'s') THEN
               IF (CCURVE(Jside,Je).ne.CCURVE(Iside,Ie)   .or.
     $            CURVE(1,Jside,Je).ne.CURVE(1,Iside,Ie)  .or.
     $            CURVE(2,Jside,Je).ne.CURVE(2,Iside,Ie)  .or.
     $            CURVE(3,Jside,Je).ne.CURVE(3,Iside,Ie)  .or.
     $            CURVE(4,Jside,Je).ne.CURVE(4,Iside,Ie)) THEN
                 IFfail = .true.
                 call qchk(ie,iside,je,jside,IF)
                 CALL PRS('ERROR: Curve side consistency failure.$')
  300            CALL PRS(
     $          'Enter element number to be changed (other=abort).$')
                 WRITE(S,310) Ie,Iside,CCURVE(Iside,Ie),
     $                        (Curve(j,Iside,Ie),j=1,4)
                 CALL PRS(S)
                 WRITE(S,310) Je,Jside,CCURVE(Jside,Je),
     $                        (Curve(j,Jside,Je),j=1,4)
                 CALL PRS(S)
  310            FORMAT(
     $           ' El:',I5,3x,'F:',I2,3x,'c: ',A3,4G13.5,'$')
cccc             CALL REI(Ke)
                 IF (Ke.eq.Ie) THEN
                    CCURVE(Iside,Ie)   = CCURVE(Jside,je)
                    CURVE (1,Iside,Ie) = CURVE (1,Jside,je)
                    CURVE (2,Iside,Ie) = CURVE (2,Jside,je)
                    CURVE (3,Iside,Ie) = CURVE (3,Jside,je)
                    CURVE (4,Iside,Ie) = CURVE (4,Jside,je)
                    CURVE (5,Iside,Ie) = CURVE (5,Jside,je)
                 ELSE IF(Ke.eq.Je) THEN
                    CCURVE(Jside,Je)   = CCURVE(Iside,Ie)
                    CURVE (1,Jside,Je) = CURVE (1,Iside,Ie)
                    CURVE (2,Jside,Je) = CURVE (2,Iside,Ie)
                    CURVE (3,Jside,Je) = CURVE (3,Iside,Ie)
                    CURVE (4,Jside,Je) = CURVE (4,Iside,Ie)
                    CURVE (5,Jside,Je) = CURVE (5,Iside,Ie)
                 ELSE
                  CALL PRS(
     $            'Are you sure you want to keep it this way?$')
cccc              CALL RES(YESNO,1)
cccc              IF (YESNO.ne.'y'.and.YESNO.ne.'Y') GOTO 300
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
  400 CONTINUE
      IF (.not.IFfail) CALL PRS('Curve side consistency check OK.$')
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine fndsida(kside,ke,iside,ie,if)
      include 'basics.inc'
      real vec1(3),vec2(3),vec3(3)
C
C     Find a side JSIDE,JEL which corresponds to ISIDE,IE and doesn't
C     have BC's set already
C
      NSIDES=NDIM*2
C
C     Define Normal to this plane (side), and find the element
C     center which also lies on this line, and is farthest from
C     the current element.
      CALL MKSIDE
      CALL GENCEN
C
C     Find Sides' Normal Vector at midpoint
C
      IC1=ISIDE
      IC2=ISIDE+1
      IF(ISIDE.EQ.4)IC2=1
C     This stuff only relevant for 3d
      IC3=IC2+4
      IC4=IC1+4
      IF (ISIDE.EQ.5) THEN
         IC1=1
         IC2=2
         IC3=3
         IC4=4
      ELSE IF(ISIDE.EQ.6) THEN
         IC1=1+4
         IC2=2+4
         IC3=3+4
         IC4=4+4
      ENDIF
      write(s,9) ie,iside,ic1,ic2,ic3,ic4
    9 format(' iesc14:',6i5,'$')
      if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
      IF (IF3D) THEN
         VEC1(1)=X(IE,IC3)-X(IE,IC1)
         VEC1(2)=Y(IE,IC3)-Y(IE,IC1)
         VEC1(3)=Z(IE,IC3)-Z(IE,IC1)
C
         VEC2(1)=X(IE,IC4)-X(IE,IC2)
         VEC2(2)=Y(IE,IC4)-Y(IE,IC2)
         VEC2(3)=Z(IE,IC4)-Z(IE,IC2)
C
         CALL CROSS(VEC3,VEC1,VEC2)
      ELSE
         VEC1(1)=X(IE,IC2)-X(IE,IC1)
         VEC1(2)=Y(IE,IC2)-Y(IE,IC1)
C
         VEC3(1)= VEC1(2)
         VEC3(2)=-VEC1(1)
         VEC3(3)=0.0
      ENDIF
      CALL NORM3D(VEC3)
c
      if (if_lattice) then          ! Make vector point in lattice direction
         d1 = dotprod(vec3,ulat1)
         d2 = dotprod(vec3,ulat2)
         d3 = dotprod(vec3,ulat3)
         if (d1.le.d2.and.d1.le.d3) then
            call copy(vec3,ulat1,3)
         elseif (d2.le.d3) then
            call copy(vec3,ulat2,3)
         else
            call copy(vec3,ulat3,3)
         endif
      endif
c
      write(s,19) vec3(1),vec3(2),vec3(3)
   19 format(' vec3:',3e13.5,'$')
      if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
c
c     EPS2=(EPSM*RCEN(IE))**2
      EPSM=1.e-4
      EPS2=EPSM**2
c
      KE=0
      KSIDE=0
      DSTMAX=0.0
c
      DO 100 JE=1,NEL
      DO 100 JSIDE=1,NSIDES
C
C        Don't find yourself.
         IF (JE.EQ.IE.AND.JSIDE.EQ.ISIDE) GOTO 100
C
C        Don't find a bc which is already defined. (to be 'E')
         IF (CBC(JSIDE,JE,IF).EQ.'E') GOTO 100
C
C        OK, is the center of this element side close to the line?
C
         VEC1(1)=SIDES(IE,ISIDE,1)-SIDES(JE,JSIDE,1)
         VEC1(2)=SIDES(IE,ISIDE,2)-SIDES(JE,JSIDE,2)
         VEC1(3)=SIDES(IE,ISIDE,3)-SIDES(JE,JSIDE,3)
         DISTP=DOTPROD(VEC1,VEC1)
         DIST2=DISTP-(DOTPROD(VEC1,VEC3))**2
         if (distp.gt.0) DIST2=DIST2/DISTP
c        write(6,11) dist2,distp,eps2,dstmax,je,jside
c  11    format(1p4e12.4,2i6,' dist2')
         IF (DIST2.LE.EPS2.AND.DISTP.GT.DSTMAX) THEN
            KE=JE
            KSIDE=JSIDE
            DSTMAX=DISTP
         ENDIF
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine qchk(ie,iside,je,jside,IF)
      include 'basics.inc'
      write(6,9) ie,iside,cbc(iside,ie,IF),
     $            (bc(j,iside,ie,IF),j=1,5)
      write(6,9) je,jside,cbc(jside,je,IF),
     $            (bc(j,jside,je,IF),j=1,5)
    9 format(' ie,side,c,b:',i5,i2,1x,a3,1x,5g11.3)
      return
      end
c-----------------------------------------------------------------------
      function i_periodic(iside,ie,ifld,bclab)
c
c     Assign periodic-auto bc
c
      include 'basics.inc'
      character*1 bclab(2)
c
      i_periodic=0
c
      call fndsida(jside,jel,iside,ie,ifld)
      if (jside.eq.0) then
        CALL prsii('For element,side:',ie,iside)
        CALL PRS('Cannot find a corresponding side, choose a new BC.$')
        i_periodic=1
        return
      endif
c
      call letbc(jside,je,ifld,bclab)
      call letbc(iside,ie,ifld,bclab)
      cbc( iside,ie,ifld) = 'P  '
      bc(1,iside,ie,ifld) = je
      bc(2,iside,ie,ifld) = jside
      cbc( jside,je,ifld) = 'P  '
      bc(1,jside,je,ifld) = ie
      bc(2,jside,je,ifld) = iside
c
      return
      end
c-----------------------------------------------------------------------
      subroutine getnormals(unx,uny,unz,ie,is,x,y,z,ld,ndim)
      real x(ld,8),y(ld,8),z(ld,8)
c
      real u(3),v(3),w(3)
c
c     Base normals on cross-product of diagonals
c
      if (is.eq.1) then
         u(1) = x(ie,6)-x(ie,1)
         u(2) = y(ie,6)-y(ie,1)
         u(3) = z(ie,6)-z(ie,1)
c
         v(1) = x(ie,5)-x(ie,2)
         v(2) = y(ie,5)-y(ie,2)
         v(3) = z(ie,5)-z(ie,2)
      elseif(is.eq.2) then
         u(1) = x(ie,7)-x(ie,2)
         u(2) = y(ie,7)-y(ie,2)
         u(3) = z(ie,7)-z(ie,2)
c
         v(1) = x(ie,6)-x(ie,3)
         v(2) = y(ie,6)-y(ie,3)
         v(3) = z(ie,6)-z(ie,3)
      elseif(is.eq.3) then
         u(1) = x(ie,8)-x(ie,3)
         u(2) = y(ie,8)-y(ie,3)
         u(3) = z(ie,8)-z(ie,3)
c
         v(1) = x(ie,7)-x(ie,4)
         v(2) = y(ie,7)-y(ie,4)
         v(3) = z(ie,7)-z(ie,4)
      elseif(is.eq.4) then
         u(1) = x(ie,5)-x(ie,4)
         u(2) = y(ie,5)-y(ie,4)
         u(3) = z(ie,5)-z(ie,4)
c
         v(1) = x(ie,8)-x(ie,1)
         v(2) = y(ie,8)-y(ie,1)
         v(3) = z(ie,8)-z(ie,1)
      elseif(is.eq.5) then
         u(1) = x(ie,4)-x(ie,2)
         u(2) = y(ie,4)-y(ie,2)
         u(3) = z(ie,4)-z(ie,2)
c
         v(1) = x(ie,3)-x(ie,1)
         v(2) = y(ie,3)-y(ie,1)
         v(3) = z(ie,3)-z(ie,1)
      elseif(is.eq.6) then
         u(1) = x(ie,7)-x(ie,5)
         u(2) = y(ie,7)-y(ie,5)
         u(3) = z(ie,7)-z(ie,5)
c
         v(1) = x(ie,8)-x(ie,6)
         v(2) = y(ie,8)-y(ie,6)
         v(3) = z(ie,8)-z(ie,6)
      endif
c
      call vcross_normal(w,u,v)
      unx = w(1)
      uny = w(2)
      unz = w(3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine vcross_normal(u,v,w)
C
C     Compute a Cartesian vector cross product.
C
      real u(3),v(3),w(3)
c
      u(1) = v(2)*w(3) - v(3)*w(2)
      u(2) = v(3)*w(1) - v(1)*w(3)
      u(3) = v(1)*w(2) - v(2)*w(1)
c
c     Normalize
c
      r2 = u(1)*u(1) + u(2)*u(2) + u(3)*u(3)
      if (r2.gt.0.) then
         r2 = sqrt(r2)
         u(1) = u(1)/r2
         u(2) = u(2)/r2
         u(3) = u(3)/r2
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_cent
      include 'basics.inc'
c
c     Dump element centroids, w/ element numbers for positional sorting
c                                   pff 4/8/99
c
      return

      call gencen
c
      open(unit=53,file='elcent.dat')
      write(53,*) nel
      do ie=1,nel
         write(53,53) zcen(ie),ycen(ie),xcen(ie),ie
      enddo
   53 format(3f18.8,i9)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc1
      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for hemisphere/plate problem.  pff 12/22/98
c
c     1st, set periodic bcs in x-direction
c
      bclab(2)='$'
      nsides = 2*ndim
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ') then
            if (ccurve(is,ie).eq.'s') then
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'T  '
               bc( 1,is,ie,2) = 1.
            elseif (ccurve(is,ie).eq.'C'.and.is.le.4) then
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'T  '
               bc( 1,is,ie,2) = 1.
            else
c
c              Straight side --- which way is it pointing?
c
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unx).ge.0.9) then
                  ierr = i_periodic(is,ie,1,bclab)
                  ierr = i_periodic(is,ie,2,bclab)
               elseif (abs(uny).ge.0.9) then
                  cbc  (is,ie,1) = 'SYM'
                  cbc  (is,ie,2) = 'I  '
               elseif (unz.ge.0.9) then
                  cbc  (is,ie,1) = 'SYM'
                  cbc  (is,ie,2) = 'I  '
               elseif (unz.le.-0.9) then
                  cbc  (is,ie,1) = 'W  '
                  cbc  (is,ie,2) = 'T  '
                  bc( 1,is,ie,2) = 0.
               endif
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc2
      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for AV graft problem.  pff 6/23/01
c
      ntot = nx*ny*nz*nel
      nsides = ndim*2
      ncrnrs = ndim**2
c
      zmax = z(1,1)
      zmin = z(1,1)
      do ic=1,ncrnrs
      do ie=1,nel
         zmax = max(z(ie,ic),zmax)
         zmin = min(z(ie,ic),zmin)
      enddo
      enddo
      dz = 0.05*(zmax-zmin)
      zdmax = zmax - dz
      zdmin = zmin + dz
c
c
      bclab(2)='$'
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ') then
            cbc  (is,ie,1) = 'W  '
            cbc  (is,ie,2) = 'I  '
            call rzero(bc(1,is,ie,1),5)
            call rzero(bc(1,is,ie,2),5)
c
            if (zcen(ie).gt.zdmax) then
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unz).ge.0.9) then
                  cbc  (is,ie,1) = 'v  '
                  cbc  (is,ie,2) = 't  '
               endif
            elseif (zcen(ie).lt.zdmin) then
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unz).ge.0.9) then
                  cbc  (is,ie,1) = 'O  '
                  cbc  (is,ie,2) = 'O  '
               endif
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc
c
c     Base boundary conditions on face number
c
c     grep " 3 " b.3 | wc 560          3360         41440
c     grep " 5 " b.3 | wc 96           576          7104
c     grep " 6 " b.3 | wc 48           288          3552
c

      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for AV graft problem.  pff 6/23/01
c
      ntot = nx*ny*nz*nel
      nsides = ndim*2
      ncrnrs = ndim**2
c
      zmax = z(1,1)
      zmin = z(1,1)
      do ic=1,ncrnrs
      do ie=1,nel
         zmax = max(z(ie,ic),zmax)
         zmin = min(z(ie,ic),zmin)
      enddo
      enddo
      dz = 0.05*(zmax-zmin)
      zdmax = zmax - dz
      zdmin = zmin + dz
c
c
      bclab(2)='$'
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ' .and.
     $       cbc(is,ie,1).ne.'J  ' .and.
     $       cbc(is,ie,1).ne.'SP ') then
            call rzero(bc(1,is,ie,1),5)
            call rzero(bc(1,is,ie,2),5)
            if (is.eq.6) then
               cbc  (is,ie,1) = 'v  '
               cbc  (is,ie,2) = 't  '
            elseif (is.eq.5) then
               cbc  (is,ie,1) = 'O  '
               cbc  (is,ie,2) = 'O  '
            else
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'I  '
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod
c
      tol=1.e-1
      call autoperiod1(nfail,tol)
      write(6,*) 'NFAIL:',nfail,tol
      if (nfail.eq.0) return
      do k=1,10
         tol = .7*tol
         call autoperiod2(nfail,tol)
         write(6,*) 'NFAIL:',nfail,tol
         if (nfail.eq.0) return
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod1(nfail,tol)
      include 'basics.inc'
      character*3 cb,ck
c
c     Set any remaining candidate bcs to P
c
      if (.not.if_lattice) call get_lattice_per_bc

      call mkside   ! moved, 9/1/09 pff
      call gencen

      nsides = 2*ndim
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'   ') then
C           Automatic selection of Periodic bc
            call fndsidb(js,je,is,ie,if,tol)
            if (js.ne.0) then
               CALL LETBC(JS,JE,IF,BCLAB)
               CALL LETBC(IS,IE,IF,BCLAB)
               CBC( IS,IE,IF) = 'P  '
               call rzero(bc(1,is,ie,if),5)
               call rzero(bc(1,js,je,if),5)
               BC(1,IS,IE,IF) = JE
               BC(2,IS,IE,IF) = JS
               CBC( JS,JE,IF) = 'P  '
               BC(1,JS,JE,IF) = IE
               BC(2,JS,JE,IF) = IS
               if (nel.le.1000.or.mod(ie,100).eq.0)
     $         write(6,1) if,ie,is,je,js
    1          format('autoper:',5i8)
            endif
         endif
      enddo
      enddo
      enddo
c
c     consistency check
c
      nfail = 0
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'P  ') then
C           check for consistency
            je = bc(1,is,ie,if)
            js = bc(2,is,ie,if)
            ke = bc(1,js,je,if)
            ks = bc(2,js,je,if)
            ck = cbc(js,je,if)
            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) then
               write(6,*) ie,is,je,js,ke,ks,ck,' fail',nfail
               cbc(is,ie,if)='   '
               nfail = nfail + 1
               write(6,11) (sides(ie,is,k),k=1,3),ie,is,' side i'
               write(6,11) (sides(je,js,k),k=1,3),je,js,' side j'
               write(6,11) (sides(ke,ks,k),k=1,3),ke,ks,' side k'
            endif
  11        format(1p3e14.5,i8,i3,a7)
         endif
      enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod2(nfail,tol)
      include 'basics.inc'
      character*3 cb,ck
c
c     Set any remaining candidate bcs to P
c
      if (.not.if_lattice) call get_lattice_per_bc

      call mkside   ! moved, 9/1/09 pff
      call gencen

      nsides = 2*ndim
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'   '.or.cb.eq.'P  ') then
C           Automatic selection of Periodic bc
            call fndsidb(js,je,is,ie,if,tol)
            if (js.ne.0) then
               CALL LETBC(JS,JE,IF,BCLAB)
               CALL LETBC(IS,IE,IF,BCLAB)
               CBC( IS,IE,IF) = 'P  '
               call rzero(bc(1,is,ie,if),5)
               call rzero(bc(1,js,je,if),5)
               BC(1,IS,IE,IF) = JE
               BC(2,IS,IE,IF) = JS
               CBC( JS,JE,IF) = 'P  '
               BC(1,JS,JE,IF) = IE
               BC(2,JS,JE,IF) = IS
               if (nel.le.1000.or.mod(ie,100).eq.0)
     $         write(6,1) if,ie,is,je,js
    1          format('autoper:',5i8)
            endif
         endif
      enddo
      enddo
      enddo
c
c     consistency check
c
      nfail = 0
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'P  ') then
C           check for consistency
            je = bc(1,is,ie,if)
            js = bc(2,is,ie,if)
            ke = bc(1,js,je,if)
            ks = bc(2,js,je,if)
            ck = cbc(js,je,if)
            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) cbc(is,ie,if)='   '
            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) nfail = nfail + 1
         endif
      enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fndsidb(kside,ke,iside,ie,if,tol)
      include 'basics.inc'
      real vec1(3),vec2(3),vec3(3),vecj(3)
      integer icrn(4,6)
      save    icrn
      data    icrn / 1 , 2 , 6 , 5 
     $             , 2 , 3 , 7 , 6 
     $             , 3 , 4 , 8 , 7 
     $             , 4 , 1 , 5 , 8 
     $             , 1 , 2 , 3 , 4 
     $             , 5 , 6 , 7 , 8 /
C
C     Find a side JSIDE,JEL which corresponds to ISIDE,IE and doesn't
C     have BC's set already
C
      NSIDES=NDIM*2
C
C     Define Normal to this plane (side), and find the element
C     center which also lies on this line, and is farthest from
C     the current element.
c
c     CALL MKSIDE   ! moved, 9/1/09 pff
c     CALL GENCEN
C
C     Find Sides' Normal Vector at midpoint

      ic1 = icrn(1,iside)
      ic2 = icrn(2,iside)
      ic3 = icrn(3,iside)
      ic4 = icrn(4,iside)

c     write(s,9) ie,iside,ic1,ic2,ic3,ic4
c   9 format(' iesc14:',6i5,'$')
c     if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)

      IF (IF3D) THEN
         VEC1(1)=X(IE,IC3)-X(IE,IC1)
         VEC1(2)=Y(IE,IC3)-Y(IE,IC1)
         VEC1(3)=Z(IE,IC3)-Z(IE,IC1)

         VEC2(1)=X(IE,IC4)-X(IE,IC2)
         VEC2(2)=Y(IE,IC4)-Y(IE,IC2)
         VEC2(3)=Z(IE,IC4)-Z(IE,IC2)

         CALL CROSS(VEC3,VEC1,VEC2)
         CALL NORM3D(VEC3)
      ELSE

         VEC1(1)=X(IE,IC2)-X(IE,IC1)
         VEC1(2)=Y(IE,IC2)-Y(IE,IC1)
C
         VEC3(1)= VEC1(2)
         VEC3(2)=-VEC1(1)
         VEC3(3)=0.0
         CALL NORM3D(VEC3)
      ENDIF
c
      if (if_lattice) then          ! Make vector point in lattice direction
         d1 = dotprod(vec3,ulat1)
         d2 = dotprod(vec3,ulat2)
         d3 = dotprod(vec3,ulat3)
         if (d1.le.d2.and.d1.le.d3) then
            call copy(vec3,ulat1,3)
         elseif (d2.le.d3) then
            call copy(vec3,ulat2,3)
         else
            call copy(vec3,ulat3,3)
         endif
      endif
c
c     write(s,19) vec3(1),vec3(2),vec3(3)
c  19 format(' vec3:',3e13.5,'$')
c     if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
c
c     EPS2=(EPSM*RCEN(IE))**2
      EPSM=tol
      EPS2=EPSM**2
      D2MN=10*EPS2
c
      KE=0
      KSIDE=0
      DSTMAX=0.0
c
      DO 100 JE=1,NEL
      DO 100 JSIDE=1,NSIDES
C
C        Don't find yourself.
         IF (JE.EQ.IE.AND.JSIDE.EQ.ISIDE) GOTO 100
C
C        Don't find a bc which is already defined. (to be 'E')
         IF (CBC(JSIDE,JE,IF).EQ.'E  ') GOTO 100
c
c        Check if dot product of surface normals are O(1)
c

         ic1 = icrn(1,jside)
         ic2 = icrn(2,jside)
         ic3 = icrn(3,jside)
         ic4 = icrn(4,jside)

         vec1(1)=x(je,ic3)-x(je,ic1)
         vec1(2)=y(je,ic3)-y(je,ic1)
         vec1(3)=z(je,ic3)-z(je,ic1)

         vec2(1)=x(je,ic4)-x(je,ic2)
         vec2(2)=y(je,ic4)-y(je,ic2)
         vec2(3)=z(je,ic4)-z(je,ic2)

         call cross (vecj,vec1,vec2)
         call norm3d(vecj)
         dotnorm = dotprod(vec3,vecj)
         if (abs(dotnorm).lt.0.9) goto 100
C
C        OK, is the center of this element side close to the line?
C
         VEC1(1)=SIDES(IE,ISIDE,1)-SIDES(JE,JSIDE,1)
         VEC1(2)=SIDES(IE,ISIDE,2)-SIDES(JE,JSIDE,2)
         VEC1(3)=SIDES(IE,ISIDE,3)-SIDES(JE,JSIDE,3)
         DISTP=DOTPROD(VEC1,VEC1)
         dis13=(dotprod(vec1,vec3))**2
         dist2=distp-dis13
         if (distp.gt.0) DIST2=DIST2/DISTP
c        write(6,11) dist2,distp,eps2,dstmax,je,jside
c  11    format(1p4e12.4,2i6,' dist2')
c        IF (DIST2.LE.EPS2.AND.DISTP.GT.DSTMAX) THEN
         IF (DIST2.LE.D2MN) THEN
            KE    = JE
            KSIDE = JSIDE
            DSTMAX= DISTP
c     if ((69.le.ie.and.ie.le.76).or.(9114.le.ie.and.ie.le.9120)) then
c     write(6,12) distp,dist2,(sides(ie,iside,k),k=1,3),ie,iside,' si'
c     write(6,12) dis13,d2mn ,(sides(je,jside,k),k=1,3),je,jside,' sj'
c     write(6,12) dis13,d2mn ,(vec1(k),k=1,3),je,jside,'vc1'
c     write(6,12) dis13,dotnorm ,(vecj(k),k=1,3),je,jside,'vcj'
c  12 format(1p5e14.5,i8,i3,1x,a3)
c     endif

            D2MN  = DIST2

         ENDIF
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
