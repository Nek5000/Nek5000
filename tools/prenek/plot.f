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
C     Beginning of PLOT subroutines, common to preprocessor and postprocessor
c-----------------------------------------------------------------------
      SUBROUTINE DATA
C     This subroutine replaces data statements
#     include "basics.inc"
      COMMON /TTSCAL/ TX(16)
C
      OPEN(UNIT=13,STATUS='SCRATCH')
C
C     Initialize Character Strings to blank.
C
      CALL BLANK(ARG,50)
      LINE = ' '
      SESION = ' '
      OLDVER = ' '
      RSTV = ' '
      RSTT = ' '
      X13 = ' '
      PLANE2 = ' '
      VCOND = ' '
      QUANTY = ' '
      COMPON = ' '
      SCALPT = ' '
      VECTPT = ' '
      XYATTR = ' '
      SUATTR = ' '
      PLFORM = ' '
      CBC1 = ' '
      DARROW = ' '
      TURBMOD = ' '
      CRVTYP = ' '
      ANS = ' '
      button = ' '
      CHOICE = ' '
C
C     Initialize transformation matrix
C     Default (pre- or post-) is 2D. For 2-d, use overhead perspective
      ISOBJ=1
      THETA =  90.
      PHI   = -90.
      IROT  = 0
      PI=4.0*ATAN(1.0)
      THETAr=THETA*PI/180.0
      PHIr  =PHI  *PI/180.0
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
C     Set clipping window

      WFR= 1.3
      IFFULL=.FALSE.
      IF (IFPOST) IFFULL=.TRUE.
      WT = YPHY(.995)
      WB = YPHY(.005)
      WL = XPHY(.005)
      WR = XPHY(.995)
      IF (IFFULL) WR = XPHY(WFR)
C
C     Assign EB's numbering scheme to PF's scheme.
C
      EFACE(1)=4
      EFACE(2)=2
      EFACE(3)=1
      EFACE(4)=3
      EFACE(5)=5
      EFACE(6)=6
C
C     Assign inverse of EB's numbering scheme to PF's scheme.
C
      EFACE1(1)=3
      EFACE1(2)=2
      EFACE1(3)=4
      EFACE1(4)=1
      EFACE1(5)=5
      EFACE1(6)=6
C
C     Save image data initialization
C
      IDRAWD=0
      IDRCNT=0
C
      N=6*NELM*(MAXFLD+1)*3
      CALL BLANK(CBC,N)
C     Equation types
      NKTONV=2
      VNEKTON=2.61
      WRITE(S,'(A16,F5.1)') ' NEKTON Version ',VNEKTON
      CALL PRS(S//'$')
      LOCLIN=1
      IFGRID = .TRUE.
      IFGRDC = .TRUE.
      IFGRDH = .FALSE.
      IFGRDP = .FALSE.
      IFSTRS = .FALSE.
      IFHEAT = .FALSE.
      IFTRAN = .TRUE.
      IFCHAR = .FALSE.
      IFNAV  = .TRUE.
      IFFLOW = .TRUE.
      IFAXIS = .FALSE.
      IF3D   = .FALSE.
      IFCEIL = .FALSE.
      IFERROR= .FALSE.
      IFLEARN= .FALSE.
      IFDEMO = .FALSE.
      IFDRAX = .TRUE.
      IFVMS  = .FALSE.
      ifrevbk  = .FALSE.
      if_lattice  = .FALSE.

      ifmid  = .false.   ! Allow 'm' ccurve
      ifcstd = .true.    ! Use std ccurve options


      DO 2 I=1,11
         IFADVC(I) = .FALSE.
         IFTMSH(I) = .TRUE.
2     CONTINUE
      IFTMSH(1)=.FALSE.
C     For the default heat transfer only, default IFADVC is false.
C      IFADVC(2)=.TRUE.
      NDIM=3
C
C     DEFAULT OCODES
      IFMOVB=.FALSE.
      IFXYH=.FALSE.
      IFXYO=.FALSE.
      IFVO =.TRUE.
      IFTO =.TRUE.
      IFPO =.TRUE.
      IFTGO=.FALSE.
C
      TEXTSW(1,1)='NONE'
      TEXTSW(1,2)='TURBMOD'
      IPSCO=1
      DO 18 IPSCAL=1,9
         IFPSCO(IPSCAL)=.TRUE.
         WRITE(PSNAME(IPSCAL),'(''PS'',1X,I1,1X)')IPSCAL
18    CONTINUE
      XFAC=1
      YFAC=1
c      IF(NKTONV.EQ.1)THEN
         OCODE(1) = 'XX'
         OCODE(2) = 'YY'
         OCODE(3) = 'UU'
         OCODE(4) = 'VV'
         OCODE(5) = 'TT'
         OCODE(6) = 'PP'
         OCODE(7) = 'TX'
         OCODE(8) = 'TY'
c      ENDIF
      DO 20 I=11,20
         OCODE(I)= '  '
20    CONTINUE
C
C     Objects
      NVOBJS = 0
      NSOBJS = 0
      NEOBJS = 0
      NPOBJS = 0

      MAXLET=0
      ILETAP=MAXLET+96-32

      DO 30 I=1,500
         PARAM(I)=0.0
30    CONTINUE
      DO 40 IEL=1,NELM
         DO 40 IEDGE=1,12
            CCURVE(IEDGE,IEL)=' '
            DO 40 I=1,6
               CURVE(I,IEDGE,IEL)=0.0
40    CONTINUE
      DO 34 IF=1,10
        DO 34 IG=-5,10
           MATYPE(IG,IF)=0
           DO 34 IPROP=1,3
            VPROP(IG,IF,IPROP)='C'
            CPROP(IG,IF,IPROP)=0.0
34    CONTINUE
      IGRP=0
C
      RSTV =' '
      RSTT =' '
      IRSTV=-1
      IRSTT=-1
      IRSTIM=-1
      RSTIM=0.0
      DERIV='NO'
C
      ISEGDUM=30
C
C     TEMPERATURE SCALE
      TX(1)=-0.0004
      TX(2)= .071
      TX(3)=0.143
      TX(4)=0.214
      TX(5)=0.286
      TX(6)=0.357
      TX(7)=0.429
      TX(8)=0.5
      TX(9)=0.571
      TX(10)=0.642
      TX(11)=0.710
      TX(12)=0.786
      TX(13)=0.857
      TX(14)=0.929
      TX(15)=1.0004
      TX(16)=1.071

C     Stuff to be written out in fortran statements
      DO 35 I=1,15
         INITC(I)='C Default'
35    CONTINUE
      DO 39 I=1,15
         DRIVC(I)='C'
39    CONTINUE

      call setgraph(ifgraf)

C     Data that used to be initialized in postnek and postnek2
      SARROW=1.0
      DARROW='HIGH'
      HFISH=0.3
      NXGRID=9
      NXBAND=10
      NXCONT=15
      cont_lev=2.0
      IPLTYP=2
      NPTPRO=100
      NCSEGS=10
      VELTYP=1.0
      LOCLIN=1
C
      IDUMP=0

      IPLANE=1
      JPLANE=1
      KPLANE=1
      ILEVEL=1
      PLANE='Z-PLANE'

      DO 45 I=1,NELM
         IFPLOT(I)=.true.
45    CONTINUE
      IINTEG=1
      NINTEG=0
C
C
      RETURN
      END
C     These 4 functions transform coordinates from screen to world and back
c-----------------------------------------------------------------------
      FUNCTION XSCR(XPHYS)
#     include "basics.inc"
      XSCR = (XPHYS-XZERO) /XFAC
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION YSCR(YPHYS)
#     include "basics.inc"
      YSCR = (YPHYS-YZERO) /YFAC
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION XPHY(XSCREN  )
#     include "basics.inc"
      XPHY = XZERO + XSCREN  * XFAC
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION YPHY(YSCREN  )
#     include "basics.inc"
      YPHY = YZERO + YSCREN  * YFAC
      RETURN
      END
      subroutine ginGRD(X1)
C     Sets Coarseness of Grid on tablet and screen. 0 < x1 < 1
      INCLUDE 'devices.inc'
      IF(X1.LT.0.0 .OR. X1.GT.1.0) THEN
         CALL PRS('Error: Grid must be between 0 and 1 (fraction of$')
         CALL PRS
     $   ('Screen height).  Screen gridding couldnt be changed$')
         RETURN
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DRWBOX(X1,Y1,X2,Y2,IC)
c     draw a box
      CALL COLOR(IC)
      CALL MOVEC(X1,Y1)
      CALL DRAWC(X2,Y1)
      CALL DRAWC(X2,Y2)
      CALL DRAWC(X1,Y2)
      CALL DRAWC(X1,Y1)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE XDOT(XX,YY)
#     include "basics.inc"
c     draws cross around x,y
      CALL MOVEC(XX-.005*XFAC,YY-.005*YFAC)
      CALL DRAWC(XX+.005*XFAC,YY+.005*YFAC)
      CALL MOVEC(XX+.005*XFAC,YY-.005*YFAC)
      CALL DRAWC(XX-.005*XFAC,YY+.005*YFAC)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DDUMMY(IEDGE,IEL)
C     Used to help evade a compiler bug in the Titan
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DASHSC(X1,Y1,X2,Y2,N)
#     include "basics.inc"
c     draws a dashed line from pt 1 to pt2, N dashes.
      IF (N.LE.0) RETURN
      DELTX=(X2-X1)/(2*N-1)
      DELTY=(Y2-Y1)/(2*N-1)
      DELTX2=2*DELTX
      DELTY2=2*DELTY
      XX1=X1
      YY1=Y1
      XX2=XX1+DELTX
      YY2=YY1+DELTY
      DO 10 I=1,N
         CALL MOVESC(XX1,YY1)
         CALL DRAWSC(XX2,YY2)
         XX1=XX1+DELTX2
         YY1=YY1+DELTY2
         XX2=XX2+DELTX2
         YY2=YY2+DELTY2
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ARROW(X,Y,DXBAR,DYBAR,SARROW)
      COMMON/SCALE/XFAC,YFAC,XZERO,YZERO
      DX=DXBAR*SARROW
      DY=DYBAR*SARROW
      IF(DX.EQ.0.0 .AND. DY.EQ.0.0) RETURN
      CALL MOVEC(X,Y)
      CALL DRAWC(X+DX,Y+DY)
      CALL DRAWC(X+0.85*DX-.1*DY,Y+0.85*DY+.1*DX*yfac/xfac)
      CALL MOVEC(X+DX,Y+DY)
      CALL DRAWC(X+0.85*DX+0.1*DY,Y+0.85*DY-0.1*DX*yfac/xfac)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ARROW4(X,Y,Z,DX,DY,DZ)
      COMMON/SCALE/XFAC,YFAC,XZERO,YZERO
C     Plot a X pointed arrow
      CALL COLOR(12)
      IF(DX.EQ.0.0 .AND. DY.EQ.0.0 .AND.DZ.EQ.0.0) RETURN
      CALL PENW(1)
      CALL COLOR(1)
      CALL MOVE3(X,Y,Z)
      X1=X+DX
      Y1=Y+DY
      Z1=Z+DZ
      CALL DRAW3(X1,Y1,Z1)
C     Wierd arrowheads
      CALL DRAW3(X+0.85*DX-.1*DY ,Y+0.85*DY+0.1*DX*yfac/xfac,Z+DZ*0.85)
      CALL MOVE3(X+DX,Y+DY,Z+DZ)
      CALL DRAW3(X+0.85*DX+0.1*DY,Y+0.85*DY-0.1*DX*yfac/xfac,Z+DZ*0.85)
      IF(DX.EQ.0.0 .AND. DY.EQ.0.0) THEN
C        Kludge to get arrowhead on vertical arrow
         CALL MOVE3(X,Y,Z+DZ)
         CALL DRAW3(X,Y+0.1*Dz,Z+DZ*0.85)
         CALL MOVE3(X,Y,Z+DZ)
         CALL DRAW3(X,Y-0.1*Dz,Z+DZ*0.85)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DRAW3(X3,Y3,Z3)
      X=XISO(X3,Y3,Z3)
      Y=YISO(X3,Y3,Z3)
      CALL DRAWC(X,Y)
c     CALL DRAWSC(XSCR(X),YSCR(Y))
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE MOVE3(X3,Y3,Z3)
      X=XISO(X3,Y3,Z3)
      Y=YISO(X3,Y3,Z3)
      CALL MOVEC(X,Y)
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION XISO(X3,Y3,Z3)
#     include "basics.inc"
C     Deals in PHYSICAL (World) Coordinates.
      REAL X3D(3)
      IF (IROT.EQ.0) THEN
         X3D(1)=X3
         X3D(2)=Y3
         X3D(3)=Z3
      ELSE
         X3D(1)=Z3
         X3D(2)=X3
         X3D(3)=Y3
      ENDIF
      XISO=DOTPROD(X3D,XHOBS)
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION YISO(X3,Y3,Z3)
#     include "basics.inc"
      REAL X3D(3)
      IF (IROT.EQ.0) THEN
         X3D(1)=X3
         X3D(2)=Y3
         X3D(3)=Z3
      ELSE
         X3D(1)=Z3
         X3D(2)=X3
         X3D(3)=Y3
      ENDIF
      YISO=DOTPROD(X3D,YHOBS)
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION ZISO(X3,Y3,Z3)
#     include "basics.inc"
      REAL X3D(3)
      IF (IROT.EQ.0) THEN
         X3D(1)=X3
         X3D(2)=Y3
         X3D(3)=Z3
      ELSE
         X3D(1)=Z3
         X3D(2)=X3
         X3D(3)=Y3
      ENDIF
      ZISO=DOTPROD(X3D,VHOBS)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNISO(X3D,Y3D,Z3D,XI,YI,ZI)
#     include "basics.inc"
C     Given XI and YI (as mouse inputs, and knowing Z3D (from the current plane)
C     Calculate X3D and Y3D (and incidentally, ZISO)
      RATIO=-YHOBS(1)/XHOBS(1)
      Y3D=(RATIO*XI + YI - RATIO*Z3D*XHOBS(3) - Z3D*YHOBS(3))
     $  / (YHOBS(2)+RATIO*XHOBS(2))
      X3D=(XI-Y3D*XHOBS(2)-Z3D*XHOBS(3))/XHOBS(1)
      ZI = X3D*VHOBS(1) + Y3D*VHOBS(2) + Z3D*VHOBS(3)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNISO2(X3D,Y3D,Z3D,XI,YI,ZI)
#     include "basics.inc"
C     Given XI and YI (as mouse inputs, and knowing y3D (from the current plane)
C     Calculate X3D and Z3D (and incidentally, ZISO)
      RATIO=-YHOBS(1)/XHOBS(1)
      Z3D=(RATIO*XI + YI - Y3D*(RATIO*XHOBS(2) + YHOBS(2)) )
     $  / (YHOBS(3)+RATIO*XHOBS(3))
      X3D=(XI-Y3D*XHOBS(2)-Z3D*XHOBS(3))/XHOBS(1)
      ZI = X3D*VHOBS(1) + Y3D*VHOBS(2) + Z3D*VHOBS(3)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE UNISO3(X3D,Y3D,Z3D,XI,YI,ZI)
#     include "basics.inc"
C     Given XI and YI (as mouse inputs, and knowing X3D (from the current plane)
C     Calculate Y3D and Z3D (and incidentally, ZISO)
      RATIO=-YHOBS(2)/XHOBS(2)
      Z3D=(RATIO*XI + YI - X3D*(RATIO*XHOBS(1) + YHOBS(1)) )
     $  / (YHOBS(3)+RATIO*XHOBS(3))
      Y3D=(XI-X3D*XHOBS(1)-Z3D*XHOBS(3))/XHOBS(2)
      ZI = X3D*VHOBS(1) + Y3D*VHOBS(2) + Z3D*VHOBS(3)
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION DOTPROD(X,Y)
      DIMENSION X(3),Y(3)
      DOTPROD=0.0
      DO 1 I=1,3
         DOTPROD=DOTPROD+X(I)*Y(I)
1     CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE G3WRIT(X,Y,Z,SIZE,STRING)
      CHARACTER STRING(100)
      CALL GWRITE(XISO(X,Y,Z),YISO(X,Y,Z),SIZE,STRING)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DRELEV(IDRAW,IFADE,FLAG)
C     Color-coded floors (with matching sides of isometric drawing??!!
C     # of elements on floor next to elevator (like the red sox??!!)
#     include "basics.inc"
      CHARACTER*5 FLAG

      IF(.NOT. IF3D) RETURN

C     2nd arg: 0 for no fade; negative for bottom only (for floor mods)
C     Draw Elevator
      ELEC =0.95
      ELEW =0.02
      ELEB =0.1
      ELEDY=0.1
      WLINE=1/600.
      DO 1 I=1,2
C        First fade old Elevator to gray (i=1), then draw(i=2)
         IF(I.EQ.1)ITYPE=IFADE
         IF(I.EQ.2)ITYPE=IDRAW
         IF(I.EQ.1)CALL COLOR(15)
         IF(I.EQ.2)CALL COLOR(1)
         IF(I.EQ.2.AND.FLAG.EQ.'BLINK')CALL COLOR(5)
         IF(ITYPE.NE.0)THEN
            IF(ITYPE.GT.0)THEN
               DO 31 ILINE=-2,2
                  CALL MOVESC(ELEC+WLINE*ILINE,ELEB+ELEDY*(ITYPE-1))
                  CALL DRAWSC(ELEC+WLINE*ILINE,ELEB+ELEDY*(ITYPE  ))
31             CONTINUE
            ENDIF
C           Do this part only if ITYPE was < 0
C           But Skip this part for blinkers on sides 1-4
            IF(ITYPE.LT.0 .OR. FLAG.NE.'BLINK')THEN
               IF(ITYPE.LT.0) ITYPE= -ITYPE
               IF(FLAG.EQ.'POST ') CALL COLOR(15)
               CALL MOVESC(ELEC-ELEW,ELEB+ELEDY*(ITYPE-1))
               CALL DRAWSC(ELEC+ELEW,ELEB+ELEDY*(ITYPE-1))
            ENDIF
         ENDIF
1     CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SORTZ
#     include "basics.inc"
C
      COMMON /CTMPS/  Wk(NELM,6) , IND(NELM,6)
C     Fcorns has the corners, in ccw order, corresponding to each of 6 faces.
      INTEGER FCORNS (4,6)
      save    FCORNS
      DATA FCORNS / 1,2,6,5,
     $              2,3,7,6,
     $              3,4,8,7,
     $              4,1,5,8,
     $              1,2,3,4,
     $              8,7,6,5 /
C
      NUMBER=0
      DO 250 IFACE=1,6
         DO 250 IEL=1,NEL
             XFA=0.0
             YFA=0.0
             ZFA=0.0
C            FIND X,Y,Z at center of each face
             DO 20 ICORN=1,4
                xfa=xfa+x(fcorns(icorn,iface),iel)
                yfa=yfa+y(fcorns(icorn,iface),iel)
                zfa=zfa+z(fcorns(icorn,iface),iel)
20           CONTINUE
             NUMBER=NUMBER+1
             xfa=0.25*xfa
             yfa=0.25*yfa
             zfa=0.25*zfa
             ZDEPTH(NUMBER,1)=ZISO(XFA,YFA,ZFA)
             IZDPTH(NUMBER,1)=IEL
             IZDPTH(NUMBER,2)=IFACE
C            ZDEPTH(IEL,IFACE)=ZISO(XFA,YFA,ZFA)
C            ZLEFT(IEL,IFACE)=.TRUE.
250   continue
C     Sort so that they are output in the right order
      NFCTOT=6*NEL
      CALL SORT(ZDEPTH,IND,NFCTOT)
      CALL ISWAP(IZDPTH(1,1),Wk,IND,NFCTOT)
      CALL ISWAP(IZDPTH(1,2),Wk,IND,NFCTOT)
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SETVUE
#     include "basics.inc"
C
      CALL PRS('Enter viewing direction from which $')
      CALL PRS('the observer''s eye sees the data$')
10    CALL PRS('Enter Altitude above x-y plane (in degrees):$')
      CALL PRSRS('Default <cr> =$',THETA,' degrees$')
      CALL RES(LINE,70)
      IF(IFLEARN)WRITE(3,'(A70)')LINE
      IF(LINE.EQ.' ')THEN
      ELSE
         CALL READER(THETA,IERR)
         IF(IERR.NE.0) GO TO 20
c         IF(THETA.LT.-90.0 .OR. THETA .GT.90.0) THEN
c            PRINT*,' -90.0 < THETA < 90.0   Try again.'
c            GO TO 10
c         ENDIF
      ENDIF
20    CALL PRS('Enter Azumuth (in degrees):$')
      CALL PRS('0 degrees= from above x-axis; 90=from above y axis$')
      CALL PRSRS('Default <cr> =$',PHI,' degrees$')
      CALL RES(LINE,70)
      IF(IFLEARN)WRITE(3,'(A70)')LINE
      IF(LINE.EQ.' ')THEN
      ELSE
         CALL READER(PHI  ,IERR)
         IF(IERR.NE.0) GO TO 20
      ENDIF
C     Now set viewing transformation
      THETAr=THETA*PI/180.0
      PHIr  =PHI  *PI/180.0
C
C       Direction toward eye (normal vector)
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
C     Rescale screen-world factors so it fits on screen
      IF (.NOT.IFZOOM)CALL RESCAL
      RETURN
      END
C     Beginning of TEKPLOT preprocessor and postprocessor subroutines *******
c-----------------------------------------------------------------------
      SUBROUTINE MENU(XMOUSE,YMOUSE,BUTTON,MTITLE)
C       Segments:  ; Cover:13;  Menu: 14; Color bar:15
#     include "basics.inc"
      CHARACTER IT,CH,MTITLE*(*)
      CHARACTER*26 OLDITM(20)
      LOGICAL NEWMNU
      SAVE IFIRST
      DATA IFIRST /0/
c
      LOGICAL ifgrdt
c
c     Turn off ifgrid!  (pff 4/25/97)
c
      nc13=max(nchoic+1,13)
      nc14=max(nchoic+2,14)
      nc15=max(nchoic+2,15)
      ifgrdt = ifgrid
      if (ifpost) ifgrid = .false.
      if (mtitle.eq.'CENTRAL') ifgrid = .false.
C
        IF (IFGRAF) THEN
C
5          CONTINUE
           IF(IFIRST.NE.1)THEN
              IFIRST=1
C             Draw dummy color bar on 15
              CALL DRCOVR(nc15)
C             Draw cover first time around.
C             DRCOVR physically draws cover on nc13
              CALL DRCOVR(nc13)
           ENDIF
           IF(IFNOSEG)THEN
C             DRAW MENU (Always)
              CALL DRCOVR(nc13)
              NEWMNU=.TRUE.
           ELSE
C             DRAW MENU (ONLY if NEW!!!)
              NEWMNU=.FALSE.
              DO 20 I=1,20
                 IF(OLDITM(I).NE.ITEM(I)) NEWMNU=.TRUE.
20            CONTINUE
           ENDIF
C          !!?? Patch - REDRAW ALTER MENU ALWAYS;
C          There is scuzzy numbers on screen
           IF(ITEM(2).EQ.'ALTER NUMERICAL PARAMETERS') NEWMNU=.TRUE.
           IF(NEWMNU)THEN
              DO 25 I=1,20
25               OLDITM(I)=ITEM(I)
C             CALL GSWRIT(MTITLE)
C             Drmenu always does the full draw (doesn't just turn on and off)
              CALL DRMENU(mtitle)
           ELSE
C             Turn on menu
              IF(MTITLE.NE.'NOCOVER')THEN
                 CALL SGVIS(14,0)
                 CALL SGVIS(14,1)
              ENDIF
           ENDIF
130        CONTINUE
           CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
           XSCM=XSCR(XMOUSE)
           YSCM=YSCR(YMOUSE)
C          Special kludge for building
           IF(XSCR(XMOUSE).LT.1.0
     $        .AND. ( ITEM(1).EQ.'ADD    ELEMENT'
     $         .OR.   ITEM(1).EQ.'END    ELEMENTS')) THEN
                      ICHOIC=1
                      CHOICE='ADD    ELEMENT'
                      CH=CHOICE
                      CALL PRS(' $')
                      CALL PRS(' $')
                      CALL PRS(CHOICE//'$')
                      ifgrid = ifgrdt
                      RETURN
           ELSE IF(XSCR(XMOUSE).GE.XLMEN
     $         .AND. XSCR(XMOUSE).LE.XRMEN ) THEN
               DO 135 IBOX = 1,nCHOIC
                  IF(YSCR(YMOUSE).GE.YBS(IBOX) .AND. YSCR(YMOUSE)
     $                .LE.YTS(IBOX))THEN
                      ICHOIC=IBOX
                      CHOICE=ITEM(ICHOIC)
                      CH=CHOICE
                      IF(MTITLE.NE.'NOECHO')THEN
                         CALL PRS(' $')
                         CALL PRS(' '//CHOICE//'$')
                      ENDIF
C                     Need to turn off centered menu immediately for grid.
                      IF(IFCEN)CALL SGVIS(nc14,0)
C                     Turn on cover
                      IF(MTITLE.NE.'NOCOVER'.AND. .NOT.IFCEN)THEN
                         CALL SGVIS(nc13,0)
                         CALL SGVIS(nc13,1)
                      ENDIF
                      ifgrid = ifgrdt
                      RETURN
                  ENDIF
135            CONTINUE
          ENDIF
          IF (BUTTON.EQ.'RIGHT') then
             CHOICE=' '
             ifgrid = ifgrdt
             return
          endif
C         Here we have a point outside the menu area.  Redraw menu.
          CALL DRMENU(mtitle)
          CALL BEEP
          CALL PRS('Choose one of the items in the menu$')
          GO TO 130
C
      ELSE  ! ifgraf turned off
C
C        NON-GRAPHICS DEVICE INPUT
C
 1130    do 1135 ibox = 1,nchoic
c           call pris (ibox,item(ibox)//'$')
            call pris4(ibox,item(ibox)//'$')
 1135    continue
         call prs('Choose item:$')
         call rei(ichoice)
         if (ichoice.lt.1.or.ichoice.gt.nchoic) then
            call prs('Choose one of the items in the menu$')
            goto 1130
         else
            choice = item(ichoice)
            ch     = choice
            button = 'LEFT'
            xmouse = 1.04
            ymouse = 0.00
         endif
       endif

       ifgrid = ifgrdt
       return

       end
c-----------------------------------------------------------------------
       SUBROUTINE DRMENU(mtitle)
#     include "basics.inc"
       CHARACTER*(*) MTITLE
       CHARACTER*21 TITLE
       CHARACTER*101 BOXLAB
      SAVE IFIRST
      DATA IFIRST /0/
C      Draws menu onto segment 14
       nc13=max(nchoic+1,13)
       nc14=max(nchoic+2,14)
       nc15=max(nchoic+2,15)
       denom=25.
       IF(IFCEN)THEN
          XLMEN = 0.25
          XRMEN = 0.75
          denom = 20.
       ELSE
          XLMEN = 1.003
          XRMEN = 1.297
          CALL SGVIS(nc13,0)
          CALL SGVIS(nc13,1)
       ENDIF
       CALL CLSSEG
       if(ifirst.eq.1)CALL DELSEG(nc14)
       CALL OPNSEG(nc14)
       ifirst=1
       IF(.NOT.IFCEN)NBOX=NCHOIC+1
       IF(     IFCEN)NBOX=NCHOIC
       DO 10 IBOX = 1,NBOX
           IF(IFCEN)THEN
              YBS(IBOX) = 0.89-(IBOX  )/denom
              YTS(IBOX) = 0.89-(IBOX-1)/denom
           ELSE
              YBS(IBOX) = 0.69-(IBOX  )/denom
              YTS(IBOX) = 0.69-(IBOX-1)/denom
           ENDIF
C          Draw Red box below last item so that black spot will not appear
           IF(IBOX.EQ.NCHOIC+1)YBS(IBOX) = 0.69-(nc13-1)/denom
cpff 2/99  IF(IBOX.EQ.NCHOIC+1.AND.IFPOST)YBS(IBOX) = 0.11
           IF(IBOX.EQ.NCHOIC+1.AND.IFPOST)YBS(IBOX) = 0.03
C
           CALL COLOR(1)
C          BLUE
           CALL FILLP(-4)
           IF(IBOX.EQ.NCHOIC+1)THEN
              IF(IFPOST) THEN
                 CALL FILLP(-13)
                 CALL COLOR( 13)
              ELSE
                 CALL FILLP(-2)
                 CALL color(2)
              ENDIF
           ENDIF
           CALL BEGINB (XPHY(XLMEN),YPHY(YBS(IBOX)))
           CALL MOVESC (XLMEN,YBS(IBOX))
           CALL DRAWSC (XRMEN,YBS(IBOX))
           CALL DRAWSC (XRMEN,YTS(IBOX))
           IF(IBOX.EQ.NCHOIC+1)CALL COLOR(1)
C          make top line white
           CALL DRAWSC (XLMEN,YTS(IBOX))
           IF(.NOT. IFPOST .AND. IBOX.EQ.NCHOIC+1)CALL COLOR(2)
           IF(      IFPOST .AND. IBOX.EQ.NCHOIC+1)CALL COLOR(13)
           CALL DRAWSC (XLMEN,YBS(IBOX))
           CALL ENDP
           CALL COLOR(1)
C
           XCH=XLMEN+.0122
           YCH=(YBS(IBOX)+YTS(IBOX))/2 - 0.005
           BOXLAB=ITEM(IBOX)
           BOXLAB(26:26)='$'
           if (ifgraf) then
             IF(IBOX.LE.NCHOIC)CALL GSWRIT(XCH,YCH,1.0,BOXLAB)
             IF(IBOX.LE.NCHOIC)write(6,'(i2,2x,A26)') ibox,item(ibox)
           endif
10    CONTINUE
C     Draw menu title
      CALL FILLP(-0)
      CALL BEGINP (XPHY(.01),YPHY(.95))
      CALL MOVESC (.01,.95)
      CALL MOVESC (.30,.95)
      CALL MOVESC (.30,.999)
      CALL MOVESC (.01,.999)
      CALL MOVESC (.01,.95)
      CALL ENDP
      TITLE=MTITLE
      IF(MTITLE.EQ.'NOCOVER')THEN
C        Either build or bound
         IF(ITEM(1).EQ.'ADD    ELEMENT' .OR. 
     $      ITEM(1).EQ.'MODIFY ELEMENT' )THEN
            TITLE='BUILD'
         ELSE
            TITLE='BOUNDARY CONDITIONS'
         ENDIF
      ENDIF
      IF(MTITLE.EQ.'NOECHO')THEN
         TITLE='SET EQUATION TYPE'
      ENDIF
      TITLE(21:21)='$'
      CALL GSWRIT(.02,.965,1.0,TITLE)
C
      CALL CLSSEG
      CALL OPNSEG(10)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DRCOVR(ICOVER)
#     include "basics.inc"
C
      CALL CLSSEG
c      CALL DELSEG(ICOVER)
      CALL OPNSEG(ICOVER)
      IF(IFPOST) THEN
          CALL FILLP(-0)
      ELSE
          CALL FILLP(-2)
          call color(3)
      ENDIF
C      Draws cover onto segment 113
      IF(IFCEN.AND.IFNOSEG)THEN
          XLMEN = 0.25
          XRMEN = 0.76
          ytop=0.9
          ybot=0.05
          CALL FILLP(-0)
      ELSE
          XLMEN = 1.003
          XRMEN = 1.297
          ytop=0.69
cpff 2/99 ybot=0.69-12./20.
          ybot=0.69-12./25.
      ENDIF
c
      CALL BEGINP (XPHY(XLMEN),YPHY(YBot))
      CALL MOVESC (XLMEN,YBot)
      CALL MOVESC (XRMEN,YBot)
      CALL MOVESC (XRMEN,YTop)
      CALL MOVESC (XLMEN,YTop)
      CALL MOVESC (XLMEN,YBot)
      CALL ENDP
C     Now draw cover for menu title
      CALL FILLP(-0)
      CALL BEGINP (XPHY(.01),YPHY(.95))
      CALL MOVESC (.01,.95)
      CALL MOVESC (.30,.95)
      CALL MOVESC (.30,.999)
      CALL MOVESC (.01,.999)
      CALL MOVESC (.01,.95)
      CALL ENDP
c
      CALL CLSSEG
      CALL OPNSEG(10)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE DRBAR(TX,NTX)
c     Draws Rectangles for temperature scale onto segment 15
#     include "basics.inc"
      CHARACTER*15 C
      REAL TX(16)
C     Draw colorbar first time around
      CALL CLSSEG
      CALL DELSEG(15)
      CALL OPNSEG(15)
      Do 20 icolor=14,1,-1
         XL=1.2
         XR=XL  +1.0/27.0
         YB=.12+ 1.0/27.0 * ICOLOR
         YT=YB + 1.0/27.0
         XC=(XL+XR)/2.
         YC=(YT+YB)/2.
         CALL FILLP(-(ICOLOR+1))
         CALL BEGINB(XPHY( XL),YPHY( YB))
         CALL MOVE  (XPHY( XL),YPHY( YB))
         CALL MOVE  (XPHY( XL),YPHY( YT))
         CALL MOVE  (XPHY( XR),YPHY( YT))
         CALL MOVE  (XPHY( XR),YPHY( YB))
         CALL MOVE  (XPHY( XL),YPHY( YB))
         CALL ENDP
20    CONTINUE
C     Draw Numbers
      DO 2 L=1,NTX
          YY=L/27.0 + YB-0.05
          WRITE(C,'(F6.3,''$'')')TX(L)
          XX=XL-.11
C          WRITE(C,'(G10.3,''$'')')TX(L)
C          XX=XL+.05
          CALL GWRITE(xphy(XX),yphy(YY),1.0,C)
2     CONTINUE
      CALL CLSSEG
      CALL OPNSEG(10)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE RDMENU(XMOUSE,YMOUSE)
#     include "basics.inc"
      character ch

      IF(XPHY(XMOUSE).GE.XLMEN.AND. XPHY(XMOUSE).LE.XRMEN ) THEN
            DO 135 IBOX = 1,nCHOIC
               IF(YPHY(YMOUSE).GE.YBS(IBOX) .AND. YPHY(YMOUSE)
     $             .LE.YTS(IBOX))THEN
                   ICHOIC=IBOX
                   cHoice=item(ichoic)
                   ch=choice
                   RETURN
               ENDIF
135         CONTINUE
      ENDIF
      CALL PRS('Not on menu item$')
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE KEYPAD(TEMP)
#     include "basics.inc"

      call rer(temp)

      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE HARD
      INCLUDE 'devices.inc'
      IFHARD0=.TRUE.
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE NOHARD
      INCLUDE 'devices.inc'
      IFHARD0=.FALSE.
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine getpts(npoint,cspace,iels,isid,xcrved,ycrved)
#     include "basics.inc"

      real cspace(npoint),xcrved(npoint),ycrved(npoint)
      real car(4,4),cpts(3,4),crvs(2,4),zu(3),wk(nxm*nxm),jx(10*nxm*nxm)

      integer npo
      save    npo
      data    npo / 0 /

      save car,zu,jx !     Cardinal spline stuff
      data car/
     +        -0.5, 1.5, -1.5, 0.5,
     +         1.0,-2.5,  2.0,-0.5,
     +        -0.5, 0.0,  0.5, 0.0,
     +         0.0, 1.0,  0.0, 0.0/
      data zu / 0.0 , 0.5 , 1.0 /

      pt1x = x(isid,iels)
      pt1y = y(isid,iels)
      if (isid.eq.4) then
         pt2x = x(1,iels)
         pt2y = y(1,iels)
      elseif (isid.eq.8) then
         pt2x = x(5,iels)
         pt2y = y(5,iels)
      else
         pt2x = x(isid+1,iels)
         pt2y = y(isid+1,iels)
      endif
      IF (CCURVE(ISID,IELS).EQ.' ' .OR.
     $    CCURVE(ISID,IELS).EQ.'M' .OR.
     $    CCURVE(ISID,IELS).EQ.'s'     )THEN
C        Straight line
         DO 10 IX=1,NPOINT
            XCRVED(IX) = PT1X + (PT2X-PT1X)* CSPACE(IX)
            YCRVED(IX) = PT1Y + (PT2Y-PT1Y)* CSPACE(IX)
   10    CONTINUE
      ELSE IF (CCURVE(ISID,IELS).EQ.'C' .OR.
     $        CCURVE(ISID,IELS).EQ.'L' .OR.
     $        CCURVE(ISID,IELS).EQ.'W' ) THEN
C        Find slope of perpendicular
         RADIUS=CURVE(1,ISID,IELS)
         GAP=SQRT( (PT1X-PT2X)**2 + (PT1Y-PT2Y)**2 )
         IF(ABS(2.0*RADIUS).LE.GAP*1.00001) THEN
            IF(RADIUS.GE.0.0) RADIUS=   1.00001 * GAP/2.0
            IF(RADIUS.LT.0.0) RADIUS= - 1.00001 * GAP/2.0
         ENDIF
         XS = PT2Y-PT1Y
         YS = -(PT2X-PT1X)
C        Make length Radius
         XYS=SQRT(XS**2+YS**2)
C        Find Center
         DTHETA = ABS(ASIN(GAP/2./RADIUS))
         PT12X = (PT1X + PT2X)/2
         PT12Y = (PT1Y + PT2Y)/2
         XCENN=PT12X - XS/XYS * RADIUS*COS(DTHETA)
         YCENN=PT12Y - YS/XYS * RADIUS*COS(DTHETA)
         THETA0 = ATAN2((PT12Y-YCENN),(PT12X-XCENN))
C
         IF(CCURVE (ISID,IELS).EQ.'L') THEN
              CALL PRS(' Enter X stretching factor (a/b of ellipse):$')
              CALL PRS('              =1 Original Circle$')
              CALL PRS('              >1 Wider Than Original Circle$')
              CALL PRS('              <1 Thinner Than Original Circle$')
              call rer(ECCENT)
              CALL PRS(
     $        'Enter rotation angle (CCW) in degrees (-45<Rot<45):$')
              call rer(rot)
         ENDIF
         DO 20 IX=1,NPOINT
            R=CSPACE(IX)*2.0-1.0
            IF(RADIUS.LT.0.0)            R=-R
            IF(CCURVE (ISID,IELS).EQ.'W') THEN
               XCRVED(IX) = XCENN + (PT2X-PT1X)*R/2
               YCRVED(IX)=CURVE(1,ISID,IELS)*SIN(CURVE(2,ISID,IELS)
     $         *XCRVED(IX)) + CURVE(3,ISID,IELS)
            ELSE
               XCRVED(IX) = XCENN + ABS(RADIUS) * COS(THETA0 + R*DTHETA)
               YCRVED(IX) = YCENN + ABS(RADIUS) * SIN(THETA0 + R*DTHETA)
            ENDIF
            IF(CCURVE (ISID,IELS).EQ.'L') THEN
C              Correct for the fact that it's an ellipse, not a circle
               XCRVED(IX) = XCENN + (XCRVED(IX) - XCENN) * ECCENT
C              Rotate ellipse
               XROT = XCENN + (XCRVED(IX)-XCENN)*COS(ROT*PI/180.)
     $                      - (YCRVED(IX)-YCENN)*SIN(ROT*PI/180.)
               YROT = YCENN + (XCRVED(IX)-XCENN)*SIN(ROT*PI/180.)
     $                      + (YCRVED(IX)-YCENN)*COS(ROT*PI/180.)
               XCRVED(IX) = XROT
               YCRVED(IX) = YROT
            ENDIF
   20    CONTINUE
      ELSE IF(CCURVE(ISID,IELS).EQ.'S')THEN
C         On split,
C         Cw  corner: keep 1&2  3 is split pt; 4 is what used to be 3
C         CCw corner: keep 3&4  2 is split pt; 1 is what used to be 2
          CPTS(1,1)=CURVE(1,ISID,IELS)
          CPTS(2,1)=CURVE(2,ISID,IELS)
          CPTS(1,2)=PT1X
          CPTS(2,2)=PT1Y
          CPTS(1,3)=PT2X
          CPTS(2,3)=PT2Y
          CPTS(1,4)=CURVE(3,ISID,IELS)
          CPTS(2,4)=CURVE(4,ISID,IELS)
          DO 50 IDIR=1,2
             DO 45 I=1,4
                CRVS(IDIR,I) = 0.0
                DO 40 J=1,4
                   CRVS(IDIR,I) = CRVS(IDIR,I) +
     $             CAR(j,i)*CPTS(IDIR,J)
   40           CONTINUE
   45        CONTINUE
   50     CONTINUE
          DO 70 IU=1,NPOINT
             UU=CSPACE(IU)
             XX = CRVS(1,4)
             YY = CRVS(2,4)
             DO 60 I = 1,3
                XX = XX + UU**(4-I) * CRVS(1,I)
                YY = YY + UU**(4-I) * CRVS(2,I)
   60         CONTINUE
              XCRVED(IU)=XX
              YCRVED(IU)=YY
   70     CONTINUE
      elseif (ccurve(isid,iels).eq.'O') then
C        O - object defined by points
C        o - object defined by function
C
C        Find points closeto straight line
         IOBJ   = CURVE(1,ISID,IELS)
C zero -- 2D kludge
         ZCRV=ZPLANE
         DO 90 IX=1,NPOINT
            XCRVED(IX) = PT1X + (PT2X-PT1X)* CSPACE(IX)
            YCRVED(IX) = PT1Y + (PT2Y-PT1Y)* CSPACE(IX)
            call latchob(xt,yt,xcrved(ix),ycrved(ix),zcrv
     $                                     ,dist2,ko,io,iobj)
            XCRVED(IX) = xt
            YCRVED(IX) = yt
   90    CONTINUE

      elseif (ccurve(isid,iels).eq.'m') then ! midside node

c        write(6,*) 'in getpts:',iels,npoint
c        write(6,*) 'in getpts:',(cspace(k),k=1,npoint)

         zu(1) = 0    ! Uniform input for midside node
         zu(2) = 0.5
         zu(3) = 1
         nzu   = 3
         if (npoint.ne.npo) call gen_int_gz(jx,wk,cspace,npoint,zu,nzu)
         npo = npoint
         
         do k=1,2
            if (k.eq.1) then
               wk(1) = pt1x
               wk(2) = curve(1,isid,iels)
               wk(3) = pt2x
               call mxm(jx,npoint,wk,3,xcrved,1)
            else
               wk(1) = pt1y
               wk(2) = curve(2,isid,iels)
               wk(3) = pt2y
               call mxm(jx,npoint,wk,3,ycrved,1)
            endif
         enddo

c        write(6,*)
c        write(68,*)
c        m=npoint
c        do k=1,npoint
c         write(6,1) iels,isid,k
c    $     ,cspace(k),xcrved(k),ycrved(k),(jx(k+j*m),j=0,2)
c         write(68,1) iels,isid,k
c    $     ,cspace(k),xcrved(k),ycrved(k),(jx(k+j*m),j=0,2)
c  1      format(i9,2i4,1p6e12.4,' xyc')
c        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE DRAWED(IEL,IEDGE,IFLIP)
C     DRAWs EDge.  IFLIP CAuses to draw from end to beginning
C     POSITIVE MEANS CCW; ON VERTICAL STRUTS 9-12, MEANS UPWARD
C     Draws only.  No move or fill.
      DIMENSION XISM(8),YISM(8),CSPACE(100),XCRVED(100),YCRVED(100)
#     include "basics.inc"
C
      ZERO=0.0
      NPOINT=10
      DO 18 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
18    CONTINUE
C     Patch on drised to make drawing simple edge complicated.
      DO 7 IC=1,8
         xism(ic)=x(ic,iel)
         yism(ic)=y(ic,iel)
7     CONTINUE
      IF(IEDGE.GT.8)THEN
C        Vertical strut
         IF(IFLIP.EQ.1)THEN
            CALL MOVEC(XISO(XISM(IEDGE-8),YISM(IEDGE-8),ZERO)
     $               ,YISO(XISM(IEDGE-8),YISM(IEDGE-8),ZERO))
            CALL DRAWC(XISO(XISM(IEDGE-4),YISM(IEDGE-4),ZERO)
     $               ,YISO(XISM(IEDGE-4),YISM(IEDGE-4),ZERO))
         ELSE
            CALL MOVEC(XISO(XISM(IEDGE-4),YISM(IEDGE-4),ZERO)
     $               ,YISO(XISM(IEDGE-4),YISM(IEDGE-4),ZERO))
            CALL DRAWC(XISO(XISM(IEDGE-8),YISM(IEDGE-8),ZERO)
     $               ,YISO(XISM(IEDGE-8),YISM(IEDGE-8),ZERO))
         ENDIF
      ELSE IF(CCURVE(IEDGE,IEL).EQ.' '   .or.
     $        CCURVE(IEDGE,IEL).EQ.'M' ) then
C        Draw straight side
         IF(IFLIP.EQ.1)THEN
            IC=IEDGE+1
            IF(IC.EQ.5)IC=1
            IF(IC.EQ.9)IC=5
         ELSE
            IC=IEDGE
         ENDIF
         CALL MOVEC(XISO(XISM(IEDGE),YISM(IEDGE),ZERO)
     $             ,YISO(XISM(IEDGE),YISM(IEDGE),ZERO))
         CALL DRAWC(XISO(XISM(IC),YISM(IC),ZERO)
     $             ,YISO(XISM(IC),YISM(IC),ZERO))
      ELSE
C        Draw curved side
         CALL GETPTS(NPOINT,CSPACE,IEL,IEDGE,XCRVED,YCRVED)
         IF(IFLIP.EQ.1)THEN
            IBEGIN=1
            IEND  =NPOINT
         ELSE
            IBEGIN=NPOINT
            IEND  =1
         ENDIF
         I=IBEGIN
         CALL MOVEC(XISO(XCRVED(I),YCRVED(I),ZERO)
     $             ,YISO(XCRVED(I),YCRVED(I),ZERO))
         DO 118 I=IBEGIN,IEND,IFLIP
            CALL DRAWC(XISO(XCRVED(I),YCRVED(I),ZERO)
     $                ,YISO(XCRVED(I),YCRVED(I),ZERO))

118      CONTINUE
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SETENV
#     include "basics.inc"
      CHARACTER*10 ARGS(5)
C     Default Environmant
C     Look for nekenv.dat
      IERR=0
      IF(IFVMS)THEN
C        VMS
C        First try Home directory
         CALL OPENF(12,'SYS$LOGIN:DEFAULTS.NEK','old',3,ierr)
C        If no file was found, then look to NEKTON directory
         IF(IERR.NE.0)CALL OPENF(12,'NEKTON:DEFAULTS.NEK','old',3,ierr)
      ELSE
C        Unix
         CALL SYSTEM('cp $HOME/.nekdefaults tmp.nekdefaults')
         CALL OPENF(12,'tmp.nekdefaults','old',3,ierr)
         CALL SYSTEM('date > tmp.date')
         CALL OPENF(14,'tmp.date','old',3,ierr)
      ENDIF
      DATE=' '
      IFPOSTS =.TRUE.
      IF(IERR.EQ.0)THEN
C        We succesfully got a date; put it into the string date
         READ(14,'(A28)')DATE
C pff 7/29/92 photo op.....
c        call blank(date,28)
C        We sucessfully opened a defaults environment file, let's look inside
         DO 271 I=1,1000
            LINE='                                   '//
     $           '                                   '
            READ(12,'(A70)',ERR=13,END=13)LINE
            CALL PARSE(LINE,ARGS,IA)
            IFLASE  =.FALSE.
c           IFPOSTS =.FALSE.
c
c           hardwire postscript as default... (pff 10/14/97)
            IFPOSTS =.TRUE.
            IF(ARGS(1).EQ.'PRINTER')THEN
               IF(ARGS(2).EQ.'POSTSCRIPT')IFPOSTS=.TRUE.
               IF(ARGS(2).EQ.'QMS       ')IFLASE =.TRUE.
               CALL PRS('Printer: '//ARGS(2)//'$')
            ELSE IF(ARGS(1).EQ.'MOUSE')THEN
               IF(ARGS(2).EQ.'TEK3BUTTON')IF3BUT=.TRUE.
               IF(ARGS(2).EQ.'TEK_TABLET')IF3BUT=.FALSE.
               CALL PRS('Mouse: '//ARGS(2)//'$')
            ENDIF
271      CONTINUE
13       CONTINUE
         CLOSE(UNIT=12)
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE INITQMS
C     Initializes QMS Printer (QUIK mode graphics)
      WRITE(45,103)
103   FORMAT(' ^PY^-',/,' ^IGV^PW03')
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE INITPS
C     Initializes POSTSCRIPT Printer
C
      WRITE(45,'(a13)')'%!PSAdobe-1.0'
      WRITE(45,*)' %%Creator: ssr'
      WRITE(45,*)' %%Title:'
      WRITE(45,*)' %%CreationDate: 02/02/1988 22:53:04 '
      WRITE(45,*)' %%Pages: (atend) '
      WRITE(45,*)' %%DocumentFonts: (atend) '
      WRITE(45,*)' %%BoundingBox: (atend)    '
      WRITE(45,*)' %%EndComments              '
      WRITE(45,*)' /saveobj save def           '
      WRITE(45,*)' /center {/yarea exch def /xarea exch def newpath'
      WRITE(45,*)' clippath pathbbox'
      WRITE(45,*)
     $' /ury exch def /urx exch def /lly exch def /llx exch def'
      WRITE(45,*)'  /xoff urx llx sub xarea sub 2 div llx add def'
      WRITE(45,*)'  /yoff ury lly sub yarea sub 2 div lly add def '
      WRITE(45,*)'  xoff yoff translate} def     '
C
c
c
c
c...12/23/92/pff       WRITE(45,*)
c...12/23/92/pff      $' /draw {lineto currentpoint stroke newpath moveto} def'
c...12/23/92/pff       WRITE(45,*)' /move {newpath moveto} def                 '
c...12/23/92/pff 
      WRITE(45,*)
     $' /d {lineto currentpoint stroke newpath moveto} def'
      WRITE(45,*)' /m {newpath moveto} def                 '
c...12/23/92/pff 
c
c
c
      WRITE(45,*)' /font {exch findfont exch scalefont setfont} def '

      WRITE(45,*)'/background {/yarea exch def /xarea exch def gsave'
      WRITE(45,*)' newpath 0 0 moveto 0 yarea lineto xarea yarea lineto'
      WRITE(45,*)' xarea 0 lineto closepath eofill grestore} def'

      WRITE(45,*)'/defaults {setlinewidth setlinecap setlinejoin'
      WRITE(45,*)
     $'setmiterlimit 0 array 0 setdash newpath 0 0 moveto} def'
      WRITE(45,*)' /dot { /y exch def  /x exch def  /s 0.1 def'
      WRITE(45,*)' 		x       y       move'
      WRITE(45,*)' 		x       y s add draw'
      WRITE(45,*)' 		x s add y s add draw'
      WRITE(45,*)' 		x s add y       draw'
      WRITE(45,*)' 		x       y       draw } def'
      WRITE(45,*)'/Times-Roman findfont'
      WRITE(45,*)'12 scalefont'
      write(45,*)'setfont'
      WRITE(45,*)' %%EndProlog                             '
      WRITE(45,*)' %%Page: 1 1                              '
c     WRITE(45,*)
c    $' 575.28 743.803 center 0 0 0 setrgbcolor 10 0 0 0.072 defaults '
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VDRAW(XXIS,YYIS,N)
      DIMENSION XXIS(N),YYIS(N)
      write(6,*) 'move it A',n
      CALL MOVEC(XXIS(1),YYIS(1))
      write(6,*) 'move it B',n
      DO 100 I=2,N
         CALL DRAWC(XXIS(I),YYIS(I))
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TRNGL1(XY1,XY2,XY3)
C
C     Draw the triangle requested by the vector XYZ -
C
C     XYZ(3,3)
C     V(3) - V ranges from 0 to 1
C
#     include "basics.inc"
      DIMENSION XY1(4) ,XY2(4) ,XY3(4)
      DIMENSION XYT1(4),XYT2(4),XYT3(4)
      DIMENSION PT1(3),PT2(3)
C
C     Check if all vertices have the same color - if yes, plot.
C
      IC1=ICOLOR1(XY1(4))
      IC2=ICOLOR1(XY2(4))
      IC3=ICOLOR1(XY3(4))
      IF (IC1.EQ.IC2.AND.IC1.EQ.IC3) THEN
         CALL TRNGL2(XY1,XY2,XY3,IC1)
         RETURN
      ENDIF
C
C     Map to canonical 1-2-3 triangle (via Cyclic shift).
C
      IF (IC1.EQ.IC2) THEN
         CALL COPY(XYT1,XY1,4)
         CALL COPY(XYT2,XY2,4)
         CALL COPY(XYT3,XY3,4)
      ELSEIF (IC1.EQ.IC3) THEN
         CALL COPY(XYT1,XY3,4)
         CALL COPY(XYT2,XY1,4)
         CALL COPY(XYT3,XY2,4)
      ELSE
         CALL COPY(XYT1,XY2,4)
         CALL COPY(XYT2,XY3,4)
         CALL COPY(XYT3,XY1,4)
      ENDIF
C
C     Linear interpolate along legs connected to vertex 3 
C     to find crossover point.
C
      CVAL1=BCOLOR(XYT1(4),XYT3(4))
      cmean = (xyt1(4)+xyt3(4))/2.0
      crad  = abs(xyt1(4)-xyt3(4))/2.0
      ctest = abs(cmean-cval1)
      if (ctest.gt.crad) then
         write(6,*) 'warning: cval1:',xyt1(4),xyt3(4),cval1
         write(6,*) 'cmean,crad,ctest:',cmean,crad,ctest
      endif
      D1VAL=(CVAL1-XYT1(4))/(XYT3(4)-XYT1(4))
      CVAL2=BCOLOR(XYT2(4),XYT3(4))
      D2VAL=(CVAL2-XYT2(4))/(XYT3(4)-XYT2(4))
      PT1(1) = XYT1(1) + D1VAL*(XYT3(1)-XYT1(1))
      PT1(2) = XYT1(2) + D1VAL*(XYT3(2)-XYT1(2))
      PT1(3) = XYT1(3) + D1VAL*(XYT3(3)-XYT1(3))
      PT2(1) = XYT2(1) + D2VAL*(XYT3(1)-XYT2(1))
      PT2(2) = XYT2(2) + D2VAL*(XYT3(2)-XYT2(2))
      PT2(3) = XYT2(3) + D2VAL*(XYT3(3)-XYT2(3))
C
C     Plot 3 new iso-color triangles 
C     (again, note preservation of cyclic symmetry, for normals)
C
      IC1=ICOLOR1(XYT1(4))
      IC2=ICOLOR1(XYT2(4))
      IC3=ICOLOR1(XYT3(4))
      CALL TRNGL2(XYT1,XYT2,PT1 ,IC1)
      CALL TRNGL2(PT1 ,XYT2,PT2 ,IC2)
      CALL TRNGL2(PT1 ,PT2 ,XYT3,IC3)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE TRNGL2(XY1,XY2,XY3,IC)
#     include "basics.inc"
      LOGICAL IFBKFC
      DIMENSION XY1(3),XY2(3),XY3(3)
C
      UUU1=XISO(XY1(1),XY1(2),XY1(3))
      UUU2=XISO(XY2(1),XY2(2),XY2(3))
      UUU3=XISO(XY3(1),XY3(2),XY3(3))
      VVV1=YISO(XY1(1),XY1(2),XY1(3))
      VVV2=YISO(XY2(1),XY2(2),XY2(3))
      VVV3=YISO(XY3(1),XY3(2),XY3(3))
C
C     Backface culling?
C
      IFBKFC=.FALSE.
      IF (IFBKFC) THEN
         UUUU1=UUU2-UUU1
         UUUU2=UUU3-UUU1
         VVVV1=VVV2-VVV1
         VVVV2=VVV3-VVV1
         VCCCC=UUUU1*VVVV2-UUUU2*VVVV1
         IF (VCCCC.LE.0.0) RETURN
      ENDIF
C
      CALL FILLP(-(IC))
      CALL COLOR( (IC))
      CALL BEGINP(UUU1,VVV1)
      CALL DRAWC (UUU2,VVV2)
      CALL DRAWC (UUU3,VVV3)
      CALL DRAWC (UUU1,VVV1)
      CALL ENDP
C
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION ICOLOR1(VAL)
      INTEGER NCOLOR,ICALLD
      SAVE NCOLOR,ICALLD
      DATA NCOLOR,ICALLD /14,0/
      real valmax
      save valmax
      data valmax /0.0/
C
C     Return a integer value between 2 and 16 (inclusive) according
C     to VAL which is between 0 and 1.
C
      icalld=icalld+1
c     if (val.gt.valmax) then
c        write(6,*) 'max color:',val,icalld
c        valmax=val
c     endif
      VALTMP = MAX(0.001,VAL)
      VALTMP = MIN(0.999,VALTMP)
      VALTMP = FLOAT(NCOLOR)*VALTMP
      ICOLOR1= 2+INT(VALTMP)
      if (icolor1.lt.2) write(6,*) 'icolor1:',icolor1,val,valtmp
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION BCOLOR(VAL1,VAL2)
C
C     Return a the break point value corresponding to a change in color,
C     bracketed by VAL1 and VAL2.
C
      REAL BPS(20)
      SAVE BPS
      INTEGER NCOLOR
      SAVE NCOLOR
      INTEGER IFIRST
      SAVE IFIRST
      DATA IFIRST /0/
      DATA NCOLOR /14/
      DATA BPS    /20*0.0/
C
      IF (IFIRST.EQ.0) THEN
         IFIRST=1
         DBP=1.0/FLOAT(NCOLOR)
         DO 10 I=1,NCOLOR+1
            BPS(I)=FLOAT(I-1)*DBP
   10    CONTINUE
      ENDIF
C
C     Use seqential search:
C
      VALMIN=MIN(VAL1,VAL2)
      VALMAX=MAX(VAL1,VAL2)
      DO 20 I=1,NCOLOR+1
         BCOLOR=BPS(I)
         IF (VALMIN.LE.BPS(I).AND.VALMAX.GE.BPS(I)) RETURN 
   20 CONTINUE
      BCOLOR=(VAL1+VAL2)/2.0
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SAVEBOX(XY1,XY2,XY3,XY4,IC)
      DIMENSION XY1(3),XY2(3),XY3(3),XY4(3)
      PARAMETER (NXGM=25)
      COMMON /BOXBOX/ XYB1(3),XYB2(3),XYB3(3,NXGM),XYB4(3,NXGM)
      COMMON /INTBOX/ ICALLD ,ICB
      INTEGER ICALL0
      SAVE    ICALL0
      DATA    ICALL0 /0/
C
      IF (ICALL0.EQ.0) THEN
         ICALLD=0
         ICALL0=1
      ENDIF
C
C     A new box? 
C
      ICALLD=ICALLD+1
      IF (ICALLD.EQ.1) THEN
         CALL COPY(XYB1,XY1,3)
         CALL COPY(XYB2,XY2,3)
         CALL COPY(XYB3,XY3,3)
         CALL COPY(XYB4,XY4,3)
         ICB=IC
         RETURN
      ENDIF
C
C     A different color?
C
      IF (IC.NE.ICB.OR.ICALLD.EQ.NXGM) THEN 
         CALL FLSHBOX
         ICALLD=ICALLD+1
         CALL COPY(XYB1,XY1,3)
         CALL COPY(XYB2,XY2,3)
         CALL COPY(XYB3,XY3,3)
         CALL COPY(XYB4,XY4,3)
         ICB=IC
         RETURN
      ENDIF
C
C     Extend the current box.
C
      CALL COPY(XYB3(1,ICALLD),XY3,3)
      CALL COPY(XYB4(1,ICALLD),XY4,3)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE FLSHBOX
      PARAMETER (NXGM=25)
      COMMON /BOXBOX/ XYB1(3),XYB2(3),XYB3(3,NXGM),XYB4(3,NXGM)
      COMMON /INTBOX/ ICALLD ,ICB
C
C     Flush the box
C
      IF (ICALLD.EQ.0) RETURN
      CALL FILLP(-(ICB))
      CALL COLOR( (ICB))
      CALL BEGINP(XISO(XYB2(1),XYB2(2),XYB2(3))
     $           ,YISO(XYB2(1),XYB2(2),XYB2(3)))
      DO 10 I=1,ICALLD
         CALL DRAWC(XISO(XYB3(1,I),XYB3(2,I),XYB3(3,I))
     $             ,YISO(XYB3(1,I),XYB3(2,I),XYB3(3,I)))
   10 CONTINUE
C
      DO 20 I=ICALLD,1,-1
         CALL DRAWC(XISO(XYB4(1,I),XYB4(2,I),XYB4(3,I))
     $             ,YISO(XYB4(1,I),XYB4(2,I),XYB4(3,I)))
   20 CONTINUE
C
      CALL DRAWC(XISO(XYB1(1),XYB1(2),XYB1(3))
     $          ,YISO(XYB1(1),XYB1(2),XYB1(3)))
      CALL ENDP
      ICALLD=0
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE BOX1(XY1,XY2,XY3,XY4,IC)
#     include "basics.inc"
      LOGICAL IFBKFC
      DIMENSION XY1(3),XY2(3),XY3(3),XY4(3)
C
      UUU1=XISO(XY1(1),XY1(2),XY1(3))
      UUU2=XISO(XY2(1),XY2(2),XY2(3))
      UUU3=XISO(XY3(1),XY3(2),XY3(3))
      UUU4=XISO(XY4(1),XY4(2),XY4(3))
      VVV1=YISO(XY1(1),XY1(2),XY1(3))
      VVV2=YISO(XY2(1),XY2(2),XY2(3))
      VVV3=YISO(XY3(1),XY3(2),XY3(3))
      VVV4=YISO(XY4(1),XY4(2),XY4(3))
C
C     Backface culling?
C
      IFBKFC=.FALSE.
      IF (IFBKFC) THEN
         UUUU1=UUU2-UUU1
         UUUU2=UUU3-UUU1
         VVVV1=VVV2-VVV1
         VVVV2=VVV3-VVV1
         VCCCC=UUUU1*VVVV2-UUUU2*VVVV1
         IF (VCCCC.LE.0.0) RETURN
      ENDIF
C
      CALL FILLP(-(IC))
      CALL COLOR( (IC))
      CALL BEGINP(UUU1,VVV1)
      CALL DRAWC (UUU2,VVV2)
      CALL DRAWC (UUU3,VVV3)
      CALL DRAWC (UUU4,VVV4)
      CALL DRAWC (UUU1,VVV1)
      CALL ENDP
C
      RETURN
      END
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFMENU(XMOUSE,YMOUSE)
C
C     Returns TRUE if MOUSE is in the menu area.
C
#     include "basics.inc"
      IERR=0
      J=0
      YBkey = YPHY(0.7 + J*3./40.)
      XSC=XSCR(XMOUSE)
      YSC=YSCR(YMOUSE)
      I=INT((XSC-     1     )*10*4./3.+1)
      J=INT((YSC-YSCR(YBKEY))*10*4./3.+1)
      IF (I.LT.1.OR.I.GT.4) THEN 
         IFMENU=.FALSE.
      ELSE
         IFMENU=.TRUE.
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VDRAW3(XYZ,N)
      DIMENSION XYZ(4,N)
C
      ICLR=ICOLOR1(XYZ(4,1))
      CALL COLOR(ICLR)
C
      CALL MOVE3( XYZ(1,1),XYZ(2,1),XYZ(3,1) )
      DO 100 I=2,N
         CALL DRAW3( XYZ(1,I),XYZ(2,I),XYZ(3,I) )
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE VDRAW3A(X,Y,Z,N)
      DIMENSION X(N),Y(N),Z(N)
C
      CALL MOVE3( X(1),Y(1),Z(1) )
      DO 100 I=2,N
         CALL DRAW3( X(I),Y(I),Z(I) )
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE ARROW3(X,Y,Z,DXBAR,DYBAR,DZBAR,PLANE,SARROW)
      COMMON/SCALE/XFAC,YFAC,XZERO,YZERO
      CHARACTER PLANE*15
C     Now draw shadow on plane at Zbase (without +DZ)
C     Sort of pale orange in heat map
C     FIX COLOR(WHITE USUALLY, YELLOW FOR AXIS.  <FIX 3-D ARROWHEAD)!!??
      DX=DXBAR*SARROW
      DY=DYBAR*SARROW
      DZ=DZBAR*SARROW
      IF(DX.EQ.0.0 .AND. DY.EQ.0.0 .AND.DZ.EQ.0.0) RETURN
C     Draw base at Z=Zplane (or whichever plane you're plotting on)
      CALL PENW(1)
      CALL COLOR(12)
      CALL MOVE3(X,Y,Z)
      IF(PLANE.EQ.'X-PLANE')CALL DRAW3(X   ,Y+DY,Z+DZ)
      IF(PLANE.EQ.'Y-PLANE')CALL DRAW3(X+DX,Y   ,Z+DZ)
      IF(PLANE.EQ.'Z-PLANE')CALL DRAW3(X+DX,Y+DY,Z   )
      CALL PENW(3)
C     Don't put arrowhead on hard copy (OR ANYWHERE ELSE)
C     Draw 3-dimensional arrow
      CALL COLOR(1)
      CALL MOVE3(X,Y,Z)
      CALL DRAW3(X+DX,Y+DY,Z+DZ)
C     Wierd arrowheads
      CALL DRAW3(X+0.85*DX-.1*DY ,Y+0.85*DY+0.1*DX*yfac/xfac,Z+DZ*0.85)
      CALL MOVE3(X+DX,Y+DY,Z+DZ)
      CALL DRAW3(X+0.85*DX+0.1*DY,Y+0.85*DY-0.1*DX*yfac/xfac,Z+DZ*0.85)
      IF(DX.EQ.0.0 .AND. DY.EQ.0.0) THEN
C        Kludge to get arrowhead on vertical arrow
         CALL MOVE3(X,Y,Z+DZ)
         CALL DRAW3(X,Y+0.1*Dz,Z+DZ*0.85)
         CALL MOVE3(X,Y,Z+DZ)
         CALL DRAW3(X,Y-0.1*Dz,Z+DZ*0.85)
      ENDIF
      RETURN
      END
c    End of PLOT subroutines, common to preprocessor and postprocessor*******
c
c-----------------------------------------------------------------------
      subroutine draw_circ(x3,y3,z3,r)
c     draws a circle around XX,YY,ZZ
c
      xm=xiso(x3-r,y3,z3)
      ym=yiso(x3-r,y3,z3)
      xp=xiso(x3+r,y3,z3)
      yp=yiso(x3+r,y3,z3)
      dx = (xp-xm)**2 + (yp-ym)**2
c
      xm=xiso(x3,y3-r,z3)
      ym=yiso(x3,y3-r,z3)
      xp=xiso(x3,y3+r,z3)
      yp=yiso(x3,y3+r,z3)
      dy = (xp-xm)**2 + (yp-ym)**2
c
      xm=xiso(x3,y3,z3-r)
      ym=yiso(x3,y3,z3-r)
      xp=xiso(x3,y3,z3+r)
      yp=yiso(x3,y3,z3+r)
      dz = (xp-xm)**2 + (yp-ym)**2
c
      dr = max(dx,dy)
      dr = max(dr,dz)
      dr = sqrt(dr)
c
      xc=xiso(x3,y3,z3)
      yc=yiso(x3,y3,z3)
      write(6,1) xc,yc,x3,y3,z3,dr,r
    1 format(8f10.5)
c
      x = xc+dr
      y = yc
      call movec(x,y)
c
      n = 100
      pi2n = 8.*atan(1.0)/n
      do i=1,n
         theta = i*pi2n
         x = xc + dr*cos(theta)
         y = xc + dr*sin(theta)
         call drawc(x,y)
      enddo
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE ADJVUE(ifadjust)
#     include "basics.inc"
      logical ifadjust
      character*1 adj(2)
      character*3 adj2
      equivalence (adj2,adj(2))
C
      CALL PRS('Enter h,j,k,l or <cr> to return$')
      call blank(adj,4)
      call res  (adj,4)
c
      ifadjust=.true.
      IF(adj(1).EQ.' ')THEN
        ifadjust=.false.
        return
      ENDIF
c
c
c     Get magnitude
c
      idelta = 1
      if (adj2.ne.' ') then
         read(adj(2),1,err=2)  jdelta
         idelta = jdelta
      endif
    1 format(i1)
    2 continue
c
c     theta:   angle above horizon in degrees
c     phi:     azimuthal angle from x axis
c
      if (adj(1).eq.'h') phi   = phi   + idelta
      if (adj(1).eq.'j') theta = theta - idelta
      if (adj(1).eq.'k') theta = theta + idelta
      if (adj(1).eq.'l') phi   = phi   - idelta
c
C     Now set viewing transformation
      THETAr=THETA*PI/180.0
      PHIr  =PHI  *PI/180.0
C
C     Direction toward eye (normal vector)
      VHOBS(1)=COS(THETAR)*COS(PHIR)
      VHOBS(2)=COS(THETAR)*SIN(PHIR)
      VHOBS(3)=SIN(THETAR)
C     Perpendicular on x-y plane
      XHOBS(1)=-1.0*SIN(PHIR)
      XHOBS(2)=COS(PHIR)
      XHOBS(3)=0.0
C     Perpendicular to above two
      YHOBS(1)=-1.0*SIN(THETAR)*COS(PHIR)
      YHOBS(2)=-1.0*SIN(THETAR)*SIN(PHIR)
      YHOBS(3)=COS(THETAR)
C     Rescale screen-world factors so it fits on screen
      IF (.NOT.IFZOOM)CALL RESCAL
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE XDIAMD(XX,YY)
#     include "basics.inc"
c     draws a diamond around XX,YY
      CALL MOVEC(XX+.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWC(XX-.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWC(XX-.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWC(XX+.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWC(XX+.0025*XFAC ,YY+.0025*YFAC)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE SDIAMD(XX,YY)
#     include "basics.inc"
c     draws a diamond around XX,YY
      CALL MOVESC(XX+.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWSC(XX-.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWSC(XX-.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWSC(XX+.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWSC(XX+.0025*XFAC ,YY+.0025*YFAC)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE RDIAMD(XX,YY)
#     include "basics.inc"
c     draws a diamond around XX,YY
      XLFAC=1.0
      YLFAC=1.0
      CALL MOVESC(XX+.0025*XLFAC ,YY+.0025*YLFAC)
      CALL DRAWSC(XX-.0025*XLFAC ,YY+.0025*YLFAC)
      CALL DRAWSC(XX-.0025*XLFAC ,YY-.0025*YLFAC)
      CALL DRAWSC(XX+.0025*XLFAC ,YY-.0025*YLFAC)
      CALL DRAWSC(XX+.0025*XLFAC ,YY+.0025*YLFAC)
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine diam2(xmouse,ymouse,iclr)
#     include "basics.inc"
c     draws a diamond around xmouse,ymouse
      call color(iclr)
      call movesc(xscr(xmouse)-.015,yscr(ymouse))
      call drawsc(xscr(xmouse),yscr(ymouse)-.015)
      call drawsc(xscr(xmouse)+.015,yscr(ymouse))
      call drawsc(xscr(xmouse),yscr(ymouse)+.015)
      call drawsc(xscr(xmouse)-.015,yscr(ymouse))
      return
      end
c-----------------------------------------------------------------------
      subroutine xdiam2(xmouse,ymouse,iclr)
#     include "basics.inc"
c     draws a plus around xmouse,ymouse
      call color(iclr)
      call movesc(xscr(xmouse)-.015,yscr(ymouse))
      call drawsc(xscr(xmouse)+.015,yscr(ymouse))
      call movesc(xscr(xmouse),yscr(ymouse)-.015)
      call drawsc(xscr(xmouse),yscr(ymouse)+.015)
      return
      end
c-----------------------------------------------------------------------
