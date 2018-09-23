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
      subroutine setgrd
#     include "basics.inc"
      CHARACTER*26 CCART,CPOLAR,CGRIDX,CGRIDY,CGRIDN,CGRIDR,BLNK
      CHARACTER*26 CGRDX,CGRDY,COBJCT,CHEX
      CHARACTER*1  BLNK1(26)
      EQUIVALENCE (BLNK1,BLNK)
      LOGICAL IFCHNG
      DATA BLNK1 /26*' '/
C
C     Set up grids for polar or cartesian based input
C
      ifobjs = .true.
      nobjs  = max(nobjs,1)
      ONE = 1.0
      PI = 4.0*ATAN(ONE)
      IFCHNG=.FALSE.

      if (ngrid.eq.0) then   ! Set defaults for polar grid
         ngrid =8
         grida = 2.0*pi/float(ngrid)
         gridr = griddx
      endif

 1000 CONTINUE

      CCART  = BLNK
      CHEX   = BLNK
      CPOLAR = BLNK
      CGRIDX = BLNK
      CGRIDY = BLNK
      CGRIDN = BLNK
      CGRIDR = BLNK
      CGRDX  = BLNK
      CGRDY  = BLNK
      COBJCT = BLNK
      WRITE(CCART , 2) IFGRDC
      WRITE(CHEX  , 3) IFGRDH ! Hex
      WRITE(CPOLAR, 4) IFGRDP
      WRITE(CGRIDX, 5) GRIDDX
c     WRITE(CGRIDY, 6) GRIDDY
      WRITE(CGRIDN, 7) NGRID
      WRITE(CGRIDR, 8) GRIDR
      WRITE(CGRDX , 9) GRIDXP
      WRITE(CGRDY ,10) GRIDYP
      if (ifobjs) write(cobjct,11) (ifobjg(i),i=1,nobjs)
    2 FORMAT('CARTESIAN',L4)
    3 FORMAT('HEXAGONAL',L4)
    4 FORMAT('POLAR    ',L4)
    5 FORMAT('DX =',G11.4)
    6 FORMAT('DY =',G11.4)
    7 FORMAT('# DIV =',I3)
    8 FORMAT('dR =',G11.4)
    9 FORMAT('Xo=',2G11.4)
   10 FORMAT('Yo=',2G11.4)
   11 FORMAT('OBJ: ',20L1)
C
      ITEM(1) ='MAIN MENU'
      ITEM(2) =CCART
      ITEM(3) =CHEX
      ITEM(4) =CPOLAR
      ITEM(5) =CGRIDX
c     ITEM(6) =CGRIDY
      ITEM(6) =CGRIDN
      ITEM(7) =CGRIDR
      ITEM(8) =CGRDX
      ITEM(9) =CGRDY
      ITEM(10)=COBJCT
      NCHOIC = 9
      IF (IFOBJS) NCHOIC=NCHOIC+1
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'SET GRID')
C
      IF (CHOICE.EQ.'MAIN MENU') THEN
         IF (IFCHNG) THEN
            CALL REFRESH
            CALL DRGRID
            DO 50 IEL=1,NEL
               CALL DRAWEL(IEL)
   50       CONTINUE
         ENDIF
         return
      ELSE
         IFCHNG=.TRUE.
      ENDIF
      IF (CHOICE.EQ.CCART)  IFGRDC=.NOT.IFGRDC
      IF (CHOICE.EQ.CHEX)   IFGRDH=.NOT.IFGRDH
      IF (CHOICE.EQ.CPOLAR) IFGRDP=.NOT.IFGRDP
      IF (IFGRDP.AND.NGRID.EQ.0) THEN
         WRITE(S,401)
         CALL PRS(S//'$')
         call rer(RGRID)
         NGRID = IFIX(RGRID)
         NGRID = MAX(NGRID,3)
         GRIDA = 2.0*PI/FLOAT(NGRID)
      ENDIF
      IF (IFGRDP.AND.GRIDR.EQ.0.0) THEN
         WRITE(S,501)
         CALL PRS(S//'$')
         call rer(GRIDR)
      ENDIF
      IF (CHOICE.EQ.CGRIDX) THEN
         WRITE(S,301)
         CALL PRS(S//'$')
  301    FORMAT(2X,'Type new value of DX:')
         call rer(GRIDDX)
C        Don't change GRID !
         GRIDT =GRIDDX/XFAC
C        GRIDDX=XFAC*GRIDT  kludge for now
         GRIDDY=YFAC*GRIDT
         GRIDDY=GRIDDX
      ENDIF
      IF (CHOICE.EQ.CGRIDY) THEN
         WRITE(S,302)
         CALL PRS(S//'$')
  302    FORMAT(2X,'Type new value of DY:')
         call rer(GRIDDY)
C        Don't change GRID !
         GRIDT =GRIDDY/XFAC
         GRIDDX=GRIDDY
c        GRIDDX=XFAC*GRIDT
C        GRIDDY=YFAC*GRIDT  kludge for now
      ENDIF
      IF (CHOICE.EQ.CGRIDN) THEN
         WRITE(S,401)
         CALL PRS(S//'$')
  401    FORMAT(2X,'Input number of radial divisions for polar grid.')
         call rer(RGRID)
         NGRID = IFIX(RGRID)
         NGRID = MAX(NGRID,3)
         GRIDA = PI/FLOAT(NGRID)
      ENDIF
      IF (CHOICE.EQ.CGRIDR) THEN
  500    WRITE(S,501)
         CALL PRS(S//'$')
  501    FORMAT(2X,'Input value of dR for polar grid.')
         call rer(GRIDR)
         IF (GRIDR.LE.0.0) GOTO 500
      ENDIF
      IF (CHOICE.EQ.CGRDX) THEN
         WRITE(S,601)
         CALL PRS(S//'$')
  601    FORMAT(2X,'Input X-ordinate for center of polar grid.')
         call rer(GRIDXP)
      ENDIF
      IF (CHOICE.EQ.CGRDY) THEN
         WRITE(S,602)
         CALL PRS(S//'$')
  602    FORMAT(2X,'Input Y-ordinate for center of polar grid.')
         call rer(GRIDYP)
      ENDIF
      IF (CHOICE.EQ.COBJCT) THEN
         if (nobjs.gt.1) then
            WRITE(S,701)
            CALL PRS(S//'$')
  701       FORMAT(1x,
     $      'Input object # to change (ALL: 0=F,99=T,-1=toggle)')
            call rei(iobj)
            if (1.le.iobj.and.iobj.le.nobjs) then
               ifobjg(iobj)=.not.ifobjg(iobj)
            elseif (iobj.eq.-1.or.iobj.eq.0.or.iobj.ge.99) then
               do i=1,nobjs
                  if (iobj.eq.-1) ifobjg(i)=.not.ifobjg(i)
                  if (iobj.eq. 0) ifobjg(i)=.false.
                  if (iobj.ge.99) ifobjg(i)=.true.
               enddo
            endif
         ELSE
            IFOBJG(1)=.NOT.IFOBJG(1)
         ENDIF
      ENDIF
      GOTO 1000
      END
c-----------------------------------------------------------------------
      subroutine drgrid
#     include "basics.inc"
c
      common /chexb/ hexscale(8,2),hexbase(8,2)
c
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
c
      CALL GINGRD(GRID)
      CALL COLOR(14)
c     write(6,*) 'ifgrdc,grid:',ifgrdc,grid
      IF (IFGRDC .AND. GRID .GT. .0099) THEN
C        PUT GRAPH PAPER ON SCREEN
         IF (GRIDDX.EQ.0.OR.GRIDDX.NE.GRIDDY) THEN
C           Then we're in Prenek, first time to draw grid
c             write(6,*) 'weve got the wrong one:',griddx,griddy
              DO 2 I=1,IFIX(1/GRID - 1)
                 XX=I*GRID * .9963
C                Don't write grid lines in bottom 1/10 of screen.
                 IF (XX.GT. 0.11) THEN
                    CALL MOVESC(0.01,XX)
                    CALL DRAWSC(0.99,XX)
                 ENDIF
                 CALL MOVESC(XX,0.11)
                 CALL DRAWSC(XX,0.99)
    2         CONTINUE
         ELSE
C           We use GRIDDX which measures the Grid in physical dimensions
C           and which can be shifted with the zero
            ZT = YPHY(.995)
            ZB = YPHY(.005)
            ZL = XPHY(.005)
            ZR = XPHY(.995)
            IXLEFT=INT(ZL/GRIDDX)
            IXRGHT=INT(ZR/GRIDDX)
            IF (ZL.GT.0) IXLEFT=IXLEFT+1
            IF (ZR.LT.0) IXRGHT=IXRGHT-1
c         write(6,*) 'weve got the right one:',griddx,griddy,iffull
            DO 3 I=IXLEFT,IXRGHT
               XX=FLOAT(I)*GRIDDX+GRIDXP
               CALL MOVEC(XX,ZT)
               CALL DRAWC(XX,ZB)
    3       CONTINUE
            IYBOTT=INT(ZB/GRIDDY)
            IYTOP =INT(ZT/GRIDDY)
            IF (ZB.GT.0) IYBOTT=IYBOTT+1
            IF (ZT.LT.0) IYTOP =IYTOP -1
            DO 4 I=IYBOTT,IYTOP
               YY=FLOAT(I)*GRIDDY+GRIDYP
               CALL MOVEC(ZL,YY)
               CALL DRAWC(ZR,YY)
    4       CONTINUE
         ENDIF
      ENDIF
C
      IF (IFGRDP) THEN
         WFR= 1.3
         WT = YPHY(.995)
         WB = YPHY(.005)
         WL = XPHY(.005)
c        WR = XPHY(.995)
         WR = XPHY(1.1)
         IF (IFFULL) WR = XPHY(WFR)
C        set range for radius based on windows in physical coordinates
         RADMIN = 0.0
         RADMAX = 0.0
         RADIUS = SQRT(WL**2 + WB**2)
         RADMAX = MAX(RADMAX,RADIUS)
         RADIUS = SQRT(WR**2 + WB**2)
         RADMAX = MAX(RADMAX,RADIUS)
         RADIUS = SQRT(WR**2 + WT**2)
         RADMAX = MAX(RADMAX,RADIUS)
         RADIUS = SQRT(WL**2 + WT**2)
         RADMAX = MAX(RADMAX,RADIUS)
C        draw cirular segments at each radius
         NRAD  = IFIX(RADMAX/GRIDR) + 1
         NSEGS = 44
         DANG  = 2.0*PI/FLOAT(NSEGS)
         DO 220 IRAD = 1,NRAD
            RADIUS = GRIDR*FLOAT(IRAD)
            ANGLE = 0.0
            XPHY0 = RADIUS*COS(ANGLE)+GRIDXP
            YPHY0 = RADIUS*SIN(ANGLE)+GRIDYP
            CALL MOVEC(XPHY0,YPHY0)
            DO 200 ISEG = 1,NSEGS
               ANGLE = DANG*FLOAT(ISEG)
               XPHY1 = RADIUS*COS(ANGLE)+GRIDXP
               YPHY1 = RADIUS*SIN(ANGLE)+GRIDYP
               CALL DRAWC(XPHY1,YPHY1)
  200       CONTINUE
  220    CONTINUE
C        draw radial segments at each angle
         XPHY0 = GRIDXP
         YPHY0 = GRIDYP
         DANG  = 2.0*PI/FLOAT(NGRID)
         DO 320 IRAD = 1,NGRID
            ANGLE = DANG*FLOAT(IRAD)
            XPHY1 = RADMAX*COS(ANGLE)+GRIDXP
            YPHY1 = RADMAX*SIN(ANGLE)+GRIDYP
            CALL MOVEC(XPHY0,YPHY0)
            CALL DRAWC(XPHY1,YPHY1)
  320    CONTINUE
      ENDIF
c
      if (ifgrdh) then
         WFR= 1.3
         WT = YPHY(.995)
         WB = YPHY(.005)
         WL = XPHY(.005)
         WR = XPHY(1.25)
c        WR = XPHY(.995)
         IF (IFFULL) WR = XPHY(WFR)
c
         dxhex = 3*griddx
         dyhex = 3.
         dyhex = sqrt(dyhex)*griddx
         imin = wl/dxhex
         imax = wr/dxhex
         jmin = wb/dyhex
         jmax = wt/dyhex
         call cmult2(hexscale,hexbase,griddx,16)
         do j=jmin,jmax
            j2 = j + 2*abs(jmin)
c           if (mod(j2,2).eq.0) then
               do i=imin,imax
                  call cmult2(hexscale,hexbase,griddx,16)
                  xi = i*dxhex
                  yj = j*dyhex
                  call cadd  (hexscale(1,1),xi,8)
                  call cadd  (hexscale(1,2),yj,8)
                  call vdraw (hexscale(1,1),hexscale(1,2),8)
               enddo
c           endif
         enddo
      endif
c
C   DRAW AXIS
C! Axisymmetric
      IF(IFAXIS) THEN
              CENTRL=YSCR(0.0)
              CENTRx=XSCR(0.0)
              IF(CENTRL.LT.0.95 .AND. CENTRL .GT. 0.0) THEN
                CALL COLOR(10)
                CALL MOVESC(0.03,CENTRL)
                CALL DRAWSC(0.97,CENTRL)
                CALL GSWRIT(0.94,CENTRL+.01,1.0,'C/L$')
                IF(CENTRX.LT.0.95 .AND. CENTRX .GT. 0.0) THEN
                  CALL MOVESC(CENTRX,CENTRL-.02)
                  CALL DRAWSC(CENTRX,CENTRL+.02)
                ENDIF
              ENDIF
      ELSE
C! Planar 2-D
              CENTRX=XSCR(0.0)
              CENTRY=YSCR(0.0)
              IF(CENTRX.GT.0.0 .AND.CENTRY.GT.0.0 .AND.
     $           CENTRX.LE.0.85.AND.CENTRY.LE.0.85) THEN
C               Draw Axis
                CALL COLOR(10)
                CALL MOVE(0.0,0.0)
                CALL DRAW(0.0,0.1*YFAC)
                CALL MOVE(0.0,0.0)
                CALL DRAW(0.1*XFAC,0.0)
                CALL GWRITE(0.1 *XFAC ,0.01*YFAC,1.0,'X$')
                CALL GWRITE(0.01*XFAC,0.1 *YFAC,1.0,'Y$')
              ENDIF
      ENDIF

      call drwobjs

      return
      end
c-----------------------------------------------------------------------
      subroutine refresh
#     include "basics.inc"
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
      CALL NOHARD
      CALL CLSSEG
      CALL OPNSEG(10)
C     WRITE ON SCREEN
C     Draw Keypad
c     CALL DRKEY
C     KEYPAD OFF; COVER ON
      CALL SGVIS(12,0)
      CALL SGVIS(11,1)
      return
      END
c-----------------------------------------------------------------------
      subroutine movec(x1,y1)
C
C     Clip image denoted by vectors in { (QX,QY),i=1,NBUF } to be
C     within window specified by WT,WB,WL,WR.
C
      include 'devices.inc'
      LOGICAL IFTMP
C
C     put this extra move command in because it doesn't hurt.
      IFTMP=IFHARD
      IFHARD=.FALSE.
      CALL MOVE(X1,Y1)
      IFHARD=IFTMP
      XMOVEC = X1
      YMOVEC = Y1
      return
      END
c-----------------------------------------------------------------------
      subroutine drawc(x1,y1)
C
C     Clip image denoted by vectors in { (QX,QY),i=1,NBUF } to be
C     within window specified by WT,WB,WL,WR.  Note number of
C     vectors will possibly be reduced if many of the vectors lie
C     completely outside the window.
C
      COMMON /CLIP/ XMOVEC,YMOVEC,WT,WB,WL,WR
      REAL NX(4),NY(4),FX(4),FY(4)
      LOGICAL IN,LASTIN
      REAL XLAST,YLAST
      SAVE XLAST,YLAST
      DATA XLAST,YLAST /10.0E+04,10.0E+04/
c
c     compute inward facing normals for CYBECK algorithm
c
      nx(1) =  0
      nx(2) = -1
      nx(3) =  0
      nx(4) =  1
      ny(1) =  1
      ny(2) =  0
      ny(3) = -1
      ny(4) =  0
c
      fx(1) = wl
      fx(2) = wr
      fx(3) = wr
      fx(4) = wl
      fy(1) = wb
      fy(2) = wb
      fy(3) = wt
      fy(4) = wt
C
C     Cyrus-Beck algorithm arrays are all set.
C
      X = XMOVEC
      Y = YMOVEC
      IN     = .FALSE.
      IF ((X.GE.WL.AND.X.LE.WR).AND.(Y.GE.WB.AND.Y.LE.WT)) IN=.TRUE.
      IF (IN.AND.(X.NE.XLAST.OR.Y.NE.YLAST)) CALL MOVE(X,Y)
C
      XO = XMOVEC
      YO = YMOVEC
      X  = X1
      Y  = Y1
C
C     Check if vector has zero length
C
      DOT = (X-XO)**2 + (Y-YO)**2
      IF (DOT.LT.1.0E-11) GOTO 1000
C
C     Check if both points are in
C
      LASTIN = IN
      IN     = .FALSE.
      IF ((X.GE.WL.AND.X.LE.WR).AND.(Y.GE.WB.AND.Y.LE.WT)) IN=.TRUE.
      IF (IN.AND.LASTIN) THEN
         CALL DRAW(X,Y)
         XLAST=X
         YLAST=Y
         GOTO 1000
      ENDIF
C
C     Check if both points are trivially out
C
      IF ( (X.LT.WL.AND.XO.LT.WL) .OR.
     $     (X.GT.WR.AND.XO.GT.WR) .OR.
     $     (Y.GT.WT.AND.YO.GT.WT) .OR.
     $     (Y.LT.WB.AND.YO.LT.WB) ) GOTO 1000
C
C     We may have a partial line segment.
C
C     Employ the Cyrus-Beck code from George Anagnostou
C
      dx = x-xo
      dy = y-yo
      tl = 0.0
      tu = 1.0
C
C     Loop over the 4 window edges
C
      do 505 ii = 1,4
         wx = xo - fx(ii)
         wy = yo - fy(ii)
         ddotn = dx*nx(ii)+dy*ny(ii)
         wdotn = wx*nx(ii)+wy*ny(ii)
         if (abs(ddotn).lt.1.0e-06) goto 502
         t = -wdotn/ddotn
         if (ddotn .gt. 0.0) goto 501
         if (t     .lt. 0.0) goto 505
         tu = min(t,tu)
         goto 505
c
c        looking for the lower limit
c
  501    if (t.gt.1.0) goto 505
         tl = max(t,tl)
         goto 505
  502    if (wdotn.lt.0.0) goto 1000
  505 continue
      if (  tl .gt. tu  ) goto 1000
c
c     If we arrived here, then we have clipped a segment.
c
      IF (tl.gt.0.0) THEN
C
C        We need an extra move (UP) command
C
         xp = xo + (x-xo)*tl
         yp = yo + (y-yo)*tl
         call move(xp,yp)
      ENDIF
C
C     New point is within box
C
      xp = xo + (x-xo)*tu
      yp = yo + (y-yo)*tu
      call draw(xp,yp)
      xlast=xp
      ylast=yp
c
 1000 CONTINUE
C
C     Clean up.
C
      XMOVEC = X1
      YMOVEC = Y1
      return
      END
c-----------------------------------------------------------------------
      subroutine filter(xmouse,ymouse,xscreen,yscreen)
#     include "basics.inc"
c
      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) then
         icalld = icalld + 1
         one = 1.
         pi  = 4.*atan(one)
         pi3 = pi/3.
         call rzero(xhex,6)
         call rzero(yhex,6)
c
         do i=0,5
            angle = pi + i*pi3
            k     = 1+i
            xhex(k) = cos(angle)
            yhex(k) = sin(angle)
         enddo
      endif

      write(6,*) iobjct,' FILTER IOBJ CAA'

      xmse0 = xphy(xscreen)
      ymse0 = yphy(yscreen)
      zmse0 = zplane

      xmouse=xmse0 ! default return values
      ymouse=ymse0
      zmouse=zmse0

c     Check objects first, then background grids.

      iobjct = 0
      write(6,*) iobjct,' FILTER IOBJ AAA',griddx,griddy

      tolobj = (griddx**2+griddy**2)/10.0
      do 10 iobj=1,nobjs
         if (ifobjg(iobj)) then ! grab object?
            call latchob(xms1,yms1,xmse0,ymse0,zmse0,dist2,k,i,iobj)
            write(6,9) k,i,iobj,xms1,yms1,xmse0,ymse0,dist2
   9        format(3i5,1p5e12.4,' OBJECT')
            if (dist2.lt.tolobj) then
               write(6,*) dist2,tolobj,iobjct,' filter obj'
               iobjct=iobj
               tolobj=dist2
               xmouse=xms1
               ymouse=yms1
               zmouse=zmse0
c              zmouse = sqrt(-1./(iobj-1))
            endif
         endif
         write(6,*) iobjct,' FILTER IOBJ BAA',iobj
   10 continue
      if (iobjct.gt.0) return

      if (ifgrdh) then
         dxhx = 3*griddx
         dyhx = 3
         dyhx = sqrt(dyhx)*griddx
c
         xphyc  = xphy(xscreen)-gridxp
         yphyc  = yphy(yscreen)-gridyp
         i      = (xphyc+.5*dxhx)/dxhx
         j      = (yphyc+.5*dyhx)/dyhx
         d2mn   = 1.e8*dxhx*dxhx
         kmin   = 1
         do k=1,6
            xtmp = i*dxhx + griddx*xhex(k)
            ytmp = j*dyhx + griddx*yhex(k)
            dst2 = (xphyc-xtmp)**2 +(yphyc-ytmp)**2
            if (dst2.lt.d2mn) then
               d2mn = dst2
               kmin = k
               xmouse  = xtmp + gridxp
               ymouse  = ytmp + gridyp
               zmouse  = zmse0
            endif
            write(6,*) i,'xms',xhex(k),xtmp,k
            write(6,*) j,'yms',yhex(k),ytmp,dst2
         enddo
c
         if (ifgrdp) then !           polar grid
            radius = sqrt(xphyc**2+yphyc**2)
            radius = anint(radius/gridr)*gridr
            if (radius.gt.0) then
               angle = atan2(yphyc,xphyc)
               angle = anint(angle/grida)*grida
            else
               angle = 0.0
            endif
            xphyp = radius*cos(angle)+gridxp
            yphyp = radius*sin(angle)+gridyp
            dpolar2 = (xphyp -xmse0)**2 + (yphyp -ymse0)**2
            dhex2   = (xmouse-xmse0)**2 + (ymouse-ymse0)**2
c           write(6,*) xmse0,xmouse,xphyp,' Xmouse'
c           write(6,*) ymse0,ymouse,yphyp,' Ymouse'
            if (dpolar2.lt.dhex2) then
               xmouse  = xphyp
               ymouse  = yphyp
            endif
         endif
c
      elseif (ifgrdc.or.ifgrdp) then
         IF (IFGRDC) THEN
C           cartesian grid
            IF (GRIDDX.EQ.0.0.OR.GRIDDY.EQ.0.0) THEN
               XSCRNC = ANINT(XSCREEN/GRID)*GRID
               YSCRNC = ANINT(YSCREEN/GRID)*GRID
               XPHYC  = XPHY(XSCRNC)
               YPHYC  = YPHY(YSCRNC)
               XMSE0  = XPHYC
               YMSE0  = YPHYC
            ELSE
               XPHYC  = XPHY(XSCREEN)-GRIDXP
               YPHYC  = YPHY(YSCREEN)-GRIDYP
               XPHYC  = ANINT(XPHYC/GRIDDX)*GRIDDX+GRIDXP
               YPHYC  = ANINT(YPHYC/GRIDDY)*GRIDDY+GRIDYP
               XMSE0  = XPHYC
               YMSE0  = YPHYC
            ENDIF
         ENDIF
         IF (IFGRDP) THEN
C           polar grid
            XPHYP  = XPHY(XSCREEN)-GRIDXP
            YPHYP  = YPHY(YSCREEN)-GRIDYP
            RADIUS = SQRT(XPHYP**2+YPHYP**2)
            RADIUS = ANINT(RADIUS/GRIDR)*GRIDR
            IF (RADIUS.GT.0) THEN
               ANGLE = ATAN2(YPHYP,XPHYP)
               ANGLE = ANINT(ANGLE/GRIDA)*GRIDA
            ELSE
               ANGLE = 0.0
            ENDIF
            XPHYP = RADIUS*COS(ANGLE)+GRIDXP
            YPHYP = RADIUS*SIN(ANGLE)+GRIDYP
            XMSE0 = XPHYP
            YMSE0 = YPHYP
         ENDIF
         IF (IFGRDC.AND.IFGRDP) THEN
C           check for the closest
            RADC = (XPHYC-XPHY0)**2 + (YPHYC-YPHY0)**2
            RADP = (XPHYP-XPHY0)**2 + (YPHYP-YPHY0)**2
            IF (RADC.LT.RADP) THEN
               XMSE0 = XPHYC
               YMSE0 = YPHYC
            ELSE
               XMSE0 = XPHYP
               YMSE0 = YPHYP
            ENDIF
         ENDIF
         XMOUSE=XMSE0
         YMOUSE=YMSE0
         ZMOUSE=ZMSE0
      endif
C     grid filtering complete
      return
      end
c-----------------------------------------------------------------------
      subroutine getobjs
#     include "basics.inc"

      open(unit=39,file='obj.dat')

      nobjs = 0
      l=1
      do iobj=1,mobj
         read(39,*,end=99,err=99) npts(iobj)
         n = abs(npts(iobj))  ! < 0 --> closed object
         do k=1,n
            read(39,*,end=99,err=99) xobj(l),yobj(l)
            l=l+1
            if (l.gt.lobj) then
               call prsii ('ERROR: increase lobj in basics.inc$',k,n)
               close(39)
               return
            endif
         enddo
         nobjs = nobjs+1
         ccobjs(nobjs) = 'o'
         ifobjs        = .true.
         ifobjg(nobjs) = .true.  ! Turn on object, as default
         ifgrdc        = .false. ! Turn off Cartesian, as default
      enddo
      close(39)
      call redraw_mesh
      call prsis('Found$',nobjs,' objects in obj.dat file.$')
      return

   99 continue
      close(39)
      call redraw_mesh
      call prsis('Found$',nobjs,' objects.$')
      return
      end
c-----------------------------------------------------------------------
      subroutine drwobjs
#     include "basics.inc"

      write(6,*) 'drwobjs: ',nobjs

      do i=1,nobjs
         if (ifobjg(i)) then
            call drwobj(i)
         else
            call drwobj_dash(i)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drwobj_dash(iobj)
#     include "basics.inc"

      ioff=1
      do i=1,iobj-1
         ioff=ioff+npts(i)
      enddo

      ndash = 64
      do ipass=1,5
         n_per_dash = npts(i)/(2*ndash)
         if (n_per_dash.gt.0) goto 10
         ndash = ndash/2
      enddo
      return  ! Not enough points to draw dashed line

   10 continue

      iclro = mod1(iobj,7)+5
      call color(iclro)

      imax = ioff-1+npts(iobj)
      i    = ioff
      do k=1,2*ndash
         ilast = i+n_per_dash
         if (ilast.gt.imax) n_per_dash = imax-i
         call vdraw(xobj(i),yobj(i),n_per_dash)
         i = i+2*n_per_dash
         if (i.ge.imax) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine drwobj(iobj)
#     include "basics.inc"

      ioff=1
      do 100 i=1,iobj-1
         ioff=ioff+npts(i)
  100 continue

      iclro = mod1(iobj,7)+5
      call color(iclro)
      call vdraw(xobj(ioff),yobj(ioff),npts(iobj))

      if (npts(iobj).le.60) then
         do i=0,npts(iobj)-1
            call diam2(xobj(ioff+i),yobj(ioff+i),iclro)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine latchob(x1,y1,x0,y0,z0,dist2,k,i,iobj)
#     include "basics.inc"

      dist2 =1.e23
      x1 = x0
      y1 = y0

      if (.not.ifobjg(iobj)) return

      ioff=1
      do i=1,iobj-1
         ioff=ioff+npts(i)
      enddo

      k = 0
      do i=ioff,ioff+npts(iobj)
         k = k+1
         dist=(x0-xobj(i))**2+(y0-yobj(i))**2
         if (dist.lt.dist2) then
            x1=xobj(i)
            y1=yobj(i)
            dist2=dist
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setzoom
#     include "basics.inc"
      LOGICAL IFTMP
C
      IFTMP =IFGRID
      IFGRID=.FALSE.
      CALL PRS('Input with mouse 2 points on edge of new frame.$')
      CALL PRS('Menu area for keybd input.$')
c     CALL PRS('Menu area to restore to original size.$')
C
      CALL MOUSE(XMSE1,YMSE1,BUTTON)
   10 CONTINUE
      IF (XSCR(XMSE1).LT.1.0) THEN
         CALL COLOR(1)
         CALL MOVESC(XSCR(XMSE1)-.015,YSCR(YMSE1))
         CALL DRAWSC(XSCR(XMSE1)+.015,YSCR(YMSE1))
         CALL MOVESC(XSCR(XMSE1),YSCR(YMSE1)-.015)
         CALL DRAWSC(XSCR(XMSE1),YSCR(YMSE1)+.015)
         CALL MOUSE(XMSE2,YMSE2,BUTTON)
         IF (XSCR(XMSE2).LT.1.0) THEN
            CALL COLOR(1)
            CALL MOVESC(XSCR(XMSE2)-.015,YSCR(YMSE2))
            CALL DRAWSC(XSCR(XMSE2)+.015,YSCR(YMSE2))
            CALL MOVESC(XSCR(XMSE2),YSCR(YMSE2)-.015)
            CALL DRAWSC(XSCR(XMSE2),YSCR(YMSE2)+.015)
C           draw a box around frame area
            DELTX=ABS(XMSE1-XMSE2)
            DELTY=ABS(YMSE1-YMSE2)
            DELT2=MAX(DELTX,DELTY)/2.0
            XZ0=(XMSE1+XMSE2)/2.0
            YZ0=(YMSE1+YMSE2)/2.0
            CALL MOVEC(XZ0-DELT2,YZ0-DELT2)
            CALL DRAWC(XZ0+DELT2,YZ0-DELT2)
            CALL DRAWC(XZ0+DELT2,YZ0+DELT2)
            CALL DRAWC(XZ0-DELT2,YZ0+DELT2)
            CALL DRAWC(XZ0-DELT2,YZ0-DELT2)
C
c           CALL PRS('Click in menu area if OK, or re-enter points.$')
            CALL PRS('Click in menu area for keybd or re-enter pts.$')
            CALL MOUSE(XMSE3,YMSE3,BUTTON)
            IF (XSCR(XMSE3).GE.1.0) THEN
               IFZOOM=.TRUE.
               DELTX=ABS(XMSE1-XMSE2)
               DELTY=ABS(YMSE1-YMSE2)
               DELT =MAX(DELTX,DELTY)
               XFAC =DELT
               YFAC =XFAC
               write(6,*) 'z:xfyf',xfac,yfac
               XZERO=(XMSE1+XMSE2)/2 - XFAC/2
               YZERO=(YMSE1+YMSE2)/2 - YFAC/2
C              reset clipping window
               WFR=1.3
               WT = YPHY(.995)
               WB = YPHY(.005)
               WL = XPHY(.005)
               WR = XPHY(.995)
               IF (IFPOST.AND.IFFULL) THEN
                  WR = XPHY(WFR)
c                 XFAC=XFAC/WFR
c                 YFAC=YFAC/WFR
c                 XZERO=WFR*( (XMSE1+XMSE2)/2 - XFAC/2)
c                 YZERO=WFR*( (YMSE1+YMSE2)/2 - YFAC/2)
               ENDIF
c              write(6,*) 'zoom:',xfac,xzero,yzero,wt,wb,wl,wr
c              write(7,*) 'zoom:',xfac,xzero,yzero,wt,wb,wl,wr
            ELSE
               CALL COLOR(0)
               CALL MOVESC(XSCR(XMSE1)-.015,YSCR(YMSE1))
               CALL DRAWSC(XSCR(XMSE1)+.015,YSCR(YMSE1))
               CALL MOVESC(XSCR(XMSE1),YSCR(YMSE1)-.015)
               CALL DRAWSC(XSCR(XMSE1),YSCR(YMSE1)+.015)
               CALL MOVESC(XSCR(XMSE2)-.015,YSCR(YMSE2))
               CALL DRAWSC(XSCR(XMSE2)+.015,YSCR(YMSE2))
               CALL MOVESC(XSCR(XMSE2),YSCR(YMSE2)-.015)
               CALL DRAWSC(XSCR(XMSE2),YSCR(YMSE2)+.015)
C              un-draw a box around frame area
               CALL MOVEC(XZ0-DELT2,YZ0-DELT2)
               CALL DRAWC(XZ0+DELT2,YZ0-DELT2)
               CALL DRAWC(XZ0+DELT2,YZ0+DELT2)
               CALL DRAWC(XZ0-DELT2,YZ0+DELT2)
               CALL DRAWC(XZ0-DELT2,YZ0-DELT2)
               XMSE1=XMSE3
               YMSE1=YMSE3
               CALL PRS(
     $         'Click 2nd point on edge of new frame, or in$')
               CALL PRS(
     $         'menu area to abort and restore to original size.$')
               GOTO 10
            ENDIF
         ELSE
            IFZOOM=.FALSE.
            CALL RESCAL
         ENDIF
      ELSE
c
c        1st click was in menu area... therefore use keybd input
c
         CALL PRS('Input two X,Y pairs (x1,y1,x2,y2)$')
         CALL PRS(' to bracket zoom area.$')
         CALL PRS('Input all zeros to restore screen.$')
         CALL PRS('Input 1 1 * * to cancel operation.$')
         CALL PRS('Input 2 2 s * to to set scale zoom by s.$')
         CALL PRS('Input 2 -2 s * to to set scale y zoom by s.$')
         CALL PRS('Input 3 3 sx sy to to scale x by sx, y by sy.$')
         CALL RERRRR(x1,y1,x2,y2)
c  
         if (x1.eq.x2.and.y1.eq.y2) then
            if (x1.eq.0.0.and.y1.eq.0.0) then
c             Restore
               IFZOOM=.FALSE.
               CALL RESCAL
            elseif (x1.eq.1.0.and.y1.eq.1.0) then
               CALL PRS('Aborting zoom request.$')
            endif
         else if (x1.eq.2. .and. y1.eq.2.) then
c           Scale current coordinates by scale
            scale = 1./x2
            IFZOOM=.TRUE.
            delt =xfac*scale
            xzero=xzero + xfac/2 - delt/2
            yzero=yzero + yfac/2 - delt/2
            xfac =delt
            yfac =xfac
            write(6,*) 'z:xfyf',xfac,yfac
c
c           reset clipping window
            WFR=1.3
            WT = YPHY(.995)
            WB = YPHY(.005)
            WL = XPHY(.005)
            WR = XPHY(.995)
            IF (IFPOST.AND.IFFULL) THEN
               WR = XPHY(WFR)
            ENDIF
         else if (x1.eq.2. .and. y1.eq.-2.) then
c           Scale current y coordinates by scale
            scale = 1./x2
            IFZOOM=.TRUE.
            delt =yfac*scale
            yzero=yzero + yfac/2 - delt/2
            yfac =delt
            write(6,*) 'z:xfyf',xfac,yfac
c
c           reset clipping window
            WFR=1.3
            WT = YPHY(.995)
            WB = YPHY(.005)
            WL = XPHY(.005)
            WR = XPHY(.995)
            IF (IFPOST.AND.IFFULL) THEN
               WR = XPHY(WFR)
            ENDIF
         else if (x1.eq.3. .and. y1.eq.3.) then
c           Scale x coordinates by x2, y by y2
            scalx = 1./x2
            scaly = 1./y2
            IFZOOM=.TRUE.
            deltx =xfac*scalx
            delty =yfac*scaly
            xzero=xzero + xfac/2 - deltx/2
            yzero=yzero + yfac/2 - delty/2
            xfac =deltx
            yfac =delty
            write(6,*) 'z:xfyf',xfac,yfac
c
c           reset clipping window
            WFR=1.3
            WT = YPHY(.995)
            WB = YPHY(.005)
            WL = XPHY(.005)
            WR = XPHY(.995)
            IF (IFPOST.AND.IFFULL) THEN
               WR = XPHY(WFR)
            ENDIF
         else
c
c        Got two (x,y) pairs
c
            xmse1 = x1
            ymse1 = y1
            xmse2 = x2
            ymse2 = y2
c
            CALL COLOR(1)
            CALL MOVESC(XSCR(XMSE1)-.015,YSCR(YMSE1))
            CALL DRAWSC(XSCR(XMSE1)+.015,YSCR(YMSE1))
            CALL MOVESC(XSCR(XMSE1),YSCR(YMSE1)-.015)
            CALL DRAWSC(XSCR(XMSE1),YSCR(YMSE1)+.015)
c           CALL MOUSE(XMSE2,YMSE2,BUTTON)
            CALL COLOR(1)
            CALL MOVESC(XSCR(XMSE2)-.015,YSCR(YMSE2))
            CALL DRAWSC(XSCR(XMSE2)+.015,YSCR(YMSE2))
            CALL MOVESC(XSCR(XMSE2),YSCR(YMSE2)-.015)
            CALL DRAWSC(XSCR(XMSE2),YSCR(YMSE2)+.015)
C           draw a box around frame area
            DELTX=ABS(XMSE1-XMSE2)
            DELTY=ABS(YMSE1-YMSE2)
            DELT2=MAX(DELTX,DELTY)/2.0
            XZ0=(XMSE1+XMSE2)/2.0
            YZ0=(YMSE1+YMSE2)/2.0
            CALL MOVEC(XZ0-DELT2,YZ0-DELT2)
            CALL DRAWC(XZ0+DELT2,YZ0-DELT2)
            CALL DRAWC(XZ0+DELT2,YZ0+DELT2)
            CALL DRAWC(XZ0-DELT2,YZ0+DELT2)
            CALL DRAWC(XZ0-DELT2,YZ0-DELT2)
C
c           CALL PRS('Click in menu area if OK, or re-enter points.$')
c           CALL MOUSE(XMSE3,YMSE3,BUTTON)
cccccc      IF (XSCR(XMSE3).GE.1.0) THEN
               IFZOOM=.TRUE.
               DELTX=ABS(XMSE1-XMSE2)
               DELTY=ABS(YMSE1-YMSE2)
               DELT =MAX(DELTX,DELTY)
               XFAC =DELT
               YFAC =XFAC
               write(6,*) 'z:xfyf',xfac,yfac
               XZERO=(XMSE1+XMSE2)/2 - XFAC/2
               YZERO=(YMSE1+YMSE2)/2 - YFAC/2
C              reset clipping window
               WFR=1.3
               WT = YPHY(.995)
               WB = YPHY(.005)
               WL = XPHY(.005)
               WR = XPHY(.995)
               IF (IFPOST.AND.IFFULL) THEN
                  WR = XPHY(WFR)
               ENDIF
ccccccc     endif
         endif
c.old    IFZOOM=.FALSE.
c.old    CALL RESCAL
      ENDIF
      IFGRID=IFTMP
      return
      END
c-----------------------------------------------------------------------
      subroutine setrans
#     include "basics.inc"
      LOGICAL IFTMP
C
      IFTMP =IFGRID
      IFGRID=.TRUE.
      CALL DRGRID
      CALL PRS('Input with mouse 2 points on edge of new frame.$')
      CALL PRS('Menu area to exit.$')
C
      CALL MOUSE(XMSE1,YMSE1,BUTTON)
   10 CONTINUE
      IF (XSCR(XMSE1).LT.1.0) THEN
         CALL COLOR(1)
         CALL MOVESC(XSCR(XMSE1)-.015,YSCR(YMSE1))
         CALL DRAWSC(XSCR(XMSE1)+.015,YSCR(YMSE1))
         CALL MOVESC(XSCR(XMSE1),YSCR(YMSE1)-.015)
         CALL DRAWSC(XSCR(XMSE1),YSCR(YMSE1)+.015)
         CALL MOUSE(XMSE2,YMSE2,BUTTON)
         IF (XSCR(XMSE2).LT.1.0) THEN
            CALL COLOR(1)
            CALL MOVESC(XSCR(XMSE2)-.015,YSCR(YMSE2))
            CALL DRAWSC(XSCR(XMSE2)+.015,YSCR(YMSE2))
            CALL MOVESC(XSCR(XMSE2),YSCR(YMSE2)-.015)
            CALL DRAWSC(XSCR(XMSE2),YSCR(YMSE2)+.015)
C
            CALL PRS('Click in menu area if OK, or re-enter points.$')
            CALL MOUSE(XMSE3,YMSE3,BUTTON)
            IF (XSCR(XMSE3).GE.1.0) THEN
               IFZOOM=.TRUE.
               DELTX=ABS(XMSE1-XMSE2)
               DELTY=ABS(YMSE1-YMSE2)
               XZERO=XZERO+XMSE1-XMSE2
               YZERO=YZERO+YMSE1-YMSE2
C              reset clipping window
               WFR=1.3
               WT = YPHY(.995)
               WB = YPHY(.005)
               WL = XPHY(.005)
               WR = XPHY(.995)
               IF (IFFULL) THEN
                  WR = XPHY(WFR)
               ENDIF
            ELSE
               CALL COLOR(0)
               CALL MOVESC(XSCR(XMSE1)-.015,YSCR(YMSE1))
               CALL DRAWSC(XSCR(XMSE1)+.015,YSCR(YMSE1))
               CALL MOVESC(XSCR(XMSE1),YSCR(YMSE1)-.015)
               CALL DRAWSC(XSCR(XMSE1),YSCR(YMSE1)+.015)
               CALL MOVESC(XSCR(XMSE2)-.015,YSCR(YMSE2))
               CALL DRAWSC(XSCR(XMSE2)+.015,YSCR(YMSE2))
               CALL MOVESC(XSCR(XMSE2),YSCR(YMSE2)-.015)
               CALL DRAWSC(XSCR(XMSE2),YSCR(YMSE2)+.015)
               XMSE1=XMSE3
               YMSE1=YMSE3
               CALL PRS(
     $         'Click 2nd point in the domain, or in$')
               CALL PRS(
     $         'menu area to return.$')
               GOTO 10
            ENDIF
         ENDIF
      ENDIF
      IFGRID=IFTMP
      return
      END
c-----------------------------------------------------------------------
      subroutine rescal
#     include "basics.inc"
c
      if (ifzoom) return
c
      XMAX = -1.0E20
      YMAX = -1.0E20
      XMIN =  1.0E20
      YMIN =  1.0E20
      nc = 4
      if (if3d) nc = 8
      DO 50 IEL=1,NEL
         DO 50 IC=1,nc
         xmax = max(xmax,xiso(x(ic,iel),y(ic,iel),z(ic,iel)) )
         xmin = min(xmin,xiso(x(ic,iel),y(ic,iel),z(ic,iel)) )
         ymax = max(ymax,yiso(x(ic,iel),y(ic,iel),z(ic,iel)) )
         ymin = min(ymin,yiso(x(ic,iel),y(ic,iel),z(ic,iel)) )
50    CONTINUE
      IF(IFPOST)THEN
         IF(IF3D)THEN
C           Fill screen with object
            XFAC = MAX( (XMAX-XMIN) , (YMAX-YMIN)) *1.3
            YFAC = XFAC
C           Center on Screen
            XZERO = ((XMAX+XMIN)/2) - XFAC/2
            YZERO = ((YMAX+YMIN)/2) - YFAC/2
            IF (IFFULL) THEN
C              Fill screen with object
               DELTX=(XMAX-XMIN)/WFR
               XFAC = MAX( DELTX , (YMAX-YMIN)) *1.025
               YFAC = XFAC
C              Center on Screen
               XZERO = WFR*( ((XMAX+XMIN)/2) - XFAC/2)
C
C              quickie fix 12-5-91 pff
C
               XFAC=XFACO*5.0
               YFAC=YFACO*5.0
               XZERO=XZEROO
               YZERO=YZEROO
               YZERO =       ((YMAX+YMIN)/2) - YFAC/2
               write(6,*) 'r1xfyf',xfac,yfac
C
C
C
            ENDIF
         ELSE
C        For 2-d, leave xfac,yfac alone if THETA=90 and PHI=-90
            IF (THETA.EQ.90.0.AND.PHI.EQ.-90.0.AND.IROT.EQ.0) THEN
               XFAC=XFACO
               YFAC=YFACO
               XZERO=XZEROO
               YZERO=YZEROO
            ELSE
C              Fill screen with object
               XFAC = MAX( (XMAX-XMIN) , (YMAX-YMIN)) *1.3
               YFAC = XFAC
C              Center on Screen
               XZERO = ((XMAX+XMIN)/2) - XFAC/2
               YZERO = ((YMAX+YMIN)/2) - YFAC/2
            ENDIF
         ENDIF
      ELSE
         XFAC=XFACO
         YFAC=YFACO
c              XFAC=XFACO*5.0
c              YFAC=YFACO*5.0
c              write(6,*) 'r2xfyf',xfac,yfac
         XZERO=XZEROO
         YZERO=YZEROO
C        Put up in right upper corner  (obsolete, draw icon seperately)
c        XFAC = MAX( (XMAX-XMIN) , (YMAX-YMIN)) * 10.0
c        YFAC = XFAC
C        Center on Screen
c        XZERO = ((XMAX+XMIN)/2) - XFAC*.9
c        YZERO = ((YMAX+YMIN)/2) - YFAC*.9
      ENDIF
C     Reset clipping window
      WFR= 1.3
      WT = YPHY(.995)
      WB = YPHY(.005)
      WL = XPHY(.005)
      WR = XPHY(.995)
      IF (IFFULL) WR = XPHY(WFR)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine ellgen(x,y,n,a,b)
c
c     Generate a set of points on an ellipse (a,b)
c
      real x(0:1),y(0:1)
      real t(0:800)
c
      call prs('Input a,b (0 or <0 = quit) :$')
      call rerr(a,b)
      if (a.le.0 .or. b.le.0) then
         call prs('returning from ellgen, nothing done$')
         return
      endif
c
c     Recursively bisect the ellipse to get constant relative curvature
c
      b2i = 1.0/(b*b)
      a2i = 1.0/(a*a)
      x(0) = -a
      y(0) = 0.0
c
      dx  = .1
      dx4 = dx/4.
      do i=1,12
         x(i) = x(i-1)+dx4
         y(i) = b*sqrt(1.0-a2i*x(i)**2)
c        write(36,*) i,x(i),y(i)
      enddo
      do i=11,1000
         x(i) = x(i-1)+dx
         if (x(i).ge.0.0) goto 10
         y(i) = b*sqrt(1.0-a2i*x(i)**2)
c        write(36,*) i,x(i),y(i)
      enddo
   10 continue
      x(i) = 0.0
      y(i) =  b
c
      nnew = i
      do irecur = 1,7
         j = 0
         do i = 1,nnew
            call genpair(t(j),x(i-1),y(i-1),x(i),y(i),a,b)
c        write(39,6) j,t(j),t(j+1),y(i-1),y(i),x(i-1),x(i),i,irecur,a,b
            j=j+2
         enddo
    6    format(i3,6f9.5,2i5,2f5.2)
         nold = nnew
         nnew = j
         t(nnew) = y(nold)
c        write(38,*) 'this is nnew',nnew,nold,irecur,j
         call copy(y,t,nnew+1)
         do j=1,nnew-1
            disc = 1.0 - b2i*y(j)**2
            if (disc.le.0.0) then
               x(j) = 0.0
            else
               x(j) = -a*sqrt(1.0-b2i*y(j)**2)
            endif
c           write(38,8) j,x(j),y(j),b2i,disc,t(j)
    8    format(i5,5g13.5)
         enddo
         x(0) =    -a
         y(0) =    0.0
         x(nnew) = 0.0
         y(nnew) =  b
            j=0
c           write(38,8) j,x(j),y(j),b2i,disc
            j=nnew
c           write(38,8) j,x(j),y(j),b2i,disc
c
         if (nnew.gt.200) goto 20
      enddo
   20 continue
c
c     Now copy to  rh plane
c
      n = nnew
c     n = 2*nnew
c     do i=0,nnew-1
c        write(6,*) i,x(i),y(i)
c        x(n-i) = -x(i)
c        y(n-i) =  y(i)
c     enddo
c
      n = n+1
      write(6,*) 'Num obj segments:',n
      return
      end
c-----------------------------------------------------------------------
      subroutine genpair(ynew,x0,y0,x1,y1,a,b)
c
c     Generate new endpoints for bisected secants
c
      real ynew(0:1)
c
      ynew(0) = y0
c
c     Slope of perpendicular bisector
c
      if (y1.ne.y0) then
         d = -(x1-x0)/(y1-y0)
      else
         d = -9.e9
      endif
      e = 0.5*(y0+y1 - d*(x0+x1))
c
      ab2 = (a/b)**2
      AA = (1.+ab2*d**2)
      BB = 2.*d*e*ab2
      CC = ab2*(e**2 - (b**2))
c
c     write(60,*) 'a,b',a,b,aa,bb,cc
c     write(60,*) 'x,y',x0,y0,x1,y1
c     write(60,*) 'd,e',d,e,ab2
      if (BB.gt.0) then
        xt = 0.5*(-BB-sqrt(BB**2-4.*AA*CC))/AA
      else
        xt = 2.0*CC/(-BB+sqrt(BB**2-4.*AA*CC))
      endif
      xt = -abs(xt)
      ynew(1) = abs(d*xt+e)
c     write(60,*) 'x,t',xt,ynew(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine redraw_mesh
#     include "basics.inc"
c
      call refresh
      call drmenu('NOCOVER')
      call drgrid

      nelcap_l = min(nelcap,nel)

      do ie=1,nelcap_l
         call drawel(ie)
      enddo

      call drwobjs         !  draw object on top of mesh, for visibility

      if (if3d) then       !  redraw all the isometric elements.
         call sortel
         do i=1,nelcap_l
            call drawis(isrt(i))
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine redraw_mesh_small
#     include "basics.inc"
c
      call refresh
      call drmenu('NOCOVER')
      call drgrid

      nelcap_l = min(nelcap,nel)

      do ie=1,nelcap_l
         call drawel(ie)
      enddo

      call drwobjs         !  draw object on top of mesh, for visibility


C     Now redraw all the isometric elements.
      call sortel
      do i=1,nelcap_l
         call drawis(isrt(i))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine setobj
#     include "basics.inc"

      call getobjs
      return

    1 item(1) = 'MAIN MENU'
      item(2) = 'BEZIER SET'
      nchoic  = 2
      call menu(xmouse,ymouse,button,'SET OBJECT')
 

      if (choice.eq.'MAIN MENU') then
         call redraw_mesh
         return
      endif

      if (choice.eq.'BEZIER SET') call bezier

      goto 1

      return
      end
c-----------------------------------------------------------------------
      subroutine color_test
#     include "basics.inc"

      parameter (mm=100)
      common /ctmp0/ xbz(0:mm),ybz(0:mm)

      do iclr=0,15
         write(6,*) iclr,' iclr before'
         call mouse_diam(xmouse,ymouse,button,iclr)
         write(6,*) iclr,' iclr after'
      enddo

c     0=black
c     1=white
c     2=red
c     3=yellow
c     4=blank?
c     5=blue
c     6=red
c     7=yellow
c     8=pale green
c     9=bright green
c     0=light blue
c    11=dark blue
c    12=magenta
c    13=red
c    14=gray
c    15= ???

      return
      end
c-----------------------------------------------------------------------
      subroutine bezier
#     include "basics.inc"

      i0   = 1
      nbez = nbez+1 
      if (nbez.gt.1) i0=ibez(nbez-1)

      call prs('Click 4 control points or menu area to end.$')

      ib = i0
      call mouse_diam(xybez(1,ib),xybez(2,ib),button,10) ! light blue
      if (xybez(1,ib).gt.xphy(1.0)) then  ! menu area
         nbez = nbez-1
         return
      endif

      do i=1,lbez
         i0 = ib

         do j=1,3  ! Get next 3 Bezier points

            ib = ib+1
            call mouse_diam(xybez(1,ib),xybez(2,ib),button,10) ! light blue

            if (xybez(1,ib).gt.xphy(1.0)) then  ! menu area
               ibez(nbez) = i0+1
               if (i.eq.1) nbez = nbez-1        ! aborting
               return
            endif

            if (ib.gt.lxybez) then
               call prsii
     $         ('Exceeded max Bezier points.  Returning.$',ib,lxybez)
               nbez = nbez-1
               return
            endif

         enddo

         j1 = ib
         j0 = ib-3
         call dr_bezier(j0,j1)

         call prs('Click 3 more points to continue, menu to end.$')

      enddo

c   1 item(1) = 'UP MENU'
c     item(2) = 'DELETE SPLINE'
c     nchoic  = 2
c     call menu(xmouse,ymouse,button,'BEZIER')

      return
      end
c-----------------------------------------------------------------------
      subroutine mouse_diam(xmouse,ymouse,button,iclr)
#     include "basics.inc"

      call mouse(xmouse,ymouse,button)
      if (xmouse.lt.xphy(1.0)) call diam2(xmouse,ymouse,iclr)

      return
      end
c-----------------------------------------------------------------------
      subroutine dr_bezier(j0,j1)
#     include "basics.inc"

      parameter (mm=50)
      common /ctmp0/ xbz(0:mm),ybz(0:mm)

      ds = 1./mm

      do i0=j0,j1-1,3

         x0 = xybez(1,i0  )
         x1 = xybez(1,i0+1)
         x2 = xybez(1,i0+2)
         x3 = xybez(1,i0+3)

         y0 = xybez(2,i0  )
         y1 = xybez(2,i0+1)
         y2 = xybez(2,i0+2)
         y3 = xybez(2,i0+3)

         do i=0,mm
            s1 = i*ds
            t1 = (mm-i)*ds
            c0 = t1*t1*t1
            c1 = 3*t1*t1*s1
            c2 = 3*t1*s1*s1
            c3 = s1*s1*s1
            xbz(i) = c0*x0+c1*x1+c2*x2+c3*x3
            ybz(i) = c0*y0+c1*y1+c2*y2+c3*y3
         enddo

         call color(8)                 ! pale green line
         call vdraw (xbz,ybz,mm+1)
         call diam2(x0,y0,7)           ! yellow endpoints
         call diam2(x3,y3,7)

         write(6,*) j0,x0,y0,' x0'
         write(6,*) j1,x3,y3,' x3'

      enddo

      return
      end
c-----------------------------------------------------------------------
