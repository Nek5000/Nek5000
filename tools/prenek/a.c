c-----------------------------------------------------------------------
      subroutine sphmesh
      include 'basics.inc'
      common /ctmp0/ sphctr(3),xcs(4,24),ycs(4,24),zcs(4,24)
      character*1 SHELL,HEMI,YESNO
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      logical ifpshp
c
      real xlat0(3)
C
C     Build either spherical or hemi-spherical meshes.
C
C     Modified to include stretching in X for prolate spheroird. pff 7/30/92
C
C
c     CALL PRS('A 6 or 24 element mesh?$')
c     CALL REI(MESH)
      MESH = 24
C
      CALL PRS
     $('Hemi, whole sphere, hex-pkd hemi, tet/dia/lat? (H/W/X/T/D/L):$')
      CALL RES(HEMI,1)
      CALL CAPIT(HEMI,1)
c
      if (hemi.eq.'t'.or.hemi.eq.'T') then
c
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_tet     (dlat,sphrad,xlat0)
         return
c
      elseif (hemi.eq.'d'.or.hemi.eq.'D') then
c
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_dia     (dlat,sphrad,xlat0)
         return
c
      elseif (hemi.eq.'l'.or.hemi.eq.'L') then
c
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_lat     (dlat,sphrad,xlat0)
         return
c
      endif
C
      CALL PRS('Enter the (X,Y,Z) coordinates of the center:$')
      CALL RERRR(SPHCTR(1),SPHCTR(2),SPHCTR(3))
c
      IF (HEMI.EQ.'H') THEN
         NELSPH=12
         IF (MESH.EQ.6) NELSPH=5
      ELSEIF (HEMI.EQ.'X') THEN
         NELSPH=12
         IF (MESH.EQ.6) NELSPH=5
      ELSE
c        CALL PRS('A 6 or 24 element mesh?$')
c        CALL REI(MESH)
         mesh = 24
         NELSPH=MESH
      ENDIF
      NLSPH4=4*NELSPH
C
C     PROLATE SPHEROID QUERY:
C
      IFpSPH=.FALSE.
      RATIO=1.0
      CALL PRS  ('Prolate spheroid? (Y/N):$')
      CALL RES  (YESNO,1)
      CALL CAPIT(YESNO,1)
      IF (YESNO.EQ.'Y') THEN
         IFpSPH=.TRUE.
         CALL PRS('Enter ratio:$')
         CALL RER (RATIO)
         CALL PRSR('The ratio is:$',RATIO)
      ENDIF
C
C---------------------------------------------------------------------------
C     Begin building shell sequence, working from the inner most to outer.
C---------------------------------------------------------------------------
C
      DO 6000 ISHLL=1,1000
C
         CALL PRS(
     $  'Enter S or C for a spherical or cartesian layer (E=exit):$')
         CALL RES(SHELL,1)
         CALL CAPIT(SHELL,1)
         IF (SHELL.EQ.'E') GOTO 9000
C
         IF (SHELL.EQ.'S') THEN
            CALL PRS('Enter radius:$')
            CALL RER(RADIUS)
            CALL SPHERE(XCS,YCS,ZCS,HEMI,MESH,RADIUS)
            CALL TRANS2(XCS,YCS,ZCS,SPHCTR,NLSPH4)
         ELSE
            CALL PRS(
     $     'Enter minimum distance from center to edge of the box:$')
            CALL RER(RADIUS)
            CALL CRTBOX(XCS,YCS,ZCS,HEMI,MESH,RADIUS)
            CALL TRANS2(XCS,YCS,ZCS,SPHCTR,NLSPH4)
         ENDIF
C
C        Update the elements
C
         IF (ISHLL.GT.1) NEL=NEL+NELSPH
         DO 1000 IE=1,NELSPH
C
            IEL=NEL+IE
            DO 101 I=1,4
               X(IEL,I)=XCS(I,IE)*RATIO
               Y(IEL,I)=YCS(I,IE)
               Z(IEL,I)=ZCS(I,IE)
  101       CONTINUE

            IF (SHELL.EQ.'S') THEN
               CCURVE(5,IEL)='s'
               IF (IFpSPH) CCURVE(5,IEL)='p'
               CURVE(1,5,IEL)=SPHCTR(1)
               CURVE(2,5,IEL)=SPHCTR(2)
               CURVE(3,5,IEL)=SPHCTR(3)
               CURVE(4,5,IEL)=RADIUS
               CURVE(5,5,IEL)=RATIO
            ENDIF
            IF (ISHLL.GT.1) THEN
               IEL=NEL-NELSPH+IE
               DO 201 I=1,4
                  J=I+4
                  X(IEL,J)=XCS(I,IE)*RATIO
                  Y(IEL,J)=YCS(I,IE)
                  Z(IEL,J)=ZCS(I,IE)
  201          CONTINUE
               IF (SHELL.EQ.'S') THEN
                  CCURVE(6,IEL)='s'
                  IF (IFpSPH) CCURVE(6,IEL)='p'
                  CURVE(1,6,IEL)=SPHCTR(1)
                  CURVE(2,6,IEL)=SPHCTR(2)
                  CURVE(3,6,IEL)=SPHCTR(3)
                  CURVE(4,6,IEL)=RADIUS
                  CURVE(5,6,IEL)=RATIO
               ENDIF
               NUMAPT(IEL)=ISHLL-1
               ilet = mod1(ie,52)
               LETAPT(IEL)=alphabet(ilet)
            ENDIF
 1000    CONTINUE
 6000 CONTINUE
C
 9000 CONTINUE

C     Recount the number of curved sides
C
      NCURVE=0
      DO 9001 IE=1,NEL
      DO 9001 IEDGE=1,8
         IF (CCURVE(IEDGE,IE).NE.' ') THEN
            NCURVE=NCURVE+1
            WRITE(6,*) 'Curve:',IE,IEDGE,CCURVE(IEDGE,IE)
         ENDIF
 9001 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine sphere(xcs,ycs,zcs,hemi,mesh,radius)
      DIMENSION XCS(4,24),YCS(4,24),ZCS(4,24)
      character*1 HEMI
C
      ONE=1.0
      PI2=2.0*ATAN(ONE)
      PI =2.0*PI2
      RAD2=RADIUS/SQRT(2.0)
      RAD3=RADIUS/SQRT(3.0)
C
      IF (MESH.ne.24) THEN
         call prs('Sorry, no cubic gnomonics at this time.$')
         mesh=24
      endif
c
      IF (MESH.EQ.24) THEN
C
C        Form octant first, then replicate.
C
         XCS(1,1)=RADIUS
         YCS(1,1)=0.0
         ZCS(1,1)=0.0
C
         XCS(2,1)=RAD2
         YCS(2,1)=RAD2
         ZCS(2,1)=0.0
C
         XCS(3,1)=RAD3
         YCS(3,1)=RAD3
         ZCS(3,1)=RAD3
C
         XCS(4,1)=RAD2
         YCS(4,1)=0.0
         ZCS(4,1)=RAD2
C
         XCS(1,2)=RAD2
         YCS(1,2)=RAD2
         ZCS(1,2)=0.0
C
         XCS(2,2)=0.0
         YCS(2,2)=RADIUS
         ZCS(2,2)=0.0
C
         XCS(3,2)=0.0
         YCS(3,2)=RAD2
         ZCS(3,2)=RAD2
C
         XCS(4,2)=RAD3
         YCS(4,2)=RAD3
         ZCS(4,2)=RAD3
C
         XCS(1,3)=0.0
         YCS(1,3)=0.0
         ZCS(1,3)=RADIUS
C
         XCS(2,3)=RAD2
         YCS(2,3)=0.0
         ZCS(2,3)=RAD2
C
         XCS(3,3)=RAD3
         YCS(3,3)=RAD3
         ZCS(3,3)=RAD3
C
         XCS(4,3)=0.0
         YCS(4,3)=RAD2
         ZCS(4,3)=RAD2
C
C        Replicate octant
C
         CALL COPY(XCS(1,4),XCS(1,1),12)
         CALL COPY(YCS(1,4),YCS(1,1),12)
         CALL COPY(ZCS(1,4),ZCS(1,1),12)
         CALL ROTAT2(XCS(1,4),YCS(1,4),ZCS(1,4),12,'Z',PI2)
         IF (HEMI.EQ.'X') THEN
c
c           Modify vertices to yield hexagonal box
c
            rad32 = 0.5*radius*sqrt(3.)
            rad37 = radius*sqrt(3./7.)
            rad17 = radius*sqrt(1./7.)
            rad5  = 0.5*radius
c
            XCS(2,1)=rad5
            YCS(2,1)=RAD32
            ZCS(2,1)=0.0
C
            XCS(3,1)=RAD17
            YCS(3,1)=RAD37
            ZCS(3,1)=RAD37
C
            XCS(1,2)=rad5
            YCS(1,2)=RAD32
            ZCS(1,2)=0.0
C
            XCS(4,2)=RAD17
            YCS(4,2)=RAD37
            ZCS(4,2)=RAD37
C
            XCS(3,3)=RAD17
            YCS(3,3)=RAD37
            ZCS(3,3)=RAD37
c
            XCS(2,4)= -rad5
            YCS(2,4)=  rad32
            ZCS(2,4)=  0.0
C
            XCS(3,4)= -RAD17
            YCS(3,4)=  RAD37
            ZCS(3,4)=  RAD37
C
            XCS(1,5)= -rad5
            YCS(1,5)=  rad32
            ZCS(1,5)=  0.0
C
            XCS(4,5)= -RAD17
            YCS(4,5)=  RAD37
            ZCS(4,5)=  RAD37
C
            XCS(3,6)= -RAD17
            YCS(3,6)=  RAD37
            ZCS(3,6)=  RAD37
         ENDIF
C
C
C        Replicate quadrant
C
         CALL COPY(XCS(1,7),XCS(1,1),24)
         CALL COPY(YCS(1,7),YCS(1,1),24)
         CALL COPY(ZCS(1,7),ZCS(1,1),24)
         CALL ROTAT2(XCS(1,7),YCS(1,7),ZCS(1,7),24,'Z',PI)
         CALL ROUNDER(XCS,48)
         CALL ROUNDER(YCS,48)
         CALL ROUNDER(ZCS,48)
         IF (HEMI.EQ.'W') THEN
C
C           Replicate hemisphere
C
            CALL COPY(XCS(1,13),XCS(1,1),48)
            CALL COPY(YCS(1,13),YCS(1,1),48)
            CALL COPY(ZCS(1,13),ZCS(1,1),48)
            CALL ROTAT2(XCS(1,13),YCS(1,13),ZCS(1,13),48,'X',PI)
            CALL ROUNDER(XCS(1,13),48)
            CALL ROUNDER(YCS(1,13),48)
            CALL ROUNDER(ZCS(1,13),48)
         ENDIF
c     ELSE
C
C        We mesh the 6 element configuration here
C
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine crtbox(xcs,ycs,zcs,hemi,mesh,radius)
      DIMENSION XCS(4,24),YCS(4,24),ZCS(4,24)
      character*1 HEMI
C
      ONE=1.0
      PI2=2.0*ATAN(ONE)
      PI =2.0*PI2
      RAD2=RADIUS
      RAD3=RADIUS
C
      IF (MESH.EQ.24) THEN
C
C        Form octant first, then replicate.
C
         XCS(1,1)=RADIUS
         YCS(1,1)=0.0
         ZCS(1,1)=0.0
C
         XCS(2,1)=RAD2
         YCS(2,1)=RAD2
         ZCS(2,1)=0.0
C
         XCS(3,1)=RAD3
         YCS(3,1)=RAD3
         ZCS(3,1)=RAD3
C
         XCS(4,1)=RAD2
         YCS(4,1)=0.0
         ZCS(4,1)=RAD2
C
         XCS(1,2)=RAD2
         YCS(1,2)=RAD2
         ZCS(1,2)=0.0
C
         XCS(2,2)=0.0
         YCS(2,2)=RADIUS
         ZCS(2,2)=0.0
C
         XCS(3,2)=0.0
         YCS(3,2)=RAD2
         ZCS(3,2)=RAD2
C
         XCS(4,2)=RAD3
         YCS(4,2)=RAD3
         ZCS(4,2)=RAD3
C
         XCS(1,3)=0.0
         YCS(1,3)=0.0
         ZCS(1,3)=RADIUS
C
         XCS(2,3)=RAD2
         YCS(2,3)=0.0
         ZCS(2,3)=RAD2
C
         XCS(3,3)=RAD3
         YCS(3,3)=RAD3
         ZCS(3,3)=RAD3
C
         XCS(4,3)=0.0
         YCS(4,3)=RAD2
         ZCS(4,3)=RAD2
C
C        Replicate octant
C
         CALL COPY(XCS(1,4),XCS(1,1),12)
         CALL COPY(YCS(1,4),YCS(1,1),12)
         CALL COPY(ZCS(1,4),ZCS(1,1),12)
         CALL ROTAT2(XCS(1,4),YCS(1,4),ZCS(1,4),12,'Z',PI2)
c
         IF (HEMI.EQ.'X') THEN
c
c           Modify vertices to yield hexagonal box
c
            rad13 =    radius/sqrt(3.)
            rad23 = 2.*radius/sqrt(3.)
c
            XCS(1,1)=rad23
            XCS(2,1)=rad13
            XCS(3,1)=rad13
            XCS(4,1)=rad23
C
            XCS(1,2)=rad13
            XCS(4,2)=rad13
C
            XCS(2,3)=rad23
            XCS(3,3)=rad13
c
            XCS(2,4)= -rad13
            XCS(3,4)= -rad13
C
            XCS(1,5)= -rad13
            XCS(2,5)= -rad23
            XCS(3,5)= -rad23
            XCS(4,5)= -rad13
C
            XCS(3,6)= -rad13
            XCS(4,6)= -rad23
C
         ENDIF
C
C
C        Replicate quadrant
C
         CALL COPY(XCS(1,7),XCS(1,1),24)
         CALL COPY(YCS(1,7),YCS(1,1),24)
         CALL COPY(ZCS(1,7),ZCS(1,1),24)
         CALL ROTAT2(XCS(1,7),YCS(1,7),ZCS(1,7),24,'Z',PI)
         CALL ROUNDER(XCS,48)
         CALL ROUNDER(YCS,48)
         CALL ROUNDER(ZCS,48)
         IF (HEMI.EQ.'W') THEN
C
C           Replicate hemisphere
C
            CALL COPY(XCS(1,13),XCS(1,1),48)
            CALL COPY(YCS(1,13),YCS(1,1),48)
            CALL COPY(ZCS(1,13),ZCS(1,1),48)
            CALL ROTAT2(XCS(1,13),YCS(1,13),ZCS(1,13),48,'X',PI)
            CALL ROUNDER(XCS(1,13),48)
            CALL ROUNDER(YCS(1,13),48)
            CALL ROUNDER(ZCS(1,13),48)
         ENDIF
c     ELSE
C
C        We mesh the 6 element configuration here
C
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine rotat2(x,y,z,n,dir,angle)
C
      DIMENSION X(1),Y(1),Z(1)
      character*1 DIR
      DIMENSION ROTAT(3,3)
C
      CALL RZERO(ROTAT,9)
      COSANG=COS(ANGLE)
      SINANG=SIN(ANGLE)
C
      IF (DIR.EQ.'X') THEN
         ROTAT(1,1) =   1.0
         ROTAT(2,2) =   COSANG
         ROTAT(3,3) =   COSANG
         ROTAT(2,3) = - SINANG
         ROTAT(3,2) =   SINANG
      ENDIF
C
      IF (DIR.EQ.'Y') THEN
         ROTAT(1,1) =   COSANG
         ROTAT(2,2) =   1.0
         ROTAT(3,3) =   COSANG
         ROTAT(1,3) =   SINANG
         ROTAT(3,1) = - SINANG
      ENDIF
C
      IF (DIR.EQ.'Z') THEN
         ROTAT(1,1) =   COSANG
         ROTAT(2,2) =   COSANG
         ROTAT(3,3) =   1.0
         ROTAT(1,2) = - SINANG
         ROTAT(2,1) =   SINANG
      ENDIF
C
      DO 100 I=1,N
         XP=ROTAT(1,1)*X(I)+ROTAT(1,2)*Y(I)+ROTAT(1,3)*Z(I)
         YP=ROTAT(2,1)*X(I)+ROTAT(2,2)*Y(I)+ROTAT(2,3)*Z(I)
         ZP=ROTAT(3,1)*X(I)+ROTAT(3,2)*Y(I)+ROTAT(3,3)*Z(I)
         X(I)=XP
         Y(I)=YP
         Z(I)=ZP
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine rounder(x,n)
      DIMENSION X(1)
C
C     Try to Round X to fractional integer - eg .05 .1 .15 - if it's within 10-6
C
      DO 100 I=1,N
         EPS=1.0E-5
         XTMP=20.0*X(I)
         XTMP2=XTMP+0.5
         ITMP=INT(XTMP2)
         XTMP2=FLOAT(ITMP)
         IF (ABS(XTMP-XTMP2).LT.EPS) X(I)=XTMP2/20.0
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine trans2(x,y,z,xyzoff,n)
      DIMENSION X(1),Y(1),Z(1)
      DIMENSION XYZOFF(3)
C
      DO 10 I=1,N
         X(I)=X(I)+XYZOFF(1)
         Y(I)=Y(I)+XYZOFF(2)
         Z(I)=Z(I)+XYZOFF(3)
   10 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine genmesh
      include 'basics.inc'
C
      common /ctmp0/ sphctr(3),xcs(4,24),ycs(4,24),zcs(4,24)
      DIMENSION CUBE(3,8),CANON(3,8),XCTR(15),YCTR(15),XYZCTR(3)
      DIMENSION ANG(2,15)
      DIMENSION IANG(2,15)
      SAVE ANG,XCTR,YCTR,CANON
      SAVE IANG
C
      DATA CANON/
     $ -1.,-1.,-1.,    1.,-1.,-1.,    1., 1.,-1.,   -1., 1.,-1.  ,
     $ -1.,-1., 1.,    1.,-1., 1.,    1., 1., 1.,   -1., 1., 1.  /
C
      DATA XCTR/ 1.0 , 3.0 , 5.0 , 7.0 , 9.0
     $         , 1.0 , 3.0 , 5.0 , 7.0 , 9.0
     $         , 1.0 , 3.0 , 5.0 , 7.0 , 9.0 /
C
      DATA YCTR/ 1.0 , 1.0 , 1.0 , 1.0 , 1.0  
     $         , 3.0 , 3.0 , 3.0 , 3.0 , 3.0  
     $         , 5.0 , 5.0 , 5.0 , 5.0 , 5.0 /
C
      DATA IANG/ 0,0 , 0,0 , 0,1 , 0,3 , 0,6
     $         , 2,0 , 0,0 , 2,0 , 0,1 , 2,2
     $         , 2,0 , 0,2 , 2,4 , 0,7 , 2,10 /
C
C
      ONE=1.0
      PI2=2.0*ATAN(ONE)
      PI =2.0*PI2
      DO 10 I=1,30
         ANG(I,1)=PI2*FLOAT(IANG(I,1))
   10 CONTINUE
C
      DO 1000 ILVEL=1,3
         ILEVEL=ILVEL
         write(6,*) 'generating level',ilevel
         ZCTR=1.0+2.0*FLOAT(ILEVEL-1)
         DO 100 III=1,15
            CALL COPY(CUBE,CANON,24)
            IF (ILEVEL.EQ.2) CALL ROTATE(CUBE,8,'Y',PI2)
            IF (ILEVEL.EQ.3) CALL ROTATE(CUBE,8,'Z',PI2)
            CALL ROTATE(CUBE,8,'Z',ANG(1,III))
            CALL ROTATE(CUBE,8,'X',ANG(2,III))
            XYZCTR(1)=XCTR(III)
            XYZCTR(2)=YCTR(III)
            XYZCTR(3)=ZCTR
            CALL TRANSL(CUBE,XYZCTR,8)
            CALL GENELE(CUBE)
         write(6,*) ilvel,iii,'element',nel,' created'
         write(8,*) ilvel,iii,'element',nel,' created'
  100    CONTINUE
 1000 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine genele(xyz)
      DIMENSION XYZ(3,8)
      include 'basics.inc'
C
      NEL=NEL+1
      DO 10 I=1,8
         X(NEL,I)=XYZ(1,I)
         Y(NEL,I)=XYZ(2,I)
         Z(NEL,I)=XYZ(3,I)
   10 CONTINUE
      NUMAPT(NEL)=ILEVEL
      LETAPT(NEL)='A'
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate(xyz,n,dir,angle)
C
      DIMENSION XYZ(3,1)
      character*1 DIR
      DIMENSION ROTAT(3,3)
C
      CALL RZERO(ROTAT,9)
      COSANG=COS(ANGLE)
      SINANG=SIN(ANGLE)
C
      IF (DIR.EQ.'X') THEN
         ROTAT(1,1) =   1.0
         ROTAT(2,2) =   COSANG
         ROTAT(3,3) =   COSANG
         ROTAT(2,3) = - SINANG
         ROTAT(3,2) =   SINANG
      ENDIF
C
      IF (DIR.EQ.'Y') THEN
         ROTAT(1,1) =   COSANG
         ROTAT(2,2) =   1.0
         ROTAT(3,3) =   COSANG
         ROTAT(1,3) =   SINANG
         ROTAT(3,1) = - SINANG
      ENDIF
C
      IF (DIR.EQ.'Z') THEN
         ROTAT(1,1) =   COSANG
         ROTAT(2,2) =   COSANG
         ROTAT(3,3) =   1.0
         ROTAT(1,2) = - SINANG
         ROTAT(2,1) =   SINANG
      ENDIF
C
      DO 100 I=1,N
         XP=ROTAT(1,1)*XYZ(1,I)+ROTAT(1,2)*XYZ(2,I)+ROTAT(1,3)*XYZ(3,I)
         YP=ROTAT(2,1)*XYZ(1,I)+ROTAT(2,2)*XYZ(2,I)+ROTAT(2,3)*XYZ(3,I)
         ZP=ROTAT(3,1)*XYZ(1,I)+ROTAT(3,2)*XYZ(2,I)+ROTAT(3,3)*XYZ(3,I)
         XYZ(1,I)=XP
         XYZ(2,I)=YP
         XYZ(3,I)=ZP
  100 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine transl(xyz,xyzoff,n)
      DIMENSION XYZ(3,1)
      DIMENSION XYZOFF(3)
C
      DO 10 I=1,N
         XYZ(1,I)=XYZ(1,I)+XYZOFF(1)
         XYZ(2,I)=XYZ(2,I)+XYZOFF(2)
         XYZ(3,I)=XYZ(3,I)+XYZOFF(3)
   10 CONTINUE
      return
      end
c-----------------------------------------------------------------------
      subroutine get_lattice_0  (d,r,x)
c
      real x(3)
c
      call prs ('Input lattice spacing, d, and radius r < d/2:$')
      call rerr(d,r)
c
      call prs ('Input coordinates for 1st sphere:$')
      call rerrr(x(1),x(2),x(3))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_lattice(x,d,x0)
c
c     Get vertices for 8-noded 3D rombus that forms a periodic
c     cell for hexagonally close packed spheres
c
c     x(:,1)=x0, and x(:,2)=x0+d
c
      real x(3,8),x0(3)
      real t(3,4)
c
      call get_vert_tet(t,d,x0) ! Base tet
c
      z1 = t(3,4)
c
      call copy(x,t,3*3)
      call copy(x(1,5),t(1,4),3)
c
      x(1,4) = x(1,3) + d
      x(2,4) = x(2,3)
      x(3,4) = x(3,3)
c
      x(1,6) = x(1,5) + d
      x(2,6) = x(2,5)
      x(3,6) = x(3,5)
c
      dy = x(2,3) - x(2,1)
c
      x(1,7) = x(1,2)
      x(2,7) = x(2,5) + dy
      x(3,7) = x(3,5)
c
      x(1,8) = x(1,7) + d
      x(2,8) = x(2,7)
      x(3,8) = x(3,7)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_dia(x,d,x0)
c
c     Get vertices for a diamond with lattice spacing d
c
c     Equilateral diamond is oriented as follows:
c
c
c
c         Y                 F
c                          / \                 
c         ^            c ./ . \. b
c         |             ./     \.              
c         |            D---------E
c         |               .   .                
c         |                 a
c         +-------> X
c
c
c
c     Here, x1=a, x2=b, etc., and "a=x0" is the prescribed input.
c
c     Lower case letters indicate vertices on the z=0 plane.
c
c     Upper case on the upper Z-plane
c
c     Separation (e.g., E-D) is d.
c
      real x0(3)
      real xs(3),x(3,6),p(3,8)
c
      call copy(xs,x0,3)
      xs(1) = x0(1) - d
      call get_vert_lattice(p,d,xs)
c
      call copy(x(1,1),p(1,2),3)
      call copy(x(1,2),p(1,4),3)
      call copy(x(1,3),p(1,3),3)
      call copy(x(1,4),p(1,5),3)
      call copy(x(1,5),p(1,6),3)
      call copy(x(1,6),p(1,7),3)
c
      write(6,*)
      do k=1,6
         write(6,1) k,(x(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   dia')
c
      return
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_tet(x,d,x0)
c
c     Get vertices for a tet with lattice spacing d
c
c     x(:,1) is presumed known, and x(:,2) is on the "x-axis"
c
      real x0(3)
      real x(3,4)
c
      x(1,1) = 0                    +  x0(1)
      x(2,1) = 0                    +  x0(2)
      x(3,1) = 0                    +  x0(3)
c
      x(1,2) = d                    +  x0(1)
      x(2,2) = 0                    +  x0(2)
      x(3,2) = 0                    +  x0(3)
c
      x3     = 3
      x(1,3) = d*.5                 +  x0(1)
      x(2,3) = d*.5*sqrt(x3)        +  x0(2)
      x(3,3) = 0                    +  x0(3)
c
      one    = 1.
      pi6    = 4.*atan(one)/6.
      y3     = .5*tan(pi6)
      z3     = 1 - .25 - y3*y3
      z3     = sqrt(z3)
c
      x(1,4) = d*.5                 +  x0(1)
      x(2,4) = d*y3                 +  x0(2)
      x(3,4) = d*z3                 +  x0(3)
c
      write(6,*)
      do k=1,4
         write(6,1) k,(x(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   tets')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sdl_element_dia(v,xdia,rad)
c
c     Update SEM element data for diamond saddle
c
      include 'basics.inc'
      real v(3,8,24),xdia(3,6)
c
      integer e
      integer e2pfv(8)
      save    e2pfv
      data    e2pfv / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      k = 0
      do idia=1,6
      do itel=1,4
         k = k+1
         e = nel+k
         do i=1,8
            j      = e2pfv(i)
            x(e,i) = v(1,j,k)
            y(e,i) = v(2,j,k)
            z(e,i) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xdia(1,idia)
         curve (2,5,e) = xdia(2,idia)
         curve (3,5,e) = xdia(3,idia)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
c
         cbc   (5,e,1) = 'v  '
         cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature
c
         call rzero(bc(1,5,e,1),5)
         call rzero(bc(1,5,e,2),5)
c
         ilet      = mod1(e,52)
         letapt(e) = alphabet(ilet)
         numapt(e) = 1
c
      enddo
      enddo
c
      nel = nel+24
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sdl_element_tet(v,xtet,rad)
c
c     Update SEM element data for tet saddle
c
      include 'basics.inc'
      real v(3,8,12),xtet(3,4)
c
      integer e
      integer e2pfv(8)
      save    e2pfv
      data    e2pfv / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      k = 0
      do itet=1,4
      do itel=1,3
         k = k+1
         e = nel+k
         do i=1,8
            j      = e2pfv(i)
            x(e,i) = v(1,j,k)
            y(e,i) = v(2,j,k)
            z(e,i) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xtet(1,itet)
         curve (2,5,e) = xtet(2,itet)
         curve (3,5,e) = xtet(3,itet)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
c
         cbc   (5,e,1) = 'v  '
         cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature
c
         call rzero(bc(1,5,e,1),5)
         call rzero(bc(1,5,e,2),5)
c
         ilet      = mod1(e,52)
         letapt(e) = alphabet(ilet)
         numapt(e) = 1
c
      enddo
      enddo
c
      nel = nel+12
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_tet_saddle(v,x,d,r)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular tet with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,12),x(3,4)
      integer e
c
c
c
c     12 elements
c
      e = 0
      call tet_saddle_element(v(1,e+1),x(1,1),x(1,2),x(1,3),x(1,4),r)
      call tet_saddle_element(v(1,e+2),x(1,1),x(1,4),x(1,2),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,1),x(1,3),x(1,4),x(1,2),r)
c
      e = e+3  ! Need a-cyclic permutation here
      call tet_saddle_element(v(1,e+1),x(1,2),x(1,4),x(1,3),x(1,1),r)
      call tet_saddle_element(v(1,e+2),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,2),x(1,3),x(1,1),x(1,4),r)
c
      e = e+3
      call tet_saddle_element(v(1,e+1),x(1,3),x(1,4),x(1,1),x(1,2),r)
      call tet_saddle_element(v(1,e+2),x(1,3),x(1,2),x(1,4),x(1,1),r)
      call tet_saddle_element(v(1,e+3),x(1,3),x(1,1),x(1,2),x(1,4),r)
c
      e = e+3  ! Need a-cyclic permutation here
      call tet_saddle_element(v(1,e+1),x(1,4),x(1,3),x(1,2),x(1,1),r)
      call tet_saddle_element(v(1,e+2),x(1,4),x(1,1),x(1,3),x(1,2),r)
      call tet_saddle_element(v(1,e+3),x(1,4),x(1,2),x(1,1),x(1,3),r)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine tet_saddle_element(v,x1,x2,x3,x4,r)
c
c     Build a single element, situated on sph1, connecting to sph2
c
      real v(3,8),x1(3),x2(3),x3(3),x4(3)
c
      do i=1,3 
         v(i,5) = (x1(i)+x2(i)+x3(i)+x4(i))/4
         v(i,6) = (x1(i)+x2(i)+x4(i))/3
         v(i,7) = (x1(i)+x2(i)+x3(i))/3
         v(i,8) = (x1(i)+x2(i))/2
      enddo
c
      do k=1,4 
         call sub3 (v(1,k),v(1,k+4),x1,3)
         vnrm = vlsc2(v(1,k),v(1,k),3)
         vnrm = r/sqrt(vnrm)
         call cmult(v(1,k),vnrm,3)
         call add2 (v(1,k),x1,3)
      enddo
c
      write(6,*)
      do k=1,8
         write(6,1) k,(v(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   saddle')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dia_sdle(v,x1,x2,x3,x4,x5,x6,r)
c
c     Build a single element, situated on sph1, connecting to sph2
c
      real v(3,8),x1(3),x2(3),x3(3),x4(3),x5(3),x6(3)
c
      do i=1,3 
         v(i,5) = (x1(i)+x2(i)+x3(i)+x4(i)+x5(i)+x6(i))/6
         v(i,6) = (x1(i)+x2(i)+x5(i))/3
         v(i,7) = (x1(i)+x2(i)+x3(i))/3
         v(i,8) = (x1(i)+x2(i))/2
      enddo
c
      do k=1,4 
         call sub3 (v(1,k),v(1,k+4),x1,3)
         vnrm = vlsc2(v(1,k),v(1,k),3)
         vnrm = r/sqrt(vnrm)
         call cmult(v(1,k),vnrm,3)
         call add2 (v(1,k),x1,3)
      enddo
c
      write(6,*)
      do k=1,8
         write(6,1) k,(v(i,k),i=1,3),r
      enddo
      write(6,*)
    1 format(i3,1p4e12.4,' dia saddle')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_dia_saddle(v,x,d,r,x0)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular diamond with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,24),x(3,6),x0(3)
      integer e,face
c
c
c     24 elements, 4 for each of 6 nodes
c
      call dia_sdle(v(1, 1),x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),x(1,6),r)
      call dia_sdle(v(1, 2),x(1,1),x(1,5),x(1,2),x(1,3),x(1,4),x(1,6),r)
      call dia_sdle(v(1, 3),x(1,1),x(1,4),x(1,5),x(1,2),x(1,3),x(1,6),r)
      call dia_sdle(v(1, 4),x(1,1),x(1,3),x(1,4),x(1,5),x(1,2),x(1,6),r)
c
      call dia_sdle(v(1, 5),x(1,2),x(1,3),x(1,1),x(1,5),x(1,6),x(1,4),r)
      call dia_sdle(v(1, 6),x(1,2),x(1,6),x(1,3),x(1,1),x(1,5),x(1,4),r)
      call dia_sdle(v(1, 7),x(1,2),x(1,5),x(1,6),x(1,3),x(1,1),x(1,4),r)
      call dia_sdle(v(1, 8),x(1,2),x(1,1),x(1,5),x(1,6),x(1,3),x(1,4),r)
c
      call dia_sdle(v(1, 9),x(1,3),x(1,4),x(1,1),x(1,2),x(1,6),x(1,5),r)
      call dia_sdle(v(1,10),x(1,3),x(1,6),x(1,4),x(1,1),x(1,2),x(1,5),r)
      call dia_sdle(v(1,11),x(1,3),x(1,2),x(1,6),x(1,4),x(1,1),x(1,5),r)
      call dia_sdle(v(1,12),x(1,3),x(1,1),x(1,2),x(1,6),x(1,4),x(1,5),r)
c
      call dia_sdle(v(1,13),x(1,4),x(1,5),x(1,1),x(1,3),x(1,6),x(1,2),r)
      call dia_sdle(v(1,14),x(1,4),x(1,6),x(1,5),x(1,1),x(1,3),x(1,2),r)
      call dia_sdle(v(1,15),x(1,4),x(1,3),x(1,6),x(1,5),x(1,1),x(1,2),r)
      call dia_sdle(v(1,16),x(1,4),x(1,1),x(1,3),x(1,6),x(1,5),x(1,2),r)
c
      call dia_sdle(v(1,17),x(1,5),x(1,2),x(1,1),x(1,4),x(1,6),x(1,3),r)
      call dia_sdle(v(1,18),x(1,5),x(1,6),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call dia_sdle(v(1,19),x(1,5),x(1,4),x(1,6),x(1,2),x(1,1),x(1,3),r)
      call dia_sdle(v(1,20),x(1,5),x(1,1),x(1,4),x(1,6),x(1,2),x(1,3),r)
c
      call dia_sdle(v(1,21),x(1,6),x(1,2),x(1,5),x(1,4),x(1,3),x(1,1),r)
      call dia_sdle(v(1,22),x(1,6),x(1,3),x(1,2),x(1,5),x(1,4),x(1,1),r)
      call dia_sdle(v(1,23),x(1,6),x(1,4),x(1,3),x(1,2),x(1,5),x(1,1),r)
      call dia_sdle(v(1,24),x(1,6),x(1,5),x(1,4),x(1,3),x(1,2),x(1,1),r)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_tet(d,rad,x0)
c
c     Build the void between 4 spheres on a regular tet lattice
c
      real x0(3)
      real v(3,8,12),x(3,4)
c
      call get_vert_tet     (x,d,x0)
      call build_tet_saddle (v,x,d,rad)
      call sdl_element_tet  (v,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_dia(d,rad,x0)
c
c     Build the void between 6 spheres on a regular diamond lattice
c
      real x0(3)
      real v(3,8,24),x(3,6)
c
      call get_vert_dia     (x,d,x0)
      call build_dia_saddle (v,x,d,rad)
      call sdl_element_dia  (v,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_lat(d,rad,x0)
c
c     Build the void between 8 spheres on a regular diamond lattice
c
      include 'basics.inc'
      real x0(3)
      real v(3,8,48),p(3,8),q(3,8)  ! 48 elements total (12 + 24 + 12)
c
      call get_vert_lattice (p,d,x0)
c
      call copy(q(1,1),p(1,1),3)  ! 1st Tet
      call copy(q(1,2),p(1,2),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call build_tet_saddle(v(1,1, 1),q,d,rad)
      call sdl_element_tet (v(1,1, 1),q,rad)
c
      if_lattice = .true.                   ! Store unit lattice vectors
      call sub3   (ulat1,q(1,2),q(1,1),3)   ! -- needed for periodic bcs
      call norm3d (ulat1)
      call sub3   (ulat2,q(1,3),q(1,1),3)
      call norm3d (ulat2)
      call sub3   (ulat3,q(1,4),q(1,1),3)
      call norm3d (ulat3)
c
      call copy(q(1,1),p(1,2),3)  ! Diamond
      call copy(q(1,2),p(1,4),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call copy(q(1,5),p(1,6),3)
      call copy(q(1,6),p(1,7),3)
      call build_dia_saddle(v(1,1,13),q,d,rad)
      call sdl_element_dia (v(1,1,13),q,rad)
c
      call copy(q(1,1),p(1,4),3)  ! 1st Tet
      call copy(q(1,2),p(1,6),3)
      call copy(q(1,3),p(1,8),3)
      call copy(q(1,4),p(1,7),3)
      call build_tet_saddle(v(1,1,37),q,d,rad)
      call sdl_element_tet (v(1,1,37),q,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_lattice_per_bc
c
      include 'basics.inc'
c
      real x0(3),tlat(3,4)
c
      call prs ('Is this a lattice for hex-close-packed spheres?:$')
      call res(ans,1)
c
      if (ans.eq.'n'.or.ans.eq.'N') return
      call rzero(x0,3)
      dlat = 1.
      call get_vert_tet(t,d,x0) ! Base tet
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_lattice(x,d,x0)
c
c     Get vertices for 8-noded 3D rombus that forms a periodic
c     cell for hexagonally close packed spheres
c
c     x(:,1)=x0, and x(:,2)=x0+d
c
      real x(3,8),x0(3)
      real t(3,4)
c
      call get_vert_tet(t,d,x0) ! Base tet
c
      z1 = t(3,4)
c
      call copy(x,t,3*3)
      call copy(x(1,5),t(1,4),3)
c
      x(1,4) = x(1,3) + d
      x(2,4) = x(2,3)
      x(3,4) = x(3,3)
c
      x(1,6) = x(1,5) + d
      x(2,6) = x(2,5)
      x(3,6) = x(3,5)
c
      dy = x(2,3) - x(2,1)
c
      x(1,7) = x(1,2)
      x(2,7) = x(2,5) + dy
      x(3,7) = x(3,5)
c
      x(1,8) = x(1,7) + d
      x(2,8) = x(2,7)
      x(3,8) = x(3,7)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_dia(x,d,x0)
c
c     Get vertices for a diamond with lattice spacing d
c
c     Equilateral diamond is oriented as follows:
c
c
c
c         Y                 F
c                          / \                 
c         ^            c ./ . \. b
c         |             ./     \.              
c         |            D---------E
c         |               .   .                
c         |                 a
c         +-------> X
c
c
c
c     Here, x1=a, x2=b, etc., and "a=x0" is the prescribed input.
c
c     Lower case letters indicate vertices on the z=0 plane.
c
c     Upper case on the upper Z-plane
c
c     Separation (e.g., E-D) is d.
c
      real x0(3)
      real xs(3),x(3,6),p(3,8)
c
      call copy(xs,x0,3)
      xs(1) = x0(1) - d
      call get_vert_lattice(p,d,xs)
c
      call copy(x(1,1),p(1,2),3)
      call copy(x(1,2),p(1,4),3)
      call copy(x(1,3),p(1,3),3)
      call copy(x(1,4),p(1,5),3)
      call copy(x(1,5),p(1,6),3)
      call copy(x(1,6),p(1,7),3)
c
      write(6,*)
      do k=1,6
         write(6,1) k,(x(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   dia')
c
      return
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_tet(x,d,x0)
c
c     Get vertices for a tet with lattice spacing d
c
c     x(:,1) is presumed known, and x(:,2) is on the "x-axis"
c
      real x0(3)
      real x(3,4)
c
      x(1,1) = 0                    +  x0(1)
      x(2,1) = 0                    +  x0(2)
      x(3,1) = 0                    +  x0(3)
c
      x(1,2) = d                    +  x0(1)
      x(2,2) = 0                    +  x0(2)
      x(3,2) = 0                    +  x0(3)
c
      x3     = 3
      x(1,3) = d*.5                 +  x0(1)
      x(2,3) = d*.5*sqrt(x3)        +  x0(2)
      x(3,3) = 0                    +  x0(3)
c
      one    = 1.
      pi6    = 4.*atan(one)/6.
      y3     = .5*tan(pi6)
      z3     = 1 - .25 - y3*y3
      z3     = sqrt(z3)
c
      x(1,4) = d*.5                 +  x0(1)
      x(2,4) = d*y3                 +  x0(2)
      x(3,4) = d*z3                 +  x0(3)
c
      write(6,*)
      do k=1,4
         write(6,1) k,(x(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   tets')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sdl_element_dia(v,xdia,rad)
c
c     Update SEM element data for diamond saddle
c
      include 'basics.inc'
      real v(3,8,24),xdia(3,6)
c
      integer e
      integer e2pfv(8)
      save    e2pfv
      data    e2pfv / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      k = 0
      do idia=1,6
      do itel=1,4
         k = k+1
         e = nel+k
         do i=1,8
            j      = e2pfv(i)
            x(e,i) = v(1,j,k)
            y(e,i) = v(2,j,k)
            z(e,i) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xdia(1,idia)
         curve (2,5,e) = xdia(2,idia)
         curve (3,5,e) = xdia(3,idia)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
c
         cbc   (5,e,1) = 'v  '
         cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature
c
         call rzero(bc(1,5,e,1),5)
         call rzero(bc(1,5,e,2),5)
c
         ilet      = mod1(e,52)
         letapt(e) = alphabet(ilet)
         numapt(e) = 1
c
      enddo
      enddo
c
      nel = nel+24
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sdl_element_tet(v,xtet,rad)
c
c     Update SEM element data for tet saddle
c
      include 'basics.inc'
      real v(3,8,12),xtet(3,4)
c
      integer e
      integer e2pfv(8)
      save    e2pfv
      data    e2pfv / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      k = 0
      do itet=1,4
      do itel=1,3
         k = k+1
         e = nel+k
         do i=1,8
            j      = e2pfv(i)
            x(e,i) = v(1,j,k)
            y(e,i) = v(2,j,k)
            z(e,i) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xtet(1,itet)
         curve (2,5,e) = xtet(2,itet)
         curve (3,5,e) = xtet(3,itet)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
c
         cbc   (5,e,1) = 'v  '
         cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature
c
         call rzero(bc(1,5,e,1),5)
         call rzero(bc(1,5,e,2),5)
c
         ilet      = mod1(e,52)
         letapt(e) = alphabet(ilet)
         numapt(e) = 1
c
      enddo
      enddo
c
      nel = nel+12
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_tet_saddle(v,x,d,r)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular tet with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,12),x(3,4)
      integer e
c
c
c
c     12 elements
c
      e = 0
      call tet_saddle_element(v(1,e+1),x(1,1),x(1,2),x(1,3),x(1,4),r)
      call tet_saddle_element(v(1,e+2),x(1,1),x(1,4),x(1,2),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,1),x(1,3),x(1,4),x(1,2),r)
c
      e = e+3  ! Need a-cyclic permutation here
      call tet_saddle_element(v(1,e+1),x(1,2),x(1,4),x(1,3),x(1,1),r)
      call tet_saddle_element(v(1,e+2),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,2),x(1,3),x(1,1),x(1,4),r)
c
      e = e+3
      call tet_saddle_element(v(1,e+1),x(1,3),x(1,4),x(1,1),x(1,2),r)
      call tet_saddle_element(v(1,e+2),x(1,3),x(1,2),x(1,4),x(1,1),r)
      call tet_saddle_element(v(1,e+3),x(1,3),x(1,1),x(1,2),x(1,4),r)
c
      e = e+3  ! Need a-cyclic permutation here
      call tet_saddle_element(v(1,e+1),x(1,4),x(1,3),x(1,2),x(1,1),r)
      call tet_saddle_element(v(1,e+2),x(1,4),x(1,1),x(1,3),x(1,2),r)
      call tet_saddle_element(v(1,e+3),x(1,4),x(1,2),x(1,1),x(1,3),r)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine tet_saddle_element(v,x1,x2,x3,x4,r)
c
c     Build a single element, situated on sph1, connecting to sph2
c
      real v(3,8),x1(3),x2(3),x3(3),x4(3)
c
      do i=1,3 
         v(i,5) = (x1(i)+x2(i)+x3(i)+x4(i))/4
         v(i,6) = (x1(i)+x2(i)+x4(i))/3
         v(i,7) = (x1(i)+x2(i)+x3(i))/3
         v(i,8) = (x1(i)+x2(i))/2
      enddo
c
      do k=1,4 
         call sub3 (v(1,k),v(1,k+4),x1,3)
         vnrm = vlsc2(v(1,k),v(1,k),3)
         vnrm = r/sqrt(vnrm)
         call cmult(v(1,k),vnrm,3)
         call add2 (v(1,k),x1,3)
      enddo
c
      write(6,*)
      do k=1,8
         write(6,1) k,(v(i,k),i=1,3)
      enddo
      write(6,*)
    1 format(i3,1p3e12.4,'   saddle')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dia_sdle(v,x1,x2,x3,x4,x5,x6,r)
c
c     Build a single element, situated on sph1, connecting to sph2
c
      real v(3,8),x1(3),x2(3),x3(3),x4(3),x5(3),x6(3)
c
      do i=1,3 
         v(i,5) = (x1(i)+x2(i)+x3(i)+x4(i)+x5(i)+x6(i))/6
         v(i,6) = (x1(i)+x2(i)+x5(i))/3
         v(i,7) = (x1(i)+x2(i)+x3(i))/3
         v(i,8) = (x1(i)+x2(i))/2
      enddo
c
      do k=1,4 
         call sub3 (v(1,k),v(1,k+4),x1,3)
         vnrm = vlsc2(v(1,k),v(1,k),3)
         vnrm = r/sqrt(vnrm)
         call cmult(v(1,k),vnrm,3)
         call add2 (v(1,k),x1,3)
      enddo
c
      write(6,*)
      do k=1,8
         write(6,1) k,(v(i,k),i=1,3),r
      enddo
      write(6,*)
    1 format(i3,1p4e12.4,' dia saddle')
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_dia_saddle(v,x,d,r,x0)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular diamond with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,24),x(3,6),x0(3)
      integer e,face
c
c
c     24 elements, 4 for each of 6 nodes
c
      call dia_sdle(v(1, 1),x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),x(1,6),r)
      call dia_sdle(v(1, 2),x(1,1),x(1,5),x(1,2),x(1,3),x(1,4),x(1,6),r)
      call dia_sdle(v(1, 3),x(1,1),x(1,4),x(1,5),x(1,2),x(1,3),x(1,6),r)
      call dia_sdle(v(1, 4),x(1,1),x(1,3),x(1,4),x(1,5),x(1,2),x(1,6),r)
c
      call dia_sdle(v(1, 5),x(1,2),x(1,3),x(1,1),x(1,5),x(1,6),x(1,4),r)
      call dia_sdle(v(1, 6),x(1,2),x(1,6),x(1,3),x(1,1),x(1,5),x(1,4),r)
      call dia_sdle(v(1, 7),x(1,2),x(1,5),x(1,6),x(1,3),x(1,1),x(1,4),r)
      call dia_sdle(v(1, 8),x(1,2),x(1,1),x(1,5),x(1,6),x(1,3),x(1,4),r)
c
      call dia_sdle(v(1, 9),x(1,3),x(1,4),x(1,1),x(1,2),x(1,6),x(1,5),r)
      call dia_sdle(v(1,10),x(1,3),x(1,6),x(1,4),x(1,1),x(1,2),x(1,5),r)
      call dia_sdle(v(1,11),x(1,3),x(1,2),x(1,6),x(1,4),x(1,1),x(1,5),r)
      call dia_sdle(v(1,12),x(1,3),x(1,1),x(1,2),x(1,6),x(1,4),x(1,5),r)
c
      call dia_sdle(v(1,13),x(1,4),x(1,5),x(1,1),x(1,3),x(1,6),x(1,2),r)
      call dia_sdle(v(1,14),x(1,4),x(1,6),x(1,5),x(1,1),x(1,3),x(1,2),r)
      call dia_sdle(v(1,15),x(1,4),x(1,3),x(1,6),x(1,5),x(1,1),x(1,2),r)
      call dia_sdle(v(1,16),x(1,4),x(1,1),x(1,3),x(1,6),x(1,5),x(1,2),r)
c
      call dia_sdle(v(1,17),x(1,5),x(1,2),x(1,1),x(1,4),x(1,6),x(1,3),r)
      call dia_sdle(v(1,18),x(1,5),x(1,6),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call dia_sdle(v(1,19),x(1,5),x(1,4),x(1,6),x(1,2),x(1,1),x(1,3),r)
      call dia_sdle(v(1,20),x(1,5),x(1,1),x(1,4),x(1,6),x(1,2),x(1,3),r)
c
      call dia_sdle(v(1,21),x(1,6),x(1,2),x(1,5),x(1,4),x(1,3),x(1,1),r)
      call dia_sdle(v(1,22),x(1,6),x(1,3),x(1,2),x(1,5),x(1,4),x(1,1),r)
      call dia_sdle(v(1,23),x(1,6),x(1,4),x(1,3),x(1,2),x(1,5),x(1,1),r)
      call dia_sdle(v(1,24),x(1,6),x(1,5),x(1,4),x(1,3),x(1,2),x(1,1),r)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_tet(d,rad,x0)
c
c     Build the void between 4 spheres on a regular tet lattice
c
      real x0(3)
      real v(3,8,12),x(3,4)
c
      call get_vert_tet     (x,d,x0)
      call build_tet_saddle (v,x,d,rad)
      call sdl_element_tet  (v,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_dia(d,rad,x0)
c
c     Build the void between 6 spheres on a regular diamond lattice
c
      real x0(3)
      real v(3,8,24),x(3,6)
c
      call get_vert_dia     (x,d,x0)
      call build_dia_saddle (v,x,d,rad)
      call sdl_element_dia  (v,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_lat(d,rad,x0)
c
c     Build the void between 8 spheres on a regular diamond lattice
c
      include 'basics.inc'
      real x0(3)
      real v(3,8,48),p(3,8),q(3,8)  ! 48 elements total (12 + 24 + 12)
c
      call get_vert_lattice (p,d,x0)
c
      call copy(q(1,1),p(1,1),3)  ! 1st Tet
      call copy(q(1,2),p(1,2),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call build_tet_saddle(v(1,1, 1),q,d,rad)
      call sdl_element_tet (v(1,1, 1),q,rad)
c
      if_lattice = .true.                   ! Store unit lattice vectors
      call sub3   (ulat1,q(1,2),q(1,1),3)   ! -- needed for periodic bcs
      call norm3d (ulat1)
      call sub3   (ulat2,q(1,3),q(1,1),3)
      call norm3d (ulat2)
      call sub3   (ulat3,q(1,4),q(1,1),3)
      call norm3d (ulat3)
c
      call copy(q(1,1),p(1,2),3)  ! Diamond
      call copy(q(1,2),p(1,4),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call copy(q(1,5),p(1,6),3)
      call copy(q(1,6),p(1,7),3)
      call build_dia_saddle(v(1,1,13),q,d,rad)
      call sdl_element_dia (v(1,1,13),q,rad)
c
      call copy(q(1,1),p(1,4),3)  ! 1st Tet
      call copy(q(1,2),p(1,6),3)
      call copy(q(1,3),p(1,8),3)
      call copy(q(1,4),p(1,7),3)
      call build_tet_saddle(v(1,1,37),q,d,rad)
      call sdl_element_tet (v(1,1,37),q,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_lattice_per_bc
c
c     Build the void between 8 spheres on a regular diamond lattice
c
      include 'basics.inc'
      real x0(3)
      real v(3,8,48),p(3,8),q(3,8)  ! 48 elements total (12 + 24 + 12)
c
      call get_vert_lattice (p,d,x0)
c
      call copy(q(1,1),p(1,1),3)  ! 1st Tet
      call copy(q(1,2),p(1,2),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call build_tet_saddle(v(1,1, 1),q,d,rad)
      call sdl_element_tet (v(1,1, 1),q,rad)
c
      if_lattice = .true.                   ! Store unit lattice vectors
      call sub3   (ulat1,q(1,2),q(1,1),3)   ! -- needed for periodic bcs
      call norm3d (ulat1)
      call sub3   (ulat2,q(1,3),q(1,1),3)
      call norm3d (ulat2)
      call sub3   (ulat3,q(1,4),q(1,1),3)
      call norm3d (ulat3)
c
      call copy(q(1,1),p(1,2),3)  ! Diamond
      call copy(q(1,2),p(1,4),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call copy(q(1,5),p(1,6),3)
      call copy(q(1,6),p(1,7),3)
      call build_dia_saddle(v(1,1,13),q,d,rad)
      call sdl_element_dia (v(1,1,13),q,rad)
c
      call copy(q(1,1),p(1,4),3)  ! 1st Tet
      call copy(q(1,2),p(1,6),3)
      call copy(q(1,3),p(1,8),3)
      call copy(q(1,4),p(1,7),3)
      call build_tet_saddle(v(1,1,37),q,d,rad)
      call sdl_element_tet (v(1,1,37),q,rad)
c
      return
      end
c-----------------------------------------------------------------------
c
c     Build the void between 8 spheres on a regular diamond lattice
c
      include 'basics.inc'
      real x0(3)
      real v(3,8,48),p(3,8),q(3,8)  ! 48 elements total (12 + 24 + 12)
c
      call get_vert_lattice (p,d,x0)
c
      call copy(q(1,1),p(1,1),3)  ! 1st Tet
      call copy(q(1,2),p(1,2),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call build_tet_saddle(v(1,1, 1),q,d,rad)
      call sdl_element_tet (v(1,1, 1),q,rad)
c
