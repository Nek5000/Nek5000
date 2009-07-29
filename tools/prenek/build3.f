c-----------------------------------------------------------------------
      SUBROUTINE CYLMESH
C
C     Build cylindrical meshes, in fashion similar to spherical meshes
C
c     7/25/04 pff
c
      INCLUDE 'basics.inc'
      COMMON /CTMP0/ SPHCTR(3),XCS(4,24),YCS(4,24),ZCS(4,24)
      CHARACTER*1 SHELL,HEMI,YESNO
      CHARACTER*1 alphabet(52)
      CHARACTER*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
C
c     CALL PRS('A 6 or 24 element mesh?$')
c     CALL REI(MESH)
      MESH = 24
C
      CALL PRS('Enter the (X,Y,Z) coordinates of the center:$')
      CALL RERRR(SPHCTR(1),SPHCTR(2),SPHCTR(3))
C
      CALL PRS('Hemisphere, whole sphere, or hex-pkd hemi? (H/W/X):$')
      CALL RES(HEMI,1)
      CALL CAPIT(HEMI,1)
      IF (HEMI.EQ.'H') THEN
         NELSPH=12
         IF (MESH.EQ.6) NELSPH=5
      ELSEIF (HEMI.EQ.'X') THEN
         NELSPH=12
         IF (MESH.EQ.6) NELSPH=5
      ELSE
         CALL PRS('A 6 or 24 element mesh?$')
         CALL REI(MESH)
         NELSPH=MESH
      ENDIF
      NLSPH4=4*NELSPH
C
C---------------------------------------------------------------------------
C     Begin building shell sequence, working from the inner most to outer.
C---------------------------------------------------------------------------
C
      DO 6000 ISHLL=1,1000
C
         CALL PRS(
     $  'Enter C or B for a cylinder or box layer (E=exit):$')
         CALL RES(SHELL,1)
         CALL CAPIT(SHELL,1)
         IF (SHELL.EQ.'E') GOTO 9000
C
         IF (SHELL.EQ.'C') THEN
            CALL PRS('Enter radius:$')
            CALL RER(RADIUS)
            CALL CYLSHELL(XCS,YCS,ZCS,RADIUS)
            CALL TRANS2(XCS,YCS,ZCS,CYLCTR,NLSPH4)
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
               X(IEL,I)=XCS(I,IE)
               Y(IEL,I)=YCS(I,IE)
               Z(IEL,I)=ZCS(I,IE)
  101       CONTINUE

            IF (SHELL.EQ.'C') THEN
               CCURVE(4,IEL)='C'      ! Circle
               CURVE(1,4,IEL)=RADIUS
               CURVE(2,4,IEL)=0
               CURVE(3,4,IEL)=0
               CURVE(4,4,IEL)=0
               CURVE(5,4,IEL)=0
               if (if3d) then
                  CCURVE(8,IEL)='C'      ! Circle
                  CURVE(1,8,IEL)=RADIUS
                  CURVE(2,8,IEL)=0
                  CURVE(3,8,IEL)=0
                  CURVE(4,8,IEL)=0
                  CURVE(5,8,IEL)=0
               endif
            ENDIF
            IF (ISHLL.GT.1) THEN
               IEL=NEL-NELSPH+IE
               DO 201 I=1,4
                  J=I+4
                  X(IEL,J)=XCS(I,IE)
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
      RETURN
      END
      SUBROUTINE SPHERE(XCS,YCS,ZCS,HEMI,MESH,RADIUS)
      DIMENSION XCS(4,24),YCS(4,24),ZCS(4,24)
      CHARACTER*1 HEMI
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
      RETURN
      END
      SUBROUTINE CRTBOX(XCS,YCS,ZCS,HEMI,MESH,RADIUS)
      DIMENSION XCS(4,24),YCS(4,24),ZCS(4,24)
      CHARACTER*1 HEMI
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
      RETURN
      END
      SUBROUTINE ROTAT2(X,Y,Z,N,DIR,ANGLE)
C
      DIMENSION X(1),Y(1),Z(1)
      CHARACTER*1 DIR
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
      RETURN
      END
      SUBROUTINE ROUNDER(X,N)
      DIMENSION X(1)
C
C     Try to Round X to fractional integer - eg .05 .1 .15 - if it's within 10-6
C
      DO 100 I=1,N
         EPS=1.0E-4
         XTMP=20.0*X(I)
         XTMP2=XTMP+0.5
         ITMP=INT(XTMP2)
         XTMP2=FLOAT(ITMP)
         IF (ABS(XTMP-XTMP2).LT.EPS) X(I)=.05*itmp
  100 CONTINUE
      RETURN
      END
      SUBROUTINE TRANS2(X,Y,Z,XYZOFF,N)
      DIMENSION X(1),Y(1),Z(1)
      DIMENSION XYZOFF(3)
C
      DO 10 I=1,N
         X(I)=X(I)+XYZOFF(1)
         Y(I)=Y(I)+XYZOFF(2)
         Z(I)=Z(I)+XYZOFF(3)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE GENMESH
      INCLUDE 'basics.inc'
C
      COMMON /CTMP0/ SPHCTR(3),XCS(4,24),YCS(4,24),ZCS(4,24)
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
      RETURN
      END
      SUBROUTINE GENELE(XYZ)
      DIMENSION XYZ(3,8)
      INCLUDE 'basics.inc'
C
      NEL=NEL+1
      DO 10 I=1,8
         X(NEL,I)=XYZ(1,I)
         Y(NEL,I)=XYZ(2,I)
         Z(NEL,I)=XYZ(3,I)
   10 CONTINUE
      NUMAPT(NEL)=ILEVEL
      LETAPT(NEL)='A'
      RETURN
      END
      SUBROUTINE ROTATE(XYZ,N,DIR,ANGLE)
C
      DIMENSION XYZ(3,1)
      CHARACTER*1 DIR
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
      RETURN
      END
      SUBROUTINE TRANSL(XYZ,XYZOFF,N)
      DIMENSION XYZ(3,1)
      DIMENSION XYZOFF(3)
C
      DO 10 I=1,N
         XYZ(1,I)=XYZ(1,I)+XYZOFF(1)
         XYZ(2,I)=XYZ(2,I)+XYZOFF(2)
         XYZ(3,I)=XYZ(3,I)+XYZOFF(3)
   10 CONTINUE
      RETURN
      END
