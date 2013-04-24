c-----------------------------------------------------------------------
cc   TODO:
c
c .check (nelm-1) and nelm usage to ensure safe_haven is ok
c
c .need arc_circle transform with:
c
c    .arclength preservation (a bit of calculus)
c
c    .need midside node adjust so that midside nodes stay centered
c       -- would the "deviation" (d) ccurve idea work here?
c       -- It should, if we compute and preserve the relative
c          orthonormal vectors so that the deflection is still
c          orthogonal --- _and_ the magnitude change of the normal
c          vector coincides with that of the chord; but only if
c          the magnitude change is computed according to the curvature
c          rather than the length.  Curvature is easy in this case
c          because 'm' implies a parabola!  So, a little algebra
c          will save this.
c
c          What about the sphere ?  That should be OK.  Ditto for 'C'
c
c          .So, idea for 'm' is the following:
c
c               - Pass v1,m,v2 and v1' v2' into a routine that returns m'.
c
c               - Easy!
c
c-----------------------------------------------------------------------
      subroutine mid_3d(xm,x0,x1)  ! xm = (x0+x1)/2

      real xm(3),x0(3),x1(3)

      call add3 (xm,x0,x1,3)
      call cmult(xm,0.5,3)

      return
      end
c-----------------------------------------------------------------------
      subroutine midside_convert_all ! Convert all elements to midside
      include 'basics.inc'
      integer e
      do e=1,nel
         call fix_m_curve(e) ! Convert all elements to midside
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine midside_clean(xm,x0,x1,linear) 

c     clean up midside node so it's projection is midpoint of (x0,x1)

      real xm(3),x0(3),x1(3)

      real v(3),d(3),n(3),xh(3),length,lp

      call sub3(v,x1,x0,3)            ! v  = x1-x0
      call normalize(v,length,3)      ! d := d/||d||, disp = ||d||
      call mid_3d(xh,x0,x1)           ! xh = (x0+x1)/2
      call sub3(d,xm,xh,3)            ! d  = xm-xh
      call normalize(d,disp,3)        ! d := d/||d||, disp = ||d||
      call vcross_normal(o,sine,v,d)  ! o  = v x d

      eps  = 2.e-5
      test = min(disp,sine)
      linear = 0
      if (test.lt.eps*length) then    ! point is nearly colinear with v
         linear = 1
         call mid_3d(xm,x0,x1)        ! set xm = (x0+x1)/2 and return
         return
      endif

c     p = projection of xm onto v

      call sub3   (d0,xm,x0,3)        ! d0 = xm-x0
      lp = dot    (d0,v,3)
      call add3s2 (p,x0,v,lp,3)       ! p = x0 + lp*v

      call sub3     (d,xm,p,3)        ! d = xm-p
      call normalize(d,disp,3)        ! d := d/||d||, disp = ||d||

      amp = disp/(lp*(length-lp))     ! a = y/(x*(l-x))
      call add3s2 (xm,xh,d,amp,3)     ! xm = xh + d*amp

      return
      end
c-----------------------------------------------------------------------
      subroutine midside_refine(xmp,x0p,x1p,xmi,x0,x1,linear)
      real xmp(3),x0p(3),x1p(3),xmi(3),x0(3),x1(3)

c     Compute new midpoint, xmp, assuming that (x0,x1) 
c     has been moved to (x0p,x1p).
c
c     Assumption is that displacement will be in same plane and
c     that magnitude of Curvature is preserved as v is shrunk.
c
c     This strategy is good for refinement but not for general
c     affine transformations.   pff 3/30/13

      real xm(3),o(3)
      real xh (3),v (3),n (3),l
      real xhp(3),vp(3),np(3),lp


      call copy          (xm,xmi,3)     ! Ensure quality input data
      call midside_clean (xm,x0,x1,linear) 
      if (linear.eq.1) then
         call mid_3d (xmp,x0p,x1p)      ! Midpoint of original vector
         return
      endif

      call mid_3d        (xh,x0,x1)     ! Midpoint of original vector
      call sub3          (n,xm,xh,3)    ! Deviation of original vector: n=xm-xh
      call normalize     (n,alpha,3)    ! alpha := | n |

      call sub3          (v ,x1,x0,3)   ! v = x1-x0
      call normalize     (v ,l ,3)      ! l = original length
      call sub3          (v ,x1p,x0p,3) ! v = x1-x0
      call normalize     (vp,lp,3)      ! l = original length

      call vcross_normal (o,sine,v,d)   ! bi-normal

      a = 4*alpha/(l*l)                 ! curvature

      alphap = a*(lp*lp)/4              ! amplitude of new curve
      call vcross_normal (np,sine,o,vp) ! normal for new midpoint

      call mid_3d        (xhp,x0p,x1p)  ! midpoint of original vector

      call add3s2(xmp,xhp,np,alphap,3)  ! xm' = xh' + n'*amp

      return
      end
c-----------------------------------------------------------------------
      subroutine sphmesh
      include 'basics.inc'
      common /ctmp0/ sphctr(3),xcs(4,24),ycs(4,24),zcs(4,24)
      common /ctmpr/ radii(100)
      character*1 SHELL,HEMI,YESNO
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      logical ifpsph,if_sph_ctr
c
      real xlat0(3)
C
C     Build either spherical or hemi-spherical meshes.
C
C     Modified to include stretching in X for prolate spheroird. pff 7/30/92
C
c     CALL PRS('A 6 or 24 element mesh?$')
c     CALL REI(MESH)
c
c        call get_lattice_0  (dlat,sphrad,xlat0)
c        call saddle_fcc     (dlat,sphrad,xlat0)
c        return
c
      MESH = 24

      call prs('H/W sph/hex-pkd hemi/tet/dia/lat/fcc/r12/cap:$')
      call prs('(H/W/X/T/D/L/F/R/C)?$')
      call res(hemi,1)
      call capit(hemi,1)

      if (hemi.eq.'C'.or.hemi.eq.'c') then

         radii(1)=1.0
         radii(2)=1.5
         radii(3)=2.0
         radii(4)=0.5  ! thickness of pipe
         radii(5)=0.25 ! 4x reduction factor
         nr       = 3
         call sc_make_sphere_cap(radii,nr)
         return

      elseif (hemi.eq.'2') then
         radii(1)=0.8
         radii(2)=1.2
         radii(3)=1.4
         radii(4)=1.5
         radii(5)=1.6
         radii(6)=1.65
         radii(7)=2.0
         radii(8)=0.5  ! thickness of pipe
         radii(9)=0.25 ! 4x reduction factor
         nr       = 7

         call sc_make_sphere_cap(radii,nr)
         return

      elseif (hemi.eq.'t'.or.hemi.eq.'T') then
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_tet     (dlat,sphrad,xlat0,if_lat_sph_cent)
         return

      elseif (hemi.eq.'d'.or.hemi.eq.'D') then
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_dia     (dlat,sphrad,xlat0,if_lat_sph_cent)
         return
c
      elseif (hemi.eq.'l'.or.hemi.eq.'L') then
         call get_lattice_0  (dlat,sphrad,xlat0)
         call saddle_lat     (dlat,sphrad,xlat0)
         return

      elseif (hemi.eq.'f'.or.hemi.eq.'F') then

c        call get_lattice_0  (dlat,sphrad,xlat0)
c        call saddle_fcc     (dlat,sphrad,xlat0)

         call fcc2  ! Build fcc mesh
         return

      elseif (hemi.eq.'r'.or.hemi.eq.'R') then
         dlat = 1.
         call rzero(xlat0,3)
         call rhombic_dodec  (dlat,xlat0,if_sph_ctr)
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
               x(i,iel)=xcs(i,ie)*ratio
               y(i,iel)=ycs(i,ie)
               z(i,iel)=zcs(i,ie)
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
                  j=i+4
                  x(j,iel)=xcs(i,ie)*ratio
                  y(j,iel)=ycs(i,ie)
                  z(j,iel)=zcs(i,ie)
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
      DO 9001 IEDGE=1,12
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
      nel=nel+1
      do 10 i=1,8
         x(i,nel)=xyz(1,i)
         y(i,nel)=xyz(2,i)
         z(i,nel)=xyz(3,i)
   10 continue
      numapt(nel)=ilevel
      letapt(nel)='A'
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
      call prs ('Input coordinates for 1st sphere or lattice ctr:$')
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
      subroutine permute_dia(x,p)
c
c     permute vertices for diamond
c
c     Input / Output diamonds oriented as follows:
c
c
c         Y                 6
c                          / \                 
c         ^            3 ./ . \. 1               INPUT
c         |             ./     \.              
c         |            2---------4
c         |               .   .                
c         |                 5
c         +-------> X
c
c
c         Y                 F
c                          / \                 
c         ^            c ./ . \. b               OUTPUT
c         |             ./     \.              
c         |            D---------E
c         |               .   .                
c         |                 a
c         +-------> X
c
c
c
c
      real x(3,6),p(3,6)
      character*1 blah
c
      call copy(x(1,2),p(1,1),3)
      call copy(x(1,4),p(1,2),3)
      call copy(x(1,3),p(1,3),3)
      call copy(x(1,5),p(1,4),3)
      call copy(x(1,1),p(1,5),3)
      call copy(x(1,6),p(1,6),3)
c
      call outmat(p,3,6,'before',1)
      call outmat(x,3,6,'permut',1)
      call prs('LOOK at window$')
      call res(blah,1)
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
      end
c-----------------------------------------------------------------------
      subroutine get_vert_fcc_lat(p,d,c0)
c
c     Get vertices for an fcc lattice with spacing d, ctr. c0:
c
c
      real p(3,8,2),c0(3)
      real x0(3)
c
      d2 = d/2
      do i=1,3
         x0(i) = c0(i) - d2
      enddo
c
      l = 0
      do k=0,1
      do j=0,1
      do i=0,1
         l = l+1
         p(1,l,1) = x0(1) + d*i
         p(2,l,1) = x0(2) + d*j
         p(3,l,1) = x0(3) + d*k
      enddo
      enddo
      enddo
c
      do k=1,6
         call copy(p(1,k,2),c0,3)
      enddo
c
      l = 0
      do k=1,3
      do i=-1,1,2
         l = l+1
         p(k,l,2) = c0(k) + d2*i
      enddo
      enddo
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
            x(i,e) = v(1,j,k)
            y(i,e) = v(2,j,k)
            z(i,e) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xdia(1,idia)
         curve (2,5,e) = xdia(2,idia)
         curve (3,5,e) = xdia(3,idia)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
c
         if (if_lat_sph_cent) call rzero(curve(1,5,e),3)
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
            x(i,e) = v(1,j,k)
            y(i,e) = v(2,j,k)
            z(i,e) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = 's'
         curve (1,5,e) = xtet(1,itet)
         curve (2,5,e) = xtet(2,itet)
         curve (3,5,e) = xtet(3,itet)
         curve (4,5,e) = rad
         curve (5,5,e) = 0.
         if (if_lat_sph_cent) call rzero(curve(1,5,e),3)
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
      subroutine build_tet_saddle(v,xi,d,r,if_sph_ctr)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular tet with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,12),xi(3,4)
      logical if_sph_ctr
c
      real x(3,4)
      integer e
      character*80 s
c
      call copy(x,xi,12)
      vol = tet_chk(xi)
      if (vol.lt.0) then
         write(6,*) 'in build.tet.saddle',vol
         call copy(x(1,1),xi(1,2),3)
         call copy(x(1,2),xi(1,1),3)
      endif
c
      c1 = ( x(1,1)+x(1,2)+x(1,3)+x(1,4) )/4
      c2 = ( x(2,1)+x(2,2)+x(2,3)+x(2,4) )/4
      c3 = ( x(3,1)+x(3,2)+x(3,3)+x(3,4) )/4
c
      write(s,1) c1,c2,c3,d,r
    1 format('Tet Ctr:',3g12.4,', Latt., R:',2g12.4,'$')
      call prs(s)
c
c     12 elements
c
      e = 0
      if (if_sph_ctr) call shift_xyz(x,x(1,1),4)
      call tet_saddle_element(v(1,e+1),x(1,1),x(1,2),x(1,3),x(1,4),r)
      call tet_saddle_element(v(1,e+2),x(1,1),x(1,4),x(1,2),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,1),x(1,3),x(1,4),x(1,2),r)
c
      e = e+3  ! Need a-cyclic permutation here
      if (if_sph_ctr) call shift_xyz(x,x(1,2),4)
      call tet_saddle_element(v(1,e+1),x(1,2),x(1,4),x(1,3),x(1,1),r)
      call tet_saddle_element(v(1,e+2),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call tet_saddle_element(v(1,e+3),x(1,2),x(1,3),x(1,1),x(1,4),r)
c
      e = e+3
      if (if_sph_ctr) call shift_xyz(x,x(1,3),4)
      call tet_saddle_element(v(1,e+1),x(1,3),x(1,4),x(1,1),x(1,2),r)
      call tet_saddle_element(v(1,e+2),x(1,3),x(1,2),x(1,4),x(1,1),r)
      call tet_saddle_element(v(1,e+3),x(1,3),x(1,1),x(1,2),x(1,4),r)
c
      e = e+3  ! Need a-cyclic permutation here
      if (if_sph_ctr) call shift_xyz(x,x(1,4),4)
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
      real v(3,8),x1(3),x2(3),x3(3),x4(3),w(3,8)
      logical ifnonsym
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
      ifnonsym = .true.
      if (ifnonsym) then  ! correction for non-regular tet
         call copy(w(1,6),v(1,6),3)
         call copy(w(1,7),v(1,7),3)
         call non_reg_tet_mod(v(1,7),x1,x2,x3,r)
         call non_reg_tet_mod(v(1,6),x1,x2,x4,r)
c
         do k=2,3 
            h  = .5
            call add2 (w(1,k+4),v(1,k+4),3)
            call cmult(w(1,k+4),h,3)
            call sub3 (v(1,k),w(1,k+4),x1,3)
            vnrm = vlsc2(v(1,k),v(1,k),3)
            vnrm = r/sqrt(vnrm)
            call cmult(v(1,k),vnrm,3)
            call add2 (v(1,k),x1,3)
         enddo
      endif
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
      subroutine build_dia_saddle(v,xi,d,r,if_sph_ctr)
c
c     Build a saddle to fill the void between 4 spheres of radius r
c     situated at the vertices of a regular diamond with distance d.
c
c     For now, assume that sphere 1 and 2 are separated along the x-axis.
c
      real v(3*8,24),xi(3,6)
      logical if_sph_ctr
c
      real x(3,6)
      integer e,face
      character*80 s
c
      call copy(x,xi,18)
c
      c1 = ( x(1,1)+x(1,2)+x(1,3)+x(1,4)+x(1,5)+x(1,6) )/6
      c2 = ( x(2,1)+x(2,2)+x(2,3)+x(2,4)+x(2,5)+x(2,6) )/6
      c3 = ( x(3,1)+x(3,2)+x(3,3)+x(3,4)+x(3,5)+x(3,6) )/6
c
      write(s,1) c1,c2,c3,d,r
    1 format('DIA Ctr:',3g12.4,', Latt., R:',2g12.4,'$')
      call prs(s)
c
c
c     24 elements, 4 for each of 6 nodes
c
      if (if_sph_ctr) call shift_xyz(x,x(1,1),6)
      call dia_sdle(v(1, 1),x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),x(1,6),r)
      call dia_sdle(v(1, 2),x(1,1),x(1,5),x(1,2),x(1,3),x(1,4),x(1,6),r)
      call dia_sdle(v(1, 3),x(1,1),x(1,4),x(1,5),x(1,2),x(1,3),x(1,6),r)
      call dia_sdle(v(1, 4),x(1,1),x(1,3),x(1,4),x(1,5),x(1,2),x(1,6),r)
c
      if (if_sph_ctr) call shift_xyz(x,x(1,2),6)
      call dia_sdle(v(1, 5),x(1,2),x(1,3),x(1,1),x(1,5),x(1,6),x(1,4),r)
      call dia_sdle(v(1, 6),x(1,2),x(1,6),x(1,3),x(1,1),x(1,5),x(1,4),r)
      call dia_sdle(v(1, 7),x(1,2),x(1,5),x(1,6),x(1,3),x(1,1),x(1,4),r)
      call dia_sdle(v(1, 8),x(1,2),x(1,1),x(1,5),x(1,6),x(1,3),x(1,4),r)
c
      if (if_sph_ctr) call shift_xyz(x,x(1,3),6)
      call dia_sdle(v(1, 9),x(1,3),x(1,4),x(1,1),x(1,2),x(1,6),x(1,5),r)
      call dia_sdle(v(1,10),x(1,3),x(1,6),x(1,4),x(1,1),x(1,2),x(1,5),r)
      call dia_sdle(v(1,11),x(1,3),x(1,2),x(1,6),x(1,4),x(1,1),x(1,5),r)
      call dia_sdle(v(1,12),x(1,3),x(1,1),x(1,2),x(1,6),x(1,4),x(1,5),r)
c
      if (if_sph_ctr) call shift_xyz(x,x(1,4),6)
      call dia_sdle(v(1,13),x(1,4),x(1,5),x(1,1),x(1,3),x(1,6),x(1,2),r)
      call dia_sdle(v(1,14),x(1,4),x(1,6),x(1,5),x(1,1),x(1,3),x(1,2),r)
      call dia_sdle(v(1,15),x(1,4),x(1,3),x(1,6),x(1,5),x(1,1),x(1,2),r)
      call dia_sdle(v(1,16),x(1,4),x(1,1),x(1,3),x(1,6),x(1,5),x(1,2),r)
c
      if (if_sph_ctr) call shift_xyz(x,x(1,5),6)
      call dia_sdle(v(1,17),x(1,5),x(1,2),x(1,1),x(1,4),x(1,6),x(1,3),r)
      call dia_sdle(v(1,18),x(1,5),x(1,6),x(1,2),x(1,1),x(1,4),x(1,3),r)
      call dia_sdle(v(1,19),x(1,5),x(1,4),x(1,6),x(1,2),x(1,1),x(1,3),r)
      call dia_sdle(v(1,20),x(1,5),x(1,1),x(1,4),x(1,6),x(1,2),x(1,3),r)
c
      if (if_sph_ctr) call shift_xyz(x,x(1,6),6)
      call dia_sdle(v(1,21),x(1,6),x(1,2),x(1,5),x(1,4),x(1,3),x(1,1),r)
      call dia_sdle(v(1,22),x(1,6),x(1,3),x(1,2),x(1,5),x(1,4),x(1,1),r)
      call dia_sdle(v(1,23),x(1,6),x(1,4),x(1,3),x(1,2),x(1,5),x(1,1),r)
      call dia_sdle(v(1,24),x(1,6),x(1,5),x(1,4),x(1,3),x(1,2),x(1,1),r)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_tet(d,rad,x0,if_sph_ctr)
c
c     Build the void between 4 spheres on a regular tet lattice
c
      real x0(3)
      real v(3,8,12),x(3,4)
      logical if_sph_ctr
c
      call get_vert_tet     (x,d,x0)
      call build_tet_saddle (v,x,d,rad,if_sph_ctr)
      call sdl_element_tet  (v,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_dia(d,rad,x0,if_sph_ctr)
c
c     Build the void between 6 spheres on a regular diamond lattice
c
      real x0(3)
      real v(3,8,24),x(3,6)
      logical if_sph_ctr
c
      call get_vert_dia     (x,d,x0)
      call build_dia_saddle (v,x,d,rad,if_sph_ctr)
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
      real a0(3)
c
      call rzero(a0,3)  ! put first sphere at (0,0,0)
      call get_vert_lattice (p,d,a0)
c
      do i=1,8                  ! center point of lattice
         a0(1) = a0(1) + p(1,i)/8
         a0(2) = a0(2) + p(2,i)/8
         a0(3) = a0(3) + p(3,i)/8
      enddo
      call sub2(a0,x0,3)
      call shift_xyz(p,a0,8)    ! p = p - a0
c
      call prs ('Enter s/v for sphere- or void-centric mesh:$')
      call res(ans,1)
      if_lat_sph_cent = .true.  ! Build sphere-centric mesh
      if (ans.eq.'v'.or.ans.eq.'V') if_lat_sph_cent = .false.
c
      call copy(q(1,1),p(1,1),3)  ! 1st Tet
      call copy(q(1,2),p(1,2),3)
      call copy(q(1,3),p(1,3),3)
      call copy(q(1,4),p(1,5),3)
      call build_tet_saddle(v(1,1, 1),q,d,rad,if_lat_sph_cent)
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
      call build_dia_saddle(v(1,1,13),q,d,rad,if_lat_sph_cent)
      call sdl_element_dia (v(1,1,13),q,rad)
c
      call copy(q(1,1),p(1,4),3)  ! 1st Tet
      call copy(q(1,2),p(1,6),3)
      call copy(q(1,3),p(1,8),3)
      call copy(q(1,4),p(1,7),3)
      call build_tet_saddle(v(1,1,37),q,d,rad,if_lat_sph_cent)
      call sdl_element_tet (v(1,1,37),q,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_lattice_per_bc
c
c     This routine will define the unit-lattice direction vectors
c     (ulat_k, k=1,...,3) that are used to determine coincidence
c     of periodic boundaries.
c
c
      include 'basics.inc'
c
      real x0(3),q(3,4)
c
      call prs ('Is this a lattice for hex-close-packed spheres?:$')
      call res(ans,1)
c
      if (ans.eq.'n'.or.ans.eq.'N') return
      call rzero(x0,3)
      dlat = 1.
      call get_vert_tet(q,dlat,x0) ! Base tet
c
      if_lattice = .true.                   ! Store unit lattice vectors
      call sub3   (ulat1,q(1,2),q(1,1),3)   ! -- needed for periodic bcs
      call norm3d (ulat1)
      call sub3   (ulat2,q(1,3),q(1,1),3)
      call norm3d (ulat2)
      call sub3   (ulat3,q(1,4),q(1,1),3)
      call norm3d (ulat3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine shift_xyz(p,ai,n)
c
      real p(3,n),ai(3)
      common /ashift/ a(3)
c
      a(1) = ai(1) ! Avoid aliasing
      a(2) = ai(2)
      a(3) = ai(3)
c
      do i=1,n
         p(1,i) = p(1,i) - a(1)
         p(2,i) = p(2,i) - a(2)
         p(3,i) = p(3,i) - a(3)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat(a,m,n,name5,k)
      real a(m,n)
      character*5 name5
c
      write(6 ,1) name5,k
      write(66,1) name5,k
      n10 = min(n,10)
      do i=1,m
         write(6 ,2) k,name5,(a(i,j),j=1,n10)
         write(66,2) k,name5,(a(i,j),j=1,n10)
      enddo
    1 format(/,'MATRIX: ',a5,i9)
    2 format(i4,1x,a5,1p10e11.3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sdl_element_dodec(v,xrhm,rad)
c
c     Update SEM element data for dodec
c
      include 'basics.inc'
      real v(3,8,4),xrhm(3,4)
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
      do irhm=1,4
         k = k+1
         e = nel+k
         do i=1,8
            j      = e2pfv(i)
            x(i,e) = v(1,j,k)
            y(i,e) = v(2,j,k)
            z(i,e) = v(3,j,k)
         enddo
c
         ccurve(5,e)   = ' '
         call rzero(curve(1,1,e),6)
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
c
      nel = nel+4
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rhombic_dodec(d,x0,if_sph_ctr)
c
      logical if_sph_ctr
c
c
c     x(:,1)=x0, and x(:,2)=x0+d
c
      real x(3,8),x0(3)
      real t(3,4),x00(3),v(3,4),ve(3,8,4)
c
      call get_vert_tet(t,d,x0) ! Base tet
c
c     Get unit vectors emanating from origin of tet
c
      call rzero(x00,3)
      do k=1,4
      do i=1,3
         x00(i) = x00(i) + t(i,k)
      enddo
      enddo
      scale = 0.25
      call cmult(x00,scale,3)
c
      do k=1,4
         do i=1,3
            v(i,k) = t(i,k)-x00(i)
         enddo
         v(3,k) = -v(3,k)   ! Need this for right-handedness of elements
         call norm3d (v(1,k))
      enddo
c
      call rzero(ve,3*8*4)
c
      call rhomboid(ve(1,1,1),v(1,4),v(1,3),v(1,2))
      call rhomboid(ve(1,1,2),v(1,3),v(1,4),v(1,1))
      call rhomboid(ve(1,1,3),v(1,2),v(1,1),v(1,4))
      call rhomboid(ve(1,1,4),v(1,1),v(1,2),v(1,3))
c
      call sdl_element_dodec  (ve,x,rad)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rhomboid(ve,v1,v2,v4)
c
      real ve(3,0:7),v1(3),v2(3),v4(3)
c
      vol = triple_prod(v1,v2,v4) ! Check for positive volume
c
      call rzero(ve,24)
      call add3(ve(1,1),ve(1,0),v1,3)
      call add3(ve(1,2),ve(1,0),v2,3)
      call add3(ve(1,3),ve(1,1),v2,3)
c
c
      call add3(ve(1,4),ve(1,0),v4,3)
      call add3(ve(1,5),ve(1,1),v4,3)
      call add3(ve(1,6),ve(1,2),v4,3)
      call add3(ve(1,7),ve(1,3),v4,3)
c
      return
      end
c-----------------------------------------------------------------------
      function tet_chk(q)
c
c     Check for positive volume (i.e., right-handed tet)
c
      real q(3,0:3)
      real v(3,3)
c
c
      do i=1,3
         call sub3(v(1,i),q(1,i),q(1,0),3)
      enddo
c
      vol = triple_prod(v(1,1),v(1,2),v(1,3))
      if (vol.le.0) write(6,*) 'in tet chk'
c
      tet_chk = vol
c
      return
      end
c-----------------------------------------------------------------------
      function triple_prod(v1,v2,v3)
c
      real v1(3),v2(3),v3(3)
      real w1(3)
c
c     vol = triple_prod(v1,v2,v3) ! Check for positive volume
c
      call cross(w1,v1,v2)
      vol = dot(w1,v3,3)
c
      if (vol.le.0) then
         write(6,1) vol
         write(6,2) 'V1:',(v1(k),k=1,3)
         write(6,2) 'V2:',(v2(k),k=1,3)
         write(6,2) 'V3:',(v3(k),k=1,3)
      endif
    1 format(/,'Nonpositive volume:',1pe12.4)
    2 format(a3,1p3e14.4)
c
      triple_prod = vol
c
      return
      end
c-----------------------------------------------------------------------
      subroutine vswap2(x,y,n)
      real x(1),y(1)
c
      do i=1,n
         t    = x(i)
         x(i) = y(i)
         y(i) = t
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine edge_tet_sdl(v,nt,p,i0,i1,k0,k1,d,rad)
c
c     Build the void for an fcc lattice
c
      real p(3,8,2) ! sphere centers 8 + 6
      real q(3,6)   ! permuted diamond sphere centers 6
c
      real v(3,8,288)  ! 288 elements total
c
c
c     Edge tets 3 x 4 = 12 total
c
      call copy(q(1,1),p(1,i0,1),3)  ! Tet 1
      call copy(q(1,2),p(1,i1,1),3)
      call copy(q(1,3),p(1,k0,2),3)
      call copy(q(1,4),p(1,k1,2),3)
c
      vol = tet_chk(q)
      if (vol.lt.0) write(6,*) 'in edge_tet_sdl',vol
      if (vol.lt.0) call vswap2(q(1,1),q(1,2),3)
c
      call build_tet_saddle(v(1,1,nt),q,d,rad,.false.)
      call sdl_element_tet (v(1,1,nt),q,rad)
      nt = nt + 12
c
      return
      end
c-----------------------------------------------------------------------
      subroutine crnr_tet_sdl(v,nt,p,i0,k0,k1,k2,d,rad)
c
c     Build the void for an fcc lattice
c
      real p(3,8,2) ! sphere centers 8 + 6
      real q(3,6)   ! permuted diamond sphere centers 6
c
      real v(3,8,288)  ! 288 elements total
c
c
c     Corner tets: 8 total
c
      call copy(q(1,1),p(1,i0,1),3)  ! Tet 1
      call copy(q(1,2),p(1,k0,2),3)
      call copy(q(1,3),p(1,k1,2),3)
      call copy(q(1,4),p(1,k2,2),3)
c
      vol = tet_chk(q)
      if (vol.lt.0) call vswap2(q(1,1),q(1,2),3)
c
      call build_tet_saddle(v(1,1,nt),q,d,rad,.false.)
      call sdl_element_tet (v(1,1,nt),q,rad)
      nt = nt + 12
c
      return
      end
c-----------------------------------------------------------------------
      subroutine saddle_fcc(d,radi,x0)
c
c     Build the void for an fcc lattice
c
      include 'basics.inc'
      real x0(3)
      real p(3,8,2) ! sphere centers 8 + 6
      real q(3,6)   ! permuted diamond sphere centers 6
c
      real v(3,8,288)  ! 48 elements total (12 + 24 + 12)
      real a0(3)
c
      two = 2.
      rad = radi/sqrt(two)
c
      call get_vert_fcc_lat (p,d,x0)
c
      nt = 1
c
c     Edge tets 3 x 4 = 12 total
c
      call edge_tet_sdl(v,nt,p,1,2,3,5,d,rad)
c     return
      call edge_tet_sdl(v,nt,p,3,4,4,5,d,rad)
      call edge_tet_sdl(v,nt,p,5,6,3,6,d,rad)
      call edge_tet_sdl(v,nt,p,7,8,4,6,d,rad)
c
      call edge_tet_sdl(v,nt,p,1,3,1,5,d,rad)
      call edge_tet_sdl(v,nt,p,2,4,2,5,d,rad)
      call edge_tet_sdl(v,nt,p,5,7,1,6,d,rad)
      call edge_tet_sdl(v,nt,p,6,8,2,6,d,rad)
c
      call edge_tet_sdl(v,nt,p,1,5,1,3,d,rad)
      call edge_tet_sdl(v,nt,p,2,6,2,3,d,rad)
      call edge_tet_sdl(v,nt,p,3,7,1,4,d,rad)
      call edge_tet_sdl(v,nt,p,4,8,2,4,d,rad)
c
c
c     BUILD CORNER TETS  (8 total)
c
c
      call crnr_tet_sdl(v,nt,p,1,1,3,5,d,rad)
      call crnr_tet_sdl(v,nt,p,2,2,3,5,d,rad)
      call crnr_tet_sdl(v,nt,p,3,1,4,5,d,rad)
      call crnr_tet_sdl(v,nt,p,4,2,4,5,d,rad)
      call crnr_tet_sdl(v,nt,p,5,1,3,6,d,rad)
      call crnr_tet_sdl(v,nt,p,6,2,3,6,d,rad)
      call crnr_tet_sdl(v,nt,p,7,1,4,6,d,rad)
      call crnr_tet_sdl(v,nt,p,8,2,4,6,d,rad)
c
      call permute_dia      (q,p(1,1,2))  ! Diamond
      call build_dia_saddle (v(1,1,nt),q,d,rad,if_sph_ctr)
      call sdl_element_dia  (v(1,1,nt),q,rad)
      nt = nt+24
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine add3s2(x,y,z,c,n)
      real x(1),y(1),z(1)
      do i=1,n
         x(i) = c * ( y(i) + z(i) )
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)

c
c     Return face index, preprocessor notation
c

      kx1=1
      ky1=1
      kz1=1
      kx2=nx1
      ky2=ny1
      kz2=nz1

      if (iface.eq.1) ky2=1
      if (iface.eq.2) kx1=nx1
      if (iface.eq.3) ky1=ny1
      if (iface.eq.4) kx2=1
      if (iface.eq.5) kz2=1
      if (iface.eq.6) kz1=nz1

      return
      end
c-----------------------------------------------------------------------
      subroutine sfacind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)

c
c     Return face index, symmetric notation
c

      kx1=1
      ky1=1
      kz1=1
      kx2=nx1
      ky2=ny1
      kz2=nz1

      if (iface.eq.1) kx2=1
      if (iface.eq.2) kx1=nx1
      if (iface.eq.3) ky2=1
      if (iface.eq.4) ky1=ny1
      if (iface.eq.5) kz2=1
      if (iface.eq.6) kz1=nz1

      return
      end
c-----------------------------------------------------------------------
      subroutine non_reg_tet_mod(v,x1,x2,x3,r)
      real v(3),x1(3),x2(3),x3(3)
      real w(3),x(3,3),s(3),a(3,3),xr(3,3),sl(3),t(3)
c
      call copy(x(1,1),x1,3)
      call copy(x(1,2),x2,3)
      call copy(x(1,3),x3,3)
c
      h = .5
      call add3s2(a(1,1),x2,x3,h,3)
      call add3s2(a(1,2),x3,x1,h,3)
      call add3s2(a(1,3),x1,x2,h,3)
c
      call copy(w,v,3)

      do i=1,3
         call sub3(s,a(1,i),x(1,i),3)
         sli = sqrt(vlsc2(s,s,3)) - r
         slr = sli - r
         call sub3(t,v,x(1,i),3)
         sll   = sqrt(vlsc2(t,t,3)) - r
         ds    = 0.5*slr - sll
c
         ds  = -ds/sli
         call add2s2(w,s,ds,3)
      enddo
c
      call copy(v,w,3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sph_intersect(xi,x0,x1,xc,r)
      real xi(3),x0(3),x1(3),xc(3),r
      real p0(3),p1(3),dx(3)

c     xi:   output point where [x0,x1] intersects sphere (xc,r)
c
c     x0    is assumed inside the sphere
c     x1    is assumed outside the sphere
c

      call sub3(p0,x0,xc,3)   ! subtract center
      call sub3(p1,x1,xc,3)

      scale = 1./r
      call cmult(p0,scale,3)  ! scale out radius
      call cmult(p1,scale,3)

      call sub3(dx,p1,p0,3)   ! compute unit vector pointing x0 --> x1
      scale = 1./dot(dx,dx,3)
      scale = sqrt(scale)
      call cmult(dx,scale,3)

      c = dot(p0,p0,3) - 1.   ! find intersection w/ unit sphere
      b = 2.*dot(p0,dx,3)
      d = b*b-4*c
      s = sqrt(d)
      ap = .5*(-b + s)        ! take positive root
      call copy  (xi,p0,3)
      call add2s2(xi,dx,ap,3)

      scale = r
      call cmult (xi,scale,3)  ! rescale by radius
      call add2  (xi,xc   ,3)  ! add back sphere center

c     write(6,*)
c     write(6,*) 'this is rad:',r
c     call outmat(x0,1,3,'xsph0',1)
c     call outmat(x1,1,3,'xsph1',2)
c     call outmat(xi,1,3,'xsph ',3)

      return
      end
c-----------------------------------------------------------------------
      subroutine nw_sphmesh

      call sph_wall_elmt
      call redraw_mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine sph_wall_elmt
c
c     Place a sphere near a wall;   10/15/06  pff

      include 'basics.inc'
      real xvi(3),zvi(5),tsph(3)

c     tsph -- translation of sphere from base position

      real xt(3,25,6,4)  ! base coordinates
      real x0  (3,3)       ! sphere centers
      real rads(0:3)       ! sphere centers

      character*80 fname

      nelo = nel

      call prs('Input file name containing r0,h0,xvi,zvi,tsph:$')
      call blank(fname,80)
      call res  (fname,80)

      open(unit=80,file=fname)
      read(80,*) r0,h0,(xvi(k),k=1,3),(zvi(k),k=1,5),(tsph(k),k=1,3)
      read(80,*) s1,s2,s3
      close(unit=80)

c     r0 = 0.5
c     h0 = 0.505
   
c     xvi(1) = 0
c     xvi(2) = 0.4
c     xvi(3) = 0.8
   
c     zvi(1) = 0
c     zvi(2) = 0.2
c     zvi(3) = 0.4
c     zvi(4) = 0.8
c     zvi(5) = 1.3

c     tsph(1) = 0.
c     tsph(2) = 0.
c     tsph(3) = 0.

      rm = min(xvi(3),zvi(5)-h0)
      dm = rm - r0

c     r1 = r0 + 0.10*dm          ! 1st shell
c     r2 = r0 + 0.65*dm          ! 2nd shell
c     r3 = r0 + 1.00*dm          ! 3rd shell

      r1 = r0 + s1*dm          ! 1st shell
      r2 = r0 + s2*dm          ! 2nd shell
      r3 = r0 + s3*dm          ! 3rd shell

      rads(0) = r0
      rads(1) = r1
      rads(2) = r2
      rads(3) = r3


      ntp  = 5
      call sph_wall_elmt1(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)
      call sph_wall_elmt2(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)
      call sph_wall_elmt3(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)

      call translate_sub_mesh(nelo,nel,tsph(1),tsph(2),tsph(3))
      write(6,*) 'this is nel:',nelo,nel

      return
      end
c-----------------------------------------------------------------------
      subroutine sph_w_gt_fc(xt,ntp,nlev,pface)
c
c     Place a sphere near a wall;   10/15/06  pff
c
c     Top & Bottom parts
c

      real xt(3,ntp,ntp,6,nlev)  ! base coordinates

      integer pface,f

      if (pface.eq.6) then

         k = ntp
         do ilev = 1,nlev
         do idir = 1,3

            f = 1
            i = ntp
            do j=1,ntp
               xt(idir,i,j,6,ilev) = xt(idir,j,k,f,ilev)
            enddo

            f = 2
            j = ntp
            l = 0
            do i=ntp,1,-1
               l = l+1
               xt(idir,i,j,6,ilev) = xt(idir,l,k,f,ilev)
            enddo
      
            f = 3
            i = 1
            l = 0
            do j=ntp,1,-1
               l = l+1
               xt(idir,i,j,6,ilev) = xt(idir,l,k,f,ilev)
            enddo
      
            f = 4
            j = 1
            do i=1,ntp
               xt(idir,i,j,6,ilev) = xt(idir,i,k,f,ilev)
            enddo
      
         enddo
         enddo

      elseif (pface.eq.5) then

         k = 1
         do ilev = 1,nlev
         do idir = 1,3

            f = 1
            i = ntp
            l = 0
            do j=ntp,1,-1
               l = l+1
               xt(idir,i,j,5,ilev) = xt(idir,l,k,f,ilev)
            enddo

            f = 2
            j = 1
            l = 0
            do i=ntp,1,-1
               l = l+1
               xt(idir,i,j,5,ilev) = xt(idir,l,k,f,ilev)
            enddo
      
            f = 3
            i = 1
            do j=1,ntp
               xt(idir,i,j,5,ilev) = xt(idir,j,k,f,ilev)
            enddo
      
            f = 4
            j = ntp
            do i=1,ntp
               xt(idir,i,j,5,ilev) = xt(idir,i,k,f,ilev)
            enddo
      
         enddo
         enddo
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine sph_wall_elmt1(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)
c
c     Place a sphere near a wall;   10/15/06  pff
c
      include 'basics.inc'

      real r0,h0           ! radius and height of sphere from wall
      real xvi(3),zvi(5)   ! local coords of box vertices
      real rads(0:3)       ! location of sphere centers
      real tsph(3)         ! x-y location of sphere, after translation

      integer e,f
      integer pf2ev(8)
      save    pf2ev
      data    pf2ev / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      real xt(3,25,6,4)  ! base coordinates
      real x0(3,3)
      real x2(3)
      real udata(3)

      integer pface


      r1 = rads(1)             ! 1st shell
      r2 = rads(2)             ! 2nd shell
      r3 = rads(3)             ! 3rd shell

      udata(1) = r0            ! this data passed to prenek as part
      udata(2) = h0            ! of the 'u' curve side parameter list
      udata(3) = r1

      h1 = .5*(h0-r0) + r1
      h1 = max(h1,h0)          ! outer sphere not closer then inner
      h2 = h1
      h3 = h2                  ! 3rd shell


      call rzero(x0,3*3)
      x0(3,1) = h0
      x0(3,2) = h1
      x0(3,3) = h2
      call rzero(x2,3)
      x2(1)   = -r2*1
      x2(3)   = h2
      r22     = r2*2

      nt2 = ntp/2 + 1

      call rzero(xt,3*ntp*ntp*6*4)

      pface = 1

      dy = r3/2.
      dz = r3/2.

      m = 0
      do k=1,ntp
      do j=1,ntp

         kk = k-nt2
         jj = j-nt2
         m  = m+1
         xt(1,m,pface,4) = r3
         xt(2,m,pface,4) = jj*dy
         xt(3,m,pface,4) = kk*dz + h3

         do kl = 1,3

            if (kl.eq.1) call sph_intersect
     $        (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x0(1,kl),r0)
            if (kl.eq.2) call user_s
     $        (xt(1,m,pface,kl),xt(1,m,pface,4),udata,1,5)
            if (kl.eq.3) call sph_intersect
     $        (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x2      ,r22)

            xt(3,m,pface,kl) = max(xt(3,m,pface,kl),0.)
            if (kl.ge.3.and.k.eq.1) xt(3,m,pface,kl) = 0.

         enddo

c        adjust x-pos on outer shell to user position, post-projection
         xt(1,m,pface,4) = xvi(nt2)
         ja = abs(jj)+1
         if (jj.lt.0) xt(2,m,pface,4) = -xvi(ja)
         if (jj.gt.0) xt(2,m,pface,4) =  xvi(ja)
         if (jj.eq.0) xt(2,m,pface,4) =  0

c        adjust z-level on outer shell to user position, post-projection
         xt(3,m,pface,4) = zvi(k)

c        floor z-level to zero
         xt(3,m,pface,4) = max(xt(3,m,pface,4),0.)
         if (k.eq.1) xt(3,m,pface,4) = 0.


      enddo
      enddo

      do f=2,4   !  save tensor-product brick

         if (f.eq.2) then     ! ang = 90
            ca = 0
            sa = 1
         elseif (f.eq.3) then ! ang = 180
            ca = -1
            sa = 0
         elseif (f.eq.4) then ! ang = 270
            ca = 0
            sa = -1
         endif

         do k = 1,4
         do m=1,ntp*ntp
               xx = xt(1,m,pface,k)
               yy = xt(2,m,pface,k)
               zz = xt(3,m,pface,k)
               xt(1,m,f,k) = ca*xx-sa*yy
               xt(2,m,f,k) = sa*xx+ca*yy
               xt(3,m,f,k) =    zz
         enddo
         enddo
      enddo

      
      nsh = 0
      do klev =1,3   ! fill all shells w/ solid

         m = 0
         do ke=1,ntp-1
         do je=1,ntp-1

            nsh = nsh+1
            e   = nel + nsh

            j = 0
            do kk=0,1
            do jj=0,1
            do ii=0,1

               m = je + ii + ntp*(ke-1 + jj)
               k = klev + kk

               j      = j+1
               i      = pf2ev(j)
               x(i,e) = xt(1,m,pface,k)
               y(i,e) = xt(2,m,pface,k)
               z(i,e) = xt(3,m,pface,k)

            enddo
            enddo
            enddo

            call rzero(curve(1,1,e),72)

            if (klev.eq.1) then
               ccurve(  5,e) = 's'
               curve (1,5,e) = x0(1,1)
               curve (2,5,e) = x0(2,1)
               curve (3,5,e) = x0(3,1)
               curve (4,5,e) = r0

               ccurve(  6,e) = 'u'
               curve (1,6,e) = udata(1)
               curve (2,6,e) = udata(2)
               curve (3,6,e) = udata(3)
               curve (4,6,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,6,e) = tsph(2)   ! y-loc of translated sphere
            elseif (klev.eq.2) then
               ccurve(  5,e) = 'u'
               curve (1,5,e) = udata(1)
               curve (2,5,e) = udata(2)
               curve (3,5,e) = udata(3)
               curve (4,5,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,5,e) = tsph(2)   ! y-loc of translated sphere
            endif

            do f=1,6
               cbc(f,e,1) = 'v  '
            enddo

            cbc   (5,e,1) = 'v  '
            cbc   (6,e,1) = 'v  '
            cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature

            call rzero(bc(1,5,e,1),5)
            call rzero(bc(1,5,e,2),5)

            ilet      = mod1(e,52)
            letapt(e) = alphabet(ilet)
            numapt(e) = 1
         enddo
         enddo
c
      enddo


      call copy_sub_mesh    (nel+1,nel+nsh,nel+nsh+1)
      call rotate_submesh_2d(nel+nsh+1,nel+2*nsh,90.)
      nsh = 2*nsh
      call copy_sub_mesh    (nel+1,nel+nsh,nel+nsh+1)
      call rotate_submesh_2d(nel+nsh+1,nel+2*nsh,180.)
      nsh = 2*nsh

      nel = nel+nsh

      return
      end
c-----------------------------------------------------------------------
      subroutine sph_wall_elmt2(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)
c
c     Place a sphere near a wall;   10/15/06  pff
c
c     Top & Bottom parts
c
      include 'basics.inc'

      real r0,h0           ! radius and height of sphere from wall
      real xvi(3),zvi(5)   ! local coords of box vertices
      real rads(0:3)       ! location of sphere centers
      real tsph(3)         ! x-y location of sphere, after translation

      integer e,f
      integer pf2ev(8)
      save    pf2ev
      data    pf2ev / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      real xt(3,ntp*ntp,6,4)  ! base coordinates
      real x0(3,3)
      real x2(3)

      integer pface
      real udata(3)

      r1 = rads(1)             ! 1st shell
      r2 = rads(2)             ! 2nd shell
      r3 = rads(3)             ! 3rd shell

      h1 = .5*(h0-r0) + r1
      h1 = max(h1,h0)          ! outer sphere not closer then inner
      h2 = h1
      h3 = h2                  ! 3rd shell

      udata(1) = r0
      udata(2) = h0
      udata(3) = r1

      call rzero(x2,3)
      x2(3)   = h2-r2
      r22     = r2*2

      nt2 = ntp/2 + 1

      pface = 6

      dx = r3/2.
      dy = r3/2.

      m = 0
      do j=1,ntp
      do i=1,ntp

         jj = j-nt2
         ii = i-nt2
         m  = m+1
         xt(1,m,pface,4) = ii*dx
         xt(2,m,pface,4) = jj*dy
         xt(3,m,pface,4) = h3 + r3

         do kl = 1,3
            if (kl.eq.1) call sph_intersect
     $        (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x0(1,kl),r0)
            if (kl.eq.2) call user_s
     $        (xt(1,m,pface,kl),xt(1,m,pface,4),udata,1,5)
c           if (kl.eq.2) call sph_intersect
c    $        (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x0(1,kl),r1)
            if (kl.eq.3) call sph_intersect
     $        (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x2      ,r22)
         enddo

c        adjust xy-pos on outer shell to user position, post-projection

         if (ii.gt.0) then
            xt(1,m,pface,4) = xvi(ii+1)
         elseif (ii.lt.0) then
            ia = abs(ii)
            xt(1,m,pface,4) = -xvi(ia+1)
         else
            xt(1,m,pface,4) = 0
         endif

         if (jj.gt.0) then
            xt(2,m,pface,4) = xvi(jj+1)
         elseif (jj.lt.0) then
            ja = abs(jj)
            xt(2,m,pface,4) = -xvi(ja+1)
         else
            xt(2,m,pface,4) = 0
         endif
         xt(3,m,pface,4) = zvi(ntp)

      enddo
      enddo

      nlev = 4
      call sph_w_gt_fc(xt,ntp,nlev,pface)

      nsh = 0
      do klev =1,3   ! fill all shells w/ solid

         m = 0
         do ke=1,ntp-1
         do je=1,ntp-1

            nsh = nsh+1
            e   = nel + nsh

            j = 0
            do kk=0,1
            do jj=0,1
            do ii=0,1

               m = je + ii + ntp*(ke-1 + jj)
               k = klev + kk

               j      = j+1
               i      = pf2ev(j)
               x(i,e) = xt(1,m,pface,k)
               y(i,e) = xt(2,m,pface,k)
               z(i,e) = xt(3,m,pface,k)

            enddo
            enddo
            enddo

            call rzero(curve(1,1,e),72)

            if (klev.eq.1) then
               ccurve(  5,e) = 's'
               curve (1,5,e) = x0(1,1)
               curve (2,5,e) = x0(2,1)
               curve (3,5,e) = x0(3,1)
               curve (4,5,e) = r0

               ccurve(  6,e) = 'u'
               curve (1,6,e) = udata(1)
               curve (2,6,e) = udata(2)
               curve (3,6,e) = udata(3)
               curve (4,6,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,6,e) = tsph(2)   ! y-loc of translated sphere
            elseif (klev.eq.2) then
               ccurve(  5,e) = 'u'
               curve (1,5,e) = udata(1)
               curve (2,5,e) = udata(2)
               curve (3,5,e) = udata(3)
               curve (4,5,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,5,e) = tsph(2)   ! y-loc of translated sphere
            endif

            do f=1,6
               cbc(f,e,1) = 'v  '
            enddo

            cbc   (5,e,1) = 'v  '
            cbc   (6,e,1) = 'v  '
            cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature

            call rzero(bc(1,5,e,1),5)
            call rzero(bc(1,5,e,2),5)

            ilet      = mod1(e,52)
            letapt(e) = alphabet(ilet)
            numapt(e) = 1
         enddo
         enddo
c
      enddo

      nel = nel+nsh

      return
      end
c-----------------------------------------------------------------------
      subroutine sph_wall_elmt3(r0,h0,xvi,zvi,xt,x0,rads,ntp,tsph)
c
c     Place a sphere near a wall;   10/15/06  pff
c
c     Bottom parts
c
      include 'basics.inc'

      real r0,h0           ! radius and height of sphere from wall
      real xvi(3),zvi(5)   ! local coords of box vertices
      real rads(0:3)       ! location of sphere centers
      real tsph(3)         ! x-y location of sphere, after translation

      integer e,f
      integer pf2ev(8)
      save    pf2ev
      data    pf2ev / 1 , 2 ,4 ,3 ,5 ,6 ,8 ,7 /
c
      character*1 alphabet(52)
      character*26 alpha(2)
      equivalence (alphabet,alpha)
      save         alpha
      data         alpha 
     $      /'abcdefghijklmnopqrstuvwxyz','ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      real xt(3,ntp*ntp,6,4)  ! base coordinates
      real x0(3,3)
      real x2(3)

      integer pface
      real    udata(3)

      r1 = rads(1)             ! 1st shell
      r2 = rads(2)             ! 2nd shell
      r3 = rads(3)             ! 3rd shell

      h1 = .5*(h0-r0) + r1
      h1 = max(h1,h0)          ! outer sphere not closer then inner
      h2 = h1
      h3 = h2                  ! 3rd shell

      udata(1) = r0
      udata(2) = h0
      udata(3) = r1

      nt2 = ntp/2 + 1

      pface = 5

      dx  = r3/2.  ! to correspond to other sections
      dy  = r3/2.

      dxx = r2/3.
      dyy = r2/3.

      m = 0
      do j=1,ntp
      do i=1,ntp

         jj = nt2-j
         ii = i-nt2
         m  = m+1
         xt(1,m,pface,3) = ii*dxx
         xt(2,m,pface,3) = jj*dyy
         xt(3,m,pface,3) = 0

         xt(1,m,pface,4) = ii*dx
         xt(2,m,pface,4) = jj*dy
         xt(3,m,pface,4) = h3 - r3

         ia = abs(ii)
         ja = abs(jj)

         do kl = 1,2
            if (ia.le.1 .and. ja.le.1) then
               if (kl.eq.1) call sph_intersect
     $           (xt(1,m,pface,kl),x0,xt(1,m,pface,3),x0(1,kl),r0)
c              if (kl.eq.2) call sph_intersect
c    $           (xt(1,m,pface,kl),x0,xt(1,m,pface,3),x0(1,kl),r1)
               if (kl.eq.2) call user_s
     $           (xt(1,m,pface,kl),xt(1,m,pface,3),udata,1,5)
            else
               if (kl.eq.1) call sph_intersect
     $           (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x0(1,kl),r0)
c              if (kl.eq.2) call sph_intersect
c    $           (xt(1,m,pface,kl),x0,xt(1,m,pface,4),x0(1,kl),r1)
               if (kl.eq.2) call user_s
     $           (xt(1,m,pface,kl),xt(1,m,pface,4),udata,1,5)
            endif
         enddo

      enddo
      enddo

      nlev = 3
      call sph_w_gt_fc(xt,ntp,nlev,pface)

      nsh = 0
      do klev =1,2   ! fill all shells w/ solid

         m = 0
         do ke=1,ntp-1
         do je=1,ntp-1

            nsh = nsh+1
            e   = nel + nsh

            j = 0
            do kk=0,1
            do jj=0,1
            do ii=0,1

               m = je + ii + ntp*(ke-1 + jj)
               k = klev + kk

               j      = j+1
               i      = pf2ev(j)
               x(i,e) = xt(1,m,pface,k)
               y(i,e) = xt(2,m,pface,k)
               z(i,e) = xt(3,m,pface,k)

            enddo
            enddo
            enddo

            call rzero(curve(1,1,e),72)

            if (klev.eq.1) then
               ccurve(  5,e) = 's'
               curve (1,5,e) = x0(1,1)
               curve (2,5,e) = x0(2,1)
               curve (3,5,e) = x0(3,1)
               curve (4,5,e) = r0

               ccurve(  6,e) = 'u'
               curve (1,6,e) = udata(1)
               curve (2,6,e) = udata(2)
               curve (3,6,e) = udata(3)
               curve (4,6,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,6,e) = tsph(2)   ! y-loc of translated sphere
            elseif (klev.eq.2) then
               ccurve(  5,e) = 'u'
               curve (1,5,e) = udata(1)
               curve (2,5,e) = udata(2)
               curve (3,5,e) = udata(3)
               curve (4,5,e) = tsph(1)   ! x-loc of translated sphere
               curve (5,5,e) = tsph(2)   ! y-loc of translated sphere
            endif

            do f=1,6
               cbc(f,e,1) = 'v  '
            enddo

            cbc   (5,e,1) = 'v  '
            cbc   (6,e,1) = 'v  '
            cbc   (5,e,2) = 'f  '  ! flux bc as default for Temperature

            call rzero(bc(1,5,e,1),5)
            call rzero(bc(1,5,e,2),5)

            ilet      = mod1(e,52)
            letapt(e) = alphabet(ilet)
            numapt(e) = 1
         enddo
         enddo
c
      enddo

      nel = nel+nsh

      return
      end
c-----------------------------------------------------------------------
      subroutine sph_combine_d(xi,r0,x0,x1,xs1)
      real xi(3),x0(3),x1(3),xs1(0:3)

      real t1(3),t2(3),un(3),d1(3),d2(3)

c     xi:   output point where [x0,x1] intersects sphere (xc,r)
c
c     x0    is assumed inside the sphere
c     x1    is assumed outside the sphere
c

      call sph_intersect(t1,x0,x1,xs1(1),xs1(0))

      h0 = x0(3)

      call copy(xi,t1,3)
      if (x1(3).ge.h0) return

      call sub3(t2,t1,x0,3)
      call normalize(t2,alpha,3)

   
      call rzero(un,3) 
      un(3) = -1.
      cos_th = dot(t2,un,3)

      r = 0.5*(r0 + h0/cos_th)

      call add2s1(t2,x0,r,3)
      

c     Blend

      call sub3(d1,t1,x0,3)
      d12 = vlsc2(d1,d1,3)
      d12 = sqrt(d12)

      call sub3(d2,t2,x0,3)
      d22 = vlsc2(d2,d2,3)
      d22 = sqrt(d22)

      cos2   = cos_th**2
      sin2   = 1-cos2
      dm = min(d12,d22)
      w1 = dm/d12
      w2 = dm/d22
      a  = 12.
      w1 = sin2*(w1**a)
      w2 = cos2*(w2**a)
      ww = w1+w2
      w1 = w1/ww
      w2 = w2/ww
      do i=1,3
         xi(i) = w1*t1(i) + w2*t2(i)
      enddo
      write(6,1) (xi(k),k=1,3),(x1(k),k=1,3)
      write(29,1) (xi(k),k=1,3)
1     format(2(2x,3f11.5))

      return
      end
c-----------------------------------------------------------------------
      subroutine user_s(xu,xf,xo,eg,face)

      real xu(3),xf(3),xo(3)
      integer eg,face

      real x0(3),x1(0:3)    ! sphere center
      save x0,x1
      data x0,x1 / 7*0./

      
c     r0 = .5      ! sphere def'n
c     h0 = .505

      r0 = xo(1)   ! hack into user_s interface
      h0 = xo(2)
      r1 = xo(3)

      write(6,*) r0,h0,r1,' R0'

      x0 (1) = 0.
      x0 (2) = 0.
      x0 (3) = h0

c     r1 = .53
      h1 = r1 + (h0-r0)/2.
      h1 = h0
      x1(0) = r1
      x1(1) = 0.
      x1(2) = 0.
      x1(3) = h1

c     write(6,1) (xf(k),k=1,3)
c   1 format(6f11.3)

      call sph_combine_d(xu,r0,x0,xf,x1)

      return
      end
c-----------------------------------------------------------------------
      subroutine permute_xyz (x,y,z,n)
      real x(n),y(n),z(n)

      do i=1,n

         tx=x(i)
         ty=y(i)
         tz=z(i)

         y(i) = tx
         z(i) = ty
         x(i) = tz

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine el_convert (e,xr,yr,zr)
      include 'basics.inc'
      integer e
      real xr(8),yr(8),zr(8)

      write(6,*)
      do i=1,8
         x(i,e) = xr(i)
         y(i,e) = yr(i)
         z(i,e) = zr(i)
c        write(6,1) e,xr(i),yr(i),zr(i),i
c  1     format(i8,3f9.4,i3,' ec')
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine el_convert_2to3(e,xr,yr,zr)
      include 'basics.inc'
      integer e
      real xr(8),yr(8),zr(8)

      write(6,*)
      do i=1,8
         i2d = mod1(i,4)
         x(i,e) = xr(i2d)
         y(i,e) = yr(i2d)
         z(i,e) = zr(i)
c        write(6,1) e,xr(i),yr(i),zr(i),i
c  1     format(i8,3f9.4,i3,' ec2')
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_central_shell(r1,r2,e)

      include 'basics.inc'
      integer e
      common /xyzr/ xr(100),yr(100),zr(100)

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      rr = sqrt(r1*r2) ! geometric mean
      rr = .5*(sqrt(r1*r2)+r1)
      x1 = r1
      x2 = rr/s2
      x3 = rr/s3

      xr(2)=x1
      xr(3)=x2
      xr(6)=x2
      xr(7)=x3
      yr(3)=x2
      yr(4)=x1
      yr(7)=x3
      yr(8)=x2
      zr(5)=x1
      zr(6)=x2
      zr(7)=x3
      zr(8)=x2

      call el_convert(e,xr,yr,zr)

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_transition_shell(r1,r2,e) ! build 3 elements, use rotation

      include 'basics.inc'
      integer e,e1,e2
      common /xyzr/ xr(100),yr(100),zr(100)

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      rr = sqrt(r1*r2) ! geometric mean
      rr = .5*(sqrt(r1*r2)+r1)
      p1 = r1
      p2 = rr/s2
      p3 = rr/s3

      q1 = r2
      q2 = r2/s2
      q3 = r2/s3


      zr(1)=p1
      zr(2)=p2
      zr(3)=p3
      zr(4)=p2
      zr(5)=q1
      zr(6)=q2
      zr(7)=q3
      zr(8)=q2

      xr(2)=p2
      xr(3)=p3
      xr(6)=q2
      xr(7)=q3

      yr(3)=p3
      yr(4)=p2
      yr(7)=q3
      yr(8)=q2

      call el_convert  (e,xr,yr,zr)

      e1 = e+1
      e2 = e+2

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e1,xr,yr,zr)

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e2,xr,yr,zr)

      ccurve(6,e ) = 's'
      ccurve(6,e1) = 's'
      ccurve(6,e2) = 's'

      call rzero(curve(1,1,e),3*72)
      curve(4,6,e ) = r2
      curve(4,6,e1) = r2
      curve(4,6,e2) = r2

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_std_shell(r1,r2,e) ! build 3 elements, use rotation

      include 'basics.inc'
      integer e,e1,e2
      common /xyzr/ xr(100),yr(100),zr(100)

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      p1 = r1
      p2 = r1/s2
      p3 = r1/s3

      q1 = r2
      q2 = r2/s2
      q3 = r2/s3

      zr(1)=p1
      zr(2)=p2
      zr(3)=p3
      zr(4)=p2
      zr(5)=q1
      zr(6)=q2
      zr(7)=q3
      zr(8)=q2

      xr(2)=p2
      xr(3)=p3
      xr(6)=q2
      xr(7)=q3

      yr(3)=p3
      yr(4)=p2
      yr(7)=q3
      yr(8)=q2

      call el_convert  (e,xr,yr,zr)

      e1 = e+1
      e2 = e+2

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e1,xr,yr,zr)

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e2,xr,yr,zr)

      ccurve(5,e ) = 's'
      ccurve(5,e1) = 's'
      ccurve(5,e2) = 's'
      ccurve(6,e ) = 's'
      ccurve(6,e1) = 's'
      ccurve(6,e2) = 's'

      call rzero(curve(1,1,e),3*72)
      curve(4,5,e ) = r1
      curve(4,5,e1) = r1
      curve(4,5,e2) = r1
      curve(4,6,e ) = r2
      curve(4,6,e1) = r2
      curve(4,6,e2) = r2

      letapt(e) = 'A'
      numapt(e) = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_cartesian_shell(r1,r2,e,ifsix) ! build 7 elements, use rotation

      include 'basics.inc'
      integer e,e1

      logical ifsix

      call sc_cartesian_shell_a(r1,r2,e,ifsix)
      e1 = e+3

      call sc_cartesian_shell_b(r1,r2,e1,ifsix)
      e1 = e1+3

      call sc_cartesian_shell_c(r1,r2,e1,ifsix)

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_cartesian_shell_a(r1,r2,e,ifsix) ! build 3 elements, use rotation

      include 'basics.inc'
      logical ifsix
      integer e,e1,e2
      common /xyzr/ xr(100),yr(100),zr(100)

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      p1 = r1
      p2 = r1/s2
      p3 = r1/s3

      x1 = r2
      x2 = r2/2   ! by fiat
      if (ifsix) x2 = 2.*r2/3   ! by fiat

      zr(1)=p1
      zr(2)=p2
      zr(3)=p3
      zr(4)=p2
      zr(5)=r2
      zr(6)=r2
      zr(7)=r2
      zr(8)=r2

      xr(2)=p2
      xr(3)=p3
      xr(6)=x2
      xr(7)=x2

      yr(3)=p3
      yr(4)=p2
      yr(7)=x2
      yr(8)=x2

      call el_convert  (e,xr,yr,zr)

      e1 = e+1
      e2 = e+2

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e1,xr,yr,zr)

      call permute_xyz    (xr,yr,zr,8)
      call el_convert     (e2,xr,yr,zr)

      ccurve(5,e ) = 's'
      ccurve(5,e1) = 's'
      ccurve(5,e2) = 's'

      call rzero(curve(1,1,e),3*72)
      curve(4,5,e ) = r1
      curve(4,5,e1) = r1
      curve(4,5,e2) = r1

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_cartesian_shell_b(r1,r2,e,ifsix) ! build 3 elements, use rotation

      include 'basics.inc'
      logical ifsix
      integer e,e1,e2
      common /xyzr/ xr(100),yr(100),zr(100)

c---------- the stuff below for midside node extraction -----------------
      parameter      (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation



      e2 = e-2  ! 2nd element in sc_cartesian_shell_a, below fiducial element
      ie = 3    ! edge number for elment "below" fiducial element

      call genxyz_e (xp,yp,zp,e2,nxm,nym,nzm) ! fill x27(.,e2)
      xm = x27(eindx(ie),e2)
      ym = y27(eindx(ie),e2)
      zm = z27(eindx(ie),e2)

c---------- the stuff above for midside node extraction -----------------

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      p1 = r1
      p2 = r1/s2
      p3 = r1/s3

      x1 = r2
      x2 = r2/2   ! by fiat
      if (ifsix) x2 = 2.*r2/3   ! by fiat

      xr(1)=p2    !             ^
      xr(2)=x1    !             |            |                          
      xr(3)=x1    !             |\           |                         
      xr(4)=p3    !             |  \         |                    
      xr(5)=x2    !             |    p3------+                    
      xr(6)=x1    !             |     \      |                      
      xr(7)=x1    !             |      \     |                       
      xr(8)=x2    !             O------------+---------------> x        
                  !                    p2    x1
      yr(3)=x2
      yr(4)=p3
      yr(7)=x2
      yr(8)=x2

      zr(1)=p2
      zr(2)=x2
      zr(3)=x2
      zr(4)=p3
      zr(5)=x1
      zr(6)=x1
      zr(7)=x1
      zr(8)=x1

      xr(9)=xm ! midside node, edge 1
      yr(9)=ym
      zr(9)=zm

c     do k=1,9
c        write(6,*) xr(k),yr(k),zr(k),'  MIDSIDE ?'
c     enddo
c     stop

      call rzero(curve(1,1,e),3*72)
      do e1=e,e+2
         call el_convert  (e1,xr,yr,zr)
         ccurve(4,e1) = 'm'
         curve (1,4,e1) = xr(9)
         curve (2,4,e1) = yr(9)
         curve (3,4,e1) = zr(9)
         call permute_xyz    (xr,yr,zr,9)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_cartesian_shell_c(r1,r2,e,ifsix) 

c     build one element, use rotation

      include 'basics.inc'
      logical ifsix

      integer e,e1
      common /xyzr/ xr(100),yr(100),zr(100)

      call rzero(xr,8)
      call rzero(yr,8)
      call rzero(zr,8)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      p1 = r1
      p2 = r1/s2
      p3 = r1/s3

      x1 = r2
      x2 = r2/2   ! by fiat
      if (ifsix) x2 = 2.*r2/3   ! by fiat

      xr(1)=p3
      xr(2)=x1
      xr(3)=x1
      xr(4)=x2
      xr(5)=x2
      xr(6)=x1
      xr(7)=x1
      xr(8)=x2

      yr(1)=p3
      yr(2)=x2
      yr(3)=x1
      yr(4)=x1
      yr(5)=x2
      yr(6)=x2
      yr(7)=x1
      yr(8)=x1

      zr(1)=p3
      zr(2)=x2
      zr(3)=x2
      zr(4)=x2
      zr(5)=x1
      zr(6)=x1
      zr(7)=x1
      zr(8)=x1

      call rzero(curve(1,1,e),72)
      do e1=e,e+2
         call el_convert     (e1,xr,yr,zr)
         call permute_xyz    (xr,yr,zr,8)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_pipe_section (radii,nr,e,ifsix)

      include 'basics.inc'
      logical ifsix

      integer e,e1,e2
      real radii(nr)

      common /xyzr/ xr(100),yr(100),zr(100)

      e1 = e

      call rzero(xr,8)
      call rzero(yr,8)

      z0 = -radii(nr+1)
      z1 = 0
      call cfill(zr(1),z0,4)
      call cfill(zr(5),z1,4)

      s2=sqrt(2.0)
      s3=sqrt(3.0)

      r1 = radii(1)
      r2 = radii(2)

      rr = sqrt(r1*r2) ! geometric mean
      rr = .5*(sqrt(r1*r2)+r1)
      x1 = r1
      x2 = rr/s2

      xr(2) = x1
      xr(3) = x2
      yr(3) = x2
      yr(4) = x1

      call el_convert_2to3(e1,xr,yr,zr)
      e1 = e1+1

      do i=2,nr-1  ! work from inner to outer radii: 1st & last are special

         r1 = radii(i-1)
         r2 = radii(i)
         rr = r1
         if (i.eq.2) rr=.5*(sqrt(r1*r2)+r1)

         p1 = r1
         p2 = rr/s2
         q1 = r2
         q2 = r2/s2

         xr(1) = p1
         xr(2) = q1
         xr(3) = q2
         xr(4) = p2

         yr(1) =  0
         yr(2) =  0
         yr(3) = q2
         yr(4) = p2

         call el_convert_2to3(e1,xr,yr,zr)

         if (i.gt.2) then
            ccurve(  4,e1) = 'C'
            curve (1,4,e1) = -r1
            ccurve(  8,e1) = 'C'
            curve (1,8,e1) = -r1
         endif
         ccurve(  2,e1) = 'C'
         curve (1,2,e1) = r2
         ccurve(  6,e1) = 'C'
         curve (1,6,e1) = r2

         call convert_top_to_mid(e1)

         e1 = e1+1

         xr(1) = p2
         xr(2) = q2
         xr(3) =  0
         xr(4) =  0

         yr(1) = p2
         yr(2) = q2
         yr(3) = q1
         yr(4) = p1

         call el_convert_2to3(e1,xr,yr,zr)

         if (i.gt.2) then
            ccurve(  4,e1) = 'C'
            curve (1,4,e1) = -r1
            ccurve(  8,e1) = 'C'
            curve (1,8,e1) = -r1
         endif
         ccurve(  2,e1) = 'C'
         curve (1,2,e1) = r2
         ccurve(  6,e1) = 'C'
         curve (1,6,e1) = r2

         call convert_top_to_mid(e1)

         e1 = e1+1

      enddo

      i  = nr
      r1 = radii(i-1)
      r2 = radii(i)
      x1 = r2
      x2 = r2/2   ! by fiat
      if (ifsix) x2 = 2.*r2/3   ! by fiat

      p1 = r1
      p2 = r1/s2

      xr(1) = r1
      xr(2) = r2
      xr(3) = r2
      xr(4) = p2
      yr(1) =  0
      yr(2) =  0
      yr(3) = x2
      yr(4) = p2
      call el_convert_2to3(e1,xr,yr,zr)
      ccurve(  4,e1) = 'C'
      curve (1,4,e1) = -r1
      ccurve(  8,e1) = 'C'
      curve (1,8,e1) = -r1
      call convert_top_to_mid(e1)
      e1 = e1+1

      xr(1) = p2
      xr(2) = x2
      xr(3) =  0
      xr(4) =  0
      yr(1) = p2
      yr(2) = r2
      yr(3) = r2
      yr(4) = r1
      call el_convert_2to3(e1,xr,yr,zr)
      ccurve(  4,e1) = 'C'
      curve (1,4,e1) = -r1
      ccurve(  8,e1) = 'C'
      curve (1,8,e1) = -r1
      call convert_top_to_mid(e1)
      e1 = e1+1

      xr(1) = p2
      xr(2) = x1
      xr(3) = x1
      xr(4) = x2
      yr(1) = p2
      yr(2) = x2
      yr(3) = x1
      yr(4) = x1
      call el_convert_2to3(e1,xr,yr,zr)
      e1 = e1+1

      return
      end
c-----------------------------------------------------------------------
      subroutine convert_top_to_mid(e)

      include 'basics.inc'
      integer e,e1,e2
      common /xyzr/ xr(100),yr(100),zr(100)

      parameter      (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      call genxyz_e (xp,yp,zp,e,nxm,nym,nzm) ! fill x27(.,e2)

      do ie=5,8
         if (ccurve(ie,e).eq.'C') then
            ccurve(  ie,e) = 'm'
            curve (1,ie,e) = x27(eindx(ie),e)
            curve (2,ie,e) = y27(eindx(ie),e)
            curve (3,ie,e) = z27(eindx(ie),e)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rep_rotate_mesh2
      include 'basics.inc'
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
      call do_rep_rotate_mesh(1,nel,angle_deg,nrep,axisr)

      return
      end
c-----------------------------------------------------------------------
      subroutine do_rep_rotate_meshg(e0,e1,angle_deg,nrep,normal)

c     Replicate & Rotate about general normal vector

      include 'basics.inc'
      integer e,e0,e1,e2,e3
      real normal(3)


      e3 = nel
      do i=1,nrep

         angle_deg_i = i*angle_deg

         e2 = e3+1
         e3 = e2 + (e1-e0)

         call copy_sub_mesh(e0,e1,e2)

         do e=e2,e3
            call rotate_el_vec(e,normal,angle_deg_i)
         enddo

      enddo
      nel = e3

      return
      end
c-----------------------------------------------------------------------
      subroutine do_rep_rotate_mesh(e0,e1,angle_deg,nrep,axisr)
      include 'basics.inc'
      integer e0,e1,e2,e3
      character*1 axisr


      e3 = nel
      do i=1,nrep

         e2 = e3+1
         e3 = e2 + (e1-e0)

         call copy_sub_mesh(e0,e1,e2)
         angle_deg_i = i*angle_deg
         call rotate_submesh_3d(e2,e3,angle_deg_i,axisr)

      enddo

      nel = e3

      return
      end
c-----------------------------------------------------------------------
      subroutine octquad_split_submesh(e0,e1,ifoct)
      include 'basics.inc'
      common /splitt/ enew(nelm),ind(nelm)
      dimension liste(8)
      integer e,e0,e1,en

      logical ifoct

      noct = 4
      if (ifoct.and.if3d) noct=8

      neln=nel + (1+e1-e0)*noct
      nelm1=nelm-2
      if (neln.gt.nelm1) then
         call prs('Number of elements after OctSplit operation$')
         write(s,51) neln,nelm1
         call prs(s)
   51    format(2x,'(',i9,') would be greater than the allowed maxium '
     $            ,'(',i9,').$')
         call prs('Aborting octquad_split.$')
         return
      endif

c     Renumber ALL elements, so that low numbered elements will remain
c     low numbered.  This will be achieved by assigning IE+0.1 to the
c     new element numbers, and then sorting the list.

      neln=nel
      call rint(enew,nel)

      do e=e0,e1
         en = 0
         do ioct = 2,noct
            neln=neln+1
c           enew(neln) = e + 0.1*ioct
            enew(neln) = e + 0.1
            en = en+1
            liste(en) = neln
         enddo
         liste(Noct) = Neln+1

         if (ifoct) then
            call octsplite(e,liste) 
         else
            call qsplite(e,liste) 
         endif
      enddo

      write(6,*) e0,e1,nel,noct,' e0e1'
c     do i=1,neln,8
c        write(6,8) (enew(j),j=i,i+7)
c  8     format(8f8.1,' enew')
c     enddo

c     Generate new sub-elements as a result of the oct-split action.

      nnew = neln-nel
      nel  = neln

C     Elements are all set. Sort.

      call sort   (enew,ind,nel)
      call swapel (ind,nel,nelm1)
      call curcnt
c     call vertadj


      write(s,300) nnew,nel
  300 format(i11,' new elements generated in octspl. NEL=',i11,'$')
      call prs(s)

      call gencen

      return
      end
c-----------------------------------------------------------------------
      subroutine mesh4to3(p9,idir)

      include 'basics.inc'
      real p9(3,9)
      integer ptr9(4,4),e
      save    ptr9
      data    ptr9 / 5,4,1,2 , 8,7,4,5 , 9,8,5,6 , 6,5,2,3 /

      if (idir.eq.1) then ! mesh to p9
         do e=1,4
         do i=1,4
            j=i+4
            k=ptr9(i,e)
            p9(1,k) = x(j,e)
            p9(2,k) = y(j,e)
            p9(3,k) = z(j,e)
            write(30,1) k,e,i,p9(1,k),p9(2,k),p9(3,k)
 1          format(3i5,3f12.5)
         enddo
         enddo
      else
         do e=1,4
         do i=1,4
            j=i+4
            k=ptr9(i,e)
            x(j,e) = p9(1,k)
            y(j,e) = p9(2,k)
            z(j,e) = p9(3,k)
         enddo
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine proj2plane(x,x0,p1,p2,p3,tol)
      real x(3),x0(3),p1(3),p2(3),p3(3)
      real nh(3),v0(3),v2(3),v3(3)


      call sub3(v2,p2,p1,3)
      call sub3(v3,p3,p1,3)
      call vcross_normal(nh,sine,v2,v3)

      call sub3(v0,x0,p1,3)

      alpha = dot(v0,nh,3)
      if (alpha.lt.0) call chsign(nh,3)

      d2plane=(x(1)-p1(1))*nh(1)+(x(2)-p1(2))*nh(2)+(x(3)-p1(3))*nh(3)

      if (d2plane.lt.tol) then
         x(1) = x(1) - d2plane*nh(1)
         x(2) = x(2) - d2plane*nh(2)
         x(3) = x(3) - d2plane*nh(3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine proj2line(x,p0,p1)
      real x(3),p0(3),p1(3),dx(3),t(3)

      call sub3  (t,p1,p0,3)
      call norm3d(t)                 ! Unit tangent vector

      call sub3(dx,x,p0,3)           ! Difference between x and base point

      alpha = dot(dx,t,3)            ! Inner product of dx and t

      do i=1,3
         x(i) = p0(i) + alpha*t(i)   ! Projection onto line [p0,p1]
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine build_corner_mesh(p9,nrefine) ! Here, we build the corner mesh

      include 'basics.inc'
      integer e,f,e0,e1,e2,e3

      real pp(3,3,2),pl(3,2),p0(3),p9(3,9)
      save pp,pl,p0
      data pp / 0.,0.,1. , 1.,0.,1. , 0.,1.,0.
     $        , 0.,0.,1. , 0.,1.,1. , 1.,0.,0. /
      data pl / 0.,0.,1. , 1.,1.,0. /
      integer pc(4)
      save    pc
      data    pc / 9,8,5,6 /


      e2=nel+1

      call copyel(nel,e2)

      do f=1,6
         ccurve(f,e2)=' '
      enddo

      call rone (x(1,e2),8)
      call rone (y(1,e2),8)
      call rzero(z(1,e2),8)

      do i=1,4
         k=pc(i)
         x(i,e2)=p9(1,k)
         y(i,e2)=p9(2,k)
         z(i,e2)=p9(3,k)
      enddo

      x(8,e2) = 1 - z(1,e2)  ! Exploit symmetries
      x(7,e2) = x(4,e2)

      y(6,e2) = 1 - z(1,e2)
      y(7,e2) = y(2,e2)

      nel = nel+1

      do k=1,nrefine
         call octquad_split_submesh(e2,nel,.true.) ! Oct refine base mesh
      enddo

c     Now copy and rotate sub-block

      e3 = nel
      call rone(p0,3)
      call do_rep_rotate_meshg(e2,e3,120.,1,p0)

      call translate_sub_mesh(e2,nel,-.5,-.5,-.5)
      call do_rep_rotate_mesh(e2,nel,180.,1,'Z')
      call translate_sub_mesh(e2,nel,0.5,0.5,0.5)

      return
      end
c-----------------------------------------------------------------------
      subroutine build_center_mesh(p9,nrefine) ! Build the center tet mesh

      include 'basics.inc'
      integer e,f,e0,e1,e2,e3,e4

      real p0(3),p9(3,9),p5(3,5),pr(3,5),a(3,3),normal(3)
      integer pc(4)
      save    pc
      data    pc / 5,4,1,2 /
      save    normal
      data    normal / -1. , -1. , 1. /

      call cfill(p5,0.5,15)
      do i=1,4
         k=pc(i)
         p5(1,i)=p9(1,k)-p9(1,5)  ! Shift defining tet face to origin
         p5(2,i)=p9(2,k)-p9(2,5)
         p5(3,i)=p9(3,k)-p9(3,5)
      enddo
      call sub2(p5(1,5),p9(1,5),3)
      
      call gen_rotate_mat_3d(a,normal, 120.)

      e0=nel+1

      call copyel(nel,e0)

      do f=1,6
         ccurve(f,e0)=' '
         cbc(f,e0,1) = 'V'
      enddo

      do i=1,4
         x(i,e0)=p5(1,i)
         y(i,e0)=p5(2,i)
         z(i,e0)=p5(3,i)
      enddo
      x(7,e0)=p5(1,5)
      y(7,e0)=p5(2,5)
      z(7,e0)=p5(3,5)

      call mxm(a,3,p5,3,pr,4)
      x(6,e0)=pr(1,3)
      y(6,e0)=pr(2,3)
      z(6,e0)=pr(3,3)
      call mxm(a,3,pr,3,p5,4)
      x(5,e0)=p5(1,4)
      y(5,e0)=p5(2,4)
      z(5,e0)=p5(3,4)
      x(8,e0)=p5(1,3)
      y(8,e0)=p5(2,3)
      z(8,e0)=p5(3,3)

      do i=1,8                          ! Translate back
         x(i,e0) = x(i,e0) + p9(1,5)
         y(i,e0) = y(i,e0) + p9(2,5)
         z(i,e0) = z(i,e0) + p9(3,5)
      enddo


      nel = nel+1
      do k=1,nrefine
         call octquad_split_submesh(e0,nel,.true.) ! Oct refine base mesh
      enddo

c     Now copy and rotate sub-block
      e1 = nel
      call rone(p0,3)
      call do_rep_rotate_meshg(e0,e1,120.,2,p0)

      call translate_sub_mesh(1,nel,-.5,-.5,-.5)
      p0(2) = -1
      call do_rep_rotate_meshg(e0,e1,120.,1,p0)
      call translate_sub_mesh(1,nel,0.5,0.5,0.5)

      return
      end
c-----------------------------------------------------------------------
      subroutine skin_fcc_spheres(rt,r0,nrefine) 

      include 'basics.inc'
      integer e,f,e0,e1,e2,e3,e4
      real p0(3)

      nel0 = nel ! Track baseline number of elements
      e0   = nel0+1

      call sc_std_shell(rt,r0,e0)     ! Standard spherical shell
      nel = nel+3                     ! Keep only the first elements

      e1  = nel

      write(6,*) e0,e1,nrefine,' in skin',nel,rt,r0

c     call prexit  ! Definitely messed up

      do e=e0,e1
         do f=1,4
            cbc(f,e,1) = 'SYM'
            cbc(f,e,1) = 'O  '
            cbc(f,e,2) = 'I  '
         enddo
         cbc(5,e,1) = 'W  '
         cbc(5,e,2) = 't  '
         cbc(6,e,1) = 'v  '
         cbc(6,e,2) = 't  '
      enddo
      e = e0                  ! Identify first element
      cbc(6,e,1) = 'O  '
      cbc(6,e,2) = 'O  '

      call rotate_submesh_3d (e0,e1,180.,'X')
      call rotate_submesh_3d (e0,e1,-90.,'Z')
      call translate_sub_mesh(e0,e1,1.0,1.0,1.0)

      call octquad_split_submesh(e0,e1,.false.) ! Quad refine base mesh

      call translate_sub_mesh (e0,nel,-.5,-.5,-.5)
      call do_rep_rotate_mesh (e0,nel,180.,1,'X')
      call do_rep_rotate_mesh (e0,nel,180.,1,'Z')
      call translate_sub_mesh (e0,nel,0.5,0.5,0.5)

      if (nrefine.ge.1)
     $    call octquad_split_submesh(e0,nel,.false.) ! Quad refine base mesh
      if (nrefine.ge.2)
     $    call octquad_split_submesh(e0,nel,.true.)  ! Oct refine base mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine fcc_x8 ! copy unit 1/2-cell to 4x-cell

      include 'basics.inc'

      call do_rep_rotate_mesh (1,nel, 90.,3,'X')
      call do_rep_rotate_mesh (1,nel,180.,1,'Z')

      return
      end
c-----------------------------------------------------------------------
      subroutine fcc_base(p9,rt,r0,r1,nrefine)  ! Build fcc mesh
      include 'basics.inc'

      real p9(3,9)

      integer e,f,e0,e1,e2

      real pp(3,3,2),pl(3,2),p0(3)
      save pp,pl,p0
      data pp / 0.,0.,1. , 1.,0.,1. , 0.,1.,0.
     $        , 0.,0.,1. , 0.,1.,1. , 1.,0.,0. /
      data pl / 0.,0.,1. , 1.,1.,0. /


c     NOTES:
c
c        Target R = .5*6 cm / 4.59619 = .652714
c

      nel0 = nel ! Track baseline number of elements
      e0   = nel0+1

      call sc_std_shell(r0,r1,1)     ! Standard spherical shell
      nel = 1                        ! Keep only the first elements

      e1  = nel

      do e=e0,e1
         do f=1,4
            cbc(f,e,1) = 'SYM'
            cbc(f,e,1) = 'O  '
            cbc(f,e,2) = 'I  '
         enddo
         cbc(5,e,1) = 'W  '
         cbc(5,e,2) = 't  '
         cbc(6,e,1) = 'v  '
         cbc(6,e,2) = 't  '
      enddo
      e = e0                  ! Identify first element
      cbc(6,e,1) = 'O  '
      cbc(6,e,2) = 'O  '

      call rotate_submesh_3d (e0,e1,180.,'X')
      call rotate_submesh_3d (e0,e1,-90.,'Z')
      call translate_sub_mesh(e0,e1,1.0,1.0,1.0)

      call octquad_split_submesh(e0,e1,.false.) ! Quad refine base mesh

      do e=1,4
         ccurve(6,e)=' '
      enddo

      call rone(p0,3)
      call mesh4to3(p9,1) ! map 4 elements to 3 x 3 array
      call proj2plane(p9(1,2),p0,pp(1,1,2),pp(1,2,2),pp(1,3,2),1.e20)
      call proj2plane(p9(1,3),p0,pp(1,1,2),pp(1,2,2),pp(1,3,2),1.e20)
      call proj2plane(p9(1,6),p0,pp(1,1,2),pp(1,2,2),pp(1,3,2),1.e20)
      call proj2plane(p9(1,4),p0,pp(1,1,1),pp(1,2,1),pp(1,3,1),1.e20)
      call proj2plane(p9(1,7),p0,pp(1,1,1),pp(1,2,1),pp(1,3,1),1.e20)
      call proj2plane(p9(1,8),p0,pp(1,1,1),pp(1,2,1),pp(1,3,1),1.e20)
      call proj2line (p9(1,5),pl(1,1),pl(1,2))
      call mesh4to3(p9,2) ! map 3x3 array to 4 elements

      call do_rep_rotate_meshg(1,nel,120.,2,p0)
      call translate_sub_mesh (1,nel,-.5,-.5,-.5)
      call do_rep_rotate_mesh (1,nel,180.,1,'X')
      call do_rep_rotate_mesh (1,nel,180.,1,'Z')
      call translate_sub_mesh (1,nel,0.5,0.5,0.5)


      if (nrefine.ge.1)
     $   call octquad_split_submesh(1,nel,.false.) ! Quad refine base mesh

      do k=2,nrefine
         call octquad_split_submesh(1,nel,.true.)  ! Oct refine base mesh
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fcc2  ! Build fcc mesh

      include 'basics.inc'
      integer e,f,e0,e1,e2

      real pp(3,3,2),pl(3,2),p0(3)
      save pp,pl,p0
      data pp / 0.,0.,1. , 1.,0.,1. , 0.,1.,0.
     $        , 0.,0.,1. , 0.,1.,1. , 1.,0.,0. /
      data pl / 0.,0.,1. , 1.,1.,0. /

      common /fccpts/ p9(3,9)



c     NOTES:  Target R = .5*6 cm / 4.59619 = .652714

      rt = 0.652714  ! Target radius
      r0 = 0.678
      r1 = 0.735

      nrefine = 3  ! NEL = 143360   ( untested )
      nrefine = 2  ! NEL =  28000   ( n ~ 10 million for N=8 )

      call fcc_base           (p9,rt,r0,r1,nrefine)  ! Build fcc mesh
      call build_corner_mesh           (p9,nrefine)
      call build_center_mesh           (p9,nrefine)  ! Remaining tet-mesh
      call skin_fcc_spheres         (rt,r0,nrefine) 
      call fcc_x8                                    ! copy to 4x-cell

      call fcc_set_bcs(rt)

      return
      end
c-----------------------------------------------------------------------
      subroutine fcc_set_bcs(r0)  ! BCs for fcc mesh
      include 'basics.inc'

      integer e,f

      tol = 1.e-4

      do e=1,nel
      do f=1,6
         cbc(f,e,1)='   '
         cbc(f,e,2)='   '

         dr=abs(curve(4,f,e)-r0)

         if (ccurve(f,e).eq.'s'.and.dr.lt.tol) then
            cbc(f,e,1) = 'W  '
            cbc(f,e,2) = 'f  '
         endif

      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_set_bc_etc ! Extraneous spherical cap (SC) cleanup
      include 'basics.inc'
      integer e,f

      do e=1,nel
      do f=1,6
         cbc(f,e,1) = 'v  '
         cbc(f,e,2) = 't  '
      enddo
      enddo

      do e=1,nel
         numapt(e) = 1
         letapt(e) = 'A'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_el_to_safe_haven(e_safe) ! Cache element list
      include 'basics.inc'

      integer e,e0,e1,e_safe

      common /csafe/ ne_safe

      e_safe = (nelm-3) - ne_safe - nel ! This keeps nelm open

      ne_safe = ne_safe + nel       ! Track number being saved

      call copy_sub_mesh(1,nel,e_safe) ! Copy current element list

      return
      end
c-----------------------------------------------------------------------
      subroutine free_safe_haven(e_safe) ! Free (ALL) cached space
      include 'basics.inc'

      integer e,e0,e1,e_safe

      common /csafe/ ne_safe

      ne_safe = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_make_sphere_cap(radii,nr)
      include 'basics.inc'

      real radii(nr)

      integer e,e0,e1,f
      logical ifsix

      common /xyzr/ xr(100),yr(100),zr(100)

      parameter      (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3) ! for genxyze_e

      ifsix = .false.
      if (nr.eq.7) ifsix = .true.

      e0 = nel
      e  = e0 + 1

      call sc_central_shell    (radii(1),radii(2),e)
      e = e+1

      call sc_transition_shell (radii(1),radii(2),e)
      e = e+3

      do i=3,nr-1  ! work from inner to outer radii: 1st & last are special
         call sc_std_shell     (radii(i-1),radii(i),e)
         e = e+3
      enddo
      do e1=e0+1,e
      do f=1,6
         cbc(f,e1,1) = 'v  '
         cbc(f,e1,2) = 't  '
      enddo
      enddo


      call sc_cartesian_shell  (radii(nr-1),radii(nr),e,ifsix)
      do e1=e+1,e+7
      do f=1,6
         cbc(f,e1,1) = 'W  '
         cbc(f,e1,2) = 'I  '
      enddo
      enddo

      e = e+7

      call prs('WARNING: Geometry inconsistent!!$')
      call prs('You must call fix_geom from usrdat2!!$')
c     do i=nel+1,e-1
c        call genxyz_e (xp,yp,zp,i,nxm,nym,nzm) ! fill x27(.,e2)
c        call fix_m_curve(i) ! convert cap to midside nodes
c     enddo

      ne_cap = e-1  ! number of elements in spherical cap

      call sc_pipe_section  (radii,nr,e,ifsix) ! Add pipe
      e   = e+6 + 2*(nr-3)
      nel = e-1
      call sc_set_bc_etc

      if (ifsix) then
         call sc_6x6_modification(radii,nr,ne_cap) ! 3x3 faces
         call sc_set_bc_etc
c        call prexit
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sc_6x6_modification(radii,nr,ne_cap) ! 3x3 faces
      include 'basics.inc'

      common /cisplit/ isplit(nelm)

      real radii(nr)

      integer e,e0,e1,f,e0s
      integer e_pipe_0,e_pipe_1

      write(6,*) nr,(radii(k),k=1,nr),' Radii 6x6'

      ne_orig = nel        ! Number of elements before refinement
      ne_pipe = nel-ne_cap ! Number of elements in pipe section

      call copy_el_to_safe_haven(e0s)         ! Cache existing element list


!- - - Refine pipe section - - - - - - - - - - - - - - - - - - - - - - -

      e_pipe_0 = e0s + ne_cap                 ! Pointer to pipe section
      e_pipe_1 = e_pipe_0 + (ne_pipe-1)
      call copy_sub_mesh(e_pipe_0,e_pipe_1,1) ! Retrieve pipe
      
      nel = ne_pipe

c     call midside_convert_all ! Convert all elements to midside?

      call mark (1,1,0.5,isplit)  ! Zipper pipe element 1, side 1
      call split(0.5,isplit)

      call mark (1,2,0.5,isplit)  ! Zipper pipe element 1, side 2
      call split(0.5,isplit)

      return
      end
c-----------------------------------------------------------------------
