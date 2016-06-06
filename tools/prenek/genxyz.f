      SUBROUTINE GENXYZ (XML,YML,ZML,IE,NXL,NYL,NZL)
C
      INCLUDE 'basics.inc'
C
C     Note : CTMP1 is used in this format in several subsequent routines
C
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMP1/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      COMMON /PFFLG/  ILGRNG,ISCND,IENTRP,INEWTX
C
      DIMENSION XML(NXL,NYL,NZL,1),YML(NXL,NYL,NZL,1),ZML(NXL,NYL,NZL,1)
      DIMENSION XCB(2,2,2),YCB(2,2,2),ZCB(2,2,2)
C
      CHARACTER*1 CCV
C
      DIMENSION INDX(8)
      SAVE      INDX
      DATA      INDX  / 1, 2, 4, 3, 5, 6, 8, 7 /
C
C     Initialize geometry arrays
C
      ILGRNG=0
      NXYZ = NXL*NYL*NZL
      CALL RZERO(XML(1,1,1,IE),NXYZ)
      CALL RZERO(YML(1,1,1,IE),NXYZ)
      CALL RZERO(ZML(1,1,1,IE),NXYZ)
C
C   Preprocessor Corner notation:      Symmetric Corner notation:
C
C           4+-----+3    ^ s                    3+-----+4    ^ s
C           /     /|     |                      /     /|     |
C          /     / |     |                     /     / |     |
C        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
C         |     | /     /                     |     | /     /
C         |     |/     /                      |     |/     /
C        5+-----+6    t                      5+-----+6    t
C
      INDX(1)=1
      INDX(2)=2
      INDX(3)=4
      INDX(4)=3
      INDX(5)=5
      INDX(6)=6
      INDX(7)=8
      INDX(8)=7
      NDIM2 = 2**NDIM
C
      CALL LEGEND(ZGML,WORK,NX)
C
      DO 10 IX=1,NXL
         H(IX,1,1)=(1.0-ZGML(IX,1))*0.5
         H(IX,1,2)=(1.0+ZGML(IX,1))*0.5
   10 CONTINUE
      DO 20 IY=1,NYL
         H(IY,2,1)=(1.0-ZGML(IY,1))*0.5
         H(IY,2,2)=(1.0+ZGML(IY,1))*0.5
   20 CONTINUE
      IF (IF3D) THEN
         DO 30 IZ=1,NZL
            H(IZ,3,1)=(1.0-ZGML(IZ,1))*0.5
            H(IZ,3,2)=(1.0+ZGML(IZ,1))*0.5
   30    CONTINUE
      ELSE
         CALL RONE(H(1,3,1),NZL)
         CALL RONE(H(1,3,2),NZL)
      ENDIF
C
      do 50 ix=1,ndim2
         i=indx(ix)
         xcb(ix,1,1)=x(i,ie)
         ycb(ix,1,1)=y(i,ie)
         zcb(ix,1,1)=z(i,ie)
   50 continue
C
C     Map R-S-T space into physical X-Y-Z space.
C
      IZTMAX = NDIM-1
      DO 200 IZT=1,IZTMAX
      DO 200 IYT=1,2
      DO 200 IXT=1,2
C
      DO 200 IZ=1,NZL
      DO 200 IY=1,NYL
         HH = H(IY,2,IYT)*H(IZ,3,IZT)
         DO 100 IX=1,NXL
            HHH = H(IX,1,IXT)*HH
            XML(IX,IY,IZ,IE)=XML(IX,IY,IZ,IE)+HHH*XCB(IXT,IYT,IZT)
            YML(IX,IY,IZ,IE)=YML(IX,IY,IZ,IE)+HHH*YCB(IXT,IYT,IZT)
            ZML(IX,IY,IZ,IE)=ZML(IX,IY,IZ,IE)+HHH*ZCB(IXT,IYT,IZT)
  100    CONTINUE
  200 CONTINUE
C
C     Deform surfaces - general 3D deformations
C                     - extruded geometry deformations
C
      NFACES = 2*NDIM
      DO 1000 IFACE=1,NFACES
        CCV = CCURVE(IFACE,IE)
        IF (CCV.EQ.'s') 
     $     CALL SPHSRF(XML,YML,ZML,IFACE,IE,NXL,NYL,NZL,WORK) 
        IF (CCV.EQ.'e') 
     $     CALL gensrf(XML,YML,ZML,IFACE,IE,NXL,NYL,NZL,zgml)
 1000 CONTINUE
C
c     DO 2000 ISID=1,8
c       CCV = CCURVE(ISID,IE)
c       IF (CCV.EQ.'C') 
c    $  CALL ARCSRF(XML,YML,ZML,NXL,NYL,NZL,IE,ISID)
C        IF (CCV.EQ.'W')
c2000 CONTINUE
C
      RETURN
      END
      SUBROUTINE SPHSRF(XML,YML,ZML,IFCE,IE,MX,MY,MZ,XYSRF) 
C
C     5 Aug 1988 19:29:52 
C
C     Program to generate spherical shell elements for NEKTON
C     input.  Paul F. Fischer
C
      INCLUDE 'basics.inc'
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      DIMENSION XML(MX,MY,MZ,1),YML(MX,MY,MZ,1),ZML(MX,MY,MZ,1)
      DIMENSION XYSRF(3,MX,MZ)
C
      COMMON /CTMP1/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      COMMON /CTMP0/ XCV(3,2,2),VN1(3),VN2(3)
     $              ,X1(3),X2(3),X3(3),DX(3),XCC(8),YCC(8),ZCC(8)
      DIMENSION IOPP(3),MXX(3)
C
C     Stuff which is normally in genxyz
C
      DO 11 IX=1,MX
         H(IX,1,1)=(1.0-ZGML(IX,1))*0.5
         H(IX,1,2)=(1.0+ZGML(IX,1))*0.5
   11 CONTINUE
      DO 21 IY=1,MY
         H(IY,2,1)=(1.0-ZGML(IY,1))*0.5
         H(IY,2,2)=(1.0+ZGML(IY,1))*0.5
   21 CONTINUE
      IF (IF3D) THEN
         DO 31 IZ=1,MZ
            H(IZ,3,1)=(1.0-ZGML(IZ,1))*0.5
            H(IZ,3,2)=(1.0+ZGML(IZ,1))*0.5
   31    CONTINUE
      ELSE
         CALL RONE(H(1,3,1),MZ)
         CALL RONE(H(1,3,2),MZ)
      ENDIF
C
C     Determine geometric parameters
C
      MXM1 = MX-1
      MYM1 = MY-1
      MXY  = MX*MZ
      MXY3 = 3*MX*MZ
      XCTR   = CURVE(1,IFCE,IE)
      YCTR   = CURVE(2,IFCE,IE)
      ZCTR   = CURVE(3,IFCE,IE)
      RADIUS = CURVE(4,IFCE,IE)
      IFACE  = EFACE1(IFCE)
C
C     Generate (normalized) corner vectors XCV(1,i,j):
C
      do 10 i=1,8
         xcc(i)=x(i,ie)
         ycc(i)=y(i,ie)
         zcc(i)=z(i,ie)
   10 continue
      CALL CRN3D(XCV,XCC,YCC,ZCC,CURVE(1,IFCE,IE),IFACE)
C
C     Generate edge vectors on the sphere RR=1.0, 
C     for (r,s) = (-1,*),(1,*),(*,-1),(*,1)      
C
      CALL EDG3D(XYSRF,XCV(1,1,1),XCV(1,1,2), 1, 1, 1,MY,MX,MY)
      CALL EDG3D(XYSRF,XCV(1,2,1),XCV(1,2,2),MX,MX, 1,MY,MX,MY)
      CALL EDG3D(XYSRF,XCV(1,1,1),XCV(1,2,1), 1,MX, 1, 1,MX,MY)
      CALL EDG3D(XYSRF,XCV(1,1,2),XCV(1,2,2), 1,MX,MY,MY,MX,MY)
C
C     Generate intersection vectors for (i,j)
C
      DO 200 J=2,MYM1
         CALL CROSS(VN1,XYSRF(1,1,J),XYSRF(1,MX,J))
         DO 200 I=2,MXM1
            CALL CROSS(VN2,XYSRF(1,I,1),XYSRF(1,I,MY))
            IF (IFACE.EQ.1.OR.IFACE.EQ.4.OR.IFACE.EQ.5) THEN
               CALL CROSS(XYSRF(1,I,J),VN2,VN1)
            ELSE
               CALL CROSS(XYSRF(1,I,J),VN1,VN2)
            ENDIF
  200 CONTINUE
C
C     Normalize all vectors to the unit sphere.
C
      DO 300 I=1,MXY
         CALL NORM3D(XYSRF(1,I,1))
  300 CONTINUE
C
C     Scale by actual radius
C
      CALL CMULT(XYSRF,RADIUS,MXY3)
C
C     Add back the sphere center offset
C
      DO 400 I=1,MXY
         XYSRF(1,I,1)=XYSRF(1,I,1)+XCTR
         XYSRF(2,I,1)=XYSRF(2,I,1)+YCTR
         XYSRF(3,I,1)=XYSRF(3,I,1)+ZCTR
  400 CONTINUE
C
C
C     Transpose data, if necessary
C
      IF (IFACE.EQ.1.OR.IFACE.EQ.4.OR.IFACE.EQ.5) THEN
         DO 500 J=1  ,MY
         DO 500 I=J+1,MX
            TMP=XYSRF(1,I,J)
            XYSRF(1,I,J)=XYSRF(1,J,I)
            XYSRF(1,J,I)=TMP
            TMP=XYSRF(2,I,J)
            XYSRF(2,I,J)=XYSRF(2,J,I)
            XYSRF(2,J,I)=TMP
            TMP=XYSRF(3,I,J)
            XYSRF(3,I,J)=XYSRF(3,J,I)
            XYSRF(3,J,I)=TMP
  500    CONTINUE
      ENDIF
C
C     Compute surface deflection and perturbation due to face IFACE
C
      CALL DSSET(MX,MY,MZ)
      JS1    = SKPDAT(IFACE,1)
      JF1    = SKPDAT(IFACE,2)
      JSKIP1 = SKPDAT(IFACE,3)
      JS2    = SKPDAT(IFACE,4)
      JF2    = SKPDAT(IFACE,5)
      JSKIP2 = SKPDAT(IFACE,6)
C
      IOPP(1) = MX-1
      IOPP(2) = MX*(MY-1)
      IOPP(3) = MX*MY*(MZ-1)
      MXX(1)  = MX
      MXX(2)  = MY
      MXX(3)  = MZ
      IDIR    = 2*MOD(IFACE,2) - 1
      IFC2    = (IFACE+1)/2
      I=0
      DO 700 J2=JS2,JF2,JSKIP2
      DO 700 J1=JS1,JF1,JSKIP1
         I=I+1
         JOPP = J1 + IOPP(IFC2)*IDIR
         X2(1) = XML(J1,J2,1,IE)
         X2(2) = YML(J1,J2,1,IE)
         X2(3) = ZML(J1,J2,1,IE)
C
         DX(1) = XYSRF(1,I,1)-X2(1)
         DX(2) = XYSRF(2,I,1)-X2(2)
         DX(3) = XYSRF(3,I,1)-X2(3)
C
         MXS = MXX(IFC2)
         JOFF = (J1-JOPP)/(MXS-1)
         DO 600 IX = 2,MXS
            J = JOPP + JOFF*(IX-1)
            ZETA = 0.5*(ZGML(IX,IFC2) + 1.0)
            ZETA = 0.5*(ZGML(IX,1) + 1.0)
            XML(J,J2,1,IE) = XML(J,J2,1,IE)+DX(1)*ZETA
            YML(J,J2,1,IE) = YML(J,J2,1,IE)+DX(2)*ZETA
            ZML(J,J2,1,IE) = ZML(J,J2,1,IE)+DX(3)*ZETA
  600    CONTINUE
  700 CONTINUE
C
      RETURN
      END
      SUBROUTINE EDG3D(XYSRF,X1,X2,I1,I2,J1,J2,MX,MY)
C
C     Generate XYZ vector along an edge of a surface.
C
      INCLUDE 'basics.inc'
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMP1/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
C
      DIMENSION XYSRF(3,MX,MY)
      DIMENSION X1(3),X2(3)
      REAL U1(3),U2(3),B(3)
C
C     Normalize incoming vectors
C
      CALL COPY (U1,X1,3)
      CALL COPY (U2,X2,3)
      CALL NORM3D (U1)
      CALL NORM3D (U2)
C
C     Find normal to the plane and tangent to the curve.
C
      CALL CROSS(VN,X1,X2)
      CALL CROSS( B,VN,X1)
      CALL NORM3D (VN)
      CALL NORM3D (B)
C
      CTHETA = DOT(U1,U2,3)
      THETA  = ACOS(CTHETA)
C
      IJ = 0
      DO 200 J=J1,J2
      DO 200 I=I1,I2
         IJ = IJ + 1
         THETAP = 0.5*THETA*(ZGML(IJ,1)+1.0)
         CTP = COS(THETAP)
         STP = SIN(THETAP)
         DO 200 IV = 1,3
            XYSRF(IV,I,J) = CTP*U1(IV) + STP*B(IV)
  200 CONTINUE
      RETURN
      END
      REAL FUNCTION DOT(V1,V2,N)
C
C     Compute Cartesian vector dot product.
C
      DIMENSION V1(N),V2(N)
C
      SUM = 0
      DO 100 I=1,N
         SUM = SUM + V1(I)*V2(I)
  100 CONTINUE
      DOT = SUM
      RETURN
      END
      SUBROUTINE CRN3D(XCV,XC,YC,ZC,CURVE,IFACE)
      DIMENSION XCV(3,2,2),XC(8),YC(8),ZC(8),CURVE(4)
      DIMENSION INDX(8)
      SAVE      INDX
      DATA      INDX  / 1, 2, 4, 3, 5, 6, 8, 7 /
      DIMENSION INDVTX(4,6)
      SAVE      INDVTX
      DATA      INDVTX  / 1,5,3,7 , 2,4,6,8 , 1,2,5,6  
     $                  , 3,7,4,8 , 1,3,2,4 , 5,6,7,8 /
C
      EPS    = 1.0E-5
      XCTR   = CURVE(1)
      YCTR   = CURVE(2)
      ZCTR   = CURVE(3)
      RADIUS = CURVE(4)
C
      DO 10 I=1,4
         J=INDVTX(I,IFACE)
         K=INDX(J)
         XCV(1,I,1)=XC(K)-XCTR
         XCV(2,I,1)=YC(K)-YCTR
         XCV(3,I,1)=ZC(K)-ZCTR
   10 CONTINUE
C
C     Check to ensure that these points are indeed on the sphere.
C
      IF (RADIUS.LE.0.0) THEN
         WRITE(6,20) XCTR,YCTR,ZCTR,IFACE
  20     FORMAT(5X,'ERROR: Sphere of radius zero requested.'
     $       ,/,5X,'EXITING in CRN3Da',3E12.4,I3)
         CALL EXIT
      ELSE
         DO 40 I=1,4
            RADT=XCV(1,I,1)**2+XCV(2,I,1)**2+XCV(3,I,1)**2
            RADT=SQRT(RADT)
            TEST=ABS(RADT-RADIUS)/RADIUS
            IF (TEST.GT.EPS) THEN
             WRITE(6,30) 
     $       RADT,RADIUS,XCV(1,I,1),XCV(2,I,1),XCV(3,I,1)
   30        FORMAT(5X,'ERROR: Element vertex not on requested sphere.'
     $           ,/,5X,'EXITING in CRN3Db',5E12.4)
             CALL EXIT
            ENDIF
   40    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE DSSET(MX,MY,MZ)
C
C     Set up arrays IXCN,ESKIP,SKPDAT,NEDG,NOFFST for new MX,MY,MZ
C
      INCLUDE 'basics.inc'
      INTEGER MXO,MYO,MZO
      SAVE    MXO,MYO,MZO
      DATA    MXO,MYO,MZO /3*0/
C
C     Check if element surface counters are already set from last call...
C
      IF (MXO.EQ.MX.AND.MYO.EQ.MY.AND.MZO.EQ.MZ) RETURN
C
C     else, proceed....
C
      MXO = MX
      MYO = MY
      MZO = MZ
C
C     Establish corner to elemental node number mappings
C
C     Assign indices for direct stiffness summation of arbitrary faces.
C
C
C     Y-Z Planes (Faces 1 and 2)
C
      SKPDAT(1,1)=1
      SKPDAT(1,2)=MX*(MY-1)+1
      SKPDAT(1,3)=MX
      SKPDAT(1,4)=1
      SKPDAT(1,5)=MY*(MZ-1)+1
      SKPDAT(1,6)=MY
C
      SKPDAT(2,1)=1             + (MX-1)
      SKPDAT(2,2)=MX*(MY-1)+1   + (MX-1)
      SKPDAT(2,3)=MX
      SKPDAT(2,4)=1
      SKPDAT(2,5)=MY*(MZ-1)+1
      SKPDAT(2,6)=MY
C
C     X-Z Planes (Faces 3 and 4)
C
      SKPDAT(3,1)=1
      SKPDAT(3,2)=MX
      SKPDAT(3,3)=1
      SKPDAT(3,4)=1
      SKPDAT(3,5)=MY*(MZ-1)+1
      SKPDAT(3,6)=MY
C
      SKPDAT(4,1)=1           + MX*(MY-1)
      SKPDAT(4,2)=MX          + MX*(MY-1)
      SKPDAT(4,3)=1
      SKPDAT(4,4)=1
      SKPDAT(4,5)=MY*(MZ-1)+1
      SKPDAT(4,6)=MY
C
C     X-Y Planes (Faces 5 and 6)
C
      SKPDAT(5,1)=1
      SKPDAT(5,2)=MX
      SKPDAT(5,3)=1
      SKPDAT(5,4)=1
      SKPDAT(5,5)=MY
      SKPDAT(5,6)=1
C
      SKPDAT(6,1)=1           + MX*MY*(MZ-1)
      SKPDAT(6,2)=MX          + MX*MY*(MZ-1)
      SKPDAT(6,3)=1
      SKPDAT(6,4)=1
      SKPDAT(6,5)=MY
      SKPDAT(6,6)=1
C
      RETURN
      END
      SUBROUTINE RONE(X,N)
      DIMENSION X(1)
      DO 10 I=1,N
         X(I)=1.0
   10 CONTINUE
      RETURN
      END
      subroutine gensrf(XML,YML,ZML,IFCE,IE,MX,MY,MZ,zgml)
C
C     9 Mar 1994 
C
C     Program to generate surface deformations for NEKTON
C     input.  Paul F. Fischer
C
      INCLUDE 'basics.inc'
C
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      DIMENSION XML(MX,MY,MZ,1),YML(MX,MY,MZ,1),ZML(MX,MY,MZ,1)
     $             ,ZGML(LX1,3)
C
      COMMON /CTMP0/ dummy0
      COMMON /CTMP1/ dummy1
      real IOPP(3),MXX(3),X0(3),DX(3)
C
C
C     Algorithm:  .Project original point onto surface S
C                 .Apply Gordon Hall to vector of points between x_s and
C                  opposite face
C
C
      CALL DSSET(MX,MY,MZ)
C
      IFACE  = EFACE1(IFCE)
      JS1    = SKPDAT(IFACE,1)
      JF1    = SKPDAT(IFACE,2)
      JSKIP1 = SKPDAT(IFACE,3)
      JS2    = SKPDAT(IFACE,4)
      JF2    = SKPDAT(IFACE,5)
      JSKIP2 = SKPDAT(IFACE,6)
C
      IOPP(1) = MX-1
      IOPP(2) = MX*(MY-1)
      IOPP(3) = MX*MY*(MZ-1)
      MXX(1)  = MX
      MXX(2)  = MY
      MXX(3)  = MZ
      IDIR    = 2*MOD(IFACE,2) - 1
      IFC2    = (IFACE+1)/2
      I=0
C
C     Find a characteristic length scale for initializing secant method
C
      x0(1) = xml(js1,js2,1,ie)
      x0(2) = yml(js1,js2,1,ie)
      x0(3) = zml(js1,js2,1,ie)
      rmin  = 1.0e16
c
      DO 100 J2=JS2,JF2,JSKIP2
      DO 100 J1=JS1,JF1,JSKIP1
         if (j1.ne.js1.or.j2.ne.js2) then
            r2 = (x0(1) - xml(j1,j2,1,ie))**2 
     $         + (x0(2) - yml(j1,j2,1,ie))**2 
     $         + (x0(3) - zml(j1,j2,1,ie))**2 
            rmin = min(r2,rmin)
         endif
  100 CONTINUE
      dxc = 0.05*sqrt(rmin)
C
C     Project each point on this surface onto curved surface
C
      DO 300 J2=JS2,JF2,JSKIP2
      DO 300 J1=JS1,JF1,JSKIP1
         I=I+1
         JOPP = J1 + IOPP(IFC2)*IDIR
         X0(1) = XML(J1,J2,1,IE)
         X0(2) = YML(J1,J2,1,IE)
         X0(3) = ZML(J1,J2,1,IE)
C
         call prjects(x0,dxc,curve(1,ifce,ie),ccurve(ifce,ie))
         DX(1) = X0(1)-xml(j1,j2,1,ie)
         DX(2) = X0(2)-yml(j1,j2,1,ie)
         DX(3) = X0(3)-zml(j1,j2,1,ie)
         MXS = MXX(IFC2)
         JOFF = (J1-JOPP)/(MXS-1)
         DO 200 IX = 2,MXS
            J = JOPP + JOFF*(IX-1)
            ZETA = 0.5*(ZGML(IX,1) + 1.0)
            XML(J,J2,1,IE) = XML(J,J2,1,IE)+DX(1)*ZETA
            YML(J,J2,1,IE) = YML(J,J2,1,IE)+DX(2)*ZETA
            ZML(J,J2,1,IE) = ZML(J,J2,1,IE)+DX(3)*ZETA
  200    CONTINUE
  300 CONTINUE
C
      RETURN
      END
      subroutine prjects(x0,dxc,c,cc)
c
c     Project the point x0 onto surface described by characteristics
c     given in the array c and cc.  
c
c     dxc - characteristic length scale used to estimate gradient.
c
      real x0(3)
      real c(5)
      character*1 cc
      real x1(3)
      logical if3d
c
      if3d = .true.
      if (dxc.le.0.0) then
         write(6,*) 'invalid dxc',dxc,x0
         write(6,*) 'Abandoning prjects'
         return
      endif
c
      call copy(x1,x0,3)
      R0 = ressrf(x0,c,cc)
      if (r0.eq.0) return
c
c     Must at least use ctr differencing to capture symmetry!
c
      x1(1) = x0(1) - dxc
      R1 = ressrf(x1,c,cc)
      x1(1) = x0(1) + dxc
      R2 = ressrf(x1,c,cc)
      x1(1) = x0(1)
      Rx = 0.5*(R2-R1)/dxc
c
      x1(2) = x0(2) - dxc
      R1 = ressrf(x1,c,cc)/dxc
      x1(2) = x0(2) + dxc
      R2 = ressrf(x1,c,cc)/dxc
      x1(2) = x0(2)
      Ry = 0.5*(R2-R1)/dxc
c
      if (if3d) then
         x1(3) = x0(3) - dxc
         R1 = ressrf(x1,c,cc)/dxc
         x1(3) = x0(3) + dxc
         R2 = ressrf(x1,c,cc)/dxc
         Rz = 0.5*(R2-R1)/dxc
      endif
      Rnorm2 = Rx**2 + Ry**2 + Rz**2
      alpha  = - R0/Rnorm2
c
c     Apply secant method:  Use an initial segment twice expected length
c
      x1(1) = x0(1) + 2.0*Rx * alpha
      x1(2) = x0(2) + 2.0*Ry * alpha
      x1(3) = x0(3) + 2.0*Rz * alpha
      call srfind(x1,x0,c,cc)
c
c     write(38,38) x1,x0,cc,c(2),c(3)
c  38 format(3f9.6,1x,3f9.6,1x,a1,2f6.2)
c
      call copy(x0,x1,3)
      return
      end
      subroutine srfind(x1,x0,c,cc)
      real x1(3),x0(3)
      real c(5)
      character*1 cc
c
c     Find point on line segment that intersects the ellipsoid 
c     specified by:
c                       (x/a)**2 + (y/b)**2 + (z/b)**2 = 1
c
c
c     Algorithm:  4 rounds of secant  x_k+1 = x_k - f/f'
c
      a0 = 0.0
      a1 = 1.0
      r0 = ressrf(x0,c,cc)
      dx = x1(1) - x0(1)
      dy = x1(2) - x0(2)
      dz = x1(3) - x0(3)
c     write(6,*) 'dxyz',dx,dy,dz
c     write(6,*) 'cc  ',x0,cc,c(2),c(3)
      do 10 i=1,9
         r1 = ressrf(x1,c,cc)
         if (r1.ne.r0) then
            da = r1*(a1-a0)/(r1-r0)
            r0 = r1
            a0 = a1
            a1 = a1 - da
         endif
         x1(1) = x0(1) + a1*dx
         x1(2) = x0(2) + a1*dy
         x1(3) = x0(3) + a1*dz
   10 continue
c     write(6,*) ' r1',r1,r0,a1
      return
      end
      function ressrf(x,c,cc)
      real x(3) 
      real c(5)
      character*1 cc
c
      ressrf = 0.0
      if (cc.eq.'e') then
         a = c(2)
         b = c(3)
         ressrf = 1.0 - (x(1)/a)**2 - (x(2)/b)**2 - (x(3)/b)**2
         return
      endif
c
      return
      end
