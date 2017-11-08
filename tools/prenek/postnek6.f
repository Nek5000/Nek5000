      SUBROUTINE SETQUAD 
C 
C     Set up the quadrant points and associated element-element links
C     for plotting r-s-t surfaces in 3D and for establishing general
C     zipper planes.        pff 11-20-90  21:30 and  1-8-91  
C
C     (Note: this algorithm is inherently CONFORMING, as constant r-s-t
C            planes have little significance in a non-conforming
C            discretization where the continuum assumption is not honored.)
C
#     include "basics.inc"
C
C     Algorithm:  Set up indices, e.g:
C                                                     ^ s
C         These points define planes                  |          
C         which propagate from one         +----------|---------+
C         element to the next, as          |          |         |
C         shown near pt. 3:                |         4+         |
C                                          |          |         |
C                                          |          |         |
C                                          +--1+------------+2--+--> r
C                                          |          |         |
C                                          |          |         |
C                                          | - - - - 3+ - - - - |
C                                          |          |         |
C                                          +--------------------+
C
C     Set up the connectivity list by performing the K^2 loop which
C     searches for coincident nodes, thus identifying an element-element
C     interaction.
C
CXXX   
CXXX   recall, that normals will always be defined positive if they point
CXXX   away from the element center...,,, thus, once the topology of the 
CXXX   normal (in or out) is established for one element, it will carry through
CXXX   to the rest of the elements.
CXXX   
CXXX   2)  we can have a pseudo XYZ plane for now which takes as the default
CXXX   normal to be x y or z at the given point.
CXXX   
CXXX   3)  when determining which  element contains the originating point, just
CXXX       use the first found point, and choose the closest quadpoint with an
CXXX       appropriately oriented normal vector (in some distance/normal norm).
CXXX   
C
      DIMENSION XYZ(3,8) 
      DIMENSION INEIGH(6,6)
      DATA      INEIGH/  1, -1, -1,  1,  1, -1
     $                , -1,  1,  1, -1, -1,  1
     $                , -1,  1,  1, -1, -1,  1
     $                ,  1, -1, -1,  1,  1, -1
     $                ,  1, -1, -1,  1,  1, -1
     $                , -1,  1,  1, -1, -1,  1 /
C
      IF (.NOT.IF3D) RETURN
      CALL GENCEN
      NFACES=6
      EPSM = 1.0E-10
      NEL5=120*NELM
      CALL RZERO(XYZQ,NEL5)
C     Define quadpoints:
      DO 100 IE=1,NEL
         CALL MAPXC(XYZ,IE)
C
C     Plane 1
         XYZQ(1,1,1,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,5))+XYZ(1,2)+XYZ(1,6))
         XYZQ(2,1,1,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,5))+XYZ(2,2)+XYZ(2,6))
         XYZQ(3,1,1,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,5))+XYZ(3,2)+XYZ(3,6))
C
         XYZQ(1,2,1,IE)=.125*(3.0*(XYZ(1,3)+XYZ(1,7))+XYZ(1,4)+XYZ(1,8))
         XYZQ(2,2,1,IE)=.125*(3.0*(XYZ(2,3)+XYZ(2,7))+XYZ(2,4)+XYZ(2,8))
         XYZQ(3,2,1,IE)=.125*(3.0*(XYZ(3,3)+XYZ(3,7))+XYZ(3,4)+XYZ(3,8))
C
         XYZQ(1,3,1,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,3))+XYZ(1,2)+XYZ(1,4))
         XYZQ(2,3,1,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,3))+XYZ(2,2)+XYZ(2,4))
         XYZQ(3,3,1,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,3))+XYZ(3,2)+XYZ(3,4))
C
         XYZQ(1,4,1,IE)=.125*(3.0*(XYZ(1,5)+XYZ(1,7))+XYZ(1,6)+XYZ(1,8))
         XYZQ(2,4,1,IE)=.125*(3.0*(XYZ(2,5)+XYZ(2,7))+XYZ(2,6)+XYZ(2,8))
         XYZQ(3,4,1,IE)=.125*(3.0*(XYZ(3,5)+XYZ(3,7))+XYZ(3,6)+XYZ(3,8))
C
C     Plane 2
         XYZQ(1,1,2,IE)=.125*(3.0*(XYZ(1,2)+XYZ(1,6))+XYZ(1,1)+XYZ(1,5))
         XYZQ(2,1,2,IE)=.125*(3.0*(XYZ(2,2)+XYZ(2,6))+XYZ(2,1)+XYZ(2,5))
         XYZQ(3,1,2,IE)=.125*(3.0*(XYZ(3,2)+XYZ(3,6))+XYZ(3,1)+XYZ(3,5))
C
         XYZQ(1,2,2,IE)=.125*(3.0*(XYZ(1,4)+XYZ(1,8))+XYZ(1,3)+XYZ(1,7))
         XYZQ(2,2,2,IE)=.125*(3.0*(XYZ(2,4)+XYZ(2,8))+XYZ(2,3)+XYZ(2,7))
         XYZQ(3,2,2,IE)=.125*(3.0*(XYZ(3,4)+XYZ(3,8))+XYZ(3,3)+XYZ(3,7))
C
         XYZQ(1,3,2,IE)=.125*(3.0*(XYZ(1,2)+XYZ(1,4))+XYZ(1,1)+XYZ(1,3))
         XYZQ(2,3,2,IE)=.125*(3.0*(XYZ(2,2)+XYZ(2,4))+XYZ(2,1)+XYZ(2,3))
         XYZQ(3,3,2,IE)=.125*(3.0*(XYZ(3,2)+XYZ(3,4))+XYZ(3,1)+XYZ(3,3))
C
         XYZQ(1,4,2,IE)=.125*(3.0*(XYZ(1,6)+XYZ(1,8))+XYZ(1,5)+XYZ(1,7))
         XYZQ(2,4,2,IE)=.125*(3.0*(XYZ(2,6)+XYZ(2,8))+XYZ(2,5)+XYZ(2,7))
         XYZQ(3,4,2,IE)=.125*(3.0*(XYZ(3,6)+XYZ(3,8))+XYZ(3,5)+XYZ(3,7))
C
C     Plane 5
         XYZQ(1,1,5,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,3))+XYZ(1,5)+XYZ(1,7))
         XYZQ(2,1,5,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,3))+XYZ(2,5)+XYZ(2,7))
         XYZQ(3,1,5,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,3))+XYZ(3,5)+XYZ(3,7))
C
         XYZQ(1,2,5,IE)=.125*(3.0*(XYZ(1,2)+XYZ(1,4))+XYZ(1,6)+XYZ(1,8))
         XYZQ(2,2,5,IE)=.125*(3.0*(XYZ(2,2)+XYZ(2,4))+XYZ(2,6)+XYZ(2,8))
         XYZQ(3,2,5,IE)=.125*(3.0*(XYZ(3,2)+XYZ(3,4))+XYZ(3,6)+XYZ(3,8))
C
         XYZQ(1,3,5,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,2))+XYZ(1,5)+XYZ(1,6))
         XYZQ(2,3,5,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,2))+XYZ(2,5)+XYZ(2,6))
         XYZQ(3,3,5,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,2))+XYZ(3,5)+XYZ(3,6))
C
         XYZQ(1,4,5,IE)=.125*(3.0*(XYZ(1,3)+XYZ(1,4))+XYZ(1,7)+XYZ(1,8))
         XYZQ(2,4,5,IE)=.125*(3.0*(XYZ(2,3)+XYZ(2,4))+XYZ(2,7)+XYZ(2,8))
         XYZQ(3,4,5,IE)=.125*(3.0*(XYZ(3,3)+XYZ(3,4))+XYZ(3,7)+XYZ(3,8))
C
C     Plane 6
         XYZQ(1,1,6,IE)=.125*(3.0*(XYZ(1,5)+XYZ(1,7))+XYZ(1,1)+XYZ(1,3))
         XYZQ(2,1,6,IE)=.125*(3.0*(XYZ(2,5)+XYZ(2,7))+XYZ(2,1)+XYZ(2,3))
         XYZQ(3,1,6,IE)=.125*(3.0*(XYZ(3,5)+XYZ(3,7))+XYZ(3,1)+XYZ(3,3))
C
         XYZQ(1,2,6,IE)=.125*(3.0*(XYZ(1,6)+XYZ(1,8))+XYZ(1,2)+XYZ(1,4))
         XYZQ(2,2,6,IE)=.125*(3.0*(XYZ(2,6)+XYZ(2,8))+XYZ(2,2)+XYZ(2,4))
         XYZQ(3,2,6,IE)=.125*(3.0*(XYZ(3,6)+XYZ(3,8))+XYZ(3,2)+XYZ(3,4))
C
         XYZQ(1,3,6,IE)=.125*(3.0*(XYZ(1,5)+XYZ(1,6))+XYZ(1,1)+XYZ(1,2))
         XYZQ(2,3,6,IE)=.125*(3.0*(XYZ(2,5)+XYZ(2,6))+XYZ(2,1)+XYZ(2,2))
         XYZQ(3,3,6,IE)=.125*(3.0*(XYZ(3,5)+XYZ(3,6))+XYZ(3,1)+XYZ(3,2))
C
         XYZQ(1,4,6,IE)=.125*(3.0*(XYZ(1,7)+XYZ(1,8))+XYZ(1,3)+XYZ(1,4))
         XYZQ(2,4,6,IE)=.125*(3.0*(XYZ(2,7)+XYZ(2,8))+XYZ(2,3)+XYZ(2,4))
         XYZQ(3,4,6,IE)=.125*(3.0*(XYZ(3,7)+XYZ(3,8))+XYZ(3,3)+XYZ(3,4))
C
C     Plane 3
         XYZQ(1,1,3,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,5))+XYZ(1,3)+XYZ(1,7))
         XYZQ(2,1,3,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,5))+XYZ(2,3)+XYZ(2,7))
         XYZQ(3,1,3,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,5))+XYZ(3,3)+XYZ(3,7))
C
         XYZQ(1,2,3,IE)=.125*(3.0*(XYZ(1,2)+XYZ(1,6))+XYZ(1,4)+XYZ(1,8))
         XYZQ(2,2,3,IE)=.125*(3.0*(XYZ(2,2)+XYZ(2,6))+XYZ(2,4)+XYZ(2,8))
         XYZQ(3,2,3,IE)=.125*(3.0*(XYZ(3,2)+XYZ(3,6))+XYZ(3,4)+XYZ(3,8))
C
         XYZQ(1,3,3,IE)=.125*(3.0*(XYZ(1,1)+XYZ(1,2))+XYZ(1,3)+XYZ(1,4))
         XYZQ(2,3,3,IE)=.125*(3.0*(XYZ(2,1)+XYZ(2,2))+XYZ(2,3)+XYZ(2,4))
         XYZQ(3,3,3,IE)=.125*(3.0*(XYZ(3,1)+XYZ(3,2))+XYZ(3,3)+XYZ(3,4))
C
         XYZQ(1,4,3,IE)=.125*(3.0*(XYZ(1,5)+XYZ(1,6))+XYZ(1,7)+XYZ(1,8))
         XYZQ(2,4,3,IE)=.125*(3.0*(XYZ(2,5)+XYZ(2,6))+XYZ(2,7)+XYZ(2,8))
         XYZQ(3,4,3,IE)=.125*(3.0*(XYZ(3,5)+XYZ(3,6))+XYZ(3,7)+XYZ(3,8))
C
C     Plane 4
         XYZQ(1,1,4,IE)=.125*(3.0*(XYZ(1,3)+XYZ(1,7))+XYZ(1,1)+XYZ(1,5))
         XYZQ(2,1,4,IE)=.125*(3.0*(XYZ(2,3)+XYZ(2,7))+XYZ(2,1)+XYZ(2,5))
         XYZQ(3,1,4,IE)=.125*(3.0*(XYZ(3,3)+XYZ(3,7))+XYZ(3,1)+XYZ(3,5))
C
         XYZQ(1,2,4,IE)=.125*(3.0*(XYZ(1,4)+XYZ(1,8))+XYZ(1,2)+XYZ(1,6))
         XYZQ(2,2,4,IE)=.125*(3.0*(XYZ(2,4)+XYZ(2,8))+XYZ(2,2)+XYZ(2,6))
         XYZQ(3,2,4,IE)=.125*(3.0*(XYZ(3,4)+XYZ(3,8))+XYZ(3,2)+XYZ(3,6))
C
         XYZQ(1,3,4,IE)=.125*(3.0*(XYZ(1,3)+XYZ(1,4))+XYZ(1,1)+XYZ(1,2))
         XYZQ(2,3,4,IE)=.125*(3.0*(XYZ(2,3)+XYZ(2,4))+XYZ(2,1)+XYZ(2,2))
         XYZQ(3,3,4,IE)=.125*(3.0*(XYZ(3,3)+XYZ(3,4))+XYZ(3,1)+XYZ(3,2))
C
         XYZQ(1,4,4,IE)=.125*(3.0*(XYZ(1,7)+XYZ(1,8))+XYZ(1,5)+XYZ(1,6))
         XYZQ(2,4,4,IE)=.125*(3.0*(XYZ(2,7)+XYZ(2,8))+XYZ(2,5)+XYZ(2,6))
         XYZQ(3,4,4,IE)=.125*(3.0*(XYZ(3,7)+XYZ(3,8))+XYZ(3,5)+XYZ(3,6))
C
  100 CONTINUE
C
C     Find Neighbor-Neighbor connections:
C
      DO 500 IE=1,NEL
      DO 500 IFACE=1,NFACES
      DO 500 IQPT=1,4
C
         NEIGHI=INT(XYZQ(4,IQPT,IFACE,IE))
         IF (NEIGHI.EQ.0) THEN
            IQP1=IQPT+1
            IQP1=MOD1(IQP1,4)
            EPSI=(XYZQ(1,IQPT,IFACE,IE)-XYZQ(1,IQP1,IFACE,IE))**2
     $          +(XYZQ(2,IQPT,IFACE,IE)-XYZQ(2,IQP1,IFACE,IE))**2
     $          +(XYZQ(3,IQPT,IFACE,IE)-XYZQ(3,IQP1,IFACE,IE))**2
C
c           DO 204 JE=1,NEL
c   pff: 4/13/93
            DO 204 JE=ie,Nel
            IF (JE.EQ.IE) GOTO 202
            DIST=(XCEN(IE)-XCEN(JE))**2+(YCEN(IE)-YCEN(JE))**2
     $          +(ZCEN(IE)-ZCEN(JE))**2
            DST2=(RCEN(IE)+RCEN(JE))**2
c           write(6,101) ie,je,dst2,dist,rcen(ie),rcen(je)
  101       format(2x,'dist:',2i3,4f12.7)
            IF (DST2.LT.DIST) GOTO 202
            DO 200 JFACE=1,NFACES
            DO 200 JQPT=1,4
c   pff: 4/13/93
            NEIGHJ=INT(XYZQ(4,IQPT,IFACE,IE))
            IF (NEIGHJ.ne.0) GOTO 200
c   pff: 4/13/93
               JQP1=JQPT+1
               JQP1=MOD1(JQP1,4)
               EPSJ=(XYZQ(1,JQPT,JFACE,JE)-XYZQ(1,JQP1,JFACE,JE))**2
     $             +(XYZQ(2,JQPT,JFACE,JE)-XYZQ(2,JQP1,JFACE,JE))**2
     $             +(XYZQ(3,JQPT,JFACE,JE)-XYZQ(3,JQP1,JFACE,JE))**2
               EPS=EPSM*(EPSJ+EPSI)
               EPS=MAX(EPS,EPSM)
               DIST=(XYZQ(1,IQPT,IFACE,IE)-XYZQ(1,JQPT,JFACE,JE))**2
     $             +(XYZQ(2,IQPT,IFACE,IE)-XYZQ(2,JQPT,JFACE,JE))**2
     $             +(XYZQ(3,IQPT,IFACE,IE)-XYZQ(3,JQPT,JFACE,JE))**2
               IF (DIST.LE.EPS) THEN
C                 We found our match, make a note of it 
C                 and move to the next point.
                  NEIGHI=24*(JE-1)+4*(JFACE-1)+JQPT
                  NEIGHJ=24*(IE-1)+4*(IFACE-1)+IQPT
                  INRM=INEIGH(IFACE,JFACE)
                  XYZQ(4,NEIGHJ,1,1)=FLOAT(INRM*NEIGHI)
                  XYZQ(4,NEIGHI,1,1)=FLOAT(INRM*NEIGHJ)
                  GOTO 300
               ENDIF
  200       CONTINUE
  202       CONTINUE
  204       CONTINUE
  300       CONTINUE
         ENDIF
  500 CONTINUE
C
C     Diagnostic (1-8-91)
C
c     CALL DRQUAD
C
      RETURN
      END
      SUBROUTINE DRQUAD
#     include "basics.inc"
C
      NFACES=6
      DO 600 IE=1,NEL
      DO 600 IFACE=1,NFACES
      DO 600 IQPT=1,4
         ICLR=2*IFACE
c        ICLR=MOD1(ICLR,14)
c        ICLR=5
         CALL COLOR(ICLR)
         CALL XDIAM3(XYZQ(1,IQPT,IFACE,IE))
  600 CONTINUE
C
      RETURN
      END
      SUBROUTINE XDIAM3(XYZ)
c     draws a diamond around XX,YY
#     include "basics.inc"
C
      DIMENSION XYZ(3)
      XX=XISO(XYZ(1),XYZ(2),XYZ(3))
      YY=YISO(XYZ(1),XYZ(2),XYZ(3))
      CALL MOVEC(XX+.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWC(XX-.0025*XFAC ,YY+.0025*YFAC)
      CALL DRAWC(XX-.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWC(XX+.0025*XFAC ,YY-.0025*YFAC)
      CALL DRAWC(XX+.0025*XFAC ,YY+.0025*YFAC)
      RETURN
      END
      SUBROUTINE MAPXC(XYZ,IE)
C     Map the element corner points (x,y,z) to standard vector.
#     include "basics.inc"
      DIMENSION INDX(8)
      DIMENSION XYZ(3,8)
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
      DO 50 IX=1,NDIM2
         I=INDX(IX)
         xyz(1,ix)=x(i,ie)
         xyz(2,ix)=y(i,ie)
         xyz(3,ix)=z(i,ie)
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DECOD(IX,IY,IZ,IE,IXYZ,NX,NY,NZ)
C
      NXY  = NX*NY
      NXYZ = NXY*NZ
C
      IE  = (IXYZ + NXYZ - 1)/NXYZ
      IXY = IXYZ-(IE-1)*NXYZ
      IZ  = (IXY + NXY - 1)/NXY
      IXX = IXY-(IZ-1)*NXY
      IY  = (IXX + NX - 1)/NX
      IX  = IXX-(IY-1)*NX
C
C     Check
C
      ICHECK = (IE-1)*NXYZ + (IZ-1)*NXY + (IY-1)*NX + IX
      IF (ICHECK.NE.IXYZ) THEN
         WRITE(6,*) 'hey, what gives here?  icheck,ixyz:',icheck,ixyz
      ENDIF
      RETURN
      END
      SUBROUTINE SETSRFP(SRFQAL,NQUAL,IFLD)
      CHARACTER*1 SRFQAL(1)
      CHARACTER*1 CB1,CB
C
#     include "basics.inc"
C
C     Select a list of planes according to the specified boundary
C     qualifications given in SRFQAL.  If the number of qualifications,
C     NQUAL, is negative, select all boundaries NOT meeting the 
C     qualifications.
C
C     Note:  Modification 10/2/92 -
C            Set INRM such that the "+" direction associated with the
C            plane has the same direction as the closest fundamental
C            unit vector (e_1,e_2, or e_3).    pff.
C
C
c     NLSTP=0                      pff 12/21/96
      NFACES=2*NDIM
C
      IF (NQUAL.GT.0) THEN
C
         DO 100 IQ=1,NQUAL
         CB1=SRFQAL(IQ)
         DO 100 IE=1,NEL
         DO 100 IFACE=1,NFACES
            CB=CBC(IFACE,IE,IFLD)
            IF (CB.EQ.CB1) THEN
               NLSTP=NLSTP+1
               IJKPLN=EFACE1(IFACE)
               IX=1
               IF (MOD(IJKPLN,2).EQ.0) IX=NX
               LIST(NLSTP)=6*NX*(IE-1)+NX*(IJKPLN-1)+IX
               CALL SIGNPL(LIST(NLSTP))
            ENDIF
  100    CONTINUE
C
      ELSE
C
         NQ=-NQUAL
         DO 300 IE=1,NEL
         DO 300 IFACE=1,NFACES
            CB=CBC(IFACE,IE,IFLD)
            DO 200 IQ=1,NQ
               CB1=SRFQAL(IQ)
C              if it's one of these, don't plot.
               IF (CB1.EQ.CB) GOTO 300
  200       CONTINUE
C
C           Since the b.c. isn't one of the one's specified, we plot
C           this surface....
            NLSTP=NLSTP+1
            IJKPLN=EFACE1(IFACE)
            IX=1
            IF (MOD(IJKPLN,2).EQ.0) IX=NX
            LIST(NLSTP)=6*NX*(IE-1)+NX*(IJKPLN-1)+IX
            CALL SIGNPL(LIST(NLSTP))
  300    CONTINUE
C
      ENDIF
      RETURN
      END
      SUBROUTINE SORTL
#     include "basics.inc"
      INCLUDE 'basicsp.inc'
      COMMON /CTMPS/ IND(NELM,6),WK(NELM,6)
C
      NXY =NX*NY
      NXYZ=NZ*NXY
      DO 100 I=1,NLSTP
C        Get an XYZ point close to the center of the requested sub-plane.
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         IX=(NX+1)/2
         IY=(NY+1)/2
         IZ=(NZ+1)/2
         IF (IPLN.LE.2) THEN
            IX=IPLANE
         ELSEIF (IPLN.LE.4) THEN
            IY=IPLANE
         ELSE
            IZ=IPLANE
         ENDIF
         IEOFF=NXYZ*(IE-1)+NXY*(IZ-1)+NX*(IY-1)+IX
         ZDEPTH(I,1)=ZISO(XP(IEOFF),YP(IEOFF),ZP(IEOFF))
100   CONTINUE
      CALL SORT(ZDEPTH,   IND,NLSTP)
      CALL ISWAP(LIST ,WK,IND,NLSTP)
C
      RETURN
      END
      SUBROUTINE SIGNPL(Lplane)
C
C     Find the closest component of (1,1,1) to the normal associated
C     with this plane, and assign the sign of the dot product of that
C     component with the normal to Lplane.
C     11/92 pff.
C
#     include "basics.inc"
      INCLUDE 'basicsp.inc'
      DIMENSION VEC0(3),VEC1(3),VEC2(3),VEC3(3)
C
      CALL DECOD(IPLANE,IPLN,IE,IDUM,Lplane,NX,6,NELM)
C
      I=(1+NX)/2
      J=(1+NY)/2
      K=(1+NZ)/2
C
      Nxyz=Nx*Ny*Nz
      Nxy =Nx*Ny
C
      IF (IPLN.LE.2) THEN
         I =IPLANE
         ieoff1 = nxyz*(ie-1)+nxy*(k-1)+nx*(j  )+i
         ieoff2 = nxyz*(ie-1)+nxy*(k  )+nx*(j-1)+i
      ELSEIF (IPLN.LE.4) THEN
         J =IPLANE
         ieoff1 = nxyz*(ie-1)+nxy*(k-1)+nx*(j-1)+i+1
         ieoff2 = nxyz*(ie-1)+nxy*(k  )+nx*(j-1)+i
      ELSE
         K =IPLANE
         ieoff1 = nxyz*(ie-1)+nxy*(k-1)+nx*(j-1)+i+1
         ieoff2 = nxyz*(ie-1)+nxy*(k-1)+nx*(j  )+i
      ENDIF
      ieoff  = nxyz*(ie-1)+nxy*(k-1)+nx*(j-1)+i
C
      VEC0(1)=Xp(ieoff)
      VEC0(2)=Yp(ieoff)
      VEC0(3)=Zp(ieoff)
C
      VEC1(1)=Xp(ieoff1)
      VEC1(2)=Yp(ieoff1)
      VEC1(3)=Zp(ieoff1)
C
      VEC2(1)=Xp(ieoff2)
      VEC2(2)=Yp(ieoff2)
      VEC2(3)=Zp(ieoff2)
C
      CALL SUB3  (VEC1,VEC1,VEC0,3)
      CALL SUB3  (VEC2,VEC2,VEC0,3)
      CALL CROSS (VEC0,VEC1,VEC2)
      CALL NORM3D(VEC0)
C
      Vmax=0.0
      DO 10 I=1,3
         IF (ABS(VEC0(i)).GT.Vmax) THEN
             Vmax=ABS(VEC0(i))
             Imax=i
         ENDIF
   10 CONTINUE
      IF (VEC0(Imax).LT.0.0) Lplane=-Lplane
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine gencen
#     include "basics.inc"

      integer e

      common /ctmp2/ xp(nxm,nxm,nxm),yp(nxm,nxm,nxm),zp(nxm,nxm,nxm)

      if (mod(nxm,2).eq.0) then
         write(6,*) 'ERROR: Recompile with nxm odd in basics.inc'
         call prexit(0)
      endif

      nxh = (nxm+1)/2
      nh1 = nxh-1

      write(6,*) 'inside gencen ',nel,if3d,nxh

C     Generate the element centers
      do e=1,nel
c        call genxyz_e (xp,yp,zp,e,3,3,3)
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

c        call out27(x27(1,e),y27(1,e),z27(1,e),e,'genc')

      enddo

C     Compute the maximum radius from the center

      call rzero(rcen,nel) 
      if (if3d) then
         do 300 e=1,nel
         do 300 j=1,8
            rad=(x(j,e)-xcen(e))**2 + (y(j,e)-ycen(e))**2
     $         +(z(j,e)-zcen(e))**2
            rcen(e)=max(rcen(e),rad)
  300    continue
      else
         do 400 e=1,nel
         do 400 j=1,4
            rad=(x(j,e)-xcen(e))**2 + (y(j,e)-ycen(e))**2
            rcen(e)=max(rcen(e),rad)
  400    CONTINUE
      ENDIF
      call vsqrt(rcen,nel)

      write(6,*) 'done gencen ',nel

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_neigh
#     include "basics.inc"
c
      call izero(neighb,6*nelm)
      call mkside
      call gencen
      nsides = 2*ndim
      do je=1   ,nel
      do ie=je+1,nel
         d2 = sqrt( ( xcen(ie)-xcen(je) )**2 +
     $              ( ycen(ie)-ycen(je) )**2 +
     $              ( zcen(ie)-zcen(je) )**2 )
         if (d2.gt.(rcen(ie)+rcen(je))) goto 360

         delta = .01 * min(rcen(ie),rcen(je))
C     
         do iside=1,nsides
         do jside=1,nsides
            deltax = abs(sides(ie,iside,1)-sides(je,jside,1))
            deltay = abs(sides(ie,iside,2)-sides(je,jside,2))
            deltaz = abs(sides(ie,iside,3)-sides(je,jside,3))
            if (deltax .lt. delta .and.
     $          deltay .lt. delta .and.
     $          deltaz .lt. delta ) then
                  neighb(iside,ie) = je
                  neighb(jside,je) = ie
            endif
         enddo
         enddo
  360    continue
      enddo
      enddo
c
c     do ie=1,nel
c        write(6,*) ie,(neighb(j,ie),j=1,4),' neig'
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine out27(x,y,z,e,name4)
      real x(3,3,3),y(3,3,3),z(3,3,3)
      integer e
      character*4 name4
      integer icalld
      save    icalld
      data    icalld /0/

      return

      icalld = icalld+1

      do k=1,3
      do j=1,3
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),i,i,j,k,e,name4,i=1,3)
      enddo
      enddo


      write(6,*)
      write(6,*) i,j,k,e,icalld,'XYZ27 ',name4,i,ncurve
   3  format(1p3e12.4,4i3,i7,' xyz27 ',a4)

c     X-lines

      do k=1,3,2
      do j=1,3
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),i,i,j,k,e,name4,i=1,3)
      enddo
      enddo

      do k=2,2
      do j=1,3,2
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),i,i,j,k,e,name4,i=1,3)
      enddo
      enddo

c     Y-lines

      do k=1,3,2
      do i=1,3
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),j,i,j,k,e,name4,j=1,3)
      enddo
      enddo

      do k=2,2
      do i=1,3,2
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),j,i,j,k,e,name4,j=1,3)
      enddo
      enddo

c     Y-lines

      do j=1,3,2
      do i=1,3
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),k,i,j,k,e,name4,k=1,3)
      enddo
      enddo

      do j=2,2
      do i=1,3,2
         write(6,*)
         write(6,3) (x(i,j,k),y(i,j,k),z(i,j,k),k,i,j,k,e,name4,k=1,3)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
