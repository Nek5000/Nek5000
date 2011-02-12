      SUBROUTINE ZIP2
      INCLUDE 'basics.inc'
      DIMENSION ENEW(NELM),IND(NELM)
C
C     This value can be changed later.
C
      NX=9
      NY=NX
      NZ=NX
      CALL LEGEND(ZPTS,WGHT,NZ)
C
C     Set RSTP planes will find the set of dividing planes for the zipper.
C
      CALL SETQUAD 
      CALL SETRST2(0)
C
C     Check on the number and numbering of the new elements to be generated.
C
      NELN=NEL+NLSTP
      NELM1=NELM-1
      IF (NELN.GT.NELM1) THEN
C
         WRITE(S,10) 
         CALL PRS(S)
         WRITE(S,11) NELN,NELM1
         CALL PRS(S)
         WRITE(S,12) 
         CALL PRS(S)
   10  FORMAT(2X,'WARNING: Number of elements after zipper operation$')
   11  FORMAT(2X,'(',I5,') would be greater than the allowed maxium ('
     $                ,I5,').$')
   12    FORMAT(2X,'Aborting zipper.$')
         RETURN
      ENDIF
C
C     Renumber ALL elements, so that low numbered elements will remain
C     low numbered.  This will be achieved by assigning IE+0.5 to the
C     new element numbers, and then sorting the list.
C
      DO 100 IE=1,NEL
         ENEW(IE)=IE
  100 CONTINUE
      DO 101 I=1,NLSTP
         NELN=NEL+I
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         ENEW(NELN)=FLOAT(IE)+0.5
  101 CONTINUE
C
C     Generate new sub-elements as a result of the zipper action.
C
      S0=0.5
      S1=1.0-S0
      DO 200 I=1,NLSTP
         NELN=NEL+I
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         CALL SPLITE(NELN,IE,I,S1)
  200 CONTINUE
      NEL=NELN
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
C
C     Exit
C
      WRITE(S,300) NLSTP,NEL
  300 FORMAT(I5,' new elements generated. NEL=',I5,'$')
      CALL PRS(S)
C
      CALL GENCEN
C
      RETURN
      END
      SUBROUTINE SPLITE(JE,IE,ILIST,FRAC)
C
C     A subroutine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
C
      IF (FRAC.LE.0.0.OR.FRAC.GE.1.0) RETURN
C
C     Begin with copying the essential information
C
      CALL COPYEL(IE,JE)
C
C     Compute XYZ on high definition mesh
C
      CALL GENXYZ (XP,YP,ZP,IE,1,NX,NY,NZ)
C
C     Are we taking an R, S or T cut?
C
      FRACS=2.0*FRAC-1.0
      LISTA=ABS(LIST(ILIST))
      CALL DECOD(IPLANE,IPLN,IIE,IDUM,LISTA,NX,6,NELM)
C
      IF (IPLN.LE.2) THEN
C        X-plane
         I0=-1
         RRL(1)=FRACS
         DO 200 JCRN=-1,1,2
            RRL(3)=FLOAT(JCRN)
            DO 100 ICRN=-1,1,2
               RRL(2)=FLOAT(ICRN)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
               I0=I0+2
               I1=I0+1
               X(IE,INV(I1))=XVAL
               Y(IE,INV(I1))=YVAL
               Z(IE,INV(I1))=ZVAL
               X(JE,INV(I0))=XVAL
               Y(JE,INV(I0))=YVAL
               Z(JE,INV(I0))=ZVAL
  100       CONTINUE
  200    CONTINUE
C
C        Don't forget to eliminate Curve side data on new faces!
C
         JFAC1=EFACE(1)
         JFAC2=EFACE(2)
         IF (CCURVE(JFAC1,IE).EQ.'s') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'s') CCURVE(JFAC2,IE)=' '
         IF (CCURVE(JFAC1,IE).EQ.'C') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'C') CCURVE(JFAC2,IE)=' '
         JFAC1=JFAC1+4
         JFAC2=JFAC2+4
         IF (CCURVE(JFAC1,IE).EQ.'C') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'C') CCURVE(JFAC2,IE)=' '
C
C
      ELSEIF (IPLN.LE.4) THEN
C        Y-plane
         RRL(2)=FRACS
         DO 210 JCRN=-1,1,2
            RRL(3)=FLOAT(JCRN)
            DO 110 ICRN=-1,1,2
               RRL(1)=FLOAT(ICRN)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
C              II,JJ = 1,2
               II=1+(ICRN+1)/2
               JJ=1+(JCRN+1)/2
C              I0=1,2,5,6; I1=3,4,7,8
               I0=II+4*(JJ-1)
               I1=I0+2
               X(IE,INV(I1))=XVAL
               Y(IE,INV(I1))=YVAL
               Z(IE,INV(I1))=ZVAL
               X(JE,INV(I0))=XVAL
               Y(JE,INV(I0))=YVAL
               Z(JE,INV(I0))=ZVAL
C
  110       CONTINUE
  210    CONTINUE
C
         JFAC1=EFACE(3)
         JFAC2=EFACE(4)
         IF (CCURVE(JFAC1,IE).EQ.'s') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'s') CCURVE(JFAC2,IE)=' '
         IF (CCURVE(JFAC1,IE).EQ.'C') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'C') CCURVE(JFAC2,IE)=' '
         JFAC1=JFAC1+4
         JFAC2=JFAC2+4
         IF (CCURVE(JFAC1,IE).EQ.'C') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'C') CCURVE(JFAC2,IE)=' '
C
      ELSE
C        Z-plane
         I0=0
         RRL(3)=FRACS
         DO 220 JCRN=-1,1,2
            RRL(2)=FLOAT(JCRN)
            DO 120 ICRN=-1,1,2
               RRL(1)=FLOAT(ICRN)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
               I0=I0+1
               I1=I0+4
               X(IE,INV(I1))=XVAL
               Y(IE,INV(I1))=YVAL
               Z(IE,INV(I1))=ZVAL
               X(JE,INV(I0))=XVAL
               Y(JE,INV(I0))=YVAL
               Z(JE,INV(I0))=ZVAL
C
  120       CONTINUE
  220    CONTINUE
         JFAC1=EFACE(5)
         JFAC2=EFACE(6)
         IF (CCURVE(JFAC1,IE).EQ.'s') CCURVE(JFAC1,JE)=' '
         IF (CCURVE(JFAC2,IE).EQ.'s') CCURVE(JFAC2,IE)=' '
C
      ENDIF
C
C     Delete periodic bcs and undo curved sides
C
      RETURN
      END
      SUBROUTINE SWAPEL(IND,N,NELM1)
      DIMENSION IND(1)
C***
C***  SWAP ELEMENTS  NOTE:  Index tells us that the 1st entry
C***                        can be found in slot IND(1), etc.
C***
C
      DO 100 I=1,N
C
C        Put element IE into slot I
C
         IE=IND(I)
         IF (IE.NE.I) THEN
C
C           Move element I to the end
            CALL COPYEL(I,NELM1)
C           Move correct element to Ith location
            CALL COPYEL(IE,I  )
C           Copy temporary data to old slot.
            CALL COPYEL(NELM1,IE)
C           Update IND:
            DO 20 J=1,N
               IF (IND(J).EQ.I) THEN
                   IND(J)=IE
                   GOTO 21
               ENDIF
   20       CONTINUE
   21       CONTINUE
         ENDIF
  100    CONTINUE
      RETURN
      END
      SUBROUTINE GENXYZ (XML,YML,ZML,IE,IEL,NXL,NYL,NZL)
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
      CALL RZERO(XML(1,1,1,IEL),NXYZ)
      CALL RZERO(YML(1,1,1,IEL),NXYZ)
      CALL RZERO(ZML(1,1,1,IEL),NXYZ)
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
      DO 50 IX=1,NDIM2
         I=INDX(IX)
         XCB(IX,1,1)=X(IE,I)
         YCB(IX,1,1)=Y(IE,I)
         ZCB(IX,1,1)=Z(IE,I)
   50 CONTINUE
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
            XML(IX,IY,IZ,IEL)=XML(IX,IY,IZ,IEL)+HHH*XCB(IXT,IYT,IZT)
            YML(IX,IY,IZ,IEL)=YML(IX,IY,IZ,IEL)+HHH*YCB(IXT,IYT,IZT)
            ZML(IX,IY,IZ,IEL)=ZML(IX,IY,IZ,IEL)+HHH*ZCB(IXT,IYT,IZT)
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
     $     CALL SPHSRF(XML,YML,ZML,IFACE,IE,IEL,NXL,NYL,NZL,WORK) 
 1000 CONTINUE
C
      DO 2000 ISID=1,8
        CCV = CCURVE(ISID,IE)
        IF (CCV.EQ.'C') 
     $  CALL ARCSRF(XML,YML,ZML,NXL,NYL,NZL,IE,IEL,ISID)
 2000 CONTINUE
C
      RETURN
      END
      SUBROUTINE SPHSRF(XML,YML,ZML,IFCE,IE,IEL,MX,MY,MZ,XYSRF) 
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
      DO 10 I=1,8
         XCC(I)=X(IE,I)
         YCC(I)=Y(IE,I)
         ZCC(I)=Z(IE,I)
   10 CONTINUE
      CALL CRN3D(XCV,XCC,YCC,ZCC,CURVE(1,IFCE,IE),IFACE,IE)
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
      DELT    = 0.0
      I=0
      DO 700 J2=JS2,JF2,JSKIP2
      DO 700 J1=JS1,JF1,JSKIP1
         I=I+1
         JOPP = J1 + IOPP(IFC2)*IDIR
         X2(1) = XML(J1,J2,1,IEL)
         X2(2) = YML(J1,J2,1,IEL)
         X2(3) = ZML(J1,J2,1,IEL)
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
            XML(J,J2,1,IEL) = XML(J,J2,1,IEL)+DX(1)*ZETA
            YML(J,J2,1,IEL) = YML(J,J2,1,IEL)+DX(2)*ZETA
            ZML(J,J2,1,IEL) = ZML(J,J2,1,IEL)+DX(3)*ZETA
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
      REAL U1(3),U2(3),VN(3),B(3)
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
      SUBROUTINE CRN3D(XCV,XC,YC,ZC,CURVE,IFACE,IE)
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
  20     FORMAT(5X,'ERROR: Sphere of radius zero requested.',/,
     $         ,5X,'EXITING in CRN3De',3E12.4,I3)
         CALL EXIT
      ELSE
         DO 40 I=1,4
            RADT=XCV(1,I,1)**2+XCV(2,I,1)**2+XCV(3,I,1)**2
            RADT=SQRT(RADT)
            TEST=ABS(RADT-RADIUS)/RADIUS
            IF (TEST.GT.EPS) THEN
             WRITE(6,30) 
     $       RADT,RADIUS,XCV(1,I,1),XCV(2,I,1),XCV(3,I,1)
   30        FORMAT(5X,'ERROR: vertex not on requested sphere B.'
     $           ,/,5X,'EXITING in CRN3Df',5F12.7)
             WRITE(6,31) IE,IFACE,XCTR,YCTR,ZCTR
   31        FORMAT(5X,'IE,IF,XYZCTR:',2I4,3F12.7)
             WRITE(6,32) (xc(j),yc(j),zc(j),j=1,8)
   32        FORMAT(3f12.7)
             CALL EXIT
            ENDIF
   40    CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE DSSET(MX,MY,MZ)
C
C     Set up arrays ESKIP,SKPDAT,NEDG,NOFFST for new MX,MY,MZ
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
      SUBROUTINE ARCSRF(XML,YML,ZML,NXL,NYL,NZL,IE,IEL,ISID)
      INCLUDE 'basics.inc'
C
C     ....note..... CTMP1 is used in this format in several subsequent routines
C
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMP1/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      DIMENSION XML(NXL,NYL,NZL,1),YML(NXL,NYL,NZL,1),ZML(NXL,NYL,NZL,1)
      LOGICAL IFGLJ
C
      IFGLJ = .FALSE.
c     IF (IFAXIS .AND. IFRZER(IE) .AND. (ISID.EQ.2 .OR. ISID.EQ.4)) 
c    $IFGLJ = .TRUE.
C
      PT1X  = X(IE,ISID)
      PT1Y  = Y(IE,ISID)
      IF(ISID.EQ.4) THEN
         PT2X = X(IE,1)
         PT2Y = Y(IE,1)
      ELSE IF(ISID.EQ.8) THEN
         PT2X = X(IE,5)
         PT2Y = Y(IE,5)
      ELSE
         PT2X = X(IE,ISID+1)
         PT2Y = Y(IE,ISID+1)
      ENDIF
C
C     Find slope of perpendicular
      RADIUS=CURVE(1,ISID,IE)
      GAP=SQRT( (PT1X-PT2X)**2 + (PT1Y-PT2Y)**2 )
      IF (ABS(2.0*RADIUS).LE.GAP*1.00001) THEN
         WRITE(6,10) RADIUS,ISID,IE,GAP
   10    FORMAT(//,2X,'ERROR: Too small a radius (',G11.3
     $  ,') specified for side',I2,' of element',I4,':  '
     $  ,G11.3,/,2X,'ABORTING during mesh generation.')
         CALL EXIT
      ENDIF
      XS = PT2Y-PT1Y
      YS = PT1X-PT2X
C     Make length Radius
      XYS=SQRT(XS**2+YS**2)
C     Find Center
      DTHETA = ABS(ASIN(0.5*GAP/RADIUS))
      PT12X  = (PT1X + PT2X)/2.0
      PT12Y  = (PT1Y + PT2Y)/2.0
      XCENN  = PT12X - XS/XYS * RADIUS*COS(DTHETA)
      YCENN  = PT12Y - YS/XYS * RADIUS*COS(DTHETA)
      THETA0 = ATAN2((PT12Y-YCENN),(PT12X-XCENN))
      IF (IFGLJ) THEN
         FAC    = SIGN(1.0,RADIUS)
         THETA1 = THETA0 - FAC*DTHETA
         THETA2 = THETA0 + FAC*DTHETA
      ENDIF
C     Compute perturbation of geometry
      ISID1 = MOD1(ISID,4)
      IF (IFGLJ) THEN
         I1 = ISID/2
         I2 = 2 - ISID/4
         DO 15 IY=1,NYL
           ANG  = H(IY,2,I1)*THETA1 + H(IY,2,I2)*THETA2
           XCRVED(IY)=XCENN + ABS(RADIUS)*COS(ANG)
     $                      - (H(IY,2,I1)*PT1X + H(IY,2,I2)*PT2X)
           YCRVED(IY)=YCENN + ABS(RADIUS) * SIN(ANG)
     $                      - (H(IY,2,I1)*PT1Y + H(IY,2,I2)*PT2Y)
   15    CONTINUE
      ELSE
         DO 20 IX=1,NXL
            IXT=IX
            IF (ISID1.GT.2) IXT=NXL+1-IX
            R=ZGML(IX,1)
            IF (RADIUS.LT.0.0) R=-R
            XCRVED(IXT) = XCENN + ABS(RADIUS) * COS(THETA0 + R*DTHETA)
     $                          - ( H(IX,1,1)*PT1X + H(IX,1,2)*PT2X )
            YCRVED(IXT) = YCENN + ABS(RADIUS) * SIN(THETA0 + R*DTHETA)
     $                          - ( H(IX,1,1)*PT1Y + H(IX,1,2)*PT2Y )
   20    CONTINUE
      ENDIF
C     Points all set, add perturbation to current mesh.
      ISID1 = MOD1(ISID,4)
      ISID1 = EFACE1(ISID1)
      IZT = (ISID-1)/4+1
      IYT = ISID1-2
      IXT = ISID1
      IF (ISID1.LE.2) THEN
         CALL ADDTNSR(XML(1,1,1,IEL),H(1,1,IXT),XCRVED,H(1,3,IZT)
     $               ,NXL,NYL,NZL)
         CALL ADDTNSR(YML(1,1,1,IEL),H(1,1,IXT),YCRVED,H(1,3,IZT)
     $               ,NXL,NYL,NZL)
      ELSE
         CALL ADDTNSR(XML(1,1,1,IEL),XCRVED,H(1,2,IYT),H(1,3,IZT)
     $               ,NXL,NYL,NZL)
         CALL ADDTNSR(YML(1,1,1,IEL),YCRVED,H(1,2,IYT),H(1,3,IZT)
     $               ,NXL,NYL,NZL)
      ENDIF
      RETURN
      END
      SUBROUTINE ADDTNSR(S,H1,H2,H3,NX,NY,NZ)
C
C     Map and add to S a tensor product form of the three functions H1,H2,H3.
C     This is a single element routine used for deforming geometry.
C
      DIMENSION H1(1),H2(1),H3(1)
      DIMENSION S(NX,NY,NZ)
C
      DO 200 IZ=1,NZ
      DO 200 IY=1,NY
         HH = H2(IY)*H3(IZ)
         DO 100 IX=1,NX
            S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
  100    CONTINUE
  200 CONTINUE
      RETURN
      END
      SUBROUTINE SETRST2(IVEC)
C     Set up list of candidate RST planes
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3)
      COMMON /PFFLG/  ILGRNG,ISCND,IENTRP,INEWTX
      DIMENSION VEC(3),VEC1(3),VEC2(3),VEC3(3)
      CHARACTER*1 YESNO
      LOGICAL IFANY,IFTMP
C
   10 CALL PRS
     $('Enter a point near RST surface to be plotted.$') 
      CALL RERRR(XPT1,YPT1,ZPT1)
C
C     Find the element closest to this point
C
      CALL FINDE(IE,XPT1,YPT1,ZPT1)
      WRITE(S,15) IE,XCEN(IE),YCEN(IE),ZCEN(IE)
   15 FORMAT(2X,'Found element number',I5,3f9.4,'$')
      CALL PRS(S)
C
C     Exception handling:
      IF (IE.EQ.0) THEN
         WRITE(S,20) XPT1,YPT1,ZPT1
         CALL PRS(S)
   20    FORMAT(2X,'Point (',G10.3,',',G10.3,',',G10.3,') is not in '
     $            ,'the domain.$')
         CALL PRS('  Continue (y/n)?$')
         CALL RES(YESNO,1)
         IF (YESNO.EQ.'Y'.OR.YESNO.EQ.'y') GOTO 10
         RETURN
      ENDIF
C
C     Now that element is found, find what type of plane is closest: R,S, or T?
C
      CALL RZERO(VEC,3)
      IF (IVEC.GE.1.AND.IVEC.LE.3) THEN
         VEC(IVEC)=1.0
      ELSE
         CALL PRS
     $('Enter approximate normal vector of the surface at this point.$')
         CALL PRS
     $('  (x_n,y_n,z_n):$') 
         CALL RERRR(VEC(1),VEC(2),VEC(3))
         CALL NORM3D(VEC)
      ENDIF
C
C     Does this correspond most closely to an R, S, or T plane?
C
      RLMAX=0.0
      DO 30 IPLN=1,6 
         CALL SUB3(VEC1,XYZQ(1,2,IPLN,IE),XYZQ(1,1,IPLN,IE),3)
         CALL SUB3(VEC2,XYZQ(1,4,IPLN,IE),XYZQ(1,3,IPLN,IE),3)
         CALL CROSS(VEC3,VEC1,VEC2)
         CALL NORM3D(VEC3)
         RLNGTH = DOTPROD(VEC,VEC3)
         RLNGTH = ABS(RLNGTH)
         IF (RLNGTH.GT.RLMAX) THEN
            IJKPLN  = IPLN
            RLMAX   = RLNGTH
            DXORD=VEC3(1) 
            DYORD=VEC3(2) 
            DZORD=VEC3(3) 
            write(6,*) ' ipln:',ipln,rlmax
         ENDIF
   30 CONTINUE
      VEC3(1)=DXORD
      VEC3(2)=DYORD
      VEC3(3)=DZORD
C
C     Diagnostics
c     IFTMP=IFHARD
c     IFHARD=.FALSE.
c     CALL COLOR(5)
c     CALL ARROW4S(XPT1,YPT1,ZPT1,DXORD,DYORD,DZORD)
c     CALL COLOR(1)
C
      DXORD=VEC(1)
      DYORD=VEC(2)
      DZORD=VEC(3)
c     CALL ARROW4S(XPT1,YPT1,ZPT1,DXORD,DYORD,DZORD)
      IFHARD=IFTMP
      write(S,40) IE,IJKPLN
      CALL PRS(S)
   40 FORMAT(2X,'IE=',I5,' PLANE=',I1,'.$')
C
C     Find the plane corresponding to the chosen point and normal.
C
      NXY =NX*NY
      NXYZ=NX*NY*NZ
      DMIN=10.0E15
      CALL GENXYZ(XP,YP,ZP,IE,1,NX,NY,NZ)
      DO 100 I=1,NXYZ
         DIST=(XP(I)-XPT1)**2+(YP(I)-YPT1)**2+
     $        (ZP(I)-ZPT1)**2
         IF (DIST.LT.DMIN) THEN
            IJKMIN=I
            DMIN=DIST
         ENDIF
  100 CONTINUE
      CALL DECOD(IXP,IYP,IZP,IDUM,IJKMIN,NX,NY,NZ)
      IF (IJKPLN.LE.2) IPLAN=IXP
      IF (IJKPLN.GT.2) IPLAN=IYP
      IF (IJKPLN.GT.4) IPLAN=IZP
C     Set IJKPLN according to actual value of IPLAN
      NXT=2*IPLAN-1
      IF (NXT.LT.NX) IJKPLN=2*((IJKPLN-1)/2)+1
      IF (NXT.GT.NX) IJKPLN=2*((IJKPLN-1)/2)+2
      write(S,40) IE,IJKPLN
      CALL PRS(S)
C
C===============================================
C     Start generating LIST of plotting planes
C===============================================
C
C     Initialize pointers to zero
      DO 105 I=1,24*NELM
         XYZQ(5,I,1,1)=0.0
  105 CONTINUE
C
      NLSTP=1
      LIST(NLSTP)=6*NX*(IE-1)+NX*(IJKPLN-1)+IPLAN
      IMD=MOD(IJKPLN,2)
C
C     Set the sign of this plane according to requested normal vector:
C
      RLNGTH = DOTPROD(VEC,VEC3)
      INRM=1
      IF (RLNGTH.LT.0.0) INRM=-1
      LIST(NLSTP)=INRM*LIST(NLSTP)
C
      DO 110 I=1,4
         NEIGHI=INT(XYZQ(4,I,IJKPLN,IE))
         IF (NEIGHI.EQ.0) THEN
            XYZQ(5,I,IJKPLN,IE)=2.0
         ELSE
            RN=FLOAT(NEIGHI)
            XYZQ(5,I,IJKPLN,IE)=SIGN(1.0,RN)
         ENDIF
  110 CONTINUE
C
C=======================================
C     Find all adjoining planes
C=======================================
C
      ICOUNT=0
  200 CONTINUE
      IFANY=.FALSE.
      ICOUNT=ICOUNT+1
      DO 500 IE=1,NEL
      DO 500 IPLN=1,6
      DO 500 IQ=1,4
         IF (ABS(XYZQ(5,IQ,IPLN,IE)).EQ.1.0) THEN
            IFANY=.TRUE.
            NEIGHI=INT(XYZQ(4,IQ,IPLN,IE))
            NEIGHA=ABS(NEIGHI)
            IF (XYZQ(5,NEIGHA,1,1).EQ.0.0) THEN
C              The neighboring plane has Not been set.
               CALL DECOD(JQ,JPLN,JE,IDUM,NEIGHA,4,6,NEL)
               DO 400 J=1,3
                  JJ=JQ+J
                  JJ=MOD1(JJ,4)
                  NEIGHJ=INT(XYZQ(4,JJ,JPLN,JE))
                  IF (NEIGHJ.EQ.0) THEN
                     XYZQ(5,JJ,JPLN,JE)=2.0
                  ELSE
                     RN=FLOAT(NEIGHJ)
                     XYZQ(5,JJ,JPLN,JE)=
     $                   XYZQ(5,IQ,IPLN,IE)*SIGN(1.0,RN)

                  ENDIF
  400          CONTINUE
C
C              New neighbor points all set, now add neighbor plane to LIST.
               JMD=MOD(JPLN,2)
               IF (IMD.EQ.JMD) THEN
                  JPLAN=IPLAN
               ELSE
                  JPLAN=NX+1-IPLAN
               ENDIF
               NLSTP=NLSTP+1
               LIST(NLSTP)=6*NX*(JE-1)+NX*(JPLN-1)+JPLAN
               LIST(NLSTP)=LIST(NLSTP)*SIGN(1.0,XYZQ(5,IQ,IPLN,IE))
C              Flip according to initial plane
               LIST(NLSTP)=INRM*LIST(NLSTP)
               write(S,401) NLSTP,JE,JPLN,JPLAN,LIST(NLSTP)
  401          format(' Found',i3,'th plane: JE,JP,PLANE:',4I6,'$')
               CALL PRS(S)
C
C              Remove pointers and exit
               XYZQ(5,NEIGHA,1,1)=2.0
            ENDIF
            XYZQ(5,IQ,IPLN,IE)=2.0
         ENDIF
  500 CONTINUE
      IF (IFANY) GOTO 200
C
C     { LIST is complete, sort by zbuff. (? or perhaps after rotation?)}
C     LIST is complete, set normals according to the following:
C
C     If plane is a Y-plane, i.e. 3 or 4, then it must be flipped (?).
C
      RETURN
      END
      SUBROUTINE FINDE(IE,XPT1,YPT1,ZPT1)
      INCLUDE 'basics.inc'
C
C     Find element which is closest to the point xpt1
C
      DSTMIN=10.0E10
      IF (IF3D) THEN
         DO 100 IEL=1,NEL
            DIST = (XCEN(IEL)-XPT1)**2 
     $           + (YCEN(IEL)-YPT1)**2 
     $           + (ZCEN(IEL)-ZPT1)**2 
            IF (DIST.LT.DSTMIN) THEN
               DSTMIN=DIST
               IEMIN=IEL
            ENDIF
  100    CONTINUE
      ELSE
         DO 200 IEL=1,NEL
            DIST = (XCEN(IEL)-XPT1)**2 
     $           + (YCEN(IEL)-YPT1)**2 
            IF (DIST.LT.DSTMIN) THEN
               DSTMIN=DIST
               IEMIN=IEL
            ENDIF
  200    CONTINUE
      ENDIF
      IE=IEMIN
      RETURN
      END
      SUBROUTINE EVALSC( X0 , SCAL , RRL , INEW )
C
C     Evaluate a scalar, SCAL, at position RRL and return the result in X0.
C
      INCLUDE 'basics.inc'
      DIMENSION SCAL(1)
      DIMENSION RRL(3)
      COMMON  /CTMP0q/ HR(NXM),HS(NXM),HT(NXM)
     $               ,HHH(NXM,NXM,NXM)
C
      REAL    ICALLD
      SAVE    ICALLD
      REAL    RLXOLD,RLYOLD,RLZOLD
      SAVE    RLXOLD,RLYOLD,RLZOLD
      DATA    ICALLD /0/
      DATA    RLXOLD,RLYOLD,RLZOLD/3*0.0/
C
      SUM = 0.0
C
      IF (RRL(1).NE.RLXOLD.OR.RRL(2).NE.RLYOLD.OR.RRL(3).NE.RLZOLD
     $   .OR.ICALLD.EQ.0) THEN
         ICALLD = 1
         RLXOLD = RRL(1)
         RLYOLD = RRL(2)
         RLZOLD = RRL(3)
         CALL HLEGN(HR,RRL(1),ZPTS,NX)
         CALL HLEGN(HS,RRL(2),ZPTS,NY)
         IF (IF3D) THEN
            CALL HLEGN(HT,RRL(3),ZPTS,NZ)
            IJK=0
            DO 100 K=1,NZ
            DO 100 J=1,NY
            DO 100 I=1,NX
               IJK=IJK+1
               HHH(IJK,1,1)=HR(I)*HS(J)*HT(K)
               SUM=SUM+SCAL(IJK)*HHH(IJK,1,1)
  100       CONTINUE
         ELSE
            IJK=0
            DO 200 J=1,NY
            DO 200 I=1,NX
               IJK=IJK+1
               HHH(IJK,1,1)=HR(I)*HS(J)
               SUM=SUM+SCAL(IJK)*HHH(IJK,1,1)
 200       CONTINUE
         ENDIF
      ELSE
C
C     Reevaluate a new scalar at an old point.
C
         NXYZ=NX*NY*NZ
         DO 300 I=1,NXYZ
            SUM=SUM+SCAL(I)*HHH(I,1,1)
  300    CONTINUE
      ENDIF
      X0 = SUM
      RETURN
      END
      SUBROUTINE CROSS(V1,V2,V3)
C
C     Compute Cartesian vector cross product.
C
      DIMENSION V1(3),V2(3),V3(3)
C
      V1(1) = V2(2)*V3(3) - V2(3)*V3(2)
      V1(2) = V2(3)*V3(1) - V2(1)*V3(3)
      V1(3) = V2(1)*V3(2) - V2(2)*V3(1)
C
      RETURN
      END
      SUBROUTINE NORM3D(V1)
C
C     Compute Cartesian vector dot product.
C
      DIMENSION V1(3)
C
      VLNGTH = DOTPROD(V1,V1)
      VLNGTH = SQRT(VLNGTH)
      V1(1) = V1(1) / VLNGTH
      V1(2) = V1(2) / VLNGTH
      V1(3) = V1(3) / VLNGTH
C
      RETURN
      END
      SUBROUTINE RENUM
C
C     Renumber elements so that the lowest numbered elements are
C     in the user specified box.
C
      INCLUDE 'basics.inc'
      DIMENSION IEIN(NELM),IEOUT(NELM) 
      NELM1=NELM-1
C 
C     Find out which element are to be renumbered:
      CALL PRS(
     $    'Enter (with mouse) 2 points in element to be deleted,$')
      CALL PRS(
     $    'or, 2 points framing a box containing elements.$')
      CALL PRS(
     $    'Enter in menu area to abort DELETE ELEMENT operation.$')
      IFTMP =IFGRID
      IFGRID=.FALSE.
  120 CONTINUE
      CALL PRS('Enter 1st point:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF (XMOUSE.GT.XPHY(1.0)) THEN
C        look for a keypad input
         CALL PRS('Aborting renumber operation.$')
         CALL BEEP
         IFGRID=IFTMP 
         RETURN
      ELSE
         CALL PRS('Enter 2nd point:$')
         CALL MOUSE(XMOUS2,YMOUS2,BUTTON)
         IF (XMOUS2.GT.XPHY(1.0)) THEN
          CALL PRS('Aborting renumber operation.$')
          CALL BEEP
          IFGRID=IFTMP 
          RETURN
         ENDIF
      ENDIF
C
C     We successfully inputted 2 points in the build area
C     Generate element centers
      CALL GENCEN
      XMAX=MAX(XMOUSE,XMOUS2)
      XMIN=MIN(XMOUSE,XMOUS2)
      YMAX=MAX(YMOUSE,YMOUS2)
      YMIN=MIN(YMOUSE,YMOUS2)
C
C     Check box to see if it contains any elements, and renumber them.
C
C     Count the number which are in
      NUMIN=0
      CALL DRWBOX(XMIN,YMIN,XMAX,YMAX,1)
      DO 100 IIEL=1,NEL
         IF (XMIN.LE.XCEN(IIEL) .AND.
     $       XCEN(IIEL).LE.XMAX .AND.
     $       YMIN.LE.YCEN(IIEL) .AND.
     $       YCEN(IIEL).LE.YMAX )         NUMIN=NUMIN+1
  100 CONTINUE
      WRITE(S,101) NUMIN,NEL
  101 FORMAT('Renumbering',I5,' out of',I5,' elements.$')
      CALL PRS(S)
      IF (NUMIN.GT.0) THEN
C        renumber the elements which are inside to be less than NUMIN
         J=0
         K=0
         DO 200 IIEL=1,NEL
            IF (XMIN.LE.XCEN(IIEL) .AND.
     $          XCEN(IIEL).LE.XMAX .AND.
     $          YMIN.LE.YCEN(IIEL) .AND.
     $          YCEN(IIEL).LE.YMAX ) THEN
                IF (IIEL.GT.NUMIN) THEN
C               We've got one which is in, but too high.
                   J=J+1
                   IEIN(J)=IIEL
                ENDIF
             ELSE
                IF (IIEL.LT.NUMIN) THEN
C               We've got one which is out, but too low.
                   K=K+1
                   IEOUT(K)=IIEL
                ENDIF
             ENDIF
  200    CONTINUE
      WRITE(S,201) J,K
  201 FORMAT('Found',I5,' elements in, and',I5,' elements out.$')
      CALL PRS(S)
      IF (J.NE.K) RETURN
C
C        Swap the elements
C
         DO 300 I=1,J
            IE=IEIN(I)
            JE=IEOUT(I)
            CALL COPYEL(IE,NELM1)
            CALL COPYEL(JE,IE)
            CALL COPYEL(NELM1,JE)
  300    CONTINUE
C
      ENDIF
      RETURN
      END
