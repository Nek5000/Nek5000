c-----------------------------------------------------------------------
      subroutine octspl
      INCLUDE 'basics.inc'
      COMMON /SPLITT/ ENEW(NELM),IND(NELM)
      DIMENSION XPT(3),VEC(3),LISTE(8)
c
      logical ifoct,ifhex,ifzsp
C
C     This value can be changed later.
C
      NX=15
      NX=MIN(NX,NXM)
      NY=NX
      NZ=NX
      CALL LEGEND(ZPTS,WGHT,NZ)
C
C     Check on the number and numbering of the new elements to be generated.
C
      Noct = 4
      ifoct = .true.
      ifhex = .false.
      ifzsp = .false.
c
      IF (IF3D) then
         CALL PRS('Oct, quad, hex, z, multi-split? (o/q/h/z/m)$')
         CALL RES(ANS,1)
         ifoct = .false.
         ifhex = .false.
         ifzsp = .false.
         IF (ANS.eq.'o'.or.ANS.eq.'O') then
            Noct = 8
            ifoct = .true.
         elseif (ans.eq.'h'.or.ans.eq.'H') then
            Noct = 6
            ifhex = .true.
         elseif (ans.eq.'z'.or.ans.eq.'Z') then
            Noct = 2
            ifzsp = .true.
         elseif (ans.eq.'m'.or.ans.eq.'M') then
            call msplit
            return
         endif
      else
         CALL PRS('quad or multi-split? (q/m)$')
         CALL RES(ANS,1)
         if (ans.eq.'m'.or.ans.eq.'M') then
            call msplit
            return
         endif
      endif
c
      NELN=NEL*Noct
      NELM1=NELM-1
      IF (NELN.GT.NELM1) THEN
C
         CALL PRS(
     $ 'WARNING: Number of elements after OctSplit operation$')
C
         WRITE(S,51) NELN,NELM1
         CALL PRS(S)
         WRITE(S,52) 
         CALL PRS(S)
   51  FORMAT(2X,'(',I5,') would be greater than the allowed maxium ('
     $                ,I5,').$')
   52    FORMAT(2X,'Aborting zipper.$')
         return
      ELSE
C
         CALL PRS(
     $ 'WARNING: Number of elements after OctSplit operation$')
         WRITE(S,61) NELN
         CALL PRS(S)
   61  FORMAT(2X,'will be (',I5
     $          ,'). Are you sure you want to split (y/n)?$')
         CALL RES(ANS,1)
         IF (ANS.ne.'y'.and.ANS.ne.'Y') return
      ENDIF
C
C     Renumber ALL elements, so that low numbered elements will remain
C     low numbered.  This will be achieved by assigning IE+0.1 to the
C     new element numbers, and then sorting the list.
C
      Neln=Nel
      DO 101 IE=1,NEL
         ENEW(IE)=IE
         Ke = 0
         DO 100 Ioct = 2,Noct
            Neln=Neln+1
            Enew(Neln) = float(ie) + 0.1
            Ke = Ke+1
            liste(ke) = Neln
  100    continue
         liste(Noct) = Neln+1
         if (ifoct) then
            call octsplite(ie,liste) 
         elseif (ifzsp) then
            call zsplite(ie,liste) 
         elseif (ifhex) then
            call hsplite(ie,liste) 
         else
            call qsplite(ie,liste) 
         endif
  101 CONTINUE
C
C     Generate new sub-elements as a result of the oct-split action.
C
      Nnew = Neln-Nel
      NEL=NELN
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
      CALL CURCNT
      CALL VERTADJ
C
C     Exit
C
      WRITE(S,300) Nnew,NEL
  300 FORMAT(I5,' new elements generated. NEL=',I5,'$')
      CALL PRS(S)
C
      CALL GENCEN
C
      return
      END
c-----------------------------------------------------------------------
      subroutine wind3d
C
C     Generate a set of corner window refinements in 3D
C
      INCLUDE 'basics.inc'
      DIMENSION ENEW(NELM),IND(NELM)
      DIMENSION XPT(3),VEC(3)
      COMMON /wierd/ LIST1(NELM),LIST2(NELM)
      DIMENSION VEC1(3),VEC2(3),VEC3(3)
      CHARACTER*1 CBCNXT
C
C     This value can be changed later.
C
      NX=15
      NX=MIN(NX,NXM)
      NY=NX
      NZ=NX
      CALL LEGEND(ZPTS,WGHT,NZ)
      CALL SETQUAD 
C
C     Set RST planes will find the set of dividing planes for the zipper.
C
C     First, get the point of interest
      WRITE(S,10)
   10 FORMAT(' Enter a point near the corner to be refined.$')
      CALL GETPT3(XPT,IE,S)
      IF (IE.EQ.0) return
C
C     Second, find the defining plane
      WRITE(S,11)
      CALL PRS(S)
   11 FORMAT(' Enter normal vector defining the surface$')
      CALL RERRR(VEC(1),VEC(2),VEC(3))
C
C     This normal, plus the specified point can be used to define
C     two normals which determine "splitting planes".  The intersection
C     of these splitting planes yields the requisite corner refinement:
C
C        +---------------------+
C        | .                   |      (The specified normal would be in/out
C        |   .                 |       of the page in this figure.)
C        |     .         ^     |
C        |       .       | v2  |
C        |         +-----------| splitting plane 1
C        |         |           |
C        |    <----|           |
C        |     v1  |     X     |          X is the specified point
C        |         |           |
C        +---------------------+
C                  ^
C                  |
C                  splitting plane 2
C
C
C     Therefore, we need to call SETRST2 twice, with vectors
C     v1 and v2 which are orthogonal to the "given" vector, vec.
C
      CALL FNDPLN(IJKPLN,IPLAN,VEC3,VEC,XPT,IE)
      CALL COPY(VEC,VEC3,3)
      IPLN0=IJKPLN
      IE0  =IE
      write(6,*) 'ie0:',ie0,ipln0
C
C     For a given IJKPLN (taking values 1-6), we can find the orthogonal
C     planes by adding 2, and taking mod1(ijkpln+2,6).  NOTE that MOD1 is
C     like mod, except that MOD1(6,6)=6 instead of 0.
C
C     Compute v1:
C
      IPLN=IJKPLN+2
      IPLN=MOD1(IPLN,6)
C
      CALL SUB3(VEC1,XYZQ(1,2,IPLN,IE),XYZQ(1,1,IPLN,IE),3)
      CALL SUB3(VEC2,XYZQ(1,4,IPLN,IE),XYZQ(1,3,IPLN,IE),3)
      CALL CROSS(VEC3,VEC1,VEC2)
C     Ok, now let's find all corresponding planes
      CALL SETRST2(XPT,VEC3)
C     Save the list.
      CALL ICOPY(LIST1,LIST,NLSTP)
      NLSTP1=NLSTP
C
C     v2 is easy, given v1 (VEC3), and VEC:
C
      CALL CROSS(VEC2,VEC3,VEC)
C     Ok, now let's find all corresponding planes for v2
      CALL SETRST2(XPT,VEC2)
      CALL ICOPY(LIST2,LIST,NLSTP)
      NLSTP2=NLSTP
C
C     The elements we want lie in the intersection of LIST1 and LIST2
C
      NLSTP =0
      DO 30 IL2=1,NLSTP2
         LISTA=ABS(LIST2(IL2))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         DO 20 IL1=1,NLSTP1
            LISTB=ABS(LIST1(IL1))
            CALL DECOD(IPLANE,JPLN,JE,IDUM,LISTB,NX,6,NELM)
            IF (IE.EQ.JE) THEN
C              Save information about both planes, as they define the split.
               NLSTP=NLSTP+1
               LIST(NLSTP)=6*6*(IE-1)+6*(IPLN-1)+JPLN
      write(6,*) 'ieijp0:',ie,ipln,jpln,nlstp
               GOTO 21
            ENDIF
   20    CONTINUE
   21    CONTINUE
   30 CONTINUE
C
C     Unfortunately, even this list is still too large, so we
C     restrict the elements to the ones which are indeed adjacent
C     to the set we want:
C
C     Set up E-E boundary conditions
      IFLD=2
      IF (IFFLOW) IFLD=1
c     IF (ifld.lt.10) goto 917
      CALL SETEBC(IFLD)
C
C     Find list of adjacent elements in the direction we want
      NLSTP1=0
      IELCUR=IE0
      NLSTP1=NLSTP1+1
      LIST1(NLSTP1)=IELCUR
      DO 120 IDIR=1,2
         IELCUR=IE0
         IFCCUR=2*((IPLN0-1)/2)+IDIR
         IFCCUR=EFACE(IFCCUR)
         DO 100 I=1,NLSTP
            CBCNXT=CBC(IFCCUR,IELCUR,IFLD)
            write(6,*) 'cbc',cbcnxt,ifccur,ielcur,idir,ifld
            IF (CBCNXT.NE.'E') GOTO 101
C           Else, we found another element
            IFCNXT=INT(BC(2,IFCCUR,IELCUR,IFLD))
            IELCUR=INT(BC(1,IFCCUR,IELCUR,IFLD))
C
C           Find next face
            IFCNXT=EFACE1(IFCNXT)
            IFCCUR=IFCNXT+2*MOD(IFCNXT,2)-1
            IFCCUR=EFACE(IFCCUR)
C
C           Bump up list count
            NLSTP1=NLSTP1+1
            LIST1(NLSTP1)=IELCUR
  100    CONTINUE
  101    CONTINUE
  120 CONTINUE
C
C     Compare E-E
      CALL ICOPY(LIST2,LIST,NLSTP)
      NLSTP2=NLSTP
C
C     Again, the elements we want lie in the intersection of LIST1 and LIST2
C
      NLSTP =0
      DO 230 IL2=1,NLSTP2
         LISTA=ABS(LIST2(IL2))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         DO 220 IL1=1,NLSTP1
            JE=LIST1(IL1)
            IF (IE.EQ.JE) THEN
C              Save this element in the list
               NLSTP=NLSTP+1
               LIST(NLSTP)=LIST2(IL2)
               GOTO 221
            ENDIF
  220    CONTINUE
  221    CONTINUE
  230 CONTINUE
C
C     Check on the number and numbering of the new elements to be generated.
C
  917 continue
      NELN=NEL+2*NLSTP
      NELM1=NELM-1
      IF (NELN.GT.NELM1) THEN
C
         WRITE(S,150) 
         CALL PRS(S)
         WRITE(S,151) NELN,NELM1
         CALL PRS(S)
         WRITE(S,152) 
         CALL PRS(S)
  150  FORMAT(2X,'WARNING: Number of elements after zipper operation$')
  151  FORMAT(2X,'(',I5,') would be greater than the allowed maxium ('
     $                ,I5,').$')
  152    FORMAT(2X,'Aborting corner refine operation.$')
         return
      ENDIF
C
C     Renumber ALL elements, so that low numbered elements will remain
C     low numbered.  This will be achieved by assigning IE+1/3,2/3 to the
C     new element numbers, and then sorting the list.
C
      DO 400 IE=1,NEL
         ENEW(IE)=IE
  400 CONTINUE
      NELN=NEL
      DO 401 I=1,NLSTP
         LISTA=ABS(LIST(I))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         NELN=NELN+1
         ENEW(NELN)=FLOAT(IE)+0.3333
         NELN=NELN+1
         ENEW(NELN)=FLOAT(IE)+0.6666
  401 CONTINUE
C
C     Generate new sub-elements as a result of the zipper action.
C
C     Third, find the ratio to be split
      WRITE(S,501)
      CALL PRS(S)
  501 FORMAT(' Enter ratio (0=small corner element, negative-abort).$')
      CALL RER(S0)
      DO 500 I=1,NLSTP
         NELN=NEL+2*I-1
         LISTA=ABS(LIST(I))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         WRITE(S,502) IE,NELN
         CALL PRS(S)
         CALL WIN3D(NELN,IE,I,S0)
  500 CONTINUE
  502 FORMAT(' Generating window for element',2I5,'.$')
      NEL=NELN+1
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
C
C     Exit
C
      NLSTP2=NLSTP*2
      WRITE(S,900) NLSTP2,NEL
  900 FORMAT(I5,' new elements generated. NEL=',I5,'$')
      CALL PRS(S)
C
      CALL GENCEN
C
      return
      END
c-----------------------------------------------------------------------
      subroutine zip2
      INCLUDE 'basics.inc'
      DIMENSION ENEW(NELM),IND(NELM)
      DIMENSION XPT(3),VEC(3)
C
C     This value can be changed later.
C
      NX=15
      NX=MIN(NX,NXM)
      NY=NX
      NZ=NX
      CALL LEGEND(ZPTS,WGHT,NZ)
      CALL SETQUAD 
C
C     Set RSTP planes will find the set of dividing planes for the zipper.
C
C     First, get the point of interest
      WRITE(S,10)
   10 FORMAT(' Enter a point near the dividing surface.$')
      CALL GETPT3(XPT,IE,S)
      IF (IE.EQ.0) return
C
C     Second, find the defining plane
      WRITE(S,11)
      CALL PRS(S)
      WRITE(S,12)
      CALL PRS(S)
   11 FORMAT(' Enter the components of a normal vector defining the$')
   12 FORMAT(' surface, e.g.,  1 0 0.$')
      CALL RERRR(VEC(1),VEC(2),VEC(3))
C
C     Ok, now let's find all corresponding planes
      CALL SETRST2(XPT,VEC)
C
C     Check on the number and numbering of the new elements to be generated.
C
      NELN=NEL+NLSTP
      NELM1=NELM-1
      IF (NELN.GT.NELM1) THEN
C
         WRITE(S,50) 
         CALL PRS(S)
         WRITE(S,51) NELN,NELM1
         CALL PRS(S)
         WRITE(S,52) 
         CALL PRS(S)
   50  FORMAT(2X,'WARNING: Number of elements after zipper operation$')
   51  FORMAT(2X,'(',I5,') would be greater than the allowed maxium ('
     $                ,I5,').$')
   52    FORMAT(2X,'Aborting zipper.$')
         return
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
      I0 = 0
      DO 200 I=1,NLSTP
         NELN=NEL+I
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         CALL SPLITE(NELN,IE,I,I0,S1)
  200 CONTINUE
      NEL=NELN
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
      CALL CURCNT
C
C     Exit
C
      WRITE(S,300) NLSTP,NEL
  300 FORMAT(I5,' new elements generated. NEL=',I5,'$')
      CALL PRS(S)
C
      CALL GENCEN
C
      return
      END
c-----------------------------------------------------------------------
      subroutine swapel(ind,n,nelm1)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine genxyz (xml,yml,zml,ie,iel,nxl,nyl,nzl)
C
      INCLUDE 'basics.inc'
C
C     Note : CTMPg is used in this format in several subsequent routines
C
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMPg/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
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
      SAVE      NXO
      DATA      NXO /0/
C
C     Initialize geometry arrays
C
      IF (NXL.NE.NXO) THEN
        NXO = NXL
C
C    Preprocessor Corner notation:      Symmetric Corner notation:
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
C
         CALL LEGEND(ZGML,WORK,NXL)
C
         DO 10 IX=1,NXL
            H(IX,1,1)=(1.0-ZGML(IX,1))*0.5
            H(IX,1,2)=(1.0+ZGML(IX,1))*0.5
   10    CONTINUE
         DO 20 IY=1,NYL
            H(IY,2,1)=(1.0-ZGML(IY,1))*0.5
            H(IY,2,2)=(1.0+ZGML(IY,1))*0.5
   20    CONTINUE
         IF (IF3D) THEN
            DO 30 IZ=1,NZL
               H(IZ,3,1)=(1.0-ZGML(IZ,1))*0.5
               H(IZ,3,2)=(1.0+ZGML(IZ,1))*0.5
   30       CONTINUE
         ELSE
            CALL RONE(H(1,3,1),NZL)
            CALL RONE(H(1,3,2),NZL)
         ENDIF
      ENDIF
C
C
      NXYZ = NXL*NYL*NZL
      CALL RZERO(XML(1,1,1,IEL),NXYZ)
      CALL RZERO(YML(1,1,1,IEL),NXYZ)
      CALL RZERO(ZML(1,1,1,IEL),NXYZ)
C
      NDIM2 = 2**NDIM
      DO 50 IX=1,NDIM2
         I=INDX(IX)
         XCB(IX,1,1)=X(IE,I)
         YCB(IX,1,1)=Y(IE,I)
         ZCB(IX,1,1)=Z(IE,I)
c     write(6,*) 'xycb:',xcb(ix,1,1),ycb(ix,1,1),zcb(ix,1,1),ix,ie
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
        IF (CCV.EQ.'p') 
     $     CALL pSPHSRF(XML,YML,ZML,IFACE,IE,IEL,NXL,NYL,NZL,WORK) 
 1000 CONTINUE
C
      DO 2000 ISID=1,8
        CCV = CCURVE(ISID,IE)
        IF (CCV.EQ.'C') 
     $  CALL ARCSRF(XML,YML,ZML,NXL,NYL,NZL,IE,IEL,ISID)
 2000 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine sphsrf(xml,yml,zml,ifce,ie,iel,mx,my,mz,xysrf) 
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
      COMMON /CTMPg/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      COMMON /CTMP0/ XCV(3,2,2),VN1(3),VN2(3)
     $     ,X1(3),X2(3),X3(3),DX(3),XCC(8),YCC(8),ZCC(8),ZHMT(237)
C HMT change to make consisteny w/build2.f
C      COMMON /CTMP0/ XCV(3,2,2),VN1(3),VN2(3)
C     $              ,X1(3),X2(3),X3(3),DX(3),XCC(8),YCC(8),ZCC(8)
C
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
      return
      END
c-----------------------------------------------------------------------
      subroutine edg3d(xysrf,x1,x2,i1,i2,j1,j2,mx,my)
C
C     Generate XYZ vector along an edge of a surface.
C
      INCLUDE 'basics.inc'
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMPg/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
C
      DIMENSION XYSRF(3,MX,MY)
      DIMENSION X1(3),X2(3)
      REAL U1(3),U2(3),UN(3),B(3)
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
      CALL CROSS(UN,X1,X2)
      CALL CROSS( B,UN,X1)
      CALL NORM3D (UN)
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
      return
      END
c-----------------------------------------------------------------------
      real function dot(v1,v2,n)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine crn3d(xcv,xc,yc,zc,curve,iface,ie)
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
     $       ,/,5X,'EXITING in CRN3Dc',3E12.4,I3)
         CALL EXITT
      ELSE
         DO 40 I=1,4
            RADT=XCV(1,I,1)**2+XCV(2,I,1)**2+XCV(3,I,1)**2
            RADT=SQRT(RADT)
            TEST=ABS(RADT-RADIUS)/RADIUS
            IF (TEST.GT.EPS) THEN
             WRITE(6,30) 
     $       RADT,RADIUS,XCV(1,I,1),XCV(2,I,1),XCV(3,I,1)
   30        FORMAT(5X,'ERROR: Element vertex not on requested sphere.'
     $           ,/,5X,'EXITING in CRN3Dd',5F12.7)
             WRITE(6,31) IE,IFACE,XCTR,YCTR,ZCTR
   31        FORMAT(5X,'IE,IF,XYZCTR:',2I4,3F12.7)
             WRITE(6,32) (xc(j),yc(j),zc(j),j=1,8)
   32        FORMAT(3f12.7)
             CALL EXITT
            ENDIF
   40    CONTINUE
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine dsset(mx,my,mz)
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
      IF (MXO.EQ.MX.AND.MYO.EQ.MY.AND.MZO.EQ.MZ) return
C
C     else, proceed....
C
      MXO = MX
      MYO = MY
      MZO = MZ
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
      return
      END
c-----------------------------------------------------------------------
      subroutine rone(x,n)
      DIMENSION X(1)
      DO 10 I=1,N
         X(I)=1.0
   10 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine arcsrf(xml,yml,zml,nxl,nyl,nzl,ie,iel,isid)
      INCLUDE 'basics.inc'
C
C     ....note..... CTMPg is used in this format in several subsequent routines
C
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      COMMON /CTMPg/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
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
         CALL EXITT
      ENDIF
      XS = PT2Y-PT1Y
      YS = PT1X-PT2X
C     Make length Radius
      XYS=SQRT(XS**2+YS**2)
      write(6,977) ie,isid,pt1x,pt2x,pt1y,pt2y,xys,xs,ys
  977 format('arc',i3,i3,7f10.4)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine addtnsr(s,h1,h2,h3,nx,ny,nz)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine setrst2(xpt,vec)
C     Set up list of candidate RST planes
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3)
      COMMON /PFFLG/  ILGRNG,ISCND,IENTRP,INEWTX
      DIMENSION XPT(3),VEC(3)
      DIMENSION VEC1(3),VEC2(3),VEC3(3)
      CHARACTER*1 YESNO
      LOGICAL IFANY,IFTMP
C
C======================================================================
C     Find out the exact plane which is defined by the 
C     incoming point and vector.
C======================================================================
C
C     Find the element closest to this point
C
      CALL FINDE(IE,XPT(1),XPT(2),XPT(3))
C
C     Does this correspond most closely to an R, S, or T plane?
C
      CALL FNDPLN(IJKPLN,IPLAN,VEC3,VEC,XPT,IE)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine finde(ie,xpt1,ypt1,zpt1)
      INCLUDE 'basics.inc'
C
C     Find element which is closest to the point xpt1
C
      CALL GENCEN
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
      return
      END
c-----------------------------------------------------------------------
      subroutine evalsc( x0 , scal , rrl , inew )
C
C     Evaluate a scalar, SCAL, at position RRL and return the result in X0.
C
      INCLUDE 'basics.inc'
      DIMENSION SCAL(1)
      DIMENSION RRL(3)
      COMMON  /qTMP0q/ HR(NXM),HS(NXM),HT(NXM)
     $               ,HHH(NXM,NXM,NXM)
      common /cevals/ zptc(nxm),wptc(nxm)
C
      integer icalld,nxold
      save    icalld,nxold
      data    icalld /0/
      data    nxold  /0/
c
      REAL    RLXOLD,RLYOLD,RLZOLD
      SAVE    RLXOLD,RLYOLD,RLZOLD
      DATA    RLXOLD,RLYOLD,RLZOLD/3*1.e20/
C
      SUM = 0.0
C
      if (nx.ne.nxold) then
         nxold = nx
         call zwgll(zptc,wptc,nx)
         call copy (zpts,zptc,nx)
      endif
c
      IF (RRL(1).NE.RLXOLD.OR.RRL(2).NE.RLYOLD.OR.RRL(3).NE.RLZOLD
     $   .OR.ICALLD.EQ.0) THEN
         ICALLD = 1
         RLXOLD = RRL(1)
         RLYOLD = RRL(2)
         RLZOLD = RRL(3)
         CALL HLEGN(HR,RRL(1),zptc,NX)
         CALL HLEGN(HS,RRL(2),zptc,NY)

c        write(6,*) 'this is rr1:',rrl(1),nx
c        call outmat(zptc,1,nx,'zptc ',nx)
c        call outmat(hr,1,nx,'hr   ',nx)
c        write(6,*) 'this is rr2:',rrl(2),ny
c        call outmat(hs,1,ny,'hs   ',ny)

         IF (IF3D) THEN
            CALL HLEGN(HT,RRL(3),zptc,NZ)
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
c              write(6,1) ijk,sum,scal(ijk),hr(i),hs(j),hhh(ijk,1,1)
c  1           format(i5,5f12.7,' sum 1')
 200       CONTINUE
         ENDIF
      ELSE
C
C     Reevaluate a new scalar at an old point.
C
         NXYZ=NX*NY*NZ
         DO 300 I=1,NXYZ
            SUM=SUM+SCAL(I)*HHH(I,1,1)
c           write(6,*) i,sum,scal(i),hhh(i,1,1),' sum 2'
  300    CONTINUE
      ENDIF
      X0 = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine cross(v1,v2,v3)
C
C     Compute Cartesian vector cross product.
C
      DIMENSION V1(3),V2(3),V3(3)
C
      V1(1) = V2(2)*V3(3) - V2(3)*V3(2)
      V1(2) = V2(3)*V3(1) - V2(1)*V3(3)
      V1(3) = V2(1)*V3(2) - V2(2)*V3(1)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine norm3d(v1)
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
      return
      END
c-----------------------------------------------------------------------
      subroutine renum
C
C     Renumber elements so that the lowest numbered elements are
C     in the user specified box.
C
      INCLUDE 'basics.inc'
      DIMENSION IEIN(NELM),IEOUT(NELM) 
      logical iftmp
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
         return
      ELSE
         CALL PRS('Enter 2nd point:$')
         CALL MOUSE(XMOUS2,YMOUS2,BUTTON)
         IF (XMOUS2.GT.XPHY(1.0)) THEN
          CALL PRS('Aborting renumber operation.$')
          CALL BEEP
          IFGRID=IFTMP 
          return
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
      IF (J.NE.K) return
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
      return
      END
c-----------------------------------------------------------------------
      subroutine getpt3(xpt,ie,text)
C     Set up list of candidate RST planes
      INCLUDE 'basics.inc'
      CHARACTER*1 YESNO
      CHARACTER*80 TEXT
      DIMENSION XPT(3)
C
      
   10 CALL PRS(TEXT)
      CALL RERRR(XPT(1),XPT(2),XPT(3))
C
C     Find the element closest to this point
C
      CALL FINDE(IE,XPT(1),XPT(2),XPT(3))
C
C     Exception handling:
      IF (IE.EQ.0) THEN
         WRITE(S,20) XPT(1),XPT(2),XPT(3)
         CALL PRS(S)
   20    FORMAT(2X,'Point (',G10.3,',',G10.3,',',G10.3,') is not in '
     $            ,'the domain.$')
         CALL PRS('  Continue (y/n)?$')
         CALL RES(YESNO,1)
         IF (YESNO.EQ.'Y'.OR.YESNO.EQ.'y') GOTO 10
      ELSE
         WRITE(S,30) IE,XCEN(IE),YCEN(IE),ZCEN(IE)
   30    FORMAT(2X,'Found element number',I5,3f9.4,'$')
         CALL PRS(S)
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine fndpln(ijkpln,iplan,vec3,vec,xpt,ie)
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      DIMENSION VEC(3) ,VEC3(3),XPT(3)
      DIMENSION VEC1(3),VEC2(3)
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
         ENDIF
   30 CONTINUE
      VEC3(1)=DXORD
      VEC3(2)=DYORD
      VEC3(3)=DZORD
C
      DXORD=VEC(1)
      DYORD=VEC(2)
      DZORD=VEC(3)
C
C     Find the plane corresponding to the chosen point and normal.
C
      NXY =NX*NY
      NXYZ=NX*NY*NZ
      DMIN=10.0E15
      CALL GENXYZ(XP,YP,ZP,IE,1,NX,NY,NZ)
      DO 100 I=1,NXYZ
         DIST=(XP(I)-XPT(1))**2+(YP(I)-XPT(2))**2+
     $        (ZP(I)-XPT(3))**2
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
C
      return
      END
c-----------------------------------------------------------------------
      subroutine win3d(je,ie,ilist,frac)
C
C     A routine which will perform a corner refinement.
C
C     This normal, plus the specified point can be used to define
C     two normals which determine "splitting planes".  The intersection
C     of these splitting planes yields the requisite corner refinement:
C
C        +---------------------+
C        | .                   |
C        |   .                 |
C        |     .       JE      |
C        |       .             |
C        |         +-----------| splitting plane 1
C        |         |           |
C        |   KE    |           |
C        |         |    IE     |
C        |         |           |
C        +---------------------+
C                  ^
C                  |
C                  splitting plane 2
C
C
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      INTEGER   ITRANS(3,6,6)
      CHARACTER CTRANS(3,6,6)
      INTEGER INV(8)
      SAVE    INV,ITRANS,CTRANS
      DATA    INV /1,2,4,3,5,6,8,7/
      DATA    ITRANS /108*0/
      DATA    CTRANS /108*' '/
C
      IF (FRAC.LE.0.0.OR.FRAC.GE.1.0) return
C
c     ITRANS(1,1,4)=Z+1
c     ITRANS(1,1,5)=Z+1;Y+1
c     ITRANS(1,1,6)=X+1
c     ITRANS(1,2,3)=Z-1
c     ITRANS(1,2,4)=Z+2
c     ITRANS(1,2,5)=Z+1;Y+1;X+1
c     ITRANS(1,2,6)=X+1;Z-1
c     ITRANS(1,3,5)=Y+1
c     ITRANS(1,3,6)=Y+1;Z-1
c     ITRANS(1,4,5)=Y-1;Z-1
c     ITRANS(1,4,6)=Y-1

      ITRANS(1,1,4)= +1
      ITRANS(1,1,5)= +1
      ITRANS(2,1,5)= +1
      ITRANS(1,1,6)= +1
      ITRANS(1,2,3)= -1
      ITRANS(1,2,4)= +2
      ITRANS(1,2,5)= -1
      ITRANS(2,2,5)= -1
C
      ITRANS(1,2,6)= +1
      ITRANS(2,2,6)= -1
      ITRANS(1,3,5)= +1
      ITRANS(1,3,6)= +1
      ITRANS(2,3,6)= -1
C
c     ITRANS(1,4,5)= -1
c     ITRANS(2,4,5)= -1
c     ITRANS(1,4,6)= -1
      ITRANS(1,4,5)= +2
      ITRANS(2,4,5)= -1
      ITRANS(1,4,6)= +1
      ITRANS(2,4,6)= -1
C
      CTRANS(1,1,4)='Z'
      CTRANS(1,1,5)='Z'
      CTRANS(2,1,5)='Y'
      CTRANS(1,1,6)='X'
C
      CTRANS(1,2,3)='Z'
      CTRANS(1,2,4)='Z'
      CTRANS(1,2,5)='X'
      CTRANS(2,2,5)='Z'
C
      CTRANS(1,2,6)='X'
      CTRANS(2,2,6)='Z'
C
      CTRANS(1,3,5)='Y'
      CTRANS(1,3,6)='Y'
      CTRANS(2,3,6)='Z'
C
c     CTRANS(1,4,5)='Y'
c     CTRANS(2,4,5)='Z'
c     CTRANS(1,4,6)='Y'
      CTRANS(1,4,5)='X'
      CTRANS(2,4,5)='Y'
      CTRANS(1,4,6)='X'
      CTRANS(2,4,6)='Y'
C
      FRACS=2.0*FRAC-1.0
      LISTA=ABS(LIST(ILIST))
      CALL DECOD(JPLN,IPLN,IDUM1,IDUM2,LISTA,6,6,NELM)
      KE=JE+1
      write(6,*) 'ieijp:',ie,ipln,jpln
C
C     Highlight the element to be split
C
      CALL GENXYZ (XP,YP,ZP,IE,1,NX,NY,NZ)
      CALL HILITE(IE,XP,YP,ZP,5)
C
C     Begin with copying the essential information
C
C     Note COPYEL copies IE ==> JE
      CALL COPYEL(IE,JE)
      CALL COPYEL(IE,KE)
C                                                   +----+
C                                                   |ke /|
C     Make rotation(s) to get into canonical form:  |__+ |
C                                                   |  | |
C                                               ie  +----+  je
      I1=MIN(IPLN,JPLN)
      I2=MAX(IPLN,JPLN)
      DO 10 IROT=1,3
         CALL ROTELM(IE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM(JE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM(KE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
   10 CONTINUE
C     Compute XYZ on high definition mesh
C             (Note, because of this, ROTELM must rotate CURVE data as well.)
      CALL GENXYZ (XP,YP,ZP,IE,1,NX,NY,NZ)
C
C     Generate element IE
C
C     First part, IE,KE are lower r-halves, JE is upper r-half.
      I1=-3
      DO 100 KCRN=-1,1,2
         RRL(3)=FLOAT(KCRN)
         I1=I1+4
         I2=I1+1
         I3=I1+2
         I4=I1+3
C
C
C       Notation:   +-------+      
C                   |     .'|      ^ Y
C                   C---B'  |      |
C                   |   |   |      |
C                   +---A---+      +---->  X
C
C     "A" point
         RRL(1)=FRACS
         RRL(2)=-1.0
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv1:',xval,yval,zval
C
         X(IE,I2)=XVAL
         Y(IE,I2)=YVAL
         Z(IE,I2)=ZVAL
C
         X(JE,I1)=XVAL
         Y(JE,I1)=YVAL
         Z(JE,I1)=ZVAL
C
C     "B" point
         RRL(1)=FRACS
         RRL(2)=FRACS
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv2:',xval,yval,zval
C
         X(IE,I3)=XVAL
         Y(IE,I3)=YVAL
         Z(IE,I3)=ZVAL
C
         X(JE,I4)=XVAL
         Y(JE,I4)=YVAL
         Z(JE,I4)=ZVAL
C
         X(KE,I2)=XVAL
         Y(KE,I2)=YVAL
         Z(KE,I2)=ZVAL
C
C     "C" point
         RRL(1)=-1.0
         RRL(2)=FRACS
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv3:',xval,yval,zval
C
         X(IE,I4)=XVAL
         Y(IE,I4)=YVAL
         Z(IE,I4)=ZVAL
C
         X(KE,I1)=XVAL
         Y(KE,I1)=YVAL
         Z(KE,I1)=ZVAL
C
  100 CONTINUE
C
C     Don't forget to eliminate Curve side data on new faces!
C
      JFAC1=EFACE(1)
      JFAC2=EFACE(2)
      JFAC3=EFACE(3)
      JFAC4=EFACE(4)
C
      CCURVE(JFAC2,IE)=' '
      CCURVE(JFAC4,IE)=' '
      CCURVE(JFAC1,JE)=' '
      CCURVE(JFAC4,JE)=' '
      CCURVE(JFAC2,JE)=' '
      CCURVE(JFAC3,JE)=' '
C
C     Delete periodic bcs .... will do later, MC?
C
C     Undo rotations via inverse sequence
      I1=MIN(IPLN,JPLN)
      I2=MAX(IPLN,JPLN)
      DO 510 IROT=3,1,-1
         CALL ROTELM1(IE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM1(JE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM1(KE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
  510 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine rotelm(je,irot,crot)
      CHARACTER*1 CROT
C
C     If there's no rotation to be done, we return.
      IF (IROT.EQ.0) return
C
C     Else, we construct the necessary rotation as a series of unit
C     rotations about the specified axis.
C
      NROT=IROT+8
      NROT=MOD(NROT,4)
C
      DO 100 I=1,NROT
         CALL ROTEL(JE,CROT)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine rotelm1(je,irot,crot)
      CHARACTER*1 CROT
C
C     Inverse rotation: if there's no rotation to be done, we return.
      IF (IROT.EQ.0) return
C
C     Else, we construct the necessary rotation as a series of unit
C     rotations about the specified axis.
C
      NROT=IROT+8
      NROT=MOD(NROT,4)
      IF (NROT.EQ.3) THEN
         NROT=1
      ELSEIF (NROT.EQ.1) THEN
         NROT=3
      ENDIF
C
      DO 100 I=1,NROT
         CALL ROTEL(JE,CROT)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine rotel(ie,crot)
C
C     Perform a unit rotation about the specified axis.
C
      INCLUDE 'basics.inc'
      CHARACTER*1 CROT,CCTMP
      DIMENSION TMP(6)
      INTEGER FCRN2 (4,6)
      DATA    FCRN2 / 1,4,8,5,
     $                2,3,7,6,
     $                1,5,6,2,
     $                4,8,7,3,
     $                1,2,3,4,
     $                5,6,7,8 /
      INTEGER FFAC (4,3)
      DATA    FFAC / 1,5,3,6,
     $               4,6,2,5,
     $               1,2,3,4/
C
C     'X' axis
C
      IF (CROT.EQ.'X') THEN
         IFC1=1
         IFC2=2
         IFAC=1
      ELSEIF (CROT.EQ.'Y') THEN
         IFC1=3
         IFC2=4
         IFAC=2
      ELSE
         IFC1=5
         IFC2=6
         IFAC=3
      ENDIF
C
      X1=X(IE,FCRN2(4,IFC1))
      Y1=Y(IE,FCRN2(4,IFC1))
      Z1=Z(IE,FCRN2(4,IFC1))
      X2=X(IE,FCRN2(4,IFC2))
      Y2=Y(IE,FCRN2(4,IFC2))
      Z2=Z(IE,FCRN2(4,IFC2))
      CCTMP=CCURVE(FFAC(4,IFAC),IE)
      CALL COPY(TMP,CURVE(1,FFAC(4,IFAC),IE),6)
      DO 10 I=3,1,-1
         I1=I+1
         X(IE,FCRN2(I1,IFC1))=X(IE,FCRN2(I,IFC1))
         Y(IE,FCRN2(I1,IFC1))=Y(IE,FCRN2(I,IFC1))
         Z(IE,FCRN2(I1,IFC1))=Z(IE,FCRN2(I,IFC1))
         X(IE,FCRN2(I1,IFC2))=X(IE,FCRN2(I,IFC2))
         Y(IE,FCRN2(I1,IFC2))=Y(IE,FCRN2(I,IFC2))
         Z(IE,FCRN2(I1,IFC2))=Z(IE,FCRN2(I,IFC2))
C
         CCURVE(FFAC(I1,IFAC),IE)=CCURVE(FFAC(I,IFAC),IE)
         CALL COPY(CURVE(1,FFAC(I1,IFAC),IE)
     $            ,CURVE(1,FFAC(I,IFAC),IE),6)
   10 CONTINUE
      X(IE,FCRN2(1,IFC1))=X1
      Y(IE,FCRN2(1,IFC1))=Y1
      Z(IE,FCRN2(1,IFC1))=Z1
      X(IE,FCRN2(1,IFC2))=X2
      Y(IE,FCRN2(1,IFC2))=Y2
      Z(IE,FCRN2(1,IFC2))=Z2
      CCURVE(FFAC(1,IFAC),IE)=CCTMP
      CALL COPY(CURVE(1,FFAC(1,IFAC),IE),TMP,6)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine hilite(ie,xp,yp,zp,iclr)
C
C     Draw ISOMETRIC view of general 3D object.
C
      INCLUDE 'basics.inc'
      CHARACTER*1 CB1
      INTEGER ICALLD,NCSOLD,NCSEG1,NCSGM1
      INTEGER ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
      SAVE    ICEDG,IEDGFC,ICFACE
      SAVE    ICALLD,NCSOLD,NCSEG1,NCSGM1
      DIMENSION XP(1),YP(1),ZP(1)
C
      COMMON /HLINE/ HIK(NXM,50),XXIS(50),YYIS(50),ZZIS(50),
     $               XK(50),YK(50),ZK(50)
      COMMON /ILINE/ IC(8),NSKIPX(3)
      COMMON /LLINE/ IFEDG(12,2,NELM)
      LOGICAL IFEDG,IFCURV
      INTEGER NXOLD,NYOLD,NZOLD
      SAVE    NXOLD,NYOLD,NZOLD
      DATA    NXOLD,NYOLD,NZOLD /0,0,0/
C
      DATA    ICALLD,NCSOLD,NCSEG1,NCSGM1 /4*0/
      DATA    IEDGFC /  5,7,9,11,  6,8,10,12,
     $                  1,3,9,10,  2,4,11,12,
     $                  1,2,5,6,   3,4,7,8    /
      DATA    ICEDG / 1,2,1,   3,4,1,   5,6,1,   7,8,1,
     $                1,3,2,   2,4,2,   5,7,2,   6,8,2,
     $                1,5,3,   2,6,3,   3,7,3,   4,8,3,
c      -2D-  (NOTE:  Edges default to faces in 2D! Not consistent, but...)
     $                1,3,2,   2,4,2,   1,2,1,   3,4,1 /
      DATA    ICFACE/ 1,3,5,7, 2,4,6,8,
     $                1,2,5,6, 3,4,7,8,
     $                1,2,3,4, 5,6,7,8,
c      -2D-
     $                1,3,0,0, 2,4,0,0,
     $                1,2,0,0, 3,4,0,0  /
C
C
C     Sort the faces according to Z-location
C     for the color fill plotting routines
c     IF (IF3D) CALL SORTZ
C
      NFACES=2*NDIM
      NEDGE=4
      IF (IF3D) NEDGE=12
      IF(ICALLD.EQ.0)THEN
         ICALLD=1
C
C        Set up logicals determining whether edge is to be
C        displayed or not.
C
C        Assign EB's numbering scheme to PF's scheme.
C
         EFACE(1)=4
         EFACE(2)=2
         EFACE(3)=1
         EFACE(4)=3
         EFACE(5)=5
         EFACE(6)=6
C
C        Assign inverse of EB's numbering scheme to PF's scheme.
C
         EFACE1(1)=3
         EFACE1(2)=2
         EFACE1(3)=4
         EFACE1(4)=1
         EFACE1(5)=5
         EFACE1(6)=6
C
      ENDIF
      IF (NX.NE.NXOLD.OR.NY.NE.NYOLD.OR.NZ.NE.NZOLD) THEN
         NXOLD=NX
         NYOLD=NY
         NZOLD=NZ
         NXYZ = NX*NY*NZ
         NSKIPX(1) = 1
         NSKIPX(2) = NX
         NSKIPX(3) = NX*NY
C
C        Fill corners - 1 through 8.
C
         IC(1) =                                 0
         IC(2) =                             (NX-1)
         IC(3) =                 NX*(NY-1)
         IC(4) =                 NX*(NY-1) + (NX-1)
         IC(5) =  NX*NY*(NZ-1)
         IC(6) =  NX*NY*(NZ-1)             + (NX-1)
         IC(7) =  NX*NY*(NZ-1) + NX*(NY-1)
         IC(8) =  NX*NY*(NZ-1) + NX*(NY-1) + (NX-1)
         IF (NCSOLD.NE.NCSEGS) THEN
C           Set up curve line segment evaluation array (exclude endpoints) :
            NCSOLD = NCSEGS
            NCSEG1 = NCSEGS + 1
            NCSGM1 = NCSEGS - 1
            DR1 = 2.0 / FLOAT(NCSEGS)
            DO 400 K=1,NCSGM1
               RR1 = -1.0 + DR1*FLOAT(K)
               CALL HLEGN(HIK(1,K),RR1,ZPTS,NX)
  400       CONTINUE
         ENDIF
      ENDIF
C
C  ==============================================
C     Preliminaries all set, draw the element
C  ==============================================
C
C
      CALL COLOR(ICLR)
      CALL PENW(1)
C
C     Draw element edges (for now, we do the 4X redundant work)
C
         IEOFF=1
C
C        Check if we can take a short cut in plotting the element bdry.
         IFCURV=.FALSE.
         DO 500 IEDG = 1,8
            IF (CCURVE(IEDG,IE).NE.' ') IFCURV=.TRUE.
  500    CONTINUE
C
C        Edges
C
         DO 800 IEDG = 1,NEDGE
            IEDG0 = IEDG + 12
            IF (IF3D) IEDG0 = IEDG
            ICE1 = ICEDG(1,IEDG0)
            ICE2 = ICEDG(2,IEDG0)
            NSKIP = NSKIPX(ICEDG(3,IEDG0))
            IC0 = IEOFF+IC(ICE1)
            IC1 = IEOFF+IC(ICE2)
            XK(     1) = XP(IC0)
            YK(     1) = YP(IC0)
            ZK(     1) = ZP(IC0)
            XK(NCSEG1) = XP(IC1)
            YK(NCSEG1) = YP(IC1)
            ZK(NCSEG1) = ZP(IC1)
            IF (IFCURV) THEN
               CALL EVLINE(XK(2),XP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL EVLINE(YK(2),YP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL EVLINE(ZK(2),ZP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL VDRAW3A(XK,YK,ZK,NCSEG1)
            ELSE
               XK(2) = XP(IC1)
               YK(2) = YP(IC1)
               ZK(2) = ZP(IC1)
               CALL VDRAW3A(XK,YK,ZK,2)
            ENDIF
  800    CONTINUE
 1000 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine evline(uv,uvp,hik,n1,n1d,npt,nskip)
      INCLUDE 'basics.inc'
C
C     Evaluate vector uv <== uvp, @ k=1,npts
C
      REAL UV(NPT),UVP(1)
      REAL HIK(N1D,NPT)
      IF (NPT.LE.N1) THEN
         DO 100 K=1,NPT
            UV(K) = 0.0
            ISKIP = 1 - NSKIP
            DO 100 I=1,N1
               ISKIP = ISKIP + NSKIP
               UV(K) = UV(K) + HIK(I,K)*UVP(ISKIP)
  100    CONTINUE
      ELSE
         DO 200 K=1,NPT
            UV(K)=0.0
  200    CONTINUE
         ISKIP = 1 - NSKIP
         DO 300 I=1,N1
            ISKIP = ISKIP + NSKIP
            DO 300 K=1,NPT
               UV(K) = UV(K) + HIK(I,K)*UVP(ISKIP)
  300    CONTINUE
      ENDIF
      return
      END
c-----------------------------------------------------------------------
      subroutine setebc(ifld)
      INCLUDE 'basics.inc'
C     Find Sides' Midpoints
      CALL MKSIDE
      CALL GENCEN
C     Find Sides' Overlaps (used for finding boundary sides)
C!! Check here for Irrational B.c.'s: unrequited connectivities;
C     Warn if physical boundaries internally.  This check is not exhaustive!!??
C     ZERO OUT OLD ELEMENTAL CONNECTIVITY
      DO 300 IE=1,NEL
      DO 300 ISIDE=1,NSIDES
         IF(  CBC( ISIDE, IE,IFLD).EQ.'E')THEN
              CBC(   ISIDE, IE,IFLD)=' '
              BC (1, ISIDE, IE,IFLD)= 0
              BC (2, ISIDE, IE,IFLD)= 0
              BC (3, ISIDE, IE,IFLD)= 0
         ENDIF
  300 CONTINUE
C
C     Loop over all elements to find element-element interactions
C
      DO 500 IE=1,NEL
       DELTA1 = .001*RCEN(IE)
       DO 500 JE=1,NEL
C
         DIST2=(XCEN(IE)-XCEN(JE))**2+(YCEN(IE)-YCEN(JE))**2
     $        +(ZCEN(IE)-ZCEN(JE))**2
         DISTR=(RCEN(IE)+RCEN(JE))**2
         IF (DIST2.LE.DISTR) THEN
C           There is a potential interaction
C
             DELTA2 = .001*RCEN(JE)
             DELTA  = MIN(DELTA1,DELTA2)
             DO 400 ISIDE=1,NSIDES
             DO 400 JSIDE=1,NSIDES
                IF(IE.EQ.JE.AND.ISIDE.EQ.JSIDE) GOTO 400
                DELTAX = ABS(SIDES(IE,ISIDE,1)-SIDES(JE,JSIDE,1))
                DELTAY = ABS(SIDES(IE,ISIDE,2)-SIDES(JE,JSIDE,2))
                DELTAZ = ABS(SIDES(IE,ISIDE,3)-SIDES(JE,JSIDE,3))
                IF(DELTAX .LT. DELTA .AND.
     $             DELTAY .LT. DELTA .AND.
     $             DELTAZ .LT. DELTA ) THEN
C                  BC Array used to define neighboring elements
C                  For want of better notation, 
C                 'E' means elemental (internal) bdry
C                  1st two reals are element & side of neighbor; 
C                  3rd is orientation
C
C                  Re-do internal connectivity
C                  "Normal" internal side.
                   CBC3=CBC(ISIDE,IE,IFLD)
                   IF(CBC3(3:3).NE.'I'.AND.CBC3(3:3).NE.'i')THEN
C                    Overlapping edges not Internal B.C.'s are elemental b.c.'s
                     CBC(ISIDE,IE,IFLD)='E'
                     CBC(JSIDE,JE,IFLD)='E'
                   ENDIF
                   IORIEN = 0
                   BC (1,ISIDE,IE,IFLD) = JE
                   BC (2,ISIDE,IE,IFLD) = JSIDE
                   BC (3,ISIDE,IE,IFLD) = IORIEN
                   BC (1,JSIDE,JE,IFLD) = IE
                   BC (2,JSIDE,JE,IFLD) = ISIDE
                   BC (3,JSIDE,JE,IFLD) = IORIEN
                ENDIF
  400       CONTINUE
         ENDIF
  500 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine psphsrf(xml,yml,zml,ifce,ie,iel,mx,my,mz,xysrf) 
C
C     5 Aug 1988 19:29:52 
C     Program to generate spherical shell elements for NEKTON
C     input.  Paul F. Fischer
C
C     Updated 7/30/92 to compute prolate spheroids
C
      INCLUDE 'basics.inc'
      PARAMETER (LX1=NXM,LY1=NYM,LZ1=NZM)
      DIMENSION XML(MX,MY,MZ,1),YML(MX,MY,MZ,1),ZML(MX,MY,MZ,1)
      DIMENSION XYSRF(3,MX,MZ)
C
      COMMON /CTMPg/ H(LX1,3,2),XCRVED(LX1),YCRVED(LY1),ZCRVED(LZ1)
     $             , ZGML(LX1,3),WORK(3,LX1,LZ1)
      COMMON /CTMP0/ XCV(3,2,2),VN1(3),VN2(3)
     $     ,X1(3),X2(3),X3(3),DX(3),XCC(8),YCC(8),ZCC(8),ZHMT(237)
C HMT change to make consisteny w/build2.f
C      COMMON /CTMP0/ XCV(3,2,2),VN1(3),VN2(3)
C     $              ,X1(3),X2(3),X3(3),DX(3),XCC(8),YCC(8),ZCC(8)
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
C
C     For now, we support only stretching in X, for lack of better
C     interface... the code is general however.
C
      pRATIO = CURVE(5,IFCE,IE)
      pRATIi = 1.0/pRATIO
      IFACE  = EFACE1(IFCE)
C
C     First, we map to unstretched coordinates:
C
C       ^      +-----------------+              ^      +-----+
C     y |      |                 |           y' |      |     |
C       |      |                 |              |      |     |
C       +--->  +-----------------+              +--->  +-----+
C         x                                      x'
C
C     Generate edge vectors on the sphere RR=1.0, 
C
C     Generate (normalized) corner vectors XCV(1,i,j),
C     in unstretched coordinates:
C
C       ^      +-----------------+              ^      +-----+
C     y |      |                 |           y' |      |     |
C       |      |                 |              |      |     |
C       +--->  +-----------------+              +--->  +-----+
C         x                                      x'
C
      DO 10 I=1,8
         XCC(I)=X(IE,I)*pRATIi
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
         X2(1) = XML(J1,J2,1,IEL)*pRATIi
         X2(2) = YML(J1,J2,1,IEL)
         X2(3) = ZML(J1,J2,1,IEL)
C
         DX(1) = (XYSRF(1,I,1)-X2(1))*pRATIO
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
      return
      END
c-----------------------------------------------------------------------
      subroutine hilites(ie,iside)
C
C     High-light just the side of the element
C
      INCLUDE 'basics.inc'
C
C     Flash element side
C     Hilight side with Blinking lights
      CALL COLOR(2)
      IF(IFBWGKS)CALL COLOR(0)
      IF(ISIDE.EQ.5.OR.ISIDE.EQ.6) THEN
C        Flash mesh sides
         IF(ISIDE.EQ.5)THEN
            IED1=1
            IED2=4
         ELSE IF(ISIDE.EQ.6)THEN
            IED1=1
            IED2=4
         ENDIF
         CALL MOVE(X(IE, IED1),Y(IE, IED1))
         DO 220 IEDGE=IED1,IED2
            CALL DRAWED(IE,IEDGE,1)
220      CONTINUE
      ELSE
C        SIDES 1-4  (SAME AS EDGES IN THIS CASE)
c        write(6,*) 'ie,is:',ie,iside
c        CALL MOVE(X(IE, ISIDE),Y(IE, ISIDE))
         CALL DRAWED(IE,ISIDE,1)
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine shift
      INCLUDE 'basics.inc'
C
    1 CONTINUE
c
      nchoic = 0
      nchoic = nchoic+1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Redraw mesh'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Center'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Shift X'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Shift Y'
      IF (IF3D) THEN
         nchoic = nchoic+1
         ITEM(nchoic)    =             'Shift Z'
      ENDIF
c
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
c
      n  = nel*(2**ndim)
      xmean = glsum(x,8*nelm)/n
      ymean = glsum(y,8*nelm)/n
      zmean = glsum(z,8*nelm)/n
      if (if3d) then
         call prsrrr('XYZ mean:$',xmean,ymean,zmean)
      else
         call prsrr ('XY  mean:$',xmean,ymean)
      endif
c
      IF (CHOICE.EQ.'UP MENU') return
      IF (CHOICE.EQ.'Redraw mesh') then
         call redraw_mesh
      elseif (choice.eq.'Center') then
c
         xsep = glmin(x,8*nelm) - 1.e5
         ysep = glmin(y,8*nelm) - 1.e5
         zsep = glmin(z,8*nelm) - 1.e5
c
         xshift = -xmean
         yshift = -ymean
         zshift = -zmean
c
         call shifter(xshift,xsep,'>',x,'X')
         call shifter(yshift,ysep,'>',y,'Y')
         if (if3d) call shifter(zshift,zsep,'>',z,'Z')
c
      elseif (choice.eq.'Shift X') then
         CALL PRS(
     $   'Input X-location separating shifted section.$')
         CALL RER(Xsep)
c        CALL DRAWLINE(Xsep,Ymax,Xsep,Ymin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired shift section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to shift in X-direction$')
         CALL RER(Xshift)
         CALL Shifter(Xshift,Xsep,ANS,X,'X')
      ELSEIF (CHOICE.EQ.'Shift Y') THEN
         CALL PRS(
     $   'Input Y-location separating shifted section.$')
         CALL RER(Ysep)
c        CALL DRAWLINE(Ysep,Xmax,Ysep,Xmin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired shift section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to shift in Y-direction$')
         CALL RER(Yshift)
         CALL Shifter(Yshift,Ysep,ANS,Y,'Y')
      ELSEIF (CHOICE.EQ.'Shift Z') THEN
         CALL PRS(
     $   'Input Z-location separating shifted section.$')
         CALL RER(Zsep)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired shift section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to shift in Z-direction$')
         CALL RER(Zshift)
         CALL Shifter(Zshift,Zsep,ANS,Z,'Z')
      ENDIF
      GOTO 1
      END
c-----------------------------------------------------------------------
      subroutine shifter(Shift,Sep,DIR,pts,coord)
      INCLUDE 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
C
      Nvts = 4
      IF (IF3D) Nvts=8
      Nkshift = 0
C
      IF (DIR.eq.'>') THEN
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
   10       CONTINUE
c           IF ((IFG.AND.IFL).and.Shift.lt.0.0) THEN ! cmt out; 9/29/05
C
C              If an element straddles the Separator, we have to
C              ensure that the shift operation doesn't "invert" the
C              element
C
c              IF (pmin+shift.le.pmax) THEN
c                 CALL PRS(
c    $           'Error:  Attempt to shrink element to zero length$')
c                 CALL PRS(' Smax     Smin    Shift $')
c                 CALL PRRR(pmax,pmin,Shift)
c                 CALL PRS('Aborting shift operation$')
c                 return
c              ENDIF
c           ENDIF
            DO 20 j=1,Nvts
               IF (pts(k,j).ge.sep) pts(k,j)=pts(k,j)+shift
   20       CONTINUE
            if (if3d.and.pmin.ge.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  100   CONTINUE
      ELSE
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
  110       CONTINUE
c           IF ((IFG.AND.IFL).and.Shift.gt.0.0) THEN
C
C              If an element straddles the Separator, we have to
C              ensure that the shift operation doesn't "invert" the
C              element
C
c              IF (pmax+shift.le.pmin) THEN
c                 CALL PRS(
c    $           'Error:  Attempt to shrink element to zero length$')
c                 CALL PRS(' Smax     Smin    Shift $')
c                 CALL PRRR(pmax,pmin,Shift)
c                 CALL PRS('Aborting shift operation$')
c                 return
c              ENDIF
c           ENDIF
            DO 120 j=1,Nvts
               IF (pts(k,j).le.sep) pts(k,j)=pts(k,j)+shift
  120       CONTINUE
            if (if3d.and.pmax.le.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  200    CONTINUE
      ENDIF
C
      WRITE(S,500) Nkshift
  500 FORMAT(' Shifted',I5,' elements.$')
      CALL PRS(S)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine curcnt
      INCLUDE 'basics.inc'
C
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
      END
c-----------------------------------------------------------------------
      subroutine splite(je,ie,ilist,ipln,frac)
C
C     A routine which will split an element IE and assign the 
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
c
      integer icalld,ieb1,ieb2
      save    icalld,ieb1,ieb2
      data    icalld,ieb1,ieb2 /0,0,0/
C
      IF (FRAC.LE.0.0.OR.FRAC.GE.1.0) return
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
      IF (IPLN.EQ.0) THEN
         LISTA=ABS(LIST(ILIST))
         CALL DECOD(IPLANE,IPLN,IIE,IDUM,LISTA,NX,6,NELM)
      ENDIF
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
C        Don't forget Curve side data on new faces!
         JFAC1=EFACE(1)
         JFAC2=EFACE(2)
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
         JFAC1=EFACE(3)
         JFAC2=EFACE(4)
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
  120       CONTINUE
  220    CONTINUE
         JFAC1=EFACE(5)
         JFAC2=EFACE(6)
      ENDIF
C
C     Curve side fix up
C
      IF (CCURVE(JFAC1,IE).EQ.'s'  .and.
     $    CCURVE(JFAC2,IE).EQ.'s') THEN
          R = 0.5 * ( CURVE(4,JFAC1,IE) + CURVE(4,JFAC2,IE) )
          CCURVE(JFAC1,JE)='s'
          CCURVE(JFAC2,IE)='s'
          CURVE(4,JFAC2,IE) = R
          CURVE(4,JFAC1,JE) = R
          CALL COPY(CURVE(1,JFAC1,JE),CURVE(1,JFAC1,JE),3)
          CALL COPY(CURVE(1,JFAC2,IE),CURVE(1,JFAC1,IE),3)
      ELSEIF (CCURVE(JFAC1,IE).EQ.'s'  .or.
     $        CCURVE(JFAC2,IE).EQ.'s') THEN
          R1i = 0.0
          R2i = 0.0
          IF (CCURVE(JFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,JFAC1,IE)
          IF (CCURVE(JFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,JFAC2,IE)
          R  = 2.0/(R1i+R2i)
C
c         CCURVE(JFAC1,JE)='s'
c         CCURVE(JFAC2,IE)='s'
          CCURVE(JFAC1,JE)=' '
          CCURVE(JFAC2,IE)=' '
          CURVE(4,JFAC2,IE) = R
          CURVE(4,JFAC1,JE) = R
          CALL COPY(CURVE(1,JFAC1,JE),CURVE(1,JFAC1,JE),3)
          CALL COPY(CURVE(1,JFAC2,IE),CURVE(1,JFAC1,IE),3)
      ENDIF
C
C     Check for cylinder  (Note this scheme fails for convex-convex case!)
C                         pff 4/10/93
C
      IF (IPLN.LE.4) THEN
         DO 300 Ilev=1,2
           IF (CCURVE(JFAC1,IE).EQ.'C'.and.
     $         CCURVE(JFAC2,IE).EQ.'C') THEN
               R1 = CURVE(1,JFAC1,IE)
               R2 = CURVE(1,JFAC2,IE)
               R = 0.5*(R1-R2)
               CURVE(1,JFAC2,IE) = -R
               CURVE(1,JFAC1,JE) =  R
               CCURVE(JFAC1,JE)='C'
               CCURVE(JFAC2,IE)='C'
           ELSEIF (CCURVE(JFAC1,IE).EQ.'C'.  or.
     $             CCURVE(JFAC2,IE).EQ.'C') THEN
               R1i = 0.0
               R2i = 0.0
               IF (CCURVE(JFAC1,IE).EQ.'C') R1i = 1.0/CURVE(1,JFAC1,IE)
               IF (CCURVE(JFAC2,IE).EQ.'C') R2i = 1.0/CURVE(1,JFAC2,IE)
               R  = 2.0/(R1i-R2i)
               CURVE(1,JFAC2,IE) = -R
               CURVE(1,JFAC1,JE) =  R
               CCURVE(JFAC1,JE)=' '
               CCURVE(JFAC2,IE)=' '
           ENDIF
           JFAC1=JFAC1+4
           JFAC2=JFAC2+4
  300    CONTINUE
      ENDIF
C
C     Delete periodic bcs and undo curved sides
C
      return
      END
c-----------------------------------------------------------------------
      subroutine octsplite(ie,liste)
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
     $              ,XVAL(-1:1,-1:1,-1:1),YVAL(-1:1,-1:1,-1:1)
     $              ,ZVAL(-1:1,-1:1,-1:1)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
      DIMENSION LISTE(8)
C
C     Compute XYZ on high definition mesh
C
      CALL GENXYZ (XP,YP,ZP,IE,1,NX,NY,NZ)
C
C     Evaluate x,y,z on 27 pt stencil
C
      if (if3d) then
c
      DO 160 It= -1,1,1
        rrl(3)=float(it)
        DO 160 Is= -1,1,1
          rrl(2)=float(is)
          DO 160 Ir= -1,1,1
            rrl(1)=float(ir)
            CALL EVALSC(XVAL(ir,is,it),XP,RRL,1)
            CALL EVALSC(YVAL(ir,is,it),YP,RRL,0)
            CALL EVALSC(ZVAL(ir,is,it),ZP,RRL,0)
  160 CONTINUE
C
C     Eight quadrants
C
      Ke = 0
      DO 860 Iz=-1,0,1
        DO 840 Iy=-1,0,1
          DO 820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          Ke = Ke+1
          Je = LISTE(Ke)
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
          DO 760 Jz=0,1
            jz0 = Jz+Iz
            DO 760 Jy=0,1
              jy0 = Jy+Iy
              DO 760 Jx=0,1
                jx0 = Jx+Ix
                iv = iv+1
                X(JE,INV(Iv))=XVAL(jx0,jy0,jz0)
                Y(JE,INV(Iv))=YVAL(jx0,jy0,jz0)
                Z(JE,INV(Iv))=ZVAL(jx0,jy0,jz0)
  760     CONTINUE
C
          DO 400 Ifce=1,3
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           IF (Ifce.eq.1) THEN
              JFAC = EFACE(1-Ix)
           ELSEIF (Ifce.eq.2) THEN
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           IF (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') THEN
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           ELSEIF (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') THEN
             R1i = 0.0
             R2i = 0.0
             IF (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             IF (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           IF (Ifce.LE.2) THEN
              DO 300 Ilev=1,2
                IF (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  IF (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  IF (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
                ELSEIF (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1i = 0.0
                  R2i = 0.0
                  IF (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  IF (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
                ENDIF
                JFAC =JFAC +4
                IFAC1=IFAC1+4
                IFAC2=IFAC2+4
  300         CONTINUE
           ENDIF
  400     CONTINUE
  820     CONTINUE
  840   CONTINUE
  860 CONTINUE
C
C     Copy last element generated back onto IE
C
      CALL COPYEL(jE,iE)
c
      else
c
        rrl(3)=0.0
        DO 2160 Is= -1,1,1
          rrl(2)=float(is)
          DO 2160 Ir= -1,1,1
            rrl(1)=float(ir)
            CALL EVALSC(XVAL(ir,is,1),XP,RRL,1)
            CALL EVALSC(YVAL(ir,is,1),YP,RRL,0)
 2160 CONTINUE
C
C     Four quadrants
C
      Ke = 0
        DO 2840 Iy=-1,0,1
          DO 2820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          Ke = Ke+1
          Je = LISTE(Ke)
C
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
            DO 2760 Jy=0,1
              jy0 = Jy+Iy
              DO 2760 Jx=0,1
                jx0 = Jx+Ix
                iv = iv+1
                X(JE,INV(Iv))=XVAL(jx0,jy0,1)
                Y(JE,INV(Iv))=YVAL(jx0,jy0,1)
                Z(JE,INV(Iv))=ZVAL(jx0,jy0,1)
 2760     CONTINUE
C
          DO 2400 Ifce=1,2
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           IF (Ifce.eq.1) THEN
              JFAC = EFACE(1-Ix)
           ELSEIF (Ifce.eq.2) THEN
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           IF (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') THEN
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           ELSEIF (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') THEN
             R1i = 0.0
             R2i = 0.0
             IF (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             IF (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           IF (Ifce.LE.2) THEN
              ilev = 1
              IF (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  IF (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  IF (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
              ELSEIF (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1i = 0.0
                  R2i = 0.0
                  IF (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  IF (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
              ENDIF
           ENDIF
 2400     CONTINUE
 2820     CONTINUE
 2840   CONTINUE
 2860 CONTINUE
C
C     Copy last element generated back onto IE
C
      CALL COPYEL(jE,iE)
c
      endif
C
      return
      END
c-----------------------------------------------------------------------
      subroutine qsplite(ie,liste)
C
C     This routine is a hack-copy of octsplite.     (pff 12/21/98)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      INCLUDE 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
     $              ,XVAL(-1:1,-1:1,-1:1),YVAL(-1:1,-1:1,-1:1)
     $              ,ZVAL(-1:1,-1:1,-1:1)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
      DIMENSION LISTE(8)
C
C     Compute XYZ on high definition mesh
C
      CALL GENXYZ (XP,YP,ZP,IE,1,NX,NY,NZ)
C
C     Evaluate x,y,z on 27 pt stencil
C
      if (if3d) then
c
      DO 160 It= -1,1,1
        rrl(3)=float(it)
        DO 160 Is= -1,1,1
          rrl(2)=float(is)
          DO 160 Ir= -1,1,1
            rrl(1)=float(ir)
            CALL EVALSC(XVAL(ir,is,it),XP,RRL,1)
            CALL EVALSC(YVAL(ir,is,it),YP,RRL,0)
            CALL EVALSC(ZVAL(ir,is,it),ZP,RRL,0)
  160 CONTINUE
C
C     Eight quadrants
C
      Ke = 0
c     DO 860 Iz=-1,0,1
      DO 860 Iz= -1,-1,1
        DO 840 Iy=-1,0,1
          DO 820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          Ke = Ke+1
          Je = LISTE(Ke)
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
c         DO 760 Jz=0,1
          DO 760 Jz=0,2,2
            jz0 = Jz+Iz
            DO 760 Jy=0,1
              jy0 = Jy+Iy
              DO 760 Jx=0,1
                jx0 = Jx+Ix
                iv = iv+1
                X(JE,INV(Iv))=XVAL(jx0,jy0,jz0)
                Y(JE,INV(Iv))=YVAL(jx0,jy0,jz0)
                Z(JE,INV(Iv))=ZVAL(jx0,jy0,jz0)
  760     CONTINUE
C
          DO 400 Ifce=1,2
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           IF (Ifce.eq.1) THEN
              JFAC = EFACE(1-Ix)
           ELSEIF (Ifce.eq.2) THEN
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           IF (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') THEN
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           ELSEIF (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') THEN
             R1i = 0.0
             R2i = 0.0
             IF (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             IF (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           IF (Ifce.LE.2) THEN
              DO 300 Ilev=1,2
                IF (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  IF (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  IF (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
                ELSEIF (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1i = 0.0
                  R2i = 0.0
                  IF (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  IF (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
                ENDIF
                JFAC =JFAC +4
                IFAC1=IFAC1+4
                IFAC2=IFAC2+4
  300         CONTINUE
           ENDIF
  400     CONTINUE
  820     CONTINUE
  840   CONTINUE
  860 CONTINUE
C
C     Copy last element generated back onto IE
C
      CALL COPYEL(jE,iE)
c
      else
c
        rrl(3)=0.0
        DO 2160 Is= -1,1,1
          rrl(2)=float(is)
          DO 2160 Ir= -1,1,1
            rrl(1)=float(ir)
            CALL EVALSC(XVAL(ir,is,1),XP,RRL,1)
            CALL EVALSC(YVAL(ir,is,1),YP,RRL,0)
 2160 CONTINUE
C
C     Four quadrants
C
      Ke = 0
        DO 2840 Iy=-1,0,1
          DO 2820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          Ke = Ke+1
          Je = LISTE(Ke)
C
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
            DO 2760 Jy=0,1
              jy0 = Jy+Iy
              DO 2760 Jx=0,1
                jx0 = Jx+Ix
                iv = iv+1
                X(JE,INV(Iv))=XVAL(jx0,jy0,1)
                Y(JE,INV(Iv))=YVAL(jx0,jy0,1)
                Z(JE,INV(Iv))=ZVAL(jx0,jy0,1)
 2760     CONTINUE
C
          DO 2400 Ifce=1,2
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           IF (Ifce.eq.1) THEN
              JFAC = EFACE(1-Ix)
           ELSEIF (Ifce.eq.2) THEN
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           IF (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') THEN
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           ELSEIF (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') THEN
             R1i = 0.0
             R2i = 0.0
             IF (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             IF (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           IF (Ifce.LE.2) THEN
              ilev = 1
              IF (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  IF (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  IF (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
              ELSEIF (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') THEN
                  R1i = 0.0
                  R2i = 0.0
                  IF (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  IF (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
              ENDIF
           ENDIF
 2400     CONTINUE
 2820     CONTINUE
 2840   CONTINUE
 2860 CONTINUE
C
C     Copy last element generated back onto IE
C
      CALL COPYEL(jE,iE)
c
      endif
C
      return
      END
c-----------------------------------------------------------------------
      subroutine shift2
c
c     This routine is like "shift" except that the shifted points
c     are moved a scaled distance, proportional to the distance from
c     the selected point.
c
c
      INCLUDE 'basics.inc'
C
    1 CONTINUE
c
      nchoic = 0
      nchoic = nchoic+1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Redraw mesh'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch X'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch Y'
      IF (IF3D) THEN
         nchoic = nchoic+1
         ITEM(nchoic)    =             'Stretch Z'
      ENDIF
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch R'
c
C     Menu's all set, prompt for user input:
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER')
c
C
      IF (CHOICE.EQ.'UP MENU') return
      IF (CHOICE.EQ.'Redraw mesh') then
         call redraw_mesh
      ELSEIF (CHOICE.EQ.'Stretch R') THEN
         call stretch_rad
      ELSEIF (CHOICE.EQ.'Stretch X') THEN
         CALL PRS(
     $   'Input X-location separating shifted section.$')
         CALL RER(Xsep)
c        CALL DRAWLINE(Xsep,Ymax,Xsep,Ymin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in X-direction$')
         CALL RER(Xshift)
         CALL Shifter2(Xshift,Xsep,ANS,X,'X')
      ELSEIF (CHOICE.EQ.'Stretch Y') THEN
         CALL PRS(
     $   'Input Y-location separating shifted section.$')
         CALL RER(Ysep)
c        CALL DRAWLINE(Ysep,Xmax,Ysep,Xmin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in Y-direction$')
         CALL RER(Yshift)
         CALL Shifter2(Yshift,Ysep,ANS,Y,'Y')
      ELSEIF (CHOICE.EQ.'Stretch Z') THEN
         CALL PRS(
     $   'Input Z-location separating shifted section.$')
         CALL RER(Zsep)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         IF (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in Z-direction$')
         CALL RER(Zshift)
         CALL Shifter2(Zshift,Zsep,ANS,Z,'Z')
      ENDIF
      GOTO 1
      END
c-----------------------------------------------------------------------
      subroutine redraw_mesh
      include 'basics.inc'
c
      CALL REFRESH
      CALL DRMENU('NOCOVER')
      CALL DRGRID
      DO 170 IEL=1,NEL
         CALL DRAWEL(IEL)
 170  CONTINUE
c     
C     Now redraw all the isometric elements.
      call sortel
      do i=1,nel
         call drawis(isrt(i))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine shifter2a(pmin,pmax,gain,Shift,Sep,DIR,pts,coord)
c
c     In this pass, we just figure out range of geometry in shifted 
c     section
c
      INCLUDE 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
C
      CALL PRS('Input Density Gain, G :$')
      CALL RER(gain)
      if (gain.eq.0) gain=1.
      if (gain.lt.0) gain=-1./gain
c
      Nvts = 4
      IF (IF3D) Nvts=8
      Nkshift = 0
C
      pmin= 9.99e15
      pmax=-9.99e15
c
      DO K=1,NEL
      DO j=1,Nvts
         pmin=min(pmin,pts(k,j))
         pmax=max(pmax,pts(k,j))
      enddo
      enddo
c
      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2b(Shift,Sep,DIR,pts,coord)
c
c     Standard linear stretch -- unity gain.
c
      INCLUDE 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
c
      Nvts = 4
      IF (IF3D) Nvts=8
      Nkshift = 0
C
      IF (DIR.eq.'>') THEN
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
   10       CONTINUE
            DO 20 j=1,Nvts
               IF (pts(k,j).ge.sep) pts(k,j)=shift*(pts(k,j)-sep) + sep
   20       CONTINUE
            if (if3d.and.pmin.ge.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
c                 if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
c                 if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
c                 if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  100   CONTINUE
      ELSE
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
  110       CONTINUE
            DO 120 j=1,Nvts
c              IF (pts(k,j).le.sep) pts(k,j)=shift*(pts(k,j)-sep) + sep
               IF (pts(k,j).le.sep) then
                   pold = pts(k,j)
                   pts(k,j)=shift*(pts(k,j)-sep) + sep
                   pnew = pts(k,j)
                   write(6,*) 'stretch',pold,pnew,sep,shift,j,dir,coord
               endif
  120       CONTINUE
            if (if3d.and.pmax.le.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
c                 if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
c                 if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
c                 if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  200    CONTINUE
      ENDIF
C
      WRITE(S,500) Nkshift
  500 FORMAT(' Shifted',I5,' elements.$')
      CALL PRS(S)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2c(Shift,Sep,DIR,pts,coord,gain,qmax,qmin)
c
c     Geometric gain
c
      INCLUDE 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
C
c
      Nvts = 4
      IF (IF3D) Nvts=8
      Nkshift = 0
c
c     We only get here if gain > 0,  different than 1
c
      g = log(gain)
      c = shift*g/(gain-1.)
c
      IF (DIR.eq.'>') THEN
         dx  = qmax - sep
         gdx = g/dx
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
   10       CONTINUE
            DO 20 j=1,Nvts
               if (pts(k,j).ge.sep) then
                  argu = gdx*(pts(k,j)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  write(6,1) shift,gain,g,c,dx,gdx,pts(k,j),xn
    1             format('s:',1p8e12.3)
                  pts(k,j)=xn
               endif
   20       CONTINUE
            if (if3d.and.pmin.ge.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
c                 if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
c                 if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
c                 if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  100   CONTINUE
      ELSE
         dx  = qmin - sep
         gdx = g/dx
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,Nvts
               IF (pts(k,j).ge.sep) THEN
                  IFG=.TRUE.
                  pmin=min(pmin,pts(k,j))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(k,j))
               ENDIF
  110       CONTINUE
            DO 120 j=1,Nvts
               if (pts(k,j).le.sep) then
                  argu = gdx*(pts(k,j)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  pts(k,j)=xn
               endif
  120       CONTINUE
            if (if3d.and.pmax.le.sep) then
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
c                 if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
c                 if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
c                 if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            endif
            Nkshift=Nkshift+1
  200    CONTINUE
      ENDIF
C
      WRITE(S,500) Nkshift
  500 FORMAT(' Shifted',I5,' elements.$')
      CALL PRS(S)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2(Shift,Sep,DIR,pts,coord)
c
c     Standard linear stretch -- unity gain.
c
      INCLUDE 'basics.inc'
      DIMENSION pts(nelm,8)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
c
      call shifter2a(qmin,qmax,gain,Shift,Sep,DIR,pts,coord)
      write(6,*) 'qmnx gn:',qmin,qmax,gain
c
      if (abs(gain-1.) .gt. 1.e-5) then
         call shifter2c(Shift,Sep,DIR,pts,coord,gain,qmax,qmin)
      else
         call shifter2b(Shift,Sep,DIR,pts,coord)
      endif
c
      return
      END
c-----------------------------------------------------------------------
      subroutine hsplite(ie,liste)
C
C     This routine is a hack-copy of octsplite to generate hex
c     decompositions of the square     (pff 8/10/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
      integer inv(8)
      save    inv
      data    inv /1,2,4,3,5,6,8,7/
      dimension liste(8)
c
c     Evaluate x,y,z on new stencil
c
c
c     3---a-------4
c     |   |     / |
c     a---a   /   |
c     |    \0\    | 
c     |   /   a---+
c     | /     |   |
c     1-------+---2
c
c
      integer qmap(4,6) ! quad map
      save    qmap
      data    qmap /   1,  2,  6,  7   ! These orderings will
     $             ,   2,  3,  7,  8   ! preserve the existing BCs
     $             ,   7,  8,  6, 11   !
     $             ,   1,  6,  4,  5   ! Except for periodic, of course
     $             ,   4,  5,  9, 10   !
     $             ,   5,  6, 10, 11 /
c
c     9---0-------1
c     |   |     / |
c     4---5   /   |
c     |    \6\    | 
c     |   /   7---8
c     | /     |   |
c     1-------2---3
c
c
      parameter (aaa = .20)
      parameter (aam = -aaa)
      real rs(2,11)
      save rs
      data rs /  -1.0 , -1.0   !  1
     $        ,   aaa , -1.0   !  2
     $        ,   1.0 , -1.0   !  3
     $        ,  -1.0 ,  aaa   !  4
     $        ,   aam ,  aaa   !  5
     $        ,   0.0 ,  0.0   !  6
     $        ,   aaa ,  aam   !  7
     $        ,   1.0 ,  aam   !  8
     $        ,  -1.0 ,  1.0   !  9
     $        ,   aam ,  1.0   ! 10
     $        ,   1.0 ,  1.0 / ! 11
      real xval(3,11,2)
      integer e,v,q
c
c     Compute XYZ on high definition mesh
c
      call genxyz (xp,yp,zp,ie,1,nx,ny,nz)
c
      iz1 = 1
      if (if3d) iz1 = 2
c
      do k=1,iz1
      do v=1,11
         rrl(1) = rs(1,v)
         rrl(2) = rs(2,v)
         rrl(3) = 2*k-3
         call evalsc(xval(1,v,k),xp,rrl,1)
         call evalsc(xval(2,v,k),yp,rrl,0)
         call evalsc(xval(3,v,k),zp,rrl,0)
      enddo
      enddo
c
C
C     six elements
C
      Ke = 0
      do e=1,6
C
C        Begin by copying the essential information
C
         Ke = Ke+1
         Je = LISTE(Ke)
         call copyel(ie,je)
C
C        Eight vertices
C
         v = 0
         do kv=1,2
         do jv=1,4
            q = qmap(jv,e)
            v = v+1
            x(je,inv(v))=xval(1,q,kv)
            y(je,inv(v))=xval(2,q,kv)
            z(je,inv(v))=xval(3,q,kv)
         enddo
         enddo
C
      enddo
C
C     Copy last element generated back onto IE
C
      call copyel(Je,Ie)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine zsplite(ie,liste)
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
      integer inv(8)
      save    inv
      data    inv /1,2,4,3,5,6,8,7/
      dimension liste(8)
c
      real xval(3,12)
      integer e,v,q
c
c     Compute XYZ on high definition mesh
c
      call prsiii('This is nxyz:$',nx,ny,nz)
      call genxyz (xp,yp,zp,ie,1,nx,ny,nz)
c
      l=0
      do k=-1,1
      do j=-1,1,2
      do i=-1,1,2
         l = l+1
         rrl(1) = i
         rrl(2) = j
         rrl(3) = k
         call evalsc(xval(1,l),xp,rrl,1)
         call evalsc(xval(2,l),yp,rrl,0)
         call evalsc(xval(3,l),zp,rrl,0)
      enddo
      enddo
      enddo
c
c
c     two elements
c
      Ke = 0
      kv = 0
      do e=1,2
C
C        Begin by copying the essential information
C
         Ke = Ke+1
         Je = LISTE(Ke)
         call copyel(ie,je)
C
C        Eight vertices
C
         do v=1,8
            x(je,inv(v))=xval(1,v+kv)
            y(je,inv(v))=xval(2,v+kv)
            z(je,inv(v))=xval(3,v+kv)
         enddo
         kv = 4
C
      enddo
C
C     Copy last element generated back onto IE
C
      call copyel(Je,Ie)
c
      return
      end
c-----------------------------------------------------------------------
      function  ie_click(prompt)  ! element that is clicked upon
      include 'basics.inc'
      character*80 prompt
C
      call prs(prompt)
      call mouse(xmouse,ymouse,button)
      IF(XSCR(XMOUSE).GT.1.0 .AND. YSCR(YMOUSE).GT.0.62) THEN
C        apparently is trying to use the keypad
         CALL PRS('** Entering Coordinates of element. **$')
         CALL PRS('Enter X-coordinate with keypad:$')
         CALL KEYPAD(XMOUSE)
         CALL PRS('Now enter Y-coordinate with keypad:$')
         CALL KEYPAD(YMOUSE)
      ENDIF
      RMIN=1.0E10
      DO 100 IEL=1,NEL
         RAD=SQRT( (XMOUSE-XCEN(IEL))**2 + (YMOUSE-YCEN(IEL))**2 )
         IF(RAD.LT.RMIN .AND. NUMAPT(IEL).EQ.ILEVEL)THEN
            RMIN=RAD
            ie_click=IEL
         ENDIF
100   CONTINUE
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msplit
C
C     This routine is a hack-copy of octsplite to generate 
c     multi-element decompositions of the square     (pff 9/28/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
      include 'basics.inc'
c
      ie = ie_click('Click on element to refine:$')
      if (if3d) then
         call prs('Enter number of partitions in r,s,t (>0):')
         call reiii(nxsp,nysp,nzsp)
      else
         call prs('Enter number of partitions in r,s (>0):')
         call reii(nxsp,nysp)
         nzsp = 1
      endif
      call msplite(ie,nxsp,nysp,nzsp)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msplite(ie,nxsp,nysp,nzsp)
C
C     This routine is a hack-copy of octsplite to generate 
c     multi-element decompositions of the square     (pff 9/28/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),
     $               rrl(3),xval(3,0:1,0:1,0:1)
      integer inv(8)
      save    inv
      data    inv /1,2,4,3,5,6,8,7/
      integer e,v
c
      call genxyz (xp,yp,zp,ie,1,nx,ny,nz) ! high definition mesh
c
      dzt = 2./nzsp
      dyt = 2./nysp
      dxt = 2./nxsp
c
      k1 = 0
      if (if3d) k1=1
c
      e  = 0
      do k=1,nzsp
         rz = -1. + (k-1)*dzt
         do j=1,nysp
            ry = -1. + (j-1)*dyt
c
            do kk=0,k1
            do jj=0,1
               rrl(1) = -1
               rrl(2) =  ry + jj*dyt
               rrl(3) =  rz + kk*dzt
               call evalsc(xval(1,1,jj,kk),xp,rrl,1)
               call evalsc(xval(2,1,jj,kk),yp,rrl,1)
               if (if3d) call evalsc(xval(3,1,jj,kk),zp,rrl,1)
            enddo
            enddo
c
            do i=1,nxsp
c
               do kk=0,k1
               do jj=0,1
                  xval(1,0,jj,kk) = xval(1,1,jj,kk)
                  xval(2,0,jj,kk) = xval(2,1,jj,kk)
                  xval(3,0,jj,kk) = xval(3,1,jj,kk)
                  rrl(1) = -1  + i*dxt
                  rrl(2) =  ry + jj*dyt
                  rrl(3) =  rz + kk*dzt
                  call evalsc(xval(1,1,jj,kk),xp,rrl,1)
                  call evalsc(xval(2,1,jj,kk),yp,rrl,1)
                  if (if3d) call evalsc(xval(3,1,jj,kk),zp,rrl,1)
               enddo
               enddo
c
               e  = e+1
               je = nel+e
               call copyel(ie,je)  ! copy ie to je
C
C              Eight vertices
C
               v = 0
               do kk=0,k1
               do jj=0,1
               do ii=0,1
                  v = v+1
                  x(je,inv(v))=xval(1,ii,jj,kk)
                  y(je,inv(v))=xval(2,ii,jj,kk)
                  z(je,inv(v))=xval(3,ii,jj,kk)
         write(6,9) ii,jj,kk,v,inv(v),je,x(je,inv(v)),y(je,inv(v))
    9    format(5i4,i7,1p2e13.5,' XY')
               enddo
               enddo
               enddo

               call fix_curve(je,i,j,k,nxsp,nysp,nzsp)

               call drawel(je)
 
            enddo
         enddo
      enddo
C
C     Copy last element generated back onto IE
C
      call copyel(Je,Ie)
      nel = nel + (e-1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fix_curve(e,ii,jj,kk,nxsp,nysp,nzsp)
      include 'basics.inc'

      integer c_pair(2,4)
      save    c_pair
      data    c_pair / 4,2 , 8,6 , 1,3 , 5,7 /

      integer e
      character*1 c0,c1

      nskip=2
      if (if3d) nskip=1
c
c     Check "r" curves (edges 4 and 2)
c
         
      do i=1,4,nskip

         i0 = c_pair(1,i)
         i1 = c_pair(2,i)
         c0 = ccurve(i0,e)
         c1 = ccurve(i1,e)

         if (i.le.2) then
            r0 = ii-1
            mm = ii
            nn = nxsp
         else
            r0 = jj-1
            mm = jj
            nn = nysp
         endif
         dr = 1./nn
         r0 = r0/nn
         r1 = r0 + dr

c        write(6,1) r0,ccurve(i0,e),curve(1,i0,e),' r0 b4'
c        write(6,1) r1,ccurve(i1,e),curve(1,i1,e),' r1 b4'

         if (c0.eq.'C' .and. c1.eq.'C') then
            curve0 = abs(curve(1,i0,e))
            curve1 = abs(curve(1,i1,e))
            cnew_0 = (1.-r0)*curve0 + r0*curve1
            cnew_1 = (1.-r1)*curve0 + r1*curve1
            if (curve(1,i0,e).lt.0) cnew_0 = -cnew_0
            if (curve(1,i1,e).lt.0) cnew_1 = -cnew_1
            curve(1,i0,e) = cnew_0
            curve(1,i1,e) = cnew_1
         elseif (c0.eq.'C') then
            c_old         = curve(1,i0,e)
            curve(1,i0,e) = c_old/(1.-r0)
            if (mm.lt.nn) then
               ccurve(i1,e) = 'C'
               curve(1,i1,e) = -c_old/(1.-r1)
            endif
         elseif (c1.eq.'C') then
            c_old         = curve(1,i1,e)
            curve(1,i1,e) = c_old/r1
            if (mm.gt.1) then
               ccurve(i0,e) = 'C'
               curve(1,i0,e) = -c_old/r0
            endif
         endif

c        write(6,1) r0,ccurve(i0,e),curve(1,i0,e),' r0 af'
c        write(6,1) r1,ccurve(i1,e),curve(1,i1,e),' r1 af'
c  1     format(f12.4,2x,a1,2x,f12.4,a6)

      enddo

c     write(6,*) 'continue?'
c     read (5,*) c0

      return
      end
c-----------------------------------------------------------------------
