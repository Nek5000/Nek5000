      subroutine jacobian
      include 'basics.inc'
      common  
C
      nx3 = 2
      ny3 = 2
      nz3 = 2
      if (ndim.eq.2) nz3=1
c
      call zwgll(z3,w3,nx3)
      call dgll (dxm3,dxmt3,
c



      NXY3  = NX3*NY3
      NYZ3  = NY3*NZ3
      NXYZ3 = NX3*NY3*NZ3
      NTOT3 = NXYZ3*NELT
      NXYZ1 = NX1*NY1*NZ1
      NTOT1 = NXYZ1*NELT
C
C
C     Compute isoparametric partials.
C
      IF (NDIM.EQ.2) THEN
C
C     Two-dimensional case
C
      DO 200 ie=1,NELT
C
C     Use the appropriate derivative- and interpolation operator in
C     the y-direction (= radial direction if axisymmetric).
C
      IF (IFAXIS) THEN
         NY33   = NY3*NY3
         IF (IFRZER(ie)) THEN
            CALL COPY (DYTM3,DATM3,NY33)
         ELSE
            CALL COPY (DYTM3,DCTM3,NY33)
         ENDIF
      ENDIF
C
      CALL MXM(DXM3,NX3,XM3(1,1,1,ie),NX3,XRM3(1,1,1,ie),NY3)
      CALL MXM(DXM3,NX3,YM3(1,1,1,ie),NX3,YRM3(1,1,1,ie),NY3)
      CALL MXM(XM3(1,1,1,ie),NX3,DYTM3,NY3,XSM3(1,1,1,ie),NY3)
      CALL MXM(YM3(1,1,1,ie),NX3,DYTM3,NY3,YSM3(1,1,1,ie),NY3)
C
 200  CONTINUE
C
      CALL RZERO   (JACM3,NTOT3)
      CALL ADDCOL3 (JACM3,XRM3,YSM3,NTOT3)
      CALL SUBCOL3 (JACM3,XSM3,YRM3,NTOT3)
C
      ELSE
C
C     Three-dimensional case
C
      DO 300 ie=1,NELT
C
      CALL MXM(DXM3,NX3,XM3(1,1,1,ie),NX3,XRM3(1,1,1,ie),NYZ3)
      CALL MXM(DXM3,NX3,YM3(1,1,1,ie),NX3,YRM3(1,1,1,ie),NYZ3)
      CALL MXM(DXM3,NX3,ZM3(1,1,1,ie),NX3,ZRM3(1,1,1,ie),NYZ3)
C
      DO 310 IZ=1,NZ3
      CALL MXM(XM3(1,1,IZ,ie),NX3,DYTM3,NY3,XSM3(1,1,IZ,ie),NY3)
      CALL MXM(YM3(1,1,IZ,ie),NX3,DYTM3,NY3,YSM3(1,1,IZ,ie),NY3)
      CALL MXM(ZM3(1,1,IZ,ie),NX3,DYTM3,NY3,ZSM3(1,1,IZ,ie),NY3)
 310  CONTINUE
C
      CALL MXM(XM3(1,1,1,ie),NXY3,DZTM3,NZ3,XTM3(1,1,1,ie),NZ3)
      CALL MXM(YM3(1,1,1,ie),NXY3,DZTM3,NZ3,YTM3(1,1,1,ie),NZ3)
      CALL MXM(ZM3(1,1,1,ie),NXY3,DZTM3,NZ3,ZTM3(1,1,1,ie),NZ3)
C
 300  CONTINUE
C
      CALL RZERO   (JACM3,NTOT3)
      CALL ADDCOL4 (JACM3,XRM3,YSM3,ZTM3,NTOT3)
      CALL ADDCOL4 (JACM3,XTM3,YRM3,ZSM3,NTOT3)
      CALL ADDCOL4 (JACM3,XSM3,YTM3,ZRM3,NTOT3)
      CALL SUBCOL4 (JACM3,XRM3,YTM3,ZSM3,NTOT3)
      CALL SUBCOL4 (JACM3,XSM3,YRM3,ZTM3,NTOT3)
      CALL SUBCOL4 (JACM3,XTM3,YSM3,ZRM3,NTOT3)
C
      ENDIF
C
      DO 400 ie=1,NELT
         CALL CHKJAC(JACM3(1,1,1,ie),NXYZ3,ie,xm3,ym3,zm3,ndim)
 400  CONTINUE
C
      RETURN
      END
