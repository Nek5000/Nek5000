      subroutine plan4

C     Splitting scheme A.G. Tomboulides et al.
c     Journal of Sci.Comp.,Vol. 12, No. 2, 1998
c
C     NOTE: QTL denotes the so called thermal
c           divergence and has to be provided
c           by an external subroutine e.g qthermal
c
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOP'
C
      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELV)
     $ ,             RES2  (LX1,LY1,LZ1,LELV)
     $ ,             RES3  (LX1,LY1,LZ1,LELV)
     $ ,             DV1   (LX1,LY1,LZ1,LELV)
     $ ,             DV2   (LX1,LY1,LZ1,LELV)
     $ ,             DV3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCREV/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP

      REAL DVC (LX1,LY1,LZ1,LELV), DFC(LX1,LY1,LZ1,LELV)
      REAL DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
c
      INTYPE = -1
      NTOT1  = NX1*NY1*NZ1*NELV

      CALL MAKEF     ! nonlinear contributions, bfx, bfy, bfz
      CALL LAGVEL

      CALL BCDIRVC    (VX,VY,VZ,v1mask,v2mask,v3mask)
C
C     first, compute pressure

      call copy     (h1,vdiff,ntot1)
      call rzero    (h2,ntot1)

      call crespsp  (respr)
C
      call invers2  (h1,vtrans,ntot1)

      call ctolspl (tolspl,respr)

c     CALL HMHOLTZ ('PRES',DPR,RESPR,H1,H2,PMASK,VMULT,
c    $              IMESH,TOLSPL,NMXH,1)

      napprox(1) = laxt
      call hsolve  ('PRES',dpr,respr,h1,h2 
     $                    ,pmask,vmult
     $                    ,imesh,tolspl,nmxh,1
     $                    ,approx,napprox,binvm1)


      CALL ADD2    (PR,DPR,NTOT1)
      CALL ZAVER1  (PR)
C
C     Compute velocity
C
      CALL CRESVSP (RES1,RES2,RES3)

c      CALL SETHLM  (H1,H2,INTYPE)

      CALL OPHINV  (DV1,DV2,DV3,RES1,RES2,RES3,H1,H2,TOLHV,NMXH)
c
      CALL OPADD2  (VX,VY,VZ,DV1,DV2,DV3)

c
c     Below is just for diagnostics...
c
c calculate Divergence norms of new VX,VY,VZ
      CALL OPDIV   (DVC,VX,VY,VZ)
      CALL DSSUM   (DVC,NX1,NY1,NZ1)
      CALL COL2    (DVC,BINVM1,NTOT1)

      CALL COL3    (DV1,DVC,BM1,NTOT1)
      DIV1 = GLSUM (DV1,NTOT1)/VOLVM1

      CALL COL3    (DV2,DVC,DVC,NTOT1)
      CALL COL2    (DV2,BM1   ,NTOT1)
      DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
      DIV2 = SQRT  (DIV2)
c
c calculate Divergence difference norms

      CALL SUB3    (DFC,DVC,QTL,NTOT1)
      CALL COL3    (DV1,DFC,BM1,NTOT1)
      DIF1 = GLSUM (DV1,NTOT1)/VOLVM1
  
      CALL COL3    (DV2,DFC,DFC,NTOT1)
      CALL COL2    (DV2,BM1   ,NTOT1)
      DIF2 = GLSUM (DV2,NTOT1)/VOLVM1
      DIF2 = SQRT  (DIF2)

      CALL COL3    (DV1,QTL,BM1,NTOT1)
      QTL1 = GLSUM (DV1,NTOT1)/VOLVM1
  
      CALL COL3    (DV2,QTL,QTL,NTOT1)
      CALL COL2    (DV2,BM1   ,NTOT1)
      QTL2 = GLSUM (DV2,NTOT1)/VOLVM1
      QTL2 = SQRT  (QTL2)


      IF (NID .EQ. 0) WRITE(6,'(13X,A,1p2e13.4)')
     &                      'L1/L2 DIV(V)    :',DIV1,DIV2

      IF (NID .EQ. 0) WRITE(6,'(13X,A,1p2e13.4)') 
     &                      'L1/L2 QTL       :',QTL1,QTL2

      IF (NID .EQ. 0) WRITE(6,'(13X,A,1p2e13.4)')
     &                      'L1/L2 DIV(V)-QTL:',DIF1,DIF2
      

      RETURN
      END

c-----------------------------------------------------------------------
      SUBROUTINE CRESPSP (RESPR)
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the pressure
C     SPLIT SCHEME
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      REAL           RESPR (LX2,LY2,LZ2,LELV)
c
      COMMON /SCRNS/ TA1   (LX1,LY1,LZ1,LELV)
     $ ,             TA2   (LX1,LY1,LZ1,LELV)
     $ ,             TA3   (LX1,LY1,LZ1,LELV)
     $ ,             WA1   (LX1,LY1,LZ1,LELV)
     $ ,             WA2   (LX1,LY1,LZ1,LELV)
     $ ,             WA3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRMG/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
c
      CHARACTER CB*3
C
      NTOT1  = NX1*NY1*NZ1*NELV
      NTOT2  = NX2*NY2*NZ2*NELV
      NXYZ1  = NX1*NY1*NZ1
      NFACES = 2*NDIM

      call rzero(respr,ntot1)
      call rzero(ta1  ,ntot1)
      call rzero(ta2  ,ntot1)
      call rzero(ta3  ,ntot1)
      call rzero(wa1  ,ntot1)
      call rzero(wa2  ,ntot1)
      call rzero(wa3  ,ntot1)
c
c     Add -mu*curl(omega) to Dirichet boundaries
c
c     *** VISCOUS TERMS ***
c
      call op_curl  (ta1,ta2,ta3,vx ,vy ,vz ,.true. ,w1,w2)
      if(IFAXIS) then  
         CALL COL2 (TA2, OMASK,NTOT1)
         CALL COL2 (TA3, OMASK,NTOT1)
      endif

      call op_curl  (wa1,wa2,wa3,ta1,ta2,ta3,.true.,w1,w2)
      if(IFAXIS) then  
c         CALL COL2  (WA1, OMASK,NTOT1)
         CALL COL2  (WA2, OMASK,NTOT1)
         CALL COL2  (WA3, OMASK,NTOT1)
      endif

      call opcolv   (wa1,wa2,wa3,bm1)
c
      call opgrad   (ta1,ta2,ta3,QTL)
      if(IFAXIS) then  
         CALL COL2  (ta2, OMASK,ntot1)
         CALL COL2  (ta3, OMASK,ntot1)
      endif

      scale = -4./3. 
      call opadd2cm (wa1,wa2,wa3,ta1,ta2,ta3,scale)
      call invcol3  (w1,vdiff,vtrans,ntot1)
c
      call opcolv   (wa1,wa2,wa3,w1)
c
c     *** PRESSURE TERM ***
c
C                                                  solve for delta PP
      CALL INVERS2 (TA1,VTRANS,NTOT1)
      CALL RZERO   (TA2,NTOT1)
      CALL AXHELM  (RESPR,PR,TA1,TA2,IMESH,1)
      CALL CHSIGN  (RESPR,NTOT2)

c
c     *** NONLINEAR TERMS ***
C                                                  x-component
      call invcol3 (TA1,bfx,vtrans,ntot1)
      CALL SUB2    (TA1,WA1,NTOT1)
      CALL DSSUM   (TA1,NX1,NY1,NZ1)
      CALL COL2    (TA1,BINVM1,NTOT1)
      CALL CDTP    (TA2,TA1,RXM2,SXM2,TXM2,1)
      CALL ADD2    (RESPR,TA2,NTOT2)
C                                                  y-component
      call invcol3 (ta1,bfy,vtrans,ntot1)
      CALL SUB2    (TA1,WA2,NTOT1)
      CALL DSSUM   (TA1,NX1,NY1,NZ1)
      CALL COL2    (TA1,BINVM1,NTOT1)
      CALL CDTP    (TA2,TA1,RYM2,SYM2,TYM2,2)
      CALL ADD2    (RESPR,TA2,NTOT2)
C                                                  z-component
      IF (NDIM.EQ.3) THEN
         call invcol3 (ta1,bfz,vtrans,ntot1)
         CALL SUB2    (TA1,WA3,NTOT1)
         CALL DSSUM   (TA1,NX1,NY1,NZ1)
         CALL COL2    (TA1,BINVM1,NTOT1)
         CALL CDTP    (TA2,TA1,RZM2,SZM2,TZM2,3)
         CALL ADD2    (RESPR,TA2,NTOT2)
      ENDIF
C                                        add thermal divergence
      DTBD = BD(1)/DT
      call admcol3(respr,QTL,bm1,dtbd,ntot1)
 
C                                                 surface terms
      DO 300 IFC=1,NFACES
         CALL RZERO  (TA1,NTOT1)
         CALL RZERO  (TA2,NTOT1)
         IF (NDIM.EQ.3)
     $   CALL RZERO  (TA3,NTOT1)
         DO 100 IEL=1,NELV
            CB = CBC(IFC,IEL,IFIELD)
            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v') THEN
               CALL FACCL3 
     $         (TA1(1,1,1,IEL),VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3 
     $         (TA2(1,1,1,IEL),VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
               IF (NDIM.EQ.3) 
     $          CALL FACCL3 
     $         (TA3(1,1,1,IEL),VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
            ENDIF
            CALL ADD2   (TA1(1,1,1,IEL),TA2(1,1,1,IEL),NXYZ1)
            IF (NDIM.EQ.3)
     $      CALL ADD2   (TA1(1,1,1,IEL),TA3(1,1,1,IEL),NXYZ1)
            CALL FACCL2 (TA1(1,1,1,IEL),AREA(1,1,IFC,IEL),IFC)
  100    CONTINUE
         CALL CMULT(TA1,dtbd,NTOT1)
         CALL SUB2 (RESPR,TA1,NTOT1)
  300 CONTINUE

C     Assure that the residual is orthogonal to (1,1,...,1)T 
C     (only if all Dirichlet b.c.)
C
      CALL ORTHO (RESPR)

C
      RETURN
      END

c----------------------------------------------------------------------
      SUBROUTINE CRESVSP (RESV1,RESV2,RESV3)
C----------------------------------------------------------------------
C
C     Compute the residual for the velocity - SPLIT SCHEME.
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      REAL           RESV1 (LX1,LY1,LZ1,LELV)
      REAL           RESV2 (LX1,LY1,LZ1,LELV)
      REAL           RESV3 (LX1,LY1,LZ1,LELV)
      COMMON /SCRUZ/ TA1   (LX1,LY1,LZ1,LELV)
     $ ,             TA2   (LX1,LY1,LZ1,LELV)
     $ ,             TA3   (LX1,LY1,LZ1,LELV)
     $ ,             TA4   (LX1,LY1,LZ1,LELV)
      COMMON /SCREV/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)


      NTOT = NX1*NY1*NZ1*NELV
      INTYPE = -1

      CALL SETHLM  (H1,H2,INTYPE)

      CALL OPHX    (RESV1,RESV2,RESV3,VX,VY,VZ,H1,H2)
      CALL OPCHSGN (RESV1,RESV2,RESV3)

      scale = -1./3.
      call col3    (ta4,vdiff,qtl,ntot)
      call add2s1  (ta4,pr,scale,ntot)    
      call opgrad  (ta1,ta2,ta3,TA4)
      if(IFAXIS) then
         CALL COL2 (TA2, OMASK,NTOT)
         CALL COL2 (TA3, OMASK,NTOT)
      endif
c
      call opsub2  (resv1,resv2,resv3,ta1,ta2,ta3)
      call opadd2  (resv1,resv2,resv3,bfx,bfy,bfz)
C
      RETURN
      END

c-----------------------------------------------------------------------
      subroutine op_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)
c
      include 'SIZE'
      include 'TOTAL'
c
      real duax(lx1), ta(lx1,ly1,lz1,lelv)

      logical ifavg
c
      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)
c
      ntot  = nx1*ny1*nz1*nelv
      nxyz  = nx1*ny1*nz1
c     work1=dw/dy ; work2=dv/dz
        call dudxyz(work1,u3,rym1,sym1,tym1,jacm1,1,2)
        if (if3d) then
           call dudxyz(work2,u2,rzm1,szm1,tzm1,jacm1,1,3)
           call sub3(w1,work1,work2,ntot)
        else
           call copy(w1,work1,ntot)

           if(ifaxis) then
              call copy (ta,u3,ntot)
              do iel = 1,nelv
                if(IFRZER(iel)) then
                  call rzero (ta(1,1,1,iel),nx1)
                  call MXM   (ta(1,1,1,iel),nx1,DATM1,ny1,duax,1)
                  call copy  (ta(1,1,1,iel),duax,nx1)
                endif
                call col2    (ta(1,1,1,iel),yinvm1(1,1,1,iel),nxyz)
              enddo
              call add2      (w1,ta,ntot)
           endif

        endif
c     work1=du/dz ; work2=dw/dx
        if (if3d) then
           call dudxyz(work1,u1,rzm1,szm1,tzm1,jacm1,1,3)
           call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
           call sub3(w2,work1,work2,ntot)
        else
           call rzero (work1,ntot)
           call dudxyz(work2,u3,rxm1,sxm1,txm1,jacm1,1,1)
           call sub3(w2,work1,work2,ntot)
        endif
c     work1=dv/dx ; work2=du/dy
        call dudxyz(work1,u2,rxm1,sxm1,txm1,jacm1,1,1)
        call dudxyz(work2,u1,rym1,sym1,tym1,jacm1,1,2)
        call sub3(w3,work1,work2,ntot)
c
c    Avg at bndry
c
      if (ifavg) then
         ifielt = ifield
         ifield = 1
c
         call col2  (w1,bm1,ntot)
         call dssum (w1,nx1,ny1,nz1)
         call col2  (w1,binvm1,ntot)
c
         call col2  (w2,bm1,ntot)
         call dssum (w2,nx1,ny1,nz1)
         call col2  (w2,binvm1,ntot)
c
         call col2  (w3,bm1,ntot)
         call dssum (w3,nx1,ny1,nz1)
         call col2  (w3,binvm1,ntot)
c
         ifield = ifielt
      endif
c
      return
      end

      SUBROUTINE OPADD2CM (A1,A2,A3,B1,B2,B3,C)
      INCLUDE 'SIZE'
      REAL A1(1),A2(1),A3(1),B1(1),B2(1),B3(1),C
      NTOT1=NX1*NY1*NZ1*NELV
      if (ndim.eq.3) then
         do i=1,ntot1
            a1(i) = a1(i) + b1(i)*c
            a2(i) = a2(i) + b2(i)*c
            a3(i) = a3(i) + b3(i)*c
         enddo
      else
         do i=1,ntot1
            a1(i) = a1(i) + b1(i)*c
            a2(i) = a2(i) + b2(i)*c
         enddo
      endif
      RETURN
      END

      SUBROUTINE ADMCOL3(A,B,C,D,N)
      REAL A(1),B(1),C(1),D
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'admcol3'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
C
      DO 100 I=1,N
         A(I)=A(I)+B(I)*C(I)*D
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
