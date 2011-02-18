c-----------------------------------------------------------------------
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
      INCLUDE 'CTIMER'
C
      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELV)
     $ ,             RES2  (LX1,LY1,LZ1,LELV)
     $ ,             RES3  (LX1,LY1,LZ1,LELV)
     $ ,             DV1   (LX1,LY1,LZ1,LELV)
     $ ,             DV2   (LX1,LY1,LZ1,LELV)
     $ ,             DV3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      COMMON /SCRSF/ VEXT  (LX1*LY1*LZ1*LELV,3)
      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP

      REAL DVC (LX1,LY1,LZ1,LELV), DFC(LX1,LY1,LZ1,LELV)
      REAL DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
c
      INTYPE = -1
      NTOT1  = NX1*NY1*NZ1*NELV

      ! add user defined divergence to qtl 
      call add2 (qtl,usrdiv,ntot1)

      CALL V_EXTRAP(vext)

      ! compute explicit contributions bfx,bfy,bfz 
      CALL MAKEF 
      CALL LAGVEL

      ! split viscosity into explicit/implicit part
      if (ifexplvis) call split_vis

      ! extrapolate velocity

      ! mask Dirichlet boundaries
      CALL BCDIRVC  (VX,VY,VZ,v1mask,v2mask,v3mask) 

C     first, compute pressure
#ifndef NOTIMER
      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld
      etime1=dnekclock()
#endif

      call crespsp  (respr)
      call invers2  (h1,vtrans,ntot1)
      call rzero    (h2,ntot1)
      call ctolspl  (tolspl,respr)
      napprox(1) = laxt
      call hsolve   ('PRES',dpr,respr,h1,h2 
     $                     ,pmask,vmult
     $                     ,imesh,tolspl,nmxh,1
     $                     ,approx,napprox,binvm1)
      call add2    (pr,dpr,ntot1)
      call ortho   (pr)
#ifndef NOTIMER
      tpres=tpres+(dnekclock()-etime1)
#endif

C     Compute velocity
      call cresvsp (res1,res2,res3,h1,h2)
      call ophinv  (dv1,dv2,dv3,res1,res2,res3,h1,h2,tolhv,nmxh)
      call opadd2  (vx,vy,vz,dv1,dv2,dv3)

      if (ifexplvis) call redo_split_vis

c Below is just for diagnostics...

c     Calculate Divergence norms of new VX,VY,VZ
      CALL OPDIV   (DVC,VX,VY,VZ)
      CALL DSSUM   (DVC,NX1,NY1,NZ1)
      CALL COL2    (DVC,BINVM1,NTOT1)

      CALL COL3    (DV1,DVC,BM1,NTOT1)
      DIV1 = GLSUM (DV1,NTOT1)/VOLVM1

      CALL COL3    (DV2,DVC,DVC,NTOT1)
      CALL COL2    (DV2,BM1   ,NTOT1)
      DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
      DIV2 = SQRT  (DIV2)
c     Calculate Divergence difference norms
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

      IF (NID.EQ.0) THEN
         WRITE(6,'(15X,A,1p2e13.4)')
     &      'L1/L2 DIV(V)    :',DIV1,DIV2
         WRITE(6,'(15X,A,1p2e13.4)') 
     &      'L1/L2 QTL       :',QTL1,QTL2
         WRITE(6,'(15X,A,1p2e13.4)')
     &      'L1/L2 DIV(V)-QTL:',DIF1,DIF2
         IF (DIF2.GT.0.1) WRITE(6,'(15X,A)') 
     &          'WARNING: DIV(V)-QTL too large!'
      ENDIF
 
 
      return
      END

c-----------------------------------------------------------------------
      subroutine crespsp (respr)

C     Compute startresidual/right-hand-side in the pressure

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      REAL           RESPR (LX2*LY2*LZ2*LELV)
c
      COMMON /SCRNS/ TA1   (LX1*LY1*LZ1,LELV)
     $ ,             TA2   (LX1*LY1*LZ1,LELV)
     $ ,             TA3   (LX1*LY1*LZ1,LELV)
     $ ,             WA1   (LX1*LY1*LZ1*LELV)
     $ ,             WA2   (LX1*LY1*LZ1*LELV)
     $ ,             WA3   (LX1*LY1*LZ1*LELV)
      COMMON /SCRMG/ W1    (LX1*LY1*LZ1*LELV)
     $ ,             W2    (LX1*LY1*LZ1*LELV)
      COMMON /SCRSF/ VEXT  (LX1*LY1*LZ1*LELV,3)

      CHARACTER CB*3
      
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NXYZ1*NELV
      NFACES = 2*NDIM

c     -mu*curl(curl(v))
      call op_curl (ta1,ta2,ta3,vext(1,1),vext(1,2),vext(1,3),
     &              .true.,w1,w2)
      if(IFAXIS) then  
         CALL COL2 (TA2, OMASK,NTOT1)
         CALL COL2 (TA3, OMASK,NTOT1)
      endif
      call op_curl  (wa1,wa2,wa3,ta1,ta2,ta3,.true.,w1,w2)
      if(IFAXIS) then  
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
      call opcolv   (wa1,wa2,wa3,w1)

c     add old pressure term because we solve for delta p 
      CALL INVERS2 (TA1,VTRANS,NTOT1)
      CALL RZERO   (TA2,NTOT1)
      CALL AXHELM  (RESPR,PR,TA1,TA2,IMESH,1)
      CALL CHSIGN  (RESPR,NTOT1)

c     add explicit (NONLINEAR) terms 
      n = nx1*ny1*nz1*nelv
      do i=1,n
         ta1(i,1) = bfx(i,1,1,1)/vtrans(i,1,1,1,1)-wa1(i)
         ta2(i,1) = bfy(i,1,1,1)/vtrans(i,1,1,1,1)-wa2(i)
         ta3(i,1) = bfz(i,1,1,1)/vtrans(i,1,1,1,1)-wa3(i)
      enddo
      call opdssum (ta1,ta2,ta3)
      do i=1,n
         ta1(i,1) = ta1(i,1)*binvm1(i,1,1,1)
         ta2(i,1) = ta2(i,1)*binvm1(i,1,1,1)
         ta3(i,1) = ta3(i,1)*binvm1(i,1,1,1)
      enddo
      if (if3d) then
         call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
         call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
         call cdtp    (wa3,ta3,rzm2,szm2,tzm2,1)
         do i=1,n
            respr(i) = respr(i)+wa1(i)+wa2(i)+wa3(i)
         enddo
      else
         call cdtp    (wa1,ta1,rxm2,sxm2,txm2,1)
         call cdtp    (wa2,ta2,rym2,sym2,tym2,1)
         do i=1,n
            respr(i) = respr(i)+wa1(i)+wa2(i)
         enddo
      endif

C     add thermal divergence
      dtbd = BD(1)/DT
      call admcol3(respr,QTL,bm1,dtbd,ntot1)
 
C     surface terms
      DO 300 IFC=1,NFACES
         CALL RZERO  (TA1,NTOT1)
         CALL RZERO  (TA2,NTOT1)
         IF (NDIM.EQ.3)
     $   CALL RZERO  (TA3,NTOT1)
         DO 100 IEL=1,NELV
            CB = CBC(IFC,IEL,IFIELD)
            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v') THEN
               CALL FACCL3 
     $         (TA1(1,IEL),VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3 
     $         (TA2(1,IEL),VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
               IF (NDIM.EQ.3) 
     $          CALL FACCL3 
     $         (TA3(1,IEL),VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
            ENDIF
            CALL ADD2   (TA1(1,IEL),TA2(1,IEL),NXYZ1)
            IF (NDIM.EQ.3)
     $      CALL ADD2   (TA1(1,IEL),TA3(1,IEL),NXYZ1)
            CALL FACCL2 (TA1(1,IEL),AREA(1,1,IFC,IEL),IFC)
  100    CONTINUE
         CALL CMULT(TA1,dtbd,NTOT1)
         CALL SUB2 (RESPR,TA1,NTOT1)
  300 CONTINUE

C     Assure that the residual is orthogonal to (1,1,...,1)T 
C     (only if all Dirichlet b.c.)
      CALL ORTHO (RESPR)

      return
      END
c----------------------------------------------------------------------
      subroutine cresvsp (resv1,resv2,resv3,h1,h2)

C     Compute the residual for the velocity

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      real resv1(lx1,ly1,lz1,lelv)
     $   , resv2(lx1,ly1,lz1,lelv)
     $   , resv3(lx1,ly1,lz1,lelv)
     $   , h1   (lx1,ly1,lz1,lelv)
     $   , h2   (lx1,ly1,lz1,lelv)

      COMMON /SCRUZ/ TA1   (LX1,LY1,LZ1,LELV)
     $ ,             TA2   (LX1,LY1,LZ1,LELV)
     $ ,             TA3   (LX1,LY1,LZ1,LELV)
     $ ,             TA4   (LX1,LY1,LZ1,LELV)

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
      return
      end

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
c     if (ifavg) then
      if (ifavg .and. .not. ifcyclic) then

         ifielt = ifield
         ifield = 1
       
         call opcolv  (w1,w2,w3,bm1)
         call opdssum (w1,w2,w3)
         call opcolv  (w1,w2,w3,binvm1)

         ifield = ifielt

      endif
c
      return
      end

c-----------------------------------------------------------------------
      subroutine opadd2cm (a1,a2,a3,b1,b2,b3,c)
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
      return
      END

c-----------------------------------------------------------------------
      subroutine v_extrap(vext)
c
c     extrapolate velocity
c
      include 'SIZE'
      include 'TOTAL'
   
      real vext(lx1*ly1*lz1*lelv,1) 

      NTOT = NX1*NY1*NZ1*NELV

      AB0 = AB(1)
      AB1 = AB(2)
      AB2 = AB(3)

c     call copy(vext(1,1),vx,ntot)
c     call copy(vext(1,2),vy,ntot)
c     call copy(vext(1,3),vz,ntot)
c     return

      call add3s2(vext(1,1),vx,vxlag,ab0,ab1,ntot)
      call add3s2(vext(1,2),vy,vylag,ab0,ab1,ntot)
      if(if3d) call add3s2(vext(1,3),vz,vzlag,ab0,ab1,ntot)

      if(nab.eq.3) then
        call add2s2(vext(1,1),vxlag(1,1,1,1,2),ab2,ntot)
        call add2s2(vext(1,2),vylag(1,1,1,1,2),ab2,ntot)
        if(if3d) call add2s2(vext(1,3),vzlag(1,1,1,1,2),ab2,ntot)
      endif

      return
      end  
c-----------------------------------------------------------------------
      subroutine split_vis

c     Split viscosity into a constant implicit (VDIFF) and variable 
c     explicit (VDIFF_E) part.

      include 'SIZE'
      include 'TOTAL'

      n = nx1*ny1*nz1*nelv

      dnu_star = -nu_star
      call cadd2 (vdiff_e,vdiff,dnu_star,n)   ! set explicit part

      call cfill (vdiff,nu_star,n)            ! set implicit part

      return
      end
c-----------------------------------------------------------------------
      subroutine redo_split_vis     !     Redo split viscosity

      include 'SIZE'
      include 'TOTAL'

      n = nx1*ny1*nz1*nelv
      call add2(vdiff,vdiff_e,n) ! sum up explicit and implicit part

      return
      end
c-----------------------------------------------------------------------
