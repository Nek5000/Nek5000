c-----------------------------------------------------------------------
      subroutine plan4_acc (igeom)

C     Splitting scheme A.G. Tomboulides et al.
c     Journal of Sci.Comp.,Vol. 12, No. 2, 1998
c
C     NOTE: QTL denotes the so called thermal
c           divergence and has to be provided
c           by userqtl.
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      COMMON /SCRNS/ RES1  (LX1,LY1,LZ1,LELV)
     $ ,             RES2  (LX1,LY1,LZ1,LELV)
     $ ,             RES3  (LX1,LY1,LZ1,LELV)
     $ ,             DV1   (LX1,LY1,LZ1,LELV)
     $ ,             DV2   (LX1,LY1,LZ1,LELV)
     $ ,             DV3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)

      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP

      REAL DVC (LX1,LY1,LZ1,LELV), DFC(LX1,LY1,LZ1,LELV)
      REAL DIV1, DIV2, DIF1, DIF2, QTL1, QTL2

      INTEGER e

      ifxyo = .true.

      if (icalld.eq.0) tpres=0.0
      icalld=icalld+1
      npres=icalld

      intype = -1
      ntot1  = nx1*ny1*nz1*nelv
      n      = ntot1

      if (igeom.eq.1) then

         call plan4_acc_data_copyin
         call plan4_acc_update_device
         call hsmg_acc_update_device

         call makef_acc

         call sumab_acc(vx_e,vx,vxlag,n,ab,nab)
         call sumab_acc(vy_e,vy,vylag,n,ab,nab)
         call sumab_acc(vz_e,vz,vzlag,n,ab,nab)

!$acc    update host(vx_e,vy_e,vz_e)

      else
         ! add user defined divergence to qtl 
         call add2_acc (qtl,usrdiv,n)

         call lagvel_acc

!$acc    update host(vx,vy,vz,v1mask,v2mask,v3mask)
         call bcdirvc  (vx,vy,vz,v1mask,v2mask,v3mask) 
!$acc    update device(vx,vy,vz)


c        first, compute pressure
         if (icalld.eq.0) tpres=0.0
         icalld=icalld+1
         npres=icalld
         etime1=dnekclock()

         if (istep.eq.1) then
!$acc       enter data copyin(h1,h2,respr,pmask,res1,res2,res3)
!$acc       enter data copyin(dv1,dv2,dv3)
         endif

         call crespsp_acc(respr)

         call invers2_acc (h1,vtrans,n)
         call rzero_acc   (h2,n)
         call ctolspl_acc(tolspl,respr)

         call dssum     (respr,nx1,ny1,nz1)
         call col2_acc  (respr,pmask,n)

         iter=nmxh
         call hmh_gmres (respr,h1,h2,vmult,iter)

         call add2_acc (pr,respr,n)
         call ortho_acc(pr)

         tpres=tpres+(dnekclock()-etime1)

         call cresvsp_acc(res1,res2,res3,h1,h2)

         ! ophinv_pr starts here
         call dssum(res1,nx1,ny1,nz1)
         call dssum(res2,nx1,ny1,nz1)
         call dssum(res3,nx1,ny1,nz1)

         call col2_acc(res1,v1mask,ntot1)
         call col2_acc(res2,v2mask,ntot1)
         call col2_acc(res3,v3mask,ntot1)

         call chktcg1_acc(tol1,res1,h1,h2,v1mask,vmult,imesh,1)
         call chktcg1_acc(tol2,res1,h1,h2,v2mask,vmult,imesh,2)
         call chktcg1_acc(tol3,res1,h1,h2,v3mask,vmult,imesh,3)

         call cggo_acc(dv1,res1,h1,h2,v1mask,vmult,imesh,tol1,nmxh,1,
     $                 binvm1)
         call cggo_acc(dv2,res2,h1,h2,v2mask,vmult,imesh,tol2,nmxh,2,
     $                 binvm1)
         call cggo_acc(dv3,res3,h1,h2,v3mask,vmult,imesh,tol3,nmxh,3,
     $                 binvm1)
         ! ophinv_pr ends here

c below gives correct values in iterations
c but printed values are wierd  L1/L2 DIV(V) 6.9034-310   6.9034-310  

         call add2_acc  (vx,dv1,n)      
         call add2_acc  (vy,dv2,n)
         call add2_acc  (vz,dv3,n)

         call plan4_acc_update_host

         IF (NIO.EQ.0) THEN
            WRITE(6,'(13X,A,1p2e13.4)')
     &         'L1/L2 DIV(V)        ',DIV1,DIV2
            WRITE(6,'(13X,A,1p2e13.4)') 
     &         'L1/L2 QTL           ',QTL1,QTL2
            WRITE(6,'(13X,A,1p2e13.4)')
     &         'L1/L2 DIV(V)-QTL    ',DIF1,DIF2
            IF (DIF2.GT.0.1) WRITE(6,'(13X,A)') 
     &         'WARNING: DIV(V)-QTL too large!'
         ENDIF
 
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine plan4 (igeom)

C     Splitting scheme A.G. Tomboulides et al.
c     Journal of Sci.Comp.,Vol. 12, No. 2, 1998
c
C     NOTE: QTL denotes the so called thermal
c           divergence and has to be provided
c           by userqtl.
c
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'MVGEOM'
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
 
      REAL           DPR   (LX2,LY2,LZ2,LELV)
      EQUIVALENCE   (DPR,DV1)
      LOGICAL        IFSTSP

      REAL DVC (LX1,LY1,LZ1,LELV), DFC(LX1,LY1,LZ1,LELV)
      REAL DIV1, DIV2, DIF1, DIF2, QTL1, QTL2
c
      INTYPE = -1
      NTOT1  = NX1*NY1*NZ1*NELV


      call plan4_acc_data_copyin()


      if (igeom.eq.1) then

         ! compute explicit contributions bfx,bfy,bfz 
         call makef 

         call sumab(vx_e,vx,vxlag,ntot1,ab,nab)
         call sumab(vy_e,vy,vylag,ntot1,ab,nab)
         if (if3d) call sumab(vz_e,vz,vzlag,ntot1,ab,nab)

      else

         if(iflomach) call opcolv(bfx,bfy,bfz,vtrans)

         ! add user defined divergence to qtl 
         call add2 (qtl,usrdiv,ntot1)

         ! split viscosity into explicit/implicit part
         if (ifexplvis) call split_vis

         call lagvel

         ! mask Dirichlet boundaries
         call bcdirvc  (vx,vy,vz,v1mask,v2mask,v3mask) 

C        first, compute pressure

         if (icalld.eq.0) tpres=0.0
         icalld=icalld+1
         npres=icalld
         etime1=dnekclock()

         call crespsp  (respr)
         call invers2  (h1,vtrans,ntot1)
         call rzero    (h2,ntot1)
         call ctolspl  (tolspl,respr)
         napproxp(1) = laxtp
         call hsolve   ('PRES',dpr,respr,h1,h2 
     $                        ,pmask,vmult
     $                        ,imesh,tolspl,nmxh,1
     $                        ,approxp,napproxp,binvm1)
         call add2    (pr,dpr,ntot1)
         call ortho   (pr)

         tpres=tpres+(dnekclock()-etime1)

C        Compute velocity
         call cresvsp (res1,res2,res3,h1,h2)
c        call ophinv_pr(dv1,dv2,dv3,res1,res2,res3,h1,h2,tolhv,nmxh)
         call ophinv_pr_debug(dv1,dv2,dv3,res1,res2,res3,
     $                        h1,h2,tolhv,nmxh)
         call opadd2  (vx,vy,vz,dv1,dv2,dv3)

         if (ifexplvis) call redo_split_vis

c Below is just for diagnostics...

c        Calculate Divergence norms of new VX,VY,VZ
         CALL OPDIV   (DVC,VX,VY,VZ)
         CALL DSSUM   (DVC,NX1,NY1,NZ1)
         CALL COL2    (DVC,BINVM1,NTOT1)

         CALL COL3    (DV1,DVC,BM1,NTOT1)
         DIV1 = GLSUM (DV1,NTOT1)/VOLVM1

         CALL COL3    (DV2,DVC,DVC,NTOT1)
         CALL COL2    (DV2,BM1   ,NTOT1)
         DIV2 = GLSUM (DV2,NTOT1)/VOLVM1
         DIV2 = SQRT  (DIV2)
c        Calculate Divergence difference norms
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

         IF (NIO.EQ.0) THEN
            WRITE(6,'(13X,A,1p2e13.4)')
     &         'L1/L2 DIV(V)        ',DIV1,DIV2
            WRITE(6,'(13X,A,1p2e13.4)') 
     &         'L1/L2 QTL           ',QTL1,QTL2
            WRITE(6,'(13X,A,1p2e13.4)')
     &         'L1/L2 DIV(V)-QTL    ',DIF1,DIF2
            IF (DIF2.GT.0.1) WRITE(6,'(13X,A)') 
     &         'WARNING: DIV(V)-QTL too large!'
         ENDIF
 
      endif
 
      return
      END

c-----------------------------------------------------------------------
      subroutine crespsp (respr)

C     Compute startresidual/right-hand-side in the pressure

      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      REAL           RESPR (LX1*LY1*LZ1,LELV)
c
      COMMON /SCRNS/ TA1   (LX1*LY1*LZ1,LELV)
     $ ,             TA2   (LX1*LY1*LZ1,LELV)
     $ ,             TA3   (LX1*LY1*LZ1,LELV)
     $ ,             WA1   (LX1*LY1*LZ1*LELV)
     $ ,             WA2   (LX1*LY1*LZ1*LELV)
     $ ,             WA3   (LX1*LY1*LZ1*LELV)
      COMMON /SCRMG/ W1    (LX1*LY1*LZ1,LELV)
     $ ,             W2    (LX1*LY1*LZ1,LELV)
     $ ,             W3    (LX1*LY1*LZ1,LELV)

      common /scruz/         sij (lx1*ly1*lz1,6,lelv)
      parameter (lr=lx1*ly1*lz1)
      common /scrvz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)

      CHARACTER CB*3
      
      NXYZ1  = NX1*NY1*NZ1
      NTOT1  = NXYZ1*NELV
      NFACES = 2*NDIM

c      call lagvel

c     -mu*curl(curl(v))
      call op_curl (ta1,ta2,ta3,vx_e,vy_e,vz_e,.true.,w1,w2)
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

c compute stress tensor for ifstrs formulation - variable viscosity Pn-Pn
      if (ifstrs) then
         call opgrad   (ta1,ta2,ta3,vdiff)
         call invcol2  (ta1,vdiff,ntot1)
         call invcol2  (ta2,vdiff,ntot1)
         call invcol2  (ta3,vdiff,ntot1)

         nij = 3
         if (if3d.or.ifaxis) nij=6

         call comp_sij   (sij,nij,vx_e,vy_e,vz_e,
     &                    ur,us,ut,vr,vs,vt,wr,ws,wt)
         call col_mu_sij (w1,w2,w3,ta1,ta2,ta3,sij,nij)

         call opcolv   (ta1,ta2,ta3,QTL)
         scale2 = -2./3. 
         call opadd2cm (w1,w2,w3,ta1,ta2,ta3,scale2)
         call opsub2   (wa1,wa2,wa3,w1,w2,w3)

      endif

      call invcol3  (w1,vdiff,vtrans,ntot1)
      call opcolv   (wa1,wa2,wa3,w1)

c     add old pressure term because we solve for delta p 
      call invers2 (ta1,vtrans,ntot1)
      call rzero   (ta2,ntot1)

      call bcdirsc (pr)

      call axhelm  (respr,pr,ta1,ta2,imesh,1)
      call chsign  (respr,ntot1)

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
         call cdtp    (wa1,ta1,rxm1,sxm1,txm1,1)
         call cdtp    (wa2,ta2,rym1,sym1,tym1,1)
         call cdtp    (wa3,ta3,rzm1,szm1,tzm1,1)
         do i=1,n
            respr(i,1) = respr(i,1)+wa1(i)+wa2(i)+wa3(i)
         enddo
      else
         call cdtp    (wa1,ta1,rxm1,sxm1,txm1,1)
         call cdtp    (wa2,ta2,rym1,sym1,tym1,1)
         do i=1,n
            respr(i,1) = respr(i,1)+wa1(i)+wa2(i)
         enddo
      endif

C     add thermal divergence
      dtbd = BD(1)/DT
      call admcol3(respr,QTL,bm1,dtbd,ntot1)
 
C     surface terms
      DO 100 IEL=1,NELV
         DO 300 IFC=1,NFACES
            CALL RZERO  (W1(1,IEL),NXYZ1)
            CALL RZERO  (W2(1,IEL),NXYZ1)
            IF (NDIM.EQ.3)
     $      CALL RZERO  (W3(1,IEL),NXYZ1)
            CB = CBC(IFC,IEL,IFIELD)
            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v'.or.
     $         cb.eq.'MV '.or.cb.eq.'mv ') then
               CALL FACCL3
     $         (W1(1,IEL),VX(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3
     $         (W2(1,IEL),VY(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
               IF (NDIM.EQ.3)
     $          CALL FACCL3
     $         (W3(1,IEL),VZ(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
            ELSE IF (CB(1:3).EQ.'SYM') THEN
               CALL FACCL3
     $         (W1(1,IEL),TA1(1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3
     $         (W2(1,IEL),TA2(1,IEL),UNY(1,1,IFC,IEL),IFC)
               IF (NDIM.EQ.3)
     $          CALL FACCL3
     $         (W3(1,IEL),TA3(1,IEL),UNZ(1,1,IFC,IEL),IFC)
            ENDIF
            CALL ADD2   (W1(1,IEL),W2(1,IEL),NXYZ1)
            IF (NDIM.EQ.3)
     $      CALL ADD2   (W1(1,IEL),W3(1,IEL),NXYZ1)
            CALL FACCL2 (W1(1,IEL),AREA(1,1,IFC,IEL),IFC)
            IF (CB(1:1).EQ.'V'.OR.CB(1:1).EQ.'v'.or.
     $         cb.eq.'MV '.or.cb.eq.'mv ') then
              CALL CMULT(W1(1,IEL),dtbd,NXYZ1)
            endif
            CALL SUB2 (RESPR(1,IEL),W1(1,IEL),NXYZ1)
  300    CONTINUE
  100 CONTINUE

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

c     call outpost(resv1,resv2,resv3,h1,h2,'gv_')
      CALL SETHLM  (H1,H2,INTYPE)
c     call outpost(vx,vy,vz,h1,h2,'gv_')

      CALL OPHX    (RESV1,RESV2,RESV3,VX,VY,VZ,H1,H2)
c     call outpost(resv1,resv2,resv3,h1,h2,'gv_')
      CALL OPCHSGN (RESV1,RESV2,RESV3)

      scale = -1./3.
      if (ifstrs) scale =  2./3.

      call col3(ta4,vdiff,qtl,ntot)

      do i=1,ntot/nelv
c        write (6,*) 'pr=',pr(i,1,1,1)
      enddo

c     call outpost(ta4,vdiff,qtl,pr,t,'gv_')
      call add2s1 (ta4,pr,scale,ntot)    

      call outpost(resv1,resv2,resv3,ta4,h2,'gcv')
c     call opgrad (ta1,ta2,ta3,ta4)
c     do i=1,ntot/nelv
c        write (6,*) 'ta4=',ta4(i,1,1,1)
c     enddo
      call wgradm1(ta1,ta2,ta3,ta4,nelv)
c     call exitti('exit after wgradm1 in cresvsp$',1)
c     call outpost(ta1,ta2,ta3,ta4,h2,'gcv')

c     do i=1,ntot/nelv
c        write (6,*) 'ta1=',ta1(i,1,1,1)
c        write (6,*) 'ta2=',ta2(i,1,1,1)
c        write (6,*) 'ta3=',ta3(i,1,1,1)
c     enddo

      if(IFAXIS) then
         CALL COL2 (TA2, OMASK,NTOT)
         CALL COL2 (TA3, OMASK,NTOT)
      endif

      call opsub2  (resv1,resv2,resv3,ta1,ta2,ta3)
      call outpost(resv1,resv2,resv3,h1,h2,'gv_')
      call opadd2  (resv1,resv2,resv3,bfx,bfy,bfz)
      call outpost(resv1,resv2,resv3,h1,h2,'gv_')
      return
      end

c-----------------------------------------------------------------------
      subroutine op_curl(w1,w2,w3,u1,u2,u3,ifavg,work1,work2)

      include 'SIZE'
      include 'TOTAL'

      real duax(lx1), ta(lx1,ly1,lz1,lelv)

      logical ifavg

      real w1(1),w2(1),w3(1),work1(1),work2(1),u1(1),u2(1),u3(1)

      ntot  = nx1*ny1*nz1*nelv
      nxyz  = nx1*ny1*nz1

ccc!$acc update host   (u1(1:ntot),u2(1:ntot),u3(1:ntot))

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
                  call mxm   (ta(1,1,1,iel),nx1,DATM1,ny1,duax,1)
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

      if (ifavg .and. .not. ifcyclic) then ! average between elements

         ifielt = ifield
         ifield = 1

         call opcolv(w1,w2,w3,bm1)

!$acc update device(w1(1:ntot),w2(1:ntot),w3(1:ntot))
         call dssum   (w1,nx1,ny1,nz1)
         call dssum   (w2,nx1,ny1,nz1)
         call dssum   (w3,nx1,ny1,nz1)
!$acc update host(w1(1:ntot),w2(1:ntot),w3(1:ntot))
         call opcolv(w1,w2,w3,binvm1)

         ifield = ifielt

      endif


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
      subroutine col_mu_sij(w1,w2,w3,ta1,ta2,ta3,sij,nij)

      include 'SIZE'
      include 'TOTAL'

      parameter (lr=lx1*ly1*lz1)
      real w1 (lr,1),w2 (lr,1),w3 (lr,1)
      real ta1(lr,1),ta2(lr,1),ta3(lr,1)
      real sij(lr,nij,1)
      integer e

      nxyz1 = nx1*ny1*nz1

      if (if3d.or.ifaxis) then
        do e=1,nelv
          call vdot3 (w1(1,e),ta1(1,e),ta2(1,e),ta3(1,e)
     $               ,sij(1,1,e),sij(1,4,e),sij(1,6,e),nxyz1)
          call vdot3 (w2(1,e),ta1(1,e),ta2(1,e),ta3(1,e)
     $                  ,sij(1,4,e),sij(1,2,e),sij(1,5,e),nxyz1)
          call vdot3 (w3(1,e),ta1(1,e),ta2(1,e),ta3(1,e)
     $                  ,sij(1,6,e),sij(1,5,e),sij(1,3,e),nxyz1)
        enddo

      else

        do e=1,nelv
          call vdot2 (w1(1,e),ta1(1,e),ta2(1,e)
     $               ,sij(1,1,e),sij(1,3,e),nxyz1)
          call vdot2 (w2,ta1(1,e),ta2(1,e)
     $               ,sij(1,3,e),sij(1,2,e),nxyz1)
          call rzero (w3(1,e),nxyz1)
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sumab(v,vv,vvlag,ntot,ab_,nab_)
c
c     sum up AB/BDF contributions 
c
      include 'SIZE'

      real vvlag(lx1*ly1*lz1*lelv,*)
      real ab_(*)

      ab0 = ab_(1)
      ab1 = ab_(2)
      ab2 = ab_(3)

      call add3s2(v,vv,vvlag(1,1),ab0,ab1,ntot)
      if(nab_.eq.3) call add2s2(v,vvlag(1,2),ab2,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine userqtl_scig
c
c     Compute the thermal divergence QTL for an ideal single component gas 
c     QTL := div(v) = -1/rho * Drho/Dt
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      common /scrns/ w1(lx1,ly1,lz1,lelt)
     $              ,w2(lx1,ly1,lz1,lelt)
     $              ,w3(lx1,ly1,lz1,lelt)
     $              ,tx(lx1,ly1,lz1,lelt)
     $              ,ty(lx1,ly1,lz1,lelt)
     $              ,tz(lx1,ly1,lz1,lelt)

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelv

      ifld_save = ifield

c - - Assemble RHS of T-eqn
      ifield=2

      call makeuq
      call copy(qtl,bq,ntot)

      ifield=1     !set right gs handle (QTL is only defined on the velocity mesh)
      call opgrad  (tx,ty,tz,t)
      call opdssum (tx,ty,tz)
      call opcolv  (tx,ty,tz,binvm1)
      call opcolv  (tx,ty,tz,vdiff(1,1,1,1,2))
      call opdiv   (w2,tx,ty,tz)

      call add2    (qtl,w2,ntot)
      call dssum   (qtl,nx1,ny1,nz1)
      call col2    (qtl,binvm1,ntot)

      ! QTL = T_RHS/(rho*cp**T)
      call col3    (w2,vtrans(1,1,1,1,2),t,ntot)
      call invcol2 (qtl,w2,ntot)

      dp0thdt = 0.0
      if (ifdp0dt) then
         dd = (1.0 - gamma0)/gamma0
         call rone(w1,ntot)
         call cmult(w1,dd,ntot)

         call invcol3(w2,vtrans(1,1,1,1,2),vtrans,ntot)
         call invcol2(w1,w2,ntot)

         call cadd(w1,1.0,ntot)
         call copy(w2,w1,ntot)
         call col2(w1,bm1,ntot)
  
         p0alph1 = p0th / glsum(w1,ntot)
  
         call copy   (w1,qtl,ntot)
         call col2   (w1,bm1,ntot)

         termQ = glsum(w1,ntot)
         if (ifcvfun) then
            termV = glcflux(vx,vy,vz)
         else
            termV = glcflux(vx_e,vy_e,vz_e)
         endif
         dp0thdt = p0alph1*(termQ - termV)

         dd =-dp0thdt/p0th
         call cmult(w2,dd,ntot)
         call add2 (qtl,w2,ntot)
      endif

      ifield = ifld_save

      return
      end
c-----------------------------------------------------------------------
      subroutine qthermal
c
c     generic qtl wrapper
c
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'

      ntot = nx1*ny1*nz1*nelv

      call rzero(qtl,ntot)
      call userqtl()

      return
      end
c-----------------------------------------------------------------------
      subroutine sumab_acc(v,vv,vvlag,ntot,ab_,nab_)
c
c     sum up AB/BDF contributions 
c
      include 'SIZE'

      real v(lx1*ly1*lz1*lelv)
      real vv(lx1*ly1*lz1*lelv)
      real vvlag(lx1*ly1*lz1*lelv,2)
      real ab_(10)

      ab0 = ab_(1)
      ab1 = ab_(2)
      ab2 = ab_(3)

c     INLINED:
c     call add3s2_acc(v,vv,vvlag(1,1),ab0,ab1,ab2,ntot)
c     call add2s2_acc(v,vvlag(1,2),ab2,ntot)
!$acc parallel present(v,vv,vvlag)
      do i=1,ntot
         v(i) = ab0*vv(i) + ab1*vvlag(i,1) + ab2*vvlag(i,2)
      enddo
!$acc end parallel

      return
      end
c-----------------------------------------------------------------------
      subroutine crespsp_face_update(respr,ta1,ta2,ta3,dtbd,w1,w2,w3)
c
      include 'SIZE'
      include 'TOTAL'

      parameter(lz=lx1*ly1*lz1)
      real respr(lz,lelt),ta1(lz,lelt),ta2(lz,lelt),ta3(lz,lelt)
      real w1(lz),w2(lz),w3(lz)

      character*1 c1
      character*3 cb
      integer e,f

      nfaces = 2*ndim

      do e=1,nelv
      do f=1,nfaces

         call rzero  (w1,lz)
         call rzero  (w2,lz)
         if (ldim.eq.3) call rzero  (w3,lz)

         cb = cbc(f,e,ifield)
         c1 = cbc(f,e,ifield)

         if (c1.eq.'V'.or.c1.eq.'v'.or.cb.eq.'MV '.or.cb.eq.'mv ') then
            call faccl3(w1,vx(1,1,1,e),unx(1,1,f,e),f)
            call faccl3(w2,vy(1,1,1,e),uny(1,1,f,e),f)
            if (ldim.eq.3)
     $      call faccl3(w3,vz(1,1,1,e),unz(1,1,f,e),f)
         elseif (cb.eq.'SYM') then
            call faccl3(w1,ta1(1,e),unx(1,1,f,e),f)
            call faccl3(w2,ta2(1,e),uny(1,1,f,e),f)
            if (ldim.eq.3)
     $       call faccl3(w3,ta3(1,e),unz(1,1,f,e),f)
         endif

         call add2   (w1,w2,lz)
         if (ldim.eq.3) call add2   (w1,w3,lz)

         call faccl2 (w1,area(1,1,f,e),f)

         if (c1.eq.'V'.or.c1.eq.'v'.or.cb.eq.'MV '.or.cb.eq.'mv ') 
     $     call cmult(w1,dtbd,lz)
         
         call sub2 (respr(1,e),w1,lz)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine crespsp_acc(respr)

C     Compute startresidual/right-hand-side in the pressure

      include 'SIZE'
      include 'TOTAL'

      real           respr (lx1,ly1,lz1,lelv)
c
      real           ta1   (lx1,ly1,lz1,lelv)
     $ ,             ta2   (lx1,ly1,lz1,lelv)
     $ ,             ta3   (lx1,ly1,lz1,lelv)
     $ ,             wa1   (lx1*ly1*lz1*lelv)
     $ ,             wa2   (lx1*ly1*lz1*lelv)
     $ ,             wa3   (lx1*ly1*lz1*lelv)
     $ ,             wa4   (lx1*ly1*lz1*lelv*ldimt1)
      real           w1    (lx1,ly1,lz1,lelv)
     $ ,             w2    (lx1,ly1,lz1,lelv)
     $ ,             w3    (lx1,ly1,lz1,lelv)

      character cb*3,c1*1
      integer e,f

!$acc routine(facind) seq
      
      nxyz1  = lx1*ly1*lz1
      ntot1  = nxyz1*nelv
      nfaces = 2*ndim

!$acc enter data create(ta1,ta2,ta3,wa1,wa2,wa3,wa4,w1,w2,w3)
c     -mu*curl(curl(v))

c!$acc update host(ta1,ta2,ta3,vx_e,vy_e,vz_e)

c     write (6,*) 'synchronize'
c     do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'v_e',vx_e(i)
c        write (6,*) 'v_e',vy_e(i)
c        write (6,*) 'v_e',vz_e(i)
c     enddo
c     stop
c     call op_curl(ta1,ta2,ta3,vx_e,vy_e,vz_e,.true.,w1,w2)
c     call op_curl(wa1,wa2,wa3,ta1,ta2,ta3,.true.,w1,w2)
c!$acc update device(wa1,wa2,wa3)

      call op_curl_acc(ta1,ta2,ta3,vx_e,vy_e,vz_e,.true.,w1,w2)
      call op_curl_acc(wa1,wa2,wa3,ta1,ta2,ta3,.true.,w1,w2)

c     write (6,*) 'synchronize'

c!$acc update host(ta1,ta2,ta3,wa1,wa2,wa3)

c     do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'taa1 sync',i
c        write (6,*) 'taa1',wa1(i)
c        write (6,*) 'taa1',wa2(i)
c        write (6,*) 'taa1',wa3(i)
c     enddo

      call opcolv_acc(wa1,wa2,wa3,bm1)
      call wgradm1_acc(ta1,ta2,ta3,qtl,nelv)

      scale = -4./3. 
      call opadd2cm_acc(wa1,wa2,wa3,ta1,ta2,ta3,scale)

c!$acc update host(wa1,wa2,wa3)
c     do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'wa iii',i
c        write (6,*) 'wa1=',wa1(i)
c        write (6,*) 'wa2=',wa2(i)
c        write (6,*) 'wa3=',wa3(i)
c     enddo 

c compute stress tensor for ifstrs formulation - variable viscosity Pn-Pn
      if (ifstrs) 
     $   call exitti('ifstrs not yet support on gpu$',nelv)

      call invcol3_acc(w1,vdiff,vtrans,ntot1)

c!$acc update host(w1,wa1,wa2,wa3)
c
c      do i=1,lx1*ly1*lz1*nelv
c        write (6,*) 'w1=',w1(i,1,1,1)
c      enddo

      call opcolv_acc   (wa1,wa2,wa3,w1)

c     add old pressure term because we solve for delta p 

      call invers2_acc (ta1,vtrans,ntot1)
      call rzero_acc   (ta2,ntot1)

!$acc update host  (pr)
      call bcdirsc (pr)
!$acc update device(pr)

c!$acc update host(pr,ta1,ta2,respr)
c
c      call outpost(pr,ta1,ta2,respr,t,'wpa')

      call axhelm_acc_debug(respr,pr,ta1,ta2,1,1)

c      call outpost(respr,pr,ta1,ta2,t,'woa')

      call chsign_acc  (respr,ntot1)

c      call outpost(bfx,bfy,bfz,pr,t,'wb2')

c     add explicit (NONLINEAR) terms 

      n = nx1*ny1*nz1*nelv
!$acc parallel loop present(ta1,ta2,ta3,wa1,wa2,wa3)
      do i=1,n
         ta1(i,1,1,1) = bfx(i,1,1,1)/vtrans(i,1,1,1,1)-wa1(i)
         ta2(i,1,1,1) = bfy(i,1,1,1)/vtrans(i,1,1,1,1)-wa2(i)
         ta3(i,1,1,1) = bfz(i,1,1,1)/vtrans(i,1,1,1,1)-wa3(i)
      enddo
!$acc end parallel
c!$acc end paranvllel



c!$acc update host(ta1,ta2,ta3)
c      do i=1,n
c         write (6,*) 'taa1',ta1(i,1,1,1)
c         write (6,*) 'taa2',ta2(i,1,1,1)
c         write (6,*) 'taa3',ta3(i,1,1,1)
c      enddo
c      call outpost(ta1,ta2,ta3,pr,t,'wc0')

      call dssum (ta1,nx1,ny1,nz1)
      call dssum (ta2,nx1,ny1,nz1)
      call dssum (ta3,nx1,ny1,nz1)

c     call outpost(ta1,ta2,ta3,vtrans,t,'gc_')
c     call outpost(bm1,binvm1,wa1,wa2,wa3,'gc_')

!$acc parallel loop present(ta1,ta2,ta3,binvm1)
      do i=1,n
         ta1(i,1,1,1) = ta1(i,1,1,1)*binvm1(i,1,1,1)
         ta2(i,1,1,1) = ta2(i,1,1,1)*binvm1(i,1,1,1)
         ta3(i,1,1,1) = ta3(i,1,1,1)*binvm1(i,1,1,1)
      enddo
!$acc end parallel

      dtbd = bd(1)/dt  !! FOR NOW, no QTL support (pff, 7/31/17)
c     call admcol3(respr,qtl,bm1,dtbd,ntot1)

!$acc update host(ta1,ta2,ta3)
c******************************************
c TODO: still done on the cpu

c!$acc parallel loop gang
c!$acc& private(w1,w2,w3)
c!$acc parallel loop
c!$acc& present(ta1,ta2,ta3,w3m1,rxm1,rym1,rzm1)
c!$acc& present(sxm1,sym1,szm1,txm1,tym1,tzm1)
c!$acc& present(dxtm1,bm1,qtl,respr)

      do e=1,nelv
c!$acc  loop vector
       do i=1,lx1*ly1*lz1
         w1(i,1,1,e) = (rxm1(i,1,1,e)*ta1(i,1,1,e) ! Jacobian
     $               +rym1(i,1,1,e)*ta2(i,1,1,e) ! included
     $               +rzm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
         w2(i,1,1,e) = (sxm1(i,1,1,e)*ta1(i,1,1,e)
     $               +sym1(i,1,1,e)*ta2(i,1,1,e)
     $               +szm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
         w3(i,1,1,e) = (txm1(i,1,1,e)*ta1(i,1,1,e)
     $               +tym1(i,1,1,e)*ta2(i,1,1,e)
     $               +tzm1(i,1,1,e)*ta3(i,1,1,e))*w3m1(i,1,1)
       enddo

c!$acc  loop vector collapse(3)
       do k=1,nz1
       do j=1,ny1
       do i=1,nx1
          t1 = 0.0
c!$acc     loop seq
          do l=1,nx1
             t1 = t1 + dxm1(l,i)*w1(l,j,k,e) ! D^T
     $               + dxm1(l,j)*w2(i,l,k,e)
     $               + dxm1(l,k)*w3(i,j,l,e)
          enddo
          respr(i,j,k,e) = respr(i,j,k,e) + t1
     $                   + dtbd*bm1(i,j,k,e)*qtl(i,j,k,e)
       enddo
       enddo
       enddo

      enddo
c!$acc end parallel
c******************************************

c     call outpost(respr,wa1,wa2,wa3,t,'wc2')
c     stop

c!$acc update host(wa1,wa2,wa3,ta1,ta2,ta3,respr)
c      if (if3d) then
c         call cdtp    (wa1,ta1,rxm1,sxm1,txm1,1)
c         call cdtp    (wa2,ta2,rym1,sym1,tym1,1)
c         call cdtp    (wa3,ta3,rzm1,szm1,tzm1,1)
c         do i=1,n
c            respr(i,1,1,1) = respr(i,1,1,1)+wa1(i)+wa2(i)+wa3(i)
c         enddo
c      else
c         call cdtp    (wa1,ta1,rxm1,sxm1,txm1,1)
c         call cdtp    (wa2,ta2,rym1,sym1,tym1,1)
c         do i=1,n
c            respr(i,1,1,1) = respr(i,1,1,1)+wa1(i)+wa2(i)
c         enddo
c      endif
c!$acc update device(respr)

      dtbd = bd(1)/dt  !! FOR NOW, no QTL support (pff, 7/31/17)
c     call admcol3_acc(respr,qtl,bm1,dtbd,ntot1)

      call crespsp_face_update(respr,ta1,ta2,ta3,dtbd,w1,w2,w3)
c     call outpost(respr,wa1,wa2,wa3,t,'wc3')
c     stop
!$acc update device(respr)

C     Orthogonalize to (1,1,...,1)T for all-Dirichlet case

      call ortho_acc (respr)

!$acc update host(respr)

!$acc exit data

      return
      end
c----------------------------------------------------------------------
      subroutine invcol1_acc(a,n)
      real a(n)
!$acc parallel loop present(a)
      do i=1,n
         a(i)=1./a(i)
      enddo
!$acc end parallel
      return
      end
c-----------------------------------------------------------------------
      subroutine chktcg1_acc2(tol,res,h1,h2,mask,mult,imesh,isd)
C-------------------------------------------------------------------
C
C     Check that the tolerances are not too small for the CG-solver.
C     Important when calling the CG-solver (Gauss-Lobatto mesh) with
C     zero Neumann b.c.
C
C-------------------------------------------------------------------
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'EIGEN'

c     common  /cprint/ ifprint
c     logical          ifprint

      real w1   (lx1,ly1,lz1,lelt),
     $     w2   (lx1,ly1,lz1,lelt)

      real res  (lx1,ly1,lz1,lelt)
      real h1   (lx1,ly1,lz1,lelt)
      real h2   (lx1,ly1,lz1,lelt)
      real mult (lx1,ly1,lz1,lelt)
      real mask (lx1,ly1,lz1,lelt)

!$acc data create(w1,w2) present(res,h1,h2,mult,mask)

      if (eigaa.ne.0.) then
         acondno = eigga/eigaa
      else
         acondno = 10.
      endif

c     Single or double precision???

      delta = 1.e-9
      x     = 1.+delta
      y     = 1.
      diff  = abs(x-y)
      if (diff.eq.0.) eps = 1.e-6
      if (diff.gt.0.) eps = 1.e-13

      if (imesh.eq.1) then
          nl  = nelv
          vol = volvm1
      elseif (imesh.eq.2) then
          nl  = nelt
          vol = voltm1
      endif

      ntot1 = nx1*ny1*nz1*nl
      call copy (w1,res,ntot1)

      if (imesh.eq.1) then
         call col3_acc(w2,binvm1,w1,ntot1)
         rinit  = sqrt(glsc3_acc(w2,w1,mult,ntot1)/volvm1)
      else
         call col3_acc(w2,bintm1,w1,ntot1)
         rinit  = sqrt(glsc3_acc(w2,w1,mult,ntot1)/voltm1)
      endif

c     write (6,*) 'rinit=',rinit
c     stop

      rmin   = eps*rinit

c     if (tol.lt.rmin) then
c        if (nio.eq.0.and.ifprint)
c    $   write (6,*) 'New CG1-tolerance (RINIT*epsm) = ',rmin,tol
c        tol = rmin
c     endif

      call rone_acc(w1,ntot1)
      bcneu1 = glsc3_acc(w1,mask,mult,ntot1)
      bcneu2 = glsc3_acc(w1,w1  ,mult,ntot1)
      bctest = abs(bcneu1-bcneu2)

      call axhelm_acc(w2,w1,h1,h2,imesh,isd)
      call col2_acc(w2,w2,ntot1)
      call col2_acc(w2,bm1,ntot1)
      bcrob  = sqrt(glsum_acc(w2,ntot1)/vol)

      if ((bctest .lt. .1).and.(bcrob.lt.(eps*acondno))) then
c         OTR = GLSC3 (W1,RES,MULT,NTOT1)
         tolmin = rinit*eps*10.
         if (tol .lt. tolmin) then
             tol = tolmin
c            if (nio.eq.0.and.ifprint)
c    $       write(6,*) 'New CG1-tolerance (Neumann) = ',tolmin
         endif
      endif

!$acc end data

      return
      end
