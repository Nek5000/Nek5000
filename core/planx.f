      SUBROUTINE PLAN3 (IGEOM)
C-----------------------------------------------------------------------
C
C     Compute pressure and velocity using consistent approximation spaces.     
C     Operator splitting technique.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'EIGEN'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      COMMON /SCRNS/  RESV1 (LX1,LY1,LZ1,LELV)
     $ ,              RESV2 (LX1,LY1,LZ1,LELV)
     $ ,              RESV3 (LX1,LY1,LZ1,LELV)
     $ ,              DV1   (LX1,LY1,LZ1,LELV)
     $ ,              DV2   (LX1,LY1,LZ1,LELV)
     $ ,              DV3   (LX1,LY1,LZ1,LELV)
      COMMON /SCRVH/  H1    (LX1,LY1,LZ1,LELV)
     $ ,              H2    (LX1,LY1,LZ1,LELV)
C
      IF (IGEOM.EQ.1) THEN
C
C        Old geometry
C
         CALL MAKEF
C
      ELSE
C
C        New geometry, new b.c.
C
         INTYPE = -1
         CALL SETHLM  (H1,H2,INTYPE)
         CALL CRESVIF (RESV1,RESV2,RESV3,H1,H2)
         mstep = abs(param(94))
         if (param(94).ne.0. .and. istep.ge.mstep) then
           CALL OPHINVpr(DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
c          CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
         else
           CALL OPHINV  (DV1,DV2,DV3,RESV1,RESV2,RESV3,H1,H2,TOLHV,NMXH)
         endif
         CALL OPADD2  (VX,VY,VZ,DV1,DV2,DV3)
c
c        Default Filtering
c
c        alpha_filt = 0.05
c        if (param(103).ne.0.) alpha_filt=param(103)
c        call q_filter(alpha_filt)
c
c        CALL SSNORMD (DV1,DV2,DV3)
c
         CALL INCOMPR
C
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE LAGPRES 
C--------------------------------------------------------------------
C
C     Keep old pressure values
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
C
      IF (NBDINP.EQ.3) THEN
         NTOT2 = NX2*NY2*NZ2*NELV
         CALL COPY (PRLAG,PR,NTOT2)
      ENDIF
      RETURN
      END
C
      SUBROUTINE CRESVIF (RESV1,RESV2,RESV3,H1,H2)
C---------------------------------------------------------------------
C
C     Compute startresidual/right-hand-side in the velocity solver
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      REAL           RESV1 (LX1,LY1,LZ1,1)
      REAL           RESV2 (LX1,LY1,LZ1,1)
      REAL           RESV3 (LX1,LY1,LZ1,1)
      REAL           H1    (LX1,LY1,LZ1,1)
      REAL           H2    (LX1,LY1,LZ1,1)
      COMMON /SCRUZ/ W1    (LX1,LY1,LZ1,LELV)
     $ ,             W2    (LX1,LY1,LZ1,LELV)
     $ ,             W3    (LX1,LY1,LZ1,LELV)
     $ ,             PREXTR(LX2,LY2,LZ2,LELV)
C
      NTOT1 = NX1*NY1*NZ1*NELV
      NTOT2 = NX2*NY2*NZ2*NELV
      CALL LAGVEL 
      CALL BCDIRVC (VX,VY,VZ,v1mask,v2mask,v3mask)
      IF (IFSTRS)  CALL BCNEUTR
C
      CALL EXTRAPP (PREXTR)
      CALL OPGRADT (RESV1,RESV2,RESV3,PREXTR)
      CALL OPADD2  (RESV1,RESV2,RESV3,BFX,BFY,BFZ)
      CALL OPHX    (W1,W2,W3,VX,VY,VZ,H1,H2)
      CALL OPSUB2  (RESV1,RESV2,RESV3,W1,W2,W3)
C
      RETURN
      END
C
      SUBROUTINE EXTRAPP (PREXTR)
C--------------------------------------------------------------------
C
C     Pressure extrapolation
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      COMMON /CTMP0/ DPR (LX2,LY2,LZ2,LELV)
      REAL        PREXTR (LX2,LY2,LZ2,LELV)
C
      NTOT2 = NX2*NY2*NZ2*NELV
      IF (NBD.EQ.1.OR.NBD.EQ.2) THEN
         CALL COPY (PREXTR,PR,NTOT2)
      ELSEIF (NBD.EQ.3) THEN
         CONST = DTLAG(1)/DTLAG(2)
         CALL SUB3 (DPR,PR,PRLAG,NTOT2)
         CALL CMULT(DPR,CONST,NTOT2)
         CALL ADD3 (PREXTR,PR,DPR,NTOT2)
      ELSEIF (NBD.GT.3) THEN
         WRITE (6,*) 'Pressure extrapolation cannot be completed'
         WRITE (6,*) 'Try a lower-order temporal scheme'
         call exitt
      ENDIF
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine ophinvpr(ot1,ot2,ot3,in1,in2,in3,h1,h2,tolh,nmxi)
C
C     OT = (H1*A+H2*B)-1 * IN  (implicit)
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ORTHOV'
      INCLUDE 'TSTEP'
      INCLUDE 'SOLN'
c
      REAL OT1 (LX1,LY1,LZ1,1)
      REAL OT2 (LX1,LY1,LZ1,1)
      REAL OT3 (LX1,LY1,LZ1,1)
      REAL IN1 (LX1,LY1,LZ1,1)
      REAL IN2 (LX1,LY1,LZ1,1)
      REAL IN3 (LX1,LY1,LZ1,1)
      REAL H1  (LX1,LY1,LZ1,1)
      REAL H2  (LX1,LY1,LZ1,1)
c
c
      IMESH = 1
C
      IF (IFSTRS) THEN
         MATMOD = 0
         CALL HMHZSF  ('NOMG',OT1,OT2,OT3,IN1,IN2,IN3,H1,H2,
     $                  V1MASK,V2MASK,V3MASK,VMULT,
     $                  TOLH,NMXi,MATMOD)
      ELSE
         CALL hmzpf2 ('VELX',OT1,IN1,H1,H2,V1MASK,VMULT,
     $                                   IMESH,TOLH,NMXi,1)
         CALL hmzpf2 ('VELY',OT2,IN2,H1,H2,V2MASK,VMULT,
     $                                   IMESH,TOLH,NMXi,2)
         IF (NDIM.EQ.3) 
     $   CALL hmzpf2 ('VELZ',OT3,IN3,H1,H2,V3MASK,VMULT,
     $                                   IMESH,TOLH,NMXi,3)
      ENDIF
C
      return
      end
c-----------------------------------------------------------------------
      subroutine hmzpf2(nm,u,rhs,h1,h2,mask,mult,imsh,tol,mxit,isd)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      character*4 nm
c
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           RHS  (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)
C
      ntot1 = nx1*ny1*nz1*nelv
      if (imsh.eq.2) ntot1 = nx1*ny1*nz1*nelt
c
      call col2   (rhs,mask,ntot1)
      call dssum  (rhs,nx1,ny1,nz1)
      call projh2 (rhs,h1,h2,mult,mask,isd)
      call hmhzpf (nm,u,rhs,h1,h2,mask,mult,imsh,tol,mxit,isd,binvm1)
      call gensh2 (u,h1,h2,mult,mask,isd)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine projh2(v1,h1,h2,vml,vmask,isd)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOV'
C
      real v1(1),h1(1),h2(1),vml(1),vmask(1)
      real work(mxprev)
C
      integer icalld
      save    icalld
      data    icalld/0/
C
      ntot1=nx1*ny1*nz1*nelv
      n1 = 1
      if (icalld.eq.0) then
C        First call, we have no vectors to orthogonalize against.
         do j=1,ndim
            Nprev(j)=0
         enddo
         Mprev=param(93)
         Mprev=min(Mprev,Mxprev)
         if (mprev.eq.0) mprev = mxprev
         if (nid.eq.0) write(6,*) 'this is mprev:',mprev,mxprev
      endif
C
C     Diag to see how much reduction in the residual is attained.
C
      ALPHA1 = glsc23(v1,binvm1,vml,ntot1)
      ALPHA1 = sqrt(alpha1/volvm1)
      if (icalld.eq.0.and.nid.eq.0) 
     $     write(6,*) 'alpha1:',alpha1,volvm1,ntot1
      if (icalld.eq.0.and.nid.eq.0) 
     $     write(6,*) 'binvm1:',binvm1(1,1,1,1),vml(1),v1(1)
C
C     Update rhs's if matrix has changed
C
      call updrhsh2(h1,h2,vml,vmask,isd)
C
C     Perform Gram-Schmidt for previous soln's
C
      call rzero(alpha,mxprev)
      ioff = 1
      do i=1,Nprev(isd)
         ALPHA(i) = vlsc3(v1,sln(ioff,isd),vml,NTOT1)
         ioff = ioff + ntot1
      enddo
      IF (Nprev(isd).GT.0) then
         call gop(alpha,WORK,'+  ',Nprev(isd))
         call cmult2(vbar(n1,isd),sln(1,isd),alpha(1),ntot1)
C
         do i=2,Nprev(isd)
            ioff = ntot1*(i-1)+1
            call add2s2(vbar(n1,isd),sln(ioff,isd),alpha(i),ntot1)
         enddo
c        alphmn = glmin(vbar(n1,isd),ntot1)
c        alphmx = glmax(vbar(n1,isd),ntot1)
         call axhelm  (vnew(n1,isd),vbar(n1,isd),H1,H2,1,1)
c        alp1mn = glmin(vnew(n1,isd),ntot1)
c        alp1mx = glmax(vnew(n1,isd),ntot1)
         call col2    (vnew(n1,isd),vmask,ntot1)
c        alp2mn = glmin(vnew(n1,isd),ntot1)
c        alp2mx = glmax(vnew(n1,isd),ntot1)
         call dssum   (vnew(n1,isd),nx1,ny1,nz1)
c        alp3mn = glmin(vnew(n1,isd),ntot1)
c        alp3mx = glmax(vnew(n1,isd),ntot1)
         call sub2    (v1,vnew(n1,isd),ntot1)
      endif
c     if (nid.eq.0) write(6,90) istep,alphmn,alphmx
c    $              ,alp1mn,alp1mx,alp2mn,alp2mx,alp3mn,alp3mx
   90 format(i4,1p8e11.3,' xx')
C ................................................................
C     Diag.
      ALPHA2 = glsc23(v1,binvm1,vml,ntot1)
      ALPHA2 = sqrt(alpha2/volvm1)
      alphmn = glmax(v1,ntot1)
      alphmx = glmin(v1,ntot1)
      ratio  = alpha1/alpha2
      n10=min(10,Nprev(isd))
      if (nid.eq.0) write(6,10) istep,alpha1,alpha2,ratio,Nprev(isd)
c     if (nid.eq.0) write(6,10) istep,alphmn,alphmx,ratio,Nprev(isd)
   10 FORMAT(I6,1p3e12.4,i6,' alph1x')
      IF (NID.EQ.0) WRITE(6,11) ISTEP,Nprev(isd),(ALPHA(I),I=1,n10)
   11 FORMAT(I6,' halpha',i4,10(1p10e12.4,/,17x))
C     Diag.
C ................................................................
c
      icalld=icalld+1
      return
      end
c-----------------------------------------------------------------------
      subroutine gensh2(v1,h1,h2,vml,vmask,isd)
C
C     Reconstruct the solution to the original problem by adding back
C     the previous solutions
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOV'
      real v1(1),h1(1),h2(1),vml(1),vmask(1)
C
      ntot1=nx1*ny1*nz1*nelv
      n1 = 1
C
C     First, save current solution
C
      call copy (vnew(n1,isd),v1,ntot1)
C
C     Reconstruct solution
C
      call add2(v1,vbar(n1,isd),ntot1)
C
C     Update the set of <sln>
C
      call updtseth2(v1,h1,h2,vml,vmask,isd)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine updtseth2(v1,h1,h2,vml,vmask,isd)
C
C     Update the set of rhs's and the corresponding p-set:
C
C        . Standard case is to add P_new, and RHS_new = E*P_new
C
C        . However, when Nprev=Mprev (max. allowed), we throw out
C          the old set, and set P_1 = P, RHS_1=E*P_1
C
C        . Other schemes are possible, e.g., let's save a bunch of
C          old vectors, perhaps chosen wisely via P.O.D.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOV'
      real v1(1),h1(1),h2(1),vml(1),vmask(1)
C
      ntot1=nx1*ny1*nz1*nelv
      n1 = 1
C     
      IF (Nprev(isd).EQ.Mprev) THEN
         CALL COPY(vnew(n1,isd),v1,ntot1)
         Nprev(isd)=0
      ENDIF
C
C     Increment solution set
      Nprev(isd) = Nprev(isd)+1
      ioff = ntot1*(Nprev(isd)-1)+1
      CALL COPY(sln(ioff,isd),vnew(n1,isd),NTOT1)
c
C     Orthogonalize rhs against previous rhs and normalize
      call hconj2(Nprev(isd),H1,H2,vml,vmask,isd)
c
c     Save last sol'n
      CALL COPY(vnew(n1,isd),v1,ntot1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hconj2(kprev,h1,h2,vml,vmask,isd)
C
C     Orthonormalize the last saved vector against vector set
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ORTHOV'
      real h1(1),h2(1),vml(1),vmask(1)
      real work(mxprev)
C
      ntot1=nx1*ny1*nz1*nelv
      n1 = 1
C
      i1 = (kprev-1)*ntot1 + 1
      call axhelm  (vbar(n1,isd),sln(i1,isd),H1,H2,1,1)
      call col2    (vbar(n1,isd),vmask,ntot1)
      call dssum   (vbar(n1,isd),nx1,ny1,nz1)
C
C     Compute part of the norm
C
      call col2     (vbar(1,isd),vml   ,ntot1)
      alphad = glsc2(vbar(1,isd),sln(i1,isd),ntot1)
C
C     Gram-Schmidt
      Kprev1=Kprev-1
      DO 10 I=1,Kprev1
         ioff = (i-1)*ntot1 + 1
         ALPHA(I) = VLSC2(vbar(1,isd),sln(ioff,isd),ntot1)
   10 CONTINUE
      IF (Kprev1.GT.0) call gop(alpha,WORK,'+  ',Kprev1)
C
      DO 20 I=1,Kprev1
         alpham = -alpha(i)
         ioff = (i-1)*ntot1 + 1
         CALL ADD2S2(sln(i1,isd),sln(ioff,isd),alpham,ntot1)
         Alphad = Alphad - alpha(i)**2
   20 CONTINUE
C
C    .Normalize new element in P~
C
      ALPHAd = 1.0/SQRT(ALPHAd)
      CALL CMULT(sln(i1,isd),alphad,ntot1)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine updrhsh2(h1,h2,vml,vmask,isd)
C
C     Update rhs's if A-matrix has changed
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'MASS'
      INCLUDE 'ORTHOV'
      INCLUDE 'TSTEP'
C
C
      REAL vml(1),h1(1),h2(1),vmask(1)
C
      real    dtold
      save    dtold
      data    dtold/0.0/
c
C     First, we have to decide if the E matrix has changed.
C
      if (dt.eq.dtold) return
      dtold = dt
      nprev(isd) = 0
      return
c
      DO 100 Iprev=1,Nprev(isd)
C        Orthogonalize this rhs w.r.t. previous rhs's
         call hconj2(Iprev,H1,H2,vml,vmask,isd)
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
