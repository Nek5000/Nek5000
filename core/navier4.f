c-----------------------------------------------------------------------
      subroutine setrhs(p,h1,h2,h2inv)
C
C     Project rhs onto best fit in the "E" norm.
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      logical ifdump
      save    ifdump
      data    ifdump /.false./
C
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev
      REAL ALPHA,WORK
C
C
      integer icalld
      save    icalld
      data    icalld/0/
C
C     First call, we have no vectors to orthogonalize against.
      IF (ICALLD.EQ.0) THEN
         icalld=icalld+1
         Nprev=0
         Mprev=param(93)
         Mprev=min(Mprev,Mxprev)
      ENDIF
C
C     Diag to see how much reduction in the residual is attained.
C
      NTOT2  = NX2*NY2*NZ2*NELV
      ALPHA1 = GLSC3(p,p,bm2inv,NTOT2)
      if (alpha1.gt.0) ALPHA1 = sqrt(alpha1/volvm2)
C
C     Update rhs's if E-matrix has changed
C
      CALL UPDRHSE(P,H1,H2,H2INV,ierr)
      if (ierr.eq.1) Nprev=0
C
C     Perform Gram-Schmidt for previous rhs's.
C
      DO 10 I=1,Nprev
         ALPHA(i) = VLSC2(P,RHS(1,i),NTOT2)
   10 CONTINUE
C
      IF (Nprev.GT.0) CALL gop(alpha,WORK,'+  ',Nprev)
C
      CALL RZERO(Pbar,NTOT2)
      DO 20 I=1,Nprev
         alphas = alpha(i)
         CALL ADD2S2(Pbar,RHS(1,i),alphas,NTOT2)
   20 CONTINUE
C
      if (Nprev.gt.0) then
         INTETYPE = 1
         CALL CDABDTP(Pnew,Pbar,H1,H2,H2INV,INTETYPE)
         CALL SUB2   (P,Pnew,NTOT2)
C    ................................................................
C      Diag.
         ALPHA2 = GLSC3(p,p,bm2inv,NTOT2)
         if (alpha2.gt.0) ALPHA2 = sqrt(alpha2/volvm2)
         ratio  = alpha1/alpha2
         n10=min(10,nprev)
         IF (NID.EQ.0) WRITE(6,11)ISTEP,Nprev,(ALPHA(I),I=1,N10)
         IF (NID.EQ.0) WRITE(6,12) ISTEP,nprev,ALPHA1,ALPHA2,ratio
   11    FORMAT(2I5,' alpha:',1p10e12.4)
   12    FORMAT(I6,i4,1p3e12.4,' alph12')
C
c        if (alpha1.gt.0.001 .and. .not.ifdump) then
c           IF (NID.EQ.0) WRITE(6,*) 'alph1 large ... '
c           if (istep.gt.10) then 
c              IF (NID.EQ.0) WRITE(6,*) ' ... dumping'
c              call prepost (.true.,'   ')
c              ifdump = .true.
c           else
c              IF (NID.EQ.0) WRITE(6,*) ' ... doing nothing'
c           endif
c        endif
c        if (alpha1.gt.0.1.and.istep.gt.10) then
c           IF (NID.EQ.0) WRITE(6,*) 'alph1 too large ... aborting'
c           call prepost (.true.,'   ')
c           call exitt
c        endif
C    ................................................................
      endif
C
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine gensoln(p,h1,h2,h2inv)
C
C     Reconstruct the solution to the original problem by adding back
C     the previous solutions
C     know the soln.
C
      include 'SIZE'
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev
      REAL ALPHA,WORK

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      NTOT2=NX2*NY2*NZ2*NELV
C
C     First, save current solution
C
      CALL COPY (Pnew,P,NTOT2)
C
C     Reconstruct solution
C
      CALL ADD2(P,Pbar,NTOT2)
C
C     Update the set of <p,rhs>
C
      CALL UPDTSET(P,H1,H2,H2INV,ierr)
      if (ierr.eq.1) Nprev = 0
c
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine updtset(p,h1,h2,h2inv,IERR)
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
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev

      REAL ALPHA,WORK

      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      NTOT2=NX2*NY2*NZ2*NELV
C
      IF (Nprev.EQ.Mprev) THEN
         CALL COPY(Pnew,P,NTOT2)
         Nprev=0
      ENDIF
C
C     Increment solution set
      Nprev = Nprev+1
C
      CALL COPY   (RHS(1,Nprev),Pnew,NTOT2)
C
C     Orthogonalize rhs against previous rhs and normalize
C
      CALL ECONJ (Nprev,H1,H2,H2INV,ierr)
c     CALL ECHECK(Nprev,H1,H2,H2INV,INTETYPE)
C
c     Save last sol'n
      CALL COPY(Pnew,P,NTOT2)
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine econj(kprev,h1,h2,h2inv,ierr)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2),Pbrr(ltot2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev
      REAL ALPHA,WORK
      real ALPHAd
C
C
      ierr  = 0
      NTOT2 = NX2*NY2*NZ2*NELV
      INTETYPE=1
C
C     Gram Schmidt, w re-orthogonalization
C
      npass=1
      if (abs(param(105)).eq.2) npass=2
      do ipass=1,npass
c
         CALL CDABDTP(Pbrr,RHS(1,Kprev),H1,H2,H2INV,INTETYPE)
C
C        Compute part of the norm
         Alphad = GLSC2(RHS(1,Kprev),Pbrr,NTOT2)
C
C        Gram-Schmidt
         Kprev1=Kprev-1
         DO 10 I=1,Kprev1
            ALPHA(I) = VLSC2(Pbrr,RHS(1,i),NTOT2)
   10    CONTINUE
         IF (Kprev1.GT.0) CALL gop(alpha,WORK,'+  ',Kprev1)
C
         DO 20 I=1,Kprev1
            alpham = -alpha(i)
            CALL ADD2S2(RHS(1,Kprev),RHS (1,i),alpham,NTOT2)
            Alphad = Alphad - alpha(i)**2
   20    CONTINUE
      enddo
C
C    .Normalize new element in P~
C
      if (ALPHAd.le.0.0) then
         write(6,*) 'ERROR:  alphad .le. 0 in ECONJ',alphad,Kprev
         ierr = 1
         return
      endif
      ALPHAd = 1.0/SQRT(ALPHAd)
      ALPHAN = Alphad
      CALL CMULT(RHS (1,Kprev),alphan,NTOT2)
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine chkptol
C--------------------------------------------------------------------
C
C     Check pressure tolerance for transient case.
C
C     pff 6/20/92
C     This routine has been modified for diagnostic purposes only.
C     It can be replaced with the standard nekton version.
C 
C--------------------------------------------------------------------
      include 'SIZE'
      include 'SOLN'
      include 'MASS'
      include 'INPUT'
      include 'TSTEP'
      COMMON /CTOLPR/ DIVEX
      COMMON /CPRINT/ IFPRINT
      LOGICAL         IFPRINT
C
      COMMON /SCRUZ/  DIVV (LX2,LY2,LZ2,LELV)
     $ ,              BDIVV(LX2,LY2,LZ2,LELV)
C
      if (ifsplit) return
      IF (param(102).eq.0.and.(TOLPDF.NE.0. .OR. ISTEP.LE.5)) RETURN
      five = 5.0
      if (param(102).ne.0.0) five=param(102)
C
      NTOT2 = NX2*NY2*NZ2*NELV
      if (ifield.eq.1) then     ! avo: sub arguments?
         CALL OPDIV (BDIVV,VX,VY,VZ)
      else
         CALL OPDIV (BDIVV,BX,BY,BZ)
      endif
      CALL COL3 (DIVV,BDIVV,BM2INV,NTOT2)
      DNORM = SQRT(GLSC2(DIVV,BDIVV,NTOT2)/VOLVM2) 
C
      if (nid.eq.0) WRITE (6,*) istep,' DNORM, DIVEX',DNORM,DIVEX
C
c     IF (istep.gt.10.and.DNORM.GT.(1.01*DIVEX).AND.DIVEX.GT.0.) then
c        if (DNORM.gt.1e-8) then
c           if (nid.eq.0) WRITE(6,*) 'DNORM-DIVEX div. ... aborting'
c           call prepost (.true.,'   ')
c           call exitt
c        else
c           if (nid.eq.0) WRITE(6,*) 'DNORM-DIVEX div. ... small'
c        endif
c     endif
C
c     IF (DNORM.GT.(1.2*DIVEX).AND.DIVEX.GT.0.) TOLPDF = 5.*DNORM
      IF (istep.gt.5.and.tolpdf.eq.0.0.and.
     $    DNORM.GT.(1.2*DIVEX).AND.DIVEX.GT.0.) 
     $     TOLPDF = FIVE*DNORM
C
      RETURN
      END
      FUNCTION VLSC3(X,Y,B,N)
C
C     local inner product, with weight
C
      DIMENSION X(1),Y(1),B(1)
      REAL DT
C
      include 'OPCTR'
C
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + dfloat(isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + dfloat(isbcnt)
C
      DT = 0.0
      DO 10 I=1,N
         T = X(I)*Y(I)*B(I)
         DT = DT+T
 10   CONTINUE
      T=DT
      VLSC3 = T
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine updrhse(p,h1,h2,h2inv,ierr)
C
C     Update rhs's if E-matrix has changed
C
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
C
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev
      COMMON /ORTHOL/ IFNEWE
      REAL ALPHA,WORK
      LOGICAL IFNEWE
C
C
      REAL             P    (LX2,LY2,LZ2,LELV)
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      integer icalld
      save    icalld
      data    icalld/0/
      NTOT2=NX2*NY2*NZ2*NELV
C
C
C     First, we have to decide if the E matrix has changed.
C
      IF (icalld.eq.0) THEN
         icalld=1
         DTlast=DT
      ENDIF
C
      IFNEWE=.FALSE.
      IF (IFMVBD) THEN
         IFNEWE=.TRUE.
         CALL INVERS2(bm2inv,bm2,Ntot2)
      ELSEIF (DTlast.ne.DT) THEN
         IFNEWE=.TRUE.
         DTlast=DT
      ENDIF
      IF (IFNEWE.and.nid.eq.0) write(6,*) 'reorthogo:',nprev
C
C     
C     
C     Next, we reconstruct a new rhs set.
C     
      IF (IFNEWE) THEN
c
c        new idea...
         if (nprev.gt.0) nprev=1
         call copy(rhs,pnew,ntot2)
c
         Nprevt = Nprev
         DO 100 Iprev=1,Nprevt
C           Orthogonalize this rhs w.r.t. previous rhs's
            CALL ECONJ (Iprev,H1,H2,H2INV,ierr)
            if (ierr.eq.1) then
               Nprev = 0
               return
            endif
  100    CONTINUE
C
      ENDIF
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine echeck(kprev,h1,h2,h2inv,intetype)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
C
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
C
      REAL             H1   (LX1,LY1,LZ1,LELV)
      REAL             H2   (LX1,LY1,LZ1,LELV)
      REAL             H2INV(LX1,LY1,LZ1,LELV)
C
      PARAMETER (LTOT2=LX2*LY2*LZ2*LELV)
      COMMON /ORTHOV/ RHS(LTOT2,MXPREV)
      COMMON /ORTHOX/ Pbar(LTOT2),Pnew(LTOT2)
      COMMON /ORTHOS/ ALPHA(Mxprev), WORK(Mxprev), ALPHAN, DTLAST
      COMMON /ORTHOI/ Nprev,Mprev
      REAL ALPHA,WORK,GLSC2
      REAL Alphad
C
C
      NTOT2=NX2*NY2*NZ2*NELV
C
C     Compute part of the norm
C
      do 20 j=1,kprev
         CALL CDABDTP(Pbar,rhs(1,j),H1,H2,H2INV,INTETYPE)
         do 10 i=1,kprev
            Alphad = GLSC2(RHS(1,i),Pbar,NTOT2)
            Alphas = alphad
            if (nid.eq.0) then
               write(6,5) i,j,alphad,alphas,istep,kprev
    5          format(' E-check:',2i4,e16.8,g12.5,i6,i4)
            endif
   10    continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
c     THE ROUTINES BELOW ARE THE NEW Helmholtz projectors
c-----------------------------------------------------------------------
      subroutine projh(r,h1,h2,bi,vml,vmk,approx,napprox,wl,ws,name4)
C
C     Orthogonalize the rhs wrt previous rhs's for which we already
C     know the soln.
c
c     Input:   r         -- residual
c              h1,h2     -- Helmholtz arrays
c              bi        -- inverse mass matrix
c              vml,vmk   -- multiplicity and mask arrays
c              approx    -- approximation space
c              napprox   -- (1) = max vecs,  (2) = current number of vecs
c              wl        -- large work array of size lx1*ly1*lz1*nelv
c              ws        -- small work array of size 2*max vecs
c
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
c
      parameter(lt=lx1*ly1*lz1*lelt)
C
      real r(1),h1(1),h2(1),vml(1),vmk(1)
      real bi(1)
      real wl(1),ws(1)
      real approx(lt,0:1)
      integer napprox(2)
      character*4 name4
c
      n_max = napprox(1)
      n_sav = napprox(2)
      if (n_sav.eq.0) return
      nel =nelfld(ifield)
      ntot=nx1*ny1*nz1*nel

      vol = voltm1
      if (nel.eq.nelv) vol = volvm1
c
c     Diag to see how much reduction in the residual is attained.
c
      alpha1 = glsc23(r,bi,vml,ntot)
      if (alpha1.gt.0) alpha1 = sqrt(alpha1/vol)
c
c     Update approximation space if dt has changed
      call updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)
c
c
c     Perform Gram-Schmidt for previous soln's
c
      do i=1,n_sav
         ws(i) = vlsc3(r,approx(1,i),vml,ntot)
      enddo
      call gop    (ws,ws(n_sav+1),'+  ',n_sav)
c
      call cmult2   (approx(1,0),approx(1,1),ws(1),ntot)
      do i=2,n_sav
         call add2s2(approx(1,0),approx(1,i),ws(i),ntot)
      enddo
c
      call axhelm  (wl,approx(1,0),h1,h2,1,1)
      call col2    (wl,vmk,ntot)
      call dssum   (wl,nx1,ny1,nz1)
      call sub2    (r ,wl,ntot)
c ................................................................
c   Diag.
      alpha2 = glsc23(r,bi,vml,ntot)
      if (alpha2.gt.0) alpha2 = sqrt(alpha2/vol)
      ratio  = alpha1/alpha2
      n10=min(10,n_sav)
c
      if (nid.eq.0) write(6,10) istep,name4,alpha1,alpha2,ratio,n_sav
   10 format(4X,I7,4x,a4,' alph1n',1p3e12.4,i6)
c
      if (nid.eq.0) write(6,11) istep,name4,n_sav,(ws(i),i=1,n10)
   11 format(4X,I7,4x,a4,' halpha',i6,10(1p10e12.4,/,17x))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gensh(v1,h1,h2,vml,vmk,approx,napprox,wl,ws,name4)
c
c     Reconstruct the solution to the original problem by adding back
c     the previous solutions
c
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
c
      common /iterhm/ niterhm
c
      parameter(lt=lx1*ly1*lz1*lelt)
      real v1(1),h1(1),h2(1),vml(1),vmk(1)
      real wl(1),ws(1)
      real approx(lt,0:1)
      integer napprox(2)
      character*4 name4
c
      n_max = napprox(1)
      n_sav = napprox(2)
      ntot=nx1*ny1*nz1*nelfld(ifield)
c
c     Reconstruct solution and save current du
c
      if (n_sav.lt.n_max) then
c
         if (niterhm.gt.0) then      ! new vector not in space
            n_sav = n_sav+1
            call copy(approx(1,n_sav),v1,ntot)
            call add2(v1,approx(1,0),ntot)
c           orthogonalize rhs against previous rhs and normalize
            call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)

c           if (ierr.ne.0) n_sav = n_sav-1
            if (ierr.ne.0) n_sav = 0

         else

            call add2(v1,approx(1,0),ntot)

         endif
      else
         n_sav = 1
         call add2(v1,approx(1,0),ntot)
         call copy(approx(1,n_sav),v1,ntot)
c        normalize
         call hconj(approx,n_sav,h1,h2,vml,vmk,ws,name4,ierr)
         if (ierr.ne.0) n_sav = 0
      endif

      napprox(2)=n_sav

      return
      end
c-----------------------------------------------------------------------
      subroutine hconj(approx,k,h1,h2,vml,vmk,ws,name4,ierr)
c
c     Orthonormalize the kth vector against vector set
c
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'PARALLEL'
      include 'TSTEP'
c
      parameter  (lt=lx1*ly1*lz1*lelt)
      real approx(lt,0:1),h1(1),h2(1),vml(1),vmk(1),ws(1)
      character*4 name4
c
      ierr=0
      ntot=nx1*ny1*nz1*nelfld(ifield)
c
      call axhelm  (approx(1,0),approx(1,k),h1,h2,1,1)
      call col2    (approx(1,0),vmk,ntot)
      call dssum   (approx(1,0),nx1,ny1,nz1)
      call col2    (approx(1,0),vml        ,ntot)
c
c     Compute part of the norm   (Note:  a(0) already scaled by vml)
c
      alpha = glsc2(approx(1,0),approx(1,k),ntot)
      alph1 = alpha
c
c     Gram-Schmidt
c
      km1=k-1
      do i=1,km1
         ws(i) = vlsc2(approx(1,0),approx(1,i),ntot)
      enddo
      if (km1.GT.0) call gop(ws,ws(k),'+  ',km1)
c
      do i=1,km1
         alpham = -ws(i)
         call add2s2(approx(1,k),approx(1,i),alpham,ntot)
         alpha = alpha - ws(i)**2
      enddo
c
c    .Normalize new element in approximation space
c
      eps = 1.e-7
      if (wdsize.eq.8) eps = 1.e-15
      ratio = alpha/alph1
c
      if (ratio.le.0) then
         ierr=1
         if (nid.eq.0) write(6,12) istep,name4,k,alpha,alph1
   12    format(I6,1x,a4,' alpha b4 sqrt:',i4,1p2e12.4)
      elseif (ratio.le.eps) then
         ierr=2
         if (nid.eq.0) write(6,12) istep,name4,k,alpha,alph1
      else
         ierr=0
         alpha = 1.0/sqrt(alpha)
         call cmult(approx(1,k),alpha,ntot)
      endif
c
      if (ierr.ne.0) then
         call axhelm  (approx(1,0),approx(1,k),h1,h2,1,1)
         call col2    (approx(1,0),vmk,ntot)
         call dssum   (approx(1,0),nx1,ny1,nz1)
         call col2    (approx(1,0),vml        ,ntot)
c
c        Compute part of the norm   (Note:  a(0) already scaled by vml)
c
         alpha = glsc2(approx(1,0),approx(1,k),ntot)
         if (nid.eq.0) write(6,12) istep,name4,k,alpha,alph1
         if (alpha.le.0) then
            ierr=3
            if (nid.eq.0) write(6,12) istep,name4,k,alpha,alph1
            return
         endif
         alpha = 1.0/sqrt(alpha)
         call cmult(approx(1,k),alpha,ntot)
         ierr = 0
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine updrhsh(approx,napprox,h1,h2,vml,vmk,ws,name4)
c
c     Reorthogonalize approx if dt has changed
c
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
c
c
      parameter  (lt=lx1*ly1*lz1*lelt)
      real approx(lt,0:1),h1(1),h2(1),vml(1),vmk(1),ws(1)
      integer napprox(2)
      character*4 name4
c
      logical ifupdate
      logical ifnewdt
      save    ifnewdt
      data    ifnewdt /.false./
c
      character*4 name_old
      save    name_old
      data    name_old/'DMIR'/
c
      real    dtold
      save    dtold
      data    dtold/0.0/
c
c     First, we have to decide if the dt has changed.
c
      ifupdate = .false.
      if (dt.ne.dtold) then
         dtold    = dt
         name_old = name4
         ifnewdt  = .true.
         ifupdate = .true.
      elseif (ifnewdt) then
         if (name4.eq.name_old) then
            ifnewdt = .false.
         else
            ifupdate = .true.
         endif
      endif
      if (ifvarp(ifield)) ifupdate = .true.
      if (iflomach)       ifupdate = .true.

      if (ifupdate) then    ! reorthogonalize 
         n_sav = napprox(2)
         l     = 1
         do k=1,n_sav
c           Orthogonalize kth vector against {v_1,...,v_k-1}
            if (k.ne.l) then
               ntot = nx1*ny1*nz1*nelfld(ifield)
               call copy(approx(1,l),approx(1,k),ntot)
            endif
            call hconj(approx,l,h1,h2,vml,vmk,ws,name4,ierr)
            if (ierr.eq.0) l=l+1
         enddo
         napprox(2)=min(l,n_sav)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hmhzpf(name,u,r,h1,h2,mask,mult,imesh,tli,maxit,isd,bi)
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'FDMH1'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           MASK (LX1,LY1,LZ1,1)
      REAL           MULT (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      COMMON /CTMP0/ W1   (LX1,LY1,LZ1,LELT)
     $ ,             W2   (LX1,LY1,LZ1,LELT)
c
      etime1=dnekclock()
c
      IF (IMESH.EQ.1) NTOT = NX1*NY1*NZ1*NELV
      IF (IMESH.EQ.2) NTOT = NX1*NY1*NZ1*NELT
c
      tol = tli
      if (param(22).ne.0) tol = abs(param(22))
      CALL CHKTCG1 (TOL,R,H1,H2,MASK,MULT,IMESH,ISD)
c
c
c     Set flags for overlapping Schwarz preconditioner (pff 11/12/98)
c
                          kfldfdm = -1
c     if (name.eq.'TEMP') kfldfdm =  0
c     if (name.eq.'VELX') kfldfdm =  1
c     if (name.eq.'VELY') kfldfdm =  2
c     if (name.eq.'VELZ') kfldfdm =  3
      if (name.eq.'PRES') kfldfdm =  ndim+1
c
      call cggo
     $      (u,r,h1,h2,mask,mult,imesh,tol,maxit,isd,bi,name)
      thmhz=thmhz+(dnekclock()-etime1)
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hsolve(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd
     $                 ,approx,napprox,bi)
c
c     Either std. Helmholtz solve, or a projection + Helmholtz solve
c
      include 'SIZE'
      include 'INPUT'
      include 'MASS'
      include 'TSTEP'
      include 'FDMH1'
      include 'CTIMER'
c
      CHARACTER*4    NAME
      REAL           U    (LX1,LY1,LZ1,1)
      REAL           R    (LX1,LY1,LZ1,1)
      REAL           H1   (LX1,LY1,LZ1,1)
      REAL           H2   (LX1,LY1,LZ1,1)
      REAL           vmk  (LX1,LY1,LZ1,1)
      REAL           vml  (LX1,LY1,LZ1,1)
      REAL           bi   (LX1,LY1,LZ1,1)
      REAL           approx (1)
      integer        napprox(1)
      common /ctmp2/ w1   (lx1,ly1,lz1,lelt)
      common /ctmp3/ w2   (2+2*mxprev)

      logical ifstdh
      character*4  cname

      call chcopy(cname,name,4)
      call capit (cname,4)

      ifstdh = .true.

      if (.not.ifflow) ifstdh = .false.

      if (param(95).ne.0.and.istep.gt.param(95)) then
         if (cname.eq.'PRES') ifstdh = .false.
      elseif (param(94).ne.0.and.istep.gt.param(94)) then
         ifstdh = .false.
      endif

      if (param(93).eq.0) ifstdh = .true.

      if (ifstdh) then
         call hmholtz(name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd)
      else

         n = nx1*ny1*nz1*nelfld(ifield)

         call dssum  (r,nx1,ny1,nz1)
         call col2   (r,vmk,n)
         call projh  (r,h1,h2,bi,vml,vmk,approx,napprox,w1,w2,name)
         call hmhzpf (name,u,r,h1,h2,vmk,vml,imsh,tol,maxit,isd,bi)
         call gensh  (u,h1,h2,vml,vmk,approx,napprox,w1,w2,name)

      endif

      return
      end
c-----------------------------------------------------------------------
