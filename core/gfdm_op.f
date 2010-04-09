c-----------------------------------------------------------------------
      subroutine gfdm_ops
c
c     Initialize Fast Diagonalization Method for x-y-z tensor product 
c     geometry
c
c     Assumes that  p116 = Nelx, p117=Nely, p118=Nelz
c
c     Also, assumes that elements are *ordered lexicographically*
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'ZPER'

      call gfdm_chk_size

      mflds = 1
      if (ifmhd) mflds = 2
c
      ii = 1          ! we'll need to change this for multiple flds
c
      do mfld=1,mflds
c
        call gfdm_set_prs_op(mfld)
        l = ngfdm_p(1)
        m = ngfdm_p(2)
        n = ngfdm_p(3)
c
        i = mlp(1,mfld)
        j = mlp(2,mfld)
        k = mlp(3,mfld)
c
c       write(6,*) 'mfld:',mfld,l,m,n,i,j,k,eigp(1)
c
        call gfdm_set_diagp
     $   (wavep(ii),tpn1,mcex,eigp(i),l,eigp(j),m,eigp(k),n)
c
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_prs_op(mfld)
c
c     Set up global 1D pressure operators for general penciled 
c     partition of 2D or 3D tensor-product domain, to be used
c     with fast-diagonalization method solver
c
c     Based on exact E operator, w/ correct velocity bc's in each direction.
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
c
      common /fastdr/ lx(lelg_sm),ly(lelg_sm),lz(lelg_sm)
      real            lx,      ly,      lz
c
      parameter (lbw=7*lx1*ly1*lz1*lelv-1)
      common /SCRNS/ bw(0:lbw)
c
      character*1  cb0(3),cbn(3)
c
c     Set up array of element pointers
c
      nelq = nelx*nely*nelz
      if (nelq.ne.nelgv.or.nelgv.gt.lelg_sm) then
         write(6,1) nid,nelq,nelv,nelx,nely,nelz,nelgv,lelg_sm
    1    format(i8,' problem in set_fast_eprec, nelv?',7i8)
         call exitt
      endif
c
      if(nelq.gt.lbw) then
         write(6,2) nid,nelq,nelv,nelx,nely,nelz,lbw
    2    format(i8,' problem in set_fast_eprec, lbw?',7i8)
         call exitt
      endif
c
      call gfdm_set_bc    (cb0,cbn,mfld)        !  Find tp-box bdry conds.
      call gfdm_set_geom  (bw,nelx,nely,nelz)   !  Get 1D elem't lengths
      call gfdm_set_genwz (nx1,nx2)             !  Compute weights in WFZ
c
c
c     Generate eigenvector and eigenvalue arrays for 1D E-operators
c
      l = 1             !   X-direction
      m = 1
      if (mfld.eq.2) then
         do k=1,ndim
            l = l+ngfdm_p(k)
            m = m+ngfdm_p(k)**2
         enddo
      endif
c
      msp(1,mfld) = m
      mlp(1,mfld) = l
c
      call set_1d_e_mat(sp(m),eigp(l),mx,lx,nelx,nx1,nx2
     $                 ,cb0(1),cbn(1),dglgt,wglgt,wgl,bw,lbw)
      call transpose(spt(m),mx,sp(m),mx)
c
      l = l+mx          !   Y-direction
      m = m+mx*mx
      msp(2,mfld) = m
      mlp(2,mfld) = l
      call set_1d_e_mat(sp(m),eigp(l),my,ly,nely,nx1,nx2
     $                 ,cb0(2),cbn(2),dglgt,wglgt,wgl,bw,lbw)
      call transpose(spt(m),my,sp(m),my)
c
      if (if3d) then
         l = l+my       !   Z-direction
         m = m+my*my
         msp(3,mfld) = m
         mlp(3,mfld) = l
         call set_1d_e_mat(sp(m),eigp(l),mz,lz,nelz,nx1,nx2
     $                    ,cb0(3),cbn(3),dglgt,wglgt,wgl,bw,lbw)
         call transpose(spt(m),mz,sp(m),mz)
      else
         mz=1
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_1d_e_mat(eigv,eigs,n,length,nel,nxv,nxp
     $                          ,cb0,cbn,dglgt,wglgt,wgl,bw,lbw)
c
c     Construct and diagonalize 1D E-matrix (pressure op)
c
c
      real    eigv(1),eigs(1)     ! output eigenvec and eigenvalues
      integer n                   ! output size of eigv (nxn) and eigs(n)
c
      real    length(nel)         ! input array containing length of els.
      integer nxv,nxp             ! # els, # vel pts/el, # pres pts/el.
      character*1 cb0,cbn         ! bc at each end of domain
      real dglgt(1),wglgt(1),wgll ! 1D deriv and interp ops (nxp to nxv)
      real bw(0:lbw)                ! big work array
c
c
c= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c     Build  one-dimensional E operator
c= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
      n   = nel*nxp
      n1  = nel*(nxv-1)
      nn  = n*n
      n2  = n*(n1+5)
      nbw = nn+n2+n*n1
      if (nbw.gt.lbw) then
         write(6,*) 
     $   'ABORT. Insufficient space in set_1d_e_mat(bw).',nbw,lbw
         call exitt
      endif
c
      call build_D1d_d(eigv,bw,dglgt,wglgt,wgl,eigs,length,nel,nxv,nxp
     $               ,bw(nn),bw(nn+n2),cb0,cbn)
c
c     Compute eigenvalues and eigenvectors of 1D operator
c     call outmat(eigv,n,n,'E1D   ',nn)
      call generalev  (eigv,bw,eigs,n,bw(nn))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine build_D1d_d(DI,II,dgg1,igg1,wm1,b1d,length,nel,nxv,nxp
     $                      ,wi,ti,cb0,cbn)
c
c     Build  1-dimensional multi-element D   matrix
c                                         I
c
c
c            T         -1    T                             
c     D  := D    Q    B     Q    D
c      I     GG   I    I     I    GG                       
c
c
c
c            T         -1    T                             
c     I  := I    Q    B     Q    I
c      I     GG   I    I     I    GG                       
c
c
c
c    Here,
c
c
c       
c    I   := K1-block-diagonal matrix of interpolation operators from M2 to M1.
c     GG
c
c       
c    D   := K1-block-diagonal matrix of derivative operators from M2 to M1.
c     GG
c
c     T 
c    Q   := K1 one-dimensional spectral element direct stiffness summation op.
c     I 
c
c                                                             T
c    B   := K1-block-diagonal (dssumed) 1D mass-matrix (B  = Q  B  Q   )
c     I                                                  I    I  L  I
c
c
c     Return the matrices, I  and D  , ready for a gen'd eigenv. solver.
c                           I      I
c
      real di(1),ii(1),wi(1),ti(1),dgg1(1),igg1(1),wm1(1),b1d(1)
      real length(nel)
      character*1 cb0,cbn
c
c
c     Number of active points on mesh 1 depends on the boundary
c     conditions at z0 and zn:
c
c     If periodic:                                   n1 = nel*(nxv-1)  := n1p
c 
c     If Dirichlet at both ends,                     n1 = n1p - 1
c     (However, we retain both, and use mask...)
c 
c     If Neumann (outflow) at both ends,             n1 = n1p + 1
c 
c     If Dirichlet at one end outflow at other,      n1 = n1p
c
c
c *** NOTE that for symmetric boundary conditions ('S' := SYM ) the number 
c     of degrees-of-freedom for the velocity field in the D_GG term is n1p-1,
c     while the number of degrees-of-freedom in the velocity field for the 
c     I_GG term is n1p+1.   This is because the I_GG term operates on the
c     other *velocity* components (those tangential to the symmetry plane)
c     for which Neumann boundary conditions are appropriate.
c
      n  = nel*nxp
      n1 = nel*(nxv-1) + 1
      if (cb0.eq.'P') n1 = n1-1
c
c
c     Perform direct summation of 1D mass matrices
c
      i = 1
      call rzero(b1d,n1+1)
      do ie=1,nel
c        build B
         s = 0.5*length(ie)
         call add2s2(b1d(i),wm1,s,nxv)
         i = i + (nxv-1)
      enddo
      if (cb0.eq.'P') b1d( 1) = b1d( 1)+b1d(n1+1)
      call invcol1(b1d,n1)
c
c     Mask velocity mesh at endpoints of 1D operator
c
      if (cb0.eq.'D') b1d( 1) = 0.
      if (cbn.eq.'D') b1d(n1) = 0.
c
c
c     Perform direct summation of 1D I    matrices
c                                     GG
      i = 1
      j = 1
      call rzero(wi,n*n1)
      do ie=1,nel
         s = 0.5*length(ie)
         call add2s2mat2p(WI,n1,n,i,j,igg1,nxv,s,nxv,nxp)
         i = i + (nxv-1)
         j = j + (nxp)
      enddo
c                 -1    T
c     Set ti = ( B   W )    
c                  I  I
      do i=1,n1
      do j=1,n
         ij = i + n1*(j-1)
         ji = j + n*(i-1)
         ti( ji ) = b1d(i)*wi(ij)
      enddo
      enddo
      call mxm(ti,n,wi,n1,ii,n)
c
c
c     Mask velocity mesh at endpoints of 1D operator for SYM+E case
c
      if (cb0.eq.'S') b1d( 1) = 0.
      if (cbn.eq.'S') b1d(n1) = 0.
c
c
c                   T
c     Compute w  = Q  D      (Note:  w   is (n1 x n) ,  n := (nxp)*nel )
c              I    I  GG             I
c
c
c
c
c     Perform direct summation of 1D D    matrices
c                                     GG
c
      i = 1
      j = 1
      call rzero(wi,n*n1)
      do ie=1,nel
         s  = 1.0
         call add2s2mat2p(WI,n1,n,i,j,dgg1,nxv,s,nxv,nxp)
         i = i + (nxv-1)
         j = j +  nxp
      enddo
c
c                 -1    T
c     Set ti = ( B   W )    
c                  I  I
      do i=1,n1
      do j=1,n
         ij = i + n1*(j-1)
         ji = j + n*(i-1)
         ti( ji ) = b1d(i)*wi(ij)
      enddo
      enddo
      call mxm(ti,n,wi,n1,di,n)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_diagp(eigi,tpn,nn,eigx,l,eigy,m,eigz,n)
c
c     Set up diagonal for pressure operator inversion
c     via Fast Diag. Meth.
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
c
      real    eigi(nn)
      integer tpn(nn)
      real    eigx(l),eigy(m),eigz(n)
c
      eps = 1.e-12
      if (wdsize.eq.4) eps = 1.e-6
      ni  = l
      nj  = m
      nk  = n
      nij = ni*nj
c
      eigxmin = glmin(eigx,l)
      eigymin = glmin(eigy,m)
      eigzmin = glmin(eigz,n)
      if (nid.eq.0) write(6,11) eigxmin,eigymin,eigzmin,l,m,n
      eigxmax = glmax(eigx,l)
      eigymax = glmax(eigy,m)
      eigzmax = glmax(eigz,n)
      if (nid.eq.0) write(6,12) eigxmax,eigymax,eigzmax,l,m,n
   11 format(1p3e12.4,' gfdm eigmins',3i6)
   12 format(1p3e12.4,' gfdm eigmaxs',3i6)
      epx = eps*abs(eigxmax)
      epy = eps*abs(eigymax)
      epz = eps*abs(eigzmax)
      if (if3d) then
         do ii=1,nn
            ijk = tpn(ii)
            i   = mod1(ijk,ni)
            k   = 1+(ijk-1)/nij
            j   = 1+(ijk-1)/ni
            j   = mod1(j,nj)
            if (eigx(i).lt.epx.and.eigy(j).lt.epy.and.eigz(k).lt.epz) 
     $                                                           then
               eigi(ii) = 0.0
               write(6,3) i,j,k,eigx(i),eigy(j),eigz(k),eps
   3           format(3i5,1p4e10.2,' zero eigenvalue, diagp')
            else
               eig3d     = ( eigx(i) + eigy(j) + eigz(k) )
               eigi(ii) = 1./eig3d
            endif
         enddo
      else
         do ii=1,nn
            ij = tpn(ii)
            i  = mod1(ij,ni)
            j  = 1+(ij-1)/ni
c           write(6,4) ij,i,j,ii,eigx(i),eigy(j),eps
c  4        format(4i8,1p3e10.2,' EIG, diagp')
            if (eigx(i).lt.eps.and.eigy(j).lt.eps) then
               eigi(ii) = 0.0
               write(6,2) i,j,eigx(i),eigy(j),eps
   2           format(2i5,1p3e10.2,' zero eigenvalue, diagp')
            else
               eig2d     = ( eigx(i) + eigy(j) )
               eigi(ii) = 1./eig2d
            endif
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_geom(work,melx,mely,melz)
c
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'ZPER'
      real work(nelx,nely,nelz)
c
      common /fastdr/ lx(lelg_sm),ly(lelg_sm),lz(lelg_sm)
      real            lx,      ly,      lz
      integer e,ex,ey,ez,eg
c
      call rzero(lx,nelgt)
      call rzero(ly,nelgt)
      call rzero(lz,nelgt)
      do e=1,nelv
         eg = lglel(e)
         lx(eg) = xm1(nx1,1,1,e) - xm1(1,1,1,e)
         ly(eg) = ym1(1,ny1,1,e) - ym1(1,1,1,e)
         lz(eg) = zm1(1,1,nz1,e) - zm1(1,1,1,e)
      enddo
c
      call gop(lx,work,'+  ',nelgt)
      call gop(ly,work,'+  ',nelgt)
      call gop(lz,work,'+  ',nelgt)
c
      call copy(work,ly,nelgt)
      do ey=1,nely
         ly(ey) = work(1,ey,1)
      enddo
c
      call copy(work,lz,nelgt)
      do ez=1,nelz
         lz(ez) = work(1,1,ez)
      enddo
c
c     Store xgtp,ygtp, and zgtp (gtp = global tensor product)
c
      n = nx1*ny1*nz1*nelv
      xgtp(0) = glmin(xm1,n)
      do ex=1,nelx
         xgtp(ex) = xgtp(ex-1) + lx(ex)
      enddo
c
      ygtp(0) = glmin(ym1,n)
      do ey=1,nely
         ygtp(ey) = ygtp(ey-1) + ly(ey)
      enddo
c
      zgtp(0) = glmin(zm1,n)
      do ez=1,nelz
         zgtp(ez) = zgtp(ez-1) + lz(ez)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_genwz(nx,nxp)
C
C     GENERATE
C
C            - DERIVATIVE OPERATORS
C            - INTERPOLATION OPERATORS
C            - WEIGHTS
C            - COLLOCATION POINTS
C
C     ASSOCIATED WITH THE
C
C            - GAUSS-LOBATTO LEGENDRE MESH (SUFFIX M1/M2/M3)
C            - GAUSS LEGENDRE         MESH (SUFFIX M2)
C            - GAUSS-LOBATTO JACOBI   MESH (SUFFIX M1/M2/M3)
C
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'WZF'
c
      parameter (lt=lx1*ly1*lz1*lelv)
      common /SCRNS/ bw(0:7*lt-1)
C
C
C
C     Gauss-Lobatto Legendre mesh (suffix M1)
C     Generate collocation points and weights
c
c
      call zwgll (zgl , wgl , nx  )
c
c     Pressure mesh weights and points
c
      if(ifsplit)then
         call zwgll (zgp , wgp , nxp )
      else
         call zwgl  (zgp , wgp , nxp )
      endif
c
      call iglm  (iggl,igglt,zgp,zgl,nxp,nx,nxp,nx)
      call igllm (wglg,wglgt,zgl,zgp,nx,nxp,nx,nxp)
      call dgllgl(dglg,dglgt,zgl,zgp,wglg,nx,nxp,nx,nxp)
c
c
c     perform row-scaling of dglg and wglg
c
      call row_mult (dglg,wgp,nxp,nx)
      call transpose(dglgt,nx,dglg,nxp)
c
      call row_mult (wglg,wgp,nxp,nx)
      call transpose(wglgt,nx,wglg,nxp)
c
c
c     Generate 1st derivative matrices
      call dgll (d1,d1t,zgl,nx,nx)
c
c
c     Generate 2nd-order derivative matrices 
      call A1D (B1iA1,B1iA1t,d2,b2p,d1,d1t,wgl,nx)
c
c
c     Compute diagonal of Laplacian 
      call set_diaga(da,dat,wgl,d1,bw,nx)
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_bc(cb0,cbn,mfld)
c
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
c
      character*1 cb0(3),cbn(3)
      integer     icb(3,2),iwk(3,2)
      character*3 c
c
      character*1 cbi(6)
      save        cbi
      data        cbi   / 'N','D','S','P','n','V' /
c
      if (ifheat) ifld=2
      if (ifflow) ifld=1
      if (mfld.eq.2) ifld=ifldmhd
      nface = 2*ndim
c
      call izero(icb,6)
c
      ieg=1                              ! Find bc's on element # 1
      jnid = gllnid(ieg)
      if (jnid.eq.nid) then
         ie = gllel(ieg)
         do jf=1,ndim
            iface = eface(2*jf-1)        ! faces 1,3,5 (pf notation)
            c   = cbc(iface,ie,ifld)
            icb(jf,1) = 1                !cb = 'N'
            if (c.eq.'v  '.or.c.eq.'V  ' .or.c.eq.'W  ') icb(jf,1)=2  ! cb0='D'
            if (c.eq.'SYM'                             ) icb(jf,1)=3  ! cb0='S'
            if (c.eq.'P  '                             ) icb(jf,1)=4  ! cb0='P'
            if (c.eq.'DVN'                             ) icb(jf,1)=5  ! cb0='P'
            if (c.eq.'DIV'                             ) icb(jf,1)=6  ! cb0='P'
         enddo
      endif
c
      ieg=nelgv                          ! Find bc's on element # NELGV
      jnid = gllnid(ieg)
      if (jnid.eq.nid) then
         ie = gllel(ieg)
         do jf=1,ndim
            iface = eface(2*jf)          ! faces 2,4,6 (pf notation)
            c   = cbc(iface,ie,ifld)
            icb(jf,2) = 1                !cb = 'N'
            if (c.eq.'v  '.or.c.eq.'V  ' .or.c.eq.'W  ') icb(jf,2)=2  ! cb0='D'
            if (c.eq.'SYM'                             ) icb(jf,2)=3  ! cb0='S'
            if (c.eq.'P  '                             ) icb(jf,2)=4  ! cb0='P'
            if (c.eq.'DVN'                             ) icb(jf,2)=5  ! cb0='P'
            if (c.eq.'DIV'                             ) icb(jf,2)=6  ! cb0='P'
         enddo
      endif
c
      call igop(icb,iwk,'+  ',6)         ! Communicate to all processors
c
c
c     Convert back to character strings
c
      do i=1,ndim
         cb0(i) = cbi(icb(i,1))
         cbn(i) = cbi(icb(i,2))
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine add2s2mat2p(a,ma,na,i1,j1,b,ldb,s,m,n)
c
c     Add two sub-matrices with periodic wrap, scaling the 2nd matrix
c
      real a(ma,na),b(ldb,1)
c
      do j=1,n
         ja     = mod1(j1+j-1,na)
         do i=1,m
            ia       = mod1(i1+i-1,ma)
            a(ia,ja) = a(ia,ja) + s*b(i,j)
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine overflow_chk(n_req,n_avail,var,sub)
c
c     Check for buffer overflow
c
      character*3 var
      character*6 sub
c
      nid = mynode()
c
      if (n_req.gt.n_avail) then
         write(6,9) nid,n_req,n_avail,nid,var,sub
    9    format(i7,' ERROR: requested array space (',i9
     $            ,') exceeds allocated amount (',i9,').'
     $            ,/,i7,' ABORTING.',3x,a3,2x,a6,' overflow_chk')
         call exitt
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine solveMp(z,r,n,w,nza)
c
c     Preconditioner
c
      real z(n,nza),r(n,nza),w(n,nza)
      call copy(z,r,n*nza)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine map12q(r2,r1,rt)
c
      include 'SIZE'
      include 'IXYZ'
c
      real r2(lx2,ly2),r1(lx1,ly1),rt(lx2,lx1)
c
      call mxm(ixm12 ,nx2,r1,nx1,rt,nx1)
      call mxm(rt,nx2,ixm12t,nx1,r2,nx2)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine cgpa(x,b,r,p,z,w,niter,tolin)
      return
      end
c-----------------------------------------------------------------------
      subroutine row_mult (A,B,n1,n2)
      real a(n1,n2),b(n1)
c
      do i=1,n1
         c = b(i)
         do j=1,n2
            a(i,j) = c*a(i,j)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine a1d (A1,A1t,d2,b2,d1,d1t,wgl,nx)
c
c     Generate 1D stiffness matrix, nodal basis
c
      real A1(nx,nx),A1t(nx,nx)
      real d2(nx,nx),b2 (nx,nx)
      real d1(nx,nx),d1t(nx,nx)
      real wgl(nx)
c
      do i=1,nx
         do j=1,nx
            b2(i,j) = wgl(i)*d1(i,j)
         enddo
      enddo
      call mxm(d1t,nx,b2,nx,d2,nx)
c
c     Make nx x nx mass matrix
      call rzero(b2,nx*nx)
      do i=1,nx
         b2(i,i) = wgl(i)
      enddo
c
c              -1
c     Compute B   A, for tensor product evaluation w/ B factored out.
c
      do i=1,nx
         s = 1./wgl(i)
         do j=1,nx
            A1(i,j) = s*d2(i,j)
         enddo
      enddo
      call transpose(A1t,nx,A1,nx)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_diagA (da,dat,b,d,w,nx)
c
c     Generate diagonal of discrete Laplacian
c
      real da(nx,nx),dat(nx,nx),b(nx),d(nx,nx),w(nx)
c
      do i=1,nx
         w(i) = 0
         do k=1,nx
            w(i) = w(i) + b(k)*d(k,i)*d(k,i)
         enddo
      enddo
c
c     Recall:   A2d = By X Ax + Ay X Bx
c
      do j=1,nx
         do i=1,nx
            da(i,j) = b(j)*w(i)
         enddo
      enddo
      call transpose(dat,nx,da,nx)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_chk_size
c
c     Check arrays for fast tensor product solver
c
c     Assumes p116 = Nelx, p117=Nely, p118=Nelz
c         and elements are *ordered lexicographically*
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'ZPER'

      if (nelx.gt.lelx .or.
     $    nely.gt.lely .or.
     $    nelz.gt.lelz ) then
         if (nid.eq.0) write(6,1) nelx,nely,nelz,lelx,lely,lelz
    1    format('gfdm_chk_size fail:',6i6)
         if (nid.eq.0) write(6,*) 'Change lelx...lelz in SIZEu.'
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
