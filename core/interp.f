C=======================================================================

      subroutine interpl_gen(ust,xst,u,cat)

      include 'SIZE'
      include 'INPUT'
      real ust,xst(ndim),u(nx1,ny1,nz1,nelv)
      
      real t, Jm
      common /interpl_gen_wk/ t(ldim),Jm(lx1,ldim)
      
      real xlst(3)
      save xlst
      data xlst / 3*1.0e27 /

      real dist
      integer ie, cat
      save    ie

      logical init
      data init /.false./
      save init
      
      if(.not.init) then
         init = .true.
         call bnd_hash_setup
c        call intp_test
      endif

c     Check for quick re-evaluation:

      if (.not.if3d) then
        if( xst(1).eq.xlst(1)  .and.
     $      xst(2).eq.xlst(2) ) then
               call intp_evalpt(u(1,1,1,ie),Jm,ust)
         endif
      else
         if( xst(1).eq.xlst(1)  .and.
     $       xst(2).eq.xlst(2)  .and.
     $       xst(3).eq.xlst(3) )  then
               call intp_evalpt(u(1,1,1,ie),Jm,ust)
             return
         endif
      endif

      if (.not. if3d) then
        call copy(xlst,xst,2)
        xlst(3) = 0.0e0
      else
        call copy(xlst,xst,3)
      endif
      
c     print *, 'interpl_gen called',(xst(k),k=1,3)
      call intp_hash(xst,ie,t,Jm,dist,cat)

      tol = 1.e-13
      if (wdsize.eq.4) tol = 5.e-6

      xnrm = xst(1)**2 + xst(2)**2
      if (if3d) xnrm = xnrm + xst(3)**2 
      if (xnrm.gt.1) tol = tol * sqrt(xnrm)   ! relative tolerance

      if (dist.gt.tol) then
c     if(cat.eq.1 .or. (cat.eq.2 .and. dist.gt.tol)) then
c     if(cat.eq.1 .or. (cat.eq.2 .and. dist.gt.1.e-20)) then
c         print *, 'interpolation failed on:'
c         print *, 'xst=', (xst(id),id=1,ndim)
c         print *, 'cat,dist', cat, dist
c         print *, 'ie,r,s,t=', ie, (t(id),id=1,ndim)
c         call exitt ! sorry, didn't exactly find it
c                    ! although it could be a boundary point (cat=2)
      endif
      call intp_evalpt(u(1,1,1,ie),Jm,ust)
      
      end
      
C=======================================================================
      
      subroutine intp_hash(x,iep,t,Jm,dist,cat)

c----------------------------------------------------------
c     Find (r,s,t,ie) coords for given (x,y,z)
c
c     Input:
c        x       : (x,y,z?) coordinates
c
c     Output:
c        iep  : element point (possibly) found in
c        t    : tentative (r,s,t?) coords for point
c        Jm   : interpolation weights for t (ala intp_wts)
c        dist : distance from x(r,s,t,iep) to input x
c        cat  : 1 - point rejected by bounding box test
c                     for each element
c               2 - best we did was converge to boundary
c                     of one of the elements
c               3 - point found inside an element
c
c----------------------------------------------------------

      include 'SIZE'
      include 'BND'

      real x(ndim),t(ndim),Jm(nx1,ndim),dist
      integer iep,cat
      integer d,i,ii,j,k

      cat=1
      ii=0
c      print *, 'x', (x(id),id=1,ndim)
      do d=ndim,1,-1
         i = bnd_fac(d)*(x(d)-bnd_all(1,d))
c         if(i.lt.0 .or. i.ge.bnd_hn) 
c     $      write(6,*) 'cat',i,d,ii,(x(kk),kk=1,3),bnd_hn
         if(i.lt.0 .or. i.ge.bnd_hn) return
         ii=bnd_hn*ii+i
      enddo
      j=bnd_hash(ii+1)
      k=bnd_hash(ii+2)-j
c      if(k.le.0) write(6,*) 'cat',j,k,ii,(x(kk),kk=1,3)
      if(k.le.0) return
      call intp_findpt(x,bnd_hash(j),k,iep,t,Jm,dist,cat)

      end

C=======================================================================

      subroutine intp_simp(x,iep,t,Jm,dist,cat)

c----------------------------------------------------------
c     Find (r,s,t,ie) coords for given (x,y,z) by searching
c     sequentially through all elements on processor
c
c     Input:
c        x       : (x,y,z?) coordinates
c
c     Output:
c        iep  : element point (possibly) found in
c        t    : tentative (r,s,t?) coords for point
c        Jm   : interpolation weights for t (ala intp_wts)
c        dist : distance from x(r,s,t,iep) to input x
c        cat  : 1 - point rejected by bounding box test
c                     for each element
c               2 - best we did was converge to boundary
c                     of one of the elements
c               3 - point found inside an element
c
c----------------------------------------------------------

      include 'SIZE'
      
      real x(ndim),t(ndim),Jm(nx1,ndim),dist
      integer iep,cat
      
      integer iel
      common /intpsimp/ iel(lelv)
      
      logical ifflag
      data ifflag /.false./
      save ifflag
      
      if(.not.ifflag) then
         ifflag = .true.
         do i=1,nelv
            iel(i)=i
         enddo
      endif
      call intp_findpt(x,iel,nelv,iep,t,Jm,dist,cat)

      end

C=======================================================================

      subroutine intp_findpt(x,iel,nl,iep,t,Jm,dist,cat)

c----------------------------------------------------------
c     Find (r,s,t,ie) coords for given (x,y,z) by searching
c     sequentially through given list of elements
c
c     Input:
c        x       : (x,y,z?) coordinates
c        iel(nl) : list of elements to look through
c
c     Output:
c        iep  : element point (possibly) found in
c        t    : tentative (r,s,t?) coords for point
c        Jm   : interpolation weights for t (ala intp_wts)
c        dist : distance from x(r,s,t,iep) to input x
c        cat  : 1 - point rejected by bounding box test
c                     for each element
c               2 - best we did was converge to boundary
c                     of one of the elements
c               3 - point found inside an element
c
c----------------------------------------------------------

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'BND'

      real x(ndim),t(ndim),dist,Jm(nx1,ndim)
      integer nl, iel(nl), iep, cat

      integer i, ie, nsuc
      real tt(ldim),tJm(lx1,ldim),tdist,tx(ldim) !tentative quantities
      real tr                                    !real temporary


      tol = 1.0e-13
      if (wdsize.eq.4) tol = 5.0e-6

      cat = 1
c     print *, 'ln', nl
c     print *, 'fi',x(1),x(2),x(3)
c     "goto 100" functions as f90's "cycle" here
      do i = 1,nl
         ie = iel(i)
         ! bounding box test
         if(     x(1).lt.bnd_el(1,1,ie)
     $      .or. x(1).gt.bnd_el(2,1,ie)
     $      .or. x(2).lt.bnd_el(1,2,ie)
     $      .or. x(2).gt.bnd_el(2,2,ie)) goto 100
         if(if3d) then
         if(     x(3).lt.bnd_el(1,3,ie)
     $      .or. x(3).gt.bnd_el(2,3,ie)) goto 100
         endif
c         print *, 'in', ie
         ! for added robustness, should do a search to find
         ! initial guess for Newton's method;
         ! for now, we'll just start at center of element
         tt(1)=0
         tt(2)=0
         if(if3d) tt(3)=0
c        print *, 'XI', iel(i), i, (x(k),k=1,ndim)
         call intp_newt(iel(i),x,tt,nsuc)
c        Newton's method really ought to have converged
         if(nsuc.eq.2) then
            if (.not.if3d) then
             write(6,'(2A,1p2E14.6,A,I6)')
     &      'Newton''s method did not converge in intp_findpt ',
     &      'at point ', (x(k),k=1,ndim),' expected to be in elm ',ie
            else
             write(6,'(2A,1p3E14.6,A,I6)')
     &      'Newton''s method did not converge in intp_findpt ',
     &      'at point ', (x(k),k=1,ndim),' expected to be in elm ',ie
            endif
c           call exitt
         endif

         call intp_wts(tt,tJm)
         nxyz = nx1*ny1*nz1
         ioff = nxyz*(ie-1) + 1
         call intp_evalpt(xm1(ioff,1,1,1),tJm,tx(1))
         call intp_evalpt(ym1(ioff,1,1,1),tJm,tx(2))
         if(if3d) call intp_evalpt(zm1(ioff,1,1,1),tJm,tx(3))

c         print *, 'tt', (tt(id),id=1,3)
c         print *, 'tx', (tx(id),id=1,3)

         tr = tx(1)-x(1)
         tdist = tr*tr
         tr = tx(2)-x(2)
         tdist = tdist + tr*tr
         if(if3d) then
            tr = tx(3)-x(3)
            tdist = tdist + tr*tr
         endif

         zdist = sqrt(tdist)
c        write(6,*) 'zdist:',zdist,(x(k),k=1,3)
         if (zdist.le.tol) then ! cat. 3 (pff 12/4/05)
            iep = ie
            call copy(t,tt,ndim)
            call copy(Jm,tJm,nx1*ndim)
            cat = 3
            dist = zdist
            return
         endif

         if(nsuc.ne.0 .and. cat.ne.1
     $      .and. tdist.ge.dist) goto 100
         ! ok, this is the best point so far
         iep = ie
         call copy(t,tt,ndim)
         call copy(Jm,tJm,nx1*ndim)
         dist  = tdist
         if(nsuc.eq.0) then ! in fact, it is category 3
            cat = 3
            dist = sqrt(dist)
            return
         endif
         cat=2  ! just the best border point so far
100      continue
      enddo
      dist = sqrt(dist)

      end

C=======================================================================

      subroutine intp_fullwts(t,Jm,Dm)

c----------------------------------------------------------
c     Interpolation and Derivative Weights Computation
c
c     Input:
c        t(:) : (r,s,t) or (r,s) cooridnates
c     Output:
c        Jm(i,d) := h_i(t(d))
c        Dm(i,d) := h_i'(t(d))
c        where h_i is the ith Lagrangian basis function
c
c----------------------------------------------------------

      include 'SIZE'
      include 'BND' ! for bnd_z, the GLL nodes

      real t(ndim), Jm(nx1,ndim), Dm(nx1,ndim)
      real work
      common /bndwk1/ work(lx1,5)
      integer d

      do d=1,ndim
         call fdwquick01(bnd_z,nx1,t(d),1,
     $                    Jm(1,d),Dm(1,d),bnd_fdw,work)
      enddo
      end

C=======================================================================

      subroutine intp_wts(t,Jm)

c----------------------------------------------------------
c     Interpolation Weights Computation
c
c     Input:
c        t(:) : (r,s,t) or (r,s) cooridnates
c     Output:
c        Jm(i,d) := h_i(t(d))
c        where h_i is the ith Lagrangian basis function
c
c----------------------------------------------------------

      include 'SIZE'
      include 'BND' ! for bnd_z, the GLL nodes

      real t(ndim), Jm(nx1,ndim)
      real work
      common /bndwk1/ work(lx1,3)
      integer d

      do d=1,ndim
         call fdwquick0(bnd_z,nx1,t(d),1,Jm(1,d),bnd_fdw,work)
      enddo
      end

C=======================================================================

      subroutine intp_fulleval(uv,Jm,Dm,u,up)

c----------------------------------------------------------
c     Polynomial + Gradient Evaluation at a Point
c
c     Output:
c        u := u(r,s,t)
c        up(i) := (du/dx_i) (r,s,t)
c     Input:
c        uv : Lagrangian components of polynomial
c        Jm, Dm:  Interpolation, derivative weights
c                 computed from (r,s,t) by intp_wts
c
c----------------------------------------------------------

      include 'SIZE'
      include 'INPUT'

      real uv(nx1,ny1,nz1),u,up(ndim)
      real Jm(nx1,ndim),Dm(nx1,ndim)
      real u1, u1x, u2, u2x, u2y
      common /bndwk1/ u1(ly1,lz1), u1x(ly1,lz1),
     $                u2(lz1), u2x(lz1), u2y(lz1)

      ! Intermediates:
      ! u1[|x] (j,k) := [u|du/dx] (r, z_j, z_k)
      ! u2[|x|y] (k) := [u|du/dx|du/dy] (r, s, z_k)
      call mxm(Jm(1,1),1,uv ,nx1,u1 ,ny1*nz1)
      call mxm(Dm(1,1),1,uv ,nx1,u1x,ny1*nz1)
      call mxm(Jm(1,2),1,u1 ,ny1,u2 ,nz1)
      call mxm(Jm(1,2),1,u1x,ny1,u2x,nz1)
      call mxm(Dm(1,2),1,u1 ,ny1,u2y,nz1)
      if(if3d) then
         call mxm(Jm(1,3),1,u2 ,nz1,u    ,1)
         call mxm(Jm(1,3),1,u2x,nz1,up(1),1)
         call mxm(Jm(1,3),1,u2y,nz1,up(2),1)
         call mxm(Dm(1,3),1,u2 ,nz1,up(3),1)
      else
         u = u2(1)
         up(1) = u2x(1)
         up(2) = u2y(1)
      endif
      end

C=======================================================================

      subroutine intp_evalpt(uv,Jm,u)

c----------------------------------------------------------
c     Polynomial Evaluation at a Point
c
c     Output:
c        u := u(r,s,t)
c     Input:
c        uv : Lagrangian components of polynomial
c        Jm : Interpolation weights
c                 computed from (r,s,t) by intp_wts
c
c----------------------------------------------------------

      include 'SIZE'
      include 'INPUT'

      real uv(nx1,ny1,nz1),u
      real Jm(nx1,ndim)
      real u1, u2
      common /bndwk1/ u1(ly1,lz1), u2(lz1)

      ! Intermediates:
      ! u1(j,k) := u(r, z_j, z_k)
      ! u2(k)   := u(r, s,   z_k)
      call mxm(Jm(1,1),1,uv,nx1,u1,ny1*nz1)
      call mxm(Jm(1,2),1,u1,ny1,u2,nz1)
      if(if3d) then
         call mxm(Jm(1,3),1,u2,nz1,u,1)
      else
         u = u2(1)
      endif
      end

C=======================================================================

      subroutine intp_newt(ie,x,t,success)

c----------------------------------------------------------
c     Newton's Method to find (r,s,t?) corresponding to
c        given (x,y,z?) inside a particular element
c
c     Input:
c        ie : element number
c        x  : (x,y,z?) coordinates
c        t  : initial guess for (r,s,t?) coordinates
c             --- if none available pass (0,0,0?)
c                 as no check is done
c     Output:
c        t  : possibly converged (r,s,t?) coords
c        success :  0 - converged to internal point
c                   1 - converged to boundary point
c                   2 - did not converge
c
c----------------------------------------------------------

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      real x(ndim), t(ndim)
      integer ie, success
      
      real Jm,Dm,dist,ot,jact,cx
      common /bndwk2/ Jm(lx1,ldim),Dm(lx1,ldim),dist
     $               ,ot(ldim),jact(ldim,ldim),cx(ldim)

      integer maxstep
      parameter (maxstep=100)

      integer jacr(ldim), jacc(ldim), d, step
      real tol

      tol = 1.0e-13
      if (wdsize.eq.4) tol = 5.0e-6

      xnrm = x(1)**2 + x(2)**2
      if (if3d) xnrm = xnrm + x(3)**2 
      if (xnrm.gt.1) tol = tol * sqrt(xnrm)   ! relative tolerance

      
      dist = 2.0
      step = 0
      do while (dist.gt.tol .and. step.lt.maxstep)
         step = step + 1
         ! find cx = x(t) and jact = x'(t)^T
         call intp_fullwts(t,Jm,Dm)
         nxyz = nx1*ny1*nz1
         ioff = nxyz*(ie-1) + 1
         call intp_fulleval(
     $         xm1(ioff,1,1,1),Jm,Dm,cx(1),jact(1,1))
         call intp_fulleval(
     $         ym1(ioff,1,1,1),Jm,Dm,cx(2),jact(1,2))
         if(if3d) call intp_fulleval(
     $         zm1(ioff,1,1,1),Jm,Dm,cx(3),jact(1,3))
         ! now x(t+dt) ~= x(t) + x'(t) dt,
         !     x ~= cx + jact^T dt
         !     dt ~= (jact^T)^{-1} (cx - x)
c        print *, 'CI', ie, step, (cx(k),k=1,ndim)
         do d=1,ndim
            cx(d) = x(d)-cx(d)    ! set cx := cx - x
            ot(d) = t(d)
         enddo
         if (if3d) then                             ! cef
           call ludecomp(jact,3,jacr,jacc)
           call lufbsub_t(jact,3,jacr,jacc,cx,t)
         else                                       ! cef
           call ludecomp(jact,2,jacr,jacc)          ! cef
           call lufbsub_t(jact,2,jacr,jacc,cx,t)    ! cef
         endif                                      ! cef
         ! now t = (jact^T)^{-1} cx    (really dt)
         success = 0 ! "found" state
         dist = 0.0
         scal = 1.0
         if (step.ge.10.and.mod(step,5).eq.0) scal = 0.5
         do d=1,ndim
            t(d) = ot(d) + scal*t(d)  ! under-relax; pff; 12/4/05
            ! clip to [-1,1]
            if(t(d).le.-1) then
               t(d) = -1
               success = 1 ! border state
            else if(t(d).ge.1) then
               t(d) = 1
               success = 1 ! border state
            endif
            dist = dist + abs(t(d)-ot(d))
         enddo
c        print *, 'NI', step, dist, (t (d),d=1,ndim)
c        print *, 'XI', step, dist, (cx(d),d=1,ndim)
      enddo
c     print *, 'Newton iterations:', step,tol,wdsize
c     print *, 'XI', step, dist, (x(d),d=1,ndim)
      if(step.ge.maxstep) success=2 ! failure
      end

C=======================================================================

      subroutine ludecomp(A,n,r,c)

c----------------------------------------------------------
c     LU Decomposition   with full pivoting
c
c     Input:
c        A : n by n matrix
c     Computes the LU Decomposition in place
c     Output:
c        r and c are the row and column permutations
c        let Ap be the input matrix A permuted according
c        to r and c:
c           Ap(i,j) = A(r(i),c(j))    !input A
c        let B be the output matrix A permuted according
c        to r and c:
c           B(i,j) = A(r(i),c(j))     !output A (LU decomp)
c        let L and U be given by
c           L(i,j) = 1      if i=j
c                  = 0      if i<j
c                  = B(i,j) if i>j
c           U(i,j) = B(i,j) if i<=j
c                  = 0      if i>j
c        then Ap = L*U
c
c----------------------------------------------------------

      implicit none

      integer n
      real A(n,n)
      integer r(n), c(n)

      integer i,j,k, ii,jj, temp
      real m, mm

      do i=1,n
         r(i)=i
         c(i)=i
      enddo
      do k=1,n-1
         !find maximum element for pivot
         ii = k
         jj = k
         mm = A(r(ii),c(jj))
         do i=k,n
            do j=k,n
               m = A(r(i),c(j))
               if(abs(m).gt.abs(mm)) then
                  mm = m
                  ii = i
                  jj = j
               endif
            enddo
         enddo
         !swap maximum element into place
         temp = r(ii)
         r(ii) = r(k)
         r(k) = temp
         temp = c(jj)
         c(jj) = c(k)
         c(k) = temp
         !now mm = A(r(k),c(k)) is our pivot
         if(mm.eq.0) return  ! singular matrix
         mm = 1/mm  ! factor out the division
         !proceed with LU decomposition
         do i=k+1,n
            A(r(i),c(k)) = A(r(i),c(k)) * mm
         enddo
         do j=k+1,n
            do i=k+1,n
               A(r(i),c(j)) = A(r(i),c(j)) 
     $                      - A(r(i),c(k))*A(r(k),c(j))
            enddo
         enddo
      enddo
      end

C=======================================================================

      subroutine lufbsub(A,n,r,c,b,x)

c----------------------------------------------------------
c     LU Decomposition   forward & backward substitution
c
c     Once ludecomp(A,n,r,c) has been called,
c     lufbsub(A,n,r,c,b,x) may be called multiple times
c     to yield solutions to
c        A x = b
c     where A is the matrix passed to ludecomp
c
c     b is overwritten
c
c----------------------------------------------------------

      implicit none

      integer n
      real A(n,n), b(n), x(n)
      integer r(n), c(n)
      
      integer i,j,k

      do k=1,n-1
         do i=k+1,n
            b(r(i)) = b(r(i))-A(r(i),c(k))*b(r(k))
         enddo
      enddo
      do j=n,1,-1
         x(c(j)) = b(r(j))/A(r(j),c(j))
         do i=1,j-1
            b(r(i)) = b(r(i))-A(r(i),c(j))*x(c(j))
         enddo
      enddo
      end

C=======================================================================

      subroutine lufbsub_t(A,n,r,c,b,x)

c----------------------------------------------------------
c     LU Decomposition   forward & backward substitution
c                        for the transpose system
c
c     Once ludecomp(A,n,r,c) has been called,
c     lufbsub_t(A,n,r,c,b,x) may be called multiple times
c     to yield solutions to
c        A' x = b        (note the transpose)
c     where A is the matrix passed to ludecomp
c
c     b is overwritten
c
c----------------------------------------------------------

      implicit none

      integer n
      real A(n,n), b(n), x(n)
      integer r(n), c(n)
      
      integer i,j,k

      do k=1,n
         b(c(k)) = b(c(k)) / A(r(k),c(k))
         do i=k+1,n
            b(c(i)) = b(c(i))-A(r(k),c(i))*b(c(k))
         enddo
      enddo
      do j=n,1,-1
         x(r(j)) = b(c(j))
         do i=1,j-1
            b(c(i)) = b(c(i))-A(r(j),c(i))*x(r(j))
         enddo
      enddo
      end

C=======================================================================

      subroutine fdweights01(z,n,x,m,Jm,Dm,work)

c----------------------------------------------------------
c     Finite Difference Weights for 0th and 1st derivative
c
c     Given:  n nodes  z ,...,z    and m nodes  x ,...,x
c                       1      n                 1      m
c
c     Computes:  weights Jm(m,n) and Dm(m,n) such that
c              n                            n
c             __                           __
c     f(x ) = >   Jm   f(z )  and f'(x ) = >   Dm   f(z )
c        i    --    ij    j           i    --    ij    j
c             j=1                          j=1
c
c     whenever f is a polynomial of degree at most n-1
c
c     User must provide work array with space for 6*n reals
c----------------------------------------------------------

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), Dm(m,n), work(n,6)

      call fdweights01_hlpr(z,n,x,m,Jm,Dm,
     $       work(1,1), work(1,2), work(1,3),
     $       work(1,4), work(1,5), work(1,6))      
      
      end

C=======================================================================
      subroutine fdweights01_hlpr(z,n,x,m,Jm,Dm,w,d,u,v,up,vp)

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), Dm(m,n)
      real w(n), d(n), u(n), v(n), up(n), vp(n)
      
      integer i,j
      
      call fdwquick_init(z,n,w)
      u(1) = 1.0
      v(n) = 1.0
      up(1) = 0.0
      vp(n) = 0.0
      do i=1,m
         do j=1,n
            d(j) = x(i)-z(j)
         enddo
         do j=2,n
            u(j) = d(j-1)*u(j-1)
            up(j) = d(j-1)*up(j-1) + u(j-1)
         enddo
         do j=n-1,1,-1
            v(j) = d(j+1)*v(j+1)
            vp(j) = d(j+1)*vp(j+1) + v(j+1)
         enddo
         do j=1,n
            Jm(i,j) = w(j)*u(j)*v(j)
            Dm(i,j) = w(j)*(up(j)*v(j)+u(j)*vp(j))
         enddo
      enddo
      
      end

C=======================================================================

      subroutine fdweights0(z,n,x,m,Jm,work)

c----------------------------------------------------------
c     Finite Difference Weights for 0th derivative
c
c     Given:  n nodes  z ,...,z    and m nodes  x ,...,x
c                       1      n                 1      m
c
c     Computes:  weights Jm(m,n)
c              n
c             __
c     f(x ) = >   Jm   f(z )
c        i    --    ij    j
c             j=1
c
c     whenever f is a polynomial of degree at most n-1
c
c     User must provide work array with space for 4*n reals
c----------------------------------------------------------

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), work(n,4)

      call fdweights0_hlpr(z,n,x,m,Jm,
     $       work(1,1), work(1,2),
     $       work(1,3), work(1,4))
      
      end

C=======================================================================

      subroutine fdweights0_hlpr(z,n,x,m,Jm,w,d,u,v)

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n)
      real w(n), d(n), u(n), v(n)
      
      integer i,j

      call fdwquick_init(z,n,w)
      u(1) = 1.0
      v(n) = 1.0
      do i=1,m
         do j=1,n
            d(j) = x(i)-z(j)
         enddo
         do j=2,n
            u(j) = d(j-1)*u(j-1)
         enddo
         do j=n-1,1,-1
            v(j) = d(j+1)*v(j+1)
         enddo
         do j=1,n
            Jm(i,j) = w(j)*u(j)*v(j)
         enddo
      enddo
      
      end

C=======================================================================

      subroutine fdwquick_init(z,n,w)

C-----------------------------------------------------
C  Compute the quantities 
C               
C           N
C         ------
C    w_i = |  |  (Z_i - Z_j))
C          |  |
C       j=1, j \ne i
C
C-----------------------------------------------------

      implicit none

      integer n
      real z(n), w(n)
      integer i,j
      
      do i=1,n
         w(i) = 1.0
         do j=1,n
            if(i.ne.j) w(i)=w(i)*(z(i)-z(j))
         enddo
         w(i) = 1.0 / w(i)
      enddo

      end

C=======================================================================

      subroutine fdwquick01(z,n,x,m,Jm,Dm,w,work)

c----------------------------------------------------------
c     fdwquick01 is identical to fdweights01,
c         except that it is assumed that
c         fdwquick_init(z,n,w) has been called,
c         and work need only hold 5*n reals (not 6*n)
c----------------------------------------------------------

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), Dm(m,n), w(n), work(n,5)

      call fdwquick01_hlpr(z,n,x,m,Jm,Dm,w,
     $                  work(1,1), work(1,2),
     $       work(1,3), work(1,4), work(1,5))
      
      end

C=======================================================================

      subroutine fdwquick01_hlpr(z,n,x,m,Jm,Dm,w,d,u,v,up,vp)

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), Dm(m,n)
      real w(n), d(n), u(n), v(n), up(n), vp(n)
      
      integer i,j
      
      u(1) = 1.0
      v(n) = 1.0
      up(1) = 0.0
      vp(n) = 0.0
      do i=1,m
         do j=1,n
            d(j) = x(i)-z(j)
         enddo
         do j=2,n
            u(j) = d(j-1)*u(j-1)
            up(j) = d(j-1)*up(j-1) + u(j-1)
         enddo
         do j=n-1,1,-1
            v(j) = d(j+1)*v(j+1)
            vp(j) = d(j+1)*vp(j+1) + v(j+1)
         enddo
         do j=1,n
            Jm(i,j) = w(j)*u(j)*v(j)
            Dm(i,j) = w(j)*(up(j)*v(j)+u(j)*vp(j))
         enddo
      enddo
      
      end

C=======================================================================

      subroutine fdwquick0(z,n,x,m,Jm,w,work)

c----------------------------------------------------------
c     fdwquick0 is identical to fdweights0,
c         except that it is assumed that
c         fdwquick_init(z,n,w) has been called,
c         and work need only hold 3*n reals (not 4*n)
c----------------------------------------------------------

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), w(n), work(n,3)

      call fdwquick0_hlpr(z,n,x,m,Jm,w,
     $       work(1,1), work(1,2), work(1,3))
      
      end

C=======================================================================

      subroutine fdwquick0_hlpr(z,n,x,m,Jm,w,d,u,v)

      implicit none

      integer n,m
      real z(n), x(m), Jm(m,n), Dm(m,n)
      real w(n), d(n), u(n), v(n), up(n), vp(n)
      
      integer i,j
      
      u(1) = 1.0
      v(n) = 1.0
      do i=1,m
         do j=1,n
            d(j) = x(i)-z(j)
         enddo
         do j=2,n
            u(j) = d(j-1)*u(j-1)
         enddo
         do j=n-1,1,-1
            v(j) = d(j+1)*v(j+1)
         enddo
         do j=1,n
            Jm(i,j) = w(j)*u(j)*v(j)
         enddo
      enddo
      
      end
