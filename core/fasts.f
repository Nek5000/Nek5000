c-----------------------------------------------------------------------
      subroutine local_solves_fdm(u,v)
c
c     Given an input vector v, this returns the additive Schwarz solution
c
c                       -1
c                  u = M   v
c
c                      T     -1  
c     where M = sum ( R_i A_i    R_i )
c                i
c
c     The R_i's are simply index_set restriction operators.
c
c     The local solves are performed with the fast diagonalization method.
c
      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'PARALLEL'
c
      include 'CTIMER'
c
      real u(lx2,ly2,lz2,lelv),v(lx2,ly2,lz2,lelv)
      common /scrpre/ v1(lx1,ly1,lz1,lelv)
     $               ,w1(lx1,ly1,lz1),w2(lx1,ly1,lz1)
      common /scrover/ ar(lelv)
c
      if (icalld.eq.0) tsolv=0.0
      icalld=icalld+1
      nsolv=icalld
c
      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv
c
c     Fill interiors
      iz1 = 0
      if (if3d) iz1=1
      call rzero(v1,ntot1)
      do ie=1,nelv
         do iz=1,nz2
         do iy=1,ny2
         do ix=1,nx2
            v1(ix+1,iy+1,iz+iz1,ie) = v(ix,iy,iz,ie)
         enddo
         enddo
         enddo
      enddo
      call dface_ext    (v1)
      call dssum        (v1,nx1,ny1,nz1)
      call dface_add1si (v1,-1.)
c
c     Now solve each subdomain problem:
c
      etime1=dclock()
      do ie = 1,nelv
         call fastdm1(v1(1,1,1,ie),ie,w1,w2)
      enddo
      tsolv=tsolv+dclock()-etime1
c
c     Exchange/add elemental solutions
c
      call s_face_to_int (v1,-1.)
      call dssum         (v1,nx1,ny1,nz1)
      call s_face_to_int (v1, 1.)
      if(param(42).eq.0) call do_weight_op(v1)
c
c     Map back to pressure grid (extract interior values)
c
      do ie=1,nelv
         do iz=1,nz2
         do iy=1,ny2
         do ix=1,nx2
            u(ix,iy,iz,ie) = v1(ix+1,iy+1,iz+iz1,ie)
         enddo
         enddo
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fastdm1(R,ie,w1,w2)
c
c     Fast diagonalization solver for FEM on mesh 1
c
      include 'SIZE'
      parameter (lxx=lx1*lx1)
      common /fastd/  sr(lxx,2,lelv),ss(lxx,2,lelv),st(lxx,2,lelv)
     $             ,  df(lx1*ly1*lz1,lelv)
c
      parameter (lxyz = lx1*ly1*lz1)
c
      real r(1),w1(1),w2(1)
c
c
c      T
c     S  r
      call tensr3 (w1,nx1,r ,nx1,sr(1,2,ie),ss(1,1,ie),st(1,1,ie),w2)
c
c     
c      -1 T
c     D  S  r
c
      call col2   (w1,df(1,ie),lxyz)
c
c
c        -1 T
c     S D  S  r
c
      call tensr3 (r ,nx1,w1,nx1,sr(1,1,ie),ss(1,2,ie),st(1,2,ie),w2)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine tensr3(v,nv,u,nu,A,Bt,Ct,w)
C
C     -  Tensor product application of v = (C x B x A) u
C        NOTE -- the transpose of B & C must be input, rather than B & C.
C
C     -  scratch arrays: w(nu*nu*nv)
C
C
      include 'SIZE'
      include 'INPUT'
      real v(nv,nv,nv),u(nu,nu,nu)
      real A(1),B(1),C(1)
      real w(1)

      if (nu.gt.nv) then
         write(6,*) nid,nu,nv,' ERROR in tensr3. Contact P.Fischer.'
         write(6,*) nid,nu,nv,' Memory problem.'
         call exitt
      endif

      if (if3d) then
         nuv = nu*nv
         nvv = nv*nv
         call mxm(A,nv,u,nu,v,nu*nu)
         k=1
         l=1
         do iz=1,nu
            call mxm(v(k,1,1),nv,Bt,nu,w(l),nv)
            k=k+nuv
            l=l+nvv
         enddo
         call mxm(w,nvv,Ct,nu,v,nv)
      else
         call mxm(A,nv,u,nu,w,nu)
         call mxm(w,nv,Bt,nu,v,nv)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine s_face_to_int(x,c)
c
c     Scale face and add to interior of element     
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelv
c
        if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,2    ,iz,ie) = c*x(ix,1  ,iz,ie) + x(ix,2    ,iz,ie)
            x(ix,ny1-1,iz,ie) = c*x(ix,ny1,iz,ie) + x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(2    ,iy,iz,ie) = c*x(1  ,iy,iz,ie) + x(2    ,iy,iz,ie)
            x(nx1-1,iy,iz,ie) = c*x(nx1,iy,iz,ie) + x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,2    ,ie) = c*x(ix,iy,1  ,ie) + x(ix,iy,2    ,ie)
            x(ix,iy,nz1-1,ie) = c*x(ix,iy,nz1,ie) + x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c         2D
          do ix=2,nx1-1
            x(ix,2    ,1,ie) = c*x(ix,1  ,1,ie) + x(ix,2    ,1,ie)
            x(ix,ny1-1,1,ie) = c*x(ix,ny1,1,ie) + x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(2    ,iy,1,ie) = c*x(1  ,iy,1,ie) + x(2    ,iy,1,ie)
            x(nx1-1,iy,1,ie) = c*x(nx1,iy,1,ie) + x(nx1-1,iy,1,ie)
          enddo
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_ext(x)
c     Extend interior to face of element     
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelv
c
         if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,1  ,iz,ie) = x(ix,2    ,iz,ie)
            x(ix,ny1,iz,ie) = x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(1  ,iy,iz,ie) = x(2    ,iy,iz,ie)
            x(nx1,iy,iz,ie) = x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,1  ,ie) = x(ix,iy,2    ,ie)
            x(ix,iy,nz1,ie) = x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c
          do ix=2,nx1-1
            x(ix,1  ,1,ie) = x(ix,2    ,1,ie)
            x(ix,ny1,1,ie) = x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(1  ,iy,1,ie) = x(2    ,iy,1,ie)
            x(nx1,iy,1,ie) = x(nx1-1,iy,1,ie)
          enddo
c
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine dface_add1si(x,c)
c     Scale interior and add to face of element     
c
      include 'SIZE'
      include 'INPUT'
      real x(nx1,ny1,nz1,1)
c
      do ie=1,nelv
c
        if (if3d) then
c
          do iz=2,nz1-1
          do ix=2,nx1-1
            x(ix,1  ,iz,ie) = x(ix,1  ,iz,ie) + c*x(ix,2    ,iz,ie)
            x(ix,ny1,iz,ie) = x(ix,ny1,iz,ie) + c*x(ix,ny1-1,iz,ie)
          enddo
          enddo
c
          do iz=2,nz1-1
          do iy=2,ny1-1
            x(1  ,iy,iz,ie) = x(1  ,iy,iz,ie) + c*x(2    ,iy,iz,ie)
            x(nx1,iy,iz,ie) = x(nx1,iy,iz,ie) + c*x(nx1-1,iy,iz,ie)
          enddo
          enddo
c
          do iy=2,ny1-1
          do ix=2,nx1-1
            x(ix,iy,1  ,ie) = x(ix,iy,1  ,ie) + c*x(ix,iy,2    ,ie)
            x(ix,iy,nz1,ie) = x(ix,iy,nz1,ie) + c*x(ix,iy,nz1-1,ie)
          enddo
          enddo
c
        else
c
c         2D
c
          do ix=2,nx1-1
            x(ix,1  ,1,ie) = x(ix,1  ,1,ie) + c*x(ix,2    ,1,ie)
            x(ix,ny1,1,ie) = x(ix,ny1,1,ie) + c*x(ix,ny1-1,1,ie)
          enddo
          do iy=2,ny1-1
            x(1  ,iy,1,ie) = x(1  ,iy,1,ie) + c*x(2    ,iy,1,ie)
            x(nx1,iy,1,ie) = x(nx1,iy,1,ie) + c*x(nx1-1,iy,1,ie)
          enddo
c
        endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine init_weight_op

      include 'SIZE'
      include 'INPUT'
      common /swaplengths/ l(lx1,ly1,lz1,lelv)
      common /weightop/ w(lx2,lz2,2,3,lelv)
      real l
      real w

      n=nx2+1
      if (if3d) then
         do ie=1,nelv
            call rzero(l(1,1,1,ie),nx1*ny1*nz1)
            do k=2,n
               do j=2,n
                  l(2,j,k,ie)=1
                  l(n,j,k,ie)=1
               enddo
            enddo
            do k=2,n
               do i=2,n
                  l(i,2,k,ie)=1
                  l(i,n,k,ie)=1
               enddo
            enddo
            do j=2,n
               do i=2,n
                  l(i,j,2,ie)=1
                  l(i,j,n,ie)=1
               enddo
            enddo
         enddo
      else
         do ie=1,nelv
            call rzero(l(1,1,1,ie),nx1*ny1*nz1)
            do j=2,n
               l(2,j,1,ie)=1
               l(n,j,1,ie)=1
            enddo
            do i=2,n
               l(i,2,1,ie)=1
               l(i,n,1,ie)=1
            enddo
         enddo
      endif
      
      call dface_ext(l)
      call dssum(l,nx1,ny1,nz1)
      call dface_add1si(l,-1.)
      call s_face_to_int(l,1.)
c     l now holds the count matrix C on the outer pressure nodes
      if (if3d) then
         do ie=1,nelv
            do k=1,nz2
               do j=1,ny2
                  w(j,k,1,1,ie)=1.0/l(2,j+1,k+1,ie)
                  w(j,k,2,1,ie)=1.0/l(n,j+1,k+1,ie)
               enddo
            enddo
            do k=1,nz2
               do i=1,nx2
                  w(i,k,1,2,ie)=1.0/l(i+1,2,k+1,ie)
                  w(i,k,2,2,ie)=1.0/l(i+1,n,k+1,ie)
               enddo
            enddo
            do j=1,ny2
               do i=1,nx2
                  w(i,j,1,3,ie)=1.0/l(i+1,j+1,2,ie)
                  w(i,j,2,3,ie)=1.0/l(i+1,j+1,n,ie)
               enddo
            enddo
         enddo
      else
         do ie=1,nelv
            do j=1,ny2
               w(j,1,1,1,ie)=1.0/l(2,j+1,1,ie)
               w(j,1,2,1,ie)=1.0/l(n,j+1,1,ie)
            enddo
            do i=1,nx2
               w(i,1,1,2,ie)=1.0/l(i+1,2,1,ie)
               w(i,1,2,2,ie)=1.0/l(i+1,n,1,ie)
            enddo
         enddo
      endif
      end
c-----------------------------------------------------------------------
      subroutine do_weight_op(x)
      include 'SIZE'
      include 'INPUT'
      common /weightop/ w(lx2,lz2,2,3,lelv)
      real w

      real x(0:nx1-1,0:ny1-1,0:nz1-1,1)

      if (if3d) then
         do ie=1,nelv
            do k=1,nz2
               do j=1,ny2
                  x(  1,j,k,ie)=w(j,k,1,1,ie)*x(  1,j,k,ie)
                  x(nx2,j,k,ie)=w(j,k,2,1,ie)*x(nx2,j,k,ie)
               enddo
            enddo
            do k=1,nz2
               do i=2,nx2-1
                  x(i,  1,k,ie)=w(i,k,1,2,ie)*x(i,  1,k,ie)
                  x(i,ny2,k,ie)=w(i,k,2,2,ie)*x(i,ny2,k,ie)
               enddo
            enddo
            do j=2,ny2-1
               do i=2,nx2-1
                  x(i,j,  1,ie)=w(i,j,1,3,ie)*x(i,j,  1,ie)
                  x(i,j,nz2,ie)=w(i,j,2,3,ie)*x(i,j,nz2,ie)
               enddo
            enddo
         enddo
      else
         do ie=1,nelv
            do j=1,ny2
               x(  1,j,0,ie)=w(j,1,1,1,ie)*x(  1,j,0,ie)
               x(nx2,j,0,ie)=w(j,1,2,1,ie)*x(nx2,j,0,ie)
            enddo
            do i=2,nx2-1
               x(i,  1,0,ie)=w(i,1,1,2,ie)*x(i,  1,0,ie)
               x(i,ny2,0,ie)=w(i,1,2,2,ie)*x(i,ny2,0,ie)
            enddo
         enddo
      endif
      end
