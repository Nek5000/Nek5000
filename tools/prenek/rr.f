c-----------------------------------------------------------------------
      program amat_test
c
c     Create a stiffness matrix for the model FEM 
c     problem discussed on 4/29/97
c
      parameter (ntmax=200)
      common /c1/    x(3,ntmax),y(3,ntmax)
      common /c2/    f(3,ntmax)
c
      common /c3/    bf(3,ntmax),u(3,ntmax),au(3,ntmax)
      common /c4/    vx(3,ntmax),vy(3,ntmax),w(3,ntmax)
c
      parameter (nabd = ntmax*40)
      common /c5/    abd(nabd),acrs(ntmax*10)
c
      common /c6/    b(3,ntmax),e2v(3,ntmax)
      integer        b         ,e2v
      common /c7/    ia(ntmax),ja(10*ntmax),order(3*ntmax/2)
      integer        order
      common /c8/    il(ntmax),jl(nabd)
c
      real*4 etime,eee(2)
      eee(1) =0.
      eee(2) =0.
c
      npass = 5
      do ipass = 1,npass
         write(6,*) 'this is ipass',ipass,npass
c
         if (ipass.eq.1) then
            call init_tri   (x,y,b,nt,f,abd)
         else
            call ref_triangles(x,y,b,nt)
         endif
         if (nt.gt.ntmax) then
            write(6,*) 'increase ntmax from:',ntmax,' to',nt
            call exit
         endif
c
c        Copy x,y, and b to work arrays so they're not destroyed
c
         call  copy(abd(nt*12 + 1),x,3*nt)
         call  copy(abd(nt*15 + 1),y,3*nt)
         call icopy(abd(nt*18 + 1),b,3*nt)
         call set_e2v (e2v,n,m
     $                   ,abd(nt*12+1),abd(nt*15+1),abd(nt*18+1),nt
     $                   ,f,abd(1),abd(nt*3+1),abd(nt*6+1),abd(nt*9+1))
         call find_lda (lda,n,m,e2v,b,nt,.false.)
c
c
         if (ipass.lt.0) then
c           Postscript plot
            io = 20+ipass
            call plot_mesh (x,y,3,nt,io)
            call plot_vert (x,y,e2v,3,nt,io)
            call close_plot(io)
         endif
c
c
c        Set up nested dissection ordering for sparse solver
c
c
         ndim = 2
         nve  = ndim+1
c
         call izero  (order,n)
         call rsb_fem(order,il,ia,ja,e2v,nve,nt,abd,nabd,ndim)
         do i=1,n
            order(i) = n+1 - order(i)
         enddo
c
c        The output of the rsb ordering is given in "order()".
c        It is the array of global index numbers to which each
c        e2v() value should be redefined to have near optimal
c        fill in the Cholesky decomposition...
c
c        Renumber e2v
c
         do it=1,nt
         do iv=1,nve
            e2v(iv,it) = order(e2v(iv,it))
         enddo
         enddo
c
c        Renumber bc's to be last
c
         call bc_relabel(e2v,b,nve,nt,n,abd)
         call find_lda  (lda,n,m,e2v,b,nt,.true.)
c
         if (ipass.lt.0) then
            io = 80+ipass
            call plot_mesh (x,y,3,nt,io)
            call plot_vert (x,y,e2v,3,nt,io)
            call close_plot(io)
         endif
c
c
         ldan = lda*n
         if (ldan.gt.nabd) then
            write(6,*) 'ABORT.  Increase nabd from',nabd,' to >',ldan
            call exit
         endif
c
c
c        Solve this problem using sparse matrix solver
c
         m = 3*nt+20
         call gen_graph_l(ia,ja,nr,e2v,b,nt
     $                 ,abd,abd(m),abd(2*m),abd(3*m),m)
c        call outjmat(ja,ia,n,'Acrs')
         call gen_acrs (acrs,ia,ja,nr,x,y,e2v,b,nt)
c
c
         call gen_abd(abd,lda,n,x,y,b,e2v,nt)
c
         call set_ubdry (u,x,y,b,nt)
         call set_Abdry (au,u,n,x,y,b,e2v,nt)
         call set_bf    (bf,  n,x,y,b,e2v,nt)
         call sub2      (bf,au,n)
c
         call copy      (au,bf,n)
         call copy      (f ,u ,nt*3)
c
         t1 = etime(eee)
c        call out_mat   (abd,lda,n,'Fabd',6)
         call solve_fem (u,bf,n,abd,lda,b,e2v,nt)
c        call out_mat   (abd,lda,n,'Fabd',6)
         t2 = etime(eee)
         db = t2-t1
c
         ldl = nabd
         t1 = etime(eee)
         call solve_spr (f,au,abd,il,jl,ldl,acrs,ia,ja,nr,b,e2v,nt,w)
         t2 = etime(eee)
         ds = t2-t1
         call sub2      (au,bf,n)
         dmax = glamax(au,n)
         write(6,8) 'max dif:',n,dmax,ds,db
   8     format(a8,i9,1p3e12.4)
c
c
c        ioc = 40+ipass
c        call do_contour(u,x,y,nt,ioc)
c
c        call velocity(vx,vy,u,x,y,e2v,nt,abd,f,m)
c
c        iou = 50+ipass
c        iov = 60+ipass
c
c        call do_contour(vx,x,y,nt,iou)
c        call do_contour(vy,x,y,nt,iov)
c
      enddo
c
      stop
      end
c-----------------------------------------------------------------------
      subroutine gen_abd(abd,lda,n,x,y,b,e2v,nt)
c
c     Generate banded stiffness matrix  Abd.
c
c     Note:  the "1's" in the next 2 lines can also be n and nt, resp.
      real abd(lda,1),x(3,1),y(3,1)
      integer b(3,1),e2v(3,1)
c
c     Initialize abd:
      call rzero (abd, lda*n)
c
      do k=1,nt
         Ak4 = 4.0*tri_area(x(1,k),y(1,k))
         do j=1,3
            if (b(j,k).le.0) then
c
c              in interior for j...
               jj = e2v(j,k)
               j1 = mod(j  ,3) + 1
               j2 = mod(j+1,3) + 1
               dxj = x(j1,k) - x(j2,k)
               dyj = y(j1,k) - y(j2,k)
c
               do i=1,3
                  if (b(i,k).le.0) then
c                    in interior for i...
c
                     ii = e2v(i,k)
                     if (jj.ge.ii) then
c                       only compute upper-half of A, due to symmetry...
c
                        i1 = mod(i  ,3) + 1
                        i2 = mod(i+1,3) + 1
                        dxi = x(i1,k) - x(i2,k)
                        dyi = y(i1,k) - y(i2,k)
c
                        call bnd_index(ib,jb,ii,jj,lda)
                        dabd =  (dxi*dxj + dyi*dyj)/Ak4
                        abd(ib,jb)=abd(ib,jb) + dabd
                     endif
                  endif
               enddo
            endif
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine out_mat(a,m,n,name,io)
      character*4 name
      real a(m,n)
c
c
      n8 = min(10,n)
c
      write(io,1) name,m,n
    1 format('Matrix: ',a4,'  m =',i6,'  n =',i6)
c
      m8 = min(m,30)
      m8 = m
      do i=1,m8
         write(io,2) (a(i,j),j=1,n8)
      enddo
    2 format(10f9.6)
c   2 format(1p10e12.5)
      return
      end
c-----------------------------------------------------------------------
      subroutine bnd_index(ib,jb,ii,jj,lda)
c
      jb = jj
      ib = lda + ii - jj
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outv(x,n,name)
      real x(1)
      character*3 name
c
      n10 = min(n,10)
      write(6,1) name,(x(i),i=1,n10)
    1 format(a3,1p10e12.4)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_e2v (e2v,n,m,x,y,b,nt
     $                      ,ind,init_addr,ifblock,n_in_block,work)
      integer e2v(1),b(1),init_addr(1),n_in_block(1),ind(1)
      real    x(1),y(1),work(1)
      logical ifblock(1)
c
      nv = 3*nt
      do i=1,nv
         init_addr(i) = i
         ifblock  (i) = .false.
      enddo
c
      call set_tol(tol,x,y,nt)
c
c     Preprocess geometry via slight rotation
c
      call vrotate(x,y,nv,-.10 )
c
      ifblock(1)    = .true.
      n_in_block(1) = nv
      nblocks       = 1
      do idim = 0,2
         ib = 1
         do iblock = 1,nblocks
            nb = n_in_block(iblock)
            if (idim.eq.0) then
c              sort by bc's
               call irank (b(ib),ind,nb)
            elseif (idim.eq.1) then
               call  rank (x(ib),ind,nb)
            elseif (idim.eq.2) then
               call  rank (y(ib),ind,nb)
            endif
c
            call iswap (b(ib)         ,ind,nb,work)
            call  swap (x(ib)         ,ind,nb,work)
            call  swap (y(ib)         ,ind,nb,work)
            call iswap (init_addr(ib) ,ind,nb,work)
c
            ib = ib + n_in_block(iblock)
         enddo
c
c        Data now sorted within each block.  Now scan for start
c        of new blocks
c
         if (idim.eq.0) then
            do i=2,nv
               if (b(i).ne.b(i-1)) ifblock(i)=.true.
            enddo
         elseif (idim.eq.1) then
            do i=2,nv
               if (abs(x(i)-x(i-1)).gt.tol) ifblock(i)=.true.
            enddo
         elseif (idim.eq.2) then
            do i=2,nv
               if (abs(y(i)-y(i-1)).gt.tol) ifblock(i)=.true.
            enddo
         endif
c
c        Now, count number of blocks and number in each block
c
         nblocks = 1
         ilast   = 1
         do i=2,nv
            if (ifblock(i)) then
               n_in_block(nblocks) = i-ilast
               nblocks             = nblocks+1
               ilast               = i
            endif
         enddo
         n_in_block(nblocks) = nv+1 - ilast
c
         if (nv.le.200) write(6,6) (ifblock(i),i=1,nv)
    6    format(80l1)
c
      enddo
c
c     Data all sorted... now count number of blocks, number in each block.
c
      nblocks = 1
      ilast   = 1
      e2v(init_addr(1)) = nblocks
      do i=2,nv
         if (ifblock(i)) then
            n_in_block(nblocks) = i-ilast
            nblocks             = nblocks+1
            ilast               = i
         endif
c
c        Assign e2v
c
         e2v(init_addr(i)) = nblocks
c
      enddo
      m = nblocks
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_tol(tol,x,y,nt)
      real x(3,1),y(3,1)
c
      xmin = x(1,1)
      ymin = y(1,1)
      xmax = x(1,1)
      ymax = y(1,1)
      do i=1,3*nt
         xmin = min(xmin,x(i,1))
         ymin = min(ymin,y(i,1))
         xmax = max(xmax,x(i,1))
         ymax = max(ymax,y(i,1))
      enddo
      d2m = (xmax-xmin)**2 + (ymax-ymin)**2
c
      do it=1,nt
         do iv=1,3
            iv1 = mod(iv,3) + 1
            d2  = (x(iv1,it)-x(iv,it))**2
     $          + (y(iv1,it)-y(iv,it))**2
            d2m = min(d2,d2m)
         enddo
      enddo
      tol = .001*d2m
      return
      end
c-----------------------------------------------------------------------
      subroutine bc_relabel(e2v,b,nv,nt,n,iw)
      integer e2v(nv,1),b(nv,1),iw(0:1)
c
c     Relabel boundary nodes last  (boundary nodes are denoted
c     by a positive integer in b(iv,ie) )
c
      do it=1,nt
         do iv=1,3
            if (b(iv,it).gt.0) e2v(iv,it) = e2v(iv,it) + n
         enddo
      enddo
c
c     Compress this list of integers to eliminate gaps
c
      call compress_int(e2v,nv*nt,iw,iw(nv*nt))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine find_lda (lda,n,m,e2v,b,nt,ifbdry)
      integer e2v(3,1),b(3,1)
      logical ifbdry
c
c     Compute lda according to max |I-J| over all triangles
c
      lda = 0
      nn  = 0
      do it=1,nt
         do iv=1,3
            iv1 = mod(iv,3) + 1
            if (ifbdry) then
               if (b(iv1,it).le.0 .and. b(iv,it).le.0) then
                  ldt = abs(e2v(iv1,it) - e2v(iv,it))
                  lda = max(lda,ldt)
                  lda = max(lda,ldt)
                  nn  = max(nn ,e2v(iv,it))
               endif
            else
               ldt = abs(e2v(iv1,it) - e2v(iv,it))
               lda = max(lda,ldt)
               lda = max(lda,ldt)
               nn  = max(nn ,e2v(iv,it))
            endif
         enddo
      enddo
      lda = lda + 1
      n   = nn
c
      write(6,*) 'this is lda:',lda,n,nt,ifbdry
c
      return
      end
c-----------------------------------------------------------------------
      function func(x,y)
      func = 0.0
      return
      end
c-----------------------------------------------------------------------
      subroutine set_bf(bf,n,x,y,b,e2v,nt)
      real bf(1)
      real x(3,1),y(3,1)
      integer b(3,1),e2v(3,1)
c
      real    bl(3,3)
      save    bl
      data    bl / 2.,1.,1.,1.,2.,1.,1.,1.,2. /
c
c     First call, divide each entry by 12
c
      if (bl(1,1).eq.2.) then
         do i=1,9
            bl(i,1) = bl(i,1)/12.
         enddo
      endif
c
c     Compute bf = B*f
c
      call rzero(bf,n)
c
      do k=1,nt
         ak = tri_area(x(1,k),y(1,k))
         do j=1,3
            jj = e2v(j,k)
            do i=1,3
               if (b(i,k).le.0) then
                  ii = e2v(i,k)
                  bf(ii) = bf(ii) + ak*bl(i,j)*func(x(j,k),y(j,k))
               endif
            enddo
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_Abdry (au,u,n,x,y,b,e2v,nt)
c
c     Generate A*u on the boundary
c
      real    au(1)
      real    u(3,1),x(3,1),y(3,1)
      integer b(3,1),e2v(3,1)
c
      call rzero(au,n)
c
      do k=1,nt
         Ak4 = 4.0*tri_area(x(1,k),y(1,k))
         do j=1,3
            if (b(j,k).gt.0) then
c              On boundary for j
c
               jj = e2v(j,k)
               j1 = mod(j  ,3) + 1
               j2 = mod(j+1,3) + 1
               dxj = x(j1,k) - x(j2,k)
               dyj = y(j1,k) - y(j2,k)
c
               do i=1,3
                  if (b(i,k).le.0) then
c                    In interior for i...
c
                     ii = e2v(i,k)
c
                     i1 = mod(i  ,3) + 1
                     i2 = mod(i+1,3) + 1
                     dxi = x(i1,k) - x(i2,k)
                     dyi = y(i1,k) - y(i2,k)
c
                     da     = (dxi*dxj + dyi*dyj)/Ak4
                     au(ii) = au(ii) + da*u(j,k)
                  endif
               enddo
            endif
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_ubdry (u,x,y,b,nt)
c
c     Generate A*u on the boundary
c
c     Here, we generate values for the streamfunction (u) for
c     flow past a cylinder
c
c
      real    u(3,1),x(3,1),y(3,1)
      integer b(3,1)
c
      call rzero(u,3*nt)
c
      do k=1,nt
         do i=1,3
            if (abs(b(i,k)).eq.1) then
c              On circle
               u(i,k) = 0.
            elseif (b(i,k).ne.0) then
               u(i,k) = y(i,k)
            endif
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine init_tri(x,y,b,nt,xp,yp)
      real    x(3,1),y(3,1),xp(1),yp(1)
      integer b(3,1)
c
c     This routine initializes a cylinder of radius r=0.5
c
      nc =2
      r = 0.5
c
      one = 1.
      pi2 = 2.*atan(one)
      dtheta = pi2/nc
      do i=0,nc
         theta = i*dtheta
         xp(i+1) = r*cos(theta)
         yp(i+1) = r*sin(theta)
      enddo
c
      xp(4) = 2.50*r
      xp(5) = 1.50*r
      xp(6) = xp(2)
      xp(7) = 0.00
c
      yp(4) = 0.00
      yp(5) = yp(2)
      yp(6) = 1.50*r
      yp(7) = 2.50*r
c
      xp(8) = xp(4)
      yp(8) = yp(7)
c
      x(1,1) = xp(1)
      x(2,1) = xp(5)
      x(3,1) = xp(2)
c
      x(1,2) = xp(2)
      x(2,2) = xp(6)
      x(3,2) = xp(3)
c
      x(1,3) = xp(1)
      x(2,3) = xp(4)
      x(3,3) = xp(5)
c
      x(1,4) = xp(2)
      x(2,4) = xp(5)
      x(3,4) = xp(6)
c
      x(1,5) = xp(3)
      x(2,5) = xp(6)
      x(3,5) = xp(7)
c
      x(1,6) = xp(5)
      x(2,6) = xp(4)
      x(3,6) = xp(8)
c
      x(1,7) = xp(6)
      x(2,7) = xp(5)
      x(3,7) = xp(8)
c
      x(1,8) = xp(7)
      x(2,8) = xp(6)
      x(3,8) = xp(8)
c
      y(1,1) = yp(1)
      y(2,1) = yp(5)
      y(3,1) = yp(2)
c
      y(1,2) = yp(2)
      y(2,2) = yp(6)
      y(3,2) = yp(3)
c
      y(1,3) = yp(1)
      y(2,3) = yp(4)
      y(3,3) = yp(5)
c
      y(1,4) = yp(2)
      y(2,4) = yp(5)
      y(3,4) = yp(6)
c
      y(1,5) = yp(3)
      y(2,5) = yp(6)
      y(3,5) = yp(7)
c
      y(1,6) = yp(5)
      y(2,6) = yp(4)
      y(3,6) = yp(8)
c
      y(1,7) = yp(6)
      y(2,7) = yp(5)
      y(3,7) = yp(8)
c
      y(1,8) = yp(7)
      y(2,8) = yp(6)
      y(3,8) = yp(8)
c
      nt = 8
c     return
c
c     Expand box several times 
c
      r_b = xp(4)
c     fac = 2.0/r_b
      fac = 1.7/r_b
      do ib=1,2
         call qbox(x(1,nt+1),y(1,nt+1),nt,r_b,fac,xp,yp)
      enddo
c
c     copy and rotate by pi/2 new geometry
c
      call copy  (x(1,nt+1),x,3*nt)
      call copy  (y(1,nt+1),y,3*nt)
      call vrotate(x(1,nt+1),y(1,nt+1),3*nt,pi2)
      nt = 2*nt
c
      call izero(b,3*nt)
c
c     Set boundary flags on domain edge
c
      xmax = 0.99*glmax(x,3*nt)
c
      do iv=1,3*nt
         if ( abs(x(iv,1)).ge.xmax .or. abs(y(iv,1)).ge.xmax )
     $      b(iv,1) = 2
         if ( abs(y(iv,1)).le.0.001) b(iv,1) = 2
      enddo
c
c     Set boundary flags on circle
c
      rr = 1.01*r*r
      do iv=1,3*nt
         r2 = x(iv,1)**2 + y(iv,1)**2
         if (r2.le.rr) b(iv,1) = 1
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine refine_tri(xo,yo,bo,xn,yn,bn)
c
c     Refine the input triangle (x,y) into 4 triangles.  
c     Original triangle becomes center triangle.
c
c     x ,y  - vertices of input triangle  (becomes ctr. tri. on exit).
c     xn,yn - vertices (3 pair) of new triangles.
c     
      real    xo(3)  ,yo(3)
      real    xn(3,3),yn(3,3)
      integer bo(3)  ,bn(3,3)
c
c     First 3 new triangles...
c
      do k=1,3
c
c        1st vertex
c
         xn(1,k) = xo(k)
         yn(1,k) = yo(k)
         bn(1,k) = bo(k)
c
c        2nd and 3rd vertices
c
         do iv=2,3
            k0 = k + (iv-2)
            k1 = mod(k0,3) + 1
            if (abs(bo(k)).eq.1.and.abs(bo(k1)).eq.1) then
c
c             both points are on curve, therefore project midpoint...
c
              call proj_mid(xn(iv,k),yn(iv,k),xo(k),yo(k),xo(k1),yo(k1))
              bn(iv,k) = bo(k)
            else
c
c             interior point....
c
              xn(iv,k) = 0.5*(xo(k) + xo(k1))
              yn(iv,k) = 0.5*(yo(k) + yo(k1))
              if (bo(k).eq.0.or.bo(k1).eq.0) then
                 bn(iv,k) = 0
              elseif (abs(bo(k)).gt.abs(bo(k1))) then
                 bn(iv,k) = bo(k)
              else
                 bn(iv,k) = bo(k1)
              endif
c
            endif
         enddo
      enddo
c
c     Replace original triangle with center triangle!
c
c        4th center triangle.  Vertices are 2nd point of each triangle
c        4th triangle overwrites the incoming "old" triangle...
c
      do iv=1,3
         xo(iv) = xn(2,iv)
         yo(iv) = yn(2,iv)
         bo(iv) = bn(2,iv)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ircopy(r,i,n)
      real r(1)
      integer i(1)
      do k=1,n
         r(k) = i(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine velocity(vx,vy,psi,x,y,e2v,nt,ux,uy,n)
c
c     Compute ux =  dpsi/dy
c             uy = -dpsi/dx
c
c     Use formulas from Celia & Gray (3.5.21), (3.5.24) on p. 158.
c
      real    vx(3,1),vy(3,1),psi(3,1),x(3,1),y(3,1),ux(1),uy(1)
      integer e2v(3,1)
c
      call rzero(ux,n)
      call rzero(uy,n)
      call rzero(vx,n)
c
      do k=1,nt
         ak   = tri_area(x(1,k),y(1,k))
         ak2i = 0.5/ak
         bl3  = ak/3.0
         uxk = 0.
         uyk = 0.
         do i=1,3
            i2 = mod(i  ,3) + 1
            i3 = mod(i+1,3) + 1
            uxk = uxk + psi(i,k)*(x(i3,k)-x(i2,k))*ak2i
            uyk = uyk - psi(i,k)*(y(i2,k)-y(i3,k))*ak2i
         enddo
         uxkm = max(uxkm,abs(uxk))
         uykm = max(uykm,abs(uyk))
c
c        Add contribution from this element to each vertex
         do i=1,3
            ii = e2v(i,k)
            ux(ii)   = ux(ii)   + uxk*bl3
            uy(ii)   = uy(ii)   + uyk*bl3
c
c           Compute "lumped" mass matrix, temporarily stored in vx.
            vx(ii,1) = vx(ii,1) + bl3
c
         enddo
      enddo
c
      write(6,*) 'this is vx_max:',uxkm,uykm,nt
c
c
c     Divide through by lumped mass matrix
c
      do i=1,n
         uxo = ux(i)
         uyo = uy(i)
         ux(i) = ux(i)/vx(i,1)
         uy(i) = uy(i)/vx(i,1)
      enddo
c
c     map to "local" (triangle-based) reference
c
      do k=1,nt
         do i=1,3
            vx(i,k) = ux(e2v(i,k))
            vy(i,k) = uy(e2v(i,k))
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rsb_fem(order,list,ia,ja,e2v,nve,nel,bw,ldw,ndim)
c
c     Input: a list of elements with nve vertices per element
c
c     Output:  a sorted list of global vertices indicating optimal
c              nested dissection firing order
c     
c     bw is assumed to be a *big* work array.. perhaps 125*n in 3d, 25*n in 2d
c
c
      integer order(1),list(1),e2v(nve,1)
      integer ia(1),ja(1)
      real    bw(0:1)
      integer ptr(100),len(100),ptr2
c
c     Initialize list of elements
c
      do k=1,nel
         list(k) = k
      enddo
c
c
c     Set up dual of mesh graph given local-to-global pointer, e2v:
c
      ncn = (16/nve)**ndim * nel
      i0=0
      i1=i0 + nve*nel
      i2=i1 + nve*nel
      i3=i2 + ncn
      i4=i3 + ncn
      i5=i4 + ncn
c
      call set_dual_graph(ia,ja,list,1,nel,e2v,nve,ndim
     $                   ,bw(i0),bw(i1),bw(i2),bw(i3),bw(i4),bw(i5))
c    $                   ,gv    ,ge    ,iel   ,jel   ,ind   ,iw    )
c
c     Initialize Lanczos vector
c
      call init_lanczos(bw(1),bw(nel+1),ia,ja,nel)
c
c     Begin recursive bi-section
c
      nnz = 0
      levmax = 2*nel
      iptr = 1
      ptr(iptr) = 1
      len(iptr) = nel
      ldb = ldw-3*nel
      do lev=1,levmax
         call partition_list(list(ptr(iptr)),ia,ja,bw(ptr(iptr))
     $             ,ptr(iptr),len(iptr)
     $             ,ptr2,len2,e2v,nve,nd,bw(nel+1),bw(2*nel+1),ldb)
         call color_vertices(order,nnz,list(ptr(iptr))
     $                      ,len(iptr),len2,e2v,nve)
         if (len2.eq.0) then
            iptr = iptr-1
         else
            iptr = iptr+1
            ptr(iptr) = ptr2
            len(iptr) = len2
         endif
         if (iptr.eq.0) goto 10
      enddo
      write(6,*) 'exited loop',iptr
   10 continue
      write(6,*) 'went to 10 ',lev,nel
c
c     Reset order
c
      return
      end
c-----------------------------------------------------------------------
      subroutine partition_list(list,ia,ja,v1,ptr,n,pt2,ln2,e2v,nv,nd,v2
     $                         ,bw,ldw)
      integer list(1),e2v(nv,1),ptr,n,pt2,ln2
      real v1(1),v2(1),bw(0:1)
c
c     Input:   list[ptr:ptr+len]  - determines active d.o.f.s
c
c     Output:  list[ptr:ptr+len],list[pt2,pt2:ln2]
c              list is partitioned into two sets of length len & ln2
c             
c
c     Check for quick return
c
      if (n.eq.1) then
         ln2 = 0
         return
      elseif (n.eq.2) then
         n   = 1
         ln2 = 1
         pt2 = ptr+n
         return
      elseif (n.eq.3) then
c        Check for weird disconnected case
         call diag_check(idiag,bw,bw(n),list,n,ia,ja)
         if (idiag.eq.1) then
            n   = 1
            ln2 = 2
            pt2 = ptr+n
            return
         elseif (idiag.eq.3) then
            n   = 2
            ln2 = 1
            pt2 = ptr+n
            return
         elseif (idiag.eq.2) then
            n   = 2
            ln2 = 1
            pt2 = ptr+n
            itmp    = list(2)
            list(2) = list(3)
            list(3) = itmp
            return
         endif
      endif
c
c     RSB on dual graph
c
c
      ldr   = ldw - 8*n
      niter = ldr/n
      niter = min(niter,100)
      niter = min(niter,n)
      call perturb(v1,n,0.03)
      call ortho1 (v1,n)
      call lanczos1 (bw(7*n),n,bw(5*n),bw(6*n),niter,v1,list,ia,ja
     $              ,bw(0*n),bw(1*n),bw(2*n),bw(3*n),ierr)
      write(6,*) 'NITER',niter,n
c
c     Construct fiedler vectors 1 & 2
c
      if (ierr.eq.0) then 
         mn = 7+niter
         call lanczos2(v1,v2,bw(7*n),n,bw(mn*n),bw(5*n),bw(6*n),niter)
      else
c        p corresponds to a zero eigenvector of A, use it as Fiedler vector
         call copy(v1,bw(3*n),n)
         call copy(v2,bw(3*n),n)
         scale = vlsc2(v2,v2,n)
         scale = 1./sqrt(scale)
         call  cmult(v2,scale,n)
      endif
c
c     Sort list according to Fiedler vector
c
      call  rank(v1  ,bw(n),n)
      call  swap(v1  ,bw(n),n,bw(2*n))
c
c     Copy 2nd Fiedler vector to first for good initial guess...
c     In addition, perturb it slightly so that Lanczos will require
c     more than O(1) iterations....
c
      call  copy  (v1  ,v2  ,n)
      call  swap  (v1  ,bw(n),n,bw(2*n))
c
c     Swap list according to Fiedler vector
c
      call iswap  (list,bw(n),n,bw(2*n))
      m   = n
      ln2 = n/2
      n   = n-ln2
      pt2 = ptr+n
c
c     Sort sublists so each is in ascending vertex order
c
      call  irank(list,bw(m),n)
      call  iswap(list,bw(m),n,bw(2*m))
      call   swap(v1  ,bw(m),n,bw(2*m))
c
      call  irank(list(n+1),bw(m),ln2)
      call  iswap(list(n+1),bw(m),ln2,bw(2*m))
      call   swap(v1  (n+1),bw(m),ln2,bw(2*m))
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lanczos1(rr,n,diag,upper,niter,f,list,ia,ja
     $                   ,r,z,w,p,ierr)
c
      real    rr(n,0:1),diag(1),upper(1),f(1),r(1),p(1),w(1),z(1)
      integer list(1),ia(1),ja(1)
c
      call rzero(diag ,niter)
      call rzero(upper,niter)
      pap  = 1.0
      ierr = 0
c
c     set machine tolerances
c
      one = 1.
      eps = 1.e-28
      if (one+eps .eq. one) eps = 1.e-13
      if (one+eps .eq. one) eps = 1.e-6
      eps2 = eps*eps
      eps1 = sqrt(eps)
c
      rtz1=1.0
c
      call copy  (r,f,n)
      call ortho1(r,n)
      rtr = vlsc2(r,r,n)
      rlim2 = rtr*eps**2
      rnorm = sqrt(rtr)
c
      iter = 0
      rnrmi = 1./rnorm
      call cmult2 (rr(1,iter),r,rnrmi,n)
c
      write(6 ,6) n,rnorm
c
c     call graph_out(w,p,list,n,ia,ja)
c
      niter=min(n,niter)
      miter=min(n,niter)
      do 1000 iter=1,miter
c
c        Invert preconditioner here
         call copy(z,r,n)
c
         rtz2=rtz1
         rtz1=vlsc2(r,z,n)
c
         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0
c
         call add2s1(p,z,beta,n)
         call ortho1(p,n)
c
         call ax(w,p,list,n,ia,ja)
c
c        Save p^Ap for eigenvalue estimates
         pap_old = pap
         pap=vlsc2(w,p,n)
         if (iter.eq.1) pap0 = pap
         if (pap.le.eps2.or.
     $      (pap/pap_old.le.eps1.and.iter.gt.1)) then
            write(6 ,6) iter,rnorm,pap_old,pap
            niter = iter-1
            ierr  = 1
            return
         endif
         if (pap/pap0.gt.1.e6) then
            write(6 ,6) iter,rnorm,pap_old,pap
c           call vout(p,n,0,' p  ')
c           call vout(w,n,0,'Ap  ')
            call graph_out(w,r,list,n,ia,ja)
c
c           write(6,*) 'quit in lanczos1',pap,iter,n
c           call exit
c
            niter = iter-1
            ierr  = 1
            return
         endif
c
         alpha=rtz1/pap
         alphm=-alpha
c        call add2s2(x,p,alpha,n)    No solution needed for Lanczos
         call add2s2(r,w,alphm,n)
         call ortho1(r,n)
c
         rtr   = vlsc2(r,r,n)
         rnorm = sqrt(rtr)
         rnrmi = 1./rnorm
         call cmult2 (rr(1,iter),r,rnrmi,n)
c        write(6,6) iter,rnorm,alpha,beta,pap
         if (iter.le.10.or.mod(iter,10).eq.0) write(6,6) iter,rnorm
     $                                              ,alpha,beta,pap
    6    format('cg:',i4,1p4e12.4)
c
c        Generate tridiagonal matrix for Lanczos scheme 
         if (iter.eq.1) then
            diag(iter) = pap/rtz1
         else
            diag(iter)    = (beta**2 * pap_old + pap ) / rtz1
            upper(iter-1) = -beta * pap_old / sqrt(rtz2 * rtz1)
         endif
         if (rtr.le.rlim2) goto 1001
 1000 continue
      iter = iter-1
 1001 continue
c
c     niter = iter-1
      niter = iter
c     Check for "complete" convergence...
c     if(rnorm.lt.eps) niter = niter-1
      write(6,6) iter,rnorm,rlim2,rtr
c
c     Call eigenvalue routine for Lanczos scheme:
c
c     call calcz(diag,upper,niter,dmax,dmin,ev)  Done outside
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ax(v,u,list,n,ia,ja)
c
c     This routine computes v = Lu, where L is the graph Laplacian of A,
c     and u and v are sparse vectors with entries given in "list".
c
c
      real v(1),u(1)
      integer list(1),ia(1),ja(1)
c
      do i=1,n
         ii = list(i)
         v(i) = sparse_dot_glap (ia,ja,ii,u,list,n)
      enddo
      return
      end
c-----------------------------------------------------------------------
      function sparse_dot_glap (ia,ja,i,b,jb,nb)
c
c     Compute dot product of row i of graph laplacian with vector b
c     1...k.
c
      real b(1)
      integer ia(1),ja(1),jb(1)
      integer ca,cb,pa,pb
c
      na = ia(i+1)-ia(i)
      kmin = na+nb
c
      pa = ia(i)
      pb = 1
c
      ca = ja(pa)
      cb = jb(pb)
c
      ma = 1
      mb = 1
      noff_diag = 0
      diag      = 0.
      dot       = 0.
      do jj=1,kmin
         if (ca.lt.cb.and.ma.le.na) then
            ma = ma+1
            pa = pa+1
            ca = ja(pa)
         elseif (cb.lt.ca.and.mb.le.nb) then
            mb = mb+1
            pb = pb+1
            cb = jb(pb)
         elseif (ma.gt.na.or.mb.gt.nb) then
            goto 1
         else
            if (ja(pa).eq.i) then
               diag = b(pb)
            else
               dot = dot - b(pb)
               noff_diag = noff_diag+1
            endif
            ma = ma+1
            mb = mb+1
            pa = pa+1
            pb = pb+1
            ca = ja(pa)
            cb = jb(pb)
         endif
      enddo
    1 continue
      sparse_dot_glap = dot + noff_diag*diag
      return
      end
c-----------------------------------------------------------------------
      subroutine lanczos2(v1,v2,rr,n,ev,diag,upper,m)
c
c     Compute two smallest eigenvector estimates 
c
      real    v1(1),v2(1),rr(n,1),ev(m,1),diag(1),upper(1)
c
      call calcz(diag,upper,m,dmax,dmin,ev)
c
c     Find dmin1 & dmin2
c
      dmin1 = diag(1)+diag(m)
      dmin2 = diag(2)+diag(m)
c
      do i=1,m
         if (diag(i).lt.dmin1) then
            i1 = i
            dmin1 = diag(i)
         endif
      enddo
c
      do i=1,m
         if (diag(i).lt.dmin2.and.i.ne.i1) then
            i2 = i
            dmin2 = diag(i)
         endif
      enddo

c
      call rzero(v1,n)
      call rzero(v2,n)
      do j=1,m
         call add2s2(v1,rr(1,j),ev(j,i1),n)
         call add2s2(v2,rr(1,j),ev(j,i2),n)
         s  = glsum(rr(1,j),n)
         sv = glsum(v1,n)
c        write(6,3) n,i1,i2,j,' ev:',s,sv,v1(n),rr(n,j),ev(j,1),ev(j,2)
      enddo
   3  format(4i4,a4,1p8e12.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ortho1(p,n)
c
c     Orthogonalize wrt the 1 vector
c
      real p(1)
c
      s = 0.
      do i=1,n
         s = s + p(i)
      enddo
      s = s/n
      do i=1,n
         p(i) = p(i) - s
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_dual_graph(ia,ja,list,ptr,nel,e2v,nve,nd
     $                         ,gv,ge,iel,jel,ind,iw)
      integer ja(1),ia(1),list(1),e2v(nve,1),ptr,nel
      integer gv(1),ge(1),iel(1),jel(1),ind(1),iw(1)
c
c     Work arrays:
c
c     integer gv(nve*nel),ge(nve*nel)
c     integer jel(ncan),iel(ncan)    ncan ~ nel*(8/nve)**nd
c     integer ind(ncan),iw (ncan)    ncan ~ nel*(8/nve)**nd
c
c
c     Construct element centroid-based graph (dual of vertex graph)
c     Graph is based upon shared edges (faces in 3d)
c
      do i=1,nel
      enddo
c
      nvg = 0
      do ipt=0,nel-1
         i  = ptr+ipt
         ie = list(i)
         do il=1,nve
            nvg = nvg+1
            gv(nvg) = e2v(il,ie)
            ge(nvg) = ie
         enddo
      enddo
c
c     Sort by global vertex number
c
      call irank(gv,ind,nvg)
      call iswap(gv,ind,nvg,iw)
      call iswap(ge,ind,nvg,iw)
c
c
c     Form candidate ie,je pairings
c
      ncan  = 0
      jlast = 1
      do iv=1,nvg
         do j=jlast+1,nvg
c           scan for change in gv
            if (gv(j).ne.gv(jlast)) goto 20
         enddo
         j=nvg+1
   20    continue
         n_in_seg = j-jlast
         call irank(ge(jlast),ind,n_in_seg)
         call iswap(ge(jlast),ind,n_in_seg,iw)
c
         do j2=jlast,j-1
         do j1=jlast,j-1
            if (j1.ne.j2) then
               ncan = ncan+1
               iel(ncan) = ge(j2)
               jel(ncan) = ge(j1)
            endif
         enddo
         enddo
c
         jlast = j
         if (jlast.gt.nvg) goto 30
c
      enddo
   30 continue
c
c
c     Sort candidate pairings, keep only if a number is paired nd 
c     times or more
c
      call irank(iel,ind,ncan)
      call iswap(iel,ind,ncan,iw)
      call iswap(jel,ind,ncan,iw)
c
      nnz  =0
      call izero(ia,nel)
c
      jlast=1
      do ie=1,nel
c
c        Add "ie" to graph...  (self-referential)
c
         nnz     = nnz+1
         ja(nnz) = ie
         ia(ie)  = ia(ie)+1
c
         do j=jlast+1,ncan
c           scan for change in iel
            if (iel(j).ne.iel(jlast)) goto 40
         enddo
         j=ncan+1
   40    continue
c
         n_in_seg = j-jlast
         call irank(jel(jlast),ind,n_in_seg)
         call iswap(jel(jlast),ind,n_in_seg,iw)
c
c        Scan candidate pairing for more than nd matches
c
         nmatch=1
         do i=1,n_in_seg
            if (jel(jlast+i).eq.jel(jlast+i-1).and.i.lt.n_in_seg) then
               nmatch = nmatch+1
            else
c
c              if (nmatch.ge.nd) then   (In order to get connected graphs, I had to go to the
c                                        below.  Actually, it's more correct!  (1/3/98 pff).
               if (nmatch.ge.1) then
                  nnz     = nnz+1
                  ja(nnz) = jel(jlast+i-1)
                  ia(ie)  = ia(ie)+1
               endif
               nmatch = 1
            endif
         enddo
         jlast = j
      enddo
c
c     Reset ia to be a pointer rather than a length counter
c
      do i=nel+1,2,-1
         ia(i) = ia(i-1)
      enddo
      ia(1) = 1
      do i=2,nel+1
         ia(i) = ia(i)+ia(i-1)
      enddo
c
c     Sort columns within each row
c
      do i=1,nel
         na = ia(i+1) - ia(i)
         call irank(ja(ia(i)),ind,na)
         call iswap(ja(ia(i)),ind,na,iw)
      enddo
c
c
      do i=1,nel
         nn = ia(i+1) - ia(i)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine iout(x,n,num,name4)
      integer x(1)
      character*4 name4
      m = min(n,40)
      write(6,1) num,n,name4,(x(k),k=1,m)
    1 format(2i4,1x,a4,26i3,5(/,13x,26i3))
      return
      end
c-----------------------------------------------------------------------
      subroutine vout(x,n,num,name4)
      real x(1)
      character*4 name4
      m = min(n,40)
      write(6,1) num,n,name4,(x(k),k=1,m)
    1 format(2i4,1x,a4,16f6.3,5(/,13x,16f6.3))
      return
      end
c-----------------------------------------------------------------------
      subroutine color_vertices(order,nnz,list,n1,n2,e2v,nv)
      integer order(1),list(1),e2v(nv,1)

c
c     Set flag on list1
c
      do il=1,n1
        ie = list(il)
        do iv=1,nv
           ig=e2v(iv,ie)
           if (order(ig).eq.0) order(ig)=-1
        enddo
      enddo
c
      if (n2.gt.0) then
c       Compare flag for list2
        do il=1,n2
           ie = list(n1+il)
           do iv=1,nv
              ig=e2v(iv,ie)
              if (order(ig).eq.-1) then
                 nnz = nnz+1
                 order(ig)=nnz
              endif
           enddo
         enddo
      else
c       check flag on list1
        do il=1,n1
           ie = list(il)
           do iv=1,nv
              ig=e2v(iv,ie)
              if (order(ig).eq.-1) then
                 nnz = nnz+1
                 order(ig)=nnz
              endif
           enddo
         enddo
      endif
c
c     UnSet flag on list1
c
      do il=1,n1
        ie = list(il)
        do iv=1,nv
           ig=e2v(iv,ie)
           if (order(ig).eq.-1) order(ig)=0
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine graph_out(w,p,list,n,ia,ja)
      real w(1),p(1)
      integer list(1),ia(1),ja(1)
      call rzero(p,n)
c
      write(6,1) n,'axn ',(list(k),k=1,n)
    1 format('list',i4,1x,a4,26i5,5(/,13x,26i5))
c
      do i=1,n
         p(i) = 1.
         call ax(w,p,list,n,ia,ja)
         call gout(w,n,i,'axn ')
         call bout(w,n,i,'axn ')
         p(i) = 0.
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gout(x,n,num,name4)
      real x(1)
      character*4 name4
      m = min(n,40)
      write(6,2) num,n,name4,(x(k),k=1,m)
    2 format(2i4,1x,a4,26f5.1,5(/,13x,26f5.1))
      return
      end
c-----------------------------------------------------------------------
      subroutine bout(x,n,num,name4)
      real x(1)
      character*4 name4
c
      write(n,2) (x(k),k=1,n)
    2 format(1pe20.12)
      return
      end
c-----------------------------------------------------------------------
      subroutine a_trid(v,u,n)
      real v(1),u(1)
c
      do i=2,n-1
         v(i) = -( u(i-1) + u(i+1) - 2.*u(i))
      enddo
      v(1) = -(u(  2) - 1.*u(1))
      v(n) = -(u(n-1) - 1.*u(n))
      return
      end
c-----------------------------------------------------------------------
      subroutine perturb(v,n,amp)
      real v(1)
c
c     Perturb v() with pseudo random variable.
c
      a = amp*(glmax(v,n)-glmin(v,n))
c
      do i=1,n
         x = i*1.e5
c        v(i) = v(i) + a*cos(x)
c        v(i) = v(i)*(1. + a*cos(x))
         v(i) = v(i)*(1. + a*cos(x)) + .1*a*sin(x)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine diag_check(id,w,p,list,n,ia,ja)
c
c     Check for first zero on diagonal of graph laplacian
c
      real w(1),p(1)
      integer list(1),ia(1),ja(1)
      call rzero(p,n)
c
      id = 0
      do i=1,n
         p(i) = 1.
         call ax(w,p,list,n,ia,ja)
         p(i) = 0.
         if (w(i).eq.0) then
            id = i
            return
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine compress_int(list,n,ind,w)
      integer list(1),ind(1),w(1)
c
c     Input:   list() -- a sequence of integers, possibly w/ repeats
c     Output:  list() -- entries in list_out preserve the rank of
c                        entries in list_in, but list is now contiguous
c
      call irank(list,ind,n)
      call iswap(list,ind,n,w)
c
      last = list(1)
      nn = 1
      w(1) = nn
c
      do i=2,n
         if (list(i).ne.last) then
            last = list(i)
            nn = nn+1
         endif
         w(i) = nn
      enddo
c
      call icopy   (list,w,n)
      call iunswap (list,ind,n,w)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_rc(ir,ic,nl,e2v,b,nt,ind,w,type)
c
c     Generate (r,c) pairs for sparse matrix.
c
c     On input,  nl is the amt. of space allocated for ir & ic.
c     On return, nl is the amt. of space required  for ir & ic.
c
      integer ir(1),ic(1)
      integer e2v(3,1),b(3,1)
      integer ind(1),w(1)
      character*1 type
c
c
c     First, get (r,c) pairs... upper triangle only, no diagonal.
c
      l = 0
      do k=1,nt
         do j=1,3
            if (b(j,k).le.0) then
c              In interior of Omega...
               jj = e2v(j,k)
               do i=j+1,3
                  if (b(i,k).le.0) then
c                    In interior of Omega...
                     ii = e2v(i,k)
                     l  = l+1
                     if (l.gt.nl) then
                        write(6,*) 'Error increase nl in gen_rc.'
                        write(6,*) 'ABORT  ijk:',i,j,k,nt,nl
                        call exit
                     endif
c
                     if (type.eq.'u') then
c                       Format for storing upper half
                        ir(l) = min(ii,jj)
                        ic(l) = max(ii,jj)
                     else
c                       Format for storing lower half
                        ir(l) = max(ii,jj)
                        ic(l) = min(ii,jj)
                     endif
                  endif
               enddo
            endif
         enddo
      enddo
      nl = l
c
c     Now sort and compress out redundant entries
c
      call irank(ir,ind,nl)
      call iswap(ir,ind,nl,w)
      call iswap(ic,ind,nl,w)
c
c     Now for each row...
c
      istart = 1
      irs    = ir(istart)
      do i=2,nl
         if (ir(i).ne.irs) then
            nrs = i-istart
            call irank(ic(istart),ind,nrs)
            call iswap(ic(istart),ind,nrs,w)
            istart = i
            irs    = ir(istart)
         elseif (i.eq.nl) then
            nrs = i-istart+1
            call irank(ic(istart),ind,nrs)
            call iswap(ic(istart),ind,nrs,w)
         endif
      enddo
c
c     Eliminate redundant pairs
c
      l = 1
      do i=2,nl
         if (ir(i).ne.ir(l) .or. ic(i).ne.ic(l) ) then
            l = l+1
            ir(l) = ir(i)
            ic(l) = ic(i)
         endif
      enddo
      nl = l
c     do l=1,nl
c        write(6,*) 'irc:',l,ir(l),ic(l)
c     enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_graph_l(ia,ja,nr,e2v,b,nt,ir,ic,ind,w,ldw)
c
c     Generate ia() and ja() for lower half stiffness matrix
c
      integer ia(1),ja(1)
      integer e2v(3,1),b(3,1)
      integer ir(1),ic(1),ind(1),w(1)
c
c
c     First, get (r,c) pairs... lower triangle only, no diagonal.
c
      nl = ldw
      call gen_rc(ir,ic,nl,e2v,b,nt,ind,w,'l')
c
c     Convert to crs format
c
      nr    = 1
      ia(1) = 1
      nc    = 0
      l     = 1
c
      nrow1 = iglmax(ir,nl)
      nrow2 = iglmax(ic,nl)
      nrows = max(nrow1,nrow2)
      do i=1,nl+nrows
         if (ir(l).ne.nr) then
c           New row: put diagonal element in ja() of current row
            nc     = nc+1
            ja(nc) = nr
c           new row -- 
            nr     = nr+1
            ia(nr) = nc+1
         endif
         if (ir(l).eq.nr) then
            nc     = nc+1
            ja(nc) = ic(l)
            l      = l+1
            if (l.gt.nl) goto 3
         endif
      enddo
    3 continue
    5 format(a4,12i6)
c
c     Last row: put diagonal element in ja(n)
      nc       = nc+1
      ja(nc)   = nr
      ia(nr+1) = nc+1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_graph_u(ia,ja,nr,e2v,b,nt,ir,ic,ind,w,ldw)
c
c     Generate ia() and ja() for augmented stiffness matrix
c
      integer ia(1),ja(1)
      integer e2v(3,1),b(3,1)
      integer ir(1),ic(1),ind(1),w(1)
c
c
c     First, get (r,c) pairs... upper triangle only, no diagonal.
c
      nl = ldw
      call gen_rc(ir,ic,nl,e2v,b,nt,ind,w,'u')
c
c     Convert to crs format
c
      nr    = 1
      ia(1) = 1
      nc    = 1
      ja(1) = 1
c
      do i=1,nl
         if (ir(i).ne.nr) then
c           new row -- 
            nr     = nr+1
            ia(nr) = nc+1
c           put diagonal element in ja()
            nc     = nc+1
            ja(nc) = nr
         endif
         nc     = nc+1
         ja(nc) = ic(i)
      enddo
C
c     Last row: put diagonal element in ja(n)
      nr     = nr+1
      ia(nr) = nc+1
c     put diagonal element in ja(n)
      nc     = nc+1
      ja(nc) = nr
c
      ia(nr+1) = nc+1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_acrs (a,ia,ja,nr,x,y,e2v,b,nt)
c
c     Generate CRS based stiffness matrix
c
c     Note:  the "1's" in the next 2 lines can also be n and nt, resp.
      real a(1),x(3,1),y(3,1)
      integer ia(1),ja(1)
      integer b(3,1),e2v(3,1)
c
c     Initialize a:
      nnz = ia(nr+1)-ia(1)
      call rzero (a,nnz)
c
      do k=1,nt
         Ak4 = 4.0*tri_area(x(1,k),y(1,k))
         do j=1,3
            if (b(j,k).le.0) then
c
c              in interior for j...
               jj = e2v(j,k)
               j1 = mod(j  ,3) + 1
               j2 = mod(j+1,3) + 1
               dxj = x(j1,k) - x(j2,k)
               dyj = y(j1,k) - y(j2,k)
c
               do i=1,3
                  if (b(i,k).le.0) then
c                    in interior for i...
c
                     ii = e2v(i,k)
                     if (jj.le.ii) then
c                       only compute lower-half of A, due to symmetry...
c
                        i1 = mod(i  ,3) + 1
                        i2 = mod(i+1,3) + 1
                        dxi = x(i1,k) - x(i2,k)
                        dyi = y(i1,k) - y(i2,k)
                        da  =  (dxi*dxj + dyi*dyj)/Ak4
c
                        call acrs_index(icrs,ii,jj,ia,ja)
                        a(icrs)=a(icrs) + da
                     endif
                  endif
               enddo
            endif
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine acrs_index(icrs,ii,jj,ia,ja)
      integer ia(1),ja(1)
c
      icrs = 0
c
      i0 = ia(ii)
      i1 = ia(ii+1)
c
      do i=i0,i1-1
         if (ja(i).eq.jj) then
            icrs = i
            return
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outjmat(ja,ia,n,name4)
      integer ja(1),ia(1)
      character*4 name4
c
      write(6,1) name4,n
    1 format(/,'integer matrix: ',a4,2i5)
      do i=1,n
         i0  = ia(i)
         i1  = ia(i+1) - 1
         write(6,2) i,(ja(j),j=i0,i1)
    2    format(i3,3x,40i4)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine init_lanczos(lanc,rc,ia,ja,n)
c
c     Use reverse Cuthil-McKee to generate initial guess for Fiedler vector
c
      real    lanc(1)
      integer rc(2,0:1)
      integer ja(1),ia(1)
      integer sm_bandwidth_reduction
c
c     Fill array for hmt's RCM code:
c
c
      l = 0
      do i=1,n
         i0 = ia(i)
         i1 = ia(i+1)-1
         do j=i0,i1
c           Store strict upper half
            if (ja(j).gt.i) then
               l = l+1
               rc(1,l) = i
               rc(2,l) = ja(j)
            endif
         enddo
      enddo
      rc(1,0) = n
      rc(2,0) = l
c
      itype = -1
      ierr = sm_bandwidth_reduction(rc,lanc,itype)
      call icopy(rc,lanc,n)
c
c     This is the way A should be reordered:
c
c     do j=1,n
c     do i=1,n
c        ii = rc(i,0)
c        jj = rc(j,0)
c        a_min_bw(i,j) = a(ii,jj)
c     enddo
c     enddo
c
c     Thus, row rc(i) moves to row i...
c     So, we initialize the Lanczos (Fiedler) vector
c     according to this permutation
c
      do i=1,n
         lanc(rc(i,0)) = i
      enddo
c
      return
      end
c-----------------------------------------------------------------------
