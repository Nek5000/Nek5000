      subroutine hmh_gmres(res,h1,h2,wt,iter)
c-----------------------------------------------------------------------
c
c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.
c
c     
c
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'GMRES'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx1*ly1*lz1*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             wt   (lx1,ly1,lz1,lelv)
c
      real /scrmg/ wp(lx1,ly1,lz1,lelv)
c
      real y(lgmres)
      real alpha, l, temp
      integer j,m
c
      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac
c
      real*8 etime1,dnekclock
c
      n = nx1*ny1*nz1*nelv

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml,mu,bm1,binvm1,nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif
c
      etime1 = dnekclock()
      divex = 0.
      iter  = 0
      m = lgmres
c
      call chktcg1(tolps,res,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      iconv = 0
      call rzero(x,n)
      do while(iconv.eq.0.and.iter.lt.200)
         if(iter.eq.0) then
                                          !      -1
            call col3(r,ml,res,n)         ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r,res,n)                  ! r = res
            call ax    (w,x,h1,h2)                ! w = A x
            call add2s2(r,w,-1.,n)                ! r = r - w
                                                  !      -1
            call col2(r,ml,n)                     ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc3(r,r,wt,n))         ! gamma  = \/ (r,r) 
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,n)             ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3(w,mu,v(1,j),n)              ! w  = U   v
                                                  !           j
                                                  !       -1
            call hsmg_solve_h1 (z(1,j),w)         ! z  = M   w
                                                  !  j        
     
            call ax  (w,z(1,j),h1,h2)             ! w = A x
                                                  !        j
     
                                                  !      -1
            call col2(w,ml,n)                     ! w = L   w
            !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  !  i,j       i
                                               
               call add2s2(w,v(1,i),-h(i,j),n)! w = w - h    v
                                                  !          i,j  i
            enddo
            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                 !            ______
            alpha = sqrt(glsc3(w,w,wt,n))        ! alpha =  \/ (w,w)
            rnorm = 0.
            if(alpha.eq.0.) goto 900  !converged
            l = sqrt(h(j,j)*h(j,j)+alpha*alpha)
            temp = 1./l
            c(j) = h(j,j) * temp
            s(j) = alpha  * temp
            h(j,j) = l
            gamma(j+1) = -s(j) * gamma(j)
            gamma(j)   =  c(j) * gamma(j)

c           call outmat(h,m,j,' h    ',j)
            
            rnorm = abs(gamma(j+1))*norm_fac
            ratio = rnorm/div0
            if (ifprint.and.nid.eq.0) 
     $         write (6,66) iter,tolpss,rnorm,div0,ratio,istep
   66       format(i5,1p4e12.5,i8,' Divergence')

            if (rnorm .lt. tolpss) goto 900  !converged
            if (j.eq.m) goto 1000 !not converged, restart
            
            temp = 1./alpha
            call cmult2(v(1,j+1),w,temp,n)   ! v    = w / alpha
                                             !  j+1            
         enddo
  900    iconv = 1
 1000    continue
         !back substitution
         !     -1
         !c = H   gamma
         do k=j,1,-1
            temp = gamma(k)
            do i=j,k+1,-1
               temp = temp - h(k,i)*c(i)
            enddo
            c(k) = temp/h(k,k)
         enddo
         !sum up Arnoldi vectors
         do i=1,j
            call add2s2(x,z(1,i),c(i),n)     ! x = x + c  z
                                             !          i  i
         enddo
c        if(iconv.eq.1) call dbg_write(x,nx1,ny1,nz1,nelv,'esol',3)
      enddo
 9000 continue
c
      divex = rnorm
      call copy(res,x,n)
c
      if (ifvcor) then
         xaver = glsc2(bm1,res,n)/volvm1
         call cadd(res,-xaver,n)
      endif
c
      etime1 = dnekclock()-etime1
      if (nid.eq.0) write(6,9999) istep,iter,divex,tolpss,div0,etime1
c     call flush_hack
 9999 format(I10,' U-Press gmres: ',I6,1p4E13.4)
19999 format(I10,' U-Press 1.e-5: ',I6,1p4E13.4)
c
c
      return
      end

      subroutine uzawa_gmres_split(l,u,b,binv,n)
      integer n
      real l(n),u(n),b(n),binv(n)
      integer i
      do i=1,n
         l(i)=sqrt(binv(i))
         u(i)=sqrt(b(i))
         if(abs(u(i)*l(i)-1.0).gt.1e-13) print *, i, u(i)*l(i)
      enddo
      return
      end

      subroutine uzawa_gmres_temp(a,b,n)
      integer n
      real a(n),b(n)
      integer i
      do i=1,n
         a(i)=sqrt(b(i))*a(i)
      enddo
      return
      end
      
