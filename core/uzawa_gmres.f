      subroutine uzawa_gmres(res,h1,h2,h2inv,intype,iter)
c-----------------------------------------------------------------------
c
c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)
c
c-----------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'GMRES'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx2*ly2*lz2*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             h2inv(lx1,ly1,lz1,lelv)

      real wp(lx2,ly2,lz2,lelv)

      common /ctmp0/ wk1(lgmres),wk2(lgmres)

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
      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml,mu,bm2,bm2inv,nx2*ny2*nz2*nelv)
         norm_fac = 1./sqrt(volvm2)
      endif
c
      etime1 = dnekclock()
      divex = 0.
      iter  = 0
      m = lgmres
c
      call chktcg2(tolps,res,iconv)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
c     if (param(21).lt.0) tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      ntot2  = nx2*ny2*nz2*nelv
c
      iconv = 0
      call rzero(x,ntot2)
      do while(iconv.eq.0.and.iter.lt.200)
         if(iter.eq.0) then
                                          !      -1
            call col3(r,ml,res,ntot2)     ! r = L  res
c            call copy(r,res,ntot2)
         else
            !update residual
            call copy(r,res,ntot2)                ! r = res
            call cdabdtp(w,x,h1,h2,h2inv,intype)  ! w = A x
            call add2s2(r,w,-1.,ntot2)            ! r = r - w
                                                  !      -1
            call col2(r,ml,ntot2)                 ! r = L   r
         endif
                                                  !            ______
         gamma(1) = sqrt(glsc2(r,r,ntot2))        ! gamma  = \/ (r,r) 
                                                  !      1
         if(iter.eq.0) then
            div0 = gamma(1)*norm_fac
            if (param(21).lt.0) tolpss=abs(param(21))*div0
         endif

         !check for lucky convergence
         rnorm = 0.
         if(gamma(1) .eq. 0.) goto 9000
         temp = 1./gamma(1)
         call cmult2(v(1,1),r,temp,ntot2)         ! v  = r / gamma
                                                  !  1            1
         do j=1,m
            iter = iter+1
                                                  !       -1
            call col3(w,mu,v(1,j),ntot2)          ! w  = U   v
                                                  !           j
            if(param(43).eq.1) then
               call uzprec(z(1,j),w,h1,h2,intype,wp)
            else                                  !       -1
               call hsmg_solve(z(1,j),w)          ! z  = M   w
            endif                                 !  j        
     
            call cdabdtp(w,z(1,j),                ! w = A z
     $                   h1,h2,h2inv,intype)      !        j
     
                                                  !      -1
            call col2(w,ml,ntot2)                 ! w = L   w
            !modified Gram-Schmidt

c           do i=1,j
c
c              h(i,j)=glsc2(w,v(1,i),ntot2)       ! h    = (w,v )
c                                                 !  i,j       i
c              call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
c                                                 !          i,j  i
c           enddo

c           1-PASS GS, 1st pass:
            do i=1,j
               h(i,j)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
            enddo                                 !  i,j       i
            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs
            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
            enddo                                 !          i,j  i


            !apply Givens rotations to new column
            do i=1,j-1
               temp = h(i,j)                   
               h(i  ,j)=  c(i)*temp + s(i)*h(i+1,j)  
               h(i+1,j)= -s(i)*temp + c(i)*h(i+1,j)
            enddo
                                                  !            ______
            alpha = sqrt(glsc2(w,w,ntot2))        ! alpha =  \/ (w,w)
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
            call cmult2(v(1,j+1),w,temp,ntot2)   ! v    = w / alpha
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
            call add2s2(x,z(1,i),c(i),ntot2)     ! x = x + c  z
                                                 !          i  i
         enddo
c        if(iconv.eq.1) call dbg_write(x,nx2,ny2,nz2,nelv,'esol',3)
      enddo
 9000 continue
c
      divex = rnorm
c     iter = iter - 1
c
c     DIAGNOSTICS
c      call copy(w,x,ntot2)
c      if (ifvcor) then
c         xaver = glsc2(bm2,w,ntot2)/volvm2
c         call cadd(w,-xaver,ntot2)
c      endif
c      call copy(r,res,ntot2) !r = res
c      call cdabdtp(r,w,h1,h2,h2inv,intype)  ! r = A w
c      do i=1,ntot2
c         r(i) = res(i) - r(i)               ! r = res - r
c      enddo
c      call uzawa_gmres_temp(r,bm2inv,ntot2)
c                                               !            ______
c      gamma(1) = sqrt(glsc2(r,r,ntot2)/volvm2) ! gamma  = \/ (r,r) 
c                                               !      1
c      print *, 'GMRES end resid:',gamma(1)
c     END DIAGNOSTICS
      call copy(res,x,ntot2)
c
      if (ifvcor) then
         xaver = glsc2(bm2,res,ntot2)/volvm2
         call cadd(res,-xaver,ntot2)
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
