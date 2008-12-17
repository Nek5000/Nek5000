c-----------------------------------------------------------------------
      subroutine uzawa_gmres(res,h1,h2,h2inv,intype,iter)

c     Solve the pressure equation by right-preconditioned 
c     GMRES iteration.
c     intype =  0  (steady)
c     intype =  1  (explicit)
c     intype = -1  (implicit)

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

      common /scrmg/    wp (lx2,ly2,lz2,lelv)

      common /ctmp0/   wk1(lgmres),wk2(lgmres)
      common /cgmres1/ y(lgmres)

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
      do while(iconv.eq.0.and.iter.lt.100)
         if(iter.eq.0) then
                                                  !      -1
            call col3(r,ml,res,ntot2)             ! r = L  res
c           call copy(r,res,ntot2)
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
c              h(i,j)=glsc2(w,v(1,i),ntot2)       ! h    = (w,v )
c                                                 !  i,j       i
c              call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
c           enddo                                 !          i,j  i


c           2-PASS GS, 1st pass:

            do i=1,j
               h(i,j)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
            enddo                                 !  i,j       i

            call gop(h(1,j),wk1,'+  ',j)          ! sum over P procs

            do i=1,j
               call add2s2(w,v(1,i),-h(i,j),ntot2)! w = w - h    v
            enddo                                 !          i,j  i


c           2-PASS GS, 2nd pass:
c
c           do i=1,j
c              wk1(i)=vlsc2(w,v(1,i),ntot2)       ! h    = (w,v )
c           enddo                                 !  i,j       i
c                                                 !
c           call gop(wk1,wk2,'+  ',j)             ! sum over P procs
c
c           do i=1,j
c              call add2s2(w,v(1,i),-wk1(i),ntot2)! w = w - h    v
c              h(i,j) = h(i,j) + wk1(i)           !          i,j  i
c           enddo

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

c            call outmat(h,m,j,' h    ',j)
            
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

c-----------------------------------------------------------------------
      subroutine uzawa_gmres_split(l,u,b,binv,n)
      integer n
      real l(n),u(n),b(n),binv(n)
      integer i
      do i=1,n
c        l(i)=sqrt(binv(i))
c        u(i)=sqrt(b(i))

c        u(i)=sqrt(b(i))
c        l(i)=1./u(i)

c        l(i)=sqrt(binv(i))
         l(i)=1.
         u(i)=1./l(i)


c        if(abs(u(i)*l(i)-1.0).gt.1e-13)write(6,*) i,u(i)*l(i),' gmr_sp'
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine uzawa_gmres_temp(a,b,n)
      integer n
      real a(n),b(n)
      integer i
      do i=1,n
         a(i)=sqrt(b(i))*a(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,x,h1,h2,n)
      include 'SIZE'
      include 'TOTAL'

c
c     w = A*x for pressure iteration
c

      integer n
      real w(n),x(n),h1(n),h2(n)

      imsh = 1
      isd  = 1
      call axhelm (w,x,h1,h2,imsh,isd)
      call dssum  (w,nx1,ny1,nz1)
      call col2   (w,pmask,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine hmh_gmres(res,h1,h2,wt,iter)

c     Solve the Helmholtz equation by right-preconditioned 
c     GMRES iteration.

     
      include 'SIZE'
      include 'TOTAL'
      include 'FDMH1'
      include 'GMRES'
      common  /ctolpr/ divex
      common  /cprint/ ifprint
      logical          ifprint
      real             res  (lx1*ly1*lz1*lelv)
      real             h1   (lx1,ly1,lz1,lelv)
      real             h2   (lx1,ly1,lz1,lelv)
      real             wt   (lx1,ly1,lz1,lelv)

      common /scrcg/ d(lx1*ly1*lz1*lelv),wk(lx1*ly1*lz1*lelv)

      common /cgmres1/ y(lgmres)
      real alpha, l, temp
      integer j,m

      logical iflag
      save    iflag
      data    iflag /.false./
      real    norm_fac
      save    norm_fac

      real*8 etime1,dnekclock


      n = nx1*ny1*nz1*nelv

      etime1 = dnekclock()
      divex = 0.
      iter  = 0
      m     = lgmres

      if(.not.iflag) then
         iflag=.true.
         call uzawa_gmres_split(ml,mu,bm1,binvm1,nx1*ny1*nz1*nelv)
         norm_fac = 1./sqrt(volvm1)
      endif

      call set_fdm_prec_h1b(d,h1,h2,nelv)
      if (ifvcor) smean = -1./glsum(vmult,n)

      call chktcg1(tolps,res,h1,h2,pmask,vmult,1,1)
      if (param(21).gt.0.and.tolps.gt.abs(param(21))) 
     $   tolps = abs(param(21))
      if (istep.eq.0) tolps = 1.e-4
      tolpss = tolps
c
      iconv = 0
      call rzero(x,n)
      do while(iconv.eq.0.and.iter.lt.500)
         if(iter.eq.0) then               !      -1
            call col3(r,ml,res,n)         ! r = L  res
c           call copy(r,res,n)
         else
            !update residual
            call copy  (r,res,n)                  ! r = res
            call ax    (w,x,h1,h2,n)              ! w = A x
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
c           call hsmg_solve(z(1,j),w)             ! z  = M   w
                                                  !  j        


c . . . . . Overlapping Schwarz + coarse-grid . . . . . . .
            kfldfdm = ndim+1
            call fdm_h1
     $           (z(1,j),w,d,pmask,vmult,nelv,ktype(1,1,kfldfdm),wk)
            call crs_solve_h1 (wk,w)  ! Currently, crs grd only for P
            call add2         (z(1,j),wk,n)
            if (ifvcor) then
               rmean = smean*glsc2(z(1,j),vmult,n)
               call cadd(z(1,j),rmean,n)
            endif
c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 

     
            call ax  (w,z(1,j),h1,h2,n)           ! w = A z
                                                  !        j
     
                                                  !      -1
            call col2(w,ml,n)                     ! w = L   w
            !modified Gram-Schmidt
            do i=1,j
               h(i,j)=glsc3(w,v(1,i),wt,n)        ! h    = (w,v )
                                                  !  i,j       i
                                               
               call add2s2(w,v(1,i),-h(i,j),n)    ! w = w - h    v
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
         enddo                               !          i  i
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

      return
      end
c-----------------------------------------------------------------------
