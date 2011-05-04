c Some relevant parameters
c
c param(41):
c     0 - use additive SEMG
c     1 - use hybrid SEMG (not yet working... but coming soon!)
c
c param(42):   navier0.f, fasts.f
c     0 - use GMRES for iterative solver, use non-symmetric weighting
c     1 - use PCG for iterative solver, do not use weighting
c
c param(43):   uzawa_gmres.f, navier6.f
c     0 - use additive multilevel scheme (requires param(42).eq.0)
c     1 - use original 2 level scheme
c
c param(44):   fast3d.f, navier6.f
c     0 - base top-level additive Schwarz on restrictions of E
c     1 - base top-level additive Schwarz on restrictions of A
c
c----------------------------------------------------------------------
      subroutine hsmg_setup()
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'SEMHAT'
      include 'TSTEP'

      integer nf,nc,nr
      integer nx,ny,nz

      mg_fld = 1
      if (ifield.gt.1) mg_fld = 2
      if (ifield.eq.1) call hsmg_index_0 ! initialize index sets

c     set up the nx values for each level of multigrid
      call hsmg_setup_mg_nx

c     set up the spectral element hat matrices
      call hsmg_setup_semhat

      call hsmg_setup_intp

c     set up the direct stiffness summation handles
      call hsmg_setup_dssum

c     set up restriction weight matrices
c     and the boundary condition masks
      call hsmg_setup_wtmask

c     set up the fast diagonalization method
      call hsmg_setup_fdm
      call hsmg_setup_schwarz_wt(.false.)

c     set up the solver
      call hsmg_setup_solve
c     call hsmg_setup_dbg
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_semhat
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      include 'SEMHAT'
      integer n,l
c     generate the SEM hat matrices for each level
c     top level
      n = mg_nx(mg_lmax)
      call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
      call copy(mg_zh(1,mg_lmax),zgl,n-1) !top level based on gl points
      mg_nh(mg_lmax)=n-1
      mg_nhz(mg_lmax)=n-1
      if(.not.if3d) mg_nhz(mg_lmax)=1
c     lower levels
      do l=1,mg_lmax-1
         n = mg_nx(l)
         if(n.gt.1) then
            call semhat(ah,bh,ch,dh,zh,dph,jph,bgl,zgl,dgl,jgl,n,wh)
            call copy(mg_ah(1,l),ah,(n+1)*(n+1))
            call copy(mg_bh(1,l),bh,n+1)
            call copy(mg_dh(1,l),dh,(n+1)*(n+1))
            call transpose(mg_dht(1,l),n+1,dh,n+1)
            call copy(mg_zh(1,l),zh,n+1)
         else
            mg_zh(1,l) = -1.
            mg_zh(2,l) =  1.
         endif
         mg_nh(l)=n+1
         mg_nhz(l)=mg_nz(l)+1
      enddo
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_intp
      include 'SIZE'
      include 'HSMG'
      include 'SEMHAT'
      integer l,nf,nc
c     set up interpolation matrices
      do l=1,mg_lmax-1
         nf=mg_nh(l+1)
         nc=mg_nh(l)
         call hsmg_setup_intpm(
     $           mg_jh(1,l),mg_zh(1,l+1),mg_zh(1,l),nf,nc)
         call transpose(mg_jht(1,l),nc,mg_jh(1,l),nf)
c        call outmat(mg_jh(1,l),nf,nc,'  jh  ',l)
      enddo
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_intpm(jh,zf,zc,nf,nc)
      integer nf,nc
      real jh(nf,nc),zf(1),zc(1)
      include 'SIZE'
      real w(2*lx1+2)
      do i=1,nf
         call fd_weights_full(zf(i),zc,nc-1,1,w)
         do j=1,nc
            jh(i,j)=w(j)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_dssum
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      common /c_is1/ glo_num(lx1*ly1*lz1*lelv)
      common /ivrtx/ vertex ((2**ldim)*lelt)

      integer*8 glo_num
      integer vertex
      integer nx,ny,nz
      integer l
      
c     set up direct stiffness summation for each level
      ncrnr = 2**ndim
      call get_vert

c++   write(6,*) mg_fld,' mgfld in hsmg_setup_dssum'

      do l=1,mg_lmax-1
         nx=mg_nh(l)
         ny=mg_nh(l)
         nz=mg_nhz(l)
         call setupds(mg_gsh_handle(l,mg_fld),nx,ny,nz
     $                ,nelv,nelgv,vertex,glo_num)
         nx=nx+2
         ny=ny+2
         nz=nz+2
         if(.not.if3d) nz=1
         call setupds(mg_gsh_schwarz_handle(l,mg_fld),nx,ny,nz
     $                ,nelv,nelgv,vertex,glo_num)
      enddo
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_wtmask
      include 'SIZE'
      include 'HSMG'
      integer i,l
      i = mg_mask_index(mg_lmax,mg_fld-1)
      do l=1,mg_lmax-1
         mg_rstr_wt_index(l,mg_fld)=i
         mg_mask_index   (l,mg_fld)=i
         i=i+mg_nh(l)*mg_nhz(l)*2*ndim*nelv
         if(i .gt. lmgs*lmg_rwt*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'parameter lmg_rwt too small',i,itmp,lmg_rwt
            call exitt
         endif
         call hsmg_setup_rstr_wt(
     $           mg_rstr_wt(mg_rstr_wt_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
         call hsmg_setup_mask(
     $           mg_mask(mg_mask_index(l,mg_fld))
     $          ,mg_nh(l),mg_nh(l),mg_nhz(l),l,mg_work)
      enddo
      mg_mask_index(l,mg_fld)=i
      end
c----------------------------------------------------------------------
      subroutine hsmg_intp(uf,uc,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      call hsmg_tnsr(uf,mg_nh(l+1),uc,mg_nh(l),mg_jh(1,l),mg_jht(1,l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_rstr(uc,uf,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      if(l.ne.mg_lmax-1)
     $   call hsmg_do_wt(uf,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call hsmg_tnsr(uc,mg_nh(l),uf,mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      call hsmg_dssum(uc,l)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_rstr_no_dssum(uc,uf,l) ! l is coarse level
      real uf(1),uc(1)
      integer l
      include 'SIZE'
      include 'HSMG'
      if(l.ne.mg_lmax-1)
     $   call hsmg_do_wt(uf,mg_rstr_wt(mg_rstr_wt_index(l+1,mg_fld))
     $                     ,mg_nh(l+1),mg_nh(l+1),mg_nhz(l+1))
      call hsmg_tnsr(uc,mg_nh(l),uf,mg_nh(l+1),mg_jht(1,l),mg_jh(1,l))
      return
      end
c----------------------------------------------------------------------
c     computes
c     v = [A (x) A] u      or
c     v = [A (x) A (x) A] u 
      subroutine hsmg_tnsr(v,nv,u,nu,A,At)
      integer nv,nu
      real v(1),u(1),A(1),At(1)
      include 'SIZE'
      include 'INPUT'
      if (.not. if3d) then
         call hsmg_tnsr2d(v,nv,u,nu,A,At)
      else
         call hsmg_tnsr3d(v,nv,u,nu,A,At,At)
      endif
      return
      end
c----------------------------------------------------------------------
c     computes
c              T
c     v = A u B
      subroutine hsmg_tnsr2d(v,nv,u,nu,A,Bt)
      integer nv,nu
      real v(nv*nv,nelv),u(nu*nu,nelv),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work(lx1*lx1)
      integer ie
      do ie=1,nelv
         call mxm(A,nv,u(1,ie),nu,work,nu)
         call mxm(work,nv,Bt,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine hsmg_tnsr3d(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv,nelv),u(nu*nu*nu,nelv),A(1),Bt(1),Ct(1)
      include 'SIZE'
      common /hsmgw/ work(0:lx1*ly1*lz1-1),work2(0:lx1*ly1*lz1-1)
      integer ie, i
      do ie=1,nelv
         call mxm(A,nv,u(1,ie),nu,work,nu*nu)
         do i=0,nu-1
            call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
         enddo
         call mxm(work2,nv*nv,Ct,nu,v(1,ie),nv)
      enddo
      return
      end
c----------------------------------------------------------------------
c     computes
c              T
c     v = A u B
      subroutine hsmg_tnsr2d_el(v,nv,u,nu,A,Bt)
      integer nv,nu
      real v(nv*nv),u(nu*nu),A(1),Bt(1)
      include 'SIZE'
      common /hsmgw/ work(lx1*lx1)
c
      call mxm(A,nv,u,nu,work,nu)
      call mxm(work,nv,Bt,nu,v,nv)
c
      return
      end
c----------------------------------------------------------------------
c     computes
c              
c     v = [C (x) B (x) A] u
      subroutine hsmg_tnsr3d_el(v,nv,u,nu,A,Bt,Ct)
      integer nv,nu
      real v(nv*nv*nv),u(nu*nu*nu),A(1),Bt(1),Ct(1)
      include 'SIZE'
      common /hsmgw/ work(0:lx1*ly1*lz1-1),work2(0:lx1*ly1*lz1-1)
      integer i
c
      call mxm(A,nv,u,nu,work,nu*nu)
      do i=0,nu-1
         call mxm(work(nv*nu*i),nv,Bt,nu,work2(nv*nv*i),nv)
      enddo
      call mxm(work2,nv*nv,Ct,nu,v,nv)
c
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_dssum(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'

      if (ifsync) call gsync()
#ifndef NOTIMER
      etime1=dnekclock()
#endif
      call gs_op(mg_gsh_handle(l,mg_fld),u,1,1,0)
#ifndef NOTIMER
      tdadd =tdadd + dnekclock()-etime1
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_dsprod(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'

      if (ifsync) call gsync()

      call gs_op(mg_gsh_handle(l,mg_fld),u,1,2,0)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_dssum(u,l)
      include 'SIZE'
      include 'HSMG'
      include 'CTIMER'

      if (ifsync) call gsync()
#ifndef NOTIMER
      etime1=dnekclock()
#endif
      call gs_op(mg_gsh_schwarz_handle(l,mg_fld),u,1,1,0)
#ifndef NOTIMER
      tdadd =tdadd + dnekclock()-etime1
#endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer l1,l2,nx,ny,nz
      real arr1(nx,ny,nz,nelv),arr2(nx,ny,nz,nelv)
      real f1,f2
      
      integer i,j,k,ie,i0,i1
      i0=2
      i1=nx-1
      
      if(.not.if3d) then
         do ie=1,nelv
            do j=i0,i1
               arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie)
     $                             +f2*arr2(l2+1 ,j,1,ie)
               arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie)
     $                             +f2*arr2(nx-l2,j,1,ie)
            enddo
            do i=i0,i1
               arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie)
     $                             +f2*arr2(i,l2+1 ,1,ie)
               arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie)
     $                             +f2*arr2(i,nx-l2,1,ie)
            enddo
         enddo
      else
         do ie=1,nelv
            do k=i0,i1
            do j=i0,i1
               arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie)
     $                             +f2*arr2(l2+1 ,j,k,ie)
               arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie)
     $                             +f2*arr2(nx-l2,j,k,ie)
            enddo
            enddo
            do k=i0,i1
            do i=i0,i1
               arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie)
     $                             +f2*arr2(i,l2+1 ,k,ie)
               arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie)
     $                             +f2*arr2(i,nx-l2,k,ie)
            enddo
            enddo
            do j=i0,i1
            do i=i0,i1
               arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie)
     $                             +f2*arr2(i,j,l2+1 ,ie)
               arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie)
     $                             +f2*arr2(i,j,nx-l2,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz(e,r,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      real e(1),r(1)
      integer l
      integer enx,eny,enz
      integer i

      real zero,one,onem
      zero =  0
      one  =  1
      onem = -1

c     apply mask (zeros Dirichlet nodes)
      !!!!! uncommenting
      call hsmg_do_wt(r,mg_mask(mg_mask_index(l,mg_fld)),
     $                mg_nh(l),mg_nh(l),mg_nhz(l))
      
c     go to extended size array (room for overlap)      
      if(.not.if3d) then
         call hsmg_schwarz_toext2d(mg_work,r,mg_nh(l))
      else
         call hsmg_schwarz_toext3d(mg_work,r,mg_nh(l))
      endif
      
      enx=mg_nh(l)+2
      eny=mg_nh(l)+2
      enz=mg_nh(l)+2
      if(.not.if3d) enz=1
      i = enx*eny*enz*nelv+1
c     exchange interior nodes
c
      call hsmg_extrude(mg_work,0,zero,mg_work,2,one,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work,l)
      call hsmg_extrude(mg_work,0,one ,mg_work,2,onem,enx,eny,enz)
c     do the local solves
      call hsmg_fdm(mg_work(i),mg_work,l)
c     sum overlap region (border excluded)
      call hsmg_extrude(mg_work,0,zero,mg_work(i),0,one ,enx,eny,enz)
      call hsmg_schwarz_dssum(mg_work(i),l)
      call hsmg_extrude(mg_work(i),0,one ,mg_work,0,onem,enx,eny,enz)
      call hsmg_extrude(mg_work(i),2,one,mg_work(i),0,one,enx,eny,enz)
c     go back to regular size array
      if(.not.if3d) then
         call hsmg_schwarz_toreg2d(e,mg_work(i),mg_nh(l))
      else
         call hsmg_schwarz_toreg3d(e,mg_work(i),mg_nh(l))
      endif
c     sum border nodes
      call hsmg_dssum(e,l)
c     apply mask (zeros Dirichlet nodes)
      !!!!!! changing r to e
      call hsmg_do_wt(e,mg_mask(mg_mask_index(l,mg_fld)),
     $                mg_nh(l),mg_nh(l),mg_nhz(l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toext2d(a,b,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,nelv),b(n,n,nelv)
      
      integer i,j,ie
c      call rzero(a,(n+2)*(n+2)*nelv)
      do ie=1,nelv
         do i=0,n+1
            a(i,0,ie)=0.
         enddo
         do j=1,n
            a(0  ,j,ie)=0.
            do i=1,n
               a(i,j,ie)=b(i,j,ie)
            enddo
            a(n+1,j,ie)=0.
         enddo
         do i=0,n+1
            a(i,n+1,ie)=0.
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toext3d(a,b,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
      
      integer i,j,k,ie
      call rzero(a,(n+2)*(n+2)*(n+2)*nelv)
      do ie=1,nelv
      do k=1,n
      do j=1,n
      do i=1,n
         a(i,j,k,ie)=b(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toreg2d(b,a,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,nelv),b(n,n,nelv)
      
      integer i,j,ie
      do ie=1,nelv
      do j=1,n
      do i=1,n
         b(i,j,ie)=a(i,j,ie)
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_toreg3d(b,a,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
      
      integer i,j,k,ie
      do ie=1,nelv
      do k=1,n
      do j=1,n
      do i=1,n
         b(i,j,k,ie)=a(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fdm()
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      
      integer l,i,j,nl
      i = mg_fast_s_index(mg_lmax,mg_fld-1)
      j = mg_fast_d_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax-1
         mg_fast_s_index(l,mg_fld)=i
         nl = mg_nh(l)+2
         i=i+nl*nl*2*ndim*nelv
         if(i .gt. lmg_fasts*2*ldim*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fasts too small',i,itmp,lmg_fasts,l
            call exitt
         endif
         mg_fast_d_index(l,mg_fld)=j
         j=j+(nl**ndim)*nelv
         if(j .gt. lmg_fastd*lelv) then
            itmp = i/(2*ldim*lelv)
            write(6,*) 'lmg_fastd too small',i,itmp,lmg_fastd,l
            call exitt
         endif
         call hsmg_setup_fast(
     $             mg_fast_s(mg_fast_s_index(l,mg_fld))
     $            ,mg_fast_d(mg_fast_d_index(l,mg_fld))
     $            ,mg_nh(l)+2,mg_ah(1,l),mg_bh(1,l),mg_nx(l))
      enddo
      mg_fast_s_index(l,mg_fld)=i
      mg_fast_d_index(l,mg_fld)=j
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast(s,d,nl,ah,bh,n)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      real s(nl*nl,2,ndim,nelv)
      real d(nl**ndim,nelv)
      real ah(1),bh(1)
      common /ctmpf/  lr(2*lx1),ls(2*lx1),lt(2*lx1)
     $              , llr(lelt),lls(lelt),llt(lelt)
     $              , lmr(lelt),lms(lelt),lmt(lelt)
     $              , lrr(lelt),lrs(lelt),lrt(lelt)
      real lr ,ls ,lt
      real llr,lls,llt
      real lmr,lms,lmt
      real lrr,lrs,lrt
      
      integer i,j,k
      integer ie,il,nr,ns,nt
      integer lbr,rbr,lbs,rbs,lbt,rbt
      real eps,diag
      
      ierr = 0
      do ie=1,nelv
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,2,ierr)
         nr=nl
         ns=nl
         nt=nl
         call hsmg_setup_fast1d(s(1,1,1,ie),lr,nr,lbr,rbr
     $            ,llr(ie),lmr(ie),lrr(ie),ah,bh,n,ie)
         call hsmg_setup_fast1d(s(1,1,2,ie),ls,ns,lbs,rbs
     $            ,lls(ie),lms(ie),lrs(ie),ah,bh,n,ie)
         if(if3d) call hsmg_setup_fast1d(s(1,1,3,ie),lt,nt,lbt,rbt
     $                     ,llt(ie),lmt(ie),lrt(ie),ah,bh,n,ie)
         il=1
         if(.not.if3d) then
            eps = 1.e-5*(vlmax(lr(2),nr-2) + vlmax(ls(2),ns-2))
            do j=1,ns
            do i=1,nr
               diag = lr(i)+ls(j)
               if (diag.gt.eps) then
                  d(il,ie) = 1.0/diag
               else
c                 write(6,2) ie,'Reset Eig in hsmg setup fast:',i,j,l
c    $                         ,eps,diag,lr(i),ls(j)
    2             format(i6,1x,a21,3i5,1p4e12.4)
                  d(il,ie) = 0.0
               endif
               il=il+1
            enddo
            enddo
         else
            eps = 1.e-5 * (vlmax(lr(2),nr-2)
     $                  + vlmax(ls(2),ns-2) + vlmax(lt(2),nt-2))
            do k=1,nt
            do j=1,ns
            do i=1,nr
               diag = lr(i)+ls(j)+lt(k)
               if (diag.gt.eps) then
                  d(il,ie) = 1.0/diag
               else
c                 write(6,3) ie,'Reset Eig in hsmg setup fast:',i,j,k,l
c    $                         ,eps,diag,lr(i),ls(j),lt(k)
    3             format(i6,1x,a21,4i5,1p5e12.4)
                  d(il,ie) = 0.0
               endif
               il=il+1
            enddo
            enddo
            enddo
         endif
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('A INVALID BC FOUND in genfast$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d(s,lam,nl,lbc,rbc,ll,lm,lr,ah,bh,n,ie)
      integer nl,lbc,rbc,n
      real s(nl,nl,2),lam(nl),ll,lm,lr
      real ah(0:n,0:n),bh(0:n)
      
      include 'SIZE'
      real b(2*lx1*lx1),w(2*lx1*lx1)
      
      call hsmg_setup_fast1d_a(s,lbc,rbc,ll,lm,lr,ah,n)
      call hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
      
c     if (nid.eq.0) write(6,*) 'THIS is generalev call',nl,lbc
      call generalev(s,b,lam,nl,w)
      if(lbc.gt.0) call row_zero(s,nl,nl,1)
      if(lbc.eq.1) call row_zero(s,nl,nl,2)
      if(rbc.gt.0) call row_zero(s,nl,nl,nl)
      if(rbc.eq.1) call row_zero(s,nl,nl,nl-1)
      
      call transpose(s(1,1,2),nl,s,nl)
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_a(a,lbc,rbc,ll,lm,lr,ah,n)
      integer lbc,rbc,n
      real a(0:n+2,0:n+2),ll,lm,lr
      real ah(0:n,0:n)
      
      real fac
      integer i,j,i0,i1
      i0=0
      if(lbc.eq.1) i0=1
      i1=n
      if(rbc.eq.1) i1=n-1
      
      call rzero(a,(n+3)*(n+3))
      fac = 2.0/lm
      a(1,1)=1.0
      a(n+1,n+1)=1.0
      do j=i0,i1
         do i=i0,i1
            a(i+1,j+1)=fac*ah(i,j)
         enddo
      enddo
      if(lbc.eq.0) then
         fac = 2.0/ll
         a(0,0)=fac*ah(n-1,n-1)
         a(1,0)=fac*ah(n  ,n-1)
         a(0,1)=fac*ah(n-1,n  )
         a(1,1)=a(1,1)+fac*ah(n  ,n  )
      else
         a(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 2.0/lr
         a(n+1,n+1)=a(n+1,n+1)+fac*ah(0,0)
         a(n+2,n+1)=fac*ah(1,0)
         a(n+1,n+2)=fac*ah(0,1)
         a(n+2,n+2)=fac*ah(1,1)
      else
         a(n+2,n+2)=1.0
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_fast1d_b(b,lbc,rbc,ll,lm,lr,bh,n)
      integer lbc,rbc,n
      real b(0:n+2,0:n+2),ll,lm,lr
      real bh(0:n)
      
      real fac
      integer i,j,i0,i1
      i0=0
      if(lbc.eq.1) i0=1
      i1=n
      if(rbc.eq.1) i1=n-1
      
      call rzero(b,(n+3)*(n+3))
      fac = 0.5*lm
      b(1,1)=1.0
      b(n+1,n+1)=1.0
      do i=i0,i1
         b(i+1,i+1)=fac*bh(i)
      enddo
      if(lbc.eq.0) then
         fac = 0.5*ll
         b(0,0)=fac*bh(n-1)
         b(1,1)=b(1,1)+fac*bh(n  )
      else
         b(0,0)=1.0
      endif
      if(rbc.eq.0) then
         fac = 0.5*lr
         b(n+1,n+1)=b(n+1,n+1)+fac*bh(0)
         b(n+2,n+2)=fac*bh(1)
      else
         b(n+2,n+2)=1.0
      endif
      return
      end
c----------------------------------------------------------------------
c     clobbers r
      subroutine hsmg_fdm(e,r,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      call hsmg_do_fast(e,r,
     $      mg_fast_s(mg_fast_s_index(l,mg_fld)),
     $      mg_fast_d(mg_fast_d_index(l,mg_fld)),
     $      mg_nh(l)+2)
      return
      end
c----------------------------------------------------------------------
c     clobbers r
      subroutine hsmg_do_fast(e,r,s,d,nl)
      include 'SIZE'
      include 'INPUT'
      real e(nl**ndim,nelv)
      real r(nl**ndim,nelv)
      real s(nl*nl,2,ndim,nelv)
      real d(nl**ndim,nelv)
      
      integer ie,nn,i
      nn=nl**ndim
      if(.not.if3d) then
         do ie=1,nelv
            call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,2,1,ie),s(1,1,2,ie))
            do i=1,nn
               r(i,ie)=d(i,ie)*e(i,ie)
            enddo
            call hsmg_tnsr2d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,1,1,ie),s(1,2,2,ie))
         enddo
      else
         do ie=1,nelv
            call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,2,1,ie),s(1,1,2,ie),s(1,1,3,ie))
            do i=1,nn
               r(i,ie)=d(i,ie)*e(i,ie)
            enddo
            call hsmg_tnsr3d_el(e(1,ie),nl,r(1,ie),nl
     $                         ,s(1,1,1,ie),s(1,2,2,ie),s(1,2,3,ie))
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
c     u = wt .* u
      subroutine hsmg_do_wt(u,wt,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz
      real u(nx,ny,nz,nelv)
      real wt(nx,nz,2,ndim,nelv)
      
      integer ie
      if (.not. if3d) then
         do ie=1,nelv
            do j=1,ny
               u( 1,j,1,ie)=u( 1,j,1,ie)*wt(j,1,1,1,ie)
               u(nx,j,1,ie)=u(nx,j,1,ie)*wt(j,1,2,1,ie)
            enddo
            do i=2,nx-1
               u(i, 1,1,ie)=u(i, 1,1,ie)*wt(i,1,1,2,ie)
               u(i,ny,1,ie)=u(i,ny,1,ie)*wt(i,1,2,2,ie)
            enddo
         enddo
      else
         do ie=1,nelv
            do k=1,nz
            do j=1,ny
               u( 1,j,k,ie)=u( 1,j,k,ie)*wt(j,k,1,1,ie)
               u(nx,j,k,ie)=u(nx,j,k,ie)*wt(j,k,2,1,ie)
            enddo
            enddo
            do k=1,nz
            do i=2,nx-1
               u(i, 1,k,ie)=u(i, 1,k,ie)*wt(i,k,1,2,ie)
               u(i,ny,k,ie)=u(i,ny,k,ie)*wt(i,k,2,2,ie)
            enddo
            enddo
            do j=2,ny-1
            do i=2,nx-1
               u(i,j, 1,ie)=u(i,j, 1,ie)*wt(i,j,1,3,ie)
               u(i,j,nz,ie)=u(i,j,nz,ie)*wt(i,j,2,3,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_rstr_wt(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelv)
      real wt(nx,nz,2,ndim,nelv)
      
      integer ie
      !init border nodes to 1
      call rzero(w,nx*ny*nz*nelv)
c     print *, 'Setup rstr wt: ',nx,ny,nz,nelv
      if (.not.if3d) then
         do ie=1,nelv
            do i=1,nx
               w(i,1,1,ie)=1.0
               w(i,ny,1,ie)=1.0
            enddo
            do j=1,ny
               w(1,j,1,ie)=1.0
               w(nx,j,1,ie)=1.0
            enddo
         enddo
      else
         do ie=1,nelv
            do j=1,ny
            do i=1,nx
               w(i,j,1,ie)=1.0
               w(i,j,nz,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=1.0
               w(i,ny,k,ie)=1.0
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=1.0
               w(nx,j,k,ie)=1.0
            enddo
            enddo
         enddo
      endif
      call hsmg_dssum(w,l)
      !invert the count w to get the weight wt
      if (.not. if3d) then
         do ie=1,nelv
            do j=1,ny
               wt(j,1,1,1,ie)=1.0/w(1,j,1,ie)
               wt(j,1,2,1,ie)=1.0/w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=1.0/w(i,1,1,ie)
               wt(i,1,2,2,ie)=1.0/w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelv
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=1.0/w(1,j,k,ie)
               wt(j,k,2,1,ie)=1.0/w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=1.0/w(i,1,k,ie)
               wt(i,k,2,2,ie)=1.0/w(i,ny,k,ie)
            enddo
            enddo
            do j=1,ny
            do i=1,nx
               wt(i,j,1,3,ie)=1.0/w(i,j,1,ie)
               wt(i,j,2,3,ie)=1.0/w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_mask(wt,nx,ny,nz,l,w)
      include 'SIZE'
      include 'INPUT'
      integer nx,ny,nz,l
      real w(nx,ny,nz,nelv)
      real wt(nx,nz,2,ndim,nelv)
      
      integer ie
      integer lbr,rbr,lbs,rbs,lbt,rbt
c     init everything to 1
      do ie=1,nelv
      do k=1,nz
      do j=1,ny
      do i=1,nx
         w(i,j,k,ie)=1.0
      enddo
      enddo
      enddo
      enddo

c     set dirichlet nodes to zero
      ierr = 0
      do ie=1,nelv
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,2,ierr)
         if(lbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(1,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbr.eq.1) then
            do k=1,nz
            do j=1,ny
               w(nx,j,k,ie)=0.0
            enddo
            enddo
         endif
         if(lbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,1,k,ie)=0.0
            enddo
            enddo
         endif
         if(rbs.eq.1) then
            do k=1,nz
            do i=1,nx
               w(i,ny,k,ie)=0.0
            enddo
            enddo
         endif
         if(if3d) then
            if(lbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,1,ie)=0.0
               enddo
               enddo
            endif
            if(rbt.eq.1) then
               do j=1,ny
               do i=1,nx
                  w(i,j,nz,ie)=0.0
               enddo
               enddo
            endif
         endif
      enddo
c     do direct stiffness multiply

      call hsmg_dsprod(w,l)


c     store weight
      if (.not. if3d) then
         do ie=1,nelv
            do j=1,ny
               wt(j,1,1,1,ie)=w(1,j,1,ie)
               wt(j,1,2,1,ie)=w(nx,j,1,ie)
            enddo
            do i=1,nx
               wt(i,1,1,2,ie)=w(i,1,1,ie)
               wt(i,1,2,2,ie)=w(i,ny,1,ie)
            enddo
         enddo
      else
         do ie=1,nelv
            do k=1,nz
            do j=1,ny
               wt(j,k,1,1,ie)=w(1,j,k,ie)
               wt(j,k,2,1,ie)=w(nx,j,k,ie)
            enddo
            enddo
            do k=1,nz
            do i=1,nx
               wt(i,k,1,2,ie)=w(i,1,k,ie)
               wt(i,k,2,2,ie)=w(i,ny,k,ie)
            enddo
            enddo
            do k=1,nz
            do j=1,ny
               wt(j,k,1,3,ie)=w(i,j,1,ie)
               wt(j,k,2,3,ie)=w(i,j,nz,ie)
            enddo
            enddo
         enddo
      endif

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('B INVALID BC FOUND in genfast$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_schwarz_wt(ifsqrt)
      logical ifsqrt
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      
      integer l,i,nl,nlz

      i = mg_schwarz_wt_index(mg_lmax,mg_fld-1)
      do l=2,mg_lmax-1
         mg_schwarz_wt_index(l,mg_fld)=i
         nl = mg_nh(l)
         nlz = mg_nh(l)
         if(.not.if3d) nlz=1
         i=i+nl*nlz*4*ndim*nelv
         if(i .gt. lmg_swt*4*ldim*lelv) then
            itmp = i/(4*ldim*lelv)
            write(6,*) 'lmg_swt too small',i,itmp,lmg_swt,l
            call exitt
         endif
         if(.not.if3d) call hsmg_setup_schwarz_wt2d(
     $       mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld))
     $      ,nl,mg_worke,ifsqrt)
c        if(if3d) write(6,*) mg_schwarz_wt_index(l,mg_fld),l,'SCHWARZ'
         if(if3d) call hsmg_setup_schwarz_wt3d(
     $       mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld))
     $      ,nl,mg_worke,ifsqrt)
      enddo
      mg_schwarz_wt_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_schwarz_wt2d(wt,n,work,ifsqrt)
      logical ifsqrt
      include 'SIZE'
      integer n
      real wt(n,4,2,nelv)
      real work(n,n)
      
      integer ie,i,j
      integer lbr,rbr,lbs,rbs,lbt,rbt
      ierr = 0
      do ie=1,nelv
         call rzero(work,n*n)
         do j=1,n
            work(1,j)=1.0
            work(2,j)=1.0
            work(n-1,j)=1.0
            work(n,j)=1.0
         enddo
         do i=1,n
            work(i,1)=1.0
            work(i,2)=1.0
            work(i,n-1)=1.0
            work(i,n)=1.0
         enddo
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,2,ierr)
         if(lbr.eq.0) then
            do j=1,n
               work(1,j)=work(1,j)+1.0
               work(2,j)=work(2,j)+1.0
            enddo
         endif
         if(rbr.eq.0) then
            do j=1,n
               work(n-1,j)=work(n-1,j)+1.0
               work(n,j)=work(n,j)+1.0
            enddo
         endif
         if(lbs.eq.0) then
            do i=1,n
               work(i,1)=work(i,1)+1.0
               work(i,2)=work(i,2)+1.0
            enddo
         endif
         if(rbs.eq.0) then
            do i=1,n
               work(i,n-1)=work(i,n-1)+1.0
               work(i,n)=work(i,n)+1.0
            enddo
         endif
         do j=1,n
            wt(j,1,1,ie)=1.0/work(1,j)
            wt(j,2,1,ie)=1.0/work(2,j)
            wt(j,3,1,ie)=1.0/work(n-1,j)
            wt(j,4,1,ie)=1.0/work(n,j)
         enddo
         do i=1,n
            wt(i,1,2,ie)=1.0/work(i,1)
            wt(i,2,2,ie)=1.0/work(i,2)
            wt(i,3,2,ie)=1.0/work(i,n-1)
            wt(i,4,2,ie)=1.0/work(i,n)
         enddo
         if(ifsqrt) then
            do ii=1,2
            do j=1,4
            do i=1,n
               wt(i,j,ii,ie)=sqrt(wt(i,j,ii,ie))
            enddo
            enddo
            enddo
         endif
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('C INVALID BC FOUND in genfast$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_schwarz_wt3d(wt,n,work,ifsqrt)
      logical ifsqrt
      include 'SIZE'
      integer n
      real wt(n,n,4,3,nelv)
      real work(n,n,n)
      
      integer ie,i,j,k
      integer lbr,rbr,lbs,rbs,lbt,rbt

      ierr = 0
      do ie=1,nelv
         do k=1,n
         do j=1,n
            work(1,j,k)=1.0
            work(2,j,k)=1.0
            work(n-1,j,k)=1.0
            work(n,j,k)=1.0
         enddo
         enddo
         do k=1,n
         do i=1,n
            work(i,1,k)=1.0
            work(i,2,k)=1.0
            work(i,n-1,k)=1.0
            work(i,n,k)=1.0
         enddo
         enddo
         do j=1,n
         do i=1,n
            work(i,j,1)=1.0
            work(i,j,2)=1.0
            work(i,j,n-1)=1.0
            work(i,j,n)=1.0
         enddo
         enddo
         call get_fast_bc(lbr,rbr,lbs,rbs,lbt,rbt,ie,2,ierr)
         if(lbr.eq.0) then
            do k=1,n
            do j=1,n
               work(1,j,k)=work(1,j,k)+1.0
               work(2,j,k)=work(2,j,k)+1.0
            enddo
            enddo
         endif
         if(rbr.eq.0) then
            do k=1,n
            do j=1,n
               work(n-1,j,k)=work(n-1,j,k)+1.0
               work(n,j,k)=work(n,j,k)+1.0
            enddo
            enddo
         endif
         if(lbs.eq.0) then
            do k=1,n
            do i=1,n
               work(i,1,k)=work(i,1,k)+1.0
               work(i,2,k)=work(i,2,k)+1.0
            enddo
            enddo
         endif
         if(rbs.eq.0) then
            do k=1,n
            do i=1,n
               work(i,n-1,k)=work(i,n-1,k)+1.0
               work(i,n,k)=work(i,n,k)+1.0
            enddo
            enddo
         endif
         if(lbt.eq.0) then
            do j=1,n
            do i=1,n
               work(i,j,1)=work(i,j,1)+1.0
               work(i,j,2)=work(i,j,2)+1.0
            enddo
            enddo
         endif
         if(rbt.eq.0) then
            do j=1,n
            do i=1,n
               work(i,j,n-1)=work(i,j,n-1)+1.0
               work(i,j,n)=work(i,j,n)+1.0
            enddo
            enddo
         endif
         do k=1,n
         do j=1,n
            wt(j,k,1,1,ie)=1.0/work(1,j,k)
            wt(j,k,2,1,ie)=1.0/work(2,j,k)
            wt(j,k,3,1,ie)=1.0/work(n-1,j,k)
            wt(j,k,4,1,ie)=1.0/work(n,j,k)
         enddo
         enddo
         do k=1,n
         do i=1,n
            wt(i,k,1,2,ie)=1.0/work(i,1,k)
            wt(i,k,2,2,ie)=1.0/work(i,2,k)
            wt(i,k,3,2,ie)=1.0/work(i,n-1,k)
            wt(i,k,4,2,ie)=1.0/work(i,n,k)
         enddo
         enddo
         do j=1,n
         do i=1,n
            wt(i,j,1,3,ie)=1.0/work(i,j,1)
            wt(i,j,2,3,ie)=1.0/work(i,j,2)
            wt(i,j,3,3,ie)=1.0/work(i,j,n-1)
            wt(i,j,4,3,ie)=1.0/work(i,j,n)
         enddo
         enddo
         if(ifsqrt) then
            do ii=1,3
            do k=1,4
            do j=1,4
            do i=1,n
               wt(i,j,k,ii,ie)=sqrt(wt(i,j,k,ii,ie))
            enddo
            enddo
            enddo
            enddo
         endif
      enddo

      ierrmx = iglmax(ierr,1)
      if (ierrmx.gt.0) then
         if (ierr.gt.0) write(6,*) nid,ierr,' BC FAIL'
         call exitti('D INVALID BC FOUND in genfast$',ierrmx)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt(e,l)
      include 'SIZE'
      include 'INPUT'
      include 'HSMG'
      
      if(.not.if3d) call hsmg_schwarz_wt2d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      if(if3d) call hsmg_schwarz_wt3d(
     $    e,mg_schwarz_wt(mg_schwarz_wt_index(l,mg_fld)),mg_nh(l))
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt2d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,nelv)
      real wt(n,4,2,nelv)
      
      integer ie,i,j
      do ie=1,nelv
         do j=1,n
            e(1  ,j,ie)=e(1  ,j,ie)*wt(j,1,1,ie)
            e(2  ,j,ie)=e(2  ,j,ie)*wt(j,2,1,ie)
            e(n-1,j,ie)=e(n-1,j,ie)*wt(j,3,1,ie)
            e(n  ,j,ie)=e(n  ,j,ie)*wt(j,4,1,ie)
         enddo
         do i=3,n-2
            e(i,1  ,ie)=e(i,1  ,ie)*wt(i,1,2,ie)
            e(i,2  ,ie)=e(i,2  ,ie)*wt(i,2,2,ie)
            e(i,n-1,ie)=e(i,n-1,ie)*wt(i,3,2,ie)
            e(i,n  ,ie)=e(i,n  ,ie)*wt(i,4,2,ie)
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_schwarz_wt3d(e,wt,n)
      include 'SIZE'
      integer n
      real e(n,n,n,nelv)
      real wt(n,n,4,3,nelv)
      
      integer ie,i,j,k
      do ie=1,nelv
         do k=1,n
         do j=1,n
            e(1  ,j,k,ie)=e(1  ,j,k,ie)*wt(j,k,1,1,ie)
            e(2  ,j,k,ie)=e(2  ,j,k,ie)*wt(j,k,2,1,ie)
            e(n-1,j,k,ie)=e(n-1,j,k,ie)*wt(j,k,3,1,ie)
            e(n  ,j,k,ie)=e(n  ,j,k,ie)*wt(j,k,4,1,ie)
         enddo
         enddo
         do k=1,n
         do i=3,n-2
            e(i,1  ,k,ie)=e(i,1  ,k,ie)*wt(i,k,1,2,ie)
            e(i,2  ,k,ie)=e(i,2  ,k,ie)*wt(i,k,2,2,ie)
            e(i,n-1,k,ie)=e(i,n-1,k,ie)*wt(i,k,3,2,ie)
            e(i,n  ,k,ie)=e(i,n  ,k,ie)*wt(i,k,4,2,ie)
         enddo
         enddo
         do j=3,n-2
         do i=3,n-2
            e(i,j,1  ,ie)=e(i,j,1  ,ie)*wt(i,j,1,3,ie)
            e(i,j,2  ,ie)=e(i,j,2  ,ie)*wt(i,j,2,3,ie)
            e(i,j,n-1,ie)=e(i,j,n-1,ie)*wt(i,j,3,3,ie)
            e(i,j,n  ,ie)=e(i,j,n  ,ie)*wt(i,j,4,3,ie)
         enddo
         enddo
      enddo
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_coarse_solve(e,r)
      include 'SIZE'
      include 'DOMAIN'
      include 'ESOLV'
      include 'GEOM'
      include 'SOLN'
      include 'PARALLEL'
      include 'HSMG'
      include 'CTIMER'
      include 'INPUT'
      include 'TSTEP'
      real e(1),r(1)
c
      integer n_crs_tot
      save    n_crs_tot
      data    n_crs_tot /0/
c
      if (icalld.eq.0) then ! timer info
         ncrsl=0
         tcrsl=0.0
      endif
      icalld = 1

      if (ifsync) call gsync()

      ncrsl  = ncrsl  + 1
#ifndef NOTIMER
      etime1=dnekclock()
#endif

      call crs_solve(xxth(ifield),e,r)

#ifndef NOTIMER
      tcrsl=tcrsl+dnekclock()-etime1
#endif

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_solve
      include 'SIZE'
      include 'HSMG'
      
      integer l,i,nl,nlz
      i = mg_solve_index(mg_lmax+1,mg_fld-1)
      do l=1,mg_lmax
         mg_solve_index(l,mg_fld)=i
         i=i+mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
         if(i .gt. lmg_solve*lelv) then
            itmp = i/lelv
            write(6,*) 'lmg_solve too small',i,itmp,lmg_solve,l
            call exitt
         endif
      enddo
      mg_solve_index(l,mg_fld)=i

      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_solve(e,r)
      real e(1),r(1)
      include 'SIZE'
      include 'HSMG'
      include 'GEOM'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'
      
      common /quick/ ecrs  (2)  ! quick work array
     $             , ecrs2 (2)  ! quick work array
c     common /quick/ ecrs  (lx2*ly2*lz2*lelv)  ! quick work array
c    $             , ecrs2 (lx2*ly2*lz2*lelv)  ! quick work array

      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv),
     $               h2    (lx1,ly1,lz1,lelv)
      
      integer ilstep,iter
      save    ilstep,iter
      data    ilstep,iter /0,0/

      real    rhoavg,copt(2),copw(2)
      save    rhoavg,copt1,copt2
      data    rhoavg,copt1,copt2 /3*1./  ! Default copt = 1 for additive

      integer l,nt
      integer*8 ntotg,nxyz2

      logical if_hybrid

      mg_fld = 1
      if (ifield.gt.1) mg_fld = 2

      if (istep.ne.ilstep) then
         ilstep = istep
         ntot1  = nx1*ny1*nz1*nelv
         rhoavg = glsc2(vtrans,bm1,ntot1)/volvm1
      endif

      n = nx2*ny2*nz2*nelv
c     call copy(e,r,n)
c     return
 
      if (icalld.eq.0) then

         tddsl=0.0
         nddsl=0

         icalld = 1
         taaaa = 0
         tbbbb = 0
         tcccc = 0
         tdddd = 0
         teeee = 0
      endif

      nddsl  = nddsl  + 1
#ifndef NOTIMER
      etime1 = dnekclock()
#endif
      
c     n = nx2*ny2*nz2*nelv
c     rmax = glmax(r,n)
c     if (nid.eq.0) write(6,*) istep,n,rmax,' rmax1'
       
      iter = iter + 1      

      l = mg_lmax
      nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
      ! e := W M        r
      !         Schwarz
      time_0 = dnekclock()
      call local_solves_fdm(e,r)

      time_1 = dnekclock()

c     if (param(41).eq.1) if_hybrid = .true.
      if_hybrid = .false.

      if (if_hybrid) then
         ! w := E e
         rbd1dt = rhoavg*bd(1)/dt ! Assumes constant density!!!
         call cdabdtp(mg_work2,e,h1,h2,h2inv,1)
         call cmult  (mg_work2,rbd1dt,nt)
         time_2 = dnekclock()
         if (istep.eq.1) then
            copt(1)  = vlsc2(r       ,mg_work2,nt)
            copt(2)  = vlsc2(mg_work2,mg_work2,nt)
            call gop(copt,copw,'+  ', 2)
            copt(1)  = copt(1)/copt(2)
            avg2     = 1./iter
            avg1     = 1.-avg2
            copt1    = avg1*copt1 + avg2*copt(1)
            if(nid.eq.0)write(6,1)istep,iter,rbd1dt,copt(1),copt1,'cpt1'
    1       format(2i6,1p3e14.5,2x,a4)
         endif
         ! w := r - w
         do i = 1,nt
            mg_work2(i) = r(i) - copt1*mg_work2(i)
            e       (i) = copt1*e(i)
            ecrs2   (i) = mg_work2(i)
         enddo

      else   ! Additive
         ! w := r - w
         do i = 1,nt
            mg_work2(i) = r(i)
         enddo
         time_2 = dnekclock()
      endif
 
      do l = mg_lmax-1,2,-1

c        rmax = glmax(mg_work2,nt)
c        if (nid.eq.0) write(6,*) l,nt,rmax,' rmax2'

         nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
         !          T
         ! r   :=  J w
         !  l         
         call hsmg_rstr(mg_solve_r(mg_solve_index(l,mg_fld)),mg_work2,l)

         ! w  := r
         !        l
         call copy(mg_work2,mg_solve_r(mg_solve_index(l,mg_fld)),nt)
         ! e  := M        w
         !  l     Schwarz  
         call hsmg_schwarz(
     $          mg_solve_e(mg_solve_index(l,mg_fld)),mg_work2,l)

         ! e  := W e
         !  l       l
         call hsmg_schwarz_wt(mg_solve_e(mg_solve_index(l,mg_fld)),l)

c        call exitti('quit in mg$',l)

         ! w  := r  - w
         !        l
         do i = 0,nt-1
            mg_work2(i+1) = mg_solve_r(mg_solve_index(l,mg_fld)+i)
     $         !-alpha*mg_work2(i+1)
         enddo
      enddo

      call hsmg_rstr_no_dssum(
     $   mg_solve_r(mg_solve_index(1,mg_fld)),mg_work2,1)

      nzw = ndim-1

      call hsmg_do_wt(mg_solve_r(mg_solve_index(1,mg_fld)),
     $                mg_mask(mg_mask_index(1,mg_fld)),2,2,nzw)

      !        -1
      ! e  := A   r
      !  1         1
      call hsmg_coarse_solve(mg_solve_e(mg_solve_index(1,mg_fld)),
     $                       mg_solve_r(mg_solve_index(1,mg_fld)))

      call hsmg_do_wt(mg_solve_e(mg_solve_index(1,mg_fld)),
     $                mg_mask(mg_mask_index(1,mg_fld)),2,2,nzw)
      time_3 = dnekclock()
      do l = 2,mg_lmax-1
         nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
         ! w   :=  J e
         !            l-1
         call hsmg_intp
     $      (mg_work2,mg_solve_e(mg_solve_index(l-1,mg_fld)),l-1)

         ! e   :=  e  + w
         !  l       l
         do i = 0,nt-1
            mg_solve_e(mg_solve_index(l,mg_fld)+i) =
     $        + mg_solve_e(mg_solve_index(l,mg_fld)+i) + mg_work2(i+1)
         enddo
      enddo
      l = mg_lmax
      nt = mg_nh(l)*mg_nh(l)*mg_nhz(l)*nelv
      ! w   :=  J e
      !            m-1

      call hsmg_intp(mg_work2,
     $   mg_solve_e(mg_solve_index(l-1,mg_fld)),l-1)

      if (if_hybrid.and.istep.eq.1) then
         ! ecrs := E e_c
         call cdabdtp(ecrs,mg_work2,h1,h2,h2inv,1)
         call cmult  (ecrs,rbd1dt,nt)
         copt(1)  = vlsc2(ecrs2,ecrs,nt)
         copt(2)  = vlsc2(ecrs ,ecrs,nt)
         call gop(copt,copw,'+  ', 2)
         copt(1)  = copt(1)/copt(2)
         avg2     = 1./iter
         avg1     = 1.-avg2
         copt2    = avg1*copt2 + avg2*copt(1)
         if(nid.eq.0)write(6,1)istep,iter,rbd1dt,copt(1),copt2,'cpt2'
      endif
      ! e := e + w

      do i = 1,nt
         e(i) = e(i) + copt2*mg_work2(i)
      enddo
      time_4 = dnekclock()
c     print *, 'Did an MG iteration'
c
      taaaa = taaaa + (time_1 - time_0)
      tbbbb = tbbbb + (time_2 - time_1)
      tcccc = tcccc + (time_3 - time_2)
      tdddd = tdddd + (time_4 - time_3)
      teeee = teeee + (time_4 - time_0)
c
c     A typical time breakdown:
c
c  1.3540E+01  5.4390E+01  1.1440E+01  1.2199E+00  8.0590E+01 HSMG time
c
c  ==>  54/80 = 67 % of preconditioner time is in residual evaluation!
c
      call ortho (e)

#ifndef NOTIMER
      tddsl  = tddsl + ( dnekclock()-etime1 )
#endif


      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_setup_mg_nx()
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'HSMG'
      include 'SEMHAT'
      real w(lx1)
      integer nf,nc,nr
      integer nx,ny,nz

      integer mgn2(10)
      save    mgn2
      data    mgn2 / 1, 2, 2, 2, 2, 3, 3, 5, 5, 5/
c     data    mgn2 / 1, 2, 3, 4, 5, 6, 7, 8, 9, 0

c     if (param(82).eq.0) param(82)=2  ! nek default
c     if (np.eq.1)        param(82)=2  ! single proc. too slow
      p82 = 2                          ! potentially variable nxc
      mg_lmax = 3
c     mg_lmax = 4
      if (lx1.eq.4) mg_lmax = 2
c     if (param(79).ne.0) mg_lmax = param(79)
      mgnx1    = p82-1 !1
      mg_nx(1) = mgnx1
      mg_ny(1) = mgnx1
      mg_nz(1) = mgnx1
      if (.not.if3d) mg_nz(1) = 0 

      mgnx2 = 2*(lx2/4) + 1
      if (lx1.eq.5) mgnx2 = 3
c     if (lx1.eq.6) mgnx2 = 3
      if (lx1.le.10) mgnx2 = mgn2(nx1)
      mg_nx(2) = mgnx2
      mg_ny(2) = mgnx2
      mg_nz(2) = mgnx2
      if (.not.if3d) mg_nz(2) = 0 

      mg_nx(3) = mgnx2+1
      mg_ny(3) = mgnx2+1
      mg_nz(3) = mgnx2+1
      if (.not.if3d) mg_nz(3) = 0 

      mg_nx(mg_lmax) = lx1-1
      mg_ny(mg_lmax) = ly1-1
      mg_nz(mg_lmax) = lz1-1

      if (nid.eq.0) write(*,*) 'mg_nx:',(mg_nx(i),i=1,mg_lmax)
      if (nid.eq.0) write(*,*) 'mg_ny:',(mg_ny(i),i=1,mg_lmax)
      if (nid.eq.0) write(*,*) 'mg_nz:',(mg_nz(i),i=1,mg_lmax)
      
      return
      end
c----------------------------------------------------------------------
      subroutine hsmg_index_0 ! initialize index sets
      include 'SIZE'
      include 'HSMG'

      n = lmgn*(lmgs+1)

      call izero( mg_rstr_wt_index      , n )
      call izero( mg_mask_index         , n )
      call izero( mg_solve_index        , n )
      call izero( mg_fast_s_index       , n )
      call izero( mg_fast_d_index       , n )
      call izero( mg_schwarz_wt_index   , n )
      
      return
      end
c----------------------------------------------------------------------
      subroutine outfldn (x,n,txt10,ichk) ! writes into unit=40+ifiled
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(n,n,1,lelt)
      character*10 txt10
c
      integer idum,e
      save    idum
      data    idum /3/
      if (idum.lt.0)   return
      m = 40 + ifield                 ! unit #
c
C
      mtot = n*n*nelv
      if (n.gt.7.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(m,116) txt10,k,ie,xmin,xmax,ichk,time
            write(m,117) 
            do j=n,1,-1
              if (n.eq.2) write(m,102) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.3) write(m,103) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.4) write(m,104) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.5) write(m,105) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.6) write(m,106) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.7) write(m,107) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.8) write(m,108) ((x(i,j,k,e+l),i=1,n),e=1,ne)
            enddo
         enddo
      enddo

C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldn0 (x,n,txt10,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(n,n,1,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /3/
      if (idum.lt.0) return
c
C
      mtot = n*n*nelv
      if (n.gt.7.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(6,116) txt10,k,ie,xmin,xmax,ichk,time
            write(6,117) 
            do j=n,1,-1
              if (n.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,n),e=1,ne)
              if (n.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,n),e=1,ne)
            enddo
         enddo
      enddo

C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
c     if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outflda (x,n,txt10,ichk) ! writes into unit=p130+ifiled
      INCLUDE 'SIZE'                      ! or into std. output for p130<9
      INCLUDE 'TSTEP'                     ! truncated below eps=p131
      INCLUDE 'INPUT'                     ! param(130)
      real x(1)
      character*10 txt10                  ! note: n is not used
c     parameter (eps=1.e-7)
C
      p130 = param(130)
      eps  = param(131)
      if (p130.le.0)    return
      m    = 6
      if (p130.gt.9)  m = p130 + ifield

      ntot = nx1*ny1*nz1*nelfld(ifield)

      xmin = glmin(x,ntot)
      xmax = glmax(x,ntot)
      xavg = glsum(x,ntot)/ntot

      if (abs(xavg).lt.eps) xavg = 0.     ! truncation

      if (nid.eq.0) write(m,10) txt10,ichk,ntot,xavg,xmin,xmax

   10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g14.6)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldan(x,n,txt10,ichk) ! writes x(1:n) into unit=p130+ifiled
      INCLUDE 'SIZE'                      ! or into std. output for 0<p130<9
      INCLUDE 'TSTEP'                     ! truncated below eps=p131
      INCLUDE 'INPUT'
      real x(1)
      character*10 txt10
c     parameter (eps=1.e-7)
C
      p130 = param(130)
      eps  = param(131)
      if (p130.le.0)    return
      m    = 6
      if (p130.gt.9)  m = p130 + ifield

      ntot = n

      xmin = glmin(x,ntot)
      xmax = glmax(x,ntot)
      xavg = glsum(x,ntot)/ntot

      if (abs(xavg).lt.eps) xavg = 0.     ! truncation

      if (nid.eq.0) write(m,10) txt10,ichk,ntot,xavg,xmin,xmax

   10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g11.3)
c  10 format(3X,a10,2i8,' pts, avg,min,max = ',1p3g14.6)
c
      return
      end
c-----------------------------------------------------------------------
