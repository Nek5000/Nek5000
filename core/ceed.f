c-----------------------------------------------------------------------
      subroutine setup_ceed(ctx,qdata,q,u1,u2,u3,v1,v2,ierr)
      integer q,ierr
      real*8 ctx(1)
      real*8 u1(1)
      real*8 u2(1)
      real*8 u3(1)
      real*8 qdata(6,q)
      real*8 jacmq
      real*8 scl

      do i=1,q
        a11=u2(i+q*0)
        a21=u2(i+q*1)
        a31=u2(i+q*2)

        a12=u2(i+q*3)
        a22=u2(i+q*4)
        a32=u2(i+q*5)

        a13=u2(i+q*6)
        a23=u2(i+q*7)
        a33=u2(i+q*8)

        jacmq = a11*(a22*a33-a23*a32)
     $        - a12*(a21*a33-a23*a31)
     $        + a13*(a21*a32-a22*a31)

        g11 = (a22*a33-a23*a32)/jacmq
        g12 = (a13*a32-a33*a12)/jacmq
        g13 = (a12*a23-a22*a13)/jacmq

        g21 = (a23*a31-a21*a33)/jacmq
        g22 = (a11*a33-a31*a13)/jacmq
        g23 = (a13*a21-a23*a11)/jacmq

        g31 = (a21*a32-a22*a31)/jacmq
        g32 = (a12*a31-a32*a11)/jacmq
        g33 = (a11*a22-a21*a12)/jacmq

        scl = abs(u3(i)*jacmq)

        qdata(1,i) = scl*(g11*g11+g12*g12+g13*g13) ! Grr
        qdata(2,i) = scl*(g11*g21+g12*g22+g13*g23) ! Grs
        qdata(3,i) = scl*(g11*g31+g12*g32+g13*g33) ! Grt
        qdata(4,i) = scl*(g21*g21+g22*g22+g23*g23) ! Gss
        qdata(5,i) = scl*(g21*g31+g22*g32+g23*g33) ! Gst
        qdata(6,i) = scl*(g31*g31+g32*g32+g33*g33) ! Gtt

      enddo

      ierr=0
      end
c-----------------------------------------------------------------------
      subroutine diffusion_ceed(ctx,qdata,q,u1,u2,u3,v1,v2,ierr)
      integer q,ierr
      real*8 ctx
      real*8 u1(1)
      real*8 u2(1)
      real*8 u3(1)
      real*8 v1(1)
      real*8 v2(1)
      real*8 qdata(6,q)

      write(6,*) 'Thilina is in CEED'
      do i=1,q
        v2(i+0*q)=
     $     qdata(1,i)*u2(i)+qdata(2,i)*u2(i+q)+qdata(3,i)*u2(i+2*q)
        v2(i+1*q)=
     $     qdata(2,i)*u2(i)+qdata(4,i)*u2(i+q)+qdata(5,i)*u2(i+2*q)
        v2(i+2*q)=
     $     qdata(3,i)*u2(i)+qdata(5,i)*u2(i+q)+qdata(6,i)*u2(i+2*q)
      enddo

      ierr=0
      end
c-----------------------------------------------------------------------
      subroutine axhm1_ceed(pap,ap1,p1,h1,h2,ceed,op_diffusion,
     $  vec_ap1,
     $  vec_p1,vec_qdata)

      include 'ceedf.h'

c     Vector conjugate gradient matvec for solution of uncoupled
c     Helmholtz equations

      include 'SIZE'
      include 'TOTAL'

      parameter (lzq=lx1)
      parameter (lx=lx1*ly1*lz1,lg=3+3*(ldim-2),lq=lzq**ldim)

      real ap1(lx,lelt)
      real  p1(lx,lelt)
      integer ceed,op_diffusion,vec_p1,vec_ap1,vec_qdata

      common /ctmp1/ ur,us,ut,wk
      real ur(lq),us(lq),ut(lq),wk(lq)

      integer e

!$acc data present (ap1(lx,lelt))
!$acc&     present (p1 (lx,lelt))

      call ceedvectorsetarray(vec_p1,ceed_mem_host,ceed_use_pointer,
     $  p1,err)
      call ceedoperatorapply(op_diffusion,vec_qdata,vec_p1,vec_ap1,
     $  ceed_request_immediate,err)
      call ceedvectorgetarray(vec_ap1,ceed_mem_host,ap1,err)

      pap=0.

      do e=1,nelt
         do i=1,lx
           pap=pap+p1(i,e)*ap1(i,e)
         enddo
      enddo

!$acc end data
      return
      end
c-----------------------------------------------------------------------
      subroutine init_axhm1_ceed(spec,ceed,op_diffusion,op_setup,
     $  vec_p1,vec_ap1,vec_qdata)

      include 'SIZE'
      include 'TOTAL'
      include 'ceedf.h'

      parameter (lzq=lx1)
      parameter (lx=lx1*ly1*lz1,lg=3+3*(ldim-2),lq=lzq**ldim)

c     spec is given, everything else is setup by this routine
      integer ceed,op_diffusion,op_setup
      integer vec_p1,vec_ap1,vec_qdata,vec_coords
      character*64 spec

      integer*8 ndof
      integer err
      integer p,q,ncomp,edof,ldof
      integer basisu,basisx
      integer erstrctu,erstrctx
      integer qdata
      integer qf_diffusion,qf_setup
      integer ii,i,e,ngeo,n
      integer identity(lelt*(lx1**ldim))

      real*8 coords(ldim*lx*lelt)
!$acc enter data copyin(coords(ldim*lx*nelt))

      external diffusion_ceed,setup_ceed

c     Init ceed library
      call ceedinit(trim(spec)//char(0),ceed,err)

      n      = nx1*ny1*nz1*nelt
      nzq    = nx1

c     Create ceed basis for mesh and computation
      p=nx1
      q=p
      ncomp=1
      call ceedbasiscreatetensorh1lagrange(ceed,ldim,3*ncomp,p,q,
     $  ceed_gauss_lobatto,basisx,err)
      call ceedbasiscreatetensorh1lagrange(ceed,ldim,ncomp,p,q,
     $  ceed_gauss_lobatto,basisu,err)

      ncount=0
      do i=1,nelt
      do j=1,lx1**ldim
        ncount = ncount+1
        identity(ncount)=ncount-1
      enddo
      enddo
c     Create ceed element restrictions for mesh and computation
      edof=nx1**ldim
      ldof=edof*nelt*ncomp
      call ceedelemrestrictioncreate(ceed,nelt,edof,ldof,
     $  ceed_mem_host,ceed_use_pointer,identity,
     $  erstrctx,err)

      edof=nx1**ldim
      ldof=edof*nelt*ncomp
      call ceedelemrestrictioncreate(ceed,nelt,edof,ldof,
     $  ceed_mem_host,ceed_use_pointer,identity,
     $  erstrctu,err)

c     Create ceed qfunctions for setupf and diffusionf
      ngeo=(ldim*(ldim+1))/2
      call ceedqfunctioncreateinterior(ceed,1,ncomp,ngeo*8,
     $  ior(ceed_eval_grad,ceed_eval_weight),ceed_eval_none,setupf,
     $  __FILE__
     $  //':setup_ceed'//char(0),qf_setup,err)

      call ceedqfunctioncreateinterior(ceed,1,ncomp,ngeo*8,
     $  ceed_eval_grad,ceed_eval_grad,diffusionf,
     $  __FILE__
     $  //':diffusion_ceed'//char(0),qf_diffusion,err)

c     Create a ceed operator
      call ceedoperatorcreate(ceed,erstrctx,basisx,qf_setup,
     $  ceed_null,ceed_null,op_setup,err)
      call ceedoperatorcreate(ceed,erstrctu,basisu,qf_diffusion,
     $  ceed_null,ceed_null,op_diffusion,err)

c     Create ceed vectors
      call ceedvectorcreate(ceed,ldof,vec_p1,err)
      call ceedvectorcreate(ceed,ldof,vec_ap1,err)
      call ceedoperatorgetqdata(op_setup,vec_qdata,err)

!$acc data present (xm1(nx1,ny1,nz1,nelt))
!$acc&     present (ym1(nx1,ny1,nz1,nelt))
!$acc&     present (zm1(nx1,ny1,nz1,nelt))

!$acc kernels
      ii=0
      do j=0,nelt-1
      do i=1,lx
        ii=ii+1
        coords(i+0*lx+3*j*lx)=xm1(ii,1,1,1)
        coords(i+1*lx+3*j*lx)=ym1(ii,1,1,1)
        coords(i+2*lx+3*j*lx)=zm1(ii,1,1,1)
      enddo
      enddo
!$acc end kernels

      call ceedvectorcreate(ceed,3*n,vec_coords,err)
      call ceedvectorsetarray(vec_coords,ceed_mem_host,
     $  ceed_use_pointer,coords,err)
      call ceedoperatorapply(op_setup,vec_qdata,vec_coords,ceed_null,
     $  ceed_request_immediate,err)

!$acc end data

      end
c-----------------------------------------------------------------------
      subroutine test_ceed(ctx,qdata,q,u1,u2,u3,v1,v2,ierr)
      integer q,ierr
      real*8 ctx(1)
      real*8 u1(1)
      real*8 u2(1)
      real*8 u3(1)
      real*8 qdata(q)
      real*8 jacmq
      real*8 scl

      do i=1,q
        qdata(i)=9.0 
      enddo

      ierr=0
      end
c-----------------------------------------------------------------------
      subroutine openacc_test_ceed()

      include 'SIZE'
      include 'TOTAL'
      include 'ceedf.h'

      parameter (lzq=lx1)
      parameter (lx=lx1*ly1*lz1,lg=3+3*(ldim-2),lq=lzq**ldim)

c     spec is given, everything else is setup by this routine
      integer ceed,op_test
      integer vec_qdata

      integer*8 ndof
      integer err
      integer p,q,ncomp,edof,ldof
      integer basisx
      integer erstrctx
      integer qdata
      integer qf_test
      integer ii,i,e,ngeo,n
      integer identity(lelt*(lx1**ldim))

      real*8 qio(lx*lelt)

      external test_ceed

c     Init ceed library
      call ceedinit('/gpu/occa',ceed,err)

      n      = nx1*ny1*nz1*nelt
      nzq    = nx1

c     Create ceed basis for mesh and computation
      p=nx1
      q=p
      ncomp=1
      call ceedbasiscreatetensorh1lagrange(ceed,ldim,ncomp,p,q,
     $  ceed_gauss_lobatto,basisx,err)

      ncount=0
      do i=1,nelt
      do j=1,lx1**ldim
        ncount = ncount+1
        identity(ncount)=ncount-1
        qio(ncount)=ncount-1
      enddo
      enddo

!$acc enter data copyin(qio(1:ncount))

c     Create ceed element restrictions for mesh and computation
      edof=nx1**ldim
      ldof=edof*nelt*ncomp
      call ceedelemrestrictioncreate(ceed,nelt,edof,ldof,
     $  ceed_mem_host,ceed_use_pointer,identity,
     $  erstrctx,err)

c     Create ceed qfunctions for setupf and diffusionf
      call ceedqfunctioncreateinterior(ceed,1,ncomp,8,
     $  ceed_eval_none,ceed_eval_none,test_ceed,
     $  __FILE__
     $  //':test_ceed'//char(0),qf_test,err)

c     Create a ceed operator
      call ceedoperatorcreate(ceed,erstrctx,basisx,qf_test,
     $  ceed_null,ceed_null,op_test,err)

c     Create ceed vectors
      call ceedoperatorgetqdata(op_test,vec_qdata,err)

c     call ceedvectorcreate(ceed,3*n,vec_coords,err)
      call ceedvectorsetarray(vec_qdata,ceed_mem_device,
     $  ceed_use_pointer,qio,err)
      call ceedoperatorapply(op_test,vec_qdata,ceed_null,ceed_null,
     $  ceed_request_immediate,err)

!$acc exit data copyout(qio(1:ncount))

      end
c-----------------------------------------------------------------------
