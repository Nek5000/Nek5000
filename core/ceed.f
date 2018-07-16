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

      return
      end
c-----------------------------------------------------------------------
