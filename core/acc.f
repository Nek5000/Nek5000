
      subroutine acc_copy_all_in()
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'DXYZ'
      include 'GEOM'
      include 'GMRES'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'


      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv),
     $     h2    (lx1,ly1,lz1,lelv)
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
!$ACC ENTER DATA COPYIN(work,work2)
!$ACC ENTER DATA COPYIN(mg_mask,mg_imask,pmask)
!$ACC ENTER DATA COPYIN(tmult,vmult)
!$ACC ENTER DATA CREATE(e,w,r)
!$ACC ENTER DATA COPYIN(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC ENTER DATA COPYIN(mg_work,mg_fast_s,mg_fast_d)
!$ACC ENTER DATA COPYIN(g1m1,g2m1,g3m1,g4m1,g5m1,g6m1,dxm1,dxtm1)
!$ACC ENTER DATA COPYIN(h_gmres,w_gmres,v_gmres,z_gmres)
!$ACC ENTER DATA COPYIN(c_gmres,s_gmres,x_gmres,gamma_gmres)
!$ACC ENTER DATA COPYIN(r_gmres)
!$ACC ENTER DATA COPYIN(ml_gmres,mu_gmres)
!$ACC ENTER DATA COPYIN(h1,h2)


      return
      end

      subroutine acc_copy_all_out()
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'DXYZ'
      include 'GEOM'
      include 'GMRES'
      include 'INPUT'
      include 'MASS'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
      include 'PARALLEL'


      common /scrhi/ h2inv (lx1,ly1,lz1,lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv),
     $     h2    (lx1,ly1,lz1,lelv)
      parameter (lt=lx1*ly1*lz1*lelt)
      common /scrmg/ e(2*lt),w(lt),r(lt)
      parameter (lwk=(lx1+2)*(ly1+2)*(lz1+2))
      common /hsmgw/ work(0:lwk-1),work2(0:lwk-1)
!$ACC EXIT DATA DELETE(work,work2)
!$ACC EXIT DATA DELETE(mg_mask,mg_imask,pmask)
!$ACC EXIT DATA DELETE(e,w,r)
!$ACC EXIT DATA DELETE(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC EXIT DATA DELETE(mg_work,mg_fast_s,mg_fast_d)
!$ACC EXIT DATA DELETE(g1m1,g2m1,g3m1,g4m1,g5m1,g6m1,dxm1,dxtm1)
!$ACC EXIT DATA COPYOUT(h_gmres,w_gmres,v_gmres,z_gmres)
!$ACC EXIT DATA COPYOUT(c_gmres,s_gmres,x_gmres,gamma_gmres)
!$ACC EXIT DATA COPYOUT(r_gmres)
!$ACC EXIT DATA COPYOUT(ml_gmres,mu_gmres)
!$ACC EXIT DATA COPYOUT(h1,h2)

      return
      end

      subroutine print_acc(a,n,text,mode)
      include 'SIZE'
      integer n,i,mode
      real    a(n)
      character*4 text

#ifdef _OPENACC
!$ACC update host(a) if (mode.eq.1)  !copies the data from device to host
#endif
c     if (nid.eq.0) write(6,*) text,(a(i),i=300,n) ! prints data on host
      if (nid.eq.0) then
         do i=1,n
            write(6,10) text,i,a(i)
         enddo
 10      format(a4,i300,1p1e15.4)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine print_acc_int(a,n,text,mode)
      include 'SIZE'
      integer n,i,mode
      integer a(n)
      character*4 text

#ifdef _OPENACC
!$ACC update host(a) if (mode.eq.1)  !copies the data from device to !host
#endif
      if (nid.eq.0) write(6,*) text,(a(i),i=1,n) ! prints data on host

      return
      end
c-----------------------------------------------------------------------
      subroutine global_grad3(d,u,u1,u2,u3)

      include 'SIZE'
      integer i,j,k,l,e
      real d (lx1,lx1)
      real u (lx1,ly1,lz1,lelt)
      real u1(lx1,ly1,lz1,lelt)
      real u2(lx1,ly1,lz1,lelt)
      real u3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (d(nx1,nx1))
!$ACC&     PRESENT (u(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u1(nx1,ny1,nz1,nelt),u2(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u3(nx1,ny1,nz1,nelt))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    PRIVATE(tmpu1,tmpu2,tmpu3)
!dir$ NOBLOCKING
      do e=1,nelt
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  tmpu1 = 0.0
                  tmpu2 = 0.0
                  tmpu3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpu1 = tmpu1 + d(i,l)*u(l,j,k,e)
                     tmpu2 = tmpu2 + d(j,l)*u(i,l,k,e)
                     tmpu3 = tmpu3 + d(k,l)*u(i,j,l,e)
                  enddo
                  u1(i,j,k,e) = tmpu1
                  u2(i,j,k,e) = tmpu2
                  u3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end
c---------------------------------------------------------------
      subroutine global_div3(d,u1,u2,u3,v1,v2,v3)

      include 'SIZE'
      integer i,j,k,l,e
      real d (lx1,ly1)
      real u1(lx1,ly1,lz1,lelt)
      real u2(lx1,ly1,lz1,lelt)
      real u3(lx1,ly1,lz1,lelt)
      real v1(lx1,ly1,lz1,lelt)
      real v2(lx1,ly1,lz1,lelt)
      real v3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (d (nx1,nx1))
!$ACC&     PRESENT (u1(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (u2(nx1,ny1,nz1,nelt),u3(nx1,ny1,nz1,nelt))
!$ACC$     PRESENT (v1(nx1,ny1,nz1,nelt),v2(nx1,ny1,nz1,nelt))
!$ACC&     PRESENT (v3(nx1,ny1,nz1,nelt))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    PRIVATE(tmpu1,tmpu2,tmpu3)
!dir$ NOBLOCKING
      do e=1,nelt
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  tmpu1 = 0.0
                  tmpu2 = 0.0
                  tmpu3 = 0.0
!$ACC LOOP SEQ
                  do l=1,nx1
                     tmpu1 = tmpu1 + d(i,l)*u1(l,j,k,e)
                     tmpu2 = tmpu2 + d(j,l)*u2(i,l,k,e)
                     tmpu3 = tmpu3 + d(k,l)*u3(i,j,l,e)
                  enddo
                  v1(i,j,k,e) = tmpu1
                  v2(i,j,k,e) = tmpu2
                  v3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end

