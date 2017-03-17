
      subroutine acc_copy_all_in()
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'GEOM'
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
!$ACC ENTER DATA COPYIN(mg_imask)
!$ACC ENTER DATA COPYIN(tmult,vmult)
!$ACC ENTER DATA CREATE(e,w,r)
!$ACC ENTER DATA COPYIN(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC ENTER DATA COPYIN(mg_work,mg_fast_s,mg_fast_d)

      return
      end

      subroutine acc_copy_all_out()
      include 'SIZE'
      include 'HSMG'            ! Same array space as HSMG
      include 'GEOM'
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
!$ACC EXIT DATA DELETE(mg_imask)
!$ACC EXIT DATA DELETE(e,w,r)
!$ACC EXIT DATA DELETE(mg_jht,mg_jh,mg_rstr_wt,mg_schwarz_wt)
!$ACC EXIT DATA DELETE(mg_work,mg_mask,mg_fast_s,mg_fast_d)

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
      subroutine GLO_GRAD(D,U,U1,U2,U3)

      include 'SIZE'
      integer i,j,k,l,e
      real D (lx1,lx1)
      real U (lx1,ly1,lz1,lelt)
      real U1(lx1,ly1,lz1,lelt)
      real U2(lx1,ly1,lz1,lelt)
      real U3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (D(NX1,NX1))
!$ACC&     PRESENT (U(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U1(NX1,NY1,NZ1,NELT),U2(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U3(NX1,NY1,NZ1,NELT))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    private(tmpu1,tmpu2,tmpu3)
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
                     tmpu1 = tmpu1 + D(i,l)*U(l,j,k,e)
                     tmpu2 = tmpu2 + D(j,l)*U(i,l,k,e)
                     tmpu3 = tmpu3 + D(k,l)*U(i,j,l,e)
                  enddo
                  U1(i,j,k,e) = tmpu1
                  U2(i,j,k,e) = tmpu2
                  U3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end
c---------------------------------------------------------------
      subroutine GLO_2GRAD(D,U1,U2,U3,V1,V2,V3)

      include 'SIZE'
      integer i,j,k,l,e
      real D (lx1,ly1)
      real U1(lx1,ly1,lz1,lelt)
      real U2(lx1,ly1,lz1,lelt)
      real U3(lx1,ly1,lz1,lelt)
      real V1(lx1,ly1,lz1,lelt)
      real V2(lx1,ly1,lz1,lelt)
      real V3(lx1,ly1,lz1,lelt)
      real tmpu1,tmpu2,tmpu3

!$ACC DATA PRESENT (D (NX1,NX1))
!$ACC&     PRESENT (U1(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (U2(NX1,NY1,NZ1,NELT),U3(NX1,NY1,NZ1,NELT))
!$ACC$     PRESENT (V1(NX1,NY1,NZ1,NELT),V2(NX1,NY1,NZ1,NELT))
!$ACC&     PRESENT (V3(NX1,NY1,NZ1,NELT))
!$ACC PARALLEL LOOP COLLAPSE(4) GANG WORKER VECTOR
!$ACC&    private(tmpu1,tmpu2,tmpu3)
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
                     tmpu1 = tmpu1 + D(i,l)*U1(l,j,k,e)
                     tmpu2 = tmpu2 + D(j,l)*U2(i,l,k,e)
                     tmpu3 = tmpu3 + D(k,l)*U3(i,j,l,e)
                  enddo
                  V1(i,j,k,e) = tmpu1
                  V2(i,j,k,e) = tmpu2
                  V3(i,j,k,e) = tmpu3
               enddo
            enddo
         enddo
      enddo
!$ACC END LOOP
!$ACC END DATA
      return
      end

