
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
 10      format(a4,i,1p1e15.4)
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

