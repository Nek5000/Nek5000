      SUBROUTINE ESOLVER (RES,H1,H2,H2INV,INTYPE)
      use scruz_mod
C---------------------------------------------------------------------
C
C     Choose E-solver
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
      include 'CTIMER'
c
      REAL RES   (LX2,LY2,LZ2,LELV)
      REAL H1    (LX1,LY1,LZ1,LELV)
      REAL H2    (LX1,LY1,LZ1,LELV)
      REAL H2INV (LX1,LY1,LZ1,LELV)
      real, pointer :: wk1(:), wk2(:), wk3(:)

      real kwave2

      ioff = 1
      ltot2 = lx2*ly2*lz2*lelv
      wk1(1:lx2*ly2*lz2*lelv) => cb_scruz(ioff : ioff + ltot2 - 1)
      ioff = ioff + ltot2
      wk2(1:lx2*ly2*lz2*lelv) => cb_scruz(ioff : ioff + ltot2 - 1)
      ioff = ioff + ltot2
      wk3(1:lx2*ly2*lz2*lelv) => cb_scruz(ioff : ioff + ltot2 - 1)

      call ortho(res) !Ensure that residual is orthogonal to null space

      neslv=neslv+1
      etime1=dnekclock()

      if (.not. ifsplit) then
         if (param(42).eq.1) then
            CALL UZAWA (RES,H1,H2,H2INV,INTYPE,ICG)
         else
            call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
         endif
      else
         WRITE(6,*) 'ERROR: E-solver does not exist PnPn'
         CALL EXITT
      ENDIF

      teslv=teslv+(dnekclock()-etime1)

      RETURN
      END
c-----------------------------------------------------------------------
      subroutine dmp_map(imap)
      use scruz_mod
      use ivrtx_mod
      use, intrinsic :: iso_c_binding, only : c_loc, c_f_pointer
c
c     Dump map file and element center point
c
      include 'SIZE'
      include 'TOTAL'

      real,    pointer :: xbar(:,:)
      integer, pointer :: ibar(:)
      integer*8, pointer :: vertex(:)
      integer imap(nelgt)

      integer e,eg

      xbar(1:ldim,1:lelt) => cb_scruz(1 : ldim*lelt)
      vertex(1:(2**ldim)*lelt) => cb_ivrtx(1:(2**ldim)*lelt)
      call c_f_pointer(c_loc(cb_scruz(ldim*lelt+1)), ibar, [lelt])

      nxb = (lx1+1)/2
      nyb = (ly1+1)/2
      nzb = (lz1+1)/2
      
      do e=1,nelt
         xbar(ldim,e) = zm1(nxb,nyb,nzb,e)
         xbar(1   ,e) = xm1(nxb,nyb,nzb,e)
         xbar(2   ,e) = ym1(nxb,nyb,nzb,e)
         eg           = lglel(e)
         ibar(e)      = imap(eg)
      enddo
      call p_outvec_ir(ibar,xbar,ldim,'mpxyz.dat')

      return
      end
c-----------------------------------------------------------------------
      subroutine p_outvec_ir(ia,a,lda,name9)
      include 'SIZE'
      include 'TOTAL'

      integer ia(1)
      real    a(lda,1)
      character*9 name9

      parameter (lbuf=50)
      common /scbuf/ buf(lbuf)
      integer ibuf(10),e,eg
      equivalence (buf,ibuf)

      if (nid.eq.0) then
         open(unit=49,file=name9)
         write(6,*) 'Opening ',name9,' in p_outveci. lda=',lda
      endif

      len = wdsize*(lda+1)
      dum = 0.

      do eg=1,nelgt

         mid   = gllnid(eg)
         e     = gllel (eg)
         mtype = 2000+e

         if (nid.eq.0) then
            if (mid.eq.0) then
               call icopy(buf(1),ia(e),1)
               call  copy(buf(2),a(1,e),lda)
            else
               call csend (mtype,dum,wdsize,mid,nullpid)
               call crecv2 (mtype,buf,len,mid)
            endif
            write(49,49) mid,ibuf(1),(buf(k+1),k=1,lda)
   49       format(2i12,1p3e16.7)
         elseif (nid.eq.mid) then
            call icopy(buf(1),ia(e),1)
            call  copy(buf(2),a(1,e),lda)
            call crecv2 (mtype,dum,wdsize,0)
            call csend (mtype,buf,len,node0,nullpid)
         endif
      enddo

      if (nid.eq.0) then
         close(49)
         write(6,*) 'Done writing to ',name9,' p_outveci.'
      endif

      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
