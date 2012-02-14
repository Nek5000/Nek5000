      SUBROUTINE ESOLVER (RES,H1,H2,H2INV,INTYPE)
C---------------------------------------------------------------------
C
C     Choose E-solver
C
C--------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'ESOLV'
      INCLUDE 'INPUT'
C
      REAL RES   (LX2,LY2,LZ2,LELV)
      REAL H1    (LX1,LY1,LZ1,LELV)
      REAL H2    (LX1,LY1,LZ1,LELV)
      REAL H2INV (LX1,LY1,LZ1,LELV)
      common /scruz/ wk1(lx2*ly2*lz2*lelv)
     $             , wk2(lx2*ly2*lz2*lelv)
     $             , wk3(lx2*ly2*lz2*lelv)
c
      include 'CTIMER'
c
#ifndef NOTIMER
      if (icalld.eq.0) teslv=0.0
      icalld=icalld+1
      neslv=icalld
      etime1=dnekclock()
#endif
C
c     write(6,*) solver_type,' solver type',iesolv
      IF (IESOLV.EQ.1) THEN
         if (solver_type.eq.'fdm') then
            ntot2 = nx2*ny2*nz2*nelv
            call gfdm_pres_solv  (wk1,res,wk2,wk3)
            call copy            (res,wk1,ntot2)
         else
            if (param(42).eq.1.or.solver_type.eq.'pdm') then
               CALL UZAWA (RES,H1,H2,H2INV,INTYPE,ICG)
            else
               call uzawa_gmres(res,h1,h2,h2inv,intype,icg)
            endif
         endif
      ELSE
         WRITE(6,*) 'ERROR: E-solver does not exist',iesolv
         WRITE(6,*) 'Stop in ESOLVER'
         CALL EXITT
      ENDIF
#ifndef NOTIMER
      teslv=teslv+(dnekclock()-etime1)
#endif
      RETURN
      END
      SUBROUTINE ESTRAT
C---------------------------------------------------------------------------
C
C     Decide strategy for E-solver
C
C---------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
C
      IESOLV = 1
      if (ifsplit) iesolv=0
c
      solver_type='itr'
      if (param(116).ne.0) solver_type='fdm'
c     if (param(90).ne.0)  solver_type='itn'
c
c     The following change recognizes that geometry is logically 
c     tensor-product, but deformed:  pdm = Preconditioner is fdm
c
      if (param(59).ne.0.and.solver_type.eq.'fdm') solver_type='pdm'

      IF (ISTEP.LT.2.AND.NID.EQ.0) WRITE(6,10) IESOLV,solver_type
   10 FORMAT(2X,'E-solver strategy: ',I2,1X,A)



C
      RETURN
      END
      SUBROUTINE EINIT
C-----------------------------------------------------------------------------
C
C     Initialize E-solver
C
C-----------------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'ESOLV'
      COMMON /SCRHI/  H2INV (LX1,LY1,LZ1,LELV)
      LOGICAL IFNEWP
C
      CALL ESTRAT 
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine dmp_map(imap)
c
c     Dump map file and element center point
c
      include 'SIZE'
      include 'TOTAL'

      common /ivrtx/ vertex ((2**ldim)*lelt)
      common /scruz/ xbar(ldim,lelt),ibar(lelt)
      integer vertex
      integer imap(nelgt)

      integer e,eg

      nxb = (nx1+1)/2
      nyb = (ny1+1)/2
      nzb = (nz1+1)/2
      
      do e=1,nelt
         xbar(ndim,e) = zm1(nxb,nyb,nzb,e)
         xbar(1   ,e) = xm1(nxb,nyb,nzb,e)
         xbar(2   ,e) = ym1(nxb,nyb,nzb,e)
         eg           = lglel(e)
         ibar(e)      = imap(eg)
      enddo
      call p_outvec_ir(ibar,xbar,ndim,'mpxyz.dat')

      return
      end
c-----------------------------------------------------------------------
      subroutine p_outvec_ir(ia,a,lda,name9)
      integer ia(1)
      real    a(lda,1)
      character*9 name9

      include 'SIZE'
      include 'TOTAL'

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
         mtype = 2000+eg

         if (nid.eq.0) then
            if (mid.eq.0) then
               call icopy(buf(1),ia(e),1)
               call  copy(buf(2),a(1,e),lda)
            else
               call csend (mtype,dum,wdsize,mid,nullpid)
               call crecv (mtype,buf,len)
            endif
            write(49,49) mid,ibuf(1),(buf(k+1),k=1,lda)
   49       format(2i12,1p3e16.7)
         elseif (nid.eq.mid) then
            call icopy(buf(1),ia(e),1)
            call  copy(buf(2),a(1,e),lda)
            call crecv (mtype,dum,wdsize)
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

