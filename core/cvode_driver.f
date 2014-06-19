#ifndef CVODE
      subroutine cv_setsize(n_in,nfld_last)

      include 'SIZE'
      if(nid.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine cv_init

      include 'SIZE'
      if(nid.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine cdscal_cvode(igeom)

      include 'SIZE'
      if(nid.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
#else
      subroutine cv_setsize(n_in,nfld_last)

      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      integer LongIntSize
#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
      data  LongIntSize / 8 /
#else
      integer*4 iout,iout_old,ipar
      data LongIntSize / 4 /
#endif
      ! cvode will not allocate these arrays
      common /cv_rout/ rout(6),rpar(1)
      common /cv_iout/ iout(21),iout_old(21),ipar(1)

      integer sizeOfLongInt
      external sizeOfLongInt

      integer*8 nn,i8glsum

      if (sizeOfLongInt() .ne. LongIntSize) then
        if(nid.eq.0) write(6,*)
     &  'ABORT cv_setsize(): invalid long int size! ', 
     &   sizeOfLongInt(), LongIntSize
        call exitt
      endif

      nxyz = nx1*ny1*nz1

      cv_nfld = nfld_last

      ! set local ODE size
      ipar(1) = n_in
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         ipar(1) = ipar(1) + ntot
      enddo
      ! determine global ODE size 
      nn = ipar(1)
      cv_nglobal = i8glsum(nn,1)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_init
c
c     Initialize CVODE
c
c     Note:  In contrast to the default CVODE version
c            cv_nglobal is defined as 'long long int'
c 
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      ! use TLAG as scratch space (be careful!)
      common /VPTSOL/  dummy(6*lx1*ly1*lz1*lelv)
     &                ,y0(lx1*ly1*lz1*lelt*ldimt+1)
#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
#else
      integer*4 iout,iout_old,ipar
#endif
      ! cvode will not allocate these arrays
      integer cvcomm
      common /cv_rout/ rout(6),rpar(1)
      common /cv_iout/ iout(21),iout_old(21),ipar(1),cvcomm

      real*8 etime1

      etime1=dnekclock_sync()

      nxyz = nx1*ny1*nz1
      ifcvodeinit   = .false.
      cv_time       = time

      if(nio.eq.0) write(*,*) 'initializing ODE integrator ...'

      if(cv_nfld.lt.2) then
        if(nid.eq.0) write(6,*)
     &  'ABORT cv_init(): invalid cv_nfld index!', cv_nfld
        call exitt
      endif

      if(ldimt.lt.2) then
        if(nid.eq.0) write(6,*)
     &  'ABORT cv_init(): ldimt has to be >= 2!', ldimt
        call exitt
      endif

      ! set dirichlet BCs
      ! note: these values will never change in time (ydot=0)
      !       time varying dirichlet bcs are not supported for now
      do ifield=2,cv_nfld
        call bcdirsc (t(1,1,1,1,ifield-1))
      enddo

      ! set solver tolerances
      cv_rtol = tolrel
      cv_atol = tolabs

      ! initialize vector module
      call create_comm(cvcomm)
#ifdef MPI
      call fnvinitp(cvcomm, 1, ipar(1), cv_nglobal, ier)
#else
      call fnvinits(1, ipar(1), ier)
#endif
      if (ier.ne.0) then
        write(*,'(a,i3)') 'ABORT cv_init(): fnvinitp ier=', ier
        call exitt
      endif

      ! set initial integrator state
      call cv_pack_sol(y0)

      ! initialize main solver
      call fcvmalloc(cv_time, y0, meth, itmeth, iatol,
     &               cv_rtol, cv_atol, iout, rout, ipar, rpar, ier)
      if (ier.ne.0) then
        write(*,'(a,i3)') 'ABORT cv_init(): fcvmalloc ier=', ier
        call exitt
      endif

      ! initialize linear solver
      call fcvspGMR(0,igstype,cv_maxl,cv_delt,ier)
      if (ier .ne. 0) then
        write(*,'(a,i3)') 'ABORT cv_init(): fcvspgmr ier=', ier
        call exitt
      endif

      ! limit time integration order
c      call fcvsetiin('MAX_ORD',3,ier)
c      if (ier .ne. 0) then
c        write(*,'(a,i3)') 'ABORT cv_init(): fcvsetiin ier=', ier
c        call exitt
c      endif

       ! limit maximum integration step size
c      call fcvsetrin('MAX_STEP',dt,ier)
c      if (ier .ne. 0) then
c        write(*,'(a,i3)') 'ABORT cv_init(): fcvsetrin ier=', ier
c        call exitt
c      endif

      ! check for user-supplied Jacobian routine
      if(abs(PARAM(16)).ge.3) then 
        call FCVSPILSSETJAC(1, ier)
        if(nio.eq.0) then 
          write(6,*) '  user-supplied Jacobian enabled'
          write(6,'(A,4x,e7.1)') '   DQ perturbation scaling factor :',
     &                            CV_SIGS
        endif
      endif

      etime1 = dnekclock_sync() - etime1
      if(nio.eq.0) then
        write(6,'(A,i11)')     '   degrees of freedom            : ',
     &                         cv_nglobal
        write(6,'(A,7x,i4)')   '   krylov dimension              : ',
     &                         cv_maxl
        write(6,'(A,6x,f5.2)') '   linear convergence factor     : ',
     &                         cv_delt
        write(6,'(A,7x,i4)')   '   last field index to integrate : ',
     &                         cv_nfld
        write(6,'(A,g15.3,A,/)') ' done :: initializing ODE integrator',
     &                         etime1, ' sec'
      endif

      ifcvodeinit = .true.

      return
      end
c----------------------------------------------------------------------
      subroutine cdscal_cvode(igeom)
c
c     Top level driver for the ODE integrator CVODE
c     webpage: https://computation.llnl.gov/casc/sundials 
c
c     Integrate the IVP d/dt[y] = f(y(t),t); y(t=t0) := f0
c     using BDF(stiff) or AM(non-stiff).
c     f denotes the RHS function and is evaluated in fcvfun().
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      ! use TLAG as scratch space (be careful!)
      common /VPTSOL/  dummy(6*lx1*ly1*lz1*lelv)
     &                ,y   (lx1*ly1*lz1*lelt*ldimt)
     &                ,vx_ (lx1,ly1,lz1,lelv)
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)

#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
#else
      integer*4 iout,iout_old,ipar
#endif
      common /cv_iout/ iout(21),iout_old(21),ipar(1)

      real nfe_avg,nli_nni_avg,nli_nni
      save nfe_avg,nli_nni_avg

      data nfe_avg / 0 /
      data nli_nni_avg / 0 /

      real nni_sum,nli_sum,nfe_sum
      save nni_sum,nli_sum,nfe_sum
      data nni_sum / 0 /
      data nli_sum / 0 / 
      data nfe_sum / 0 /

      nxyz = nx1*ny1*nz1

      if (.not.ifcvodeinit) then
        write(6,*) 'ABORT: cv_init() was not called!'
        call exitt
      endif 

      ! save old output values
      do i=1,21
         iout_old(i) = iout(i)
      enddo

      ! recompute fine grid velocities
      if(ifchar) call set_convect_new(vxd,vyd,vzd,vx,vy,vz)

      cv_time=0.0
      call store_vel(vx_,vy_,vz_)
      call fcvode(time,cv_time,y,itask,ier)
      call set_vel(vx_,vy_,vz_)
      call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
 
      ! copy the cvode solution (y) back into internal array (t)
      call cv_unpack_sol(y)

      ! print solver statistics 
      if (nid.eq.0) then
         nstp        = iout(3)-iout_old(3) 
         nfe         = iout(4)-iout_old(4)+iout(20)-iout_old(20)
         nni         = iout(7)-iout_old(7)  
         nli         = iout(20)-iout_old(20)
         nfe_sum     = nfe_sum + nfe
         nfe_avg     = (1.-1./istep)*nfe_avg + float(nfe)/istep
         nli_sum     = nli_sum + nli
         nni_sum     = nni_sum + nni
         nli_nni     = float(nli)/max(1,nni)
         nli_nni_avg = (1.-1./istep)*nli_nni_avg + nli_nni/istep
         
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,f5.1)')
     &    'nsteps: ',   nstp     ,
     &    'nfe:    ',   nfe      ,
     &    '<nfe>:    ', nfe_avg
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,f5.1)')  
     &    'nni:    ',   nni      ,
     &    'nli:    ',   nli      ,
     &    '<nli/nni>:', nli_nni_avg
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,i5)') 
     &    'ncfn:   '  ,iout(6)-iout_old(6),
     &    'ncfl:   '  ,iout(21)-iout_old(21),  
     &    'netf:     ',iout(5)-iout_old(5) 
      endif

      if (time.ne.cv_time .or. ier.lt.0) then
         if (nid.eq.0) then
            write(*,*) 'ABORT: integration error tout=',cv_time
            write(*,'(a,i3)') 'fcvode returned ier =', ier
         endif
         call exitt
      endif

      return
      end
c----------------------------------------------------------------------
      SUBROUTINE FCVJTIMES (V,FJV,TT,Y,FY,H,IPAR,RPAR,WORK,IER)
c
c     Compute Jacobian Vetor product FJV
c
      REAL V(1), FJV(1), TT, Y(1), FY(1), H, RPAR(1), WORK(1)

      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'CVODE'

#ifdef LONGINT8
      integer*8 ipar(1),neql
#else
      integer*4 ipar(1),neql
#endif
      logical ifcdjv
      data    ifcdjv / .false. /

      common /VPTSOL/  dummy(6*lx1*ly1*lz1*lelv)
     &                ,FJV_(lx1*ly1*lz1*lelt*ldimt)

      if(abs(PARAM(16)).eq.3.1) ifcdjv = .true.

      ! set local size
      neql = ipar(1)

      ! compute weighted rms norm ||v||
      sum = 0.0
      do i = 1,neql
         EWT = 1./(cv_rtol*abs(Y(i)) + cv_atol)   
         dnorm = V(i)*EWT
         sum = sum + dnorm*dnorm
      enddo
      sum = sqrt(glsum(sum,1)/cv_nglobal)
     
      ! set perturbation sig to 1/||v||
      sig =  1./sum

      ! scale perturbation with user supplied constant
      sig = cv_sigs * sig

      ! set FJV = f(t, y + sigs*v/||v||)
      do i = 1,NEQL
         WORK(i) = Y(i) + sig*V(i)
      enddo
      call FCVFUN(TT,WORK,FJV,IPAR,RPAR,IER)

      if(ifcdjv) then ! approximate Jacobian by 2nd-order fd quotient
        do i = 1,NEQL
           WORK(i) = Y(i) - sig*V(i)
        enddo
        call FCVFUN(TT,WORK,FJV_,IPAR,RPAR,IER)
         siginv = 1./(2*sig)
        do i = 1,NEQL
c           FJV(i) = FJV(i)*siginv - FJV_(i)*siginv
           FJV(i) = (FJV(i) - FJV_(i))*siginv
        enddo
      else ! approximate Jacobian by 1st-order fd quotient 
        siginv = 1./sig
        do i = 1,NEQL
c           FJV(i) = (FJV(i) - FY(i))*siginv
           FJV(i) = FJV(i)*siginv - FY(i)*siginv
        enddo
      endif

      ier = 0

      return
      end
c----------------------------------------------------------------------
      subroutine store_vel(vx_,vy_,vz_)
 
      include 'SIZE'
      include 'TOTAL'
 
      real vx_(1),vy_(1),vz_(1)
 
      ntot  = nx1*ny1*nz1*nelv
 
      ! save velocities
      call copy(vx_,vx,ntot)
      call copy(vy_,vy,ntot)
      if (if3d) call copy(vz_,vz,ntot)
 
      return
      end
c----------------------------------------------------------------------
      subroutine set_vel(vx_,vy_,vz_)
 
      include 'SIZE'
      include 'TOTAL'
 
      real vx_(1),vy_(1),vz_(1)
 
      ntot  = nx1*ny1*nz1*nelv
 
      ! save velocities
      call copy(vx,vx_,ntot)
      call copy(vy,vy_,ntot)
      if (if3d) call copy(vz,vz_,ntot)
 
      return
      end
c----------------------------------------------------------------------
      subroutine update_vel(time_)
c
c     extrapolated volocity at t=time_
c
      include 'SIZE'
      include 'TOTAL'
      ! use TLAG as scratch space (this can get tricky!)
      common /VPTSOL/  ttmp(6*lx1*ly1*lz1*lelv)
     &                ,y   (lx1*ly1*lz1*lelt*ldimt)
     &                ,vx_ (lx1,ly1,lz1,lelt)
     &                ,vy_ (lx1,ly1,lz1,lelt)
     &                ,vz_ (lx1,ly1,lz1,lelt)

      real dtlag_(3),ab_(3)
      real timel
      save timel
      data timel /-1.0/
      real tmp(lx1,ly1,lz1,lelt)

      if (time_ .eq. timel) return

c      if(nio.eq.0) write(6,*) 'recompute extrapolated velocity'
      timel = time_

      ntot = nx1*ny1*nz1*nelv
      ! restore velocities from previous time step 
      call copy(vx,vx_,ntot)
      call copy(vy,vy_,ntot)
      if(if3d) call copy(vz,vz_,ntot)

      ! compute extrapolate velocities
      dt_ = time_ - (time-dt)
      dtlag_(1) = dt_
      dtlag_(2) = dtlag(2)
      dtlag_(3) = dtlag(3)
      N = 3
      if(istep.le.2) N = istep
      call rzero(ab_,3)
      call setabbd (ab_,dtlag_,N,NBD)
      ab0 = ab_(1)
      ab1 = ab_(2)
      ab2 = ab_(3)
      call add3s2(tmp,vxlag(1,1,1,1,1),vxlag(1,1,1,1,2),ab1,ab2,ntot)
      call add2s1(vx,tmp,ab0,ntot)
      call add3s2(tmp,vylag(1,1,1,1,1),vylag(1,1,1,1,2),ab1,ab2,ntot)
      call add2s1(vy,tmp,ab0,ntot)
      if (if3d) then
         call add3s2(tmp,vzlag(1,1,1,1,1),vzlag(1,1,1,1,2),ab1,ab2,ntot)
         call add2s1(vz,tmp,ab0,ntot)
      endif
 
      ! update fine grid velocity
      if (param(99).gt.0) call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
 
      return
      end
c----------------------------------------------------------------------
      subroutine fcvfun (cv_time, y, ydot, ipar, rpar, ier)
c
c     Compute RHS function f (allocated within cvode)
c     NOTE: working array is ydot, do not change! 
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real cv_time,y(1),ydot(1),rpar(1)
#ifdef LONGINT8
      integer*8 ipar(1)
#else
      integer*4 ipar(1)
#endif
      data timel / -1 /
      save timel

      nxyz = nx1*ny1*nz1

      if (cv_time.ne.timel) then
        timel = cv_time
        if(nio.eq.0) write(6,10) cv_time
  10                 format(14X,'substepping t=',1pE14.7)
      endif

      call cv_unpack_sol(y)

      ! extrapolate velocity using CVODE's internal time
      if(abs(PARAM(16)).eq.3) call update_vel(cv_time)

      ntflds = 0
      j = 1
      do ifield=2,cv_nfld
         ntot = nxyz*nelfld(ifield)
         if (iftmsh(ifield)) ntflds = ntflds + 1
         call makeq
         call invcol3(ydot(j),bq(1,1,1,1,ifield-1),
     &               vtrans(1,1,1,1,ifield),ntot)
         j = j + ntot
      enddo

      j = 1
      if(ntflds.gt.0) then
        do ifield = 2,cv_nfld
           call dssum (ydot(j),nx1,ny1,nz1)
           j = j + nxyz*nelfld(ifield)
        enddo
      else ! same gs handle; use vector dssum
        ntot = nxyz*nelv
        if(cv_nfld-1.gt.0) 
     &    call nvec_dssum (ydot,ntot,cv_nfld-1,gsh_fld(1))
      endif

      j = 1
      do ifield=2,cv_nfld
         ntot = nxyz*nelfld(ifield)
         if (iftmsh(ifield)) then
            call col2(ydot(j),bintm1,ntot)
         else
            call col2(ydot(j),binvm1,ntot)
         endif
         j = j + ntot
      enddo

      call add_fcvfun_usr(ydot,j)

      j = 1
      do ifield=2,cv_nfld
         ntot = nxyz*nelfld(ifield)
         call col2(ydot(j),tmask(1,1,1,1,ifield-1),ntot)
         j = j + ntot
      enddo

      ier = 0

      return
      end

#endif
