#ifndef CVODE
      subroutine cv_setsize(n_in,nfld_last)

      include 'SIZE'
      if(nid.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
 
      subroutine cv_init

      include 'SIZE'
      if(nid.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end

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

#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
#else
      integer*4 iout,iout_old,ipar
#endif
      ! cvode will not allocate these arrays
      common /cv_rout/ rout(6),rpar(1)
      common /cv_iout/ iout(21),iout_old(21),ipar(1)

      integer*8 nn,i8glsum

      nxyz = nx1*ny1*nz1

      cv_nfld = nfld_last

      ! set local ODE size
      ipar(1) = n_in
      ipar(1) = ipar(1) + nxyz*nelfld(2)
      do i=3,cv_nfld
         ntot = nxyz*nelfld(i)
         ipar(1) = ipar(1) + ntot
      enddo
      ! determine global ODE size 
      nn = ipar(1)
      cv_nglobal = i8glsum(nn,1)

      return
      end

      subroutine cv_init
c
c     Initialize CVODE solver
c
c     Note:  In contrast to the original CVODE version we
c            define cv_nglobal as 'long long int'
c            (on most platforms 8bytes)
c            YOU HAVE TO USE A PATCHED CVODE VERSION
c 
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      ! use TLAG as scratch space
      common /VPTSOL/  dummy(6*lx1*ly1*lz1*lelv)
     &                ,y0(lx1*ly1*lz1*lelt*ldimt)
#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
#else
      integer*4 iout,iout_old,ipar
#endif
      ! cvode will not allocate these arrays
      common /cv_rout/ rout(6),rpar(1)
      common /cv_iout/ iout(21),iout_old(21),ipar(1),cvcomm

      nxyz = nx1*ny1*nz1
      ifcvodeinit   = .false.
      time          = time - dt 
      cv_time       = time

      if(nid.eq.0) write(*,*) 'initializing ODE integrator ...'

      if(cv_nfld.lt.1) then
        if(nid.eq.0) write(6,*)
     &  'ABORT cv_init(): cv_nfld<1!'
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

      ! check for user-supplied Jacobian routine
      if(abs(PARAM(16)).eq.3) then 
        call FCVSPILSSETJAC(1, ier)
        if(nid.eq.0) then 
          write(6,*) '  user-supplied Jacobian enabled'
          write(6,'(A,4x,e7.1)') '   DQ perturbation scaling factor :',
     &                            CV_SIGS
        endif
      endif

      time        = time + dt
      ifcvodeinit = .true.

      if(nid.eq.0) then
        write(6,'(A,i11)')     '   degrees of freedom            : ',
     &                         cv_nglobal
        write(6,'(A,7x,i4)')   '   krylov dimension              : ',
     &                         cv_maxl
        write(6,'(A,7x,f4.1)') '   linear convergence factor     : ',
     &                         cv_delt
        write(6,*) 'done :: initializing ODE integrator'
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cdscal_cvode(igeom)
c
c     Top level driver for CVODE 
c     Integrate the IVP d/dt Y = f(y(t)) using BDF.
c     f denotes the RHS vector and is computed in fcvfun().
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      ! use TLAG as scratch space (this can get tricky!)
      common /VPTSOL/  dummy(6*lx1*ly1*lz1*lelv)
     &                ,y   (lx1*ly1*lz1*lelt*ldimt)
     &                ,vx_ (lx1,ly1,lz1,lelv)
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)
     &                ,vxd_(lxd,lyd,lzd,lelv)
     &                ,vyd_(lxd,lyd,lzd,lelv)
     &                ,vzd_(lxd,lyd,lzd,lelv)

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

      nxyz = nx1*ny1*nz1

      if(istep.eq.1 .or. .not.ifcvodeinit) call cv_init

      if (.not.ifcvodeinit) then
        write(6,*) 'ABORT: cv_init() was not called!'
        call exitt
      endif 

      ! save old solver output
      do i=1,21
         iout_old(i) = iout(i)
      enddo

      ! recompute fine grid velocities 
      if(ifchar) call set_convect_new(vxd,vyd,vzd,vx,vy,vz)

      ! calling cvode solver
      cv_time=0.0
      call store_vel(vx_,vy_,vz_,vxd_,vyd_,vzd_)
      call fcvode(time,cv_time,y,itask,ier)
      call restore_vel(vx_,vy_,vz_,vxd_,vyd_,vzd_)
 
      ! copy the cvode solution (y) back to the internal nek array (t)
      call cv_unpack_sol(y)

      ! print some solver statistics 
      if (nid.eq.0) then
         nfe         = iout(4)-iout_old(4)+iout(20)-iout_old(20)
         nfe_avg     = (1.-1./istep)*nfe_avg + nfe/float(istep)
         nli_nni     = real(iout(20)-iout_old(20))
         nli_nni     = nli_nni / (real(iout(7)-iout_old(7)) + 1e-50)
         nli_nni_avg = (1.-1./istep)*nli_nni_avg + 
     &                 nli_nni/float(istep)
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,f5.1)')
     &    'nsteps: ',iout(3)-iout_old(3) ,
     &    'nfe:    ', nfe                ,
     &    '<nfe>:    ', nfe_avg
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,f5.1)')  
     &    'nni:    ',iout(7)-iout_old(7)  ,
     &    'nli:    ',iout(20)-iout_old(20),
     &    '<nli/nni>:', nli_nni_avg
         write(*,'(13x,a8,i8,3x,a8,i8,3x,a10,i5)') 
     &    'ncfn:   ',iout(6)-iout_old(6)  ,
     &    'ncfl:   ',iout(21)-iout_old(21),  
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
c----------------------------------------------------------------------
      SUBROUTINE FCVJTIMES (V,FJV,TT,Y,FY,H,IPAR,RPAR,WORK,IER)
c
c     Compute Jacobian Vetor product FJV
c
      REAL V(1), FJV(1), TT, Y(1), FY(1), H, RPAR(1), WORK(1)

      INCLUDE 'SIZE'
      INCLUDE 'CVODE'

#ifdef LONGINT8
      integer*8 ipar(1),neql
#else
      integer*4 ipar(1),neql
#endif

      ! set local size
      neql = ipar(1)

      ! compute weighted rms norm ||v||
      sum = 0.0
      do i = 1,NEQL
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

      ! approximate Jacobian by a finite difference quotient 
      siginv = 1./sig
      do i = 1,NEQL
         FJV(i) = (FJV(i) - FY(i))*siginv
      enddo

      ier = 0

      return
      end
c----------------------------------------------------------------------
      subroutine store_vel(vx_,vy_,vz_,vxd_,vyd_,vzd_)
 
      include 'SIZE'
      include 'TOTAL'
 
      real vx_(1),vy_(1),vz_(1)
      real vxd_(1),vyd_(1),vzd_(1)
 
      ntot  = nx1*ny1*nz1*nelv
      ntotd = nxd*nyd*nzd*nelv
 
      ! save velocities
      call copy(vx_,vx,ntot)
      call copy(vy_,vy,ntot)
      if (if3d) call copy(vz_,vz,ntot)
      if(param(99).gt.0) then
        call copy(vxd_,vxd,ntotd)
        call copy(vyd_,vyd,ntotd)
        if (if3d) call copy(vzd_,vzd,ntotd)
      endif
 
      return
      end
c----------------------------------------------------------------------
      subroutine restore_vel(vx_,vy_,vz_,vxd_,vyd_,vzd_)
 
      include 'SIZE'
      include 'TOTAL'
 
      real vx_(1),vy_(1),vz_(1)
      real vxd_(1),vyd_(1),vzd_(1)
 
      ntot  = nx1*ny1*nz1*nelv
      ntotd = nxd*nyd*nzd*nelv
 
      ! save velocities
      call copy(vx,vx_,ntot)
      call copy(vy,vy_,ntot)
      if (if3d) call copy(vz,vz_,ntot)
      if(param(99).gt.0) then
        call copy(vxd,vxd_,ntotd)
        call copy(vyd,vyd_,ntotd)
        if (if3d) call copy(vzd,vzd_,ntotd)
      endif
 
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
      data timel /0.0/
      real tmp(lx1,ly1,lz1,lelt)

      if (abs(time_-timel).gt.1e-14) then
        timel = time_
      else ! in cache, don't need to evaluate again 
        return
      endif

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
 
      ! compute dealiased velocity
      if (param(99).gt.0) call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
 
      return
      end

      subroutine fcvfun (cv_time, y, ydot, ipar, rpar, ier)
c
c     Compute RHS vector ydot (allocated within cvode)
c     NOTE: working array is bq, do not change! 
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'
c      include 'CHEMISTRY'

      real           cv_time,y(1),ydot(1),rpar(1)
      real           h2(lx1,ly1,lz1,lelt)
      common /scrns/ ta(lx1,ly1,lz1,lelt)
     $              ,tb(lx1,ly1,lz1,lelt)
#ifdef LONGINT8
      integer*8 ipar(1)
#else
      integer*4 ipar(1)
#endif
      save etime_rhs,etime_makeq,etime_dssum,etime_axhelm,etime_setprop
      save timel
      data timel /0.0/

      if (abs(cv_time-timel).gt.1e-13 .or. timel.eq.0.0) then
        timel = cv_time
        icalld = 1
        if(nid.eq.0) write(6,10) cv_time
  10                 format(14X,'substepping t=',1pE14.7)
      endif

#ifdef RHS_TIMING
      etime_rhs0 = dnekclock()
      etime_rhs    = 0.0
      etime_makeq  = 0.0
      etime_setprop= 0.0
      etime_dssum  = 0.0
      etime_axhelm = 0.0
#endif

      nxyz = nx1*ny1*nz1
      isd = 1

      call cv_unpack_sol(y)

#ifdef RHS_TIMING
      etime_tmp = dnekclock()
#endif
      ! calculate thermodynamic and transport properties
      ! including kinetics 
      call setprop
#ifdef RHS_TIMING
      etime_setprop = etime_setprop + dnekclock() - etime_tmp
#endif
 
      ntot = nxyz*nelt
      call rzero (h2,ntot)

      if(abs(PARAM(16)).eq.3)  call update_vel(cv_time)

      j = 0
      do ifield=2,cv_nfld
         ntot = nxyz*nelfld(ifield)
         if (.not.iftmsh(ifield)) imesh = 1
         if (     iftmsh(ifield)) imesh = 2
#ifdef RHS_TIMING
         etime_tmp = dnekclock()
#endif
         call makeq
#ifdef RHS_TIMING
         etime_makeq = etime_makeq + dnekclock() - etime_tmp
#endif

#ifdef RHS_TIMING
         etime_tmp = dnekclock()
#endif
         call axhelm  (ta,t(1,1,1,1,ifield-1),vdiff(1,1,1,1,ifield),
     &                 h2,imesh,isd)
#ifdef RHS_TIMING
         etime_axhelm = etime_axhelm + dnekclock()-etime_tmp
#endif
         ! apply Neumann boundary conditions
         call bcneusc (tb,1)                  

         if (iftmsh(ifield)) then
            do i=1,ntot
               ydot(i+j) = ( bq(i,1,1,1,ifield-1)
     $                   -   ta(i,1,1,1) + tb(i,1,1,1) )
     $                   * bintm1(i,1,1,1)*tmask(i,1,1,1,ifield-1)
     $                   / vtrans(i,1,1,1,ifield)
            enddo
         else
            do i=1,ntot
               ydot(i+j) = ( bq(i,1,1,1,ifield-1)
     $                   -   ta(i,1,1,1) + tb(i,1,1,1) )
     $                   * binvm1(i,1,1,1)*tmask(i,1,1,1,ifield-1)
     $                   / vtrans(i,1,1,1,ifield)
            enddo
         endif
         j = j + ntot
      enddo

#ifdef RHS_TIMING
      etime_tmp = dnekclock()
#endif
      if(nelfld(2).ne.nelfld(1)) then
        ifield = 2 ! set right gs-handle
        call dssum (ydot,nx1,ny1,nz1)
        ntott = nxyz*nelt
        ntot  = nxyz*nelv
        if(cv_nfld-2.gt.0) 
     &    call nvec_dssum (ydot(ntott+1),ntot,cv_nfld-2,gsh_fld(1))
#ifdef RHS_TIMING
        etime_dssum = etime_dssum + dnekclock()-etime_tmp
#endif
      else
        ntot = nxyz*nelv
        if(cv_nfld-1.gt.0) 
     &    call nvec_dssum (ydot,ntot,cv_nfld-1,gsh_fld(1))
#ifdef RHS_TIMING
        etime_dssum = etime_dssum + dnekclock()-etime_tmp
#endif
      endif 

#ifdef RHS_TIMING
      if(nid.eq.0) then
         write(*,*) 'etime_setprop = ',etime_setprop
         write(*,*) 'etime_makeq   = ',etime_makeq
         write(*,*) 'etime_axhelm  = ',etime_axhelm
         write(*,*) 'etime_dssum   = ',etime_dssum
         etime_rhs = etime_rhs + dnekclock()- etime_rhs0
         write(*,*) 'etime_rhs     = ',etime_rhs
      endif
#endif

      ier = 0

      return
      end

#endif
