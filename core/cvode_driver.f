#ifndef CVODE
      subroutine cv_setsize(nin)

      include 'SIZE'
      if(nio.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine cv_init

      include 'SIZE'
      if(nio.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
c----------------------------------------------------------------------
      subroutine cdscal_cvode

      include 'SIZE'
      if(nio.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
#else
      subroutine cv_setsize(nin)

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
      integer cvcomm
      common /cv_rout/ rout(6),rpar(1)
      common /cv_iout/ iout(21),iout_old(21),ipar(1),cvcomm

      integer sizeOfLongInt
      external sizeOfLongInt

      integer*8 nn,i8glsum

      if (sizeOfLongInt() .ne. LongIntSize) then
        if(nio.eq.0) write(6,*)
     &  'ABORT cv_setsize(): invalid long int size! ', 
     &   sizeOfLongInt(), LongIntSize
        call exitt
      endif

      nxyz = nx1*ny1*nz1

      ! set local ODE size
      nadd = nin
      if (ifdp0dt) nadd = nadd + 1
      ipar(1) = nadd 
      do i = 2,nfield
         if (ifcvfld(i)) then
           ntot = nxyz*nelfld(i)
           ipar(1) = ipar(1) + ntot
         endif
      enddo

      ! check array size is large enough
      if (ipar(1) .gt. cv_lysize) then
        if(nio.eq.0) write(6,*)
     &  'ABORT cv_setsize(): workspace for cvode too small! ' 
        call exitt
      endif 

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

      common /CVWRK1/  y0(cv_lysize) 

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

      meth = 2 ! 1:non-stiff / 2:stiff

      etime1=dnekclock_sync()

      nxyz = nx1*ny1*nz1
      ifcvodeinit   = .false.

      if(nio.eq.0) write(*,*) 'Initializing CVODE ...'

      ! set solver tolerances
      cv_atol = param(18)
      cv_rtol = param(19) 

      if (cv_rtol.le.0 .or. cv_atol.le.0) then
        write(*,'(a,i3)') 'ABORT cv_init(): invalid tolerances!'
        call exitt
      endif

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

      call cvpack(y0)

      ! initialize main solver
      call fcvmalloc(time, y0, meth, itmeth, iatol,
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

c      ! check for user-supplied Jacobian routine
c      if(abs(PARAM(16)).ge.3) then 
c        call FCVSPILSSETJAC(1, ier)
c        if(nio.eq.0) then 
c          write(6,*) '  user-supplied Jacobian enabled'
c          write(6,'(A,4x,e8.1)') '   DQ perturbation scaling factor :',
c     &                            CV_SIGS
c        endif
c      endif

      etime1 = dnekclock_sync() - etime1
      if(nio.eq.0) then
        write(6,'(A,i11)')       '   degrees of freedom             : ',
     &                         cv_nglobal
        write(6,'(A,2(1pe8.1))') '   rel/abs tolerances             : ',
     &                         cv_rtol, cv_atol
        write(6,'(A,7x,i4)')     '   krylov dimension               : ',
     &                         cv_maxl
        write(6,'(A,6x,f5.2)')   '   linear convergence factor      : ',
     &                         cv_delt
        write(6,'(A,7x,30L2)')   '   fields to integrate            : ',
     &                         ifcvfld
        write(6,'(A,g15.3,A,/)') ' done :: initializing ODE integrator',
     &                         etime1, ' sec'
      endif

      ifcvodeinit = .true.

      return
      end
c----------------------------------------------------------------------
      subroutine cdscal_cvode
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

      common /CVWRK1/  y(cv_lysize) 

      common /CVWRK2/  vx_ (lx1,ly1,lz1,lelv) 
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)
     &                ,xm1_(lx1,ly1,lz1,lelv)
     &                ,ym1_(lx1,ly1,lz1,lelv)
     &                ,zm1_(lx1,ly1,lz1,lelv)
     &                ,wx_ (lx1,ly1,lz1,lelv)
     &                ,wy_ (lx1,ly1,lz1,lelv)
     &                ,wz_ (lx1,ly1,lz1,lelv)

#ifdef LONGINT8
      integer*8 iout,iout_old,ipar
#else
      integer*4 iout,iout_old,ipar
#endif
      integer cvcomm
      common /cv_iout/ iout(21),iout_old(21),ipar(1),cvcomm

      real nfe_avg,nli_nni_avg,nli_nni
      save nfe_avg,nli_nni_avg

      data nfe_avg / 0 /
      data nli_nni_avg / 0 /

      real nni_sum,nli_sum,nfe_sum
      save nni_sum,nli_sum,nfe_sum
      data nni_sum / 0 /
      data nli_sum / 0 / 
      data nfe_sum / 0 /


      itask = 1 ! fixed, at least for now

      nxyz = nx1*ny1*nz1
      ntot = nxyz * nelv

      if (.not.ifcvodeinit) then
        write(6,*) 'ABORT: cv_init() was not called!'
        call exitt
      endif 

      ! save old output values
      if (itask.eq.1) then
        do i=1,21
           iout_old(i) = iout(i)
        enddo
      endif

      ! save coord and mesh vel
      if(ifmvbd) then
        call copy(xm1_,xm1,ntot)           
        call copy(ym1_,ym1,ntot)           
        if (if3d) call copy(zm1_,zm1,ntot) 
        call copy(wx_,wx,ntot)           
        call copy(wy_,wy,ntot)           
        if (if3d) call copy(wz_,wz,ntot) 
      endif

      ! save velocities                  
      call copy(vx_,vx,ntot)             
      call copy(vy_,vy,ntot)            
      if (if3d) call copy(vz_,vz,ntot)  

      ! call solver
      time_ = time
      if (itask.eq.3) then
        if(istep.gt.1) then 
          call cvpack(y)
          call fcvreinit(timef,y,iatol,cv_rtol,cv_atol,ier)
          if (ier .ne. 0) then
            write(*,'(a,i3)') 'ABORT cv_init(): fcvsetrin ier=', ier
            call exitt
          endif
          cv_timel = 0
        endif
        call fcvsetrin('STOP_TIME',time,ier)
      else
        call fcvsetrin('MAX_STEP',2*dt,ier)
      endif
      call fcvode(time,tout,y,itask,ier)
      time = time_

      ! restore velocities               
      call copy(vx,vx_,ntot)             
      call copy(vy,vy_,ntot)             
      if (if3d) call copy(vz,vz_,ntot) 
      if (param(99).gt.0) call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
 
      if(ifmvbd) then        
        ! restore coord and mesh vel 
        call copy(xm1,xm1_,ntot)           
        call copy(ym1,ym1_,ntot)          
        if (if3d) call copy(zm1,zm1_,ntot) 
        call copy(wx,wx_,ntot)           
        call copy(wy,wy_,ntot)           
        if (if3d) call copy(wz,wz_,ntot) 

        call cv_eval_geom                      
      endif

      ! copy the cvode solution (y) back into internal array (t)
      call cvunpack(y)

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

      if (ier.lt.0) then
         if (nid.eq.0) then
            write(*,'(a,i3)') 'ABORT: fcvode returned ier =', ier
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

      common /CVWRK1/  FJV_(cv_lysize) 

      ifcdjv = .false. ! fixed for now

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
      subroutine cv_upd_v
c
      include 'SIZE'
      include 'TSTEP'
      include 'DEALIAS'
      include 'SOLN'
      include 'INPUT'
      include 'CVODE'

      common /CVWRK2/  vx_ (lx1,ly1,lz1,lelv) 
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)
     &                ,xm1_(lx1,ly1,lz1,lelv)
     &                ,ym1_(lx1,ly1,lz1,lelv)
     &                ,zm1_(lx1,ly1,lz1,lelv)
     &                ,wx_ (lx1,ly1,lz1,lelv)
     &                ,wy_ (lx1,ly1,lz1,lelv)
     &                ,wz_ (lx1,ly1,lz1,lelv)

      ntot = nx1*ny1*nz1*nelv

      call sumab(vx,vx_,vxlag,ntot,cv_ab,nbd)
      call sumab(vy,vy_,vylag,ntot,cv_ab,nbd)
      if(if3d) call sumab(vz,vz_,vzlag,ntot,cv_ab,nbd)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_upd_w
c
      include 'SIZE'
      include 'TSTEP'
      include 'MVGEOM'
      include 'INPUT'
      include 'CVODE'

      common /CVWRK2/  vx_ (lx1,ly1,lz1,lelv) 
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)
     &                ,xm1_(lx1,ly1,lz1,lelv)
     &                ,ym1_(lx1,ly1,lz1,lelv)
     &                ,zm1_(lx1,ly1,lz1,lelv)
     &                ,wx_ (lx1,ly1,lz1,lelv)
     &                ,wy_ (lx1,ly1,lz1,lelv)
     &                ,wz_ (lx1,ly1,lz1,lelv)

      ntot = nx1*ny1*nz1*nelv

      call sumab(wx,wx_,wxlag,ntot,cv_ab,nbd)
      call sumab(wy,wy_,wylag,ntot,cv_ab,nbd)
      if(if3d) call sumab(wz,wz_,wzlag,ntot,cv_ab,nbd)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_upd_coor
c
      include 'SIZE'
      include 'TSTEP'
      include 'GEOM'
      include 'MVGEOM'
      include 'INPUT'
      include 'CVODE'

      common /CVWRK2/  vx_ (lx1,ly1,lz1,lelv) 
     &                ,vy_ (lx1,ly1,lz1,lelv)
     &                ,vz_ (lx1,ly1,lz1,lelv)
     &                ,xm1_(lx1,ly1,lz1,lelv)
     &                ,ym1_(lx1,ly1,lz1,lelv)
     &                ,zm1_(lx1,ly1,lz1,lelv)
     &                ,wx_ (lx1,ly1,lz1,lelv)
     &                ,wy_ (lx1,ly1,lz1,lelv)
     &                ,wz_ (lx1,ly1,lz1,lelv)

      COMMON /SCRSF/ dtmp(lx1*ly1*lz1*lelv)

      ntot = nx1*ny1*nz1*nelv

      call sumab(dtmp,wx_,wxlag,ntot,cv_abmsh,nab)
      call add3 (xm1,xm1_,dtmp,ntot)

      call sumab(dtmp,wy_,wylag,ntot,cv_abmsh,nab)
      call add3 (ym1,ym1_,dtmp,ntot)

      if(if3d) then
        call sumab(dtmp,wz_,wzlag,ntot,cv_abmsh,nab)
        call add3 (zm1,zm1_,dtmp,ntot)
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine fcvfun (time_, y, ydot, ipar, rpar, ier)
c
c     Compute RHS function f (called within cvode)
c     NOTE: output is ydot, don't touch y! 
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real time_,y(1),ydot(1),rpar(1)

      real ta(lx1,ly1,lz1,lelt)
      integer ntf      

#ifdef LONGINT8
      integer*8 ipar(1)
#else
      integer*4 ipar(1)
#endif

      time = time_ ! set nek time to internal cvode time    
      nxyz = nx1*ny1*nz1
      ntotv= nxyz*nelv
       
      if (time.ne.cv_timel) then
         call cv_settime     
 
         if(nio.eq.0) write(6,10) istep,time,time-cv_timel
  10       format(4x,i7,2x,'t=',1pE14.7,'  stepsize=',1pE13.4)

         call cv_upd_v
         if (ifmvbd) then
            call cv_upd_coor 
            call cv_eval_geom
            call cv_upd_w
         endif 
         cv_timel = time          
      endif

      call cvunpack(y)          

      j = 1
      ntf = 0
      do ifield=2,nfield
         if (ifcvfld(ifield)) then
           ntot = nxyz*nelfld(ifield)
c           call filter_s(t(1,1,1,1,ifield-1))
           call vprops
           call makeq
           if (iftmsh(ifield)) then                                
              ntf = 1                                          
              call col3(ta, vtrans(1,1,1,1,ifield),bm1, ntot)      
              call dssum(ta, nx1, ny1, nz1)                        
              call col2(ta, bintm1,ntot)                           
              call invcol3(ydot(j), bq(1,1,1,1,ifield-1), ta, ntot) 
           else                                                    
              call invcol3(ydot(j),bq(1,1,1,1,ifield-1),        
     &                     vtrans(1,1,1,1,ifield),ntot)           
           endif     
           j = j + ntot
         endif
      enddo

      j = 1
      if(ntf.eq.1) then
        do ifield = 2,nfield
           if (ifcvfld(ifield)) then                                
             if (.not.iftmsh(ifield)) call dssum (ydot(j),nx1,ny1,nz1)
             j = j + nxyz*nelfld(ifield)
           endif
        enddo
      else ! if no field is a tmesh, use vector dssum
        call nvec_dssum (ydot,ntotv,cv_nfld,gsh_fld(1))
      endif

      j = 1
      do ifield = 2,nfield
         if (ifcvfld(ifield)) then                                
           ntot = nxyz*nelfld(ifield)
           if (.not.iftmsh(ifield)) call col2(ydot(j),binvm1,ntot)
           j = j + ntot
         endif
      enddo

      if (ifdp0dt) then
         call qthermal(.false.,.true.,tlag) ! compute dp0thdt using the 
                                            ! cached source terms stored in tlag
         ydot(j) = dp0thdt
         j = j + 1

         ! add contribution to temp-eqn (assumed to be stored in ifield=2)
         dd = gamma0
         dd = (dd - 1.)/dd
         xfacr= dd * dp0thdt
         do i = 1,ntotv
            ydot(i) = ydot(i) + xfacr/vtrans(i,1,1,1,2)
         enddo
      endif

      call fcvfun_usr(ydot,j)

      j = 1 ! mask out dirichlet points 
      do ifield = 2,nfield
         if (ifcvfld(ifield)) then                                
           ntot = nxyz*nelfld(ifield)
           call col2(ydot(j),tmask(1,1,1,1,ifield-1),ntot)
           j = j + ntot
         endif
      enddo

      ier = 0

      return
      end

#endif
c----------------------------------------------------------------------
      subroutine cv_eval_geom

      call glmapm1
      call geodat1
      call geom2
      call volume
      call setinvm
      call setdef

      return
      end

c----------------------------------------------------------------------
      subroutine cv_settime

      include 'SIZE'
      include 'TSTEP'
      include 'CVODE'

      cv_dtNek = time - timef    

      cv_dtlag(1) = cv_dtNek 
      cv_dtlag(2) = dtlag(2)
      cv_dtlag(3) = dtlag(3)

      call rzero(cv_abmsh,3)
      call setabbd (cv_abmsh,cv_dtlag,nab,1) ! why is nabmsh wrong, use nab for now
      do i = 1,3
         cv_abmsh(i) = cv_dtNek*cv_abmsh(i) 
      enddo

      call rzero(cv_ab,3)
      call setabbd (cv_ab,cv_dtlag,nbd,nbd)

      return
      end
c----------------------------------------------------------------------
      subroutine cvunpack(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do ifield = 2,nfield
         if (ifcvfld(ifield)) then
           ntot = nxyz*nelfld(ifield)
           call copy   (t(1,1,1,1,ifield-1),y(j),ntot)
           call bcdirsc(t(1,1,1,1,ifield-1))
           j = j + ntot
         endif
      enddo

      if (ifdp0dt) then
         p0th = y(j)
         j = j + 1 
      endif

      call cvunpack_usr(y,j) 

      return
      end
c----------------------------------------------------------------------
      subroutine cvpack(y)
c
c     copy the internal nek array (t) to the cvode solution (y)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do ifield = 2,nfield
         if (ifcvfld(ifield)) then
           ntot = nxyz*nelfld(ifield)
           call copy (y(j),t(1,1,1,1,ifield-1),ntot)
           j = j + ntot
         endif
      enddo

      if (ifdp0dt) then
         y(j) = p0th
         j = j + 1
      endif

      call cvpack_usr(y,j)

      return
      end
