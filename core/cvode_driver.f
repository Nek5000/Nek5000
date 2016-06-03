#ifndef CVODE
      subroutine cv_setsize(n_in,nfld_last)

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
      subroutine cdscal_cvode(igeom)

      include 'SIZE'
      if(nio.eq.0) write(6,*) 'ABORT: not compiled with CVODE support!'
      call exitt

      return
      end
#else
      subroutine cv_setsize(i1fld,ilfld,n_aux)

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

      cv_1fld = i1fld
      cv_lfld = ilfld

      ! set local ODE size
      ipar(1) = n_aux
      do i=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(i)
         ipar(1) = ipar(1) + ntot
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

      cv_ystart = 1
      cv_1fld   = 2 

      etime1=dnekclock_sync()

      nxyz = nx1*ny1*nz1
      ifcvodeinit   = .false.
      cv_time       = time

      if(nio.eq.0) write(*,*) 'initializing ODE integrator ...'

      if(cv_lfld-cv_1fld .lt. 0) then
        if(nid.eq.0) write(6,*)
     &  'ABORT cv_init(): no fields to integrate!', cv_1fld, cv_lfld
        call exitt
      endif

c      ! check to see if lagged arrarys are large enough
c      if(lorder.lt.3) then
c        if(nid.eq.0) write(6,*)
c     &  'ABORT cv_init(): lorder has to be >= 3!', lorder
c        call exitt
c      endif
c      if(ldimt.lt.1) then
c        if(nid.eq.0) write(6,*)
c     &  'ABORT cv_init(): ldimt has to be >= 3!', ldimt
c        call exitt
c      endif

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

      ! enforce dirichlet bcs, time varying bcs are not supported for now
      do ifield=cv_1fld,cv_lfld
        call bcdirsc (t(1,1,1,1,ifield-1))
      enddo
      call cv_pack_sol(y0(cv_ystart))

      ! initialize main solver
      call fcvmalloc(cv_time, y0(cv_ystart), meth, itmeth, iatol,
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

      ! limit intneral integration step size
      ! to ensure extrapolation v,w is accurate 
      call fcvsetrin('MAX_STEP',dt,ier)
      if (ier .ne. 0) then
        write(*,'(a,i3)') 'ABORT cv_init(): fcvsetrin ier=', ier
        call exitt
      endif

      ! check for user-supplied Jacobian routine
      if(abs(PARAM(16)).ge.3) then 
        call FCVSPILSSETJAC(1, ier)
        if(nio.eq.0) then 
          write(6,*) '  user-supplied Jacobian enabled'
          write(6,'(A,4x,e8.1)') '   DQ perturbation scaling factor :',
     &                            CV_SIGS
        endif
      endif

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
        write(6,'(A,7x,i4)')     '   first field index to integrate : ',
     &                         cv_1fld
        write(6,'(A,7x,i4)')     '   last field index to integrate  : ',
     &                         cv_lfld
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

      nxyz = nx1*ny1*nz1
      ntot = nxyz * nelv

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
      call fcvode(time,cv_time,y(cv_ystart),itask,ier)


      ! restore velocities               
      call copy(vx,vx_,ntot)             
      call copy(vy,vy_,ntot)             
      if (if3d) call copy(vz,vz_,ntot) 
      call set_convect_new(vxd,vyd,vzd,vx,vy,vz)
 
      ! restore coord and mesh vel 
      if(ifmvbd) then         
        call copy(xm1,xm1_,ntot)           
        call copy(ym1,ym1_,ntot)          
        if (if3d) call copy(zm1,zm1_,ntot) 
        call copy(wx,wx_,ntot)           
        call copy(wy,wy_,ntot)           
        if (if3d) call copy(wz,wz_,ntot) 
        call cv_eval_geom                      
      endif

      ! copy the cvode solution (y) back into internal array (t)
      call cv_unpack_sol(y(cv_ystart))

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

      common /CVWRK1/  FJV_(cv_lysize) 

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
      subroutine cv_upd_v
c
      include 'SIZE'
      include 'TSTEP'
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
      subroutine fcvfun (cv_time, y, ydot, ipar, rpar, ier)
c
c     Compute RHS function f (allocated within cvode)
c     NOTE: working array is ydot, do not change! 
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real cv_time,y(1),ydot(1),rpar(1)

      real ta(lx1,ly1,lz1,lelt)
      integer ntf,ntflds(nfield-1)      

#ifdef LONGINT8
      integer*8 ipar(1)
#else
      integer*4 ipar(1)
#endif
      data timel / -1 /
      save timel

      nxyz = nx1*ny1*nz1
      ntotv= nxyz*nelv                  
      call izero(ntflds,nfield-1)       
      call cv_unpack_sol(y)          
        
      if (cv_time.ne.timel) then
         call cv_settime(cv_time)     
 
         if(nio.eq.0) write(6,10) cv_time,timel,cv_time-timel
  10       format(14X,'substepping t=',1pE14.7,1pE14.7,1pE14.7)

         call cv_upd_v
         if (ifmvbd) then
            call cv_upd_coor 
            call cv_eval_geom
            call cv_upd_w
         endif 
 
         timel = cv_time          
      endif

      j = 1
      ntf = 0
      do ifield=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(ifield)
         call makeq
         if (iftmsh(ifield)) then                                
            ntflds(ifield-1) = 1                                 
            ntf = 1                                          
            call col3(ta, vtrans(1,1,1,1,ifield),bm1, ntot)      
            call dssum(ta, nx1, ny1, nz1)                        
            call col2(ta, bintm1,ntot)                           
            call invcol3(ydot(j), bq(1,1,1,1,ifield-1), ta, ntot) 
         else                                                    
            call invcol3(ydot(j),bq(1,1,1,1,ifield-1),        
     &                   vtrans(1,1,1,1,ifield),ntot)           
         endif     
         j = j + ntot
      enddo


      j = 1
      if(ntf.eq.1) then
        do ifield = cv_1fld,cv_lfld
           if (ntflds(ifield-1).eq.0) call dssum (ydot(j),nx1,ny1,nz1)
           j = j + nxyz*nelfld(ifield)
        enddo
      else ! if no field is a tmesh, use vector dssum
        n    = cv_lfld-cv_1fld + 1
        call nvec_dssum (ydot,ntotv,n,gsh_fld(1))
      endif

      j = 1
      do ifield=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(ifield)
         if (ntflds(ifield-1).eq.0) call col2(ydot(j),binvm1,ntot)
         j = j + ntot
      enddo

      call add_fcvfun_usr(ydot,j)

      j = 1 ! enforce dirichlet bcs, no time varying bcs support for now
      do ifield=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(ifield)
         call col2(ydot(j),tmask(1,1,1,1,ifield-1),ntot)
         j = j + ntot
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
      subroutine cv_settime(cv_time)

      include 'SIZE'
      include 'TSTEP'
      include 'CVODE'

      cv_dtNek = cv_time - (time-dt)    
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
