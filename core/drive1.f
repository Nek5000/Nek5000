c-----------------------------------------------------------------------
      subroutine nek_init(intracomm)

      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

C     used scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs automatically
c
c      COMMON /CTMP1/ DUMMY1(LCTMP1)
c      COMMON /CTMP0/ DUMMY0(LCTMP0)
c
c      COMMON /SCRNS/ DUMMY2(LX1,LY1,LZ1,LELT,7)
c      COMMON /SCRUZ/ DUMMY3(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCREV/ DUMMY4(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRVH/ DUMMY5(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRMG/ DUMMY6(LX1,LY1,LZ1,LELT,4)
c      COMMON /SCRCH/ DUMMY7(LX1,LY1,LZ1,LELT,2)
c      COMMON /SCRSF/ DUMMY8(LX1,LY1,LZ1,LELT,3)
c      COMMON /SCRCG/ DUMM10(LX1,LY1,LZ1,LELT,1)
  
      real kwave2
      real*8 t0, tpp

      logical ifemati,ifsync_

      call get_session_info(intracomm)
      nio = -1
      if(nid.eq.0) nio=0
      
C     Initalize Nek (MPI stuff, word sizes, ...)
c     call iniproc (initalized in get_session_info)

      etimes = dnekclock()
      istep  = 0
      tpp    = 0.0

      call opcount(1)

C     Initialize and set default values.
      call initdim
      call initdat
      call files

C     Read .rea +map file
      etime1 = dnekclock()
      call readat

      ifsync_ = ifsync
      ifsync = .true.

C     Initialize some variables
      call setvar  

c     Map BCs
      if (ifmoab) then
#ifdef MOAB
        call nekMOAB_bcs
#endif
      endif

c     Check for zero steps
      instep=1
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

C     Setup domain topology  
      igeom = 2
      call setup_topo

C     Compute GLL stuff (points, weights, derivate op, ...)
      call genwz

C     Initalize io unit
      call io_init

C     Set size for CVODE solver
      if(ifcvode .and. nsteps.gt.0) call cv_setsize(0,nfield)

C     USRDAT
      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

C     generate geometry (called after usrdat in case something changed)
      call gengeom (igeom)

      if (ifmvbd) call setup_mesh_dssum ! Set up dssum for mesh (needs geom)

C     USRDAT2
      if(nio.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call geom_reset(1)    ! recompute Jacobians, etc.
      call vrdsmsh          ! verify mesh topology

      call echopar ! echo back the parameter stack
      call setlog  ! Initalize logical flags

C     Zero out masks corresponding to Dirichlet boundary points.
      call bcmask

C     Need eigenvalues to set tolerances in presolve (SETICS)
      if (fintim.ne.0.0.or.nsteps.ne.0) call geneig (igeom)

C     Verify mesh topology
      call vrdsmsh

C     Pressure solver initialization (uses "SOLN" space as scratch)
      if (ifflow.and.(fintime.ne.0.or.nsteps.ne.0)) then
         call estrat
         if (iftran.and.solver_type.eq.'itr') then
            call set_overlap
         elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
            ifemati = .true.
            kwave2  = 0.0
            if (ifsplit) ifemati = .false.
            call gfdm_init(nx2,ny2,nz2,ifemati,kwave2)
         elseif (solver_type.eq.'25D') then
            call g25d_init
         endif
      endif

C     Initialize optional plugin
      call init_plugin

C     USRDAT3
      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

C     Set initial conditions + compute field properties
      call setics
      call setprop

C     USRCHK
      if(instep.ne.0) then
        if(nio.eq.0) write(6,*) 'call userchk'
        call userchk
        if(nio.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif

C     Initialize CVODE
      if(ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow

      call in_situ_init()

C     Initalize timers to ZERO
      call time00
      call opcount(2)

      etims0 = dnekclock_sync()
      IF (NIO.EQ.0) THEN
        WRITE (6,*) ' '
        IF (TIME.NE.0.0) WRITE (6,'(A,E14.7)') ' Initial time:',TIME
        WRITE (6,'(A,g13.5,A)') 
     &              ' Initialization successfully completed ',
     &              etims0-etimes, ' sec'
      ENDIF

      ifsync = ifsync_ ! restore initial value

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'

      real*4 papi_mflops
      integer*8 papi_flops

      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      itime = 0
      if(ifsync) isyc=1
      itime = 1
      call nek_comm_settings(isyc,itime)

      call nek_comm_startstat()

      istep  = 0
      msteps = 1

      do kstep=1,nsteps,msteps
         call nek__multi_advance(kstep,msteps)
         call userchk
         call prepost (.false.,'his')
         call in_situ_check()
         if (lastep .eq. 1) goto 1001
      enddo
 1001 lastep=1


      call nek_comm_settings(isyc,0)

      call comment

c     check for post-processing mode
      if (instep.eq.0) then
         nsteps=0
         istep=0
         if(nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,*) 'done :: userchk'
         call prepost (.true.,'his')
      else
         if (nio.eq.0) write(6,'(/,A,/)') 
     $      'end of time-step loop' 
      endif


      RETURN
      END
c-----------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      call nekgsync
      IF (IFTRAN) CALL SETTIME
      if (ifmhd ) call cfl_check
      CALL SETSOLV
      CALL COMMENT

      if (ifsplit) then   ! PN/PN formulation

         igeom = 1
         if (ifheat)          call heat     (igeom)
         call setprop
         call qthermal
         igeom = 1
         if (ifflow)          call fluid    (igeom)
         if (param(103).gt.0) call q_filter(param(103))
         call setup_convect (2) ! Save convective velocity _after_ filter

      else                ! PN-2/PN-2 formulation

         call setprop
         do igeom=1,ngeom

            if (igeom.gt.2) call userchk_set_xfer

            if (ifgeom) then
               call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
               if (ifheat)      call heat     (igeom)
                                call induct   (igeom)
            elseif (ifpert) then
               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
            else  ! std. nek case
               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
            endif

            if (igeom.eq.ngeom.and.param(103).gt.0) 
     $          call q_filter(param(103))

            call setup_convect (igeom) ! Save convective velocity _after_ filter

         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'OPCTR'

      if(instep.ne.0)  call runstat
      if(xxth(1).gt.0) call crs_stats(xxth(1))

   
      call in_situ_end()
      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance(kstep,msteps)

      include 'SIZE'
      include 'TOTAL'

      do i=1,msteps
         istep = istep+i
         call nek_advance
      enddo

      return
      end
c-----------------------------------------------------------------------
