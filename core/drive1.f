c-----------------------------------------------------------------------
      subroutine nek_mem_init
c
c     Allocate all runtime memory blocks (the module init routines that
c     replaced the old static common blocks). Lives here in the CORE
c     archive so BOTH the standalone program (drive.f) and library callers
c     (who call nek_init directly) get the allocations. Guarded so it runs
c     only once.
c
      use parallel_mod, only : parallel_init => init
      use scrns_mod, only : scrns_init => init
      use ctmp0_mod, only : ctmp0_init => init
      use ctmp1_mod, only : ctmp1_init => init
      use soln_mod, only : soln_init => init
      use orthov_mod, only : orthov_init => init
      use screv_mod, only : screv_init => init
      use scrcg_mod, only : scrcg_init => init
      use scrch_mod, only : scrch_init => init
      use scrmg_mod, only : scrmg_init => init
      use scrsf_mod, only : scrsf_init => init
      use scruz_mod, only : scruz_init => init
      use scrvh_mod, only : scrvh_init => init
      use noncon_mod, only : noncon_init => init
      use neknek_mod, only : neknek_init => init
      use mass_mod, only : mass_init => init
      use geom_mod, only : geom_init => init
      use dealias_mod, only : dealias_init => init
      use cvode_mod, only : cvode_init => init
      use adjoint_mod, only : adjoint_init => init
      use cvflow_a_mod, only : cvflow_a_init => init
      use scrdg_mod, only : scrdg_init => init
      use c_is1_mod, only : c_is1_init => init
      use fastd_mod, only : fastd_init => init
      use scrxxti_mod, only : scrxxti_init => init
      use fdmh1_mod, only : fdmh1_init => init
      use scrpre_mod, only : scrpre_init => init
      use swaplengths_mod, only : swaplengths_init => init
      use weightop_mod, only : weightop_init => init
      use scrhi_mod, only : scrhi_init => init
      use input_mod, only : input_init => init
      use topol_mod, only : topol_init => init
      use scrct_mod, only : scrct_init => init
      use ctmpf_mod, only : ctmpf_init => init
      use scrxxt_mod, only : scrxxt_init => init
      use scrpr2_mod, only : scrpr2_init => init
      use fastg_mod, only : fastg_init => init
      use ivrtx_mod, only : ivrtx_init => init
      use domain_mod, only : domain_init => init

      logical icalld
      save    icalld
      data    icalld /.false./

      if (icalld) return
      icalld = .true.

      call parallel_init()
      call scrns_init()
      call ctmp0_init()
      call ctmp1_init()
      call soln_init()
c     hsmg_mod: lazy in hsmg.f
c     mvgeom_mod: lazy in bdry.f/connect1.f/plan5.f/cvode_driver.f
c     vproj_mod: lazy in induct.f's ophinv
c     orthot_mod: lazy in conduct.f/perturb.f
      call orthov_init()
      call screv_init()
      call scrcg_init()
      call scrch_init()
      call scrmg_init()
      call scrsf_init()
      call scruz_init()
      call scrvh_init()
c     orthostrs_mod: lazy in subs1.f's hmhzsf
c     orthop_mod: lazy in plan4.f
      call noncon_init()
      call neknek_init() ! tiny -- just nfld_neknek, touched every run
c     the rest of neknek_mod: lazy in multimesh.f's neknek_setup
      call mass_init()
c     gmres_mod: lazy in gmres.f/hmholtz.f
      call geom_init()
      call dealias_init()
      call cvode_init()
c     avg_mod: lazy in navier5.f's avg_all
      call adjoint_init()
c     cbplan_vol_ms_mod: lazy in multimesh.f's plan3_vol_ms/plan4_vol_ms
      call cvflow_a_init()
      call scrdg_init()
c     outtmp_mod: lazy in prepost.f/postpro.f
c     orthox_mod: lazy in induct.f's setrhsp/gensolnp
c     cvflow_nn_mod: lazy in multimesh.f's plan3_vol_ms/plan4_vol_ms
      call c_is1_init()
      call fastd_init()
      call scrxxti_init()
      call fdmh1_init()
      call scrpre_init()
      call swaplengths_init()
      call weightop_init()
      call scrhi_init()
      call input_init()
      call topol_init()
      call scrct_init()
      call ctmpf_init()
      call scrxxt_init()
      call scrpr2_init()
      call fastg_init()
      call ivrtx_init()
      call domain_init()

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_init(comm)
      use c_is1_mod
      use ivrtx_mod
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
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

      integer comm
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
  
      common /rdump/ ntdump

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest 

      integer*8, pointer :: glo_num(:,:)
      integer*8, pointer :: vertex(:)
      integer*8 ngv

c     allocate all runtime memory blocks (once) before any cb_* is used.
      call nek_mem_init

      glo_num(1:lx1*ly1*lz1,1:lelt) => cb_c_is1(1:lx1*ly1*lz1*lelt)
      vertex(1:(2**ldim)*lelt) => cb_ivrtx(1:(2**ldim)*lelt)

      ! set word size for REAL
      wdsize = sizeof(rtest)
      ! set word size for INTEGER
      isize = sizeof(itest)
      ! set word size for INTEGER*8
      isize8 = sizeof(itest8) 
      ! set word size for LOGICAL
      lsize = sizeof(ltest) 
      ! set word size for CHARACTER
      csize = sizeof(ctest)

      call setupcomm(comm,newcomm,newcommg,'','')
      intracomm   = newcomm   ! within a session
      nekcomm     = newcomm
      iglobalcomm = newcommg  ! across all sessions
      call iniproc()

      if (nid.eq.nio) call printHeader

      etimes = dnekclock()
      istep  = 0

      call opcount(1)

      call initdim         ! Initialize / set default values.
      call initdat
      call files

      call readat          ! Read .rea +map file

      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1/lx2/lx3/lxd: ',lx1,lx2,lx3,lxd
 12      format(1X,A,4I12)
         write(6,*)
      endif

      call setvar          ! Initialize most variables

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology  

      call genwz           ! Compute GLL points, weights, etc.

      if(nio.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat' 

      call gengeom(igeom)  ! Generate geometry, after usrdat 

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      if(nio.eq.0) write(6,*) 'call usrdat2'
      do iref=1,nhref
         call h_refine_usrdat2(hrefcuts(iref))
         call fix_geom
      enddo

      call usrdat2
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat2' 

      call count_bdry   ! count the number of faces with assigned BCs
      call fix_geom

      call chk_axis         ! verify axisymmetric mesh/BC requirements
      call vrdsmsh          ! verify mesh topology
      call mesh_check(ifjac0_abort,2,0) ! check mesh and print metrics

      call setlog(.true.)   ! Initalize logical flags

      if (ifneknekc) call neknek_setup

      call bcmask  ! Set BC masks for Dirichlet boundaries.

      if (fintim.ne.0.0 .or. nsteps.ne.0) 
     $   call geneig(igeom) ! eigvals for tolerances

      call dg_setup ! Setup DG, if dg flag is set.

      if (ifflow.and.iftran) then ! Init pressure solver 
         if (fintim.ne.0 .or. nsteps.ne.0) call prinit
      endif

      if(ifcvode) call cv_setsize

      if(nio.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nio.eq.0) write(6,'(A,/)') ' done :: usrdat3'

      call setics
      call setprop

      if (instep.ne.0) then
         if (ifneknekc) call neknek_exchange
         if (ifneknekc) call chk_outflow

         if (nio.eq.0) write(6,*) 'call userchk'
         call userchk
         if(nio.eq.0) write(6,'(A,/)') ' done :: userchk' 
      endif
      call mesh_check(.true.,1,1) ! check mesh for possible changes from setics or userchk
      call setprop      ! call again because input has changed in userchk

      if (ifcvode .and. nsteps.gt.0) call cv_init

      call comment
      call sstest (isss) 

      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow
      p0thn = p0th

      call in_situ_init()

      call time00       !     Initalize timers to ZERO
      call opcount(2)

      ntdump=0
      if (timeio.ne.0.0) ntdump = int( time/timeio )

      tinit = dnekclock_sync() - etimes
      if (nio.eq.0) then
        write (6,*) ' '
        if (time.ne.0.0) write (6,'(a,e14.7)') ' Initial time:',time
        write (6,'(a,g13.5,a)') 
     &     ' Initialization successfully completed ', tinit, ' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'

      call nekgsync()

      if (instep.eq.0) then
        if(nid.eq.0) write(6,'(/,A,/,A,/)') 
     &     ' nsteps=0 -> skip time loop',
     &     ' running solver in post processing mode'
      else
        if(nio.eq.0) write(6,'(/,A,/)') 'Starting time loop ...'
      endif

      isyc  = 0
      if(ifsync) isyc=1
      itime = 0
#ifdef TIMER
      itime = 1
#endif

      ! start measurements
      dtmp = dnekgflops()

      istep  = 0
      msteps = 1

      irstat = int(param(120))

      do kstep=1,nsteps,msteps
         call nek__multi_advance(kstep,msteps)
         if(kstep.ge.nsteps) lastep = 1
         call check_ioinfo  
         call set_outfld
         etime1 = dnekclock()
         call userchk
         tuchk = tuchk + dnekclock()-etime1
         call prepost (ifoutfld,'his')
         call in_situ_check()
         if (mod(kstep,irstat).eq.0 .and. lastep.eq.0) call runstat 
         if (lastep .eq. 1) goto 1001
      enddo
 1001 lastep=1

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

      ntot = lx1*ly1*lz1*nelv

      call nekgsync

      call setup_convect(2) ! Save conv vel

      if (iftran) call settime
      if (ifmhd ) call cfl_check
      call setsolv
      call comment

      if (ifsplit) then   ! PN/PN formulation

         do igeom=1,ngeom

         if (ifneknekc .and. igeom.gt.2) then
            if (ifneknekm.and.igeom.eq.3) call neknek_setup
            call neknek_exchange
         endif

         ! call here before we overwrite wx 
         if (ifheat .and. ifcvode) call heat_cvode (igeom)   

         if (ifgeom) then
            call gengeom (igeom)
            call geneig  (igeom)
         endif

         if (ifheat) call heat (igeom)

         if (igeom.eq.2) then  
            call setprop
            call rzero(qtl,ntot)
            if (iflomach) call qthermal
         endif

         if (ifflow)          call fluid    (igeom)
         if (ifmvbd)          call meshv    (igeom)
         if (igeom.eq.ngeom.and.filterType.eq.1)
     $                        call q_filter(param(103))

         enddo

      else                ! PN-2/PN-2 formulation
         call setprop
         do igeom=1,ngeom

            if (ifneknekc .and. igeom.gt.2) then
              if (ifneknekm.and.igeom.eq.3) call neknek_setup
              call neknek_exchange
            endif

            ! call here before we overwrite wx 
            if (ifheat .and. ifcvode) call heat_cvode (igeom)   

            if (ifgeom) then
               if (.not.ifrich) call gengeom (igeom)
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
            if (igeom.eq.ngeom.and.filterType.eq.1)
     $         call q_filter(param(103))
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TOTAL'
      include 'DPROCMAP'
      include 'RESTART'
      include 'mpif.h'

      if(instep.ne.0) call runstat

c      if (ifstrs) then
c         call crs_free(xxth_strs) 
c      else
c         call crs_free(xxth(1))
c      endif

#ifdef DPROCMAP
#ifdef MPI
      call MPI_Win_free(dProcmapH, ierr)
#endif
#endif 
 
c     rsH (the restart RMA window) is now created+freed per restart inside mfi,
c     so there is nothing to free here. commrs (the communicator dup) is kept
c     for possible reuse; MPI_Finalize releases it.

      call in_situ_end()
      call exitt0()

      return
      end
c-----------------------------------------------------------------------
      subroutine nek__multi_advance(kstep,msteps)

      include 'SIZE'
      include 'TOTAL'

      do i=1,msteps
         istep = istep+i
         call nek_advance

         if (ifneknekc) then 
            call neknek_exchange
            call bcopy
            call chk_outflow
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
