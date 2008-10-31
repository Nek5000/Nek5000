      subroutine NEKTON
C
C       VERSION 2.7
C
C------------------------------------------------------------------------
C
C       3-D Isoparametric Legendre Spectral Element Solver.
C
C
C
C       N a v i e r     S t o k e s     S o l v e r 
C
C
C       LAPLACIAN FORMULATION
C
C       The program solves the equations 
C
C       (1) ro ( dv/dt + (v*grad)v ) = mu*div(grad(v)) - grad(p) + f 
C       (2) div(v) = 0,
C
C       subject to periodic, Dirichlet or Neumann boundary 
C       conmonditions for the velocity v.
C
C       STRESS FORMULATION
C
C       (1) ro*( dv/dt + (v*grad)v ) = mu*( div grad v + grad div v )
C                                      - grad p + f 
C       (2) div(v) = 0,
C
C       subject to velocity and traction boundary conditions.
C
C
C       In both formulations, the velocity v is solved on a 
C       Gauss-Legendre Lobatto spectral element mesh, while the
C       pressure is solved on a Gauss-Legendre mesh (staggered mesh).
C
C
C
C       P a s s i v e     S c a l a r     S o l v e r
C
C
C       The program solves the equation
C
C       (3)  rhocp ( dT/dt + v*grad T ) = k div(grad T) + q 
C
C       subject to Dirichlet and Neumann boundary conditions
C       for the passive scalar T.
C
C--------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
      include 'DEALIAS'
      include 'DOMAIN'
      include 'ZPER'
c
      include 'OPCTR'
      include 'CTIMER'

C     Declare scratch arrays
C     NOTE: no initial declaration needed. Linker will take 
c           care about the size of the CBs
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
  
      REAL e, oe
      integer WDS
      REAL*8 t0,tp

      call iniproc !  processor initialization 
      if (nid.eq.0) write(6,*) 'Number of Processors ::',np

      TIME0  = dnekclock()
      etimes = dnekclock()
      ISTEP  = 0
      tpp    = 0.0

      call opcount(1)

      call initdim
C
C     Data initialization
C
      call initdat
      call files
      t0 = dnekclock()

      call readat  ! Read processor map, followed by data.
      if (nid.eq.0) write(6,*) 'readat time ::',dnekclock()-t0,
     &                         ' seconds'

      call setvar  ! initialize some variables
      call echopar ! echo back the parameter stack

c     Check for zero steps

      instep=1
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0


C     Geometry initialization

      igeom = 2

      call connect
      call genwz
      call usrdat
      call gengeom (igeom)
      call usrdat2
      call geom_reset(1)    ! recompute Jacobians, etc.

      if (nid.eq.0) write(6,*) 'NELGV/NELGT/NX1:',nelgv,nelgt,nx1
      if (nid.eq.0) write(6,*) nid,' call vrdsms',instep

      call vrdsmsh  ! verify mesh topology


C     Field initialization
      if (nid.eq.0) write(6,*) nid,' call setlog',instep
      call setlog

      call bcmask
C
C     Need eigenvalues to set tolerances in presolve (SETICS)
C
      if (nid.eq.0) write(6,*) 'this is ifflow:',ifflow,nsteps
      if (fintim.ne.0.0.or.nsteps.ne.0) call geneig (igeom)
      call vrdsmsh
C
C     Solver initialization  (NOTE:  Uses "SOLN" space as scratch...)

      if (ifflow.and.nsteps.gt.0) then

         if (ifsplit) then
            call set_up_h1_crs
         else

            call estrat

            if (nid.eq.0) 
     $         write(6,*) nid,'estrat',iesolv,nsteps,fintim,solver_type

            if (fintim.ne.0.0.or.nsteps.ne.0) then
               if (iftran.and.solver_type.eq.'itr') then
                  call set_overlap
               elseif (solver_type.eq.'fdm'.or.solver_type.eq.'pdm')then
                  call gfdm_init
               elseif (solver_type.eq.'25D') then
                  call g25d_init
               endif
            endif

         endif
      endif


      call usrdat3

C     The properties are set if PRESOLVE is used in SETICS,
C     otherwise they are set in the beginning of the time stepping loop

      call setics

      CALL SETPROP

      CALL USERCHK
      CALL COMMENT
      CALL SSTEST (ISSS) 

      CALL TIME00
      CALL opcount(2)
      CALL dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow

      DO 1000 ISTEP=1,NSTEPS

         IF (IFTRAN) CALL SETTIME
         if (ifmhd ) call cfl_check

         CALL SETSOLV
         CALL COMMENT

         if (ifsplit) then

            if (ifheat)      call heat     (0)

                             call setprop
                             call qthermal

            if (ifflow)      call fluid    (0)

         else

           call setprop
           do igeom=1,2

              if (ifgeom) then
                 call gengeom (igeom)
                 call geneig  (igeom)
              endif

              if (ifmhd) then
                             call induct   (igeom)
                 if (ifheat) call heat     (igeom)
              else
                 if (ifheat)      call heat     (igeom)
                 if (ifflow)      call fluid    (igeom)
                 if (ifmvbd)      call meshv    (igeom)
              endif

              if (ifpert) then
                 if (ifflow) call fluidp   (igeom)
                 if (ifheat) call heatp    (igeom)
              endif

           enddo

         endif

         if (.not.ifmhd) then  ! (filter in induct.f for ifmhd)
             if (param(103).gt.0) alpha_filt=param(103)
             if (param(103).gt.0) call q_filter(alpha_filt)
         endif

         call prepost (.false.,'his')
         call userchk

         if (lastep .eq. 1) goto 1001
 1000 CONTINUE
 1001 CONTINUE
c
      call opcount(3)
      call timeout
C
C----------------------------------------------------------------------
C     Time stepping loop end
C----------------------------------------------------------------------
C
c     call prepost (.true.,'his')
c
      if (instep.eq.0) then
         lastep=1
         t0 = dnekclock()
         call prepost (.true.,'his')
         tpp = tpp + (dnekclock()-t0)
         nsteps=0
         call userchk
      endif
C
      if (nid.eq.0) then
         write(6,*) 'prepost time ::',tpp,' seconds'
      endif
C
      CALL COMMENT
      CALL DIAGNOS
      call crs_stats(xxth)
      call exitt

      return
      end
