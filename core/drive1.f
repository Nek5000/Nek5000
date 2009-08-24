      subroutine nek_init
C--------------------------------------------------------------------------

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
  
      REAL e, oe
      REAL*8 t0, tpp
      common /drive1f/ e, oe, t0, tpp
      
C     Initalize Nek (MPI stuff, word sizes, ...)
      call iniproc

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

C     Initialize some variables
      call setvar  

c     Check for zero steps
      instep=1
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

C     Geometry initialization
      igeom = 2
      call connect
      call genwz

C     Initalize io unit
      call io_init

C     USRDAT
      if(nid.eq.0) write(6,*) 'call usrdat'
      call usrdat
      if(nid.eq.0) then 
        write(6,*) 'done :: usrdat'
        write(6,*) ' '
      endif

C     Reinitalize geometry (in case it was changed in usrdat)
      call gengeom (igeom)

C     USRDAT2
      if(nid.eq.0) write(6,*) 'call usrdat2'
      call usrdat2
      if(nid.eq.0) then 
        write(6,*) 'done :: usrdat2'
        write(6,*) ' '
      endif

      call geom_reset(1)    ! recompute Jacobians, etc.
      call vrdsmsh          ! verify mesh topology

      call echopar ! echo back the parameter stack
      call setlog  ! Initalize logical flags

C     Zero out masks corresponding to Dirichlet boundary points.
      call bcmask

C     Need eigenvalues to set tolerances in presolve (SETICS)
      if (fintim.ne.0.0.or.nsteps.ne.0) call geneig (igeom)
      call vrdsmsh

C     Pressure solver initialization  (NOTE:  Uses "SOLN" space as scratch...)
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

C     USRDAT3
      if(nid.eq.0) write(6,*) 'call usrdat3'
      call usrdat3
      if(nid.eq.0) then 
        write(6,*) 'done :: usrdat3'
        write(6,*) ' '
      endif

C     Set initial conditions + compute field properties
      call setics
      call setprop

C     USRCHK
      if(nid.eq.0) write(6,*) 'call userchk'
      call userchk
      if(nid.eq.0) then 
        write(6,*) 'done :: userchk'
        write(6,*) ' '
      endif

      call comment
      call sstest (isss) 

C     Initalize timers to ZERO
      call time00
      call opcount(2)
      call dofcnt

      jp = 0  ! Set perturbation field count to 0 for baseline flow

      call gsync()

      etims0 = dnekclock()
      IF (NID.EQ.0) THEN
        WRITE (6,*) ' '
        IF (TIME.NE.0.0) WRITE (6,'(A,E14.7)') ' Initial time:',TIME
        WRITE (6,'(A,g13.5,A)') 
     &              ' Initialization successfully completed ',
     &              etims0-etimes, ' sec'
      ENDIF

      return
      end
C--------------------------------------------------------------------------
      subroutine nek_solve

      include 'SIZE'
      include 'TSTEP'
      include 'INPUT'
      include 'CTIMER'

      IF (NSTEPS.EQ.0) then
          if (nid.eq.0) then
             write(6,*) ' '
             write(6,*) '0 time steps -> skip time loop'
          endif
          RETURN
      ENDIF

      IF (NID.EQ.0) THEN
        WRITE (6,*) ' '
        WRITE (6,*) 'Starting time loop ...'
        WRITE (6,*) ' '
      ENDIF

      DO ISTEP=1,NSTEPS
         call nek_advance
         call userchk
         call prepost (.false.,'his')
         if (lastep .eq. 1) goto 1001
      ENDDO

 1001 RETURN
      END
C--------------------------------------------------------------------------
      subroutine nek_advance

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      common /cgeom/ igeom

      IF (IFTRAN) CALL SETTIME
      if (ifmhd ) call cfl_check
      CALL SETSOLV
      CALL COMMENT
         
      if (ifsplit) then   ! PN/PN formulation
         igeom = 1
         if (ifheat)      call heat     (igeom)
         call setprop
         call qthermal
         igeom = 1
         if (ifflow)      call fluid    (igeom)
                          call setup_convect (2)
      else                ! PN-2/PN-2 formulation
         call setprop
         do igeom=1,2

            if (ifgeom) then
               call gengeom (igeom)
               call geneig  (igeom)
            endif

            if (ifmhd) then
                                call induct   (igeom)
               if (ifheat)      call heat     (igeom)

            elseif (ifpert) then

               if (ifbase.and.ifheat)  call heat          (igeom)
               if (ifbase.and.ifflow)  call fluid         (igeom)
               if (ifflow)             call fluidp        (igeom)
               if (ifheat)             call heatp         (igeom)
                                       call setup_convect (igeom)
            else  ! std. nek case

               if (ifheat)             call heat          (igeom)
               if (ifflow)             call fluid         (igeom)
               if (ifmvbd)             call meshv         (igeom)
                                       call setup_convect (igeom)
            endif
         enddo
      endif
      
      if (.not.ifmhd) then      ! (filter in induct.f for ifmhd)
         if (param(103).gt.0) alpha_filt=param(103)
         if (param(103).gt.0) call q_filter(alpha_filt)
      endif

      return
      end
C--------------------------------------------------------------------------
      subroutine nek_end

      include 'SIZE'
      include 'TSTEP'
      include 'PARALLEL'
      include 'OPCTR'

c      write(6,*) nid,etimes

      if (nid.eq.0) then
         write(6,*) ' '
         write(6,*) 'call nek_end'
         write(6,*) ' '
      endif

c     check for zero steps
      if (instep.eq.0) then
         lastep=1
         call prepost (.true.,'his')
         nsteps=0
         call userchk
      else
#ifndef NOTIMER
         if(nid.eq.0) write(6,*) 'runtime statistics:'
         call opcount(3)
    
         call timeout
#endif
         CALL COMMENT
         CALL DIAGNOS
         if(xxth.gt.0) call crs_stats(xxth)
      endif
   
      return
      end
