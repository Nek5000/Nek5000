c-----------------------------------------------------------------------
      subroutine readat_new
C
C     Read in run parameters from .par file
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ZPER'
c
      logical ifbswap

      call setDefaultParam

      call open_bin_file(ifbswap) ! this will also read the header

      if(nid.eq.0) call par_read(ierr)
      call bcast(ierr,isize)
      if(ierr .ne. 0) call exitt

      call bcastParam
      call chkParam
      if (.not.ifgtp) call mapelpr  ! read .map file, est. gllnid, etc.

      call bin_rd1(ifbswap) 

      call nekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine setDefaultParam
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ZPER'
c
      call rzero(param,200)
      call rzero(uparam,20)

      param(1)  = 1    ! density
      param(16) = 1    ! helmholz for temperature and scalars
      param(21) = 1e-6 ! pressure tolerance
      param(22) = 1e-8 ! helmholz tolerance

      param(26) = 0.5  ! max Courant number
      param(27) = 2    ! 2nd order in time

      param(65) = 1    ! just one i/o node
      param(66) = 6    ! write in binary
      param(67) = 6    ! read in binary
      param(93) = 20   ! number of vectors for projection

      param(94) = 0    ! turn on projection for helmholz solves
      param(95) = 5    ! turn on projection for pressure solve
      param(99) = 4    ! dealising
c
      iftmsh(0) = .false. 
      iftmsh(1) = .false. 
      do i=1,ldimt
         iftmsh(1+i) = .false.
      enddo

      ifxyo = .true.
      ifvo  = .false.
      ifpo  = .false.
      ifto  = .false.   
      do i=1,ldimt1
         ifpsco(i) = .false.
      enddo 

      ifadvc(1) = .true.  
      do i=1,ldimt1
         ifadvc(i+1) = .true.  
      enddo 

      ifflow    = .false.
      ifheat    = .false.  
      iftran    = .true.   
      ifaxis    = .false.  
      ifstrs    = .false. 
      iflomach  = .false. 
      ifmvbd    = .false.
      ifchar    = .false.  
      ifmhd     = .false. 
      ifuservp  = .false.  
      ifcyclic  = .false.
      ifusermv  = .false.
      ifmgrid   = .false.
      ifessr    = .false.
      ifreguo   = .false.
      ifbase    = .true.   
      ifpert    = .false. 
      ifaziv    = .false. 
      ifmoab    = .false.  
      ifcvode   = .false.

      ifgtp     = .false.
      ifdg      = .false.
      ifsync    = .false.  
      ifanls    = .false.  
      ifcoup    = .false.  
      ifvcoup   = .false.  
      ifexplvis = .false.
      ifcons    = .false.   

      ifmodel   = .false. 
      ifkeps    = .false.
      ifschclob = .false. 

      ifcmt     = .false. 
      ifvisc    = .false. 
      iffltr    = .false.

      call izero(matype,16*ldimt1)
      call rzero(cpgrp ,48*ldimt1)

      call blank (hcode ,11*lhis)
      call izero (lochis, 4*lhis)

      call blank (initc,15*132)

      return
      end
c-----------------------------------------------------------------------
      subroutine par_read(ierr)
c
c     parse .par file and set run parameters
c
c     still missing:
c     - field specific solver for temp+scalars
c     - field specific tolerances for temp+scalars 
c     - abs/rel tolerance support for helmholz
c     - mhd support

      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ZPER'

      character*132 c_out,txt
 
      call finiparser_load(parfle,ierr)
      if(ierr .ne. 0) return

      call par_verify(ierr)
      if(ierr .ne. 0) return

c set parameters
      call finiparser_getString(c_out,'general:stopAt',ifnd)
      call capit(c_out,132)
      if (index(c_out,'ENDTIME') .gt. 0) then
         call finiparser_getDbl(d_out,'general:endTime',ifnd)
         if(ifnd .eq. 1) param(10) = d_out 
      endif

      call finiparser_getDbl(d_out,'general:numSteps',ifnd)
      if(ifnd .eq. 1) param(11) = d_out 

      call finiparser_getDbl(d_out,'general:dt',ifnd)
      if(ifnd .eq. 1) param(12) = d_out

      param(12) = -1*abs(param(12))
      call finiparser_getBool(i_out,'general:variableDt',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 1) param(12) = abs(param(12)) 

      d_out = param(15)
      call finiparser_getDbl(d_out,'general:writeInterval',ifnd)
      call finiparser_getString(c_out,'general:writeControl',ifnd)
      call capit(c_out,132)
      if (index(c_out,'RUNTIME') .gt. 0) then
         param(14) = d_out
      else
         param(14) = 0
         param(15) = d_out
      endif

      ! overrule abs/rel tol
      param(20) = -1
      call finiparser_getDbl(d_out,'temperature:residualTol',ifnd)
      if(ifnd .eq. 1) param(20) = d_out 
      call finiparser_getDbl(d_out,'pressure:residualTol',ifnd)
      if(ifnd .eq. 1) param(21) = d_out 
      call finiparser_getDbl(d_out,'velocity:residualTol',ifnd)
      if(ifnd .eq. 1) param(22) = d_out 

      j = 0
      do i = 1,99
         write(txt,"('scalar',i2.2)") i
         call finiparser_find(i_out,txt,ifnd)
         if (ifnd .eq. 1) then
            j = j + 1 
            ifpso(i) = .true.
         endif
      enddo
      param(23) = j 

      call finiparser_getDbl(d_out,'temperature:relativeTol',ifnd)
      if(ifnd .eq. 1) param(24) = d_out 
      call finiparser_getDbl(d_out,'temperature:absoluteTol',ifnd)
      if(ifnd .eq. 1) param(25) = d_out 

      call finiparser_getString(c_out,'temperature:solver',ifnd)
      call capit(c_out,132)
      if (index(c_out,'CVODE') .gt. 0) then
         param(16) = 2
      endif

      call finiparser_getDbl(d_out,'general:maxCFL',ifnd)
      if(ifnd .eq. 1) param(26) = d_out

      call finiparser_getDbl(d_out,'general:tOrder',ifnd)
      if(ifnd .eq. 1) param(27) = d_out 

      call finiparser_getDbl(d_out,'magnetic:viscosity',ifnd)
      if(ifnd .eq. 1) param(29) = d_out 
      if(param(29).lt.0.0) param(29) = -1.0/param(29)

      call finiparser_getDbl(d_out,'general:perturbationModes',ifnd)
      if(ifnd .eq. 1) param(31) = d_out 

      call finiparser_getString(c_out,'pressure:preconditioner',ifnd)
      call capit(c_out,132)
      if (index(c_out,'SCHWARZ') .gt. 0) param(43) = 1

      call finiparser_getBool(i_out,'general:write8Byte',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 1) param(63) = 1 

      call finiparser_getDbl(d_out,'general:writeNParallelFiles',ifnd)
      if(ifnd .eq. 1) param(65) = d_out 

      call finiparser_getBool(i_out,'velocity:residualProj',ifnd) ! for all helmholz solves
      if(ifnd .eq. 1) then
        param(94) = 0
        if(i_out .eq. 1) param(94) = 5 
      endif

      call finiparser_getBool(i_out,'pressure:residualProj',ifnd)
      if(ifnd .eq. 1) then
        param(95) = 0
        if(i_out .eq. 1) param(95) = 5 
      endif

      call finiparser_getBool(i_out,'general:dealiasing',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 0) param(99) = -1 

      call finiparser_getBool(i_out,'general:filtering',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 1) then
        call finiparser_getDbl(d_out,'general:filterWeight',ifnd)
        if(ifnd .eq. 1) param(103) = d_out 
        call finiparser_getDbl(d_out,'general:addFilterModes',ifnd)
        if(ifnd .eq. 1) param(101) = d_out 
      endif

      do i = 1,20
         call blank(txt,132)
         write(txt,"('general:userParam',i2.2)") i
         call finiparser_getDbl(d_out,txt,ifnd)
         if(ifnd .eq. 1) uparam(i) = d_out
      enddo

c set logical flags
      call finiparser_getString(c_out,'general:timeStepper',ifnd)
      call capit(c_out,132)

      if (index(c_out,'CHAR') .gt. 0) then
         ifchar = .true.
      else if (index(c_out,'STEADY') .gt. 0) then
         iftran = .false.
      endif

      call finiparser_find(i_out,'velocity',ifnd)
      if(ifnd .eq. 1) then
        ifflow = .true.
        ifvo   = .true.
        ifpo   = .true.
      endif

      call finiparser_find(i_out,'temperature',ifnd)
      if(ifnd .eq. 1) then
        ifheat = .true.
        ifto   = .true.
      endif

      call finiparser_getBool(i_out,'mesh:motion',ifnd)
      if(ifnd .eq. 1 .and. i_out .eq. 1) then
        ifmvbd = .true.
        call finiparser_getString(c_out,'mesh:meshVelocity',ifnd)
        call capit(c_out,132)
        if (index(c_out,'USER') .gt. 0) ifusermv = .true.
      endif

      call finiparser_getBool(i_out,'problemType:axiSymmetry',ifnd)
      if(ifnd .eq. 1) then
        ifaxis = .false.
        if(i_out .eq. 1) ifaxis = .true.
      endif

      call finiparser_getBool(i_out,'problemType:swirl',ifnd)
      if(ifnd .eq. 1) then
        ifaziv = .false.
        if(i_out .eq. 1) ifaziv = .true.
      endif

      call finiparser_getBool(i_out,'problemType:cyclicBoundaries',ifnd)
      if(ifnd .eq. 1) then
        ifcyclic = .false.
        if(i_out .eq. 1) ifcyclic = .true.
      endif

      call finiparser_getBool(i_out,'problemType:perturbations',ifnd)
      if(ifnd .eq. 1) then
        ifpert = .false.
        if(i_out .eq. 1) ifpert = .true.
      endif

      call finiparser_getBool(i_out,'problemType:solveBaseFlow',ifnd)
      if(ifnd .eq. 1) then
        ifbase = .false.
        if(i_out .eq. 1) ifbase = .true.
      endif

      call finiparser_getBool(i_out,'problemType:lowMachNumber',ifnd)
      if(ifnd .eq. 1) then
        iflomach = .false.
        if(i_out .eq. 1) iflomach = .true.
      endif
  
      call finiparser_getBool(i_out,
     &                        'problemType:stressFormulation',ifnd)
      if(ifnd .eq. 1) then
        ifstrs = .false.
        if(i_out .eq. 1) ifstrs = .true.
      endif

      call finiparser_getBool(i_out,'problemType:userProperties',ifnd)
      if(ifnd .eq. 1) then
        ifuservp = .false.
        if(i_out .eq. 1) ifuservp = .true.
      endif

c set advection
      call finiparser_getBool(i_out,'velocity:advection',ifnd)
      if(ifnd .eq. 1) then
        ifadvc(1) = .false.
        if(i_out .eq. 1) ifadvc(1) = .true.
      endif
 
      call finiparser_getBool(i_out,'temperature:advection',ifnd)
      if(ifnd .eq. 1) then
        ifadvc(2) = .false.
        if(i_out .eq. 1) ifadvc(2) = .true.
      endif

      do i = 1,ldimt-1
         write(txt,"('scalar',i2.2)") i
         call finiparser_getBool(i_out,txt // ':advection',ifnd)
         if(ifnd .eq. 1) then
           ifadvc(i+2) = .false.
           if(i_out .eq. 1) ifadvc(i+2) = .true.
         endif
      enddo

c set mesh-field mapping
conjugateHeatTransfer
      call finiparser_getBool(i_out,'temperature:conjugateHeatTransfer',
     &                        ifnd)
      if(ifnd .eq. 1) then
        iftmsh(2) = .false.
        if(i_out .eq. 1) iftmsh(2) = .true.
      endif

c set output flags
      call finiparser_getBool(i_out,'mesh:writeToFieldFile',ifnd) 
      if(ifnd .eq. 1) then
        ifxyo = .false.
        if(i_out .eq. 1) ifxyo = .true.
      endif

      call finiparser_getBool(i_out,'velocity:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifvo = .false.
        if(i_out .eq. 1) ifvo = .true.
      endif

      call finiparser_getBool(i_out,'pressure:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifpo = .false.
        if(i_out .eq. 1) ifpo = .true.
      endif

      call finiparser_getBool(i_out,'temperature:writeToFieldFile',ifnd)
      if(ifnd .eq. 1) then
        ifto = .false.
        if(i_out .eq. 1) ifto = .true.
      endif

      do i = 1,ldimt-1
         write(txt,"('scalar',i2.2)") i
         call finiparser_getBool(i_out,txt // ':writeToFieldFile',ifnd)
         if(ifnd .eq. 1) then
           ifpsco(i) = .false.
           if(i_out .eq. 1) ifpsco(i) = .true.
         endif
      enddo

c set properties
      call finiparser_getDbl(d_out,'velocity:viscosity',ifnd)
      if(ifnd .eq. 1) cpfld(1,1) = d_out 
      if (cpfld(1,1) .lt.0.0) cpfld(1,1)  = -1.0/cpfld(1,1)
      call finiparser_getDbl(d_out,'velocity:density',ifnd)
      if(ifnd .eq. 1) cpfld(1,2) = d_out 

      call finiparser_getDbl(d_out,'temperature:conductivity',ifnd)
      if(ifnd .eq. 1) cpfld(2,1) = d_out 
      if (cpfld(2,1) .lt.0.0) cpfld(2,1)  = -1.0/cpfld(2,1)
      call finiparser_getDbl(d_out,'temperature:rhoCp',ifnd)
      if(ifnd .eq. 1) cpfld(2,2) = d_out 

      do i = 1,ldimt-1
         write(txt,"('scalar',i2.2)") i
         call finiparser_getDbl(d_out, txt // ':conductivity',ifnd)
         if(ifnd .eq. 1) cpfld(2+i,1) = d_out 
         if(cpfld(2+i,1) .lt.0.0) cpfld(2+i,1)  = -1.0/cpfld(2+i,1)
         call finiparser_getDbl(d_out, txt // ':rhoCp',ifnd)
         if(ifnd .eq. 1) cpfld(2+i,2) = d_out 
      enddo

c set restart options
      call finiparser_findTokens('general:startFrom', ',' , ifnd)
      do i = 1,min(ifnd,15)
         call finiparser_getToken(initc(i),i)
         if(index(initc(i),'0') .eq. 1) call blank(initc(i),132)
      enddo

      call finiparser_dump()
      call finiparser_free()

      return
      end
c-----------------------------------------------------------------------
      subroutine bcastParam
C
C     Broadcast run parameters to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ZPER'

      call bcast(param , 200*wdsize)
      call bcast(uparam, 20*wdsize)

      call bcast(ifchar , lsize)
      call bcast(iftran  , lsize)
      call bcast(ifflow  , lsize)
      call bcast(ifheat  , lsize)     
      call bcast(iflomach, lsize)
      call bcast(ifstrs  , lsize)
      call bcast(ifmvbd  , lsize)
      call bcast(ifusermv, lsize)
      call bcast(ifaxis  , lsize)
      call bcast(ifcyclic, lsize)
      call bcast(ifmhd   , lsize)
      call bcast(ifuservp, lsize)
      call bcast(ifpert, lsize)
      call bcast(ifbase, lsize)
      call bcast(ifmoab, lsize)
      call bcast(ifaziv, lsize)

      call bcast(ifadvc, ldimt1*lsize)
      call bcast(iftmsh, (ldimt1+1)*lsize)

      call bcast(cpfld, 3*ldimt1*wdsize)

      call bcast(ifxyo , lsize)
      call bcast(ifvo  , lsize)
      call bcast(ifpo  , lsize)
      call bcast(ifto  , lsize)
      call bcast(ifpsco, ldimt1*lsize)

      call bcast(initc, 15*132*csize) 

c set some internals 
      if (ndim.eq.3) if3d=.true.
      if (ndim.ne.3) if3d=.false.
      if (ndim.lt.0) ifgtp = .true.     ! domain is a global tensor product
      if (ifsplit) ifmgrid   = .true.

      param(1) = cpfld(1,2)
      param(2) = cpfld(1,1)
      param(7) = cpfld(2,2)
      param(8) = cpfld(2,1)

      npscal=int(param(23)) ! number of scalar sections
      npscl1=npscal+1
      npscl2=npscal+2
 
      npert = abs(param(31)) 

      if (if3d) ifaxis = .false.
      if (ifaxis) param(99) = 3             ! For now, at least.

      ifldmhd = npscal + 3
      if (ifmhd) then
         ifessr = .true.
         ifchar = .false.
         npscl1 = npscl1 + 1
         cpfld(ifldmhd,1) = param(29)  ! magnetic viscosity
         cpfld(ifldmhd,2) = param( 1)  ! magnetic rho same as for fluid
      endif

      if (ifaxis.and..not.ifsplit) then ! Use standard Schwarz/PCG solver
         ifmgrid   = .false.
         param(42) = 1.00000  !  p042 0=gmres/1=pcg
         param(43) = 1.00000  !  p043 0=semg/1=schwarz
         param(44) = 1.00000  !  p044 0=E-based/1=A-based prec.
      endif

      if (.not.iftran) then
         if (ifflow.and.ifsplit) then
            iftran=.true. ! no steady Stokes support for Pn/Pn
         else
            param(11) = 1.0
            param(12) = 1.0
            param(19) = 0.0
         endif
      endif
c
c     Check here for global fast diagonalization method or z-homogeneity.
c     This is here because it influence the mesh read, which follows.
      nelx   = abs(param(116))   ! check for global tensor-product structure
      nely   = abs(param(117))
      nelz   = abs(param(118))
      n_o    = 0

      if (n_o.eq.0) then
         ifzper=.false.
         ifgfdm=.false.
         if (nelz.gt.0) ifzper=.true.
         if (nelx.gt.0) ifgfdm=.true.
         if (nelx.gt.0) ifzper=.false.
      endif

      if (ifdg) param(59)=1  ! No fast operator eval for DG

      return
      END
c-----------------------------------------------------------------------
      subroutine chkParam
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'PARALLEL'
      INCLUDE 'CTIMER'
      INCLUDE 'ZPER'
c
      neltmx=np*lelt
      nelvmx=np*lelv

      neltmx=min(neltmx,lelg)
      nelvmx=min(nelvmx,lelg)

      nelgt = iglmax(nelgt,1)
      nelgv = iglmax(nelgv,1)

      if (npscal+1.gt.ldimt) then
         if(nid.eq.0) then
           write(6,21) ldimt,npscal+1
   21      format(//,2x,'Error: Nek has been compiled'
     $             /,2x,'       for',i4,' scalars. This run'
     $             /,2x,'       requires that ldimt be set to',i4,'.')
         endif
         call exitt
      endif

      if (nelgt.gt.neltmx.or.nelgv.gt.nelvmx) then
         if (nid.eq.0) then
          lelt_needed = nelgt/np
          if (mod(nelgt,np).ne.0) lelt_needed = lelt_needed + 1 
          write(6,82) lelt,lelg,lelt_needed,np,nelgt
   82         format(//,2X,'ABORT: Problem size too large!'
     $         ,/,2X
     $         ,/,2X,'This solver has been compiled for:'
     $         ,/,2X,'   number of elements/proc  (lelt):',i12
     $         ,/,2X,'   total number of elements (lelg):',i12
     $         ,/,2X
     $         ,/,2X,'Recompile with the following SIZE  parameters:'
     $         ,/,2X,'   lelt >= ',i12,'  for np = ',i12
     $         ,/,2X,'   lelg >= ',i12,/)
c           write(6,*)'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
c           write(6,*)'help:',lelt,lelv,lelgv
         endif
         call exitt
      endif

      if(nelgt.gt.nelgt_max) then
        if(nid.eq.0) write(6,*)
     $               'ABORT: Total number of elements too large!',
     $               '       nel_max = ', nelgt_max 
        call exitt
      endif

      if (nelt.gt.lelt) then
        write(6,'(A,3I12)') 'ABORT: nelt>lelt!', nid, nelt, lelt
        call exitt
      endif

      if (ifmvbd) then
         if (lx1.ne.lx1m.or.ly1.ne.ly1m.or.lz1.ne.lz1m) 
     $    call exitti
     $    ('Mesh motion requires lx1m=lx1 etc. in SIZE . $',lx1m)
      endif

      IF(NDIM.NE.LDIM) THEN
         IF(NID.EQ.0) THEN
           WRITE(6,10) LDIM,NDIM
   10      FORMAT(//,2X,'ERROR: Nek has been compiled'
     $             /,2X,'       for spatial dimension equal to',I2,'.'
     $             /,2X,'       The data file has dimension',I2,'.')
         ENDIF
         call exitt
      ENDIF

      IF (NPSCL1.GT.LDIMT .AND. IFMHD) THEN
         if(nid.eq.0) then
           WRITE(6,22) LDIMT,NPSCL1
   22      FORMAT(/s,2X,'Error: Nek has been compiled'
     $             /,2X,'       for',I4,' scalars.  A MHD run'
     $             /,2X,'       requires that LDIMT be set to',I4,'.')
         endif
         call exitt
      ENDIF

      if (if3d) then
         if (ly1.ne.lx1.or.lz1.ne.lx1) then
            if (nid.eq.0) write(6,13) lx1,ly1,lz1
   13       format('ERROR: lx1,ly1,lz1:',3i5,' must be equal for 3D')
            call exitt
         endif
         if (ly2.ne.lx2.or.lz2.ne.lx2) then
            if (nid.eq.0) write(6,14) lx2,ly2,lz2
   14       format('ERROR: lx2,ly2,lz2:',3i5,' must be equal for 3D')
            call exitt
         endif
      else
         if (ly1.ne.lx1.or.lz1.ne.1) then
            if (nid.eq.0) write(6,12) lx1,ly1,lz1
   12       format('ERROR: ',3i5,' must have lx1=ly1; lz1=1, for 2D')
            call exitt
         endif
         if (ly2.ne.lx2.or.lz2.ne.1) then
            if (nid.eq.0) write(6,11) lx2,ly2,lz2
   11       format('ERROR: ',3i5,' must have lx2=ly2; lz2=1, for 2D')
            call exitt
         endif
      endif

      if (lgmres.lt.5 .and. param(42).eq.0) then
         if(nid.eq.0) write(6,*)
     $   'WARNING: lgmres might be too low!'
      endif


      if (ifsplit) then
         if (lx1.ne.lx2) then
            if (nid.eq.0) write(6,43) lx1,lx2
   43    format('ERROR: lx1,lx2:',2i4,' must be equal for IFSPLIT=T')
            call exitt
         endif
      else
         if (lx2.lt.lx1-2) then
            if (nid.eq.0) write(6,44) lx1,lx2
   44    format('ERROR: lx1,lx2:',2i4,' lx2 must be lx-2 for IFSPLIT=F')
           call exitt
         endif
      endif

      if (ifmoab .and..not. ifsplit) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: MOAB in Pn-Pn-2 is not supported'
         call exitt
      endif

      ktest = (lx1-lx1m) + (ly1-ly1m) + (lz1-lz1m)
      if (ifstrs.and.ktest.ne.0) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: Stress formulation requires lx1m=lx1, etc. in SIZE'
         call exitt
      endif

      if (ifgfdm.and.ifsplit) call exitti
     $  ('ERROR: FDM (p116>0) requires lx2=lx1-2 in SIZE$',lx2)

      if (ifgfdm.and.lfdm.eq.0) call exitti
     $  ('ERROR: FDM requires lfdm=1 in SIZE$',lfdm)

      if (ifsplit .and. ifmhd) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: MHD in Pn-Pn is not supported'
         call exitt
      endif

      if (ifmhd .and. lbx1.ne.lx1) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: For MHD, need lbx1=lx1, etc.; Change SIZE '
         call exitt
      endif

      if (ifpert .and. lpx1.ne.lx1) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: For Lyapunov, need lpx1=lx1, etc.; Change SIZE '
      endif

      if (iflomach .and. .not.ifsplit) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: For lowMach, need lx2=lx1, etc.; Change SIZE '
         call exitt
      endif

      if (iflomach .and. .not.ifheat) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT For lowMach, need to solve for temperature too!'
         call exitt
      endif

      if (ifchar .and. param(99).lt.0) then
        if (nid.eq.0) write(6,*) 
     &     'ABORT: Characteristic scheme needs dealiasing!'
        call exitt
      endif

      if (ifchar.and.(nelgv.ne.nelgt)) call exitti(
     $ 'ABORT: Characteristics not supported w/ conj. ht transfer$',1)

      if (param(99).gt.-1 .and. (lxd.lt.lx1 .or. lyd.lt.ly1 .or.
     &   lzd.lt.lz1)) then
         if(nid.eq.0) write(6,*)
     &   'ABORT: Dealiasing space too small; Check lxd,lyd,lzd in SIZE '
         call exitt
      endif

      if (ifmoab) then
         if(nid.eq.0) write(6,*) 
     $   'ABORT: MOAB is not supported!'
         call exitt
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine par_verify(ierr)

      INCLUDE 'PARDICT'

      character*132 val,key

      character*132 txt
      character*1   tx1(132)
      equivalence   (tx1,txt)

      ierr = 0

      call finiparser_getDictEntries(n)
      do i = 1,n
         call finiparser_getPair(key,val,i,ifnd)
         call capit(key,132)
         do j = 1,PARDICT_NKEYS ! do we find the key in the par-dictionary 
            if(index(pardictkey(j),key).eq.1) goto 10      

            is = index(key,'SCALAR')
            if(is .eq. 1) then
              call chcopy(txt,key,132)
              call chcopy(tx1(is+6),'%%',2) 
              if(index(pardictkey(j),txt).eq.1) goto 10
            endif
         enddo
         write(6,*) 'ERROR: Par file contains unknown key ', key
         ierr = ierr + 1
   10 enddo

      return
      end
