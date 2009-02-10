c-----------------------------------------------------------------------
      subroutine setics
C-----------------------------------------------------------------------
C
C     Set initial conditions.
C
C-----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'DEALIAS'
      INCLUDE 'INPUT'
      INCLUDE 'IXYZ'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'MASS'
      INCLUDE 'MVGEOM'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
c
      logical  iffort(  ldimt1,0:lpert)
     $       , ifrest(0:ldimt1,0:lpert)
     $       , ifprsl(  ldimt1,0:lpert)
c
      LOGICAL  IFANYP
      common /rdump/ ntdump
      common /inelr/ nelrr
      common /ctmp1/ work(lx1,ly1,lz1,lelv)
     $ ,             ta1 (lx2,ly1,lz1)
     $ ,             ta2 (lx2,ly2,lz1)


      real psmax(LDIMT)

      if(nid.eq.0) then
        write(6,*) 'set initial conditions'
      endif
C
C     Initialize all fields:
C


      nxyz2=nx2*ny2*nz2
      ntot2=nxyz2*nelv
      nxyz1=nx1*ny1*nz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1
c
c
      CALL RZERO(VX,NTOTT)
      CALL RZERO(VY,NTOTT)
      CALL RZERO(VZ,NTOTT)
      CALL RZERO(PR,nxyz2*nelt)
      DO 10 IFLD=1,LDIMT
         CALL RZERO(T(1,1,1,1,IFLD),NTOTT)
   10 CONTINUE

      jp = 0                  ! set counter for perturbation analysis

      irst = param(46)        ! for lee's restart (rarely used)
      if (irst.gt.0) then
         if (param(99).eq.4) call set_convect_new(vxd,vyd,vzd)
         return
      endif


c     If moving geometry then add a perturbation to the
c     mesh coordinates (see Subroutine INIGEOM)

      IF (IFMVBD) CALL PTBGEOM      
C
C     Find out what type of i.c. is requested
C     Current options: 
C
C     (1) - User specified fortran function (default is zero i.c.)
C     (2) - Restart from file(s)
C     (3) - Activate pre-solver => steady diffusion / steady Stokes
C
C     If option (2) is requested, also return with the name of the
C     restart file(s) together with the associated dump number
C
      call slogic (iffort,ifrest,ifprsl,nfiles)
C
C     Set up proper initial values for turbulence model arrays
C
      IF (IFMODEL) CALL PRETMIC
C
C      ***** TEMPERATURE AND PASSIVE SCALARS ******
C
C     Check if any pre-solv necessary for temperature/passive scalars
C
      IFANYP = .FALSE.
      DO 100 IFIELD=2,NFIELD
         IF (IFPRSL(IFIELD,jp)) THEN
            IF (NID.EQ.0) WRITE(6,101) IFIELD
            IFANYP = .TRUE.
         ENDIF
  100 CONTINUE
  101 FORMAT(2X,'Using PRESOLVE option for field',I2,'.')
C
C
C     If any pre-solv, do pre-solv for all temperatur/passive scalar fields
C
      IF (IFANYP) CALL PRSOLVT
C
C     Fortran function initial conditions for temp/pass. scalars.
C
      MAXFLD = NFIELD
      IF (IFMODEL.AND.IFKEPS) MAXFLD = NFIELD-2
      if (ifmhd) maxfld = npscal+3
c
      jp = 0
      do 200 ifield=2,maxfld
         if (iffort(ifield,jp)) then
         if (nid.eq.0) write(6,*) 'call nekuic for ifld ', ifield
            call nekuic
         endif
 200  continue
c
      if (ifpert) then
         ifield=2
         do jp=1,npert
         if (nid.eq.0) write(6,*) 'nekuicP',ifield,jp,iffort(ifield,jp)
            if (iffort(ifield,jp)) call nekuic
         enddo
      endif
      jp = 0
C     
C
C     Restart files
C
      call restart(nfiles)
C
C
C      ***** VELOCITY ******
C
C     (If restarting for V, we're done,
C     ...else, do pre-solv for fluid if requested.)
C
      IFIELD = 1
      IF (IFPRSL(IFIELD,jp)) CALL PRSOLVV
C
C
C     Fortran function initial conditions for velocity.
C
      ifield = 1
      if (iffort(ifield,jp)) then
         if (nid.eq.0) write(6,*) 'call nekuic for vel  '
         call nekuic
      endif
c
      if (ifpert) then
         ifield=1
         do jp=1,npert
            if (iffort(ifield,jp)) call nekuic
            if (nid.eq.0) write(6,*) 'ic vel pert:',iffort(1,jp),jp
         enddo
      endif
      jp = 0
c
      ntotv = nx1*ny1*nz1*nelv
C
C     Fortran function initial conditions for turbulence k-e model
C
      if (ifmodel .and. ifkeps) then
         mfldt = nfield - 1
         do 300 ifield=mfldt,nfield
            if (iffort(ifield,jp)) call nekuic
 300     continue
      endif
C
C     Initial mesh velocities
C
      IF (IFMVBD) CALL OPCOPY (WX,WY,WZ,VX,VY,VZ)
      IF (IFMVBD.AND..NOT.IFREST(0,jp)) CALL MESHV (2)
C
C     Compute additional initial values for turbulence model arrays
C     based on I.C.
C
      IF (IFMODEL) CALL POSTMIC
C
C     If convection-diffusion of a passive scalar with a fixed velocity field,
C     make sure to fill up lagged arrays since this will not be done in
C     the time-stepping procedure (no flow calculation) (01/18/91 -EMR).
C
      IF (.NOT.IFFLOW.AND.IFHEAT) THEN
         ITEST=0
         DO 400 IFIELD=2,NFIELD
            IF (IFADVC(IFIELD)) ITEST=1
 400     CONTINUE
         IF (ITEST.EQ.1) THEN
            NBDMAX = 3
            NBDSAV = NBDINP
            NBDINP = NBDMAX
            DO 500 I=1,NBDMAX
               CALL LAGVEL
 500        CONTINUE
            NBDINP = NBDSAV
         ENDIF
      ENDIF
C     
C     Ensure that all processors have the same time as node 0.
C
      IF (NID.NE.0) TIME=0.0
      TIME=GLSUM(TIME,1)
      NTDUMP=0
      IF (TIMEIO.NE.0.0) NTDUMP = INT( TIME/TIMEIO )
C
C     Ensure that initial field is continuous!
C
      NXYZ1=NX1*NY1*NZ1
      NTOTT=NELT*NXYZ1
      NTOTV=NELV*NXYZ1
      NTOTG=NELGV*NXYZ1
c
C     first.. a test...
      ifield = 2
      if (ifflow) ifield = 1
      call rone(work,ntotv)
      ifield = 1
      CALL DSSUM(work,NX1,NY1,NZ1)
      CALL COL2(work,VMULT,NTOTV)
      rtot  = glsum(work,ntotv)
      rtotv = ntotg
      rdif  = (rtot-rtotv)/rtotv
c
      if (rdif.gt.0.0) then
         if (nid.eq.0) write(*,*) 'Abort: dssum has failed!'
         call exitt
      endif

      vxmax = glamax(vx,ntotv)
      vymax = glamax(vy,ntotv)
      vzmax = glamax(vz,ntotv)
      prmax = glamax(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glamax(t ,ntot)

      do i=1,NPSCAL
         ntot = nx1*ny1*nz1*nelfld(i+2)
         psmax(i) = glamax(T(1,1,1,1,i+1),ntot)
      enddo

c 
      small=1.0E-20
      ifldsave = ifield
      if (vxmax.eq.0.0) call perturb(vx,1,small)
      if (vymax.eq.0.0) call perturb(vy,1,small)
      if (vzmax.eq.0.0) call perturb(vz,1,small)
      if (prmax.eq.0.0) call perturb(pr,1,small)
      if (ttmax.eq.0.0) call perturb(t ,2,small)
c
      do i=1,NPSCAL
         ntot = nxyz1*nelfld(i+2)
         if(psmax(i).eq.0) call perturb(t(1,1,1,1,1+i),i+2,small)
      enddo
      ifield = ifldsave
    
      if(ifflow) then
         ifield = 1
         call opdssum(vx,vy,vz)
         call col2 (vx,vmult,ntotv)
         call col2 (vy,vmult,ntotv)
         call col2 (vz,vmult,ntotv)
         if (ifsplit) call dsavg(pr)  ! continuous pressure
      endif

      if (ifmhd) then
         ifield = ifldmhd
         call opdssum(bx,by,bz)
         call col2 (bx,vmult,ntotv)
         call col2 (by,vmult,ntotv)
         call col2 (bz,vmult,ntotv)
      endif
c
      if (ifsplit) then
         call dssum(pr,nx1,ny1,nz1)
         call col2 (pr,vmult,ntotv)
      endif
c
      if (ifheat) then
         ifield = 2
         call dssum(t ,nx1,ny1,nz1)
         call col2 (t ,tmult,ntott)
         do ifield=3,nfield
            call dssum(t(1,1,1,1,i-1),nx1,ny1,nz1)
            if(iftmsh(ifield)) then
              call col2 (t(1,1,1,1,i-1),tmult,ntott)
            else
              call col2 (t(1,1,1,1,i-1),vmult,ntotv)
            endif
         enddo
      endif
c
      if (ifpert) then
         do jp=1,npert
            ifield = 1
            call opdssum(vxp(1,jp),vyp(1,jp),vzp(1,jp))
            call col2 (vxp(1,jp),vmult,ntotv)
            call col2 (vyp(1,jp),vmult,ntotv)
            call col2 (vzp(1,jp),vmult,ntotv)
            ifield = 2
            call dssum(tp(1,1,jp),nx1,ny1,nz1)
            call col2 (tp(1,1,jp),tmult,ntotv)
c           note... must be updated for addl pass. scal's. pff 4/26/04
            vxmax = glamax(vxp(1,jp),ntotv)
            vymax = glamax(vyp(1,jp),ntotv)
            if (nid.eq.0) write(6,111) jp,vxmax,vymax
  111       format(i5,1p2e12.4,' max pert vel')
         enddo
      endif
      jp = 0

C print min values
      xxmax = glmin(xm1,ntott)
      yymax = glmin(ym1,ntott)
      zzmax = glmin(zm1,ntott)

      vxmax = glmin(vx,ntotv)
      vymax = glmin(vy,ntotv)
      vzmax = glmin(vz,ntotv)
      prmax = glmin(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glmin(t ,ntott)

      do i=1,NPSCAL
         ntot = nxyz1*nelfld(i+2)
         psmax(i) = glmin(T(1,1,1,1,i+1),ntot)
      enddo

      if (nid.eq.0) then
         write(6,19) xxmax,yymax,zzmax
   19    format(' xyz min  ',5g13.5)
      endif
      if (nid.eq.0) then
         write(6,20) vxmax,vymax,vzmax,prmax,ttmax
   20    format(' uvwpt min',5g13.5)
      endif
      if (npscal.gt.0) then
         if (nid.eq.0) write(6,21) (psmax(i),i=1,NPSCAL)
   21    format(' PS min   ',50g13.5)
      endif

c print max values
      xxmax = glmax(xm1,ntott)
      yymax = glmax(ym1,ntott)
      zzmax = glmax(zm1,ntott)

      vxmax = glmax(vx,ntotv)
      vymax = glmax(vy,ntotv)
      vzmax = glmax(vz,ntotv)
      prmax = glmax(pr,ntot2)

      ntot = nxyz1*nelfld(2)
      ttmax = glmax(t ,ntott)

      do i=1,NPSCAL
         ntot = nxyz1*nelfld(i+2)
         psmax(i) = glmax(T(1,1,1,1,i+1),ntot)
      enddo

      if (nid.eq.0) then
         write(6,16) xxmax,yymax,zzmax
   16    format(' xyz max  ',5g13.5)
      endif

      if (nid.eq.0) then
         write(6,17) vxmax,vymax,vzmax,prmax,ttmax
   17    format(' uvwpt max',5g13.5)
      endif

      if (npscal.gt.0) then
         if (nid.eq.0)  then
            write(6,18) (psmax(i),i=1,NPSCAL)
   18       format(' PS max   ',50g13.5)
         endif
      endif


      if (ifrest(0,jp)) then !  mesh has been read in.
         if (nid.eq.0) write(6,*) 'Restart: recompute geom. factors.'
         call geom_reset(1)  !  recompute geometric factors
      endif

c     ! save velocity on fine mesh for dealiasing
      if (param(99).eq.4) call set_convect_new(vxd,vyd,vzd)

      if(nid.eq.0) then
        write(6,*) 'done :: set initial conditions'
        write(6,*) ' '
      endif

      return
      end
C            
c-----------------------------------------------------------------------
      subroutine slogic (iffort,ifrest,ifprsl,nfiles)
C---------------------------------------------------------------------
C
C     Set up logicals for initial conditions.
C
C---------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
c
      logical  iffort(  ldimt1,0:lpert)
     $       , ifrest(0:ldimt1,0:lpert)
     $       , ifprsl(  ldimt1,0:lpert)
c
      character*80 line,fname,cdum
      character*2  s2
      character*1  line1(80)
      equivalence (line1,line)
C
C     Default is user specified fortran function (=0 if not specified)
C
      nfldt = nfield
      if (ifmhd) nfldt = nfield+1

      do jp=0,npert
         ifrest(0,jp) = .false.
         do ifld=1,nfldt
            iffort(ifld,jp) = .true.
            ifrest(ifld,jp) = .false.
            ifprsl(ifld,jp) = .false.
         enddo
      enddo

      jp = 0
      nfiles=0
C
C     Check for Presolve options     
C
      DO 1000 ILINE=1,15 
         LINE=INITC(ILINE)
         CALL CAPIT(LINE,80)
         IF (INDX1(LINE,'PRESOLV',7).NE.0) THEN
C           found a presolve request
            CALL BLANK(INITC(ILINE),80)
            CALL LJUST(LINE)
            CALL CSPLIT(CDUM,LINE,' ',1)
C
            IF (LTRUNC(LINE,80).EQ.0) THEN
               IF (NID.EQ.0) WRITE(6,700)
  700          FORMAT(/,2X,'Presolve options: ALL')
C              default - all fields are presolved.
               DO 800 IFIELD=1,nfldt
                  ifprsl(ifield,jp) = .true.
                  iffort(ifield,jp) = .false.
  800          CONTINUE
            ELSE
C           check line for arguments
C
               LL=LTRUNC(LINE,80)
               IF (NID.EQ.0) WRITE(6,810) (LINE1(L),L=1,LL)
  810          FORMAT(/,2X,'Presolve options: ',80A1)
C
               IF (INDX2(LINE,'U',1).NE.0) THEN
                  ifprsl(1,jp) = .true.
                  iffort(1,jp) = .false.
               ENDIF
C
               IF (INDX2(LINE,'T',1).NE.0) THEN
                  ifprsl(2,jp) = .true.
                  iffort(2,jp) = .false.
               ENDIF
C
               DO 900 IFIELD=3,NPSCAL+2
                  IP=IFIELD-2
                  WRITE(S2,901) IP
                  IF (INDX2(LINE,S2,2).NE.0) THEN
                     ifprsl(ifield,jp) = .true.
                     iffort(ifield,jp) = .false.
                  ENDIF
  900          CONTINUE
  901          FORMAT('P',I1)
            ENDIF
         ENDIF
 1000    CONTINUE
C
C     Check for restart options
C
      jp = 0
      DO 2000 ILINE=1,15
         if (ifpert) jp=iline-1
         LINE=INITC(ILINE)
         IF (LTRUNC(LINE,80).NE.0) THEN
C           found a filename
            NFILES=NFILES+1
            INITC(NFILES)=LINE
C
            IF (NID.EQ.0.AND.NFILES.EQ.1) WRITE(6,1010) LINE
 1010       FORMAT(1X,'Checking restart options: ',A80)
c            IF (NID.EQ.0) WRITE(6,'(A80)') LINE
C
C           Parse restart options
 
            call sioflag(ndumps,fname,line)

            IF (IFGETX) THEN
               IFREST(0,jp) = .TRUE.
            ENDIF
            IF (IFGETU) THEN
               iffort(1,jp) = .false.
               ifprsl(1,jp) = .false.
               ifrest(1,jp) = .true.
            ENDIF
            IF (IFGETT) THEN
               iffort(2,jp) = .false.
               ifprsl(2,jp) = .false.
               ifrest(2,jp) = .true.
            ENDIF
            DO 1900 IFIELD=3,nfldt
c              write(6,*) 'ifgetps:',(ifgtps(k),k=1,nfield)
               IF (IFGTPS(IFIELD-2)) THEN
                  iffort(ifield,jp) = .false.
                  ifprsl(ifield,jp) = .false.
                  ifrest(ifield,jp) = .true.
               ENDIF
 1900       CONTINUE
         ENDIF
 2000 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine restart(nfiles)
C----------------------------------------------------------------------
C
C     (1) Open restart file(s)
C     (2) Check previous spatial discretization 
C     (3) Map (K1,N1) => (K2,N2) if necessary
C
C     nfiles > 1 has several implications:
C
C     i.   For std. run, data is taken from last file in list, unless
C          explicitly specified in argument list of filename
C
C     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
C          subsequent files are for B-field or perturbation fields
C
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'RESTART'

      common /inelr/ nelrr

      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZT=LX1*LY1*LZ1*LELT)
      PARAMETER (LPSC9=LDIMT+9)
c
      COMMON /SCRNS/ SDUMP(LXYZT,7)
      integer mesg(40)
c
C     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 9.
      COMMON /CTMP1/ TDUMP(LXYZR,LPSC9)
      real*4         tdump
c
      REAL SDMP2(LXYZT,LDIMT)
c
c     cdump comes in via PARALLEL (->TOTAL)
c     COMMON /CBLELG/ CDUMP(LELG)
c
      CHARACTER*30 EXCODER
      CHARACTER*1  EXCODER1(30)
      EQUIVALENCE (EXCODER,EXCODER1)


      CHARACTER*80 FNAME
      CHARACTER*1  FNAME1(80)
      EQUIVALENCE (FNAME1,FNAME)
C
      INTEGER      HNAMI (30)
      CHARACTER*80 HNAME
      CHARACTER*1  HNAME1(80)
      EQUIVALENCE (HNAME,HNAME1)
      EQUIVALENCE (HNAME,HNAMI )

      CHARACTER*80 header
C
C     Local logical flags to determine whether to copy data or not.
C
      logical ifok,iffmat
      integer iposx,iposz,iposu,iposw,iposp,ipost,ipsps(ldimt1)
C
      logical ifbytsw, if_byte_swap_test
      real*4   bytetest
c
      REAL AXISM1 (LX1,LY1)
      REAL AXISM2 (LX2,LY2)
c
      ifok=.false.
      ifbytsw = .false.

      if(nfiles.lt.1) return
      if(nid.eq.0) write(6,*) 'Reading restart data'

      if (param(67).lt.1.0) then
         ! ascii 
         iffmat=.true.
      else
         ! binary
         iffmat=.false.
      endif

      if (param(67).eq.6.0) then
         do ifile=1,nfiles
            call sioflag(ndumps,fname,initc(ifile))
            call  mlti_file_input(fname)
         enddo
         return
      endif
c
      DO 6000 IFILE=1,NFILES
C
C        Loop over the requested restart files
C
         call sioflag(ndumps,fname,initc(ifile))
C
C
C       Standard Case:   Requested inputs are in the file.
C
        if (nid.eq.0) then

          if (iffmat) then
            open (unit=91,file=fname,status='old',err=500)
          else
            len= ltrunc(fname,79)
            call izero (hnami,20)
            call chcopy(hname1,fname,len)
c           test for presence of file
            open (unit=91,file=hname
     $           ,form='unformatted',status='old',err=500)
            close(UNIT=91)
            call byte_open(hname)
          ENDIF
          ifok = .true.
        endif

  500   continue
        call lbcast(ifok)
        if (.not.ifok) goto 5000
         
c
         ndumps = 1
C
C        Only NODE 0 reads from the disk.
C
         DO 1000 IDUMP=1,NDUMPS
C
            IF (NID.EQ.0) THEN
                ! read header
               if (iffmat) then
                 if(mod(param(67),1.0).eq.0) then ! old header format
                   if(nelgt.lt.10000) then
                     read(91,91,err=1500,end=1500)
     $              neltr,nxr,nyr,nzr,rstime,istepr,(excoder1(i),i=1,30)
   91                format(4i4,1x,g13.4,i5,1x,30a1)
                   else
                     read(91,92,err=1500,end=1500)
     $              neltr,nxr,nyr,nzr,rstime,istepr,(excoder1(i),i=1,30)
   92                format(i10,3i4,1P1e18.9,i9,1x,30a1)
                   endif
                 else                          ! new head format
                   read(91,'(A80)',err=1500,end=1500) header
                   read(header,*)
     &                  neltr,nxr,nyr,nzr,rstime,istepr,excoder
                 endif
               else
                 if(mod(param(67),1.0).eq.0) then  ! old header format
                   call byte_read(hnami,20)
                   icase = 2
                   if (nelgt.lt.10000) icase = 1
                   ipass = 0
   93              continue  ! test each possible case  UGLY (7/31/07)
                   if(ipass.lt.2) then
                     ipass = ipass+1
                     if(icase.eq.1) then
                       read(hname,'(4i4,1x,g13.4,i5,1x,30a1)',
     $                     err=94,end=94) 
     $                     neltr,nxr,nyr,nzr,rstime,istepr,
     $                     (excoder1(i),i=1,30)
                       goto 95
                     else
                       read(hname,'(i10,3i4,1P1e18.9,i9,1x,30a1)',
     $                     err=94,end=94)
     $                     neltr,nxr,nyr,nzr,rstime,istepr,
     $                     (excoder1(i),i=1,30)
                       goto 95
                     endif
   94                icase = 3-icase  !  toggle: 2-->1  1-->2
                     goto 93
                   else
                     goto 1500 ! failed
                   endif
   95              continue
                 else                         ! new head format
                   call byte_read(header,20)
                   read(header,*)
     &             neltr,nxr,nyr,nzr,rstime,istepr,excoder
                  endif
                  call byte_read(bytetest,1)
                  ifbytsw = if_byte_swap_test(bytetest)
               endif
               mesg(1) = neltr
               mesg(2) = nxr
               mesg(3) = nyr
               mesg(4) = nzr
               write(6,*)  'Read mode: ', param(67)
               write(6,333)'neltr,nxr,nyr,nzr: ', neltr,nxr,nyr,nzr
  333          format(A,i8,3i4)
               call chcopy(mesg(5),excoder1,20)
               len  = 14*isize
            endif
c
            IF (IDUMP.EQ.1) THEN
               len  = 14*isize
               call bcast(mesg,len)
               neltr = mesg(1)
               nxr   = mesg(2)
               nyr   = mesg(3)
               nzr   = mesg(4)
               call   chcopy(excoder1,mesg(5),20)
c
               call lbcast(ifbytsw)
C
C              Bounds checking on mapped data.
               IF (NXR.GT.LXR) THEN
                  WRITE(6,20) NXR,NX1
   20             FORMAT(//,2X,
     $            'ABORT:  Attempt to map from',I3,
     $            ' to N=',I3,'.',/,2X,
     $            'NEK5000 currently supports mapping from N+6 or less.'
     $            ,/,2X,'Increase N or LXR in IC.FOR.')
                  CALL EMERXIT
               ENDIF
C
C
C              Figure out position of data in file "IFILE"
C
               NOUTS=0
               IPOSX=0
               IPOSY=0
               IPOSZ=0
               IPOSU=0
               IPOSV=0
               IPOSW=0
               IPOSP=0
               IPOST=0
               DO 40 I=1,NPSCAL
                  IPSPS(I)=0
   40          CONTINUE

               IPS = 0
               NPS = 0
               DO 50 I=1, 30
                  IF (EXCODER1(I).EQ.'X') THEN
                     NOUTS=NOUTS + 1
                     IPOSX=NOUTS
                     NOUTS=NOUTS+1
                     IPOSY=NOUTS
                     IF (IF3D) THEN
                        NOUTS=NOUTS + 1
                        IPOSZ=NOUTS
                     ENDIF
                  ENDIF
                  IF (EXCODER1(I).EQ.'U') THEN
                     NOUTS=NOUTS + 1
                     IPOSU=NOUTS
                     NOUTS=NOUTS+1
                     IPOSV=NOUTS
                     IF (IF3D) THEN
                        NOUTS=NOUTS + 1
                        IPOSW=NOUTS
                     ENDIF
                  ENDIF
                  IF (EXCODER1(I).EQ.'P') THEN
                     NOUTS=NOUTS + 1
                     IPOSP=NOUTS
                  ENDIF
                  IF (EXCODER1(I).EQ.'T') THEN
                     NOUTS=NOUTS + 1
                     IPOST=NOUTS
                  ENDIF
                  IF(mod(param(67),1.0).eq.0.0) THEN
                    IF (EXCODER1(I).EQ.'1') THEN
                       NOUTS=NOUTS + 1
                       IPSPS(1)=NOUTS
                    ENDIF
                    IF (EXCODER1(I).EQ.'2') THEN
                       NOUTS=NOUTS + 1
                       IPSPS(2)=NOUTS
                    ENDIF
                    IF (EXCODER1(I).EQ.'3') THEN
                       NOUTS=NOUTS + 1
                       IPSPS(3)=NOUTS
                    ENDIF
                    IF (EXCODER1(I).EQ.'4') THEN
                       NOUTS=NOUTS + 1
                       IPSPS(4)=NOUTS
                    ENDIF
                  ELSE
                    IF(EXCODER1(I).EQ.'S') THEN
                       READ(EXCODER1(I+1),'(I1)') NPS1
                       READ(EXCODER1(I+2),'(I1)') NPS0
                       NPS=10*NPS1 + NPS0 
                       DO IS = 1, NPS
                         NOUTS=NOUTS + 1
                         IPSPS(IS)=NOUTS
                       ENDDO
                       GOTO 50
                    ENDIF
                  ENDIF
   50          CONTINUE
 
               IF (NPS.GT.(LDIMT-1)) THEN
                  IF (NID.EQ.0) THEN 
                    WRITE(*,'(A)') 
     &               'ERROR: restart file has a NSPCAL > LDIMT'
                    WRITE(*,'(A,I2)') 
     &               'Change LDIMT in SIZE'
                  ENDIF
                  CALL EXITT
               ENDIF

               LNAME=LTRUNC(FNAME,80)
               IF (NID.EQ.0) WRITE(6,61) (FNAME1(I),I=1,LNAME)
               IF (NID.EQ.0) WRITE(6,62) 
     $             IPOSU,IPOSV,IPOSW,IPOSP,IPOST,NPS,NOUTS
   61          FORMAT(/,2X,'Restarting from file ',80A1)
   62          FORMAT(2X,'Columns for restart data U,V,W,P,T,S,N: ',7I4)
C
C              Make sure the requested data is present in this file....
               IF (IPOSX.EQ.0) IFGETX=.FALSE.
               IF (IPOSY.EQ.0) IFGETX=.FALSE.
               IF (IPOSZ.EQ.0) IFGETZ=.FALSE.
               IF (IPOSU.EQ.0) IFGETU=.FALSE.
               IF (IPOSV.EQ.0) IFGETU=.FALSE.
               IF (IPOSW.EQ.0) IFGETW=.FALSE.
               IF (IPOSP.EQ.0) IFGETP=.FALSE.
               IF (IPOST.EQ.0) IFGETT=.FALSE.
               DO 65 I=2,NPSCAL
                  IF (IPSPS(I).EQ.0) IFGTPS(I)=.FALSE.
   65          CONTINUE
C
C              End of restart file header evaluation.
C
            ENDIF
C
C           Read the error estimators
C           not supported at the moment => just do dummy reading
C
            IF(NID.EQ.0)THEN
               if (iffmat)
     &             READ(91,'(6G11.4)',END=1500)(CDUMP,I=1,NELTR)
            ENDIF
C
C           Read the current dump, double buffer so that we can
C           fit the data on a distributed memory machine,
C           and so we won't have to read the restart file twice
C           in case of an incomplete data file.
C
            NXYZR = NXR*NYR*NZR
C
C           Read the data
C
            nelrr = min(nelgt,neltr) ! # of elements to _really_read
                                     ! why not just neltr?
            do 200 ieg=1,nelrr
               ifok = .false.
               IF (NID.EQ.0) THEN
                 IF (MOD(IEG,100).EQ.1) WRITE(6,*) 'Reading',IEG
                 IF (iffmat) THEN
                    READ(91,*,ERR=1500,END=1500)
     $              ((tdump(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
                 ELSE
                    do ii=1,nouts
                       call byte_read(tdump(1,II),nxyzr)
                    enddo
                 ENDIF
                 IFOK=.TRUE.
               ENDIF
C
C              Notify other processors that we've read the data OK.
C
               call lbcast(ifok)
               IF (.NOT.IFOK) GOTO 1600
C
C              MAPDMP maps data from NXR to NX1
C              (and sends data to the appropriate processor.)
C
C              The buffer SDUMP is used so that if an incomplete dump
C              file is found (e.g. due to UNIX io buffering!), then
C              the previous read data stored in VX,VY,.., is not corrupted.
C
               IF (IFGETX) CALL MAPDMP
     $         (SDUMP(1,1),TDUMP(1,IPOSX),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETX) CALL MAPDMP
     $         (SDUMP(1,2),TDUMP(1,IPOSY),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETZ) CALL MAPDMP
     $         (SDUMP(1,3),TDUMP(1,IPOSZ),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETU) CALL MAPDMP
     $         (SDUMP(1,4),TDUMP(1,IPOSU),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETU) CALL MAPDMP
     $         (SDUMP(1,5),TDUMP(1,IPOSV),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETW) CALL MAPDMP
     $         (SDUMP(1,6),TDUMP(1,IPOSW),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETP) CALL MAPDMP
     $         (SDUMP(1,7),TDUMP(1,IPOSP),IEG,NXR,NYR,NZR,ifbytsw)
               IF (IFGETT) CALL MAPDMP
     $         (SDMP2(1,1),TDUMP(1,IPOST),IEG,NXR,NYR,NZR,ifbytsw)

C              passive scalars
               DO 100 IPS=1,NPSCAL
                  IF (IFGTPS(IPS)) CALL MAPDMP(SDMP2(1,IPS+1)
     $               ,TDUMP(1,IPSPS(IPS)),IEG,NXR,NYR,NZR,ifbytsw)
  100          CONTINUE
 
  200       CONTINUE               
C
C           Successfully read a complete field, store it.
C
            nerr = 0              ! Count number of elements rec'd by nid
            do ieg=1,nelrr
               mid = gllnid(ieg)
               if (mid.eq.nid) nerr = nerr+1
            enddo

            nxyz2=nx2*ny2*nz2
            nxyz1=nx1*ny1*nz1
            ntott=nerr*nxyz1
            ntotv=nerr*nxyz1   ! Problem for differing Vel. and Temp. counts!
                               ! for now we read nelt dataset
c
            if (ifmhd.and.ifile.eq.2) then
               if (ifgetu) call copy(bx,sdump(1,4),ntott)
               if (ifgetu) call copy(by,sdump(1,5),ntott)
               if (ifgetw) call copy(bz,sdump(1,6),ntott)
               if (ifgetp) then
                 if (nid.eq.0) write(6,*) 'getting restart pressure'
                 if (ifsplit) then
                    call copy( pm,sdump(1,7),ntotv)
                 else
                  do iel=1,nelv
                    iiel = (iel-1)*nxyz1+1
                    call map12 (pm(1,1,1,iel),sdump(iiel,7),iel)
                  enddo
                 endif
               endif
               if (ifaxis.and.ifgett) 
     $            call copy(t(1,1,1,1,2),sdmp2(1,1),ntott)
            elseif (ifpert.and.ifile.ge.2) then
               j=ifile-1  ! pointer to perturbation field
               if (ifgetu) call copy(vxp(1,j),sdump(1,4),ntotv)
               if (ifgetu) call copy(vyp(1,j),sdump(1,5),ntotv)
               if (ifgetw) call copy(vzp(1,j),sdump(1,6),ntotv)
               if (ifgetp) then
                  if (nid.eq.0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                     call copy(prp(1,j),sdump(1,7),ntotv)
                  else
                     do ie=1,nelv
                        ie1 = (ie-1)*nxyz1+1
                        ie2 = (ie-1)*nxyz2+1
                        call map12 (prp(ie2,j),sdump(ie1,7),ie)
                     enddo
                  endif
               endif
               if (ifgett) call copy(tp(1,1,j),sdmp2(1,1),ntott)
C              passive scalars
               do ips=1,NPSCAL
                  if (ifgtps(ips))
     $            call copy(tp(1,ips+1,j),sdmp2(1,ips+1),ntott)
               enddo
c
            else  ! Std. Case
               if (ifgetx) call copy(xm1,sdump(1,1),ntott)
               if (ifgetx) call copy(ym1,sdump(1,2),ntott)
               if (ifgetz) call copy(zm1,sdump(1,3),ntott)
               if (ifgetu) call copy(vx ,sdump(1,4),ntotv)
               if (ifgetu) call copy(vy ,sdump(1,5),ntotv)
               if (ifgetw) call copy(vz ,sdump(1,6),ntotv)
               if (ifgetp) then
                  if (nid.eq.0) write(6,*) 'getting restart pressure'
                  if (ifsplit) then
                     call copy(pr,sdump(1,7),ntotv)
                  else
                     do iel=1,nelv
                        iiel = (iel-1)*nxyz1+1
                        call map12 (pr(1,1,1,iel),sdump(iiel,7),iel)
                     enddo
                  endif
               endif
               if (ifgett) call copy(t,sdmp2(1,1),ntott)
C              passive scalars
               do i=1,NPSCAL
                  if (ifgtps(i))
     $            call copy(t(1,1,1,1,i+1),sdmp2(1,i+1),ntott)
               enddo
c
               IF (IFAXIS) THEN

               DO IEL=1,NELV
               IF(IFRZER(IEL)) THEN
                 IF(IFGETX) THEN
                   CALL MXM   (XM1(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY  (XM1(1,1,1,IEL),AXISM1,NX1*NY1)
                   CALL MXM   (YM1(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY  (YM1(1,1,1,IEL),AXISM1,NX1*NY1)
                 ENDIF
                 IF(IFGETU) THEN
                   CALL MXM    (VX(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY   (VX(1,1,1,IEL),AXISM1,NX1*NY1)
                   CALL MXM    (VY(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY   (VY(1,1,1,IEL),AXISM1,NX1*NY1)
                 ENDIF
                 IF(IFGETW) THEN
                   CALL MXM    (VZ(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY   (VZ(1,1,1,IEL),AXISM1,NX1*NY1)
                 ENDIF
                 IF(IFGETP) THEN
                   CALL MXM    (PR(1,1,1,IEL),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY   (PR(1,1,1,IEL),AXISM1,NX1*NY1)
                 ENDIF
                 IF(IFGETT) THEN
                   CALL MXM  (T (1,1,1,IEL,1),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY (T (1,1,1,IEL,1),AXISM1,NX1*NY1)
                 ENDIF
                 DO IPS=1,NPSCAL
                  IS1 = IPS + 1
                  IF(IFGTPS(IPS)) THEN
                   CALL MXM (T(1,1,1,IEL,IS1),NX1,IATLJ1,NY1,AXISM1,NY1)
                   CALL COPY(T(1,1,1,IEL,IS1),AXISM1,NX1*NY1)
                  ENDIF
                 ENDDO
               ENDIF
               ENDDO
               ENDIF
c
               if (ifgtim) time=rstime
            endif

 1000    CONTINUE
         GOTO 1600
C
C           Else, end of file found - notify other processors.
 1500       CONTINUE
            IFOK=.FALSE.
            call lbcast(ifok)
C
 1600    CONTINUE
C
         IF (IDUMP.EQ.1.AND.NID.EQ.0) THEN
            WRITE(6,1700) FNAME
            WRITE(6,1701) IEG,IXYZ
            WRITE(6,1702) 
     $            ((TDUMP(JXYZ,II),II=1,NOUTS),JXYZ=IXYZ-1,IXYZ)
 1700       FORMAT(5X,'WARNING:  No data read in for file ',A80)
 1701       FORMAT(5X,'Failed on  element',I4,',  point',I5,'.')
 1702       FORMAT(5X,'Last read dump:',/,5G15.7)
            write(6,*) nid,'call exitt 1702a',idump
            call exitt
         ELSEIF (IDUMP.EQ.1) THEN
            write(6,*) nid,'call exitt 1702b',idump
            call exitt
         ELSE
            IDUMP=IDUMP-1
            IF (NID.EQ.0) WRITE(6,1800) IDUMP
 1800       FORMAT(2X,'Successfully read data from dump number',I3,'.')
         ENDIF
         if (iffmat) then
            if (nid.eq.0) close(unit=91)
         else
            if (nid.eq.0) call byte_close()
         endif
         GOTO 6000
C
C        Can't open file...
 5000    CONTINUE
         IF (NID.EQ.0) WRITE(6,5001) FNAME 
 5001    FORMAT(2X,'   *******   ERROR   *******    '
     $       ,/,2X,'   *******   ERROR   *******    '
     $       ,/,2X,'   Could not open restart file:'
     $       ,/,A80
     $      ,//,2X,'Quitting in routine RESTART.')
         CLOSE(UNIT=91)
         call exitt
 5002    CONTINUE
         IF (NID.EQ.0) WRITE(6,5001) HNAME 
         call exitt
C
C
C     End of IFILE loop
 6000 CONTINUE
C
      return
      end
C
c-----------------------------------------------------------------------
      subroutine sioflag(ndumps,fname,rsopts)
C
C     Set IO flags according to Restart Options File, RSOPTS
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'RESTART'
      INCLUDE 'TSTEP'
C
      CHARACTER*80 RSOPTS,FNAME
      CHARACTER*2  S2
      LOGICAL IFGTRL
C
C     Scratch variables..
      LOGICAL IFDEFT,IFANYC
      CHARACTER*80 RSOPT     ,LINE
      CHARACTER*1  RSOPT1(80),LINE1(80)
      EQUIVALENCE (RSOPT1,RSOPT)
      EQUIVALENCE (LINE1,LINE)
C
C     Parse filename
C
C        CSPLIT splits S1 into two parts, delimited by S2.  
C        S1 returns with 2nd part of S1.  CSPLIT returns 1st part.
C
      RSOPT=RSOPTS
      CALL LJUST(RSOPT)
      CALL CSPLIT(FNAME,RSOPT,' ',1)
C     check fname for user supplied extension.
      IF (INDX1(FNAME,'.',1).EQ.0) THEN
         LEN=LTRUNC(FNAME,80)
         LEN1=LEN+1
         LEN4=LEN+4
         FNAME(LEN1:LEN4)='.fld'
      ENDIF
C
C     Parse restart options
C
C        set default flags
C
      IFGETX=.FALSE.
      IFGETZ=.FALSE.
      IFGETU=.FALSE.
      IFGETW=.FALSE.
      IFGETP=.FALSE.
      IFGETT=.FALSE.
      DO 100 I=1,LDIMT1
         IFGTPS(I)=.FALSE.
  100 CONTINUE
      IFGTIM=.TRUE.
      NDUMPS=0
C
C     Check for default case - just a filename given, no i/o options specified
C
      IFDEFT=.TRUE.
C
C     Parse file for i/o options and/or dump number
C
      CALL CAPIT(RSOPT,80)


      IF (LTRUNC(RSOPT,80).NE.0) THEN
C
C        Check for explicit specification of restart TIME.
C
         ITO=INDX2(RSOPT,'TIME',4)
         IFGTIM=.TRUE.
         IF (ITO.NE.0) THEN
C           user has specified the time explicitly.
            IT1=INDX2(RSOPT,'=',1)
            IT8=80-IT1
            CALL BLANK(LINE,80)
            CALL CHCOPY(LINE,RSOPT1(IT1),IT8)
            IF (IFGTRL(TTIME,LINE)) THEN
               IFGTIM=.FALSE.
               TIME=TTIME
            ENDIF
C           remove the user specified time from the RS options line.
            ITA=80-ITO+1
            CALL BLANK(RSOPT1(ITO),ITA)
            CALL LJUST(LINE)
            IT1=INDX1(LINE,' ',1)
            ITB=80-IT1+1
            CALL CHCOPY(RSOPT1(ITO),LINE1(IT1),ITB)
         ENDIF
C
C        Parse field specifications.
C
         IXO=INDX2(RSOPT,'X',1)
         IF (IXO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETX=.TRUE.
            IF (IF3D) IFGETZ=.TRUE.
         ENDIF
C
         IVO=INDX2(RSOPT,'U',1)
         IF (IVO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETU=.TRUE.
            IF (IF3D) IFGETW=.TRUE.
         ENDIF
C
         IPO=INDX2(RSOPT,'P',1)
         IF (IPO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETP=.TRUE.
         ENDIF
C
         ITO=INDX2(RSOPT,'T',1)
         IF (ITO.NE.0) THEN
            IFDEFT=.FALSE.
            IFGETT=.TRUE.
         ENDIF
C
         DO 300 I=2,NPSCAL
            WRITE (S2,301) I
            IPO=INDX2(RSOPT,S2,2)
            IF (IPO.NE.0) THEN
               IFDEFT=.FALSE.
               IFGTPS(I)=.TRUE.
            ENDIF
  300    CONTINUE
  301    FORMAT('P',I1)
C
C        Get number of dumps from remainder of user supplied line.
C
         IF (IFGTRL(TDUMPS,RSOPT)) NDUMPS=INT(TDUMPS)
      ENDIF
C
C     If no fields were explicitly specified, assume getting all fields. 
C
      IF (IFDEFT) THEN
         IF (IFXYO) THEN
            IFGETX=.TRUE.
            IF (IF3D) IFGETZ=.TRUE.
         ENDIF
         IFANYC=.FALSE.
         DO 400 I=1,NFIELD
            IF (IFADVC(I)) IFANYC=.TRUE.
  400    CONTINUE
         IF (IFFLOW.OR.IFANYC) THEN
            IFGETU=.TRUE.
            IF (IF3D) IFGETW=.TRUE.
         ENDIF
         IF (IFFLOW) IFGETP=.TRUE.
         IF (IFHEAT) IFGETT=.TRUE.
         DO 410 I=1,NPSCAL
            IFGTPS(I)=.TRUE.
  410    CONTINUE
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mapdmp(sdump,tdump,ieg,nxr,nyr,nzr,if_byte_sw)
C----------------------------------------------------------------------
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'PARALLEL'
C
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
C
      REAL   SDUMP(LXYZ1,LELT)
      REAL*4 TDUMP(LXYZR)
c
      logical if_byte_sw
c
      NXYZ=NX1*NY1*NZ1
      NXYR=NXR*NYR*NZR
C
C     Serial processor code:
C
      IF (NP.EQ.1) THEN
C
         IF (if_byte_sw)           CALL byte_reverse(TDUMP,NXYR)
         IF (NXR.EQ.NX1.AND.NYR.EQ.NY1.AND.NZR.EQ.NZ1) THEN
            CALL COPY4r(SDUMP(1,IEG),TDUMP,NXYZ)
         ELSE
C           do the map    (assumes that NX=NY=NZ, or NX=NY, NZ=1)
            call mapab4r(sdump(1,ieg),tdump,nxr,1)
         ENDIF
C
       ELSE
C
C     Parallel code - send data to appropriate processor and map.
C
         JNID=GLLNID(IEG)
         MTYPE=3333+IEG
         LEN=4*NXYR
         LE1=4
         IF (NID.EQ.0.AND.JNID.NE.0) THEN
c           hand-shake
            CALL CSEND(MTYPE,TDUMP,LE1,JNID,NULLPID)
            CALL CRECV(MTYPE,dummy,LE1)
            CALL CSEND(MTYPE,TDUMP,LEN,JNID,NULLPID)
         ELSEIF (NID.NE.0.AND.JNID.EQ.NID) THEN
C           Receive data from node 0
            CALL CRECV(MTYPE,dummy,LE1)
            CALL CSEND(MTYPE,TDUMP,LE1,0,NULLPID)
            CALL CRECV(MTYPE,TDUMP,LEN)
         ENDIF
C
C        If the data is targeted for this processor, then map 
C        to appropriate element.
C
         IF (JNID.EQ.NID) THEN
            IE=GLLEL(IEG)
            IF (if_byte_sw)           CALL byte_reverse(TDUMP,NXYR)
            IF (NXR.EQ.NX1.AND.NYR.EQ.NY1.AND.NZR.EQ.NZ1) THEN
               CALL COPY4r(SDUMP(1,IE),TDUMP,NXYZ)
            ELSE
               call mapab4r(sdump(1,ie),tdump,nxr,1)
            ENDIF
         ENDIF
C
C        End of parallel distribution/map routine.
C
      ENDIF
      return
      END
c-----------------------------------------------------------------------
      subroutine mapab(x,y,nxr,nel)
C---------------------------------------------------------------
C
C     Interpolate Y(NXR,NYR,NZR,NEL) to X(NX1,NY1,NZ1,NEL)
C     (assumes that NXR=NYR=NZR, or NXR=NYR, NZR=1)
C---------------------------------------------------------------
C
      INCLUDE 'SIZE'
      INCLUDE 'IXYZ'
      INCLUDE 'WZ'
C
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      DIMENSION X(NX1,NY1,NZ1,NEL)
      DIMENSION Y(NXR,NXR,NXR,NEL)

      common /ctmpab/ xa(lxyzr)      ,xb(lx1,ly1,lzr) ,xc(lxyzr)
     $              , ires(lxr,lxr)  ,itres(lxr,lxr)
     $              , zgmr(lxr)      ,wgtr(lxr)
      real ires,itres

      INTEGER NOLD
      SAVE    NOLD
      DATA    NOLD /0/
C
      NZR = NXR
      IF(NZ1.EQ.1) NZR=1
      NYZR = NXR*NZR
      NXY1 = NX1*NY1
C
      IF (NXR.NE.NOLD) THEN
         NOLD=NXR
         CALL ZWGLL   (ZGMR,WGTR,NXR)
         CALL IGLLM   (IRES,ITRES,ZGMR,ZGM1,NXR,NX1,NXR,NX1)      
         IF (NID.EQ.0) WRITE(6,10) NXR,NX1
   10       FORMAT(2X,'Mapping restart data from Nold=',I2
     $               ,' to Nnew=',I2,'.')
      ENDIF
C
      DO 1000 IE=1,NEL
         CALL MXM (IRES,NX1,Y(1,1,1,IE),NXR,XA,NYZR)
         DO 100 IZ=1,NZR
            IZOFF = 1 + (IZ-1)*NX1*NXR
            CALL MXM (XA(IZOFF),NX1,ITRES,NXR,XB(1,1,IZ),NY1)
  100    CONTINUE
         IF (NDIM.EQ.3) THEN
            CALL MXM (XB,NXY1,ITRES,NZR,X(1,1,1,IE),NZ1)
         ELSE
            CALL COPY(X(1,1,1,IE),XB,NXY1)
         ENDIF
 1000 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine mapab4R(x,y,nxr,nel)
C---------------------------------------------------------------
C
C     Interpolate Y(NXR,NYR,NZR,NEL) to X(NX1,NY1,NZ1,NEL)
C     (assumes that NXR=NYR=NZR, or NXR=NYR, NZR=1)
c
c     Input:  real*4,  Output:  default precision
c
C---------------------------------------------------------------
C
      INCLUDE 'SIZE'
      INCLUDE 'IXYZ'
      INCLUDE 'WZ'
C
      PARAMETER (LXR=LX1+6)
      PARAMETER (LYR=LY1+6)
      PARAMETER (LZR=LZ1+6)
      PARAMETER (LXYZR=LXR*LYR*LZR)
      PARAMETER (LXYZ1=LX1*LY1*LZ1)
      REAL*4 X(NX1,NY1,NZ1,NEL)
      REAL   Y(NXR,NXR,NXR,NEL)

      common /ctmpab/ xa(lxyzr)      ,xb(lx1,ly1,lzr) ,xc(lxyzr)
     $              , ires(lxr,lxr)  ,itres(lxr,lxr)
     $              , zgmr(lxr)      ,wgtr(lxr)
      real ires,itres

      INTEGER NOLD
      SAVE    NOLD
      DATA    NOLD /0/

      NZR = NXR
      IF(NZ1.EQ.1) NZR=1
      NYZR = NXR*NZR
      NXY1 = NX1*NY1
      nxyzr = nxr*nxr*nzr
C
      IF (NXR.NE.NOLD) THEN
         NOLD=NXR
         CALL ZWGLL   (ZGMR,WGTR,NXR)
         CALL IGLLM   (IRES,ITRES,ZGMR,ZGM1,NXR,NX1,NXR,NX1)      
         IF (NID.EQ.0) WRITE(6,10) NXR,NX1
   10       FORMAT(2X,'Mapping restart data from Nold=',I2
     $               ,' to Nnew=',I2,'.')
      ENDIF
C
      DO 1000 IE=1,NEL
         call copy4r(xc,y(1,1,1,ie),nxyzr)
         CALL MXM (IRES,NX1,xc,NXR,XA,NYZR)
         DO 100 IZ=1,NZR
            IZOFF = 1 + (IZ-1)*NX1*NXR
            CALL MXM (XA(IZOFF),NX1,ITRES,NXR,XB(1,1,IZ),NY1)
  100    CONTINUE
         IF (NDIM.EQ.3) THEN
            CALL MXM (XB,NXY1,ITRES,NZR,X(1,1,1,IE),NZ1)
         ELSE
            CALL COPY(X(1,1,1,IE),XB,NXY1)
         ENDIF
 1000 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX132(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=80-L2+1
      INDX132=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX132=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      INDX1=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            return
         ENDIF
  100 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX2(S1,S2,L2)
C
C     INDX2 is returned with the location of S2 in S1 (0 if not found)
C     S1     is returned with 1st occurance of S2 removed.
C
      CHARACTER*1 S1(80),S2(80)
C
      I1=INDX1(S1,S2,L2)
C
      IF (I1.NE.0) THEN
C
         N1=80-L2
         DO 100 I=I1,N1
            I2=I+L2
C           remove the 1st occurance of S2 from S1.
            S1(I)=S1(I2)
  100    CONTINUE
         N2=N1+1
         DO 200 I=N2,80
            S1(I)=' '
  200    CONTINUE
      ENDIF
C
      INDX2=I1
      return
      END
c-----------------------------------------------------------------------
      subroutine csplit(s0,s1,s2,l0)
      CHARACTER*80 S0,S1,S2
C     split string S1 into two parts, delimited by S2.
C
      I2=INDX2(S1,S2,L0)
      IF (I2.EQ.0) return
C
      I1=I2-1
      CALL BLANK(S0,80)
      S0(1:I1)=S1(1:I1)
      CALL LSHFT(S1,I2)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine lshft(string,ipt)
C     shift string from IPT to the left
C     INPUT : "abcde......    test    "
C     OUTPUT: "e......    test        "     if ipt.eq.5
      CHARACTER*1 STRING(80)
C
      DO 20 J=1,81-IPT
         IJ=IPT+J-1
         STRING(J)=STRING(IJ)
   20 CONTINUE
      DO 30 J=82-IPT,80
         STRING(J)=' '
   30 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ljust(string)
C     left justify string
      CHARACTER*1 STRING(80)
C
      IF (STRING(1).NE.' ') return
C
      DO 100 I=2,80
C
         IF (STRING(I).NE.' ') THEN
            DO 20 J=1,81-I
               IJ=I+J-1
               STRING(J)=STRING(IJ)
   20       CONTINUE
            DO 30 J=82-I,80
               STRING(J)=' '
   30       CONTINUE
            return
         ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine chknorm (ifzero)
C----------------------------------------------------------------------
C
C     Check if trivial user specified initial conditions
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      LOGICAL IFZERO
C
      IFZERO = .TRUE.
C
      IF (IFFLOW) THEN
         IFIELD = 1
         IMESH  = 1
         CALL UNORM
         IF (VNRML8.GT.0.) IFZERO = .FALSE.
      ENDIF
      IF (IFHEAT) THEN
         DO 100 IFIELD=2,NFIELD
            IMESH = 1
            IF (IFTMSH(IFIELD)) IMESH = 2
            CALL UNORM
            IF (TNRML8(IFIELD).GT.0.) IFZERO = .FALSE.
 100     CONTINUE
      ENDIF
c
      return
      END
C
c-----------------------------------------------------------------------
      subroutine prsolvt
C----------------------------------------------------------------------
C
C     Use steady state solution as initial condition 
C     for temperatur/passive scalar
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TSTEP'
      LOGICAL  IFSAV1,IFSAV2(LDIMT1)
C
      IF (NID.EQ.0) WRITE(6,*) ' '
      IF (NID.EQ.0) WRITE(6,*) 'Conduction pre-solver activated'
C
C     Set logical IFTRAN to false (steady state)
C     Save logicals for convection
C     Turn convection off
C
      IFSAV1 = IFTRAN
      IFTRAN = .FALSE.
      DO 100 IFIELD=2,NFIELD
         IFSAV2(IFIELD) = IFADVC(IFIELD)
         IFADVC(IFIELD) = .FALSE.
 100  CONTINUE
C
      CALL SETPROP
      CALL SETSOLV
C
      IF(NID.EQ.0)WRITE(6,*)'Steady conduction/passive scalar problem'
C
      DO 200 IGEOM=1,2
         CALL HEAT (IGEOM)
 200  CONTINUE
C
C     Set IFTRAN to true again
C     Turn convection on again
C
      IFTRAN = IFSAV1
      DO 300 IFIELD=2,NFIELD
         IFADVC(IFIELD) = IFSAV2(IFIELD)
 300  CONTINUE
C
      return
      END
C
c-----------------------------------------------------------------------
      subroutine prsolvv
C----------------------------------------------------------------------
C
C     Use steady Stokes solution as initial condition 
C     for flow problem
C
C----------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      LOGICAL  IFSAV1,IFSAV2
C
      IF (NID.EQ.0) WRITE(6,*) ' '
      IF (NID.EQ.0) WRITE(6,*) 'Velocity pre-solver activated'
C
C     Initialize velocity to some non-trivial RHS to avoid FP trap in i860.
C
      IF (PARAM(60).NE.0.0) THEN
         SMALL=10.0E-10
         CALL CFILL(VX,SMALL,NTOTV)
         CALL CFILL(VY,SMALL,NTOTV)
         CALL CFILL(VZ,SMALL,NTOTV)
      ENDIF
C
C     Set logical IFTRAN to false (steady state)
C     Save logicals for convection
C     Turn convection off
C
      IF (IFSPLIT) THEN
        WRITE(6,10)
   10   FORMAT(
     $ /,2X,'ERROR: Steady Stokes Flow initial condition cannot'
     $,/,2X,'       be computed using the splitting formulation.'
     $,/,2X,'       Either compute using UZAWA, or remove PRESOLVE'
     $,/,2X,'       request for velocity.'
     $,/,2X
     $,/,2X,'       ABORTING IN PRSOLVV.')
        CALL EXITT
      ENDIF
C
      IFSAV1 = IFTRAN
      IFSAV2 = IFADVC(IFIELD)
      IFTRAN = .FALSE.
      IFADVC(IFIELD) = .FALSE.
C
      CALL SETPROP
      CALL SETSOLV
      IF (IFNATC) GTHETA = GTHETA+10.
C
      IF (NID.EQ.0) WRITE (6,*) 'Steady Stokes problem'
      DO 100 IGEOM=1,2
         IF (.NOT.IFSPLIT) CALL FLUID (IGEOM)
 100  CONTINUE
      IF (IFNATC) GTHETA = GTHETA-10.
C
C     Set IFTRAN to true again
C     Turn convection on again
C
      IFTRAN = IFSAV1
      IFADVC(IFIELD) = IFSAV2
C
      return
      END
C
c-----------------------------------------------------------------------
      subroutine nekuic
C------------------------------------------------------------------
C
C     User specified fortran function (=0 if not specified)
C
C------------------------------------------------------------------
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      INCLUDE 'TURBO'
      INCLUDE 'PARALLEL'
      INCLUDE 'NEKUSE'
C
      NEL   = NELFLD(IFIELD)
C
      IF (IFMODEL .AND. IFKEPS .AND. IFIELD.EQ.IFLDK) THEN
C
      DO 100 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 100 K=1,NZ1
         DO 100 J=1,NY1
         DO 100 I=1,NX1
            CALL NEKASGN (I,J,K,IEL)
            CALL USERIC  (I,J,K,IEG)
            T(I,J,K,IEL,IFIELD-1) = TURBK
 100  CONTINUE
C
      ELSEIF (IFMODEL .AND. IFKEPS .AND. IFIELD.EQ.IFLDE) THEN
C
      DO 200 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 200 K=1,NZ1
         DO 200 J=1,NY1
         DO 200 I=1,NX1
            CALL NEKASGN (I,J,K,IEL)
            CALL USERIC  (I,J,K,IEG)
            T(I,J,K,IEL,IFIELD-1) = TURBE
 200  CONTINUE
C
      ELSE
C
      DO 300 IEL=1,NEL
         IEG = LGLEL(IEL,NODE)
         DO 300 K=1,NZ1
         DO 300 J=1,NY1
         DO 300 I=1,NX1
           CALL NEKASGN (I,J,K,IEL)
           CALL USERIC  (I,J,K,IEG)
           if (jp.eq.0) then
             IF (IFIELD.EQ.1) THEN
               VX(I,J,K,IEL) = UX
               VY(I,J,K,IEL) = UY
               VZ(I,J,K,IEL) = UZ
             ELSEIF (IFIELD.EQ.ifldmhd) THEN
               BX(I,J,K,IEL) = UX
               BY(I,J,K,IEL) = UY
               BZ(I,J,K,IEL) = UZ
             ELSE
               T(I,J,K,IEL,IFIELD-1) = TEMP
             ENDIF
           else
             ijke = i+nx1*((j-1)+ny1*((k-1) + nz1*(iel-1)))
             IF (IFIELD.EQ.1) THEN
               VXP(IJKE,jp) = UX
               VYP(IJKE,jp) = UY
               VZP(IJKE,jp) = UZ
             ELSE
               TP(IJKE,IFIELD-1,jp) = TEMP
             ENDIF
           endif

 300  CONTINUE
C
      ENDIF
c
      return
      END
c-----------------------------------------------------------------------
      subroutine capit(lettrs,n)
C     Capitalizes string of length n
      CHARACTER LETTRS(N)
C
      DO 5 I=1,N
         INT=ICHAR(LETTRS(I))
         IF(INT.GE.97 .AND. INT.LE.122) THEN
            INT=INT-32
            LETTRS(I)=CHAR(INT)
         ENDIF
5     CONTINUE
      return
      END
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFGTRL(VALUE,LINE)
C
C     Read VALUE from LINE and set IFGTRL to .TRUE. if successful,
C                                  IFGTRL to .FALSE. otherwise.
C
C     This complicated function is necessary thanks to the Ardent,
C     which won't allow free formatted reads (*) from internal strings!
C
      CHARACTER*80 LINE
      CHARACTER*80 WORK
      CHARACTER*8  FMAT
C
C     Note that the format Fn.0 is appropriate for fields of type:
C          34   34.0  34.0e+00
C     The only difficulty would be with '34' but since we identify
C     the field width exactly, there is no problem.
C
      IFGTRL=.FALSE.
      VALUE=0.0
C
      WORK=LINE
      CALL LJUST(WORK)
      IFLDW=INDX1(WORK,' ',1)-1
C
      IF (IFLDW.GT.0) THEN
         WRITE(FMAT,10) IFLDW
   10    FORMAT('(F',I3.3,'.0)')
         READ(WORK,FMAT,ERR=100,END=100) TVAL
         VALUE=TVAL
         IFGTRL=.TRUE.
         return
      ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      LOGICAL FUNCTION IFGTIL(IVALUE,LINE)
C
C     Read IVALUE from LINE and set IFGTIL to .TRUE. if successful,
C                                   IFGTIL to .FALSE. otherwise.
C
C     This complicated function is necessary thanks to the Ardent,
C     which won't allow free formatted reads (*) from internal strings!
C
      CHARACTER*80 LINE
      CHARACTER*80 WORK
      CHARACTER*8  FMAT
C
      IFGTIL=.FALSE.
      IVALUE=0
C
      WORK=LINE
      CALL LJUST(WORK)
      IFLDW=INDX1(WORK,' ',1)-1
C
      IF (IFLDW.GT.0) THEN
         WRITE(FMAT,10) IFLDW
   10    FORMAT('(F',I3.3,'.0)')
         READ(WORK,FMAT,ERR=100,END=100) TVAL
         IVALUE=INT(TVAL)
         IFGTIL=.TRUE.
         return
      ENDIF
C
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine perturb(tt,ifld,eps)
      include 'SIZE'
      include 'TOTAL'
c
      real tt(1)
      integer ifld

      ifield = ifld

      n = nx1*ny1*nz1*nelfld(ifield)
      call vcospf(tt,bm1,n)
      call cmult(tt,eps,n)
      call dssum(tt,nx1,ny1,nz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine vcospf(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = cos(1000.*y(i))
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine vbyte_swap(x,n)
      character*1 x(0:3,1),tmp0,tmp1
      character*1 in (0:3), out(0:3)
      real*4      in4     , out4
      equivalence (in ,in4 )
      equivalence (out,out4)
c
      do i=1,n
         do j=0,3
            in (j) = x(j,i)
         enddo
         tmp0   = x(0,i)
         tmp1   = x(1,i)
         x(0,i) = x(3,i)
         x(1,i) = x(2,i)
         x(2,i) = tmp1
         x(3,i) = tmp0
         do j=0,3
            out(j) = x(j,i)
         enddo
         write(6,*) 'swap:',i,in4,out4
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)
      include 'SIZE'
c
      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.
c
      test2 = bytetest
      call byte_reverse(test2,1)
      if (nid.eq.0) 
     $   write(6,*) 'byte swap:',if_byte_swap_test,bytetest,test2
      return
      end
c-----------------------------------------------------------------------
      subroutine geom_reset(icall)
C
C     Generate geometry data
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'GEOM'
      INCLUDE 'SOLN'
      INCLUDE 'TSTEP'
      include 'WZ'
c
      COMMON /scruz/ XM3 (LX1,LY1,LZ1,LELT)
     $ ,             YM3 (LX1,LY1,LZ1,LELT)
     $ ,             ZM3 (LX1,LY1,LZ1,LELT)
C
c
      if(nid.eq.0) write(6,*) 'regenerate geomerty data',icall

      ntot = nx1*ny1*nz1*nelt
c
      if (lx3.eq.lx1) then
         call copy(xm3,xm1,ntot)
         call copy(ym3,ym1,ntot)
         call copy(zm3,zm1,ntot)
      else
         call map13_all(xm3,xm1)
         call map13_all(ym3,ym1)
         if (if3d) call map13_all(zm3,zm1)
      endif
c
      CALL GEOM1 (XM3,YM3,ZM3)
      CALL GEOM2
      CALL UPDMSYS (1)
      CALL VOLUME
      CALL SETINVM
      CALL SETDEF
      CALL SFASTAX
c
      if(nid.eq.0) then
        write(6,*) 'done :: regenerate geomerty data',icall
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dsavg(u)
c
c
      include 'SIZE'
      include 'TOTAL'
      real u(lx1,ly1,lz1,lelt)
c
c
c     Take direct stiffness avg of u
c
c
      ifieldo = ifield
      if (ifflow) then
         ifield = 1
         ntot = nx1*ny1*nz1*nelv
         call dssum(u,nx1,ny1,nz1)
         call col2 (u,vmult,ntot)
      else
         ifield = 2
         ntot = nx1*ny1*nz1*nelt
         call dssum(u,nx1,ny1,nz1)
         call col2 (u,tmult,ntot)
      endif
      ifield = ifieldo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine map13_all(x3,x1)
c
      include 'SIZE'
      include 'TOTAL'
c
      real x3(lx3,ly3,lz3,lelt)
      real x1(lx1,ly1,lz1,lelt)
c
      integer e
c
      do e=1,nelt
         call map13 (x3(1,1,1,e),x1(1,1,1,e),e)
      enddo
c
      return
      end
c-----------------------------------------------------------------------

c     Parallel file i/o reader  (param(66) = 6)
c
c     6/14/2006  pff
c
c     File format:

c     132 byte header  (containing nel_block, etc.)
c     4-byte swap test
c     nel_block 4-byte integers, indicating elements contained in file
c     nxyz . nel_block vec data
c     nxyz . nel_block vec data
c     nxyz . nel_block scalar data


c-----------------------------------------------------------------------
      subroutine mfi_gets(u,wk,lwk,idummy)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1)
      real*4 wk(lwk)
      integer e,eg,msg_id(lelt)
      logical idummy

      common /scruz/ w2(4*lx1*ly1*lz1*lelt)

      call gsync() ! clear outstanding message queues.

      num_recv  = nxr*nyr*nzr*nelv
      num_avail = lwk*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_gets a')

      num_read  = nxr*nyr*nzr
      num_avail = 4*lx1*ly1*lz1*lelv*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_gets b')

      len   = nxr*nyr*nzr*wdsizr
      nxyzr = nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      l = 1
      if (np.gt.1) then
         do e=1,nelt
            eg = lglel(e,node)
            msg_id(e) = irecv(eg,wk(l),len)
            l = l+nxyzw
         enddo
      endif

      if (nid.eq.pid0.and.np.gt.1) then ! open file fid0
         do e=1,nelr
            jnid = gllnid(er(e))
            call byte_read(w2,nxyzw)
            call csend(er(e),w2,len,jnid,0)  ! blocking send
         enddo
      elseif (np.eq.1) then
         l = 1
         do e=1,nelr
            call byte_read(wk(l),nxyzw)
            l = l+nxyzw
         enddo
      endif

      if (if_byte_sw.and.wdsizr.eq.8) 
     &   write(6,*) 'ABORT: no support for 64-bit and byte_swap'
      if (if_byte_sw.and.wdsizr.eq.8) call exitt

      if (idummy) goto 100

      l = 1
      do e=1,nelt
         if (np.gt.1) call msgwait(msg_id(e))
         if (if_byte_sw) call byte_reverse(wk(l),nxyzr)
         if (nxr.eq.nx1.and.nyr.eq.ny1.and.nzr.eq.nz1) then
            if (wdsizr.eq.4) then         ! COPY
               call copy4r(u(1,e),wk(l),nxyzr)
            else
               call copy  (u(1,e),wk(l),nxyzr)
            endif
         else                             ! INTERPOLATE
            if (wdsizr.eq.4) then
               call mapab4r(u(1,e),wk(l),nxr,1)
            else
               call mapab  (u(1,e),wk(l),nxr,1)
            endif
         endif
         l = l+nxyzw
      enddo

 100  call gsync() ! clear outstanding message queues.

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_getv(u,v,w,wk,lwk,idummy)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'RESTART'

      real u(lx1*ly1*lz1,1),v(lx1*ly1*lz1,1),w(lx1*ly1*lz1,1)
      real*4 wk(lwk)
      integer e,eg,msg_id(lelt)
      logical idummy

      common /scruz/ w2(4*lx1*ly1*lz1*lelt)

      call gsync() ! clear outstanding message queues.

      num_recv  = ndim*nxr*nyr*nzr*nelv
      num_avail = lwk*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv a')

      num_read  = ndim*nxr*nyr*nzr
      num_avail = 4*lx1*ly1*lz1*lelv*wdsize/wdsizr
      call lim_chk(num_recv,num_avail,'     ','     ','mfi_getv b')

      len   = ndim*nxr*nyr*nzr*wdsizr
      nxyzr = ndim*nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzr = 2*nxyzr

      if (np.gt.1) then
         l = 1
         do e=1,nelt
            eg = lglel(e,node)
            msg_id(e) = irecv(eg,wk(l),len)
            l = l+nxyzr
         enddo
      endif

      if (nid.eq.pid0.and.np.gt.1) then ! open file fid0
         do e=1,nelr
            jnid = gllnid(er(e))
            call byte_read(w2,nxyzr)
            call csend(er(e),w2,len,jnid,0)  ! blocking send
         enddo
      elseif (np.eq.1) then
         l = 1
         do e=1,nelr
            call byte_read(wk(l),nxyzr)
            l = l+nxyzr
         enddo
      endif

      if (if_byte_sw.and.wdsizr.eq.8) write(6,*) 'ERROR mfi_getv a'
      if (if_byte_sw.and.wdsizr.eq.8) call exitt

      nxyzr = nxr*nyr*nzr
      nxyzv = ndim*nxr*nyr*nzr
      nxyzw = nxr*nyr*nzr
      if (wdsizr.eq.8) nxyzw = 2*nxyzw

      if (idummy) goto 100

      l = 1
      do e=1,nelt

         if (np.gt.1) call msgwait(msg_id(e))

         if (if_byte_sw) call byte_reverse(wk(l),nxyzv)

         if (nxr.eq.nx1.and.nyr.eq.ny1.and.nzr.eq.nz1) then
            if (wdsizr.eq.4) then         ! COPY
               call copy4r(u(1,e),wk(l        ),nxyzr)
               call copy4r(v(1,e),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy4r(w(1,e),wk(l+2*nxyzw),nxyzr)
            else
               call copy  (u(1,e),wk(l        ),nxyzr)
               call copy  (v(1,e),wk(l+  nxyzw),nxyzr)
               if (if3d) 
     $         call copy  (w(1,e),wk(l+2*nxyzw),nxyzr)
            endif
         else                             ! INTERPOLATE
            if (wdsizr.eq.4) then
               call mapab4r(u(1,e),wk(l        ),nxr,1)
               call mapab4r(v(1,e),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab4r(w(1,e),wk(l+2*nxyzw),nxr,1)
            else
               call mapab  (u(1,e),wk(l        ),nxr,1)
               call mapab  (v(1,e),wk(l+  nxyzw),nxr,1)
               if (if3d) 
     $         call mapab  (w(1,e),wk(l+2*nxyzw),nxr,1)
            endif
         endif
         l = l+ndim*nxyzw
      enddo

 100  call gsync() ! clear outstanding message queues.

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_hdr(hdr)
      include 'SIZE'

      character*132 hdr

      if (indx1(hdr,'#std',4).eq.1) then
          call parse_std_hdr(hdr)
      else
         if (nid.eq.0) write(6,80) hdr
         if (nid.eq.0) write(6,80) 'NONSTD HDR, parse_hdr, abort.'
  80     format(a132)
         call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr(hdr)
      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      character*132 hdr

      character*4 dummy

      read(hdr,*,err=99) dummy
     $         ,  wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         ,  ifiler,nfiler
     $         ,  rdcode      ! 74+20=94

c          NOTE:  This will be extended to general case in future.
c                 For now, what you see in file is what you get.
      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,npscal
         ifgtpsr(k) = .false.
      enddo

      ifgtim = .true.  ! this is the default

      NPS=0
      do i=1,10 
         if (rdcode1(i).eq.'X') ifgetxr = .true.
         if (rdcode1(i).eq.'U') ifgetur = .true.
         if (rdcode1(i).eq.'P') ifgetpr = .true.
         if (rdcode1(i).eq.'T') ifgettr = .true.
         if (npscal.gt.0 .and. rdcode1(i).eq.'S') then
            read(rdcode1(i+1),'(I1)') NPS1
            read(rdcode1(i+2),'(I1)') NPS0
            NPS = 10*NPS1+NPS0
            do k=1,NPS
               ifgtpsr(k) = .true.
            enddo
            ! nothing will follow
            GOTO 50
         endif
      enddo
  
 50   if (nps.gt.(ldimt-1)) then
         if (nid.eq.0) then 
           write(*,'(A)') 
     &      'ERROR: restart file has a NSPCAL > LDIMT'
           write(*,'(A,I2)') 
     &      'Change LDIMT in SIZE'
         endif
         call exitt
      endif

      return

   99 continue   !  If we got here, then the May 2008 variant of std hdr
                 !  failed and we may have an older input file.

      call parse_std_hdr_2006(hdr,rdcode)  ! try the original header format

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_std_hdr_2006(hdr,rlcode)
      include 'SIZE'
      include 'INPUT'
      include 'RESTART'

      character*132 hdr
      character*1 rlcode(20)

c                4  7  10  13   23    33    53    62     68     74
      read(hdr,1) wdsizr,nxr,nyr,nzr,nelr,nelgr,timer,istpr
     $         , ifiler,nfiler
     $         , (rlcode(k),k=1,20)                   ! 74+20=94
    1 format(4x,i2,3i3,2i10,e20.13,i9,2i6,20a1)

c     Assign read conditions, according to rdcode
c          NOTE:  This will be extended to general case in future.
c                 For now, what you see in file is what you get.

      ifgetxr = .false.
      ifgetur = .false.
      ifgetpr = .false.
      ifgettr = .false.
      do k=1,npscal
         ifgtpsr(k) = .false.
      enddo

      ifgtim = .true.  ! this is the default

      if (rlcode(1).eq.'X') ifgetxr = .true.
      if (rlcode(2).eq.'U') ifgetur = .true.
      if (rlcode(3).eq.'P') ifgetpr = .true.
      if (rlcode(4).eq.'T') ifgettr = .true.
      do k=1,npscal
         if (rlcode(4+k).ne.' ') ifgtpsr(k) = .true.
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mfi_get_hdr0(hdr,fname)
c
      include 'SIZE'
      include 'PARALLEL'
      character*132 hdr
      character*80  fname

      if (nid.eq.0) then
         call mbyte_open(fname,0) ! open  blah000.fldnn
         call blank     (hdr,132)
         call byte_read (hdr, 33) ! 4 x 33 = 132
         call byte_close()
      endif
      call bcast(hdr,132)

      return
      end
c-----------------------------------------------------------------------
      subroutine mlti_file_input(fname)
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
c          subsequent files are for B-field or perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      character*132 hdr
      character*80  fname


      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)
      integer e

      call mfi_get_hdr0(hdr,fname)
      call parse_hdr   (hdr)

C     Bounds checking on mapped data.
      IF (NXR.GT.LX1+6) THEN
         if(nid.eq.0) WRITE(6,20) NXR,NX1
   20    FORMAT(//,2X,
     $   'ABORT:  Attempt to map from',I3,
     $   ' to N=',I3,'.',/,2X,
     $   'NEK5000 currently supports mapping from N+6 or less.')
         CALL EXITT
      ENDIF

      call set_pid(fname)! determine reader nodes; open files

      if (ifgetxr) then      ! if available
         if (ifgetx) then
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.false.)
         else                ! skip
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.true.)
         endif
      endif

      if (ifgetur) then
         if (ifgetu) then
            if (ifmhd.and.ifile.eq.2) then
               call mfi_getv(bx,by,bz,wk,lwk,.false.)
            else
               call mfi_getv(vx,vy,vz,wk,lwk,.false.)
            endif
         else
            call mfi_getv(vx,vy,vz,wk,lwk,.true.)
         endif
      endif

      if (ifgetpr) then
         if (ifgetp) then
            call mfi_gets(pm1,wk,lwk,.false.)
            if (ifmhd.and.ifile.eq.2) then
               do e=1,nelv
                  call map12 (pm(1,1,1,e),pm1(1,e),e)
               enddo
            elseif (ifsplit) then
               call copy (pr,pm1,nx2*ny2*nz2*nelv)
            else
               do e=1,nelv
                  call map12 (pr(1,1,1,e),pm1(1,e),e)
               enddo
            endif
         else
            call mfi_gets(pm1,wk,lwk,.true.)
         endif
      endif

      if (ifgettr) then
         if (ifgett) then
            call mfi_gets(t,wk,lwk,.false.)
         else
            call mfi_gets(t,wk,lwk,.true.)
         endif
      endif

      do k=1,npscal
         if (ifgtpsr(k)) then
            if (ifgtps(k)) then
               call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.false.)
            else
               call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.true.)
            endif
         endif
      enddo

      if (ifgtim) time = timer

      if (nid.eq.pid0) call byte_close()

      return
      end
c-----------------------------------------------------------------------
      subroutine mbyte_open(hname,fid) ! open  blah000.fldnn
      include 'SIZE'
      include 'TSTEP'
c
      integer fid
      character*80 hname

      character*6  six,fmt,s6
      save         six
      data         six / "??????" /

      character*80 fname
      character*1  fname1(80)
      equivalence (fname1,fname)

      integer      iname(20)
      equivalence (iname,fname)

      call  izero (iname,20)
      len = ltrunc(hname,80)
      call chcopy (fname,hname,len)

      do ipass=1,2      ! 2nd pass, in case 1 file/directory
         do k=6,1,-1
            i1 = indx1(fname,six,k)
            if (i1.ne.0) then
               write(fmt,1) k,k
    1          format('(i',i1,'.',i1,')')
               write(s6,fmt) fid
               call chcopy(fname1(i1),s6,k)
               goto 10
            endif
         enddo
   10    continue
      enddo
      
      write(6,6) nid,istep,(fname1(k),k=1,len)
    6 format(2i8,' OPEN: ',80a1)

      call byte_open(fname)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_pid(hname)  ! determine which nodes are readers
      character*80 hname

      include 'SIZE'
      include 'PARALLEL'
      include 'RESTART'

      integer stride
      character*132 hdr
      logical if_byte_swap_test
      real*4 bytetest

      stride = np / nfiler
c      write(6,*) nid,' STRIDE:',stride,np,nfiler

      if (stride.lt.1) then
         write(6,*) nfiler,np,'  TOO MANY FILES, set_pid abort'
         call exitt
      endif

      if (mod(nid,stride).eq.0) then
         pid0 = nid
         pid1 = nid + stride
         fid0 = nid / stride

         call mbyte_open(hname,fid0) ! open  blah000.fldnn

         call blank     (hdr,132)
         call byte_read (hdr, 33) ! 4 x 33 = 132
         call parse_hdr (hdr)     ! Re-parse hdr, for nelr info

         call byte_read (bytetest,1)
         if_byte_sw = if_byte_swap_test(bytetest)

         call byte_read (er,nelr) ! Get list of elements to be read
         if (if_byte_sw) call byte_reverse(er,nelr)

      endif

      call lbcast(if_byte_sw)     ! broadcast byte swap from node 0

      return
      end
c-----------------------------------------------------------------------
