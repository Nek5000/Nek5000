c-----------------------------------------------------------------------
      subroutine readat
C
C     Read in data supplied by preprocessor and
C     (eventually) echo check.
C
C     NOTE:  This routine has been broken up into several submodules
C            for ease of comprehension and in anticipation of future
C            restructuring of the NEKTON input.
C
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'ZPER'
C
      if(nid.eq.0) write(6,*) 'read .rea file'

      OPEN (UNIT=9,FILE=REAFLE,STATUS='OLD')

C     Read Parameters
      CALL RDPARAM

C     Read Mesh Data and Group ID

      read(9,*)  ! xfac,yfac,xzero,yzero
      read(9,*)
      read(9,*) nelgs,ndim,nelgv
      nelgt = abs(nelgs)

      if (nid.eq.0) then
         write(6,22) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,22) 'lx1  /lx2  /lx3 :',lx1,lx2,lx3
 22      format(1X,A,4I9)
         write(6,*) ''
      endif

      ifgtp = .false.
      if (ndim.lt.0) ifgtp = .true.     ! domain is a global tensor product

      if ((ifgfdm .or. ifgtp) .and. iand(np,np-1).ne.0) then
         write(6,*)
     $   'For GFDM or GTP, number of processors need to be 2^k'
         call exitt
      endif
c
      call chk_nel  ! make certain sufficient array sizes

      if(nid.eq.0) write(6,*) 'read .map file'
      call mapelpr  ! read .map file, est. gllnid, etc.
      if(nid.eq.0) then
        write(6,*) 'done :: read .map file'
        write(6,*) ''
      endif

      if (nelgs.lt.0) then

         ! new binary reader (nid.eq.0 reads the re2 file)
         ! and sends the data to the other processors
         call bin_rd1

      else

         if(nid.eq.0) write(6,*) 
     &      'ABORT: ASCII no longer supported, use .re2 file!'
         call exitt

         maxrd = 32               ! max # procs to read at once
         mread = (np-1)/maxrd+1   ! mod param
         iread = 0                ! mod param
         x     = 0
         do i=0,np-1,maxrd
            call gsync()
            if (mod(nid,mread).eq.iread) then
               if (ifgtp) then
                  call genbox
               else
                  call rdmesh
                  call rdcurve !  Curved side data
                  call rdbdry  !  Boundary Conditions
               endif
            endif
            iread = iread + 1
         enddo

      endif

C     Read Initial Conditions / Drive Force
      CALL RDICDF

C     Read materials property data
      CALL RDMATP

C     Read history data
      CALL RDHIST

C     Read output specs
      CALL RDOUT

C     Read objects
      CALL RDOBJ
C
C     End of input data, close read file.
C
      CLOSE(UNIT=9)
      if (.not.IFFMTIN) CLOSE(UNIT=8)

      if(nid.eq.0) then
        write(6,*) 'done :: read .rea file'
        write(6,*) ''
      endif
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdparam
C
C     .Read in parameters supplied by preprocessor and
C      (eventually) echo check.
C
C     .Broadcast run parameters to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'ZPER'
C
      READ(9,*,ERR=400)
      READ(9,*,ERR=400) VNEKTON
      NKTONV=(VNEKTON+0.5)
      VNEKMIN=2.5
      IF(VNEKTON.LT.VNEKMIN)THEN
         PRINT*,' Error: This NEKTON Solver Requires a .rea file'
         PRINT*,' from prenek version ',VNEKMIN,' or higher'
         PRINT*,' Please run the session through the preprocessor'
         PRINT*,' to bring the .rea file up to date.'
         call exitt
      ENDIF
      READ(9,*,ERR=400) NDIM
c     error check
      IF(NDIM.NE.LDIM)THEN
         WRITE(6,10) LDIM,NDIM
   10       FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $              /,2X,'       for spatial dimension equal to',I2,'.'
     $              /,2X,'       The data file has dimension',I2,'.')
         call exitt
      ENDIF
      IF (NDIM.EQ.3) IF3D=.TRUE.
      IF (NDIM.NE.3) IF3D=.FALSE.


      if (if3d) then
         if (ly1.ne.lx1.or.lz1.ne.lx1) then
            if (nid.eq.0) write(6,13) lx1,ly1,lz1
   13       format('ERROR: lx1,ly1,lz1:',3i5,' must be equal for 3D')
            call exitt
         endif
      else
         if (ly1.ne.lx1.or.lz1.ne.1) then
            if (nid.eq.0) write(6,12) lx1,ly1,lz1
   12       format('ERROR: ',3i5,' must have lx1=ly1; lz1=1, for 2D')
            call exitt
         endif
      endif


      READ(9,*,ERR=400) NPARAM
      DO 20 I=1,NPARAM
         READ(9,*,ERR=400)PARAM(I)
   20 CONTINUE
C
C     Check to see if code is compiled with sufficient number of fields.
C
      ifmhd  = .false.
      ifessr = .false.
      if (param(29).ne.0.) ifmhd  = .true.
      if (ifmhd)           ifessr = .true.

      if (ifmvbd .and. ifsplit) then
         write(6,*) 
     $   'Moving boundary in Pn-Pn is not supported'
         call exitt
      endif

      if (ifmhd .and. lbx1.ne.lx1) then
         write(6,*) 
     $   'For MHD, need lbx1=lx1, etc.; Change SIZEu & recompile'
         call exitt
      endif

      ifpert = .false.
      if (param(31).ne.0.) ifpert = .true.
      npert = abs(param(31)) 
c
      if (ifpert .and. lpx1.ne.lx1) then
         write(6,*) 
     $   'For Lyapunov, need lpx1=lx1, etc.; Change SIZEu & recompile'
      endif

c
      NPSCAL=INT(PARAM(23))
      NPSCL1=NPSCAL+1
      if (ifmhd) npscl1 = npscl1 + 1
      NPSCL2=NPSCAL+2
      IF (NPSCL1.GT.LDIMT) THEN
         WRITE(6,21) LDIMT,NPSCL1
   21    FORMAT(//,2X,'Error: This NEKTON Solver has been compiled'
     $           /,2X,'       for',I3,' passive scalars.  This run'
     $           /,2X,'       requires that LDIMT be set to',I3,'.')
         call exitt
      ENDIF
   

c     dealiasing handling
      if (param(99).lt.0) then
         param(99) = -1       ! No  dealiasing 
      else
         param(99) = 4        ! default
      endif

      if (param(99).gt.-1 .and. lxd.le.lx1) then
          if(nid.eq.0) write(6,*) 
     &    'ABORT: LXD=<LX1, change LXD and recompile!'
          call exitt
      endif


c     I/O format handling
      if (param(67).lt.0) then
         param(67) = 0        ! ASCII
 
      else
         param(67) = 6        ! binary is default
      endif

      if (param(66).lt.0) then
         param(66) = 0        ! ASCII
 
      else
         param(66) = 6        ! binary is default
      endif

c   
c     Read in the passive scalar conduct and rhocp's:
c
c     fluid
c                 .viscosity is PARAM(2)
c                 .if it is negative, it indicates that Re has been input
c                 .therefore, redefine PARAM(2) = -1.0/PARAM(2)
c
      if(param(2) .lt.0.0) param(2)  = -1.0/param(2)
      if(param(8) .lt.0.0) param(8)  = -1.0/param(8)
      if(param(29).lt.0.0) param(29) = -1.0/param(29)
C
      CPFLD(1,1)=PARAM(2)
      CPFLD(1,2)=PARAM(1)
C     temperature
      CPFLD(2,1)=PARAM(8)
      CPFLD(2,2)=PARAM(7)
      CPFLD(2,3)=PARAM(9)
c
c     passive scalars
c
      READ(9,*,ERR=400) NSKIP
      IF (NPSCAL.GE.1 .AND. .NOT.IFUSERVP) THEN
         READ(9,*,ERR=400)(CPFLD(I,1),I=3,NPSCL2)
         IF(NPSCL2.LT.9)READ(9,*)
         READ(9,*,ERR=400)(CPFLD(I,2),I=3,NPSCL2)
         IF(NPSCL2.LT.9)READ(9,*)
         do i=3,npscl2
            if (cpfld(i,1).lt.0) cpfld(i,1) = -1./cpfld(i,1)
            if (cpfld(i,2).lt.0) cpfld(i,2) = -1./cpfld(i,2)
         enddo
      ELSE
         DO 25 I=1,NSKIP
            READ(9,*,ERR=500)
   25       CONTINUE
      ENDIF
      ifldmhd = npscal + 3
      if (ifmhd) then
         cpfld(ifldmhd,1) = param(29)  ! magnetic viscosity
         cpfld(ifldmhd,2) = param( 1)  ! magnetic rho same as for fluid
      endif
C
C     Read logical equation type descriptors....
C
      READ(9,*,ERR=500) NLOGIC
      READ(9,*,ERR=500) IFFLOW
      READ(9,*,ERR=500) IFHEAT
      READ(9,*,ERR=500) IFTRAN
C     IFADVC(1) is IFNAV; IFADVC(2) istemperature; rest are 9 passive scalars
      READ(9,*,ERR=500) (IFADVC(I),I=1,NPSCL2)
      READ(9,*,ERR=500) (IFTMSH(I),I=0,NPSCL2)
c
      IFAXIS   = .false.
      IFSTRS   = .false.
      IFLOMACH = .false.
      IFMGRID  = .false.
      IFMODEL  = .false.
      IFKEPS   = .false.
      IFMVBD   = .false.
      IFCHAR   = .false.
      ifanls   = .FALSE.
      ifanl2   = .FALSE.
c
      IF (NLOGIC.GE.6 ) READ(9,*,ERR=500) IFAXIS
      if (if3d) ifaxis = .false.
c
      IF (NLOGIC.GE.7 ) READ(9,*,ERR=500) IFSTRS
      IF (NLOGIC.GE.8 ) READ(9,*,ERR=500) IFLOMACH
      IF (NLOGIC.GE.9 ) READ(9,*,ERR=500) IFMGRID
      IF (NLOGIC.GE.10) READ(9,*,ERR=500) IFMODEL
      IF (NLOGIC.GE.11) READ(9,*,ERR=500) IFKEPS
      IF (NLOGIC.GE.12) READ(9,*,ERR=500) IFMVBD
      IF (NLOGIC.GE.13) READ(9,*,ERR=500) IFCHAR
      IF (NLOGIC.GE.14) READ(9,*,ERR=500) ifanls
      IF (NLOGIC.GE.15) READ(9,*,ERR=500) ifanl2
      MAXLOG=15
      IF(NLOGIC.GT.MAXLOG)THEN
         DO 30 IL=1,NLOGIC-MAXLOG
            READ(9,*,ERR=500)
   30    CONTINUE
      ENDIF
C
C     Set up default time dependent coefficients - NSTEPS,DT.
C
      IF (.NOT.IFTRAN) THEN
         PARAM(11) = 1.0
         PARAM(12) = 1.0
         PARAM(19) = 0.0
      ENDIF
c     IF (PARAM(22).EQ.0.0) PARAM(22)=1.0E-5
c
c     Check here for global fast diagonalization method or z-homogeneity.
c     This is here because it influence the mesh read, which follows.
c
c
      nelx   = abs(param(116))   ! check for global tensor-product structure
      nely   = abs(param(117))
      nelz   = abs(param(118))
      n_o    = 0
c     n_o    = abs(param(90))
c
      if (n_o.eq.0) then
         ifzper=.false.
         ifgfdm=.false.
         if (nelz.gt.0) ifzper=.true.
         if (nelx.gt.0) ifgfdm=.true.
         if (nelx.gt.0) ifzper=.false.
      endif


      if (iflomach .and. .not.ifsplit) then
         write(6,*) 
     $   'For low Mach, need lx2=lx1, etc.; Change SIZEu & recompile'
         call exitt
      endif

      if (iflomach .and. .not.ifheat) then
         write(6,*) 
     $   'Low Mach requires ifheat=true; Set IFHEAT true in .rea file'
         call exitt
      endif

      if (ifsplit)         ifchar = .false.   ! For now, at least.
      if (ifmhd)           ifchar = .false.   ! For now, at least.


C
C     End of PARAMETER read for all processors.
C
      return
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401)
  401 FORMAT(2X,'ERROR READING PARAMETER DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
      call exitt
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING LOGICAL DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDPARAM.')
      call exitt
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdmesh
C
C     .Read number of elements
C
C     .Construct sequential element-processor partition according
C      to number of elements and processors
C
C     .Selectively read mesh (defined by element vertices, and group numbers)
C      on each processor
C
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      character*1 adum
      real    dum(4)


c     Read elemental mesh data, formatted
      iffmtin = .true.

      call chk_nel

      NSIDES=NDIM*2
      DO 40 IEG=1,NELGT
         IF (GLLNID(IEG).EQ.NID) THEN
            IEL=GLLEL(IEG)
            READ(9,30,ERR=500,END=500) IGROUP(IEL)
   30       FORMAT(43X,I5)
C           Read Corner data
            IF(NDIM.EQ.2)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
                              call rzero (zc(1 ,iel)     ,4)
            ELSE IF(NDIM.EQ.3)THEN
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=1,4)
               READ(9,*,ERR=500,END=500) (XC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (YC(IC,IEL),IC=5,8)
               READ(9,*,ERR=500,END=500) (ZC(IC,IEL),IC=5,8)
            ENDIF
         ELSE
C           Skip over this data for element NOT on this processor
            READ(9,41,ERR=500,END=500) ADUM
C           Read Corner data
            IF(NDIM.EQ.2)THEN
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
            ELSE IF(NDIM.EQ.3)THEN
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
               READ(9,41,ERR=500,END=500) ADUM
            ENDIF
         ENDIF
   40 CONTINUE
   41 FORMAT(A1)
C
C     End of mesh read.
C
      return
C
C     Error handling:
C
  400 CONTINUE
      WRITE(6,401) 
  401 FORMAT(2X,'ERROR READING SCALE FACTORS, CHECK READ FILE'
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt
C
  500 CONTINUE
      WRITE(6,501) IEG
  501 FORMAT(2X,'ERROR READING MESH DATA NEAR ELEMENT',I7
     $    ,/,2X,'ABORTING IN ROUTINE RDMESH.')
      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine rdcurve
C
C     .Read curve side data
C
C     .Disperse curve side data to all processors according 
C      to sequential partition scheme
C
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      CHARACTER*1 ANS
C
C
C
      IF (IFFMTIN) THEN
C
C     Read formatted curve side data 
C
      READ(9,*)
      READ(9,*)NCURVE
      CALL RZERO(CURVE ,48*LELT)
      CALL BLANK(CCURVE, 8*LELT)
      IF (NCURVE.GT.0) THEN
         DO 50 ICURVE=1,NCURVE
            IF (NELGT.LT.1000) THEN
               READ(9,60,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSEIF (NELGT.LT.1000000) THEN
               READ(9,61,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSE
               READ(9,62,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ENDIF
   60       FORMAT(I3,I3 ,5G14.6,1X,A1)
   61       FORMAT(I2,I6 ,5G14.6,1X,A1)
   62       FORMAT(I2,I10,5G14.6,1X,A1)

            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               CURVE (1,IEDG,IEL)=R1
               CURVE (2,IEDG,IEL)=R2
               CURVE (3,IEDG,IEL)=R3
               CURVE (4,IEDG,IEL)=R4
               CURVE (5,IEDG,IEL)=R5
               CCURVE(  IEDG,IEL)=ANS
            ENDIF
   50    CONTINUE
      ENDIF
      return
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
      return
C
      ELSE
C
C     Read unformatted curve side data 
C
      READ(8) NCURVE
      CALL RZERO(CURVE ,48*LELT)
      CALL BLANK(CCURVE, 8*LELT)
      IF (NCURVE.GT.0) THEN
         DO 1050 ICURVE=1,NCURVE
            READ(8,ERR=1500,END=1500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               CURVE (1,IEDG,IEL)=R1
               CURVE (2,IEDG,IEL)=R2
               CURVE (3,IEDG,IEL)=R3
               CURVE (4,IEDG,IEL)=R4
               CURVE (5,IEDG,IEL)=R5
               CCURVE(  IEDG,IEL)=ANS
            ENDIF
 1050    CONTINUE
      ENDIF
      return
C
C     Error handling:
C
 1500 CONTINUE
      WRITE(6,1501)
 1501 FORMAT(2X,'ERROR READING unformatted CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      call exitt
C
      return
      ENDIF
      END
c-----------------------------------------------------------------------
      subroutine rdbdry
C
C     .Read Boundary Conditions (and connectivity data)
C
C     .Disperse boundary condition data to all processors 
C      according to sequential partition scheme
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      CHARACTER CBC1*1,CBC3*3,CHTEMP*1,CHTMP3*3
      EQUIVALENCE (CHTEMP,CHTMP3)
      character*80 string
C
C     Set up TEMPORARY value for NFIELD - NFLDT
C
      NFLDT = 1
      IF (IFHEAT) NFLDT=2+NPSCAL
      if (ifmhd ) nfldt=2+npscal+1
      NBCS      = NFLDT
      IBCS      = 2
      IF (IFFLOW) IBCS = 1
      NSIDES    = 2*NDIM
C
C     Read boundary conditions for all fields
C
      LCBC=18*LELT*(LDIMT1 + 1)
      LRBC=30*LELT*(LDIMT1 + 1)
      CALL RZERO(BC ,LRBC)
      CALL BLANK(CBC,LCBC)
C
C-----------------------------------------------------------------
C  Formatted Reads
C-----------------------------------------------------------------
C
      IF (IFFMTIN) THEN
C
      READ(9,*,ERR=500,END=500)  !   ***** BOUNDARY CONDITIONS *****
      ibcnew = 1
      DO 100 IFIELD=ibcnew,NBCS  !     DO 100 IFIELD=IBCS,NBCS
        NEL=NELGT
        if (.not.iftmsh(ifield)) nel = nelgv
C       Fluid and/or thermal
        read(9,81) string        !  ***** FLUID   BOUNDARY CONDITIONS *****
        call capit(string,80)

c       write(6,*) 'reading BC:',ifield,ibcs,nbcs
c       write(6,81) string
c       if1 = indx1(string,'NO ',3)
c       write(6,*) if1,' if NO.  quit.',ifield,ibcs,nbcs
c       write(6,*) ifield,iftmsh(ifield),nel,' iftmsh'
c       call exitt


        if (indx1(string,'NO ',3).eq.0) then ! we have acitve bc info
C
         IF(VNEKTON .LE. 2.52) NBCREA = 3
         IF(VNEKTON .GE. 2.55) NBCREA = 5
C
         DO 80 IEG=1,NEL
         DO 80 ISIDE=1,NSIDES
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               IF (NELGT.LT.1000) THEN
                  READ(9,50,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
c                 write(6,50)
c    $            CHTEMP,
c    $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
c    $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
   50             FORMAT(A1,A3,2I3,5G14.6)
               ELSEIF (NELGT.LT.100000) THEN
                  READ(9,51,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
   51             FORMAT(A1,A3,I5,I1,5G14.6)
               ELSE
                  READ(9,*,ERR=500,END=500)    
     $            CHTEMP,
     $            CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $            (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
               ENDIF
C              Mesh B.C.'s in 1st column of 1st field
               IF (CHTEMP.NE.' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
C              check for fortran function as denoted by lower case bc's:
               CBC1=CBC(ISIDE,IEL,IFIELD)
               CBC3=CBC(ISIDE,IEL,IFIELD)
               ICBC1=ICHAR(CBC1)
               IF (ICBC1.GE.97.AND.ICBC1.LE.122) THEN
                  IF(CBC3(3:3).NE.'i')NLINES=BC(1,ISIDE,IEL,IFIELD)
                  IF(CBC3(3:3).EQ.'i')NLINES=BC(4,ISIDE,IEL,IFIELD)
                  DO 60 I=1,NLINES
   60             READ(9,*,ERR=500,END=500)
               ENDIF
            ELSE
               READ(9,*,ERR=500,END=500)   cbc1  ! dummy read, pff 4/28/05
            ENDIF
   80    CONTINUE
        endif
   81   format(a80)
  100 CONTINUE
C
C     END OF BC READ
C
C     Check for dummy line:  "NO THERMAL B.C.'S"
      IF (NFLDT.EQ.1) READ(9,*,ERR=500,END=500)
C
      return
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501) IFIELD,IEG
  501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I6
     $    ,/,2X,'ABORTING IN ROUTINE RDBDRY.')
      call exitt
      return
C
C
      ELSE
C
C-----------------------------------------------------------------
C  UNformatted Reads
C-----------------------------------------------------------------
C
c     READ(8,ERR=500,END=500)
      DO 1100 IFIELD=IBCS,NBCS
         NEL=NELGT
C        Fluid and/or thermal
         NBCREA = 5
C
         DO 1080 IEG=1,NEL
         DO 1080 ISIDE=1,NSIDES
            IF (GLLNID(IEG).EQ.NID) THEN
               IEL=GLLEL(IEG)
               READ(8,ERR=1500,END=1500)    
     $         CHTMP3,
     $         CBC(ISIDE,IEL,IFIELD),ID1,ID2,
     $         (BC(II,ISIDE,IEL,IFIELD),II=1,NBCREA)
C
C              Mesh B.C.'s in 1st column of 1st field
               IF (CHTEMP.NE.' ') CBC(ISIDE,IEL,0)(1:1)= CHTEMP
C              check for fortran function as denoted by lower case bc's:
            ELSE
               IEL=1
               READ(8,ERR=1500,END=1500) CHTMP3,
     $         CBCS(ISIDE,IEL),ID1,ID2,(BCS(II,ISIDE,IEL),II=1,NBCREA)
C              check for fortran function as denoted by lower case bc's:
            ENDIF
 1080    CONTINUE
 1100 CONTINUE
C
C     END OF BC READ
C
      return
C
C     Error handling:
C
 1500 CONTINUE
      WRITE(6,1501) IFIELD,IEG
 1501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I6
     $    ,/,2X,'(unformatted) ABORTING IN ROUTINE RDBDRY.')
      call exitt
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdicdf
C
C     .Read Initial Conditions / Drive Force
C
C     .Broadcast ICFILE to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER*80 LINE
      LOGICAL      IFGTIL
C
C     Read Initial Conditions/Restart Files
      CALL BLANK(INITC,1200)
      READ(9,80,ERR=200,END=200) LINE
      IF (INDX1(LINE,'RESTART',7).NE.0) THEN
         IF (.NOT.IFGTIL(NSKIP,LINE)) GOTO 200
c        READ(LINE,*,ERR=200,END=200) NSKIP
         DO 50 I=1,NSKIP
            READ(9,80,ERR=200,END=200) INITC(I)
   50    CONTINUE
         READ(9,80,ERR=200,END=200) LINE
      ENDIF
   80 FORMAT(A80)
      IF (.NOT.IFGTIL(NSKIP,LINE)) GOTO 200
c     READ(LINE,*,ERR=200,END=200)NSKIP
      DO 100 I=1,NSKIP
         READ(9,80,ERR=200,END=200) LINE
  100 CONTINUE
C     Read drive force data
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NSKIP
      DO 110 I=1,NSKIP
        READ(9,80,ERR=200,END=200) LINE
  110 CONTINUE
      return
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING INITIAL CONDITION/DRIVE FORCE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDICDF.')
      call exitt
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdmatp
C
C     .Read materials property data
C
C     .Disperse material properties to all processors according 
C      to sequential partition scheme
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      CHARACTER*80 LINE
C
      CALL IZERO(MATYPE,16*LDIMT1)
      CALL RZERO(CPGRP ,48*LDIMT1)
C
C     Read material property data
C
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NSKIP
      READ(9,*,ERR=200,END=200) NPACKS
      DO 100 IIG=1,NPACKS
         IFVPS=.TRUE.
         READ(9,*)IGRP,IFLD,ITYPE
         MATYPE(IGRP,IFLD)=ITYPE
         DO 100 IPROP=1,3
            IF(ITYPE.EQ.1) READ(9,* ) CPGRP(IGRP,IFLD,IPROP)
            IF(ITYPE.EQ.2) READ(9,80) LINE
   80 FORMAT(A80)
  100 CONTINUE
C
      return
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING MATERIAL PROPERTIES DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDMATP.')
      call exitt
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdhist
C
C     .Read history data
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
C
      CALL BLANK (HCODE ,11*lhis)
      CALL IZERO (LOCHIS, 4*lhis)
C
C     Read history data
      READ (9,*)
      READ (9,*,ERR=200,END=200) NHIS
      if (nhis.gt.lhis) then
         write(6,*) nid,' Too many history pts. RESET LHIS.',nhis,lhis
         call exitt
      endif
c
C     HCODE(10) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
      if (nhis.gt.0) then
         do i=1,nhis
            if (nelgt.lt.100000) then
               read(9,130,err=200,end=200)
     $         (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
  130          format(1x,11a1,1x,4i5)
            else
               read(9,131,err=200,end=200)
     $         (hcode(ii,i),ii=1,11),(lochis(i2,i),i2=1,4)
  131          format(1x,11a1,1x,3i5,i10)
            endif
c
c           threshold lochis locations to allow easy specification of "NX,NY,NZ"
c           pff 1/7/97
c
            if (hcode(10,i).eq.'H') then
               lochis(1,i) = min(lochis(1,i),nx1)
               lochis(2,i) = min(lochis(2,i),ny1)
               lochis(3,i) = min(lochis(3,i),nz1)
c
c              if lochis_k = -1, set it to nxk/2   pff 8/21/03
c
               if (lochis(1,i).eq.-1) lochis(1,i) = (nx1+1)/2
               if (lochis(2,i).eq.-1) lochis(2,i) = (ny1+1)/2
               if (lochis(3,i).eq.-1) lochis(3,i) = (nz1+1)/2
            endif
         enddo
      endif
C
      return
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING HISTORY DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDHIST.')
      call exitt
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdout
C
C     .Read output specs
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
C
C     Read output specs
      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOUTS
      READ(9,*,ERR=200,END=200) IFXYO
      READ(9,*,ERR=200,END=200) IFVO
      READ(9,*,ERR=200,END=200) IFPO
      READ(9,*,ERR=200,END=200) IFTO
      READ(9,*,ERR=200,END=200) IFBO   !  IFTGO
      READ(9,*,ERR=200,END=200) IPSCO
      IF (IPSCO.GT.0) THEN
         IF (IPSCO.GT.NPSCAL) GOTO 200
         DO 120 I=1,IPSCO
            READ(9,*,ERR=200,END=200) IFPSCO(I)
  120    CONTINUE
      ENDIF
C
      return
C
C     Error handling:
C
  200 CONTINUE
      WRITE(6,201)
  201 FORMAT(2X,'ERROR READING OUTPUT SPECIFICATION DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDOUT.')
      call exitt
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rdobj
C
C     .Read objects
C
C     .Broadcast to all processors
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
C
C     Default if no data is read No Objects
C
      NOBJ=0

      READ(9,*,ERR=200,END=200)
      READ(9,*,ERR=200,END=200) NOBJ
C
      IF(NOBJ.GT.MAXOBJ)THEN
         write(6,*) nid,'ERROR, too many objects:',nobj,maxobj
         call exitt
      ENDIF
C
      DO 10 IOBJ = 1,NOBJ
         READ(9,*,ERR=200,END=200) NMEMBER(IOBJ)
         IF(NMEMBER(IOBJ).GT.MAXMBR)THEN
            PRINT*,'ERROR: Too many members in object ',IOBJ
            call exitt
         ENDIF
         DO 5 MEMBER=1,NMEMBER(IOBJ)
            READ(9,*,ERR=200,END=200) OBJECT(IOBJ,MEMBER,1),
     $                                OBJECT(IOBJ,MEMBER,2)
    5    CONTINUE
   10 CONTINUE
      if (nid.eq.0) write(6,*) nobj,' objects found,'
     $                             ,(nmember(k),k=1,nobj)
C
      return
C
C     Error handling:  For old versions, default to no objects
C
  200 CONTINUE
      NOBJ=0
C
      return
      END
c-----------------------------------------------------------------------
      subroutine vrdsmsh
C
C=====================================================================
C     Verify that mesh and dssum are properly defined by performing
C        a direct stiffness operation on the X,Y and Z coordinates.
C     Note that periodic faces are not checked here.
C=====================================================================
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      COMMON /SCRNS/ TA(LX1,LY1,LZ1,LELT),TB(LX1,LY1,LZ1,LELT)
     $           ,QMASK(LX1,LY1,LZ1,LELT),tmp(2)
      CHARACTER*3 CB

c      call  vrdsmshx  ! verify mesh topology

      if(nid.eq.0) write(*,*) 'verify mesh topology'

      IERR      = 0
      EPS       = 1.0e-04
      EPS       = 1.0e-03
      IFIELD    = 1
      IF (IFHEAT) IFIELD = 2
      NXYZ1     = NX1*NY1*NZ1
      NTOT      = NX1*NY1*NZ1*NELT
      NFACES    = 2*NDIM

      xmx = glmax(xm1,ntot)
      xmn = glmin(xm1,ntot)
      ymx = glmax(ym1,ntot)
      ymn = glmin(ym1,ntot)
      zmx = glmax(zm1,ntot)
      zmn = glmin(zm1,ntot)
      if (nid.eq.0) write(6,*) xmn,xmx,' Xrange'
      if (nid.eq.0) write(6,*) ymn,ymx,' Yrange'
      if (nid.eq.0) write(6,*) zmn,zmx,' Zrange'
c     return
C
C     First check - use 1/Multiplicity
C
      IF (IFHEAT) THEN 
         CALL COPY(TA,TMULT,NTOT)
      ELSE
         CALL COPY(TA,VMULT,NTOT)
      ENDIF
c
c     write(6,1) 
c    $(nid,'tab4',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c   1 format(i3,a4,i3,16f5.2)
c
      CALL DSSUM(TA,NX1,NY1,NZ1)
c
c     write(6,1) 
c    $(nid,'taaf',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c
      CALL RONE (TB,NTOT)
      CALL SUB2 (TB,TA,NTOT)
      DO 1000 IE=1,NELT
      IEG=LGLEL(IE,NODE)
      DO 1000 IZ=1,NZ1
      DO 1000 IY=1,NY1
      DO 1000 IX=1,NX1
         IF (ABS(TB(IX,IY,IZ,IE)).GT.EPS ) THEN
            WRITE(6,1005) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TA(IX,IY,IZ,IE),eps
c           WRITE(7,1005) IX,IY,IZ,IEG
c    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
c    $      ,QMASK(IX,IY,IZ,IE)
 1005       FORMAT(2X,'WARNING: DSSUM problem at:',/
     $            ,1X,'I,J,K,IE:',3I5,i9,/
     $            ,2X,'Near X =',3G16.8,', d:',2G16.8)
            IERR=4
         ENDIF
 1000 CONTINUE
C
C     Set up QMASK quickly to annihilate checks on periodic bc's
C
      CALL RONE(QMASK,NTOT)
      DO 100 IEL=1,NELT
      DO 100 IFACE=1,NFACES
         CB =CBC(IFACE,IEL,IFIELD)
         IF (CB.EQ.'P  ') 
     $         CALL FACEV(QMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
  100 CONTINUE
      CALL DSOP(QMASK,'MUL',NX1,NY1,NZ1)

c      xxmin = glmin(xm1,ntot)
c      yymin = glmin(ym1,ntot)
c      zzmin = glmin(zm1,ntot)
c      xxmax = glmax(xm1,ntot)
c      yymax = glmax(ym1,ntot)
c      zzmax = glmax(zm1,ntot)
c      if (nid.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
c    7 format('xyz minmx2:',6g13.5)



C
C     X-component
C
      CALL COPY(TA,XM1,NTOT)
      CALL COPY(TB,XM1,NTOT)
      CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
      CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
      CALL SUB2(TA,XM1,NTOT)
      CALL SUB2(TB,XM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1100 IE=1,NELT
         XSCMAX = VLMAX(XM1(1,1,1,IE),NXYZ1)
         XSCMIN = VLMIN(XM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(XSCMAX-XSCMIN)
         SCAL2=ABS(XSCMAX)
         SCAL3=ABS(XSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         XSCALE = 1./SCAL1
         IEG=LGLEL(IE,NODE)
         DO 1100 IZ=1,NZ1
         DO 1100 IY=1,NY1
         DO 1100 IX=1,NX1
         if (abs(ta(ix,iy,iz,ie)*xscale).gt.eps .or.
     $       abs(tb(ix,iy,iz,ie)*xscale).gt.eps ) then
            write(6,1105) ix,iy,iz,ieg
     $      ,xm1(ix,iy,iz,ie),ym1(ix,iy,iz,ie),zm1(ix,iy,iz,ie)
     $      ,tb(ix,iy,iz,ie),ta(ix,iy,iz,ie),XSCALE
 1105       format(1x,'WARNING1 Element mesh mismatch at:',/
     $            ,1x,'i,j,k,ie:',3i5,I9,/
     $            ,1X,'Near X =',3G16.8,', d:',3G16.8)
            ierr=1
            call exitt
         endif
 1100 CONTINUE
C
C     Y-component
C
      CALL COPY(TA,YM1,NTOT)
      CALL COPY(TB,YM1,NTOT)
      CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
      CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
      CALL SUB2(TA,YM1,NTOT)
      CALL SUB2(TB,YM1,NTOT)
      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)
      DO 1200 IE=1,NELT
         YSCMAX = VLMAX(YM1(1,1,1,IE),NXYZ1)
         YSCMIN = VLMIN(YM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(YSCMAX-YSCMIN)
         SCAL2=ABS(YSCMAX)
         SCAL3=ABS(YSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         YSCALE = 1./SCAL1
         IEG=LGLEL(IE,NODE)
         DO 1200 IZ=1,NZ1
         DO 1200 IY=1,NY1
         DO 1200 IX=1,NX1
         IF (ABS(TA(IX,IY,IZ,IE)*YSCALE).GT.EPS .OR.
     $       ABS(TB(IX,IY,IZ,IE)*YSCALE).GT.EPS ) THEN
            WRITE(6,1205) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),yscale
 1205       FORMAT(1X,'WARNING2 Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',3I5,i9,/
     $            ,1X,'Near Y =',3G16.8,', d:',3G16.8)
            IERR=2
            call exitt
         ENDIF
 1200 CONTINUE
C
C     Z-component
C
      IF (IF3D) THEN
       CALL COPY(TA,ZM1,NTOT)
       CALL COPY(TB,ZM1,NTOT)
       CALL DSOP(TA,'MIN',NX1,NY1,NZ1)
       CALL DSOP(TB,'MAX',NX1,NY1,NZ1)
       CALL SUB2(TA,ZM1,NTOT)
       CALL SUB2(TB,ZM1,NTOT)
       CALL COL2(TA,QMASK,NTOT)
       CALL COL2(TB,QMASK,NTOT)
       DO 1300 IE=1,NELT
          ZSCMAX = VLMAX(ZM1(1,1,1,IE),NXYZ1)
          ZSCMIN = VLMIN(ZM1(1,1,1,IE),NXYZ1)
          SCAL1=ABS(ZSCMAX-ZSCMIN)
          SCAL2=ABS(ZSCMAX)
          SCAL3=ABS(ZSCMIN)
          SCAL1=MAX(SCAL1,SCAL2)
          SCAL1=MAX(SCAL1,SCAL3)
          ZSCALE = 1./SCAL1
          IEG=LGLEL(IE,NODE)
          DO 1300 IZ=1,NZ1
          DO 1300 IY=1,NY1
          DO 1300 IX=1,NX1
          IF (ABS(TA(IX,IY,IZ,IE)*ZSCALE).GT.EPS .OR.
     $        ABS(TB(IX,IY,IZ,IE)*ZSCALE).GT.EPS ) THEN
           WRITE(6,1305) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE),zscale
 1305       FORMAT(1X,'WARNING3 Element mesh mismatch at:',/
     $            ,1X,'I,J,K,IE:',3I5,i9,/
     $            ,1X,'Near Z =',3G16.8,', d:',3G16.8)
            IERR=3
            call exitt
          ENDIF
 1300  CONTINUE
      ENDIF
C
      IF (IERR.gt.0) THEN
         WRITE(6,1400) 
 1400    FORMAT
     $   (' Mesh consistency check failed.  EXITING in VRDSMSH.')
            call exitt
      ENDIF
C
      tmp(1)=ierr
      CALL GOP(tmp,tmp(2),'M  ',1)
      IF (tmp(1).ge.4.0) THEN
         WRITE(6,1400) 
     $   (' Mesh consistency check failed.  EXITING in VRDSMSH.')
         call exitt
      ENDIF
C
      if(nid.eq.0) then
        write(6,*) 'done :: verify mesh topology'
        write(6,*) ''
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine vrdsmshx  ! verify mesh topology
C
C=====================================================================
C     Verify that mesh and dssum are properly defined by performing
C        a direct stiffness operation on the X,Y and Z coordinates.
C     Note that periodic faces are not checked here.
C=====================================================================
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      common /scrns/ tc(lx1,ly1,lz1,lelt),td(lx1,ly1,lz1,lelt)
     $             , ta(lx1,ly1,lz1,lelt),tb(lx1,ly1,lz1,lelt)
     $             , qmask(lx1,ly1,lz1,lelt)
      CHARACTER*3 CB
C
      IERR      = 0
      EPS       = 1.0e-04
      EPS       = 1.0e-03
      IFIELD    = 1
      IF (IFHEAT) IFIELD = 2
      NXYZ1     = NX1*NY1*NZ1
      NTOT      = NX1*NY1*NZ1*NELT
      NFACES    = 2*NDIM

      xmx = glmax(xm1,ntot)
      xmn = glmin(xm1,ntot)
      ymx = glmax(ym1,ntot)
      ymn = glmin(ym1,ntot)
      zmx = glmax(zm1,ntot)
      zmn = glmin(zm1,ntot)
      if (nid.eq.0) write(6,*) xmn,xmx,' Xrange'
      if (nid.eq.0) write(6,*) ymn,ymx,' Yrange'
      if (nid.eq.0) write(6,*) zmn,zmx,' Zrange'
c     return
C
C     First check - use 1/Multiplicity
C
      IF (IFHEAT) THEN 
         CALL COPY(TA,TMULT,NTOT)
      ELSE
         CALL COPY(TA,VMULT,NTOT)
      ENDIF
c
c     write(6,1) 
c    $(nid,'tab4',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c   1 format(i3,a4,i3,16f5.2)
c
      CALL DSSUM(TA,NX1,NY1,NZ1)
c
c     write(6,1) 
c    $(nid,'taaf',lglel(ie,node),(ta(k,1,1,ie),k=1,nx1*ny1),ie=1,nelt)
c
      CALL RONE (TB,NTOT)
      CALL SUB2 (TB,TA,NTOT)
      DO 1000 IE=1,NELT
      IEG=LGLEL(IE,NODE)
      DO 1000 IZ=1,NZ1
      DO 1000 IY=1,NY1
      DO 1000 IX=1,NX1
         IF (ABS(TB(IX,IY,IZ,IE)).GT.EPS ) THEN
            WRITE(6,1005) IX,IY,IZ,IEG
     $      ,XM1(IX,IY,IZ,IE),YM1(IX,IY,IZ,IE),ZM1(IX,IY,IZ,IE)
     $      ,TA(IX,IY,IZ,IE),eps
c           WRITE(7,1005) IX,IY,IZ,IEG
c    $      ,XM1(IX,IY,IZ,IE),TB(IX,IY,IZ,IE),TA(IX,IY,IZ,IE)
c    $      ,QMASK(IX,IY,IZ,IE)
 1005       FORMAT(2X,'WARNING: DSSUM problem at:',/
     $            ,1X,'I,J,K,IE:',3I5,i9,/
     $            ,2X,'Near X =',3G16.8,', d:',2G16.8)
            IERR=4
         ENDIF
 1000 CONTINUE
C
C     Set up QMASK quickly to annihilate checks on periodic bc's
C
      CALL RONE(QMASK,NTOT)
      DO 100 IEL=1,NELT
      DO 100 IFACE=1,NFACES
         CB =CBC(IFACE,IEL,IFIELD)
         IF (CB.EQ.'P  ') 
     $         CALL FACEV(QMASK,IEL,IFACE,0.0,NX1,NY1,NZ1)
  100 CONTINUE
      call dsop(QMASK,'MUL',NX1,NY1,NZ1)

      xxmin = glmin(xm1,ntot)
      yymin = glmin(ym1,ntot)
      zzmin = glmin(zm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ym1,ntot)
      zzmax = glmax(zm1,ntot)
      if (nid.eq.0) write(6,7) xxmin,yymin,zzmin,xxmax,yymax,zzmax
    7 format('xyz minmx2:',6g13.5)


C
C     X-component
C
      call copy(ta,xm1,ntot)
      call copy(tb,xm1,ntot)
      call dsop(ta,'min',nx1,ny1,nz1)
      call dsop(tb,'max',nx1,ny1,nz1)

      call copy(tc,xm1,ntot)
      call copy(td,xm1,ntot)
      call dsop(tc,'min',nx1,ny1,nz1)
      call dsop(td,'max',nx1,ny1,nz1)

      xxmin = glmin(xm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nid.eq.0) write(6,9) xxmin,yymin,zzmin,xxmax,yymax,zzmax
    9 format('xyz minmx3:',6g13.5)

      CALL SUB2(TA,XM1,NTOT)
      CALL SUB2(TB,XM1,NTOT)

      xxmin = glmin(qmask,ntot)
      xxmax = glmax(qmask,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nid.eq.0) write(6,19) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   19 format('xyz minmx4:',6g13.5)

      CALL COL2(TA,QMASK,NTOT)
      CALL COL2(TB,QMASK,NTOT)

      xxmin = glmin(qmask,ntot)
      xxmax = glmax(qmask,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nid.eq.0) write(6,29) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   29 format('xyz minmx5:',6g13.5)

      DO 1100 IE=1,NELT
         XSCMAX = VLMAX(XM1(1,1,1,IE),NXYZ1)
         XSCMIN = VLMIN(XM1(1,1,1,IE),NXYZ1)
         SCAL1=ABS(XSCMAX-XSCMIN)
         SCAL2=ABS(XSCMAX)
         SCAL3=ABS(XSCMIN)
         SCAL1=MAX(SCAL1,SCAL2)
         SCAL1=MAX(SCAL1,SCAL3)
         XSCALE = 1./SCAL1
         IEG=LGLEL(IE,NODE)
         DO 1100 IZ=1,NZ1
         DO 1100 IY=1,NY1
         DO 1100 IX=1,NX1
         if (abs(ta(ix,iy,iz,ie)*xscale).gt.eps .or.
     $       abs(tb(ix,iy,iz,ie)*xscale).gt.eps ) then
            write(6,1105) nid,ix,iy,iz,ie,ieg
     $      ,xm1(ix,iy,iz,ie),tc(ix,iy,iz,ie),td(ix,iy,iz,ie)
     $      ,ym1(ix,iy,iz,ie),zm1(ix,iy,iz,ie)
     $      ,ta(ix,iy,iz,ie),tb(ix,iy,iz,ie),xscale
     $      ,qmask(ix,iy,iz,ie)
 1105       format(i4.4,1x,'ie:',3i3,i4,i6,1p9e11.3)
c1105       format(i4.4,1x,'ie:',3i3,i6,1p9e11.3)
            ierr=1
            goto 1101
         endif
 1100 CONTINUE
 1101 CONTINUE

      xxmin = glmin(xm1,ntot)
      xxmax = glmax(xm1,ntot)
      yymax = glmax(ta ,ntot)
      yymin = glmin(ta ,ntot)
      zzmin = glmin(tb ,ntot)
      zzmax = glmax(tb ,ntot)
      if (nid.eq.0) write(6,39) xxmin,yymin,zzmin,xxmax,yymax,zzmax
   39 format('xyz minmx5:',6g13.5)

c     ifvo = .true.
c     ifpo = .false.
c     ifto = .true.
c     call outpost(xm1,ta,tb,pr,qmask,'   ')
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine rotat2(xyz,angle,npts)
C
C     Rotate NPTS through ANGLE (in two directions IF3D).
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZ(3,1)
      COMMON /CTMP0/ RMTRX(3,3),RX(3,3),RZ(3,3),XYZN(3,10)
C
      SINA=SIN(ANGLE)
      COSA=COS(ANGLE)
      CALL RZERO(RX,9)
      CALL RZERO(RZ,9)
      RX(1,1)=COSA
      RX(2,2)=COSA
      RX(1,2)=SINA
      RX(2,1)=-SINA
      RX(3,3)=1.0
      IF (IF3D) THEN
         RZ(1,1)=COSA
         RZ(3,3)=COSA
         RZ(1,3)=SINA
         RZ(3,1)=-SINA
         RZ(2,2)=1.0
      ELSE
         RZ(1,1)=1.0
         RZ(2,2)=1.0
         RZ(3,3)=1.0
      ENDIF
      CALL MXM(RX,3,RZ,3,RMTRX,3)
C
C     Strip mine mxms in chunks of 10:
      DO 100 I=1,NPTS-10,10
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,10)
         CALL COPY(XYZ(1,I),XYZN,30)
  100 CONTINUE
      N10=MOD1(NPTS,10)
      I=NPTS-N10+1
      CALL RZERO(XYZN,30)
      IF (N10.GT.0) THEN
         CALL MXM(RMTRX,3,XYZ(1,I),3,XYZN,N10)
         CALL COPY(XYZ(1,I),XYZN,3*N10)
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine scale(xyzl,nl)
C
C     Rescale XYZL such that the mean value of IXX=IYY=IZZ for each element.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZL(3,8,LELT)
      COMMON /CTMP0/ VO(LELT),XYZI(3,LELT),CG(3,LELT)
     $              ,TI(6),WORK(6)
C
C     Compute volumes -
C
      CALL VOLUME2(VO,XYZL,NL)
      VTOT=GLSUM (VO,NL)
C
C     Compute (weighted) average inertia for each element.
C
      NCRNR=2**NDIM
      CALL RZERO(TI,6)
      DO 100 IL=1,NL
         VO0 = VO(IL)/VTOT
         CALL INRTIA(XYZI(1,IL),CG(1,IL),XYZL(1,1,IL),NCRNR,1)
         TI(1)=TI(1)+XYZI(1,IL)*VO0
         TI(2)=TI(2)+XYZI(2,IL)*VO0
         TI(3)=TI(3)+XYZI(3,IL)*VO0
         TI(4)=TI(4)+CG(1,IL)  *VO0
         TI(5)=TI(5)+CG(2,IL)  *VO0
         TI(6)=TI(6)+CG(3,IL)  *VO0
  100 CONTINUE
      CALL GOP(TI,WORK,'+  ',6)
      XI  =SQRT(TI(1))
      YI  =SQRT(TI(2))
      ZI  =1.0
      IF (IF3D) ZI=SQRT(TI(3))
C
C     Rescale ( & shift to a nearly mean zero )
C
      DO 200 IL=1,NL
      DO 200 IC=1,NCRNR
         XYZL(1,IC,IL)=(XYZL(1,IC,IL)-TI(4))/XI
         XYZL(2,IC,IL)=(XYZL(2,IC,IL)-TI(5))/YI
         XYZL(3,IC,IL)=(XYZL(3,IC,IL)-TI(6))/ZI
  200 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine inrtia(xyzi,cg,xyzl,n,itype)
C
C     Compute cg and inertia for a collection of unit point masses.
C     This is a global (multiprocessor) operation, only IF itype=2.
C
      DIMENSION XYZI(3),CG(3),XYZL(3,1)
      DIMENSION TI(4),WORK(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      TI(4)=N
      DO 100 I=1,N
         TI(1)=TI(1)+XYZL(1,I)
         TI(2)=TI(2)+XYZL(2,I)
         TI(3)=TI(3)+XYZL(3,I)
  100 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',4)
      IF (TI(4).EQ.0.0) TI(4)=1.0
      CG(1)=TI(1)/TI(4)
      CG(2)=TI(2)/TI(4)
      CG(3)=TI(3)/TI(4)
C
      TI(1)=0.0
      TI(2)=0.0
      TI(3)=0.0
      DO 200 I=1,N
         TI(1)=TI(1)+( XYZL(1,I)-CG(1) )**2
         TI(2)=TI(2)+( XYZL(2,I)-CG(2) )**2
         TI(3)=TI(3)+( XYZL(3,I)-CG(3) )**2
  200 CONTINUE
      IF (ITYPE.EQ.2) CALL GOP(TI,WORK,'+  ',3)
      TI(1)=TI(1)/TI(4)
      TI(2)=TI(2)/TI(4)
      TI(3)=TI(3)/TI(4)
      IF (ITYPE.EQ.2) THEN
C        std. def'n of inertia.
         XYZI(1)=TI(2)+TI(3)
         XYZI(2)=TI(3)+TI(1)
         XYZI(3)=TI(1)+TI(2)
      ELSE
         XYZI(1)=TI(1)
         XYZI(2)=TI(2)
         XYZI(3)=TI(3)
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine volume2(vol,xyz,n)
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      DIMENSION XYZ(3,2,2,2,1)
      DIMENSION VOL(1)
C
      DO 1000 IE=1,N
         VOL(IE)=0.0
         IF (IF3D) THEN        
           DO 20 K=1,2
           DO 20 J=1,2
           DO 20 I=1,2
              VOL1 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL2 = (XYZ(1,2,J,K,IE)-XYZ(1,1,J,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL3 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,2,J,K,IE)-XYZ(2,1,J,K,IE))
     $             * (XYZ(3,I,J,2,IE)-XYZ(3,I,J,1,IE))
              VOL4 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,I,2,K,IE)-XYZ(3,I,1,K,IE))
              VOL5 = (XYZ(1,I,2,K,IE)-XYZ(1,I,1,K,IE))
     $             * (XYZ(2,I,J,2,IE)-XYZ(2,I,J,1,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL6 = (XYZ(1,I,J,2,IE)-XYZ(1,I,J,1,IE))
     $             * (XYZ(2,I,2,K,IE)-XYZ(2,I,1,K,IE))
     $             * (XYZ(3,2,J,K,IE)-XYZ(3,1,J,K,IE))
              VOL(IE) = VOL(IE)+VOL1+VOL2+VOL3+VOL4+VOL5+VOL6
   20      CONTINUE
           VOL(IE)=VOL(IE)/8.0
         ELSE
C     2-D:
            DO 40 J=1,2
            DO 40 I=1,2
              VOL1 = (XYZ(1,2,J,1,IE)-XYZ(1,1,J,1,IE))
     $             * (XYZ(2,I,2,1,IE)-XYZ(2,I,1,1,IE))
              VOL3 = (XYZ(1,I,2,1,IE)-XYZ(1,I,1,1,IE))
     $             * (XYZ(2,2,J,1,IE)-XYZ(2,1,J,1,IE))
              VOL(IE)=VOL(IE)+VOL1+VOL3
   40      CONTINUE
           VOL(IE)=VOL(IE)/4.0
         ENDIF
         VOL(IE)=ABS(VOL(IE))
 1000 CONTINUE
C
      return
      END
c-----------------------------------------------------------------------
      subroutine findcg(cg,xyz,n)
C
C     Compute cg for N elements.
C
      INCLUDE 'SIZE'
      DIMENSION CG(3,1),XYZ(3,8,1)
C
      NCRNR=2**NDIM
      CALL RZERO(CG,3*N)
      DO 100 I =1,N
      DO 100 IC=1,NCRNR
         CG(1,I)=CG(1,I)+XYZ(1,IC,I)
         CG(2,I)=CG(2,I)+XYZ(2,IC,I)
         CG(3,I)=CG(3,I)+XYZ(3,IC,I)
  100 CONTINUE
      TMP=1.0/(NCRNR)
      CALL CMULT(CG,TMP,3*N)
      return
      END
c-----------------------------------------------------------------------
      subroutine divide(list1,list2,nl1,nl2,ifok,list,nl,xyzi,cg,WGT)
C
C     Divide the elements associated with this subdomain according to
C     the direction having the smallest moment of inertia (the "long"
C     direction).
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'TSTEP'
C
      DIMENSION LIST(LELT),LIST1(LELT),LIST2(LELT)
      DIMENSION XYZI(3),CG(3,LELT),wgt(1)
      COMMON /CTMP0/ XCG(LELT),YCG(LELT),ZCG(LELT)
      REAL IXX,IYY,IZZ
      INTEGER WORK(2),WRK2(2)
      LOGICAL IFOK
C
C     Choose "long" direction:
C
      IXX=XYZI(1)
      IYY=XYZI(2)
      IZZ=XYZI(3)
      IF (IF3D) THEN
         IF (IXX.LE.IYY.AND.IXX.LE.IZZ) THEN
            DO 104 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
               ZCG(IE)=CG(3,IE)
  104       CONTINUE
         ELSEIF (IYY.LE.IXX.AND.IYY.LE.IZZ) THEN
            DO 106 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(3,IE)
               ZCG(IE)=CG(1,IE)
  106       CONTINUE
         ELSEIF (IZZ.LE.IXX.AND.IZZ.LE.IYY) THEN
            DO 108 IE=1,NL
               XCG(IE)=CG(3,IE)
               YCG(IE)=CG(1,IE)
               ZCG(IE)=CG(2,IE)
  108       CONTINUE
         ENDIF
      ELSE
C     2-D:
         IF (IXX.LE.IYY) THEN
            DO 114 IE=1,NL
               XCG(IE)=CG(1,IE)
               YCG(IE)=CG(2,IE)
  114       CONTINUE
         ELSE
            DO 116 IE=1,NL
               XCG(IE)=CG(2,IE)
               YCG(IE)=CG(1,IE)
  116       CONTINUE
         ENDIF
      ENDIF
      call col2(xcg,wgt,nl)
      call col2(ycg,wgt,nl)
      call col2(zcg,wgt,nl)
C
C     Find median value of CG to determine dividing point:
C
      XM=FMDIAN(XCG,NL,IFOK)
      YM=FMDIAN(YCG,NL,IFOK)
      ZM=0.0
      IF (IF3D) ZM=FMDIAN(ZCG,NL,IFOK)
C
C     Diagnostics
C
      IF (.NOT.IFOK) THEN
         WRITE(6,130) NID,NL,XM,YM,ZM
         DO 120 IL=1,NL
            WRITE(6,135) NID,IL,XCG(IL),YCG(IL),ZCG(IL)
  120    CONTINUE
  130    FORMAT(I3,'DIVIDE: NL,XM,YM,ZM',I3,3F12.5)
  135    FORMAT(I3,'DIVIDE: NID,IL,XC,YC,ZCG',I4,3F12.5)
      ENDIF
C
C=============================================================
C     Divide LIST into LIST1 (XCG < XM) and LIST2 (XCG>XM).
C=============================================================
C
      NL1=0
      NL2=0
      DO 200 IE=1,NL
         IF (XCG(IE).LT.XM) THEN
            NL1=NL1+1
            LIST1(NL1)=LIST(IE)
         ENDIF
         IF (XCG(IE).GT.XM) THEN
            NL2=NL2+1
            LIST2(NL2)=LIST(IE)
         ENDIF
         IF (XCG(IE).EQ.XM) THEN
C
C           We have to look at the other directions to arrive at
C           a unique subdivision algortithm.
C
C
C           More Diagnostics
C
            IF (.NOT.IFOK) WRITE(6,201) NID,IE,XCG(IE),XM
  201    FORMAT(I3,'DIVIDE: IE,XCG,XM:',I4,3F12.5)
C
            IF (YCG(IE).LT.YM) THEN
               NL1=NL1+1
               LIST1(NL1)=LIST(IE)
            ENDIF
            IF (YCG(IE).GT.YM) THEN
               NL2=NL2+1
               LIST2(NL2)=LIST(IE)
            ENDIF
            IF (YCG(IE).EQ.YM) THEN
C              look at 3rd direction.  
               IF (IF3D .AND. ZCG(IE).LT.ZM) THEN
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ELSE IF (IF3D .AND. ZCG(IE).GT.ZM) THEN
                  NL2=NL2+1
                  LIST2(NL2)=LIST(IE)
               ELSE 
C                 for 2- or 3-D intdeterminate case:
                  NL1=NL1+1
                  LIST1(NL1)=LIST(IE)
               ENDIF
            ENDIF
C
         ENDIF
  200 CONTINUE
C
C     Check for an even distribution (i.e. - not different by
C     more than 1):
C
      IFOK=.TRUE.  
      WORK(1)=NL1 
      WORK(2)=NL2
      CALL IGOP(WORK,WRK2,'+  ',2)
      IF (ABS(WORK(1)-WORK(2)).GT.1) IFOK=.FALSE.
C
      return
      END
c-----------------------------------------------------------------------
      subroutine gray(ipg,ipg1,itmp,n)
      DIMENSION IPG(1),IPG1(1),ITMP(1)
      INTEGER CD
C
      CD=LOG2(N)
      NPR = 1
      IPG(1) = 1
      DO 100 I=1,CD
         N1  = NPR+1
         N2  = NPR*2
         CALL ICOPY(IPG(N1),IPG(1),NPR)
         CALL IADD (IPG(N1),NPR,NPR)
         CALL IFLIP(IPG(N1),NPR)
         NPR = 2*NPR
  100 CONTINUE
      N1=NPR+1
c     IPG(0)=IPG(NPR)
c     IPG(N1)=IPG(1)
C
C     Compute the inverse of the gray code mapping.  In other words,
C     given the ith entry on the gray code list, IPG(i) returns the
C     appropriate Node number (=nid+1).  Given a node number, NODE,
C     IPG1(NODE) returns the location in the gray code (i.e. ring) map.
C
      CALL ICOPY(ITMP,IPG ,NPR)
      CALL ISORT(ITMP,IPG1,NPR)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine bin_rd1  ! read mesh, curve, and bc info

      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      logical ifbswap

      if(nid.eq.0) write(6,*) 'read .re2 file'
   
      etime1 = dnekclock()

                  ibc = 2
      if (ifflow) ibc = 1

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal
      if (ifmhd ) nfldt = 2+npscal+1

c
c     If p32 = 0.1, there will be no bcs read in
c
      if (param(32).gt.0) nfldt = ibc + param(32)-1


      lcbc=18*lelt*(ldimt1 + 1)
      call blank(cbc,lcbc)

      call open_bin_file (ifbswap)

      if (nid.eq.0) write(6,*)    '  read mesh '
      call bin_rd1_mesh  (ifbswap)   ! version 1 of binary reader
      if (nid.eq.0) write(6,*) '  read curved sides '
      call bin_rd1_curve (ifbswap)

      do ifield = ibc,nfldt
         if (nid.eq.0) write(6,*) '  read bc ifld',ifield
         call bin_rd1_bc (cbc(1,1,ifield),bc(1,1,1,ifield),ifbswap)
      enddo

      call close_bin_file

      if(nid.eq.0) then
        write(6,*) 'done :: read .re2 file'
        write(6,*) ''
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,e,ifbswap)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'
      logical ifbswap

      integer e,eg,buf(0:30)

      nwds = 1 + ndim*(2**ndim) ! group + 2x4 for 2d, 3x8 for 3d

      if (ifbswap) call byte_reverse(buf,nwds)


      igroup(e) = buf(0)

      if (if3d) then
         call copy4r(xc(1,e),buf( 1),8)
         call copy4r(yc(1,e),buf( 9),8)
         call copy4r(zc(1,e),buf(17),8)
      else
         call copy4r(xc(1,e),buf( 1),4)
         call copy4r(yc(1,e),buf( 5),4)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_curve(buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      integer e,eg,f,buf(30)

      eg = buf(1)
      e  = gllel(eg)
      f  = buf(2)

      call copy4r( curve(1,f,e),buf(3),5)
      call chcopy(ccurve(f,e)  ,buf(8),1)

c     write(6,1) eg,e,f,(curve(k,f,e),k=1,5),ccurve(f,e)
c   1 format(2i7,i3,5f10.3,1x,a1,'ccurve')

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbl,bl,buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      character*3 cbl(6,lelt)
      real         bl(5,6,lelt)

      integer e,eg,f,buf(30)

      eg = buf(1)
      e  = gllel(eg)
      f  = buf(2)

      call copy4r ( bl(1,f,e),buf(3),5)
      call chcopy (cbl(  f,e),buf(8),3)

c     write(6,1) eg,e,f,cbl(f,e),' CBC',nid
c  1  format(2i8,i4,2x,a3,a4,i8)

      return
      end
c-----------------------------------------------------------------------
      subroutine bin_rd1_mesh(ifbswap)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'
      logical ifbswap

      integer e,eg,buf(30)

      nwds = 1 + ndim*(2**ndim) ! group + 2x4 for 2d, 3x8 for 3d
      len  = 4*nwds             ! 4 bytes / wd

      if (nwds.gt.30.or.isize.gt.4) then
         write(6,*) nid,' Error in bin_rd1_mesh: buf size',nwds,isize
         call exitt
      endif

      call gsync()

      do eg=1,nelgt             ! sync NOT needed here

         mid = gllnid(eg)
         e   = gllel (eg)

         if (nid.eq.0.and.mod(eg,1000).eq.0) write(6,*) eg,' mesh read'
         if (mid.ne.nid.and.nid.eq.0) then              ! read & send

            call byte_read  (buf,nwds)
            call csend      (eg,buf,len,mid,0)

         elseif (mid.eq.nid.and.nid.ne.0) then          ! recv & process

            call crecv      (eg,buf,len)
            call buf_to_xyz (buf,e,ifbswap)

         elseif (mid.eq.nid.and.nid.eq.0) then          ! read & process

            call byte_read  (buf,nwds)
            call buf_to_xyz (buf,e,ifbswap)

         endif

      enddo

      call gsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine bin_rd1_curve (ifbswap) ! v. 1 of curve side reader

      include 'SIZE'
      include 'TOTAL'
      logical ifbswap

      integer e,eg,buf(30)

      nwds = 2 + 1 + 5   ! eg + iside + ccurve + curve(6,:,:) !only 5 in rea
      len  = 4*nwds      ! 4 bytes / wd

      if (nwds.gt.30.or.isize.gt.4) then
         write(6,*)nid,' Error in bin_rd1_curve: buf size',nwds,isize
         call exitt
      endif

      call gsync()

      if (nid.eq.0) then  ! read & send/process

         call byte_read(ncurve,1)
         if (ifbswap) call byte_reverse(ncurve,1)

         do k=1,ncurve
            call byte_read(buf,nwds)
            if (ifbswap) call byte_reverse(buf,nwds-1) ! last is char
            eg  = buf(1)
            mid = gllnid(eg)
            if (mid.eq.0) then
               call buf_to_curve(buf)
            else
               call csend(mid,buf,len,mid,0)
            endif
         enddo
         call buf_close_out  ! notify all procs: no more data

      else               ! wait for data from node 0

         ncurve_mx = 8*nelt
         do k=1,ncurve_mx

            call crecv(nid,buf,len)

            if (buf(1).eq.0) then
               goto 99
            else
               call buf_to_curve(buf)
            endif
            
         enddo
   99    call buf_close_out

      endif

      call gsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine bin_rd1_bc (cbl,bl,ifbswap) ! v. 1 of bc reader

      include 'SIZE'
      include 'TOTAL'
      logical ifbswap

      character*3 cbl(6,lelt)
      real         bl(5,6,lelt)

      integer e,eg,buf(30)

      nwds = 2 + 1 + 5   ! eg + iside + cbc + bc(5,:,:)
      len  = 4*nwds      ! 4 bytes / wd

      if (nwds.gt.30.or.isize.gt.4) then
         write(6,*) nid,' Error in bin_rd1_bc: buf size',nwds,isize
         call exitt
      endif

      do e=1,nelt   ! fill up cbc w/ default
      do k=1,6
         cbl(k,e) = 'E  '
      enddo
      enddo

      call gsync()

      if (nid.eq.0) then  ! read & send/process

         call byte_read(nbc_max,1)
         if (ifbswap) call byte_reverse(nbc_max,1) ! last is char
         do k=1,nbc_max

c           write(6,*) k,' dobc1 ',nbc_max
            call byte_read(buf,nwds)
            if (ifbswap) call byte_reverse(buf,nwds-1) ! last is char

            eg  = buf(1)
            mid = gllnid(eg)
c           write(6,*) k,' dobc3 ',eg,mid

            if (mid.eq.0) then
               call buf_to_bc(cbl,bl,buf)
            else
c              write(6,*) mid,' sendbc1 ',eg
               call csend(mid,buf,len,mid,0)
c              write(6,*) mid,' sendbc2 ',eg
            endif

c           write(6,*) k,' dobc2 ',nbc_max,eg
         enddo
c        write(6,*) mid,' bclose ',eg,nbc_max
         call buf_close_outv ! notify all procs: no more data

      else               ! wait for data from node 0

         nbc_max = 2*ndim*nelt
         do k=1,nbc_max

c           write(6,*) nid,' recvbc1',k
            call crecv(nid,buf,len)
c           write(6,*) nid,' recvbc2',k,buf(1)

            if (buf(1).eq.0) then
               goto 99
            else
               call buf_to_bc(cbl,bl,buf)
            endif
            
         enddo
   99    call buf_close_outv

      endif

      call gsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_close_outv  ! this is the stupid O(P) formulation

      include 'SIZE'
      include 'PARALLEL'
      integer*4 zero

      len  = 4
      zero = 0
c     write(6,*) nid,' bufclose'
      if (nid.eq.0) then
         do mid=1,np-1
            call csend(mid,zero,len,mid,0)
c           write(6,*) mid,' sendclose'
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_close_out  ! this is the stupid O(P) formulation

      include 'SIZE'
      include 'PARALLEL'
      integer*4 zero

      len  = 4
      zero = 0
      if (nid.eq.0) then
         do mid=1,np-1
            call csend(mid,zero,len,mid,0)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine open_bin_file(ifbswap) ! open file & chk for byteswap

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap,if_byte_swap_test

      integer fnami (33)
      character*132 fname
      equivalence (fname,fnami)

      character*80 hdr
      character*5 version
      real*4      test

      if (nid.eq.0) then
         call izero(fnami,33)
         m = indx132(re2fle,' ',1)-1
         call chcopy(fname,re2fle,m)
   
         call byte_open(fname)
         call byte_read(hdr,20)

         read (hdr,1) version,nelgt,ndum,nelgv
    1    format(a5,i9,i3,i9)

         call byte_read(test,1)
         ifbswap = if_byte_swap_test(test)

      endif


      call lbcast(ifbswap)     ! broadcast byte swap from node 0

      nelgv = iglmax(nelgv,1)
      nelgt = iglmax(nelgt,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine close_bin_file

      include 'SIZE'
      include 'TOTAL'

      if (nid.eq.0) call byte_close()

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_xyz
      include 'SIZE'
      include 'TOTAL'
      integer e,f,eg

      do e=1,nelt
         eg = lglel(e,node)
         write(6,1) nid,eg,e,(cbc(f,e,1),f=1,6)
      enddo
    1 format(3i8,6(1x,a3),'  cbc')

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_nel
      include 'SIZE'
      include 'TOTAL'

      neltmx=np*lelt
      nelvmx=np*lelv

      neltmx=min(neltmx,lelg)
      nelvmx=min(nelvmx,lelg)

      if (nelgt.gt.neltmx.or.nelgv.gt.nelvmx) then
        if (nid.eq.0) write(6,*)'help:',lp,np,nelvmx,nelgv,neltmx,nelgt
        if (nid.eq.0) write(6,*)'help:',lelt,lelv,lelgv
        if (np.eq.1) write(6,11) lelt,lelg,nelgv,nelgt
        if (np.gt.1) write(6,12) lelt,lelg,nelgv,nelgt
   11       format(//,2X,'Error: This NEKTON Solver has been compiled'
     $             ,/,2X,'       for',2i9,' elements.'
     $             ,/,2X,'       The data file has dimensions',2i9,'.')
   12       format(//,2X,'Error: This NEKTON Solver has been compiled'
     $             ,/,2X,'       for',2i9,' elements per processor.'
     $             ,/,2X,'       The data file has dimensions',2i9,'.'
     $             ,/,6X,'   Rerun with more processors or recompile.')
         call exitt
      endif

      if (nelt.gt.lelt) then
        write(6,'(A,3I9)') 'Error: nelt>lelt!', nid, nelt, lelt
        call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
