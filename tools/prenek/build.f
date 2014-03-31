c-----------------------------------------------------------------------
C     
C     NEKTON 2.6  2/8/90
C     
C     Copyright (C) 1990, by the 
C     
C     Massachusetts Institute of Technology  and Nektonics, Inc.
C     
C     All Rights Reserved
C     
C     This program is a licenced product of MIT and Nektonics, Inc.,  and it is 
C     not to be disclosed to others, copied, distributed, or displayed 
C     without prior authorization.
C     
C     
c-----------------------------------------------------------------------
      subroutine build

      include 'basics.inc'

      if (ifconj_merge) then

         call build0             ! Read in fluid mesh
         nelv = nel

         call imp_mesh(.false.)  ! Get thermal mesh
         nelt  = nel
         ncond = nelt-nelv

         call set_igroup(nelv,nelt)

         call build2             ! Set bcs

         nelt  = nel
         ncond = nelt-nelv

      else          ! Std. menu-driven mesh construction

         call build0
         call build1

         nelt = nel
         nelv = nel
         call build2

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine build0

c     Menu-based module that prompts the user to input corners.
     
      include 'basics.inc'
      dimension icrvs(4)
      character key,string*6,leto,char1
      logical iftmp
      common /splitt/ enew(nelm),ind(nelm)

      if (ifnoseg) call drcovr(13)
     
      pi=4.0*atan(1.0)
      nx=nxm
      call legend(zpts,wght,nx)

      IF (NX.GT.NXM) THEN
         WRITE(S,104) NX,NXM
 104     FORMAT(' Warning, current N exceeds NXM, resetting from',I3
     $        ,' to',I3,'.$')
         CALL PRS(S)
         NX=NXM
      ENDIF
      NY=NX
      NZ=1
      IF(IF3D)NZ=NX
      IF(NKTONV.EQ.1) THEN
         DO 10 I=1,NX
            XXPTS(I) = -COS(PI*FLOAT(I-1)/FLOAT(NX-1))
            YYPTS(I) = -COS(PI*FLOAT(I-1)/FLOAT(NX-1))
 10      CONTINUE
      ELSE
C     Need legendre points for mesh drawn on screen
         CALL LEGEND(XXPTS,WGHT,NX)
         CALL LEGEND(YYPTS,WGHT,NY)
      ENDIF
      NSIDES=4
      IFCEN=.FALSE.
      IF(IF3D)NSIDES=6

      IF(CHOICE.EQ.'BUILD FROM FILE'.OR.
     $     CHOICE.EQ.'IMPORT UNIVERSAL FILE')THEN
         if (choice.eq.'BUILD FROM FILE') then
           call readat
           close(unit=9)
           if (.not.if3d) call chk_right_hand(nel)
         endif
C     Set grid stuff based on old data
         GRIDDX=XFAC*GRID
         GRIDDY=YFAC*GRID
C     Set clipping window based on scale factors
         WT = YPHY(.995)
         WB = YPHY(.005)
         WL = XPHY(.005)
         WR = XPHY(.995)
C     Save scale factors for UnZoom function
         XFACO=XFAC
         YFACO=YFAC
         XZEROO=XZERO
         YZEROO=YZERO
C     Now set stuff based on data from readat
         do 20 ie=1,nel
            height(numapt(ie))=z(5,ie)-z(1,ie)
            xcen(ie)=(x(1,ie)+x(2,ie)+x(3,ie)+x(4,ie))/4.0
            ycen(ie)=(y(1,ie)+y(2,ie)+y(3,ie)+y(4,ie))/4.0
 20      continue
C     Display Mesh (or at least 1st floor)
         CALL DRGRID
         NLEVEL=0
C     Sort out the elements so that they will be drawn correctly
         CALL SORTEL
         DO 30 IEL=1,NEL
            CALL DRAWIS(ISRT(IEL))
            IF(NUMAPT(IEL).EQ.1)CALL DRAWEL(IEL)
            IF(NUMAPT(IEL).GT.NLEVEL)NLEVEL=NUMAPT(IEL)
 30      CONTINUE
C     
         IF(NDIM.EQ.2)NLEVEL=1
C     Display Elevator 1st floor hilighted
         IF(NLEVEL.GT.1)THEN
            DO 40 I=NLEVEL,2,-1
               CALL DRELEV(I-1,I,'     ')
 40         CONTINUE
         ENDIF
         CALL SCROLL(5)
         CALL GNABLE
         ITEM(1)='ACCEPT MESH'
         ITEM(2)='REVIEW/MODIFY'

c        goto 331  ! quick hack       <---------------------------!
c        ITEM(3)='Edit Mesh'                                      !
c        NCHOIC =  2                                              !
c        CALL MENU(XMOUSE,YMOUSE,BUTTON,'ACCEPT/REVIEW')          !
c        IF(IFNOSEG)CALL DRCOVR(13)                               !
c        IF(CHOICE.EQ.'ACCEPT MESH')THEN                          !
c           goto 330                                              !
c        ELSE IF(CHOICE.EQ.'Edit Mesh')THEN                       !
c           call mesh_edit                                        !
c        ELSE IF(CHOICE.EQ.'REVIEW/MODIFY')THEN                   !
c           prepare to modify floor by floor                      !
c        ENDIF                                                    !
c 331    continue   ! quick hack       <--------------------------!


      ELSE
C     Interactive Input
         CALL SCROLL(5)
         CALL GNABLE
         CALL SETSCL
      ENDIF

      return
      end
c-----------------------------------------------------------------------
      subroutine build1

c     Menu-based module that prompts the user to input corners.
     
      include 'basics.inc'
      dimension icrvs(4)
      character key,string*6,leto,char1
      logical iftmp
      common /splitt/ enew(nelm),ind(nelm)


      if (ifmerge) return

C     Just in case it didn't get set in setscl
      IFGRID=.TRUE.
      IF(NLEVEL.EQ.0)NLEVEL=1
      CALL PRS('                    *** BUILD MENU ***$')
C     
C     
      IFCEN=.FALSE.
      IF(.NOT.IFREAD)CALL DRGRID
      IF(.NOT.IFREAD)NEL=0
C     Start at First level
      ILEVEL=1
      IF(IF3D) THEN
         IF(IFREAD)THEN
C     Put in bogus height of first level
         ELSE
C     For interactive sesion, Demand Height of First Level
            CALL PRSI('Please Enter Height of Level$',ILEVEL)
 50         call rer(HEIGHT(ILEVEL))
            IF(HEIGHT(ILEVEL).LE.0.0)THEN
               CALL PRS('HEIGHT must be a positive number!$')
               goto 50
            ENDIF
         ENDIF
C     2nd arg: 0 for no fade; negative for bottom only (for floor modify)
         CALL DRELEV(ILEVEL,0,'     ')
      ENDIF
C     Draw FIRST Level
      IF(NEL.GT.0) THEN
         DO 60 IEL=1,NEL
C     CALL DRAWIS(ISRT(IEL))
C     IF(NUMAPT(IEL).EQ.ILEVEL)CALL DRAWEL(IEL)
            IF (NUMAPT(IEL).EQ.ILEVEL)
     $           MAXLET=MAX(MAXLET,ICHAR(LETAPT(IEL)))
 60      CONTINUE
      ENDIF
      ILETAP=MAXLET+96-32
C     ! ??!!
C     
 1000 CONTINUE

      call gencen
      call mkside


C***  BIG DECISION POINT  ****

      call set_main_menu
      call menu(xmouse,ymouse,button,'NOCOVER') ! Prompt for user input.
     
      IF (CHOICE.EQ.'ADD    ELEMENT') THEN
         if(nel.ge.nelm-3)then
            call prs('Too many elements. Increase nelm in basics.inc.$')
            goto 1000
         endif

         ne1=nel+1
         numapt(ne1)=ilevel
         if (iletap.le.122) letapt(ne1)=char(iletap)
         if (iletap.gt.122) letapt(ne1)=' '

         call prs('Enter element corners. Use mouse,$')
         call prs('or click menu area to type x-y coordinate$')

C     Turn on Keypad

         call newel(xmouse,ymouse,button,ierr)
         nel = nel+1

         if (ierr.eq.1) then
            call drawel(-nel)
            nel=nel-1
         endif

      ELSE IF(CHOICE.EQ.'CEILING') THEN
C        Make sure everything is drawn on ceiling.  And that
C        MODEL and CURVE know about it, too
         DO 80 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 80      CONTINUE
         IFCEIL=.TRUE.
         DO 90 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
 90      CONTINUE
         CALL DRELEV(-(ILEVEL+1),ILEVEL,'     ')
      ELSE IF(CHOICE.EQ.'OFF CEILING') THEN
         DO 100 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 100     CONTINUE
         IFCEIL=.FALSE.
         DO 110 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(I)
 110     CONTINUE
         CALL DRELEV(ILEVEL,-(ILEVEL+1),'     ')
      ELSE IF(CHOICE.EQ.'IMPORT MESH')THEN
         call imp_mesh(.true.)
         if (.not.if3d) call chk_right_hand(nel)
         call redraw_mesh_small
      ELSE IF(CHOICE.EQ.'REFLECT MESH')THEN
         CALL REFLECT_MESH
      ELSE IF(CHOICE.EQ.'IMPORT VTK MESH')THEN
         call imp_mesh_vtk
         if (.not.if3d) call chk_right_hand(nel)
      ELSE IF(CHOICE.EQ.'IMPORT vtx MESH')THEN
         call imp_mesh_vtx
         if (.not.if3d) call chk_right_hand(nel)
      ELSE IF(CHOICE.EQ.'CURVE SIDES')THEN
         CALL CURVES
      ELSE IF(CHOICE.EQ.'MODIFY ELEMENT')THEN
C     Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
         ELSE
            CALL MODEL(NEL)
         ENDIF
      ELSE IF(CHOICE.EQ.'GLOBAL REFINE')THEN
C     Normal Modify
         IF(NEL.EQ.0)THEN
            CALL PRS('ERROR: No elements to modify$')
         ELSE
C     Only floor of elevator hilighted during modify
            CALL DRELEV(-ILEVEL,ILEVEL,'     ')
            CALL GLOMOD
            CALL DRELEV(ILEVEL,0,'     ')
         ENDIF
      elseif (choice.eq.'DELETE ELEMENT') then
         call delete
      elseif (choice.eq.'REDRAW ISOMETRIC') then ! redraw isometric elements.
         call sortel
         do 160 i=1,nel
            call drawis(isrt(i))
 160     continue
      ELSE IF(CHOICE.EQ.'REDRAW MESH')THEN
         call redraw_mesh
      ELSE IF(CHOICE.EQ.'ZOOM')THEN
         call setzoom
         call redraw_mesh
      ELSE IF(CHOICE.EQ.'SET GRID')THEN
         call setgrd
      ELSE IF(CHOICE.EQ.'Edit Mesh')THEN
         call mesh_edit
      ELSE IF(CHOICE.EQ.'DEFINE OBJECT')THEN
         CALL SETOBJ
      ELSE IF(CHOICE.EQ.'END    ELEMENTS')THEN
C      WHAT ELSE TO DO WHEN 2-D PROBLEM?
         IF(NEL.EQ.0) THEN
            CALL PRS('ERROR: Can''t "END ELEMENTS" when no elements$')
            CALL PRS('are input. ^C if you want to give up$')
         ELSE
            NLEVEL=0
            DO 180 I=1,NEL
               NLEVEL=MAX(NLEVEL,NUMAPT(I))
 180        CONTINUE
            goto 320
         ENDIF
      ELSE IF(CHOICE.EQ.'UP   LEVEL') THEN
C     Make Checks
         NTHISL=0
         NNEXTL=0
C     !!?? maybe 96-32??  also, sometimes up gives wrong ele
         MAXLET=0
         DO 190 I=1,NEL
            IF(NUMAPT(I).EQ.ILEVEL  )NTHISL=NTHISL+1
            IF(NUMAPT(I).EQ.ILEVEL+1)NNEXTL=NNEXTL+1
            IF(NUMAPT(I).EQ.ILEVEL  )MAXLET=
     $           MAX(MAXLET,ICHAR(LETAPT(I)))
 190     CONTINUE
         IF(NTHISL.EQ.0)THEN
            CALL PRSI('*** ERROR *** No elements on level$',ILEVEL)
            CALL PRS('Can''t make level above empty one$')
         ENDIF
         IF(NNEXTL.EQ.0)THEN
            CALL PRS('** New Level **  Enter Height:$')
            CALL PRS('Negative to Abort start of new level$')
            call rer(HEIGHT(ILEVEL+1))
C     !!?? HOW TO MODIFY HEIGHT ONCE IT'S IN?
            IF(HEIGHT(ILEVEL+1) .LE. 0.0) THEN
C     Abort level change
               CALL PRSI('Aborting level change. Still on $',ilevel)
               HEIGHT(ILEVEL+1)=0
            ELSE IF(HEIGHT(ILEVEL+1) .GT. 0.0) THEN
C     Really go up one
               CALL PRSIS('Using Ceiling of Level$',ILEVEL,
     $              '  Mesh as default$')
               DO 200 I=1,NEL !     Erase old mesh
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 200           CONTINUE
               ILEVEL=ILEVEL+1
               CALL DRELEV(ILEVEL,ILEVEL-1,'     ')
               ILETAP=MAXLET
               NLEVEL=NLEVEL+1
            ENDIF
            IELNOW=NEL
            DO 250 I=1,IELNOW
               IF(NUMAPT(I).EQ.ILEVEL-1)THEN
                  NEL=NEL+1
                  CALL COPYEL(I,NEL)
                  NUMAPT(NEL)=ILEVEL
C     Correct Stuff affected if walls are not vertical
                  XCEN(NEL)=0.0
                  YCEN(NEL)=0.0
C     Remove periodic b.c.'s
                  DO 210 IF=1,NFLDS
                     DO 210 ISIDE=1,NSIDES
                        IF(  CBC( ISIDE, NEL,IF).EQ.'P  ')THEN
                           CBC(   ISIDE, NEL,IF)='   '
                           BC (1, ISIDE, NEL,IF)= 0
                           BC (2, ISIDE, NEL,IF)= 0
                           BC (3, ISIDE, NEL,IF)= 0
                           ibc(   ISIDE, NEL,IF)= 0
                        ENDIF
 210                 CONTINUE
                     DO 230 IC=1,4
                        x(ic,nel)=x(ic+4,i)
                        y(ic,nel)=y(ic+4,i)
                        z(ic,nel)=z(ic+4,i)
                        SIDES(NEL,IC,3)=SIDES(NEL,IC,3)+
     $                       (HEIGHT(ILEVEL-1)+HEIGHT(ILEVEL))/2.
                        xcen(nel)=xcen(nel)+x(ic,nel)/4.
                        ycen(nel)=ycen(nel)+y(ic,nel)/4.
C     Curved side stuff
                        CCURVE(IC,NEL)=CCURVE(IC+4,I)
                        DO 220 II=1,6
                           CURVE(II,IC,NEL)=CURVE(II,IC+4,I)
 220                    CONTINUE
 230                 CONTINUE
C     RAISE Z OF CEILING OF NEW ELEMENTS
                     DO 240 IC=5,8
                        z(ic,nel)=z(ic,nel)+height(ilevel)
 240                 CONTINUE
                     SIDES(NEL,5,3)=SIDES(NEL,5,3)+HEIGHT(ILEVEL-1)
                     SIDES(NEL,6,3)=SIDES(NEL,6,3)+HEIGHT(ILEVEL)
                     CALL DRAWEL(NEL)
                  ENDIF
 250           CONTINUE
C     Sort out the elements so that they will be drawn correctly
               CALL SORTEL
               DO 260 I=1,NEL
                  IF(NUMAPT(ISRT(I)).EQ.ILEVEL)CALL DRAWIS(ISRT(I))
 260           CONTINUE
            ELSE
C     Already have elements on next level. Erase old mesh& draw new
               CALL DRELEV(ILEVEL+1,ILEVEL,'     ')
               DO 270 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 270           CONTINUE
               ILEVEL=ILEVEL+1
               DO 280 I=1,NEL
                  IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
 280           CONTINUE
            ENDIF
         ELSE IF(CHOICE.EQ.'DOWN LEVEL')THEN
C     Make Checks
            IF(ILEVEL.EQ.1)THEN
               CALL PRS('You already are on Level 1.  You cannot$')
               CALL PRS('go down any further.  Choose something else.$')
            ENDIF
            NTHISL=0
            NDOWNL=0
            MAXLET=0
            DO 290 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL  )NTHISL=NTHISL+1
               IF(NUMAPT(I).EQ.ILEVEL-1)NDOWNL=NDOWNL+1
               IF(NUMAPT(I).EQ.ILEVEL-1)MAXLET=
     $              MAX(MAXLET,ICHAR(LETAPT(I)))
 290        CONTINUE
C     Go down one level.  Erase old mesh& draw new
            CALL DRELEV(ILEVEL-1,ILEVEL,'     ')
            DO 300 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL(-I)
 300        CONTINUE
            ILEVEL=ILEVEL-1
            DO 310 I=1,NEL
               IF(NUMAPT(I).EQ.ILEVEL) CALL DRAWEL( I)
 310        CONTINUE
         ELSE
            CALL PRS('CHOICE:'//CHOICE//'NOT IN MENU$')
         ENDIF
         goto 1000
 320     CONTINUE
 330     CONTINUE

      return   ! End of menu-driven query
      end
c-----------------------------------------------------------------------
      subroutine build2
     
c     Establish connectivity, set global node numbers, and
c         construct elemental mesh. 

c     Call routines that set boundary conditions.
     
      include 'basics.inc'
      dimension icrvs(4)
      character key,string*6,leto,char1
      logical iftmp
      common /splitt/ enew(nelm),ind(nelm)
      real znew(nelm)
      equivalence (znew,enew)
      integer e


C     Now, cover menu
         CALL SGVIS(13,0)
         CALL SGVIS(13,1)
         nelv= nel
C     
C     Finished with gridlatch
         IFGRID=.FALSE.
 330     CONTINUE
C     Variable properties put in here
         CALL VPROPS
C     JUMP IN HERE IF YOU DON'T WANT TO MODIFY MESH YOU JUST READ IN
C     Find Sides' Midpoints
         CALL MKSIDE
C     Find Sides' Overlaps (used for finding boundary sides)
C!!   Check here for Irrational B.c.'s: unrequited connectivities;
C     Warn if physical boundaries internally.  This check is not exhaustive!!??

      ifld0=1
      ifld1=nflds

      if (ifconj_merge) then  ! Update THERMAL BCS  ONLY
         ifld0=2
         ifld1=2
      endif

      nsides = 2*ndim
      do 335 ifld=ifld0,ifld1 !  ZERO OUT OLD ELEMENTAL CONNECTIVITY
      do 335 iel=1,nel
      do 335 iside=1,nsides
         IF (CBC( ISIDE, IEL,Ifld).EQ.'E  ') then
c        IF (CBC( ISIDE, IEL,Ifld).EQ.'E  ' .or.
c    $       CBC( ISIDE, IEL,Ifld).EQ.'P  ' .or.
c    $       CBC( ISIDE, IEL,Ifld).EQ.'SP ' .or.
c    $       CBC( ISIDE, IEL,Ifld).EQ.'J  ') then
            CBC(   ISIDE, IEL,Ifld)='   '
            BC (1, ISIDE, IEL,Ifld)= 0
            BC (2, ISIDE, IEL,Ifld)= 0
            BC (3, ISIDE, IEL,Ifld)= 0
            ibc(   ISIDE, IEL,Ifld)= 0
         endif
 335  continue
C     
C     Make internal side comparisons
C     
      icount=1
      do jcount=1,8
         if (icount*100.gt.nel) goto 345
         icount = icount*10
      enddo
  345 continue


c     Find 'E-E' boundary pairs in n log n time

      call cell_cell_connectivity

c     All 'E-E' boundaries are set.  Check for nonconforming (SP-J)

      call check_spj


c     All internal boundaries are set.  Query remaining bcs

      call gingrd(0.0)
      call scroll(7)

      call bound
      call bound
     
      nsides = 2*ndim
      ifmvbd=.false.
      do 381 Ifld=1,NFLDS
      do 381 IEL=1,NEL
      do 381 ISIDE=1,NSIDES
         CHAR1 = CBC(ISIDE,IEL,Ifld)
         if (char1.eq.'m'.or.char1.eq.'M') then
            ifmvbd=.true.
            goto 382
         endif
 381  continue
 382  continue
      if (ifmvbd) ifxyo = .true.
     
      if (.not.ifmerge) call mesgen
     
      return
      end
c-----------------------------------------------------------------------
      subroutine readat
      include 'basics.inc'
      logical iffold,ifhold
      character chtemp*3
      character*80 string
      real*8 bc8(5)

      NLINF=0
      NLINP=0
      NLINR=0
     
C     Read in all the build,bound,curve,history,output,initcond data
C     Read Dummy Parameters
      READ(9,*,ERR=33)
      READ(9,*,ERR=33)VNEKOLD
      nktold=VNEKOLD
      READ(9,*,ERR=33)
      READ(9,*,ERR=33)NDUM
      do 188 IDUM=1,NDUM
         read(9,*,err=33)xxx
         if(idum.eq.23)npsold=xxx
 188  CONTINUE
C     Read Passive scalar data
      read(9,*,err=33)NSKIP
      IF(NSKIP.GT.0) THEN
         do 177 I=1,NSKIP
            read(9,*,err=33)
 177     CONTINUE
      endif
C     
C     READ DUMMY LOGICALS
      READ(9,*,ERR=33)NLOGIC
      do 189 I=1,NLOGIC
         if(i.eq.1)read(9,*,err=33)iffold
         if(i.eq.2)read(9,*,err=33)ifhold
         if(I.ne.1.and.i.ne.2)READ(9,*,ERR=33)
 189  CONTINUE
      NoLDS=1
      IF(ifhold)NoLDS=2+NPSold
C     
C     Read Elemental Mesh data
      READ(9,*,ERR=33,end=34)XFAC,YFAC,XZERO,YZERO
      READ(9,*,ERR=33)
      READ(9,*,ERR=33,END=33)NEL,NDIM

      if (nel.ge.nelm) then
         call prsii('NELM too small in basics.inc.$',nel,nelm)
         call prs  ('Recompile prenek.  ABORT$')
      endif
     
      iffmtin = .true.
      IF (NEL.lt.0) iffmtin = .false.
      IF (NEL.lt.0) NEL = -NEL
c     
      IF(IF3D.AND.NDIM.EQ.2)THEN
         CALL PRS
     $        ('WARNING: PARAMETER SET TO 2D BUT DATA IS FOR 3D$')
         CALL PRS('PLEASE ENTER Z-HEIGHT: $')
         CALL RER( ZHGHT)
      ENDIF
      IF(.NOT.IF3D.AND.NDIM.EQ.3) THEN
         CALL PRS(
     $        'ERROR: PARAMETER SET TO 3-D BUT DATA IS FOR 2-D$')
      ENDIF
      NCOND=0
      nel50=nel/50
      if (nel.lt.100) nel50=500
      do 98 IEL=1,NEL

c        READ(9,'(20X,I4,4X,I3,A1,11x,i5)',ERR=33,END=33)
c    $        IDUM,NUMAPT(IEL),LETAPT(IEL),IGROUP(IEL)

         call blank(string,80)
         read(9,80,err=33,end=33) string
   80    format(a80)
c        write(6,*) iel,' ',string

         call parse_e_hdr
     $       (idum,numapt(iel),letapt(iel),igroup(iel),string)


         IF(IGROUP(IEL).GT.0) NCOND=NCOND+1
         IF (NEL.LT.100.or.mod(iel,nel50).eq.0) THEN
            WRITE(S,'(A20,I8,A4,I5,A1,A2)',ERR=33)
     $           ' Reading Element',
     $           idum,' [',NUMAPT(IEL),LETAPT(IEL),']$'
            CALL PRS(S)
         ENDIF
         IF(NDIM.EQ.2)THEN
            read(9,*,err=33,end=33)(x(ic,iel),ic=1,4)
            read(9,*,err=33,end=33)(y(ic,iel),ic=1,4)
            IF (IF3D) THEN
               do 197 IC=1,4
                  I4 = IC+4
                  z(ic,iel) = 0.0
                  x(i4,iel) = x(ic,iel)
                  y(i4,iel) = y(ic,iel)
                  z(i4,iel) = zhght
 197           CONTINUE
            ENDIF
         ELSE IF(NDIM.EQ.3)THEN
            read(9,*,err=33,end=33)(x(ic,iel),ic=1,4)
            read(9,*,err=33,end=33)(y(ic,iel),ic=1,4)
            read(9,*,err=33,end=33)(z(ic,iel),ic=1,4)
            read(9,*,err=33,end=33)(x(ic,iel),ic=5,8)
            read(9,*,err=33,end=33)(y(ic,iel),ic=5,8)
            read(9,*,err=33,end=33)(z(ic,iel),ic=5,8)
C     For Flat floors only.  On write will flatten floors.
         ENDIF
 98   CONTINUE
      nelv=nel-ncond
C     Read curved side data
      READ(9,*,ERR=57,END=57)
      READ(9,*,ERR=57,END=57)NCURVE
      IF(NCURVE.GT.0)THEN
         do 19 I=1,NCURVE

            if (nel.lt.1000) then
               READ(9,'(I3,I3,5G14.6,1X,A1)',ERR=57,END=57)
     $              IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            elseif (nel.lt.1 000 000) then
               READ(9,'(I2,I6,5G14.6,1X,A1)',ERR=57,END=57)
     $              IEDGE,IEL,R1,R2,R3,R4,R5,ANS
            else
               read(9,'(i2,i12,5g14.6,1x,a1)',err=158,end=158)
     $              iedge,iel,r1,r2,r3,r4,r5,ans
            endif
            goto 160
  158       continue
              read(9,'(i2,i12,5g14.6,1x,a1)',err=57,end=57) ! old format
     $              iedge,iel,r1,r2,r3,r4,r5,ans            ! pff, feb 2013
  160       continue

            CALL DDUMMY(IEDGE,IEL)
            CURVE(1,IEDGE,IEL)=R1
            CURVE(2,IEDGE,IEL)=R2
            CURVE(3,IEDGE,IEL)=R3
            CURVE(4,IEDGE,IEL)=R4
            CURVE(5,IEDGE,IEL)=R5
            CCURVE( IEDGE,IEL)=ANS
 19      CONTINUE
      ENDIF

C     Read Boundary Conditions (and connectivity data)

      if (iffmtin) READ(9,*,ERR=44,END=44)
      do 90 IFLD=1,NoLDS
C     Fluid and/or thermal
C     !!?? nelv different from nel??
         if (iffmtin) READ(9,*,ERR=44,END=44)
         IF( (IFLD.EQ.1.AND.IFFOLD) .OR. IFLD.GT.1)THEN
C     FIX UP FOR WHICH OF FIELDS TO BE USED
            do 88 IEL=1,NEL
               do 88 ISIDE=1,NSIDES
C     !Fix to a4,i2 when you make cbc character*4
                  IF(VNEKOLD .LE. 2.5) NBCREA = 3
                  IF(VNEKOLD .GE. 2.6) NBCREA = 5
                  IF (NEL.LT.1000.and.iffmtin) THEN
                     READ(9,'(1X,A3,1x,I2,I3,5G14.6)',ERR=44,END=44)
     $                  CBC(ISIDE,IEL,IFLD),ID,ID,(bc8(ii),ii=1,nbcrea)
                  ELSEIF (iffmtin.and.nel.lt. 100 000) then
                     READ(9,'(1X,A3,I5,I1,5G14.6)',ERR=441,END=441)
     $                  CBC(ISIDE,IEL,IFLD),ID,ID,(bc8(ii),ii=1,nbcrea)
                  ELSEIF (iffmtin.and.nel.lt. 1 000 000) then
                     READ(9,'(1X,A3,6x,5G14.6)',ERR=441,END=441)
     $                  CBC(ISIDE,IEL,IFLD),(bc8(ii),ii=1,nbcrea)
                  ELSEIF (iffmtin.and.nel.lt.2 000 000 000) then
                     read(9,'(1x,a3,i12,5g18.11)',err=2443,end=2443)
     $                  cbc(iside,iel,ifld),id,(bc8(ii),ii=1,nbcrea)
                  ELSE
                     READ(8,ERR=443,END=443) chtmp3,
     $                  CBC(ISIDE,IEL,IFLD),ID,ID,(bc8(ii),ii=1,nbcrea)
                  ENDIF
                  goto 2444
 2443             continue
                     read(9,'(1x,a3,9x,i1,5g14.6)',err=442,end=442) ! old format
     $                  cbc(iside,iel,ifld),id,(bc8(ii),ii=1,nbcrea)
 2444             continue
                  do ii=1,nbcrea
                     bc(ii,iside,iel,ifld)=bc8(ii)
                  enddo
                  ibc(iside,iel,ifld)=bc8(1)
                  CBC1=CBC(ISIDE,IEL,IFLD)
                  IF(CBC1.EQ.'M'.OR.CBC1.EQ.'m')THEN
                     IFMOVB=.TRUE.
                     IFGNGM=.TRUE.
                  ENDIF
                  ICBC1=ICHAR(CBC1)
                  IF(ICBC1.GE.97 .AND. ICBC1.LE.122)THEN
                     CBC3 = CBC(ISIDE,IEL,IFLD)
c
c                    pff  no nlines!!
                     NLINES=0
                     do k=1,5
                        BC(k,ISIDE,IEL,IFLD) = 0
                     enddo
                  ENDIF
 88            CONTINUE
               ENDIF
 90         CONTINUE
            IF(NFLDS.EQ.1.and.iffmtin)READ(9,*,err=45,end=45)

C           Read Initial Conditions
            IF(VNEKOLD.LT.2.5)THEN
               CALL PRS('***** OLD VERSION RESTART DATA ****$')
               CALL PRS(
     $              'Please Re-enter IC''s from the options menu$')
               READ(9,*,ERR=50,END=50)
               do 21 I=1,7
 21               READ(9,'(A80)',ERR=50,END=50)
               ELSE
C     Check for 2.6 style read file - separate restart/fortran sections.
                  READ(9,'(A70)',ERR=50,END=50) LINE
                  IF (INDX1(LINE,'RESTART',7).NE.0) THEN
                     REWIND(13)
                     WRITE(13,'(A70)')LINE
                     REWIND(13)
                     READ (13,*,ERR=50,END=50) NLINR
                     REWIND(13)
                     NSKIP=NLINR
                     do 201 I=1,NSKIP
                        READ(9,'(A80)')INITP(I)
                        LINE=INITP(I)
                        CALL CAPIT(LINE,70)
                        IF (INDX1(LINE,'PRESOLV',7).NE.0) THEN
                           NLINP=1
                           NLINR=NLINR-1
                        ENDIF
 201                 CONTINUE
                     READ(9,'(A70)',ERR=50,END=50) LINE
                  ENDIF
C     Read fortran initial condition data.
                  REWIND(13)
                  WRITE(13,'(A70)')LINE
                  REWIND(13)
                  READ(13,*,ERR=50,END=50)NLINF
                  REWIND(13)
                  do 22 I=1,NLINF
                     READ(9,'(A80)')INITC(I)
 22               CONTINUE
               ENDIF
C     Read drive force data
               READ(9,*,ERR=50,END=50)
               READ(9,*,ERR=50,END=50)nskip
               do 23 i=1,nskip
 23               READ(9,'(A80)',ERR=50,END=50)DRIVC(I)
C     READ CONDUCTION ELEMENT DATA
                  READ(9,*,ERR=60,END=60)
                  NFLDSC=0
                  do 117 IFLD=1,NFLDS
                     IF(IFTMSH(IFLD))NFLDSC=NFLDSC+1
 117              CONTINUE
                  READ(9,*,ERR=60,END=60)NSKIP
                  READ(9,*,ERR=60,END=60)NPACKS
                  IF(NPACKS.GT.0)THEN
                     do 218 IIG=1,NPACKS
                        READ(9,*)IGRP,IF,ITYPE
                        MATYPE(IGRP,IF)=ITYPE
                        do 218 IPROP=1,3
                           IF(ITYPE.EQ.1) READ(9,*      ) 
     $                          CPROP(IGRP,IF,IPROP)
                           IF(ITYPE.EQ.2) READ(9,'(A80)') 
     $                          VPROP(IGRP,IF,IPROP)
 218                    CONTINUE
                     ENDIF
C     
C     Read history data
                     READ(9,*,ERR=50,END=50)
                     READ(9,*,ERR=50,END=50)NHIS
C     HCODE(11) IS WHETHER IT IS HISTORY, STREAKLINE, PARTICLE, ETC.
                     IF(NHIS.GT.0)THEN
                        do 51 I=1,NHIS
                           READ(9,'(1X,11A1,1X,4I5)',ERR=50,END=50)
     $                          (HCODE(II,I),II=1,11),
     $                          (LOCHIS(II,I),II=1,4)
 51                     CONTINUE
                     ENDIF
C     Read output specs
                     READ(9,*,ERR=50,END=50)
                     READ(9,*,ERR=50,END=50)NOUTS
                     READ(9,*,ERR=50,END=50)IFXYO
                     READ(9,*,ERR=50,END=50)IFVO
                     READ(9,*,ERR=50,END=50)IFPO
                     READ(9,*,ERR=50,END=50)IFTO
                     READ(9,*,ERR=50,END=50)IFTGO
                     READ(9,*,ERR=50,END=50)IPSCO
                     IF(IPSCO .GT.0)THEN
                        do 1221 I=1,IPSCO
                           READ(9,'(1X,L1,1X,A5)',ERR=50,END=50) 
     $                          IFPSCO(I),PSNAME(I)
 1221                   CONTINUE
                     ENDIF
C     OBJECT SPECIFICATION DATA
                     READ(9,*,ERR=50,END=50)
                     READ(9,*,ERR=50,END=50)NSOBJS
                     IF(NSOBJS .GT. 0)THEN
                        do 62 IOBJ=1,NSOBJS
                           READ(9,'(4X,I4,35X,A20)',ERR=50,END=50)
     $                          NFACE(IOBJ),SOBJ(IOBJ)
                           do 63 IFACE=1,NFACE(IOBJ)
                              READ(9,*,ERR=50,END=50)
     $                             ILSURF(1,IOBJ,IFACE),
     $                             ILSURF(2,IOBJ,IFACE)
 63                        CONTINUE
 62                     CONTINUE
                     ENDIF
                     READ(9,*,ERR=50,END=50)NVOBJS
                     READ(9,*,ERR=50,END=50)NEOBJS
                     READ(9,*,ERR=50,END=50)NPOBJS
C     
                     return
 33                  CALL PRS
     $                    ('ERR: MESH DATA MISSING OR CORRUPTED.$')
                     CALL BEEP
                     CALL PRS('<CR> TO CONTINUE$')
                     CALL RES(LINE,0)
                     CALL DEVEX
                     CALL EXITT
 34                  CALL PRS('ERR: FILE CONTAINS ONLY PARAM DATA.$')
                     CALL PRS('ERR: MESH DATA MISSING$')
                     CALL BEEP
                     CALL PRS('<CR> TO CONTINUE$')
                     CALL RES(LINE,0)
                     CALL DEVEX
                     CALL EXITT
 35                  CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $                    '  CONTINUEING WITH MESH ONLY.$')
                     return
 44                  WRITE(S,'(1X,A39,I4,A6,I3)')
     $                    '0. Err reading B.C. data for elm',IEL,
     $                    ' side ',ISIDE
                     CALL PRS(S//'$')
c                    kzm = 1/(iel-iside)
                     return
 441                 WRITE(S,'(1X,A39,I4,A6,I3)')
     $                    'A. Err reading B.C. data for elm',IEL,
     $                    ' side ',ISIDE
                     CALL PRS(S//'$')
c                    kzm = 1/(iel-iside)
                     return
 442                 WRITE(S,'(1X,A39,I4,A6,I3)')
     $                    'B. Err reading B.C. data for elm',IEL,
     $                    ' side ',ISIDE
                     CALL PRS(S//'$')
c                    kzm = 1/(iel-iside)
                     return
 443                 WRITE(S,'(1X,A39,I4,A6,I3)')
     $                    'C. Err reading B.C. data for elm',IEL,
     $                    ' side ',ISIDE
                     CALL PRS(S//'$')
c                    kzm = 1/(iel-iside)
                     return

 45                  CALL PRS('FILE DOES NOT CONTAIN VALID B.C. DATA.'//
     $                    '  CONTINUEING WITH MESH ONLY.$')
                     return
 40                  CALL PRS(
     $                    'CURVED SIDES NOT READ CORRECTLY. PLS CHECK$')
                     CALL PRS('SIDES IN CURVED SIDE MENU$')
 50                  CALL PRS
     $                    ('OPTIONS MAY NOT HAVE BEEN READ CORRECTLY.$')
                     CALL PRS('PLS CHECK THEM IN OPTIONS MENU$')
                     return
 57                  CALL PRS
     $                    ('CURVED SIDES HAVE NOT BEEN READ CORRECTLY$')
                     CALL PRS('PLS CHECK THEM.$')
                     return
 60                  CALL PRS
     $                    ('CONDUCTION ELEMENTS MAY NOT HAVE BEEN$')
                     CALL PRS('READ CORRECTLY$')
                     CALL PRS('PLEASE CHECK THEM IN OPTIONS MENU$')
 70                  CALL PRS('NO OBJECTS$')
      return
      END
c-----------------------------------------------------------------------
      subroutine setscl
C     Sets scale factors
C     
      include 'basics.inc'
      LOGICAL IFMENU
      REAL XSC(4),YSC(4)
      GRID=PARAM(18)
 3    CONTINUE
      CALL CLEAR
      CALL REFRESH
      CALL DRGRID
c     CALL DRMENU('NOCOVER')
      CALL PRS(' SET SCALE FACTORS $')
c     NCHOIC =  3
      ITEM(1)='MOUSE'
      ITEM(2)='TYPE'
      ITEM(3)='USE NORMALIZED'
      CALL PRS(
     $     'Scale Factors: Input thru Mouse, type, or use normalized?$')
      IFCEN=.TRUE.
C     
C     10/7/92 pff .... always use MOUSE
C     
c     CALL MENU(XMOUSE,YMOUSE,BUTTON,'SCALE FACTOR')
c     IF(IFNOSEG)CALL DRCOVR(13)
      CHOICE=ITEM(1)
c     
c     
c     
      IFCEN=.FALSE.
C     Begin using grid latch filter
      IFGRID=.TRUE.
      IF(CHOICE.EQ.'TYPE') THEN
         CALL PRS('Enter Scale Factors Xfac,Yfac:$')
         CALL PRS(
     $        'They are the lengths that correspond to 1 scr height$')
         CALL RERR(XFAC,YFAC)
         IF(IFLEARN)WRITE(3,*)XFAC,YFAC
         XFAC=ABS(XFAC)
         YFAC=ABS(YFAC)
C     
         CALL PRS
     $        ('Enter the coordinates of Left Bottom corner of screen$')
         CALL RERR(XZERO,YZERO)
         IF(IFLEARN)WRITE(3,*)XZERO,YZERO
C     Mark Zero
         CALL COLOR(7)
      ELSE IF(CHOICE.EQ.'MOUSE') THEN
 17      CALL PRS('Push left button with mouse at 2 different'//
     $        ' x (horizontal)$')
         CALL PRS('locations$')
         do 18 I=1,3
            IF(I.EQ.2)CALL PRS('2nd x location$')
            IF(I.EQ.4)CALL PRS('2nd y location$')
            IF(I.EQ.3)THEN
               CALL PRS('Enter corresponding x length$')
               CALL GINGRD(0.0)
               call rer(XLEN)
               CALL GINGRD(GRID)
c              CALL PRS('Push left button with mouse at 2'//
c    $              ' different y (vertical) locations$')
c              CALL PRS('or click in menu area for equal X-Y scaling.$')
            ENDIF
            if (i.lt.3) CALL MOUSE(XSC(I),YSC(I),BUTTON)
C     
C     Check for keypad entry?
C     
c           IF (IFMENU(XSC(I),YSC(I)).AND.I.LE.2) THEN
c              CALL PRS('Please enter two points in build area.$')
c              goto 17
c           ENDIF
c           IF (IFMENU(XSC(I),YSC(I)).AND.I.GT.2) THEN
            IF (I.GT.2) THEN
               CALL PRS('Assuming equal scaling in X and Y.$')
               YSC(4)=XSC(2)
               YSC(3)=XSC(1)
               YLEN=XLEN
               ICOUNT=I-1
               goto 19
            ENDIF
C     
            CALL COLOR(7)
            CALL MOVE(XSC(I)-.02,YSC(I))
            CALL DRAW(XSC(I)+.02,YSC(I))
            CALL MOVE(XSC(I),YSC(I)-.02)
            CALL DRAW(XSC(I),YSC(I)+.02)
            ICOUNT=I
 18      CONTINUE
         CALL PRS('Enter in corresponding y length$')
         CALL GINGRD(0.0)
         call rer(YLEN)
         CALL GINGRD(GRID)
 19      CONTINUE
         IF(XSC(2).EQ.XSC(1) .OR. YSC(4).EQ.YSC(3))THEN
            CALL PRS('ERROR: MUST BE NONZERO GAP.  START OVER$')
            GO TO 3
         ENDIF
C     
C     Compute and check scale factors
C     
         XFACT=ABS(XLEN/(XSC(2)-XSC(1)))
         YFACT=ABS(YLEN/(YSC(4)-YSC(3)))
         IF(XFACT.EQ.0.0 .OR. YFACT.EQ.0.0)THEN
            CALL PRS('ERROR: MUST BE NONZERO SCALE FACTORS '//
     $           'START OVER$')
            GO TO 3
         ENDIF
C     
         do 20 I=1,ICOUNT
C     Erase markers
            CALL COLOR(0)
            CALL MOVEC(XSC(I)-.02,YSC(I))
            CALL DRAWC(XSC(I)+.02,YSC(I))
            CALL MOVEC(XSC(I),YSC(I)-.02)
            CALL DRAWC(XSC(I),YSC(I)+.02)
 20      CONTINUE
C     
C===============================================================
C     Set Origin
C===============================================================
C     
         CALL PRS (' $')
         CALL PRS
     $        ('Push left button with mouse at known point in x,y$')
         CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
C     
         CALL PRS('Enter X coordinate of point:$')
         CALL GINGRD(0.0)
         call rer(XPHYS)
         CALL PRS('Enter Y coordinate of point:$')
         call rer(YPHYS)
         CALL GINGRD(GRID)
C     
         XFAC  = ABS(XLEN/(XSC(2)-XSC(1)))
         YFAC  = ABS(YLEN/(YSC(4)-YSC(3)))
         XZERO = XPHYS - XFAC * XMOUSE
         YZERO = YPHYS - YFAC * YMOUSE
C===============================================================
      ELSE IF(CHOICE.EQ.'USE NORMALIZED') THEN
         XFAC=  2.0
         YFAC=  2.0
         XZERO =-1.
         YZERO =-1.
      ELSE
         CALL PRS('Error Reading Input$')
         GO TO 3
      ENDIF
C     
C===============================================================
C     Set grid stuff
C===============================================================
      GRIDDX=XFAC*GRID
      GRIDDY=YFAC*GRID
C     Set clipping window based on scale factors
      WT = YPHY(.995)
      WB = YPHY(.005)
      WL = XPHY(.005)
      WR = XPHY(.995)
C     Save scale factors for UnZoom function
      XFACO=XFAC
      YFACO=YFAC
      XZEROO=XZERO
      YZEROO=YZERO
C     
      CALL REFRESH
      CALL DRGRID
      return
      end

c-----------------------------------------------------------------------
      subroutine imp_mesh(ifquery_displace)

      include 'basics.inc'
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)

      logical ifquery_displace,ifdisplace
      logical iflow,iheat

      real    xyzbox(6)

      call prs('input name of new .rea file$')
      call blank(fname,70)
      call res  (fname,70)

      ifdisplace = .false.
      if (nel.gt.0.and.ifquery_displace) then
       call prs('Would you like to displace existing elements in box?$')
       call res  (ans,1)
       if (ans.eq.'y'.or.ans.eq.'Y') ifdisplace = .true.
      endif

      if (indx1(fname,'.rea',4).eq.0) then !  Append .rea, if not present
         len = ltrunc(fname,70)
         call chcopy (fnam1(len+1),'.rea',4)
      endif
      open(unit=47,file=fname,status='old',err=1000)

c     Successfully opened file, scan until 'MESH DATA' is found
      call readscan('DIMENSION',9,47)

      do ipass=1,2              ! Read parameters and passive scalar data
         read(47,*) mpar
         do k=1,mpar
            read(47,'(a1)') ans ! dummy read
         enddo
      enddo

      read(47,*) nlogic
      call get_flow_heat(iflow,iheat,nlogic,47) ! iflow/iheat _local_ logicals

      call readscan('MESH DATA',9,47)
      read (47,*) nelin,ndimn
      neln = nelin + nel

      if (ndimn.ne.ndim) then
         call prs('Dimension must match dimension of current session.$')
         call prs('Returning.$')
         close(47)
         return
      endif

      if (neln.ge.nelm) then
         call prs('Sorry, number of elements would exceed nelm.$')
         call prs('Returning.$')
         call prii(neln,nelm)
         close(47)
         return
      endif

c     Read geometry

      nelo = nel
      nels = nel+1

      ierr=imp_geom(x(1,nels),y(1,nels),z(1,nels),nelm
     $             ,numapt(nels),letapt(nels),igroup(nels)
     $             ,ndim,nelin,47)
      if (ierr.ne.0) then
         call prs('Error reading geometry... returning.$')
         close(47)
         return
      endif
c
      ierr=imp_curv(nc,ccurve,curve,ndim,nelin,nel,47)
      if (ierr.ne.0) then
         call prs('Error reading curve side info... returning.$')
         close(47)
         return
      endif
c
c
c     Read BC info
c
c
      read(47,3) d
    3 format(a3)
c
      ifld0 = 1
      if (.not.iflow) ifld0=2
      write(6,*) 'IFLD:',ifld0,nflds,iflow,iheat
c
      do ifld=ifld0,nflds
         if (.not.iflow) read(47,*) ans  ! dummy read
         ierr=imp_bc(cbc(1,nels,ifld),bc(1,1,nels,ifld),ibc(1,nels,ifld)
     $      ,ndim,nelin,nel,47)
         if (ierr.ne.0) then
            call prsii('nelin,ifld:$',nelin,ifld)
            call prs('Error reading boundary conditions. Returning.$')
            close(47)
            return
         endif
      enddo
c
c     ALL DONE
c
      close(47)
c
      write(6,*) 'This is nel,ncurve old:',nel,ncurve
      nel = neln
      ncurve = ncurve + nc
      write(6,*) 'This is nel,ncurve new:',nel,ncurve
c
      call gencen
c
c
      if (ifdisplace) then
c
         ncrnr = 2**ndim
c
         xyzbox(1) = 9.e21
         xyzbox(2) =-9.e21
         xyzbox(3) = 9.e21
         xyzbox(4) =-9.e21
         xyzbox(5) = 9.e21
         xyzbox(6) =-9.e21
c
         do ie=nels,nel   ! Find box containing imported mesh
         do i=1,ncrnr
            xyzbox(1) = min(xyzbox(1),x(i,ie))
            xyzbox(2) = max(xyzbox(2),x(i,ie))
            xyzbox(3) = min(xyzbox(3),y(i,ie))
            xyzbox(4) = max(xyzbox(4),y(i,ie))
            xyzbox(5) = min(xyzbox(5),z(i,ie))
            xyzbox(6) = max(xyzbox(6),z(i,ie))
         enddo
         enddo
c
         write(s,110) 'X',xyzbox(1),xyzbox(2)
         call prs(s)
         write(s,110) 'Y',xyzbox(3),xyzbox(4)
         call prs(s)
         if (if3d) write(s,110) 'Z',xyzbox(5),xyzbox(6)
         if (if3d) call prs(s)
  110    format(' ',a1,' min/max of new elements:',2g15.7,'$')
c
         call substitute_el(xyzbox,nelo)
      endif
c
      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      close(47)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine readscan(key,len_key,io)
c
c     Read a file until "key" is found or eof is found.
c
c
      character*80 key
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
   80 format(a80)
   81 format(80a1)
c
      do i=1,90000
         call blank(string,80)
         read (io,80,end=100,err=100) string
         len = ltrunc(string,80)
         if (indx1(string,key,len_key).ne.0) return
      enddo
  100 continue
      return
      end
c-----------------------------------------------------------------------
      function imp_geom(x,y,z,ld,numapt,letapt,igroup,ndim,nel,io)
      real        x(8,ld),y(8,ld),z(8,ld)
      integer     numapt(1),igroup(1)
      character*1 letapt(1)
c
      character*1 c
c
c
c
      if (ndim.eq.3) then
         do ie=1,nel
c
c           read(io,1) c
c   1       format(a1)
c
c           read(io,'(20x,i4,4x,i3,a1,11x,i5)',err=9,end=9)
c    $      idum,numapt(ie),letapt(ie),igroup(ie)
            read(io,*)
            numapt(ie) = 1
            letapt(ie) = 'A'
            igroup(ie) = 0
c
            read(io,*,end=9,err=9) (x(k,ie),k=1,4)
            read(io,*,end=9,err=9) (y(k,ie),k=1,4)
            read(io,*,end=9,err=9) (z(k,ie),k=1,4)
c
            read(io,*,end=9,err=9) (x(k,ie),k=5,8)
            read(io,*,end=9,err=9) (y(k,ie),k=5,8)
            read(io,*,end=9,err=9) (z(k,ie),k=5,8)
c
         enddo
      else
         do ie=1,nel
c
c           read(io,1) c
c
c           read(io,'(20x,i4,4x,i3,a1,11x,i5)',err=9,end=9)
c    $      idum,numapt(ie),letapt(ie),igroup(ie)
            read(io,*)
            numapt(ie) = 1
            letapt(ie) = 'A'
            igroup(ie) = 0
c
            read(io,*,end=9,err=9) (x(k,ie),k=1,4)
            read(io,*,end=9,err=9) (y(k,ie),k=1,4)
         enddo
      endif
c
      imp_geom = 0
      return
c
    9 continue
      imp_geom = 1
      return
      end
c-----------------------------------------------------------------------
      function imp_curv(nc,cc,c,ndim,nel,nelo,io)
      character*1 cc(12,1)
      real        c(6,12,1)
c
      character*1 d
      logical iffmtin
c
c
c
      iffmtin = .true.
c
      read(io,1) d
    1 format(a1)
c
      read(io,*) nc
      write(6,*) 'Found',nc,' curve sides.'
 
      do i=1,nc

         if (nel.lt.1000.AND.iffmtin) then
            read(io,'(2i3,5g14.6,1x,a1)',err=9,end=9)
     $                    iedge,ie,r1,r2,r3,r4,r5,d
         elseif (nel.lt.1 000 000.and.iffmtin) THEN
            read(io,'(i2,i6,5g14.6,1x,a1)',err=9,end=9)
     $                    iedge,ie,r1,r2,r3,r4,r5,d
         elseif (iffmtin) THEN
            read(io,'(i2,i12,5g14.6,1x,a1)',err=9,end=9)
     $                    iedge,ie,r1,r2,r3,r4,r5,d
         else
            read(io,err=9,end=9) iedge,ie,r1,r2,r3,r4,r5,d
         endif

         ie = nelo + ie
         c(1,iedge,ie)=r1
         c(2,iedge,ie)=r2
         c(3,iedge,ie)=r3
         c(4,iedge,ie)=r4
         c(5,iedge,ie)=r5
         c(6,iedge,ie)=0.
         cc( iedge,ie)=d
      enddo
c
      imp_curv = 0
      return
c
    9 continue
      imp_curv = 1
      return
      end
c-----------------------------------------------------------------------
      function imp_bc(cbc,bc,ibc,ndim,neln,nel,io)
      character*3  cbc   (6,1)
      real          bc (5,6,1)
      integer      ibc (6,1)
      real*8       bc8 (5)

      character*80 a80
      character*1  a81(1)
      equivalence(a81,a80)

      imp_bc = 0

      call blank (a80,80)
      read(io,80) a80
   80 format(a80)

      call capit(a80,80)
      len = ltrunc(a80,80)
      a81(len+1) = '$'
      call prs(a80)
      write(6,80) a80

      i1 = indx1(a80,'NO',2)
c     write(6,*) i1,' indx1 '
      if (i1.ne.0) return
c     if (indx1(a80,'NO',2).ne.0) return
c
      nsides = 2*ndim
      nbcrea = 5

      do ie=1,neln
      do iside = 1,nsides
c        write(6,*) ie,iside,neln,' READING BC'
         if (neln.lt.1000) then
            read(io,'(1x,a3,1X,2i3,5g14.6)',err=9,end=9)
     $      cbc(iside,ie),id,jd,(bc8(ii),ii=1,nbcrea)
         elseif (neln.lt.1 000 000) then
            read(io,'(1x,a3,i6,5g14.6)',err=9,end=9)
     $      cbc(iside,ie),id,(bc8(ii),ii=1,nbcrea)
         else
            read(io,'(1x,a3,i12,5g18.11)',err=191,end=191)
     $      cbc(iside,ie),id,(bc8(ii),ii=1,nbcrea)
         endif
         goto 192

  191    continue
            if (iside.eq.1.and.mod(ie,100 000).eq.0) write(6,*)
     $      'Warning: reading bc in old format! Check results!'
            write(6,*) 'Warning: input bc in old format! Check results!'
            read(io,'(1x,a3,i6,5g14.6)',err=9,end=9)
     $      cbc(iside,ie),id,(bc8(ii),ii=1,nbcrea)
  192    continue

         do ii=1,nbcrea
            bc(ii,iside,ie) = bc8(ii)
            ibc(iside,ie)   = bc8(1)   ! High precision for P bc
         enddo
c        Adjust periodic boundary conditions
         if (cbc(iside,ie).eq.'P  ') bc(1,iside,ie) = bc(1,iside,ie)+nel
         if (cbc(iside,ie).eq.'P  ') ibc(iside,ie) =  ibc(iside,ie)+nel
      enddo
      enddo
c
      imp_bc = 0
      return
c
    9 continue
      write(6,*) 'made it to:',ie,iside
      imp_bc = 1
      return
      end
c
c-----------------------------------------------------------------------
      subroutine new_out(fname,lname,gv,ncrnr,xyz,ind)
c
c     Read a file until "key" is found or eof is found.
c
      include 'basics.inc'

      parameter (maxv = 8*nelm)
      real      xp(maxv),yp(maxv),zp(maxv)
      common /c_xyz/ xp,yp,zp
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      integer gv(ncrnr,1),ind(1)
      real    xyz(1)
c
      integer ecrnr(8)
      save    ecrnr
      data    ecrnr  / 1,2,4,3,5,6,8,7 /
c
   80 format(a80)
   81 format(80a1)
c
      call blank(string,80)
      call chcopy(string,fname,lname)
      call chcopy(string1(lname+1),'.map',4)
c
c     Get current vertex map info
c
      open(unit=20,file=string,err=999)
      read(20,*) neln
      do ie=1,neln
         read(20,*,end=998,err=998) idum,(gv(k,ie),k=1,ncrnr)
      enddo
      close(unit=20)
c
c
c     Sort data and gridpoints by global vertex number
c
      l = 1
      do ie = 1,neln
         do k=1,ncrnr
            j     = ecrnr(k)
            xp(l) = x(j,ie)
            yp(l) = y(j,ie)
            zp(l) = z(j,ie)
            l=l+1
         enddo
      enddo
c
      nv = ncrnr*neln
      call irank(gv,ind,nv)
      call swap (xp,xyz,ind,nv)
      call swap (yp,xyz,ind,nv)
      call swap (zp,xyz,ind,nv)
c
c     Compress lists
c
      nnv=1
      do i=2,nv
         if (gv(i,1).ne.gv(nnv,1)) then
            nnv       = nnv+1
            gv(nnv,1) = gv(i,1)
            xp(nnv)   = xp(i)
            yp(nnv)   = yp(i)
            zp(nnv)   = zp(i)
         endif
      enddo
      nv = nnv
c
c     Output vertex info
c
      l=1
      if (ndim.eq.3) then
         do i=1,nv
            xyz(l  ) = xp(i)
            xyz(l+1) = yp(i)
            xyz(l+2) = zp(i)
            l=l+3
         enddo
      else
         do i=1,nv
            xyz(l  ) = xp(i)
            xyz(l+1) = yp(i)
            l=l+2
         enddo
      endif
      npts_o = nv*ndim
c
      call chcopy(string1(lname+1),'.vtx',4)
      open(unit=11,file=string,err=999)
      write(11,*) npts_o
      if (ndim.eq.3) then
         write(11,33) (xyz(i),i=1,npts_o)
      else
         write(11,32) (xyz(i),i=1,npts_o)
      endif
   32 format(1p2e15.7)
   33 format(1p3e15.7)
      return
c
  998 continue
      call prs('end of .map file reached')
      return
c
  999 continue
      call prs('Could not open .map file.')
      return
      end
c-----------------------------------------------------------------------
      function i_finds(list,ig,n)
c     Find if there is a match of ig in *sorted* list
      integer i_finds,ig,list(1),n
      integer hi,lo,m
c
      lo = 1
      hi = n
      icount  = 0
      i_finds = 0
c
    1 continue 
      m = (hi+lo)/2
      if (ig.eq.list(m)) then
         i_finds = m
         return
      elseif (ig.lt.list(m)) then
         hi = m
      else
         lo = max(m,lo+1)
      endif
      icount = icount+1
      if (icount.gt.n.or.lo.gt.hi) then
         write(6,*) 'ERROR IN i_finds:',ig,n,hi,lo,m,' ABORT'
         j = i_findu(list,ig,n)
         i_finds = j
         write(6,*) 'result of i_findu:',j,list(j)
         call exitt
         return
      else
         goto 1
      endif
      return
      end
c-----------------------------------------------------------------------
      function i_findu(list,ig,n)
c     Find if there is a match of ig in *unsorted* list
      integer i_findu,ig,list(1),n
      integer hi,lo,m
c
      i_findu = 1
c
      do m=1,n
         if (ig.eq.list(m)) then
            i_findu = m
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      function i_findu_zero(list,ig,n)
c     Find if there is a match of ig in *unsorted* list	starting from 0
      integer i_findu,ig,list(1),n
      integer hi,lo,m
c
      i_findu_zero = 0
c
      do m=1,n
         if (ig.eq.list(m)) then
            i_findu_zero = m-1
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_vtx
c
c     Read vtk-like unstructured mesh format
c
      include 'basics.inc'
c
      parameter (maxv = 8*nelm)
      integer        vv(8,nelm),vnum(maxv)
      common /c_vtx/ vv,vnum

      real      xp(maxv),yp(maxv),zp(maxv)
      common /c_xyz/ xp,yp,zp

      integer      kcell(8),indx(8)
      save         indx
      data         indx  / 1,2,4,3,5,6,8,7 /
C
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)
c
      call prsi('Current number of elements:$',nel)
      call prs('input name of new vertex/cell file$')
      call blank(fname,70)
      call res  (fname,70)
      ikill = indx1(fname,' ',1)
      nkill = 70-ikill+1
      call blank(fnam1(ikill),nkill)
c
      open(unit=47,file=fname,status='old',err=1000)
c
c     Read vertex info
c
      read(47,*) nvtx 
c
      if (nvtx.gt.maxv) then
         nvtx = nvtx+1
         call prsii('Too many pts. Inc. maxv to:$',nvtx,nelm)
         close(47)
         return
      endif
c
      if (if3d) then
        do i=1,nvtx
c        read(47,*) vnum(i),xp(i),yp(i),zp(i)
         read(47,*) xp(i),yp(i),zp(i)
         vnum(i) = i
        enddo
      else
        do i=1,nvtx
c        read(47,*) vnum(i),xp(i),yp(i)
         read(47,*) xp(i),yp(i)
         vnum(i) = i
        enddo
      endif
      call prsis('Found$',nvtx,' points.$')
      xmin = glmin(xp,nvtx)
      xmax = glmax(xp,nvtx)
      ymin = glmin(yp,nvtx)
      ymax = glmax(yp,nvtx)
      call prsrr('xmin xmax:$',xmin,xmax)
      call prsrr('ymin ymax:$',ymin,ymax)
c
c
c     Read cell data
c
      read(47,*) ncell
c
      if (ncell+nel.gt.nelm) then
         write(6,*) ncell,nel,nelm,' ncell,nel,nelm'
         ncell = ncell+nel
         call prsii('Too many elements. Increase nelm to:$',ncell,nelm)
         close(47)
         return
      endif

      npt=2**ndim
      mcell = 0
      do ie=1,ncell
         read(47,*) icell,(kcell(k),k=1,npt)
         if (icell.eq.npt) then
            mcell = mcell+1
            do k=1,npt
               kcell(k) = kcell(k) + 1   ! for vtk, need to add 1
               vv(k,ie) = kcell(k)

c              rrr = .05
c              ii = vv(k,ie)
c              call draw_circ(xp(ii),yp(ii),zp(ii),rrr)
c              write(6,*) ie,k,vv(k,ie),kcell(k),nvtx
c              vv(k,ie) = i_finds(vnum,kcell(k),nvtx)

            enddo
         endif

c        call prsi('continue? $',ie)
c        call res(ans,1)
         numapt(ie)=1
         letapt(ie)='A'
      enddo
  444 ncell = mcell

      call prs('continue?$')
      call res(ans,1)

c
c     At this point, we have the cell data in the "new" format.
c
c     Now we chuck it in favor of the worthless nekton format.
c
c     4/17/01 -- first, try to figure out the local vertex ordering
c
c
c     call prs('Flip for Gambit? (Yes=Gambit, No=Seung)$')
c     call res(ans,1)
      ans='y'
c
c
      if (ans.eq.'y'.or.ans.eq.'Y') then
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               x(k,nel)=xp(l)
               y(k,nel)=yp(l)
               z(k,nel)=zp(l)
            enddo
         enddo
      else                     ! flip
         do ie=1,ncell
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               j = indx(k)
               x(j,nel)=xp(l)
               y(j,nel)=yp(l)
               z(j,nel)=zp(l)
            enddo
         enddo
      endif
c     do ie=1,ncell
c     do k=1,npt
c        write(6,*) k,ie,x(k,ie),y(k,ie)
c     enddo
c     enddo
      if_unstructure=.true.
      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      call prs(fname)
      close(47)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine imp_mesh_vtk
c
c     Read vtk unstructured mesh format
c
      include 'basics.inc'
c
      parameter (maxv = 8*nelm)
      integer        vv(8,nelm),vnum(0:maxv-1)
      common /c_vtx/ vv,vnum

      real      xp(0:maxv-1),yp(0:maxv-1),zp(0:maxv-1)
      common /c_xyz/ xp,yp,zp

      integer      kcell(8),indx(8)
      save         indx
      data         indx  / 1,2,4,3,5,6,8,7 /
C
      character*3 d
      character*1  fnam1(70)
      character*70 fname
      equivalence (fname,fnam1)
c
      character*1  s81(80)
      character*80 s80
      equivalence (s80,s81)
c
      call prsi('Current number of elements:$',nel)
      call prs('input name of vtk file$')
      call blank(fname,70)
      call res  (fname,70)
      ikill = indx1(fname,' ',1)
      nkill = 70-ikill+1
      call blank(fnam1(ikill),nkill)
c
      open(unit=47,file=fname,status='old',err=1000)
c
c     Successfully opened file, scan until 'MESH DATA' is found
c
c     Read vertex info
c
      read(47,*) 
      read(47,*) 
      read(47,*) 
      read(47,*) 
 
      call blank(s80,80)
      read(47,80) s80
      ib = indx1(s80,' ',1)
      open (unit=48,file='my_vtk.tmp')
      write(48,81) (s81(k),k=ib+1,80)
      rewind (unit=48)
      read (48,*) nvtx
      close (unit=48)
c
      if (nvtx.gt.maxv) then
         nvtx = nvtx+1
         call prsii('Too many pts. Inc. maxv to:$',nvtx,nelm)
         close(47)
         return
      endif
c
      if (if3d) then
        do i=0,nvtx-1
         read(47,*) xp(i),yp(i),zp(i)
         write(6,*) xp(i),yp(i),zp(i),i
         vnum(i) = i
        enddo
      else
        do i=0,nvtx-1
         read(47,*) xp(i),yp(i)
         vnum(i) = i
        enddo
      endif
c
c
c     Read cell data
c
      read(47,*)
      call blank(s80,80)
      read(47,80) s80
      write(6,*) s80
      ib = indx1(s80,' ',1)
      open (unit=48,file='my_vtk.tmp')
      write(48,81) (s81(k),k=ib+1,80)
      rewind (unit=48)
      read (48,*) ncell
      close (unit=48)
c
c
      if (ncell+nel.gt.nelm) then
         ncell = ncell+nel
         call prsii('Too many elements. Increase nelm to:$',ncell,nelm)
         close(47)
         return
      endif

      nelo = 0
      ncrd = 0 ! # of cell lines read
      npt=2**ndim
      do ie=1,ncell
         call blank(s80,80)
         read(47,80) s80
         read(s80,1) icell_type
   1     format(i1)
         if (icell_type.eq.npt) then   ! =4 for quad/2D, =8 for hex/3D
            ib = indx1(s80,' ',1)
            open (unit=48,file='my_vtk.tmp')
            write(48,81) (s81(k),k=ib+1,80)
            rewind (unit=48)
            read (48,*) (kcell(k),k=1,npt)
            close (unit=48)

            ncrd = ncrd + 1

            do k=1,npt
               vv(k,ncrd) = i_findu_zero(vnum,kcell(k),nvtx)

               l=vv(k,ncrd)
               write(6,*) k,ie,ncrd,l,xp(l),yp(l),'  cell numbers'
            enddo
         endif
      enddo

c     Check that only CELL_TYPES 1,3,8,9,11,12 are in the .vtk file (i.e. point,line,quad,hex)?
 
      close (unit=47)
c
c     At this point, we have the cell data in the "new" format.
c
c     Now we chuck it in favor of the worthless nekton format.
c
c     4/17/01 -- first, try to figure out the local vertex ordering
c
c
c     call prs('Flip for Gambit? (Yes=Gambit, No=Seung)$')
c     call res(ans,1)
      ans = 'y'
c
c
      if (ans.eq.'y'.or.ans.eq.'Y') then
c        do ie=1,ncell
         do ie=1,ncrd
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               x(k,nel)=xp(l)
               y(k,nel)=yp(l)
               z(k,nel)=zp(l)
            enddo
         enddo
      else                     ! don't flip
         do ie=1,ncrd
            nel = nel+1
            do k=1,npt
               l = vv(k,ie)
               j = indx(k)
               x(j,nel)=xp(l)
               y(j,nel)=yp(l)
               z(j,nel)=zp(l)
            enddo
         enddo
      endif
      if_unstructure=.true.

      call set_std_bcs(nelo)


      return
c
 1000 continue
      call prs('Unable to open file.  Returning.$')
      call prs(fname)
      close(47)
c
   80 format(a80)
   81 format(80a1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine check_spj
c
c     Check for split (parent) - join (child) nonconforming interfaces.
c
      include 'basics.inc'
      logical ifclose,ifok
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /
c
c     Currently designed only for 2D, NO CURVE nonconf.
c
      if (ndim.eq.3) return


      ifld0=1
      ifld1=nflds

      if (ifconj_merge) then  ! Update THERMAL BCS  ONLY
         ifld0=2
         ifld1=2
      endif


      nsides = 2*ndim
      do kf=ifld0,ifld1
      do ie=1,nel
       if (maskel(ie,kf).gt.0) then
        do is=1,nsides
c         write(6,*) 'CBC 1:',cbc(is,ie,kf),ie,is
          if (cbc(is,ie,kf).ne.'E  ') then   ! n^2 loop time!
            do je=1,nel
             if (maskel(je,kf).gt.0) then
              do js=1,nsides
               if (cbc(js,je,kf).ne.'E  '.and.ie.ne.je) then
                  if (ifclose(ie,is,je,js)) then
c                    Currently designed only for 2D, NO CURVE nonconf.
                     if (if3d) then
                       call chk_segs(ifok,icond,ie,is,je,js)
                       if (ifok) then
                          cbc(is,ie,kf)='SP '
                          bc (5,is,ie,kf)=bc (5,is,ie,kf)+1
                          cbc(js,je,kf)='J  '
                          ibc(js,je,kf)=ie
                          bc (1,js,je,kf)=ie
                          bc (2,js,je,kf)=is
c                         write(6,*) 'Match:',ie,is,je,js,icond
c                       write(6,*) 'Match:',cbc(is,ie,kf),cbc(js,je,kf)
                       else
                       endif
                     else
                       call chk_seg(g1,e1,g2,e2,ie,is,je,js)
                       if (abs(e1).le.0.001          .and.
     $                     abs(e2).le.0.001          .and.
     $                     -.01.le. g1.and.g1.le.1.01 .and.
     $                     -.01.le. g2.and.g2.le.1.01) then  ! segment
                          cbc(is,ie,kf)='SP '
                          bc (5,is,ie,kf)=bc (5,is,ie,kf)+1
                          cbc(js,je,kf)='J  '
                          ibc(js,je,kf)=ie
                          bc (1,js,je,kf)=ie
                          bc (2,js,je,kf)=is
                          bc (3,js,je,kf)=g1
                          bc (4,js,je,kf)=g2
                       endif
                     endif
                  endif
               endif
              enddo
             endif
            enddo
         endif
c        write(6,*) 'CBC 2:',cbc(is,ie,kf),ie,is
        enddo
       endif
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      logical function ifclose(ie,is,je,js)
c
c
      include 'basics.inc'
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /

c
      ifclose = .true.
c
      nv=2**(ndim-1)
      dz=0.
      do jv=1,nv
         dx = sides(ie,is,1)-x(jvs(jv,js),je)
         dy = sides(ie,is,2)-y(jvs(jv,js),je)
         if (if3d) dz = sides(ie,is,3)-z(jvs(jv,js),je)
         ds = dx*dx+dy*dy+dz*dz
         if (ds.le.rcen(ie)*rcen(ie)) return
      enddo
      ifclose = .false.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_seg(g1,e1,g2,e2,ie,is,je,js)
c
c     Upon return,  (g1,e1) and (g2,e2) are the transformed
c                   coordinates of the two vertices associated with
c                   side (je,js) in the coordinate system defined
c                   by the two vertices associated with side (ie,is).
c
c     For example, if side (je,js) is *coincident* with (ie,is), then either
c                  (g1,e1) = (0.,0.) and  (g2,e2) = (1.,0.),  or
c                  (g1,e1) = (1.,0.) and  (g2,e2) = (0.,0.).
c
c
      include 'basics.inc'
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,5,6 , 2,3,6,7 , 3,4,7,8
     $            , 4,1,8,5 , 1,2,3,4 , 5,6,7,8 /
c
c
c     find the (x,y) coordinates of the two endpoints of side (ie,is)
      x1 = x(jvs(1,is),ie)
      y1 = y(jvs(1,is),ie)
      x2 = x(jvs(2,is),ie)
      y2 = y(jvs(2,is),ie)
c
c     find the (x,y) coordinates of the two endpoints of side (je,js)
      x3 = x(jvs(1,js),je)
      y3 = y(jvs(1,js),je)
      x4 = x(jvs(2,js),je)
      y4 = y(jvs(2,js),je)
c
c     use these endpoints to set up the transformation.
c     the 3x3 matrix is not explicitly formed, just the closed 
c     form result.
c
      x2_x1 = x2 - x1
      y2_y1 = y2 - y1
      hyp2 = ( x2_x1*x2_x1 + y2_y1*y2_y1 )
c
c
      scale  = 1./hyp2
      g1 = scale*( x3*x2_x1 + y3*y2_y1 - x1*x2_x1 - y1*y2_y1 )
      e1 = scale*(-x3*y2_y1 + y3*x2_x1 + x1*y2_y1 - y1*x2_x1 )
      g2 = scale*( x4*x2_x1 + y4*y2_y1 - x1*x2_x1 - y1*y2_y1 )
      e2 = scale*(-x4*y2_y1 + y4*x2_x1 + x1*y2_y1 - y1*x2_x1 )
c
c     write(8,1) ie,is,je,js,g1,e1,g2,e2
c   1 format(4i4,1p4e12.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_segs(ifok,icond,ie,is,je,js)
c
c     3D routine for determining if face(je,js) is a subset of
c     face(ie,is)
c
c     Currently predicated on assumption that faces are PLANAR.
c
c
      include 'basics.inc'
c
      logical ifok,ifxoeq0
c
      integer jvs(4,6)
      save    jvs
      data    jvs / 1,2,6,5 , 2,3,7,6 , 3,4,8,7
     $            , 4,1,5,8 , 4,3,2,1 , 5,6,7,8 /  ! EB cclockws ordering
c
      real n(3),nhat(3),xx(3,4),xj(3),xo(3),d1(3),d2(3),norm
      real vnt(3),vv(3),ee(3),op(3),xjp(3)
c
c     First, get canonical coordinates of plane containing face i
c
c     if (ie.eq.1.and.je.gt.1.and.is.eq.6.and.js.eq.5)
c    $    write(6,*) 'ie,is,je,js:',ie,is,je,js
c     write(6,*) 'Continue? (any key)'
c     read (5,*) ans
c
      do i=1,4
         xx(1,i) = x(jvs(i,is),ie)
         xx(2,i) = y(jvs(i,is),ie)
         xx(3,i) = z(jvs(i,is),ie)
      enddo
c     write(6,1) (i,(xx(k,i),k=1,3),' xx',i=1,4)
    1 format(/,4(i4,3f10.4,a3,/))
c
      do k=1,3
        xo(k) = 0.
        do i=1,4
           xo(k) = xo(k)+xx(k,i)
        enddo
        xo(k) = 0.25*xo(k)
      enddo
c     write(6,1) i,(xo(k),k=1,3),' xo'
c
      call sub3(d1,xx(1,3),xx(1,1),3) ! vectors connecting diagonally
      call sub3(d2,xx(1,2),xx(1,4),3) ! opposed vertices of face i
      call vcross_normal(nhat,sine,d2,d1)
c     write(6,1) i,(nhat(k),k=1,3),' nh'
      d1n = dotprod(d1,d1)
      d2n = dotprod(d2,d2)
      clen = max(d1n,d2n)
      clen = sqrt(clen)
c
c     Generate basis vectors for plane i
c
      call norm3d(d1)
      call vcross_normal(d2,sine,nhat,d1)
c     write(6,1) i,(d1(k),k=1,3),' d1'
c     write(6,1) i,(d2(k),k=1,3),' d2'
c
c     Now for each vertex of face j, check if it lies on the plane of face i
c     and if it's on the edge of or interior to face i.
c
      ifok  = .false.
      icond = 0
      do j=1,4
         xj(1) = x(jvs(j,js),je)
         xj(2) = y(jvs(j,js),je)
         xj(3) = z(jvs(j,js),je)
c        write(6,1) j,(xj(k),k=1,3),' xj'
c
c        map to canonical plane of face i
c
         call sub3(xjp,xj,xo,3)
         xp1 = dotprod(d1  ,xjp)
         xp2 = dotprod(d2  ,xjp)
         xp3 = dotprod(nhat,xjp)
         xjp(1) = xp1
         xjp(2) = xp2
         xjp(3) = xp3
c        write(6,1) j,(xjp(k),k=1,3),'xjp'
c
c        Now check for any disqualifications 
c
         icond = 2
         epsr  = 1.e-4
         rel_err = abs(xjp(3))/clen
c        write(6,*) j,rel_err,xjp(3),clen,'rel_err'
         if (rel_err.gt.epsr) return
c
c        Vertex is co-planar with face i, check if interior to face i
c
         epsn = 1.e-4*clen
         sign = 0.
         do i=1,4
            i1 = i+1
            i1 = mod1(i1,4)
            call sub3(ee,xx(1,i1),xx(1,i),3)
            call sub3(vv,xj(1)   ,xx(1,i),3)
c
            vnorm = dotprod(vv,vv)
            if (vnorm.gt.0) vnorm = sqrt(vnorm)
c
            if (vnorm.gt.epsn) then
               call norm3d(vv)
               call norm3d(ee)
               scal = 1.0 - abs(dotprod(vv,ee))  ! check for co-linearity
               if (scal.gt.epsr) then
                  call cross(vnt,ee,vv)
                  dotp = dotprod(vnt,nhat)
                  if (sign.eq.0) then
                     sign = dotp
                  elseif (sign*dotp.lt.0) then
                     icond=3+i
                     return
                  endif
               endif
            else
               goto 11
            endif
         enddo
   11    continue
      enddo
      ifok  = .true.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine chk_right_hand(nl)
      include 'basics.inc'
c
      real xyz(2,4)
c
      do iel=1,nl
         do 210 I=1,4
            XYZ(1,I)= x(i,iel)
            XYZ(2,I)= y(i,iel)
  210    CONTINUE
C
C        CRSS2D(A,B,O) = (A-O) X (B-O)
C                      (note that the notation for corner number here
C                       differs slightly from that used in NEKTON the code.)
         C1=CRSS2D(XYZ(1,2),XYZ(1,4),XYZ(1,1))
         C2=CRSS2D(XYZ(1,3),XYZ(1,1),XYZ(1,2))
         C3=CRSS2D(XYZ(1,4),XYZ(1,2),XYZ(1,3))
         C4=CRSS2D(XYZ(1,1),XYZ(1,3),XYZ(1,4))
C
         IERR=0
         IF (C1.LE.0.0.AND.C2.LE.0.0.AND.
     $       C3.LE.0.0.AND.C4.LE.0.0 ) THEN
C            cyclic permutation (counter clock-wise):  reverse
             call prsi('Reversing Element: $',iel)
             x(2,iel) = xyz(1,4)
             y(2,iel) = xyz(2,4)
             x(4,iel) = xyz(1,2)
             y(4,iel) = xyz(2,2)
             IF (IF3D) THEN
                do 400 I=1,4
                   x(i+4,iel)=x(i,iel)
                   y(i+4,iel)=y(i,iel)
  400           CONTINUE
             ENDIF
         ELSEIF (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $           C3.LE.0.0.OR.C4.LE.0.0 ) THEN
             CALL PRSi('ERROR in entering element.  Re-enter.$',iel)
             IERR=1
         ENDIF
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine special_delete
      include 'basics.inc'
c
      integer dflag(nelm),e,emin,ecount
c
      open(unit=49,file='cyls.dat',status='old',err=999)
c
      read(49,*) ncyl
c
c     Delete all elements within radius
c
c
      call izero(dflag,nel)
c
      do k=1,ncyl
         read(49,*) rad,x0,y0
         r1 = .99*rad
         r1 = r1*r1
         do e=1,nel
         do i=1,4
            dx = x0-x(i,e)
            dy = y0-y(i,e)
            r2 = dx*dx + dy*dy
            if (r2.lt.r1) dflag(e) = 1
         enddo
         enddo
      enddo
      close(49)
c
      emin = nel+1
      do e=1,nel
         if (dflag(e).ne.0) then
            emin=e
            goto 10
         endif
      enddo
   10 continue
c
      ecount = emin-1
      do e=emin+1,nel
         if (dflag(e).eq.0) then
            ecount = ecount + 1
            call copyel(e,ecount)
         endif
      enddo
      ndel = nel - ecount
      call prsis('Special delete of$',ndel,' elements.$')
      nel  = ecount
c
  999 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine list_delete
      include 'basics.inc'
      character*80 dname

      common /idelt/ dflag(nelm)
      integer dflag,e,emin,ecount

      call prs('Enter name of file with elements to delete:$')
      call blank(dname,80)
      call res(dname,80)
      open(unit=49,file=dname,status='old',err=999)


      call izero(dflag,nel)

      emin = nel+1
      do e=1,nel
         read(49,*,end=10) i_delete
         if (i_delete.gt.0.and.i_delete.le.nel) then
            dflag(i_delete) = 1         ! Mark for deletion
            emin = min(emin,i_delete)
         endif
      enddo
   10 continue

      close(49)

      ecount = emin-1
      do e=emin+1,nel
         if (dflag(e).eq.0) then
            ecount = ecount + 1
            call copyel(e,ecount)
         endif
      enddo
      ndel = nel - ecount
      call prsis('List delete of$',ndel,' elements.$')
      nel  = ecount

  999 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine special_delete_cyls
      include 'basics.inc'
c
      integer dflag(nelm),e,emin,ecount
c
      open(unit=49,file='cyls.dat',status='old',err=999)
c
      read(49,*) ncyl
c
c     Delete all elements within radius
c
c
      call izero(dflag,nel)
c
      do k=1,ncyl
         read(49,*) rad,x0,y0
         r1 = .99*rad
         r1 = r1*r1
         do e=1,nel
         do i=1,4
            dx = x0-x(i,e)
            dy = y0-y(i,e)
            r2 = dx*dx + dy*dy
            if (r2.lt.r1) dflag(e) = 1
         enddo
         enddo
      enddo
      close(49)
c
      emin = nel+1
      do e=1,nel
         if (dflag(e).ne.0) then
            emin=e
            goto 10
         endif
      enddo
   10 continue
c
      ecount = emin-1
      do e=emin+1,nel
         if (dflag(e).eq.0) then
            ecount = ecount + 1
            call copyel(e,ecount)
         endif
      enddo
      ndel = nel - ecount
      call prsis('Special delete of$',ndel,' elements.$')
      nel  = ecount
c
  999 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine set_std_bcs(nelo)
      include 'basics.inc'
      integer e

      do e=nelo+1,nel
         cbc(1,e,1) = 'W  '
         cbc(2,e,1) = 'W  '
         cbc(3,e,1) = 'W  '
         cbc(4,e,1) = 'W  '
         cbc(5,e,1) = 'v  '
         cbc(6,e,1) = 'O  '
         numapt(e) = 1
         letapt(e) = 'A'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine parse_e_hdr(iel,numapt,letapt,igroup,string)
c
c     Replaces:
c
c        READ(9,'(20X,I4,4X,I3,A1,11x,i5)',ERR=33,END=33)
c    $        IDUM,NUMAPT(IEL),LETAPT(IEL),IGROUP(IEL)
c
c
c     A "typical" format (but this sometimes changes, hence this routine):
c
c           ELEMENT    1 [    1a]    GROUP     0
c

      character*1  letapt
      character*80 string

      character*80 s,t
      character*1   s1(80),t1(80)
      equivalence  (s1,s),(t1,t)


      numapt = 1     ! Std. defaults
      letapt = 'A'
      igroup = 0

      i0 = 1
      call blank (s1,80)
      call chcopy(s1,string,80)

c     4 items: iel, numapt, letapt, igroup

      i1 = nindx1(s1(i0),' ',1)+i0-1  ! 1st non-blank location
      i2 =  indx1(s1(i1),' ',1)+i1-1  ! 2nd blank location
      i1 = nindx1(s1(i2),' ',1)+i2-1  ! 2nd non-blank location
      i2 =  indx1(s1(i1),' ',1)+i1-1  ! 3rd blank location
      n2 =  i2-i1
      call blank (t1,80)
      call chcopy(t,s1(i1),n2)
      read(t,*) iel

c           ELEMENT    1 [    1a]    GROUP     0

      i1 = indx1(s1(i0),'[',1)+1  ! Locate "["
      i2 = indx1(s1(i0),']',1)-1  ! Locate "]"
      n2 =  i2-i1
      call blank (t1,80)
      call chcopy(t,s1(i1),n2)
      read(t,*) numapt
      letapt = s1(i2)

      len = ltrunc(s,80)         ! Get last integer in the string
      do i1=len,1,-1
         if (s1(i1).eq.' ') goto 10
      enddo
   10 continue

      call blank (t1,80)
      call chcopy(t,s1(i1),n2)
      read(t,*) igroup

      return
      end
c-----------------------------------------------------------------------
      subroutine set_main_menu

c***  BIG DECISION POINT  ****

      include 'basics.inc'

      nchoic = 0

      if (if3d) then

         if (ifceil) then !  Restrict choices for person on CEILING

            ITEM(1)='OFF CEILING'
            ITEM(2)='MODIFY ELEMENT'
            ITEM(3)='CURVE SIDES'
            nchoic=3

         else

            nchoic = nchoic+1
            ITEM(nchoic)        =   'END    ELEMENTS'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'MODIFY ELEMENT'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'GLOBAL REFINE'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'CURVE SIDES'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'DELETE ELEMENT'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'ZOOM'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'SET GRID'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'DEFINE OBJECT'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'REDRAW MESH'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'IMPORT MESH'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'IMPORT vtx MESH'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'REFLECT MESH '
            nchoic = nchoic+1
            ITEM(nchoic)        =   'UP   LEVEL'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'DOWN LEVEL'
            nchoic = nchoic+1
            ITEM(nchoic)        =   'CEILING'
         endif

      else  !  2-D

         nchoic = nchoic+1
         ITEM(nchoic)       =       'END    ELEMENTS'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'MODIFY ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'GLOBAL REFINE'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'CURVE SIDES'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'DELETE ELEMENT'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'ZOOM'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'SET GRID'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'DEFINE OBJECT'
         nchoic = nchoic+1
c        ITEM(nchoic)       =       'Edit Mesh'
c        nchoic = nchoic+1
         ITEM(nchoic)       =       'REDRAW MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'IMPORT MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'IMPORT VTK MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'IMPORT vtx MESH'
         nchoic = nchoic+1
         ITEM(nchoic)       =       'REFLECT MESH '
      ENDIF
     
      return   ! End of menu-driven query
      end
c-----------------------------------------------------------------------
      subroutine cell_cell_connectivity

      include 'basics.inc'


      parameter(lpts=8*nelm)
      common /arrayi/ i_n(lpts)   , j_n(4*lpts)
     $              , j_o(4*lpts)
     $              , ee(lpts)   , cell(lpts)
      integer cell,ee

      common /arrayr/ wk (lpts),dx(4*lpts)

      integer clist_i(lpts),clist_j(4*lpts)
      equivalence (clist_i,i_n)
      equivalence (clist_j,j_n)

      integer e,f,fe

      call chker

      nv = 2**ndim

      call set_d2a(dx,x,y,z,nv,nel,ndim)

      call makecell  (cell,nv,nel,dx,ndim,wk,i_n,j_n,j_o)
      call chker2    (cell,nv,nel)


      nf  = 2*ndim
      nv  = 2**ndim
      nic = lpts
      njc = lpts
      itype = 1   ! return structured cell pointer 
      call cell2v(i0,i1,clist_i,nic,clist_j,njc,cell,nv,nel,itype,j_o) 
c     fe=1/(nf-nv)
c     stop

      call find_ee(ee,nf,nel,clist_i,clist_j,cell,nv,ndim)

      do e=1,nel   ! This needs to be fixed for CHT case.
      do f=1,nf
         fe = f + (e-1)*nf
         if (ee(fe).ne.0) then
            je = ee(fe)
            call get_ee_face(jf,je,f,e,ee,nf)
            nfld0=1
            nfld1=nflds
            if(ifconj_merge) then
              nfld0=2
              nfld1=2
            endif
            nflda = nfld0
            if (e.gt.nelv) nflda=2
            do i=nflda,nfld1
               cbc(f,e,i)   = 'E  '
               ibc(f,e,i)   = je
               bc (1,f,e,i) = je
               bc (2,f,e,i) = jf
            enddo
         endif
      enddo
      enddo

      if (nelt.gt.nel) then
         call prs('FIX CHT CASE FOR find_ee. ABORT.$')
         call prexit
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine find_ee(ee,nf,nel,clist_i,clist_j,cell,nv,ndim)
c
c     For each cell, use vertex connectivity to establish cell-cell connections.
c
c     A cell is connected if 2^(ndim-1) vertices are shared with another cell.
c
c     We also scan for pathological cases where the number of connections in
c     3D is only 3 and not 4.
c
c     Algorithm is as follows:
c
c     1. For each vertex, collect cell numbers sharing that vertex.
c
c     2. For each cell face collect lists of cells sharing vertices on face,
c        and tally these.
c
c     3. Look for high tally counts (e.g., 3 or 4 in the 3D case).
c

      integer ee(nf,nel)
      integer clist_i(1),clist_j(1),cell(nv,nel)
      integer e,f,e_share,share(200),wk(200)

      integer flist_3d(4,6)  ! use preprocessor numbering for f AND for v
      save    flist_3d
      data    flist_3d / 1,2,5,6 , 2,3,6,7 , 3,4,7,8 
     $                 , 1,4,5,8 , 1,2,3,4 , 5,6,7,8 /

      integer flist_2d(2,4)  ! use preprocessor numbering for f AND for v
      save    flist_2d
      data    flist_2d / 1,2 , 2,3 , 3,4 , 4,1 /


      nface = 2*ndim
      nvf   = 2**(ndim-1)

      n  = nv*nel
      i0 = iglmin(cell,n)
      i1 = iglmax(cell,n)

      nfail = 0

      do e=1,nel
      do f=1,nface

         ee(f,e) = 0
         n_shared = 0
         do ii=1,nvf      ! Get tally for this face

            iv = flist_3d(ii,f)
            if (ndim.eq.2) iv = flist_2d(ii,f)

            i  = cell(iv,e)
            j0 = clist_i(i)    ! CSR pointer
            j1 = clist_i(i+1)-1

            do jj=j0,j1
               je=clist_j(jj) ! je is an element shared by element e on face f
               if (je.ne.e) then
                  n_shared = n_shared+1
                  share(n_shared) = je
               endif
            enddo
         enddo

         call isortg(share,wk,n_shared)

         ne_3     = 0
         ne_share = 0
         icount   = 0
         ilast    = 0


         do i=1,n_shared
            if (share(i).eq.ilast) then
               icount = icount+1
               if (icount.eq.3) then  ! Consistency check for 3D
                  ne_3 = ne_3+1
               elseif (icount.eq.nvf) then
                  e_share = share(i)
                  ne_share = ne_share+1
               endif
            else
               ilast  = share(i)
               icount = 1
            endif
         enddo
         if (ne_share.eq.1) then ! OK
            if (ndim.eq.2.or.ne_3.eq.1) then
               ee(f,e) = e_share
            elseif (ndim.eq.3.and.ne_3.gt.1) then
               nfail = nfail+1
               write(6,3) 'A',e,f,ne_3,n_shared,ne_share
    3          format('Consistency failure ',a1,i9,i3,3i8)
            endif
         elseif (ne_share.gt.1) then
               nfail = nfail+1
               write(6,3) 'B',e,f,ne_3,n_shared,ne_share
         elseif (ne_share.lt.ne_3) then
               nfail = nfail+1
               write(6,3) 'C',e,f,ne_3,n_shared,ne_share
         else
            ee(f,e) = 0
         endif

      enddo
      enddo

      if (nfail.eq.0) return
      write(6,*) 'FAIL in find_ee: nfail=',nfail
      call prexit

      end
c-----------------------------------------------------------------------
      subroutine unique_vertex2(cell,dx,ndim,nel,q,ind,ninseg,ifseg,wk)
 
      integer cell(1),ind(1),ninseg(1)
      logical ifseg(1)
      real dx(0:ndim,1)

      integer e
      real dxt(4),t1(4),t2(4),wk(0:ndim,1)


      nv = 2**ndim
      n = nv*nel

      write(6,*) 'start locglob_lexico:',nv,nel,n,q

      qq = q*q  ! Square of relative tolerance

      do i=1,n
         cell(i)   = i
         ifseg (i) = .false.
      enddo

c     Sort by directions

      lda         = 1+ndim
      nseg        = 1
      ifseg(1)    = .true.
      ninseg(1)   = n

      do ipass=1,ndim   ! Multiple passes eliminates false positives
      do j=1,ndim       ! Sort within each segment
         write(6,*) 'locglob:',j,nseg,n
         i =1
         j1=j+1
         do iseg=1,nseg
            call tuple_sort(dx(0,i),lda,ninseg(iseg),j1,1,ind,dxt) !key=j1
            call iswap_ip  (cell(i),ind,ninseg(iseg)) ! Swap position
            i  =   i + ninseg(iseg)
         enddo
 
c        q=0.0010   ! Smaller is better
c        q=0.2      ! But generous is more forgiving for bad meshes!

         do i=2,n
           if ((dx(j,i)-dx(j,i-1))**2.gt.qq*min(dx(0,i),dx(0,i-1)))
     $        ifseg(i)=.true.
         enddo

         nseg = 0              !  Count up number of different segments
         do i=1,n
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
      enddo
c
c     Assign global node numbers (now sorted lexigraphically!)
c
      ig  = 0
      ic  = 0
      icm = 0
      do i=1,n
        ic = ic+1     ! count number of instances at present ig

        if (ifseg(i)) then
           ig=ig+1
           icm = max(ic,icm)
           ic  = 0
        endif
        ind(cell(i)) = ig
      enddo
      nglb = ig

c     Unshuffle geometry:
      call tuple_swapt_ip(dx,lda,n,cell,t1,t2)

c     Reassign cell to hold global index numbering

      call icopy     (cell,ind,n)
      write(6,*) 'Performing unique_vertex2 self_chk',n
      call self_chk  (cell,nv,nel,0)       ! check for not self-ptg.


      write(6,6) nseg,nglb,n,icm
    6 format(' done locglob_lexico:',4i10)


      return
      end
c-----------------------------------------------------------------------
      subroutine makecell(cell,nv,nel,dx,ndim,wk,i_n,j_n,j_o)

      integer      cell(1),i_n(1),j_n(1),j_o(1)
      real         dx  (1)
      real         wk  (1)
      integer e,f


      logical ifbinary,ifbswap
      integer buf(30)

      integer eface(6)  ! return Nekton preprocessor face ordering
      save    eface
      data    eface / 4 , 2 , 1 , 3 , 5 , 6 /
         
      q = 0.02
      q = 0.002

c     Compress vertices based on coordinates
      call unique_vertex2(cell,dx,ndim,nel,q,i_n,j_n,j_o,wk)

      nv   = 2**ndim
      n    = nel*nv

c     write(6,*) 'THIS is NEL,NV:',nel,nv,n

      call iranku    (cell,j_n,n)
      write(6,*) 'Performing makecell self_chk',n
      call self_chk  (cell,nv,nel,32)       ! check for not self-ptg.

      return
      end
c-----------------------------------------------------------------------
      subroutine set_d2a(dx,x,y,z,nv,nel,ndim)
      real dx(0:ndim,nv,nel)
      real x(8,nel),y(8,nel),z(8,nel) ! Hardcoded to 8 in prenek :(

      integer e

      if (ndim.eq.3) then
         do e=1,nel
         do i=1,nv
            dx(1,i,e) = x(i,e)
            dx(2,i,e) = y(i,e)
            dx(3,i,e) = z(i,e)
         enddo
         enddo
      else
         do e=1,nel
         do i=1,nv
            dx(1,i,e) = x(i,e)
            dx(2,i,e) = y(i,e)
         enddo
         enddo
      endif

      call set_d2(dx,nv,nel,ndim)  ! Fill up the "0" location

      return
      end
c-----------------------------------------------------------------------
      subroutine set_d2(dx,nv,nel,ndim)
      real dx(0:ndim,nv,nel)

      integer neigh(3,8),e
      save    neigh
      data    neigh / 2,4,5 , 1,3,6 , 2,4,7 , 1,3,8   ! symm. ordering
     $              , 1,6,8 , 2,5,7 , 3,6,8 , 4,5,7 /

c     data    neigh / 2,3,5 , 1,4,6 , 1,4,7 , 2,3,8   ! symm. ordering
c    $              , 1,6,7 , 2,5,8 , 3,5,8 , 4,6,7 / ! used for genmap?

      b = 1.e22

      do e = 1,nel
      do i = 1,nv

         dx(0,i,e) = b

         do k = 1,ndim

            n   = neigh(k,i)
            d2l = 0.

            do j=1,ndim
               d2l = d2l + (dx(j,n,e)-dx(j,i,e))**2
            enddo

            dx(0,i,e) = min (dx(0,i,e),d2l)
c           write(6,9) i,e,' dx',(dx(j,i,e),j=0,ndim)
c  9        format(2i6,a3,1p4e11.3)

         enddo

      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine iranku(a,p,n)
c
c     Return r = rank of a() in place of a(), where r =< n
c     so that a() is the number of unique pts sorted such that
c     it is in the original spot
c     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
c
c
      integer a(1),p(1)
      integer rank

      call isortg(a,p,n)

      rank = 1
      last = a(1)
      a(1) = 1

      do k=2,n
         next = a(k)
         if (next.ne.last) rank=rank+1
         a(k) = rank
         last = next
      enddo

      call iswapt_ip(a,p,n)  ! restore a() to orginal location

      return
      end
c-----------------------------------------------------------------------
      subroutine tuple_swapt_ip(x,m,n,p,t1,t2)
      real x(m,n),t1(m),t2(m)
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            call copy(t1,x(1,loop_start),m)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt(1)
               elseif (next.eq.loop_start) then
                  call copy(x(1,next),t1,m)
                  p(next) = -p(next)
                  goto 10
               else
                  call copy(t2,x(1,next),m)
                  call copy(x(1,next),t1,m)
                  call copy(t1,t2       ,m)
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine iswapt_ip(x,p,n)
      integer x(1),t1,t2
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt(1)
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine self_chk(cell,nv,nel,flag)     ! check for not self-ptg.
      integer cell(nv,nel),flag
      integer e


      n10 = min(nel,10)

      do e=1,nel
      do i=1,nv
         do j=i+1,nv
            if (cell(i,e).eq.cell(j,e)) then

               write(6,*)
               call outmati(cell(1,e),2,4,'SELF!!',e,flag)

               write(6,*)
               write(6,*) 'ABORT: SELF-CHK ',i,j,e,flag
               write(6,*) 'Try to tighten the mesh tolerance!' 

               call outmatti  (cell,nv,n10,'slfchk',nel,flag)

               call prexit
               call exitt(flag)

            endif
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outmatti(u,m,n,name6,nid,ic)
      integer u(m,n)
      character*6 name6
      character*1 adum
c
c     Print out copies of a global matrix, transpose
c
         write(6,1) nid,m,n,name6,nid,ic
   1     format(3i6,'  Matrix T:',2x,a6,2i9)

         m8 = min(m,10)
         do j=1,n
            write(6,2) j,name6,(u(i,j),i=1,m8)
         enddo
   2     format(i8,1x,a6,10i9)

      return
      end
c-----------------------------------------------------------------------
      subroutine outmati(u,m,n,name6,nid,ic)
      integer u(m,n)
      character*6 name6
      character*1 adum
c
c     return
c
c     Print out copies of a global matrix
c
      n200 = min(n,200)
      if (m.gt.1) then
         write(6,1) nid,m,n,name6
   1     format(3i6,'  Matrix:',2x,a6)
         do i=1,m
            write(6,2) i,name6,(u(i,j),j=1,n200)
         enddo
   2     format(i8,1x,a6,20(10i9,/,10x))
      else
         write(6,3) nid,n,name6,(u(1,j),j=1,n200)
   3     format(2i8,1x,a6,20(10i9,/,10x))
      endif
      if (ic.eq.0) then
         write(6,*) 'cont: ',name6,nid,'  ??'
c        read (5,*) adum
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine outmatx (u,m,n,name6,nid)
      real u(m,n)
      character*6 name6
c
c     return
c
c     Print out copies of a global matrix
c
      write(6,1) nid,m,n,name6
      n8 = min(8,n)
   1  format(3i6,'  Matrix:',2x,a6)
      do i=1,m
         write(6,2) i,name6,(u(i,j),j=1,n8)
      enddo
   2  format(i3,1x,a6,1p8e12.4)
c  2  format(i3,1x,a6,20(1p8e12.4,/,10x))

      return
      end
c-----------------------------------------------------------------------
      subroutine cell2v(i0,i1,ic,nic,jc,njc,cell,nv,ncell,type,wk)
c
c     Generate a list of cells that point to each numbered
c     vertex in [i0,i1].  ic/jc is the csr-formated list.
c
c     Cells (aka elements) have nv local verties.
c
      integer ic(1),jc(1),cell(nv,ncell),wk(1),type
c
      i0 = iglmin(cell,nv*ncell)
      i1 = iglmax(cell,nv*ncell)

      if (i1-i0.gt.nic) then
         write(6,1) i0,i1,nic,nv,ncell,njc
    1    format(' ERROR: nic too small in cell2v:',6i10)
         i0 = 0 ! error return code
         i1 = 0 ! error return code
         call prexit
      endif

      call cell2v1(ic,i0,i1,jc,njc,cell,nv,ncell,type,wk)

      return
      end
c-----------------------------------------------------------------------
      subroutine cell2v1(ic,i0,i1,jc,njc,cell,nv,ncell,type,wk)
c
c     Generate a list of cells that point to each numbered
c     vertex in [i0,i1].  ic/jc is the csr-formated list.
c
c     Cells (aka elements) have nv local verties.
c
      integer ic(i0:i1),jc(1),cell(nv,ncell),wk(i0:i1),type

c     Step 1:  count number of cells per vertex
c  the number of cells each vertex 'touches'

      nic = i1-i0+1
      call izero(wk,nic)
      do j=1,nv*ncell
         i = cell(j,1)
         wk(i) = wk(i) + 1
      enddo
c  wk() is an array refering to each vertex, and holds the num of cells
c     Step 2:  Generate csr pointer

      ic(1) = 1
      do i=i0,i1
         ic(i+1) = ic(i) + wk(i)
         wk(i)   = 0
      enddo

c     Step 3:  fill jc()

      if (type.eq.1) then ! return cell
         do k=1,ncell
         do j=1,nv
            i = cell(j,k)
            jc(ic(i)+wk(i)) = k  ! This is the cell number
            wk(i) = wk(i) + 1    ! This is number filled for ith vertex
         enddo
         enddo
         
      else                ! return structured pointer
         do j=1,nv*ncell
            i = cell(j,1)
            jc(ic(i)+wk(i)) = j  ! This is the structured pointer
            wk(i) = wk(i) + 1    ! This is number filled for ith vertex
         enddo
      endif

c     Diagnostics

      write(6,*)
c     do i=i0,i1
c        j0 = ic(i)
c        j1 = ic(i+1)-1
c        write(64,1) i,(jc(j),j=j0,j1)
c  1     format(i8,' c2v:',20i5)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine isortg(a,ind,n)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(1),ind(1)
      integer aa
C
      dO 10 j=1,n
         ind(j)=j
   10 continue
C
      if (n.le.1) return
      L=n/2+1
      ir=n
  100 continue
         if (l.gt.1) then
            l=l-1
            aa  = a  (l)
            ii  = ind(l)
         else
                 aa =   a(ir)
                 ii = ind(ir)
              a(ir) =   a( 1)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
                 a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if ( a(j).lt.a(j+1) ) j=j+1
            endif
            if (aa.lt.a(j)) then
                 a(i) = a(j)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         goto 200
         endif
           a(i) = aa
         ind(i) = ii
      goto 100
      end
c-----------------------------------------------------------------------
      subroutine get_ee_face(jf,je,f,e,ee,nf)

      integer f,e,ee(nf,1)

      do jf=1,nf
         ie=ee(jf,je)
         if (ie.eq.e) return
      enddo
      write(6,*)
      write(6,*) 'ERROR: inconsistent ee pairing:'
      do i=1,nf
         write(6,1) i,e,f,ee(i,e),je,ee(i,je)
      enddo
    1 format(6i8,' i,e,f,ee(i,e),je,ee(i,je) ERR')

      return
      end
c-----------------------------------------------------------------------
      subroutine chker

      include 'basics.inc'

c     character*1 ans

      return

      call prs('Write or Read xyz pts? (r/w):$')
      call res(ans,1)

      open(unit=65,file='xyz.dat')
      if (ans.eq.'r') then
         ntot = 8*nel
         do i=1,ntot
            read(65,*) je,jv,x(i,1),y(i,1),z(i,1)
         enddo
      else
         do je=1,nel
         do jv=1,8
            write(65,1) je,jv,x(jv,je),y(jv,je),z(jv,je)
    1       format(i9,i3,1p3e14.6)
         enddo
         enddo
      endif
      close(65)

      return
      end
c-----------------------------------------------------------------------
      subroutine chker2(cell,nv,nel)

      integer cell(nv,nel)

      character*1 ans

      return
      call prs('Write or Read cell dat? (r/w):$')
      call res(ans,1)

      open(unit=65,file='cell.dat')
      if (ans.eq.'r') then
         do ie=1,nel
            read(65,*) je,(cell(jv,ie),jv=1,8)
         enddo
      else
         do je=1,nel
            write(65,1) je,(cell(jv,je),jv=1,8)
    1       format(i6,8i9)
         enddo
      endif
      close(65)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_flow_heat(ifflow,ifheat,nlogic,io)
      logical ifflow,ifheat
      character*132 string
      do i=1,nlogic
         read(io,132) string
  132    format(a132)
         call capit(string,132)
         if (indx1(string,'IFFLOW' ,6).gt.0) then
              read(string,*) IFFLOW
         elseif (indx1(string,'IFHEAT' ,6).gt.0) then 
              read(string,*) IFHEAT
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine delete  ! Query to delete elements

      include 'basics.inc'
      logical iftmp

      if (nel.eq.0) call prs('ERROR: No elements to delete$')
      if (nel.eq.0) return

C     Find out which element to delete

      if (.not.ifgraf) then
         ielmin = 1
         call drawel(-ielmin)
         call drawis(-ielmin)
         call delel(ielmin)
         return
      endif

      call prs('Enter (w/mouse) 2 points in element to be deleted,$')
      call prs('or, 2 points framing a box containing elements.$')
      call prs('Enter in menu area to abort DELETE ELEMENT op.$')

      iftmp =ifgrid
      ifgrid=.false.

      call prs('Enter 1st point:$')
      call mouse(xmouse,ymouse,button)
      if (xmouse.gt.xphy(1.0)) then ! look for a keypad input
         ifgrid=iftmp 
         call beep
         call prs('Abort DELETE? (y/l=list):$')
         call res(ans,1)
         if (ans.eq.'y'.or.ans.eq.'Y') return
         if (ans.eq.'l'.or.ans.eq.'L') call list_delete
         if (ans.eq.'s'.or.ans.eq.'S') call special_delete_cyls
         return
      else
         call prs('Enter 2nd point:$')
         call mouse(xmous2,ymous2,button)
         if (xmous2.gt.xphy(1.0)) then
            ifgrid=iftmp 
            call beep
            call prs('Abort DELETE? (y/l=list):$')
            call res(ans,1)
            if (ans.eq.'y'.or.ans.eq.'Y') return
            if (ans.eq.'l'.or.ans.eq.'L') call list_delete
            if (ans.eq.'s'.or.ans.eq.'S') call special_delete_cyls
            return
         endif
      endif
C     
C     We successfully input 2 points in the build area
C     1) count number of centroids in the bounding box
C     2) if ncntrd=0, find nearest element and anihilate it.
C     
      xmax=max(xmouse,xmous2) ! Box
      xmin=min(xmouse,xmous2)
      ymax=max(ymouse,ymous2)
      ymin=min(ymouse,ymous2)
C     
C     Check box to see if it contains any elements, and delete them.
C     
C     This is a gross N^2 algorithm, but it beats entering them all
C     by hand....
C     
      numdel=0
 124  continue
      nelt=nel
      DO 125 JEL=1,NELT
         if (jel.le.nel.and.xmin.le.xcen(jel) .and.xcen(jel).le.xmax
     $       .and.ymin.le.ycen(jel) .and.ycen(jel).le.ymax ) then
            if (numdel.eq.0) call drwbox(xmin,ymin,xmax,ymax,1)
            call delel(jel)
            numdel=numdel+1
            goto 124
         endif
 125  continue
      if (numdel.gt.0) then !     redraw the mesh
         call refresh
         call drmenu('NOCOVER')
         call drgrid
         do 128 iel=1,nel
            call drawel(iel)
 128     continue
         ifgrid=iftmp 
         return
      else ! Look for closest element (standard delete element option)
         xmouse=(xmin+xmax)/2.0
         ymouse=(ymin+ymax)/2.0
         rmin=1.0e20
         do 130 jel=1,nel
            if (.not.if3d .or. numapt(jel).eq.ilevel) then
               r2= (xcen(jel)-xmouse)**2 + (ycen(jel)-ymouse)**2
               if(r2.lt.rmin) then
                  rmin  = r2
                  ielmin= jel
               endif
            endif
 130     continue
      endif

      if (if3d) then ! check if it's OK to delete
         ieq=0
         igt=0
         do 140 i=1,nel
            if (numapt(i).eq.ilevel  ) ieq=ieq+1
            if (numapt(i).eq.ilevel+1) igt=igt+1
 140     continue

         if (ieq.eq.1) then ! Last element on floor
            if (igt.gt.0) then
               CALL PRS('**ERROR** You cannot delete all the$')
               CALL PRS('elements on a level when there are $')
               CALL PRS('elements above$')
               ifgrid=iftmp 
               return
            elseif (ilevel.gt.1) then
               CALL PRS(' ** EMPTY LEVEL ** Please ADD ELEMENT$')
               CALL PRS('Go DOWN LEVEL to continue$')
               nlevel=nlevel-1
            endif
         endif
      endif
      call prsi('Deleting Element $',ielmin) ! Now, black out old element
      call drawel(-ielmin)
      call drawis(-ielmin)
      call delel(ielmin)

c     call sortel    ! Sort elements for isometric
c     do 150 i=1,nel ! Redraw all the isometric elements.
C        call drawis(isrt(i))
c150  continue

      if(.not.if3d .and.ielmin.le.nel) call drawel(ielmin)
      ifgrid=iftmp 

      return
      end
c-----------------------------------------------------------------------
      subroutine set_igroup(nelv_i,nelt_i)
      include 'basics.inc'

      ncond=nelt_i - nelv_i ! recalculate ncond
      do i=nelv_i+1,nelt_i
         igroup(i) = 1
      enddo

      return
      end
c-----------------------------------------------------------------------
