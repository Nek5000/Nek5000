c-----------------------------------------------------------------------
      subroutine getfld(newdump,jerr,if_call_setw,if_call_coef)
C----------------------------------------------------------------------
C
C     Modified from routine RESTART, 5-2-24 pff
C     (1) Open restart file(s)
C     (2) Check previous spatial discretization 
C     (3) Map (K1,N1) => (K2,N2) if necessary
C
C----------------------------------------------------------------------
      INCLUDE 'basics.inc'
      INCLUDE 'basicsp.inc'
c
      common /regularize/ regdir
      character*1  regdir
c
      PARAMETER (LXR=NXM)
      PARAMETER (LYR=NYM)
      PARAMETER (LZR=NZM)
      PARAMETER (LXYZR=LXR*LYR*LZR)
C     note, this usage of CTMP1 will be less than elsewhere if nelt ~> 9.
      PARAMETER (LPSC9=MAXFLD+3)
      COMMON /CTMP1/ XDUMP(LXYZR,LPSC9)
      CHARACTER*4    CDUMP(LXYZR,LPSC9)
      EQUIVALENCE   (CDUMP,XDUMP)
      CHARACTER*2  EXCODER(15)
      CHARACTER*80 FNAME
      CHARACTER*1  FNAME1(80)
      EQUIVALENCE (FNAME1,FNAME)
c
      logical if_call_setw
      logical if_call_coef
c
      CHARACTER*80 hname
      CHARACTER*1  hname1(80)
      EQUIVALENCE (hname1,hname)
c
      real*4 bytetest
      real*8 bytetest8
c
      character*1  fname_c(80)
      integer      fname_i(20)
      equivalence (fname_i,fname_c)
      character*2 s2
      character*3 s3
      character*4 s4
c
      INTEGER IDUMP
c
      integer icalld,nelhdr
      save    icalld,nelhdr
      data    icalld /0/
c
c     For "SETWRK":                  pff 2/28/98
      common /sfldwrk/ odumpw,odump
      integer odumpw,odump
c
      logical if_byte_sw, if_byte_swap_test, if_byte_swap_test8
C
C     Local logical flags to determine whether to copy data or not.
C
      LOGICAL IFGETX,IFGETZ,IFGETU,IFGETW,IFGETP,IFGETT,IFGTPS(MAXFLD)
     $       ,IFGTIM,IFOK,IFFMAT
      INTEGER IPOSX,IPOSZ,IPOSU,IPOSW,IPOSP,IPOST,IPSPS(MAXFLD)
      INTEGER rdsize
C
      jerr = 0
      close(unit=24,err=1)
    1 continue
      close(unit=27,err=11)
   11 continue
c
      if (icalld.eq.0) then
         nelhdr = nel
         odump  = -99
      endif
      icalld = icalld+1
c
c     if (newdump.eq.odump) return
      odump  = newdump
      odumpw = newdump
      newdumps=newdump
      NID=0
      nelv=nel
      nelt=nel
      nelgt=nel
      IFOK=.FALSE.
      IFFMAT=.TRUE.
      IF (PARAM(66).EQ.1.0.or.param(66).eq.3) IFFMAT=.FALSE.
      IF (PARAM(66).EQ.4.0.or.param(66).eq.5) IFFMAT=.FALSE.
c
      len = ltrunc(session_name,80)
      call blank(fname,80)
      call chcopy(fname,session_name,len)
      call chcopy(fname1(len+1),'.fld',4)
c
      call blank(hname,80)
      call chcopy(hname,session_name,len)
      call chcopy(hname1(len+1),'.fhd',4)
c
c     Append numerical suffix
c
      len = ltrunc(fname,80)
      if (newdump.ne.0) then
         if (newdump.le.99) then
            write(s2,92) newdump
   92       format(i2.2)
            call chcopy(fname1(len+1),s2,2)
            call chcopy(hname1(len+1),s2,2)
         elseif (newdump.le.999) then
            write(s3,93) newdump
   93       format(i3.3)
            call chcopy(fname1(len+1),s3,3)
            call chcopy(hname1(len+1),s3,3)
         else
            write(s4,94) newdump
   94       format(i4.4)
            call chcopy(fname1(len+1),s4,4)
            call chcopy(hname1(len+1),s4,4)
         endif
      endif
      len = ltrunc(fname,80)
c
c     prepare file name for "byte_open" call 
c
      call izero(fname_i,10)
      call chcopy(fname_c,fname,len)


      IFGETX=.TRUE.
      IFGETZ=.TRUE.
      IFGETU=.TRUE.
      IFGETW=.TRUE.
      IFGETP=.TRUE.
      IFGETT=.TRUE.
      DO 8 I=1,MAXFLD-4
         IFGTPS(I)=.TRUE.
    8 CONTINUE

C     Standard Case:   Requested inputs are in the file.

      write(6,*) 
      write(6,*) 'this is iffmat:',iffmat, param(66)

      call chcopy(fld_name,fname,80)
      write(6,80) session_name
      write(6,80) fname
   80 format(a80)
      write(6,*) 
c
      IF (IFFMAT) THEN
         write(6,*) 'opening file:',len,' ',(fname1(k),k=1,len)
         write(6,*) 'opening file:',len,' ',fname
         open (UNIT=24,FILE=FNAME,ERR=9,status='OLD')
         write(6,*) 'opened  file:',len,' ',fname
      elseif (param(66).eq.2) then
         open (UNIT=24,FILE=FNAME
     $     ,FORM='UNFORMATTED',ERR=9,status='OLD')
      elseif (param(66).eq.3) then
         write(6,*) 'opening file:',len,' ',(fname1(k),k=1,len)
         write(6,*) 'opening file:',len,' ',(hname1(k),k=1,len)
         open (UNIT=27,FILE=hname,ERR=9,status='OLD')
         open (UNIT=24,FILE=FNAME,ERR=9,status='OLD')
         close(unit=24)
         call byte_open(fname_c)
      elseif (param(66).eq.4.or.param(66).eq.5) then
         write(6,*) 'opening file:',len,' ',(fname1(k),k=1,len)
         open (UNIT=24,FILE=FNAME,ERR=9,status='OLD')
         close(unit=24)
         call byte_open(fname_c)
      ENDIF
      GOTO 10
c
c     Check file
c
    9 CONTINUE
c
c        Error condition
c
         write(6,*) 'this is param(66):',param(66)
         len = ltrunc(fname,80)
         fname1(len+1) = '$'
         call prs('WARNING:  Could not open file.$')
         call prs(fname)
         ierr = 1
         return
   10 continue
c
C
      rdsize = 4
      ndumps_read = newdumps
      if (param(65).ne.0) ndumps_read = 1
      ndumps_read = 1                        !   11/04/02
c
      DO 1000 JDUMP=1,ndumps_read
C
         IF (NID.EQ.0.OR.JDUMP.EQ.1) THEN
           call blank(hname,80)
           IF (IFFMAT) THEN
            READ(24,80,ERR=1500,END=1500) hname
           ELSEIF(param(66).eq.1) then
            READ(24,ERR=1500,END=1500) hname
           ELSEIF(param(66).eq.3) then
            READ(27,80,ERR=1500,END=1500) hname
           ELSEIF(param(66).eq.4.or.param(66).eq.5) then
            call byte_read(hname,20)
           endif
           write(6,*) ' got hname',nelgt,neltr,nel,nelhdr
           write(6,80) hname
c
           nxr=nx
           nyr=ny
           nzr=nz
           neltr=nel
           if (indx1(hname,'*',1).ne.0) then
c              hack to get around overflow in file header...
               read(hname,'(16x,1X,G13.4,5x,1X,10A2)',ERR=1501,END=1501)
     $         RSTIME,(EXCODER(I),I=1,10)
           else
            if (nelhdr.le.9999) then
            READ(hname,'(4I4,1X,G13.4,I5,1X,10A2,i2)',ERR=1502,END=1502)
     $        neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10),RDSIZE
            else
              read
     $        (hname,'(i10,3i4,1P1e18.9,i9,1x,10a2)',ERR=1503,END=1503)
     $        neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
              rdsize = 4
            endif
           endif
           IF (PARAM(68).EQ.2.0) 
     $      WRITE(29,'(4I4,1X,G13.4,I5,1X,10A2)')
     $      neltr,nxr,nyr,nzr,rstime,istepr,(excoder(i),i=1,10)
           write(6,*) ' done read from hname',nxr,nyr,nzr
c
           time=rstime
           if (rdsize.eq.0) rdsize=4
           write(6,*) 'this is rdsize:',rdsize
         ENDIF
c
c        Reset spatial dimension according to contents of .fld file
c
         param(20) = nxr
         nx = nxr
         ny = nxr
         nz = 1
         if (nx.gt.nxm.or.ny.gt.nym.or.nz.gt.nzm) then
            write(6,*) 'ERROR - nxm,nym,nzm too small.',nxm,nym,nzm
            write(6,*) 'Increase in basics.inc and recompile:',nx,ny,nz
            call exitt
         endif
c
         if (ndim.eq.3) nz = nxr
         neltm=maxpts/(nx*ny*nz)
         neltr=min(neltr,neltm)
         nel  =min(neltr,nel)
         write(6,*) neltm,maxpts,nxr,nel,' NELTM '
         IF (JDUMP.EQ.1) THEN
C
C           Find out the position of the data via the header info.
C
C           Assumptions about position of data:
C
C           .If X is found then Y (Z) are implied
C           .If U is found then V (W) are implied
C           .Column data is space or , delimited.
C
C
C           Figure out position of data in file "IFILE"
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
   40       CONTINUE
C
c           if (idump.eq.1) write(6,41) (excoder(j),j=1,10),idump
            write(6,41) (excoder(j),j=1,10),idump
   41       format(' excoder:',2x,10a2,i4)
            DO 50 I=1,10
               IF (EXCODER(I).EQ.'X') THEN
                  NOUTS=NOUTS + 1
                  IPOSX=NOUTS
                  IFGETX = .TRUE.
               elseif (EXCODER(I).EQ.'Y') THEN
                  NOUTS=NOUTS+1
                  IPOSY=NOUTS
                  IFGETX = .TRUE.
               elseif (EXCODER(I).EQ.'Z') THEN
                  NOUTS=NOUTS + 1
                  IPOSZ=NOUTS
                  IFGETZ = .TRUE.
               elseif (EXCODER(I).EQ.'U') THEN
                  NOUTS=NOUTS + 1
                  IPOSU=NOUTS
                  NOUTS=NOUTS+1
                  IPOSV=NOUTS
                  IFGETU = .TRUE.
                  IF (IF3D) THEN
                     NOUTS=NOUTS + 1
                     IPOSW=NOUTS
                     IFGETW = .TRUE.
                  ENDIF
                  IFFLOW=.TRUE.
               elseif (EXCODER(I).EQ.'P') THEN
                  NOUTS=NOUTS + 1
                  IPOSP=NOUTS
                  IFGETP=.TRUE.
               elseif (EXCODER(I).EQ.'T') THEN
                  NOUTS=NOUTS + 1
                  IPOST=NOUTS
                  IFGETT=.TRUE.
                  IFHEAT=.TRUE.
               elseif (EXCODER(I).EQ.'V') THEN
                  IFFLOW=.TRUE.
               elseif (EXCODER(I).NE.' ') THEN
                  write(6,*) 'Excoder:', i,excoder(i),nouts
                  read(excoder(i),111) ii
  111             format(i1)
c
                  NOUTS=NOUTS + 1
                  IPSPS(ii) = nouts
                  IFGTPS(ii)=.true.
               ENDIF
         write(6,*) 
     $  'IFGNGM,IFGETX bbb',ifgngm,ifgetx,i,excoder(i),iposx,iposy
   50       CONTINUE
C
            LNAME=LTRUNC(FNAME,80)
            IF (NID.EQ.0) WRITE(6,61) (FNAME1(I),I=1,LNAME)
            IF (NID.EQ.0) WRITE(6,62) IPOSU,IPOSV,IPOSW,IPOST,NOUTS
   61       FORMAT(/,2X,'Reading from file ',80A1)
   62       FORMAT(2X,'Columns for data, U,V,W,T,N: ',5I4)
            write(6,*) 'IFGNGM,IFGETX ccc',ifgngm,ifgetx,iposx,iposy
C
C           Make sure the requested data is present in this file....
            IF (IPOSX.EQ.0) IFGETX=.FALSE.
            IF (IPOSY.EQ.0) IFGETX=.FALSE.
            IF (IPOSZ.EQ.0) IFGETZ=.FALSE.
            IF (IPOSU.EQ.0) IFGETU=.FALSE.
            IF (IPOSW.EQ.0) IFGETW=.FALSE.
            IF (IPOSP.EQ.0) IFGETP=.FALSE.
            IF (IPOST.EQ.0) IFGETT=.FALSE.
                            IFGTIM=.TRUE.
            DO 65 I=2,NPSCAL
               IF (IPSPS(I).EQ.0) IFGTPS(I)=.FALSE.
   65       CONTINUE
            write(6,*) 'IFGNGM,IFGETX ddd',ifgngm,ifgetx,iposx,iposy
C
C           End of restart file header evaluation.
C
         ENDIF
C
C        Read the error estimators
C
         if_byte_sw = .false.
         IF(NID.EQ.0)THEN
            IF (IFFMAT) THEN
             read(24,'(6g11.4)',end=1504)(cerror(1,i),i=1,neltr)
             IF (PARAM(68).EQ.2.0) 
     $       write(29,'(6g11.4)')(cerror(1,i),i=1,neltr)
            ELSEIF(param(66).eq.1) then
             read(24,end=1505)(cerror(1,i),i=1,neltr)
            ELSEIF(rdsize.eq.4) then
c            read test pattern to determine if byte-swapping is req'd
             call byte_read(bytetest,1)
             if_byte_sw = if_byte_swap_test(bytetest)
            ELSEIF(rdsize.eq.8) then
             call byte_read(bytetest8,2)
             if_byte_sw = if_byte_swap_test8(bytetest8)
            ENDIF
         ENDIF
            write(6,*) 'IFGNGM,IFGETX 1dd',ifgngm,ifgetx,iposx,iposy
C
C        Read the current dump, double buffer so that we can
C        fit the data on a distributed memory machine,
C        and so we won't have to read the restart file twice
C        in case of an incomplete data file.
C
         NXYZR = NXR*NYR*NZR
C
C        Read the data
C
         write(6,*) 'reading data',idump,nouts,sdump(1,1),rdsize

         i100 = 10
         if (neltr/i100.gt.100) i100 = i100*10
         if (neltr/i100.gt.100) i100 = i100*10
         if (neltr/i100.gt.100) i100 = i100*10
         if (neltr/i100.gt.100) i100 = i100*10
         if (neltr/i100.gt.100) i100 = i100*10
         if (neltr/i100.gt.100) i100 = i100*10

         do 200 ieg=1,neltr
            IF (NID.EQ.0) THEN
c
              IF (MOD(IEG,i100).EQ.0 .or. ieg.eq.1 ) 
     $        WRITE(6,*) 'Reading',IEG,nouts,nxyzr,neltr
c
              IF (PARAM(66).EQ.0.0) THEN
                 READ(24,*,ERR=1506,END=1506)
     $           ((XDUMP(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
              ELSEIF (PARAM(66).EQ.2.0) THEN
                 READ(24,201,ERR=1507,END=1507)
     $           ((CDUMP(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
  201            FORMAT(20A4)
C
C
                 DO 202 II=1,NOUTS
                    IF (PARAM(68).EQ.2.0) THEN
C                   .... fix up for Delta
                       CALL VRNVERTq(XDUMP(1,II),NXYZR)
                    ELSE
                       CALL VRNVERT(XDUMP(1,II),NXYZR)
                    ENDIF
  202            CONTINUE
C
                 IF (PARAM(68).EQ.2.0) THEN
C                   .... fix up for Delta
                    DO 2021 II=1,NOUTS
 2021               CALL VCNVERT(XDUMP(1,II),NXYZR)
                    WRITE(29,201)
     $             ((XDUMP(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
                 ENDIF
              ELSEIF(param(66).ge.3.and.rdsize.eq.4) then
                 do ii=1,nouts
                    call byte_read(xdump(1,ii),nxyzr)
                 enddo
              ELSEIF(param(66).ge.3.and.rdsize.eq.8) then
                 call prs('8 byte reads not yet supported in getfld')
c                do ii=1,nouts
c                   call byte_read(xdump(1,ii),nxyzr)
c                enddo
              ELSE
                 READ(24,ERR=1508,END=1508)
     $           ((XDUMP(IXYZ,II),II=1,NOUTS),IXYZ=1,NXYZR)
              ENDIF
              IFOK=.TRUE.
            ENDIF
C
C           Notify other processors that we've read the data OK.
C
c           CALL LBCAST(IFOK)
            IF (.NOT.IFOK) GOTO 1600
C
C           MAPDMP maps data from NXR to NX
C           (and sends data to the appropriate processor.)
C
C           The buffer SDUMP is used so that if an incomplete dump
C           file is found (e.g. due to UNIX io buffering!), then
C           the previous read data stored in VX,VY,.., is not corrupted.
C
c     write(6,*) ieg,ifgetx,ifgetz,iposx,iposz
c    $        ,xdump(1,iposx),xdump(1,iposy),xdump(1,iposz)
            if (ieg.le.nelt) then
              IF (IFGETX) CALL MAPDMP
     $        (SDUMP(1,1),XDUMP(1,IPOSX),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETX) CALL MAPDMP
     $        (SDUMP(1,2),XDUMP(1,IPOSY),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETZ) CALL MAPDMP
     $        (SDUMP(1,3),XDUMP(1,IPOSZ),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETU) CALL MAPDMP
     $        (SDUMP(1,4),XDUMP(1,IPOSU),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETU) CALL MAPDMP
     $        (SDUMP(1,5),XDUMP(1,IPOSV),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETW) CALL MAPDMP
     $        (SDUMP(1,6),XDUMP(1,IPOSW),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETP) CALL MAPDMP
     $        (SDUMP(1,7),XDUMP(1,IPOSP),IEG,NXR,NYR,NZR,if_byte_sw)
              IF (IFGETT) CALL MAPDMP
     $        (SDUMP(1,8),XDUMP(1,IPOST),IEG,NXR,NYR,NZR,if_byte_sw)
              do i=1,npscal
                 i8 = i+8
                 IF (ifgtps(i)) CALL MAPDMP
     $              (SDUMP(1,i8),XDUMP(1,ipsps(i))
     $              ,IEG,NXR,NYR,NZR,if_byte_sw)
              enddo
C             passive scalars
              NP4=MIN(NPSCAL,4)
              NP4=NPSCAL-4
              NP4=MIN(NP4,4)
            ENDIF
  200      CONTINUE               
            write(6,*) 'IFGNGM,IFGETX 2dd',ifgngm,ifgetx,iposx,iposy
C
C        Successfully read a complete field, store it.
C
C        General geometry?
         IF (IFGETX) IFGNGM=.TRUE.
         write(6,*) 'IFGNGM,IFGETX aaa',ifgngm,ifgetx
C        passive scalars
         NP4=MIN(4,NPSCAL)
         NP4=NPSCAL-4
         NP4=MIN(NP4,4)
         TIME=RSTIME
         ISTEP = ISTEPR
 1000 CONTINUE
      GOTO 1600
C
C     Else, end of file found - notify other processors.
 1500    write(s,1591) jdump,ieg,'A'
         goto 1599
 1501    write(s,1591) jdump,ieg,'B'
         goto 1599
 1502    write(s,1591) jdump,ieg,'C'
         goto 1599
 1503    write(s,1591) jdump,ieg,'D'
         goto 1599
 1504    write(s,1591) jdump,ieg,'E'
         goto 1599
 1505    write(s,1591) jdump,ieg,'F'
         goto 1599
 1506    write(s,1591) jdump,ieg,'G'
         goto 1599
 1507    write(s,1591) jdump,ieg,'H'
         goto 1599
 1508    write(s,1591) jdump,ieg,'I'
         goto 1599
 1591    format
     $('Unable to read dump',I4,' near element',I8,1x,a1,'$')

 1599    ifok=.false.
         call prs(s)

 1600 CONTINUE
C
      write(6,*) 'IFGNGM,IFGETX 3dd',ifgngm,ifgetx,iposx,iposy
      write(6,*) 'IFGNGM,IFGETX 3de',nel,neltr
      IF (JDUMP.EQ.1.AND.NID.EQ.0) THEN
         WRITE(6,1700) FNAME
         WRITE(6,1701) IEG,IXYZ
         WRITE(6,1702) 
     $         ((XDUMP(JXYZ,II),II=1,NOUTS),JXYZ=IXYZ-1,IXYZ)
 1700    FORMAT(5X,'WARNING:  No data read in for file ',A80)
 1701    FORMAT(5X,'Failed on  element',I4,',  point',I5,'.')
 1702    FORMAT(5X,'Last read dump:',/,5G15.7)
      ELSE
         IDUMP=JDUMP-1
         NDUMPS=IDUMP
         if (newdump.ne.1) idump=newdump
C
         IF (NID.EQ.0) WRITE(6,1800) IDUMP
 1800    FORMAT(2X,'Successfully read data from dump number',I3,'.')
         TIME=RSTIME
         ISTEP = ISTEPR
         IDSTEP(IDUMP)=ISTEP
         DMPTIM(IDUMP)=TIME
      ENDIF
c
      IF (param(66).eq.3.) then
         close(unit=27)
         call byte_close()
      elseif (param(66).eq.4 .or. param(66).eq.5) then
         call byte_close()
      else
         close(unit=24)
      endif
            write(6,*) 'IFGNGM,IFGETX 4dd',ifgngm,ifgetx,iposx,iposy
c
      GOTO 6000
C
C     Can't open file...
 5000 CONTINUE
      IF (NID.EQ.0) WRITE(6,5001) FNAME 
 5001 FORMAT(2X,'   *******   WARNING   *******    '
     $    ,/,2X,'   *******   WARNING   *******    '
     $    ,/,2X,'   Could not open fld file:'
     $    ,/,A80
     $   ,//,2X,'ASSUMING DEFAULT INITIAL CONDITIONS.')
c     call exitt
c
      IF (param(66).eq.3.) then
         close(unit=27)
         call byte_close()
      elseif (param(66).eq.4 .or. param(66).eq.5) then
         call byte_close()
      else
         close(unit=24,err=5009)
      endif
c
 5009 CONTINUE
C
C
C     End of IFILE loop
 6000 CONTINUE
C
c
      if (ifregz) then
         call prs('regularize$')
         write(6,*) 'call mapz: ',regdir
         call mapz(u,regdir)
         call mapz(v,regdir)
         call mapz(w,regdir)
         call mapz(p,regdir)
         call mapz(t,regdir)
      endif
C
            write(6,*) 'IFGNGM,IFGETX 5dd',ifgngm,ifgetx,iposx,iposy
      write(6,66) 
     $  ifavgupt,if_call_setw,if_call_coef,ifgetx,ifgngm,ifregz
   66 format('ifavgupt,call_setw,call_coef,ifgetx,ifgngm',7l4)
c
      if (ifavgupt)                      call avg_uvwpt_regular
            write(6,*) 'IFGNGM,IFGETX 6dd',ifgngm,ifgetx,iposx,iposy
      if (ifgetx)                        call reset_xc
            write(6,*) 'IFGNGM,IFGETX 7dd',ifgngm,ifgetx,iposx,iposy
      if (ifgetx.and.if_call_coef)       call coef
            write(6,*) 'IFGNGM,IFGETX 8dd',ifgngm,ifgetx,iposx,iposy
c     if (ifgngm)                        call coef
            write(6,*) 'IFGNGM,IFGETX 9dd',ifgngm,ifgetx,iposx,iposy
      if (if_call_setw)                  call setwrk(.true.)
            write(6,*) 'IFGNGM,IFGETX 0dd',ifgngm,ifgetx,iposx,iposy
c
c     call prs('CALLING avg_cylinder.$')
c     call avg_cylinder
C
            write(6,*) 'IFGNGM,IFGETX add',ifgngm,ifgetx,iposx,iposy
      return
      end
C
c-----------------------------------------------------------------------
      SUBROUTINE MAPDMP(YD,XD,IEG,NXR,NYR,NZR,if_byte_sw)
C----------------------------------------------------------------------
      INCLUDE 'basics.inc'
      INCLUDE 'basicsp.inc'
C
      PARAMETER (LXYZ1=NXM*NYM*NZM)
      PARAMETER (LXR=NXM)
      PARAMETER (LYR=NYM)
      PARAMETER (LZR=NZM)
      PARAMETER (LXYZR=LXR*LYR*LZR)
c
      logical if_byte_sw
C
      dimension yd(nxr,nyr,nzr,nelm)
      dimension xd(lxyzr)
      NXYZ=NX*NY*NZ
      NXYR=NXR*NYR*NZR
C
C
      CALL COPY(YD(1,1,1,IEG),XD,NXYZ)
      if (if_byte_sw) call byte_reverse(yd(1,1,1,ieg),nxyz)
c
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE VRNVERTq(A,N)
      REAL*4 A(1)
      ioldexp=-2
      A(1)=REVERTq(A(1),ioldexp,izero)
      DO 100 I=2,N
         i1=i-1
         A(I)=REVERTq(A(I),inewexp,izero)
         if (izero.ne.0.or.mod(i1,1000).eq.0) then
            a(i1)=a(i1)*10.0**(inewexp-ioldexp)
         endif
         ioldexp=inewexp
  100 CONTINUE
      return
      END
      FUNCTION REVERTq(T,ioldexp,izero)
C
      REAL*4 REVERT1
      CHARACTER*4 T
      CHARACTER*12 VC
      CHARACTER*1  V1(12)
      EQUIVALENCE (V1,VC)
C
      INTEGER IC
      CHARACTER*1  C(4)
      EQUIVALENCE (IC,C)
C
      CHARACTER*4  WC
      CHARACTER*1  W1(4)
      EQUIVALENCE (W1,WC)
C
      CHARACTER*1 ALPSRT(4,0:63)
      CHARACTER*4 ALPSR4(0:63)
      SAVE        ALPH64,ALPSRT
      INTEGER INDEX(0:63),INTALP(0:63),LOG64
      INTEGER IWk  (0:63)
      SAVE    INDEX,INTALP,LOG64
      EQUIVALENCE (INTALP,ALPSRT)
      EQUIVALENCE (INTALP,ALPSR4)
C
      CHARACTER*1 ALPH64(0:63)
      INTEGER ICALLD
      SAVE    ICALLD
      DATA    ICALLD /0/
      DATA        ALPH64 
     $  /'1','2','3','4','5','6','7','8','9','0'
     $   ,'a','b','c','d','e','f','g','h','i','j'
     $   ,'k','l','m','n','o','p','q','r','s','t'
     $   ,'u','v','w','x','y','z'
     $   ,'A','B','C','D','E','F','G','H','I','J'
     $   ,'K','L','M','N','O','P','Q','R','S','T'
     $   ,'U','V','W','X','Y','Z','+','-'/
C
C     On first call, set up sorted list of Alpha values so
C     that we can use binary chop to evaluate the integer
C     values.
C
      IF (ICALLD.EQ.0) THEN
         ICALLD=1
         CALL BLANK(ALPSRT,256)
         DO 10 I=0,63
            ALPSRT(4,I)=ALPH64(I)
   10    CONTINUE
         CALL ISORT(INTALP,INDEX,64)      ! isort is now a sort, not a rank
c        CALL ISWAP(INTALP,Iwk,INDEX,64)  ! pff 12/2/04
         DO 20 I=0,63
   20    INDEX(I)=INDEX(I)-1
         LOG64=7
      ENDIF
C
C     Copy T to W
C
      WC=T
      CALL BLANK(C,4)
C
C     Begin divide and conquer search
C
      ILO=0
      IHI=63
      C(4)=W1(1)
      DO 100 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IHUN=INDEX(I)
c           write(6,103) ihun,i,ilog,ilo,ihi
  103 format(' ihun',6i9)
            GOTO 101
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  100 CONTINUE
  101 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(2)
      DO 200 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            ITEN=INDEX(I)
c           write(6,203) iten,i,ilog,ilo,ihi
  203 format(' iten',6i9)
            GOTO 201
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  200 CONTINUE
  201 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(3)
      DO 300 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IONE=INDEX(I)
c           write(6,303) ione,i,ilog,ilo,ihi
  303 format(' ione',6i9)
            GOTO 301
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  300 CONTINUE
  301 CONTINUE
C
      ILO=0
      IHI=63
      C(4)=W1(4)
      DO 400 ILOG=0,LOG64
         I=(ILO+IHI)/2
         IF (IC.EQ.INTALP(I)) THEN
            IEXP=INDEX(I)
c           write(6,403) iexp,i,ilog,ilo,ihi
  403 format(' iexp',6i9)
            GOTO 401
         ELSEIF (IC.LT.INTALP(I)) THEN
            IHI=I-1
         ELSE
            ILO=I+1
         ENDIF
  400 CONTINUE
  401 CONTINUE
C
c     MANTIS=4096*IHUN+64*ITEN+IONE
      MANTIS=4096*IHUN+64*ITEN+IONE-131072
      IEXP  =IEXP-31
      IF (MANTIS.LT.0) THEN
         MANTIS=ABS(MANTIS)
         WRITE(VC,1001) MANTIS,IEXP
      ELSE
         WRITE(VC,1002) MANTIS,IEXP
      ENDIF
 1001 FORMAT('-0.',I5.5,'E',I3.2)
 1002 FORMAT(' 0.',I5.5,'E',I3.2)
C
      READ(VC,2001) W
 2001 FORMAT(E12.5)
 2002 FORMAT(2X,A12)
C
      izero=1
      if (w.eq.0.0) izero=0
      ioldexp=iexp
      REVERTq=W
      return
      END
c-----------------------------------------------------------------------
      SUBROUTINE VRNVERT(A,N)
      REAL*4 A(1)
      DO 100 I=1,N
         A(I)=REVERT(A(I))
c        write(6,*) a(i),i
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine copy_byte_sw(y,x,n)
c
      character*1 y(0:3,1),x(0:3,1)
c
      do i=1,n
         do j=0,3
            y(3-j,i) = x(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test8(bytetest)
c
      real*8 bytetest
      real*8 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      if_byte_swap_test8 = .false.
      if (bytetest.ne.test_pattern) if_byte_swap_test8 = .true.
c
      write(6,*) 'Byte swap8:',if_byte_swap_test8,bytetest
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)
c
      real*4 bytetest
      real*4 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      eps          = 0.0001
      if_byte_swap_test = .true.
      if (abs(bytetest-test_pattern).lt.eps) if_byte_swap_test = .false.
c
      write(6,*) 'Byte swap:',if_byte_swap_test,bytetest
      return
      end
c-----------------------------------------------------------------------
      subroutine openfld(iunit,file,SUFFCH,SUFFI,status,iffmat,ierr)
      CHARACTER*(*) FILE,STATUS
      CHARACTER*30 NEWFILE,LFILE
      CHARACTER*80 zfile
      character*1  lfile1(30)
      character*2  s2
      character*3  s3
      equivalence (lfile1,lfile)
      character*4 suffch
      integer     suffi
      logical iffmat
C
C     NULL out local filename character strings
      call blank(lfile  ,30)
      call blank(newfile,30)
C
      LFILE = FILE
      IF (STATUS.EQ.'NEW' .OR. STATUS.EQ.'new') THEN
         NEWFILE=FILE
         do 10 i=30,1,-1
            if(ichar(newfile(I:I)).ne.0.and.newfile(I:I).ne.' ') then
               newfile(I+1:I+1) = '~'
               GO TO 20
            endif
   10    continue
   20    continue
         print*,'renaming ',file,' ',newfile
         MVERR=RENAME(FILE,NEWFILE)
      ENDIF
      if (suffi.ne.0) then
c
c        Append suffix to file name
c
c        len1=ltrunc(lfile,30)+1
         len1=indx1(lfile,suffch,4)+4
         if (suffi.lt.100) then
            write(s2,32) suffi
            call chcopy(lfile1(len1),s2,2)
   32       format(i2.2)
         else
            write(s3,33) suffi
            call chcopy(lfile1(len1),s3,3)
   33       format(i3.3)
         endif
      endif
c
      len=ltrunc(lfile,30)
      call izero(zfile,20)
      call chcopy(zfile,lfile,len)
c
      if (iffmat) then
         write(6,*) 'Trying to open fld file:',zfile,status,iffmat
         write(6,*) 'Trying to open fld file:',lfile,status,iffmat
c        OPEN(UNIT=IUNIT,FILE=zfile,STATUS=STATUS
         OPEN(UNIT=IUNIT,FILE=lfile,STATUS=STATUS
     $       ,FORM='FORMATTED',ERR=1)
         IERR=0
c        IFFMAT=.TRUE.
         write(6,*) 'Successfully opened file: '
     $             ,lfile,' Format:',iffmat
         return
      endif
c
    1 continue
      if (.not.iffmat) then
         write(6,*) 'Trying to open unfmt fld file:',zfile,status
         OPEN(UNIT=IUNIT,FILE=zfile,STATUS=STATUS
     $      ,FORM='UNFORMATTED',ERR=2)
         IERR=0
c        IFFMAT=.FALSE.
         write(6,*) 'Successfully opened file: '
     $             ,lfile,' Format:',iffmat
         return
      endif
c
    2 continue
      IERR=1
c     IFFMAT=.TRUE.
      return
      END
c-----------------------------------------------------------------------
      subroutine reset_xc
      INCLUDE 'basics.inc'
      INCLUDE 'basicsp.inc'
c
      integer ied(8)
      save    ied
      data    ied / 1,2,4,3 , 5,6,8,7 /
c
      do ie=1,nel
         l=0
         if (if3d) then
            do k=1,nz,nz-1
            do j=1,ny,ny-1
            do i=1,nx,nx-1
               l=l+1
               m=i+nx*(j-1) + nx*ny*(k-1) + nx*ny*nz*(ie-1)
               x(ie,ied(l)) = xp(m)
               y(ie,ied(l)) = yp(m)
               z(ie,ied(l)) = zp(m)
            enddo
            enddo
            enddo
         else
            do j=1,ny,ny-1
            do i=1,nx,nx-1
               l=l+1
               m=i+nx*(j-1) + nx*ny*(ie-1)
c
c              The following is to foil the stupid optimizer
c
               if (j.gt.90) write(6,*) 'm:',m,i,j,nx,ny,ie,l
               x(ie,ied(l)) = xp(m)
               y(ie,ied(l)) = yp(m)
            enddo
            enddo
         endif
      enddo
c
      call gencen
C
      return
      end
c-----------------------------------------------------------------------
