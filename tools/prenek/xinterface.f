      PROGRAM MAINF
      CALL MAINC
      STOP
      END
      SUBROUTINE READER(VALUE,IERR)
#     include "basics.inc"
      IERR=0
      REWIND(13)
      WRITE(13,'(A70)')LINE        
      REWIND(13)
      IF(LINE.EQ.' ')THEN
         VALUE=0
         RETURN
      ENDIF
      READ(13,*,ERR=2,END=2)VALUE
      REWIND(13)
      RETURN
2     IERR=1
      RETURN
      END

      SUBROUTINE OPENF(IUNIT,FILE,STATUS,ITYPE,IERR)
      CHARACTER*(*) FILE,STATUS
      CHARACTER*30 NEWFILE,LFILE
      CHARACTER NULL
C
C     NULL out local filename character strings
      NULL = CHAR(0)
      DO 5 I=1,30
         LFILE  (I:I) = NULL
         NEWFILE(I:I) = NULL
   5  CONTINUE
C
      LFILE = FILE
      if (status.eq.'new' .or. status.eq.'NEW') then
         newfile=file
         do 10 i=30,1,-1
            if (ichar(newfile(I:I)).ne.0.and.newfile(I:I).ne.' ') then
               newfile(I+1:I+1) = '~'
               goto 20
            endif
 10      continue
 20      continue
c        print*,'renaming ',file,' ',newfile
c        mverr=rename(file,newfile)
         open(unit=iunit,file=lfile,status='unknown',err=1)
      else
         open(unit=iunit,file=lfile,status=status,err=1)
      endif

        ierr=0
        return

1       ierr=1
        return

        end

        SUBROUTINE FINDFL(SESION,OLDVER)
        RETURN
        END

        SUBROUTINE MOUSE(XMOUSE,YMOUSE,BUTTON)
#     include "basics.inc"
	INTEGER*4 RECORD(3),SIZE,PET
	INTEGER ERRIND,MODE,ESW,XFORM,LDR,RETSIZ,ISTAT,DEVNO
	REAL LOCX,LOCY,EAREA(4)
        COMMON/MOUSEFAC/ WINDOW
C
C                                                Universal
C       BUTTON only recorded in interactive mode, not from journal
C       BUTTON=' ' for typed in input or when reading from journal
C       BUTTON='LEFT','RIGHT', or 'MIDDLE'
C	NEW! IMPROVED!!! Green button means "EXIT" end, etc
C       Return XMOUSE,YMOUSE, BUTTON; backspace and write x,y in journal file
20      CONTINUE
c
C       Look for button input on SUN
97      CONTINUE

        IF(.NOT.IFDEMO)call XWAITM(IBUTTON,xm,ym)
C       IF(     IFDEMO)CALL XWAITM(IBUTTON,XM,YM)
        IF(     IFDEMO)READ(55,*)IBUTTON,XM,YM
c       IF(     IFDEMO)READ(5,*)IBUTTON,XM,YM

        IF(IFLEARN)WRITE(3,*)IBUTTON,XM,YM
        if(ibutton.eq.0) then
           CALL PRS('Still waiting for mouse input$')
           go to 97
	endif	
        IF(IBUTTON.EQ.1)BUTTON='LEFT'
        IF(IBUTTON.EQ.2)BUTTON='MIDDLE'
        IF(IBUTTON.EQ.3)BUTTON='RIGHT'
	window=1.3
        XSCREEN=(XM/WINDOWW*WINDOW) 
        YSCREEN=(YM/(WINDOWH-180.))-0.05
C     OLD VERSION -- coordinate transformation is almost
C                    device dependent!
C        YSCREEN=((YM-12)/1000.0*WINDOW) + 0.015
        IF(XSCREEN .EQ. OLDXSCREEN .AND.
     $     YSCREEN .EQ. OLDYSCREEN .AND. 
     $     CHOICE  .EQ.'ACCEPT CURRENT SWITHCES') THEN
           CALL PRS('Error reading mouse input; please move mouse'//
     $              'slightly and re-enter$')
           GO TO 97
        ENDIF
        OLDXSCREEN = XSCREEN
        OLDYSCREEN = YSCREEN
        IF(XSCREEN .LT. 0.0 .OR. XSCREEN .GT. 1.3  .OR.
     $     YSCREEN .LT. 0.0 .OR. YSCREEN .GT. 1.0) THEN
           CALL PRS('Error reading mouse input; please re-enter$')
           GO TO 97
        ENDIF
        IF(XSCREEN .EQ. 0.0 .AND. YSCREEN .EQ. 1.5e-2) THEN
           CALL PRS('Mouse not activated; please move mouse slightly'//
     $          'and re-enter$')
           GO TO 97
        ENDIF

c       filter x,y thru grid
        grid=param(18)
        if (grid.eq.0) grid=0.05
        if (grid.lt.0) grid=-1./grid
        if (grid.lt. .0099) ifgrdc=.false.

C       Filter Only if building and outside build menu area.
        iobjct=0
        if(xscreen.lt.1.0 .and.ifgrid.and.button.ne.'RIGHT') then
           call filter(xmouse,ymouse,xscreen,yscreen)
        elseif (xscreen.lt.1.0) then                  ! Button = right
           call filter(xmouse,ymouse,xscreen,yscreen) ! to latch to element
           xmouse=xphy(xscreen)
           ymouse=yphy(yscreen)
        else
           xmouse=xphy(xscreen)
           ymouse=yphy(yscreen)
        endif

        IF(BUTTON.EQ.'MIDDLE')THEN
           CALL PRSR('Cursor Coordinates: x: $',XMOUSE)
           CALL PRSR('Cursor Coordinates: y: $',YMOUSE)
           GO TO 20
        ELSE IF(BUTTON.EQ.'LEFT'.OR.BUTTON.EQ.'RIGHT')then
        ENDIF

        RETURN
        END
        
      SUBROUTINE MOVE(X,Y)
      CALL MOVESC(XSCR(X),YSCR(Y))
      RETURN
      END

      SUBROUTINE DRAW(X,Y)
      CALL DRAWSC(XSCR(X),YSCR(Y))
      RETURN
      END

      SUBROUTINE MOVESC(X,Y)
C     Moves with "Pen up", no line Drawn  IN SCREEN COORDINATES
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      COMMON /FILL/ ICFILL,IFILL,NFILL,XFILL(210),YFILL(210)

      IF(IFILL.NE.0) THEN
C          Currently Filling region defined by moves and draws
           IF(NFILL.GT.210)THEN
C             Forgot somewhere to shut off fill
              IFILL=0
              NFILL=0
              RETURN
           ENDIF
           NFILL = NFILL+1
           XFILL(NFILL) = SCOORDX(X)
           YFILL(NFILL) = SCOORDY(Y)
      ENDIF
      CALL XMOVE(SCOORDX(X),SCOORDY(Y))
      IF(IFPOSTS.AND.IFHARD.AND.IFHARD0)
     $WRITE(45,'(2F10.3,A3)')REAL(550-Y*500),REAL(X*500+50),' m '
c12/23/92/pff $WRITE(45,'(2F10.3,A6)')REAL(550-Y*500),REAL(X*500+50),' move '
      RETURN
      END

      SUBROUTINE DRAWSC(X,Y)
C
C     Draws lines for various devices with "Pen down" IN SCREEN COORDINATES
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      COMMON /FILL/ ICFILL,IFILL,NFILL,XFILL(210),YFILL(210)
C
      IF(IFILL.NE.0) THEN
C          Currently Filling region defined by moves and draws
           IF(NFILL.GT.210)THEN
C             Forgot somewhere to shut off fill
              IFILL=0
              NFILL=0
              RETURN
           ENDIF
           NFILL = NFILL+1
           XFILL(NFILL) = SCOORDX(X)
           YFILL(NFILL) = SCOORDY(Y)
      ENDIF
      CALL XDRAWTO(SCOORDX(X),SCOORDY(Y))
      IF(IFPOSTS.AND.IFHARD.AND.IFHARD0)
     $WRITE(45,'(2F9.3,A3)')REAL(550-Y*500),REAL(X*500+50),' d '
c12/23/92/pff $WRITE(45,'(2F9.3,A6)')REAL(550-Y*500),REAL(X*500+50),' draw '
      RETURN
      END

      SUBROUTINE BEGINP(X1,Y1)
C     Begins panel for fill  NO BOUNDARY DISPLAY
      LOGICAL IFTMP 
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      COMMON /FILL/ ICFILL,IFILL,NFILL,XFILL(210),YFILL(210)
C
      X=(X1-XZERO)/XFAC
      Y=(Y1-YZERO)/YFAC
C
      ICFILL=0
      IFILL = 1
      NFILL = 0
      CALL MOVEC(X1,Y1)
      IFTMP=IFHARD
      IFHARD=.FALSE.
      CALL MOVE(X1,Y1)
      IFHARD=IFTMP
      RETURN
      END

      SUBROUTINE BEGINB(X1,Y1)
C     Begins panel for fill WITH BOUNDARY DISPLAY
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      COMMON /FILL/ ICFILL,IFILL,NFILL,XFILL(210),YFILL(210)
      LOGICAL IFTMP 
      X=(X1-XZERO)/XFAC
      Y=(Y1-YZERO)/YFAC
      ICFILL=1
      IFILL = 1
      NFILL = 0
      CALL MOVEC(X1,Y1)
      IFTMP=IFHARD
      IFHARD=.FALSE.
      CALL MOVE(X1,Y1)
      IFHARD=IFTMP
      RETURN
      END

      subroutine gnable
      return
      end

      SUBROUTINE GINDIS
      RETURN
      END

      subroutine sgvis(iseg,ivisibility)
      RETURN
      END

      SUBROUTINE HEJECT
C     EJECTS page from hard copy laser
      INCLUDE 'devices.inc'
      IF(IFHARD.AND.IFHARD0) THEN
         IF(IFPOSTS)WRITE(45,*)' showpage'
         IF(IFPOSTS)CALL PENW(3)
         IF(IFLASE)WRITE(45,*)' ^,'
      ENDIF
      RETURN
      END

      SUBROUTINE CLSSEG
      RETURN
      END

      SUBROUTINE OPNSEG(ISEG)
      RETURN
      END

      SUBROUTINE DELSEG(ISEG)
      RETURN
      END

      SUBROUTINE PENW(IWIDTH)
C     Sets Pen width of laser. IWIDTH should be an odd integer 1 to 31
C     Default is 3.
      INCLUDE 'devices.inc'
      INTEGER IWDOLD
      SAVE    IWDOLD
      DATA    IWDOLD /-99999/ 
      IF (IWIDTH.EQ.IWDOLD) RETURN
C
      IWDOLD=IWIDTH
      IF(IFHARD.AND.IFHARD0) THEN
        IF(IFPOSTS)
     $    WRITE(45,'(f7.2,a14)')REAL(IWIDTH/15.),' setlinewidth '
        IF(IFLASE)WRITE(45,100) IWIDTH
      ENDIF
100   FORMAT(' ^PW',I2.2)
      RETURN
      END

      SUBROUTINE SCROLL(ILINES)
      RETURN
      END

      SUBROUTINE RUBBER
      RETURN
      END

      SUBROUTINE OFFRUB
      RETURN
      END

      SUBROUTINE GWRITE(X,Y,SIZE,STRING)
      CHARACTER STRING(100)
      COMMON /CLIP/ XMOVEC,YMOVEC,WT,WB,WL,WR
      CALL GSWRIT(XSCR(X),YSCR(Y),SIZE,STRING)     
c     if ((wl.lt.x.and.x.lt.wr) .and. (wb.lt.y.and.y.lt.wt)) 
c    $CALL GSWRIT(XSCR(X),YSCR(Y),SIZE,STRING)     
      RETURN
      END      

      SUBROUTINE GSWRIT(X,Y,SIZE,STRING)
      CHARACTER STRING(100),DIR,STRING2(100)
      COMMON /FONT/ IFONT,IDIR
      INCLUDE 'devices.inc'
      DATA IFONT,IDIR /204,1/
      NULL=CHAR(0)
      DO 1 I=1,100
          STRING2(i)=STRING(I)
          IF(STRING(I).EQ.'$') THEN
                STRING2(i)=NULL
                NCHARS=I-1
                GO TO 2
          ENDIF
1     CONTINUE
      CALL PRS('PUT DOLLAR SIGN AFTER LAST CHARACTER!$')
2     CONTINUE
      CALL XPUTXT(SCOORDX(X),SCOORDY(Y),SIZE,STRING2,NCHARS)
      IF(IFPOSTS.AND.IFHARD.AND.IFHARD0) THEN
          CALL MOVESC(x,y)
          write(45,*)' 90 rotate'
          write(45,'(80A1)') '(',(STRING(I),I=1,NCHARS),')',' ','s',
     $    'h','o','w'
          write(45,*)'-90 rotate'
      ENDIF
5     RETURN
      END

      SUBROUTINE DEVINI(DEVICE)
#     include "basics.inc"
      
      IFNOSEG= .TRUE.

      IFXWIN = .TRUE.

      CALL GETWIN(XSIZE,YSIZE)
      WINDOWW=XSIZE
      WINDOWH=YSIZE
      CALL XPLANES(NPLANES)
      IFBWGKS = .FALSE.
      IF(NPLANES .EQ. 1) IFBWGKS = .TRUE.
      RETURN
      END

      SUBROUTINE DEVEX
      INCLUDE 'devices.inc'
C
C     Set Screen to scroll
      IF(IFPOSTS)WRITE(45,*)' showpage'
      IF(IFPOSTS)       CLOSE(UNIT=45)
      CALL CLOSEW
      RETURN
      END
      SUBROUTINE ENDP
      INCLUDE 'devices.inc'
      COMMON /FILL/ ICFILL,IFILL,NFILL,XFILL(210),YFILL(210)
      IF (NFILL.GT.210) THEN
         NFILL=210 
         XFILL(NFILL)=XFILL(1)
         YFILL(NFILL)=YFILL(1)
      ENDIF
      CALL XPOLY(XFILL,YFILL,NFILL,ICFILL)
      IFILL = 0
      NDRAW = NFILL
      NFILL = 0
      RETURN
      END

      SUBROUTINE FILLP(ICOLOR)
C     SELECTS FILL panel
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      CALL XSETFILL(IABS(ICOLOR))
      RETURN
      END


      SUBROUTINE HMT_FILLP(ICOLOR)
C     SELECTS FILL panel
      CHARACTER STRING*5
      INCLUDE 'devices.inc'
      CALL HMT_XSETFILL(IABS(ICOLOR))
      RETURN
      END

      SUBROUTINE CLEAR
C     clears screen (graphics portion)
      CALL XCLEAR
      RETURN
      END

      SUBROUTINE CLEARA
C     clears screen
      CALL XCLEARA
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE COLOR(ICOLOR)
C     SETS COLOR OF LINES TO BE
C     WHICHEVER DRAW IS TO BE THIS COLOR
      EQUIVALENCE(CCOLOR,ICDUM)
      INCLUDE 'devices.inc'
      CHARACTER*3 COLORS(0:200),CCOLOR
      CHARACTER STRING*3
C
C     Save color for use with IFSAVE DRAW option.
C

c     write (6,*) 'COLOR :: ICOLOR = ',ICOLOR

      IF (ICOLOR.EQ.LSTCLR) RETURN
      LSTCLR=ICOLOR
C
      ICDUM = ICOLOR
      INDEX=ICOLOR
      IF(IFBWGKS .AND. ICOLOR .NE.0) INDEX=1
      NCOLOR=15
C
      COLORS(0)='BLACK'
      COLORS(1)='WHITE'
      COLORS(2)='RED'
      COLORS(3)='GREEN'
      COLORS(4)='BLUE'
      COLORS(5)='LTBLUE'
      COLORS(6)='PURPLE'
      COLORS(7)='YELLOW'
      COLORS(8)='ORANGE'
      COLORS(9)='GREEN'
      COLORS(10)='LTGREEN'
      COLORS(11)='MEDBLUE'
      COLORS(12)='VIOLET'
      COLORS(13)='PINK'
      COLORS(14)='GRAY'
      COLORS(15)='LGRAY'
C
      
      DO 1 I=0,NCOLOR
          IF(CCOLOR.EQ.COLORS(I)) INDEX=I
1     CONTINUE

      CALL XSETPEN(INDEX)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE HMT_COLOR(ICOLOR)
C     SETS COLOR OF LINES TO BE
C     WHICHEVER DRAW IS TO BE THIS COLOR
      EQUIVALENCE(CCOLOR,ICDUM)
      INCLUDE 'devices.inc'
      CHARACTER*3 COLORS(0:200),CCOLOR
      CHARACTER STRING*3
C
C     Save color for use with IFSAVE DRAW option.
C

C HMT TRACE
C      write (6,*) 'COLOR :: ICOLOR = ',ICOLOR
C
      CALL XSETPEN(ICOLOR)
      RETURN
C
      IF (ICOLOR.EQ.LSTCLR) RETURN
      LSTCLR=ICOLOR
C
      ICDUM = ICOLOR
      INDEX=ICOLOR
      IF(IFBWGKS .AND. ICOLOR .NE.0) INDEX=1
      NCOLOR=15
C
      COLORS(0)='BLACK'
      COLORS(1)='WHITE'
      COLORS(2)='RED'
      COLORS(3)='GREEN'
      COLORS(4)='BLUE'
      COLORS(5)='LTBLUE'
      COLORS(6)='PURPLE'
      COLORS(7)='YELLOW'
      COLORS(8)='ORANGE'
      COLORS(9)='GREEN'
      COLORS(10)='LTGREEN'
      COLORS(11)='MEDBLUE'
      COLORS(12)='VIOLET'
      COLORS(13)='PINK'
      COLORS(14)='GRAY'
      COLORS(15)='LGRAY'
C
      
      DO 1 I=0,NCOLOR
          IF(CCOLOR.EQ.COLORS(I)) INDEX=I
1     CONTINUE

      CALL XSETPEN(INDEX)
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE CLRMAP(TYPE)
      CHARACTER*4 TYPE
      REAL RED(0:15),GREEN(0:15),BLUE(0:15)
      REAL MESHMP(4,0:15),HEATMP(4,0:15)
C
C      DATA MESHMP / 0,      0,       0,       0,
C     $            1,      61440,   61440,   61440,
C     $            2,      61440,   0,       0,
C     $            3,      0,       61440,   0,
C     $            4,      0,       0,       61440,
C     $            5,      12288,   40960,   53248,
C     $            6,      40960,   12288,   53248,
C    $            7,      61440,   61440,   0,
C     $            8,      61440,   32748,   0,
C     $            9,      0,       61440,   0,
C     $            10,     0,       61440,   32768,
C     $            11,     12288,   12288,   53248,
C     $            12,     20480,   12288,   20480,
C     $            13,     61440,   0,       61440,
C     $            14,     20480,   20480,   20480,
C     $            15,     45056,   45056,   45056  /
C                         RED      GREEN    BLUE

C      DATA HEATMP / 0,      0,       0,       0,
C     $            1,      61440,   61440,   61440,
C     $            2,      61440,   0,       0,
C     $            3,      0,       61440,   0,
C     $            4,      0,       0,       61440,
C     $            5,      12288,   40960,   53248,
C     $            6,      40960,   12288,   53248,
C     $            7,      61440,   61440,   0,
C     $            8,      61440,   32748,   0,
C     $            9,      0,       61440,   0,
C     $            10,     0,       61440,   32768,
C     $            11,     12288,   12288,   53248,
C     $            12,     20480,   12288,   20480,
C     $            13,     61440,   0,       61440,
C     $            14,     20480,   20480,   20480,
C     $            15,     45056,   45056,   45056  /
C                         RED      GREEN    BLUE

      
      DATA MESHMP /    0.,      .1,    .1,   .1,
     $                 1.,      1.0,    1.0,   1.0,
     $                  2.,      1.0,    0.0,   0.0,
     $                  3.,      1.0,    1.0,	0.0,
     $                  4.,      0.5,    0.5,	0.5,
     $                  5.,      0.0,    .0,     1.0,
     $                  6.,      1.0,    0.0,	0.7,
     $                  7.,      1.0,    1.0,	0.0,
     $                  8.,      0.7,    1.0,	0.0,
     $                  9.,      0.2,    1.0,	0.2,
     $                  10.,     0.0,    1.0,	1.0,
     $                  11.,     0.0,    0.3,	1.0,
     $                  12.,     1.0,    0.0,	1.0,
     $                  13.,     1.0,    0.0,	0.2,
     $                  14.,     0.6,    0.6,	0.6,
     $                  15.,     0.5,    0.5,	0.7 /
C                               RED     Green   BLUE
C
C
      DATA HEATMP /     0.0,     0.0,    0.0,	0.0,
     $                  1.,      1.0,    1.0,   1.0,
     $                  2.,      0.0,    0.0,   1.0,
     $                  3.,      0.0,    0.0,	1.0,
     $                  4.,      0.0,    0.1,	1.0,
     $                  5.,      0.1,    0.2,	0.9,
     $                  6.,      0.2,    0.7,	0.8,
     $                  7.,      0.4,    0.9,	0.6,
     $                  8.,      0.6,    1.0,	0.4,
     $                  9.,      0.7,    1.0,	0.3,
     $                  10.,     0.8,    0.9,	0.2,
     $                  11.,     0.9,    0.7,	0.1,
     $                  12.,     1.0,    0.2,	0.0,
     $                  13.,     1.0,    0.1,	0.0,
     $                  14.,     1.0,    0.0,	0.0,
     $                  15.,     1.0,    0.0,	0.0 /
C                               RED     Green   BLUE
C

      IF(TYPE.EQ.'HEAT') THEN
          DO 10 i=0,15
             RED   (i)= HEATMP(2,I)
             GREEN (i)= HEATMP(3,I)
             BLUE  (i)= HEATMP(4,I)
10        CONTINUE
          DO 11 i=2,15
             RED(i)  = (i -2)/13.
             blue(i) = (15-i)/13.
             green(i)= amin1(blue(i),red(i)) *2.0
 11          CONTINUE
      ENDIF
      IF(TYPE.EQ.'MESH') THEN
          DO 20 i=0,15
             RED   (i)= MESHMP(2,I)
             GREEN (i)= MESHMP(3,I)
             BLUE  (i)= MESHMP(4,I)
20        CONTINUE
      ENDIF
C     scale colormap up by 61440
      DO 30 I=0,15
         RED(I)=RED(I)*61440.0
         GREEN(I)=GREEN(I)*61440.0
         BLUE(I)=BLUE(I)*61440.0
 30   CONTINUE
C     white text
      RED(1)=61440.0
      GREEN(1)=61440.0
      BLUE(1)=61440.0
      CALL HMT_XCMAP(RED,GREEN,BLUE,15)
      CALL COLOR(1)
      RETURN
      END

      FUNCTION SCOORDX(WORLDX)
#     include "basics.inc"

      SCOORDX=WORLDX/1.3*WINDOWW
      RETURN
      END

      FUNCTION SCOORDY(WORLDY)
#     include "basics.inc"
      SCOORDY=ABS(WORLDY-1.0)*(WINDOWH - 180.0)
      RETURN
      END
      SUBROUTINE OPENFLD(IUNIT,FILE,suffix,STATUS,IFFMAT,IERR)
      CHARACTER*(*) FILE,STATUS
      CHARACTER*30 NEWFILE,LFILE
      character*1  lfile1(30)
      character*2  s2
      character*3  s3
      equivalence (lfile1,lfile)
      integer suffix
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
c        print*,'renaming ',file,' ',newfile
c        MVERR=RENAME(FILE,NEWFILE)
      ENDIF
      if (suffix.ne.0) then
c
c        Append suffix to file name
c
c        len1=ltrunc(lfile,30)+1
         len1=indx1(lfile,'.fld',4)+4
         if (suffix.lt.100) then
            write(s2,32) suffix
            call chcopy(lfile1(len1),s2,2)
   32       format(i2.2)
         else
            write(s3,33) suffix
            call chcopy(lfile1(len1),s3,3)
   33       format(i3.3)
         endif
      endif
c
      if (iffmat) then
         write(6,*) 'Trying to open fld file:',LFILE,status,iffmat
         OPEN(UNIT=IUNIT,FILE=LFILE,STATUS=STATUS
     $       ,FORM='FORMATTED',ERR=1)
         IERR=0
c        IFFMAT=.TRUE.
         write(6,*) 'Successfully opened file: '
     $             ,lfile,' Format:',iffmat
         RETURN
      endif
c
    1 continue
      if (.not.iffmat) then
         write(6,*) 'Trying to open unfmt fld file:',LFILE,status
         OPEN(UNIT=IUNIT,FILE=LFILE,STATUS=STATUS
     $      ,FORM='UNFORMATTED',ERR=2)
         IERR=0
c        IFFMAT=.FALSE.
         write(6,*) 'Successfully opened file: '
     $             ,lfile,' Format:',iffmat
         RETURN
      endif
c
    2 continue
      IERR=1
c     IFFMAT=.TRUE.
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine mouse2(xmouse,ymouse,button,ifout)
#     include "basics.inc"
      integer*4 record(3),size,pet
      integer errind,mode,esw,xform,ldr,retsiz,istat,devno
      real locx,locy,earea(4)
      common/mousefac/ window
      integer e
      logical ifout

   20 continue

      if (.not.ifdemo) call xwaitm(ibutton,XM,YM)
c     if (     ifdemo) call xwaitm(ibutton,xm,ym)
      if (     ifdemo) read(55,*)ibutton,xm,ym
C     if (     ifdemo) read(5,*)ibutton,xm,ym
      if(iflearn)write(3,*)ibutton,xm,ym

      if(ibutton.eq.0) then
         call prs('Still waiting for mouse input$')
         goto 20
      endif

      if (ibutton.eq.1) button='LEFT'
      if (ibutton.eq.2) button='MIDDLE'
      if (ibutton.eq.3) button='RIGHT'

      window=1.3
      xscreen=(xm/windoww*window) 
      yscreen=(ym/(windowh-180.))-0.05
c     OLD VERSION -- coordinate transformation is almost device dependent!
c     yscreen=((ym-12)/1000.0*window) + 0.015
      if(xscreen .eq. oldxscreen .and.
     $   yscreen .eq. oldyscreen .and. 
     $   choice  .eq.'accept current swithces') then
         call prs('Error reading mouse input; please move mouse'//
     $            'slightly and re-enter$')
         goto 20
      endif

      oldxscreen = xscreen
      oldyscreen = yscreen

      if(xscreen .lt. 0.0 .or. xscreen .gt. 1.3  .or.
     $     yscreen .lt. 0.0 .or. yscreen .gt. 1.0) then
         call prs('Error reading mouse input; please re-enter$')
         goto 20
      endif

      if(xscreen .eq. 0.0 .and. yscreen .eq. 1.5E-2) then
         CALL PRS('Mouse not activated; please move mouse slightly'//
     $        'and re-enter$')
         goto 20
      endif

c     filter x,y thru grid
      grid=param(18)
      if (grid.eq.0) grid=0.05
      if (grid.lt.0) grid=-1./grid
      if (grid.lt. .0099) ifgrdc=.false.


      xmouse=xphy(xscreen)
      ymouse=yphy(yscreen)

c     Filter only if outside build menu area.
      if (xscreen.lt.1.0.and.button.ne.'RIGHT') then
         call filter(xmouse,ymouse,xscreen,yscreen)
      elseif (button.eq.'RIGHT') then ! Latch to closest point
         rmin = 1.e20
         do e=1,nel
         do i=1,4
            r=(x(i,e)-xmouse)**2+(y(i,e)-ymouse)**2
            if (r.lt.rmin) then
               rmin  = r
               iemin = e
               icmin = i
            endif
         enddo
         enddo
         xmouse = x(icmin,iemin)
         ymouse = y(icmin,iemin)
      else
         xmouse=xphy(xscreen)
         ymouse=yphy(yscreen)
      endif

      if(button.eq.'MIDDLE')then
         CALL PRSR('Cursor Coordinates: x: $',XMOUSE)
         CALL PRSR('Cursor Coordinates: y: $',YMOUSE)
         goto 20
      endif

      if (ifout) call prsrr('x,y:$',xmouse,ymouse)

      return
      end
c-----------------------------------------------------------------------
