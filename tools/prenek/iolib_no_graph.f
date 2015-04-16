C------------------------------------------------------------------------------
C
C                          NEKTON 2.6  2/8/90
C
C			Copyright (C) 1990, by the
C
C		Massachusetts Institute of Technology  and Nektonics, Inc.
C
C All Rights Reserved
C
C This program is a licenced product of MIT and Nektonics, Inc.,  and it is
C not to be disclosed to others, copied, distributed, or displayed
C without prior authorization.
C
C------------------------------------------------------------------------------
C
C
C     I/O Subroutine Library for reading & writing to screen device
C     There are two layers of I/O routines.  This upper layer has
C     about 10 routines for handling different formats.  Each of them
C     has one call to the lower layer.
C
C     Output Subroutine Convention:
C     PRS            Prints String                  call PRS   (S)
C     PRI            Prints Integer                 call PRI   (I)
C     PRR            Prints Real                    call PRR   (R)
C     PRRR           Prints Real,Real               call PRRR  (R,R)
C     PRII           Prints Integer,Integer         call PRII  (I,I)
C     PRIS           Prints Integer,String          call PRSI  (I,S)
C     PRSI           Prints String, Integer         call PRSI  (S,I)
C     PRSR           Prints String, Real            call PRSR  (S,R)
C     PRSIR          Prints String, Integer, Real   call PRSIR (S,I,R)
C     PRSIS          Prints String, Integer, String call PRSIS (S,I,S)
C     PRSRS          Prints String, Real,    String call PRSRS (S,R,S)
C
C     All strings must have $ in them to signal termination
C
C
C     Input Subroutine Convention:
C     RES           Reads String of Length N       call RES   (S,N)
C     REI           Reads Integer                  call REI   (I)
C     RER           Reads Real                     call RER   (R)
C     RERR          Reads Real, Real               call RERR  (R,R)
C     RERRR         Reads 3 Real                   call RERRR (R,R,R)
C     RERRRR        Reads 4 Real                   call RERRRR(R,R,R,R)
C
C
C     The lower layer (not seen by the applications code) contains two
C     routines, GETS and PUTS, which perform I/O with the device.
C
C     Driver Subroutine Calling Convention:
C     GETS         Reads  String of 130 chars        call GETS(S)
C     PUTS         Prints String of  N chars         call PUTS(S,N)

c-----------------------------------------------------------------------
      subroutine prs(s)
      include 'devices.inc'

      integer icalld
      save    icalld
      data    icalld /0/
      
      character s(*)
      character*1 s1(132),s2(132)
      character*132 s131,s132
      equivalence (s131,s1)
      equivalence (s132,s2)
      save s2
      data s2 /132*' '/
      integer scount
      save    scount
      data    scount /0/

      if (icalld.eq.0) open(unit=87,file='pretex.jou')
      icalld=1

      do 1 i=1,133
         if(s(i).eq.'$')then
            nchars=i-1
            goto 2
         endif
 1    continue

      call putsold('I/O Error: No String Terminator sent to PRS',43)
      write(6 ,*)  'I/O Error: No String Terminator sent to PRS'
      write(87,*)  'I/O Error: No String Terminator sent to PRS'

      ri=-i
      nchars=sqrt(ri)

      nchars = 132
      call chcopy (s1,s,nchars)
      write(87,3) (s1(j),j=1,nchars)

 2    continue

      if (ifgraf) call putsold  (s   ,nchars)
      call chcopy (s1,s,nchars)

      if (s131.eq.s132) then  ! Check for infinite loop
         scount = scount+1
         if (scount.gt.10) then
            write(6,*) 'ABORT: too many repeat calls to prs.'
            call exitt
         endif
      else
         scount = 0
      endif
      call chcopy (s2,s1,132)

      write(6 ,3) (s1(j),j=1,nchars)
      write(87,3) (s1(j),j=1,nchars)

 3    format(132a1)
      return
      end
c-----------------------------------------------------------------------
      subroutine pri(i)
      character S*132
C
      S=' '
      WRITE(S,'(I13)',ERR=1)I
      S(9:9)='$'
      call PRS(S)
 1    return
      end
C
c-----------------------------------------------------------------------
      subroutine prr(r)
      character S*132
C
      S=' '
      WRITE(S,'(G14.6)',ERR=1)R
      S(15:15)='$'
      call PRS(S)
 1    return
      end
C
C
c-----------------------------------------------------------------------
      subroutine prii(i1,i2)
      character S*132
C
      S=' '
      WRITE(S,'(2I13)',ERR=1)I1,I2
      S(17:17)='$'
      call PRS(S)
 1    return
      end
C
c-----------------------------------------------------------------------
      subroutine prrrr(r1,r2,r3)
      character S*132
C
      S=' '
      WRITE(S,'(3G13.9)',ERR=1)R1,R2,R3
      S(55:55)='$'
      call PRS(S)
 1    return
      end
C
c-----------------------------------------------------------------------
      subroutine prrr(r1,r2)
      character S*132
C
      S=' '
      WRITE(S,'(2G14.6)',ERR=1)R1,R2
      S(29:29)='$'
      call PRS(S)
 1    return
      end
C
c-----------------------------------------------------------------------
      subroutine prsi(s,i)
      character S(*),SS*132
C
      DO 1 J=1,70
         IF(S(J).EQ.'$')THEN
            NC=J-1
            goto 2
         ENDIF
         SS(J:J)=S(J)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSI',44)
      call putsold(s,132)
      return
 2    CONTINUE
C
      WRITE(SS(NC+1:NC+13),'(I13)',ERR=13)I
      SS(NC+14:NC+14)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSI',39)
      return
      end
C
c-----------------------------------------------------------------------
      subroutine pris(i,s)
      character s(*),ss*146

      call blank(ss,146)
      write(ss(1:13),'(i13)',err=13)i
      do 1 j=1,133
         if(s(j).eq.'$')then
            nc=j-1
            go to 2
         endif
         ss(j+14:j+14)=s(j)
 1    continue
      call putsold('I/O Error: No String Terminator sent to PRSI',44)

 2    continue
      ss(nc+14:nc+14)='$'
      call prs(ss)
      return

 13   call putsold('I/O Error: cant write to string in PRIS',39)
      return
      end
c-----------------------------------------------------------------------
      subroutine pris4(i,s)
      character s(*),ss*136

      call blank(ss,136)
      write(ss(1:4),'(i4)',err=8)i
      do 1 j=1,133
         if(s(j).eq.'$')then
            nc=j-1
            go to 2
         endif
         ss(j+5:j+5)=s(j)
 1    continue
      call putsold('I/O Error: No String Terminator sent to PRSI',44)

 2    continue
      ss(nc+5:nc+5)='$'
      call prs(ss)
      return

 8    call putsold('I/O Error: cant write to string in PRIS',39)
      return
      end
c-----------------------------------------------------------------------
      subroutine pris8(i,s)
      character s(*),ss*140

      call blank(ss,140)
      write(ss(1:8),'(i8)',err=8)i
      do 1 j=1,133
         if(s(j).eq.'$')then
            nc=j-1
            go to 2
         endif
         ss(j+9:j+9)=s(j)
 1    continue
      call putsold('I/O Error: No String Terminator sent to PRSI',44)

 2    continue
      ss(nc+9:nc+9)='$'
      call prs(ss)
      return

 8    call putsold('I/O Error: cant write to string in PRIS',39)
      return
      end
c-----------------------------------------------------------------------
      subroutine prsrrr(s,R1,R2,R3)
      character S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            goto 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSRR',44)
      call putsold(s,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+42),'(1p3e14.6)',ERR=13)r1,r2,r3
      SS(NC+43:NC+43)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSRR',39)
      return
      end
C-----------------------------------------------------------------------
      subroutine prsrr(s,R1,R2)
      character S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            goto 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSRR',44)
      call putsold(s,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+28),'(1p2e14.6)',ERR=13)r1,r2
      SS(NC+29:NC+29)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSRR',39)
      return
      end
C-----------------------------------------------------------------------
      subroutine prsr(s,r)
      character S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            goto 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSR',44)
      call putsold(s,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+14),'(G14.6)',ERR=13)R
      SS(NC+15:NC+15)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSR',39)
      return
      end
C-----------------------------------------------------------------------
      subroutine prsii(s,i,j)
      character S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            goto 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSii',44)
      call putsold(s,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+20),'(2I10)',ERR=13)I,J
      SS(NC+21:NC+21)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSii',39)
      return
      end
C
c-----------------------------------------------------------------------
      subroutine prsir(s,i,r)
      character S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            goto 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSir',44)
      call putsold(s,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+20),'(I6,G14.6)',ERR=13)I,R
      SS(NC+21:NC+21)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSir',39)
      return
      end
c-----------------------------------------------------------------------
      subroutine prsis(s1,i,s2)
      character S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 J=1,70
         IF(S1(J).EQ.'$')THEN
            NC=J-1
            goto 2
         ENDIF
         SS(J:J)=S1(J)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSIS',45)
      call putsold(sS,70)
      return
 2    CONTINUE
      NC = NC + 1
      WRITE(SS(NC+1:NC+8),'(I8)',ERR=13)I
C
      ISTART = NC+10
      IS2=0
      DO 11 J=ISTART,132
         IS2=IS2+1
         IF(S2(IS2).EQ.'$')THEN
            NC=J-1
            goto 12
         ENDIF
         SS(J:J)=S2(IS2)
 11   CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSIS',45)
      call putsold(sS,70)
      return
 12   CONTINUE
C
      SS(NC+1:NC+1)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSIS',40)
      return
      end
c-----------------------------------------------------------------------
      subroutine prsisi(s1,i1,s2,i2)
      character S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 J=1,70
         IF(S1(J).EQ.'$')THEN
            NC=J-1
            goto 2
         ENDIF
         SS(J:J)=S1(J)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSISI',46)
      call putsold(sS,70)
      return
 2    CONTINUE
      NC = NC + 1
      WRITE(SS(NC+1:NC+8),'(I8)',ERR=13) I1
C
      ISTART = NC+10
      IS2=0
      DO 11 J=ISTART,132
         IS2=IS2+1
         IF(S2(IS2).EQ.'$')THEN
            NC=J-1
            goto 12
         ENDIF
         SS(J:J)=S2(IS2)
 11   CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSISI',46)
      call putsold(sS,70)
      return
 12   CONTINUE
C
      NC = NC + 1
      WRITE(SS(NC+1:NC+8),'(I8)',ERR=13) I2
      NC = NC+13
      SS(NC+1:NC+1)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSISI',41)
      return
      end
c-----------------------------------------------------------------------
      subroutine prsrs(s1,r,s2)
      character S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 I=1,60
         IF(S1(I).EQ.'$')THEN
            NC=I-1
            goto 2
         ENDIF
         SS(I:I)=S1(I)
 1    CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSRS',45)
      call putsold(ss,70)
      return
 2    CONTINUE
C
C     Pad a few blanks
      NC=NC+1
      WRITE(SS(NC+1:NC+14),'(G14.6)',ERR=13)R
C
      ISTART = NC + 16
      IS2 = 0
      DO 11 I=ISTART,132
         IS2 = IS2 + 1
         IF(S2(IS2).EQ.'$')THEN
            NC=I-1
            goto 12
         ENDIF
         SS(I:I)=S2(IS2)
 11   CONTINUE
      call putsold('I/O Error: No String Terminator sent to PRSRS',45)
      call putsold(ss,70)
      return
 12   CONTINUE
      SS(NC+1:NC+1)='$'
      call PRS(SS)
      return
 13   call putsold('I/O Error: cant write to string in PRSRS',40)
      return
      end
C
C
C
C                                          Input Routine Section
C
C
c-----------------------------------------------------------------------
      subroutine res(s1,n)
      include 'devices.inc'
      character S1(*),S(132)

C     All text input comes thru here
C

      IF (.NOT.IFDEMO) THEN
         call GETSOLD(S)
      ELSE
         READ(55,132) S
  132    FORMAT(132A1)
      ENDIF
C
c     call GETSOLD(S)
      DO 1 I=1,N
         S1(I) = S(I)
 1    CONTINUE
      return
      end
C
c-----------------------------------------------------------------------
      subroutine reiii(I,J,K)
C     Read Integer
      character*132 S
      integer count
      count=0
C
 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13) i,j,k
      write (6,*) i,j,k
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter 3 integer Values$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop
      end

c-----------------------------------------------------------------------
      subroutine reii(i,J) !     Read Integer
      character*132 S
      integer count
      count=0

 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)I,j
      write (6,*) i,j
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter 2 integer Values$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end
c
c-----------------------------------------------------------------------
      subroutine rei(i) ! read Integer
      character*132 s
      integer count
      count=0

 1    call res(s,132)
      rewind(13)
      write (13,'(a132)')s
      rewind(13)
      read  (13,*,err=13,end=13)i
      write (6,*) i
      rewind(13)
      return

 13   call prs('Error reading input.  Enter Integer Value$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end

c-----------------------------------------------------------------------
      subroutine rer(r1)
C     Read Real
      character*132 S
      integer count
      count=0
C
 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1
      write (6,*) r1
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter Real    Value$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end
C
c-----------------------------------------------------------------------
      subroutine rerr(r1,r2)
C     Read Real
      character*132 S
      integer count
      count=0
C
 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2
      write (6,*) r1,r2
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter 2 Real Values$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end
C
c-----------------------------------------------------------------------
      subroutine rerrr(r1,r2,r3)
C     Read Real
      character*132 S
      integer count
      count=0
C
 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2,R3
      write (6,*) r1,r2,r3
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter 3 Real Values$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end
C
c-----------------------------------------------------------------------
      subroutine rerrrr(r1,r2,r3,r4)
C     Read Real
      character*132 S
      integer count
      count=0
C
 1    call RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2,R3,R4
      write (6,*) r1,r2,r3,r4
      REWIND(13)
      return
 13   call PRS('Error reading input.  Enter 4 Real Values$')

      write(6,*) s
      count=count+1
      if (count.lt.10) goto 1
      i=1/(count-10)
      stop

      end
c-----------------------------------------------------------------------
      subroutine putsold(s,nchars)
C     putsold is the one device-dependent output routine.
C     It Goes in Tekplot.f (Or Xinterface.f)  This is the Tek version
C     It displays a string on the output display device
      character S(*)
      WRITE(6,'(1X,132A1)',ERR=1)(S(I),I=1,NCHARS)
 1    return
      end
c-----------------------------------------------------------------------
      subroutine getsold(s)
C     GETSOLD is the one device-dependent input routine.
C     It Goes in Tekplot.f (Or Xinterface.f)  This is the Tek version
C     It returns an 132 character string entered from the input device
C
      character S*132
C
      S=' '
      READ(5,'(A132)',ERR=1,END=1)S
 1    return
      end
C-----------------------------------------------------------------------
      subroutine prsiii(s,i1,I2,I3)
      character S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            goto 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      call puts('I/O Error: No String Terminator sent to PRSiii',44)
      call puts(S,132)
      return
 2    CONTINUE
      WRITE(SS(NC+1:NC+30),'(3I10)',ERR=13)i1,i2,i3
      SS(NC+31:NC+31)='$'
      call PRS(SS)
      return
 13   call puts('I/O Error: cant write to string in PRSii',39)
      return
      end
C-----------------------------------------------------------------------
      subroutine priii(i1,i2,i3)
      character S*132
      S=' '
      WRITE(S,'(3I13)',ERR=1)I1,I2,I3
      S(25:25)='$'
      call PRS(S)
 1    return
      end
c-----------------------------------------------------------------------
      subroutine setgraph(ifgraf)
      logical ifgraf
      ifgraf = .false.
      return
      end
c-----------------------------------------------------------------------
