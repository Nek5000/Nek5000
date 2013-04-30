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
C     PRS            Prints String                  CALL PRS   (S)
C     PRI            Prints Integer                 CALL PRI   (I)
C     PRR            Prints Real                    CALL PRR   (R)
C     PRRR           Prints Real,Real               CALL PRRR  (R,R)
C     PRII           Prints Integer,Integer         CALL PRII  (I,I)
C     PRIS           Prints Integer,String          CALL PRSI  (I,S)
C     PRSI           Prints String, Integer         CALL PRSI  (S,I)
C     PRSR           Prints String, Real            CALL PRSR  (S,R)
C     PRSIR          Prints String, Integer, Real   CALL PRSIR (S,I,R)
C     PRSIS          Prints String, Integer, String CALL PRSIS (S,I,S)
C     PRSRS          Prints String, Real,    String CALL PRSRS (S,R,S)
C
C     All strings must have $ in them to signal termination
C
C
C     Input Subroutine Convention:
C     RES           Reads String of Length N       CALL RES   (S,N)
C     REI           Reads Integer                  CALL REI   (I)
C     RER           Reads Real                     CALL RER   (R)
C     RERR          Reads Real, Real               CALL RERR  (R,R)
C     RERRR         Reads 3 Real                   CALL RERRR (R,R,R)
C     RERRRR        Reads 4 Real                   CALL RERRRR(R,R,R,R)
C
C
C     The lower layer (not seen by the applications code) contains two
C     routines, GETS and PUTS, which perform I/O with the device.
C
C     Driver Subroutine Calling Convention:
C     GETS         Reads  String of 130 chars        CALL GETS(S)
C     PUTS         Prints String of  N chars        CALL PUTS(S,N)

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
      SUBROUTINE PRI(I)
      CHARACTER S*132
C
      S=' '
      WRITE(S,'(I13)',ERR=1)I
      S(9:9)='$'
      CALL PRS(S)
 1    RETURN
      END
C
      SUBROUTINE PRR(R)
      CHARACTER S*132
C
      S=' '
      WRITE(S,'(G14.6)',ERR=1)R
      S(15:15)='$'
      CALL PRS(S)
 1    RETURN
      END
C
C
      SUBROUTINE PRII(I1,I2)
      CHARACTER S*132
C
      S=' '
      WRITE(S,'(2I13)',ERR=1)I1,I2
      S(17:17)='$'
      CALL PRS(S)
 1    RETURN
      END
C
      SUBROUTINE PRRRR(R1,R2,R3)
      CHARACTER S*132
C
      S=' '
      WRITE(S,'(3G13.9)',ERR=1)R1,R2,R3
      S(55:55)='$'
      CALL PRS(S)
 1    RETURN
      END
C
      SUBROUTINE PRRR(R1,R2)
      CHARACTER S*132
C
      S=' '
      WRITE(S,'(2G14.6)',ERR=1)R1,R2
      S(29:29)='$'
      CALL PRS(S)
 1    RETURN
      END
C
      SUBROUTINE PRSI(S,I)
      CHARACTER S(*),SS*132
C
      DO 1 J=1,70
         IF(S(J).EQ.'$')THEN
            NC=J-1
            GO TO 2
         ENDIF
         SS(J:J)=S(J)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSI',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
C
      WRITE(SS(NC+1:NC+13),'(I13)',ERR=13)I
      SS(NC+14:NC+14)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSI',39)
      RETURN
      END
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
      SUBROUTINE PRSRRR(S,r1,r2,r3)
      CHARACTER S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            GO TO 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSRR',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+42),'(1p3e14.6)',ERR=13)r1,r2,r3
      SS(NC+43:NC+43)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSRR',39)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PRSRR(S,r1,r2)
      CHARACTER S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            GO TO 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSRR',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+28),'(1p2e14.6)',ERR=13)r1,r2
      SS(NC+29:NC+29)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSRR',39)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PRSR(S,R)
      CHARACTER S(*),SS*132
C
      DO 1 I=1,60
         IF(S(I).EQ.'$')THEN
            NC=I-1
            GO TO 2
         ENDIF
         SS(I:I)=S(I)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSR',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+14),'(G14.6)',ERR=13)R
      SS(NC+15:NC+15)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSR',39)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PRSII(S,I,J)
      CHARACTER S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            GO TO 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSii',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+20),'(2I10)',ERR=13)I,J
      SS(NC+21:NC+21)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSii',39)
      RETURN
      END
C
      SUBROUTINE PRSIR(S,I,R)
      CHARACTER S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            GO TO 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSir',44)
      CALL PUTSOLD(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+20),'(I6,G14.6)',ERR=13)I,R
      SS(NC+21:NC+21)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSir',39)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PRSIS(S1,I,S2)
      CHARACTER S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 J=1,70
         IF(S1(J).EQ.'$')THEN
            NC=J-1
            GO TO 2
         ENDIF
         SS(J:J)=S1(J)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSIS',45)
      CALL PUTSOLD(SS,70)
      RETURN
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
            GO TO 12
         ENDIF
         SS(J:J)=S2(IS2)
 11   CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSIS',45)
      CALL PUTSOLD(SS,70)
      RETURN
 12   CONTINUE
C
      SS(NC+1:NC+1)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSIS',40)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PRSISI(S1,I1,S2,I2)
      CHARACTER S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 J=1,70
         IF(S1(J).EQ.'$')THEN
            NC=J-1
            GO TO 2
         ENDIF
         SS(J:J)=S1(J)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSISI',46)
      CALL PUTSOLD(SS,70)
      RETURN
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
            GO TO 12
         ENDIF
         SS(J:J)=S2(IS2)
 11   CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSISI',46)
      CALL PUTSOLD(SS,70)
      RETURN
 12   CONTINUE
C
      NC = NC + 1
      WRITE(SS(NC+1:NC+8),'(I8)',ERR=13) I2
      NC = NC+13
      SS(NC+1:NC+1)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSISI',41)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE PRSRS(S1,R,S2)
      CHARACTER S1(*),S2(*),SS*132
C
      SS = ' '
      DO 1 I=1,60
         IF(S1(I).EQ.'$')THEN
            NC=I-1
            GO TO 2
         ENDIF
         SS(I:I)=S1(I)
 1    CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSRS',45)
      CALL PUTSOLD(SS,70)
      RETURN
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
            GO TO 12
         ENDIF
         SS(I:I)=S2(IS2)
 11   CONTINUE
      CALL PUTSOLD('I/O Error: No String Terminator sent to PRSRS',45)
      CALL PUTSOLD(SS,70)
      RETURN
 12   CONTINUE
      SS(NC+1:NC+1)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTSOLD('I/O Error: cant write to string in PRSRS',40)
      RETURN
      END
C
C
C
C                                          Input Routine Section
C
C
      SUBROUTINE RES(S1,N)
      include 'devices.inc'
      CHARACTER S1(*),S(132)
C
C     All text input comes thru here
C
C
      IF (.NOT.IFDEMO) THEN
         CALL GETSOLD(S)
      ELSE
         READ(55,132) S
  132    FORMAT(132A1)
      ENDIF
C
c     CALL GETSOLD(S)
      DO 1 I=1,N
         S1(I) = S(I)
 1    CONTINUE
      RETURN
      END
C
      SUBROUTINE REiii(i,j,k)
C     Read Integer
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13) i,j,k
      write (6,*) i,j,k
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter 3 integer Values$')
      GO TO 1
      END
C
      SUBROUTINE REIi(I,j)
C     Read Integer
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)I,j
      write (6,*) i,j
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter 2 integer Values$')
      GO TO 1
      END
c
      SUBROUTINE REI(I)
C     Read Integer
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)I
      write (6,*) i
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter Integer Value$')
      GO TO 1
      END
C
      SUBROUTINE RER(R1)
C     Read Real
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1
      write (6,*) r1
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter Real    Value$')
      GO TO 1
      END
C
      SUBROUTINE RERR(R1,R2)
C     Read Real
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2
      write (6,*) r1,r2
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter 2 Real Values$')
      GO TO 1
      END
C
      SUBROUTINE RERRR(R1,R2,R3)
C     Read Real
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2,R3
      write (6,*) r1,r2,r3
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter 3 Real Values$')
      GO TO 1
      END
C
      SUBROUTINE RERRRR(R1,R2,R3,R4)
C     Read Real
      CHARACTER*132 S
C
 1    CALL RES(S,132)
      REWIND(13)
      WRITE (13,'(A132)')S
      REWIND(13)
      READ  (13,*,ERR=13,END=13)R1,R2,R3,R4
      write (6,*) r1,r2,r3,r4
      REWIND(13)
      RETURN
 13   CALL PRS('Error reading input.  Enter 4 Real Values$')
      GO TO 1
      END
C
      SUBROUTINE PUTSOLD(S,NCHARS)
C     PUTSOLD is the one device-dependent output subroutine.
C     It Goes in Tekplot.f (Or Xinterface.f)  This is the Tek version
C     It displays a string on the output display device
      CHARACTER S(*)
      WRITE(6,'(1X,132A1)',ERR=1)(S(I),I=1,NCHARS)
 1    RETURN
      END
C
      SUBROUTINE GETSOLD(S)
C     GETSOLD is the one device-dependent input subroutine.
C     It Goes in Tekplot.f (Or Xinterface.f)  This is the Tek version
C     It returns an 132 character string entered from the input device
C
      CHARACTER S*132
C
      S=' '
      READ(5,'(A132)',ERR=1,END=1)S
 1    RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PRSIII(S,I1,i2,i3)
      CHARACTER S(*),SS*132
C
      DO 1 k=1,60
         IF(S(k).EQ.'$')THEN
            NC=k-1
            GO TO 2
         ENDIF
         SS(k:k)=S(k)
 1    CONTINUE
      CALL PUTS('I/O Error: No String Terminator sent to PRSiii',44)
      CALL PUTS(S,132)
      RETURN
 2    CONTINUE
      WRITE(SS(NC+1:NC+30),'(3I10)',ERR=13)i1,i2,i3
      SS(NC+31:NC+31)='$'
      CALL PRS(SS)
      RETURN
 13   CALL PUTS('I/O Error: cant write to string in PRSii',39)
      RETURN
      END
C-----------------------------------------------------------------------
      SUBROUTINE PRIII(I1,I2,I3)
      CHARACTER S*132
      S=' '
      WRITE(S,'(3I13)',ERR=1)I1,I2,I3
      S(25:25)='$'
      CALL PRS(S)
 1    RETURN
      END
c-----------------------------------------------------------------------
      subroutine setgraph(ifgraf)
      logical ifgraf
      ifgraf = .false.
      return
      end
c-----------------------------------------------------------------------
