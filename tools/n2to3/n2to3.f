c----------------------------------------------------------------------
      program trans
c
c     Stretch 2d data set to 3d       Paul F. Fischer  (pff@cfm.brown.edu)  9/6/95.
c     Updated to include circle sides                  (fischer@mcs.anl.gov) 3/29/99.
c
c
c     This code takes as input a 2D nekton file, blah.rea, and an optional
c     .fld file, blah.fld  (asci format -- param(67)=0 -- only). 
c
c     It generates a correpsonding translated 3D mesh of height Z-max, 
c     with an equi-partitioned number of levels in Z, given by nlev.
c     Both these are input by the user.
c
c     The code generates periodic-bc's in the Z-direction, all other bc's
c     are preserved.
c
c     Run n2to3.   Enter blah,  enter out, enter nlev and zmax.
c     
c
      character*80 file,fname
      character*1  file1(80)
      equivalence (file1,file)
      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)

      parameter(nelm=9999)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      logical ifflow,ifheat
      common /arraz/ zmin,zmax,dz(nelm)

      common /arral/ ifcirc
      logical ifcirc

      integer e,en

      write(6,*)
      write(6,*) 'This is the code that establishes proper'
      write(6,*) 'element connectivities, as of Oct., 2006.'
      write(6,*)
c
c
c     for workstation:
      in = 5
c
c     for delta:
c     open(unit=7,file='indat',status='old',err=1999)
c     in = 7
c1999 continue
c
c     Get file name
      write(6,*) 'Input old (source) file name:'      
      call blank(file,80)
      read(in,80) file
      len = ltrunc(file,80)
   80 format(a80)
c
c     Get file name
      write(6,*) 'Input new (output) file name:'      
      call blank(fout,80)
      read(in,80) fout
      lou = ltrunc(fout,80)
c
      nlev = 1
      write(6,*) 
     $'input number of levels: (1, 2, 3,...; < 0 for circular sweep.):'
      read(in,*) nlev

      ifcirc = .false.
      if (nlev.lt.0) ifcirc=.true.
      nlev = abs(nlev)

      write(6,*) 'input z min:'
      read(in,*) zmin
      write(6,*) 'input z max:'
      read(in,*) zmax
      write(6,*) 
     $  'input gain (0=custom,1=uniform,other=geometric spacing):'
      read(in,*) gain
      if (gain.gt.0) then
         dz(1) = 1.
         zsum  = dz(1)
         do i = 2,nlev
            dz(i) = gain*dz(i-1)
            zsum  = zsum + dz(i)
         enddo
         ratio = (zmax-zmin)/zsum
         do i=1,nlev
            dz(i) = ratio*dz(i)
         enddo
      else ! Custom distribution
         write(6,*) 'input file containing ASCENDING z-values:'
         read (5,80) fname
         open(unit=61,file=fname,status='old')
         nlev = 0
         read(61,*) z0
         zmin = z0
         do i=1,nelm
            read(61,*,end=88) z1
            zmax = z1
            nlev = nlev + 1
            dz(nlev) = z1-z0
            z0 = z1
         enddo
   88    continue
         write(6,*) 'Found',nlev,' levels, z=',zmin,' to',zmax
      endif
c
c     stretch .rea data 
c
      call chcopy(file1(len+1),'.rea',4)
      call chcopy(fout1(lou+1),'.rea',4)
c
      open(unit=10, file=file)
      open(unit=11, file=fout)
      call rea23(dz,zmin,neln)
c
      write(6,*)
      write(6,6) neln,(fout1(k),k=1,len+4)
    6 format(i8,' elements written to ',40a1)
c
      close (unit=10)
      close (unit=11)
c
      write(6,*)
      write(6,6) neln,(fout1(k),k=1,len+4)
c
      close (unit=10)
      close (unit=11)
  999 continue
      stop
 9999 continue
      write(6,*) 'could not find file "indat". abort.'
      stop
      end
c-----------------------------------------------------------------------
      subroutine rea23(dzi,zmin,neln)
      real dzi(1)
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
      character*6  string21
      equivalence (string21,string1(19))
      character*4  string30
      equivalence (string30,string1(29))

      character*1 abc(52)
      character*52 abcs
      equivalence (abc,abcs)
      data abcs /'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvw'/
      save abc

      character*3  cb5,cb6
      character*3  a1
      character*3  option
      logical      ifcem,ifper, ifpec,ifpmc,ifpml

c     Nekton stuff
      parameter(nelm=9999)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      logical ifflow,ifheat

      real xc(4),yc(4),zc(4)

      common /arral/ ifcirc
      logical ifcirc

      integer e,en


      call rzero(x,4*nelm)
      call rzero(y,4*nelm)
      call rzero(bc,30*nelm)
      call blank(cbc,18*nelm)
c
c
c     Choose BC for Z direction
c
    
      write(6,*) 'This is for CEM: yes or no:'
      read (5,1) option
    1 format(a3)

      ifcem =.false.
      if (option.eq.'yes')   ifcem = .true.

      if (ifcem) then
      write(6,*) 'Enter Z (5) boundary condition (P,PEC,PMC,PML):'
      read (5,1) a1           
      call blank(cb5,3)
      write(cb5,3) a1
    3 format(a3)
      
      ifper = .false.
      if (a1.eq.'P  ')   ifper = .true.
      if (a1.eq.'PEC')   ifpec = .true.
      if (a1.eq.'PMC')   ifpmc = .true.
      if (a1.eq.'PML')   ifpml = .true.

      if (.not.ifper) then
         write(6,*) 'Enter Z (6) boundary condition (PEC,PMC,PML):'
         read (5,1) a1
         call blank(cb6,3)
         write(cb6,3) a1
      endif

      else

      write(6,*) 'Enter Z (5) boundary condition (P,v,O):'
      read (5,1) a1
      call blank(cb5,3)
      write(cb5,3) a1

      ifper = .false.
      if (a1.eq.'P  ')   ifper = .true.

      if (.not.ifper) then
         write(6,*) 'Enter Z (6) boundary condition (v,O):'
         read (5,1) a1
         call blank(cb6,3)
         write(cb6,3) a1
      endif

      endif

      write(6,*) 'this is cbz: ',cb5, cb6,' <--- '

c   
c
c     Read parameters from old .rea file
c
      call readwrite(string,'NEKTON',6)
      read(10,80) string
      write(11,'(2x,a18)') ' 3 DIMENSIONAL RUN'

      call readwrite(string,'NPSCAL',6)
      read (string,*) rpscal
      npscal = (rpscal + 0.01)

      ifflow = .false.
      call readwrite(string,'IFFLOW',6)
      if (indx1(string,' T ',3).ne.0) ifflow = .true.

      ifheat = .false.
      call readwrite(string,'IFHEAT',6)
      if (indx1(string,' T ',3).ne.0) ifheat = .true.

      call readwrite(string,'XFAC,YFAC',9)
      read(10,80) string
      write(11,10)
   10 format(1x,
     $'**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8')
c
c     Read mesh data
c
      read (10,*) nel,ndim
      ndim3 = 3
      if (nel.gt.nelm) then
         write(6,*) 'ABORT:  increase nelm in n2to3.f',nel,nelm
         call exit
      endif
      neln = nlev*nel
      write(11,11) neln,ndim3,neln
   11 format(3i10,11x,'NEL,NDIM,NELV')
c
c     Read & write xy data 
      dz = dzi(1)
      z0 = zmin
      z1 = z0+dz
      do e=1,nel      
         read (10,80) string
c        ca(e) = string1(32)
         ie52 = mod1(e,52)
         ca(e) = abc(ie52)
         len = ltrunc(string,80)
         read (10,*)   (x(k,e),k=1,4)
         read (10,*)   (y(k,e),k=1,4)
      enddo
      len = ltrunc(string,80)

      igroup = 0       ! Later, this will reflect the true group

      z1 = zmin
      do ilev=1,nlev

         z0 = z1
         z1 = z1+dzi(ilev)

         do e=1,nel      

            en = e + nel*(ilev-1)
            call out_e_hdr(en,ilev,ca(e),igroup,11)

            if (ifcirc) then ! Sweep in circular arc

               call sweep_circ(xc,yc,zc,x(1,e),y(1,e),4,z0)             ! z0=theta
               write(11,90)   (xc(k),k=1,4),(yc(k),k=1,4),(zc(k),k=1,4)

               call sweep_circ(xc,yc,zc,x(1,e),y(1,e),4,z1)             ! z1=theta
               write(11,90)   (xc(k),k=1,4),(yc(k),k=1,4),(zc(k),k=1,4)

            else   ! Translate data in z direction (std n2to3)

               write(11,90)  (x(k,e),k=1,4),(y(k,e),k=1,4),(z0,k=1,4)
               write(11,90)  (x(k,e),k=1,4),(y(k,e),k=1,4),(z1,k=1,4)

            endif

   90       format(4e14.6)

         enddo
      enddo


c     Curve sides (added 3/26/99   pff)

      call rdcurve

      if (ifcirc) then
         call newcurve
      else
         call out_curve
      endif

      read (10,80) string
      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c
      if (.not.ifflow) then
         read (10,80) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
      endif
c
      nbc = 0
      if (ifflow) nbc=nbc+1
      if (ifheat) nbc=nbc+1
                  nbc=nbc+npscal

      write(6,*) ifheat,npscal,nbc,' ifheat'

c     Read and write boundary conditions

      do ibc=1,nbc

         read (10,80) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
c
         do e = 1,nel
            if (nel.lt.1000) then
               do  k = 1,4
                  read (10,20) cbc(k,e),id,jd,(bc(j,k,e),j=1,5)
               enddo
            else
               do  k = 1,4
                  read (10,21) cbc(k,e),id,jd,(bc(j,k,e),j=1,5)
               enddo
            endif

            call rzero(bc(1,5,e),5)
            if     (ifper) then
               bc(1,5,e) = e+(nlev-1)*nel
               bc(2,5,e) = 6
               cbc(5,e) = 'P  '
            else
               cbc(5,e) = cb5
            endif

            call rzero(bc(1,6,e),5)
            bc(1,6,e) = e+nel
            bc(2,6,e) = 5
            cbc(6,e)  = 'E  '

            if (nlev.eq.1) then
               if     (ifper) then
                  bc(1,6,e) = e
                  cbc(6,e)  = 'P  '
               else
                  cbc(6,e) = cb6
                  call rzero(bc(1,6,e),5)
               endif
            endif

            do k=1,6
               if (neln.lt.1000) then
                  write(11,20) cbc(k,e),e,k,(bc(j,k,e),j=1,5)
               elseif (neln.lt.100000) then
                  write(11,21) cbc(k,e),id,k,(bc(j,k,e),j=1,5)
               else
                  write(11,22) cbc(k,e),id,k,(bc(j,k,e),j=1,5)
               endif
            enddo
            cbc(5,e) = 'E  '
            bc (1,5,e) = e   !  -nel
c
         enddo

         do ilev = 2,nlev
            do e = 1,nel
               id = e + nel*(ilev-1)
c              Periodic bc's on Z plane
               if (ilev.eq.nlev) then
                 ! bcs for final level
                 if     (ifper) then
                   cbc(6,e) = 'P  '
                   bc(1,6,e) = e
                 else
                   cbc(6,e) = cb6
                   call rzero(bc(1,6,e),5)
                 endif
               else
                 ! bcs for intermediate level
                 cbc(6,e) = 'E  '
                 call rzero(bc(1,6,e),5)
                 bc(1,6,e) = id+nel
                 bc(2,6,e) = 5
               endif

               cbc(5,e) = 'E  '
               call rzero(bc(1,5,e),5)
               bc(1,5,e) = id-nel
               bc(2,5,e) = 6

               do  k = 1,4
                 if (cbc(k,e).eq.'P  ') bc(1,k,e) = bc(1,k,e)+nel
                 if (cbc(k,e).eq.'E  ') bc(1,k,e) = bc(1,k,e)+nel
               enddo
               do  k = 1,6
                  if (neln.lt.1000) then
                     write(11,20) cbc(k,e),id,k,(bc(j,k,e),j=1,5)
                  elseif (neln.lt.100000) then
                     write(11,21) cbc(k,e),id,k,(bc(j,k,e),j=1,5)
                  else
                     write(11,22) cbc(k,e),id,k,(bc(j,k,e),j=1,5)
                  endif
               enddo
            enddo
         enddo
c
         if (cb5.eq.'v  ') cb5='t  '
         if (cb6.eq.'v  ') cb6='t  '
c
      enddo
c
      call readwrite(string,'endendend',9)
c
   20 FORMAT(1x,A3,2I3,5G14.7)
   21 FORMAT(1x,A3,i5,i1,5G14.7)
   22 FORMAT(1x,A3,i6,i1,5G14.7)

c
   80 format(a80)
   81 format(80a1)
      return
      end
c-----------------------------------------------------------------------
      subroutine blank(s,n)
      character*1 s(1)
      do i=1,n
        s(i)=' '
      enddo
      return
      end
c-----------------------------------------------------------------------
      function ltrunc(s,n)
      character*1 s(1)
      ltrunc = 0
      do j=n,1,-1
         if (s(j).ne.' ') then
            ltrunc = j 
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      INDX1=0
      IF (N1.LT.1) RETURN
C
      DO 300 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            RETURN
         ENDIF
300   CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine readwrite(sout,key,nk)

      character*80 sout,key

      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
      character*1 skey(80)

   80 format(a80)
   81 format(80a1)

      call blank (skey,80)
      call chcopy(skey,key,nk)

      do i=1,90000      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)

c        write(6,*) string,(skey(j),j=1,nk)

         call chcopy  (sout,string,80)
         call capit(sout,80)
         call capit(skey,nk)
         if (indx1(sout,skey,nk).ne.0) return

      enddo
  100 continue
      return
      end
c-----------------------------------------------------------------------
      function glmax(a,n)
      real a(1)
      tmax=-99.0E20
      do 100 i=1,n
         tmax=max(tmax,a(i))
  100 continue
      glmax=tmax
      return
      end
c-----------------------------------------------------------------------
      function glmin(a,n)
      real a(1)
      tmin=99.0E20
      do 100 i=1,n
         tmin=min(tmin,a(i))
  100 continue
      glmin = tmin
      return
      end
c-----------------------------------------------------------------------
      subroutine addc(x,c,n)
      real x(1)
      do i=1,n
         x(i) = x(i)+c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cmult(x,c,n)
      real x(1)
      do i=1,n
         x(i) = x(i)*c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(x,y,n)
      character*1 x(1),y(1)
      do i=1,n
         x(i) = y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rzero(x,n)
      real x(1)
      do i=1,n
         x(i) = 0.0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rdcurve
C
      parameter(nelm=9999)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      character*1 ans
      integer e


C     .Read formatted curve side data 

      read(10,*)
      read(10,*) ncurve
      call rzero(curve ,72*NELM)
      call blank(ccurve,12*NELM)

      IF (NCURVE.GT.0) THEN
         DO 50 ICURVE=1,NCURVE
            IF (nel.LT.1000) THEN
               READ(10,60,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSE
               READ(10,61,ERR=500,END=500) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ENDIF
            CURVE (1,IEDG,IEG)=R1
            CURVE (2,IEDG,IEG)=R2
            CURVE (3,IEDG,IEG)=R3
            CURVE (4,IEDG,IEG)=R4
            CURVE (5,IEDG,IEG)=R5
            CCURVE(  IEDG,IEG)=ANS
   50    CONTINUE
   60    FORMAT(I3,I3,5G14.6,1X,A1)
   61    FORMAT(I2,I6,5G14.6,1X,A1)
      ENDIF
      RETURN
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      CALL EXIT
      RETURN
C
      END
c-----------------------------------------------------------------------
      subroutine out_curve
C
      parameter(nelm=9999)
      common /arraz/ zmin,zmax,dz(nelm)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      CHARACTER*1 ANS
c
C     .Read formatted curve side data 
C
      neln = nlev*nel
      ncun = 2*nlev*ncurve
      write(11,11)
      write(11,12) ncun
   11 format(' ***** CURVED SIDE DATA *****')
   12 format(i8
     $ ,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')

      ilev = 0
      r31 = zmin
      r30 = r31
      r31 = r31 + dz(ilev+1)

      IF (NCURVE.GT.0) THEN
         DO 55 ilev=1,nlev
           DO 50 ieg=1,nel
           DO 50 iedg=1,8
            r1=CURVE (1,IEDG,IEG)
            r2=CURVE (2,IEDG,IEG)
            r3=CURVE (3,IEDG,IEG)
            r4=CURVE (4,IEDG,IEG)
            r5=CURVE (5,IEDG,IEG)
            ans=CCURVE( IEDG,IEG)
            ie =ieg + nel*(ilev-1)
            if (ans.eq.'C') then
               IF (neln.LT.1000) THEN
                  write(11,60) IEDG,IE,R1,R2,R3,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,60) IED4,IE,R1,R2,R3,R4,R5,ANS
               ELSEIF (neln.lt.1000000) then
                  write(11,61) IEDG,IE,R1,R2,R3,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,61) IED4,IE,R1,R2,R3,R4,R5,ANS
               ELSE
                  write(11,62) IEDG,IE,R1,R2,R3,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,62) IED4,IE,R1,R2,R3,R4,R5,ANS
               ENDIF
            elseif (ans.eq.'m') then
               IF (neln.LT.1000) THEN
                  write(11,60) IEDG,IE,R1,R2,R30,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,60) IED4,IE,R1,R2,R31,R4,R5,ANS
               ELSEIF (neln.lt.1000000) then
                  write(11,61) IEDG,IE,R1,R2,R30,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,61) IED4,IE,R1,R2,R31,R4,R5,ANS
               ELSE
                  write(11,62) IEDG,IE,R1,R2,R30,R4,R5,ANS
                  ied4 = iedg+4
                  write(11,62) IED4,IE,R1,R2,R31,R4,R5,ANS
               ENDIF
            endif
   50      CONTINUE
           r30 = r31
           r31 = r31 + dz(ilev+1)
   55    CONTINUE
   60    FORMAT(I3,I3,5G14.6,1X,A1)
   61    FORMAT(I2,I6,5G14.6,1X,A1)
   62    FORMAT(I1,I7,5G14.6,1X,A1)

      ENDIF
      RETURN
C
C     Error handling:
C
  500 CONTINUE
      WRITE(6,501)
  501 FORMAT(2X,'ERROR READING CURVE SIDE DATA'
     $    ,/,2X,'ABORTING IN ROUTINE RDCURVE.')
      CALL EXIT
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION MOD1(I,N)
C
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C
      MOD1=0
      IF (I.EQ.0) THEN
         RETURN
      ENDIF
      IF (N.EQ.0) THEN
         WRITE(6,*) 
     $  'WARNING:  Attempt to take MOD(I,0) in FUNCTION MOD1.'
         RETURN
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine sweep_circ(xc,yc,zc,x,y,n,theta)             ! z0=theta

      real xc(n),yc(n),zc(n),x(n),y(n)

c     project x() and y() into theta plane parallel to z-axis

c     Intitial convention:   x=r, y=z

      one = 1.
      pi  = 4.*atan(one)
      th  = -pi*theta/180.  ! Flip sweep direction to preserve right-handedness

      do i=1,n
         rad   = x(i)
         xc(i) = rad*cos(th)
         yc(i) = rad*sin(th)
         zc(i) = y(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine out_e_hdr(e,ilev,a,igroup,io)
      integer e
      character*1 a

      jlev = ilev
      do k=1,9
         if (jlev.le.99999) goto 1
         jleve = jlev / 10
      enddo
    1 continue

      write(io,10) e,jlev,a,igroup
   10 format(5x,'ELEMENT',i12,' [',i5,a1,']    GROUP',i5)

      return
      end
c-----------------------------------------------------------------------
      subroutine newcurve

      parameter(nelm=9999)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      logical ifflow,ifheat
      common /arraz/ zmin,zmax,dz(nelm)

      common /arral/ ifcirc
      logical ifcirc

      integer e,en,edge,edg1


      neln = nel*nlev


      do ipass=1,2     ! 1st pass, count number of nontrivial curves
       zo = 0.
       z1 = zmin
       do ilev=1,nlev

        z0 = z1
        z1 = z1+dz(ilev)

        do e=1,nel      

          en = e + nel*(ilev-1)

          do edge=1,4
            if (ccurve(edge,e).ne.' ') then
              call get_midside(xm,ym,edge,e)
              call sweep_circ (xc,yc,zc,xm,ym,1,z0)  ! z0=theta
              call outcurve(xc,yc,zc,zo,zo,'m',edge,en,neln,ipass)
            endif
          enddo

          do edge=5,8
            edg1 = edge-4
            if (ccurve(edg1,e).ne.' ') then
              call get_midside(xm,ym,edg1,e)
              call sweep_circ (xc,yc,zc,xm,ym,1,z1)  ! z1=theta
              call outcurve(xc,yc,zc,zo,zo,'m',edge,en,neln,ipass)
            endif
          enddo

          do edge=9,12
            zh = .5*(z0+z1)
            i  = edge-8
            call sweep_circ (xc,yc,zc,x(i,e),y(i,e),1,zh)  ! zh=theta
            call outcurve   (xc,yc,zc,zo,zo,'m',edge,en,neln,ipass)
          enddo

        enddo
       enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outcurve(r1,r2,r3,r4,r5,cc,edge,e,nel,ipass)
      character*1 cc

      integer ncurve,ilast
      save    ncurve,ilast
      data    ncurve,ilast  /0,0/

      integer edge,e

      if (ilast.eq.1.and.ipass.eq.2) then

          write(11,11)
          write(11,12) ncurve
   11     format(' ***** CURVED SIDE DATA *****')
   12     format(i8
     $      ,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')

      endif

      ilast = ipass

      if (ipass.eq.1.and.cc.ne.' ') then

         ncurve = ncurve+1

      elseif (ipass.eq.2.and.cc.ne.' ') then

         if (nel.lt.1000) then

            write(11,60) edge,e,r1,r2,r3,r4,r5,cc
   60       format(i3,i3,5g14.6,1x,a1)

         elseif (nel.LT.1000000) then

            write(11,61) edge,e,r1,r2,r3,r4,r5,cc
   61       format(i2,i6,5g14.6,1x,a1)

         else

            write(11,62) edge,e,r1,r2,r3,r4,r5,cc
   62       format(i1,i7,5g14.6,1x,a1)

         endif

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_midside(xm,ym,edge,e) ! Works only for edge=1,...,4
      integer edge,e

      parameter(nelm=9999)
      common /array/ x(4,nelm),y(4,nelm),bc(5,6,nelm),curve(6,12,nelm)
      common /arrai/ nlev,nel,ncurve,npscal
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6,nelm)
      character*1 ccurve(12,nelm)
      character*1 ca(nelm)
      logical ifflow,ifheat

      real l

      i0 = edge
      i1 = edge+1
      i1 = mod1(i1,4)
      x0 = x(i0,e)
      y0 = y(i0,e)
      x1 = x(i1,e)
      y1 = y(i1,e)

      ym = .5*(y0+y1)  ! Default, straight-sided values
      xm = .5*(x0+x1)

      if (ccurve(edge,e).eq.'m') then

         xm = curve(1,edge,e)
         ym = curve(2,edge,e)

      elseif (ccurve(edge,e).eq.'C') then

         dx = x1-x0
         dy = y1-y0

         l   = dx*dx + dy*dy
         if (l.gt.0) l=sqrt(l)/2.
         rad = curve(1,edge,e)  ! signed radius since (x0,x1) is directed graph

         if (l.ge.abs(rad)) then
            write(6,1) e,edge,rad
    1       format('Error: invalid curve side radius, el/edge:'
     $            ,i9,i3,1pe12.5)
            return
         endif

         delta = 1.-(l/rad)**2
         delta = rad*(1-sqrt(delta)) ! signed offset from chord

         xm    = xm + .5*delta*(y1-y0)/l
         ym    = ym - .5*delta*(x1-x0)/l

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine capit(lettrs,n)  ! capitalizes string of length n
      character lettrs(n)

      do i=1,n
         int=ichar(lettrs(i))
         if(int.ge.97 .and. int.le.122) then
            int=int-32
            lettrs(i)=char(int)
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
