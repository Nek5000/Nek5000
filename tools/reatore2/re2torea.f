c-----------------------------------------------------------------------
      program re2torea

c     this program will take a nek all-ascii rea file and
c     extract geometry, bcs, and curve data to a new binary re2 file, plus
c     an ascii rea file for just the parameters

      
      character*80 fin
      character*1  fin1(80)
      equivalence (fin1,fin)

      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)

      character*80 fbin
      character*1  fbin1(80)
      equivalence (fbin1,fbin)

      character*80 string
      character*1  string1(80)
      equivalence (string,string1)
c
      character*80 sstring

      parameter(nelm=8888)
      common /array/ x(8),y(8),z(8),bc(5,6),curve(6,8)
      common /arrai/ nel,ncurve
      common /arrac/ cbc,ccurve,ca
      character*3 cbc(6)
      integer id(6), jd(6)
      character*1 ccurve(8,nelm)
      character*1 ca(nelm)
      logical ifflow,ifheat

      character*80 hdr
      real*4 test
      data   test  / 6.54321 /


      
      integer ignor,igroup
      character cignor
c     
c
      call rzero(x,8)
c
c
c     for workstation:
      in = 5
c
c     for delta:
c     open(unit=7,fin='indat',status='old',err=1999)
c     in = 7
c     1999 continue
c
c     Get file name
      write(6,*) 'Input old (source) file name:'      
      call blank(fin,80)
      read(in,80) fin
      fin = 'test'
      len = ltrunc(fin,80)
   80 format(a80)
c
c     Get file name
      write(6,*) 'Input new (output) file name:'      
      call blank(fout,80)
      read(in,80) fout
      fout = 'back'
      fbin = fin
      lou = ltrunc(fout,80)

      call chcopy(fin1(len+1),'.rea',4)
      call chcopy(fout1(lou+1),'.rea',4)
      call chcopy(fbin1(len+1),'.re2\0',5)
c
      print *, fin, fbin, fout
      open(unit=10, file=fin)
      open(unit=11, file=fout)
      call byte_open(fbin)


      call scanout(sstring,'PRESOLVE',8,10,11)

      write(11,10)
   10 format(1x,
     $'**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8')

      call blank(hdr,80)
      call byte_read(hdr,20)
      read (hdr,1) nelv,ndim,nel

 1    format('#v001',i9,i3,i9,' this is the hdr')

      write(11,11) nel, ndim, nelv
 11   format(3i6,11x,'NEL,NDIM,NELV')



      call byte_read(test,1)

      print *, nel, ndim, nelv

      do ie=1,nel      

         call byte_read(igroup, 1)
         print *, ie, igroup
         do ixyz = 1,8
            call byte_read(x(ixyz),1)
            call byte_read(y(ixyz),1)
            call byte_read(z(ixyz),1)
         enddo

         write (11,12) ie, ie, 'a', igroup
         write (11,*)   (x(k),k=1,4)
         write (11,*)   (y(k),k=1,4)
         write (11,*)   (z(k),k=1,4)

         write (11,*)   (x(k),k=5,8)
         write (11,*)   (y(k),k=5,8)
         write (11,*)   (z(k),k=5,8)

      enddo

 12   format(
     $     '            ELEMENT',i5,' [',i5,a1,']    GROUP     ',i1)
      
      call byte_read(ncurve,1)

      write(11,13)
      write(11,14) ncurve
 13   format(' ***** CURVED SIDE DATA *****')
 14   format(i8
     $     ,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')

      IF (NCURVE.GT.0) THEN
         DO 50 ICURVE=1,NCURVE
            call byte_read(iedg, 1)
            print *, iedg, 'iedg'
            call byte_read(ieg, 1)
            call byte_read(r1, 1)
            call byte_read(r2, 1)
            call byte_read(r3, 1)
            call byte_read(r4, 1)
            call byte_read(r5, 1)
            call byte_read(ans, 1)

            IF (nel.LT.1000) THEN
               write(11,60) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ELSE
               write(11,61) IEDG,IEG,R1,R2,R3,R4,R5,ANS
            ENDIF
 50      CONTINUE

 60      FORMAT(I3,I3,5G14.6,1X,A1)
 61      FORMAT(I2,I6,5G14.6,1X,A1)
         
      endif


c     boundary conditions

      call byte_read(string,20) 
      write(11,80) string
c     fluid boundary conditions?

      call byte_read(string,20)
      write (11,80) string

      if (string1(9).eq. 'F') then
         do ie = 1,nel
            do k = 1,6
               call byte_read(cbc(k), 1)
               call byte_read(id(k), 1)
               call byte_read(jd(k), 1)
               do j = 1, 5
                  call byte_read(bc(j,k), 1)
               enddo
            enddo

            if (nel.lt.1000) then
               do  k = 1,6
                  write (11,20) cbc(k),id(k),jd(k),(bc(j,k),j=1,5)
               enddo
            else
               do  k = 1,6
                  write (11,21) cbc(k),id(k),jd(k),(bc(j,k),j=1,5)
               enddo
            endif
               

         enddo
      endif


      call byte_read(string,20)
      write (11,80) string
c     thermal boundary conditions?

      if (string1(9).eq. 'T') then

         do ie = 1,nel
            do k = 1,6
               call byte_read(cbc(k), 1)
               call byte_read(id(k), 1)
               call byte_read(jd(k), 1)
               do j = 1, 5
                  call byte_read(bc(j,k), 1)
               enddo
            enddo
            if (nel.lt.1000) then
               do  k = 1,6
                  write (11,20) cbc(k),id(k),jd(k),(bc(j,k),j=1,5)
               enddo
            else
               do  k = 1,6
                  write (11,21) cbc(k),id(k),jd(k),(bc(j,k),j=1,5)
               enddo
            endif
               

         enddo
      else
         print *, 'no thermal boundary conditions'
      endif
   20 FORMAT(1x,A3,2I3,5G14.6)
   21 FORMAT(1x,A3,i5,i1,5G14.6)


      write(11,80) sstring      
      call scanout(string,'endendend',9,10,11)
c
cc      write(6,*)
cc      write(6,6) neln,(fout1(k),k=1,len+4)
cc    6 format(i4,' elements written to ',40a1)
c
      close (unit=10)
      close (unit=11)
      call byte_close (fbin)

  999 continue
      stop
 9999 continue
      write(6,*) 'could not find file "indat". abort.'
      stop
      end
c-----------------------------------------------------------------------
