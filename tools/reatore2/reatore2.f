c-----------------------------------------------------------------------
      program trans

c     this program will take a nek all-ascii rea file and
c     extract geometry, bcs, and curve data to a new binary re2 file, plus
c     an ascii rea file for just the parameters

      
      character*4 ans
      real*4      buf(10)

      character*80 file,fout,fbout,string
      character*1  file1(80),fout1(80),fbout1(80),string1(80)
      equivalence (file1,file)
      equivalence (fout1,fout)
      equivalence (fbout1,fbout)
      equivalence (string,string1)

      parameter(lelt=900000)
      common /array/ x(8),y(8),z(8),bc(5,6,lelt),curve(6,8)
      common /arrai/ nel,ncurve
      common /arrac/ cbc
      character*3 cbc(6,lelt)
      integer e,f
      character*1 ccurve(4)
      logical ifflow,ifheat

      character*80 hdr
      real*4 test
      data   test  / 6.54321 /

      write(6,*) 'Input old (source) file name:'      
      call blank(file,80)
      read(5,80) file
      len = ltrunc(file,80)
   80 format(a80)
   81 format(80a1)

      write(6,*) 'Input new (output) file name:'      
      call blank(fout,80)
      read(5,80) fout
      fbout = fout
      lou = ltrunc(fout,80)

      call chcopy(file1(len+1),'.rea',4)
      call chcopy(fout1(lou+1),'.rea',4)
      call chcopy(fbout1(lou+1),'.re2\0',5)

      open(unit=10, file=file)
      open(unit=11, file=fout)
      call byte_open(fbout)

      call scanout(string,'MESH DATA',9,10,11)
      call lineout(11,string,80)

      read (10,*)  nel, ndim, nelv
      nels = -nel
      write(11,11) nels, ndim, nelv
   11 format(i9,i3,i9,1x,'NEL,NDIM,NELV')


      call blank(hdr,80)
      write(hdr,1) nel,ndim,nelv
    1 format('#v001',i9,i3,i9,' this is the hdr')
      call byte_write(hdr,20)   ! assumes byte_open() already issued
      call byte_write(test,1)   ! write the endian discriminator

      write(6,*) nel, ndim, nelv
      do ie=1,nel      
         if (mod(ie,10000).eq.0) write(6,*) ie,nel,' mesh'
csk         read (10,12) igroup
csk   12    format(43x,i5)
         read(10,*) 
         igroup = 0
         call byte_write(igroup, 1)

         if (ndim.eq.3) then
            read (10,*)   (x(k),k=1,4)
            read (10,*)   (y(k),k=1,4)
            read (10,*)   (z(k),k=1,4)
            read (10,*)   (x(k),k=5,8)
            read (10,*)   (y(k),k=5,8)
            read (10,*)   (z(k),k=5,8)

            call byte_write(x,8)
            call byte_write(y,8)
            call byte_write(z,8)
         else
            read (10,*)   (x(k),k=1,4)
            read (10,*)   (y(k),k=1,4)
            call byte_write(x,4)
            call byte_write(y,4)
         endif
      enddo

      read(10,*)
      read(10,*) ncurve
      call byte_write(ncurve,1)
      call blank(ccurve,4)
      write(*,*) ncurve
      if (ncurve.gt.0) then
         do icurve=1,ncurve
            if (mod(icurve,10000).eq.0) write(6,*) icurve, 
     &                                  ncurve, ' curve'
            if (nel.lt.1000) then
               read(10,60) f,e,(buf(k),k=1,5),ccurve(1)
            elseif (nel.lt.1000000) then
               read(10,61) f,e,(buf(k),k=1,5),ccurve(1)
            else
               read(10,62) f,e,(buf(k),k=1,5),ccurve(1)
            endif
            call byte_write(e     ,1)
            call byte_write(f     ,1)
            call byte_write(buf   ,5)
            call byte_write(ccurve,1)
         enddo
   60    format(i3,i3,5g14.6,1x,a1)
   61    format(i2,i6,5g14.6,1x,a1)
   62    FORMAT(I1,i7,5G14.6,1X,A1)
      endif
         
      read (10,80) string ! ***** BOUNDARY CONDITIONS *****

      nface = 2*ndim
      do kpass = 1,50   ! for now, at most 50 fields
         read (10,80) string
         if (indx1(string,'BOUN',4).ne.0) then ! we might have bcs
            if (indx1(string,'NO ',3).eq.0)  then ! we have bcs, read and count
               nbc = 0
               do e=1,nel
               do f=1,nface
                  if (nel.lt.1000) then
                     read(10,20) cbc(f,e),(bc(k,f,e),k=1,5)
                  elseif (nel.lt.100000) then
                     read(10,21) cbc(f,e),(bc(k,f,e),k=1,5)
                  else
                     read(10,22) cbc(f,e),(bc(k,f,e),k=1,5)
                  endif
                  if (cbc(f,e).ne.'E  ') nbc = nbc+1
               enddo
               enddo
   20          FORMAT(1x,A3,6x,5G14.7)
   21          FORMAT(1x,A3,6x,5G14.7)
   22          FORMAT(1x,A3,7x,5G14.7)

               write(6,*) kpass,nbc,' Number of bcs'
               call byte_write(nbc,1)

               do e=1,nel
               do f=1,nface
                  if (cbc(f,e).ne.'E  ') then
                     call icopy      (buf(1),e,1)
                     call icopy      (buf(2),f,1)
                     call copy       (buf(3),bc(1,f,e),5)
                     call blank      (buf(8),4)
                     call chcopy     (buf(8),cbc(f,e),3)
                     call byte_write (buf,8)
                  endif
               enddo
               enddo

            endif
         else
            goto 51
         endif

      enddo
   51 continue

      len = ltrunc(string,80)
      write(11,81) (string1(k),k=1,len)
c     write(6 ,81) (string1(k),k=1,len)

      call scanout(string,'end',3,10,11)

      close (unit=10)
      close (unit=11)
      call byte_close (fbout)

  999 continue
      stop
 9999 continue
      write(6,*) 'could not find file "indat". abort.'
      stop
      end
c-----------------------------------------------------------------------
