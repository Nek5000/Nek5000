c=======================================================================
      program convert

      integer e
      integer nfld
      integer maxps,nxmax
      parameter(nout_max=30,nmax=20)
       

      CHARACTER*1  EXCODE1(10)
      CHARACTER*10 EXCODE
      EQUIVALENCE (EXCODE,EXCODE1)


      CHARACTER*80 FNAME
      CHARACTER*1  FNAME1(80)
      EQUIVALENCE (FNAME,FNAME1)
   

      logical      if_byte_sw

      INTEGER      HNAMI (79),len
      CHARACTER*80 HNAME
      CHARACTER*1  HNAME1(80)
      EQUIVALENCE (HNAME,HNAME1)
      EQUIVALENCE (HNAME,HNAMI )

      real*4 tdump(nmax*nmax*nmax,nout_max)
      real time
      integer istep,ist,iposps(nout_max)
      data ist /0/

      write(*,*) 'Enter the fld binary filename to convert:'
      read(*,*) fname
c
      ierr = 1
      len= ltrunc(fname,79)
      call izero (hnami,79)
      call chcopy(hname1,fname,len)
c           test for presence of file
      open (unit=91,file=hname
     &      ,form='unformatted',status='old',err=500)
      close(UNIT=91)
      call byte_open(hname)

      call get_hdr     (nel,nx,ny,nz,nfld,time,istep,excode1,if_byte_sw)
      nx1  = nx
      ny1  = ny
      nz1  = nz
      nelv = nel
      nelt = nel

      write(*,*) 'Reading fld file ', fname
      write(*,'(I6,X,3I3)') nel,nx,ny,nz
      ndim = 2
      if (nz1.gt.1) ndim = 3

      if3d = .false.
      if (nz1.gt.1) if3d = .true.
      write(*,*) 'if3d', if3d

      iposx = -1
      iposy = -1
      iposz = -1
      iposu = -1
      iposv = -1
      iposw = -1
      iposp = -1
      ipost = -1
      do j = 1, maxps
        iposps(j) = -1
      enddo

C     Figure out position of data in file
      write(*,*) 'scaning EXCODER ...'
      write(*,*) 'readed EXCODER=', excode
      nouts = 0
      do j = 1,10
            if(excode1(j).eq.'X') then
              nouts=nouts+1
              iposx=nouts
              nouts=nouts+1
              iposy=nouts
              write(6,3001) 'found X', iposx
              if(id3d) then
                nouts=nouts+1
                iposz=nouts
              endif
            endif
            if (excode1(j).eq.'U') then
               nouts=nouts+1
               iposu=nouts
               nouts=nouts+1
               iposv=nouts
               write(6,3001) 'found U', iposu
               if (if3d) then
                  nouts=nouts+1
                  iposw=nouts
               endif
            endif
            if (excode1(j).eq.'P') then
               nouts=nouts+1
               iposp=nouts
               write(6,3001) 'found P', iposp
            endif
            if (excode1(j).eq.'T') then
               nouts=nouts+1
               ipost=nouts
               write(6,3001) 'found T', ipost
            endif
            if (excode1(j).eq.'S') then
               READ(EXCODE1(j+1),'(I1)') NPS1
               READ(EXCODE1(j+2),'(I1)') NPS0
               NPS=10*NPS1 + NPS0
               DO IS = 1, NPS
                  NOUTS=NOUTS + 1
                  iposps(IS)=NOUTS
               ENDDO
               write(6,3001) 'npscal', nps
               goto 51
            endif
      enddo

  51  nfld=nouts
      write(*,'(A,I3)') 'total number of fields: ', nfld
      if(nfld.eq.0) then
        write(*,*) 'ERROR: no fields found!'
        call exitt
      endif

      nxyz = nx*ny*nz

      fname = 'asc_' // hname
      len= ltrunc(fname,79)
      call izero (hnami,79)
      call chcopy(hname1,fname,len)

      open(55,file=hname)
 
      ! write header to ASCII file
      call blank(fname,80)
      write(fname,'(i10,1x,i2,1x,i2,1x,i2,1x,1p1e18.9,1x,i9,1x,a)') 
     &      nel,nx,ny,nz,time,istep,excode
      write(55,'(A)') fname


C     Read fld File
      do 101 e=1,nel
             if (mod(e,2000).eq.0) 
     &          write(6,*)'Converting element ',e,' out of ',nel
             do ii=1,nfld
                call byte_read(tdump(1,II),nxyz)
                if(if_byte_sw) call byte_reverse(tdump(1,II),nxyz)
             enddo
             ! write data to ASCII file         
             do m=1,nxyz 
                write(55,'(1p50E15.7)') 
     &               (tdump(m,ii),ii=1,nfld)
             enddo
  101 continue

      ierr = 0
      close(55)
      call byte_close

      write(*,*)
      write(*,*) hname
      write(*,*) 'ASCII file written'

 500  write(*,'(A)') 'ERROR: cannot read fld file!'
      ierr =1

 3001 format(A,I2)

      end

c=======================================================================

      subroutine get_hdr(nel,nx,ny,nz,nfld,time,istep,excode,if_byte_sw)

      character*10 excode

      character*80 s80
      character*1  s81(80)
      equivalence (s80,s81)

      logical if_byte_sw, if_byte_swap_test

      logical if3d

      real bytetest

      call byte_read(s80,20) !  read 80 character header
      read(s80,*,end=101)
     $    nel,nx,ny,nz,time,istep,excode
 101  continue

      ! small or big endian?
      call byte_read(bytetest,1)
      if(bytetest .gt. 6.5 .and. bytetest .lt. 6.6) then
        if_byte_sw = .false.
      else
        write(*,'(A)') 'Data format seems to be big endian'
        if_byte_sw = .true.
      endif

      return
      end

c=======================================================================
      subroutine copy4r(a,b,n)
      real   a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end

