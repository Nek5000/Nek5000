c=======================================================================
      subroutine get_fld_data(fname,ierr)

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MISC'

      integer e
      integer*8 ioff
      integer nfld 

      CHARACTER*20 FRMAT
      CHARACTER  EXCODE(10)
      CHARACTER*80 FNAME
      CHARACTER*1  FNAME1(80)

      logical      if_byte_sw

      INTEGER      HNAMI (20),len
      CHARACTER*80 HNAME
      CHARACTER*1  HNAME1(80)
      EQUIVALENCE (HNAME,HNAME1)
      EQUIVALENCE (HNAME,HNAMI )

      real*4 tdump(lx1*ly1*lz1,20)
      integer ist
      data ist /0/
c
      ierr = 1

      len= ltrunc(fname,79)
      call izero (hnami,20)
      call chcopy(hname1,fname,len)
c           test for presence of file
      open (unit=91,file=hname
     &      ,form='unformatted',status='old',err=500)
      close(UNIT=91)
      call byte_open(hname)
c      call byte_open('test.fld' // CHAR(0))



      call get_hdr     (nel,nx,ny,nz,nfld,time,excode,if_byte_sw)
      nx1  = nx
      ny1  = ny
      nz1  = nz
      nelv = nel
      nelt = nel

      write(*,'A,A') 'Reading fld file ', fname
      write(*,'I6,X,3I3') nel,nx,ny,nz

      if (nx1.gt.lx1.or.nz1.gt.lz1.or.nelv.gt.lelv) then
         write(6,*) 'Increase lx1,ly1,lz1,lelv,lelt in SIZE'
         write(6,*) 'Old:',lx1,ly1,lz1,lelv,lelt
         write(6,*) 'New:',nx1,ny1,nz1,nelv,nelt
         call exitt
      endif

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
      write(*,*) 'EXCODER=', (excode(j),j=1,10)
      write(*,*) 'scaning EXCODER ...', excode(j)
      nouts = 0
      do j = 1,10
            if(excode(j).eq.'X') then
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
            if (excode(j).eq.'U') then
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
            if (excode(j).eq.'P') then
               nouts=nouts+1
               iposp=nouts
               write(6,3001) 'found P', iposp
            endif
            if (excode(j).eq.'T') then
               nouts=nouts+1
               ipost=nouts
               write(6,3001) 'found T', ipost
            endif
            if (excode(j).eq.'S') then
               READ(EXCODE(j+1),'(I1)') NPS1
               READ(EXCODE(j+2),'(I1)') NPS0
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
C     Read fld File
      ioff = 1
      do 101 e=1,nel
             if (mod(e,2000).eq.0) 
     &          write(6,*)'Reading element ',e,' out of ',nel
             do ii=1,nfld
                call byte_read(tdump(1,II),nxyz)
             enddo
 
             IF(IPOSX.GE.1) CALL COPY4r(XM1(ioff),TDUMP(1,IPOSX),nxyz)
             IF(IPOSY.GE.1) CALL COPY4r(YM1(ioff),TDUMP(1,IPOSY),nxyz)
             IF(IPOSZ.GE.1) CALL COPY4r(ZM1(ioff),TDUMP(1,IPOSZ),nxyz)

             IF(IPOSU.GE.1) CALL COPY4r(VX (ioff),TDUMP(1,IPOSU),nxyz)
             IF(IPOSV.GE.1) CALL COPY4r(VY (ioff),TDUMP(1,IPOSV),nxyz)
             IF(IPOSW.GE.1) CALL COPY4r(VZ (ioff),TDUMP(1,IPOSW),nxyz)

             IF(IPOSP.GE.1) CALL COPY4r(PR (ioff),TDUMP(1,IPOSP),nxyz)
             IF(IPOST.GE.1) CALL COPY4r(T(ioff,1),TDUMP(1,IPOST),nxyz)
             do is = 1, ist
                IF(IPOSPS(is).GE.1) 
     &             CALL COPY4r(T(ioff,is+1),TDUMP(1,IPOSPS(is)),nxyz)
             enddo
             ioff = ioff + nxyz
  101 continue

      ierr = 0
      call byte_close
      return

 500  write(*,'(A)') 'ERROR: cannot read fld file!'
      ierr =1

 3001 format(A,I2)

      end

c=======================================================================

      subroutine get_hdr(nel,nx,ny,nz,nfld,time,excode,if_byte_sw)

      character*10 excode

      common /c80/ s80
      character*80 s80
      character*1  s81(80)
      equivalence (s80,s81)

      logical if_byte_sw, if_byte_swap_test, if3d
      common /byte_key/ bytetest

      call byte_read(s80,20) !  get 80 character header
      read(s80,*,end=101)
     $    nel,nx,ny,nz,time,istep,excode
 101  continue

      call byte_read(bytetest,1)
      if_byte_sw = if_byte_swap_test(bytetest)

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

