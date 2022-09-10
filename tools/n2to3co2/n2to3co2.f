c----------------------------------------------------------------------
      program n2to3co2
#     include "SIZE"
      character*80 file
      character*1  file1(80)
      equivalence (file1,file)
      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)

      logical ifflow,ifheat,ifcht

      common /arral/ ifcirc
      logical ifcirc

      integer e

      write(6,*)
      write(6,*) 'This is the code that establishes proper'
      write(6,*) 'element connectivities, as of May., 2020.'
      write(6,*)
 
c     for workstation:
      in = 5
 
c     for delta:
c     open(unit=7,file='indat',status='old',err=1999)
c     in = 7
c1999 continue

 
c     Get file name (input)
      write(6,*) 'Input .con/.co2 name to extrude:'      
      call blank(file,80)
      read(in,80) file
      len = ltrunc(file,80)
   80 format(a80)
 
c     Get file name (output)
      write(6,*) 'Input .con/.co2 output name'      
      call blank(fout,80)
      read(in,80) fout
      lou = ltrunc(fout,80)

c     Get output type
      write(6,*) 'Input 0:ASCII or 1:BINARY for output file'
      read(in,*) itype
      if (itype.ne.0.AND.itype.ne.1) call exitti
     $  ('Error: invalid itype ',itype)
 
c     Get nlevels
      nlev = 1
      write(6,*) 
     $'input number of levels: (1, 2, 3,...; < 0 for circular sweep.):'
      read(in,*) nlev

      ifcirc = .false.
      if (nlev.lt.0) ifcirc=.true.
      if (ifcirc) then
        write(6,*) '[TODO] circular mode is not support yet'
        call exitt
      endif
      nlev = abs(nlev)

c     Find input type
      call search_file(file,'.con',ifile)
      if (ifile.gt.0) then
        intype=0
      else 
        call search_file(file,'.co2',ifile)
        if (ifile.gt.0) then
          intype=1
        else
          write(6,'(A)') 'Error no .con / .co2 file found!'
          call exitt
        endif
      endif

c     Read input: 2D con/co2
      call read_2dcon(file,nel,nelv,nvtx2d,intype,ifcht) ! the data is saved into global variables

c     Open output file
      if (itype.eq.0) then
        call chcopy(fout1(lou+1),'.con',4)
        open(unit=11, file=fout)
      else
        call chcopy(fout1(lou+1),'.co2',4)
        call byte_open(fout,ierr)
        if(ierr.gt.0) call exitti('fail to open file',ierr)
      endif

c     Write output 
      call con23(neln,nvtx2d,nvtx3d,itype,ifcht)
      write(6,6) nvtx3d,neln,(fout1(k),k=1,lou+4)
    6 format(i16,' vertices and ',i16,' elements written to ',40a1)

c     Close output files
      if(itype.eq.0) then
        close (unit=11)
      else
        call byte_close(ierr)
        if(ierr.gt.0) call exitti('fail to close file',ierr)
      endif

      stop
      end
c-----------------------------------------------------------------------
      subroutine con23(neln,nvtx2d,nvtx3d,itype,ifcht)
c     input nlev,nel,nvtx2d,itype,icon2d
c     output neln,nvtx3d
#     include "SIZE"

      character*132 hdr
      character*5   version
      real*4 test
      data   test  / 6.54321 /
 
c     Nekton stuff
      logical ifflow,ifheat,ifmhd,ifcht

      common /arral/ ifcirc
      logical ifcirc

      logical ifper
      integer e,iper

      neln = nlev*nel
      nelnv= nlev*nelv
 
      ! Choose BC for Z direction
      write(6,*)'Enter Z boundary condition (1=periodic 0=others)'
      read(5,*) iper
      ifper = .false.
      if (iper.eq.1) ifper = .true.
      if (ifper.and.nlev.lt.3) then
         write(6,*) 'NOTE: nlev < 3 not allowed with periodic bcs'
         write(6,*) 'nlev =',nlev
         write(6,*) 'ABORT'
         stop
      endif


      if(itype.eq.1) then    !co2
        ! write header
        call blank(hdr,132)
        if (neln.lt.10000000)  then
          write(hdr,1) '#v001',neln,nelnv,8
        else   
          write(hdr,2) '#v002',neln,nelnv,8
    2     format(a5,i16,i16,i5) 
        endif   

        write(6,*) 'hdr:', hdr
        call byte_write(hdr,132/4,ierr)
        call byte_write(test,1,ierr) ! write the endian discriminator
        if(ierr.gt.0) call exitti('fail to write header',ierr)
        ! write co2
        call dump_co2(ifper,nvtx2d,nvtx3d,ifcht)
      else
        version = '#v001'
        write(11,1) version,neln,nelnv,8
        write(6,1) version,neln,nelnv,8
        ! write con
        call dump_con(ifper,nvtx2d,nvtx3d,ifcht)
      endif

    1 format(a5,3i12) ! header

      return
      end
c-----------------------------------------------------------------------
      subroutine dump_con(ifper,nvtx2d,nvtx3d,ifcht)
c     input: iper, icon2d, nvtx2d
c     output: nvtx3d
#     include "SIZE"
      
      logical ifper,ifcht
      integer e

      write(6,*)'Z(5)-Z(6) Periodic=',ifper

      ivtx=0
      do i=1,nlev-1
        do j=1,nelv
          e = j + (i-1)*nelv
          write(11,2) e
     $     ,(icon2d(ii,j)+ivtx,ii=1,4)
     $     ,(icon2d(ii,j)+ivtx+nvtx2d,ii=1,4)
        enddo
        ivtx=ivtx+nvtx2d
      enddo

      if (ifper) then
        i=nlev
        do j=1,nelv
          e = j + (i-1)*nelv
          write(11,2) e
     $     ,(icon2d(ii,j)+ivtx,ii=1,4)
     $     ,(icon2d(ii,j),ii=1,4)
        enddo
        ivtx=ivtx+nvtx2d
      else
        i=nlev
        do j=1,nelv
          e = j + (i-1)*nelv
          write(11,2) e
     $     ,(icon2d(ii,j)+ivtx,ii=1,4)
     $     ,(icon2d(ii,j)+ivtx+nvtx2d,ii=1,4)
        enddo
        ivtx=ivtx+nvtx2d*2
      endif

      if(ifcht)then
        ndel=nel-nelv
        ivtx=0
        do i=1,nlev-1
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
            write(11,2) e
     $       ,(icon2d(ii,j2d)+ivtx,ii=1,4)
     $       ,(icon2d(ii,j2d)+ivtx+nvtx2d,ii=1,4)
          enddo
          ivtx=ivtx+nvtx2d
        enddo

        if (ifper) then
          i=nlev
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
            write(11,2) e
     $       ,(icon2d(ii,j2d)+ivtx,ii=1,4)
     $       ,(icon2d(ii,j2d),ii=1,4)
          enddo
          ivtx=ivtx+nvtx2d
        else
          i=nlev
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
            write(11,2) e
     $       ,(icon2d(ii,j2d)+ivtx,ii=1,4)
     $       ,(icon2d(ii,j2d)+ivtx+nvtx2d,ii=1,4)
          enddo
          ivtx=ivtx+nvtx2d*2
        endif

      endif

      nvtx3d=ivtx

    2 format(9i16)
      return
      end
c-----------------------------------------------------------------------
      subroutine dump_co2(ifper,nvtx2d,nvtx3d,ifcht)
c     input: iper, icon2d, nvtx2d
c     output: nvtx3d
#     include "SIZE"

      logical ifper,ifcht
      integer iwrk(8+1),e

      write(6,*)'Z(5)-Z(6) Periodic=',ifper

      ierr=0

      ivtx=0
      do i=1,nlev-1
        do j=1,nelv
          e = j + (i-1)*nelv

          iwrk(1)=e
          do k=1,4
            iwrk(k+1) = icon2d(k,j)+ivtx
            iwrk(k+5) = icon2d(k,j)+ivtx+nvtx2d
          enddo
          call byte_write(iwrk,9,ierr)
        enddo
        ivtx=ivtx+nvtx2d
      enddo

      if (ifper) then
        i=nlev
        do j=1,nelv
          e = j + (i-1)*nelv

          iwrk(1)=e
          do k=1,4
            iwrk(k+1) = icon2d(k,j)+ivtx
            iwrk(k+5) = icon2d(k,j)
          enddo
          call byte_write(iwrk,9,ierr)
        enddo
        ivtx=ivtx+nvtx2d

      else
        i=nlev
        do j=1,nelv
          e = j + (i-1)*nelv

          iwrk(1)=e
          do k=1,4
            iwrk(k+1) = icon2d(k,j)+ivtx
            iwrk(k+5) = icon2d(k,j)+ivtx+nvtx2d
          enddo
          call byte_write(iwrk,9,ierr)
        enddo
        ivtx=ivtx+nvtx2d*2
      endif


      if(ifcht)then

        ivtx=0
        ndel=nel-nelv
        do i=1,nlev-1
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
  
            iwrk(1)=e
            do k=1,4
              iwrk(k+1) = icon2d(k,j2d)+ivtx
              iwrk(k+5) = icon2d(k,j2d)+ivtx+nvtx2d
            enddo
            call byte_write(iwrk,9,ierr)
          enddo
          ivtx=ivtx+nvtx2d
        enddo
  
        if (ifper) then
          i=nlev
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
  
            iwrk(1)=e
            do k=1,4
              iwrk(k+1) = icon2d(k,j2d)+ivtx
              iwrk(k+5) = icon2d(k,j2d)
            enddo
            call byte_write(iwrk,9,ierr)
          enddo
          ivtx=ivtx+nvtx2d
  
        else
          i=nlev
          do j=1,ndel
            e = j + (i-1)*ndel + nelv*nlev
            j2d=j+nelv
  
            iwrk(1)=e
            do k=1,4
              iwrk(k+1) = icon2d(k,j2d)+ivtx
              iwrk(k+5) = icon2d(k,j2d)+ivtx+nvtx2d
            enddo
            call byte_write(iwrk,9,ierr)
          enddo
          ivtx=ivtx+nvtx2d*2
        endif

      endif

      if(ierr.gt.0) call exitti('write error',ierr)

      nvtx3d=ivtx

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
      integer function indx1(s1,s2,l2)
      character*80 s1,s2

      n1=80-l2+1
      indx1=0
      if (n1.lt.1) return

      do 300 i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx1=i
            return
         endif
300   continue

      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      integer a(1), b(1)
      do i = 1,n
         a(i) = b(i)
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
      subroutine izero(x,n)
      integer x(1)
      do i=1,n
         x(i) = 0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(name,ie)
      character*40 name
      write(6,*) name, ie
      stop
      end
c-----------------------------------------------------------------------
      subroutine exitt
      call exit
      return
      end
c-----------------------------------------------------------------------
      subroutine exitrr(stringi,r1,r2)
      character*1 stringi(132)
      character*1 stringo(132)
      character*26 s26

      call blank  (stringo,132)
      call chcopy (stringo,stringi,132)
      len = indx1 (stringo,'$',1)
      write(s26,26) r1,r2
   26 format(1p2e13.4)
      call chcopy(stringo(len),s26,26)

      if (nid.eq.0) write(6,1) (stringo(k),k=1,len+25)
      if (nid.eq.0) write(6,*)
    1 format(/,'EXIT: ',132a1)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      integer function ivlmax(vec,n)
      integer vec(1),tmax
      if (n.eq.0) then
         ivlmax=0
         return
      endif
      TMAX =-8888888
      do i=1,n
         TMAX = MAX(TMAX,VEC(I))
      enddo
      Ivlmax = tmax
      return
      end
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)
c
      real*4 bytetest,test2
      real*4 test_pattern
      save   test_pattern
c
      test_pattern = 6.54321
      eps          = 0.00020
      etest        = abs(test_pattern-bytetest)
      if_byte_swap_test = .true.
      if (etest.le.eps) if_byte_swap_test = .false.

      ierr  = 0
      test2 = bytetest
      call byte_reverse(test2,1,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error with byte_reverse in if_byte_swap_test ',ierr)
c     write(6,*) 'Byte swap:',if_byte_swap_test,bytetest,test2

      return
      end
c-----------------------------------------------------------------------
      subroutine isort(a,ind,n)
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
         GOTO 200
         endif
           a(i) = aa
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine chk_vtx(con,nel,nvtx,ierr)
c     This will check uniq(con) = 1:nvtx_max, so max(vtx)=nvtx
      integer con(nel*4),ind(nel*4)
      integer iwrk(nel*4)

      ierr=10

      call icopy(iwrk,con,nel*4)
      call isort(iwrk,ind,nel*4)

      inow=1
      do i=1,nel*4
        if (iwrk(i).gt.inow) then
          if(iwrk(i).gt.inow+1) then
            ierr=1
            write(*,*)'chk_vtx fail, missing uniq vtx',inow+1
            return
          else
            inow=inow+1
          endif
        endif
        nvtx=inow
      enddo

      if(nvtx.ne.iwrk(nel*4)) then
        ierr=2
        write(6,*)'chk_vtx is not passed',nvtx,iwrk(nel*4)
        return
      endif

      ierr=0
      return
      end
c-----------------------------------------------------------------------
      subroutine search_file(fname,fext,ierr)
c
      logical ifexist
c
      character*1  fname(1)
      character*80 file
      character*1  file1(80)
      equivalence (file1,file)
      character*4  fext

      ifexist = .false.

c     Get file name
      len = ltrunc(fname,80)
      call chcopy(file1(1),fname,80)
      call chcopy(file1(len+1),fext,4)
      inquire(file=file, exist=ifexist)

      ierr=0
      if(ifexist) then
        ierr=1
        write(6,*)'Find file: ',file
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_2dco2(fname,nelgti,nelgvi,ifcht)
#     include "SIZE"

      character*1  fname(1)
      character*80 file
      character*1  file1(80)
      equivalence (file1,file)

      character*132 hdr
      character*5 version
      logical ifbswap,if_byte_swap_test
      real*4  test
      integer iwrk(1+8),nelgti,nelgvi
      logical ifcht

      ierr=0

      ! Open co2
      len = ltrunc(fname,80)
      call chcopy(file1(1),fname,80)
      call chcopy(file1(len+1),'.co2',4)
      write(6,*) 'reading ', file
      call byte_open(file,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error opening file in open_bin_file ',ierr)

      ! Read header
      call byte_read(hdr,sizeof(hdr)/4,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error reading header in open_bin_file ',ierr)
      read (hdr,*) version,nelgti,nelgvi,nv
      write(6,1) version,nelgti,nelgvi,nv

      ! Check before reading more
      wdsizi=4
      if(version.eq.'#v001') wdsizi=8
      if(nelgti.gt.nelxym) call exitti
     $  ('ABORT: input mesh is too large ',nelxym)
      if(nv.ne.4) call exitti
     $  ('ABORT: input co2 is not 2D ',nv)

      ifcht=.false.
      if(nelgti.ne.nelgvi) ifcht=.true.

      call byte_read(test,1,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error reading test number in open_bin_file ',ierr)
      ifbswap = if_byte_swap_test(test)

      ! Read connectivity
      do e=1,nelgti 
        call byte_read(iwrk,nv+1,ierr)
        call icopy(icon2d(1,e),iwrk(2),nv)
      enddo

      ! Close
      call byte_close(ierr)
      if(ierr.gt.0) call exitti('closing file',ierr)

    1 format(a5,3i12) ! header
      return
      end
c-----------------------------------------------------------------------
      subroutine read_2dcon(fname,nelgti,nelgvi,nvtx2d,intype,ifcht)

#     include "SIZE"
      character*1  fname(1)
      character*80 file
      character*1  file1(80)
      equivalence (file1,file)

      character*5 version
      integer wdsizi
      integer iwrk(1+8),nelgti,nelgvi,nvtx2d
      logical ifcht

      if(intype.eq.0) then
        ! Open con
        len = ltrunc(fname,80)
        call chcopy(file1(1),fname,80)
        call chcopy(file1(len+1),'.con',4)
        write(6,*) 'reading ', file
        open (unit=29,file=file)
        read(29,1) version,nelgti,nelgvi,nv
        write(6,1) version,nelgti,nelgvi,nv

        ! Check before reading more
        wdsizi=4
        if(version.eq.'#v001')wdsizi=8
        if(nelgti.gt.nelxym) call exitti
     $  ('ABORT:  increase nelxym in con2to3/SIZE ',nelxym)
        if(nv.ne.4) call exitti
     $  ('ABORT:  input con is not 2D ',nv)

        ifcht=.false.
        if(nelgti.ne.nelgvi) ifcht=.true.

        ! Read connectivity
        do e=1,nelgti
          read(29,2)(iwrk(ii),ii=1,nv+1)
          call icopy(icon2d(1,e),iwrk(2),nv)
        enddo
        close(unit=29)

      else

        call read_2dco2(fname,nelgti,nelgvi,ifcht)

      endif

      ! Get nvtx
      ierr=0
      call chk_vtx(icon2d,nelgti,nvtx2d,ierr)
      if(ierr.gt.0) call exitti
     $  ('ABORT: chk_vtx failed',ierr)
      write(6,*)'found',nvtx2d,' vertices in 2d' 

    1 format(a5,3i12) ! header
    2 format(5i12)
      return
      end
c-----------------------------------------------------------------------

