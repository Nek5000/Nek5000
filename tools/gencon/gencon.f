      program gencon

#     include "SIZE"

      parameter(lpts=8*lelm)
      integer cell
      common /carrayi/ cell (lpts) 

      real w1,w2,w3,w4
      common /carrayw/ w1   (lpts) , w2   (lpts)
     $               , w3   (lpts) , w4   (lpts)

      real        w14(4*lpts)
      equivalence (w1,w14)

      common /arrayr/  dx(4*lpts)

      common /carrayr/ bc(5*6*lelm)
      real*8 bc
      common /carrayc/ cbc(6,lelm)
      character*3      cbc

      logical ifconn,is_connected,face_conn

      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) then
         wdsize=8
      endif
 
      call makemesh  (cell,nelv,nelt,irnk,dx,cbc,bc,ndim,w14)
c                                    irnk is # unique points

      nfc = 2*ndim
      nv  = 2**ndim

      nic = lpts
      njc = lpts
      if (ndim.eq.3) then
         face_conn = .false.
         do i =1,5
            if (.not.face_conn) then
               call face_chk(face_conn,cell,nv,nelt,nic,njc,w1,w2,w3)
            else
               goto 15
            endif
         enddo 
         write(6,*) "WARNING:Missing Face Connection Not Resolved"
         call exitt(1) 
      endif
  15  continue

      call periodic_vtx(cell,nv,nelt,irnk,dx,ndim,cbc,bc,nfc,w14,w1)
     
c      npts = nv*nelt
c      call iranku       (cell,nrnk,npts,w1)
c      call self_chk     (cell,nv,nelt,1)    ! check for not self-ptg.

      call dmpfile(cell,nv,nelt,nelv,irnk)

      end
c-----------------------------------------------------------------------
      subroutine makemesh(cell,nelv,nelt,irnk,dx,cbc,bc,ndim,wk)

c     read nekton .rea file and make a mesh

#     include "SIZE" 

      integer      cell(1)
      character*3  cbc (6,1)
      real*8       bc  (5,6,1)
      real         dx  (1)
      real         wk  (1)
      integer e,f

      character*3 cbt(6)
      real*8      bt(5,6)

      parameter(lpts=8*lelm)

      common /arrayi/ i_n(lpts) , j_n(4*lpts)
     $                          , j_o(4*lpts)

      logical ifbinary,ifbswap
      integer buf(30)

      integer eface(6)  ! return Nekton preprocessor face ordering
      save    eface
      data    eface / 4 , 2 , 1 , 3 , 5 , 6 /
         
      io = 10
      ifco2 = .false.

      call getreafile('Input .rea / .re2 name:$',ifbinary,io,ierr)
      if (ierr.gt.0) then 
         write(6,'(A)') 'Error no .rea / .re2 file found!'
         call exitt(1) 
      endif

      write(6,'(A)') 'Input mesh tolerance (default 0.2):'
      write(6,'(A,A)') 'NOTE: smaller is better, but generous is more ',
     &                 'forgiving for bad meshes.'

      read(5,'(f7.2)') qin
      if(qin.gt.0) then
        q = qin
      else
        write(6,'(A,2f7.2)') ' using default value'
        q = 0.2
      endif
  
      call cscan_dxyz (dx,nelt,nelv,ndim,ifbinary,ifbswap)

      ierr = 0
      if (ifbinary) then
         ifco2 = .true.
         ! skip curved side data
         if(wdsizi.eq.8) then 
           call byte_read(rcurve,2,ierr)
           if (ifbswap) call byte_reverse8(rcurve,2,ierr)
           ncurve = rcurve
           if(ierr.ne.0) call exitti
     $         ('Error reading ncurve in makemesh ',ierr)
           do k = 1,ncurve
              call byte_read(buf,16,ierr)
              if(ierr.ne.0) call exitti
     $         ('Error reading curve data in makemesh ',ierr)
           enddo
         else
           call byte_read(ncurve,1,ierr)
           if (ifbswap) call byte_reverse(ncurve,1,ierr)
           if(ierr.ne.0) call exitti
     $         ('Error reading ncurve in makemesh ',ierr)
           do k = 1,ncurve
              call byte_read(buf,8,ierr)
              if(ierr.ne.0) call exitti
     $         ('Error reading curve data in makemesh ',ierr)
           enddo
         endif
c        For current version, only need the fluid bcs.
c        Also, if we start to have fluid/thermal domains with differing
c        periodic bc topologies, then we'll need something completely
c        different (two distinct numberings).
c
c        In fact, we need the thermal bcs because of periodicity... 
c        This is not pretty --- there are many cases to be considered.
c
c        For now, we default to solid topology if nelt > nelv

         call rd_bc_bin(cbc,bc,nelv,nelt,ifbswap)

         call byte_close(ierr)
         if(ierr.ne.0) call exitti
     $       ('Error closing file in makemesh ',ierr)
      else
         call cscan_bcs   (cbc,bc,nelv,nelt,ndim)
      endif
c     call outbc(cbc,bc,nelt,ndim,' CBC 1')
      close (unit=10)  ! close .rea file

      nface = 2*ndim
      do e=1,nelt !  SWAP TO PREPROCESSOR NOTATION
         call chcopy(cbt,cbc(1,e)  ,nface*3)
         call copy8 ( bt, bc(1,1,e),nface*5)
         do f=1,nface
            call copy8 ( bc(1,f,e), bt(1,eface(f)),5)
            call chcopy(cbc(  f,e),cbt(  eface(f)),3)
         enddo
      enddo
c     call outbc(cbc,bc,nelt,ndim,' CBC 2')


c     Compress vertices based on coordinates
      call unique_vertex2(cell,dx,ndim,nelt,q,i_n,j_n,j_o,wk)

      nv   = 2**ndim
      npts = nelt*nv

      call iranku    (cell,irnk,npts,i_n)
      call self_chk  (cell,nv,nelt,32)       ! check for not self-ptg.
      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(name,ie)
      character*40 name
      write(6,*) name, ie
      stop
      end
c-----------------------------------------------------------------------
      subroutine exitt(ie)
      write(6,*) 'exit status', ie
      stop
      end
c-----------------------------------------------------------------------
      subroutine cscan_dxyz (dx,nelt,nelv,ndim,ifbinary,ifbswap)
c
c     Scan for xyz data, read it, and set characteristic length, d2
c

#     include "SIZE" 
     
      character*80 string
c
      real dx(1)
      real x(8),y(8),z(8)
      integer e,buf(50)

      integer h2s(8) ! hypercube to strange ordering
      save    h2s
      data    h2s / 1,2,4,3,5,6,8,7 /

      logical ifbinary,ifbswap

      ifbswap  = .false.

      write(6,*) 'reading mesh data ...'

      if (.not. ifbinary) then
         call cscan(string,'MESH DATA',9)
         read (10,*) nelt,ndim,nelv
      endif
       
      if (nelt.lt.0 .or. ifbinary) then
         ifbinary = .true.
         call open_bin_file(ifbswap,nelgtr,ndimr,nelgvr,wdsizi)
         if(wdsize.eq.4.and.wdsizi.eq.8) then
             write(6,*) "Double Precision .rea not supported ",
     $                  "in Single Precision mode, compile with -r8"
             call exitt(wdsize)
         endif
         nelt = nelgtr
         ndim = ndimr
         nelv = nelgvr
         nwds = (1 + ndim*(2**ndim))*(wdsizi/4) ! group + 2x4 for 2d, 3x8 for 3d
      endif

c      write(6,*) nelt,ndim,nelv,ifbinary, ' nelt,ndim,nelv,ifre2 '

      if (nelt.gt.lelm) then
        write(6,*) 'Abort: number of elements too large', nelt
        write(6,*) 'change MAXNEL and recompile' 
        call exitt(1)
      endif

      b = 1.e22
      l = 1

      ierr = 0
      if (ndim.eq.3) then
         do e=1,nelt
            if(ifbinary) then
              call byte_read(buf,nwds,ierr)
              if(ierr.ne.0) goto 100
              call buf_to_xyz(buf,x,y,z,e,ifbswap,ndim,wdsizi)
            else 
              read (10,80) string
              read (10,*)   (x(k),k=1,4)
              read (10,*)   (y(k),k=1,4)
              read (10,*)   (z(k),k=1,4)
              read (10,*)   (x(k),k=5,8)
              read (10,*)   (y(k),k=5,8)
              read (10,*)   (z(k),k=5,8)
            endif
c           write(6,*) e
c           do ii=1,8
c          write(6,*) x(ii),y(ii),z(ii)
c           enddo
            do k=1,8
               dx(l+0) = b
               dx(l+1) = x(h2s(k))
               dx(l+2) = y(h2s(k))
               dx(l+3) = z(h2s(k))
               l = l + (ndim+1)
            enddo
         enddo
      else
         do e=1,nelt
            if(ifbinary) then
              call byte_read(buf,nwds,ierr)
              if(ierr.ne.0) goto 100
              call buf_to_xyz(buf,x,y,z,e,ifbswap,ndim,wdsizi)
            else
              read (10,80) string
              read (10,*)   (x(k),k=1,4)
              read (10,*)   (y(k),k=1,4)
            endif
            do k=1,4
               dx(l+0) = b
               dx(l+1) = x(h2s(k))
               dx(l+2) = y(h2s(k))
               l = l + (ndim+1)
            enddo
         enddo
      endif
   80 format(a80)
 
      nvrt = 2**ndim
      call set_d2(dx,nvrt,nelt,ndim)
 
      return
 100  write(6,*) "Error reading xyz byte data in scan_dxyz. Abort"
      call exitt(ierr)
      return
      end
c-----------------------------------------------------------------------
      subroutine set_d2(dx,nvrt,nel,ndim)
      real dx(0:ndim,nvrt,nel)

      integer neigh(3,8),e
      save    neigh
      data    neigh / 2,3,5 , 1,4,6 , 1,4,7 , 2,3,8   ! symm. ordering
     $              , 1,6,7 , 2,5,8 , 3,5,8 , 4,6,7 /

      b = 1.e22

      do e = 1,nel
      do i = 1,nvrt

         dx(0,i,e) = b

         do k = 1,ndim

            n   = neigh(k,i)
            d2l = 0.

            do j=1,ndim
               d2l = d2l + (dx(j,n,e)-dx(j,i,e))**2
            enddo

            dx(0,i,e) = min (dx(0,i,e),d2l)
c           write(6,9) i,e,' dx',(dx(j,i,e),j=0,ndim)
c  9        format(2i6,a3,1p4e11.3)

         enddo

      enddo
      enddo
c     call exitt(5)

      return
      end
c-----------------------------------------------------------------------
      subroutine cscan_bcs  (cbc,bc,nelv,nelt,ndim)
c
c     Scan for cbc data and read it
c
      character*3  cbc(6,nelt)
      real*8       bc (5,6,nelt)
      character*80 string

      integer e

      call cscan(string,'BOUNDARY',8) ! for now, fluid only

      npass = 1
      if (nelt.gt.nelv) npass = 2

      do kpass = 1,npass

         read (10,80) string
   80    format(a80)

         ifield = kpass
         if (indx1(string,'NO ',3).ne.0) then
            ifield = ifield+1
            call cscan(string,'BOUNDARY',8) ! then, temp only
         endif

         nel=nelv
         if (kpass.eq.2) nel=nelt
         call rd_bc(cbc,bc,nel,ndim,ifield,10)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rd_bc(cbc,bc,nel,ndim,ifield,io)

c     .Read Boundary Conditions (and connectivity data)

      character*3 cbc(6,nel)
      real*8      bc(5,6,nel)
      integer e,f

c     write(6,*) 'inside rd_bc: ',nel,ndim,ifield,io

      nbcrea = 5
      nface  = 2*ndim
      do e=1,nel
      do f=1,nface
         if (nel.lt.1000) then
            read(io,50,err=510,end=600)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
   50       format(a1,a3,2i3,5g14.6)
         elseif (nel.lt.100 000) then
            read(io,51,err=520,end=600)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
   51       format(a1,a3,i5,i1,5g14.6)
         elseif (nel.lt.1 000 000) then
            read(io,52,err=530,end=600)    
     $      cbc(f,e),id1,(bc(ii,f,e),ii=1,nbcrea)
   52       format(1x,a3,i6,5g14.6)
         else
            read(io,53,err=540,end=600)    
     $      cbc(f,e),id1,(bc(ii,f,e),ii=1,nbcrea)
   53       format(1x,a3,i12,5g18.11)
         endif
c        write(6,*) e,f,' ',cbc(f,e),' BC IN?',nel,ifield
      enddo
      enddo

      return
C
C     Error handling:
C
  500 format(2x,'ERROR: error reading ',i4,2i12,/,
     $       2x,'aborting ',a3,' in routine rdbdry.')
  510 write(6,500) ifield,e,nel,'510'
      call exitt(ifield)
      return

  520 write(6,500) ifield,e,nel,'520'
      call exitt(ifield)
      return

  530 write(6,500) ifield,e,nel,'530'
      call exitt(ifield)
      return

  540 write(6,500) ifield,e,nel,'540'
      call exitt(ifield)
      return

      return

  600 continue
      write(6,601) ifield,e,nel
  601 FORMAT(2X,'ERROR: end of file',i4,2i12,/,
     $       2X,'ABORTING 600 IN ROUTINE RDBDRY.')
      call exitt(ifield)
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
C
      n1=80-l2+1
      indx1=0
      if (n1.lt.1) return
C
      do 300 i=1,n1
         i2=i+l2-1
         if (s1(i:i2).eq.s2(1:l2)) then
            indx1=i
            return
         endif
300   continue
c
      return
      end
c-----------------------------------------------------------------------
      subroutine cscan(sout,key,nk)
      character*80 sout,key
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      do i=1,99999999      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         call chcopy(sout,string,80)
c        write (6,*) string
         if (indx1(string,key,nk).ne.0) return
      enddo
  100 continue
c
   80 format(a80)
      return
      end
c
c-----------------------------------------------------------------------
      subroutine readwrite(sout,key,nk)
      character*80 sout,key
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      do i=1,90000      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
         if (indx1(string,key,nk).ne.0) return
      enddo
  100 continue
c
   80 format(a80)
   81 format(80a1)
      return
      end
c
c-----------------------------------------------------------------------
      subroutine readwrite2(sout,key1,nk1,key2,nk2)
      character*80 sout,key1,key2
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      do i=1,90000      
         call blank(string,80)
         read (10,80,end=100,err=100) string
         len = ltrunc(string,80)
         write(11,81) (string1(k),k=1,len)
c        write(6 ,81) (string1(k),k=1,len)
         if (indx1(string,key1,nk1).ne.0) return
         if (indx1(string,key2,nk2).ne.0) return
      enddo
  100 continue
   80 format(a80)
   81 format(80a1)
      return
      end
c
c-----------------------------------------------------------------------
      integer function log2(k)
      rk=(k)
      rlog=log10(rk)
      rlog2=log10(2.0)
      rlog=rlog/rlog2  + 1.e-6  !  + 0.5  ! don't round up!
      log2=int(rlog)
      return
      end
c-----------------------------------------------------------------------
      function iglmin(a,n)
      integer a(1),tmin
      tmin=999999999
      do 100 i=1,n
         tmin=min(tmin,a(i))
  100 continue
      iglmin=tmin
      return
      end
c-----------------------------------------------------------------------
      function ivlmax(a,n)
      integer a(1),tmax
      tmax=-999999999
      do 100 i=1,n
         tmax=max(tmax,a(i))
  100 continue
      ivlmax=tmax
      return
      end
c-----------------------------------------------------------------------
      function iglmax(a,n)
      integer a(1),tmax
      tmax=-999999999
      do 100 i=1,n
         tmax=max(tmax,a(i))
  100 continue
      iglmax=tmax
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
      subroutine icopy(x,y,n)
      integer x(1),y(1)
      do i=1,n
         x(i) = y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cmult2(x,y,c,n)
      real x(1),y(1)
      do i=1,n
         x(i) = c*y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy8(x,y,n)
      real*8 x(1),y(1)
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
      subroutine cfill(x,c,n)
      real x(1)
      do i=1,n
         x(i) = c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine rint(x,n)
      real x(1)
      do i=1,n
         x(i) = i
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine jjnt(x,n)
      integer x(1)
      do i=1,n
         x(i) = i
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ifill(x,c,n)
      integer x(1),c
      do i=1,n
         x(i) = c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine icadd(x,c,n)
      integer x(1),c
      do i=1,n
         x(i) = x(i)+c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cadd(x,c,n)
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
      subroutine getreafile(prompt,ifbinary,io,ierr)
c
      character*1 prompt(1)
      logical ifbinary
c
      common /sess/ session
      character*80 session

      character*80 file
      character*1  file1(80)
      equivalence (file1,file)

      ierr = 0
      ifbinary = .false.

c     Get file name
      len = indx1(prompt,'$',1) - 1
      write(6,81) (prompt(k),k=1,len)
   81 format(80a1)
      call blank(session,80)
      read(5,80) session
   80 format(a80)

      if (session.eq.'-1') then
         io = -1
         ierr = 1
         return
      else
         call chcopy(file,session,80)
         len = ltrunc(file,80)
         call chcopy(file1(len+1),'.rea',4)
         open(unit=io, file=file, status='old', iostat=ierr)

         if (ierr.gt.0) then
            call chcopy(file,session,80)
            len = ltrunc(file,80)
            call chcopy(file1(len+1),'.re2',4)
            inquire(file=file, exist=ifbinary)
            if(ifbinary) ierr = 0
         endif
      endif

      write(6,*) 'reading ', file

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
      subroutine iswap_ip(x,p,n)
      integer x(1),xstart
      integer p(1)
c
c     In-place permutation: x' = x(p)
c
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            xstart     = x(k)
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! iswap_ip problem.',j,k,n,next
                  call exitt(0)
               elseif (next.eq.loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine iswapt_ip(x,p,n)
      integer x(1),t1,t2
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            t1         = x(loop_start)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt(1)
               elseif (next.eq.loop_start) then
                  x(next) = t1
                  p(next) = -p(next)
                  goto 10
               else
                  t2      =  x(next)
                  x(next) =  t1
                  t1      =  t2
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine sort(a,p,n)
c
c     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
c
c     returns sorted list a(), and permutation vector p()
c
      real    a(1),q
      integer p(1)
C
      if (n.le.1) return
      DO 10 J=1,N
         p(j)=j
   10 continue
C
      L=n/2+1
      ir=n
  100 CONTINUE
         IF (l.gt.1) THEN
            l=l-1
            indx=p(l)
            q=a(indx)
         ELSE
            indx=p(ir)
            q=a(indx)
            p(ir)=p(1)
            ir=ir-1
            if (ir.eq.1) then
               p(1)=indx
               return
            endif
         endif
         i=l
         j=l+l
  200    CONTINUE
         IF (J.le.IR) THEN
            IF (J.lt.IR) THEN
               IF ( A(p(j)).lt.A(p(j+1)) ) j=j+1
            endif
            IF (q.lt.A(p(j))) THEN
               p(I)=p(J)
               I=J
               J=J+J
            ELSE
               J=IR+1
            endif
         goto 200
         endif
         p(I)=indx
      goto 100
      end
c-----------------------------------------------------------------------
      subroutine csr_sort_colj(ia,ja,n,ind)
      integer ia(0:1),ja(1)
c
      do i=1,n
         j0 = ia(i-1)
         nj = ia(i)-j0
         call isort(ja(j0),ind,nj)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_csrmati(ia,ja,n,name9)
      integer ia(0:1),ja(1)
c
      character*9 name9
      character*4 s(33)
c
      write(6,1) name9,n
    1 format(/,'CSR Mat:',a9,3x,'n =',i9,/)
c
      mj = 11
      do i=1,n
         j0 = ia(i-1)
         j1 = ia(i)-1
         nj = j1-j0 + 1
         j1 = min(j1,j0+mj)
         if (nj.gt.0) write(6,26) i,nj,(ja(k),k=j0,j1)
      enddo
   26 format(2i4,': ',26i4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine iranku(a,rank,n,p)
c
c     Return r = rank of a() in place of a(), where r =< n
c     so that a() is the number of unique pts sorted such that
c     it is in the original spot
c     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
c
c
      integer a(1),p(1)
      integer rank

      write(6,*) 'compressing verticies'

      call isort(a,p,n)

      rank = 1
      last = a(1)
      a(1) = 1

      do k=2,n
         next = a(k)
         if (next.ne.last) rank=rank+1
         a(k) = rank
         last = next
      enddo

      call iswapt_ip(a,p,n)  ! restore a() to orginal location

      return
      end
c-----------------------------------------------------------------------
      subroutine dmpfile(cell, nv, nelt, nelv, nuids)

#     include "SIZE"

      common /sess/ session
      character*80 session
      character*80 fname
      character*1  fnam1(80)
      equivalence (fnam1,fname)

      integer cell(nv,nelt)

      integer d2,e,p0
      common /arrayi2/ iwrk((1+8)*lelm)

      character*132 hdr
      character*5   version
      real*4 test
      data   test  / 6.54321 /

      version = '#v001'
c      ifco2 = .false. ! force ASCII for debugging 
 
      len = ltrunc(session,80)
      call chcopy(fname,session,80)
      if (ifco2) then
         call chcopy(fnam1(len+1),'.co2',4)
      else
         call chcopy(fnam1(len+1),'.con',4)
      endif

      if (nv.ne.2)  then 
         write(6,'(A,A)') 'writing ', fname
         if (ifco2) then
            call byte_open(fname,ierr)
            call blank(hdr,132)
            write(hdr,1) version,nelt,nelv,nv
    1       format(a5,3i12)
            call byte_write(hdr,132/4,ierr)
            call byte_write(test,1,ierr) ! write the endian discriminator
         else
            open (unit=29,file=fname)
            write(29,1) version,nelt,nelv,nv
        endif
      endif

      if (ifco2) then
         do e=1,nelt
            iwrk(1) = e
            call icopy(iwrk(2),cell(1,e),nv)
            call byte_write(iwrk,nv+1,ierr)
         enddo
      else
         do e=1,nelt
            write(29,2) e,(cell(k,e),k=1,nv)
    2       format(9i12)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ituple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      integer a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_ialtb
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
c           aa  = a  (l)
            call icopy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call icopy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call icopy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call icopy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
               if (iftuple_ialtb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
            if (iftuple_ialtb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call icopy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call icopy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      subroutine tuple_sort(a,lda,n,key,nkey,ind,aa)
C
C     Use Heap Sort (p 231 Num. Rec., 1st Ed.)
C
      real a(lda,n),aa(lda)
      integer ind(1),key(nkey)
      logical iftuple_altb
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
c           aa  = a  (l)
            call copy(aa,a(1,l),lda)
            ii  = ind(l)
         else
c           aa =   a(ir)
            call copy(aa,a(1,ir),lda)
            ii = ind(ir)
c           a(ir) =   a( 1)
            call copy(a(1,ir),a(1,1),lda)
            ind(ir) = ind( 1)
            ir=ir-1
            if (ir.eq.1) then
c              a(1) = aa
               call copy(a(1,1),aa,lda)
               ind(1) = ii
               return
            endif
         endif
         i=l
         j=l+l
  200    continue
         if (j.le.ir) then
            if (j.lt.ir) then
c              if ( a(j).lt.a(j+1) ) j=j+1
               if (iftuple_altb(a(1,j),a(1,j+1),key,nkey)) j=j+1
            endif
c           if (aa.lt.a(j)) then
            if (iftuple_altb(aa,a(1,j),key,nkey)) then
c              a(i) = a(j)
               call copy(a(1,i),a(1,j),lda)
               ind(i) = ind(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
         GOTO 200
         endif
c        a(i) = aa
         call copy(a(1,i),aa,lda)
         ind(i) = ii
      GOTO 100
      end
c-----------------------------------------------------------------------
      logical function iftuple_ialtb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_ialtb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_ialtb = .false.
            return
         endif
      enddo
      iftuple_ialtb = .false.
      return
      end
c-----------------------------------------------------------------------
      logical function iftuple_altb(a,b,key,nkey)
      real a(1),b(1)
      integer key(nkey)
c
      do i=1,nkey
         k=key(i)
         if (a(k).lt.b(k)) then
            iftuple_altb = .true.
            return
         elseif (a(k).gt.b(k)) then
            iftuple_altb = .false.
            return
         endif
      enddo
      iftuple_altb = .false.
      return
      end
c-----------------------------------------------------------------------
      logical function iftuple_ianeb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      iftuple_ianeb = .false.
      do i=1,nkey
         k=key(i)
         if (a(k).ne.b(k)) then
            iftuple_ianeb = .true.
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      logical function iftuple_iaeqb(a,b,key,nkey)
      integer a(1),b(1)
      integer key(nkey)
c
      iftuple_iaeqb = .true.
      do i=1,nkey
         k=key(i)
         if (a(k).ne.b(k)) then
            iftuple_iaeqb = .false.
            return
         endif
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine tuple_swapt_ip(x,m,n,p,t1,t2)
      real x(m,n),t1(m),t2(m)
      integer p(1)
c
c     In-place permutation: x'(p) = x
c

      do k=1,n
         if (p(k).gt.0) then   ! not swapped
            loop_start = k
            next       = p(loop_start)
            call copy(t1,x(1,loop_start),m)
            do j=1,n
               if (next.lt.0) then
                  write(6,*) 'Hey! iswapt_ip problem.',j,k,n,next
                  call exitt(1)
               elseif (next.eq.loop_start) then
                  call copy(x(1,next),t1,m)
                  p(next) = -p(next)
                  goto 10
               else
                  call copy(t2,x(1,next),m)
                  call copy(x(1,next),t1,m)
                  call copy(t1,t2       ,m)
                  nextp   =  p(next)
                  p(next) = -p(next)
                  next    =  nextp
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine periodic_vtx
     $               (cell,nv,nel,irnk,dx,ndim,cbc,bc,nfc,iper,jmin)
c
c     NOTE:  .genmap() uses the symmetric vertex/face ordering.
c            .cbc() is ordered in symm. fashion
c            .Contents of cbc(), however, refer to nekton preprocessor
c             ordering, _not_ to the symmetric ordering.
c
c
c     Reassign cell() pointers to account for periodic connections.
c
      integer cell(nv,nel),iper(ndim,1),jmin(1)
      real dx(0:ndim,nv,nel)

      character*3      cbc(  6,nel)
      real*8           bc (5,6,nel)

      character*3      cb,cj
      integer e,f,v
      integer ipair(2,4)

      integer h2s(8) ! hypercube to strange ordering
      save    h2s
      data    h2s / 1,2,4,3,5,6,8,7 /

      integer eface(6)  ! return Nekton preprocessor face ordering
      save    eface
      data    eface / 4 , 2 , 1 , 3 , 5 , 6 /

      integer efaci(6)  ! return symmetric face ordering
      save    efaci
      data    efaci / 3 , 2 , 4 , 1 , 5 , 6 /

      logical foundp

      write(6,*) 'start periodic vtx:',nel,irnk

      call izero(iper,ndim*irnk)   ! Zero out periodic indicator flag
      call jjnt (jmin,irnk)        ! Initial permutation = identity
      nvf = nv/2                   ! # vertices/face = 1/2 # vertices/cell

      nmn = irnk
      nmx = 0
      foundp = .false.
      do e=1,nel
      do f=1,nfc
         cb = cbc(f,e)
c        write(6,*) cb,e,f,' cb'
         if (cb.eq.'P  ') then
           foundp = .true.
           je = abs(bc(1,f,e))
           jf = bc(2,f,e)
           jf = efaci(jf)

           cj = cbc(jf,je)
           ke = abs(bc(1,jf,je))
           kf = bc(2,jf,je)
           kf = efaci(kf)

c          write(26,26)   e,f,cb,je,jf,cj,ke,kf,nel
c  26      format(i9,i2,1x,a3,i9,i2,1x,a3,i9,i2)

           if (bc(1,f,e).gt.0 .and. bc(1,jf,je).gt.0) then
              if (ke.ne.e .or. kf.ne.f .or. cj.ne.'P  ') then
               write(6,*)
               write(6,*) 'abort: PERIODIC MISMATCH 1:'
               write(6,6)   e,f,cb,' ie '
               write(6,6) je,jf,cj,' je '
               write(6,6) ke,kf,cj,' ke '
               write(6,*)
    6          format(i12,i3,1x,a3,1x,a4)
               call exitt(9)
              endif

              call find_connctd_pairs
     $                 (jmin,nvf,e,f,je,jf,cell,nv,dx,ndim,nel)
            

c             bc(1, f, e) = -bc(1, f, e) ! indicate that 
c             bc(1,jf,je) = -bc(1,jf,je) ! pairing is done !!!

           elseif (bc(1,f,e)*bc(1,jf,je).le.0) then
              write(6,*)
              write(6,*) 'abort: PERIODIC MISMATCH 2:'
              write(6,6)   e,f,cb,' ie '
              write(6,6) je,jf,cj,' je '
              write(6,*) bc(1,f,e),bc(1,jf,je),' bc'
              write(6,*)
              call exitt(8)
           endif
         endif
      enddo
      enddo

      if (.not.foundp) goto 99

c     Okay -- we now have the updated pointers, time to update
c     the cell pointers
      npts = nv*nel
      do i=1,npts
         cell(i,1) = jmin(cell(i,1))  ! permuted identity
      enddo
      call iranku      (cell,irnk,npts,iper) ! compress cell list
      call self_chk    (cell,nv,nel,6)       ! check for not self-ptg.

c     Reset bc array
      do e=1,nel
      do f=1,nfc
         if (cbc(f,e).eq.'P  ') bc(1,f,e) = abs(bc(1,f,e))
      enddo
      enddo

 99   write(6,*) 'done periodic vtx', foundp

      return
      end
c-----------------------------------------------------------------------
      subroutine find_connctd_pairs
     $   (jmin,nvf,e,f,je,jf,cell,nv,dx,ndim,nel)


      integer jmin(1),cell(nv,1)
      real dx(0:ndim,nv,1)

      real x0(0:3,4),x1(0:3,4)
      real z0(0:3,4),z1(0:3,4)

      integer e,f,shift,smin

      integer vface(4,6)  ! circulant vertices on symm. faces, 3D
      save    vface       ! order ctr-clkws when looking at face
      data    vface / 1,5,7,3 , 2,4,8,6 , 1,2,6,5 , 3,7,8,4
     $              , 1,3,4,2 , 5,6,8,7 /
      integer wface(4,6)  ! circulant vertices on symm. faces, 2D
      save    wface       ! 
      data    wface / 3,1,0,0 , 2,4,0,0 , 1,2,0,0 , 4,3,0,0
     $              , 0,0,0,0 , 0,0,0,0 /

      integer icalld
      save    icalld
      data    icalld /0/

      icalld = icalld+1

      if (ndim.eq.2) call icopy(vface,wface,24)

      nvf = nv/2     ! # vertices/face = 1/2 # vertices/cell
      
c     write(6,4) e,f,je,jf,nv,nvf,ndim,' FACE'
c   4 format(i9,i3,i9,i3,3i3,a5)

      do i=1,nvf     ! Grab geometry for P-P face pair

         j0 = vface (i ,f)
         call copy  (x0(0,i),dx(0,j0, e),ndim+1)

         j1 = vface (i ,jf)
         call copy  (x1(0,i),dx(0,j1,je),ndim+1)

      enddo
      call copy(z0,x0,16)   ! For failure diagnosis
      call copy(z1,x1,16)   ! For failure diagnosis


      x0m = 0.
      do k=1,ndim           ! Subtract off mean of faces
         x0a = 0.
         x1a = 0.
         do i=1,nvf
            x0a = x0a + x0(k,i)
            x1a = x1a + x1(k,i)
         enddo
         do i=1,nvf
            x0m = max(x0m,abs(x0(k,i)))
            x0m = max(x0m,abs(x1(k,i)))
            x0(k,i) = x0(k,i) - x0a/nvf
            x1(k,i) = x1(k,i) - x1a/nvf
         enddo
      enddo

c     call outmat(x0,4,4,'  x0  ', e)
c     call outmat(x1,4,4,'  x1  ',je)

      d2min = 1.e22
      do shift=0,nvf-1
         d2 = 0.
         do i=1,nvf
            j=i+shift
            if (j.gt.nvf) j=j-nvf
            j=nvf+1-j              ! go backward for j !
c           write(6,5)
            do k=1,ndim
               d2 = d2 + (x0(k,i)-x1(k,j))**2
c              write(6,5) shift,i,j,k,x0(k,i),x1(k,j),d2
c   5          format(4i4,1p3e12.4,'  d2')
            enddo
         enddo
         if (d2.lt.d2min) then
            smin  = shift
            d2min = d2
         endif
      enddo
      shift = smin
      if (d2min.gt.0) d2min = sqrt(d2min)
      eps = 1.e-4
      eps = 1.e-3
      tol = eps*x0m
      if (d2min.gt.tol) then
c         call outmat(x0,4,4,'  x0  ', e)
c         call outmat(x1,4,4,'  x1  ',je)
c         call outmat(z0,4,4,'  z0  ', e)
c         call outmat(z1,4,4,'  z1  ',je)
         write(6,6) e , f,shift,eps,x0m
         write(6,6) je,jf,    i,tol,d2min
    6    format(i12,i2,i3,1p2e16.8,' abort: FACE MATCH FAIL')
         call exitt(0)
      endif

c      if (nel.le.100000.or.mod(icalld,1000).eq.0)
c     $   write(6,7) e,f,i,shift,d2min,icalld
c    7    format(i12,i2,2i3,1p1e16.8,i9,' shift')

      do i=1,nvf

         j=i+shift
         if (j.gt.nvf) j=j-nvf
         j=nvf+1-j              ! go backward for j to make faces match

         iv = vface(i, f)
         jv = vface(j,jf)

         ic = cell(iv, e)
         jc = cell(jv,je)

         ijmin = jmin(ic)
         jjmin = jmin(jc)

         jmin(ic) = min(ijmin,jjmin)
         jmin(jc) = min(ijmin,jjmin)

      enddo

      return
      end
c-----------------------------------------------------------------------
      function ipush(stack,val)
      integer stack(0:1),val

      n    = stack(0)     
      if (n.ge.0) then
         n        = n+1
         stack(n) = val
         stack(0) = n
         ipush    = n
      else
         ipush    = -1
         write(6,*) 'ipush: negative stack count'
      endif

      return
      end
c-----------------------------------------------------------------------
      function ipop(nout,stack)
      integer stack(0:1)

      n    = stack(0)     
      if (n.gt.0) then
         ipop     = stack(n)
         nout     = n-1
         stack(0) = nout
      else
         nout     = 0
         ipop     = 0
c        write(6,*) 'ipop: stack empty'
      endif

      return
      end
c-----------------------------------------------------------------------
      logical function is_connected(jactive,n0,ia,ja,n,jstack)

      integer ia(0:n),ja(1),jactive(n),jstack(0:n)

      character*1 adum

c     if ia     < 0, then already on list
c     if ja(ia) < 0, then row is already processed


      call izero(jactive,n)
      call izero(jstack ,n+1)

      m = ia(n)

      icurr      = 1
      jactive(1) = ia(icurr-1)

      do i=1,m             ! move through all entries in array

         j0 = ia     (icurr-1)
         j1 = ia     (icurr)
         jj = jactive(icurr)

         do j=jj,j1-1
            inext = ja(j)
            if (1.le.inext.and.inext.le.n) then  ! range check
               if (jactive(inext).eq.0) then
                  nstack = ipush(jstack,icurr)   ! push the return point
                  jactive(icurr) = j+1
                  jactive(inext) = ia(inext-1)
                  goto 10
               endif
            endif
         enddo

         jactive(icurr) = j1                     ! we've exhausted this row
         inext = ipop(nstack,jstack)

   10    icurr = inext
         if (icurr.eq.0) goto 100

      enddo
  100 continue

      n0           = 0
      is_connected = .true.

      do i=1,n
         if (jactive(i).eq.0) then
            n0           = n0+1
            is_connected = .false.
         endif
      enddo


c     write(6,*) is_connected,n,nstack,' is connected?'
c     if (.not.is_connected) then
c        n8 = min(n,18)
c        write(6,8) n,'jact:',(jactive(k),k=1,n8)
c        write(6,8) n,'ia:  ',(ia(k),k=0,n8-1)
c 8      format(i4,1x,a5,1x,18i5)
c        n8 = ia(n)-ia(0)-1
c        n8 = min(n8,18)
c        write(6,8) n8,'ja:  ',(ja(k),k=1,n8)
c     endif
c     call exitt(8)
c     read (5,*) adum

      return
      end
c-----------------------------------------------------------------------
      subroutine ident(a,n)
      real a(n,n)
      call rzero(a,n*n)
      do i=1,n
         a(i,i) = 1.0
      enddo
      return
      end
c=======================================================================
      subroutine rzero(a,n)
      real  a(1)
      do i = 1, n
         a(i) = 0.0
      enddo
      return
      end
c=======================================================================
      function glsc2(x,y,n)
      real x(1), y(1)
      s = x(1)*y(1)
      do i=2,n
         s = s+x(i)*y(i)
      enddo
      glsc2 = s
      return
      end

c=======================================================================
      subroutine add2s2(a,b,c,n)
      real  a(1),b(1)
      do i = 1, n
         a(i) = a(i) + c*b(i)
      enddo
      return
      end
c=======================================================================
      subroutine add2s1(a,b,c,n)
      real  a(1),b(1)
      do i = 1, n
         a(i) = c*a(i) + b(i)
      enddo
      return
      end

c=======================================================================
      subroutine reverse_p (p,n)          ! lightly load node 0
      integer p(n),pmax
c
      pmax = p(1)
      do i=1,n
         pmax = max(pmax,p(i))
      enddo
c      write(6,*) 'pmax:',pmax,n

      do i=1,n
         p(i) = pmax+1 - p(i)  ! range of p is [1:pmax]
      enddo
      return
      end

c-----------------------------------------------------------------------
      subroutine unique_vertex2(cell,dx,ndim,nel,q,ind,ninseg,ifseg,wk)
 
      integer cell(1),ind(1),ninseg(1)
      logical ifseg(1)
      real dx(0:ndim,1)

      integer e
      real dxt(4),t1(4),t2(4),wk(0:ndim,1)


      nvtx = 2**ndim
      n = nvtx*nel

      write(6,*) 'start locglob_lexico:',nvtx,nel,n,q

      qq = q*q  ! Square of relative tolerance

      do i=1,n
         cell(i)   = i
         ifseg (i) = .false.
      enddo

c     Sort by directions

      lda         = 1+ndim
      nseg        = 1
      ifseg(1)    = .true.
      ninseg(1)   = n

      do ipass=1,ndim   ! Multiple passes eliminates false positives
      do j=1,ndim       ! Sort within each segment
         write(6,*) 'locglob:',j,nseg,n
         i =1
         j1=j+1
         do iseg=1,nseg
            call tuple_sort(dx(0,i),lda,ninseg(iseg),j1,1,ind,dxt) !key=j1
            call iswap_ip  (cell(i),ind,ninseg(iseg)) ! Swap position
            i  =   i + ninseg(iseg)
         enddo
 
         do i=2,n
           if ((dx(j,i)-dx(j,i-1))**2.gt.qq*min(dx(0,i),dx(0,i-1)))
     $        ifseg(i)=.true.
         enddo

         nseg = 0              !  Count up number of different segments
         do i=1,n
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
      enddo
c
c     Assign global node numbers (now sorted lexigraphically!)
c
      ig  = 0
      ic  = 0
      icm = 0
      do i=1,n
        ic = ic+1     ! count number of instances at present ig

        if (ifseg(i)) then
           ig=ig+1
           icm = max(ic,icm)
           ic  = 0
        endif
        ind(cell(i)) = ig
      enddo
      nglb = ig

c     Unshuffle geometry:
       call tuple_swapt_ip(dx,lda,n,cell,t1,t2)

c     Reassign cell to hold global index numbering

      call icopy     (cell,ind,n)

c      call self_chk  (cell,nvtx,nel,0)       ! check for not self-ptg.


c     Reassign geometry to match global index numbering
c     Retain the geometry that is associated with the smallest bounding radius

c     call copy          (wk,dx,lda*n)
c     call izero         (ind,nglb)      ! Flag to see if value already assigned

c     do i=1,n
c        ig = cell(i)
c        if (ind(ig).eq.0) then
c           call copy(dx(0,ig),wk(0,i),lda) ! lda = ndim+1 values ->global vtx
c           ind(ig) = 1
c        elseif (wk(0,i).lt.dx(0,ig)) then
c           call copy(dx(0,ig),wk(0,i),lda) ! lda = ndim+1 values ->global vtx
c        endif
c     enddo
      write(6,6) nseg,nglb,n,icm
    6 format(' done locglob_lexico:',4i12)


      return
      end
c-----------------------------------------------------------------------
      subroutine outcell(cell,dx,n,ndim)

      integer cell(n)
      real dx(0:ndim,n)
    
      do i = 1,n
      write(6,6) cell(i),i,dx(1,i),dx(2,i)
      enddo
   6  format(i8,i3,3g15.6)
   
      return
      end
c-----------------------------------------------------------------------
      subroutine face_chk(face_conn,cell,
     $                        nv,ncell,nic,njc,v2cind,v2c,w3)
c     checks that every face has proper connectivity
c        i.e. if a 3D face has 3 vertices connected a another, 
c             it must have a fourth
c     vc2    : vertex to cell pointer
c     vc2ind : indicies for the vertex to cell pointer
c     cell   : global vertex numbering

      integer vface(4,6)  ! symm. vertices ordered on symm. faces
      save    vface
      data    vface / 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8
     $              , 1,2,3,4 , 5,6,7,8 /

      logical face_conn,is_ok,n_is_3
      integer cell(nv,ncell),v2c(1),v2cind(1),w3(200),share(200)
      integer f,v,vfail,v2

      ndim = 3
      itype = 1          ! vertex to cell pointers !
      call cell2v(i0,i1,v2cind,nic,v2c,njc,cell,nv,ncell,itype,w3) 

      nface = ndim*2     ! number of faces per element
      nfv   = 2*(ndim-1) ! number of vertices in each face

      face_conn = .true.      

      do i = 1,ncell
      do f = 1,nface
         n_shared = 0     
         do v = 1,nfv
            iv1 = cell(vface(v,f),i)   
            j0  = v2cind(iv1)
            j1  = v2cind(iv1+1)-1
            do j=j0,j1
               je = v2c(j)
               if (je.ne.i) then
                  n_shared = n_shared+1
                  share(n_shared) = je
               endif
            enddo
         enddo

         call isort(share,w3,n_shared)

         ilast = 0
         icount = 0
         n_is_3 = .false.

         do k = 1,n_shared
            if (share(k).eq.ilast) then
               icount = icount+1
               if(icount.eq.3) then
                  n_is_3 = .true.
               elseif(icount.eq.nfv) then
                  n_is_3 = .false.
               endif
            else
              icount = 1
              ilast  = share(k)
              if(n_is_3) then  !FAIL
                 n_fail = share(k-1)
                 write(6,6) f,i,n_fail,k,n_shared
    6            format('MISSING FACE CONNECTION  a',i3,2i10,i4,i3)
                 write(6,*) (share(kk),kk=1,n_shared)
                 goto 10
              endif
            endif
         enddo
         if (n_is_3) then !FAIL
            n_fail = share(k-1)
            write(6,*) 'MISSING FACE CONNECTION  b',i,f,n_fail,k
            goto 10
         endif
      enddo
      enddo

      return

   10 continue
      do vfail = 1,nfv
         iv1 = cell(vface(vfail,f),i)
         j0    = v2cind(iv1)
         j1    = v2cind(iv1+1)-1
         is_ok = .false.
         do j  = j0,j1          ! Checks that pt_iv1 connects to k 
            if (v2c(j).eq.n_fail) is_ok = .true.
         enddo
         if(.not.is_ok) goto 20 ! Then this pt is our failed one
      enddo
 
  20  continue
      vfail = vface(vfail,f)             ! Failed vertex, on element i
      call find_v2(v2,i,n_fail,cell,ncell,nv,v2cind,v2c)
      call fix_geom(cell,i,vfail,n_fail,v2,v2cind,v2c,ncell,nv)
      face_conn = .false.

      return
      end
c-----------------------------------------------------------------------
      subroutine find_v2(v2,e1,e2,cell,nel,nv,ind,jc)
c
c     Find the vertex(v2) on e2 that (v1,e1) should be connected to      
c
      integer cell(nv,nel),ind(1),jc(1)
      integer e1,e2,v2
      integer f,v,vf,ncount
      logical is_ok
 
      integer vface(4,6)  ! symm. vertices ordered on symm. faces
      save    vface
      data    vface / 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8
     $              , 1,2,3,4 , 5,6,7,8 /

      do f = 1,6
         ncount = 0
         do v = 1,4
            iv2 = cell(vface(v,f),e2)   
            j0  = ind(iv2)
            j1  = ind(iv2+1)-1
            is_ok = .false.
            do j=j0,j1
               if(jc(j).eq.e1) then
                 ncount = ncount +1
                 is_ok = .true.
               endif
            enddo
            if(.not.is_ok)vf = v
         enddo
         if(ncount.eq.3) goto 10
      enddo
  10  continue
      v2 = vface(vf,f)   

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_geom(cell,e1,v1,e2,v2,ind,jc,nel,nv)
c
c     Changes the wrong point (v1,e1) to the correct point (v2,e2)
c     Then checks that no other (v1,e1) exist in the geom.
c
      integer cell(nv,nel),ind(1),jc(1)
      integer e1,v1,e2,v2,wrong_cell
      

      iv1 = cell(v1,e1)
      iv2 = cell(v2,e2)
      ivt = min(iv1,iv2)         ! what will now be the "correct" cell
      ivf = max(iv1,iv2)         ! what will now be the "incorrect" cell

      j0 = ind(ivf)
      j1 = ind(ivf)
      do j = j0,j1
         ke = jc(j)
         do k = 1,nv
            if (cell(k,ke).eq.ivf) cell(k,ke)=ivt
         enddo
      enddo

      
      return
      end
c-----------------------------------------------------------------------
      function mod1(i,n)
C
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C
      mod1=0
      if (i.eq.0) return
      if (n.eq.0) then
         write(6,*)
     $  'WARNING:  Attempt to take MOD(I,0) in FUNCTION MOD1.'
         return
      endif
      ii = i+n-1
      mod1 = mod(ii,n)+1
      return
      end
c-----------------------------------------------------------------------
      subroutine self_chk(cell,nv,nel,flag)     ! check for not self-ptg.
      integer cell(nv,nel),flag
      integer e

      write(6,*) 'performing self_chk'

      do e=1,nel
      do i=1,nv
         do j=i+1,nv
            if (cell(i,e).eq.cell(j,e)) then

               write(6,*)
               write(6,*) 'ABORT: SELF-CHK ',i,j,e,flag
               write(6,*) 'Try to tighten the mesh tolerance!' 

               call exitt(flag)

            endif
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rd_bc_bin(cbc,bc,nelv,nelt,ifbswap)

c     .Read Boundary Conditions (and connectivity data)

#     include "SIZE"

      character*3 cbc(6,lelm)
      real*8      bc (5,6,lelm)
      logical     ifbswap
      
      integer e,f,buf(30)

      npass = 1
      if (nelt.gt.nelv) npass = 2   ! default to thermal topology (for now)

      ierr = 0
      do kpass = 1,npass

         do e=1,nelt   ! fill up cbc w/ default
         do k=1,6
            cbc(k,e) = 'E  '
         enddo
         enddo

         nwds =(2 + 1 + 5)*(wdsizi/4)   ! eg + iside + cbc + bc(5,:,:)

         if(wdsizi.eq.8) then
            call byte_read(rbc_max,2,ierr)
            if (ifbswap) call byte_reverse8(rbc_max,2,ierr) 
            nbc_max = rbc_max
            do k=1,nbc_max
c              write(6,*) k,' dobc1 ',nbc_max
               call byte_read(buf,nwds,ierr)
               n8wds=nwds/2
               if (ifbswap) call byte_reverse8(buf,nwds-2,ierr) ! last is char
               if(ierr.ne.0) call exitti
     &              ('Error reading byte bcs ',ierr)
               call buf_to_bc(cbc,bc,buf,nelt)
            enddo
         else
            call byte_read(nbc_max,1,ierr)
            if (ifbswap) call byte_reverse(nbc_max,1,ierr) 
            do k=1,nbc_max
c              write(6,*) k,' dobc1 ',nbc_max
               call byte_read(buf,nwds,ierr)
               if (ifbswap) call byte_reverse(buf,nwds-1,ierr) ! last is char
               if(ierr.ne.0) call exitti
     &              ('Error reading byte bcs ',ierr)
               call buf_to_bc(cbc,bc,buf,nelt)
            enddo
         endif

      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbl,bl,buf,nelt)  ! version 1 of binary reader

#     include "SIZE"

      character*3 cbl(6,nelt)
      real*8      bl(5,6,nelt)
      integer     buf(30)

      integer e,eg,f

      if(wdsizi.eq.8) then
         call copyi4(e,buf(1),1) !1-2
         call copyi4(f,buf(3),1) !3-4
         call copy  (bl(1,f,e),buf(5),5) !5--14
         call chcopy(cbl( f,e),buf(15),3)!15-16
c     the following is not needed since buf(5)-buf(6) are already real
c         if(nelt.ge.1000000.and.cbl(f,e).eq.'P  ') 
c     $   call copyi4(bl(1,f,e),buf(5),1) !Integer assign connecting P element
      else
         e  = buf(1)
         f  = buf(2)
         call copy48r ( bl(1,f,e),buf(3),5)
         call chcopy  (cbl(  f,e),buf(8),3)
         if(nelt.ge.1000000.and.cbl(f,e).eq.'P  ') 
     $     bl(1,f,e) = buf(3) ! Integer assign of connecting periodic element
      endif


c      write(6,1) e,f,cbl(f,e),(bl(k,f,e),k=1,5),' CBC'
c  1   format(i8,i4,2x,a3,5f8.3,1x,a4)

      return
      end

c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,xc,yc,zc,e,ifbswap,ndim,wdsizi)  
c     version 1 of binary

      logical ifbswap
      integer e,eg,buf(0:49),wdsizi

      real xc(8),yc(8),zc(8)   !these are *8

      nwds = (1 + ndim*(2**ndim)) ! group + 2x4 for 2d, 3x8 for 3d

      ierr = 0

      if (ifbswap.and.wdsizi.eq.8)     then
          nwds=nwds*2
          call byte_reverse8(buf,nwds,ierr)
      elseif (ifbswap.and.wdsizi.eq.4) then 
          call byte_reverse (buf,nwds,ierr)
      endif

      if (ierr.ne.0) call exitti
     $   ('Error byte_reverse in buf_to_xy ',ierr)

      if(wdsizi.eq.8) then
         call copyi4(igroup,buf(0),1) !0-1
         if (ndim.eq.3) then 
            call copy  (xc,buf( 2),8) !2 --17
            call copy  (yc,buf(18),8) !18--33
            call copy  (zc,buf(34),8) !34--49
         else
            call copy  (xc,buf( 2),4) !2 --9
            call copy  (yc,buf(10),4) !10--17
         endif
      else 
         igroup = buf(0)
         if (ndim.eq.3) then
            call copy4r(xc,buf( 1),8)
            call copy4r(yc,buf( 9),8)
            call copy4r(zc,buf(17),8)
         else
            call copy4r(xc,buf( 1),4)
            call copy4r(yc,buf( 5),4)
         endif
      endif

      return
      end

      subroutine copyI4(a,b,n)
      integer A(1)
      REAL   B(1)
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      end
c-----------------------------------------------------------------------
      subroutine copy4(a,b,n)
      real*4 a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy4r(a,b,n)
      real   a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine copy48r(a,b,n)
      real*8 a(1)
      real*4 b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cell2v(i0,i1,ic,nic,jc,njc,cell,nv,ncell,type,wk)
c
c     Generate a list of cells that point to each numbered
c     vertex in [i0,i1].  ic/jc is the csr-formated list.
c
c     Cells (aka elements) have nv local verties.
c
      integer ic(1),jc(1),cell(nv,ncell),wk(1),type
c
      i0 = iglmin(cell,nv*ncell)
      i1 = iglmax(cell,nv*ncell)

      if (i1-i0.gt.nic) then
         write(6,1) i0,i1,nic,nv,ncell,njc
    1    format(' ERROR: nic too small in cell2v:',6i12)
         i0 = 0 ! error return code
         i1 = 0 ! error return code
         return
      endif

      call cell2v1(ic,i0,i1,jc,njc,cell,nv,ncell,type,wk)

      return
      end
c-----------------------------------------------------------------------
      subroutine cell2v1(ic,i0,i1,jc,njc,cell,nv,ncell,type,wk)
c
c     Generate a list of cells that point to each numbered
c     vertex in [i0,i1].  ic/jc is the csr-formated list.
c
c     Cells (aka elements) have nv local verties.
c
      integer ic(i0:i1),jc(1),cell(nv,ncell),wk(i0:i1),type

c     Step 1:  count number of cells per vertex
c  the number of cells each vertex 'touches'

      nic = i1-i0+1
      call izero(wk,nic)
      do j=1,nv*ncell
         i = cell(j,1)
         wk(i) = wk(i) + 1
      enddo
c  wk() is an array refering to each vertex, and holds the num of cells
c     Step 2:  Generate csr pointer

      ic(1) = 1
      do i=i0,i1
         ic(i+1) = ic(i) + wk(i)
         wk(i)   = 0
      enddo

c     Step 3:  fill jc()

      if (type.eq.1) then ! return cell
         do k=1,ncell
         do j=1,nv
            i = cell(j,k)
            jc(ic(i)+wk(i)) = k  ! This is the cell number
            wk(i) = wk(i) + 1    ! This is number filled for ith vertex
         enddo
         enddo
         
      else                ! return structured pointer
         do j=1,nv*ncell
            i = cell(j,1)
            jc(ic(i)+wk(i)) = j  ! This is the structured pointer
            wk(i) = wk(i) + 1    ! This is number filled for ith vertex
         enddo
      endif

c     Diagnostics

c     write(6,*)
c     do i=i0,i1
c        j0 = ic(i)
c        j1 = ic(i+1)-1
c        write(6,1) i,(jc(j),j=j0,j1)
c  1     format(i8,' c2v:',20i5)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine open_bin_file(ifbswap,nelgt,ndim,nelgv,wdsizi) 
c     open file & chk for byteswap & 8byte reals

      logical ifbswap,if_byte_swap_test

      CHARACTER*132 NAME
      CHARACTER*1  NAM1(132)
      EQUIVALENCE  (NAME,NAM1)

      integer fnami (33)
      character*132 fname,re2fle
      equivalence (fname,fnami)

      character*80 hdr
      character*5 version
      real*4      test

      character*1 re2(4)
      character*4 re24
      equivalence (re2,re24)
      DATA re24   /'.re2'       /

      common /sess/ session
      character*80 session

      integer wdsizi

      call izero  (fnami,33)

      len = ltrunc(session,80)
      call chcopy (nam1,session,80)
      len = len + 1
      call chcopy (nam1(len),re2,4)
      len = len + 3
      call chcopy (fname,nam1,len)

      ierr=0
      call byte_open(fname,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error opening file in open_bin_file ',ierr)
      call byte_read(hdr,20,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error reading header in open_bin_file ',ierr)
c      write(6,80) hdr
c   80 format(a80)

      read (hdr,1) version,nelgt,ndim,nelgv
    1 format(a5,i9,i3,i9)
      wdsizi=4
      if(version.eq.'#v002')wdsizi=8
      if(version.eq.'#v003')wdsizi=8

      call byte_read(test,1,ierr)
      if(ierr.ne.0) call exitti
     $  ('Error reading test number in open_bin_file ',ierr)
      ifbswap = if_byte_swap_test(test)

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
 
