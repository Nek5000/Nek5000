c-----------------------------------------------------------------------
c
c  This code no longer relies on Metis to do the partitioning.
c
c  For large problems ( nel > 1e6 ), compile with -mcmodel=medium and
c  change parameter lelm (maximum number of elements)
c
c
c  genmap() uses the symmetric vertex ordering
c
c      3 ----- 4
c      |       |
c      |       |
c      |       |
c      1 ----- 2
c
c  which can be extended to an arbitrary number of space dimensions
c  (like 3).
c
c
c
c  7/27/07 -- add self-connected check
c  7/27/07 -- verify that don't double-check periodic faces
c
c  1/13/08 -- added support for conjugate heat transfer;
c             Main features:
c
c             .recursive bisection only on fluid graph
c             .remaining solid elements distributed round-robin/greedy
c
c             Will need to later modify for case nel_mhd > nelv.
c
c-----------------------------------------------------------------------
      program genmap

c     read nekton .rea file and make a .map file


      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)
      common /carrayi/ cell (lpts) , pmap (lpts)
     $               , order(lpts) , elist(lpts)
      common /carrayw/ w1   (lpts) , w2   (lpts)
     $               , w3   (lpts) , w4   (lpts)
     $               , w5   (lpts)
      integer     w13(3*lpts)
      equivalence (w1,w13)

      common /arrayr/  dx(4*lpts)

      common /carrayr/ bc(5*6*lelm)
      common /carrayc/ cbc(6,lelm)
      character*3      cbc

      logical ifconn,is_connected

      integer     cell,pmap,order,elist,w1,w2,w3,w4,depth

      call makemesh  (cell,nelv,nelt,irnk,dx,cbc,bc,ndim)
c                                    irnk is # unique points


      nfc = 2*ndim
      nv  = 2**ndim

      nic = lpts
      njc = lpts
      itype = 2   ! return structured vertex (not cell) pointer
      call cell2v(i0,i1,w1,nic,w2,njc,cell,nv,nelt,itype,w3) 
c     call out_cell(cell,nv,nelt)
c     call exitt

c     Determine number of outflow points and order them last
      call izero (order,irnk)
      call set_outflow
     $     (no,order,mo,cell,nv,nelt,irnk,cbc,nfc,w1,i0,i1,w2)
c     call out_cell(cell,nv,nelt)


c     Find all periodic connections, based on cbc info.
      call periodic_vtx
     $               (cell,nv,nelt,irnk,dx,ndim,cbc,bc,nfc,w13,w4)
c     call out_cell(cell,nv,nelt)
c     call exitt(1)

c     Recursive bisection of element graph; reverse-order interface points
      call rec_bisect (elist,pmap,order,mo,cell,nv,nelv,ndim
     $                                              ,w1,w2,w3,w4,w5)

c     Clean up 
      call isort     (elist,w1,nelv)
      call iswap_ip  (pmap ,w1,nelv)


      if (nelt.gt.nelv) then
         itype = 2   ! return structured vertex (not cell) pointer
         call cell2v(i0,i1,w1,nic,w2,njc,cell,nv,nelt,itype,w3) 

         call greedy  (elist,pmap,cell,nv,nelv,nelt,ndim
     $                                              ,w1,w2,w3,w4,w5)

c        Clean up 
         call isort     (elist,w1,nelt)
         call iswap_ip  (pmap ,w1,nelt)

      endif



      npts = nv*nelt
      call iranku       (cell,nrnk,npts,w1)
      call self_chk     (cell,nv,nelt,1)    ! check for not self-ptg.

      call fill_order   (order,mo,cell,nv,nelt)
      call assign_order (cell,nv,nelt,order)


      call iranku       (cell,nrnk,npts,w1) ! make cell numbering contiguous
      call self_chk     (cell,nv,nelt,2)    ! check for not self-ptg.
      call reverse_p    (pmap,nelt)         ! lightly load node 0

c     Output to .map file:
      noutflow    = no    ! for now - no outflow bcs
      call out_mapfile (pmap,nelt,cell,nv,nrnk,noutflow)

c     call out_geofile (dx,ndim,nv,nelt,pmap,39)
c     call out_geofile2(dx,ndim,nv,nelt,cell,nrnk)

c     call outmati(pmap,13,9,'pmap  ',nelt,1)
c     open(unit=22,file='p.dat')
c     write(22,1) (pmap(k),k=1,nelt)
c   1 format(i9)
c     close(unit=22)

      stop
      end
c-----------------------------------------------------------------------
      subroutine makemesh(cell,nelv,nelt,irnk,dx,cbc,bc,ndim)

c     read nekton .rea file and make a mesh

      integer      cell(1)
      character*3  cbc (6,1)
      real         bc  (5,6,1)
      real         dx  (1)
      integer e,f

      character*3 cbt(6)
      real        bt(5,6)

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)

      common /arrayi/ i_n(lpts) , j_n(4*lpts)
     $              , i_o(lpts) , j_o(4*lpts)

      logical ifbinary,ifbswap
      integer buf(30)

      integer eface(6)  ! return Nekton preprocessor face ordering
      save    eface
      data    eface / 4 , 2 , 1 , 3 , 5 , 6 /
         
      call getfile2('Input (.rea) file name:$','.rea$',10)
      call cscan_dxyz (dx,nelt,nelv,ndim,ifbinary,ifbswap)

      if (ifbinary) then
         ! skip curved side data
         call byte_read(ncurve,1)
         if (ifbswap) call byte_reverse(ncurve,1)
         do k = 1,ncurve
            call byte_read(buf,8)
         enddo

c        For current version of genmap, only need the fluid bcs.
c        Later, for more complex MHD configs, we'll need fluid + induct.

c        Also, if we start to have fluid/thermal domains with differing
c        periodic bc topologies, then we'll need something completely
c        different (two distinct numberings).
c
c        In fact, we need the thermal bcs because of periodicity... 
c        This is not pretty --- there are many cases to be considered.
c
c        For now, we default to solid topology if nelt > nelv
c


         call rd_bc_bin(cbc,bc,nelv,nelt,ifbswap)

         call byte_close()
      else
         call cscan_bcs   (cbc,bc,nelv,nelt,ndim)
      endif
c     call outbc(cbc,bc,nelt,ndim,' CBC 1')
      close (unit=10)  ! close .rea file

      nface = 2*ndim
      do e=1,nelt !  SWAP TO PREPROCESSOR NOTATION
         call chcopy(cbt,cbc(1,e)  ,nface*3)
         call copy  ( bt, bc(1,1,e),nface*5)
         do f=1,nface
            call copy  ( bc(1,f,e), bt(1,eface(f)),5)
            call chcopy(cbc(  f,e),cbt(  eface(f)),3)
         enddo
      enddo
c     call outbc(cbc,bc,nelt,ndim,' CBC 2')


c     Compress vertices based on coordinates
      call unique_vertex2(cell,dx,ndim,nelt,i_n,j_n,j_o)

      nv   = 2**ndim
      npts = nelt*nv
      call iranku    (cell,irnk,npts,i_n)
      call self_chk  (cell,nv,nelt,3)       ! check for not self-ptg.

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt(ie)
      write(6,*)
      write(6,*) ie,' quit'
      ke = 2*ie
      ff = 1./(ke-ie-ie)
      stop
      end
c-----------------------------------------------------------------------
      subroutine cscan_dxyz (dx,nelt,nelv,ndim,ifbinary,ifbswap)

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !

c
c     Scan for xyz data, read it, and set characteristic length, d2
c
      character*80 string
c
      real dx(1)
      real x(8),y(8),z(8)
      integer e,buf(30)

      integer h2s(8) ! hypercube to strange ordering
      save    h2s
      data    h2s / 1,2,4,3,5,6,8,7 /

      logical ifbinary,ifbswap


      ifbinary = .false.
      ifbswap  = .false.

      call cscan(string,'MESH DATA',9)
      read (10,*) nelt,ndim,nelv
       
      if (nelt.lt.0) then
         ifbinary = .true.
         nelt = abs(nelt)
         call open_bin_file(ifbswap)
         nwds = 1 + ndim*(2**ndim) ! group + 2x4 for 2d, 3x8 for 3d
      endif

      write(6,*) nelt,ndim,nelv,ifbinary, ' nelt,ndim,nelv,ifre2 '

      if (nelt.gt.lelm) then 
         write(6,*) 'ABORT: NELT>LELM, modify LELM and recompile'
         call exitt(1)
      endif

      b = 1.e22
      l = 1
      if (ndim.eq.3) then
         do e=1,nelt
            if(ifbinary) then
              call byte_read(buf,nwds)
              call buf_to_xyz(buf,x,y,z,e,ifbswap,ndim)
            else 
              read (10,80) string
              read (10,*)   (x(k),k=1,4)
              read (10,*)   (y(k),k=1,4)
              read (10,*)   (z(k),k=1,4)
              read (10,*)   (x(k),k=5,8)
              read (10,*)   (y(k),k=5,8)
              read (10,*)   (z(k),k=5,8)
            endif
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
              call byte_read(buf,nwds)
              call buf_to_xyz(buf,x,y,z,e,ifbswap,ndim)
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
c
      nvrt = 2**ndim
      call set_d2(dx,nvrt,nelt,ndim)
c
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
      real         bc (5,6,nelt)
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
         if (ipass.eq.2) nel=nelt
         call rd_bc(cbc,bc,nel,ndim,ifield,10)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rd_bc(cbc,bc,nel,ndim,ifield,io)

c     .Read Boundary Conditions (and connectivity data)

      character*3 cbc(6,nel)
      real        bc(5,6,nel)
      integer e,f
C
      nbcrea = 5
      nface  = 2*ndim
      do e=1,nel
      do f=1,nface
         if (nel.lt.1000) then
            read(io,50,err=500,end=500)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
   50       format(a1,a3,2i3,5g14.7)
         elseif (nel.lt.100 000) then
            read(io,51,err=500,end=500)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
   51       format(a1,a3,i5,i1,5g14.7)
         else
            read(io,*,err=500,end=500)    
     $      cbc(f,e),id1,(bc(ii,f,e),ii=1,nbcrea)
         endif
c        write(6,*) e,f,' ',cbc(f,e),' BC IN?'
      enddo
      enddo

      return
C
C     Error handling:
C
  500 continue
      write(6,501) ifield,e
  501 FORMAT(2X,'ERROR READING BOUNDARY CONDITIONS FOR FIELD',I4,I6
     $    ,/,2X,'ABORTING IN ROUTINE RDBDRY.')
      call exitt(ifield)
      return
C
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
      subroutine getfile2(prompt,suffix,io)
c
      character*1 prompt(1),suffix(1)
c
      common /sess/ session
      character*80 session

      character*80 file
      character*1  file1(80)
      equivalence (file1,file)
      character*80 fout
      character*1  fout1(80)
      equivalence (fout1,fout)
c
c     Get file name

      len = indx1(prompt,'$',1) - 1
      write(6,81) (prompt(k),k=1,len)
   81 format(80a1)

      call blank(session,80)
      read(5,80) session
   80 format(a80)
      call chcopy(file,session,80)
      len = ltrunc(file,80)
      lsf = indx1 (suffix,'$',1) - 1
      call chcopy(file1(len+1),suffix,lsf)

      open(unit=io, file=file)

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
c
c     Use Heap Sort (p 233 Num. Rec.), 5/26/93 pff.
c
c
      integer a(1),p(1)
      integer rank

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
      subroutine part_dual(pmap,nmap,c,nv,nel,irnk)
c
c     Patition map into two subdomains, based on dual of mesh graph
c

      integer pmap(1),nmap(1),c(nv,nel)
      integer etype,edgecut,e

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)
      parameter(mm  =50)

      common /arrayr1/ f(lelm),r(lelm),p(lelm),w(lelm),rr(lelm,mm)
     $               , ev(mm*mm),d(mm),u(mm)

      common /arrayi2/ jdual(lpts) , vdual(lpts)
     $               , idual(lelm) 
      common /arrayi3/ iv2c (lpts) , jv2c (lpts)
      common /arrayi4/ wk(lpts+2*lelm)
      integer jdual,vdual,idual,wk
      integer list(lpts)
      equivalence (list ,jv2c)
      logical ifconn,is_connected

      if (nv.eq.4) then   ! 2D
         ndim  = 2
         etype = 4      ! Quadrilateral
      else
         ndim  = 3
         etype = 3      ! Hexahedral
      endif


c     write(6,*) 'part dual:',nv,nel,irnk
c     call out_cell(c,nv,nel)


      nv  = 2**ndim
      call c2c(idual,jdual,vdual,ni,nj,c,nv,nel,iv2c,jv2c,lpts,wk)
c     call outaij(idual,jdual,nel,'jdual ',9)
      ifconn = is_connected(list,n0,idual,jdual,nel,wk)

      if (ifconn) then
         m = mm
         call spec_bis(pmap,idual,jdual,vdual
     $                ,nel,d,u,f,r,p,w,rr,ev,m,ndim)

      else

         write(6,*) 'not connected',n0,nel

         n1 = nel/2
         n2 = nel - n1
         i1 = 0
         i2 = 0
         do i=1,nel
            if (i1.lt.n1.and.i2.lt.n2) then
               if (list(i).ne.0) then
                  i1=i1+1
                  pmap(i)=1
               else
                  i2=i2+1
                  pmap(i)=2
               endif
            elseif (i1.lt.n1) then
               pmap(i)=1
               i1 = i1+1
            else
               pmap(i)=2
               i2 = i2+1
            endif
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rec_bisect (elist,pmap,order,mo,cell,nv,nel,ndim
     $                      ,ia,da,w1,w2,w3)

      integer pmap(1),order(1),cell(nv,nel),elist(1),ia(0:1),da(0:1)
      integer w1(1),w2(1),w3(1)
      integer e,p,v,depth

      integer pmaps(100),elists(100)

      do e=1,nel
         elist(e) = e
      enddo

      ia(0) = 1
      ia(1) = nel+1
      da(0) = 0

      max_depth = log2(nel)
      nel_even  = 2**max_depth

      if (nel.eq.1) pmap(1) = 1

      p  = 0
      l  = 1
      d  = 0
 
      depth = 0
      do i=1,nel
         j0 = ia(l-1)
         j1 = ia(l)-1
         n  = ia(l) - ia(l-1)
c        write(6,8) i,depth,max_depth,l,j0,j1,n1,n2,p,' d2   '
         if (da(l-1).lt.max_depth.and.l.gt.0) then

            da(l-1) = da(l-1) + 1
            da(l  ) = da(l-1)
            depth   = da(l-1)

            call bipart_sort
     $        (n1,n2,pmap(j0),order,mo,elist(j0),n,cell,nv,p,w1,w2,w3)
            write(6,8) i,depth,max_depth,l,j0,j1,n1,n2,p,' DEPTH'
    8       format(i9,2i5,6i8,a6)

c           write(6,18) 'A',(pmap(k),k=1,nel) 
c 18        format(a1,'pmap:',32i3)
c           call outmati(ia,1,l+1,'idepth',depth,1)
c           call outmati(da,1,l+1,'ddepth',depth,1)

            depth = da(l-1)+1
            if (depth.le.max_depth) then
               ia(l)   = j0 + n1
               ia(l+1) = j1 + 1
               l       = l+1
            else
               p       = p+2   !  increase base processor offset
               l       = l-1   !  go back in list
            endif
         endif

c        call checker(elist,pmap,nel,ndim,i)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fill_order(order,mo,cell,nv,nel)
      integer order(1),cell(nv,nel)
      integer e,v

c     Fill in remaining separator sets
      do e=1,nel
      do v=1,nv
         i = cell(v,e)
         if (order(i).eq.0) then
            mo = mo+1
            order(i) = mo
         endif
         order(i) = -abs(order(i)) ! set flag
      enddo
      enddo

c     Reverse separator set ordering
      do e=1,nel
      do v=1,nv
         i = cell(v,e)
         if (order(i).lt.0) order(i) = mo+1+order(i)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ident_sep( order, mo, elist, cell, nv, n1, n2 )
c
c     order separator nodes
c
c
c     Output:
c
c     order - updated vertex separator ordering
c     mo    - current max(order)
c
c     Input:
c
c     elist - list of active elements
c     cell  - list of vertices for each element (in global addr. space)
c     mo    - current max(order)
c

      integer elist(1),cell(nv,1),order(1)
      integer e,v

c--- diagnostic use only -----------------
                                         !
      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE
      common /arrayr/  dx(0:3,8,lelm)    !
                                         !
      integer icalld,kj(10000),ke(10000) !
      save    icalld                     !
      data    icalld /0/                 !
                                         !
                                         !
c--- diagnostic use only -----------------

      mod = mo

      do k=1,n1         ! Set flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.0) order(i) = -1
         enddo
      enddo

      m0=0
      do k=n1+1,n1+n2         ! Check flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-1) then
                mo = mo+1
                order(i) = mo
                m0 = m0+1
                kj(m0) = j
                ke(m0) = e
            endif
         enddo
      enddo

      nsep = mo-mod
      if (nsep.gt.0) icalld=icalld+1

c     write(6,6) nsep,mo,mod,icalld
c   6 format(4i11,' nsep')

      do i=1,nsep
         x=dx(1,kj(i),ke(i))
         y=dx(2,kj(i),ke(i))
         z=dx(3,kj(i),ke(i))
         write(9,4) x,y,z,icalld
      enddo
    4 format(1p3e12.4,i9)

      do k=1,n1         ! Unset flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-1) order(i) = 0
         enddo
      enddo

      call out_order(order,mo,elist,cell,nv,n1,n2)
c     call exitt(5)


      return
      end
c-----------------------------------------------------------------------
      subroutine count_sep( count, order, nsep, elist, cell, nv, n1, n2)
c
c     order separator nodes
c
c
c     Output:
c
c     count - number of separators on each element
c
c     Input:
c
c     elist - list of active elements
c     cell  - list of vertices for each element (in global addr. space)
c

      integer count(1),elist(1),cell(nv,1),order(1)
      integer e,v

      nel = n1+n2
      call izero (count,nel)

      do k=1,n1         ! Set flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.0) order(i) = -1
         enddo
      enddo

      do k=n1+1,nel           ! Check flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-1) order(i) = -2
         enddo
      enddo

      do k=1,nel              ! Check flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-2) count(k) = count(k)+1
         enddo
      enddo

      nsep = 0
      do k=1,n1         ! Unset flags and count separators
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-2) nsep = nsep+1
            if (order(i).lt.0)  order(i) = 0
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bipart_sort
     $   (n1,n2,pmap,order,mo,elist,nel,cell,nv,p,c,w1,w2)
c
c     Bipartition elist, order separator nodes, add p to partition 
c
c
c     Output:
c
c     n1,n2 - number of elements on each 1/2
c
c     pmap  - sorted processor map
c     order - updated vertex separator ordering
c     mo    - current max(order)
c
c     elist - sorted by processor map
c
c     Input:
c
c     elist - list of active elements
c     nel   - number of active elements
c     cell  - list of vertices for each element (in global addr. space)
c     nv    - 8 for 3D, 4 for 2D
c     p     - current processor offset
c
c     Work arrays:
c
c     c     - nv*nel
c     w1    - nv*nel
c     w2    - nv*nel
c

      integer pmap(nel),order(1),elist(nel),cell(nv,1),p,c(nv,1)
      integer w1(1),w2(1),wk2(2)
      integer e,key(2)

      call isort (elist,w1,nel)
      do k=1,nel   ! extract cell sublist
         e = elist(k)
         call icopy(c(1,k),cell(1,e),nv)
      enddo

      npts = nv*nel
      call iranku   (c,irnk,npts,w1)
      call self_chk (c,nv,nel,4)       ! check for not self-ptg.

      call part_dual(pmap,w1,c,nv,nel,irnk)  ! pmap contains processor map

      k = 1
      do i=1,nel
         c(k  ,1) = pmap (i)
         c(k+1,1) = elist(i)
         k = k+2
      enddo
      key(1) = 1
      key(2) = 2
      nkey   = 2
      call ituple_sort(c,2,nel,key,nkey,w1,wk2)
      k = 1
      do i=1,nel
         pmap (i) = c(k  ,1)
         elist(i) = c(k+1,1)
         k = k+2
      enddo

      n1 = nel     ! Count number in each 1/2
      do i=2,nel
         if (pmap(i).ne.pmap(1)) then
            n1 = i-1
            goto 10
         endif
      enddo
   10 continue
      n2 = nel - n1
c     write(6,*) n1,n2,' n1,n2 raw'

      if (abs(n2-n1).gt.1) then ! rebalance load
         if (n2.gt.n1) then ! swap

            j1 = n1+1
            j2 = n2+1
            call icopy(w1,elist(j1),n2)
            call icopy(elist(j2),elist,n1)
            call icopy(elist,w1,n2)

            call icopy(w1,pmap(j1),n2)
            call icopy(pmap(j2),pmap,n1)
            call icopy(pmap,w1,n2)

            m1 = n1
            n1 = n2
            n2 = m1

         endif
         ndif = (n1-n2)/2
         do i=1,ndif
            call count_sep( w1, order, nsep, elist, cell, nv, n1, n2 )
            j0 = 1
            j1 = n1+1
            call isort   ( w1    ,w2,n1)
            call iswap_ip( elist ,w2,n1)
            pmap(n1) = pmap(nel)
            n1 = n1-1
            n2 = n2+1
            mcount = iglmax(w1,nel)
c           if (mod(i,4).eq.0 .or.    i.le.3)
c    $         write(6,*) ndif,n1,n2,mcount,' Sep count',nsep
         enddo
      endif

      j0=1      ! Reset pmap and shift by current offset
      j1=n1+1
      call ifill(pmap(j0),2,n1)   ! Since n1 > n2, let's lightly load node 0
      call ifill(pmap(j1),1,n2)   ! by assigning node 0 to shorter stack.
      call icadd(pmap,p,nel)
c     maxp = iglmax(pmap,nel)
c     write(6,7) maxp,n1,n2,nel
c   7 format('bipart sort: ',4i12)


      call part_clean( order, nsep, elist, cell, nv, n1, n2, w1, w2)

      call count_sep( w1, order, nsep, elist, cell, nv, n1, n2 )

      mcount = iglmax(w1,nel)
c     write(6,*) n1,n2,mcount,' sepsep',nsep

      call ident_sep( order, mo, elist, cell, nv, n1, n2 )

      return
      end
c-----------------------------------------------------------------------
      subroutine part_clean(order,nsep,elist,cell,nv,n1,n2,count,wk)
c
c     Clean up partition by swapping
c
      integer count(1),elist(1),cell(nv,1),order(1),wk(1)
      integer e,v

      nel  = n1+n2
      snel = nel
      snel = sqrt(snel)
      nels = snel + 5      ! check lots of options for small cases
      nels = min(nels,nel)
      
      nvmx = (3*nv)/4

      do i=1,nels

         call count_sep(count, order, nsep, elist, cell, nv, n1, n2)

         j0 = 1
         j1 = n1+1

         call isort   ( count(j0),wk,n1)
         call iswap_ip( elist(j0),wk,n1)
         m1 = iglmax  ( count(j0),n1)

         call isort   ( count(j1),wk,n2)
         call iswap_ip( elist(j1),wk,n2)
         m2 = iglmax  ( count(j1),n2)

         mx = max(m1,m2)

c        write(6,6) i,nels,n1,n2,m1,m2,nsep,' CLEAN'
         if (mx.eq.0) then
c           write(6,6) i,nels,n1,n2,m1,m2,nsep,' CLEAN'
c  6        format(7i7,a6)
c           call out_cell2(elist,cell,order,nv,nel)
c           call exitt(0)
            return
         elseif (mx.le.nvmx) then
            return
         else          ! swap
c           write(6,6) i,nels,n1,n2,m1,m2,nsep,' SWAP'
            e = elist(n1)
            elist(n1)  = elist(nel)
            elist(nel) = e
         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine out_cell2(elist,c,o,nv,nel)
c
      integer elist(1),c(nv,nel),o(1)
      integer e,etype

      etype = 3                   ! Hexahedral
      if (nv.eq.4) etype = 4      ! Quadrilateral

      open (unit=11,file='t.dat') ! Output results
      write(11,11) nel,etype
      do l=1,nel
         e = elist(l)
         write(11,11) (c(k,e),k=1,nv)
      enddo
   11 format(8i9)
      close(11)

      write(6 ,10) nel,etype,' CELL '
c     do l=1,nel
c        e = elist(l)
c        write(6 ,12) e,(c(k,e),k=1,nv)
c     enddo

      write(6 ,10) nel,etype,' SEPRT'
c     do l=1,nel
c        e = elist(l)
c        write(6,12) e,(o(c(k,e)),k=1,nv)
c     enddo

   10 format(/,2i9,a6)
   12 format(i9,'e:',8i7)

      return
      end
c-----------------------------------------------------------------------
      subroutine out_cell(c,nv,nel)

      integer c(nv,nel)
      integer e,etype

      etype = 3                   ! Hexahedral
      if (nv.eq.4) etype = 4      ! Quadrilateral

      write(6 ,10) nel,etype

      open (unit=11,file='t.dat') ! Output results
      write(11,11) nel,etype
      do e=1,nel
         write(11,11) (c(k,e),k=1,nv)
         write(6 ,12) e,(c(k,e),k=1,nv)
      enddo

c     write(6,*) 'part_dual cont?'
c     read (5,*) xx

   10 format(/,2i9,' CELL')
   11 format(8i9)
   12 format(i9,'e:',8i7)
      close(11)

      return
      end
c-----------------------------------------------------------------------
      subroutine outmati(u,m,n,name6,nid,ic)
      integer u(m,n)
      character*6 name6
      character*1 adum
c
c     return
c
c     Print out copies of a global matrix
c
      if (m.gt.1) then
         write(6,1) nid,m,n,name6
   1     format(3i6,'  Matrix:',2x,a6)
         do i=1,m
            write(6,2) i,name6,(u(i,j),j=1,n)
         enddo
   2     format(i8,1x,a6,20(10i9,/,10x))
      else
         write(6,3) nid,n,name6,(u(1,j),j=1,n)
   3     format(2i8,1x,a6,20(10i9,/,10x))
      endif
      if (ic.eq.0) then
         write(6,*) 'cont: ',name6,nid,'  ??'
c        read (5,*) adum
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat (u,m,n,name6,nid)
      real u(m,n)
      character*6 name6
c
c     return
c
c     Print out copies of a global matrix
c
      write(6,1) nid,m,n,name6
      n8 = min(8,n)
   1  format(3i6,'  Matrix:',2x,a6)
      do i=1,m
         write(6,2) i,name6,(u(i,j),j=1,n8)
      enddo
   2  format(i3,1x,a6,1p8e12.4)
c  2  format(i3,1x,a6,20(1p8e12.4,/,10x))
      return
      end
c-----------------------------------------------------------------------
      subroutine out_mapfile (pmap,nel,cell,nv,nrnk,noutflow)
      integer pmap(nel),cell(nv,nel)
      integer depth,d2

      depth            = log2(nel)
      d2               = 2**depth
      npts             = nel*nv
      call dmp_mapfile (pmap,nel,depth,cell,nv,nrnk,npts,noutflow)

      return
      end
c-----------------------------------------------------------------------
      subroutine dmp_mapfile (pmap,nel,depth,cell,nv,nrnk,npts,noutflow)

      common /sess/ session
      character*80 session
      character*80 fname
      character*1  fnam1(80)
      equivalence (fnam1,fname)

      integer pmap(nel),depth,cell(nv,nel)
      integer d2,e,p0

      d2 = 2**depth

      write(6,*) 'DEPTH:',depth,d2,nel,nrnk,npts,noutflow

      nactive = nrnk - noutflow

      len = ltrunc(session,80)
      call chcopy(fname,session,80)
      call chcopy(fnam1(len+1),'.map',4)
      open (unit=29,file=fname)

      write(29,1) nel,nactive,depth,d2,npts,nrnk,noutflow
    1 format(9i9)

      do e=1,nel
         p0 = pmap(e)-1
         write(29,1) p0,(cell(k,e),k=1,nv)
      enddo

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
      subroutine set_outflow(no,order,mo,cell,nv,nel,nrnk,cbc,nfc
     $                      ,ic,i0,i1,jc)
c
c     Order outflow nodes last
c
      integer cell(nv,nel),order(1)
      character*3      cbc(6,nel)

      integer ic(i0:i1),jc(1)

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)
      common /arrayi2/ face (3*lpts) , elist(lelm) , ind  (lpts)
      integer face,elist

      integer     e,f,out_vtx,out_vtm
      character*3 cb
      logical     ifoutflow

      integer vface(4,6)  ! symm. vertices ordered on symm. faces
      save    vface
      data    vface / 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8
     $              , 1,2,3,4 , 5,6,7,8 /

      ifoutflow = .false.

      nvf  = nv/2
      mvtx = nel*nv

      do e=1,nel
      do f=1,nfc
         cb = cbc(f,e)
c        write(6,*) cb,e,f,' cb'
         if (cb.eq.'O  ') then
            do i=1,nvf
               j = vface(i,f)
               cell(j,e) = cell(j,e) + mvtx
            enddo
            ifoutflow = .true.
         endif
      enddo
      enddo

c     Make cells consistent if you have an 'O  ' next to 'W  ' (say)
      do j=1,nv*nel
         i = cell(j,1)
         if (i.gt.mvtx) then  ! Outflow node, make attached nodes same
            k0 = ic(i)
            k1 = ic(i+1)-1
            do k=k0,k1
               jj = jc(k)
               cell(jj,1) = i
            enddo
         endif
      enddo

      mo = 0
      no = 0
      write(6,*) no,mo,out_vtm,ifoutflow,' OUTFLOW'
      if (.not. ifoutflow) return

      npts             = nel*nv
      call iranku      (cell,nrnk,npts,ind)
      call self_chk    (cell,nv,nel,5)       ! check for not self-ptg.

      out_vtm = mvtx + 1
      do e=1,nel  ! Determine number of unique outflow pts
      do f=1,nfc
         cb = cbc(f,e)
         if (cb.eq.'O  ') then
            do i=1,nvf
               j = vface(i,f)
               out_vtx = cell(j,e)
               out_vtm = min(out_vtm,out_vtx)
            enddo
         endif
      enddo
      enddo

      no = nrnk - out_vtm + 1

      mo = 0
      do i=out_vtm,nrnk
         mo = mo+1
         order(i) = mo
      enddo
      write(6,*) no,mo,out_vtm,ifoutflow,' OUTFLOW'

      return
      end
c-----------------------------------------------------------------------
      subroutine assign_order(cell,nv,nel,order)
c
c     Order nodes by "order"
c
      integer cell(nv,nel),order(1)

      integer     e,f

      do e=1,nel
      do k=1,nv
         i = cell(k,e)
         j = order(i)
         cell(k,e) = j
      enddo
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
      real             bc (5,6,nel)

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

      nvf = nv/2                   ! # vertices/face = 1/2 # vertices/cell

      call izero(iper,ndim*irnk)   ! Zero out periodic indicator flag

      call jjnt (jmin,irnk)        ! Initial permutation = identity

      nmn = irnk
      nmx = 0
      do e=1,nel
      do f=1,nfc
         cb = cbc(f,e)
c        write(6,*) cb,e,f,' cb'
         if (cb.eq.'P  ') then
           je = abs(bc(1,f,e))
           jf = bc(2,f,e)
           jf = efaci(jf)

           cj = cbc(jf,je)
           ke = abs(bc(1,jf,je))
           kf = bc(2,jf,je)
           kf = efaci(kf)
            
           if (bc(1,f,e).gt.0 .and. bc(1,jf,je).gt.0) then
              if (ke.ne.e .or. kf.ne.f .or. cj.ne.'P  ') then
               write(6,*)
               write(6,*) 'abort: PERIODIC MISMATCH 1:'
               write(6,6)   e,f,cb,' ie '
               write(6,6) je,jf,cj,' je '
               write(6,6) ke,kf,cj,' ke '
               write(6,*)
    6          format(i9,i3,1x,a3,1x,a4)
               call exitt(9)
              endif

              call find_connctd_pairs
     $                 (jmin,nvf,e,f,je,jf,cell,nv,dx,ndim)
            

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

c
c     Okay -- we now have the updated pointers, time to update
c     the cell pointers
c
      npts = nv*nel
      do i=1,npts
         cell(i,1) = jmin(cell(i,1))  ! permuted identity
      enddo
      call iranku      (cell,irnk,npts,iper) ! compress cell list
      call self_chk    (cell,nv,nel,6)       ! check for not self-ptg.


c
c     Reset bc array
c
      do e=1,nel
      do f=1,nfc
         if (cbc(f,e).eq.'P  ') bc(1,f,e) = abs(bc(1,f,e))
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine find_connctd_pairs(jmin,nvf,e,f,je,jf,cell,nv,dx,ndim)


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
            x0m = max(x0m,x0(k,i))
            x0m = max(x0m,x0(1,i))
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
         call outmat(x0,4,4,'  x0  ', e)
         call outmat(x1,4,4,'  x1  ',je)
         call outmat(z0,4,4,'  z0  ', e)
         call outmat(z1,4,4,'  z1  ',je)
         write(6,6) e , f,shift,eps,x0m
         write(6,6) je,jf,    i,tol,d2min
   6     format(i8,i2,i3,1p2e16.8,' abort: FACE MATCH FAIL')
         call exitt(0)
      endif

      write(6,7) e,f,i,shift,d2min
   7  format(i8,i2,2i3,1p1e16.8,' shift')

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
      subroutine out_order(order,mo,elist,cell,nv,n1,n2)
c
c
c         +-----+
c         |     |
c         |     |
c         +-----+
c
c
c
c
c     Output orders for an 8x8 array of elements

      integer elist(1),cell(nv,16,16),order(1)
      integer e,v
      character*1 a(7,4,16,16)

      integer h2s(8) ! hypercube to strange ordering
      save    h2s
      data    h2s / 1,2,4,3,5,6,8,7 /

      return

      call set_a(a)

      do l=1,n1+n2

         e = elist(l)
         if (l.le.n1) a(4,3,e,1) = 'X'
         if (l.gt.n1) a(4,3,e,1) = 'O'

         k=0
         do j=0,1
         do i=0,1
            k = k+1
            kstupid = h2s(k)
            v = cell(kstupid,e,1)
            if (i.eq.0.and.j.eq.0) write(a(1,1,e,1),1) order(v)
            if (i.eq.1.and.j.eq.0) write(a(5,1,e,1),1) order(v)
            if (i.eq.0.and.j.eq.1) write(a(1,4,e,1),1) order(v)
            if (i.eq.1.and.j.eq.1) write(a(5,4,e,1),1) order(v)
   1        format(i2)
         enddo
         enddo

      enddo
      write(6,*) 'n12:',n1,n2,nv
      call out_a(a)

      return
      end
c-----------------------------------------------------------------------
      subroutine out_a(a)
      character*1 a(7,4,16,16)
c
      do je=8,1,-1
      do j =4,1,-1
         write(6,1) ((a(i,j,ie,je),i=1,7),ie=1,8)
      enddo
      enddo
   1  format(56a1)

      write(6,*)
      write(6,*) 'continue ?'
      read (5,*) dumm
      
      return
      end
c-----------------------------------------------------------------------
      subroutine set_a(a)
      character*1 a(7,4,16,16)
c
c         +-----+
c         |     |
c         |     |
c         +-----+
c
      call blank(a,7*4*64)
c
      do je=1,8
      do ie=1,8
         write(a(1,1,ie,je),1)
         write(a(1,2,ie,je),2)
         write(a(1,3,ie,je),2)
         write(a(1,4,ie,je),1)
      enddo
      enddo
   1  format('+-----+')
   2  format('|     |')

      return
      end
c-----------------------------------------------------------------------
      subroutine outaij(ia,ja,n,name6,key)

      integer ia(0:n),ja(1)
      character*6 name6

c     if ia     < 0, then already on list
c     if ja(ia) < 0, then row is already processed

      write(6,*)
      write(6,*) ' OUT Aij ',name6,n,key

      do i=1,n
         j0 = abs(ia(i-1))
         j1 = abs(ia(i))-1
         m  = j1-j0 + 1
         m  = min(16,m)
         jm = j0 + m - 1
         write(6,1) i,j0,j1,'aij:',(ja(j),j=j0,jm)
      enddo
   1  format(3i4,1x,a4,1x,16i4)

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
c     call outmati(stack,1,n+1,'stack:',val,1)

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
c     call outmati(stack,1,n+1,'popst:',ipop,1)

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
      subroutine spec_bis(pmap,ia,ja,va,n,d,u,f,r,p,w,rr,ev,m,ndim)

c     n = dimension of A
c     m = max # iterations

      real d(m),u(m),f(n),r(n),p(n),w(n),rr(n,m)
      integer pmap(n),ia(0:n),ja(1),va(1)

      if (n.lt.3) then
         call sbisect (pmap,f,p,w,n)
         return
      endif

      call rint(f,n)
      n2 = n/2
      r1 = 1000.*n
      do i=1,n2              ! bias in 0-0-1 direction
         f(i) = f(i) + r1
      enddo
      call ortho1(f,n)
      ftf = glsc2(f,f,n)
      fnm = 1./sqrt(ftf)
      call cmult(f,fnm,n)
      
c     rn = n
c     rd = rn**(1./ndim)

c     one = 1.               ! For diagnostics only
c     pi  = 4.*atan(one)
c     i   = 0
c     do jj=1,8
c     do ii=1,8
c        i = i+1
c        a = pi*(ii-1)/7
c        f(i) = cos(a)
c     enddo
c     enddo

      npass = 50
      do k=1,npass
         niter = m
         call glanczos(rr,n,d,u,niter,f,ia,ja,va,r,p,w)
         call lanczos2(f,rr,n,ev,d,u,niter)
         if (niter.lt.m) goto 100
      enddo
  100 continue

      call sbisect (pmap,f,p,w,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine sbisect(pmap,f,p,w,n)
      real f(n)
      integer pmap(n),p(n),w(n)

      if (n.gt.2) then
         call sort     (f,w,n)
         call jjnt     (p,n)
         call iswap_ip (p,w,n)
      else
         call jjnt     (p,n)
      endif

      n2 = n/2
      if (n.eq.1) n2 = 1

      do i=1,n2
         pmap(p(i)) = 1
      enddo

      do i=n2+1,n
         pmap(p(i)) = 2
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine glanczos(rr,n,diag,upper,niter,f,ia,ja,va,r,p,w)
c
c     Lanczos applied to graph Laplacian
c
      real    rr(n,1),diag(1),upper(1),f(1),r(1),p(1),w(1)
      integer ia(1),ja(1),va(1)
      real one,eps
c
      call rzero(diag ,niter)
      call rzero(upper,niter)
      pap = 0.0
c
c     set machine tolerances
c
      one = 1.
      eps = 1.e-28
      if (one+eps .eq. one) eps = 1.e-13
      if (one+eps .eq. one) eps = 1.e-5
      eps = 1.e-5
c
      rtz1=1.0
c
      call copy  (r ,f,n)
      call ortho1(r,n)
      rtr   = glsc2(r,r,n)
      rnorm = sqrt(rtr)
      rtol  = rnorm*eps
      rni   = 1./rnorm
      call cmult2  (rr,r,rni,n)
      iter = 0
c
      do 1000 iter=1,niter

         rtz2=rtz1
         rtz1=rtr

         beta = rtz1/rtz2
         if (iter.eq.1) beta=0.0

         call add2s1(p,r,beta,n)
         call ortho1(p,n)
         call ax(w,p,ia,ja,va,n)
c
c        Save p^Ap for eigenvalue estimates
         pap_old = pap
         pap=glsc2(w,p,n)
         alpha=rtz1/pap
         alphm=-alpha
         call add2s2(r,w,alphm,n)
c        call outmat (r,16,16,'resid ',iter)
c
         rtr = glsc2(r,r,n)
         if (rtr.le.0) goto 1001
         rnorm = sqrt(rtr)
         rni   = 1./rnorm
         call cmult2  (rr(1,iter+1),r,rni,n)
c        write(6 ,6) iter,n,rnorm,rtol,alpha,beta,pap
c        write(60,6) iter,n,rnorm,rtol,alpha,beta,pap
c        if (iter.le.10.or.mod(iter,10).eq.0) write(6,6) iter,n,rnorm
c
c        Generate tridiagonal matrix for Lanczos scheme 
         if (iter.eq.1) then
            diag(iter) = pap/rtz1
         else
            diag(iter)    = (beta**2 * pap_old + pap ) / rtz1
            upper(iter-1) = -beta * pap_old / sqrt(rtz2 * rtz1)
         endif
         if (rnorm.le.rtol) goto 1001
 1000 continue
      iter = iter-1
 1001 continue
c
      niter = iter
      write(6,6) iter,n,rnorm,rtol
    6 format(i4,i10,' cg:',1p6e12.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ax(y,x,ia,ja,va,n)

c     This routine computes y = Ax, where A is the graph Laplacian 

      real y(1),x(1)
      integer ia(0:1),ja(1),va(1)
c
      do i=1,n
         j0 = ia(i-1)
         j1 = ia(i) - 1
         m  = j1-j0 + 1
         y(i) = 0
         do j=j0,j1
            y(i) = y(i) + va(j)*x(ja(j))
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine ortho1(p,n)
c
c     Orthogonalize wrt the 1 vector
c
      real p(1)
c
      s = 0.
      do i=1,n
         s = s + p(i)
      enddo
      s = s/n
      do i=1,n
         p(i) = p(i) - s
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lanczos2(v1,rr,n,ev,diag,upper,m)
c
c     Compute two smallest eigenvector estimates 
c
      real    v1(1),rr(n,1),ev(m,1),diag(1),upper(1)

      call calcvec (diag,upper,v1,v1(m+1),m,imin,imax,ev)
      call rzero(v1,n)
      do j=1,m
         call add2s2(v1,rr(1,j),ev(j,imin),n)
      enddo

c     call outmat(ev(1,imin),1,m,'evmin ',imin)
c     call outmat(v1,16,16,'rsbeig',imin)

      return
      end
c-----------------------------------------------------------------------
      subroutine calcvec (diag,upper,d,e,n,imin,imax,z)
c
c     This is the numerical recipes tridiagonal eigenvalue/vector solver.
c
c     -  Eigen vectors are returned in z(n,n)
c
c
      real diag(n),upper(n)
      real d(n),e(n),z(n,n)
c
      call copy (d,diag ,n)
      call copy (e,upper,n)
      call ident(z,n)
c
      do 15 l=1,n
         iter = 0
c
    1    do 12 m=l,n-1
            dd = abs( d(m) ) + abs( d(m+1) )
            if ( abs(e(m)) + dd .eq. dd ) goto 2
   12    continue
c
         m = n
    2    if ( m .ne. l ) then
c
            if ( iter .eq. 30 ) then
               write (6,*) 'too many iterations'
               return
            endif
c
            iter = iter + 1
            g = ( d(l+1) - d(l) ) / ( 2.0 * e(l) )
            r = sqrt( g**2 + 1.0 )
c
c    sign is defined as a(2) * abs( a(1) )
c
            g = d(m) - d(l) + e(l)/(g+sign(r,g))
            s = 1.0
            c = 1.0
            p = 0.0
c
            do 14 i = m-1,l,-1
               f = s * e(i)
               b = c * e(i)
               if ( abs(f) .ge. abs(g) ) then
                  c = g/f
                  r = sqrt( c**2 + 1.0 )
                  e(i+1) = f*r
                  s = 1.0/r
                  c = c*s
               else
                  s = f/g
                  r = sqrt( s**2 + 1.0 )
                  e(i+1) = g*r
                  c = 1.0 / r
                  s = s * c
               endif
c
               g = d(i+1) - p
               r = ( d(i) - g ) * s + 2.0 * c * b
               p = s * r
               d(i+1) = g + p
               g = c*r - b
c      ...     find eigenvectors ... (new, 11/19/94, pff, p.363. Num.Rec.I.)
               do 13 k=1,n
                  f = z(k,i+1)
                  z(k,i+1) = s*z(k,i)+c*f
                  z(k,i  ) = c*z(k,i)-s*f
   13          continue
c      ...     end of eigenvector section ... 
   14       continue
c
            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0
            goto 1
         endif
c
   15 continue
c
      imin = 1
      imax = 1
      dmin = d(imin)
      dmax = d(imax)
      do 40 i = 1 , n
        if (d(i).lt.dmin) then
           dmin = d(i)
           imin = i
        endif
        if (d(i).gt.dmax) then
           dmax = d(i)
           imax = i
        endif
c       write(6,*) i,imin,imax,d(i),' iminx'
   40 continue
c     write(6,41) 'eig:',(d(i),i=1,n)
   41 format(a4,2x,5(1p10g12.4,/))
c
c     Orthonormalize eigenvectors
c
      n10 = min(n,10)
      do ko=1,n
         do ki=1,n
            e(ki) = glsc2(z(1,ki),z(1,ko),n)
            if (e(ki).ne.0.0) e(ki) = sqrt(abs(e(ki)))
         enddo
c        write(6,9) d(ko),(e(ki),ki=1,n10)
c   9    format(1pe12.4,' e:',1p10e12.4)
         scale = 1.0/e(ko)
         call cmult(z(1,ko),scale,n)
      enddo
c
      return
      end
c=======================================================================
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
      write(6,*) 'pmax:',pmax,n

      do i=1,n
         p(i) = pmax+1 - p(i)  ! range of p is [1:pmax]
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine unique_vertex2(cell,dx,ndim,nel,ind,ninseg,ifseg)
c
      integer cell(1),ind(1),ninseg(1)
      logical ifseg(1)
      real dx(0:ndim,nel)

      integer e
      real dxt(4),t1(4),t2(4)

      nvtx = 2**ndim
      n = nvtx*nel

      do i=1,n
         cell(i)   = i
         ifseg (i) = .false.
      enddo

c
c     Sort by directions
c
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
            call tuple_sort(dx(0,i),lda,ninseg(iseg),j1,1,ind,dxt) ! key=j1
            call iswap_ip  (cell(i),ind,ninseg(iseg)) ! Swap position 
            i  =   i + ninseg(iseg)
         enddo
c
c        q=0.0010   ! Smaller is better
         q=0.2      ! But generous is more forgiving for bad meshes!
         do i=2,n
           if ((dx(j,i)-dx(j,i-1))**2.gt.q*min(dx(0,i),dx(0,i-1)))
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
c     Unshuffle geometry
c
      call tuple_swapt_ip(dx,lda,n,cell,t1,t2)
c
c     Assign global node numbers (now sorted lexigraphically!)
c
      ig = 0
      do i=1,n
         if (ifseg(i)) ig=ig+1
         ind(cell(i)) = ig
      enddo
      nglb = ig
      call icopy(cell,ind,n)

      write(6,6) nseg,nglb,n
    6 format('done locglob_lexico:',3i9)

c     nv = 2**ndim
c     call out_cell(cell,nv,nel)
c     call exitt(1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_geofile (dx,ndim,nv,nel,pmap,io)
      integer pmap(nel),p0,e
      real dx(0:ndim,nv,nel)

      do e=1,nel
         p0 = pmap(e)-1
         x=0
         y=0
         z=0
         do i=1,nv
            x=x+dx(1,i,e)
            y=y+dx(2,i,e)
            z=z+dx(3,i,e)
         enddo
         x=x/nv
         y=y/nv
         z=z/nv
         write(io,1) p0,x,y,z
      enddo
    1 format(i9,1p3e12.4)

      return
      end
c-----------------------------------------------------------------------
      subroutine self_chk(cell,nv,nel,flag)     ! check for not self-ptg.
      integer cell(nv,nel),flag
      integer e
      
      do e=1,nel
      do i=1,nv
         do j=i+1,nv
            if (cell(i,e).eq.cell(j,e)) then
               write(6,*)
               call outmati(cell(1,e),2,4,'SELF!!',e,flag)
               write(6,*)
               write(6,*) 'ABORT SELF-CHK:',i,j,e,flag
               call exitt(flag)
            endif
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine out_geofile2(dx,ndim,nv,nel,cell,nrnk)
      real dx(0:ndim,nv,nel)
      integer cell(nv,nel)

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)
      common /carrayw/ w1   (lpts) , w2   (lpts)
     $               , w3   (lpts) , w4   (lpts)
     $               , w5   (lpts)

      integer e,v,emax,vmax

      call rzero(w4,nrnk)
      do e=1,nel                     ! Get global vertex multiplicities
      do v=1,nv
         i = cell(v,e)
         if (i.gt.nrnk) then
            write(6,1) e,v,i,nrnk,(dx(k,v,e),k=1,3)
    1       format(i9,i3,2i9,1p3e12.4,' i>nrnk! ERROR!')
         else
            w4(i)=w4(i)+1.
         endif
      enddo
      enddo

      do i=1,nrnk
         if (w4(i).gt.0) then
            w4(i)=1./w4(i)
         else
            write(6,*) i,' detected blank index in geofile2'
         endif
      enddo

      call rzero(w1,nrnk)
      call rzero(w2,nrnk)
      call rzero(w3,nrnk)
      do e=1,nel
      do v=1,nv
         i = cell(v,e)
         if (i.le.nrnk) then           ! average global vertex coords
            w1(i)=w1(i)+dx(1,v,e)*w4(i)
            w2(i)=w2(i)+dx(2,v,e)*w4(i)
            w3(i)=w3(i)+dx(3,v,e)*w4(i)
         endif
      enddo
      enddo

      dmax = 0.
      do e=1,nel
      do v=1,nv
         i = cell(v,e)
         if (i.le.nrnk) then    ! check for global/local Euclidian variance
            ddx = dx(1,v,e)-w1(i)
            ddy = dx(2,v,e)-w2(i)
            ddz = dx(3,v,e)-w3(i)
            dd2 = ddx*ddx + ddy*ddy
            if (ndim.eq.3) dd2 = dd2 + ddz*ddz
            if (dd2.ge.dmax) then
               dmax = dd2
               imax = i
               emax = e
               vmax = v
            endif
         endif
      enddo
      enddo
      if (dmax.gt.0) dmax = sqrt(dmax)
      write(6,3) imax,emax,vmax,dmax
      write(6,4) w1(imax),w2(imax),w3(imax),' global xyz'
      write(6,4) (dx(k,vmax,emax),k=1,3)   ,' local  xyz'
    3 format(2i9,i3,1pe12.4,'dmax xyz')
    4 format(1p3e14.6,1x,a11)

c
c     Write coordinates to a file
c
      write(6 ,*) 'Dumping vertex coordinates to unit 12'
      write(12,*) nrnk

      if (ndim.eq.3) then
       do i=1,nrnk
         if (mod(i,10000).eq.0) write(6,6) i,w1(i),w2(i),w3(i)
         write(12,5) w1(i),w2(i),w3(i)
       enddo
      else
       do i=1,nrnk
         if (mod(i,10000).eq.0) write(6,6) i,w1(i),w2(i)
         write(12,5) w1(i),w2(i)
       enddo
      endif

    5 format(1p3e14.6)
    6 format(i9,1x,1p3e14.6)

      return
      end
c-----------------------------------------------------------------------
      subroutine checker(in_elist,in_pmap,nel,ndim,ii)
      integer in_elist(1),in_pmap(1)

      parameter(lelm=1 000 000)  ! DO GLOBAL REPLACE FOR THIS EVERYWHERE !
      parameter(lpts=8*lelm)

      common /qarrayi/ pmap (lpts) , elist(lpts) , w1(lpts)
      integer pmap,elist,w1

      common /arrayr/  dx(4*lpts)

      call outmati(in_pmap ,2,9,'inpmap',nel,1)
      call outmati(in_elist,2,9,'ielist',nel,1)

      call icopy(elist,in_elist,nel)
      call icopy(pmap ,in_pmap ,nel)

      call isort     (elist,w1,nel)
      call iswap_ip  (pmap ,w1,nel)

      call outmati(pmap ,2,9,'s pmap',nel,1)

      nv = 2**ndim
      io = 50+ii
      call out_geofile (dx,ndim,nv,nel,pmap,io)

      return
      end
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
      subroutine rd_bc_bin(cbc,bc,nelv,nelt,ifbswap)

c     .Read Boundary Conditions (and connectivity data)

      parameter (lelm=1 000 000)

      character*3 cbc(6,lelm)
      real        bc (5,6,lelm)
      logical     ifbswap
      
      integer e,f,buf(30)

      npass = 1
      if (nelt.gt.nelv) npass = 2   ! default to thermal topology (for now)

      do kpass = 1,npass

         do e=1,lelm   ! fill up cbc w/ default
         do k=1,6
            cbc(k,e) = 'E  '
         enddo
         enddo

         nwds = 2 + 1 + 5   ! eg + iside + cbc + bc(5,:,:)

         call byte_read(nbc_max,1)
         if (ifbswap) call byte_reverse(nbc_max,1) ! last is char
         do k=1,nbc_max
c           write(6,*) k,' dobc1 ',nbc_max
            call byte_read(buf,nwds)
            if (ifbswap) call byte_reverse(buf,nwds-1) ! last is char
            call buf_to_bc(cbc,bc,buf)
         enddo

      enddo

      return
      end


c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbl,bl,buf)    ! version 1 of binary reader

      parameter(lelm=1 000 000)

      character*3 cbl(6,lelm)
      real        bl(5,6,lelm)
      integer     buf(30)

      integer e,eg,f

      e  = buf(1)
      f  = buf(2)

      call copy4r ( bl(1,f,e),buf(3),5)
      call chcopy (cbl(  f,e),buf(8),3)

c      write(6,1) e,f,cbl(f,e),(bl(k,f,e),k=1,5),' CBC'
c  1   format(i8,i4,2x,a3,5f8.3,1x,a4)

      return
      end

c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,xc,yc,zc,e,ifbswap,ndim)  ! version 1 of binary

      logical ifbswap
      integer e,eg,buf(0:30)

      real xc(8),yc(8),zc(8)

      nwds = 1 + ndim*(2**ndim) ! group + 2x4 for 2d, 3x8 for 3d

      if (ifbswap) call byte_reverse(buf,nwds)

      igroup = buf(0)

      if (ndim.eq.3) then
         call copy4r(xc,buf( 1),8)
         call copy4r(yc,buf( 9),8)
         call copy4r(zc,buf(17),8)
      else
         call copy4r(xc,buf( 1),4)
         call copy4r(yc,buf( 5),4)
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine open_bin_file(ifbswap) ! open file & chk for byteswap

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

      call izero  (fnami,33)

      len = ltrunc(session,80)
      call chcopy (nam1,session,80)
      len = len + 1
      call chcopy (nam1(len),re2,4)
      len = len + 3
      call chcopy (fname,nam1,len)

      call byte_open(fname)
      call byte_read(hdr,20)

      read (hdr,1) version,nelgt,ndum,nelgv
    1 format(a5,i9,i3,i9)

      call byte_read(test,1)
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
c
      test2 = bytetest
      call byte_reverse(test2,1)
      write(6,*) 'Byte swap:',if_byte_swap_test,bytetest,test2

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

      INTEGER FUNCTION INDX132(S1,S2,L2)
      CHARACTER*132 S1,S2
C
      N1=80-L2+1
      INDX132=0
      IF (N1.LT.1) return
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX132=I
            return
         ENDIF
  100 CONTINUE
C
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
    1    format(' ERROR: nic too small in cell2v:',6i10)
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

      nic = i1-i0+1
      call izero(wk,nic)
      do j=1,nv*ncell
         i = cell(j,1)
         wk(i) = wk(i) + 1
      enddo

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
      subroutine greedy (elist,pmap,cell,nv,nelv,nelt,ndim
     $                                              ,ip,ep,ic,jc,ne)
c
c     Distribute elements e > nelv to processors p \in [1,pmax]
c
c     Two pass strategy:
c
c     1st pass:  fill underloaded processors
c
c     2nd pass:  assign any T element to any connected processor p
c                for which ne(p) < nep_max
c
c
c     Required resources:
c
c     ne(p), p(e), pmax, nep_max
c

      integer elist(1),pmap(1)
     $      , cell (nv,1)
     $      , ip(0:1),ep(1)          ! These are 
     $      , ic(1),jc(1)            ! work arrays 
     $      , ne(0:nelv)             ! upon input 

      integer e,p,pnext,pmax


      pmax = iglmax(pmap,nelv)  ! pmap is sorted by e
      call izero(ne,pmax+1)
      do e=1,nelv
         ne(pmap(e)) = ne(pmap(e)) + 1
      enddo
      nepf_min = iglmin(ne(1),pmax)   ! min number of elem/proc, fluid
      nepf_max = iglmax(ne(1),pmax)   ! max number of elem/proc, fluid

      if (nepf_max-nepf_min.gt.1) then
         write(6,*) 'ERROR: nepfmax/min:',nepf_max,nepf_min,nelv,nelt
         call exitt(11)
      endif

      nep_max = nelt/pmax
      if (nep_max*pmax .lt. nelt) nep_max = nep_max + 1

c
c     Build processor-to-element map
c
      call build_proc_el(ip,ep,ne,pmap,pmax,nelv)

c
c     First loop, use greedy to fill underloaded processors
c
      call izero(ne,pmax+1)
      do p=1,pmax
       if (ne(p).lt.nepf_max) then  !  [ip,ep: csr proc-el list]
          j0=ip(p)
          j1=ip(p+1)-1
          do j=j0,j1   ! loop over all elements e on proc p
             e = ep(j)
             do i=1,nv
                iv=cell(i,e)   ! Find elements ke attached to e
                k0 = ic(iv)
                k1 = ic(iv+1)-1
                do k=k0,k1
                   kv = jc(k)  !   This is the structured vertex pointer
                   ke = 1 + (kv-1)/nv      !  = neighboring element
                   if (pmap(ke).le.0) then ! Found an unattached element
                      pmap(ke) = p
                      ne  (p ) = ne(p)+1
                      goto 10
                   endif
                enddo
             enddo
          enddo
       endif
   10  continue  ! next processor
      enddo
c
c     Repeat first loop, using _anything_ to fill underloaded processors
c
      ke_min = nelv
      do p=1,pmax
       if (ne(p).lt.nepf_max) then
          do ke=ke_min+1,nelt
             if (pmap(ke).le.0) then
                pmap(ke) = p
                ne  (p ) = ne(p)+1
                ke_min = ke
                goto 20
             endif
          enddo
       endif
   20  continue
      enddo

c
c     Now, use greedy to distribute remaining elements
c
c     (First version of this code, we revert ot round robin...)
c

      pnext = 0
      do e=nelv+1,nelt
         if (pmap(e).le.0) then
            pnext = pnext+1
            if (pnext.gt.pmax) pnext = 1
            pmap(e) = pnext
            ne  (pnext) = ne(pnext) + 1
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine build_proc_el(ip,ep,ne,pmap,pmax,nel)
      integer ip(0:1),ep(1),ne(0:1),pmap(1),pmax
      integer e,p

c     Step 1:  count number of elems per proc
      call izero(ne,pmax+1) 
      do e=1,nel
         p=pmap(e)   ! pmap is 1-based,  a "0" --> not specified
         ne(p) = ne(p) + 1
      enddo

c     Step 2:  Generate csr pointer

      ip(0) = 1
      do p=0,pmax
         ip(p+1) = ip(p) + ne(p)
         ne(p)   = 0
      enddo

c     Step 3:  fill ep()

      do e=1,nel
         p = pmap(e)
         ep(ip(p) + ne(p)) = e
         ne(p)             = ne(p) + 1 ! This is number filled for proc p
      enddo

      return
      end
c-----------------------------------------------------------------------
      function indxi(list,item,n)
      integer list(1)

      indxi = 0
      do j=1,n
         if (list(j).eq.item) then
            indxi = j
            return
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c2c(ic2c,jc2c,vc2c,ni,nj,cell,nv,nel,iv2c,jv2c,nn,wk)

c     build the c-to-c connectivity matrix in CSR format

      integer ic2c(1),jc2c(1),vc2c(1),cell(nv,nel)
      integer iv2c(1),jv2c(1),wk(1)
      integer e,type,row_sum

c     Step 1:  build the v-to-c connectivity matrix in CSR format

      nic = nn
      njc = nn
      type = 1 ! return vertex-to-cell pointer in (iv2c,jv2c) csr pair
      call cell2v(i0,i1,iv2c,nic,jv2c,njc,cell,nv,nel,type,wk)
      nic = i1-i0+1
c     call outaij(iv2c,jv2c,nic,'iv2c  ',type)


c     Step 2: strip mine

      n_in_strip = 20
      nstrips    = nel / n_in_strip + 1

      k0 = 1                        ! offsets for work array
      k1 = k0 + n_in_strip + 1
      k2 = k1 + n_in_strip + 1

      ie0 = 1
      ic2c(ie0) = 1

      do istrip=1,nstrips
         ie1   = ie0 + n_in_strip - 1
         ie1   = min(ie1,nel)
         ncell = ie1 - ie0 + 1
         call c2cs(ic2c,jc2c,vc2c
     $            ,iv2c,jv2c,i0,i1
     $            ,cell,nv,ie0,ie1
     $            ,wk(k0),wk(k1),wk(k2))
         ie0 = ie1+1
      enddo

c
c     Convert vc2c to graph Laplacian
c
      do e=1,nel
         j0 = ic2c(e)
         j1 = ic2c(e+1)-1
         nj = j1-j0 + 1
         ii = indxi(jc2c(j0),e,nj)
         row_sum = 0
         do k=1,nj
            if (k.ne.ii) then
               j = j0+k-1
               row_sum=row_sum + vc2c(j)
               vc2c(j)=-vc2c(j)           ! off-diag elements < 0
            endif
         enddo
         j = j0+ii-1
         vc2c(j)=row_sum  ! diagonal entry of graph Laplacian = row_sum

         if (ii.eq.0) then
            write(6,*) 'c2c graph, did not find self?',e,j0,j1
            call outmati(jc2c(j0),1,nj,'c2cslf',e,1)
         endif

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine c2cs(ic2c,jc2c,vc2c
     $               ,iv2c,jv2c,k0,k1
     $               ,cell,nv,ie0,ie1
     $               ,nee,ic2t,jc2t)

      integer ic2c(1),jc2c(1),vc2c(1),cell(nv,1)
      integer iv2c(k0:k1),jv2c(1)

      integer nee(ie0:ie1),ic2t(ie0:ie1+1),jc2t(1)


      ne = ie1-ie0 + 1
      call izero(nee(ie0),ne)

      do ie=ie0,ie1              ! Counting phase
      do iv=1,nv
         jv = cell(iv,ie)
         j0 = iv2c(jv)
         j1 = iv2c(jv+1)-1
c        write(6,*) ie,iv,jv,j0,j1,nee(ie),' nee A'
         do jj=j0,j1
            nee(ie) = nee(ie)+1  ! this count too high; compressed below
         enddo
      enddo
      enddo
c     call outmati(nee,1,ne,'nee  1',ne,1)

      ic2t(ie0) = 1
      do ie=ie0,ie1
         ic2t(ie+1) = ic2t(ie) + nee(ie)
      enddo

      j0 = ic2t(ie0)
      j1 = ic2t(ie1)+1
      nj = j1-j0 
      call izero(jc2t(j0),nj)
      call izero(nee(ie0),ne)
  
      do ie=ie0,ie1          ! First filling phase
         i0 = ic2t(ie)
         do iv=1,nv
            jv = cell(iv,ie)
            j0 = iv2c(jv)
            j1 = iv2c(jv+1)-1
            do jj=j0,j1
               je       = jv2c(jj)
               ii       = indxi(jc2t(i0),je,nee(ie))
               if (ii.eq.0) then
                  ji       = i0 + nee(ie)
                  jc2t(ji) = je
                  nee(ie)  = nee(ie)+1     ! # unique elements connected to ie
               endif
            enddo
         enddo
      enddo

      do ie=ie0,ie1
         ic2c(ie+1) = ic2c(ie) + nee(ie)
      enddo
      call izero(nee(ie0),ne)

      do ie=ie0,ie1          ! Compression + 2nd filling phase
         i0 = ic2c(ie)
         do iv=1,nv
            jv = cell(iv,ie)
            j0 = iv2c(jv)
            j1 = iv2c(jv+1)-1
            do jj=j0,j1
               je       = jv2c(jj)
               ii       = indxi(jc2c(i0),je,nee(ie))
               if (ii.eq.0) then
                  ji       = i0 + nee(ie)
                  jc2c(ji) = je
                  vc2c(ji) =  1            ! # connections for pair (ie,je)
                  nee(ie)  = nee(ie)+1     ! # unique elements connected to ie
               else
                  ji       = i0 + ii - 1
                  vc2c(ji) = vc2c(ji) + 1  ! # connections for pair (ie,je)
               endif
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine outbc(cbc,bc,nel,ndim,name6)
      character*3 cbc(6,1)
      real         bc(5,6,1)
      character*6 name6
      integer e,f

      nface = 2*ndim
      do e=1,nel
      do f=1,nface
         write(6,1) e,f,cbc(f,e),(bc(k,f,e),k=1,5),name6
      enddo
      enddo
    1 format(i8,i4,2x,a3,5f8.3,1x,a6)
      return
      end
c-----------------------------------------------------------------------
