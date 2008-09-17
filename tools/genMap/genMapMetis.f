c-----------------------------------------------------------------------
c
c  BEWARE!  This code relies on Metis to do the partitioning.
c
c  Unfortunately, the Metis authors chose to use a nonsymmetric
c  local node numbering for quad and hex elements of the form:
c
c
c      4 ----- 3
c      |       |
c      |       |
c      |       |
c      1 ----- 2
c
c  genmap() uses the symmetric form
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
c-----------------------------------------------------------------------
      program genmap

c     read nekton .rea file and make a .map file


      parameter(lelm=1000000)
      parameter(lpts=8*lelm)

      common /carrayi/ cell (lpts) , pmap (lpts)
     $               , order(lpts) , elist(lpts)
     $               , w1   (lpts) , w2   (lpts)
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


      call makemesh  (cell,nel,irnk,dx,cbc,bc,ndim) ! irnk is # unique points
      nv  = 2**ndim

!     Find all periodic connections, based on cbc info.
      nfc = 2*ndim
      call periodic_vtx
     $               (cell,nv,nel,irnk,dx,ndim,cbc,bc,nfc,w13,w4)


!     Determine number of outflow points and order them last

      call izero (order,irnk)
      call set_outflow(no,order,mo,cell,nv,nel,irnk,cbc,nfc)



!     Recursive bisection of element graph; reverse-order interface points

      call rec_bisect (elist,pmap,order,mo,cell,nv,nel,ndim
     $                                              ,w1,w2,w3,w4,w5)
      call isort     (elist,w1,nel)
      call iswap_ip  (pmap ,w1,nel)

      npts = nv*nel
      call iranku       (cell,nrnk,npts,w1)
      call assign_order (cell,nv,nel,order)


      call iranku       (cell,nrnk,npts,w1) ! make cell numbering contiguous
      call reverse_p    (pmap,nel)          ! lightly load node 0

c     Output to .map file:
      noutflow    = no    ! for now - no outflow bcs
      call out_mapfile (pmap,nel,cell,nv,nrnk,noutflow)


c     call outmati(pmap,13,9,'pmat  ',nel,1)
      open(unit=22,file='p.dat')
      write(22,1) (pmap(k),k=1,nel)
    1 format(i9)
      close(unit=22)

      stop
      end
c-----------------------------------------------------------------------
      subroutine makemesh(cell,nel,irnk,dx,cbc,bc,ndim)

c     read nekton .rea file and make a metis mesh file

      integer     cell(1)
      character*3 cbc (1)
      real         bc (1)
      real         dx (1)

      parameter(lelm=1000000)
      parameter(lpts=8*lelm)

      common /arrayi/ i_n(lpts) , j_n(lpts)
     $              , i_o(lpts) , j_o(lpts)


      call getfile2('Input (.rea) file name:$','.rea$',10)
      call cscan_dxyz  (dx ,nel,ndim)

      nface = 2*ndim
      call cscan_bcs   (cbc,bc,nface,nel,ndim)

      close (unit=10)

c     Compress vertices based on coordinates
      ln = lpts-1
      lo = lpts-1
      call unique_vertex(cell,i_n,j_n,ln,i_o,j_o,lo,dx,nel,ndim)

      nv   = 2**ndim
      npts = nel*nv
      call iranku(cell,irnk,npts,i_n)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt(ie)
      write(6,*) ie,' quit'
      stop
      end
c-----------------------------------------------------------------------
      subroutine cscan_dxyz (dx,nel,ndim)
c
c     Scan for xyz data, read it, and set characteristic length, d2
c
      character*80 string
c
      real dx(1)
      real x(8),y(8),z(8)
      integer e

      integer h2s(8) ! hypercube to strange ordering
      save    h2s
      data    h2s / 1,2,4,3,5,6,8,7 /

      call cscan(string,'MESH DATA',9)
      read (10,*) nel,ndim
      write(6,*) nel,ndim, ' nel,ndim '

      b = 1.e22
      l = 1
      if (ndim.eq.3) then
         do e=1,nel
            read (10,80) string
            read (10,*)   (x(k),k=1,4)
            read (10,*)   (y(k),k=1,4)
            read (10,*)   (z(k),k=1,4)
            read (10,*)   (x(k),k=5,8)
            read (10,*)   (y(k),k=5,8)
            read (10,*)   (z(k),k=5,8)
            do k=1,8
               dx(l+0) = b
               dx(l+1) = x(h2s(k))
               dx(l+2) = y(h2s(k))
               dx(l+3) = z(h2s(k))
               l = l + (ndim+1)
            enddo
         enddo
      else
         do e=1,nel
            read (10,80) string
            read (10,*)   (x(k),k=1,4)
            read (10,*)   (y(k),k=1,4)
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
      call set_d2(dx,nvrt,nel,ndim)
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
      subroutine cscan_bcs  (cbc,bc,nface,nel,ndim)
c
c     Scan for cbc data and read it
c
      character*3  cbc(nface,nel)
      real         bc (5,nface,nel)
      character*80 string

      integer e

      call cscan(string,'BOUNDARY',8) ! for now, fluid only

      read (10,80) string
   80 format(a80)

      ifield = 1
      if (indx1(string,'NO ',3).ne.0) then
         ifield = ifield+1
         call cscan(string,'BOUNDARY',8) ! then, temp only
      endif


      call rd_bc(cbc,bc,nface,nel,ndim,ifield,10)

      return
      end
c-----------------------------------------------------------------------
      subroutine rd_bc(cbc,bc,nfc,nel,ndim,ifield,io)

c     .Read Boundary Conditions (and connectivity data)

      character*3 cbc(nfc,nel)
      real        bc(5,nfc,nel)

      character*3 cbt(nfc)
      real        bt(5,nfc)

      integer eface(6)  ! return Nekton preprocessor face ordering
      save    eface
      data    eface / 4 , 2 , 1 , 3 , 5 , 6 /
c
c
      integer e,f
C
      nbcrea = 5
      do e=1,nel
      do f=1,nfc
         if (nel.lt.1000) then
            read(io,50,err=500,end=500)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
c           write(6,50)
c    $      chtemp,
c    $      cbc(f,e),id1,id2,
c    $      (bc(ii,f,e),ii=1,nbcrea)
   50       format(a1,a3,2i3,5g14.7)
         elseif (nel.lt.100000) then
            read(io,51,err=500,end=500)    
     $      chtemp,
     $      cbc(f,e),id1,id2,
     $      (bc(ii,f,e),ii=1,nbcrea)
   51       format(a1,a3,i5,i1,5g14.7)
         elseif (nel.lt.1000000) then
            read(io,52,err=500,end=500)    
     $      chtemp,
     $      cbc(f,e),id1,
     $      (bc(ii,f,e),ii=1,nbcrea)
   52       format(a1,a3,i6,i1,5g14.7)
         endif
      enddo
      enddo

      do e=1,nel
         call chcopy(cbt,cbc(1,e)  ,nfc*3)
         call copy  ( bt, bc(1,1,e),nfc*5)
         do f=1,nfc
            call copy  ( bc(1,f,e), bt(1,eface(f)),5)
            call chcopy(cbc(  f,e),cbt(  eface(f)),3)
         enddo
         call chcopy(cbt,cbc(1,e),nfc*3)
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
      subroutine unique_vertex(cell,i_n,j_n,ln,i_o,j_o,lo,dx,nel,ndim)
c
      integer i_n(0:ln),j_n(ln),i_o(0:lo),j_o(lo),cell(1)
      real dx(0:ndim,1)

      nvtx = 2**ndim
      do i=1,nvtx*nel
         cell(i) = i
      enddo
      ndx = nvtx*nel

      nv = 2**ndim
c     call out_cell(cell,nv,nel)
c     call exitt(4)


      call build_ohash(i_n,nni,j_n,nnj,ln,i_o,noi,j_o,noj,lo,nh
     $                                                ,dx,cell,ndx,ndim)

      call csr_sort_colj(i_n,j_n,nni,cell)
      call csr_sort_colj(i_o,j_o,noi,cell)

c     call out_csrmati(i_n,j_n,nni,'mat nonov')
c     call out_csrmati(i_o,j_o,noi,'mat    ov')

      nrow = nni

      do i=1,nvtx*nel  ! initialize index set
         cell(i) = -i
      enddo

      do k=1,nrow  ! compare and condense each index set, k=1...,nrow

        jn0 = i_n(k-1)
        jn1 = i_n(k)
        nn  = jn1-jn0

        jo0 = i_o(k-1)
        jo1 = i_o(k)
        no  = jo1-jo0

c       if (mod(k,100).eq.0.or.k.le.10) 
c    $     write(6,*) k,jn0,jn1,jo0,jo1,' condense'

        eps = .10 ! condense anything smaller than .01 x local spacing
        call comp_condense(cell,j_n(jn0),nn,j_o(jo0),no,dx,ndim,eps)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_condense (indx,i_n,nn,i_o,no,dx,ndim,eps)
c
c     compare & condense positions in nonoverlapping (I_n) 
c     and overlapping (I_o) index sets
c
      integer indx(1),i_n(nn),i_o(no),a1,a2
      real    dx(0:ndim,1) ! dx(0,:) = local scale, dx(1,:) = x(:)
c
      eps2 = eps**2

c     write(6,*) nn,no,' nn',eps2

      do ii=1,nn
        i  = i_n(ii)
        i1 = indx(i)
        a1 = abs(i1)
        if (i1.lt.0) then   !  not condensed
          do jj=1,no        !  compare x_i vs. x_j, j in I_o{}
            j  = i_o(jj)
            i2 = indx(j)
            a2 = abs(i2)
            d2 = min(dx(0,i),dx(0,j))
            e2 = eps2*d2
            s2 = 0.
            do k=1,ndim
              s2 = s2 + (dx(k,i)-dx(k,j))**2
            enddo
c           write (6,8) i,j,i1,i2,' COMP',d2,s2,e2
c  8        format( 4i5,a5,1p3e9.2)
            if (s2.lt.e2) then   ! found a match
              if (i2.gt.0) then       ! jj has already been set
                indx(i) = i2
              else
                i0 = min(a1,a2)
                indx(i) = i0
                indx(j) = i0
              endif
            endif
          enddo
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine build_ohash(i_n,nni,j_n,nnj,ln,i_o,noi,j_o,noj,lo,nh
     $                                                ,dx,indx,ndx,ndim)
c
c     Build two hash tables, one based on nonoverlapping domains, 
c     one based on overlapping domains.
c
c     Hash tables returned in csr format with (nh**ndim) rows
c
      integer i_n(0:ln),j_n(ln),i_o(0:lo),j_o(lo),indx(1)
      real dx(0:ndim,1)
      real dmn(0:3),dmx(0:3)

      integer stride
c
      bn =  1.e22
      bx = -1.e22
      call cfill(dmn,bn,4)
      call cfill(dmx,bx,4)

      do i=1,ndx     ! get ranges
         j=indx(i)
         do k=0,ndim
            dmn(k) = min(dmn(k),dx(k,j))
            dmx(k) = max(dmx(k),dx(k,j))
         enddo
      enddo

      zdx  = ndx
      zdx  = zdx**(1./ndim)
      nh   = zdx/3 + 1
c     nh   = 1
      nh   = max(nh,1)
      nrow = nh**ndim

      dmn(0) = sqrt(dmn(0)) ! convert to distance
      dmx(0) = sqrt(dmx(0))

c
c     Set overlapping hash table, csr format
c
      stride = 4
      nho    = stride*nh

c     Build fine hash table, in order to generate overlap
      if (ndim.eq.3) then
         call build_hash3(i_n,nni,j_n,nnj,nho,dmn,dmx,dx,indx,ndx)
      else
         call build_hash2(i_n,nni,j_n,nnj,nho,dmn,dmx,dx,indx,ndx)
      endif
c     call out_csrmati(i_n,j_n,nni,'mat   ov1')


c     concatenate rows of (i_n,j_n) to form (i_o,j_o)

      i_o(0) = 1
      nptr   = 1
      if (ndim.eq.3) then
         do k=1,nh
         do j=1,nh
         do i=1,nh

            iii      = i + nh*( (j-1) + nh*(k-1) )
            i_o(iii) = i_o(iii-1)

            i0 = max(stride*(i-1),1)
            i1 = min(stride*(i+1),stride*nh)
            j0 = max(stride*(j-1),1)
            j1 = min(stride*(j+1),stride*nh)
            k0 = max(stride*(k-1),1)
            k1 = min(stride*(k+1),stride*nh)

            do kk=k0,k1
            do jj=j0,j1
            do ii=i0,i1
               irow     = ii + nho*( (jj-1) + nho*(kk-1) )
               jj0      = i_n(irow-1)
               jj1      = i_n(irow)
               njj      = jj1 - jj0
               if (njj.gt.0) then
                  call icopy(j_o(nptr),j_n(jj0),njj)
                  i_o(iii) = i_o(iii) + njj
                  nptr     = nptr     + njj
               endif
            enddo
            enddo
            enddo
         enddo
         enddo
         enddo
      else              ! 2D
         do j=1,nh
         do i=1,nh

            iii      = i + nh*(j-1)
            i_o(iii) = i_o(iii-1)

            i0 = max(stride*(i-1),1)
            i1 = min(stride*(i+1),stride*nh)
            j0 = max(stride*(j-1),1)
            j1 = min(stride*(j+1),stride*nh)

            do jj=j0,j1
            do ii=i0,i1
               irow     = ii + nho*(jj-1)
               jj0      = i_n(irow-1)
               jj1      = i_n(irow)
               njj      = jj1 - jj0
               if (njj.gt.0) then
                  call icopy(j_o(nptr),j_n(jj0),njj)
                  i_o(iii) = i_o(iii) + njj
                  nptr     = nptr     + njj
               endif
            enddo
            enddo
         enddo
         enddo
      endif

      noi = nrow
      noj = nptr-1

c     Set nonoverlapping hash table, csr format
      if (ndim.eq.3) then
         call build_hash3(i_n,nni,j_n,nnj,nh,dmn,dmx,dx,indx,ndx)
      else
         call build_hash2(i_n,nni,j_n,nnj,nh,dmn,dmx,dx,indx,ndx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine build_hash3(ia,ni,ja,nj,nh,dmn,dmx,dx,indx,ndx)
c
c     build an ( nh x nh x nh ) hash table in csr format
c
      integer ia(0:1),ja(1),indx(1)
      real dx(0:3,1)                 ! scale + x + y + z
      real dmn(0:3),dmx(0:3)

      ddx = (dmx(1)-dmn(1))/nh
      ddy = (dmx(2)-dmn(2))/nh
      ddz = (dmx(3)-dmn(3))/nh

      x0  = dmn(1)
      y0  = dmn(2)
      z0  = dmn(3)

      nrow = nh**3
      call izero(ia(1),nrow)

      do ii=1,ndx     ! loop 1, count number in each row
         i = indx(ii)
         x = dx(1,i)
         y = dx(2,i)
         z = dx(3,i)

         sx = (x - x0)/ddx
         ix = 1 + sx
         ix = min(ix,nh)
         ix = max(ix,1)

         sy = (y - y0)/ddy
         iy = 1 + sy
         iy = min(iy,nh)
         iy = max(iy,1)

         sz = (z - z0)/ddz
         iz = 1 + sz
         iz = min(iz,nh)
         iz = max(iz,1)

         irow     = ix + nh * ( (iy-1) + nh*(iz-1) )
         ia(irow) = ia(irow) + 1
      enddo

      ia(0) = 1      ! convert column count to pointers
      do i=1,nrow
         ia(i) = ia(i-1) + ia(i)
      enddo
      nnzero = ia(nrow)-1
      call izero(ja,nnzero)  ! DANGER, should check size

      do ii=1,ndx     ! loop 2, fill column pointers
         i = indx(ii)
         x = dx(1,i)
         y = dx(2,i)
         z = dx(3,i)

         sx = (x - x0)/ddx
         ix = 1 + sx
         ix = min(ix,nh)
         ix = max(ix,1)

         sy = (y - y0)/ddy
         iy = 1 + sy
         iy = min(iy,nh)
         iy = max(iy,1)

         sz = (z - z0)/ddz
         iz = 1 + sz
         iz = min(iz,nh)
         iz = max(iz,1)

         irow     = ix + nh * ( (iy-1) + nh*(iz-1) )
         j0       = ia(irow-1)
         j1       = ia(irow)
         do j=j0,j1-1
            if (ja(j).eq.0) then
               ja(j) = i            ! put indx in csr column
               goto 100
            endif
         enddo
  100    continue
      enddo

      ni = nrow
      nj = nnzero

      return
      end
c-----------------------------------------------------------------------
      subroutine build_hash2(ia,ni,ja,nj,nh,dmn,dmx,dx,indx,ndx)
c
c     build an ( nh x nh ) hash table in csr format
c
      integer ia(0:1),ja(1),indx(1)
      real dx(0:2,1)                 ! scale + x + y
      real dmn(0:2),dmx(0:2)

      ddx = (dmx(1)-dmn(1))/nh
      ddy = (dmx(2)-dmn(2))/nh

      x0  = dmn(1)
      y0  = dmn(2)

      nrow = nh**2
      call izero(ia(1),nrow)


      do ii=1,ndx     ! loop 1, count number in each row
         i = indx(ii)
         x = dx(1,i)
         y = dx(2,i)

         sx = (x - x0)/ddx
         ix = 1 + sx
         ix = min(ix,nh)
         ix = max(ix,1)

         sy = (y - y0)/ddy
         iy = 1 + sy
         iy = min(iy,nh)
         iy = max(iy,1)

         irow     = ix + nh * (iy-1)
         ia(irow) = ia(irow) + 1
      enddo

      ia(0) = 1      ! convert column count to pointers
      do i=1,nrow
         ia(i) = ia(i-1) + ia(i)
      enddo
      nnzero = ia(nrow)-1
      call izero(ja,nnzero)  ! DANGER, should check size

      do ii=1,ndx     ! loop 2, fill column pointers
         i = indx(ii)
         x = dx(1,i)
         y = dx(2,i)

         sx = (x - x0)/ddx
         ix = 1 + sx
         ix = min(ix,nh)
         ix = max(ix,1)

         sy = (y - y0)/ddy
         iy = 1 + sy
         iy = min(iy,nh)
         iy = max(iy,1)

         irow     = ix + nh * (iy-1)
         j0       = ia(irow-1)
         j1       = ia(irow)
         do j=j0,j1-1
            if (ja(j).eq.0) then
               ja(j) = i            ! put indx in csr column
               goto 100
            endif
         enddo
  100    continue
      enddo

      ni = nrow
      nj = nnzero

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
      function iglmax(a,n)
      integer a(1),tmax
      tmax=-9999999
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

      parameter(lelm=1000000)
      parameter(lpts=8*lelm)
      parameter(mm  =50)

      common /arrayr1/ f(lelm),r(lelm),p(lelm),w(lelm),rr(lelm,mm)
     $               , ev(mm*mm),d(mm),u(mm)

      common /arrayi2/ jdual(lpts) , face (3*lpts)
     $               , idual(lelm) , list (lpts) , list2 (0:lpts)
     $               , elist(lelm) 
      integer jdual,face,idual,list,elist
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
      nfc = 2*ndim
      nvf = 2**(ndim-1)
      call jjnt(elist,nel)
      call build_dualj(idual,jdual,nfc,c,nv,nel,elist,face,nvf,list)
c     call outaij(idual,jdual,nel,'jdual ',9)
      ifconn = is_connected(idual,jdual,nel,list,list2)

      if (ifconn) then
c        nflag = 1
c        np    = 2
c        call METIS_PartMeshDual_f
c    $      (nel,irnk,c,etype,nflag,np,edgecut,pmap,nmap)
         m = mm
         call spec_bis(pmap,idual,jdual,nel,d,u,f,r,p,w,rr,ev,m,ndim)
      else
         n1 = nel/2
         n2 = nel - n1
c        write(6,*) 'not connected',n1,n2,nel
         do i=1,n
            pmap(i)=1
            if (i.gt.n1) pmap(i)=2
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
     $         (n1,n2,pmap(j0),order,mo,elist(j0),n,cell,nv,p,w1,w2,w3)
            write(6,8) i,depth,max_depth,l,j0,j1,n1,n2,p,' DEPTH'
            call outmati(ia,1,l+1,'idepth',depth,1)
            call outmati(da,1,l+1,'ddepth',depth,1)
    8       format(i4,2i3,6i8,a6)

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
      enddo

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

      do k=1,n1         ! Set flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.0) order(i) = -1
         enddo
      enddo

      do k=n1+1,n1+n2         ! Check flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-1) then
                mo = mo+1
                order(i) = mo
            endif
         enddo
      enddo

      do k=1,n1         ! Unset flags
         e = elist(k)
         do j=1,nv
            i = cell(j,e)
            if (order(i).eq.-1) order(i) = 0
         enddo
      enddo

      call out_order(order,mo,elist,cell,nv,n1,n2)

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
      integer w1(1),w2(1)
      integer e,key(2)

      call isort (elist,w1,nel)
      do k=1,nel   ! extract cell sublist
         e = elist(k)
         call icopy(c(1,k),cell(1,e),nv)
      enddo

      npts = nv*nel
      call iranku   (c,irnk,npts,w1)

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
      call ituple_sort(c,2,nel,key,nkey,w1,pmap)
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
      do l=1,nel
         e = elist(l)
c        write(6 ,12) e,(c(k,e),k=1,nv)
      enddo

      write(6 ,10) nel,etype,' SEPRT'
      do l=1,nel
         e = elist(l)
c        write(6,12) e,(o(c(k,e)),k=1,nv)
      enddo

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
c        write(6 ,12) e,(c(k,e),k=1,nv)
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
      return
c
c     Print out copies of a global matrix
c
      if (m.gt.1) then
         write(6,1) nid,m,n,name6
   1     format(3i6,'  Matrix:',2x,a6)
         do i=1,m
            write(6,2) i,name6,(u(i,j),j=1,n)
         enddo
   2     format(i3,1x,a6,20(20i5,/,10x))
      else
         write(6,3) nid,n,name6,(u(1,j),j=1,n)
   3     format(2i3,1x,a6,20(20i5,/,10x))
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
      return
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
      subroutine get_faces(face,nvf,nfc,cell,nv,nel,elist)
c
c     Count number of elements sharing each vertex
c
      integer face(nvf+2,nfc,nel),cell(nv,nel),elist(nel)
      integer e,f

      integer vface(4,6)  ! symm. vertices ordered on symm. faces
      save    vface
      data    vface / 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8
     $              , 1,2,3,4 , 5,6,7,8 /


      do l=1,nel
         e = elist(l)
         do f=1,nfc
            do i=1,nvf
               j = vface(i,f)
               face(i,f,l) = cell(j,e)
            enddo
            face(nvf+1,f,l) = l ! return local element number
            face(nvf+2,f,l) = f
         enddo
      enddo
c     call exitt(7)

      return
      end
c-----------------------------------------------------------------------
      subroutine build_dual(dual,nfc,cell,nv,nel,elist,face,nvf,ind)
c
      integer dual(nfc,nel),cell(nv,nel),elist(1)
c
      integer face(nvf+2,nfc,nel),ind(nfc,nel)
      integer e,f,key(4),w4(4),rank,last(8),next(8)
      logical iftuple_iaeqb,is_same

c     call out_cell(cell,nv,nel)
      call get_faces(face,nvf,nfc,cell,nv,nel,elist)

      call izero(dual,nfc*nel)
      do e=1,nel
      do f=1,nfc
         call isort (face(1,f,e),ind,nvf) ! sort local vertices
      enddo
      enddo

      call jjnt(key,4)
      nkey = nvf
      nfcs = nfc*nel
      nvf2 = nvf+2
      call ituple_sort(face,nvf2,nfcs,key,nkey,ind,w4)


      nvf2 = nvf + 2               !     find matched pairs
      rank = 1
      call icopy(last,face(1,1,1),nvf2)

      i = 0
      do e=1,nel
      do f=1,nfc
         i = i+1
         if (i.gt.1) then
            is_same = iftuple_iaeqb(face(1,f,e),last,key,nkey)
            if (is_same) then
               rank = rank+1
               je   = face(nvf+1,f,e)
               jf   = face(nvf+2,f,e)
               ke   = last(nvf+1)
               kf   = last(nvf+2)
               if (ie.eq.je) then
                  write(6,*) 'dual: found self:',je,jf,kf,rank
               endif
               if (rank.gt.2) then
                  write(6,*) 'dual: high rank:',je,jf,kf,rank
               endif
               dual (jf,je) = ke
               dual (kf,ke) = je
            else
               rank = 1
            endif
            call icopy(last,face(1,i,1),nvf2)
         endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine build_dualj(idual,jdual,nfc,cell,nv,nel,elist
     $                                                   ,face,nvf,dual)

      integer idual(0:nel),jdual(nfc*nel),cell(nv,nel),dual(nfc,nel)
      integer elist(1),face(nvf+2,nfc,nel)
      integer e

      call build_dual(dual,nfc,cell,nv,nel,elist,face,nvf,jdual)

c     Convert jdual to csr format

      l = 0
      idual(0) = 1
      do e=1,nel
         do j=1,nfc
            if (dual(j,e).ne.0) then
               l = l+1
               jdual(l) = dual(j,e)
            endif
         enddo
         idual(e) = l+1
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine set_outflow(no,order,mo,cell,nv,nel,nrnk,cbc,nfc)
c
c     Order outflow nodes last
c
      integer cell(nv,nel),order(1)
      character*3      cbc(nfc,nel)

      parameter(lelm=1000000)
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

      mo = 0
      no = 0
      write(6,*) no,mo,out_vtm,ifoutflow,' OUTFLOW'
      if (.not. ifoutflow) return

      npts             = nel*nv
      call iranku      (cell,nrnk,npts,ind)

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

      character*3      cbc(  nfc,nel)
      real             bc (5,nfc,nel)

      character*3      cb
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

      nmn = irnk
      nmx = 0
      do e=1,nel
      do f=1,nfc
         cb = cbc(f,e)
c        write(6,*) cb,e,f,' cb'
         if (cb.eq.'P  ') then
            je = bc(1,f,e)
            jf = bc(2,f,e)
            jf = efaci(jf)
            call find_connctd_pairs
     $                 (ipair,nvf,e,f,je,jf,cell,nv,dx,ndim)
            
            do k=1,nvf   ! for each cnctd pr, store cross pointer
            do ip=1,2
               jp = 3-ip

               i = ipair(ip,k)
               j = ipair(jp,k)

               nmn = min(i,nmn) ! find active range
               nmn = min(j,nmn)
               nmx = max(i,nmx)
               nmx = max(j,nmx)

               do l=1,ndim  ! at most ndim periodic connections
                  if (iper(l,i).eq.0) then
                     iper(l,i) = j
                     goto 30
                  endif
               enddo
   30          continue

            enddo
            enddo
         endif
      enddo
      enddo

      call jjnt(jmin,irnk)

      do ipass=1,50 ! should be more than enough
         nchange = 0
         do k=nmn,nmx
            do l=1,ndim
               j = iper(l,k)
               if (j.ne.0) then
                  if (jmin(j).lt.jmin(k)) then
                     nchange = nchange+1
                     jmin(k) = jmin(j)
                  endif
               else
                  goto 100
               endif
            enddo
  100       continue
         enddo
         write(6,101) ipass,nchange,nmn,nmx
  101    format(4i9,' ip,nc,mn,mx,periodic')
         if (nchange.eq.0) goto 200
      enddo
  200 continue

c
c     Okay -- we now have the updated pointers, time to update
c     the cell pointers
c

      do e=1,nel
      do v=1,nv
         j = cell(v,e)
         j = min(j,jmin(j))
         cell(v,e) = j
      enddo
      enddo

      npts = nv*nel
      call iranku   (cell,irnk,npts,iper) ! compress cell list

      return
      end
c-----------------------------------------------------------------------
      subroutine find_connctd_pairs(ipair,nvf,e,f,je,jf,cell,nv,dx,ndim)


      integer ipair(2,nvf)
      integer cell(nv,1)
      real dx(0:ndim,nv,1)

      real x0(0:3,4),x1(0:3,4)

      integer e,f
      integer tpair(2,4)

      integer vface(4,6)  ! symm. vertices ordered on symm. faces
      save    vface
      data    vface / 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8
     $              , 1,2,3,4 , 5,6,7,8 /


      nvf = nv/2     ! # vertices/face = 1/2 # vertices/cell

      do i=1,nvf     ! Grab geometry for P-P face pair

c        call outmat(dx(0,j0, e),1,4,' dx0  ', e)
c        call outmat(dx(0,j1,je),1,4,' dx1  ',je)

         j0 = vface (i ,f)
         i0 = cell  (j0,e)
         call copy  (x0(0,i),dx(0,j0, e),ndim+1)
         tpair(1,i) = i0

         j1 = vface (i ,jf)
         i1 = cell  (j1,je)
         call copy  (x1(0,i),dx(0,j1,je),ndim+1)
         tpair(2,i) = i1

      enddo

      do k=1,ndim           ! Subtract off mean of faces
         x0a = 0.
         x1a = 0.
         do i=1,nvf
            x0a = x0a + x0(k,i)/nvf
            x1a = x1a + x1(k,i)/nvf
         enddo
         do i=1,nvf
            x0(k,i) = x0(k,i) - x0a
            x1(k,i) = x1(k,i) - x1a
         enddo
      enddo

      do i=1,4

         d2min = 1.e22
         do j=1,4
            d2 = 0.
            do k=1,ndim
               d2 = d2 + (x0(k,i)-x1(k,j))**2
            enddo
            if (d2.lt.d2min) then
               jmin  = j
               d2min = d2
            endif
         enddo

         ipair(1,i) = tpair(1,i)
         ipair(2,i) = tpair(2,jmin)
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
c     call out_a(a)

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

      do i=1,n
         j0 = abs(ia(i-1))
         j1 = abs(ia(i))
         m  = j1-j0
         call outmati(ja(j0),1,m,name6,i,1)
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
      logical function is_connected(ia,ja,n,jactive,jstack)

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

         iii = 0
         do j=jj,j1-1
            inext = ja(j)
            if (1.le.inext.and.inext.le.n) then  ! range check
               if (jactive(inext).eq.0) then
                  nstack = ipush(jstack,icurr)   ! push the return point
                  jactive(icurr) = j+1
                  jactive(inext) = ia(inext-1)
                  iii = 1
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

      is_connected = .true.
      do i=1,n
         if (jactive(i).eq.0) is_connected = .false.
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
      subroutine spec_bis(pmap,ia,ja,n,d,u,f,r,p,w,rr,ev,m,ndim)

c     n = dimension of A
c     m = max # iterations

      real d(m),u(m),f(n),r(n),p(n),w(n),rr(n,m)
      integer pmap(n),ia(0:n),ja(1)

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

      npass = 3
      do k=1,npass
         niter = m
         call glanczos(rr,n,d,u,niter,f,ia,ja,r,p,w)
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
      subroutine glanczos(rr,n,diag,upper,niter,f,ia,ja,r,p,w)
c
c     Lanczos applied to graph Laplacian
c
      real    rr(n,1),diag(1),upper(1),f(1),r(1),p(1),w(1)
      integer ia(1),ja(1)
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
         call ax(w,p,ia,ja,n)
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
      write(6,6) iter,n,rnorm,rtol,rtr
    6 format(i4,i6,' cg:',1p6e12.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ax(y,x,ia,ja,n)

c     This routine computes y = Ax, where A is the graph Laplacian 

      real y(1),x(1)
      integer ia(0:1),ja(1)
c
      do i=1,n
         j0 = ia(i-1)
         j1 = ia(i) - 1
         m  = j1-j0 + 1
         y(i) = m*x(i)
         do j=j0,j1
            y(i) = y(i) - x(ja(j))
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
c=======================================================================
