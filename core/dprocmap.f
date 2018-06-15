      subroutine dProcmapInit()

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      integer   disp_unit
      integer*8 winsize

#ifdef MPI
      disp_unit = ISIZE
      winsize = disp_unit*size(dProcmapWin)
      call MPI_Win_create(dProcmapWin,winsize,disp_unit,MPI_INFO_NULL,
     $                    nekcomm,dProcmapH,ierr)

      if (ierr .ne. 0 ) call exitti('MPI_Win_allocate failed!$',0)
#endif
      
      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapPut(ibuf,lbuf,ioff,ieg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer ibuf(lbuf)
      integer*8 disp

      if (lbuf.lt.l .or. lbuf.gt.3)
     $   call exitti('invalid lbuf!',lbuf)

#ifdef MPI
      call dProcMapFind(iloc,nids,ieg)
      disp = 3*(iloc-1) + ioff

      call mpi_win_lock(MPI_LOCK_EXCLUSIVE,nids,0,dProcmapH,ierr)
csk      call mpi_win_lock(MPI_LOCK_SHARED,nids,0,dProcmapH,ierr)
      call mpi_put(ibuf,lbuf,MPI_INTEGER,nids,disp,lbuf,MPI_INTEGER,
     $             dProcmapH,ierr)
      call mpi_win_unlock(nids,dProcmapH,ierr)
#else
      call icopy(dProcmapWin(3*(ieg-1) + ioff + 1),ibuf,lbuf)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapGet(ibuf,ieg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer ibuf(3)

      integer*8 disp

      ! cache for local and remote elements
      parameter (lcr = 128)
      parameter (lc = lelt+lcr+8-mod(lelt+lcr,8)) ! multiple of 8
      integer   cache(lc,3)
      save      cache

      save icalld
      data icalld /0/

      save iran
      parameter(im = 6075, ia = 106, ic = 1283)

      if (icalld .eq. 0) then
         call ifill(cache,-1,size(cache))
         icalld = 1
      endif

      ii = lsearch_ur(cache(1,3),lc,ieg)
      if (ii.gt.lc) call exitti('lsearch_ur returns invalid index$',ii)
      if (ii.gt.0 .and. ii.ne.lelt+lcr) then ! cache hit
c         write(6,*) nid, 'cache hit ', 'ieg:', ieg
         ibuf(1) = cache(ii,1)
         ibuf(2) = cache(ii,2)
      else
#ifdef MPI
         call dProcmapFind(il,nidt,ieg)
         disp = 3*(il-1)
         call mpi_win_lock(MPI_LOCK_SHARED,nidt,0,dProcmapH,ierr)
         call mpi_get(ibuf,3,MPI_INTEGER,nidt,disp,3,MPI_INTEGER,
     $                dProcmapH,ierr)
         call mpi_win_unlock(nidt,dProcmapH,ierr)
#else
         call icopy(ibuf,dProcmapWin(3*(ieg-1) + 1),3)
#endif
         if (dProcmapCache) then
            ii = ibuf(1)
            if (ibuf(2).ne.nid) then
               iran = mod(iran*ia+ic,im)
               ii = lelt + (lcr*iran)/im + 1 ! randomize array location 
            endif
            cache(ii,1) = ibuf(1)
            cache(ii,2) = ibuf(2)
            cache(ii,3) = ieg
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcMapFind(il,nids,ieg)

      include 'SIZE'
      include 'PARALLEL'

      ! distribute array in blocks across ranks
      nstar = nelgt/np
      nids = (ieg-1)/nstar
      il = ieg - nids * nstar
      if (ieg .gt. np*nstar) then
         nids = mod(ieg,np) - 1
         il = nstar + 1
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function lsearch_ur(a,n,k)

      integer a(n), n, k
      parameter(lvec = 8) ! unroll factor

      lsearch_ur = 0
      do i = 1,n,lvec
         do j = 0,lvec-1
            if (a(i+j).eq.k) lsearch_ur = i + j
         enddo
         if (lsearch_ur.gt.0) goto 10
      enddo

10    continue
      end
