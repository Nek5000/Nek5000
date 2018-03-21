      subroutine dProcmapInit()

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      integer disp_unit
      integer*8 winsize, winptr

#ifdef MPI
c      call MPI_Type_Extent(MPI_INTEGER,disp_unit,ierr)
      disp_unit = ISIZE

      winsize = 2*lelt*disp_unit
      call MPI_Win_allocate(winsize,disp_unit,MPI_INFO_NULL,
     $                      nekcomm,winptr,dProcmapH,ierr)

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

#ifdef MPI
      call dProcMapFind(iloc,nids,ieg)
      disp = 2*(iloc-1) + ioff-1

      call mpi_win_lock(MPI_LOCK_EXCLUSIVE,nids,0,dProcmapH,ierr)
      call mpi_put(ibuf,lbuf,MPI_INTEGER,nids,disp,1,MPI_INTEGER,
     $             dProcmapH,ierr)
      call mpi_win_unlock(nids,dProcmapH,ierr)
#else
      call icopy(dProcmapWin(2*(ieg-1) + ioff),ibuf,lbuf)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapGet(ibuf,ieg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer ibuf(2)

      integer*8 disp

      parameter (lc = 128) ! cache size for remote entries
      parameter (lcache = lelt+lc+8-mod(lelt+lc,8)) ! multiple of 8
      integer   cache(lcache,3)
      save      cache

      save icalld
      data icalld /0/

      save iran
      parameter(im = 6075, ia = 106, ic = 1283)

      if (icalld .eq. 0) then
         do i = 1,lelt+lc
            cache(i,1) = -1
         enddo
         icalld = 1
      endif

      ii = lsearch_ur(cache,lcache,ieg)
      if (ii.gt.0 .and. ii.ne.lelt+lc) then ! cache hit
c         write(6,*) nid, 'cache hit ', 'ieg:', ieg
         ibuf(1) = cache(ii,2)
         ibuf(2) = cache(ii,3)
      else
#ifdef MPI
         call dProcmapFind(il,nidt,ieg)
         disp = 2*(il-1)
         call mpi_win_lock(MPI_LOCK_SHARED,nidt,0,dProcmapH,ierr)
         call mpi_get(ibuf,2,MPI_INTEGER,nidt,disp,2,MPI_INTEGER,
     $                dProcmapH,ierr)
         call mpi_win_unlock(nidt,dProcmapH,ierr)
#else
         ibuf(1) = dProcmapWin(2*(ieg-1) + 1)
         ibuf(2) = dProcmapWin(2*(ieg-1) + 2)
#endif
         if (dProcmapCache) then
            ii = ibuf(1)
            if (ibuf(2).ne.nid) then
               iran = mod(iran*ia+ic,im)
               ii = lelt + (lc*iran)/im + 1 ! randomize array location 
            endif
            cache(ii,1) = ieg
            cache(ii,2) = ibuf(1)
            cache(ii,3) = ibuf(2)
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
      integer function lsearch_ur(a, n, k)

      integer a(n), n, k

      parameter(lvec=8) ! unroll factor

      slsearch = 0
      do i = 1,n,lvec
         do j = 0,lvec-1
            if (a(i+j).eq.k) slsearch = i + j
         enddo
         if (slsearch.gt.0) goto 10
      enddo

10    continue
      end
