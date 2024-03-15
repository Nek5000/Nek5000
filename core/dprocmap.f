#ifdef DPROCMAP
c-----------------------------------------------------------------------
c
c     gllnid and gllel are stored as a distributed array ordered by 
c     the global element index. Access is provided by two
c     functions gllnid() and gllel(). Each ranks holds a local cache
c     for its local and some remote elements.
c
c-----------------------------------------------------------------------
      subroutine dProcmapInit()

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      integer   disp_unit
      integer*8 wsize

#ifdef MPI
      disp_unit = ISIZE
      wsize  = disp_unit*size(dProcmapWin)
      call mpi_comm_dup(nekcomm,commproc,ierr)
      call MPI_Win_create(dProcmapWin,
     $                    wsize,
     $                    disp_unit,
     $                    MPI_INFO_NULL,
     $                    nekcomm,dProcmapH,ierr)

      if (ierr .ne. 0 ) call exitti('MPI_Win_allocate failed!$',0)
#endif
 
      call dProcMapClearCache()
      dProcmapCache = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapPut(ibuf,lbuf,ioff,eg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer ibuf(lbuf),eg
      integer*8 disp

      if (lbuf.lt.1 .or. lbuf.gt.2)
     $   call exitti('invalid lbuf!',lbuf)

#ifdef MPI
      call dProcMapFind(iloc,nids,eg)
      disp = 2*(iloc-1) + ioff

      call mpi_win_lock(MPI_LOCK_EXCLUSIVE,nids,0,dProcmapH,ierr)
      call mpi_put(ibuf,lbuf,MPI_INTEGER,nids,disp,lbuf,MPI_INTEGER,
     $             dProcmapH,ierr)
      call mpi_win_unlock(nids,dProcmapH,ierr)
#else
      call icopy(dProcmapWin(2*(eg-1) + ioff + 1),ibuf,lbuf)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapGet(ibuf,eg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer ibuf(2),eg

      integer*8 disp

      save icalld
      data icalld /0/

      integer ind(3*lcs+lcu)
      common /ctmpx/ ind

      integer n_in_sort,n_in_unsort
      common  /mycache/ n_in_sort,n_in_unsort

      if (icalld .eq. 0) then
         icalld = 1
         call ifill( cache,-1,size( cache))
         call ifill(ucache,-9,size(ucache))
         n_in_sort   = 0
         n_in_unsort = 0
      endif

      ii = ibsearch  (cache(1,3),n_in_sort  ,eg)
      if (ii.gt.0) then  ! cache hit
c        write(6,*) nid,'a cache hit ', 'eg:', eg, n_in_sort, lcs
         ibuf(1) = cache(ii,1)
         ibuf(2) = cache(ii,2)
         return
      endif

      ii = lsearch_ur(ucache(1,3),n_in_unsort,eg)
      if (ii.gt.lcu)call exitti('lsearch_ur returns invalid index$',nid)

      if (ii.gt.0) then ! ucache hit
c        write(6,*) nid,'Aucache hit ', 'eg:', eg, n_in_unsort, lcu
         ibuf(1) = ucache(ii,1)
         ibuf(2) = ucache(ii,2)
         return
      else                                   ! cache miss
c       write(6,*) nid,'A cache miss eg:',eg,n_in_sort,n_in_unsort
#ifdef MPI
         call dProcmapFind(il,nidt,eg)
         disp = 2*(il-1)
         call mpi_win_lock(MPI_LOCK_SHARED,nidt,0,dProcmapH,ierr)
         call mpi_get(ibuf,2,MPI_INTEGER,nidt,disp,2,MPI_INTEGER,
     $                dProcmapH,ierr)
         call mpi_win_unlock(nidt,dProcmapH,ierr)
#else
         call icopy(ibuf,dProcmapWin(2*(eg-1) + 1),2)
#endif
         if (dProcmapCache) then

            n_in_unsort = n_in_unsort + 1      ! Add current ref to unsorted bin
            ucache(n_in_unsort,1) = ibuf(1)
            ucache(n_in_unsort,2) = ibuf(2)
            ucache(n_in_unsort,3) = eg

            if (n_in_unsort.eq.lcu) then ! merge with sorted list

               if (n_in_sort.ge.lcs-lcu) then ! decimate cached data
                  j=0
                  do i=1,n_in_sort,2
                     j=j+1
                     cache(j,1) = cache(i,1)
                     cache(j,2) = cache(i,2)
                     cache(j,3) = cache(i,3)
                  enddo
                  n_in_sort = j
               endif


c              Sort unsorted part of cache

               ipt = n_in_unsort+1
               call isort(ucache(1,3),ind,n_in_unsort)
               call iswap(ucache(1,2),ind,n_in_unsort,ind(ipt))
               call iswap(ucache(1,1),ind,n_in_unsort,ind(ipt))

c              Merge unsorted part of with sorted part

               key = 3
               call merge_tuple
     $             (ind,lcs,mc,3,cache,lcs,n_in_sort
     $                         ,ucache,lcu,n_in_unsort,key)

               call icopy(cache(1,1),ind(1+0*lcs),mc)
               call icopy(cache(1,2),ind(1+1*lcs),mc)
               call icopy(cache(1,3),ind(1+2*lcs),mc)

               n_in_sort   = mc      ! Update counters
               n_in_unsort = 0

            endif

         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcMapFind(il,nids,eg)

      include 'SIZE'
      include 'PARALLEL'

      integer eg

      nstar = nelgt/np
      nids = (eg-1)/nstar
      il = eg - nids * nstar
      if (eg .gt. np*nstar) then
         nids = mod(eg,np) - 1
         il = nstar + 1
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function ibsearch(a,n,k)      !     P. 88-89, NUMERICAL RECIPES
      include 'SIZE'
      include 'PARALLEL'
      integer a(n),n,k

      ibsearch = 0

      if (n.eq.0) return

      ilo=1
      ihi=n

        lcount = 0
    1   if ((ihi-ilo).gt.1) then
           lcount = lcount + 1
           i=(ihi+ilo)/2
           if (a(i).eq.k) goto 10
           if (a(i).gt.k) then
              ihi=i
           else
              ilo=i
           endif
           goto 1
        endif

      return

   10 ibsearch = i

      return
      end
c-----------------------------------------------------------------------
      integer function lsearch_ur(a,n,k)

      integer a(n), n, k
      parameter(lvec = 8) ! unroll factor

      lsearch_ur = 0
      ipt        = 0

      if (n.eq.0) return

      if (nvec.gt.4*lvect) then
         do i = 1,n-lvec,lvec
            do j = 0,lvec-1
               ipt = i+j
               if (a(ipt).eq.k) lsearch_ur = ipt
            enddo
            if (lsearch_ur.gt.0) return
         enddo
      endif


      do j = ipt+1,n
         if (a(j).eq.k) lsearch_ur = j
         if (lsearch_ur.gt.0) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      integer function gllnid(eg)

      include 'SIZE'
      include 'mpif.h'

      integer egl, nidl, icalld
      save    egl, nidl, icalld
      data    egl, nidl, icalld  /0,0,0/

      integer ibuf(2),eg

      include 'DPROCMAP'

      integer n_in_sort,n_in_unsort
      common  /mycache/ n_in_sort,n_in_unsort


      icalld = icalld+1

      if (eg.eq.0) egl = 0

      if (eg.eq.egl) then
         ibuf(2) = nidl
         goto 100
      endif
      call dProcmapGet(ibuf,eg)

c     iout = 60+nid
c     write(iout,*) nid,'gllnid: ',icalld,eg,ibuf(2),ibuf(1)
c     write(6,*)    nid,'gllnid: ',icalld,eg,ibuf(2),ibuf(1)
c    $                            ,n_in_unsort,n_in_sort,lcu,lcs

 100  egl    = eg
      nidl   = ibuf(2)
      gllnid = ibuf(2)


      end
c-----------------------------------------------------------------------
      integer function gllel(eg)

      include 'mpif.h'

      integer egl, iell, eg
      save    egl, iell
      data    egl, iell /0,0/

      integer ibuf(2)

      if (eg.eq.0) egl = 0

      if (eg.eq.egl) then
         ibuf(1) = iell
         goto 100
      endif
      call dProcmapGet(ibuf,eg)

 100  egl   = eg
      iell  = ibuf(1)
      gllel = ibuf(1)

      end
c-----------------------------------------------------------------------
      subroutine dProcMapClearCache 

      include 'SIZE'
      include 'PARALLEL'
      include 'DPROCMAP'

      integer n_in_sort,n_in_unsort
      common  /mycache/ n_in_sort,n_in_unsort

c     write(6,*) nid,' Clearing dpcache:', n_in_sort,n_in_unsort

      call ifill( cache,-1,size(cache) )
      call ifill(ucache,-1,size(ucache))

      itmp = gllnid(0) ! reset last element cache
      itmp = gllel(0)  ! reset last element cache

      n_in_sort   = 0
      n_in_unsort = 0

      end
c-----------------------------------------------------------------------
      subroutine merge_tuple(c,ldc,mc,nc,a,lda,ma,b,ldb,mb,key)

      include 'SIZE'
      include 'PARALLEL'

c     c(mc,nc) = merge(a(ma,nc),b(mb,nc))

      integer c(ldc,nc),a(lda,nc),b(ldb,nc)

      mc = 0  ! Number of merged entries [ c ] = [ a ] + [ b ]

      if (ma.eq.0.and.mb.eq.0) return

      if (mb.eq.0) then
         do k=1,nc
         do i=1,ma
            c(i,k) = a(i,k)
         enddo
         enddo
         mc = ma
         return
      endif

      if (ma.eq.0) then
         do k=1,nc
         do i=1,mb
            c(i,k) = b(i,k)
         enddo
         enddo
         mc = mb
         return
      endif


      i=1
      j=1
      k=1
      do while (i.le.ma.and.j.le.mb)
         if (a(i,key).lt.b(j,key)) then
            do l=1,nc
               c(k,l) = a(i,l)
            enddo
            i=i+1
         else
            do l=1,nc
               c(k,l) = b(j,l)
            enddo
            j=j+1
         endif

         if (k.eq.1) then
            k=k+1
         elseif (c(k,key).gt.c(k-1,key)) then ! Avoid repeats
            k=k+1
         endif
      enddo

      do while (i.le.ma)
         do l=1,nc
            c(k,l) = a(i,l)
         enddo
         i=i+1
         k=k+1
      enddo

      do while (j.le.mb)
         do l=1,nc
            c(k,l) = b(j,l)
         enddo
         j=j+1
         k=k+1
      enddo

      mc = k-1 ! Length of merged list

      return
      end
c-----------------------------------------------------------------------
#endif
