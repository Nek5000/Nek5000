      module vrthov_mod
      include 'SIZE'            ! = use size_mod: lelt for the /mfi_hs/ bound
      implicit none

c     restart buffers (dedicated, problem-size independent):
c       cb_rdbuf (rdbuf) : file-order read buffer for one batch
c       cb_cmbuf (cmbuf) : CR/RMA redistribution buffer (also the RMA window)
c     both sized to (at least) lrdbuf_bytes/4 words (16 MB) in vrthov_reserve()
      real*4,  allocatable, target :: cb_rdbuf(:)
      real*4,  allocatable, target :: cb_cmbuf(:)

c     per-element floor (worst single element, real*4 words): FP64 + mapab
c     source-order margin (nxr<=lx1+6). one element must fit either buffer.
      integer, parameter :: lrst_elem = 2*(lx1+6)*(ly1+6)*(lz1+6)

c     fixed byte budget per buffer (16 MB); word count = lrdbuf_bytes/4
      integer, parameter :: lrdbuf_bytes = 16*1024*1024 ! 16 MB

c     lrst_idst: RMA idstage (round_rma) length = max read-batch; mfi caps
c            nrst_rd <= lrst_idst ('d' check). Sibling of lrst_hs (both = lelt).
      integer, parameter :: lrst_idst = lelt

c     redist handshake work arrays, bound by lrst_hs=lelt:
c       send-side: nb<=lelt elems/batch;
c       recv-side (it): #rows n1 <= recvcnt <= nelt_hr0 <= lelt (indep. of np)
c       actual usage << lelt (depends on #dests/#sources per rank)
      integer, parameter :: lrst_hs = lelt
      integer kv,ord,ioff,dstlist,cnt,boff,it,ndest
      common /mfi_hs/ kv(2,lrst_hs),ord(lrst_hs),ioff(lrst_hs+1),
     $               dstlist(lrst_hs),cnt(lrst_hs),boff(lrst_hs),
     $               it(3,lrst_hs),ndest

      contains

c     vrthov_reserve: ensure buffers hold >= max(lrdbuf_bytes/4, need) words
c     grown only if `need` (one elem) exceeds the budget (rare, high order).
c     MUST NOT run while cb_cmbuf backs a live window.
      subroutine vrthov_reserve(need)
         implicit none
         integer, optional :: need
         integer ns, ierr

         ns = lrdbuf_bytes / 4        ! default = 16MB budget, 4 bytes per word
         if (present(need)) then
            if (need.gt.ns) ns = need      ! increase only if needed
         endif
         if (mod(ns,2).ne.0) ns = ns + 1   ! keep even (rzero mcm/2)

         if (allocated(cb_rdbuf)) then
            if (size(cb_rdbuf).ge.ns) return ! big enough: no realloc
            deallocate(cb_rdbuf)             ! grow -> both together
         endif
         if (allocated(cb_cmbuf)) deallocate(cb_cmbuf)

         allocate(cb_rdbuf(ns), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_rdbuf$',ns)
         cb_rdbuf = 0
         allocate(cb_cmbuf(ns), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_cmbuf$',ns)
         cb_cmbuf = 0
      end subroutine vrthov_reserve

c     escape hatch: release both buffers. free any MPI window over cb_cmbuf first.
      subroutine vrthov_free
         implicit none
         if (allocated(cb_rdbuf)) deallocate(cb_rdbuf)
         if (allocated(cb_cmbuf))     deallocate(cb_cmbuf)
      end subroutine vrthov_free
      end module vrthov_mod
