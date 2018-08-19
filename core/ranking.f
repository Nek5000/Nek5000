      subroutine isort8(a, ind, n)
      integer*8 a(1), aa
      integer ind(1)

      do 10 j = 1, n
         ind(j) = j
   10 continue

      if (n .le. 1) return
      L = n / 2 + 1
      ir = n
  100 continue
         if (l .gt. 1) then
            l = l - 1
            aa = a(l)
            ii = ind(l)
         else
            aa = a(ir)
            ii = ind(ir)
            a(ir) = a(1)
            ind(ir) = ind(1)
            ir = ir - 1
            if (ir .eq. 1) then
               a(1) = aa
               ind(1) = ii
               return
            endif
         endif
         i = l
         j = l + l
  200    continue
         if (j .le. ir) then
            if (j .lt. ir) then
               if (a(j) .lt. a(j + 1)) j = j + 1
            endif
            if (aa .lt. a(j)) then
               a(i) = a(j)
               ind(i) = ind(j)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
         goto 200
         endif
         a(i) = aa
         ind(i) = ii
      goto 100

      return
      end
c-----------------------------------------------------------------------
      subroutine iswap8_ip(x, p, n)
      integer*8 x(1), xstart
      integer p(1)

      do k = 1, n
         if (p(k) .gt. 0) then ! Not swapped
            xstart = x(k)
            loop_start = k
            last = k
            do j = k, n
               next = p(last)
               if (next .lt. 0) then
                  write(6, *) 'Hey! iswap_ip problem.', j, k, n, next
                  call exitt
               elseif (next .eq. loop_start) then
                  x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  x(last) = x(next)
                  p(last) = -p(last)
                  last = next
               endif
            enddo
   10       continue
         endif
      enddo

      do k = 1, n
         p(k) = -p(k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine isortcols(a, m, key, col, ind, n)
      integer a(m, 1), col(m)
      integer ind(1)

      do 10 j = 1, n
         ind(j) = j
   10 continue

      if (n .le. 1) return
      L = n / 2 + 1
      ir = n
  100 continue
         if (l .gt. 1) then
            l = l - 1
            call icopy(col, a(1, l), m)
            ii = ind(l)
         else
            call icopy(col, a(1, ir), m)
            call icopy(a(1, ir), a(1, 1), m)
            ii = ind(ir)
            ind(ir) = ind( 1)
            ir = ir - 1
            if (ir .eq. 1) then
               call icopy(a(1, ir), col, m)
               ind(1) = ii
               return
            endif
         endif
         i = l
         j = l + l
  200    continue
         if (j .le. ir) then
            if (j .lt. ir) then
               if (a(key, j) .lt. a(key, j + 1)) j = j + 1
            endif
            if (col(key) .lt. a(key, j)) then
               call icopy(a(1, i), a(1, j), m)
               ind(i) = ind(j)
               i = j
               j = j + j
            else
               j = ir + 1
            endif
         goto 200
         endif
         call icopy(a(1, i), col, m)
         ind(i) = ii
      goto 100
      end
c-----------------------------------------------------------------------
      subroutine i8cadd(a, const, n)
      integer*8 a(n), const

      do i = 1, n
         a(i) = a(i) + const
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine glob_rank_i8(rank, n, ind, wk, ldw)
      include 'SIZE'
      include 'PARALLEL'

      parameter (lt = lx1 * ly1 * lz1 * lelv)
      common /nekmpi/ mid, mp, nekcomm, nekgroup, nekreal

      integer*8 rank(n), icount8, input_max, rlast
      integer ind(ldw), wk(2, ldw), col(2)

      integer*4 cr_hndl
      save cr_hndl
      data cr_hndl /-999/

      real dummy

      ! Check if we can swap int*8 as "real".  (Ok in double precision.)
      if (wdsize .lt. 8) then
         call exitti('ABORT: glob_rank_i8 wdsize<8.$', ldw)
      endif

      if (cr_hndl .eq. -999) then
         call fgslib_crystal_setup(cr_hndl, nekcomm, np) 
      endif

      input_max = i8glmax(rank, n)                    ! Data range for binning


      do i = 1, n
         wk(1, i) = (np * (rank(i) - 1)) / input_max  ! Target processor
         wk(2, i) = i                                 ! Local position
      enddo

      nmax = ldw ! Bound on # of returns. (From declaration of wk)
      call fgslib_crystal_tuple_transfer 
     &        (cr_hndl, n, nmax, wk, 2, rank, 1, dummy, 0, 1)
      nmx = iglmax(n, 1)

      if (nmx .gt. ldw) then
         call exitti('ABORT: glob_rank_i8 nmax>ldw.$', ldw)
      endif

      call isort8(rank, ind, n)   ! Sort by original global number 
      call iswap8_ip(wk, ind, n)  ! Swap 2x32 bit ints as 64 bit ints

      icount = 0                  ! Establish local ranking
      rlast = wk(2, 1) - 1
      do i = 1, n
         if (rlast .ne. rank(i)) icount = icount + 1 ! Bump counter indicating
         if (rlast .ne. rank(i)) rlast = rank(i)     ! number of unique entries
         rank(i) = icount
      enddo

      icount8 = icount            ! Establish global ranking
      icount8 = i8gl_running_sum(icount8) - icount
      call i8cadd(rank, icount8, n)

      call fgslib_crystal_tuple_transfer !  Restore rank-ordered lists
     &        (cr_hndl, n, nmax, wk, 2, rank, 1, dummy, 0, 1)

      nmx = iglmax(n, 1)

      if (nmx .gt. ldw) then
         call exitti('ABORT: glob_rank_i8 nmax>ldw.$', ldw)
      endif

      call isortcols(wk, 2, 2, col, ind, n) ! Sort by original global number 
      call iswap8_ip(rank, ind, n)          ! Swap new rank back to original
                                            ! location
      return
      end
c-----------------------------------------------------------------------
      subroutine matrix_distribution(mat_dist)
      include 'SIZE'
      include 'TOTAL'

      common /nekmpi/ mid, mp, nekcomm, nekgroup, nekreal

      parameter (lt = lx1 * ly1 * lz1 * lelv)
      integer*8 mat_dist(lt), icount8, ngv
      integer gsl, icount

      common /ivrtx/ vertex ((2 ** ldim) * lelt)
      integer vertex

      common /ctmp0/ iwk(2, 2 * lt), ind(2 * lt)

      nxyz = nx1 * ny1 * nz1
      n = nxyz * nelv

      call set_vert(mat_dist, ngv, nx1, nelv, vertex, .true.) ! std glo_num
      call fgslib_gs_setup(gsl, mat_dist, n, nekcomm, np)     ! Need centers > 0.

      icount8 = n                              ! Establish new numbering
      icount8 = i8gl_running_sum(icount8) - n  ! for rank 0 to P-1

      do i = 1, n
         if (pmask(i, 1, 1, 1) > 0.0) then
            mat_dist(i) = i + icount8
         else
            mat_dist(i) = -1
         endif
      enddo

      call fgslib_gs_op(gsl, mat_dist, 3, 3, 0)        ! make shared id's unique
      call glob_rank_i8(mat_dist, n, ind, iwk, 2 * lt) ! Compress global index set

      return
      end
c-----------------------------------------------------------------------

