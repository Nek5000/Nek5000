c-----------------------------------------------------------------------
      SUBROUTINE MXM(A,N1,B,N2,C,N3)
C-----------------------------------------------------------------------
C
C     NOTE -- this code has been set up with the "mxmf3" routine
c             referenced in memtime.f.   On most machines, the f2
c             and f3 versions give the same performance (f2 is the
c             nekton standard).  On the t3e, f3 is noticeably faster.
c             pff  10/5/98
C
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C----------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
c
      include 'SIZE'
      include 'OPCTR'
      include 'TOTAL'
c
      integer mxmid
      integer aligned
      integer n3m4, n3r4

      call aligned3(a,b,c,aligned)

      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'mxm   '
      endif
      isbcnt = n1*n3*(2*n2-1)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
c

      mxmid = 0;
C      call checktime1(n1,n2,n3,mxmid,a,b,c)
      if (n2.eq.2) then
         call mxm44_2(a,n1,b,n2,c,n3)
      else if ((aligned.eq.1) .and.
     $         (n1.ge.8) .and. (n2.ge.8) .and. (n3.ge.8) .and.
     $         (modulo(n1,2).eq.0) .and. (modulo(n2,2).eq.0) ) then
         if (modulo(n3,4).eq.0) then
            call mxm44(a,n1,b,n2,c,n3)
         else
            call mxm44_uneven(a,n1,b,n2,c,n3)
         endif
      else if((aligned.eq.1) .and.
     $        (modulo(n1,6).eq.0) .and. (modulo(n3,6).eq.0) .and.
     $        (n2.ge.4) .and. (modulo(n2,2).eq.0) ) then
         call mxm3(a,n1,b,n2,c,n3)
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
C      call checktime2
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine amxm(a,n1,b,n2,c,n3)
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
c
      include 'SIZE'
      include 'OPCTR'
      include 'TOTAL'
c
      integer opct
      opct = opct + n1*n3*(2*n2-1)
c
      if (n2.eq.2) then
         call mxm44_2(a,n1,b,n2,c,n3)
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine fmxm(a,n1,b,n2,c,n3)
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
c
      include 'SIZE'
      include 'OPCTR'
      include 'TOTAL'
c
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'fmxm  '
      endif
      isbcnt = n1*n3*(2*n2-1)
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
c
      if (n2.eq.2) then
         call mxm44_2(a,n1,b,n2,c,n3)
      else
         call mxm44_0(a,n1,b,n2,c,n3)
      endif
c
      return
      end
c
c-----------------------------------------------------------------------
      subroutine mxm44_0(a, m, b, k, c, n)
*      subroutine matmul44(m, n, k, a, lda, b, ldb, c, ldc)
*      real*8 a(lda,k), b(ldb,n), c(ldc,n)
      real a(m,k), b(k,n), c(m,n)
      real s11, s12, s13, s14, s21, s22, s23, s24
      real s31, s32, s33, s34, s41, s42, s43, s44
c
c matrix multiply with a 4x4 pencil 
c

      mresid = iand(m,3) 
      nresid = iand(n,3) 
      m1 = m - mresid + 1
      n1 = n - nresid + 1

      do i=1,m-mresid,4
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s41 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s42 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          s43 = 0.0d0
          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0
          s44 = 0.0d0
          do l=1,k
            s11 = s11 + a(i,l)*b(l,j)
            s12 = s12 + a(i,l)*b(l,j+1)
            s13 = s13 + a(i,l)*b(l,j+2)
            s14 = s14 + a(i,l)*b(l,j+3)

            s21 = s21 + a(i+1,l)*b(l,j)
            s22 = s22 + a(i+1,l)*b(l,j+1)
            s23 = s23 + a(i+1,l)*b(l,j+2)
            s24 = s24 + a(i+1,l)*b(l,j+3)

            s31 = s31 + a(i+2,l)*b(l,j)
            s32 = s32 + a(i+2,l)*b(l,j+1)
            s33 = s33 + a(i+2,l)*b(l,j+2)
            s34 = s34 + a(i+2,l)*b(l,j+3)

            s41 = s41 + a(i+3,l)*b(l,j)
            s42 = s42 + a(i+3,l)*b(l,j+1)
            s43 = s43 + a(i+3,l)*b(l,j+2)
            s44 = s44 + a(i+3,l)*b(l,j+3)
          enddo
          c(i,j)     = s11 
          c(i,j+1)   = s12 
          c(i,j+2)   = s13
          c(i,j+3)   = s14

          c(i+1,j)   = s21 
          c(i+2,j)   = s31 
          c(i+3,j)   = s41 

          c(i+1,j+1) = s22
          c(i+2,j+1) = s32
          c(i+3,j+1) = s42

          c(i+1,j+2) = s23
          c(i+2,j+2) = s33
          c(i+3,j+2) = s43

          c(i+1,j+3) = s24
          c(i+2,j+3) = s34
          c(i+3,j+3) = s44
        enddo
* Residual when n is not multiple of 4
        if (nresid .ne. 0) then
          if (nresid .eq. 1) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,n)
              s21 = s21 + a(i+1,l)*b(l,n)
              s31 = s31 + a(i+2,l)*b(l,n)
              s41 = s41 + a(i+3,l)*b(l,n)
            enddo
            c(i,n)     = s11 
            c(i+1,n)   = s21 
            c(i+2,n)   = s31 
            c(i+3,n)   = s41 
          elseif (nresid .eq. 2) then
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
            enddo
            c(i,j)     = s11 
            c(i,j+1)   = s12

            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 

            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
          else
            s11 = 0.0d0
            s21 = 0.0d0
            s31 = 0.0d0
            s41 = 0.0d0
            s12 = 0.0d0
            s22 = 0.0d0
            s32 = 0.0d0
            s42 = 0.0d0
            s13 = 0.0d0
            s23 = 0.0d0
            s33 = 0.0d0
            s43 = 0.0d0
            do l=1,k
              s11 = s11 + a(i,l)*b(l,j)
              s12 = s12 + a(i,l)*b(l,j+1)
              s13 = s13 + a(i,l)*b(l,j+2)

              s21 = s21 + a(i+1,l)*b(l,j)
              s22 = s22 + a(i+1,l)*b(l,j+1)
              s23 = s23 + a(i+1,l)*b(l,j+2)

              s31 = s31 + a(i+2,l)*b(l,j)
              s32 = s32 + a(i+2,l)*b(l,j+1)
              s33 = s33 + a(i+2,l)*b(l,j+2)

              s41 = s41 + a(i+3,l)*b(l,j)
              s42 = s42 + a(i+3,l)*b(l,j+1)
              s43 = s43 + a(i+3,l)*b(l,j+2)
            enddo
            c(i,j)     = s11 
            c(i+1,j)   = s21 
            c(i+2,j)   = s31 
            c(i+3,j)   = s41 
            c(i,j+1)   = s12 
            c(i+1,j+1) = s22
            c(i+2,j+1) = s32
            c(i+3,j+1) = s42
            c(i,j+2)   = s13
            c(i+1,j+2) = s23
            c(i+2,j+2) = s33
            c(i+3,j+2) = s43
          endif
        endif
      enddo

* Residual when m is not multiple of 4
      if (mresid .eq. 0) then
        return
      elseif (mresid .eq. 1) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,j)
            s12 = s12 + a(m,l)*b(l,j+1)
            s13 = s13 + a(m,l)*b(l,j+2)
            s14 = s14 + a(m,l)*b(l,j+3)
          enddo
          c(m,j)     = s11 
          c(m,j+1)   = s12 
          c(m,j+2)   = s13
          c(m,j+3)   = s14
        enddo
* mresid is 1, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n)
          enddo
          c(m,n) = s11
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s12 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-1)
            s12 = s12 + a(m,l)*b(l,n)
          enddo
          c(m,n-1) = s11
          c(m,n) = s12
          return
        else
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          do l=1,k
            s11 = s11 + a(m,l)*b(l,n-2)
            s12 = s12 + a(m,l)*b(l,n-1)
            s13 = s13 + a(m,l)*b(l,n)
          enddo
          c(m,n-2) = s11
          c(m,n-1) = s12
          c(m,n) = s13
          return
        endif          
      elseif (mresid .eq. 2) then
        do j=1,n-nresid,4
          s11 = 0.0d0
          s12 = 0.0d0
          s13 = 0.0d0
          s14 = 0.0d0
          s21 = 0.0d0
          s22 = 0.0d0
          s23 = 0.0d0
          s24 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,j)
            s12 = s12 + a(m-1,l)*b(l,j+1)
            s13 = s13 + a(m-1,l)*b(l,j+2)
            s14 = s14 + a(m-1,l)*b(l,j+3)

            s21 = s21 + a(m,l)*b(l,j)
            s22 = s22 + a(m,l)*b(l,j+1)
            s23 = s23 + a(m,l)*b(l,j+2)
            s24 = s24 + a(m,l)*b(l,j+3)
          enddo
          c(m-1,j)   = s11 
          c(m-1,j+1) = s12 
          c(m-1,j+2) = s13
          c(m-1,j+3) = s14
          c(m,j)     = s21
          c(m,j+1)   = s22 
          c(m,j+2)   = s23
          c(m,j+3)   = s24
        enddo
* mresid is 2, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n)
          enddo
          c(m-1,n) = s11
          c(m,n) = s21
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-1)
            s12 = s12 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-1)
            s22 = s22 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-1) = s11
          c(m-1,n)   = s12
          c(m,n-1)   = s21
          c(m,n)     = s22
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-1,l)*b(l,n-2)
            s12 = s12 + a(m-1,l)*b(l,n-1)
            s13 = s13 + a(m-1,l)*b(l,n)
            s21 = s21 + a(m,l)*b(l,n-2)
            s22 = s22 + a(m,l)*b(l,n-1)
            s23 = s23 + a(m,l)*b(l,n)
          enddo
          c(m-1,n-2) = s11
          c(m-1,n-1) = s12
          c(m-1,n)   = s13
          c(m,n-2)   = s21
          c(m,n-1)   = s22
          c(m,n)     = s23
          return
        endif
      else
* mresid is 3
        do j=1,n-nresid,4
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0

          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0

          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0

          s14 = 0.0d0
          s24 = 0.0d0
          s34 = 0.0d0

          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,j)
            s12 = s12 + a(m-2,l)*b(l,j+1)
            s13 = s13 + a(m-2,l)*b(l,j+2)
            s14 = s14 + a(m-2,l)*b(l,j+3)

            s21 = s21 + a(m-1,l)*b(l,j)
            s22 = s22 + a(m-1,l)*b(l,j+1)
            s23 = s23 + a(m-1,l)*b(l,j+2)
            s24 = s24 + a(m-1,l)*b(l,j+3)

            s31 = s31 + a(m,l)*b(l,j)
            s32 = s32 + a(m,l)*b(l,j+1)
            s33 = s33 + a(m,l)*b(l,j+2)
            s34 = s34 + a(m,l)*b(l,j+3)
          enddo
          c(m-2,j)   = s11 
          c(m-2,j+1) = s12 
          c(m-2,j+2) = s13
          c(m-2,j+3) = s14

          c(m-1,j)   = s21 
          c(m-1,j+1) = s22
          c(m-1,j+2) = s23
          c(m-1,j+3) = s24

          c(m,j)     = s31 
          c(m,j+1)   = s32
          c(m,j+2)   = s33
          c(m,j+3)   = s34
        enddo
* mresid is 3, check nresid
        if (nresid .eq. 0) then
          return
        elseif (nresid .eq. 1) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n)
          enddo
          c(m-2,n) = s11
          c(m-1,n) = s21
          c(m,n) = s31
          return
        elseif (nresid .eq. 2) then
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-1)
            s12 = s12 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-1)
            s22 = s22 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-1)
            s32 = s32 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-1) = s11
          c(m-2,n)   = s12
          c(m-1,n-1) = s21
          c(m-1,n)   = s22
          c(m,n-1)   = s31
          c(m,n)     = s32
          return
        else
          s11 = 0.0d0
          s21 = 0.0d0
          s31 = 0.0d0
          s12 = 0.0d0
          s22 = 0.0d0
          s32 = 0.0d0
          s13 = 0.0d0
          s23 = 0.0d0
          s33 = 0.0d0
          do l=1,k
            s11 = s11 + a(m-2,l)*b(l,n-2)
            s12 = s12 + a(m-2,l)*b(l,n-1)
            s13 = s13 + a(m-2,l)*b(l,n)
            s21 = s21 + a(m-1,l)*b(l,n-2)
            s22 = s22 + a(m-1,l)*b(l,n-1)
            s23 = s23 + a(m-1,l)*b(l,n)
            s31 = s31 + a(m,l)*b(l,n-2)
            s32 = s32 + a(m,l)*b(l,n-1)
            s33 = s33 + a(m,l)*b(l,n)
          enddo
          c(m-2,n-2) = s11
          c(m-2,n-1) = s12
          c(m-2,n)   = s13
          c(m-1,n-2) = s21
          c(m-1,n-1) = s22
          c(m-1,n)   = s23
          c(m,n-2)   = s31
          c(m,n-1)   = s32
          c(m,n)     = s33
          return
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mxm44_2(a, m, b, k, c, n)
      real a(m,2), b(2,n), c(m,n)

      nresid = iand(n,3) 
      n1 = n - nresid + 1

      do j=1,n-nresid,4
         do i=1,m
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
            c(i,j+1) = a(i,1)*b(1,j+1)
     $             + a(i,2)*b(2,j+1)
            c(i,j+2) = a(i,1)*b(1,j+2)
     $             + a(i,2)*b(2,j+2)
            c(i,j+3) = a(i,1)*b(1,j+3)
     $             + a(i,2)*b(2,j+3)
         enddo
      enddo
      if (nresid .eq. 0) then
        return
      elseif (nresid .eq. 1) then
         do i=1,m
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      elseif (nresid .eq. 2) then
         do i=1,m
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      else
         do i=1,m
            c(i,n-2) = a(i,1)*b(1,n-2)
     $             + a(i,2)*b(2,n-2)
            c(i,n-1) = a(i,1)*b(1,n-1)
     $             + a(i,2)*b(2,n-1)
            c(i,n) = a(i,1)*b(1,n)
     $             + a(i,2)*b(2,n)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
