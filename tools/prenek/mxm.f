      SUBROUTINE MXM(A,N1,B,N2,C,N3)
C----------------------------------------------------------------------
C
C     Matrix-vector product routine. 
C     NOTE: Use assembly coded routine if available.
C
C---------------------------------------------------------------------
      REAL A(N1,N2),B(N2,N3),C(N1,N3)
C
      if (n2.le.8) then
         if (n2.eq.1) then
            call mxm1(a,n1,b,n2,c,n3)
         elseif (n2.eq.2) then
            call mxm2(a,n1,b,n2,c,n3)
         elseif (n2.eq.3) then
            call mxm3(a,n1,b,n2,c,n3)
         elseif (n2.eq.4) then
            call mxm4(a,n1,b,n2,c,n3)
         elseif (n2.eq.5) then
            call mxm5(a,n1,b,n2,c,n3)
         elseif (n2.eq.6) then
            call mxm6(a,n1,b,n2,c,n3)
         elseif (n2.eq.7) then
            call mxm7(a,n1,b,n2,c,n3)
         else
            call mxm8(a,n1,b,n2,c,n3)
         endif
      elseif (n2.le.16) then
         if (n2.eq.9) then
            call mxm9(a,n1,b,n2,c,n3)
         elseif (n2.eq.10) then
            call mxm10(a,n1,b,n2,c,n3)
         elseif (n2.eq.11) then
            call mxm11(a,n1,b,n2,c,n3)
         elseif (n2.eq.12) then
            call mxm12(a,n1,b,n2,c,n3)
         elseif (n2.eq.13) then
            call mxm13(a,n1,b,n2,c,n3)
         elseif (n2.eq.14) then
            call mxm14(a,n1,b,n2,c,n3)
         elseif (n2.eq.15) then
            call mxm15(a,n1,b,n2,c,n3)
         else
            call mxm16(a,n1,b,n2,c,n3)
         endif
      else
         N0=N1*N3
         DO 10 I=1,N0
            C(I,1)=0.
 10      CONTINUE
         DO 100 J=1,N3
         DO 100 K=1,N2
         BB=B(K,J)
         DO 100 I=1,N1
            C(I,J)=C(I,J)+A(I,K)*BB
 100     CONTINUE
      ENDIF
      RETURN
      END
c
      subroutine mxm1(a,n1,b,n2,c,n3)
c
      real a(n1,1),b(1,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
         enddo
      enddo
      return
      end
      subroutine mxm2(a,n1,b,n2,c,n3)
c
      real a(n1,2),b(2,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
         enddo
      enddo
      return
      end
      subroutine mxm3(a,n1,b,n2,c,n3)
c
      real a(n1,3),b(3,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
         enddo
      enddo
      return
      end
      subroutine mxm4(a,n1,b,n2,c,n3)
c
      real a(n1,4),b(4,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
         enddo
      enddo
      return
      end
      subroutine mxm5(a,n1,b,n2,c,n3)
c
      real a(n1,5),b(5,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
         enddo
      enddo
      return
      end
      subroutine mxm6(a,n1,b,n2,c,n3)
c
      real a(n1,6),b(6,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
         enddo
      enddo
      return
      end
      subroutine mxm7(a,n1,b,n2,c,n3)
c
      real a(n1,7),b(7,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
         enddo
      enddo
      return
      end
      subroutine mxm8(a,n1,b,n2,c,n3)
c
      real a(n1,8),b(8,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
         enddo
      enddo
      return
      end
      subroutine mxm9(a,n1,b,n2,c,n3)
c
      real a(n1,9),b(9,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
         enddo
      enddo
      return
      end
      subroutine mxm10(a,n1,b,n2,c,n3)
c
      real a(n1,10),b(10,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
         enddo
      enddo
      return
      end
      subroutine mxm11(a,n1,b,n2,c,n3)
c
      real a(n1,11),b(11,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
         enddo
      enddo
      return
      end
      subroutine mxm12(a,n1,b,n2,c,n3)
c
      real a(n1,12),b(12,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
         enddo
      enddo
      return
      end
      subroutine mxm13(a,n1,b,n2,c,n3)
c
      real a(n1,13),b(13,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
         enddo
      enddo
      return
      end
      subroutine mxm14(a,n1,b,n2,c,n3)
c
      real a(n1,14),b(14,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
         enddo
      enddo
      return
      end
      subroutine mxm15(a,n1,b,n2,c,n3)
c
      real a(n1,15),b(15,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
         enddo
      enddo
      return
      end
      subroutine mxm16(a,n1,b,n2,c,n3)
c
      real a(n1,16),b(16,n3),c(n1,n3)
c
      do j=1,n3
         do i=1,n1
            c(i,j) = a(i,1)*b(1,j)
     $             + a(i,2)*b(2,j)
     $             + a(i,3)*b(3,j)
     $             + a(i,4)*b(4,j)
     $             + a(i,5)*b(5,j)
     $             + a(i,6)*b(6,j)
     $             + a(i,7)*b(7,j)
     $             + a(i,8)*b(8,j)
     $             + a(i,9)*b(9,j)
     $             + a(i,10)*b(10,j)
     $             + a(i,11)*b(11,j)
     $             + a(i,12)*b(12,j)
     $             + a(i,13)*b(13,j)
     $             + a(i,14)*b(14,j)
     $             + a(i,15)*b(15,j)
     $             + a(i,16)*b(16,j)
         enddo
      enddo
      return
      end
