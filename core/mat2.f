c-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      real a(1),b(1)
      include 'OPCTR'

#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)*b(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine col2c(a,b,c,n)
      real a(1),b(1),c

      do i=1,n
         a(i)=a(i)*b(i)*c
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)
      real a(1),b(1),c(1)
      include 'OPCTR'

#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'col3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

!xbm* unroll (10)
      do i=1,n
         a(i)=b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      real a(1),b(1)
      include 'OPCTR'

#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD2  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add3(a,b,c,n)
      real a(1),b(1),c(1)
      include 'OPCTR'

#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'ADD3  '
      endif
      isbcnt = N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

!xbm* unroll (10)
      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)
      real a(1),b(1),c(1)
      include 'OPCTR'

#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'addcl3'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif

!xbm* unroll (10)
      do i=1,n
         a(i)=a(i)+b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2s1(a,b,c1,n)
      real a(1),b(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s1'
      endif
      isbcnt = 2*N
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      DO 100 I=1,N
        A(I)=C1*A(I)+B(I)
  100 CONTINUE
      return
      END
C
c-----------------------------------------------------------------------
      subroutine add2s2(a,b,c1,n)
      real a(1),b(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add2s2'
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      DO 100 I=1,N
        A(I)=A(I)+C1*B(I)
  100 CONTINUE
      return
      END
C
c-----------------------------------------------------------------------
      subroutine add3s2(a,b,c,c1,c2,n)
      real a(1),b(1),c(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add3s2'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      DO 100 I=1,N
        A(I)=C1*B(I)+C2*C(I)
  100 CONTINUE
      return
      END
C
c-----------------------------------------------------------------------
      subroutine add4(a,b,c,d,n)
      REAL A(1),B(1),C(1),D(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'add4  '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      DO 100 I=1,N
         A(I)=B(I)+C(I)+D(I)
 100  CONTINUE
      return
      END
      real function vlsc2(x,y,n)
      REAL X(1),Y(1)
      include 'SIZE'
      include 'OPCTR'
      include 'PARALLEL'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      s = 0.
      do i=1,n
         s = s + x(i)*y(i)
      enddo
      vlsc2=s
      return
      end
c-----------------------------------------------------------------------
      real function vlsc21(x,y,n)
      real x(1),y(1)
      include 'SIZE'
      include 'OPCTR'
      include 'PARALLEL'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'VLSC21'
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      s = 0.
      do i=1,n
         s = s + x(i)*x(i)*y(i)
      enddo
      vlsc21=s
      return
      end
c-----------------------------------------------------------------------
      function glsc3(a,b,mult,n)
C
C     Perform inner-product in double precision
C
      REAL A(1),B(1),MULT(1)
      REAL TMP,WORK(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc3 '
      endif
      isbcnt = 3*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C
      TMP = 0.0
      DO 10 I=1,N
         TMP = TMP + A(I)*B(I)*MULT(I)
 10   CONTINUE
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC3 = TMP
      return
      END
C
c-----------------------------------------------------------------------
      function glsc2(x,y,n)
C
C     Perform inner-product in double precision
C
      real x(1), y(1)
      real tmp,work(1)
C
      include 'OPCTR'
C
#ifndef NOTIMER
      if (isclld.eq.0) then
          isclld=1
          nrout=nrout+1
          myrout=nrout
          rname(myrout) = 'glsc2 '
      endif
      isbcnt = 2*n
      dct(myrout) = dct(myrout) + (isbcnt)
      ncall(myrout) = ncall(myrout) + 1
      dcount      =      dcount + (isbcnt)
#endif
C

      tmp=0.0
      do 10 i=1,n
         tmp = tmp+ x(i)*y(i)
   10 continue
      CALL GOP(TMP,WORK,'+  ',1)
      GLSC2 = TMP
      return
      END
c-----------------------------------------------------------------------
      function glsc23(x,y,z,n)
c
C     Perform inner-product  x*x*y*z
c
      real x(1), y(1),z(1)
      real tmp,work(1)

      ds = 0.0
      do 10 i=1,n
         ds=ds+x(i)*x(i)*y(i)*z(i)
   10 continue
      tmp=ds
      call gop(tmp,work,'+  ',1)
      glsc23 = tmp
      return
      end
c-----------------------------------------------------------------------
