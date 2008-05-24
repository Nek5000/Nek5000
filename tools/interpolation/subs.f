c-----------------------------------------------------------------------
      subroutine exitt
      stop
      end
c-----------------------------------------------------------------------
      subroutine add2(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = x(i) + y(i)
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
      subroutine sub2(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = x(i) - y(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2sq(x,y,n)
      real x(1),y(1)
      do i=1,n
         x(i) = x(i) + y(i)*y(i)
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
      subroutine rzero(x,n)
      real x(1)
      do i=1,n
         x(i) = 0.
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cmult(x,c,n)
      real x(1),c
      do i=1,n
         x(i) = x(i) * c
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine blank(x,n)
      character*1 x(1)
      do i=1,n
         x(i) = ' '
      enddo
      return
      end
c-----------------------------------------------------------------------
      function glmin(x,n)
      real x(1)
c
      s = x(1)
c
      do i=1,n
         s = min(s,x(i))
      enddo
      glmin = s
      return
      end
c-----------------------------------------------------------------------
      function glmax(x,n)
      real x(1)
c
      s = x(1)
c
      do i=1,n
         s = max(s,x(i))
      enddo
      glmax = s
      return
      end
c-----------------------------------------------------------------------
      FUNCTION LTRUNC(STRING,L)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/
C
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CHCOPY(A,B,N)
      CHARACTER*1 A(1), B(1)
      DO 100 I = 1, N
 100     A(I) = B(I)
      RETURN
      END
c-----------------------------------------------------------------------
      logical function if_byte_swap_test(bytetest)

      real*4 bytetest
      real*4 test_pattern
      save   test_pattern

      integer icalld
      save    icalld
      data    icalld /0/

      eps          = 0.00021
      test_pattern = 6.54321

      if_byte_swap_test = .false.
      if (abs(bytetestn-test_pattern).gt.eps) if_byte_swap_test = .true.


      if (abs(bytetestn-test_pattern).gt.eps) then
          call byte_reverse(bytetest,1)
c          write(6,*)'Byte swap2',if_byte_swap_test,bytetest,test_pattern
      endif

c      if (icalld.eq.0) write(6,*)'Byte swap:',if_byte_swap_test,bytetest
      icalld = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine outl(x,n,name)
      character*3 name
      real x(1)
c
      n1 = min(n,6)
      write(6,1) name,(x(k),k=1,n1)
   1  format(a3,1p6e13.5)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine stats(x,n,name)
      character*3 name
      real x(n,1)
      real ymax(5),ymin(5)
c
      do k=1,5
         ymax(k) = glmax(x(1,k),n)
         ymin(k) = glmin(x(1,k),n)
      enddo
c
      write(6,*)
      write(6,1) name,' min ',(ymin(k),k=1,5)
      write(6,1) name,' max ',(ymax(k),k=1,5)
c
   1  format(a3,a5,1p6e13.5)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dnorm (d2,di,x1,x2,n)
      real x1(n),x2(n)
c
      d2 = 0.
      di = 0.
      do i=1,n
         dif = x1(i)-x2(i)
         di  = max(di,abs(dif))
         d2  = d2 + dif*dif
      enddo
      if (n.gt.0.and.d2.gt.0) d2 = sqrt(d2)/n
      return
      end
c-----------------------------------------------------------------------
      subroutine outmat(a,m,n,name6,ie)
      real a(m,n)
      character*6 name6
c
      return
      write(6,*) 
      write(6,*) ie,' matrix: ',name6,m,n
      n12 = min(n,12)
      do i=1,m
         write(6,6) ie,name6,(a(i,j),j=1,n12)
      enddo
    6 format(i3,1x,a6,12f9.5)
      write(6,*) 
      return
      end
c-----------------------------------------------------------------------
      subroutine lex2sem_fld
     $  (v,nr,ns,nt,nelx,nely,nelz,nfld,if3d,u,mx,my,mz)
c
      logical if3d
c
      real v(0:nr,0:ns,0:nt,nfld,1)
      real u(mx,my,mz,nfld)
      integer e,ex,ey,ez
c
      l = 0
c
      do kfld=1,nfld
       if (if3d) then
         kk=2
         do ez=1,nelz
            kk=kk-1
            do k=0,nt
               jj=2
               do ey=1,nely
                  jj=jj-1
                  do j=0,ns
                     ii=2
                     do ex=1,nelx
                        ii=ii-1
                        do i=0,nr
                           e=ex+nelx*((ey-1)+nelz*(ez-1))
                           v(i,j,k,kfld,e) = u(ii,jj,kk,kfld)
                           ii=ii+1
                        enddo
                     enddo
                     jj=jj+1
                  enddo
               enddo
               kk=kk+1
            enddo
         enddo
       else
         jj=2
         do ey=1,nely
            jj = jj-1
            do j=0,ns
               ii=2
               do ex=1,nelx
                  ii = ii-1
                  do i=0,nr
                     e=ex+nelx*(ey-1)
                     v(i,j,0,kfld,e) = u(ii,jj,1,kfld)
                     ii=ii+1
                  enddo
               enddo
               jj=jj+1
            enddo
         enddo
       endif
c
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outfld_fm
     $    (uw,nx1,ny1,nz1,nelx,nely,nelz,nfld,if3d)
c
      real uw(nx1*ny1*nz1,nfld,1)
      integer e
c
      common /cexcod/ excode(10)
      character*2     excode
c
      istep = 0
      time  = 0
c
      nelt = nelx*nely
      if (if3d) nelt = nelx*nely*nelz
c
      write(11,'(4i4,1PE14.7,i5,1x,15a2,1x,a12)')
     $ nelt,nx1,ny1,nz1,time,istep,(excode(i),i=1,10)
c
      cdrror = 0.
      write(11,'(6g11.4)')(cdrror,i=1,nelt)
c
      nxyz = nx1*ny1*nz1
      do e=1,nelt
      do i=1,nxyz
         write(11,2) (uw(i,k,e),k=1,nfld)
      enddo
      enddo
    2 format(1p6e15.7)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_header(s80)
c
      common /iold/ nelx,nely,nelz,nx,ny,nz
      common /inew/ nenx,neny,nenz,nx1,ny1,nz1
c
      character*80 s80
c
      common /cexcod/ excode(10)
      character*2     excode
c
c
      nelot = nenx*neny*nenz
      if (nelot.lt.10000) then
         read (s80,'(4i4,1p1e14.7,i5,1x,15a2,1x,a12)')
     $   melot,mx1,my1,mz1,time,istep,(excode(i),i=1,10)
      else
         read (s80,'(i10,3i4,1p1e18.9,i9,1x,15a2)')
     $   melot,mx1,my1,mz1,time,istep,(excode(i),i=1,10)
      endif
c
      call blank(s80,80)
      nelnt = nenx*neny*nenz
      if (nelnt.lt.10000) then
         write(s80,'(4i4,1p1e14.7,i5,1x,15a2,1x,a12)')
     $   nelnt,nx1,ny1,nz1,time,istep,(excode(i),i=1,10),
     $   ' 4 NELT,NX,NY,N'
      else
         write(s80,'(i10,3i4,1p1e18.9,i9,1x,15a2)')
     $   nelnt,nx1,ny1,nz1,time,istep,(excode(i),i=1,10)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_data(x,n,fname)
C
      real*4 x(1)
      character*40 fname
      common /c80/ s80
      character*80 s80
c
c
      real*4 bytetest
      bytetest = 6.54321
c
c
c     Open file
c
      call byte_open(fname)
c
c     Write new header
c
      call gen_header(s80)
      call byte_write(s80,20)
      call byte_write(bytetest,1)
c
      call byte_write(x,n)
c
      call byte_close()
c
      return
      end
c-----------------------------------------------------------------------
