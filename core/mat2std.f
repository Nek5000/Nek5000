c-----------------------------------------------------------------------
      subroutine col2(a,b,n)
      dimension a(1),b(1)
      include 'OPCTR'
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


      do i=1,n
         a(i)=a(i)*b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine col3(a,b,c,n)
      dimension a(1),b(1),c(1)
      include 'OPCTR'
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


      do i=1,n
         a(i)=b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add2(a,b,n)
      dimension a(1),b(1)
      include 'OPCTR'
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

      do i=1,n
         a(i)=a(i)+b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine add3(a,b,c,n)
      dimension a(1),b(1),c(1)
      include 'OPCTR'
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


      do i=1,n
         a(i)=b(i)+c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine addcol3(a,b,c,n)
      dimension a(1),b(1),c(1)
      include 'OPCTR'
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

      do i=1,n
         a(i)=a(i) +  b(i)*c(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
