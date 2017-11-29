      subroutine cvunpack(w1,w1p,y)
c
c     copy the internal cvode vector y to nek array w1 and w1p 
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real w1(lx1,ly1,lz1,lelt,*)
      real w1p
      real y(*)

      nxyz = lx1*ly1*lz1

      j = 1
        do ifield = 2,nfield
           if (ifcvfld(ifield)) then
              ntot = nxyz*nelfld(ifield)
              call copy   (w1(1,1,1,1,ifield-1),y(j),ntot)
              call bcdirsc(w1(1,1,1,1,ifield-1)) ! restore dirichlet bcs
              j = j + ntot
           endif
        enddo

      if (ifdp0dt) then
         w1p = y(j)
         j = j + 1 
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine cvpack(y,w1,w1p,ifrhs)
c
c     copy the nek array w1 and w1p to the internal cvode vector y
c     note: assumes temperature is stored in ifield=2 (only for ifdp0dt)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(*)
      real w1(lx1,ly1,lz1,lelt,*)
      real w1p
      logical ifrhs

      common /scrsf/ dtmp(lx1,ly1,lz1,lelt)


      nxyz = lx1*ly1*lz1

      if (ifrhs .and. ifdp0dt) then
         call qthermal ! computes dp0thdt  
         dd = (gamma0 - 1.)/gamma0
         dd = dd * dp0thdt
         ntot = nxyz*nelfld(2)
         call invers2(dtmp,vtrans(1,1,1,1,2),ntot)
         call cmult(dtmp,dd,ntot)
         call add2 (w1,dtmp,ntot)
      endif

      j = 1
      do ifield = 2,nfield
         if (ifcvfld(ifield)) then
            ntot = nxyz*nelfld(ifield)
            call copy (y(j),w1(1,1,1,1,ifield-1),ntot)
            if (ifrhs) call col2(y(j),tmask(1,1,1,1,ifield-1),ntot)
            j = j + ntot
         endif
      enddo

      if (ifdp0dt) then
         y(j) = w1p
         if (ifrhs) y(j) = w1p
      endif

      return
      end

