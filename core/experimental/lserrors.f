c-----------------------------------------------------------------------
      real function ls_l1norm(ifld,exact)
      include 'SIZE'
      include 'TOTAL'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      integer ntot,i
      integer ifld

      ntot = lx1*ly1*lz1*nelv

      call sub3(lstemp,t(1,1,1,1,ifld-1),exact,ntot)

      do i=1,ntot
        lstemp(i,1,1,1) = abs(lstemp(i,1,1,1))
      enddo

      call col2(lstemp,bm1,ntot)

      ls_l1norm = glsum(lstemp,ntot)

      return
      end
c-----------------------------------------------------------------------
      real function ls_l2norm(ifld,exact)
      include 'SIZE'
      include 'TOTAL'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      integer ntot,i
      integer ifld

      ntot = lx1*ly1*lz1*nelv

      call sub3(lstemp,t(1,1,1,1,ifld-1),exact,ntot)

      ls_l2norm = glsc3(lstemp,lstemp,bm1,ntot)

      ls_l2norm = sqrt(ls_l2norm)

      return
      end
c-----------------------------------------------------------------------
      real function enclosedVol(ifld,direc)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ntot,i
      real glsum

      real vol(lx1,ly1,lz1,lelt)
      integer direc,ifld

      ntot = lx1*ly1*lz1*nelv

      if(direc.lt.0)then
        call copy(vol,t(1,1,1,1,ifld-1),ntot)
        call cmult(vol,-1.0,ntot)
        call cadd(vol,1.0,ntot)
      else
        call copy(vol,t(1,1,1,1,ifld-1),ntot)
      endif

      call col2(vol,bm1,ntot)

      enclosedVol = glsum(vol,ntot)

      return
      end
