c-----------------------------------------------------------------------
      real function ls_shapeerr(ifld,exact)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      real exact(1)

      integer ntot,i
      integer ifld

      common /ls_err_arrs2/ lstemp2(lx1,ly1,lz1,lelv),
     $                      lstemp3(lx1,ly1,lz1,lelv) 
      real lstemp2,lstemp3

      real vol,vole
      real enclosedVol
      real glmax,glmin
      real tmax,tmin

      ntot = lx1*ly1*lz1*nelv

      call sharpHeaviside(t(1,1,1,1,ifld-1),lstemp2)

      call sharpHeaviside(exact,lstemp3)

      vole = enclosedVol(lstemp3,1)

      do i=1,ntot
        lstemp(i,1,1,1) = abs(lstemp2(i,1,1,1)-lstemp3(i,1,1,1))
      enddo

      vol = enclosedVol(lstemp,1)

      ls_shapeerr = vol/vole

      return
      end
c-----------------------------------------------------------------------
      real function ls_volerr(ifld,exact,initial)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      integer ntot,i
      integer ifld
      real initial(1), exact(1)

      real vol,voli,vole
      real enclosedVol

      ntot = lx1*ly1*lz1*nelv

      vol = enclosedVol(t(1,1,1,1,ifld-1),1)

      voli = enclosedVol(initial,1)

      vole = enclosedVol(exact,1)

      ls_volerr = (vol - voli)/vole

      return
      end
c-----------------------------------------------------------------------
      real function ls_relerr(ifld,exact)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      real exact(1)

      integer ntot,i
      integer ifld

      real ls_l1norm
      real glsum, tmp

      ntot = lx1*ly1*lz1*nelv

      call sub3(lstemp,exact,bm1,ntot)

      tmp = glsum(lstemp,ntot)

      ls_relerr = ls_l1norm(ifld,exact)/tmp

      return
      end
c-----------------------------------------------------------------------
      real function ls_l1norm(ifld,exact)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      real exact(1)

      integer ntot,i
      integer ifld
      real glsum

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
      implicit none
      include 'SIZE'
      include 'TOTAL'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      real exact(1)

      integer ntot,i
      integer ifld
      real glsc3

      ntot = lx1*ly1*lz1*nelv

      call sub3(lstemp,t(1,1,1,1,ifld-1),exact,ntot)

      ls_l2norm = glsc3(lstemp,lstemp,bm1,ntot)

      ls_l2norm = sqrt(ls_l2norm)

      return
      end
c-----------------------------------------------------------------------
      real function enclosedVol(psi,direc)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      integer ntot,i
      real glsum

      integer direc
      real psi(1)

      ntot = lx1*ly1*lz1*nelv

      if(direc.lt.0)then
        call copy(lstemp,psi,ntot)
        call cmult(lstemp,-1.0,ntot)
        call cadd(lstemp,1.0,ntot)
      else
        call copy(lstemp,psi,ntot)
      endif

      call col2(lstemp,bm1,ntot)

      enclosedVol = glsum(lstemp,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine sharpHeaviside(psi,heavi)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ntot,i
      real psi(1),heavi(1)

      ntot = lx1*ly1*lz1*nelv

      call copy(heavi,psi,ntot)

      call cadd(heavi,-0.5,ntot)

      do i=1,ntot
        if(heavi(i).ge.0.0)then
          heavi(i) = 1.0
        elseif(heavi(i).lt.0.0)then
          heavi(i) = 0.0
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine printVol(ifld,direc)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ifld

      integer icalld
      save icalld
      data icalld /0/

      real voli 
      save voli
      integer direc

      real vol,enclosedVol

      vol = enclosedVol(t(1,1,1,1,ifld-1),direc)

      if(icalld.eq.0)then
        voli = vol
        icalld = 1
      endif

      if(nio.eq.0)then 
        write(*,*)"Volume (initial/current/ratio):",voli,vol,vol/voli
      endif

      return
      end
c-----------------------------------------------------------------------
      real function enclosedVolInt(psi,a,direc)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      common /ls_err_arrs/ lstemp(lx1,ly1,lz1,lelv)
      real lstemp

      integer ntot,i
      real glsum

      integer direc
      real psi(1),a(1)

      ntot = lx1*ly1*lz1*nelv

      if(direc.lt.0)then
        call copy(lstemp,psi,ntot)
        call cmult(lstemp,-1.0,ntot)
        call cadd(lstemp,1.0,ntot)
      else
        call copy(lstemp,psi,ntot)
      endif

      call col2(lstemp,bm1,ntot)
      call col2(lstemp,a,ntot)

      enclosedVolInt = glsum(lstemp,ntot)

      return
      end

