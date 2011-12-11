      subroutine add_fcvfun_usr(ydot,j)
c
      real ydot(1)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_unpack_sol(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy(t(1,1,1,1,i-1),y(j),ntot)
         j = j + ntot
      enddo

      return
      end
c----------------------------------------------------------------------
      subroutine cv_pack_sol(y)
c
c     copy the cvode solution (y) back to the internal nek array (t)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=2,cv_nfld
         ntot = nxyz*nelfld(i)
         call copy (y(j),t(1,1,1,1,i-1),ntot)
         j = j + ntot
      enddo

      return
      end
