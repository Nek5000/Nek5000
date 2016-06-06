      subroutine add_fcvfun_usr(ydot,jstart)
c                                               
c     plug-in specific contributions to rhs     
c                                               
      include 'SIZE'                            
      include 'TOTAL'                          
                                               
      real ydot(1)
                                                
      ntotv = nx1*ny1*nz1*nelv                

      if (ifvcor .and. iflomach) then

         dd = gamma0   ! CVref/CPref ! Note CVref denotes the inverse CPref  
         dd = (dd - 1.)/dd
         xfacr= dd * dp0thdt

         do i = 1,ntotv
            ydot(i) = ydot(i) + xfacr/vtrans(i,1,1,1,2)
         enddo

      else
          dp0thdt = 0.0
      endif

      if (ifvcor .and. iflomach) ydot(jstart) = dp0thdt

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
      do i=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(i)
         call copy(t(1,1,1,1,i-1),y(j),ntot)
         j = j + ntot
      enddo

      if (ifvcor .and. iflomach) p0th = y(j)

      return
      end
c----------------------------------------------------------------------
      subroutine cv_pack_sol(y)
c
c     copy the internal nek array (t) to the cvode solution (y)
c
      include 'SIZE'
      include 'TOTAL'
      include 'CVODE'

      real y(1)

      nxyz = nx1*ny1*nz1

      j = 1
      do i=cv_1fld,cv_lfld
         ntot = nxyz*nelfld(i)
         call copy (y(j),t(1,1,1,1,i-1),ntot)
         j = j + ntot
      enddo

      if (ifvcor .and. iflomach) y(j) = p0th

      return
      end
