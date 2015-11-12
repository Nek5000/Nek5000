!--------------------------------------------------------------------
! JH061914 propose name change to intpdiff since that is what is in
!          here.
!--------------------------------------------------------------------

      subroutine compute_gradients(e)
      include 'SIZE'
      include 'INPUT'
      include 'DXYZ'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      parameter (ldd=lxd*lyd*lzd)
      common /ctmp1/ ur(ldd),us(ldd),ut(ldd),ju(ldd),ud(ldd),tu(ldd)
      real ju

      integer eq,e

!     !  Compute d/dx, d/dy and d/dz of all the cons vars

      nxy1  = nx1*ny1
      nyz1  = ny1*nz1
      nxyz1 = nx1*ny1*nz1
!     m0 = nx1-1

      do eq=1,toteq

         if (if3d) then
! assumes dxm1=dym1=dzm1, nx1=ny1=nz1
!           call local_grad3(ur,us,ut,u(1,1,1,eq,e),m0,1,dxm1,dxm1t)
! eh, why do that?
            call mxm(dxm1,nx1,u(1,1,1,eq,e),nx1,ur,nyz1)
            do iz=1,nz1
               il=(iz-1)*nxy1+1
               call mxm(u(1,1,iz,eq,e),nx1,dytm1,ny1,us(il),ny1)
            enddo
            call mxm(u(1,1,1,eq,e),nxy1,dztm1,nz1,ut,nz1)
            do i=1,nxyz1
               gradu(i,1,1,eq,1) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)+
     >                                         sxm1(i,1,1,e)*us(i)+
     >                                         txm1(i,1,1,e)*ut(i))
            enddo
            do i=1,nxyz1
               gradu(i,1,1,eq,2) = jacmi(i,e)*(rym1(i,1,1,e)*ur(i)+
     >                                         sym1(i,1,1,e)*us(i)+
     >                                         tym1(i,1,1,e)*ut(i))
            enddo
            do i=1,nxyz1
               gradu(i,1,1,eq,3) = jacmi(i,e)*(rzm1(i,1,1,e)*ur(i)+
     >                                         szm1(i,1,1,e)*us(i)+
     >                                         tzm1(i,1,1,e)*ut(i))
            enddo

         else

!           call local_grad2(ur,us   ,u(1,1,1,eq,e),m0,1,dxm1,dxm1t)
            call mxm (dxm1,nx1,u(1,1,1,eq,e),nx1,ur,nyz1)
            call mxm (u(1,1,1,eq,e),nx1,dytm1,ny1,us,ny1)
            do i=1,nxyz1
               gradu(i,1,1,eq,1) = jacmi(i,e)*(rxm1(i,1,1,e)*ur(i)+
     >                                         sxm1(i,1,1,e)*us(i))
            enddo
            do i=1,nxyz1
               gradu(i,1,1,eq,2) = jacmi(i,e)*(rym1(i,1,1,e)*ur(i)+
     >                                         sym1(i,1,1,e)*us(i))
            enddo

         endif

      enddo ! equation loop

      return
      end

!-----------------------------------------------------------------------

      subroutine set_rxgll

!-----------------------------------------------------------------------
! JH091015 legacy. minimizes changes to flux_div from the bad old days
!-----------------------------------------------------------------------

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP' ! for istep
      include 'CMTDATA'
      include 'WZ'
      integer e

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
      if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep = istep

      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd

      if (if3d) then
         do e=1,nelt
            call copy(rx(1,1,e),rxm1(1,1,1,e),nxyz1)
            call copy(rx(1,2,e),rym1(1,1,1,e),nxyz1)
            call copy(rx(1,3,e),rzm1(1,1,1,e),nxyz1)
            call copy(rx(1,4,e),sxm1(1,1,1,e),nxyz1)
            call copy(rx(1,5,e),sym1(1,1,1,e),nxyz1)
            call copy(rx(1,6,e),szm1(1,1,1,e),nxyz1)
            call copy(rx(1,7,e),txm1(1,1,1,e),nxyz1)
            call copy(rx(1,8,e),tym1(1,1,1,e),nxyz1)
            call copy(rx(1,9,e),tzm1(1,1,1,e),nxyz1)
         enddo
      else
         do e=1,nelt
            call copy(rx(1,1,e),rxm1(1,1,1,e),nxyz1)
            call copy(rx(1,2,e),rym1(1,1,1,e),nxyz1)
            call copy(rx(1,3,e),sxm1(1,1,1,e),nxyz1)
            call copy(rx(1,4,e),sym1(1,1,1,e),nxyz1)
         enddo
      endif

      return
      end
