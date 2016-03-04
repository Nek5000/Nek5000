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

      subroutine set_dealias_face

!-----------------------------------------------------------------------
! JH111815 needed for face Jacobian and surface integrals
!-----------------------------------------------------------------------

      include 'SIZE'
      include 'INPUT' ! for if3d
      include 'GEOM'  ! for ifgeom
      include 'TSTEP' ! for istep
      include 'WZ'    ! for wxm1
      include 'DG'    ! for facewz

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
      if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep = istep

      call zwgl(zptf,wgtf,nxd)

      if (if3d) then
         k=0
         do j=1,ny1
         do i=1,nx1
            k=k+1
            wghtc(k)=wxm1(i)*wzm1(j)
         enddo
         enddo
         k=0
         do j=1,nyd
         do i=1,nxd
            k=k+1
            wghtf(k)=wgtf(i)*wgtf(j)
         enddo
         enddo
      else
         call copy(wghtc,wxm1,nx1)
         call copy(wghtf,wgtf,nxd)
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine set_alias_rx(istp)
! note that set_alias_rx will be called only when nxd = nx1
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
c     include 'TSTEP' ! for istep

      common /dealias1/ zd(lx1),wd(lx1)
      integer e

      integer ilsp
      save    ilsp
      data    ilsp /-1/

      if (.not.ifgeom.and.ilsp.gt.1) return  ! already computed
      if (ifgeom.and.ilsp.eq.istp)  return  ! already computed
      ilsp = istp
      nxyz = nx1*ny1*nz1
      call zwgll (zd,wd,nx1)  

      if (if3d)then
         do e=1,nelv

            call copy(rx(1,1,e),rxm1(1,1,1,e),nxyz) 
            call copy(rx(1,2,e),rym1(1,1,1,e),nxyz) 
            call copy(rx(1,3,e),rzm1(1,1,1,e),nxyz) 
            call copy(rx(1,4,e),sxm1(1,1,1,e),nxyz) 
            call copy(rx(1,5,e),sym1(1,1,1,e),nxyz) 
            call copy(rx(1,6,e),szm1(1,1,1,e),nxyz) 
            call copy(rx(1,7,e),txm1(1,1,1,e),nxyz) 
            call copy(rx(1,8,e),tym1(1,1,1,e),nxyz) 
            call copy(rx(1,9,e),tzm1(1,1,1,e),nxyz) 

            l = 0
            do k=1,nz1
            do j=1,ny1
            do i=1,nx1
               l = l+1
               w = wd(i)*wd(j)*wd(k)
               do ii=1,9
                  rx(l,ii,e) = w*rx(l,ii,e)
               enddo
            enddo
            enddo
            enddo
            enddo
      else

         do e=1,nelv

c           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call copy(rx(1,1,e),rxm1(1,1,1,e),nxyz) 
            call copy(rx(1,2,e),rym1(1,1,1,e),nxyz) 
            call copy(rx(1,3,e),sxm1(1,1,1,e),nxyz) 
            call copy(rx(1,4,e),sym1(1,1,1,e),nxyz) 

            l = 0
            do j=1,ny1
            do i=1,nx1
               l = l+1
               w = wd(i)*wd(j)
               do ii=1,4
                  rx(l,ii,e) = w*rx(l,ii,e)
               enddo
            enddo
            enddo
         enddo

      endif

      return
      end
!-----------------------------------------------------------------------
