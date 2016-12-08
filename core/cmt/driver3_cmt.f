      subroutine compute_primitive_vars
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'CMTDATA'
      include 'SOLN'
      include 'DEALIAS' ! until we are comfortable with setup_convect

      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ energy(lx1,ly1,lz1),scr(lx1,ly1,lz1)
      integer e, eq

      nxyz= nx1*ny1*nz1
      ntot=nxyz*nelt

      do e=1,nelt
         call invcol3(vx(1,1,1,e),u(1,1,1,irpu,e),u(1,1,1,irg,e),nxyz)
         call invcol3(vy(1,1,1,e),u(1,1,1,irpv,e),u(1,1,1,irg,e),nxyz)
!        if (if3d)
         call invcol3(vz(1,1,1,e),u(1,1,1,irpw,e),u(1,1,1,irg,e),nxyz)
! first kinetic energy
         if (if3d) then
            call vdot3(scr,
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),
     >             u(1,1,1,irpu,e),u(1,1,1,irpv,e),u(1,1,1,irpw,e),nxyz)
         else
            call vdot2(scr,u(1,1,1,irpu,e),u(1,1,1,irpv,e),
     >                     u(1,1,1,irpu,e),u(1,1,1,irpv,e),nxyz)
         endif
         call invcol2(scr,u(1,1,1,irg,e),nxyz)
         call cmult(scr,0.5,nxyz)
! then to internal energy
         call sub3(energy,u(1,1,1,iret,e),scr,nxyz)
! now mass-specific
         call invcol2(energy,u(1,1,1,irg,e),nxyz)
! don't forget to get density where it belongs
         call invcol3(vtrans(1,1,1,e,irho),u(1,1,1,irg,e),phig(1,1,1,e),
     >                nxyz)
         call tdstate(e,energy)
      enddo

! setup_convect has the desired effect
! if IFPART=F
! if IFCHAR=F
! if IFCONS=T
! if igeom .ne. 1
! if param(99) .ge. 0
!-----------------------------------------------------------------------
!     call setup_convect(0)
!-----------------------------------------------------------------------
! to make life easier until we master this stuff and harmonize even better with nek,
! I'm including 'DEALIAS' and calling set_convect_cons here
      if (nxd.gt.nx1) then
         call set_convect_cons (vxd,vyd,vzd,vx,vy,vz)
      else
         call copy(vxd,vx,ntot) 
         call copy(vyd,vy,ntot) 
         call copy(vzd,vz,ntot) 
      endif

      return
      end
!-----------------------------------------------------------------------

      subroutine tdstate(e,energy)!,energy)
c compute the gas properties. We will have option to add real gas models
c We have perfect gas law. Cvg is stored full field
      include 'SIZE'
      include 'CMTDATA'
      include 'SOLN'
      include 'PARALLEL'
      include 'NEKUSE'
      integer   e,eg
      real energy(nx1,ny1,nz1)

      eg = lglel(e)
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         call nekasgn(i,j,k,e)
         e_internal=energy(i,j,k) !nekasgn should do this, but can't
         call userEOS(i,j,k,eg)
         vtrans(i,j,k,e,icp)= cp*rho
         vtrans(i,j,k,e,icv)= cv*rho
         t(i,j,k,e,1)       = temp
         pr(i,j,k,e)        = pres
         csound(i,j,k,e)    = asnd
      enddo
      enddo
      enddo

      return
      end
