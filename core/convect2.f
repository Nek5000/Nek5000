c-----------------------------------------------------------------------
      subroutine set_char_mask(mask,u,v,w) ! mask for hyperbolic system 

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'TSTEP'
      integer msk(0:1)
      character      cb*3
      parameter (lxyz1=lx1*ly1*lz1)
      common /ctmp1/ work(lxyz1,lelt)

      real mask(lxyz1,1),u(lxyz1,1),v(lxyz1,1),w(lxyz1,1)

      integer e,f

      nfaces= 2*ldim
      ntot1 = lx1*ly1*lz1*nelv
      call rzero (work,ntot1)
      call rone  (mask,NTOT1)

      ifldv = 1
      do 100 e=1,nelv
      do 100 f=1,nfaces
         cb=cbc(f,e,ifldv)
         if (cb(1:1).eq.'v' .or. cb(1:1).eq.'V' .or.
     $       cb.eq.'mv ' .or. cb.eq.'MV ') then

           call faccl3 (work(1,e),u(1,e),unx(1,1,f,e),f)
           call faddcl3(work(1,e),v(1,e),uny(1,1,f,e),f)
           if (if3d) 
     $     call faddcl3(work(1,e),w(1,e),unz(1,1,f,e),f)

           call fcaver (vaver,work,e,f)

           if (vaver.lt.0) call facev (mask,e,f,0.0,lx1,ly1,lz1)
         endif
         if (cb(1:2).eq.'ws' .or. cb(1:2).eq.'WS') 
     $   call facev (mask,e,f,0.0,lx1,ly1,lz1)
 100  continue
      call dsop(mask,'MUL',lx1,ly1,lz1)

      return
      end
c-----------------------------------------------------------------------
      subroutine set_dealias_rx
C
C     Eulerian scheme, add convection term to forcing function
C     at current time step.
C
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'TSTEP' ! for istep

      common /dealias1/ zd(lxd),wd(lxd)
      integer e

      integer ilstep
      save    ilstep
      data    ilstep /-1/

      if (.not.ifgeom.and.ilstep.gt.1) return  ! already computed
      if (ifgeom.and.ilstep.eq.istep)  return  ! already computed
      ilstep = istep

      nxyz1 = lx1*ly1*lz1
      nxyzd = lxd*lyd*lzd

      call zwgl (zd,wd,lxd)  ! zwgl -- NOT zwgll!

      if (if3d) then
c
         do e=1,nelv

c           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),rzm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sxm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,5,e),sym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,6,e),szm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,7,e),txm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,8,e),tym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,9,e),tzm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd

            l = 0
            do k=1,lzd
            do j=1,lyd
            do i=1,lxd
               l = l+1
               w = wd(i)*wd(j)*wd(k)
               do ii=1,9
                  rx(l,ii,e) = w*rx(l,ii,e)
               enddo
            enddo
            enddo
            enddo
         enddo

      else ! 2D
c
         do e=1,nelv

c           Interpolate z+ and z- into fine mesh, translate to r-s-t coords

            call intp_rstd(rx(1,1,e),rxm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,2,e),rym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,3,e),sxm1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd
            call intp_rstd(rx(1,4,e),sym1(1,1,1,e),lx1,lxd,if3d,0) ! 0 --> fwd

            l = 0
            do j=1,lyd
            do i=1,lxd
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
c-----------------------------------------------------------------------

