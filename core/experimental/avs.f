      real function avs_vdiff(ix,iy,iz,e,ifld,c1,c2,ncut)
c
c FE/SE methods, even with linear stabilization, are known to exhibit 
c oscillations around discontinuities. In order to reduce (or even eliminate) 
c these oscillations we introduce an artificial vicosity in a controlled way. 
c Note, this can spoil the accuracy of the resulting method for smooth
c functions and so, it should be active only around the shocks/discontinuities 
c (to reduce oscillations) and vanishes in smooth regions (to keep accuracy).
c
c The artificial viscosity is defined as
c
c     nu_a   = min(nu_max, c1 * h**2 * R) 
c     where nu_max = c2 * h * vmax 
c
c and is inspired by the work of Nazarov and co-workers,
c where a similar method is proposed. Instead of computing an entropy or
c equation residual, we use a simple gradient approach to detect spatial
c irregularities.  
c
c c1, c2 and ncut a user tunable control parameters. 
c First, set ncut to 2. Then, set c1 to a huge value make c2
c (0.5 is always a good first guess) to be as small as possible.
c Once c2 is chosen, set c2 = 1.0 and reduce/increase as much 
c possible/required.
c
c Known limitations:
c - scalars only (ifield > 1)
c
      include 'SIZE'
      include 'TOTAL'

      integer ix, iy, iz, e, ifld
      real c1, c2

      integer res_mode
      parameter (res_mode=1)

      parameter (lt=lx1*ly1*lz1*lelt)
      common /avs_rcb/  visc(lx1,ly1,lz1,lelt),
     $                  h0(lt),uf(lt),r(lt),
     $                  tx(lt),ty(lt),tz(lt)

      parameter (lm=40)
      parameter (lm2=lm*lm)
      real hpf_filter(lm2)
      real hpf_op(lx1,lx1)
      save hpf_op

      real uinf

      save    ibuild
      data    ibuild / 0 /

      real tla
      data tla    /-1./
      save tla

      if (time .eq. tla) then
         avs_vdiff = visc(ix,iy,iz,e)
         return
      endif

      ifldt  = ifield
      ifield = ifld

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      dinv = 1./ldim
      do i = 1,n
         h0(i) = bm1(i,1,1,1)**dinv
         vmax  = sqrt(vx(i,1,1,1)**2 + vy(i,1,1,1)**2 + vz(i,1,1,1)**2)
         visc(i,1,1,1) = c2*h0(i)*vmax  ! max possible viscosity
      enddo

      ! compute residual
      if (res_mode.eq.1) then
         if (ibuild.eq.0) then
           call hpf_trns_fcn(hpf_filter,ncut)
           call build_hpf_mat(hpf_op,hpf_filter,.false.)
           ibuild = ibuild + 1
         endif
         call build_hpf_fld(uf,t(1,1,1,1,ifield-1),hpf_op,lx1,lz1)
         call opgrad(tx,ty,tz,uf)
         do i = 1,n
            binv = 1/bm1(i,1,1,1)
            r(i) = tx(i)*tx(i) + ty(i)*ty(i) + tz(i)*tz(i)
            r(i) = r(i) * binv*binv
         enddo
          uinf = 1.
      endif

      ! evaluate arificial viscosity
      do i = 1,n
        ent_visc = c1 * h0(i)*h0(i) * r(i) *uinf
        visc(i,1,1,1) = min(visc(i,1,1,1),ent_visc)
      enddo

      ! make it piecewise constant across elements
      do ie = 1,nelv
        vmax=vlmax(visc(1,1,1,ie),nxyz)
        call cfill(visc(1,1,1,ie),vmax,nxyz)
      enddo

      if (loglevel.gt.2) then
        vismx = glamax(visc,n)
        vismn = glamin(visc,n)
        visav = glsc2(visc,bm1,n)/volvm1
        if (nio.eq.0) write(6,10) time,vismx,vismn,visav
 10     format(1p4e12.4,' AVM')
      endif

      ifield = ifldt

      tla = time
      avs_vdiff = visc(ix,iy,iz,e)

      return
      end
