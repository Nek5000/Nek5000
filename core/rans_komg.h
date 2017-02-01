      real function rans_komg_nut(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      real timel
      save timel
      data timel /-1.0/

      if(time.ne.timel) then
         call rans_komg_compute
         timel = time
      endif

      rans_komg_nut = nut(ix,iy,iz,iel)
      
      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_getCoeffs(coeffs_in)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      real coeffs_in(1)

c      if(ncoeffs_in.lt.ncoeffs) 
c     $  call exitti('dim of user provided komg coeffs array 
c     $               should be >=$',ncoeffs)

      do i=1,ncoeffs
         coeffs_in(i) = coeffs(i)
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      real function rans_komg_kSrc(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_komg_kSrc = kSrc(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_komg_omgSrc(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_komg_omgSrc = omgSrc(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_compute
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)
      common /scrns/ u_x(lxyz),u_y(lxyz),u_z(lxyz)
     $             , v_x(lxyz),v_y(lxyz),v_z(lxyz)
     $             , w_x(lxyz),w_y(lxyz),w_z(lxyz)
     $             , k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $             , o_x(lxyz),o_y(lxyz),o_z(lxyz)
      real k_x,k_y,k_z

      real           tempv(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,nu_t,mu,k_t

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

!====turbulence constants==========

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigma_k      = coeffs( 2)
        sigma_omega  = coeffs( 3)

c Low Reynolds number correction constants (see Eq. 29 in notes)
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta_i       = coeffs( 6)
        alp0_str     = coeffs( 7)

c Production of omega constants
        alp_inf      = coeffs( 8)
        alpha_0      = coeffs( 9)
        r_w          = coeffs(10)

c Dissipation of K constants
        betainf_str  = coeffs(11)
        r_b          = coeffs(12)
        akk          = coeffs(13)

c Dissipation of omega constants
c         beta_i = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

!================================

!compute over all elements

      do e=1,nelv
c     Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c  http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
!computation of G_k term later
       if(if3d)then
         call gradm11(u_x,u_y,u_z,vx,e)
         call gradm11(v_x,v_y,v_z,vy,e)
         call gradm11(w_x,w_y,w_z,vz,e)
         do i=1,lxyz
            g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $           +   (u_y(i) + v_x(i))**2
     $           +   (u_z(i) + w_x(i))**2
     $           +   (w_y(i) + v_z(i))**2
         enddo
       else 
         call gradm11(u_x,u_y,u_z,vx,e)
         call gradm11(v_x,v_y,v_z,vy,e)
         do i=1,lxyz
            g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $           +   (u_y(i) + v_x(i))**2
         enddo
       endif

       call gradm11(k_x,k_y,k_z,t(1,1,1,1,ifld_k    -1),e)
       call gradm11(o_x,o_y,o_z,t(1,1,1,1,ifld_omega-1),e)

! solve for omega_pert
        call add2 (o_x, dfdx_omegb(1,1,1,e), lxyz)
        call add2 (o_y, dfdy_omegb(1,1,1,e), lxyz)
        call add2 (o_z, dfdz_omegb(1,1,1,e), lxyz)

 
         do i=1,lxyz

c           rho = param(1) ! vtrans(i,1,1,e,1)
           rho = vtrans(i,1,1,e,1)
           mu  = param(2) ! vdiff (i,1,1,e,1)

!======= limits for k, omega and t	 
          if(t(i,1,1,e,ifld_omega-1).lt.0.0)
     $     t(i,1,1,e,ifld_omega-1) = 0.01*abs(t(i,1,1,e,ifld_omega-1))

          if(t(i,1,1,e,ifld_k-1).lt.0.0)
     $     t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))

c          if(t(i,1,1,e,1).lt.0.9*tref1)
c     $       t(i,1,1,e,1) = 0.9*tref1

c          if(t(i,1,1,e,1).gt.3*tref1)then
c          write(6,*)'High temp','in element ',e,' at istep = ',istep
c          call exitti('High Temp$',uxmax)
c!	  stop
c          endif

!       if(t(i,1,1,e,ifld_k-1).gt.200)t(i,1,1,e,ifld_k-1) = 200
!       if(t(i,1,1,e,ifld_omega-1).gt.1.2*omeg_max)
!     $  t(i,1,1,e,ifld_omega-1) = 1.2*omeg_max

!=======all limits off	 

! solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

!See equations from eqns_k_omega1.pdf from Eq. (3) onwards
!Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
!no source terms Sk or S_w are added

          re_t    = rho * k /(mu * omega + tiny) 

          alp_str = alpinf_str * (alp0_str + (re_t/r_k))
     $                                         / (1+(re_t/r_k)) 
          nu_t    = alp_str*k/(omega + tiny)
!nu_t is kinematic turbulent viscosity
!units of k = m2/s2, units of omega = 1/s, nu units = m2/s
!set limit for nu_t
!            nu_t    = max(tiny,nu_t)
!	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

            betai_str = betainf_str * (akk + (re_t/r_b)**4)
     $                / (1 + (re_t/r_b)**4)
            if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
            else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
            endif


            if (xk.le.0)then
              f_beta_str = 1
            else
              f_beta_str = (1 + 680.0*xk*xk)/(1 + 400.0*xk*xk)
            endif
            Y_k = rho * betai_str * f_beta_str * k * omega
 
c betai_str = beta_star in 12.5.15 for incompressible flow
 
            G_k = rho*nu_t*g(i)
!g(i) is S**2 in equation sheet
            kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega
            alpha = (alp_inf/alp_str) * 
     $            ((alpha_0 + (re_t/r_w))/(1 + (re_t/r_w)))

            G_w = alpha*omega*G_k/(k+tiny)

c Compute dissipation of omega
            beta = beta_i
c no compressibility correction M < 0.25

            omg(1,1) = 0.0
            omg(1,2) = 0.5*(u_y(i) - v_x(i))
            if(if3d)then
              omg(1,3) = 0.5*(u_z(i) - w_x(i))
            else
              omg(1,3) = 0
            endif    
            omg(2,1) = -omg(1,2)
            omg(2,2) = 0.0
            if(if3d)then
               omg(2,3) = 0.5*(v_z(i) - w_y(i))
            else
               omg(2,3) = 0.0
            endif
            omg(3,1) = -omg(1,3)
            omg(3,2) = -omg(2,3)
            omg(3,3) = 0.0 
c ----------
            s(1,1)  = u_x (i)
            s(1,2)  = 0.5*(v_x(i) + u_y(i))
            if(if3d)then
                s(1,3)  = 0.5*(w_x(i) + u_z(i))
            else
                s(1,3) = 0.0
            endif
            s(2,1)  = 0.5*(u_y(i) + v_x(i))
            s(2,2)  = v_y(i) 
            if(if3d)then
               s(2,3)  = 0.5*(w_y(i) + v_z(i))
            else
               s(2,3) = 0.0
            endif
            if(if3d)then
               s(3,1)  = 0.5*(u_z(i) + w_x(i))
               s(3,2)  = 0.5*(v_z(i) + w_y(i))
               s(3,3)  = w_z(i) 
            else
               s(3,1) = 0
               s(3,2) = 0
               s(3,3) = 0
            endif
   
            sum_xx = 0
            do kk = 1,3
               do jj = 1,3
                  do ii = 1,3
                     sum_xx = sum_xx + omg(ii,jj)*omg(jj,kk)*s(kk,ii)
                  enddo
               enddo
            enddo

            x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
            f_b = (1 + 70*x_w)/(1 + 80*x_w)

            Y_w = rho*beta*f_b * omega * omega


            omgSrc(i,1,1,e) = G_w - Y_w

            nut  (i,1,1,e)   = nu_t
c            dif_k  (i,1,1,e)   = rho*nu_t/sigma_k
c            dif_omega(i,1,1,e) = rho*nu_t/sigma_omega
c            dif_u  (i,1,1,e)   = rho*nu_t
c            dif_t  (i,1,1,e)   = rho*nu_t/Pr_t

!turbulent Prandtl number = 0.85 fluent manual


c            src_t  (i,1,1,e) = 0.0 ! (mu+rho*nu_t)*g(i)

! viscous dissipation term to be added to the energy equation
c            src_u  (i,1,1,e,1) = -2.*k_x(i)/3.   ! No rho here because 
c            src_u  (i,1,1,e,2) = -2.*k_y(i)/3.   ! Nek will later multiply
c            if(if3d)then
c               src_u  (i,1,1,e,3) = -2.*k_z(i)/3.   ! by rho!
c            else
c               src_u  (i,1,1,e,3) = 0.0
c            endif

         enddo
!enddo for lxyz

! solve for omega_pert

c      goto 111
      sigma_omega1 = 1.0/sigma_omega
      call copy   (mu_omeg,nut(1,1,1,e),lxyz)
      call cmult  (mu_omeg,sigma_omega1,lxyz)
      call cadd   (mu_omeg,mu          ,lxyz)
      call col3   (extra_src_omega,mu_omeg , delsqf_omegb(1,1,1,e),lxyz)
      call gradm11(mu_omegx,mu_omegy,mu_omegz,nut,e)
      call cmult  (mu_omegx,sigma_omega1,lxyz)
      call cmult  (mu_omegy,sigma_omega1,lxyz)
      call cmult  (mu_omegz,sigma_omega1,lxyz)

! need to multiply VX, VY, VZ by variable rho (vtrans) for compressible 
      call sub3   (tempv, mu_omegx, VX(1,1,1,e),lxyz)
      call addcol3(extra_src_omega,tempv   ,   dfdx_omegb(1,1,1,e),lxyz)
      call sub3   (tempv, mu_omegy, VY(1,1,1,e),lxyz)
      call addcol3(extra_src_omega,tempv   ,   dfdy_omegb(1,1,1,e),lxyz)
      call sub3   (tempv, mu_omegz, VZ(1,1,1,e),lxyz)
      call addcol3(extra_src_omega,tempv   ,   dfdz_omegb(1,1,1,e),lxyz)

      call add2   (omgSrc(1,1,1,e), extra_src_omega            ,lxyz)
c 111  continue
c      if(nid.eq.0) then
c        write(*,*) 'Skipping extra source omega e ', e
c      endif

      enddo
!enddo for elements (1,nelv)

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_init(ifld_k_in,ifld_omega_in,ifcoeffs
     $                         ,coeffs_in)
c
c     Initialize values ifld_omega & ifld_k for RANS k-omega turbulence
c     modeling
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      real coeffs_in(1)
      logical ifcoeffs

      character*3 bcw

      ifld_k     = ifld_k_in
      ifld_omega = ifld_omega_in
      if (ifld_omega.gt.ldimt1 .or. ifld_k.gt.ldimt1) 
     $  call exitti('nflds gt ldimt1, recompile with ldimt1 > ',
     $  ifld_omega)

! specify k-omega model coefficients

c      if(ncoeffs_in.lt.ncoeffs) 
c     $  call exitti('dim of user provided komg coeffs array 
c     $               should be >=$',ncoeffs)

      if(ifcoeffs) then
         do i=1,ncoeffs
            coeffs(i) =coeffs_in(i)
         enddo
      else
         call rans_komg_set_defaultcoeffs
      endif

! solve for omega_pert
      bcw = 'W  '
      if(nid.eq.0) write(*,*) 'BC for distance ', bcw
      call cheap_dist(ywd,1,bcw)
      call rans_komg_omegabase

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_omegabase
c
c     Compute Omega base solution which diverges at walls
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      integer i,j,k,e
      real    nu0

      real dudx(lx1,ly1,lz1,lelv), dudy(lx1,ly1,lz1,lelv)
     $   , dudz(lx1,ly1,lz1,lelv), temt(lx1,ly1,lz1,lelv)


      real kv_min,omeg_max

      omeg_max   = coeffs(15)
      beta0      = coeffs(6)
      nu0  = param(2)
      ntot1 = nx1*ny1*nz1*nelv

      Cfcon = 6.0 * nu0 / beta0
      expn  = -2.0

      call gradm1 (dfdx_omegb,dfdy_omegb,dfdz_omegb,   ywd)
      call opcolv (dfdx_omegb,dfdy_omegb,dfdz_omegb,   bm1)
      call opdssum(dfdx_omegb,dfdy_omegb,dfdz_omegb)
      call opcolv (dfdx_omegb,dfdy_omegb,dfdz_omegb,binvm1)

      call gradm1 (dudx,      dudy      ,dudz,  dfdx_omegb)
      call copy   (delsqf_omegb, dudx, ntot1)
      call gradm1 (dudx,      dudy      ,dudz,  dfdy_omegb)
      call add2   (delsqf_omegb, dudy, ntot1)
      call gradm1 (dudx,      dudy      ,dudz,  dfdz_omegb)
      call add2   (delsqf_omegb, dudz, ntot1)

!      call opgrad (dfdx_omegb,dfdy_omegb,dfdz_omegb,   ywd)
!      call opdssum(dfdx_omegb,dfdy_omegb,dfdz_omegb)
!      call opcolv (dfdx_omegb,dfdy_omegb,dfdz_omegb,binvm1)

!      call opdiv  (delsqf_omegb,dfdx_omegb,dfdy_omegb,dfdz_omegb)
!      call dssum  (delsqf_omegb,nx1,ny1,nz1)
!      call col2   (delsqf_omegb,binvm1,ntot1)

      call vdot3  (delfsq_omegb,dfdx_omegb,dfdy_omegb,dfdz_omegb 
     $                         ,dfdx_omegb,dfdy_omegb,dfdz_omegb,ntot1)
 
      do e = 1,nelv
      do k = 1,nz1
      do j = 1,ny1
      do i = 1,nx1

         yw   = ywd(i,j,k,e)       ! 1.0 - abs(y)
         ywmin=sqrt(Cfcon/omeg_max)
         if(yw.gt.ywmin) then
            ywm2 = 1.0/(yw*yw)
            ywm3 = ywm2/yw
         else
            ywm2 = omeg_max / Cfcon
            ywm3 = ywm2/ywmin
         endif
         ywm4 = ywm2*ywm2

         f_omegb     (i,j,k,e) =        Cfcon * ywm2
         delfpart              = expn * Cfcon * ywm3
         dfdx_omegb  (i,j,k,e) = delfpart * dfdx_omegb(i,j,k,e) 
         dfdy_omegb  (i,j,k,e) = delfpart * dfdy_omegb(i,j,k,e) 
         dfdz_omegb  (i,j,k,e) = delfpart * dfdz_omegb(i,j,k,e) 
         delsqfpart            = expn * (expn - 1.0) * Cfcon * ywm4
         delsqf_omegb(i,j,k,e) = delsqfpart  * delfsq_omegb(i,j,k,e)
     $                         + delfpart    * delsqf_omegb(i,j,k,e)
c         t(i,j,k,e,1) = f_omegb(i,j,k,e)
c         t(i,j,k,e,2) = ywd(i,j,k,e)
c!         t(i,j,k,e,2) = delfsq_omegb(i,j,k,e)
c         t(i,j,k,e,3) = delsqf_omegb(i,j,k,e)

      enddo
      enddo
      enddo
      enddo

!      call outpost2(dfdx_omegb,dfdy_omegb,dfdz_omegb,pr ,t,3,'yom')
!      stop

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_set_defaultcoeffs
c
      include 'SIZE'
      include 'RANS_KOMG'

!=====various problem-specific turbulence constants
c omeg_max = value of omega on the walls
c kv_min = value of K on the walls
c Pr_t is the turbulent prandtl number

c Turbulent viscosity constants
        Pr_t         = 0.85
        coeffs( 1)   = Pr_t
        sigma_k      = 2.0
        coeffs( 2)   = sigma_k
        sigma_omega  = 2.0
        coeffs( 3)   = sigma_omega

c Low Reynolds number correction constants (see Eq. 29 in notes)
        alpinf_str   = 1.0
        coeffs( 4)   = alpinf_str
        r_k          = 6.0
        coeffs( 5)   = r_k
        beta_i       = 0.072
        coeffs( 6)   = beta_i
        alp0_str     = beta_i/3.0
        coeffs( 7)   = alp0_str

c Production of omega constants
        alp_inf      = 0.52
        coeffs( 8)   = alp_inf
        alpha_0      = 1.0/9.0
        coeffs( 9)   = alpha_0
        r_w          = 2.95
        coeffs(10)   = r_w

c Dissipation of K constants
        betainf_str  = 0.09
        coeffs(11)   = betainf_str
        r_b          = 8.0
        coeffs(12)   = r_b
        akk          = 4.0/15.0
        coeffs(13)   = akk

c Dissipation of omega constants
c         beta_i = defined earlier	     

        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 200000.0
        coeffs(15)   = omeg_max
        tiny         = 1.0e-20
        coeffs(16)   = tiny

        return
        end
