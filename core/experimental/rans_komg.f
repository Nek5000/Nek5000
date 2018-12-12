      real function rans_komg_mut(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      real timel
      save timel
      data timel /-1.0/

      if(time.ne.timel) then
         if(ifrans_komg_stndrd)    call rans_komg_stndrd_compute
         if(ifrans_komg_lowRe)     call rans_komg_lowRe_compute
         if(ifrans_komgSST_stndrd) call rans_komgSST_stndrd_compute
         if(ifrans_komgSST_lowRe)  call rans_komgSST_lowRe_compute
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_compute_noreg
         if(ifrans_komg_lowRe_noreg)call rans_komg_lowRe_compute_noreg
         timel = time
      endif

      rans_komg_mut = mut(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_komg_mutsk(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_komg_mutsk = mutsk(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_komg_mutso(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_komg_mutso = mutso(ix,iy,iz,iel)

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
      subroutine rans_komg_stndrd_compute
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigma_k      = coeffs( 2)
        sigma_omega  = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta_0       = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        betainf_str  = coeffs( 8)
        alp_inf      = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv
c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d) 
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom (i)= 1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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

! calculate del k * del omega / omega

          if(if3d)then
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
          else
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
          endif

          alp_str = alpinf_str 
          mu_t    = rho * alp_str*k/(omega + tiny) ! should multiply these by rho!!!
          rhoalpk(i)= rho*alp_str*k

          factr = 1.0
          sigma_omega1 = 1.0/sigma_omega
          rhoalpfr(i) = expn * rho * alp_str * factr * omwom(i)
     $                                           * sigma_omega1

c nu_t is kinematic turbulent viscosity
c units of k = m2/s2, units of omega = 1/s, nu units = m2/s
c set limit for nu_t
c           nu_t    = max(tiny,nu_t)
c	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

          betai_str = betainf_str

          if (xk3.le.0)then
            f_beta_str = 1.0
          else
            f_beta_str = (1.0 + 680.0*xk3*xk3)/(1.0 + 400.0*xk3*xk3)
          endif
          Y_k = rho * betai_str * f_beta_str * k * omega
 
c betai_str = beta_star in 12.5.15 for incompressible flow
 
          G_k = mu_t*g(i)
c g(i) is S**2 in equation sheet
          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega
          alpha = (alp_inf/alp_str)

          G_w = alpha*alp_str*rho*g(i)

c Compute dissipation of omega
          beta = beta_0
c no compressibility correction M < 0.25

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
          if(if3d) f_b = (1.0 + 70.0*x_w)/(1.0 + 80.0*x_w)

          Y_w = rho*beta*f_b * omega * omega

          omgSrc(i,1,1,e) = G_w - Y_w

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

c solve for omega_pert

c add mut*delsqf
        sigma_omega1 = 1.0/sigma_omega
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigma_omega1,lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call cmult  (tempv,mu,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
        expm = -expn
        call col3   (term2,omwom,delfsq_omegb(1,1,1,e),   lxyz)
        call col2   (term2           ,t(1,1,1,e,ifld_k-1),lxyz)
        call cmult  (term2,expm,lxyz)
        call add3   (tempv, term1, term2, lxyz)

c  form 3: -(omw/om) k (del_yw/yw) \del_omp/omw
        call col3   (term3     ,omwom,t(1,1,1,e,ifld_k-1),lxyz)
        call invcol2(term3     ,f_omegb(1,1,1,e)         ,lxyz)
        call col3   (term4,  dfdx_omegb(1,1,1,e),omp_x,   lxyz)
        call addcol3(term4,  dfdy_omegb(1,1,1,e),omp_y,   lxyz)
        if(if3d)
     $  call addcol3(term4,  dfdz_omegb(1,1,1,e),omp_z,   lxyz)
        call subcol3(tempv, term3, term4, lxyz)

        call addcol3(extra_src_omega, rhoalpfr, tempv,    lxyz)

c add rho v * del_omw
        call col3   (tempv,dfdx_omegb(1,1,1,e),VX(1,1,1,e),lxyz)
        call addcol3(tempv,dfdy_omegb(1,1,1,e),VY(1,1,1,e),lxyz)
        if(if3d)
     $  call addcol3(tempv,dfdz_omegb(1,1,1,e),VZ(1,1,1,e),lxyz)
        call col2   (tempv,vtrans(1,1,1,e,1),lxyz)
        call admcol3(extra_src_omega,tempv,f_omegb(1,1,1,e),expm,lxyz)

        call add2   (omgSrc(1,1,1,e), extra_src_omega           ,lxyz)

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_lowRe_compute
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigma_k      = coeffs( 2)
        sigma_omega  = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta_0       = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        betainf_str  = coeffs( 8)
        alp_inf      = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv
c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d) 
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom (i)= 1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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

! calculate del k * del omega / omega

          if(if3d)then
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
          else
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
          endif

          re_t    = rho * k /(mu * omega + tiny) 
          alp_str = alpinf_str * (alp0_str + (re_t/r_k))
     $                                         / (1+(re_t/r_k)) 
          mu_t    = rho * alp_str*k/(omega + tiny) ! should multiply these by rho!!!
          rhoalpk(i)= rho*alp_str*k

          rtld = re_t / r_k
          funcr = (1.0 - alp0_str) / (alp0_str + rtld)/(1.0 + rtld)
          factr = (1.0 + rtld * funcr)
          sigma_omega1 = 1.0/sigma_omega
          rhoalpfr(i) = expn * rho * alp_str * factr * omwom(i)
     $                                           * sigma_omega1

c nu_t is kinematic turbulent viscosity
c units of k = m2/s2, units of omega = 1/s, nu units = m2/s
c set limit for nu_t
c           nu_t    = max(tiny,nu_t)
c	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

          betai_str = betainf_str * (akk + (re_t/r_b)**4)
     $              / (1.0 + (re_t/r_b)**4)


          if (xk3.le.0)then
            f_beta_str = 1.0
          else
            f_beta_str = (1.0 + 680.0*xk3*xk3)/(1.0 + 400.0*xk3*xk3)
          endif
          Y_k = rho * betai_str * f_beta_str * k * omega
 
c betai_str = beta_star in 12.5.15 for incompressible flow
 
          G_k = mu_t*g(i)
c g(i) is S**2 in equation sheet
          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega
          alpha = (alp_inf/alp_str) * 
     $          ((alpha_0 + (re_t/r_w))/(1.0 + (re_t/r_w)))

          G_w = alpha*alp_str*rho*g(i)

c Compute dissipation of omega
          beta = beta_0
c no compressibility correction M < 0.25

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
          if(if3d) f_b = (1.0 + 70.0*x_w)/(1.0 + 80.0*x_w)

          Y_w = rho*beta*f_b * omega * omega

          omgSrc(i,1,1,e) = G_w - Y_w

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

c solve for omega_pert

c add mut*delsqf
        sigma_omega1 = 1.0/sigma_omega
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigma_omega1,lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call cmult  (tempv,mu,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
        expm = -expn
        call col3   (term2,omwom,delfsq_omegb(1,1,1,e),   lxyz)
        call col2   (term2           ,t(1,1,1,e,ifld_k-1),lxyz)
        call cmult  (term2,expm,lxyz)
        call add3   (tempv, term1, term2, lxyz)

c  form 3: -(omw/om) k (del_yw/yw) \del_omp/omw
        call col3   (term3     ,omwom,t(1,1,1,e,ifld_k-1),lxyz)
        call invcol2(term3     ,f_omegb(1,1,1,e)         ,lxyz)
        call col3   (term4,  dfdx_omegb(1,1,1,e),omp_x,   lxyz)
        call addcol3(term4,  dfdy_omegb(1,1,1,e),omp_y,   lxyz)
        if(if3d)
     $  call addcol3(term4,  dfdz_omegb(1,1,1,e),omp_z,   lxyz)
        call subcol3(tempv, term3, term4, lxyz)

        call addcol3(extra_src_omega, rhoalpfr, tempv,    lxyz)

c add rho v * del_omw
        call col3   (tempv,dfdx_omegb(1,1,1,e),VX(1,1,1,e),lxyz)
        call addcol3(tempv,dfdy_omegb(1,1,1,e),VY(1,1,1,e),lxyz)
        if(if3d)
     $  call addcol3(tempv,dfdz_omegb(1,1,1,e),VZ(1,1,1,e),lxyz)
        call col2   (tempv,vtrans(1,1,1,e,1),lxyz)
        call admcol3(extra_src_omega,tempv,f_omegb(1,1,1,e),expm,lxyz)

        call add2   (omgSrc(1,1,1,e), extra_src_omega           ,lxyz)

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komgSST_stndrd_compute
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c====turbulence constants==========

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigk1        = coeffs( 2)
        sigom1       = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta1        = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        beta_str     = coeffs( 8)
        gamma1       = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv

c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d) 
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom (i)= 1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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


          sum_Om = 0
          sum_St = 0
          do jj = 1,3
          do ii = 1,3
             sum_Om = sum_Om + omg(ii,jj)*omg(ii,jj)
             sum_St = sum_St + s  (ii,jj)*s  (ii,jj)
          enddo              
          enddo               

          St_magn = sqrt(2.0*sum_St)
          Om_magn = sqrt(2.0*sum_Om)

! calculate del k * del omega / omega

          if(if3d)then
          xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                          / (omega+tiny)
          else
          xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                          / (omega+tiny)
          endif

          rhoalpk(i)= rho*k
          rhoalpfr(i) = expn * rho * omwom(i) * sigom1

! calculate F2 based on arg2

          yw   = ywd  (i,1,1,e)
          ywm1 = ywdm1(i,1,1,e)
          ywm2 = ywm1*ywm1
          arg2_1 =      sqrt(abs(k)) * ywm1 / omega / beta_str
          arg2_2 =          500.0*nu * ywm2 / omega
          arg2   = max(2.0*arg2_1, arg2_2)
          Fun2   = tanh(arg2 * arg2)

! calculate F1 based on arg1

          argCD  =     2.0 * rho * sigom2 * xk
          CDkom  = max(argCD, 1.0e-10)
          arg1_1 = max(    arg2_1, arg2_2)
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / CDkom
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

! calculate mu_t 
          argn_1 = alp1*(omega + tiny)
          argn_2 = Fun2*St_magn ! this can also be St_magn

          mu_t = 0.0
          denom = max(argn_1, argn_2)
          if(yw.ne.0) mu_t = rho * alp1 * k / denom 

! Compute Y_k = dissipation of k

          Y_k = rho * beta_str * k * omega

! Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)

          G_k = mu_t*g(i)
          G_k = min(G_k, 10.0*Y_k)

! Compute Source term for k 

          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega

          beta  = Fun1 * beta1  + (1.0 - Fun1) * beta2
          gamma = Fun1 * gamma1 + (1.0 - Fun1) * gamma2
          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

! Compute production of omega

          G_w = rho * gamma * g(i)

! Compute dissipation of omega

          Y_w = rho * beta * omega * omega

! Compute additional SST term for omega

          S_w = (1.0 - Fun1) * argCD

! Compute Source term for omega

          omgSrc(i,1,1,e) = G_w - Y_w + S_w
          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t * sigk
          mutso(i,1,1,e)   = mu_t * sigom

        enddo

c solve for omega_pert

c add mut*delsqf
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigom1, lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call cmult  (tempv,mu,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
        expm = -expn
        call col3   (term2,omwom,delfsq_omegb(1,1,1,e),   lxyz)
        call col2   (term2           ,t(1,1,1,e,ifld_k-1),lxyz)
        call cmult  (term2,expm,lxyz)
        call add3   (tempv, term1, term2, lxyz)

c  form 3: -(omw/om) k (del_yw/yw) \del_omp/omw
        call col3   (term3     ,omwom,t(1,1,1,e,ifld_k-1),lxyz)
        call invcol2(term3     ,f_omegb(1,1,1,e)         ,lxyz)
        call col3   (term4,  dfdx_omegb(1,1,1,e),omp_x,   lxyz)
        call addcol3(term4,  dfdy_omegb(1,1,1,e),omp_y,   lxyz)
        if(if3d)
     $  call addcol3(term4,  dfdz_omegb(1,1,1,e),omp_z,   lxyz)
        call subcol3(tempv, term3, term4, lxyz)

        call addcol3(extra_src_omega, rhoalpfr, tempv,    lxyz)

c add rho v * del_omw
        call col3   (tempv,dfdx_omegb(1,1,1,e),VX(1,1,1,e),lxyz)
        call addcol3(tempv,dfdy_omegb(1,1,1,e),VY(1,1,1,e),lxyz)
        if(if3d)
     $  call addcol3(tempv,dfdz_omegb(1,1,1,e),VZ(1,1,1,e),lxyz)
        call col2   (tempv,vtrans(1,1,1,e,1),lxyz)
        call admcol3(extra_src_omega,tempv,f_omegb(1,1,1,e),expm,lxyz)

        call add2   (omgSrc(1,1,1,e), extra_src_omega           ,lxyz)

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komgSST_lowRe_compute
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c====turbulence constants==========

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigk1        = coeffs( 2)
        sigom1       = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta1        = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        beta_str     = coeffs( 8)
        alp_inf      = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv

c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d) 
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom (i)= 1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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


          sum_Om = 0
          sum_St = 0
          do jj = 1,3
          do ii = 1,3
             sum_Om = sum_Om + omg(ii,jj)*omg(ii,jj)
             sum_St = sum_St + s  (ii,jj)*s  (ii,jj)
          enddo              
          enddo               

          St_magn = sqrt(2.0*sum_St)
          Om_magn = sqrt(2.0*sum_Om)

          sum_xx = 0
          do kk = 1,3
             do jj = 1,3
                do ii = 1,3
                   sum_xx = sum_xx + omg(ii,jj)*omg(jj,kk)*s(kk,ii)
                enddo
             enddo
          enddo

! calculate del k * del omega / omega

          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega+tiny)
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega+tiny)
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
          endif

c ------------ begin low Re part of k-omega -------------------

          re_t    = rho * k /(mu * omega + tiny) 
          alp_str = alpinf_str * (alp0_str + (re_t/r_k))
     $                                         / (1+(re_t/r_k)) 
          rhoalpk(i)= rho*alp_str*k

          rtld = re_t / r_k
          funcr = (1.0 - alp0_str) / (alp0_str + rtld)/(1.0 + rtld)
          factr = (1.0 + rtld * funcr)
          rhoalpfr(i) = expn * rho * alp_str * factr * omwom(i) * sigom1

c nu_t is kinematic turbulent viscosity
c units of k = m2/s2, units of omega = 1/s, nu units = m2/s
c set limit for nu_t
c           nu_t    = max(tiny,nu_t)
c	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

          betai_str = beta_str * (akk + (re_t/r_b)**4)
     $              / (1.0 + (re_t/r_b)**4)

          if (xk3.le.0)then
            f_beta_str = 1.0
          else
            f_beta_str = (1.0 + 680.0*xk3*xk3)/(1.0 + 400.0*xk3*xk3)
          endif

c Compute production of omega
          alpha = (alp_inf/alp_str) *
     $          ((alpha_0 + (re_t/r_w))/(1.0 + (re_t/r_w)))

          gammai = alpha*alp_str

c Compute dissipation of omega

          x_w = abs((sum_xx)/(beta_str*omega + tiny)**3)
          f_b = 1.0
          if(if3d) f_b = (1.0 + 70.0*x_w)/(1.0 + 80.0*x_w)
          betai= beta1 * f_b

c ------------ end   low Re part of k-omega -------------------

! calculate F2 based on arg2

          yw   = ywd  (i,1,1,e)
          ywm1 = ywdm1(i,1,1,e)
          ywm2 = ywm1*ywm1
          arg2_1 =      sqrt(abs(k)) * ywm1 / omega / beta_str
          arg2_2 =          500.0*nu * ywm2 / omega
          arg2   = max(2.0*arg2_1, arg2_2)
          Fun2   = tanh(arg2 * arg2)

! calculate F1 based on arg1

          argCD  =     2.0 * rho * sigom2 * xk
          CDkom  = max(argCD, 1.0e-10)
          arg1_1 = max(    arg2_1, arg2_2)
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / CDkom
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

! calculate mu_t 
          argn_1 = alp1*(omega + tiny) / alp_str
          argn_2 = Fun2*St_magn ! this can also be St_magn

          mu_t = 0.0
          denom = max(argn_1, argn_2)
          if(yw.ne.0) mu_t = rho * alp1 * k / denom 

! Compute Y_k = dissipation of k

c          Y_k = rho * beta_str * k * omega
          Y_k = rho * betai_str * f_beta_str * k * omega

! Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)

          G_k = mu_t*g(i)
          G_k = min(G_k, 10.0*Y_k)

! Compute Source term for k 

          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega

          beta  = Fun1 * betai  + (1.0 - Fun1) * beta2
          gamma = Fun1 * gammai + (1.0 - Fun1) * gamma2
          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

! Compute production of omega

          G_w = rho * gamma * g(i)

! Compute dissipation of omega

          Y_w = rho * beta * omega * omega

! Compute additional SST term for omega

          S_w = (1.0 - Fun1) * argCD

! Compute Source term for omega

          omgSrc(i,1,1,e) = G_w - Y_w + S_w
          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t * sigk
          mutso(i,1,1,e)   = mu_t * sigom

        enddo

c solve for omega_pert

c add mut*delsqf
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigom1, lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call cmult  (tempv,mu,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
        expm = -expn
        call col3   (term2,omwom,delfsq_omegb(1,1,1,e),   lxyz)
        call col2   (term2           ,t(1,1,1,e,ifld_k-1),lxyz)
        call cmult  (term2,expm,lxyz)
        call add3   (tempv, term1, term2, lxyz)

c  form 3: -(omw/om) k (del_yw/yw) \del_omp/omw
        call col3   (term3     ,omwom,t(1,1,1,e,ifld_k-1),lxyz)
        call invcol2(term3     ,f_omegb(1,1,1,e)         ,lxyz)
        call col3   (term4,  dfdx_omegb(1,1,1,e),omp_x,   lxyz)
        call addcol3(term4,  dfdy_omegb(1,1,1,e),omp_y,   lxyz)
        if(if3d)
     $  call addcol3(term4,  dfdz_omegb(1,1,1,e),omp_z,   lxyz)
        call subcol3(tempv, term3, term4, lxyz)

        call addcol3(extra_src_omega, rhoalpfr, tempv,    lxyz)

c add rho v * del_omw
        call col3   (tempv,dfdx_omegb(1,1,1,e),VX(1,1,1,e),lxyz)
        call addcol3(tempv,dfdy_omegb(1,1,1,e),VY(1,1,1,e),lxyz)
        if(if3d)
     $  call addcol3(tempv,dfdz_omegb(1,1,1,e),VZ(1,1,1,e),lxyz)
        call col2   (tempv,vtrans(1,1,1,e,1),lxyz)
        call admcol3(extra_src_omega,tempv,f_omegb(1,1,1,e),expm,lxyz)

        call add2   (omgSrc(1,1,1,e), extra_src_omega           ,lxyz)

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_stndrd_compute_noreg
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigma_k      = coeffs( 2)
        sigma_omega  = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta_0       = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        betainf_str  = coeffs( 8)
        alp_inf      = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv
c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)
          o_y(i)= omp_y(i)
          if(if3d) 
     $    o_z(i)= omp_z(i)

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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

! calculate del k * del omega / omega

          if(if3d)then
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
          else
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
          endif

          alp_str = alpinf_str 
          mu_t    = rho * alp_str*k/(omega + tiny) ! should multiply these by rho!!!

c nu_t is kinematic turbulent viscosity
c units of k = m2/s2, units of omega = 1/s, nu units = m2/s
c set limit for nu_t
c           nu_t    = max(tiny,nu_t)
c	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

          betai_str = betainf_str

          if (xk3.le.0)then
            f_beta_str = 1.0
          else
            f_beta_str = (1.0 + 680.0*xk3*xk3)/(1.0 + 400.0*xk3*xk3)
          endif
          Y_k = rho * betai_str * f_beta_str * k * omega
 
c betai_str = beta_star in 12.5.15 for incompressible flow
 
          G_k = mu_t*g(i)
c g(i) is S**2 in equation sheet
          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega
          alpha = (alp_inf/alp_str)

          G_w = alpha*alp_str*rho*g(i)

c Compute dissipation of omega
          beta = beta_0
c no compressibility correction M < 0.25

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
          if(if3d) f_b = (1.0 + 70.0*x_w)/(1.0 + 80.0*x_w)

          Y_w = rho*beta*f_b * omega * omega

          omgSrc(i,1,1,e) = G_w - Y_w

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_lowRe_compute_noreg
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

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)
      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,term5(lxyz)
      common /scruz/ g(lxyz)

      real kv_min,omeg_max

      integer e
      real omg(3,3),s(3,3)
      real k,mu_t,nu_t,mu,nu

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

c Turbulent viscosity constants
        Pr_t         = coeffs( 1)
        sigma_k      = coeffs( 2)
        sigma_omega  = coeffs( 3)

c Low Reynolds number correction constants
        alpinf_str   = coeffs( 4)
        r_k          = coeffs( 5)
        beta_0       = coeffs( 6)
        alp0_str     = coeffs( 7)

c Dissipation of K constants
        betainf_str  = coeffs( 8)
        alp_inf      = coeffs( 9)
        r_b          = coeffs(10)
        akk          = coeffs(11)

c Production of omega constants
        alpha_0      = coeffs(12)
        r_w          = coeffs(13)

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)
c================================

      nome_neg = 0
      nkey_neg = 0
      xome_neg = 0.0
      xkey_neg = 0.0

      do e=1,nelv
c Evaluate G from eq. 4.2.15 cms.gre.ac.uk used in 
c http://staffweb.cms.gre.ac.uk/~ct02/research/thesis/node54.html
c computation of G_k term later
        if(if3d)then
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          call gradm11(w_x,w_y,w_z,vz,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2 + w_z(i)**2)
     $            +   (u_y(i) + v_x(i))**2
     $            +   (u_z(i) + w_x(i))**2
     $            +   (w_y(i) + v_z(i))**2
          enddo
        else 
          call gradm11(u_x,u_y,u_z,vx,e)
          call gradm11(v_x,v_y,v_z,vy,e)
          do i=1,lxyz
             g(i) = 2*(u_x(i)**2 + v_y(i)**2)
     $            +   (u_y(i) + v_x(i))**2
          enddo
        endif

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

c solve for omega_pert

c ---------------------
c        call check_omwall_behavior
c ---------------------
        do i=1,lxyz

          rho = param(1) ! vtrans(i,1,1,e,1)
          mu  = param(2) ! vdiff (i,1,1,e,1)
          nu  = mu/rho

c limits for k, omega

c solve for omega_pert
          omega   = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          iflim_omeg = 0 ! limit omega^{prime} 1 limit omega_total

          if    (iflim_omeg.eq.0) then
            if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c             write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
              xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
              nome_neg = nome_neg + 1
c             write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
              t(i,1,1,e,ifld_omega-1) =0.01*abs(t(i,1,1,e,ifld_omega-1))
              omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
            endif

            if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c             write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
              xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
              nkey_neg = nkey_neg + 1
c             write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
              t(i,1,1,e,ifld_k-1) = 0.01*abs(t(i,1,1,e,ifld_k-1))
              k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
            endif

          elseif(iflim_omeg.eq.1) then

            if(omega.lt.0.0) then
              write(*,*) 'OMEG tot is neg', omega
              omega = 0.01*abs(omega)
            endif

            if(k.lt.0.0) then
               write(*,*) 'K  is neg', k
               k = 0.01*abs(k)
            endif

          endif

          expn = -2.0 
          o_x(i)= omp_x(i)
          o_y(i)= omp_y(i)
          if(if3d) 
     $    o_z(i)= omp_z(i)

c See equations from eqns_k_omega1.pdf from Eq. (3) onwards
c Eq.(1) and (2) in eqns_k_omega1.pdf are the governing equations
c no source terms Sk or S_w are added

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

! calculate del k * del omega / omega

          if(if3d)then
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
     $                                            / (omega**3+tiny)
          else
            xk3 = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
     $                                            / (omega**3+tiny)
          endif

          re_t    = rho * k /(mu * omega + tiny) 
          alp_str = alpinf_str * (alp0_str + (re_t/r_k))
     $                                         / (1+(re_t/r_k)) 
          mu_t    = rho * alp_str*k/(omega + tiny) ! should multiply these by rho!!!

c nu_t is kinematic turbulent viscosity
c units of k = m2/s2, units of omega = 1/s, nu units = m2/s
c set limit for nu_t
c           nu_t    = max(tiny,nu_t)
c	    if(nu_t.gt.5000*mu)nu_t = 5000*mu

          betai_str = betainf_str * (akk + (re_t/r_b)**4)
     $              / (1.0 + (re_t/r_b)**4)


          if (xk3.le.0)then
            f_beta_str = 1.0
          else
            f_beta_str = (1.0 + 680.0*xk3*xk3)/(1 + 400.0*xk3*xk3)
          endif
          Y_k = rho * betai_str * f_beta_str * k * omega
 
c betai_str = beta_star in 12.5.15 for incompressible flow
 
          G_k = mu_t*g(i)
c g(i) is S**2 in equation sheet
          kSrc  (i,1,1,e) = G_k - Y_k

c Compute production of omega
          alpha = (alp_inf/alp_str) * 
     $          ((alpha_0 + (re_t/r_w))/(1 + (re_t/r_w)))

          G_w = alpha*alp_str*rho*g(i)

c Compute dissipation of omega
          beta = beta_0
c no compressibility correction M < 0.25

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
          if(if3d) f_b = (1.0 + 70.0*x_w)/(1.0 + 80.0*x_w)

          Y_w = rho*beta*f_b * omega * omega

          omgSrc(i,1,1,e) = G_w - Y_w

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      nome_neg =iglsum(nome_neg,1)
      nkey_neg =iglsum(nkey_neg,1)
      xome_neg = glmin(xome_neg,1)
      xkey_neg = glmin(xkey_neg,1)
c      write(*,*) 'Neg Omega ', nome_neg, xome_neg
c      write(*,*) 'Neg Key   ', nkey_neg, xkey_neg

      if(nid.eq.0) then
        if(nome_neg.gt.0) write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nkey_neg.gt.0) write(*,*) 'Neg Key   ', nkey_neg, xkey_neg
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_init(ifld_k_in,ifld_omega_in,ifcoeffs
     $                       ,coeffs_in,wall_id,ywd_in,model_id)
c
c     Initialize values ifld_omega & ifld_k for RANS k-omega turbulence
c     modeling
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      real w1,w2,w3,w4,w5
      common /SCRNS/
     & w1(lx1*ly1*lz1*lelv)
     &,w2(lx1*ly1*lz1*lelv)
     &,w3(lx1*ly1*lz1*lelv)
     &,w4(lx1*ly1*lz1*lelv)
     &,w5(lx1*ly1*lz1*lelv)

      integer n,wall_id
      real coeffs_in(1),ywd_in(1)
      logical ifcoeffs

      character*3 bcw
      character*32 mname(6)

      data mname
     &/'regularized standard k-omega    '
     &,'regularized low-Re k-omega      '
     &,'regularized standard k-omega SST'
     &,'regularized low-Re k-omega SST  '
     &,'non-regularized standard k-omega'
     &,'non-regularized low-Re komega   '/

      n=nx1*ny1*nz1*nelv

      ifrans_komg_stndrd       = .FALSE.
      ifrans_komg_lowRe        = .FALSE.
      ifrans_komgSST_stndrd    = .FALSE.
      ifrans_komgSST_lowRe     = .FALSE.
      ifrans_komg_stndrd_noreg = .FALSE.
      ifrans_komg_lowRe_noreg  = .FALSE.
      if(model_id .eq.0) ifrans_komg_stndrd       = .TRUE.
      if(model_id .eq.1) ifrans_komg_lowRe        = .TRUE.
      if(model_id .eq.2) ifrans_komgSST_stndrd    = .TRUE.
      if(model_id .eq.3) ifrans_komgSST_lowRe     = .TRUE.
      if(model_id .eq.4) ifrans_komg_stndrd_noreg = .TRUE.
      if(model_id .eq.5) ifrans_komg_lowRe_noreg  = .TRUE.

      if(nid.eq.0) write(*,'(a,1x,a)') 
     &                      'Using model:',mname(model_id+1)
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
         if(ifrans_komg_stndrd .or. ifrans_komg_lowRe) 
     $                    call rans_komg_set_defaultcoeffs
         if(ifrans_komgSST_stndrd .or. ifrans_komgSST_lowRe) 
     $                    call rans_komgSST_set_defaultcoeffs
         if(ifrans_komg_stndrd_noreg .or. ifrans_komg_lowRe_noreg)
     $                    call rans_komg_set_defaultcoeffs
      endif

c solve for omega_pert
      if(wall_id.eq.0) then
        write(*,*) 'Using user supplied wall distance'
        call copy(ywd,ywd_in,n)
      else
        bcw    = 'W  '
        ifld   = 1
        if(nid.eq.0) write(*,*) 'BC for distance ', bcw
        if(wall_id.eq.1) call cheap_dist(ywd,ifld,bcw)
        if(wall_id.eq.2) call distf(ywd,ifld,bcw,w1,w2,w3,w4,w5)
        call copy(ywd_in,ywd,n)
      endif

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

      parameter (nprof_max = 10000)
      integer i,j,k,e
      real    nu

c     real dudx(lx1,ly1,lz1,lelv), dudy(lx1,ly1,lz1,lelv)
c    $   , dudz(lx1,ly1,lz1,lelv), temt(lx1,ly1,lz1,lelv)

      real kv_min,omeg_max

      omeg_max   = coeffs(15)
      beta0      = coeffs(6)
      nu  = param(2)/param(1)
      ntot1 = nx1*ny1*nz1*nelv

      betainf_str = coeffs(11)
      yw_min = glmin(ywd,ntot1)

      Cfcon = 6.0 * nu / beta0 ! 2.0 * nu0 / betainf_str ! 
      expn  = -2.0
c     write(*,*) 'Cf, beta, ywd_min is ', Cfcon, betainf_str, yw_min

      call gradm1 (dfdx_omegb,dfdy_omegb,dfdz_omegb,   ywd)
      call opcolv (dfdx_omegb,dfdy_omegb,dfdz_omegb,   bm1)
      call opdssum(dfdx_omegb,dfdy_omegb,dfdz_omegb)
      call opcolv (dfdx_omegb,dfdy_omegb,dfdz_omegb,binvm1)

c     call gradm1 (dudx,      dudy      ,dudz,  dfdx_omegb)
c     call copy   (delsqf_omegb, dudx, ntot1)
c     call gradm1 (dudx,      dudy      ,dudz,  dfdy_omegb)
c     call add2   (delsqf_omegb, dudy, ntot1)
c     call gradm1 (dudx,      dudy      ,dudz,  dfdz_omegb)
c     call add2   (delsqf_omegb, dudz, ntot1)

      toll = 1.0e-08

c     call vdot3  (delfsq_omegb,dfdx_omegb,dfdy_omegb,dfdz_omegb 
c    $                         ,dfdx_omegb,dfdy_omegb,dfdz_omegb,ntot1)

      call rzero  (delsqf_omegb, ntot1)
      call rone   (delfsq_omegb, ntot1)

      do e = 1,nelv
      do k = 1,nz1
      do j = 1,ny1
      do i = 1,nx1

         ieg = lglel(e)
         yw   = ywd(i,j,k,e)       ! 1.0 - abs(y)
         ywmin=sqrt(Cfcon/omeg_max)
         if(yw.gt.ywmin) then
            ywm1 = 1.0 /yw
         else
c           if(yw.ne.0) write(*,'(a,3G14.7,4I5)') 
c    $                   'ywmin and yw ',ywmin,yw,omeg_max, i, j, k, ieg
            ywm1 = 1.0 /ywmin
         endif
         ywdm1(i,j,k,e) = ywm1
         ywm2 = ywm1*ywm1
         ywm3 = ywm2*ywm1
         ywm4 = ywm2*ywm2

         f_omegb     (i,j,k,e) =        Cfcon * ywm2
         delfpart              = expn * ywm1
         dfdx_omegb  (i,j,k,e) = dfdx_omegb(i,j,k,e)  * ywm1
         dfdy_omegb  (i,j,k,e) = dfdy_omegb(i,j,k,e)  * ywm1
         dfdz_omegb  (i,j,k,e) = dfdz_omegb(i,j,k,e)  * ywm1
         delsqfpart            = expn * (expn - 1.0)  * ywm2
         delsqf_omegb(i,j,k,e) =(delsqfpart  * delfsq_omegb(i,j,k,e)
     $                         + delfpart    * delsqf_omegb(i,j,k,e))
         delfsq_omegb(i,j,k,e) = delfsq_omegb(i,j,k,e)* ywm2

      enddo
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_set_defaultcoeffs
c
      include 'SIZE'
      include 'RANS_KOMG'

c ====various problem-specific turbulence constants
c omeg_max = value of omega on the walls
c kv_min = value of K on the walls
c Pr_t is the turbulent prandtl number

        vkappa = 0.41

c Turbulent viscosity constants
        Pr_t         = 0.85
        coeffs( 1)   = Pr_t
        sigma_k      = 2.0
        coeffs( 2)   = sigma_k
        sigma_omega  = 2.0
        coeffs( 3)   = sigma_omega

c Low Reynolds number correction constants

c Production of K constants
        alpinf_str   = 1.0
        coeffs( 4)   = alpinf_str
        r_k          = 6.0
        coeffs( 5)   = r_k
        beta_0       = 0.072 ! should be 0.075 for SST
        coeffs( 6)   = beta_0
        alp0_str     = beta_0/3.0
        coeffs( 7)   = alp0_str

c Dissipation of K constants
        betainf_str  = 0.09
        coeffs( 8)   = betainf_str
        alp_inf      = 0.52
c        alp_inf      = beta_0/betainf_str 
c     $               - vkappa**2/sqrt(betainf_str)/sigma_omega ! should be 0.52 for k-omega
        coeffs( 9)   = alp_inf
        r_b          = 8.0
        coeffs(10)   = r_b
        akk          = 4.0/15.0
        coeffs(11)   = akk

c Production of omega constants
        alpha_0      = 1.0/9.0
        coeffs(12)   = alpha_0
        r_w          = 2.95
        coeffs(13)   = r_w

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 2.0e8 ! 400.0 
        coeffs(15)   = omeg_max
        tiny         = 1.0e-20
        coeffs(16)   = tiny

        return
        end
c-----------------------------------------------------------------------
      subroutine rans_komgSST_set_defaultcoeffs
c
      include 'SIZE'
      include 'RANS_KOMG'

c ====various problem-specific turbulence constants
c omeg_max = value of omega on the walls
c kv_min = value of K on the walls
c Pr_t is the turbulent prandtl number

        vkappa = 0.41

c Turbulent viscosity constants
        Pr_t         = 0.85
        coeffs( 1)   = Pr_t
        sigk1        = 0.5 ! should be 0.85 for SST
        coeffs( 2)   = sigk1
        sigom1       = 0.5
        coeffs( 3)   = sigom1

c Low Reynolds number correction constants

c Production of K constants
        alpinf_str   = 1.0
        coeffs( 4)   = alpinf_str
        r_k          = 6.0
        coeffs( 5)   = r_k
        beta_0       = 0.072 ! should be 0.075 for SST
        coeffs( 6)   = beta_0
        alp0_str     = beta_0/3.0
        coeffs( 7)   = alp0_str

c Dissipation of K constants
        betainf_str  = 0.09
        coeffs( 8)   = betainf_str
        alp_inf      = beta_0/betainf_str 
     $               - sigom1 * vkappa**2/sqrt(betainf_str) ! should be 0.52 for k-omega
c        alp_inf      = 0.52
        coeffs( 9)   = alp_inf
        r_b          = 8.0
        coeffs(10)   = r_b
        akk          = 4.0/15.0
        coeffs(11)   = akk

c Production of omega constants
        alpha_0      = 1.0/9.0
        coeffs(12)   = alpha_0
        r_w          = 2.95
        coeffs(13)   = r_w

c Dissipation of omega constants
c         beta_0 = defined earlier	     

        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 2.0e8 ! 400.0 
        coeffs(15)   = omeg_max
        tiny         = 1.0e-20
        coeffs(16)   = tiny

c additional SST and k and epsilon constants
        alp1         = 0.31
        coeffs(17)   = alp1
        beta2        = 0.0828
        coeffs(18)   = beta2
        sigk2        = 1.0
        coeffs(19)   = sigk2
        sigom2       = 0.856
        coeffs(20)   = sigom2
        gamma2       = beta2/betainf_str 
     $               - sigom2 * vkappa**2/sqrt(betainf_str) ! should be 0.44 for k-epsilon
c        gamma2       = 0.44
        coeffs(21)   = gamma2

        return
        end
c-----------------------------------------------------------------------
