      real function rans_mut(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      save ifldla
      data ifldla /ldimt1/ 

      if(ix*iy*iz*iel.eq.1 .and. ifield.le.ifldla) then
         if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'updating rans_mut'
         if(ifrans_komg_stndrd)      call rans_komg_stndrd_eddy
         if(ifrans_komg_lowRe)       call rans_komg_lowRe_eddy
         if(ifrans_komgSST_stndrd)   call rans_komgSST_stndrd_eddy
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_eddy_noreg
         if(ifrans_ktau_stndrd)      call rans_ktau_stndrd_eddy
         if(ifrans_ktau_lowRe)       call rans_ktau_lowRe_eddy
         if(ifrans_ktauSST_stndrd)   call rans_ktauSST_stndrd_eddy
      endif

      ifldla = ifield
      rans_mut = mut(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_turbPrandtl
      include 'RANS_KOMG'

      rans_turbPrandtl = coeffs(1)

      return
      end
c-----------------------------------------------------------------------
      real function rans_mutsk(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_mutsk = mutsk(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_mutso(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      rans_mutso = mutso(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_getCoeffs(coeffs_in)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      real coeffs_in(1)

      do i=1,ncoeffs
         coeffs_in(i) = coeffs(i)
      enddo
      
      return
      end
c-----------------------------------------------------------------------
      real function rans_kSrc(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      logical ifevalsrc
      data ifevalsrc /.true./
      common /komgifsrc/ ifevalsrc

      if(ix*iy*iz*iel.eq.1 .and. ifevalsrc) then
         if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'updating rans_src'
         if(ifrans_komg_stndrd)      call rans_komg_stndrd_compute
         if(ifrans_komg_lowRe)       call rans_komg_lowRe_compute
         if(ifrans_komgSST_stndrd)   call rans_komgSST_stndrd_compute
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_compute_noreg
         if(ifrans_ktau_stndrd)      call rans_ktau_stndrd_compute
         if(ifrans_ktau_lowRe)       call rans_ktau_lowRe_compute
         if(ifrans_ktauSST_stndrd)   call rans_ktauSST_stndrd_compute
         ifevalsrc = .false.
      endif

      if(ifld_k.gt.ifld_omega) ifevalsrc = .true.
      rans_kSrc = kSrc(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_omgSrc(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      logical ifevalsrc
      data ifevalsrc /.true./
      common /komgifsrc/ ifevalsrc

      if(ix*iy*iz*iel.eq.1 .and. ifevalsrc) then
         if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'updating rans_src'
         if(ifrans_komg_stndrd)      call rans_komg_stndrd_compute
         if(ifrans_komg_lowRe)       call rans_komg_lowRe_compute
         if(ifrans_komgSST_stndrd)   call rans_komgSST_stndrd_compute
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_compute_noreg
         if(ifrans_ktau_stndrd)      call rans_ktau_stndrd_compute
         if(ifrans_ktau_lowRe)       call rans_ktau_lowRe_compute
         if(ifrans_ktauSST_stndrd)   call rans_ktauSST_stndrd_compute
         ifevalsrc = .false.
      endif

      if(ifld_omega.gt.ifld_k) ifevalsrc = .true.
      rans_omgSrc = omgSrc(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_kDiag(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      logical ifevalsrc
      data ifevalsrc /.true./
      common /komgifsrc/ ifevalsrc

      if(ix*iy*iz*iel.eq.1 .and. ifevalsrc) then
         if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'updating rans_src'
         if(ifrans_komg_stndrd)      call rans_komg_stndrd_compute
         if(ifrans_komg_lowRe)       call rans_komg_lowRe_compute
         if(ifrans_komgSST_stndrd)   call rans_komgSST_stndrd_compute
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_compute_noreg
         if(ifrans_ktau_stndrd)      call rans_ktau_stndrd_compute
         if(ifrans_ktau_lowRe)       call rans_ktau_lowRe_compute
         if(ifrans_ktauSST_stndrd)   call rans_ktauSST_stndrd_compute
         ifevalsrc = .false.
      endif

      if(ifld_k.gt.ifld_omega) ifevalsrc = .true.
      rans_kDiag = kDiag(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      real function rans_omgDiag(ix,iy,iz,iel)
      include 'SIZE'
      include 'TSTEP'
      include 'RANS_KOMG'

      logical ifevalsrc
      data ifevalsrc /.true./
      common /komgifsrc/ ifevalsrc

      if(ix*iy*iz*iel.eq.1 .and. ifevalsrc) then
         if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'updating rans_src'
         if(ifrans_komg_stndrd)      call rans_komg_stndrd_compute
         if(ifrans_komg_lowRe)       call rans_komg_lowRe_compute
         if(ifrans_komgSST_stndrd)   call rans_komgSST_stndrd_compute
         if(ifrans_komg_stndrd_noreg)call rans_komg_stndrd_compute_noreg
         if(ifrans_ktau_stndrd)      call rans_ktau_stndrd_compute
         if(ifrans_ktau_lowRe)       call rans_ktau_lowRe_compute
         if(ifrans_ktauSST_stndrd)   call rans_ktauSST_stndrd_compute
         ifevalsrc = .false.
      endif

      if(ifld_omega.gt.ifld_k) ifevalsrc = .true.
      rans_omgDiag = omgDiag(ix,iy,iz,iel)

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_init(ifld_k_in,ifld_omega_in,ifcoeffs
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

      integer n,wall_id,ifld_mx
      real coeffs_in(1),ywd_in(1)
      logical ifcoeffs,ifransD

      character*3 bcw
      character*36 mname(7)

      data mname
     &/'regularized standard k-omega        '
     &,'regularized low-Re k-omega          '
     &,'regularized standard k-omega SST    '
     &,'non-regularized standard k-omega    '
     &,'standard k-tau                      '
     &,'low-Re   k-tau                      '
     &,'standard k-tau SST                  '/

      n=lx1*ly1*lz1*nelv

      if(nid.eq.0) write(6,*) 'initialize RANS model'

      if(iflomach) then
        if(nid.eq.0) write(6,*)
     &         "ERROR: RANS NOT YET SUPPORTED WITH LOW MACH FORMULATION"
        call exitt
      endif

      ifrans_komg_stndrd       = .FALSE.
      ifrans_komg_lowRe        = .FALSE.
      ifrans_komgSST_stndrd    = .FALSE.
      ifrans_komg_stndrd_noreg = .FALSE.
      ifrans_ktau_stndrd       = .FALSE.
      ifrans_ktau_lowRe        = .FALSE.
      ifrans_ktauSST_stndrd    = .FALSE.
      if(model_id .eq.0) ifrans_komg_stndrd          = .TRUE.
      if(model_id .eq.1) ifrans_komg_lowRe           = .TRUE.
      if(model_id .eq.2) ifrans_komgSST_stndrd       = .TRUE.
      if(model_id .eq.3) ifrans_komg_stndrd_noreg    = .TRUE.
      if(model_id .eq.4) ifrans_ktau_stndrd          = .TRUE.
      if(model_id .eq.5) ifrans_ktau_lowRe           = .TRUE.
      if(model_id .eq.6) ifrans_ktauSST_stndrd       = .TRUE.

c split diagonal of the source term into implicit, by Sigfried
      ifrans_diag=.TRUE.

      if(nid.eq.0) write(*,'(a,a)')
     &                      '  model: ',mname(model_id+1)
      if(nio.eq.0) write(*,*)
     &                      '  ifrans_diag: ',ifrans_diag
      ifld_k     = ifld_k_in
      ifld_omega = ifld_omega_in
      ifld_mx=max(ifld_k,ifld_omega)
      if (ifld_mx.gt.ldimt1)
     $  call exitti('nflds gt ldimt+1, recompile with ldimt > ',
     $  ifld_mx+1)

! specify k-omega model coefficients


      if(ifcoeffs) then
        if(ncoeffs_in.lt.ncoeffs) call exitti(
     $   'dim of user provided komg coeffs array should be >=$',ncoeffs)
        do i=1,ncoeffs
          coeffs(i) =coeffs_in(i)
        enddo
      else
        if(ifrans_komg_stndrd .or. ifrans_komg_lowRe .or.
     $  ifrans_komg_stndrd_noreg .or. ifrans_ktau_stndrd .or.
c     $  ifrans_ktau_lowRe) call rans_komg_set_defaultcoeffs
     $  ifrans_ktau_lowRe) call rans_komg2006_set_defaultcoeffs
        if(ifrans_komgSST_stndrd .or. ifrans_ktauSST_stndrd)
     $                               call rans_komgSST_set_defaultcoeffs
      endif

c setup wall distance
      if(wall_id.eq.0) then
        if(nid.eq.0) write(6,*) ' user supplied wall distance'
        call copy(ywd,ywd_in,n)
      else
        bcw    = 'W  '
        ifld   = 1
        if(nid.eq.0) write(6,*) 'BC for distance , w_id ',bcw, wall_id
        if(wall_id.eq.1) call cheap_dist(ywd,ifld,bcw)
        if(wall_id.eq.2) call distf(ywd,ifld,bcw,w1,w2,w3,w4,w5)
        call copy(ywd_in,ywd,n)
      endif

c set cbc array for k and omega/tau (need to revise for wall-functions)
      do 10 ie = 1,nelv
      do 10 ifc = 1,2*ndim
        bcw=cbc(ifc,ie,1)
        cbc(ifc,ie,ifld_k)=bcw
        cbc(ifc,ie,ifld_omega)=bcw
        if(bcw.eq.'W  '.or.bcw.eq.'v  ') then
          cbc(ifc,ie,ifld_k)='t  '
          cbc(ifc,ie,ifld_omega)='t  '
        elseif(bcw.eq.'SYM'.or.bcw.eq.'O  '.or.bcw.eq.'o  ') then
          cbc(ifc,ie,ifld_k)='I  '
          cbc(ifc,ie,ifld_omega)='I  '
        endif
  10  continue

c solve for omega_wall & setup molecular viscosity
      call rans_komg_omegabase
      call cfill(mul,cpfld(1,1),n)

      if(nid.eq.0) write(6,*) 'done :: init RANS'

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

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,o_x(lxyz),o_y(lxyz),o_z(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)

      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg ! limit omega^{prime}
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy (div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          omp   = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k-1)     ! from previous timestep
          omw   = f_omegb(i,1,1,e)        ! omega wall
          omega = omp + omw               ! total omega

          expn = -2.0
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d)
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom(i) =  1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))
          sum_xx  =      OiOjSk (i,e)

! calculate del k * del omega / omega
          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
          endif

          alp_str   = alpinf_str 
          betai_str = betainf_str

          rhoalpk(i)= rho*alp_str*k

          factr = 1.0
          sigma_omega1 = 1.0/sigma_omega
          rhoalpfr(i) = expn * rho * alp_str * factr * omwom(i)
     $                                           * sigma_omega1

          f_beta_str = 1.0
          sigd       = sigd_min

          if (xk.gt.0)then
            xk3 = xk/(omega**3+tiny)
            f_beta_str = (1.0 + fb_c1st*xk3*xk3)/(1.0 + fb_c2st*xk3*xk3)
            sigd       = sigd_max
          endif

c calculate mu_t
          mu_t    = rho * alp_str*k/(omega + tiny)
          mu_k    = rho * alp_str  /(omega + tiny)

c Compute Y_k = dissipation of k
          Y_k = rho * betai_str * f_beta_str * omega

c Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho*k + mu_t*div(i) )*extra_prod
          G_p =             ( rho   + mu_k*div(i) )*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - ( Y_k + G_p ) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          alpha = (alp_inf/alp_str)

c         G_w = alpha*alp_str*rho*(g(i)-extra_prod*(omega+div(i)))! the full term
          G_w = alpha*alp_str*rho*(g(i)-extra_prod*(omw  +div(i)))! the explicit term
          G_p = alpha*alp_str*rho*(     extra_prod               )! the implicit term

c Compute dissipation of omega
          beta = beta_0

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
c         if(if3d) f_b = (1.0 + fb_c1*x_w)/(1.0 + fb_c2*x_w)

          Y_w1= rho*beta*f_b * omw * omw
          Y_w2= rho*beta*f_b *(2.*omw + omp)
          Y_w = rho*beta*f_b * omega * omega

c Compute extra source term of omega
          S_w = rho * sigd * xk / (omega+tiny)

          if (ifrans_diag) then
            omgSrc(i,1,1,e) = G_w - Y_w1 + S_w
            omgDiag(i,1,1,e)= Y_w2 + G_p
          else
            omgSrc(i,1,1,e) = G_w - Y_w + S_w  - G_p*omp
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

c solve for omega_pert
        expm = -expn
c add mut*delsqf
        sigma_omega1 = 1.0/sigma_omega
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigma_omega1,lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf ! This will cancel out Y_w at wall
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call col2   (tempv,mul,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
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

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,o_x(lxyz),o_y(lxyz),o_z(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)

      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg ! limit omega^{prime}
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy(div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          omp   = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
          omw   = f_omegb(i,1,1,e)        ! omega wall
          omega = omp + omw               ! total omega

          if    (iflim_omeg.eq.1) then

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
          omwom(i) =  1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))
          sum_xx  =      OiOjSk (i,e)

c calculate del k * del omega / omega
          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
          endif

          re_t    = rho * k /(mul(i,1,1,e) * omega + tiny) 
          alp_str  = alpinf_str  * (alp0_str + (re_t/r_k))
     $                                   / (1.+(re_t/r_k))

          rhoalpk(i)= rho*alp_str*k

          rtld = re_t / r_k
          funcr = (1.0 - alp0_str) / (alp0_str + rtld)/(1.0 + rtld)
          factr = (1.0 + rtld * funcr)
          sigma_omega1 = 1.0/sigma_omega
          rhoalpfr(i) = expn * rho * alp_str * factr * omwom(i)
     $                                           * sigma_omega1

          betai_str= betainf_str *   (akk + (re_t/r_b)**4)
     $                           /   (1.0 + (re_t/r_b)**4)

          f_beta_str = 1.0
          sigd       = sigd_min

          if (xk.gt.0)then
            xk3 = xk/(omega**3+tiny)
            f_beta_str = (1.0 + fb_c1st*xk3*xk3)/(1.0 + fb_c2st*xk3*xk3)
            sigd       = sigd_max
          endif

c calculate mu_t
          mu_t    = rho * alp_str*k/(omega + tiny)
          mu_k    = rho * alp_str  /(omega + tiny)

c Compute Y_k = dissipation of k
          Y_k = rho * betai_str * f_beta_str * omega

c Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho * k + mu_t * div(i) ) * extra_prod
          G_p =             ( rho     + mu_k * div(i) ) * extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - ( Y_k + G_p ) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          alpha = (alp_inf/alp_str) *
     $          ( (alpha_0 + (re_t/r_w))/(1.0 + (re_t/r_w)) )

c         G_w = alpha*alp_str*rho*(g(i)-extra_prod*(omega+div(i)))
          G_w = alpha*alp_str*rho*(g(i)-extra_prod*(  omw+div(i)))
          G_p = alpha*alp_str*rho*(     extra_prod               )

c Compute dissipation of omega
          beta = beta_0

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
c         if(if3d) f_b = (1.0 + fb_c1*x_w)/(1.0 + fb_c2*x_w)

          Y_w1= rho*beta*f_b * omw * omw
          Y_w2= rho*beta*f_b *(2.*omw + omp)
          Y_w = rho*beta*f_b * omega * omega

c Compute extra source term of omega
          S_w = rho * sigd * xk / (omega+tiny)

          if (ifrans_diag) then
            omgSrc(i,1,1,e) = G_w - Y_w1+ S_w
            omgDiag(i,1,1,e)= Y_w2 + G_p
          else
            omgSrc(i,1,1,e) = G_w - Y_w + S_w - G_p*omp
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

c solve for omega_pert
        expm = -expn
c add mut*delsqf
        sigma_omega1 = 1.0/sigma_omega
        call copy   (mu_omeg,rhoalpk,lxyz)
        call cmult  (mu_omeg,sigma_omega1,lxyz)
        call col4   (extra_src_omega,mu_omeg ,omwom
     $                                ,delsqf_omegb(1,1,1,e),lxyz)

c add mu*delsqf ! This will cancel out Y_w at wall
        call copy   (tempv,delsqf_omegb(1,1,1,e),lxyz)
        call col2   (tempv,mul,lxyz)
        call col2   (tempv,f_omegb(1,1,1,e),lxyz)
        call add2   (extra_src_omega, tempv,lxyz)

c Form (1/sigma_w) del_mut * del_omw
c  form 1: (del_yw/yw) del_k
        call col3   (term1,dfdx_omegb(1,1,1,e),k_x,lxyz)
        call addcol3(term1,dfdy_omegb(1,1,1,e),k_y,lxyz)
        if(if3d)
     $  call addcol3(term1,dfdz_omegb(1,1,1,e),k_z,lxyz)

c  form 2: 2(omw/om) (del_omw/yw)^2
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

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,o_x(lxyz),o_y(lxyz),o_z(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)

      real           tempv(lxyz), rhoalpk (lxyz)
     $              ,omwom(lxyz), rhoalpfr(lxyz)

      real           term1(lxyz), term2   (lxyz)
     $              ,term3(lxyz), term4   (lxyz)
     $              ,g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

      real mu_omeg(lxyz), mu_omegx(lxyz), mu_omegy(lxyz), mu_omegz(lxyz)
      real extra_src_omega(lxyz)

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy(div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          mu    = mul(i,1,1,e)
          nu    = mu/rho
          omp   = t(i,1,1,e,ifld_omega-1) ! Current k & omega prime values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep
          omw   = f_omegb(i,1,1,e)        ! omega wall
          omega = omp + omw               ! total omega

          expn = -2.0
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d)
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          omwom(i) =  1.0/(1.0+t(i,1,1,e,ifld_omega-1)/f_omegb(i,1,1,e))

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))

c calculate del k * del omega / omega
          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
          endif

          rhoalpk(i)= rho*k
          rhoalpfr(i) = expn * rho * omwom(i) * sigom1

c calculate F2 based on arg2
          yw     = ywd  (i,1,1,e)
          ywm1   = ywdm1(i,1,1,e)
          ywm2   = ywm1*ywm1
          arg2_1 =     sqrt(k) * ywm1 / omega / beta_str
          arg2_2 =          500.0*nu * ywm2 / omega
          arg2   = 2.0*arg2_1
          argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
          if(2.0*argF2 .le. 1.0) arg2   = arg2_2
          Fun2   = tanh(arg2 * arg2)

c calculate F1 based on arg1
          tinySST= 1.0e-10
          arg1_1 = arg2_1
          if(    argF2 .le. 1.0) arg1_1   = arg2_2
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
          argF1  = tinySST * omega /(2.0 * rho * sigom2)
          if(xk .gt. argF1) arg1_2 =   2.0 * k * omega * ywm2 / xk
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

c calculate mu_t
          mu_t   = rho * k/(omega + tiny)
          mu_k   = rho    /(omega + tiny)
          argn   = Fun2*St_magn ! this can also be Om_magn
          if(omega.le.argn/alp1) then
             mu_t   = rho * alp1 * k/argn
             mu_k   = rho * alp1    /argn
             denom  = argn/ alp1
          else
             denom  = omega
          endif

c Compute Y_k = dissipation of k
          Y_k = rho * beta_str * omega

c Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)! - ( rho*k + mu_t*div(i) )*extra_prod
          G_p =              ( rho   + mu_k*div(i) )*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k 
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - (G_p + Y_k) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          beta  = Fun1 * beta1  + (1.0 - Fun1) * beta2
          gamma = Fun1 * gamma1 + (1.0 - Fun1) * gamma2
          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

          G_w = rho * gamma * (g(i)-(div(i)+denom)*extra_prod)

c Compute dissipation of omega
          Y_w1= rho*beta* omw * omw
          Y_w2= rho*beta*(2.*omw + omp)
          Y_w = rho*beta* omega * omega

c Compute additional SST term for omega
          S_w = 2.0 * rho * sigom2 * (1.0 - Fun1) * xk / (omega+tiny)

          if (ifrans_diag) then
            omgSrc(i,1,1,e) = G_w - Y_w1+ S_w
            omgDiag(i,1,1,e)= Y_w2
          else
            omgSrc(i,1,1,e) = G_w - Y_w + S_w
            omgDiag(i,1,1,e)= 0.0
          endif

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

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,o_x(lxyz),o_y(lxyz),o_z(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

        vkappa = 0.41

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

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg_noreg
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy(div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          mu    = mul(i,1,1,e) 
          nu    = mu/rho
          omega = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k-1)     ! from previous timestep

          expn = -2.0
          o_x(i)= omp_x(i)
          o_y(i)= omp_y(i)
          if(if3d) o_z(i)= omp_z(i)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))
          sum_xx  =      OiOjSk (i,e)

c calculate del k * del omega / omega
          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
          endif

          alp_str   = alpinf_str 
          betai_str = betainf_str

          f_beta_str = 1.0
          sigd       = sigd_min

          if (xk.gt.0)then
            xk3 = xk/(omega**3+tiny)
            f_beta_str = (1.0 + fb_c1st*xk3*xk3)/(1.0 + fb_c2st*xk3*xk3)
            sigd       = sigd_max
          endif

c calculate mu_t
          mu_t    = rho * alp_str*k/(omega + tiny)
          mu_k    = rho * alp_str  /(omega + tiny)

          yw   = ywd  (i,1,1,e)
          toll = 1.0e-08
          if(yw.le.toll) then
            veddy = vkappa*nu*yplus
            mu_t  = rho * veddy
          endif

c Compute Y_k = dissipation of k
          Y_k = rho * betai_str * f_beta_str * omega

c Compute G_k = production of k
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho*k + mu_t*div(i) )*extra_prod
          G_p = (rho + mu_k*div(i))*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k
          else
            kSrc  (i,1,1,e) = G_k - Y_k * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          alpha = (alp_inf/alp_str)

c         G_w = alpha*alp_str*rho*(g(i)-(omega+div(i))*extra_prod)
          G_w = alpha*alp_str*rho*(g(i)-(      div(i))*extra_prod)
          G_p = alpha*alp_str*rho*(                    extra_prod)

c Compute dissipation of omega
          beta = beta_0

          x_w = abs((sum_xx)/(betainf_str*omega + tiny)**3)
          f_b = 1.0
c         if(if3d) f_b = (1.0 + fb_c1*x_w)/(1.0 + fb_c2*x_w)

c         Y_w = rho*beta*f_b * omega * omega
          Y_wp= rho*beta*f_b * omega

c Compute extra source term of omega
          S_w = rho * sigd * xk / (omega+tiny)

          if (ifrans_diag) then
            omgSrc(i,1,1,e) = G_w + S_w
            omgDiag(i,1,1,e)= Y_wp + G_p
          else
            omgSrc(i,1,1,e) = G_w + S_w - (Y_wp + G_p) * omega
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktau_stndrd_compute
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,t_x(lxyz),t_y(lxyz),t_z(lxyz)
     $              ,tau_x(lxyz), tau_y(lxyz), tau_z(lxyz)

      real           tausq(lxyz,lelv)
     $              ,tsq_x(lxyz), tsq_y(lxyz), tsq_z(lxyz)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv

      mu_min    = edd_frac_free*param(2)

      call limit_ktau ! check for negative values
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)
      call sqrt_tau(tausq,t(1,1,1,1,ifld_omega-1),ntot)
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy(div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(tau_x,tau_y,tau_z,t(1,1,1,1,ifld_omega-1),e)
        call gradm11(tsq_x,tsq_y,tsq_z,tausq                  ,e)

        do i=1,lxyz

          rho = vtrans(i,1,1,e,1)
          mu  = mul(i,1,1,e)
          tau = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k   = t(i,1,1,e,ifld_k-1)     ! from previous timestep

          t_x(i)= tau_x(i)
          t_y(i)= tau_y(i)
          if(if3d) t_z(i)= tau_z(i)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))
          sum_xx  =      OiOjSk (i,e)

! calculate del k * del tau   / tau  
          if(if3d)then
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i) + k_z(i)*t_z(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i) + t_z(i)*t_z(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i)+tsq_z(i)*tsq_z(i))
          else
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i))
          endif

          alp_str   = alpinf_str 
          betai_str = betainf_str

          f_beta_str = 1.0
          sigd       = sigd_min

          if (xk.gt.0)then
            xk3 = xk * tau
            f_beta_str = (1.0 + fb_c1st*xk3*xk3)/(1.0 + fb_c2st*xk3*xk3)
            sigd       = sigd_max
          endif

c calculate mu_t
          mu_t = rho * alp_str * k * tau    ! eddy viscosity
          mu_k = rho * alp_str *     tau    ! eddy viscosity without k
          mu_tp= rho * alp_str * k          ! eddy viscosity without tau

c Limit source terms in far field
          yw   = ywd  (i,1,1,e)
          Rfact= 1.
          if( mu_t.lt.mu_min .and .yw.gt.ywlim) Rfact= mu_t/mu_min

c Compute Y_k = dissipation of k
c         Y_k = 0.
c         if(tau.gt.tiny) Y_k = rho * betai_str * f_beta_str / tau
          Y_k = rho * betai_str * f_beta_str / (tau + tiny)

c Compute G_k = production of k 
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho*k + mu_t*div(i) )*extra_prod
          G_p =             ( rho   + mu_k*div(i) )*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - ( Y_k + G_p ) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of tau
          alpha = (alp_inf/alp_str)
          gamm  = alpha*alp_str

          G_p = rho*gamm*Rfact*(tau*g(i)-(1.+div(i)*tau)*extra_prod)

c Compute dissipation of tau
          beta = beta_0

          x_w = abs((sum_xx)*(tau/betainf_str)**3)
          f_b = 1.0
c         if(if3d) f_b = (1.0 + fb_c1*x_w)/(1.0 + fb_c2*x_w)

          Y_w = rho*beta*f_b * Rfact

          S_tau = 8.0*mu    *xtq * Rfact
          S_taup= 8.0*mu_tp *xtq * Rfact/sigma_omega

c Compute extra source term of tau
c         S_w =-rho * sigd * xk * tau * Rfact
          S_wp= rho * sigd * xk *       Rfact

c Compute Source term for tau
          if(ifrans_diag) then
            if(tau.le.tiny) then
              omgSrc(i,1,1,e) = Y_w - S_tau
              omgDiag(i,1,1,e)= G_p + S_taup + S_wp
            else
              omgSrc(i,1,1,e) = Y_w
              omgDiag(i,1,1,e)= G_p + S_taup + S_wp + S_tau/tau
            endif
          else
            omgSrc(i,1,1,e) = Y_w - S_tau - (G_p + S_taup + S_wp) * tau
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktau_lowRe_compute
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,t_x(lxyz),t_y(lxyz),t_z(lxyz)
     $              ,tau_x(lxyz), tau_y(lxyz), tau_z(lxyz)

      real           tausq(lxyz,lelv)
     $              ,tsq_x(lxyz), tsq_y(lxyz), tsq_z(lxyz)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv

      mu_min    = edd_frac_free*param(2)

      call limit_ktau ! check for negative values
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)
      call sqrt_tau(tausq,t(1,1,1,1,ifld_omega-1),ntot)
      call rzero(div,lxyz)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        if(iflomach) call copy(div,DivQ(1,e),lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(tau_x,tau_y,tau_z,t(1,1,1,1,ifld_omega-1),e)
        call gradm11(tsq_x,tsq_y,tsq_z,tausq                  ,e)

        do i=1,lxyz

          rho = vtrans(i,1,1,e,1)
          mu  = mul(i,1,1,e)
          tau = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k   = t(i,1,1,e,ifld_k-1)     ! from previous timestep

          t_x(i)= tau_x(i)
          t_y(i)= tau_y(i)
          if(if3d) t_z(i)= tau_z(i)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))
          sum_xx  =      OiOjSk (i,e)

c calculate del k * del tau   / tau  
          if(if3d)then
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i) + k_z(i)*t_z(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i) + t_z(i)*t_z(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i)+tsq_z(i)*tsq_z(i))
          else
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i))
          endif

          re_t     = rho * k * tau / mu

          alp_str  = alpinf_str  * (alp0_str + (re_t/r_k))
     $                                   / (1.+(re_t/r_k))
  
          betai_str= betainf_str *   (akk + (re_t/r_b)**4)
     $                           /   (1.0 + (re_t/r_b)**4)

          f_beta_str = 1.0
          sigd       = sigd_min

          if (xk.gt.0)then
            xk3 = xk * tau
            f_beta_str = (1.0 + fb_c1st*xk3*xk3)/(1.0 + fb_c2st*xk3*xk3)
            sigd       = sigd_max
          endif

c calculate mu_t
          mu_t = rho * alp_str * k * tau    ! eddy viscosity
          mu_k = rho * alp_str *     tau    ! eddy viscosity without k
          mu_tp= rho * alp_str * k          ! eddy viscosity without tau

c Limit source terms in far field
          yw   = ywd  (i,1,1,e)
          Rfact= 1.
          if( mu_t.lt.mu_min .and .yw.gt.ywlim) Rfact= mu_t/mu_min

c Compute Y_k = dissipation of k
          Y_k = 0.
          if(tau.gt.0) Y_k = rho * betai_str * f_beta_str / tau

c Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho*k + mu_t*div(i) )*extra_prod
          G_p =             ( rho   + mu_k*div(i) )*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - ( Y_k + G_p ) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          alpha = (alp_inf/alp_str) *
     $          ( (alpha_0 + (re_t/r_w))/(1.0 + (re_t/r_w)) )

          gamm  = alpha*alp_str

          G_p = rho*gamm*Rfact*(tau*g(i)-(1.+div(i)*tau)*extra_prod)

c Compute dissipation of tau
          beta = beta_0

          x_w = abs((sum_xx)*(tau/betainf_str)**3)
          f_b = 1.0
c         if(if3d) f_b = (1.0 + fb_c1*x_w)/(1.0 + fb_c2*x_w)

          Y_w = rho*beta*f_b * Rfact

c Compute extra source term of tau
c         S_w =-rho * sigd * xk * tau * Rfact
          S_wp= rho * sigd * xk *       Rfact

c Compute Source term for tau
          S_tau = 8.0*mu    *xtq * Rfact
          S_taup= 8.0*mu_tp *xtq * Rfact/sigma_omega

          if(ifrans_diag) then
            if(tau.le.tiny) then
              omgSrc(i,1,1,e) = Y_w - S_tau
              omgDiag(i,1,1,e)= G_p + S_taup + S_wp
            else
              omgSrc(i,1,1,e) = Y_w
              omgDiag(i,1,1,e)= G_p + S_taup + S_wp + S_tau/tau
            endif
          else
            omgSrc(i,1,1,e) = Y_w - S_tau - (G_p + S_taup + S_wp) * tau
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktauSST_stndrd_compute
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,t_x(lxyz),t_y(lxyz),t_z(lxyz)
     $              ,tau_x(lxyz), tau_y(lxyz), tau_z(lxyz)

      real           tausq(lxyz,lelv)
     $              ,tsq_x(lxyz), tsq_y(lxyz), tsq_z(lxyz)

      real           g    (lxyz), div     (lxyz)

      real     tempR1(lx1,ly1,lz1,lelv)
      real     tempR2(lx1,ly1,lz1,lelv)
      real     tempR3(lx1,ly1,lz1,lelv)
      real     tempR4(lx1,ly1,lz1,lelv)
      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      ntot = nx1*ny1*nz1*nelv

      mu_min    = edd_frac_free*param(2)

      call limit_ktau ! check for negative values
      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)
      call sqrt_tau(tausq,t(1,1,1,1,ifld_omega-1),ntot)

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)
        call copy   (div, DivQ   (1,e),       lxyz)
        if(.not.iflomach) call rzero  (div,   lxyz)

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(tau_x,tau_y,tau_z,t(1,1,1,1,ifld_omega-1),e)
        call gradm11(tsq_x,tsq_y,tsq_z,tausq                  ,e)

        do i=1,lxyz

          rho = vtrans(i,1,1,e,1)
          mu  = mul(i,1,1,e)
          nu  = mu/rho
          tau = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k   = t(i,1,1,e,ifld_k-1)     ! from previous timestep

          t_x(i)= tau_x(i)
          t_y(i)= tau_y(i)
          if(if3d) t_z(i)= tau_z(i)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))

c calculate del k * del tau   / tau  
          if(if3d)then
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i) + k_z(i)*t_z(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i) + t_z(i)*t_z(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i)+tsq_z(i)*tsq_z(i))
          else
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i))
            xtq= (tsq_x(i)*tsq_x(i)+tsq_y(i)*tsq_y(i))
          endif

c calculate F2 based on arg2
          yw     = ywd  (i,1,1,e)
          ywm1   = ywdm1(i,1,1,e)
          ywm2   = ywm1*ywm1
          arg2_1 =     sqrt(k) * ywm1 * tau / beta_str
          arg2_2 =    500.0*nu * ywm2 * tau
          arg2   = 2.0*arg2_1
          argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
          if(2.0*argF2 .le. 1.0) arg2   = arg2_2
          Fun2   = tanh(arg2 * arg2)

c calculate F1 based on arg1
          tinySST= 1.0e-10
          arg1_1 = arg2_1
          if(    argF2 .le. 1.0) arg1_1   = arg2_2
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
          argF1  = tinySST * tau /(2.0 * rho * sigom2)
          if(xk .gt. argF1) arg1_2 =   2.0 * k * tau * ywm2 / xk
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

c calculate mu_t
          mu_t   = rho * k * tau
          mu_k   = rho *     tau
          argn   = Fun2*St_magn ! this can also be Om_magn
          if(alp1.le.(argn*tau)) then
             mu_t   = 0.0
             mu_k   = 0.0
             if(argn.ne.0.) then
               mu_t   = rho * alp1 * k/argn
               mu_k   = rho * alp1    /argn
             endif
             denom  = argn/ alp1
          else
             denom  = 0.
             if(tau.ne.0.) denom  = 1.0/tau
          endif

          yw   = ywd  (i,1,1,e)
          Rfact= 1.
          if( mu_t.lt.mu_min .and .yw.gt.ywlim) Rfact= mu_t/mu_min  ! limit source terms in far field

c Compute Y_k = dissipation of k
          Y_k = 0.
          if(tau.gt.0) Y_k = rho * beta_str / tau

c Compute G_k = production of  k and limit it to 10*Y_k (the dissipation of k)
          extra_prod = twothird*div(i)

          G_k = mu_t*g(i)!- ( rho*k + mu_t*div(i) )*extra_prod
          G_p =             ( rho   + mu_k*div(i) )*extra_prod

c Compute Source term for k
          if (ifrans_diag) then
            kSrc  (i,1,1,e) = G_k
            kDiag (i,1,1,e) = Y_k + G_p
          else
            kSrc  (i,1,1,e) = G_k - ( Y_k + G_p ) * k
            kDiag (i,1,1,e) = 0.0
          endif

c Compute production of omega
          beta  = Fun1 * beta1  + (1.0 - Fun1) * beta2
          gamma = Fun1 * gamma1 + (1.0 - Fun1) * gamma2
          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

          G_w  = rho*tau *gamma*(g(i)-(div(i)+denom)*extra_prod)
     $          *Rfact

c Compute dissipation of omega
          Y_w = rho * beta * Rfact

c Compute additional SST term for tau
c         S_w =-2.0 * rho * sigom2 * (1.0 - Fun1) * xk * tau * Rfact
          S_wp= 2.0 * rho * sigom2 * (1.0 - Fun1) * xk *       Rfact

c Compute Source term for omega
          S_tau = 8.0*mu   *xtq * Rfact
          S_taup= 8.0*rho*k*xtq * Rfact*sigom

          if(ifrans_diag) then
            if(tau.le.tiny) then
              omgSrc(i,1,1,e) = Y_w - S_tau
              omgDiag(i,1,1,e)= G_w + S_taup + S_wp
            else
              omgSrc(i,1,1,e) = Y_w
              omgDiag(i,1,1,e)= G_w + S_taup + S_wp + S_tau/tau
            endif
          else
            omgSrc(i,1,1,e) = Y_w - S_tau - (G_w + S_taup + S_wp) * tau
            omgDiag(i,1,1,e)= 0.0
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t * sigk
          mutso(i,1,1,e)   = mu_t * sigom

        enddo

      enddo

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

      integer ix,iy,iz,e ! had to change index k to iz to avoid conflict with real k (TKE) in RANS_KOMG

c     real dudx(lx1,ly1,lz1,lelv), dudy(lx1,ly1,lz1,lelv)
c    $   , dudz(lx1,ly1,lz1,lelv), temt(lx1,ly1,lz1,lelv)

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

      do e  = 1,nelv
      do iz = 1,lz1
      do iy = 1,ly1
      do ix = 1,lx1

         ieg = lglel(e)
         yw   = ywd(ix,iy,iz,e)       ! 1.0 - abs(y)
         ywmin=sqrt(Cfcon/omeg_max)
         if(yw.gt.ywmin) then
            ywm1 = 1.0 /yw
         else
           if(yw.ne.0) write(*,'(a,3G14.7,4I5)') 
     $                   'ywmin and yw ',ywmin,yw,omeg_max,ix,iy,iz,ieg
            ywm1 = 1.0 /ywmin
         endif
         ywdm1(ix,iy,iz,e) = ywm1
         ywm2 = ywm1*ywm1
         ywm3 = ywm2*ywm1
         ywm4 = ywm2*ywm2

         f_omegb     (ix,iy,iz,e) =        Cfcon * ywm2
         delfpart              = expn * ywm1
         dfdx_omegb  (ix,iy,iz,e) = dfdx_omegb(ix,iy,iz,e)  * ywm1
         dfdy_omegb  (ix,iy,iz,e) = dfdy_omegb(ix,iy,iz,e)  * ywm1
         dfdz_omegb  (ix,iy,iz,e) = dfdz_omegb(ix,iy,iz,e)  * ywm1
         delsqfpart            = expn * (expn - 1.0)  * ywm2
         delsqf_omegb(ix,iy,iz,e)=(delsqfpart*delfsq_omegb(ix,iy,iz,e)
     $                          + delfpart   *delsqf_omegb(ix,iy,iz,e))
         delfsq_omegb(ix,iy,iz,e) = delfsq_omegb(ix,iy,iz,e)* ywm2

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
        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 2.0e10 ! 400.0 Lan
        coeffs(15)   = omeg_max
        tiny         = 1.e-8
        coeffs(16)   = tiny

c additional constants 
        fb_c1        = 70.0
        coeffs(17)   = fb_c1
        fb_c2        = 80.0
        coeffs(18)   = fb_c2
        fb_c1st      = 680.0
        coeffs(19)   = fb_c1st
        fb_c2st      = 400.0
        coeffs(20)   = fb_c2st

        sigd_min     = 0.0
        coeffs(21)   = sigd_min
        sigd_max     = 0.0
        coeffs(22)   = sigd_max
        Clim         = 7.0/8.0
        coeffs(23)   = Clim      

c constants related to limiting source terms or mu_t
        Hlen         = 1.0
        coeffs(24)   = Hlen
        ywlim        = 0.5
        coeffs(25)   = ywlim
        edd_frac_free= 0.01
        coeffs(26)   = edd_frac_free
        tke_frac_free= 1.e-6
        coeffs(27)   = tke_frac_free

c yplus boundary related to wall functions
        yplus        = 100.
        coeffs(28)   = yplus

        if(nid.eq.0) write(*,*) 'Using kw98 coeffs'

        return
        end
c-----------------------------------------------------------------------
      subroutine rans_komg2006_set_defaultcoeffs
c
      include 'SIZE'
      include 'RANS_KOMG'

c ====various problem-specific turbulence constants
c omeg_max = value of omega on the walls
c kv_min = value of K on the walls
c Pr_t is the turbulent prandtl number

        logical if_cfl3d
        if_cfl3d = .false.

        vkappa = 0.41

c Turbulent viscosity constants
        Pr_t         = 0.85
        coeffs( 1)   = Pr_t
        sigma_k      = 1.0/0.6
        coeffs( 2)   = sigma_k
        sigma_omega  = 2.0
        coeffs( 3)   = sigma_omega

c Low Reynolds number correction constants

c Production of K constants
        alpinf_str   = 1.0
        coeffs( 4)   = alpinf_str
        r_k          = 6.0
        coeffs( 5)   = r_k
        beta_0       = 0.0708 ! should be 0.075 for SST
        if(if_cfl3d) 
     $  beta_0       = 0.075
        coeffs( 6)   = beta_0
        alp0_str     = beta_0/3.0
        coeffs( 7)   = alp0_str

c Dissipation of K constants
        betainf_str  = 0.09
        coeffs( 8)   = betainf_str
        alp_inf      = 0.52
        if(if_cfl3d) 
     $  alp_inf      = beta_0/betainf_str 
     $               - vkappa**2/sqrt(betainf_str)/sigma_omega ! should be 0.52 for k-omega
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
        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 2.0e10 ! 400.0 Lan
        coeffs(15)   = omeg_max
        tiny         = 1.e-8
        coeffs(16)   = tiny

c additional constants 
        fb_c1        = 85.0
        coeffs(17)   = fb_c1
        fb_c2        = 100.0
        coeffs(18)   = fb_c2
        fb_c1st      = 400.0
        coeffs(19)   = fb_c1st
        fb_c2st      = 400.0
        coeffs(20)   = fb_c2st

        sigd_min     = 0.0
        coeffs(21)   = sigd_min
        sigd_max     = 1.0/8.0
        coeffs(22)   = sigd_max
        Clim         = 7.0/8.0
        coeffs(23)   = Clim     

c constants related to limiting source terms or mu_t
        Hlen         = 1.0
        coeffs(24)   = Hlen
        ywlim        = 0.5
        coeffs(25)   = ywlim
        edd_frac_free= 0.01
        coeffs(26)   = edd_frac_free
        tke_frac_free= 1.e-6
        coeffs(27)   = tke_frac_free
        
c yplus boundary related to wall functions
        yplus        = 100.
        coeffs(28)   = yplus

        if(nid.eq.0) then
          if(if_cfl3d) then
            write(*,*) 'Using kw06_cfl3d coeffs'
          else
            write(*,*) 'Using kw06 coeffs'
          endif
          write(*,*) 'beta_0, alp_inf ', beta_0, alp_inf
        endif

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
        kv_min       = 0.0
        coeffs(14)   = kv_min
        omeg_max     = 2.0e8 ! 400.0 
        coeffs(15)   = omeg_max
        tiny         = 1.e-8
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
        coeff22      = 0.0
        coeffs(22)   = coeff22
        coeff23      = 0.0
        coeffs(23)   = coeff23

c constants related to limiting source terms or mu_t
        Hlen         = 1.0
        coeffs(24)   = Hlen
        ywlim        = 0.5
        coeffs(25)   = ywlim
        edd_frac_free= 0.01
        coeffs(26)   = edd_frac_free
        tke_frac_free= 1.e-6
        coeffs(27)   = tke_frac_free
        
c yplus boundary related to wall functions
        yplus        = 100.
        coeffs(28)   = yplus

        return
        end
c-----------------------------------------------------------------------
      subroutine rans_komg_stndrd_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg ! limit omega^{prime}

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

c calculate mu_t
          alp_str   = alpinf_str 

          mu_t    = rho * alp_str*k/(omega + tiny)

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_lowRe_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg ! limit omega^{prime}

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

c calculate mu_t
          re_t    = rho * k /(mul(i,1,1,e) * omega + tiny) 
          alp_str  = alpinf_str  * (alp0_str + (re_t/r_k))
     $                                   / (1.+(re_t/r_k))

          mu_t    = rho * alp_str*k/(omega + tiny)

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komgSST_stndrd_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,o_x(lxyz),o_y(lxyz),o_z(lxyz)
     $              ,omp_x(lxyz), omp_y(lxyz), omp_z(lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg ! limit omega^{prime}

      do e=1,nelv

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(omp_x,omp_y,omp_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          mu    = mul(i,1,1,e)
          nu    = mu/rho
          omega = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          expn = -2.0
          o_x(i)= omp_x(i)+expn * dfdx_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          o_y(i)= omp_y(i)+expn * dfdy_omegb(i,1,1,e) *f_omegb(i,1,1,e)
          if(if3d)
     $    o_z(i)= omp_z(i)+expn * dfdz_omegb(i,1,1,e) *f_omegb(i,1,1,e)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))

c calculate del k * del omega / omega
          if(if3d)then
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i) + k_z(i)*o_z(i))
          else
            xk = (k_x(i)*o_x(i) + k_y(i)*o_y(i))
          endif

c calculate F2 based on arg2
          yw     = ywd  (i,1,1,e)
          ywm1   = ywdm1(i,1,1,e)
          ywm2   = ywm1*ywm1
          arg2_1 =     sqrt(k) * ywm1 / omega / beta_str
          arg2_2 =          500.0*nu * ywm2 / omega
          arg2   = 2.0*arg2_1
          argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
          if(2.0*argF2 .le. 1.0) arg2   = arg2_2
          Fun2   = tanh(arg2 * arg2)

c calculate F1 based on arg1
          tinySST= 1.0e-10
          arg1_1 = arg2_1
          if(    argF2 .le. 1.0) arg1_1   = arg2_2
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
          argF1  = tinySST * omega /(2.0 * rho * sigom2)
          if(xk .gt. argF1) arg1_2 =   2.0 * k * omega * ywm2 / xk
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

c calculate mu_t
          mu_t   = rho * k/(omega + tiny)
          argn   = Fun2*St_magn ! this can also be Om_magn
          if(omega.le.argn/alp1) then
             mu_t   = rho * alp1 * k/argn
             denom  = argn/ alp1
          else
             denom  = omega
          endif

          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t * sigk
          mutso(i,1,1,e)   = mu_t * sigom

        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_komg_stndrd_eddy_noreg
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

      real kwall, kappa, kwallo

        vkappa = 0.41

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_komg_noreg ! limit omega^{prime}

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)

        do i=1,lxyz

          rho   = vtrans(i,1,1,e,1)
          mu    = mul(i,1,1,e)
          nu    = mu/rho
          omega = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k     = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

c calculate mu_t
          alp_str = alpinf_str 

          mu_t    = rho * alp_str*k/(omega + tiny)

          yw   = ywd  (i,1,1,e)
          toll = 1.0e-08
          if(yw.le.toll) then
            veddy = vkappa*nu*yplus
            mu_t  = rho * veddy
          endif

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktau_stndrd_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_ktau !check for negative values

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)

        do i=1,lxyz

          rho     = vtrans(i,1,1,e,1)
          tau     = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

c calculate mu_t
          alp_str   = alpinf_str 
          mu_t = rho * alp_str * k * tau    ! eddy viscosity

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktau_lowRe_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           g    (lxyz), div     (lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional constants 
        fb_c1        = coeffs(17)
        fb_c2        = coeffs(18)
        fb_c1st      = coeffs(19)
        fb_c2st      = coeffs(20)

        sigd_min     = coeffs(21)
        sigd_max     = coeffs(22)
        Clim         = coeffs(23)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_ktau !check for negative values

      do e=1,nelv

        call copy   (g,   St_mag2(1,e),       lxyz)
c       call copy   (g,   Om_mag2(1,e),       lxyz)

        do i=1,lxyz

          rho = vtrans(i,1,1,e,1)
          mu  = mul(i,1,1,e)
          tau     = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

c calculate mu_t
          re_t     = rho * k * tau / mu
          alp_str  = alpinf_str  * (alp0_str + (re_t/r_k))
     $                                   / (1.+(re_t/r_k))

          mu_t = rho * alp_str * k * tau    ! eddy viscosity

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t / sigma_k
          mutso(i,1,1,e)   = mu_t / sigma_omega
        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rans_ktauSST_stndrd_eddy
c
c     Compute RANS source terms and diffusivities on an 
c     element-by-element basis
c
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      parameter (lxyz=lx1*ly1*lz1)

      real           k_x(lxyz),k_y(lxyz),k_z(lxyz)
     $              ,t_x(lxyz),t_y(lxyz),t_z(lxyz)
     $              ,tau_x(lxyz), tau_y(lxyz), tau_z(lxyz)

      common /storesom/ St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      integer e

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
        kv_min       = coeffs(14)
        omeg_max     = coeffs(15)
        tiny         = coeffs(16)

c additional SST and k and epsilon constants
        alp1         = coeffs(17)
        beta2        = coeffs(18)
        sigk2        = coeffs(19)
        sigom2       = coeffs(20)
        gamma2       = coeffs(21)

c constants related to limiting source terms or mu_t
        Hlen         = coeffs(24)
        ywlim        = coeffs(25)
        edd_frac_free= coeffs(26)
        tke_frac_free= coeffs(27)

c yplus boundary related to wall functions
        yplus        = coeffs(28)
c================================

      call comp_StOm (St_mag2, Om_mag2, OiOjSk, DivQ)

      mu_min    = edd_frac_free*param(2)

      call limit_ktau !check for negative values

      do e=1,nelv

        call gradm11(k_x,  k_y,  k_z,  t(1,1,1,1,ifld_k    -1),e)
        call gradm11(tau_x,tau_y,tau_z,t(1,1,1,1,ifld_omega-1),e)

        do i=1,lxyz

          rho = vtrans(i,1,1,e,1)
          mu  = mul(i,1,1,e)
          nu  = mu/rho

          tau     = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          t_x(i)= tau_x(i)
          t_y(i)= tau_y(i)
          if(if3d) t_z(i)= tau_z(i)

          St_magn = sqrt(St_mag2(i,e))
          Om_magn = sqrt(Om_mag2(i,e))

c calculate del k * del tau   / tau  
          if(if3d)then
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i) + k_z(i)*t_z(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i) + t_z(i)*t_z(i))
          else
            xk =-(k_x(i)*t_x(i) + k_y(i)*t_y(i))
            xt = (t_x(i)*t_x(i) + t_y(i)*t_y(i))
          endif

c calculate F2 based on arg2
          yw     = ywd  (i,1,1,e)
          ywm1   = ywdm1(i,1,1,e)
          ywm2   = ywm1*ywm1
          arg2_1 =     sqrt(k) * ywm1 * tau / beta_str
          arg2_2 =    500.0*nu * ywm2 * tau
          arg2   = 2.0*arg2_1
          argF2  =     sqrt(k)*yw/(500.0*nu*beta_str)
          if(2.0*argF2 .le. 1.0) arg2   = arg2_2
          Fun2   = tanh(arg2 * arg2)

c calculate F1 based on arg1
          tinySST= 1.0e-10
          arg1_1 = arg2_1
          if(    argF2 .le. 1.0) arg1_1   = arg2_2
          arg1_2 =     4.0 * rho * sigom2 * k * ywm2 / tinySST
          argF1  = tinySST * tau /(2.0 * rho * sigom2)
          if(xk .gt. argF1) arg1_2 =   2.0 * k * tau * ywm2 / xk
          arg1   = min(    arg1_1, arg1_2)
          Fun1   = tanh(arg1 * arg1 * arg1 * arg1)

c calculate mu_t
          mu_t   = rho * k * tau
          argn   = Fun2*St_magn ! this can also be Om_magn
          if(alp1.le.(argn*tau)) then
             mu_t   = 0.0
             if(argn.ne.0.) mu_t   = rho * alp1 * k/argn
             denom  = argn/ alp1
          else
             denom  = 0.
             if(tau.ne.0.) denom  = 1.0/tau
          endif

          sigk  = Fun1 * sigk1  + (1.0 - Fun1) * sigk2
          sigom = Fun1 * sigom1 + (1.0 - Fun1) * sigom2

          mut  (i,1,1,e)   = mu_t
          mutsk(i,1,1,e)   = mu_t * sigk
          mutso(i,1,1,e)   = mu_t * sigom

        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_StOm(St_mag2,Om_mag2,OiOjSk,DivQ)
c
c     Compute the square of the magnitude of the stress and rotation tensors
c     St_mag2=2Sij*Sij=S'ij*S'ij/2, S'ij=2Sij, Sij=(dui/dxj+duj/dxi)/2
c     Om_mag2=2Oij*Oij=O'ij*O'ij/2, O'ij=OSij, Oij=(dui/dxj-duj/dxi)/2
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e

      parameter(lt  =lx1*ly1*lz1*lelv)
      parameter(lxyz=lx1*ly1*lz1)

      common /scruz/    sij  (lx1*ly1*lz1,6,lelv)
     $                , oij  (lx1*ly1*lz1,lelv,3)

      real              St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)
     $                , DivQ   (lx1*ly1*lz1,lelv)

      logical iflmc, ifdss
c
      thqrt    = 0.75

      iflmc = .false.! to add lowmach correction to 2Sij (2Sij - 2Q/3*delta_ij) 
      ifdss = .true. ! to dssum sij and oij
      nij = 3
      if (if3d.or.ifaxis) nij=6

      if(nid.eq.0 .and. loglevel.gt.2) write(*,*) 'updating StOm '
      call comp_sij_oij     (sij, oij, nij, vx, vy, vz, iflmc, ifdss) ! S'ij=2Sij, O'ij=2Oij
      call comp_sij_oij_mag2(sij, oij, nij, St_mag2, Om_mag2, OiOjSk) ! St_mag2=2S^2=S'^2/2
c                                                                     ! Om_mag2=2O^2=O'^2/2
      do e = 1, nelv
         call add3   (DivQ(1,e), sij(1,1,e), sij(1,2,e), lxyz)
         if(if3d.or.ifaxis)
     $   call add2   (DivQ(1,e), sij(1,3,e),             lxyz)
         if(iflmc)     call cmult(DivQ(1,e), thqrt,      lxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_sij_oij(sij,oij,nij,u,v,w,iflmc,ifdss)
c                                       du_i       du_j
c     Compute the stress tensor S_ij := ----   +   ----
c                                       du_j       du_i
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e
      logical iflmc, ifdss
c
      parameter(lt  =lx1*ly1*lz1*lelv)
      parameter(lxyz=lx1*ly1*lz1)

      real sij(lx1*ly1*lz1,nij,lelv), oij(lx1*ly1*lz1,lelv,3)
      real u  (lx1*ly1*lz1,lelv)
      real v  (lx1*ly1*lz1,lelv)
      real w  (lx1*ly1*lz1,lelv)

      real ur(lxyz),us(lxyz),ut(lxyz)
     $    ,vr(lxyz),vs(lxyz),vt(lxyz)
     $    ,wr(lxyz),ws(lxyz),wt(lxyz)

      real j ! Inverse Jacobian

      n    = lx1-1      ! Polynomial degree
      nxyz = lx1*ly1*lz1

      if (if3d) then     ! 3D CASE
       do e=1,nelv
        call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
        call local_grad3(vr,vs,vt,v,N,e,dxm1,dxtm1)
        call local_grad3(wr,ws,wt,w,N,e,dxm1,dxtm1)

        do i=1,nxyz

         j = jacmi(i,e)

         sij(i,1,e) = j*  ! du/dx + du/dx
     $   2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e))

         sij(i,2,e) = j*  ! dv/dy + dv/dy
     $   2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e)+vt(i)*tym1(i,1,1,e))

         sij(i,3,e) = j*  ! dw/dz + dw/dz
     $   2*(wr(i)*rzm1(i,1,1,e)+ws(i)*szm1(i,1,1,e)+wt(i)*tzm1(i,1,1,e))

         sij(i,4,e) = j*  ! du/dy + dv/dx
     $   (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e) +
     $    vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)+vt(i)*txm1(i,1,1,e) )
         oij(i,e,3) = j*  ! dv/dx - du/dy
     $   (vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)+vt(i)*txm1(i,1,1,e) 
     $   -ur(i)*rym1(i,1,1,e)-us(i)*sym1(i,1,1,e)-ut(i)*tym1(i,1,1,e) )

         sij(i,5,e) = j*  ! dv/dz + dw/dy
     $   (wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e)+wt(i)*tym1(i,1,1,e) +
     $    vr(i)*rzm1(i,1,1,e)+vs(i)*szm1(i,1,1,e)+vt(i)*tzm1(i,1,1,e) )
         oij(i,e,1) = j*  ! dw/dy - dv/dz
     $   (wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e)+wt(i)*tym1(i,1,1,e)  
     $   -vr(i)*rzm1(i,1,1,e)-vs(i)*szm1(i,1,1,e)-vt(i)*tzm1(i,1,1,e) )

         sij(i,6,e) = j*  ! du/dz + dw/dx
     $   (ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e) +
     $    wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e)+wt(i)*txm1(i,1,1,e) )
         oij(i,e,2) = j*  ! du/dz - dw/dx
     $   (ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e)  
     $   -wr(i)*rxm1(i,1,1,e)-ws(i)*sxm1(i,1,1,e)-wt(i)*txm1(i,1,1,e) )


        enddo
       enddo

      elseif (ifaxis) then  ! AXISYMMETRIC CASE  

c
c        Notation:                       ( 2  x  Acheson, p. 353)
c                     Cylindrical
c            Nek5k    Coordinates
c
c              x          z
c              y          r
c              z          theta
c

         do e=1,nelv
            call setaxdy ( ifrzer(e) )  ! change dytm1 if on-axis
            call local_grad2(ur,us,u,N,e,dxm1,dytm1)
            call local_grad2(vr,vs,v,N,e,dxm1,dytm1)
            call local_grad2(wr,ws,w,N,e,dxm1,dytm1)

            do i=1,nxyz
               j = jacmi(i,e)
               r = ym1(i,1,1,e)                              ! Cyl. Coord:

               sij(i,1,e) = j*  ! du/dx + du/dx              ! e_zz
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))

               sij(i,2,e) = j*  ! dv/dy + dv/dy              ! e_rr
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))

               if (r.gt.0) then                              ! e_@@
                  sij(i,3,e) = 2.*v(i,e)/r  ! v / r  ! corrected AT: factor of 2, 10/30/18
               else
                  sij(i,3,e) = j*  ! L'Hopital's rule: e_@@ = 2dv/dr
     $            2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))
               endif

               sij(i,4,e) = j*  ! du/dy + dv/dx             ! e_zr
     $            ( ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $              vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )
               oij(i,e,3) = j*  ! dv/dx - du/dy
     $            ( vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)  
     $             -ur(i)*rym1(i,1,1,e)-us(i)*sym1(i,1,1,e) )

               if (r.gt.0) then                             ! e_r@
                  sij(i,5,e) = j*  ! dw/dy 
     $              ( wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e) )
     $              - w(i,e) / r
                  oij(i,e,1) = j*  ! dw/dy 
     $              ( wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e) )
     $              + w(i,e) / r
               else
                  sij(i,5,e) = 0
                  oij(i,e,2) = j*  ! L'Hopital's rule: e_r@ = 2dw/dr
     $            2*(wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e))
               endif

               sij(i,6,e) = j*  ! dw/dx                     ! e_@z
     $            ( wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e) )
               oij(i,e,2) = j*  !-dw/dx
     $            (-wr(i)*rxm1(i,1,1,e)-ws(i)*sxm1(i,1,1,e) )

            enddo
         enddo

      else              ! 2D CASE

         do e=1,nelv
            call local_grad2(ur,us,u,N,e,dxm1,dxtm1)
            call local_grad2(vr,vs,v,N,e,dxm1,dxtm1)

            do i=1,nxyz
               j = jacmi(i,e)

               sij(i,1,e) = j*  ! du/dx + du/dx
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))
               oij(i,e,2) = 0.

               sij(i,2,e) = j*  ! dv/dy + dv/dy
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))
               oij(i,e,3) = 0.

               sij(i,3,e) = j*  ! du/dy + dv/dx
     $           (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $            vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )
               oij(i,e,1) = j*  ! dv/dx - du/dy
     $           (vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)  
     $           -ur(i)*rym1(i,1,1,e)-us(i)*sym1(i,1,1,e) )

            enddo
         enddo
      endif

      if(iflomach .and. iflmc) then

        onethird = 1./3.
        if(if3d .or. ifaxis) then
          do e=1,nelv
            do i=1,nxyz
             trS = sij(i,1,e) + sij(i,2,e) + sij(i,3,e) ! 2(du/dx + dv/dy + dw/dz or v/r) in axisym
             sij(i,1,e) = sij(i,1,e) - onethird*trS     ! 2S - (2/3)Q = S'-(1/3)tr(S')
            enddo
          enddo
        else
          do e=1,nelv
            do i=1,nxyz
             trS = sij(i,1,e) + sij(i,2,e)              ! 2(du/dx + dv/dy)
             sij(i,1,e) = sij(i,1,e) - onethird*trS     ! 2S - (2/3)Q = S'-(1/3)tr(S')
            enddo
           enddo
        endif

      endif

      ifielt = ifield
      ifield = 1
      if(ifdss) call dssum_sij_oij(sij, oij, nij)
      ifield = ifielt

      return
      end
c-----------------------------------------------------------------------
      subroutine dssum_sij_oij(sij,oij,nij)
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e
c
      parameter(lt  =lx1*ly1*lz1*lelv)
      parameter(lxyz=lx1*ly1*lz1)

      real sij(lx1*ly1*lz1,nij,lelv), oij(lx1*ly1*lz1,lelv,3)

      common /scrsij/  work1(lx1*ly1*lz1,lelv), work2(lx1*ly1*lz1,lelv)
     $                                        , work3(lx1*ly1*lz1,lelv)

      nxyz = lx1*ly1*lz1
      ntot = nxyz*nelv

      if (if3d) then
         call opcolv  (oij(1,1,1), oij(1,1,2), oij(1,1,3),   bm1)
         call opdssum (oij(1,1,1), oij(1,1,2), oij(1,1,3)       )
         call opcolv  (oij(1,1,1), oij(1,1,2), oij(1,1,3),binvm1)
         do e=1,nelv
            call copy (work1(1,e), sij(1,1,e),nxyz)
            call copy (work2(1,e), sij(1,2,e),nxyz)
            call copy (work3(1,e), sij(1,3,e),nxyz)
         enddo
         call opcolv  (work1, work2, work3,   bm1)
         call opdssum (work1, work2, work3       )
         call opcolv  (work1, work2, work3,binvm1)
         do e=1,nelv
            call copy (sij(1,1,e),work1(1,e), nxyz)
            call copy (sij(1,2,e),work2(1,e), nxyz)
            call copy (sij(1,3,e),work3(1,e), nxyz)
         enddo
         do e=1,nelv
            call copy (work1(1,e), sij(1,4,e),nxyz)
            call copy (work2(1,e), sij(1,5,e),nxyz)
            call copy (work3(1,e), sij(1,6,e),nxyz)
         enddo
         call opcolv  (work1, work2, work3,   bm1)
         call opdssum (work1, work2, work3       )
         call opcolv  (work1, work2, work3,binvm1)
         do e=1,nelv
            call copy (sij(1,4,e),work1(1,e), nxyz)
            call copy (sij(1,5,e),work2(1,e), nxyz)
            call copy (sij(1,6,e),work3(1,e), nxyz)
         enddo
      elseif(ifaxis) then
         call col2    (oij(1,1,1),bm1   ,ntot)
         call dssum   (oij(1,1,1),lx1,ly1,lz1)
         call col2    (oij(1,1,1),binvm1,ntot)
         call col2    (oij(1,1,2),bm1   ,ntot)
         call dssum   (oij(1,1,2),lx1,ly1,lz1)
         call col2    (oij(1,1,2),binvm1,ntot)
         call col2    (oij(1,1,3),bm1   ,ntot)
         call dssum   (oij(1,1,3),lx1,ly1,lz1)
         call col2    (oij(1,1,3),binvm1,ntot)
         do e=1,nelv
            call copy (work1(1,e), sij(1,1,e),nxyz)
            call copy (work2(1,e), sij(1,2,e),nxyz)
            call copy (work3(1,e), sij(1,3,e),nxyz)
         enddo
         call col2    (work1     ,bm1   ,ntot)
         call dssum   (work1     ,lx1,ly1,lz1)
         call col2    (work1     ,binvm1,ntot)
         call col2    (work2     ,bm1   ,ntot)
         call dssum   (work2     ,lx1,ly1,lz1)
         call col2    (work2     ,binvm1,ntot)
         call col2    (work3     ,bm1   ,ntot)
         call dssum   (work3     ,lx1,ly1,lz1)
         call col2    (work3     ,binvm1,ntot)
         do e=1,nelv
            call copy (sij(1,1,e),work1(1,e), nxyz)
            call copy (sij(1,2,e),work2(1,e), nxyz)
            call copy (sij(1,3,e),work3(1,e), nxyz)
         enddo
         do e=1,nelv
            call copy (work1(1,e), sij(1,4,e),nxyz)
            call copy (work2(1,e), sij(1,5,e),nxyz)
            call copy (work3(1,e), sij(1,6,e),nxyz)
         enddo
         call col2    (work1     ,bm1   ,ntot)
         call dssum   (work1     ,lx1,ly1,lz1)
         call col2    (work1     ,binvm1,ntot)
         call col2    (work2     ,bm1   ,ntot)
         call dssum   (work2     ,lx1,ly1,lz1)
         call col2    (work2     ,binvm1,ntot)
         call col2    (work3     ,bm1   ,ntot)
         call dssum   (work3     ,lx1,ly1,lz1)
         call col2    (work3     ,binvm1,ntot)
         do e=1,nelv
            call copy (sij(1,4,e),work1(1,e), nxyz)
            call copy (sij(1,5,e),work2(1,e), nxyz)
            call copy (sij(1,6,e),work3(1,e), nxyz)
         enddo
      else
         call col2    (oij(1,1,1),bm1   ,ntot)
         call dssum   (oij(1,1,1),lx1,ly1,lz1)
         call col2    (oij(1,1,1),binvm1,ntot)
         do e=1,nelv
            call copy (work1(1,e), sij(1,1,e),nxyz)
            call copy (work2(1,e), sij(1,2,e),nxyz)
            call copy (work3(1,e), sij(1,3,e),nxyz)
         enddo
c         call opcolv  (work1, work2, work3,   bm1)
c         call opdssum (work1, work2, work3       )
c         call opcolv  (work1, work2, work3,binvm1)
         call col2    (work1     ,bm1   ,ntot)
         call dssum   (work1     ,lx1,ly1,lz1)
         call col2    (work1     ,binvm1,ntot)
         call col2    (work2     ,bm1   ,ntot)
         call dssum   (work2     ,lx1,ly1,lz1)
         call col2    (work2     ,binvm1,ntot)
         call col2    (work3     ,bm1   ,ntot)
         call dssum   (work3     ,lx1,ly1,lz1)
         call col2    (work3     ,binvm1,ntot)
         do e=1,nelv
            call copy (sij(1,1,e),work1(1,e), nxyz)
            call copy (sij(1,2,e),work2(1,e), nxyz)
            call copy (sij(1,3,e),work3(1,e), nxyz)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_sij_oij_mag2(sij,oij,nij,St_mag2,Om_mag2,OiOjSk)
c
c     Compute the square of the magnitude of the stress and rotation tensors
c     St_mag2=2Sij*Sij=S'ij*S'ij/2, S'ij=2Sij, Sij=(dui/dxj+duj/dxi)/2
c     Om_mag2=2Oij*Oij=O'ij*O'ij/2, O'ij=OSij, Oij=(dui/dxj-duj/dxi)/2
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e
c
      parameter(lt  =lx1*ly1*lz1*lelv)
      parameter(lxyz=lx1*ly1*lz1)

      real sij(lx1*ly1*lz1,nij,lelv), oij(lx1*ly1*lz1,lelv,3)

      real              St_mag2(lx1*ly1*lz1,lelv)
     $                , Om_mag2(lx1*ly1*lz1,lelv)
     $                , OiOjSk (lx1*ly1*lz1,lelv)

      common /scrsij/  work1(lx1*ly1*lz1,lelv), work2(lx1*ly1*lz1,lelv)
     $                                        , work3(lx1*ly1*lz1,lelv)
      real tmp1(lxyz), tmp2(lxyz), tmp3(lxyz)

      nxyz    = lx1*ly1*lz1
      ntot    = nxyz*nelv
      two     = 2.
      onehalf = 0.5
      onethird= 1./3.
      oneeight= 1./8.
      beta    = onehalf-onethird

      if (if3d .or. ifaxis) then     ! 3D CASE or axisymmetric

       do e=1,nelv
          call    col3 (St_mag2(1,e), sij(1,1,e), sij(1,1,e),nxyz)
          call addcol3 (St_mag2(1,e), sij(1,2,e), sij(1,2,e),nxyz)
          call addcol3 (St_mag2(1,e), sij(1,3,e), sij(1,3,e),nxyz)
          call    col3 (work1  (1,e), sij(1,4,e), sij(1,4,e),nxyz)
          call addcol3 (work1  (1,e), sij(1,5,e), sij(1,5,e),nxyz)
          call addcol3 (work1  (1,e), sij(1,6,e), sij(1,6,e),nxyz)
          call  add2s2 (St_mag2(1,e), work1(1,e), two,       nxyz)
c
          call    col3 (work1  (1,e), oij(1,e,1), oij(1,e,1),nxyz)
          call    col3 (work2  (1,e), oij(1,e,2), oij(1,e,2),nxyz)
          call    col3 (work3  (1,e), oij(1,e,3), oij(1,e,3),nxyz)
          call    add4 (Om_mag2(1,e), work1(1,e), work2(1,e)
     $                                          , work3(1,e),nxyz)

          call    add3 (tmp1, work2(1,e), work3(1,e), nxyz)
          call    add3 (tmp2, work1(1,e), work3(1,e), nxyz)
          call    add3 (tmp3, work1(1,e), work2(1,e), nxyz)
          call    col3 (OiOjSk(1,e),tmp1, sij(1,1,e), nxyz)
          call addcol3 (OiOjSk(1,e),tmp2, sij(1,2,e), nxyz)
          call addcol3 (OiOjSk(1,e),tmp3, sij(1,3,e), nxyz)

          call    col4 (tmp1, oij(1,e,1), oij(1,e,2), sij(1,4,e), nxyz)
          call addcol4 (tmp1, oij(1,e,2), oij(1,e,3), sij(1,5,e), nxyz)
          call addcol4 (tmp1, oij(1,e,1), oij(1,e,3), sij(1,6,e), nxyz)
          call add2s2  (OiOjSk(1,e), tmp1,-two,                   nxyz)

       enddo
           
      else              ! 2D CASE

       do e=1,nelv
          call    col3 (St_mag2(1,e), sij(1,1,e), sij(1,1,e),nxyz)
          call addcol3 (St_mag2(1,e), sij(1,2,e), sij(1,2,e),nxyz)
          call    col3 (work1  (1,e), sij(1,3,e), sij(1,3,e),nxyz)
          call  add2s2 (St_mag2(1,e), work1(1,e), two,       nxyz)
c
          call    col3 (Om_mag2(1,e), oij(1,e,1), oij(1,e,1),nxyz)
          call   rzero (OiOjSk (1,e),                        nxyz)
       enddo

      endif

      call  cmult (St_mag2, onehalf, ntot) ! St_mag2=2*Sij*Sij=S'ij*S'ij/2

      return
      end
c-----------------------------------------------------------------------
      subroutine sqrt_tau(tausq,tinput,n)
      implicit real(a-h,o-z)
      include 'SIZE'
      include 'PARALLEL'
c
      real tausq(n), tinput(n)

      do i=1,n
         tau     = tinput(i)               ! Current tau    values
         tausq(i) = sqrt(tau)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine limit_komg
      implicit real(a-h,o-z)
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
c
      integer e

      nxyz       = lx1*ly1*lz1
      ntot       = nxyz*nelv
      nome_neg   = 0
      nkey_neg   = 0
      xome_neg   = 0.
      xkey_neg   = 0.
      frac       = 0.01

c      if(nid.eq.0) write(*,*) 'loglevel is ', loglevel

c limits for k, omega

      do e=1,nelv
      do i=1,nxyz

          omega   = t(i,1,1,e,ifld_omega-1) + f_omegb(i,1,1,e) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c           write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
            xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
            nome_neg = nome_neg + 1
c           write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
            t(i,1,1,e,ifld_omega-1)=frac*abs(t(i,1,1,e,ifld_omega-1))
          endif

          if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c           write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
            xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
            nkey_neg = nkey_neg + 1
c           write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
            t(i,1,1,e,ifld_k-1)=frac*abs(t(i,1,1,e,ifld_k-1))
          endif

      enddo
      enddo

c      if(loglevel.gt.2) then
        nome_neg =iglsum(nome_neg,1)
        nkey_neg =iglsum(nkey_neg,1)
        xome_neg = glmin(xome_neg,1)
        xkey_neg = glmin(xkey_neg,1)

        if(nid.eq.0 .and. nome_neg.gt.0)
     $    write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nid.eq.0 .and. nkey_neg.gt.0)
     $    write(*,*) 'Neg TKE   ', nkey_neg, xkey_neg
c      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine limit_komg_noreg
      implicit real(a-h,o-z)
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
c
      integer e

      nxyz       = lx1*ly1*lz1
      ntot       = nxyz*nelv
      nome_neg   = 0
      nkey_neg   = 0
      xome_neg   = 0.
      xkey_neg   = 0.
      frac       = 0.01

c      if(nid.eq.0) write(*,*) 'loglevel is ', loglevel

c limits for k, omega

      do e=1,nelv
      do i=1,nxyz

          omega   = t(i,1,1,e,ifld_omega-1) ! Current k & omega values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c           write(*,*) 'Zero OMEG ', t(i,1,1,e,ifld_omega-1)
            xome_neg = min(xome_neg,t(i,1,1,e,ifld_omega-1))
            nome_neg = nome_neg + 1
c           write(*,*) 'Neg  OMEG ', nome_neg, xome_neg
            t(i,1,1,e,ifld_omega-1)=frac*abs(t(i,1,1,e,ifld_omega-1))
          endif

          if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c           write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
            xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
            nkey_neg = nkey_neg + 1
c           write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
            t(i,1,1,e,ifld_k-1)=frac*abs(t(i,1,1,e,ifld_k-1))
          endif

      enddo
      enddo

c      if(loglevel.gt.2) then
        nome_neg =iglsum(nome_neg,1)
        nkey_neg =iglsum(nkey_neg,1)
        xome_neg = glmin(xome_neg,1)
        xkey_neg = glmin(xkey_neg,1)

        if(nid.eq.0 .and. nome_neg.gt.0)
     $    write(*,*) 'Neg Omega ', nome_neg, xome_neg
        if(nid.eq.0 .and. nkey_neg.gt.0)
     $    write(*,*) 'Neg TKE   ', nkey_neg, xkey_neg
c      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine limit_ktau
      implicit real(a-h,o-z)
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
c
      integer e

      nxyz       = lx1*ly1*lz1
      ntot       = nxyz*nelv
      ntau_neg   = 0
      nkey_neg   = 0
      xtau_neg   = 0.
      xkey_neg   = 0.
      frac       = 1.

c limits for k, omega

      do e=1,nelv
      do i=1,nxyz

c limits for k, tau

          tau     = t(i,1,1,e,ifld_omega-1) ! Current k & tau    values
          k       = t(i,1,1,e,ifld_k  -1)   ! from previous timestep

          if(t(i,1,1,e,ifld_omega-1).lt.0.0) then
c           write(*,*) 'Zero TAU  ', t(i,1,1,e,ifld_omega-1)
            xtau_neg = min(xtau_neg,t(i,1,1,e,ifld_omega-1))
            ntau_neg = ntau_neg + 1
c           write(*,*) 'Neg  TAU  ', ntau_neg, xtau_neg
            t(i,1,1,e,ifld_omega-1)=frac*abs(t(i,1,1,e,ifld_omega-1))
          endif

          if(t(i,1,1,e,ifld_k-1).lt.0.0) then
c           write(*,*) 'Zero K    ', t(i,1,1,e,ifld_k-1)
            xkey_neg = min(xkey_neg,t(i,1,1,e,ifld_k-1))
            nkey_neg = nkey_neg + 1
c           write(*,*) 'Neg  KEY  ', nkey_neg, xkey_neg
            t(i,1,1,e,ifld_k-1)=frac*abs(t(i,1,1,e,ifld_k-1))
          endif

      enddo
      enddo

c      if(loglevel.gt.2) then
        ntau_neg =iglsum(ntau_neg,1)
        nkey_neg =iglsum(nkey_neg,1)
        xtau_neg = glmin(xtau_neg,1)
        xkey_neg = glmin(xkey_neg,1)

        if(nid.eq.0 .and. ntau_neg.gt.0)
     $    write(*,*) 'Neg Tau   ', ntau_neg, xtau_neg
        if(nid.eq.0 .and. nkey_neg.gt.0)
     $    write(*,*) 'Neg TKE   ', nkey_neg, xkey_neg
c      endif

      return
      end
c-----------------------------------------------------------------------
