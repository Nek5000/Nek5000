c---------------------------------------------------------------------
      subroutine rans_init_wf(wid)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'
      
      integer wid
      
      real w1,w2,w3,w4,w5
      common /SCRNS/
     & w1(lx1*ly1*lz1*lelv)
     &,w2(lx1*ly1*lz1*lelv)
     &,w3(lx1*ly1*lz1*lelv)
     &,w4(lx1*ly1*lz1*lelv)
     &,w5(lx1*ly1*lz1*lelv)

      character*3 bcw

      integer ifc,ie
      
      common /gradywd/ ywdx(lx1,ly1,lz1,lelv),ywdy(lx1,ly1,lz1,lelv)
     $                ,ywdz(lx1,ly1,lz1,lelv),ywdc(lx1,ly1,lz1,lelv)
      real ywdx, ywdy, ywdz, ywdc

c     fix corners
      if (ifstrs) call fixmska1 (v1mask,v2mask,v3mask)
      call fixcorners('shl','W  ')
      call fixcorners('shl','A  ')
      
c     set cbc array for k and omega/tau
      do 10 ie = 1,nelv
      do 10 ifc = 1,2*ndim
         bcw=cbc(ifc,ie,1)
         if(bcw.eq.'W  '.or.bcw.eq.'v  ') then
            cbc(ifc,ie,ifld_k)='t  '
            cbc(ifc,ie,ifld_omega)='t  '
         elseif(bcw.eq.'shl') then
            cbc(ifc,ie,ifld_k)='f  '
            cbc(ifc,ie,ifld_omega)='f  '
         elseif(bcw.eq.'SYM'.or.bcw.eq.'O  '
     $           .or.bcw.eq.'o  ' .or. bcw.eq.'A  ') then
            cbc(ifc,ie,ifld_k)='I  '
            cbc(ifc,ie,ifld_omega)='I  '
         elseif(bcw.eq.'E  ')then
            cbc(ifc,ie,ifld_k)='E  '
            cbc(ifc,ie,ifld_omega)='E  '
         elseif(bcw.eq.'P  ')then
            cbc(ifc,ie,ifld_k)='P  '
            cbc(ifc,ie,ifld_omega)='P  '
         endif
 10   continue      

      if(wid.eq.2)then
         call distf2(ywd,1,'shl','W  ',w1,w2,w3,w4,w5)
      endif

      call get_gradywd
      
      return
      end
c---------------------------------------------------------------------
      subroutine ktau_wf(ix,iy,iz,iside,e,ifpwf)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      integer ix,iy,iz,iside,e
      real tw1,tw2,tauwall,flux_tau,tkewall

      real wd
      common /walldist/ wd(lx1,ly1,lz1,lelv)
      
      common /gradywd/ ywdx(lx1,ly1,lz1,lelv),ywdy(lx1,ly1,lz1,lelv)
     $                ,ywdz(lx1,ly1,lz1,lelv),ywdc(lx1,ly1,lz1,lelv)
      real ywdx, ywdy, ywdz, ywdc
      
      real ydx,ydy,ydz,ydn
      real usn(3)
      real sgnydn

      character*3 cbv3, cbt3, cbk3, cbo3

      logical ifpwf

      integer icalld
      save icalld
      data icalld /0/

      if(icalld.eq.0)then
         icalld = 1
         call init_wf_param
      endif
      
      cbv3 = cbc(iside,e,1)
      cbt3 = cbc(iside,e,2)
      cbk3 = cbc(iside,e,3)
      cbo3 = cbc(iside,e,4)

      call getSnormal  (usn,ix,iy,iz,iside,e)
      ydx = ywdx(ix,iy,iz,e)
      ydy = ywdy(ix,iy,iz,e)
      ydz = ywdz(ix,iy,iz,e)
      if(if3d) then
         ydn= usn(1)*ydx+ usn(2)*ydy+ usn(3)*ydz
      else
         ydn= usn(1)*ydx+ usn(2)*ydy
      endif

      sgnydn = sign(1.,ydn)

      if(cbv3.eq.'shl')then
         if(.not.ifpwf)then
            call standard_ktau_wf(ix,iy,iz,iside,e,tw1,tw2,
     $           tkewall,tauwall,flux_tau)
         else
            call pcorrected_ktau_wf(ix,iy,iz,iside,e,tw1,tw2,
     $           tkewall,tauwall,flux_tau)
         endif
         flux_tau = flux_tau*sgnydn
      endif

      if(ifield.eq.1)then
         if(cbv3.eq.'shl')then
            trn = 0
            tr1 = -tw1
            tr2 = -tw2
         endif
      elseif(ifield.eq.2)then
         temp = 0.0
      elseif(ifield.eq.3)then
         if(cbk3.eq.'t  ' .and. cbv3.eq.'shl') then
            temp= tkewall
         elseif(cbk3.eq.'t  ' .and. cbv3.eq.'W  ') then
            temp= 0.
         elseif(cbk3.eq.'f  ' .and. cbv3.eq.'shl') then
            flux= 0.0
         endif
      elseif (ifield.eq.4) then
         if(cbo3.eq.'t  ' .and. cbv3.eq.'shl') then
            temp = tauwall 
         elseif(cbo3.eq.'t  ' .and. cbv3.eq.'W  ') then
            temp= 0.
         elseif(cbo3.eq.'f  ' .and. cbv3.eq.'shl') then
            flux = flux_tau
        endif
      endif
      
      
      return
      end
c---------------------------------------------------------------------
      subroutine standard_ktau_wf(ix,iy,iz,iside,e,tw1,tw2,
     $     tke,tau,flux_tau)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer ix,iy,iz,e,iside
      real tw1,tw2,tau,flux_tau

      real tsn(3),bsn(3),usn(3)
      real tol,tolmax
      real visc,dens
      real kw,tauw,u_k,tke,omega
      
      real unormal,utx,uty,utz,uw,ut1,ut2
      real utau,utau1,utau2

      real ukstar,veddy,factro,u1plusc
            
      real yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      common /wfparam/ yplusc, Econ, kappa, sCmu, Ccon, alp, bet

      real velx,vely,velz
      
      tol = 1.e-6
      tolmax = 1.e5
      visc = cpfld(1,1)
      dens = cpfld(1,2)

      kw       = abs(ps(1))
      tauw     = abs(ps(2))
      u_k      = sqrt(sCmu*kw)

!     Get the tangent and bi-tangent vectors
      call getangent   (tsn,ix,iy,iz,iside,e)
      call getbitangent(bsn,ix,iy,iz,iside,e)

!     Get the surface normal
      call getSnormal (usn,ix,iy,iz,iside,e)

!     Get the tangential velocity
!     Velocity from previous time step should be tangential
!     Just in case it is not:
      velx = vx(ix,iy,iz,e)
      vely = vy(ix,iy,iz,e)
      velz = vz(ix,iy,iz,e)
      if(if3d)then
         unormal = velx*usn(1)+vely*usn(2)+velz*usn(3)
      else
         unormal = velx*usn(1)+vely*usn(2)
      endif
      utx = velx-unormal*usn(1)
      uty = vely-unormal*usn(2)
      utz = velz-unormal*usn(3)

      ut1=tsn(1)*utx+tsn(2)*uty
      ut2=0.0
      if(if3d) then
         ut1=ut1+tsn(3)*utz
         ut2=bsn(1)*utx+bsn(2)*uty+bsn(3)*utz
      endif
      uw=sqrt(ut1*ut1+ut2*ut2)

      u1plusc = (1./kappa)*log(Econ*yplusc)
      
      utau2 = u_k
      utau1 = uw/u1plusc
      utau = max(utau1,utau2)
            
      tw1 = 0.
      tw2 = 0.
      
      tw1 = (ut1/u1plusc)*utau*dens
      tw2 = (ut2/u1plusc)*utau*dens

      ukstar = utau
            
      veddy = kappa*visc*yplusc !*ukstar/uc
      factro= 0.5 + visc/veddy 

      if(ukstar.ne.0.)then
c         flux_tau = kappa*ukstar*tauw*factro
         flux_tau = kw*tauw*kappa*sCmu*factro/ukstar
      else
         flux_tau = 0.
      endif
      
      tke = ukstar**2/sCmu
      if(veddy.eq.0. .or. veddy.ne.veddy)then
         omega = 0.0
      else
         omega = dens*tke/veddy
      endif
      if(omega.lt.1.)then
         tau = 1e-6
      else
         tau = 1./omega
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine pcorrected_ktau_wf(ix,iy,iz,iside,e,tw1,
     $     tw2,tke,tau,flux_tau)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer ix,iy,iz,e,iside
      real tw1,tw2,tke,flux_tau,tau,omega

      real tsn(3),bsn(3),usn(3)
      real tol,tolmax
      real visc,dens
      real kw,tauw,u_k

      common /pgrads/ dpdx(lx1,ly1,lz1,lelt),dpdy(lx1,ly1,lz1,lelt),
     $     dpdz(lx1,ly1,lz1,lelt)
      real dpdx,dpdy,dpdz
      real dpx,dpy,dpz,dpn
      real dptx,dpty,dptz

      real dpt1,dpt2,dpw,up,up3,ut3
      real unormal,utx,uty,utz,uw,ut1,ut2
      real cosphi,cospsi
      real utau,uc

      real Ccon1, Ccon2
      real u1plusc,u2plusc
      real u1t1,u1t2,ukstar
      real veddy,factro

      real utau1,utau2
      real velx,vely,velz
      
      integer method_ut2,i
            
      real utold(lx1,ly1,lz1,lelv)
      common /oldu/ utold

      real yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      common /wfparam/ yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      
      integer icalld
      save icalld
      data icalld /0/

      method_ut2 = 1 !tomboulides
      method_ut2 = 2 !saini
      
      tol = 1.e-6
      tolmax = 1.e5
      visc = cpfld(1,1)
      dens = cpfld(1,2)

      kw       = abs(ps(1))
      tauw     = abs(ps(2))
      u_k      = sqrt(sCmu*kw)

      if(icalld.eq.0)then
         do i=1,lx1*ly1*lz1*nelt
            utold(i,1,1,1) = sqrt(t(i,1,1,1,2)*sCmu)
         enddo
         icalld = 1
      endif

!     Get the tangent and bi-tangent vectors
      call getangent   (tsn,ix,iy,iz,iside,e)
      call getbitangent(bsn,ix,iy,iz,iside,e)

!     Get the surface normal
      call getSnormal (usn,ix,iy,iz,iside,e)

!     Get the tangent of pressure gradient 
      dpx = dpdx(ix,iy,iz,e)
      dpy = dpdy(ix,iy,iz,e)
      dpz = dpdz(ix,iy,iz,e)
      if(if3d)then !Pressure normal
         dpn = dpx*usn(1)+dpy*usn(2)+dpz*usn(3)
      else
         dpn = dpx*usn(1)+dpy*usn(2)
      endif
      dptx = dpx-dpn*usn(1)
      dpty = dpy-dpn*usn(2)
      dptz = dpz-dpn*usn(3)

!     Get the dot product of press tangent with
!     tangent and bi-tangent vectors
      dpt1 = tsn(1)*dptx+tsn(2)*dpty
      dpt2 = 0.0      
      if(if3d)then
         dpt1 = dpt1+tsn(3)*dptz
         dpt2 = bsn(1)*dptx+bsn(2)*dpty+bsn(3)*dptz
      endif
      dpw = sqrt(dpt1*dpt1+dpt2*dpt2)
      if(dpw.gt.tol)then
         dpt1 = dpt1/dpw
         dpt2 = dpt2/dpw
      else
         dpt1 = 0.
         dpt2 = 0.
         dpw = 0.
      endif

!     Get the pressure velocity scale
      up3 = visc*dpw/dens
      ut3 = utold(ix,iy,iz,e)**3.0
      if(up3/ut3 .gt. 0.001)then
         up = up3**(1./3.)
      else
         up = 0.0
         dpt1 = 0.0
         dpt2 = 0.0
         dpw = 0.0
      endif
      
!     Get the tangential velocity
!     Velocity from previous time step should be tangential
!     Just in case it is not:
      velx = vx(ix,iy,iz,e)
      vely = vy(ix,iy,iz,e)
      velz = vz(ix,iy,iz,e)
      if(if3d)then
         unormal = velx*usn(1)+vely*usn(2)+velz*usn(3)
      else
         unormal = velx*usn(1)+vely*usn(2)
      endif
      utx = velx-unormal*usn(1)
      uty = vely-unormal*usn(2)
      utz = velz-unormal*usn(3)
      
      ut1=tsn(1)*utx+tsn(2)*uty
      ut2=0.0
      if(if3d) then
         ut1=ut1+tsn(3)*utz
         ut2=bsn(1)*utx+bsn(2)*uty+bsn(3)*utz
      endif
      uw=sqrt(ut1*ut1+ut2*ut2)

!     Get utau
      if(uw.ne.0.)then
         cosphi = (dpt1*ut1+dpt2*ut2)/uw
      else
         cosphi = sign(1.,dpt1*ut1+dpt2*ut2)
      endif

      if(.not.if3d)then
         if(cosphi .gt. 0.0)then
            cosphi = 1.0
         else
            cosphi = -1.0
         endif
      endif
      if(cosphi.gt.1.0)cosphi=1.
      if(cosphi.lt.-1.)cosphi=-1.
      
      if(method_ut2.eq.2)then
         utau2 = utold(ix,iy,iz,e)
         call newton_utau(utau2,up,cosphi,uw,kw,if3d,nid)
         utold(ix,iy,iz,e) = utau2 
      elseif(method_ut2.eq.1)then
         utau2 = sqrt(u_k**2.+(up*alp*kappa)**2.
     $        -2.*alp*up*kappa*u_k*cosphi)
c        if(cosphi.eq. 1.) utau2 = abs(u_k-alp*kappa*up)
c        if(cosphi.eq.-1.) utau2 =     u_k+alp*kappa*up
      endif
      
      call finducut(uc,utau1,up,uw,u_k,yplusc,kappa,Ccon,alp,bet,cosphi)

      utau = max(utau1,utau2)

!     Combined velocity scale
      uc = utau+up

      Ccon1   = 0.
      u1plusc = 0.
      if(utau.ne.0. .and. uc.ne.0.) then
         Ccon1   = (utau/uc)*(log(utau/uc)/kappa+ Ccon)
         u1plusc = (utau/uc)* log(  yplusc)/kappa+ Ccon1
      endif
      Ccon2   = 0.
      u2plusc = 0.
      if(up   .ne.0. .and. uc.ne.0.) then
         Ccon2   = (up   /uc)*(log(up   /uc)*alp  + bet )
         u2plusc = (up   /uc)* log(  yplusc)*alp  + Ccon2
      endif

      u1t1 = ut1 - uc*dpt1*u2plusc
      u1t2 = ut2 - uc*dpt2*u2plusc

      if(sqrt(u1t1*u1t1+u1t2*u1t2).ne.0.)then
         cospsi = (u1t1*dpt1+u1t2*dpt2)/sqrt(u1t1*u1t1+u1t2*u1t2)
      else
         cospsi = sign(1.,u1t1*dpt1+u1t2*dpt2)
      endif
      if(cospsi.gt.1.0)cospsi=1.
      if(cospsi.lt.-1.)cospsi=-1.
      
      tw1 = 0.
      tw2 = 0.
      if(uc.ne.0. .and. u1plusc .ne.0.)then
         tw1 = (u1t1/u1plusc)*(utau/uc)*utau*dens
     $        + up**2.*dpt1*(up/uc)*yplusc
         if(if3d)tw2 = (u1t2/u1plusc)*(utau/uc)*utau*dens
     $        + up**2.*dpt2*(up/uc)*yplusc
      endif
      
      ukstar = utau**2. + (up*alp*kappa)**2.
     $     + 2.*utau*up*alp*kappa*cospsi
      ukstar = sqrt(max(0.,ukstar))
      
      veddy = kappa*visc*yplusc*ukstar/uc
      if(veddy.le.tol .or. veddy.ne.veddy) then
         factro= 0.5
      else
         factro= 0.5 + visc/veddy 
      endif

      if(ukstar.gt.tol)then
c         flux_tau = kappa*ukstar*tauw*factro
         flux_tau = kw*tauw*kappa*sCmu*factro/ukstar
      else
         flux_tau = 0.
      endif

      if(tw1.ne.tw1 .or. flux_tau.ne.flux_tau)then
         write(*,*)tw1,flux_tau,utau
         call exit(1)
      endif

      tke = ukstar**2/sCmu
      if(veddy.eq.0. .or. veddy.ne.veddy)then
         omega = 0.0
      else
         omega = dens*tke/veddy
      endif
      if(omega.lt.1.)then
         tau = 0.
      else
         tau = 1./omega
      endif

      return
      end
c---------------------------------------------------------------------
      real function futau(utau,up,cosphi,u,k,if3d)
      implicit none

      real utau,up,cosphi,u,k
      
      logical if3d
            
      real Ccon1,Ccon2
      real u1plusc, u2plusc
      real uc,cospsi,llim

      real yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      common /wfparam/ yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      
!     The lower limit of utau is known
!     Close to zero value would also work
      llim = abs(-up*alp*kappa+sqrt(k*sCmu))
      utau = max(llim,utau)
      uc = utau+up
            
      Ccon1   = 0.
      u1plusc = 0.
      if(utau.ne.0. .and. uc.ne.0.) then
         Ccon1   = (utau/uc)*(log(utau/uc)/kappa+ Ccon)
         u1plusc = (utau/uc)* log(  yplusc)/kappa+ Ccon1
      endif
      Ccon2   = 0.
      u2plusc = 0.
      if(up   .ne.0. .and. uc.ne.0.) then
         Ccon2   = (up   /uc)*(log(up   /uc)*alp  + bet )
         u2plusc = (up   /uc)* log(  yplusc)*alp  + Ccon2
      endif

      cospsi = 1.0
      if(utau.ne.0. .and. uc.ne.0.)then
         cospsi = (u*cosphi - uc*u2plusc)/(u1plusc*uc)
      endif

      if(cospsi.gt.1.0)cospsi = 1.
      if(cospsi.lt.-1.0)cospsi = -1.
            
      futau = utau**2. + (up*alp*kappa)**2.
     $     + 2.*up*alp*kappa*utau*cospsi - k*sCmu

      return
      end
c---------------------------------------------------------------------          
      subroutine newton_utau(utau,up,cosphi,u,k,if3d,nid)
      implicit none

      real utau,cospsi,up,cosphi,u,k

      logical if3d
     
      real tol
      integer  maxiter,i,nid

      real uc,utauplus,upplus
      real f0,f1,f2
      
      real ut0,ut1,futau
      real cpsi1,cpsi2,sgncos

      real yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      common /wfparam/ yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      
      tol = 1e-8
      maxiter = 20

      if(up.eq.0.)then
         utau = sqrt(k*sCmu)
         return
      else
c     Initial guesses for utau
         ut1 = utau 
         ut0 = ut1 + 0.1 
         
         do i=1,maxiter
            f0 = futau(ut0,up,cosphi,u,k,if3d)
            f1 = futau(ut1,up,cosphi,u,k,if3d)

            if(abs(f1-f0).lt.tol)then
               utau = ut0
               exit
            else
               utau = ut0 - (ut1-ut0)*f0/(f1-f0)
               f2 = futau(utau,up,cosphi,u,k,if3d)
               ut0 = ut1
               ut1 = utau
               if(abs(f2).lt.tol)then
                  exit
               endif
            endif
         enddo
      endif

      return
      end
c---------------------------------------------------------------------
      subroutine finducut(uc,utau,up,ut,uk,yplusc,kappa,C,alp,bet,cosf)
c      implicit real (a-h, o-z)
      implicit none

      real uc,ut,up,yplusc,kappa,C,tol,eps,uc1,ucnp1,error,utau
      real uplust, uplusp, cosf, eru
      real func, fder, quant
      real alp, bet
      real u1 , u2 , uk, du1 , du2 
      integer iter, niter
      parameter (niter = 20)
      parameter (tol = 1.e-12, eps=1.e-4)
      real duu2, diffu

      if(up.eq.0) then
         uplust = (    (log(yplusc)             )/kappa+ C  )
         uplusp = 0.
         utau = ut/uplust
         uc   = utau
         goto 20
      endif

      uc = up + 0.01

        do iter = 1,niter
           utau   = uc-up
           uplust = (    (log(yplusc)+log(utau/uc))/kappa+ C  )
           uplusp = (alp*(log(yplusc)+log(up  /uc))      + bet)

           u1     = utau * uplust
           u2     = up   * uplusp
           du1    = uplust + up/uc/kappa
           du2    =-alp*up/uc 

           func  =   u1**2 - u2**2 + 2.*ut*u2*cosf - ut**2
           fder  =2.*u1*du1-2.*(u2-ut*cosf)*du2

           ucnp1 = uc - func/fder
           error = ucnp1-uc

           if(ucnp1.le.0.) then
c              write(*,*) iter,ut,up,uc,ucnp1,utau,cosf, 'Negative1 uc '
              uc = 1.
              goto 12
           endif

c           write(*,'(I4,10G14.7)') iter, ut, up, utau, uc, ucnp1
c     $                           , func, fder ,uplust, uplusp, error

           if(abs(error).le.tol) goto 10
           uc  = ucnp1
 12        continue
        enddo

        if(iter.ge.niter) then
          utau = 0.
          uc   = 0.
c          write(*,*)'No convergence1 for'
c     $            ,ut,up,uc,utau,func,fder,ucnp1,cosf
c             call exitt
          return
        endif
        uc  = ucnp1
 10     continue
        utau = uc-up
 20     continue
        if(cosf.eq. 1.) duu2 = abs(ut-u2)
        if(cosf.eq.-1.) duu2 =    (ut+u2)
        diffu = u1 - duu2
c        write(*,'(I4,10G14.7)') iter, ut, up, utau, uc
c     $                  , uplust, uplusp, error, cosf, diffu

      return
      end
c---------------------------------------------------------------------      
      subroutine fixcorners(cbtype1,cbtype2)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'

      common /scruz/ rmlt(lx1,ly1,lz1,lelv),runx(lx1,ly1,lz1,lelv)
     $              ,runy(lx1,ly1,lz1,lelv),runz(lx1,ly1,lz1,lelv)
     $              ,rt1x(lx1,ly1,lz1,lelv),rt1y(lx1,ly1,lz1,lelv)
     $              ,rt1z(lx1,ly1,lz1,lelv),rt2x(lx1,ly1,lz1,lelv)
     $              ,rt2y(lx1,ly1,lz1,lelv),rt2z(lx1,ly1,lz1,lelv)

      character*3 cb,cbtype1,cbtype2

      nxyz1= lx1*ly1*lz1
      ntot1= nxyz1*nelv
      nfaces = 2*ldim
      tol  = 1.e-01

      call rzero  (rmlt,    ntot1)
      call oprzero(runx,runy,runz)
      call oprzero(rt1x,rt1y,rt1z)
      call oprzero(rt2x,rt2y,rt2z)

c      write(*,*) 'element faces from fixmask2'
      do 1000 iel=1,nelv
      ieg = lglel(iel)
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
              ia =ia + 1
              rmlt(ix,iy,iz,iel)=rmlt(ix,iy,iz,iel)+1.
              runx(ix,iy,iz,iel)=runx(ix,iy,iz,iel)+unx(ia,1,iface,iel)
              runy(ix,iy,iz,iel)=runy(ix,iy,iz,iel)+uny(ia,1,iface,iel)
              rt1x(ix,iy,iz,iel)=rt1x(ix,iy,iz,iel)+t1x(ia,1,iface,iel)
              rt1y(ix,iy,iz,iel)=rt1y(ix,iy,iz,iel)+t1y(ia,1,iface,iel)
              rt2x(ix,iy,iz,iel)=rt2x(ix,iy,iz,iel)+t2x(ia,1,iface,iel)
              rt2y(ix,iy,iz,iel)=rt2y(ix,iy,iz,iel)+t2y(ia,1,iface,iel)
              if(if3d) then
               runz(ix,iy,iz,iel)=runz(ix,iy,iz,iel)+unz(ia,1,iface,iel)
               rt1z(ix,iy,iz,iel)=rt1z(ix,iy,iz,iel)+t1z(ia,1,iface,iel)
               rt2z(ix,iy,iz,iel)=rt2z(ix,iy,iz,iel)+t2z(ia,1,iface,iel)
              endif
 10         continue
         endif
 100  continue
 1000 continue

      call dssum  (rmlt,nx1,ny1,nz1)
      call opdssum(runx, runy, runz)
      call opdssum(rt1x, rt1y, rt1z)
      call opdssum(rt2x, rt2y, rt2z)

      do 2000 iel=1,nelv
      ieg = lglel(iel)
      do 200 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 20 iz=kz1,kz2
            do 20 iy=ky1,ky2
            do 20 ix=kx1,kx2
               ia =ia + 1
               amul = rmlt(ix,iy,iz,iel)
               runxs= runx(ix,iy,iz,iel)/amul
               runys= runy(ix,iy,iz,iel)/amul
               runzs= 0.0
               rt1xs= rt1x(ix,iy,iz,iel)/amul
               rt1ys= rt1y(ix,iy,iz,iel)/amul
               rt1zs= 0.0
               rt2xs= rt2x(ix,iy,iz,iel)/amul
               rt2ys= rt2y(ix,iy,iz,iel)/amul
               rt2zs= 0.0
               if(if3d) then
                  runzs= runz(ix,iy,iz,iel)/amul
                  rt1zs= rt1z(ix,iy,iz,iel)/amul
                  rt2zs= rt2z(ix,iy,iz,iel)/amul
               endif
               unmag = sqrt(runxs*runxs+runys*runys+runzs*runzs)
               t1mag = sqrt(rt1xs*rt1xs+rt1ys*rt1ys+rt1zs*rt1zs)
               t2mag = sqrt(rt2xs*rt2xs+rt2ys*rt2ys+rt2zs*rt2zs)

               if((1.0-abs(unmag)).ge.tol) then
c                 write(*,'(4(1X,A),3I5,2(2X,G14.7))') 'converting BC '
c     $            ,cb, ' to ','W  ', ieg, iface, ia, unmag, amul

                 cbc(iface,iel,1) = 'W  '
               endif

 20         continue
         endif
 200  continue
 2000 continue

      return
      end
c---------------------------------------------------------------------
      subroutine fixcorners_bid(bid1,bid2)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer bid1,bid2,bid,n,iel,ifc,bidmn,bidmx
      integer i,i0,i1,j,j0,j1,k,k0,k1

      real w1,w2,term
      common /scrns/ w1(lx1,ly1,lz1,lelt)
     &              ,w2(lx1,ly1,lz1,lelt)

      n=lx1*ly1*lz1*nelv

      bidmn=min(bid1,bid2)
      bidmx=max(bid1,bid2)

      call rzero(w1,n)
      call rzero(w2,n)
      do iel = 1, nelv
      do ifc = 1, 2*ldim
        bid = BoundaryID(ifc,iel)
        if(bid.gt.0) then
          term = bid
          call cadd_face(ifc,iel,w1,term)
          call cadd_face(ifc,iel,w2,1.0)
        endif
      enddo
      enddo

      call dssum(w1,lx1,ly1,lz1)
      call dssum(w2,lx1,ly1,lz1)

      do iel = 1,nelv
      do ifc = 1,2*ldim
        bid = BoundaryID(ifc,iel)
        if(bid.gt.0) then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do k=k0,k1
          do j=j0,j1
          do i=i0,i1
            term = w1(i,j,k,iel)/w2(i,j,k,iel)
            if(term.gt.bidmn.and.term.lt.bidmx) cbc(ifc,iel,1)='W  '
          enddo
          enddo
          enddo
        endif
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine getangent(st,ix,iy,iz,iside,e)

c     calculate surface normal

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'

      real st(3)
      integer e,f

      f = eface1(iside)

      if (1.le.f.and.f.le.2) then     ! "r face"
         st(1) = t1x(iy,iz,iside,e)
         st(2) = t1y(iy,iz,iside,e)
         st(3) = t1z(iy,iz,iside,e)
      elseif (3.le.f.and.f.le.4) then ! "s face"
         st(1) = t1x(ix,iz,iside,e)
         st(2) = t1y(ix,iz,iside,e)
         st(3) = t1z(ix,iz,iside,e)
      elseif (5.le.f.and.f.le.6) then ! "t face"
         st(1) = t1x(ix,iy,iside,e)
         st(2) = t1y(ix,iy,iside,e)
         st(3) = t1z(ix,iy,iside,e)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine getbitangent(sb,ix,iy,iz,iside,e)

c     calculate surface normal

      include 'SIZE'
      include 'GEOM'
      include 'TOPOL'

      real sb(3)
      integer e,f

      f = eface1(iside)

      if (1.le.f.and.f.le.2) then     ! "r face"
         sb(1) = t2x(iy,iz,iside,e)
         sb(2) = t2y(iy,iz,iside,e)
         sb(3) = t2z(iy,iz,iside,e)
      elseif (3.le.f.and.f.le.4) then ! "s face"
         sb(1) = t2x(ix,iz,iside,e)
         sb(2) = t2y(ix,iz,iside,e)
         sb(3) = t2z(ix,iz,iside,e)
      elseif (5.le.f.and.f.le.6) then ! "t face"
         sb(1) = t2x(ix,iy,iside,e)
         sb(2) = t2y(ix,iy,iside,e)
         sb(3) = t2z(ix,iy,iside,e)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fixmask_old(c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      common /ctmp0/ im1(lx1,ly1,lz1),im2(lx1,ly1,lz1)
      integer e,f,val,im1,im2

      character*3 cb

      n = lx1*ly1*lz1

      nface = 2*ldim

      do e=1,nelv
         call izero (im1,n)
         call izero (im2,n)
         do f=1,nface
            cb  = cbc (f,e,1)
            if (cb.eq.'SYM')  call add_iface_e(im1,f,1,lx1,ly1,lz1)
c            if (cb.eq.'SYM')  call iface_e(im2,f,2,lx1,ly1,lz1)
         enddo
c         call icol2(im2,im1,n)

         tolr = 1.e-8
         k = 1
         do j=1,ly1,ly1-1
         do i=1,lx1,lx1-1
            xx = xm1(i,j,k,e) 
            yy = ym1(i,j,k,e) 
            rr = sqrt(xx*xx+yy*yy)
c            if  ( im2(i,j,k) .eq. 2) then  ! corner of SYM & 'SYM' faces
            if  ( im1(i,j,k) .eq. 2) then  ! corner of SYM & 'SYM' faces
               c1mask(i,j,k,e) = 0.
               c2mask(i,j,k,e) = 0.
            elseif(rr.le.tolr) then
               c1mask(i,j,k,e) = 0.
               c2mask(i,j,k,e) = 0.
            endif
         enddo
         enddo
      enddo

c      do e=1,nelv
c         if ( v1mask(1,1,1,e).eq.0 .and.
c     $        v1mask(2,1,1,e).eq.0 .and.
c     $        v1mask(1,2,1,e).eq.0 ) v2mask(1,1,1,e)=0.
c
c         if ( v1mask(nx1  ,1,1,e).eq.0 .and.
c     $        v1mask(nx1-1,1,1,e).eq.0 .and.
c     $        v1mask(nx1  ,2,1,e).eq.0 ) v2mask(nx1,1,1,e)=0.
c
c         if ( v1mask(1,ny1  ,1,e).eq.0 .and.
c     $        v1mask(1,ny1-1,1,e).eq.0 .and.
c     $        v1mask(2,ny1  ,1,e).eq.0 ) v2mask(1,ny1,1,e)=0.
c
c         if ( v1mask(nx1,ny1  ,1,e).eq.0 .and.
c     $        v1mask(nx1,ny1-1,1,e).eq.0 .and.
c     $        v1mask(nx1-1,ny1,1,e).eq.0 ) v2mask(nx1,ny1,1,e)=0.
c      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine fixmska1(c1mask,c2mask,c3mask)

c     fixes masks for A/shl face corners

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      common /ctmp0/ im1(lx1,ly1,lz1),im2(lx1,ly1,lz1)
      integer e,f,val,im1,im2

      character*3 cb

      n = lx1*ly1*lz1

      nface = 2*ldim

      do e=1,nelv
         call izero (im1,n)
         call izero (im2,n)
         do f=1,nface
            cb  = cbc (f,e,1)
            if (cb.eq.'shl')  call iface_e(im1,f,1,lx1,ly1,lz1)
            if (cb.eq.'A  ')  call iface_e(im2,f,2,lx1,ly1,lz1)
         enddo
         call icol2(im2,im1,n)

         k = 1
         do j=1,ly1,ly1-1
         do i=1,lx1,lx1-1
            if  ( im2(i,j,k) .eq. 2) then  ! corner of shl & 'A  ' faces
               c1mask(i,j,k,e) = 0.
               c2mask(i,j,k,e) = 0.
            endif
         enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
      subroutine add_iface_e(a,iface,val,nx,ny,nz)

C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.

      include 'SIZE'
      integer a(nx,ny,nz),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz)=a(ix,iy,iz)+val
 100      continue
      return
      end

c-----------------------------------------------------------------------
      subroutine un_face_e(u,iface,val,nx,ny,nz)

C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.

      include 'SIZE'
      real a(nx,ny,nz),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz)=a(ix,iy,iz)+val
 100      continue
      return
      end

c-----------------------------------------------------------------------
      subroutine get_gradywd
      include 'SIZE'
      include 'TOTAL'
      include 'RANS_KOMG'

      common /gradywd/ ywdx(lx1,ly1,lz1,lelv),ywdy(lx1,ly1,lz1,lelv)
     $                ,ywdz(lx1,ly1,lz1,lelv),ywdc(lx1,ly1,lz1,lelv)

      ntot = lx1*ly1*lz1*lelv

      call gradm1 (ywdx,ywdy,ywdz,   ywd)
      call opcolv (ywdx,ywdy,ywdz,   bm1)
      call opdssum(ywdx,ywdy,ywdz)
      call opcolv (ywdx,ywdy,ywdz,binvm1)
      call copy   (ywdc,ywd ,ntot)

      return
      end

c-----------------------------------------------------------------------
      subroutine fixmask(cbtype,c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      real   arr(lx1,ly1,lz1,lelv),ar1(lx1,ly1,lz1,lelv)
     $      ,ar2(lx1,ly1,lz1,lelv),ar3(lx1,ly1,lz1,lelv)

      character*3 cb,cbtype

      nxyz1= lx1*ly1*lz1
      ntot1= nxyz1*nelv
      nfaces = 2*ldim
      tol  = 1.e-07

      call rzero  (arr,  ntot1)
      call oprzero(ar1,ar2,ar3)

      write(*,*) 'element faces from fixmask'
      do 1000 iel=1,nelv
      ieg = lglel(iel)
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
               ia =ia + 1
               arr(ix,iy,iz,iel)=arr(ix,iy,iz,iel)+1.
               ar1(ix,iy,iz,iel)=ar1(ix,iy,iz,iel)+unx(ia,1,iface,iel)
               ar2(ix,iy,iz,iel)=ar2(ix,iy,iz,iel)+uny(ia,1,iface,iel)
               if(if3d)
     $         ar3(ix,iy,iz,iel)=ar3(ix,iy,iz,iel)+unz(ia,1,iface,iel)
 10                     continue
         endif
 100       continue
 1000       continue

      call dssum  (arr,nx1,ny1,nz1)
      call opdssum(ar1,  ar2,  ar3)

      do 2000 iel=1,nelv
      ieg = lglel(iel)
      do 200 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 20 iz=kz1,kz2
            do 20 iy=ky1,ky2
            do 20 ix=kx1,kx2
               ia =ia + 1
               amul = arr(ix,iy,iz,iel)
               ar1s = ar1(ix,iy,iz,iel)/amul
               ar2s = ar2(ix,iy,iz,iel)/amul
               ar3s = 0.0

               if(if3d)
     $         ar3s = ar3(ix,iy,iz,iel)/amul
               unmag = sqrt(ar1s*ar1s+ar2s*ar2s+ar3s*ar3s)

               if((1.0-abs(unmag)).ge.tol) then
                 write(*,'(A,3I5,2(2X,G14.7))') 'normal vector '
     $                             , ieg, iface, ia, unmag, amul
                 if(amul.eq.2.) then
                    c1mask(ix,iy,iz,iel) = 0.
                    c2mask(ix,iy,iz,iel) = 0.
                 endif
                 if(amul.eq.3.) then
                    c1mask(ix,iy,iz,iel) = 0.
                    c2mask(ix,iy,iz,iel) = 0.
                    c3mask(ix,iy,iz,iel) = 0.
                 endif
               endif
 20                     continue
         endif
 200       continue
 2000       continue

      return
      end
c-----------------------------------------------------------------------
      subroutine fixmask2(cbtype1,cbtype2,c1mask,c2mask,c3mask)

c     fixes masks for SYM face corners

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'PARALLEL'

      real   c1mask(lx1,ly1,lz1,1)
     $      ,c2mask(lx1,ly1,lz1,1)
     $      ,c3mask(lx1,ly1,lz1,1)

      common /scruz/ rmlt(lx1,ly1,lz1,lelv),runx(lx1,ly1,lz1,lelv)
     $              ,runy(lx1,ly1,lz1,lelv),runz(lx1,ly1,lz1,lelv)
     $              ,rt1x(lx1,ly1,lz1,lelv),rt1y(lx1,ly1,lz1,lelv)
     $              ,rt1z(lx1,ly1,lz1,lelv),rt2x(lx1,ly1,lz1,lelv)
     $              ,rt2y(lx1,ly1,lz1,lelv),rt2z(lx1,ly1,lz1,lelv)

      character*3 cb,cbtype1,cbtype2

      nxyz1= lx1*ly1*lz1
      ntot1= nxyz1*nelv
      nfaces = 2*ldim
      tol  = 1.e-07

      call rzero  (rmlt,    ntot1)
      call oprzero(runx,runy,runz)
      call oprzero(rt1x,rt1y,rt1z)
      call oprzero(rt2x,rt2y,rt2z)

      write(*,*) 'element faces from fixmask2'
      do 1000 iel=1,nelv
      ieg = lglel(iel)
      do 100 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 10 iz=kz1,kz2
            do 10 iy=ky1,ky2
            do 10 ix=kx1,kx2
              ia =ia + 1
              rmlt(ix,iy,iz,iel)=rmlt(ix,iy,iz,iel)+1.
              runx(ix,iy,iz,iel)=runx(ix,iy,iz,iel)+unx(ia,1,iface,iel)
              runy(ix,iy,iz,iel)=runy(ix,iy,iz,iel)+uny(ia,1,iface,iel)
              rt1x(ix,iy,iz,iel)=rt1x(ix,iy,iz,iel)+t1x(ia,1,iface,iel)
              rt1y(ix,iy,iz,iel)=rt1y(ix,iy,iz,iel)+t1y(ia,1,iface,iel)
              rt2x(ix,iy,iz,iel)=rt2x(ix,iy,iz,iel)+t2x(ia,1,iface,iel)
              rt2y(ix,iy,iz,iel)=rt2y(ix,iy,iz,iel)+t2y(ia,1,iface,iel)
              if(if3d) then
               runz(ix,iy,iz,iel)=runz(ix,iy,iz,iel)+unz(ia,1,iface,iel)
               rt1z(ix,iy,iz,iel)=rt1z(ix,iy,iz,iel)+t1z(ia,1,iface,iel)
               rt2z(ix,iy,iz,iel)=rt2z(ix,iy,iz,iel)+t2z(ia,1,iface,iel)
              endif
 10                    continue
         endif
 100       continue
 1000       continue

      call dssum  (rmlt,nx1,ny1,nz1)
      call opdssum(runx, runy, runz)
      call opdssum(rt1x, rt1y, rt1z)
      call opdssum(rt2x, rt2y, rt2z)

      do 2000 iel=1,nelv
      ieg = lglel(iel)
      do 200 iface=1,nfaces
         cb = cbc(iface,iel,1)
         if (cb.eq.cbtype1 .or. cb.eq.cbtype2) then
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            ia = 0
            do 20 iz=kz1,kz2
            do 20 iy=ky1,ky2
            do 20 ix=kx1,kx2
               ia =ia + 1
               amul = rmlt(ix,iy,iz,iel)
               runxs= runx(ix,iy,iz,iel)/amul
               runys= runy(ix,iy,iz,iel)/amul
               runzs= 0.0
               rt1xs= rt1x(ix,iy,iz,iel)/amul
               rt1ys= rt1y(ix,iy,iz,iel)/amul
               rt1zs= 0.0
               rt2xs= rt2x(ix,iy,iz,iel)/amul
               rt2ys= rt2y(ix,iy,iz,iel)/amul
               rt2zs= 0.0

               if(if3d) then
               runzs= runz(ix,iy,iz,iel)/amul
               rt1zs= rt1z(ix,iy,iz,iel)/amul
               rt2zs= rt2z(ix,iy,iz,iel)/amul
               endif
               unmag = sqrt(runxs*runxs+runys*runys+runzs*runzs)
               t1mag = sqrt(rt1xs*rt1xs+rt1ys*rt1ys+rt1zs*rt1zs)
               t2mag = sqrt(rt2xs*rt2xs+rt2ys*rt2ys+rt2zs*rt2zs)

               if((1.0-abs(unmag)).ge.tol) then
                 write(*,'(3(1X,A),3I5,2(2X,G14.7))') 'normal vector '
     $            ,cbtype1, cbtype2, ieg, iface, ia, unmag, amul
                 if    (if3d) then
                   if(amul.eq.2.) then
                    if    ((1.0-abs(t1mag)).ge.tol) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                    elseif((1.0-abs(t2mag)).ge.tol) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c3mask(ix,iy,iz,iel) = 0.
                    endif
                   endif
                   if(amul.eq.3.) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                      c3mask(ix,iy,iz,iel) = 0.
                   endif
                 else
                   if(amul.eq.2.) then
                      c1mask(ix,iy,iz,iel) = 0.
                      c2mask(ix,iy,iz,iel) = 0.
                   endif
                 endif
               endif

 20                     continue
         endif
 200       continue
 2000       continue

      return
      end
c-----------------------------------------------------------------------
      subroutine distf2(d,ifld,b1,b2,dmin,emin,xn,yn,zn)

c     Generate a distance function to boundary with bc "b1" or "b2".
c     This approach does not yet work with periodic boundary conditions.

c     INPUT:  ifld - field type for which distance function is to be found.
c             ifld = 1 for velocity
c             ifld = 2 for temperature, etc.

c     OUTPUT: d = distance to nearest boundary with boundary condition "b"

c     Work arrays:  dmin,emin,xn,yn,zn

      include 'SIZE'
      include 'GEOM'       ! Coordinates
      include 'INPUT'      ! cbc()
      include 'TSTEP'      ! nelfld
      include 'PARALLEL'   ! gather-scatter handle for field "ifld"

      real d(lx1,ly1,lz1,lelt)
      character*3 b1, b2

      real dmin(lx1,ly1,lz1,lelt),emin(lx1,ly1,lz1,lelt)
      real xn(lx1,ly1,lz1,lelt),yn(lx1,ly1,lz1,lelt)
      real zn(lx1,ly1,lz1,lelt)


      integer e,eg,f

      nel = nelfld(ifld)
      n = lx1*ly1*lz1*nel

      call domain_size(xmin,xmax,ymin,ymax,zmin,zmax)

      xmn = min(xmin,ymin)
      xmx = max(xmax,ymax)
      if (if3d) xmn = min(xmn ,zmin)
      if (if3d) xmx = max(xmx ,zmax)

      big = 10*(xmx-xmn)
      call cfill (d,big,n)

      call opcopy(xn,yn,zn,xm1,ym1,zm1)

      nface = 2*ldim
      do e=1,nel     ! Set d=0 on walls
      do f=1,nface
         if (cbc(f,e,1).eq.b1 .or. cbc(f,e,1).eq.b2)
     $              call facev(d,e,f,0.,lx1,ly1,lz1)
      enddo
      enddo

      nxyz = lx1*ly1*lz1

      do ipass=1,10000
         dmax    = 0
         nchange = 0
         do e=1,nel
            do k=1,lz1
            do j=1,ly1
            do i=1,lx1
              i0=max(  1,i-1)
              j0=max(  1,j-1)
              k0=max(  1,k-1)
              i1=min(lx1,i+1)
              j1=min(ly1,j+1)
              k1=min(lz1,k+1)
              do kk=k0,k1
              do jj=j0,j1
              do ii=i0,i1

               dself  = d(i,j,k,e)
               dneigh = d(ii,jj,kk,e)
               if (dneigh.lt.dself) then  ! check neighbor's nearest point
                  d2 = (xm1(i,j,k,e)-xn(ii,jj,kk,e))**2
     $               + (ym1(i,j,k,e)-yn(ii,jj,kk,e))**2
                  if (if3d) d2 = d2 + (zm1(i,j,k,e)-zn(ii,jj,kk,e))**2
                  if (d2.gt.0) d2 = sqrt(d2)
                  if (d2.lt.dself) then
                    nchange = nchange+1
                    d (i,j,k,e) = d2
                    xn(i,j,k,e) = xn(ii,jj,kk,e)
                    yn(i,j,k,e) = yn(ii,jj,kk,e)
                    zn(i,j,k,e) = zn(ii,jj,kk,e)
                    dmax = max(dmax,d(i,j,k,e))
                  endif
               endif
              enddo
              enddo
              enddo
            enddo
            enddo
            enddo
            re = lglel(e)
            call cfill(emin(1,1,1,e),re,nxyz)
            call copy (dmin(1,1,1,e),d(1,1,1,e),nxyz)
         enddo
         nchange = iglsum(nchange,1)
         call fgslib_gs_op(gsh_fld(ifld),dmin,1,3,0) ! min over all elements
         nchange2=0
         do e=1,nel
         do i=1,nxyz
          if (dmin(i,1,1,e).ne.d(i,1,1,e)) then
             nchange2 = nchange2+1
             emin(i,1,1,e) = 0  ! Flag
          endif
         enddo
         enddo
         call copy(d,dmin,n)                !   Ensure updated distance
         nchange2 = iglsum(nchange2,1)
         nchange  = nchange + nchange2
         call fgslib_gs_op(gsh_fld(ifld),emin,1,4,0) ! max over all elements
         do e=1,nel    ! Propagate nearest wall points
         do i=1,nxyz
          eg = emin(i,1,1,e)
          if (eg.ne.lglel(e)) then
             xn(i,1,1,e) = 0
             yn(i,1,1,e) = 0
             zn(i,1,1,e) = 0
          endif
         enddo
         enddo
         call fgslib_gs_op(gsh_fld(ifld),xn,1,1,0) !   Sum over all elements to
         call fgslib_gs_op(gsh_fld(ifld),yn,1,1,0) !   convey nearest point
         call fgslib_gs_op(gsh_fld(ifld),zn,1,1,0) !   to shared neighbor.
         dmax = glmax(dmax,1)
         if (nio.eq.0) write(6,1) ipass,nchange,dmax
 1           format(i9,i12,1pe12.4,' max wall distance 2')
         if (nchange.eq.0) goto 1000
      enddo
 1000  continue
c     wgt = 0.3
c     call filter_d2(d,lx1,lz1,wgt,.true.)
      return
      end
c-----------------------------------------------------------------------
      subroutine finducut1(uc,utau,up,ut,uk,yplusc,kappa,C,alp,bet,cosf)
c      implicit real (a-h, o-z)
      implicit none

      real uc,ut,up,yplusc,kappa,C,tol,eps,uc1,ucnp1,error,utau
      real uplust, uplusp, cosf, eru
      real func, fder, quant
      real alp, bet
      real u1 , u2 , uk, du1 , du2 
      integer iter, niter
      parameter (niter = 20)
      parameter (tol = 1.e-12, eps=1.e-4)

      if(up.eq.0) then
         uplust = (    (log(yplusc)             )/kappa+ C  )
         uplusp = 0.
         utau = ut/uplust
         uc   = utau
         goto 10
      endif
      uc = up + 0.01

        do iter = 1,niter
           utau   = uc-up
           uplust = (    (log(yplusc)+log(utau/uc))/kappa+ C  )
           uplusp = (alp*(log(yplusc)+log(up  /uc))      + bet)

           u1     = utau * uplust
           u2     = up   * uplusp
           du1    = uplust + up/uc/kappa
           du2    =-alp*up/uc 

           func  =   u1 - u2 - ut
           fder  =  du1 -du2

           ucnp1 = uc - func/fder
           error = ucnp1-uc

           if(ucnp1.le.0.) then
c              write(*,*) iter,ut,up,uc,ucnp1,utau,cosf, 'Negative1 uc '
              uc = 1.
              goto 12
           endif

c           write(*,'(I4,10G14.7)') iter, ut, up, utau, uc, ucnp1
c     $                           , func, fder ,uplust, uplusp, error

           if(abs(error).le.tol) goto 10
           uc  = ucnp1
 12        continue
        enddo

        if(iter.ge.niter) then
           utau = 0.
           uc   = up
           return
c        write(*,*)'No convergence1 for'
c     $          ,ut,up,uc,utau,func,fder,ucnp1,cosf
c           call exitt
        endif
        uc  = ucnp1
 10     continue
        utau = uc-up
c 20     write(*,'(A,I4,10G14.7)') 'Converged ', iter, ut, up, utau, uc
c     $                         , uplust, uplusp, error

      return
      end

c-----------------------------------------------------------------------
      subroutine finducut2(uc,utau,up,ut,uk,yplusc,kappa,C,alp,bet,cosf)
c      implicit real (a-h, o-z)
      implicit none

      real uc,ut,up,yplusc,kappa,C,tol,eps,uc1,ucnp1,error,utau
      real uplust, uplusp, cosf, eru
      real func, fder, quant
      real alp, bet
      real u1 , u2 , uk, du1 , du2 
      integer iter, niter
      parameter (niter = 20)
      parameter (tol = 1.e-12, eps=1.e-4)

      if(up.eq.0) then
         uplust = (    (log(yplusc)             )/kappa+ C  )
         uplusp = 0.
         utau = ut/uplust
         uc   = utau
         iter = 0
         goto 30
      endif
      uc = up + 0.01

        do iter = 1,niter
           utau   = uc-up
           uplust = (    (log(yplusc)+log(utau/uc))/kappa+ C  )
           uplusp = (alp*(log(yplusc)+log(up  /uc))      + bet)

           u1     = utau * uplust
           u2     = up   * uplusp
           du1    = uplust + up/uc/kappa
           du2    =-alp*up/uc 

           func  =   u1 - u2 + ut
           fder  =  du1 -du2

           ucnp1 = uc - func/fder
           error = ucnp1-uc

           if(ucnp1.le.0.) then
c              write(*,*) iter,ut,up,uc,ucnp1,utau,cosf, 'Negative1 uc '
              uc = 1.
              goto 12
           endif

c           write(*,'(I4,10G14.7)') iter, ut, up, utau, uc, ucnp1
c     $                           , func, fder ,uplust, uplusp, error

           if(abs(error).le.tol) goto 30
           uc  = ucnp1
 12        continue
        enddo

        if(iter.ge.niter) then
c        write(*,*)'No convergence1 for'
c     $          ,ut,up,uc,utau,func,fder,ucnp1,cosf
           goto 40
        endif
        uc  = ucnp1
 10     continue
        utau = uc-up
        return

 40     continue
        uc = up + 0.01
 
        do iter = 1,niter
           utau   = uc-up
           uplust = (    (log(yplusc)+log(utau/uc))/kappa+ C  )
           uplusp = (alp*(log(yplusc)+log(up  /uc))      + bet)

           u1     = utau * uplust
           u2     = up   * uplusp
           du1    = uplust + up/uc/kappa
           du2    =-alp*up/uc 

           func  =   u1 + u2 - ut
           fder  =  du1 +du2

           ucnp1 = uc - func/fder
           error = ucnp1-uc

           if(ucnp1.le.0.) then
c              write(*,*) iter,ut,up,uc,ucnp1,utau,cosf, 'Negative1 uc '
              uc = 1.
              goto 32
           endif

c           write(*,'(I4,10G14.7)') iter, ut, up, utau, uc, ucnp1
c     $                           , func, fder ,uplust, uplusp, error

           if(abs(error).le.tol) goto 30
           uc  = ucnp1
 32        continue
        enddo

        if(iter.ge.niter) then
           utau = 0.
           uc   = up
           cosf = -1.
           return
        write(*,*)'No convergence1 for'
     $          ,ut,up,uc,utau,func,fder,ucnp1,cosf
c           call exitt
        endif
        uc  = ucnp1
 30     continue
        utau = uc-up
c 20     write(*,'(A,I4,10G14.7)') 'Converged ', iter, ut, up, utau, uc
c     $                         , uplust, uplusp, error

      return
      end

c-----------------------------------------------------------------------
      subroutine finducuk(uc,utau,up,ut,uk,yplusc,kappa,C,alp,bet,cosf)
c      implicit real (a-h, o-z)
      implicit none

      real uc,ut,up,yplusc,kappa,C,tol,eps,uc1,ucnp1,error,utau
      real uplust, uplusp, cosf, eru
      real func, fder, quant
      real alp, bet
      real u1 , u2 , uk, du1 , du2 
      integer iter, niter
      parameter (niter = 20)
      parameter (tol = 1.e-12, eps=1.e-4)

      if(up.eq.0) then
         utau = uk
         uc   = utau
         goto 20
      endif
      uc = up + 0.01

        do iter = 1,niter
           utau   = uc-up
           uplust = (    (log(yplusc)+log(utau/uc))/kappa+ C  )
           uplusp = (alp*(log(yplusc)+log(up  /uc))      + bet)

           u1     = utau * uplust
           u2     = up   * uplusp
           du1    = uplust + up/uc/kappa
           du2    =-alp*up/uc 

           func  =   (ut - u2*cosf)*uk + (u2 - ut*cosf)*alp*kappa*up
     $                   - u1*utau
           fder  =  -(uk*cosf-alp*kappa*up)*du2
     $                - (2.*uplust + up/uc/kappa)*utau

           ucnp1 = uc - func/fder
           error = ucnp1-uc

           if(ucnp1.le.0.) then
c              write(*,*) iter,ut,up,uc,ucnp1,utau,cosf, 'Negative2 uc '
              uc = 1.
              goto 12
           endif

c           write(*,'(I4,10G14.7)') iter, ut, up, utau, uc, ucnp1
c     $                           , func, fder ,uplust, uplusp, error

           if(abs(error).le.tol) goto 10
           uc  = ucnp1
 12        continue
        enddo

        if(iter.ge.niter) then
          uc  = 0.
          utau= 0.
c          write(*,*)'No convergence2 for'
c     $            ,ut,up,uc,utau,func,fder,ucnp1,cosf
          return
c             call exitt
        endif
        uc  = ucnp1
 10     continue
        utau = uc-up
 20     continue
c        write(*,'(I4,10G14.7)') iter, ut, up, utau, uc
c     $                         , uplust, uplusp, error

      return
      end

c-----------------------------------------------------------------------
      subroutine getpgrads
      include 'SIZE'
      include 'TOTAL'

      common /pgrads/ dpdx(lx1,ly1,lz1,lelt),dpdy(lx1,ly1,lz1,lelt),
     $     dpdz(lx1,ly1,lz1,lelt)
      real dpdx,dpdy,dpdz
      
      call gradm1(dpdx,dpdy,dpdz,pr)
      call opcolv(dpdx,dpdy,dpdz,bm1)
      call opdssum(dpdx,dpdy,dpdz)
      call opcolv(dpdx,dpdy,dpdz,binvm1)

      return
      end
c-----------------------------------------------------------------------
      subroutine init_wf_param
      implicit none
      
      real yplusc, Econ, kappa, sCmu, Ccon, alp, bet
      common /wfparam/ yplusc, Econ, kappa, sCmu, Ccon, alp, bet

      yplusc = 30.0
      Econ = 9.0
      kappa = 0.4
      sCmu = 0.3
      Ccon = log(Econ)/kappa
      alp = 5.0
      bet = 8.0

      return
      end
c-----------------------------------------------------------------------      
      subroutine cadd_face(ifc,iel,phi,cc)
      implicit none
      include 'SIZE'

      integer ifc,iel
      real cc, phi(lx1,ly1,lz1,*),ph0

      integer i0,i1,j0,j1,k0,k1,i,j,k

      call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
      do 100 k=k0,k1
      do 100 j=j0,j1
      do 100 i=i0,i1
        ph0 = phi(i,j,k,iel)
        phi(i,j,k,iel) = ph0+cc
 100  continue

      return
      end
C-----------------------------------------------------------------------
