C----------------------------------------------------------------------     
      include "experimental/lshmholtz.f"
C----------------------------------------------------------------------     
      subroutine ls_init(ifld_cls_in, ifld_clsr_in,
     $                   ifld_tls_in, ifld_tlsr_in,
     $                   eps_in, ifdebug, ifixCLSbdry_in)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer ifld_cls_in, ifld_clsr_in
      integer ifld_tls_in, ifld_tlsr_in
      integer ifdebug, ifixCLSbdry_in
      real eps_in


      common /ellength/ dxmax, dxmin
      real dxmax, dxmin

      real deltael, dxave

      real dt_cls_in, dt_tls_in

      integer nsteps_cls_in, nsteps_tls_in

      real nfac

      !get the element lengths
      dxave = deltael(1,1,1,1)

      !Based on unit velocity and shortest element
      !Based on experiment /4 factor gives CFL~0.6
      dt_tls_in = dxmin / lx1 / 4.0

      !Characteristics must travel nfac times largest element
      nfac = 6.0
      nsteps_tls_in = floor(dxmax * nfac /dt_tls_in)

      dt_cls_in = 0.5 * dt_tls_in
      nfac = 0.1
      nsteps_cls_in = floor(dxmax * nfac / dt_cls_in)

      if(nio.eq.0)then 
        write(*,*) "dt - CLSR, TLSR:",dt_cls_in,dt_tls_in
        write(*,*) "nsteps - CLSR, TLSR:",nsteps_cls_in, nsteps_tls_in
      endif

      call ls_init2(nsteps_cls_in, nsteps_tls_in,
     $              eps_in, dt_cls_in, dt_tls_in,
     $              ifld_cls_in, ifld_clsr_in,
     $              ifld_tls_in, ifld_tlsr_in,
     $              ifdebug, ifixCLSbdry_in)

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_init2(nsteps_cls_in,nsteps_tls_in,
     $                   eps_in,dt_cls_in,dt_tls_in,
     $                   ifld_cls_in,ifld_clsr_in,
     $                   ifld_tls_in,ifld_tlsr_in,
     $                   ifdebug,ifixCLSbdry_in)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real eps_in,dt_cls_in,dt_tls_in
      integer nsteps_cls_in,nsteps_tls_in
      integer ifld_cls_in, ifld_tls_in
      integer ifld_clsr_in, ifld_tlsr_in
      integer ntot, ifdebug
      integer ifixCLSbdry_in
      
      ! multiple of element length
      eps_cls = eps_in
      
      nsteps_cls = nsteps_cls_in
      nsteps_tls = nsteps_tls_in

      ifld_cls = ifld_cls_in
      ifld_clsr = ifld_clsr_in
      ifld_tls = ifld_tls_in
      ifld_tlsr = ifld_tlsr_in

      dt_cls = dt_cls_in
      dt_tls = dt_tls_in

      ifls_debug = ifdebug
      ifixCLSbdry = ifixCLSbdry_in

      ntot = lx1*ly1*lz1*nelv
      !no natural BCs for LS fields
      !need to set tmasks to one
      !also turn off internal solvers
      if(ifld_cls_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_cls_in-1),ntot)
      endif
      if(ifld_clsr_in.ne.0)then 
        ! call rone(tmask(1,1,1,1,ifld_clsr_in-1),ntot)
        idpss(ifld_clsr_in-1) = -1
      endif
      if(ifld_tls_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_tls_in-1),ntot)
      endif
      if(ifld_tlsr_in.ne.0)then 
        call rone(tmask(1,1,1,1,ifld_tlsr_in-1),ntot)
        idpss(ifld_tlsr_in-1) = -1
      endif

      if(nio.eq.0)write(*,*)"Initialized Level-Set"
      
      if(nio.eq.0)write(*,*)"Debug mode",ifls_debug

      call fixcorners('shl','W  ')

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_drive(ifld,sgntype)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer istep_save, ifld_save
      real dt_save, time_save, timef_save
      integer i,ntot,ifld
      integer nsteps_in
      real dtlag_save(10)
      integer nbdinp_save
      integer sgntype

      signtype = sgntype 

      if(ifld.ne.ifld_clsr .and. ifld.ne.ifld_tlsr)then
        if(nio.eq.0)then 
          write(*,*)"Solvers only for re-distancing fields"
        endif
        call exit(1)
      endif

      ntot = lx1*ly1*lz1*nelv

      !replace internal istep
      istep_save = ISTEP
      dt_save = dt
      time_save = time
      timef_save = timef
      ifld_save = ifield
      do i=1,10
        dtlag_save(i) = dtlag(i)
      enddo
      nbdinp_save = NBDINP

      ISTEP = 0
      if(ifld.eq.ifld_clsr) dt = dt_cls
      if(ifld.eq.ifld_tlsr) dt = dt_tls
      time = 0.0
      ifield = ifld
      NBDINP = 2

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*) "Field", ifield
        write(*,*) "istep", istep
        write(*,*) "time step", dt_cls
        write(*,*) "Max iteration count", nsteps_cls
      endif

      nsteps_in = 0
      if(ifld.eq.ifld_clsr) nsteps_in = nsteps_cls
      if(ifld.eq.ifld_tlsr) nsteps_in = nsteps_tls

      do i=1,nsteps_in
        ! if(ifld.eq.ifld_tlsr)then
        !   dt = (0.5*(1.0 + tanh(2.0*PI*(3.*i/nsteps_in-0.5))))*dt_tls
        ! endif
        istep = istep + 1
        call ls_advance
      enddo

      !replace back
      ISTEP = istep_save
      dt = dt_save
      time = time_save
      timef = timef_save
      ifield = ifld_save
      do i=1,10
        dtlag(i) = dtlag_save(i)
      enddo
      NBDINP = nbdinp_save

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_advance
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer igeom

      call nekgsync()

      call settime_cls

      call setprop_cls

      if(ifls_debug.eq.1 .and. nio.eq.0)then 
        write(*,*)"ngeom: ",ngeom
      endif

      if (.not.iftmsh(ifield)) imesh = 1
      if (     iftmsh(ifield)) imesh = 2

      do igeom = 1,ngeom
        call unorm
        !diffusion array must be filled out here - pending
        call settolt
        call cdcls(igeom)
      enddo


      return
      end
C----------------------------------------------------------------------     
      subroutine settime_cls
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'LVLSET'

      integer irst, ilag

      irst = param(46)

      do ILAG=10,2,-1
        DTLAG(ILAG) = DTLAG(ILAG-1)
      enddo

      CALL setdt_cls
      DTLAG(1) = DT
      IF (ISTEP.EQ.1 .and. irst.le.0) DTLAG(2) = DT

      TIMEF    = TIME
      TIME     = TIME+DT

      CALL SETORDBD
      if (irst.gt.0) nbd = nbdinp
      CALL RZERO (BD,10)
      CALL SETBD (BD,DTLAG,NBD)
      if (PARAM(27).lt.0) then
        NAB = NBDINP
      else
        NAB = 3
      endif
      IF (ISTEP.lt.NAB.and.irst.le.0) NAB = ISTEP
      CALL RZERO   (AB,10)
      CALL SETABBD (AB,DTLAG,NAB,NBD)

      if(ifls_debug.eq.1 .and. nio.eq.0)then
        write(*,*)"BDF/EXT order",nbd,nab,irst
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine setdt_cls
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real cfl
      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt)
      real ta,tb

      integer i,ntot
      real signls

      ntot = lx1*ly1*lz1*nelv

      !This is redundant
      if(ifield.eq.ifld_clsr)then
        if(istep.eq.1)call cls_normals(clsnx,clsny,clsnz,ifld_tls)
      elseif(ifield.eq.ifld_tlsr)then
        call cls_normals(clsnx,clsny,clsnz,ifld_tlsr)
        do i=1,ntot
          tb(i,1,1,1) = signls(i,1,1,1)
        enddo
        call col2(clsnx,tb,ntot)
        call col2(clsny,tb,ntot)
        if(if3d)call col2(clsnz,tb,ntot)
      endif

      call compute_cfl(cfl,clsnx,clsny,clsnz,dt)

      if(nio.eq.0 .and. istep.eq.1)then
        write(*,*)"CFL: ",cfl
      endif
      ! worry about adjust dt based on CFL later

      return
      end
C----------------------------------------------------------------------     
      subroutine cdcls(igeom)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      include 'ORTHOT'
      include 'SVV'
      include 'AVM'

      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt) 
      real ta,tb

      COMMON /SCRVH/ H1(LX1,LY1,LZ1,LELT),
     $               H2(lx1,ly1,lz1,lelt) 

      real h1,h2

      integer igeom,n,iter
      logical ifconv
      integer ifld1,isd,intype

      n = lx1*ly1*lz1*nelv

      ifld1 = ifield-1
      napproxt(1,ifld1) = laxtt


      if(igeom.eq.1)then
        call makeq_cls
        call lagscal
      else
        if(ifield.eq.ifld_clsr)then
          write(name4t,'(A4)')"CLSR"
        elseif(ifield.eq.ifld_tlsr)then
          write(name4t,'(A4)')"TLSR"
        endif

        if((ifsvv(ifield).and.ifupwindsvv(ifield)).or. 
     $       (ifavm(ifield).and.ifupwindavm(ifield)))then
              call setUpwindSVVAVM(clsnx,clsny,clsnz)
        endif

        isd = 1
        do iter=1,nmxnl
          intype = 0
          if(iftran) intype = -1
          call sethlm_ls(h1,h2,intype)
          call bcneusc(ta,-1)
          call add2(h2,ta,n)
          !following is divergence term
          call add2 (h2,adq(1,1,1,1,ifield-1),n)
          call bcdirsc(t(1,1,1,1,ifield-1))
          call axhelm_cls(ta,t(1,1,1,1,ifield-1),h1,h2,imesh,isd) 
          ! call axhelm_cls2(ta,t(1,1,1,1,ifield-1),h1,h2,imesh,isd) 
          call sub3(tb,bq(1,1,1,1,ifield-1),ta,n)
          call bcneusc(ta,1)
          call add2(tb,ta,n)

          if(ifield.eq.ifld_tlsr)call constrainTLSR(0)
          call hsolve_cls(name4t,ta,tb,h1,h2,
     $                tmask(1,1,1,1,ifield-1),
     $                tmult(1,1,1,1,ifield-1),
     $                imesh,tolht(ifield),nmxt(ifield-1),1,
     $                approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

          call add2(t(1,1,1,1,ifield-1),ta,n)
          if(ifield.eq.ifld_tlsr) call constrainTLSR(1)
          call cvgnlps (ifconv)
          if (ifconv) exit
        enddo
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine setprop_cls
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      include 'AVM'

      integer n,i
      real deltael
      real avm_vdiff

      n = lx1*ly1*lz1*lelv

      call cfill(VTRANS(1,1,1,1,ifield),1.0,n)

      if(ifield.eq.ifld_clsr)then
        do i=1,n
          VDIFF(i,1,1,1,ifield) = deltael(i,1,1,1)*eps_cls/4.0
        enddo
      elseif(ifield.eq.ifld_tlsr)then
        call cfill(vdiff(1,1,1,1,ifield),1e-10,n)
      endif

      if(ifavm(ifield))then
        do i=1,n
          avm_diff(i,1,1,1) = avm_vdiff(i,1,1,1,clsnx,clsny,clsnz)   
        enddo
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine vector_to_rst(ux,uy,uz,ur,us,ut)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real ux(1),uy(1),uz(1)
      real ur(1),us(1),ut(1)

      integer i,ntot
      real xr,yr,zr,xs,ys,zs,xt,yt,zt

      integer icalld
      save icalld
      data icalld /0/

      common /lvlset_x/ xrm1(lx1,ly1,lz1,lelt),
     $                  yrm1(lx1,ly1,lz1,lelt),
     $                  zrm1(lx1,ly1,lz1,lelt),
     $                  xsm1(lx1,ly1,lz1,lelt),
     $                  ysm1(lx1,ly1,lz1,lelt),
     $                  zsm1(lx1,ly1,lz1,lelt),
     $                  xtm1(lx1,ly1,lz1,lelt),
     $                  ytm1(lx1,ly1,lz1,lelt),
     $                  ztm1(lx1,ly1,lz1,lelt)
      real xrm1,yrm1,zrm1 
      real xsm1,ysm1,zsm1
      real xtm1,ytm1,ztm1

      if(icalld.eq.0)then
        call XYZRST(XRM1,YRM1,ZRM1,
     $               XSM1,YSM1,ZSM1,
     $                XTM1,YTM1,ZTM1,ifaxis)
        icalld = 1
      endif
      ntot = lx1*ly1*lz1*nelv

      if(if3d)then
        do i=1,ntot
          xr = xrm1(i,1,1,1)
          yr = yrm1(i,1,1,1)
          zr = zrm1(i,1,1,1)
          xs = xsm1(i,1,1,1)
          ys = ysm1(i,1,1,1)
          zs = zsm1(i,1,1,1)
          xt = xtm1(i,1,1,1)
          yt = ytm1(i,1,1,1)
          zt = ztm1(i,1,1,1)
          ur(i) = xr*ux(i) + yr*uy(i) + zr*uz(i)
          us(i) = xs*ux(i) + ys*uy(i) + zs*uz(i)
          ut(i) = xt*ux(i) + yt*uy(i) + zt*uz(i)
        enddo
      else
        do i=1,ntot
          xr = xrm1(i,1,1,1)
          yr = yrm1(i,1,1,1)
          xs = xsm1(i,1,1,1)
          ys = ysm1(i,1,1,1)
          ur(i) = xr*ux(i) + yr*uy(i)
          us(i) = xs*ux(i) + ys*uy(i)
        enddo
      endif

      return
      end
C----------------------------------------------------------------------     
      subroutine makeq_cls
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ntot

      common /lsscratch2/ du (lx1,ly1,lz1,lelt)
      real du

      integer i
      real signls

      ntot = lx1*ly1*lz1*nelv

      !this is to add divergence and forcing term
      call makeq_aux

      !convop
      if(ifield.eq.ifld_clsr)then
        call conv_clsr(du,t(1,1,1,1,ifield-1))
      elseif(ifield.eq.ifld_tlsr)then
        call conv_tlsr(du,t(1,1,1,1,ifield-1))
      endif

      do i=1,ntot
        bq(i,1,1,1,ifield-1) = bq(i,1,1,1,ifield-1)
     $          -bm1(i,1,1,1)*du(i,1,1,1)*vtrans(i,1,1,1,ifield)
      enddo

      call makeabq

      call makebdq

      return
      end
c---------------------------------------------------------------
      real function deltael(ix,iy,iz,iel)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ix,iy,iz,iel

      real dx(lx1,ly1,lz1,lelt)
      save dx 

      integer icalld
      data icalld /0/
      save icalld

      integer nedge
      parameter(nedge = 4 + 8*(ldim-2))
      real ledg(nedge)

      integer nxyz,n,ie
      real dd,dinv
      real dxmin_e
      real dxmax_e
      real dist_xyzc
      real vlmax
      real glsum, dxsum, glmax, glmin
      integer iglsum

      real delta_save
      save delta_save

      common /ellength/ dxmax, dxmin
      real dxmax, dxmin

      real dxave

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      if (icalld.eq.0 .or. ifmvbd) then
         dinv = 1./ldim
         do ie = 1,nelv
           ledg(1) = dist_xyzc(1,2,ie)
           ledg(2) = dist_xyzc(1,4,ie)
           ledg(3) = dist_xyzc(2,3,ie)
           ledg(4) = dist_xyzc(4,3,ie)
           if (ndim.eq.3) then
             ledg(5)  = dist_xyzc(1,5,ie)
             ledg(6)  = dist_xyzc(2,6,ie)
             ledg(7)  = dist_xyzc(4,8,ie)
             ledg(8)  = dist_xyzc(3,7,ie)

             ledg(9)  = dist_xyzc(5,6,ie)
             ledg(10) = dist_xyzc(5,8,ie)
             ledg(11) = dist_xyzc(8,7,ie)
             ledg(12) = dist_xyzc(6,7,ie)
           endif
            dd = vlmax(ledg,nedge)
            ! dd = dxmax_e(ie)
            call cfill(dx(1,1,1,ie),dd,nxyz) 
         enddo

         dxsum = glsum(dx,n)
         dxmax = glmax(dx,n)
         dxmin = glmin(dx,n)
        
         delta_save = dxsum/iglsum(n,1)

         dxave = delta_save

         if(nio.eq.0)then
           write(*,*)"Max/min/avg el length",dxmax,dxmin,dxave
         endif
         icalld = 1
      endif

      deltael = delta_save !dx(ix,iy,iz,iel)

      return
      end 
c---------------------------------------------------------------
      real function heaviside(ix,iy,iz,iel,phi,epsin)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ix,iy,iz,iel

      real eps, deltael, phi, epsin

      if(epsin.eq.0.0)then
        eps = deltael(ix,iy,iz,iel)*eps_cls
      else
        eps = deltael(ix,iy,iz,iel)*epsin
      endif
      !this factor (/4.0) is introduced so that eps_cls=1
      !gives heaviside transition roughly equal to the
      !element edge
      eps = eps/4.0
      heaviside = 0.5*(tanh(phi/(2.0*eps))+1.0)

      !It looks like CLSR equation is really designed for
      !the tanh heaviside profile rather than the below function
      !Be also wary of the diffusion coefficient for the CLSR
      ! if(phi.ge.eps)then
      !   heaviside = 1.0
      ! elseif(phi.le.-eps)then
      !   heaviside = 0.0
      ! elseif(abs(phi).lt.eps)then
      !   heaviside = 0.5*(1.0 + phi/eps + (1./PI)*sin(PI*phi/eps))
      ! endif

      return
      end
c---------------------------------------------------------------
      subroutine cls_normals(cnx,cny,cnz,ifld)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real cnx(1)
      real cny(1)
      real cnz(1)

      common /cls_norm_temp/ cmag(lx1*ly1*lz1*lelv)
      real cmag

      integer ntot,ifld,i

      ntot = lx1*ly1*lz1*nelv

      !must be calc from TLS field
      call gradm1(cnx,cny,cnz,t(1,1,1,1,ifld-1))
      call opcolv(cnx,cny,cnz,bm1)
      call opdssum(cnx,cny,cnz)
      call opcolv(cnx,cny,cnz,binvm1)

      call col3(cmag,cnx,cnx,ntot)
      call addcol3(cmag,cny,cny,ntot)
      if(if3d) call addcol3(cmag,cnz,cnz,ntot)
      call vsqrt(cmag,ntot)

      do i=1,ntot
        if(cmag(i).gt.1e-12)then
          cnx(i) = cnx(i)/cmag(i)
          cny(i) = cny(i)/cmag(i)
          if(if3d)cnz(i) = cnz(i)/cmag(i)
        else
          cnx(i) = 0.0
          cny(i) = 0.0
          cnz(i) = 0.0
        endif
      enddo

      return
      end
c---------------------------------------------------------------
      subroutine lsmonitor(u,aname)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real u(1)
      character*5 aname
      
      integer n

      real norm, amin, amax
      real gl2norm, glmax, glmin

      n = lx1*ly1*lz1*nelv

      norm = gl2norm(u,n)
      
      amin = glmin(u,n)

      amax = glmax(u,n)

      if(nio.eq.0)then
        write(6,1000)aname," norm, min, max:",
     $   norm,amin,amax
      endif

1000  format(a,10x,a,1p3E13.4)

      return
      end
c---------------------------------------------------------------
      subroutine axhelm_cls2(au,u,helm1,helm2,imsh,isd)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real au(lx1,ly1,lz1,1)
      real u(lx1,ly1,lz1,1)
      real helm1(lx1,ly1,lz1,1)
      real helm2(lx1,ly1,lz1,1)

      integer imsh,isd

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $  ,             DUDS  (LX1,LY1,LZ1)
     $  ,             DUDT  (LX1,LY1,LZ1)
     $  ,             TMP1  (LX1,LY1,LZ1)
     $  ,             TMP2  (LX1,LY1,LZ1)
     $  ,             TMP3  (LX1,LY1,LZ1)
      real dudr,duds,dudt,tmp1,tmp2,tmp3

      real tm1(lx1,ly1,lz1)
      real tm2(lx1,ly1,lz1)
      real tm3(lx1,ly1,lz1)
      equivalence (dudr,tm1),(duds,tm2),(dudt,tm3)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer ntot,nxy,nyz,nxz,nxyz
      integer ie,iz

      common /ltmp/ ur(lx1,ly1,lz1),
     $              us(lx1,ly1,lz1),
     $              ut(lx1,ly1,lz1)
      real ur,us,ut

      real tmp(lx1,ly1,lz1,lelt)
      real tmpx(lx1,ly1,lz1,lelt)
      real tmpy(lx1,ly1,lz1,lelt)
      real tmpz(lx1,ly1,lz1,lelt)

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NELV

      call rzero(au,ntot)

      if(ifield.eq.ifld_clsr)then
        call gradm1(tmpx,tmpy,tmpz,u)
        call opcolv(tmpx,tmpy,tmpz,bm1)
        call opdssum(tmpx,tmpy,tmpz)
        call opcolv(tmpx,tmpy,tmpz,binvm1)

        if(if3d)then
          call vdot3(tmp,tmpx,tmpy,tmpz,clsnx,clsny,clsnz,ntot)
        else
          call vdot2(tmp,tmpx,tmpy,clsnx,clsny,ntot)
        endif

        do ie = 1,nelv
          call vdot2(tmp1,rxm1(1,1,1,ie),rym1(1,1,1,ie),
     $           clsnx(1,1,1,ie),clsny(1,1,1,ie),nxyz)
          call vdot2(tmp2,sxm1(1,1,1,ie),sym1(1,1,1,ie),
     $           clsnx(1,1,1,ie),clsny(1,1,1,ie),nxyz)
          call col2(tmp1,w3m1,nxyz)
          call col2(tmp2,w3m1,nxyz)

          call col2(tmp1,tmp(1,1,1,ie),nxyz)
          call col2(tmp2,tmp(1,1,1,ie),nxyz)

          call col2(tmp1,helm1(1,1,1,ie),nxyz)
          call col2(tmp2,helm1(1,1,1,ie),nxyz)

          call mxm(dxtm1,lx1,tmp1,lx1,tm1,nyz)
          call mxm(tmp2,lx1,dym1,ly1,tm2,ly1)

          call add2(au(1,1,1,ie),tm1,nxyz)
          call add2(au(1,1,1,ie),tm2,nxyz)
        enddo
      endif

      if(ifavm(ifield))call axhelm_avm(au,u,imsh,isd)

      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield)
      if(ifls_debug.eq.1) call lsmonitor(au,'Diff ')

      return
      end
c---------------------------------------------------------------
      subroutine axhelm_cls(au,u,helm1,helm2,imsh,isd)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real au(lx1,ly1,lz1,1)
      real u(lx1,ly1,lz1,1)
      real helm1(lx1,ly1,lz1,1)
      real helm2(lx1,ly1,lz1,1)

      integer imsh,isd

      COMMON /CTMP1/ DUDR  (LX1,LY1,LZ1)
     $  ,             DUDS  (LX1,LY1,LZ1)
     $  ,             DUDT  (LX1,LY1,LZ1)
     $  ,             TMP1  (LX1,LY1,LZ1)
     $  ,             TMP2  (LX1,LY1,LZ1)
     $  ,             TMP3  (LX1,LY1,LZ1)
      real dudr,duds,dudt,tmp1,tmp2,tmp3

      real tm1(lx1,ly1,lz1)
      real tm2(lx1,ly1,lz1)
      real tm3(lx1,ly1,lz1)
      equivalence (dudr,tm1),(duds,tm2),(dudt,tm3)

      COMMON /FASTMD/ IFDFRM(LELT), IFFAST(LELT), IFH2, IFSOLV
      LOGICAL IFDFRM, IFFAST, IFH2, IFSOLV

      integer ntot,nxy,nyz,nxz,nxyz
      integer e,iz

      common /lsaxhelm/ tmpx(lx1,ly1,lz1,lelv),
     $                  tmpy(lx1,ly1,lz1,lelv), 
     $                  tmpz(lx1,ly1,lz1,lelv),
     $                  tmp(lx1,ly1,lz1,lelv)

      real tmpx, tmpy, tmpz, tmp
      real psi
      integer i 

      NXY=lx1*ly1
      NYZ=ly1*lz1
      NXZ=lx1*lz1
      NXYZ=lx1*ly1*lz1
      NTOT=NXYZ*NELV

      call rzero(au,ntot)

      if(ifield.eq.ifld_clsr)then
        call gradm1(tmpx,tmpy,tmpz,u)
        call opcolv(tmpx,tmpy,tmpz,bm1)
        call opdssum(tmpx,tmpy,tmpz)
        call opcolv(tmpx,tmpy,tmpz,binvm1)

        if(if3d)then
          call vdot3(tmp,tmpx,tmpy,tmpz,clsnx,clsny,clsnz,ntot)
        else
          call vdot2(tmp,tmpx,tmpy,clsnx,clsny,ntot)
        endif

        call col2(tmp,helm1,ntot)
        do i=1,ntot
         psi = u(i,1,1,1) 
         tmp(i,1,1,1) = tmp(i,1,1,1)- psi*(1.0-psi)
        enddo

        do e=1,nelv
          if(.not.if3d)then
            call col3(tmp1,g1m1(1,1,1,e),clsnr(1,1,1,e),nxyz)
            call col3(tmp2,g2m1(1,1,1,e),clsns(1,1,1,e),nxyz)
            if (ifdfrm(e)) then
              call addcol3 (tmp1,clsns(1,1,1,e),g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,clsnr(1,1,1,e),g4m1(1,1,1,e),nxyz)
            endif
            call col2 (tmp1,tmp(1,1,1,e),nxyz)
            call col2 (tmp2,tmp(1,1,1,e),nxyz)
            call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
            call mxm  (tmp2,lx1,dym1,ly1,tm2,ly1)
            call add2 (au(1,1,1,e),tm1,nxyz)
            call add2 (au(1,1,1,e),tm2,nxyz)
          else
            call col3(tmp1,g1m1(1,1,1,e),clsnr(1,1,1,e),nxyz)
            call col3(tmp2,g2m1(1,1,1,e),clsns(1,1,1,e),nxyz)
            call col3(tmp3,g3m1(1,1,1,e),clsnt(1,1,1,e),nxyz)
            if (ifdfrm(e)) then
              call addcol3 (tmp1,clsns(1,1,1,e),g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp1,clsnt(1,1,1,e),g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,clsnr(1,1,1,e),g4m1(1,1,1,e),nxyz)
              call addcol3 (tmp2,clsnt(1,1,1,e),g6m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,clsnr(1,1,1,e),g5m1(1,1,1,e),nxyz)
              call addcol3 (tmp3,clsns(1,1,1,e),g6m1(1,1,1,e),nxyz)
            endif
            call col2 (tmp1,tmp(1,1,1,e),nxyz)
            call col2 (tmp2,tmp(1,1,1,e),nxyz)
            call col2 (tmp3,tmp(1,1,1,e),nxyz)
            call mxm  (dxtm1,lx1,tmp1,lx1,tm1,nyz)
            do iz=1,lz1
              call mxm(tmp2(1,1,iz),lx1,dym1,ly1,tm2(1,1,iz),ly1)
            enddo
            call mxm  (tmp3,nxy,dzm1,lz1,tm3,lz1)
            call add2 (au(1,1,1,e),tm1,nxyz)
            call add2 (au(1,1,1,e),tm2,nxyz)
            call add2 (au(1,1,1,e),tm3,nxyz)
          endif
        enddo
      endif

      if(ifavm(ifield))call axhelm_avm(au,u,imsh,isd)

      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield)
      if(ifls_debug.eq.1) call lsmonitor(au,'Diff ')

      return
      end
c---------------------------------------------------------------
      real function signls(ix,iy,iz,ie)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ix,iy,iz,ie
      real deltael,phi,eps

      if(signtype.eq.1)then
        phi = t(ix,iy,iz,ie,ifld_tlsr-1)
        eps = deltael(ix,iy,iz,ie) * eps_cls

        signls = tanh(phi/(2.0 * eps))
      else
      !The TLSR works better with below definition
      !Therefore do not use ifld_tls to define the sign function
      !for TLS re-distancing
        signls = (t(ix,iy,iz,ie,ifld_cls-1)-0.5)*2.0
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine sethlm_ls (h1,h2,intloc)
 
c     Set the variable property arrays H1 and H2
c     in the Helmholtz equation.
c     (associated with variable IFIELD)
c     INTLOC =      integration type

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      include 'LVLSET'

      real h1(1),h2(1)

      nel   = nelfld(ifield)
      ntot1 = lx1*ly1*lz1*nel

      if (iftran) then
         dtbd = bd(1)/dt
         call copy  (h1,vdiff (1,1,1,1,ifield),ntot1)
         if (intloc.eq.0) then
            call rzero (h2,ntot1)
         else
            call cmult2 (h2,vtrans(1,1,1,1,ifield),dtbd,ntot1)
         endif

c        if (ifield.eq.1 .and. ifanls) then   ! this should be replaced
c           const = 2.                        ! with a correct stress
c           call cmult (h1,const,ntot1)       ! formulation
c        endif

      ELSE
         CALL COPY  (H1,VDIFF (1,1,1,1,IFIELD),NTOT1)
         CALL RZERO (H2,NTOT1)
      endif

      if(ifsvv(ifield))then 
        call setmu_svv(t(1,1,1,1,ifield-1),clsnx,clsny,clsnz)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine conv_clsr(du,u)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt)
      real ta,tb

      real u(1),du(1)
      integer ntot

      ntot = lx1*ly1*lz1*nelv

      if(istep.eq.1)then
        !Convert normal vector to rst space
        !Note that normals do not change over re-dist steps
        call cls_normals(clsnx,clsny,clsnz,ifld_tls)
        call vector_to_rst(clsnx,clsny,clsnz,
     $                   clsnr,clsns,clsnt)
      endif

      !(1-psi)*psi

      call rzero(du,ntot)
      return

      call copy(tb,u,ntot)
      call cmult(tb,-1.0,ntot)
      call cadd(tb,1.0,ntot)
      call col2(tb,u,ntot)

      call rzero(ta,ntot)
      call convect_cons(ta,tb,.false.,clsnx,clsny,clsnz,.false.)

      call invcol3(du,ta,bm1,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine conv_tlsr(du,u)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      common /lsscratch/ ta(lx1,ly1,lz1,lelt),
     $                   tb(lx1,ly1,lz1,lelt)
      real ta,tb

      real u(1),du(1)
      integer ntot,i 
      real signls

      ntot = lx1*ly1*lz1*nelv


      !Need to change this to ifld_tlsr later
      !there might exist a novel better solution for this
      !maybe narrow band?
      !the local divergence should give a max of what this should be
      call cls_normals(clsnx,clsny,clsnz,ifld_tlsr)
      call bdry_tlsr_fix(clsnx,clsny,clsnz)

      do i=1,ntot
        tb(i,1,1,1) = signls(i,1,1,1)
      enddo
      call col2(clsnx,tb,ntot)
      call col2(clsny,tb,ntot)
      if(if3d)call col2(clsnz,tb,ntot)

      call convect_new(ta,t(1,1,1,1,ifld_tlsr-1),.false.,
     $                    clsnx,clsny,clsnz,.false.)  

      call invcol3(du,ta,bm1,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine bdry_tlsr_fix(cx,cy,cz)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real cx(lx1,ly1,lz1,1)
      real cy(lx1,ly1,lz1,1)
      real cz(lx1,ly1,lz1,1)

      integer ie,ifc
      integer kx1,kx2,ky1,ky2,kz1,kz2
      real usn(3),udot
      real velx,vely,velz
      character*3 cb
      integer ix,iy,iz

      do ie=1,nelv
        do ifc=1,2*ndim
          cb = cbc(ifc,ie,1)
          if(cb.eq.'O  ' .or. cb.eq.'o  '
     $       .or. cb.eq.'W  '.or. cb.eq.'SYM'
     $        .or. cb.eq.'shl')then

            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,lx1,ly1,lz1,ifc)

            do iz=kz1,kz2
              do iy=ky1,ky2
                do ix=kx1,kx2
                  call getSnormal(usn,ix,iy,iz,ifc,ie)
                  velx = cx(ix,iy,iz,ie)
                  vely = cy(ix,iy,iz,ie)
                  velz = cz(ix,iy,iz,ie)
                  udot = velx*usn(1)+vely*usn(2)
                  if(if3d)udot = udot+velz*usn(3)
                  if(udot .lt. 0.0)then
                    cx(ix,iy,iz,ie)=-cx(ix,iy,iz,ie)
                    cy(ix,iy,iz,ie)=-cy(ix,iy,iz,ie)
                    if(if3d)cz(ix,iy,iz,ie)=-cz(ix,iy,iz,ie)
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine constrainTLSR(op)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      !prevents sign being flipped across the interface
      !with subsequent TLSR iterations
      integer i,ntot
      real sgn,phi
      real sgnold

      integer op
      common /constraintls/ tlsr(lx1,ly1,lz1,lelt)
      real tlsr

      ntot = lx1*ly1*lz1*nelt

      if(op.eq.0)then
        call copy(tlsr,t(1,1,1,1,ifld_tlsr-1),ntot)
      else
        do i=1,ntot
          phi = (t(i,1,1,1,ifld_cls-1)-0.5)*2.0
          sgn = sign(1.,phi)
          phi = tlsr(i,1,1,1)
          sgnold = sign(1.,phi)
          if(sgnold*sgn.lt.0.0) then !Do not update
            t(i,1,1,1,ifld_tlsr-1) = tlsr(i,1,1,1)
          endif
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine surfacetension(ix,iy,iz,e,gamm,sfx,sfy,sfz)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ix,iy,iz,e
      real sfx,sfy,sfz
      common /sforce/ stx(lx1,ly1,lz1,lelv),
     $                sty(lx1,ly1,lz1,lelv), 
     $                stz(lx1,ly1,lz1,lelv), 
     $                curv(lx1,ly1,lz1,lelv),
     $                delta(lx1,ly1,lz1,lelv)

      real stx,sty,stz,curv,delta
      real gamm

      integer ntot

      ntot = lx1*ly1*lz1*nelv

      if(ix*iy*iz*e.eq.1)then

        call deltals(t(1,1,1,1,ifld_tls-1),delta)
        call col2(delta,bm1,ntot)
        call dssum(delta,lx1,ly1,lz1)
        call col2(delta,binvm1,ntot)

        call cls_normals(clsnx,clsny,clsnz,ifld_tls)
        call opdiv(curv,clsnx,clsny,clsnz)
        call dssum(curv,lx1,ly1,lz1)
        call col2(curv,binvm1,ntot)

        call col4(stx,curv,delta,clsnx,ntot)
        call col4(sty,curv,delta,clsny,ntot)
        if(if3d) call col4(stz,curv,delta,clsnz,ntot)

        call cmult(stx,gamm,ntot)
        call cmult(sty,gamm,ntot)
        if(if3d)call cmult(stz,gamm,ntot)
      endif

      !BEWARE: forces are internally multiplied by density
      !So you need to account for that here. 
      !sfx,sfy,sfz MUST be accelerations
      sfx = stx(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)
      sfy = sty(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)
      if(if3d) sfz = stz(ix,iy,iz,e) / vtrans(ix,iy,iz,e,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine deltals(phi,delta)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ntot,i
      real phi(1),delta(1)
      real eps
      real deltael

      ntot = lx1*ly1*lz1*nelv

      do i=1,ntot
        eps = deltael(i,1,1,1)*eps_cls
        if(abs(phi(i)).gt.eps)then
          delta(i) = 0.0
        else
          delta(i) = 0.5*(1.0+cos(PI*phi(i)/eps))/eps
        endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine LS_default_driver(ftlsr,fclsr)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      include 'CTIMER'

      integer ftlsr       !freq of TLS re-distancing
      integer fclsr       !freq of CLS re-initialization

      integer ntot

      integer nclsr
      save nclsr
      data nclsr /0/

      ntot = lx1*ly1*lz1*nelt

      ifcoupledls = .true.

      !re-distancing TLS every n steps
      if(mod(istep,ftlsr).eq.0)then
        call copy(t(1,1,1,1,ifld_tlsr-1),t(1,1,1,1,ifld_cls-1),ntot)

        call cadd(t(1,1,1,1,ifld_tlsr-1),-0.5,ntot)
        
        call cmult(t(1,1,1,1,ifld_tlsr-1),0.01,ntot)

        call ls_drive(ifld_tlsr,1)

        call copy(t(1,1,1,1,ifld_tls-1),t(1,1,1,1,ifld_tlsr-1),ntot)

        nclsr = 0

        ireset_ls = 0
      endif

      call LS_CLS_driver(nclsr, fclsr)

      return
      end
c-----------------------------------------------------------------------
      subroutine LS_CLS_driver(nclsr,fclsr)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'
      include 'CTIMER'

      integer fclsr, nclsr

      integer ntot

      ntot = lx1*ly1*lz1*nelt

      if(mod(nclsr,fclsr).eq.0)then
        call copy(t(1,1,1,1,ifld_clsr-1),t(1,1,1,1,ifld_cls-1),ntot)
        call ls_drive(ifld_clsr)
        call copy(t(1,1,1,1,ifld_cls-1),t(1,1,1,1,ifld_clsr-1),ntot)
        if(ifixCLSbdry .eq. 1)call fixCLSbdryI
        ireset_ls = 0
      endif

      nclsr = nclsr + 1

      return
      end
c-----------------------------------------------------------------------
      real function enclosedVol(direc)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ntot,i
      real glsum

      real vol(lx1,ly1,lz1,lelt)
      integer direc

      ntot = lx1*ly1*lz1*lelv

      if(direc.lt.0)then
        call copy(vol,t(1,1,1,1,ifld_cls-1),ntot)
        call cmult(vol,-1.0,ntot)
        call cadd(vol,1.0,ntot)
      else
        call copy(vol,t(1,1,1,1,ifld_cls-1),ntot)
      endif

      !Convert to step heaviside function
      do i=1,ntot
        if(vol(i,1,1,1) .ge. 0.5)then
          vol(i,1,1,1) = 1.0
        else
          vol(i,1,1,1) = 0.0
        endif
      enddo

      call col2(vol,bm1,ntot)

      enclosedVol = glsum(vol,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine limit_cls(ix,iy,iz,e)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ix,iy,iz,e
      real psi

      if(ifield.eq.ifld_clsr .or. ifield.eq.ifld_cls)then
        psi = t(ix,iy,iz,e,ifield-1)
        if(psi.lt.0.0) psi = 0.0
        if(psi.gt.1.0) psi = 1.0
        t(ix,iy,iz,e,ifield-1) = psi
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine traction_ls(ix,iy,iz,e,iside)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      integer ix,iy,iz,e,iside

      real usn(3),tsn(3),bsn(3)
      real velx,vely,velz
      real utx,uty,utz
      real ut1,ut2,uw
      real yplus,uplus,utau
      real Econ,kappa
      real tw1,tw2
      real unormal,rho

      if(cbc(iside,e,1).eq.'shl')then
        rho = vtrans(ix,iy,iz,e,1)

        yplus = 1.0
        Econ = 9.0
        kappa = 0.41

        call getSnormal(usn,ix,iy,iz,iside,e)
        call getangent(tsn,ix,iy,iz,iside,e)
        call getbitangent(bsn,ix,iy,iz,iside,e)

        !Get the tangential velocity
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
        uw = sqrt(ut1*ut1+ut2*ut2)

        uplus = yplus !(1./kappa)*log(Econ*yplus)

        utau = uw/uplus

        tw1 = 0.0
        tw2 = 0.0

        tw1 = (ut1/uplus)*utau*rho
        tw2 = (ut2/uplus)*utau*rho

        trn = 0.0
        tr1 = -tw1
        tr2 = -tw2
      endif

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
      subroutine fixCLSbdryI
c     This routine will avoid spurious interfaces at the I boundary of CLSR
c     field
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer ie, ifc, ix, iy, iz
      integer kx1,kx2,ky1,ky2,kz1,kz2
      character*3 cb
      real deltael, eps, phi

      do ie=1,nelv
        do ifc=1,2*ndim
          cb = cbc(ifc,ie,ifld_clsr)
          if(cb.eq.'I  ')then
            CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,lx1,ly1,lz1,ifc)
            do iz=kz1,kz2
              do iy=ky1,ky2
                do ix=kx1,kx2
                  eps = deltael(ix,iy,iz,ie) * eps_cls
                  phi = t(ix,iy,iz,ie,ifld_tls-1)
                  if(phi .gt. 2.0*eps)then
                    t(ix,iy,iz,ie,ifld_cls-1) = 1.0
                  elseif(phi .lt. -2.0*eps)then
                    t(ix,iy,iz,ie,ifld_cls-1) = 0.0
                  endif
                enddo
              enddo
            enddo
          endif
        enddo
      enddo

      return
      end
