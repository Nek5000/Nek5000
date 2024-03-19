C----------------------------------------------------------------------     
      include "experimental/lshmholtz.f"
C----------------------------------------------------------------------     
      subroutine ls_init(nsteps_cls_in,nsteps_tls_in,
     $                   eps_in,dt_in,
     $                   ifld_cls_in,ifld_clsr_in,
     $                   ifld_tls_in,ifld_tlsr_in,
     $                   ifdebug)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      real eps_in,dt_in
      integer nsteps_cls_in,nsteps_tls_in
      integer ifld_cls_in, ifld_tls_in
      integer ifld_clsr_in, ifld_tlsr_in
      integer ntot, ifdebug
      
      ! multiple of element length
      eps_cls = eps_in
      
      nsteps_cls = nsteps_cls_in
      nsteps_tls = nsteps_tls_in

      ifld_cls = ifld_cls_in
      ifld_clsr = ifld_clsr_in
      ifld_tls = ifld_tls_in
      ifld_tlsr = ifld_tlsr_in

      dt_cls = dt_in

      ifls_debug = ifdebug

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

      return
      end
C----------------------------------------------------------------------     
      subroutine ls_drive(ifld)
      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'LVLSET'

      integer istep_save, ifld_save
      real dt_save, time_save, timef_save
      integer i,ntot,ifld
      integer nsteps_in
      real dtlag_save(10)

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

      ISTEP = 0
      dt = dt_cls
      time = 0.0
      ifield = ifld

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
        !   dt = (0.5*(1.0 + tanh(2.0*PI*(2.*i/nsteps_in-0.5))))*dt_cls
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

        if((ifsvv(ifield-1).and.ifupwindsvv(ifield-1)).or. 
     $       (ifavm(ifield-1).and.ifupwindavm(ifield-1)))then
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

          call hsolve_cls(name4t,ta,tb,h1,h2,
     $                tmask(1,1,1,1,ifield-1),
     $                tmult(1,1,1,1,ifield-1),
     $                imesh,tolht(ifield),nmxt(ifield-1),1,
     $                approxt(1,0,ifld1),napproxt(1,ifld1),binvm1)

          call add2(t(1,1,1,1,ifield-1),ta,n)
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
          VDIFF(i,1,1,1,ifield) = deltael(i,1,1,1) * eps_cls
        enddo
      elseif(ifield.eq.ifld_tlsr)then
        call cfill(vdiff(1,1,1,1,ifield),1e-10,n)
      endif

      if(ifavm(ifield-1))then
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

      integer ix,iy,iz,iel

      real dx(lx1,ly1,lz1,lelt)
      save dx 

      integer icalld
      data icalld /0/
      save icalld

      integer nxyz,n,ie
      real dd,dinv
      real dxmin_e
      real dxmax_e

      nxyz = nx1*ny1*nz1
      n    = nxyz*nelv

      if (icalld.eq.0 .or. ifmvbd) then
         dinv = 1./ldim
         do ie = 1,nelv
            dd = dxmax_e(ie)
            call cfill(dx(1,1,1,ie),dd,nxyz) 
         enddo
         icalld = 1
      endif

      deltael = dx(ix,iy,iz,iel)

      return
      end 
c---------------------------------------------------------------
      real function heaviside(ix,iy,iz,iel,phi,epsin)
      include 'SIZE'
      include 'LVLSET'

      integer ix,iy,iz,iel

      real eps, deltael, phi, epsin

      if(epsin.eq.0.0)then
        eps = deltael(ix,iy,iz,iel)*eps_cls
      else
        eps = deltael(ix,iy,iz,iel)*epsin
      endif
      heaviside = 0.5*(tanh(phi/(2.0*eps))+1.0)

      return
      end
c---------------------------------------------------------------
      subroutine cls_normals(cnx,cny,cnz,ifld)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real cnx(lx1,ly1,lz1,1)
      real cny(lx1,ly1,lz1,1)
      real cnz(lx1,ly1,lz1,1)

      real cmag(lx1,ly1,lz1,lelv)

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
        if(cmag(i,1,1,1).gt.0.)then
          cnx(i,1,1,1) = cnx(i,1,1,1)/cmag(i,1,1,1)
          cny(i,1,1,1) = cny(i,1,1,1)/cmag(i,1,1,1)
          if(if3d)cnx(i,1,1,1) = cnx(i,1,1,1)/cmag(i,1,1,1)
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

      if(ifavm(ifield-1))call axhelm_avm(au,u,imsh,isd)

      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield-1))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield-1)
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
            call col3(tmp1,g2m1(1,1,1,e),clsns(1,1,1,e),nxyz)
            call col3(tmp1,g3m1(1,1,1,e),clsnt(1,1,1,e),nxyz)
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

      if(ifavm(ifield-1))call axhelm_avm(au,u,imsh,isd)

      call addcol4 (au,helm2,bm1,u,ntot)

      if(ifsvv(ifield-1))call axhelm_svv(au,u,imsh,isd)
      !lets worry about axisymmetry later

      if(ifls_debug.eq.1 .and. nio.eq.0)
     $ write(*,*)"SVV status",ifsvv(ifield-1)
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

      phi = t(ix,iy,iz,ie,ifld_tls-1)
      eps = deltael(ix,iy,iz,ie) * eps_cls

      signls = tanh(phi/(2.0 * eps))

      !The TLSR works better with below definition
      !Therefore do not use ifld_tls to define the sign function
      !for TLS re-distancing
      signls = (t(ix,iy,iz,ie,ifld_cls-1)-0.5)*2.0

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

      if(ifsvv(ifield-1))then 
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
      do i=1,ntot
        tb(i,1,1,1) = signls(i,1,1,1)
      enddo
      call col2(clsnx,tb,ntot)
      call col2(clsny,tb,ntot)
      if(if3d)call col2(clsnz,tb,ntot)

      call bdry_tlsr_fix(clsnx,clsny,clsnz)
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
          if(cb.eq.'O  ' .or. cb.eq.'o  ')then

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
