c-----------------------------------------------------------------------
c     Routines for multidomain (neknek) simulations.
c
c     References:
c
c     Mittal, Ketan, Som Dutta, and Paul Fischer. "Nonconforming
c     Schwarz-spectral element methods for incompressible flow." 
c     Computers & Fluids (2019): 104237.
c
c     "A spectrally accurate method for overlapping grid solution of
c     incompressible Navier–Stokes equations" Brandon E. Merrill,
c     Yulia T. Peet, Paul F. Fischer, and James W. Lottes, J. Comp. Phys.
c     307 (2016) 60-93.
c
c     "Stability analysis of interface temporal discretization in grid
c      overlapping methods," Y. Peet, P.F. Fischer, SIAM J. Numer. Anal.
c      50 (6) (2012) 3375–3401.
c-----------------------------------------------------------------------
      subroutine setup_neknek_wts

      include 'SIZE'
      include 'TOTAL'

      ntot1 = lx1*ly1*lz1*nelt
      call get_ms_dist()
      call getupf()
      call col3(bm1ms,bm1,upf,ntot1)

      return
      end
c-------------------------------------------------------------
      subroutine neknek_setup

      include 'SIZE'
      include 'TOTAL'

      integer icalld
      save    icalld
      data    icalld  /0/

      integer nfld_neknek
      common /inbc/ nfld_neknek

      if (icalld.eq.0.and.nid.eq.0) write(6,*) 'setup neknek'

      if (nsessmax.eq.1) 
     $  call exitti('set nsessmax > 1 in SIZE!$',nsessmax)

      if (icalld.eq.0) then
         ! just in case we call setup from usrdat2 
         call fix_geom
         call geom_reset(1)

         call set_intflag
         call neknekmv
         if (nid.eq.0) write(6,*) 'session id:', idsess
         if (nid.eq.0) write(6,*) 'extrapolation order:', ninter
         if (nid.eq.0) write(6,*) 'nfld_neknek:', nfld_neknek

         nfld_min = iglmin_ms(nfld_neknek,1)
         nfld_max = iglmax_ms(nfld_neknek,1)
         if (nfld_min .ne. nfld_max) then
            nfld_neknek = nfld_min 
            if (nid.eq.0) write(6,*)
     $         'WARNING: reset nfld_neknek to ', nfld_neknek
         endif
      endif

      call setup_int_neknek()  !sets up findpts handle

      call setup_neknek_wts    !sets up integration weights
      
      call exchange_points()   !find donor elements

      if (icalld.eq.0) then
        if(nio.eq.0) write(6,'(A,/)') ' done :: setup neknek'
        icalld = 1
      endif

      return
      end
c-------------------------------------------------------------
      subroutine set_intflag 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      character*3 cb
      character*2 cb2
      equivalence (cb2,cb)
      integer j,e,f

c     Set interpolation flag: points with boundary condition = 'int' 
c     get intflag=1, with 'inp' get intflag=2
c
c     Boundary conditions are changed back to 'v' or 't'.

      nfaces = 2*ldim
      
      nflag=nelt*nfaces
      call izero(intflag,nflag)

      do j=1,nfield
         nel = nelfld(j)
      do e=1,nel
      do f=1,nfaces
         cb=cbc(f,e,j)
         if (cb2.eq.'in') then
            intflag(f,e)=1
            if (j.ge.2) cbc(f,e,j)='t  '
            if (j.eq.1) cbc(f,e,j)='v  '
            if (cb.eq.'inp') then 
               cbc(f,e,j)='o  ' ! Pressure
               intflag(f,e) = 2
            endif
         endif
      enddo
      enddo
      enddo

c      zero out valint
       do i=1,nfld_neknek
         call rzero(valint(1,1,1,1,i),lx1*ly1*lz1*nelt)
       enddo

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      integer k,i,n
      real tvec(3),ccoeff(3)

      if (.not.ifneknekc.or.istep.eq.1) return

      n    = lx1*ly1*lz1*nelt

      do k=1,nfld_neknek
         call copy(bdrylg(1,k,2),bdrylg(1,k,1),n)
         call copy(bdrylg(1,k,1),bdrylg(1,k,0),n)
         call copy(bdrylg(1,k,0),valint(1,1,1,1,k),n)
      enddo

c     Order of extrpolation is contolled by the parameter NINTER contained 
c     in NEKNEK. First order interface extrapolation, NINTER=1 (time lagging) 
c     is activated. It is unconditionally stable.  If you want to use 
c     higher-order interface extrapolation schemes, you need to increase 
c     ngeom to 3-5 for scheme to be stable.

      ichk = 0
      if (iffxdt) ichk = 1
      ichkg = iglmin_ms(ichk,1)
      if (ichkg.eq.1) then 
        if (NINTER.eq.1.or.istep.le.2) then
         c0=1.
         c1=0.
         c2=0.
         else if (NINTER.eq.2.or.istep.eq.3) then
           c0=2.
           c1=-1.
           c2=0.
         else 
           c0=3.
           c1=-3.
           c2=1.
        endif
      else
         tvec(1) = time-dtlag(1)
         tvec(2) = tvec(1)-dtlag(2)
         tvec(3) = tvec(2)-dtlag(3)
         if (NINTER.eq.1.or.istep.le.2) then
           c0=1.
           c1=0.
           c2=0.
         elseif (NINTER.eq.2.or.istep.eq.3) then
           call fd_weights_full(time,tvec,1,0,ccoeff)
           c0=ccoeff(1)
           c1=ccoeff(2)    
           c2=0. 
         else
           call fd_weights_full(time,tvec,2,0,ccoeff)
           c0=ccoeff(1)
           c1=ccoeff(2) 
           c2=ccoeff(3)
         endif
      endif

      do k=1,nfld_neknek
      do i=1,n
         valint(i,1,1,1,k) = 
     $      c0*bdrylg(i,k,0)+c1*bdrylg(i,k,1)+c2*bdrylg(i,k,2)
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------
      subroutine chk_outflow ! Assign neighbor velocity to outflow if
                             ! characteristics are going the wrong way.

      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      integer e,eg,f

      nface = 2*ldim
      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,1).eq.'o  ') then
           eg = lglel(e)
           call facind (i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,f)
           l=0
           do k=k0,k1
           do j=j0,j1
           do i=i0,i1
              l=l+1
              vo=vx(i,j,k,e)*unx(l,1,f,e)
     $          +vy(i,j,k,e)*uny(l,1,f,e)
     $          +vz(i,j,k,e)*unz(l,1,f,e)
              if (vo.lt.0) then            ! We have inflow
                 cbu = cbc(f,e,1)
                 call userbc(i,j,k,f,eg)
                 vx(i,j,k,e) = ux
                 vy(i,j,k,e) = uy
                 vz(i,j,k,e) = uz
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
      subroutine neknekmv()
      include 'SIZE'
      include 'TOTAL'
      integer imove

      imove=1
      if (ifmvbd) imove=0
      iglmove = iglmin_ms(imove,1)

      if (iglmove.eq.0) then
         ifneknekm=.true.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_int_neknek()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      include 'mpif.h'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

c     Get total number of processors and number of p
      npall = 0
      do i=1,nsessions
       npall = npall+npsess(i-1)
      enddo

c     Setup findpts    
      tol     = 5e-13
      npt_max = 128
      nxf     = 2*lx1 ! fine mesh for bb-test
      nyf     = 2*ly1
      nzf     = 2*lz1
      bb_t    = 0.01 ! relative size to expand bounding boxes by
      ntot    = lx1*ly1*lz1*nelt

      if (istep.gt.0) call fgslib_findptsms_free(fpth_ms)
      call fgslib_findptsms_setup(fpth_ms,mpi_comm_world,npall,ldim,
     &                            xm1,ym1,zm1,lx1,ly1,lz1,
     &                            nelt,nxf,nyf,nzf,bb_t,ntot,ntot,
     &                            npt_max,tol,idsess,distfint)

      return
      end
c-----------------------------------------------------------------------
      subroutine exchange_points()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer jsend(nmaxl_nn)
      common /exchr/ rsend(ldim*nmaxl_nn)
      real    dist(nmaxl_nn)
      integer isid_nn(nmaxl_nn)

c     Look for boundary points with Diriclet b.c. (candidates for
c     interpolation)

      ifield = 1
      if (ifheat)  ifield   = 2

      nfaces = 2*ldim
      nel    = nelfld(ifield)
      nxy    = lx1*ly1
      nxyz   = nxy*lz1
      ntot   = nxyz*nel
      call izero(imask,ntot)

c     Setup arrays of x,y,zs to send to findpts and indices of boundary 
c     points in iList
      ip = 0
      do ie=1,nel
      do iface=1,nfaces
         if (intflag(iface,ie).gt.0) then
            call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ip=ip+1
               idx = (ie-1)*nxyz+(iz-1)*nxy+(iy-1)*lx1+ix
               iList(ip) = idx 
               isid_nn(ip) = idsess
               rsend(ldim*(ip-1)+1)=xm1(ix,iy,iz,ie)
               rsend(ldim*(ip-1)+2)=ym1(ix,iy,iz,ie)
               if (if3d) 
     $         rsend(ldim*(ip-1)+3)=zm1(ix,iy,iz,ie)

               if (ip.gt.nmaxl_nn) then
                  write(6,*) nid,
     &            ' ABORT: nbp (current ip) too large',ip,nmaxl_nn
                  call exitt
               endif
               imask(idx,1,1,1)=1
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
      nbp = ip
      nbpg = iglsum(nbp,1)

c     JL's routine to find which points these procs are on
      call fgslib_findptsms(fpth_ms,rcode,1,
     &                      proc,1,
     &                      elid,1,
     &                      rst,ldim,
     &                      dist,1,
     &                      rsend(1),ldim,
     &                      rsend(2),ldim,
     &                      rsend(3),ldim,
     &                      isid_nn,1,0,
     &                      nbp)

      ip     = 0
      ierror = 0
c     Make sure all points were found
      do i=1,nbp
        if (rcode(i).eq.2) then
           ierror = 1
        else 
           ip = ip+1
        endif
      enddo

      ipg = iglsum(ip,1)
      if (nid.eq.0) write(6,*) idsess,nbpg,ibpg,
     $                    ' Interdomain boundary points'
      npoints_nn = nbp

      ierror = iglmax_ms(ierror,1)
      if (ierror.eq.1) call exitti('Not all interdomain boundary points
     $                              were found in domain: $',idsess)
     
 
      return
      end
c-----------------------------------------------------------------------
      subroutine neknek_exchange
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'CTIMER'

      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)

      if (nio.eq.0) write(6,98) 
     $   ' Multidomain data exchange ... ', nfld_neknek
 98   format(12x,a,i3)

      etime0 = dnekclock_sync()
      call neknekgsync()
      etime1 = dnekclock()

      call mappr(pm1,pr,wk1,wk2)  ! Map pressure to pm1 
      nv = lx1*ly1*lz1*nelv
      nt = lx1*ly1*lz1*nelt

c     Interpolate using findpts_eval
      call field_eval(fieldout(1,1),1,vx)
      call field_eval(fieldout(1,2),1,vy)
      if (ldim.eq.3) call field_eval(fieldout(1,ldim),1,vz)
      call field_eval(fieldout(1,ldim+1),1,pm1)
      if (nfld_neknek.gt.ldim+1) then 
        do i=ldim+2,nfld_neknek
          call field_eval(fieldout(1,i),1,t(1,1,1,1,i-ldim-1))
        enddo
      endif
         
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
      do i=1,npoints_nn
       idx = iList(i)
       do ifld=1,nfld_neknek
         valint(idx,1,1,1,ifld)=fieldout(i,ifld)
       enddo
      enddo

      call nekgsync()
      etime = dnekclock() - etime1
      tsync = etime1 - etime0

      if (nio.eq.0) write(6,99) istep,
     $              '  done :: Multidomain data exchange',
     $              etime, etime+tsync
 99   format(i11,a,1p2e13.4)

      return
      end
c--------------------------------------------------------------------------
      subroutine field_eval(fieldout,fieldstride,fieldin)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real fieldout(1)
      real fieldin(1)
      integer fieldstride

      call fgslib_findptsms_eval(fpth_ms,fieldout,fieldstride,
     &                           rcode,1,proc,1,elid,1,rst,ldim,
     &                           npoints_nn,fieldin)

      return
      end
c--------------------------------------------------------------------------
      subroutine neknek_xfer_fld(u,ui)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real fieldout(nmaxl_nn)
      real u(1),ui(1)

      call field_eval(fieldout,1,u)
c     Now we transfer fieldout to ui at the interdomain boundary points
      do i=1,npoints_nn
        idx = iList(i)
        ui(idx)=fieldout(i)
      enddo

      return
      end
C--------------------------------------------------------------------------
      subroutine fix_surface_flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      integer e,f
      common /ctmp1/ work(lx1*ly1*lz1*lelt)
      integer itchk
      common /idumochk/ itchk
      integer icalld
      save    icalld
      data    icalld /0/
c     assume that this routine is called at the end of bcdirvc
c     where all the boundary condition data has been read in for 
c     velocity.
      if (icalld.eq.0) then
       itchk = 0
       do e=1,nelv
       do f=1,2*ldim
         if (cbc(f,e,1).eq.'o  '.or.cbc(f,e,1).eq.'O  ') then
           itchk = 1
         endif
       enddo
       enddo
       itchk = iglmax(itchk,1)
       icalld = 1
      endif

      if (itchk.eq.1) return

      dqg=0
      aqg=0
      do e=1,nelv
      do f=1,2*ldim
         if (cbc(f,e,1).eq.'v  '.or.cbc(f,e,1).eq.'V  ') then
            call surface_flux_area(dq,aq,vx,vy,vz,e,f,work)
            dqg = dqg+dq
            if (intflag(f,e).eq.1) aqg = aqg+aq
         endif
      enddo
      enddo
      dqg=glsum(dqg,1) ! sum over all processors for this session
      aqg=glsum(aqg,1) ! sum over all processors for this session
      gamma = 0.
      if (aqg.gt.0) gamma = -dqg/aqg
      if (nid.eq.0) write(6,104) idsess,istep,time,dqg,aqg,gamma
 104  format(i4,i10,1p4e13.4,' NekNek_bdry_flux')

      do e=1,nelv
      do f=1,2*ldim
        if (intflag(f,e).eq.1) then
          call facind (i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,f)
          l=0
          do k=k0,k1
          do j=j0,j1
          do i=i0,i1
            l=l+1
            vx(i,j,k,e) = vx(i,j,k,e) + gamma*unx(l,1,f,e)
            vy(i,j,k,e) = vy(i,j,k,e) + gamma*uny(l,1,f,e)
            if (ldim.eq.3) 
     $      vz(i,j,k,e) = vz(i,j,k,e) + gamma*unz(l,1,f,e)
          enddo
          enddo
          enddo
        endif
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine surface_flux_area(dq,aq,qx,qy,qz,e,f,w)
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      parameter (l=lx1*ly1*lz1)
      real qx(l,1),qy(l,1),qz(l,1),w(lx1,ly1,lz1)
      integer e,f

      call           faccl3  (w,qx(1,e),unx(1,1,f,e),f)
      call           faddcl3 (w,qy(1,e),uny(1,1,f,e),f)
      if (if3d) call faddcl3 (w,qz(1,e),unz(1,1,f,e),f)
      call dsset(lx1,ly1,lz1)

      iface  = eface1(f)
      js1    = skpdat(1,iface)
      jf1    = skpdat(2,iface)
      jskip1 = skpdat(3,iface)
      js2    = skpdat(4,iface)
      jf2    = skpdat(5,iface)
      jskip2 = skpdat(6,iface)

      dq = 0
      aq = 0
      i  = 0

      do 100 j2=js2,jf2,jskip2
      do 100 j1=js1,jf1,jskip1
          i = i+1
         dq = dq + area(i,1,f,e)*w(j1,j2,1)
         aq = aq + area(i,1,f,e)
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine rescale_x_ms (x,x0,x1)
      include 'SIZE'
      real x(1)

      n = lx1*ly1*lz1*nelt
      xmin = glmin(x,n)
      xmax = glmax(x,n)
      xming = glmin_ms(x,n)
      xmaxg = glmax_ms(x,n)

      if (xmax.le.xmin) return

      scale = (x1-x0)/(xmaxg-xming)
      x0n   = x0 + scale*(xmin-xming)

      do i=1,n
         x(i) = x0n + scale*(x(i)-xmin)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine vol_flow_ms
c
c
c     Adust flow volume at end of time step to keep flow rate fixed by
c     adding an appropriate multiple of the linear solution to the Stokes
c     problem arising from a unit forcing in the X-direction.  This assumes
c     that the flow rate in the X-direction is to be fixed (as opposed to Y-
c     or Z-) *and* that the periodic boundary conditions in the X-direction
c     occur at the extreme left and right ends of the mesh.
c
c     pff 6/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
c     Swap the comments on these two lines if you don't want to fix the
c     flow rate for periodic-in-X (or Z) flow problems.
c
      parameter (kx1=lx1,ky1=ly1,kz1=lz1,kx2=lx2,ky2=ly2,kz2=lz2)
c
      common /cvflow_a/ vxc(kx1,ky1,kz1,lelv)
     $                , vyc(kx1,ky1,kz1,lelv)
     $                , vzc(kx1,ky1,kz1,lelv)
     $                , prc(kx2,ky2,kz2,lelv)
     $                , vdc(kx1*ky1*kz1*lelv,2)
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv
c
      real bd_vflow,dt_vflow
      save bd_vflow,dt_vflow
      data bd_vflow,dt_vflow /-99.,-99./

      logical ifcomp

      ntot1 = lx1*ly1*lz1*nelv
      ntot2 = lx2*ly2*lz2*nelv

      if (param(55).eq.0.) return
      if (kx1.eq.1) then
         write(6,*) 'ABORT. Recompile vol_flow with kx1=lx1, etc.'
         call exitt
      endif

      icvflow   = 1                                  ! Default flow dir. = X
      if (param(54).ne.0) icvflow = abs(param(54))
      iavflow   = 0                                  ! Determine flow rate
      if (param(54).lt.0) iavflow = 1                ! from mean velocity
      flow_rate = param(55)

      chv(1) = 'X'
      chv(2) = 'Y'
      chv(3) = 'Z'

c     If either dt or the backwards difference coefficient change,
c     then recompute base flow solution corresponding to unit forcing:

      ifcomp = .false.
      if (dt.ne.dt_vflow.or.bd(1).ne.bd_vflow.or.ifmvbd) ifcomp=.true.
      if (.not.ifcomp) then
         ifcomp=.true.
         do i=1,ntot1
            if (vdiff (i,1,1,1,1).ne.vdc(i,1)) goto 20
            if (vtrans(i,1,1,1,1).ne.vdc(i,2)) goto 20
         enddo
         ifcomp=.false.  ! If here, then vdiff/vtrans unchanged.
   20    continue
      endif

      call copy(vdc(1,1),vdiff (1,1,1,1,1),ntot1)
      call copy(vdc(1,2),vtrans(1,1,1,1,1),ntot1)
      dt_vflow = dt
      bd_vflow = bd(1)

      if (ifcomp) call compute_vol_soln_ms(vxc,vyc,vzc,prc)

      if (icvflow.eq.1)
     $       current_flow=glsc2_ms(vx,bm1ms,ntot1)/domain_length  ! for X
      if (icvflow.eq.2)
     $       current_flow=glsc2_ms(vy,bm1ms,ntot1)/domain_length  ! for Y
      if (icvflow.eq.3)
     $       current_flow=glsc2_ms(vz,bm1ms,ntot1)/domain_length  ! for Z
      volvm1ms = glsum_ms(bm1ms,ntot1)

      if (iavflow.eq.1) then
         xsec = volvm1ms / domain_length
         flow_rate = param(55)*xsec
      endif

      delta_flow = flow_rate-current_flow

c     Note, this scale factor corresponds to FFX, provided FFX has
c     not also been specified in userf.   If ffx is also specified
c     in userf then the true FFX is given by ffx_userf + scale.

      scale = delta_flow/base_flow
      scale_vf(icvflow) = scale
      if (nio.eq.0) write(6,1) istep,chv(icvflow)
     $   ,time,scale,delta_flow,current_flow,flow_rate
    1    format(i10,'  volflow ',a1,11x,1p5e12.4)

      call add2s2(vx,vxc,scale,ntot1)
      call add2s2(vy,vyc,scale,ntot1)
      call add2s2(vz,vzc,scale,ntot1)
      call add2s2(pr,prc,scale,ntot2)

      return
      end
c-----------------------------------------------------------------------
      subroutine compute_vol_soln_ms(vxc,vyc,vzc,prc)
c
c     Compute the solution to the time-dependent Stokes problem
c     with unit forcing, and find associated flow rate.
c
c     pff 2/28/98
c
      include 'SIZE'
      include 'TOTAL'
c
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
c
      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)
      common /cvflow_i/ icvflow,iavflow
      common /cvflow_c/ chv(3)
      character*1 chv
 
      integer icalld
      save    icalld
      data    icalld/0/
 
      ntot1 = lx1*ly1*lz1*nelv
      if (icalld.eq.0) then
         icalld=icalld+1
         xlmin = glmin_ms(xm1,ntot1)
         xlmax = glmax_ms(xm1,ntot1)
         ylmin = glmin_ms(ym1,ntot1)          !  for Y!
         ylmax = glmax_ms(ym1,ntot1)
         zlmin = glmin_ms(zm1,ntot1)          !  for Z!
         zlmax = glmax_ms(zm1,ntot1)
 
         if (icvflow.eq.1) domain_length = xlmax - xlmin
         if (icvflow.eq.2) domain_length = ylmax - ylmin
         if (icvflow.eq.3) domain_length = zlmax - zlmin
      endif
 
      if (ifsplit) then
         call plan4_vol_ms(vxc,vyc,vzc,prc)
      else
         call plan3_vol_ms(vxc,vyc,vzc,prc)
      endif
c
c     Compute base flow rate
c 
      if (icvflow.eq.1)
     $       base_flow = glsc2_ms(vxc,bm1ms,ntot1)/domain_length
      if (icvflow.eq.2)
     $       base_flow = glsc2_ms(vyc,bm1ms,ntot1)/domain_length
      if (icvflow.eq.3)
     $       base_flow = glsc2_ms(vzc,bm1ms,ntot1)/domain_length
c
      if (nio.eq.0 .and. loglevel.gt.0) write(6,1) 
     $   istep,chv(icvflow),base_flow,domain_length,flow_rate
    1    format(i11,'  basflow ',a1,11x,1p3e13.4)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine plan3_vol_ms(vxc,vyc,vzc,prc)
c
c     Compute pressure and velocity using fractional step method.
c     (PLAN3).
c
      include 'SIZE'
      include 'TOTAL'
 
      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)
 
      COMMON /SCRNS/ rw1   (LX1,LY1,LZ1,LELV)
     $ ,             rw2   (LX1,LY1,LZ1,LELV)
     $ ,             rw3   (LX1,LY1,LZ1,LELV)
     $ ,             dv1   (LX1,LY1,LZ1,LELV)
     $ ,             dv2   (LX1,LY1,LZ1,LELV)
     $ ,             dv3   (LX1,LY1,LZ1,LELV)
     $ ,             RESPR (LX2,LY2,LZ2,LELV)
      COMMON /SCRVH/ H1    (LX1,LY1,LZ1,LELV)
     $ ,             H2    (LX1,LY1,LZ1,LELV)
      COMMON /SCRHI/ H2INV (LX1,LY1,LZ1,LELV)
      common /cvflow_i/ icvflow,iavflow

      real vxcbc(lx1,ly1,lz1,lelv)
      real vycbc(lx1,ly1,lz1,lelv)
      real vzcbc(lx1,ly1,lz1,lelv)
      real vxcp (lx1,ly1,lz1,lelv)
      real vycp (lx1,ly1,lz1,lelv)
      real vzcp (lx1,ly1,lz1,lelv)
      real resbc(lx1*ly1*lz1*lelv,ldim)

      common /cvflow_nn/ vxcbc,vycbc,vzcbc,vxcp,vycp,vzcp,resbc
 
 
c     Compute velocity, 1st part 
      n  = lx1*ly1*lz1*nelv
      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      ifield = 1
 
      call opzero(vxcbc,vycbc,vzcbc)
      call opzero(vxc,vyc,vzc)

      ngeompv = 20
      do ictr = 1,ngeompv
        if (icvflow.eq.1) then
         call copy     (rw1,bm1,ntot1)
         call rzero    (rw2,ntot1)
         call rzero    (rw3,ntot1)
        elseif (icvflow.eq.2) then
         call rzero    (rw1,ntot1)
         call copy     (rw2,bm1,ntot1)
         call rzero    (rw3,ntot1)
        else
         call rzero    (rw1,ntot1)        ! Z-flow!
         call rzero    (rw2,ntot1)        ! Z-flow!
         call copy     (rw3,bm1,ntot1)    ! Z-flow!
        endif

        if (ictr.eq.1) then
          intype = -1
          call sethlm   (h1,h2,intype)
          call ophinv   (vxc,vyc,vzc,rw1,rw2,rw3,h1,h2,tolhv,nmxh)
          call ssnormd  (vxc,vyc,vzc)
        else
          intype = -1
          call sethlm   (h1,h2,intype)

          call opcopy(vxcp,vycp,vzcp,vxc,vyc,vzc)

          call neknek_xfer_fld(vxc,vxcbc)
          call neknek_xfer_fld(vyc,vycbc)
          if (ldim.eq.3) call neknek_xfer_fld(vzc,vzcbc)
  
          call ophx(resbc(1,1),resbc(1,2),resbc(1,3),
     $             vxcbc,vycbc,vzcbc,h1,h2)

          call opsub2(rw1,rw2,rw3,resbc(1,1),resbc(1,2),resbc(1,3))
          call ophinv(vxc,vyc,vzc,rw1,rw2,rw3,h1,h2,tolhv,nmxh)
          call opadd2(vxc,vyc,vzc,vxcbc,vycbc,vzcbc)
          call ssnormd  (vxc,vyc,vzc)
        endif
c
c     Compute pressure  (from "incompr")
c
        intype = 1
        dtinv  = 1./dt
 
        call rzero   (h1,ntot1)
        call copy    (h2,vtrans(1,1,1,1,ifield),ntot1)
        call cmult   (h2,dtinv,ntot1)
        call invers2 (h2inv,h2,ntot1)
        call opdiv   (respr,vxc,vyc,vzc)
        call chsign  (respr,ntot2)
        call ortho   (respr)
 
 
c     Set istep=0 so that h1/h2 will be re-initialized in eprec
        i_tmp = istep
        istep = 0
        call esolver (respr,h1,h2,h2inv,intype)
        istep = i_tmp
 
        call opgradt (rw1,rw2,rw3,respr)
        call opbinv  (dv1,dv2,dv3,rw1,rw2,rw3,h2inv)
        call opadd2  (vxc,vyc,vzc,dv1,dv2,dv3)
 
        call cmult2  (prc,respr,bd(1),ntot2)

        call opsub2(vxcp,vycp,vzcp,vxc,vyc,vzc)
        dvxmax = glamax_ms(vxcp,ntot1)
        dvymax = glamax_ms(vycp,ntot1)
        dvzmax = glamax_ms(vzcp,ntot1)
         if (nio.eq.0)
     $      write(6,'(i2,i8,i4,1p4e13.4,a11)') idsess,istep,ictr,time,
     $      dvxmax,dvymax,dvzmax,' del-vol-vxy'
        call neknekgsync()
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine plan4_vol_ms(vxc,vyc,vzc,prc)

c     Compute pressure and velocity using fractional step method.
c     (Tombo splitting scheme).
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      real vxc(lx1,ly1,lz1,lelv)
     $   , vyc(lx1,ly1,lz1,lelv)
     $   , vzc(lx1,ly1,lz1,lelv)
     $   , prc(lx2,ly2,lz2,lelv)

      common /scrns/ resv1 (lx1,ly1,lz1,lelv)
     $ ,             resv2 (lx1,ly1,lz1,lelv)
     $ ,             resv3 (lx1,ly1,lz1,lelv)
     $ ,             respr (lx2*ly2*lz2,lelv)
     $ ,             TA1 (lx1*ly1*lz1*lelv)
     $ ,             TA2 (lx1*ly1*lz1*lelv)
     $ ,             TA3 (lx1*ly1*lz1*lelv)
     $ ,             WA1 (lx1*ly1*lz1*lelv)
     $ ,             WA2 (lx1*ly1*lz1*lelv)
     $ ,             WA3 (lx1*ly1*lz1*lelv)
      common /scrvh/ h1    (lx1,ly1,lz1,lelv)
     $ ,             h2    (lx1,ly1,lz1,lelv)
      COMMON /SCRMG/ W1    (LX1*LY1*LZ1,LELV)
     $ ,             W2    (LX1*LY1*LZ1,LELV)
     $ ,             W3    (LX1*LY1*LZ1,LELV)

      common /cvflow_i/ icvflow,iavflow

      real vxcbc(lx1,ly1,lz1,lelv)
      real vycbc(lx1,ly1,lz1,lelv)
      real vzcbc(lx1,ly1,lz1,lelv)
      real vxcp (lx1,ly1,lz1,lelv)
      real vycp (lx1,ly1,lz1,lelv)
      real vzcp (lx1,ly1,lz1,lelv)
      real resbc(lx1*ly1*lz1*lelv,ldim)

      common /cvflow_nn/ vxcbc,vycbc,vzcbc,vxcp,vycp,vzcp,resbc

      CHARACTER CB*3


      ntot1  = lx1*ly1*lz1*nelv
      ntot2  = lx2*ly2*lz2*nelv
      nxyz1  = lx1*ly1*lz1

      ngeompv = 20
      do ictr = 1,ngeompv
        call invers2  (h1,vtrans,ntot1)
        call rzero    (h2,       ntot1)
        if (ictr.eq.1) then
          call opzero(vxc,vyc,vzc)
        else 
          call neknek_xfer_fld(vxc,vxcbc)
          call neknek_xfer_fld(vyc,vycbc)
          if (ldim.eq.3) call neknek_xfer_fld(vzc,vzcbc)
        endif

c       Compute pressure 
        if (icvflow.eq.1) call cdtp(respr,h1,rxm2,sxm2,txm2,1)
        if (icvflow.eq.2) call cdtp(respr,h1,rym2,sym2,tym2,1)
        if (icvflow.eq.3) call cdtp(respr,h1,rzm2,szm2,tzm2,1)

        dtbd = BD(1)/DT

C     surface terms
        DO 100 IEL=1,NELV
          DO 300 IFC=1,2*ldim
            CALL RZERO  (W1(1,IEL),nxyz1)
            CALL RZERO  (W2(1,IEL),nxyz1)
            IF (ldim.EQ.3)
     $      CALL RZERO  (W3(1,IEL),nxyz1)
            CB = CBC(IFC,IEL,IFIELD)
            IF (intflag(ifc,iel).eq.1) then
               CALL FACCL3
     $         (W1(1,IEL),vxcbc(1,1,1,IEL),UNX(1,1,IFC,IEL),IFC)
               CALL FACCL3
     $         (W2(1,IEL),vycbc(1,1,1,IEL),UNY(1,1,IFC,IEL),IFC)
               IF (ldim.EQ.3)
     $          CALL FACCL3
     $         (W3(1,IEL),vzcbc(1,1,1,IEL),UNZ(1,1,IFC,IEL),IFC)
            ENDIF
            CALL ADD2   (W1(1,IEL),W2(1,IEL),nxyz1)
            IF (ldim.EQ.3)
     $      CALL ADD2   (W1(1,IEL),W3(1,IEL),nxyz1)
            CALL FACCL2 (W1(1,IEL),AREA(1,1,IFC,IEL),IFC)
            IF (intflag(ifc,iel).eq.1) then
              CALL CMULT(W1(1,IEL),dtbd,nxyz1)
            endif
            CALL SUB2 (RESPR(1,IEL),W1(1,IEL),nxyz1)
  300     CONTINUE
  100   CONTINUE


        call ortho    (respr)
        call ctolspl  (tolspl,respr)

        call hmholtz  ('PRES',prc,respr,h1,h2,pmask,vmult,
     $                             imesh,tolspl,nmxh,1)
        call ortho    (prc)

C     Compute velocity
        call opgrad   (resv1,resv2,resv3,prc)
        if (ifaxis) call col2 (resv2,omask,ntot1)
        call opchsgn  (resv1,resv2,resv3)

        if (icvflow.eq.1) call add2col2(resv1,v1mask,bm1,ntot1) ! add forcing
        if (icvflow.eq.2) call add2col2(resv2,v2mask,bm1,ntot1)
        if (icvflow.eq.3) call add2col2(resv3,v3mask,bm1,ntot1)
        if (ifexplvis) call split_vis ! split viscosity into exp/imp part

        call opcopy(vxcp,vycp,vzcp,vxc,vyc,vzc)

        intype = -1
        call sethlm(h1,h2,intype)
        call ophx(resbc(1,1),resbc(1,2),resbc(1,3),
     $             vxcbc,vycbc,vzcbc,h1,h2)
        call opsub2(resv1,resv2,resv3,resbc(1,1),resbc(1,2),resbc(1,3))
        call ophinv(vxc,vyc,vzc,resv1,resv2,resv3,h1,h2,tolhv,nmxh)
        call opadd2(vxc,vyc,vzc,vxcbc,vycbc,vzcbc)

        call opsub2(vxcp,vycp,vzcp,vxc,vyc,vzc)
        dvxmax = glamax_ms(vxcp,ntot1)
        dvymax = glamax_ms(vycp,ntot1)
        dvzmax = glamax_ms(vzcp,ntot1)
        if (ifexplvis) call redo_split_vis ! restore vdiff
        if (nid.eq.0)
     $    write(6,'(i2,i8,i4,1p4e13.4,a11)') idsess,istep,ictr,time,
     $    dvxmax,dvymax,dvzmax,' del-vol-vxy'

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_ms_dist()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      do ie=1,nelv
      do ifac=1,2*ldim
         if (intflag(ifac,ie).gt.0) cbc(ifac,ie,1) = 'int'
      enddo
      enddo

      call cheap_dist(distfint,1,'int')
      call dsavg(distfint)

      do ie=1,nelv
      do ifac=1,2*ldim
         if (intflag(ifac,ie).eq.1) cbc(ifac,ie,1) = 'v  '
         if (intflag(ifac,ie).eq.2) cbc(ifac,ie,1) = 'o  '
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine getupf()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'mpif.h'
      parameter (ltot=lx1*ly1*lz1*lelt)
      parameter (nxyz=lx1*ly1*lz1     )
      real    wtglls(nxyz,lelt,0:nsessmax-1),rsend(ltot*ldim),
     &        dist_all(ltot),rst_all(ltot*ldim),disti_all(ltot)
      integer rcode_all(ltot),elid_all(ltot),proc_all(ltot),
     &        rsid_nn(ltot)
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      npt = lx1*ly1*lz1*nelt
      do i=1,npt
        wtglls(i,1,idsess) = distfint(i,1,1,1)
      enddo

      do ids=1,nsessions-1
         idcheck = mod(idsess+ids,nsessions)
         call ifill(rsid_nn,idcheck,npt)
         call fgslib_findptsms(fpth_ms,rcode_all,1,proc_all,1,
     &                         elid_all,1,rst_all,ldim,dist_all,1,
     &                         xm1,1,ym1,1,zm1,1,rsid_nn,1,1,npt)
         call fgslib_findptsms_eval(fpth_ms,disti_all,1,
     &      rcode_all,1,proc_all,1,elid_all,1,rst_all,ldim,npt,distfint)

         do i=1,npt
           icd = rcode_all(i)
           dst = dist_all(i)
           idx = i
           if (icd.eq.2) then !this point not found
             wtglls(idx,1,idcheck) = 0.
           else
             wtglls(idx,1,idcheck) = disti_all(i)
           endif
         enddo
      enddo
      do i=1,npt
         prod = 0.
         do j=0,nsessions-1
           prod = prod+wtglls(i,1,j)
         enddo
         upval = wtglls(i,1,idsess)/prod
         upf(i,1,1,1) = upval
      enddo

      return
      end
c----------------------------------------------------------------------
