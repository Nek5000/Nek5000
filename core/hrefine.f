c-----------------------------------------------------------------------
      subroutine h_refine_usrdat2(ncut) ! interface to oct-refine code
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      common /c_is1/ glo_num(lxyz*lelt)
      integer*8 glo_num

      ncut_o = ncut
      nelv_o = nelv
      nelt_o = nelt

      call h_refine(glo_num,ncut_o)

c     call h_refine_fld(vx,nelv_o,ncut_o)
c     call h_refine_fld(vy,nelv_o,ncut_o)
c     call h_refine_fld(vz,nelv_o,ncut_o)
c
c     call h_refine_fld(t,nelt_o,ncut_o)

c     call outpost(xm1,ym1,zm1,pr,t,'   ')

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine(glo_num,ncut)
c
c     Refine mesh including derivated vars
c                                   ncut = 1 --> do nothing
c     Here we do an "ncut" refine:  ncut = 2 --> oct-refine (8x number of elements)
c                                   ncut = 3 --> 27x number of elements
c                                   ncut = 4 --> 64x number of elements

      include 'SIZE'
      include 'TOTAL'
      include 'SCRCT'  ! For xyz() array

      integer*8 glo_num(ncut+1,ncut+1,ncut*(ldim-2)+1,lelt)
      integer e,eg,egn,el,en,er,es,et

      parameter(lxyz=lx1*ly1*lz1,mxmin=512,mxnew=max(mxmin,lelt))
      common /scrns/ x0(lxyz,mxnew),y0(lxyz,mxnew),z0(lxyz,mxnew)
     $             , pc(lx1*lx1,mxnew),pt(lx1*lx1,mxnew)
      real x0,y0,z0,pc,pt

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex
      integer*8 ngv

      integer ibuf(2)
      integer iwork(lelt)

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      save    isym2pre
      data    isym2pre / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      nvrt = ncut+1
      nblk = ncut**ldim
      call lim_chk(nblk,mxnew,'nblk ','mxnew',' h_refine ')

      if (ncut.lt.2) then
         if (nio.eq.0) write(6,11) ncut
 11      format('h-refine: ncut < 2, do nothing!  ncut=',I3)
         return
      endif

      if (ifaxis) then
         call exitti('h-refine does not support ifaxis=T$',ncut)
      endif

      nelt_new = nblk*nelt+1 ! +1 for temporary storage
      nelgt_new = nblk*nelgt
      call lim_chk(nelt_new,lelt,'n_new','lelt ',' h_refine ')
      call lim_chk(nelgt_new,lelg,'ngnew','lelg ',' h_refine ')

      if (nio.eq.0) write(6,12) nblk
 12      format('h-refine: split each element into',I12)

      call h_refine_set_interp_mat(lx1,ncut,pc,pt,.true.)

      call set_vert(glo_num,ngv,nvrt,nelt,vertex,.true.) ! Get new vertex set

      do e=nelt,1,-1  ! REPLICATE EACH ELEMENT, working backward

        eg = lglel(e)

        call elcopy(lelt,e) ! Save current element

        call get_rst_m1_vec(x0,y0,z0,ncut
     $           ,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),pc,pt)


        el = 0
        en = nblk*(e -1)
        egn= nblk*(eg-1)

        net=ncut
        if (ldim.eq.2) net=1
        do et=1,net          ! ncut^ldim reps
        do es=1,ncut
        do er=1,ncut

          el = el+1
          en = en+1
          egn= egn+1

          call elcopy(en,lelt) ! Copy saved element e to en

          call copy(xm1(1,1,1,en),x0(1,el),lxyz)
          call copy(ym1(1,1,1,en),y0(1,el),lxyz)
          call copy(zm1(1,1,1,en),z0(1,el),lxyz)

          k0=1 + (et-1)
          j0=1 + (es-1)
          i0=1 + (er-1)

          nk = 1
          if (ldim.eq.2) nk=0

          lc = 0
          do kc=0,nk
          do jc=0,1
          do ic=0,1
             lc=lc+1
             vertex(lc,en)=glo_num(i0+ic,j0+jc,k0+kc,e)

             ii = 1 + ic*(lx1-1)
             jj = 1 + jc*(ly1-1)
             kk = 1 + kc*(lz1-1)
             ipre=isym2pre(lc)
             xc(ipre,en) = xm1(ii,jj,kk,en)
             yc(ipre,en) = ym1(ii,jj,kk,en)
             zc(ipre,en) = zm1(ii,jj,kk,en)

             xyz(1,ipre,en)=xc(ipre,en)
             xyz(2,ipre,en)=yc(ipre,en)
             if (ldim.eq.3) xyz(3,ipre,en)=zc(ipre,en)

          enddo
          enddo
          enddo

          if (er.gt.1)    call fczero(4,en)  ! r-boundaries (face=4,2)
          if (er.lt.ncut) call fczero(2,en)
          if (es.gt.1)    call fczero(1,en)  ! s-boundaries (face=1,3)
          if (es.lt.ncut) call fczero(3,en)
          if (ldim.eq.3) then
             if (et.gt.1)    call fczero(5,en)  ! t-boundaries (face=5,6)
             if (et.lt.ncut) call fczero(6,en)
          endif

          ! remove curves as it's not updated correspondingly
          call rzero(curve(1,1,en),72)
          call blank(ccurve(1,en),12)

        enddo
        enddo
        enddo

      enddo

      ! Update elements count
      nelt0 = nelt

      nelt  = nblk*nelt
      nelv  = nblk*nelv
      nelgt = nblk*nelgt
      nelgv = nblk*nelgv
      do ifld=0,nfield
         nelfld(ifld)=nblk*nelfld(ifld)
      enddo

      ! Update lglel, gllel, gllnid
#if defined(DPROCMAP)
      call dProcMapClearCache()
#else
      call izero(gllnid,lelg)
      call izero(gllel,lelg)
#endif

      do e=nelt0,1,-1
        eg = lglel(e)
        en = nblk*(e -1)
        egn= nblk*(eg-1)

        do el=1,ncut**ldim
          en = en+1
          egn= egn+1

          lglel(en)  = egn
#if defined(DPROCMAP)
          ibuf(1) = en
          ibuf(2) = nid
          call dProcmapPut(ibuf,2,0,egn)
#else
          gllel(egn) = en
          gllnid(egn) = nid
#endif
        enddo
      enddo

#if !defined(DPROCMAP)
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllnid(k),iwork,'+  ',m)
         if (m.gt.0) call igop(gllel(k) ,iwork,'+  ',m)
         k = k+m
      enddo
#endif

      call nek_init_2  ! Reset some of the key arrays, etc.

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_mat(jmat,zo,no,zi,ni,s,t)
      real jmat(no,ni),zo(no),zi(ni),s(ni),t(ni)

      call rone(s,ni)
      do i=1,ni                ! Compute denominators ("alpha_i")
         do j=1,i-1
            s(i) = s(i)*(zi(i)-zi(j))
         enddo
         do j=i+1,ni
            s(i) = s(i)*(zi(i)-zi(j))
         enddo
      enddo

      do j=1,ni
         a_j = 1./s(j)
         do i=1,no
            jmat(i,j) = a_j
         enddo
      enddo

      call rone(s,ni)
      call rone(t,ni)
      do i=1,no
         x=zo(i)
         do j=2,ni
            jm   = ni+1-j
            s(j) = s(j-1 )*(x-zi(j-1))
            t(jm)= t(jm+1)*(x-zi(jm+1))
         enddo
         do j=1,ni
            jmat(i,j) = jmat(i,j)*s(j)*t(j)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_set_interp_mat(nx,ncut,pc,pt,ifrecomp)
c     interp mat for h-refine, GLL
c        nx: npts in 1 direction
c        ncut: new nel in 1 direction
      implicit none
      include 'SIZE'

      logical ifrecomp
      integer nx,ncut
      real pc(nx*nx,ncut) ! Interpolation matrix
      real pt(nx*nx,ncut) ! Transpose

      integer k, i, ncut_save, nx_save
      real z_out, z_hat, r0, dr,  wk, wk2
      common /qtmp0/ z_out(lx1),z_hat(lx1),wk(lx1*lx1),wk2(lx1*lx1)

      save ncut_save, nx_save
      data ncut_save /0/
      data nx_save /0/

      if (ncut.le.1) call exitti('invalid ncut in set_interp_mat$',ncut)
      if (nx.lt.2) call exitti('invalid nx in set_interp_mat$',nx)
      call lim_chk(nx,lx1,'nx   ','lx1  ',' set_intp ')

      if (ifrecomp) then ! recompute
         ncut_save = 0
         nx_save = 0
      endif

      if (ncut.ne.ncut_save.OR.nx.ne.nx_save) then

        call zwgll(z_hat,wk,nx)

        dr = 2./ncut
        do k=1,ncut
          r0 = -1. + (k-1)*dr
          do i=1,nx
            z_out(i) = r0 + dr*(z_hat(i)+1)/2.
          enddo
          call interp_mat(pc(1,k),z_out,nx,z_hat,nx,wk,wk2)
          call transpose (pt(1,k),nx,pc(1,k),nx)
        enddo

        ncut_save = ncut
        nx_save = nx

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_rst_m1_vec(xb,yb,zb,ncut,x1,y1,z1,pc,pt)

c       call get_rst_m1(x0,y0,z0,ncut
c    $                 ,xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),pc,pt)

      include 'SIZE'
      include 'TOTAL'

      real xb(lx1*ly1*lz1,ncut*ncut*ncut)
     $    ,yb(lx1*ly1*lz1,ncut*ncut*ncut)
     $    ,zb(lx1*ly1*lz1,ncut*ncut*ncut)

      real x1(lx1*ly1*lz1),y1(lx1*ly1*lz1),z1(lx1*ly1*lz1)
      real pc(lx1*lx1,ncut) ! Interpolation matrix
      real pt(lx1*lx1,ncut) ! Interpolation matrix

      common /qtmp0/ z_out(lx1),z_hat(lx1),wk(lx1*lx1),wk2(lx1*lx1)

      if (ldim.eq.3) then

       lc=0
       do kc=1,ncut
       do jc=1,ncut
       do ic=1,ncut
          lc=lc+1
          call tensr3(xb(1,lc),lx1,x1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
          call tensr3(yb(1,lc),lx1,y1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
          call tensr3(zb(1,lc),lx1,z1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
       enddo
       enddo
       enddo

      else ! 2D

       lc=0
       do jc=1,ncut
       do ic=1,ncut
          lc=lc+1
          call tensr3(xb(1,lc),lx1,x1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
          call tensr3(yb(1,lc),lx1,y1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
       enddo
       enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_rst_m1_fld(ub,ncut,u1,pc,pt)

      include 'SIZE'
      include 'TOTAL'

      real ub(lx1*ly1*lz1,ncut*ncut*ncut)

      real u1(lx1*ly1*lz1)
      real pc(lx1*lx1,ncut) ! Interpolation matrix
      real pt(lx1*lx1,ncut) ! Interpolation matrix

      common /qtmp0/ z_out(lx1),z_hat(lx1),wk(lx1*lx1),wk2(lx1*lx1)

      kcut=ncut
      if (ldim.eq.2) kcut=1

      lc=0
      do kc=1,kcut
      do jc=1,ncut
      do ic=1,ncut
         lc=lc+1
         call tensr3(ub(1,lc),lx1,u1,lx1,pc(1,ic),pt(1,jc),pt(1,kc),wk)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine elcopy(en,e)
      include 'SIZE'
      include 'TOTAL'

      integer e,en,eg

      lxyz = lx1*ly1*lz1

      call copy(xm1(1,1,1,en),xm1(1,1,1,e),lxyz)
      call copy(ym1(1,1,1,en),ym1(1,1,1,e),lxyz)
      call copy(zm1(1,1,1,en),zm1(1,1,1,e),lxyz)

      do kf=0,ldimt1
         call chcopy(cbc(1,en,kf),cbc(1,e,kf),18)
         call copy  (bc(1,1,en,kf),bc(1,1,e,kf),30)
      enddo
      call icopy (boundaryID(1,en),boundaryID(1,e),6)
      call icopy (boundaryIDt(1,en),boundaryIDt(1,e),6)

      call copy  (curve(1,1,en),curve(1,1,e),72)
      call chcopy(ccurve(1,en),ccurve(1,e),12)

      return
      end
c-----------------------------------------------------------------------
      subroutine fczero(f,e)
      include 'SIZE'
      include 'TOTAL'

      integer f,e

      do kf=0,nfield
          cbc(f,e,kf)='E  '
          call rzero(bc(1,f,e,kf),5)
      enddo
      boundaryID(f,e)=0
      boundaryIDt(f,e)=0

      return
      end
c-----------------------------------------------------------------------
c     Input: global element id from coarse mesh
c     Output: id of the first child  element in fine mesh
      integer function ie_map_o2r(eg,nblk)
      integer eg, nblk
      ie_map_o2r = (eg - 1) * nblk + 1
      return
      end
c-----------------------------------------------------------------------
c     Input: global element id from fine mesh
c     Output: id of the parents element in coarse mesh
      integer function ie_map_r2o(egn,nblk)
      integer egn, nblk
      ie_map_r2o = (egn - 1) / nblk + 1
      return
      end
c-----------------------------------------------------------------------
      subroutine nek_init_2
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1/lx2/lx3/lxd: ',lx1,lx2,lx3,lxd
 12      format(1X,A,4I12)
         write(6,*)
      endif

      igeom = 2
      call setup_topo      ! Setup domain topology

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      return
      end
c-----------------------------------------------------------------------
c     h refine + restart
c-----------------------------------------------------------------------
      subroutine h_refine_copy(u,nel,ncut)
      include 'SIZE'
      real u(lx1,ly1,lz1,lelt)
      real ubak(lx1,ly1,lz1,lelt)
      integer ncut, nblk

      nblk = ncut**ldim
      nxyz = lx1*ly1*lz1
      call rzero(ubak,nxyz*lelt)
      call copy(ubak,u,nxyz*lelt)
      call rzero(u,nxyz*lelt)

      do ie=1,nel
        ien = ie_map_o2r(ie,nblk)
        call copy(u(1,1,1,ie),ubak(1,1,1,ien),nxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_fld(u,nel,ncut)
c     apply one round of refinement to a field
      include 'SIZE'

      real u(lx1,ly1,lz1,lelt)
      integer e,eg,egn,el,en,er,es,et

      parameter(lxyz=lx1*ly1*lz1,mxmin=512,mxnew=max(mxmin,lelt))
      common /scrns/ x0(lxyz,mxnew),y0(lxyz,mxnew),z0(lxyz,mxnew)
     $             , pc(lx1*lx1,mxnew),pt(lx1*lx1,mxnew)
      real x0,y0,z0,pc,pt

      nblk = ncut**ldim
      call lim_chk(nblk,mxnew,'nblk ','mxnew',' h_refine_fld ')

      nel_new = nblk*nel+1 ! +1 for temporary storage
      call lim_chk(nel_new,lelt,'n_new','lelt ',' h_refine_fld ')

      nelg_new = iglsum(nblk*nel,1)
      call lim_chk(nelg_new,lelg,'ngnew','lelg ',' h_refine_fld ')

      do e=nel,1,-1  ! REPLICATE EACH ELEMENT, working backward

        call get_rst_m1_fld(x0,ncut,u(1,1,1,e),pc,pt)

        do el=1,ncut**ldim

          en = nblk*(e-1) + el
          call copy(u(1,1,1,en),x0(1,el),lxyz)

        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_remap_elem(refine, refineSize)
c     restart, modify er to re-distribute elements
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      integer refineSize
      integer refine(refineSize)

      ! check number of elements
      ncut_total = 1
      do iref=1,refineSize
        ncut_total = ncut_total * refine(iref)
      enddo
      nblk_total = ncut_total**ldim
      if (nio.eq.0) write(*,31) ncut_total,nblk_total
  31  format(3x,'mfi:href rs map_e ncut/nblk:',2(1I8))
      if (refineSize.eq.0.OR.ncut_total.lt.2) return

      ierr = 0
      if (nid.eq.0) then
        if (nelgr*nblk_total.ne.nelgt) ierr = 1
      endif
      ierr = iglsum(ierr,1)
      if (ierr.ne.0)
     $  call exitti('h_refine_map nel mismstched$',nelgr)

      if (np.gt.1) then
      do iref=1,refineSize
        ncut = refine(iref)
        nblk = ncut**ldim
        do i=1,nelr                       ! go through elem in file
          iegr = er(i)                    ! global id of this elem
          iegnr = ie_map_o2r(iegr,nblk)   ! map to the global id of the refined elem
          er(i) = iegnr                   ! put it back to the list
        enddo
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_readfld(xm1_,ym1_,zm1_,vx_,vy_,vz_
     $                           ,pm1_,t_,ps_, refine, refineSize)
c     restart, refine fields after readfld
      implicit none
      include 'SIZE'
      include 'INPUT' ! ifaxis
      include 'PARALLEL' ! np
      include 'RESTART'

      real xm1_(lx1,ly1,lz1,*), ym1_(lx1,ly1,lz1,*), zm1_(lx1,ly1,lz1,*)
      real vx_ (lx1,ly1,lz1,*), vy_ (lx1,ly1,lz1,*), vz_ (lx1,ly1,lz1,*)
      real pm1_(lx1,ly1,lz1,*)
      real t_  (lx1,ly1,lz1,*)
      real ps_ (lx1,ly1,lz1,lelt,*)

      integer refineSize
      integer refine(refineSize)
      integer iref,k,ncut,nblk,ncut_total,nblk_total,nelt0

      integer lxyz,mxmin,mxnew
      parameter(lxyz=lx1*ly1*lz1,mxmin=512,mxnew=max(mxmin,lelt))
      common /scrns/ x0(lxyz,mxnew),y0(lxyz,mxnew),z0(lxyz,mxnew)
     $             , pc(lx1*lx1,mxnew),pt(lx1*lx1,mxnew)
      real x0,y0,z0,pc,pt

      ncut_total = 1
      do iref=1,refineSize
        ncut_total = ncut_total * refine(iref)
      enddo
      nblk_total = ncut_total**ldim
      if (nio.eq.0) write(*,31) ncut_total,nblk_total
  31  format(3x,'mfi:href rs ref_e ncut/nblk:',2(1I8))
      if (refineSize.eq.0.OR.ncut_total.lt.2) return

      if (ifaxis) then
         call exitti('h-refine does not support ifaxis=T$',ncut)
      endif

      if (np.gt.1) then
        nelt0 = nelt
        do iref=refineSize,1,-1 ! reverse copy
          ncut = refine(iref)
          nblk = ncut**ldim
          nelt0 = nelt0 / nblk

          if (ifgetxr.AND.ifgetx) then
            call h_refine_copy(xm1_,nelt0,ncut)
            call h_refine_copy(ym1_,nelt0,ncut)
            call h_refine_copy(zm1_,nelt0,ncut)
          endif
          if (ifgetur.AND.ifgetu) then
            call h_refine_copy(vx_,nelt0,ncut)
            call h_refine_copy(vy_,nelt0,ncut)
            call h_refine_copy(vz_,nelt0,ncut)
          endif
          if (ifgetpr.AND.ifgetp) then
            call h_refine_copy(pm1_,nelt0,ncut)
          endif
          if (ifgettr.AND.ifgett) then
            call h_refine_copy(t_,nelt0,ncut)
          endif
          do k=1,npsr
            if (ifgtpsr(k).AND.ifgtps(k))then
              call h_refine_copy(ps_(1,1,1,1,k),nelt0,ncut)
            endif
          enddo
        enddo
      endif

      nelt0 = nelt / nblk_total
      do iref=1,refineSize
        ncut = refine(iref)
        nblk = ncut**ldim

        call h_refine_set_interp_mat(lx1,ncut,pc,pt,.true.)

        if (ifgetxr.AND.ifgetx) then
          call h_refine_fld(xm1_,nelt0,ncut)
          call h_refine_fld(ym1_,nelt0,ncut)
          call h_refine_fld(zm1_,nelt0,ncut)
        endif
        if (ifgetur.AND.ifgetu) then
          call h_refine_fld(vx_,nelt0,ncut)
          call h_refine_fld(vy_,nelt0,ncut)
          call h_refine_fld(vz_,nelt0,ncut)
        endif
        if (ifgetpr.AND.ifgetp) then
          call h_refine_fld(pm1_,nelt0,ncut)
        endif
        if (ifgettr.AND.ifgett) then
          call h_refine_fld(t_,nelt0,ncut)
        endif
        do k=1,npsr
          if (ifgtpsr(k).AND.ifgtps(k))then
            call h_refine_fld(ps_(1,1,1,1,k),nelt0,ncut)
          endif
        enddo
        nelt0 = nelt0 * nblk
      enddo

      return
      end
c-----------------------------------------------------------------------
c     Extra interface to recover original mesh info, for hMG
c-----------------------------------------------------------------------
      subroutine h_refine_r2o_nel(nelv_o,nelt_o,ncut)
      implicit none
      include 'SIZE'
      integer nelv_o,nelt_o,ncut,nblk

      nblk = ncut**ldim

      if (mod(nelv,nblk).ne.0)
     $   call exitti('ref_r2o nelv not divisible$',nblk)
      if (mod(nelt,nblk).ne.0)
     $   call exitti('ref_r2o nelt not divisible$',nblk)

      nelv_o = nelv / nblk
      nelt_o = nelt / nblk

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_r2o_vertex(vtxo,vtxr,nelo,ncut)
      implicit none
      include 'SIZE'
      integer*8 vtxo(2**ldim,1), vtxr(2**ldim,1)
      integer nelo, ncut, nblk
      integer e, en, e1,e2,e3,e4,e5,e6,e7,e8

      nblk = ncut**ldim

      e1 = 1
      e2 = ncut
      e3 = (ncut-1)*ncut+1
      e4 = ncut*ncut
      e5 = e1 + ncut*ncut
      e6 = e2 + ncut*ncut
      e7 = e3 + ncut*ncut
      e8 = e4 + ncut*ncut

      do e=1,nelo
         en = (e-1) * nblk
         vtxo(1,e) = vtxr(1,en+e1)
         vtxo(2,e) = vtxr(2,en+e2)
         vtxo(3,e) = vtxr(3,en+e3)
         vtxo(4,e) = vtxr(4,en+e4)

         if (ldim.eq.3) then
            vtxo(5,e) = vtxr(5,en+e5)
            vtxo(6,e) = vtxr(6,en+e6)
            vtxo(7,e) = vtxr(7,en+e7)
            vtxo(8,e) = vtxr(8,en+e8)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_r2o_cbc(CBCo,CBCr,nelo,ncut)
      implicit none
      include 'SIZE'
      integer e,el,er,nelo,ncut,kcut,nblk
      integer ic,jc,kc
      character*3 CBCo(6,1), CBCr(6,1)

      nblk = ncut**ldim

      kcut = ncut
      if (ldim.eq.2) kcut = 1

      do e=1,nelo

         el = 0
         do kc=1,kcut
         do jc=1,ncut
         do ic=1,ncut
            el = el + 1
            er = (e-1) * nblk + el

            if (ic.eq.1)         CBCo(4,e) = CBCr(4,er)
            if (ic.eq.ncut)      CBCo(2,e) = CBCr(2,er)
            if (jc.eq.1)         CBCo(1,e) = CBCr(1,er)
            if (jc.eq.ncut)      CBCo(3,e) = CBCr(3,er)
            if (ldim.eq.3) then
               if (kc.eq.1)      CBCo(5,e) = CBCr(5,er)
               if (kc.eq.ncut)   CBCo(6,e) = CBCr(6,er)
            endif
         enddo
         enddo
         enddo

      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine hrefcuts_i2c(cout) ! convert integer to base 62 alphabet
      implicit none
      include 'SIZE'
      include 'INPUT'

      character*4 cout, ctmp
      character*1 c4(4)
      equivalence (c4,ctmp)

c     Base 62 alphabet
      character*62 B62a
      character*1  B62b(62)
      equivalence (B62a,B62b)
      save B62a

      integer i, nround, ncut, ierr

      B62a = '0123456789'
     &    // 'abcdefghijklmnopqrstuvwxyz'
     &    // 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

      ierr = 0
      cout = '0000'
      if (nhref.eq.0) return

      call blank(ctmp,4)

      nround = min(nhref,4)
      do i=1,nround
         ncut = hrefcuts(i)
         if (ncut.gt.61.OR.ncut.lt.1) then
            if (nio.eq.0) write(6,20) ncut
            ierr = 1
         else
            c4(i) = B62b(ncut+1)
         endif
      enddo

      if (nhref.gt.4) then
         if (nio.eq.0) write(6,21) nhref
         ierr = 2
      endif

      cout = ctmp

   20 format('WARN: href i2c only support ncut in [1,61] ncut=',i3)
   21 format('WARN: href i2c only support upto 4 rounds nhref=',i3)

      if (ierr.ne.0) then
        if(nio.eq.0)write(*,*) 'WARN: href i2c fallback to empty string'
        cout = '    '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hrefcuts_c2i(cin) ! convert string to int list
      implicit none
      include 'SIZE'
      include 'RESTART'

      character*4 cin
      character*4 cout, ctmp
      character*1 c4(4)
      equivalence (c4,ctmp)

c     Base 61 alphabet
      character*62 B62
      parameter (B62 = '0123456789'
     &              // 'abcdefghijklmnopqrstuvwxyz'
     &              // 'ABCDEFGHIJKLMNOPQRSTUVWXYZ')

      integer i, nround, ncut, pos, ierr

      if (cin.eq.'    ') return ! do nothing for legacy format to allow overwrite

      nhrefrs = 0
      call izero(hrefcutsrs, lhref)
      if (cin.eq.'0000') return ! empty schedule

      ierr = 0
      call chcopy(ctmp,cin,4)

c     convert all characters into integers
      nround = min(4, lhref)
      do i=1,nround
         pos = index(B62, c4(i))
         if (c4(i).eq.' ') then
            hrefcutsrs(i) = 0
         elseif (pos.gt.0) then
            hrefcutsrs(i) = pos - 1
         else
            hrefcutsrs(i) = 0
            goto 90 ! err
         endif
      enddo

c     use first zero to determine the length
      nhrefrs = nround
      do i=1,nround
         if (hrefcutsrs(i).eq.0) then
            nhrefrs = i - 1
            goto 40
         endif
      enddo
   40 continue

c     zero out the rest
      do i=1,lhref
         if (i.gt.nhrefrs) hrefcutsrs(i) = 0
      enddo
      return

c     zero everything when error
   90 continue
      if (nio.eq.0)
     $  write(*,*)'WARN: href c2i invalid fmt, zero hrefcutsrs|',cin,'|'

      nhrefrs = 0
      call izero(hrefcutsrs, lhref)

      return
      end
c-----------------------------------------------------------------------
      subroutine hrefcuts_chkdiff
c
c     Input:
c        hrefcuts:      h-refine schedule from INPUT, e.g. par
c        hrefcutsrs:    h-refine schedule from RESTART, e.g., checkpoint hdr
c     This subroutine return the ordered diffeference
c        hrefcutsrs = hrefcuts \ hrefcutsrs
c     which is the extra refinement on the top of checkpoint to match simulation
c
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'RESTART'
      include 'PARALLEL' ! nelgt

      integer nblk, nblk_rs, ncut, ncut_rs, i, j, ierr
      integer nelgr0, nelgt0

      if (nhref.eq.0) return

      if (nio.eq.0) then
         write(*,*)'href schdule, sim: ', (hrefcuts(i),i=1,nhref)
         write(*,*)'href schdule, fld: ', (hrefcutsrs(i),i=1,nhrefrs)
      endif

      ncut_rs = 1
      do i=1,nhrefrs
         ncut_rs = ncut_rs * hrefcutsrs(i)
      enddo
      nblk_rs = ncut_rs**ldim

      ncut = 1
      do i=1,nhref
         ncut = ncut * hrefcuts(i)
      enddo
      nblk = ncut**ldim

      nelgr0 = nelgr / nblk_rs
      nelgt0 = nelgt / nblk

      if (nelgr0.ne.nelgt0)
     $  call exitti('href schdule diff initial mesh mismatched$',nelgr0)

      if (nhrefrs.gt.nhref) then
         ierr = nhrefrs
         goto 90
      endif

      ierr = 0
      do i=1,nhrefrs
         if (hrefcuts(i).ne.hrefcutsrs(i)) then
            ierr = i
            goto 90
         endif
      enddo

      ! set new schedule to be applied
      call izero(hrefcutsrs, lhref)
      j = 0
      do i = nhrefrs+1,nhref
         j = j + 1
         hrefcutsrs(j) = hrefcuts(i)
      enddo
      nhrefrs = j;

      if (nhrefrs.GT.lhref) call exitti('nhref rs > lhref$',nhrefrs)

      if (nio.eq.0) then
         write(*,*)'href schdule, dif: ', (hrefcutsrs(i),i=1,nhrefrs)
      endif

      return

   90 continue
      call exitti('href rs schedule is not a subset$',ierr)

      return
      end
c-----------------------------------------------------------------------
