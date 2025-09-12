c-----------------------------------------------------------------------
      subroutine usrdat2_oct(ncut) ! interface to oct-refine code

c     Note that lelt and lelg need to be LARGE enough

      include 'SIZE'
      include 'TOTAL'

      parameter(lxyz=lx1*ly1*lz1)
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
      subroutine set_interp_mat(ncut,pc,pt)
      include 'SIZE'
      include 'TOTAL'

      real pc(lx1*lx1,ncut) ! Interpolation matrix
      real pt(lx1*lx1,ncut)

      common /qtmp0/ z_out(lx1),wk(lx1*lx1),wk2(lx1*lx1)

      save ncut_save
      data ncut_save /0/

      if (ncut.le.1)
     $  call exitti('invalid ncut in set_interp_mat$',ncut)

      if (ncut.ne.ncut_save) then

        dr = 2./ncut
        do k=1,ncut
          r0 = -1. + (k-1)*dr
          do i=1,lx1
            z_out(i) = r0 + dr*(zgm1(i,1)+1)/2.
          enddo
          call interp_mat(pc(1,k),z_out,lx1,zgm1,lx1,wk,wk2)
          call transpose (pt(1,k),lx1,pc(1,k),lx1)
        enddo

        ncut_save = ncut

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

      common /qtmp0/ z_out(lx1),wk(lx1*lx1),wk2(lx1*lx1)

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

      common /qtmp0/ z_out(lx1),wk(lx1*lx1),wk2(lx1*lx1)

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

c     ADD MORE ITEMS, as needed

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
      subroutine nek_init_2
c
      include 'SIZE'
      include 'TOTAL'
      include 'DOMAIN'
c
      include 'OPCTR'
      include 'CTIMER'

      real kwave2
      logical ifemati

      real rtest
      integer itest
      integer*8 itest8
      character ctest
      logical ltest

      common /c_is1/ glo_num(lx1 * ly1 * lz1, lelt)
      common /ivrtx/ vertex((2 ** ldim) * lelt)
      integer*8 glo_num, ngv
      integer*8 vertex

      if (nio.eq.0) then
         write(6,12) 'nelgt/nelgv/lelt:',nelgt,nelgv,lelt
         write(6,12) 'lx1/lx2/lx3/lxd: ',lx1,lx2,lx3,lxd
 12      format(1X,A,4I12)
         write(6,*)
      endif

      instep=1             ! Check for zero steps
      if (nsteps.eq.0 .and. fintim.eq.0.) instep=0

      igeom = 2
      call setup_topo      ! Setup domain topology

      if (ifmvbd) call setup_mesh_dssum ! Set mesh dssum (needs geom)

      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine(glo_num,ncut)

c     Here we do an "ncut" refine:  ncut = 2 --> oct-refine (8x number of elements)
c                                   ncut = 3 --> 27x number of elements
c                                   ncut = 4 --> 64x number of elements

      include 'SIZE'
      include 'TOTAL'
      include 'SCRCT'  ! For xyz() array

      integer*8 glo_num(ncut+1,ncut+1,ncut*(ldim-2)+1,lelt)
      integer e,eg,egn,el,en,er,es,et

      parameter(lxyz=lx1*ly1*lz1,mxnew=500)
      common /qcrmg/ x0(lxyz,mxnew),y0(lxyz,mxnew),z0(lxyz,mxnew)
     $             , pc(lx1*lx1,mxnew),pt(lx1*lx1,mxnew)
      real pc,pt

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex
      integer*8 ngv

      integer ibuf(2)
      integer iwork(lelt)

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      save    isym2pre
      data    isym2pre / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      nblk = ncut**ldim
      nnew = nblk - 1
      lxyc = 2**ldim
      nvrt = ncut+1

      if (nio.eq.0) write(6,12) nblk
 12      format('h-refine: split each element into',I12)
c     CHECK limit sizes

      call lim_chk(nblk,mxnew,'nblk ','mxnew',' h_refine ')

      nelt_new = nblk*nelt+1 ! +1 for temporary storage
      call lim_chk(nelt_new,lelt,'n_new','lelt ',' h_refine ')

      nelgt_new = nblk*nelgt
      call lim_chk(nelgt_new,lelg,'ngnew','lelg ',' h_refine ')

      call set_interp_mat(ncut,pc,pt)

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

c         Note - we need to update xc,yc,zc
c         Note - we also should extract midside nodes.

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
          gllnid(eg) = nid
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
      subroutine h_refine_fld(u,nel,ncut)

c     Interpolate field onto the refined mesh 1

      include 'SIZE'

      real u(lx1,ly1,lz1,lelt)
      integer e,eg,egn,el,en,er,es,et

      parameter(lxyz=lx1*ly1*lz1,mxnew=500)
      common /qcrmg/ x0(lxyz,mxnew),y0(lxyz,mxnew),z0(lxyz,mxnew)
     $             , pc(lx1*lx1,mxnew),pt(lx1*lx1,mxnew)
      real pc,pt

      nblk = ncut**ldim

c     CHECK limit sizes

      call lim_chk(nblk,mxnew,'nblk ','mxnew',' h_refine_fld ')

      nel_new = nblk*nel+1 ! +1 for temporary storage
      call lim_chk(nel_new,lelt,'n_new','lelt ',' h_refine_fld ')

      nelg_new = iglsum(nblk*nel,1)
      call lim_chk(nelg_new,lelg,'ngnew','lelg ',' h_refine_fld ')

      call set_interp_mat(ncut,pc,pt)

      do e=nel,1,-1  ! REPLICATE EACH ELEMENT, working backward

        call get_rst_m1_fld(x0,ncut,u(1,1,1,e),pc,pt)

        en = nblk*(e - 1)

        do el=1,ncut**ldim

          en = en+1
          call copy(u(1,1,1,en),x0(1,el),lxyz)

        enddo

      enddo

      return
      end
c-----------------------------------------------------------------------
c     Input: global element id from coarse mesh
c     Output: id of the corresponding root-element in fine mesh
      integer function ie_map_c2f(eg,nblk)
      integer eg, nblk

c     nblk = ncut**ldim
      ie_map_c2f = (eg - 1) * nblk + 1
      return
      end
c-----------------------------------------------------------------------
c     In: global element id from fine mesh
c     Output: id of the corresponding root-element in coarse mesh
      integer function ie_map_f2c(egn,nblk)
      integer egn, nblk

      ie_map_f2c = (egn - 1) / nblk + 1
      return
      end
c-----------------------------------------------------------------------
      subroutine h_refine_copy(u,nel,ncut)
      include 'SIZE'
      real u(lx1,ly1,lz1,lelt)
      real ubak(lx1,ly1,lz1,lelt)
      integer ncut, nblk

      integer icalld
      save icalld
      data icalld / 0 /

      nblk = ncut**ldim
      nxyz = lx1*ly1*lz1
      call copy(ubak,u,nxyz*lelt)

      do ie=1,nel
        ien = ie_map_c2f(ie,nblk)
        call copy(u(1,1,1,ie),ubak(1,1,1,ien),nxyz)
      enddo

      icalld = icalld + 1

      return
      end
c-----------------------------------------------------------------------
      subroutine refine_map_elements(refine, refineSize)
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

      ierr = 0
      if (nid.eq.0) then
        if (nelgr*nblk_total.ne.nelgt) ierr = 1
      endif
      ierr = iglsum(ierr,1)
      if (ierr.ne.0)
     $  call exitti('nekf_refine_map nel mismstched$',nelgr)

      if (np.gt.1) then
      do iref=1,refineSize
        ncut = refine(iref)
        nblk = ncut**ldim
        do i=1,nelr                       ! go through elem in file
          iegr = er(i)                    ! global id of thie elem
          iegnr = ie_map_c2f(iegr,nblk)   ! map to the global id of the refined elem
          er(i) = iegnr                   ! put it back to the list
        enddo
      enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine refine_readfld(xm1_,ym1_,zm1_,vx_,vy_,vz_
     $                         ,pm1_,t_,ps_, refine, refineSize)
c      implicit none
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

      real xm1_(lx1,ly1,lz1,*), ym1_(lx1,ly1,lz1,*), zm1_(lx1,ly1,lz1,*)
      real vx_ (lx1,ly1,lz1,*), vy_ (lx1,ly1,lz1,*), vz_ (lx1,ly1,lz1,*)
      real pm1_(lx1,ly1,lz1,*)
      real t_  (lx1,ly1,lz1,*)
      real ps_ (lx1,ly1,lz1,lelt,*)

      integer refineSize
      integer refine(refineSize)

      ncut_total = 1
      do iref=1,refineSize
        ncut_total = ncut_total * refine(iref)
      enddo
      nblk_total = ncut_total**ldim

      if (np.gt.1) then
        nelt0 = nelt
        do iref=refineSize,1,-1 ! reverse copy
          ncut = refine(iref)
          nblk = ncut**ldim
          nelt0 = nelt0 / nblk

          if (ifgetxr) then
            call h_refine_copy(xm1_,nelt0,ncut)
            call h_refine_copy(ym1_,nelt0,ncut)
            call h_refine_copy(zm1_,nelt0,ncut)
          endif
          if (ifgetur) then
            call h_refine_copy(vx_,nelt0,ncut)
            call h_refine_copy(vy_,nelt0,ncut)
            call h_refine_copy(vz_,nelt0,ncut)
          endif
          if (ifgetpr) then
            call h_refine_copy(pm1_,nelt0,ncut)
          endif
          if (ifgettr) then
            call h_refine_copy(t_,nelt0,ncut)
          endif
          do k=1,npsr
            call h_refine_copy(ps_(1,1,1,1,k),nelt0,ncut)
          enddo
        enddo
      endif

      nelt0 = nelt / nblk_total
      do iref=1,refineSize
        ncut = refine(iref)
        nblk = ncut**ldim
        if (ifgetxr) then
          call h_refine_fld(xm1_,nelt0,ncut)
          call h_refine_fld(ym1_,nelt0,ncut)
          call h_refine_fld(zm1_,nelt0,ncut)
        endif
        if (ifgetur) then
          call h_refine_fld(vx_,nelt0,ncut)
          call h_refine_fld(vy_,nelt0,ncut)
          call h_refine_fld(vz_,nelt0,ncut)
        endif
        if (ifgetpr) then
          call h_refine_fld(pm1_,nelt0,ncut)
        endif
        if (ifgettr) then
          call h_refine_fld(t_,nelt0,ncut)
        endif
        do k=1,npsr
          call h_refine_fld(ps_(1,1,1,1,k),nelt0,ncut)
        enddo
        nelt0 = nelt0 * nblk
      enddo

      return
      end
c-----------------------------------------------------------------------
