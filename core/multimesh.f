c-----------------------------------------------------------------------
c     Routines for multidomain (neknek) simulations.
c
c     References:
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
      subroutine multimesh_create

      include 'SIZE'
      include 'TOTAL'
      real dxf,dyf,dzf

      integer icalld
      save    icalld
      data    icalld  /0/

      integer nfld_neknek
      common /inbc/ nfld_neknek

      call neknekgsync()
c     Do some sanity checks - just once at setup
c     Set interpolation flag: points with bc = 'int' get intflag=1. 
c     Boundary conditions are changed back to 'v' or 't'.

      if (icalld.eq.0) then
         nfld_neknek = ldim+nfield
         call nekneksanchk(1)
         call set_intflag
         call neknekmv()
         icalld = icalld + 1
      endif 
      call neknekgsync()

c     Figure out the displacement for the first mesh 
      call setup_int_neknek(dxf,dyf,dzf)  !sets up interpolation for 2 meshes

c     exchange_points2 finds the processor and element number at
c     comm_world level and displaces the 1st mesh back
      call exchange_points2(dxf,dyf,dzf)

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
c     get intflag=1. 
c
c     Boundary conditions are changed back to 'v' or 't'.

      ifield = 1
      if (ifheat) ifield = 2


      nfaces = 2*ldim
      nel    = nelfld(ifield)
      
      nflag=nel*nfaces
      call izero(intflag,nflag)

      do j=1,ifield
      do e=1,nel
      do f=1,nfaces
         cb=cbc(f,e,j)
         if (cb2.eq.'in') then
            intflag(f,e)=1
            if (j.eq.2) cbc(f,e,j)='t  '
            if (j.eq.1) cbc(f,e,j)='v  '
c           if (cb.eq.'inp') cbc(f,e,ifield)='on ' ! Pressure
c            if (cb.eq.'inp') cbc(f,e,ifield)='o  ' ! Pressure
            if (cb.eq.'inp') cbc(f,e,j)='o  ' ! Pressure
         endif
      enddo
      enddo
      enddo

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      integer k,i,n

      if (.not.ifneknekc) return

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
c     ngeom to ngeom=3-5 for scheme to be stable.


      if (NINTER.eq.1.or.istep.eq.0) then
       c0=1.
       c1=0.
       c2=0.
       else if (NINTER.eq.2.or.istep.eq.1) then
         c0=2.
         c1=-1.
         c2=0.
       else 
         c0=3.
         c1=-3.
         c2=1.
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
      subroutine chk_outflow_short ! Assign neighbor velocity to outflow
c
c     This is just an experimental routine for PnPn only...
c
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      integer e,eg,f

      n = lx1*ly1*lz1*nelt
      ipfld=ldim
c     if (ifsplit) ipfld=ldim+1
      ipfld=ldim+1                ! Now for both split and nonsplit (12/14/15)
      itfld=ipfld+1

      do i=1,n ! The below has not been checked for ifheat=.true., pff 6/27/15
         if (imask(i,1,1,1).eq.1) then
            vx(i,1,1,1) = valint(i,1,1,1,1)
            vy(i,1,1,1) = valint(i,1,1,1,2)
            if (if3d)    vz(i,1,1,1)  = valint(i,1,1,1,3)
            if (ifsplit) pr(i,1,1,1)  = valint(i,1,1,1,ipfld)
            if (ifheat)  t(i,1,1,1,1) = valint(i,1,1,1,itfld)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
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
      call neknekgsync()

      iglmove = ms_iglmin(imove,1)

      if (iglmove.eq.0) then
         ifneknekm=.true.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_int_neknek(dxf,dyf,dzf)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      include 'mpif.h'

      real dx1,dy1,dz1,dxf,dyf,dzf,mx_glob,mn_glob
      integer i,j,k,n,ntot2,npall
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
c     THIS ROUTINE DISPLACES THE FIRST MESH AND SETUPS THE FINDPTS
c     THE MESH IS DISPLACED BACK TO ORIGINAL POSITION IN EXCH_POINTS2

c     Get total number of processors and number of p
      npall = 0
      call neknekgsync()
      do i=1,nsessions
       npall = npall+npsess(i-1)
      enddo

c     Get diamter of the domain
      mx_glob = ms_glmax(xm1,lx1*ly1*lz1*nelt)
      mn_glob = ms_glmin(xm1,lx1*ly1*lz1*nelt)
      dx1 = mx_glob-mn_glob

      dxf = 10.+dx1
      dyf = 0.
      dzf = 0.

c     Displace MESH 1
      ntot = lx1*ly1*lz1*nelt
      if (idsess.eq.0) then
         call cadd(xm1,-dxf,ntot)
      endif

      call neknekgsync()

c     Setup findpts    
      tol     = 5e-13
      npt_max = 256
      nxf     = 2*lx1 ! fine mesh for bb-test
      nyf     = 2*ly1
      nzf     = 2*lz1
      bb_t    = 0.01 ! relative size to expand bounding boxes by

      if (istep.gt.1) call fgslib_findpts_free(inth_multi2)
      call fgslib_findpts_setup(inth_multi2,mpi_comm_world,npall,ldim,
     &                          xm1,ym1,zm1,lx1,ly1,lz1,
     &                          nelt,nxf,nyf,nzf,bb_t,ntot,ntot,
     &                          npt_max,tol)

      return
      end
c-----------------------------------------------------------------------
      subroutine exchange_points2(dxf,dyf,dzf)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      integer i,j,k,n
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer jsend(nmaxl_nn)
      common /exchr/ rsend(ldim*nmaxl_nn)
      integer rcode_all(nmaxl_nn),elid_all(nmaxl_nn)
      integer proc_all(nmaxl_nn)
      real    dist_all(nmaxl_nn)
      real    rst_all(nmaxl_nn*ldim)
      integer e,ip,iface,nel,nfaces,ix,iy,iz
      integer kx1,kx2,ky1,ky2,kz1,kz2,idx,nxyz,nxy
      integer icalld
      save    icalld
      data    icalld /0/

c     Look for boundary points with Diriclet b.c. (candidates for
c     interpolation)

      ifield = 1
      if (ifheat)  ifield   = 2

      nfaces = 2*ldim
      nel    = nelfld(ifield)

      nxyz  = lx1*ly1*lz1
      ntot  = nxyz*nel
      nxy = lx1*ly1
      call izero(imask,ntot)

c     Setup arrays of x,y,zs to send to findpts and indices of boundary 
c     points in jsend
      ip = 0
      if (idsess.eq.0) then
        dxf = -dxf
      endif
      do e=1,nel
      do iface=1,nfaces
         if (intflag(iface,e).eq.1) then
            call facind (kx1,kx2,ky1,ky2,kz1,kz2,lx1,ly1,lz1,iface)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               call nekasgn (ix,iy,iz,e)
               ip=ip+1
               idx = (e-1)*nxyz+(iz-1)*nxy+(iy-1)*lx1+ix
               jsend(ip) = idx 
               if (if3d) then
                 rsend(ldim*(ip-1)+1)=x-dxf
                 rsend(ldim*(ip-1)+2)=y
                 rsend(ldim*(ip-1)+3)=z
               else
                 rsend(ldim*(ip-1)+1)=x-dxf
                 rsend(ldim*(ip-1)+2)=y
               endif

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

      call neknekgsync()

c     JL's routine to find which points these procs are on
      call fgslib_findpts(inth_multi2,rcode_all,1,
     &                    proc_all,1,
     &                    elid_all,1,
     &                    rst_all,ldim,
     &                    dist_all,1,
     &                    rsend(1),ldim,
     &                    rsend(2),ldim,
     &                    rsend(3),ldim,nbp)

      call neknekgsync()

c     Move mesh 1 back to its original position
      if (idsess.eq.0) then
        call cadd(xm1,-dxf,lx1*ly1*lz1*nelt)
      endif

      ip=0
      icount=0
      ierror=0

c     Make sure rcode_all is fine
      do 200 i=1,nbp

      if (rcode_all(i).lt.2) then

        if (rcode_all(i).eq.1.and.dist_all(i).gt.1e-02) then
           if (ldim.eq.2) write(6,*)
     &     'WARNING: point on boundary or outside the mesh xy[z]d^2: '
           if (ldim.eq.3) write(6,*)
     &     'WARNING: point on boundary or outside the mesh xy[z]d^2: '
           goto 200
         endif
         ip=ip+1
         rcode(ip) = rcode_all(i)
         elid(ip)  = elid_all(i)
         proc(ip)  = proc_all(i)
         do j=1,ldim
           rst(ldim*(ip-1)+j)   = rst_all(ldim*(i-1)+j)
         enddo
         iList(1,ip) = jsend(i)

      endif  !  rcode_all

 200  continue

      ipg = iglsum(ip,1)
      nbpg = iglsum(nbp,1)
      if (nid.eq.0) write(6,*) 
     $      idsess,ipg,nbpg,'int pts bcs found findpt'
      npoints_nn = ip

c     zero out valint
      do i=1,nfld_neknek
       call rzero(valint(1,1,1,1,i),lx1*ly1*lz1*nelt)
      enddo

      ierror = ms_iglmax(ierror,1)
      if (ierror.eq.1) call exitt
 
      call neknekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine xfer_bcs_neknek
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      include 'CTIMER'

      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)
      integer nv,nt,i,j,k,n,ie,ix,iy,iz,idx,ifld

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
        do i=ldim+2,nfld_neknek  !do all passive scalars
          call field_eval(fieldout(1,i),1,t(1,1,1,1,i-ldim-1))
        enddo
      endif
         
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
       do i=1,npoints_nn
        idx = iList(1,i)
        do ifld=1,nfld_neknek
          valint(idx,1,1,1,ifld)=fieldout(i,ifld)
        enddo
       enddo

      call nekgsync()
      etime = dnekclock() - etime1
      tsync = etime1 - etime0

      if (nio.eq.0) write(6,99) istep,
     $              '  Multidomain data exchange done',
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

c     Used for findpts_eval of various fields
      call fgslib_findpts_eval(inth_multi2,fieldout,fieldstride,
     &                         rcode,1,
     &                         proc,1,
     &                         elid,1,
     &                         rst,ldim,npoints_nn,
     &                         fieldin)

      return
      end
c--------------------------------------------------------------------------
      subroutine nekneksanchk(flag)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      integer flag,nintvh(2),nintv,ninth,e,f,i,j,k
      character*3 cb
      character*2 cb2
      equivalence (cb2,cb)
c     Some sanity checks for neknek 

      if (nfld_neknek.gt.nfldmax_nn) then
        call exitti('Error: nfld_neknek > nfldmax:$',idsess)
      endif

      if (nfld_neknek.eq.0)
     $ call exitti('Error: set nfld_neknek in usrdat. Session:$',idsess)

      if (nfld_neknek.lt.ldim) then
        if (nid.eq.0) write(6,*) 'Warning: Not all velocities are 
     $      being interpolated between sessions'
      endif       

      if (ifsplit) then
       if (nfld_neknek.lt.ldim+1) then
        if (nid.eq.0) write(6,*) 'Warning: Pressure is not being 
     $ interpolated.'
       endif       
      endif

      call izero(nintvh,2)
      ifield = 1
      if (ifheat) ifield = 2

      do i=1,ifield
      do e=1,nelt
      do f=1,2*ldim
         cb=cbc(f,e,i)
         if (cb2.eq.'in') then
           nintvh(i) = nintvh(i)+1
         endif
      enddo
      enddo
      enddo

      nintv = iglsum(nintvh(1),1)
      ninth = iglsum(nintvh(2),1)

      if (ifield.eq.2) then
      if (nfld_neknek.lt.ldim+2) then
      if (nid.eq.0) then
      write(6,*) 'Warning: Temperature might not be
     $ interpolated during the run. Check nfld_neknek in usrdat'
      if (nintv.eq.ninth) write(6,*) 'set nfld_neknek to atleast ldim+2'
      if (nintv.eq.ninth) call exitt
      endif
      endif
      endif

      call neknekgsync()

      if (nid.eq.0) write(6,*) idsess,nfld_neknek,nfldmax_nn,
     $  ifflow,ifheat,'Neknek log'
c     idsess - session number
c     nfld_neknek - fields to interpolate
      if (nid.eq.0) write(6,*) ngeom,ninter,'Neknek log ngeom ninter' 
      return
      end
C--------------------------------------------------------------------------
      subroutine neknek_xfer_fld(u,ui)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real fieldout(nmaxl_nn,nfldmax_nn)
      real u(1),ui(1)
      integer nv,nt

cccc  Exchanges field u between the two neknek sessions and copies it 
cccc  to ui
c     Interpolate using findpts_eval
      call field_eval(fieldout(1,1),1,u)
cccc
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
       do i=1,npoints_nn
        idx = iList(1,i)
        ui(idx)=fieldout(i,1)
       enddo
      call neknekgsync()

      return
      end
C--------------------------------------------------------------------------
