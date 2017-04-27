c-----------------------------------------------------------------------
      subroutine get_session_info(intracomm)
      include 'mpif.h'
      include 'SIZE'
      include 'GLOBALCOM' 
      include 'TSTEP' 
      include 'INPUT' 
      
      common /happycallflag/ icall
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
      integer nid_global_root(0:nsessmax-1)
      character*132 session_mult(0:nsessmax-1), path_mult(0:nsessmax-1)

      logical ifhigh

C     nsessmax = upper limit for number of sessions
C     nfldmax_nn  = max number of fields to be interpolated
C     nmaxl_nn = max number of points at the boundary


C     Read from a file: number of sessions (nsessions), 
C     session name (SESSION_MULT(n-1)),
C     number of pocessors in each session (npsess(n-1)) 
C     and path (PATH_MULT(n-1))


      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
      if ( mpi_is_initialized .eq. 0 ) call mpi_init (ierr)

      call mpi_comm_size(mpi_comm_world,np_global,ierr)
      call mpi_comm_rank(mpi_comm_world,nid_global,ierr)

      nid=nid_global
      nekcomm=mpi_comm_world

      ierr = 0
      if (nid_global.eq.0) then
         open (unit=8,file='SESSION.NAME',status='old',err=24)
         read(8,*) nsessions
         do n=0,nsessions-1
            call blank(session_mult(n),132)
            call blank(path_mult(n)   ,132)
            read(8,10) session_mult(n)
            read(8,10) path_mult(n)
            read(8,*)  npsess(n)
         enddo
 10      format(a132)
         close(unit=8)
         goto 23
 24      ierr = 1
      endif
 23   continue
      
      call err_chk(ierr,' Cannot open SESSION.NAME!$')

      call bcast(nsessions,4)

      if (nsessions.gt.2) 
     &  call exitti('More than 2 sessions are currently 
     &  not supported!$',1)

      do n=0,nsessions-1
         call bcast(npsess(n),4)
         call bcast(session_mult(n),132)
         call bcast(path_mult(n),132)
      enddo

      npall=0
      do n=0,nsessions-1
         npall=npall+npsess(n)
      enddo
     
C     Check if number of processors in each session is consistent 
C     with the total number of processors

      if (npall.ne.np_global) 
     &  call exitti('Wrong number of processors!$',1)

C     Assign key for splitting into multiple groups

      nid_global_root_next=0
      do n=0,nsessions-1
         nid_global_root(n)=nid_global_root_next
         nid_global_root_next=nid_global_root(n)+npsess(n)
         if (nid_global.ge.nid_global_root(n).and.
     &   nid_global.lt.nid_global_root_next) 
     &    idsess=n
      enddo

      call mpi_comm_split(mpi_comm_world,idsess,nid,intracomm,ierr)
 
      session = session_mult(idsess)
      path    = path_mult   (idsess)

C     Intercommunications set up only for 2 sessions 

      if (nsessions.gt.1) then

         if (idsess.eq.0) idsess_neighbor=1
         if (idsess.eq.1) idsess_neighbor=0
 
         call mpi_intercomm_create(intracomm,0,mpi_comm_world, 
     &     nid_global_root(idsess_neighbor), 10,intercomm,ierr)

         np_neighbor=npsess(idsess_neighbor)
      
         call iniproc(intracomm)

         ifhigh=.true.
         call mpi_intercomm_merge(intercomm, ifhigh, iglobalcomm, ierr)
      
         ifneknek   = .true.
         ifneknekm  = .false.

         ninter = 1 ! Initialize NEKNEK interface extrapolation order to 1.

         icall = 0  ! Emergency exit call flag

      endif 

      return
      end
C-----------------------------------------------------------------------
      subroutine multimesh_create

      include 'SIZE'
      include 'TOTAL'
      real dxf,dyf,dzf

      integer icalld
      save    icalld
      data    icalld  /0/
c   Do some sanity checks - just once at setup
      call nekneksanchk(1)
C     Set interpolation flag: points with bc = 'int' get intflag=1. 
C     Boundary conditions are changed back to 'v' or 't'.

      if (icalld.eq.0) then
         call set_intflag
         call neknekmv()
         icalld = icalld + 1
      endif 
      call neknekgsync()

c   Figure out the displacement for the first mesh 
      call setup_int_neknek(dxf,dyf,dzf)  !sets up interpolation for 2 meshes

c    exchange_points2 finds the processor and element number at
c    comm_world level and displaces the 1st mesh back
      call exchange_points2(dxf,dyf,dzf)
      

      return
      end
C-------------------------------------------------------------
      subroutine set_intflag 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      character*3 cb
      character*2 cb2
      equivalence (cb2,cb)
      integer j,e,f

C     Set interpolation flag: points with boundary condition = 'int' 
c     get intflag=1. 
c
C     Boundary conditions are changed back to 'v' or 't'.

      ifield = 1
      if (ifheat) ifield = 2


      nfaces = 2*ndim
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
      subroutine userchk_set_xfer
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real l2,linf
      character*3 which_field(nfldmax_nn)

c     nfld_neknek is the number of fields to interpolate.
c     nfld_neknek = ndim+1 for velocities+pressure 
c     nfld_neknek = ndim+2 for velocities+pressure+temperature

      which_field(1)='vx'
      which_field(2)='vy'
      which_field(3)='vz'
      which_field(ndim+1)='pr'
      if (nfld_neknek.gt.ndim+1) which_field(ndim+2)='t'
c
c     Special conditions set for flow-poro coupling
      if(nfld_neknek.eq.1) then
         if(session.eq.'flow') which_field(1)='pr'
         if(session.eq.'poro') which_field(1)='vy'
      endif

      call neknekgsync()
      if (nsessions.gt.1) call get_values2(which_field)

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      n    = nx1*ny1*nz1*nelt

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
       c0=1
       c1=0
       c2=0
       else if (NINTER.eq.2.or.istep.eq.1) then
         c0=2
         c1=-1
         c2=0
       else 
         c0=3
         c1=-3
         c2=1
      endif
     
      do k=1,nfld_neknek
      do i=1,n
         ubc(i,1,1,1,k) = 
     $      c0*bdrylg(i,k,0)+c1*bdrylg(i,k,1)+c2*bdrylg(i,k,2)
      enddo
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine setintercomm(nekcommtrue,nptrue) 
      include 'SIZE' 
      include 'GLOBALCOM'
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal 
      
      nekcommtrue=nekcomm
      nekcomm=iglobalcomm

      nptrue=np
      np=np_global

      return
      end
c-----------------------------------------------------------------------
      subroutine unsetintercomm(nekcommtrue,nptrue)
      include 'SIZE' 
      include 'GLOBALCOM'
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal 

      nekcomm=nekcommtrue
      np=nptrue

      return
      end
c-----------------------------------------------------------------------
      function uglmin(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglmin=glmin(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c-----------------------------------------------------------------------
      function uglamax(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglamax=glamax(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c------------------------------------------------------------------------
      function uglmax(a,n)
      real a(1)

      call happy_check(1)
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      uglmax=glmax(a,n)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

      return
      end
c-----------------------------------------------------------------------
      subroutine neknekgsync()
      include 'SIZE' 
      include 'GLOBALCOM'

      call happy_check(1)
      call mpi_barrier(intercomm,ierr)
      return
      end
c------------------------------------------------------------------------
      subroutine happy_check(ihappy)
      include 'SIZE'
      include 'TOTAL'
      common /happycallflag/ icall

c     Happy check
      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      iglhappy=iglmin(ihappy,1)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original
      if (ihappy.eq.1.and.iglhappy.eq.0) then
         if (nid.eq.0) then
         write (6,*) '       '
         write (6,'(A,1i7,A,1e13.5)') 
     $   ' Emergency exit due to the other session:',
     $     ISTEP,'   time =',TIME
         write (6,*)   
         endif
         icall=1
       call exitt
      endif 

      return
      end
c------------------------------------------------------------------------
      subroutine chk_outflow_short ! Assign neighbor velocity to outflow

c
c     This is just an experimental routine for PnPn only...
c


      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'NEKNEK'
      integer e,eg,f

      n = nx1*ny1*nz1*nelt
      ipfld=ndim
c     if (ifsplit) ipfld=ndim+1
      ipfld=ndim+1                ! Now for both split and nonsplit (12/14/15)
      itfld=ipfld+1

      do i=1,n ! The below has not been checked for ifheat=.true., pff 6/27/15
         if (imask(i,1,1,1).eq.1) then
            vx(i,1,1,1) = ubc(i,1,1,1,1)
            vy(i,1,1,1) = ubc(i,1,1,1,2)
            if (if3d)    vz(i,1,1,1)  = ubc(i,1,1,1,3)
            if (ifsplit) pr(i,1,1,1)  = ubc(i,1,1,1,ipfld)
            if (ifheat)  t(i,1,1,1,1) = ubc(i,1,1,1,itfld)
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

      nface = 2*ndim
      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,1).eq.'o  ') then
           eg = lglel(e)
           call facind (i0,i1,j0,j1,k0,k1,nx1,ny1,nz1,f)
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

      call setintercomm(nekcommtrue,nptrue)    ! nekcomm=iglobalcomml
      iglmove=iglmin(imove,1)
      call unsetintercomm(nekcommtrue,nptrue)  ! nekcomm=nekcomm_original

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

ccccc
c     Get total number of processors and number of p
      npall = 0
      call neknekgsync()
      do i=1,nsessions
       npall = npall+npsess(i-1)
      enddo

ccccc
c     Get diamter of the domain
      call neknekgsync()
      mx_glob=uglmax(xm1,lx1*ly1*lz1*nelt)
      mn_glob=uglmin(xm1,lx1*ly1*lz1*nelt)
      dx1 = mx_glob-mn_glob
      call neknekgsync()

      dxf = 10.+dx1
      dyf = 0.
      dzf = 0.
ccccc
c     Displace MESH 1
      ntot = lx1*ly1*lz1*nelt
      if (idsess.eq.0) then
         call cadd(xm1,-dxf,ntot)
      endif

      call neknekgsync()
ccccc
c     Setup findpts    
      tol     = 1e-13
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by

      call findpts_setup(inth_multi2,mpi_comm_world,npall,ndim,
     &                   xm1,ym1,zm1,nx1,ny1,nz1,
     &                   nelt,nxf,nyf,nzf,bb_t,ntot,ntot,
     &                   npt_max,tol)

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

cccc
c     Look for boundary points with Diriclet b.c. (candidates for
c     interpolation)

      ifield = 1
      if (ifheat)  ifield   = 2

      nfaces = 2*ndim
      nel    = nelfld(ifield)

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nel
      nxy = lx1*ly1
      call izero(imask,ntot)

cccc
c     Setup arrays of x,y,zs to send to findpts and indices of boundary 
c     points in jsend
      ip = 0
      if (idsess.eq.0) then
        dxf = -dxf
      endif
      do e=1,nel
      do iface=1,nfaces
         if (intflag(iface,e).eq.1) then
            call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
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

cccc
c     JL's routine to find which points these procs are on
      call findpts(inth_multi2,rcode_all,1,
     &             proc_all,1,
     &             elid_all,1,
     &             rst_all,ndim,
     &             dist_all,1,
     &             rsend(1),ndim,
     &             rsend(2),ndim,
     &             rsend(3),ndim,nbp)

      call neknekgsync()
cccc
c     Move mesh 1 back to its original position
      if (idsess.eq.0) then
        call cadd(xm1,-dxf,lx1*ly1*lz1*nelt)
      endif

      ip=0
      icount=0
      ierror=0
cccc
c     Make sure rcode_all is fine
      do 200 i=1,nbp

      if (rcode_all(i).lt.2) then

        if (rcode_all(i).eq.1.and.dist_all(i).gt.1e-02) then
           if (ndim.eq.2) write(6,*)
     &     'WARNING: point on boundary or outside the mesh xy[z]d^2: '
           if (ndim.eq.3) write(6,*)
     &     'WARNING: point on boundary or outside the mesh xy[z]d^2: '
           goto 200
         endif
         ip=ip+1
         rcode(ip) = rcode_all(i)
         elid(ip)  = elid_all(i)
         proc(ip)  = proc_all(i)
         do j=1,ndim
           rst(ndim*(ip-1)+j)   = rst_all(ndim*(i-1)+j)
         enddo
         iList(1,ip) = jsend(i)

      endif  !  rcode_all

 200  continue

      npoints_nn = ip

      call iglmax(ierror,1)
      if (ierror.eq.1) call exitt
 
      call neknekgsync()

      return
      end
c-----------------------------------------------------------------------
      subroutine get_values2(which_field)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'

      parameter (lt=lx1*ly1*lz1*lelt,lxyz=lx1*ly1*lz1)
      common /scrcg/ pm1(lt),wk1(lxyz),wk2(lxyz)

      character*3 which_field(nfld_neknek)
      real fieldout(nmaxl_nn,nfldmax_nn)
      real field(lx1*ly1*lz1*lelt)
      integer nv,nt,i,j,k,n,ie,ix,iy,iz,idx,ifld

      call mappr(pm1,pr,wk1,wk2)  ! Map pressure to pm1 

      nv = nx1*ny1*nz1*nelv
      nt = nx1*ny1*nz1*nelt


cccc
c     Interpolate using findpts_eval
      do ifld=1,nfld_neknek
        if (which_field(ifld).eq.'vx') call copy(field,vx ,nt)
        if (which_field(ifld).eq.'vy') call copy(field,vy ,nt)
        if (which_field(ifld).eq.'vz') call copy(field,vz ,nt)
        if (which_field(ifld).eq.'pr') call copy(field,pm1,nt)
        if (which_field(ifld).eq.'t' ) call copy(field,t  ,nt)

        call field_eval(fieldout(1,ifld),1,field)
      enddo
      call neknekgsync()
         
cccc
c     Now we can transfer this information to valint array from which
c     the information will go to the boundary points
       do i=1,npoints_nn
        idx = iList(1,i)
        do ifld=1,nfld_neknek
          valint(idx,1,1,1,ifld)=fieldout(i,ifld)
        enddo
       enddo

      call neknekgsync()
      return
      end
C--------------------------------------------------------------------------
      subroutine field_eval(fieldout,fieldstride,fieldin)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKNEK'
      real fieldout(1)
      real fieldin(1)
      integer fieldstride
cccc
c     Used for findpts_eval of various fields
      call findpts_eval(inth_multi2,fieldout,fieldstride,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,npoints_nn,
     &                     fieldin)

      return
      end
C--------------------------------------------------------------------------
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

      if (nfld_neknek.lt.ndim) then
        if (nid.eq.0) write(6,*) 'Warning: Not all velocities are 
     $      being interpolated between sessions'
      endif       

      if (ifsplit) then
       if (nfld_neknek.lt.ndim+1) then
        if (nid.eq.0) write(6,*) 'Warning: Pressure is not being 
     $ interpolated.'
       endif       
      endif

      call izero(nintvh,2)
      ifield = 1
      if (ifheat) ifield = 2

      do i=1,ifield
      do e=1,nelt
      do f=1,2*ndim
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
      if (nfld_neknek.lt.ndim+2) then
      if (nid.eq.0) then
      write(6,*) 'Warning: Temperature might not be
     $ interpolated during the run. Check nfld_neknek in usrdat'
      if (nintv.eq.ninth) write(6,*) 'set nfld_neknek to atleast ndim+2'
      if (nintv.eq.ninth) call exitt
      endif
      endif
      endif

      call neknekgsync()

      if (nid.eq.0) write(6,*) idsess,nfld_neknek,nfldmax_nn,
     $  ifflow,ifheat,'Neknek log'
c     idsess - session number
c     nfld_neknek - fields to interpolate
      return
      end
C--------------------------------------------------------------------------
