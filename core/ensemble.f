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

      call iniproc(intracomm)

      ifneknek   = .true.
      ifneknekm  = .false.


      return
      end
c-----------------------------------------------------------------------
      subroutine savg(ua,u,n,hndl)

c     Compute the average of quantity u() across sessions

      include 'SIZE'
      include 'TOTAL'
      include 'GLOBALCOM'

      real u,ua,weight
      integer hndl

      call copy(ua,u,n)
      call msgsop(ua,'sum',hndl)   ! Sum ua() across sessions

      weight = 1./nsessions
c      weight = 1.0

      call cmult(ua,weight,n)   ! Scale with 1./nsession for sess avg

      return
      end
c-----------------------------------------------------------------------
      subroutine msgsop_get_hndl(hndl,nel,lx,ly,lz)

c     Get a multi-session gsop handle

      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'

      integer e,eg

      common /c_is1/ glo_num(lx1*ly1*lz1*lelt)
      integer*8 glo_num
      integer hndl

      n = nel*lx*ly*lz

      do e=1,nel
         eg = lglel(e)
         do k=1,lz      
         do j=1,ly
         do i=1,lx
            ii = i + lx*(j-1) + lx*ly*(k-1) + lx*ly*lz*(e-1)
            glo_num(ii) = i + lx*(j-1) + lx*ly*(k-1) + 
     $                                   lx*ly*lz*(eg-1)
         enddo
         enddo
         enddo
      enddo


      call mpi_comm_size(mpi_comm_world,np_global,ierr)
      call fgslib_gs_setup(hndl,glo_num,n,mpi_comm_world,np_global)

      return
      end
c-----------------------------------------------------------------------
      subroutine msgsop(u,op,hndl)

c     multi-session version of gsop

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TSTEP'
      include  'CTIMER'

      real u(1)
      character*3 op

      if(ifsync) call nekgsync()

      if (op.eq.'+  ' .or. op.eq.'sum' .or. op.eq.'SUM')
     &   call fgslib_gs_op(hndl,u,1,1,0)

      if (op.eq.'*  ' .or. op.eq.'mul' .or. op.eq.'MUL')
     &   call fgslib_gs_op(hndl,u,1,2,0)

      if (op.eq.'m  ' .or. op.eq.'min' .or. op.eq.'mna' 
     &                .or. op.eq.'MIN' .or. op.eq.'MNA')
     &   call fgslib_gs_op(hndl,u,1,3,0)

      if (op.eq.'M  ' .or. op.eq.'max' .or. op.eq.'mxa'
     &                .or. op.eq.'MAX' .or. op.eq.'MXA')
     &   call fgslib_gs_op(hndl,u,1,4,0)

      return
      end
c------------------------------------------------------------------------
      subroutine ms_plan_tensr_op(ua,u,dir_hndl,ms_hndl,nelx,nely,nelz,
     $                            nel,nx,ifld,idx,op)
      include 'SIZE'
      include 'TOTAL'
      real u (lx1,ly1,lz1,lelt)
      real ua(lx1,ly1,lz1,lelt)
      real ub(lx1,ly1,lz1,lelt)
      integer dir_hndl,e,ex,ey,ez,eg,ms_hndl
      character*3 op

      n = nel*(nx**ldim)

      if (dir_hndl.eq.0.or.idx.lt.0) then
       idir = abs(idx)
       if (idir.eq.1) call set_gs_xavg_hndl(dir_hndl,nelx,nelyz,ifld)
       if (idir.eq.2) call set_gs_yavg_hndl(dir_hndl,nelx,nely,nelz,
     $                                                          ifld)
       if (idir.eq.3) call set_gs_zavg_hndl(dir_hndl,nelxy,ifld)
      endif
      if (ms_hndl.eq.0) then
        call msgsop_get_hndl(ms_hndl,nel,nx,ny,nz)
      endif

      if (op.eq.'ave') then
         call savg(ub,u,n,ms_hndl)
         call plan_tensr_avg(ua,ub,dir_hndl,ifld)
      elseif (op.eq.'*  ' .or. op.eq.'mul' .or. op.eq.'MUL') then
         call copy(ua,u,n)
         call msgsop(ua,op,ms_hndl)
         call fgslib_gs_op(dir_hndl,ua,1,2,0)
      elseif (op.eq.'m  ' .or. op.eq.'min' .or. op.eq.'mna'
     &   .or. op.eq.'MIN' .or. op.eq.'MNA') then
         call copy(ua,u,n)
         call msgsop(ua,op,ms_hndl)
         call fgslib_gs_op(dir_hndl,ua,1,3,0)
      elseif (op.eq.'M  ' .or. op.eq.'max' .or. op.eq.'mxa'
     &   .or. op.eq.'MAX' .or. op.eq.'MXA') then
         call copy(ua,u,n)
         call msgsop(ua,op,ms_hndl)
         call fgslib_gs_op(dir_hndl,ua,1,4,0)
      else
         if (nid.eq.0) write(6,*) 'Please enter a valid operation'
      endif


      return
      end
c------------------------------------------------------------------------
      subroutine multimesh_create

c     Dummy for ensemble

      return
      end
c------------------------------------------------------------------------
      subroutine userchk_set_xfer

c     Dummy for ensemble

      return
      end
c------------------------------------------------------------------------
      subroutine bcopy

c     Dummy for ensemble

      return
      end
C---------------------------------------------------------------------
      subroutine setintercomm(nekcommtrue,nptrue) 

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
      subroutine unsetintercomm(nekcommtrue,nptrue)

c     Dummy for ensemble

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
      subroutine happy_check(ihappy)

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_outflow ! Assign neighbor velocity to outflow if
                             ! characteristics are going the wrong way.

c     Dummy for ensemble

      return
      end
c-----------------------------------------------------------------------
