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
C-----------------------------------------------------------------------
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
