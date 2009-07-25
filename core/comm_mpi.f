c-----------------------------------------------------------------------
      subroutine iniproc
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal

      logical flag

      call mpi_initialized(mpi_is_initialized, ierr) !  Initialize MPI
      if ( mpi_is_initialized .eq. 0 ) then
         call mpi_init (ierr)
      endif

      ! create communicator
      call init_nek_comm
      np  = np_
      nid = nid_

      if(nid.eq.0) call printHeader

      ! check upper tag size limit
      call mpi_attr_get(MPI_COMM_WORLD,MPI_TAG_UB,nval,flag,ierr)
      if (nval.lt.(10000+max(lp,lelg))) then
         if(nid.eq.0) write(6,*) 'ABORT: MPI_TAG_UB too small!'
         call exitt
      endif

      IF (NP.GT.LP) THEN
         WRITE(6,*) 
     $   'ERROR: Code compiled for a max of',LP,' processors.'
         WRITE(6,*) 
     $   'Recompile with LP =',NP,' or run with fewer processors.'
         WRITE(6,*) 
     $   'Aborting in routine INIPROC.'
         call exitt
      endif

      ! check word size for REAL
      wdsize=4
      eps=1.0e-12
      oneeps = 1.0+eps
      if (oneeps.ne.1.0) then
         wdsize=8
      else
         if(nid.eq.0) 
     &     write(6,*) 'ABORT: single precision mode not supported!'
         call exitt
      endif
      nekreal = mpi_real
      if (wdsize.eq.8) nekreal = mpi_double_precision

      ifdblas = .false.
      if (wdsize.eq.8) ifdblas = .true.

      ! set word size for INTEGER
      ! HARDCODED since there is no secure way to detect an int overflow
      isize = 4
c
c
      PID = 0
      NULLPID=0
      NODE0=0
      NODE= NID+1

      if (nid.eq.0) then 
         write(6,*) 'Number of processors:',np
         WRITE(6,*) 'REAL    wdsize      :',WDSIZE
         WRITE(6,*) 'INTEGER wdsize      :',ISIZE
      endif

      call crystal_new(cr_h,nekcomm,np)  ! set cr handle to new instance

      return
      end
c-----------------------------------------------------------------------
      subroutine init_nek_comm
      include 'mpif.h'
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
C
      call nek_comm   ! set up nekton specific communicator
c
      nid_  = mynode()
      np_   = numnodes()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'CTIMER'
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
      real x(n), w(n)
      character*3 op
c
      if (ifsync) call gsync()

#ifndef NOTIMER
      if (icalld.eq.0) then
        tgop =0.0d0
        ngop =0
        icalld=1
      endif
      ngop = ngop + 1
      etime1=dnekclock()
#endif
c
      if (op.eq.'+  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
         call mpi_allreduce (x,w,n,nekreal,mpi_prod,nekcomm,ierr)
      else
         write(6,*) nid,' OP ',op,' not supported.  ABORT in GOP.'
         call exitt
      endif

      call copy(x,w,n)

#ifndef NOTIMER
      tgop =tgop +(dnekclock()-etime1)
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine igop( x, w, op, n)
c
c     Global vector commutative operation
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      common /ctimel/ ifsync
      logical ifsync

      integer x(n), w(n)
      character*3 op

      if (ifsync) call gsync()

      if     (op.eq.'+  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_sum ,nekcomm,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_max ,nekcomm,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_min ,nekcomm,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_allreduce (x,w,n,mpi_integer,mpi_prod,nekcomm,ierr)
      else
        write(6,*) nid,' OP ',op,' not supported.  ABORT in igop.'
        call exitt
      endif

      call icopy(x,w,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine csend(mtype,buf,len,jnid,jpid)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_send (buf,len,mpi_byte,jnid,mtype,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine crecv(mtype,buf,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $              ,jnid,mtype,nekcomm,status,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exitt
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine crecv3(mtype,buf,len,lenm)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
C
      real*4 buf(1)
      len = lenm
      jnid = mpi_any_source

      call mpi_recv (buf,len,mpi_byte
     $            ,jnid,mtype,nekcomm,status,ierr)
      call mpi_get_count (status,mpi_byte,len,ierr)
c
      if (len.gt.lenm) then 
          write(6,*) nid,'long message in mpi_crecv:',len,lenm
          call exit
      endif
c
      return
      end
c-----------------------------------------------------------------------
      integer function numnodes()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_comm_size (nekcomm, numnodes , ierr)

      return
      end
c-----------------------------------------------------------------------
      integer function mynode()
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer myid

      call mpi_comm_rank (nekcomm, myid, ierr)
      mynode = myid

      return
      end
c-----------------------------------------------------------------------
      real*8 function dnekclock()
      include 'mpif.h'
c
      dnekclock=mpi_wtime ()
c
      return
      end
c-----------------------------------------------------------------------
      subroutine lbcast(ifif)
C
C     Broadcast logical variable to all processors.
C
      include 'SIZE'
      include 'PARALLEL'
      include 'mpif.h'

      logical ifif

      if (np.eq.1) return

      item=0
      if (ifif) item=1
      call bcast(item,isize)
      ifif=.false.
      if (item.eq.1) ifif=.true.

      return
      end
c-----------------------------------------------------------------------
      subroutine bcast(buf,len)
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      real*4 buf(1)

      call mpi_bcast (buf,len,mpi_byte,0,nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine nek_comm
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_comm_group (mpi_comm_world,nekgroup,ierr)
      call mpi_comm_create (mpi_comm_world,nekgroup,nekcomm,ierr)
      call mpi_group_free (nekgroup,ierr)
c     write(6,*) 'nekcomm:',nekcomm

      return
      end
c-----------------------------------------------------------------------
      function isend(msgtag,x,len,jnid,jpid)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call mpi_isend (x,len,mpi_byte,jnid,msgtag
     $       ,nekcomm,imsg,ierr)
      isend = imsg
c     write(6,*) nid,' isend:',imsg,msgtag,len,jnid,(x(k),k=1,len/4)
c
      return
      end
c-----------------------------------------------------------------------
      function irecv(msgtag,x,len)
c
c     Note: len in bytes
c
      integer x(1)
C
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
C
      call mpi_irecv (x,len,mpi_byte,mpi_any_source,msgtag
     $       ,nekcomm,imsg,ierr)
      irecv = imsg
c     write(6,*) nid,' irecv:',imsg,msgtag,len
c
c
      return
      end
c-----------------------------------------------------------------------
      subroutine msgwait(imsg)
c
      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
c
c     write(6,*) nid,' msgwait:',imsg
c
      call mpi_wait (imsg,status,ierr)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gsync()

      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      call mpi_barrier(nekcomm,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine exitti(string,idata)
      character*1 string(80)
      character*11 s11
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'

      len = indx1(string,'$',1)
      write(s11,11) idata
   11 format(1x,i10)
      call chcopy(string(len),s11,11)

      if (nid.eq.0) write(6,1) (string(k),k=1,len+10)
    1 format('EXIT: ',80a1)

      call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine exitt
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'mpif.h'
c
      call gsync()

      tstop = dnekclock()
      ttotal= tstop-etimes

      if (nid.eq.0) then
         write(6,*) ' '
         write(6,'(A)') 'call exitt: dying ...'
         write(6,*) ' '
         call print_stack()
         write(6,*) ' '
         write(6,'(3(A,1g13.5,A,/)') 
     &      'total elapsed time         : ',ttotal, ' sec',
     &      'total solve time incl. I/O : ',ttime , ' sec',
     &      'time/timestep              : ',ttime/max(istep,1), ' sec'
      endif 
      call flush_io

      call mpi_finalize (ierr)
      call exit

c     z = -nx1
c     z = sqrt(z)
c     y = 1./(nx1-lx1)
c     y = 0.*y
c     a = 1./y
c     b = 1./y
c     write(6,*) 'quittin3',z,b
      call exit


      return
      end
c-----------------------------------------------------------------------

      subroutine printHeader

      INCLUDE 'HEADER'

      return
      end

c-----------------------------------------------------------------------
