c*********************************************************************72
      subroutine mpi_scan(data1, data2, n, datatype,
     &  operation, comm, ierror )

      data2 = data1

      return
      end

c*********************************************************************72
      subroutine mpi_abort ( comm, errorcode, ierror )

c*********************************************************************72
c
cc MPI_ABORT shuts down the processes in a given communicator.
c
      implicit none

      integer comm
      integer errorcode
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_ABORT:'
      write ( *, '(a,i12)' ) 
     &  '  Shut down with error code = ', errorcode

      stop
      end
      subroutine mpi_allgather ( data1, nsend, sendtype, data2, 
     &  nrecv, recvtype, comm, ierror )

c*********************************************************************72
c
cc MPI_ALLGATHER gathers data from all the processes in a communicator.
c
      implicit none

      include "mpi_dummy.h"

      integer nsend

      integer comm
      integer data1(nsend)
      integer data2(nsend)
      integer ierror
      integer nrecv
      integer recvtype
      integer sendtype

      ierror = MPI_SUCCESS

      if ( sendtype .eq. mpi_double_precision ) then
        call mpi_copy_double_precision ( data1, data2, nsend, ierror )
      else if ( sendtype .eq. mpi_integer ) then
        call mpi_copy_integer ( data1, data2, nsend, ierror )
      else if ( sendtype .eq. mpi_real ) then
        call mpi_copy_real ( data1, data2, nsend, ierror )
      else
        ierror = MPI_FAILURE
      end if

      return
      end
      subroutine mpi_allgatherv ( data1, nsend, sendtype,
     &  data2, nrecv, ndispls, recvtype, comm, ierror )

c*********************************************************************72
c
cc MPI_ALLGATHERV gathers data from all the processes in a communicator.
c
      implicit none

      include "mpi_dummy.h"

      integer nsend

      integer comm
      integer data1(nsend)
      integer data2(nsend)
      integer ierror
      integer ndispls
      integer nrecv
      integer recvtype
      integer sendtype

      ierror = MPI_SUCCESS

      if ( sendtype .eq. mpi_double_precision ) then
        call mpi_copy_double_precision ( data1, data2, nsend, ierror )
      else if ( sendtype .eq. mpi_integer ) then
        call mpi_copy_integer ( data1, data2, nsend, ierror )
      else if ( sendtype .eq. mpi_real ) then
        call mpi_copy_real ( data1, data2, nsend, ierror )
      else
        ierror = MPI_FAILURE
      end if

      return
      end
      subroutine mpi_allreduce ( data1, data2, n, datatype,
     &  operation, comm, ierror )

c*********************************************************************72
c
cc MPI_ALLREDUCE carries out a reduction operation.
c
      implicit none

      include "mpi_dummy.h"

      integer n

      integer comm
      integer data1(n)
      integer data2(n)
      integer datatype
      integer ierror
      integer operation

      ierror = MPI_SUCCESS

      if ( datatype .eq. mpi_double_precision ) then

        call mpi_reduce_double_precision ( 
     &    data1, data2, n, operation, ierror )

      else if ( datatype .eq. mpi_integer ) then

        call mpi_reduce_integer ( 
     &    data1, data2, n, operation, ierror )

      else if ( datatype .eq. mpi_real ) then

        call mpi_reduce_real ( 
     &    data1, data2, n, operation, ierror )

      else

        ierror = MPI_FAILURE

      end if

      return
      end
      subroutine mpi_barrier ( comm, ierror )

c*********************************************************************72
c
cc MPI_BARRIER forces processes within a communicator to wait together.
c
      implicit none

      integer comm
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      return
      end
      subroutine mpi_bcast ( data, n, datatype, node, comm, ierror )

c*********************************************************************72
c
cc MPI_BCAST broadcasts data from one process to all others.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )
      integer node

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_bsend ( data, n, datatype, iproc, itag,
     &  comm, ierror )

c*********************************************************************72
c
cc MPI_BSEND sends data from one process to another, using buffering.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_BSEND - Error!'
      write ( *, '(a)' )  '  Should not send message to self.'

      return
      end
      subroutine mpi_cart_create ( comm, ndims, dims, periods,
     &  reorder, comm_cart, ierror )

c*********************************************************************72
c
cc MPI_CART_CREATE creates a communicator for a Cartesian topology.
c
      implicit none

      integer ndims

      integer comm
      integer comm_cart
      integer dims(*)
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )
      logical periods(*)
      logical reorder

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_cart_get ( comm, ndims, dims, periods,
     &  coords, ierror )

c*********************************************************************72
c
cc MPI_CART_GET returns the "Cartesian coordinates" of the calling process.
c
      implicit none

      integer ndims

      integer comm
      integer coords(*)
      integer dims(*)
      integer i
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )
      logical periods(*)

      ierror = MPI_SUCCESS

      do i = 1, ndims
        coords(i) = 0
      end do

      return
      end
      subroutine mpi_cart_shift ( comm, idir, idisp, isource, 
     &  idest, ierror )

c*********************************************************************72
c
cc MPI_CART_SHIFT finds the destination and source for Cartesian shifts.
c
      implicit none

      integer comm
      integer idest
      integer idir
      integer idisp
      integer ierror
      integer isource
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS
      isource = 0
      idest = 0

      return
      end
      subroutine mpi_comm_dup ( comm, comm_out, ierror )

c*********************************************************************72
c
cc MPI_COMM_DUP duplicates a communicator.
c
      implicit none

      integer comm
      integer comm_out
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS
      comm_out = comm

      return
      end
      subroutine mpi_comm_free ( comm, ierror )

c*********************************************************************72
c
cc MPI_COMM_FREE "frees" a communicator.
c
      implicit none

      integer comm
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_comm_rank ( comm, me, ierror )

c*********************************************************************72
c
cc MPI_COMM_RANK reports the rank of the calling process.
c
      implicit none

      integer comm
      integer ierror
      integer me
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS
      me = 0

      return
      end
      subroutine mpi_comm_size ( comm, nprocs, ierror )

c*********************************************************************72
c
cc MPI_COMM_SIZE reports the number of processes in a communicator.
c
      implicit none

      integer comm
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )
      integer nprocs

      ierror = MPI_SUCCESS
      nprocs = 1

      return
      end
      subroutine mpi_comm_split ( comm, icolor, ikey, comm_new,
     &  ierror )

c*********************************************************************72
c
cc MPI_COMM_SPLIT splits up a communicator based on a key.
c
      implicit none

      integer comm
      integer comm_new
      integer icolor
      integer ierror
      integer ikey
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_copy_double_precision ( data1, data2, n, ierror )

c*********************************************************************72
c
cc MPI_COPY_DOUBLE copies a double precision vector.
c
      implicit none

      integer n

      double precision data1(n)
      double precision data2(n)
      integer i
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      do i = 1, n
        data2(i) = data1(i)
      end do

      return
      end
      subroutine mpi_copy_integer ( data1, data2, n, ierror )

c*********************************************************************72
c
cc MPI_COPY_INTEGER copies an integer vector.
c
      implicit none

      integer n

      integer data1(n)
      integer data2(n)
      integer i
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      do i = 1, n
        data2(i) = data1(i)
      end do

      return
      end
      subroutine mpi_copy_real ( data1, data2, n, ierror )

c*********************************************************************72
c
      implicit none

      integer n

      real data1(n)
      real data2(n)
      integer i
      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      do i = 1, n
        data2(i) = data1(i)
      end do

      return
      end
      subroutine mpi_finalize ( ierror )

c*********************************************************************72
c
cc MPI_FINALIZE shuts down the MPI library.
c
      implicit none

      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_get_count ( istatus, datatype, icount, ierror )

c*********************************************************************72
c
cc MPI_GET_COUNT reports the actual number of items transmitted.
c
      implicit none

      integer datatype
      integer icount
      integer ierror
      integer istatus
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_GET_COUNT - Error!'
      write ( *, '(a)' ) '  Should not query message from self.'

      return
      end
      subroutine mpi_init ( ierror )

c*********************************************************************72
c
cc MPI_INIT initializes the MPI library.
c
      implicit none

      integer ierror
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_SUCCESS

      return
      end
      subroutine mpi_irecv ( data, n, datatype, iproc, itag,
     &  comm, irequest, ierror )

c*********************************************************************72
c
cc MPI_IRECV receives data from another process.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer irequest
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_IRECV - Error!'
      write ( *, '(a)' ) '  Should not recv message from self.'

      return
      end
      subroutine mpi_isend ( data, n, datatype, iproc, itag,
     &  comm, request, ierror )

c*********************************************************************72
c
cc MPI_ISEND sends data from one process to another using nonblocking transmission.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )
      integer request

      request = 0
      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_ISEND - Error!'
      write ( *, '(a)' )  '  Should not send message to self.'

      return
      end
      subroutine mpi_recv ( data, n, datatype, iproc, itag,
     &  comm, istatus, ierror )

c*********************************************************************72
c
cc MPI_RECV receives data from another process within a communicator.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer istatus
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_RECV - Error!'
      write ( *, '(a)' ) '  Should not recv message from self.'

      return
      end
      subroutine mpi_reduce ( data1, data2, n, datatype, operation,
     &  receiver, comm, ierror )

c*********************************************************************72
c
cc MPI_REDUCE carries out a reduction operation.
c
      implicit none

      include "mpi_dummy.h"

      integer n

      integer comm
      integer data1(n)
      integer data2
      integer datatype
      integer ierror
      integer operation
      integer receiver

      ierror = MPI_SUCCESS

      if ( datatype .eq. mpi_double_precision ) then

        call mpi_reduce_double_precision ( 
     &    data1, data2, n, operation, ierror )

      else if ( datatype .eq. mpi_integer ) then

        call mpi_reduce_integer ( 
     &    data1, data2, n, operation, ierror )

      else if ( datatype .eq. mpi_real ) then

        call mpi_reduce_real ( 
     &    data1, data2, n, operation, ierror )

      else

        ierror = MPI_FAILURE

      end if

      return
      end
      subroutine mpi_reduce_double_precision ( 
     &  data1, data2, n, operation, ierror )

c*********************************************************************72
c
cc MPI_REDUCE_DOUBLE_PRECISION carries out a reduction operation on double precision values.
c
      implicit none

      include "mpi_dummy.h"

      integer n

      double precision data1(n)
      double precision data2(n)
      integer i
      integer ierror
      integer operation


      ierror = MPI_SUCCESS

      do i = 1, n
        data2(i) = data1(i)
      end do

      return
      end

      subroutine mpi_reduce_integer ( 
     &  data1, data2, n, operation, ierror )

c*********************************************************************72
c
      implicit none

      include "mpi_dummy.h"

      integer n

      integer data1(n)
      integer data2(n)
      integer i
      integer ierror
      integer operation

      ierror = MPI_SUCCESS

      do i = 1, n
         data2(i) = data1(i)
      end do

      ierror = MPI_FAILURE

      return
      end
      subroutine mpi_reduce_real ( 
     &  data1, data2, n, operation, ierror )

c*********************************************************************72
c
cc MPI_REDUCE_REAL carries out a reduction operation on reals.
c
c  Discussion:
c
      implicit none

      include "mpi_dummy.h"

      integer n

      real data1(n)
      real data2(n)
      integer i
      integer ierror
      integer operation

      ierror = MPI_SUCCESS

        do i = 1, n
          data2(i) = data1(i)
        end do

      return
      end
      subroutine mpi_reduce_scatter ( data1, data2, n, datatype,
     &  operation, comm, ierror )

c*********************************************************************72
c
cc MPI_REDUCE_SCATTER collects a message of the same length from each process.
c
      implicit none

      include "mpi_dummy.h"

      integer n

      integer comm
      integer data1(n)
      integer data2(n)
      integer datatype
      integer ierror
      integer operation

      ierror = MPI_SUCCESS

      if ( datatype .eq. mpi_double_precision ) then
        call mpi_copy_double_precision ( data1, data2, n, ierror )
      else if ( datatype .eq. mpi_integer ) then
        call mpi_copy_integer ( data1, data2, n, ierror )
      else if ( datatype .eq. mpi_real ) then
        call mpi_copy_real ( data1, data2, n, ierror )
      else
        ierror = MPI_FAILURE
      end if

      return
      end
      subroutine mpi_rsend ( data, n, datatype, iproc, itag,
     &  comm, ierror )

c*********************************************************************72
c
cc MPI_RSEND "ready sends" data from one process to another.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_RSEND - Error!'
      write ( *, '(a)' ) '  Should not send message to self.'

      return
      end
      subroutine mpi_send ( data, n, datatype, iproc, itag,
     &  comm, ierror )

c*********************************************************************72
c
cc MPI_SEND sends data from one process to another.
c
      implicit none

      integer n

      integer comm
      integer data(n)
      integer datatype
      integer ierror
      integer iproc
      integer itag
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_SEND - Error!'
      write ( *, '(a)' )  '  Should not send message to self.'

      return
      end
      subroutine mpi_wait ( irequest, istatus, ierror )

c*********************************************************************72
c
cc MPI_WAIT waits for an I/O request to complete.
c
      implicit none

      integer ierror
      integer irequest
      integer istatus
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_WAIT - Error!'
      write ( *, '(a)' ) '  Should not wait on message from self.'

      return
      end
      subroutine mpi_waitall ( icount, irequest, istatus, ierror )

c*********************************************************************72
c
cc MPI_WAITALL waits until all I/O requests have completed.
c
      implicit none

      integer icount
      integer ierror
      integer irequest
      integer istatus
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_WAITALL - Error!'
      write ( *, '(a)' ) '  Should not wait on message from self.'

      return
      end
      subroutine mpi_waitany ( icount, array_of_requests, index, 
     &  istatus, ierror )

c*********************************************************************72
c
cc MPI_WAITANY waits until one I/O requests has completed.
c
      implicit none

      integer array_of_requests(*)
      integer icount
      integer ierror
      integer index
      integer istatus
      integer MPI_FAILURE
      parameter ( MPI_FAILURE = 1 )
      integer MPI_SUCCESS
      parameter ( MPI_SUCCESS = 0 )

      ierror = MPI_FAILURE

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MPI_WAITANY - Error!'
      write ( *, '(a)' ) '  Should not wait on message from self.'

      return
      end
      function mpi_wtick ( )

c*********************************************************************72
c
cc MPI_WTICK returns the time between clock ticks.
c
      implicit none
      
      double precision mpi_wtick
      
      mpi_wtick = 1.0D+00
      
      return
      end
      function mpi_wtime ( )

c*********************************************************************72
c
cc MPI_WTIME returns the elapsed wall clock time.
c
      implicit none

      real*8 mpi_wtime
      real*4 a(2),etime
      a(1)=0.0
      a(2)=0.0
      mpi_wtime = etime(a)

      return
      end

      subroutine mpi_initialized(mpi_is_initialized, ierr)

      mpi_is_initialized = 0 
      ierr = 0

      return
      end

      subroutine mpi_comm_create(icomm,igroup,icommd,ierr)

      icommd = 1

      return
      end

      subroutine mpi_comm_group(icomm,igroup,ierr)

      igroup = 1
      ierr = 0

      return
      end

      subroutine mpi_group_free

      return
      end

      subroutine mpi_attr_get(icomm,ikey,ival,iflag,ierr)
 
      logical iflag

      ival =  9 999 999  ! dummy
 
      return
      end
uuuu
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

      ! set word size for REAL
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

      ! set word size for LOGICAL
      lsize = 4

      ! set word size for CHARACTER
      csize = 1
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

      call crystal_setup(cr_h,nekcomm,np)  ! set cr handle to new instance

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
          call exitt
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
         write(6,'(3(A,1g13.5,A,/))') 
     &      'total elapsed time         : ',ttotal, ' sec',
     &      'total solve time incl. I/O : ',ttime , ' sec',
     &      'time/timestep              : ',ttime/max(istep,1), ' sec'
      endif 
      call flush_io

      call mpi_finalize (ierr)
      call exit(0)

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
      function igl_running_sum(in)
c
c     Global vector commutative operation using spanning tree.
c

      include 'mpif.h'
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
      integer status(mpi_status_size)
      integer x,w,r


      x = in  ! running sum
      w = in  ! working buff
      r = 0   ! recv buff

#ifndef USERMPICOLL
      call mpi_scan(x,r,1,mpi_integer,mpi_sum,nekcomm,ierr)
      igl_running_sum = r
#else
      log2p = log2(np)
      mp    = 2**log2p
      lim   = log2P
      if (mp.ne.np) lim = log2P+1

      do l=1,lim
         mtype = l
         jid   = 2**(l-1)
         jid   = xor(nid,jid)   ! Butterfly, not recursive double

         if (jid.lt.np) then
            call mpi_irecv (r,1,mpi_integer,mpi_any_source,mtype
     $                                            ,nekcomm,msg,ierr)
            call mpi_send  (w,1,mpi_integer,jid,mtype,nekcomm,ierr)
            call mpi_wait  (msg,status,ierr)
            w = w+r
            if (nid.gt.jid) x = x+r
         endif
c        write(6,1) l,nid,jid,r,w,x,'summer'
c   1    format(2i6,'nid',4i6,1x,a6)
      enddo

      igl_running_sum = x
c     write(6,2) nid,in,x,'running sum'
c   2 format(3i9,1x,a6)
#endif

      return
      end

