      subroutine byte_sync_mpi(mpi_fh)

#ifdef MPIIO
      include 'mpif.h'
      call MPI_file_sync(mpi_fh,ierr)
#endif

      return
      end
C--------------------------------------------------------------------------
      subroutine byte_open_mpi(fname,mpi_fh,ifro,ierr)

      include 'SIZE'
      include 'RESTART'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

#ifdef MPIIO
      include 'mpif.h'

      character*132 fname
      logical ifro 

      imode = MPI_MODE_WRONLY+MPI_MODE_CREATE
      if(ifro) then
        imode = MPI_MODE_RDONLY 
      endif

      if(nid.eq.pid0 .or. nid.eq.pid0r) then
c        write(*,*) nid, 'call MPI_file_open',fname
        call MPI_file_open(nekcomm,fname,imode,
     &                     MPI_INFO_NULL,mpi_fh,ierr)
      endif
#else
      write(6,*) 'byte_open_mpi: No MPI-IO support!'
      ierr=1
      return
#endif
      ierr=0
      return
      end
C--------------------------------------------------------------------------
      subroutine byte_read_mpi(buf,icount,iorank,mpi_fh,ierr)

      include 'SIZE'
      include 'RESTART'

#ifdef MPIIO
      include 'mpif.h'

      real*4 buf(1)          ! buffer

      if(nid.eq.pid0 .or. nid.eq.pid0r) then
        iout = icount ! icount is in 4-byte words
        if(iorank.ge.0 .and. nid.ne.iorank) iout = 0
c        write(*,*) 'byte_read_mpi', nid, iout/4
        call MPI_file_read_all(mpi_fh,buf,iout,MPI_REAL,
     &                           MPI_STATUS_IGNORE,ierr)
      endif
#else
      write(6,*) 'byte_read_mpi: No MPI-IO support!'
      ierr=1
      return
#endif
     
      ierr=0

      return
      end
C--------------------------------------------------------------------------
      subroutine byte_write_mpi(buf,icount,iorank,mpi_fh,ierr)

      include 'SIZE'
      include 'RESTART'

#ifdef MPIIO
      include 'mpif.h'

      real*4 buf(1)          ! buffer

      if(nid.eq.pid0 .or. nid.eq.pid0r) then
        iout = icount ! icount is in 4-byte words
        if(iorank.ge.0 .and. nid.ne.iorank) iout = 0
c        write(*,*) 'byte_write', nid, iout/4
        call MPI_file_write_all(mpi_fh,buf,iout,MPI_REAL,
     &                          MPI_STATUS_IGNORE,ierr)
      endif
#else
      write(6,*) 'byte_write_mpi: No MPI-IO support!'
      ierr=1
      return
#endif
      ierr=0
      return
      end
C--------------------------------------------------------------------------
      subroutine byte_close_mpi(mpi_fh,ierr)

      include 'SIZE'
      include 'RESTART'

#ifdef MPIIO
      include 'mpif.h'
      if(nid.eq.pid0 .or. nid.eq.pid0r) then
        call MPI_file_close(mpi_fh,ierr)
      endif
#else
      if(nio.eq.0) write(6,*) 'byte_close_mpi: No MPI-IO support!'
      ierr=1
      return
#endif

      return
      end
C--------------------------------------------------------------------------
      subroutine byte_set_view(ioff_in,mpi_fh)

      include 'SIZE'
      include 'RESTART'

#ifdef MPIIO
      include 'mpif.h'
      integer*8 ioff_in
    
      if(nid.eq.pid0 .or. nid.eq.pid0r) then
         if(ioff_in.lt.0) then
           write(6,*) 'byte_set_view: offset<0!'
           call exitt
         endif
c         write(*,*) 'dataoffset', nid, ioff_in
         call MPI_file_set_view(mpi_fh,ioff_in,MPI_BYTE,MPI_BYTE,
     &                          'native',MPI_INFO_NULL,ierr)
      endif
#endif

      return
      end
