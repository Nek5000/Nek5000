c=======================================================================
      program int_tp
c
c     Interpolate from a tensor-product SE geometry to a regular tp array
c     NOTE:  COMPILE in real*4 mode only, because of byte_read/write.
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'

      integer mode
      character*80 fname

      call proc_init

      write(*,'A') 'Enter the fld filename:'
      read(*,*) fname

      call get_fld_data (fname,ierr)
      if (ierr.ne.0) call exitt

      ! interpolate data for a bunch of given points
      call interpolate

      end

c=======================================================================

      subroutine proc_init

      include 'SIZE'
      include 'INPUT'
      include 'GEOM'

      one = 1.0
      eps = 1.e-12

      wdsize = 8
      if (one+eps.eq.1.0) wdsize = 4

      return
      end

c=======================================================================
      subroutine interpolate
c=======================================================================
c
c Example: Resampling VX,Vy,VZ of an existing fld-file 
c          on a NELX^3 point cube using spectral
c          interpolation. The interpolation is top-level parallelized using
c          simple OpenMP pragmas.
c          The interpolated data is saved in a binary file called
c          int.dat 
c-----------------------------------------------------------------------

      include 'SIZE'

      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MISC'

      real value(10)
      real point(3),counter

      data point / 0.,0.,0. /

      integer npoints
      parameter(npoints=128)
 
      real*4 tdump(6,npoints,npoints)

      real dx,dy,dz
      integer nx,ny,nz
      integer average
       
      ! define cube point mesh
      nelx=npoints
      nely=nelx
      nelz=nelx
      dx=1.0/(nelx-1)
      dy=dx
      dz=dx

      ! open binary output file
      call byte_open("int.dat" // CHAR(0))

      do iz=1,nelz
         point(3) = float(iz-1)*dz
         do iy=1,nely
            point(2) = float(iy-1)*dy
            do ix=1,nelx
               point(1) = float(ix-1)*dx
               ! do interpolation
               call interpl_gen(value(1),point,VX)
               call interpl_gen(value(2),point,VY)
               call interpl_gen(value(3),point,VZ)
               ! save results in dumping array
               call copy4(tdump(1,ix,iy),point,3)
               call copy4(tdump(4,ix,iy),value,3)
            enddo
         enddo
         call byte_write(tdump,6*nelx*nely)
         write(*,*) 100.0*float(iz)/float(nelz), '% done'
      enddo

      call byte_close
  
      return
      end

c=======================================================================

      subroutine copy4(a,b,n)
      real*4   a(1)
      real     b(1)
      do i = 1, n
         a(i) = b(i)
      enddo
      return
      end

