h1mg_solve
**********

h1mg_schwarz
~~~~~~~~~~~~

h1mg_schwarz_part1
==================

hsmg_schwarz_toext3d
--------------------

In this subroutine, ``rzero`` was inlined to avoid re-architecting ``rzero`` for OpenACC.  If we
rearchitected ``rzero`` to run on GPU, we would need to account for the case where: ``rzero`` was
called on the host, expected to zero-out ``A`` on the host but the argument ``A`` was present also
present on the device.  

.. todo:: 
  Opportunity exists to zero out fewer array elements in hsmg_schwarz_toext3d

Rearchitecting ``hsmg_schwarz_toext3d`` was done with the intent of running h1mg_solve on GPU.  The
former subroutine is also called from ``hsmg_solve``.  This *could* cause an issue if
``hsmg_solve`` is expected to run on CPU during the same execution, since we have a hard PRESENT on
the data

.. code-block:: fortranfixed

      subroutine hsmg_schwarz_toext3d(a,b,n)
      include 'SIZE'
      integer n
      real a(0:n+1,0:n+1,0:n+1,nelv),b(n,n,n,nelv)
      
      integer i,j,k,ie

      call rzero(a,(n+2)*(n+2)*(n+2)*nelv)

      do ie=1,nelv
      do k=1,n
      do j=1,n
      do i=1,n
         a(i,j,k,ie)=b(i,j,k,ie)
      enddo
      enddo
      enddo
      enddo
      return
      end

hsmg_extrude
------------

``hsmg_extrude`` was rewritten for the call from h1mg_schwarz_part_1.  It might be problematic if
called from elsewhere.

Outer loop over elements was parallelized over gangs.  Inner loops over i,j will not fill a vector
width if we collapse them alone (``COLLAPSE(2)``).  The opportunity exists to split the loop over
elements to allow ``COLLAPSE(3)``, but that would change serial exectution.  Is this okay?

.. todo::
  Possible loop fission in hsmg_extrude

Loop over elements was implemented with PRESENT_OR_COPY(arr1) to prevent it from crashing during
setup (h1mg_setup_schwarz_wt_1). This could be slow (could be unnecessary to run on GPU).  Options:
ifpresent check, different subroutines for host/devcie

Also called during hsmg_schwarz; same situation. 


.. code-block:: fortranfixed

      subroutine hsmg_extrude(arr1,l1,f1,arr2,l2,f2,nx,ny,nz)
      include 'SIZE'
      include 'INPUT'
      integer l1,l2,nx,ny,nz
      real arr1(nx,ny,nz,nelv),arr2(nx,ny,nz,nelv)
      real f1,f2
      
      integer i,j,k,ie,i0,i1
      i0=2
      i1=nx-1
      
      if(.not.if3d) then
         do ie=1,nelv
            do j=i0,i1
               arr1(l1+1 ,j,1,ie) = f1*arr1(l1+1 ,j,1,ie)
     $                             +f2*arr2(l2+1 ,j,1,ie)
               arr1(nx-l1,j,1,ie) = f1*arr1(nx-l1,j,1,ie)
     $                             +f2*arr2(nx-l2,j,1,ie)
            enddo
            do i=i0,i1
               arr1(i,l1+1 ,1,ie) = f1*arr1(i,l1+1 ,1,ie)
     $                             +f2*arr2(i,l2+1 ,1,ie)
               arr1(i,ny-l1,1,ie) = f1*arr1(i,ny-l1,1,ie)
     $                             +f2*arr2(i,nx-l2,1,ie)
            enddo
         enddo
      else
         do ie=1,nelv
            do k=i0,i1
            do j=i0,i1
               arr1(l1+1 ,j,k,ie) = f1*arr1(l1+1 ,j,k,ie)
     $                             +f2*arr2(l2+1 ,j,k,ie)
               arr1(nx-l1,j,k,ie) = f1*arr1(nx-l1,j,k,ie)
     $                             +f2*arr2(nx-l2,j,k,ie)
            enddo
            enddo
            do k=i0,i1
            do i=i0,i1
               arr1(i,l1+1 ,k,ie) = f1*arr1(i,l1+1 ,k,ie)
     $                             +f2*arr2(i,l2+1 ,k,ie)
               arr1(i,nx-l1,k,ie) = f1*arr1(i,nx-l1,k,ie)
     $                             +f2*arr2(i,nx-l2,k,ie)
            enddo
            enddo
            do j=i0,i1
            do i=i0,i1
               arr1(i,j,l1+1 ,ie) = f1*arr1(i,j,l1+1 ,ie)
     $                             +f2*arr2(i,j,l2+1 ,ie)
               arr1(i,j,nx-l1,ie) = f1*arr1(i,j,nx-l1,ie)
     $                             +f2*arr2(i,j,nx-l2,ie)
            enddo
            enddo
         enddo
      endif
      return
      end

hsschwarz_dssum
---------------

gs_op is called from gslib, which implements ifpresent checks, etc.  should work.  

hsmg_fdm
--------

Using #ifdef 

Downside: ``hsmg_solve`` call tree would crash in call to hsmg_do_fast_acc (similarly to
hsmg_schwarz_toext3d) because of hard PRESENT.  Would actually crash in hsmg_schwarz_toext3d first.  

hsmg_do_fast
............

A separate GPU subroutine (``hsmg_do_fast_acc``) was implmented so mxm could be inlined (compiler
issues with ACC ROUTINE).  Also avoids reimplementing the many CPU optimized mxm routines.
``hsmg_do_fast_acc`` is naive for GPU and slower on CPU.  

Downside: If hsmg_do_fast_acc is called on CPU during 

.. todo::
  Next steps for optimizing mxm on GPU?  Think about caching, as analogy on CPU

hsmg_schwarz_toreg3d
--------------------

This is called during setup (h1mg_setup_schwarz_wt_1) so we need PRESENT_OR_COPY.  Could be slow on
setup, not needed to run on GPU during setup, see above comment.

hsmg_schwarz_wt
===============

hsmg_schwarz_wt3d
-----------------

Could split loop over elements to allow COLLAPSE(3), but that would change serial execution.  

Doing hard PRESENT()

.. todo::
  Assess whether we can change loops

cmult
=====

Inlined cmult and used ACC LOOP.  However, since ``n`` is not known during the data declaration,
the PRESENT(e) clause might crash.  

copy
~~~~

Inlining this to avoid CPU/GPU if/else's. Hard PRESENT(r,rhs) might crash because size of r and rhs
aren't known.  Might need to explicitly declare size of r and rhs.

.. todo::

  At this point, the restriction operations could be more efficient on CPU, based on Matt's
  experience on NekCEM.  However, to execute the following routines on CPU, we would need to test
  for presence of data on GPU (right now, we have a hard PRESENT, assume we only run on GPU, and
  allow it to crash).  To overcome this problem:
  * Use acc_if_present within one version of subroutine, i.e:
      if data is on GPU
        run on GPU w/o copy
      else
        run on CPU w/o copy
      endif
    (assumes that data we want is on CPU)
  * Use two versions of subroutine
  Both of these involve significant changes to the codebase.  

h1mg_rstr
~~~~~~~~

h1mg_azm
~~~~~~~~

.. todo::
  We'll need to refactor this if_hybrid.  
  Note: Paul says this is always false!  No need to redo this...

hsmg_do_wt
==========

Another situtation where we could COLLAPSE(3) if we split loop over elements.  Another situation
that could crash b/c hard PRESENT.

hsmg_tsnr1
==========

hsmg_tnsr1_3d
-------------

Also implemented an alternate ACC subroutine, hsmg_tnsr1_3d_acc, with a naive inline/unrolling of
mxm.  hsmg_tnsr1_3d_acc will be less efficient on CPU.  Into reformulating to allow COLLAPSE(3) or
COLLAPSE(4).  

hsmg_coarse_solve
~~~~~~~~~~~~~~~~~

.. todo::
  A bunch of optimizations to study here!

Right now, only this coarse solve is done on the CPU.  On NekCEM, Matt saw better performance when
he was also doing the some restriction operations on the CPU.  We will attempt to redo some of his
performance studies.

Also, consider how much of the array ``e`` to UPDATE.

hsmg_intp
~~~~~~~~

hsmg_tnsr
`````````

hsmg_tnsr3d
:::::::::::

Also implemented an alternate ACC subroutine, hsmg_tnsr3d_acc, with a naive inline/unrolling of
mxm.  hsmg_tnsr1_3d_acc will be less efficient on CPU.  Into reformulating to allow COLLAPSE(3) or
COLLAPSE(4).  Does a hard PRESENT check

dsavg
~~~~~

Implemented dsavg_acc, which inlined col2 (wanted to parallelize loop but didn't want to change
math.f) and called dssum w/o modificatiod.  Hard PRESENT on vmult and tmult, which don't change
after initialization.





