      module neknek_mod

      implicit none
c
c     Multimesh variables
c
      integer, allocatable, target :: cb_intflag(:)
      integer, allocatable, target :: cb_intmask(:)
      real,    allocatable, target :: cb_valmask(:)
      integer, allocatable, target :: cb_cgeom(:)
      integer, allocatable, target :: cb_inbc(:)
      real,    allocatable, target :: cb_mybd(:)
      real,    allocatable, target :: cb_multipts_r(:)
      integer, allocatable, target :: cb_multipts_i(:)
      integer, allocatable, target :: cb_intp_h_nn(:)

      integer, pointer :: intflag(:,:)

      integer, pointer :: imask(:,:,:,:)

      real, pointer :: valint(:,:,:,:,:)

      integer, pointer :: igeom

      integer, pointer :: nfld_neknek

      real, pointer :: bdrylg(:,:,:)

      real, pointer :: rst(:)

      integer, pointer :: rcode(:), elid(:)
      integer, pointer :: proc(:), ilist(:,:), npoints_nn

      integer, pointer :: inth_multi2

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_intflag(2*ldim*lelt), stat=ierr)     ! intflag
         if (ierr.ne.0) call exitti('alloc cb_intflag$',ierr)
         cb_intflag = 0

         allocate(cb_intmask(lx1*ly1*lz1*lelt), stat=ierr) ! imask
         if (ierr.ne.0) call exitti('alloc cb_intmask$',ierr)
         cb_intmask = 0

         allocate(cb_valmask(lx1*ly1*lz1*lelt*nfldmax_nn),
     $            stat=ierr)                              ! valint
         if (ierr.ne.0) call exitti('alloc cb_valmask$',ierr)
         cb_valmask = 0

         allocate(cb_cgeom(1), stat=ierr)                 ! igeom
         if (ierr.ne.0) call exitti('alloc cb_cgeom$',ierr)
         cb_cgeom = 0

         allocate(cb_inbc(1), stat=ierr)                  ! nfld_neknek
         if (ierr.ne.0) call exitti('alloc cb_inbc$',ierr)
         cb_inbc = 0

         allocate(cb_mybd(lx1*ly1*lz1*lelt*nfldmax_nn*3),
     $            stat=ierr)                              ! bdrylg
         if (ierr.ne.0) call exitti('alloc cb_mybd$',ierr)
         cb_mybd = 0

         allocate(cb_multipts_r(nmaxl_nn*ldim), stat=ierr) ! rst
         if (ierr.ne.0) call exitti('alloc cb_multipts_r$',ierr)
         cb_multipts_r = 0

         allocate(cb_multipts_i(
     $        nmaxl_nn                 ! rcode
     $      + nmaxl_nn                 ! elid
     $      + nmaxl_nn                 ! proc
     $      + 1*nmaxl_nn               ! ilist
     $      + 1), stat=ierr)           ! npoints_nn
         if (ierr.ne.0) call exitti('alloc cb_multipts_i$',ierr)
         cb_multipts_i = 0

         allocate(cb_intp_h_nn(1), stat=ierr)             ! inth_multi2
         if (ierr.ne.0) call exitti('alloc cb_intp_h_nn$',ierr)
         cb_intp_h_nn = 0

c        Group 1: /intflag/
         intflag(1:2*ldim,1:lelt) => cb_intflag(1 : 2*ldim*lelt)

c        Group 2: /intmask/
         imask(1:lx1,1:ly1,1:lz1,1:lelt) =>
     $         cb_intmask(1 : lx1*ly1*lz1*lelt)

c        Group 3: /valmask/
         valint(1:lx1,1:ly1,1:lz1,1:lelt,1:nfldmax_nn) =>
     $         cb_valmask(1 : lx1*ly1*lz1*lelt*nfldmax_nn)

c        Group 4: /cgeom/
         igeom => cb_cgeom(1)

c        Group 5: /inbc/
         nfld_neknek => cb_inbc(1)

c        Group 6: /mybd/
         bdrylg(1:lx1*ly1*lz1*lelt,1:nfldmax_nn,0:2) =>
     $         cb_mybd(1 : lx1*ly1*lz1*lelt*nfldmax_nn*3)

c        Group 7: /multipts_r/
         rst(1:nmaxl_nn*ldim) => cb_multipts_r(1 : nmaxl_nn*ldim)

c        Group 8: /multipts_i/
         ioff = 1
         rcode(1:nmaxl_nn) => cb_multipts_i(ioff : ioff+nmaxl_nn-1)
         ioff = ioff + nmaxl_nn
         elid(1:nmaxl_nn) => cb_multipts_i(ioff : ioff+nmaxl_nn-1)
         ioff = ioff + nmaxl_nn
         proc(1:nmaxl_nn) => cb_multipts_i(ioff : ioff+nmaxl_nn-1)
         ioff = ioff + nmaxl_nn
         ilist(1:1,1:nmaxl_nn) =>
     $         cb_multipts_i(ioff : ioff + 1*nmaxl_nn - 1)
         ioff = ioff + 1*nmaxl_nn
         npoints_nn => cb_multipts_i(ioff)

c        Group 9: /intp_h_nn/
         inth_multi2 => cb_intp_h_nn(1)

      end subroutine init
      end module neknek_mod
