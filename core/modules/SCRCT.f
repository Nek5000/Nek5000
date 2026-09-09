      module scrct_mod

      implicit none
C
C     Scratch arrays used in CONNECT and associated subroutines.
C
C     ring pass arrays  (was COMMON /SCRMG/)
C
      real, allocatable, target :: cb_scrmg_ct(:)
      real, pointer :: rmxs(:),rmax(:)
     $               , xcg(:),ycg(:),zcg(:)
     $               , xgs(:),ygs(:),zgs(:)
     $               , xml(:,:,:,:),xms(:,:,:,:)
     $               , yml(:,:,:,:),yms(:,:,:,:)
     $               , zml(:,:,:,:),zms(:,:,:,:)
C
C     was COMMON /SCREV/
C
      real, allocatable, target :: cb_screv_ct(:)
      real, pointer :: side(:,:,:), sides(:,:,:)
C
C     was COMMON /CTMP1/
C
      real, allocatable, target :: cb_ctmp1_ct(:)
      real, pointer :: flag(:,:,:,:),tmp2(:,:,:,:)
     $               , lmult(:,:,:,:),bcs(:,:,:),xyz(:,:,:)
C
C     nested dissection arrays  (was COMMON /SCRVH/, /SCRCH/)
C
      character*3, allocatable, target :: cbcs(:,:)

      integer, allocatable, target :: ibrnch(:), nbrnch(:), list(:)
     $                              , list1(:), list2(:)
      logical, allocatable, target :: ifcnst(:,:)
C
C     was: real xyzl(3,8,lelt),cg(3,lelt); equivalence(xyzl,xms),
C     (cg,xgs) -- xyzl/cg alias the start of xms/xgs's storage
C
      real, pointer :: xyzl(:,:,:), cg(:,:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff
         integer n_scrmg

         n_scrmg = 8*lelt + 54*lzl*lelt

         allocate(cb_scrmg_ct(n_scrmg), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_scrmg_ct$',ierr)
         cb_scrmg_ct = 0
         allocate(cb_screv_ct(2*4*6*lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_screv_ct$',ierr)
         cb_screv_ct = 0
         allocate(cb_ctmp1_ct(3*3*3*lzl*lelt + 5*6*lelt + 3*8*lelt),
     $            stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_ctmp1_ct$',ierr)
         cb_ctmp1_ct = 0
         allocate(cbcs(6,lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cbcs$',ierr)
         cbcs = ''
         allocate(ibrnch(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ibrnch$',ierr)
         ibrnch = 0
         allocate(nbrnch(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc nbrnch$',ierr)
         nbrnch = 0
         allocate(list(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc list$',ierr)
         list = 0
         allocate(list1(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc list1$',ierr)
         list1 = 0
         allocate(list2(lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc list2$',ierr)
         list2 = 0
         allocate(ifcnst(6,lelt), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ifcnst$',ierr)
         ifcnst = .false.

         ioff = 1
         rmxs(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         rmax(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         xcg(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         ycg(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         zcg(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         xgs(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         ygs(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         zgs(1:lelt) => cb_scrmg_ct(ioff : ioff+lelt-1)
         ioff = ioff + lelt
         xml(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         xms(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         yml(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         yms(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         zml(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         zms(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_scrmg_ct(ioff : ioff+9*lzl*lelt-1)

      ioff = 1
      side(1:4,1:6,1:lelt) => cb_screv_ct(ioff : ioff + 4*6*lelt - 1)
      ioff = ioff + 4*6*lelt
      sides(1:4,1:6,1:lelt) => cb_screv_ct(ioff : ioff + 4*6*lelt - 1)

         ioff = 1
         flag(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_ctmp1_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         tmp2(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_ctmp1_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         lmult(1:3,1:3,1:lzl,1:lelt) =>
     $      cb_ctmp1_ct(ioff : ioff+9*lzl*lelt-1)
         ioff = ioff + 9*lzl*lelt
         bcs(1:5,1:6,1:lelt) => cb_ctmp1_ct(ioff : ioff+5*6*lelt-1)
         ioff = ioff + 5*6*lelt
         xyz(1:3,1:8,1:lelt) => cb_ctmp1_ct(ioff : ioff+3*8*lelt-1)

         xyzl(1:3,1:8,1:lelt) => cb_scrmg_ct(
     $      8*lelt+9*lzl*lelt+1 : 8*lelt+9*lzl*lelt+3*8*lelt)
         cg(1:3,1:lelt) => cb_scrmg_ct(5*lelt+1 : 5*lelt+3*lelt)

      end subroutine init
      end module scrct_mod
