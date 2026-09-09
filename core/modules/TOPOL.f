      module topol_mod
      include 'SIZE'
c
c     Arrays for direct stiffness summation
c
      integer nomlis(2,3),nmlinv(6),group(6),skpdat(6,6),eface(6)
     $       ,eface1(6)
      common /cfaces/ nomlis,nmlinv,group,skpdat,eface,eface1

      integer eskip(-12:12,3),nedg(3),ncmp
     $       ,ixcn(8),noffst(3,0:ldimt1)
     $       ,maxmlt,nspmax(0:ldimt1)
     $       ,ngspcn(0:ldimt1),ngsped(3,0:ldimt1)
      common /cedges/ eskip,nedg,ncmp,ixcn,noffst,maxmlt,nspmax
     $               ,ngspcn,ngsped

      integer, allocatable, target ::
     $     numscn(:,:),numsed(:,:)
     $   , gcnnum(:,:,:),lcnnum(:,:,:)
     $   , gednum(:,:,:),lednum(:,:,:)
     $   , gedtyp(:,:,:)
     $   , ngcomm(:,:)

      integer iedge(20),iedgef(2,4,6,0:1)
     $       ,indx(8),invedg(27)
      common /edges/ iedge,iedgef,indx,invedg

      integer iedgfc(4,6)
      DATA    IEDGFC /  5,7,9,11,  6,8,10,12,
     $                  1,3,9,10,  2,4,11,12,
     $                  1,2,5,6,   3,4,7,8    /

      integer icedg(3,16)
      DATA    ICEDG / 1,2,1,   3,4,1,   5,6,1,   7,8,1,
     $                1,3,2,   2,4,2,   5,7,2,   6,8,2,
     $                1,5,3,   2,6,3,   3,7,3,   4,8,3,
C      -2D-
     $                1,2,1,   3,4,1,   1,3,2,   2,4,2 /

      integer icface(4,10)
      DATA    ICFACE/ 1,3,5,7, 2,4,6,8,
     $                1,2,5,6, 3,4,7,8,
     $                1,2,3,4, 5,6,7,8,
C      -2D-
     $                1,3,0,0, 2,4,0,0,
     $                1,2,0,0, 3,4,0,0  /
C
      contains

      subroutine init
         implicit none
         integer ierr

         allocate(numscn(lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc numscn$',ierr)
         numscn = 0
         allocate(numsed(lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc numsed$',ierr)
         numsed = 0
         allocate(gcnnum(8,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc gcnnum$',ierr)
         gcnnum = 0
         allocate(lcnnum(8,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc lcnnum$',ierr)
         lcnnum = 0
         allocate(gednum(12,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc gednum$',ierr)
         gednum = 0
         allocate(lednum(12,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc lednum$',ierr)
         lednum = 0
         allocate(gedtyp(12,lelt,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc gedtyp$',ierr)
         gedtyp = 0
         allocate(ngcomm(2,0:ldimt1), stat=ierr)
         if (ierr.ne.0) call exitti('alloc ngcomm$',ierr)
         ngcomm = 0

      end subroutine init
      end module topol_mod
