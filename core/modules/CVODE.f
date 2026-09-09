      module cvode_mod

      use size_mod
      implicit none
c
c     Variables for the method of characteristics
c
      integer, parameter :: cv_lysize = lcvx1*lcvy1*lcvz1*lcvelt*ldimt

      integer igstype,itmeth
      parameter(
     &     itmeth  = 2,         ! newton iter
     &     igstype = 1          ! gs 1: modified 2: classical
     &         )

      integer,   allocatable, target :: cb_icvode(:)
      integer*8, allocatable, target :: cb_ilcvode(:)
      logical,   allocatable, target :: cb_lcvode(:)
      real,      allocatable, target :: cb_rcvode(:)
      real,      allocatable, target :: cb_cvrstat(:)
      integer*8, allocatable, target :: cb_cvistat(:)

      integer, pointer :: cv_nfld, cv_iatol
      integer, pointer :: cv_maxl, cv_itask, cv_ipretype

      integer*8, pointer :: cv_nlocal, cv_nglobal

      logical, pointer :: ifcvodeinit, ifdqj, ifcvfun

      real, pointer :: cv_atol(:)
      real, pointer :: cv_dtlag(:), cv_abmsh(:), cv_ab(:), cv_bd(:)
      real, pointer :: cv_rtol, cv_sigs, cv_delt
      real, pointer :: cv_time, cv_timel, cv_dtnek, cv_dtmax

      real, pointer :: nfe_avg, nli_nni_avg

      integer*8, pointer :: cv_istep
      integer*8, pointer :: iout_save(:)

      contains

      subroutine init
         implicit none

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_icvode(5), stat=ierr)
         if (ierr.ne.0) call exitti('alloc cb_icvode$',ierr)
         cb_icvode = 0
                          ! cv_nfld,cv_iatol,cv_maxl,cv_itask,cv_ipretype

         allocate(cb_ilcvode(2), stat=ierr)    ! cv_nlocal, cv_nglobal
         if (ierr.ne.0) call exitti('alloc cb_ilcvode$',ierr)
         cb_ilcvode = 0_8

         allocate(cb_lcvode(3), stat=ierr)   ! ifcvodeinit,ifdqj,ifcvfun
         if (ierr.ne.0) call exitti('alloc cb_lcvode$',ierr)
         cb_lcvode = .false.

         allocate(cb_rcvode(
     $        cv_lysize                ! cv_atol
     $      + 3                        ! cv_dtlag
     $      + 3                        ! cv_abmsh
     $      + 3                        ! cv_ab
     $      + 4                        ! cv_bd
     $      + 1                        ! cv_rtol
     $      + 1                        ! cv_sigs
     $      + 1                        ! cv_delt
     $      + 1                        ! cv_time
     $      + 1                        ! cv_timel
     $      + 1                        ! cv_dtnek
     $      + 1), stat=ierr)           ! cv_dtmax
         if (ierr.ne.0) call exitti('alloc cb_rcvode$',ierr)
         cb_rcvode = 0

         allocate(cb_cvrstat(2), stat=ierr)   ! nfe_avg, nli_nni_avg
         if (ierr.ne.0) call exitti('alloc cb_cvrstat$',ierr)
         cb_cvrstat = 0

         allocate(cb_cvistat(1 + 21), stat=ierr) ! cv_istep, iout_save
         if (ierr.ne.0) call exitti('alloc cb_cvistat$',ierr)
         cb_cvistat = 0_8

c        Group 1: /icvode/
         ioff = 1
         cv_nfld => cb_icvode(ioff)
         ioff = ioff + 1
         cv_iatol => cb_icvode(ioff)
         ioff = ioff + 1
         cv_maxl => cb_icvode(ioff)
         ioff = ioff + 1
         cv_itask => cb_icvode(ioff)
         ioff = ioff + 1
         cv_ipretype => cb_icvode(ioff)

c        Group 2: /ilcvode/
         ioff = 1
         cv_nlocal => cb_ilcvode(ioff)
         ioff = ioff + 1
         cv_nglobal => cb_ilcvode(ioff)

c        Group 3: /lcvode/
         ioff = 1
         ifcvodeinit => cb_lcvode(ioff)
         ioff = ioff + 1
         ifdqj => cb_lcvode(ioff)
         ioff = ioff + 1
         ifcvfun => cb_lcvode(ioff)

c        Group 4: /rcvode/
         ioff = 1
         cv_atol(1:cv_lysize) =>
     $         cb_rcvode(ioff : ioff + cv_lysize - 1)
         ioff = ioff + cv_lysize
         cv_dtlag(1:3) => cb_rcvode(ioff : ioff + 3 - 1)
         ioff = ioff + 3
         cv_abmsh(1:3) => cb_rcvode(ioff : ioff + 3 - 1)
         ioff = ioff + 3
         cv_ab(1:3) => cb_rcvode(ioff : ioff + 3 - 1)
         ioff = ioff + 3
         cv_bd(1:4) => cb_rcvode(ioff : ioff + 4 - 1)
         ioff = ioff + 4
         cv_rtol => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_sigs => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_delt => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_time => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_timel => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_dtnek => cb_rcvode(ioff)
         ioff = ioff + 1
         cv_dtmax => cb_rcvode(ioff)

c        Group 5: /cvrstat/
         ioff = 1
         nfe_avg => cb_cvrstat(ioff)
         ioff = ioff + 1
         nli_nni_avg => cb_cvrstat(ioff)

c        Group 6: /cvistat/
         ioff = 1
         cv_istep => cb_cvistat(ioff)
         ioff = ioff + 1
         iout_save(1:21) => cb_cvistat(ioff : ioff + 21 - 1)

      end subroutine init
      end module cvode_mod
