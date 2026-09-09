      module avg_mod

      implicit none
c
c     Averaged fields
c
      real, allocatable, target :: cb_avgcmnr(:)
      real, allocatable, target :: cb_chkavg(:)
      real, allocatable, target :: cb_chkrms(:)

      real, pointer :: atime, timel

      real, pointer :: uavg(:,:,:,:), vavg(:,:,:,:), wavg(:,:,:,:)
      real, pointer :: tavg(:,:,:,:,:)
      real, pointer :: pavg(:,:,:,:)

      real, pointer :: urms(:,:,:,:), vrms(:,:,:,:), wrms(:,:,:,:)
      real, pointer :: trms(:,:,:,:,:)
      real, pointer :: prms(:,:,:,:)
      real, pointer :: vwms(:,:,:,:), wums(:,:,:,:), uvms(:,:,:,:)

      contains

      subroutine init
         use size_mod
         implicit none

         integer ierr, ioff

c        --- allocate backing arrays ---

         allocate(cb_avgcmnr(2), stat=ierr)          ! atime, timel
         if (ierr.ne.0) call exitti('alloc cb_avgcmnr$',ierr)
         cb_avgcmnr = 0

         allocate(cb_chkavg(
     $        ax1*ay1*az1*lelt         ! uavg
     $      + ax1*ay1*az1*lelt         ! vavg
     $      + ax1*ay1*az1*lelt         ! wavg
     $      + ax1*ay1*az1*lelt*ldimt   ! tavg
     $      + ax2*ay2*az2*lelt), stat=ierr) ! pavg
         if (ierr.ne.0) call exitti('alloc cb_chkavg$',ierr)
         cb_chkavg = 0

         allocate(cb_chkrms(
     $        ax1*ay1*az1*lelt         ! urms
     $      + ax1*ay1*az1*lelt         ! vrms
     $      + ax1*ay1*az1*lelt         ! wrms
     $      + ax1*ay1*az1*lelt*ldimt   ! trms
     $      + ax2*ay2*az2*lelt         ! prms
     $      + ax1*ay1*az1*lelt         ! vwms
     $      + ax1*ay1*az1*lelt         ! wums
     $      + ax1*ay1*az1*lelt), stat=ierr) ! uvms
         if (ierr.ne.0) call exitti('alloc cb_chkrms$',ierr)
         cb_chkrms = 0

c        Group 1: /avgcmnr/
         ioff = 1
         atime => cb_avgcmnr(ioff)
         ioff = ioff + 1
         timel => cb_avgcmnr(ioff)

c        Group 2: /chkavg/
         ioff = 1
         uavg(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkavg(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         vavg(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkavg(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         wavg(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkavg(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         tavg(1:ax1,1:ay1,1:az1,1:lelt,1:ldimt) =>
     $         cb_chkavg(ioff : ioff + ax1*ay1*az1*lelt*ldimt - 1)
         ioff = ioff + ax1*ay1*az1*lelt*ldimt
         pavg(1:ax2,1:ay2,1:az2,1:lelt) =>
     $         cb_chkavg(ioff : ioff + ax2*ay2*az2*lelt - 1)

c        Group 3: /chkrms/
         ioff = 1
         urms(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         vrms(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         wrms(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         trms(1:ax1,1:ay1,1:az1,1:lelt,1:ldimt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt*ldimt - 1)
         ioff = ioff + ax1*ay1*az1*lelt*ldimt
         prms(1:ax2,1:ay2,1:az2,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax2*ay2*az2*lelt - 1)
         ioff = ioff + ax2*ay2*az2*lelt
         vwms(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         wums(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)
         ioff = ioff + ax1*ay1*az1*lelt
         uvms(1:ax1,1:ay1,1:az1,1:lelt) =>
     $         cb_chkrms(ioff : ioff + ax1*ay1*az1*lelt - 1)

      end subroutine init
      end module avg_mod
