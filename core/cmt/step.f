      subroutine setdtcmt
!--------------------------------------------------------------
! JH072914 now summed instead of maximized for compressible flow
! JH082015 why was I doing that again?
!     someday compute new timestep based on CFL-condition. follow
!     setdtc in nek/subs1.f very carefully.
!--------------------------------------------------------------
      include 'SIZE'
      include 'GEOM'
      include 'MVGEOM'
      include 'MASS'
      include 'INPUT'
      include 'SOLN'
      include 'TSTEP'
      COMMON /CMTGASPROP/  CSOUND(lx1,ly1,lz1,lelcmt)
      real                 csound

!--------------------------------------------------------------
! YOU REALLY PROBABLY WANT YOUR OWN DAMN SCRATCH SPACE FOR THIS
      common /scrsf/ u(lx1,ly1,lz1,lelcmt) ! mind if I borrow this?
     $ ,             v(lx1,ly1,lz1,lelcmt) ! as long as the mesh
     $ ,             w(lx1,ly1,lz1,lelcmt) ! doesn't move
! YOU REALLY PROBABLY WANT YOUR OWN DAMN SCRATCH SPACE FOR THAT
!--------------------------------------------------------------

      common /udxmax/ umax
      REAL VCOUR
      SAVE VCOUR

      NTOT   = NX1*NY1*NZ1*NELV
      NTOTL  = LX1*LY1*LZ1*lelcmt
      NTOTD  = NTOTL*NDIM
      COLD   = COURNO
      CMAX   = 1.2*CTARG
      CMIN   = 0.8*CTARG
      do i=1,ntot
         u(i,1,1,1)=abs(vx(i,1,1,1))+csound(i,1,1,1)
         v(i,1,1,1)=abs(vy(i,1,1,1))+csound(i,1,1,1)
         w(i,1,1,1)=abs(vz(i,1,1,1))+csound(i,1,1,1)
      enddo
      CALL CUMAX (u,v,w,UMAX)
      COURNO = DT*UMAX
      VOLD   = VCOUR
      VCOUR  = UMAX
      if (nio.eq.0) WRITE(6,100)ISTEP,TIME,DT,COURNO
 100  FORMAT('CMT ',I7,', t=',1pE14.7,', DT=',1pE14.7
     $,', C=',0pF7.3)
!     Zero DT
!
!     IF (DT .EQ. 0.0) THEN
!
!        IF (UMAX .NE. 0.0) THEN
!           DT = CTARG/UMAX
!           VCOUR = UMAX
      return
      end
