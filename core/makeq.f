c-----------------------------------------------------------------------
      subroutine makeq
C
C     Generate forcing function for the solution of a passive scalar.
C     !! NOTE: Do not change the content of the array BQ until the current
C              time step is completed 
C
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      logical  ifturb,if_conv_std
c
      if_conv_std = .true.
      if (ifmhd.and.ifaxis) if_conv_std = .false. ! conv. treated in induct.f
C
      CALL WHATFLD (IFTURB)
C
      IF (IFTURB)                                          CALL MAKETQ
      IF (.NOT.IFTURB                   .and.if_conv_std)  CALL MAKEUQ
      IF (IFADVC(IFIELD).AND..NOT.IFCHAR.and.if_conv_std)  CALL CONVAB
      IF (IFMVBD)                                          CALL ADMESHT
      IF (IFTRAN.AND..NOT.IFCVODE)                         CALL MAKEABQ
      IF ((IFTRAN.AND..NOT.IFCHAR.AND..NOT.IFCVODE).OR.
     $    (IFTRAN.AND..NOT.IFADVC(IFIELD).AND.IFCHAR.AND..NOT.IFCVODE))
     $                                                     CALL MAKEBDQ
      IF (IFADVC(IFIELD).AND.IFCHAR.AND..NOT.IFMVBD)       CALL CONVCH  
C
      return
      end
C

