      subroutine setprop
C------------------------------------------------------------------------
C
C     Set variable property arrays
C
C------------------------------------------------------------------------
      include 'SIZE'
      include 'TOTAL'
C
C     Caution: 2nd and 3rd strainrate invariants residing in scratch
C              common /SCREV/ are used in STNRINV and NEKASGN
C
      common /screv/ sii (lx1,ly1,lz1,lelt),siii(lx1,ly1,lz1,lelt)

      NXYZ1 = NX1*NY1*NZ1
      MFIELD=2
      IF (IFFLOW) MFIELD=1
      nfldt = nfield
      if (ifmhd) nfldt = nfield+1
      DO IFIELD=MFIELD,nfldt
         IF (IFSTRS .AND. IFIELD.EQ.1) CALL STNRINV ! expensive !

         CALL VPROPS

         nel = nelfld(ifield)
         vol = volfld(ifield)
         ntot1 = nxyz1*nel

         avdiff(ifield) = glsc2 (bm1,vdiff (1,1,1,1,ifield),ntot1)/vol
         avtran(ifield) = glsc2 (bm1,vtrans(1,1,1,1,ifield),ntot1)/vol

      ENDDO

C
      RETURN
      END

