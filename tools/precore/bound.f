C------------------------------------------------------------------------------
C
C                          NEKTON 2.6  2/8/90
C
C			Copyright (C) 1990, by the 
C
C		Massachusetts Institute of Technology  and Nektonics, Inc.
C
C All Rights Reserved
C
C This program is a licenced product of MIT and Nektonics, Inc.,  and it is 
C not to be disclosed to others, copied, distributed, or displayed 
C without prior authorization.
C
C------------------------------------------------------------------------------
C
c-----------------------------------------------------------------------
      subroutine overlap(iel,iside,ielo,isideo)
C
      include 'basics.inc'

C     Find closest element, side we want to duplicate, If there is no overlap,
C     returns original element and side
C
      RMIN=10000.0 * XFAC
      RMAX=    0.0
C
      IELO   = IEL
      ISIDEO = ISIDE
      DO 3 IIEL=1,NEL
         DO 1 IISIDE=1,NSIDES
            R=    (
     $       (SIDES(IEL,ISIDE,1)-SIDES(IIEL,IISIDE,1))**2
     $      +(SIDES(IEL,ISIDE,2)-SIDES(IIEL,IISIDE,2))**2
     $      +(SIDES(IEL,ISIDE,3)-SIDES(IIEL,IISIDE,3))**2
     $      )
            IF(IEL.NE.IIEL.OR.IISIDE.NE.ISIDE)THEN
               IF(R.GT.RMAX) RMAX = R
               IF(R.LT.RMIN) THEN
                  RMIN   = R
                  ISIDEO = IISIDE
                  IELO   = IIEL
               ENDIF
            ENDIF
1        CONTINUE
3     CONTINUE
C
      IF(RMIN .EQ. 10000.0 * XFAC)THEN
         CALL PRS('ERROR FINDING OVERLAPPING SIDE FOR B.C.$')
      ELSE IF (RMIN .GE. RMAX/100.)THEN
         CALL PRS('ERROR FINDING OVERLAPPING SIDE FOR B.C.$')
         CALL PRS('NO ADJACENT ELEMENT$')
         IELO   = IEL
         ISIDEO = ISIDE
      ENDIF
C
      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine INFLOW(IEL,ISIDE,IF,cbcI)
C! REAL KLUDGE USED NOT TO CORRUPT LINE  !!??DO THE POINTERS GET READ IN RIGHT?
      include 'basics.inc'
      CHARACTER*3 cbcI
C
C     Don't read inline stuff anymore... pff 8/30/93.
C
      if (nel.gt.0) return
C
      cbc3 = cbcI
      IF(cbcI .NE. 'on ')THEN
         CALL PRS(' Enter boundary conditions for this element $')
         CALL PRS
     $   (' as a function of X,Y [Z]& TIME in the style of a fortran$')
         CALL PRS(' statement.  Start in column 1.  Use one line for $')
      ENDIF
      IF(cbcI.EQ.'s  ')THEN
C        stress (Traction)
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRSIS(' each of the $',NLINES,
     $   ' TRACTION components$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TRY=SIN(X)$')
         IF(IF3d)CALL PRS(' TRZ=1.0$')
      ELSE IF(cbcI.EQ.'sl ')THEN
C        stress (Traction)
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRSIS(' each of the $',NLINES,
     $   ' TRACTION components$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRN=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
         IF(IF3D)CALL PRS(' TR2=SIN(X)$')
      ELSE IF(cbcI.EQ.'sh ')THEN
C        shear
         IF(.NOT.IF3D)THEN
            NLINES=2
            CALL PRS(' each of the 2 shear (TRaction) components $')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TRY=    Y  +  3 * TIME)$')
         ELSE IF(IF3D)THEN
            NLINES=3
            CALL PRS(' each of the 3 shear (TRaction) componnts$')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TRX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TRY=    Y  +  3 * TIME)$')
            CALL PRS(' TRZ=    Z  +  3 * TIME)$')
         ENDIF
      ELSE IF(cbcI.EQ.'shl')THEN
C        shear local
         IF(.NOT.IF3D)THEN
            NLINES=1
            CALL PRS(' shear (TRaction 1) $')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
         ELSE IF(IF3D)THEN
            NLINES=2
            CALL PRS(' each of the 2 shear (TRaction) componnts$')
            CALL PRS(' EXAMPLE:$')
            CALL PRS(' TR1=SIN(X) + EXP(-TIME)$')
            CALL PRS(' TR2=X      + EXP( TIME)$')
         ENDIF
      ELSE IF(cbcI.EQ.'sl ' .OR. cbcI .EQ. 'ms ')THEN
C        stress
         IF(.NOT.IF3D)NLINES=3
         IF(     IF3D)NLINES=4
         CALL PRSIS(' each of the $',NLINES,
     $   ' stress (TRACTION) components $')
         CALL PRS(' And Surface Tension$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TRN=SIN(X) + EXP(-TIME)$')
         CALL PRS(' TR1=SIN(X)$')
         IF(IF3D)CALL PRS(' TR2=1.0$')
         CALL PRS(' SIGMA=1.2 * TEMP$')
      ELSE IF(cbcI .EQ. 'msi')THEN
C        Fluid Layers
         NLINES = 1
         CALL PRS(' the surface tension coefficient$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' SIGMA = 1.2 * TEMP$')
      ELSE IF(cbcI .EQ. 'on ')THEN
C        Outflow
         NLINES = 1
         CALL PRS(' Enter exit pressure (TRactioN)function$')
         CALL PRS(' EXAMPLE:  TRN = 2. * TIME $')
      ELSE IF(cbcI .EQ. 'mf ')THEN
C        Melting Front
         NLINES = 2
         CALL PRS(' the freezing temp and the product Rho*Latent heat$')
         CALL PRS(' EXAMPLE:$')
         CALL PRS(' TEMP??? What is it Lee? = 1.2$')
         CALL PRS(' RHOL =What is it LEE??? 3 + $')
      ELSE IF(cbcI(1:1).EQ.'f')THEN
C        flux
         NLINES=1
         CALL PRS(' the flux.  EXAMPLE: FLUX=1.0$')
      ELSE IF(cbcI(1:1).EQ.'t')THEN
C        temperature
         NLINES=1
         CALL PRS(' TEMP.  Example:$')
         CALL PRS(' TEMP = 0.667 * (1.0-Y**2) +EXP(TIME)$')
         CALL PRS(' Enter statement For TEMP:$')
      ELSE IF(cbcI(1:1).EQ.'v'.OR.cbcI(1:1).EQ.'m')THEN
C        velocity
         IF(.NOT.IF3D)NLINES=2
         IF(     IF3D)NLINES=3
         CALL PRS(' each of the $',NLINES,' velocity components as $')
         CALL PRS(' functions of space and time.  Example:$')
         IF(cbcI .EQ. 'vl') THEN
            CALL PRS(' UN=SIN(X) + EXP(-TIME)$')
            CALL PRS(' U1=SIN(Y) + EXP(-TIME)$')
            IF(IF3D)
     $      CALL PRS(' U2 = 1 - TIME$')
         ELSE
            CALL PRS(' UX=SIN(X) + EXP(-TIME)$')
            CALL PRS(' UY=1$')
            IF(IF3D)
     $      CALL PRS(' UZ = 0.0$')
         ENDIF
      ELSE IF(cbcI(1:1).EQ.'c')THEN
C        convection
         NLINES=2
         CALL PRS(' Tinf & one line for HC.  Example:$')
         CALL PRS(' TINF=SIN(X) + EXP(TIME)$')
         CALL PRS(' HC=5.0$')
      ELSE IF(cbcI(1:1).EQ.'r')THEN
         CALL PRS(
     $   ' The HRAD (Sigma*F) and TINF (the absolute Temperature).$')
         CALL PRS('  EXAMPLE: HRAD = 1.0e-4$')
         NLINES = 2
      ELSE
         CALL PRS('New BC Type.$')
         CALL PRS('Enter BC in the Following three lines $')
         CALL PRS('Ignoring Prompts for components$')
         NLINES=3
      ENDIF
      IF(cbc3(3:3) .EQ. 'i')THEN
         BC(4,ISIDE,IEL,IF)=NLINES
         BC(5,ISIDE,IEL,IF)=LOCLIN
      ELSE
         BC(1,ISIDE,IEL,IF)=NLINES
         BC(2,ISIDE,IEL,IF)=LOCLIN
      ENDIF
      DO 100 I=1,NLINES
         IF(NLINES.GT.1 .AND. cbcI .eq. 'ms ')THEN
            IF(I.EQ.1)CALL PRS(' Normal Component:$')
            IF(IF3d)THEN
               IF(I.EQ.2)CALL PRS(' #1 Tangential Component:$')
               IF(I.EQ.3)CALL PRS(' #2 Tangential Component:$')
               IF(I.EQ.4)CALL PRS(' Surface Tension SIGMA:$')
            ELSE
               IF(I.EQ.2)CALL PRS(' #1 Tangential Component:$')
               IF(I.EQ.3)CALL PRS(' Surface Tension SIGMA:$')
            ENDIF
         ELSE IF(NLINES.GT.1 .AND. cbcI .eq. 'shl')THEN
            IF(I.EQ.1)CALL PRS(' 1st component:$')
            IF(I.EQ.2)CALL PRS(' 2nd component:$')
         ELSE IF(NLINES.NE.1 .AND.(cbcI(2:2).EQ.'l'
     $   .OR.cbcI(3:3).eq.'l'))THEN
            IF(I.EQ.1)CALL PRS(' Enter Normal component:$')
            IF(I.EQ.2.AND.NLINES.EQ.2)
     $      CALL PRS(' Enter Tangential component:$')
            IF(I.EQ.2.AND.NLINES.EQ.3)
     $      CALL PRS(' Enter Tangential component #1:$')
            IF(I.EQ.3)CALL PRS(' Enter Tangential component #2:$')
         ELSE IF(NLINES.NE.1 .AND.cbcI(1:1).EQ.'r'.AND.I.EQ.2)THEN
            CALL PRS(' Now TINF.   EXAMPLE: TINF = 1.0e-4$')
         ELSE IF(NLINES.NE.1 .AND.cbcI(1:1).NE.'c')THEN
            IF(I.EQ.1)CALL PRS(' Enter X-component:$')
            IF(I.EQ.2)CALL PRS(' Enter Y-component:$')
            IF(I.EQ.3)CALL PRS(' Enter Z-component:$')
         ENDIF
 99      CALL RES(INBC(LOCLIN),60)
         IF(IFLEARN)WRITE(3,'(A60)')INBC(loclin)
         IF(INBC(LOCLIN) .EQ. ' ')THEN
            CALL PRS('Please enter a non-blank fortran line:$')
            GO TO 99
         ENDIF
         LOCLIN=LOCLIN+1
100   CONTINUE
C
      IF(cbc3(3:3) .EQ. 'i')THEN
C        EXPORT B.C. FOR INTERNAL Boundary
         CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
         cbc(  ISIDEO,IELO,IF) =cbc(  ISIDE,IEL,IF)
         IF(cbc3 .EQ. 'msi') cbc(  ISIDEO,IELO,IF) = 'mpi'
         BC (4,ISIDEO,IELO,IF) = BC(4,ISIDE,IEL,IF)
         DO 121 I=LOCLIN-1,LOCLIN-1+NLINES
            INBC(I+NLINES)=INBC(I)
121      CONTINUE
          LOCLIN=LOCLIN+NLINES
          BC(5,ISIDEO,IELO,IF) = LOCLIN
      ENDIF

      RETURN
      END
C
c-----------------------------------------------------------------------
      subroutine letbc(f,e,ifld,bclab)
      include 'basics.inc'
      character bclab,bclab2*2
      integer e,f

      bclab2(1:1)=bclab
      bclab2(2:2)='$'
      yshift= 0.0
      if (f.eq.5) yshift=-yfac*0.05
      if (f.eq.6) yshift= yfac*0.02
      xshift= xfac*(-3.0/200.) +(ifld-1)*xfac/50.0
      xlab=sides(e,f,1)+xshift
      ylab=sides(e,f,2)+yshift
C     call prsii(' B.C.?$',e,f)

      if (bclab.eq.' ') then ! Draw empty box
         dx=xfac/50.
         dy=yfac/40.
         call fillp(-15)
         call beginp(xlab,ylab)
         call draw(xlab+dx,ylab   )
         call draw(xlab+dx,ylab+dy)
         call draw(xlab   ,ylab+dy)
         call draw(xlab   ,ylab   )
         call endp
      else ! Normal Write
         if(cwrite.ne.0. .and.bclab2 .ne.'E$')
     $   call gwrite(xlab+xfac/300.,ylab+yfac/300.,1.0,bclab2)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine chkbcs
      include 'basics.inc'
      character*3 cbctmp
      character*1 yesno
      logical iffaiL
C
      IF=1
      do 100 iel=1,nelv
         DO 100 ISIDE=1,NSIDES
            IF( cbc(ISIDE,IEL,IF).EQ.'MSI'
     $      .OR.cbc(ISIDE,IEL,IF).EQ.'msi')THEN
               CALL OVERLAP(IEL,ISIDE,IELO,ISIDEO)
               IF(cbc (ISIDEO,IELO,IF).NE. 'MPI'
     $        .and.cbc(ISIDEO,IELO,IF).NE. 'mpi')THEN
                 CALL BEEP
                 CALL PRS('Inconsistency in Fluid Layer B.C.s$')
                 CALL PRS('Attempting Correction....$')
                 IF(cbc(ISIDE,IEL,IF).EQ.'MSI')
     $              cbc(ISIDEO,IELO,IF)='MPI'
                 IF(cbc(ISIDE,IEL,IF).EQ.'msi')
     $              cbc(ISIDEO,IELO,IF)='mpi'
               ENDIF
               IF(IGROUP(IELO) .EQ. IGROUP(IEL))THEN
                  CALL BEEP
                  CALL PRS('Inconsistency in Fluid Layer B.C.s$')
                  CALL PRII(IEL,IELO)
                  CALL PRS
     $            ('Are 2 elements of same fluid that straddle a$')
                  CALL PRS
     $            ('Fluid Layer.  Different fluids (groups) are reqd$')
               ENDIF
               IF(IGROUP(IELO) .GT. IGROUP(IEL))THEN
C                 Swap
                  cbcTMP              = cbc(ISIDEO,IELO,IF)
                  cbc(ISIDEO,IELO,IF) = cbc(ISIDE ,IEL ,IF)
                  cbc(ISIDE ,IEL ,IF) = cbcTMP
                  CALL PRS('Normal Fluid layer B.C. swap$')
               ENDIF
            ENDIF
  100 CONTINUE
C
C     Check for curve side consistency
C
  101 CONTINUE
      if (ifflow) ifld=1
      if (ifheat) ifld=2
C
      iffail = .false.
      do 400 ie=1,nel
      do 400 iside=1,nsides
         IF (cbc(iside,ie,ifld).eq.'E  '   .or.
     $       cbc(iside,ie,ifld).eq.'P  ' )   then
 
             je    = bc(1,iside,ie,ifld)
             jside = bc(2,iside,ie,ifld)
             je    = ibc (iside,ie,ifld) ! high-precision periodic bc
 
c           Check edges
            IF (ISIDE.le.4) THEN
             Do 250 ilev=0,4,4
               ied = iside + ilev
               IF (CCURVE(ied,ie).eq.'C') THEN
                jed = jside+ilev
                IF (CCURVE(jed,je).ne.CCURVE(ied,ie)    .or.
     $             -CURVE(1,jed,je).ne.CURVE(1,ied,ie)) THEN
                 call qchk(ie,iside,je,jside,IF)
                 iffail = .true.
                 CALL PRS('ERROR: Curve side consistency failure A.$')
  200            CALL PRS(
     $          'Enter element number to be changed (other=abort).$')
                 WRITE(S,210) ie,ied,CCURVE(ied,ie),Curve(1,ied,ie)
                 CALL PRS(S)
                 WRITE(S,210) je,jed,CCURVE(jed,je),Curve(1,jed,je)
                 CALL PRS(S)
  210            FORMAT(
     $           ' El:',I10,3x,'Edge:',I2,3x,'C: ',A3,' Rad:',G14.6,'$')
cccc             CALL REI(Ke)
                 IF (Ke.eq.ie) THEN
                    CCURVE(ied,ie)   = CCURVE(jed,je)
                    CURVE (1,ied,ie) = -CURVE (1,jed,je)
                 ELSE IF(Ke.eq.je) THEN
                    CCURVE(jed,je)   = CCURVE(ied,ie)
                    CURVE (1,jed,je) = -CURVE (1,ied,ie)
                 ELSE
                  CALL PRS(
     $            'Are you sure you want to keep it this way (A)?$')
                  je =jed
c                 call prexit
cccc              CALL RES(YESNO,1)
cccc              IF (YESNO.ne.'y'.and.YESNO.ne.'Y') GOTO 200
                 ENDIF
                ENDIF
               ENDIF
  250        CONTINUE
            ENDIF
C
            if (ccurve(iside,ie).eq.'s') then
               tol   = 1.e-4
               diffc =
     $           dif_rel_inf_norm(curve(1,jside,je),curve(1,iside,ie),4)
               if (ccurve(jside,je).ne.ccurve(iside,ie)   .or.
     $            diffc.gt.tol) then
                 iffail = .true.
                 call qchk(ie,iside,je,jside,if)
                 call prs('ERROR: Curve side consistency failure B.$')
  300            call prs(
     $          'Enter element number to be changed (other=abort).$')
                 write(s,310) ie,iside,CCURVE(iside,ie),
     $                        (curve(j,iside,ie),j=1,4)
                 call prs(s)
                 write(s,310) je,jside,CCURVE(jside,je),
     $                        (curve(j,jside,je),j=1,4)
                 call prs(s)
  310            format(
     $           ' El:',I10,3x,'F:',I2,3x,'c: ',A3,4G13.5,'$')
cccc             CALL REI(Ke)
                 IF (Ke.eq.ie) THEN
                    CCURVE(iside,ie)   = CCURVE(jside,je)
                    CURVE (1,iside,ie) = CURVE (1,jside,je)
                    CURVE (2,iside,ie) = CURVE (2,jside,je)
                    CURVE (3,iside,ie) = CURVE (3,jside,je)
                    CURVE (4,iside,ie) = CURVE (4,jside,je)
                    CURVE (5,iside,ie) = CURVE (5,jside,je)
                 ELSE IF(Ke.eq.je) THEN
                    CCURVE(jside,je)   = CCURVE(iside,ie)
                    CURVE (1,jside,je) = CURVE (1,iside,ie)
                    CURVE (2,jside,je) = CURVE (2,iside,ie)
                    CURVE (3,jside,je) = CURVE (3,iside,ie)
                    CURVE (4,jside,je) = CURVE (4,iside,ie)
                    CURVE (5,jside,je) = CURVE (5,iside,ie)
                 ELSE
                  CALL PRS(
     $            'Are you sure you want to keep it this way? (B)$')
cccc              CALL RES(YESNO,1)
cccc              IF (YESNO.ne.'y'.and.YESNO.ne.'Y') GOTO 300
                 ENDIF
               ENDIF
            ENDIF
         ENDIF
  400 CONTINUE
      if (.not.iffail) CALL PRS('Curve side consistency check OK.$')
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine fndsida(kside,ke,iside,ie,if)
      include 'basics.inc'
      real vec1(3),vec2(3),vec3(3)
C
C     Find a side JSIDE,JEL which corresponds to ISIDE,IE and doesn't
C     have BC's set already
C
      NSIDES=NDIM*2
C
C     Define Normal to this plane (side), and find the element
C     center which also lies on this line, and is farthest from
C     the current element.
      CALL MKSIDE
      CALL GENCEN

      ifail = 0
      ifail = 2397


C     Find Sides' Normal Vector at midpoint

      IC1=ISIDE
      IC2=ISIDE+1
      IF(ISIDE.EQ.4)IC2=1
C     This stuff only relevant for 3d
      IC3=IC2+4
      IC4=IC1+4
      IF (ISIDE.EQ.5) THEN
         IC1=1
         IC2=2
         IC3=3
         IC4=4
      ELSE IF(ISIDE.EQ.6) THEN
         IC1=1+4
         IC2=2+4
         IC3=3+4
         IC4=4+4
      ENDIF
      write(s,9) ie,iside,ic1,ic2,ic3,ic4
    9 format(' iesc14:',6i5,'$')
      if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
      IF (IF3D) THEN
         vec1(1)=x(ic3,ie)-x(ic1,ie)
         vec1(2)=y(ic3,ie)-y(ic1,ie)
         vec1(3)=z(ic3,ie)-z(ic1,ie)
         vec2(1)=x(ic4,ie)-x(ic2,ie)
         vec2(2)=y(ic4,ie)-y(ic2,ie)
         vec2(3)=z(ic4,ie)-z(ic2,ie)
         CALL CROSS(VEC3,VEC1,VEC2,IE)
      ELSE
         VEC1(1)=X(IC2,IE)-X(IC1,IE)
         VEC1(2)=Y(IC2,IE)-Y(IC1,IE)
C
         VEC3(1)= VEC1(2)
         VEC3(2)=-VEC1(1)
         VEC3(3)=0.0
      ENDIF
      CALL NORM3D(VEC3)
c
      if (if_lattice) then          ! Make vector point in lattice direction
         d1 = dotprod(vec3,ulat1)
         d2 = dotprod(vec3,ulat2)
         d3 = dotprod(vec3,ulat3)
         if (d1.le.d2.and.d1.le.d3) then
            call copy(vec3,ulat1,3)
         elseif (d2.le.d3) then
            call copy(vec3,ulat2,3)
         else
            call copy(vec3,ulat3,3)
         endif
      endif
c
      if (ie.eq.ifail) write(6,19) vec3(1),vec3(2),vec3(3)
      write(s,19) vec3(1),vec3(2),vec3(3)
   19 format(' vec3:',3e13.5,'$')
      if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
c
c     EPS2=(EPSM*RCEN(IE))**2
      EPSM=1.e-4
      EPS2=EPSM**2
c
      KE=0
      KSIDE=0
      DSTMAX=0.0
c
      DO 100 JE=1,NEL
      DO 100 JSIDE=1,NSIDES
C
C        Don't find yourself.
         IF (JE.EQ.IE.AND.JSIDE.EQ.ISIDE) GOTO 100
C
C        Don't find a bc which is already defined. (to be 'E')
         IF (cbc(JSIDE,JE,IF).EQ.'E') GOTO 100
C
C        OK, is the center of this element side close to the line?
C
         VEC1(1)=SIDES(IE,ISIDE,1)-SIDES(JE,JSIDE,1)
         VEC1(2)=SIDES(IE,ISIDE,2)-SIDES(JE,JSIDE,2)
         VEC1(3)=SIDES(IE,ISIDE,3)-SIDES(JE,JSIDE,3)
         DISTP=DOTPROD(VEC1,VEC1)
         DIST2=DISTP-(DOTPROD(VEC1,VEC3))**2
         if (distp.gt.0) DIST2=DIST2/DISTP
         if (ie.eq.ifail) write(6,11) dist2,distp,eps2,dstmax,je,jside
   11    format(1p4e12.4,2i6,' dist2')
         IF (DIST2.LE.EPS2.AND.DISTP.GT.DSTMAX) THEN
            KE=JE
            KSIDE=JSIDE
            DSTMAX=DISTP
         ENDIF
  100 CONTINUE
      write(6,*) ie,ke,' pbc'
      return
      end
c-----------------------------------------------------------------------
      subroutine qchk(ie,iside,je,jside,ifld)
      include 'basics.inc'
      write(6,9) ie,iside,cbc(iside,ie,ifld),
     $            (bc(j,iside,ie,ifld),j=1,5)
      write(6,9) je,jside,cbc(jside,je,ifld),
     $            (bc(j,jside,je,ifld),j=1,5)
    9 format(' ie,side,c,b:',i5,i2,1x,a3,1x,5g11.3)
      return
      end
c-----------------------------------------------------------------------
      function i_periodic(iside,ie,ifld,bclab)
c
c     Assign periodic-auto bc
c
      include 'basics.inc'
      character*1 bclab(2)
c
      i_periodic=0
c
      call fndsida(jside,jel,iside,ie,ifld)
      if (jside.eq.0) then
        call prsii('For element,side:',ie,iside)
        call prs('Cannot find a corresponding side, choose a new BC.$')
        i_periodic=1
        return
      endif
c
      call letbc(jside,je,ifld,bclab)
      call letbc(iside,ie,ifld,bclab)
      cbc (iside,ie,ifld) = 'P  '
      ibc (iside,ie,ifld) = je       ! high-precision periodic bc
      bc(1,iside,ie,ifld) = je
      bc(2,iside,ie,ifld) = jside

      cbc (jside,je,ifld) = 'P  '
      ibc (jside,je,ifld) = ie       ! high-precision periodic bc
      bc(1,jside,je,ifld) = ie
      bc(2,jside,je,ifld) = iside
 
      return
      end
c-----------------------------------------------------------------------
      subroutine getnormals(unx,uny,unz,ie,is,x,y,z,ld,ndim)
      real x(ld,8),y(ld,8),z(ld,8)
c
      real u(3),v(3),w(3)
c
c     Base normals on cross-product of diagonals
c
      if (is.eq.1) then
         u(1) = x(6,ie)-x(1,ie)
         u(2) = y(6,ie)-y(1,ie)
         u(3) = z(6,ie)-z(1,ie)
         v(1) = x(5,ie)-x(2,ie)
         v(2) = y(5,ie)-y(2,ie)
         v(3) = z(5,ie)-z(2,ie)
      elseif(is.eq.2) then
         u(1) = x(7,ie)-x(2,ie)
         u(2) = y(7,ie)-y(2,ie)
         u(3) = z(7,ie)-z(2,ie)
         v(1) = x(6,ie)-x(3,ie)
         v(2) = y(6,ie)-y(3,ie)
         v(3) = z(6,ie)-z(3,ie)
      elseif(is.eq.3) then
         u(1) = x(8,ie)-x(3,ie)
         u(2) = y(8,ie)-y(3,ie)
         u(3) = z(8,ie)-z(3,ie)
         v(1) = x(7,ie)-x(4,ie)
         v(2) = y(7,ie)-y(4,ie)
         v(3) = z(7,ie)-z(4,ie)
      elseif(is.eq.4) then
         u(1) = x(5,ie)-x(4,ie)
         u(2) = y(5,ie)-y(4,ie)
         u(3) = z(5,ie)-z(4,ie)
         v(1) = x(8,ie)-x(1,ie)
         v(2) = y(8,ie)-y(1,ie)
         v(3) = z(8,ie)-z(1,ie)
      elseif(is.eq.5) then
         u(1) = x(4,ie)-x(2,ie)
         u(2) = y(4,ie)-y(2,ie)
         u(3) = z(4,ie)-z(2,ie)
         v(1) = x(3,ie)-x(1,ie)
         v(2) = y(3,ie)-y(1,ie)
         v(3) = z(3,ie)-z(1,ie)
      elseif(is.eq.6) then
         u(1) = x(7,ie)-x(5,ie)
         u(2) = y(7,ie)-y(5,ie)
         u(3) = z(7,ie)-z(5,ie)
         v(1) = x(8,ie)-x(6,ie)
         v(2) = y(8,ie)-y(6,ie)
         v(3) = z(8,ie)-z(6,ie)
      endif
      call vcross_normal(w,sine,u,v)
      unx = w(1)
      uny = w(2)
      unz = w(3)

      return
      end
c-----------------------------------------------------------------------
      subroutine vcross_normal(u,r2,v,w)
C
C     Compute a Cartesian vector cross product.
C
      real u(3),v(3),w(3)
c
      u(1) = v(2)*w(3) - v(3)*w(2)
      u(2) = v(3)*w(1) - v(1)*w(3)
      u(3) = v(1)*w(2) - v(2)*w(1)
c
c     Normalize
c
      r2 = u(1)*u(1) + u(2)*u(2) + u(3)*u(3)
      if (r2.gt.0.) then
         r2 = sqrt(r2)
         u(1) = u(1)/r2
         u(2) = u(2)/r2
         u(3) = u(3)/r2
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine out_cent
      include 'basics.inc'
c
c     Dump element centroids, w/ element numbers for positional sorting
c                                   pff 4/8/99
c
      return

      call gencen
c
      open(unit=53,file='elcent.dat')
      write(53,*) nel
      do ie=1,nel
         write(53,53) zcen(ie),ycen(ie),xcen(ie),ie
      enddo
   53 format(3f18.8,i9)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc1
      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for hemisphere/plate problem.  pff 12/22/98
c
c     1st, set periodic bcs in x-direction
c
      bclab(2)='$'
      nsides = 2*ndim
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ') then
            if (ccurve(is,ie).eq.'s') then
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'T  '
               bc( 1,is,ie,2) = 1.
            elseif (ccurve(is,ie).eq.'C'.and.is.le.4) then
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'T  '
               bc( 1,is,ie,2) = 1.
            else
c
c              Straight side --- which way is it pointing?
c
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unx).ge.0.9) then
                  ierr = i_periodic(is,ie,1,bclab)
                  ierr = i_periodic(is,ie,2,bclab)
               elseif (abs(uny).ge.0.9) then
                  cbc  (is,ie,1) = 'SYM'
                  cbc  (is,ie,2) = 'I  '
               elseif (unz.ge.0.9) then
                  cbc  (is,ie,1) = 'SYM'
                  cbc  (is,ie,2) = 'I  '
               elseif (unz.le.-0.9) then
                  cbc  (is,ie,1) = 'W  '
                  cbc  (is,ie,2) = 'T  '
                  bc( 1,is,ie,2) = 0.
               endif
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc2
      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for AV graft problem.  pff 6/23/01
c
      ntot = nx*ny*nz*nel
      nsides = ndim*2
      ncrnrs = ndim**2
c
      zmax = z(1,1)
      zmin = z(1,1)
      do ic=1,ncrnrs
      do ie=1,nel
         zmax = max(z(ic,ie),zmax)
         zmin = min(z(ic,ie),zmin)
      enddo
      enddo
      dz = 0.05*(zmax-zmin)
      zdmax = zmax - dz
      zdmin = zmin + dz
c
c
      bclab(2)='$'
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ') then
            cbc  (is,ie,1) = 'W  '
            cbc  (is,ie,2) = 'I  '
            call rzero(bc(1,is,ie,1),5)
            call rzero(bc(1,is,ie,2),5)
c
            if (zcen(ie).gt.zdmax) then
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unz).ge.0.9) then
                  cbc  (is,ie,1) = 'v  '
                  cbc  (is,ie,2) = 't  '
               endif
            elseif (zcen(ie).lt.zdmin) then
               call getnormals(unx,uny,unz,ie,is,x,y,z,nelm,ndim)
               if (abs(unz).ge.0.9) then
                  cbc  (is,ie,1) = 'O  '
                  cbc  (is,ie,2) = 'O  '
               endif
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine user_bc
c
c     Base boundary conditions on face number
c
c     grep " 3 " b.3 | wc 560          3360         41440
c     grep " 5 " b.3 | wc 96           576          7104
c     grep " 6 " b.3 | wc 48           288          3552
c

      include 'basics.inc'
      character*1 bclab(2)
c
c     set all non E-E bc's
c
c     This is currently set up for AV graft problem.  pff 6/23/01
c
      ntot = nx*ny*nz*nel
      nsides = ndim*2
      ncrnrs = ndim**2
c
      zmax = z(1,1)
      zmin = z(1,1)
      do ic=1,ncrnrs
      do ie=1,nel
         zmax = max(z(ic,ie),zmax)
         zmin = min(z(ic,ie),zmin)
      enddo
      enddo
      dz = 0.05*(zmax-zmin)
      zdmax = zmax - dz
      zdmin = zmin + dz
c
c
      bclab(2)='$'
      do ie=1,nel
      do is=1,nsides
         if (cbc(is,ie,1).ne.'E  ' .and.
     $       cbc(is,ie,1).ne.'J  ' .and.
     $       cbc(is,ie,1).ne.'SP ') then
            call rzero(bc(1,is,ie,1),5)
            call rzero(bc(1,is,ie,2),5)
            if (is.eq.6) then
               cbc  (is,ie,1) = 'v  '
               cbc  (is,ie,2) = 't  '
            elseif (is.eq.5) then
               cbc  (is,ie,1) = 'O  '
               cbc  (is,ie,2) = 'O  '
            else
               cbc  (is,ie,1) = 'W  '
               cbc  (is,ie,2) = 'I  '
            endif
         endif
c        write(6,1) is,ie,cbc(is,ie,1),unx,uny,unz
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine set_match_p(je,jf,ie,if)
      include 'basics.inc'

      if1=1
      if(.not. ifflow) if1=2

      ifld0=if1
      ifld1=nflds

      do ifld=ifld0,ifld1

         cbc (  if,ie,ifld) = 'P  '
         ibc (  if,ie,ifld) = je    ! high-precision pointer
         bc  (1,if,ie,ifld) = je
         bc  (2,if,ie,ifld) = jf

         cbc (  jf,je,ifld) = 'P  '
         ibc (  jf,je,ifld) = ie    ! high-precision pointer
         bc  (1,jf,je,ifld) = ie
         bc  (2,jf,je,ifld) = if

      enddo


      return
      end
c-----------------------------------------------------------------------
      function get_dist_ctr(x,y,z,f,ndim) ! dist to centroid

      real x(3,3,3),y(3,3,3),z(3,3,3)
      integer f

      nz = 1 + 2*(ndim-2)  ! 1 or 3
      call facind(i0,i1,j0,j1,k0,k1,3,3,nz,f)

      d2m = 1.e22
      do ipass=1,2
         l=0
         do k=k0,k1
         do j=j0,j1
         do i=i0,i1
            l=l+1
            if (ipass.eq.1.and.l.eq.5) then
               xc=x(i,j,k)
               yc=y(i,j,k)
               zc=z(i,j,k)
            elseif (ipass.eq.2.and.l.ne.5) then
               dx = xc-x(i,j,k)
               dy = yc-y(i,j,k)
               dz = zc-z(i,j,k)
               d2 = dx*dx + dy*dy + dz*dz
               d2m = min(d2,d2m)
            endif
         enddo
         enddo
         enddo

      enddo

      if (d2m.gt.0) d2m = sqrt(d2m)
      get_dist_ctr = d2m

      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiodz(dir,tol,cb_clear,xyz,m) ! dir=X,Y, or Z.
      include 'basics.inc'
      common /domainr/ scal(nelm),list_e(2*nelm)
      character*1 dir
      character*3 cb,ck,cb_clear
      real xyz(27,nelm)
      integer e,f,en,ep,fn,fp

c
c     Set any remaining candidate bcs to P
c
      n = 27*nel
      xmax=glmax(xyz,n)
      xmin=glmin(xyz,n)
      dx  =xmax-xmin

      nmin=0
      nmax=0
      dxmn=0
      dxmx=0
      do e=1,nel
         xmn=vlmin(xyz(1,e),27)
         xmx=vlmax(xyz(1,e),27)
         xnx=xmx-xmn
         if (abs(xmn-xmin).le.tol*dx) then
            dxmn = dxmn + xnx
            nmin = nmin + 1
         endif
         if (abs(xmx-xmax).le.tol*dx) then
            dxmx = dxmx + xnx
            nmax = nmax + 1
         endif
      enddo
      if (nmin.gt.0) dxmn=dxmn/nmin ! avg element thickness at bottom
      if (nmax.gt.0) dxmx=dxmx/nmax ! avg element thickness at top
      toln = 0.1*dxmn
      tolx = 0.1*dxmx

      n_ext=0                 ! Number of extrema
      call izero(list_e,nel)   ! List of extrema
      do ipass=1,2
        do e=1,nel
         xmn=vlmin(xyz(1,e),27)
         xmx=vlmax(xyz(1,e),27)
         if (ipass.eq.1.and.abs(xmn-xmin).le.toln) then
            n_ext = n_ext + 1
            list_e(n_ext) = -e
         elseif (ipass.eq.2.and.abs(xmx-xmax).le.tolx) then
            n_ext = n_ext + 1
            list_e(n_ext) = e
         endif
        enddo
        if (ipass.eq.1) nneg = n_ext
      enddo
      npos = n_ext-nneg

      do ie=1,n_ext           ! get dist to centroid
         e=abs(list_e(ie))
         f=6
         if (list_e(ie).lt.0) f=5
         scal(ie) = get_dist_ctr(x27(1,e),y27(1,e),z27(1,e),f,ndim)
      enddo
      
      
      do in=1,nneg

         en = abs(list_e(in))

         fn = 5
         jn = 5
         fp = 6
         jp = 23


         dxymn = 1.e20
         ipm   = 0
         do ip=nneg+1,n_ext
            ep=list_e(ip)
            scale = min(scal(ep),scal(en))
            dxy = (x27(jp,ep)-x27(jn,en))**2+(y27(jp,ep)-y27(jn,en))**2
            if (dxy.lt.dxymn) then
               ipm=ip
               dxymn=dxy
            endif
         enddo

         ip = ipm
         ep = list_e(ip)

         scale = .01*min(scal(ip),scal(in))
         if (dxymn.lt.scale*scale) call set_match_p(en,fn,ep,fp)

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         if (nel.gt.0) then                     ! Start Diagnostics
           if (dxymn.gt.0) dxymn = sqrt(dxymn)
           if (dxymn.lt.scale) then
              write(77,77) en,fn,ep,fp,x27(jn,en),y27(jn,en),z27(jn,en)
     $                                ,x27(jp,ep),y27(jp,ep),z27(jp,ep)
     $                                ,dxymn
           else
              write(78,77) en,fn,ep,fp,x27(jn,en),y27(jn,en),z27(jn,en)
     $                                ,x27(jp,ep),y27(jp,ep),z27(jp,ep)
     $                                ,dxymn
              list_e(1) = 1/ep
           endif
   77      format(i8,i2,i8,i2,1p7e10.2)
         endif                                  ! End Diagnostics
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine bound ! Set boundary conditions
      include 'basics.inc'
      character heatbc,velbc,key,string*5,ele(4),bclab(2),reply,ques(2)
     $         ,reply2,ipapt,bcchoice*26,mode*15,cmenu*15,cc*3
      logical iftmp
      integer icalld,e,f
      save    icalld
      data icalld / 0 /
      data cmenu /' '/

c     call out_cent ! dump element centroids, w/ element numbers for positional sorting

      call semi_init_bcs(maxlev,needbc,if1)

      call first_bc_menu(needbc,if1)

      write(6,*) 'bound A: ',choice,needbc,nlevel

      if (choice.eq.'ACCEPT B.C.''s') return

      call prs('                   *** BOUNDARY CONDITION MENU ***$')
      if (if3d) then !     Set up gray elevator
         if (maxlev.ge.2) then
            do 50 i=maxlev,2,-1
               call drelev(i-1,i,'     ')
50          continue
         else
            call drelev(1,0,'     ')
         endif
      endif

      do ifld=if1,nflds 
         call set_bc_fld(ifld)  ! Set BCs for each field
      enddo

      do if=if1,nflds
         call period_check(if) ! Check and reset stray period bcs  , pff  4/7/99
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine semi_init_bcs(maxlev,needbc,if1)
      include 'basics.inc'
      character heatbc,velbc,key,string*5,ele(4),bclab(2),reply,ques(2)
     $         ,reply2,ipapt,bcchoice*26,mode*15,cmenu*15,cc*3
      logical iftmp
      integer fcorns (4,6),icalld,e,f
      save    fcorns,icalld
      data icalld / 0 /
      data cmenu /' '/

      ifgrid=.false. ! Turn off grid latch for improved mouse resolution.

      needbc=0
      maxlev=0

      if1=1
      if (.not. ifflow) if1=2

      ifld0=if1
      ifld1=nflds
      if (ifconj_merge) ifld0=2
      if (ifconj_merge) ifld1=2

      if (nelv.ne.nel) then
         do ifld=ifld0,ifld1 ! Blank out B.C.'s for nonexistent elements.
         do e=nelv+1,nel
         do f=1,nsides
            if(.not.iftmsh(ifld)) cbc(f,e,ifld)='   '
         enddo
         enddo
         enddo
      endif

      do ifld=ifld0,ifld1
        if (     iftmsh(ifld)) nnel=nel
        if (.not.iftmsh(ifld)) nnel=nelv
        do e=1,nnel
           maxlev=max(numapt(e),maxlev)
         enddo
      enddo

      do ifld=ifld0,ifld1
        if (     iftmsh(ifld)) nnel=nel
        if (.not.iftmsh(ifld)) nnel=nelv
        do e=1,nnel
        do f=1,nsides
             if (cbc(f,e,ifld).eq.'   ') then
                needbc=1
                goto 11
             endif
         enddo
         enddo
      enddo
   11 continue

      if (nlevel.eq.1) then ! Display B.C.'s already here (for 2-d case)
       do ifld=ifld0,ifld1
       do e=1,nel
       do f=1,nsides
         if (cbc(f,e,ifld).ne.'   ') call letbc(f,e,ifld,cbc(f,e,ifld))
       enddo
       enddo
       enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine first_bc_menu(needbc,if1)

      include 'basics.inc'

      if (needbc.eq.1) then ! Ya gotta have b.c.'s
         call prs('              *** MIDWAY BREAK MENU ***$')
         call prs('ABORT to write build data to file and exit$')
         call prs('SET B.C''S to continue$')
         item(1)='SET BCs'
         item(2)='ABORT'
         nchoic =  2
         call menu(xmouse,ymouse,button,'MIDWAY BREAK')
         if (choice.eq.'ABORT') then
             call mesgen
             call prexit
         endif
      else
         item(1)='ACCEPT B.C.''s'
         item(2)='REVIEW/MODIFY'
         nchoic = 2
         if (ifconj_merge) then
            choice = item(1)
         else
            call menu(xmouse,ymouse,button,'ACCEPT/REVIEW')
         endif
         if (choice.eq.'ACCEPT B.C.''s') then

            call chkbcs
 
            do ifld=if1,nflds
               call period_check(ifld)
            enddo
 
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rehighlight(e,f,jlevel)
      include 'basics.inc'
      integer e,f
c
c     Re-highlight the current element boundaries
c
      call color(2)
      if (f.eq.5.or.f.eq.6) then
         if (f.eq.5) ilev=jlevel
         if (f.eq.6) ilev=jlevel+1
         call movesc(elec-elew,eleb+eledy*(ilev-1))
         call drawsc(elec+elew,eleb+eledy*(ilev-1))

c        Flash mesh sides
         if (f.eq.5) then
            ied1=1
            ied2=4
         elseif (f.eq.6) then
            ied1=1
            ied2=4
         endif
         call movec(x(ied1,e),y(ied1,e))
         do iedge=ied1,ied2
            call drawed(e,iedge,1)
         enddo
      else
         call movec(x(f,e),y(f,e)) ! sides 1-4  (Same as edges in this case)
         call drawed(e,f,1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine rezoom
      include 'basics.inc'
      integer e

      call setzoom
      call refresh
      call drmenu('NOCOVER')
      call drgrid

      do e=1,nel
         call drawel(e)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine black_out_level
      include 'basics.inc'

      call color(0)
      call fillp(0)
      call beginp(xphy(0.0),yphy(0.0))
      call draw  (xphy(.90),yphy(0.0))
      call draw  (xphy(.90),yphy(.99))
      call draw  (xphy(0.0),yphy(.99))
      call draw  (xphy(0.0),yphy(0.0))
      call endp

      return
      end
c-----------------------------------------------------------------------
      subroutine set_up_level(num_on_this_level,ifld,jlevel)
      include 'basics.inc'
      integer e,f

      call black_out_level  ! Black out whole level

      if (icalld.eq.0) then ! Make isometric drawing
         do e=1,nel
            call drawis(isrt(e))
         enddo
         icalld=1
      endif

c     Draw mesh OR WHOLE ELEMENT?  ??Delete previous mesh(+bclabels?)
      call color(1)
      num_on_this_level=0
      do 300 e=1,nel
         if (numapt(e).eq.jlevel) then
            num_on_this_level=num_on_this_level+1
            call drawel(e)
            do 200 f=1,nsides
               if (cbc(f,e,ifld).ne.'   ')
     $         call letbc(f,e,ifld,cbc(f,e,ifld))
200         continue
         endif
300   continue
      if (ndim.eq.3) call drelev(jlevel,jlevel-1,'     ')

      return
      end
c-----------------------------------------------------------------------
      subroutine check_end_level(icontinue,jlevel,ifld)
      include 'basics.inc'
      character heatbc,velbc,key,string*5,ele(4),bclab(2),reply,ques(2)
     $         ,reply2,ipapt,bcchoice*26,mode*15,cmenu*15,cc*3
      logical iftmp
      integer e,f,ee,ff,edup,fdup,ep,fp,eper,fper

      do 1500 e=1,nel
         if (numapt(e).eq.jlevel) then  ! Not 100% clear we need this check
             itype=0
             do 1200 f=1,nsides
                cc = cbc(f,e,ifld)
                if (cc.eq.'VL'.or.cc.eq.'vl'.or.
     $              cc.eq.'VL'.or.cc.eq.'vl') then ! at least one BC
                   if (itype.eq.1) then            ! in local coords.
                      call prs('Error: You cannot mix local and'
     $                   //' cartesian coordinates in B.C.''s in '
     $                   //'the same element.$')
                      write(s,'(1X,A7,I11,A2,I3,A2,A2)')
     $                   'Element',e,' [',NUMAPT(e),LETAPT(e),']$'
                      call prs(s)
                      call prs(' has such a mixture.$')
                      call prs('Please MODIFY velocity and stress'
     $                 //' B.C.''s to consistent coordinates.$')
                      goto 1010
                   else
                      itype=2
                   endif
                endif
                if (cc.eq.'V'.or.cc.eq.'v'.or.
     $              cc.eq.'V'.or.cc.eq.'v') then ! at least one BC global
                   if (itype.eq.2) then ! There is a conflict
                      call prs('Error: You cannot mix local and'
     $                   //' cartesian coordinates in B.C.''s in '
     $                   //'the same element.$')
                      write(s,'(1x,a7,i11,a2,i3,a2,a2)')
     $                   'Element',e,' [',NUMAPT(e),
     $                   LETAPT(e),']$'
                      call prs(s)
                      call prs('has such a mixture.$')
                      call prs('Please MODIFY velocity '
     $                   //'and stress B.C.''s to consistent '
     $                   //'coordinates.$')
                      goto 1010
                   else
                      itype=1
                   endif
                endif
1200         continue
         endif
1500  continue
      icontinue = 0
      return

1010  continue
      icontinue = 1
      return

      return
      end
c-----------------------------------------------------------------------
      subroutine flash_side_a (ifade,iflash,e,f,fcorns,jlevel,ifld)
      include 'basics.inc'
      integer e,f
      integer fcorns (4,6)

      character*2 ques
      save ques
      data ques /'?$'/


      ff=f+1 ! Side which is on Boundary and Not periodic with Previous el
      if (ff.eq.5) ff = 1
      iflash=jlevel
      if (f.eq.5) iflash= -jlevel
      if (f.eq.6) iflash=-(jlevel+1)
      ifade=iflash
      call drelev(iflash,ifade,'BLINK')

      call color(2)  !  Flash element side ;  Hilight side with Blinking lights
      if (ifbwgks) call color(0)
      call hilites(e,f)
      call letbc(f,e,ifld,ques)
     
      if (if3d) then ! FLASH SIDE ON ISOMETRIC
         do icorn=0,4
            if (icorn.eq.0) ic=fcorns(4    ,f)
            if (icorn.ne.0) ic=fcorns(icorn,f)
            xisom=xphy(xscr(x(ic,e))/5.0 + 0.8)-z(ic,e)/20.
            yisom=yphy(yscr(y(ic,e))/5.0 + 0.8)-z(ic,e)/20.
            if (icorn.eq.0) call movec(xisom,yisom)
            if (icorn.ne.0) call drawc(xisom,yisom)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine flash_element_b(ifade,iflash,e,f,jlevel,ifld,fcorns)
      include 'basics.inc'
      integer e,f
      integer fcorns (4,6)
      character bclab(2)

      call color(3)
      ff=f+1
      if (ff.eq.5) ff = 1
      call drelev(iflash,ifade,'     ') ! Turn Off Der Blinkenlights

      if (f.eq.5.or.f.eq.6) then ! Un-Flash Elevator at element floor or ceiling
         if (f.eq.5) ilev=jlevel
         if (f.eq.6) ilev=jlevel+1
         call movesc(elec-elew,eleb+eledy*(ilev-1))
         call drawsc(elec+elew,eleb+eledy*(ilev-1))
         if (f.eq.5) then ! Un-Flash mesh sides
            ied1=1
            ied2=4
         elseif (f.eq.6) then
            ied1=1
            ied2=4
         endif
         call movec(x(ied1,e),y(ied1,e))
         do 230 iedge=ied1,ied2
            call drawed(e,iedge,1)
230      continue
      else
         if (if3d) then ! Un-Flash at elevator at walls of current level
            call movesc (elec,eleb+eledy*(jlevel-1))
            call drawsc (elec,eleb+eledy*(jlevel  ))
         endif
         call movec(x( f,e),y( f,e)) ! Un-Flash mesh side
c        call draw (x(ff,e),y(ff,e))
         call drawed(e,f,1)
      endif
      if (if3d) then !  un-flash side on isometric
         do 233 icorn=0,4
             if (icorn.eq.0) ic=fcorns(4    ,f)
             if (icorn.ne.0) ic=fcorns(icorn,f)
             xisom=xphy(xscr(x(ic,e))/5.0 + 0.8)-z(ic,e)/20.
             yisom=yphy(yscr(y(ic,e))/5.0 + 0.8)-z(ic,e)/20.
             if (icorn.eq.0) call movec(xisom,yisom)
             if (icorn.ne.0) call drawc(xisom,yisom)
233      continue
      endif

      call letbc(f,e,ifld,bclab) ! Label Boundary Side

      return
      end
c-----------------------------------------------------------------------
      subroutine set_bc_fld(ifld)  ! Set BCs for each field
      include 'basics.inc'
      character heatbc,velbc,key,string*5,ele(4),bclab(2),reply,ques(2)
     $         ,reply2,ipapt,bcchoice*26,mode*15,cmenu*15,cc*3
      logical iftmp
      integer e,f,ee,ff,edup,fdup,ep,fp,eper,fper,eo,fo
      integer fcorns (4,6),icalld
      save    fcorns,icalld
      data icalld / 0 /
      data cmenu /' '/
      data fcorns / 1,2,6,5,
     $              2,3,7,6,
     $              3,4,8,7,
     $              4,1,5,8,
     $              1,2,3,4,
     $              5,6,7,8 /

      bclab(2)='$'

      nsides=4
      if (if3d) nsides=6

      write(6,*) ilevel,nlevel,num_on_this_level,ifld,' level A'

      if (ifld.eq.1) then
        call prs('Begin Inputting Fluid Boundary Conditions$')
      elseif (ifld.eq.2) then
        call prs('Begin Inputting Thermal Boundary Conditions$')
      else
        call prs('Begin Inputting PASSIVE SCALAR Boundary Conditions$')
        nn=ifld-2
        call prsi('PASSIVE SCALAR #$',nn)
      endif
c
      do 1000 ilevel=1,nlevel
1010     continue
         if (nlevel.gt.1) then
            call set_up_level(num_on_this_level,ifld,ilevel)
            write(6,*) ilevel,nlevel,num_on_this_level,ifld,' level'
            if (num_on_this_level.eq.0) goto 1999
         endif

800      continue !    Go thru normal stuff
         mode='NORMAL'
         do 900 f=1,nsides
         do 900 e=1,nel
            IF (numapt(e) .ne. ilevel )   goto 900
            if (cbc(f,e,ifld).ne.'   ')     goto 900
            if (maskel(e,ifld) .eq. 0   )   goto 900

            call flash_side_a (ifade,iflash,e,f,fcorns,jlevel,ifld)

113         continue
c           Error Reading Input
            FLUX=0.0
            TEMP=0.0
            VX  =0.0
            VY  =0.0
            if (mode.NE.'SKIP MENU') call prsii(' B.C.>$',e,f)
            IDUP=0
c           SET UP MENU choiceS
            if (ifld.eq.1) then !  FLUID B.C.'s
               item( 1)='PERIODIC-AUTO'
               item( 2)='PERIODIC'
               item( 3)='WALL'
               item( 4)='VELOCITY'
               item( 5)='OUTFLOW'
               item( 6)='OUTFLOW/N'
               item( 7)='VELOCITY (LOCAL)'
               item( 8)='SYMMETRY'
               item( 9)='MOVING WALL'
               item(10)='MELTING FRONT'
c              item(11)='QUASI-SYMMETRY'
               nchoic=10

               if (IFAXIS) then
                  nchoic=nchoic+1
                  item(nchoic)='AXIS'
               endif
               if (IFSTRS) then
                  nchoic=nchoic+1
                  item(nchoic)='STRESS B.C.s'
               endif
            elseif (ifld.GE.2) then ! THERMAL BCs:  USE FOR PASSICVE SCALARS
               item(1)='PERIODIC-AUTO'
               item(2)='PERIODIC'
               item(3)='INSULATED'
               item(4)='FLUX'
               item(5)='TEMP' !  Convection IS now ready for v2
               item(6)='CONVECTION'
               nchoic=6
               if (IFFLOW.or.IFADVC(ifld)) then
                  nchoic=nchoic+1
                  item(nchoic)='OUTFLOW'
               endif
               if (IFAXIS) then
                  nchoic=nchoic+1
                  item(nchoic)='AXIS'
               endif
c              nchoic=nchoic+1
c              item(nchoic)='RADIATION'
            endif
            nchoic=nchoic+1
            item(nchoic)='SET ENTIRE LEVEL'
            nchoic=nchoic+1
            item(nchoic)='ZOOM'

            if (CMENU.eq.'INTERNAL') then
c              Start over with short menu
               if (ifld.eq.1) then
                  nchoic=1
                  item(nchoic)='FLUID LAYERS'
               elseif (ifld.eq.2) then
c                  Fill up solid & liquid front when melting front is
c                  chosen
c                  nchoic       =1
c                  item(nchoic)='SOLID FRONT'
c                  nchoic=nchoic+1
c                  item(nchoic)='LIQUID FRONT'
c
c                  nchoic=nchoic+1   No internal radiation yet
c                  item(nchoic)='RADIATION'
                  nchoic=1
                  item(nchoic)='MELTING FRONT'
c                  call prs('No internal boundaries for temperature$')
c                  goto 313
               elseif (ifld.GT.2) then
                  call prs(
     $               'No internal boundaries for passive scalars$')
                  goto 313
               endif
            endif

c           We already got info for first value; skip menu now.
            if (mode.eq.'GET FIRST') mode='SKIP MENU'
63          if (mode.NE.'SKIP MENU') then
               if (ifaxis.and.abs(sides(e,f,2)).lt.yfac/1000.) then
                  call prs('**Axis B.C. Strongly Recommended here**$')
                  call beep
               endif
               call menu (xmouse,ymouse,button,'NOCOVER')
            endif

            if (choice.eq.'STRESS B.C.s') then ! New B.C.'s From Lee Ho
               item(1)='MAIN B.C. MENU'
               item(2)='SHEAR'
               item(3)='SHEAR    (LOCAL)'
               item(4)='FREE SURFACE'
               item(5)='TRACTION'
               item(6)='TRACTION (LOCAL)'
c              item(7)='FLUID LAYERS' ! Catch at internal b.c. menu
               nchoic=7
               call menu (xmouse,ymouse,button,'NOCOVER')
               if (choice.eq.'MAIN B.C. MENU') goto 113
            endif
c           End of patch

           if (mode.eq.'GET FIRST'.AND.choice.eq.'PERIODIC') then
              call prs('*** ERROR ***  YOU CANNOT SET WHOLE LEVEL$')
              call prs('TO BE PERIODIC; CHOOSE A DIFFERENT B.C.$')
               goto 63
            elseif (mode.eq.'GET FIRST'.AND.
     $             choice.eq.'PERIODIC-AUTO') then
               call autoperiod
            endif

            if(choice.eq.'ZOOM') then
               call rezoom
               call rehighlight(e,f,ilevel)
               goto 63
            endif

            if(choice.eq.'SET ENTIRE LEVEL') then ! Get b.c., then copy for rest of loops
               mode='GET FIRST'
               call blank(choice,26)
               edup=e
               fdup=f
               call prs('Choose B.C.''s for all remaining sides '//
     $            'in entire level$')
               goto 63
            endif

            if (mode.eq.'SKIP MENU')BUTTON='RIGHT'
            if (BUTTON.eq.'RIGHT') then
               bcchoice='DUPLICATE'
               if (XSCR(XMOUSE).GE.1.0.AND.mode.NE.'SKIP MENU') then
                  call prs('Right button not used in red menu area.$')
                  call prs('Try again.$')
                  if (mode.NE.'NORMAL')call prs(
     $               'Global mode reset.  Start over.$')
                  mode='NORMAL'
                  goto 113
               endif
               if (mode.NE.'SKIP MENU') then
                 xp=xmouse
                 yp=ymouse
                 rmin=10000.0 ! Find closest element, side that we want to duplicate
                 do ee=1,nel 
                    if (numapt(ee).eq.ilevel) then ! Only those on same floor
                       do ff=1,nsides
                          r = ((xp-sides(ee,ff,1))**2
     $                      +  (yp-sides(ee,ff,2))**2)
                          if (r.lt.rmin) then
                             rmin = r
                             fdup = ff
                             edup = ee
                          endif
                       enddo
                    endif
                 enddo
                 if (fdup.eq.5.or.fdup.eq.6) then ! check IF it's above or below center
                    fdup=5
                    if (YP.GT.SIDES(edup,fdup,2)) fdup=6
                 endif
               endif
               if (cbc(fdup,edup,ifld).eq.'   ') then
                  call prsi ('B.C. not set for fdup$',fdup,
     $                   'of element$',edup)
                  goto 113
               endif
               cbc3=cbc(fdup,edup,ifld)
               if (cbc3.eq.'E  ') then
                  write(s,'(A5,I4,A12,i11,A1)')'Side ',fdup,
     $               ' of element ',edup,'$'
                  call prs(s)
                  call prs('is internal side.  I won''t duplicate.$')
                  goto 113
               elseif  (cbc3.eq.'P  ') then
                  write(s,'(A5,I4,A12,i11,A1)')'Side ',fdup,
     $               ' of element ',edup,'$'
                  call prs(s)
                  call prs('is Periodic side.  I won''t duplicate.$')
                  goto 113
               elseif (cbc3(3:3).eq.'i'.or.cbc3(3:3).eq.'I') then
                  write(s,'(A5,I4,A12,i11,A1)')'Side ',fdup,
     $               ' of element ',edup,'$'
                  call prs(s)
                  call prs('is Internal side.  I won''t duplicate.$')
                  goto 113
               endif

c              Need flag to know if b.c. is fortran function.
c              Special handling puts new fortran in next available lines
               cbc3 = cbc(fdup,edup,ifld)
               icbc = ICHAR(cbc3(1:1))
               if (icbc.GE.97.AND.icbc.LE.122) then ! lower case signifies function
                  if (cbc3(3:3).eq.'i  ') then ! Special storage locations
                     nlines=bc(4,fdup,edup,ifld)
                     line1 =bc(5,fdup,edup,ifld)
                     BC(4,f,e,ifld)=nlines
                     BC(5,f,e,ifld)=loclin
                  else
                     nlines=bc(1,fdup,edup,ifld)
                     line1 =bc(2,fdup,edup,ifld)
                     BC(1,f,e,ifld)=nlines
                     BC(2,f,e,ifld)=loclin
                  endif
                  cbc (f,e,ifld)=cbc( fdup,edup,ifld)
                  do 777 i=1,nlines
                     inbc(loclin)=inbc(line1+i-1)
                     loclin=loclin+1
777               continue
               else
                  call copy(bc(1,f,e,ifld),bc(1,fdup,edup,ifld),5)
                  cbc (f,e,ifld)=cbc( fdup,edup,ifld)
                  ibc (f,e,ifld)=ibc( fdup,edup,ifld)
               endif
               bclab(1) = cbc(fdup,edup,ifld)
            elseif (BUTTON.eq.'LEFT') then
               bcchoice=choice
               bclab(1)      = bcchoice(1:1)
               cbc(f,e,ifld) = bcchoice(1:1)
               if (choice.eq.'VELOCITY'        ) cbc(f,e,ifld)='V  '
               if (choice.eq.'VELOCITY (LOCAL)') cbc(f,e,ifld)='VL '
               if (choice.eq.'MOVING WALL     ') cbc(f,e,ifld)='mv '
               if (choice.eq.'STRESS'          ) cbc(f,e,ifld)='S  '
               if (choice.eq.'STRESS   (LOCAL)') cbc(f,e,ifld)='SL '
               if (choice.eq.'SHEAR'           ) cbc(f,e,ifld)='SH '
               if (choice.eq.'SHEAR    (LOCAL)') cbc(f,e,ifld)='SHL'
               if (choice.eq.'SYMMETRY'        ) cbc(f,e,ifld)='SYM'
               if (choice.eq.'QUASI-SYMMETRY'  ) cbc(f,e,ifld)='QSM'
               if (choice.eq.'FREE SURFACE'    ) cbc(f,e,ifld)='MS '
               if (choice.eq.'OUTFLOW/N'       ) cbc(f,e,ifld)='ON '
               if (choice.eq.'TRACTION'        ) cbc(f,e,ifld)='S  '
               if (choice.eq.'TRACTION (LOCAL)') cbc(f,e,ifld)='SL '
               if (choice.eq.'FLUID LAYERS'    ) cbc(f,e,ifld)='MSI'
               if (choice.eq.'MELTING FRONT'   ) cbc(f,e,ifld)='MF '
               if (choice.eq.'SOLID FRONT'     ) cbc(f,e,ifld)='MCI'
               if (choice.eq.'LIQUID FRONT'    ) cbc(f,e,ifld)='MLI'

              cbc3=cbc(f,e,ifld)
              bclab(1)= cbc3(1:1)

              if (choice.eq.'TEMP'         .or.
     $            choice.eq.'FLUX'            .or.
     $            choice.eq.'VELOCITY'        .or.
     $            choice.eq.'VELOCITY (LOCAL)'.or.
     $            choice.eq.'STRESS'          .or.
     $            choice.eq.'STRESS   (LOCAL)'.or.
     $            choice.eq.'SHEAR'           .or.
     $            choice.eq.'SHEAR    (LOCAL)'.or.
     $            choice.eq.'MOVING WALL     '.or.
     $            choice.eq.'CONVECTION'      .or.
     $            choice.eq.'OUTFLOW'         .or.
     $            choice.eq.'OUTFLOW/N'       .or.
     $            choice.eq.'TRACTION'        .or.
     $            choice.eq.'TRACTION (LOCAL)'.or.
     $            choice.eq.'RADIATION'       .or.
     $            choice.eq.'FREE SURFACE'    .or.
     $            choice.eq.'FLUID LAYERS'    .or.
     $            choice.eq.'MELTING FRONT'   .or.
     $            choice.eq.'SOLID FRONT'     .or.
     $            choice.eq.'LIQUID FRONT'    ) then
                  if (choice.eq.'MOVING WALL') then   ! only have fortran function
                     choice='FORTRAN FUNCTION'
                  elseif (choice.eq.'OUTFLOW') then   ! only have constant
                     choice='CONSTANT'
                  elseif (choice.eq.'OUTFLOW/N') then ! only have constant
                     choice='CONSTANT'
                  elseif (choice.eq.'MELTING FRONT') then ! only have constant
                     choice='CONSTANT'
                  else
                     call prs(' Choose Format for Relevant Parameters$')
                     nchoic=2
                     item(1)='CONSTANT'
                     item(2)='FORTRAN FUNCTION'
                     call menu (xmouse,ymouse,button,'FORTRAN FUNCTION')
                  endif
                  if (choice.eq.'FORTRAN FUNCTION') then ! Make lower case (? pff, 2013)
                   do i=1,3
                      iii=ichar(cbc3(i:i))+32
                      if (iii.ge.97.and.iii.le.122) cbc3(i:i)=char(iii)
                   enddo
                   cbc(f,e,ifld)=cbc3
                   bclab(1)          = cbc3(1:1)
                   call inflow(e,f,ifld,cbc(f,e,ifld))
                   call rzero(bc(1,f,e,ifld),5)
                   ibc(f,e,ifld)=0
                  elseif (choice.eq.'CONSTANT') then
                     if (bcchoice.eq.'SHEAR    (LOCAL)') then
                       if (if3d) then
                         call prs(' Type #1-comp of Shear:$')
                         call rer(BC(2,f,e,ifld))
                         call prs(' Type #2-comp of Shear:$')
                         call rer(BC(3,f,e,ifld))
                       else
                         call prs(' Type Shear:$')
                         call rer(BC(2,f,e,ifld))
                       endif
                     elseif (bcchoice(10:16).eq.'(LOCAL)') then
                        call prs(' Type:$')
                        call prs(' NORMAL-'//bcchoice//'$')
                        call rer(BC(1,f,e,ifld))
                        if (if3d) then
                           call prs(' #1 TANGENTIAL '//bcchoice//'$')
                           call rer(BC(2,f,e,ifld))
                           call prs(' #2 TANGENTIAL '//bcchoice//'$')
                           call rer(BC(3,f,e,ifld))
                        else
                           call prs(' TANGENTIAL '//bcchoice//'$')
                           call rer(BC(2,f,e,ifld))
                        endif
                     elseif (bcchoice.eq.'TEMP') then
                        call prs('Type TEMP:$')
                        call rer(BC(1,f,e,ifld))
                     elseif (bcchoice.eq.'FLUX') then
                        call prs('Type FLUX:$')
                        call rer(BC(1,f,e,ifld))
                     elseif (bcchoice.eq.'OUTFLOW/N' .or.
     $                          bcchoice.eq.'OUTFLOW') then
                        if (.NOT.IFSTRS) then
                           call prs(' Exit pressure set to zero$')
                           BC(1,f,e,ifld) = 0.0
                        else
                           call prs(' Type Exit Pressure:$')
                           call rer(BC(1,f,e,ifld))
                        endif
                     elseif (bcchoice.eq.'FLUID LAYERS') then
                      write(s,'('' Type Surface Tension:'')')
                      call prs(s//'$')
                      call rer(BC(4,f,e,ifld))
                      call overlap(e,f,eo,fo)
                      cbc(  fo,eo,ifld) = 'MPI'
                      BC(4,fo,eo,ifld)=BC(4,f,e,ifld)
                     elseif (bcchoice.eq.'SHEAR') then
                       call prs(' Type X-comp of Shear:$')
                       call rer(BC(1,f,e,ifld))
                       call prs(' Type Y-comp of Shear:$')
                       call rer(BC(2,f,e,ifld))
                       if (if3d) then
                         call prs(' Type Z-comp of Shear:$')
                         call rer(BC(3,f,e,ifld))
                       endif
                     elseif (bcchoice.eq.'CONVECTION') then
                       call prs(' Heat Transfer Coefficient (h)>$')
                       call rer(BC(2,f,e,ifld)) ! 2nd storage Contains h
                       call prs(' Temperature at Infinity (Tinf)>$')
                       call rer(BC(1,f,e,ifld))
                     elseif (bcchoice.eq.'FREE SURFACE') then
                       call prs(' Pressure of Gas (Normal Traction)>$')
                       call rer(BC(1,f,e,ifld))
                       call prs(' Shear (Tangential Traction)>$')
                       call rer(BC(2,f,e,ifld))
                       if (if3d) then
                          call prs(' #2 Component of Shear '//
     $                       '(Tangential Traction)>$')
                          call rer(BC(3,f,e,ifld))
                       endif
                       call prs(' Surface tension Coefficient>$')
                       call rer(BC(4,f,e,ifld))
                     elseif (bcchoice.eq.'RADIATION') then
                       call prs(' Temperature at Infinity (Tinf)>$')
                       call rer(BC(1,f,e,ifld))
                       call prs(' Type Product of Emissivity$')
                       call prs(' and Boltzmanns Constant >$')
                       call rer(BC(2,f,e,ifld))
                     elseif (bcchoice.eq.'SOLID FRONT'
     $                  .or.    bcchoice.eq.'LIQUID FRONT') then
                       call prs(' Freezing Temperature>$')
                       call rer(BC(4,f,e,ifld))
                       call prs(' Rho * Latent Heat>$')
                       call rer(BC(5,f,e,ifld))
                     elseif (bcchoice.eq.'MELTING FRONT') then
c                      Fill up both sides of TEMPERATURE B.C.
                       call prs(' Freezing Temperature>$')
                       call rer(BC(4,f,e,2))
                       call prs(' Rho * Latent Heat>$')
                       call rer(bc(5,f,e,2))
                       call overlap(e,f,eo,fo) ! Export to overlapping edge
                       BC (4,fo,eo,2)=BC (4,f,e,2)
                       BC (5,fo,eo,2)=BC (5,f,e,2)
                       if (ee.gt.e) then
                          cbc(f ,e ,2)='MLI'
                          cbc(fo,eo,2)='MCI'
                       else
                          cbc(f ,e ,2)='MCI'
                          cbc(fo,eo,2)='MLI'
                       endif
                     else
c                       B.C. requiring 3 x,y,z components
                        write(s,'('' Type:'')')
                        call prs(s//'$')
                        write(s,'('' X-'',A26)')bcchoice
                        call prs(s//'$')
                        call rer(BC(1,f,e,ifld))
                        write(s,'('' Y-'',A26)')bcchoice
                        call prs(s//'$')
                        call rer(BC(2,f,e,ifld))
                        if (if3d) then
                           write(s,'('' Z-'',A26)')bcchoice
                           call prs(s//'$')
                           call rer(BC(3,f,e,ifld))
                        endif
                     endif
                  endif

               elseif (bcchoice.eq.'PERIODIC-AUTO') then ! Automatic selection of Periodic bc

                  call fndsida(jside,jel,f,e,ifld)
                  if (jside.eq.0) then
                   call prs(
     $            'Cannot find a corresponding side, choose a new BC.$')
                   goto 63
                  else
                    call letbc(jside,jel,ifld,bclab)
                    call letbc(f,e,ifld,bclab)

                    cbc( f,e,ifld) = 'P  '
                    bc(1,f,e,ifld) = jel
                    bc(2,f,e,ifld) = jside

                    cbc( jside,jel,ifld) = 'P  '
                    bc(1,jside,jel,ifld) = e
                    bc(2,jside,jel,ifld) = f
                    ibc(jside,jel,ifld)  = e
                    ibc(f,e,ifld)  = jel
                  endif

               elseif (bcchoice.eq.'PERIODIC') then ! Periodic B.C.
c                 !!?? In general, check if old b.c. was a periodic or
c                 !!?? B.c. before writing over it!
                  call prs('Choose corresponding periodic with mouse.$')
                  call prs('Use red menu area if side on other level.$')
124               write(s,'('' Periodic Side>'')')
                  call prs(s//'$')
116               call mouse(xmouse,ymouse,button) ! look for mouse input for periodic side
                  xp=xmouse
                  yp=ymouse
                  if (button.ne.'LEFT') then
                   call prs('Use left button to choose periodic side.$')
                   goto 124
                  endif
                 iother = 0
                 if (xscr(xmouse).gt.1.0) then
                  iother = 1
                  call prs('Type Level:$')
                  call rei(iplev)
                  if (iflearn) write(3,*) IPLEV
                  call prs('Type Apartment letter (with correct case)$')
                  call res(ipapt,1)
                  if (IFLEARN) write(3,*) IPAPT
                  call prs('Type Side number:$')
                  call rei(ipside)
                  if (IFLEARN)write(3,*) IPSIDE
                  fp = ipside ! Find Closest Side, Use data from other level
                  ep = 1
c                 do ee=1,nel
c                    if (NUMAPT(ee).eq.IPLEV .AND.LETAPT(ee).eq.IPAPT) ep=ee
c                 enddo
                  call prsi('Displaying BC label From level$',iplev)
                  call prs('On this level.  Do not be alarmed.$')
                 else !  Look on this level
                  rmin=1.e22
                  do 105 ee=1,nel
                  do 105 ff=1,nsides
                     if (numapt(ee).ne.ilevel) goto 105
                     r=((xp-sides(ee,ff,1))**2+(yp-sides(ee,ff,2))**2)
                     if (r.lt.rmin) then
                        rmin=r
                        fp = ff
                        ep = ee
                     endif
105               continue
                  if (fp.eq.5.or.fp.eq.6) then ! check IF it's above or below center
                     fp=5
                     if (yp.gt.sides(ep,fp,2)) fp=6
                  endif
                  call prsi('fp=$',fp)
               endif
               if (ep.eq.e .and. fp.eq.f) then
                  call prs('** ERROR **  A side cannot be periodic$')
                  call prs('with itself !  Choose the (different) $')
                  call prs('side that the flashing side is $')
                  call prs('periodic with.$')
                  goto 124
               endif
               write(s,110) ep,fp,e,f
               call prs(s//'$')
110            format('   ',2I7,'  IS PERIODIC WITH ',2I7)
c              Now Enter Information into Periodicity Array
c              !?? WHATS THE ORIENTATION??
               IORIEN=1
               call letbc(fp,ep,ifld,bclab)
               cbc( fp,ep,ifld) = 'P  '
               BC(1,f ,e ,ifld) = ep
               BC(2,f ,e ,ifld) = fp
c              BC(3,f ,e ,ifld) = IORIEN
               BC(1,fp,ep,ifld) = e
               BC(2,fp,ep,ifld) = f
c              BC(3,fp,ep,ifld) = IORIEN
               ibc(f ,e ,ifld)  = ep    ! high-precision periodic bc
               ibc(fp,ep,ifld)  = e
            endif
c           Post-pick warnings
            if (bcchoice.eq.'AXIS') then ! Axisymmetric B.C.
               if(abs(y(iside,e)/yfac).gt.0.01)
     $            call prs('** WARNING ** AXIS B.C. FOR Y=0 ONLY!$')
               if (f .NE. 1) then
                  call prs('*** Error!!! Gauss-Jacobi Mesh$')
                  call prs('allows only side 1 to be on axis.$')
                  call prsI('Please delete element$',e)
                  call prs('and Re-enter with axial side first.$')
                  call prs('Hit <cr> to continue$')
                  call res(line,0)
               endif
            endif
           endif

           call flash_element_b(ifade,iflash,e,f,ilevel,ifld,fcorns)

900   continue
c
c==================================================================
c     Entering MODIFY BC MENU area (apparently)
c==================================================================
c
313   continue
      item(1)='END  LEVEL'
      item(2)='MODIFY B.C.'
      item(3)='SHOW   B.C.'
      item(4)='INTERNAL B.C.'
      nchoic=3
      if (ifld.eq.1 .or. ifld.eq.2 .and. .not. ifflow) nchoic=4

      nchoic=nchoic+1
      item(nchoic)='ZOOM'
      call menu (xmouse,ymouse,button,'MODIFY B.C.')

      cmenu = '   '
      if (choice.eq.'INTERNAL B.C.') cmenu = 'INTERNAL'

      if (choice.eq.'MODIFY B.C.' .or. choice.eq.'SHOW   B.C.'
     $     .or.choice.eq.'INTERNAL B.C.') then
          call prs('Enter side with mouse$')
          call mouse(xp,yp,button) ! Find closest BOUNDARY side
          rmin=1.e20
          do ee=1,nel
          do ff=1,nsides
             if (ff.eq.5.and.yp.gt.sides(ee,ff,2)) goto 1113 ! above or below
             if (ff.eq.6.and.yp.lt.sides(ee,ff,2)) goto 1113 ! center- 5 or 6
             if (numapt(ee).eq.ilevel .and. (cbc(ff,ee,ifld).ne.'E  '
     $          .or.choice.eq.'INTERNAL B.C.')) then
                r = ((xp-sides(ee,ff,1))**2 + (yp-sides(ee,ff,2))**2)
                if (r.lt.rmin) then
                  rmin=r
                  isimod = ff
                  ielmod = ee
                endif
             endif
1113         continue
          enddo
          enddo

          if (choice.eq.'SHOW   B.C.') then ! Show   B.C.
             write(s,'(A7,i11,A6,I3,A1)')'Element',ielmod,
     $          '  Side',isimod,'$'
             call prs(s)
             write (s,'(a5,5g12.5,a1)') cbc(isimod,ielmod,ifld)
     $                            ,(bc(k,isimod,ielmod,ifld),k=1,5),'$'
             call prs(s)
             goto 313
          elseif (choice.eq.'MODIFY B.C.' .or.
     $            choice.eq.'INTERNAL B.C.') then
c            Reset to zero and send back into loop.  Loop will demand bc
             if (cbc(isimod,ielmod,ifld).eq.'P  ') then ! zero out the periodic side
                eper=bc(1,isimod,ielmod,ifld)
                eper=ibc(isimod,ielmod,ifld)
                fper=bc(2,isimod,ielmod,ifld)
                call letbc(fper,eper,ifld,'   ')
                cbc( fper,eper,ifld)='   '
                call rzero(bc(1,fper,eper,ifld),5)
                ibc (fper,eper,ifld)=0 ! high-precision periodic bc
             endif
             call letbc(isimod,ielmod,ifld,'   ')
             cbc( isimod,ielmod,ifld)='   '
             if (choice.ne.'INTERNAL B.C.') then
c               Retain old elemental connectivity for internal B.C.'s
                call rzero(bc(1,isimod,ielmod,ifld),5)
             endif
             goto 800
          endif
       elseif (choice.eq. 'END  LEVEL') then
          call check_end_level(icontinue,ilevel,ifld)
          if (icontinue.gt.0) goto 1010
       elseif (choice.eq. 'ZOOM') then
          call rezoom
          goto 313
       endif

1000  continue ! Just keep on going in loop (to 1000)
1999  continue
      ilevel=nlevel

      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod
      include 'basics.inc'
      character*3 cb_clear
 
      tol=1.e-1
      m  = 9
      if (if3d) m = 27 

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      goto 1  ! BELOW IS WORKING ONLY FOR Z AT THE MOMENT  (2/19/13,pff)

      if (if3d) then  ! Now, supported only for 3D:

       call prs('Is this for X-, Y-, or Z-extrema (x,y,z or other) ?$')
       call res(ans,1)
       call capit(ans,1)

       if (ans.eq.'X'.or.ans.eq.'Y'.or.ans.eq.'Z') then
         call mkside
         call gencen
         cb_clear = 'O  '
         if (ans.eq.'X') call autoperiodz(ans,tol,cb_clear,x27,m)
         if (ans.eq.'Y') call autoperiodz(ans,tol,cb_clear,y27,m)
         if (ans.eq.'Z') call autoperiodz(ans,tol,cb_clear,z27,m)
         return
       endif
      endif

    1 continue
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      call autoperiod1(nfail,tol)      ! Standard periodicity search
      if (nfail.eq.0) return

      write(6,*) 'NFAIL:',nfail,tol
      do k=1,10
         tol = .7*tol
         call autoperiod2(nfail,tol)
         write(6,*) 'NFAIL:',nfail,tol
         if (nfail.eq.0) return
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod1(nfail,tol)
      include 'basics.inc'
      character*3 cb,ck
c
c     Set any remaining candidate bcs to P
c
      if (.not.if_lattice) call get_lattice_per_bc

      call mkside   ! moved, 9/1/09 pff
      call gencen

      nsides = 2*ndim
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'   ') then
C           Automatic selection of Periodic bc
            call fndsidb(js,je,is,ie,if,tol)
            if (js.ne.0) then

               call letbc(js,je,if,bclab)
               call letbc(is,ie,if,bclab)

               call rzero(bc(1,is,ie,if),5)
               call rzero(bc(1,js,je,if),5)

               cbc (is,ie,if) = 'P  '
               ibc (is,ie,if) = je       ! high-precision periodic bc
               bc(1,is,ie,if) = je
               bc(2,is,ie,if) = js

               cbc (js,je,if) = 'P  '
               ibc (js,je,if) = ie       ! high-precision periodic bc
               bc(1,js,je,if) = ie
               bc(2,js,je,if) = is


               if (nel.le.1000.or.mod(ie,100).eq.0)
     $         write(6,1) if,ie,is,je,js
    1          format('autoper:',5i8)
            endif
         endif
      enddo
      enddo
      enddo
c
c     consistency check
c
      nfail = 0
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'P  ') then
C           check for consistency
            je = bc(1,is,ie,if)
            js = bc(2,is,ie,if)
            ke = bc(1,js,je,if)
            ks = bc(2,js,je,if)

            je = ibc (is,ie,if)     ! high-precision periodic bc
            ke = ibc (js,je,if)     ! high-precision periodic bc

            ck = cbc (js,je,if)

            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) then
               write(6,*) ie,is,je,js,ke,ks,ck,' fail',nfail
               cbc(is,ie,if)='   '
               nfail = nfail + 1
               write(6,11) (sides(ie,is,k),k=1,3),ie,is,' side i'
               write(6,11) (sides(je,js,k),k=1,3),je,js,' side j'
               write(6,11) (sides(ke,ks,k),k=1,3),ke,ks,' side k'
               sides(1,1,1) = -nel
               sides(1,1,1) = sqrt(sides(1,1,1))
            endif
  11        format(1p3e14.5,i8,i3,a7)
         endif
      enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine autoperiod2(nfail,tol)
      include 'basics.inc'
      character*3 cb,ck,bclab
c
c     Set any remaining candidate bcs to P
c
      if (.not.if_lattice) call get_lattice_per_bc

      call mkside   ! moved, 9/1/09 pff
      call gencen

      nsides = 2*ndim
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'   '.or.cb.eq.'P  ') then
C           Automatic selection of Periodic bc
            call fndsidb(js,je,is,ie,if,tol)
            if (js.ne.0) then
c              call letbc(js,je,if,bclab)
c              call letbc(is,ie,if,bclab)
               cbc( IS,IE,IF) = 'P  '
               call rzero(bc(1,is,ie,if),5)
               call rzero(bc(1,js,je,if),5)
               BC(1,IS,IE,IF) = JE
               BC(2,IS,IE,IF) = JS
               cbc( JS,JE,IF) = 'P  '
               BC(1,JS,JE,IF) = IE
               BC(2,JS,JE,IF) = IS
               ibc (is,ie,if) = je  ! high-precision periodic bc
               ibc (js,je,if) = ie  ! high-precision periodic bc

               if (nel.le.1000.or.mod(ie,100).eq.0)
     $         write(6,1) if,ie,is,je,js
    1          format('autoper:',5i8)
            endif
         endif
      enddo
      enddo
      enddo
c
c     consistency check
c
      nfail = 0
      do if=1,nflds
      do ie=1,nel
      do is=1,nsides
         cb = cbc(is,ie,if)
         if (cb.eq.'P  ') then
C           check for consistency
            je = bc(1,is,ie,if)
            js = bc(2,is,ie,if)
            ke = bc(1,js,je,if)
            ks = bc(2,js,je,if)
            ck = cbc(js,je,if)

            je = ibc (is,ie,if)     ! high-precision periodic bc
            ke = ibc (js,je,if)     ! high-precision periodic bc

            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) cbc(is,ie,if)='   '
            if (ck.ne.'P  '.or.ke.ne.ie.or.ks.ne.is) nfail = nfail + 1
         endif
      enddo
      enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_vec_side(vec3,iside,ie)
      include 'basics.inc'
      real vec1(3),vec2(3),vec3(3)

      ic1 = icrn(1,iside)
      ic2 = icrn(2,iside)
      ic3 = icrn(3,iside)
      ic4 = icrn(4,iside)

c     write(s,9) ie,iside,ic1,ic2,ic3,ic4
c   9 format(' iesc14:',6i5,'$')
c     if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)

      IF (IF3D) THEN
         vec1(1)=x(ic3,ie)-x(ic1,ie)
         vec1(2)=y(ic3,ie)-y(ic1,ie)
         vec1(3)=z(ic3,ie)-z(ic1,ie)
         vec2(1)=x(ic4,ie)-x(ic2,ie)
         vec2(2)=y(ic4,ie)-y(ic2,ie)
         vec2(3)=z(ic4,ie)-z(ic2,ie)
         call cross(vec3,vec1,vec2)
         call norm3d(vec3)
      else
         vec1(1)=x(ic2,ie)-x(ic1,ie)
         vec1(2)=y(ic2,ie)-y(ic1,ie)
         vec3(1)= vec1(2)
         vec3(2)=-vec1(1)
         vec3(3)=0.0
         call norm3d(vec3)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine fndsidb(kside,ke,iside,ie,if,tol)
      include 'basics.inc'
      real vec1(3),vec2(3),vec3(3),vecj(3)
C
C     Find a side JSIDE,JEL which corresponds to ISIDE,IE and doesn't
C     have BC's set already
C
      NSIDES=NDIM*2
C
C     Define Normal to this plane (side), and find the element
C     center which also lies on this line, and is farthest from
C     the current element.
c
c     CALL MKSIDE   ! moved, 9/1/09 pff
c     CALL GENCEN
C
C     Find Sides' Normal Vector at midpoint

      call get_vec_side(vec3,iside,ie)

      ic1 = icrn(1,iside)
      ic2 = icrn(2,iside)
      ic3 = icrn(3,iside)
      ic4 = icrn(4,iside)

c     write(s,9) ie,iside,ic1,ic2,ic3,ic4
c   9 format(' iesc14:',6i5,'$')
c     if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)

      IF (IF3D) THEN
         vec1(1)=x(ic3,ie)-x(ic1,ie)
         vec1(2)=y(ic3,ie)-y(ic1,ie)
         vec1(3)=z(ic3,ie)-z(ic1,ie)
         vec2(1)=x(ic4,ie)-x(ic2,ie)
         vec2(2)=y(ic4,ie)-y(ic2,ie)
         vec2(3)=z(ic4,ie)-z(ic2,ie)
         call cross(vec3,vec1,vec2)
         call norm3d(vec3)
      else
         vec1(1)=x(ic2,ie)-x(ic1,ie)
         vec1(2)=y(ic2,ie)-y(ic1,ie)
         vec3(1)= vec1(2)
         vec3(2)=-vec1(1)
         vec3(3)=0.0
         call norm3d(vec3)
      ENDIF
c
      if (if_lattice) then          ! Make vector point in lattice direction
         d1 = dotprod(vec3,ulat1)
         d2 = dotprod(vec3,ulat2)
         d3 = dotprod(vec3,ulat3)
         if (d1.le.d2.and.d1.le.d3) then
            call copy(vec3,ulat1,3)
         elseif (d2.le.d3) then
            call copy(vec3,ulat2,3)
         else
            call copy(vec3,ulat3,3)
         endif
      endif
c
c     write(s,19) vec3(1),vec3(2),vec3(3)
c  19 format(' vec3:',3e13.5,'$')
c     if (nel.lt.500.or.mod(ie,100).eq.0) call prs(s)
c
c     EPS2=(EPSM*RCEN(IE))**2
      EPSM=tol
      EPS2=EPSM**2
      D2MN=10*EPS2
c
      KE=0
      KSIDE=0
      DSTMAX=0.0
c
      do 100 je=1,nel
      do 100 jside=1,nsides
C
C        Don't find yourself.
         if (je.eq.ie.and.jside.eq.iside) goto 100
C
C        Don't find a bc which is already defined. (to be 'E')
         if (cbc(jside,je,if).eq.'E  ') goto 100
c
c        Check if dot product of surface normals are O(1)
c
         call get_vec_side(vecj,jside,je)

C        OK, is the center of this element side close to the line?

         vec1(1)=sides(ie,iside,1)-sides(je,jside,1)
         vec1(2)=sides(ie,iside,2)-sides(je,jside,2)
         vec1(3)=sides(ie,iside,3)-sides(je,jside,3)
         distp=dotprod(vec1,vec1)
         dis13=(dotprod(vec1,vec3))**2
         dist2=distp-dis13
         if (distp.gt.0) DIST2=DIST2/DISTP
c        write(6,11) dist2,distp,eps2,dstmax,je,jside
c  11    format(1p4e12.4,2i6,' dist2')
c        if (dist2.le.eps2.and.distp.gt.dstmax) then
         if (dist2.le.d2mn) then
            ke    = je
            kside = jside
            dstmax= distp
c     write(6,12) distp,dist2,(sides(ie,iside,k),k=1,3),ie,iside,' si'
c     write(6,12) dis13,d2mn ,(sides(je,jside,k),k=1,3),je,jside,' sj'
c     write(6,12) dis13,d2mn ,(vec1(k),k=1,3),je,jside,'vc1'
c     write(6,12) dis13,dotnorm ,(vecj(k),k=1,3),je,jside,'vcj'
c  12 format(1p5e14.5,i8,i3,1x,a3)

            d2mn  = dist2

         endif
  100 continue

      return
      end
c-----------------------------------------------------------------------
      function dif_rel_inf_norm(x,y,n)
      real x(n),y(n)

      dif_rel_inf_norm = 0.

      dmax = 0.
      xmax = 0.
      ymax = 0.

      do i=1,n
         dif = abs(x(i)-y(i))
         dmax = max(dmax,dif)
         xmax = max(xmax,abs(x(i)))
         ymax = max(ymax,abs(y(i)))
      enddo

      den = max(xmax,ymax)
      if (den.gt.0) dif_rel_inf_norm = dmax / den

      return
      end
c-----------------------------------------------------------------------
