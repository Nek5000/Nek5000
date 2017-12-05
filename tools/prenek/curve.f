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
      subroutine curves
#     include "basics.inc"
      integer icalld,e
      save    icalld
      data    icalld /0/
      logical iftmp
C
      iftmp=ifgrid
      ifgrid=.false.

1     continue
      ITEM(1)='BUILD MENU'
      ITEM(2)='STRAIGHTEN CURVE'
      ITEM(3)='MAKE CIRCLE'
      ITEM(4)='Convert to Midside Nodes'
      ITEM(5)='Autosphere'
c     ITEM(5)='Tile with hexagons'
c     ITEM(6)='TRANSITION HEXAGONS'
      ITEM(6)='Refine hexagons'
c     ITEM(7)='MAKE SPLINE'
      ITEM(7)='Arc-Circle Transform'
c     ITEM(5)='MAKE SINE WAVE'
C     ITEM(6)='FORTRAN FUNCTION'
      ITEM(8)='Convert Midside to Circle'
      NCHOIC=8
      nchoic=nchoic+1
      ITEM(nchoic)='Circular-Cartesian'
      if (if3d) then
         nchoic=nchoic+1
         ITEM(nchoic)='SPHERICAL MESH'
         nchoic=nchoic+1
         ITEM(nchoic)='STRETCH SPHERE'
         nchoic=nchoic+1
c        ITEM(nchoic)='GEN TEST MESH'
         ITEM(nchoic)='NEAR-WALL SPHERE'
c     else
c        nchoic=nchoic+1
c        item(nchoic)='CIRC MESH'
      endif
      nchoic=nchoic+1
      ITEM(nchoic)='RENUMBER ELEMENTS'
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'CURVE SIDES')
3013  CONTINUE
      IF(CHOICE.EQ.'MAKE SPLINE')THEN
          CALL GETEDG(ISID,IELS)
C         GET 2 Auxiliary control points
          CALL COLOR(1)
          call XDOT(x(isid,iels),y(isid,iels))
          CALL PRS(
     $    'Enter First Control point (adjacent to hilighted corner):$')
          CALL MOUSE(XMOUS1,YMOUS1,BUTTON)
          CALL COLOR(8)
          CALL xdot(XMOUS1,YMOUS1)
          CALL PRS(
     $    'Enter Second Control point (adjacent to other corner):$')
          CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
          CALL COLOR(8)
          CALL xdot(XMOUSE,YMOUSE)
          CALL COLOR(1)
          CALL DRAWEL(-IELS)
          CCURVE (ISID,IELS)='S'
          CURVE(1,ISID,IELS)=XMOUS1
          CURVE(2,ISID,IELS)=YMOUS1
          CURVE(3,ISID,IELS)=XMOUSE
          CURVE(4,ISID,IELS)=YMOUSE
C         Export copies curve ICURVE to all overlapping sides
          CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'RENUMBER ELEMENTS') THEN
           call renum_special
           GOTO 1
      ELSE IF(CHOICE.EQ.'GEN TEST MESH') THEN
           CALL GENMESH
           GOTO 1
      ELSE IF(CHOICE.EQ.'Clean up vertices') THEN
           CALL VERTADJ
           GOTO 1
      ELSE IF(CHOICE.EQ.'Refine hexagons') THEN
           call hexagon_refine
           GOTO 1
      ELSE IF(CHOICE.EQ.'TRANSITION HEXAGONS') THEN
           call hex_transition_3d
           GOTO 1 
      ELSE IF(CHOICE.EQ.'Tile with hexagons') THEN
           call hexagon_tile
           GOTO 1
      ELSE IF(CHOICE.EQ.'STRETCH SPHERE') THEN
           call stretch_sphere
           GOTO 1
c     ELSE IF(CHOICE.EQ.'CIRC MESH') THEN
c          call circ_mesh 
c          GOTO 1
      ITEM(7)='Arc-Circle Transform'
      ELSEIF (CHOICE.EQ.'Arc-Circle Transform') THEN
           call make_arc_circ
           goto 1
      ELSE IF(CHOICE.EQ.'SPHERICAL MESH') THEN
           CALL SPHMESH
           GOTO 1
      ELSE IF(CHOICE.EQ.'NEAR-WALL SPHERE') THEN
           call nw_sphmesh
           goto 1
      ELSE IF(CHOICE.EQ.'STRAIGHTEN CURVE')THEN
          call getedg(isid,iels)
          call drawel(-iels)
          call rzero(curve(1,isid,iels),6)
          ccurve (isid,iels)=' '
          write(s,'(1x,a14,i3,i4)')' Straightening',isid,iels
          call prs(s//'$')
          call export(isid,iels) ! Export to overlapping sides
      ELSE IF(CHOICE.EQ.'Convert to Midside Nodes') THEN
          call gencen
          do e=1,nel
             call fix_m_curve(e)
          enddo
          ifmid  = .true.
          ifcstd = .true.
          goto 1

      elseif (choice.eq.'Convert Midside to Circle') then
          call convert_m_to_c_all ! This works only for 2D at present

      elseif (choice.eq.'Autosphere') then

          if (if3d) then
             call filled_sphere_oct ! Builds sphere or hemisphere of diameter 1
          else
             call filled_cyl_quadrant ! Builds filled cylinder quadrant
          endif

          call redraw_mesh
          goto 1

      elseif (choice.eq.'Circular-Cartesian') then
             call cirmesh
          goto 1

      ELSE IF(CHOICE.EQ.'MAKE CIRCLE')THEN
            CALL GETEDG(ISID,IELS)
            xcent=0.25*vlsum(x(1,iels),4)
            ycent=0.25*vlsum(y(1,iels),4)
            CALL GWRITE(XCENT,YCENT,1.0,'*$')
            CALL PRS(' Enter Circle Radius for side of element *$')
            CALL PRS(' >0 for convex element, <0 for concave$')
            call rer(radius)
c           CALL KEYPADm(RADIUS,button,xmouse,ymouse)
c           if (button.eq.'RIGHT') then
c              clicking near an element 
c              call getside(jel,jside,xmouse,ymouse)
c              RADIUS =  -CURVE(1,jside,jel)
c           endif
            CALL DRAWEL(-IELS)
            CURVE(1,ISID,IELS)=RADIUS
            CCURVE(ISID,IELS)='C'
C
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'FORTRAN FUNCTION')THEN
            CALL GETEDG(ISID,IELS)
            CALL DRAWEL(-IELS)
            CCURVE(ISID,IELS)='F'
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'MAKE SINE WAVE')THEN
            CALL GETEDG(ISID,IELS)
            CALL PRS(' **WARNING** sine wave cannot be rotated.$')
            CALL PRS(' **WARNING** sine wave cannot be refined.$')
            CALL PRS('It must be of the form:   Y=Yo + A*sin(alpha x)$')
            CALL PRS
     $      ('First Enter the amplitude  A: (Negative to abort)$')
            call rer(CURVE(1,ISID,IELS))
            IF(CURVE(1,ISID,IELS).LE.0.0)GO TO 1
            CALL PRS('Enter the wave number  alpha:$')
            call rer(CURVE(2,ISID,IELS))
            CALL PRS('Enter the Y-offset Yo:$')
            call rer(CURVE(3,ISID,IELS))
            CALL DRAWEL(-IELS)
            CCURVE(ISID,IELS)='W'
C
            CALL EXPORT(ISID,IELS)
      ELSE IF(CHOICE.EQ.'BUILD MENU')THEN
            NCURVE=0
            DO 114 IE=1,NEL
            DO 114 IEDGE=1,12
               IF (CCURVE(IEDGE,IE).NE.' ') THEN
                  NCURVE=NCURVE+1
               ENDIF
114         CONTINUE
            IFGRID=IFTMP
            RETURN
      ENDIF
C     After modifiacation of floor, check if there are any elements in floors
C     above.  If not, modify ceiling also.
      IF(.NOT.IFCEIL .AND.IF3D)THEN
         DO 11 I=1,IEL
            IF(NUMAPT(I).GT.ILEVEL) THEN
               CALL PRS(' Presence of elements on higher '//
     $         'floors inhibits modifying ceiling also$')
               RETURN
            ENDIF
11       CONTINUE
C        Since we got here, there must be no elements on upper floors
C        Copy modifications.
         DO 12 I=1,6
             CURVE(I,ISID+4,IELS)=CURVE(I,ISID,IELS)
12        CONTINUE
          CCURVE(ISID+4,IELS)=CCURVE(ISID,IELS)
C         Export copies curve ICURVE to all overlapping sides
          CALL EXPORT(ISID+4,IELS)
      ENDIF
      GO TO 1
      END
      SUBROUTINE MESGEN
#     include "basics.inc"
      DIMENSION DXCHEB(20),DYCHEB(20),XCHEB(20),YCHEB(20)
      REAL CSPACE(100),XCRVED(100),YCRVED(100),XFRAC(20)
      REAL XPTSEL(NXM,NYM),YPTSEL(NXM,NYM)
C     Make mesh and draw on screen
      CALL PRSIS(' Generating Mesh Points for $',NEL,' Elements$')
      IF(NEL .GT. NELM-2 ) CALL PRS(' ***** ERROR ******* $')
      IF(NEL .GT. NELM-2 ) CALL PRSIS('EXCEEDED MAXIMUM OF$',NELM,
     $' ELEMENTS$')
      call color(8)
      IF (NX.GT.NXM) THEN
         WRITE(S,104) NX,NXM
  104    FORMAT(' Warning, current N exceeds NXM, resetting from',I3
     $         ,' to',I3,'.$')
         CALL PRS(S)
         NX=NXM
         NY=NYM
         NZ=NZM
      ENDIF
      DO 105 I=1,NX
         XFRAC(I)=XXPTS(I)/2.0 + 0.5
105   CONTINUE
      DO 108 IEL = 1,NEL
         DO 106 I = 1,NX
           DO 106 J = 1,NY
C           Calculate original xpts based on straight sides
            xpts(i,j,iel) =  .25* (
     $      x(1,iel)*(xxpts(i)-1.0)*(yypts(j)-1.0)-
     $      x(2,iel)*(xxpts(i)+1.0)*(yypts(j)-1.0)+
     $      x(3,iel)*(xxpts(i)+1.0)*(yypts(j)+1.0)-
     $      x(4,iel)*(xxpts(i)-1.0)*(yypts(j)+1.0))
            ypts(i,j,iel) =  .25* (
     $      y(1,iel)*(xxpts(i)-1.0)*(yypts(j)-1.0)-
     $      y(2,iel)*(xxpts(i)+1.0)*(yypts(j)-1.0)+
     $      y(3,iel)*(xxpts(i)+1.0)*(yypts(j)+1.0)-
     $      y(4,iel)*(xxpts(i)-1.0)*(yypts(j)+1.0))
C           Initial values for temporary array
            xptsel(i,j) = xpts(i,j,iel)
            yptsel(i,j) = ypts(i,j,iel)
106       CONTINUE
           ICRV(IEL)=0
C          For 2-d only
           DO 165 IE=1,4
             IF(CCURVE(IE,IEL).NE.' ')ICRV(IEL)=ICRV(IEL)+IE**2
165        CONTINUE
           IF(ICRV(IEL).NE.0)THEN
C             Shift internal points to account for curviness
              DO 107 ISID=1,4
                IF(CCURVE(ISID,IEL).NE.' ')THEN
                  NPOINT=NX
                  CALL GETPTS(NPOINT,XFRAC,IEL,ISID,XCRVED,YCRVED)
                  DO 117 I = 1,NX
                    IF(ISID.EQ.1)THEN
                      II=I
                      JJ=1
                    ELSE IF(ISID.EQ.2)THEN
                      II=NX
                      JJ=I
                    ELSE IF(ISID.EQ.3)THEN
                      II=NX-I+1
                      JJ=NX
                    ELSE IF(ISID.EQ.4)THEN
                      II=1
                      JJ=NX-I+1
                     ENDIF
                    DXCHEB(I)=XCRVED(I)-XPTS(II,JJ,IEL)
                    DYCHEB(I)=YCRVED(I)-YPTS(II,JJ,IEL)
117               CONTINUE
                  DO 120 I=1,NX
                    DO 120 J=1,NY
                      IF(ISID.EQ.2)III=J
                      IF(ISID.EQ.4)III=NX-J+1
                      IF(ISID.EQ.1)III=I
                      IF(ISID.EQ.3)III=NX-I+1
                      IF(ISID.EQ.1) FAC = XFRAC(NX-J+1)
                      IF(ISID.EQ.2) FAC = XFRAC(I)
                      IF(ISID.EQ.3) FAC = XFRAC(J)
                      IF(ISID.EQ.4) FAC = XFRAC(NX-I+1)
C
                      XPTSEL(I,J)=XPTSEL(I,J)+DXCHEB(III)*FAC
                      YPTSEL(I,J)=YPTSEL(I,J)+DYCHEB(III)*FAC
120               CONTINUE
                ENDIF
107           CONTINUE
C             Now that the points in the temporary array are moved, copy to xpts
              DO 109 I=1,NX
                DO 109 J=1,NY
                  XPTS(I,J,IEL)=XPTSEL(I,J)
                  YPTS(I,J,IEL)=YPTSEL(I,J)
109           CONTINUE
            endif
108       CONTINUE
C
        DXMIN=1.0E10
        DYMIN=1.0E10
        DO 119 IEL=1,NEL
              DX=SQRT((XPTS(2,1,IEL)-XPTS(1,1,IEL))**2
     $               +(YPTS(2,1,IEL)-YPTS(1,1,IEL))**2)
              DY=SQRT((XPTS(1,2,IEL)-XPTS(1,1,IEL))**2
     $               +(YPTS(1,2,IEL)-YPTS(1,1,IEL))**2)
              DXMIN=AMIN1(DXMIN,DX)
              DYMIN=AMIN1(DYMIN,DY)
119      CONTINUE
C
        CALL HARD
        IFHARD=.TRUE.
        call color(8)
        if (nel.lt.1000) then
         DO 112 IEL = 1,NEL
C         Draw only those on this floor
          IF(NUMAPT(IEL).EQ.ILEVEL)THEN
            IF(ICRV(IEL).EQ.0)THEN
               DO 209 I = 1,NX
                 CALL MOVEC(XPTS(I,1 ,IEL),YPTS(I,1 ,IEL))
                 CALL DRAWC(XPTS(I,NY,IEL),YPTS(I,NY,IEL))
209            CONTINUE
               DO 211 J = 1,NY
                 CALL MOVEC(XPTS(1 ,J,IEL),YPTS(1 ,J,IEL))
                 CALL DRAWC(XPTS(NX,J,IEL),YPTS(NX,J,IEL))
211            CONTINUE
            ELSE
C              Draw it the hard, slow, curvy way
               DO 1109 J=1,NY-1
                   DO 1109 I = 1,NX
                     IF(CCURVE(2,IEL).NE.' '.AND. I.EQ.NX) GO TO 1108
                     IF(CCURVE(4,IEL).NE.' '.AND. I.EQ. 1) GO TO 1108
                     CALL MOVEC(XPTS(I,J  ,IEL),YPTS(I,J  ,IEL))
                     CALL DRAWC(XPTS(I,J+1,IEL),YPTS(I,J+1,IEL))
1108                CONTINUE
1109            CONTINUE
               DO 111 I=1,NX-1
                   DO 111 J = 1,NY
                     IF(CCURVE(1,IEL).NE.' '.AND. J.EQ. 1) GO TO 110
                     IF(CCURVE(3,IEL).NE.' '.AND. J.EQ.NY) GO TO 110
                     CALL MOVEC(XPTS(I  ,J,IEL),YPTS(I  ,J,IEL))
                     CALL DRAWC(XPTS(I+1,J,IEL),YPTS(I+1,J,IEL))
110                  CONTINUE
111            CONTINUE
C              Draw curved sides for hardcopy
               DO 130 IEDGE=1,4
                  IF(CCURVE(IEDGE,IEL).NE.' ')THEN
                     call movec(x(iedge,iel),y(iedge,iel))
                     CALL DRAWED(IEL,IEDGE,1)
                  ENDIF
130            CONTINUE
            ENDIF
          ENDIF
112       continue
         ENDIF
         CALL NOHARD
        IFHARD=.FALSE.
         NCURVE=0
         DO 114 IEL=1,NEL
            DO 114 IEDGE=1,12
               IF(CCURVE(IEDGE,IEL).NE.' ') NCURVE=NCURVE+1
114      CONTINUE
         RETURN
      END
c-----------------------------------------------------------------------
      subroutine getedg(isid,iels)
#     include "basics.inc"

      integer eindx(12),v  ! index of 12 edges into 3x3x3 tensor
      save    eindx        ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation
      integer count


      count=0
1     count=count+1
      if (count.gt.5) then
         call prs('Returning.$')
         return
      endif

      call prs(' Enter element side.$')
      call mouse(xmouse,ymouse,button)
      xh=xmouse
      yh=ymouse
      if (xscr(xmouse).gt.1.0) then
         call prs
     $   (' **ERROR**: Expecting you to define side.  Define side$')
         call prs(' by hitting mouse button in mesh area.  Try again.$')
         goto 1
      endif

      zh=0.0
      if (ilevel.gt.1) then
         do 10 i=1,ilevel-1
            zh=zh+height(i)
10       continue
      endif
      if (ifceil) zh=zh+height(ilevel)

      ii1=1
      ii2=4
      if (ifceil) ii1=5
      if (ifceil) ii2=8

      rmin=1.0e10
      do 50 iel=1,nel
      do 40 iedge=ii1,ii2

c        v=eindx(iedge)    ! This is problematic if curve is messed up!
c        dxh=xh-x27(v,iel) ! So, use mean vertex values, as defined 
c        dyh=yh-y27(v,iel) ! below.
c        dzh=zh-z27(v,iel)

         iv=mod1(iedge+1,4)
         xmid=(x(iedge,iel)+x(iv,iel))/2
         ymid=(y(iedge,iel)+y(iv,iel))/2
         zmid=(z(iedge,iel)+z(iv,iel))/2
         dxh=xh-xmid
         dyh=yh-ymid
         dzh=zh-zmid

         r=sqrt(dxh**2+dyh**2+dzh**2)

         if(r.lt.rmin) then
            rmin=r
            iels=iel
            isid=iedge
         endif

   40 continue
   50 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine export(iedge,iel) ! Problematic if edges messed up
#     include "basics.inc"
      integer iover(2,20),v1,v2,v3,v4,vi,vj
      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

C     Takes characteristics of curve ICURVE and copies it to any overlapping
C     edges
C
C     Use apartment letters to find if there are any downstairs ceilings
C     What was wrong with going up to roof level?
C     First Check for overlapping side
C
C     Find all elements and sides that get changed
C     Put original side in IOVER 1
      NOVER = 1
      IOVER(1,NOVER)=IEDGE
      IOVER(2,NOVER)=IEL

      v1=eindx(1)
      v2=eindx(3)
      v3=eindx(2)
      v4=eindx(4)
      vi=eindx(iedge)

      eps1=(x27(v2,iel)-x27(v1,iel))**2
     $    +(y27(v2,iel)-y27(v1,iel))**2
     $    +(z27(v2,iel)-z27(v1,iel))**2

      eps2=(x27(v4,iel)-x27(v3,iel))**2
     $    +(y27(v4,iel)-y27(v3,iel))**2
     $    +(z27(v4,iel)-z27(v3,iel))**2
      epsi = min(eps1,eps2)

      do 10 jel=1,nel

         eps1=(x27(v2,jel)-x27(v1,jel))**2
     $       +(y27(v2,jel)-y27(v1,jel))**2
     $       +(z27(v2,jel)-z27(v1,jel))**2

         eps2=(x27(v4,jel)-x27(v3,jel))**2
     $       +(y27(v4,jel)-y27(v3,jel))**2
     $       +(z27(v4,jel)-z27(v3,jel))**2
         epsj = min(eps1,eps2)
         eps  = 0.02*min(epsi,epsj)

         do 10 iedg=1,nedges
            vj=eindx(iedg)
            r = (x27(vi,iel)-x27(vj,jel))**2
     $        + (y27(vi,iel)-y27(vj,jel))**2
     $        + (z27(vi,iel)-z27(vj,jel))**2
            if(r.lt.eps) then
               if(iedg.ne.iedge .or. iel.ne.jel)then
                  nover = nover+1
                  iover(1,nover)=iedg
                  iover(2,nover)=jel
c     write(s,11) nover,iedg,jel,iedge,iel,r,eps,epsi,epsj
c     call prs(s)
   11 format('jiedgel,r,epsij:',5i3,4e10.3,'$')
               endif
            endif
10    continue
      IF(NOVER.GE.2)THEN
C        Do this only for overlapping elements on this level, not original
         DO 20 I=2,NOVER
            IELS=IOVER(2,I)
            IF(NUMAPT(IELS).EQ.NUMAPT(IEL)) CALL DRAWEL(-IELS)
20       CONTINUE
         DO 30 I=2,NOVER
            IEDG=IOVER(1,I)
            IELS=IOVER(2,I)
            DO 25 II=1,6
               CURVE(II,IEDG,IELS)=CURVE(II,IEDGE,IEL)
25          CONTINUE
            CCURVE (IEDG,IELS)=CCURVE( IEDGE,IEL)
C           Convex/concave for circle (except on different levels)
C            IF(NUMAPT(IEL).EQ.NUMAPT(IELS))THEN
            IF(LETAPT(IEL).NE.LETAPT(IELS))THEN
               IF(CCURVE(IEDG,IELS).EQ.'C')
     $         CURVE(1,IEDG,IELS)=-CURVE(1,IEDGE,IEL)
               IF(CCURVE(IEDG,IELS).EQ.'S')THEN
C              Spline on same level switch auxiliary points
                  CURVE(1,IEDG,IELS)=CURVE(3,IEDGE,IEL)
                  CURVE(2,IEDG,IELS)=CURVE(4,IEDGE,IEL)
                  CURVE(3,IEDG,IELS)=CURVE(1,IEDGE,IEL)
                  CURVE(4,IEDG,IELS)=CURVE(2,IEDGE,IEL)
               ENDIF
            ENDIF
30       CONTINUE
      ENDIF
      write(s,11) nover,iedg,iels,iedge,iel,r,eps,epsi,epsj
      call prs(s)
      DO 40 I=1,NOVER
         IELS=IOVER(2,I)
         IF(NUMAPT(IELS).EQ.NUMAPT(IEL)) CALL DRAWEL(IELS)
40    CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine stretch_sphere
c
c     Stretch all points in spherical shell
c
#     include "basics.inc"
C
      call prs(' Input sphere center (x,y,z):$')
      call rerrr(x0,y0,z0)
c
      call prs(' Input inner cut-off radius, r0:$')
      call rer(r0)
c
      call prs(' Input outer cut-off radius, r1:$')
      call rer(r1)
c
      call prs(' Input expansion factor ( =< 0 to abort):$')
      call rer(sfact)
      if (sfact.le.0) return
C
C     Stretch for hemisphere problem
C
C     1)  For all elements with z < z0, just stretch in X-Y
C
C     2)  For all elements with z > z0, stretch in X-Y and Z
C
C
C     Take care of pts first
C
      nvts=8
c     sfact = 1.25
c     z0 = 0.0346
c     rr = 0.5*0.5 + .0001
c     re = 0.5 + .0001
c
      z00= -9.e9
      z0 = 0.
      rr = 0.
      re = 0.
C
      r02 = r0*r0
      r12 = r1*r1
      do ie=1,nel
      do i=1,nvts
         xs = x(i,ie)-x0
         ys = y(i,ie)-y0
         zs = z(i,ie)-z0
         r2 = xs*xs + ys*ys + zs*zs
         if (r02.lt.r2 .and. r2 .lt. r12) then
            xs = sfact*xs
            ys = sfact*ys
            zs = sfact*zs
c
            x(i,ie) = x0 + xs
            y(i,ie) = y0 + ys
            z(i,ie) = z0 + zs
         endif
      enddo
      enddo
C
C     Take care of curved sides
C
      eps = .01*r1
      do ie = 1,nel
      do is = 1,12
c
         if (ccurve(is,ie).eq.'s') then
c           spherical side, x0  = curve(1,is,ie)
c           spherical side, y0  = curve(2,is,ie)
c           spherical side, z0  = curve(3,is,ie)
c           spherical side, rad = curve(4,is,ie)
c
            xs  = curve(1,is,ie)-x0
            ys  = curve(2,is,ie)-y0
            zs  = curve(3,is,ie)-z0
            re2 = xs*xs + ys*ys + zs*zs
            rad = curve(4,is,ie)
            if (re2.lt.eps .and. r0 .lt.rad .and. rad.le. r1) 
     $         curve(4,is,ie) = sfact*curve(4,is,ie)
         endif
c
c        if (ccurve(is,ie).eq.'C') then
c           cylindrical side, rad = curve(1,is,ie)
c           if (abs(curve(1,is,ie)).gt.re)
c    $         curve(1,is,ie) = sfact*curve(1,is,ie)
c        endif
c
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine find_long_edges ! Find long edge of a 2D element
c
#     include "basics.inc"
      integer e,en,et
c
      do e=1,nel
         long_edge(e) = long_edge_e(e)
      enddo
      return
      end
c-----------------------------------------------------------------------
      function long_edge_e(e) ! Find long edge of a 2D element
c
#     include "basics.inc"
      integer e,en,et
c
c
      dx = x(1,e)-x(4,e)
      dy = y(1,e)-y(4,e)
      dm = dx*dx + dy*dy
      long_edge_e = 4
c
      do i=1,3
         dx = x(i+1,e)-x(i,e)
         dy = y(i+1,e)-y(i,e)
         d2 = dx*dx + dy*dy
         if (d2.gt.dm) then
            long_edge_e = i
            dm = d2
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_tile
c
c     Tile a 2D region w/ hexagons (2D only, pff 10/10/05)
c
#     include "basics.inc"
      real hexscale(8,2),hexbase(8,2)
      save hexscale,hexbase
c
      integer e
      integer icalld
      save    icalld
      data    icalld  /0/
c
      if (icalld.eq.0) then
         icalld = icalld + 1
         one = 1.
         pi  = 4.*atan(one)
         pi3 = pi/3.
         call rzero(hexbase,16)
         hexbase(1,1) = -2.
c
         do i=0,5
            angle = pi + i*pi3
            k     = 2+i
            hexbase(k,1) = cos(angle)
            hexbase(k,2) = sin(angle)
         enddo
         hexbase(8,1) = hexbase(2,1) ! close loop
         hexbase(8,2) = hexbase(2,2)
      endif
      call outmat(hexbase(2,1),1,8,'hxb x',1)
      call outmat(hexbase(2,2),1,8,'hxb y',1)
c
C
      CALL PRS('Enter (w/mouse) 2 points for framing box.$')

c     IFTMP =IFGRID
c     IFGRID=.FALSE.
 
      call prs('Enter 1st point:$')
         call mouse(xmouse,ymouse,button)
         xmouse = min(xmouse,xphy(1.))
         xmouse = max(xmouse,xphy(0.))
         ymouse = min(ymouse,yphy(1.))
         ymouse = max(ymouse,yphy(0.))

      call prs('Enter 2nd point:$')
         call mouse(xmous2,ymous2,button)
         xmous2 = min(xmous2,xphy(1.))
         xmous2 = max(xmous2,xphy(0.))
         ymous2 = min(ymous2,yphy(1.))
         ymous2 = max(ymous2,yphy(0.))
C     
C     We successfully input 2 points in the build area
C     1) count number of centroids in the bounding box
C     2) if ncntrd=0, find nearest element and anihilate it.
C     
C     Box
      xmax=max(xmouse,xmous2)
      xmin=min(xmouse,xmous2)
      ymax=max(ymouse,ymous2)
      ymin=min(ymouse,ymous2)
      call drwbox(xmin,ymin,xmax,ymax,1)
c
      dxhex = 3*griddx
      dyhex = 3.
      dyhex = sqrt(dyhex)*griddx
      imin = xmin/dxhex
      imax = xmax/dxhex
      jmin = ymin/dyhex
      jmax = ymax/dyhex
      e = nel
      write(6,*) imin,imax,jmin,jmax,' imnx',dxhex
      write(6,*) xmin,xmax,ymin,ymax,' xmnx',dyhex
      do j=jmin,jmax
      do i=imin,imax
         xshift = 0
         yshift = 0
         do k=1,2
            xi = i*dxhex + xshift
            yj = j*dyhex + yshift
            call cmult2(hexscale,hexbase,griddx,16)
            call cadd  (hexscale(1,1),xi,8)
            call cadd  (hexscale(1,2),yj,8)
c
            e = e + 1
            do iv=1,4
               x(iv,e) = hexscale(1+iv,1) ! kth hex, lower half
               y(iv,e) = hexscale(1+iv,2)
            enddo
            long_edge(e) = 3
            call drawel(e)
c
            e = e + 1
            do iv=1,4
               x(iv,e) = hexscale(4+iv,1) ! kth hex, upper half
               y(iv,e) = hexscale(4+iv,2)
            enddo
            long_edge(e) = 1
            call drawel(e)
c
            xshift = dxhex/2
            yshift = dyhex/2
         enddo
      enddo
      enddo
      nel = e
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine hex_transition_3d_e(el,en)  
c
c     Refine 3D hex prisms in transition zone, i.e. the zone between  
c     coarse and fine 3D hex meshes
c
#     include "basics.inc"
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)

c     Element el is partitioned into 5 elements: a, b, c, d, e. Element a
c     replaces original element el. 
c
      parameter (m3=-3,m2=-2,n2=2,n3=3)
      real    rsvx(-3:3) 
      integer rptrx(8,5), rptry(8,5), rptrz(8,5) ! 8 pts, 5 elmts, 3 dir.
      save rptr
c
c     Coordinates for each of a, b, c, d and e (in x, y and z) 
      data rptrx / m2,n2,n3,m3,m2,n2,n3,m3
     $            ,-1,m2,m3,-1,-1,m2,m3,-1
     $            ,m3,n3, 1,-1,m3,n3, 1,-1
     $            ,n2, 1, 1,n3,n2, 1, 1,n3
     $            ,-1, 1, 1,-1,m2,n2,n3,m3  /
c      
      data rptry / -1,-1, 0, 0,-1,-1, 0, 0
     $            ,-1,-1, 0, 1,-1,-1, 0, 1
     $            , 0, 0, 1, 1, 0, 0, 1, 1
     $            ,-1,-1, 1, 0,-1,-1, 1, 0
     $            ,-1,-1, 1, 1,-1,-1, 0, 0  /
c
      data rptrz /  0, 0,n3,n3, 1, 1, 1, 1 
     $            ,-1, 0,n3,-1, 1, 1, 1, 1
     $            ,n3,n3,-1,-1, 1, 1, 1, 1
     $            , 0,-1,-1,n3, 1, 1, 1, 1
     $            ,-1,-1,-1,-1, 0, 0,n3,n3  /

      integer el,en,et
c
      rsvx( 0) = 0
      do i=1,3
        rsvx( i) =  1./i
        rsvx(-i) = -1./i
      end do
c
      do k=0,3
        call copyel(el,en+k)
      enddo
c
      nx = nxm
      ny = nym
      nz = nzm
      call genxyz_e (xp,yp,zp,el,nx,ny,nz)
c
      rrl(3) = 0 
      do k=1,5    ! 5 new elements
        et = en + (k-2)
        if (k.eq.1) et = el
        do i=1,8  ! 8 points
          rrl(1) = rsvx(rptrx(i,k))
          rrl(2) = rsvx(rptry(i,k))
          rrl(3) = rsvx(rptrz(i,k))
          call evalsc (x(i,et),xp,rrl,1)
          call evalsc (y(i,et),yp,rrl,0)
          call evalsc (z(i,et),zp,rrl,0)
        enddo

        if (k.eq.3) then
           ccurve(1,et) = ' '
           ccurve(5,et) = ' '
        endif

        call drawel(et)
      enddo
c
      nel = nel + 4
c
      return
      end
c-----------------------------------------------------------------------              
      subroutine hex_transition_3d ! Refine 3D hexes in transition zone
c
#     include "basics.inc"
      integer e,en,et
c
      integer icalld
      save icalld
      data icalld  /0/
c
      call find_long_edges
      call gencen
c
      nelo = nel
      do e=1,nelo
         en = nel+1
         call rotate_element(e,long_edge(e))
         call hex_transition_3d_e(e,en) ! refine 3D hexes in a zone
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_refine_e(e,en) ! refine 2D hexes in a zone
c
#     include "basics.inc"
      common /chexb/ hexscale(8,2),hexbase(8,2)
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
c
      parameter (m3=-3,m2=-2,n2=2,n3=3)
      real    rsvx(-3:3)
      integer rptr(4,2,4,4)
      save rptr
      data rptr / m2,n2,n3,m3,-1,-1, 0, 0
     $          , -1,m2,m3,-1,-1,-1, 0, 1
     $          , n2, 1, 1,n3,-1,-1, 1, 0
     $          , m3,n3, 1,-1, 0, 0, 1, 1
c
     $          ,  0, 1, 1, 0,m3,m2,n2,n3
     $          , -1, 1, 1, 0,-1,-1,m2,m3
     $          ,  0, 1, 1,-1,n3,n2, 1, 1
     $          , -1, 0, 0,-1,-1,m3,n3, 1
c
     $          , m3,n3,n2,m2, 0, 0, 1, 1
     $          , n3, 1, 1,n2, 0,-1, 1, 1
     $          , -1,m3,m2,-1,-1, 0, 1, 1
     $          , -1, 1,n3,m3,-1,-1, 0, 0
c
     $          , -1, 0, 0,-1,m2,m3,n3,n2
     $          , -1, 0, 1,-1,n2,n3, 1, 1
     $          , -1, 1, 0,-1,-1,-1,m3,m2
     $          ,  0, 1, 1, 0,m3,-1, 1,n3   /

      integer e,en,et
c
      rsvx( 0) = 0
      do i=1,3
         rsvx( i) =  1./i
         rsvx(-i) = -1./i
      enddo
c
      do k=0,2
         call copyel(e,en+k)
      enddo
c
      if (long_edge(e).eq.1) then
         long_edge(en  ) = 4
         long_edge(en+1) = 2
         long_edge(en+2) = 3
      elseif (long_edge(e).eq.2) then
         long_edge(en  ) = 1
         long_edge(en+1) = 3
         long_edge(en+2) = 4
      elseif (long_edge(e).eq.3) then
         long_edge(en  ) = 2
         long_edge(en+1) = 4
         long_edge(en+2) = 1
      elseif (long_edge(e).eq.4) then
         long_edge(en  ) = 3
         long_edge(en+1) = 1
         long_edge(en+2) = 2
      endif
c
      do j=0,maxfld
         cbc(long_edge(e),en+2,j) = ' '
      enddo
c
      nx = nxm
      ny = nym
      nz = 1
      call genxyz_e (xp,yp,zp,e,nx,ny,nz)
c
      rrl(3) = 0.
      do k=1,4     ! 4 new elements
         et = en + (k-2)
         if (k.eq.1) et = e
         do i=1,4  ! 4 points
            rrl(1) = rsvx(rptr(i,1,k,long_edge(e)))
            rrl(2) = rsvx(rptr(i,2,k,long_edge(e)))
            call evalsc (x(i,et),xp,rrl,1)
            call evalsc (y(i,et),yp,rrl,0)
                         z(i,et) = 0
         enddo
         call drawel(et)
      enddo
      nel = nel + 3
c
      return
      end
c-----------------------------------------------------------------------
      subroutine hexagon_refine ! refine 2D hexes in a zone
c
#     include "basics.inc"
      integer e,en,et
c
      integer elist(nelm),marked(nelm)
c
      integer icalld
      save    icalld
      data    icalld  /0/
c
      call get_el2(nlist,elist)
c
c     This is ugly (potentially ON^2), but effective
c
      k  = 0
      kk = 0
   10 continue
c
      call find_long_edges
      call gencen
      call gen_neigh
c
      k  = k+1
      kk = kk+1
      nnew = 0
      do e=1,nel
         if (elist(e).ne.0) then ! pending request
            l  = long_edge(e)
            en = neighb(l,e)
            write(6,*) nnew,en,e,l,' Enew'
            if (en.gt.0) then
               ln = long_edge(en)
               et = neighb(ln,en)
               if (et.eq.e) then ! reciprocal pair, ok to refine
                  elist(e ) = -1
                  elist(en) = -1
               elseif (elist(en).eq.0) then
                  nnew = nnew+1
                  elist(en) = 1
               endif
            else
               elist(e) = -1
            endif
         endif
      enddo
      write(6,*) 'KPASS:',k,nnew,kk
      if (nnew.gt.0) goto 10
      k  = 0
c
      nelo = nel
      do e=1,nelo
         if (elist(e).lt.0) then
            en = nel+1
            call hexagon_refine_e(e,en) ! refine 2D hexes in a zone
            elist(e) = 0
         endif
      enddo
c
      icount = 0
      do e=1,nelo
         if (elist(e).ne.0) icount = icount+1
      enddo
      if (icount.gt.0) goto 10
c
      call redraw_mesh
c
      return
      end
c-----------------------------------------------------------------------
      subroutine rotate_element(e,kface)
c
c     Rotate element e so that face kface is mapped to face 1
c
c     We assume that kface  is in [1,4]
c
#     include "basics.inc"
      integer e,et
c
      if (kface.eq.1) return
      if (kface.lt.1.or.kface.gt.4) then
         call prsii('Cannot rotate element. kface out of range.$'
     $              ,e,kface)
         return
      endif
c
      et = nelm-1
      call copyel(e,et)
c
      do it=1,4
c
         i = it-(kface-1)
         if (i.lt.1) i = i + 4
         x(i  ,e) = x(it,et)
         y(i  ,e) = y(it,et)
         z(i  ,e) = z(it,et)
         x(i+4,e) = x(it+4,et)
         y(i+4,e) = y(it+4,et)
         z(i+4,e) = z(it+4,et)
c
         do ifld=0,maxfld
            cbc(i,e,ifld) = cbc(it,et,ifld)
            call copy(bc(1,i,e,ifld),bc(1,it,et,ifld),5)
         enddo
         ibc(i,e,ifld) = ibc(it,et,ifld)

         do ii=1,4
            sides(e,i,ii)=sides(et,it,ii)
         enddo

         ccurve(i  ,e) = ccurve(it  ,et)
         call copy(curve(1,i,e),curve(1,it,et),5)
 
         ccurve(i+4,e) = ccurve(it+4,et)
         call copy(curve(1,i+4,e),curve(1,it+4,et),5)

      enddo

      do k=5,6
         if (ccurve(k,et).eq.'s') then
            ccurve(k,e) = ccurve(k,et)
            call copy(curve(1,k,e),curve(1,k,et),5)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine get_el2(nlist,elist)
c
c     Grab a set of elements
c
#     include "basics.inc"
      integer elist(1),e,tlist(nelm)
c
      call get_els(nlist,tlist)
c
c     convert list to flags
c
      call izero  (elist,nelm)
      do i=1,nlist
         e = tlist(i)
         elist(e) = 1
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_els(nlist,elist)
c
c     Grab a set of elements
c
#     include "basics.inc"
      integer elist(1),e
      logical iftmp
C
      CALL PRS('Enter (w/mouse) 2 points for framing box.$')

      IFTMP =IFGRID
      IFGRID=.FALSE.
 
      call prs('Enter 1st point:$')
      call mouse(xmouse,ymouse,button)
      if (xmouse.lt.xphy(1.)) then
         xmouse = min(xmouse,xphy(1.))
         xmouse = max(xmouse,xphy(0.))
         ymouse = min(ymouse,yphy(1.))
         ymouse = max(ymouse,yphy(0.))

         call prs('Enter 2nd point:$')
            call mouse(xmous2,ymous2,button)
            xmous2 = min(xmous2,xphy(1.))
            xmous2 = max(xmous2,xphy(0.))
            ymous2 = min(ymous2,yphy(1.))
            ymous2 = max(ymous2,yphy(0.))
C     
C        We successfully input 2 points in the build area
C        1) count number of centroids in the bounding box
C        2) if ncntrd=0, find nearest element and anihilate it.
C     
C        Box
         xmax=max(xmouse,xmous2)
         xmin=min(xmouse,xmous2)
         ymax=max(ymouse,ymous2)
         ymin=min(ymouse,ymous2)
         call drwbox(xmin,ymin,xmax,ymax,1)
      else
         call prs('Enter point (drop mouse!) (99,99 for file input):$')
         call rerr(xmin,ymin)
         xmax = xmin
         ymax = ymin
         if (ymin.eq.99.and.ymax.eq.99) then
            call get_els_bc(nlist,elist)
            return
         endif
      endif
c
      nlist = 0
      do e=1,nel
         if ( xmin.le.xcen(e) .and. xcen(e).le.xmax .and.
     $        ymin.le.ycen(e) .and. ycen(e).le.ymax ) then
            nlist = nlist + 1
            elist(nlist) = e
         endif
      enddo
      if (nlist.eq.0) then ! find closest element
         xm = .5*(xmax+xmin)
         ym = .5*(ymax+ymin)
         dx = xm-xcen(1)
         dy = ym-ycen(1)
         dm = dx*dx + dy*dy
         im = 1
         do e=1,nel
            dx = xm-xcen(e)
            dy = ym-ycen(e)
            d2 = dx*dx + dy*dy
            if (d2.lt.dm) then
               dm = d2
               im = e
            endif
         enddo
         nlist = 1
         elist(1) = im
      endif
c
      IFGRID=IFTMP
c
      return
      end
c-----------------------------------------------------------------------
      subroutine get_els_bc(nlist,elist)
c
c     Grab a set of elements
c
#     include "basics.inc"
      integer elist(1),e
c
      ifld = 2
      if (ifflow) ifld = 1
c
      nlist = 0
      do e=1,nel
         do k=1,4
            if (cbc(k,e,ifld).eq.'v'.or.cbc(k,e,ifld).eq.'t') then
               nlist = nlist+1
               elist(nlist) = e
               goto 10
            endif
         enddo
   10    continue
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine fix_m_curve_old(e)

c     Assign / repair curve-side info for edges

#     include "basics.inc"
      integer e

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation


      nedge = 4 + 8*(ndim-2)
      do ic=1,nedge
         ccurve (ic,e) = 'm'   ! Midside node (now default??)
         curve(1,ic,e) = x27(eindx(ic),e)
         curve(2,ic,e) = y27(eindx(ic),e)
         curve(3,ic,e) = z27(eindx(ic),e)
      enddo
c     if (e.eq.2) call out27(x27(1,e),y27(1,e),z27(1,e),e,'fixm')
c     if (e.eq.2) write(6,*) 'stop in fix_m_curve',e
c     if (e.eq.2) call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_m_curve(e)

c     Assign / repair curve-side info for edges

#     include "basics.inc"

      real xyz(3,3)

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1
     $           , 19,20,21,   21,24,27,   27,26,25,   25,22,19
     $           ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /
      real len
      integer e,edge

      tol   = 1.e-4
      tol2  = tol**2
      nedge = 4 + 8*(ndim-2)

      do i=1,nedge
         do j=1,3
            xyz(1,j)=x27(e3(j,i),e)
            xyz(2,j)=y27(e3(j,i),e)
            xyz(3,j)=z27(e3(j,i),e)
         enddo

c        write(6,*)
c        write(65,*)
c        write(6,6) (e,i,xyz(1,k),xyz(2,k),k=1,3)
c        write(65,6) (e,i,xyz(1,k),xyz(2,k),k=1,3)

         len = 0.
         h   = 0.
         do j=1,ndim
            xmid = .5*(xyz(j,1)+xyz(j,3))
            h    = h   + (xyz(j,2)-xmid)**2
            len  = len + (xyz(j,3)-xyz(j,1))**2
         enddo
         ht = tol2*len
         if (h.gt.ht) then
            ccurve(i,e) = 'm'
            call copy(curve(1,i,e),xyz(1,2),ndim)
         else
            ccurve(i,e) = ' '
            call rzero(curve(1,i,e),ndim)
         endif

c        write(6,*) 'ccurve: ',i,e,' ',ccurve(i,e),' ',h,ht
c        write(6,6) (e,i,xyz(1,k),xyz(2,k),k=1,3)

c        write(66,*)
c        write(66,6) (e,i,xyz(1,k),xyz(2,k),k=1,3)
c  6     format(i9,i4,1p2e18.7)
         

      enddo

      call drawel(e)
c     call prs('continue?$')
c     call res(ans,1)

c     if (e.eq.2) call out27(x27(1,e),y27(1,e),z27(1,e),e,'fixm')
c     if (e.eq.2) write(6,*) 'stop in fix_m_curve',e
c     if (e.eq.2) call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine s2vec(x,y,z,u,v,w,r)  ! Project onto r*S2

      x = u
      y = v
      z = w

      zl = x*x+y*y+z*z

      if (zl.gt.0) then
         zl = sqrt(zl)
         x  = x*r/zl
         y  = y*r/zl
         z  = z*r/zl
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine make_top(xl,yl,zl,r)
#     include "basics.inc"

      real xl(3,3,3),yl(3,3,3),zl(3,3,3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      integer e,f

      e = nel+1
      call copyel(nel,e)

      call copy (x(1,e),x(5,nel),4)
      call copy (y(1,e),y(5,nel),4)
      call copy (z(1,e),z(5,nel),4)
      call s2vec(x(5,e),y(5,e),z(5,e),0.,0.,1.,r)  ! Project onto r*S2
      call s2vec(x(6,e),y(6,e),z(6,e),1.,0.,1.,r)
      call s2vec(x(7,e),y(7,e),z(7,e),1.,1.,1.,r)
      call s2vec(x(8,e),y(8,e),z(8,e),0.,1.,1.,r)

      call blank (ccurve(1,e),12)
      call rzero (curve(1,1,e),72)

      call chcopy(ccurve(1,e) ,ccurve(5,nel) , 4)  ! Top 4 edges of nel -->
      call copy  (curve(1,1,e),curve(1,5,nel),24)  ! bottom 4 of e

      ccurve(6,e)  = 's'   ! sphere
      curve(4,6,e) = r     ! radius r

c
c     Clean up a few details
c
      letapt(e)='a'

      do ifld=1,maxfld
      do f=5,6
         cbc(f,e,ifld) = 'W  '    ! totally arbitrary default
         if (f.eq.5.or.f.eq.4.or.f.eq.5) cbc(f,e,ifld) = '   '
         call rzero(bc(1,f,e,ifld),5)
         ibc(f,e,ifld) = 0
      enddo
      enddo

      nel = e

      return
      end
c-----------------------------------------------------------------------
      subroutine make_core(xl,yl,zl,r)
#     include "basics.inc"

      real xl(3,3,3),yl(3,3,3),zl(3,3,3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      integer e,f

      r = 0.5
      a = 0.120   ! delta (1,0,0) radius
      b = .90*a   ! delta (1,1,0) radius
      c = .85*a   ! delta (1,1,1) radius

      ra = r-a
      rb = r-2*b+a  ! .5*(rb+ra) = .5*(r+r-2*b+a-a) = r-b
      rc = r-2*c+a

      e = nel+1

      call rzero(x(1,e),8)
      call rzero(y(1,e),8)
      call rzero(z(1,e),8)

      call blank(ccurve(1,e),12)
      call rzero(curve(1,1,e),72)

      call s2vec(x(1,e),y(1,e),z(1,e),0.,0.,0.,0.)  ! Project onto ra*S2
      call s2vec(x(2,e),y(2,e),z(2,e),1.,0.,0.,ra)  ! Project onto ra*S2
      call s2vec(x(3,e),y(3,e),z(3,e),1.,1.,0.,rb)  ! Project onto ra*S2
      call s2vec(x(4,e),y(4,e),z(4,e),0.,1.,0.,ra)  ! Project onto ra*S2
      call s2vec(x(5,e),y(5,e),z(5,e),0.,0.,1.,ra)  ! Project onto ra*S2
      call s2vec(x(6,e),y(6,e),z(6,e),1.,0.,1.,rb)  ! Project onto ra*S2
      call s2vec(x(7,e),y(7,e),z(7,e),1.,1.,1.,rc)  ! Project onto ra*S2
      call s2vec(x(8,e),y(8,e),z(8,e),0.,1.,1.,rb)  ! Project onto ra*S2

c     call outmat(x(1,e),2,4,'x(e) ',e)
c     call outmat(y(1,e),2,4,'y(e) ',e)
c     call outmat(z(1,e),2,4,'z(e) ',e)

      call xyzlin (xl,yl,zl,3,3,3,e)

c     Curve the interior connector edges

      ccurve( 2,e) = 'm'
      ccurve( 3,e) = 'm'
      ccurve( 5,e) = 'm'
      ccurve( 6,e) = 'm'
      ccurve( 7,e) = 'm'
      ccurve( 8,e) = 'm'
      ccurve(10,e) = 'm'
      ccurve(11,e) = 'm'
      ccurve(12,e) = 'm'

      w0 = 1.10 ! weight on curve
      w1 = 1-w0
      wp = 1 + .2*a/r

      nedge = 12
      do k=1,nedge
         if (ccurve(k,e).eq.'m') then
            j = eindx(k)
            call s2vec(xs,ys,zs,xl(j,1,1),yl(j,1,1),zl(j,1,1),r )
c           curve(1,k,e) = w0*xs + w1*xl(j,1,1)
c           curve(2,k,e) = w0*ys + w1*yl(j,1,1)
c           curve(3,k,e) = w0*zs + w1*zl(j,1,1)
            curve(1,k,e) = wp*xl(j,1,1)
            curve(2,k,e) = wp*yl(j,1,1)
            curve(3,k,e) = wp*zl(j,1,1)
         endif
      enddo

      call linquad(xl,yl,zl,3,3,3,e)

c
c     Clean up a few details
c
      letapt(e)='a'

      do ifld=1,maxfld
      do f=1,6
         cbc(f,e,ifld) = 'W  '    ! totally arbitrary default
         if (f.eq.1.or.f.eq.4.or.f.eq.5) cbc(f,e,ifld) = 'SYM'
         call rzero(bc(1,f,e,ifld),5)
         ibc(f,e,ifld) = 0
      enddo
      enddo

      nel = e

      return
      end
c-----------------------------------------------------------------------
      subroutine filled_sphere_oct ! Builds sphere or hemisphere of diameter 1
#     include "basics.inc"

      common /sphoct/ sphctr(3),normal(3),xl(27),yl(27),zl(27)
      real normal
      integer e

      call make_core(xl,yl,zl,r)

      call make_top (xl,yl,zl,r)

c     Copy top and rotate about (1,1,1):

      call rone(normal,3)
      do k=1,2
         e = nel+k
         call copyel(nel,e)
         angle = 120*k
         call rotate_el_vec(e,normal,angle)
      enddo
      nel = e

c     call prs('Enter number of elements radial direction:$')
c     call rei(ne_rad)
c     if (ne_rad.le.2) return
c     if (ne_rad.le.4) then
c        nro = 1
c        nrc = ne_rad - nro
c     else
c        nrc = (ne_rad+.1)/1.5
c        nro = ne_rad-nrc
c        nro = max(2,nro)
c        nrc = ne_rad-nro
c     endif
c     ro  = .8
c     ro  = .6

      call prs('Enter number of elements in core and shell:$')
      call reii(nrc,nro)
      call prs('Enter ratio in outer shell:$')
      call rer(ro)

      r1  = 1
      e   = 1
      call msplite(e,nrc,nrc,nrc,r1,r1,r1)    ! Uniform refinement of core

      do e=2,4
         call msplite(e,nrc,nrc,nro,r1,r1,ro) ! Refine outer shell
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine make_top_2d(xl,yl,zl,r)
#     include "basics.inc"

      real xl(3,3,3),yl(3,3,3),zl(3,3,3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      integer e,f

      do ipass=1,2
        e = nel+ipass
        call copyel(nel,e)
        if (ipass.eq.1) then
          x(1,e) = x(3,nel)
          x(2,e) = x(2,nel)
          y(1,e) = y(3,nel)
          y(2,e) = y(2,nel)
          call s2vec(x(3,e),y(3,e),z(3,e),1.,0.,0.,r)  ! Project onto r*S2
          call s2vec(x(4,e),y(4,e),z(4,e),1.,1.,0.,r)  ! Project onto r*S2
        else
          x(1,e) = x(4,nel)
          x(2,e) = x(3,nel)
          y(1,e) = y(4,nel)
          y(2,e) = y(3,nel)
          call s2vec(x(3,e),y(3,e),z(3,e),1.,1.,0.,r)  ! Project onto r*S2
          call s2vec(x(4,e),y(4,e),z(4,e),0.,1.,0.,r)  ! Project onto r*S2
        endif

        call blank (ccurve(1,e),12)
        call rzero (curve(1,1,e),72)

        ccurve(1,e) = ccurve(1+ipass,nel)
        call copy  (curve(1,1,e),curve(1,1+ipass,nel),6)

        f=3
        ccurve(f,e)  = 'C'   ! circle
        curve(1,f,e) = r     ! radius r

c
c       Clean up a few details
c
        letapt(e)='a'

        do ifld=1,maxfld
           f=3
           cbc(f,e,ifld) = 'W  '    ! totally arbitrary default
           call rzero(bc(1,f,e,ifld),5)
           ibc(f,e,ifld) = 0
        enddo
      enddo

      nel = e

      return
      end
c-----------------------------------------------------------------------
      subroutine make_core_2d(xl,yl,zl,r)
#     include "basics.inc"

      real xl(3,3,3),yl(3,3,3),zl(3,3,3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      integer e,f

      r = 0.5
      a = 0.120   ! (1,0,0) radius
      b = .90*a   ! (1,1,0) radius
      c = .85*a   ! (1,1,1) radius

      ra = r-a
      rb = r-2*b+a  ! .5*(rb+ra) = .5*(r+r-2*b+a-a) = r-b
      rc = r-2*c+a

      e = nel+1

      call rzero(x(1,e),8)
      call rzero(y(1,e),8)
      call rzero(z(1,e),8)

      call blank(ccurve(1,e),12)
      call rzero(curve(1,1,e),72)

      call s2vec(x(1,e),y(1,e),z(1,e),0.,0.,0.,0.)  ! Project onto ra*S2
      call s2vec(x(2,e),y(2,e),z(2,e),1.,0.,0.,ra)  ! Project onto ra*S2
      call s2vec(x(3,e),y(3,e),z(3,e),1.,1.,0.,rb)  ! Project onto ra*S2
      call s2vec(x(4,e),y(4,e),z(4,e),0.,1.,0.,ra)  ! Project onto ra*S2

      call xyzlin (xl,yl,zl,3,3,3,e)

c     Curve the interior connector edges
      ccurve( 2,e) = 'm'
      ccurve( 3,e) = 'm'

      w0 = 1.10 ! weight on curve
      w1 = 1-w0

      nedge = 4
      do k=1,nedge
         if (ccurve(k,e).eq.'m') then
            j = eindx(k)
            call s2vec(xs,ys,zs,xl(j,1,1),yl(j,1,1),zl(j,1,1),ra)
            curve(1,k,e) = w0*xs + w1*xl(j,1,1)
            curve(2,k,e) = w0*ys + w1*yl(j,1,1)
         endif
      enddo

      call linquad(xl,yl,zl,3,3,3,e)

c
c     Clean up a few details
c
      letapt(e)='a'

      do ifld=1,maxfld
      do f=1,4
         cbc(f,e,ifld) = 'W  '    ! totally arbitrary default
         if (f.eq.1.or.f.eq.4.or.f.eq.5) cbc(f,e,ifld) = 'SYM'
         call rzero(bc(1,f,e,ifld),5)
         ibc(f,e,ifld) = 0
      enddo
      enddo

      nel = e

      return
      end
c-----------------------------------------------------------------------
      subroutine filled_cyl_quadrant ! Builds filled cylinder quadrant
#     include "basics.inc"

      common /sphoct/ sphctr(3),normal(3),xl(27),yl(27),zl(27)
      real normal
      integer e

      call make_core_2d (xl,yl,zl,r)
      call make_top_2d  (xl,yl,zl,r)

      call vertadj

c     call prs('Enter number of elements radial direction:$')
c     call rei(ne_rad)
c     if (ne_rad.le.2) return
c     if (ne_rad.le.4) then
c        nro = 1
c        nrc = ne_rad - nro
c     else
c        nrc = (ne_rad+.1)/1.5
c        nro = ne_rad-nrc
c        nro = max(2,nro)
c        nrc = ne_rad-nro
c     endif
c     ro  = .8
c     ro  = .6

      call prs('Enter number of elements in core and shell:$')
      call reii(nrc,nro)
      call prs('Enter ratio in outer shell:$')
      call rer(ro)

      r1  = 1
      e   = 1
      call msplite(e,nrc,nrc,1,r1,r1,0.)    ! Uniform refinement of core

      do e=2,3
         call msplite(e,nrc,nro,1,r1,ro,0.) ! Refine outer shell
      enddo

      call vertadj

      return
      end
c-----------------------------------------------------------------------
      function blend_circ_in_box(x,y,r0,x0i,x1i,y0i,y1i)
      real p(2),v(2,5)
      real o(2)
      save o
      data o / 0. , 0. /

      blend_circ_in_box = 1.

      r=x*x+y*y
      if (r.le.r0*r0) return

      ey = 1.e-6*(y1i-y0i)   ! Put a slight tolerance on box size
      ex = 1.e-6*(x1i-x0i)
      x0 = x0i-ex
      x1 = x1i+ex
      y0 = y0i-ey
      y1 = y1i+ey

      r = sqrt(r)
      t = atan2(y,x)

      p(1)=x
      p(2)=y

      do i=1,5
         v(1,i)=x0
         v(2,i)=y0
      enddo
      v(1,2)=x1
      v(1,3)=x1
      v(2,3)=y1
      v(2,4)=y1
      call cmult(v,2.0,10) ! Make domain bigger for "in_triangle" check

      ks=0
      do k=1,4
         in = in_triangle(p,o,v(1,k),v(1,k+1)) ! triangle: [ o vk vk1 ]
         if (in.gt.0) then
           if (k.eq.1) then       !  Lower y boundary
              yc = r0*sin(t)
              b  = (y-y0)/(yc-y0)
           elseif (k.eq.2) then   ! Right x boundary
              xc = r0*cos(t)
              b  = (x-x1)/(xc-x1)
           elseif (k.eq.3) then   ! Upper y boundary
              yc = r0*sin(t)
              b  = (y-y1)/(yc-y1)
           else                   !  Left x boundary
              xc = r0*cos(t) 
              b  = (x-x0)/(xc-x0) 
           endif
           blend_circ_in_box = max(b,0.)
           return
         endif
      enddo
      t = sqrt(-t)
      t = sqrt(-t)

      return
      end
c-----------------------------------------------------------------------
      function in_triangle(p,a,b,c)
      real p(2),a(2),b(2),c(2)
      real v0(2),v1(2),v2(2),invdenom

      call sub3(v0,c,a,2) ! v0 = c-a
      call sub3(v1,b,a,2) ! v1 = b-a
      call sub3(v2,p,a,2) ! v2 = p-a

      d00 = vlsc2(v0,v0,2)
      d01 = vlsc2(v0,v1,2)
      d02 = vlsc2(v0,v2,2)
      d11 = vlsc2(v1,v1,2)
      d12 = vlsc2(v1,v2,2)

c     Compte barycentric coordinates

      invdenom = 1./(d00*d11-d01*d01)
      u=(d11*d02-d01*d12)*invdenom
      v=(d00*d12-d01*d02)*invdenom

      in_triangle = 0
      if (u.ge.0.and.v.ge.0.and.(u+v).lt.1) in_triangle = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine x27_to_e
#     include "basics.inc"

      integer e

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      nedge = 4 + 8*(ndim-2)

      do e=1,nel

         call q_to_neklin  (x(1,e),1,x27(1,e),if3d)
         call q_to_neklin  (y(1,e),1,y27(1,e),if3d)
         call q_to_neklin  (z(1,e),1,z27(1,e),if3d)

         do kk=1,nedge
            ccurve(kk,e)='m'
            jj = eindx(kk)
            curve(1,kk,e) = x27(jj,e)
            curve(2,kk,e) = y27(jj,e)
            curve(3,kk,e) = z27(jj,e)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine std_transform(t0,t1)
#     include "basics.inc"

c     Transforms existing 2D or 3D rectangular domain to a
c     circular arc segment
c
c     Here, we judge the length at y=0 -- this arc length will be
c     preserved over the desired theta range (t0,t1)
c
c     This routine uses x27 and converts all circles and spheres to midside node.

c     Results are best if incoming domain is close to desired dimensions.

      call gencen  ! Generate x27

      n = 27*nel
      xmin=glmin(x27,n)
      xmax=glmax(x27,n)
      xbar=(xmax+xmin)/2
      ymin=glmin(y27,n)
      ymax=glmax(y27,n)
      ybar=(ymax+ymin)/2
      zmin=glmin(z27,n)
      zmax=glmax(z27,n)
      zbar=(zmax+zmin)/2

      arclength = xmax-xmin
      r1 = arclength/(t1-t0)
      write(6,*) xmin,xmax, ' xmin,xmax'
      write(6,*) t0,t1,' t0,t1'
      write(6,*) arclength,r1,' arclngth,r1'


      do i=1,n

         xx = x27(i,1)
         yy = y27(i,1)
         zz = z27(i,1) ! / (i-1)

         rr = yy + r1
         th = xx / r1

         x27(i,1) = rr*sin(th)
         y27(i,1) = rr*cos(th)

      enddo

      call x27_to_e  ! Converts all edges to 'm'
      call redraw_mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine get_x27max(xmin,xmax,ymin,ymax,zmin,zmax)
#     include "basics.inc"

      n = 27*nel
      xmin=glmin(x27,n)
      xmax=glmax(x27,n)
      ymin=glmin(y27,n)
      ymax=glmax(y27,n)
      zmin=glmin(z27,n)
      zmax=glmax(z27,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine make_arc_circ
#     include "basics.inc"
      integer e

c     Transforms existing 2D or 3D rectangular domain to a circular arc


      call gencen
      do e=1,nel
         call fix_m_curve(e)
      enddo
      ifmid  = .true.
      ifcstd = .true.

      call get_x27max(xmin,xmax,ymin,ymax,zmin,zmax)

      call prs('Standard (0) or Circle preserving (1)?$')
      call rei(icirc)

c     call prs('Enter protected radius:$')
c     call rer(r0)

c     call prs('Enter desired angle (deg.):$')
c     call rer(dth)
      dth = 30
      dth = pi*dth/180.

      t0  = -dth/2
      t1  =  dth/2

      y0  = 1.9
      y1  = 6.
      r1  = (y0+y1)/2

      arclength_b = r1*(t1-t0)
      arclength_a = xmax-xmin
      r1          = 0.5*arclength_a/(t1-t0)

      write(6,*) arclength_a,arclength_b,' arclength a,b'
      write(6,*) xmin,xmax, ' xmin,xmax'
      write(6,*) t0,t1,' t0,t1'
      write(6,*) arclength,r1,' arclength,r1'

      if (icirc.eq.0) then ! Standard Map
         call std_transform(t0,t1)
      else
         call circle_in_arc_transform_0(r1)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine circle_in_arc_transform_0(r1)
#     include "basics.inc"

c     NOTE: Must call gencen first!

c     Transforms existing 2D or 3D rectangular domain
c     to a circular arc segment.
c
c     Origin, center of circle, is mapped to radius r1 and angle t=0.
c
c     Domain is transformed to arc segment [t0,t1], with angle t
c     measured _clockwise_ from the y axis, and major radii spanning [y0,y1].

c     This routine uses x27 and converts all edges to midside node.

c     Results are best if incoming domain is close to desired dimensions.

c     call gencen  ! Generate x27

      open(76,file='x.x')

      do i=1,n
         xx = x27(i,1)
         yy = y27(i,1)

         rad2  = xx*xx + yy*yy       ! Circle-preserving transform
         ycirc = yy - .5*xx*xx/r1
         argx  = rad2 - ycirc**2
         xcirc = sqrt(argx)
         if (xx.lt.0) xcirc = -xcirc
         ycirc = ycirc + r1

         x27(i,1) = xcirc
         y27(i,1) = ycirc

         rad = sqrt(rad2)
         write(76,1) xx,yy,xcirc,ycirc,rad,r1
   1     format(1p6e12.4)
      enddo
      close(76)

      call x27_to_e  ! Converts all edges to 'm'
      call redraw_mesh

      call get_x27max(xmin,xmax,ymin,ymax,zmin,zmax)
      dx  = xmax-xmin
      r1n = ymax

      theta = 180*2*atan2(dx,r1n)/pi
      thet0 = 180*2*atan2(dx,r1 )/pi

      write(6,*) xmin,xmax,' new xmin,xmax'
      write(6,*) r1n,r1,' new r1n,r1'
      write(6,*) theta,thet0,' new theta,thet0'

      return
      end
c-----------------------------------------------------------------------
