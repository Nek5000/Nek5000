c-----------------------------------------------------------------------
      subroutine newel(xmouse,ymouse,button,ierr) ! add new element

#     include "basics.inc"
      dimension iobjs(8)
      integer e

      e = nel+1

      do icorn=1,4

         goto 100 ! else, error
   90       call prs('Error reading input.  To input corner, use $')
            call prs('mouse$')
  100    continue

         if (icorn.eq.1) call prs(' Corner I  >$')
         if (icorn.eq.2) call prs(' Corner II >$')
         if (icorn.eq.3) call prs(' Corner III>$')
         if (icorn.eq.4) call prs(' Corner IV >$')

         if (xmouse.le.xphy(1.0) .and. icorn.eq.1) then
c           Special Kludge for "ADD ELEMENT"
c           First element corner already input
         else
            xold=xmouse
            yold=ymouse
            call mouse(xmouse,ymouse,button)
            if (xmouse.eq.xold.and.ymouse.eq.yold) goto 90
         endif

         if (xmouse.gt.xphy(1.0)) then ! look for a keypad input
            call prs('Switching to keyboard inputs$')
            call sgvis(11,0)
            call sgvis(12,1)
            call prs('Enter x and y with keyboard:$')
            call rerr(xmouse,ymouse)
            button='LEFT'
            call sgvis(12,0) ! turn off keypad hilite
            call sgvis(11,1)
         endif

         write(6,*) icorn,e,iobjct,' IOBJ AAA'
         if (button.eq.'LEFT') then ! std. input:  new unique corner
            if (iobjct.eq.0) xmouse=round(xmouse) ! round input value - 7-25-90 pff
            if (iobjct.eq.0) ymouse=round(ymouse)
         elseif (button.eq.'RIGHT') then ! latch to closest vertex
            call blatch(xmouse,ymouse,icmin,iemin)
            write(6,*) xmouse,ymouse,icmin,iemin,' COOR A'
            write(6,*) x(icmin,iemin),y(icmin,iemin),e,' COOR B'
            write(6,*) icorn,e,iobjct,' IOBJ BBB'
         else
            call prs('Error reading mouse input$')
            call beep
            goto 90     ! Go back to query
         endif
         write(6,*) icorn,e,iobjct,' IOBJ CCC'

         x(icorn,e)=xmouse
         y(icorn,e)=ymouse
         iobjs(icorn)=iobjct ! identifies if this corner latched to object

         if (if3d) x(icorn+4,e)=x(icorn,e)
         if (if3d) y(icorn+4,e)=y(icorn,e)

         call color(3)
         call fillp(-15)

         xs=x(icorn,e)
         ys=y(icorn,e)
         call prsrr('xs,ys:$',xs,ys)
         call prsrr('xs1ys1$',xs1,ys1)
         call prsiii('e,crn:$',e,icorn,in)
 
         if (icorn.eq.1) then
            call beginb(xs,ys)
            xs1=xs
            ys1=ys
            if (in.eq.5)call rubber
         else
            call drawc(xs,ys)
         endif

         if (icorn.eq.4) then ! CLOSE OUT

            call offrub ! Turn off rubberbanding
            call drawc(xs1,ys1)
            call endp

            call clockw_chk(e,ierr) ! counter-clockwise check
            if (ierr.ne.0) then
                call prsi('ERROR in entering element.  Re-enter.$',e)
                return
            endif

         endif
      enddo

      xcen(e)=(x(1,e)+x(2,e)+x(3,e)+x(4,e))/4.  ! Label Corners
      ycen(e)=(y(1,e)+y(2,e)+y(3,e)+y(4,e))/4.
      iletap=iletap+1
      if (iletap.le.122)letapt(e)=char(iletap)
      if (iletap.gt.122)letapt(e)=' '

      if (if3d) then ! New elements have their Z-coordinate put in here
         zceil=0.0
         do 620 i=1,numapt(e)
            zceil=zceil+height(i)
  620    continue
         zfloor=zceil-height(numapt(e))
         do 630 ic=1,8
           if (ic.le.4)z(ic,e)=zfloor
           if (ic.gt.4)z(ic,e)=zceil
  630    continue
      endif

      call set_auto_curve(e,iobjs) ! curve sides for polar coords or objects

      call mkside_e(e) ! Define sides' midpoints
      call drawel(e)   ! Element all set , draw it.

      if (if3d) then
         call sortel
         do i=1,e    ! Figure out which to draw
            if (isrt(i).eq.e) ibegin=i
         enddo
         do i=ibegin,e
            call drawis(isrt(i))
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_auto_curve(e,iobjs) ! set curve sides

#     include "basics.inc"
      dimension iobjs(8)
      integer e


      call blank(ccurve(1,e)  ,12)
      call rzero( curve(1,1,e),72)

      call chk_polar    (e)
      call chk_obj      (e,iobjs)
      call chk_neighbor (e)

      call drawel(e)

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_polar (e)
#     include "basics.inc"
      logical ifcurv
      integer e

      if (.not.ifgrdp) return

      radusc = 0.0   ! Characteristic length
      do 800 ic=1,4
         dx = x(ic,e)-gridxp
         dy = y(ic,e)-gridyp
         radusc = radusc + sqrt(dx**2+dy**2)
  800 continue
      radusc = radusc/4.0

c     check for any curve sides, if yes, delete straight sided image.
c     (two steps required to delete proper shaped element)

      ifcurv=.false.
      do ic=1,4
         jc = mod1(ic+1,4)
         dx0=x(ic,e)-gridxp
         dy0=y(ic,e)-gridyp
         dx1=x(jc,e)-gridxp
         dy1=y(jc,e)-gridyp
         radus0 = sqrt(dx0**2+dy0**2)
         radus1 = sqrt(dx1**2+dy1**2)
         drad = (radus1-radus0)/(radus1+radus0)
         if (abs(drad).lt.0.0001) ifcurv=.true.
      enddo

      if (ifcurv) then
         call drawel(-e)
C        Now curve side and update element image.
         do ic=1,4
            jc = mod1(ic+1,4)
            if (mod(ic,4).eq.0) jc = ic-3
            dx0=x(ic,e)-gridxp
            dy0=y(ic,e)-gridyp
            dx1=x(jc,e)-gridxp
            dy1=y(jc,e)-gridyp
            radus0 = sqrt(dx0**2+dy0**2)
            radus1 = sqrt(dx1**2+dy1**2)
            drad = (radus1-radus0)/(radus1+radus0)

            if (abs(drad).lt.0.0001) then ! curve side
               if (radus0.lt.radusc) radus0 = -radus0
               curve(1,ic,e)=radus0
               ccurve(ic,e)='C'
               if (if3d) then
                  ic4 = ic+4
                  curve(1,ic4,e)=radus0
                  ccurve(ic4,e)='C'
               endif
            endif

         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_obj (e,iobjs)
#     include "basics.inc"
      integer e,iobjs(4)
      logical ifcurv

      xmn=x(1,e)
      ymn=y(1,e)
      xmx=x(1,e)
      ymx=y(1,e)

      do i=1,4
         xmn=min(xmn,x(i,e))
         ymn=min(ymn,y(i,e))
         xmx=max(xmx,x(i,e))
         ymx=max(ymx,y(i,e))
      enddo
      diam = (xmx-xmn)+(ymx-ymn)

      do i=1,4
         tol  = 0.01*diam
         iobjs(i)=0
         do iobj=1,nobjs
            if (ifobjg(iobj)) then ! Active object
               call latchob(xo,yo,x(i,e),y(i,e),0.,dist2,ko,io,iobj)
               if (dist2.gt.0) dist2=sqrt(dist2)
               if (dist2.lt.tol) then
                  iobjs(i)=iobj
                  tol = dist2
               endif
          write(6,8) i,iobj,dist2,xo    ,yo    ,' dist2,xo,yo',iobjs(i)
          write(6,8) i,iobj,tol  ,x(i,e),y(i,e),' tol  ,xe,ye',iobjs(i)
   8           format(2i4,1p3e12.4,a12,i4)
            endif
         enddo
      enddo


      ifcurv=.false.

      do ic=1,4
         jc = mod1(ic+1,4)
         if (iobjs(ic).eq.iobjs(jc).and.iobjs(ic).ne.0) ifcurv=.true.
         write(6,*) 'IOBJ1:  ',ic,jc,iobjs(ic),iobjs(jc),ifcurv
      enddo
      if (.not.ifcurv) return


      do ic=1,4
         jc = mod1(ic+1,4)
         iobjc=iobjs(ic)
         write(6,*) 'iobjs:',ic,iobjc,iobjs(jc),' ',ccobjs(iobjc)

         if (iobjc.eq.iobjs(jc).and.iobjc.ne.0) then

            if (ccobjs(iobjc).eq.'o') then   ! Std. object
               xm=.5*(x(ic,e)+x(jc,e))
               ym=.5*(y(ic,e)+y(jc,e))
               zm=0.
               z1=0.
               call latchob(x1,y1,xm,ym,zm,dist2,ko,io,iobjc)
               rad=rad_circ(x(ic,e),x1,x(jc,e),y(ic,e),y1,y(jc,e))
               if (abs(rad).gt.0) then
c                 curve(1,ic,e)=rad         ! Fit circle
c                 ccurve(ic,e)='C'
                  curve(1,ic,e)=x1          ! Fit midside node
                  curve(2,ic,e)=y1
                  curve(3,ic,e)=z1
                  curve(5,ic,e)=iobjc       ! Save object number
                  ccurve(ic,e)='m'
               endif

            else                            ! Special object, like circle

               ccurve(ic,e)=ccobjs(iobjc)
               call copy(curve(1,ic,e),cobjs(1,iobjc),6)
               if (if3d) then
                  ic4 = ic+4
                  ccurve(ic4,e)=ccurve(ic,e)
                  call copy(curve(1,ic4,e),curve(1,ic,e),6)
               endif
            endif
         endif
         write(6,*)ccurve(ic,e),e,ic,(curve(k,ic,e),k=1,2),' chk_obj'
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine chk_neighbor(e) ! Check for shared curved sides
#     include "basics.inc"
      integer e
      logical ifcurv

      do ic = 1,4
         jc = mod1(ic+1,4)
         xm = 0.5*(x(ic,e)+x(jc,e))
         ym = 0.5*(y(ic,e)+y(jc,e))
         delta =  (x(ic,e)-x(jc,e))**2 + (y(ic,e)-y(jc,e))**2
         call getside(je,js,xm,ym)

         if (je.ne.e.and.je.ne.0) then
          js1 = mod1(js+1,4)
          xj = 0.5*(x(js,je)+x(js1,je))
          yj = 0.5*(y(js,je)+y(js1,je))
          delta2 = ((xm-xj)**2 + (ym-yj)**2)/delta
          if (delta2.lt.0.002) then
c              matched sides, give curve side attributes to new
c              element.... pff 3/27/94:  only "C" curve sides for now.
            if (ccurve(js,je).eq.'C') then
               ccurve(ic,e)='C'
               curve(1,ic,e) =  -curve(1,js,je)
               ifcurv = .true.
            elseif (ccurve(js,je).eq.'m') then
               ccurve(ic,e)='m'
               call copy(curve(1,ic,e),curve(1,js,je),3)
               ifcurv = .true.
            endif
          endif
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine clockw_chk(e,ierr) ! counter-clockwise check
#     include "basics.inc"
      integer e
      logical ifcurv
      real xyz(2,4)

      do i=1,4
         xyz(1,i)=x(i,e)
         xyz(2,i)=y(i,e)
      enddo

      c1=crss2d(xyz(1,2),xyz(1,4),xyz(1,1)) ! crss2d(a,b,o) := (a-o) x (b-o)
      c2=crss2d(xyz(1,3),xyz(1,1),xyz(1,2)) ! Note that the notation for corner 
      c3=crss2d(xyz(1,4),xyz(1,2),xyz(1,3)) ! numbers here differs slightly from 
      c4=crss2d(xyz(1,1),xyz(1,3),xyz(1,4)) ! that used in NEKTON the code.

      ierr=0
      if (c1.le.0.0.and.c2.le.0.0.and.
     $    c3.le.0.0.and.c4.le.0.0 ) then    ! cyclic permutation (cclock-wise):
          x(2,e) = xyz(1,4)                 ! reverse
          y(2,e) = xyz(2,4)
          x(4,e) = xyz(1,2)
          y(4,e) = xyz(2,2)
          if (if3d) then
             do i=1,4
                x(i+4,e)=x(i,e)
                y(i+4,e)=y(i,e)
             enddo
          endif
      elseif (c1.le.0.0.or.c2.le.0.0.or.   ! Some, but not all, < 0
     $        c3.le.0.0.or.c4.le.0.0 ) then
          ierr=1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine getside(je,js,xp,yp) ! Find closest element/side to duplicate
#     include "basics.inc"
      integer e,f

      rmin=1.e22

      je = 0

      nsides = 2*ndim
      if (ndim.eq.3) then
         do e=1,nel
            if (numapt(e).eq.ilevel) then ! only try those on same floor
               do f=1,nsides
                  r=    ((xp-sides(e,f,1))**2
     $            +      (yp-sides(e,f,2))**2)
                  if (r.lt.rmin) then
                     rmin=r
                     js = f
                     je = e
                  endif
               enddo
            endif
         enddo
      else                 ! 2D
         do e=1,nel
         do f=1,nsides
            r = (xp-sides(e,f,1))**2 + (yp-sides(e,f,2))**2
            if (r.lt.rmin) then
               rmin=r
               js = f
               je = e
            endif
c           write(6,1) e,f,r,rmin,(sides(e,f,k),k=1,2),xp,yp
c 1         format(i5,i2,1p6e12.4,' sides')
         enddo
         enddo
      endif
c     write(6,2) je,js,rmin,(sides(je,js,k),k=1,2),xp,yp
c 2   format(i5,i2,1p5e12.4,' sides2')

      return
      end
c-----------------------------------------------------------------------
      subroutine delelq(idel)
C     Deletes element from mesh.  But really, it copies a null element to
C     the one to be deleted.
#     include "basics.inc"
C
c     CALL BLANK(S,80)
c     WRITE(S,10,ERR=20) IDEL,NEL
c  10 FORMAT(' Deleting element',I5,' of',I5,'.$')
c     CALL PRS(S)
   20 CONTINUE
C
      if (idel.ne.nel) CALL COPYEL(NEL,IDEL)
      CALL COPYEL(NELM,NEL)
C     If the last element was a conduction element, reduce NCOND
      IF(MASKEL(NEL,1).EQ.0) NCOND=NCOND-1
      NEL=NEL-1
C     ILETAP=ILETAP-1
C
C     Recount the number of curved sides
C
      NCURVE=0
      DO 100 IE=1,NEL
      DO 100 IEDGE=1,12
         IF(CCURVE(IEDGE,IE).NE.' ') NCURVE=NCURVE+1
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine delel(idel)
C     Deletes element from mesh.  But really, it copies a null element to
C     the one to be deleted.
#     include "basics.inc"
C
      CALL BLANK(S,80)
      WRITE(S,10,ERR=20) IDEL,NEL
   10 FORMAT(' Deleting element',I5,' of',I5,'.$')
      CALL PRS(S)
   20 CONTINUE
C
      if (idel.ne.nel) CALL COPYEL(NEL,IDEL)
      CALL COPYEL(NELM,NEL)
C     If the last element was a conduction element, reduce NCOND
      IF(MASKEL(NEL,1).EQ.0) NCOND=NCOND-1
      NEL=NEL-1
C     ILETAP=ILETAP-1
C
C     Recount the number of curved sides
C
      NCURVE=0
      DO 100 IE=1,NEL
      DO 100 IEDGE=1,12
         IF(CCURVE(IEDGE,IE).NE.' ') NCURVE=NCURVE+1
  100 CONTINUE
C
      return
      end
c-----------------------------------------------------------------------
      subroutine copyel(isrc,ides) ! Copy everything related to element

      call copyel0   (isrc,ides)
      call copyelbc  (isrc,ides)
      call copyel2   (isrc,ides)

      return
      end
c-----------------------------------------------------------------------
      subroutine copyel_p(isrc,ides) ! Copy & preserve periodic bc

      call copyel0   (isrc,ides)
      call copyelbc_p(isrc,ides)
      call copyel2   (isrc,ides)

      return
      end
c-----------------------------------------------------------------------
      subroutine copyelbc_p(isrc,ides) ! Copies everything related to element
#     include "basics.inc"
      character key,string*6

      nshift = ides-isrc

      do 100 ifld=0,maxfld
      do 100 i=1,6
         cbc(i,ides,ifld) = cbc(i,isrc,ifld)
         call copy(bc(1,i,ides,ifld),bc(1,i,isrc,ifld),5)
         ibc(i,ides,ifld) = ibc(i,isrc,ifld)
         if (cbc(i,ides,ifld).eq.'P  ') 
     $       bc(1,i,ides,ifld) = bc(1,i,ides,ifld)+nshift
         if (cbc(i,ides,ifld).eq.'P  ') 
     $       ibc(i,ides,ifld) = ibc(i,ides,ifld)+nshift
  100 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine copyelbc(isrc,ides) ! Copies everything related to element
#     include "basics.inc"
      character key,string*6

      do 100 ifld=0,maxfld
      do 100 i=1,6
         cbc(i,ides,ifld) = cbc(i,isrc,ifld)
c        if (cbc(i,ides,ifld).eq.'P  ')call prsii('Reset P BC:$',i,ides)
c        if (cbc(i,ides,ifld).eq.'P  ') cbc(i,ides,ifld)='   '
         do j=1,5
            bc(j,i,ides,ifld) = bc(j,i,isrc,ifld)
         enddo
         ibc(i,ides,ifld) = ibc(i,isrc,ifld)
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine copyel0(isrc,ides) ! Copies everything related to element
#     include "basics.inc"
      character key,string*6

c     if (ides.ge.nelm-2) then
      if (ides.ge.nelm-1) then
        call prs
     $  ('***warning*** no additional elements allowed: arrays full$')
        write(6,*) 'ides,nelm:',ides,nelm
        return
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine copyel2(isrc,ides) ! Copies everything related to element
#     include "basics.inc"
      character key,string*6
      do 200 i =1,6
      do 200 ii=1,4
            sides(ides,i,ii)=sides(isrc,i,ii)
  200 continue
C
      DO 300 IC=1,8
         x(ic,ides)=x(ic,isrc)
         y(ic,ides)=y(ic,isrc)
         z(ic,ides)=z(ic,isrc)
  300 CONTINUE
      IGROUP(IDES)=IGROUP(ISRC)
C
C     Correct NCOND
C
      IF(MASKEL(IDES,1).EQ.0 .AND. MASKEL(ISRC,1).EQ.1)NCOND=NCOND-1
      IF(MASKEL(IDES,1).EQ.1 .AND. MASKEL(ISRC,1).EQ.0)NCOND=NCOND+1
      DO 500 IF=1,MAXFLD
         MASKEL(IDES,IF)=MASKEL(ISRC,IF)
  500 CONTINUE
      ISRT  (IDES) = ISRT (ISRC)
      ICRV  (IDES) = ICRV (ISRC)
      XCEN  (IDES) = XCEN (ISRC)
      YCEN  (IDES) = YCEN (ISRC)
      ZCEN  (IDES) = ZCEN (ISRC)
      RCEN  (IDES) = RCEN (ISRC)
      NUMAPT(IDES) =NUMAPT(ISRC)
      LETAPT(IDES) =LETAPT(ISRC)

      call copy  (x27(1,ides),x27(1,isrc),27)
      call copy  (y27(1,ides),y27(1,isrc),27)
      call copy  (z27(1,ides),z27(1,isrc),27)

      call chcopy(ccurve(1,ides),ccurve(1,isrc)  ,12)
      call copy  (curve(1,1,ides),curve(1,1,isrc),72)
c     call copy  (edges(1,idir,ides)=edges(iedge,idir,isrc)

      IF(NHIS.GT.0)THEN
C        Delete any integral sides
         IDEL=0
         DO 700 I=1,NHIS
            IF(HCODE(10,I).EQ.'I')IDEL=1
            IF(HCODE(10,I).EQ.'I')HCODE(10,I)=' '
  700    CONTINUE
         IF(IDEL.NE.0)THEN
            CALL PRS('***RESETTING INTEGRAL QUANTITY***$')
            CALL PRS(
     $   ' YOU MUST REDEFINE SIDES TO GET INTEGRAL FLUX, LIFT,OR DRAG$')
         ENDIF
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine model(iel)
C     Modifies element by moving point in mesh
C     Modified neighbor points iff they had been latched to point in questio
#     include "basics.inc"
      dimension ielmov(nelm)
      real xcheck(8),ycheck(8)
      character key,string*6
      logical iftmp
      integer e

      if (.not.if3d) then
         call model2d
         return
      endif

C
C     FIND CLOSEST ELEMENT AND CORNER
C     Erase whole isometric surface
c
      CALL PRS
     $('Push left button close to the corner you want changed.$')
      IFTMP =IFGRID
      IFGRID=.FALSE.
1     CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      IF(XSCR(XMOUSE).GT.1.0) THEN ! apparently trying to enter
         call prs('Type X-Y coordinates:$')
         call rerr(xmouse,ymouse)
      ELSE IF(BUTTON.EQ.'RIGHT')THEN
C        Latch to closest element vertex
         call blatch(xmouse,ymouse,icmn,iemn)
      ENDIF
      XM=XMOUSE
      YM=YMOUSE
      IFGRID=IFTMP
C     If the areas from the triangles formed by each side and the mouse point
C     are all positive (ccw) we must be inside the element!
      if (button.ne.'RIGHT') then
       IELMV=0
       DO 2 I=1,IEL
C        Skip apartments on other floors
         IF(.NOT. IF3D .OR. NUMAPT(I).EQ.ILEVEL)THEN
            IF(IFCEIL)THEN
               A1 = (x(i,5)-XM)*(y(i,6)-YM)-(x(i,6)-XM)*(y(i,5)-YM)
               A2 = (x(i,6)-XM)*(y(i,7)-YM)-(x(i,7)-XM)*(y(i,6)-YM)
               A3 = (x(i,7)-XM)*(y(i,8)-YM)-(x(i,8)-XM)*(y(i,7)-YM)
               A4 = (x(i,8)-XM)*(y(i,5)-YM)-(x(i,5)-XM)*(y(i,8)-YM)
            ELSE
               A1 = (x(i,1)-XM)*(y(i,2)-YM)-(x(i,2)-XM)*(y(i,1)-YM)
               A2 = (x(i,2)-XM)*(y(i,3)-YM)-(x(i,3)-XM)*(y(i,2)-YM)
               A3 = (x(i,3)-XM)*(y(i,4)-YM)-(x(i,4)-XM)*(y(i,3)-YM)
               A4 = (x(i,4)-XM)*(y(i,1)-YM)-(x(i,1)-XM)*(y(i,4)-YM)
            ENDIF
4           CONTINUE
C           We're inside the element
            IF(A1.GE.0.0.AND.A2.GE.0.0.AND.A3.GE.0.0.AND.A4.GE.0.0)
     $      IELMV=I
         ENDIF
2      CONTINUE
      else
       ielmv = iemn
      endif
      if (ielmv.eq.0) then
c        If it still = 0 then the point entered must have been outside all
c        elements.
         call prs
     $   ('Enter a point INSIDE one of the elements in order to$')
         call prs('move one of the element''s corners.  Try again.$')
         goto 1
      else
         call prsi('Moving element #: ',ielmv)
      endif
C
      IF(IFCEIL)THEN
         IC1=5
         IC2=8
      ELSE
         IC1=1
         IC2=4
      ENDIF
C     FIND CLOSEST CORNER
      rmin=1.0e8
      DO 3 I=IC1,IC2
         r=sqrt((x(i,ielmv)-xmouse)**2+(y(i,ielmv)-ymouse)**2)
         IF(R.LT.RMIN) THEN
            RMIN=R
            ICOMV=I
         ENDIF
3     CONTINUE
      xpicked=x(icomv,ielmv)
      ypicked=y(icomv,ielmv)

      CALL PRS
     $('Enter new point to which element corner is to be moved:$')
      CALL MOUSE(XMOVED,YMOVED,BUTTON)
      IF(XSCR(XMoved).GT.1.0) THEN
C        apparently trying to enter
         call prs('Type X-Y coordinates:$')
         call rerr(xmoved,ymoved)
      ELSEIF (BUTTON.EQ.'RIGHT') THEN
C        Latch to closest element vertex
         call blatch(xmoved,ymoved,icmn,iemn)
      ENDIF
C
C     HAVE IELMV,ICOMV  , NOW MOVE THE APPROPRIATE POINTS AND REDRAW ELEMENTS
C
C     First loop checks legality, second moves and erases, third redraws
      CALL DRAWIS(-IELMV)
      do 10 iloop=1,3
       do 10 iiel=1,iel
         if(iloop.eq.1) ielmov(iiel)=0
         DO 10 IICORN=IC1,IC2
C           Only move corners on same floor
            if(numapt(iiel).eq.numapt(ielmv)) then
C              Check floor for same global #'s  NO: Floors automatically
C              Modify overlapping corners; ceilings move 1 at a time.
C              Skip this check on 3rd loop; the redraw will know which to do
               IF(ILOOP.NE.3)THEN
C                    Only move corners close to the one picked
                     if (abs(x(iicorn,iiel)-xpicked).gt.xfac/100.
     $               .or.abs(y(iicorn,iiel)-ypicked).gt.yfac/100.)
     $               GO TO 9
C                     call prsii('moving$',iiel,iicorn)
               ENDIF
               IF(ILOOP.EQ.1)THEN
C
C                 First check if move is legal (has 4 acute <'s for corners
                  xtest=x(iicorn,iiel)
                  ytest=y(iicorn,iiel)
                  DO 8 I=IC1,IC2
                     xcheck(i)=x(i,iiel)
                     ycheck(i)=y(i,iiel)
8                 CONTINUE
                  XCHECK(IICORN)=XMOVED
                  YCHECK(IICORN)=YMOVED
C                 Check that the 4 cross products (areas) of each adjacent
C                 side pair is positive (angles between 0 and 180 degrees)
                  DO 109 ICC=1,4
                     I1=MOD(ICC-1,4)+1
                     I2=MOD(ICC  ,4)+1
                     I3=MOD(ICC+1,4)+1
                     IF (IFCEIL) THEN
                        I1=I1+4
                        I2=I2+4
                        I3=I3+4
                     ENDIF
                     X1=XCHECK(I2)-XCHECK(I1)
                     X2=XCHECK(I3)-XCHECK(I2)
                     Y1=YCHECK(I2)-YCHECK(I1)
                     Y2=YCHECK(I3)-YCHECK(I2)
                     AREA=X1*Y2-X2*Y1
                     IF(AREA.LE.0.0) THEN
                        CALL PRSIS(
     $                  '**ERROR** Element$',iiel,'Cannot be modified$')
                        CALL PRSIS(
     $                  'Angle at corner $',I2,' Would have been $')
                        CALL PRS('illegal (>180 degrees)$')
                        CALL PRS('No Elements were modified.$' )
                        DO 52 I=1,NEL
c                           CALL DRAWIS(ISRT (I))
52                      CONTINUE
                        return
                     ENDIF
109                 CONTINUE

               else if(iloop.eq.2)then 
                  call drawel(-iiel)  ! erase old element lines
                  call drawis(-iiel)
C                 Move corners
                  ielmov(iiel)=1
                  if (ifceil) then ! First move points on floor of upper level
                    do 18 e=1,iel
                     movede=0
                     do 16 ic=1,4
                       if (numapt(e).eq.ilevel+1) then
                         if  (abs(x(ic,e)-x(iicorn,iiel)).lt.xfac/100.
     $                   .and.abs(y(ic,e)-y(iicorn,iiel)).lt.yfac/100.)
     $                   then
c                           call prs('Moving upper$')
                            movede=1
                            x(ic,e)=xmoved
                            y(ic,e)=ymoved
                         endif
                       endif
16                    continue
                      if (movede.eq.1) call drawis(-iel)
18                  continue
                  elseif (.not.ifceil) then ! move pts on ceiling of lower level
                     if (ilevel.gt.1) then
                      do 14 e=1,iel
                       movede=0
                       do 12 ic=5,8
                         if(numapt(e).eq.ilevel-1)then
                           if(abs(x(ic,e)-x(iicorn,iiel)).lt.xfac/100.
     $                   .and.abs(y(ic,e)-y(iicorn,iiel)).lt.yfac/100.)
     $                     then
                              x(ic,e)=xmoved
                              y(ic,e)=ymoved
                              movede=1
                           endif
                         endif
12                      continue
                        if (movede.eq.1) call drawis(-iel)
14                     continue
                     endif
C                    Here we also move the points on the ceiling of the element,
C                    But only if there are no elements on any levels above
C                    the current one. (to avoid complications)
                     do 11 e=1,iel
                        if(numapt(e).gt.ilevel) then
                           CALL PRS('Presence of elements on higher '//
     $                     'floors inhibits modifying ceiling also$')
                           go to 15
                        endif
11                   continue
C                    Since we got here, there must be no elements on upper floo
                     x(iicorn+4,iiel)=xmoved
                     Y(IICORN+4,IIEL)=YMOVED
                  ENDIF
                  if(ifceil.and.iicorn.eq.5)then
C                    Manually erase old ceiling
                     CALL color(0)
                     call movec(x(8,iiel),y(8,iiel))
                     DO 114 IC=5,8
                        call drawc(x(ic,iiel),y(ic,iiel))
114                  CONTINUE
                     CALL color(10)
                  ENDIF
C                 Now move points on this level
15                x(iicorn,iiel)=xmoved
                  y(iicorn,iiel)=ymoved
c                  call prs('iiel,iicorn,xmoved$')
c                  call prs(iiel,iicorn,xmoved
C                 Move Center
                  if (.not.ifceil) then
                     xcen(iiel)=.25*vlsum(x(1,iiel),4)
                     ycen(iiel)=.25*vlsum(y(1,iiel),4)
                  endif
               ELSE IF(ILOOP.EQ.3)THEN
C                 Now, draw new elements for elements that were modified
                  IF(IICORN.NE.IC2)GO TO 9
                  if(ielmov(iiel).ne.0)call drawel(iiel)
               ENDIF
            endif
9           CONTINUE
10    CONTINUE
C     Draw ISOMETRICALLY ONLY ELEMENT MOVED
      CALL DRAWIS(IELMV)
      return
      end
c-----------------------------------------------------------------------
      subroutine blatch2(xm,ym,cm,em,ec,vc,nc)! close vertex
#     include "basics.inc"
      integer cm,em,ec(nel),vc(nel)
      integer e,v,v0,v1

      parameter (n3 = nxm*nym*nzm)
      common /ctmp2/ xp(n3),yp(n3),zp(n3),rrl(3)
      common /ctmp0/ erad(nelm)

      xc=xm
      yc=ym

      v0=1
      v1=4
      if (ifceil) v0=5
      if (ifceil) v1=8

      rm = 1.0e15

      do e=1,nel

         exmn = x(1,e)
         eymn = y(1,e)

         exmx = x(1,e)
         eymx = y(1,e)

         do v=v0,v1
            xt=x(v,e)
            yt=y(v,e)
            r=(xt-xc)**2+(yt-yc)**2
   
            if (r.lt.rm) then
               rm = r
               em = e
               cm = v
            endif

            exmn = min(x(v,e),exmn)
            eymn = min(y(v,e),eymn)

            exmx = max(x(v,e),exmx)
            eymx = max(y(v,e),eymx)

         enddo

         erad(e) = (exmx-exmn)+(eymx-eymn)
c        write(6,1) e,em,cm,xt,yt,xc,yc,rm
c  1     format(2i5,i3,1p5e12.4,' rmin2')

      enddo

      xm=x(cm,em)
      ym=y(cm,em)

      eps = 1.e-4
      nc  = 0    ! Count how many elements close to this vertex

      do e=1,nel
      do v=v0,v1
         xt=x(v,e)
         yt=y(v,e)
         r=(xt-xm)**2+(yt-ym)**2
         if (r.gt.0) r = sqrt(r)

         if (r.lt.eps*erad(e)) then
            nc = nc + 1
            ec(nc) = e
            vc(nc) = v

c           call genxyz_e (xp,yp,zp,e,nxm,nym,nzm)
c           call hilite   (e,xp,yp,zp,5)

c           call prsii('CONTINUE ?$',e,nc)
c           call res(ans,1)

         endif

c        write(6,*) v,e,nc,r,erad(e),' ERAD'

      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine model2d

C     Modifies element by moving point in mesh
C     Modified neighbor points iff they had been latched to point in questio

#     include "basics.inc"

      common /ctmps/ ec(nelm),vc(nelm)  ! Vertices close to moved point
      integer ec,vc,cm,em

      real xcheck(8),ycheck(8)
      character key,string*6
      logical iftmp
      integer e,v

      iftmp =ifgrid
      ifgrid=.false.

      call prs
     $('Push button close to the corner you want changed.$')

      call mouse  (xm,ym,button)
c     write(6,*) xm,ym,'  Corner clicked A'
      call blatch2(xm,ym,cm,em,ec,vc,nc)! close vertex
c     write(6,*) xm,ym,'  Corner latched, B'

      x0=xm
      y0=ym
      ifgrid=iftmp

      call prs
     $('Enter point to which vertex is to be moved:$')
      call mouse(x1,y1,button)
c     call diam2(x1,y1,4)

      if(xscr(x1).gt.1.0) then ! apparently trying to enter
         call prs('Type X-Y coordinates:$')
         call rerr(x1,y1)
      elseif (button.eq.'RIGHT') then
         call blatch(x1,y1,icmn,iemn) ! latch to closest vertex
      endif

c     call diam2(x1,y1,13)
c     write(6,*) x1,y1,'  TARGET GRID POINT'
c     call prs('CONTINUE ?$')
c     call res(ans,1)

      k    = nel+1
      kerr = 0

      do i=1,nc  ! Check only those elements that are attaced to vertex
         e = ec(i)
         v = vc(i)
         call copyel(e,k)   ! e --> k

c        call diam2(x(v,e),y(v,e),0)
c        write(6,*) v,e,x(v,e),y(v,e),x1,y1,'  ELEMENT POINT'
c        call prs('CONTINUE ?$')
c        call res(ans,1)

         dx = x1 - x(v,k)
         dy = y1 - y(v,k)
         dz = 0
         x(v,k) = x1
         y(v,k) = y1

         do jside=v-1,v   ! Move midside nodes
            j = jside
            if (j.lt.1) j=4

c           write(6,5) j,e,k,ccurve(j,k),(curve(jj,j,k),jj=1,5)

            if (ccurve(j,k).eq.'m') then
               xm = curve(1,j,k) + 0.5*dx
               ym = curve(2,j,k) + 0.5*dy
               zm = curve(3,j,k) + 0.5*dz

               if (curve(5,j,k).ne.0) then   ! Latch to existing obj
                  iobj=abs(curve(5,j,k))
                  call latchob(xo,yo,xm,ym,zm,dist2,ko,io,iobj)
                  xm = xo
                  ym = yo
               endif

               curve(1,j,k) = xm
               curve(2,j,k) = ym
               curve(3,j,k) = zm

            endif
         enddo

         call chkjac_e(ierr,k)
         if (ierr.eq.0) then
            call copyel(k,e)   ! k --> e
            call prsii('Moving element/vertex: $',e,v)
         else
            kerr = ierr
            call prsii('Error modifying element:$',e,ierr)
         endif
      enddo

      call gencen
      call redraw_mesh

      return
      end
c-----------------------------------------------------------------------
      subroutine chkjac_e(ierr,e)
      integer e

      ierr = 0

      return
      end
c-----------------------------------------------------------------------
      subroutine drised(iel,iedge,iflip)
C     DRaws ISometric EDge.  IFLIP CAuses to draw from end to beginning
C     POSITIVE MEANS CCW; ON VERTICAL STRUTS 9-12, MEANS UPWARD
C     Draws only.  No move or fill.
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
#     include "basics.inc"
C

      IF(IEDGE.GT.8)THEN
C        Vertical strut
         IC=IEDGE-4
         xisom(ic)=xphy(xscr(x(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         yisom(ic)=yphy(yscr(y(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         ic=iedge-8
         xisom(ic)=xphy(xscr(x(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         yisom(ic)=yphy(yscr(y(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         IF(IFLIP.EQ.1)THEN
            CALL MOVEC(XISOM(IEDGE-8),YISOM(IEDGE-8))
            CALL DRAWC(XISOM(IEDGE-4),YISOM(IEDGE-4))
         ELSE
            CALL MOVEC(XISOM(IEDGE-4),YISOM(IEDGE-4))
            CALL DRAWC(XISOM(IEDGE-8),YISOM(IEDGE-8))
         ENDIF
      ELSE IF(CCURVE(IEDGE,IEL).EQ.' ')THEN
C        Draw straight side
         IF(IFLIP.EQ.1)THEN
            IC=IEDGE+1
            IF(IC.EQ.5)IC=1
            IF(IC.EQ.9)IC=5
         ELSE
            IC=IEDGE
         ENDIF
         xisom(ic)=xphy(xscr(x(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         yisom(ic)=yphy(yscr(y(ic,iel))/5.0 + 0.8)-z(ic,iel)/20.
         call drawc(xisom(ic),yisom(ic))
      ELSE
C        Draw curved side
         NPOINT=10
         DO 18 I=1,NPOINT
            CSPACE(I)=(I-1.0)/(NPOINT-1.0)
18       CONTINUE
         CALL GETPTS(NPOINT,CSPACE,IEL,IEDGE,XCRVED,YCRVED)
         IF(IFLIP.EQ.1)THEN
            IBEGIN=1
            IEND  =NPOINT
         ELSE
            IBEGIN=NPOINT
            IEND  =1
         ENDIF
         DO 118 I=IBEGIN,IEND,IFLIP
            xi=xphy(xscr(xcrved(i))/5.0 + 0.8)-z(iedge,iel)/20.
            yi=yphy(yscr(ycrved(i))/5.0 + 0.8)-z(iedge,iel)/20.
            CALL DRAWC(XI,YI)
118      CONTINUE
      ENDIF
      return
      end
c-----------------------------------------------------------------------
      subroutine drawis(iel)
#     include "basics.inc"
      dimension xisom(8),yisom(8),cspace(100),xcrved(100),ycrved(100)
      iiel=iabs(iel)
      if (.not.if3d)   return
      if (nel.gt.1000) return
C        Now draw isometric view  (??! RESCALE??)
         IF(IEL.GT.0)call color(10)
         IF(IEL.le.0)call color(0)
         DO 7 IC=1,8
          xisom(ic)=xphy(xscr(x(ic,iiel))/5.0 + 0.8)-z(ic,iiel)/20.
          yisom(ic)=yphy(yscr(y(ic,iiel))/5.0 + 0.8)-z(ic,iiel)/20.
7        CONTINUE
         CALL BEGINB(XISOM(1),YISOM(1))
         DO 8 IEDGE=1,4
            call drised(iiel,iedge,1)
8        CONTINUE
         CALL ENDP
C        Now draw side panels
         IF(IEL.GT.0)CALL fillp(-14)
         IF(IEL.LE.0)CALL fillp(0)
         CALL BEGINB(XISOM(2),YISOM(2))
         call drised(iiel, 2, 1)
         call drised(iiel,11, 1)
         call drised(iiel, 6,-1)
         call drised(iiel,10,-1)
         call endp
C
         call beginb(xisom(4),yisom(4))
         call drised(iiel,12, 1)
         call drised(iiel, 7,-1)
         call drised(iiel,11,-1)
         call drised(iiel, 3, 1)
         call endp
C
C        Draw Ceiling panel
         IF(IEL.GT.0)CALL fillp(-15)
         IF(IEL.LE.0)CALL fillp(0)
         CALL BEGINB(XISOM(5),YISOM(5))
         DO 9 I=5,8
           call drised(iiel,i, 1)
9        CONTINUE
         CALL ENDP
      return
      end
c-----------------------------------------------------------------------
      subroutine drawel(iel)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
#     include "basics.inc"
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)

      if (.not.ifgraf) then
        if (iel.eq.(nelcap+1))
     $   call prsi('Showing only nelcap elements of$',nel)
        if (iel.gt.nelcap) return
      endif


C     Now, draw new elements for elements that were modified
      iiel=iabs(iel)
C

c     write (6,*) 'DRAWEL: iel = ',IEL

      IFSPHR=.FALSE.
      DO 7 IFACE=5,6
         IF (CCURVE(IFACE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
c
         xyzctr(1,1)=x(1,ie)+x(4,ie)+x(5,ie)+x(8,ie)
         xyzctr(1,2)=x(2,ie)+x(3,ie)+x(6,ie)+x(7,ie)
         xyzctr(1,3)=x(1,ie)+x(2,ie)+x(5,ie)+x(6,ie)
         xyzctr(1,4)=x(3,ie)+x(4,ie)+x(7,ie)+x(8,ie)
         xyzctr(1,5)=x(1,ie)+x(2,ie)+x(3,ie)+x(4,ie)
         xyzctr(1,6)=x(5,ie)+x(6,ie)+x(7,ie)+x(8,ie)
         xyzctr(2,1)=y(1,ie)+y(4,ie)+y(5,ie)+y(8,ie)
         xyzctr(2,2)=y(2,ie)+y(3,ie)+y(6,ie)+y(7,ie)
         xyzctr(2,3)=y(1,ie)+y(2,ie)+y(5,ie)+y(6,ie)
         xyzctr(2,4)=y(3,ie)+y(4,ie)+y(7,ie)+y(8,ie)
         xyzctr(2,5)=y(1,ie)+y(2,ie)+y(3,ie)+y(4,ie)
         xyzctr(2,6)=y(5,ie)+y(6,ie)+y(7,ie)+y(8,ie)
         xyzctr(3,1)=z(1,ie)+z(4,ie)+z(5,ie)+z(8,ie)
         xyzctr(3,2)=z(2,ie)+z(3,ie)+z(6,ie)+z(7,ie)
         xyzctr(3,3)=z(1,ie)+z(2,ie)+z(5,ie)+z(6,ie)
         xyzctr(3,4)=z(3,ie)+z(4,ie)+z(7,ie)+z(8,ie)
         xyzctr(3,5)=z(1,ie)+z(2,ie)+z(3,ie)+z(4,ie)
         xyzctr(3,6)=z(5,ie)+z(6,ie)+z(7,ie)+z(8,ie)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            if (j.eq.1) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.2) then
               call beginb(xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call endp
            endif
            if (j.eq.3) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call endp
            endif
            if (j.eq.4) then
               call beginb(xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.5) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call endp
            endif
            if (j.eq.6) then
               call beginb(xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
   72    continue
      else
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(-14)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
               call fillp(-15)
c              i15 = mod(iel,15)+1
c              call fillp(i15)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(-2)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp( -(IGROUP(IEL)+2) )
            ELSE
               CALL fillp( -(IGROUP(IEL)+3) )
            ENDIF
            IF(IFBWGKS)call fillp(0)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(0)
      ENDIF

      IC=4
      IF(IFCEIL)IC=8

      call beginb(x(ic,iiel),y(ic,iiel))

      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          ycenter=ycenter+y(ic,iiel)/4.
          xcenter=xcenter+x(ic,iiel)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          if(ccurve(iedge,iiel).eq.' ')then
             call drawc(x(ic,iiel),y(ic,iiel))
          ELSE
C            Draw curved side
             call getpts(npoint,cspace,iiel,iedge,xcrved,ycrved)
             DO 118 I=1,NPOINT
                CALL DRAWC(XCRVED(I),YCRVED(I))
118          CONTINUE
          ENDIF
6     CONTINUE
      CALL ENDP


      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
         WRITE(STRING,'(I3)')IEL
         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
C HMT color TRACE - seems to be the outline
      call color(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine mkside_e(e) ! Define sides' midpoints
#     include "basics.inc"
      integer e,f

      do f=1,nsides
         ic1=f
         ic2=f+1
         if (f.eq.4) ic2=1

         if (if3d) then ! This stuff only relevant for 3d

            ic3=ic1+4
            ic4=ic2+4
            if (f.eq.5) then
               ic1=1
               ic2=2
               ic3=3
               ic4=4
            elseif (f.eq.6) then
               ic1=1+4
               ic2=2+4
               ic3=3+4
               ic4=4+4
            endif
            xs =( x(ic1,iel)+x(ic2,iel)
     $         +  x(ic3,iel)+x(ic4,iel) )/4.
            ys =( y(ic1,iel)+y(ic2,iel)
     $         +  y(ic3,iel)+y(ic4,iel) )/4.
            zs =( z(ic1,iel)+z(ic2,iel)
     $         +  z(ic3,iel)+z(ic4,iel) )/4.
         else
            xs =( x(ic1,iel)+x(ic2,iel) )/2.
            ys =( y(ic1,iel)+y(ic2,iel) )/2.
            zs = 0.0
         endif

         sides (iel,f,1)=xs
         sides (iel,f,2)=ys
         sides (iel,f,3)=zs

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine mkside
#     include "basics.inc"
C
C     Find Sides' Midpoints
C
      DO 25 IEL=1,NEL
      DO 25 ISIDE=1,NSIDES
         IC1=ISIDE
         IC2=ISIDE+1
         IF(ISIDE.EQ.4)IC2=1
C        This stuff only relevant for 3d
         IC3=IC1+4
         IC4=IC2+4
         IF (ISIDE.EQ.5) THEN
            IC1=1
            IC2=2
            IC3=3
            IC4=4
         ELSEIF (ISIDE.EQ.6) THEN
            IC1=1+4
            IC2=2+4
            IC3=3+4
            IC4=4+4
         ENDIF
         if (if3d) then
            xs =( x(ic1,iel)+x(ic2,iel)
     $         +  x(ic3,iel)+x(ic4,iel) )/4.
            ys =( y(ic1,iel)+y(ic2,iel)
     $         +  y(ic3,iel)+y(ic4,iel) )/4.
            zs =( z(ic1,iel)+z(ic2,iel)
     $         +  z(ic3,iel)+z(ic4,iel) )/4.
         else
            xs =( x(ic1,iel)+x(ic2,iel) )/2.
            ys =( y(ic1,iel)+y(ic2,iel) )/2.
            zs = 0.0
         endif
         sides (iel,iside,1)=xs
         sides (iel,iside,2)=ys
         sides (iel,iside,3)=zs
25    continue

      return
      end
c-----------------------------------------------------------------------
      function crss2d(xy1,xy2,xy0)
      real xy1(2),xy2(2),xy0(2)

      v1x=xy1(1)-xy0(1)
      v2x=xy2(1)-xy0(1)
      v1y=xy1(2)-xy0(2)
      v2y=xy2(2)-xy0(2)
      crss2d = v1x*v2y - v1y*v2x

      return
      end
c-----------------------------------------------------------------------
      function round(x)

C     Try to Round X to fractional integer - eg .05 .1 .15 - if it's within 10-6

      integer icalld
      save    icalld
      data    icalld /0/

      round = x

      ax = abs(x)
      if (ax.lt.0.0499) return
   
      eps=1.0e-4
      eps=max(eps,ax*5.e-6)

      xtmp  = 20.0*ax + 0.5
      itmp  = xtmp
      round = .05*itmp
      if (x.lt.0) round = -round
      diff  = abs(round-x)

      if (diff.gt.eps) then
         round=x
      elseif (diff.gt.0) then
         icalld = icalld+1
c        write(6 ,*) icalld,round,x,diff,' round'
c        write(88,*) icalld,round,x,diff
         if (mod(icalld,50).eq.0)write(6,*) icalld,round,x,diff,' round'
      endif

      return
      end

c-----------------------------------------------------------------------
      subroutine blatch(xmouse,ymouse,icmin,ielmin) ! Latch to closest vertex
#     include "basics.inc"
      integer e

      xc=xmouse
      yc=ymouse
      rmin = 1.e15

      do e=1,nel    !  Find closest corner (NOT IN SAME ELEMENT)
      do iicorn=1,4
         icorn=iicorn
         if (ifceil) icorn=icorn+4
         xt=x(icorn,e)
         yt=y(icorn,e)
         r=(xt-xc)**2+(yt-yc)**2
         if (r.lt.rmin) then
            rmin  = r
            ielmin= e
            icmin = icorn
         endif
c        write(6,1) e,icorn,ielmin,icmin,xt,yt,xc,yc,rmin
c  1     format(i5,i3,i5,i3,1p5e12.4,' blatch')
      enddo
      enddo

      xmouse=x(icmin,ielmin)
      ymouse=y(icmin,ielmin)

      return
      end
c-----------------------------------------------------------------------
      subroutine sortel
#     include "basics.inc"
C     Sorts elements according to their visibility, i.e., ISRT (1) is behind
C     all the others and gets drawn first; element ISRT (NEL) is in front
C     and is most visible.
C     Elements with lowest zplane have lowqest sorted index.
C     Given the same zplane, elements with lowest xcen+ycen are lowest.
C
      INTEGER IND(NELM)
      EQUIVALENCE (IND,ZDEPTH(1,3))
C
      DO 1 IE=1,NEL
         ZDEPTH(IE,1) = 100.0*ZCEN(IE) + XCEN(IE)+YCEN(IE)
         ISRT(IE) = IE
1     CONTINUE
C
      CALL SORT(ZDEPTH,IND,NEL)
      CALL ISWAP(ISRT,ZDEPTH,IND,NEL)
C
      return
      end
c-----------------------------------------------------------------------
      subroutine fillpf(iCOLOR)
C     SELECTS FILL panel
      CHARACTER STRING*5
      include 'devices.inc'
c
      CALL PRS('Input color:$')
      CALL REI(ic)
      CALL XSETFILL(IC)
c
      return
      end
C
C----------------------------------------------------------------------
C
      subroutine drawel_COLOR(iel,ICOLOR)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
#     include "basics.inc"
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      iiel=iabs(iel)
C
      IFSPHR=.FALSE.
      DO 7 IFACE=5,6
         IF (CCURVE(IFACE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         xyzctr(1,1)=x(1,ie)+x(4,ie)+x(5,ie)+x(8,ie)
         xyzctr(1,2)=x(2,ie)+x(3,ie)+x(6,ie)+x(7,ie)
         xyzctr(1,3)=x(1,ie)+x(2,ie)+x(5,ie)+x(6,ie)
         xyzctr(1,4)=x(3,ie)+x(4,ie)+x(7,ie)+x(8,ie)
         xyzctr(1,5)=x(1,ie)+x(2,ie)+x(3,ie)+x(4,ie)
         xyzctr(1,6)=x(5,ie)+x(6,ie)+x(7,ie)+x(8,ie)
         xyzctr(2,1)=y(1,ie)+y(4,ie)+y(5,ie)+y(8,ie)
         xyzctr(2,2)=y(2,ie)+y(3,ie)+y(6,ie)+y(7,ie)
         xyzctr(2,3)=y(1,ie)+y(2,ie)+y(5,ie)+y(6,ie)
         xyzctr(2,4)=y(3,ie)+y(4,ie)+y(7,ie)+y(8,ie)
         xyzctr(2,5)=y(1,ie)+y(2,ie)+y(3,ie)+y(4,ie)
         xyzctr(2,6)=y(5,ie)+y(6,ie)+y(7,ie)+y(8,ie)
         xyzctr(3,1)=z(1,ie)+z(4,ie)+z(5,ie)+z(8,ie)
         xyzctr(3,2)=z(2,ie)+z(3,ie)+z(6,ie)+z(7,ie)
         xyzctr(3,3)=z(1,ie)+z(2,ie)+z(5,ie)+z(6,ie)
         xyzctr(3,4)=z(3,ie)+z(4,ie)+z(7,ie)+z(8,ie)
         xyzctr(3,5)=z(1,ie)+z(2,ie)+z(3,ie)+z(4,ie)
         xyzctr(3,6)=z(5,ie)+z(6,ie)+z(7,ie)+z(8,ie)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            if (j.eq.1) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.2) then
               call beginb(xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call endp
            endif
            if (j.eq.3) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call endp
            endif
            if (j.eq.4) then
               call beginb(xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.5) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call endp
            endif
            if (j.eq.6) then
               call beginb(xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
   72    continue

      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
C              call fillp(icolor)
               i15 = mod(iel,15)+1
               call fillp(icolor)
C               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8

      call beginb(x(ic,iiel),y(ic,iiel))

      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          ycenter=ycenter+y(ic,iiel)/4.
          xcenter=xcenter+x(ic,iiel)/4.
          IEDGE=IC-1
          IF(IC.EQ.1)IEDGE=4
          IF(IC.EQ.5)IEDGE=8
          if(ccurve(iedge,iiel).eq.' ')then
             call drawc(x(ic,iiel),y(ic,iiel))
          ELSE
C            Draw curved side
             call getpts(npoint,cspace,iiel,iedge,xcrved,ycrved)
             do 118 i=1,npoint
                call drawc(xcrved(i),ycrved(i))
118          continue
          endif
6     continue
      call endp

      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
         WRITE(STRING,'(I3)')IEL
         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
      call color(1)
      return
      end
C
C----------------------------------------------------------------------
C
      subroutine hmt_drawel_COLOR(iel,ICOLOR)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
#     include "basics.inc"
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      iiel=iabs(iel)
C
      IFSPHR=.FALSE.
      DO 7 IFACE=5,6
         IF (CCURVE(IFACE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         xyzctr(1,1)=x(1,ie)+x(4,ie)+x(5,ie)+x(8,ie)
         xyzctr(1,2)=x(2,ie)+x(3,ie)+x(6,ie)+x(7,ie)
         xyzctr(1,3)=x(1,ie)+x(2,ie)+x(5,ie)+x(6,ie)
         xyzctr(1,4)=x(3,ie)+x(4,ie)+x(7,ie)+x(8,ie)
         xyzctr(1,5)=x(1,ie)+x(2,ie)+x(3,ie)+x(4,ie)
         xyzctr(1,6)=x(5,ie)+x(6,ie)+x(7,ie)+x(8,ie)
         xyzctr(2,1)=y(1,ie)+y(4,ie)+y(5,ie)+y(8,ie)
         xyzctr(2,2)=y(2,ie)+y(3,ie)+y(6,ie)+y(7,ie)
         xyzctr(2,3)=y(1,ie)+y(2,ie)+y(5,ie)+y(6,ie)
         xyzctr(2,4)=y(3,ie)+y(4,ie)+y(7,ie)+y(8,ie)
         xyzctr(2,5)=y(1,ie)+y(2,ie)+y(3,ie)+y(4,ie)
         xyzctr(2,6)=y(5,ie)+y(6,ie)+y(7,ie)+y(8,ie)
         xyzctr(3,1)=z(1,ie)+z(4,ie)+z(5,ie)+z(8,ie)
         xyzctr(3,2)=z(2,ie)+z(3,ie)+z(6,ie)+z(7,ie)
         xyzctr(3,3)=z(1,ie)+z(2,ie)+z(5,ie)+z(6,ie)
         xyzctr(3,4)=z(3,ie)+z(4,ie)+z(7,ie)+z(8,ie)
         xyzctr(3,5)=z(1,ie)+z(2,ie)+z(3,ie)+z(4,ie)
         xyzctr(3,6)=z(5,ie)+z(6,ie)+z(7,ie)+z(8,ie)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            if (j.eq.1) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.2) then
               call beginb(xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call endp
            endif
            if (j.eq.3) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call endp
            endif
            if (j.eq.4) then
               call beginb(xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.5) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call endp
            endif
            if (j.eq.6) then
               call beginb(xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
   72    continue
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
            call color(10)
            IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
C              call fillp(icolor)
               i15 = mod(iel,15)+1
C HMT               call fillp(icolor)
               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      IC=4
      IF(IFCEIL)IC=8

      call beginb(X(ic,iiel),Y(ic,iiel))

      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      DO 6 IC=ICBEG,ICEND
          ccenter=xcenter+x(ic,iiel)/4.
          ycenter=ycenter+y(ic,iiel)/4.
          iedge=ic-1
          if(ic.eq.1)iedge=4
          if(ic.eq.5)iedge=8
          if(ccurve(iedge,iiel).eq.' ')then
             call drawc(x(ic,iiel),y(ic,iiel))
          else
C            Draw curved side
             call getpts(npoint,cspace,iiel,iedge,xcrved,ycrved)
             do 118 i=1,npoint
                call drawc(xcrved(i),ycrved(i))
118          continue
          endif
6     CONTINUE
      CALL ENDP
C
C HMT
C      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
C HMT
C         WRITE(STRING,'(I3)')IEL
C         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
C         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
C         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
C      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
      call color(1)
      return
      end
C
C----------------------------------------------------------------------
C
      subroutine hmt_drawel_COLOR_NLN(iel,ICOLOR)
C     IF ELEMENT NUMBER IS NEGATIVE, ERASE ELEMENT
      integer icolor
#     include "basics.inc"
      CHARACTER STRING*6
      DIMENSION XISOM(8),YISOM(8),CSPACE(100),XCRVED(100),YCRVED(100)
      LOGICAL IFSPHR
      DIMENSION ZBUFF(6),XYZCTR(3,6)
      DIMENSION IND(6)
C     Now, draw new elements for elements that were modified
      iiel=iabs(iel)
C
      IFSPHR=.FALSE.
      DO 7 IEDGE=5,6
         IF (CCURVE(IEDGE,IEL).EQ.'s') IFSPHR=.TRUE.
    7 CONTINUE
      IF (IFSPHR) THEN
         IE  =IABS(IEL)
C
C        Draw spherical mesh
C
         xyzctr(1,1)=x(1,ie)+x(4,ie)+x(5,ie)+x(8,ie)
         xyzctr(1,2)=x(2,ie)+x(3,ie)+x(6,ie)+x(7,ie)
         xyzctr(1,3)=x(1,ie)+x(2,ie)+x(5,ie)+x(6,ie)
         xyzctr(1,4)=x(3,ie)+x(4,ie)+x(7,ie)+x(8,ie)
         xyzctr(1,5)=x(1,ie)+x(2,ie)+x(3,ie)+x(4,ie)
         xyzctr(1,6)=x(5,ie)+x(6,ie)+x(7,ie)+x(8,ie)
         xyzctr(2,1)=y(1,ie)+y(4,ie)+y(5,ie)+y(8,ie)
         xyzctr(2,2)=y(2,ie)+y(3,ie)+y(6,ie)+y(7,ie)
         xyzctr(2,3)=y(1,ie)+y(2,ie)+y(5,ie)+y(6,ie)
         xyzctr(2,4)=y(3,ie)+y(4,ie)+y(7,ie)+y(8,ie)
         xyzctr(2,5)=y(1,ie)+y(2,ie)+y(3,ie)+y(4,ie)
         xyzctr(2,6)=y(5,ie)+y(6,ie)+y(7,ie)+y(8,ie)
         xyzctr(3,1)=z(1,ie)+z(4,ie)+z(5,ie)+z(8,ie)
         xyzctr(3,2)=z(2,ie)+z(3,ie)+z(6,ie)+z(7,ie)
         xyzctr(3,3)=z(1,ie)+z(2,ie)+z(5,ie)+z(6,ie)
         xyzctr(3,4)=z(3,ie)+z(4,ie)+z(7,ie)+z(8,ie)
         xyzctr(3,5)=z(1,ie)+z(2,ie)+z(3,ie)+z(4,ie)
         xyzctr(3,6)=z(5,ie)+z(6,ie)+z(7,ie)+z(8,ie)
         TMP=.25
         CALL CMULT(XYZCTR,TMP,18)
         DO 71 I=1,6
            ZBUFF(I)=ZISO(XYZCTR(1,I),XYZCTR(2,I),XYZCTR(3,I))
   71    CONTINUE
         CALL SORT(ZBUFF,IND,6)
         DO 72 I=1,6
            J=IND(I)
            if (j.eq.1) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.2) then
               call beginb(xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call endp
            endif
            if (j.eq.3) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call endp
            endif
            if (j.eq.4) then
               call beginb(xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
            if (j.eq.5) then
               call beginb(xiso(x(1,ie),y(1,ie),z(1,ie)) 
     $                    ,yiso(x(1,ie),y(1,ie),z(1,ie)))
               call drawc (xiso(x(2,ie),y(2,ie),z(2,ie)) 
     $                    ,yiso(x(2,ie),y(2,ie),z(2,ie)))
               call drawc (xiso(x(3,ie),y(3,ie),z(3,ie)) 
     $                    ,yiso(x(3,ie),y(3,ie),z(3,ie)))
               call drawc (xiso(x(4,ie),y(4,ie),z(4,ie)) 
     $                    ,yiso(x(4,ie),y(4,ie),z(4,ie)))
               call endp
            endif
            if (j.eq.6) then
               call beginb(xiso(x(5,ie),y(5,ie),z(5,ie)) 
     $                    ,yiso(x(5,ie),y(5,ie),z(5,ie)))
               call drawc (xiso(x(6,ie),y(6,ie),z(6,ie)) 
     $                    ,yiso(x(6,ie),y(6,ie),z(6,ie)))
               call drawc (xiso(x(7,ie),y(7,ie),z(7,ie)) 
     $                    ,yiso(x(7,ie),y(7,ie),z(7,ie)))
               call drawc (xiso(x(8,ie),y(8,ie),z(8,ie)) 
     $                    ,yiso(x(8,ie),y(8,ie),z(8,ie)))
               call endp
            endif
   72    continue
      ELSE
C
C     Std draw:
C
      npoint=10
      DO 8 I=1,NPOINT
         CSPACE(I)=(I-1.0)/(NPOINT-1.0)
8     CONTINUE
      IF(IEL.GT.0)THEN
C        Normal Draw, unless it is a conduction element
C           Draw conduction elements red
C HMT Looking for outline
C         call color(icolor)
         call hmt_color(icolor)
         IF(IGROUP(IEL).LT.0)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.0)THEN
C              NORMAL
c              call fillp(icolor)
               i15 = mod(iel,15)+1
C HMT HACK - found where he fills elms
C               write (6,*) 'icolor = ',icolor
C               write (6,*) 'i15 = ',i15
               call hmt_fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.1)THEN
               call fillp(icolor)
            ELSE IF(IGROUP(IEL).EQ.2)THEN
               CALL fillp(icolor)
            ELSE
               CALL fillp(icolor)
            ENDIF
            IF(IFBWGKS)call fillp(icolor)
      ELSE
C        fill black (i.e., erase)
         CALL color(0)
         CALL fillp(icolor)
      ENDIF
C     One more Kludge: if IEL is .GT. 10,000 then draw outline only
      ic=4
      if (ifceil) ic=8
      call beginb(x(ic,iiel),y(ic,iiel))

      XCENTER=0.0
      YCENTER=0.0
      IF(IFCEIL)THEN
         ICBEG=5
         ICEND=8
      ELSE
         ICBEG=1
         ICEND=4
      ENDIF
      do 6 ic=icbeg,icend
          ycenter=ycenter+y(ic,iiel)/4.
          xcenter=xcenter+x(ic,iiel)/4.
          iedge=ic-1
          if(ic.eq.1)iedge=4
          if(ic.eq.5)iedge=8
          if(ccurve(iedge,iiel).eq.' ')then
             call drawc(x(ic,iiel),y(ic,iiel))
          else
C            Draw curved side
             call getpts(npoint,cspace,iiel,iedge,xcrved,ycrved)
             do 118 i=1,npoint
                call drawc(xcrved(i),ycrved(i))
118          continue
          endif
6     continue
      CALL ENDP
C
      IF(IEL.GT.0)THEN
C        LABEL Element Center
c        IF(IF3D)     WRITE(STRING,'(I3,A1)')NUMAPT(IEL),LETAPT(IEL)
c        IF(.NOT.IF3D)WRITE(STRING,'(I3)')IEL
C HMT Turn numbers off
C         WRITE(STRING,'(I3)')IEL
C         IF(IFCEIL)WRITE(STRING(5:5),'(A1)')'C'
C         STRING(6:6)='$'
c        IF(CWRITE.NE.0.)CALL GWRITE(XCENTER-.055*XFAC,
c    $   YCENTER-.02*YFAC,1.0,STRING)
C         IF(CWRITE.NE.0.)CALL GWRITE(XCENTER,YCENTER,1.0,STRING)
      ENDIF
C
C     End of spherical - regular choice
C
      ENDIF
C      call color(1)
      return
      end
c-----------------------------------------------------------------------
      subroutine flipel(ieg,fplane)
c
c     This routine flips about x, y, or z plane
c
#     include "basics.inc"
      character key,string*6
      character*1 fplane
      character*3 cbt
c
      integer efac(6)
      save    efac
      data    efac  / 4,2,1,3,5,6 /
c
      logical ifswap(8)
c
      i1 = efac(3)
      i2 = efac(4)
      if (if3d) then
         i1 = efac(5)
         i2 = efac(6)
      endif
c
      DO IFLD=0,MAXFLD
         cbt              = CBC(i1,ieg,IFLD)
         CBC(i1,ieg,IFLD) = CBC(i2,ieg,IFLD)
         CBC(i2,ieg,IFLD) = cbt
         if (CBC(i1,ieg,IFLD).eq.'P  ') CBC(i1,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'P  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'E  ') CBC(i1,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'E  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'SP ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i1,ieg,IFLD).eq.'J  ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'SP ') CBC(i2,ieg,IFLD) = '   '
         if (CBC(i2,ieg,IFLD).eq.'J  ') CBC(i2,ieg,IFLD) = '   '
         do j=1,5
            Bt                = BC(J,i1,ieg,IFLD)
            BC(J,i1,ieg,IFLD) = BC(J,i2,ieg,IFLD)
            BC(J,i2,ieg,IFLD) = bt
         enddo
      enddo
C
      DO II=1,4
         st              =sides(ieg,i1,II)
         sides(ieg,i1,II)=sides(ieg,i2,II)
         sides(ieg,i2,II)=st
      enddo
C
      if (if3d) then
         do ic=1,4
            i4 = ic+4
            xt       =x(ic,ieg)
            x(ic,ieg)=x(i4,ieg)
            x(i4,ieg)=xt
            yt       =y(ic,ieg)
            y(ic,ieg)=y(i4,ieg)
            y(i4,ieg)=yt
            zt       =z(ic,ieg)
            z(ic,ieg)=z(i4,ieg)
            z(i4,ieg)=zt
         enddo
      else
         xt      =x(1,ieg)
         x(1,ieg)=x(4,ieg)
         x(4,ieg)=xt
         xt      =x(2,ieg)
         x(2,ieg)=x(3,ieg)
         x(3,ieg)=xt
         yt      =y(1,ieg)
         y(1,ieg)=y(4,ieg)
         y(4,ieg)=yt
         yt      =y(2,ieg)
         y(2,ieg)=y(3,ieg)
         y(3,ieg)=yt
      endif
c
      if (fplane.eq.'x') then
         idir = 1
         xcen(ieg) = -xcen(ieg)
         do ic=1,8
            x(ic,ieg) = -x(ic,ieg)
         enddo
      elseif (fplane.eq.'y') then
         idir = 2
         ycen(ieg) = -ycen(ieg)
         do ic=1,8
            y(ic,ieg) = -y(ic,ieg)
         enddo
      else
         idir = 3
         zcen(ieg) = -zcen(ieg)
         do ic=1,8
            z(ic,ieg) = -z(ic,ieg)
         enddo
      endif
C
C     Correct NCOND
C
      IF(MASKEL(ieg,1).EQ.0 .AND. MASKEL(ieg,1).EQ.1)NCOND=NCOND-1
      IF(MASKEL(ieg,1).EQ.1 .AND. MASKEL(ieg,1).EQ.0)NCOND=NCOND+1
      DO 500 IF=1,MAXFLD
         MASKEL(ieg,IF)=MASKEL(ieg,IF)
  500 CONTINUE
c     ISRT  (ieg) = ISRT (ieg)
c     ICRV  (ieg) = ICRV (ieg)
c     NUMAPT(ieg) =NUMAPT(ieg)
c     LETAPT(ieg) =LETAPT(ieg)
      if (if3d) then
         do iface=1,6
            if (ccurve(iface,ieg).eq.'s') 
     $      curve(idir,iface,ieg)=-curve(idir,iface,ieg)
         enddo
         if (ccurve(5,ieg).eq.'s'.or.ccurve(6,ieg).eq.'s') then
            cbt           = ccurve(6,ieg)
            ccurve(6,ieg) = ccurve(5,ieg)
            ccurve(5,ieg) = cbt
            do i=1,5
               ct             = curve(i,6,ieg)
               curve(i,6,ieg) = curve(i,5,ieg)
               curve(i,5,ieg) = ct
            enddo
         endif
c
         do iedge=1,8
            ifswap(iedge)=.false.
         enddo
c
         do iedge=1,4
            iedg4=iedge+4
            if (ccurve(iedge,ieg).eq.'C'.and. .not.ifswap(iedge)) then
c
               ifswap(iedge) = .true.
               ifswap(iedg4) = .true.
c
               cbt           = ccurve(6,ieg)
               ccurve(6,ieg) = ccurve(5,ieg)
               ccurve(5,ieg) = cbt
               do i=1,5
                  ct                 = curve(i,iedg4,ieg)
                  curve(i,iedg4,ieg) = curve(i,iedge,ieg)
                  curve(i,iedge,ieg) = ct
               enddo
            endif
         enddo
c
         do iedge=5,8
            iedg4=iedge-4
            if (ccurve(iedge,ieg).eq.'C'.and. .not.ifswap(iedge)) then
c
               ifswap(iedge) = .true.
               ifswap(iedg4) = .true.
c
               cbt           = ccurve(6,ieg)
               ccurve(6,ieg) = ccurve(5,ieg)
               ccurve(5,ieg) = cbt
               do i=1,5
                  ct                 = curve(i,iedg4,ieg)
                  curve(i,iedg4,ieg) = curve(i,iedge,ieg)
                  curve(i,iedge,ieg) = ct
               enddo
            endif
         enddo
      else   ! 2D:  swap 1 & 3
         cbt           = ccurve(1,ieg)
         ccurve(1,ieg) = ccurve(3,ieg)
         ccurve(3,ieg) = cbt
c
         do i=1,5
            ct             = curve(i,1,ieg)
            curve(i,1,ieg) = curve(i,3,ieg)
            curve(i,3,ieg) = ct
         enddo
      endif

      if (fplane.eq.'x') i=1
      if (fplane.eq.'y') i=2
      if (fplane.eq.'z') i=3
      nedge = 4 + 8*(ndim-2)
      do k=1,nedge
         if (ccurve(k,ieg).eq.'m') curve(i,k,ieg) = -curve(i,k,ieg)
      enddo


      return
      end
c-----------------------------------------------------------------------
      subroutine reflect_mesh
c
c     Reflect about x-, y-, or z-plane
c
#     include "basics.inc"
      character*1 fplane
c
1     CONTINUE
      ITEM(1)='BUILD MENU'
      ITEM(2)='REFLECT ABOUT X=0'
      ITEM(3)='REFLECT ABOUT Y=0'
      ITEM(4)='REFLECT ABOUT Z=0'
      NCHOIC=3
      if (if3d) NCHOIC=4
c
      CALL MENU(XMOUSE,YMOUSE,BUTTON,'Reflect Mesh')
c
3013  CONTINUE
      IF (CHOICE.EQ.'BUILD MENU') return
      IF (CHOICE.EQ.'REFLECT ABOUT X=0') fplane = 'x'
      IF (CHOICE.EQ.'REFLECT ABOUT Y=0') fplane = 'y'
      IF (CHOICE.EQ.'REFLECT ABOUT Z=0') fplane = 'z'
c
      do ieg=1,nel
c        write(6,*) fplane,'  flip: ',ieg
         call flipel(ieg,fplane)
      enddo
c
      goto 1
c
      return
      end
c-----------------------------------------------------------------------
      subroutine substitute_el(xyzbox,nelold)
C     Delete elements in xyzbox, provied e < nelold + 1
#     include "basics.inc"
      real xyzbox(6)
c
      common /cdell/ ifkeep(nelm)
      logical ifkeep
      integer e
C
      do e=1,nel
         ifkeep(e) = .true.
      enddo
c
      do e=1,nelold
         if (if3d) then
            if (xyzbox(1).le.xcen(e).and.xcen(e).le.xyzbox(2) .and.
     $          xyzbox(3).le.ycen(e).and.ycen(e).le.xyzbox(4) .and.
     $          xyzbox(5).le.zcen(e).and.zcen(e).le.xyzbox(6)) then
                ifkeep(e)=.false.
            endif
         else
            if (xyzbox(1).le.xcen(e).and.xcen(e).le.xyzbox(2)  .and.
     $          xyzbox(3).le.ycen(e).and.ycen(e).le.xyzbox(4)) then
                ifkeep(e)=.false.
            endif
         endif
c        if (ifkeep(e)) write(6,*) 'KEEP:',e
         if (.not.ifkeep(e)) write(6,*) 'NO KEEP:',e,if3d
     $                                   ,xcen(e),ycen(e),zcen(e)
      enddo
c
c     Now, shift all elements down, to compress out deleted elements
c
      jskip = 0
      je    = 0
      do e=1,nel
         if (ifkeep(e)) then
            je = je+1
            if (je.lt.e) then
               call copyel(e,je)
C              If the last element was a conduction element, reduce NCOND
               if(maskel(je,1).eq.0) ncond=ncond-1
            endif
         else
            write(s,10) e,nel
   10       format(' Deleting element',I9,' of',I9,'.$')
            call prs(s)
         endif
      enddo
      nel = je
C
C     Recount the number of curved sides
C
      ncurve=0
      do 100 ie=1,nel
      do 100 iedge=1,12
         if(ccurve(iedge,ie).ne.' ') ncurve=ncurve+1
  100 continue
C
      return
      end
c-----------------------------------------------------------------------
      function rad_circ(x0i,x1i,x2i,y0i,y1i,y2i)
      real*8 x0,x1,x2,y0,y1,y2

      real*8 rad_circ8,max_ratio
      real*8 dx1,dx2,dx3,dy1,dy2,dy3,aaa,bbb,area

      rad_circ = 0

      x0 = x0i ! real*4 --> real*8 conversion
      x1 = x1i
      x2 = x2i
      y0 = y0i
      y1 = y1i
      y2 = y2i

c     Directed arc:

      dx1  = x1-x0
      dx2  = x2-x0
      dx3  = x2-x1
      dy1  = y1-y0
      dy2  = y2-y0
      dy3  = y2-y1
      aaa  = dx1*dx1 + dx2*dx2 + dx3*dx3 + dy1*dy1 + dy2*dy2 + dy3*dy3
      bbb  = (dx1*dx1+dy1*dy1) * (dx2*dx2+dy2*dy2) * (dx3*dx3+dy3*dy3)
      area = 0.5*(dx1*dy2 - dx2*dy1)

      tol  = 1.e-5
      tol  = 1.e-9
c     write(6,*) x0,y0,' x0'
c     write(6,*) x1,y1,' x1'
c     write(6,*) x2,y2,' x2'
      if (abs(area).lt.tol*aaa) return  ! nearly colinear
      rad_circ8 = 0.
      if (bbb.gt.0) rad_circ8 = 0.25*sqrt(bbb)/area

c     x20 = dx2**2 + dy2**2
c     if (x20.gt.0) x20=sqrt(x20)
c     write(6 ,1) 'area:',area,aaa,bbb,rad_circ8,x20
c     write(86,1) 'area:',area,aaa,bbb,rad_circ8,x20
c   1 format(a5,1p5e12.4)


c     if (x20.gt.0) x20=sqrt(x20)
c     max_ratio = 20.
c     if (abs(rad_circ).gt.max_ratio*x20) rad_circ=0

      rad_circ = rad_circ8

      return
      end
c-----------------------------------------------------------------------
      subroutine side_to_obj(iobj,xo,yo,e,side)
#     include "basics.inc"

c     Is this side attached to object iobj ?

      integer e,side,v


      tolobj = 1.e-4*(griddx**2+griddy**2)

      do iobj=1,nobjs
         if (ifobjg(iobj)) then ! Active object
            nv = 0              ! Number of vertices on this object
            do j=side,side+1
               v = j
               if (v.gt.4) v=1
               call latchob(xo,yo,x(v,e),y(v,e),0.,dist2,ko,io,iobj)
               if (dist2.lt.tolobj) nv = nv+1
               if (nv.eq.2) return
            enddo
         endif
      enddo

      iobj = 0   ! Nothing found

      return
      end
c-----------------------------------------------------------------------
      subroutine convert_m_to_c_e(f,e,rad_mx)
#     include "basics.inc"
      integer f,e

      if (ccurve(f,e).ne.'m') return


      ic=f
      jc=ic+1
      if (ic.eq.4) jc=1
      if (ic.eq.8) jc=5

      xm = curve(1,f,e)
      ym = curve(2,f,e)
      zm = curve(3,f,e)
      write(6,*) f,e,' ',ccurve(f,e),(curve(k,f,e),k=1,2),' curve'

      rad = rad_circ(x(ic,e),xm,x(jc,e),y(ic,e),ym,y(jc,e))

      ccurve(f,e)=' '
      call rzero(curve(1,f,e),6)

      if (rad.eq.0.0 .or. abs(rad).gt.rad_mx) return

      ccurve (f,e)='C'
      curve  (1,f,e)=rad

      return
      end
c-----------------------------------------------------------------------
      subroutine convert_m_to_c_all ! This works only for 2D at present

#     include "basics.inc"
      integer edge,e

      call prs('Input maxium radius:$')
      call rer(rad_mx)
      
      do e=1,nel
      do edge=1,4*(ndim-1)
         call convert_m_to_c_e(edge,e,rad_mx)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine convert_m_to_c_allr(rmx) ! This works only for 2D at present

#     include "basics.inc"
      integer edge,e

      do e=1,nel
      do edge=1,4*(ndim-1)
         call convert_m_to_c_e(edge,e,rmx)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
