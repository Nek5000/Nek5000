c-----------------------------------------------------------------------
      subroutine octspl
      include 'basics.inc'
      COMMON /SPLITT/ ENEW(NELM),IND(NELM)
      DIMENSION XPT(3),VEC(3),LISTE(8)
c
      logical ifoct,ifhex,ifzsp

      call set_zgml (nxm)

C     Check on the number and numbering of the new elements to be generated.

      Noct = 4
      ifoct = .true.
      ifhex = .false.
      ifzsp = .false.

      if (if3d) then
         call prs('Oct, quad, hex, z, multi-split? (o/q/h/z/m)$')
         call res(ANS,1)
         ifoct = .false.
         ifhex = .false.
         ifzsp = .false.
         if (ANS.eq.'o'.or.ANS.eq.'O') then
            Noct = 8
            ifoct = .true.
         elseif (ans.eq.'h'.or.ans.eq.'H') then
            Noct = 6
            ifhex = .true.
         elseif (ans.eq.'z'.or.ans.eq.'Z') then
            Noct = 2
            ifzsp = .true.
         elseif (ans.eq.'m'.or.ans.eq.'M') then
            call msplit
            return
         endif
      else
         call prs('quad or multi-split? (q/m)$')
         call res(ANS,1)
         if (ans.eq.'m'.or.ans.eq.'M') then
            call msplit
            return
         endif
      endif
c
      neln=nel*noct
      nelm1=nelm-2
      if (neln.gt.nelm1) then
       n1   =nelm-2
       call prs('WARNING: Number of elements after OctSplit operation$')
       call prsii('would be greater than the allowed maximum:$',neln,n1)
       call prs('Aborting zipper.$')
       return
      else
       call prsi('Number of elements after OctSplit operation:$',neln)
       call prs ('Are you sure you want to split (y/n)?$')
       call res(ans,1)
       if (ans.ne.'Y'.and.ans.ne.'y') return
      endif

c     Renumber ALL elements, so that low numbered elements will remain
c     low numbered.  This will be achieved by assigning IE+0.1 to the
c     new element numbers, and then sorting the list.

      Neln=Nel
      DO 101 IE=1,NEL
         ENEW(IE)=IE
         ke = 0
         DO 100 Ioct = 2,Noct
            Neln=Neln+1
            Enew(Neln) = float(ie) + 0.1
            ke = ke+1
            liste(ke) = Neln
  100    continue
         liste(Noct) = Neln+1
         if (ifoct) then
c           call octsplite(ie,liste)  ! Standard oct/quad split
            call octsplitn(ie,liste)  ! New oct/quad split
         elseif (ifzsp) then
            call zsplite(ie,liste) 
         elseif (ifhex) then
            call hsplite(ie,liste) 
         else
            call qsplite(ie,liste) 
         endif
  101 continue

C     Generate new sub-elements as a result of the oct-split action.
      Nnew = Neln-Nel
      nel  = neln
      call sort   (enew,ind,nel)     !     Elements are all set. Sort.
      call swapel (ind,nel,nelm1)
      call curcnt
      call vertadj

      write(s,300) Nnew,NEL
  300 format(i11,' new elements generated in octspl. NEL=',i11,'$')
      call prs(s)

      call gencen

      return
      end
c-----------------------------------------------------------------------
      subroutine wind3d
C
C     Generate a set of corner window refinements in 3D
C
      include 'basics.inc'
      common /cenewi/ enew(nelm),ind(nelm)
      DIMENSION XPT(3),VEC(3)
      common /wierd/ list1(nelm),list2(nelm)
      DIMENSION VEC1(3),VEC2(3),VEC3(3)
      CHARACTER*1 CBCNXT
C
C
      NX=NXM
      NY=NX
      NZ=NX
      call set_zgml(nxm)

      call setquad 
C
C     Set RST planes will find the set of dividing planes for the zipper.
C
C     First, get the point of interest
      WRITE(S,10)
   10 FORMAT(' Enter a point near the corner to be refined.$')
      CALL GETPT3(XPT,IE,S)
      if (IE.EQ.0) return
C
C     Second, find the defining plane
      WRITE(S,11)
      CALL PRS(S)
   11 FORMAT(' Enter normal vector defining the surface$')
      CALL RERRR(VEC(1),VEC(2),VEC(3))
C
C     This normal, plus the specified point can be used to define
C     two normals which determine "splitting planes".  The intersection
C     of these splitting planes yields the requisite corner refinement:
C
C        +---------------------+
C        | .                   |      (The specified normal would be in/out
C        |   .                 |       of the page in this figure.)
C        |     .         ^     |
C        |       .       | v2  |
C        |         +-----------| splitting plane 1
C        |         |           |
C        |    <----|           |
C        |     v1  |     X     |          X is the specified point
C        |         |           |
C        +---------------------+
C                  ^
C                  |
C                  splitting plane 2
C
C
C     Therefore, we need to call SETRST2 twice, with vectors
C     v1 and v2 which are orthogonal to the "given" vector, vec.
C
      CALL FNDPLN(IJKPLN,IPLAN,VEC3,VEC,XPT,IE)
      CALL COPY(VEC,VEC3,3)
      IPLN0=IJKPLN
      IE0  =IE
c     write(6,*) 'ie0:',ie0,ipln0
C
C     For a given IJKPLN (taking values 1-6), we can find the orthogonal
C     planes by adding 2, and taking mod1(ijkpln+2,6).  NOTE that MOD1 is
C     like mod, except that MOD1(6,6)=6 instead of 0.
C
C     Compute v1:
C
      IPLN=IJKPLN+2
      IPLN=MOD1(IPLN,6)
C
      CALL SUB3(VEC1,XYZQ(1,2,IPLN,IE),XYZQ(1,1,IPLN,IE),3)
      CALL SUB3(VEC2,XYZQ(1,4,IPLN,IE),XYZQ(1,3,IPLN,IE),3)
      CALL CROSS(VEC3,VEC1,VEC2)
C     Ok, now let's find all corresponding planes
      CALL SETRST2(XPT,VEC3)
C     Save the list.
      CALL ICOPY(LIST1,LIST,NLSTP)
      NLSTP1=NLSTP
C
C     v2 is easy, given v1 (VEC3), and VEC:
C
      CALL CROSS(VEC2,VEC3,VEC)
C     Ok, now let's find all corresponding planes for v2
      CALL SETRST2(XPT,VEC2)
      CALL ICOPY(LIST2,LIST,NLSTP)
      NLSTP2=NLSTP
C
C     The elements we want lie in the intersection of LIST1 and LIST2
C
      NLSTP =0
      DO 30 IL2=1,NLSTP2
         LISTA=ABS(LIST2(IL2))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         DO 20 IL1=1,NLSTP1
            LISTB=ABS(LIST1(IL1))
            CALL DECOD(IPLANE,JPLN,JE,IDUM,LISTB,NX,6,NELM)
            if (IE.EQ.JE) then
C              Save information about both planes, as they define the split.
               NLSTP=NLSTP+1
               LIST(NLSTP)=6*6*(IE-1)+6*(IPLN-1)+JPLN
      write(6,*) 'ieijp0:',ie,ipln,jpln,nlstp
               GOTO 21
            ENDIF
   20    continue
   21    continue
   30 continue
C
C     Unfortunately, even this list is still too large, so we
C     restrict the elements to the ones which are indeed adjacent
C     to the set we want:
C
C     Set up E-E boundary conditions
      IFLD=2
      if (IFFLOW) IFLD=1
c     if (ifld.lt.10) goto 917
      CALL SETEBC(IFLD)
C
C     Find list of adjacent elements in the direction we want
      NLSTP1=0
      IELCUR=IE0
      NLSTP1=NLSTP1+1
      LIST1(NLSTP1)=IELCUR
      DO 120 IDIR=1,2
         IELCUR=IE0
         IFCCUR=2*((IPLN0-1)/2)+IDIR
         IFCCUR=EFACE(IFCCUR)
         DO 100 I=1,NLSTP
            CBCNXT=CBC(IFCCUR,IELCUR,IFLD)
            write(6,*) 'cbc',cbcnxt,ifccur,ielcur,idir,ifld
            if (CBCNXT.NE.'E') GOTO 101
C           Else, we found another element
            IFCNXT=INT(BC(2,IFCCUR,IELCUR,IFLD))
            IELCUR=INT(BC(1,IFCCUR,IELCUR,IFLD))
C
C           Find next face
            IFCNXT=EFACE1(IFCNXT)
            IFCCUR=IFCNXT+2*MOD(IFCNXT,2)-1
            IFCCUR=EFACE(IFCCUR)
C
C           Bump up list count
            NLSTP1=NLSTP1+1
            LIST1(NLSTP1)=IELCUR
  100    continue
  101    continue
  120 continue
C
C     Compare E-E
      CALL ICOPY(LIST2,LIST,NLSTP)
      NLSTP2=NLSTP
C
C     Again, the elements we want lie in the intersection of LIST1 and LIST2
C
      NLSTP =0
      DO 230 IL2=1,NLSTP2
         LISTA=ABS(LIST2(IL2))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         DO 220 IL1=1,NLSTP1
            JE=LIST1(IL1)
            if (IE.EQ.JE) then
C              Save this element in the list
               NLSTP=NLSTP+1
               LIST(NLSTP)=LIST2(IL2)
               GOTO 221
            ENDIF
  220    continue
  221    continue
  230 continue
C
C     Check on the number and numbering of the new elements to be generated.
C
  917 continue
      NELN=NEL+2*NLSTP
      NELM1=NELM-1
      if (NELN.GT.NELM1) then
C
         WRITE(S,150) 
         CALL PRS(S)
         WRITE(S,151) NELN,NELM1
         CALL PRS(S)
         WRITE(S,152) 
         CALL PRS(S)
  150  FORMAT(2X,'WARNING: Number of elements after zipper operation$')
  151  FORMAT(2X,'(',i11,') would be greater than the allowed maxium ('
     $                ,i11,').$')
  152    FORMAT(2X,'Aborting corner refine operation.$')
         return
      ENDIF
C
C     Renumber ALL elements, so that low numbered elements will remain
C     low numbered.  This will be achieved by assigning IE+1/3,2/3 to the
C     new element numbers, and then sorting the list.
C
      DO 400 IE=1,NEL
         ENEW(IE)=IE
  400 continue
      NELN=NEL
      DO 401 I=1,NLSTP
         LISTA=ABS(LIST(I))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         NELN=NELN+1
         ENEW(NELN)=FLOAT(IE)+0.3333
         NELN=NELN+1
         ENEW(NELN)=FLOAT(IE)+0.6666
  401 continue
C
C     Generate new sub-elements as a result of the zipper action.
C
C     Third, find the ratio to be split
      WRITE(S,501)
      CALL PRS(S)
  501 FORMAT(' Enter ratio (0=small corner element, negative-abort).$')
      CALL RER(S0)
      DO 500 I=1,NLSTP
         NELN=NEL+2*I-1
         LISTA=ABS(LIST(I))
         CALL DECOD(JPLN,IPLN,IE,IDUM,LISTA,6,6,NELM)
         WRITE(S,502) IE,NELN
         CALL PRS(S)
         CALL WIN3D(NELN,IE,I,S0)
  500 continue
  502 FORMAT(' Generating window for element',2i11,'.$')
      NEL=NELN+1
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
C
C     Exit
C
      NLSTP2=NLSTP*2
      WRITE(S,900) NLSTP2,NEL
  900 FORMAT(i11,' new elements generated in win3d. NEL=',i11,'$')
      CALL PRS(S)
C
      call gencen
C
      return
      END
c-----------------------------------------------------------------------
      subroutine zip2
      include 'basics.inc'
      DIMENSION ENEW(NELM),IND(NELM)
      DIMENSION XPT(3),VEC(3)
C
C     This value can be changed later.
C
      NX=NXM
      NY=NX
      NZ=NX
      call set_zgml (nx)

      CALL SETQUAD 
C
C     Set RSTP planes will find the set of dividing planes for the zipper.
C
C     First, get the point of interest
      WRITE(S,10)
   10 FORMAT(' Enter a point near the dividing surface.$')
      CALL GETPT3(XPT,IE,S)
      if (IE.EQ.0) return
C
C     Second, find the defining plane
      WRITE(S,11)
      CALL PRS(S)
      WRITE(S,12)
      CALL PRS(S)
   11 FORMAT(' Enter the components of a normal vector defining the$')
   12 FORMAT(' surface, e.g.,  1 0 0.$')
      CALL RERRR(VEC(1),VEC(2),VEC(3))
C
C     Ok, now let's find all corresponding planes
      CALL SETRST2(XPT,VEC)
C
C     Check on the number and numbering of the new elements to be generated.
C
      NELN=NEL+NLSTP
      NELM1=NELM-1
      if (NELN.GT.NELM1) then
C
         WRITE(S,50) 
         CALL PRS(S)
         WRITE(S,51) NELN,NELM1
         CALL PRS(S)
         WRITE(S,52) 
         CALL PRS(S)
   50  FORMAT(2X,'WARNING: Number of elements after zipper operation$')
   51  FORMAT(2X,'(',i11,') would be greater than the allowed maxium ('
     $                ,i11,').$')
   52    FORMAT(2X,'Aborting zipper.$')
         return
      ENDIF
C
C     Renumber ALL elements, so that low numbered elements will remain
C     low numbered.  This will be achieved by assigning IE+0.5 to the
C     new element numbers, and then sorting the list.
C
      DO 100 IE=1,NEL
         ENEW(IE)=IE
  100 continue
      DO 101 I=1,NLSTP
         NELN=NEL+I
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         ENEW(NELN)=FLOAT(IE)+0.5
  101 continue
C
C     Generate new sub-elements as a result of the zipper action.
C
      S0=0.5
      S1=1.0-S0
      I0 = 0
      DO 200 I=1,NLSTP
         NELN=NEL+I
         LISTA=ABS(LIST(I))
         CALL DECOD(IPLANE,IPLN,IE,IDUM,LISTA,NX,6,NELM)
         CALL SPLITE(NELN,IE,I,I0,S1)
  200 continue
      NEL=NELN
C
C     Elements are all set. Sort.
C
      CALL SORT(ENEW,IND,NEL)
      CALL SWAPEL(IND,NEL,NELM1)
      CALL CURCNT
C
C     Exit
C
      WRITE(S,300) NLSTP,NEL
  300 FORMAT(i11,' new elements generated in zip2. NEL=',i11,'$')
      CALL PRS(S)
C
      call gencen
C
      return
      END
c-----------------------------------------------------------------------
      subroutine swapel(ind,n,nelm1)
      DIMENSION IND(1)
C***
C***  SWAP ELEMENTS  NOTE:  Index tells us that the 1st entry
C***                        can be found in slot IND(1), etc.
C***
C
      DO 100 I=1,N
C
C        Put element IE into slot I
C
         IE=IND(I)
         if (IE.NE.I) then
C
            call copyel(i,nelm1)  ! Move element I to the end
            call copyel(ie,i  )   ! Move correct element to Ith location
            call copyel(nelm1,ie) ! Copy temporary data to old slot.
C           Update IND:
            DO 20 J=1,N
               if (IND(J).EQ.I) then
                   IND(J)=IE
                   GOTO 21
               ENDIF
   20       continue
   21       continue
         ENDIF
  100    continue
      return
      END
c-----------------------------------------------------------------------
      subroutine gh_face_extend_2d(x,zg,n,gh_type,e,v)
c
c     Extend 2D faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and faces
c
c
      real x(n,n)
      real zg(n)
      real e(n,n)
      real v(n,n)
      integer gh_type
c
c     Build vertex interpolant
c
      ntot=n*n
      call rzero(v,ntot)
      do jj=1,n,n-1
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            si     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j) = v(i,j) + si*sj*x(ii,jj)
         enddo
         enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif


c     Extend 4 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do jj=1,n,n-1
         do j=1,n
         do i=1,n
            hj     = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j) = e(i,j) + hj*(x(i,jj)-v(i,jj))
         enddo
         enddo
      enddo
c
c     y-edges
c
      do ii=1,n,n-1
         do j=1,n
         do i=1,n
            hi     = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            e(i,j) = e(i,j) + hi*(x(ii,j)-v(ii,j))
         enddo
         enddo
      enddo

      call add3(x,e,v,ntot)

      return
      end
c-----------------------------------------------------------------------
      subroutine gh_face_extend(x,zg,n,gh_type,e,v)
c
c     Extend faces into interior via gordon hall
c
c     gh_type:  1 - vertex only
c               2 - vertex and edges
c               3 - vertex, edges, and faces
c
c
      real x(n,n,n)
      real zg(n)
      real e(n,n,n)
      real v(n,n,n)
      integer gh_type
c
c     Build vertex interpolant
c
      ntot=n*n*n
      call rzero(v,ntot)
      do kk=1,n,n-1
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            si       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            sj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            sk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + si*sj*sk*x(ii,jj,kk)
         enddo
         enddo
         enddo
      enddo
      enddo
      enddo
      if (gh_type.eq.1) then
         call copy(x,v,ntot)
         return
      endif
c
c
c     Extend 12 edges
      call rzero(e,ntot)
c
c     x-edges
c
      do kk=1,n,n-1
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hj*hk*(x(i,jj,kk)-v(i,jj,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     y-edges
c
      do kk=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hk*(x(ii,j,kk)-v(ii,j,kk))
         enddo
         enddo
         enddo
      enddo
      enddo
c
c     z-edges
c
      do jj=1,n,n-1
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            e(i,j,k) = e(i,j,k) + hi*hj*(x(ii,jj,k)-v(ii,jj,k))
         enddo
         enddo
         enddo
      enddo
      enddo
c
      call add2(e,v,ntot)
c
      if (gh_type.eq.2) then
         call copy(x,e,ntot)
         return
      endif
c
c     Extend faces
c
      call rzero(v,ntot)
c
c     x-edges
c
      do ii=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hi       = 0.5*((n-ii)*(1-zg(i))+(ii-1)*(1+zg(i)))/(n-1)
            v(i,j,k) = v(i,j,k) + hi*(x(ii,j,k)-e(ii,j,k))
         enddo
         enddo
         enddo
      enddo
c
c     y-edges
c
      do jj=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hj       = 0.5*((n-jj)*(1-zg(j))+(jj-1)*(1+zg(j)))/(n-1)
            v(i,j,k) = v(i,j,k) + hj*(x(i,jj,k)-e(i,jj,k))
         enddo
         enddo
         enddo
      enddo
c
c     z-edges
c
      do kk=1,n,n-1
         do k=1,n
         do j=1,n
         do i=1,n
            hk       = 0.5*((n-kk)*(1-zg(k))+(kk-1)*(1+zg(k)))/(n-1)
            v(i,j,k) = v(i,j,k) + hk*(x(i,j,kk)-e(i,j,kk))
         enddo
         enddo
         enddo
      enddo
c
      call add2(v,e,ntot)
      call copy(x,v,ntot)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine sphsrf_e(xml,yml,zml,ifce,e,mx,my,mz,xysrf) 
C
C     5 Aug 1988 19:29:52 
C
C     Program to generate spherical shell elements for NEKTON
C     input.  Paul F. Fischer
C
      include 'basics.inc'
      real xml(mx,my,mz),yml(mx,my,mz),zml(mx,my,mz),xysrf(3,mx,mz)
      integer e

      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)
      common /ctmp0/ xcv(3,2,2),vn1(3),vn2(3)
     $     ,x1(3),x2(3),x3(3),dx(3),xcc(8),ycc(8),zcc(8),zhmt(237)

      integer iopp(3),mxx(3)

      call set_hlin(h,mx) ! linear to mx mapping

C     Determine geometric parameters
      mxm1 = mx-1
      mym1 = my-1
      mxy  = mx*mz
      mxy3 = 3*mx*mz
      xctr   = curve(1,ifce,e)
      yctr   = curve(2,ifce,e)
      zctr   = curve(3,ifce,e)
      radius = curve(4,ifce,e)
      iface  = eface1(ifce)

C     Generate (normalized) corner vectors XCV(1,i,j):
      do 10 i=1,8
         xcc(i)=x(i,e)
         ycc(i)=y(i,e)
         zcc(i)=z(i,e)
   10 continue
      call crn3d(xcv,xcc,ycc,zcc,curve(1,ifce,e),iface,e)

C     Generate edge vectors on the sphere RR=1.0, 
C     for (r,s) = (-1,*),(1,*),(*,-1),(*,1)      
      call edg3d(xysrf,xcv(1,1,1),xcv(1,1,2), 1, 1, 1,my,mx,my)
      call edg3d(xysrf,xcv(1,2,1),xcv(1,2,2),mx,mx, 1,my,mx,my)
      call edg3d(xysrf,xcv(1,1,1),xcv(1,2,1), 1,mx, 1, 1,mx,my)
      call edg3d(xysrf,xcv(1,1,2),xcv(1,2,2), 1,mx,my,my,mx,my)

C     Generate intersection vectors for (i,j)
      do 200 j=2,mym1
         call cross(vn1,xysrf(1,1,j),xysrf(1,mx,j))
         DO 200 I=2,MXM1
            call cross(vn2,xysrf(1,i,1),xysrf(1,i,my))
            if (iface.eq.1.or.iface.eq.4.or.iface.eq.5) then
               call cross(xysrf(1,i,j),vn2,vn1)
            else
               call cross(xysrf(1,i,j),vn1,vn2)
            endif
  200 continue

C     Normalize all vectors to the unit sphere.
      do 300 i=1,mxy
         call norm3d(xysrf(1,i,1))
  300 continue

      call cmult(xysrf,radius,mxy3) ! Scale by actual radius
C
C     Add back the sphere center offset
C
      do 400 i=1,mxy
         xysrf(1,i,1)=xysrf(1,i,1)+xctr
         xysrf(2,i,1)=xysrf(2,i,1)+yctr
         xysrf(3,i,1)=xysrf(3,i,1)+zctr
  400 continue
C
C
C     Transpose data, if necessary
C
      if (iface.eq.1.or.iface.eq.4.or.iface.eq.5) then
         do 500 j=1  ,my
         do 500 i=j+1,mx
            tmp=xysrf(1,i,j)
            xysrf(1,i,j)=xysrf(1,j,i)
            xysrf(1,j,i)=tmp
            tmp=xysrf(2,i,j)
            xysrf(2,i,j)=xysrf(2,j,i)
            xysrf(2,j,i)=tmp
            tmp=xysrf(3,i,j)
            xysrf(3,i,j)=xysrf(3,j,i)
            xysrf(3,j,i)=tmp
  500    continue
      endif
C
C     Compute surface deflection and perturbation due to face IFACE
      call dsset2(js1,jf1,jskip1,js2,jf2,jskip2,iface,mx,my,mz)

      iopp(1) = mx-1
      iopp(2) = mx*(my-1)
      iopp(3) = mx*my*(mz-1)
      mxx(1)  = mx
      mxx(2)  = my
      mxx(3)  = mz
      idir    = 2*mod(iface,2) - 1
      ifc2    = (iface+1)/2
      delt    = 0.0
      i=0
      do 700 j2=js2,jf2,jskip2
      do 700 j1=js1,jf1,jskip1
         i=i+1
         jopp = j1 + iopp(ifc2)*idir
         x2(1) = xml(j1,j2,1)
         x2(2) = yml(j1,j2,1)
         x2(3) = zml(j1,j2,1)

         dx(1) = xysrf(1,i,1)-x2(1)
         dx(2) = xysrf(2,i,1)-x2(2)
         dx(3) = xysrf(3,i,1)-x2(3)

         mxs = mxx(ifc2)
         joff = (j1-jopp)/(mxs-1)
         do 600 ix = 2,mxs
            j = jopp + joff*(ix-1)
            zeta = 0.5*(zgml(ix) + 1.0)
            xml(j,j2,1) = xml(j,j2,1)+dx(1)*zeta
            yml(j,j2,1) = yml(j,j2,1)+dx(2)*zeta
            zml(j,j2,1) = zml(j,j2,1)+dx(3)*zeta
  600    continue
  700 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine edg3d(xysrf,x1,x2,i1,i2,j1,j2,mx,my)
C
C     Generate XYZ vector along an edge of a surface.
C
      include 'basics.inc'
      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)
C
      DIMENSION XYSRF(3,MX,MY)
      DIMENSION X1(3),X2(3)
      REAL U1(3),U2(3),UN(3),B(3)
C
C     Normalize incoming vectors
C
      CALL COPY (U1,X1,3)
      CALL COPY (U2,X2,3)
      CALL NORM3D (U1)
      CALL NORM3D (U2)
C
C     Find normal to the plane and tangent to the curve.
C
      CALL CROSS(UN,X1,X2)
      CALL CROSS( B,UN,X1)
      CALL NORM3D (UN)
      CALL NORM3D (B)
C
      CTHETA = DOT(U1,U2,3)
      THETA  = ACOS(CTHETA)
C
      IJ = 0
      DO 200 J=J1,J2
      DO 200 I=I1,I2
         IJ = IJ + 1
         THETAP = 0.5*THETA*(ZGML(IJ)+1.0)
         CTP = COS(THETAP)
         STP = SIN(THETAP)
         DO 200 IV = 1,3
            XYSRF(IV,I,J) = CTP*U1(IV) + STP*B(IV)
  200 continue
      return
      END
c-----------------------------------------------------------------------
      real function dot(v1,v2,n)
C
C     Compute Cartesian vector dot product.
C
      DIMENSION V1(N),V2(N)
C
      SUM = 0
      DO 100 I=1,N
         SUM = SUM + V1(I)*V2(I)
  100 continue
      DOT = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine crn3d(xcv,xc,yc,zc,curve,iface,ie)
      DIMENSION XCV(3,2,2),XC(8),YC(8),ZC(8),CURVE(4)
      DIMENSION INDX(8)
      SAVE      INDX
      DATA      INDX  / 1, 2, 4, 3, 5, 6, 8, 7 /
      DIMENSION INDVTX(4,6)
      SAVE      INDVTX
      DATA      INDVTX  / 1,5,3,7 , 2,4,6,8 , 1,2,5,6  
     $                  , 3,7,4,8 , 1,3,2,4 , 5,6,7,8 /
C
      EPS    = 1.0E-5
      XCTR   = CURVE(1)
      YCTR   = CURVE(2)
      ZCTR   = CURVE(3)
      RADIUS = CURVE(4)
C
      DO 10 I=1,4
         J=INDVTX(I,IFACE)
         K=INDX(J)
         XCV(1,I,1)=XC(K)-XCTR
         XCV(2,I,1)=YC(K)-YCTR
         XCV(3,I,1)=ZC(K)-ZCTR
   10 continue
C
C     Check to ensure that these points are indeed on the sphere.
C
      if (RADIUS.LE.0.0) then
         WRITE(6,20) XCTR,YCTR,ZCTR,IFACE
  20     FORMAT(5X,'ERROR: Sphere of radius zero requested.'
     $       ,/,5X,'EXITING in CRN3Dc',3E12.4,I3)
         CALL EXITT
      ELSE
         DO 40 I=1,4
            RADT=XCV(1,I,1)**2+XCV(2,I,1)**2+XCV(3,I,1)**2
            RADT=SQRT(RADT)
            TEST=ABS(RADT-RADIUS)/RADIUS
            if (TEST.GT.EPS) then
             WRITE(6,30) 
     $       RADT,RADIUS,XCV(1,I,1),XCV(2,I,1),XCV(3,I,1)
   30        FORMAT(5X,'ERROR: Element vertex not on requested sphere.'
     $           ,/,5X,'REPAIRING in CRN3Dd',5F12.7)
c    $           ,/,5X,'EXITING in CRN3Dd',5F12.7)
             WRITE(6,31) IE,IFACE,XCTR,YCTR,ZCTR
   31        FORMAT(5X,'IE,IF,XYZCTR:',2I4,3F12.7)
             WRITE(6,32) (xc(j),yc(j),zc(j),j=1,8)
   32        FORMAT(3f12.7)

c            CALL EXITT

             xvec = radius*xcv(1,i,1)/radt
             yvec = radius*xcv(2,i,1)/radt
             zvec = radius*xcv(3,i,1)/radt

             j=indvtx(i,iface)
             k=indx(j)
             xc(k)=xctr + xvec
             yc(k)=yctr + yvec
             zc(k)=zctr + zvec

            endif

   40    continue
      ENDIF

      do I=1,4
         j=indvtx(i,iface)
         k=indx(j)
         xcv(1,i,1)=xc(k)-xctr
         xcv(2,i,1)=yc(k)-yctr
         xcv(3,i,1)=zc(k)-zctr
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dsset2(js1,jf1,jskip1,js2,jf2,jskip2,iface,mx,my,mz)
      include 'basics.inc'

      call dsset(mx,my,mz)

      js1    = skpdat(iface,1)
      jf1    = skpdat(iface,2)
      jskip1 = skpdat(iface,3)
      js2    = skpdat(iface,4)
      jf2    = skpdat(iface,5)
      jskip2 = skpdat(iface,6)

      return
      end
c-----------------------------------------------------------------------
      subroutine dsset(mx,my,mz)
C
C     Set up arrays ESKIP,SKPDAT,NEDG,NOFFST for new MX,MY,MZ
C
      include 'basics.inc'
      INTEGER MXO,MYO,MZO
      SAVE    MXO,MYO,MZO
      DATA    MXO,MYO,MZO /3*0/
C
C     Check if element surface counters are already set from last call...
C
      if (MXO.EQ.MX.AND.MYO.EQ.MY.AND.MZO.EQ.MZ) return
C
C     else, proceed....
C
      MXO = MX
      MYO = MY
      MZO = MZ
C
C     Assign indices for direct stiffness summation of arbitrary faces.
C
C
C     Y-Z Planes (Faces 1 and 2)
C
      SKPDAT(1,1)=1
      SKPDAT(1,2)=MX*(MY-1)+1
      SKPDAT(1,3)=MX
      SKPDAT(1,4)=1
      SKPDAT(1,5)=MY*(MZ-1)+1
      SKPDAT(1,6)=MY
C
      SKPDAT(2,1)=1             + (MX-1)
      SKPDAT(2,2)=MX*(MY-1)+1   + (MX-1)
      SKPDAT(2,3)=MX
      SKPDAT(2,4)=1
      SKPDAT(2,5)=MY*(MZ-1)+1
      SKPDAT(2,6)=MY
C
C     X-Z Planes (Faces 3 and 4)
C
      SKPDAT(3,1)=1
      SKPDAT(3,2)=MX
      SKPDAT(3,3)=1
      SKPDAT(3,4)=1
      SKPDAT(3,5)=MY*(MZ-1)+1
      SKPDAT(3,6)=MY
C
      SKPDAT(4,1)=1           + MX*(MY-1)
      SKPDAT(4,2)=MX          + MX*(MY-1)
      SKPDAT(4,3)=1
      SKPDAT(4,4)=1
      SKPDAT(4,5)=MY*(MZ-1)+1
      SKPDAT(4,6)=MY
C
C     X-Y Planes (Faces 5 and 6)
C
      SKPDAT(5,1)=1
      SKPDAT(5,2)=MX
      SKPDAT(5,3)=1
      SKPDAT(5,4)=1
      SKPDAT(5,5)=MY
      SKPDAT(5,6)=1
C
      SKPDAT(6,1)=1           + MX*MY*(MZ-1)
      SKPDAT(6,2)=MX          + MX*MY*(MZ-1)
      SKPDAT(6,3)=1
      SKPDAT(6,4)=1
      SKPDAT(6,5)=MY
      SKPDAT(6,6)=1
C
      return
      END
c-----------------------------------------------------------------------
      subroutine rone(x,n)
      real x(1)
      do 10 i=1,n
         x(i)=1.0
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine arcsrf_e(xml,yml,zml,nxl,nyl,nzl,e,isid)
      include 'basics.inc'
      real xml(nxl,nyl,nzl),yml(nxl,nyl,nzl),zml(nxl,nyl,nzl)
      integer e

C     ctmpg is used in this format in several subsequent routines
      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)
      logical ifglj

      call set_hlin(h,nxl) ! linear to mx mapping

      ifglj = .false.
C     if (ifaxis .and. ifrzer(ie) .and. (isid.eq.2 .or. isid.eq.4)) 
C    $ifglj = .true.

      pt1x  = x(isid,e)
      pt1y  = y(isid,e)
      if(isid.eq.4) then
         pt2x = x(1,e)
         pt2y = y(1,e)
      else if(isid.eq.8) then
         pt2x = x(5,e)
         pt2y = y(5,e)
      else
         pt2x = x(isid+1,e)
         pt2y = y(isid+1,e)
      endif

C     Find slope of perpendicular
      radius=curve(1,isid,e)
      gap=sqrt( (pt1x-pt2x)**2 + (pt1y-pt2y)**2 )
      if (abs(2.0*radius).le.gap*1.00001) then
         write(6,10) radius,isid,e,gap
   10    format(//,2x,'ERROR: Too small a radius (',g11.3
     $  ,') specified for side',I2,' of element',I4,':  '
     $  ,g11.3,/,2x,'ABORTING during mesh generation.')
        radius = 1.5*gap
        curve(1,isid,e) = radius
      endif
      xs = pt2y-pt1y
      ys = pt1x-pt2x
C     Make length Radius
      xys=sqrt(xs**2+ys**2)
c     write(6,977) e,isid,pt1x,pt2x,pt1y,pt2y,xys,xs,ys
c 977 format('arc',i3,i3,7f10.4)
C     Find Center
      dtheta = abs(asin(0.5*gap/radius))
      pt12x  = (pt1x + pt2x)/2.0
      pt12y  = (pt1y + pt2y)/2.0
      xcenn  = pt12x - xs/xys * radius*cos(dtheta)
      ycenn  = pt12y - ys/xys * radius*cos(dtheta)
      theta0 = atan2((pt12y-ycenn),(pt12x-xcenn))
      if (ifglj) then
         fac    = sign(1.0,radius)
         theta1 = theta0 - fac*dtheta
         theta2 = theta0 + fac*dtheta
      endif
C     Compute perturbation of geometry
      ISID1 = MOD1(ISID,4)
      myflag = 0
      if (ifglj) myflag=99
      if (ifglj) then
         i1 = isid/2
         i2 = 2 - isid/4
         do 15 iy=1,nyl
           ang  = h(iy,2,i1)*theta1 + h(iy,2,i2)*theta2
           xcrved(iy)=xcenn + abs(radius)*cos(ang)
     $                      - (h(iy,2,i1)*pt1x + h(iy,2,i2)*pt2x)
           ycrved(iy)=ycenn + abs(radius) * sin(ang)
     $                      - (h(iy,2,i1)*pt1y + h(iy,2,i2)*pt2y)
   15    continue
      else
         do 20 ix=1,nxl
            ixt=ix
            if (isid1.gt.2) ixt=nxl+1-ix
            r=zgml(ix)
            if (radius.lt.0.0) r=-r
            xcrved(ixt) = xcenn + abs(radius) * cos(theta0 + r*dtheta)
     $                          - ( h(ix,1,1)*pt1x + h(ix,1,2)*pt2x )
            ycrved(ixt) = ycenn + abs(radius) * sin(theta0 + r*dtheta)
     $                          - ( h(ix,1,1)*pt1y + h(ix,1,2)*pt2y )
c           if (ix.eq.1) then
c            call color(15)
c            call move3(xcrved(ixt),ycrved(ixt),0.)
c            write(6,*) ixt,xcrved(ixt),ycrved(ixt),' xcrv',radius,isid
c           else
c            call draw3(xcrved(ixt),ycrved(ixt),0.)
c            write(6,*) ixt,xcrved(ixt),ycrved(ixt),' xcrv',radius,isid
c           endif
   20    continue
      endif
c     call prsi('KEEP going?$',myflag)
c     call res(ans,1)

c     call color(15)
c     call move3(xcrved,ycrved,0.)
c     do ix=2,nxl
c        call draw3(xcrved(ix),ycrved(ix),0.)
c        write(6,*) xcrved(ix),ycrved(ix),' xcrv',radius,isid
c     enddo
c     call prsi('keep going?$',myflag)
c     call res(ans,1)


C     Points all set, add perturbation to current mesh.
      isid1 = mod1(isid,4)
      isid1 = eface1(isid1)
      izt = (isid-1)/4+1
      iyt = isid1-2
      ixt = isid1
      if (isid1.le.2) then
         call addtnsr(xml,h(1,1,ixt),xcrved,h(1,3,izt),nxl,nyl,nzl)
         call addtnsr(yml,h(1,1,ixt),ycrved,h(1,3,izt),nxl,nyl,nzl)
      else
         call addtnsr(xml,xcrved,h(1,2,iyt),h(1,3,izt),nxl,nyl,nzl)
         call addtnsr(yml,ycrved,h(1,2,iyt),h(1,3,izt),nxl,nyl,nzl)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine addtnsr(s,h1,h2,h3,nx,ny,nz)
C
C     Map and add to S a tensor product form of the three functions H1,H2,H3.
C     This is a single element routine used for deforming geometry.
C
      DIMENSION H1(1),H2(1),H3(1)
      DIMENSION S(NX,NY,NZ)
C
      DO 200 IZ=1,NZ
      DO 200 IY=1,NY
         HH = H2(IY)*H3(IZ)
         DO 100 IX=1,NX
            S(IX,IY,IZ)=S(IX,IY,IZ)+HH*H1(IX)
  100    continue
  200 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine setrst2(xpt,vec)
C     Set up list of candidate RST planes
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3)
      COMMON /PFFLG/  ILGRNG,ISCND,IENTRP,INEWTX
      DIMENSION XPT(3),VEC(3)
      DIMENSION VEC1(3),VEC2(3),VEC3(3)
      CHARACTER*1 YESNO
      LOGICAL IFANY,IFTMP
C
C======================================================================
C     Find out the exact plane which is defined by the 
C     incoming point and vector.
C======================================================================
C
C     Find the element closest to this point
C
      call finde(ie,xpt(1),xpt(2),xpt(3))
C
C     Does this correspond most closely to an R, S, or T plane?
C
      CALL FNDPLN(IJKPLN,IPLAN,VEC3,VEC,XPT,IE)
C
C===============================================
C     Start generating LIST of plotting planes
C===============================================
C
C     Initialize pointers to zero
      DO 105 I=1,24*NELM
         XYZQ(5,I,1,1)=0.0
  105 continue
C
      NLSTP=1
      LIST(NLSTP)=6*NX*(IE-1)+NX*(IJKPLN-1)+IPLAN
      IMD=MOD(IJKPLN,2)
C
C     Set the sign of this plane according to requested normal vector:
C
      RLNGTH = DOTPROD(VEC,VEC3)
      INRM=1
      if (RLNGTH.LT.0.0) INRM=-1
      LIST(NLSTP)=INRM*LIST(NLSTP)
C
      DO 110 I=1,4
         NEIGHI=INT(XYZQ(4,I,IJKPLN,IE))
         if (NEIGHI.EQ.0) then
            XYZQ(5,I,IJKPLN,IE)=2.0
         ELSE
            RN=FLOAT(NEIGHI)
            XYZQ(5,I,IJKPLN,IE)=SIGN(1.0,RN)
         ENDIF
  110 continue
C
C=======================================
C     Find all adjoining planes
C=======================================
C
      ICOUNT=0
  200 continue
      IFANY=.FALSE.
      ICOUNT=ICOUNT+1
      DO 500 IE=1,NEL
      DO 500 IPLN=1,6
      DO 500 IQ=1,4
         if (ABS(XYZQ(5,IQ,IPLN,IE)).EQ.1.0) then
            IFANY=.TRUE.
            NEIGHI=INT(XYZQ(4,IQ,IPLN,IE))
            NEIGHA=ABS(NEIGHI)
            if (XYZQ(5,NEIGHA,1,1).EQ.0.0) then
C              The neighboring plane has Not been set.
               CALL DECOD(JQ,JPLN,JE,IDUM,NEIGHA,4,6,NEL)
               DO 400 J=1,3
                  JJ=JQ+J
                  JJ=MOD1(JJ,4)
                  NEIGHJ=INT(XYZQ(4,JJ,JPLN,JE))
                  if (NEIGHJ.EQ.0) then
                     XYZQ(5,JJ,JPLN,JE)=2.0
                  ELSE
                     RN=FLOAT(NEIGHJ)
                     XYZQ(5,JJ,JPLN,JE)=
     $                   XYZQ(5,IQ,IPLN,IE)*SIGN(1.0,RN)

                  ENDIF
  400          continue
C
C              New neighbor points all set, now add neighbor plane to LIST.
               JMD=MOD(JPLN,2)
               if (IMD.EQ.JMD) then
                  JPLAN=IPLAN
               ELSE
                  JPLAN=NX+1-IPLAN
               ENDIF
               NLSTP=NLSTP+1
               LIST(NLSTP)=6*NX*(JE-1)+NX*(JPLN-1)+JPLAN
               LIST(NLSTP)=LIST(NLSTP)*SIGN(1.0,XYZQ(5,IQ,IPLN,IE))
C              Flip according to initial plane
               LIST(NLSTP)=INRM*LIST(NLSTP)
C
C              Remove pointers and exit
               XYZQ(5,NEIGHA,1,1)=2.0
            ENDIF
            XYZQ(5,IQ,IPLN,IE)=2.0
         ENDIF
  500 continue
      if (IFANY) GOTO 200
C
C     { LIST is complete, sort by zbuff. (? or perhaps after rotation?)}
C     LIST is complete, set normals according to the following:
C
C     If plane is a Y-plane, i.e. 3 or 4, then it must be flipped (?).
C
      return
      END
c-----------------------------------------------------------------------
      subroutine finde(ie,xpt1,ypt1,zpt1)
      include 'basics.inc'
C
C     Find element which is closest to the point xpt1
C
      call gencen
      DSTMIN=10.0E10
      if (IF3D) then
         DO 100 IEL=1,NEL
            DIST = (XCEN(IEL)-XPT1)**2 
     $           + (YCEN(IEL)-YPT1)**2 
     $           + (ZCEN(IEL)-ZPT1)**2 
            if (DIST.LT.DSTMIN) then
               DSTMIN=DIST
               IEMIN=IEL
            ENDIF
  100    continue
      ELSE
         DO 200 IEL=1,NEL
            DIST = (XCEN(IEL)-XPT1)**2 
     $           + (YCEN(IEL)-YPT1)**2 
            if (DIST.LT.DSTMIN) then
               DSTMIN=DIST
               IEMIN=IEL
            ENDIF
  200    continue
      ENDIF
      IE=IEMIN
      return
      END
c-----------------------------------------------------------------------
      subroutine evalsc( x0 , scal , rrl , inew )
C
C     Evaluate a scalar, SCAL, at position RRL and return the result in X0.
C
      include 'basics.inc'
      DIMENSION SCAL(1)
      DIMENSION RRL(3)
      COMMON  /qTMP0q/ HR(NXM),HS(NXM),HT(NXM)
     $               ,HHH(NXM,NXM,NXM)

      integer icalld,nxold
      save    icalld,nxold
      data    icalld /0/
      data    nxold  /0/
c
      REAL    RLXOLD,RLYOLD,RLZOLD
      SAVE    RLXOLD,RLYOLD,RLZOLD
      DATA    RLXOLD,RLYOLD,RLZOLD/3*1.e20/
C

      call set_zgml(nxm)
      
      SUM = 0.0
      if (RRL(1).NE.RLXOLD.OR.RRL(2).NE.RLYOLD.OR.RRL(3).NE.RLZOLD
     $   .OR.ICALLD.EQ.0) then
         ICALLD = 1
         RLXOLD = RRL(1)
         RLYOLD = RRL(2)
         RLZOLD = RRL(3)
         call hlegn(hr,rrl(1),zgml,nx)
         call hlegn(hs,rrl(2),zgml,ny)

c        write(6,*) 'this is rr1:',rrl(1),nx
c        call outmat(zgml,1,nx,'zgml ',nx)
c        call outmat(hr,1,nx,'hr   ',nx)
c        write(6,*) 'this is rr2:',rrl(2),ny
c        call outmat(hs,1,ny,'hs   ',ny)

         if (IF3D) then
            call hlegn(ht,rrl(3),zgml,nz)
            IJK=0
            DO 100 K=1,NZ
            DO 100 J=1,NY
            DO 100 I=1,NX
               IJK=IJK+1
               HHH(IJK,1,1)=HR(I)*HS(J)*HT(K)
               SUM=SUM+SCAL(IJK)*HHH(IJK,1,1)
  100       continue
         ELSE
            IJK=0
            DO 200 J=1,NY
            DO 200 I=1,NX
               IJK=IJK+1
               HHH(IJK,1,1)=HR(I)*HS(J)
               SUM=SUM+SCAL(IJK)*HHH(IJK,1,1)
c              write(6,1) ijk,sum,scal(ijk),hr(i),hs(j),hhh(ijk,1,1)
c  1           format(i5,5f12.7,' sum 1')
 200       continue
         ENDIF
      ELSE
C
C     Reevaluate a new scalar at an old point.
C
         NXYZ=NX*NY*NZ
         DO 300 I=1,NXYZ
            SUM=SUM+SCAL(I)*HHH(I,1,1)
c           write(6,*) i,sum,scal(i),hhh(i,1,1),' sum 2'
  300    continue
      ENDIF
      X0 = SUM
      return
      END
c-----------------------------------------------------------------------
      subroutine chsign(x,n)

      real x(n)

      do i=1,n
         x(i) = -x(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cross(v1,v2,v3) ! Cartesian vector cross product.

      real v1(3),v2(3),v3(3)

      v1(1) = v2(2)*v3(3) - v2(3)*v3(2)
      v1(2) = v2(3)*v3(1) - v2(1)*v3(3)
      v1(3) = v2(1)*v3(2) - v2(2)*v3(1)

      return
      end
c-----------------------------------------------------------------------
      subroutine norm3d(v1)
C
C     Compute Cartesian vector dot product.
C
      DIMENSION V1(3)
C
      VLNGTH = DOTPROD(V1,V1)
      VLNGTH = SQRT(VLNGTH)
      V1(1) = V1(1) / VLNGTH
      V1(2) = V1(2) / VLNGTH
      V1(3) = V1(3) / VLNGTH
C
      return
      END
c-----------------------------------------------------------------------
      subroutine renum
C
C     Renumber elements so that the lowest numbered elements are
C     in the user specified box.
C
      include 'basics.inc'
      DIMENSION IEIN(NELM),IEOUT(NELM) 
      logical iftmp
      NELM1=NELM-1
C 
C     Find out which element are to be renumbered:
      CALL PRS(
     $    'Enter (with mouse) 2 points in element to be deleted,$')
      CALL PRS(
     $    'or, 2 points framing a box containing elements.$')
      CALL PRS(
     $    'Enter in menu area to abort RENUMBER operation.$')
      IFTMP =IFGRID
      IFGRID=.FALSE.
  120 continue
      CALL PRS('Enter 1st point:$')
      CALL MOUSE(XMOUSE,YMOUSE,BUTTON)
      if (XMOUSE.GT.XPHY(1.0)) then
C        look for a keypad input
         CALL PRS('Aborting renumber operation.$')
         CALL BEEP
         IFGRID=IFTMP 
         return
      ELSE
         CALL PRS('Enter 2nd point:$')
         CALL MOUSE(XMOUS2,YMOUS2,BUTTON)
         if (XMOUS2.GT.XPHY(1.0)) then
          CALL PRS('Aborting renumber operation.$')
          CALL BEEP
          IFGRID=IFTMP 
          return
         ENDIF
      ENDIF
C
C     We successfully input 2 points in the build area
C     Generate element centers
      call gencen
      XMAX=MAX(XMOUSE,XMOUS2)
      XMIN=MIN(XMOUSE,XMOUS2)
      YMAX=MAX(YMOUSE,YMOUS2)
      YMIN=MIN(YMOUSE,YMOUS2)
C
C     Check box to see if it contains any elements, and renumber them.
C
C     Count the number which are in
      NUMIN=0
      CALL DRWBOX(XMIN,YMIN,XMAX,YMAX,1)
      DO 100 IIEL=1,NEL
         if (XMIN.LE.XCEN(IIEL) .AND.
     $       XCEN(IIEL).LE.XMAX .AND.
     $       YMIN.LE.YCEN(IIEL) .AND.
     $       YCEN(IIEL).LE.YMAX )         NUMIN=NUMIN+1
  100 continue
      WRITE(S,101) NUMIN,NEL
  101 FORMAT('Renumbering',i11,' out of',i11,' elements.$')
      CALL PRS(S)
      if (NUMIN.GT.0) then
C        renumber the elements which are inside to be less than NUMIN
         J=0
         K=0
         DO 200 IIEL=1,NEL
            if (XMIN.LE.XCEN(IIEL) .AND.
     $          XCEN(IIEL).LE.XMAX .AND.
     $          YMIN.LE.YCEN(IIEL) .AND.
     $          YCEN(IIEL).LE.YMAX ) then
                if (IIEL.GT.NUMIN) then
C               We've got one which is in, but too high.
                   J=J+1
                   IEIN(J)=IIEL
                ENDIF
             ELSE
                if (IIEL.LT.NUMIN) then
C               We've got one which is out, but too low.
                   K=K+1
                   IEOUT(K)=IIEL
                ENDIF
             ENDIF
  200    continue
      WRITE(S,201) J,K
  201 FORMAT('Found',i11,' elements in, and',i11,' elements out.$')
      CALL PRS(S)
      if (J.NE.K) return
C
C        Swap the elements
C
         DO 300 I=1,J
            IE=IEIN(I)
            JE=IEOUT(I)
            CALL COPYEL(IE,NELM1)
            CALL COPYEL(JE,IE)
            CALL COPYEL(NELM1,JE)
  300    continue
C
      ENDIF

      call redraw_mesh

      return
      END
c-----------------------------------------------------------------------
      subroutine getpt3(xpt,ie,text)
C     Set up list of candidate RST planes
      include 'basics.inc'
      CHARACTER*1 YESNO
      CHARACTER*80 TEXT
      DIMENSION XPT(3)
C
      
   10 CALL PRS(TEXT)
      CALL RERRR(XPT(1),XPT(2),XPT(3))
C
C     Find the element closest to this point
C
      call finde(ie,xpt(1),xpt(2),xpt(3))
C
C     Exception handling:
      if (IE.EQ.0) then
         WRITE(S,20) XPT(1),XPT(2),XPT(3)
         CALL PRS(S)
   20    FORMAT(2X,'Point (',G10.3,',',G10.3,',',G10.3,') is not in '
     $            ,'the domain.$')
         CALL PRS('  Continue (y/n)?$')
         CALL RES(YESNO,1)
         if (YESNO.EQ.'Y'.OR.YESNO.EQ.'y') GOTO 10
      ELSE
         WRITE(S,30) IE,XCEN(IE),YCEN(IE),ZCEN(IE)
   30    FORMAT(2X,'Found element number',i11,3f9.4,'$')
         CALL PRS(S)
      ENDIF
C
      return
      END
c-----------------------------------------------------------------------
      subroutine fndpln(ijkpln,iplan,vec3,vec,xpt,ie)
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      DIMENSION VEC(3) ,VEC3(3),XPT(3)
      DIMENSION VEC1(3),VEC2(3)
C
      RLMAX=0.0
      DO 30 IPLN=1,6 
         CALL SUB3(VEC1,XYZQ(1,2,IPLN,IE),XYZQ(1,1,IPLN,IE),3)
         CALL SUB3(VEC2,XYZQ(1,4,IPLN,IE),XYZQ(1,3,IPLN,IE),3)
         CALL CROSS(VEC3,VEC1,VEC2)
         CALL NORM3D(VEC3)
         RLNGTH = DOTPROD(VEC,VEC3)
         RLNGTH = ABS(RLNGTH)
         if (RLNGTH.GT.RLMAX) then
            IJKPLN  = IPLN
            RLMAX   = RLNGTH
            DXORD=VEC3(1) 
            DYORD=VEC3(2) 
            DZORD=VEC3(3) 
         ENDIF
   30 continue
      VEC3(1)=DXORD
      VEC3(2)=DYORD
      VEC3(3)=DZORD
C
      DXORD=VEC(1)
      DYORD=VEC(2)
      DZORD=VEC(3)
C
C     Find the plane corresponding to the chosen point and normal.
C
      nxy =nx*ny
      nxyz=nx*ny*nz
      dmin=10.0e15
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
      do 100 i=1,nxyz
         dist=(xp(i)-xpt(1))**2+(yp(i)-xpt(2))**2+
     $        (zp(i)-xpt(3))**2
         if (dist.lt.dmin) then
            ijkmin=i
            dmin=dist
         endif
  100 continue
      call decod(ixp,iyp,izp,idum,ijkmin,nx,ny,nz)
      if (ijkpln.le.2) iplan=ixp
      if (ijkpln.gt.2) iplan=iyp
      if (ijkpln.gt.4) iplan=izp
C     Set IJKPLN according to actual value of IPLAN
      nxt=2*iplan-1
      if (nxt.lt.nx) ijkpln=2*((ijkpln-1)/2)+1
      if (nxt.gt.nx) ijkpln=2*((ijkpln-1)/2)+2
C
      return
      end
c-----------------------------------------------------------------------
      subroutine win3d(je,ie,ilist,frac)
C
C     A routine which will perform a corner refinement.
C
C     This normal, plus the specified point can be used to define
C     two normals which determine "splitting planes".  The intersection
C     of these splitting planes yields the requisite corner refinement:
C
C        +---------------------+
C        | .                   |
C        |   .                 |
C        |     .       JE      |
C        |       .             |
C        |         +-----------| splitting plane 1
C        |         |           |
C        |   KE    |           |
C        |         |    IE     |
C        |         |           |
C        +---------------------+
C                  ^
C                  |
C                  splitting plane 2
C
C
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      INTEGER   ITRANS(3,6,6)
      CHARACTER CTRANS(3,6,6)
      INTEGER INV(8)
      SAVE    INV,ITRANS,CTRANS
      DATA    INV /1,2,4,3,5,6,8,7/
      DATA    ITRANS /108*0/
      DATA    CTRANS /108*' '/
C
      if (FRAC.LE.0.0.OR.FRAC.GE.1.0) return
C
c     ITRANS(1,1,4)=Z+1
c     ITRANS(1,1,5)=Z+1;Y+1
c     ITRANS(1,1,6)=X+1
c     ITRANS(1,2,3)=Z-1
c     ITRANS(1,2,4)=Z+2
c     ITRANS(1,2,5)=Z+1;Y+1;X+1
c     ITRANS(1,2,6)=X+1;Z-1
c     ITRANS(1,3,5)=Y+1
c     ITRANS(1,3,6)=Y+1;Z-1
c     ITRANS(1,4,5)=Y-1;Z-1
c     ITRANS(1,4,6)=Y-1

      ITRANS(1,1,4)= +1
      ITRANS(1,1,5)= +1
      ITRANS(2,1,5)= +1
      ITRANS(1,1,6)= +1
      ITRANS(1,2,3)= -1
      ITRANS(1,2,4)= +2
      ITRANS(1,2,5)= -1
      ITRANS(2,2,5)= -1
C
      ITRANS(1,2,6)= +1
      ITRANS(2,2,6)= -1
      ITRANS(1,3,5)= +1
      ITRANS(1,3,6)= +1
      ITRANS(2,3,6)= -1
C
c     ITRANS(1,4,5)= -1
c     ITRANS(2,4,5)= -1
c     ITRANS(1,4,6)= -1
      ITRANS(1,4,5)= +2
      ITRANS(2,4,5)= -1
      ITRANS(1,4,6)= +1
      ITRANS(2,4,6)= -1
C
      CTRANS(1,1,4)='Z'
      CTRANS(1,1,5)='Z'
      CTRANS(2,1,5)='Y'
      CTRANS(1,1,6)='X'
C
      CTRANS(1,2,3)='Z'
      CTRANS(1,2,4)='Z'
      CTRANS(1,2,5)='X'
      CTRANS(2,2,5)='Z'
C
      CTRANS(1,2,6)='X'
      CTRANS(2,2,6)='Z'
C
      CTRANS(1,3,5)='Y'
      CTRANS(1,3,6)='Y'
      CTRANS(2,3,6)='Z'
C
c     CTRANS(1,4,5)='Y'
c     CTRANS(2,4,5)='Z'
c     CTRANS(1,4,6)='Y'
      CTRANS(1,4,5)='X'
      CTRANS(2,4,5)='Y'
      CTRANS(1,4,6)='X'
      CTRANS(2,4,6)='Y'
C
      FRACS=2.0*FRAC-1.0
      LISTA=ABS(LIST(ILIST))
      CALL DECOD(JPLN,IPLN,IDUM1,IDUM2,LISTA,6,6,NELM)
      KE=JE+1
      write(6,*) 'ieijp:',ie,ipln,jpln
C
C     Highlight the element to be split
C
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
      call hilite(ie,xp,yp,zp,5)
C
C     Begin with copying the essential information
C
C     Note COPYEL copies IE ==> JE
      CALL COPYEL(IE,JE)
      CALL COPYEL(IE,KE)
C                                                   +----+
C                                                   |ke /|
C     Make rotation(s) to get into canonical form:  |__+ |
C                                                   |  | |
C                                               ie  +----+  je
      I1=MIN(IPLN,JPLN)
      I2=MAX(IPLN,JPLN)
      DO 10 IROT=1,3
         CALL ROTELM(IE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM(JE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM(KE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
   10 continue
C     Compute XYZ on high definition mesh
C             (Note, because of this, ROTELM must rotate CURVE data as well.)
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
C
C     Generate element IE
C
C     First part, IE,KE are lower r-halves, JE is upper r-half.
      I1=-3
      DO 100 KCRN=-1,1,2
         RRL(3)=FLOAT(KCRN)
         I1=I1+4
         I2=I1+1
         I3=I1+2
         I4=I1+3
C
C
C       Notation:   +-------+      
C                   |     .'|      ^ Y
C                   C---B'  |      |
C                   |   |   |      |
C                   +---A---+      +---->  X
C
C     "A" point
         RRL(1)=FRACS
         RRL(2)=-1.0
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv1:',xval,yval,zval

         x(i2,ie)=xval
         y(i2,ie)=yval
         z(i2,ie)=zval

         x(i1,je)=xval
         y(i1,je)=yval
         z(i1,je)=zval
C
C     "B" point
         RRL(1)=FRACS
         RRL(2)=FRACS
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv2:',xval,yval,zval
C
         x(i3,ie)=xval
         y(i3,ie)=yval
         z(i3,ie)=zval
C
         x(i4,je)=xval
         y(i4,je)=yval
         z(i4,je)=zval
C
         x(i2,ke)=xval
         y(i2,ke)=yval
         z(i2,ke)=zval
C
C     "C" point
         RRL(1)=-1.0
         RRL(2)=FRACS
         CALL EVALSC(XVAL,XP,RRL,1)
         CALL EVALSC(YVAL,YP,RRL,0)
         CALL EVALSC(ZVAL,ZP,RRL,0)
         write(6,*) 'xyzv3:',xval,yval,zval
C
         x(i4,ie)=xval
         y(i4,ie)=yval
         z(i4,ie)=zval
c
         x(i1,ke)=xval
         y(i1,ke)=yval
         z(i1,ke)=zval
C
  100 continue
C
C     Don't forget to eliminate Curve side data on new faces!
C
      JFAC1=EFACE(1)
      JFAC2=EFACE(2)
      JFAC3=EFACE(3)
      JFAC4=EFACE(4)
C
      CCURVE(JFAC2,IE)=' '
      CCURVE(JFAC4,IE)=' '
      CCURVE(JFAC1,JE)=' '
      CCURVE(JFAC4,JE)=' '
      CCURVE(JFAC2,JE)=' '
      CCURVE(JFAC3,JE)=' '
C
C     Delete periodic bcs .... will do later, MC?
C
C     Undo rotations via inverse sequence
      I1=MIN(IPLN,JPLN)
      I2=MAX(IPLN,JPLN)
      DO 510 IROT=3,1,-1
         CALL ROTELM1(IE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM1(JE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
         CALL ROTELM1(KE,ITRANS(IROT,I1,I2),CTRANS(IROT,I1,I2))
  510 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine rotelm(je,irot,crot)
      CHARACTER*1 CROT
C
C     If there's no rotation to be done, we return.
      if (IROT.EQ.0) return
C
C     Else, we construct the necessary rotation as a series of unit
C     rotations about the specified axis.
C
      NROT=IROT+8
      NROT=MOD(NROT,4)
C
      DO 100 I=1,NROT
         CALL ROTEL(JE,CROT)
  100 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine rotelm1(je,irot,crot)
      CHARACTER*1 CROT
C
C     Inverse rotation: if there's no rotation to be done, we return.
      if (IROT.EQ.0) return
C
C     Else, we construct the necessary rotation as a series of unit
C     rotations about the specified axis.
C
      NROT=IROT+8
      NROT=MOD(NROT,4)
      if (NROT.EQ.3) then
         NROT=1
      elseif (NROT.EQ.1) then
         NROT=3
      ENDIF
C
      DO 100 I=1,NROT
         CALL ROTEL(JE,CROT)
  100 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine rotel(ie,crot)
C
C     Perform a unit rotation about the specified axis.
C
      include 'basics.inc'
      CHARACTER*1 CROT,CCTMP
      DIMENSION TMP(6)
      INTEGER FCRN2 (4,6)
      DATA    FCRN2 / 1,4,8,5,
     $                2,3,7,6,
     $                1,5,6,2,
     $                4,8,7,3,
     $                1,2,3,4,
     $                5,6,7,8 /
      INTEGER FFAC (4,3)
      DATA    FFAC / 1,5,3,6,
     $               4,6,2,5,
     $               1,2,3,4/
C
C     'X' axis
C
      if (CROT.EQ.'X') then
         IFC1=1
         IFC2=2
         IFAC=1
      elseif (CROT.EQ.'Y') then
         IFC1=3
         IFC2=4
         IFAC=2
      ELSE
         IFC1=5
         IFC2=6
         IFAC=3
      ENDIF
C
      x1=x(fcrn2(4,ifc1),ie)
      y1=y(fcrn2(4,ifc1),ie)
      z1=z(fcrn2(4,ifc1),ie)
      x2=x(fcrn2(4,ifc2),ie)
      y2=y(fcrn2(4,ifc2),ie)
      z2=z(fcrn2(4,ifc2),ie)
      CCTMP=CCURVE(FFAC(4,IFAC),IE)
      CALL COPY(TMP,CURVE(1,FFAC(4,IFAC),IE),6)
      DO 10 I=3,1,-1
         I1=I+1
         x(fcrn2(i1,ifc1),ie)=x(fcrn2(i,ifc1),ie)
         y(fcrn2(i1,ifc1),ie)=y(fcrn2(i,ifc1),ie)
         z(fcrn2(i1,ifc1),ie)=z(fcrn2(i,ifc1),ie)
         x(fcrn2(i1,ifc2),ie)=x(fcrn2(i,ifc2),ie)
         y(fcrn2(i1,ifc2),ie)=y(fcrn2(i,ifc2),ie)
         z(fcrn2(i1,ifc2),ie)=z(fcrn2(i,ifc2),ie)
C
         CCURVE(FFAC(I1,IFAC),IE)=CCURVE(FFAC(I,IFAC),IE)
         CALL COPY(CURVE(1,FFAC(I1,IFAC),IE)
     $            ,CURVE(1,FFAC(I,IFAC),IE),6)
   10 continue
      x(fcrn2(1,ifc1),ie)=x1
      y(fcrn2(1,ifc1),ie)=y1
      z(fcrn2(1,ifc1),ie)=z1
      x(fcrn2(1,ifc2),ie)=x2
      y(fcrn2(1,ifc2),ie)=y2
      z(fcrn2(1,ifc2),ie)=z2
      CCURVE(FFAC(1,IFAC),IE)=CCTMP
      CALL COPY(CURVE(1,FFAC(1,IFAC),IE),TMP,6)
C
      return
      END
c-----------------------------------------------------------------------
      subroutine hilite(ie,xp,yp,zp,iclr)
C
C     Draw ISOMETRIC view of general 3D object.
C
      include 'basics.inc'
      CHARACTER*1 CB1
      INTEGER ICALLD,NCSOLD,NCSEG1,NCSGM1
      INTEGER ICEDG(3,16),IEDGFC(4,6),ICFACE(4,10)
      SAVE    ICEDG,IEDGFC,ICFACE
      SAVE    ICALLD,NCSOLD,NCSEG1,NCSGM1
      DIMENSION XP(1),YP(1),ZP(1)
C
      COMMON /HLINE/ HIK(NXM,50),XXIS(50),YYIS(50),ZZIS(50),
     $               XK(50),YK(50),ZK(50)
      COMMON /ILINE/ IC(8),NSKIPX(3)
      COMMON /LLINE/ IFEDG(12,2,NELM)
      LOGICAL IFEDG,IFCURV
      INTEGER NXOLD,NYOLD,NZOLD
      SAVE    NXOLD,NYOLD,NZOLD
      DATA    NXOLD,NYOLD,NZOLD /0,0,0/
C
      DATA    ICALLD,NCSOLD,NCSEG1,NCSGM1 /4*0/
      DATA    IEDGFC /  5,7,9,11,  6,8,10,12,
     $                  1,3,9,10,  2,4,11,12,
     $                  1,2,5,6,   3,4,7,8    /
      DATA    ICEDG / 1,2,1,   3,4,1,   5,6,1,   7,8,1,
     $                1,3,2,   2,4,2,   5,7,2,   6,8,2,
     $                1,5,3,   2,6,3,   3,7,3,   4,8,3,
c      -2D-  (NOTE:  Edges default to faces in 2D! Not consistent, but...)
     $                1,3,2,   2,4,2,   1,2,1,   3,4,1 /
      DATA    ICFACE/ 1,3,5,7, 2,4,6,8,
     $                1,2,5,6, 3,4,7,8,
     $                1,2,3,4, 5,6,7,8,
c      -2D-
     $                1,3,0,0, 2,4,0,0,
     $                1,2,0,0, 3,4,0,0  /
C
C
C     Sort the faces according to Z-location
C     for the color fill plotting routines
c     if (IF3D) CALL SORTZ
C
      NFACES=2*NDIM
      NEDGE=4
      if (IF3D) NEDGE=12
      IF(ICALLD.EQ.0)THEN
         ICALLD=1
C
C        Set up logicals determining whether edge is to be
C        displayed or not.
C
C        Assign EB's numbering scheme to PF's scheme.
C
         EFACE(1)=4
         EFACE(2)=2
         EFACE(3)=1
         EFACE(4)=3
         EFACE(5)=5
         EFACE(6)=6
C
C        Assign inverse of EB's numbering scheme to PF's scheme.
C
         EFACE1(1)=3
         EFACE1(2)=2
         EFACE1(3)=4
         EFACE1(4)=1
         EFACE1(5)=5
         EFACE1(6)=6
C
      ENDIF
      if (NX.NE.NXOLD.OR.NY.NE.NYOLD.OR.NZ.NE.NZOLD) then
         NXOLD=NX
         NYOLD=NY
         NZOLD=NZ
         NXYZ = NX*NY*NZ
         NSKIPX(1) = 1
         NSKIPX(2) = NX
         NSKIPX(3) = NX*NY
C
C        Fill corners - 1 through 8.
C
         IC(1) =                                 0
         IC(2) =                             (NX-1)
         IC(3) =                 NX*(NY-1)
         IC(4) =                 NX*(NY-1) + (NX-1)
         IC(5) =  NX*NY*(NZ-1)
         IC(6) =  NX*NY*(NZ-1)             + (NX-1)
         IC(7) =  NX*NY*(NZ-1) + NX*(NY-1)
         IC(8) =  NX*NY*(NZ-1) + NX*(NY-1) + (NX-1)
         if (NCSOLD.NE.NCSEGS) then
C           Set up curve line segment evaluation array (exclude endpoints) :
            NCSOLD = NCSEGS
            NCSEG1 = NCSEGS + 1
            NCSGM1 = NCSEGS - 1
            DR1 = 2.0 / FLOAT(NCSEGS)
            DO 400 K=1,NCSGM1
               RR1 = -1.0 + DR1*FLOAT(K)
               call hlegn(hik(1,k),rr1,zgml,nx)
  400       continue
         ENDIF
      ENDIF
C
C  ==============================================
C     Preliminaries all set, draw the element
C  ==============================================
C
C
      CALL COLOR(ICLR)
      CALL PENW(1)
C
C     Draw element edges (for now, we do the 4X redundant work)
C
         IEOFF=1
C
C        Check if we can take a short cut in plotting the element bdry.
         IFCURV=.FALSE.
         DO 500 IEDG = 1,12
            if (CCURVE(IEDG,IE).NE.' ') IFCURV=.TRUE.
  500    continue
C
C        Edges
C
         DO 800 IEDG = 1,NEDGE
            IEDG0 = IEDG + 12
            if (IF3D) IEDG0 = IEDG
            ICE1 = ICEDG(1,IEDG0)
            ICE2 = ICEDG(2,IEDG0)
            NSKIP = NSKIPX(ICEDG(3,IEDG0))
            IC0 = IEOFF+IC(ICE1)
            IC1 = IEOFF+IC(ICE2)
            XK(     1) = XP(IC0)
            YK(     1) = YP(IC0)
            ZK(     1) = ZP(IC0)
            XK(NCSEG1) = XP(IC1)
            YK(NCSEG1) = YP(IC1)
            ZK(NCSEG1) = ZP(IC1)
            if (IFCURV) then
               CALL EVLINE(XK(2),XP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL EVLINE(YK(2),YP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL EVLINE(ZK(2),ZP(IC0),HIK,NX,NXM,NCSGM1,NSKIP)
               CALL VDRAW3A(XK,YK,ZK,NCSEG1)
            ELSE
               XK(2) = XP(IC1)
               YK(2) = YP(IC1)
               ZK(2) = ZP(IC1)
               CALL VDRAW3A(XK,YK,ZK,2)
            ENDIF
  800    continue
 1000 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine evline(uv,uvp,hik,n1,n1d,npt,nskip)
      include 'basics.inc'
C
C     Evaluate vector uv <== uvp, @ k=1,npts
C
      REAL UV(NPT),UVP(1)
      REAL HIK(N1D,NPT)
      if (NPT.LE.N1) then
         DO 100 K=1,NPT
            UV(K) = 0.0
            ISKIP = 1 - NSKIP
            DO 100 I=1,N1
               ISKIP = ISKIP + NSKIP
               UV(K) = UV(K) + HIK(I,K)*UVP(ISKIP)
  100    continue
      ELSE
         DO 200 K=1,NPT
            UV(K)=0.0
  200    continue
         ISKIP = 1 - NSKIP
         DO 300 I=1,N1
            ISKIP = ISKIP + NSKIP
            DO 300 K=1,NPT
               UV(K) = UV(K) + HIK(I,K)*UVP(ISKIP)
  300    continue
      ENDIF
      return
      END
c-----------------------------------------------------------------------
      subroutine setebc(ifld)
      include 'basics.inc'
C     Find Sides' Midpoints
      call mkside
      call gencen
C     Find Sides' Overlaps (used for finding boundary sides)
C!! Check here for Irrational B.c.'s: unrequited connectivities;
C     Warn if physical boundaries internally.  This check is not exhaustive!!??
C     ZERO OUT OLD ELEMENTAL CONNECTIVITY
      DO 300 IE=1,NEL
      DO 300 ISIDE=1,NSIDES
         IF(  CBC( ISIDE, IE,IFLD).EQ.'E')THEN
              CBC(   ISIDE, IE,IFLD)=' '
              ibc(   iside, ie,ifld)= 0
              BC (1, ISIDE, IE,IFLD)= 0
              BC (2, ISIDE, IE,IFLD)= 0
              BC (3, ISIDE, IE,IFLD)= 0
         ENDIF
  300 continue
C
C     Loop over all elements to find element-element interactions
C
      DO 500 IE=1,NEL
       DELTA1 = .001*RCEN(IE)
       DO 500 JE=1,NEL
C
         DIST2=(XCEN(IE)-XCEN(JE))**2+(YCEN(IE)-YCEN(JE))**2
     $        +(ZCEN(IE)-ZCEN(JE))**2
         DISTR=(RCEN(IE)+RCEN(JE))**2
         if (DIST2.LE.DISTR) then
C           There is a potential interaction
C
             DELTA2 = .001*RCEN(JE)
             DELTA  = MIN(DELTA1,DELTA2)
             DO 400 ISIDE=1,NSIDES
             DO 400 JSIDE=1,NSIDES
                IF(IE.EQ.JE.AND.ISIDE.EQ.JSIDE) GOTO 400
                DELTAX = ABS(SIDES(IE,ISIDE,1)-SIDES(JE,JSIDE,1))
                DELTAY = ABS(SIDES(IE,ISIDE,2)-SIDES(JE,JSIDE,2))
                DELTAZ = ABS(SIDES(IE,ISIDE,3)-SIDES(JE,JSIDE,3))
                IF(DELTAX .LT. DELTA .AND.
     $             DELTAY .LT. DELTA .AND.
     $             DELTAZ .LT. DELTA ) then
C                  BC Array used to define neighboring elements
C                  For want of better notation, 
C                 'E' means elemental (internal) bdry
C                  1st two reals are element & side of neighbor; 
C                  3rd is orientation
C
C                  Re-do internal connectivity
C                  "Normal" internal side.
                   CBC3=CBC(ISIDE,IE,IFLD)
                   IF(CBC3(3:3).NE.'I'.AND.CBC3(3:3).NE.'i')THEN
C                    Overlapping edges not Internal B.C.'s are elemental b.c.'s
                     CBC(ISIDE,IE,IFLD)='E'
                     CBC(JSIDE,JE,IFLD)='E'
                   ENDIF
                   IORIEN = 0
                   ibc(iside,ie,ifld)   = je
                   BC (1,ISIDE,IE,IFLD) = JE
                   BC (2,ISIDE,IE,IFLD) = JSIDE
                   BC (3,ISIDE,IE,IFLD) = IORIEN

                   ibc(jside,je,ifld)   = ie
                   BC (1,JSIDE,JE,IFLD) = IE
                   BC (2,JSIDE,JE,IFLD) = ISIDE
                   BC (3,JSIDE,JE,IFLD) = IORIEN
                ENDIF
  400       continue
         ENDIF
  500 continue
      return
      END
c-----------------------------------------------------------------------
      subroutine psphsrf_e(xml,yml,zml,ifce,e,mx,my,mz,xysrf) 
C
C     5 Aug 1988 19:29:52 
C     Program to generate spherical shell elements for NEKTON
C     input.  Paul F. Fischer
C
C     Updated 7/30/92 to compute prolate spheroids
C
      include 'basics.inc'
      real xml(mx,my,mz),yml(mx,my,mz),zml(mx,my,mz)
      real xysrf(3,mx,mz)
      integer e

      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)
      common /ctmp0/ xcv(3,2,2),vn1(3),vn2(3)
     $     ,x1(3),x2(3),x3(3),dx(3),xcc(8),ycc(8),zcc(8),zhmt(237)
      dimension iopp(3),mxx(3)

      call set_hlin(h,mx) ! linear to mx mapping

C     Determine geometric parameters
      MXM1 = MX-1
      MYM1 = MY-1
      MXY  = MX*MZ
      MXY3 = 3*MX*MZ
      xctr   = curve(1,ifce,e)
      yctr   = curve(2,ifce,e)
      zctr   = curve(3,ifce,e)
      radius = curve(4,ifce,e)
C
C     For now, we support only stretching in X, for lack of better
C     interface... the code is general however.
C
      Pratio = curve(5,ifce,e)
      PratiI = 1.0/Pratio
      iface  = eface1(ifce)
C
C     First, we map to unstretched coordinates:
C
C       ^      +-----------------+              ^      +-----+
C     y |      |                 |           y' |      |     |
C       |      |                 |              |      |     |
C       +--->  +-----------------+              +--->  +-----+
C         x                                      x'
C
C     Generate edge vectors on the sphere RR=1.0, 
C
C     Generate (normalized) corner vectors XCV(1,i,j),
C     in unstretched coordinates:
C
C       ^      +-----------------+              ^      +-----+
C     y |      |                 |           y' |      |     |
C       |      |                 |              |      |     |
C       +--->  +-----------------+              +--->  +-----+
C         x                                      x'
C
      do 10 i=1,8
         xcc(i)=x(i,e)*PratiI
         ycc(i)=y(i,e)
         zcc(i)=z(i,e)
   10 continue
      call crn3d(xcv,xcc,ycc,zcc,curve(1,ifce,e),iface,e)
C
C     Generate edge vectors on the sphere RR=1.0, 
C     for (r,s) = (-1,*),(1,*),(*,-1),(*,1)      
C
      call edg3d(xysrf,xcv(1,1,1),xcv(1,1,2), 1, 1, 1,my,mx,my)
      call edg3d(xysrf,xcv(1,2,1),xcv(1,2,2),mx,mx, 1,my,mx,my)
      call edg3d(xysrf,xcv(1,1,1),xcv(1,2,1), 1,mx, 1, 1,mx,my)
      call edg3d(xysrf,xcv(1,1,2),xcv(1,2,2), 1,mx,my,my,mx,my)
C
C     Generate intersection vectors for (i,j)
C
      do 200 j=2,mym1
         call cross(vn1,xysrf(1,1,j),xysrf(1,mx,j))
         do 200 i=2,mxm1
            call cross(vn2,xysrf(1,i,1),xysrf(1,i,my))
            if (iface.eq.1.or.iface.eq.4.or.iface.eq.5) then
               call cross(xysrf(1,i,j),vn2,vn1)
            else
               call cross(xysrf(1,i,j),vn1,vn2)
            endif
  200 continue
C
C     Normalize all vectors to the unit sphere.
C
      do 300 i=1,mxy
         call norm3d(xysrf(1,i,1))
  300 continue
C
C     Scale by actual radius
C
      call cmult(xysrf,radius,mxy3)
C
C     Add back the sphere center offset
C
      do 400 i=1,mxy
         xysrf(1,i,1)=xysrf(1,i,1)+xctr
         xysrf(2,i,1)=xysrf(2,i,1)+yctr
         xysrf(3,i,1)=xysrf(3,i,1)+zctr
  400 continue
C
C
C     Transpose data, if necessary
C
      if (iface.eq.1.or.iface.eq.4.or.iface.eq.5) then
         do 500 j=1  ,my
         do 500 i=j+1,mx
            tmp=xysrf(1,i,j)
            xysrf(1,i,j)=xysrf(1,j,i)
            xysrf(1,j,i)=tmp
            tmp=xysrf(2,i,j)
            xysrf(2,i,j)=xysrf(2,j,i)
            xysrf(2,j,i)=tmp
            tmp=xysrf(3,i,j)
            xysrf(3,i,j)=xysrf(3,j,i)
            xysrf(3,j,i)=tmp
  500    continue
      endif
C
C     Compute surface deflection and perturbation due to face IFACE
      call dsset2(js1,jf1,jskip1,js2,jf2,jskip2,iface,mx,my,mz)

      iopp(1) = mx-1
      iopp(2) = mx*(my-1)
      iopp(3) = mx*my*(mz-1)
      mxx(1)  = mx
      mxx(2)  = my
      mxx(3)  = mz
      idir    = 2*mod(iface,2) - 1
      ifc2    = (iface+1)/2
      delt    = 0.0
      i=0
      do 700 j2=js2,jf2,jskip2
      do 700 j1=js1,jf1,jskip1
         i=i+1
         jopp = j1 + iopp(ifc2)*idir
         x2(1) = xml(j1,j2,1)*PratiI
         x2(2) = yml(j1,j2,1)
         x2(3) = zml(j1,j2,1)
c
         dx(1) = (xysrf(1,i,1)-x2(1))*Pratio
         dx(2) = xysrf(2,i,1)-x2(2)
         dx(3) = xysrf(3,i,1)-x2(3)
c
         mxs = mxx(ifc2)
         joff = (j1-jopp)/(mxs-1)
         do 600 ix = 2,mxs
            j = jopp + joff*(ix-1)
            zeta = 0.5*(zgml(ix) + 1.0)
            xml(j,j2,1) = xml(j,j2,1)+dx(1)*zeta
            yml(j,j2,1) = yml(j,j2,1)+dx(2)*zeta
            zml(j,j2,1) = zml(j,j2,1)+dx(3)*zeta
  600    continue
  700 continue
C
      return
      end
c-----------------------------------------------------------------------
      subroutine hilites(ie,iside)
C
C     High-light just the side of the element
C
      include 'basics.inc'
C
C     Flash element side
C     Hilight side with Blinking lights
      CALL COLOR(2)
      IF(IFBWGKS)CALL COLOR(0)
      IF(ISIDE.EQ.5.OR.ISIDE.EQ.6) then
C        Flash mesh sides
         IF(ISIDE.EQ.5)THEN
            IED1=1
            IED2=4
         ELSE IF(ISIDE.EQ.6)THEN
            IED1=1
            IED2=4
         ENDIF
         call move(x(ied1,ie),y(ied1,ie))
         DO 220 IEDGE=IED1,IED2
            CALL DRAWED(IE,IEDGE,1)
220      continue
      else
c        SIDES 1-4  (SAME AS EDGES IN THIS CASE)
c        write(6,*) 'ie,is:',ie,iside
c        call move(x(iside,ie),y(iside,ie))
         call drawed(ie,iside,1)
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine shift
      include 'basics.inc'

      real x0(3),p1(3),p2(3),p3(3)

    1 continue           !    Set up menu and query

      nchoic = 0
      nchoic = nchoic+1
      item(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      item(nchoic)       =             'Redraw mesh'
      nchoic = nchoic+1
      item(nchoic)       =             'Center'
      nchoic = nchoic+1
      item(nchoic)       =             'Shift X'
      nchoic = nchoic+1
      item(nchoic)       =             'Shift Y'
      if (if3d) then
         nchoic = nchoic+1
         item(nchoic)    =             'Shift Z'
         nchoic = nchoic+1
         item(nchoic)    =             'Project on Plane'
      endif

      call menu(xmouse,ymouse,button,'NOCOVER') ! Prompt for user input

      n  = nel*(2**ndim)
      xmean = glsum(x,8*nel)/n
      ymean = glsum(y,8*nel)/n
      zmean = glsum(z,8*nel)/n
      if (if3d) then
         call prsrrr('XYZ mean:$',xmean,ymean,zmean)
      else
         call prsrr ('XY  mean:$',xmean,ymean)
      endif

      if (choice.eq.'UP MENU') return
      if (choice.eq.'Redraw mesh') then

         call redraw_mesh

      elseif (choice.eq.'Center') then

         xsep = glmin(x,8*nelm) - 1.e5
         ysep = glmin(y,8*nelm) - 1.e5
         zsep = glmin(z,8*nelm) - 1.e5

         xshift = -xmean
         yshift = -ymean
         zshift = -zmean

         call shifter(xshift,xsep,'>',x,'X')
         call shifter(yshift,ysep,'>',y,'Y')
         if (if3d) call shifter(zshift,zsep,'>',z,'Z')

      elseif (choice.eq.'Shift X') then
         call prs(
     $   'Input X-location separating shifted section.$')
         call rer(Xsep)
c        call drawline(Xsep,Ymax,Xsep,Ymin)
         call prs(
     $   'Input "<" or ">" to indicate desired shift section.$')
         call prs('("=" implies abort.)$')
         call res(ans,1)
         if (ANS.eq.'=') return
         call prs('Input amount to shift in X-direction$')
         call rer(Xshift)
         call shifter(Xshift,Xsep,ANS,X,'X')
      elseif (choice.EQ.'Shift Y') then
         call prs(
     $   'Input Y-location separating shifted section.$')
         call rer(Ysep)
c        call drawline(Ysep,Xmax,Ysep,Xmin)
         call prs(
     $   'Input "<" or ">" to indicate desired shift section.$')
         call prs('("=" implies abort.)$')
         call res(ANS,1)
         if (ANS.eq.'=') return
         call prs('Input amount to shift in Y-direction$')
         call rer(Yshift)
         call shifter(Yshift,Ysep,ANS,Y,'Y')

      elseif (choice.EQ.'Shift Z') then
         call prs(
     $   'Input Z-location separating shifted section.$')
         call rer(Zsep)
         call prs(
     $   'Input "<" or ">" to indicate desired shift section.$')
         call prs('("=" implies abort.)$')
         call res(ANS,1)
         if (ANS.eq.'=') return
         call prs('Input amount to shift in Z-direction$')
         call rer(Zshift)
         call Shifter(Zshift,Zsep,ANS,Z,'Z')

      elseif (choice.EQ.'Project on Plane') then
         call prs('Enter plane Point 1 (x,y,z):$')
         call rerrr(p1(1),p1(2),p1(3))
         call prs('Enter plane Point 2 (x,y,z):$')
         call rerrr(p2(1),p2(2),p2(3))
         call prs('Enter plane Point 3 (x,y,z):$')
         call rerrr(p3(1),p3(2),p3(3))
         call prs('Enter point on side to be preserved (x,y,z):$')
         call rerrr(x0(1),x0(2),x0(3))
         call prs('Enter tolerance (e.g., 0):$')
         call rer(tol)
         call plane_project_mesh(x0,p1,p2,p3,1,nel,tol)
      endif
      goto 1
      end
c-----------------------------------------------------------------------
      subroutine shifter(Shift,Sep,DIR,pts,coord)
      include 'basics.inc'
      real pts(8,nelm)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
C
      nvts = 4
      if (IF3D) nvts=8
      nkshift = 0
C
      if (DIR.eq.'>') then
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
   10       continue
c           if ((IFG.AND.IFL).and.Shift.lt.0.0) then ! cmt out; 9/29/05
C
C              If an element straddles the Separator, we have to
C              ensure that the shift operation doesn't "invert" the
C              element
C
c              if (pmin+shift.le.pmax) then
c                 CALL PRS(
c    $           'Error:  Attempt to shrink element to zero length$')
c                 CALL PRS(' Smax     Smin    Shift $')
c                 CALL PRRR(pmax,pmin,Shift)
c                 CALL PRS('Aborting shift operation$')
c                 return
c              ENDIF
c           ENDIF
            DO 20 j=1,nvts
               if (pts(j,k).ge.sep) pts(j,k)=pts(j,k)+shift
   20       continue
            if (if3d.and.pmin.ge.sep) then
              do l=1,12
                if (ccurve(l,k).eq.'m') then ! midside node
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            elseif (pmin.ge.sep) then  ! 2D
              do l=1,4
                if (ccurve(l,k).eq.'m') then ! midside node
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                endif
              enddo
            endif
            nkshift=nkshift+1
  100   continue
      ELSE
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
  110       continue
c           if ((IFG.AND.IFL).and.Shift.gt.0.0) then
C
C              If an element straddles the Separator, we have to
C              ensure that the shift operation doesn't "invert" the
C              element
C
c              if (pmax+shift.le.pmin) then
c                 CALL PRS(
c    $           'Error:  Attempt to shrink element to zero length$')
c                 CALL PRS(' Smax     Smin    Shift $')
c                 CALL PRRR(pmax,pmin,Shift)
c                 CALL PRS('Aborting shift operation$')
c                 return
c              ENDIF
c           ENDIF
            DO 120 j=1,nvts
               if (pts(j,k).le.sep) pts(j,k)=pts(j,k)+shift
  120       continue
            if (if3d.and.pmax.le.sep) then
              do l=1,12
                if (ccurve(l,k).eq.'m') then ! midside node
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
              do l=1,6
                if (ccurve(l,k).eq.'s') then
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                  if (coord.eq.'Z') curve(3,l,k)=curve(3,l,k)+shift
                endif
              enddo
            elseif (pmax.le.sep) then  ! 2D
              do l=1,4
                if (ccurve(l,k).eq.'m') then ! midside node
                  if (coord.eq.'X') curve(1,l,k)=curve(1,l,k)+shift
                  if (coord.eq.'Y') curve(2,l,k)=curve(2,l,k)+shift
                endif
              enddo
            endif
            nkshift=nkshift+1
  200    continue
      ENDIF

      WRITE(S,500) nkshift
  500 FORMAT(' Shifter',i11,' elements.$')
      CALL PRS(S)

      return
      END
c-----------------------------------------------------------------------
      subroutine curcnt  !     Recount the number of curved sides
      include 'basics.inc'

      ncurve=0
      do 9001 ie=1,nel
      do 9001 iedge=1,12
         if (ccurve(iedge,ie).ne.' ') then
            ncurve=ncurve+1
C           write(6,*) 'curve:',ie,iedge,ccurve(iedge,ie)
         endif
 9001 continue

      return
      end
c-----------------------------------------------------------------------
      subroutine splite(je,ie,ilist,ipln,frac)
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
c
      integer icalld,ieb1,ieb2
      save    icalld,ieb1,ieb2
      data    icalld,ieb1,ieb2 /0,0,0/
C
      if (frac.le.0.0.or.frac.ge.1.0) RETURN
C
C     Begin with copying the essential information
C
      call copyel(ie,je)
C
C     Compute XYZ on high definition mesh
C
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
C
C     Are we taking an R, S or T cut?
C
      FRACS=2.0*FRAC-1.0
      if (IPLN.EQ.0) then
         LISTA=ABS(LIST(ILIST))
         CALL DECOD(IPLANE,IPLN,IIE,IDUM,LISTA,NX,6,NELM)
      ENDIF
C
      if (IPLN.LE.2) then
C        X-plane
         I0=-1
         RRL(1)=FRACS
         DO 200 JCRN=-1,1,2
            RRL(3)=FLOAT(JCRN)
            DO 100 ic=-1,1,2
               RRL(2)=FLOAT(ic)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
               I0=I0+2
               I1=I0+1
               x(inv(i1),ie)=xval
               y(inv(i1),ie)=yval
               z(inv(i1),ie)=zval
               x(inv(i0),je)=xval
               y(inv(i0),je)=yval
               z(inv(i0),je)=zval
  100       continue
  200    continue
C
C        Don't forget Curve side data on new faces!
         JFAC1=EFACE(1)
         JFAC2=EFACE(2)
C
      elseif (IPLN.LE.4) then
C        Y-plane
         RRL(2)=FRACS
         DO 210 JCRN=-1,1,2
            RRL(3)=FLOAT(JCRN)
            DO 110 ic=-1,1,2
               RRL(1)=FLOAT(ic)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
C              II,JJ = 1,2
               II=1+(ic+1)/2
               JJ=1+(JCRN+1)/2
C              I0=1,2,5,6; I1=3,4,7,8
               I0=II+4*(JJ-1)
               I1=I0+2
               x(inv(i1),ie)=xval
               y(inv(i1),ie)=yval
               z(inv(i1),ie)=zval
               x(inv(i0),je)=xval
               y(inv(i0),je)=yval
               z(inv(i0),je)=zval
C
  110       continue
  210    continue
         JFAC1=EFACE(3)
         JFAC2=EFACE(4)
C
      ELSE
C        Z-plane
         I0=0
         RRL(3)=FRACS
         DO 220 JCRN=-1,1,2
            RRL(2)=FLOAT(JCRN)
            DO 120 ic=-1,1,2
               RRL(1)=FLOAT(ic)
               CALL EVALSC(XVAL,XP,RRL,1)
               CALL EVALSC(YVAL,YP,RRL,0)
               CALL EVALSC(ZVAL,ZP,RRL,0)
               I0=I0+1
               I1=I0+4
               x(inv(i1),ie)=xval
               y(inv(i1),ie)=yval
               z(inv(i1),ie)=zval
               x(inv(i0),je)=xval
               y(inv(i0),je)=yval
               z(inv(i0),je)=zval
  120       continue
  220    continue
         JFAC1=EFACE(5)
         JFAC2=EFACE(6)
      ENDIF
C
C     Curve side fix up
C
      if (CCURVE(JFAC1,IE).EQ.'O') CCURVE(JFAC1,JE)=' '
      if (CCURVE(JFAC2,IE).EQ.'O') CCURVE(JFAC2,JE)=' '

      if (CCURVE(JFAC1,IE).EQ.'s'  .and.
     $    CCURVE(JFAC2,IE).EQ.'s') then
          R = 0.5 * ( CURVE(4,JFAC1,IE) + CURVE(4,JFAC2,IE) )
          CCURVE(JFAC1,JE)='s'
          CCURVE(JFAC2,IE)='s'
          CURVE(4,JFAC2,IE) = R
          CURVE(4,JFAC1,JE) = R
          CALL COPY(CURVE(1,JFAC1,JE),CURVE(1,JFAC1,JE),3)
          CALL COPY(CURVE(1,JFAC2,IE),CURVE(1,JFAC1,IE),3)
      elseif (CCURVE(JFAC1,IE).EQ.'s'  .or.
     $        CCURVE(JFAC2,IE).EQ.'s') then
          R1i = 0.0
          R2i = 0.0
          if (CCURVE(JFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,JFAC1,IE)
          if (CCURVE(JFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,JFAC2,IE)
          R  = 2.0/(R1i+R2i)
C
c         CCURVE(JFAC1,JE)='s'
c         CCURVE(JFAC2,IE)='s'
          CCURVE(JFAC1,JE)=' '
          CCURVE(JFAC2,IE)=' '
          CURVE(4,JFAC2,IE) = R
          CURVE(4,JFAC1,JE) = R
          CALL COPY(CURVE(1,JFAC1,JE),CURVE(1,JFAC1,JE),3)
          CALL COPY(CURVE(1,JFAC2,IE),CURVE(1,JFAC1,IE),3)
      ENDIF
C
C     Check for cylinder  (Note this scheme fails for convex-convex case!)
C                         pff 4/10/93
C
      if (IPLN.LE.4) then
         DO 300 Ilev=1,2
           if (CCURVE(JFAC1,IE).EQ.'C'.and.
     $         CCURVE(JFAC2,IE).EQ.'C') then
               R1 = CURVE(1,JFAC1,IE)
               R2 = CURVE(1,JFAC2,IE)
               R = 0.5*(R1-R2)
               CURVE(1,JFAC2,IE) = -R
               CURVE(1,JFAC1,JE) =  R
               CCURVE(JFAC1,JE)='C'
               CCURVE(JFAC2,IE)='C'
           elseif (CCURVE(JFAC1,IE).EQ.'C'.  or.
     $             CCURVE(JFAC2,IE).EQ.'C') then
               R1i = 0.0
               R2i = 0.0
               if (CCURVE(JFAC1,IE).EQ.'C') R1i = 1.0/CURVE(1,JFAC1,IE)
               if (CCURVE(JFAC2,IE).EQ.'C') R2i = 1.0/CURVE(1,JFAC2,IE)
               R  = 2.0/(R1i-R2i)
               CURVE(1,JFAC2,IE) = -R
               CURVE(1,JFAC1,JE) =  R
               CCURVE(JFAC1,JE)=' '
               CCURVE(JFAC2,IE)=' '
           ENDIF
           JFAC1=JFAC1+4
           JFAC2=JFAC2+4
  300    continue
      ENDIF
C
C     Delete periodic bcs and undo curved sides
C
      return
      END
c-----------------------------------------------------------------------
      subroutine qsplite(ie,liste)
C
C     This routine is a hack-copy of octsplite.     (pff 12/21/98)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
     $              ,XVAL(-1:1,-1:1,-1:1),YVAL(-1:1,-1:1,-1:1)
     $              ,ZVAL(-1:1,-1:1,-1:1)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
      DIMENSION LISTE(8)
C
C     Compute XYZ on high definition mesh
C
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)

      if (if3d) then
c
        do 160 it= -1,1,1 !     Evaluate x,y,z on 27 pt stencil
        do 160 is= -1,1,1
        do 160 ir= -1,1,1
           rrl(3)=float(it)
           rrl(2)=float(is)
           rrl(1)=float(ir)
           call evalsc(xval(ir,is,it),xp,rrl,1)
           call evalsc(yval(ir,is,it),yp,rrl,0)
           call evalsc(zval(ir,is,it),zp,rrl,0)
  160 continue
C
C     Eight quadrants
C
      ke = 0
c     DO 860 Iz=-1,0,1
      DO 860 Iz= -1,-1,1
        DO 840 Iy=-1,0,1
          DO 820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          ke = ke+1
          Je = LISTE(ke)
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
c         DO 760 Jz=0,1
          DO 760 Jz=0,2,2
            jz0 = Jz+Iz
            DO 760 Jy=0,1
              jy0 = Jy+Iy
              DO 760 Jx=0,1
                jx0 = Jx+Ix
                iv = iv+1
                x(inv(iv),je)=xval(jx0,jy0,jz0)
                y(inv(iv),je)=yval(jx0,jy0,jz0)
                z(inv(iv),je)=zval(jx0,jy0,jz0)
  760     continue
C
          DO 400 Ifce=1,2
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           if (Ifce.eq.1) then
              JFAC = EFACE(1-Ix)
           elseif (Ifce.eq.2) then
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           if (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') then
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           elseif (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') then
             R1i = 0.0
             R2i = 0.0
             if (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             if (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           if (Ifce.LE.2) then
              DO 300 Ilev=1,2
                if (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') then
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  if (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  if (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
                elseif (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') then
                  R1i = 0.0
                  R2i = 0.0
                  if (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  if (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
                ENDIF
                JFAC =JFAC +4
                IFAC1=IFAC1+4
                IFAC2=IFAC2+4
  300         continue
           ENDIF
  400     continue
  820     continue
  840   continue
  860 continue


      else  ! 2D

        DO 2160 Is= -1,1,1
        DO 2160 Ir= -1,1,1
           rrl(3)=0.0
           rrl(2)=float(is)
           rrl(1)=float(ir)
           CALL EVALSC(XVAL(ir,is,1),XP,RRL,1)
           CALL EVALSC(YVAL(ir,is,1),YP,RRL,0)
 2160 continue
C
C     Four quadrants
C
      ke = 0
        DO 2840 Iy=-1,0,1
          DO 2820 Ix=-1,0,1
C
C         Begin with copying the essential information
C
          ke = ke+1
          Je = LISTE(ke)
C
          CALL COPYEL(IE,JE)
C
C         Eight vertices
C
          iv = 0
            do 2760 jy=0,1
              jy0 = jy+iy
              do 2760 jx=0,1
                jx0 = jx+ix
                iv = iv+1
                x(inv(iv),je)=xval(jx0,jy0,1)
                y(inv(iv),je)=yval(jx0,jy0,1)
                z(inv(iv),je)=zval(jx0,jy0,1)
 2760     continue
C
          DO 2400 Ifce=1,2
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           if (Ifce.eq.1) then
              JFAC = EFACE(1-Ix)
           elseif (Ifce.eq.2) then
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C        Curve side fix up - only internal faces have to be changed
C
           if (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') then
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           elseif (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') then
             R1i = 0.0
             R2i = 0.0
             if (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             if (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C
C          Check for cylinder  
C          (Note this scheme fails for convex-convex case!)
C          pff 4/10/93
C
           if (Ifce.LE.2) then
              ilev = 1
              if (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') then
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  if (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  if (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
              elseif (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') then
                  R1i = 0.0
                  R2i = 0.0
                  if (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  if (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
              ENDIF
           ENDIF
 2400     continue
 2820     continue
 2840   continue
 2860 continue



      endif

      call copyel(je,ie) !     Copy last element generated back onto IE

      return
      end
c-----------------------------------------------------------------------
      subroutine shift2

c     This routine is like "shift" except that the shifted points
c     are moved a scaled distance, proportional to the distance from
c     the selected point.

      include 'basics.inc'

    1 continue

      nchoic = 0
      nchoic = nchoic+1
      ITEM(nchoic)       =             'UP MENU'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Redraw mesh'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch X'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch Y'
      if (IF3D) then
         nchoic = nchoic+1
         ITEM(nchoic)    =             'Stretch Z'
      ENDIF
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch R'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch theta'
c     nchoic = nchoic+1
c     ITEM(nchoic)       =             'Stretch inside circle'
      nchoic = nchoic+1
      ITEM(nchoic)       =             'Stretch outside circle'


      CALL MENU(XMOUSE,YMOUSE,BUTTON,'NOCOVER') ! Prompt for user input.


      if (choice.eq.'UP MENU') return

      if (choice.eq.'Redraw mesh') then
         call redraw_mesh
      elseif (choice.eq.'Stretch outside circle') then
         call stretch_outside_circ
      elseif (choice.eq.'Stretch theta') then
         call stretch_theta
      elseif (choice.eq.'Stretch R') then
         call stretch_rad
      elseif (choice.eq.'Stretch X') then
         CALL PRS('Input X-location separating shifted section.$')
         CALL RER(Xsep)
c        CALL DRAWLINE(Xsep,Ymax,Xsep,Ymin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         if (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in X-direction$')
         CALL RER(Xshift)
         CALL Shifter2(Xshift,Xsep,ANS,X,'X')
      elseif (choice.eq.'Stretch Y') then
         CALL PRS(
     $   'Input Y-location separating shifted section.$')
         CALL RER(Ysep)
c        CALL DRAWLINE(Ysep,Xmax,Ysep,Xmin)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         if (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in Y-direction$')
         CALL RER(Yshift)
         CALL Shifter2(Yshift,Ysep,ANS,Y,'Y')
      elseif (choice.eq.'Stretch Z') then
         CALL PRS(
     $   'Input Z-location separating shifted section.$')
         CALL RER(Zsep)
         CALL PRS(
     $   'Input "<" or ">" to indicate desired stretch section.$')
         CALL PRS('("=" implies abort.)$')
         CALL RES(ANS,1)
         if (ANS.eq.'=') return
         CALL PRS('Input amount to stretch in Z-direction$')
         CALL RER(Zshift)
         CALL Shifter2(Zshift,Zsep,ANS,Z,'Z')
      endif
      goto 1
      end
c-----------------------------------------------------------------------
      subroutine shifter2a(pmin,pmax,gain,Shift,Sep,DIR,pts,coord)
c
c     In this pass, we just figure out range of geometry in shifted 
c     section
c
      include 'basics.inc'
      real pts(8,nelm)
      CHARACTER*1 DIR,coord
      LOGICAL IFG,IFL
C
      CALL PRS('Input Density Gain, G :$')
      CALL RER(gain)
      if (gain.eq.0) gain=1.
      if (gain.lt.0) gain=-1./gain
c
      nvts = 4
      if (IF3D) nvts=8
      nkshift = 0
C
      pmin= 9.99e15
      pmax=-9.99e15
c
      DO K=1,NEL
      DO j=1,nvts
         pmin=min(pmin,pts(j,k))
         pmax=max(pmax,pts(j,k))
      enddo
      enddo
c
      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2b(Shift,Sep,DIR,pts,coord)
c
c     Standard linear stretch -- unity gain.
c
      include 'basics.inc'
      real pts(8,nelm)
      character*1 dir,coord
      logical ifg,ifl

      nvts = 4
      if (if3d) nvts=8
      nkshift = 0

      if (coord.eq.'X') i=1
      if (coord.eq.'Y') i=2
      if (coord.eq.'Z') i=3
      nedge = 4 + 8*(ndim-2)

      if (DIR.eq.'>') then
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
   10       continue
            do 20 j=1,nvts
               if (pts(j,k).ge.sep) pts(j,k)=shift*(pts(j,k)-sep) + sep
   20       continue

            if (pmin.ge.sep) then
              do l=1,nedge
                if (ccurve(l,k).eq.'m') then ! midside node
                  curve(i,l,k)=shift*(curve(i,l,k)-sep)+sep
                elseif (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
                endif
              enddo
            endif
            nkshift=nkshift+1
  100   continue
      ELSE
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
  110       continue
            DO 120 j=1,nvts
c              if (pts(j,k).le.sep) pts(j,k)=shift*(pts(j,k)-sep) + sep
               if (pts(j,k).le.sep) then
                   pold = pts(j,k)
                   pts(j,k)=shift*(pts(j,k)-sep) + sep
                   pnew = pts(j,k)
                   write(6,*) 'stretch',pold,pnew,sep,shift,j,dir,coord
               endif
  120       continue

            if (pmax.le.sep) then
              do l=1,nedge
                if (ccurve(l,k).eq.'m') then ! midside node
                  curve(i,l,k)=shift*(curve(i,l,k)-sep)+sep
                elseif (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
                endif
              enddo
            endif
            nkshift=nkshift+1
  200    continue
      ENDIF

      WRITE(S,500) nkshift
  500 FORMAT(' Shifter2b',i11,' elements.$')
      CALL PRS(S)

      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2c(Shift,Sep,DIR,pts,coord,gain,qmax,qmin)

c     Geometric gain

      include 'basics.inc'
      real pts(8,nelm)
      character*1 dir,coord
      logical ifg,ifl

      nvts = 4
      if (IF3D) nvts=8
      nkshift = 0
c
c     We only get here if gain > 0,  different than 1
c
      g = log(gain)
      c = shift*g/(gain-1.)

      nedge = 4 + 8*(ndim-2)
      if (coord.eq.'X') i=1
      if (coord.eq.'Y') i=2
      if (coord.eq.'Z') i=3


      if (DIR.eq.'>') then
         dx  = qmax - sep
         gdx = g/dx
         DO 100 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 10 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
   10       continue
            DO 20 j=1,nvts
               if (pts(j,k).ge.sep) then
                  argu = gdx*(pts(j,k)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  write(6,1) shift,gain,g,c,dx,gdx,pts(j,k),xn
    1             format('s:',1p8e12.3)
                  pts(j,k)=xn
               endif
   20       continue
            if (pmin.ge.sep) then
              do l=1,nedge
                if (ccurve(l,k).eq.'m') then ! midside node
                  argu = gdx*(curve(i,l,k)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  curve(i,l,k)=xn
                elseif (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
                endif
              enddo
            endif
            nkshift=nkshift+1
  100   continue
      ELSE
         dx  = qmin - sep
         gdx = g/dx
         DO 200 K=1,NEL
            IFG=.FALSE.
            IFL=.FALSE.
            pmin= 9.99e15
            pmax=-9.99e15
            DO 110 j=1,nvts
               if (pts(j,k).ge.sep) then
                  IFG=.TRUE.
                  pmin=min(pmin,pts(j,k))
               ELSE
                  IFL=.TRUE.
                  pmax=max(pmax,pts(j,k))
               ENDIF
  110       continue
            DO 120 j=1,nvts
               if (pts(j,k).le.sep) then
                  argu = gdx*(pts(j,k)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  pts(j,k)=xn
               endif
  120       continue
            if (pmax.le.sep) then
              do l=1,nedge
                if (ccurve(l,k).eq.'m') then ! midside node
                  argu = gdx*(curve(i,l,k)-sep)
                  xn   = sep + c*dx*(exp(argu)-1.)/g
                  curve(i,l,k)=xn
                elseif (ccurve(l,k).eq.'s') then
                  call prs('need to fix spherical+shifter2, pff$')
                endif
              enddo
            endif
            nkshift=nkshift+1
  200    continue
      ENDIF

      write(s,500) nkshift
  500 format(' Shifter2c',i11,' elements.$')
      call prs(s)

      return
      END
c-----------------------------------------------------------------------
      subroutine shifter2(Shift,Sep,DIR,pts,coord)
c
c     Standard linear stretch -- unity gain.
c
      include 'basics.inc'
      real pts(8,nelm)
      character*1 dir,coord
      logical ifg,ifl
c
      call shifter2a(qmin,qmax,gain,Shift,Sep,DIR,pts,coord)
      write(6,*) 'qmnx gn:',qmin,qmax,gain
c
      if (abs(gain-1.) .gt. 1.e-5) then
         call shifter2c(Shift,Sep,DIR,pts,coord,gain,qmax,qmin)
      else
         call shifter2b(Shift,Sep,DIR,pts,coord)
      endif
c
      return
      END
c-----------------------------------------------------------------------
      subroutine hsplite(ie,liste)
C
C     This routine is a hack-copy of octsplite to generate hex
c     decompositions of the square     (pff 8/10/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
      integer inv(8)
      save    inv
      data    inv /1,2,4,3,5,6,8,7/
      dimension liste(8)
c
c     Evaluate x,y,z on new stencil
c
c
c     3---a-------4
c     |   |     / |
c     a---a   /   |
c     |    \0\    | 
c     |   /   a---+
c     | /     |   |
c     1-------+---2
c
c
      integer qmap(4,6) ! quad map
      save    qmap
      data    qmap /   1,  2,  6,  7   ! These orderings will
     $             ,   2,  3,  7,  8   ! preserve the existing BCs
     $             ,   7,  8,  6, 11   !
     $             ,   1,  6,  4,  5   ! Except for periodic, of course
     $             ,   4,  5,  9, 10   !
     $             ,   5,  6, 10, 11 /
c
c     9---0-------1
c     |   |     / |
c     4---5   /   |
c     |    \6\    | 
c     |   /   7---8
c     | /     |   |
c     1-------2---3
c
c
      parameter (aaa = .20)
      parameter (aam = -aaa)
      real rs(2,11)
      save rs
      data rs /  -1.0 , -1.0   !  1
     $        ,   aaa , -1.0   !  2
     $        ,   1.0 , -1.0   !  3
     $        ,  -1.0 ,  aaa   !  4
     $        ,   aam ,  aaa   !  5
     $        ,   0.0 ,  0.0   !  6
     $        ,   aaa ,  aam   !  7
     $        ,   1.0 ,  aam   !  8
     $        ,  -1.0 ,  1.0   !  9
     $        ,   aam ,  1.0   ! 10
     $        ,   1.0 ,  1.0 / ! 11
      real xval(3,11,2)
      integer e,v,q
c
c     Compute XYZ on high definition mesh
c
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
c
      iz1 = 1
      if (if3d) iz1 = 2
c
      do k=1,iz1
      do v=1,11
         rrl(1) = rs(1,v)
         rrl(2) = rs(2,v)
         rrl(3) = 2*k-3
         call evalsc(xval(1,v,k),xp,rrl,1)
         call evalsc(xval(2,v,k),yp,rrl,0)
         call evalsc(xval(3,v,k),zp,rrl,0)
      enddo
      enddo
c
C
C     six elements
C
      ke = 0
      do e=1,6
C
C        Begin by copying the essential information
C
         ke = ke+1
         Je = LISTE(ke)
         call copyel(ie,je)
C
C        Eight vertices
C
         v = 0
         do kv=1,2
         do jv=1,4
            q = qmap(jv,e)
            v = v+1
            x(inv(v),je)=xval(1,q,kv)
            y(inv(v),je)=xval(2,q,kv)
            z(inv(v),je)=xval(3,q,kv)
         enddo
         enddo
C
      enddo
C
C     Copy last element generated back onto IE
C
      call copyel(Je,Ie)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine zsplite(ie,liste)
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      parameter (nxm3=nxm*nym*nzm)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),rrl(3)
      integer inv(8)
      save    inv
      data    inv /1,2,4,3,5,6,8,7/
      dimension liste(8)
c
      real xval(3,12)
      integer e,v,q
c
c     Compute XYZ on high definition mesh
c
      call prsiii('This is nxyz:$',nx,ny,nz)
      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)
c
      l=0
      do k=-1,1
      do j=-1,1,2
      do i=-1,1,2
         l = l+1
         rrl(1) = i
         rrl(2) = j
         rrl(3) = k
         call evalsc(xval(1,l),xp,rrl,1)
         call evalsc(xval(2,l),yp,rrl,0)
         call evalsc(xval(3,l),zp,rrl,0)
      enddo
      enddo
      enddo

      ke = 0
      kv = 0
      do e=1,2 ! two elements

         ke = ke+1
         je = liste(ke)
         call copyel(ie,je) ! copy essential information

         do v=1,8                       !  Eight vertices
            x(inv(v),je)=xval(1,v+kv)
            y(inv(v),je)=xval(2,v+kv)
            z(inv(v),je)=xval(3,v+kv)
         enddo
         kv = 4
      enddo

      call copyel(Je,Ie) ! Copy last element generated back onto IE

      return
      end
c-----------------------------------------------------------------------
      function ie_click(prompt)  ! element that is clicked upon
      include 'basics.inc'
      character*80 prompt
C
1     call prs(prompt)
      call mouse(xmouse,ymouse,button)
c     if(xscr(xmouse).gt.1.0 .and. yscr(ymouse).gt.0.62) then
      if(xscr(xmouse).gt.1.0) then   ! apparently trying to use keypad
         call prs('Choosing all elements? (y/n)$')
         call res(ans,1)
         ie_click = 0
         if (ans.eq.'y'.or.ans.eq.'Y') then
            ie_click=-1
            call prs('All elements selected.$')
            return
         elseif (ans.eq.'n'.or.ans.eq.'N') then
            ie_click=0
            call prs('No elements selected.$')
            return
         else
            goto 1
         endif
      endif

      rmin=1.0e10
      do iel=1,nel
         rad=sqrt( (xmouse-xcen(iel))**2 + (ymouse-ycen(iel))**2 )
         if (if3d.and.rad.lt.rmin .and. numapt(iel).eq.ilevel) then
            rmin=rad
            ie_click=iel
         elseif (rad.lt.rmin)then
            rmin=rad
            ie_click=iel
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine msplit
C
C     This routine is a hack-copy of octsplite to generate 
c     multi-element decompositions of the square     (pff 9/28/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
      include 'basics.inc'
c
      ie = ie_click('Click on element to refine:$')
      if (ie.eq.0) return

      if (if3d) then

         call prs('Enter number of partitions in r,s,t (>0):$')
         call reiii(nxsp,nysp,nzsp)

         call prs  ('Enter ratios for r,s,t (1=uniform):$')
         call rerrr(rax,ray,raz)

      else
         call prs('Enter number of partitions in r,s (>0):$')
         call reii(nxsp,nysp)

         call prs('Enter ratios for r,s (1=uniform):$')
         call rerr(rax,ray)

         nzsp    = 1
         raz     = 1.

      endif

      write(6,*)
      write(6,9) ie,(x(jj,ie),jj=1,4),' xorg'
      write(6,9) ie,(y(jj,ie),jj=1,4),' yorg'
    9 format(i5,4f14.5,a5)

      if (ie.gt.0) then
         call msplite(ie,nxsp,nysp,nzsp,rax,ray,raz)
      else
         nelo = nel
         do ie=1,nelo
            call msplite(ie,nxsp,nysp,nzsp,rax,ray,raz)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_ratios(rrz,nzspi,razi)

      real rrz(1),razi(0:1)

      nzsp = abs(nzspi)
      if (nzspi.lt.0) then ! custom ratio
         rdz = razi(nzsp)-razi(0)
         if (rdz.gt.0) rdz = 2./rdz
         do k=1,nzsp+1
            rrz(k) = -1 + rdz*(razi(k)-razi(0))
         enddo
      else
         do k=1,nzsp+1
            rrz(k) = get_ratio(k,nzsp,razi)
         enddo
      endif
c     call outmat(rrz,1,nzsp+1,'rrz  ',nzspi)
 
      return
      end
c-----------------------------------------------------------------------
      function get_ratio(i,n,ratio)
c
c     Compute x_i \in [-1,1], where
c
c         dx_i := x_i+1 - x_i,  i=1,...,n
c
c     forms a geometric progression defined by ratio
c
c
      if (i.eq.1) then

         get_ratio = -1.

      else

         d=1
         s=0
         do j=2,n+1
            s = s+d
            d = d*ratio
         enddo

         d = 2/s   ! rescale initial value of d
         s = 0

         do j=2,i
            s = s+d
            d = d*ratio
         enddo

         get_ratio = s-1.

      endif

c     write(6,*) i,n,ratio,get_ratio

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_curve(e,i,j,k,nxsp,nysp,nzsp,rrx,rry,rrz,xl,yl,zl)

c     Assign / repair curve-side info for edges

      include 'basics.inc'
      real rrx(1),rry(1),rrz(1)
      real xl(3,3,3),yl(3,3,3),zl(3,3,3)

      integer e

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation

      character m12(12),cct(12)
      save      m12
      data      m12  / 12 * 'm' /


      call fix_c_curve(e,i,j,k,nxsp,nysp,nzsp,rrx,rry,rrz)
      call fix_s_curve(e,i,j,k,nxsp,nysp,nzsp,rrx,rry,rrz)


      call chcopy(cct,m12,12)

      if (if3d) then                      !  NOTE THE SERIOUS CONFLICT:
         if (ccurve(6,e).eq.'s') then     !
            cct(6) = 's'                  !  cct(5) = 's'
            call blank(cct(5),4)          !
         endif                            !  and 
         if (ccurve(5,e).eq.'s') then     !
            cct(5) = 's'                  !  cct(5) = 'm' are compatible, but
            call blank(cct,4)             !  not allowed because of naming conv!
         endif
      endif

      nedge = 4 + 8*(ndim-2)
      do kk=1,nedge
         if (ccurve(kk,e).eq.'C') then
c           do nothing
         elseif (cct(kk).eq.'m') then
            ccurve(kk,e)='m'
            jj = eindx(kk)
            curve(1,kk,e) = xl(jj,1,1)
            curve(2,kk,e) = yl(jj,1,1)
            curve(3,kk,e) = zl(jj,1,1)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_c_curve(e,ii,jj,kk,nxsp,nysp,nzsp,rrx,rry,rrz)
      include 'basics.inc'
      real rrx(1),rry(1),rrz(1)

      integer c_pair(2,4)
      save    c_pair
      data    c_pair / 4,2 , 8,6 , 1,3 , 5,7 /

      integer e
      character*1 c0,c1

      nskip=2
      if (if3d) nskip=1
c
c     Check "r" curves (edges 4 and 2)
c
         
      do i=1,4,nskip

         i0 = c_pair(1,i)
         i1 = c_pair(2,i)
         c0 = ccurve(i0,e)
         c1 = ccurve(i1,e)

         if (i.le.2) then
            r0 = rrx(ii  )
            r1 = rrx(ii+1)
            mm = ii
            nn = nxsp
         else
            r0 = rry(jj  )
            r1 = rry(jj+1)
            mm = jj
            nn = nysp
         endif
         r0 = 0.5*(r0+1.)
         r1 = 0.5*(r1+1.)

c        write(6,1) r0,ccurve(i0,e),curve(1,i0,e),' r0 b4'
c        write(6,1) r1,ccurve(i1,e),curve(1,i1,e),' r1 b4'

         if (c0.eq.'C' .and. c1.eq.'C') then
            curve0 = abs(curve(1,i0,e))
            curve1 = abs(curve(1,i1,e))
            cnew_0 = (1.-r0)*curve0 + r0*curve1
            cnew_1 = (1.-r1)*curve0 + r1*curve1
            if (curve(1,i0,e).lt.0) cnew_0 = -cnew_0
            if (curve(1,i1,e).lt.0) cnew_1 = -cnew_1
            curve(1,i0,e) = cnew_0
            curve(1,i1,e) = cnew_1
         elseif (c0.eq.'C') then
            c_old         = curve(1,i0,e)
            curve(1,i0,e) = c_old/(1.-r0)
            if (mm.lt.nn) then
               ccurve(i1,e) = 'C'
               curve(1,i1,e) = -c_old/(1.-r1)
            endif
         elseif (c1.eq.'C') then
            c_old         = curve(1,i1,e)
            curve(1,i1,e) = c_old/r1
            if (mm.gt.1) then
               ccurve(i0,e) = 'C'
               curve(1,i0,e) = -c_old/r0
            endif
         endif

c        write(6,1) r0,ccurve(i0,e),curve(1,i0,e),' r0 af'
c        write(6,1) r1,ccurve(i1,e),curve(1,i1,e),' r1 af'
c  1     format(f12.4,2x,a1,2x,f12.4,a6)

      enddo

c     write(6,*) 'continue?'
c     read (5,*) c0

      return
      end
c-----------------------------------------------------------------------
      subroutine fix_s_curve(e,ii,jj,kk,nxsp,nysp,nzsp,rrx,rry,rrz)
      include 'basics.inc'

      real rrx(1),rry(1),rrz(1)

      integer e
      character*1 c0,c1

c     Check "s" curves faces 5 and 6

      i0 = 5
      i1 = 6
      c0 = ccurve(i0,e)
      c1 = ccurve(i1,e)

      r0 = rrz(kk  )
      r1 = rrz(kk+1)
      r0 = 0.5*(r0+1.) ! ratio on [0,1]
      r1 = 0.5*(r1+1.)

c     write(6,1) r0,curve(1,i0,e),ccurve(i0,e),' r0 b4'
c     write(6,1) r1,curve(1,i1,e),ccurve(i1,e),' r1 b4'
   1  format(1p2e13.5,1x,a1,1x,a6)

      if (c0.eq.'s' .and. c1.eq.'s') then
         do i=1,4  ! average all 4 quantities:  radius, ctr_x, ctr_y, ctr_z
            curve0 = abs(curve(i,i0,e))
            curve1 = abs(curve(i,i1,e))
            cnew_0 = (1.-r0)*curve0 + r0*curve1
            cnew_1 = (1.-r1)*curve0 + r1*curve1
            curve(i,i0,e) = cnew_0
            curve(i,i1,e) = cnew_1
         enddo
      elseif (c0.eq.'s') then

         if (kk.gt.1) ccurve(i0,e) = ' '  ! turn off spherical projection

c        c_old         = curve(1,i0,e)
c        curve(1,i0,e) = c_old/(1.-r0)
c        if (mm.lt.nn) then
c           ccurve(i1,e) = 's'
c           curve(1,i1,e) = -c_old/(1.-r1)
c        endif

      elseif (c1.eq.'s') then

         if (kk.lt.nzsp) ccurve(i1,e) = ' '  ! turn off spherical projection

c        c_old         = curve(1,i1,e)
c        curve(1,i1,e) = c_old/r1
c        if (mm.gt.1) then
c           ccurve(i0,e) = 's'
c           curve(1,i0,e) = -c_old/r0
c        endif

      endif

c     Check "a" curve on edge 9

      i0 = 9
      c0 = ccurve(i0,e)
      if ((ii.gt.1.or.jj.gt.1).and.c0.eq.'a') ccurve(i0,e) = ' '

      return
      end
c-----------------------------------------------------------------------
      subroutine arcsph_e(xml,yml,zml,nxl,nyl,nzl,e,isid)
      include 'basics.inc'
      real xml(nxl,nyl,nzl),yml(nxl,nyl,nzl),zml(nxl,nyl,nzl)
      integer e

c     ctmpg is used in this format in several subsequent routines
      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)

      real xcrv(3,lx1),osph(3),v1(3),v2(3)

c     Edg2v stores start point for each edge, in hypercube notation
      integer edg2v(12)  ! vertices in symm notation, edge in prep not.
      save    edg2v
      data    edg2v / 0,1,2,0, 4,5,6,4, 0,1,3,4 /

      integer eskip(12)  ! vertices in symm notation, edge in prep not.
      save    eskip
      data    eskip / 1,2,1,2, 1,2,1,2, 3,3,3,3 /

      integer eindx(8)
      save    eindx
      data    eindx  / 1, 2, 4, 3, 5, 6, 8, 7 /

      integer estride(3),estart(8)

      estride(1) = 1
      estride(2) = nxl
      estride(3) = nxl*nyl

      l=0
      do k=0,1
      do j=0,1
      do i=0,1
         l=l+1
         estart(l) = i + j*nxl + k*nxl*nyl
      enddo
      enddo
      enddo

      istart = edg2v  (isid)
      istart = estart (istart)
      istrid = estride(isid)

      i=istart
      do k=1,nxl   ! this inveighs against nzl^=nxl
         i=i+istrid
         xcrv(1,k) = xml(i,1,1)
         xcrv(2,k) = yml(i,1,1)
         xcrv(3,k) = zml(i,1,1)
c        write(79,1) i,xcrv(1,i),xcrv(2,i),xcrv(3,i)
      enddo

      osph(1) = curve(1,isid,e)
      osph(2) = curve(2,isid,e)
      osph(3) = curve(3,isid,e)
      radius  = curve(4,isid,e)

      call sub3  (v1,xcrv(1,1  ),osph,3)
      call sub3  (v2,xcrv(1,nxl),osph,3)
      call edg3d (xcrv,v1,v2,1,nxl,1,1,nxl,1)

      call cmult (xcrv,radius,3*nxl)

      do i=1,nxl       !     Add sphere center offset
         xcrv(1,i)=xcrv(1,i)+osph(1)
         xcrv(2,i)=xcrv(2,i)+osph(2)
         xcrv(3,i)=xcrv(3,i)+osph(3)
c        write(78,1) i,xcrv(1,i),xcrv(2,i),xcrv(3,i)
c  1     format(i4,1p3e13.5)
      enddo

      i=istart
      do k=1,nxl   ! map back to edge
         i=i+istrid
         xml(i,1,1) = xcrv(1,k)
         yml(i,1,1) = xcrv(2,k)
         zml(i,1,1) = xcrv(3,k)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine el_swap_ip(p,n,nmax)
      integer p(1),xstart

      xstart = nmax
c
c     In-place permutation: x' = x(p)
c
      do k=1,n
         if (p(k).gt.0) then   ! not swapped
c           xstart     = x(k)
            call copyel(k,xstart)  ! k --> xstart
            loop_start = k
            last       = k
            do j=k,n
               next    = p(last)
               if (next.lt.0) then
                  write(6,*) 'Hey! swap_ip problem.',j,k,n,next
                  return
c                 call exitt
               elseif (next.eq.loop_start) then
                  call copyel(xstart,last)  ! xstart --> last
c                 x(last) = xstart
                  p(last) = -p(last)
                  goto 10
               else
                  call copyel(next,last)  ! next --> last
c                 x(last) = x(next)
                  p(last) = -p(last)
                  last    = next
               endif
            enddo
   10       continue
         endif
      enddo
c
      do k=1,n
         p(k) = -p(k)
      enddo

      call el_bc_swap_ip(p,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine el_bc_swap_ip(p,n)
      integer p(1),xstart
      include 'basics.inc'
      common /wierd/ list1(nelm),list2(nelm)
      integer e,f,e_old,e_new
      character*3 cb

      nfaces = 2*ndim

      call icopy    (list1,p,n)
      call iswap_ip (list1,p,n)

      do ifld=1,nflds
      do e=1,n
      do f=1,nfaces
         cb=cbc(f,e,ifld)
         if (cb.eq.'P  ') then
            e_old = bc(1,f,e,ifld) ! old element number
            e_old = ibc(f,e,ifld)  ! old element number
            e_new = list1(e_old)
            bc(1,f,e,ifld) = e_new
            ibc(f,e,ifld)  = e_new
         endif
      enddo
      enddo
      enddo
         
      return
      end
c-----------------------------------------------------------------------
      subroutine renum_special
C
C     Renumber elements 
C
      include 'basics.inc'
      common /cenewi/ enew(nelm),ind(nelm)
      integer enew,e,v

      nsolid = 0

      call refresh
      call drmenu('NOCOVER')
      call drgrid


      eps = -.01
      do e=1,nel
         enew(e) = e
         dmin    = dist_solid(e)
         if (dmin.lt.eps) then
            enew(e) = e + 2*nel ! solid element
            nsolid = nsolid+1
            call drawel(e)
            if (mod(nsolid,10).eq.0) then
               call prsis('Number found =$',nsolid,'. Continue?$')
c              call res(ans,1)
            endif
         endif
      enddo
      call prsis('Found $',nsolid,' solid elements.$')

      call irank     (enew,ind,nel)
      call el_swap_ip(ind ,nel,nelm-1)

      return
      end
c-----------------------------------------------------------------------
      function dist_solid(e) ! distance of element from solid

      include 'basics.inc'
      integer e,v

      nv = 2**ndim

      dmin = 1.e22

      do v=1,nv
         dist = dist_special(x(v,e),y(v,e),z(v,e))
         dmin = min(dmin,dist)
      enddo

      dist_solid = dmin

      return
      end
c-----------------------------------------------------------------------
      function dist_special(x,y,z)

      integer ncyl
      save    ncyl
      data    ncyl /0/

      common /cyldata/ cyl(3,100)

      if (ncyl.eq.0) then
         open(unit=88,file='cyl.in',status='old',err=20)
         read(88,*) ncyl
         do i=1,ncyl
            read(88,*) cyl(1,i),cyl(2,i),cyl(3,i)
c           write(6,*) cyl(1,i),cyl(2,i),cyl(3,i),' cyl ',i
         enddo
         close(88)
         goto 40

   20    continue
         call prs('Input number of cylinders (e.g., 2):$')
         call rei(ncyl)
         do i=1,ncyl
            call prsis('x0,y0,rad for cyl $',i,' ?$')
            call rerrr(cyl(1,i),cyl(2,i),cyl(3,i))
         enddo

   40    continue
      endif
      dmin = 1.e22
      do i=1,ncyl
c        write(6,*) cyl(1,i),cyl(2,i),cyl(3,i),' cylb ',i
         x0  = cyl(1,i)
         y0  = cyl(2,i)
         rad = cyl(3,i)
         dist = dcircle(x,y,x0,y0,rad) ! negative --> inside circle
         dmin = min(dmin,dist)
c        write(6,*) cyl(1,i),cyl(2,i),cyl(3,i),' cyld ',dmin
      enddo
      dist_special = dmin

      return
      end
c-----------------------------------------------------------------------
      function dcircle(x,y,xc,yc,r)

      d=(x-xc)**2+(y-yc)**2
      if (d.gt.0) d=sqrt(d)
      dcircle = d-r

      return
      end
c-----------------------------------------------------------------------
      subroutine map_m_to_n(a,na,b,nb,if3d,w,ldw)
c
c     Input:   b
c     Output:  a
c
      real a(1),b(1),w(1)
      logical if3d
c
      parameter(lx=50)
      real za(lx),zb(lx)
c
      real iba(lx*lx),ibat(lx*lx)
      save iba,ibat
c
      integer nao,nbo
      save    nao,nbo
      data    nao,nbo  / -9, -9/
c
      if (na.gt.lx.or.nb.gt.lx) then
         write(6,*)'ERROR: increase lx in map_m_to_n to max:',na,nb
         call exitt
      endif
c
      if (na.ne.nao  .or.   nb.ne.nbo) then
         nao = na
         nbo = nb
         call zwgll(za,w,na)
         call zwgll(zb,w,nb)
         call igllm(iba,ibat,zb,za,nb,na,nb,na)
      endif

      call specmpn(a,na,b,nb,iba,ibat,if3d,w,ldw)

      return
      end
c
c-----------------------------------------------------------------------
      subroutine specmpn(b,nb,a,na,ba,ab,if3d,w,ldw)

c     -  Spectral interpolation from A to B via tensor products
c     -  scratch arrays: w(na*na*nb + nb*nb*na)

c     5/3/00  -- this routine replaces specmp in navier1.f, which
c                has a potential memory problem



      real b(nb,nb,nb),a(na,na,na),ba(1),ab(1),w(ldw)
      logical if3d

      ltest = na*nb
      if (if3d) ltest = na*na*nb + nb*na*na
      if (ldw.lt.ltest) then
         write(6,*) 'ERROR specmp:',ldw,ltest,if3d
         call exitt
      endif

      if (if3d) then
         nab = na*nb
         nbb = nb*nb
         call mxm(ba,nb,a,na,w,na*na)
         k=1
         l=na*na*nb + 1
         do iz=1,na
            call mxm(w(k),nb,ab,na,w(l),nb)
            k=k+nab
            l=l+nbb
         enddo
         l=na*na*nb + 1
         call mxm(w(l),nbb,ab,na,b,nb)
      else
         call mxm(ba,nb,a,na,w,na)
         call mxm(w,nb,ab,na,b,nb)
      endif
      return
      end
c
c-----------------------------------------------------------------------
      subroutine tensr3(v,n1,n2,n3,u,m1,m2,m3,A,Bt,Ct,w,ldw,if3d)
C
C     -  Tensor product application of v = (C x B x A) u
C        NOTE -- the transpose of B & C must be input, rather than B & C.
C
C     -  scratch arrays: w(ldw)

c     Nominal declarations:
c        real v(n1,n2,n3),u(m1,m2,m3)
c        real a(n1,m1),bt(m2,n2),ct(m3,n3)

      real v(1),u(1),a(1),bt(1),ct(1),w(ldw)

      logical if3d

      ndw = n1*m2*m3 + n1*n2*m3
      if (ndw.gt.ldw) then
         write(6,*) ndw,ldw,' ERROR in tensr3. Memory problem.'
         write(6,*) n1,n2,n3
         write(6,*) m1,m2,m3,if3d
         call exitt
      endif

      if (if3d) then
         call mxm(a,n1,u,m1,w,m2*m3)
         i=1
         j=1 + n1*m2*m3
         do k=1,m3
            call mxm(w(i),n1,bt,m2,w(j),n2)
            i=i+n1*m2
            j=j+n1*n2
         enddo
         call mxm(w(i),n1*n2,ct,m3,v,n3)
      else
         call mxm(a,n1,u ,m1,w,m2)
         call mxm(w,n1,bt,m2,v,n2)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine gen_int_gz(j,jt,g,n,z,m)

c     Generate interpolater from m z points to n g points

c        j   = interpolation matrix, mapping from z to g
c        jt  = transpose of interpolation matrix
c        m   = number of points on z grid
c        n   = number of points on g grid

      real j(n,m),jt(m,n),g(n),z(m)

c     call outmat(z,1,m,'zptsA',n)


      mpoly  = m-1
      do i=1,n
         call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
      enddo

      call transpose_r(j,n,jt,m)

      return
      end
c-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
c
c     This routine evaluates the derivative based on all points
c     in the stencils.  It is more memory efficient than "fd_weights"
c
c     This set of routines comes from the appendix of 
c     A Practical Guide to Pseudospectral Methods, B. Fornberg
c     Cambridge Univ. Press, 1996.   (pff)
c
c     Input parameters:
c       xx -- point at wich the approximations are to be accurate
c       x  -- array of x-ordinates:   x(0:n)
c       n  -- polynomial degree of interpolant (# of points := n+1)
c       m  -- highest order of derivative to be approxxmated at xi
c
c     Output:
c       c  -- set of coefficients c(0:n,0:m).
c             c(j,k) is to be applied at x(j) when
c             the kth derivative is approxxmated by a 
c             stencil extending over x(0),x(1),...x(n).
c
c
      real x(0:n),c(0:n,0:m)
c
      c1       = 1.
      c4       = x(0) - xx
c
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
c
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo
c     call outmat(c,n+1,m+1,'fdw',n)
      return
      end
c-----------------------------------------------------------------------
      subroutine transpose_r(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
c
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine set_zgml (nxi)
      include 'basics.inc'
      integer nxo
      save    nxo
      data    nxo / 0 /

      nxl = min(nxi,nxm)

      if (nxo.ne.nxl) call legend(zgml,wght,nxl)

c     write(6,*) nxo,nxl,nxi,nxm,' POINTS'
c     call outmat(zgml,1,nxm,'setzg',nxl)

      nxo=nxl
      return
      end
c-----------------------------------------------------------------------
      subroutine msplite(ie,nxspi,nyspi,nzspi,raxi,rayi,razi)
C
C     This routine is a hack-copy of octsplite to generate 
c     multi-element decompositions of the square     (pff 9/28/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
      include 'basics.inc'

      real raxi(0:1),rayi(0:1),razi(0:1)

      parameter (nxm3=nxm*nym*nzm,ldw=2*nxm3,maxsp=100)
      common /ctmp2/ xp(nxm3),yp(nxm3),zp(nxm3),wk(ldw)
      common /ctmp0/ jx(nxm*3),jyt(nym*3),jzt(nzm*3)
     $             , rrx(0:maxsp),rry(0:maxsp),rrz(0:maxsp)
     $             , rxt(3),ryt(3),rzt(3)
     $             , xq(27),yq(27),zq(27)
      real jx,jyt,jzt

      integer eindx(12),e
      save    eindx
      data    eindx /  2 ,  6 ,  8 ,  4  ! Index of 12 edges into 3x3x3 tensor.
     $              , 20 , 24 , 26 , 22  ! Follows preproc. vtx notation
     $              , 10 , 12 , 18 , 16  /


c     ifmid  = .false.
c     ifcstd = .true.
c     write(6,*) nx,ny,nz,nxm,' NXM!'

      call genxyz_e (xp,yp,zp,ie,nx,ny,nz) ! high definition mesh

      nxsp = abs(nxspi)
      nysp = abs(nyspi)
      nzsp = abs(nzspi)

      call get_ratios(rrx,nxspi,raxi)
      call get_ratios(rry,nyspi,rayi)
      call get_ratios(rrz,nzspi,razi)

c     call outmat(rrx,1,nxsp+1,'rrx  ',nxspi)
c     call outmat(rry,1,nysp+1,'rry  ',nyspi)
c     call outmat(rrz,1,nzsp+1,'rrz  ',nzspi)
      call rzero(zq,27) ! for 2D

      call set_zgml (nx)

      e  = 0
      do k=1,nzsp

        if (if3d) then
          rzt(1)=rrz(k-1)
          rzt(3)=rrz(k  )
          rzt(2)=(rzt(1)+rzt(3))/2.
          call gen_int_gz(wk,jzt,rzt,3,zgml,nz)
        endif

        do j=1,nysp

          ryt(1)=rry(j-1)
          ryt(3)=rry(j  )
          ryt(2)=(ryt(1)+ryt(3))/2.
          call gen_int_gz(wk,jyt,ryt,3,zgml,ny)

          do i=1,nxsp

            rxt(1)=rrx(i-1)
            rxt(3)=rrx(i  )
            rxt(2)=(rxt(1)+rxt(3))/2.
            call gen_int_gz(jx,wk,rxt,3,zgml,nx)

            e  = e+1
            je = nel+e
            call copyel(ie,je)  ! copy ie to je
   
            if (if3d) then
              call tensr3(xq,3,3,3,xp,nx,ny,nz,jx,jyt,jzt,wk,ldw,if3d)
              call tensr3(yq,3,3,3,yp,nx,ny,nz,jx,jyt,jzt,wk,ldw,if3d)
              call tensr3(zq,3,3,3,zp,nx,ny,nz,jx,jyt,jzt,wk,ldw,if3d)
            else    ! 2D
              call tensr3(xq,3,3,3,xp,nx,ny,nz,jx,jyt,jzt,wk,ldw,if3d)
              call tensr3(yq,3,3,3,yp,nx,ny,nz,jx,jyt,jzt,wk,ldw,if3d)
            endif

            call q_to_neklin(x(1,je),1,xq,if3d)
            call q_to_neklin(y(1,je),1,yq,if3d)
            call q_to_neklin(z(1,je),1,zq,if3d)

c           write(6,*)
c           write(6,9) ie,(x(jj,ie),jj=1,4),' xold'
c           write(6,9) ie,(y(jj,ie),jj=1,4),' yold'
c           write(6,*)
c           write(6,9) je,(x(jj,je),jj=1,4),' xnew'
c           write(6,9) je,(y(jj,je),jj=1,4),' ynew'
c   9       format(i5,4f14.5,a5)

            call fix_curve(je,i,j,k,nxsp,nysp,nzsp,rrx,rry,rrz,xq,yq,zq)

            if (je.lt.500) call drawel(je)
 
          enddo
        enddo
      enddo
C
C     Copy last element generated back onto IE
C
      call copyel(je,ie)
      nel = nel + (e-1)
 
      return
      end
c-----------------------------------------------------------------------
      subroutine linquad(xl,yl,zl,nxl,nyl,nzl,e)

      include 'basics.inc'

      real xl(nxl*nyl*nzl),yl(nxl*nyl*nzl),zl(nxl*nyl*nzl)
      integer e

      logical ifmid_e

      nedge = 4 + 8*(ndim-2)

      ifmid_e = .false.
      do k=1,nedge
         if (ccurve(k,e).eq.'m') ifmid_e = .true.
      enddo

      if (ifmid_e) then
         call xyzquad(xl,yl,zl,nxl,nyl,nzl,e)
      else
         call xyzlin (xl,yl,zl,nxl,nyl,nzl,e)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine xyzlin(xl,yl,zl,nxl,nyl,nzl,e)
c     Generate bi- or trilinear mesh

      include 'basics.inc'

      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
      integer e

c   Preprocessor Corner notation:      Symmetric Corner notation:
c
c           4+-----+3    ^ s                    3+-----+4    ^ s
c           /     /|     |                      /     /|     |
c          /     / |     |                     /     / |     |
c        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r
c         |     | /     /                     |     | /     /
c         |     |/     /                      |     |/     /
c        5+-----+6    t                      5+-----+6    t

      integer indx(8)
      save    indx
      data    indx / 1,2,4,3,5,6,8,7 /

      parameter (ldw=2*lx1*ly1*lz1)
      common /ctmp0/ xcb(2,2,2),ycb(2,2,2),zcb(2,2,2),w(ldw)

c     Note : ctmpg is used in this format in several subsequent routines
      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)


      ndim2 = 2**ndim

      call set_zgml (nxl)

      do ix=1,nxl
         h(ix,1,1)=(1.0-zgml(ix))*0.5
         h(ix,1,2)=(1.0+zgml(ix))*0.5
      enddo
      do iy=1,nyl
         h(iy,2,1)=(1.0-zgml(iy))*0.5
         h(iy,2,2)=(1.0+zgml(iy))*0.5
      enddo
      if (if3d) then
         do iz=1,nzl
            h(iz,3,1)=(1.0-zgml(iz))*0.5
            h(iz,3,2)=(1.0+zgml(iz))*0.5
         enddo
      else
         call rone(h(1,3,1),nzl)
         call rone(h(1,3,2),nzl)
      endif

      do ix=1,ndim2
         i=indx(ix)
         xcb(ix,1,1)=x(i,e)  !  nek:   xc(i,e)
         ycb(ix,1,1)=y(i,e)  !  nek:   yc(i,e)
         zcb(ix,1,1)=z(i,e)  !  nek:   zc(i,e)
      enddo

c     Map R-S-T space into physical X-Y-Z space.

      ! NOTE:  Assumes nxl=nyl=nzl !

      call map_m_to_n(xl,nxl,xcb,2,if3d,w,ldw)
      call map_m_to_n(yl,nxl,ycb,2,if3d,w,ldw)
      call map_m_to_n(zl,nxl,zcb,2,if3d,w,ldw)

      return
      end
c-----------------------------------------------------------------------
      subroutine xyzquad(xl,yl,zl,nxl,nyl,nzl,e)
c     Generate bi- or triquadratic mesh

      include 'basics.inc'

      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)
      real xq(27),yq(27),zq(27)
      integer e

      parameter (ldw=lx1*ly1*lz1)
      common /ctmp0/ w(ldw,2),zg(3)

      integer eindx(12)  ! index of 12 edges into 3x3x3 tensor
      save    eindx      ! Follows preprocessor notation..
      data    eindx /  2 ,  6 ,  8 ,  4
     $              , 20 , 24 , 26 , 22
     $              , 10 , 12 , 18 , 16  /  ! preproc. vtx notation


      call xyzlin(xq,yq,zq,3,3,3,e) ! map bilin to bi- or triquadratic

      nedge = 4 + 8*(ndim-2)

      do k=1,nedge
         if (ccurve(k,e).eq.'m') then
            j = eindx(k)
            xq(j) = curve(1,k,e)
            yq(j) = curve(2,k,e)
            zq(j) = curve(3,k,e)
         endif
      enddo

      zg(1) = -1
      zg(2) =  0 
      zg(3) =  1

      if (if3d) then
         call gh_face_extend(xq,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
         call gh_face_extend(yq,zg,3,2,w(1,1),w(1,2))
         call gh_face_extend(zq,zg,3,2,w(1,1),w(1,2))

c        Map R-S-T space into physical X-Y-Z space.
         ! NOTE:  Assumes nxl=nyl=nzl !
         call map_m_to_n(xl,nxl,xq,3,if3d,w,ldw)
         call map_m_to_n(yl,nxl,yq,3,if3d,w,ldw)
         call map_m_to_n(zl,nxl,zq,3,if3d,w,ldw)

      else

         call gh_face_extend_2d(xq,zg,3,2,w(1,1),w(1,2)) ! 2 --> edge extend
         call gh_face_extend_2d(yq,zg,3,2,w(1,1),w(1,2))

         call map_m_to_n(xl,nxl,xq,3,if3d,w,ldw)
         call map_m_to_n(yl,nxl,yq,3,if3d,w,ldw)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine genxyz_e (xl,yl,zl,e,nxl,nyl,nzl)

      include 'basics.inc'

      real xl(nxl,nyl,nzl),yl(nxl,nyl,nzl),zl(nxl,nyl,nzl)

C     ctmpg is used in this format in several subsequent routines
      common /ctmpg/ h(lx1,3,2),xcrved(lx1),ycrved(ly1),zcrved(lz1)
     $             , work(3,lx1,lz1)

      character*1 ccv

      integer e

      call set_zgml (nxl)
c     Initialize geometry arrays with bi- triquadratic deformations
      call linquad(xl,yl,zl,nxl,nyl,nzl,e)
      call drawqel(xl,yl,nxl,nyl,e,4)

c     Deform surfaces - general 3D deformations
c                     - extruded geometry deformations
      call set_zgml (nxl)

      nfaces = 2*ndim
      do iface=1,nfaces
        ccv = ccurve(iface,e)
        if (ccv.eq.'s') 
     $     call sphsrf_e(xl,yl,zl,iface,e,nxl,nyl,nzl,work) 
      enddo

      do isid=1,8
        ccv = ccurve(isid,e)
        if (ccv.eq.'C') call arcsrf_e(xl,yl,zl,nxl,nyl,nzl,e,isid)
c       write(6,*) ccv,isid,e,nxl,' ccv'
c       if (ccv.eq.'C') call drawqel(xl,yl,nxl,nyl,e,14)
      enddo

      if (mod(nxl,2).eq.1) then
         nxh = nxl/2
         l=0
         do k=1,nzl,nxh
         do j=1,nyl,nxh
         do i=1,nxl,nxh
            l = l+1
            x27(l,e) = xl(i,j,k)
            y27(l,e) = yl(i,j,k)
            z27(l,e) = zl(i,j,k)
c           write(6,1) x27(l,e),y27(l,e),z27(l,e),l,i,j,k,e,'x27'
  1         format(1p3e12.4,4i4,i6,a4)
         enddo
         enddo
         enddo
      else
         if (e.eq.1) write(6,*) 'Warning: nxl even in genxyz',nxl
         if (e.eq.1) call exitt
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_hlin(h,mx) ! linear to mx mapping
      include 'basics.inc'
      real h(mx,3,2)

      call set_zgml (mx)

      do 11 ix=1,mx
         h(ix,1,1)=(1.0-zgml(ix))*0.5
         h(ix,1,2)=(1.0+zgml(ix))*0.5
   11 continue
      do 21 iy=1,mx
         h(iy,2,1)=(1.0-zgml(iy))*0.5
         h(iy,2,2)=(1.0+zgml(iy))*0.5
   21 continue
      if (if3d) then
         do 31 iz=1,mx
            h(iz,3,1)=(1.0-zgml(iz))*0.5
            h(iz,3,2)=(1.0+zgml(iz))*0.5
   31    continue
      else
         call rone(h(1,3,1),mx)
         call rone(h(1,3,2),mx)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine q_to_neklin(x,ldx,xq,if3d)
      real x(ldx,8),xq(27)
      logical if3d
      integer inv(8),v
      save    inv
      data    inv /1,2,4,3,5,6,8,7/

      if (if3d) then

        v=0
        do kk=1,3,2
        do jj=1,3,2
        do ii=1,3,2
           v   = v+1
           ijk = ii + 3*(jj-1) + 9*(kk-1)
           x(1,inv(v))=xq(ijk)
        enddo
        enddo
        enddo

      else    ! 2D

        v=0
        do jj=1,3,2
        do ii=1,3,2
           v   = v+1
           ijk = ii + 3*(jj-1)
           x(1,inv(v))=xq(ijk)
        enddo
        enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine plane_project_mesh(x0,p1,p2,p3,e0,e1,tol)
      include 'basics.inc'
      real x0(3),p1(3),p2(3),p3(3)
      real nh(3),v0(3),v2(3),v3(3)
      integer e,e0,e1


      call sub3(v2,p2,p1,3)
      call sub3(v3,p3,p1,3)
      call vcross_normal(nh,v2,v3)

      call sub3(v0,x0,p1,3)

      alpha = dot(v0,nh,3)
      if (alpha.lt.0) call chsign(nh,3)

      do e=e0,e1
        do i=1,8
           d2plane = (x(i,e)-p1(1))*nh(1)
     $             + (y(i,e)-p1(2))*nh(2)
     $             + (z(i,e)-p1(3))*nh(3)
           if (d2plane.lt.tol) then
              x(i,e) = x(i,e) - d2plane*nh(1)
              y(i,e) = y(i,e) - d2plane*nh(2)
              z(i,e) = z(i,e) - d2plane*nh(3)
           endif
        enddo

        do i=1,12 ! Check midside nodes
           if (ccurve(i,e).eq.'m') then
              d2plane = (curve(1,i,e)-p1(1))*nh(1)
     $                + (curve(2,i,e)-p1(2))*nh(2)
     $                + (curve(3,i,e)-p1(3))*nh(3)
              if (d2plane.lt.tol) then ! PROJECT
                 curve(1,i,e) = curve(1,i,e)-d2plane*nh(1)
                 curve(2,i,e) = curve(2,i,e)-d2plane*nh(2)
                 curve(3,i,e) = curve(3,i,e)-d2plane*nh(3)
              endif
           endif
        enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine draw3v(x3,y3,z3)
      x=xiso(x3,y3,z3)
      y=yiso(x3,y3,z3)
      write(6,1) x,y,x3,y3,z3
    1 format(1p5e11.3,' draw3')
      call drawc(x,y)
C     call drawsc(xscr(x),yscr(y))
      return
      end
c-----------------------------------------------------------------------
      subroutine move3v(x3,y3,z3)
      x=xiso(x3,y3,z3)
      y=yiso(x3,y3,z3)
      write(6,1) x,y,x3,y3,z3
    1 format(1p5e11.3,' move3')
      call movec(x,y)
      return
      end
c-----------------------------------------------------------------------
      subroutine drawqel(xl,yl,nxl,nyl,e,ic)
      real xl(nxl,nyl),yl(nxl,nyl)
      integer e,ic
      character*1 ans

      return

      call color(ic)

      write(6,*) nxl,nyl,e,ic,' nxl,nyl,e,ic' 
      do i=1,nxl
         call move3(xl(i,1),yl(i,1),0.)
         do j=2,nyl
            call draw3(xl(i,j),yl(i,j),0.)
         enddo
      enddo

      do j=1,nyl
         call move3(xl(1,j),yl(1,j),0.)
         do i=2,nxl
            call draw3(xl(i,j),yl(i,j),0.)
         enddo
      enddo

      call prsi('Continue?$',e)
      call res(ans,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine octsplite(ie,liste)
C
C     A routine which will split an element IE and assign the 
C     appropriate fraction to element JE.   IPLN indicates the
C     type of split.   FRAC indicates the sub-division point on
C     the range [0,1].
C
      include 'basics.inc'
      PARAMETER (NXM3=NXM*NYM*NZM)
      COMMON /CTMP2/ XP(NXM3),YP(NXM3),ZP(NXM3),RRL(3)
     $              ,XVAL(-1:1,-1:1,-1:1),YVAL(-1:1,-1:1,-1:1)
     $              ,ZVAL(-1:1,-1:1,-1:1)
      INTEGER INV(8)
      SAVE    INV
      DATA    INV /1,2,4,3,5,6,8,7/
      DIMENSION LISTE(8)

      call genxyz_e (xp,yp,zp,ie,nx,ny,nz)  ! Compute XYZ on high definition mesh

      if (if3d) then       !   Evaluate x,y,z on 27 pt stencil
        do 160 It= -1,1,1
        do 160 Is= -1,1,1
        do 160 Ir= -1,1,1
          rrl(3)=float(it)
          rrl(2)=float(is)
          rrl(1)=float(ir)
          call evalsc(xval(ir,is,it),xp,rrl,1)
          call evalsc(yval(ir,is,it),yp,rrl,0)
          call evalsc(zval(ir,is,it),zp,rrl,0)
  160   continue

        ke = 0
        do 860 Iz=-1,0,1 ! Eight quadrants
        do 840 Iy=-1,0,1
        do 820 Ix=-1,0,1

          ke = ke+1
          je = liste(ke)
          call copyel(ie,je) ! copy essential information

          iv = 0
          do 760 Jz=0,1 !         Eight vertices
          do 760 Jy=0,1
          do 760 Jx=0,1
           jz0 = Jz+Iz
           jy0 = Jy+Iy
           jx0 = Jx+Ix
           iv = iv+1
           x(inv(iv),je)=xval(jx0,jy0,jz0)
           y(inv(iv),je)=yval(jx0,jy0,jz0)
           z(inv(iv),je)=zval(jx0,jy0,jz0)
  760     continue

          do 400 Ifce=1,3
           IFAC1=EFACE(2*Ifce-1)
           IFAC2=EFACE(2*Ifce  )
           if (Ifce.eq.1) then
              JFAC = EFACE(1-Ix)
           elseif (Ifce.eq.2) then
              JFAC = EFACE(3-Iy)
           ELSE
              JFAC = EFACE(5-Iz)
           ENDIF
C
C          Curve side fix up - only internal faces have to be changed
C
           if (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') then
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           elseif (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') then
             R1i = 0.0
             R2i = 0.0
             if (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             if (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF

C          Check for cylinder (Note this scheme fails for convex-convex case!) pff 4/10/93

           if (ifce.le.2) then
              do 300 Ilev=1,2
                if (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') then
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  CCURVE(JFAC,JE)  = 'C'
                  if (JFAC.eq.IFAC1) CURVE(1,JFAC,JE) =  R
                  if (JFAC.eq.IFAC2) CURVE(1,JFAC,JE) = -R
                elseif (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') then
                  R1i = 0.0
                  R2i = 0.0
                  if (CCURVE(IFAC1,IE).EQ.'C') 
     $               R1i = 1.0/CURVE(1,IFAC1,IE)
                  if (CCURVE(IFAC2,IE).EQ.'C') 
     $               R2i = 1.0/CURVE(1,IFAC2,IE)
                  R  = 2.0/(R1i-R2i)
c                 CURVE(1,JFAC,JE) =  R
                  CCURVE(JFAC,JE)=' '
                ENDIF
                JFAC =JFAC +4
                IFAC1=IFAC1+4
                IFAC2=IFAC2+4
  300         continue
           ENDIF
  400     continue
  820  continue
  840  continue
  860  continue

      else ! 2D

        do 2160 Is= -1,1,1
        do 2160 Ir= -1,1,1
           rrl(3)=0.0
           rrl(2)=float(is)
           rrl(1)=float(ir)
           call evalsc(xval(ir,is,1),xp,rrl,1)
           call evalsc(yval(ir,is,1),yp,rrl,0)
 2160  continue

      ke = 0
      do 2840 Iy=-1,0,1      ! Four quadrants
      do 2820 Ix=-1,0,1

          ke = ke+1
          je = liste(ke)

          call copyel(ie,je) ! copy essential information

          iv = 0
          do 2760 Jy=0,1     ! Four vertices
          do 2760 Jx=0,1
             jy0 = Jy+Iy
             jx0 = Jx+Ix
             iv = iv+1
             x(inv(iv),je)=xval(jx0,jy0,1)
             y(inv(iv),je)=yval(jx0,jy0,1)
             z(inv(iv),je)=zval(jx0,jy0,1)
 2760     continue
C
          do 2400 ifce=1,2
           ifac1=eface(2*Ifce-1)
           ifac2=eface(2*Ifce  )
           if (ifce.eq.1) then
              jfac = eface(1-Ix)
           elseif (Ifce.eq.2) then
              jfac = eface(3-iy)
           else
              jfac = eface(5-iz)
           endif
C
C        Curve side fix up - only internal faces have to be changed
C
           if (CCURVE(IFAC1,IE).EQ.'s'  .and.
     $         CCURVE(IFAC2,IE).EQ.'s') then
             R = 0.5 * ( CURVE(4,IFAC1,IE) + CURVE(4,IFAC2,IE) )
             CCURVE(JFAC,JE)='s'
             CURVE(4,JFAC,JE) = R
           elseif (CCURVE(IFAC1,IE).EQ.'s'  .or.
     $             CCURVE(IFAC2,IE).EQ.'s') then
             R1i = 0.0
             R2i = 0.0
             if (CCURVE(IFAC1,IE).EQ.'s') R1i = 1.0/CURVE(4,IFAC1,IE)
             if (CCURVE(IFAC2,IE).EQ.'s') R2i = 1.0/CURVE(4,IFAC2,IE)
             R  = 2.0/(R1i+R2i)
c            CCURVE(JFAC,JE)='s'
c            CURVE(4,JFAC,JE) = R
             CCURVE(JFAC,JE)=' '
           ENDIF
C          Check for cylinder (Note this scheme fails for convex-convex case!) pff 4/10/93

           if (Ifce.LE.2) then
              ilev = 1
              if (CCURVE(IFAC1,IE).EQ.'C'.and.
     $              CCURVE(IFAC2,IE).EQ.'C') then
                  R1 = CURVE(1,IFAC1,IE)
                  R2 = CURVE(1,IFAC2,IE)
                  R = 0.5*(R1-R2)
                  ccurve(jfac,je)  = 'C'
                  if (jfac.EQ.ifac1) curve(1,jfac,je) =  r
                  if (jfac.EQ.ifac2) curve(1,jfac,je) = -r
              elseif (CCURVE(IFAC1,IE).EQ.'C'.  or.
     $                  CCURVE(IFAC2,IE).EQ.'C') then
                  R1i = 0.0
                  R2i = 0.0
                  if (ccurve(ifac1,ie).eq.'C') r1i=1/curve(1,ifac1,ie)
                  if (ccurve(ifac2,ie).eq.'C') r2i=1/curve(1,ifac2,ie)
                  R  = 2.0/(R1i-R2i)
c                 curve(1,jfac,je) =  r
                  ccurve(jfac,je)=' '
              ENDIF
           ENDIF
 2400     continue
 2820   continue
 2840   continue
 2860   continue

      endif  ! End of 2D/3D decision branch

      call copyel(je,ie) ! copy last element generated back onto ie

      write(6,*) 'quit in octsplite'
      stop

      return
      end
c-----------------------------------------------------------------------
      subroutine octsplitn(ie,liste)
C
C     This routine is a hack-copy of octsplite to generate 
c     multi-element decompositions of the square     (pff 9/28/05)
c
C     It simply modifies in the X-Y plane, but doesn't refine in Z.
C
C
      include 'basics.inc'

      nxsp=2
      nysp=2
      nzsp=1
      if (if3d) nzsp=2

      rax = 1.
      ray = 1.
      raz = 1.

      call msplite(ie,nxsp,nysp,nzsp,rax,ray,raz)

      return
      end
c-----------------------------------------------------------------------
