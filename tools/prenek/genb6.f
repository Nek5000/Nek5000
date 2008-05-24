C
C     5 Aug 1994 updated  Paul F. Fischer
c     18 Jan 1996  pff
c     3/19 1996    pff
C
c     This program generates a nekton .rea file which comprises
c     a number of *tensor-product* boxes, i.e., a set of Nelx x Nely x Nelz
c     elements, where the locations of the vertices of the elements are
c     given as one-dimensional arrays of length (Nelx+1), (Nely+1), and (Nelz+1), 
c     respectively.
c
C    .For each box, the User specifies the number of elements in the 
c     x-y-z directions and the corresponding locations of the vertices
c     and the boundary conditions on the box
c
C    .INPUT FILE FORMAT:
c
c     line  0:    name of .rea file (start in col. 1, include .rea suffix)
c     line  1:    ndim              (MUST be same as in above .rea file! )
c     line  2:    Nelx, Nely, Nelz Nfld   for Box 1 
c     line  3:    x_0  x_1 ....  x_Nelx   for Box 1
c     line  4:    y_0  y_1 ....  y_Nely   for Box 1
c     line  5:    z_0  z_1 ....  z_Nelz   for Box 1
c     line  6:    bc1,bc2,bc3,bc4,bc5,bc6 for Box 1, fld 1
c     line  7:    bc1,bc2,bc3,bc4,bc5,bc6 for Box 1, fld 2  (If Nfld > 1)
c     line  8:    Nelx, Nely, Nelz Nfld   for Box 2 
c     line  9:          etc.
c
c
c    .Note that the boundary conditions are read with (6(a3,1x)) format, so
c     those lines must start in column 1 of the input file.  The rest is
c     free formatted, so any standard input format is ok.
c
c    .A typical input file for a 3D case is the following:
c
c>#
c>#========================================================
c>#    comments
c>#
c>#    This is the set of box data for full spherical configuration
c>#
c>#========================================================
c>#
c>ps3.rea
c>3                         NDIM
c>3                         NFLDS
c>box_1                     (Any string with 1st character .ne. "c" or "C")
c>-4   4  -4                (nelx,nely,nelz for Box) 
c>0 1.                      (if nelx<0, length to be equally divided)
c>0 .25 .5 .75 1.0          (if nely>0, nely+1 y-values for start/stop of els)
c>0 1.
c>P  ,P  ,P  ,P  ,P  ,P     (fluid bc's)
c>P  ,P  ,P  ,P  ,P  ,P     (temp  bc's)
c>P  ,P  ,P  ,P  ,P  ,P     (passive scalar bc's)
c>box_2                     (Any string with 1st character .ne. "c" or "C")
c>-4   4  -4                (nelx,nely,nelz for Box) 
c>1 2.                      (if nelx<0, length to be equally divided)
c>0 .25 .5 .75 1.0          (if nely>0, nely+1 y-values for start/stop of els)
c>0 1.
c>P  ,P  ,P  ,P  ,P  ,P     (fluid bc's)
c>P  ,P  ,P  ,P  ,P  ,P     (temp  bc's)
c>P  ,P  ,P  ,P  ,P  ,P     (passive scalar bc's)
c
c
c------------------------------------------------------------------------------
      program genbox
      character*80 string
      character*1  string1(80)
      equivalence (string,string1)
C
      logical     if3d,ifevenx,ifeveny,ifevenz,if_multi_seg
      character*3 cbc1,   cbc2,   cbc3,   cbc4,   cbc5,   cbc6
      real        rbc1(5),rbc2(5),rbc3(5),rbc4(5),rbc5(5),rbc6(5)
      integer eface(6)
      save    eface
      data    eface  / 4,2,1,3,5,6 /
      integer nelxyz(3)
c
      character*1  apt(52)
      character*52 apt52
      equivalence (apt52,apt)
      data apt52/'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      parameter (mbox  = 60)
      parameter (maxx  = 160)
      parameter (maxel = mbox*maxx*maxx)
      integer nlx(mbox),nly(mbox),nlz(mbox)
      real x(0:maxx,mbox),y(0:maxx,mbox),z(0:maxx,mbox)
      real xc(mbox),yc(mbox),zc(mbox)
      real curve(8,maxel)
      character*3 cbc(6,mbox,2)
      character*1 boxcirc(mbox)
c
c
c     Get the input file, which specifies the base .rea file
      call blank(string,80)
      write(6,*) 'input file name:'
      read (5,80) string
   80 format(a80)
      open (unit=7,file=string,status='old')
c
C     Read in name of previously generated NEKTON data set.
C     (Must have same dimension and number of fields as current run)
C
      call gets(string,80,iend,7)
c
      open (unit=8,file=string,status='old')
      open (unit=9,file='box.rea')
      call scanout(string,'NEL,NDIM,NELV',13,8,9)
c
c-----------------------------------------------------------------------
c     here is where the 2d/3d determination is made....
c
      if3d = .false.
      call geti1(ndim,iend,7)
      if (ndim.eq.3) if3d=.true.
c
c-----------------------------------------------------------------------
c     here is where the fluid, fluid+heat determination is made....
c
      call geti1(nfld,iend,7)
c-----------------------------------------------------------------------
C
      nel  = 0
      nbox = 0
      do ibox=1,mbox
c
         ifevenx = .false.
         ifeveny = .false.
         ifevenz = .false.
c
         call gets(boxcirc(ibox),1,iend,7)
         if (iend.eq.1) goto 99
         if (boxcirc(ibox).eq.'c'.or.boxcirc(ibox).eq.'C') 
     $                       call getr2(xc(ibox),yc(ibox),iend,7)
         if_multi_seg = .false.
         if (boxcirc(ibox).eq.'m'.or.boxcirc(ibox).eq.'M') 
     $      if_multi_seg = .true.
       
         if (if_multi_seg) then
            call get_multi_seg
     $           (nelxyz,x(0,ibox),y(0,ibox),z(0,ibox),maxx,if3d)
            nelx = nelxyz(1)
            nely = nelxyz(2)
            nelz = nelxyz(3)
            do ifld=1,nfld
               call getcv(cbc(1,ibox,ifld),3,6,iend,7)
               write(6,*) 'CBC:',(cbc(k,ibox,ifld),k=1,6),ifld,ibox
            enddo
         elseif (if3d) then
            call geti3(nelx,nely,nelz,iend,7)
c
            if (nelx.lt.0) ifevenx = .true.
            if (nely.lt.0) ifeveny = .true.
            if (nelz.lt.0) ifevenz = .true.
            nelx = abs(nelx)
            nely = abs(nely)
            nelz = abs(nelz)
c
            if (nelx.gt.maxx.or.nely.gt.maxx.or.nelz.gt.maxx) then
               write(6,*) 'Abort, increase maxx and recompile',
     $                       nelx,nely,nelz,maxx
               stop
            endif
c
            nery = nely
            if (boxcirc(ibox).eq.'C') nery = 1
            if (iend.eq.1) goto 99
            ninbox = nelx*nely*nelz
            write(6,6) ninbox,nelx,nely,nelz,ibox,nfld
c
            if (ifevenx) then
               call getr3(x0,x1,ratio,iend,7)
               if (ratio.le.0) ratio=1.
               dx = (x1-x0)/nelx
               x(0,ibox) = 0.
               do k=1,nelx
                  x(k,ibox) = x(k-1,ibox) + dx
                  dx        = dx*ratio
               enddo
               scale = (x1-x0)/x(nelx,ibox)
               do k=0,nelx
                  x(k,ibox) = x0 + scale*x(k,ibox)
               enddo
            else
               call getrv(x(0,ibox),nelx+1,iend,7)
            endif
c
            if (ifeveny) then
               call getr3(y0,y1,ratio,iend,7)
               if (ratio.le.0) ratio=1.
               dy = (y1-y0)/nely
               y(0,ibox) = 0.
               do k=1,nely
                  y(k,ibox) = y(k-1,ibox) + dy
                  dy        = dy*ratio
               enddo
               scale = (y1-y0)/y(nely,ibox)
               do k=0,nely
                  y(k,ibox) = y0 + scale*y(k,ibox)
               enddo
            else
               call getrv(y(0,ibox),nely+1,iend,7)
            endif
c
            if (ifevenz) then
               call getr3(z0,z1,ratio,iend,7)
               if (ratio.le.0) ratio=1.
               dz = (z1-z0)/nelz
               z(0,ibox) = 0.
               do k=1,nelz
                  z(k,ibox) = z(k-1,ibox) + dz
                  dz        = dz*ratio
               enddo
               scale = (z1-z0)/z(nelz,ibox)
               do k=0,nelz
                  z(k,ibox) = z0 + scale*z(k,ibox)
               enddo
            else
               call getrv(z(0,ibox),nelz+1,iend,7)
            endif
c
            do ifld=1,nfld
               call getcv(cbc(1,ibox,ifld),3,6,iend,7)
               write(6,*) 'CBC:',(cbc(k,ibox,ifld),k=1,6),ifld,ibox
            enddo
         else
            nelz = 1
            call geti2(nelx,nely,iend,7)
c
            if (nelx.lt.0) ifevenx = .true.
            if (nely.lt.0) ifeveny = .true.
            nelx = abs(nelx)
            nely = abs(nely)
            if (nelx.gt.maxx.or.nely.gt.maxx) then
               write(6,*) 'Abort, increase maxx and recompile',
     $                       nelx,nely,maxx
               stop
            endif
c
            nery = nely
            if (boxcirc(ibox).eq.'C') nery = 1
            if (iend.eq.1) goto 99
            ninbox = nelx*nely*nelz
            write(6,6) ninbox,nelx,nely,nery,ibox,nfld
c
            if (ifevenx) then
               call getr3(x0,x1,ratio,iend,7)
               if (ratio.le.0) ratio=1.
               dx = (x1-x0)/nelx
               x(0,ibox) = 0.
               do k=1,nelx
                  x(k,ibox) = x(k-1,ibox) + dx
                  dx        = dx*ratio
               enddo
               scale = (x1-x0)/x(nelx,ibox)
               do k=0,nelx
                  x(k,ibox) = x0 + scale*x(k,ibox)
               enddo
            else
               call getrv(x(0,ibox),nelx+1,iend,7)
            endif
c
            if (ifeveny) then
               call getr3(y0,y1,ratio,iend,7)
               if (ratio.le.0) ratio=1.
               dy = (y1-y0)/nely
               y(0,ibox) = 0.
               do k=1,nely
                  y(k,ibox) = y(k-1,ibox) + dy
                  dy        = dy*ratio
               enddo
               scale = (y1-y0)/y(nely,ibox)
               do k=0,nely
                  y(k,ibox) = y0 + scale*y(k,ibox)
               enddo
            else
               call getrv(y(0,ibox),nely+1,iend,7)
            endif
c
            do ifld=1,nfld
               call getcv(cbc(1,ibox,ifld),3,4,iend,7)
            enddo
         endif
         nlx(ibox) = nelx
         nly(ibox) = nely
         nlz(ibox) = nelz
         nel = nel + nelx*nely*nelz
         nbox = nbox+1
      enddo
    6 format('Reading',i6,' =',3i6,' elements for box',i3,'.')
   99 continue
c
      write(6,*) 'Beginning construction of new .rea file to box.rea'
      write(6,*) nel,' elements will be created for',nbox,' boxes.'
      write(6,*) 
c
c     Construct mesh
c
      write(9,10) nel,ndim,nel
   10 format(3i6,'           NEL,NDIM,NELV')
c
      ncurv = 0
      call rzero(curve,8*nel)
c
      if (if3d) then
        ie   = 0
        ilev = 0
        do ibox=1,nbox
          do k=1,nlz(ibox)
            z1 = z(k-1,ibox)
            z2 = z(k  ,ibox)
c
            ia   = 0
c           ilev = ilev+1
            ilev = k
            do j=1,nly(ibox)
               y1 = y(j-1,ibox)
               y2 = y(j  ,ibox)
c
               do i=1,nlx(ibox)
                  x1 = x(i-1,ibox)
                  x2 = x(i  ,ibox)
c
                  ie = ie+1
                  ia = ia+1
                  ia = mod1(ia,52)
                  write(9,11) ie,ilev,apt(ia),'0'
c
                  write(9,12) x1,x2,x2,x1
                  write(9,12) y1,y1,y2,y2
                  write(9,12) z1,z1,z1,z1
c
                  write(9,12) x1,x2,x2,x1
                  write(9,12) y1,y1,y2,y2
                  write(9,12) z2,z2,z2,z2
               enddo
            enddo
          enddo
        enddo
   12   format(4g14.6)
   11   format(
     $'            ELEMENT',i5,' [',i5,a1,']    GROUP     ',a1)
c    $'            ELEMENT',i5,' [',i5,a1,']    GROUP     0')
      else
        ilev = 1
        ie=0
        ia=0
        one    = 1.
        pi     = 4.*atan(one)
        pi180  = pi/180.
        do ibox=1,nbox
           if (boxcirc(ibox).eq.'c') then
c             CIRCLE (x = r, y = theta) .... specified in std. "box" way...
              do j=1,nly(ibox)
                 t1 = pi180*y(j-1,ibox)
                 t2 = pi180*y(j  ,ibox)
                 do i=1,nlx(ibox)
                    r1 = x(i-1,ibox)
                    r2 = x(i  ,ibox)
                    ie = ie+1
                    ia = ia+1
                    ia = mod1(ia,52)
                    x1 = xc(ibox) + r1*cos(t1)
                    x2 = xc(ibox) + r2*cos(t1)
                    x3 = xc(ibox) + r2*cos(t2)
                    x4 = xc(ibox) + r1*cos(t2)
                    y1 = yc(ibox) + r1*sin(t1)
                    y2 = yc(ibox) + r2*sin(t1)
                    y3 = yc(ibox) + r2*sin(t2)
                    y4 = yc(ibox) + r1*sin(t2)
                    curve(2,ie) =  r2
                    curve(4,ie) = -r1
c                   write(9,11) ie,ilev,apt(ia),'0'
                    write(9,11) ie,ilev,' ','0'
                    write(9,12) x1,x2,x3,x4
                    write(9,12) y1,y2,y3,y4
                    ncurv = ncurv+2
                 enddo
              enddo
           elseif (boxcirc(ibox).eq.'C') then
c             CIRCLE (x = r, y = theta) .... specified in Cool fast way...
              dt = y(1,ibox)
              t2 = y(0,ibox)
              do j=1,nly(ibox)
                 t1 = t2
                 t2 = t1 + dt
                 p1 = pi180*t1
                 p2 = pi180*t2
                 write(6,*) 'p1:',p1,t1,t2
                 do i=1,nlx(ibox)
                    r1 = x(i-1,ibox)
                    r2 = x(i  ,ibox)
                    ie = ie+1
                    ia = ia+1
                    ia = mod1(ia,52)
                    x1 = xc(ibox) + r1*cos(p1)
                    x2 = xc(ibox) + r2*cos(p1)
                    x3 = xc(ibox) + r2*cos(p2)
                    x4 = xc(ibox) + r1*cos(p2)
                    y1 = yc(ibox) + r1*sin(p1)
                    y2 = yc(ibox) + r2*sin(p1)
                    y3 = yc(ibox) + r2*sin(p2)
                    y4 = yc(ibox) + r1*sin(p2)
                    curve(2,ie) =  r2
                    curve(4,ie) = -r1
                    write(9,11) ie,ilev,apt(ia),'0'
                    write(9,12) x1,x2,x3,x4
                    write(9,12) y1,y2,y3,y4
                    ncurv = ncurv+2
                 enddo
              enddo
           else
c             BOX
              do j=1,nly(ibox)
                 y1 = y(j-1,ibox)
                 y2 = y(j  ,ibox)
                 do i=1,nlx(ibox)
                    x1 = x(i-1,ibox)
                    x2 = x(i  ,ibox)
                    ie = ie+1
                    ia = ia+1
                    ia = mod1(ia,52)
                    write(9,11) ie,ilev,apt(ia),'0'
                    write(9,12) x1,x2,x2,x1
                    write(9,12) y1,y1,y2,y2
                 enddo
              enddo
           endif
        enddo
      endif
c
c     output curve stuff and Boundary conditions
c                                     .... **for now, curved sides unsupported**
      write(9,28) ncurv
   28 format(
     $ '  ***** CURVED SIDE DATA *****',/,
     $    i6,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')
c
      zero = 0.
      if (nel.lt.1000) then
         do ie=1,nel
            do iedge = 2,4,2
               if (curve(iedge,ie).ne.0) write(9,290) 
     $            iedge,ie,curve(iedge,ie),(zero,k=1,4),'C'
            enddo
         enddo
  290    format(i3,i3,5g14.6,1x,a1)
       else
         do ie=1,nel
            do iedge = 2,4,2
               if (curve(iedge,ie).ne.0) write(9,291) 
     $            iedge,ie,curve(iedge,ie),(zero,k=1,4),'C'
            enddo
         enddo
  291    format(i2,i6,5g14.6,1x,a1)
       endif
c
      write(9,31) 
   31 format('  ***** BOUNDARY CONDITIONS *****')
c
      do ifld=1,nfld
c
         if2 = ifld-2
         if (ifld.eq.1) write(9,41)
         if (ifld.eq.2) write(9,42)
         if (ifld.eq.3) write(9,43) if2
   41    format('  ***** FLUID   BOUNDARY CONDITIONS *****')
   42    format('  ***** THERMAL BOUNDARY CONDITIONS *****')
   43    format('  ***** PASSIVE SCALAR',i2
     $                       ,'  BOUNDARY CONDITIONS *****')
c
         ie = 0
         if (if3d) then
            do ibx=1,nbox
            do iez=1,nlz(ibx)
            do iey=1,nly(ibx)
            do iex=1,nlx(ibx)
               nelx = nlx(ibx)
               nely = nly(ibx)
               nelz = nlz(ibx)
               ilftz = mod1(iez-1+nelz,nelz)
               irgtz = mod1(iez+1+nelz,nelz)
               ilfty = mod1(iey-1+nely,nely)
               irgty = mod1(iey+1+nely,nely)
               ilftx = mod1(iex-1+nelx,nelx)
               irgtx = mod1(iex+1+nelx,nelx)
c
               rbc1(1) = nely*nelx*(iez-1)+nelx*(iey-1) + ilftx
               rbc1(2) = eface(2)
               rbc2(1) = nely*nelx*(iez-1)+nelx*(iey-1) + irgtx
               rbc2(2) = eface(1)
c
               rbc3(1) = nely*nelx*(iez-1)+nelx*(ilfty-1) + iex
               rbc3(2) = eface(4)
               rbc4(1) = nely*nelx*(iez-1)+nelx*(irgty-1) + iex
               rbc4(2) = eface(3)
c
               rbc5(1) = nely*nelx*(ilftz-1)+nelx*(iey-1) + iex
               rbc5(2) = eface(6)
               rbc6(1) = nely*nelx*(irgtz-1)+nelx*(iey-1) + iex
               rbc6(2) = eface(5)
c
               if (iex.eq.1) then
                  cbc1=cbc(1,ibx,ifld)
                  if (cbc1.ne.'P  ') call rzero(rbc1,5)
               else
                  cbc1='E  '
               endif
c
               if (iex.eq.nelx) then
                  cbc2=cbc(2,ibx,ifld)
                  if (cbc2.ne.'P  ') call rzero(rbc2,5)
               else
                  cbc2='E  '
               endif
c
               if (iey.eq.1) then
                  cbc3=cbc(3,ibx,ifld)
                  if (cbc3.ne.'P  ') call rzero(rbc3,5)
               else
                  cbc3='E  '
               endif
c
               if (iey.eq.nely) then
                  cbc4=cbc(4,ibx,ifld)
                  if (cbc4.ne.'P  ') call rzero(rbc4,5)
               else
                  cbc4='E  '
               endif
c
               if (iez.eq.1) then
                  cbc5=cbc(5,ibx,ifld)
                  if (cbc5.ne.'P  ') call rzero(rbc5,5)
               else
                  cbc5='E  '
               endif
c
               if (iez.eq.nelz) then
                  cbc6=cbc(6,ibx,ifld)
                  if (cbc6.ne.'P  ') call rzero(rbc6,5)
               else
                  cbc6='E  '
               endif
c
c              output bc's in preproc. notation
c
               ie = ie+1
               if (nel.lt.1000) then
                  write(9,20) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                  write(9,20) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                  write(9,20) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                  write(9,20) cbc1,ie,eface(1),(rbc1(j),j=1,5)
                  write(9,20) cbc5,ie,eface(5),(rbc5(j),j=1,5)
                  write(9,20) cbc6,ie,eface(6),(rbc6(j),j=1,5)
               else
                  write(9,21) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                  write(9,21) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                  write(9,21) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                  write(9,21) cbc1,ie,eface(1),(rbc1(j),j=1,5)
                  write(9,21) cbc5,ie,eface(5),(rbc5(j),j=1,5)
                  write(9,21) cbc6,ie,eface(6),(rbc6(j),j=1,5)
               endif
   20          format(1x,a3,2i3,5g14.7)
   21          format(1x,a3,i5,i1,5g14.7)
            enddo
            enddo
            enddo
            enddo
         else
            do ibx=1,nbox
            do iey=1,nly(ibx)
            do iex=1,nlx(ibx)
               nelx = nlx(ibx)
               nely = nly(ibx)
               ilfty = mod1(iey-1+nely,nely)
               irgty = mod1(iey+1+nely,nely)
               ilftx = mod1(iex-1+nelx,nelx)
               irgtx = mod1(iex+1+nelx,nelx)
c
               rbc1(1) = nelx*(iey-1) + ilftx
               rbc1(2) = eface(2)
               rbc2(1) = nelx*(iey-1) + irgtx
               rbc2(2) = eface(1)
c
               rbc3(1) = nelx*(ilfty-1) + iex
               rbc3(2) = eface(4)
               rbc4(1) = nelx*(irgty-1) + iex
               rbc4(2) = eface(3)
c
               if (iex.eq.1) then
                  cbc1=cbc(1,ibx,ifld)
                  if (cbc1.ne.'P  ') call rzero(rbc1,5)
               else
                  cbc1='E  '
               endif
c
               if (iex.eq.nelx) then
                  cbc2=cbc(2,ibx,ifld)
                  if (cbc2.ne.'P  ') call rzero(rbc2,5)
               else
                  cbc2='E  '
               endif
c
               if (iey.eq.1) then
                  cbc3=cbc(3,ibx,ifld)
                  if (cbc3.ne.'P  ') call rzero(rbc3,5)
               else
                  cbc3='E  '
               endif
c
               if (iey.eq.nely) then
                  cbc4=cbc(4,ibx,ifld)
                  if (cbc4.ne.'P  ') call rzero(rbc4,5)
               else
                  cbc4='E  '
               endif
c
c              output bc's in preproc. notation
c
               ie = ie+1
               if (nel.lt.1000) then
                  write(9,20) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                  write(9,20) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                  write(9,20) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                  write(9,20) cbc1,ie,eface(1),(rbc1(j),j=1,5)
               else
                  write(9,21) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                  write(9,21) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                  write(9,21) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                  write(9,21) cbc1,ie,eface(1),(rbc1(j),j=1,5)
               endif
            enddo
            enddo
            enddo
         endif
      enddo
c
c     Scan through .rea file until end of bcs
c
      if (nfld.eq.1) call pfscan(string,'THERMAL BOUNDARY',16,8)
      if (nfld.ge.2) call pfscan(string,'PRESOLVE',8,8)
      lout = ltrunc(string1,80)
      write (9,81) (string1(j),j=1,lout)
   81 format(80a1)
c
c     Scan through and output .rea file until end of file
c
      call scanout(string,'xxxx',4,8,9)
c
c
c     generate  box.map       6/15/00
c
c     This is currently broken, at least for periodic domains, it seems.
c
c     call out_map(nelx,nely,nelz,if3d)
c
      stop
      end
c-----------------------------------------------------------------------
      subroutine scanout(string,input,len,infile,outfile)
c
c     scan through infile until input is found and output
c     all lines save that containing input to "outfile"
c
      character*80 string
      character*1 input(1)
      integer infile,outfile
c
      character*1 string1(80)
c
      do line=1,100000
         call blank(string,80)
         read (infile ,80,end=100,err=100) string
         if (indx1(string,input,len).ne.0) return
         call ccopy(string1,string,80)
         lout = ltrunc(string1,80)
         write (outfile,81) (string1(j),j=1,lout)
      enddo
   80 format(a80)
   81 format(80a1)
c
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine pfscan(string,input,len,infile)
c
c     scan through infile until input is found 
c
      character*80 string
      character*1 input(1)
      integer infile,outfile
c
      do line=1,100000
         call blank(string,80)
         read (infile ,80,end=100,err=100) string
         if (indx1(string,input,len).ne.0) return
      enddo
   80 format(a80)
c
  100 continue
      return
      end
c-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION LTRUNC(STRING,L)
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
      DATA BLNK/' '/
C
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      RETURN
      END
      INTEGER FUNCTION INDX1(S1,S2,L2)
      CHARACTER*80 S1,S2
C
      N1=80-L2+1
      INDX1=0
      IF (N1.LT.1) RETURN
C
      DO 100 I=1,N1
         I2=I+L2-1
         IF (S1(I:I2).EQ.S2(1:L2)) THEN
            INDX1=I
            RETURN
         ENDIF
  100 CONTINUE
C
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE CCOPY(A,B,N)
      CHARACTER*1 A(1), B(1)
      DO 100 I = 1, N
 100     A(I) = B(I)
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE RZERO(A,N)
      DIMENSION  A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      RETURN
      END
      FUNCTION MOD1(I,N)
C
C     Yields MOD(I,N) with the exception that if I=K*N, result is N.
C
      MOD1=0
      IF (I.EQ.0) RETURN
      IF (N.EQ.0) THEN
         WRITE(6,*)
     $  'WARNING:  Attempt to take MOD(I,0) in FUNCTION MOD1.'
         RETURN
      ENDIF
      II = I+N-1
      MOD1 = MOD(II,N)+1
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine scannocom(iend,infile)
c
c     scan through infile until "no comment" is found
c
      character*80 string
      character*1  string1(80)
      equivalence (string1,string)
c
      character*1 comment
      save        comment
      data        comment /'#'/
c
      iend = 0
      do line=1,100000
         call blank(string,80)
         read (infile ,80,end=100,err=100) string
c
         write(6,80) string
         if   (indx1(string,comment,1).ne.1) then
              open(unit=99,file='box.tmp')
              len = ltrunc(string,80)
              write(99,81) (string1(k),k=1,len)
              close(unit=99)
              return
         else
              i1 = indx1(string,comment,1)
c             write(6,*) i1,' i1', comment
         endif
c
      enddo
   80 format(a80)
   81 format(80a1)
c
  100 continue
      iend = 1
      return
      end
c-----------------------------------------------------------------------
      subroutine geti1(i,iend,io)
c
c     Get an integer from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine geti2(i1,i2,iend,io)
c
c     Get an integer from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i1,i2
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine geti3(i1,i2,i3,iend,io)
c
c     Get an integer from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) i1,i2,i3
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getr2(r1,r2,iend,io)
c
c     Get two reals from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) r1,r2
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getr3(r1,r2,r3,iend,io)
c
c     Get two reals from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,*) r1,r2,r3
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine getrv(r,n,iend,io)
      real r(n)
      parameter (big=1.e29)
c
c     Get reals from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
c
c     serious hack for f90
      call cfill(r,big,n)
c
      open(unit=99,file='box.tmp')
      read(99,*,end=1) (r(k),k=1,n)
      close(unit=99)
      return
    1 continue
      close(unit=99)
      do k=1,n
         if (r(k).eq.big) goto 15
      enddo
   15 continue
      write(6,*) 'this is k:',k,n
      read( 7,*,end=2) (r(j),j=k,n)
      return
    2 continue
      iend=1
      return
      end
c-----------------------------------------------------------------------
      subroutine getiv(r,n,iend,io)
      integer r(n)
      parameter (ibig=99999999)
c
c     Get integers from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
c
c     serious hack for f90
      call ifill(r,ibig,n)
c
      open(unit=99,file='box.tmp')
      read(99,*,end=1) (r(k),k=1,n)
      close(unit=99)
      return
    1 continue
      close(unit=99)
      do k=1,n
         if (r(k).eq.ibig) goto 15
      enddo
   15 continue
      write(6,*) 'this is k:',k,n
      read( 7,*,end=2) (r(j),j=k,n)
      return
    2 continue
      iend=1
      return
      end
c-----------------------------------------------------------------------
      subroutine getcv(c,m,n,iend,io)
      character*1 c(m,n)
      character*1 adum(80)
c
c     Get character strings, with single seperator, from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,1,end=2) (adum(k),k=1,80)
    1 format(80a1)
    2 continue
c
      i = 0
      do l=1,n
         do k=1,m
            i = i+1
            c(k,l) = adum(i)
         enddo
c
c        bump pointer for "," in input string
         i = i+1
      enddo
c
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine gets(c,n,iend,io)
      character*1 c(n)
c
c     Get character string from first uncommented line
c
      call scannocom(iend,io)
      if (iend.ne.0) return
      open(unit=99,file='box.tmp')
      read(99,1) (c(k),k=1,n)
    1 format(80a1)
      close(unit=99)
      return
      end
c-----------------------------------------------------------------------
      subroutine ifill(x,ic,n)
      integer x(1)
      do i=1,n
         x(i) = ic
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine cfill(x,c,n)
      real x(1)
      do i=1,n
         x(i) = c
      enddo
      return
      end
c-----------------------------------------------------------------------
      function ilsum(x,n)
      integer x(1)
      ilsum = 0
      do i=1,n
         ilsum = ilsum + x(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine out_map(nelx,nely,nelz,if3d)
c
      logical if3d
      integer ivv(8)
c
      open (unit=21,file='box.map',status='unknown')
c
      ndim=2
      if (if3d) ndim=3
      nel = nelx*nely*nelz
      nlv = 2**ndim
c
      write(21,*) nel,nelx,nely,nelz,ndim
c
c     output dual of lexicographical element ordering
c
      if (if3d) then
         ie = 0
         do jz=0,nelz-1
         do jy=0,nely-1
         do jx=0,nelx-1
            ie = ie+1
            l = 0
            do iz=0,1
            do iy=0,1
            do ix=0,1
               l = l+1
               ivv(l)=1+ix+2*iy+4*iz+8*jx+8*nelx*jy+8*nelx*nely*jz
            enddo
            enddo
            enddo
            write(21,2) ie,(ivv(k),k=1,8)
         enddo
         enddo
         enddo
      else
         ie = 0
         do jy=0,nely-1
         do jx=0,nelx-1
            ie = ie+1
            l = 0
            do iy=0,1
            do ix=0,1
               l = l+1
               ivv(l)=1+ix+jx+(nelx+1)*(iy+jy)
            enddo
            enddo
            write(21,2) ie,(ivv(k),k=1,4)
         enddo
         enddo
      endif
   2  format(9i9)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_multi_seg(nelxyz,x,y,z,m,if3d)
c
      integer nelxyz(3)
      real x(0:m),y(0:m),z(0:m)
      logical if3d
c
      integer nels (1000)
      real    gains(1000)
      real    xs   (0:1000)
c
c     This routine allows the user to input a complex sequence of
c     segments for each of the x, y and z directions, so that
c     genbox will generate the desired tensor product array.
c
c     Input format for each coordinat is
c
c     nsegs
c     nel_1    nel_2 ...  nel_nsegs
c     gain_1  gain_2 ... gain_nsegs
c     x_0, x_1, x_2, ...    x_nsegs
c
c     The output will be:
c
c     x(0)...x(nel_1)...x_(nel_2)...x(nel_segs),
c
c     where the segment between x(e-1) and x(e) is filled with
c     a distribution determined by gain_e.   Uniform spacing
c     corresponds to gain_e = 1, otherwise, a geometric sequence
c     is generated, with dx_i+1 = dx_i * gain_e
c
      ndim = 2
      if (if3d) ndim=3
c
      do id=1,ndim
c
         call geti1(nseg,iend,7)
         call getiv(nels ,nseg  ,iend,7)
         call getrv(xs   ,nseg+1,iend,7)
         call getrv(gains,nseg  ,iend,7)
         write(6,*) 'NSEG: ',id,nseg,(nels(k),k=1,nseg)
c
         nelxyz(id) = ilsum(nels,nseg)
         if (nelxyz(id).gt.m) then
            write(6,*) 'NEL too large:',nelxyz(id),m,id
            write(6,*) 'Abort in get_multi_seg.'
            call exit
         endif
c
         k = 0
         do is=1,nseg
            nel = nels(is)
            x0  = xs(is-1)
            x1  = xs(is)
            if (id.eq.1) call geometric_x(x(k),nel,x0,x1,gains(is))
            if (id.eq.2) call geometric_x(y(k),nel,x0,x1,gains(is))
            if (id.eq.3) call geometric_x(z(k),nel,x0,x1,gains(is))
            k = k + nel
         enddo
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine geometric_x(x,n,x0,x1,gain)
c
      real x(0:n)
c
      x(0) = 0.
      dx = (x1-x0)/n
      do i=1,n
         x(i) = x(i-1) + dx
         dx   = dx*gain
      enddo
c
      scale = (x1-x0)/(x(n)-x(0))
      do i=0,n
         x(i) = x0 + scale*x(i)
      enddo
c
      return
      end
c-----------------------------------------------------------------------
