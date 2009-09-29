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
      logical     if3d,ifevenx,ifeveny,ifevenz
      character*3 cbc1,   cbc2,   cbc3,   cbc4,   cbc5,   cbc6
      real        rbc1(5),rbc2(5),rbc3(5),rbc4(5),rbc5(5),rbc6(5)
      integer eface(6)
      save    eface
      data    eface  / 4,2,1,3,5,6 /
c
      character*1  apt(52)
      character*52 apt52
      equivalence (apt52,apt)
      data apt52/'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
      parameter (mbox  = 300)
      parameter (maxx  = 300)
      parameter (maxel = mbox*maxx*maxx)
      integer nlx(mbox),nly(mbox),nlz(mbox)
      real x(0:maxx,mbox),y(0:maxx,mbox),z(0:maxx,mbox)
      real xc(mbox),yc(mbox),zc(mbox)
      real curve(8,maxel)
      character*3 cbc(6,mbox,20)
      character*1 boxcirc(mbox)

      real*4 buf(30)
      character*80 hdr
      real*4 test
      data   test  / 6.54321 /
      logical iffo
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
      call scanout(string,'MESH DATA',9,8,9)
      read(8,*) i,i,i
      write(9,*) ' *** MESH DATA ***'
      
c-----------------------------------------------------------------------
c     here is where the 2d/3d determination is made....
      if3d = .false.
      iffo = .true.
      call geti1(ndim,iend,7)
      if(ndim.lt.0) then
        iffo = .false.
        call byte_open('box.re2\0')
      endif
      ndim = abs(ndim)
      if (ndim.eq.3) if3d = .true.
c-----------------------------------------------------------------------
c     here is where the fluid, fluid+heat determination is made....
      call geti1(nfld,iend,7)
c-----------------------------------------------------------------------
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
       
         if (if3d) then
            call geti3(nelx,nely,nelz,iend,7)
c
            if (nelx.lt.0) ifevenx = .true.
            if (nely.lt.0) ifeveny = .true.
            if (nelz.lt.0) ifevenz = .true.
            nelx = abs(nelx)
            nely = abs(nely)
            nelz = abs(nelz)

            if(nelx*nely.gt.maxel) then
              write(6,*) 'ABORT: increase maxel and recompile!'
              call exitt
            endif
c
            nery = nely
            if (boxcirc(ibox).eq.'C') nery = 1
            if (iend.eq.1) goto 99
            ninbox = nelx*nely*nelz
            write(6,6) ninbox,nelx,nely,nelz,ibox,nfld
c
            if (ifevenx) then
               call getr2(x1,x2,iend,7)
               dx = (x2-x1)/nelx
               do k=0,nelx
                  x(k,ibox) = x1 + dx*k
               enddo
            else
               call getrv(x(0,ibox),nelx+1,iend,7)
            endif
c
            if (ifeveny) then
               call getr2(y1,y2,iend,7)
               dy = (y2-y1)/nely
               do k=0,nely
                  y(k,ibox) = y1 + dy*k
               enddo
            else
               call getrv(y(0,ibox),nely+1,iend,7)
            endif
c
            if (ifevenz) then
               call getr2(z1,z2,iend,7)
               dz = (z2-z1)/nelz
               do k=0,nelz
                  z(k,ibox) = z1 + dz*k
               enddo
            else
               call getrv(z(0,ibox),nelz+1,iend,7)
            endif
c
            do ifld=1,nfld
               call getcv(cbc(1,ibox,ifld),3,6,iend,7)
            enddo
         else
            nelz = 1
            call geti2(nelx,nely,iend,7)
c
            if (nelx.lt.0) ifevenx = .true.
            if (nely.lt.0) ifeveny = .true.
            nelx = abs(nelx)
            nely = abs(nely)

            if(nelx*nely.gt.maxel) then
              write(6,*) 'ABORT: increase maxel and recompile!'
              call exitt
            endif
c
            nery = nely
            if (boxcirc(ibox).eq.'C') nery = 1
            if (iend.eq.1) goto 99
            ninbox = nelx*nely*nelz
            write(6,6) ninbox,nelx,nely,nery,ibox,nfld
c
            if (ifevenx) then
               call getr2(x1,x2,iend,7)
               dx = (x2-x1)/nelx
               do k=0,nelx
                  x(k,ibox) = x1 + dx*k
               enddo
            else
               call getrv(x(0,ibox),nelx+1,iend,7)
            endif
c            write(999,*) (x(i,1),i=0,nelx)

c
            if (ifeveny) then
               call getr2(y1,y2,iend,7)
               dy = (y2-y1)/nely
               do k=0,nely
                  y(k,ibox) = y1 + dy*k
               enddo
            else
               call getrv(y(0,ibox),nely+1,iend,7)
            endif
c            write(998,*) (y(i,1),i=0,nely)

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
    6 format('Reading',i8,' =',3i8,' elements for box',i4,'.')
   99 continue
c
      if(iffo) then
        write(6,*) 'Beginning construction of box.rea'
      else
        write(6,*) 'Beginning construction of box.rea/box.re2'
      endif
      write(6,*) nel,' elements will be created for ',nbox,' boxes.'
      write(6,*) 
c
c     Construct mesh
c
      if(iffo) then
        write(9,10) nel,ndim,nel
   10   format(3i10,'           NEL,NDIM,NELV')
      else
        write(9,10) -nel,ndim,nel
        call blank(hdr,80)
        write(hdr,111) nel,ndim,nel
  111   format('#v001',i9,i3,i9,' this is the hdr')
        call byte_write(hdr,20)   ! assumes byte_open() already issued
        call byte_write(test,1)   ! write the endian discriminator
      endif
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
                  if(iffo) then
                    write(9,11) ie,ilev,apt(ia),'0'
c
                    write(9,12) x1,x2,x2,x1
                    write(9,12) y1,y1,y2,y2
                    write(9,12) z1,z1,z1,z1
c
                    write(9,12) x1,x2,x2,x1
                    write(9,12) y1,y1,y2,y2
                    write(9,12) z2,z2,z2,z2
                  else
                    igroup = 0
                    call byte_write(igroup, 1)
                    buf(1)  = x1
                    buf(2)  = x2
                    buf(3)  = x2
                    buf(4)  = x1
                    buf(5)  = x1
                    buf(6)  = x2
                    buf(7)  = x2
                    buf(8)  = x1
                    buf(9)  = y1
                    buf(10) = y1
                    buf(11) = y2
                    buf(12) = y2
                    buf(13) = y1
                    buf(14) = y1
                    buf(15) = y2
                    buf(16) = y2
                    buf(17) = z1
                    buf(18) = z1
                    buf(19) = z1
                    buf(20) = z1
                    buf(21) = z2
                    buf(22) = z2
                    buf(23) = z2
                    buf(24) = z2
                    call byte_write(buf,24)
                  endif
               enddo
            enddo
          enddo
        enddo
   12   format(4g14.6)
   11   format(
     $'            ELEMENT',i8,' [',i5,a1,']  GROUP  ',a1)
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
                    if(iffo) then
                      write(9,11) ie,ilev,apt(ia),'0'
                      write(9,12) x1,x2,x3,x4
                      write(9,12) y1,y2,y3,y4
                    else
                      igroup = 0
                      call byte_write(igroup, 1)
                      buf(1) = x1
                      buf(2) = x2
                      buf(3) = x3
                      buf(4) = x4
                      buf(5) = y1
                      buf(6) = y2
                      buf(7) = y3
                      buf(8) = y4
                      call byte_write(buf,8)
                    endif
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
                    if(iffo) then
                      write(9,11) ie,ilev,apt(ia),'0'
                      write(9,12) x1,x2,x3,x4
                      write(9,12) y1,y2,y3,y4
                    else
                      igroup = 0
                      call byte_write(igroup, 1)
                      buf(1) = x1
                      buf(2) = x2
                      buf(3) = x3
                      buf(4) = x4
                      buf(5) = y1
                      buf(6) = y2
                      buf(7) = y3
                      buf(8) = y4
                      call byte_write(buf,8)
                    endif
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
                    if(iffo) then
                      write(9,11) ie,ilev,apt(ia),'0'
                      write(9,12) x1,x2,x3,x4
                      write(9,12) y1,y2,y3,y4
                    else
                      igroup = 0
                      call byte_write(igroup, 1)
                      buf(1) = x1
                      buf(2) = x2
                      buf(3) = x3
                      buf(4) = x4
                      buf(5) = y1
                      buf(6) = y2
                      buf(7) = y3
                      buf(8) = y4
                      call byte_write(buf,8)
                    endif
                 enddo
              enddo
           endif
        enddo
      endif
c
c     output curve stuff and Boundary conditions
      if(iffo) then
        write(9,28) ncurv
   28   format(
     $   '  ***** CURVED SIDE DATA *****',/,
     $      i6,' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE')

        zero = 0.
        if (nel.lt.1000) then
           do ie=1,nel
              do iedge = 2,4,2
                 if (curve(iedge,ie).ne.0) write(9,290) 
     $            iedge,ie,curve(iedge,ie),(zero,k=1,4),'C'
              enddo
           enddo
  290      format(i3,i3,5g14.6,1x,a1)
         else
           do ie=1,nel
              do iedge = 2,4,2
                 if (curve(iedge,ie).ne.0) write(9,291) 
     $              iedge,ie,curve(iedge,ie),(zero,k=1,4),'C'
              enddo
           enddo
  291      format(i2,i6,5g14.6,1x,a1)
         endif
      else
         call byte_write(ncurv,1)  
         do ie=1,nel
            do iedge = 2,4,2
               if (curve(iedge,ie).ne.0) then
                  if(iffo) then
                    write(9,291) 
     $                iedge,ie,curve(iedge,ie),(zero,k=1,4),'C'
                  else
                      buf(1) = ie
                      buf(2) = iedge
                      buf(3) = curve(iedge,ie)
                      buf(4) = zero
                      buf(5) = zero
                      buf(6) = zero
                      buf(7) = zero
                      buf(8) = 'C'
                      call byte_write(buf,8)
                  endif
               endif
            enddo
         enddo
      endif

c
      if(iffo) 
     &   write(9,31) 
   31    format('  ***** BOUNDARY CONDITIONS *****')
c
      do ifld=1,nfld
c
         if2 = ifld-2
         if(iffo) then
           if (ifld.eq.1) write(9,41)
           if (ifld.eq.2) write(9,42)
           if (ifld.ge.3) write(9,43) if2
   41      format('  ***** FLUID   BOUNDARY CONDITIONS *****')
   42      format('  ***** THERMAL BOUNDARY CONDITIONS *****')
   43      format('  ***** PASSIVE SCALAR',i2
     $                         ,'  BOUNDARY CONDITIONS *****')
         endif
c
         nbc = 0
         if (if3d) then
            do ipass=1,2
               ie = 0
               if(ipass.eq.2 .and. .not. iffo) then
                 call byte_write(nbc,1)
               endif
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
                  if (cbc1.ne.'P  ')  call rzero(rbc1,5)
                  nbc = nbc + 1
               else
                  cbc1='E  '
               endif
c
               if (iex.eq.nelx) then
                  cbc2=cbc(2,ibx,ifld)
                  if (cbc2.ne.'P  ') call rzero(rbc2,5)
                  nbc = nbc + 1
               else
                  cbc2='E  '
               endif
c
               if (iey.eq.1) then
                  cbc3=cbc(3,ibx,ifld)
                  if (cbc3.ne.'P  ') call rzero(rbc3,5)
                  nbc = nbc + 1
               else
                  cbc3='E  '
               endif
c
               if (iey.eq.nely) then
                  cbc4=cbc(4,ibx,ifld)
                  if (cbc4.ne.'P  ') call rzero(rbc4,5)
                  nbc = nbc + 1
               else
                  cbc4='E  '
               endif
c
               if (iez.eq.1) then
                  cbc5=cbc(5,ibx,ifld)
                  if (cbc5.ne.'P  ') call rzero(rbc5,5)
                  nbc = nbc + 1
               else
                  cbc5='E  '
               endif
c
               if (iez.eq.nelz) then
                  cbc6=cbc(6,ibx,ifld)
                  if (cbc6.ne.'P  ') call rzero(rbc6,5)
                  nbc = nbc + 1
               else
                  cbc6='E  '
               endif
c
c              output bc's in preproc. notation
c
               ie = ie+1
               if(iffo .and. ipass.eq.1) then
                 if (nel.lt.1000) then
                    write(9,20) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                    write(9,20) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                    write(9,20) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                    write(9,20) cbc1,ie,eface(1),(rbc1(j),j=1,5)
                    write(9,20) cbc5,ie,eface(5),(rbc5(j),j=1,5)
                    write(9,20) cbc6,ie,eface(6),(rbc6(j),j=1,5)
                 elseif (nel.lt.100000) then
                    write(9,21) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                    write(9,21) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                    write(9,21) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                    write(9,21) cbc1,ie,eface(1),(rbc1(j),j=1,5)
                    write(9,21) cbc5,ie,eface(5),(rbc5(j),j=1,5)
                    write(9,21) cbc6,ie,eface(6),(rbc6(j),j=1,5)
                 elseif (nel.lt.10000000) then
                    write(9,22) cbc3,ie,eface(3),(rbc3(j),j=1,5)
                    write(9,22) cbc2,ie,eface(2),(rbc2(j),j=1,5)
                    write(9,22) cbc4,ie,eface(4),(rbc4(j),j=1,5)
                    write(9,22) cbc1,ie,eface(1),(rbc1(j),j=1,5)
                    write(9,22) cbc5,ie,eface(5),(rbc5(j),j=1,5)
                    write(9,22) cbc6,ie,eface(6),(rbc6(j),j=1,5)
                 endif
   20            format(1x,a3,2i3,5g14.7)
   21            format(1x,a3,i5,i1,5g14.7)
   22            format(1x,a3,i6,i1,5g14.7)
               elseif (.not. iffo) then
                 if(ipass.eq.2) then
                   call blank(buf(1),30*4)
                   call icopy(buf(1),ie,1)
c                   call blank(buf(8),4)
 
                   if(cbc3.ne.'E  ') then 
                     call icopy(buf(2),eface(3),1)
                     call copy(buf(3),rbc3,5)
                     call chcopy(buf(8),cbc3,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc2.ne.'E  ') then 
                     call icopy(buf(2),eface(2),1)
                     call copy(buf(3),rbc2,5)
                     call chcopy(buf(8),cbc2,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc4.ne.'E  ') then 
                     call icopy(buf(2),eface(4),1)
                     call copy(buf(3),rbc4,5)
                     call chcopy(buf(8),cbc4,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc1.ne.'E  ') then 
                     call icopy(buf(2),eface(1),1)
                     call copy(buf(3),rbc1,5)
                     call chcopy(buf(8),cbc1,3)
                     call byte_write(buf,8)
                   endif 
 
                   if(cbc5.ne.'E  ') then 
                     call icopy(buf(2),eface(5),1)
                     call copy(buf(3),rbc5,5)
                     call chcopy(buf(8),cbc5,3)
                     call byte_write(buf,8)
                   endif

                   if(cbc6.ne.'E  ') then 
                     call icopy(buf(2),eface(6),1)
                     call copy(buf(3),rbc6,5)
                     call chcopy(buf(8),cbc6,3)
                     call byte_write(buf,8)
                   endif
                 endif
               endif
            enddo
            enddo
            enddo
            enddo
            enddo
         else
            do ipass=1,2
               ie = 0
               if(ipass.eq.2 .and. .not. iffo) then
                 call byte_write(nbc,1)
               endif
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
                  nbc = nbc + 1
               else
                  cbc1='E  '
               endif
c
               if (iex.eq.nelx) then
                  cbc2=cbc(2,ibx,ifld)
                  if (cbc2.ne.'P  ') call rzero(rbc2,5)
                  nbc = nbc + 1
               else
                  cbc2='E  '
               endif
c
               if (iey.eq.1) then
                  cbc3=cbc(3,ibx,ifld)
                  if (cbc3.ne.'P  ') call rzero(rbc3,5)
                  nbc = nbc + 1
               else
                  cbc3='E  '
               endif
c
               if (iey.eq.nely) then
                  cbc4=cbc(4,ibx,ifld)
                  if (cbc4.ne.'P  ') call rzero(rbc4,5)
                  nbc = nbc + 1
               else
                  cbc4='E  '
               endif
c
c              output bc's in preproc. notation
c
               ie = ie+1
               if(iffo .and. ipass.eq.1) then
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
               elseif (.not. iffo) then
                 if(ipass.eq.2) then
                   call icopy(buf(1),ie,1)
                   call blank(buf(8),4)
 
                   if(cbc3.ne.'E  ') then 
                     call icopy(buf(2),eface(3),1)
                     call copy(buf(3),rbc3,5)
                     call chcopy(buf(8),cbc3,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc2.ne.'E  ') then 
                     call icopy(buf(2),eface(2),1)
                     call copy(buf(3),rbc2,5)
                     call chcopy(buf(8),cbc2,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc4.ne.'E  ') then 
                     call icopy(buf(2),eface(4),1)
                     call copy(buf(3),rbc4,5)
                     call chcopy(buf(8),cbc4,3)
                     call byte_write(buf,8)
                   endif
 
                   if(cbc1.ne.'E  ') then 
                     call icopy(buf(2),eface(1),1)
                     call copy(buf(3),rbc1,5)
                     call chcopy(buf(8),cbc1,3)
                     call byte_write(buf,8)
                   endif 
                 endif
               endif
            enddo
            enddo
            enddo
            enddo
         endif
      enddo
c
c     Scan through .rea file until end of bcs
c
      if(iffo) then
        if (nfld.eq.1) call scan(string,'THERMAL BOUNDARY',16,8)
        if (nfld.ge.2) call scan(string,'PRESOLVE',8,8)
        lout = ltrunc(string1,80)
        write (9,81) (string1(j),j=1,lout)
   81   format(80a1)
      else
        call byte_close()
      endif
c
c     Scan through and output .rea file until end of file
c
      call scanout(string,'xxxx',4,8,9)
c
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
      subroutine scan(string,input,len,infile)
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
         if   (indx1(string,comment,1).ne.1) then
c              write(*,*) line
c              write(6,80) string
              open(unit=99,file='box.tmp')
              len = ltrunc(string,80)
              write(99,81) (string1(k),k=1,len)
              close(unit=99)
              return
         else
              i1 = indx1(string,comment,1)
c              write(6,*) i1,' i1', comment
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
      subroutine getrv(r,n,iend,io)
      real r(n)
c
c     Get reals from first uncommented line
c
      istart=1
      call scannocom(iend,io)
      if (iend.ne.0) return

      call cfill(r,1e+30,n)
 222  open(unit=99,file='box.tmp')
      read(99,*,end=1,iostat=istat) (r(k),k=istart,n)
   1  close(unit=99)

      ! check how many points we have already read
      icounter = 0
      do i=1,n
         if(r(i).ne.1e+30) icounter = icounter + 1
      enddo
 
      if(icounter.lt.n) then
        write(*,*) 'error:', icounter,n
        call scannocom(iend,io)
        istart=icounter+1
        goto 222
      endif


      return
c    1 continue
c      close(unit=99)
c      write(6,*) 'this is k:',k,n
c      read( 7,*,end=2) (r(j),j=k,n)
c      return
c    2 continue
c      iend=1
c      return
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
      subroutine cfill(a,b,n)
      DIMENSION  A(1)
C
      DO I = 1, N
        A(I) = B
      ENDDO
      
      return
      end
c-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      integer a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      character*1 a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1), b(1)
      do 100 i = 1, n
 100     a(i) = b(i)
      return
      end
c-----------------------------------------------------------------------
      subroutine exitt

      call exit 
 
      return
      end
