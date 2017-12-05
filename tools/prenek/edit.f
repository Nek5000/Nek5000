c-----------------------------------------------------------------------
      subroutine mesh_edit   ! enter 2D mesh editor
#     include "basics.inc"

      integer e
      logical ifdomdef
      save    ifdomdef
      data    ifdomdef / .false. /

      integer nelo
      save    nelo
      data    nelo / 0 /

      common /domainr/ xdom(2,100),ddom(40*nelm)
      common /domaini/ npdom,ndom,ein(5*nelm),nin
      integer ein

    1 continue

      if (nel.ne.nelo) ifdomdef = .false.
      nelo = nel

      item(1) ='UP MENU'
      item(2) ='HELP'
      item(3) ='Define Domain'
      nchoic = 3

      if (ifdomdef) then ! domain has been defined
         call dom_draw_bdry(npdom,xdom)
         item(4) ='Move Vertex'
         item(5) ='Edge Spline'
         item(6) ='Delete Elements'
         nchoic  = 6
      endif

      call menu(xmouse,ymouse,button,'Edit Mesh')  ! Prompt for input

      if (choice.eq.'UP MENU')          return
      if (choice.eq.'HELP')             call prs('No help yet.$')
      if (choice.eq.'Define Domain')    call dom_def(ifdomdef)
      if (choice.eq.'Move Vertex')      call dom_mv_vtx
     $                                       (ein,nin,ddom,xdom,npdom)
      if (choice.eq.'Edge Spline')      call dom_spline
      if (choice.eq.'Delete Elements')  call dom_delete
     $                                       (ein,nin,ddom,xdom,npdom)

      call refresh
      call drgrid
      do e=1,nel
         call drawel(e)
      enddo
      if (ifdomdef) call dom_draw_bdry(npdom,xdom)

      goto 1

      return
      end
c-----------------------------------------------------------------------
      subroutine dom_spline
      return
      end
c-----------------------------------------------------------------------
      subroutine dom_draw_bdry(npdom,xdom)
      real xdom(2,npdom)

      call move(xdom(1,1),xdom(2,1))
      do k=2,npdom
         call prsrr('BDRY:$',xdom(1,k),xdom(2,k))
         call draw(xdom(1,k),xdom(2,k))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine dom_mv_vtx(ein,nin,ddom,xdom,npdom)

#     include "basics.inc"

      integer e,ein(nin)
      real ddom(npdom-1,nin),xdom(2,npdom)

      call prsi('Select domain vertex to move.  nin=$',nin)
      call mouse2(xmouse,ymouse,button,.true.)

      call prs('Select new vertex position (rt attch to el.)$.')
      call mouse2(xmove,ymove,button,.true.)

      rmin = 1.e20
      do i=1,npdom-1
         r = (xdom(1,i)-xmouse)**2 + (xdom(2,i)-ymouse)**2
         if (r.lt.rmin) then
            imin = i
            rmin = r
         endif
      enddo

      dx = xmove-xdom(1,imin)   ! displacement
      dy = ymove-xdom(2,imin) 
      write(6,*) 'this is dx:',dx,dy,imin

      do k=1,nin
         i = mod(ein(k),5)
         e = 1 + ein(k)/5
         write(6,9) 'i,e:',i,e,k,ein(k),x(i,e),y(i,e),dx,dy,ddom(imin,k)
  9      format(a4,4i4,5f11.4)
         if (i.gt.0) then ! move vtx (i,e)
            x(i,e) = x(i,e) + ddom(imin,k)*dx
            y(i,e) = y(i,e) + ddom(imin,k)*dy
         endif
      enddo

      xdom(1,imin) = xmove
      xdom(2,imin) = ymove

      if (imin.eq.1) xdom(1,npdom) = xmove
      if (imin.eq.1) xdom(2,npdom) = ymove

      return
      end
c-----------------------------------------------------------------------
      subroutine dom_def(ifdomdef) ! grab a polygon
#     include "basics.inc"

      integer e
      logical ifdomdef
      common /domainr/ xdom(2,100),ddom(40*nelm)
      common /domaini/ npdom,ndom,ein(5*nelm),nin
      integer ein

      call prs('Click points defining convex domain; menu to close$.')

      k = 0
      do kk=1,100
         call mouse2(xmouse,ymouse,button,.true.)
         if (xscr(xmouse).le.xlmen) then

            k = k+1
c           write(line,10) k,xmouse,ymouse
   10       format('Point:',i3,1p2e12.4,'$')
            call prs(line)

            xdom(1,k) = xmouse                 ! std entry
            xdom(2,k) = ymouse

            if (k.eq.1) call move(xmouse,ymouse)
            if (k.gt.1) call draw(xmouse,ymouse)

         else
           if (k.lt.2) then
              call prs(' You have < 3 points. Delete domain? (y/n)$')
              call res(ans,1)
              if (ans.eq.'y' .or. ans.eq.'Y') then
                 ifdomdef = .false.
                 return
              endif
              call prs('Click to define domain, menu area to close$.')
           endif
           goto 20
         endif
      enddo
   20 continue

      ifdomdef = .true.
      npdom = k+1
      xdom(1,npdom) = xdom(1,1)
      xdom(2,npdom) = xdom(2,1)
      call draw(xdom(1,npdom),xdom(2,npdom))

      call dom_el_inside(ein,nin,ddom,xdom,npdom)

c     if (nin.gt.5*nelm .or. nin*(npdom-1).gt.40*nelm) then
c        call prsiii('Too many vertices:$',nin,npdom,nelm)
c        nin = 10
c        ifdomdef = .false.
c     endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dom_el_inside(ein,nin,ddom,xdom,npdom)

#     include "basics.inc"

      integer e,ein(1)
      real ddom(npdom-1,1),xdom(2,npdom)
      real ddist(0:100)
      real dnorm(0:100)

      call prs('Click points to define domain, menu area to close$.')

      blx = xdom(1,1)
      bux = xdom(1,1)
      bly = xdom(2,1)
      buy = xdom(2,1)

      do i=1,npdom                 ! bounding box
         blx = min(blx,xdom(1,i))
         bly = min(bly,xdom(2,i))
         bux = max(bux,xdom(1,i))
         buy = max(buy,xdom(2,i))
      enddo

      dlx = bux-blx
      dly = buy-bly
      blx = blx-.001*dlx
      bly = bly-.001*dly
      bux = bux+.001*dlx
      buy = buy+.001*dly

c     write(6,*) 'bb:',blx,bux,bly,buy

      epsm = -1.e-4*max(dlx,dly)
      nin = 0
      do e=1,nel
         xce = 0.
         yce = 0.
         do i=1,4
            xce = xce + x(i,e)
            yce = yce + y(i,e)
            if (blx.le.x(i,e).and.x(i,e).le.bux .and.
     $          bly.le.y(i,e).and.y(i,e).le.buy ) then ! inside bbox
                ddist(0) = dom_dist(ddist(1),xdom,npdom,x(i,e),y(i,e))
                if (ddist(0).ge.epsm) then
                   nin = nin+1
                   ein(nin) = i+5*(e-1)
                   do k=1,npdom-1
                      ddom(k,nin) = max(ddist(k),0.)
c                     write(6,*) e,i,k,nin,ddom(k,nin),ein(nin),'dst'
                   enddo
                endif
            endif
c           write(6,*) 'NIN:',nin,i,e,epsm,ddist(0)
         enddo

         xce = .25*xce  ! check centroid
         yce = .25*yce

         if (blx.le.xce.and.xce.le.bux .and.
     $       bly.le.yce.and.yce.le.buy ) then ! inside bbox
             ddist(0) = dom_dist(ddist(1),xdom,npdom,xce,yce)
             if (ddist(0).ge.epsm) then
                nin = nin+1
                ein(nin) = 0+5*(e-1)
                do k=1,npdom-1
                   ddom(k,nin) = ddist(k)
                enddo
             endif
         endif
c        write(6,*) 'NIN:',nin,i,e,epsm,ddist(0)
      enddo
c     write(6,*) 'THIS IS nin:',nin

c     normalize all distances by distance fct at bdry

      do k=1,npdom-1
        ddist(0) = dom_dist(ddist(1),xdom,npdom,xdom(1,k),xdom(2,k))
        dnorm(k) = ddist(k)
c       write(6,*) 'DNORM:',k,dnorm(k)
        if (dnorm(k).ne.0.) dnorm(k) = 1./dnorm(k)
      enddo

      do i=1,nin
      do k=1,npdom-1
         ddom(k,i) = ddom(k,i)*dnorm(k)
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      function dom_dist(dlist,xdom,npdom,x,y)

      real dlist(npdom-1),xdom(2,npdom)
      real xx(2),xc(2),dtmp(100)

      xx(1) = x
      xx(2) = y

      xc(1) = 0.           ! Assumption - centroid inside domain
      xc(2) = 0.
      do k=1,npdom-1
         xc(1)=xc(1)+xdom(1,k)
         xc(2)=xc(2)+xdom(2,k)
      enddo
      xc(1) = xc(1)/(npdom-1)
      xc(2) = xc(2)/(npdom-1)

      dist = 1.e20
      do k=1,npdom-1
         dtmp(k) = dline(xx,xdom(1,k),xdom(1,k+1),xc) ! positive > inside
         dist = min(dist,dtmp(k))
      enddo

      do k=1,npdom-1
         dlist(k) = 1.e20
         k1 = k-1
         if (k1.eq.0) k1=npdom-1
         do j=1,npdom-1
            if (j.ne.k.and.j.ne.k1) dlist(k) = min(dlist(k),dtmp(j))
         enddo
      enddo

      dom_dist = dist

      return
      end
c-----------------------------------------------------------------------
      function dline(x,x1,x2,x0)

      real x(2),x1(2),x2(2),x0(2)
      real v(2),c(2),nx,ny,nn

      call sub3(v,x2,x1,2)
      call sub3(c,x0,x1,2)

      vxc = v(1)*c(2) - v(2)*c(1)

      nx  = -vxc*v(2)
      ny  =  vxc*v(1)
      nn  =  nx**2 + ny**2

      dline = 0.

      if (nn.gt.0) then
         nn  =  sqrt(nn)
         nx  =  nx/nn
         ny  =  ny/nn
         dline = nx*(x(1)-x1(1)) + ny*(x(2)-x1(2))
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine dom_delete(ein,nin,ddom,xdom,npdom)

#     include "basics.inc"
      common /deletei/ dflag(nelm)
      integer dflag

      integer ein(nin)
      real ddom(npdom-1,nin),xdom(2,npdom)

      integer e,emin,ecount

      call izero(dflag,nel)

      emin = nelm
      do k=1,nin
         i = mod(ein(k),5)
         e = 1 + ein(k)/5
         if (i.eq.0) dflag(e)=1
         if (i.eq.0) emin = min(e,emin)
      enddo

      ndel = 0
      do e=1,nel
         ndel = ndel+dflag(e)
      enddo
      if (ndel.eq.0) return  ! no centroids in domain


      ecount = 0
      do e=1,nel
         if (dflag(e).eq.0) then
            ecount = ecount + 1                     ! keep this one
            if (e.gt.ecount) call copyel(e,ecount)  ! e --> ecount
         endif
      enddo
      call prsis('Deleted$',ndel,' elements.$')
      nel  = ecount
      call prsis('Now have$',nel,' elements.$')
c
      return
      end
c-----------------------------------------------------------------------
