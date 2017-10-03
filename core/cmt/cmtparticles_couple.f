c----------------------------------------------------------------------
      subroutine spread_props_grid
c
c     spread particle properties at fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'MASS'
      include 'TSTEP'
      include 'CMTPART'

      real   xdrange(2,3)
      common /domainrange/ xdrange
      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      real    xx,yy,zz,vol,pfx,pfy,pfz,pmass,pmassf,vcell,spl,multfc
     >       ,qgqf,rvx,rvy,rvz,rcountv(8,nelt)
      integer e

      ! begin timer
      ptdum(6) = dnekclock()

      nlxyze = lx1*ly1*lz1*lelt
      nxyze  = nx1*ny1*nz1*nelt
      call rzero(ptw,nlxyze*8)
      call rzero(rcountv,8*nelt)

c     ! for grid spreading line in finite volume.. no ghost particles
c     local mpi rank effects
      do ip=1,n
         e   = ipart(je0,ip) + 1
         xx  = rpart(jx,ip)
         yy  = rpart(jy,ip)
         zz  = rpart(jz,ip)
         spl = rpart(jspl,ip)
         rgam= rpart(jgam,ip)

         rdum = spl

         pfx = -rpart(jf0,ip)*rdum
         pfy = -rpart(jf0+1,ip)*rdum
         pfz = -rpart(jf0+2,ip)*rdum
         vol = rpart(jvol,ip)*rdum
         qgqf= -(rpart(jg0,ip) + rpart(jq0,ip))*rdum
         rvx = rpart(jv0  ,ip)*vol
         rvy = rpart(jv0+1,ip)*vol
         rvz = rpart(jv0+2,ip)*vol

         rcountv(1,e) = rcountv(1,e) + pfx
         rcountv(2,e) = rcountv(2,e) + pfy
         rcountv(3,e) = rcountv(3,e) + pfz
         rcountv(4,e) = rcountv(4,e) + vol
         rcountv(5,e) = rcountv(5,e) + qgqf
         rcountv(6,e) = rcountv(6,e) + rvx
         rcountv(7,e) = rcountv(7,e) + rvy
         rcountv(8,e) = rcountv(8,e) + rvz
      enddo

         call local_part_to_grid_fv(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
     >                             ptw(1,1,1,1,3),ptw(1,1,1,1,4),
     >                             ptw(1,1,1,1,5),ptw(1,1,1,1,6),
     >                             ptw(1,1,1,1,7),ptw(1,1,1,1,8),
     >                             rcountv)


      pi       = 4.0d+0*atan(1.0d+0)
      rbexpi   = 1./(-2.*rsig**2)

c     use if exponential mollification
      ralph    = d2chk(1)     ! assume all directions same!
      multfc   = 1./(sqrt(2.*pi)**3 * rsig**3) ! exponential
c     end use

c     use if uniform mollification
c     ralph    = d2chk(1)
c     multfc   = 1./(2.*ralph)**3              ! uniform
c     end use

      ralph2   = ralph**2



c     ! for rectangular grid spreading of gaussian.. need to change gen
c     local mpi rank effects
c     do ip=1,n
c        e   = ipart(je0,ip) + 1
c        xx  = rpart(jx,ip)
c        yy  = rpart(jy,ip)
c        zz  = rpart(jz,ip)
c        spl = rpart(jspl,ip)
c        rgam= rpart(jgam,ip)

c        rdum = multfc*spl*rgam
c        rdum1= multfc*rgam
c        rdum = multfc*spl

c        pfx = -rpart(jf0,ip)*rdum
c        pfy = -rpart(jf0+1,ip)*rdum
c        pfz = -rpart(jf0+2,ip)*rdum
c        vol = rpart(jvol,ip)*rdum
c        qgqf= -(rpart(jg0,ip) + rpart(jq0,ip))*rdum
c        rvx = rpart(jv0  ,ip)*vol
c        rvy = rpart(jv0+1,ip)*vol
c        rvz = rpart(jv0+2,ip)*vol

c        call local_part_to_grid(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
c    >                             ptw(1,1,1,1,3),ptw(1,1,1,1,4),
c    >                             ptw(1,1,1,1,5),ptw(1,1,1,1,6),
c    >                             ptw(1,1,1,1,7),ptw(1,1,1,1,8),
c    >                             pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,
c    >                             xx,yy,zz,rbexpi,
c    >                             ralph,ralph2,e)
c     enddo

c     ntmp1 = iglmax(n,1)
c     ntmp2 = iglmax(nfptsgp,1)
c     if (nid.eq.0) write(6,*) 'Passed local spreading to grid'
c    >                                                     ,ntmp1,ntmp2

c     remote mpi rank effects
c     do ip=1,nfptsgp
c        e   = iptsgp(jgpes,ip) + 1
c        xx  = rptsgp(jgpx,ip)
c        yy  = rptsgp(jgpy,ip)
c        zz  = rptsgp(jgpz,ip)
c        spl = rptsgp(jgpspl,ip)
c        rgam= rptsgp(jgpgam,ip)

c        rdum = multfc*spl*rgam
c        rdum1= multfc*rgam
c        rdum = multfc*spl

c        pfx = -rptsgp(jgpfh,ip)*rdum
c        pfy = -rptsgp(jgpfh+1,ip)*rdum
c        pfz = -rptsgp(jgpfh+2,ip)*rdum
c        vol = rptsgp(jgpvol,ip)*rdum
c        qgqf= -(rptsgp(jgpg0,ip) + rptsgp(jgpq0,ip))*rdum
c        rvx = rptsgp(jgpv0  ,ip)*vol
c        rvy = rptsgp(jgpv0+1,ip)*vol
c        rvz = rptsgp(jgpv0+2,ip)*vol

c        call remote_part_to_grid(ptw(1,1,1,1,1),ptw(1,1,1,1,2),
c    >                             ptw(1,1,1,1,3),ptw(1,1,1,1,4),
c    >                             ptw(1,1,1,1,5),ptw(1,1,1,1,6),
c    >                             ptw(1,1,1,1,7),ptw(1,1,1,1,8),
c    >                             pfx,pfy,pfz,vol,qgqf,rvx,rvy,rvz,
c    >                             xx,yy,zz,rbexpi,
c    >                             ralph,ralph2,e)
c     enddo

c     ntmp = iglsum(n,1)
c     if (nid.eq.0) write(6,*) 'Passed remote spreading to grid'

      rvfmax = 0.7
      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         if (ptw(i,j,k,ie,4) .gt. rvfmax) ptw(i,j,k,ie,4) = rvfmax
         phig(i,j,k,ie) = 1. - ptw(i,j,k,ie,4)
      enddo
      enddo
      enddo
      enddo

      ! end timer
      pttime(6) = pttime(6) + dnekclock() - ptdum(6)

      return
      end
c----------------------------------------------------------------------
      subroutine local_part_to_grid_fv(fvalgx,fvalgy,fvalgz,fvalgv,
     >                              fvalgg,fvalv1,fvalv2,fvalv3,
     >                              rcountv)
c
c     spread a local particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      integer e,er
      real    fvalgx(nx1,ny1,nz1,nelt),fvalgy(nx1,ny1,nz1,nelt),
     >        fvalgz(nx1,ny1,nz1,nelt),fvalgv(nx1,ny1,nz1,nelt),
     >        fvalgg(nx1,ny1,nz1,nelt),fvalv1(nx1,ny1,nz1,nelt),
     >        fvalv2(nx1,ny1,nz1,nelt),fvalv3(nx1,ny1,nz1,nelt),
     >        rcountv(8,nelt)

      ! begin timer
      ptdum(7) = dnekclock()


         do ie=1,nelt

            rvole = (xerange(2,1,ie) - xerange(1,1,ie))*
     >              (xerange(2,2,ie) - xerange(1,2,ie))*
     >              (xerange(2,3,ie) - xerange(1,3,ie))
            rvolei=1./rvole
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            fvalgx(i,j,k,ie) = rcountv(1,ie)*rvolei
            fvalgy(i,j,k,ie) = rcountv(2,ie)*rvolei
            fvalgz(i,j,k,ie) = rcountv(3,ie)*rvolei
            fvalgv(i,j,k,ie) = rcountv(4,ie)*rvolei
            fvalgg(i,j,k,ie) = rcountv(5,ie)*rvolei
            fvalv1(i,j,k,ie) = rcountv(6,ie)*rvolei
            fvalv2(i,j,k,ie) = rcountv(7,ie)*rvolei
            fvalv3(i,j,k,ie) = rcountv(8,ie)*rvolei
         enddo
         enddo
         enddo
         enddo
            


      ! end timer
      pttime(7) = pttime(7) + dnekclock() - ptdum(7)

      return
      end
c----------------------------------------------------------------------
      subroutine local_part_to_grid(fvalgx,fvalgy,fvalgz,fvalgv,fvalgg,
     >                              fvalv1,fvalv2,fvalv3,
     >                              pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                              ppv1,ppv2,ppv3,
     >                              xx,yy,zz,rbexpi,ralph,ralph2,e)
c
c     spread a local particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e,er
      real    fvalgx(lx1,ly1,lz1,lelt),fvalgy(lx1,ly1,lz1,lelt),
     >        fvalgz(lx1,ly1,lz1,lelt),fvalgv(lx1,ly1,lz1,lelt),
     >        fvalgg(lx1,ly1,lz1,lelt),fvalv1(lx1,ly1,lz1,lelt),
     >        fvalv2(lx1,ly1,lz1,lelt),fvalv3(lx1,ly1,lz1,lelt),
     >        pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,ppg,ppv1,ppv2,ppv3

      ! begin timer
      ptdum(7) = dnekclock()

c     this element
      call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
     >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
     >                   fvalgg(1,1,1,e),fvalv1(1,1,1,e),
     >                   fvalv2(1,1,1,e),fvalv3(1,1,1,e),
     >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)

c     faces
      do ii=1,nfacegp
         er=el_face_el_map(e,ii) + 1
         impi=el_face_proc_map(e,ii)
         if (impi .eq. nid) 
     >      call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
      enddo

c     edges
      do ii=1,nedgegp
         er=el_edge_el_map(e,ii) + 1
         impi=el_edge_proc_map(e,ii)
         if (impi .eq. nid) 
     >      call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
      enddo

c     corners
      do ii=1,ncornergp
         er=el_corner_el_map(e,ii) + 1
         impi=el_corner_proc_map(e,ii)
         if (impi .eq. nid) 
     >      call point_to_grid(fvalgx(1,1,1,er),fvalgy(1,1,1,er),
     >                   fvalgz(1,1,1,er),fvalgv(1,1,1,er),
     >                   fvalgg(1,1,1,er),fvalv1(1,1,1,er),
     >                   fvalv2(1,1,1,er),fvalv3(1,1,1,er),
     >                   xm1(1,1,1,er),ym1(1,1,1,er),zm1(1,1,1,er),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)
      enddo

      ! end timer
      pttime(7) = pttime(7) + dnekclock() - ptdum(7)

      return
      end
c----------------------------------------------------------------------
      subroutine remote_part_to_grid(fvalgx,fvalgy,fvalgz,fvalgv,fvalgg,
     >                                fvalv1,fvalv2,fvalv3,
     >                                pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                                ppv1,ppv2,ppv3,
     >                                xx,yy,zz,rbexpi,ralph,ralph2,e)
c
c     spread a remote particle property to local fluid grid points
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer e
      real    fvalgx(lx1,ly1,lz1,lelt),fvalgy(lx1,ly1,lz1,lelt),
     >        fvalgz(lx1,ly1,lz1,lelt),fvalgv(lx1,ly1,lz1,lelt),
     >        fvalgg(lx1,ly1,lz1,lelt),fvalv1(lx1,ly1,lz1,lelt),
     >        fvalv2(lx1,ly1,lz1,lelt),fvalv3(lx1,ly1,lz1,lelt),
     >        pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,ppg,
     >        ppv1,ppv2,ppv3

      ! begin timer
      ptdum(8) = dnekclock()

      call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
     >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
     >                   fvalgg(1,1,1,e),fvalv1(1,1,1,e),
     >                   fvalv2(1,1,1,e),fvalv3(1,1,1,e),
     >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
     >                   pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                   ppv1,ppv2,ppv3,
     >                   xx,yy,zz,rbexpi,
     >                   ralph,ralph2)

      ! end timer
      pttime(8) = pttime(8) + dnekclock() - ptdum(8)

      return
      end
c----------------------------------------------------------------------
      subroutine point_to_grid(gval1,gval2,gval3,gval4,gval5,
     >                      gval6,gval7,gval8,
     >                      xgd,ygd,zgd,
     >                      pvalpx,pvalpy,pvalpz,pvalpv,ppg,
     >                      ppv1,ppv2,ppv3,
     >                      xx,yy,zz,rbexpi,
     >                      ralph,ralph2)
c
c     spreads point onto grid in element e
c       gval: grid value
c       pval: particle value
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'

      integer i,j,k,e,ip
      real    gval1(lx1,ly1,lz1),gval2(lx1,ly1,lz1),gval3(lx1,ly1,lz1),
     >        gval4(lx1,ly1,lz1),gval5(lx1,ly1,lz1),gval6(lx1,ly1,lz1),
     >        gval7(lx1,ly1,lz1),gval8(lx1,ly1,lz1),
     >        xgd(lx1,ly1,lz1),ygd(lx1,ly1,lz1),zgd(lx1,ly1,lz1),
     >        pvalpx,pvalpy,pvalpz,pvalpv,pi,
     >        distx,disty,distz,xx,yy,zz,distx2,disty2,distz2,multfc,
     >        ppg,ppv1,ppv2,ppv3

      ! begin timer
      ptdum(9) = dnekclock()

c     optimized code ------------------------------------------------
c     can we skip this element?
c     rxl      = abs(xx - xgd(1,1,1))
c     rxr      = abs(xx - xgd(nx1,1,1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514
c     rxl      = abs(yy - ygd(1,1,1))
c     rxr      = abs(yy - ygd(1,ny1,1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514
c     rxl      = abs(zz - zgd(1,1,1))
c     rxr      = abs(zz - zgd(1,1,nz1))
c     if (rxl.gt.ralph .and. rxr.gt.ralph) goto 1514

      do k=1,nz1
         rtmp = 0.
         distz = zz - zgd(1,1,k) ! even element spacing only!
         distz2 = distz**2
         if (distz2 .gt. ralph2) goto 1513
      do j=1,ny1
         disty = yy - ygd(1,j,1) ! even element spacing only!
         disty2 = disty**2
c        if (disty2 .gt. ralph2) goto 1512
         rtmp1   =  distz2 + disty2
         if (rtmp1 .gt. ralph2) goto 1512
      do i=1,nx1
         distx = xx - xgd(i,1,1) ! even element spacing only!
         distx2 = distx**2
c        if (distx2 .gt. ralph2) goto 1511
         rtmp2 = rtmp1 + distx2
         if (rtmp2 .gt. ralph2) goto 1511

         rdum = rtmp2*rbexpi
         rexp = exp(rdum)
c        rexp = 1. ! change here for uniform and above in spread_prop... 

         gval1(i,j,k) = gval1(i,j,k) + pvalpx*rexp
         gval2(i,j,k) = gval2(i,j,k) + pvalpy*rexp
         gval3(i,j,k) = gval3(i,j,k) + pvalpz*rexp
         gval4(i,j,k) = gval4(i,j,k) + pvalpv*rexp
         gval5(i,j,k) = gval5(i,j,k) + ppg   *rexp
         gval6(i,j,k) = gval6(i,j,k) + ppv1  *rexp
         gval7(i,j,k) = gval7(i,j,k) + ppv2  *rexp
         gval8(i,j,k) = gval8(i,j,k) + ppv3  *rexp
 1511 continue
      enddo
 1512 continue
      enddo
 1513 continue
      enddo
 1514 continue

      ! end timer
      pttime(9) = pttime(9) + dnekclock() - ptdum(9)

      return
      end
c-----------------------------------------------------------------------
      subroutine search_nearest_neighbor
      include 'SIZE'
      include 'TOTAL'
      include 'CMTDATA'
      include 'CMTPART'
c
c     this routine implements a naive particle search. This should be
c     updated for the future with some kind of more recent algorithm for
c     nearest neighbor searching. Particles will need to search local
c     particles (in rpart and ipart) and remote particles that are
c     nearby but on different MPI ranks (rptsgp and iptsgp of length nfptsgp)
c
c     particles will check if they are within d2chk of each other
c

      integer nneigh

      ! begin timer
      ptdum(15) = dnekclock()

      d3 = 0.5*d2chk(1) ! user can change, but d2chk is robust max value
                        ! Note: 1/2*d2chk seems to work even w/outflow

c     let every particle search for itself
      do i = 1,n
         ipart(jai,i) = ipart(jpnn,i) ! for testing
         nneigh = 0
c        particles in local elements
         do j = 1,n
            if (i .ne. j) then
               pdist = abs(rpart(jx,i)-rpart(jx,j))**2  
     >                          + abs(rpart(jy,i)-rpart(jy,j))**2
     >                          + abs(rpart(jz,i)-rpart(jz,j))**2
               pdist = sqrt(pdist)
               if (pdist .gt. d3) goto 1109
               nneigh = nneigh + 1
            endif
1109        continue
         enddo

c        search list of ghost particles
         do j = 1,nfptsgp
            if (iptsgp(jgpes,j).eq. ipart(je0,i)) then ! exclude ghosts not
                                                      ! meant for this eleme
            pdist = abs(rpart(jx,i)-rptsgp(jgpx,j))**2  
     >                    + abs(rpart(jy,i)-rptsgp(jgpy,j))**2
     >                    + abs(rpart(jz,i)-rptsgp(jgpz,j))**2
            pdist = sqrt(pdist)
            if (pdist .gt. d3) goto 11092
            nneigh = nneigh + 1
            endif
11092       continue
         enddo
         ipart(jpnn,i) = nneigh
         ipart(jai,i) = ipart(jai,i) - ipart(jpnn,i) ! comptued distance
                                                     ! for testing
      enddo

      ! end timer
      pttime(15) = pttime(15) + dnekclock() - ptdum(15)

      return
      end
c----------------------------------------------------------------------
      subroutine point_to_grid_corr_init
c
c
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'CMTDATA'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      integer i,j,k,e,ip
      real    xx,yy,zz,msum

      common /point2gridc/ p2gc
      real   p2gc(lx1,ly1,lz1,lelt,4)

      if (nid.eq.0) write(6,*) 'Starting point_to_grid_corr_init'

c     local mpi rank effects
      do ie=1,nelt
         xs = xerange(1,1,ie)
         xe = xerange(2,1,ie)
         ys = xerange(1,2,ie)
         ye = xerange(2,2,ie)
         zs = xerange(1,3,ie)
         ze = xerange(2,3,ie)
         xdelta = (xe-xs)/(nx1-1)
         ydelta = (ye-ys)/(ny1-1)
         zdelta = (ze-zs)/(nz1-1)
         do k=1,nz1
            zz = zs + (k-1)*zdelta
         do j=1,ny1
            yy = ys + (j-1)*ydelta
         do i=1,nx1
            xx = xs + (i-1)*xdelta
c           p2gc(i,j,k,ie,1) = xx
c           p2gc(i,j,k,ie,2) = yy 
c           p2gc(i,j,k,ie,3) = zz 

            p2gc(i,j,k,ie,1) = xm1(i,j,k,ie)
            p2gc(i,j,k,ie,2) = ym1(i,j,k,ie) 
            p2gc(i,j,k,ie,3) = zm1(i,j,k,ie) 
         enddo
         enddo
         enddo
      enddo

      do ie=1,nelt
      do k=1,nz1
      do j=1,ny1
      do i=1,nx1
         xx = p2gc(i,j,k,ie,1) 
         yy = p2gc(i,j,k,ie,2) 
         zz = p2gc(i,j,k,ie,3) 
         call compute_gamma_grid(ie,xx,yy,zz,p2gc(i,j,k,ie,4))
      enddo
      enddo
      enddo
      enddo

      if (nid.eq.0) write(6,*) 'Ending point_to_grid_corr_init'

      return
      end
c----------------------------------------------------------------------
      subroutine compute_gamma_grid(ie,xx,yy,zz,gam_val)
c
c
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'MASS'
      include 'CMTDATA'
      include 'CMTPART'

      integer e,er
      real    msum,msum_total,pi,multfc
      real    xx,yy,zz,Lx,Ly,Lz,rdum(lx1,ly1,lz1)
      real    mesharound(81),gam_val,dumval(lx1,ly1,lz1,27)
      real    xgd(lx1,ly1,lz1),ygd(lx1,ly1,lz1),zgd(lx1,ly1,lz1)

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      xs = xerange(1,1,ie)
      xe = xerange(2,1,ie)
      ys = xerange(1,2,ie)
      ye = xerange(2,2,ie)
      zs = xerange(1,3,ie)
      ze = xerange(2,3,ie)

      Lx = xe - xs
      Ly = ye - ys
      Lz = ze - zs

      mesharound = (/ 
     >                0. ,0. ,0. , ! 1
     >                -Lx,-Ly,0. , ! 2
     >                0. ,-Ly,0. , ! 3
     >                Lx ,-Ly,0. , ! 4
     >                -Lx,0. ,0. , ! 5
     >                Lx ,0. ,0. , ! 6
     >                -Lx,Ly ,0. , ! 7
     >                0. ,Ly ,0. , ! 8
     >                Lx ,Ly ,0. , ! 9
     >                0. ,0. ,-Lz, ! 10
     >                -Lx,-Ly,-Lz, ! 11
     >                0. ,-Ly,-Lz, ! 12
     >                Lx ,-Ly,-Lz, ! 13
     >                -Lx,0. ,-Lz, ! 14
     >                Lx ,0. ,-Lz, ! 15
     >                -Lx,Ly ,-Lz, ! 16
     >                0. ,Ly ,-Lz, ! 17
     >                Lx ,Ly ,-Lz, ! 18
     >                0. ,0. ,Lz , ! 19
     >                -Lx,-Ly,Lz , ! 20
     >                0. ,-Ly,Lz , ! 21
     >                Lx ,-Ly,Lz , ! 22
     >                -Lx,0. ,Lz , ! 23
     >                Lx ,0. ,Lz , ! 24
     >                -Lx,Ly ,Lz , ! 25
     >                0. ,Ly ,Lz , ! 26
     >                Lx ,Ly ,Lz   ! 27
     >                            /)

      call rzero(dumval,lx1*ly1*lz1*27)
      call rzero(rdum,lx1*ly1*lz1)

      pi       = 4.0d+0*atan(1.0d+0)
      multfc   = 1./(sqrt(2.*pi)**3 * rsig**3)
      rbexpi   = 1./(-2.*rsig**2)

c     ralphd   = 1E10   ! dummy so it will spread everywhere
c     ralphd2   = 1E10  ! dummy so it will spread everywhere

      ralphd    = d2chk(1)     ! assume all directions same!
      ralphd2   = ralphd**2

      do iie=1,27         
         ioff = (iie-1)*3
         do k=1,ny1
         do j=1,ny1
         do i=1,nx1
            xgd(i,j,k) = xm1(i,j,k,ie) + mesharound(ioff+1)
            ygd(i,j,k) = ym1(i,j,k,ie) + mesharound(ioff+2)
            zgd(i,j,k) = zm1(i,j,k,ie) + mesharound(ioff+3)
         enddo
         enddo
         enddo
         
         rpass = 1.*multfc

c        call point_to_grid(dumval(1,1,1,iie),1.,xx,yy,zz,1.,
c    >             xgd,ygd,zgd)

         call point_to_grid(rdum(1,1,1),rdum(1,1,1),
     >                  rdum(1,1,1),dumval(1,1,1,iie),
     >                  rdum(1,1,1),rdum(1,1,1),rdum(1,1,1),rdum(1,1,1),
     >                  xgd,ygd,zgd,
     >                  1.,1.,1.,rpass,1.,1.,1.,1.,xx,yy,zz,rbexpi,
     >                  ralphd,ralphd2)

c        call point_to_grid(fvalgx(1,1,1,e),fvalgy(1,1,1,e),
c    >                   fvalgz(1,1,1,e),fvalgv(1,1,1,e),
c    >                   xm1(1,1,1,e),ym1(1,1,1,e),zm1(1,1,1,e),
c    >                   pvalpx,pvalpy,pvalpz,pvalpv,xx,yy,zz,rbexpi,
c    >                   ralph,ralph2)
c        do k=1,ny1
c        do j=1,ny1
c        do i=1,nx1
c           print *, i,j,k,xgd(i,j,k),ygd(i,j,k),zgd(i,j,k),
c    >               dumval(i,j,k,iie)
c        enddo
c        enddo
c        enddo
      enddo

      msum = 0.
      do iie=1,27
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            msum = msum + dumval(i,j,k,iie)*bm1(i,j,k,ie)
         enddo
         enddo
         enddo
      enddo
c     msum_total = glsum(msum,1)
      if (abs(msum) .lt. 1E-16) then
         gam_val = 1.
      else
         gam_val = 1./msum
      endif

      return
      end
c----------------------------------------------------------------------
      subroutine correct_spl
c
c     correct initial super particle loading
c
      include 'SIZE'
      include 'INPUT'
      include 'GEOM'
      include 'SOLN'
      include 'CMTDATA'
      include 'MASS'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real   xdrange(2,3)
      common /domainrange/ xdrange

      real rdumvol(llpart,2*3)

      do i=1,n
         rdumvol(i,1) = rpart(jvol,i)  ! particle volume
         rdumvol(i,2) = rpart(jvol1,i) ! interp vol frac @ part loc
         rdumvol(i,3) = rpart(jspl,i)  ! super part. loading
      enddo

c     call usr_particles_io(istep)

c     begin diagnostics ----
c
c     eulerian volume frac 
      nxyze = nx1*ny1*nz1*nelt
      rmu1  = glsc2(bm1,ptw(1,1,1,1,4),nxyze)
      rmu1  = rmu1/vol_distrib

c
c     lagrangian volume frac
      rmu2  = glsum(rdumvol(1,2),n)
      rmu2  = rmu2/nw
      rmin2 = glmin(rdumvol(1,2),n)
      rmax2 = glmax(rdumvol(1,2),n)

c
c     what spl mean should be
      rdumt   = glsum(rdumvol(1,1),n)
      rsplavg = phi_desire*vol_distrib/rdumt

c
c     what spl mean actually is
      rmu3  = glsum(rdumvol(1,3),n)
      rmu3  = rmu3/nw
      rmin3 = glmin(rdumvol(1,3),n)
      rmax3 = glmax(rdumvol(1,3),n)

c
c     variance and skew stuff
      do i=1,n
         rdumvol(i,4) = (rpart(jvol1,i) - rmu2)**2
         rdumvol(i,5) = (rpart(jspl,i) - rmu3)**2
      enddo

      rvar2 = glsum(rdumvol(1,4),n)
      rvar2 = rvar2/nw

      rvar3 = glsum(rdumvol(1,5),n)
      rvar3 = rvar3/nw

      if (nid.eq.0) write(6,*) '-DZ- Md,bd,Me --'
      if (nid.eq.0) write(6,*) phi_desire,rsplavg,rmu1
      if (nid.eq.0) write(6,*) '-DZ-- Ml,Sl,Minl,Maxl'
      if (nid.eq.0) write(6,*) rmu2,sqrt(rvar2),rmin2,rmax2
      if (nid.eq.0) write(6,*) '-DZ--- Mb,Sb,Minb,Maxb'
      if (nid.eq.0) write(6,*) rmu3,sqrt(rvar3),rmin3,rmax3
c     end diagnostics ----


      do ip=1,n
c        linear roll off near bed edges
c        rthresh = 2.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (0.-phi_desire)/(rxmin-rxmax)
c           phi_val = 0. + rm*(rpart(jx,ip)-rxmin)
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (phi_desire-0.)/(rxmin-rxmax)
c           phi_val = phi_desire + rm*(rpart(jx,ip)-rxmin)
c        else
c           phi_val = phi_desire
c        endif

c        tanh roll off near bed edges ! do this one if bed edge!!
c        rthresh = 4.*rleng
c        rxboavg = (rxbo(2,1)+rxbo(1,1))/2.
c        if (rpart(jx,ip) .lt. rxboavg) then
c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=0.
c           rbp=phi_desire
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rap + rd*(0.5 + 0.5*tanh(200.*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxboavg) then
c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=phi_desire
c           rbp=0.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rap + rd*(0.5 + 0.5*tanh(200.*(rpart(jx,ip)-rc)))
c        else
c           phi_val = phi_desire
c        endif

         phi_val = phi_desire ! comment out if bed and uncomment above
         rtmp = 0.30*rpart(jspl,ip)
         rxi = rtmp*(1. - rpart(jvol1,ip)/phi_val)
         rpart(jspl,ip)=rpart(jspl,ip) + rxi
         if (rpart(jspl,ip).lt.0.) rpart(jspl,ip) = 0.

c        rthresh = 2.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (1.-rsplavg)/(rxmin-rxmax)
c           rpart(jspl,ip) = 1. + rm*(rpart(jx,ip)-rxmin)

c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=1.
c           rbp=rsplavg
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           rpart(jspl,ip) = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (rsplavg-1.)/(rxmin-rxmax)
c           rpart(jspl,ip) = rsplavg + rm*(rpart(jx,ip)-rxmin)

c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=rsplavg
c           rbp=1.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           rpart(jspl,ip) = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        else
c           rtmp = 0.30*rpart(jspl,ip)
c           rxi = rtmp*(1. - rpart(jvol1,ip)/phi_desire)
c           rpart(jspl,ip)=rpart(jspl,ip) + rxi
c        endif
c        if (rpart(jspl,ip).lt.0.) rpart(jspl,ip) = 0.


c        linear roll off near bed edges
c        rthresh = 1.*rleng
c        if (rpart(jx,ip) .lt. rxbo(1,1)+rthresh) then
c           rxmin = rxbo(1,1)
c           rxmax = rxbo(1,1)+rthresh
c           rm    = (0.-phi_desire)/(rxmin-rxmax)
c           phi_val = 0. + rm*(rpart(jx,ip)-rxmin)
c        elseif (rpart(jx,ip) .gt. rxbo(2,1)-rthresh) then
c           rxmin = rxbo(2,1)-rthresh
c           rxmax = rxbo(2,1)
c           rm    = (phi_desire-0.)/(rxmin-rxmax)
c           phi_val = phi_desire + rm*(rpart(jx,ip)-rxmin)
c        else
c           phi_val = phi_desire
c        endif

c        tanh roll off near bed edges
c        rthresh = 2.*rleng
c        rxboavg = (rxbo(2,1)+rxbo(1,1))/2.
c        if (rpart(jx,ip) .lt. rxboavg) then
c           ra = rxbo(1,1)
c           rb = rxbo(1,1)+rthresh
c           rap=0.
c           rbp=phi_desire
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        elseif (rpart(jx,ip) .gt. rxboavg) then
c           ra = rxbo(2,1)-rthresh
c           rb = rxbo(2,1)
c           rap=phi_desire
c           rbp=0.
c           rc = (ra+rb)/2.
c           rd = rbp - rap
c           phi_val = rbp + rd*(0.5 + 
c    >                 0.5*tanh(100.*(rb-ra)/rc*(rpart(jx,ip)-rc)))
c        else
c           phi_val = phi_desire
c        endif



 1511 continue
      enddo


      return
      end
c----------------------------------------------------------------------
