c-----------------------------------------------------------------------
      subroutine init_interpolation
      include 'SIZE' 
      include 'INPUT' 
      include 'CMTPART' 
c
c     calculates the barycentric lagrange weights
c

c     get gll points in all directions
      call zwgll(xgll,wxgll,lx1)
      call zwgll(ygll,wygll,ly1)
      call rone(zgll,lz1)
      if(if3d) call zwgll(zgll,wzgll,lz1)
c     set all weights to ones first
      call rone(wxgll,lx1)
      call rone(wygll,ly1)
      call rone(wzgll,lz1)
c
c     copy for reduced interpolation
      nx1r = lx1
      if (red_interp.gt.0) then
         nx1r = red_interp
         ic = 0
         do j=1,lx1,2
            ic = ic + 1
            xgll(ic) = xgll(j)
            ygll(ic) = ygll(j)
            zgll(ic) = zgll(j)
         enddo
      endif

c     calc x bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wxgll(j) = wxgll(j)/(xgll(j) - xgll(k))
            endif
         enddo
      enddo
c     calc y bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wygll(j) = wygll(j)/(ygll(j) - ygll(k))
            endif
         enddo
      enddo
c     calc z bary weights
      do j=1,nx1r
         do k=1,nx1r
            if (j .NE. k) then
               wzgll(j) = wzgll(j)/(zgll(j) - zgll(k))
            endif
         enddo
      enddo

      return 
      end
c-----------------------------------------------------------------------
      subroutine init_baryinterp(x,y,z,nxyz)
c     used for 3d interpolation only
      include 'SIZE'
      include 'CMTPART'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot

      real x, y, z, repy, repz,repx,diff
      real bwgtx(lx1),bwgty(ly1),bwgtz(lz1)

      bot= 0.00
      do k=1,nx1r
         diff = z - zgll(k)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtz(k) = wzgll(k)/diff
      enddo
      do i=1,nx1r
         diff = x - xgll(i)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgtx(i) = wxgll(i)/diff
      enddo 
      do j=1,nx1r
         diff = y-ygll(j)
           if (abs(diff) .le. 1E-16) diff = sign(1E-16,diff)
         bwgty(j) = wygll(j)/diff
      enddo

      do k=1,nx1r
      do j=1,nx1r
         repdum = bwgty(j)*bwgtz(k)
      do i=1,nx1r
         rep(i,j,k) =  repdum* bwgtx(i)
         bot        =  bot + rep(i,j,k)
      enddo
      enddo
      enddo 

      do k=1,nx1r
      do j=1,nx1r
      do i=1,nx1r
         rep(i,j,k) =  rep(i,j,k)/bot
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine baryinterp(field,pofx,nxyz)
c     used for 3d interpolation only
      include 'SIZE'
      include 'CMTPART'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot
      real field(1),pofx,top

      ! begin timer
      ptdum(18) = dnekclock()

      pofx = 0.00
      if (nx1r.eq.lx1) then ! full interpolation
      do i=1,nxyz
         pofx =  pofx + rep(i,1,1)*field(i)
      enddo
      else                  ! reduced interpolation

      kk = 0 
      do k=1,nx1,2
         kk = kk + 1
         jj = 0
         ijk3 = (k-1)*nx1**2
      do j=1,nx1,2
         jj = jj + 1
         ii = 0
         ijk2 = ijk3+(j-1)*nx1
      do i=1,nx1,2
         ii   = ii + 1
         ijk1 = ijk2 + i
         pofx =  pofx + rep(ii,jj,kk)*field(ijk1)
      enddo
      enddo
      enddo
      endif

      ! end timer
      pttime(18) = pttime(18) + dnekclock() - ptdum(18)

      return
      end
c-----------------------------------------------------------------------
      subroutine triinterp(xf,yf,zf,field,x,y,z,r,s,t,ie,pval)
c     
c     used for 3d trilinear interpolation
c
      include 'SIZE'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      real field(nx1,ny1,nz1),xf(nx1,ny1,nz1),yf(nx1,ny1,nz1),
     >                        zf(nx1,ny1,nz1)
      real x,y,z,pval,c00,c01,c10,c11,c0,c1_0,c1_1,r,s,t

      ! begin timer
      ptdum(19) = dnekclock()

      rdelta = 2./(nx1-1.)
      sdelta = 2./(ny1-1.)
      tdelta = 2./(nz1-1.)

      mxx = floor((1.+r)/rdelta)+1
      myy = floor((1.+s)/sdelta)+1
      mzz = floor((1.+t)/tdelta)+1

      xd = (x - xf(mxx,myy,mzz))/(xf(mxx+1,myy,mzz)-xf(mxx,myy,mzz))
      yd = (y - yf(mxx,myy,mzz))/(yf(mxx,myy+1,mzz)-yf(mxx,myy,mzz))
      zd = (z - zf(mxx,myy,mzz))/(zf(mxx,myy,mzz+1)-zf(mxx,myy,mzz))

      c00=field(mxx,myy,mzz)*(1.-xd)+field(mxx+1,myy,mzz)*xd
      c01=field(mxx,myy,mzz+1)*(1.-xd)+field(mxx+1,myy,mzz+1)*xd
      c10=field(mxx,myy+1,mzz)*(1.-xd)+field(mxx+1,myy+1,mzz)*xd
      c11=field(mxx,myy+1,mzz+1)*(1.-xd)+field(mxx+1,myy+1,mzz+1)*xd

      c1_0 = c00*(1.-yd) + c10*yd
      c1_1 = c01*(1.-yd) + c11*yd

      pval = c1_0*(1.-zd) + c1_1*zd

      ! end timer
      pttime(19) = pttime(19) + dnekclock() - ptdum(19)

      return
      end
c-----------------------------------------------------------------------
      subroutine interp_props_part_location
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'CMTDATA'
      include 'CMTPART'
      include 'GEOM'
      common /BARRYREP/ rep, bot
      real              rep(lx1,ly1,lz1), bot

      common /point2gridc/ p2gc
      real   p2gc(lx1,ly1,lz1,lelt,4)

      ! begin timer
      ptdum(20) = dnekclock()

      nxyz = nx1*ny1*nz1
        do i=1,n
           rrdum = 1.0
           if(if3d) rrdum = rpart(jr+2,i)
c          init. barycentric interp for this particle (cost savings)
           call init_baryinterp(rpart(jr,i),rpart(jr+1,i),rrdum,nxyz)

c          interpolate fields for this particle
           ie  =  ipart(je0,i) + 1
           call baryinterp(vx(1,1,1,ie),rpart(ju0,i),nxyz)   !fluid uvel
           call baryinterp(vy(1,1,1,ie),rpart(ju0+1,i),nxyz) !fluid vvel
           if (if3d) call baryinterp(vz(1,1,1,ie),           !fluid wvel
     >                                  rpart(ju0+2,i),nxyz)

           call baryinterp(t(1,1,1,ie,1),rpart(jtempf,i),nxyz)!fluid temp
           call baryinterp(vtrans(1,1,1,ie,1),rpart(jrho,i), !fluid dens
     >                                                  nxyz) 

           call baryinterp(ptw(1,1,1,ie,4),rpart(jvol1,i),   !change!!
     >                                                  nxyz) 
c             call triinterp(xm1(1,1,1,ie),ym1(1,1,1,ie),
c    >                       zm1(1,1,1,ie),ptw(1,1,1,ie,4),
c    >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
c    >                       rpart(jr,i),rpart(jr+1,i),rpart(jr+2,i),
c    >                       ie,rpart(jvol1,i))

           call baryinterp(rhs_fluidp(1,1,1,ie,1),           !dp/dx
     >                                rpart(jDuDt,i),nxyz)
           call baryinterp(rhs_fluidp(1,1,1,ie,2),           !dp/dy
     >                                rpart(jDuDt+1,i),nxyz)
           call baryinterp(rhs_fluidp(1,1,1,ie,3),           !dp/dz
     >                                rpart(jDuDt+2,i),nxyz)

           if (two_way.ne.0) ! corrections for integration
     >     call baryinterp(p2gc(1,1,1,ie,4),rpart(jgam,i),   !change!!
     >                                                  nxyz) 
c    >        call triinterp(p2gc(1,1,1,ie,1),p2gc(1,1,1,ie,2),
c    >                       p2gc(1,1,1,ie,3),p2gc(1,1,1,ie,4),
c    >                       rpart(jx,i),rpart(jy,i),rpart(jz,i),
c    >                       rpart(jr,i),rpart(jr+1,i),rpart(jr+2,i),
c    >                       ie,rpart(jgam,i))
        enddo

      ! end timer
      pttime(20) = pttime(20) + dnekclock() - ptdum(20)

      return
      end
c----------------------------------------------------------------------
