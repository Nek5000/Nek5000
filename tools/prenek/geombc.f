c-----------------------------------------------------------------------
      subroutine set_geom_bc(ifld)
c
c     Set BCs based upon geometric information
c
      include 'basics.inc'
c
c
      call compute_adj_srf




uuu
c-----------------------------------------------------------------------
      subroutine compute_adj_srf(ia,ja,vals,x_centroid,y_centroid
     $     ,z_centroid,n,vertex,nv)
c
      INCLUDE 'basics.inc'
      INCLUDE 'basicsp.inc'
c
      integer n,nv
      integer ia(1),ja(1),vertex(1)
      real*8  vals(1),x_centroid(1),y_centroid(1),z_centroid(1)
c
      integer iseg,i_min,nc,k,ie,ic,kl  
      integer start_pt(1),start_el(1),ind(1),ninseg(1)
      real    dx(1),dy(1),dz(1)
      logical ifseg(1)
c
      equivalence (start_pt,rxm1)
      equivalence (start_el,rym1)
      equivalence (ind     ,rzm1)
      equivalence (dx      ,sxm1)
      equivalence (dy      ,sym1)
      equivalence (dz      ,szm1)
      equivalence (ninseg  ,txm1)
      equivalence (ifseg   ,tym1)
      equivalence (ww      ,tzm1)
c
      integer      va(1)
      equivalence (va      ,sxm1)
C
C
C     only works if each element is defined by exactly 8 pts!!!
C
      n = nel
      nc = 2**ndim
      nv = n*nc
c
c
c     Lexicographical sort to establish adj. array:
c     Any elements sharing associated vertices share and edge.
c
      k  = 0
      do ie=1,nel
         xmin = x(1,ie)
         xmax = x(1,ie)
         ymin = y(1,ie)
         ymax = y(1,ie)
         zmin = z(1,ie)
         zmax = z(1,ie)
         x_centroid(ie) = 0.0
         y_centroid(ie) = 0.0
         z_centroid(ie) = 0.0
         kl   = k+1
         do ic=1,nc
            k=k+1
            xp(k)=x(ic,ie)
            xmin = min(x(ic,ie),xmin)
            xmax = max(x(ic,ie),xmax)
            x_centroid(ie)=x_centroid(ie)+x(ic,ie)
            yp(k)=y(ic,ie)
            ymin = min(y(ic,ie),ymin)
            ymax = max(y(ic,ie),ymax)
            y_centroid(ie)=y_centroid(ie)+y(ic,ie)
            if (if3d) then
               zp(k)=z(ic,ie)
               zmin = min(z(ic,ie),zmin)
               zmax = max(z(ic,ie),zmax)
               z_centroid(ie)=z_centroid(ie)+z(ic,ie)
            else
               zmin          =0.
               zmax          =0.
               z_centroid(ie)=0.
            endif
            vertex(k) = 9999999
            start_pt(k) = k
            start_el(k) = ie
         enddo
         x_centroid(ie)=x_centroid(ie)/nc
         y_centroid(ie)=y_centroid(ie)/nc
         z_centroid(ie)=z_centroid(ie)/nc
c
         endif
         dxe = .001*(xmax-xmin)
         dye = .001*(ymax-ymin)
         dze = .001*(zmax-zmin)
         call cfill(dx(kl),dxe,nc)
         call cfill(dy(kl),dye,nc)
         call cfill(dz(kl),dze,nc)
      enddo
c
      if (k.ne.nv) then
         write(6,*) '# elements          = ',nel
         write(6,*) '# vertices          = ',nv
         write(6,*) 'k                   = ',k
      else
         write(6,*) '# elements          = ',nel
         write(6,*) '# vertices          = ',nv
      endif
c
c     Lexicographical sorting scheme... modifies x,y,z,dx,dy,dz
c
      call lfalse(ifseg   ,nv)
      nseg        = 1
      ifseg(1)    = .true.
      ninseg(1)   = nv
c
      do j=1,ndim
c
c       Sort within each segment
c
         ioff=1
         do iseg=1,nseg
            if  (j.eq.1) then
               call rank (xp(ioff),ind,ninseg(iseg))
            elseif  (j.eq.2) then
               call rank (yp(ioff),ind,ninseg(iseg))
            else
               call rank (zp(ioff),ind,ninseg(iseg))
            endif
            call iswap (start_pt(ioff),ww,ind,ninseg(iseg))
            call iswap (start_el(ioff),ww,ind,ninseg(iseg))
            call  swap (xp      (ioff),ww,ind,ninseg(iseg))
            call  swap (yp      (ioff),ww,ind,ninseg(iseg))
            call  swap (zp      (ioff),ww,ind,ninseg(iseg))
            call  swap (dx      (ioff),ww,ind,ninseg(iseg))
            call  swap (dy      (ioff),ww,ind,ninseg(iseg))
            call  swap (dz      (ioff),ww,ind,ninseg(iseg))
            ioff=ioff+ninseg(iseg)
         enddo
c
c       Check for jumps in current coordinate
c
         if (j.eq.1) then
            do i=2,nv
               if (abs(xp(i)-xp(i-1)).gt.min(dx(i),dx(i-1)))
     $              ifseg(i)=.true.
            enddo
         elseif (j.eq.2) then
            do i=2,nv
               if (abs(yp(i)-yp(i-1)).gt.min(dy(i),dy(i-1)))
     $              ifseg(i)=.true.
            enddo
         elseif (j.eq.3) then
            do i=2,nv
               if (abs(zp(i)-zp(i-1)).gt.min(dz(i),dz(i-1)))
     $              ifseg(i)=.true.
            enddo
         endif
c
c       Count up number of different segments
c
         nseg = 0
         do i=1,nv
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
c
c     Data now sorted lexicographically!
c
c     To sort additional incoming data:
c
c     >
c     >      do k=1,n
c     >         t_sort(k) = t_in(init_loc(k))
c     >      enddo
c     >
c
      return
      end
c-----------------------------------------------------------------------
