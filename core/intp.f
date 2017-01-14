      subroutine intp_setup(tolin)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nidd.eq.0) write(6,*) 'call intp_setup(), tol=', tol
      call findpts_setup(ih_intp,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end
c-----------------------------------------------------------------------
      subroutine intp_do(fieldout,fieldin,nfld,xp,yp,zp,n,ifpts)
c
c wrapper to interpolate input field at given points
c
c in:
c fieldin   ... source field(s)
c nfld      ... number of fields in fieldin
c xp,yp,zp  ... interpolation xp(n),yp(n),zp(n)
c n         ... dim of xp,yp,zp
c ifpts     ... find interpolation points
c ih        ... interpolation handle
c
c out:
c fieldout  ... target (interpolated) field(s) (:,1:nfld)
c
      include 'SIZE'

      real    fieldin(*),fieldout(*)
      real    xp(*),yp(*),zp(*)

      real    dist(lpts) ! squared distance
      real    rst(lpts*ldim)
      common /intp_r/ rst,dist

      integer rcode(lpts),elid(lpts),proc(lpts)
      common /intp_i/ rcode,elid,proc

      integer nn(2)
      logical ifot,ifpts

      ifot = .false.

      if(nio.eq.0) write(6,*) 'call intp_do'

      if(n.gt.lpts) then
        write(6,*)
     &   'ABORT: intpts() n>lpts, increase lelt in SIZE', n, lpts
        call exitt
      endif

      ! locate points (iel,iproc,r,s,t)
      nfail = 0
      if(ifpts) then
        if(nio.eq.0) write(6,*) 'call findpts'
        call findpts(ih_intp,rcode,1,
     &               proc,1,
     &               elid,1,
     &               rst,ndim,
     &               dist,1,
     &               xp,1,
     &               yp,1,
     &               zp,1,n)
        do in=1,n
           ! check return code
           if(rcode(in).eq.1) then
             if(dist(in).gt.1e-12) then
               nfail = nfail + 1
               if (nfail.le.5) write(6,'(a,1p4e15.7)')
     &     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',
     &     xp(in),yp(in),zp(in),dist(in)
             endif
           elseif(rcode(in).eq.2) then
             nfail = nfail + 1
             if (nfail.le.5) write(6,'(a,1p3e15.7)')
     &        ' WARNING: point not within mesh xy[z]: !',
     &        xp(in),yp(in),zp(in)
           endif
        enddo
      endif

      ! evaluate inut field at given points
      ltot = lelt*lx1*ly1*lz1
      do ifld = 1,nfld
         iin    = (ifld-1)*ltot + 1
         iout   = (ifld-1)*n + 1
         is_out = 1
         if(ifot) then ! transpose output
           iout   = ifld
           is_out = nfld
         endif
         call findpts_eval(ih_intp,fieldout(iout),is_out,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,n,
     &                     fieldin(iin))
      enddo

      nn(1) = iglsum(n,1)
      nn(2) = iglsum(nfail,1)
      if(nio.eq.0) then
        write(6,1) nn(1),nn(2)
  1     format('   total number of points = ',i12,/,'   failed = '
     &         ,i12,/,' done :: intp_do')
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_free()

      call findpts_free(ih_intp)

      return
      end

