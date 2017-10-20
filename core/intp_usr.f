c
c interpolation wrapper for usage in .usr file
c note: don't call outside
c
c-----------------------------------------------------------------------
      subroutine intp_setup(tolin)
c
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal
      common /intp/   tol
      common /intp_h/ ih_intp


      tol = tolin
      if (tolin.lt.0) tol = 5e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.01 ! relative size to expand bounding boxes by
c
      if(nidd.eq.0) write(6,*) 'call intp_setup(), tol=', tol
      call fgslib_findpts_setup(ih_intp,nekcomm,npp,ndim,
     &                          xm1,ym1,zm1,nx1,ny1,nz1,
     &                          nelt,nxf,nyf,nzf,bb_t,n,n,
     &                          npt_max,tol)
c       
      return
      end
c-----------------------------------------------------------------------
      subroutine intp_do(fldout,fldin,nfld,xp,yp,zp,n,iwk,rwk,nmax,iflp)
c
c fldout    ... output field(s) dim (*,nfld)
c fldin     ... source field(s) dim (*,nfld)
c nfld      ... number of fields
c xp,yp,zp  ... interpolation points dim (n)
c n         ... number of points
c iwk       ... integer working array dim > (nmax,3)
c rwk       ... real working array dim > (nmax,ldim+1)
c nmax      ... maximum number of points per MPI rank
c iflp      ... look-up interpolation points
c
      include 'SIZE'

      common /intp/   tol
      common /intp_h/ ih_intp

      real    fldin(*),fldout(*)
      real    xp(*),yp(*),zp(*)

      real    rwk(nmax,*)
      integer iwk(nmax,*) 

      logical iflp

      integer nn(2)
      logical ifot


      ifot = .false. ! transpose output field

      if(nio.eq.0) write(6,*) 'call intp_do'

      if(n.gt.nmax) then
        write(6,*)
     &   'ABORT: n>nmax in intp_do', n, nmax
        call exitt
      endif

      ! locate points (iel,iproc,r,s,t)
      nfail = 0
      if(iflp) then
        if(nio.eq.0) write(6,*) 'call findpts'
        call fgslib_findpts(ih_intp,
     &                      iwk(1,1),1,
     &                      iwk(1,3),1,
     &                      iwk(1,2),1,
     &                      rwk(1,2),ndim,
     &                      rwk(1,1),1,
     &                      xp,1,
     &                      yp,1,
     &                      zp,1,n)
        do in=1,n
           ! check return code
           if(iwk(in,1).eq.1) then
             if(rwk(in,1).gt.10*tol) then
               nfail = nfail + 1
               if (nfail.le.5) write(6,'(a,1p4e15.7)')
     &     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',
     &     xp(in),yp(in),zp(in),rwk(in,1)
             endif
           elseif(iwk(in,1).eq.2) then
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
         call fgslib_findpts_eval(ih_intp,fldout(iout),is_out,
     &                            iwk(1,1),1,
     &                            iwk(1,3),1,
     &                            iwk(1,2),1,
     &                            rwk(1,2),ndim,n,
     &                            fldin(iin))
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

      common /intp_h/ ih_intp


      call fgslib_findpts_free(ih_intp)

      return
      end

