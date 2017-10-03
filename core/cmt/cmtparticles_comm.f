c----------------------------------------------------------------------
      subroutine move_particles_inproc
c     Interpolate fluid velocity at current xyz points and move
c     data to the processor that owns the points.
c     Input:    n = number of points on this processor
c     Output:   n = number of points on this processor after the move
c     Code checks for n > llpart and will not move data if there
c     is insufficient room.
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
      include 'CMTPART'

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /myparth/ i_fp_hndl, i_cr_hndl

      integer icalld1
      save    icalld1
      data    icalld1 /0/

      logical partl         ! This is a dummy placeholder, used in cr()

      ! begin timer
      ptdum(23) = dnekclock()

      nl = 0                ! No logicals exchanged

      if (icalld1.eq.0) then
         tolin = 1.e-12
         if (wdsize.eq.4) tolin = 1.e-6
         call intpts_setup  (tolin,i_fp_hndl)
         call crystal_setup (i_cr_hndl,nekcomm,np)
         icalld1 = icalld1 + 1
      endif

      call particles_in_nid

      call findpts(i_fp_hndl !  stride     !   call findpts( ihndl,
     $           , ifpts(jrc,1),lif        !   $             rcode,1,
     $           , ifpts(jpt,1),lif        !   &             proc,1,
     $           , ifpts(je0,1),lif        !   &             elid,1,
     $           , rfpts(jr ,1),lrf        !   &             rst,ndim,
     $           , rfpts(jd ,1),lrf        !   &             dist,1,
     $           , rfpts(jx ,1),lrf        !   &             pts(    1),1,
     $           , rfpts(jy ,1),lrf        !   &             pts(  n+1),1,
     $           , rfpts(jz ,1),lrf ,nfpts)    !   &             pts(2*n+1),1,n)

      nmax = iglmax(n,1)
      if (nmax.gt.llpart) then
         if (nid.eq.0) write(6,1) nmax,llpart
    1    format('WARNING: Max number of particles:',
     $   i9,'.  Not moving because llpart =',i9,'.')
      else
c        copy rfpts and ifpts back into their repsected positions in rpart and ipart
         call update_findpts_info
c        Move particle info to the processor that owns each particle
c        using crystal router in log P time:

         jps = jpid1-1     ! Pointer to temporary proc id for swapping
         do i=1,n        ! Can't use jpt because it messes up particle info
            ipart(jps,i) = ipart(jpt,i)
         enddo
         call crystal_tuple_transfer(i_cr_hndl,n,llpart
     $              , ipart,ni,partl,nl,rpart,nr,jps)
c        Sort by element number - for improved local-eval performance
         call crystal_tuple_sort    (i_cr_hndl,n 
     $              , ipart,ni,partl,nl,rpart,nr,je0,1)
      endif

      ! end timer
      pttime(23) = pttime(23) + dnekclock() - ptdum(23)

      return
      end
c-----------------------------------------------------------------------
      subroutine particles_in_nid
      include 'SIZE'
      include 'CMTPART'

      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      integer icalld
      save    icalld
      data    icalld  /-1/

      ! begin timer
      ptdum(24) = dnekclock()

      icalld = icalld + 1

      nfpts = 0
      do ip = 1,n
         xloc = rpart(jx,ip)
         yloc = rpart(jy,ip)
         zloc = rpart(jz,ip)
         itest = 0
         if (nrect_assume .eq. 0) goto 1511
         do ie=1,nelt
            if (xloc.ge.xerange(1,1,ie).and.xloc.le.xerange(2,1,ie))then
            if (yloc.ge.xerange(1,2,ie).and.yloc.le.xerange(2,2,ie))then
            if (zloc.ge.xerange(1,3,ie).and.zloc.le.xerange(2,3,ie))then
                ipart(je0 ,ip) = ie-1
                if (icalld .eq. 0) ipart(je00,ip) = ie-1 ! set previous element as well
                ipart(jrc ,ip) = 0
                ipart(jpt ,ip) = nid
                rpart(jd  ,ip) = 1.0 
                rloc = -1.0 + 2.0*(xloc - xerange(1,1,ie))/
     $                 (xerange(2,1,ie)-xerange(1,1,ie))
                sloc = -1.0 + 2.0*(yloc - xerange(1,2,ie))/
     $                 (xerange(2,2,ie)-xerange(1,2,ie))
                tloc = -1.0 + 2.0*(zloc - xerange(1,3,ie))/
     $                 (xerange(2,3,ie)-xerange(1,3,ie))
                rpart(jr  ,ip) = rloc
                rpart(jr+1,ip) = sloc
                rpart(jr+2,ip) = tloc
                itest = 1
                goto 123
            endif
            endif
            endif
         enddo
 1511 continue
         if (itest.eq.0)then
            nfpts = nfpts + 1
            ifptsmap(nfpts) = ip
            call copy (rfpts(1,nfpts),rpart(1,ip),nrf) 
            call icopy(ifpts(1,nfpts),ipart(1,ip),nif) 
            if(nfpts.gt.llpart)then
               write(6,*)'Too many points crossing over ',
     $                      nfpts,llpart,nid
               call exitt
            endif
         endif
123      continue
      enddo

      ! end timer
      pttime(24) = pttime(24) + dnekclock() - ptdum(24)

      return
      end
c----------------------------------------------------------------------
      subroutine update_particle_location
c     check if particles are outside domain
c     > if bc_part = 0 then it is periodic
c     > if bc_part = -1,1 then particles are killed (outflow)
      include 'SIZE'
      include 'CMTDATA'
      include 'CMTPART'

      real  xdrange(2,3) 
      common /domainrange/ xdrange

      integer in_part(llpart), icount_p,itmp(li,llpart)
      real    rtmp(lr,llpart)

      ! begin timer
      ptdum(17) = dnekclock()

      jx0 = jx

      do i=1,n
         in_part(i) = 0
         do j=0,ndim-1
            if (rpart(jx0+j,i).lt.xdrange(1,j+1))then
c           if (bc_part(1).eq.0) then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(2,j+1) - 
     &                             abs(xdrange(1,j+1) - rpart(jx0+j,i))
                  rpart(jx1+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx1+j,i))
                  rpart(jx2+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx2+j,i))
                  rpart(jx3+j,i) = xdrange(2,j+1) +
     &                             abs(xdrange(1,j+1) - rpart(jx3+j,i))
                  goto 1512
c              elseif (bc_part(1).eq. 1) then
               elseif (((bc_part(1).ne.0) .and. (j.eq.0)) .or. ! outflow
     >                 ((bc_part(3).ne.0) .and. (j.eq.1)) .or.     
     >                 ((bc_part(5).ne.0) .and. (j.eq.2)) ) then
                  in_part(i) = -1
                  goto 1511
               endif
            endif
            if (rpart(jx0+j,i).gt.xdrange(2,j+1))then
c              if (bc_part(1).eq. 0) then
               if (((bc_part(1).eq.0) .and. (j.eq.0)) .or.   ! periodic
     >             ((bc_part(3).eq.0) .and. (j.eq.1)) .or.     
     >             ((bc_part(5).eq.0) .and. (j.eq.2)) ) then
                  rpart(jx0+j,i) = xdrange(1,j+1) +
     &                             abs(rpart(jx0+j,i) - xdrange(2,j+1))
                  rpart(jx1+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx1+j,i) - xdrange(2,j+1))
                  rpart(jx2+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx2+j,i) - xdrange(2,j+1))
                  rpart(jx3+j,i) = xdrange(1,j+1) -
     &                             abs(rpart(jx3+j,i) - xdrange(2,j+1))
                  goto 1512
c              elseif (bc_part(1).eq. 1) then
               elseif (((bc_part(1).ne.0) .and. (j.eq.0)) .or. ! outflow
     >                 ((bc_part(3).ne.0) .and. (j.eq.1)) .or.     
     >                 ((bc_part(5).ne.0) .and. (j.eq.2)) ) then
                  in_part(i) = -1
                  goto 1511
               endif
            endif
 1512 continue
         enddo
 1511 continue
      enddo

      nbc_sum = abs(bc_part(1)) + abs(bc_part(2)) + 
     >          abs(bc_part(3)) + abs(bc_part(4)) +
     >          abs(bc_part(5)) + abs(bc_part(6)) ! all periodic, don't search
      if (nbc_sum .gt. 0) then
      ic = 0
      do i=1,n
         if (in_part(i).eq.0) then
            ic = ic + 1 
            if (i .ne. ic) then
               call copy(rpart(1,ic),rpart(1,i),nr)
               call icopy(ipart(1,ic),ipart(1,i),ni)
            endif
         endif
      enddo
      n = ic
      endif

      ! end timer
      pttime(17) = pttime(17) + dnekclock() - ptdum(17)

      return
      end
c-----------------------------------------------------------------------
      subroutine update_findpts_info
      include 'SIZE'
      include 'CMTPART'

      ! begin timer
      ptdum(25) = dnekclock()

      do ifp = 1,nfpts
         call copy(rpart(1,ifptsmap(ifp)),rfpts(1,ifp),nrf)
         call icopy(ipart(1,ifptsmap(ifp)),ifpts(1,ifp),nif)
      enddo

      ! end timer
      pttime(25) = pttime(25) + dnekclock() - ptdum(25)

      return
      end
c-----------------------------------------------------------------------
      subroutine intpts_setup(tolin,ih)
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
      if(nidd.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end
c-----------------------------------------------------------------------
