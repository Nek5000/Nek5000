#ifndef NOMPIIO

      subroutine gfldr(sourcefld)
c
c     generic field file reader
c     reads sourcefld and interpolates all avaiable fields
c     onto current mesh
c
c     memory requirement: 
c     nelgs*nxs**ndim < np*(4*lelt*lx1**ldim)
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      character sourcefld*(*)

      common /scrcg/  pm1(lx1*ly1*lz1,lelv)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      character*132 hdr
      character*1   hdr1(132)
      equivalence   (hdr1,hdr)

      integer*8 dtmp8

      logical if_full_pres_tmp

      logical ifbswp, if_byte_swap_test
      real*4 bytetest


      etime_t = dnekclock_sync()
      if(nio.eq.0) write(6,*) 'call gfldr' 

      ! open source field file
      ierr = 0
      if(nid.eq.0) then
        open (90,file=sourcefld,status='old',err=100)
        close(90)
        goto 101
 100    ierr = 1
 101  endif
      call err_chk(ierr,' Cannot open source fld file!$')
      call byte_open_mpi(sourcefld,fldh_gfldr,.true.,ierr)

      ! read and parse header
      call byte_read_mpi(hdr,iHeaderSize/4,0,fldh_gfldr,ierr)
      call byte_read_mpi(bytetest,1,0,fldh_gfldr,ierr)

      call mfi_parse_hdr(hdr,ierr)
      call err_chk(ierr,' Invalid header!$')
      ifbswp = if_byte_swap_test(bytetest,ierr)
      call err_chk(ierr,' Invalid endian tag!$')

      nelgs   = nelgr
      nxs     = nxr
      nys     = nyr
      nzs     = nzr
      if(nzs.gt.1) then 
        ndims = 3
      else
        ndims = 2
      endif
      if (ifgtim) time = timer

      if_full_pres_tmp = if_full_pres
      if (wdsizr.eq.8) if_full_pres = .true.

      ! distribute elements across all ranks
      nels = nelgs/np
      do i = 0,mod(nelgs,np)-1
         if(i.eq.nid) nels = nels + 1
      enddo
      nxyzs      = nxs*nys*nzs
      dtmp8      = nels
      ntots_b    = dtmp8*nxyzs*wdsizr
      rankoff_b  = igl_running_sum(nels) - dtmp8
      rankoff_b  = rankoff_b*nxyzs*wdsizr  
      dtmp8      = nelgs
      nSizeFld_b = dtmp8*nxyzs*wdsizr
      noff0_b    = iHeaderSize + iSize + iSize*dtmp8

      ! do some checks
      if(ndims.ne.ndim) 
     $ call exitti('ndim of source does not match target!$',0)
      if(ntots_b/wdsize .gt. ltots) then
        dtmp8 = nelgs
        lelt_req = dtmp8*nxs*nys*nzs / (np*ltots/lelt)
        lelt_req = lelt_req + 1
        if(nio.eq.0) write(6,*)
     $   'ABORT: buffer too small, increase lelt > ', lelt_req
        call exitt
      endif

      ifldpos = 0
      if(ifgetxr) then
        ! read source mesh coordinates
        call gfldr_getxyz(xm1s,ym1s,zm1s,ifbswp)
        ifldpos = ndim
      else
        call exitti('source does not contain a mesh!$',0)
      endif

      ! initialize interpolation tool using source mesh
      nxf   = 2*nxs
      nyf   = 2*nys
      nzf   = 2*nzs
      nhash = nxs*nys*nzs 
      nmax  = 256

      call findpts_setup(inth_gfldr,nekcomm,np,ndim,
     &                   xm1s,ym1s,zm1s,nxs,nys,nzs,
     &                   nels,nxf,nyf,nzf,bb_t,
     &                   nhash,nhash,nmax,tol)

      ! read source fields and interpolate
      if(ifgetur .and. ifflow) then
        call gfldr_getfld(vx,vy,vz,ndim,ifldpos+1,ifbswp)
        ifldpos = ifldpos + ndim
      endif
      if(ifgetpr) then
        call gfldr_getfld(pm1,dum,dum,1,ifldpos+1,ifbswp)
        ifldpos = ifldpos + 1
        if (ifaxis) call axis_interp_ic(pm1)
        if (ifgetpr) call map_pm1_to_pr(pm1,1)
      endif
      if(ifgettr .and. ifheat) then
        call gfldr_getfld(t(1,1,1,1,1),dum,dum,1,ifldpos+1,ifbswp)
        ifldpos = ifldpos + 1
      endif
      do i = 1,ldimt-1
         if(ifgtpsr(i)) then
           call gfldr_getfld(t(1,1,1,1,i+1),dum,dum,1,ifldpos+1,ifbswp) 
           ifldpos = ifldpos + 1
         endif
      enddo

      if_full_pres = if_full_pres_tmp

      call byte_close_mpi(fldh_gfldr,ierr)
      call findpts_free(inth_gfldr)

      etime_t = dnekclock_sync() - etime_t
      if(nio.eq.0) write(6,'(A,1(1g8.2),A)')
     &                   ' done :: gfldr  ', etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_getxyz(xout,yout,zout,ifbswp)

      include 'SIZE'
      include 'GFLDR'
      include 'RESTART'

      real xout(*)
      real yout(*)
      real zout(*)
      logical ifbswp

      integer*8 ioff_b

 
      ioff_b = noff0_b + ndim*rankoff_b
      call byte_set_view(ioff_b,fldh_gfldr)

      nread = ndim*ntots_b/4
      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
      if(ifbswp) then
        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
      endif

      call gfldr_buf2vi (xout,1,bufr,ndim,wdsizr,nels,nxyzs)
      call gfldr_buf2vi (yout,2,bufr,ndim,wdsizr,nels,nxyzs)
      if(ndim.eq.3)
     $ call gfldr_buf2vi(zout,3,bufr,ndim,wdsizr,nels,nxyzs)

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_getfld(out1,out2,out3,nndim,ifldpos,ifbswp)

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'
      include 'RESTART'

      real out1(*)
      real out2(*)
      real out3(*)
      logical ifbswp

      integer*8 ioff_b

      logical ifpts

      integer icalld
      save    icalld
      data    icalld /0/


      ifpts = .false.
      if(icalld.eq.0) then
        ifpts = .true. ! find points
        icalld = 1
      endif

      ! read field data from source fld file
      ioff_b = noff0_b + (ifldpos-1)*nSizeFld_b
      ioff_b = ioff_b  + nndim*rankoff_b
      call byte_set_view(ioff_b,fldh_gfldr)
      nread = nndim*ntots_b/4
      call byte_read_mpi(bufr,nread,-1,fldh_gfldr,ierr)
      if(ifbswp) then
        if(wdsizr.eq.4) call byte_reverse (bufr,nread,ierr)
        if(wdsizr.eq.8) call byte_reverse8(bufr,nread,ierr)
      endif

      ! interpolate onto current mesh
      ntot = nx1*ny1*nz1*nelt
      call gfldr_buf2vi  (buffld,1,bufr,nndim,wdsizr,nels,nxyzs)
      call gfldr_intp    (out1,buffld,ifpts)
      if(nndim.eq.1) return

      call gfldr_buf2vi  (buffld,2,bufr,nndim,wdsizr,nels,nxyzs)
      call gfldr_intp    (out2,buffld,.false.)
      if(nndim.eq.2) return

      if(nndim.eq.3) then
        call gfldr_buf2vi(buffld,3,bufr,nndim,wdsizr,nels,nxyzs)
        call gfldr_intp  (out3,buffld,.false.)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_buf2vi(vi,index,buf,ldim,wds,nel,nxyz)

      real    vi(*)
      real*4  buf(*)
      integer wds


      do iel = 1,nel
         j = (iel-1)*nxyz
         k = (iel-1)*ldim*nxyz

         if(index.eq.2) k = k+nxyz
         if(index.eq.3) k = k+2*nxyz

         if(wds.eq.4) call copy4r(vi(j+1),buf(k+1)  ,nxyz)
         if(wds.eq.8) call copy  (vi(j+1),buf(2*k+1),nxyz)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_intp(fieldout,fieldin,iffpts)

      include 'SIZE'
      include 'GEOM'
      include 'GFLDR'

      real    fieldout(*),fieldin(*)
      logical iffpts

      integer*8 i8glsum,nfail,nfail_sum


      nfail = 0
      ntot  = nx1*ny1*nz1*nelt

      if(iffpts) then ! locate points (iel,iproc,r,s,t)

        call findpts(inth_gfldr,
     &               rcode,1,
     &               proc,1,
     &               elid,1,
     &               rst,ndim,
     &               dist,1,
     &               xm1,1,
     &               ym1,1,
     &               zm1,1,ntot)

        do i=1,ntot
           if(rcode(i).eq.1 .and. dist(i).gt.10*tol) nfail = nfail + 1
           if(rcode(i).eq.2) nfail = nfail + 1
        enddo

        nfail_sum = i8glsum(nfail,1)
        if(nio.eq.0 .and. nfail_sum.gt.0) write(6,*) 
     &   ' WARNING: Unable to find all mesh points in source fld ',
     &   nfail_sum

      endif

      ! evaluate inut field at given points
      call findpts_eval(inth_gfldr,
     &                  fieldout,1,
     &                  rcode,1,
     &                  proc,1,
     &                  elid,1,
     &                  rst,ndim,ntot,
     &                  fieldin)

      return
      end

#endif
