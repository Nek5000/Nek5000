c-----------------------------------------------------------------------
      subroutine read_re2_data(ifbswap, ifxyz, ifcur, ifbc)  ! .re2 reader

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'CTIMER'

      logical ifbswap
      logical ifxyz, ifcur, ifbc
      integer idummy(100)

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

c TODO: set as a .par file param      
      ifmpiio = .true.
#ifdef NOMPIIO
      ifmpiio = .false.
#endif

c TODO: hardcoding default 10, to set as a .par file param
      cbnodes = min(1024,np)
 
      etime0 = dnekclock_sync()

                  ibc = 2
      if (ifflow) ibc = 1

                  nfldt = 1
      if (ifheat) nfldt = 2+npscal
      if (ifmhd ) nfldt = 2+npscal+1

      ! first field to read
      if (param(33).gt.0) ibc = int(param(33))

      ! number of fields to read
      if (param(32).gt.0) then
        nfldt = ibc + int(param(32)) - 1
        nfldt = max(nfldt,1) 
        if (nelgt.gt.nelgv) nfldt = max(nfldt,2) 
      endif

      call blank(cbc,3*size(cbc))
      call rzero(bc ,size(bc))

      call fgslib_crystal_setup(cr_re2,nekcomm,np)
      call nek_file_open(nekcomm,re2fle,0,ifmpiio,cbnodes,re2_h,ierr)
      call err_chk(ierr,' Cannot open .re2 file!$')
      if (nid.eq.0) write(6,*) 'ifmpiio: ', ifmpiio
      if (nid.eq.0) write(6,*) 'calling readp_re2_*'
      call readp_re2_mesh (ifbswap,ifxyz)
      call readp_re2_curve(ifbswap,ifcur)
      do ifield = ibc,nfldt
      call readp_re2_bc(cbc(1,1,ifield),bc(1,1,1,ifield),
     &                  ifbswap,ifbc)
      enddo
      call nek_file_close(re2_h,ierr)
      call fgslib_crystal_free(cr_re2)



      etime_t = dnekclock_sync() - etime0
      if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: read .re2 file   ',
     &                   etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------
      subroutine readp_re2_mesh(ifbswap,ifread) ! version 2 of .re2 reader

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap
      logical ifread

      parameter(nrmax = lelt)             ! maximum number of records
      parameter(lrs   = 1+ldim*(2**ldim)) ! record size: group x(:,c) ...
      parameter(li    = 2*lrs+2)

      integer         bufr(li-2,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg

      nrg       = nelgt
      nr        = nelt
      irankoff  = igl_running_sum(nr) - nr
      dtmp8     = irankoff
      re2off_b  = 84 ! set initial offset (hdr + endian)
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      ! read coordinates from file
      nwds4r = nr*lrs4
      call nek_file_read(re2_h,int8(nwds4r),int8(lre2off_b),bufr,ierr)
      re2off_b = re2off_b + nrg*4*lrs4
      if (ierr.gt.0) goto 100

      if (.not.ifread) return

      if (nio.eq.0) write(6,*) 'reading mesh '

      ! pack buffer
      do i = 1,nr
         jj      = (i-1)*lrs4 + 1
         ielg    = irankoff + i ! elements are stored in global order
         vi(1,i) = gllnid(ielg)
         vi(2,i) = ielg
         call icopy(vi(3,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n   = nr
      key = 1 
      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,
     &   vl,0,vr,0,key)

      ! unpack buffer
      ierr = 0
      if (n.gt.nrmax) then
         ierr = 1
         goto 100
      endif

      do i = 1,n
         iel = gllel(vi(2,i)) 
         call icopy     (bufr,vi(3,i),lrs4)
         call buf_to_xyz(bufr,iel,ifbswap,ierr)
      enddo

 100  call err_chk(ierr,'Error reading .re2 mesh$')

      return
      end
c-----------------------------------------------------------------------
      subroutine readp_re2_curve(ifbswap,ifread)

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap
      logical ifread

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      parameter(nrmax = 12*lelt) ! maximum number of records
      parameter(lrs   = 2+1+5)   ! record size: eg iside curve(5) ccurve
      parameter(li    = 2*lrs+1)

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg
      integer*4       nrg4(2)
     
      integer*8       i8gl_running_sum 

      ! read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call nek_file_read(re2_h,int8(nwds4r),int8(lre2off_b),nrg4,ierr)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r

      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4

      if (.not.ifread) return
      if(nio.eq.0) write(6,*) 'reading curved sides '

      nwds4r = nr*lrs4
      call nek_file_read(re2_h,int8(nwds4r),int8(lre2off_b),bufr,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1,nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then 
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         ielg = bufr(jj,1)
         if(wdsizi.eq.8) call copyi4(ielg,bufr(jj,1),1)

         if(ielg.le.0 .or. ielg.gt.nelgt) goto 100
         vi(1,i) = gllnid(ielg)

         call icopy (vi(2,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1
      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,vl,0,vr,0,
     &                                   key)

      ! unpack buffer
      if(n.gt.nrmax) goto 100
      do i = 1,n
         call icopy       (bufr,vi(2,i),lrs4)
         call buf_to_curve(bufr)
      enddo

      return

 100  ierr = 1
      call err_chk(ierr,'Error reading .re2 curved data$')

      end
c-----------------------------------------------------------------------
      subroutine readp_re2_bc(cbl,bl,ifbswap,ifread)

      include 'SIZE'
      include 'TOTAL'

      character*3  cbl(  6,lelt)
      real         bl (5,6,lelt)
      logical      ifbswap
      logical      ifread

      parameter(nrmax = 6*lelt) ! maximum number of records
      parameter(lrs   = 2+1+5)  ! record size: eg iside bl(5) cbl
      parameter(li    = 2*lrs+1)

      integer         bufr(li-1,nrmax)
      common /scrns/  bufr

      integer         vi  (li  ,nrmax)
      common /ctmp1/  vi

      integer*8       lre2off_b,dtmp8
      integer*8       nrg
      integer*4       nrg4(2)

      integer*8       i8gl_running_sum 

      ! read total number of records
      nwds4r    = 1*wdsizi/4
      lre2off_b = re2off_b
      call nek_file_read(re2_h,int8(nwds4r),int8(lre2off_b),nrg4,ierr)
      call bcast(nrg4,wdsizi)
      if(ierr.gt.0) goto 100

      if(wdsizi.eq.8) then
         if(ifbswap) call byte_reverse8(nrg4,nwds4r,ierr)
         call copy(dnrg,nrg4,1)
         nrg = dnrg
      else
         if(ifbswap) call byte_reverse (nrg4,nwds4r,ierr)
         nrg = nrg4(1)
      endif
      re2off_b = re2off_b + 4*nwds4r
      if(nrg.eq.0) return

      ! read data from file
      dtmp8 = np
      nr = nrg/dtmp8
      do i = 0,mod(nrg,dtmp8)-1
         if(i.eq.nid) nr = nr + 1
      enddo
      dtmp8     = i8gl_running_sum(int(nr,8)) - nr
      lre2off_b = re2off_b + dtmp8*lrs*wdsizi
      lrs4      = lrs*wdsizi/4

      re2off_b = re2off_b + nrg*4*lrs4

      if (.not.ifread) return
      if(nio.eq.0) write(6,*) 'reading bc for ifld',ifield

      nwds4r = nr*lrs4
      call nek_file_read(re2_h,int8(nwds4r),int8(lre2off_b),bufr,ierr)
      if(ierr.gt.0) goto 100

      ! pack buffer
      do i = 1,nr
         jj = (i-1)*lrs4 + 1

         if(ifbswap) then 
           lrs4s = lrs4 - wdsizi/4 ! words to swap (last is char)
           if(wdsizi.eq.8) call byte_reverse8(bufr(jj,1),lrs4s,ierr)
           if(wdsizi.eq.4) call byte_reverse (bufr(jj,1),lrs4s,ierr)
         endif

         ielg = bufr(jj,1)
         if(wdsizi.eq.8) call copyi4(ielg,bufr(jj,1),1)

         if(ielg.le.0 .or. ielg.gt.nelgt) goto 100
         vi(1,i) = gllnid(ielg)

         call icopy (vi(2,i),bufr(jj,1),lrs4)
      enddo

      ! crystal route nr real items of size lrs to rank vi(key,1:nr)
      n    = nr
      key  = 1

      call fgslib_crystal_tuple_transfer(cr_re2,n,nrmax,vi,li,vl,0,vr,0,
     &                                   key)

      ! fill up with default
      do iel=1,nelt
      do k=1,6
         cbl(k,iel) = 'E  '
      enddo
      enddo

      ! unpack buffer
      if(n.gt.nrmax) goto 100
      do i = 1,n
         call icopy    (bufr,vi(2,i),lrs4)
         call buf_to_bc(cbl,bl,bufr)
      enddo

      return

 100  ierr = 1
      call err_chk(ierr,'Error reading .re2 boundary data$')

      end
c-----------------------------------------------------------------------
      subroutine buf_to_xyz(buf,e,ifbswap,ierr)! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'
      logical ifbswap

c      integer e,eg,buf(0:49)
      integer e,eg,buf(0:49)

      nwds = (1 + ldim*(2**ldim))*(wdsizi/4) ! group + 2x4 for 2d, 3x8 for 3d

      if     (ifbswap.and.ierr.eq.0.and.wdsizi.eq.8) then
          call byte_reverse8(buf,nwds,ierr)
      elseif (ifbswap.and.ierr.eq.0.and.wdsizi.eq.4) then
          call byte_reverse (buf,nwds,ierr)
      endif
      if(ierr.ne.0) return

      if(wdsizi.eq.8) then
         call copyi4(igroup(e),buf(0),1) !0-1
         if (ldim.eq.3) then
            call copy  (xc(1,e),buf( 2),8) !2 --17
            call copy  (yc(1,e),buf(18),8) !18--33
            call copy  (zc(1,e),buf(34),8) !34--49
         else
            call copy  (xc(1,e),buf( 2),4) !2 --9
            call copy  (yc(1,e),buf(10),4) !10--17
          endif
      else
         igroup(e) = buf(0)
         if (if3d) then
            call copy4r(xc(1,e),buf( 1),8)
            call copy4r(yc(1,e),buf( 9),8)
            call copy4r(zc(1,e),buf(17),8)
         else
            call copy4r(xc(1,e),buf( 1),4)
            call copy4r(yc(1,e),buf( 5),4)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_curve(buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2
        e  = gllel(eg)

        call copyi4(f,buf(3),1) !3-4

        call copy  ( curve(1,f,e),buf(5) ,5) !5--14
        call chcopy(ccurve(  f,e),buf(15),1)!15
      else
        eg = buf(1)
        e  = gllel(eg)
        f  = buf(2)

        call copy4r( curve(1,f,e),buf(3),5)
        call chcopy(ccurve(f,e)  ,buf(8),1)
      endif

c     write(6,1) eg,e,f,(curve(k,f,e),k=1,5),ccurve(f,e)
c   1 format(2i7,i3,5f10.3,1x,a1,'ccurve')

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_to_bc(cbl,bl,buf)    ! version 1 of binary reader

      include 'SIZE'
      include 'TOTAL'

      character*3 cbl(6,lelt)
      real         bl(5,6,lelt)

      integer e,eg,f,buf(30)

      if(wdsizi.eq.8) then
        call copyi4(eg,buf(1),1) !1-2
        e  = gllel(eg)

        call copyi4(f,buf(3),1) !3-4

        call copy  (bl(1,f,e),buf(5),5) !5--14
        call chcopy(cbl( f,e),buf(15),3)!15-16

        if(nelt.ge.1000000.and.cbl(f,e).eq.'P  ')
     $   call copyi4(bl(1,f,e),buf(5),1) !Integer assign connecting P element

      else
        eg = buf(1)
        e  = gllel(eg)
        f  = buf(2)

        call copy4r ( bl(1,f,e),buf(3),5)
        call chcopy (cbl(  f,e),buf(8),3)

        if (nelgt.ge.1000000.and.cbl(f,e).eq.'P  ')
     $     bl(1,f,e) = buf(3) ! Integer assign of connecting periodic element
      endif

c      write(6,1) eg,e,f,cbl(f,e),' CBC',nid
c  1   format(2i8,i4,2x,a3,a4,i8)

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_close_outv  ! this is the stupid O(P) formulation

      include 'SIZE'
      include 'PARALLEL'
      integer*4 zero
      real      rzero

      len   = wdsizi
      rzero = 0
      zero  = 0
c     write(6,*) nid,' bufclose'
      if (nid.eq.0) then
         do mid=1,np-1
            if(wdsizi.eq.8)call csend(mid,rzero,len,mid,0)
            if(wdsizi.eq.4)call csend(mid, zero,len,mid,0)
c           write(6,*) mid,' sendclose'
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine buf_close_out  ! this is the stupid O(P) formulation

      include 'SIZE'
      include 'PARALLEL'
      integer*4 zero
      real      rzero

c     len  = 4
      len   = wdsizi
      zero = 0
      rzero = 0
      if (nid.eq.0) then
         do mid=1,np-1
            if(wdsizi.eq.8)call csend(mid,rzero,len,mid,0)
            if(wdsizi.eq.4)call csend(mid, zero,len,mid,0)
         enddo
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_re2_hdr(ifbswap, ifverbose) ! open file & chk for byteswap

      include 'SIZE'
      include 'TOTAL'

      logical ifbswap, ifverbose
      logical if_byte_swap_test

      integer fnami (33)
      character*132 fname
      equivalence (fname,fnami)

      character*132 hdr
      character*5 version
      real*4      test

      logical iffound

      ierr=0

      if (nid.eq.0) then
         if (ifverbose) write(6,'(A,A)') ' Reading ', re2fle
         call izero(fnami,33)
         m = indx2(re2fle,132,' ',1)-1
         call chcopy(fname,re2fle,m)

         inquire(file=fname, exist=iffound)
         if(.not.iffound) ierr = 1
      endif
      call err_chk(ierr,' Cannot find re2 file!$')

      if (nid.eq.0) then
         call byte_open(fname,ierr)
         if(ierr.ne.0) goto 100
         call byte_read(hdr,20,ierr)
         if(ierr.ne.0) goto 100

         read (hdr,1) version,nelgt,ldimr,nelgv
    1    format(a5,i9,i3,i9)
 
         wdsizi = 4
         if(version.eq.'#v002') wdsizi = 8
         if(version.eq.'#v003') then
           wdsizi = 8
           param(32)=1
         endif

         call byte_read(test,1,ierr)
         if(ierr.ne.0) goto 100
         ifbswap = if_byte_swap_test(test,ierr)
         if(ierr.ne.0) goto 100
        
         call byte_close(ierr)
      endif
 
 100  call err_chk(ierr,'Error reading re2 header$')

      call bcast(wdsizi, ISIZE)
      call bcast(ifbswap,LSIZE)
      call bcast(nelgv  ,ISIZE)
      call bcast(nelgt  ,ISIZE)
      call bcast(ldimr  ,ISIZE)
      call bcast(param(32),WDSIZE)

      if(wdsize.eq.4.and.wdsizi.eq.8) 
     $   call exitti('wdsize=4 & wdsizi(re2)=8 not compatible$',wdsizi)

      return
      end
