      subroutine mapelpr()

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SCRCT'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
c
      logical ifverbm
c
      etime0 = dnekclock_sync()

      if(nio.eq.0) write(6,'(A)') ' mapping elements to processors'

      MFIELD=2
      IF (IFFLOW) MFIELD=1
      IF (IFMVBD) MFIELD=0

c     Set up TEMPORARY value for NFIELD - NFLDT
      NFLDT = 1
      IF (IFHEAT) NFLDT = 2 + NPSCAL
c
c     Distributed memory processor mapping
      IF (NP.GT.NELGT) THEN
         IF(NID.EQ.0) THEN
           WRITE(6,1000) NP,NELGT
 1000      FORMAT(2X,'ABORT: Too many processors (',I8
     $          ,') for to few elements (',I8,').'
     $          ,/,2X,'ABORTING IN MAPELPR.')
         ENDIF
         call exitt
      ENDIF

      call set_proc_map()
c
      DO 1200 IFIELD=MFIELD,NFLDT
         IF (IFTMSH(IFIELD)) THEN
            NELG(IFIELD)      = NELGT
         ELSE
            NELG(IFIELD)      = NELGV
         ENDIF
 1200 CONTINUE

C     Output the processor-element map:
      ifverbm=.true.
      if (np.gt.2000.or.nelgt.gt.40000) ifverbm=.false.
      if (loglevel .gt. 0) ifverbm=.true.

      if(ifverbm) then
        idum = 1
        if(nid.eq.0) then
           N8 = min(8,nelt)
           write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
           if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
           DO inid=1,NP-1
              mtype = inid
              call csend(mtype,idum,4,inid,0)            ! handshake
              call crecv(mtype,inelt,4)               ! nelt of other cpus
              N8 = min(8,inelt)
           ENDDO
 1310      FORMAT(' RANK',I6,' IEG',8I8)
 1315      FORMAT('     ',6X,'    ',8I8)
        else
           mtype = nid
           call crecv(mtype,idum,4)                ! hand-shake
           call csend(mtype,nelt,4,0,0)            ! nelt
           if (loglevel .gt. 0) then
              N8 = min(8,nelt)
              write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
              if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
           endif
        endif
      endif

C     Check elemental distribution
C
C      IF (IPASS.EQ.2.AND.PARAM(156).eq.9) THEN
C         NXYZ=lx1*ly1*lz1
C         DO 1400 IE=1,NELT
C            VTMP1=NODE
c            VTMP2=IE
C            CALL CFILL(VX(1,1,1,IE) ,VTMP1,NXYZ)
C            CALL CFILL(VY(1,1,1,IE) ,VTMP2,NXYZ)
C            CALL CFILL(T(1,1,1,IE,1),VTMP1,NXYZ)
C 1400    CONTINUE
C         call prepost(.true.,'   ')
C      ENDIF

      nn = iglmin(nelt,1)
      nm = iglmax(nelt,1)
      dt = dnekclock() - etime0
      if(nio.eq.0) then
        write(6,*) ' '
        write(6,*) 'element load imbalance: ',nm-nn,nn,nm
        if((nm-nn)/float(nn).gt.0.2) 
     $    write(6,*) 'WARNING: imbalance >20% !!!'
        write(6,'(A,g13.5,A,/)')  ' done :: mapping ',dt,' sec'
        write(6,*) ' '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_proc_map()
C
C     Compute element to processor distribution according to (weighted) 
C     physical distribution in an attempt to minimize exposed number of
C     element interfaces.
C
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'SCRCT'
      include 'TSTEP'
      include 'ZPER'

      dProcmapCache = .false.

      if (.not.ifgtp) then
c        rsb element to processor mapping 
         if (ifgfdm)       call gfdm_elm_to_proc(np) ! gfdm w/ .map
         call get_map
      endif

      if(ifzper.or.ifgtp) call gfdm_elm_to_proc(np) ! special processor map

c     compute global to local map (no processor info)
      iel=0
      do ieg=1,nelgt
         if (nid.eq.0) mid = gllnid(ieg)
         call bcast(mid,isize)
         if (mid.eq.nid) then
            iel = iel + 1
            lglel(iel) = ieg
            call dProcmapPut(iel,1,1,ieg)
            nelt = iel
            if (ieg.le.nelgv) nelv = iel
         endif
      enddo
      call nekgsync() ! wait for pending puts

      dProcmapCache = .true.

      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_elm_to_proc(np)
c
c
      include 'SIZE'
      include 'ZPER'
c
      common /ctmp1/  map_st(lelg_sm)
      common /vptsol/ iwork(0:lp)
      integer nelbox(3),nstride_box(3)
c
      call gfdm_set_pst(ip,is,it,nelbox,nstride_box,lx2,ly2,lz2)
c
      nep = nelbox(ip)
      nes = nelbox(is)
     
      if(nelbox(it).eq.0) nelbox(it)=1
      net = nelbox(it)
c
      nst = nes*net
      if (nst.lt.np) then
         if (nid.eq.0) 
     $   write(6,*) 'ERROR, number of elements in plane must be > np'
     $   ,nst,np,nep,nes,net
         call exitt
      endif
c
c
      call gfdm_map_2d(map_st,nes,net,iwork,np)
      call gfdm_build_global_el_map (map_st,nes,net,
     $                               nelbox,nstride_box,ip,is,it)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_map_2d(map_st,nes,net,num_el,np)
c
c     Set up a 2D processor decomposition of an NES x NET array.
c
      integer map_st(nes,net),num_el(0:np)
c
c     First a stupid dealing algorithm to determine how many
c     elements on each processor

      do i=0,np
         num_el(i) = 0
      enddo
c
      k = np-1
      do i=1,nes*net
         num_el(k) = num_el(k)+1
         k=k-1
         if (k.lt.0) k = np-1
      enddo
c
      jnid = 0
      nel_cnt = 0
      nel_cur = num_el(jnid)
      do j=1,net,2
         do i=1,nes                 ! Count down
            nel_cnt = nel_cnt + 1
            if (nel_cnt.gt.nel_cur) then
               jnid=jnid+1
               nel_cur = num_el(jnid)
               nel_cnt = 1
            endif
            map_st(i,j) = jnid
         enddo
c
         j1 = j+1
         if (j1.le.net) then
            do i=nes,1,-1                ! Count up
               nel_cnt = nel_cnt + 1
               if (nel_cnt.gt.nel_cur) then
                  jnid=jnid+1
                  nel_cur = num_el(jnid)
                  nel_cnt = 1
               endif
               map_st(i,j1) = jnid
            enddo
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_set_pst(ip,is,it,nelbox,nstride_box,nxp,nyp,nzp)
c
      include 'SIZE'
      include 'INPUT'
      include 'ZPER'
c
      integer nelbox(3),nstride_box(3)
c
      if (if3d) then
         if (param(118).lt.0) then
            ip = 3
            is = 1
            it = 2
         elseif (param(117).lt.0) then
            ip = 2
            is = 3
            it = 1
         else
            ip = 1
            is = 2
            it = 3
         endif
      else
         if (param(117).lt.0) then
            ip = 2
            is = 1
            it = 3
         else
            ip = 1
            is = 2
            it = 3
         endif
      endif
c
      pst2lex(1)=ip       ! identify x-, y- or z with primary direction
      pst2lex(2)=is
      pst2lex(3)=it
c
      lex2pst(ip)=1
      lex2pst(is)=2
      lex2pst(it)=3
c
      nelbox(1) = nelx
      nelbox(2) = nely
      nelbox(3) = nelz
c
      nstride_box(1) = 1
      nstride_box(2) = nelx
      nstride_box(3) = nelx*nely
c
      ngfdm_p(1) = nelx*nxp
      ngfdm_p(2) = nely*nyp
      ngfdm_p(3) = nelz*nzp
      write(6,*) 'ngfdm:',(ngfdm_p(k),k=1,3)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_build_global_el_map (map_st,nes,net,
     $                                     nelbox,nstride_box,ip,is,it)
c
      include 'SIZE'
      integer map_st(nes,net)
      integer nelbox(3),nstride_box(3)
c
      integer proc
c
      do jt=1,nelbox(it)
      do js=1,nelbox(is)
         proc = map_st(js,jt)
         do jp=1,nelbox(ip)
            ieg = 1 + nstride_box(ip)*(jp-1)  ! nstride_p=nes*net
     $              + nstride_box(is)*(js-1)
     $              + nstride_box(it)*(jt-1)
            call dProcmapPut(proc,1,2,ieg)
         enddo
      enddo
      enddo
      call nekgsync() ! wait pending puts
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outmati(u,m,n,name6)
      integer u(m,n)
      character*6 name6
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
c     Print out copies of a global matrix
c
      do mid=0,np-1
        call nekgsync
        if (mid.eq.nid) then
         n20 = min(n,20)
         write(6,1) nid,m,n,name6
   1     format(//,3i6,'  Matrix:',2x,a6,/)
         do i=1,m
            write(6,2) nid,name6,(u(i,j),j=1,n20)
         enddo
   2     format(i3,1x,a6,20i6)
        endif
        call nekgsync
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine outmati8(u,m,n,name6)
      integer*8 u(m,n)
      character*6 name6
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal
c
c     Print out copies of a global matrix
c
      do mid=0,np-1
        call nekgsync
        if (mid.eq.nid) then
         n20 = min(n,20)
         write(6,1) nid,m,n,name6
   1     format(//,3i6,'  Matrix:',2x,a6,/)
         do i=1,m
            write(6,2) nid,name6,(u(i,j),j=1,n20)
         enddo
   2     format(i3,1x,a6,20i6)
        endif
        call nekgsync
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine get_map

      call get_vert

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert
c
c     Distribute and assign partitions using the .map file
c
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      parameter(mdw=2+2**ldim)
      parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
      common /scrns/ wk(mdw,ndw)
      integer wk

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer vertex

      integer icalld
      save    icalld
      data    icalld  /0/

      if (icalld.gt.0) return
      icalld = 1

      ncrnr = 2**ldim
      call get_vert_map(vertex,ncrnr,wk,mdw,ndw,ifgfdm)

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_map(vertex,nlv,wk,mdw,ndw,ifgfdm)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer vertex(nlv,1)
      integer wk(mdw,ndw)
      logical ifgfdm

      logical ifbswap,if_byte_swap_test

      character*132 mapfle
      character*1   mapfle1(132)
      equivalence  (mapfle,mapfle1)

      character*4   suffix
      character*1   suffix1(4)
      equivalence  (suffix,suffix1)

      character*132 hdr
      character*5   version
      real*4        test

      logical ifma2,ifmap
      integer e,eg,eg0,eg1
      integer itmp20(20)

      integer ibuf(2)

      ierr   = 0
      ifma2  = .false.
      suffix = '.map'

      if (nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (mapfle,132)
         call chcopy(mapfle,reafle,lfname)
         call chcopy(mapfle1(lfname+1),suffix,4)
         inquire(file=mapfle, exist=ifmap)

         if (.not.ifmap) then
            suffix = '.ma2'
            call chcopy(mapfle1(lfname+1),suffix,4)
            inquire(file=mapfle, exist=ifma2)
         endif

        if(.not.ifmap .and. .not.ifma2) ierr=1 
      endif
      if(nid.eq.0) write(6,'(A,A)') ' Reading ', mapfle
      call err_chk(ierr,' Cannot find map file!$')
      call bcast(ifma2,lsize)

      if (nid.eq.0) then
         if (ifma2) then         
            call byte_open(mapfle,ierr)
            if(ierr.ne.0) goto 100

            call blank(hdr,132)
            call byte_read(hdr,132/4,ierr)
            if(ierr.ne.0) goto 100

            read (hdr,1) version,neli,nnzi
    1       format(a5,2i12)

            call byte_read(test,1,ierr)
            if(ierr.ne.0) goto 100
            ifbswap = if_byte_swap_test(test,ierr)
            if(ierr.ne.0) goto 100
         else
            open(unit=80,file=mapfle,status='old',err=100)
            read(80,1,err=100) version,neli,nnzi
         endif
      endif
      call bcast(version, 5*CSIZE)
      call bcast(neli, ISIZE)

      if (version .ne. '#v002') then
         ierr = 1
         if (nid.eq.0) write(6,*) 
     $                 'Unsupported map file version - rerun genmap!'
         goto 200
      endif
        
      npass = 1 + (neli/ndw)
      if (npass.gt.np) then
         if (nid.eq.0) write(6,*) npass,np,neli,ndw,'Error get_vert_map'
         call exitt
      endif 

      len = 4*mdw*ndw
      if (nid.gt.0.and.nid.lt.npass) msg_id=irecv(nid,wk,len)
      call nekgsync

      if (nid.eq.0) then
         eg0 = 0
         do ipass=1,npass ! loop over all "ranks"
            eg1 = min(eg0+ndw,neli)

            if (ifma2) then
               nwds = (eg1 - eg0)*(mdw-1)
               call byte_read(wk,nwds,ierr)
               if (ierr.ne.0) goto 200
               if (ifbswap) call byte_reverse(wk,nwds,ierr)

               m = eg1 - eg0
               do eg=eg1,eg0+1,-1 ! reshuffle array
                  jj = (m-1)*(mdw-1) + 1
                  call icopy(itmp20,wk(jj,1),mdw-1)
                  call icopy(wk(1,m),itmp20 ,mdw-1)
                  m = m - 1
               enddo
            else
               m = 0
               do eg=eg0+1,eg1
                  m = m + 1
                  read(80,*,err=200) (wk(k,m),k=1,mdw-1)
               enddo
            endif
            
            m = 0
            do eg=eg0+1,eg1
               m = m + 1
               wk(mdw,m) = wk(1,m) !eg
               call dProcmapPut(wk(mdw,m),1,1,eg) ! store rank-sorted array
            enddo
    
            if (ipass.lt.npass) call csend(ipass,wk,len,ipass,0) !send to ipass
            eg0 = eg1
         enddo

         ntuple = m

         if (ifma2) then
            call byte_close(ierr)
         else
            close(80)
         endif

      elseif (nid.lt.npass) then

         call msgwait(msg_id)
         ntuple = ndw

      else

         ntuple = 0

      endif

      call nekgsync()

      if (.not.ifgfdm) call assign_gllnid()

      do i = 1,ntuple
         eg = wk(mdw,i)
         wk(1,i) = gllnid(eg)
      enddo

      nelt = 0 !     Count number of elements on this processor
      nelv = 0
      do eg=1,neli
         if (nid.eq.0) mid = gllnid(eg)
         call bcast(mid,isize)
         if (mid.eq.nid) then
            if (eg.le.nelgv) nelv=nelv+1
            if (eg.le.nelgt) nelt=nelt+1
         endif
      enddo

      if (np.le.64) write(6,*) nid,nelv,nelt,nelgv,nelgt,' NELV'

c     NOW: crystal route vertex by processor id

      ntuple_sum = iglsum(ntuple,1)
      if (ntuple_sum .ne. nelgt) then
         if (nid.eq.0) write(6,*) 'Error invalid tuple sum!'
         call exitt
      endif 

      key = 1  ! processor id is in wk(1,:)
      call fgslib_crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)

c      if (.not.ifgfdm) then            ! no sorting for gfdm?
         key = mdw  ! Sort tuple list by eg
         nkey = 1
         call fgslib_crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
c      endif

      iflag = 0
      if (ntuple.ne.nelt) then
         write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
         write(6,*) 'Check that .map file and .rea file agree'
         iflag=1
      else
         do e=1,nelt
            call icopy(vertex(1,e),wk(2,e),nlv)
         enddo
      endif

      iflag = iglmax(iflag,1)
      if (iflag.gt.0) then
         do mid=0,np-1
            call nekgsync
            if (mid.eq.nid)
     $      write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FB'
            call nekgsync
         enddo
         call nekgsync
         call exitt
      endif

      return

 100  continue
      call err_chk(ierr,'Error opening or reading map header$')

 200  continue
      call err_chk(ierr,'Error while reading map file$')

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapPut(ibuf,lbuf,ioff,ieg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'

      integer ibuf(lbuf)
      integer*8 disp

      call dProcMapFind(iloc,nids,ieg)
      disp = 2*(iloc-1) + ioff-1

      call mpi_win_lock(MPI_LOCK_EXCLUSIVE,nids,0,win,ierr)
      call mpi_put(ibuf,lbuf,MPI_INTEGER,nids,disp,1,MPI_INTEGER,
     $             win,ierr)
      call mpi_win_unlock(nids,win,ierr)

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcmapGet(ibuf,ieg)

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'

      integer ibuf(2)

      integer*8 disp

      parameter (lc = 128) ! cache size for remote entries
      parameter (lcache = lelt+lc+8-mod(lelt+lc,8)) ! multiple of 8
      integer   cache(lcache,3)
      save      cache

      save icalld
      data icalld /0/

      save iran
      parameter(im = 6075, ia = 106, ic = 1283)

      if (icalld .eq. 0) then
         do i = 1,lelt+lc
            cache(i,1) = -1
         enddo
         icalld = 1
      endif

      ii = lsearch_ur(cache,lcache,ieg)
      if (ii.gt.0 .and. ii.ne.lelt+lc) then ! cache hit
c         write(6,*) nid, 'cache hit ', 'ieg:', ieg
         ibuf(1) = cache(ii,2)
         ibuf(2) = cache(ii,3)
      else
         call dProcmapFind(il,nidt,ieg)
         disp = 2*(il-1)
         call mpi_win_lock(MPI_LOCK_SHARED,nidt,0,win,ierr)
         call mpi_get(ibuf,2,MPI_INTEGER,nidt,disp,2,MPI_INTEGER,
     $                win,ierr)
         call mpi_win_unlock(nidt,win,ierr)

         if (dProcmapCache) then
            ii = ibuf(1)
            if (ibuf(2).ne.nid) then
               iran = mod(iran*ia+ic,im)
               ii = lelt + (lc*iran)/im + 1 ! randomize array location 
            endif
            cache(ii,1) = ieg
            cache(ii,2) = ibuf(1)
            cache(ii,3) = ibuf(2)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dProcMapFind(il,nids,ieg)

      include 'SIZE'
      include 'PARALLEL'

      ! distribute array in blocks across ranks
      nstar = nelgt/np
      nids = (ieg-1)/nstar
      il = ieg - nids * nstar
      if (ieg .gt. np*nstar) then
         nids = mod(ieg,np) - 1
         il = nstar + 1
      endif

      return
      end
c-----------------------------------------------------------------------
      integer function lsearch_ur(a, n, k)

      integer a(n), n, k

      parameter(lvec=8) ! unroll factor

      slsearch = 0
      do i = 1,n,lvec
         do j = 0,lvec-1
            if (a(i+j).eq.k) slsearch = i + j
         enddo
         if (slsearch.gt.0) goto 10
      enddo

10    continue
      end
c-----------------------------------------------------------------------
      integer function gllnid(ieg)

      include 'mpif.h'

      integer iegl, nidl
      save    iegl, nidl
      data    iegl, nidl /0,0/

      integer ibuf(2)

      if (ieg .eq. iegl) then
         ibuf(2) = nidl
         goto 100
      endif
      call dProcmapGet(ibuf,ieg)

 100  iegl   = ieg
      nidl   = ibuf(2)
      gllnid = ibuf(2)

      end
c-----------------------------------------------------------------------
      integer function gllel(ieg)

      include 'mpif.h'

      integer iegl, iell
      save    iegl, iell
      data    iegl, iell /0,0/

      integer ibuf(2)

      if (ieg .eq. iegl) then
         ibuf(1) = iell
         goto 100
      endif
      call dProcmapGet(ibuf,ieg)

 100  iegl  = ieg
      iell  = ibuf(1)
      gllel = ibuf(1)

      end
c-----------------------------------------------------------------------
      subroutine assign_gllnid()
c
c     pratition rank sorted element list into approximately equal-sized 
c     chunks
c
      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'

      integer ibuf(2)
      integer nel(2)

      nel(1) = nelgv
      nel(2) = nelgt - nelgv

      if (nid .eq. 0) then
         do imsh = 1,2 
            n  = nel(imsh)/np
            nn = np - mod(nel(imsh),np)
            do ir = 1,nel(imsh) ! sweep through rank sorted array
               if (ir .le. nn*n) then
                  nid_el = (ir-1)/n
               else
                  nid_el = nn + (ir - nn*n - 1)/(n+1)
               endif
 
               ! store nid for ieg
               irr = ir
               if(imsh .eq. 2) irr = nelgv + irr 
               call dProcmapGet(ibuf,irr)
               ieg = ibuf(1) 
               call dProcmapPut(nid_el,1,2,ieg)
            enddo
         enddo
      endif
      call nekgsync()

      return
      end
