c-----------------------------------------------------------------------
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
      if (loglevel .gt. 2) ifverbm=.true.

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
           if (loglevel .gt. 2) then
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
      common /ctmp0/ iwork(lelt)
c
      REAL*8 dnekclock,t0
c
      t0 = dnekclock()
c     if (.not.(ifgtp.or.ifgfdm)) then
      if (.not.ifgtp) then
c
c        rsb element to processor mapping 
c
         if (ifgfdm)       call gfdm_elm_to_proc(gllnid,np) ! gfdm w/ .map

         call get_map

      endif

      if(ifzper.or.ifgtp) call gfdm_elm_to_proc(gllnid,np) ! special processor map

c     compute global to local map (no processor info)
c
      IEL=0
      CALL IZERO(GLLEL,NELGT)
      DO IEG=1,NELGT
         IF (GLLNID(IEG).EQ.NID) THEN
            IEL = IEL + 1
            GLLEL(IEG)=IEL
            NELT = IEL
            if (ieg.le.nelgv) NELV = IEL
         ENDIF
c        write(6,*) 'map2 ieg:',ieg,nelv,nelt,nelgv,nelgt
      ENDDO
c
c     dist. global to local map to all processors
c
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllel(k),iwork,'+  ',m)
         k = k+m
      enddo
c
c     compute local to global map
c     (i.e. returns global element number given local index and proc id)
c
      do ieg=1,nelgt
         mid  =gllnid(ieg)
         ie   =gllel (ieg)
         if (mid.eq.nid) lglel(ie)=ieg
      enddo
c
c     All Done.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine gfdm_elm_to_proc(gllnid,np)
c
c
      include 'SIZE'
      include 'ZPER'
c
      integer gllnid(1)
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
      call gfdm_build_global_el_map (gllnid,map_st,nes,net
     $                                     ,nelbox,nstride_box,ip,is,it)
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
      subroutine gfdm_build_global_el_map (gllnid,map_st,nes,net
     $                                     ,nelbox,nstride_box,ip,is,it)
c
      include 'SIZE'
      integer gllnid(1)
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
            gllnid(ieg) = proc
         enddo
      enddo
      enddo
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
            read(80,*,err=100) neli,nnzi
         endif
      endif
 
      call bcast(neli, ISIZE)

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
         do ipass=1,npass
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
               if (.not.ifgfdm) gllnid(eg) = wk(1,m)  ! must still be divided
               wk(mdw,m) = eg
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

      if (.not.ifgfdm) then ! gllnid is already assigned for gfdm
        lng = isize*neli
        call bcast(gllnid,lng)
        call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch
      endif

      nelt=0 !     Count number of elements on this processor
      nelv=0
      do eg=1,neli
         if (gllnid(eg).eq.nid) then
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

      do i=1,ntuple
         eg=wk(mdw,i)
         wk(1,i)=gllnid(eg)        ! processor id for element eg
      enddo

      key = 1  ! processor id is in wk(1,:)
      call fgslib_crystal_ituple_transfer(cr_h,wk,mdw,ntuple,ndw,key)

      if (.not.ifgfdm) then            ! no sorting for gfdm?
         key = mdw  ! Sort tuple list by eg
         nkey = 1
         call fgslib_crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)
      endif

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
      subroutine assign_gllnid(gllnid,iunsort,nelgt,nelgv,np)
c
      integer gllnid(1),iunsort(1),nelgt,np 
      integer e,eg


      log2p = log2(np)
      np2   = 2**log2p
      if (np2.eq.np.and.nelgv.eq.nelgt) then   ! std power of 2 case

         npstar = ivlmax(gllnid,nelgt)+1
         nnpstr = npstar/np
         do eg=1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

         return

      elseif (np2.eq.np) then   ! std power of 2 case, conjugate heat xfer

c        Assign fluid elements
         npstar = max(np,ivlmax(gllnid,nelgv)+1)
         nnpstr = npstar/np
         do eg=1,nelgv
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

c        Assign solid elements
         nelgs  = nelgt-nelgv  ! number of solid elements
         npstar = max(np,ivlmax(gllnid(nelgv+1),nelgs)+1)
         nnpstr = npstar/np
         do eg=nelgv+1,nelgt
            gllnid(eg) = gllnid(eg)/nnpstr
         enddo

         return

      elseif (nelgv.ne.nelgt) then
         call exitti
     $       ('Conjugate heat transfer requires P=power of 2.$',np)
      endif


c  Below is the code for P a non-power of two:

c  Split the sorted gllnid array (read from .map file) 
c  into np contiguous partitions. 

c  To load balance the partitions in case of mod(nelgt,np)>0 
c  add 1 contiguous entry out of the sorted list to NODE_i 
c  where i = np-mod(nelgt,np) ... np


      nel   = nelgt/np       ! number of elements per processor
      nmod  = mod(nelgt,np)  ! bounded between 1 ... np-1
      npp   = np - nmod      ! how many paritions of size nel 
 
      ! sort gllnid  
      call isort(gllnid,iunsort,nelgt)

      ! setup partitions of size nel 
      k   = 0
      do ip = 0,npp-1
         do e = 1,nel  
            k = k + 1 
            gllnid(k) = ip
         enddo
      enddo
      ! setup partitions of size nel+1
      if(nmod.gt.0) then 
        do ip = npp,np-1
           do e = 1,nel+1  
              k = k + 1 
              gllnid(k) = ip
           enddo
        enddo 
      endif

      ! unddo sorting to restore initial ordering by
      ! global element number
      call iswapt_ip(gllnid,iunsort,nelgt)

      return
      end
