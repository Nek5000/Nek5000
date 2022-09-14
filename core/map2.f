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
      if(nio.eq.0 .and. loglevel.gt.1) write(6,'(A)') 
     $  ' partioning elements to MPI ranks'

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
     $          ,') for to few elements (',I12,').'
     $          ,/,2X,'ABORTING IN MAPELPR.')
         ENDIF
         call exitt
      ENDIF

      call set_proc_map()

      if (nelt.gt.lelt) then
         call exitti('nelt > lelt, increase lelt!$',nelt)
      endif

      DO 1200 IFIELD=MFIELD,NFLDT
         IF (IFTMSH(IFIELD)) THEN
            NELG(IFIELD)      = NELGT
         ELSE
            NELG(IFIELD)      = NELGV
         ENDIF
 1200 CONTINUE

C     Output the processor-element map:
      ifverbm=.false.
      if (loglevel .gt. 2) ifverbm=.true.

      if(ifverbm) then
        idum = 1
        if(nid.eq.0) then
           N8 = min(8,nelt)
           write(6 ,1310) node-1,(lglel(ie),ie=1,n8)
           if (NELT.GT.8) write(6 ,1315) (lglel(ie),ie=9,NELT)
           DO inid=1,NP-1
              mtype = inid
              call csend(mtype,idum,4,inid,0)         ! handshake
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

      dtmp = dnekclock_sync() - etime0
      if(nio.eq.0 .and. loglevel .gt. 1) then
        write(6,'(A,g13.5,A,/)')  ' done :: partioning ',dtmp,' sec'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine mapelpr_big()

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SCRCT'
      include 'SOLN'
      include 'TSTEP'
      include 'CTIMER'
c
      logical ifverbm
      integer ibuf(2), loc_to_glob_nid(lelt), lglelo(lelt)
c
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      etime0 = dnekclock_sync()
      if(nio.eq.0 .and. loglevel.gt.1) write(6,'(A)')
     $  ' partioning elements to MPI ranks'

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
     $          ,') for to few elements (',I12,').'
     $          ,/,2X,'ABORTING IN MAPELPR.')
         ENDIF
         call exitt
      ENDIF

      call get_vert_big(loc_to_glob_nid)

      call fgslib_crystal_setup(cr_re2,nekcomm,np)

      ! TODO: transfer vertices
c     do e = 1, nelt
c       lglelo(e) = lglel(e)
c     enddo
c     call transfer_re2_mesh(loc_to_glob_nid)

      ! TODO: bcs and curves based on loc_to_glob

      call fgslib_crystal_free(cr_re2)

#ifdef DPROCMAP
      call dProcmapInit()
      do i = 1,nelt
         ieg = lglel(i)
         if (ieg.lt.1 .or. ieg.gt.nelgt)
     $      call exitti('invalid ieg!$',ieg)
         ibuf(1) = i
         ibuf(2) = nid
         call dProcmapPut(ibuf,2,0,ieg)
      enddo
#else
      call izero(gllnid,nelgt)
      do i = 1,nelt
         ieg = lglel(i)
         gllnid(ieg) = nid
      enddo
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllnid(k),iwork,'+  ',m)
         k = k+m
      enddo
#endif

      if (nelt.gt.lelt) then
         call exitti('nelt > lelt, increase lelt!$',nelt)
      endif

      DO 1200 IFIELD=MFIELD,NFLDT
         IF (IFTMSH(IFIELD)) THEN
            NELG(IFIELD)      = NELGT
         ELSE
            NELG(IFIELD)      = NELGV
         ENDIF
 1200 CONTINUE

      dtmp = dnekclock_sync() - etime0
      if(nio.eq.0 .and. loglevel .gt. 1) then
        write(6,'(A,g13.5,A,/)')  ' done :: partioning ',dtmp,' sec'
      endif

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

      if(get_vert_called.gt.0) return

      nv = 2**ldim
      call get_vert_map(nv)

      get_vert_called = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_big(loc_to_glob_nid)
c
c     Distribute and assign partitions using the .map file
c
      include 'SIZE'
      include 'TOTAL'

      integer loc_to_glob_nid(lelt)

      if(get_vert_called.gt.0) return

      nv = 2**ldim
      call get_vert_map_big(nv, loc_to_glob_nid)

      get_vert_called = 1

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_map_big(nlv, loc_to_glob_nid)

      include 'SIZE'
      include 'TOTAL'

      parameter(mdw=2+2**ldim)
      parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
      common /scrns/ wk(mdw*ndw)
      integer*8 wk

      integer nlv, loc_to_glob_nid(lelt)
      integer     wk4(2*mdw*ndw)
      equivalence (wk4,wk)

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer*8 eid8(lelt), vtx8(lelt*2**ldim)
      integer   iwork(lelt), dest(lelt)
      common /ctmp0/ eid8, vtx8, iwork, dest

      common /scrcg/ xyz(ldim*lelt*2**ldim)

      integer cnt, algo
      integer opt_parrsb(3), opt_parmetis(10)

      logical ifbswap, ifread_con

      real tol

#if !defined(PARRSB) && !defined(PARMETIS)
#if defined(DPROCMAP)
      call exitti('DPROCMAP requires PARRSB or PARMETIS!$',0)
#else
      call read_map(vertex,nlv,wk4,mdw,ndw)
      return
#endif      
#endif

#if defined(PARRSB) || defined(PARMETIS)
      neli = nelt
      ifread_con = .true.
      call read_con(wk4,size(wk),neli,nlv,ierr)

      if (ierr.ne.0) then
        ifread_con = .false.
        tol = connectivityTol
        call find_con(wk,size(wk),tol,ierr)
        if(ierr.ne.0) then
          tol = tol / 10.0;
          call find_con(wk,size(wk),tol,ierr)
        endif
        call err_chk(ierr,' find_con failed!$')
      endif

c fluid elements
      j  = 0
      ii = 0
      cnt= 0
      do i = 1, neli
         itmp = wk(ii+1)
         if (ifread_con) itmp = wk4(ii+1)

         if (itmp .le. nelgv) then
            j = j + 1

            eid8(j) = wk(ii+1)
            call i8copy(vtx8((j-1)*nlv+1),wk(ii+2),nlv)
            if (ifread_con) then
              eid8(j) = wk4(ii+1)
              call icopy48(vtx8((j-1)*nlv+1),wk4(ii+2),nlv)
            endif

            do iv=1,nlv
              xyz(cnt+1)=xc(iv,i)
              xyz(cnt+2)=yc(iv,i)
              if(ldim.eq.3) then
                xyz(cnt+3)=zc(iv,i)
                cnt=cnt+3
              else
                cnt=cnt+2
              endif
            enddo
         endif
         ii = ii + (nlv+1)
      enddo
      neliv = j

      call fpartMeshV2(dest,vtx8,xyz,neliv,nlv,nekcomm,
     $  meshPartitioner,0,loglevel,ierr)
      call err_chk(ierr,'partMesh fluid failed!$')

      do i = 1,neliv
         lglel(i) = eid8(i)
      enddo
      call isort(lglel,iwork,neliv)

      do i = 1,neliv
         call i8copy(vertex(1,i),vtx8((iwork(i)-1)*nlv+1),nlv)
         loc_to_glob_nid(i) = dest(iwork(i))
      enddo
      nelv = neliv
      nelt = neliv

c solid elements
      cnt=0
      if (nelgt.ne.nelgv) then
         write(6,*) 'There are solid elems'
         j  = 0
         ii = 0
         do i = 1,neli
            itmp = wk(ii+1)
            if (ifread_con) itmp = wk4(ii+1)

            if (itmp .gt. nelgv) then
               j = j + 1

               eid8(j) = wk(ii+1)
               call i8copy(vtx8((j-1)*nlv+1),wk(ii+2),nlv)
               if (ifread_con) then
                 eid8(j) = wk4(ii+1)
                 call icopy48(vtx8((j-1)*nlv+1),wk4(ii+2),nlv)
               endif

               do iv=1,nlv
                 xyz(cnt+1)=xc(iv,i)
                 xyz(cnt+2)=yc(iv,i)
                 if(ldim.eq.3) then
                   xyz(cnt+3)=zc(iv,i)
                   cnt=cnt+3
                 else
                   cnt=cnt+2
                 endif
               enddo
            endif
            ii = ii + (nlv+1)
         enddo
         nelit = j

         call fpartMeshV2(dest,vtx8,xyz,nelit,nlv,nekcomm,
     $                  2,0,loglevel,ierr)
         call err_chk(ierr,'partMesh solid failed!$')

         do i = 1, nelit
            lglel(nelv + i) = eid8(i)
         enddo
         ! sort locally by global element id
         call isort(lglel(nelv+1),iwork,nelit)

         do i = 1, nelit
            call i8copy(vertex(1,nelv+i),vtx8((iwork(i)-1)*nlv+1),nlv)
            loc_to_glob_nid(nelv + i) = dest(iwork(i))
         enddo
         nelt = nelv + nelit
      endif
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_vert_map(nlv)

      include 'SIZE'
      include 'TOTAL'

      parameter(mdw=2+2**ldim)
      parameter(ndw=7*lx1*ly1*lz1*lelv/mdw)
      common /scrns/ wk(mdw*ndw)
      integer*8 wk

      integer     wk4(2*mdw*ndw)
      equivalence (wk4,wk)

      common /ivrtx/ vertex ((2**ldim),lelt)
      integer*8 vertex

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer ibuf(2)

      integer*8 eid8(lelt), vtx8(lelt*2**ldim)
      integer   iwork(lelt)
      common /ctmp0/ eid8, vtx8, iwork

      common /scrcg/ xyz(ldim*lelt*2**ldim)

      integer cnt, algo
      integer opt_parrsb(3), opt_parmetis(10)

      logical ifbswap, ifread_con

      real tol

#if !defined(PARRSB) && !defined(PARMETIS)
#if defined(DPROCMAP)
      call exitti('DPROCMAP requires PARRSB or PARMETIS!$',0)
#else
      call read_map(vertex,nlv,wk4,mdw,ndw)
      return
#endif      
#endif

#if defined(PARRSB) || defined(PARMETIS)
      neli = nelt
      ifread_con = .true.
      call read_con(wk4,size(wk),neli,nlv,ierr)

      if (ierr.ne.0) then
        ifread_con = .false.
        tol = connectivityTol
        call find_con(wk,size(wk),tol,ierr)
        if(ierr.ne.0) then
          tol = tol / 10.0;
          call find_con(wk,size(wk),tol,ierr)
        endif
        call err_chk(ierr,' find_con failed!$')
      endif

c fluid elements
      j  = 0
      ii = 0
      cnt= 0
      do i = 1,neli
         itmp = wk(ii+1)
         if (ifread_con) itmp = wk4(ii+1)

         if (itmp .le. nelgv) then
            j = j + 1

            eid8(j) = wk(ii+1)
            call i8copy(vtx8((j-1)*nlv+1),wk(ii+2),nlv)
            if (ifread_con) then
              eid8(j) = wk4(ii+1)
              call icopy48(vtx8((j-1)*nlv+1),wk4(ii+2),nlv)
            endif

            do iv=1,nlv
              xyz(cnt+1)=xc(iv,i)
              xyz(cnt+2)=yc(iv,i)
              if(ldim.eq.3) then
                xyz(cnt+3)=zc(iv,i)
                cnt=cnt+3
              else
                cnt=cnt+2
              endif
            enddo
         endif
         ii = ii + (nlv+1)
      enddo
      neliv = j

      nel = neliv
      call fpartMesh(eid8,vtx8,xyz,lelt,nel,nlv,nekcomm,
     $  meshPartitioner,0,loglevel,ierr)
      call err_chk(ierr,'partMesh fluid failed!$')

      nelv = nel
      nelt = nelv
      ierr = 0 
      if (nelv .gt. lelv) ierr = 1
      call err_chk(ierr,'nelv > lelv!$')
 
      do i = 1,nelv
         lglel(i) = eid8(i)
      enddo
      call isort(lglel,iwork,nelv)
      do i = 1,nelv
         call i8copy(vertex(1,i),vtx8((iwork(i)-1)*nlv+1),nlv)
      enddo

      cnt=0
c solid elements
      if (nelgt.ne.nelgv) then
         j  = 0
         ii = 0
         do i = 1,neli
            itmp = wk(ii+1)
            if (ifread_con) itmp = wk4(ii+1)

            if (itmp .gt. nelgv) then
               j = j + 1

               eid8(j) = wk(ii+1)
               call i8copy(vtx8((j-1)*nlv+1),wk(ii+2),nlv)
               if (ifread_con) then
                 eid8(j) = wk4(ii+1)
                 call icopy48(vtx8((j-1)*nlv+1),wk4(ii+2),nlv)
               endif

               do iv=1,nlv
                 xyz(cnt+1)=xc(iv,i)
                 xyz(cnt+2)=yc(iv,i)
                 if(ldim.eq.3) then
                   xyz(cnt+3)=zc(iv,i)
                   cnt=cnt+3
                 else
                   cnt=cnt+2
                 endif
               enddo
            endif
            ii = ii + (nlv+1)
         enddo
         nelit = j

         nel = nelit
         call fpartMesh(eid8,vtx8,xyz,lelt,nel,nlv,nekcomm,
     $                  2,0,loglevel,ierr)
         call err_chk(ierr,'partMesh solid failed!$')

         nelt = nelv + nel
         ierr = 0 
         if (nelt .gt. lelt) ierr = 1
         call err_chk(ierr,'nelt > lelt!$')
    
         do i = 1,nel
            lglel(nelv+i) = eid8(i)
         enddo
         call isort(lglel(nelv+1),iwork,nel) ! sort locally by global element id
         do i = 1,nel
            call i8copy(vertex(1,nelv+i),vtx8((iwork(i)-1)*nlv+1),nlv)
         enddo
      endif

#ifdef DPROCMAP
      call dProcMapClearCache()
      do i = 1,nelt
         ieg = lglel(i)
         if (ieg.lt.1 .or. ieg.gt.nelgt) 
     $      call exitti('invalid ieg!$',ieg)
         ibuf(1) = i
         ibuf(2) = nid
         call dProcmapPut(ibuf,2,0,ieg)
      enddo
#else
      call izero(gllnid,nelgt)
      do i = 1,nelt
         ieg = lglel(i)
         gllnid(ieg) = nid
      enddo
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllnid(k),iwork,'+  ',m)
         k = k+m
      enddo
#endif

#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine read_con(wk,nwk,nelr,nv,ierr)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer nwk,nelr,nv,ierr
      integer wk(nwk)
     

      logical ifbswap,if_byte_swap_test
      logical ifco2, ifcon

      character*132 confle
      character*1   confle1(132)
      equivalence  (confle,confle1)

      character*132 hdr
      character*5   version
      real*4        test

      integer nvi
      integer*8 nelgti,nelgvi
      integer*8 offs, offs0

      ierr = 0

      ifco2 = .false.
      ifmpiio = .true.
#ifdef NOMPIIO
      ifmpiio = .false.
#endif

      if (nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (confle,132)
         call chcopy(confle,reafle,lfname)
         call chcopy(confle1(lfname+1),'.con',4)
         inquire(file=confle, exist=ifcon)

         if (.not.ifcon) then
            call chcopy(confle1(lfname+1),'.co2',4)
            inquire(file=confle, exist=ifco2)
         endif

        if(.not.ifcon .and. .not.ifco2) ierr = 1
      endif

      call bcast(ierr,sizeof(ierr))
      if(ierr.ne.0) return

      call bcast(confle,sizeof(confle))
      if(nid.eq.0) write(6,'(A,A)') ' reading ', confle
      call bcast(ifco2,lsize)

      ierr = 0

      ! read header
      if (nid.eq.0) then
         if (ifco2) then
            call byte_open(confle,ierr)
            if(ierr.ne.0) goto 100

            call blank(hdr,sizeof(hdr))
            call byte_read(hdr,sizeof(hdr)/4,ierr)
            if(ierr.ne.0) goto 100

            read (hdr,'(a5)') version 

            if (version.eq.'#v002') then
               read (hdr,*) version,nelgti,nelgvi,nvi
            else
               read (hdr,1) version,nelgti,nelgvi,nvi
            endif
            write (6,*) ' co2 hdr:',version,nelgti,nelgvi,nvi 

            call byte_read(test,1,ierr)
            if(ierr.ne.0) goto 100
            ifbswap = if_byte_swap_test(test,ierr)
            if(ierr.ne.0) goto 100
         endif
      endif

   1  format(a5,3i12)

      call bcast(nelgti,sizeof(nelgti))
      call bcast(nelgvi,sizeof(nelgvi))
      call bcast(nvi,sizeof(nvi))
      call bcast(ifbswap,sizeof(ifbswap))

      if (nvi .ne. nv)
     $   call exitti('Number of vertices do not match!$',0)
      if (nelgti .ne. nelgt)
     $   call exitti('nelgt for mesh/con differs!$',0)
      if (nelgvi .ne. nelgv)
     $   call exitti('nelgt for mesh/con differs!$',0)

      if (ifco2 .and. ifmpiio) then
        if (nid.eq.0) call byte_close(ierr)
        call byte_open_mpi(confle,ifh,.true.,ierr)
        offs0 = sizeof(hdr) + sizeof(test)

        call lim_chk(nelr*(nvi+1),nwk,'nelr ','nwk   ','read_con  ')

        nelBr = igl_running_sum(nelr) - nelr
        offs  = offs0 + int(nelBr,8)*(nvi+1)*ISIZE

        call byte_set_view(offs,ifh)
        call byte_read_mpi(wk,(nvi+1)*nelr,-1,ifh,ierr)
        call err_chk(ierr,' Error while reading con file!$')
        call byte_close_mpi(ifh,ierr)
        if (ifbswap) call byte_reverse(wk,(nvi+1)*nelr,ierr)
      endif

      return


 100  continue
      call err_chk(ierr,'Error opening/reading con file$')
      return

      end
c-----------------------------------------------------------------------
#if defined(PARRSB)      
      subroutine find_con(wk,nwk,tol,ierr)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer nwk,ierr
      integer*8 wk(nwk)
      real tol

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      common /scrcg/ xyz(ldim*(2**ldim)*lelt)

      integer*8 eid8(4*lelt),vtx8(lelt*(2**ldim+1))
      common /ctmp0/ eid8, vtx8, iwork

      ierr = 0

      nv=2**ndim
      nf=2*ndim

      k=0
      if(ndim.eq.3) then
        do i=1,nelt
          do j=1,nv
            xyz(k+1)=xc(j,i)
            xyz(k+2)=yc(j,i)
            xyz(k+3)=zc(j,i)
            k=k+3
          enddo
        enddo
      else
        do i=1,nelt
          do j=1,nv
            xyz(k+1)=xc(j,i)
            xyz(k+2)=yc(j,i)
            k=k+2
          enddo
        enddo
      endif

      istart=igl_running_sum(nelt)-nelt

      !calculate number of periodic pairs
      npf=0
      do i=1,nelt
        do j=1,nf
          if(cbc(j,i,1).eq.'P  ') then
            eid8(4*npf+1)=i+istart
            eid8(4*npf+2)=j
            eid8(4*npf+3)=bc(1,j,i,1)
            eid8(4*npf+4)=bc(2,j,i,1)
            npf=npf+1
          endif
        enddo
      enddo

      call fparrsb_find_conn(vtx8,xyz,nelt,ndim,eid8,npf,tol,nekcomm,
     $  1,ierr)

      k=1
      l=1
      do i=1,nelt
        wk(k)=istart+i
        k=k+1
        do j=1,nv 
          wk(k)=vtx8(l)
          k=k+1
          l=l+1
        enddo
      enddo

      end
#endif
c-----------------------------------------------------------------------
      subroutine set_proc_map()
C
C     Compute element to processor distribution according to (weighted) 
C     physical distribution in an attempt to minimize exposed number of
C     element interfaces.
C
      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'DPROCMAP'

      integer ibuf(2)
      logical ifbswap
      logical ifre2

      integer iwork(lelt)

      if(nid.eq.0) inquire(file=re2fle, exist=ifre2)
      call bcast(ifre2,lsize)
      if(.not.ifre2) then
        call set_proc_map_legacy
        return
      endif

      call read_re2_hdr(ifbswap, .false.)
      nelt = nelgt/np
      do i = 1,mod(nelgt,np)
        if (np-i.eq.nid) nelt = nelt + 1
      enddo

      if (nelt .gt. lelt)
     $   call exitti('nelt > lelt!$',nelt)

      ! setup gllnid + gllel
#if defined(DPROCMAP)
      call dProcmapInit()  
#endif
      nelB = igl_running_sum(nelt) - nelt
      do i = 1,nelt
         ieg = nelB + i
         lglel(i) = ieg 
         if (ieg.lt.1 .or. ieg.gt.nelgt)
     $      call exitti('invalid ieg!$',ieg)
#if defined(DPROCMAP)
         ibuf(1) = i
         ibuf(2) = nid
         call dProcmapPut(ibuf,2,0,ieg)
#else
         gllnid(ieg) = nid
         gllel(ieg) = i 
#endif
      enddo
#if !defined(DPROCMAP)
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllnid(k),iwork,'+  ',m)
         if (m.gt.0) call igop(gllel(k) ,iwork,'+  ',m)
         k = k+m
      enddo
#endif

      ! read coord for RCB and/or connectivity
      call read_re2_data(ifbswap, .true., .false., .true.) 

      ! get element-proc mapping      
      call get_map() 

#if !defined(DPROCMAP)
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
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllel(k),iwork,'+  ',m)
         k = k+m
      enddo

      do ieg=1,nelgt
         mid  =gllnid(ieg)
         ie   =gllel (ieg)
         if (mid.eq.nid) lglel(ie)=ieg
      enddo
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine set_proc_map_legacy()
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
      common /ctmp0/ iwork(lelt)

      REAL*8 dnekclock,t0

#if defined(PARRSB) || defined(PARMETIS) || defined(DPROCMAP)
      call exitti(' DPROCMAP/PARRSB not supported for rea files$',0)
#else      
      t0 = dnekclock()
      call get_map

c     compute global to local map (no processor info)
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

c     dist. global to local map to all processors
      npass = 1 + nelgt/lelt
      k=1
      do ipass = 1,npass
         m = nelgt - k + 1
         m = min(m,lelt)
         if (m.gt.0) call igop(gllel(k),iwork,'+  ',m)
         k = k+m
      enddo

c     compute local to global map
c     (i.e. returns global element number given local index and proc id)
      do ieg=1,nelgt
         mid  =gllnid(ieg)
         ie   =gllel (ieg)
         if (mid.eq.nid) lglel(ie)=ieg
      enddo
#endif

      return
      end
c-----------------------------------------------------------------------

#ifndef DPROCMAP

c-----------------------------------------------------------------------
      subroutine read_map(vertex,nlv,wk,mdw,ndw)

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'

      integer*8 vertex(nlv,1)
      integer wk(mdw,ndw)

      logical ifbswap,if_byte_swap_test

      character*132 mapfle
      character*1   mapfle1(132)
      equivalence  (mapfle,mapfle1)

      character*132 hdr
      character*5   version
      real*4        test

      logical ifma2,ifmap
      integer e,eg,eg0,eg1
      integer itmp20(20)

      ierr = 0
      ifma2 = .false.

      if (nid.eq.0) then
         lfname = ltrunc(reafle,132) - 4
         call blank (mapfle,132)
         call chcopy(mapfle,reafle,lfname)
         call chcopy(mapfle1(lfname+1),'.map',4)
         inquire(file=mapfle, exist=ifmap)

         if (.not.ifmap) then
            call chcopy(mapfle1(lfname+1),'.ma2',4)
            inquire(file=mapfle, exist=ifma2)
         endif

        if(.not.ifmap .and. .not.ifma2) ierr = 1 
      endif
      if(nid.eq.0) write(6,'(A,A)') ' Reading ', mapfle
      call err_chk(ierr,' Cannot find map file!$')
      call bcast(ifma2,lsize)
      ierr = 0

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
               gllnid(eg) = wk(1,m)  ! must still be divided
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

      lng = isize*neli
      call bcast(gllnid,lng)
      call assign_gllnid(gllnid,gllel,nelgt,nelgv,np) ! gllel is used as scratch

      nelt=0 !     Count number of elements on this processor
      nelv=0
      do eg=1,neli
         if (gllnid(eg).eq.nid) then
            if (eg.le.nelgv) nelv=nelv+1
            if (eg.le.nelgt) nelt=nelt+1
         endif
      enddo

c      if (np.le.64) write(6,*) nid,nelv,nelt,nelgv,nelgt,' NELV'

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

      key = mdw  ! Sort tuple list by eg
      nkey = 1
      call fgslib_crystal_ituple_sort(cr_h,wk,mdw,nelt,key,nkey)

      iflag = 0
      if (ntuple.ne.nelt) then
         write(6,*) nid,ntuple,nelv,nelt,nelgt,' NELT FAIL'
         write(6,*) 'Check that .map file and .rea file agree'
         iflag=1
      else
         do e=1,nelt
            call icopy48(vertex(1,e),wk(2,e),nlv)
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
c-----------------------------------------------------------------------

#endif
