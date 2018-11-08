      program exo2nek

      use SIZE

      call read_input_name
      call exodus_read
      call convert
      call gen_re2

      end 
!-----------------------------------------------------------------------
      subroutine read_input_name

      use SIZE

      character(1)  re2nam1(80)
      character(1)  exonam1(32)
      character(32) fname

      write(6,*) 'Input (.exo) file name:'
      read (5,'(A32)') fname
      len = ltrunc(fname,32)
      
      call blank  (exonam1, 32)
      call blank  (re2nam1, 80)
      call chcopy (exonam1,fname,32)
      call chcopy (re2nam1,fname,80)
      call chcopy (exonam1(len+1) ,'.exo',4)
      call chcopy (re2nam1(len+1) ,'.re2',4)

      call blank  (exoname, 32)
      call blank  (re2name, 80)
      call chcopy (exoname,exonam1,len+4)
      call chcopy (re2name,re2nam1,len+4)

      return 
      end
!-----------------------------------------------------------------------
      subroutine exodus_read
!
!  Subroutine to read an exodusII binary file containing a mesh.
!  It uses exodus fortran binding subroutines, which depend on
!  the netcdf library for low level data access.
!
      use SIZE
      include 'exodusII.inc'

      integer exoid, cpu_ws, io_ws

      character(MXSTLN) typ, qa_record(4,10)
      character(MXLNLN) titl
      character(1)      cdum

      integer,allocatable,dimension(:) :: idblk
      integer,allocatable,dimension(:) :: num_nodes_per_elem
      integer,allocatable,dimension(:) :: num_attr    !not used
!
! open EXODUS II file
!
      cpu_ws = 8 ! use real*8 to communicate with exodus
      io_ws  = 0
      exoid  = exopen (exoname, EXREAD, cpu_ws, io_ws, vers, ierr)
      if (ierr.lt.0) then
        write(6,'(2a)') "ERROR: cannot open file ", exoname 
        STOP
      endif
      write(6,*)
      write(6,'(a32,a,f4.2)') & 
            exoname," is an EXODUSII file; version ",vers
      write(6,'(a,i2)') "I/O word size", io_ws
!
! read database parameters
!
      call exgini (exoid, titl, num_dim, num_nodes, num_elem, &
                   num_elem_blk, num_node_sets, num_side_sets, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read exodusII parameters (exgini)"
        STOP
      endif
      write (6, '(/"database parameters:"/ /        &
                 "title         = ", a81 / /        &
                 "num_dim       = ", i8 /           &
                 "num_nodes     = ", i8 /           &
                 "num_elem      = ", i8 /           &
                 "num_elem_blk  = ", i8 /           &
                 "num_side_sets = ", i8)')          &
                 titl,num_dim, num_nodes, num_elem, &
                  num_elem_blk, num_side_sets
      write (6,*)
!
! allocate some arrays    
!
      ! EXODUS:
      allocate ( idblk              (num_elem_blk)        )
      allocate ( num_nodes_per_elem (num_elem_blk)        )
      allocate ( num_attr           (num_elem_blk)        )
      allocate ( num_elem_in_block  (num_elem_blk)        )
      allocate ( num_sides_in_set   (num_side_sets)       )
      allocate ( idss               (num_side_sets)       )
      allocate ( connect            (3**num_dim*num_elem) )
      allocate ( x_exo              (3**num_dim*num_elem) )
      allocate ( y_exo              (3**num_dim*num_elem) )
      allocate ( z_exo              (3**num_dim*num_elem) )
      ! Nek5000:
      allocate ( xm1                (3,3,3,num_elem)      )
      allocate ( ym1                (3,3,3,num_elem)      )
      allocate ( zm1                (3,3,3,num_elem)      )
!
! read element block parameters
!
      call exgebi (exoid, idblk, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read block ids (exgebi)"
        STOP
      endif

      do i = 1, num_elem_blk
        call exgelb (exoid, idblk(i), typ, num_elem_in_block(i), &
                     num_nodes_per_elem(i), num_attr(i), ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3,a)') &
          "ERROR: cannot read parameters for block ",i," (exgelb)"
          STOP
        endif
        write (6, '("element block id   = ", i8,/       &
                   "element type       = ", 3x,a8,/     &
                   "num_elem_in_block  = ", i8,/        &
                   "num_nodes_per_elem = ", i8)')       &
                   idblk(i), typ, num_elem_in_block(i), &
                   num_nodes_per_elem(i)
        write(6,*)

        if (i.eq.1) then
          nvert=num_nodes_per_elem(i)
          if (num_dim.eq.2) then
            if (nvert.ne.8) then
              if (nvert.eq.9) then
                write(6,*)
                write(6,'(a)') &
                "WARNING: QUAD9 elements are not officially supported"
                write(6,'(a)') &
                "as there is no exodus standard for this element type."
                write(6,*)
              else
                write(6,*)
                write(6,'(a)') &
                "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
                STOP
              endif
            endif      
          elseif (num_dim.eq.3) then
            if (nvert.ne.20) then
              if (nvert.eq.27) then
                write(6,*)
                write(6,'(a)') &
                "WARNING: HEX27 elements are not officially supported"
                write(6,'(a)') &
                "as there is no exodus standard for this element type."
                write(6,*)
              else
                write(6,*)
                write(6,'(a)') &
                "ERROR: Only HEX20 elements are allowed in a 3D mesh!"
                STOP
              endif
            endif      
          else
          write(6,'(a,i3)') &
           "ERROR: Unknown number of dimensions! num_dim= ",num_dim
          STOP
          endif
        endif

        if (num_nodes_per_elem(i).ne.nvert) then
          write(6,*)
          write(6,'(a)') &
           "ERROR: All blocks should contain elements of the same type!"
          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
          write(6,'(a,i3)') "num_nodes_per_elem of block 1 ",nvert
          STOP
        endif
      enddo
!
! read nodal coordinates values from database
!
      call exgcor (exoid, x_exo, y_exo, z_exo, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read nodal coordinates (exgcor)"
        STOP
      endif
!
! read element connectivity
!
      iend = 0
      do 60 i = 1, num_elem_blk
        istart = iend + 1
        call exgelc (exoid, idblk(i), connect(istart), ierr)
        iend = iend+num_nodes_per_elem(i)*num_elem_in_block(i)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read elm. connectivity (exgelc)"
          STOP
        endif
60    continue
!
! read individual side sets
!
      num_sides_tot = 0
      if (num_side_sets .gt. 0) then
        call exgssi (exoid, idss, ierr)
        if (ierr.lt.0) then
          write(6,'(a)') "ERROR: cannot read SideSet ids (exgssi)"
          STOP
        endif
        maxnss=0
        do i = 1, num_side_sets
          call exgsp (exoid,idss(i),num_sides_in_set(i),idum,ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgsp)"
            STOP
          endif
          write (6, '("side set ", i2, " num_sides = ", i8)') &
                 idss(i), num_sides_in_set(i)
          num_sides_tot = num_sides_tot + num_sides_in_set(i)
          maxnss        = max(maxnss,num_sides_in_set(i))
        enddo

        ! allocate sideset arrays
        allocate (elem_list(maxnss,num_side_sets) )
        allocate (side_list(maxnss,num_side_sets) )

        do i = 1, num_side_sets
          call exgss (exoid,idss(i),elem_list(1,i),side_list(1,i),ierr)
          if (ierr.lt.0) then
            write(6,'(a,i3,a)') &
            "ERROR: cannot read parameters for SideSet No.",i," (exgss)"
            STOP
          endif
        enddo
        write (6,*)
      else
        write(6,'(a)') "WARNING: No SideSets in exodus file!"
      endif
!
! read QA records
!
      call exinq (exoid, EXQA, num_qa_rec, fdum, cdum, ierr)
      if (ierr.lt.0) then
        write(6,'(a)') "ERROR: cannot read QA records (exinq QA) "
        STOP
      elseif (ierr.gt.0) then
        write(6,'(a)') "INFO: file does not contain any QA records"
      else
        call exgqa (exoid, qa_record, ierr)
        if (ierr.lt.0) then
          write(6,'(a,i3)') "WARNING: cannot read QA records (exgqa)"
        else
          write (6, '("QA records = ")')
          if (num_qa_rec.gt.10) then
            write(6,'(2a)') &
              'WARNING: Cannot handle more than 10 QA records', &
              'Printing only the first 10...'
          else
            do i = 1, num_qa_rec
              do j = 1, 4
                write (6,'(a)') qa_record(j,i)
              enddo
            enddo
          endif
        endif
      endif

      return
      end
!-----------------------------------------------------------------
      subroutine convert
!
!  Subroutine to convert an already read exodusII mesh to a nek
!  mesh. The idea is to fill each element's node coordinates of
!  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
!  coordinates.
!
      use SIZE
      include 'exodusII.inc'

! node and face conversion (it works at least for cubit):
      integer exo_to_nek_vert3D(27)                   ! hex27 to nek numbering
      data    exo_to_nek_vert3D                     &
            / 19,  1,  7, 25, 21,  3,  9, 27, 10    &     
            ,  4, 16, 22, 20,  2,  8, 26, 12,  6    &
            , 18, 24, 14, 13, 15, 23,  5, 11, 17 /

      integer exo_to_nek_vert2D(9)                     ! quad9 to nek numbering
      data    exo_to_nek_vert2D   / 1, 3, 9, 7, 2, 6, 8, 4, 5  / 

      integer exo_to_nek_face3D(6)
      data    exo_to_nek_face3D  / 1, 5, 3, 6, 4, 2 /  ! symmetric face numbering

      integer exo_to_nek_face2D(4)
      data    exo_to_nek_face2D  / 1, 2, 3, 4 /        ! symmetric face numbering

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do iel=1,num_elem
        do ivert =1,nvert
          if (num_dim.eq.2) then
            jvert = exo_to_nek_vert2D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
          else
            jvert = exo_to_nek_vert3D(ivert)
            xm1(jvert,1,1,iel)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,iel)=y_exo(connect(nvert*(iel-1)+ivert))
            zm1(jvert,1,1,iel)=z_exo(connect(nvert*(iel-1)+ivert))
          endif
        enddo
      enddo
      deallocate(x_exo,y_exo,z_exo,connect)
      write(6,'(A)') 'done :: Converting elements '
!
! allocate and zero-out curve sides arrays
!
      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
      allocate   (curve  (2*num_dim,12,   num_elem) )
      call rzero (curve,2*num_dim*12*num_elem)
      call blank (ccurve,(4+8*(num_dim-2))*num_elem)
!
! allocate and zero-out bc arrays only if sidesets are specified
!
      if (num_side_sets.eq.0) return   

      allocate   (cbc    (2*num_dim,      num_elem) )
      allocate   (bc     (5,2*num_dim,    num_elem) ) 
      call rzero (bc,5*2*num_dim*num_elem)
      call blank (cbc,3*2*num_dim*num_elem)
!
! set bc's
!
      write(6,'(a)') ''
      write(6,'(a)') 'Converting SideSets ...'
      do iss=1,num_side_sets
        write(6,'(a)') ''
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
        do i=1,num_sides_in_set(iss) 
          iel = elem_list(i,iss)
          ifc = side_list(i,iss)
          if (num_dim.eq.2) then 
            jfc = exo_to_nek_face2D(ifc)
          else
            jfc = exo_to_nek_face3D(ifc)
          endif
          cbc(jfc,iel)   = 'EXO' ! dummy exodus bc 
          bc (5,jfc,iel) = idss(iss)
        enddo
      write(6,'(A,I2)') 'done :: Sideset ',idss(iss)
      enddo
      deallocate(elem_list,side_list)

      write(6,'(a)') ''
      write(6,'(a)') 'done :: Converting SideSets '

      return
      end
!--------------------------------------------------------------------
      subroutine gen_re2

      use SIZE

      write(6,*)
      write(6,'(A,A)') 'writing ', re2name

      call open_re2
      call write_xyz
      call write_curve
      call write_bc
      call close_re2

      return
      end
!--------------------------------------------------------------------
      subroutine open_re2

      use SIZE

      character(80) hdr


      real*4 test
      data   test  / 6.54321 /

      call byte_open(re2name,ierr)
            
!  Write the header
      call blank   (hdr,80)    
      write(hdr,1) num_elem, num_dim, num_elem
    1 format('#v003',i9,i3,i9,' this is the hdr')
      call byte_write(hdr,20,ierr)         
      call byte_write(test,1,ierr)     ! write the endian discriminator

      return
      end
!--------------------------------------------------------------------
      subroutine write_xyz

      use SIZE

      real     xx(8), yy(8), zz(8)
      real*8   rgroup, buf2(30)

      integer isym2pre(8)   ! Symmetric-to-prenek vertex ordering
      data    isym2pre    / 1 , 2 , 4 , 3 , 5 , 6 , 8 , 7 /

      nxs = 3-1  ! always nx1=ny1=nz1=3 here
      nys = 3-1
      nzs = 3-1

      igr    = 0
      rgroup = igr

      do iel=1,num_elem
        l = 0
        if (num_dim.eq.3) then
          do k=0,1
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1+k*nzs,iel)
            zz(li) = zm1(1+i*nxs,1+j*nys,1+k*nzs,iel)
          enddo
          enddo
          enddo
        else
          do j=0,1
          do i=0,1
            l=l+1
            li=isym2pre(l)
            xx(li) = xm1(1+i*nxs,1+j*nys,1,iel)
            yy(li) = ym1(1+i*nxs,1+j*nys,1,iel)
          enddo
          enddo
        endif

        call byte_write(rgroup, 2,ierr)

        if (num_dim.eq.3) then
          call copy        (buf2(1) ,xx,8)
          call copy        (buf2(9) ,yy,8)
          call copy        (buf2(17),zz,8)
          call byte_write  (buf2(1) ,16, ierr)
          call byte_write  (buf2(9) ,16, ierr)
          call byte_write  (buf2(17),16, ierr)
        else
          call copy        (buf2(1),xx,4)
          call copy        (buf2(5),yy,4)
          call byte_write  (buf2(1),8, ierr)
          call byte_write  (buf2(5),8, ierr)
        endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine write_curve

      use SIZE

      real*8     buf2(30)
      real*8     rcurve

      character(1) cc

      do iel=1,num_elem
         call gen_rea_midside_e(iel)
      enddo

      nedge = 4 + 8*(num_dim-2)
      ncurv = 0
      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') ncurv = ncurv+1
        enddo
      enddo

      rcurve = ncurv
      call byte_write(rcurve,2, ierr)

      do iel=1,num_elem
        do iedge=1,nedge
          if (ccurve(iedge,iel).ne.' ') then
            if (ccurve(iedge,iel).eq.'C') cc='C'
            if (ccurve(iedge,iel).eq.'s') cc='s'
            if (ccurve(iedge,iel).eq.'m') cc='m'
            buf2(1) = iel
            buf2(2) = iedge
            call copy       (buf2(3),curve(1,iedge,iel),5)
            call blank      (buf2(8),8)
            call chcopy     (buf2(8),cc,1)
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine write_bc
      
      use SIZE

      real*8  rbc, buf2(30)

      character(3) ch3
      character(1) chdum
      data         chdum /' '/

      if (num_side_sets.eq.0) return

      rbc = num_sides_tot
      call byte_write (rbc,2,ierr)

      do iel = 1,num_elem
        do ifc = 1,2*num_dim
          ch3 = cbc(ifc,iel)
          if (ch3.eq.'EXO') then
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (num_elem.ge.1000000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

      return
      end 
!-----------------------------------------------------------------------
      subroutine close_re2

      call byte_close (ierr)

      return
      end
!-----------------------------------------------------------------------
      subroutine gen_rea_midside_e(e)

      use SIZE

      real         len
      real         x3(27),y3(27),z3(27),xyz(3,3)
      character(1) ccrve(12)
      integer      e,edge

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1, &
                   19,20,21,   21,24,27,   27,26,25,   25,22,19, &
                    1,10,19,    3,12,21,    9,18,27,    7,16,25  /

      call chcopy(ccrve,ccurve(1,e),12)

      call map2reg                   (x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
      call map2reg                   (y3,3,ym1(1,1,1,e),1)
      if (num_dim.eq.3) call map2reg (z3,3,zm1(1,1,1,e),1)

!     Take care of spherical curved face defn
      if (ccurve(5,e).eq.'s') then
         call chcopy (ccrve(1),'ssss',4) ! face 5
         call chcopy (ccrve(5),' ',1)    ! face 5
      endif
      if (ccurve(6,e).eq.'s') then
         call chcopy (ccrve(5),'ssss',4) ! face 6
      endif

      tol   = 1.e-4
      tol2  = tol**2
      nedge = 4 + 8*(num_dim-2)

      do i=1,nedge
         if (ccrve(i).eq.' ') then
            do j=1,3
               xyz(1,j) = x3(e3(j,i))
               xyz(2,j) = y3(e3(j,i))
               xyz(3,j) = z3(e3(j,i))
            enddo
            len = 0.
            h   = 0.
            do j=1,num_dim
               xmid = .5*(xyz(j,1)+xyz(j,3))
               h    = h   + (xyz(j,2)-xmid)**2
               len  = len + (xyz(j,3)-xyz(j,1))**2
            enddo
            if (h.gt.tol2*len) ccurve(i,e) = 'm'
            if (h.gt.tol2*len) call copy (curve(1,i,e),xyz(1,2),num_dim)
         endif
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg(ur,n,u,nel)
!
!     Map scalar field u() to regular n x n x n array ur

      use SIZE

      real    ur(1), u(3*3*3,1)
      integer e

      ldr = n**num_dim

      k=1
      do e=1,nel
         if (num_dim.eq.2) call map2reg_2di_e(ur(k),n,u(1,e),3)
         if (num_dim.eq.3) call map2reg_3di_e(ur(k),n,u(1,e),3)
         k = k + ldr
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg_2di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n),uc(m,m)

      parameter (l=50)
      common /cmap2d/ j(l*l),jt(l*l),w(l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      call mxmf2(j,n,uc,m,w ,m)
      call mxmf2(w,n,jt,m,uf,n)

      return
      end
!-----------------------------------------------------------------------
      subroutine map2reg_3di_e(uf,n,uc,m) ! Fine, uniform pt

      real uf(n,n,n),uc(m,m,m)

      parameter (l=50)
      common /cmap3d/ j(l*l),jt(l*l),v(l*l*l),w(l*l*l),z(l)

      integer mo,no
      save    mo,no
      data    mo,no / 0,0 /

      if (m.ne.mo .or. n.ne.no ) then

          call zwgll (z,w,m)
          call zuni  (w,n)

          call gen_int_gz(j,jt,w,n,z,m)

      endif

      mm = m*m
      mn = m*n
      nn = n*n

      call mxmf2(j,n,uc,m,v ,mm)
      iv=1
      iw=1
      do k=1,m
         call mxmf2(v(iv),n,jt,m,w(iw),n)
         iv = iv+mn
         iw = iw+nn
      enddo
      call mxmf2(w,nn,jt,m,uf,n)

      return
      end
!-----------------------------------------------------------------------
      subroutine gen_int_gz(j,jt,g,n,z,m)

!     Generate interpolater from m z points to n g points

!        j   = interpolation matrix, mapping from z to g
!        jt  = transpose of interpolation matrix
!        m   = number of points on z grid
!        n   = number of points on g grid

      real j(n,m),jt(m,n),g(n),z(m)

      mpoly  = m-1
      do i=1,n
         call fd_weights_full(g(i),z,mpoly,0,jt(1,i))
      enddo

      call transpose(j,n,jt,m)

      return
      end
!-----------------------------------------------------------------------
      SUBROUTINE BLANK(A,N)
      CHARACTER(1) A(1)
      CHARACTER(1) BLNK
      SAVE        BLNK
      DATA        BLNK /' '/
 
      DO 10 I=1,N
         A(I)=BLNK
   10 CONTINUE
      RETURN
      END
!-----------------------------------------------------------------------
      subroutine exitt

      stop
      return
      end
!-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      real a(1),b(1)

      do i=1,n
         a(i)=b(i)
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine chcopy(a,b,n)
      CHARACTER(1) A(1), B(1)
 
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
!-----------------------------------------------------------------------
      subroutine icopy(a,b,n)
      INTEGER A(1), B(1)
 
      DO 100 I = 1, N
 100     A(I) = B(I)
      return
      END
!-----------------------------------------------------------------------
      subroutine fd_weights_full(xx,x,n,m,c)
!
!     This routine evaluates the derivative based on all points
!     in the stencils.  It is more memory efficient than "fd_weights"
!
!     This set of routines comes from the appendix of
!     A Practical Guide to Pseudospectral Methods, B. Fornberg
!     Cambridge Univ. Press, 1996.   (pff)
!
!     Input parameters:
!       xx -- point at wich the approximations are to be accurate
!       x  -- array of x-ordinates:   x(0:n)
!       n  -- polynomial degree of interpolant (# of points := n+1)
!       m  -- highest order of derivative to be approxxmated at xi
!
!     Output:
!       c  -- set of coefficients c(0:n,0:m).
!             c(j,k) is to be applied at x(j) when
!             the kth derivative is approxxmated by a
!             stencil extending over x(0),x(1),...x(n).
!
!
      real x(0:n),c(0:n,0:m)
 
      c1       = 1.
      c4       = x(0) - xx
 
      do k=0,m
      do j=0,n
         c(j,k) = 0.
      enddo
      enddo
      c(0,0) = 1.
 
      do i=1,n
         mn = min(i,m)
         c2 = 1.
         c5 = c4
         c4 = x(i)-xx
         do j=0,i-1
            c3 = x(i)-x(j)
            c2 = c2*c3
            do k=mn,1,-1
               c(i,k) = c1*(k*c(i-1,k-1)-c5*c(i-1,k))/c2
            enddo
            c(i,0) = -c1*c5*c(i-1,0)/c2
            do k=mn,1,-1
               c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
            enddo
            c(j,0) = c4*c(j,0)/c3
         enddo
         c1 = c2
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine transpose(a,lda,b,ldb)
      real a(lda,1),b(ldb,1)
 
      do j=1,ldb
         do i=1,lda
            a(i,j) = b(j,i)
         enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      function ltrunc(string,l)
      CHARACTER(1) STRING(L)
      CHARACTER(1)   BLNK
      DATA BLNK/' '/
 
      DO 100 I=L,1,-1
         L1=I
         IF (STRING(I).NE.BLNK) GOTO 200
  100 CONTINUE
      L1=0
  200 CONTINUE
      LTRUNC=L1
      return
      END
!-----------------------------------------------------------------------
      subroutine zuni(z,np)
!
!     Generate equaly spaced np points on the interval [-1:1]
!
      real z(1)

      dz = 2./(np-1)
      z(1) = -1.
      do i = 2,np-1
         z(i) = z(i-1) + dz
      enddo
      z(np) = 1.

      return
      end
!-----------------------------------------------------------------------
      subroutine rzero(a,n)
      real A(1)
      DO 100 I = 1, N
 100     A(I ) = 0.0
      return
      END

