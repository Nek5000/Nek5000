! Haomin Yuan, 6/15/2020
! include all features:
! 1, CHT mesh merge
! 2, mutiple fluid and solid exo files merge
! 3, auto correct non-right hand hex elements
! 4, setting periodicity
! 5, tet-to-hex, wedge-to-hex conversion
!---------------------------------------------------------------------------------
      program exo2nek

      use SIZE

      integer option
      character ifmore,if_stitch_interface
      logical ifmorefluid
      integer iexo1
      !logical ifmoresolid

      write(6,*) 'Special exo2nek version to merge multiple fluid/solid exo files'
      write(6,*) 'notes: for 3D mesh only'
      write(6,*) 'notes: interface surface mesh should be conformal and match spatially'

      write(6,*) 'please input estimate final total hex element number:'
      write(6,*) 'if this number is smaller than actual converted total hex element number,'
      write(6,*) 'you will encounter invalid memory error '

      read (5,*) etot_est

      ! allocate Nek5000 array
      num_dim = 3
      allocate ( xm1                (3,3,3,etot_est)      )
      allocate ( ym1                (3,3,3,etot_est)      )
      allocate ( zm1                (3,3,3,etot_est)      )

      allocate   (ccurve (4+8*(num_dim-2),etot_est) )
      allocate   (curve  (2*num_dim,12,   etot_est) )
      call rzero (curve,2*num_dim*12*etot_est)
      call blank (ccurve,(4+8*(num_dim-2))*etot_est)

      allocate   (cbc    (2*num_dim,      etot_est) )
      allocate   (bc     (5,2*num_dim,    etot_est) )
      call rzero (bc,5*2*num_dim*etot_est)
      call blank (cbc,3*2*num_dim*etot_est)
 
 
      write(6,*) 'please input number of fluid exo files:'
      read (5,*) fnexo
      allocate (f_elem_exo  (fnexo))
 
      eacc = 0
      do iexo = 1,fnexo

      eacc_old = eacc

      write(6,*) 'please input exo file:'
      call read_input_name
      call exodus_read_new
 
      write(6,*) 'Please input shifting vector:'
      read (5,*)  shiftvector(1),shiftvector(2),shiftvector(3)

      write(6,*) 'Please input option for splitting:'
      write(6,*) '1: assume pure hex20 mesh, original exo2nek'
      write(6,*) '2: assume hybrid tetra4 and wedge6 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all wedge6 elements in block 2 (or second block)'
      write(6,*) '3: assume hybrid hex8, tetra4 and wedge6 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all hex8   elements in block 2 (or second block)'
      write(6,*) '    all wedge6 elements in block 3 (or third block)'

      read (5,'(I1)') option

      if(option.EQ.1) then
      call convert_new
      elseif(option.EQ.2) then
      call split_convert1
      elseif(option.EQ.3) then
      call split_convert2
      else
       write(6,*) 'Unknown option for exo2nek, ABORT.'
      endif

      f_elem_exo(iexo) = eacc - eacc_old

      if (fnexo.gt.1)  call offset_sideset(iexo)

      call checkXYZ_min_max()
      write(6,*) 'total element now is ',eacc
      write(6,*) 'fluid exo file ',iexo,' has elements ',f_elem_exo(iexo)
  
      enddo

      etot = eacc
      eftot = etot
 
      write(6,*) 'please input number of solid exo files for CHT problem (input 0 for no solid mesh):'
      read (5,*) snexo

      if (snexo.gt.0) then
      allocate (s_elem_exo  (snexo))

      do iexo = 1,snexo
      eacc_old = eacc

      write(6,*) 'please input exo file:'
      call read_input_name
      call exodus_read_new

      write(6,*) 'Please input shifting vector (if no shift, just input 0 0 0):'
      read (5,*)  shiftvector(1),shiftvector(2),shiftvector(3)

      write(6,*) 'Please input option for splitting:'
      write(6,*) '1: assume pure hex20 mesh, original exo2nek'
      write(6,*) '2: assume hybrid tetra4 and wedge6 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all wedge6 elements in block 2 (or second block)'
      write(6,*) '3: assume hybrid hex8, tetra4 and wedge6 mesh'
      write(6,*) '    all tetra4 elements in block 1 (or first block)'
      write(6,*) '    all hex8   elements in block 2 (or second block)'
      write(6,*) '    all wedge6 elements in block 3 (or third block)'

      read (5,'(I1)') option

      if(option.EQ.1) then
      call convert_new
      elseif(option.EQ.2) then
      call split_convert1
      elseif(option.EQ.3) then
      call split_convert2
      else
       write(6,*) 'Unknown option for exo2nek, ABORT.'
      endif

      s_elem_exo(iexo) = eacc - eacc_old
      iexo1 = iexo + fnexo
      call offset_sideset(iexo1)

      call checkXYZ_min_max()
      write(6,*) 'total element now is ',eacc
      write(6,*) 'solid exo file ',iexo,' has elements ',s_elem_exo(iexo)

      enddo

     endif ! if (snexo.gt.0) then

      etot = eacc
      num_elem = etot

     if(etot.gt.etot_est) then
       write(6,*) 'ABORT'
       write(6,*) 'please increase estimate element number to ',etot
       call exit
      endif 

      call right_hand_check
      call gather_bc_info
      call setbc_new
      
      call scale_mesh

      write(6,*) 'please input output re2 file name:'
      call read_re2_name

      call gen_re2

      end 
!-----------------------------------------------------------------------
      subroutine checkXYZ_min_max
! return element bound for exo file iexo
      use SIZE

      maxxyz(1) = -1e6
      maxxyz(2) = -1e6
      maxxyz(3) = -1e6

      minxyz(1) = 1e6
      minxyz(2) = 1e6 
      minxyz(3) = 1e6 

      ntot = eacc*3*3*3

      do i = 1,ntot
        xx = xm1(i,1,1,1)        
        yy = ym1(i,1,1,1)
        zz = zm1(i,1,1,1)

        if(xx.gt.maxxyz(1)) maxxyz(1)  = xx
        if(yy.gt.maxxyz(2)) maxxyz(2)  = yy
        if(zz.gt.maxxyz(3)) maxxyz(3)  = zz

        if(xx.lt.minxyz(1)) minxyz(1)  = xx
        if(yy.lt.minxyz(2)) minxyz(2)  = yy
        if(zz.lt.minxyz(3)) minxyz(3)  = zz

      enddo

      write(6,*) 'Domain max xyz:',maxxyz(1),maxxyz(2),maxxyz(3) 
      write(6,*) 'Domain min xyz:',minxyz(1),minxyz(2),minxyz(3) 

      return
      end
!------------------------------------------------------------------------
      subroutine read_input_name

      use SIZE

!      character(1)  re2nam1(80)
      character(1)  exonam1(32)
      character(32) fname

      read (5,'(A32)') fname
      len = ltrunc(fname,32)
      
      call blank  (exonam1, 32)
!      call blank  (re2nam1, 80)
      call chcopy (exonam1,fname,32)
!      call chcopy (re2nam1,fname,80)
      call chcopy (exonam1(len+1) ,'.exo',4)
!      call chcopy (re2nam1(len+1) ,'.re2',4)

      call blank  (exoname, 32)
!      call blank  (re2name, 80)
      call chcopy (exoname,exonam1,len+4)
!      call chcopy (re2name,re2nam1,len+4)
 
      return 
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine read_re2_name

      use SIZE

      character(1)  re2nam1(80)
!      character(1)  exonam1(32)
      character(32) fname

      read (5,'(A32)') fname
      len = ltrunc(fname,32)
      
!      call blank  (exonam1, 32)
      call blank  (re2nam1, 80)
!      call chcopy (exonam1,fname,32)
      call chcopy (re2nam1,fname,80)
!      call chcopy (exonam1(len+1) ,'.exo',4)
      call chcopy (re2nam1(len+1) ,'.re2',4)

!      call blank  (exoname, 32)
      call blank  (re2name, 80)
!      call chcopy (exoname,exonam1,len+4)
      call chcopy (re2name,re2nam1,len+4)
 
      return 
      end
!-----------------------------------------------------------------------
      subroutine exodus_read_new
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

!      integer,allocatable,dimension(:) :: idblk
!      integer,allocatable,dimension(:) :: num_nodes_per_elem
!      integer,allocatable,dimension(:) :: num_attr    !not used

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
!      allocate ( xm1                (3,3,3,num_elem*8)      )
!      allocate ( ym1                (3,3,3,num_elem*8)      )
!      allocate ( zm1                (3,3,3,num_elem*8)      )

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

!        if (i.eq.1) then
!          nvert=num_nodes_per_elem(i)
!          if (num_dim.eq.2) then
!            if (nvert.ne.8) then
!              if (nvert.eq.9) then
!                write(6,*)
!                write(6,'(a)') &
!                "WARNING: QUAD9 elements are not officially supported"
!                write(6,'(a)') &
!                "as there is no exodus standard for this element type."
!                write(6,*)
!              else
!                write(6,*)
!                write(6,'(a)') &
!                "ERROR: Only QUAD8 elements are allowed in a 2D mesh!"
!                STOP
!              endif
!            endif      
!          elseif (num_dim.eq.3) then
!            if (nvert.ne.20) then
!              if (nvert.eq.27) then
!                write(6,*)
!                write(6,'(a)') &
!                "WARNING: HEX27 elements are not officially supported"
!                write(6,'(a)') &
!                "as there is no exodus standard for this element type."
!                write(6,*)
!              else
!                write(6,*)
!                write(6,'(a)') &
!                "ERROR: Only HEX20 elements are allowed in a 3D mesh!"
!                STOP
!              endif
!            endif      
!          else
!          write(6,'(a,i3)') &
!           "ERROR: Unknown number of dimensions! num_dim= ",num_dim
!          STOP
!          endif
!        endif
!
!        if (num_nodes_per_elem(i).ne.nvert) then
!          write(6,*)
!          write(6,'(a)') &
!           "ERROR: All blocks should contain elements of the same type!"
!          write(6,'(a,i3)') "num_nodes_per_elem= ",num_nodes_per_elem(i)
!          write(6,'(a,i3)') "num_nodes_per_elem of block 1 ",nvert
!          STOP
!        endif
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
! -------------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine convert_new
!
!  Subroutine to convert an already read exodusII mesh to a nek
!  mesh. The idea is to fill each element's node coordinates of
!  size lx1**3 (3D) or lx1**2 (2D) (lx1=3) with the hex27/quad9
!  coordinates.
!
      use SIZE
      include 'exodusII.inc'
      integer ssID
      logical ifinterface

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

      eacc_old = eacc
      nvert = 20
      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '
      do iel=1,num_elem
        eacc = eacc + 1
        do ivert =1,nvert
          if (num_dim.eq.2) then
            jvert = exo_to_nek_vert2D(ivert)
            xm1(jvert,1,1,eacc)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,eacc)=y_exo(connect(nvert*(iel-1)+ivert))
          else
            jvert = exo_to_nek_vert3D(ivert)
            xm1(jvert,1,1,eacc)=x_exo(connect(nvert*(iel-1)+ivert))
            ym1(jvert,1,1,eacc)=y_exo(connect(nvert*(iel-1)+ivert))
            zm1(jvert,1,1,eacc)=z_exo(connect(nvert*(iel-1)+ivert))

            xm1(jvert,1,1,eacc)=xm1(jvert,1,1,eacc) + shiftvector(1)
            ym1(jvert,1,1,eacc)=ym1(jvert,1,1,eacc) + shiftvector(2)
            zm1(jvert,1,1,eacc)=zm1(jvert,1,1,eacc) + shiftvector(3)
          endif
        enddo
      enddo

      write(6,'(A)') 'done :: Converting elements '
!
! allocate and zero-out curve sides arrays
!
!      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
!      allocate   (curve  (2*num_dim,12,   num_elem) )
!      call rzero (curve,2*num_dim*12*num_elem)
!      call blank (ccurve,(4+8*(num_dim-2))*num_elem)
!
! allocate and zero-out bc arrays only if sidesets are specified
!
!      if (num_side_sets.eq.0) return   

!      allocate   (cbc    (2*num_dim,      num_elem) )
!      allocate   (bc     (5,2*num_dim,    num_elem) ) 
!      call rzero (bc,5*2*num_dim*num_elem)
!      call blank (cbc,3*2*num_dim*num_elem)

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
           cbc(jfc,iel+eacc_old)   = 'EXO' ! dummy exodus bc
           bc (5,jfc,iel+eacc_old) = idss(iss)
           bc (4,jfc,iel+eacc_old) = iexo
        enddo

      write(6,'(A,I2)') 'done :: Sideset ',idss(iss)
      enddo

      deallocate(x_exo,y_exo,z_exo,connect)
      deallocate (elem_list,side_list)
      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
      deallocate ( num_sides_in_set   )
      deallocate ( idss               )

      write(6,'(a)') ''
      write(6,'(a)') 'done :: Converting SideSets '

      return
      end
!--------------------------------------------------------------------
!-----------------------------------------------------------------
      subroutine split_convert1
! split 1 TETRA4 to 4 HEX20
! split 1 WEDGE6 to 3 HEX20

      use SIZE
      include 'exodusII.inc'

!  in fortran arraysm higher order is at back.
      real*8 hexver(3,27,8) ! hex coordinates
      real*8 tetver(3,10),wedgever(3,15)
      real*8 ehexver(3,20)

!  exoss is used to store sideset information for all exo elements.
!  tet,hex,wedge
      integer exoss(6,num_elem) 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)
      integer ehexnumber,tetnumber,wedgenumber,bctot
      integer vert_index_exo
      integer iel_nek,iel_exo,ifc_exo

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo

      eacc_old = eacc

      call rzero_int(exoss,6*max_num_elem)
! store sideset information
      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Store SideSet information from EXO file'
        do iss=1,num_side_sets   ! loop over ss
           if(idss(iss).gt.0) then
           write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel_exo = elem_list(i,iss) ! exo element number
             ifc_exo = side_list(i,iss) ! exo element face number
             exoss(ifc_exo,iel_exo) = idss(iss)
           enddo
           endif
        enddo
      endif

      !deallocate(elem_list,side_list)

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '

!      allocate   (ccurve (4+8*(num_dim-2),num_elem*4) )
!      allocate   (curve  (2*num_dim,12,   num_elem*4) )
!      call rzero (curve,2*num_dim*12*num_elem*4)
!      call blank (ccurve,(4+8*(num_dim-2))*num_elem*4)
!
!      allocate   (cbc    (2*num_dim,      num_elem*4) )
!      allocate   (bc     (5,2*num_dim,    num_elem*4) ) 
!      call rzero (bc,5*2*num_dim*num_elem*4)
!      call blank (cbc,3*2*num_dim*num_elem*4)


!  assume block 1 (or 1st block) contains all tet4 elements
!  assume block 2 (or 2nd block) contains all wedge6 elements
      call rzero_int(hexss,6*8) 
      tetnumber  = num_elem_in_block(1)
      vert_index_exo = 0
      iel_nek = eacc_old !0

      write(6,'(A)') 'Splitting elements ... '
      do iel_exo = 1, num_elem  

!! number of elements == 4*number of tet + 3*number of wedge

!  if tet.
       if (iel_exo.le.tetnumber) then

!  read tet 4 . 
!  linear interpolate to tet10
       do ivert = 1, 4
       vert_index_exo = vert_index_exo + 1  
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(tetver(1,5),tetver(1,1),tetver(1,2))
       call average2vec(tetver(1,6),tetver(1,2),tetver(1,3))
       call average2vec(tetver(1,7),tetver(1,1),tetver(1,3))
       call average2vec(tetver(1,8),tetver(1,1),tetver(1,4))
       call average2vec(tetver(1,9),tetver(1,2),tetver(1,4))
       call average2vec(tetver(1,10),tetver(1,3),tetver(1,4))

!  assign sideset to tet elements
       call rzero_int(tetss,4)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
         tetss(ifc_exo) = 0
         tetss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif
!  given tet10 vertices, and return you four hex27 coords.
       call tettohex(hexver,tetver,hexss,tetss)

       do ihex = 1, 4
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          !write(6,*)  iel_nek
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)& 
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)&
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)&
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      else ! if wedge elements 

!  read wedge 6,
!  linear interpolate to wedge 15
       do ivert = 1,6
       vert_index_exo = vert_index_exo + 1  
       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(wedgever(1,7),wedgever(1,1), &
       wedgever(1,2))
       call average2vec(wedgever(1,8),wedgever(1,2), &
       wedgever(1,3))
       call average2vec(wedgever(1,9),wedgever(1,1), &
       wedgever(1,3))
       call average2vec(wedgever(1,10),wedgever(1,1), &
       wedgever(1,4))
       call average2vec(wedgever(1,11),wedgever(1,2), &
       wedgever(1,5))
       call average2vec(wedgever(1,12),wedgever(1,3), &
      wedgever(1,6))
       call average2vec(wedgever(1,13),wedgever(1,4), &
      wedgever(1,5))
       call average2vec(wedgever(1,14),wedgever(1,5), &
      wedgever(1,6))
       call average2vec(wedgever(1,15),wedgever(1,4), &
      wedgever(1,6))

!  assign sideset to wedge elements
       call rzero_int(wedgess,5)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

!  given wedge15 vertices, and return you 3 hex coords.
       call wedgetohex(hexver,wedgever,hexss,wedgess)
       
       do ihex = 1,3
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)& 
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)&
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      endif
   
      enddo
      write(6,'(A)') ' Done :: Splitting elements ... '

      deallocate (x_exo,y_exo,z_exo,connect)
      deallocate (elem_list,side_list)
      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
      deallocate ( num_sides_in_set   )
      deallocate ( idss               )

      write(6,*) 'Converted elements in nek:',iel_nek-eacc_old
      write(6,'(A)') 'Done :: Converting elements '
 
      !num_elem = iel_nek
	  
      !bctot = 0
      !do ihex = 1,num_elem
      !  do ifc=1,6
      !     if (cbc(ifc,ihex).eq.'EXO') then
      !        bctot = bctot + 1 
      !     endif		
      !  enddo
      !enddo
      !write(6,*) 'Converted elem sides with BC :',bctot

      return
      end
!--------------------------------------------------------------------
      subroutine split_convert2
! split 1 TETRA4 to 4 HEX20, block 1
! split 1 exo HEX8 to 8 hex20, block 2
! split 1 WEDGE6 to 6 HEX20, block 3
!
      use SIZE
      include 'exodusII.inc'

!  in fortran arrays higher order is at back.
      real*8 hexver(3,27,8) ! hex coordinates
      real*8 tetver(3,10),wedgever(3,15)
      real*8 ehexver(3,20)

!  exoss is used to store sideset information for all exo elements.
!  tet,hex,wedge
      integer exoss(6,num_elem) 

      integer tetss(4),wedgess(5),ehexss(6),hexss(6,8)
      integer ehexnumber,tetnumber,wedgenumber,bctot
      integer vert_index_exo
      integer iel_nek,iel_exo,ifc_exo

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo
      eacc_old = eacc
      call rzero_int(exoss,6*max_num_elem)
! store sideset information
      if (num_side_sets.ne.0) then
      write(6,'(a)') ''
      write(6,'(a)') 'Store SideSet information from EXO file'
        do iss=1,num_side_sets   ! loop over ss 
        write(6,'(a,i2,a)') 'Sideset ',idss(iss), ' ...'
           do i=1,num_sides_in_set(iss) ! loop over sides in ss
             iel_exo = elem_list(i,iss) ! exo element number
             ifc_exo = side_list(i,iss) ! exo element face number
             exoss(ifc_exo,iel_exo) = idss(iss)
           enddo
        enddo
      endif

      write(6,'(A)') ' '
      write(6,'(A)') 'Converting elements ... '

!      allocate   (ccurve (4+8*(num_dim-2),num_elem*8) )
!      allocate   (curve  (2*num_dim,12,   num_elem*8) )
!      call rzero (curve,2*num_dim*12*num_elem*8)
!      call blank (ccurve,(4+8*(num_dim-2))*num_elem*8)
!	  
!      allocate   (cbc    (2*num_dim,      num_elem*8) )
!      allocate   (bc     (5,2*num_dim,    num_elem*8) ) 
!      call rzero (bc,5*2*num_dim*num_elem*8)
!      call blank (cbc,3*2*num_dim*num_elem*8)

      call rzero_int(hexss,6*8) 
      tetnumber  = num_elem_in_block(1) ! all exo tet elements in block 1 (or first block)
      ehexnumber = num_elem_in_block(2) ! all exo hex elements in block 2 (or second block)
      wedgenumber = num_elem_in_block(3) ! all exo wedge elements in block 3 (or third block)

      vert_index_exo = 0
      iel_nek =  eacc_old !0

      do iel_exo = 1, num_elem  

!! number of elements == 4*number of tet + 6*number of wedge + 8*number of hex

!  block 1, if tet.
       if (iel_exo.le.tetnumber) then

!  read tet 4 . 
!  linear interpolate to tet10
       do ivert = 1, 4
       vert_index_exo = vert_index_exo + 1  
       tetver(1,ivert) = x_exo(connect(vert_index_exo))
       tetver(2,ivert) = y_exo(connect(vert_index_exo))
       tetver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(tetver(1,5),tetver(1,1),tetver(1,2))
       call average2vec(tetver(1,6),tetver(1,2),tetver(1,3))
       call average2vec(tetver(1,7),tetver(1,1),tetver(1,3))
       call average2vec(tetver(1,8),tetver(1,1),tetver(1,4))
       call average2vec(tetver(1,9),tetver(1,2),tetver(1,4))
       call average2vec(tetver(1,10),tetver(1,3),tetver(1,4))


!  assign sideset to tet elements
       call rzero_int(tetss,4)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
          tetss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

!  given tet10 vertices, and return you four hex coords.
       call tettohex(hexver,tetver,hexss,tetss) 

       do ihex = 1,4
          iel_nek = iel_nek + 1
           eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)&
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)&
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

!  block 2, if exo hex element, convert to 8 nek hex elements
       elseif (iel_exo.gt.tetnumber &
      .and.(iel_exo.le.(tetnumber+ehexnumber))) then

!  read hex8
!  linear interpolate to hex20	   
       do ivert = 1,8
       vert_index_exo = vert_index_exo + 1  
       ehexver(1,ivert) = x_exo(connect(vert_index_exo))
       ehexver(2,ivert) = y_exo(connect(vert_index_exo))
       ehexver(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(ehexver(1,9),ehexver(1,1),ehexver(1,2))
       call average2vec(ehexver(1,10),ehexver(1,2),ehexver(1,3))
       call average2vec(ehexver(1,11),ehexver(1,3),ehexver(1,4))
       call average2vec(ehexver(1,12),ehexver(1,1),ehexver(1,4))
       call average2vec(ehexver(1,13),ehexver(1,1),ehexver(1,5))
       call average2vec(ehexver(1,14),ehexver(1,2),ehexver(1,6))
       call average2vec(ehexver(1,15),ehexver(1,3),ehexver(1,7))
       call average2vec(ehexver(1,16),ehexver(1,4),ehexver(1,8))
       call average2vec(ehexver(1,17),ehexver(1,5),ehexver(1,6))
       call average2vec(ehexver(1,18),ehexver(1,6),ehexver(1,7))
       call average2vec(ehexver(1,19),ehexver(1,7),ehexver(1,8))
       call average2vec(ehexver(1,20),ehexver(1,5),ehexver(1,8))

!  assign sideset
       call rzero_int(ehexss,6)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,6
          ehexss(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

       call ehexto8hex(hexver,ehexver,hexss,ehexss) 

       do ihex = 1,8
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)& 
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)& 
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

! block 3, if WEDGE15 elements
      elseif (iel_exo.gt.(tetnumber+ehexnumber)) then
!  read wedge 6,
!  linear interpolate to wedge 15
       do ivert = 1,6
       vert_index_exo = vert_index_exo + 1  
       wedgever(1,ivert) = x_exo(connect(vert_index_exo))
       wedgever(2,ivert) = y_exo(connect(vert_index_exo))
       wedgever(3,ivert) = z_exo(connect(vert_index_exo))
       enddo

       call average2vec(wedgever(1,7),wedgever(1,1),&
       wedgever(1,2))
       call average2vec(wedgever(1,8),wedgever(1,2),&
       wedgever(1,3))
       call average2vec(wedgever(1,9),wedgever(1,1),&
       wedgever(1,3))
       call average2vec(wedgever(1,10),wedgever(1,1),&
       wedgever(1,4))
       call average2vec(wedgever(1,11),wedgever(1,2),&
       wedgever(1,5))
       call average2vec(wedgever(1,12),wedgever(1,3),&
       wedgever(1,6))
       call average2vec(wedgever(1,13),wedgever(1,4),&
       wedgever(1,5))
       call average2vec(wedgever(1,14),wedgever(1,5),&
       wedgever(1,6))
       call average2vec(wedgever(1,15),wedgever(1,4),&
       wedgever(1,6))


!  assign sideset to wedge elements
       call rzero_int(wedgess,5)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo)
        enddo
       endif

!  given wedge15 vertices, and return you 6 hex coords.
       call wedgetohex2(hexver,wedgever,hexss,wedgess)

       do ihex = 1,6
          iel_nek = iel_nek + 1
          eacc = eacc + 1
          do inekvert = 1,27
             xm1(inekvert,1,1,iel_nek) = hexver(1,inekvert,ihex)&
       + shiftvector(1)
             ym1(inekvert,1,1,iel_nek) = hexver(2,inekvert,ihex)& 
       + shiftvector(2)
             zm1(inekvert,1,1,iel_nek) = hexver(3,inekvert,ihex)& 
       + shiftvector(3)
          enddo
          do ifc=1,6
             if(hexss(ifc,ihex).gt.0) then
              bc(5,ifc,iel_nek) = hexss(ifc,ihex)
              bc(4,ifc,iel_nek) = iexo
              cbc(ifc,iel_nek)   = 'EXO' ! dummy boundary condition
             endif
          enddo
       enddo

      endif
      enddo

      write(6,*) 'Converted elements in nek:',iel_nek-eacc_old
      write(6,'(A)') 'Done :: Converting elements '

      deallocate(x_exo,y_exo,z_exo,connect)
      deallocate (elem_list,side_list)
      deallocate ( idblk )
      deallocate ( num_nodes_per_elem )
      deallocate ( num_attr           )
      deallocate ( num_elem_in_block  )
      deallocate ( num_sides_in_set   )
      deallocate ( idss               )

      !num_elem = iel_nek
      !
      !bctot = 0
      !do ihex = 1,num_elem
      !  do ifc=1,6
      !     if (cbc(ifc,ihex).eq.'EXO') then
      !        bctot = bctot + 1 
      !     endif		
      !  enddo
      !enddo
      !write(6,*) 'Converted elem sides with BC :',bctot
	  
      return
      end
!--------------------------------------------------------------------
      subroutine tettohex(hexver,tetver,hexss,tetss)
      real*8 tetver(3,10) ! tet vertices
      real*8 tetface(3,4) ! tet face center
      real*8 tetcen(3,1)  ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer tetss(4),hexss(6,4)

      do i=1,6*4
      hexss(i,1)=0
      enddo

!  get face center coords.
      call average3vec(tetface(1,1),tetver(1,1),tetver(1,2),tetver(1,4))
      call average3vec(tetface(1,2),tetver(1,2),tetver(1,3),tetver(1,4))
      call average3vec(tetface(1,3),tetver(1,1),tetver(1,3),tetver(1,4))
      call average3vec(tetface(1,4),tetver(1,1),tetver(1,2),tetver(1,3))

!  get tet vol center
      call average4vec(tetcen(1,1),tetver(1,1),tetver(1,2),&
      tetver(1,3),tetver(1,4))

!  assign coordinates to four hex.
!  hex 1
      hexss(1,1) = tetss(1)
      hexss(4,1) = tetss(3)
      hexss(5,1) = tetss(4)

      call assignvec(hex8(1,1),tetver(1,1))
      call assignvec(hex8(1,2),tetver(1,5))
      call assignvec(hex8(1,3),tetface(1,4))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetver(1,8))
      call assignvec(hex8(1,6),tetface(1,1))
      call assignvec(hex8(1,7),tetcen(1,1))
      call assignvec(hex8(1,8),tetface(1,3))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 
!  hex 2
      hexss(1,2) = tetss(1)
      hexss(2,2) = tetss(2)
      hexss(5,2) = tetss(4)

      call assignvec(hex8(1,1),tetver(1,5))
      call assignvec(hex8(1,2),tetver(1,2))
      call assignvec(hex8(1,3),tetver(1,6))
      call assignvec(hex8(1,4),tetface(1,4))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetface(1,2))
      call assignvec(hex8(1,8),tetcen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
!  hex 3
      hexss(2,3) = tetss(2)
      hexss(3,3) = tetss(3)
      hexss(5,3) = tetss(4)

      call assignvec(hex8(1,1),tetface(1,4))
      call assignvec(hex8(1,2),tetver(1,6))
      call assignvec(hex8(1,3),tetver(1,3))
      call assignvec(hex8(1,4),tetver(1,7))
      call assignvec(hex8(1,5),tetcen(1,1))
      call assignvec(hex8(1,6),tetface(1,2))
      call assignvec(hex8(1,7),tetver(1,10))
      call assignvec(hex8(1,8),tetface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

!  hex 4
      hexss(2,4) = tetss(2)
      hexss(3,4) = tetss(3)
      hexss(6,4) = tetss(1)

      call assignvec(hex8(1,1),tetcen(1,1))
      call assignvec(hex8(1,2),tetface(1,2))
      call assignvec(hex8(1,3),tetver(1,10))
      call assignvec(hex8(1,4),tetface(1,3))
      call assignvec(hex8(1,5),tetface(1,1))
      call assignvec(hex8(1,6),tetver(1,9))
      call assignvec(hex8(1,7),tetver(1,4))
      call assignvec(hex8(1,8),tetver(1,8))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))

      return
      end	  
!--------------------------------------------------------------------
      subroutine wedgetohex(hexver,wedgever,hexss,wedgess)
!  convert 1 wedge to 3 hex
      real*8 wedgever(3,15) ! tet vertices
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,4) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,4)

      do i=1,6*4
      hexss(i,1)=0
      enddo

      call average3vec(wedgeface(1,4),wedgever(1,1),wedgever(1,2),&
      wedgever(1,3))
      call average3vec(wedgeface(1,5),wedgever(1,4),wedgever(1,5),&
      wedgever(1,6))

!  assign coordinates to 3 hex.
!  hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)
      hexss(6,1) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),wedgeface(1,4))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),wedgeface(1,5))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 
!  hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)
      hexss(6,2) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),wedgeface(1,4))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),wedgeface(1,5))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
!  hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)
      hexss(6,3) = wedgess(5)

      call assignvec(hex8(1,1),wedgeface(1,4))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgeface(1,5))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

      return
      end	  
! -------------------------------------------------------------------
      subroutine wedgetohex2(hexver,wedgever,hexss,wedgess)
!  convert 1 wedge to 6 nek hex20 elements
      real*8 wedgever(3,15) ! tet vertices
      real*8 wedgeface(3,5) ! tet face center
      real*8 wedgecen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,8) ! four hex vertices in nek format.
      real*8 tempvec(3,3) ! temperary vector variable

      integer wedgess(5),hexss(6,8)

      do i=1,6*8
      hexss(i,1)=0
      enddo

      call average4vec(wedgeface(1,1),wedgever(1,1),wedgever(1,2),&
      wedgever(1,4),wedgever(1,5))	 
      call average4vec(wedgeface(1,2),wedgever(1,2),wedgever(1,3),&
      wedgever(1,5),wedgever(1,6))
      call average4vec(wedgeface(1,3),wedgever(1,1),wedgever(1,3),&
      wedgever(1,4),wedgever(1,6))
      call average3vec(wedgeface(1,4),wedgever(1,1),wedgever(1,2),&
      wedgever(1,3))
      call average3vec(wedgeface(1,5),wedgever(1,4),wedgever(1,5),&
      wedgever(1,6))
      call average2vec(wedgecen(1,1),wedgeface(1,4),wedgeface(1,5))

!  assign coordinates to 6 hex.
!  hex 1
      hexss(1,1) = wedgess(1)
      hexss(4,1) = wedgess(3)
      hexss(5,1) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,1))
      call assignvec(hex8(1,2),wedgever(1,7))
      call assignvec(hex8(1,3),wedgeface(1,4))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgever(1,10))
      call assignvec(hex8(1,6),wedgeface(1,1))
      call assignvec(hex8(1,7),wedgecen(1,1))
      call assignvec(hex8(1,8),wedgeface(1,3))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))

!  hex 2
      hexss(1,2) = wedgess(1)
      hexss(2,2) = wedgess(2)
      hexss(5,2) = wedgess(4)

      call assignvec(hex8(1,1),wedgever(1,7))
      call assignvec(hex8(1,2),wedgever(1,2))
      call assignvec(hex8(1,3),wedgever(1,8))
      call assignvec(hex8(1,4),wedgeface(1,4))
      call assignvec(hex8(1,5),wedgeface(1,1))
      call assignvec(hex8(1,6),wedgever(1,11))
      call assignvec(hex8(1,7),wedgeface(1,2))
      call assignvec(hex8(1,8),wedgecen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
!  hex 3
      hexss(2,3) = wedgess(2)
      hexss(3,3) = wedgess(3)
      hexss(5,3) = wedgess(4)

      call assignvec(hex8(1,1),wedgeface(1,4))
      call assignvec(hex8(1,2),wedgever(1,8))
      call assignvec(hex8(1,3),wedgever(1,3))
      call assignvec(hex8(1,4),wedgever(1,9))
      call assignvec(hex8(1,5),wedgecen(1,1))
      call assignvec(hex8(1,6),wedgeface(1,2))
      call assignvec(hex8(1,7),wedgever(1,12))
      call assignvec(hex8(1,8),wedgeface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

!  hex 4
      hexss(1,4) = wedgess(1)
      hexss(4,4) = wedgess(3)
      hexss(6,4) = wedgess(5)

      call assignvec(hex8(1,1),wedgever(1,10))
      call assignvec(hex8(1,2),wedgeface(1,1))
      call assignvec(hex8(1,3),wedgecen(1,1))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgever(1,4))
      call assignvec(hex8(1,6),wedgever(1,13))
      call assignvec(hex8(1,7),wedgeface(1,5))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))	 
!  hex 5
      hexss(1,5) = wedgess(1)
      hexss(2,5) = wedgess(2)
      hexss(6,5) = wedgess(5)

      call assignvec(hex8(1,1),wedgeface(1,1))
      call assignvec(hex8(1,2),wedgever(1,11))
      call assignvec(hex8(1,3),wedgeface(1,2))
      call assignvec(hex8(1,4),wedgecen(1,1))
      call assignvec(hex8(1,5),wedgever(1,13))
      call assignvec(hex8(1,6),wedgever(1,5))
      call assignvec(hex8(1,7),wedgever(1,14))
      call assignvec(hex8(1,8),wedgeface(1,5))
      call hex8tohex27(hexver(1,1,5),hex8(1,1))
!  hex 6
      hexss(2,6) = wedgess(2)
      hexss(3,6) = wedgess(3)
      hexss(6,6) = wedgess(5)

      call assignvec(hex8(1,1),wedgecen(1,1))
      call assignvec(hex8(1,2),wedgeface(1,2))
      call assignvec(hex8(1,3),wedgever(1,12))
      call assignvec(hex8(1,4),wedgeface(1,3))
      call assignvec(hex8(1,5),wedgeface(1,5))
      call assignvec(hex8(1,6),wedgever(1,14))
      call assignvec(hex8(1,7),wedgever(1,6))
      call assignvec(hex8(1,8),wedgever(1,15))
      call hex8tohex27(hexver(1,1,6),hex8(1,1))

      return
      end	  
! -------------------------------------------------------------------
      subroutine ehexto8hex(hexver,ehexver,hexss,ehexss)
!  convert exodus hex20 elements 
!  to 8 nek hex20 elements
      real*8 ehexver(3,20)
      real*8 ehexface(3,6)
      real*8 ehexcen(3,1) ! tet vol center
      real*8 hex8(3,8)
      real*8 hexver(3,27,8)
      integer ehexss(6),hexss(6,8)

      do i=1,6*8
      hexss(i,1)=0
      enddo

      call average4vec(ehexface(1,1),ehexver(1,1),ehexver(1,2),&
      ehexver(1,5),ehexver(1,6))
      call average4vec(ehexface(1,2),ehexver(1,2),ehexver(1,3),&
      ehexver(1,6),ehexver(1,7))
      call average4vec(ehexface(1,3),ehexver(1,3),ehexver(1,4),&
      ehexver(1,7),ehexver(1,8))
      call average4vec(ehexface(1,4),ehexver(1,1),ehexver(1,4),&
      ehexver(1,5),ehexver(1,8))
      call average4vec(ehexface(1,5),ehexver(1,1),ehexver(1,2),&
      ehexver(1,3),ehexver(1,4))
      call average4vec(ehexface(1,6),ehexver(1,5),ehexver(1,6),&
      ehexver(1,7),ehexver(1,8))
      call average2vec(ehexcen(1,1),ehexface(1,5),ehexface(1,6))

!  assign coordinates to 8 hex.
!  hex 1
      hexss(1,1) = ehexss(1)
      hexss(4,1) = ehexss(4)
      hexss(5,1) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,1))
      call assignvec(hex8(1,2),ehexver(1,9))
      call assignvec(hex8(1,3),ehexface(1,5))
      call assignvec(hex8(1,4),ehexver(1,12))
      call assignvec(hex8(1,5),ehexver(1,13))
      call assignvec(hex8(1,6),ehexface(1,1))
      call assignvec(hex8(1,7),ehexcen(1,1))
      call assignvec(hex8(1,8),ehexface(1,4))
      call hex8tohex27(hexver(1,1,1),hex8(1,1))	 

!  hex 2
      hexss(1,2) = ehexss(1)
      hexss(2,2) = ehexss(2)
      hexss(5,2) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,9))
      call assignvec(hex8(1,2),ehexver(1,2))
      call assignvec(hex8(1,3),ehexver(1,10))
      call assignvec(hex8(1,4),ehexface(1,5))
      call assignvec(hex8(1,5),ehexface(1,1))
      call assignvec(hex8(1,6),ehexver(1,14))
      call assignvec(hex8(1,7),ehexface(1,2))
      call assignvec(hex8(1,8),ehexcen(1,1))
      call hex8tohex27(hexver(1,1,2),hex8(1,1))
!  hex 3
      hexss(2,3) = ehexss(2)
      hexss(3,3) = ehexss(3)
      hexss(5,3) = ehexss(5)

      call assignvec(hex8(1,1),ehexface(1,5))
      call assignvec(hex8(1,2),ehexver(1,10))
      call assignvec(hex8(1,3),ehexver(1,3))
      call assignvec(hex8(1,4),ehexver(1,11))
      call assignvec(hex8(1,5),ehexcen(1,1))
      call assignvec(hex8(1,6),ehexface(1,2))
      call assignvec(hex8(1,7),ehexver(1,15))
      call assignvec(hex8(1,8),ehexface(1,3))
      call hex8tohex27(hexver(1,1,3),hex8(1,1))

!  hex 4
      hexss(3,4) = ehexss(3)
      hexss(4,4) = ehexss(4)
      hexss(5,4) = ehexss(5)

      call assignvec(hex8(1,1),ehexver(1,12))
      call assignvec(hex8(1,2),ehexface(1,5))
      call assignvec(hex8(1,3),ehexver(1,11))
      call assignvec(hex8(1,4),ehexver(1,4))
      call assignvec(hex8(1,5),ehexface(1,4))
      call assignvec(hex8(1,6),ehexcen(1,1))
      call assignvec(hex8(1,7),ehexface(1,3))
      call assignvec(hex8(1,8),ehexver(1,16))
      call hex8tohex27(hexver(1,1,4),hex8(1,1))

!  hex 5
      hexss(1,5) = ehexss(1)
      hexss(4,5) = ehexss(4)
      hexss(6,5) = ehexss(6)

      call assignvec(hex8(1,1),ehexver(1,13))
      call assignvec(hex8(1,2),ehexface(1,1))
      call assignvec(hex8(1,3),ehexcen(1,1))
      call assignvec(hex8(1,4),ehexface(1,4))
      call assignvec(hex8(1,5),ehexver(1,5))
      call assignvec(hex8(1,6),ehexver(1,17))
      call assignvec(hex8(1,7),ehexface(1,6))
      call assignvec(hex8(1,8),ehexver(1,20))
      call hex8tohex27(hexver(1,1,5),hex8(1,1))	 
!  hex 6
      hexss(1,6) = ehexss(1)
      hexss(2,6) = ehexss(2)
      hexss(6,6) = ehexss(6)

      call assignvec(hex8(1,1),ehexface(1,1))
      call assignvec(hex8(1,2),ehexver(1,14))
      call assignvec(hex8(1,3),ehexface(1,2))
      call assignvec(hex8(1,4),ehexcen(1,1))
      call assignvec(hex8(1,5),ehexver(1,17))
      call assignvec(hex8(1,6),ehexver(1,6))
      call assignvec(hex8(1,7),ehexver(1,18))
      call assignvec(hex8(1,8),ehexface(1,6))
      call hex8tohex27(hexver(1,1,6),hex8(1,1))
!  hex 7
      hexss(2,7) = ehexss(2)
      hexss(3,7) = ehexss(3)
      hexss(6,7) = ehexss(6)

      call assignvec(hex8(1,1),ehexcen(1,1))
      call assignvec(hex8(1,2),ehexface(1,2))
      call assignvec(hex8(1,3),ehexver(1,15))
      call assignvec(hex8(1,4),ehexface(1,3))
      call assignvec(hex8(1,5),ehexface(1,6))
      call assignvec(hex8(1,6),ehexver(1,18))
      call assignvec(hex8(1,7),ehexver(1,7))
      call assignvec(hex8(1,8),ehexver(1,19))
      call hex8tohex27(hexver(1,1,7),hex8(1,1))

!  hex 8
      hexss(3,8) = ehexss(3)
      hexss(4,8) = ehexss(4)
      hexss(6,8) = ehexss(6)

      call assignvec(hex8(1,1),ehexface(1,4))
      call assignvec(hex8(1,2),ehexcen(1,1))
      call assignvec(hex8(1,3),ehexface(1,3))
      call assignvec(hex8(1,4),ehexver(1,16))
      call assignvec(hex8(1,5),ehexver(1,20))
      call assignvec(hex8(1,6),ehexface(1,6))
      call assignvec(hex8(1,7),ehexver(1,19))
      call assignvec(hex8(1,8),ehexver(1,8))
      call hex8tohex27(hexver(1,1,8),hex8(1,1))	 

      return
      end
!--------------------------------------------------------------------
      subroutine hex8tohex27(hex27,hex8)
!  convert hex8 coordinates to hex27 coordinates in nek.
!  
      real*8 hex8(3,8)
      real*8 hex27(3,27)
      real*8 tempvec(3,3)

!   hex27 vert 1
      call assignvec(hex27(1,1),hex8(1,1))
!   hex27 vert 2
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,2))
      call assignvec(hex27(1,2),tempvec(1,1))
!   hex27 vert 3
      call assignvec(hex27(1,3),hex8(1,2))
!   hex27 vert 4
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,4))
      call assignvec(hex27(1,4),tempvec(1,1))
!   hex27 vert 5
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,2),hex8(1,3),hex8(1,4))
      call assignvec(hex27(1,5),tempvec(1,1))
!   hex27 vert 6
      call average2vec(tempvec(1,1),hex8(1,2),hex8(1,3))
      call assignvec(hex27(1,6),tempvec(1,1))
!   hex27 vert 7
      call assignvec(hex27(1,7),hex8(1,4))
!   hex27 vert 8
      call average2vec(tempvec(1,1),hex8(1,3),hex8(1,4))
      call assignvec(hex27(1,8),tempvec(1,1))
!   hex27 vert 9
      call assignvec(hex27(1,9),hex8(1,3))

!   hex27 vert 10
      call average2vec(tempvec(1,1),hex8(1,1),hex8(1,5))
      call assignvec(hex27(1,10),tempvec(1,1))  
!   hex27 vert 11
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,2),hex8(1,5),hex8(1,6))
      call assignvec(hex27(1,11),tempvec(1,1))
!   hex27 vert 12
      call average2vec(tempvec(1,1),hex8(1,2),hex8(1,6))
      call assignvec(hex27(1,12),tempvec(1,1))

!   hex27 vert 13
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,4),hex8(1,5),hex8(1,8))
      call assignvec(hex27(1,13),tempvec(1,1))
!   hex27 vert 14
      call average4vec(tempvec(1,1),hex8(1,1), &
      hex8(1,4),hex8(1,5),hex8(1,8))
      call average4vec(tempvec(1,2),hex8(1,2), &
      hex8(1,3),hex8(1,6),hex8(1,7))
      call average2vec(tempvec(1,3),tempvec(1,1),tempvec(1,2))
      call assignvec(hex27(1,14),tempvec(1,3))
!   hex27 vert 15
      call average4vec(tempvec(1,1),hex8(1,2), & 
      hex8(1,3),hex8(1,6),hex8(1,7))
      call assignvec(hex27(1,15),tempvec(1,1))

!   hex27 vert 16
      call average2vec(tempvec(1,1),hex8(1,4),hex8(1,8))
      call assignvec(hex27(1,16),tempvec(1,1))
!   hex27 vert 17
      call average4vec(tempvec(1,1),hex8(1,3), &
      hex8(1,4),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,17),tempvec(1,1))
!   hex27 vert 18
      call average2vec(tempvec(1,1),hex8(1,3),hex8(1,7))
      call assignvec(hex27(1,18),tempvec(1,1))

!   hex27 vert 19
      call assignvec(hex27(1,19),hex8(1,5))
!   hex27 vert 20
      call average2vec(tempvec(1,1),hex8(1,5),hex8(1,6))
      call assignvec(hex27(1,20),tempvec(1,1))	  
!   hex27 vert 21
      call assignvec(hex27(1,21),hex8(1,6))

!   hex27 vert 22
      call average2vec(tempvec(1,1),hex8(1,5),hex8(1,8))
      call assignvec(hex27(1,22),tempvec(1,1))	  
!   hex27 vert 23
      call average4vec(tempvec(1,1),hex8(1,5), &
      hex8(1,6),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,23),tempvec(1,1))
!   hex27 vert 24
      call average2vec(tempvec(1,1),hex8(1,6),hex8(1,7))
      call assignvec(hex27(1,24),tempvec(1,1))
	  
!   hex27 vert 25
      call assignvec(hex27(1,25),hex8(1,8))
!   hex27 vert 26
      call average2vec(tempvec(1,1),hex8(1,7),hex8(1,8))
      call assignvec(hex27(1,26),tempvec(1,1))	  
!   hex27 vert 27
      call assignvec(hex27(1,27),hex8(1,7))

      return
      end
! -------------------------------------------------------------------
      subroutine assignvec(a,b)
      real*8 a(3),b(3)
      a(1) = b(1)
      a(2) = b(2)
      a(3) = b(3)
      return
      end
! -------------------------------------------------------------------
      subroutine average2vec(avg,a,b)
      real*8 a(3),b(3),avg(3)      
      avg(1) = (a(1)+b(1))*0.5
      avg(2) = (a(2)+b(2))*0.5
      avg(3) = (a(3)+b(3))*0.5
      return
      end
!--------------------------------------------------------------------
      subroutine average3vec(avg,a,b,c)
      real*8 a(3),b(3),c(3),avg(3)      
      avg(1) = (a(1)+b(1)+c(1))*0.33333333
      avg(2) = (a(2)+b(2)+c(2))*0.33333333
      avg(3) = (a(3)+b(3)+c(3))*0.33333333

      return
      end
! -------------------------------------------------------------------
      subroutine average4vec(avg,a,b,c,d)
      real*8  a(3),b(3),c(3),d(3),avg(3)       
      avg(1) = (a(1)+b(1)+c(1)+d(1))*0.25
      avg(2) = (a(2)+b(2)+c(2)+d(2))*0.25
      avg(3) = (a(3)+b(3)+c(3)+d(3))*0.25
      return
      end
! -------------------------------------------------------------------
      subroutine setbc
!  set boundary condition 
!  read from a file casename.bc
      use SIZE
      parameter (npbc_max=10) ! maximum pairs of periodic boundary condition

      character*3 ubc
      integer tags(2),ibc,nbc,io,er
      integer ip,np
      integer ptags(2,npbc_max)
      real  pvecs(3,npbc_max)

      character*32 bcname
      character*1 bcnam1(32)
      equivalence(bcname,bcnam1)

      call blank (bcname,32)

      len = ltrunc(exoname,32)
      call chcopy(bcnam1,exoname,(len-3))
      call chcopy(bcnam1(len-3),'.bc' , 3)
      len = ltrunc(bcnam1,32)

      open(301,file=bcname,err=1010) ! if error, direct go to label 1010, return
      write(6,*) 'Setting boundary condition from ',bcname(:len),' file' 

      read(301,*,iostat=io) nbc
      if(io.ne.0) then
        write(6,*) bcname(:len),' file is empty '
        write(6,*) 'please set boundary condition in .usr file'
         return
      endif

       do ibc = 1,nbc
        call blank (ubc,3)
        read(301,*) tags(1),ubc
! not periodic boundary condition, direct set it up
          write(6,*) 'setting ',ubc,' to surface ',tags(1)
          do ihex = 1, num_elem
            do iface = 1,6
               if(bc(5,iface,ihex).eq.tags(1)) then
                 cbc(iface,ihex) = ubc
               endif
            enddo
          enddo
      enddo
      ip = 0
      read(301,*,iostat=io) nbc,ptol 
                       ! nbc is the number of pairs of periodic boundary condition
                       ! ptol is the tolerance used to search periodic boundary condition

      if(io.ne.0) then
         write(6,*) 'No periodic boundary condition set.'
         return
      endif

      if(nbc.gt.npbc_max) then
         write(6,*) 'ERROR: increase npbc_max to ',nbc
         write(6,*) 'and recompile exo2nek'
         return
      endif

      do ibc = 1,nbc
        ip = ip + 1
        read(301,*) ptags(1,ip),ptags(2,ip),pvecs(1,ip),pvecs(2,ip),pvecs(3,ip)
! ptags(1,ip) is surface 1 number sideset number
! ptags(2,ip) is surface 2 number sideset number
! pvecs(1,ip),pvecs(2,ip),pvecs(3,ip) is user defined projection vector
!  pvecs = 1 -> 2 
      enddo
      close(301)

      write(6,*) 'Setting periodic boundary condition' 
      np = ip
! np is the number of pairs of periodic bc now
      do ip = 1,np
! mapping surface ptags(1,ip) to surface ptags(2,ip)
        call setPeriodic(ptags(1,ip),pvecs(1,ip),ptol)
      enddo

1010  return
      end
!--------------------------------------------------------------------
      subroutine setPeriodic(ptags,pvec,ptol)
      use SIZE

      integer hex_face_node(4,6)
      data hex_face_node /1,3,21,19,3,9,27,21,7,9,27,25,1,7,25,19,1,7,9,3,19,21,27,25/

      integer parray(2,2,num_elem)
      real parea(2,num_elem)

      integer fnode(4),ifnode,ptags(2),ipe,nipe(2),nperror
      real pvec(3),ptol,fpxyz(3,2)
      real dist,distMax

      ipe = 0

! collect ihex,iface for surface ptags(1)
      do ihex = 1,num_elem
         do iface = 1,6
            if(bc(5,iface,ihex).eq.ptags(1)) then 
              ipe = ipe + 1
              parray(1,1,ipe) = ihex
              parray(2,1,ipe) = iface
            endif
         enddo
      enddo
      nipe(1) = ipe

! collect ihex,iface for surface ptags(2)
      ipe = 0
      do ihex = 1, num_elem
         do iface = 1,6
            if(bc(5,iface,ihex).eq.ptags(2)) then
                ipe = ipe + 1
                parray(1,2,ipe) = ihex
                parray(2,2,ipe) = iface
             endif
         enddo
      enddo
      nipe(2) = ipe

      write(6,*)'mapping surface',ptags(1),'with',nipe(1),'faces'
      write(6,*)'to surface',ptags(2),'with',nipe(2),'faces' 

      if(nipe(1).ne.nipe(2))  then
         write(6,*) 'EORROR, face numbers are not matching'
         write(6,*) 'periodic faces should have conformal mesh'
         return
      endif

! 1st loop, loop faces on surface 1
      do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
! get face center xyz
         call rzero(fpxyz(1,1),3)
         do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25			 
         enddo

! 2nd loop over surface 2
         distMax = 1000.0
         do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
! get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
               fnode(ifnode)=hex_face_node(ifnode,iface2)
               fpxyz(1,2) = fpxyz(1,2)+xm1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(2,2) = fpxyz(2,2)+ym1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(3,2) = fpxyz(3,2)+zm1(fnode(ifnode),1,1,ihex2)*0.25
              enddo
 
         dist = sqrt((fpxyz(1,2)-fpxyz(1,1)-pvec(1))**2+(fpxyz(2,2)-fpxyz(2,1)-pvec(2))**2+(fpxyz(3,2)-fpxyz(3,1)-pvec(3))**2)

               if (dist.lt.distMax) then 
                  distMax = dist
                  !write(6,*) distMax
                  if(distMax.le.ptol) then
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
                  bc(1,iface2,ihex2) = ihex*1.0
                  bc(2,iface2,ihex2) = iface*1.0
                  cbc(iface,ihex) = 'P  '
                  cbc(iface2,ihex2) = 'P  '
! for debug use only
!          write(6,*) ihex,iface,bc(1,iface,ihex),bc(2,iface,ihex)
!          write(6,*) ihex2,iface2,bc(1,iface2,ihex2),bc(2,iface2,ihex2)
!          write(6,*) fpxyz(1,1),fpxyz(2,1),fpxyz(3,1)
!          write(6,*) fpxyz(1,2),fpxyz(2,2),fpxyz(3,2)
!          write(6,*) dist,areadiff,parea(1,ipe),parea(2,ipe2)
                  endif
               endif
         enddo
      enddo

          nperror = 0

          write(6,*)'doing periodic check for surface',ptags(1)

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror, ' faces did not map'

          nperror = 0

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(iface.ne.iface3)) then

! for debug use only
!
!                write(6,*) 'ERROR,',ihex,iface,' map to ',ihex2,iface2
!                write(6,*) 'but,',ihex2,iface2,' map to ',ihex3,iface3
!
!		     do ifnode = 1,4
!         fnode(ifnode)=hex_face_node(ifnode,iface)
!      write(6,*)xm1(fnode(ifnode),1,1,ihex),ym1(fnode(ifnode),1,1,ihex),
!     & zm1(fnode(ifnode),1,1,ihex) 
!             enddo
!			 
!		     do ifnode = 1,4
!         fnode(ifnode)=hex_face_node(ifnode,iface2)
!      write(6,*)xm1(fnode(ifnode),1,1,ihex2),
!     & ym1(fnode(ifnode),1,1,ihex2),zm1(fnode(ifnode),1,1,ihex2) 
!             enddo
!
!		     do ifnode = 1,4
!         fnode(ifnode)=hex_face_node(ifnode,iface3)
!      write(6,*)xm1(fnode(ifnode),1,1,ihex3),
!     & ym1(fnode(ifnode),1,1,ihex3),zm1(fnode(ifnode),1,1,ihex3) 
!             enddo
                nperror = nperror + 1
             endif
          enddo		  

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong out of total ',nipe(1),' faces'
          endif

      return
      end
! -------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine setbc_new
      use SIZE
	 
      integer hex_face_node(4,6)
      data hex_face_node /1,3,21,19,3,9,27,21,7,9,27,25,1,7,25,19,1,7,9,3,19,21,27,25/
	 
      integer parray(2,2,num_elem)

      character*3 ubc
      integer tags(2),ibc,nbc,io
      integer ip,np,ipe,ipe2,nipe(2)
      integer ptags(2)
      integer fnode(4)
      real pvec(3)
      real fpxyz(3,2)
      real AB_v(3),AD_v(3),farea,product_v(3)
      real dist,distMax,ptol
 
! boundary condition summary
      write(6,*) '******************************************************'
      write(6,*) 'Boundary info summary'
      write(6,*) 'sideSet ID'
      do ibc= 1,bcNumber
      write(6,*) bcID(ibc)
      enddo
      write(6,*) '******************************************************'
 
 
      write(6,*) 'Enter number of periodic boundary surface pairs:'
      read (5,*) nbc
	  
      if(nbc.ne.0) then
      write(6,*) 'Enter search tolerance (1 as default):'
      read (5,*) ptol
      endif
	  
      if(nbc.le.0) return

      do ibc = 1,nbc 
        write(6,*) 'input surface 1 and  surface 2  sideSet ID'
        read (5,*) ptags(1),ptags(2)
        write(6,*) 'input translation vector (surface 1 -> surface 2)'
        read (5,*) pvec(1),pvec(2),pvec(3)

          ipe = 0
          do ihex = 1, num_elem
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(1)) then   
                ipe = ipe + 1
                parray(1,1,ipe) = ihex
                parray(2,1,ipe) = iface
               endif
            enddo
          enddo
          nipe(1) = ipe
	  
          ipe = 0
          do ihex = 1, num_elem
            do iface = 1,6
               if(bc(5,iface,ihex).eq.ptags(2)) then
                ipe = ipe + 1
                parray(1,2,ipe) = ihex
                parray(2,2,ipe) = iface
               endif
            enddo
          enddo
          nipe(2) = ipe

          if(nipe(1).ne.nipe(2))  then
            write(6,*)'mapping sideset ',ptags(1),'with',nipe(1),'faces'
            write(6,*)'to sideset ',ptags(2),'with',nipe(2),'faces'
            write(6,*) 'EORROR, face numbers are not matching'
          endif

! 1st loop, loop faces on surface 1
      do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
! get face center xyz
         call rzero(fpxyz(1,1),3)

         do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25
         enddo

! 2nd loop over surface 2
         distMax = ptol
         do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
! get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
               fnode(ifnode)=hex_face_node(ifnode,iface2)
               fpxyz(1,2) = fpxyz(1,2)+xm1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(2,2) = fpxyz(2,2)+ym1(fnode(ifnode),1,1,ihex2)*0.25
               fpxyz(3,2) = fpxyz(3,2)+zm1(fnode(ifnode),1,1,ihex2)*0.25
              enddo
 
       dist = sqrt((fpxyz(1,2)-fpxyz(1,1)-pvec(1))**2 &
      +(fpxyz(2,2)-fpxyz(2,1)-pvec(2))**2 &
      +(fpxyz(3,2)-fpxyz(3,1)-pvec(3))**2)

               if (dist.lt.distMax) then 
                  distMax = dist
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
               endif
         enddo
      enddo

! change. only assign periodic face at the end of loop.
! this will assign the closest face for periodicity.

     do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
         ihex2 = int(bc(1,iface,ihex))
         iface2 = int(bc(2,iface,ihex)) 
         bc(1,iface2,ihex2) = ihex*1.0
         bc(2,iface2,ihex2) = iface*1.0
         cbc(iface,ihex) = 'P  '
         cbc(iface2,ihex2) = 'P  '
     enddo

        nperror = 0
   
          write(6,*)'doing periodic check for surface',ptags(1)

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror, ' faces did not for periodicity'

          nperror = 0

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(iface.ne.iface3)) then
               nperror = nperror + 1
             endif
          enddo

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong out of total ',nipe(1),' faces'
          endif
  
          write(6,*)'doing periodic check for surface',ptags(2)

          nperror = 0

          do ipe = 1,nipe(2)
             ihex = parray(1,2,ipe)
             iface = parray(2,2,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) write(6,*) 'ERROR,',nperror, ' faces did not map for periodicity'

          nperror = 0
          do ipe = 1,nipe(2)
             ihex = parray(1,2,ipe)
             iface = parray(2,2,ipe)
             ihex2 = int(bc(1,iface,ihex))
             iface2 = int(bc(2,iface,ihex))
             ihex3 = int(bc(1,iface2,ihex2))
             iface3 = int(bc(2,iface2,ihex2))
             if ((ihex.ne.ihex3).or.(iface.ne.iface3)) then
                nperror = nperror + 1
             endif
          enddo		  

          if (nperror.gt.0) then
          write(6,*) 'ERROR,',nperror,'faces are wrong out of total ',nipe(2),' faces'
          endif
 
      enddo

      write(6,*) '******************************************************'
      write(6,*) 'Please set boundary conditions to all non-periodic boundaries'
      write(6,*) 'in .usr file usrdat2() subroutine'
      write(6,*) '******************************************************'

      return
      end
!-----------------------------------------------------------------------
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

!      integer,allocatable,dimension(:) :: idblk
!      integer,allocatable,dimension(:) :: num_nodes_per_elem
!      integer,allocatable,dimension(:) :: num_attr    !not used

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
      subroutine offset_sideset(iexo1)
! offset sideset number by iexo for element btween eacc and eacc_old
      use SIZE
      integer iexo1

      do iel= eacc_old+1,eacc
       do jfc =1,6
         if (bc(5,jfc,iel).gt.0) then
           bc(5,jfc,iel) = bc(5,jfc,iel)  + iexo1*100
         endif		 
       enddo
      enddo

      return
      end
!--------------------------------------------------------------------
      subroutine gather_bc_info
      use SIZE

      integer ibc,ibc2,nbc
      logical newbc
	  
      integer bcID2
      allocate (bcID (100)) ! assuming there is no more than 100 sidesets in total
	  
      ibc = 0
	  
      ! gather all boundary information
	  
      do iel= 1,num_elem
       do jfc =1,6

         if(ibc.eq.0) then
         if (bc(5,jfc,iel).ne.0) then
          ibc = ibc +1
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc		  
         endif
         endif 
		 
         if(ibc.gt.0) then
         if (bc(5,jfc,iel).ne.0) then
          newbc = .TRUE.
          do ibc2 =1,nbc
           if (bc(5,jfc,iel).eq.bcID(ibc2)) then
            newbc = .FALSE.
           endif
          enddo 

          if(newbc) then		
          ibc = ibc +1		  
          bcID(ibc) = bc(5,jfc,iel)
          nbc = ibc
          endif

         endif
         endif 

       enddo
      enddo
	  
      ! sorting bcID array to ascend order 
	  
       do ibc2 =1,nbc-1
        do ibc =1,nbc-1
         if( bcID(ibc).gt.bcID(ibc+1)) then
         bcID2 = bcID(ibc)
         bcID(ibc) = bcID(ibc+1)
         bcID(ibc+1) = bcID2
         endif
        enddo
       enddo
      bcNumber = nbc

      return
      end
	  
!--------------------------------------------------------------------
      subroutine scale_mesh
      use SIZE
      real xx,yy,zz,ss
      integer ntot

      write(6,*) "please input scaling factor (1 for no scale):"
      read(5,*) ss
	  
      ntot = num_elem*3*3*3

      do i = 1,ntot
        xx = xm1(i,1,1,1)        
        yy = ym1(i,1,1,1)
        zz = zm1(i,1,1,1)
  
        xm1(i,1,1,1) = xx*ss
        ym1(i,1,1,1) = yy*ss
        zm1(i,1,1,1) = zz*ss
      enddo
	  
      call checkXYZ_min_max()

      return
      end
!--------------------------------------------------------------------
      subroutine right_hand_check
! check if there is non-right hand elements (3D)
! because if mesh is from ICEM, and mirror operation is made in ICEM,
! the exported exo file will contain non-right hand elements.
! this subroutine will:
! 1. check right-hand
! 2. fix if not
      use SIZE
      logical ifnonrighthand

      do iel=1,num_elem
         call check_if_non_right(ifnonrighthand,iel)
         if (ifnonrighthand) call fix_if_non_right(iel)		 
      enddo

      return
      end
!--------------------------------------------------------------------
      subroutine check_if_non_right(ifnonrighthand,iel)
      use SIZE
      logical ifnonrighthand
      integer iel
      integer hex8_to_hex27_vertex(8)
      data hex8_to_hex27_vertex /1,3,7,9,19,21,25,27/ 
      real hex8_vertex(3,8),vec12(3),vec14(3),vec15(3)
      real vec1(3),AA,dot_prod

      do iver = 1,8
       hex8_vertex(1,iver) = xm1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(2,iver) = ym1(hex8_to_hex27_vertex(iver),1,1,iel)
       hex8_vertex(3,iver) = zm1(hex8_to_hex27_vertex(iver),1,1,iel)       
      enddo
	  
      vec12(1) = hex8_vertex(1,2) - hex8_vertex(1,1)
      vec12(2) = hex8_vertex(2,2) - hex8_vertex(2,1)
      vec12(3) = hex8_vertex(3,2) - hex8_vertex(3,1)

      vec14(1) = hex8_vertex(1,4) - hex8_vertex(1,1)
      vec14(2) = hex8_vertex(2,4) - hex8_vertex(2,1)
      vec14(3) = hex8_vertex(3,4) - hex8_vertex(3,1)

      vec15(1) = hex8_vertex(1,5) - hex8_vertex(1,1)
      vec15(2) = hex8_vertex(2,5) - hex8_vertex(2,1)
      vec15(3) = hex8_vertex(3,5) - hex8_vertex(3,1)

      call cross_product(vec12,vec14,vec1,AA) 
      dot_prod = vec1(1)*vec15(1) + vec1(2)*vec15(2) + vec1(3)*vec15(3)
  
      if(dot_prod.gt.0.0) then
       ifnonrighthand = .FALSE.
      else
       ifnonrighthand = .TRUE.
       !write(6,*) 'non-right hand element detected'
      endif  

      return
      end
!-----------------------------------------------------------------------
      subroutine cross_product(AB_v,AC_v,prod_v,area)
! calculate cross product of two vectors

      real*8 AB_v(3),AC_v(3),prod_v(3)
      real*8 area

      prod_v(1) = AB_v(2)*AC_v(3)-AB_v(3)*AC_v(2)
      prod_v(2) = AB_v(3)*AC_v(1)-AB_v(1)*AC_v(3)
      prod_v(3) = AB_v(1)*AC_v(2)-AB_v(2)*AC_v(1)

      area = prod_v(1)**2 + prod_v(2)**2+prod_v(3)**2
      area = sqrt(area)
      area = area/2.0
 
      return 
      end
!------------------------------------------------------------------------------------------
!--------------------------------------------------------------------
      subroutine fix_if_non_right(iel)
! fix non-right hand element
! 1 <->19, 3 <-> 21, 7<->25, 9 <->27 
      use SIZE
      integer iel,iver
      real*8 xm2(27),ym2(27),zm2(27)
      character(3) cbc5,cbc6
      real bc55,bc56

! swap vertex
      do iver = 1,27
       xm2(iver) = xm1(iver,1,1,iel)
       ym2(iver) = ym1(iver,1,1,iel)
       zm2(iver) = zm1(iver,1,1,iel)
      enddo

      do iver = 1,9
       xm1(iver,1,1,iel) = xm2(iver+18)
       ym1(iver,1,1,iel) = ym2(iver+18)
       zm1(iver,1,1,iel) = zm2(iver+18)
      enddo
     
      do iver = 19,27
       xm1(iver,1,1,iel) = xm2(iver-18)
       ym1(iver,1,1,iel) = ym2(iver-18)
       zm1(iver,1,1,iel) = zm2(iver-18)
      enddo

! swap face 5 <-> 6

      cbc5 = cbc(5,iel)
      bc55 = bc(5,5,iel)
 
      cbc6 = cbc(6,iel)
      bc56 = bc(5,6,iel)

      cbc(5,iel)   = cbc6
      bc (5,5,iel) = bc56

      cbc(6,iel)   = cbc5
      bc (5,6,iel) = bc55 

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
            
      num_elem = etot
      num_dim = 3

!  Write the header
      call blank   (hdr,80)    
      write(hdr,1) num_elem, num_dim, eftot
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

      !rbc = num_sides_tot
      !call byte_write (rbc,2,ierr)

      nbc   = 0
      nface = 2*num_dim
      do iel=1,eftot
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo
      rbc = nbc
      call byte_write (rbc,2, ierr)
      write(6,*) 'velocity boundary faces: ',nbc


      do iel = 1,eftot
        do ifc = 1,2*num_dim
          ch3 = cbc(ifc,iel)
          if (ch3.ne.'   ') then
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (eftot.ge.1000000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16,ierr)
          endif
        enddo
      enddo

     if(num_elem.ne.eftot) then

! writing thermal bc to all elements
      nbc   = 0
      nface = 2*num_dim
      do iel=1,num_elem
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo
      rbc = nbc
      call byte_write (rbc,2, ierr)
      write(6,*) 'thermal boundary faces: ',nbc

      do iel = 1,num_elem
        do ifc = 1,2*num_dim
          ch3 = cbc(ifc,iel)
          if (ch3.ne.'   ') then
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

      endif

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
!-----------------------------------------------------------------------
      subroutine rzero_int(a,n)
      integer A(1)
      DO 100 I = 1, N
 100     A(I) = 0
      return
      END
