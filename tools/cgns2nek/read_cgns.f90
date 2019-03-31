SUBROUTINE read_cgns
  !
  !  READ  X, Y, Z GRID POINTS AND BOUNDARY CONDITIONS FROM CGNS FILE
  !
  USE CGNS
  USE module_global
  !
  IMPLICIT NONE
  !
  INTEGER(iprec) :: icelldim,iphysdim,b,z,phys_dim,s,max_vert,num_bc
  INTEGER(iprec) :: ier,c,i,j,k,cell_dim,elem_no,ngrids,Dset,l11,n,l1,l2,l3
  INTEGER(iprec) :: ib,ie,ii,i1,i2,i3,i4,i5,i6,i7,i8,is,js,nb,jj
  INTEGER(iprec) :: gridlocation,elementconn,bc,DirichletFlag, NeumannFlag
  INTEGER(iprec) :: ifirstnode,nconns,connect_type,data_type,inc_2,axes
  INTEGER(iprec) :: fam,nfamilies,ngeo,nfambc,nnames,int_dir,perio,trans
  INTEGER(iprec) :: thread, id
  INTEGER(iprec) :: index_coord,nzone,zone_type,index_section,parent
  INTEGER(iprec) :: NormalIndex(3),NormalDataType, ndataset,ptset_type,max_bc
  INTEGER(iprec) :: omp_get_thread_num,omp_get_num_threads,sum_f,nmax
  INTEGER(iprec) :: iel,iv,nv,nel,vid,flag
  INTEGER(iprec),ALLOCATABLE :: npe(:),per_fam(:),j_used(:),sum_face(:)

  INTEGER(cgsize_t),PARAMETER :: min_v= 1
  INTEGER(cgsize_t) :: max_v, elementdatas,nvert,num_vert
  INTEGER(cgsize_t) :: nelem_s,nelem_e,normallistsize,ver_id(8),bnd_vert(4)
  INTEGER(cgsize_t),ALLOCATABLE :: npnts(:),tmp_e(:),tmp_pa(:), isize(:,:),per_val(:,:)
  INTEGER(cgsize_t),ALLOCATABLE :: volnodes(:,:)
  INTEGER(cgsize_t),ALLOCATABLE :: bound_elem(:),bc_vert(:)
  !
  REAL(high),PARAMETER :: eps=1.d-6
  REAL(high) :: diff,diffx,diffy, diffz,grad,displace,t1,t2,ta(2),time
  !
  CHARACTER(32),ALLOCATABLE :: familyname(:)
  CHARACTER(1) :: dir
  CHARACTER(32) :: degree_
  !
  perio = 0

  CALL cg_get_file_type_f(ifile,s,ier)
  IF (s /= CG_FILE_ADF) THEN
     WRITE(6,*)
     WRITE(6,*) ' ADF file format required!'
     CALL cg_error_exit_f
  ENDIF
 
  CALL cg_nbases_f(ifile,nbases,ier)
  IF (nbases /= 1) THEN
     WRITE(6,'('' Number of base nodes = '',i5)') nbases
     WRITE(6,*)
     WRITE(6,*) ' No support for more than one CGNS base nodes!'
     CALL cg_error_exit_f
  ENDIF
  !
  DO b=1,nbases
     !
     CALL cg_nzones_f (ifile,b,nzones,ier)
     CALL cg_base_read_f(ifile, b, basename, cell_dim, phys_dim, ier)
  !
  IF (nzones == 0) THEN
     WRITE(6,'('' Number of zones = '',i5)') nzones
     WRITE(6,*)
     WRITE(6,*) ' No support for more than one zone!'
     CALL cg_error_exit_f
  ENDIF
  !
  IF (cell_dim /= 3 .OR. phys_dim /= 3 ) THEN
     WRITE(6,*)
     WRITE(6,*) ' 3-dim mesh required!'
     CALL cg_error_exit_f
  ENDIF
     !
     ALLOCATE (zonetype(nzones))
     ALLOCATE (index_dim(nzones))
     ALLOCATE (isize(cell_dim, phys_dim))
     ALLOCATE (zonename(nzones))
     !
     isize = 0
     !
     !  read zone
     !
     DO z=1,nzones

        CALL cg_zone_type_f (ifile,b,z,zonetype(z),ier)
        CALL cg_index_dim_f (ifile,b,z,index_dim(z),ier)
        CALL cg_zone_read_f (ifile,b,z,zonename(z),isize,ier)

        CALL cg_ngrids_f(ifile,b,z,ngrids,ier)
        CALL cg_ncoords_f(ifile,b,z,ncoords,ier)

        ALLOCATE (coord_name(ncoords))
        !
        !  read sections
        !
        CALL cg_nsections_f(ifile,b,z,nsections,ier)

        ALLOCATE (sec_name(nsections))
        ALLOCATE (sec_type(nsections))
        ALLOCATE (nelem_start(nsections))
        ALLOCATE (nelem_end(nsections))
        ALLOCATE (nbdyelem(nsections))
        ALLOCATE (sectype(nsections))
        ALLOCATE (npe(nsections))
        ALLOCATE (elementdatasize(nsections))

        nelem_start = 0
        nelem_end = 0

        max_elem = 0
        DO s=1,nsections
           CALL cg_ElementDataSize_f(ifile,b,z,s,elementdatas,ier)
           elementdatasize(s)= elementdatas
        ENDDO
        !
        num_vert = MAXVAL(elementdatasize)
        !
        ALLOCATE (elements(nsections,num_vert))
        ALLOCATE (tmp_e(num_vert))
        ALLOCATE (tmp_Pa(num_vert))
        !
        num_bc = 0
        !
        WRITE(6,*) ' Scanning sections ...'
        DO s=1,nsections
           CALL cg_section_read_f(ifile,b,z,s,sec_name(s),sectype(s),nelem_s,nelem_e &
                ,nbdyelem(s),parent,ier)

           WRITE(6,*) s, ' Element type ', ElementTypeName(sectype(s))
           IF (sectype(s) /= HEXA_8  .AND. &
               sectype(s) /= HEXA_20 .AND. &
               sectype(s) /= HEXA_27 .AND. &
               sectype(s) /= QUAD_4  .AND. &
               sectype(s) /= QUAD_8  .AND. &
               sectype(s) /= QUAD_9) THEN 
               WRITE(6,*) ' Skipping unsupported element type! ', ElementTypeName(sectype(s))
               CYCLE
           ENDIF

           CALL cg_npe_f(sectype(s),npe(s),ier)

           CALL cg_elements_read_f(ifile,b,z,s,tmp_e,tmp_Pa,ier)

           elements(s,1:elementdatasize(s)) = tmp_e(1:elementdatasize(s))
           nelem_start(s) = nelem_s
           nelem_end(s)   = nelem_e
           !
           IF (sectype(s) == QUAD_4 .OR. sectype(s) == QUAD_8 .OR. sectype(s) == QUAD_9) THEN
              num_bc = num_bc +1
           ENDIF
        ENDDO
        !
        DEALLOCATE (tmp_e)
        DEALLOCATE (tmp_Pa)
        !
        DO c=1,ncoords
           CALL cg_coord_info_f(ifile,b,z,c,data_type,coord_name(c),ier)
        ENDDO

        max_vert = isize(1,z)
        max_v = isize(1,z)

        ALLOCATE (xv(max_vert))
        ALLOCATE (yv(max_vert))
        ALLOCATE (zv(max_vert))

        IF (data_type == RealSingle ) THEN
           ALLOCATE (xvsingle(max_vert))
           ALLOCATE (yvsingle(max_vert))
           ALLOCATE (zvsingle(max_vert))
           CALL cg_coord_read_f(ifile,b,z,coord_name(1),data_type,min_v,max_v,xvsingle,ier)
           CALL cg_coord_read_f(ifile,b,z,coord_name(2),data_type,min_v,max_v,yvsingle,ier)
           CALL cg_coord_read_f(ifile,b,z,coord_name(3),data_type,min_v,max_v,zvsingle,ier)
           xv = DBLE(xvsingle)
           yv = DBLE(yvsingle)
           zv = DBLE(zvsingle)
           DEALLOCATE (xvsingle)
           DEALLOCATE (yvsingle)
           DEALLOCATE (zvsingle)
        ELSE
           CALL cg_coord_read_f(ifile,b,z,coord_name(1),data_type,min_v,max_v,xv,ier)
           CALL cg_coord_read_f(ifile,b,z,coord_name(2),data_type,min_v,max_v,yv,ier)
           CALL cg_coord_read_f(ifile,b,z,coord_name(3),data_type,min_v,max_v,zv,ier)
        ENDIF
        !
        CALL cg_nbocos_f(ifile,b,z,nbocos,ier)
        !
        ALLOCATE (boconame(nbocos))
        ALLOCATE (bocotype(nbocos))
        ALLOCATE (normallist(nbocos))
        ALLOCATE (npnts(nbocos))
        ALLOCATE (location(nbocos))
        ALLOCATE (DatasetName(nbocos))
        ALLOCATE (BCType(nbocos))
        ALLOCATE (per_fam(nbocos))
        !
        DO bc=1,nbocos
           !
           CALL cg_boco_info_f(ifile, b, z, bc, boconame(bc), bocotype(bc), ptset_type, npnts(bc), &
                NormalIndex, NormalListSize, NormalDataType, ndataset,ier)
           !
           CALL cg_boco_gridlocation_read_f(ifile, b, z, bc, location(bc), ier)
           !
           CALL cg_nfamilies_f(ifile, b, nfamilies, ier)
           !
        ENDDO
        !
        max_bc = MAXVAL(npnts)
        !
        CALL cg_nconns_f(ifile, B, Z, nconns, ier)
        !
        ! Determine the section number for the elements and boundary faces
        nb = 0
        !
        DO s=1,nsections
           IF(sectype(s) == HEXA_8 .OR. sectype(s) == HEXA_20 .OR. sectype(s) == HEXA_27) THEN
              max_elem = nelem_end(s)-nelem_start(s)+1
              sec = s
           ENDIF
           IF(sectype(s) == QUAD_4 .OR. sectype(s) == QUAD_8 .OR. sectype(s) == QUAD_9) THEN
              nb = nb + nelem_end(s)-nelem_start(s)+1
           ENDIF
        ENDDO
        !
        ! Determine the parent connectivity list
        ! First calculate the face mid points of all elements
        !
        IF (sectype(sec) == HEXA_8) THEN
           inc= 8
        ELSEIF (sectype(sec) == HEXA_20) THEN
           inc= 20
        ELSEIF (sectype(sec) == HEXA_27) THEN
           inc= 27
        ENDIF
        !
        nel = nelem_end(sec)-nelem_start(sec)+1
        !
        nv  = inc ! only first 8 vertices taken
        !
        nvert = MAXVAL(elements(sec,:))
        !
        ALLOCATE(bound_elem(nel*8))
        ALLOCATE(volnodes(0:8,nel))
        !
        bound_elem = 0
        !
        !  Filtering the inner elements out
        !
        n = 0
        DO iel = 1,num_vert,inc
           !
           n = n +1
           ver_id(1:8) = elements(sec,iel:iel+7)
           volnodes(1:8,n) = ver_id(1:8)
           bound_elem(ver_id(1:8)) = bound_elem(ver_id(1:8))+1
           !
        ENDDO
        !
        DO iel = 1,nel
           vid = SUM(bound_elem(volnodes(1:8,iel)))
           IF (vid < 64 ) THEN
              volnodes(0,iel) = 1  
           ENDIF
        ENDDO
        !
        ! Check if infinite thin walls exists
        !
        ALLOCATE(bc_vert(max_vert))
        !
        IF (sectype(2) == QUAD_4) THEN
           inc_2= 4
        ELSEIF (sectype(2) == QUAD_8) THEN
           inc_2= 8
        ELSEIF (sectype(2) == QUAD_9) THEN
           inc_2= 9
        ENDIF
        !
        bc_vert = 0
        ver_id = 0
        flag = 0
        !
        DO bc = 1,nbocos 
           DO iel = 1,elementdatasize(1+bc)
              !
              vid = elements(1+bc,iel)
              bc_vert(vid) = bc_vert(vid)+1
              !
           ENDDO
        ENDDO
        !
        !
        DO iel = 1,max_vert
           IF (bc_vert(iel) > 4 ) flag = 1 
        ENDDO
        !
        IF (flag == 1) THEN
           !
           ! Find boundary families with common vertices
           !
           DO iel = 1,max_elem
              IF (volnodes(0,iel) == 1) CYCLE
              !
              ver_id(1:8) = volnodes(1:8,iel)
              !
              bnd_loop: DO bc = 1,nbocos 
                 DO ii = 1,elementdatasize(1+bc),inc_2
                    n = 0
                    bnd_vert(1:4) = elements(1+bc,ii:ii+3)
                    DO i = 1,4 
                       DO j = 1,8 
                          IF (bnd_vert(i) == ver_id(j) ) THEN 
                             n= n+1
                          ENDIF
                       ENDDO
                    ENDDO
                    IF (n ==4) THEN
                       volnodes(0,iel) = 1  
                       EXIT bnd_loop
                    ENDIF
                 ENDDO
              ENDDO bnd_loop
              !
           ENDDO
           !
        ENDIF
        !
        DEALLOCATE (bc_vert,bound_elem)
        !
        ALLOCATE (xf(max_elem,6))
        ALLOCATE (yf(max_elem,6))
        ALLOCATE (zf(max_elem,6))
        !
        ii = 1
        ie = 1
        !
        DO i= nelem_start(sec),nelem_end(sec)
           !
           i1= elements(sec,ii)
           i2= elements(sec,ii+1)
           i3= elements(sec,ii+2)
           i4= elements(sec,ii+3)
           i5= elements(sec,ii+4)
           i6= elements(sec,ii+5)
           i7= elements(sec,ii+6)
           i8= elements(sec,ii+7)
           !
           !       Face 1    
           xf(ie,1) = 0.25*(xv(i1)+xv(i2)+xv(i5)+xv(i6))
           yf(ie,1) = 0.25*(yv(i1)+yv(i2)+yv(i5)+yv(i6))
           zf(ie,1) = 0.25*(zv(i1)+zv(i2)+zv(i5)+zv(i6))
           !       Face 2    
           xf(ie,2) = 0.25*(xv(i2)+xv(i3)+xv(i6)+xv(i7))
           yf(ie,2) = 0.25*(yv(i2)+yv(i3)+yv(i6)+yv(i7))
           zf(ie,2) = 0.25*(zv(i2)+zv(i3)+zv(i6)+zv(i7))
           !       Face 3
           xf(ie,3) = 0.25*(xv(i3)+xv(i4)+xv(i7)+xv(i8))
           yf(ie,3) = 0.25*(yv(i3)+yv(i4)+yv(i7)+yv(i8))
           zf(ie,3) = 0.25*(zv(i3)+zv(i4)+zv(i7)+zv(i8))
           !       Face 4    
           xf(ie,4) = 0.25*(xv(i1)+xv(i4)+xv(i5)+xv(i8))
           yf(ie,4) = 0.25*(yv(i1)+yv(i4)+yv(i5)+yv(i8))
           zf(ie,4) = 0.25*(zv(i1)+zv(i4)+zv(i5)+zv(i8))
           !       Face 5
           xf(ie,5) = 0.25*(xv(i1)+xv(i2)+xv(i3)+xv(i4))
           yf(ie,5) = 0.25*(yv(i1)+yv(i2)+yv(i3)+yv(i4))
           zf(ie,5) = 0.25*(zv(i1)+zv(i2)+zv(i3)+zv(i4))
           !       Face 6
           xf(ie,6) = 0.25*(xv(i5)+xv(i8)+xv(i6)+xv(i7))
           yf(ie,6) = 0.25*(yv(i5)+yv(i8)+yv(i6)+yv(i7))
           zf(ie,6) = 0.25*(zv(i5)+zv(i8)+zv(i6)+zv(i7))
           !
           ii = ii +inc
           ie = ie + 1
        ENDDO
        !
        ! Compare face coordinates to determine parent elements
        !
        ALLOCATE (vel_type(max_elem,6))
        ALLOCATE (vel_val(max_elem,6,5))
        !
        vel_type = 'E  '
        vel_val = 0.
        !
        !
        ALLOCATE (xb(nb))
        ALLOCATE (yb(nb))
        ALLOCATE (zb(nb))
        ALLOCATE (per_val(nb,2))
        ALLOCATE (bc_fam(nb))
        !
        ib= 0
        !
        DO ii= 1,nsections
           !
           jj = 1
           !
           IF (sectype(ii) == QUAD_4) THEN
              inc_2= 4
           ELSEIF (sectype(ii) == QUAD_8) THEN
              inc_2= 8
           ELSEIF (sectype(ii) == QUAD_9) THEN
              inc_2= 9
           ELSE
              CYCLE
           ENDIF
           !
           DO j= nelem_start(ii),nelem_end(ii)
              i1= elements(ii,jj)
              i2= elements(ii,jj+1)
              i3= elements(ii,jj+2)
              i4= elements(ii,jj+3)
              ib = ib+1
              !       bc Face 
              xb(ib) = 0.25*(xv(i1)+xv(i2)+xv(i3)+xv(i4))
              yb(ib) = 0.25*(yv(i1)+yv(i2)+yv(i3)+yv(i4))
              zb(ib) = 0.25*(zv(i1)+zv(i2)+zv(i3)+zv(i4))
              bc_fam(ib) = ii-1 
              jj = jj+inc_2
           ENDDO
        ENDDO
        !
        ib = 0
        !
        ! Check if one to one relation exist of nbocos and sectype
        !
        IF (num_bc /= nbocos ) THEN
           !
           DO i= 1,nbocos
              !
              DO j= 1,npnts(i)
                 ib = ib+1
                 bc_fam(ib) = i 
              ENDDO
           ENDDO
           !
        ENDIF
        !
        ! Compare face coordinates to determine all boundary conditions
        !
        ALLOCATE (j_used(nb))
        j_used=0
        !
        WRITE(6,*) 
        WRITE(6,*) ' Number of volume elements = ', max_elem
        WRITE(6,*) ' Number of boundary face elements = ',nb
        WRITE(6,*)
        !
        ALLOCATE (sum_face(nb))
        sum_face = 0
        !
        !$OMP parallel 
        !$OMP do
        DO i= 1,max_elem
           IF (volnodes(0,i) == 0) CYCLE

           DO is= 1,6
              bc_loop: DO j= 1,nb
                 IF (j_used(j) == 1) CYCLE 
                 diff = ABS(xf(i,is)-xb(j))+ ABS(yf(i,is)-yb(j))+ ABS(zf(i,is)-zb(j))
                 IF (diff < eps) THEN
                    vel_type(i,is) = 'CGN'
                    j_used(j) = 1
                    per_val (j,1) = DBLE(i)
                    per_val (j,2) = DBLE(is)
                    vel_val (i,is,1:4) = 0.
                    vel_val (i,is,5) = bc_fam(j)
                    sum_face(j) = 1
                    EXIT bc_loop
                 ENDIF
              ENDDO bc_loop
              !
           ENDDO
        ENDDO
        !$OMP end do
        !$OMP end parallel
        !
        sum_f = SUM(sum_face)
        IF (sum_f < nb) THEN
           WRITE(6,*) 
           WRITE(6,*) '   NOT enough BC found !!!  BC = ',sum_f,' Max BC =', nb
           CALL cgnstonek_exit('       ')
        ENDIF
        !
        DO bc=1,nbocos
           l1 = INDEX (boconame(bc), ' ') - 1
           WRITE(6,100) 'ID = ', bc, 'name: ', boconame(bc)
 100       FORMAT (2x,a6,i3,2x,a6,a18)
        ENDDO
        !
        ier = 0
        WRITE(6,*)
 150    WRITE(6,'('' Specify number of periodic boundary face pairs: '')', ADVANCE = "NO")
        READ(5,'(i1)', IOSTAT = ier) perio
        IF (ier /= 0) GOTO 150
        !
        DO n=1,perio
           !
           diffx=0.
           diffy=0.
           diffz=0.
           sum_face = 0
           !
           WRITE(6,'('' Type in periodic pair ID: '')', ADVANCE = "NO")
           READ(5,'(a)') dir 
           READ(dir(1:1),*) per_fam(n)
           !
           !WRITE(6,'('' Type in if : translational=0 or rotational=1 '')')
           !READ(5,'(a)') dir 
           !READ(dir(1:1),*) int_dir
           int_dir = 0
           !
           IF (int_dir == 1) THEN
              !
 200          WRITE(6,'('' Type in rotation axis (x=1/y=2/z=3): '')', ADVANCE = "NO")
              READ(5,'(a)') dir 
              READ(dir(1:1),*) axes
              IF (trans < 1 .OR. trans > 3) GOTO 200
              !
              WRITE(6,'('' Type in rotation angle in degrees: '')', ADVANCE = "NO")
              READ(5,'(a)') degree_
              l11 = INDEX (degree_, ' ') - 1
              READ(degree_(1:l11),*) grad
           ELSE
              !
 201           WRITE(6,'('' Type in translational direction (x=1/y=2/z=3): '')', ADVANCE = "NO")
              READ(5,'(a)') dir 
              READ(dir(1:1),*) trans
              IF (trans < 1 .OR. trans > 3) GOTO 201
              !
              WRITE(6,'('' Type in translational displacement '')', ADVANCE = "NO")
              READ(5,'(a)') degree_
              l11 = INDEX (degree_, ' ') - 1
              READ(degree_(1:l11),*) displace
              !
              IF (trans == 1) THEN
                 diffx = displace
              ELSEIF (trans == 2) THEN
                 diffy = displace
              ELSEIF (trans == 3) THEN
                 diffz = displace
              ENDIF
              !
           ENDIF
           !
           !$OMP parallel do
           !
           DO i= 1,nb 
              per_loop: DO j= 1,nb
                 IF (i == j) CYCLE
                 diff=ABS(xb(i)+diffx-xb(j))+ ABS(yb(i)+diffy-yb(j))+ ABS(zb(i)+diffz-zb(j))
                 IF (diff < eps ) THEN
                    ii = per_val(i,1)
                    is = per_val(i,2)
                    jj = per_val(j,1)
                    js = per_val(j,2)
                    vel_val (ii,is,1) = per_val(j,1)
                    vel_val (ii,is,2) = per_val(j,2)
                    vel_val (jj,js,1) = per_val(i,1)
                    vel_val (jj,js,2) = per_val(i,2)
                    vel_val (ii,is,3:4) = 0.
                    vel_val (jj,js,3:4) = 0.
                    vel_type(ii,is) = 'P  '
                    vel_type(jj,js) = 'P  '
                    !
                    sum_face(j) = 2
                    EXIT per_loop
                 ENDIF
              ENDDO per_loop
              !
           ENDDO
           !
           !$OMP end parallel do
           !
           sum_f = SUM(sum_face)
           !
           ii = per_fam(n)+1
           IF (sum_f < nelem_end(ii)-nelem_start(ii)+1) THEN
              WRITE(6,*) 
              WRITE(6,*) '      NOT enough periodic elements found !!! in BC = ',ii
              CALL cgnstonek_exit('       ')
           ENDIF
           !
        ENDDO
           !
        DEALLOCATE (xb)
        DEALLOCATE (yb)
        DEALLOCATE (zb)
        !
     ENDDO
     !
  ENDDO
  !
  !  close CGNS file
  !
  CALL cg_close_f(ifile,ier)
  !
  !
  RETURN
  !
END SUBROUTINE read_cgns
