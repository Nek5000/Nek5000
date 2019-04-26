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
  INTEGER(iprec) :: ib,ie,ii,i1,i2,i3,i4,i5,i6,i7,i8,is,js,nb,jj,nn
  INTEGER(iprec) :: gridlocation,elementconn,bc,DirichletFlag, NeumannFlag
  INTEGER(iprec) :: ifirstnode,nconns,connect_type,data_type,inc_2,axes
  INTEGER(iprec) :: fam,nfamilies,ngeo,nfambc,nnames,int_dir,perio,trans
  INTEGER(iprec) :: thread, id,sec_cnt
  INTEGER(iprec) :: index_coord,nzone,zone_type,index_section,parent
  INTEGER(iprec) :: NormalIndex(3),NormalDataType, ndataset,ptset_type,max_bc
  INTEGER(iprec) :: omp_get_thread_num,omp_get_num_threads,sum_f,nmax
  INTEGER(iprec) :: iel,iv,nv,nel,vid,flag,bc_par(2)
  INTEGER(iprec),ALLOCATABLE :: npe(:),per_fam(:),j_used(:),sum_face(:)

  INTEGER(cgsize_t),PARAMETER :: min_v= 1
  INTEGER(cgsize_t) :: max_v, elementdatas,nvert,num_vert,vrt(4)
  INTEGER(cgsize_t) :: nelem_s,nelem_e,normallistsize,ver_id(8),bnd_vert(4)
  INTEGER(cgsize_t),ALLOCATABLE :: npnts(:),tmp_e(:),tmp_pa(:), isize(:,:),per_val(:,:)
  INTEGER(cgsize_t),ALLOCATABLE :: volnodes(:,:),num_point(:,:)
  INTEGER(cgsize_t),ALLOCATABLE :: bound_elem(:),bc_vert(:)
  !
  REAL(high),PARAMETER :: eps=1.d-5,PI = 3.14159265358979
  REAL(high) :: diff,diffx,diffy, diffz,grad,displace,t1,t2,ta(2),time
  REAL(high) :: n1(4,3),n2(4,3),nmean(2,3),ax,ay,az,bx,by,bz,dist(2,3),tri(5,3)
  REAL(high) :: rot_angle,rot_axe,uabs,vabs,udotv,cosn,beta,radius,rot_dir
  !
  CHARACTER(32),ALLOCATABLE :: familyname(:)
  CHARACTER(100) :: dir
  CHARACTER(32) :: degree_
  !

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
        ALLOCATE (num_point(nsections,2))

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
        nbocos = num_bc 
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
           IF (bc_vert(iel) > 6 ) flag = 1 
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
        DEALLOCATE (bc_vert)
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
        sec_cnt = 0
        !
        DO ii= 1,nsections
           !
           jj = 1
           !
           IF (sectype(ii) == QUAD_4) THEN
              inc_2= 4
              sec_cnt = sec_cnt +1
           ELSEIF (sectype(ii) == QUAD_8) THEN
              inc_2= 8
              sec_cnt = sec_cnt +1
           ELSEIF (sectype(ii) == QUAD_9) THEN
              inc_2= 9
              sec_cnt = sec_cnt +1
           ELSE
              CYCLE
           ENDIF
           !
           num_point(sec_cnt,1) = ib+1
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
           num_point(sec_cnt,2) = ib
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
        ALLOCATE (j_used(nb))
        ALLOCATE (sum_face(nb))
        sum_face = 0
        j_used=0
        !
        WRITE(6,*) 
        WRITE(6,*) ' Number of volume elements = ', max_elem
        WRITE(6,*) ' Number of boundary face elements = ',nb
        WRITE(6,*)
        !
        ! Compare face coordinates to determine all boundary conditions
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
100        FORMAT (2x,a6,i3,2x,a6,a18)
        ENDDO
        !
        ier = 0
        perio = 0
        WRITE(6,*)
150     WRITE(6,'('' Enter number of periodic boundary face pairs: '')', ADVANCE = "NO")
        READ(5,'(a)', IOSTAT = ier) dir
        dir = ADJUSTL (dir)
        l1 = INDEX(dir,' ')+1
        l2 = LEN_TRIM(dir)
        READ(dir(1:l1),*) perio 

        IF (ier /= 0) GOTO 150
        !
        DO n=1,perio
           !
           sum_face = 0
           !
           WRITE(6,*)
120        WRITE(6,'('' Specify the two IDs defining the periodic pair: '')', ADVANCE = "NO")
           WRITE(6,*)
           READ(5,'(a)') dir
           dir = ADJUSTL (dir)
           l1 = INDEX(dir,' ')+1
           l2 = LEN_TRIM(dir)

           READ(dir(1:l1),*) bc_par(1)
           READ(dir(l1:l2),*) bc_par(2) 
           !
           ! Check if number of element faces are equal
           !
           l1 = num_point(bc_par(1),2)-num_point(bc_par(1),1)+1
           l2 = num_point(bc_par(2),2)-num_point(bc_par(2),1)+1
           !
           IF (l1/=l2) THEN
              WRITE(6,*) 
              WRITE(6,*) '  Number of element faces are not equal: ',l1,l2
              WRITE(6,*) 
              GOTO 120
           ENDIF
           !
           ! Find translation vector 
           !
           trans = 0
           diffx = 0.
           diffy = 0.
           diffz = 0.
           nmean = 0.
           dist = 0.
           !
           per_loop: DO j=1,2
              !
              i1 = num_point(bc_par(j),1)
              i2 = num_point(bc_par(j),2)
              nn = 0
              !
              DO iel = 1,inc_2*l1,inc_2
                 !
                 nn = nn +1
                 vrt(1:4) = elements(1+bc_par(j),iel:iel+3)
                 !
                 CALL norm_triangle (n1,tri,vrt)
                 !
                 nmean(j,1)= nmean(j,1) +0.25*(n1(1,1)+n1(2,1)+n1(3,1)+n1(4,1))
                 nmean(j,2)= nmean(j,2) +0.25*(n1(1,2)+n1(2,2)+n1(3,2)+n1(4,2))
                 nmean(j,3)= nmean(j,3) +0.25*(n1(1,3)+n1(2,3)+n1(3,3)+n1(4,3))
                 !
                 dist(j,1) = dist(j,1) +tri(5,1)
                 dist(j,2) = dist(j,2) +tri(5,2)
                 dist(j,3) = dist(j,3) +tri(5,3)
                 !
              ENDDO
              !
              nmean(j,:)= nmean(j,:)/nn
              dist(j,:) = dist(j,:)/nn
              !
           ENDDO per_loop
           !
           ! check if they are translational periodic planes
           !
           IF ( ABS(nmean(1,1))-ABS(nmean(2,1)) < eps .AND. ABS(nmean(1,2))-ABS(nmean(2,2)) < eps &
                .AND. ABS(nmean(1,3))-ABS(nmean(2,3)) < eps) THEN
              int_dir = 0
              diffx = dist(2,1)-dist(1,1)
              diffy = dist(2,2)-dist(1,2)
              diffz = dist(2,3)-dist(1,3)
              !
              !
              ! Check if they are rotational periodic planes
              !
           ELSEIF ( ABS(nmean(1,1)-nmean(2,1)) > eps .AND. ABS(nmean(1,2)-nmean(2,2)) > eps &
                .AND. nmean(1,3)-nmean(2,3) < eps) THEN
              int_dir = 1
              rot_axe = 3 
              rot_dir = nmean(1,1)*nmean(2,2)-nmean(1,2)*nmean(2,1)
              udotv = nmean(1,1)*nmean(2,1)+nmean(1,2)*nmean(2,2)
              uabs  = SQRT(nmean(1,1)*nmean(1,1)+nmean(1,2)*nmean(1,2))
              vabs  = SQRT(nmean(2,1)*nmean(2,1)+nmean(2,2)*nmean(2,2))
              cosn  = ACOS(udotv/(uabs*vabs))
              beta = PI-cosn
           ELSEIF ( ABS(nmean(1,1)-nmean(2,1)) > eps .AND. ABS(nmean(1,2)-nmean(2,2)) < eps &
                .AND. ABS(nmean(1,3)-nmean(2,3)) > eps) THEN
              int_dir = 1
              rot_axe = 2 
              rot_dir = nmean(1,3)*nmean(2,1)-nmean(1,1)*nmean(2,3)
              udotv = nmean(1,1)*nmean(2,1)+nmean(1,3)*nmean(2,3)
              uabs  = SQRT(nmean(1,1)*nmean(1,1)+nmean(1,3)*nmean(1,3))
              vabs  = SQRT(nmean(2,1)*nmean(2,1)+nmean(2,3)*nmean(2,3))
              cosn  = ACOS(udotv/(uabs*vabs))
              beta = PI-cosn

           ELSEIF ( ABS(nmean(1,1)-nmean(2,1)) < eps .AND. ABS(nmean(1,2)-nmean(2,2)) > eps &
                .AND. ABS(nmean(1,3)-nmean(2,3)) > eps) THEN
              int_dir = 1
              rot_axe = 1 
              rot_dir = nmean(1,3)*nmean(2,2)-nmean(1,3)*nmean(2,1)
              udotv = nmean(1,3)*nmean(2,3)+nmean(1,2)*nmean(2,2)
              uabs  = SQRT(nmean(1,3)*nmean(1,3)+nmean(1,2)*nmean(1,2))
              vabs  = SQRT(nmean(2,3)*nmean(2,3)+nmean(2,2)*nmean(2,2))
              cosn  = ACOS(udotv/(uabs*vabs))
              beta = PI-cosn
           ELSE
              WRITE(6,*) 
              WRITE(6,*) '      No periodic bcs  found !!! '
              CALL cgnstonek_exit('       ')
           ENDIF
           !
           IF (int_dir == 0) THEN
              !
              !$OMP parallel do
              !
              DO i=num_point(bc_par(1),1),num_point(bc_par(1),2)
                 tra_loop: DO j= num_point(bc_par(2),1),num_point(bc_par(2),2)
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
                       sum_face(j) = 1
                       EXIT tra_loop
                    ENDIF
                 ENDDO tra_loop
                 !
              ENDDO
              !
              !$OMP end parallel do
              !
              sum_f = SUM(sum_face)
              !
              ii = per_fam(n)+1
              !
              IF (sum_f < nelem_end(1+bc_par(1))-nelem_start(1+bc_par(1))+1) THEN
                 WRITE(6,*) 
                 WRITE(6,*) '      NOT enough periodic elements found !!! in BC = ',ii
                 CALL cgnstonek_exit('       ')
              ENDIF
              !
           ELSEIF (int_dir == 1) THEN
              !
              !$OMP parallel do
              !
              DO i=num_point(bc_par(1),1),num_point(bc_par(1),2)
                 IF (rot_axe == 1) THEN
                    radius = SQRT(yb(i)**2 +zb(i)**2)
                    diffx = 0.
                    diffy = radius*(1.-COS(beta))*SIGN(dble(1.),rot_dir)
                    diffz = -radius*SIN(beta)*SIGN(dble(1.),rot_dir)
                 ELSEIF (rot_axe == 2) THEN
                    radius = SQRT(xb(i)**2 +zb(i)**2)
                    diffx = radius*(1.-COS(beta))*SIGN(dble(1.),rot_dir)
                    diffy = 0.
                    diffz = -radius*SIN(beta)*SIGN(dble(1.),rot_dir)
                 ELSEIF (rot_axe == 3) THEN
                    radius = SQRT(xb(i)**2 +yb(i)**2)
                    diffx = radius*(1.-COS(beta))*SIGN(dble(1.),rot_dir)
                    diffy = -radius*SIN(beta)*SIGN(dble(1.),rot_dir)
                    diffz = 0.
                 ENDIF

                 rot_loop: DO j= num_point(bc_par(2),1),num_point(bc_par(2),2)
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
                       sum_face(j) = 1
                       EXIT rot_loop
                    ENDIF
                 ENDDO rot_loop
                 !
              ENDDO
              !
              !$OMP end parallel do
              !
              sum_f = SUM(sum_face)
              !
              ii = per_fam(n)+1
              !
              IF (sum_f < nelem_end(bc_par(1))-nelem_start(bc_par(1))+1) THEN
                 WRITE(6,*)
                 WRITE(6,*) '      NOT enough periodic elements found !!! in BC = ',ii
                 CALL cgnstonek_exit('       ')
              ENDIF

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
!
SUBROUTINE norm_triangle (n1,tri,vrt)
  !
  !  Calculate the normal vectors of triangles
  !
  USE CGNS
  USE module_global
  !
  IMPLICIT NONE

  INTEGER(cgsize_t) :: i,vrt(4)
  !
  REAL(high) :: ax,ay,az,bx,by,bz
  REAL(high) :: tri(5,3),n1(4,3)
  !
  DO i=1,4
     tri(i,1) = xv(vrt(i))
     tri(i,2) = yv(vrt(i))
     tri(i,3) = zv(vrt(i))
  ENDDO
  !
  !  center point of the plane
  !
  tri(5,1) = 0.25*(xv(vrt(1)) +xv(vrt(2)) +xv(vrt(3)) +xv(vrt(4)))
  tri(5,2) = 0.25*(yv(vrt(1)) +yv(vrt(2)) +yv(vrt(3)) +yv(vrt(4)))
  tri(5,3) = 0.25*(zv(vrt(1)) +zv(vrt(2)) +zv(vrt(3)) +zv(vrt(4)))
  !
  ! Normalvector of plane 1 computed on 4 subtriangles
  !
  ax = tri(2,1)-tri(1,1)
  ay = tri(2,2)-tri(1,2)
  az = tri(2,3)-tri(1,3)
  !
  bx = tri(5,1)-tri(1,1)
  by = tri(5,2)-tri(1,2)
  bz = tri(5,3)-tri(1,3)
  !
  n1(1,1)= ay*bz - az*by
  n1(1,2)= az*bx - ax*bz
  n1(1,3)= ax*by - ay*bx
  !
  ax = tri(3,1)-tri(2,1)
  ay = tri(3,2)-tri(2,2)
  az = tri(3,3)-tri(2,3)
  !
  bx = tri(5,1)-tri(2,1)
  by = tri(5,2)-tri(2,2)
  bz = tri(5,3)-tri(2,3)
  !
  n1(2,1)= ay*bz - az*by
  n1(2,2)= az*bx - ax*bz
  n1(2,3)= ax*by - ay*bx
  !
  ax = tri(4,1)-tri(3,1)
  ay = tri(4,2)-tri(3,2)
  az = tri(4,3)-tri(3,3)
  !
  bx = tri(5,1)-tri(3,1)
  by = tri(5,2)-tri(3,2)
  bz = tri(5,3)-tri(3,3)
  !
  n1(3,1)= ay*bz - az*by
  n1(3,2)= az*bx - ax*bz
  n1(3,3)= ax*by - ay*bx
  !
  ax = tri(1,1)-tri(4,1)
  ay = tri(1,2)-tri(4,2)
  az = tri(1,3)-tri(4,3)
  !
  bx = tri(5,1)-tri(4,1)
  by = tri(5,2)-tri(4,2)
  bz = tri(5,3)-tri(4,3)
  !
  n1(4,1)= ay*bz - az*by
  n1(4,2)= az*bx - ax*bz
  n1(4,3)= ax*by - ay*bx
  !
  RETURN
  !
END SUBROUTINE norm_triangle
