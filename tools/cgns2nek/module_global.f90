!**************************************************************
!
!
! Mesh converter from cgns to nek format
!
!
!********************************************************************
!
MODULE module_global
  !
  IMPLICIT NONE
  !
  SAVE
  !
  !
  ! Select this line for DOUBLE precision of real data
  !
  INTEGER, PARAMETER :: high = 8
  INTEGER, PARAMETER :: iprec= 4
  INTEGER, PARAMETER :: single_p= 4
  !
  INTEGER(iprec) :: max_elem,energy,scalar,ifile
  !
  INTEGER(iprec), ALLOCATABLE :: index_zone(:),nbdyelem(:),sectype(:)
  INTEGER(iprec), ALLOCATABLE,DIMENSION (:) :: zonetype,index_dim,nelem_start,nelem_end
  INTEGER(iprec), ALLOCATABLE,DIMENSION (:) :: elementdatasize,tmp_Elements,tmp_ParentData, &
                  Normallist,location,bocotype,bc_fam
  INTEGER(iprec), ALLOCATABLE,DIMENSION (:,:) :: Elements,ParentData
  !
  INTEGER(iprec) :: inc,cgns_len,sec,bcsec,nbocos,nsections,ncoords,nzones,nbases,pre
  !
  REAL (high), ALLOCATABLE, DIMENSION (:) :: xv,yv,zv,xb,yb,zb
  REAL (high), ALLOCATABLE, DIMENSION (:,:) :: xf,yf,zf
  REAL (high), ALLOCATABLE, DIMENSION (:,:,:) :: vel_val
  REAL (high), ALLOCATABLE, DIMENSION (:, :) :: x_nodes
  REAL (single_p) :: version
  REAL (single_p), ALLOCATABLE, DIMENSION (:) :: xvsingle,yvsingle,zvsingle
  !
  CHARACTER (len=80) :: cgns_file,basename
  CHARACTER (len=32),ALLOCATABLE, DIMENSION (:) :: fam_name,zonename,coord_name
  CHARACTER (len=32),ALLOCATABLE,DIMENSION (:) ::sec_type,sec_name,connectname,donorname,boconame, &
                                            datasetName,BcType
  CHARACTER (len=10) :: project,ascii
  CHARACTER (len=100),ALLOCATABLE, DIMENSION (:, :) :: bc_values
  CHARACTER (len=3),ALLOCATABLE,DIMENSION (:,:) :: vel_type,tem_type,scal_type
  CHARACTER (len=32) :: bc_nek(25)
  !
END MODULE module_global
