module SIZE
!
!
! EXODUS related variables:
!
      character(32) exoname ! MXSTLN=32 max string length in exodus file

      real,save,allocatable,dimension(:) :: x_exo, y_exo, z_exo

      integer,save,allocatable,dimension(:)   :: num_elem_in_block, connect
      integer,save,allocatable,dimension(:)   :: num_sides_in_set, idss
      integer,save,allocatable,dimension(:,:) :: elem_list, side_list
      integer  num_dim, num_elem, num_elem_blk, nvert
      integer  num_side_sets, num_sides_tot
      save     num_dim, num_elem, num_elem_blk, nvert
      save     num_side_sets, num_sides_tot
!
!
! NEK CORE variables:
!
      real,save,allocatable,dimension(:,:,:)   ::  bc, curve
      real,save,allocatable,dimension(:,:,:,:) ::  xm1, ym1, zm1

      character(1),save,allocatable,dimension(:,:) :: ccurve
      character(3),save,allocatable,dimension(:,:) :: cbc
!
!
! .RE2 file related variables
!
      character(80) re2name
!
end module SIZE
