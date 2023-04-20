module SIZE
!
!
! EXODUS related variables:
!
      character(32) exoname ! MXSTLN=32 max string length in exodus file
	  
      character,allocatable,dimension(:,:) :: fluidexo
      character,allocatable,dimension(:,:) :: solidexo
	  
      integer,allocatable,dimension(:) :: idblk
      integer,allocatable,dimension(:) :: num_nodes_per_elem
      integer,allocatable,dimension(:) :: num_attr    !not used

      real,save,allocatable,dimension(:) :: x_exo, y_exo, z_exo

      integer,save,allocatable,dimension(:)   :: num_elem_in_block, connect
      integer,save,allocatable,dimension(:)   :: num_sides_in_set, idss
      integer,save,allocatable,dimension(:,:) :: elem_list, side_list
	  integer,save,allocatable,dimension(:,:,:)   :: parray
      
      integer*8  num_dim, num_elem, num_elem_blk, nvert
      integer*8  num_side_sets, num_sides_tot
      save     num_dim, num_elem, num_elem_blk, nvert
      save     num_side_sets, num_sides_tot
	  
      integer*8  etot,etot_est,eacc,eacc_old,eftot
      save     etot,etot_est,eacc,eacc_old,eftot

      real    shiftvector(3)
      save    shiftvector
	  
      integer  fnexo,snexo,iexo ! number of exo files
      integer,save,allocatable,dimension(:)   :: f_elem_exo,s_elem_exo ! converted nek elements number for each exo file
      integer,save,allocatable,dimension(:)   :: ssinfo

      real*8    maxxyz(3),minxyz(3)
      save    maxxyz,minxyz
	 
      integer bcNumber
      integer,save,allocatable,dimension(:)   ::  bcID
	  
      integer converting_option,converting_option_old
	  
      integer quadratic_option
	  
	  
	  
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
