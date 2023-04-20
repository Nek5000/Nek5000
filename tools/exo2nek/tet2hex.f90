!--------------------------------------------------------------------
      subroutine split_convert_new
!converting hybrid tet4-hex8-wedge6 mesh
! 
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
      integer*8 ehexnumber,tetnumber,wedgenumber,bctot
      integer*8 vert_index_exo
      integer*8 iel_nek,iel_exo,ifc_exo,iel_exo_g
	  
	  integer iblk,iblk1

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo,iel_exo_g
      eacc_old = eacc
      
      call rzero_int(exoss,6*num_elem)
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

      !call rzero_int(hexss,6*8) 
      ! write(6,*) 'flag1'
      vert_index_exo = 0
      iel_nek =  eacc_old !0
	  ! write(6,*) 'flag2'

      do iblk = 1,num_elem_blk
	  
        write(6,*) 'Converting elements in block ',idblk(iblk)
        ! detect element type of this block
        nvert=num_nodes_per_elem(iblk)
        write(6,*) 'nvert, ', nvert
        if (nvert.eq.4) then ! tet4
          do iel_exo = 1,num_elem_in_block(iblk)
           !write(6,*) 'iel_exo, ',iel_exo

         iel_exo_g = iel_exo
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif

!  read tet 4 . 
!  linear interpolate to tet10

       do ivert = 1, nvert
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
       !call rzero_int(tetss,4)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,4
          tetss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
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
          enddo

        else if (nvert.eq.8) then ! hex8
          do iel_exo = 1,num_elem_in_block(iblk)

         iel_exo_g = iel_exo
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif
		  
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
       !call rzero_int(ehexss,6)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,6
          ehexss(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

       call ehexto8hex(hexver,ehexver,hexss,ehexss) 

       do ihex = 1,8
          iel_nek = iel_nek + 1
          eacc = eacc + 1
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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

          enddo
        else if (nvert.eq.6) then ! wedge6
          do iel_exo = 1,num_elem_in_block(iblk)

         iel_exo_g = iel_exo
         if (iblk.gt.1) then
            do iblk1 = 1,iblk-1
             iel_exo_g = iel_exo_g + num_elem_in_block(iblk1)
            enddo
         endif
		  
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
       !call rzero_int(wedgess,5)
       if (num_side_sets.ne.0) then
        do ifc_exo=1,5
          wedgess(ifc_exo) = exoss(ifc_exo,iel_exo_g)
        enddo
       endif

!  given wedge15 vertices, and return you 6 hex coords.
       call wedgetohex2(hexver,wedgever,hexss,wedgess)

       do ihex = 1,6
          iel_nek = iel_nek + 1
          eacc = eacc + 1
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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
      integer*8 ehexnumber,tetnumber,wedgenumber,bctot
      integer*8 vert_index_exo
      integer*8 iel_nek,iel_exo,ifc_exo

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo

      eacc_old = eacc

      call rzero_int(exoss,6*num_elem)
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
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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
      integer*8 ehexnumber,tetnumber,wedgenumber,bctot
      integer*8 vert_index_exo
      integer*8 iel_nek,iel_exo,ifc_exo

      save ehexnumber,tetnumber,wedgenumber,bctot
      save vert_index_exo
      save iel_nek,iel_exo,ifc_exo
      eacc_old = eacc
      
      call rzero_int(exoss,6*num_elem)
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
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

		   
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
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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
        if (eacc.gt.etot_est) write(6,*) 'ERROR: please increase estimate final total hex element number' 

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