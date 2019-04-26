!-----------------------------------------------------------------------
      program gmsh2nek

      use SIZE

      integer option
      write(6,'(A)', ADVANCE = "NO") 'Enter mesh dimension: '
      read (5,'(I1)') option
 
      if (option.eq.2) then
         call read_input_name
		 if(aorb.eq.0) then
         call gmsh_read_2d_ascii
         elseif(aorb.eq.1) then
         call gmsh_read_2d_binary
         endif
         call convert_2d
         call setbc_2d
      elseif(option.eq.3) then
         call read_input_name
         if(aorb.eq.0) then
         call gmsh_read_3d_ascii
         elseif(aorb.eq.1) then
         call gmsh_read_3d_binary
         endif
         call convert_3d
         call setbc_3d
      else
        write(6,*) 'Unknown input option'
        STOP
      endif

      call deallocate_all_msh_arrays  ! deallocate_all_msh_arrays to save memory
      call gen_re2                    ! write nek mesh re2 file 

      end 
!-----------------------------------------------------------------------
      subroutine read_input_name

      use SIZE
      character*80 charline
      real A
      integer B
	  
      character(1) re2nam1(80)
      character(1) mshnam1(32)
      character(32) fname

      write(6,'(A)', ADVANCE = "NO") 'Input .msh file name: '
      read (5,'(a32)') fname
      len = ltrunc(fname,32)

      call blank  (mshnam1, 32)
      call blank  (re2nam1, 80)
      call chcopy (mshnam1,fname,32)
      call chcopy (re2nam1,fname,80)
	  
      call chcopy(mshnam1(len+1) ,'.msh' , 4)
      call chcopy(re2nam1(len+1) ,'.re2' , 4)
	  
      call blank (mshname, 32)
      call blank (re2name, 80)
      call chcopy (mshname,mshnam1,len+4)
      call chcopy (re2name,re2nam1,len+4)

      open(301,file=mshname)
      read(301,*) charline
      read(301,*) A,aorb,B
      close(301)
	  
      if ((A.ge.3.0).or.(A.lt.2.0)) then
      write(6,*) 'ERROR: invalid msh file format!'  
      STOP
      endif
	  
! aorb indicates ascii or binary file
	  
      return 
      end
!-----------------------------------------------------------------------
      subroutine gmsh_read_2d_ascii
! read .msh file (version 2, ascii format)
! 2d mesh, ascii
      use SIZE

      character*32  mshnam2
      character*1   mshnam3(32)
      character*80 charline
      character*1  charlin1(80)
      integer A,B,C,elemType

      equivalence(mshnam2,mshnam3)
      equivalence(charline,charlin1) 	  

      call chcopy(mshnam2,mshname,32)	  
      len = ltrunc(mshnam2,32)
      call chcopy(mshnam3(len+1) ,'_1' , 2)


      call blank (charline,80)
      call chcopy(charlin1(1),'cp ',3)
	
      len = ltrunc(mshname,32)
      call chcopy(charlin1(4),mshname,len)
	
      len = ltrunc(charline,80)
      call chcopy(charlin1(len+1),' ',1)

      len2 = ltrunc(mshnam2,32)
      call chcopy(charlin1(len+2),mshnam2,len2)

      call system(charline)	  

      call blank (charline,80)

      open(299,file=mshname)
      open(300,file=mshnam2)

! loop to find $PhysicalNames
      do while (.true.) 
        read(299,*) charline
        read(300,*)
        charline = trim(charline)
        if (charline.eq."$PhysicalNames") goto 1010
      enddo
! end loop to $PhysicalNames
1010  read(299,*) bcNumber ! bcNumber is number of boundaries
      read(300,*)
      !bcNumber	= bcNumber - 1 
      allocate ( bcID       (2,bcNumber))
      allocate ( bcChar     (bcNumber))
      call rzero_int(bcID,2*bcNumber)
      call blank  (bcChar, 32*bcNumber)

      ibc_a = 0
      do ibc= 1,bcNumber
      read(299,*) A,bcID(1,ibc),bcChar(ibc)
      read(300,*)
      !write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
        if(A.EQ.1) ibc_a = ibc_a + 1
      enddo
      bcNumber = ibc_a

! loop to find Nodes section
      do while (.true.) 
        read(299,*) charline
        read(300,*)
        charline = trim(charline)
        if (charline.eq."$Nodes") goto 1020
      enddo
! end loop to "$Nodes"

! read all nodes xyz
1020  read(299,*) totalNode
      read(300,*)

! now we know total node number, allocate memory size.
      allocate ( node_xyz       (3,totalNode))
      allocate ( node_line      (2,totalNode))
      allocate ( node_quad      (4,totalNode))	  
      call rzero(node_xyz,      3*totalNode)
      call rzero_int(node_line,     2*totalNode) 
      call rzero_int(node_quad,     4*totalNode)

! read all node xyz.
      do inode = 1,totalNode
      read(299,*)A,node_xyz(1,inode),node_xyz(2,inode) &
      ,node_xyz(3,inode)
      read(300,*)
      enddo
! end read all nodes xyz

      read(299,*)charline ! "$EndNodes"
      read(300,*)
      read(299,*)charline ! "$Elements"
      read(300,*)
	  
      read(299,*)totalElem
      read(300,*)

      allocate ( line_array       (5,totalElem))
      allocate ( quad_array       (11,totalElem))
      call rzero_int(line_array,      5*totalElem)
      call rzero_int(quad_array,      11*totalElem)

      totalLine = 0
      totalQuad = 0

      do iElem= 1,totalElem
      read(299,*) A,elemType

	  ! detemine element type
      if (elemType.eq.8) then ! if line3
      totalLine = totalLine + 1
      read(300,*) A,B,C, &
      line_array(1,totalLine),line_array(2,totalLine),&
      line_array(3,totalLine),line_array(4,totalLine),&
      line_array(5,totalLine)

      do inode = 1,2
      call addTo_node_line(line_array(2+inode,totalLine),totalLine)
      enddo
 
      do ibc= 1,bcNumber
        if(line_array(1,totalLine).eq.bcID(1,ibc)) then
          bcID(2,ibc) =   bcID(2,ibc) + 1
        endif
      enddo

      elseif (elemType.eq.16) then ! if quad8
      totalQuad = totalQuad + 1
      read(300,*) A,B,C, &
      quad_array(1,totalQuad),quad_array(2,totalQuad),&
      quad_array(3,totalQuad),quad_array(4,totalQuad),&
      quad_array(5,totalQuad),quad_array(6,totalQuad),&
      quad_array(7,totalQuad),quad_array(8,totalQuad),&
      quad_array(9,totalQuad),quad_array(10,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo
  
      elseif (elemType.eq.10) then ! if quad9
      totalQuad = totalQuad + 1
      read(300,*) A,B,C,&
      quad_array(1,totalQuad),quad_array(2,totalQuad),&
      quad_array(3,totalQuad),quad_array(4,totalQuad),&
      quad_array(5,totalQuad),quad_array(6,totalQuad),&
      quad_array(7,totalQuad),quad_array(8,totalQuad),&
      quad_array(9,totalQuad),quad_array(10,totalQuad),&
      quad_array(11,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo 

      else 
       ! unknown element type 
       ! only line3 and quad8/9 elements are accepted.
	    write (6,*) 'ERROR: unknown element type'
        write (6,*) 'only line3 and quad8/9 elements are accepted'
        write (6,*) 'please choose "set order 2" option to set all elements to 2nd order'
        write (6,*) 'please uncheck "save all elements" when exporting mesh'
        write (6,*) 'please see readme for more information'
        STOP
      endif

      enddo

      write (6,*) 'total node number is ', totalNode
      write (6,*) 'total line element number is ', totalLine
      write (6,*) 'total quad element number is ', totalQuad

      close(299)
      close(300)

      num_dim = 2
      num_elem = totalQuad
 
      call blank (charline,80)
      call chcopy(charlin1(1),'rm ',3)
	
      len = ltrunc(mshnam2,32)
      call chcopy(charlin1(4),mshnam2,len)
	
      call system(charline)

      return
      end
!-----------------------------------------------------------------------
      subroutine gmsh_read_2d_binary
! read .msh file (version 2, binary format)
! 2d mesh
      use SIZE

      character*1   singlechar(100)
      character*100 charline
 	  logical ifbswap
      integer idummy(100),buf(100),buf1(100)
      integer bone,nlength
      integer elem_type,num_elm_follow,num_tags
      integer fileid
	 
      fileid = 302
	  
	  ! read msh file in binary format.
      open(unit=fileid,file=mshname,access="stream",form="unformatted",status="old")
	  
      ! read two lines.
! ------------------------------------------------------------------
      call bread_line(fileid,singlechar,nlength)
      call bread_line(fileid,singlechar,nlength)

! 1. test little or big endian.	  
! if binary one,  then no need to bit swap, ifbswap = false
! if not binray one, then need to bit swap, ifbswap = true

      read(fileid) bone
      read(fileid) singlechar(1)  ! move cursor to next line
      if (bone.eq.1) then
       ifbswap = .false.
      else
       ifbswap = .true.
      endif	   

! loop to find $PhysicalNames
      do while (.true.) 
         call bread_line(fileid,singlechar,nlength)
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         if (charline.eq."$PhysicalNames") goto 1010
      enddo
! end loop to $PhysicalNames
1010  call bread_line(fileid,singlechar,nlength) !bcNumber is number of boundaries
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) bcNumber  
      !bcNumber	= bcNumber - 1 
      allocate ( bcID       (2,bcNumber))
      allocate ( bcChar     (bcNumber))
      call rzero_int(bcID,2*bcNumber)
      call blank  (bcChar, 32*bcNumber)
      ibc_a = 0
      do ibc= 1,bcNumber
      call bread_line(fileid,singlechar,nlength) !bcNumber is number of boundaries
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) A,bcID(1,ibc),bcChar(ibc)
      !write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
        if(A.EQ.1) ibc_a = ibc_a + 1
      enddo
      bcNumber = ibc_a
	  
! 2. loop lines to $Node
! loop to find Nodes section
      do while (.true.) 
         call bread_line(fileid,singlechar,nlength)
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         if (charline.eq."$Nodes") goto 1130
      enddo
! end loop to "$Nodes"
1130     call bread_line(fileid,singlechar,nlength) ! read node number
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         read(charline,*) totalNode  ! convert char to int

! allocate memory size for node array
      allocate ( node_xyz       (3,totalNode))
      allocate ( node_line      (2,totalNode))
      allocate ( node_quad      (4,totalNode))	  
      call rzero(node_xyz,      3*totalNode)
      call rzero_int(node_line,     2*totalNode) 
      call rzero_int(node_quad,     4*totalNode)

! 2. Get total node number, loop all nodes
      do inode = 1,totalNode
      read(fileid) idummy(1),node_xyz(1,inode),node_xyz(2,inode)&
      ,node_xyz(3,inode)
      if(ifbswap) then
          call endian_swap_8(node_xyz(1,inode))
          call endian_swap_8(node_xyz(2,inode))
          call endian_swap_8(node_xyz(3,inode))  
      endif
      enddo

      read(fileid) singlechar(1) ! move cursor to next line
! jump $EndNodes
      call bread_line(fileid,singlechar,nlength)
! jump $Elements		 
      call bread_line(fileid,singlechar,nlength)
! read totalElem		 
      call bread_line(fileid,singlechar,nlength)
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) totalElem

      allocate ( line_array       (5,totalElem))
      allocate ( quad_array       (11,totalElem))
      call rzero_int(line_array,      5*totalElem)
      call rzero_int(quad_array,      11*totalElem)

      totalLine = 0
      totalQuad = 0

      do while (.true.) 
      read(fileid) elem_type,num_elm_follow, num_tags

	    if (elem_type.eq.8) then ! if line3
          do iLine= 1,num_elm_follow
          totalLine = totalLine + 1
          read(fileid) idummy(1),&
          line_array(1,totalLine),line_array(2,totalLine),&
          line_array(3,totalLine),line_array(4,totalLine),&
          line_array(5,totalLine)

          if(ifbswap) then
            do i = 1,5
              call endian_swap_4(line_array(i,totalLine))
            enddo
          endif
		  
          do inode = 1,2
          call addTo_node_line(line_array(2+inode,totalLine),totalLine)
          enddo
 
          do ibc= 1,bcNumber
           if(line_array(1,totalLine).eq.bcID(1,ibc)) then
             bcID(2,ibc) =   bcID(2,ibc) + 1
           endif
          enddo

          enddo
	   elseif (elem_type.eq.16) then ! quad 8
	       do iQuad= 1,num_elm_follow
           totalQuad = totalQuad + 1
           read(fileid) idummy(1),&
           quad_array(1,totalQuad),quad_array(2,totalQuad),&
           quad_array(3,totalQuad),quad_array(4,totalQuad),&
           quad_array(5,totalQuad),quad_array(6,totalQuad),&
           quad_array(7,totalQuad),quad_array(8,totalQuad),&
           quad_array(9,totalQuad),quad_array(10,totalQuad)

           if(ifbswap) then
             do i = 1,10
               call endian_swap_4(quad_array(i,totalQuad))
             enddo
           endif

           do inode = 1,4
             call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
           enddo

           enddo
       elseif (elem_type.eq.10) then ! quad 9
	       do iQuad= 1,num_elm_follow
           totalQuad = totalQuad + 1
           read(fileid) idummy(1),&
           quad_array(1,totalQuad),quad_array(2,totalQuad),&
           quad_array(3,totalQuad),quad_array(4,totalQuad),&
           quad_array(5,totalQuad),quad_array(6,totalQuad),&
           quad_array(7,totalQuad),quad_array(8,totalQuad),&
           quad_array(9,totalQuad),quad_array(10,totalQuad),&
           quad_array(11,totalQuad)
		
           if(ifbswap) then
             do i = 1,11
               call endian_swap_4(quad_array(i,totalQuad))
             enddo
           endif

           do inode = 1,4
             call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
           enddo

           enddo
      else
       ! unknown element type 
       ! only line3 and quad8/9 elements are accepted.
	    write (6,*) 'ERRPOR: unknown element type'
        write (6,*) 'only line3 and quad8/9 elements are accepted'
        write (6,*) 'please choose "set order 2" option to set all elements to 2nd order'
        write (6,*) 'please uncheck "save all elements" when exporting mesh'
        write (6,*) 'please see readme for more information'
        STOP
      endif 

      if((totalQuad+totalLine).eq.totalElem) goto 1180
      enddo

! close file
1180  close(fileid)

      write (6,*) 'total node number is ', totalNode
      write (6,*) 'total line element number is ', totalLine
      write (6,*) 'total quad element number is ', totalQuad

      num_dim = 2
      num_elem = totalQuad
 
      return
      end
!-----------------------------------------------------------------------
      subroutine gmsh_read_3d_ascii
! read .msh file (version 2, ascii format)
! 3d mesh, ascii
      use SIZE

      character*32  mshnam2
      character*1   mshnam3(32)
      character*80 charline
      character*1  charlin1(80)
      integer A,B,C,elemType

      equivalence(mshnam2,mshnam3)
      equivalence(charline,charlin1) 	  

      call chcopy(mshnam2,mshname,32)	  
      len = ltrunc(mshnam2,32)
      call chcopy(mshnam3(len+1) ,'_1' , 2)


      call blank (charline,80)
      call chcopy(charlin1(1),'cp ',3)
	
      len = ltrunc(mshname,32)
      call chcopy(charlin1(4),mshname,len)
	
      len = ltrunc(charline,80)
      call chcopy(charlin1(len+1),' ',1)

      len2 = ltrunc(mshnam2,32)
      call chcopy(charlin1(len+2),mshnam2,len2)

      call system(charline)	  

      call blank (charline,80)

      open(299,file=mshname)
      open(300,file=mshnam2)

! loop to find $PhysicalNames
      do while (.true.) 
        read(299,*) charline
        read(300,*)
        charline = trim(charline)
        if (charline.eq."$PhysicalNames") goto 1010
      enddo
! end loop to $PhysicalNames
1010  read(299,*) bcNumber ! bcNumber is number of boundaries
      read(300,*)
      !bcNumber	= bcNumber - 1 
      allocate ( bcID       (2,bcNumber))
      allocate ( bcChar     (bcNumber))
      call rzero_int(bcID,2*bcNumber)
      call blank  (bcChar, 32*bcNumber)
      ibc_a = 0
      do ibc= 1,bcNumber
      read(299,*) A,bcID(1,ibc),bcChar(ibc)
      read(300,*)
      !write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
         if(A.EQ.2) ibc_a = ibc_a + 1
      enddo
      bcNumber = ibc_a
! loop to find Nodes section
      do while (.true.) 
        read(299,*) charline
        read(300,*)
        charline = trim(charline)
        if (charline.eq."$Nodes") goto 1020
      enddo
! end loop to "$Nodes"

! read all nodes xyz
1020  read(299,*) totalNode
      read(300,*)

! now we know total node number, allocate memory size.
      allocate ( node_xyz       (3,totalNode))
      allocate ( node_quad      (4,totalNode))	  
      allocate ( node_hex       (8,totalNode))
      call rzero(node_xyz,      3*totalNode)
      call rzero_int(node_quad,     4*totalNode) 
      call rzero_int(node_hex,      8*totalNode)

! read all node xyz.
      do inode = 1,totalNode
      read(299,*)A,node_xyz(1,inode),node_xyz(2,inode) &
      ,node_xyz(3,inode)
      read(300,*)
      enddo
! end read all nodes xyz

      read(299,*)charline ! "$EndNodes"
      read(300,*)
      read(299,*)charline ! "$Elements"
      read(300,*)
	  
      read(299,*)totalElem
      read(300,*)
! msh (version2, ascci) only tells us the total element number,
! including all quad+hex elements.
! but we do not know the specific number of quads and hexs
      allocate ( quad_array       (11,totalElem)) 	  
      allocate ( hex_array        (29,totalElem)) 
      call rzero_int(quad_array,      11*totalElem)
      call rzero_int(hex_array,       29*totalElem) 

      totalQuad = 0
      totalHex = 0

      do iElem= 1,totalElem
      read(299,*) A,elemType

	  ! detemine element type
      if (elemType.eq.16) then ! if quad8
      totalQuad = totalQuad + 1
      read(300,*) A,B,C, &
      quad_array(1,totalQuad),quad_array(2,totalQuad),&
      quad_array(3,totalQuad),quad_array(4,totalQuad),&
      quad_array(5,totalQuad),quad_array(6,totalQuad),&
      quad_array(7,totalQuad),quad_array(8,totalQuad),&
      quad_array(9,totalQuad),quad_array(10,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo
	    
      do ibc= 1,bcNumber
        if(quad_array(1,totalQuad).eq.bcID(1,ibc)) then
          bcID(2,ibc) =   bcID(2,ibc) + 1
        endif
      enddo
 
      elseif (elemType.eq.10) then ! if quad9
      totalQuad = totalQuad + 1
      read(300,*) A,B,C,&
      quad_array(1,totalQuad),quad_array(2,totalQuad),&
      quad_array(3,totalQuad),quad_array(4,totalQuad),&
      quad_array(5,totalQuad),quad_array(6,totalQuad),&
      quad_array(7,totalQuad),quad_array(8,totalQuad),&
      quad_array(9,totalQuad),quad_array(10,totalQuad),&
      quad_array(11,totalQuad)

      do inode = 1,4
      call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
      enddo 
	  
      do ibc= 1,bcNumber
        if(quad_array(1,totalQuad).eq.bcID(1,ibc)) then
          bcID(2,ibc) =   bcID(2,ibc) + 1
        endif
      enddo
	 
      elseif (elemType.eq.17) then ! if hex20
      totalHex = totalHex + 1
      read(300,*) A,B,C,&
      hex_array(1,totalHex),hex_array(2,totalHex),&
      hex_array(3,totalHex),hex_array(4,totalHex),&
      hex_array(5,totalHex),hex_array(6,totalHex),&
      hex_array(7,totalHex),hex_array(8,totalHex),&
      hex_array(9,totalHex),hex_array(10,totalHex),&
      hex_array(11,totalHex),hex_array(12,totalHex),&
      hex_array(13,totalHex),hex_array(14,totalHex),&
      hex_array(15,totalHex),hex_array(16,totalHex),&
      hex_array(17,totalHex),hex_array(18,totalHex),&
      hex_array(19,totalHex),hex_array(20,totalHex),&
      hex_array(21,totalHex),hex_array(22,totalHex)
	 
      do inode = 1,8
      call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
      enddo

      elseif (elemType.eq.12) then ! if hex27
      totalHex = totalHex + 1
      read(300,*) A,B,C,&
      hex_array(1,totalHex),hex_array(2,totalHex),&
      hex_array(3,totalHex),hex_array(4,totalHex),&
      hex_array(5,totalHex),hex_array(6,totalHex),&
      hex_array(7,totalHex),hex_array(8,totalHex),&
      hex_array(9,totalHex),hex_array(10,totalHex),&
      hex_array(11,totalHex),hex_array(12,totalHex),&
      hex_array(13,totalHex),hex_array(14,totalHex),&
      hex_array(15,totalHex),hex_array(16,totalHex),&
      hex_array(17,totalHex),hex_array(18,totalHex),&
      hex_array(19,totalHex),hex_array(20,totalHex),&
      hex_array(21,totalHex),hex_array(22,totalHex),&
      hex_array(23,totalHex),hex_array(24,totalHex),&
      hex_array(25,totalHex),hex_array(26,totalHex),&
      hex_array(27,totalHex),hex_array(28,totalHex),&
      hex_array(29,totalHex)
      
      do inode = 1,8
      call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
      enddo

      else 
       ! unknown element type 
       ! only quad8/9 and hex20/27 elements are accepted.
	    write (6,*) 'ERRPOR: unknown element type'
        write (6,*) 'only quad8/9 and hex20/27 elements are accepted'
        write (6,*) 'please choose "set order 2" option to set all elements to 2nd order'
        write (6,*) 'please uncheck "save all elements" when exporting mesh'
        write (6,*) 'please see readme for more information'
        STOP
      endif

      enddo

      write (6,*) 'total node number is ', totalNode
      write (6,*) 'total quad element number is ', totalQuad
      write (6,*) 'total hex element number is ', totalHex

      close(299)
      close(300)

      num_dim = 3
      num_elem = totalHex
 
      call blank (charline,80)
      call chcopy(charlin1(1),'rm ',3)
	
      len = ltrunc(mshnam2,32)
      call chcopy(charlin1(4),mshnam2,len)
	
      call system(charline)

      return
      end
!-----------------------------------------------------------------
      subroutine gmsh_read_3d_binary
! read .msh file (version 2, binary format)
! looking for read_re2_data for reference on how to read binary files. 
! using nek build in subroutines.
! 
      use SIZE

      character*1   singlechar(100)
      character*100 charline
 	  logical ifbswap
      integer idummy(100),buf(100),buf1(100)
      integer bone,nlength
      integer elem_type,num_elm_follow,num_tags
      integer fileid
	 
      fileid = 302
	  
	  ! read msh file in binary format.
      open(unit=fileid,file=mshname,access="stream",form="unformatted",status="old")
	  
      ! read two lines.
! ------------------------------------------------------------------
      call bread_line(fileid,singlechar,nlength)
      call bread_line(fileid,singlechar,nlength)

! 1. test little or big endian.	  
! if binary one,  then no need to bit swap, ifbswap = false
! if not binray one, then need to bit swap, ifbswap = true

      read(fileid) bone
      read(fileid) singlechar(1)  ! move cursor to next line
      if (bone.eq.1) then
       ifbswap = .false.
      else
       ifbswap = .true.
      endif	   

! loop to find $PhysicalNames
      do while (.true.) 
         call bread_line(fileid,singlechar,nlength)
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         if (charline.eq."$PhysicalNames") goto 1010
      enddo
! end loop to $PhysicalNames
1010  call bread_line(fileid,singlechar,nlength) !bcNumber is number of boundaries
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) bcNumber  
      !bcNumber	= bcNumber - 1 
      allocate ( bcID       (2,bcNumber))
      allocate ( bcChar     (bcNumber))
      call rzero_int(bcID,2*bcNumber)
      call blank  (bcChar, 32*bcNumber)
      ibc_a = 0
      do ibc= 1,bcNumber
      call bread_line(fileid,singlechar,nlength) !bcNumber is number of boundaries
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) A,bcID(1,ibc),bcChar(ibc)
      !write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
        if(A.EQ.2) ibc_a = ibc_a + 1
      enddo
      bcNumber = ibc_a
	  
! 2. loop lines to $Node
! loop to find Nodes section
      do while (.true.) 
         call bread_line(fileid,singlechar,nlength)
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         if (charline.eq."$Nodes") goto 1130
      enddo
! end loop to "$Nodes"
1130     call bread_line(fileid,singlechar,nlength) ! read node number
         call blank (charline,100)
         call chcopy(charline,singlechar,nlength-1)
         charline = trim(charline)
         read(charline,*) totalNode  ! convert char to int

! allocate memory size for node array
      allocate ( node_xyz       (3,totalNode))
      allocate ( node_quad      (4,totalNode))	  
      allocate ( node_hex       (8,totalNode))
      call rzero(node_xyz,      3*totalNode)
      call rzero_int(node_quad,     4*totalNode) 
      call rzero_int(node_hex,      8*totalNode)
	 
	 
! 2. Get total node number, loop all nodes
      do inode = 1,totalNode
      read(fileid) idummy(1),node_xyz(1,inode),node_xyz(2,inode)&
      ,node_xyz(3,inode)
      if(ifbswap) then
          call endian_swap_8(node_xyz(1,inode))
          call endian_swap_8(node_xyz(2,inode))
          call endian_swap_8(node_xyz(3,inode))  
      endif
      enddo

      read(fileid) singlechar(1) ! move cursor to next line
! jump $EndNodes
      call bread_line(fileid,singlechar,nlength)
! jump $Elements		 
      call bread_line(fileid,singlechar,nlength)
! read totalElem		 
      call bread_line(fileid,singlechar,nlength)
      call blank (charline,100)
      call chcopy(charline,singlechar,nlength-1)
      charline = trim(charline)
      read(charline,*) totalElem

! 3. Get total element number
! allocate memory size for element array
      allocate ( quad_array       (11,totalElem)) 	  
      allocate ( hex_array        (29,totalElem)) 
      call rzero_int(quad_array,      11*totalElem)
      call rzero_int(hex_array,       29*totalElem) 

      totalQuad = 0
      totalHex = 0

      do while (.true.) 
       read(fileid) elem_type,num_elm_follow, num_tags

	   if (elem_type.eq.16) then ! quad 8
	       do iQuad= 1,num_elm_follow
           totalQuad = totalQuad + 1
           read(fileid) idummy(1),&
           quad_array(1,totalQuad),quad_array(2,totalQuad),&
           quad_array(3,totalQuad),quad_array(4,totalQuad),&
           quad_array(5,totalQuad),quad_array(6,totalQuad),&
           quad_array(7,totalQuad),quad_array(8,totalQuad),&
           quad_array(9,totalQuad),quad_array(10,totalQuad)

           if(ifbswap) then
             do i = 1,10
               call endian_swap_4(quad_array(i,totalQuad))
             enddo
           endif

           do inode = 1,4
             call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
           enddo
		   
           do ibc= 1,bcNumber
             if(quad_array(1,totalQuad).eq.bcID(1,ibc)) then
             bcID(2,ibc) =   bcID(2,ibc) + 1
             endif
           enddo

           enddo
       elseif (elem_type.eq.10) then ! quad 9
	       do iQuad= 1,num_elm_follow
           totalQuad = totalQuad + 1
           read(fileid) idummy(1),&
           quad_array(1,totalQuad),quad_array(2,totalQuad),&
           quad_array(3,totalQuad),quad_array(4,totalQuad),&
           quad_array(5,totalQuad),quad_array(6,totalQuad),&
           quad_array(7,totalQuad),quad_array(8,totalQuad),&
           quad_array(9,totalQuad),quad_array(10,totalQuad),&
           quad_array(11,totalQuad)
		
           if(ifbswap) then
             do i = 1,11
               call endian_swap_4(quad_array(i,totalQuad))
             enddo
           endif

           do inode = 1,4
             call addTo_node_quad(quad_array(2+inode,totalQuad),totalQuad)
           enddo

           do ibc= 1,bcNumber
             if(quad_array(1,totalQuad).eq.bcID(1,ibc)) then
             bcID(2,ibc) =   bcID(2,ibc) + 1
             endif
           enddo

           enddo
	  elseif (elem_type.eq.17) then ! if hex20
	       do iHex= 1,num_elm_follow
           totalHex = totalHex + 1
           read(fileid) idummy(1),&
           hex_array(1,totalHex),hex_array(2,totalHex),&
           hex_array(3,totalHex),hex_array(4,totalHex),&
           hex_array(5,totalHex),hex_array(6,totalHex),&
           hex_array(7,totalHex),hex_array(8,totalHex),&
           hex_array(9,totalHex),hex_array(10,totalHex),&
           hex_array(11,totalHex),hex_array(12,totalHex),&
           hex_array(13,totalHex),hex_array(14,totalHex),&
           hex_array(15,totalHex),hex_array(16,totalHex),&
           hex_array(17,totalHex),hex_array(18,totalHex),&
           hex_array(19,totalHex),hex_array(20,totalHex),&
           hex_array(21,totalHex),hex_array(22,totalHex)
		  
           if(ifbswap) then
             do i = 1,22
               call endian_swap_4(hex_array(i,totalHex))
             enddo
           endif

	       do inode = 1,8
             call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
           enddo

           enddo
	  elseif (elem_type.eq.12) then ! if hex20
	       do iHex= 1,num_elm_follow
           totalHex = totalHex + 1
           read(fileid) idummy(1),&
           hex_array(1,totalHex),hex_array(2,totalHex),&
           hex_array(3,totalHex),hex_array(4,totalHex),&
           hex_array(5,totalHex),hex_array(6,totalHex),&
           hex_array(7,totalHex),hex_array(8,totalHex),&
           hex_array(9,totalHex),hex_array(10,totalHex),&
           hex_array(11,totalHex),hex_array(12,totalHex),&
           hex_array(13,totalHex),hex_array(14,totalHex),&
           hex_array(15,totalHex),hex_array(16,totalHex),&
           hex_array(17,totalHex),hex_array(18,totalHex),&
           hex_array(19,totalHex),hex_array(20,totalHex),&
           hex_array(21,totalHex),hex_array(22,totalHex),&
           hex_array(23,totalHex),hex_array(24,totalHex),&
           hex_array(25,totalHex),hex_array(26,totalHex),&
           hex_array(27,totalHex),hex_array(28,totalHex),&
           hex_array(29,totalHex)
		   
           if(ifbswap) then
             do i = 1,29
               call endian_swap_4(hex_array(i,totalHex))
             enddo
           endif

           do inode = 1,8
             call addTo_node_hex(hex_array(2+inode,totalHex),totalHex)
           enddo

           enddo
      else
       ! unknown element type 
       ! only quad8/9 and hex20/27 elements are accepted.
	    write (6,*) 'ERRPOR: unknown element type'
        write (6,*) 'only quad8/9 and hex20/27 elements are accepted'
        write (6,*) 'please choose "set order 2" option to set all elements to 2nd order'
        write (6,*) 'please uncheck "save all elements" when exporting mesh'
        write (6,*) 'please see readme for more information'
        STOP
      endif 

      if((totalQuad+totalHex).eq.totalElem) goto 1180
      enddo

! close file
1180  close(fileid)

      write (6,*) 'total node number is ', totalNode
      write (6,*) 'total quad element number is ', totalQuad
      write (6,*) 'total hex element number is ', totalHex

      num_dim = 3
      num_elem = totalHex

      return
      end
!-----------------------------------------------------------------------
      subroutine convert_2d
!  Subroutine to convert gmsh quad8/quad9 to nek  quad8/quad9 elements.
      use SIZE

      integer msh_to_nek_right(9)
      data    msh_to_nek_right /1,3,9,7,2,6,8,4,5/ ! RIGHT HAND SIDE ELEMENT
      integer msh_to_nek_left(9)
      data    msh_to_nek_left /3,1,7,9,2,4,8,6,5/  ! LEFT HAND SIDE ELEMENT
 
      integer quad_face_node_right(2,4)
      data quad_face_node_right /1,2,2,3,3,4,4,1/ ! RIGHT HAND SIDE ELEMENT
      integer quad_face_node_left(2,4)
      data quad_face_node_left /2,1,1,4,4,3,3,2/  ! LEFT HAND SIDE ELEMENT
 
      integer lnode(2)
      integer iQuad,iline,imshvert,inekvert
      integer ifoundline,imatch
      integer physicalTag
      integer iaddhex,addhex
      logical ifnew
 
      allocate ( xm1              (3,3,3,num_elem))
      allocate ( ym1              (3,3,3,num_elem))
      allocate ( zm1              (3,3,3,num_elem))
      allocate ( quad_line_array   (4,num_elem))
      call rzero_int(quad_line_array,  4*num_elem)

      allocate (r_or_l(num_elem))
	  call rzero_int(r_or_l,num_elem)

! need a test to figure out if element is right-hand or left-hand.
      do iQuad = 1, totalQuad
! detect right or left hand elements
      call r_or_l_detect(iQuad,r_or_l(iQuad))
      if(r_or_l(iQuad).eq.0) then !  for right hand element in gmsh mesh
        do imshvert = 1,9 
          inekvert = msh_to_nek_right(imshvert)
          xm1(inekvert,1,1,iQuad)= node_xyz(1,quad_array(imshvert+2,iQuad))
          ym1(inekvert,1,1,iQuad)= node_xyz(2,quad_array(imshvert+2,iQuad))
        enddo
      elseif(r_or_l(iQuad).eq.1) then !  for left hand element in gmsh mesh
        do imshvert = 1,9
          inekvert = msh_to_nek_left(imshvert)
          xm1(inekvert,1,1,iQuad)= node_xyz(1,quad_array(imshvert+2,iQuad))
          ym1(inekvert,1,1,iQuad)= node_xyz(2,quad_array(imshvert+2,iQuad))
        enddo
      endif
      enddo

! zero-out bc and curve sides arrays

      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
      allocate   (curve  (2*num_dim,12,   num_elem) )
      call rzero (curve,2*num_dim*12*num_elem)
      call blank (ccurve,(4+8*(num_dim-2))*num_elem)

      allocate   (cbc    (2*num_dim,      num_elem) )
      allocate   (bc     (5,2*num_dim,    num_elem) ) 
      call rzero (bc,5*2*num_dim*num_elem)
      call blank (cbc,3*2*num_dim*num_elem)

!---- search for boundaries, based on physical tags of line elements
!---- use old scheme
!---- because usually 2d mesh will not be too large.
!---- so this should be enough for 2d mesh, and it is guaranteed to work.
      do iQuad = 1,totalQuad
        do iline = 1,4
            ! obtain node id for this line on this quad.
            if(r_or_l(iQuad).eq.0) then
              do ilnode = 1,2
              lnode(ilnode)=quad_array(quad_face_node_right(ilnode,iline)+2,iQuad)
              enddo
            elseif(r_or_l(iQuad).eq.1) then
              do ilnode = 1,2
              lnode(ilnode)=quad_array(quad_face_node_left(ilnode,iline)+2,iQuad)
              enddo
            endif

            call findline(lnode,ifoundline)   ! ifoundline is the matching line number
                                              ! ifoundline is 0 if no line element find
            if(ifoundline.ne.0) then
            quad_line_array(iline,iQuad) = line_array(1,ifoundline) ! physical tag
            endif
        enddo
      enddo

! assign dummy boundary condition and id to bc array	  
      do iQuad= 1,totalQuad
        do iline = 1,4
         if((quad_line_array(iline,iQuad)).ne.0) then       ! if on boundary, with physical tag
          cbc(iline,iQuad) = 'MSH'                         ! dummy boundary condition
          bc(5,iline,iQuad) = quad_line_array(iline,iQuad) ! assign tag 
         endif
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------
      subroutine convert_3d
!
!  Subroutine to convert gmsh hex20/hex27 to nek hex20/hex27 elements.

      use SIZE

      integer msh_to_nek(27) ! gmsh hex27 node order map to nek hex27 node order
      data    msh_to_nek &
           /1,3,9,7,19,21,27,25,2,4,10,6,12, &
            8,18,16,20,22,24,26,5,11,13,15,17,23,14/
 
      integer hex_face_node(4,6)
      data hex_face_node &
           /1,2,6,5,2,3,7,6,3,4,8,7,1,4,8,5, &
            1,2,3,4,5,6,7,8/
 
      integer fnode(4),quadnode(4),fhex(32),fhex_nd(32)
      integer ihex,imshvert,inekvert,iquad,iface,inode,ifnode
      integer ifoundquad,imatch
      integer physicalTag
      integer iaddhex,addhex,iaddhex_nd,iaddhex1,iaddhex2
      logical ifnew

      allocate ( xm1              (3,3,3,num_elem))
      allocate ( ym1              (3,3,3,num_elem))
      allocate ( zm1              (3,3,3,num_elem))
      allocate ( hex_face_array   (6,num_elem))
      call rzero_int(hex_face_array,  6*num_elem)

      do ihex = 1, totalHex
        do imshvert = 1,27 ! 20
          inekvert = msh_to_nek(imshvert)
          xm1(inekvert,1,1,ihex)= node_xyz(1,hex_array(imshvert+2,ihex))
          ym1(inekvert,1,1,ihex)= node_xyz(2,hex_array(imshvert+2,ihex))
          zm1(inekvert,1,1,ihex)= node_xyz(3,hex_array(imshvert+2,ihex))
        enddo
      enddo

! zero-out bc and curve sides arrays

      allocate   (ccurve (4+8*(num_dim-2),num_elem) )
      allocate   (curve  (2*num_dim,12,   num_elem) )
      call rzero (curve,2*num_dim*12*num_elem)
      call blank (ccurve,(4+8*(num_dim-2))*num_elem)

      allocate   (cbc    (2*num_dim,      num_elem) )
      allocate   (bc     (5,2*num_dim,    num_elem) ) 
      call rzero (bc,5*2*num_dim*num_elem)
      call blank (cbc,3*2*num_dim*num_elem)

! currently, does consider converting boundary condition now.
! only need to associate quad4 to hex8 faces.
!
!---- new scheme to search for boundaries.
!---- because node_quad, and node_hex now contains information of this node.
! node_quad((1-4),inode) contains the quad ids associate with this inode
! node_hex((1-8),inode) contains the hex ids associate with this inode
! this new scheme O(N) should be much faster than old scheme O(N^2)
!
! search boundary using loop in all quads.
      do iquad = 1,totalQuad
         physicalTag = quad_array(1,iquad)
	     do inode = 1,4
          quadnode(inode) = quad_array(inode+2,iquad) ! first 4 nodes of quad
         enddo

!find all hexes that share the same nodes of this quad
         call rzero_int(fhex,32)
         call rzero_int(fhex_nd,32)
         iaddhex = 0
         do inode = 1,4
            do i = 1,8
               ! find all hex id related to the nodes in this quad.
			   ! there will be duplicated hex id 
               if(node_hex(i,quadnode(inode)).gt.0) then
               iaddhex = iaddhex + 1
               fhex(iaddhex) = node_hex(i,quadnode(inode))
               endif
            enddo
         enddo

! fhex(32) contains all hexes that share the same nodes of this quad.
! now, only need to loop over this fhex(32) to find which hex face correspond to this quad.
! however, there are duplicated hex id in fhex(32)
! eliminate duplicated hex id in fhex, and store to fhex_nd

       addhex = iaddhex
       iaddhex_nd = 1
       fhex_nd(1) = fhex(1)
       do iaddhex1 = 2,addhex  ! loop in fhex 
            ifnew = .TRUE.
            do iaddhex2 =1,iaddhex_nd ! loop in fhex_nd
			! if duplicate hex id, ifnew = false
             if(fhex(iaddhex1).eq.fhex_nd(iaddhex2)) ifnew=.FALSE.		
            enddo
            if(ifnew) then
              iaddhex_nd = iaddhex_nd + 1
              fhex_nd(iaddhex_nd) = fhex(iaddhex1)
            endif
       enddo

! look over in fhex_nd to find 
! which hex which face is corresponding to this quad
!
       do iaddhex = 1,iaddhex_nd
         ihex = fhex_nd(iaddhex) ! hex id that are related to the nodes of quad
         do iface = 1,6       ! loop all faces
           do ifnode = 1,4
           fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
           enddo
			  ! fnode is the node ids of hex face
			  ! quadnode is the node ids of quad
           imatch = 0 
           call ifquadmatch(imatch,fnode,quadnode)
           if(imatch.eq.1) then
             hex_face_array(iface,ihex) = physicalTag
             if(hex_face_array(iface,ihex).eq.1) flag1 = flag1 +1
		     goto 1100
           endif
         enddo
       enddo
!
! if no hex face is found for this quad, report error
      write(6,*) 'ERROR: cannot find hex face for quad id ',iquad
      write(6,*) 'ERROR: this should not happen, please check your mesh'
      write(6,*) 'ERROR: or your mesh exporting process in gmsh'

1100  continue
      enddo

!---- old scheme to search for boundary setup
!---- this is very slow O(N^2) scheme. 	  
!      do ihex = 1, totalHex
!        do iface = 1,6
!            ! obtain node id for this face.
!            do ifnode = 1,4
!            fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
!            enddo
!            !input fnode, return quad elements number iquad
!            call findquad(fnode,ifoundquad)   ! ifoundquad is the matching quad number
!                                              ! ifoundquad is 0 if no quad element find
!            if(ifoundquad.ne.0) then
!            hex_face_array(iface,ihex) = quad_array(1,ifoundquad) ! physical tag
!                                         !quad_array(2,ifoundquad) ! geometrical tag
!            endif
!        enddo
!      enddo

! assign dummy boundary condition and id to bc array
      do ihex = 1, totalHex
        do iface = 1,6
         if(hex_face_array(iface,ihex).ne.0) then       ! if on boundary, with physical tag
          cbc(iface,ihex) = 'MSH'                       ! dummy boundary condition
          bc(5,iface,ihex) = hex_face_array(iface,ihex) ! assign tag 
         endif
        enddo
      enddo
 
      return
      end
!!-----------------------------------------------------------------
      subroutine setbc_2d
      use SIZE
!  set boundary condition for 2d mesh
!  read from a file casename.bc

      integer msh_to_nek_right(9)
      data    msh_to_nek_right /1,3,9,7,2,6,8,4,5/ ! RIGHT HAND SIDE ELEMENT
      integer msh_to_nek_left(9)
      data    msh_to_nek_left /3,1,7,9,2,4,8,6,5/  ! LEFT HAND SIDE ELEMENT
 
      integer quad_face_node_right(2,4)
      data quad_face_node_right /1,2,2,3,3,4,4,1/ ! RIGHT HAND SIDE ELEMENT
      integer quad_face_node_left(2,4)
      data quad_face_node_left /2,1,1,4,4,3,3,2/  ! LEFT HAND SIDE ELEMENT
	 
      integer parray(2,2,totalLine)

      character*3 ubc
      integer tags(2),ibc,nbc,io
      integer ip,np,ipe,ipe2,nipe(2)
      integer ptags(2)
      integer lnode(2)
      real pvec(3)
      real fpxyz(3,2)
      real AB_v(3),AD_v(3),farea,product_v(3)
      real dist,distMax,ptol

! boundary condition summary
      write(6,*) '******************************************************'
      write(6,*) 'Boundary info summary'
      write(6,*) 'BoundaryName     BoundaryID'
      do ibc= 1,bcNumber
      write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
      enddo
      write(6,*) '******************************************************'
  
      write(6,*) 'Enter number of periodic boundary surface pairs:'
      read (5,*) nbc
	  
      if(nbc.le.0) return
  
      do ibc = 1,nbc 
        ptol = 1e-5
        write(6,*) 'input surface 1 and  surface 2  BoundaryID'
        read (5,*) ptags(1),ptags(2)
        write(6,*) 'input translation vector (surface 1 -> surface 2)'
        read (5,*) pvec(1),pvec(2),pvec(3)

          ipe = 0
          do iquad = 1,totalQuad
            do iline = 1,4
               if(bc(5,iline,iquad).eq.ptags(1)) then
                ipe = ipe + 1
                parray(1,1,ipe) = iquad
                parray(2,1,ipe) = iline
               endif
            enddo
          enddo
          nipe(1) = ipe
	  
          ipe = 0
	      do iquad = 1,totalQuad
            do iline = 1,4
               if(bc(5,iline,iquad).eq.ptags(2)) then
                ipe = ipe + 1
                parray(1,2,ipe) = iquad
                parray(2,2,ipe) = iline
               endif
            enddo
          enddo
          nipe(2) = ipe

          if(nipe(1).ne.nipe(2))  then
            write(6,*)'mapping surface',ptags(1),'with',nipe(1),'lines'
            write(6,*)'to surface',ptags(2),'with',nipe(2),'lines'
            write(6,*) 'EORROR, line numbers are not matching'
          endif

          do ipe = 1,nipe(1)
            iquad = parray(1,1,ipe)
            iline = parray(2,1,ipe)
! get face center xyz
            call rzero(fpxyz(1,1),3) 
            if(r_or_l(iquad).eq.0) then
              do ilnode = 1,2
              lnode(ilnode)=quad_array(quad_face_node_right(ilnode,iline)+2,iquad)
              enddo
            elseif(r_or_l(iQuad).eq.1) then
              do ilnode = 1,2
              lnode(ilnode)=quad_array(quad_face_node_left(ilnode,iline)+2,iquad)
              enddo
            endif
            do ilnode = 1,2
            fpxyz(1,1) = fpxyz(1,1) + node_xyz(1,lnode(ilnode))/2.0
            fpxyz(2,1) = fpxyz(2,1) + node_xyz(2,lnode(ilnode))/2.0
            fpxyz(3,1) = fpxyz(3,1) + node_xyz(3,lnode(ilnode))/2.0
            enddo

! lopp over mapped faces to find its mapping face
            distMax = 1000.0
            do ipe2 = 1,nipe(2)
              iquad2 = parray(1,2,ipe2)
              iline2 = parray(2,2,ipe2)
! get face center xyz
              call rzero(fpxyz(1,2),3)
              if(r_or_l(iquad2).eq.0) then
                do ilnode = 1,2
                lnode(ilnode)=quad_array(quad_face_node_right(ilnode,iline2)+2,iquad2)
                enddo
              elseif(r_or_l(iquad2).eq.1) then
                do ilnode = 1,2
                lnode(ilnode)=quad_array(quad_face_node_left(ilnode,iline2)+2,iquad2)
                enddo
              endif
              do ilnode = 1,2
              fpxyz(1,2) = fpxyz(1,2) + node_xyz(1,lnode(ilnode))/2.0
              fpxyz(2,2) = fpxyz(2,2) + node_xyz(2,lnode(ilnode))/2.0
              fpxyz(3,2) = fpxyz(3,2) + node_xyz(3,lnode(ilnode))/2.0
              enddo
 
              dist = sqrt((fpxyz(1,2) - fpxyz(1,1) - pvec(1))**2 &
       + (fpxyz(2,2) - fpxyz(2,1) - pvec(2))**2 &
       + (fpxyz(3,2) - fpxyz(3,1) - pvec(3))**2)
 
               if(dist.lt.distMax) then 
                  distMax = dist
                  !write(6,*) distMax
                  if(distMax.le.ptol) then
                  bc(1,iline,iquad) = iquad2*1.0
                  bc(2,iline,iquad) = iline2*1.0
                  bc(1,iline2,iquad2) = iquad*1.0
                  bc(2,iline2,iquad2) = iline*1.0
                  cbc(iline,iquad) = 'P  '
                  cbc(iline2,iquad2) = 'P  '
                  endif
               endif
             enddo
          enddo

          nperror = 0

          do ipe = 1,nipe(1)
             iquad = parray(1,1,ipe)
             iline = parray(2,1,ipe)
             if (cbc(iline,iquad).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) then
          write(6,*)'doing periodic check for surface',ptags(1)
          write(6,*) 'ERROR,',nperror,'lines did not map'
          endif

          nperror = 0

          do ipe = 1,nipe(1)
             iquad = parray(1,1,ipe)
             iline = parray(2,1,ipe)
             iquad2 = int(bc(1,iline,iquad))
             iline2 = int(bc(2,iline,iquad))
             iquad3 = int(bc(1,iline2,iquad2))
             iline3 = int(bc(2,iline2,iquad2))
             if ((iquad.ne.iquad3).or.(iline.ne.iline3)) then
                nperror = nperror + 1
             endif
          enddo		  

          if (nperror.gt.0) then
          write(6,*)'doing periodic check for surface',ptags(1)
          write(6,*) 'ERROR,',nperror,'lines are wrong',&
      'out of total ',nipe(1),' lines'
          endif
      
      enddo

      write(6,*) '******************************************************'
      write(6,*) 'Please set boundary conditions to all non-periodic boundaries'
      write(6,*) 'in .usr file usrdat2() subroutine'
      write(6,*) '******************************************************'

      return
      end
!--------------------------------------------------------------------
      subroutine setbc_3d
      use SIZE
!  set boundary condition for 3d mesh
!  read from a file casename.bc

      integer hex_face_node(4,6)
      data hex_face_node &
           /1,2,6,5,2,3,7,6,3,4,8,7,1,4,8,5, &
            1,2,3,4,5,6,7,8/
	 
      integer parray(2,2,totalQuad)

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
      write(6,*) 'BoundaryName     BoundaryID'
      do ibc= 1,bcNumber
      write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
      enddo
      write(6,*) '******************************************************'
 
 
      write(6,*) 'Enter number of periodic boundary surface pairs:'
      read (5,*) nbc
	  
      if(nbc.le.0) return

      do ibc = 1,nbc 
        ptol = 1e-5
        write(6,*) 'input surface 1 and  surface 2  BoundaryID'
        read (5,*) ptags(1),ptags(2)
        write(6,*) 'input translation vector (surface 1 -> surface 2)'
        read (5,*) pvec(1),pvec(2),pvec(3)

          ipe = 0
          do ihex = 1, totalHex
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
	      do ihex = 1, totalHex
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
            write(6,*)'maping surface',ptags(1),'with',nipe(1),'faces'
            write(6,*)'to surface',ptags(2),'with',nipe(2),'faces'
            write(6,*) 'EORROR, face numbers are not matching'
          endif

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
! get face center xyz
             call rzero(fpxyz(1,1),3)
             do ifnode = 1,4
             fnode(ifnode)=hex_array(hex_face_node(ifnode,iface)+2,ihex)
             fpxyz(1,1) = fpxyz(1,1) + node_xyz(1,fnode(ifnode))/4.0
             fpxyz(2,1) = fpxyz(2,1) + node_xyz(2,fnode(ifnode))/4.0
             fpxyz(3,1) = fpxyz(3,1) + node_xyz(3,fnode(ifnode))/4.0
             enddo

! lopp over mapped faces to find its mapping face
             distMax = 1000.0
             do ipe2 = 1,nipe(2)
                ihex2 = parray(1,2,ipe2)
                iface2 = parray(2,2,ipe2)
! get face center xyz
               call rzero(fpxyz(1,2),3)
               do ifnode = 1,4
           fnode(ifnode)=hex_array(hex_face_node(ifnode,iface2)+2,ihex2)
               fpxyz(1,2) = fpxyz(1,2) + node_xyz(1,fnode(ifnode))/4.0
               fpxyz(2,2) = fpxyz(2,2) + node_xyz(2,fnode(ifnode))/4.0
               fpxyz(3,2) = fpxyz(3,2) + node_xyz(3,fnode(ifnode))/4.0
               enddo
               
              dist = sqrt((fpxyz(1,2) - fpxyz(1,1) - pvec(1))**2 &
       + (fpxyz(2,2) - fpxyz(2,1) - pvec(2))**2 &
       + (fpxyz(3,2) - fpxyz(3,1) - pvec(3))**2)
 
               if(dist.lt.distMax) then 
                  distMax = dist
                  !write(6,*) distMax
                  if(distMax.le.ptol) then
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
                  bc(1,iface2,ihex2) = ihex*1.0
                  bc(2,iface2,ihex2) = iface*1.0
                  cbc(iface,ihex) = 'P  '
                  cbc(iface2,ihex2) = 'P  '
                  endif
               endif
             enddo
          enddo

          nperror = 0

          do ipe = 1,nipe(1)
             ihex = parray(1,1,ipe)
             iface = parray(2,1,ipe)
             if (cbc(iface,ihex).ne.'P  ') then
                  nperror = nperror +1 
             endif
          enddo
          if (nperror.gt.0) then
          write(6,*)'doing periodic check for surface',ptags(1)
          write(6,*) 'ERROR,',nperror,'faces did not map'
          endif
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
          write(6,*)'doing periodic check for surface',ptags(1)
          write(6,*) 'ERROR,',nperror,'faces are wrong',&
      'out of total ',nipe(1),' faces'
          endif
      
      enddo

	  
      write(6,*) '******************************************************'
      write(6,*) 'Please set boundary conditions to all non-periodic boundaries'
      write(6,*) 'in .usr file usrdat2() subroutine'
      write(6,*) '******************************************************'
	  
      return
      end
!-----------------------------------------------------------------------
      subroutine deallocate_all_msh_arrays
! deallocate msh file related arrays
      use SIZE


      if(num_dim.eq.2) then
      deallocate(node_xyz,node_quad,node_line)
      deallocate(quad_array,line_array,quad_line_array)
      deallocate(r_or_l)
      endif

      if(num_dim.eq.3) then
      deallocate(node_xyz,node_quad,node_hex)
      deallocate(quad_array,hex_array,hex_face_array)
      endif

      return 
      end
!------------------------------------------------------------------------------------------
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

      character*80  hdr


      real*4 test
      data   test  / 6.54321 /

      call byte_open(re2name,ierr)
            
!  Write the header
      call blank     (hdr,80)    
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

      character*1 cc

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
            call byte_write (buf2,16, ierr)
          endif
        enddo
      enddo

      return
      end
!-----------------------------------------------------------------------
      subroutine write_bc
      
      use SIZE

      real*8  rbc, buf2(30)

      character*3 ch3
      character*1 chdum
      data        chdum /' '/

      nface = 2*num_dim
      nbc   = 0

      do iel=1,num_elem
        do ifc=1,nface
          if (cbc(ifc,iel).ne.'   ')  nbc = nbc + 1
        enddo
      enddo

      rbc = nbc
      call byte_write (rbc,2, ierr)

      do iel = 1,num_elem
        do ifc = 1,nface
          ch3 = cbc(ifc,iel)
!          if (ch3.eq.'MSH') then
          if (ch3.ne.'   ') then ! setting boundary condition.
            buf2(1)=iel
            buf2(2)=ifc
            call copy   (buf2(3),bc(1,ifc,iel),5)
            call blank  (buf2(8),8)
            call chcopy (buf2(8),ch3,3)
            if (num_elem.ge.1000000) then
              ibc     = bc(1,ifc,iel)
              buf2(3) = ibc
            endif
            call byte_write (buf2,16, ierr)
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

      real        len
      real        x3(27),y3(27),z3(27),xyz(3,3)
      character*1 ccrve(12)
      integer     e,edge

      integer e3(3,12)
      save    e3
      data    e3 /  1, 2, 3,    3, 6, 9,    9, 8, 7,    7, 4, 1 &
                , 19,20,21,   21,24,27,   27,26,25,   25,22,19  &
                ,  1,10,19,    3,12,21,    9,18,27,    7,16,25 /

      call chcopy(ccrve,ccurve(1,e),12)

      call map2reg     (x3,3,xm1(1,1,1,e),1)  ! Map to 3x3x3 array
      call map2reg     (y3,3,ym1(1,1,1,e),1)
      if (num_dim.eq.3) call map2reg    (z3,3,zm1(1,1,1,e),1)

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
      CHARACTER*1 A(1)
      CHARACTER*1 BLNK
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
      CHARACTER*1 A(1), B(1)

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
      CHARACTER*1 STRING(L)
      CHARACTER*1   BLNK
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
 100     A(I) = 0.0
      return
      END
!-----------------------------------------------------------------------
      subroutine rzero_int(a,n)
      integer A(1)
      DO 100 I = 1, N
 100     A(I) = 0
      return
      END
!-----------------------------------------------------------------------
      subroutine endian_swap_4(var)
! swap 4 byte variables, typical a int, real*4
      character var_swap1(4),var_swap2(4)

      call chcopy(var_swap1,var,4)
      call chcopy(var_swap2,var,4)

      var_swap1(1) = var_swap2(4)
      var_swap1(2) = var_swap2(3)
      var_swap1(3) = var_swap2(2)
      var_swap1(4) = var_swap2(1)

      call chcopy(var,var_swap1,4)

      return
      end
!-----------------------------------------------------------------------
      subroutine endian_swap_8(var)
! swap 4 byte variables, typical a real*8
! swap 4 byte variables, typical a int, real*4
      character var_swap1(8),var_swap2(8)

      call chcopy(var_swap1,var,8)
      call chcopy(var_swap2,var,8)

      var_swap1(1) = var_swap2(8)
      var_swap1(2) = var_swap2(7)
      var_swap1(3) = var_swap2(6)
      var_swap1(4) = var_swap2(5)
      var_swap1(5) = var_swap2(4)
      var_swap1(6) = var_swap2(3)
      var_swap1(7) = var_swap2(2)
      var_swap1(8) = var_swap2(1)

      call chcopy(var,var_swap1,8)
	  
      return
      end
!-----------------------------------------------------------------------
      subroutine bread_line(fileid,singlechar,nlength)
! read binary file. read each byte, convert to char
! convert to ASCII, until new line (char(10)) symbol is met
      integer fileid,nlength
      character*1 singlechar(1)

      nlength = 0
      do nlength = 1,100	  
	  read(fileid) singlechar(nlength)
      !write(6,*) singlechar(nlength)
      if (singlechar(nlength).eq.char(10)) then
      return
      endif  
      enddo

      return 
      end
!-----------------------------------------------------------------------
      subroutine addTo_node_line(inode,iline)
      use SIZE
      integer inode,iline
      integer ilinestart
      ilinestart = 1

      do i = 1,2
        if(node_line(i,inode).eq.0) then
        ilinestart = i
        node_line(ilinestart,inode) = iline
        return
        endif
      enddo

      return 
      end
!-----------------------------------------------------------------------
      subroutine addTo_node_quad(inode,iquad)
      use SIZE
      integer inode,iquad
      integer iquadstart
      iquadstart = 1

      do i = 1,4
        if(node_quad(i,inode).eq.0) then
        iquadstart = i
        node_quad(iquadstart,inode) = iquad
        return
        endif
      enddo

      return 
      end
!-----------------------------------------------------------------
      subroutine addTo_node_hex(inode,ihex)
      use SIZE
      integer inode,ihex
      integer ihexstart
      ihexstart = 1

      do i = 1,8
        if(node_hex(i,inode).eq.0) then
        ihexstart = i
        node_hex(ihexstart,inode) = ihex
        return
        endif
      enddo

      return 
      end
!-----------------------------------------------------------------
      subroutine findline(lnode,ifoundline)
      use SIZE
      integer lnode(2),ifoundline,iline
      integer linenode(2),imatch

      ifoundline = 0
      iline= 0
      imatch = 0
! loop over all quad to find the quad has the fnode numbers.

      do iline = 1,totalLine

         do inode = 1,2
          linenode(inode) = line_array(inode+2,iline)
         enddo

         call iflinematch(imatch,lnode,linenode)

         if(imatch.eq.1) then
          ifoundline = iline
          return
         endif
      enddo

      return
      end
!--------------------------------------------------------------------
      subroutine iflinematch(imatch,lnode,linenode)
      integer lnode(2),linenode(2)
      integer imatch

      imatch = 0
  
      do ilnode = 1,2
         do ilinenode = 1,2
         if(lnode(ilnode).eq.linenode(ilinenode)) imatch=imatch+1
         enddo
      enddo
 
! imatch should equal 2 if lnode and linenode is matching.
      if(imatch.eq.2) then
       imatch = 1
      else
       imatch =0
      endif
 
      return
      end
!-----------------------------------------------------------------
      subroutine r_or_l_detect(quad,rfflag)
! detect if iquad is right-hand or left-hand, elements
      use SIZE
      integer quad, rfflag
      integer node(4)
      real vec12(3),vec14(3),cz
	  
      do inode = 1,4
         node(inode) = quad_array(inode+2,quad)
      enddo

      do i = 1,3
	     vec12(i) = node_xyz(i,node(2)) - node_xyz(i,node(1))
	     vec14(i) = node_xyz(i,node(4)) - node_xyz(i,node(1))
      enddo
	  
      cz = vec12(1)*vec14(2) - vec12(2)*vec14(1)

      if(cz.gt.0.0) rfflag = 0 ! right hand element
      if(cz.lt.0.0) rfflag = 1 ! left hand element

      return
      end
!--------------------------------------------------------------------
      subroutine findquad(fnode,ifoundquad)
      use SIZE
      integer fnode(4),ifoundquad,iquad
      integer quadnode(4),imatch

      ifoundquad = 0
      iquad = 0
      imatch = 0
! loop over all quad to find the quad has the fnode numbers.

      do iquad = 1,totalQuad

         do inode = 1,4
          quadnode(inode) = quad_array(inode+2,iquad)
         enddo

         call ifquadmatch(imatch,fnode,quadnode)

         if(imatch.eq.1) then
          ifoundquad = iquad
          return
         endif
      enddo

      return
      end
!--------------------------------------------------------------------
      subroutine ifquadmatch(imatch,fnode,quadnode)
      integer fnode(4),quadnode(4)
      integer imatch,imatch1,imatch2

      imatch = 0
      imatch1 = 0
      do ifnode = 1,4
         do iquadnode = 1,4
         if(fnode(ifnode).eq.quadnode(iquadnode)) imatch1=imatch1+1
         enddo
      enddo
 
! imatch should equal 4 if fnode and quadnode is matching.
      if(imatch1.eq.4) then
       imatch1 = 1
      else
       imatch1 =0
      endif
 
      imatch2 = 0
      do iquadnode = 1,4
         do ifnode = 1,4
         if(fnode(ifnode).eq.quadnode(iquadnode)) imatch2=imatch2+1
         enddo
      enddo
 
! imatch should equal 4 if fnode and quadnode is matching.
      if(imatch2.eq.4) then
       imatch2 = 1
      else
       imatch2 =0
      endif
	  
       imatch = imatch1*imatch2
 
      return
      end
!-----------------------------------------------------------------
