!--------------------------------------------------------------------
      subroutine set_periodicity(f_or_s)
! handle 2d 3d mesh together
! set periodicity, using bucket sorting.
! 
      use SIZE
	 
      integer quad_edge_node(2,4)
      data quad_edge_node /1,3,3,9,7,9,1,7/
	 
      integer hex_face_node(4,6)
      data hex_face_node /1,3,21,19,3,9,27,21,7,9,27,25,1,7,25,19,1,7,9,3,19,21,27,25/
      
      character*3 ubc
      integer tags(2),ibc,nbc,io
      integer ip,np,ipe,ipe2,nipe(2)
      integer ptags(2)
      integer fnode(4)
      integer mappingOption,eOffset
      real pvec(3)
      real fpxyz(3,2)
      real AB_v(3),AD_v(3),farea,product_v(3)
      real dist,distMax,ptol

      real,save,allocatable,dimension(:,:) ::  fc1,fc2
      integer,save,allocatable,dimension(:) ::  findex1,findex2	  ! ipe->index
      integer,save,allocatable,dimension(:) ::  cfindex2
	  integer,save,allocatable,dimension(:,:) :: rfindex2 	      ! index->ipe
	  integer,save,allocatable,dimension(:,:,:)   :: parray

      integer nseg(3),maxEinBucket,ix,iy,iz,index1,ipe2_bucket
      
      real xmin,xmax,ymin,ymax,zmin,zmax
      real xdiff,ydiff,zdiff,maxDiff

      real ssa(2) ! total ss area
      real ssc(3,2) ! ss center xyz
      real quadNode(3,4) ! quad node xyz
      real quadArea,length
      integer adjacent_ibucket(27),n,quickCheckPass

      integer f_or_s,ih_start,ih_end

! boundary condition summary
!      write(6,*) '******************************************************'
!      write(6,*) 'Boundary info summary'
!      write(6,*) 'sideSet ID'
!      do ibc= 1,bcNumber
!      write(6,*) bcID(ibc)
!      enddo
!      write(6,*) '******************************************************'
 
! boundary condition summary
!      write(6,*) '******************************************************'
!      write(6,*) 'Boundary info summary'
!      write(6,*) 'BoundaryName     BoundaryID'
!      do ibc= 1,bcNumber
!      write(6,*) trim(bcChar(ibc)),bcID(1,ibc)
!      enddo
!      write(6,*) '******************************************************'
 
 
      if (f_or_s.eq.1) then ! for fluid 
        write(6,*) 'For Fluid domain'
        write(6,*) 'Enter number of periodic boundary surface pairs:'
        read (5,*) nbc
      elseif(f_or_s.eq.2) then ! for solid 
        write(6,*) 'For Solid domain'
        write(6,*) 'Enter number of periodic boundary surface pairs:'
        read (5,*) nbc
      endif

       !if(nbc.ne.0) then
      !write(6,*) 'Enter relative search tolerance (recommend 1e-3 for coarse mesh, 1e-5 for refine mesh):'
      !read (5,*) ptol
      !endif
      ptol = 1e-2 !1e-3	  

      if(nbc.le.0) return
	  
      allocate ( parray (2,2,num_elem))

      do ibc = 1,nbc 
        write(6,*) 'input surface 1 and  surface 2  Boundary ID'
        read (5,*) ptags(1),ptags(2)


          ih_start = 1
          ih_end = num_elem

          if (f_or_s.eq.1) ih_end = eftot 
          if (f_or_s.eq.2) ih_start = eftot+1
		
          ipe = 0
          do ihex = ih_start, ih_end
            do iface = 1,2*num_dim
               if(bc(5,iface,ihex).eq.ptags(1)) then   
                ipe = ipe + 1
                parray(1,1,ipe) = ihex
                parray(2,1,ipe) = iface
               endif
            enddo
          enddo
          nipe(1) = ipe
	  
          ipe = 0
          do ihex = ih_start, ih_end
            do iface = 1,2*num_dim
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
            write(6,*) 'ERROR, face numbers are not matching'
            call exitt()
          endif

        mappingOption = 0
        !write(6,*) 'Enter mapping option'
        !write(6,*) '(1 for general mapping (N^2),2 for advanced mapping (N), 3 for uniform element order offset(N) ):'
        !read (5,*) mappingOption
    
        mappingOption = 2  ! only use advanced mapping
 
        if (mappingOption.eq.1) then

        write(6,*) 'general mapping(could be slow for super large mesh):'
        pvec(1) = 0.0
        pvec(2) = 0.0
        pvec(3) = 0.0

        write(6,*) 'input translation vector (surface 1 -> surface 2)'
        if (num_dim.eq.2) then
          read (5,*) pvec(1),pvec(2)
        else
          read (5,*) pvec(1),pvec(2),pvec(3)
        endif 
		
        else if (mappingOption.eq.2) then

        ! write(6,*) 'advanced mapping using bucket sorting (only for 3D mesh):'
        pvec(1) = 0.0
        pvec(2) = 0.0
        pvec(3) = 0.0

        !write(6,*) 'input translation vector (surface 1 -> surface 2)'
        !if (num_dim.eq.2) then
        !  read (5,*) pvec(1),pvec(2)
        !else
        !  read (5,*) pvec(1),pvec(2),pvec(3)
        !endif 

        else if (mappingOption.eq.3) then

        write(6,*) 'only works for extruded hex20 mesh'
        write(6,*) 'uniform element order offset approach assuming element on surface 1 + offset will find element on surface 2'
        write(6,*) 'please input uniform element order offset:'
        read(5,*) eOffset
		
        else 

        write(6,*) 'invalid option, periodicity not done properly'
        call exitt()

        endif

!================================================================================
       allocate (fc1 (3,nipe(1)))
       allocate (fc2 (3,nipe(1)))
!================================================================================
! calculate face center array fc1, fc2
! calculate sideset center xyz..... inorder to calculate translation vector...
! get sideset area, and area weighted center xyz...

       call rzero(ssa,2)
       call rzero(ssc,6)

       do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
         fpxyz(1,1) = 0.0
         fpxyz(2,1) = 0.0
         fpxyz(3,1) = 0.0

         if (num_dim.eq.2) then
           do ifnode = 1,2
             fnode(ifnode)=quad_edge_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.5
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.5
           enddo

           length = ((xm1(fnode(1),1,1,ihex)-xm1(fnode(2),1,1,ihex))**2.0+(ym1(fnode(1),1,1,ihex)-ym1(fnode(2),1,1,ihex))**2.0)**0.5 
           ssa(1) = ssa(1) + length
           ssc(1,1) = ssc(1,1) + fpxyz(1,1)*length
           ssc(2,1) = ssc(2,1) + fpxyz(2,1)*length

         else
           do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25

             
            quadNode(1,ifnode) = xm1(fnode(ifnode),1,1,ihex)
            quadNode(2,ifnode) = ym1(fnode(ifnode),1,1,ihex)
            quadNode(3,ifnode) = zm1(fnode(ifnode),1,1,ihex)

           enddo

           call get_quadArea(quadNode,quadArea)

           ssa(1) = ssa(1) + quadArea
           ssc(1,1) = ssc(1,1) + fpxyz(1,1)*quadArea
           ssc(2,1) = ssc(2,1) + fpxyz(2,1)*quadArea
           ssc(3,1) = ssc(3,1) + fpxyz(3,1)*quadArea         

         endif

         fc1(1,ipe) = fpxyz(1,1)
         fc1(2,ipe) = fpxyz(2,1)
         fc1(3,ipe) = fpxyz(3,1)
       enddo

       ssc(1,1) = ssc(1,1)/ssa(1)  ! get avg sideset center
       ssc(2,1) = ssc(2,1)/ssa(1) 
       ssc(3,1) = ssc(3,1)/ssa(1)   


       do ipe = 1,nipe(1)
         ihex = parray(1,2,ipe)
         iface = parray(2,2,ipe)
         fpxyz(1,1) = 0.0
         fpxyz(2,1) = 0.0
         fpxyz(3,1) = 0.0

         if (num_dim.eq.2) then
           do ifnode = 1,2
             fnode(ifnode)=quad_edge_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.5
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.5
           enddo

           length = ((xm1(fnode(1),1,1,ihex)-xm1(fnode(2),1,1,ihex))**2.0+(ym1(fnode(1),1,1,ihex)-ym1(fnode(2),1,1,ihex))**2.0)**0.5 
           ssa(2) = ssa(2) + length
           ssc(1,2) = ssc(1,2) + fpxyz(1,1)*length
           ssc(2,2) = ssc(2,2) + fpxyz(2,1)*length
    
         else
           do ifnode = 1,4
             fnode(ifnode)=hex_face_node(ifnode,iface)
             fpxyz(1,1) = fpxyz(1,1)+xm1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(2,1) = fpxyz(2,1)+ym1(fnode(ifnode),1,1,ihex)*0.25
             fpxyz(3,1) = fpxyz(3,1)+zm1(fnode(ifnode),1,1,ihex)*0.25

            quadNode(1,ifnode) = xm1(fnode(ifnode),1,1,ihex)
            quadNode(2,ifnode) = ym1(fnode(ifnode),1,1,ihex)
            quadNode(3,ifnode) = zm1(fnode(ifnode),1,1,ihex)

           enddo

           call get_quadArea(quadNode,quadArea)

           ssa(2) = ssa(2) + quadArea
           ssc(1,2) = ssc(1,2)  + fpxyz(1,1)*quadArea
           ssc(2,2) = ssc(2,2)  + fpxyz(2,1)*quadArea
           ssc(3,2) = ssc(3,2)  + fpxyz(3,1)*quadArea         

         endif

         fc2(1,ipe) = fpxyz(1,1)
         fc2(2,ipe) = fpxyz(2,1)
         fc2(3,ipe) = fpxyz(3,1)
       enddo

       ssc(1,2) = ssc(1,2)/ssa(2)  ! get avg sideset center
       ssc(2,2) = ssc(2,2)/ssa(2) 
       ssc(3,2) = ssc(3,2)/ssa(2)   

! calculate translate vector
! based on area-weighted sideset center...
        pvec(1) = ssc(1,2) - ssc(1,1)
        pvec(2) = ssc(2,2) - ssc(2,1)
        pvec(3) = ssc(3,2) - ssc(3,1)

        write (6,*) 'translation vector: ',pvec(1),' ',pvec(2),' ',pvec(3)


!================================================================================
!================================================================================
! preprocess fc1 and fc2 for option 1 and 2
       if ((mappingOption.eq.1) .or. (mappingOption.eq.2)) then

        write(6,*) 'offset and normalize periodic face coordinates'

! offset fc1 to match fc2	   
       do ipe = 1,nipe(1)
         fc1(1,ipe) =  fc1(1,ipe) + pvec(1)
         fc1(2,ipe) =  fc1(2,ipe) + pvec(2)
         fc1(3,ipe) =  fc1(3,ipe) + pvec(3)
       enddo

! normalize fc1 and fc2 
       xmax = -1e6
       ymax = -1e6
       zmax = -1e6
       xmin = 1e6
       ymin = 1e6
       zmin = 1e6

       do ipe = 1,nipe(1)
         xmax = max(fc1(1,ipe),xmax)
         ymax = max(fc1(2,ipe),ymax)
         zmax = max(fc1(3,ipe),zmax)
         xmin = min(fc1(1,ipe),xmin)
         ymin = min(fc1(2,ipe),ymin)
         zmin = min(fc1(3,ipe),zmin)
       enddo

       xdiff = xmax - xmin
       ydiff = ymax - ymin
       zdiff = zmax - zmin
	   
       maxDiff = max(xdiff,ydiff,zdiff)

       if (xdiff.gt.(0.01*maxDiff)) then
          do ipe = 1,nipe(1)
           fc1(1,ipe) = (fc1(1,ipe)-xmin)/xdiff
           fc2(1,ipe) = (fc2(1,ipe)-xmin)/xdiff
          enddo
       endif

       if (ydiff.gt.(0.01*maxDiff)) then
          do ipe = 1,nipe(1)
           fc1(2,ipe) = (fc1(2,ipe)-ymin)/ydiff
           fc2(2,ipe) = (fc2(2,ipe)-ymin)/ydiff
          enddo
       endif

       if (zdiff.gt.(0.01*maxDiff)) then
          do ipe = 1,nipe(1)
           fc1(3,ipe) = (fc1(3,ipe)-zmin)/zdiff
           fc2(3,ipe) = (fc2(3,ipe)-zmin)/zdiff
          enddo
       endif

       endif

!================================================================================
! deal with option 2.
! use bucket sorting to make index of each point.
! 
      if (mappingOption.eq.2) then

        write(6,*) 'bucket sorting to make index of each point'

! number of buckets in each direction 
       if (nipe(1).le.1e4 )then
       nseg(1) = 10
       nseg(2) = 10
       nseg(3) = 10
       else if ((nipe(1).gt.1e4 ).and.(nipe(1).le.1e6)) then
       nseg(1) = 20
       nseg(2) = 20
       nseg(3) = 20
       else if (nipe(1).gt.1e6) then
       nseg(1) = 50
       nseg(2) = 50
       nseg(3) = 50
       endif

       if (xdiff.le.(0.01*maxDiff)) then
        nseg(1) = 1
       endif

       if (ydiff.le.(0.01*maxDiff)) then
        nseg(2) = 1
       endif

       if (zdiff.le.(0.01*maxDiff)) then
        nseg(3) = 1
       endif

       allocate (findex1 (nipe(1)))
       allocate (findex2 (nipe(1)))

       call rzero_int2(findex1,nipe(1))
       call rzero_int2(findex2,nipe(1))

       allocate (cfindex2 (nseg(1)*nseg(2)*nseg(3)))
       call rzero_int2(cfindex2,nseg(1)*nseg(2)*nseg(3))

       do ipe = 1,nipe(1)
	 
! sorting elements into buckets	 
         ix =  int(fc1(1,ipe)*dble(nseg(1)-1))+1
         iy =  int(fc1(2,ipe)*dble(nseg(2)-1))+1
         iz =  int(fc1(3,ipe)*dble(nseg(3)-1))+1

         findex1(ipe) = ix + (iy-1)*nseg(1) + (iz-1)*nseg(1)*nseg(2)

         ix =  int(fc2(1,ipe)*dble(nseg(1)-1))+1
         iy =  int(fc2(2,ipe)*dble(nseg(2)-1))+1
         iz =  int(fc2(3,ipe)*dble(nseg(3)-1))+1

         findex2(ipe) = ix + (iy-1)*nseg(1) + (iz-1)*nseg(1)*nseg(2)

         cfindex2(findex2(ipe)) = cfindex2(findex2(ipe))+1

       enddo

       ! maxEinBucket is the max number of point in all bucket
       maxEinBucket = maxval(cfindex2)
       ! write(6,*) 'maxEinBucket=',maxEinBucket

       allocate (rfindex2 (maxEinBucket,nseg(1)*nseg(2)*nseg(3)))
       call rzero_int2(rfindex2,maxEinBucket*nseg(1)*nseg(2)*nseg(3))
       call rzero_int2(cfindex2,nseg(1)*nseg(2)*nseg(3))

       do ipe = 1,nipe(1)
          cfindex2(findex2(ipe)) = cfindex2(findex2(ipe))+1 ! is the counter of each bucket
          rfindex2(cfindex2(findex2(ipe)),findex2(ipe))=ipe ! store ipe in this bucket(findex2(ipe))
       enddo

! now, loop face1, and search in this bucket only, this saves a lot of time .... 
! to find the closest point.
     do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
         
		 ! index1 = findex1(ipe) ! bucket number

         n = 27
         call rzero_int2(adjacent_ibucket,n)

! get adjacent bucket number.
         ix =  int(fc1(1,ipe)*dble(nseg(1)-1))+1
         iy =  int(fc1(2,ipe)*dble(nseg(2)-1))+1
         iz =  int(fc1(3,ipe)*dble(nseg(3)-1))+1 

         ib = 0 ! for all 27 buckets..
         do k = 1,3
         do j = 1,3
         do i = 1,3
         
         ix1 = ix-2 + i
         iy1 = iy-2 + j
         iz1 = iz-2 + k
         
         ibucket = ix1 + (iy1-1)*nseg(1) + (iz1-1)*nseg(1)*nseg(2)
         
         if ((ix1.lt.1).or.(ix1.gt.nseg(1))) ibucket = 0
         if ((iy1.lt.1).or.(iy1.gt.nseg(2)))  ibucket = 0
         if ((iz1.lt.1).or.(iz1.gt.nseg(3)))  ibucket = 0
         
         ib = ib + 1
         adjacent_ibucket(ib) = ibucket
         
         enddo
         enddo
         enddo


         distMax = ptol

         do ib = 1,27
         ibucket = adjacent_ibucket(ib)        ! loop in all adjacent bucket
         if (ibucket.gt.0) then

         do ipe2_bucket = 1,cfindex2(ibucket)  ! loop in this bucket

          ipe2 = rfindex2(ipe2_bucket,ibucket) ! counter, and bucket number 
          ihex2 = parray(1,2,ipe2)
          iface2 = parray(2,2,ipe2)

        dist = sqrt((fc1(1,ipe)-fc2(1,ipe2))**2 &
      +(fc1(2,ipe)-fc2(2,ipe2))**2 &
      +(fc1(3,ipe)-fc2(3,ipe2))**2)

          if (dist.lt.distMax) then 
             distMax = dist
             bc(1,iface,ihex) = ihex2*1.0
             bc(2,iface,ihex) = iface2*1.0
          endif

         enddo

         endif
         enddo
		 

      enddo

      deallocate (findex1)
      deallocate (findex2)
      deallocate(cfindex2)
      deallocate(rfindex2)

      endif
!================================================================================
! only deal with option 1 and 3 here
      if ((mappingOption.eq.1) .or. (mappingOption.eq.3)) then

! 1st loop, loop faces on surface 1
      do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)

        if (mappingOption.eq.1) then ! general mapping 
! 2nd loop over surface 2
         distMax = ptol
         do ipe2 = 1,nipe(2)
         ihex2 = parray(1,2,ipe2)
         iface2 = parray(2,2,ipe2)

        dist = sqrt((fc1(1,ipe)-fc2(1,ipe2))**2 &
      +(fc1(2,ipe)-fc2(2,ipe2))**2 &
      +(fc1(2,ipe)-fc2(2,ipe2))**2)

               if (dist.lt.distMax) then 
                  distMax = dist
                  bc(1,iface,ihex) = ihex2*1.0
                  bc(2,iface,ihex) = iface2*1.0
               endif
         enddo
		 
         else if (mappingOption.eq.3) then ! uniform element order offset
		 
            ihex2 = ihex + eOffset
            iface2 = -1
            do ifc2= 1,2*num_dim
               if(bc(5,ifc2,ihex2).eq.ptags(2)) then
                 iface2 = ifc2
               endif
            enddo
            if (iface2.eq.-1) then 
            write(6,*) 'ERROR: the offset element does not have face of sideSet: ',ptags(2)
            call exitt()
            endif

            bc(1,iface,ihex) = ihex2*1.0
            bc(2,iface,ihex) = iface2*1.0
            bc(1,iface2,ihex2) = ihex*1.0
            bc(2,iface2,ihex2) = iface*1.0

        endif
  
      enddo
      endif
!================================================================================
! add a quick check here
     quickCheckPass = 0
     do ipe = 1,nipe(1)
         ihex = parray(1,1,ipe)
         iface = parray(2,1,ipe)
         ihex2 = int(bc(1,iface,ihex))
         iface2 = int(bc(2,iface,ihex)) 
         if ((ihex2.eq.0).or.(iface2.eq.0)) then
          quickCheckPass = 1
         endif
     enddo
     if (quickCheckPass.eq.1) then
      write(6,*) 'Error, quick check failed, please check your periodic faces'
     endif
!================================================================================
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
 
      deallocate (fc1)
      deallocate (fc2)
 
 
      enddo

      deallocate(parray)

      write(6,*) '******************************************************'
      write(6,*) 'Please set boundary conditions to all non-periodic boundaries'
      write(6,*) 'in .usr file usrdat2() subroutine'
      write(6,*) '******************************************************'

      return
      end
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      subroutine  get_quadArea(quadNode,quadArea)

! get quad area.
! 1st, get center xyz,
! 2nd, get area of 4 triangles...

      real quadNode(3,4),quadArea
      real*8 triArea
      real center(3)
      real*8 vec1(1),vec2(3),prod_v(3)
  
      quadArea = 0.0
      call rzero(center,3)

      do i = 1,4
      center(1) = center(1)+quadNode(1,i)/4.0
      center(2) = center(2)+quadNode(2,i)/4.0
      center(3) = center(3)+quadNode(3,i)/4.0
      enddo

      do i = 1,4
      i1 = i
      i2 = i+1
      if (i2>4) i2 = 1 

      vec1(1) =  quadNode(1,i1) - center(1)
      vec1(2) =  quadNode(2,i1) - center(2)
      vec1(3) =  quadNode(3,i1) - center(3)

      vec2(1) =  quadNode(1,i2) - center(1)
      vec2(2) =  quadNode(2,i2) - center(2)
      vec2(3) =  quadNode(3,i2) - center(3)

      call cross_product(vec1,vec2,prod_v,triArea)
      quadArea = quadArea + triArea
      enddo


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
