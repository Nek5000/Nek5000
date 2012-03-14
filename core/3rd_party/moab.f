#ifdef PTRSIZE8
#define IMESH_HANDLE integer*8
#else
#define IMESH_HANDLE integer*4
#endif

#define IMESH_NULL 0
#define IMESH_ASSERT \
      if (ierr.ne.0) call imesh_err(ierr,imeshh,'moab.f ',__LINE__)

#define IMESH_NULLSTRIP(s) s(:index(s, char(0))-1)

#define MYLOC LOC


c    We are targeting hex27 meshes coming from MOAB.
c    By default we read the MOAB mesh from 'input.h5m'
c
c    fun stuff: 
c    * "imesh" is already used in TSTEP
c    * non-standard extensions: %val and loc + Cray pointers
c
c    unsupported for now:
c    * scalable mesh read (every processor will read the mesh)
c    * conj. heat transfer
c    * periodic BCs
c-----------------------------------------------------------------------
      subroutine moab_dat
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'
      include 'GEOM'
      common /mbc/ moabbc(6,lelt) 
      integer moabbc

      call nekMOAB_load                   ! read mesh using MOAB
      call nekMOAB_get_elems              ! read material sets and establish mapping

      call chk_nel

      call nekMOAB_create_tags             ! allocate MOAB tags to reflect Nek variables

      call mapelpr                        ! create gllel mapping 
      call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
      call xml2xc                         ! fill xc,yc,zc

      call ifill(moabbc, -1, 6*lelt)
      call nekMOAB_BC(moabbc)             ! read MOAB BCs 
      
c      call nekMOAB_compute_diagnostics

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_create_tags
c
      implicit none
#include "NEKMOAB"

      integer ierr, ntot
      integer semdim(3)
      iBase_TagHandle tagh

      ntot = nx1*ny1*nz1

      call iMesh_createTagWithOptions(%VAL(imeshh), "SEM_X",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), xm1Tag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "SEM_Y",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), ym1Tag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "SEM_Z",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), zm1Tag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "VX",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), vxTag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "VY",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), vyTag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "VZ",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), vzTag, ierr)
      IMESH_ASSERT

      call iMesh_createTagWithOptions(%VAL(imeshh), "TEMP",
     1     "moab:TAG_STORAGE_TYPE=DENSE ", 
     1     %VAL(ntot), %VAL(iBase_DOUBLE), tTag, ierr)
      IMESH_ASSERT

      if (nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) then
         call iMesh_createTagWithOptions(%VAL(imeshh), "PRESS",
     1        "moab:TAG_STORAGE_TYPE=DENSE ", 
     1        %VAL(ntot), %VAL(iBase_DOUBLE), pTag, ierr)
         IMESH_ASSERT
      endif

c create a tag to store SEM dimensions, and set it on the root set
      semdim(1) = nx1
      semdim(2) = ny1
      semdim(3) = nz1
      call iMesh_createTagWithOptions(%VAL(imeshh), "SEM_DIMS", 
     1     "moab:TAG_STORAGE_TYPE=SPARSE", 
     1     %VAL(3), %VAL(iBase_INTEGER), tagh, ierr)
      IMESH_ASSERT
      call iMesh_setEntSetData(%VAL(imeshh), %VAL(rootset), %VAL(tagh), 
     1     semdim, 12, ierr)
      IMESH_ASSERT

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_load
c
c     Load "filename" into imesh/moab, store imesh handle 
c     in /nekmoab/ common block
c
      implicit none
#include "NEKMOAB"
      include 'mpif.h'
c two forms of load options, depending whether we\'re running serial or parallel      
      character*(*) parLoadOpt, serLoadOpt
      parameter(parLoadOpt=" moab:PARALLEL=READ_PART   moab:PARTITION=PA
     $RALLEL_PARTITION moab:PARALLEL_RESOLVE_SHARED_ENTS moab:PARTITION_
     $DISTRIBUTE moab:CPUTIME")
      parameter(serLoadOpt = " ")
      common /nekmpi/ nid_,np,nekcomm,nekgroup,nekreal
      integer nekcomm, nekgroup, nekreal, nid_, np

      integer ierr
      IBASE_HANDLE_T ccomm

c      !Initialize imesh and load file
      imeshh = IMESH_NULL
      call iMesh_newMesh(" ", imeshh, ierr)
      IMESH_ASSERT

      if (1 .lt. np) then
         call moab_comm_f2c(nekcomm, ccomm)
         call iMeshP_createPartitionAll(%VAL(imeshh), 
     1        %VAL(ccomm),
     1        hPartn, ierr)
         IMESH_ASSERT

         call iMesh_getRootSet(%VAL(imeshh), rootset, ierr)

         call iMeshP_loadAll(%VAL(imeshh), %VAL(hPartn),%VAL(rootset), 
     $        H5MFLE, parLoadOpt, ierr)
         IMESH_ASSERT

         partsSize = 0
         rpParts = IMESH_NULL
         call iMeshP_getLocalParts(%VAL(imeshh), %VAL(hPartn), 
     1        rpParts, partsSize, partsSize, ierr)

      else
         call iMesh_getRootSet(%VAL(imeshh), rootset, ierr)

         call iMesh_load(%VAL(imeshh), %VAL(rootset), 
     $        H5MFLE, parLoadOpt, ierr)
         IMESH_ASSERT
      endif

c initialize tag handles
      globalIdTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "GLOBAL_ID",       !/*in*/ const char* tag_name,
     $     globalIdTag,       !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT
      
      matsetTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "MATERIAL_SET", !/*in*/ const char* tag_name,
     $     matSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT

      neusetTag = 0
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "NEUMANN_SET", !/*in*/ const char* tag_name,
     $     neuSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT

      return
      end  
c-----------------------------------------------------------------------
      subroutine nekMOAB_get_elems()
c
c get fluid/solid elements and establish global id mapping
c
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'
      include 'SCRCT'
      include 'ZPER'

      include 'mpif.h'

      integer iglsum, i, ierr
      iBase_EntitySetHandle dumsets(numsts)
      IBASE_HANDLE_T valptr, setsptr
      integer dumval, dumalloc, dumsize, dumnum, ns, ilast

c get fluid, other material sets, and count elements in them
      valptr = loc(dumval)
      setsptr = loc(dumsets(1))
      dumalloc = numsts
      ilast = 0
      do i = 1, numflu+numoth
         dumval = matids(i)
         dumsize = numsts
c get the set by matset number
         call iMesh_getEntSetsByTagsRec(%VAL(imeshh), %VAL(rootset),
     $        matsetTag, valptr, %VAL(1), %VAL(1), 
     $        setsptr, dumsize, dumsize, ierr)
         if (dumsize .gt. 1) then
            call exitti('More than one material set with id ', dumval)
         elseif (dumsize .eq. 0) then
            call exitti('Could not find material set with id ', dumval)
         endif
         IMESH_ASSERT
c get the number of hexes
         call iMesh_getNumOfTopoRec(%VAL(imeshh), %VAL(dumsets(1)), 
     $        %VAL(iMesh_HEXAHEDRON), %VAL(1), dumnum, ierr)
         IMESH_ASSERT
         matsets(i) = dumsets(1)
         iestart(i) = ilast + 1
         ilast = ilast + dumnum
         iecount(i) = dumnum
c get an iterator for this set, used later
         call iMesh_initEntArrIterRec(%VAL(imeshh), %VAL(dumsets(1)),
     $        %VAL(iBase_REGION), %VAL(iMesh_HEXAHEDRON),
     $        %VAL(dumnum), %VAL(0), %VAL(1), ieiter(i), ierr)
c set total number if nec
         if (i .eq. numflu) then
            nelv = ilast
         endif
c this is if, not elseif, to handle numoth=0
         if (i .eq. numflu+numoth) then
            nelt = ilast
         endif
      enddo

c set remaining values to default values
      do i = numflu+numoth+1, numsts
         iecount(i) = -1
         iestart(i) = -1
         ieiter(i) = 0
         matsets(i) = 0
      enddo

c check local size
      if(nelv .le. 0 .or. nelv .gt. lelv) then
         print *, 'ABORT: nelv is invalid in nekmoab_proc_map'
         print *, 'nelv, lelv = ', nelv, lelv
         call exitt
      endif

c reduce to get global numbers of fluid, other elements, and check size
      nelgv = iglsum(nelv,1)
      nelgt = iglsum(nelt,1)
      if (NELGT .gt. LELT) then
         print *, 'ABORT; increase lelg ',nelgv,lelg
         call exitt
      endif

c assign GLLNID, map from gid to proc
      call izero(GLLNID, NELGV)
      do i = 1, numflu+numoth
         call nekMOAB_gllnid(matsets(i), ieiter(i), iecount(i))
      enddo

      call igop(GLLNID, GLLEL, '+  ', NELGV)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_gllnid(matset, iter, count)
      implicit none
c
c initialize global ids for hexes in this set
      iBase_EntitySetHandle matset
      iBase_EntityArrIterator iter
      integer itmp, i, j, ierr, gid, count, atend
      IBASE_HANDLE_T tag_ptr
      pointer(tag_ptr, gid(1))

#include "NEKMOAB"
      include 'PARALLEL'
c 
      call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(iter), ierr)
      IMESH_ASSERT
      i = 0
      atend = 0
      do while (atend .eq. 0)
         call iMesh_tagIterate(%VAL(imeshh), %VAL(globalIdTag), 
     $        %VAL(iter), tag_ptr, itmp, ierr)
         IMESH_ASSERT

c set the global ids to this proc
         do j = 1, itmp
            if (gid(j) .gt. nelgt) 
     $           call exitti('Global id greater than NELGT', gid(i))
            gllnid(gid(j)) = nid
         enddo
c step the iterator
         call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(iter), %VAL(itmp), 
     $        atend, ierr)
         IMESH_ASSERT
         i = i + itmp
      enddo

c assert we got the right number of elements
      if (i .lt. count) call exitti(
     $     'Wrong number of entities in region iterator', i)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_loadConn(vertex, nelgt, ncrnr)
c
c     Fill the vertex array with connectivity data from imesh
c
c     vertex(ncrnr, nelt): int array, global id of vertices in element nelt
c     nelgt: int, global number of elements
c     ncrnr: int, number of corner vertices per element (should be 8)
c
      implicit none
#include "NEKMOAB"

      integer vertex(ncrnr, *), i

c
c get corner vertex gids
      integer e_in_set, eid, j, k, nv, ierr, e_in_chunk, v_per_e
      integer gids(27)
      iBase_EntityArrIterator iter
      IBASE_HANDLE_T connect_ptr
      iBase_EntityHandle connect
      pointer (connect_ptr, connect(0:1))

      integer l2c(8)
      save    l2c
      data    l2c / 1, 2, 4, 3, 5, 6, 8, 7 /

      do i = 1, numflu+numoth
         call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), ierr)
         IMESH_ASSERT
         nv = 8
         e_in_set = 0
         eid = iestart(i)
         do while (e_in_set .lt. iecount(i))
c     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), 
     $           connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

c     for each element
            do j = 0, e_in_chunk-1
c     get vertex gids for this e
               call iMesh_getIntArrData(%VAL(imeshh), !iMesh_Instance instance,
     $              connect(j*v_per_e), %VAL(8), %VAL(globalIdTag), 
     $              loc(gids), nv, nv, ierr)
               IMESH_ASSERT
c     permute into vertex array
               do k=1, 8
                  vertex(k, eid) = gids(l2c(k))
               enddo
               eid = eid + 1
            enddo

            e_in_set = e_in_set + e_in_chunk
         enddo
      enddo

      return
      end 
c-----------------------------------------------------------------------------
      subroutine nekMOAB_loadCoord(x27, y27, z27)
c     
c     stuff the xyz coords of the 27 verts of each local element -- 
c     shared vertex coords are stored redundantly
c     
      implicit none
#include "NEKMOAB"
      real x27(27,*), y27(27,*), z27(27,*)
      IBASE_HANDLE_T connect_ptr
      iBase_EntityHandle connect
      pointer(connect_ptr, connect(0:1))
      integer i, j, k, ierr, e_in_chunk, e_in_set, v_per_e

      do i = 1, numflu+numoth
         call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), ierr)
         IMESH_ASSERT
         e_in_set = 0
         do while (e_in_set .lt. iecount(i))
c     get ptr to connectivity for this chunk
            call iMesh_connectIterate(%VAL(imeshh), %VAL(ieiter(i)), 
     $           connect_ptr, v_per_e, e_in_chunk, ierr)
            IMESH_ASSERT

c     for each element
            do j = 0, e_in_chunk-1
c     get vertex gids for this e
               do k = 1, TWENTYSEVEN
                  call iMesh_getVtxCoord(%VAL(imeshh),
     $                 %VAL(connect(j*v_per_e+k-1)), 
     $                 x27(k,j+1), y27(k,j+1), z27(k,j+1), ierr)
             IMESH_ASSERT
               enddo
            enddo

            e_in_set = e_in_set + e_in_chunk
         enddo
      enddo

      return
      end 
c-----------------------------------------------------------------------------
      subroutine nekMOAB_BC(moabbc)
c
c     Copy the boundary condition information into bcdata array
c

      implicit none
#include "NEKMOAB"
      integer moabbc(6,lelt)

      IBASE_HANDLE_T hentSet(*)
      pointer (rpentSet, hentSet)
      integer entSetSize

      integer ierr, i, tagIntData

      !Sidesets in cubit come in as entity sets with the NEUMANN_SET -- see sample file
      call iMesh_getTagHandle(%VAL(imeshh),
     $     "NEUMANN_SET", neuSetTag, ierr)
      IMESH_ASSERT

      rpentSet = IMESH_NULL
      entSetSize      = 0
      call iMesh_getEntSetsByTagsRec(%VAL(imeshh),
     $     %VAL(rootset), neuSetTag, %VAL(IMESH_NULL), %VAL(1), %VAL(0),
     $     rpentSet, entSetSize, entSetSize, ierr)
      IMESH_ASSERT

      !print *, 'H3', entSetSize

      do i=1, entSetSize
         call iMesh_getIntData(%VAL(imeshh),
     $        %VAL(hentSet(i)), %VAL(neuSetTag), tagIntData, ierr)

         if (ierr .eq. 0) then !tag was defined
            !print *, 'H32', tagIntData
            call nekMOAB_intBC(moabbc, hentSet(i), tagIntData)
         endif
      enddo

      call free(rpentSet)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_intBC(bcdata, setHandle, setId)
c
c     Internal function, don\'t call directly
c

      implicit none
#include "NEKMOAB"

      integer bcdata(6,lelt)

      iBase_EntitySetHandle setHandle

      IBASE_HANDLE_T hset
      integer setId !coming from cubit

      IBASE_HANDLE_T faces(*)
      pointer (rpfaces, faces)
      integer facesSize

      IBASE_HANDLE_T ahex(*)
      pointer (rpahex, ahex)
      integer ahexSize

      IBASE_HANDLE_T fvtx(*)
      pointer (rpfvtx, fvtx)
      integer fvtxSize

      IBASE_HANDLE_T hvtx(*)
      pointer (rphvtx, hvtx)
      integer hvtxSize

      integer ierr, i, j, elno

      integer side_no, side_offset, side_sense
      integer hexMBCNType
      integer num_sides
      integer elnos(2)
      integer numids

      numids = 0

      call iMesh_MBCNType(%VAL(iMesh_HEXAHEDRON), hexMBCNType)

      !what faces are in this set?
      facesSize = 0
      rpfaces = IMESH_NULL
      call iMesh_getEntitiesRec(%VAL(imeshh), 
     $     %VAL(setHandle), %VAL(iBase_FACE), %VAL(iMesh_QUADRILATERAL),
     $     %VAL(1), rpfaces, facesSize, facesSize, ierr)
      IMESH_ASSERT

      num_sides = 0

      do i=1, facesSize
         !get vertices defining the face
         fvtxSize = 0
         rpfvtx = IMESH_NULL
         call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), 
     $        %VAL(iBase_VERTEX), rpfvtx, fvtxSize, fvtxSize, ierr)
         IMESH_ASSERT         

         !get hexes adjacent to the face (1 or 2) 
         ! -- hopefully only returns hexes on local proc, but untested
         ahexSize = 0
         rpahex = IMESH_NULL                  
         call iMesh_getEntAdj(%VAL(imeshh), %VAL(faces(i)), 
     $        %VAL(iBase_REGION), rpahex, ahexSize, ahexSize, ierr)
         IMESH_ASSERT

         do j=1, ahexSize                                
            !get verts adjacent to the hex
            hvtxSize = 0
            rphvtx = IMESH_NULL                  
            call iMesh_getEntAdj(%VAL(imeshh), 
     $           %VAL(ahex(j)), %VAL(iBase_VERTEX), 
     $           rphvtx, hvtxSize, hvtxSize, ierr)
            IMESH_ASSERT

            !Get the side number
            call MBCN_SideNumberUlong(%VAL(rphvtx),
     $           %VAL(hexMBCNType), %VAL(rpfvtx), 
     $           %VAL(4), %VAL(2),     !4 vertices, want 2D side number
     $           side_no, side_sense, side_offset) 
           if ((side_no .lt. 0) .or. (side_no .gt. 5)) ierr = 1
           IMESH_ASSERT

           side_no = side_no + 1 !moab side number is zero based
           num_sides = num_sides + 1

           !call nekMOAB_getGlobElNo(ahex(j), elno)
           !print *, 'BLAH', elno, side_no, setId

           call nekMOAB_getElNo(ahex(j), elno)
           if (ahexSize .eq. 2) elnos(j) = elno

           if (bcdata(side_no, elno) .ne. -1) 
     $          print *, 'Warning: resetting BC, bcno, elno, sideno = ', 
     $            setId, elno, side_no 
           bcdata(side_no, elno) = setId
           numids = numids + 1

           call free(rphvtx)
         enddo

         if (ahexSize .eq. 2) 
     $        print *, 'Warning: face shared by 2 hexes: ', elnos(1), 
     $        elnos(2)

         call free(rpahex)
         call free(rpfvtx)

      enddo

      call free(rpfaces)

      print *, 'Setid, numids = ', setId, numids

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getElNo(handle, elno)
c
c     Given the imesh handle of an element (hex), 
c     return its local nek element number in elno
c     Cannot be used until the GLLEL array has been set (in map2.f)
c
      implicit none
#include "NEKMOAB"
      include 'PARALLEL'

      IBASE_HANDLE_T handle
      integer elno

      !first get global id
      call nekMOAB_getGlobElNo(handle, elno)
      !then convert to local id
      elno = GLLEL(elno)

      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getGlobElNo(handle, elno)

      implicit none
#include "NEKMOAB"
      include 'PARALLEL'

      IBASE_HANDLE_T handle
      integer elno
      
      integer ierr

      call iMesh_getIntData(%VAL(imeshh), %VAL(handle), 
     $     %VAL(globalIdTag), elno, ierr)
      IMESH_ASSERT
      
      if (GLLNID(elno) .ne. NID) then
         print *, 'Wrong proc for element; gid, proc, gllnid = ',
     $        elno, nid, gllnid(elno)
         call exitt
      endif
                          
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_bcs() ! fill the nek cbc arrays
      implicit none
c
#include "NEKMOAB"      
      common /mbc/ moabbc(6,lelt)
      integer moabbc

      integer e,f
      character*3 cbi

      integer ibcs(3), i, lcbc, nface
      data ibcs / 0, 0, 0 /

      lcbc=18*lelt*(ldimt1 + 1)
      call blank(cbc,lcbc)

      nface = 2*ndim
      do e=1,nelt
         do f=1,nface
            cbi = 'E  '
            if (moabbc(f,e) .ne. -1) then
               do i = 1, numsts
                  if (ibcsts(i) .eq. moabbc(f,e)) then
                     cbi = bctyps(i)
                     goto 100
                  endif
               enddo
 100           continue
            endif
            do i = 1, nfield
               cbc(f, e, i) = cbi
            enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine moab_geometry (xmlo,ymlo,zmlo)

      implicit none
#include "NEKMOAB"      

      real      xmlo(nx1*ny1*nz1,1)
     $        , ymlo(nx1*ny1*nz1,1)
     $        , zmlo(nx1*ny1*nz1,1)
      integer   e, nmoab

      common /tcrmg/ x27(27,lelt), y27(27,lelt), z27(27,lelt)
      real x27, y27, z27

      call nekMOAB_loadCoord(x27, y27, z27) !     Get coords from moab

c     call outmat(x27,27,8,'x27dat',nelt)
c     call outmat(y27,27,8,'y27dat',nelt)
c     call outmat(z27,27,8,'z27dat',nelt)

      nmoab = 3 !not used
      do e=1,nelt   !  Interpolate for each element
         call permute_and_map(xmlo(1,e),x27(1,e),nx1,nmoab,e)
         call permute_and_map(ymlo(1,e),y27(1,e),ny1,nmoab,e)
         if (if3d) call permute_and_map(zmlo(1,e),z27(1,e),nz1,nmoab,e)
      enddo

c     param(66) = 0
c     ifxyo = .true.
c     ifvo  = .true.
c     call outpost(xm1,ym1,zm1,pr,t,'   ')
c     write(6,*) 'DONE PERMUTE; ABORT'
c     call exitt

      return
      end
c-----------------------------------------------------------------------
      subroutine xml2xc
      
      implicit none
#include "NEKMOAB"
      include 'GEOM'
      integer i, j, k, l

      integer e

      l = 0
      if (if3d) then
         do e=1,nelt
         do k=1,nz1,nz1-1
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,1) = xm1(i,j,k,e)
            yc(l,1) = ym1(i,j,k,e)
            zc(l,1) = zm1(i,j,k,e)
         enddo
         enddo
         enddo
         enddo


      else  ! 2D
         do e=1,nelt
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,1) = xm1(i,j,1,e)
            yc(l,1) = ym1(i,j,1,e)
         enddo
         enddo
         enddo
      endif

      do e=1,nelt ! flip corners back to pre-proc notation
         dtmp = xc(3,e)
         xc(3,e) = xc(4,e)
         xc(4,e) = dtmp
         dtmp = yc(3,e)
         yc(3,e) = yc(4,e)
         yc(4,e) = dtmp
         dtmp = zc(3,e)
         zc(3,e) = zc(4,e)
         zc(4,e) = dtmp

         if(if3d) then
           dtmp = xc(7,e)
           xc(7,e) = xc(8,e)
           xc(8,e) = dtmp
           dtmp = yc(7,e)
           yc(7,e) = yc(8,e)
           yc(8,e) = dtmp
           dtmp = zc(7,e)
           zc(7,e) = zc(8,e)
           zc(8,e) = dtmp
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine permute_and_map(x,x27,nx,nmoab,e)
c
      implicit none
#include "NEKMOAB"
      integer nx, nmoab, e

      real x(1),x27(0:1)

      common /ctmp0/ z3(3),zpt(lx1)
     $             , xt(3,3,3)
     $             , wk(3*lx1*ly1*lz1)
     $             , xw(3*lx1*ly1*lz1)
      real z3, zpt, xt, wk, xw

      real interp(lx1*3),interpt(lx1*3),zl(lx1*3)
      save interp,interpt

      integer moabmap(27)
      save    moabmap
      data    moabmap
     $      /  0,  8,  1, 11, 24,  9,  3, 10,  2  
     $      , 12, 20, 13, 23, 26, 21, 15, 22, 14 
     $      ,  4, 16,  5, 19, 25, 17,  7, 18,  6  /

      integer nmlast,nxlast
      save    nmlast,nxlast
      data    nmlast,nxlast / 0,0 /

      real gh_edge(3,3,3),gh_vtx(3,3,3),zgh(3)
      save zgh
      data zgh / -1., 0., 1. /
      integer gh_type, i, j, ldw

      if (nx.ne.nxlast) then ! define interp. op
         nxlast = nx
         call zwgll (zl,interp,nx)
         call igllm (interp,interpt,zgh,zl,3,nx,3,nx)
      endif

      do i=1,3**ndim     ! currently support only 3x3x3 in moab
         j = moabmap(i)
         xt(i,1,1) = x27(j)
      enddo

      gh_type = 1 ! vertex only extension
      gh_type = 2 ! edge extension
      call gh_face_extend(xt,zgh,3,gh_type,gh_edge,gh_vtx)

c     Interpolate from 3x3x3 to (nx1 x ny1 x nz1) SEM mesh
      ldw = 3*lx1*ly1*lz1
      call map_to_crs(x,nx1,xt,3,if3d,wk,ldw)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_copy_vars()
      implicit none
#include "NEKMOAB"
      include 'GEOM'
      include 'SOLN'
      integer i, j, ierr, ntot, tmpcount, count, atend
      real tag_ptr(1)

      ntot = nx1*ny1*nz1
      do i = 1, numflu+numoth

         call iMesh_resetEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), ierr)
         IMESH_ASSERT

         atend = 0
         count = 0
         do while (atend .eq. 0)
c use the same iterator for all variables, since the elems are the same
            call nekMOAB_set_tag(ieiter(i), xm1Tag, ntot, tmpcount, xm1)
            call nekMOAB_set_tag(ieiter(i), ym1Tag, ntot, tmpcount, ym1)
            call nekMOAB_set_tag(ieiter(i), zm1Tag, ntot, tmpcount, zm1)

            call nekMOAB_set_tag(ieiter(i), vxTag, ntot, tmpcount, vx)
            call nekMOAB_set_tag(ieiter(i), vyTag, ntot, tmpcount, vy)
            call nekMOAB_set_tag(ieiter(i), vzTag, ntot, tmpcount, vz)

            call nekMOAB_set_tag(ieiter(i), tTag, ntot, tmpcount, t)

            if (nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) then
               call nekMOAB_set_tag(ieiter(i), pTag, ntot, tmpcount, pr)
            endif

c     step the iterator
            call iMesh_stepEntArrIter(%VAL(imeshh), %VAL(ieiter(i)), 
     $           %VAL(tmpcount), atend, ierr)
            IMESH_ASSERT

            count = count + tmpcount
         enddo

c double-check the total number of elements in this set
         if (count .ne. iecount(i)) then
            call exitti('Wrong no of elems iterating over matset ', 
     $           matids(i))
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_set_tag(iter, tagh, size, count, vals)
      implicit none

#include "NEKMOAB"      
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      integer ierr, i, ivals, size, count
      real vals(*), tag_vals
      pointer(tag_ptr, tag_vals(1))

      call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), 
     $     %VAL(iter), tag_ptr, count, ierr)
      IMESH_ASSERT

c set the tag vals
      ivals = size * count
      do i = 1, ivals
         tag_vals(i) = vals(i)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_compute_diagnostics()
      implicit none

#include "NEKMOAB"
      include 'GEOM'

      integer i, j, k, l

      integer e, intv
      real avg(3)

      intv = nelt / 10
      e = 1
      do while (e .lt. nelt)
         avg(1) = 0.0
         avg(2) = 0.0
         avg(3) = 0.0
         do k=1,nz1
            do j=1,ny1
               do i=1,nx1
                  avg(1) = avg(1) + xm1(i,j,k,e)
                  avg(2) = avg(2) + ym1(i,j,k,e)
                  if (if3d)
     $                 avg(3) = avg(3) + zm1(i,j,k,e)
               enddo
            enddo
         enddo
         avg(1) = avg(1) / (nx1*ny1*nz1)
         avg(2) = avg(2) / (nx1*ny1*nz1)
         avg(3) = avg(3) / (nx1*ny1*nz1)
         print *, "Average for e is ", e, avg(1), avg(2), avg(3)

         e = e + intv
      enddo

      return
      end
