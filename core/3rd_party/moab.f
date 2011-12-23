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
#include "NEKMOAB"
      include 'PARALLEL'
      include 'GEOM'
      common /mbc/ moabbc(6,lelt) 
      integer moabbc

      call nekMOAB_load                   ! read mesh using MOAB
      call nekMOAB_proc_map               ! create gllnid mapping

      nelgt = nelgv
      nelt  = nelv

      call chk_nel

      call nekMOAB_assign_tag_storage             ! allocate MOAB tags to reflect Nek variables

      call mapelpr                        ! create gllel mapping 
      call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
      call xml2xc                         ! fill xc,yc,zc

      call ifill(moabbc, 0, 6*lelt)
      call nekMOAB_BC(moabbc)             ! read MOAB BCs 
      
      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_assign_tag_storage
c
      implicit none
#include "NEKMOAB"

      integer num_hexes, ierr
      iBase_EntityHandle iter
      IBASE_HANDLE_T tagh
      integer semdim(3), num_regions

c get the number of regions/hexes in the local mesh (need that to init the iterator)
      call iMesh_getNumOfType(%VAL(imeshh), %VAL(rootset), 
     1 %VAL(iBase_REGION), num_regions, ierr)
      IMESH_ASSERT

c get an iterator over regions
      call iMesh_initEntArrIter(%VAL(imeshh), %VAL(rootset), 
     1 %VAL(iBase_REGION), %VAL(iMesh_HEXAHEDRON),%VAL(num_hexes),
     1 %VAL(0), iter, ierr) 
      IMESH_ASSERT

      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "SEM_X", rpxm1)
      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "SEM_Y", rpym1)
      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "SEM_Z", rpzm1)

      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "VEL_X", rpvx)
      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "VEL_Y", rpvy)
      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "VEL_Z", rpvz)

      call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "TEMP", rpt)

      if (nx2.eq.nx1 .and. ny2.eq.ny1 .and. nz2.eq.nz1) then
         call nekMOAB_get_tag_storage(nx1*ny1*nz1, iter, "PRESS", rpp)
      endif

c get rid of the iterator, since we no longer need it
      call iMesh_endEntArrIter(%VAL(imeshh), %VAL(iter), ierr)
      IMESH_ASSERT

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
      subroutine nekMOAB_get_tag_storage(isz, iter, tag_name, rpdata) 
c
      implicit none
#include "NEKMOAB"

      iBase_EntityHandle iter
      character*(*) tag_name
      integer*8 rpdata
      integer*4 isz

      integer ierr, count
      IBASE_HANDLE_T tagh

c create the tag
      call iMesh_createTagWithOptions(%VAL(imeshh), tag_name, 
     1 "moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     1 %VAL(isz), %VAL(iBase_DOUBLE), tagh, ierr)
      IMESH_ASSERT

c iterate over tag memory
      call iMesh_tagIterate(%VAL(imeshh), %VAL(tagh), %VAL(iter), 
     1     rpdata, count, ierr)
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
      include 'INPUT'
      include 'mpif.h'
c two forms of load options, depending whether we're running serial or parallel      
      character*(*) parLoadOpt, serLoadOpt
      parameter(parLoadOpt=" moab:PARALLEL=READ_DELETE moab:PARTITION=PA
     $RALLEL_PARTITION moab:PARALLEL_RESOLVE_SHARED_ENTS moab:PARTITION_
     $DISTRIBUTE ")
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

      return
      end  
c-----------------------------------------------------------------------
      subroutine nekMOAB_proc_map()
c
c     Fills the following common-block data from a previously-loaded 
c     imesh (handle stored in /nekmoab/ common block
c     GLLNID : (int array)  Maps global element id to the processor 
c     that it's on.
c     NELGV  : (int) global number of elements
c     
c     Gotcha: Out of GLLNID, map2.f "computes" GLLEL (map from global 
c     element number to local element number). To do this, it assumes 
c     that if A has a lower global element number than B, it also
c     has a lower local element number. 
c     Who knows if this will be true with data coming from imesh.
c
      implicit none
#include "NEKMOAB"
      include 'INPUT'
      include 'PARALLEL'
      include 'SCRCT'
      include 'ZPER'

      include 'mpif.h'


      integer globalId, lastGlobalId
      integer pstatus(*), hgid(*), etype(*)
      pointer(rppstatus, pstatus), (rphgid, hgid), (rpetype, etype)
      integer pstatusSize, hgidSize, etypeSize, etypeAlloc, i, ierr, ip, 
     $     ih, nonloc_hexes
      integer iglsum

      call iMesh_getTagHandle(%VAL(imeshh),
     $     "GLOBAL_ID",       !/*in*/ const char* tag_name,
     $     globalIdTag,       !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT

      
      rpHexes = IMESH_NULL
      hexesSize = 0
      call iMesh_getEntities(%VAL(imeshh), 
     $     %VAL(rootset), 
     $     %VAL(iBase_REGION), 
     $     %VAL(iMesh_ALL_TOPOLOGIES), 
     $     rpHexes, hexesSize, hexesSize, ierr)
      IMESH_ASSERT

c     get entity types, just as a check
      rpetype = IMESH_NULL
      etypeAlloc = 0
      call iMesh_getEntArrTopo(%VAL(imeshh), 
     $     %VAL(rpHexes), %VAL(hexesSize), rpetype, 
     $     etypeAlloc, etypeSize, ierr)
      IMESH_ASSERT
      do i = 1, etypeSize
         if (etype(i) .ne. iMesh_HEXAHEDRON) then
            print *, "Not all hex elements!"
            call exitt
         endif
      enddo

c     count the number of local entities; reuse etype array
      do i = 1, etypeSize
         etype(i) = 1
      enddo
      pstatusSize = 0
      rppstatus = IMESH_NULL
      nonloc_hexes = 0
      hgidSize = 0
      rphgid = IMESH_NULL
      call iMesh_getIntArrData(%VAL(imeshh), 
     $     %VAL(rpHexes), %VAL(hexesSize),
     $     %VAL(globalIdTag), rphgid, hgidSize, hgidSize, ierr) 
      IMESH_ASSERT

      do ip = 1, partsSize
         call iMeshP_getEntStatusArr(%VAL(imeshh), %VAL(hPartn), 
     $        %VAL(hParts(ip)), %VAL(rpHexes), %VAL(hexesSize),
     $        rppstatus, pstatusSize, pstatusSize, ierr)
         IMESH_ASSERT
         do ih = 1, hexesSize
            if (pstatus(ih) .ne. iMeshP_INTERNAL) then
               nonloc_hexes = nonloc_hexes+1
               etype(ih) = 0
            endif
         enddo
      enddo
      if (nonloc_hexes .ne. 0) then
         print *, 'Number of local, total hexes = ', 
     $        hexesSize-nonloc_hexes, hexesSize  
      endif

c     nelgv = hexSize
      nelv = hexesSize - nonloc_hexes

      if(nelv .le. 0 .or. nelv .gt. lelv) then
         print *, 'ABORT: nelv is invalid in nekmoab_proc_map'
         print *, 'nelv, lelv = ', nelv, lelv
         call exitt
      endif

c     check global total number of hexes
      nelgv = iglsum(nelv,1)
      if (NELGV .gt. LELG) then
         print *, 'ABORT; increase lelg ',nelgv,lelg
         call exitt
      endif

c     check global id space for monotonicity
c     get all global ids at once, reuse rppstatus
      call iMesh_getIntArrData(%VAL(imeshh), 
     $     %VAL(rpHexes), %VAL(hexesSize),
     $     %VAL(globalIdTag), 
     $     rppstatus, pstatusSize, pstatusSize, ierr)
      IMESH_ASSERT

      call izero(GLLNID, NELGV)
      lastGlobalId = -1
      do ih=1, hexesSize
c         print *, nid, 'local hex ', ih, 'has global_id ', pstatus(ih)

!     consistency check
         if(pstatus(ih) .lt. lastGlobalId .and. etype(ih) .eq. 1) then
            print *, "Non-monotonic global id space!"
            call exitt
         endif
         lastGlobalId = pstatus(ih)

         if(pstatus(ih) .gt. nelgv) then 
            write(6,*) 'ABORT: invalid  globalId! ', globalId
            call exitt
         endif
         GLLNID(pstatus(ih)) = nid
      end do

      call igop(GLLNID, GLLEL, '+  ', NELGV)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_matSet2(matl, setHandle, setId)
c
c     Internal function, don't call directly
c
      implicit none
#include "NEKMOAB"

      integer matl(1)

      IBASE_HANDLE_T setHandle
      integer setId !coming from cubit

      IBASE_HANDLE_T loc_hexes(*)
      pointer (rploc_hexes, loc_hexes)
      integer loc_hexesSize
      integer ierr, i, elno


      !what hexes are in this set?
      loc_hexesSize = 0
      rploc_hexes = IMESH_NULL
      call iMesh_getEntitiesRec(%VAL(imeshh), 
     $     %VAL(setHandle), 
     $     %VAL(iBase_REGION), 
     $     %VAL(iMesh_HEXAHEDRON),
     $     %VAL(1),
     $     rploc_hexes, loc_hexesSize, loc_hexesSize, ierr)
      IMESH_ASSERT

      do i=1, loc_hexesSize
         call nekMOAB_getElNo(loc_hexes(i), elno)
         matl(elno) = setId
      enddo

      call free(rploc_hexes)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_loadConn(vertex, nelgt, ncrnr)
c
c     Fill the vertex array with connectivity data from imesh
c
c     vertex(ncrnr, nelgt): int array, global id of a vertex
c     nelgt: int, global number of elements
c     ncrnr: int, number of corner vertices per element (should be 8)
c
      implicit none
#include "NEKMOAB"

      integer vertex(ncrnr, 1)

      integer globalId, ierr, i, k

      integer vtxId(TWENTYSEVEN), vtxIdSize, itmp1
      pointer (rpvtxId, itmp1)

      IBASE_HANDLE_T hvtx(TWENTYSEVEN)
      pointer (rpvtx, htmp)
      integer vtxSize
      IBASE_HANDLE_T htmp

      integer l2c(27),j
      common /cccc/ l2c

      rpvtxId = MYLOC(vtxId)
      vtxIdSize = TWENTYSEVEN

      rpvtx = MYLOC(hvtx)
      vtxSize = TWENTYSEVEN

      call get_l2c_8(l2c)  ! get cubit-to-lexicographical ordering

c      write(6,*) hexSize,nelt,nelgt,ncrnr
      call izero(vertex, nelt*ncrnr)

      do i=1, hexesSize
        call iMesh_getIntData(%VAL(imeshh), !iMesh_Instance instance,
     $        %VAL(hHexes(i)), %VAL(globalIdTag), globalId, ierr)
         IMESH_ASSERT
         
         call iMesh_getEntAdj(%VAL(imeshh), 
     $        %VAL(hHexes(i)), %VAL(iBase_VERTEX), 
     $        rpvtx, vtxSize, vtxSize, ierr)
         IMESH_ASSERT

         call iMesh_getIntArrData(%VAL(imeshh), 
     $        %VAL(rpvtx), %VAL(vtxSize), %VAL(globalIdTag), 
     $        rpvtxId, vtxIdSize, vtxIdSize, ierr)
         IMESH_ASSERT

         do k=1, ncrnr
            j = l2c(k)
c            vertex(k, globalId) = vtxId(j)
            vertex(k, i) = vtxId(j)
c            write(6,*)  i,j,k,globalid,vertex(k,globalid),' id'
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
      real x27(27,1), y27(27,1), z27(27,1)

      IBASE_HANDLE_T hvtx(TWENTYSEVEN)
      pointer (rpvtx, htmp)
      integer vtxSize
      IBASE_HANDLE_T htmp

      integer i, j, ierr

      rpvtx = MYLOC(hvtx)
      vtxSize = TWENTYSEVEN

      do i=1, hexesSize         
         call iMesh_getEntAdj(%VAL(imeshh), 
     $        %VAL(hHexes(i)), %VAL(iBase_VERTEX), 
     $        rpvtx, vtxSize, vtxSize, ierr)
         IMESH_ASSERT
         
         if(vtxSize .ne. TWENTYSEVEN) then
            print *, 'Bad mesh! Element with', vtxSize, ' nodes'
            call exitt()
         endif

         do j=1, TWENTYSEVEN
            call iMesh_getVtxCoord(%VAL(imeshh),
     $           %VAL(hvtx(j)), x27(j,i), y27(j,i), z27(j,i), ierr)
            IMESH_ASSERT
         enddo

c         print *, i, 'coords: ' 
c         print *, 'X',  x27
c         print *, 'Y',  y27
c         print *, 'Z',  z27
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

      integer moabbc(6,1)

      IBASE_HANDLE_T hentSet(*)
      pointer (rpentSet, hentSet)
      integer entSetSize

      IBASE_HANDLE_T neuSetTag
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
c     Internal function, don't call directly
c

      implicit none
#include "NEKMOAB"

      integer bcdata(6,1)

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

           !print *, 'setting', elno, side_no, setId

           bcdata(side_no, elno) = setId

           call free(rphvtx)
         enddo


         call free(rpahex)
         call free(rpfvtx)

      enddo

      call free(rpfaces)

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
      subroutine moab_to_nek_bc()  ! fill the nek cbc arrays
c
#include "NEKMOAB"      
      include 'INPUT'
      common /mbc/ moabbc(6,lelt)
      integer moabbc

      integer e,f
      character*3 cbi(ldimt1)

      lcbc=18*lelt*(ldimt1 + 1)
      call blank(cbc,lcbc)

      nface = 2*ndim
      do e=1,nelt
      do f=1,nface
         call usr_moab2nek(cbi, moabbc(f,e))
         do ifld=1,nfield
            cbc(f,e,ifld) = cbi(ifld)
c            write(6,*) cbc(f,e,ifld),e
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine moab_geometry (xmlo,ymlo,zmlo)

      implicit none
#include "NEKMOAB"      
      include 'INPUT'

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
      include 'INPUT'
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
      include 'INPUT'
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
      subroutine get_l2c_8(l2c)
      implicit none

      integer l2c_save(8),l2c(8), i
      save    l2c_save
      data    l2c_save / 0, 1, 3, 2, 4, 5, 7, 6 /

      call icopy(l2c,l2c_save,8)
      do i=1,8
         l2c(i) = l2c(i) + 1
      enddo

      return
      end
c-----------------------------------------------------------------------
c      subroutine nekMOAB_loadMaterialSets
c      include 'SIZE'
c
c#include "NEKMOAB"
c
c      !store material properties here
c      common /cmatl/ matl(lelg)
c      integer matl
c      
cc very similar to the bc-loading code
c
c      IBASE_HANDLE_T entSetHandles(*)
c      pointer (entSetHandlesPointer, entSetHandles)
c      integer entSetAllocated, entSetSize
c
c      IBASE_HANDLE_T matSetTag
c      integer ierr, i, tagIntData
c
cc      matl = -1
c
c      call iMesh_getTagHandle(%VAL(imesh),
c     $     "MATERIAL_SET", !/*in*/ const char* tag_name,
c     $     matSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
c     $     ierr)
c      IMESH_ASSERT(ierr, imesh)
c
c      entSetHandlesPointer = IMESH_NULL
c      entSetAllocated      = 0
c      call iMesh_getEntSets(%VAL(imesh),
c     $     %VAL(IMESH_NULL),     !/*in*/ const iBase_EntitySetHandle entity_set_handle,
c     $     %VAL(1),              !/*in*/ const int num_hops,
c     $     entSetHandlesPointer, !/*out*/ iBase_EntitySetHandle** contained_set_handles,
c     $     entSetAllocated,      !/*out*/ int* contained_set_handles_allocated,
c     $     entSetSize,           !/*out*/ int* contained_set_handles_size,
c     $     ierr)                 !/*out*/ int *err);
c      IMESH_ASSERT(ierr, imesh)
c
c      do i=1, entSetSize
c         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
c     $        %VAL(entSetHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
c     $        %VAL(matSetTag), !/*in*/ const iBase_TagHandle tag_handle,
c     $        tagIntData,       !/*out*/ int *out_data,
c     $        ierr)             !/*out*/ int *err);
c
c         if (ierr .eq. 0) then !tag was defined
c            call nekMOAB_matSet2(matl, lelt, 
c     $           entSetHandles(i), tagIntData)
c         endif
c      enddo
c
c      call free(entSetHandlesPointer)
c
cc      print *, 'matl'
cc      print *, matl
c
c      return
c      end 
c-----------------------------------------------------------------------
