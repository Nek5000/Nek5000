#ifdef PTRSIZE8
#define IMESH_HANDLE integer*8
#else
#define IMESH_HANDLE integer*4
#endif

#define IMESH_NULL 0
#define IMESH_ASSERT(ierr, imesh) if (ierr .ne. 0) call imesh_err(ierr, imesh, 0,__LINE__)
#define IMESH_NULLSTRIP(s) s(:index(s, char(0))-1)

#define MYLOC LOC
#define MYVAL %val

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
      include 'SIZE'
      include 'TOTAL'

      common /mbc/ moabbc(6,lelt) 
      integer moabbc

      call nekMOAB_load                   ! read mesh using MOAB
      call nekMOAB_proc_map               ! create gllnid mapping

      nelgt = nelgv
      nelt  = nelv

      call chk_nel

      call mapelpr                        ! create gllel mapping 
      call moab_geometry(xm1,ym1,zm1)     ! fill xm1,ym1,zm1
      call xml2xc                         ! fill xc,yc,zc

      call ifill(moabbc, 0, 6*lelt)
      call nekMOAB_BC(moabbc)             ! read MOAB BCs 

c      call nekMOAB_loadMaterialSets

#if 0
      do i=1,6
         call outface(xm1,1,i,nx1,ny1,nz1,'xface ')
         call outface(ym1,1,i,nx1,ny1,nz1,'yface ')
         call outface(zm1,1,i,nx1,ny1,nz1,'zface ')
      enddo
      call exitt
#endif

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_load
c
c     Load "filename" into imesh/moab, store imesh handle 
c     in /nekmoab/ common block
c
      implicit none
#include "iMesh_f.h"
#include "NEKMOAB"

      character*(*) filename
      parameter(filename='input.h5m')

      character*(*) loadOpt
c     ! only root processor will read and broadcast
c      parameter(loadOpt=";PARALLEL=READ_PART;PARTITION=PARALLEL_PARTIT
c     $ION;PARTITION_DISTRIBUTE;")
      ! all processor will read the input file and delete what is not needed
      parameter(loadOpt="moab:PARALLEL=READ_DELETE  moab:PARTITION=PARAL
     $LEL_PARTITION moab:PARTITION_DISTRIBUTE")

      character*(*) loadMode
      parameter(loadMode="PARALLEL")

      integer ierr
      IMESH_HANDLE rootst

      !Initialize imesh and load file
      imesh = IMESH_NULL
      call iMesh_newMesh(loadMode, imesh, ierr)
      IMESH_ASSERT(ierr, imesh)

      call iMesh_getRootSet(MYVAL(imesh), rootst, ierr)

      call iMesh_load(MYVAL(imesh), MYVAL(rootst), filename, 
     $     loadOpt,
     $     ierr)
      IMESH_ASSERT(ierr, imesh)

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
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'SCRCT'
      include 'ZPER'

      include 'mpif.h'

#include "iMesh_f.h"
#include "NEKMOAB"

      integer globalId, lastGlobalId

      call iMesh_getTagHandle(MYVAL(imesh),
     $     "GLOBAL_ID", !/*in*/ const char* tag_name,
     $     globalIdTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT(ierr, imesh)

      
      ifooAlloc = 0
      hexHandlesPointer = IMESH_NULL
      hexSize = 0
      call iMesh_getEntities(MYVAL(imesh), 
     $     MYVAL(IMESH_NULL), 
     $     MYVAL(iBase_REGION), 
     $     MYVAL(iMesh_ALL_TOPOLOGIES), 
     $     hexHandlesPointer, ifooAlloc, hexSize, ierr)
      IMESH_ASSERT(ierr, imesh)

c      nelgv = hexSize
      nelv = hexSize

      if(nelv .le. 0 .or. nelv .gt. lelv) then
         print *, 'ABORT: nelv is invalid in nekmoab_proc_map'
         print *, 'nelv, lelv = ', nelv, lelv
         call exitt
      endif

      nelgv = iglsum(nelv,1)
      if (NELGV .gt. LELG) then
         print *, 'ABORT; increase lelg ',nelgv,lelg
         call exitt
      endif

      call izero(GLLNID, NELGV)
      lastGlobalId = -1
      do i=1, hexSize
         call iMesh_getIntData(MYVAL(imesh), !iMesh_Instance instance,
     $        MYVAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        MYVAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        globalId,
     $        ierr) 
         IMESH_ASSERT(ierr, imesh)
c         print *, nid, 'local hex ', i, 'has global_id ', globalId

         !consistency check
         if(globalId .lt. lastGlobalId) then
            print *, "bug! this won't work"
            call exitt
         endif
         lastGlobalId = globalId

         if(globalId .gt. nelgv) then 
           write(6,*) 'ABORT: invalid  globalId! ', globalId
           call exitt
         endif
         GLLNID(globalId) = nid
      end do

      call igop(GLLNID, GLLEL, '+  ', NELGV)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_matSet2(matl, lelt, setHandle, setId)
c
c     Internal function, don't call directly
c
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer lelt
      integer matl(1)

      IMESH_HANDLE setHandle
      integer setId !coming from cubit

      IMESH_HANDLE hexes(*)
      pointer (hexesPointer, hexes)
      integer hexesSize, hexesAlloc
      integer ierr, i, elno


      !what hexes are in this set?
      hexesAlloc = 0
      hexesPointer = IMESH_NULL
      call iMesh_getEntitiesRec(MYVAL(imesh), 
     $     MYVAL(setHandle), 
     $     MYVAL(iBase_REGION), 
     $     MYVAL(iMesh_HEXAHEDRON),
     $     MYVAL(1),
     $     hexesPointer, hexesAlloc, hexesSize, ierr)
      IMESH_ASSERT(ierr, imesh)

      do i=1, hexesSize
         call nekMOAB_getElNo(hexes(i), elno)
         matl(elno) = setId
      enddo

      !ASDD call iMesh_free(hexesPointer)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_loadConn(vertex, nelgt, nelt, ncrnr)
c
c     Fill the vertex array with connectivity data from imesh
c
c     vertex(ncrnr, nelgt): int array, global id of a vertex
c     nelgt: int, global number of elements
c     ncrnr: int, number of corner vertices per element (should be 8)
c
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer nelgt, ncrnr
      integer vertex(ncrnr, 1)

      integer nelt
      integer globalId, ierr,i,k,m,n,npass,lex,ivmin,ivmax,iglmin,iglmax

      integer vtxId(TWENTYSEVEN)
      pointer (vtxIdPointer, vtxId)
      integer vtxIdSize, vtxIdAlloc
      integer vtxIdData(TWENTYSEVEN)

      IMESH_HANDLE vtxHandles(TWENTYSEVEN)
      pointer (vtxHandlesPointer, vtxHandles)
      integer vtxHandlesSize, vtxHandlesAlloc
      IMESH_HANDLE vtxHandlesData(TWENTYSEVEN)

      integer l2c(27),j
      common /cccc/ l2c

      vtxIdPointer = MYLOC(vtxIdData)
      vtxIdAlloc = TWENTYSEVEN

      vtxHandlesPointer = MYLOC(vtxHandlesData)
      vtxHandlesAlloc = TWENTYSEVEN

      call get_l2c_8(l2c)  ! get cubit-to-lexicographical ordering

c      write(6,*) hexSize,nelt,nelgt,ncrnr
      call izero(vertex, nelt*ncrnr)

      do i=1, hexSize
        call iMesh_getIntData(MYVAL(imesh), !iMesh_Instance instance,
     $        MYVAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        MYVAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        globalId,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);
         IMESH_ASSERT(ierr, imesh)
         
         call iMesh_getEntAdj(MYVAL(imesh), 
     $         MYVAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $         MYVAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $         vtxHandlesPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $         vtxHandlesAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $         vtxHandlesSize , !/*out*/ int* adj_entity_handles_size,
     $         ierr)
         IMESH_ASSERT(ierr, imesh)

         call iMesh_getIntArrData(MYVAL(imesh), !iMesh_Instance instance,
     $         MYVAL(vtxHandlesPointer), !/*in*/ const iBase_EntityHandle* entity_handles
     $         MYVAL(vtxHandlesSize), !/*in*/ const int entity_handles_size,
     $         MYVAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $         vtxIdPointer, !/*inout*/ int** tag_values,
     $         vtxIdAlloc, !/*inout*/ int* tag_values_allocated,
     $         vtxIdSize, !/*out*/ int* tag_values_size,
     $         ierr)

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

#include "iMesh_f.h"
#include "NEKMOAB"
      real x27(27,1), y27(27,1), z27(27,1)

      IMESH_HANDLE vtxHandles(TWENTYSEVEN)
      pointer (vtxHandlesPointer, vtxHandles)
      integer vtxHandlesSize, vtxHandlesAlloc
      IMESH_HANDLE vtxHandlesData(TWENTYSEVEN)

      integer i, j, ierr

      vtxHandlesPointer = MYLOC(vtxHandlesData)
      vtxHandlesAlloc = TWENTYSEVEN

      do i=1, hexSize         
         call iMesh_getEntAdj(MYVAL(imesh), 
     $                        MYVAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        MYVAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $                        vtxHandlesPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                        vtxHandlesAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                        vtxHandlesSize , !/*out*/ int* adj_entity_handles_size,
     $                        ierr)
         IMESH_ASSERT(ierr, imesh)
         
         if(vtxHandlesSize .ne. TWENTYSEVEN) then
            print *, 'Bad mesh! Element with', vtxHandlesSize, ' nodes'
            call exitt()
         endif

         do j=1, TWENTYSEVEN
            call iMesh_getVtxCoord(MYVAL(imesh),
     $                             MYVAL(vtxHandles(j)), !/*in*/ const iBase_EntityHandle vertex_handle,
     $                             x27(j,i), y27(j,i), z27(j,i), !/*out*/ double *x, /*out*/ double *y, /*out*/ double *z,
     $                             ierr)
            IMESH_ASSERT(ierr, imesh)
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

#include "iMesh_f.h"
#include "NEKMOAB"

      integer moabbc(6,1)

      IMESH_HANDLE entSetHandles(*)
      pointer (entSetHandlesPointer, entSetHandles)
      integer entSetAllocated, entSetSize

      IMESH_HANDLE neuSetTag
      integer ierr, i, tagIntData


      !Sidesets in cubit come in as entity sets with the NEUMANN_SET -- see sample file
      call iMesh_getTagHandle(MYVAL(imesh),
     $     "NEUMANN_SET", !/*in*/ const char* tag_name,
     $     neuSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT(ierr, imesh)


      entSetHandlesPointer = IMESH_NULL
      entSetAllocated      = 0
      call iMesh_getEntSets(MYVAL(imesh),
     $     MYVAL(IMESH_NULL),     !/*in*/ const iBase_EntitySetHandle entity_set_handle,
     $     MYVAL(1),              !/*in*/ const int num_hops,
     $     entSetHandlesPointer, !/*out*/ iBase_EntitySetHandle** contained_set_handles,
     $     entSetAllocated,      !/*out*/ int* contained_set_handles_allocated,
     $     entSetSize,           !/*out*/ int* contained_set_handles_size,
     $     ierr)                 !/*out*/ int *err);
      IMESH_ASSERT(ierr, imesh)

      !print *, 'H3', entSetSize

      do i=1, entSetSize
         call iMesh_getIntData(MYVAL(imesh), !iMesh_Instance instance,
     $        MYVAL(entSetHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        MYVAL(neuSetTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        tagIntData,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);

         if (ierr .eq. 0) then !tag was defined
            !print *, 'H32', tagIntData
            call nekMOAB_intBC(moabbc, entSetHandles(i), tagIntData)
         endif
      enddo

      !ASDD call iMesh_free(entSetHandlesPointer)

      return
      end 
c-----------------------------------------------------------------------
      subroutine nekMOAB_intBC(bcdata, setHandle, setId)
c
c     Internal function, don't call directly
c
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer bcdata(6,1)

      IMESH_HANDLE setHandle
      integer setId !coming from cubit

      IMESH_HANDLE faces(*)
      pointer (facesPointer, faces)
      integer facesSize, facesAlloc

      IMESH_HANDLE ahex(*)
      pointer (ahexPointer, ahex)
      integer ahexSize, ahexAlloc

      IMESH_HANDLE fvtx(*)
      pointer (fvtxPointer, fvtx)
      integer fvtxSize, fvtxAlloc

      IMESH_HANDLE hvtx(*)
      pointer (hvtxPointer, hvtx)
      integer hvtxSize, hvtxAlloc

      integer ierr, i, j, elno

      integer side_no, side_offset, side_sense
      integer hexMBCNType

      call iMesh_MBCNType(MYVAL(iMesh_HEXAHEDRON), hexMBCNType)

      !what faces are in this set?
      facesAlloc = 0
      facesPointer = IMESH_NULL
      call iMesh_getEntitiesRec(MYVAL(imesh), 
     $     MYVAL(setHandle), 
     $     MYVAL(iBase_FACE), 
     $     MYVAL(iMesh_QUADRILATERAL),
     $     MYVAL(1),
     $     facesPointer, facesAlloc, facesSize, ierr)
      IMESH_ASSERT(ierr, imesh)

      do i=1, facesSize
         !get vertices defining the face
         fvtxAlloc = 0
         fvtxPointer = IMESH_NULL
         call iMesh_getEntAdj(MYVAL(imesh),
     $                        MYVAL(faces(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        MYVAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $                        fvtxPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                        fvtxAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                        fvtxSize , !/*out*/ int* adj_entity_handles_size,
     $                        ierr)
         IMESH_ASSERT(ierr, imesh)         

         !get hexes adjacent to the face (1 or 2) 
         ! -- hopefully only returns hexes on local proc, but untested
         ahexAlloc = 0
         ahexPointer = IMESH_NULL                  
         call iMesh_getEntAdj(MYVAL(imesh), 
     $                        MYVAL(faces(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        MYVAL(iBase_REGION), !/*in*/ const int entity_type_requested,
     $                        ahexPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                        ahexAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                        ahexSize , !/*out*/ int* adj_entity_handles_size,
     $                        ierr)
         IMESH_ASSERT(ierr, imesh)

         do j=1, ahexSize                                
            !get verts adjacent to the hex
            hvtxAlloc = 0
            hvtxPointer = IMESH_NULL                  
            call iMesh_getEntAdj(MYVAL(imesh), 
     $                           MYVAL(ahex(j)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                           MYVAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $                           hvtxPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                           hvtxAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                           hvtxSize , !/*out*/ int* adj_entity_handles_size,
     $                           ierr)
            IMESH_ASSERT(ierr, imesh)

            !Get the side number
            call MBCN_SideNumberUlong(MYVAL(hvtxPointer),
     $                                MYVAL(hexMBCNType), 
     $                                MYVAL(fvtxPointer), 
     $                                MYVAL(4), !not fvtxSize, because that's 5 and crashes stuff 
                                               ! -- linear nodes are stored first and that's all that's used
     $                                MYVAL(2), !Dimension of faces
     $                    side_no, side_sense, side_offset) !output
           if ((side_no .lt. 0) .or. (side_no .gt. 5)) ierr = 1
           IMESH_ASSERT(ierr, imesh)

           side_no = side_no + 1 !moab side number is zero based

           !call nekMOAB_getGlobElNo(ahex(j), elno)
           !print *, 'BLAH', elno, side_no, setId

           call nekMOAB_getElNo(ahex(j), elno)

           !print *, 'setting', elno, side_no, setId

           bcdata(side_no, elno) = setId

           !ASDD call iMesh_free(hvtxPointer)
         enddo


         !ASDD call iMesh_free(ahexPointer)
         !ASDD call iMesh_free(fvtxPointer)

      enddo

      !ASDD call iMesh_free(facesPointer)

      return
      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getElNo(handle, elno)
c
c     Given the imesh handle of an element (hex), 
c     return its local nek element number in elno
c     Cannot be used until the GLLEL array has been set (in map2.f)
c
      include 'SIZE'
      include 'PARALLEL'
#include "NEKMOAB"

      IMESH_HANDLE handle
      integer elno

      !first get global id
      call nekMOAB_getGlobElNo(handle, elno)
      !then convert to local id
      elno = GLLEL(elno)

      end
c-----------------------------------------------------------------------
      subroutine nekMOAB_getGlobElNo(handle, elno)

      include 'SIZE'
      include 'PARALLEL'
#include "NEKMOAB"

      IMESH_HANDLE handle
      integer elno
      
      integer ierr

      call iMesh_getIntData(MYVAL(imesh), !iMesh_Instance instance,
     $        MYVAL(handle), !/*in*/ const iBase_EntityHandle entity_handle,
     $        MYVAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        elno,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);
      IMESH_ASSERT(ierr, imesh)
      
      if (GLLNID(elno) .ne. NID) then
         print *, 'bug... fixme',elno, nid, gllnid(elno)
         call exitt
      endif
                          
      end
c-----------------------------------------------------------------------
      subroutine moab_to_nek_bc()  ! fill the nek cbc arrays
c
      include 'SIZE'
      include 'TOTAL'

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
      subroutine moab_geometry (xml,yml,zml)

      include 'SIZE'
      include 'TOTAL'

      parameter(lxyz=lx1*ly1*lz1)
      real      xml(lxyz,1)
     $        , yml(lxyz,1)
     $        , zml(lxyz,1)
      integer   e

      common /tcrmg/ x27(27,lelt), y27(27,lelt), z27(27,lelt)

      call nekMOAB_loadCoord(x27, y27, z27) !     Get coords from moab

c     call outmat(x27,27,8,'x27dat',nelt)
c     call outmat(y27,27,8,'y27dat',nelt)
c     call outmat(z27,27,8,'z27dat',nelt)

      nmoab = 3 !not used
      do e=1,nelt   !  Interpolate for each element
         call permute_and_map(xml(1,e),x27(1,e),nx1,nmoab,e)
         call permute_and_map(yml(1,e),y27(1,e),ny1,nmoab,e)
         if (if3d) call permute_and_map(zml(1,e),z27(1,e),nz1,nmoab,e)
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

      include 'SIZE'
      include 'GEOM'
      include 'INPUT'

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
      include 'SIZE'
      include 'TOTAL'

      real x(1),x27(0:1)

      common /ctmp0/ z3(3),zpt(lx1)
     $             , xt(3,3,3)
     $             , wk(3*lx1*ly1*lz1)
     $             , xw(3*lx1*ly1*lz1)

      real interp(lx1*3),interpt(lx1*3),zl(lx1*3)
      save interp,interpt

      integer moabmap(27),e
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
      integer gh_type

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

      integer l2c_save(8),l2c(8)
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
c#include "iMesh_f.h"
c#include "NEKMOAB"
c
c      !store material properties here
c      common /cmatl/ matl(lelg)
c      integer matl
c      
cc very similar to the bc-loading code
c
c      IMESH_HANDLE entSetHandles(*)
c      pointer (entSetHandlesPointer, entSetHandles)
c      integer entSetAllocated, entSetSize
c
c      IMESH_HANDLE matSetTag
c      integer ierr, i, tagIntData
c
cc      matl = -1
c
c      call iMesh_getTagHandle(MYVAL(imesh),
c     $     "MATERIAL_SET", !/*in*/ const char* tag_name,
c     $     matSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
c     $     ierr)
c      IMESH_ASSERT(ierr, imesh)
c
c      entSetHandlesPointer = IMESH_NULL
c      entSetAllocated      = 0
c      call iMesh_getEntSets(MYVAL(imesh),
c     $     MYVAL(IMESH_NULL),     !/*in*/ const iBase_EntitySetHandle entity_set_handle,
c     $     MYVAL(1),              !/*in*/ const int num_hops,
c     $     entSetHandlesPointer, !/*out*/ iBase_EntitySetHandle** contained_set_handles,
c     $     entSetAllocated,      !/*out*/ int* contained_set_handles_allocated,
c     $     entSetSize,           !/*out*/ int* contained_set_handles_size,
c     $     ierr)                 !/*out*/ int *err);
c      IMESH_ASSERT(ierr, imesh)
c
c      do i=1, entSetSize
c         call iMesh_getIntData(MYVAL(imesh), !iMesh_Instance instance,
c     $        MYVAL(entSetHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
c     $        MYVAL(matSetTag), !/*in*/ const iBase_TagHandle tag_handle,
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

