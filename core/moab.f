c-----------------------------------------------------------------------
!Some convenience macros. Maybe they can end up in the iBase_f.h later on


#define IMESH_HANDLE integer*8
#define IMESH_NULL 0
#define IMESH_ASSERT(ierr, imesh) if (ierr .ne. 0) call imesh_err(ierr, imesh, __FILE__, __LINE__)
#define IMESH_NULLSTRIP(s) s(:index(s, char(0))-1)

c
c    Data shared between functions in this file is stored in 
c   /nekmoab/ common block, in file NEKMOAB
c
c    We are targeting hex27 meshes coming from cubit. 
c    Sample case is in data/framework/testmeshes/2bricks/8hex.jou


c    fun stuff: 
c    "imesh" is already used in TSTEP
c    can't use implicit none because [A-Z]* includes 
c    depend on it -- watch out for typos!
c
c    remember to free things... 
c           currently using the fortran free to free memory malloced 
c           in c++. who knows what that does.

c    used igop instead of direct call to mpi_allreduce + copy, 
c    but no equivalent for allgather. Should use nek comm 
c    instead of worldcomm




      subroutine nekMOAB_load(filename)
c     Load "filename" into imesh/moab, store imesh handle 
c     in /nekmoab/ common block
      implicit none
#include "iMesh_f.h"
#include "NEKMOAB"

      character(*), intent(in)       :: filename
      
      character(*) loadOpt

#if 0
      parameter(loadOpt="PARALLEL=BCAST_DELETE;
     $PARTITION=GEOM_DIMENSION;PARTITION_VAL=3;
     $PARTITION_DISTRIBUTE;
     $PARALLEL_RESOLVE_SHARED_ENTS")
#endif 

#if 0
      parameter(loadOpt="PARALLEL=BCAST_DELETE;
     $PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;")
#endif      

#if 0
      parameter(loadOpt="PARALLEL=BCAST_DELETE;
     $PARTITION=GEOM_DIMENSION;PARTITION_DISTRIBUTE;
     $PARALLEL_RESOLVE_SHARED_ENTS;
     $PARTITION_VAL=3")
#endif      

#if 1
      parameter(loadOpt="")
#endif      

      integer       :: ierr

      !Initialize imesh and load file
      imesh = IMESH_NULL
      call iMesh_newMesh("PARALLEL", imesh, ierr)
      IMESH_ASSERT(ierr, imesh)

      call iMesh_load(%VAL(imesh), %VAL(IMESH_NULL), filename, 
     $     loadOpt,
     $     ierr)
      IMESH_ASSERT(ierr, imesh)

      end subroutine 

      subroutine nekMOAB_proc_map()
!     Fills the following common-block data from a previously-loaded 
!     imesh (handle stored in /nekmoab/ common block
!     GLLNID : (int array)  Maps global element id to the processor 
!     that it's on.
!     NELGV  : (int) global number of elements

!     Gotcha: Out of GLLNID, map2.f "computes" GLLEL (map from global 
!     element number to local element number). To do this, it assumes 
!     that if A has a lower global element number than B, it also
!     has a lower local element number. 
!     Who knows if this will be true with data coming from imesh.

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'SOLN'
      include 'SCRCT'
      include 'ZPER'

      include 'mpif.h'

#include "iMesh_f.h"
#include "NEKMOAB"

      integer ierr
      integer fooAlloc, i

      integer globalId, lastGlobalId

      call iMesh_getTagHandle(%VAL(imesh),
     $     "GLOBAL_ID", !/*in*/ const char* tag_name,
     $     globalIdTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT(ierr, imesh)
 
      fooAlloc = 0
      hexHandlesPointer = IMESH_NULL
      call iMesh_getEntities(%VAL(imesh), 
     $     %VAL(IMESH_NULL), 
     $     %VAL(iBase_REGION), 
     $     %VAL(iMesh_ALL_TOPOLOGIES), 
     $     hexHandlesPointer, fooAlloc, hexSize, ierr)
      IMESH_ASSERT(ierr, imesh)

      NELGV = hexSize
      write(6,*) nid,' hexsize:',hexsize,imesh

      if(nelgv .eq. 0) then
         write(6,*) 'nelgv is zero in nekmoab_proc_map'
         call exitt
      endif

      call igop(NELGV, LELGWORK, '+  ', 1)

      if (NELGV .gt. LELG) then
         print *, 'increase lelg ',nelgv,lelg
         call exitt
      endif

      call izero(GLLNID, NELGV)
      lastGlobalId = -1
      do i=1, hexSize
         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        globalId,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);
         IMESH_ASSERT(ierr, imesh)
         !print *, nid, 'local hex ', i, 'has global_id ', globalId

         !consistency check
         if(globalId .lt. lastGlobalId) then
            print *, "bug! this won't work"
            call exitt
         endif
         lastGlobalId = globalId

         GLLNID(globalId) = nid
      end do

      call igop(GLLNID, LELGWORK, '+  ', NELGV)


!      print *, 'NELGV:', NELGV
!      print *, 'GLLNID:'
!      print *,  GLLNID(:NELGV)


      end subroutine nekMOAB_proc_map
c-----------------------------------------------------------------------
      subroutine nekMOAB_loadMaterialSets
      include 'SIZE'

#include "iMesh_f.h"
#include "NEKMOAB"

      integer matl(lelt) !store material properties here
      common /cmatl/ matl
      
!very similar to the bc-loading code

      IMESH_HANDLE entSetHandles(*)
      pointer (entSetHandlesPointer, entSetHandles)
      integer :: entSetAllocated, entSetSize

      IMESH_HANDLE matSetTag
      integer:: ierr, i, tagIntData

      matl = -1

      call iMesh_getTagHandle(%VAL(imesh),
     $     "MATERIAL_SET", !/*in*/ const char* tag_name,
     $     matSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT(ierr, imesh)

      entSetHandlesPointer = IMESH_NULL
      entSetAllocated      = 0
      call iMesh_getEntSets(%VAL(imesh),
     $     %VAL(IMESH_NULL),     !/*in*/ const iBase_EntitySetHandle entity_set_handle,
     $     %VAL(1),              !/*in*/ const int num_hops,
     $     entSetHandlesPointer, !/*out*/ iBase_EntitySetHandle** contained_set_handles,
     $     entSetAllocated,      !/*out*/ int* contained_set_handles_allocated,
     $     entSetSize,           !/*out*/ int* contained_set_handles_size,
     $     ierr)                 !/*out*/ int *err);
      IMESH_ASSERT(ierr, imesh)

      do i=1, entSetSize
         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(entSetHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(matSetTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        tagIntData,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);

         if (ierr .eq. 0) then !tag was defined
            call nekMOAB_matSet2(matl, lelt, 
     $           entSetHandles(i), tagIntData)
         endif
      enddo

      call free(entSetHandlesPointer)

!      print *, 'matl'
!      print *, matl

      end subroutine nekMOAB_loadMaterialSets




      subroutine nekMOAB_matSet2(matl, lelt, setHandle, setId)
      !Internal function, don't call directly
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer lelt
      integer matl(lelt)

      IMESH_HANDLE :: setHandle
      integer :: setId !coming from cubit

      IMESH_HANDLE hexes(*)
      pointer (hexesPointer, hexes)
      integer hexesSize, hexesAlloc
      integer ierr, i, elno

      !what hexes are in this set?
      hexesAlloc = 0
      hexesPointer = IMESH_NULL
      call iMesh_getEntitiesRec(%VAL(imesh), 
     $     %VAL(setHandle), 
     $     %VAL(iBase_REGION), 
     $     %VAL(iMesh_HEXAHEDRON),
     $     %VAL(1),
     $     hexesPointer, hexesAlloc, hexesSize, ierr)
      IMESH_ASSERT(ierr, imesh)

      do i=1, hexesSize
         call nekMOAB_getElNo(hexes(i), elno)
         matl(elno) = setId
      enddo

      call free(hexesPointer)

      end subroutine nekMOAB_matSet2



c-----------------------------------------------------------------------
      subroutine nekMOAB_loadConn(vertex, nelgt, ncrnr)
      !Fill the vertex array with connectivity data from imesh

      !vertex(ncrnr, nelgt): int array, global id of a vertex
      !nelgt: int, global number of elements
      !ncrnr: int, number of corner vertices per element (should be 8)
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer nelgt, ncrnr
      integer vertex(ncrnr, nelgt)


      integer lwrk
      parameter (lwrk=1000)
      common /scruz/ wrk(lwrk)
      integer wrk

      integer globalId, ierr,i,k,m,n,npass,lex,ivmin,ivmax,iglmin,iglmax

      integer vtxId(TWENTYSEVEN)
      pointer (vtxIdPointer, vtxId)
      integer vtxIdSize, vtxIdAlloc
      integer vtxIdData(TWENTYSEVEN)

      IMESH_HANDLE vtxHandles(TWENTYSEVEN)
      pointer (vtxHandlesPointer, vtxHandles)
      integer vtxHandlesSize, vtxHandlesAlloc
      integer vtxHandlesData(TWENTYSEVEN)


      integer l2c(27),j
      common /cccc/ l2c
       


      vtxIdPointer = %LOC(vtxIdData)
      vtxIdAlloc = TWENTYSEVEN

      vtxHandlesPointer = %LOC(vtxHandlesData)
      vtxHandlesAlloc = TWENTYSEVEN


      call get_l2c_8(l2c)  ! get cubit-to-lexicographical ordering

      call izero(vertex, nelgt*ncrnr)

      do i=1, hexSize
         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        globalId,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);
         IMESH_ASSERT(ierr, imesh)
         
         call iMesh_getEntAdj(%VAL(imesh), 
     $            %VAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $            %VAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $            vtxHandlesPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $            vtxHandlesAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $            vtxHandlesSize , !/*out*/ int* adj_entity_handles_size,
     $            ierr)
         IMESH_ASSERT(ierr, imesh)

         call iMesh_getIntArrData(%VAL(imesh), !iMesh_Instance instance,
     $            %VAL(vtxHandlesPointer), !/*in*/ const iBase_EntityHandle* entity_handles
     $            %VAL(vtxHandlesSize), !/*in*/ const int entity_handles_size,
     $            %VAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $            vtxIdPointer, !/*inout*/ int** tag_values,
     $            vtxIdAlloc, !/*inout*/ int* tag_values_allocated,
     $            vtxIdSize, !/*out*/ int* tag_values_size,
     $            ierr)


         do k=1, ncrnr
            j = l2c(k)
            vertex(k, globalId) = vtxId(j)
c           write(6,*)  i,j,k,globalid,vertex(k,globalid),' id'
         enddo

      enddo

      n     = nelgt*ncrnr
      ivmin = iglmin(vertex,n)
      ivmax = iglmax(vertex,n)
      write(6,*) ivmin,ivmax,' ivminA ',nelgt

      npass = n/lwrk + 1
      k     = 1
      do j=1,npass
        m     = min(lwrk,n-k+1)
        call igop(vertex(k,1), wrk, '+  ', m)
        k = k+lwrk
        if (k.gt.n) goto 10
      enddo
   10 continue
      ivmin = iglmin(vertex,n)
      ivmax = iglmax(vertex,n)
      write(6,*) ivmin,ivmax,' ivmin ',nelgt
c     call exitt


      end subroutine nekMOAB_loadConn
c-----------------------------------------------------------------------------


      subroutine nekMOAB_loadCoord(x27, y27, z27)
      !stuff the xyz coords of the 27 verts of each local element -- 
      !shared vertex coords are stored redundantly
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"
      double precision, dimension(TWENTYSEVEN, hexSize) :: x27, y27, z27


      IMESH_HANDLE vtxHandles(TWENTYSEVEN)
      pointer (vtxHandlesPointer, vtxHandles)
      integer vtxHandlesSize, vtxHandlesAlloc
      integer vtxHandlesData(TWENTYSEVEN)

      integer i, j, ierr

      !!!!!!!!!!!!!!!
      vtxHandlesPointer = %LOC(vtxHandlesData)
      vtxHandlesAlloc = TWENTYSEVEN



      do i=1, hexSize         
         !slight cut and paste from above
         call iMesh_getEntAdj(%VAL(imesh), 
     $                        %VAL(hexHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        %VAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
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
            call iMesh_getVtxCoord(%VAL(imesh),
     $                             %VAL(vtxHandles(j)), !/*in*/ const iBase_EntityHandle vertex_handle,
     $                             x27(j,i), y27(j,i), z27(j,i), !/*out*/ double *x, /*out*/ double *y, /*out*/ double *z,
     $                             ierr)
            IMESH_ASSERT(ierr, imesh)
         enddo

!         print *, i, 'coords: ' 
!         print *, x27
!         print *, y27
!         print *, z27
      enddo



      end subroutine nekMOAB_loadCoord
     
      subroutine nekMOAB_BC(bcdata)
      !Copy the boundary condition information into bcdata array
      !bcdata(s,e) is the bc value of side s of local element e -- default value is 0
      implicit none

#include "iMesh_f.h"
#include 'NEKMOAB'

      integer bcdata(SIX, hexSize)

      IMESH_HANDLE entSetHandles(*)
      pointer (entSetHandlesPointer, entSetHandles)
      integer :: entSetAllocated, entSetSize

      IMESH_HANDLE neuSetTag
      integer:: ierr, i, tagIntData

      bcdata = 0

      !Sidesets in cubit come in as entity sets with the NEUMANN_SET -- see sample file
      call iMesh_getTagHandle(%VAL(imesh),
     $     "NEUMANN_SET", !/*in*/ const char* tag_name,
     $     neuSetTag, !/*out*/ iBase_TagHandle *tag_handle, 
     $     ierr)
      IMESH_ASSERT(ierr, imesh)


      entSetHandlesPointer = IMESH_NULL
      entSetAllocated      = 0
      call iMesh_getEntSets(%VAL(imesh),
     $     %VAL(IMESH_NULL),     !/*in*/ const iBase_EntitySetHandle entity_set_handle,
     $     %VAL(1),              !/*in*/ const int num_hops,
     $     entSetHandlesPointer, !/*out*/ iBase_EntitySetHandle** contained_set_handles,
     $     entSetAllocated,      !/*out*/ int* contained_set_handles_allocated,
     $     entSetSize,           !/*out*/ int* contained_set_handles_size,
     $     ierr)                 !/*out*/ int *err);
      IMESH_ASSERT(ierr, imesh)

      do i=1, entSetSize
         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(entSetHandles(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(neuSetTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        tagIntData,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);

         if (ierr .eq. 0) then !tag was defined
            call nekMOAB_intBC(bcdata, entSetHandles(i), tagIntData)
         endif
      enddo

      call free(entSetHandlesPointer)
      

!      print *, 'bcdata'
!      print *, bcdata

      end subroutine nekMOAB_BC


      subroutine nekMOAB_intBC(bcdata, setHandle, setId)
      !Internal function, don't call directly
      implicit none

#include "iMesh_f.h"
#include "NEKMOAB"

      integer bcdata(SIX, hexSize)

      IMESH_HANDLE :: setHandle
      integer :: setId !coming from cubit

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

      call iMesh_MBCNType(%VAL(iMesh_HEXAHEDRON), hexMBCNType)

      !what faces are in this set?
      facesAlloc = 0
      facesPointer = IMESH_NULL
      call iMesh_getEntitiesRec(%VAL(imesh), 
     $     %VAL(setHandle), 
     $     %VAL(iBase_FACE), 
     $     %VAL(iMesh_QUADRILATERAL),
     $     %VAL(1),
     $     facesPointer, facesAlloc, facesSize, ierr)
      IMESH_ASSERT(ierr, imesh)

      do i=1, facesSize
         !get vertices defining the face
         fvtxAlloc = 0
         fvtxPointer = IMESH_NULL
         call iMesh_getEntAdj(%VAL(imesh),
     $                        %VAL(faces(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        %VAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $                        fvtxPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                        fvtxAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                        fvtxSize , !/*out*/ int* adj_entity_handles_size,
     $                        ierr)
         IMESH_ASSERT(ierr, imesh)         

         !get hexes adjacent to the face (1 or 2) 
         ! -- hopefully only returns hexes on local proc, but untested
         ahexAlloc = 0
         ahexPointer = IMESH_NULL                  
         call iMesh_getEntAdj(%VAL(imesh), 
     $                        %VAL(faces(i)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                        %VAL(iBase_REGION), !/*in*/ const int entity_type_requested,
     $                        ahexPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                        ahexAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                        ahexSize , !/*out*/ int* adj_entity_handles_size,
     $                        ierr)
         IMESH_ASSERT(ierr, imesh)

         do j=1, ahexSize                                
            !get verts adjacent to the hex
            hvtxAlloc = 0
            hvtxPointer = IMESH_NULL                  
            call iMesh_getEntAdj(%VAL(imesh), 
     $                           %VAL(ahex(j)), !/*in*/ const iBase_EntityHandle entity_handle,
     $                           %VAL(iBase_VERTEX), !/*in*/ const int entity_type_requested,
     $                           hvtxPointer, !/*inout*/ iBase_EntityHandle** adj_entity_handles,
     $                           hvtxAlloc, !/*inout*/ int* adj_entity_handles_allocated,
     $                           hvtxSize , !/*out*/ int* adj_entity_handles_size,
     $                           ierr)
            IMESH_ASSERT(ierr, imesh)

            !Get the side number
            call MBCN_SideNumberUlong(%VAL(hvtxPointer),
     $                                %VAL(hexMBCNType), 
     $                                %VAL(fvtxPointer), 
     $                                %VAL(4), !not fvtxSize, because that's 5 and crashes stuff 
                                               ! -- linear nodes are stored first and that's all that's used
     $                                %VAL(2), !Dimension of faces
     $                    side_no, side_sense, side_offset) !output
           if ((side_no .lt. 0) .or. (side_no .gt. 5)) ierr = 1
           IMESH_ASSERT(ierr, imesh)

           side_no = side_no + 1 !moab side number is zero based

           call nekMOAB_getElNo(ahex(j), elno)

           bcdata(side_no, elno) = setId

           call free(hvtxPointer)
         enddo


         call free(ahexPointer)
         call free(fvtxPointer)

      enddo

      call free(facesPointer)

      end subroutine nekMOAB_intBC


      subroutine nekMOAB_getElNo(handle, elno)
      !Given the imesh handle of an element (hex), 
      !   return its local nek element number in elno

      !Cannot be used until the GLLEL array has been set (in map2.f)

      include 'SIZE'
      include 'PARALLEL'
#include "NEKMOAB"
      IMESH_HANDLE handle
      integer elno
      
      integer ierr

      !first get global id
      call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(handle), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(globalIdTag), !/*in*/ const iBase_TagHandle tag_handle,
     $        elno,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);
      IMESH_ASSERT(ierr, imesh)
      
      if (GLLNID(elno) .ne. NID) then
         print *, 'bug... fixme',elno, nid, gllnid(elno)
         !call exitt
      endif
                          
      !then convert to local id
      elno = GLLEL(elno)

      end subroutine nekMOAB_getElNo

 
!code below is not used
!!!!!!!!!!!!!!!!!!!
#if 0
      subroutine printTags(imesh, entHandle)
      implicit none

#include "iMesh_f.h"
      IMESH_HANDLE :: imesh, entHandle
      IMESH_HANDLE tagHandles(*)
      pointer (tagHandlesPointer, tagHandles)
      integer :: tagsAllocated, tagsSize
      integer :: j, ierr
      integer tagintdata
      character(80) :: tempString

      tagHandlesPointer = IMESH_NULL
      tagsAllocated = 0
      call iMesh_getAllTags(%VAL(imesh), 
     $     %VAL(entHandle), !/*in*/ const iBase_EntityHandle entity_handle,                                                        
     $     tagHandlesPointer,   !/*inout*/ iBase_TagHandle** tag_handles,                                                                    
     $     tagsAllocated,       !/*inout*/ int* tag_handles_allocated,                                                                              
     $     tagsSize,            !/*out*/ int* tag_handles_size,                                                                                     
     $     ierr)
      IMESH_ASSERT(ierr, imesh)
         
      print *,  tagsAllocated, "Tags"


      do j=1, tagsSize
         call iMesh_getTagName(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(tagHandles(j)), !/*in*/ const iBase_TagHandle tag_handle,                                                                   
     $        tempstring,       !char *name,                                                                                                        
     $        ierr)
         IMESH_ASSERT(ierr, imesh)
         print *, 'tag number', j, tagHandles(j)
         print *, 'tag name ', IMESH_NULLSTRIP(tempstring)


         call iMesh_getIntData(%VAL(imesh), !iMesh_Instance instance,
     $        %VAL(entHandle), !/*in*/ const iBase_EntityHandle entity_handle,
     $        %VAL(tagHandles(j)), !/*in*/ const iBase_TagHandle tag_handle,
     $        tagIntData,       !/*out*/ int *out_data,
     $        ierr)             !/*out*/ int *err);                                                                                                                              
         IMESH_ASSERT(ierr, imesh)
         print *, 'tag value', tagIntData

      enddo


      end subroutine printTags


#endif     

#if 0
      !snippets
      call iMesh_getNumOfType(%VAL(imesh), 
     $     %VAL(IMESH_NULL),    !/*in*/ const iBase_EntitySetHandle entity_set_handle,                                                            
     $     %VAL(iBase_VERTEX),  !/*in*/ const int entity_type,                                                                                  
     $     numEls,           !/*out*/ int *num_type,                                                                                                  
     $     ierr)

      IMESH_ASSERT(ierr, imesh)
      print *, 'Num of verts: ', numEls

      call iMesh_getNumEntSets(%VAL(imesh),
     $     %VAL(IMESH_NULL), !/*in*/ const iBase_EntitySetHandle entity_set_handle,
     $     %VAL(1), !/*in*/ const int num_hops,
     $     numSets, !/*out*/ int *num_sets,
     $     ierr)

      IMESH_ASSERT(ierr, imesh)
      print *, 'Num of sets: ', numSets

#endif

c-----------------------------------------------------------------------
      subroutine moab_dat
c
      include 'SIZE'
      include 'TOTAL'

      common /cmoab/ moabbc(2*ldim,lelt)

      integer f


c     Use MOAB for :
c
c     .get element to processor map
c     .get vertices
c     .get 27-pt. hex mesh
c     .get bcs

c     call nekMOAB_load(cubfle)
c     call nekMOAB_load('pin217.h5m')
      call nekMOAB_load('input.cub')


      call nekMOAB_proc_map()

      nelgt = nelgv
      call mapelpr()
c
c     Read std nek stuff
c

      open (unit=9,file=reafle,status='old')
      call rdparam !     Parameters

      call moab_geometry (xm1,ym1,zm1,nx1,ny1,nz1)  ! fill xm1
      call xml2xc        (xm1,ym1,zm1,nx1,ny1,nz1)  ! fill xc

      call nekMOAB_BC(moabbc)
      call moab_to_nek_bc       ! fill the cbc arrays


      call nekMOAB_loadMaterialSets


#if 0
      do f=1,6
         call outface(xm1,1,f,nx1,ny1,nz1,'xface ')
         call outface(ym1,1,f,nx1,ny1,nz1,'yface ')
         call outface(zm1,1,f,nx1,ny1,nz1,'zface ')
      enddo
      call exitt
#endif

c     remaining nek stuff
c

      call rdicdf  !     Initial Conditions / Drive Force
      call rdmatp  !     Materials property data
      call rdhist  !     History data
      call rdout   !     Output specs
      call rdobj   !     Objects
      close(9)     !     End of input data, close read file.

      return
      end
c-----------------------------------------------------------------------
      subroutine set_char_bc(cbi,ifld)

      include 'SIZE'
      include 'TOTAL'

      character*3 cbi(0:100)

      if (ifld.eq.1) then
         cbi(0) = 'E  '
         cbi(1) = 'V  '
         cbi(2) = 'v  '
         cbi(3) = 'W  '
         cbi(4) = 'SYM'
         cbi(5) = 'O  '
         cbi(6) = 'o  '
      elseif (ifld.eq.2) then
         cbi(0) = 'E  '
         cbi(1) = 'T  '
         cbi(2) = 't  '
         cbi(3) = 'I  '
         cbi(4) = 'O  '
         cbi(5) = 'o  '
         cbi(6) = 'F  '
         cbi(7) = 'f  '
         cbi(8) = 'C  '
         cbi(9) = 'c  '
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine moab_to_nek_bc  ! fill the nek cbc arrays
c
      include 'SIZE'
      include 'TOTAL'

      common /cmoab/ moabbc(2*ldim,lelt)

      integer e,f
      character*3 cbi(0:100,0:ldimt1)

      ifld = 1  ! fluid only, for now
      ifld = 2  ! heat only, for now

      call set_char_bc(cbi(1,ifld),ifld)

      call rzero(bc,30*lelt*(ldimt+1))

      nface = 2*ndim
      do e=1,nelt
      do f=1,nface
         cbc(f,e,ifld) = cbi(moabbc(f,e),ifld)
!         write(6,1) e,f,moabbc(f,e),cbc(f,e,ifld)
      enddo
      enddo
   1  format(3i8,2x,a3,'  moab cbc')

      return
      end
c-----------------------------------------------------------------------
      subroutine moab_geometry (xml,yml,zml,nxl,nyl,nzl)

      include 'SIZE'
      include 'TOTAL'

      real xml(nxl*nyl*nzl,nelt)
     $   , yml(nxl*nyl*nzl,nelt)
     $   , zml(nxl*nyl*nzl,nelt)
      integer e

      common /tcrmg/ x27(27,lelt), y27(27,lelt), z27(27,lelt)

      call nekMOAB_loadCoord(x27, y27, z27) !     Get coords from moab

c     call outmat(x27,27,8,'x27dat',nelt)
c     call outmat(y27,27,8,'y27dat',nelt)
c     call outmat(z27,27,8,'z27dat',nelt)

      nmoab = 3

      do e=1,nelt   !  Interpolate for each element
         call permute_and_map(xml(1,e),x27(1,e),nx1,nmoab,e)
         call permute_and_map(yml(1,e),y27(1,e),nx1,nmoab,e)
         if (if3d) call permute_and_map(zml(1,e),z27(1,e),nx1,nmoab,e)
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
      subroutine xml2xc (xml,yml,zml,nxl,nyl,nzl)

      include 'SIZE'
      include 'INPUT'

      real xml(nxl,nyl,nzl,nelt)
     $   , yml(nxl,nyl,nzl,nelt)
     $   , zml(nxl,nyl,nzl,nelt)

      integer e

      l = 0
      if (if3d) then
         do e=1,nelt
         do k=1,nz1,nz1-1
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,1) = xml(i,j,k,e)
            yc(l,1) = yml(i,j,k,e)
            zc(l,1) = zml(i,j,k,e)
         enddo
         enddo
         enddo
         enddo
      else  ! 2D
         do e=1,nelt
         do j=1,ny1,ny1-1
         do i=1,nx1,nx1-1
            l = l+1
            xc(l,1) = xml(i,j,1,e)
            yc(l,1) = yml(i,j,1,e)
         enddo
         enddo
         enddo
      endif

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
