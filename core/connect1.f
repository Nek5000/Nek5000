c-----------------------------------------------------------------------
      subroutine connect
C
C     Parallel compatible routine to find 
C     connectivity of element structure.
C
C     On Processor 0:
C
C     .Verify right-handedness of elements.
C     .Verify element-to-element reciprocity of BC's
C     .Verify correlation between E-E BC's and physical coincidence
C     .Set rotations
C     .Determine multiplicity
C     .Set up direct stiffness summation arrays.
C
C     All Processors:
C
C     .Disperse/Receive BC and MULT temporary data read from preprocessor.
C
C
      INCLUDE 'SIZE'
      INCLUDE 'TOTAL'
      INCLUDE 'NONCON'
      INCLUDE 'ZPER'
      INCLUDE 'SCRCT'
c
      COMMON /SCRUZ/ XM3 (LX3,LY3,LZ3,LELT)
     $ ,             YM3 (LX3,LY3,LZ3,LELT)
     $ ,             ZM3 (LX3,LY3,LZ3,LELT)
C
      common /c_is1/ glo_num(1*lx1*ly1*lz1*lelv)
      common /ivrtx/ vertex ((2**ldim)*lelg)
      integer glo_num,vertex
C
C     Initialize key arrays for Direct Stiffness SUM.
C
      NXL=3
      NYL=3
      NZL=1+2*(NDIM-2)

      call initds
      call dsset (nx1,ny1,nz1)
      call setedge
C
C=================================================
C     Establish (global) domain topology
C=================================================
C
C     .Generate topologically correct mesh data.
C     .Set up element centers, face centers, etc. 
C     .Check  right handedness of elements.
C     .Check  element boundary conditions.
C     .Establish Element-Element rotations
C     .Construct the element to processor map and

      call genxyzl
      call setside
      call verify

      CALL SETCDOF 
      IF (IFAXIS            ) CALL SETRZER
      IF (IFMVBD            ) CALL CBCMESH
      IF (IFMODEL.AND.IFKEPS) CALL CBCTURB
      CALL CHKAXCB

C
C========================================================================
C     Set up multiplicity and direct stiffness arrays for each IFIELD
C========================================================================
C
      mfield=2
      if (ifflow) mfield=1
      if (ifmvbd) mfield=0

      ncrnr = 2**ndim

      if (nelgv.eq.nelgt) then
         if (ifgtp) then
            call gen_gtp_vertex (vertex, ncrnr)
         else
            call f77_get_vert
         endif
         call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)
         gsh_fld(2)=gsh_fld(1)

      else

c
c        For conjugate heat transfer, it is assumed that fluid
c        elements are listed both globally and locally with lower
c        element numbers than the solid elements.
c
c        We currently assume that there is at least one fluid elem.
c        per processor.
c

         call f77_get_vert
c        call outmati(vertex,4,nelv,'vrtx V')
         call setupds(gsh_fld(1),nx1,ny1,nz1,nelv,nelgv,vertex,glo_num)

c        call f77_get_vert  (vertex, ncrnr, nelgt, '.mp2')  !  LATER !
c        call outmati(vertex,4,nelt,'vrtx T')
         call setupds(gsh_fld(2),nx1,ny1,nz1,nelt,nelgt,vertex,glo_num)

      endif

      ntotv = nx1*ny1*nz1*nelv
      ntott = nx1*ny1*nz1*nelt

      if (ifflow) then
         ifield = 1
         call rone    (vmult,ntotv)
         call dssum   (vmult,nx1,ny1,nz1)
         call invcol1 (vmult,ntotv)
      endif
      if (ifheat) then
         ifield = 2
         call rone    (tmult,ntott)
         call dssum   (tmult,nx1,ny1,nz1)
         call invcol1 (tmult,ntott)
      endif
      if (.not.ifflow) call copy(vmult,tmult,ntott)
      if (ifmvbd)  call copy (wmult,vmult,ntott)
      do ifield=3,nfield                  ! Additional pass. scalrs.
         if (nelg(ifield).eq.nelgv) then
            gsh_fld(ifield) = gsh_fld(1)
            call copy (tmult(1,1,1,1,ifield-1),vmult,ntotv)
         else
            gsh_fld(ifield) = gsh_fld(2)
            call copy (tmult(1,1,1,1,ifield-1),tmult,ntott)
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine initds
C
C          -- Direct Stiffness Initialization Routine --
C
C     Set up required data for packing data on faces of spectral cubes.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
C     Nominal ordering for direct stiffness summation of faces
C
      J=0
      DO 5 IDIM=1,NDIM
      DO 5 IFACE=1,2
        J=J+1
         NOMLIS(IFACE,IDIM)=J
    5 CONTINUE
C
C     Assign Ed's numbering scheme to PF's scheme.
C
      EFACE(1)=4
      EFACE(2)=2
      EFACE(3)=1
      EFACE(4)=3
      EFACE(5)=5
      EFACE(6)=6
C
C     Assign inverse of Ed's numbering scheme to PF's scheme.
C
      EFACE1(1)=3
      EFACE1(2)=2
      EFACE1(3)=4
      EFACE1(4)=1
      EFACE1(5)=5
      EFACE1(6)=6
C
C     Assign group designation to each face to determine ordering of indices.
C
      GROUP(1)=0
      GROUP(2)=1
      GROUP(3)=1
      GROUP(4)=0
      GROUP(5)=0
      GROUP(6)=1
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine setedge
C
C     .Initialize EDGE arrays for face and edge specific tasks.
C
C     .NOTE: Sevaral arrays in common are initialized via 
C            BLOCKDATA EDGEC
C
C     Computed arrays: 
C
C     IEDGE  -  Minimal list of wire frame nodes.
C               Used to search for all physical
C               coincidences.
C
C
C     IEDGEF - .Ordered list of wire frame nodes 
C               associated with faces 1 through 6.
C              .Each of 4 sides of square frame stored 
C               individually so that rotations are 
C               readily handled.
C              .Two types of node orderings stored -
C               (0) is clockwise marching
C               (1) is counter-clockwise marching 
C                   for image face.
C                               
C
C     IFACE         - indicates the face number.  Two notations
C                     are currently in use:
C
C                     i) Preprocessor notation:
C                
C                                       +--------+     ^ S
C                                      /        /|     |
C                                     /    3   / |     |
C                               4--> /        /  |     |
C                                   +--------+ 2 +     +----> R
C                                   |        |  /     /
C                                   |    6   | /     /
C                                   |        |/     /
C                                   +--------+     T
C                                       1
C
C                    ii) Symmetric notation:
C                     
C                                       +--------+     ^ S
C                                      /        /|     |
C                                     /    4   / |     |
C                               1--> /        /  |     |
C                                   +--------+ 2 +     +----> R
C                                   |        |  /     /
C                                   |    6   | /     /
C                                   |        |/     /
C                                   +--------+     T
C                                       3
C
C     EFACE(IFACE)  -   Given face number IFACE in symmetric notation, 
C                       returns preprocessor notation face number.
C
C     EFACE1(IFACE) -   Given face number IFACE in preprocessor notation, 
C                       returns symmetric notation face number.
C
C  The following variables all take the symmetric 
C  notation of IFACE as arguments:
C
C     ICFACE(i,IFACE) -   Gives the 4 vertices which reside on face IFACE
C                         as depicted below, e.g. ICFACE(i,2)=2,4,6,8.
C
C                        3+-----+4    ^ Y
C                        /  2  /|     |
C     Edge 1 extends    /     / |     |
C       from vertex   7+-----+8 +2    +----> X
C       1 to 2.        |  4  | /     /
C                      |     |/     /
C                     5+-----+6    Z
C                         3
C
C     IEDGFC(i,IFACE) -   Gives the 4 edges which border the face IFACE
C                         Edge numbering is as follows:
C                            Edge = 1,2,3,4     run in +r direction
C                            Edge = 5,6,7,8     run in +s direction
C                            Edge = 9,10,11,12  run in +t direction
C
C                         Ordering of each edge is such that a monotonically
C                         increasing sequence of vertices is associated with
C                         the start point of a corresponding set of 
C                         monotonically increasing edge numbers, e.g.,
C
C     ICEDG(i,IEDGE)  -   Gives 3 variables for determining the stride along
C                         a given edge, IEDGE;  i=1 gives the starting vertex
C                                               i=2 gives the stopping vertex
C                                               i=3 gives the stride size.
C
      INCLUDE 'SIZE'
      INCLUDE 'TOPOL'
C
      COMMON /CTMP0/ ITMP(3,3,3)
      INTEGER ORDER
C
      NXL=3
      NYL=3
      NZL=1+2*(NDIM-2)
      NXY   =NXL*NYL
      NXYZ  =NXL*NYL*NZL
      NFACES=2*NDIM
C
C----------------------------------------------------------------------
C     Set up edge arrays (temporary - required only for defining DS)
C----------------------------------------------------------------------
C
C     Fill corners - 1 through 8.
C
      I3D=1
      IF (NDIM.EQ.2) I3D=0
C
      I=0
      DO 10 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 10 I2=0,1
            IY=1+(NYL-1)*I2
            DO 10 I1=0,1
               IX=1+(NXL-1)*I1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
   10 CONTINUE
C
C     Fill X-direction edges.
C
      DO 20 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 20 I2=0,1
            IY=1+(NYL-1)*I2
            DO 20 IX=2,NXL-1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
   20 CONTINUE
C
C     Fill Y-direction edges.
C
      DO 30 I3=0,I3D
         IZ=1+(NZL-1)*I3
         DO 30 I1=0,1
            IX=1+(NXL-1)*I1
            DO 30 IY=2,NYL-1
               I=I+1
               IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
   30 CONTINUE
C
C     Fill Z-direction edges.
C
      IF (NDIM.EQ.3) THEN
         DO 40 I2=0,1
            IY=1+(NYL-1)*I2
            DO 40 I1=0,1
               IX=1+(NXL-1)*I1
               DO 40 IZ=2,NZL-1
                  I=I+1
                  IEDGE(I)=IX+NXL*(IY-1)+NXY*(IZ-1)
   40    CONTINUE
      ENDIF
C
      CALL IZERO(INVEDG,27)
      DO 44 II=1,20
         IX=IEDGE(II)
         INVEDG(IX)=II
   44 CONTINUE
C
C
C     GENERAL FACE, GENERAL ROTATION EDGE NUMBERS.
C
      IF (NDIM.EQ.3) THEN
C
C        Pack 3-D edge numbering:
C   
C        Fill temporary array with local index numbers:
C
         DO 50 IX=1,NXYZ
            ITMP(IX,1,1)=IX
   50    CONTINUE
C
C        Two sets are required, the base cube and the image cube 
C        which is being summed with it.
C
         DO 1000 IMAGE=0,1
C
C        Pack edges for each face, no rotation.
C
         DO 500 IFACE=1,NFACES
            JS1    = SKPDAT(1,IFACE)
            JF1    = SKPDAT(2,IFACE)
            JSKIP1 = SKPDAT(3,IFACE)
            JS2    = SKPDAT(4,IFACE)
            JF2    = SKPDAT(5,IFACE)
            JSKIP2 = SKPDAT(6,IFACE)
C
C           Choose proper indexing order according to face type and image.
C
            ORDER = (-1)**(GROUP(IFACE)+IMAGE)
            IF (ORDER.EQ.1) THEN
C
C              Forward ordering:
C
C            +-------------+    ^ v1
C            | --------->| |    |
C            | ^    2    | |    +-->
C            | |         | |      v2
C            | |1       3| |
C            | |    4    V |
C            | |<--------- |
C            F-------------I     F is fiducial node.
C
C                                I is location of fiducial node for
C                                     image face.
C
C           Load edge 1:
C
            J=0
            J2=JS2
            DO 100 J1=JS1,JF1-JSKIP1,JSKIP1
               J=J+1
               IEDGEF(J,1,IFACE,IMAGE)=ITMP(J1,J2,1)
  100       CONTINUE
C
C           Load edge 2:
C
            J=0
            J1=JF1
            DO 200 J2=JS2,JF2-JSKIP2,JSKIP2
               J=J+1
               IEDGEF(J,2,IFACE,IMAGE)=ITMP(J1,J2,1)
  200       CONTINUE
C
C
C           Load edge 3:
C
            J=0
            J2=JF2
            DO 300 J1=JF1,JS1+JSKIP1,-JSKIP1
               J=J+1
               IEDGEF(J,3,IFACE,IMAGE)=ITMP(J1,J2,1)
  300       CONTINUE
C
C           Load edge 4:
C
            J=0
            J1=JS1
            DO 400 J2=JF2,JS2+JSKIP2,-JSKIP2
               J=J+1
               IEDGEF(J,4,IFACE,IMAGE)=ITMP(J1,J2,1)
  400       CONTINUE
C
            ELSE
C
C           Reverse ordering:
C
C            +-------------+
C            | |<--------- |       ^ v2
C            | |    2    ^ |       |
C            | |         | |    <--+
C            | |3       1| |     v1
C            | V    4    | |
C            | --------->| |
C            I-------------F     F is fiducial node.
C
C                                I is location of fiducial node for
C                                     image face.
C
C           Load edge 1:
C   
            J=0
            J1=JS1
            DO 105 J2=JS2,JF2-JSKIP2,JSKIP2
               J=J+1
               IEDGEF(J,1,IFACE,IMAGE)=ITMP(J1,J2,1)
  105       CONTINUE
C
C           Load edge 2:
C
            J=0
            J2=JF2
            DO 205 J1=JS1,JF1-JSKIP1,JSKIP1
               J=J+1
               IEDGEF(J,2,IFACE,IMAGE)=ITMP(J1,J2,1)
  205       CONTINUE
C
C           Load edge 3:
C
            J=0
            J1=JF1
            DO 305 J2=JF2,JS2+JSKIP2,-JSKIP2
                J=J+1
               IEDGEF(J,3,IFACE,IMAGE)=ITMP(J1,J2,1)
  305       CONTINUE
C
C           Load edge 4:
C
            J=0
            J2=JS2
            DO 405 J1=JF1,JS1+JSKIP1,-JSKIP1
               J=J+1
              IEDGEF(J,4,IFACE,IMAGE)=ITMP(J1,J2,1)
  405       CONTINUE
            ENDIF
C
  500    CONTINUE
 1000    CONTINUE
      ELSE
C
C        Load edge information for 2-D case
C
         IEDGEF(1,1,1,0) = NXY - NXL + 1
         IEDGEF(1,2,1,0) = 1
         IEDGEF(1,1,2,0) = NXL
         IEDGEF(1,2,2,0) = NXY
         IEDGEF(1,1,3,0) = 1
         IEDGEF(1,2,3,0) = NXL
         IEDGEF(1,1,4,0) = NXY
         IEDGEF(1,2,4,0) = NXY - NXL + 1
C
         IEDGEF(1,1,1,1) = 1
         IEDGEF(1,2,1,1) = NXY - NXL + 1
         IEDGEF(1,1,2,1) = NXY
         IEDGEF(1,2,2,1) = NXL
         IEDGEF(1,1,3,1) = NXL
         IEDGEF(1,2,3,1) = 1
         IEDGEF(1,1,4,1) = NXY - NXL + 1
         IEDGEF(1,2,4,1) = NXY
      ENDIF
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine dsset(nx,ny,nz)
C
C     Set up arrays IXCN,ESKIP,SKPDAT,NEDG,NOFFST for new NX,NY,NZ
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TOPOL'
      INTEGER NXO,NYO,NZO
      SAVE    NXO,NYO,NZO
      DATA    NXO,NYO,NZO /3*0/
C
C     Check if element surface counters are already set from last call...
C
      IF (NXO.EQ.NX.AND.NYO.EQ.NY.AND.NZO.EQ.NZ) RETURN
C
C     else, proceed....
C
      NXO = NX
      NYO = NY
      NZO = NZ
C
C     Establish corner to elemental node number mappings
C
      IC=0
      DO 10 ICZ=0,1
      DO 10 ICY=0,1
      DO 10 ICX=0,1
C       Supress vectorization to 
        IF(ICX.EQ.0)DUMMY=0        
        IF(ICX.EQ.1)DUMMY=1
        DUMMY2=DUMMY2+DUMMY
        IC=IC+1
        IXCN(IC)= 1 + (NX-1)*ICX + NX*(NY-1)*ICY + NX*NY*(NZ-1)*ICZ
   10   CONTINUE
C
C     Assign indices for direct stiffness summation of arbitrary faces.
C
C
C     Y-Z Planes (Faces 1 and 2)
C
      SKPDAT(1,1)=1
      SKPDAT(2,1)=NX*(NY-1)+1
      SKPDAT(3,1)=NX
      SKPDAT(4,1)=1
      SKPDAT(5,1)=NY*(NZ-1)+1
      SKPDAT(6,1)=NY
C
      SKPDAT(1,2)=1             + (NX-1)
      SKPDAT(2,2)=NX*(NY-1)+1   + (NX-1)
      SKPDAT(3,2)=NX
      SKPDAT(4,2)=1
      SKPDAT(5,2)=NY*(NZ-1)+1
      SKPDAT(6,2)=NY
C
C     X-Z Planes (Faces 3 and 4)
C
      SKPDAT(1,3)=1
      SKPDAT(2,3)=NX
      SKPDAT(3,3)=1
      SKPDAT(4,3)=1
      SKPDAT(5,3)=NY*(NZ-1)+1
      SKPDAT(6,3)=NY
C
      SKPDAT(1,4)=1           + NX*(NY-1)
      SKPDAT(2,4)=NX          + NX*(NY-1)
      SKPDAT(3,4)=1
      SKPDAT(4,4)=1
      SKPDAT(5,4)=NY*(NZ-1)+1
      SKPDAT(6,4)=NY
C
C     X-Y Planes (Faces 5 and 6)
C
      SKPDAT(1,5)=1
      SKPDAT(2,5)=NX
      SKPDAT(3,5)=1
      SKPDAT(4,5)=1
      SKPDAT(5,5)=NY
      SKPDAT(6,5)=1
C
      SKPDAT(1,6)=1           + NX*NY*(NZ-1)
      SKPDAT(2,6)=NX          + NX*NY*(NZ-1)
      SKPDAT(3,6)=1
      SKPDAT(4,6)=1
      SKPDAT(5,6)=NY
      SKPDAT(6,6)=1
C
C     Set up skip indices for each of the 12 edges 
C
C         Note that NXY = NX*NY even for 2-D since
C         this branch does not apply to the 2D case anyway.
C
C     ESKIP(*,1) = start location
C     ESKIP(*,2) = end 
C     ESKIP(*,3) = stride
C
      NXY=NX*NY
      ESKIP( 1,1) = IXCN(1) + 1
      ESKIP( 1,2) = IXCN(2) - 1
      ESKIP( 1,3) = 1
      ESKIP( 2,1) = IXCN(3) + 1
      ESKIP( 2,2) = IXCN(4) - 1
      ESKIP( 2,3) = 1
      ESKIP( 3,1) = IXCN(5) + 1
      ESKIP( 3,2) = IXCN(6) - 1
      ESKIP( 3,3) = 1
      ESKIP( 4,1) = IXCN(7) + 1
      ESKIP( 4,2) = IXCN(8) - 1
      ESKIP( 4,3) = 1
      ESKIP( 5,1) = IXCN(1) + NX
      ESKIP( 5,2) = IXCN(3) - NX
      ESKIP( 5,3) = NX
      ESKIP( 6,1) = IXCN(2) + NX
      ESKIP( 6,2) = IXCN(4) - NX
      ESKIP( 6,3) = NX
      ESKIP( 7,1) = IXCN(5) + NX
      ESKIP( 7,2) = IXCN(7) - NX
      ESKIP( 7,3) = NX
      ESKIP( 8,1) = IXCN(6) + NX
      ESKIP( 8,2) = IXCN(8) - NX
      ESKIP( 8,3) = NX
      ESKIP( 9,1) = IXCN(1) + NXY
      ESKIP( 9,2) = IXCN(5) - NXY
      ESKIP( 9,3) = NXY
      ESKIP(10,1) = IXCN(2) + NXY
      ESKIP(10,2) = IXCN(6) - NXY
      ESKIP(10,3) = NXY
      ESKIP(11,1) = IXCN(3) + NXY
      ESKIP(11,2) = IXCN(7) - NXY
      ESKIP(11,3) = NXY
      ESKIP(12,1) = IXCN(4) + NXY
      ESKIP(12,2) = IXCN(8) - NXY
      ESKIP(12,3) = NXY
C
C     Load reverse direction edge arrays for reverse mappings...
C
      DO 20 IED=1,12
      IEDM=-IED
      ESKIP(IEDM,1) =  ESKIP(IED,2)
      ESKIP(IEDM,2) =  ESKIP(IED,1)
      ESKIP(IEDM,3) = -ESKIP(IED,3)
   20 CONTINUE
C
C     Compute offset for global edge vector given current element
C     dimensions.....
C
C     NGSPED(ITE,ICMP) = number of global (ie, distinct) special edges
C                        of type ITE (1,2, or 3)  for field ICMP.
C
C                        ITE = 1 implies an "X" edge
C                        ITE = 2 implies an "Y" edge
C                        ITE = 3 implies an "Z" edge
C
C     Set up number of nodes along each of the 3 types of edges
C     (endpoints excluded).
C
      NEDG(1)=NX-2
      NEDG(2)=NY-2
      NEDG(3)=NZ-2
C
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine genxyzl
C
C     Generate xyz coordinates 
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'SCRCT'
      COMMON /CTMP0/ XCB(2,2,2),YCB(2,2,2),ZCB(2,2,2),H(3,3,2),INDX(8)
C
      NXL=3
      NYL=3
      NZL=1+2*(NDIM-2)
      NTOT3=NXL*NYL*NZL*NELT
C
C   Preprocessor Corner notation:      Symmetric Corner notation:
C                
C           4+-----+3    ^ s                    3+-----+4    ^ s        
C           /     /|     |                      /     /|     |          
C          /     / |     |                     /     / |     |          
C        8+-----+7 +2    +----> r            7+-----+8 +2    +----> r   
C         |     | /     /                     |     | /     /           
C         |     |/     /                      |     |/     /            
C        5+-----+6    t                      5+-----+6    t             
C
      DO 10 IX=1,NXL
         H(IX,1,1)=0.5*FLOAT(3-IX)
         H(IX,1,2)=0.5*FLOAT(IX-1)
   10 CONTINUE
      DO 20 IY=1,NYL
         H(IY,2,1)=0.5*FLOAT(3-IY)
         H(IY,2,2)=0.5*FLOAT(IY-1)
   20 CONTINUE
      DO 30 IZ=1,NZL
         H(IZ,3,1)=0.5*FLOAT(3-IZ)
         H(IZ,3,2)=0.5*FLOAT(IZ-1)
   30 CONTINUE
C
      INDX(1)=1
      INDX(2)=2
      INDX(3)=4
      INDX(4)=3
      INDX(5)=5
      INDX(6)=6
      INDX(7)=8
      INDX(8)=7
C
      CALL RZERO(XML,NTOT3)
      CALL RZERO(YML,NTOT3)
      CALL RZERO(ZML,NTOT3)
      CALL RZERO(XCB,8)
      CALL RZERO(YCB,8)
      CALL RZERO(ZCB,8)
C
      DO 5000 IE=1,NELT
C
         NDIM2 = 2**NDIM
         DO 50 IX=1,NDIM2
            I=INDX(IX)
            XCB(IX,1,1)=XC(I,IE)
            YCB(IX,1,1)=YC(I,IE)
            ZCB(IX,1,1)=ZC(I,IE)
   50    CONTINUE
C
C        Map R-S-T space into physical X-Y-Z space.
C
         DO 100 IZT=1,ndim-1
         DO 100 IYT=1,2
         DO 100 IXT=1,2
C
         DO 100 IZ=1,NZL
         DO 100 IY=1,NYL
         DO 100 IX=1,NXL
            XML(IX,IY,IZ,IE)=XML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*XCB(IXT,IYT,IZT)
            YML(IX,IY,IZ,IE)=YML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*YCB(IXT,IYT,IZT)
            ZML(IX,IY,IZ,IE)=ZML(IX,IY,IZ,IE)+
     $        H(IX,1,IXT)*H(IY,2,IYT)*H(IZ,3,IZT)*ZCB(IXT,IYT,IZT)
  100    CONTINUE
C
 5000 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verify
C
C     .Verify right-handedness of elements.
C     .Verify element-to-element reciprocity of BC's
C     .Verify correlation between E-E BC's and physical coincidence
C
      include 'SIZE'
      include 'PARALLEL'
      include 'SCRCT'
C
      call verrhe
      if (np.eq.0) then
c        Note this routine *can* be called in parallel.  It is simply
c        turned off to reduce start up time.  pff 2/7/98
         call verbcs
      elseif (nid.eq.0) then
         write(6,*) 'NOTE:  No more E-E BC verification in parallel'
      endif
C
      return
      end
c-----------------------------------------------------------------------
      subroutine setside
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'TOPOL'
      INCLUDE 'SCRCT'
C
C     SIDE(i,IFACE,IE) -  Physical (xyz) location of element side midpoint.
C                         i=1,2,3 gives x,y,z value, respectively.
C                         i=4  gives average dimension of face for setting
C                              tolerances.
C
      INDX(1)=1
      INDX(2)=2
      INDX(3)=4
      INDX(4)=3
      INDX(5)=5
      INDX(6)=6
      INDX(7)=8
      INDX(8)=7
C
C     Flip vertex array structure
C
c     write(6,*) nelv,nelt,if3d
      call rzero(xyz,24*nelt)
      if (if3d) then
         do ie=1,nelt
         do j=1,8
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
            xyz(3,ivtx,ie) = zc(j,ie)
c           write(6,1) ie,j,ivtx,xc(j,ie),yc(j,ie),zc(j,ie),' xcz'
c           write(6,1) ie,j,ivtx,(xyz(k,ivtx,ie),k=1,3),' vtx'
c   1       format(3i5,1p3e12.4,a4)
         enddo
         enddo
      else
         do ie=1,nelt
         do j=1,4
            ivtx = indx(j)
            xyz(1,ivtx,ie) = xc(j,ie)
            xyz(2,ivtx,ie) = yc(j,ie)
            xyz(3,ivtx,ie) = 0.0
         enddo
         enddo
      endif
C
C     Compute location of center and "diameter" of each element side.
C
      NFACES=NDIM*2
      NCRNR =2**(NDIM-1)
      CALL RZERO(SIDE,24*NELT)
      DO 500 ICRN=1,NCRNR
      DO 500 IFAC=1,NFACES
         IVTX = ICFACE(ICRN,IFAC)
         ICR1 = NCRNR+(ICRN-1)
         ICR1 = MOD1(ICR1,NCRNR)
         IVT1 = ICFACE(ICR1,IFAC)
         DO 400 IE=1,NELT
            DO 300 IDIM=1,NDIM
               SIDE(IDIM,IFAC,IE)=SIDE(IDIM,IFAC,IE)+XYZ(IDIM,IVTX,IE)
               SIDE(   4,IFAC,IE)=SIDE(   4,IFAC,IE)+
     $                     ( XYZ(IDIM,IVTX,IE)-XYZ(IDIM,IVT1,IE) )**2
  300       CONTINUE
            SIDE(4,IFAC,IE)=SQRT( SIDE(4,IFAC,IE) )
  400    CONTINUE
  500 CONTINUE
      AVWGHT=1.0/FLOAT(NCRNR)
      CALL CMULT(SIDE,AVWGHT,24*NELT)
C
c     call exitt
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verrhe
C
C     8 Mar 1989 21:58:26   PFF
C     Verify right-handedness of given elements. 
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      INCLUDE 'TOPOL'
      LOGICAL IFYES,IFCSTT
C
      IFCSTT=.TRUE.
      IF (.NOT.IF3D) THEN
      DO 1000 IE=1,NELT
C
C        CRSS2D(A,B,O) = (A-O) X (B-O)
C
         C1=CRSS2D(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,1,IE))
         C2=CRSS2D(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,2,IE))
         C3=CRSS2D(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,3,IE))
         C4=CRSS2D(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,4,IE))
C
         IF (C1.LE.0.0.OR.C2.LE.0.0.OR.
     $       C3.LE.0.0.OR.C4.LE.0.0 ) THEN
C
            IEG=LGLEL(IE,NODE)
            WRITE(6,800) IEG,C1,C2,C3,C4
            call exitt
  800       FORMAT(/,2X,'WARNINGa: Detected non-right-handed element.',
     $      /,2X,'Number',I8,'  C1-4:',4E12.4)
            IFCSTT=.FALSE.
C           CALL QUERY(IFYES,'Proceed                                 ')
C           IF (.NOT.IFYES) GOTO 9000
         ENDIF
 1000 CONTINUE
C
C     Else 3-D:
C
      ELSE      
      DO 2000 IE=1,NELT
C
C        VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
C
         V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
         V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
         V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
         V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
         V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
         V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
         V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
         V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
C
         IF (V1.LE.0.0.OR.V2.LE.0.0.OR.
     $       V3.LE.0.0.OR.V4.LE.0.0.OR.
     $       V5.LE.0.0.OR.V6.LE.0.0.OR.
     $       V7.LE.0.0.OR.V8.LE.0.0    ) THEN
C
            IEG=LGLEL(IE,NODE)
            WRITE(6,1800) IEG,V1,V2,V3,V4,V5,V6,V7,V8
            call exitt
 1800       FORMAT(/,2X,'WARNINGb: Detected non-right-handed element.',
     $      /,2X,'Number',I8,'  V1-8:',4E12.4
     $      /,2X,'      ',4X,'       ',4E12.4)
            IFCSTT=.FALSE.
         ENDIF
 2000 CONTINUE
      ENDIF
C
 9000 CONTINUE
C
C     Print out results from right-handed check
C
      IF (.NOT.IFCSTT) WRITE(6,2001)
C
C     Check consistency accross all processors.
C
      CALL GLLOG(IFCSTT,.FALSE.)
C
      IF (.NOT.IFCSTT) THEN
         IF (NID.EQ.0) WRITE(6,2003) NELGT
         call exitt
      ELSE
         IF (NID.EQ.0) WRITE(6,2002) NELGT
      ENDIF
C
 2001 FORMAT(//,'  Elemental geometry not right-handed, ABORTING'
     $      ,' in routine VERRHE.')
 2002 FORMAT('  Right-handed check complete for',I8,' elements. OK.')
 2003 FORMAT('  Right-handed check failed for',I8,' elements.'
     $      ,'  Exiting in routine VERRHE.')
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine verbcs
C
C     8 Mar 1989 21:58:26   PFF
C     Verify reciprocity of boundary conditions.
C
      INCLUDE 'SIZE'
      INCLUDE 'INPUT'
      INCLUDE 'PARALLEL'
      INCLUDE 'SCRCT'
      INCLUDE 'TOPOL'
      INCLUDE 'TSTEP'
C
      LOGICAL IFCLOS,IFRECP
      CHARACTER*3 CB1,CBJ,CBCHAR(8,2)
      DATA CBCHAR /'E  ','P  ','MSI','MPI','msi','mpi','MLI','MCI'
     $            ,'E  ','P  ','MPI','MSI','mpi','msi','MCI','MLI'/
      IFCLOS=.TRUE.
      IFRECP=.TRUE.
C
C     In order to check E-E connectivity, we need to know the location
C     of the faces.  Find midpoint of each face or edge by averaging.
C
C     Check Element-Element reciprocity
C
C
C     PARALLEL COMPATIBLE VERSION: 
C       .Set consistency to .FALSE. for each element-face.
C       .Pass each data set around to other processors, checking off
C        each side that is indeed consistent.
C       .Scan list at completion to verify consistency of all faces.
C
      MFIELD=2
      IF (IFFLOW) MFIELD=1
      DO 6000 IFIELD=MFIELD,NFIELD
C
         NFACES=2*NDIM
         NEL=NELFLD(IFIELD)
         NELM=NEL
         DO 10 I=1,NP
            NELM = MAX(NELM,NELP(IFIELD,I))
   10    CONTINUE
         DO 100 IFACE=1,NFACES
         DO 100 IE =1,NEL
         IFCNST(IFACE,IE)=.TRUE.
         CB1 = CBC(IFACE,IE,IFIELD)
         IF (CB1.EQ.'P  ' .OR. CB1.EQ.'E  ' .OR. 
     $       CB1.EQ.'MSI' .OR. CB1.EQ.'msi' .OR.
     $       CB1.EQ.'MPI' .OR. CB1.EQ.'mpi' .OR.
     $       CB1.EQ.'MLI' .OR. CB1.EQ.'MCI') 
     $       IFCNST(IFACE,IE)=.FALSE.
  100    CONTINUE
C
         DO 2000 IRG=1,NP
            NXTNID=IPRING(IRG)
            NXNODE=NXTNID+1
            NEL  = NELFLD(IFIELD)
            NELJ = NELP(IFIELD,NXNODE)
C
C           The following data needs to be passed around the ring
C           in order to verify E-E consistency:
C
            CALL CRING(CBCS,CBC(1,1,IFIELD),18*NELM,IRG)
            CALL RRING(BCS,BC(1,1,1,IFIELD),30*NELM,IRG)
            CALL RRING(SIDES,SIDE,          24*NELM,IRG) 
C
C           See how many IFCNST's we can set to .TRUE. on this pass.
C
            DO 1000 ICHR=1,8
            DO 1000 IE=1,NEL
            DO 1000 IFACE=1,NFACES
               IEG = LGLEL(IE,NODE)
               CB1 = CBC(IFACE,IE,IFIELD)
               IF (.NOT.IFCNST(IFACE,IE).AND.
     $            CB1.EQ.CBCHAR(ICHR,1)) THEN
                  JEG    = INT( BC(1,IFACE,IE,IFIELD) )
                  JFACE  = INT( BC(2,IFACE,IE,IFIELD) )
                  JNID   = GLLNID(JEG)
C
C                 Is the opposite face/element a member of the newly arrived
C                 set of data??
C
                  IF (JNID.EQ.NXTNID) THEN
                     JE   = GLLEL(JEG)
                     CBJ  = CBCS(JFACE,JE)
                     JBC1 = INT( BCS(1,JFACE,JE) )
                     JBC2 = INT( BCS(2,JFACE,JE) )
C
C                    Check for reciprocity of E-E bc.s
C
                     IFRECP=.TRUE.
                     IF (CBJ.NE.CBCHAR(ICHR,2) .OR. JBC1.NE.IEG .OR.
     $                  JBC2.NE.IFACE) THEN
                        IFRECP=.FALSE.
                        IF (NP.GT.1) THEN
                           WRITE(6,510) NID,IEG,IFACE,IEG,IFACE
     $                    ,CB1,JEG,JFACE,JEG,JFACE,CBJ,JBC1,JBC2
                        ELSE
                           WRITE(6,511) IEG,IFACE,IEG,IFACE
     $                    ,CB1,JEG,JFACE,JEG,JFACE,CBJ,JBC1,JBC2
                        ENDIF
  510                   FORMAT(1X,I5,' WARNING: Detected inconsistency'
     $                  ,' for element',I4,', side',I2,':'
     $                  ,/,1X,' IE  SIDE',5X,'CBC',5X,'BC1',5X,'BC2'
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3   )
  511                   FORMAT(1X,' WARNING: Detected inconsistency'
     $                  ,' for element',I4,', side',I2,':'
     $                  ,/,1X,' IE  SIDE',5X,'CBC',5X,'BC1',5X,'BC2'
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3
     $                  ,/,1X,2I5,5X,A3,5X,I3,5X,I3   )
                     ENDIF
C
C                    Check for closeness of sides (w/ semi-check for periodic)
C
                     IFC = EFACE1(IFACE)
                     JFC = EFACE1(JFACE)
                     IF (IF3D) THEN
                        DIST = 
     $                     SQRT( (SIDE(1,IFC,IE)-SIDES(1,JFC,JE))**2 +
     $                           (SIDE(2,IFC,IE)-SIDES(2,JFC,JE))**2 +
     $                           (SIDE(3,IFC,IE)-SIDES(3,JFC,JE))**2 )
                     ELSE
                        DIST = 
     $                     SQRT( (SIDE(1,IFC,IE)-SIDES(1,JFC,JE))**2 +
     $                           (SIDE(2,IFC,IE)-SIDES(2,JFC,JE))**2 )
                     ENDIF
                     EPS =.0001*MIN(SIDE(4,IFC,IE),SIDES(4,JFC,JE))
C
c                    WRITE(7,512) IE,JE,IFC,JFC,IFACE,JFACE,DIST,EPS,
c    $               SIDE(1,IFC,IE),SIDE(1,IFC,IE),SIDE(1,IFC,IE),
c    $               SIDES(1,JFC,JE),SIDES(1,JFC,JE),SIDES(1,JFC,JE)
  512                FORMAT(6I4,2E11.4,/,4X,3E14.4,/,4X,3E14.4)
C
                     IFCLOS=.TRUE.
                     IF (ICHR.NE.2 .AND. DIST.GT.EPS) THEN
                        IFCLOS=.FALSE.
                        IF (NP.GT.1) THEN
                           WRITE(6,610) NID,IEG,IFACE,DIST,EPS
     $                                     ,JEG,JFACE
  610                      FORMAT(2X,'WARNING:  Sides too far '
     $                    ,'apart for internal boundary condition.'
     $                    ,/,2X,' NID   EL  FACE  DIST      TOL.'
     $                    ,/,2X,  I4,   I6,  I5,  G10.3,   G10.3
     $                    ,/,2X,  4X,   I6,  I5                 )
                        ELSE
c                          WRITE(7,611)     IEG,IFACE,DIST,EPS
c    $                                     ,JEG,JFACE
                           WRITE(6,611)     IEG,IFACE,DIST,EPS
     $                                     ,JEG,JFACE
  611                      FORMAT(2X,'WARNING:  Sides too far '
     $                    ,'apart for internal boundary condition.'
     $                    ,/,2X,'   EL  FACE  DIST      TOL.'
     $                    ,/,2X,    I6,  I5,  G10.3,   G10.3
     $                    ,/,2X,    I6,  I5                 )
                        ENDIF
                     ENDIF
C                    consistent?
                     IF (IFRECP.AND.IFCLOS) IFCNST(IFACE,IE)=.TRUE.
C                    End of element-element match for this pass of the RING.
                  ENDIF
C                 End of consistency check for this element.
               ENDIF
 1000       CONTINUE
 2000    CONTINUE
C
C        Ring check is completed, now verify that all BC's for this field
C        are consistent:
C
         DO 4000 IFACE=1,NFACES
         DO 4000 IE =1,NEL
         IF (.NOT.IFCNST(IFACE,IE)) THEN
            IEG=LGLEL(IE,NODE)
            WRITE(6,4001) NID,IFIELD,IEG,IFACE
            call exitt
            RETURN
         ENDIF
 4000    CONTINUE
C     End of fields loop
 6000 CONTINUE
C
C     Formats.
C
 4001 FORMAT(1X,I5,1X
     $   ,'Boundary condition consistency check failed for '
     $   ,'field',I2,', element',I5,', face',I2,'.'
     $   ,/,2X,'ABORTING in routine VERBCS.')
 4002 FORMAT('  Boundary condition consistency check, OK.')
C
C     Must be OK:
C
      IF (NID.EQ.0) WRITE(6,4002)
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION VOLUM0(P1,P2,P3,P0)
C
C                           3
C     Given four points in R , (P1,P2,P3,P0), VOLUM0 returns
C     the volume enclosed by the parallelagram defined by the
C     vectors { (P1-P0),(P2-P0),(P3-P0) }.  This routine has
C     the nice feature that if the 3 vectors so defined are
C     not right-handed then the volume returned is negative.
C
      REAL P1(3),P2(3),P3(3),P0(3)
C
         U1=P1(1)-P0(1)
         U2=P1(2)-P0(2)
         U3=P1(3)-P0(3)
C
         V1=P2(1)-P0(1)
         V2=P2(2)-P0(2)
         V3=P2(3)-P0(3)
C
         W1=P3(1)-P0(1)
         W2=P3(2)-P0(2)
         W3=P3(3)-P0(3)
C
         CROSS1 = U2*V3-U3*V2
         CROSS2 = U3*V1-U1*V3
         CROSS3 = U1*V2-U2*V1
C
         VOLUM0  = W1*CROSS1 + W2*CROSS2 + W3*CROSS3
         
      RETURN
      END
c-----------------------------------------------------------------------
      FUNCTION CRSS2D(XY1,XY2,XY0)
      REAL XY1(2),XY2(2),XY0(2)
C
         V1X=XY1(1)-XY0(1)
         V2X=XY2(1)-XY0(1)
         V1Y=XY1(2)-XY0(2)
         V2Y=XY2(2)-XY0(2)
         CRSS2D = V1X*V2Y - V1Y*V2X
C
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,iface)
       KX1=1
       KY1=1
       KZ1=1
       KX2=NX1
       KY2=NY1
       KZ2=NZ1
       IF (IFACE.EQ.1) KY2=1
       IF (IFACE.EQ.2) KX1=NX1
       IF (IFACE.EQ.3) KY1=NY1
       IF (IFACE.EQ.4) KX2=1
       IF (IFACE.EQ.5) KZ2=1
       IF (IFACE.EQ.6) KZ1=NZ1
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine facev(a,ie,iface,val,nx,ny,nz)
C
C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      INCLUDE 'SIZE'
      DIMENSION A(NX,NY,NZ,LELT)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
         A(IX,IY,IZ,IE)=VAL
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine ifacev(a,ie,iface,val,nx,ny,nz)
C
C     Assign the value VAL to face(IFACE,IE) of array A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      include 'SIZE'
      integer a(nx,ny,nz,lelt),val
      call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx,ny,nz,iface)
      do 100 iz=kz1,kz2
      do 100 iy=ky1,ky2
      do 100 ix=kx1,kx2
         a(ix,iy,iz,ie)=val
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine facec(a,b,ie,iface,nx,ny,nz,nel)
C
C     Copy the face (IFACE) of B to A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      DIMENSION A(NX,NY,NZ,NEL)
      DIMENSION B(NX,NY,NZ,NEL)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
         A(IX,IY,IZ,IE)=B(IX,IY,IZ,IE)
  100 CONTINUE
      RETURN
      END
c-----------------------------------------------------------------------
      subroutine setupds(gs_handle,nx,ny,nz,nel,melg,vertex,glo_num)
      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'NONCON'
      integer gs_handle
      integer vertex(1),glo_num(1)

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      call set_vert(glo_num,ngv,nx,nel,melg,vertex,.false.)

c     Initialize gather-scatter code

      t0 = dclock()
      ntot      = nx*ny*nz*nel

      call gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

c     call gs_chkr(glo_num)

      t1 = dclock()
      et = t1-t0
c
      if (nid.eq.0) then
         write(6,1) et,nx,nel,ntot,ngv,gs_handle
    1    format('gs_init time',1pe11.4,' seconds ',i3,4i10)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine combin2(glnm1,glnm2,nglob)
c
      write(6,*) 'Hey, who called combin2???  ABORT'
      call exitt
c
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldio (x,txt10)
      INCLUDE 'SIZE'
      integer x(lx1,ly1,lz1,lelt)
      character*10 txt10
C
      do ie=1,nelv,2
         do iz=1,nz1,1
            if (iz.eq.1) write(6,106) txt10,iz,ie
            if (iz.gt.1) write(6,107) 
            i1 = ie+1
            do j=ny1,1,-1
               write(6,105) (x(i,j,iz,ie),i=1,nx1)
     $                    , (x(i,j,iz,i1),i=1,nx1)
            enddo
         enddo
      enddo
C
  107 FORMAT(' ')
  105 FORMAT(4i6,20x,4i6)
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','Plane = ',I2,'/',I2,/,
     $           5X,'       X            ')
C
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldi (x,txt10)
      INCLUDE 'SIZE'
      integer x(lx1,ly1,lz1,lelt)
      character*10 txt10
c
      character*6 s(20,20)
c
      if (lx1.ne.4 .or. nelv.gt.3) return
c
      write(6,106) txt10,ie,ie
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','elem. = ',I2,'/',I2,/,
     $           5X,'       X            ')
C
C
      call blank(s,6*20*20)
      do ie=1,3
         if (ie.eq.1) then
            jstart  = 1
            istart  = 1
            istride = 3
         else
            jstart  = 2 + lx1
            istart  = 1
            istride = 1
            if (ie.eq.2) istart = 7
         endif
c
         i=istart
         do iy=ny1,1,-1
            j=jstart
            do ix=1,nx1
               write(s(i,j),6) x(ix,iy,1,ie)
               j=j+1
            enddo
            i=i+istride
         enddo
    6    format(20i6)
      enddo
c
      do i=1,10
         write(6,7) (s(i,l),l=1,j-1)
      enddo
    7 format(20a6)
c
      write(6,*)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldr (x,txt10)
      INCLUDE 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
c
      character*6 s(20,20)
c
      if (lx1.ne.4 .or. nelv.gt.3) return
      write(6,106) txt10,ie,ie
  106 FORMAT(  /,5X,'     ^              ',/,
     $           5X,'   Y |              ',/,
     $           5X,'     |              ',A10,/,
     $           5X,'     +---->         ','elem. = ',I2,'/',I2,/,
     $           5X,'       X            ')
      call blank(s,6*20*20)
C
C
      do ie=1,3
         if (ie.eq.1) then
            jstart  = 1
            istart  = 1
            istride = 3
         else
            jstart  = 2 + lx1
            istart  = 1
            istride = 1
            if (ie.eq.2) istart = 7
         endif
c
         i=istart
         do iy=ny1,1,-1
            j=jstart
            do ix=1,nx1
               write(s(i,j),6) x(ix,iy,1,ie)
               j=j+1
            enddo
            i=i+istride
         enddo
    6    format(f6.2)
      enddo
c
      do i=1,10
         write(6,7) (s(i,l),l=1,j-1)
      enddo
    7 format(20a6)
c
      write(6,*)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine checkit(idum)
      write(6,*) 'continue?'
      read (5,*) idum
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldro (x,txt10,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(nx1,ny1,nz1,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /3/
      if (idum.lt.0) return
c
C
      mtot = nx1*ny1*nz1*nelv
      if (nx1.gt.8.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      nell = nelt
      rnel = nell
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nell-ne+1
      do ie=1,2
         ne = 0
         do k=1,1
            write(6,116) txt10,k,ie,xmin,xmax,istep,time
            write(6,117) 
            do j=ny1,1,-1
              if (nx1.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
              if (nx1.eq.8) write(6,118) ((x(i,j,k,e+l),i=1,nx1),e=0,ne)
            enddo
         enddo
      enddo
C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
  118 FORMAT(8f12.9)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldrp (x,txt10,ichk)
      INCLUDE 'SIZE'
      INCLUDE 'TSTEP'
      real x(nx2,ny2,nz2,lelt)
      character*10 txt10
c
      integer idum,e
      save idum
      data idum /3/
      if (idum.lt.0) return
c
C
      mtot = nx2*ny2*nz2*nelv
      if (nx2.gt.7.or.nelv.gt.16) return
      xmin = glmin(x,mtot)
      xmax = glmax(x,mtot)
c
      rnel = nelv
      snel = sqrt(rnel)+.1
      ne   = snel
      ne1  = nelv-ne+1
      do ie=ne1,1,-ne
         l=ie-1
         do k=1,1
            if (ie.eq.ne1) write(6,116) txt10,k,ie,xmin,xmax,istep,time
            write(6,117) 
            do j=ny2,1,-1
              if (nx2.eq.2) write(6,102) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.3) write(6,103) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.4) write(6,104) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.5) write(6,105) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.6) write(6,106) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.7) write(6,107) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
              if (nx2.eq.8) write(6,108) ((x(i,j,k,e+l),i=1,nx2),e=1,ne)
            enddo
         enddo
      enddo
C
  102 FORMAT(4(2f9.5,2x))
  103 FORMAT(4(3f9.5,2x))
  104 FORMAT(4(4f7.3,2x))
  105 FORMAT(5f9.5,10x,5f9.5)
  106 FORMAT(6f9.5,5x,6f9.5)
  107 FORMAT(7f8.4,5x,7f8.4)
  108 FORMAT(8f8.4,4x,8f8.4)
c
  116 FORMAT(  /,5X,'     ^              ',/,
     $    5X,'   Y |              ',/,
     $    5X,'     |              ',A10,/,
     $    5X,'     +---->         ','Plane = ',I2,'/',I2,2x,2e12.4,/,
     $    5X,'       X            ','Step  =',I9,f15.5)
  117 FORMAT(' ')
c
      if (ichk.eq.1.and.idum.gt.0) call checkit(idum)
      return
      end
c-----------------------------------------------------------------------
      subroutine outmatp(a,m,n,name6,ie)
      include 'SIZE'
      include 'PARALLEL'
      real a(m,n)
      character*6 name6
c
      if (nelv.gt.1) return
      do jid=0,np-1
        imax = iglmax(jid,1)
        if (jid.eq.nid) then
           write(6,*) 
           write(6,*) nid,ie,' matrix: ',name6,m,n
           n12 = min(n,12)
           do i=1,m
              write(6,6) ie,name6,(a(i,j),j=1,n12)
           enddo
    6      format(i3,1x,a6,12f9.5)
           write(6,*) 
           call flush_io()
         endif
        imax = iglmax(jid,1)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine gs_chkr(glo_num)
      include 'SIZE'
      include 'TOTAL'

      integer glo_num(lx1,ly1,lz1,lelt)
      integer e,eg

      do ipass = 1,2
         iquick = 0
         if (nid.eq.0) iquick=glo_num(ipass,1,1,1)
         iquick = iglmax(iquick,1)

         do e=1,nelv
         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
           if (glo_num(i,j,k,e).eq.iquick) then
            eg = lglel(e,node)
            write(6,1) nid,i,j,k,e,eg,iquick,ipass
     $      ,xm1(i,j,k,e),ym1(i,j,k,e),zm1(i,j,k,e)
  1         format(i6,3i4,2i6,i9,i2,1p3e12.4,' iquick')
          endif
         enddo
         enddo
         enddo
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
