#include "isosurf.h"

SimplexList *facelist;

int ntcount = 0;
int hangid  = 0;

/******************************************************************************/

void isosurf_init()
{
  facelist = ( SimplexList *) malloc( sizeof(SimplexList));
  facelist->head = NULL;
  facelist->tail = NULL;
}
/******************************************************************************/

int equal_edge( const Edge *e1, const Edge *e2 )
{
  if( (e1->connect[0] == e2->connect[0]) && (e1->connect[1] == e2->connect[1])) return 1;
  if( (e1->connect[0] == e2->connect[1]) && (e1->connect[1] == e2->connect[0])) return 1;
  return 0;
}
/******************************************************************************/
Vertex* getNewVertex()
{
  Vertex *vnew    = (Vertex *)malloc(sizeof(Vertex)); assert(vnew);
  vnew->visit     = 0;
  vnew->edgelist  = NULL;
  vnew->lid       = 0;
  return vnew;
}
/******************************************************************************/

Edge* getNewEdge( Vertex *v1, Vertex *v2)
{
  Edge *newedge = ( Edge*) malloc( sizeof(Edge));
  newedge->connect[0] = v1;
  newedge->connect[1] = v2;
  newedge->hangnode   = NULL;
  newedge->next       = NULL;
  ntcount++;
  return newedge;
}

/******************************************************************************/
Simplex *getNewSimplex()
{
  Simplex *face     = (Simplex *)malloc(sizeof(Simplex)); 
  face->numConnect  = 0;
  face->next        = NULL;
  return face;
}

/******************************************************************************/
Simplex* getNewFace( Vertex *v1, Vertex *v2, Vertex *v3)
{
  Simplex *face  = (Simplex *)malloc(sizeof(Simplex)); 
  face->numConnect = 3;
  face->connect[0] = v1;
  face->connect[1] = v2;
  face->connect[2] = v3;
  face->next       = NULL;
  return face;
}
/******************************************************************************/
void delete_edge( Edge *edge )
{
   if( edge->hangnode ) free(edge->hangnode);
   free( edge );
}

void delete_vertex(Vertex *vertex)
{
  Edge *next, *curr_edge = NULL;
  if(vertex->edgelist) {
     curr_edge = vertex->edgelist->head;
     while(curr_edge) {
           next = curr_edge->next;
           free( curr_edge);
	   curr_edge = next;
     }
     free( vertex->edgelist);
  }
  free(vertex); 
}
/******************************************************************************/
void delete_cell( Simplex *sx )
{
  sx->numConnect = 0;
  free(sx);
}
/******************************************************************************/
void delete_all_objects(Mesh *mesh)
{
  int i;
  for( i = 0; i < mesh->numNodes; i++) 
    delete_vertex( mesh->vertexArray[i] );
  free(mesh->vertexArray);

  for( i = 0; i < mesh->numElems; i++) 
    delete_cell( mesh->cellArray[i] );
  free(mesh->cellArray);
}
/******************************************************************************/

Edge *check_edge( Vertex *va, Vertex *vb)
{
  Edge edge, *curr_edge;
  edge.connect[0] = va;
  edge.connect[1] = vb;

  if( va->edgelist == NULL) return NULL;

  curr_edge = va->edgelist->head;

  while(curr_edge != NULL) {
    if(equal_edge(curr_edge, &edge) ) return curr_edge;
    curr_edge = curr_edge->next;
  }

  return NULL;
}
/******************************************************************************/

double getParam( double AQ, double BQ, double Q)
{
  long double  u = 0.0;
  if( fabs(BQ-AQ) > 1.0E-14) u = (Q-AQ)/(BQ-AQ);

  if( u < 0.0) u = 0.0;
  if( u > 1.0) u = 1.0;
       
  return u;
}
/******************************************************************************/
void setAttribute( Edge *edge, double t)
{
  int k;
  Vertex *hangnode = edge->hangnode;
  Vertex *v1 = edge->connect[0];
  Vertex *v2 = edge->connect[1];
  for( k = 0; k < 3; k++) 
    hangnode->xyzCoords[k] = (1-t)*v1->xyzCoords[k] + t*v2->xyzCoords[k];
  for( k = 0; k < MAX_ATTRIB; k++) 
    hangnode->attrib[k] = (1-t)*v1->attrib[k] + t*v2->attrib[k];
}
/******************************************************************************/
void addEdge( Vertex *vertex, Edge *edge)
{
  if( vertex->edgelist == NULL) {
    vertex->edgelist = ( EdgeList*)malloc(sizeof(EdgeList));
    vertex->edgelist->head = edge;
    vertex->edgelist->tail = edge;
    return;
  }
  vertex->edgelist->tail->next = edge;
  vertex->edgelist->tail = edge;
}
/******************************************************************************/

Vertex* VertexInterp( Vertex *nA,  Vertex *nB, int attribID, double isoVal)
{
  Vertex *vertex_on_edge, *nmin, *nmax;
  Edge *oldedge;
  Edge *newedge;
  double paramval;
  double aval, bval;

  nmin = nA < nB ? nA : nB;
  nmax = nA > nB ? nA : nB;

  oldedge = check_edge(nmin,nmax);

  if( oldedge != NULL ) return oldedge->hangnode;

  if( nmin == nA ) {
    aval =  nA->attrib[attribID]; bval =  nB->attrib[attribID];
  } else {
    aval =  nB->attrib[attribID]; bval =  nA->attrib[attribID];
  }

  paramval = getParam(aval, bval, isoVal);

  newedge  = getNewEdge(nmin, nmax);
  vertex_on_edge = getNewVertex();
  vertex_on_edge->lid = hangid++;
  newedge->hangnode = vertex_on_edge;
  setAttribute(newedge, paramval);
  addEdge(nmin,newedge);

  return vertex_on_edge;
}
/******************************************************************************/
void addNewTriangle(Vertex *v1, Vertex *v2, Vertex *v3)
{

  Simplex *face = getNewFace(v1,v2,v3);

  assert( (v1->lid != v2->lid) && (v1->lid != v3->lid ) && (v2->lid != v3->lid));

  if( facelist->head == NULL ) {
    facelist->head = face;
    facelist->tail = face;
    return;
  }
  facelist->tail->next = face;
  facelist->tail = face;
}
/******************************************************************************/
int TetLookupTable(Simplex *tet, int eindex, int attribID, double isoVal)
{

  Vertex *n0 = tet->connect[0];
  Vertex *n1 = tet->connect[1];
  Vertex *n2 = tet->connect[2];
  Vertex *n3 = tet->connect[3];

  Vertex *tn0, *tn1, *tn2;
  
  int ncount = 0;
  switch (eindex) 
    {
    case 0x00:
    case 0x0F:
      break;
    case 0x0E:
    case 0x01:
      tn0  = VertexInterp(n0, n1, attribID, isoVal);
      tn1  = VertexInterp(n0, n2, attribID, isoVal);
      tn2  = VertexInterp(n0, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2);  
      ncount = 1;
      break;
    case 0x0D:
    case 0x02:
      tn0 = VertexInterp( n1, n0, attribID, isoVal);
      tn1 = VertexInterp( n1, n3, attribID, isoVal);
      tn2 = VertexInterp( n1, n2, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 
      ncount = 1;
      break;
    case 0x0C:
    case 0x03:
      tn0 = VertexInterp( n0, n3, attribID, isoVal);
      tn1 = VertexInterp( n0, n2, attribID, isoVal);
      tn2 = VertexInterp( n1, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 

      tn0 = VertexInterp( n1, n3, attribID, isoVal);
      tn1 = VertexInterp( n1, n2, attribID, isoVal);
      tn2 = VertexInterp( n0, n2, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 
      ncount = 2;
      break;
    case 0x0B:
    case 0x04:
      tn0 = VertexInterp( n2, n0, attribID, isoVal);
      tn1 = VertexInterp( n2, n1, attribID, isoVal);
      tn2 = VertexInterp( n2, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 
      ncount = 1;
      break;
    case 0x0A:
    case 0x05:
      tn0  = VertexInterp( n0, n1, attribID, isoVal);
      tn1  = VertexInterp( n2, n3, attribID, isoVal);
      tn2  = VertexInterp( n0, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 

      tn0  = VertexInterp( n0, n1, attribID, isoVal);
      tn1  = VertexInterp( n1, n2, attribID, isoVal);
      tn2  = VertexInterp( n2, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 
      ncount = 2;
      break;
    case 0x09:
    case 0x06:
      tn0  = VertexInterp( n0, n1, attribID, isoVal);
      tn1  = VertexInterp( n1, n3, attribID, isoVal);
      tn2  = VertexInterp( n2, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2 ); 

      tn0   = VertexInterp( n0, n1, attribID, isoVal);
      tn1   = VertexInterp( n0, n2, attribID, isoVal);
      tn2   = VertexInterp( n2, n3, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2 ); 
      ncount = 2;
      break;
    case 0x07:
    case 0x08:
      tn0   = VertexInterp( n3, n0, attribID, isoVal);
      tn1   = VertexInterp( n3, n2, attribID, isoVal);
      tn2   = VertexInterp( n3, n1, attribID, isoVal);
      addNewTriangle(tn0, tn1, tn2); 
      ncount = 1;
      break;
    }

  return ncount;
}
/******************************************************************************/
int isActiveTet(Simplex *tet, int attribID, double isoVal)
{

  int eindex = 0;
  if (tet->connect[0]->attrib[attribID] < isoVal) eindex |= 1;
  if (tet->connect[1]->attrib[attribID] < isoVal) eindex |= 2;
  if (tet->connect[2]->attrib[attribID] < isoVal) eindex |= 4;
  if (tet->connect[3]->attrib[attribID] < isoVal) eindex |= 8;

  return eindex;
}
/******************************************************************************/

int SearchTet(Simplex *tet, int attribID, double isoVal)
{
  int eindex = isActiveTet( tet, attribID, isoVal);

  return TetLookupTable( tet, eindex, attribID, isoVal);
}
/******************************************************************************/

void buildTet( Simplex *tet, Vertex *v0, Vertex *v1, Vertex *v2, Vertex *v3)
{
  tet->numConnect = 4;
  tet->connect[0] = v0;
  tet->connect[1] = v1;
  tet->connect[2] = v2;
  tet->connect[3] = v3;
}
/******************************************************************************/
void DecomposeHex( Simplex *hex, Simplex *tets, int *ntets)
{
  int i;
  Vertex *connect[8];
  
  for( i = 0; i < 8; i++) connect[i] = hex->connect[i];

  buildTet( &tets[0], connect[0], connect[5], connect[4], connect[6] );
  buildTet( &tets[1], connect[0], connect[7], connect[5], connect[6] );
  buildTet( &tets[2], connect[0], connect[2], connect[7], connect[6] );
  buildTet( &tets[3], connect[0], connect[7], connect[1], connect[5] );
  buildTet( &tets[4], connect[0], connect[1], connect[7], connect[2] );
  buildTet( &tets[5], connect[3], connect[7], connect[1], connect[2] );

  *ntets = 6;
}
/******************************************************************************/
int isActiveHex( Simplex *hex, int attribID, double isoVal)
{
  int i, ntets, stat = 0;
  Simplex  tets[6];

  if( hex  == NULL ) return 0;
  DecomposeHex(hex, tets, &ntets);

  for( i = 0; i < ntets; i++) 
    if(isActiveTet(&tets[i], attribID, isoVal)) stat = 1;

  return 0;
}
/******************************************************************************/
Mesh *convert_to_tets( Mesh *hexmesh)
{
  int  i, j, numHexs, numNodes;
  Simplex *tets[6];
  Vertex *connect[8];

  Mesh *tetmesh = (Mesh *)malloc(sizeof(Mesh));

  numHexs = hexmesh->numElems;
  tetmesh->cellArray = ( Simplex**)malloc(6*numHexs*sizeof(Simplex*));

  for(i = 0; i < hexmesh->numElems; i++) {
    for( j = 0; j < 8; j++) connect[j] = hexmesh->cellArray[i]->connect[j];
    for( j = 0; j < 6; j++) tets[j] = getNewSimplex();
    buildTet( tets[0], connect[0], connect[5], connect[4], connect[6] );
    buildTet( tets[1], connect[0], connect[7], connect[5], connect[6] );
    buildTet( tets[2], connect[0], connect[2], connect[7], connect[6] );
    buildTet( tets[3], connect[0], connect[7], connect[1], connect[5] );
    buildTet( tets[4], connect[0], connect[1], connect[7], connect[2] );
    buildTet( tets[5], connect[3], connect[7], connect[1], connect[2] );
    for( j = 0; j < 6; j++) tetmesh->cellArray[6*i+j] = tets[j];
  }

  numNodes = hexmesh->numNodes;
  tetmesh->vertexArray = ( Vertex **)malloc(numNodes*sizeof(Vertex*));

  for( i = 0; i < numNodes; i++) 
    tetmesh->vertexArray[i] = hexmesh->vertexArray[i];

  tetmesh->numElems = 6*numHexs;
  tetmesh->numNodes = numNodes;

  return tetmesh;
}
/*#############################################################################*/
void SearchHex( Simplex *cell, int attribID, double isoVal)
{
  int i, ntets;
  Simplex  tets[6];
  if( cell == NULL ) return;
  DecomposeHex(cell, tets, &ntets);

  for( i = 0; i < ntets; i++) SearchTet(&tets[i], attribID, isoVal);
}
/******************************************************************************/
Vertex *binarySearch(Vertex **vertexArray, int numNodes, int key)
{
  int vid, low, high, mid;

  low = 0;
  high = numNodes;

  while (low < high) {
    mid = (low + high) / 2;
    vid = vertexArray[mid]->lid;
    if (key < vid)
      high = mid;
    else if (key > vid)
      low = mid + 1;
    else return vertexArray[mid];
  }
  return NULL;
}

/******************************************************************************/
void collect_trimesh( int *tConnect, float *txyzCoords, float *tAttribs, 
                      int *numTriElems, int *numTriNodes, int numAttribs)
{
  int i, j, k, numFaces, numNodes, vid, n0, n1, n2;
  Vertex  *vertex;
  Simplex *curr_face;

  curr_face = facelist->head;

  if( curr_face == NULL) return;

  numFaces =  0;
  while( curr_face ) {
    curr_face->connect[0]->visit = 0;
    curr_face->connect[1]->visit = 0;
    curr_face->connect[2]->visit = 0;
    numFaces++;
    curr_face = curr_face->next;
  }

  numNodes = 0;
  if( numFaces ) {
    curr_face = facelist->head;
    for( i = 0; i < numFaces; i++) {
        for(j = 0; j < 3; j++) {
            vertex = curr_face->connect[j];
	    if( vertex->visit == 0) {
	        numNodes++; vertex->visit = 1;
	    }
         }
         curr_face = curr_face->next;
     }
    *numTriElems = numFaces;
    *numTriNodes = numNodes;
  }

  if( numNodes ) {
    curr_face = facelist->head;
    vid = 0;
    for( i = 0; i < numFaces; i++) {
      for( j = 0; j < 3; j++) {
	vertex = curr_face->connect[j];
	if( vertex->visit == 1) {
	  vertex->visit = 0;
	  vertex->lid   = vid;
	  txyzCoords[3*vid+0] = vertex->xyzCoords[0];
	  txyzCoords[3*vid+1] = vertex->xyzCoords[1];
	  txyzCoords[3*vid+2] = vertex->xyzCoords[2];
	  for( k = 0;  k < numAttribs; k++) 
	    tAttribs[vid*numAttribs+k] = vertex->attrib[k];
	  vid++;
	}
      }
      n0 = curr_face->connect[0]->lid;
      n1 = curr_face->connect[1]->lid;
      n2 = curr_face->connect[2]->lid;
      assert( (n0 != n1) && (n0 != n2) && (n1 != n2));
      tConnect[3*i+0] = n0;
      tConnect[3*i+1] = n1;
      tConnect[3*i+2] = n2;
      curr_face = curr_face->next;
    }
  }

}

/******************************************************************************/
void  memestimate_( const int *hConnect, const double *xyzCoords, double *hAttrib, 
		    int *vidArray, int *nElems, int *nNodes, int *nAttribs, 
		    int *attID, double *iVal, int *numTriElems, int *numTriNodes)
{
  int i, j, numActive;

  Simplex cell;
  Vertex  *vertex;

  int numNodes   = *nNodes;
  int numElems   = *nElems;
  int numAttribs = *nAttribs;
  int attribID   = *attID;
  double isoVal  = *iVal;

  Vertex **vertexArray = (Vertex **)malloc(numNodes*sizeof(Vertex *) );
  for( i = 0; i < numNodes; i++) {
    vertex = getNewVertex();
    vertex->lid = vidArray[i];
    for(j = 0; j < numAttribs; j++)
      vertex->attrib[j] = hAttrib[numAttribs*i+j];
    vertexArray[i] = vertex;
  }
  numActive = 0;
  for( i = 0; i < numElems; i++) {
    cell.numConnect = 8;
    for( j = 0; j < 8; j++) {
      int nid  = hConnect[8*i+j];
      vertex   = binarySearch(vertexArray, numNodes, nid); assert(vertex);
      cell.connect[j] = vertex;
    }
    if( isActiveHex(&cell, attribID, isoVal) ) numActive++;
  }
  *numTriElems = 4*numActive;
  *numTriNodes = 3*numActive;

  for( i = 0; i < numNodes; i++) delete_vertex(vertexArray[i]);
  free(vertexArray);
}

/******************************************************************************/
#ifdef UPCASE
void  ISOSURFACE ( const int *hConnect, const double *xyzCoords, int *vidArray, 
                   double *hAttrib, int *nElems, int *nNodes, int *nAttribs, 
                   int *attID, double *iVal, int *tConnect,  float *txyzCoords, 
		   float *tAttrib, int *numTriElems, int *numTriNodes)
#elif UNDERSCORE
void  isosurface_( const int *hConnect, const double *xyzCoords, int *vidArray, 
                   double *hAttrib, int *nElems, int *nNodes, int *nAttribs, 
                   int *attID, double *iVal, int *tConnect,  float *txyzCoords, 
		   float *tAttrib, int *numTriElems, int *numTriNodes)
#else
void  isosurface ( const int *hConnect, const double *xyzCoords, int *vidArray, 
                   double *hAttrib, int *nElems, int *nNodes, int *nAttribs, 
                   int *attID, double *iVal, int *tConnect,  float *txyzCoords, 
		   float *tAttrib, int *numTriElems, int *numTriNodes)
#endif
{

  /*****************************************************************************
   * Input Arguments:
   * hConnect  :  Array of hex elements connectivity. n11,n12,   n18, n21, n22.....
   * xyzCoords :  Array of vertex coordinates of the(3D). x11
   * vidArray  :  Array of vertex IDs.
   * hAttrib   :  Array of attribs of vertices
   * nElems    :  Number of elements in the hConnect array.
   * nNodes    :  Number of elements in the xyzCoord, vidArray, hAttribs
   * attID     :  ID of the attrib for which isosurface is sought.
   * iVal      :  Scalar value for which isosurface is sought.
   *
   * Output Arguments:
   * tConnect    : Array of triangular elements connectivity.
   * txyzCoords  : Arrary of triangle nodes coordinates.(3D)
   * tAttribs    : Number of attributes on the nodes ( Same as Hex nodes ).
   * numTriElems : Number of triangular elements in the isosurface.
   * numTriNodes : Number of triangular nodes  in the isosurface.
   *******************************************************************************/
  int      i, j, numNodes, numElems, numAttribs, attribID;
  double   x, y, z, isoVal;
  Vertex   *vertex;
  Vertex  **vertexArray;
  Simplex **hexArray, *cell;

  isosurf_init(); 

  numNodes   = *nNodes;
  numElems   = *nElems;
  attribID   = *attID;
  isoVal     = *iVal;
  numAttribs = *nAttribs;

  vertexArray = (Vertex **)malloc(numNodes*sizeof(Vertex *) );
  hexArray    = (Simplex **)malloc(numElems*sizeof(Simplex *));

  for( i = 0; i < numNodes; i++) {
    Vertex *vertex = getNewVertex();
    vertex->lid  =  vidArray[i];
    x = xyzCoords[3*i+0];
    y = xyzCoords[3*i+1];
    z = xyzCoords[3*i+2];
    vertex->xyzCoords[0] = x;
    vertex->xyzCoords[1] = y;
    vertex->xyzCoords[2] = z;
    for( j = 0; j < numAttribs; j++)
      vertex->attrib[j] =   hAttrib[i*numAttribs+j];
    vertexArray[i] =  vertex;
  }

  for( i = 0; i < numElems; i++) {
    cell =  (Simplex*) malloc(sizeof( Simplex));
    cell->numConnect = 8;
    for( j = 0; j < 8; j++) {
      int nid    = hConnect[8*i+j];
      vertex     = binarySearch(vertexArray, numNodes, nid); assert( vertex );
      cell->connect[j] =  vertex;
    }
    hexArray[i]  = cell;
  }

  for( i = 0; i < numElems; i++) SearchHex(hexArray[i], attribID, isoVal);

  collect_trimesh(tConnect, txyzCoords, tAttrib, numTriElems, numTriNodes, numAttribs);

  for( i = 0; i < numNodes; i++) delete_vertex(vertexArray[i]);
  free(vertexArray);

  for( i = 0; i < numElems; i++) delete_cell( hexArray[i]);
  free(hexArray);
}
/******************************************************************************/
