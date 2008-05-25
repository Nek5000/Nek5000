#ifndef ISO_H
#define ISO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define MAX_ATTRIB  5
#define VTK_LINE 3
#define VTK_TRIANGLE 5
#define VTK_QUAD 9
#define VTK_TET  10
#define VTK_HEX  12

typedef struct GVertex
{
  int  lid;
  short  visit;
  double xyzCoords[3];
  double attrib[MAX_ATTRIB];
  struct ListEdge *edgelist;
} Vertex;

typedef struct GEdge
{
  Vertex* connect[2];
  Vertex* hangnode;
  struct GEdge *next;
} Edge;

typedef struct ListEdge
{
  Edge *head, *tail;
} EdgeList;

typedef struct GSimplex
{
  Vertex  *connect[8];
  short   numConnect;
  struct  GSimplex *next;
} Simplex;

typedef struct ListSimplex
{
  Simplex *head, *tail;
} SimplexList;

typedef struct
{
    Vertex  **vertexArray;
    Simplex **cellArray;
    size_t  numNodes, numElems;
} Mesh;

Vertex*  getNewVertex();
Simplex* getNewSimplex();


#endif

