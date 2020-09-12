#ifndef GEOMETRY_H
#define GEOMETRY_H

#include"macro.h"
#include"gparam.h"

typedef struct Geometry {
   long **d_nnp;      // d_nnp[r][i] = next neighbour (on the local lattice) in dir.  i of the site r
   long **d_nnm;      // d_nnm[r][i] = next neighbour (on the local lattice) in dir. -i of the site r
   int  **d_bcsitep;  // d_bcsitep[r][i] = -1 if r is on the (positive dir) bounday, else +1
   int  **d_bcsitem;  // d_bcsitem[r][i] = -1 if r is on the (negative dir) bounday, else +1

} Geometry;


// these are the functions to be used in shwitching between different indices
long (*cart_to_si)(int const * const cartcoord, GParam const * const param); // cartesian coordinates -> single index
void (*si_to_cart)(int *cartcoord, long si, GParam const * const param);     // single index -> cartesian coordinates


// general functions
void init_indexing_lex(void); // has to be called before init_geometry
void init_geometry(Geometry *geo, GParam const * const param);
void free_geometry(Geometry *geo, GParam const * const param);


// next neighbour in + direction
inline long nnp(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnp[r][i];
  }


// next neighbour in - direction
inline long nnm(Geometry const * const geo, long r, int i)
  {
  return geo->d_nnm[r][i];
  }


// return -1 if nnp[r][i] is beyond the boundary, else +1
inline int bcsitep(Geometry const * const geo, long r, int i)
  {
  return geo->d_bcsitep[r][i];
  }


// return -1 if nnm[r][i] is beyond the boundary, else +1
inline int bcsitem(Geometry const * const geo, long r, int i)
  {
  return geo->d_bcsitem[r][i];
  }


// for debug
void test_geometry(Geometry const * const geo, GParam const * const param);

//------------ these are not to be used outside geometry.c ----------------

long cart_to_lex(int const * const cartcoord, GParam const * const param);   // cartesian coordinates -> lexicographic index
void lex_to_cart(int *cartcoord, long lex, GParam const * const param);      // lexicographic index -> cartesian coordinates

#endif
