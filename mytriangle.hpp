#ifndef MYTRIANGLE_H
#define MYTRIANGLE_H

// I SEGUNETI DEFINE SONO PER COMPATIBILITA` CON triangle.c

// triangle.c line 213

#ifdef SINGLE
#define REAL_TR float
#else /* not SINGLE */
#define REAL_TR double
#endif /* not SINGLE */

// triangle.c line 315

#ifdef _WIN32
#define TRILIBRARY
#endif // _WIN32

#define VOID int

//-------------------------------------------------------------------

struct triangulateio {
  REAL_TR *pointlist;                                               /* In / out */
  REAL_TR *pointattributelist;                                      /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  REAL_TR *triangleattributelist;                                   /* In / out */
  REAL_TR *trianglearealist;                                         /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  REAL_TR *holelist;                        /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  REAL_TR *regionlist;                      /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  REAL_TR *normlist;                /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

extern "C" {
  void triangulate(char *, struct triangulateio *, struct triangulateio *,
		   struct triangulateio *);
  void trifree(VOID *memptr);
}

#endif
