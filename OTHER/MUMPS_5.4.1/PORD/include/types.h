/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: types.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains the fundamental data structures
/
******************************************************************************/

/*****************************************************************************
A macro defining the size of integers 
(modified for compatibility with MUMPS)
******************************************************************************/
#if defined(INTSIZE64) || defined(PORD_INTSIZE64)
#include <inttypes.h>
#define PORD_INT int64_t
#else
#define PORD_INT int
#endif

typedef double FLOAT;
typedef PORD_INT    options_t;
typedef FLOAT  timings_t;

/*****************************************************************************
Graph object
******************************************************************************/
typedef struct _graph {
  PORD_INT    nvtx;
  PORD_INT    nedges;
  PORD_INT    type;
  PORD_INT    totvwght;
  PORD_INT   *xadj;
  PORD_INT   *adjncy;
  PORD_INT   *vwght;
} graph_t;

/*****************************************************************************
Graph bisection object
******************************************************************************/
typedef struct _gbisect {
  graph_t      *G;
  PORD_INT     *color;
  PORD_INT      cwght[3];
} gbisect_t;

/*****************************************************************************
Domain decomposition object
******************************************************************************/
typedef struct _domdec {
  graph_t      *G;
  PORD_INT      ndom;
  PORD_INT      domwght;
  PORD_INT     *vtype;
  PORD_INT     *color;
  PORD_INT      cwght[3];
  PORD_INT     *map;
  struct _domdec *prev, *next;
} domdec_t;

/*****************************************************************************
Bipartite graph object
******************************************************************************/
typedef struct _gbipart {
  graph_t *G;
  PORD_INT     nX;
  PORD_INT     nY;
} gbipart_t;

/*****************************************************************************
Recursive nested dissection object
******************************************************************************/
typedef struct _nestdiss {
  graph_t      *G;
  PORD_INT     *map;
  PORD_INT      depth;
  PORD_INT      nvint;
  PORD_INT     *intvertex;
  PORD_INT     *intcolor;
  PORD_INT      cwght[3];
  struct _nestdiss *parent, *childB, *childW;
} nestdiss_t;

/*****************************************************************************
Multisector object
******************************************************************************/
typedef struct _multisector {
  graph_t      *G;
  PORD_INT     *stage;
  PORD_INT      nstages;
  PORD_INT      nnodes;
  PORD_INT      totmswght;
} multisector_t;

/*****************************************************************************
Elimination graph object
******************************************************************************/
typedef struct _gelim {
  graph_t      *G;
  PORD_INT      maxedges;
  PORD_INT     *len;
  PORD_INT     *elen;
  PORD_INT     *parent;
  PORD_INT     *degree;
  PORD_INT     *score;
} gelim_t;

/*****************************************************************************
Bucket structure object
******************************************************************************/
typedef struct _bucket {
  PORD_INT    maxbin, maxitem;
  PORD_INT    offset;
  PORD_INT    nobj;
  PORD_INT    minbin;
  PORD_INT   *bin;
  PORD_INT   *next;
  PORD_INT   *last;
  PORD_INT   *key;
} bucket_t;

/*****************************************************************************
Minimum priority object
******************************************************************************/
typedef struct _stageinfo stageinfo_t;
typedef struct _minprior {
  gelim_t            *Gelim;
  multisector_t      *ms;
  bucket_t           *bucket;
  stageinfo_t        *stageinfo;
  PORD_INT           *reachset;
  PORD_INT            nreach;
  PORD_INT           *auxaux;
  PORD_INT           *auxbin;
  PORD_INT           *auxtmp;
  PORD_INT            flag;
} minprior_t;
struct _stageinfo {
  PORD_INT   nstep;
  PORD_INT   welim;
  PORD_INT   nzf;
  FLOAT ops;
};

/*****************************************************************************
Elimination tree object
******************************************************************************/
typedef struct _elimtree {
  PORD_INT    nvtx;
  PORD_INT    nfronts;
  PORD_INT    root;
  PORD_INT   *ncolfactor;
  PORD_INT   *ncolupdate;
  PORD_INT   *parent;
  PORD_INT   *firstchild;
  PORD_INT   *silbings;
  PORD_INT   *vtx2front;
} elimtree_t;

/*****************************************************************************
Input matrix object
******************************************************************************/
typedef struct _inputMtx {
  PORD_INT    neqs;
  PORD_INT    nelem;
  FLOAT      *diag;
  FLOAT      *nza;
  PORD_INT   *xnza;
  PORD_INT   *nzasub;
} inputMtx_t;

/*****************************************************************************
Dense matrix object
******************************************************************************/
typedef struct _workspace workspace_t;
typedef struct _denseMtx {
  workspace_t      *ws;
  PORD_INT          front;
  PORD_INT          owned;
  PORD_INT          ncol;
  PORD_INT          nrow;
  PORD_INT          nelem;
  PORD_INT          nfloats;
  PORD_INT         *colind;
  PORD_INT         *rowind;
  PORD_INT         *collen;
  FLOAT            *entries;
  FLOAT            *mem;
  struct _denseMtx *prevMtx, *nextMtx;
} denseMtx_t;
struct _workspace {
  FLOAT          *mem;
  PORD_INT        size;
  PORD_INT        maxsize;
  PORD_INT        incr;
  denseMtx_t     *lastMtx;
};

/*****************************************************************************
Compressed subscript structure object
******************************************************************************/
typedef struct _css {
  PORD_INT    neqs;
  PORD_INT    nind;
  PORD_INT    owned;
  PORD_INT   *xnzl;
  PORD_INT   *nzlsub;
  PORD_INT   *xnzlsub;
} css_t;

/*****************************************************************************
Front subscript object
******************************************************************************/
typedef struct _frontsub {
  elimtree_t      *PTP;
  PORD_INT         nind;
  PORD_INT        *xnzf;
  PORD_INT        *nzfsub;
} frontsub_t;

/*****************************************************************************
Factor matrix object
******************************************************************************/
typedef struct _factorMtx {
  PORD_INT         nelem;
  PORD_INT        *perm;
  FLOAT           *nzl;
  css_t           *css;
  frontsub_t      *frontsub;
} factorMtx_t;

/*****************************************************************************
Mapping object
******************************************************************************/
typedef struct _groupinfo groupinfo_t;
typedef struct {
  elimtree_t       *T;
  PORD_INT          dimQ;
  PORD_INT          maxgroup;
  PORD_INT         *front2group;
  groupinfo_t      *groupinfo;
} mapping_t;
struct _groupinfo {
  FLOAT      ops;
  PORD_INT   nprocs;
  PORD_INT   nfronts;
};

/*****************************************************************************
Topology object
******************************************************************************/
typedef struct {
  PORD_INT      nprocs;
  PORD_INT      mygridId;
  PORD_INT      dimX;
  PORD_INT      dimY;
  PORD_INT      myQId;
  PORD_INT      dimQ;
  PORD_INT      *cube2grid;
#ifdef PARIX
  LinkCB_t      **link;
#endif
#ifdef MPI
  MPI_Comm      comm;
  MPI_Status    status;
#endif
} topology_t;

/*****************************************************************************
Communication buffer object
******************************************************************************/
typedef struct {
  char   *data;
  size_t len;
  size_t maxlen;
} buffer_t;

/*****************************************************************************
Bit mask object
******************************************************************************/
typedef struct {
  PORD_INT   dimQ;
  PORD_INT   maxgroup;
  PORD_INT   mygroupId;
  PORD_INT   offset;
  PORD_INT   *group;
  PORD_INT   *colbits, *colmask;
  PORD_INT   *rowbits, *rowmask;
} mask_t;


