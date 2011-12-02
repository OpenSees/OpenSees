/*
 * struct.h
 *
 * This file contains data structures for the multilevel ordering schemes
 *
 * Started 8/27/94
 * George
 *
 * $Id: struct.h,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edgedef {
  int edge;
  int ewgt;
};
typedef struct edgedef EdgeType;

/*************************************************************************
* The following data structure stores information about a graph vertex
**************************************************************************/
struct vtxdef {
 int vwgt;	/* The weight of the vertex */
 int cewgt; 	/* The edge weight being collapsed in this node */
 int nedges;	/* The number of edges. Note this is different than the weight */
 int ewgtsum;   /* The sum of the edge weights */
 EdgeType *edges;	/* List of edges */
};

typedef struct vtxdef VertexType;



/*************************************************************************
* The following data structure implements a link list of matching pairs
**************************************************************************/
struct matchdef {
  int v, u;		 /* Nodes being matched. Unmatched nodes v == u */
  struct matchdef *next; /* Next node in the match list */
};

typedef struct matchdef VertexPairType;


/*************************************************************************
* The following data structure implements a closed hash-table
**************************************************************************/
struct htabledef {
  int size;		/* The size of the hash table */
  int nelem;		/* The # of elements in the hash table */
  int *ht;		/* The hash table itself */
};

typedef struct htabledef HTableType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 int id, ed;            /* ID/ED of edges */
 int ndegrees;          /* The number of different ext-degrees */
 EdgeType *degrees;     /* List of edges */
};

typedef struct rinfodef RInfoType;



/*************************************************************************
* The following data structure is a node in a doubly-link node of graphs
* used to represent the coarsion sequence.
* Each node in this list will contain the current graph, and mapping
* information to obtain the coarser graph.
**************************************************************************/
struct crsgraphdef {
  int nvtxs;			/* The number of vertices */
  int nedges;			/* The total number of edges */
  int tvwgt;			/* The total vertex weight */
  VertexType *allvtxs;		/* The pool of vertices from where you allocate */
  VertexType **vtxs;		/* The list of vertices */
  int *cmap;			/* The coersion map for the graph */
  int *match;			/* The match array */
  struct crsgraphdef *coarser;	/* Pointer to the coarser graph */
  struct crsgraphdef *finer;	/* Pointer to the finer graph */
  int level;			/* The level of the graph in the coarsening sequence */

  /*
   * Refinement related data structures.
   * They are here because they can be used to save some 
   * computation
   */
  int *where;		/* In which partition the vertices belong */
  int *id;		/* The internal degree of a vertex */
  int *ed;		/* The external degree of a vertex */
  int pwgts[2];		/* The amount of vertex weight in the 2 partitions */
  int mincut;		/* The size of the minimum cut */

  /*
   * Boundary refinement related data structures
   */
  int nbnd;		/* Number of vertices at the boundary */
  HTableType htable;	/* The hash table of boundary vertices */

  /*
   * Arrays used in ordering 
   */
  int *label;		/* The label array, with respect to original matrix */ 

  /* 
   * K-Way refinement additional data structures
   */
  RInfoType *rinfo;
  int *kpwgts;

};

typedef struct crsgraphdef CoarseGraphType;



/*************************************************************************
* The following data structure will hold a node of the link list.
* The matrix i stored as 2 link lists, one is row wise and the other 
* is column wise.
* At this point the data structure doesn't perform numerical computations.
**************************************************************************/
struct ListNodeType {
  int pos;                /* The position along row/column */
  float val;              /* The value of the entry in the matrix */
  struct ListNodeType *link;   /* It's a link list */
};

typedef struct ListNodeType ListNodeType;


/*************************************************************************
* The following data structure is a node from a doubly link list of 
* gain lists. This is an alternate way of implementing the BucketListType
**************************************************************************/
struct GainBucketType {
  int gain;			/* The gain of this bucket */
  ListNodeType *gainlist;	/* The list of nodes at this gain */
  struct GainBucketType 
     *pbucket, 			/* The previous bucket */
     *nbucket;			/* The next bucket */
};

typedef struct GainBucketType GainBucketType;

/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct BucketListType {
  int type;			/* The type of the representation used */
  int nnodes;
  int maxnodes;
  int pgainspan, ngainspan;	/* plus and negative gain span */
  int mingain, maxgain;
  ListNodeType **buckets;
  GainBucketType *head, *tail;
  GainBucketType *core;		/* Core memory for GainBucketType */
};

typedef struct BucketListType BucketListType;



/*************************************************************************
* The following data structure is used to sort the values of evector
**************************************************************************/
struct EvecWgtType {
  double ei;
  int wgt;
  int index;
};
typedef struct EvecWgtType EvecWgtType;


/*************************************************************************
* The following data structure stores a node of the separator tree
**************************************************************************/
struct SepNodeType {
  int nvtxs;	/* The numbre of vertices in the graph */
  int li;	/* Low index of separator (inclusive) */
  int hi;	/* High index of separator (exclusive) */
  int isleaf;	/* Indicates if it is a leaf node */
  double opc;	/* The opcount for this node */
  double subopc;	/* The opcount for the subtree of this node */
};

typedef struct SepNodeType SepNodeType;


/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  int key;
  int val;
};

typedef struct KeyValueType KeyValueType;



/*************************************************************************
* The following data structure stores the various memory pools
* and parameters of the multilevel algorithm
**************************************************************************/
struct controldef {
  /* Edge Pool */
  EdgeType *edgepool;   /* The pool of edges */
  int lastedge;         /* Pointer to the first unused edge */
  int maxedges;         /* Number of edges in the pool */

  /* External Degree Pool */
  EdgeType *degrees;    /* The pool of edges */
  int lastdegree;       /* Pointer to the first unused edge */
  int maxdegrees;       /* Number of edges in the pool */

  /* Icore Pool */
  int *icore;		/* Core memory for integer arrays */
  int maxicore;		/* Maximum number of ints */
  int cicore;		/* The current free icore */

  /* Gain Pool */
  ListNodeType *gaincore;  /* Core memory for bucket gains */
  int maxgain;             /* The maximum number of gains allocated */
  int cgain;               /* The current free gain */

  /* Gain Bucket Pool */
  GainBucketType *gbcore;
  int maxbucket;
  int cbucket;

  /* Various control information */
  int dbglvl;
  int CoarsenTo;
  int MatchType;
  int InitPartType;
  int RefineType;
  int OpType;
  int IsWeighted;
  int nparts;
  float cfrac;
};
typedef struct controldef CtrlType;


