/*
 * defs.h
 *
 * This file contains constant definitions
 *
 * Started 8/27/94
 * George
 *
 * $Id: defs.h,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */


/*
 * Default values for ometis, pmetis, and kmetis
 */
#define OMETIS_CTO	100
#define OMETIS_MTYPE	MATCH_SHEM
#define OMETIS_IPTYPE	INITPART_GGPKL
#define OMETIS_RTYPE 	REFINE_BGKLR

#define PMETIS_CTO	100
#define PMETIS_MTYPE	MATCH_SHEM
#define PMETIS_IPTYPE	INITPART_GGPKL
#define PMETIS_RTYPE 	REFINE_BGKLR

#define KMETIS_CTO	2000
#define KMETIS_MTYPE	MATCH_SHEM
#define KMETIS_IPTYPE	INITPART_GGPKL
#define KMETIS_RTYPE 	REFINE_BGR


/*
 * Various constants
 */
#define MMDSWITCH	200

#define KWAY_REF_GREEDY_NITER   8
#define KWAY_REF_FM_NITER       6

#define NPARTS_FACTOR           12
#define UNBALANCE_FRACTION      1.03
#define UNBALANCE_FACTOR        1

#define ET_METIS	1
#define ET_METIS_P	2
#define ET_METIS_P_K	3
#define ET_METIS_O	4

#define PLUS_GAINSPAN	2000		/* Parameters for FM buckets */
#define NEG_GAINSPAN	500
#define MAXDEGREE	500

#define HT_EMPTY	-1		/* An empty htable entry */

#define STUPID		0		/* Type of boundary refinement */
#define SMART		1		/* Type of boundary refinement */

#define UNMATCHED	-1  	 	/* Denotes that a vertex is unmatched */
#define MAXLINE		128*1024	/* Line length during I/O */	
#define MAXIDX		((1<<30) -1)	/* Maximum value for an index */

#define INC_ORDER	1
#define DEC_ORDER	2

#define MAXLANCZOS	100		/* Maximum number of Lanczos iterations */
#define LANCZOSEPS	1e-2		/* Lanczos acuracy */

#define COARSEN_FRACTION	0.90	/* Node reduction between succesive coarsening levels */
#define LOCAL_COARSEN_FRACTION	0.85	/* Node reduction between succesive coarsening levels */
#define KWAY_COARSEN_FRACTION   0.90 

#define GGP_NTRIALS	10
#define GGGP_NTRIALS	6
#define GGPKL_NTRIALS	4


/* Match Type */
#define MATCH_RM			1
#define MATCH_HEM			2
#define MATCH_LEM			3
#define MATCH_HCM			4
#define MATCH_MHEM			5
#define MATCH_SRM			11
#define MATCH_SHEM			21
#define MATCH_SMHEM			51


/* Initial Partition Types */
#define INITPART_GGP			1
#define INITPART_GGGP			2
#define INITPART_EIG			3
#define INITPART_GGPKL			4


/* Refinement Algorithm Types */
#define REFINE_GR		1
#define REFINE_KLR		2	
#define REFINE_GKLR		3
#define REFINE_BGR 		11
#define REFINE_BKLR 		12
#define REFINE_BGKLR		13
#define REFINE_NONE		20


/* Operation Type */
#define OP_RMLB		1
#define OP_MLKP		2
#define OP_MLND		3


/* Debug Levels */
#define DBG_OUTPUT	1		/* Outputs partition/permutation vectors */
#define DBG_TIME	2		/* Print timing statistics */
#define DBG_PROGRESS   	4		/* Show the coersion progress */
#define DBG_INITCUT	8		/* Show the initial cuts */
#define DBG_FFCUT	16		/* Show First and Final cuts */
#define DBG_PARTSIZES	32		/* Shows the sizes of the final parts */
#define DBG_ORDERNODE	64		/* Shows the nodes being touched by the cut */
#define DBG_ORDERBAL	128		/* Shows the balance of nodes */
#define DBG_KWAYREF	256		/* Write the graph out */

#define DBG_ITERCUT	512		/* Show cut progress as you go */
#define DBG_CRSSTAT	1024		/* Show stat during coarsening */
#define DBG_GRAPHOUT	2048		/* Write the graph out */

#define DBG_BALCUT	4096		/* Info during the balancing phase */

