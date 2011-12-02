/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */

#ifndef __SUPERLU_SYNCH /* allow multiple inclusions */
#define __SUPERLU_SYNCH

/* Structure for the globally shared work queue */

typedef int qitem_t;

typedef struct {
    int       head, tail, count;
    qitem_t   *queue;
} queue_t;

typedef enum {
    ULOCK,       /* locked once per column */
    LLOCK,       /* locked once per supernode */
    LULOCK,      /* locked once per column in L-supernode */
    NSUPER_LOCK, /* locked once per supernode */
    SCHED_LOCK,  /* locked once per panel, if succeeded each time */
    NO_GLU_LOCKS
} lu_locks_t;

typedef enum {
    RELAXED_SNODE,
    TREE_DOMAIN,   /* domain */
    REGULAR_PANEL  /* non-domain */
} panel_t;

typedef enum {
    DONE,
    BUSY,
    CANGO,
    CANPIPE,
    UNREADY
} pipe_state_t;

typedef struct {
    panel_t      type;  /* panel type: 0 -- relaxed, also domain
			               1 -- domain
			               2 -- regular, non-domain */
    pipe_state_t state; /* one of the 5 states in which the panel can be */
    int          size;  /* in the leading column, the panel size is stored;
	                   in the other columns, the offset (negative)
		           to the leading column is stored */
    int          ukids; /* number of kids not yet finished
			 * In linear pipeline --
			 *   if ukids[firstcol] = 0 then
			 *      the panel becomes a leaf (CANGO)
			 *   if ukids[firstcol] = 1 then
			 *      the panel can be taken as CANPIPE
			 */
} pan_status_t;


/* *********************
   Function prototypes
   *********************/

#ifdef __cplusplus
extern "C" {
#endif

extern void await(volatile int *);

#ifdef __cplusplus
	   }
#endif

#endif /* __SUPERLU_SYNCH */

