/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SuperLU_MT_pxgstrf_synch.h,v $
                                                                        
                                                                        
/*
 * -- SuperLU MT routine (alpha version) --
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

#endif /* __SUPERLU_SYNCH */

