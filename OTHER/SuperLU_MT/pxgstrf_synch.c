/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <stdio.h>
#include <math.h>
#include "pdsp_defs.h"
#include "util.h"

#define SPLIT_TOP

int
ParallelInit(int n, pxgstrf_relax_t *pxgstrf_relax, 
	     pdgstrf_options_t *pdgstrf_options, 
	     pxgstrf_shared_t *pxgstrf_shared)
{
    int      *etree = pdgstrf_options->etree;
    register int w, dad, ukids, i, j, k, rs, panel_size, relax;
    register int P, w_top, do_split = 0;
    panel_t panel_type;
    int      *panel_histo = pxgstrf_shared->Gstat->panel_histo;
    register int nthr, concurrency, info;

#if ( MACH==SUN )
    register int sync_type = USYNC_THREAD;
    
    /* Set concurrency level. */
    nthr = sysconf(_SC_NPROCESSORS_ONLN);
    thr_setconcurrency(nthr);            /* number of LWPs */
    concurrency = thr_getconcurrency();

#if ( PRNTlevel==1 )    
    printf(".. CPUs %d, concurrency (#LWP) %d, P %d\n",
	   nthr, concurrency, P);
#endif

    /* Initialize mutex variables. */
    pxgstrf_shared->lu_locks = (mutex_t *) 
        SUPERLU_MALLOC(NO_GLU_LOCKS * sizeof(mutex_t));
    for (i = 0; i < NO_GLU_LOCKS; ++i)
	mutex_init(&pxgstrf_shared->lu_locks[i], sync_type, 0);

#elif ( MACH==DEC || MACH==PTHREAD )
    pxgstrf_shared->lu_locks = (pthread_mutex_t *) 
        SUPERLU_MALLOC(NO_GLU_LOCKS * sizeof(pthread_mutex_t));
    for (i = 0; i < NO_GLU_LOCKS; ++i)
	pthread_mutex_init(&pxgstrf_shared->lu_locks[i], NULL);
#else
    pxgstrf_shared->lu_locks = (mutex_t *) SUPERLU_MALLOC(NO_GLU_LOCKS * sizeof(mutex_t));
#endif    
    
#if ( PRNTlevel==1 )
    printf(".. ParallelInit() ... nprocs %2d\n", pdgstrf_options->nprocs);
#endif

    pxgstrf_shared->spin_locks = intCalloc(n);
    pxgstrf_shared->pan_status = 
        (pan_status_t *) SUPERLU_MALLOC((n+1)*sizeof(pan_status_t));
    pxgstrf_shared->fb_cols    = intMalloc(n+1);

    panel_size = pdgstrf_options->panel_size;
    relax = pdgstrf_options->relax;
    w = MAX(panel_size, relax) + 1;
    for (i = 0; i < w; ++i) panel_histo[i] = 0;
    pxgstrf_shared->num_splits = 0;
    
    if ( (info = queue_init(&pxgstrf_shared->taskq, n)) ) {
	fprintf(stderr, "ParallelInit(): %d\n", info);
	ABORT("queue_init fails.");
    }

    /* Count children of each node in the etree. */
    for (i = 0; i <= n; ++i) pxgstrf_shared->pan_status[i].ukids = 0;
    for (i = 0; i < n; ++i) {
	dad = etree[i];
	++pxgstrf_shared->pan_status[dad].ukids;
    }

    
    /* Find the panel partitions and initialize each panel's status */

#ifdef PROFILE
    num_panels = 0;
#endif

    pxgstrf_shared->tasks_remain = 0;
    rs = 1;
    w_top = panel_size/2;
    if ( w_top == 0 ) w_top = 1;
    P = 12;

    for (i = 0; i < n; ) {
	if ( pxgstrf_relax[rs].fcol == i ) {
	    w = pxgstrf_relax[rs++].size;
	    panel_type = RELAXED_SNODE;
	    pxgstrf_shared->pan_status[i].state = CANGO;
	} else {
	    w = MIN(panel_size, pxgstrf_relax[rs].fcol - i);
#ifdef SPLIT_TOP
	    if ( !do_split ) {
	  	if ( (n-i) < panel_size * P ) do_split = 1;
	    }
	    if ( do_split && w > w_top ) { /* split large panel */
	    	w = w_top;
	    	++pxgstrf_shared->num_splits;
	    }
#endif
	    for (j = i+1; j < i + w; ++j) 
		/* Do not allow panel to cross a branch point in the etree. */
		if ( pxgstrf_shared->pan_status[j].ukids > 1 ) break;
	    w = j - i;    /* j should start a new panel */
	    panel_type = REGULAR_PANEL;
	    pxgstrf_shared->pan_status[i].state = UNREADY;
#ifdef DOMAINS
	    if ( in_domain[i] == TREE_DOMAIN ) panel_type = TREE_DOMAIN;
#endif
	}

	if ( panel_type == REGULAR_PANEL ) {
	    ++pxgstrf_shared->tasks_remain;
	    /*printf("nondomain panel %6d -- %6d\n", i, i+w-1);
	    fflush(stdout);*/
	}

	ukids = k = 0;
	for (j = i; j < i + w; ++j) {
	    pxgstrf_shared->pan_status[j].size = k--;
	    pxgstrf_shared->pan_status[j].type = panel_type;
	    ukids += pxgstrf_shared->pan_status[j].ukids;
	}
	pxgstrf_shared->pan_status[i].size = w; /* leading column */
	/* only count those kids outside the panel */
	pxgstrf_shared->pan_status[i].ukids = ukids - (w-1);
	panel_histo[w]++;
	
#ifdef PROFILE
	panstat[i].size = w;
	++num_panels;
#endif
	
	pxgstrf_shared->fb_cols[i] = i;
	i += w;
    } /* for i ... */
    
    /* Dummy root */
    pxgstrf_shared->pan_status[n].size = 1;
    pxgstrf_shared->pan_status[n].state = UNREADY;

#if ( PRNTlevel==1 )
    printf(".. Split: P %d, #nondomain panels %d\n", P, pxgstrf_shared->tasks_remain);
#endif
#ifdef DOMAINS
    EnqueueDomains(&pxgstrf_shared->taskq, list_head, pxgstrf_shared);
#else
    EnqueueRelaxSnode(&pxgstrf_shared->taskq, n, pxgstrf_relax, pxgstrf_shared);
#endif
#if ( PRNTlevel==1 )
    printf(".. # tasks %d\n", pxgstrf_shared->tasks_remain);
    fflush(stdout);
#endif

#ifdef PREDICT_OPT
    /* Set up structure describing children */
    for (i = 0; i <= n; cp_firstkid[i++] = EMPTY);
    for (i = n-1; i >= 0; i--) {
	dad = etree[i];
	cp_nextkid[i] = cp_firstkid[dad];
	cp_firstkid[dad] = i;
    }
#endif

    return 0;
} /* ParallelInit */

/*
 * Free the storage used by the parallel scheduling algorithm.
 */
int ParallelFinalize(pxgstrf_shared_t *pxgstrf_shared)
{
    /* Destroy mutexes */
#if ( MACH==SUN )
    register int i;
    for (i = 0; i < NO_GLU_LOCKS; ++i)
        mutex_destroy( &pxgstrf_shared->lu_locks[i] );
#elif ( MACH==DEC || MACH==PTHREAD )
    register int i;
    for (i = 0; i < NO_GLU_LOCKS; ++i) 
        pthread_mutex_destroy( &pxgstrf_shared->lu_locks[i] );
#endif    
    
    SUPERLU_FREE ((void*)pxgstrf_shared->lu_locks);
    SUPERLU_FREE ((int*)pxgstrf_shared->spin_locks);
    SUPERLU_FREE (pxgstrf_shared->pan_status);
    SUPERLU_FREE (pxgstrf_shared->fb_cols);
    SUPERLU_FREE (pxgstrf_shared->Glu->map_in_sup);
    queue_destroy(&pxgstrf_shared->taskq);

#if ( PRNTlevel==1 )
    printf(".. # panel splittings %d\n", pxgstrf_shared->num_splits);
#endif

    return 0;
}

int queue_init(queue_t *q, int n)
{
    if ( n < 1 ) return (-1);

    q->queue = (qitem_t *) SUPERLU_MALLOC(n*sizeof(qitem_t));
    q->count = 0;
    q->head = 0;
    q->tail = 0;

    return 0;
}

int queue_destroy(queue_t *q)
{
    SUPERLU_FREE( q->queue );
    return 0;
}

/*
 * Return value: number of items in the queue
 */
int Enqueue(queue_t *q, qitem_t item)
{
    q->queue[q->tail++] = item;
    ++q->count;
    return (q->count);
}

/*
 * Return value: >= 0 number of items in the queue
 *               = -1 queue is empty
 */
int Dequeue(queue_t *q, qitem_t *item)
{
    if ( q->count <= 0 ) return EMPTY;
    
    *item = q->queue[q->head++];
    --q->count;
    return (q->count);
}

int QueryQueue(queue_t *q)
{
    register int     i;
    printf("Queue count: %d\n", q->count);
    for (i = q->head; i < q->tail; ++i)
	printf("%8d\titem %8d\n", i, q->queue[i]);

    return 0;
}

int EnqueueRelaxSnode(queue_t *q, int n, pxgstrf_relax_t *pxgstrf_relax,
		      pxgstrf_shared_t *pxgstrf_shared)
{
    register int rs, j, m;

    m = pxgstrf_relax[0].size;
    for (rs = 1; rs <= m; ++rs) {
	j = pxgstrf_relax[rs].fcol;
	q->queue[q->tail++] = j;
	q->count++;
	++pxgstrf_shared->tasks_remain;
    }
#if ( PRNTlevel==1 )    
    printf(".. EnqueueRelaxSnode(): count %d\n", q->count);
#endif
    return 0;
}

/*
 * Enqueue the initial independent domains.
 * A pair of two numbers {root, fst_desc} is added in the queue.
 */
/*int EnqueueDomains(int P, queue_t *q, struct Branch **proc_domains_h)*/
int EnqueueDomains(queue_t *q, struct Branch *list_head,
		   pxgstrf_shared_t *pxgstrf_shared)
{
    struct Branch *b, *thrash;

/*    for (pnum = 0; pnum < P; ++pnum) {
	for (b = proc_domains_h[pnum]; b != NULL; ) {*/
    b = list_head;
    while ( b ) {
	thrash = b;
	q->queue[q->tail++] = b->root;
	q->queue[q->tail++] = b->first_desc;
	q->count = q->count + 2;
	STATE ( b->root ) = CANGO;
	++pxgstrf_shared->tasks_remain;
	b = b->next;
	SUPERLU_FREE (thrash);
    }
    printf("EnqueueDomains(): count %d\n", q->count);
    return 0;
}

int NewNsuper(const int pnum, mutex_t *lock, int *data)
{
    register int i;

#ifdef PROFILE
    double t = SuperLU_timer_();
#endif    
	
#if ( MACH==SUN )
    mutex_lock(lock);
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_lock(lock);
#elif ( MACH==SGI || MACH==ORIGIN )
#pragma critical lock(lock)
#elif  ( MACH==CRAY_PVP )
#pragma _CRI guard (*lock)
#endif    
    {
      i = ++(*data);
    }
#if ( MACH==SUN )
    mutex_unlock(lock);
#elif ( MACH==DEC || MACH==PTHREAD )
    pthread_mutex_unlock(lock);
#elif ( MACH==CRAY_PVP )
#pragma _CRI endguard (*lock)
#endif    

#ifdef PROFILE
    Gstat->procstat[pnum].cs_time += SuperLU_timer_() - t;
#endif
	
    return i;
}

int lockon(int *block)
{
    while ( *block ) ; /* spin-wait */
    *block = 1;
    return 0;
}

int lockoff(int *block)
{
    *block = 0;
    return 0;
}


