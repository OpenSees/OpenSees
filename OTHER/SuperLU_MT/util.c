/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 */
#include <unistd.h>
#include <math.h>
#include "pdsp_defs.h"
#include "util.h"

void superlu_abort_and_exit(char* msg)
{
    fprintf(stderr, msg);
    exit (-1);
}

void *superlu_malloc(int size)
{
    void *buf;
    buf = (void *) malloc(size);
    return (buf);
}

void superlu_free(void *addr)
{
    free (addr);
}

/* Deallocate the structure pointing to the actual storage of the matrix. */
void
Destroy_SuperMatrix_Store(SuperMatrix *A)
{
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperMatrix_Store ...\n");
#endif
}

/* A is of type Stype==NC */
void
Destroy_CompCol_Matrix(SuperMatrix *A)
{
    SUPERLU_FREE( ((NCformat *)A->Store)->rowind );
    SUPERLU_FREE( ((NCformat *)A->Store)->colptr );
    SUPERLU_FREE( ((NCformat *)A->Store)->nzval );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_Matrix ...\n");
#endif
}

/* A is of type Stype==NCP */
void
Destroy_CompCol_Permuted(SuperMatrix *A)
{
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colbeg );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colend );
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_Permuted ...\n");
#endif
}

/* A is of type Stype==NCP */
void
Destroy_CompCol_NCP(SuperMatrix *A)
{
    SUPERLU_FREE ( ((NCPformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colbeg );
    SUPERLU_FREE ( ((NCPformat *)A->Store)->colend );
    SUPERLU_FREE ( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_CompCol_NCP ...\n");
#endif
}

/* A is of type Stype==SC */
void
Destroy_SuperNode_Matrix(SuperMatrix *A)
{
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCformat *)A->Store)->rowind_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((SCformat *)A->Store)->nzval_colptr );
    SUPERLU_FREE ( ((SCformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCformat *)A->Store)->sup_to_col );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperNode_Matrix ...\n");
#endif
}

/* A is of type Stype==SCP */
void
Destroy_SuperNode_SCP(SuperMatrix *A)
{
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->nzval_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colend );
    SUPERLU_FREE( A->Store );
#if ( PRNTlevel==1 )
    printf(".. Destroy_SuperNode_SCP ...\n");
#endif
}

/*
 * Reset repfnz[*] for the current column 
 */
void
pxgstrf_resetrep_col(const int nseg, const int *segrep, int *repfnz)
{
    register int i, irep;
    
    for (i = 0; i < nseg; ++i) {
	irep = segrep[i];
	repfnz[irep] = EMPTY;
    }
}


/*
 * Count the total number of nonzeros in factors L and U,  and in the 
 * symmetrically reduced L. 
 */
void
countnz(const int n, int *xprune, int *nnzL, int *nnzU, GlobalLU_t *Glu)
{
    register int nsuper, fsupc, i, j, nnzL0, jlen, irep;
    register int nnzsup = 0;
    register int *xsup, *xsup_end, *xlsub, *xlsub_end, *supno;
	
    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    supno     = Glu->supno;
    *nnzU     = Glu->nextu;
    nnzL0     = 0;
    *nnzL     = 0;
    nsuper    = supno[n];

    if ( n <= 0 ) return;

    /* 
     * For each supernode ...
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jlen = xlsub_end[fsupc] - xlsub[fsupc];
	nnzsup += jlen * (xsup_end[i] - fsupc);
			  
	for (j = fsupc; j < xsup_end[i]; j++) {
	    *nnzL += jlen;
	    *nnzU += j - fsupc + 1;
	    jlen--;
	}
	irep = SUPER_REP(i);
	if ( SINGLETON(supno[irep]) )
	    nnzL0 += xprune[irep] - xlsub_end[irep];
	else 
	    nnzL0 += xprune[irep] - xlsub[irep];
    }

#if ( PRNTlevel==1 )
    printf(".. # supernodes = %d\n", nsuper+1);
    printf(".. # edges in symm-reduced L = %d\n", nnzL0);
    if ( Glu->dynamic_snode_bound )
      printf(".. # NZ in LUSUP %d, dynamic bound %d, utilization %.2f\n",
	     nnzsup, Glu->nextlu, (float)nnzsup/Glu->nextlu);
    else
      printf(".. # NNZ in LUSUP %d, static bound %d, utilization %.2f\n",
	     nnzsup, Glu->nzlumax, (float)nnzsup/Glu->nzlumax);
#endif
}



/*
 * Fix up the data storage lsub for L-subscripts. It reclaims the
 * storage for the adjancency lists of the pruned graph, and applies
 * row permuation to the row subscripts of matrix $L$.
 */
void
fixupL(const int n, const int *perm_r, GlobalLU_t *Glu)
{
    register int nsuper, fsupc, nextl, i, j, jstrt;
    register int *xsup, *xsup_end, *lsub, *xlsub, *xlsub_end;

    if ( n <= 1 ) return;

    xsup      = Glu->xsup;
    xsup_end  = Glu->xsup_end;
    lsub      = Glu->lsub;
    xlsub     = Glu->xlsub;
    xlsub_end = Glu->xlsub_end;
    nsuper    = Glu->supno[n];
    nextl     = 0;
    
    /* 
     * For each supernode ...
     */
    for (i = 0; i <= nsuper; i++) {
	fsupc = xsup[i];
	jstrt = xlsub[fsupc];
	xlsub[fsupc] = nextl;
	for (j = jstrt; j < xlsub_end[fsupc]; j++) {
	    lsub[nextl] = perm_r[lsub[j]]; /* Now indexed into P*A */
	    nextl++;
  	}
	xlsub_end[fsupc] = nextl;
    }
    xlsub[n] = nextl;
#if ( PRNTlevel==1 )
    printf(".. # edges in supernodal graph of L = %d\n", nextl);
    fflush(stdout);
#endif
}

/*
 * Print all definitions to be used by CPP.
 */
int cpp_defs()
{
    printf("CPP Defs:\n");
#ifdef PRNTlevel
    printf("\tPRNTlevel=%d\n", PRNTlevel);
#endif 
#ifdef DEBUGlevel
    printf("\tDEBUGlevel=%d\n", DEBUGlevel);
#endif 
#ifdef PROFILE
    printf("\tPROFILE\n");
#endif
#ifdef PREDICT_OPT
    printf("\tPREDICT_OPT\n");
#endif
#ifdef USE_VENDOR_BLAS
    printf("\tUSE_VENDOR_BLAS\n");
#endif
#ifdef GEMV2
    printf("\tGEMV2\n");
#endif
#ifdef SCATTER_FOUND
    printf("\tSCATTER_FOUND\n");
#endif

    return 0;
}

/*
 * Compress the data storage LUSUP[*] for supernodes. It removes the
 * memory holes due to untightness of the upper bounds by A'A-supernode.
 */
void
compressSUP(const int n, GlobalLU_t *Glu)
{
    register int nextlu, i, j, jstrt;
    int *xlusup, *xlusup_end;
    double *lusup;

    if ( n <= 1 ) return;

    lusup     = Glu->lusup;
    xlusup    = Glu->xlusup;
    xlusup_end= Glu->xlusup_end;
    nextlu     = 0;
    
    for (j = 0; j < n; ++j) {
	jstrt = xlusup[j];
	xlusup[j] = nextlu;
	for (i = jstrt; i < xlusup_end[j]; ++i, ++nextlu)
	    lusup[nextlu] = lusup[i];
	xlusup_end[j] = nextlu;
    }
    xlusup[n] = nextlu;
    printf("\tcompressSUP() nextlu %d\n", nextlu);
}

int check_mem_leak(char *where)
{ 
    void *addr;
    addr = (void *)sbrk(0);
    printf("\tsbrk(0) %s: addr = %u\n", where, addr);
    return 0;
}

/*
 * Diagnostic print of segment info after pdgstrf_panel_dfs().
 */
void print_panel_seg(int n, int w, int jcol, int nseg, 
		     int *segrep, int *repfnz)
{
    int j, k;
    
    for (j = jcol; j < jcol+w; j++) {
	printf("\tcol %d:\n", j);
	for (k = 0; k < nseg; k++)
	    printf("\t\tseg %d, segrep %d, repfnz %d\n", k, 
			segrep[k], repfnz[(j-jcol)*n + segrep[k]]);
    }
}

/*
 * Allocate storage for various statistics.
 */
void
StatAlloc(const int n, const int nprocs, const int panel_size, 
	  const int relax, Gstat_t *Gstat)
{
    register int w;

    w = MAX( panel_size, relax ) + 1;
    Gstat->panel_histo = intCalloc(w);
    Gstat->utime = (double *) doubleMalloc(NPHASES);
    Gstat->ops   = (flops_t *) SUPERLU_MALLOC(NPHASES * sizeof(flops_t));
    
    if ( !(Gstat->procstat 
		= (procstat_t *) SUPERLU_MALLOC(nprocs*sizeof(procstat_t))) )
	ABORT( "SUPERLU_MALLOC failed for procstat[]" );

#if (PRNTlevel==1)
    printf(".. StatAlloc(): n %d, nprocs %d, panel_size %d, relax %d\n",
		n, nprocs, panel_size, relax);
#endif
#ifdef PROFILE    
    if ( !(panstat = (panstat_t *) SUPERLU_MALLOC(n * sizeof(panstat_t))) )
	ABORT( "SUPERLU_MALLOC failed for panstat[]" );
    panhows = intCalloc(3);
    Gstat->height = intCalloc(n+1);
    if ( !(flops_by_height = (float *) SUPERLU_MALLOC(n * sizeof(float))) )
	ABORT("SUPERLU_MALLOC failed for flops_by_height[]");
    
#endif
    
#ifdef PREDICT_OPT
    if ( !(cp_panel = (cp_panel_t *) SUPERLU_MALLOC(n * sizeof(cp_panel_t))) )
	ABORT( "SUPERLU_MALLOC failed for cp_panel[]" );
    if ( !(desc_eft = (desc_eft_t *) SUPERLU_MALLOC(n * sizeof(desc_eft_t))) )
	ABORT( "SUPERLU_MALLOC failed for desc_eft[]" );
    cp_firstkid = intMalloc(n+1);
    cp_nextkid = intMalloc(n+1);
#endif
    
}

/*
 * Initialize various statistics variables.
 */
void
StatInit(const int n, const int nprocs, Gstat_t *Gstat)
{
    register int i;
    
    for (i = 0; i < NPHASES; ++i) {
	Gstat->utime[i] = 0;
	Gstat->ops[i] = 0;
    }
    
    for (i = 0; i < nprocs; ++i) {
	Gstat->procstat[i].panels = 0;
	Gstat->procstat[i].fcops = 0.0;
	Gstat->procstat[i].skedwaits = 0;
	Gstat->procstat[i].skedtime = 0.0;
	Gstat->procstat[i].cs_time = 0.0;
	Gstat->procstat[i].spintime = 0.0;
	Gstat->procstat[i].pruned = 0;
	Gstat->procstat[i].unpruned = 0;
    }

#ifdef PROFILE    
    for (i = 0; i < n; ++i) {
	panstat[i].fctime = 0.0;
	panstat[i].flopcnt = 0.0;
	panstat[i].pipewaits = 0;
	panstat[i].spintime = 0.0;
	flops_by_height[i] = 0.0;
    }
    for (i = 0; i < 3; ++i) panhows[i] = 0;
    dom_flopcnt = 0.;
    flops_last_P_panels = 0;
#endif
    
#ifdef PREDICT_OPT
    for (i = 0; i < n; ++i)
	cp_panel[i].est = cp_panel[i].pdiv = 0;
#endif
    
#if ( PRNTlevel==1 )
    printf(".. StatInit(): n %d, nprocs %d\n", n, nprocs);
#endif
}


/* Print timings used in factorization and solve. */
void
PrintStat(Gstat_t *Gstat)
{
    double         *utime;
    flops_t        *ops;

    utime = Gstat->utime;
    ops   = Gstat->ops;
    printf("Factor time  = %8.2f\n", utime[FACT]);
    if ( utime[FACT] != 0.0 )
      printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	     ops[FACT]*1e-6/utime[FACT]);

    printf("Solve time   = %8.2f\n", utime[SOLVE]);
    if ( utime[SOLVE] != 0.0 )
      printf("Solve flops = %e\tMflops = %8.2f\n", ops[SOLVE],
	     ops[SOLVE]*1e-6/utime[SOLVE]);

}

void
StatFree(Gstat_t *Gstat)
{
    SUPERLU_FREE (Gstat->panel_histo);
    SUPERLU_FREE (Gstat->utime);
    SUPERLU_FREE (Gstat->ops);
    SUPERLU_FREE (Gstat->procstat);
    
#ifdef PROFILE
    SUPERLU_FREE (Gstat->panstat);
    SUPERLU_FREE (Gstat->panhows);
    SUPERLU_FREE (Gstat->height);
    SUPERLU_FREE (flops_by_height);
#endif

#ifdef PREDICT_OPT
    SUPERLU_FREE (Gstat->cp_panel);
    SUPERLU_FREE (Gstat->desc_eft);
    SUPERLU_FREE (Gstat->cp_firstkid);
    SUPERLU_FREE (Gstat->cp_nextkid);
#endif    

#if (PRNTlevel==1)
    printf(".. StatFree(): Free Stat variables.\n");
#endif
}

flops_t
LUFactFlops(Gstat_t *Gstat)
{
    return (Gstat->ops[FACT]);
}

flops_t
LUSolveFlops(Gstat_t *Gstat)
{
    return (Gstat->ops[SOLVE]);
}



/* 
 * Fills an integer array with a given value.
 */
void ifill(int *a, int alen, int ival)
{
    register int i;
    for (i = 0; i < alen; i++) a[i] = ival;
}



/* 
 * Get the statistics of the supernodes 
 */
#define NBUCKS 10
static 	int	max_sup_size;

void super_stats(int nsuper, int *xsup, int *xsup_end)
{
    register int nsup1 = 0;
    int          i, isize, whichb, bl, bh;
    int          bucket[NBUCKS];

    max_sup_size = 0;

    /* Histogram of the supernode sizes */
    ifill (bucket, NBUCKS, 0);

    for (i = 0; i <= nsuper; i++) {
        isize = xsup_end[i] - xsup[i];
	if ( isize == 1 ) nsup1++;
	if ( max_sup_size < isize ) max_sup_size = isize;	
        whichb = (float) isize / max_sup_size * NBUCKS;
        if (whichb >= NBUCKS) whichb = NBUCKS - 1;
        bucket[whichb]++;
    }
    
    printf("** Supernode statistics:\n\tno of supernodes = %d\n", nsuper+1);
    printf("\tmax supernode size = %d\n", max_sup_size);
    printf("\tno of size 1 supernodes = %d\n", nsup1);

    printf("\tHistogram of supernode size:\n");
    for (i = 0; i < NBUCKS; i++) {
        bl = (float) i * max_sup_size / NBUCKS;
        bh = (float) (i+1) * max_sup_size / NBUCKS;
        printf("\t%3d-%3d\t\t%d\n", bl+1, bh, bucket[i]);
    }

}

void panel_stats(int n, int max_w, int* in_domain, Gstat_t *Gstat)
{
    register int i, w;
    float *histo_flops, total;

    histo_flops = (float *) SUPERLU_MALLOC( max_w * sizeof(float) );

    for (i = 0; i < max_w; ++i) histo_flops[i] = 0;
    total = 0;
    for (i = 0; i < n; i += w) {
	w = Gstat->panstat[i].size;
	if ( in_domain[i] != TREE_DOMAIN ) {
	    histo_flops[w - 1] += Gstat->panstat[i].flopcnt;
	    total += Gstat->panstat[i].flopcnt;
	}
    }

    if ( total != 0.0 ) {
	printf("** Panel & flops distribution: nondomain flopcnt %e\n", total);
	for (i = 1; i <= max_w; i++)
	    printf("\t%d\t%d\t%e (%.2f)\n", i, Gstat->panel_histo[i],
		   histo_flops[i-1], histo_flops[i-1]/total);
    }
    SUPERLU_FREE (histo_flops);
}



float SpaSize(int n, int np, float sum_npw)
{
    return (sum_npw*8 + np*8 + n*4)/1024.;
}

float DenseSize(int n, float sum_nw)
{
    return (sum_nw*8 + n*8)/1024.;;
}


/*
 * Check whether repfnz[] == EMPTY after reset.
 */
void check_repfnz(int n, int w, int jcol, int *repfnz)
{
    int jj, k;

    for (jj = jcol; jj < jcol+w; jj++) 
	for (k = 0; k < n; k++)
	    if ( repfnz[(jj-jcol)*n + k] != EMPTY ) {
		fprintf(stderr, "col %d, repfnz_col[%d] = %d\n", jj,
			k, repfnz[(jj-jcol)*n + k]);
		ABORT("repfnz[] not empty.");
	    }
}


int PrintInt10(char *name, int len, int *x)
{
    register int i;
    
    printf("(len=%d) %s:", len, name);
    for (i = 0; i < len; ++i) {
	if ( i % 10 == 0 ) printf("\n[%4d-%4d]", i, i+9);
	printf("%6d", x[i]);
    }
    printf("\n");
    return 0;
}

/* Print a summary of the testing results. */
void
PrintSumm(char *type, int nfail, int nrun, int nerrs)
{
    if ( nfail > 0 )
	printf("%3s driver: %d out of %d tests failed to pass the threshold\n",
	       type, nfail, nrun);
    else
	printf("All tests for %3s driver passed the threshold (%6d tests run)\n", type, nrun);

    if ( nerrs > 0 )
	printf("%6d error messages recorded\n", nerrs);
}


/* Print the adjacency list for graph of L, including the pruned graph,
   graph of U, and L supernodes partition */
int PrintGLGU(int n, int *xprune, GlobalLU_t *Glu)
{
    register int nsuper = Glu->nsuper;
    PrintInt10("LSUB", Glu->xlsub_end[n-1], Glu->lsub);
    PrintInt10("XLSUB", n, Glu->xlsub);
    PrintInt10("XLSUB_END", n, Glu->xlsub_end);
    PrintInt10("XPRUNE", n, xprune);
    PrintInt10("USUB", Glu->xusub_end[n-1], Glu->usub);
    PrintInt10("XUSUB", n, Glu->xusub);
    PrintInt10("XUSUB_END", n, Glu->xusub_end);
    PrintInt10("SUPNO", n, Glu->supno);
    PrintInt10("XSUP", nsuper+1, Glu->xsup);
    PrintInt10("XSUP_END", nsuper+1, Glu->xsup_end);
    return 0;
}

#if 0
/*
 * Print the statistics of the relaxed snodes for matlab process
 */
void relax_stats(int start, int end, int step)
{
    FILE *fp;
    int i;

    fp = fopen("relax.m", "w");
    
    fprintf(fp,"relax = [\n");
    for (i = start; i <= end; i += step) fprintf(fp, "%d ", i);
    fprintf(fp, "];\n");

    fprintf(fp, "fctime = [\n");
    for (i = start; i <= end; i += step) 
	fprintf(fp, "%15.8e\n ", stat_relax[i].fctime);
    fprintf(fp, "];\n");

    fprintf(fp, "mflops = [\n");
    for (i = start; i <= end; i += step)
	fprintf(fp, "%15.8e\n ", (float)stat_relax[i].flops / 1e6);
    fprintf(fp, "];\n");

    fprintf(fp, "mnzs = [\n");
    for (i = start; i <= end; i += step)
 	fprintf(fp, "%15.8e\n ", stat_relax[i].nzs / 1e6);
    fprintf(fp, "];\n");

    fclose(fp);
}

/*
 * Obtain the distribution of time/flops/nzs on the snode size.
 */
void snode_profile(int nsuper, int *xsup)
{
    FILE *fp;
    int i, j;
    int ssize;

    if ( !(stat_snode = (stat_snode_t *) SUPERLU_MALLOC((max_sup_size+1) *
	sizeof(stat_snode_t))) ) ABORT("SUPERLU_MALLOC fails for stat_snode[].");

    for (i = 0; i <= max_sup_size; i++) {
	stat_snode[i].ncols = 0;
	stat_snode[i].flops = 0;
	stat_snode[i].nzs = 0;
	stat_snode[i].fctime = 0.0;
    }	

    for (i = 0; i <= nsuper; i++) {

	ssize = xsup[i+1] - xsup[i];   
	stat_snode[ssize].ncols += ssize;

        for (j=xsup[i]; j<xsup[i+1]; j++) { 
	    stat_snode[ssize].flops += stat_col[j].flops;	    
	    stat_snode[ssize].nzs += stat_col[j].nzs;	    
	    stat_snode[ssize].fctime += stat_col[j].fctime;	    
	}

    }

    fp = fopen("snode.m", "w");
    
    fprintf(fp, "max_sup_size = %d;\n", max_sup_size);

    fprintf(fp,"ncols = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%d ", stat_snode[i].ncols);
    fprintf(fp, "];\n");

    fprintf(fp, "fctime = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", stat_snode[i].fctime);
    fprintf(fp, "];\n");

    fprintf(fp, "mflops = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", (float) stat_snode[i].flops / 1e6);
    fprintf(fp, "];\n");

    fprintf(fp, "mnzs = [");
    for (i = 1; i <= max_sup_size; i++) 
	fprintf(fp, "%15.8e\n", (float) stat_snode[i].nzs / 1e6);
    fprintf(fp, "];\n");

    fclose(fp);

    SUPERLU_FREE (stat_snode);

}
#endif

int print_int_vec(char *what, int n, int *vec)
{
    int i;
    printf("%s\n", what);
    for (i = 0; i < n; ++i) printf("%d\t%d\n", i, vec[i]);
    return 0;
}


/*
 * Print the parallel execution statistics.
 */
int ParallelProfile(const int n, const int supers, const int panels, 
		const int procs, Gstat_t *Gstat)
{
    register int i, imax, pruned, unpruned, waits, itemp, cs_numbers;
    register float loadmax, loadtot, temp, thresh, loadprint;
    register float waittime, cs_time;
    double    *utime = Gstat->utime;
    procstat_t *pstat;
    panstat_t *pan;
    void print_flops_by_height(int, panstat_t *, int *, float *);
    
    printf("\n---- Parallel Profile Per Processor ----\n");
    printf("%4s%16s%8s%10s%10s%10s%10s%8s\n", "proc", "factops",
	   "seconds", "skedwaits", "skedtime", "CS-time",
	   "spin-time", "[%tot]");
    for (i = 0; i < procs; ++i) {
	pstat = &(Gstat->procstat[i]);
	if ( pstat->fctime != 0 ) {
	    temp = pstat->spintime/pstat->fctime*100.;
	    printf("%4d%16e%8.2f%10d%10.3f%10.3f%10.3f%8.1f\n", 
		   i, pstat->fcops, pstat->fctime, pstat->skedwaits,
		   pstat->skedtime, pstat->cs_time, pstat->spintime, temp);
	}
    }

    printf("%4s%8s%12s%14s\n",
	   "proc", "#panels", "dfs_pruned","dfs_unpruned");
    pruned = unpruned = 0;
    cs_time = 0.0;
    for (i = 0; i < procs; ++i) {
	pstat = &(Gstat->procstat[i]);
	printf("%4d%8d%12d%14d\n", i, pstat->panels,
		pstat->pruned, pstat->unpruned);
	pruned += Gstat->procstat[i].pruned;
	unpruned += Gstat->procstat[i].unpruned;
	cs_time += Gstat->procstat[i].cs_time;
    }
    temp = pruned + unpruned;
    if ( temp != 0 ) {
    	printf("%12s%26s\n", "", "--------------------");
    	printf("%12s%12d%14d%14.0f\n", "total", pruned, unpruned, temp);
    	printf("%12s%12.2f%14.2f\n", "frac.", pruned/temp, unpruned/temp);
    }

    printf("%16s%16d\n", "piped-panels", Gstat->panhows[PIPE]);
    printf("%16s%16d\n", "nonpiped-DADs", Gstat->panhows[DADPAN]);
    printf("%16s%16d\n", "nonpiped-panels", Gstat->panhows[NOPIPE]);

    /* work load distribution */
    loadmax = loadtot = Gstat->procstat[0].fcops;
    imax = 0;
    for (i = 1; i < procs; ++i) {
	temp = Gstat->procstat[i].fcops;
	loadtot += temp;
	if ( temp > loadmax ) {
	    loadmax = temp;
	    imax = i;
	}
    }
    printf("%25s%8.2f\n", "Load balance [mean/max]", loadtot/loadmax/procs);

    /* Delays due to pipelining. */
    waits = waittime = 0;
    for (i = 0; i < n; i += Gstat->panstat[i].size) { /* add up all panels */
	waits += Gstat->panstat[i].pipewaits;
	waittime += Gstat->panstat[i].spintime;
    }
    printf("%25s%8d,\tper-panel %.1f\n", "total #delays in pipeline",
	    waits, (float)waits/panels);
    temp = waittime / procs;
    printf("%25s%8.2f\t[%.2f]\n", "mean spin time per-proc", 
	   temp, temp/utime[FACT]);
    
    /* Delays due to scheduling. */
    waits = waittime = 0;
    for (i = 0; i < procs; ++i) {
	waits += Gstat->procstat[i].skedwaits;
	waittime += Gstat->procstat[i].skedtime;
    }
    printf("%25s%8d\n", "total #delays in schedule", waits);
    temp = waittime / procs;
    printf("%25s%8.2f\t[%.2f]\n", "mean sched time per-proc", 
	   temp, temp/utime[FACT]);

    /* estimated overhead in spin-locks */
#if ( MACH==CRAY_PVP )    /* measured for mutex lock/unlock on 4 cpus */
#define TMUTEX          4.42e-6
#define FLOPS_PER_LOCK  221
#elif ( MACH==SUN )
#define TMUTEX          4.36e-6
#define FLOPS_PER_LOCK  109
#elif ( MACH==SGI || MACH==ORIGIN )
#define TMUTEX          2.02e-6
#define FLOPS_PER_LOCK  364
#elif ( MACH==DEC || PTHREAD )
#define TMUTEX          2.71e-6
#define FLOPS_PER_LOCK  407
#endif
    cs_numbers = n + 3*supers + panels + procs; 
    itemp = cs_numbers * FLOPS_PER_LOCK;     /* translated to flops */
    temp = cs_numbers * TMUTEX;
    printf("mutex-lock overhead (est.) %8.2f, #locks %d, equiv. flops %e\n", 
	   temp, cs_numbers, (float) itemp);
    printf("time in critical section   %8.2f\t[%.2f]\n",
	   cs_time/procs, cs_time/procs/utime[FACT]);

    printf("\n---- Parallel Profile Per Panel ----\n");
    printf("%8s%8s%16s%8s%8s%12s%8s\n", "panel", "height",
	    "factops", "[tot%]", "msec", "spin(msec)", "Mflops");
    thresh = 0.005 * loadtot;
    loadprint = 0;
    itemp = 0;
    for (i = 0; i < n; i += Gstat->panstat[i].size) {
	pan = &(Gstat->panstat[i]);
	if ( pan->flopcnt > thresh ) {
	    loadprint += pan->flopcnt;
	    ++itemp;
	    if ( pan->fctime != 0 ) temp = pan->flopcnt/pan->fctime*1e-6;
	    printf("%4d%4d%8d%16e%8.1f%8.2f%12.2f%8.2f\n", i, pan->size,
		    Gstat->height[i], pan->flopcnt, pan->flopcnt/loadtot*100.,
		    pan->fctime*1e3, pan->spintime*1e3, temp);
	}
    }
    printf("Total panels %d,  height(T) %d, height(T)/n= %.4f\n", 
	   panels, Gstat->height[n], (float)Gstat->height[n]/n);
    printf("Printed flops %e [%.1f], printed panels %d [%.1f]\n",
	    loadprint, loadprint/loadtot*100.,
	    itemp, (float)itemp/panels);

/*    print_flops_by_height(n, panstat, height, flops_by_height);*/
	
    printf("---- End ParallelProfile().\n\n");
    fflush(stdout);
    return 0;
}

/*
 * Print the distribution of flops by the height of etree.
 */
void
print_flops_by_height(int n, panstat_t *panstat,
		      int *height, float *flops_by_height)
{
    register int i, w, ht;
    register float flops;

    for (i = 0; i < n; i += w) {
	w = panstat[i].size;
	ht = height[i];
	flops_by_height[ht] += panstat[i].flopcnt;
    }

    printf("\n%8s\t%8s\n", "height", "flops");
    ht = height[n-1]; /* root */
    for (i = 0; i <= ht; ++i) {
	flops = flops_by_height[i];
	if ( flops != 0.0 )
	    printf("%8d\t%e\n", i, flops);
    }
}

   
/*
 * Print the analysis of the optimal runtime.
 */
int
CPprofile(const int n, cp_panel_t *cp_panel, pxgstrf_shared_t *pxgstrf_shared)
{
    Gstat_t *Gstat = pxgstrf_shared->Gstat;
    register int maxpan, i, j, treecnt;
    register float eft, maxeft; /* earliest (possible) finish time */
    flops_t  *ops;

    /* Find the longest (weighted) path in the elimination forest. */
    treecnt = 0;
    maxeft = 0;
    for (i = Gstat->cp_firstkid[n]; i != EMPTY; i = Gstat->cp_nextkid[i]) {
/*	printf("Root %d, height %d\n", i, height[i]);*/
	j = (pxgstrf_shared->pan_status[i].size > 0) ? 
	  i : (i + pxgstrf_shared->pan_status[i].size);
	eft   = cp_panel[j].est + cp_panel[j].pdiv;
	if ( eft > maxeft ) {
	    maxeft = eft;
	    maxpan = j;
	}
	++treecnt;
    }
    
    ops   = Gstat->ops;
    printf("\n** Runtime prediction model: #trees %d\n", treecnt);
    printf("Last panel %d, seq-time %e, EFT %e, ideal-speedup %.2f\n",
	   maxpan, ops[FACT], maxeft, ops[FACT]/maxeft);

#if ( DEBUGlevel>=2 )
    printf("last panel %d\n", maxpan);
    for (i = 0; i < n; i += pxgstrf_shared->pan_status[i].size)
	printf("%6d %8s%e\t%8s%8.0f\n", i, "est  ", cp_panel[i].est,
	       "pdiv  ", cp_panel[i].pdiv);
#endif    
    return 0;
}


/***************************************************************
 * Utilities to print the supermatrix.
 ***************************************************************/

#define PERLINE  10
#define FMT      "%7.4f "

/* A is of type Stype==SCP */
void
Print_SuperNode_SCP(SuperMatrix *A)
{
    int i, j, c;
    int n = A->ncol;
    SCPformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int *colbeg = Astore->nzval_colbeg, *colend = Astore->nzval_colend;
    printf("SuperNode_SCP: nnz %d, nsuper %d\n", Astore->nnz, Astore->nsuper);
    printf("valL=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colbeg[j]; i < colend[j]; ++i) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	    ++c;
	}
    }
    printf("];\n");
    fflush(stdout);
    /*    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->rowind_colend );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->col_to_sup );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colbeg );
    SUPERLU_FREE ( ((SCPformat *)A->Store)->sup_to_colend );*/
}

/* A is of type NCP */
void
Print_CompCol_NC(SuperMatrix *A)
{
    int i, j, c;
    int n = A->ncol;
    NCformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int *colptr = Astore->colptr;
    printf("CompCol_NC: nnz %d\n", Astore->nnz);
    printf("valA=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colptr[j]; i < colptr[j+1]; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	}
    }
    printf("];\n");
    fflush(stdout);
}

/* A is of type NCP */
void
Print_CompCol_NCP(SuperMatrix *A)
{
    int i, j, c;
    int n = A->ncol;
    NCPformat *Astore = A->Store;
    double *nzval = Astore->nzval;
    int *colbeg = Astore->colbeg, *colend = Astore->colend;
    printf("SuperNode_NCP: nnz %d\n", Astore->nnz);
    printf("nzval[U]\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = colbeg[j]; i < colend[j]; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i]);
	}
    }
    printf("\n");
    fflush(stdout);
}

/* A is of type DN */
void
Print_Dense(SuperMatrix *A)
{
    int i, j, c;
    int m = A->nrow, n = A->ncol;
    DNformat *Astore = A->Store;
    int lda = Astore->lda;
    double *nzval = Astore->nzval;
    printf("Dense: lda %d\n", lda);
    printf("val=[\n");
    for (c = 0, j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i, ++c) {
	    if (c == PERLINE) { printf("\n"); c = 0; }
	    printf(FMT, nzval[i + j*lda]);
      }
    }
    printf("];\n");
    fflush(stdout);
}

