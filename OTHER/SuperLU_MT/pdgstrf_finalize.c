#include "pdsp_defs.h"
#include "util.h"

void
pdgstrf_finalize(pdgstrf_options_t *pdgstrf_options, SuperMatrix *AC)
{
/*
 * -- SuperLU MT routine (version 1.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * August 15, 1997
 *
 * Purpose
 * =======
 * 
 * pdgstrf_finalize() deallocates storage after factorization pdgstrf().
 *
 * Arguments
 * =========
 *
 * pdgstrf_options (input) pdgstrf_options_t*
 *        The structure contains the parameters to facilitate sparse
 *        LU factorization.
 *
 * AC (input) SuperMatrix*
 *        The original matrix with columns permuted.
 */
    SUPERLU_FREE(pdgstrf_options->etree);
    SUPERLU_FREE(pdgstrf_options->colcnt_h);
    SUPERLU_FREE(pdgstrf_options->part_super_h);
    Destroy_CompCol_Permuted(AC);
#if ( DEBUGlevel>=1 )
    printf("** pdgstrf_finalize() called\n");
#endif
}
