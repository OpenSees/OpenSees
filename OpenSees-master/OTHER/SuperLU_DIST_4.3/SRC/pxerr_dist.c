/*! @file
 * \brief
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.3) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 *
 * Modified: November 21, 1999
 *
 * </pre>
 */

#include "superlu_ddefs.h"

/* pxerbla */
void pxerr_dist(char *srname, gridinfo_t *grid, int_t info)
{
    printf("{" IFMT "," IFMT "}: On entry to %6s, parameter number " IFMT " had an illegal value\n",
	   MYROW(grid->iam, grid), MYCOL(grid->iam, grid), srname, info);

}
