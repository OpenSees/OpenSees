/*! @file
 * \brief
 *
 * <pre>
 * -- Distributed SuperLU routine (version 1.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * September 1, 1999
 * </pre>
 */

#include "superlu_ddefs.h"

void pxerbla(char *srname, gridinfo_t *grid, int_t info)
{
    printf("{" IFMT "," IFMT "}: On entry to %6s, parameter number " IFMT " had an illegal value\n",
	   MYROW(grid->iam, grid), MYCOL(grid->iam, grid), srname, info);

}
