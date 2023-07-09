/*! @file
 * \brief Chooses machine-dependent parameters for the local environment
 */
/*
 * File name:		sp_ienv.c
 * History:             Modified from lapack routine ILAENV
 */
#include "superlu_ddefs.h"
#include "machines.h"

/*! \brief

</pre>
    Purpose   
    =======   

    sp_ienv_dist() is inquired to choose machine-dependent parameters for the
    local environment. See ISPEC for a description of the parameters.   

    This version provides a set of parameters which should give good,   
    but not optimal, performance on many of the currently available   
    computers.  Users are encouraged to modify this subroutine to set   
    the tuning parameters for their particular machine using the option   
    and problem size information in the arguments.   

    Arguments   
    =========   

    ISPEC   (input) int
            Specifies the parameter to be returned as the value of SP_IENV_DIST.   
            = 1: the panel size w; a panel consists of w consecutive
	         columns of matrix A in the process of Gaussian elimination.
		 The best value depends on machine's cache characters.
            = 2: the relaxation parameter relax; if the number of
	         nodes (columns) in a subtree of the elimination tree is less
		 than relax, this subtree is considered as one supernode,
		 regardless of the their row structures.
            = 3: the maximum size for a supernode, which must be greater
                 than or equal to relaxation parameter (see case 2);
	    = 4: the minimum row dimension for 2-D blocking to be used;
	    = 5: the minimum column dimension for 2-D blocking to be used;
	    = 6: the estimated fills factor for the adjacency structures 
	         of L and U, compared with A;
	    = 7: the minimum value of the product M*N*K for a GEMM call
	         to be off-loaded to accelerator (e.g., GPU, Xeon Phi).
	    
   (SP_IENV_DIST) (output) int
            >= 0: the value of the parameter specified by ISPEC   
            < 0:  if SP_IENV_DIST = -k, the k-th argument had an illegal value.
  
    ===================================================================== 
</pre>
*/


#include <stdlib.h>
#include <stdio.h>


int_t
sp_ienv_dist(int_t ispec)
{
    // printf(" this function called\n");
    int i;

    char* ttemp;

    switch (ispec) {
#if ( MACH==CRAY_T3E )
	case 2: return (6);
	case 3: return (30);

#elif ( MACH==IBM )
	case 2: return (20);
	case 3: return (100);
#else
	case 2: 
            ttemp = getenv("NREL");
            if(ttemp)
            {
                return(atoi(ttemp));
            }
            else
            return 20;
            
	case 3: 
            ttemp = getenv("NSUP");
            if(ttemp)
            {
                return(atoi(ttemp));
            }
            else
            return 128;

#endif
        case 6: return (5);
        case 7:
	    ttemp = getenv ("N_GEMM");
	    if (ttemp) return atoi (ttemp);
	    else return 10000;

    }

    /* Invalid value for ISPEC */
    i = 1;
    xerr_dist("sp_ienv", &i);
    return 0;


} /* sp_ienv_dist */

