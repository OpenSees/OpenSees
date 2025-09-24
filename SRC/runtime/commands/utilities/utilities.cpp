#include <tcl.h>
#include "utilities/linalg.hh"

int init_g3_tcl_utils(Tcl_Interp* interp)
{ 
  const int libsize = sizeof(linalg)/sizeof(char*);
  for (int i=0; i < libsize; i++)
    Tcl_Eval(interp, linalg[i]);
  return 0;
}
