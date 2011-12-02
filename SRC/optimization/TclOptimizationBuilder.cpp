


#include "TclOptimizationBuilder.h"
static Domain *theStructuralDomain = 0;
TclOptimizationBuilder::TclOptimizationBuilder( Domain &passedSructuralDomain, Tcl_Interp *interp )
{ 
   theInterp = interp;
   theStructuralDomain	= &passedSructuralDomain;
 } 
TclOptimizationBuilder::~TclOptimizationBuilder()
{
}