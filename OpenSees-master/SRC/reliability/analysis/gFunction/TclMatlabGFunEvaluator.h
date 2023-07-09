#ifndef TclMatlabGFunEvaluator_h
#define TclMatlabGFunEvaluator_h

#include <GFunEvaluator.h>
#include <Vector.h>
#include <ReliabilityDomain.h>
#include <tcl.h>


class TclMatlabGFunEvaluator : public GFunEvaluator
{

public:
	TclMatlabGFunEvaluator(	Tcl_Interp *passedTclInterp,
						ReliabilityDomain *passedReliabilityDomain);
	~TclMatlabGFunEvaluator();

	int		evaluate_g(const Vector &passed_x);
	double	get_g();

protected:

private:
	double g;
	Tcl_Interp *theTclInterp;
	ReliabilityDomain *theReliabilityDomain;
};

#endif
