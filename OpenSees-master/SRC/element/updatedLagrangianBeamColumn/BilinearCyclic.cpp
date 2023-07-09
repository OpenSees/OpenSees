#include "BilinearCyclic.h"
#include <math.h>

BilinearCyclic::BilinearCyclic(int tag, double weight)
:CyclicModel(tag, -1), weightFactor(weight)
{

}


BilinearCyclic::~BilinearCyclic()
{
  // does nothing
}


CyclicModel *BilinearCyclic::getCopy()
{
CyclicModel *newModel = new BilinearCyclic(getTag(), weightFactor);	
	return newModel;
}

double BilinearCyclic::getTaskFactor()
{
double tfactor;
	// redundant - only for print
	if(d_curr >= 0 && !initYieldPos)
		return 1.0;
	if(d_curr  < 0 && !initYieldNeg)
		return 1.0;
	// end redundant

	if(yielding /* && fabs(d_curr) >= fabs(d_end) */)
//		return resFactor; // will eventually unload
		tfactor = cycFactor_hist;
    else
    {
    	if(f_bgn*f_end < 0) // full-cycle
    	{
    		if(contains(0.0, f_bgn, f_curr))
    			tfactor=resFactor;
    		else
			  {
				  tfactor=rationalize(d_curr, f_curr, d_end, f_end); // actually with (x0, 0.0)
				  tfactor = weightFactor*tfactor + (1 - weightFactor)*resFactor;
			  }
				  
    	}
    	else // half-cycle
		 {
			 tfactor = rationalize(d_bgn, f_bgn, d_end, f_end);
			 tfactor = weightFactor*tfactor + (1 - weightFactor)*resFactor;
         }
    }

    // opserr << *this; opserr << "\a";
	// opserr << "tFactor = " << tfactor << endln;
	// if(yielding) {cout << "-yielding-" << endl;cin.get();}

	return tfactor;
}

void BilinearCyclic::Print (OPS_Stream &s, int flag)
{
	this->CyclicModel::Print (s, flag);
	s << "+BilinearCyclic\n";
	//s << "   taskFactor = " << getTaskFactor() << endln;
	s << "----------------------------------------"
      << "----------------------------------------"
	  << endln;
}



