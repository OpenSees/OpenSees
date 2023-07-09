#ifndef LinearCyclic_H
#define LinearCyclic_H

#include <CyclicModel.h>

class LinearCyclic : public CyclicModel
{
public:
	LinearCyclic(int tag);
	~LinearCyclic();
	 void Print (OPS_Stream &, int = 0);
	 CyclicModel *getCopy();

protected:
	double getTaskFactor();
};

#endif

