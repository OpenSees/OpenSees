#ifndef BilinearCyclic_H
#define BilinearCyclic_H

#include <CyclicModel.h>

class BilinearCyclic : public CyclicModel
{
public:
	BilinearCyclic(int tag, double weight=0.9);
	~BilinearCyclic();
	 void Print (OPS_Stream &, int = 0);
	 CyclicModel *getCopy();

protected:
	double getTaskFactor();
	double weightFactor;
};

#endif

