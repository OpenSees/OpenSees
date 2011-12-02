#include "LinearCyclic.h"
#include <math.h>

LinearCyclic::LinearCyclic(int tag)
:CyclicModel(tag, -1)
{

}


LinearCyclic::~LinearCyclic()
{
  // does nothing
}


CyclicModel *LinearCyclic::getCopy()
{
CyclicModel *newModel = new LinearCyclic(getTag());	
	return newModel;
}

double LinearCyclic::getTaskFactor()
{
	return resFactor;
}

void LinearCyclic::Print (OPS_Stream &s, int flag)
{
	this->CyclicModel::Print (s, flag);
	s << "+LinearCyclic\n";
	s << "   taskFactor = " << resFactor << endln;
	s << "----------------------------------------"
      << "----------------------------------------"
	  << endln;
}



