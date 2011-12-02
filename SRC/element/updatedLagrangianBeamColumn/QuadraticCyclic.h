#ifndef QuadraticCyclic_H
#define QuadraticCyclic_H

#include <CyclicModel.h>
#include <Vector.h>
#include <Matrix.h>

class QuadraticCyclic : public CyclicModel
{
public:
	QuadraticCyclic(int tag, double wt=0.9, double qy=0.33);
	~QuadraticCyclic();
	 void Print (OPS_Stream &, int = 0);
	 CyclicModel *getCopy();

protected:
	 int createFullCycleTask();
	 int createHalfCycleTask();
	double getTaskFactor();

private:
	int solveQuad(double x1, double y1, double x2,
                  double y2, double x3, double y);
	double getQuadFactor(double x1, double y1, double dx);
	int createTask();

private:
	double a, b, c;
	double weightFactor, facty;
	double qx1, qy1, qx2, qy2, qx3, qy3;
	
	static Matrix X;
	static Vector Y, A;
};

#endif

