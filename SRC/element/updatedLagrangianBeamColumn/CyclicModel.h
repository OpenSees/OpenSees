#ifndef CyclicModel_H
#define CyclicModel_H

#include <TaggedObject.h>
#include <MovableObject.h>
#include <OPS_Globals.h>

class CyclicModel : public TaggedObject, public MovableObject
{
public:
	CyclicModel(int tag, int classTag);
	~CyclicModel();
	
	int commitState(double newResidual);
    void update(double f, double d, bool yield);
	virtual CyclicModel *getCopy()=0;
	
	virtual double getFactor();
	// for now ignore
	virtual int 	sendSelf(int commitTag, Channel &theChannel){return -1;}
    virtual int 	recvSelf(int commitTag, Channel &theChannel,
							 FEM_ObjectBroker &theBroker){return -1;}
	virtual void Print (OPS_Stream &s, int flag=0);

protected:
// virtual methods
	virtual int createFullCycleTask();
	virtual int createHalfCycleTask();
	virtual double getTaskFactor()=0;

// protected:
	int setCurrent(double f, double d);
	int dir(double x);

	int taskStatus(void);
	double rationalize(double x1, double y1, double x2, double y2);
	bool   contains(double x1, double x2, double x);

private:
	int initNewTask();

protected:
double resFactor;
double cycFactor, cycFactor_hist;
double f_hist, d_hist;
double f_curr, d_curr;
double delT_curr, delT_hist;
double f_bgn, d_bgn;
double f_end, d_end;

bool   initYieldPos, initYieldNeg;
bool   initCyc;
bool   yielding, yielding_hist;
double initFpos,initDpos;
double initFneg,initDneg;
double initFMag,initDMag;
double k_init;
double k_hist, k_curr;
double fpeakPos, fpeakNeg;
double dpeakPos, dpeakNeg;
int    state_hist, state_curr;

const static int Loading,Unloading, Crossover;
const static double Tol, delK;
};

bool OPS_addCyclicModel(CyclicModel *newComponent);
CyclicModel *OPS_getCyclicModel(int tag);
void OPS_clearAllCyclicModel(void);


#endif

