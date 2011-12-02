#include "CyclicModel.h"
#include <math.h>

const int    CyclicModel::Loading(1);
const int    CyclicModel::Unloading(2);
const int    CyclicModel::Crossover(3);
const double CyclicModel::Tol(1e-10);
const double CyclicModel::delK(0.85);

CyclicModel::CyclicModel(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag),
 resFactor(1.0),
 cycFactor(1.0), cycFactor_hist(1.0),
 f_hist(0.0), d_hist(0.0), f_curr(0.0), d_curr(0.0),
 f_bgn(0.0), d_bgn(0.0), f_end(0.0), d_end(0.0),
 initFpos(0.0),initDpos(0.0),initFneg(0.0),initDneg(0.0),
 initFMag(0.0),initDMag(0.0),
 k_init(0.0),
 delT_curr(0.0), delT_hist(0.0),
 k_hist(0.0), k_curr(0.0),
 initYieldPos(false), initYieldNeg(false), 
 initCyc(false), yielding(false), yielding_hist(false),
 fpeakPos(0.0), dpeakPos(0.0),fpeakNeg(0.0),dpeakNeg(0.0),
 state_hist(Loading), state_curr(Loading)
{

}


CyclicModel::~CyclicModel()
{
  // does nothing
}


double CyclicModel::getFactor()
{
	if(cycFactor < 0.05)
		cycFactor = 0.05;

	if(state_curr == Unloading && state_hist == Loading)
		cycFactor = resFactor;
		
	return cycFactor;
}

void CyclicModel::update(double f, double d, bool yield)
{
	if( yield && (!initYieldPos &&  !initYieldNeg))
	{
		initDMag = d;
		initFMag = f;
	}
	if( !yield && (!initYieldPos &&  !initYieldNeg))
	{
		initDMag = d;
		initFMag = f;
	}

	// non-dimentionalize d_curr and f_curr
	d_curr = d/initDMag;
	f_curr = f/initFMag;

	yielding = yield;

// check state
	if(fabs(f_curr) < fabs(f_hist) &&
	   fabs(d_curr) < fabs(d_hist) &&
	   dir(f_curr)== dir(f_hist)
	  )
		state_curr = Unloading;
	else
		state_curr = Loading;

	if(f_curr*f_hist < 0)
		state_curr = Crossover;
}

int CyclicModel::commitState(double newRes)
{
//	d_curr = d;
//	f_curr = f;

// first make sure we don't encounter any divide by 0 errors
	if (fabs(d_curr - d_hist) > Tol)
		k_curr = fabs((f_curr - f_hist)/(d_curr - d_hist));

	// check if the yield-point has been found
	if(!initYieldPos && yielding && d_curr > 0.0)
	{
		initDpos = d_curr;
		initFpos = f_curr;
		initYieldPos = true;
		if(!initYieldNeg)
			k_init = f_curr/d_curr;
	}

	if(!initYieldNeg && yielding && d_curr < 0.0)
	{
		initDneg = d_curr;
		initFneg = f_curr;
		initYieldNeg = true;
		if(!initYieldPos)
			k_init = f_curr/d_curr;
	}


// Temp approx. init values
	if(initYieldPos && !initYieldNeg)
	{
		initDneg = -1*initDpos;
		initFneg = -1*initFpos;
	}
	if(initYieldNeg && !initYieldPos)
	{
		initDpos = -1*initDneg;
		initFpos = -1*initFneg;
	}

// now update the commit values

	if(d_curr > dpeakPos)
	{
		dpeakPos = d_curr;
		fpeakPos = f_curr;
	}

	if(d_curr < 0 && fabs(d_curr) > fabs(dpeakNeg))
	{
		dpeakNeg = d_curr;
		fpeakNeg = f_curr;
	}

	setCurrent(f_curr, d_curr);
	//cout << *this;

	k_hist = k_curr;
	f_hist = f_curr;
	d_hist = d_curr;
	cycFactor_hist = cycFactor;
	resFactor = newRes;

	// opserr << " Resfactor = " << newRes << endln;
	
	state_hist = state_curr;
	yielding_hist = yielding;

	delT_hist = delT_curr;

	return 0;
}

int CyclicModel::setCurrent(double f, double d)
{
	// check if initially yielded
	if(d_curr > 0 && !initYieldPos ||
       d_curr < 0 && !initYieldNeg
	  )
	{
		cycFactor = 1.0; // resFactor; - not yielded yet
		// return cycFactor;
		return 0;
	}

	if (fabs(d - d_hist) < Tol)
	{
		state_curr = Loading;
		// opserr << "CYC1 \n"; -> don't remove
		cycFactor = cycFactor_hist;
		return 0;
	}

	// check if new full-cycle task needs to be created
	if(state_curr == Unloading && state_hist == Loading && f_curr*f_hist > 0)
	{
		if (createFullCycleTask() < 0)
		{
			opserr << "WARNING - CyclicModel::getFactor(), createFullCycleTask failed\n";
			cycFactor = resFactor;
		}
		else
			cycFactor = getTaskFactor();
		// opserr << *this;

		// return cycFactor;
		return 0;
	}

	// else check the status of the current task
	int res = taskStatus();
	if(res < 0) // task aborted, create new (half-cycle)
	{
		opserr << "Task aborted, creating new half-cycle task\n";
		//cout <<"BEFORE\n" << *this;

		if (createHalfCycleTask() < 0)
		{
			opserr << "WARNING - CyclicModel::getFactor(), createHalfCycleTask failed\n";
			cycFactor = resFactor;
		}
		else
			cycFactor = getTaskFactor();
		// opserr <<"AFTER\n" << *this;
		// opserr << "\a";

	}
	else if (res > 0) // existing task continuing
	{
		// opserr << "Task continuing\n";
		cycFactor = getTaskFactor();
	}
	else if (res == 0)
	{
		// opserr << "Task does not exit (or completed)\n";
		// cycFactor = resFactor or cyc_hist?
		cycFactor = cycFactor_hist;
	}
	else
	{
		opserr << "FATAL: CyclicModel::getFactor() - nan \n";
		opserr << "CYC3 \n";
		cycFactor = 1.0;
	}

	if(cycFactor > (1.001))
	{
		// opserr << "WARNING - cyc > Ko\n";
		// opserr << *this;
		// opserr << "\a";
		// opserr << "CYC4 \n";
		cycFactor = 1.0;
	}

	// return cycFactor;
	return 0;
}

int CyclicModel::createHalfCycleTask()
{
	int res =  initNewTask();
	delT_hist = fabs(d_hist - d_end);
	delT_curr = fabs(d_curr - d_end);
	return res;
}

int CyclicModel::createFullCycleTask()
{
	opserr << "Creating new full-cycle task\n";
	// opserr << *this;
    // opserr << "\a";
    
	initCyc = true;
	int res =  initNewTask();
	delT_hist = fabs(d_hist - d_end);
	delT_curr = fabs(d_curr - d_end);
	return res;
}

// start-point should be f_hist, d_hist
// end-point is fpeak (corresponding to dpeak)
// direction is diven by (d_curr - d_hist)
int CyclicModel::initNewTask()
{
	f_bgn = f_hist;
	d_bgn = d_hist;

	double sgn = d_curr - d_hist;

	if(sgn > 0)
	{
		f_end = fpeakPos;
		d_end = dpeakPos;

		if(!initYieldPos)
		{
			f_end = initFpos;
			d_end = initDpos;
		}
	}
	else
	{
		f_end = fpeakNeg;
		d_end = dpeakNeg;

		if(!initYieldNeg)
		{
			f_end = initFneg;
			d_end = initDneg;
		}
	}

	return 0;
}

int CyclicModel::taskStatus()
{
	if(!initCyc)
		return 0; // task does not exist

	delT_curr = fabs(d_curr - d_end);

	if(fabs(d_curr) >= fabs(d_end) && dir(d_curr)==dir(d_end))
	{
		initCyc = false;
		return 0; // task completed
	}

int status = -1; // aborted
	if( delT_curr <= delT_hist)
		status = 1; // task continuing

	return status;
}


double CyclicModel::rationalize(double x1, double y1, double x2, double y2)
{
	double tangent = fabs((y2-y1)/(x2-x1));
	double tang0   = k_init; // = 1.0

	return (tangent/tang0);
}

bool CyclicModel::contains(double x1, double x2, double x)
{
bool res = false;
	if( (x >= x1 && x <= x2) || (x <= x1 && x >= x2))
		res = true;
	return res;
}

int CyclicModel::dir(double x)
{
int sgn;
	(x>0)?sgn=1:sgn=-1;
	return sgn;
}

void CyclicModel::Print (OPS_Stream &s, int flag)
{
	s << "+CyclicModel, Tag: " << getTag() << endln;
	s << "|  f curr  = " << f_curr << ", d curr  = " << d_curr << endln;
	s << "|  f commit = " << f_hist  << ", d commit = " << d_hist  << endln;
	s << "|  state = " << state_curr << endln;
	s << "|  (1: loading, 2:unloading, 3: cross-over)\n";
	s << "|  Yielding: "; (yielding)? s << "TRUE\n" : s << "FALSE" << endln;
	s << "|  " << endln;
	s << "|  d_bgn = " << d_bgn << ", f_bgn = " << f_bgn << endln;
	s << "|  d_end = " << d_end << ", f_end = " << f_end << endln;
	s << "|  " << endln;		
	s << "|  delT curr = " << delT_curr << ", delT_hist = " << delT_hist << endln;
	s << "|  initFpos: " << initFpos << ", initDpos: " << initDpos << endln;
	s << "|  initFneg: " << initFneg << ", initDneg: " << initDneg << endln;
	s << "|  k_init  : " << k_init << endln;
	s << "|  dpeakPos: " << dpeakPos << ", fpeakPos: " << fpeakPos << endln;
	s << "|  dpeakNeg: " << dpeakNeg << ", fpeakNeg: " << fpeakNeg << endln;
	s << "|  " << endln;
	s << "|  resFactor  -> " << resFactor << endln;
	s << "|  realFactor -> " << getFactor() << endln;
}


