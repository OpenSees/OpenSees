#ifndef ExplicitDifferenceStatic_h
#define ExplicitDifferenceStatic_h


#include<TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

#define ED_VSIGN_EPS 1e-4

class ExplicitDifferenceStatic : public TransientIntegrator
{
public:
	ExplicitDifferenceStatic();
	ExplicitDifferenceStatic(double alphaM, double betaK, double betaKi, double betaKc);
	~ExplicitDifferenceStatic();                                                                //constructors and unconstructor

	                                                 

	int formEleTangent(FE_Element *theEle);

	int formNodTangent(DOF_Group *theDof);
	
	const Vector & getVel(void);    //added for Modal damping

	int domainChanged(void);
	int newStep(double deltaT);
	int update(const Vector &U);

	int commit(void);

	int formNodalUnbalance(void);

	virtual int sendSelf(int commitTag, Channel &theChannel);
	virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);


protected:
private:
	double deltaT;
	static double deltaT1;
	double alphaM;
	double betaK;
	double betaKi;
	double betaKc;

	int updateCount;
	double c2, c3;
	Vector *U, *Ut;
	Vector  *Utdotdot, *Utdotdot1;
	Vector *Udot, *Utdot, *Utdot1;

    Vector *velSignMem;   // size = number of global equations; entries in {-1,0,+1}
    Vector *prevUnbal;    // Previous step's unbalance size = number of global equations; entries in {-1,0,+1}
    double vSignEps;     // deadband for |v| below which we keep the last sign

};

#endif
