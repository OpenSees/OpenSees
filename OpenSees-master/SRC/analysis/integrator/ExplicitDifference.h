
#ifndef ExplicitDifference_h
#define ExplicitDifference_h


#include<TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class ExplicitDifference : public TransientIntegrator
{
public:
	ExplicitDifference();
	ExplicitDifference(double alphaM, double betaK, double betaKi, double betaKc);
	~ExplicitDifference();                                                                //constructors and unconstructor

	                                                 

	int formEleTangent(FE_Element *theEle);

	int formNodTangent(DOF_Group *theDof);
	
	const Vector & getVel(void);    //added for Modal damping

	int domainChanged(void);
	int newStep(double deltaT);
	int update(const Vector &U);

	int commit(void);

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

};

#endif
