// MultiLinear.h: interface for the MultiLinear class.
//
//////////////////////////////////////////////////////////////////////


#ifndef ExponReducing_h
#define ExponReducing_h

#include "PlasticHardeningMaterial.h"
#include <math.h>

class ExponReducing : public PlasticHardeningMaterial
{
public:
	ExponReducing(int tag, double kp0, double alfa);
	ExponReducing(int tag, double kp0, double alfa, double res_fact);
	virtual ~ExponReducing();
	
	double getTrialPlasticStiffness();
    PlasticHardeningMaterial *getCopy(void);
    void Print(OPS_Stream &s, int flag =0);

  private:
  double Kp0;
  double alpha;
  double resFactor;
};

#endif
