#ifndef Inelastic2DYS02_H
#define Inelastic2DYS02_H

#include "InelasticYS2DGNL.h"

class CyclicModel;

class Inelastic2DYS02 : public InelasticYS2DGNL
{
public:
  Inelastic2DYS02(int tag, double A, double E, double Iz,
		  int Nd1, int Nd2,
		  YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
		  // int cyc_type, double wt,
		  CyclicModel *cycModel,
		  double del_p_max,
		  double Alpha, double Beta, int rf_algo=-1,
		  bool islinear=false, double rho=0.0);

  ~Inelastic2DYS02();
  int commitState(void);
  int update(void);

 protected:
  void getLocalStiff(Matrix &K);

 private:
  // int cycType;
  // double wT;
  double A;
  double E;
  double Iz;
  double resFactor;
  double alfa, beta;
  double delPmax;
  double delPMaxPos, delPMaxNeg;
  CyclicModel *cModel;
};

#endif

