#ifndef Inelastic2DYS03_H
#define Inelastic2DYS03_H

#include "InelasticYS2DGNL.h"

class Inelastic2DYS03 : public InelasticYS2DGNL
{
 public:
  Inelastic2DYS03(int tag, double a_ten, double a_com, double e,
		  double iz_pos, double iz_neg, int Nd1, int Nd2,
		  YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
		  int rf_algo=-1, bool islinear=false, double rho=0);
  virtual int commitState();
  
  ~Inelastic2DYS03();
  
 protected:
  virtual void getLocalStiff(Matrix &K);
  
 private:
  double Atens, Acomp, E;
  double IzPos, IzNeg;
  Vector ndisp, ndisp_hist;
};

#endif

