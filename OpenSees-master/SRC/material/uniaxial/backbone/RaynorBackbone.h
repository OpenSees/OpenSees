

#ifndef RaynorBackbone_h
#define RaynorBackbone_h

#include <HystereticBackbone.h>
#include <Vector.h>

class RaynorBackbone : public HystereticBackbone
{
 public:
  RaynorBackbone(int tag,double es,double f1,double f2,double epsh,double epsm,double c1,double ey);
  RaynorBackbone();
  ~RaynorBackbone();
  
  double getStress(double strain);
  double getTangent(double strain);
  double getEnergy(double strain);
  
  double getYieldStrain(void);
  
  HystereticBackbone *getCopy(void);
  
  void Print(OPS_Stream &s, int flag = 0);
  
  int setVariable(char *argv);
  int getVariable(int varID, double &theValue);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
 protected:
  
 private:
  double Es;
  double fy;
  double fsu;
  double Epsilonsh;
  double Epsilonsm;
  double C1;
  double Ey;
};

#endif
