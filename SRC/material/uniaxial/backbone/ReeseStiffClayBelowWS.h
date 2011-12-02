

#ifndef ReeseStiffClayBelowWS_h
#define ReeseStiffClayBelowWS_h

#include <HystereticBackbone.h>
#include <Vector.h>

class ReeseStiffClayBelowWS : public HystereticBackbone
{
 public:
  ReeseStiffClayBelowWS(int tag,double esi,double y,double as,double pc);
  ReeseStiffClayBelowWS();
  ~ReeseStiffClayBelowWS();
  
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
  double Esi;
  double y50;
  double As;
  double Pc;
};

#endif
