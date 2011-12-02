

#include <RaynorBackbone.h>
#include <Vector.h>
#include <Channel.h>

#include <math.h>

RaynorBackbone::RaynorBackbone(int tag,double es,double f1,double f2,double epsh,double epsm,double c1,double ey):
  HystereticBackbone(tag,BACKBONE_TAG_Raynor),
  Es(es), fy(f1), fsu(f2),Epsilonsh(epsh),Epsilonsm(epsm), C1(c1), Ey(ey)
{
  if (Epsilonsm <= Epsilonsh)
    opserr << "RaynorBackbone::RaynorBackbone -- Esilonsm-Epsilonsh <= 0" << endln;

  if (fy/Es > Epsilonsh) 
    opserr << "RaynorBackbone::RaynorBackbone -- Esilony > Epsilonsh" << endln;
  
}

RaynorBackbone::RaynorBackbone():
  HystereticBackbone(0,BACKBONE_TAG_Raynor),
  Es(0.0), fy(0.0), fsu(0.0),Epsilonsh(0.0),Epsilonsm(0.0), C1(0.0), Ey(0.0)
{
}

RaynorBackbone::~RaynorBackbone()
{
  
}

double
RaynorBackbone::getTangent (double strain)
{
  double Epsilony = fy/Es;
  double fsh = fy + (Epsilonsh - Epsilony)*Ey;
  if( fabs(strain) <= Epsilony ) 
    return Es;
  else if( fabs(strain)<=Epsilonsh && fabs(strain)>Epsilony )
    return Ey;
  else if( strain<=Epsilonsm && strain>Epsilonsh )
  {
    if( C1 != 0.0 )    return (fsu-fsh)/(Epsilonsm-Epsilonsh)*C1*pow((Epsilonsm-strain)/(Epsilonsm-Epsilonsh),C1-1.0);
    else               return 0.0;
  }
  else if( strain>=-Epsilonsm && strain<-Epsilonsh )
  {
    if( C1 != 0.0 )    return (fsu-fsh)/(Epsilonsm-Epsilonsh)*C1*pow((Epsilonsm+strain)/(Epsilonsm-Epsilonsh),C1-1.0);
    else               return 0.0;
  }
  else
  {
    return 0.0001*Es;
  }
}

double
RaynorBackbone::getStress (double strain)
{
  double Epsilony = fy/Es;
  double fsh = fy + (Epsilonsh - Epsilony)*Ey;
  if( fabs(strain) <= Epsilony )
    return Es*strain;
  else if( strain> Epsilony && strain<=Epsilonsh )
    return fy+(strain-Epsilony)*Ey;
  else if( strain<-Epsilony && strain>=-Epsilonsh )
    return -fy+(strain+Epsilony)*Ey;
  else if( strain>Epsilonsh && strain<=Epsilonsm )
  {
     return fsu-(fsu-fsh)*pow((Epsilonsm-strain)/(Epsilonsm-Epsilonsh),C1);
  }
  else if( strain<-Epsilonsh && strain>=-Epsilonsm )
  {
    return -fsu+(fsu-fsh)*pow((Epsilonsm+strain)/(Epsilonsm-Epsilonsh),C1);
  }
  else if (strain < -Epsilonsm)
    return -fsu;
  else if (strain > Epsilonsm)
    return fsu;
}

double
RaynorBackbone::getEnergy (double strain)
{
  double Epsilony = fy/Es;
  double fsh = fy + (Epsilonsh - Epsilony)*Ey;
  if( fabs(strain) <= Epsilony ) 
    return Es*pow(strain,2)/2.0;

  else if( fabs(strain)<=Epsilonsh && fabs(strain)>Epsilony )
    return Es*pow(Epsilony,2)/2.0+fy*(strain-Epsilony)+Ey*pow((strain-Epsilony),2)/2.0;

  else if( fabs(strain)<=Epsilonsm && fabs(strain)>Epsilonsh )
  {
    if( C1 != 0.0 )    return Es*pow(Epsilony,2)/2.0+fy*(Epsilonsh-Epsilony)+Ey*pow((Epsilonsh-Epsilony),2)/2.0
                              +fsu*(strain-Epsilonsh)
                              +(fsu-fsh)/(Epsilonsm-Epsilonsh)/(C1+1.0)
                                *(pow((Epsilonsm-strain)/(Epsilonsm-Epsilonsh),C1+1)-1.0);
    else               return Es*pow(Epsilony,2)/2.0+fy*(Epsilonsh-Epsilony)+Ey*pow((Epsilonsh-Epsilony),2)/2.0
                              +fsh*(strain-Epsilonsh);
  }
  else
  {
    opserr << "RaynorBackbone::RaynorBackbone -- fabs(strain) > Epsilonsm" << endln;
    return 0.0;
  }
}

double
RaynorBackbone::getYieldStrain(void)
{
  return fy/Es;
}

HystereticBackbone*
RaynorBackbone::getCopy(void)
{
  RaynorBackbone *theCopy =
    new RaynorBackbone (this->getTag(), Es, fy, fsu, Epsilonsh, Epsilonsm, C1, Ey);
  
  return theCopy;
}

void
RaynorBackbone::Print(OPS_Stream &s, int flag)
{
  s << "RaynorBackbone, tag: " << this->getTag() << endln;
  s << "\tEs: " << Es << endln;
  s << "\tfy: " << fy << endln;
  s << "\tfsu: " << fsu << endln;
  s << "\tEpsilonsh: " << Epsilonsh << endln;
  s << "\tEpsilonsm: " << Epsilonsm << endln;
  s << "\tC1: " << C1 << endln;
  s << "\tEy: " << Ey << endln;
}

int
RaynorBackbone::setVariable (char *argv)
{
  return -1;
}

int
RaynorBackbone::getVariable (int varID, double &theValue)
{
  return -1;
}

int
RaynorBackbone::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(8);
  
  data(0) = this->getTag();
  data(1) = Es;
  data(2) = fy;
  data(3) = fsu;
  data(4) = Epsilonsh;
  data(5) = Epsilonsm;
  data(6) = C1;
  data(7) = Ey;
  
  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "RaynorBackbone::sendSelf -- could not send Vector" << endln;

    return res;
  }
  
  return res;
}

int
RaynorBackbone::recvSelf(int commitTag, Channel &theChannel, 
			     FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(8);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "RaynorBackbone::recvSelf -- could not receive Vector" << endln;

    return res;
  }
  
  this->setTag(int(data(0)));
  Es = data(1);
  fy = data(2);
  fsu = data(3);
  Epsilonsh = data(4);
  Epsilonsm = data(5);
  C1 = data(6);
  Ey = data(7);
  return res;
}
