
                                                                       
#include <elementAPI.h>
#include "ConcreteD.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


static int numConcreteD = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_ConcreteD)
{
  // print out some KUDO's
  if (numConcreteD == 0) {
   // opserr << "ConcreteD unaxial material - Written by Zenyong Wan Tongji University 2014.01\n";
    numConcreteD =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  //
  // parse the input line for the material parameters
  //
  
  int    iData[1];
  double dData[9];
  
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid ConcreteD tag" << endln;
    return 0;
  }
  
  numData = OPS_GetNumRemainingInputArgs();
  
  if(numData != 7 && numData != 9 )
    {
      opserr << "Invalid #args, want: uniaxialMaterial ConcreteD "<< iData[0] <<"(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)" << endln;
      return 0;
    }
  
  if(numData==7)  {
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args: uniaxialMaterial ConcreteD "<< iData[0] <<"(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)" << endln;
      return 0;
    } 
    // Parsing was successful, allocate the material
    theMaterial = new ConcreteD(iData[0], dData[0], dData[1],dData[2],
				dData[3],dData[4],dData[5],dData[6]);  
    
  }else if (numData==9){
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid #args: uniaxialMaterial ConcreteD "<< iData[0] <<"(fcr? epcr? ft? eptr? Ec? alphac? alphat? <cesp? etap?>)" << endln;
      return 0;
    } 
    // Parsing was successful, allocate the material
    theMaterial = new ConcreteD(iData[0], dData[0], dData[1],dData[2],
				dData[3],dData[4],dData[5],dData[6],dData[7],dData[8]); 
  }
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ConcreteD\n";
    return 0;
  }
  
  // return the material
  return theMaterial;
}



//Constructor
ConcreteD::ConcreteD(int tag, double fc0, double ec0,double ft0,double eptt0, double Ec0,double alphac0 ,double alphat0,double cesp0,double etap0) 
:UniaxialMaterial(tag,0), 
fcc(fc0),epcc(ec0),ft(ft0),Ec(Ec0),alphac(alphac0),alphat(alphat0),TDtp(0.0),
cesp(cesp0),etap(etap0), CLoadState(0),CStress(0.0),CStrain(0.0),
CTangent(Ec0),CDc(0.0),CDt(0.0),CDcp(0.0),CEpp(0.0),CRc(0.0),eptt(eptt0),
CRt(0.0),TLoadState(0),TStrain(0.0),TStress(0.0),TTangent(0.0),TDc(0.0),
TDt(0.0),TDcp(0.0),TEpp(0.0),TRc(0.0),TRt(0.0),CSecant(Ec0),TSecant(Ec0),CDtp(0.0)

{
}
ConcreteD::ConcreteD(int tag, double fc0, double ec0,double ft0,double eptt0, double Ec0,double alphac0 ,double alphat0) 
:UniaxialMaterial(tag,0), 
fcc(fc0),epcc(ec0),ft(ft0),Ec(Ec0),alphac(alphac0),alphat(alphat0),TDtp(0.0),
cesp(0.25),etap(1.15), CLoadState(0),CStress(0.0),CStrain(0.0),
CTangent(Ec0),CDc(0.0),CDt(0.0),CDcp(0.0),CEpp(0.0),CRc(0.0),eptt(eptt0),
CRt(0.0),TLoadState(0),TStrain(0.0),TStress(0.0),TTangent(0.0),TDc(0.0),
TDt(0.0),TDcp(0.0),TEpp(0.0),TRc(0.0),TRt(0.0),CSecant(Ec0),TSecant(Ec0),CDtp(0.0)

{
}

//Constructor
ConcreteD::ConcreteD()
:UniaxialMaterial(0, 0),
fcc(0.0),epcc(0.0),ft(0.0),Ec(0.0),alphac(0.0),alphat(0.0),TDtp(0.0),
cesp(0.0),etap(0.0), CLoadState(0),CStress(0.0),CStrain(0.0),CDtp(0.0),
CTangent(0.0),CDc(0.0),CDt(0.0),CDcp(0.0),CEpp(0.0),CRc(0.0),eptt(0.0),
CRt(0.0),TLoadState(0),TStrain(0.0),TStress(0.0),TTangent(0.0),TDc(0.0),
TDt(0.0),TDcp(0.0),TEpp(0.0),TRc(0.0),TRt(0.0),CSecant(0.0),TSecant(0.0)
{
}

//Destructor
ConcreteD::~ConcreteD()
{
}

void 
ConcreteD::envelope()
{if(TStrain>=TEpp)//Tension envelope
	{
	double xt,nt,rowt;
	TRt=(TStrain-TEpp);
	xt=TRt/eptt;
	rowt=ft/(Ec*eptt);
	nt=1/(1-rowt);
	double pDcpEp,pDtppDt;
	if(xt<1.0)
	{TDt=1-rowt*nt/(nt-1+pow(xt,nt));
	TDtp=TDt+TDc-TDt*TDc;
	pDcpEp=(nt*nt*rowt*pow(xt,(nt - 1)))/pow((nt + pow(xt,nt) - 1),2)/eptt;// partial dc partial epc
	pDtppDt=1-TDc;
	}
	else
	{TDt=1-rowt/(alphat*(xt-1)*(xt-1)+xt);
	TDtp=TDt+TDc-TDt*TDc;
	pDcpEp=(rowt*(alphat*(2*xt - 2) + 1))/pow(xt + alphat*pow(xt - 1,2),2)/eptt;}
	pDtppDt=1-TDc;
	TStress=(1-TDtp)*Ec*(TStrain-TEpp);
	TTangent=Ec*(1-TDt+pDcpEp*pDtppDt*(TEpp-TStrain));
	}
else              //compression envelope
	{
	double xc,nc,rowc;
	TRc=TStrain;
	xc=TRc/epcc;
	rowc=fcc/(Ec*epcc);
	nc=1/(1-rowc);
	double pDcpEp;
	if(xc<1.0)
	{TDc=1-rowc*nc/(nc-1+pow(xc,nc));
	TDtp=TDt+TDc-TDt*TDc;
	pDcpEp=(nc*nc*rowc*pow(xc,(nc - 1)))/pow((nc + pow(xc,nc) - 1),2)/epcc;// partial dc partial epc
	}
	else
	{TDc=1-rowc/(alphac*(xc-1)*(xc-1)+xc);
	TDtp=TDt+TDc-TDt*TDc;
	pDcpEp=(rowc*(alphac*(2*xc - 2) + 1))/pow(xc + alphac*pow(xc - 1,2),2)/epcc;}
	double Fip;
	//Fip=thetap*pow(TDc,etap)/(1+thetap*pow(TDc,etap));
	Fip=cesp*(pow(2.718,TDc*etap)-1);
	//if(Fip>TDc)
	//{Fip=TDc;}
	TEpp=Fip*TRc;
	TDcp=(TDc-Fip)/(1-Fip);
	//if(TDcp<0.0)
	//{TDcp=0.0;}
	TStress=(1-TDcp)*Ec*(TStrain-TEpp);
	TTangent=Ec*(1-TDc-pDcpEp*TStrain);
	}
}

void
ConcreteD::unload()  //unloading, reloading curve
{if(TStrain>TEpp)    //
{
	TStress=(1-TDtp)*Ec*(TStrain-TEpp);
	TTangent=(1-TDtp)*Ec;
}
else
{		
	TStress=(1-TDcp)*Ec*(TStrain-TEpp);
	TTangent=(1-TDcp)*Ec;
}
}

int 
ConcreteD::setTrialStrain(double strain, double strainRate)
{	

	TLoadState	=	CLoadState;
	TStress		=	CStress;
	TTangent	=	CTangent;
	TDc			=	CDc;
	TDt			=	CDt;
	TDcp		=	CDcp;
	TDtp		=	CDtp;
	TEpp		=	CEpp;
	TRc 		=	CRc;
	TRt	    	=	CRt;

	TStrain		=	strain;
	double dStrain = strain-CStrain;

    if (fabs(dStrain) < DBL_EPSILON)
	{return 0;}

if(TStrain>TEpp)					//Tension
{
		if(TLoadState==0)		    //On the envelope
		{
		if(dStrain>0.0)			    //Tension
		{	envelope();
			return 0;	}
		else
		{	TLoadState=1;		    //unloading, reloading
			unload();
			return 0;	}
		}else 	{	
			if(TStrain-TEpp<TRt)
			{unload();
			return 0;}
			else
			{TLoadState=0;		    //On the envelope
			envelope();	
			return 0;}
			}

}
//��ѹ
else
	{
		if(TLoadState==0)		   //On the envelope
		{
		if(dStrain<0.0)			   
		{	envelope();
			return 0;	}
		else
		{	TLoadState=1;		   //unloading, reloading
			unload();
			return 0;	}
		}else 	{	
			if(TStrain>TRc)		
			{unload();
			return 0;}
			else
			{TLoadState=0;		   //On the envelope
			envelope();	
			return 0;}
		}
	}
}


/*
int
ConcreteD::setTrial(double strain,double &stress,double &tangent,double strainRate)
{ 
	int Return_a= setTrialStrain(strain,strainRate);
	return Return_a;
}*/
/*
int 
ConcreteD::determineTrialState(double dStain)
{
	int Return_b= setTrialStrain(CStrain+dStain);
	return Return_b;
}*/
double 
ConcreteD::getStrain(void)
{
  return TStrain;
}

double 
ConcreteD::getStress(void)
{
  return TStress;
}

double 
ConcreteD::getSecant(void)
{
  if (TStrain==0)
	  return Ec;
  else
	  return TStress/TStrain;

	//return TDc;
}

double 
ConcreteD::getTangent(void)
{
  return TTangent;
}

int 
ConcreteD::commitState(void)
{
    CLoadState	=	TLoadState;
	CStrain		=	TStrain;
	CStress		=	TStress;
	CTangent	=	TTangent;
	CDc			=	TDc;
	CDt			=	TDt;
	CDcp		=	TDcp;
	CDtp		=	TDtp;
	CEpp		=	TEpp;
	CRc			=	TRc;
	CRt			=	TRt;

    return 0;
}	

int 
ConcreteD::revertToLastCommit(void)
{
	TLoadState	=	CLoadState;
	TStress		=	CStress;
	TTangent	=	CTangent;
	TDc			=	CDc;
	TDt			=	CDt;
	TDcp		=	CDcp;
	TDtp		=	CDtp;
	TEpp		=	CEpp;
	TRc 		=	CRc;
	TRt	    	=	CRt;
	
  return 0;
}


int 
ConcreteD::revertToStart(void)
{
	CLoadState	= 0;
	CStrain		= 0.0;
	CStress		= 0.0;
	CTangent	= Ec;
	CDc			= 0.0;
	CDt			= 0.0;
	CDcp		= 0.0;
	CDtp		= 0.0;
	CEpp		= 0.0;
	CRc			= 0.0;
	CRt			= 0.0;

  return 0;
}


UniaxialMaterial *
ConcreteD::getCopy(void)
{
  ConcreteD *theCopy =
    new ConcreteD(this->getTag(),fcc,epcc,ft,eptt,Ec,alphac,alphat,cesp,etap);
//History Parameters
  theCopy->CLoadState = this->CLoadState;
  theCopy->CDc        = this->CDc;
  theCopy->CDt        = this->CDt;
  theCopy->CDcp       = this->CDcp;  
  theCopy->CEpp       = this->CEpp;
  theCopy->CRc        = this->CRc;
  theCopy->CRt        = this->CRt;
  theCopy->CDtp		  = this->CDtp;
//State Parameters
  theCopy->CStress    = this->CStress;
  theCopy->CStrain    = this->CStrain;
  theCopy->CTangent   = this->CTangent;
  return theCopy;
}


int 
ConcreteD::sendSelf(int cTag, Channel &theChannel)
{
int res		= 0;
static		Vector data(12);
data(0)		=   this->getTag();
data(1)		=	CLoadState;
data(2)		=	CDc;
data(3)		=	CDt;
data(4)		=	CDcp;
data(5)		=	CEpp;
data(6)		=	CRc;
data(7)		=	CRt;
data(8)		=	CStress;
data(9)		=	CStrain;
data(10)	=	CTangent;
data(11)	=	CDtp;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ConcreteD::sendSelf() - failed to send data\n";

  return res;
}

int 
ConcreteD::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(12);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ConcreteD::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
	CLoadState	=	data(1);
	CDc			=	data(2);
	CDt			=	data(3);
	CDcp		=	data(4);
	CEpp		=	data(5);
	CRc			=	data(6);
	CRt			=	data(7);
	CStress		=	data(8);
	CStrain		=	data(9);
	CTangent	=	data(10);
	CDtp		=	data(11);
	TStrain		=	CStrain;
	TStress		=	CStress;
	TTangent	=	CTangent;

  }

  return res;
}

void 
ConcreteD::Print(OPS_Stream &s, int flag)
{
  s << "ConcreteD tag: " << this->getTag() << endln;
 // s << "  Ec: " << Ec << endln;
 // s << "  epc0: " << epc0<< endln;
 // s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
}


