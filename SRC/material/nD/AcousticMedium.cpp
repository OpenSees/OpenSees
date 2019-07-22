
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// Acoustic fulid Material
//------------------------------------------

#include <string.h>

#include <AcousticMedium.h>

#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <FEM_ObjectBroker.h>  

#include <OPS_Globals.h>
#include <elementAPI.h>
#include <MaterialResponse.h>

Matrix AcousticMedium::D(1,1);	  // global for AcousticMedium only
Vector AcousticMedium::sigma(3);	// global for AcousticMedium only
Matrix AcousticMedium::DSensitivity(1,1);	  // global for AcousticMedium only

void *
OPS_AcousticMedium(void)
{
  NDMaterial *theMaterial = 0;
  
  int numArgs = OPS_GetNumRemainingInputArgs();
  
  if (numArgs < 3) {
    printf("Want: nDMaterial AcousticMedium $tag $K $rho <$gamma>\n");
    return 0;	
  }
  
  int iData[1];
  double dData[3];
  dData[2] = 0.0;
  
  int numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    printf("WARNING invalid integer tag: nDMaterial AcousticMedium \n");
    return 0;
  }
  
  if (numArgs > 3) 
    numData = 3;
  else
    numData = 2;
  
  if (OPS_GetDouble(&numData, dData) != 0) {
    printf("WARNING invalid data: nDMaterial AcousticMedium : %d\n", iData[0]);
    return 0;
  }

  // print input
  //  printf("The material tag is %d!\n", iData[0]);
  //  printf("Kf is %g!\n", dData[0]);
  //  printf("Rho is %g!\n", dData[1]);
  //  printf("Gamma is %g!\n", dData[2]);

  theMaterial = new AcousticMedium(iData[0], dData[0], dData[1], dData[2]);
  
  return theMaterial;
}

AcousticMedium::AcousticMedium
(int tag, int classTag, double k, double r, double gamma)
  :NDMaterial(tag, classTag), Kf(k), rho(r), Gamma(gamma)
{

  //============Sensitivity 
    SHVs =0;
    parameterID =0;

}

AcousticMedium::AcousticMedium
(int tag, double k, double r, double gamma)
  :NDMaterial(tag, ND_TAG_AcousticMedium), Kf(k), rho(r), Gamma(gamma)
{
  //============Sensitivity 
    SHVs =0;
   parameterID =0; 
}

AcousticMedium::~AcousticMedium()
{
	
}

double
AcousticMedium::getRho() 
{ 
  return rho ;
}

NDMaterial*
AcousticMedium::getCopy (const char *type)
{
  AcousticMedium *theModel = 
    new AcousticMedium(this->getTag(), Kf, rho, Gamma);
  
  return theModel;
}

int
AcousticMedium::setTrialStrain (const Vector &v)
{
  epsilon = v;
  
  return 0;
}

int
AcousticMedium::setTrialStrain (const Vector &v, const Vector &rate)
{
  epsilon = v;
  
  return 0;
}

int
AcousticMedium::setTrialStrainIncr (const Vector &v)
{
  epsilon = epsilon + v;
  
  return 0;
}

int
AcousticMedium::setTrialStrainIncr (const Vector &v, const Vector &rate)
{
  epsilon = epsilon + v;
  
  return 0;
}

const Matrix&
AcousticMedium::getTangent (void)
{
  D(0,0) = Kf;

  return D;
}

const Matrix&
AcousticMedium::getDampTangentSensitivity(int gradNumber) 
{
	DSensitivity.Zero();

	if (parameterID == 1) {  

    DSensitivity(0,0) = 1;

	}
	return DSensitivity;
}

double 
AcousticMedium::getRhoSensitivity(int gradNumber) 
{ 
   double rhoSensitivity = 0 ;

   if (parameterID == 2) {  

     rhoSensitivity = 1 ;
   }
   return rhoSensitivity;
}

const Matrix&
AcousticMedium::getInitialTangent (void)
{
  return this->getTangent();
}

const Vector& 
AcousticMedium::getStress (void)
{
  sigma = epsilon/rho;
  

 //opserr<< "sigma is "<<sigma<<endln;

  return sigma;
}

const Vector&
AcousticMedium::getStrain (void)
{
  return epsilon;
}

int
AcousticMedium::commitState (void)
{
  return 0;
}

int
AcousticMedium::revertToLastCommit (void)
{
  return 0;
}

int
AcousticMedium::revertToStart (void)
{
  return 0;
}

NDMaterial*
AcousticMedium::getCopy (void)
{
  AcousticMedium *theModel = 
    new AcousticMedium(this->getTag(), Kf, rho, Gamma);
  
  return theModel;
}

const char*
AcousticMedium::getType (void) const
{
  return "AcousticMedium";
}

int
AcousticMedium::getOrder (void) const
{
  return 1;
}

int
AcousticMedium::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(4);
  
  data(0) = this->getTag();
  data(1) = Kf;
  data(2) = rho;
  data(3) = Gamma;
  
 res += theChannel.sendVector(this->getDbTag(), commitTag, data);
 if (res < 0) {
   opserr << "AcousticMedium::sendSelf -- could not send Vector\n";
   return res;
 }

 return res;
}





int
AcousticMedium::recvSelf (int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(4);
  
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
   opserr << "AcousticMedium::recvSelf -- could not recv Vector\n";
   return res;
  }
    
  this->setTag((int)data(0));
  Kf = data(1);
  rho = data(2);
  Gamma = data(3);
  
  return res;
}




Response * AcousticMedium::setResponse (const char **argv, int argc, OPS_Stream &matInformation)  {

  if (strcmp(argv[0],"sigma") == 0 )
		return new MaterialResponse(this, 1, sigma);

  else if (strcmp(argv[0],"epsilon") == 0)
		return new MaterialResponse(this, 2, epsilon);


		return 0;

};

int AcousticMedium::getResponse (int responseID, Information &matInfo)  {
	
		switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = sigma;
			return 0;

		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = epsilon;
			return 0;


		}

	return 0;

};




void
AcousticMedium::Print (OPS_Stream &s, int flag)
{
	s << "Elastic Isotropic Material Model" << endln;
	s << "\tKf:  " << Kf << endln;
	s << "\trho:  " << rho << endln;
	s << "\tGamma:  " << Gamma << endln;

	return;
}

///////////////////////// add sensitivity ///////////////////////////////////

int
AcousticMedium::setParameter(const char **argv, int argc,
				      Parameter &param)
{
  if (strcmp(argv[0],"Kf") == 0)
    return param.addObject(1, this);
  
  else if (strcmp(argv[0],"rho") == 0)
    return param.addObject(2, this);
  
  else if (strcmp(argv[0],"gamma") == 0)
    return param.addObject(3, this);

  return -1;
}

int 
AcousticMedium::updateParameter(int parameterID, Information &info)
{ 
  switch(parameterID) {
  case 1:
    Kf = info.theDouble;
    return 0;
  case 2:
    rho = info.theDouble;
    return 0;
  case 3:
    Gamma = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int AcousticMedium::activateParameter(int passedParameterID){
	parameterID = passedParameterID;
	return 0;						
	}

const Vector &  AcousticMedium::getStressSensitivity (int gradNumber, bool conditional){

	double KfSensitivity = 0.0;
	double rhoSensitivity = 0.0;
	double GammaSensitivity = 0.0;  

	if (parameterID == 1) {  
		KfSensitivity = 1.0;
	}
	else if (parameterID == 2) {   
		rhoSensitivity = 1.0;
	}
	else if (parameterID == 3) {   
		GammaSensitivity = 1.0;
	}
	else {

		// Nothing random here, but may have to return something in any case

	}

  // --- define history variables -----------
	static Vector CepsilonSensitivity(3);	CepsilonSensitivity.Zero();
	static Vector CsigmaSensitivity(3);   	CsigmaSensitivity.Zero();

	static Vector sigmaSensitivity(3);
    sigmaSensitivity.Zero();

	Vector epsilonSensitivity(3);   epsilonSensitivity.Zero();   // conditional sensitivity

	sigmaSensitivity.addVector(0.0, epsilonSensitivity, rho); 
	sigmaSensitivity.addVector(1.0, epsilon, -1.0*rhoSensitivity); 
	sigmaSensitivity.addVector(0.0, sigmaSensitivity, 1.0/rho/rho); 

	//opserr<<"sigmaSensitivity in water is "<<sigmaSensitivity<<endln;
	return sigmaSensitivity;

}

int AcousticMedium::commitSensitivity(const Vector & strainGradient, int gradNumber, int numGrads){
    
	double KfSensitivity = 0.0;
	double rhoSensitivity = 0.0;
	double GammaSensitivity = 0.0;  
 
	static Vector epsilonSensitivity(3);	epsilonSensitivity.Zero();
	static Vector sigmaSensitivity(3);   	sigmaSensitivity.Zero();

		epsilonSensitivity[0] = strainGradient[0];
		epsilonSensitivity[1] = strainGradient[1];
		epsilonSensitivity[2] = strainGradient[2];


if (parameterID == 1) {  
		KfSensitivity = 1.0;
	}
	else if (parameterID == 2) {   
		rhoSensitivity = 1.0;
	}
	else if (parameterID == 3) {   
		GammaSensitivity = 1.0;
	}
	else {

		// Nothing random here, but may have to return something in any case

	}
  // --- define history variables -----------
	static Vector CepsilonSensitivity(3);	CepsilonSensitivity.Zero();
	static Vector CsigmaSensitivity(3);   	CsigmaSensitivity.Zero();


	sigmaSensitivity.addVector(0.0, epsilonSensitivity, rho); 
	sigmaSensitivity.addVector(1.0, epsilon, -1.0*rhoSensitivity); 
	sigmaSensitivity.addVector(0.0, sigmaSensitivity, 1.0/rho/rho); 

	return 0;

}; 
