// Written: KK
// Revision: A  

#include <Concrete0.h>
#include <Vector.h>
#include <Channel.h>
#include <Parameter.h>
#include <string.h>

#include <MaterialResponse.h> // KK
#include <math.h> // CSUF
#include <float.h> // CSUF

#include <OPS_Globals.h>
#include <elementAPI.h> // CSUF

// Read input parameters and build the material // CSUF
void* OPS_Concrete0(void)
{

	// Pointer to a uniaxial material that will be returned                       
	UniaxialMaterial* theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	// Parse the script for material parameters
	if (numArgs < 6) {
		opserr << "WARNING insufficient arguments\n";
		//printCommand(argc, argv);
		opserr << "Want: uniaxialMaterial Concrete0 tag? e0c? fpc? e00? <e0t? ft?>" << endln;
		return 0;
	}

	//int iData[1];
	//double dData[9];
	int tag;
	double e0c, fpc, e00, e0t, ft;
	int mon;

	int numData = 1;
	if (OPS_GetIntInput(&numData, &tag) != 0) {
		opserr << "WARNING invalid uniaxialMaterial Concrete0 tag" << endln;
		return 0;
	}

	//numData = 9;
	if (OPS_GetDoubleInput(&numData, &e0c) != 0) {
		opserr << "WARNING invalid e0c\n";
		opserr << "uniaxialMaterial Concrete0: " << tag << endln;
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, &fpc) != 0) {
		opserr << "WARNING invalid fpc\n";
		opserr << "uniaxialMaterial Concrete0: " << tag << endln;
		return 0;
	}

	if (OPS_GetDoubleInput(&numData, &e00) != 0) {
		opserr << "WARNING invalid fpc\n";
		opserr << "uniaxialMaterial Concrete0: " << tag << endln;
		return 0;
	}

	if (numArgs == 6) {

		// Parsing was successful, allocate the material
		theMaterial = new Concrete0(tag, e0c, fpc, e00);

	}
	else if (numArgs == 8) {

		//numData = 1;
		//int mon;
		if (OPS_GetDoubleInput(&numData, &e0t) != 0) {
			opserr << "WARNING invalid e0t\n";
			opserr << "uniaxialMaterial Concrete0: " << tag << endln;
			return 0;
		}

		if (OPS_GetDoubleInput(&numData, &ft) != 0) {
			opserr << "WARNING invalid ft\n";
			opserr << "uniaxialMaterial Concrete0: " << tag << endln;
			return 0;
		}

		// Parsing was successful, allocate the material
		theMaterial = new Concrete0(tag, e0c, fpc, e00, e0t, ft);

	}
	else {

		//int gap;
		//numData = 1;

		if (OPS_GetDoubleInput(&numData, &e0t) != 0) {
			opserr << "WARNING invalid e0t\n";
			opserr << "uniaxialMaterial Concrete0: " << tag << endln;
			return 0;
		}

		if (OPS_GetDoubleInput(&numData, &ft) != 0) {
			opserr << "WARNING invalid ft\n";
			opserr << "uniaxialMaterial Concrete0: " << tag << endln;
			return 0;
		}

		if (OPS_GetIntInput(&numData, &mon) != 0) {
			opserr << "WARNING invalid uniaxialMaterial Concrete0 tag" << endln;
			return 0;
		}

		if (mon != 0 && mon != 1) {
			opserr << "WARNING invalid $mon parameter for Concrete0" << endln;
			return 0;
		}

		theMaterial = new Concrete0(tag, e0c, fpc, e00, e0t, ft, mon);
	}

	return theMaterial;
}

// No tension capacity
Concrete0::Concrete0(int tag, double E0C, double FPC, double E00)
	:UniaxialMaterial(tag,MAT_TAG_Concrete0),
	e0c(E0C), fpc(FPC), e0t(0.0), e00(E00), ft(0.0), Et(0.0), mon(0),
	Tstrain(0.0), Tstress(0.0), Ttangent(0.0),
	Cstrain(0.0), Cstress(0.0), Ctangent(0.0),
	CompStiff(0.0),
	Cracking(0), Crushing(0),
	epsnmin(0.0), signmin(0.0),
	SMALL(1e-30) // for small stiffness
{
	if (fpc > 0.0)
		fpc = -fpc;
	if (e0c > 0.0)
		e0c = -e0c;
	if (e00 > 0.0)
		e00 = -e00;

	// Points for calculating pre-peak curve
	e0c033 = 1.0/3.0*e0c;
	fpc033 = fpc*(2.0*e0c033/e0c - pow(e0c033/e0c,2.0));

	e0c066 = 3.0/4.0*e0c;
	fpc066 = fpc*(2.0*e0c066/e0c - pow(e0c066/e0c,2.0));

	Ec0 = fpc033/e0c033;
	Ec033 = (fpc066-fpc033)/(e0c066-e0c033);
	Ec066 = (fpc-fpc066)/(e0c-e0c066);
	Ec00 = (0.0-fpc)/(e00-e0c);

	Et = SMALL;

}

// with tension capacity
Concrete0::Concrete0(int tag, double E0C, double FPC, double E00, double E0T, double FT)
	:UniaxialMaterial(tag,MAT_TAG_Concrete0),
	e0c(E0C), fpc(FPC), e0t(E0T), e00(E00), ft(FT), mon(0),
	Tstrain(0.0), Tstress(0.0), Ttangent(0.0),
	Cstrain(0.0), Cstress(0.0),
	CompStiff(0.0),
	Cracking(0), Crushing(0),
	epsnmin(0.0), signmin(0.0),
	SMALL(1e-30) // for small stiffness
{
	if (fpc > 0.0)
		fpc = -fpc;
	if (e0c > 0.0)
		e0c = -e0c;
	if (e00 > 0.0)
		e00 = -e00;
	if (ft < 0.0)
		ft = -ft;
	if (e0t < 0.0)
		e0t = -e0t;

	// Points for calculating pre-peak curve
	e0c033 = 1.0/3.0*e0c;
	fpc033 = fpc*(2.0*e0c033/e0c - pow(e0c033/e0c,2.0));

	e0c066 = 3.0/4.0*e0c;
	fpc066 = fpc*(2.0*e0c066/e0c - pow(e0c066/e0c,2.0));

	Ec0 = fpc033/e0c033;
	Ec033 = (fpc066-fpc033)/(e0c066-e0c033);
	Ec066 = (fpc-fpc066)/(e0c-e0c066);
	Ec00 = (0.0-fpc)/(e00-e0c);

	Et = ft/e0t;

}

// With tension capacity and mon parameter
Concrete0::Concrete0(int tag, double E0C, double FPC, double E00, double E0T, double FT, int MON)
	:UniaxialMaterial(tag,MAT_TAG_Concrete0),
	e0c(E0C), fpc(FPC), e0t(E0T), e00(E00), ft(FT), mon(MON),
	Tstrain(0.0), Tstress(0.0), Ttangent(0.0),
	Cstrain(0.0), Cstress(0.0), Ctangent(0.0),
	CompStiff(0.0),
	Cracking(0), Crushing(0),
	epsnmin(0.0), signmin(0.0),
	SMALL(1e-30) // for small stiffness
{
	if (fpc > 0.0)
		fpc = -fpc;
	if (e0c > 0.0)
		e0c = -e0c;
	if (e00 > 0.0)
		e00 = -e00;
	if (ft < 0.0)
		ft = -ft;
	if (e0t < 0.0)
		e0t = -e0t;

	// Points for calculating pre-peak curve
	e0c033 = 1.0/3.0*e0c;
	fpc033 = fpc*(2.0*e0c033/e0c - pow(e0c033/e0c,2.0));

	e0c066 = 3.0/4.0*e0c;
	fpc066 = fpc*(2.0*e0c066/e0c - pow(e0c066/e0c,2.0));

	Ec0 = fpc033/e0c033;
	Ec033 = (fpc066-fpc033)/(e0c066-e0c033);
	Ec066 = (fpc-fpc066)/(e0c-e0c066);
	Ec00 = (0.0-fpc)/(e00-e0c);

	Et = ft/e0t;
}

// Blank constructor
Concrete0::Concrete0()
	:UniaxialMaterial(0,MAT_TAG_Concrete0),
	e0c(0.0), fpc(0.0), e0t(0.0), ft(0.0),
	Tstrain(0.0), Tstress(0.0),
	Cracking(0), Crushing(0),
	SMALL(1e-30) // for small stiffness
{

}

Concrete0::~Concrete0()
{
  // does nothing
}


int 
Concrete0::setTrial(double strain, double &stress, double &Ttangent, double strainRate)
{

  Tstrain = strain;

  determineTrialState();
 
  return 0;
}

int 
Concrete0::setTrialStrain(double strain, double strainRate)
{

  Tstrain = strain;

  determineTrialState();
 
  return 0;
}

void Concrete0::determineTrialState()
{

	if (Tstrain >= 0.0) { // tension

		if (Cracking == 1) { // concrete has cracked
			Tstress = 0.0;
			Ttangent = SMALL;

		} else {
			if (Tstrain <= e0t) {
				Tstress = Et*Tstrain;
				Ttangent = Et;

			} else {
				Tstress = 0.0;
				Ttangent = SMALL;

			}
		}

	} else { // compression

		if (mon == 1) { // Monotonic concrete in compression for uncracked panel behavior
			compEnvelope();

		} else { // Cyclic concrete in compression

			if (Crushing == 1) { // concrete has crushed
				Tstress = 0.0;
				Ttangent = SMALL;

			} else {

				if (Tstrain < epsnmin) { // go to envelope
					compEnvelope();

				} else { // unloading/reloading region
					Tstress = CompStiff*Tstrain;
					Ttangent = CompStiff;

				}

			}

		}

	}

} 


void Concrete0::compEnvelope()
{

	if (Tstrain <= e00) { // < 1.0e00
		Tstress = 0.0;
		Ttangent = SMALL;
		
	} else if (Tstrain <= e0c) { // e00 - e0c
		Tstress = fpc + Ec00*(Tstrain-e0c);
		Ttangent = Ec00;

	} else if (Tstrain <= e0c066) { // 1.0 e0c - 0.66e0c
		Tstress = fpc066 + Ec066*(Tstrain-e0c066);
		Ttangent = Ec066;

	} else if (Tstrain <= e0c033) { // 0.33e0c - 0.66e0c
		Tstress = fpc033 + Ec033*(Tstrain-e0c033);
		Ttangent = Ec033;

	} else { // 0 - 0.33e0c
		Tstress = Ec0*Tstrain;
		Ttangent = Ec0;
	}
}

double 
Concrete0::getStress(void)
{
  return Tstress;
}

double 
Concrete0::getTangent(void)
{
  return Ttangent;
}

double 
Concrete0::getInitialTangent(void)
{
    return (Et > Ec0) ? Et : Ec0;
}

int 
Concrete0::commitState(void)
{

	if (Tstrain >= 0) { // tension

		if (Cracking == 1) { // concrete has cracked
			Cstrain = Tstrain;
			Cstress = 0.0;
			Ctangent = SMALL;

		} else {
			if (Tstrain <= e0t) {
				Cstrain = Tstrain;
				Cstress = Et*Cstrain;
				Ctangent = Et;

			} else {
				Cstrain = Tstrain;
				Cstress = 0.0;
				Ctangent = SMALL;

			}
		}

	} else { // compression

		if (Tstrain < epsnmin) { // update maximum negative strain
			epsnmin = Tstrain;
			compEnvelope();
			signmin = Tstress;
			CompStiff = signmin/epsnmin;

			Cstrain = Tstrain;
			Cstress = Tstress;
			Ctangent = Ttangent;

		} else {
			Cstrain = Tstrain;
			Cstress = CompStiff*Cstrain;
			Ctangent = CompStiff;

		}

	}

	// Concrete crushed
	if (Cstrain <= e00 && Crushing == 0 && mon == 0) {
		Crushing = 1;
	}

	// Concrete cracked
	if (Cstrain >= e0t && Cracking == 0 && mon == 0) {
		Cracking = 1;
	}

	return 0;
}


int 
Concrete0::revertToLastCommit(void) // TRY TO GET HERE AND FIGURE OUT !!!
{
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;

    return 0;
}


int 
Concrete0::revertToStart(void)
{

    Tstrain = 0.0;
	Tstress = 0.0;
	Ttangent = this->getInitialTangent();

	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = this->getInitialTangent();
	
    return 0;
}

UniaxialMaterial *
Concrete0::getCopy(void)
{
    Concrete0 *theCopy = new Concrete0(this->getTag(), e0c, fpc, e00, e0t, ft, mon);
	
	theCopy -> e0c033 = e0c033;
	theCopy -> e0c066 = e0c066;
	theCopy -> fpc033 = fpc033;
	theCopy -> fpc066 = fpc066;
    
	theCopy -> Ec0 = Ec0;
	theCopy -> Ec033 = Ec033;
	theCopy -> Ec066 = Ec066;
	theCopy -> Ec00 = Ec00;
	theCopy -> Et = Et;

	theCopy -> Cstrain = Cstrain;
	theCopy -> Cstress = Cstress;
	theCopy -> CompStiff = CompStiff;

	theCopy -> Crushing = Crushing;
	theCopy -> Cracking = Cracking;
	theCopy -> Ctangent = Ctangent;

	theCopy -> epsnmin = epsnmin;
	theCopy -> signmin = signmin;

	theCopy -> SMALL = SMALL;

    return theCopy;
}

int 
Concrete0::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(25);

  data(0) = this->getTag();
  data(1) = e0c;
  data(2) = fpc;
  data(3) = e00;
  data(4) = e0t;
  data(5) = ft;
  data(6) = mon;

  data(7) = e0c033;
  data(8) = e0c066;
  data(9) = fpc033;
  data(10) = fpc066;

  data(11) = Ec0;
  data(12) = Ec033;
  data(13) = Ec066;
  data(14) = Ec00;
  data(15) = Et;

  data(16) = Cstrain;
  data(17) = Cstress;
  data(18) = Ctangent;
  data(19) = CompStiff;

  data(20) = Crushing;
  data(21) = Cracking;

  data(22) = epsnmin;
  data(23) = signmin;

  data(24) = SMALL;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "Concrete0::sendSelf() - failed to send data\n";

  return res;
}

int 
Concrete0::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(25);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "Concrete0::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag(data(0));
    e0c = data(1);
    fpc = data(2);
	e00 = data(3);
    e0t = data(4);
	ft  = data(5);
	mon = data(6);

	e0c033 = data(7);
	e0c066 = data(8);
	fpc033 = data(9);
	fpc066 = data(10);

	Ec0 = data(11);
	Ec033 = data(12);
	Ec066 = data(13);
	Ec00 = data(14);
	Et = data(15);

	Cstrain = data(16);
	Cstress = data(17);
	Ctangent = data(18);
	CompStiff = data(19);

	Crushing = data(20);
	Cracking = data(21);

	epsnmin = data(22);
	signmin = data(23);

	SMALL = data(24);

  }

  // Set trial state variables
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;

  return res;
}

void 
Concrete0::Print(OPS_Stream &s, int flag)
{
    s << "Concrete0 tag: " << this->getTag() << endln;
    s << "  e0c: " << e0c << " fpc: " << fpc << " e00: " << e00 << " e0t: " << e0t << " ft: " << ft << endln;
}

int
Concrete0::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"e0c") == 0)
    return param.addObject(1, this);

  else if (strcmp(argv[0],"fpc") == 0)
    return param.addObject(2, this);
  
  else if (strcmp(argv[0],"e00") == 0)
    return param.addObject(3, this);

  else if (strcmp(argv[0],"e0t") == 0)
    return param.addObject(4, this);
  
  else if (strcmp(argv[0],"ft") == 0)
    return param.addObject(5, this);
  
  return -1;
}
  
int 
Concrete0::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case -1:
    return -1;
  case 1:
    e0c = info.theDouble;
    return 0;
  case 2:
    fpc = info.theDouble;
    return 0;
  case 3:
    e00 = info.theDouble;
    return 0;
  case 4:
    e0t = info.theDouble;
    return 0;
  case 5:
    ft = info.theDouble;
    return 0;
  default:
    return -1;
  }

  return 0;
}

// KK
Response* Concrete0::setResponse(const char **argv, int argc,
	OPS_Stream &theOutput)
{
	Response *theResponse = 0;

	if (strcmp(argv[0],"getInputParameters") == 0) {
		Vector data1(6);
		data1.Zero();
		theResponse = new MaterialResponse(this, 100, data1);

	} else

		return this->UniaxialMaterial::setResponse(argv, argc, theOutput);

	return theResponse;
}

// KK
int Concrete0::getResponse(int responseID, Information &matInfo)
{
	if (responseID == 100) {
		matInfo.setVector(this->getInputParameters()); 

	} else

		return this->UniaxialMaterial::getResponse(responseID, matInfo);

	return 0;
}

// KK
Vector Concrete0::getInputParameters(void)
{
	Vector input_par(6); // size = max number of parameters (assigned + default)

	input_par.Zero();

	input_par(0) = this->getTag(); 
	input_par(1) = e0c; 
	input_par(2) = fpc;
	input_par(3) = e00;
	input_par(4) = e0t; 
	input_par(5) = ft; 

	return input_par;
}