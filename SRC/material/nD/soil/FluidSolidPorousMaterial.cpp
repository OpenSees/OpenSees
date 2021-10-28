// $Revision: 1.26 $
// $Date: 2010-06-12 18:07:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/FluidSolidPorousMaterial.cpp,v $
                                                                        
// Written: ZHY

//
// FluidSolidPorousMaterial.cpp
// -------------------
//
#include <math.h>
#include <stdlib.h>
#include <FluidSolidPorousMaterial.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <MultiYieldSurface.h>

#include <string.h>
#include <elementAPI.h>

int* FluidSolidPorousMaterial::loadStagex = 0;
int* FluidSolidPorousMaterial::ndmx = 0;
double* FluidSolidPorousMaterial::combinedBulkModulusx = 0;
int FluidSolidPorousMaterial::matCount = 0;
double FluidSolidPorousMaterial::pAtm = 101;

Vector FluidSolidPorousMaterial::workV3(3);
Vector FluidSolidPorousMaterial::workV6(6);
Matrix FluidSolidPorousMaterial::workM3(3,3);
Matrix FluidSolidPorousMaterial::workM6(6,6);

void* OPS_FluidSolidPorousMaterial()
{
    int tag;  double param[4];
    char * arg[] = {"nd", "soilMatTag", "combinedBulkModul", "Atmospheric pressure"};

    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < 6) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial FluidSolidPorous tag? "<< arg[0];
	opserr << "? "<< "\n";
	opserr << arg[1] << "? "<< arg[2] << "? "<< endln;
	return 0;
    }

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid FluidSolidPorous tag" << endln;
	return 0;
    }

    for (int i=3; i<6; i++)
	if (OPS_GetDoubleInput(&numdata, &param[i-3]) < 0) {
	    opserr << "WARNING invalid " << " double" << "\n";
	    opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	    return 0;
	}

    NDMaterial *soil = OPS_getNDMaterial(param[1]);
    if (soil == 0) {
	opserr << "WARNING FluidSolidPorous: couldn't get soil material ";
	opserr << "tagged: " << param[1] << "\n";
	return 0;
    }

    param[3] = 101.;
    if (argc == 7) {
	if (OPS_GetDoubleInput(&numdata, &param[3]) < 0) {
	    opserr << "WARNING invalid " << " double" << "\n";
	    opserr << "nDMaterial FluidSolidPorous: " << tag << endln;
	    return 0;
	}
    }

    return new FluidSolidPorousMaterial (tag, param[0], *soil,
					 param[2],param[3]);
}

FluidSolidPorousMaterial::FluidSolidPorousMaterial (int tag, int nd, NDMaterial &soilMat,
						    double combinedBulkModul, double atm)
 : NDMaterial(tag, ND_TAG_FluidSolidPorousMaterial)
{
  if (combinedBulkModul < 0) {
    opserr << "WARNING:FluidSolidPorousMaterial::FluidSolidPorousMaterial: combinedBulkModulus < 0" << endln;
    opserr << "Will reset to 0." <<endln;
    combinedBulkModul = 0.;
  }

  if (matCount%20 == 0) {
     int * temp1 = loadStagex;
     int * temp2 = ndmx;
     double * temp3 = combinedBulkModulusx;
     loadStagex = new int[matCount+20];
     ndmx = new int[matCount+20];
     combinedBulkModulusx = new double[matCount+20];
     for (int i=0; i<matCount; i++) {
       loadStagex[i] = temp1[i];
       ndmx[i] = temp2[i];
       combinedBulkModulusx[i] = temp3[i];
     }
     if (matCount > 0) {
       delete [] temp1; delete [] temp2; delete [] temp3;
     }
  }
  
  ndmx[matCount] = nd;
  loadStagex[matCount] = 0;  //default
  combinedBulkModulusx[matCount] = combinedBulkModul;
  matN = matCount;
  matCount ++;
  pAtm = atm;

  theSoilMaterial = soilMat.getCopy();
  theSoilCommittedStress = theSoilMaterial->getStress();
  theSoilCommittedStrain = theSoilMaterial->getStrain();

  trialExcessPressure = currentExcessPressure = 0.;
  trialVolumeStrain = currentVolumeStrain = 0.;
  initMaxPress = 0.;
  e2p = 0;
}
   

FluidSolidPorousMaterial::FluidSolidPorousMaterial () 
 : NDMaterial(0,ND_TAG_FluidSolidPorousMaterial), theSoilMaterial(0)
{
  trialExcessPressure = currentExcessPressure = 0.;
  trialVolumeStrain = currentVolumeStrain = 0.;
  initMaxPress = 0.;
  e2p = 0;
}


FluidSolidPorousMaterial::FluidSolidPorousMaterial (const FluidSolidPorousMaterial & a)
 : NDMaterial(a.getTag(),ND_TAG_FluidSolidPorousMaterial)
{
  matN = a.matN;
  theSoilMaterial = a.theSoilMaterial->getCopy();
  trialExcessPressure = a.trialExcessPressure;
  currentExcessPressure = a.currentExcessPressure;
  trialVolumeStrain = a.trialVolumeStrain;
  currentVolumeStrain = a.currentVolumeStrain;
  initMaxPress = a.initMaxPress;
  e2p = a.e2p;
}


FluidSolidPorousMaterial::~FluidSolidPorousMaterial ()
{
  if (theSoilMaterial != 0)
    delete theSoilMaterial;
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain)
{
  int ndm = ndmx[matN];
  
  if (ndm==2 && strain.Size()==3)
    trialVolumeStrain = strain[0]+strain[1];
  else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
  else {
    opserr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
    exit(-1);;
  }

  return theSoilMaterial->setTrialStrain(strain);
}


int FluidSolidPorousMaterial::setTrialStrain (const Vector &strain, const Vector &rate)
{
	int ndm = ndmx[matN];

	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = strain[0]+strain[1]+strain[2];
	else {
		opserr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);;
	}

  return theSoilMaterial->setTrialStrain(strain, rate);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain)
{
	int ndm = ndmx[matN];

	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		opserr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);;
	}

  return theSoilMaterial->setTrialStrainIncr(strain);
}


int FluidSolidPorousMaterial::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
	int ndm = ndmx[matN];

	if (ndm==2 && strain.Size()==3)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1];
	else if (ndm==3 && strain.Size()==6)
		trialVolumeStrain = currentVolumeStrain + strain[0]+strain[1]+strain[2];
	else {
		opserr << "Fatal:FluidSolidPorousMaterial:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);;
	}

  return theSoilMaterial->setTrialStrainIncr(strain, rate);
}


const Matrix & FluidSolidPorousMaterial::getTangent (void)
{
	int ndm = ndmx[matN];
	int loadStage = loadStagex[matN];
	double combinedBulkModulus = combinedBulkModulusx[matN];

	Matrix *workM = (ndm == 2) ? &workM3 : &workM6;
  
	*workM = theSoilMaterial->getTangent();

	if (loadStage != 0) { 
	  for (int i=0; i<ndm; i++) 
		  for (int j=0; j<ndm; j++) 
			  (*workM)(i,j) = (*workM)(i,j) + combinedBulkModulus;
  }

	return *workM;
}

const Matrix & FluidSolidPorousMaterial::getInitialTangent (void)
{
	int ndm = ndmx[matN];

	Matrix *workM = (ndm == 2) ? &workM3 : &workM6;
  
	*workM = theSoilMaterial->getInitialTangent();

	return *workM;
}

double FluidSolidPorousMaterial::getRho(void)
{
  return theSoilMaterial->getRho();
}

const Vector & FluidSolidPorousMaterial::getStress (void)
{
  int ndm = ndmx[matN];
  int loadStage = loadStagex[matN];
  double combinedBulkModulus = combinedBulkModulusx[matN];
  
  Vector *workV = (ndm == 2) ? &workV3 : &workV6;
  
  *workV = theSoilMaterial->getStress();
  
  if (loadStage != 0) { 
    if (e2p==0) {
      e2p = 1;
      initMaxPress = ((*workV)[0] < (*workV)[1]) ? (*workV)[0] : (*workV)[1]; 
      if (ndm == 3)
	initMaxPress = (initMaxPress < (*workV)[2]) ? initMaxPress : (*workV)[2];
    }
    trialExcessPressure = currentExcessPressure;
    trialExcessPressure += 
      (trialVolumeStrain - currentVolumeStrain) * combinedBulkModulus;
    if (trialExcessPressure > pAtm-initMaxPress) 
      trialExcessPressure = pAtm-initMaxPress;
    //if (trialExcessPressure < initMaxPress)
    //	trialExcessPressure = initMaxPress;
    for (int i=0; i<ndm; i++) 
      (*workV)[i] += trialExcessPressure;
  }
  
  return *workV;
}


int FluidSolidPorousMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 2)
    return theSoilMaterial->setParameter(argv, argc, param);
  
  int matTag = atoi(argv[1]);
  if (this->getTag() == matTag) {

    if (strcmp(argv[0],"updateMaterialStage") == 0) 
    	return param.addObject(1, this);  
    else if (strcmp(argv[0],"combinedBulkModulus") == 0) 
      return param.addObject(2, this);  
  }

  return theSoilMaterial->setParameter(argv, argc, param);
}

int FluidSolidPorousMaterial::updateParameter(int responseID, Information &info)
{
  if (responseID == 1) {
    //    opserr << "FluidSolidPorousMaterial updateMaterialStage" << info.theInt << endln;
    loadStagex[matN] = info.theInt;
  }
  else if (responseID == 2) {
    //    opserr << "FluidSolidPorousMaterial combinedBulk" << info.theDouble << endln;
    combinedBulkModulusx[matN]=info.theDouble;
  }

  return 0;
}


const Vector & FluidSolidPorousMaterial::getCommittedStress (void)
{
  return theSoilCommittedStress;
}


const Vector & FluidSolidPorousMaterial::getCommittedStrain (void)
{
  return theSoilCommittedStrain;
}


const Vector & FluidSolidPorousMaterial::getCommittedPressure (void)
{
	int ndm = ndmx[matN];

	static Vector temp(2);
	
	temp[0] = currentExcessPressure;
	temp[1] = temp[0]/initMaxPress;
    
	return temp;
}


const Vector & FluidSolidPorousMaterial::getStrain (void)
{
  return theSoilMaterial->getStrain();
}


int FluidSolidPorousMaterial::commitState (void)
{
	int loadStage = loadStagex[matN];

	currentVolumeStrain = trialVolumeStrain;
	if (loadStage != 0) 
		currentExcessPressure = trialExcessPressure;
	else
        currentExcessPressure = 0.;

	int res = theSoilMaterial->commitState();
        theSoilCommittedStress = theSoilMaterial->getStress();
	theSoilCommittedStrain = theSoilMaterial->getStrain();
	return res;
}


int FluidSolidPorousMaterial::revertToLastCommit (void)
{
	return theSoilMaterial->revertToLastCommit();
}

int FluidSolidPorousMaterial::revertToStart (void)
{
	return theSoilMaterial->revertToStart();
}

NDMaterial * FluidSolidPorousMaterial::getCopy (void)
{
  FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
  return copy;
}


NDMaterial * FluidSolidPorousMaterial::getCopy (const char *code)
{
	if (strcmp(code,"FluidSolidPorous") == 0 || strcmp(code,"PlaneStrain") == 0 ||
		strcmp(code,"ThreeDimensional") == 0) {
     FluidSolidPorousMaterial * copy = new FluidSolidPorousMaterial(*this);
	   return copy;
	}

	return 0;
}


const char * FluidSolidPorousMaterial::getType (void) const
{
	int ndm = ndmx[matN];

	return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int FluidSolidPorousMaterial::getOrder (void) const
{
	int ndm = ndmx[matN];

	return (ndm == 2) ? 3 : 6;
}


int FluidSolidPorousMaterial::sendSelf(int commitTag, Channel &theChannel)
{
	int ndm = ndmx[matN];
	int loadStage = loadStagex[matN];
	double combinedBulkModulus = combinedBulkModulusx[matN];

	int res = 0;

	static Vector data(7);
	data(0) = this->getTag();
	data(1) = ndm;
	data(2) = loadStage;
    data(3) = combinedBulkModulus;
	data(4) = currentExcessPressure;
    data(5) = currentVolumeStrain;
	data(6) = matN;

    res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
	  opserr << "FluidSolidPorousMaterial::sendSelf -- could not send Vector\n";
	     
		return res;
	}

	ID classTags(2);

	classTags(0) = theSoilMaterial->getClassTag();
	int matDbTag = theSoilMaterial->getDbTag();
	// NOTE: we do have to ensure that the material has a database
	// tag if we are sending to a database channel.
	if (matDbTag == 0) {
		matDbTag = theChannel.getDbTag();
		if (matDbTag != 0)
			theSoilMaterial->setDbTag(matDbTag);
	}
	classTags(1) = matDbTag;

	res += theChannel.sendID(this->getDbTag(), commitTag, classTags);
	if (res < 0) {
	  opserr << "WARNING FluidSolidPorousMaterial::sendSelf() - " << this->getTag() << " failed to send ID\n";
	  
	  return res;
	}

	// Finally, asks the material object to send itself
	res += theSoilMaterial->sendSelf(commitTag, theChannel);
	if (res < 0) {
	  opserr << "WARNING FluidSolidPorousMaterial::sendSelf() - " << this->getTag() << " failed to send its Material\n";
	  return res;
	}

	return res;
}




int FluidSolidPorousMaterial::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	int res = 0;

	static Vector data(7);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
	  opserr << "FluidSolidPorousMaterial::recvSelf -- could not receive Vector\n";
	  return res;
	}
    
	this->setTag((int)data(0));
	int ndm = data(1);
	int loadStage = data(2);
	double combinedBulkModulus = data(3);
	currentExcessPressure = data(4);
	currentVolumeStrain = data(5);
	matN = data(6);

	ndmx[matN] = ndm;
	loadStagex[matN] = loadStage;
	combinedBulkModulusx[matN] = combinedBulkModulus;

	// now receives the ids of its material
	ID classTags(2);

	res += theChannel.recvID(this->getDbTag(), commitTag, classTags);
	if (res < 0)  {
	  opserr << "FluidSolidPorousMaterial::recvSelf() - failed to recv ID data\n";
	  return res;
	}    
	
	int matClassTag = classTags(0);
	int matDbTag = classTags(1);
	// Check that material is of the right type; if not,
	// delete it and create a new one of the right type
	if (theSoilMaterial == 0 || theSoilMaterial->getClassTag() != matClassTag) {
	  if (theSoilMaterial != 0)
	    delete theSoilMaterial;
	  theSoilMaterial = theBroker.getNewNDMaterial(matClassTag);
	  if (theSoilMaterial == 0) {
	    opserr << "FluidSolidPorousMaterial::recvSelf() - " <<
	      "Broker could not create NDMaterial of class type" << matClassTag << endln;
	    exit(-1);
	  }
	}

	// Receive the material
	theSoilMaterial->setDbTag(matDbTag);
	res += theSoilMaterial->recvSelf(commitTag, theChannel, theBroker);
	if (res < 0) {
	  opserr << "FluidSolidPorousMaterial::recvSelf() - material failed to recv itself\n";
	  return res;
	}
	theSoilCommittedStress = theSoilMaterial->getStress();
	theSoilCommittedStress = theSoilMaterial->getStrain();

	return res;
}


Response*
FluidSolidPorousMaterial::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
    return new MaterialResponse(this, 1, this->getCommittedStress());
  
  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
    return new MaterialResponse(this, 2, this->getCommittedStrain());
  
  else if (strcmp(argv[0],"tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());
  
  else if (strcmp(argv[0],"backbone") == 0)
    return theSoilMaterial->setResponse(argv, argc, output);
  
  else if (strcmp(argv[0],"pressure") == 0)
    return new MaterialResponse(this, 5, this->getCommittedPressure());
  
  else
    return 0;
}


int FluidSolidPorousMaterial::getResponse (int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getCommittedStress());
    
  case 2:
    return matInfo.setVector(this->getCommittedStrain());
    
  case 3:
    return matInfo.setMatrix(this->getTangent());
    
  case 4:
    return theSoilMaterial->getResponse(responseID, matInfo);
    
  case 5:
    return matInfo.setVector(this->getCommittedPressure());
    
  default:
    return -1;
  }
}


void FluidSolidPorousMaterial::Print(OPS_Stream &s, int flag )
{
	s << "FluidSolidPorousMaterial" << endln;
}















