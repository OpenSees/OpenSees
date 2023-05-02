#include <CreepSection.h>

#include <Parameter.h>
#include <Response.h>
#include <DummyStream.h>

#include <classTags.h>
#include <elementAPI.h>

void *OPS_CreepSection()
{
  int numData = OPS_GetNumRemainingInputArgs();
  if (numData < 2) {
    opserr << "Insufficient arguments to CreepSection" << endln;
    return 0;
  }

  numData = 1;
  int tag;
  if (OPS_GetIntInput(&numData, &tag) < 0)
    return 0;

  int secTag;
  if (OPS_GetIntInput(&numData, &secTag) < 0)
    return 0;  

  SectionForceDeformation *section = OPS_getSectionForceDeformation(secTag);
  if (section == 0) {
    opserr << "CreepSection - section with tag " << secTag << " not found" << endln;
    return 0;
  }

  return new CreepSection(tag, *section);
}

CreepSection::CreepSection():
  SectionForceDeformation(0, SEC_TAG_CreepSection), theSection(0),
  creepFactor(0.0), shrinkage(0.0), numFibers(0), initialStrain(0)
{

}

CreepSection::CreepSection(int tag, SectionForceDeformation &section):
  SectionForceDeformation(tag, SEC_TAG_CreepSection), theSection(0),
  creepFactor(0.0), shrinkage(0.0), numFibers(0), initialStrain(0)
{
  theSection = section.getCopy();
  if (theSection == 0) {
    opserr << "CreepSection::CreepSection - failed to get copy of section" << endln;
    exit(-1);
  }
}

CreepSection::~CreepSection()
{
  if (theSection != 0)
    delete theSection;

  if (initialStrain != 0)
    delete [] initialStrain;
}

SectionForceDeformation *
CreepSection::getCopy(void)
{
  if (theSection == 0)
    return 0;
  
  CreepSection *theCopy = new CreepSection(this->getTag(), *theSection);

  theCopy->creepFactor = creepFactor;
  theCopy->shrinkage = shrinkage;
  theCopy->numFibers = numFibers;
  
  if (initialStrain != 0) {
    theCopy->initialStrain = new double[numFibers];
    for (int i = 0; i < numFibers; i++)
      theCopy->initialStrain[i] = initialStrain[i];
  }
  
  return theCopy;
}

void
CreepSection::Print(OPS_Stream &s, int flag)
{
  s << "CreepSection: " << this->getTag() << endln;
  s << "  creep factor: " << creepFactor << endln;
  s << "  shrinkage: " << shrinkage << endln;  
  if (theSection != 0) {
    s << "  wrapped section: " << theSection->getTag() << endln;
    theSection->Print(s, flag);
  }
}

int
CreepSection::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  if (strcmp(argv[0],"creepFactor") == 0 || strcmp(argv[0],"creepFactorInitial") == 0) {
    param.setValue(creepFactor);
    return param.addObject(1, this);
  }

  if (strcmp(argv[0],"creepFactorCurrent") == 0) {
    param.setValue(creepFactor);
    return param.addObject(3, this);
  }  

  if (strcmp(argv[0],"shrinkageStrain") == 0) {
    param.setValue(shrinkage);
    return param.addObject(2, this);
  }

  return -1;
}

int
CreepSection::updateParameter(int paramID, Information &info)
{
  if (paramID == 1) {
    creepFactor = info.theDouble;
    
    const char *argv[3];
    
    argv[0] = "fiberData";
    DummyStream stream;
    Response *theResponse = theSection->setResponse(argv, 1, stream);
    theResponse->getResponse();
    Information &secinfo = theResponse->getInformation();
    const Vector &theVector = *(secinfo.theVector);
    
    // Need to make sure this is fiber section
    int nFibers = theVector.Size() / 5;
    
    const char *argvParam[3];
    argvParam[0] = "fiberIndex";
    argvParam[2] = "epsInit";
    
    const char *argvResp[4];
    argvResp[0] = "fiberIndex";
    argvResp[2] = "material";
    argvResp[3] = "strain";
    char buffer[80];
    
    if (initialStrain == 0) {
      numFibers = nFibers;
      initialStrain = new double[numFibers];

      for (int i = 0; i < numFibers; i++) {
	initialStrain[i] = 0.0;
	
	sprintf(buffer,"%d",i);
	argvResp[1] = buffer;
	
	// Get the mechanical strain
	//double eps0 = theVector(4 + i*5);
	Response *theResponse = theSection->setResponse(argvResp, 4, stream);
	if (theResponse == 0)
	  continue;
	theResponse->getResponse();
	Information &secinfo = theResponse->getInformation();
	double eps0 = secinfo.theDouble;
	initialStrain[i] = eps0;
	//opserr << ' ' << i << ' ' << eps0 << endln;

	delete theResponse;
      }
    }
    
    for (int i = 0; i < numFibers; i++) {
      sprintf(buffer,"%d",i);
      argvParam[1] = buffer;

      // Get the initial strain parameter
      Parameter param;
      int ok = theSection->setParameter(argvParam,3,param);
      if (ok < 0)
	continue;
	
      // Update the initial strain parameter
      //if (eps0 < 0.0)
      param.update(-creepFactor*initialStrain[i] + shrinkage);
    }
    
    delete theResponse;
    //  param.update(creepFactor);
  }

  if (paramID == 3) {
    creepFactor = info.theDouble;
    
    const char *argv[3];
    
    argv[0] = "fiberData";
    DummyStream stream;
    Response *theResponse = theSection->setResponse(argv, 1, stream);
    theResponse->getResponse();
    Information &secinfo = theResponse->getInformation();
    const Vector &theVector = *(secinfo.theVector);
    
    // Need to make sure this is fiber section
    int nFibers = theVector.Size() / 5;
    
    const char *argvParam[3];
    argvParam[0] = "fiberIndex";
    argvParam[2] = "epsInit";
    
    const char *argvResp[4];
    argvResp[0] = "fiberIndex";
    argvResp[2] = "material";
    argvResp[3] = "strain";
    char buffer[80];
    
    numFibers = nFibers;
    double *currentStrain = new double[numFibers];

    for (int i = 0; i < numFibers; i++) {
      currentStrain[i] = 0.0;
	
      sprintf(buffer,"%d",i);
      argvResp[1] = buffer;
      
      // Get the mechanical strain
      //double eps0 = theVector(4 + i*5);
      Response *theResponse = theSection->setResponse(argvResp, 4, stream);
      if (theResponse == 0)
	continue;
      theResponse->getResponse();
      Information &secinfo = theResponse->getInformation();
      double eps0 = secinfo.theDouble;
      currentStrain[i] = eps0;
      //opserr << ' ' << i << ' ' << eps0 << endln;
      
      delete theResponse;

      sprintf(buffer,"%d",i);
      argvParam[1] = buffer;

      // Get the initial strain parameter
      Parameter param;
      int ok = theSection->setParameter(argvParam,3,param);
      if (ok < 0)
	continue;
	
      // Update the initial strain parameter
      //if (eps0 < 0.0)
      param.update(-creepFactor*currentStrain[i] + shrinkage);
    }

    delete [] currentStrain;
    delete theResponse;
    //  param.update(creepFactor);
  }  

  if (paramID == 2) {
    shrinkage = info.theDouble;
    
    const char *argv[3];
    
    argv[0] = "fiberData";
    DummyStream stream;
    Response *theResponse = theSection->setResponse(argv, 1, stream);
    theResponse->getResponse();
    Information &secinfo = theResponse->getInformation();
    const Vector &theVector = *(secinfo.theVector);
    
    // Need to make sure this is fiber section
    int numFibers = theVector.Size() / 5;
    
    const char *argvParam[3];
    argvParam[0] = "fiberIndex";
    argvParam[2] = "epsInit";
    
    const char *argvResp[4];
    argvResp[0] = "fiberIndex";
    argvResp[2] = "material";
    argvResp[3] = "strain";
    char buffer[80];
    for (int i = 0; i < numFibers; i++) {
      sprintf(buffer,"%d",i);
      argvParam[1] = buffer;
      argvResp[1] = buffer;
      
      // Get the initial strain parameter
      Parameter param;
      int ok = theSection->setParameter(argvParam,3,param);
      if (ok < 0)
	continue;
      
      // Update the initial strain parameter
      param.update(shrinkage);
    }
    
    delete theResponse;
    //  param.update(creepFactor);
  }
  
  return 0;
}
