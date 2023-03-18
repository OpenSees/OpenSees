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
  SectionForceDeformation(0, SEC_TAG_CreepSection), theSection(0), creepFactor(0.0)
{

}

CreepSection::CreepSection(int tag, SectionForceDeformation &section):
  SectionForceDeformation(tag, SEC_TAG_CreepSection), theSection(0), creepFactor(0.0)
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
}

SectionForceDeformation *
CreepSection::getCopy(void)
{
  if (theSection == 0)
    return 0;
  
  CreepSection *theCopy = new CreepSection(this->getTag(), *theSection);

  theCopy->creepFactor = creepFactor;
  
  return theCopy;
}

void
CreepSection::Print(OPS_Stream &s, int flag)
{
  s << "CreepSection: " << this->getTag() << endln;
  s << "  creep factor: " << creepFactor << endln;
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
  
  if (strcmp(argv[0],"creepFactor") == 0) {
    param.setValue(creepFactor);
    return param.addObject(1, this);
  }

  return -1;
}

int
CreepSection::updateParameter(int paramID, Information &info)
{
  if (paramID != 1)
    return 0;
  
  creepFactor = info.theDouble;

  const char *argv[3];

  argv[0] = "fiberData";
  DummyStream stream;
  Response *theResponse = theSection->setResponse(argv, 1, stream);
  theResponse->getResponse();
  Information &secinfo = theResponse->getInformation();
  const Vector &theVector = *(secinfo.theVector);

  // Need to make sure this is fiber section
  int numFibers = theVector.Size() / 5;
  
  argv[0] = "fiberIndex";
  argv[2] = "epsInit";
  char buffer[80];
  for (int i = 0; i < numFibers; i++) {
    sprintf(buffer,"%d",i);
    argv[1] = buffer;
    
    Parameter param;
    int ok = theSection->setParameter(argv,3,param);
    if (ok < 0)
      continue;

    double eps0 = theVector(4 + i*5);
    param.update(creepFactor*eps0);
  }
  
  
  //  param.update(creepFactor);
  
  return 0;
}
