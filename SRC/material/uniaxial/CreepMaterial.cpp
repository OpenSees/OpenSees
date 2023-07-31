#include <CreepMaterial.h>

#include <elementAPI.h>

void *OPS_CreepMaterial()
{
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 2) {
    opserr << "CreepMaterial - insufficient input args" << endln;
    return 0;
  }

  int tag, matTag;
  int numData = 1;
  
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "CreepMaterial - unable to read tag" << endln;
    return 0;
  }

  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "CreepMaterial - unable to read material tag" << endln;
    return 0;
  }  

  UniaxialMaterial *wrappedMaterial = OPS_getUniaxialMaterial(matTag);
  if (wrappedMaterial == 0) {
    opserr << "CreepMaterial - material with tag " << matTag << " does not exist" << endln;
    return 0;
  }

  theMaterial = new CreepMaterial(tag, *wrappedMaterial);
  
  return theMaterial;
}

CreepMaterial::CreepMaterial(int tag, UniaxialMaterial &wrapped):
  UniaxialMaterial(tag, MAT_TAG_CreepMaterial), wrappedMaterial(0)
{
  wrappedMaterial = wrapped.getCopy();
  if (wrappedMaterial == 0) {
    opserr << "CreepMaterial::CreepMaterial - unable to get copy of material" << endln;
    exit(-1);
  }

  for (int i = 0; i < maxNumSteps; i++) {
    PHI_i[i] = 0.0;
    E_i[i] = 0.0;
    DSIG_i[i] = 0.0;
    dsig_i[i] = 0.0;
    TIME_i[i] = 0.0;
    DTIME_i[i] = 0.0;
  }
}

CreepMaterial::CreepMaterial():
  UniaxialMaterial(0, MAT_TAG_CreepMaterial), wrappedMaterial(0)
{
  for (int i = 0; i < maxNumSteps; i++) {
    PHI_i[i] = 0.0;
    E_i[i] = 0.0;
    DSIG_i[i] = 0.0;
    dsig_i[i] = 0.0;
    TIME_i[i] = 0.0;
    DTIME_i[i] = 0.0;
  }
}

CreepMaterial::~CreepMaterial()
{
  if (wrappedMaterial != 0)
    delete wrappedMaterial;
}

UniaxialMaterial *
CreepMaterial::getCopy(void)
{
  if (wrappedMaterial == 0)
    return 0;
  
  CreepMaterial *theCopy = new CreepMaterial(this->getTag(), *wrappedMaterial);

  for (int i = 0; i < maxNumSteps; i++) {
    theCopy->PHI_i[i] = PHI_i[i];
    theCopy->E_i[i] = E_i[i];
    theCopy->DSIG_i[i] = DSIG_i[i];
    theCopy->dsig_i[i] = dsig_i[i];
    theCopy->TIME_i[i] = TIME_i[i];
    theCopy->DTIME_i[i] = DTIME_i[i];
  }

  return theCopy;
}

int
CreepMaterial::setTrialStrain(double strain, double strainRate)
{
  if (wrappedMaterial == 0)
    return 0;
  
  int ret = wrappedMaterial->setTrialStrain(strain, strainRate);

  return ret;
}

double
CreepMaterial::getStrain(void)
{
  if (wrappedMaterial == 0)
    return 0.0;
  
  return wrappedMaterial->getStrain();
}

double
CreepMaterial::getStress(void)
{
  if (wrappedMaterial == 0)
    return 0.0;
  
  return wrappedMaterial->getStress();
}

double
CreepMaterial::getTangent(void)
{
  if (wrappedMaterial == 0)
    return 0.0;
  
  return wrappedMaterial->getTangent();
}

double
CreepMaterial::getInitialTangent(void)
{
  if (wrappedMaterial == 0)
    return 0.0;
  
  return wrappedMaterial->getInitialTangent();
}

int
CreepMaterial::commitState(void)
{
  if (wrappedMaterial == 0)
    return 0;

  return wrappedMaterial->commitState();
}

int
CreepMaterial::revertToLastCommit(void)
{
  if (wrappedMaterial == 0)
    return 0;

  return wrappedMaterial->revertToLastCommit();
}

int
CreepMaterial::revertToStart(void)
{
  if (wrappedMaterial == 0)
    return 0;

  return wrappedMaterial->revertToStart();
}

int
CreepMaterial::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
CreepMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
CreepMaterial::Print(OPS_Stream &s, int flag)
{
  s << "CreepMaterial, tag = " << this->getTag() << endln;
  if (wrappedMaterial != 0)
    s << "\twrapped material tag = " << wrappedMaterial->getTag() << endln;
  else
    s << "\tNo wrapped material" << endln;
}

Response *
CreepMaterial::setResponse(const char **argv, int argc, OPS_Stream &theOutput)
{
  if (wrappedMaterial == 0)
    return UniaxialMaterial::setResponse(argv, argc, theOutput);
  else
    return wrappedMaterial->setResponse(argv, argc, theOutput);
}

int
CreepMaterial::getResponse(int responseID, Information &matInfo)
{
  return 0;
}
