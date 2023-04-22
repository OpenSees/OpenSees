#ifndef CreepSection_h
#define CreepSection_h

#include <SectionForceDeformation.h>
#include <Vector.h>
#include <Matrix.h>

class CreepSection : public SectionForceDeformation
{
 public:
  CreepSection();
  CreepSection(int tag, SectionForceDeformation &section);
  ~CreepSection();

  const char *getClassType(void) const {return "CreepSection";}

  int setTrialSectionDeformation(const Vector &e) {return theSection->setTrialSectionDeformation(e);}
  const Vector &getSectionDeformation(void) {return theSection->getSectionDeformation();}

  const Vector &getStressResultant(void) {return theSection->getStressResultant();}
  const Matrix &getSectionTangent(void) {return theSection->getSectionTangent();}
  const Matrix &getInitialTangent(void) {return theSection->getInitialTangent();}

  int commitState(void) {return theSection->commitState();}
  int revertToLastCommit(void) {return theSection->revertToLastCommit();}
  int revertToStart(void) {return theSection->revertToStart();}

  SectionForceDeformation *getCopy(void);
  const ID &getType(void) {return theSection->getType();}
  int getOrder(void) const {return theSection->getOrder();}

  int sendSelf(int cTag, Channel &theChannel) {return -1;}
  int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {return -1;}

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s) {return theSection->setResponse(argv, argc, s);}
  int getResponse(int responseID, Information &info) {return theSection->getResponse(responseID, info);}
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  const Vector & getStressResultantSensitivity(int gradIndex, bool conditional) {return theSection->getStressResultantSensitivity(gradIndex, conditional);}
  const Matrix & getSectionTangentSensitivity(int gradIndex) {return theSection->getSectionTangentSensitivity(gradIndex);}
  int   commitSensitivity(const Vector& dedh, int gradIndex, int numGrads) {return theSection->commitSensitivity(dedh, gradIndex, numGrads);}
  const Vector & getSectionDeformationSensitivity(int gradIndex) {return theSection->getSectionDeformationSensitivity(gradIndex);}

  double getEnergy(void) const {return theSection->getEnergy();}
  
 private:
  SectionForceDeformation *theSection;
  double creepFactor;
  double shrinkage;

  int numFibers;
  double *initialStrain;
};

#endif
