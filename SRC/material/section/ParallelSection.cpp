/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.3 $
// $Date: 2008-08-26 16:48:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ParallelSection.cpp,v $
                                                                        
                                                                        
// File: ~/section/ParallelSection.C
//
// Written: MHS
// Created: June 2000
// Revision: A
//
// Purpose: This file contains the implementation for the ParallelSection class. 

#include <stdlib.h>
#include <string.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <ParallelSection.h>
#include <MaterialResponse.h>
#include <ID.h>

#include <classTags.h>
#include <elementAPI.h>

#define maxOrder 10

// Assumes section order is less than or equal to maxOrder.
// Can increase if needed!!!
double ParallelSection::workArea[2*maxOrder*(maxOrder+1)];
int    ParallelSection::codeArea[maxOrder];

void* OPS_ParallelSection()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section Parallel tag? tag1? tag2? ..." << endln;
	return 0;
    }
 
    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid section Parallel tag" << endln;
	return 0;
    }

    int numMaterials = OPS_GetNumRemainingInputArgs();
	
    if (numMaterials == 0) {
	opserr << "WARNING no component section(s) provided\n";
	opserr << "section Parallel: " << tag << endln;
	return 0;
    }
    
    // Create an array to hold pointers to component materials
    SectionForceDeformation **theMats = new SectionForceDeformation *[numMaterials];
	
    // For each material get the tag and ensure it exists in model already
    for (int i = 0; i < numMaterials; i++) {
	int tagI;
	if (OPS_GetIntInput(&numdata, &tagI) < 0) {
	    opserr << "WARNING invalid component tag\n";
	    opserr << "section Parallel: " << tag << endln;
	    return 0;
	}
	    
	SectionForceDeformation *theMat = OPS_getSectionForceDeformation(tagI);
	    
	if (theMat == 0) {
	    opserr << "WARNING component section does not exist\n";
	    opserr << "Component section: "; 
	    opserr << "\tsection Parallel: " << tag << endln;
	    delete [] theMats;
	    return 0;
	}
	else
	    theMats[i] = theMat;
    }	
	
    // Parsing was successful, allocate the material
    SectionForceDeformation* theSection = new ParallelSection(tag, numMaterials, theMats);
	
    // Deallocate the temporary pointers
    delete [] theMats;

    return theSection;
}

// constructors:
ParallelSection::ParallelSection (int tag, int numSecs,
				  SectionForceDeformation **theSecs): 
  SectionForceDeformation(tag, SEC_TAG_Parallel), 
  numSections(numSecs), theSections(0), 
  e(0), s(0), ks(0), fs(0), order(0), theCode(0),
  otherDbTag(0)
{
  if (theSecs == 0) {
    opserr << "ParallelSection::ParallelSection -- null section array passed\n";
    exit(-1);
  }
  
  theSections = new SectionForceDeformation *[numSections];
  if (theSections == 0) {
    opserr << "ParallelSection::ParallelSection -- failed to allocate pointers\n";
    exit(-1);
  }    

  for (int i = 0; i < numSections; i++) {
    if (theSecs[i] == 0) {
      opserr << "ParallelSection::ParallelSection -- null section pointer passed\n";
      exit(-1);
    }
    
    theSections[i] = theSecs[i]->getCopy();
    if (theSections[i] == 0) {
      opserr << "ParallelSection::ParallelSection -- failed to copy section\n";
      exit(-1);
    }
  }

  this->setUpCode();
}

void
ParallelSection::setUpCode(void)
{
  order = 0;

  bool haveP = false;
  bool haveMZ = false;
  bool haveMY = false;
  bool haveVZ = false;
  bool haveVY = false;
  bool haveT = false;
  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    const ID &codei = theSections[i]->getType();
    for (int j = 0; j < orderi; j++) {
      if (codei(j) == SECTION_RESPONSE_P)
	haveP = true;
      if (codei(j) == SECTION_RESPONSE_MZ)
	haveMZ = true;
      if (codei(j) == SECTION_RESPONSE_VY)
	haveVY = true;
      if (codei(j) == SECTION_RESPONSE_MY)
	haveMY = true;
      if (codei(j) == SECTION_RESPONSE_VZ)
	haveVZ = true;
      if (codei(j) == SECTION_RESPONSE_T)
	haveT = true;
    }
  }
  if (haveP)
    order++;
  if (haveMZ)
    order++;
  if (haveVY)
    order++;
  if (haveMY)
    order++;
  if (haveVZ)
    order++;
  if (haveT)
    order++;
  
  if (order > maxOrder) {
    opserr << "ParallelSection::ParallelSection -- order too big, need to modify the #define in ParallelSection.cpp to " <<
      order << endln;
    exit(-1);
  }

  theCode = new ID(codeArea, order);
  e = new Vector(workArea, order);
  s = new Vector(&workArea[maxOrder], order);
  ks = new Matrix(&workArea[2*maxOrder], order, order);
  fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);

  if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0) {
    opserr << "ParallelSection::ParallelSection -- out of memory\n";
    exit(-1);
  }

  int codeIndex = 0;
  if (haveP)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_P;
  if (haveMZ)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_MZ;
  if (haveVY)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_VY;
  if (haveMY)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_MY;
  if (haveVZ)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_VZ;
  if (haveT)
    (*theCode)(codeIndex++) = SECTION_RESPONSE_T;

}

// constructor for blank object that recvSelf needs to be invoked upon
ParallelSection::ParallelSection():
  SectionForceDeformation(0, SEC_TAG_Parallel),
  numSections(0), theSections(0),
  e(0), s(0), ks(0), fs(0), order(0), theCode(0),
  otherDbTag(0)
{

}

// destructor:
ParallelSection::~ParallelSection()
{
   for (int i = 0; i < numSections; i++)
       if (theSections[i])
	   delete theSections[i];

   if (theSections)
       delete [] theSections;
   
   if (e != 0)
     delete e;

   if (s != 0)
     delete s;

   if (ks != 0)
     delete ks;

   if (fs != 0)
     delete fs;

   if (theCode != 0)
     delete theCode;
}

int 
ParallelSection::setTrialSectionDeformation (const Vector &def)
{
  *e = def;

  int ret = 0;

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    Vector defi(orderi);
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  defi(j) = def(k);

    ret += theSections[i]->setTrialSectionDeformation(defi);
  }

  return ret;
}

const Vector &
ParallelSection::getSectionDeformation(void)
{
  return *e;
}

const Matrix &
ParallelSection::getSectionTangent(void)
{
  ks->Zero();

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    Matrix P(orderi,order); // The lazy man's approach -- MHS
    const ID &code = theSections[i]->getType();
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  P(j,k) = 1.0;
    const Matrix &ksi = theSections[i]->getSectionTangent();
    //opserr << P << ksi;
    ks->addMatrixTripleProduct(1.0, P, ksi, 1.0);
  }

  return *ks;
}

const Matrix &
ParallelSection::getInitialTangent(void)
{
  ks->Zero();

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    Matrix P(orderi,order);  // The lazy man's approach -- MHS
    const ID &code = theSections[i]->getType();
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  P(j,k) = 1.0;

    const Matrix &ksi = theSections[i]->getInitialTangent();

    ks->addMatrixTripleProduct(1.0, P, ksi, 1.0);
  }
  
  return *ks;
}

const Vector &
ParallelSection::getStressResultant(void)
{
  s->Zero();

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    const Vector &si = theSections[i]->getStressResultant();
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  (*s)(k) += si(j);
  }

  return *s;
}

SectionForceDeformation *
ParallelSection::getCopy(void)
{
  ParallelSection *theCopy = 0;
    
  theCopy = new ParallelSection(this->getTag(), numSections, theSections);
  
  if (theCopy == 0) {
    opserr << "ParallelSection::getCopy -- failed to allocate copy\n";
    exit(-1);
  }
			  
  return theCopy;
}

const ID&
ParallelSection::getType ()
{
  return *theCode;
}

int
ParallelSection::getOrder () const
{
  return order;
}

int
ParallelSection::commitState(void)
{
  int err = 0;
    
  for (int i = 0; i < numSections; i++)
    err += theSections[i]->commitState();
  
  return err;
}

int
ParallelSection::revertToLastCommit(void)
{
  int err = 0;
  
  for (int i = 0; i < numSections; i++)
    err += theSections[i]->revertToLastCommit();
  
  return err;
}	

int
ParallelSection::revertToStart(void)
{
  int err = 0;
  
  e->Zero();

  for (int i = 0; i < numSections; i++)
    err += theSections[i]->revertToStart();
  
  return err;
}

int
ParallelSection::sendSelf(int cTag, Channel &theChannel)
{
  static ID data(3); // 3 so no conflict when one section
  data(0) = this->getDbTag();
  data(1) = numSections;

  int dbTag = this->getDbTag();
  
  if (theChannel.sendID(dbTag, cTag, data) < 0) {
    opserr << "ParallelSection::sendSelf() - failed to send data" << endln;
    return -1;
  }  

  // now create an ID containing the class tags and dbTags of all
  // the MaterialModel objects in this ParallelMaterial
  // then send each of the MaterialModel objects
  ID classTags(numSections*2);
  for (int i=0; i<numSections; i++) {
    classTags(i) = theSections[i]->getClassTag();
    int matDbTag = theSections[i]->getDbTag();
    if (matDbTag == 0) {
      matDbTag  = theChannel.getDbTag();
      if (matDbTag != 0)
	theSections[i]->setDbTag(matDbTag);
    }
    classTags(i+numSections) = matDbTag;
  }
  
  if (theChannel.sendID(dbTag, cTag, classTags) < 0) {
    opserr << "ParallelSection::sendSelf() - failed to send classTags" << endln;
    return -3;
  }
  
  for (int j=0; j<numSections; j++)
    if (theSections[j]->sendSelf(cTag, theChannel) < 0) {
      opserr << "ParallelSection::sendSelf() - failed to send section" << endln;
      return -4;
    }
  
  return 0;
}


int
ParallelSection::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID data(3);
  int dbTag = this->getDbTag();
  
  if (theChannel.recvID(dbTag, cTag, data) < 0) {
    opserr << "ParallelSection::recvSelf() - failed to receive data" << endln;
    return -1;
  }
  
  this->setTag(int(data(0)));
  int numSectionsSent = int(data(1));
  if (numSections != numSectionsSent) { 
    numSections = numSectionsSent;
    if (theSections != 0) {
      for (int i=0; i<numSections; i++)
	delete theSections[i];
      
      delete [] theSections;
    }
    
    theSections = new SectionForceDeformation *[numSections];      
    if (theSections == 0) {
      opserr << "FATAL ParallelSection::recvSelf() - ran out of memory";
      opserr << " for array of size: " << numSections << endln;
      return -2;
    }
    for (int i=0; i<numSections; i++)
      theSections[i] = 0;
  }
  
  // create and receive an ID for the classTags and dbTags of the local 
  // MaterialModel objects
  ID classTags(numSections*2);
  if (theChannel.recvID(dbTag, cTag, classTags) < 0) {
    opserr << "ParallelSection::recvSelf() - failed to receive classTags" << endln;
    return -4;
  }
  
  // now for each of the MaterialModel objects, create a new object
  // and invoke recvSelf() on it
  for (int i=0; i<numSections; i++) {
    int matClassTag = classTags(i);
    if (theSections[i] == 0 || theSections[i]->getClassTag() != matClassTag) {
      if (theSections[i] == 0)
	delete theSections[i];
      SectionForceDeformation *theSectionModel = 
	theBroker.getNewSection(matClassTag);
      if (theSectionModel != 0) {
	theSections[i] = theSectionModel;
	theSectionModel->setDbTag(classTags(i+numSections));
      }
      else {
	opserr << "FATAL ParallelSection::recvSelf() ";
	opserr << " could not get a SectionForceDeformation" << endln;
	return -5;
      }    	    
    }
    if (theSections[i]->recvSelf(cTag, theChannel, theBroker) < 0) {
      opserr << "ParallelSection::recvSelf - failed to receive section" << endln;
      return -6;
    }
  }

  this->setUpCode();
  
  return 0;  
}

void
ParallelSection::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "\nSection Parallel, tag: " << this->getTag() << endln;
        if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
            for (int i = 0; i < numSections; i++) {
                s << "\t\tSection, tag: " << endln;
                theSections[i]->Print(s, flag);
            }
        }
        else {
            for (int i = 0; i < numSections; i++)
                s << "\t\tSection, tag: " << theSections[i]->getTag() << endln;
        }
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ParallelSection\", ";
        s << "\"sections\": [";
        for (int i = 0; i < numSections-1; i++)
            s << "\"" << theSections[i]->getTag() << "\", ";
        s << "\"" << theSections[numSections - 1]->getTag() << "\"]}";
    }
}

// AddingSensitivity:BEGIN ////////////////////////////////////
int
ParallelSection::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // Check if the parameter belongs to the material (only option for now)
  if (strstr(argv[0],"section") != 0) {
    
    if (argc < 3)
      return -1;

    // Get the tag of the material
    int sectionTag = atoi(argv[1]);
    
    // Loop to find the right material
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      if (sectionTag == theSections[i]->getTag())
	ok += theSections[i]->setParameter(&argv[2], argc-2, param);

    return ok;
  } 
  else {
    int ok = 0;
    for (int i = 0; i < numSections; i++)
      ok += theSections[i]->setParameter(argv, argc, param);

    return ok;
  }
}

const Vector &
ParallelSection::getSectionDeformationSensitivity(int gradIndex)
{
  s->Zero();

  return *s;
}

const Vector &
ParallelSection::getStressResultantSensitivity(int gradIndex,
					       bool conditional)
{
  s->Zero();

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    const Vector &si = theSections[i]->getStressResultantSensitivity(gradIndex, conditional);
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  (*s)(k) += si(j);
  }
  
  return *s;
}

const Matrix &
ParallelSection::getSectionTangentSensitivity(int gradIndex)
{
  ks->Zero();
  
  return *ks;
}

int
ParallelSection::commitSensitivity(const Vector& defSens,
				   int gradIndex, int numGrads)
{
  int ret = 0;

  dedh = defSens;

  for (int i = 0; i < numSections; i++) {
    int orderi = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();
    Vector defi(orderi);
    for (int j = 0; j < orderi; j++) 
      for (int k = 0; k < order; k++) 
	if (code(j) == (*theCode)(k))
	  defi(j) = defSens(k);

    ret += theSections[i]->commitSensitivity(defi, gradIndex, numGrads);
  }
  
  return ret;
}

// AddingSensitivity:END ///////////////////////////////////

const Vector&
ParallelSection::getdedh(void)
{
  return dedh;
}
