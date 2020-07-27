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
                                                                        
// $Revision: 1.23 $
// $Date: 2009-10-02 20:48:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/SectionAggregator.cpp,v $
                                                                        
                                                                        
// File: ~/section/SectionAggregator.C
//
// Written: MHS
// Created: June 2000
// Revision: A
//
// Purpose: This file contains the implementation for the SectionAggregator class. 

#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <SectionAggregator.h>
#include <MaterialResponse.h>
#include <ID.h>
#include <FiberSection2d.h>			//by SAJalali
#include <string.h>

#include <classTags.h>
#include <elementAPI.h>
#include <vector>

void* OPS_SectionAggregator()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section Aggregator tag? uniTag1? code1? ... <-section secTag?>" << endln;
	return 0;
    }
	    
    int tag;
    int secTag;
    SectionForceDeformation *theSec = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid Aggregator tag" << endln;
	return 0;
    }

    // uni mat tags and section dofs
    std::vector<UniaxialMaterial*> theMats;
    ID codes(0, 10);
    while (OPS_GetNumRemainingInputArgs() > 1) {
	int matTag;
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    OPS_ResetCurrentInputArg(-1);
	    break;
	}

	UniaxialMaterial* mat = OPS_getUniaxialMaterial(matTag);
	    
	if (mat == 0) {
	    opserr << "WARNING uniaxial material does not exist\n";
	    opserr << "uniaxial material: " << matTag; 
	    opserr << "\nsection Aggregator: " << tag << endln;
	    return 0;
	}

	theMats.push_back(mat);
	
	const char* type = OPS_GetString();
	int code = 0;
	if (strcmp(type,"Mz") == 0) 
	    code = SECTION_RESPONSE_MZ;
	else if (strcmp(type,"P") == 0)
	    code = SECTION_RESPONSE_P;
	else if (strcmp(type,"Vy") == 0)
	    code = SECTION_RESPONSE_VY;
	else if (strcmp(type,"My") == 0)
	    code = SECTION_RESPONSE_MY;
	else if (strcmp(type,"Vz") == 0)
	    code = SECTION_RESPONSE_VZ;
	else if (strcmp(type,"T") == 0)
	    code = SECTION_RESPONSE_T;
	else {
	    opserr << "WARNING invalid code" << endln;
	    opserr << "\nsection Aggregator: " << tag << endln;
	    return 0;
	}
	codes[codes.Size()] = code;
    }

    int nMats = (int)theMats.size();
    if (nMats == 0) {
      opserr << "No material is given\n";
      return 0;
    }

    // section
    if (OPS_GetNumRemainingInputArgs() > 1) {
	const char* flag = OPS_GetString();
	if (strcmp(flag, "-section") == 0) {
	    if (OPS_GetIntInput(&numdata, &secTag) < 0) {
		opserr << "WARNING invalid Aggregator section tag" << endln;
		return 0;
	    }
	    theSec = OPS_getSectionForceDeformation(secTag);
	    if (theSec == 0) {
		opserr << "WARNING section does not exist\n";
		opserr << "section: " << secTag; 
		opserr << "\nsection Aggregator: " << tag << endln;
		return 0;
	    }
	}
    }

    if (theSec) {
	return new SectionAggregator (tag, *theSec, nMats, &theMats[0], codes);
    } else {
	return new SectionAggregator (tag, nMats, &theMats[0], codes);
    }
	
    return 0;
}

void* OPS_UniaxialSection()
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 3) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: section Uniaxial tag? 1DTag? code?" << endln;
    return 0;
  }
  
  int data[2];
  numdata = 2;
  if (OPS_GetIntInput(&numdata, data) < 0) {
    opserr << "WARNING invalid integer" << endln;
    return 0;
  }
  
  int code;
  const char* type = OPS_GetString();
  if (strcmp(type,"Mz") == 0)
    code = SECTION_RESPONSE_MZ;
  else if (strcmp(type,"P") == 0)
    code = SECTION_RESPONSE_P;
  else if (strcmp(type,"Vy") == 0)
    code = SECTION_RESPONSE_VY;
  else if (strcmp(type,"My") == 0)
    code = SECTION_RESPONSE_MY;
  else if (strcmp(type,"Vz") == 0)
    code = SECTION_RESPONSE_VZ;
  else if (strcmp(type,"T") == 0)
    code = SECTION_RESPONSE_T;
  else {
    opserr << "WARNING invalid code" << endln;
    opserr << "Uniaxial section: " << data[0] << endln;
    return 0;
  }
  
  // Retrieve the uniaxial material from the model builder
  UniaxialMaterial *theMat = OPS_getUniaxialMaterial(data[1]);
  
  if (theMat == 0) {
    opserr << "WARNING uniaxial material does not exist\n";
    opserr << "uniaxial material: " << data[0];
    opserr << "\nUniaxial section: " << data[1] << endln;
    return 0;
  }
  
  // Parsing was successful, allocate the section
  //theSection = new GenericSection1d (tag, *theMat, code);
  
  UniaxialMaterial *theMats[1];
  theMats[0] = theMat;
  ID codeID(1);
  codeID(0) = code;
  return new SectionAggregator(data[0], 1, theMats, codeID);
}

#define maxOrder 10

// Assumes section order is less than or equal to maxOrder.
// Can increase if needed!!!
double SectionAggregator::workArea[2*maxOrder*(maxOrder+1)];
int    SectionAggregator::codeArea[maxOrder];

// constructors:
SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      int numAdds, UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  e(0), s(0), ks(0), fs(0), theCode(0),
  otherDbTag(0)
{
    theSection = theSec.getCopy();
    
    if (!theSection) {
      opserr << "SectionAggregator::SectionAggregator " << tag << " -- failed to get copy of section\n";
      exit(-1);
    }

    if (!theAdds) {
      opserr << "SectionAggregator::SectionAggregator " << tag << " -- null uniaxial material array passed\n";
      exit(-1);
    }
    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions) {
      opserr << "SectionAggregator::SectionAggregator " << tag << "  -- failed to allocate pointers\n";
      exit(-1);
    }    
    int i;
    
    for (i = 0; i < numMats; i++) {
      if (!theAdds[i]) {
	opserr << "SectionAggregator::SectionAggregator " << tag << " -- null uniaxial material pointer passed\n";
	exit(-1);
      }	
      //theAdditions[i] = theAdds[i]->getCopy(this);
      theAdditions[i] = theAdds[i]->getCopy();
      
      if (!theAdditions[i]) {
	opserr << "SectionAggregator::SectionAggregator " << tag << " -- failed to copy uniaxial material\n";
	opserr << theAdds[i];
	exit(-1);
      }
    }

    int order = theSec.getOrder()+numAdds;

    if (order > maxOrder) {
      opserr << "SectionAggregator::SectionAggregator   " << tag << "  -- order too big, need to modify the #define in SectionAggregator.cpp to " <<
	order << endln;
      exit(-1);
    }

    theCode = new ID(codeArea, order);
    e = new Vector(workArea, order);
    s = new Vector(&workArea[maxOrder], order);
    ks = new Matrix(&workArea[2*maxOrder], order, order);
    fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
    matCodes = new ID(addCodes);

    if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
      opserr << "SectionAggregator::SectionAggregator   " << tag << " -- out of memory\n";
      exit(-1);
    }	
}

SectionAggregator::SectionAggregator (int tag, int numAdds,
				      UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), matCodes(0), numMats(numAdds),
  e(0), s(0), ks(0), fs(0), theCode(0),
  otherDbTag(0)
{
  if (!theAdds) {
    opserr << "SectionAggregator::SectionAggregator  " << tag << " -- null uniaxial material array passed\n";
    exit(-1);
  }

  theAdditions = new UniaxialMaterial *[numMats];

  if (!theAdditions) {
    opserr << "SectionAggregator::SectionAggregator   " << tag << " -- failed to allocate pointers\n";
    exit(-1);
  }    
    int i;
    
    for (i = 0; i < numMats; i++) {
      if (!theAdds[i]) {
	opserr << "SectionAggregator::SectionAggregator   " << tag << " -- null uniaxial material pointer passed\n";
	exit(-1);
      }		
	
      //      theAdditions[i] = theAdds[i]->getCopy(this);
      theAdditions[i] = theAdds[i]->getCopy();
      
      if (!theAdditions[i]) {
	opserr << "SectionAggregator::SectionAggregator   " << tag << " -- failed to copy uniaxial material\n";
	opserr << theAdds[i];
	exit(-1);
      }
    }

    int order = numAdds;

    if (order > maxOrder) {
      opserr << "SectionAggregator::SectionAggregator   " << tag << " -- order too big, need to modify the #define in SectionAggregator.cpp to %d\n";
      exit(-1);
    }

    theCode = new ID(codeArea, order);
    e = new Vector(workArea, order);
    s = new Vector(&workArea[maxOrder], order);
    ks = new Matrix(&workArea[2*maxOrder], order, order);
    fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
    matCodes = new ID(addCodes);

    if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
      opserr << "SectionAggregator::SectionAggregator   " << tag << " -- out of memory\n";
      exit(-1);
    }
}

SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      UniaxialMaterial &theAddition, int c) :
  SectionForceDeformation(tag, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(1),
  e(0), s(0), ks(0), fs(0), theCode(0),
  otherDbTag(0)
{
  theSection = theSec.getCopy();

  if (!theSection) {
    opserr << "SectionAggregator::SectionAggregator   " << tag << " -- failed to get copy of section\n";
    exit(-1);
  }

  theAdditions = new UniaxialMaterial *[1];
  
   theAdditions[0] = theAddition.getCopy();
  
  if (!theAdditions[0]) {
    opserr << "SectionAggregator::SectionAggregator   " << tag << " -- failed to copy uniaxial material\n";
    exit(-1);
  }
    
  matCodes = new ID(1);
  (*matCodes)(0) = c;
  
  int order = theSec.getOrder()+1;
  
  if (order > maxOrder) {
    opserr << "SectionAggregator::SectionAggregator   " << tag << " -- order too big, need to modify the #define in SectionAggregator.cpp to %d\n";
    exit(-1);
  }
  
  theCode = new ID(codeArea, order);
  e = new Vector(workArea, order);
  s = new Vector(&workArea[maxOrder], order);
  ks = new Matrix(&workArea[2*maxOrder], order, order);
  fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
  
  if (theCode == 0 || e == 0 || s == 0 || ks == 0 || fs == 0 || matCodes == 0) {
    opserr << "SectionAggregator::SectionAggregator   " << tag << " -- out of memory\n";
    exit(-1);
  }
}

// constructor for blank object that recvSelf needs to be invoked upon
SectionAggregator::SectionAggregator():
  SectionForceDeformation(0, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), matCodes(0), numMats(0), 
  e(0), s(0), ks(0), fs(0), theCode(0),
  otherDbTag(0)
{

}

// destructor:
SectionAggregator::~SectionAggregator()
{
   if (theSection)
       delete theSection;

   for (int i = 0; i < numMats; i++)
       if (theAdditions[i])
	   delete theAdditions[i];

   if (theAdditions)
       delete [] theAdditions;
   
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

   if (matCodes != 0)
     delete matCodes;
}

int SectionAggregator::setTrialSectionDeformation (const Vector &def)
{
  int ret = 0;
  int i = 0;

  int theSectionOrder = 0;

  if (theSection) {
    theSectionOrder = theSection->getOrder();
    Vector v(workArea, theSectionOrder);
    
    for (i = 0; i < theSectionOrder; i++)
      v(i) = def(i);
    
    ret = theSection->setTrialSectionDeformation(v);
  }

  int order = theSectionOrder + numMats;
  
  for ( ; i < order; i++)
    ret += theAdditions[i-theSectionOrder]->setTrialStrain(def(i));
  
  return ret;
}

const Vector &
SectionAggregator::getSectionDeformation(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &eSec = theSection->getSectionDeformation();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*e)(i) = eSec(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*e)(i) = theAdditions[i-theSectionOrder]->getStrain();

  return *e;
}

const Matrix &
SectionAggregator::getSectionTangent(void)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getSectionTangent();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getTangent();
  
  return *ks;
}

const Matrix &
SectionAggregator::getInitialTangent(void)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getInitialTangent();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getInitialTangent();
  
  return *ks;
}

const Matrix &
SectionAggregator::getSectionFlexibility(void)
{
  int i = 0;
    
  int theSectionOrder = 0;

  // Zero before assembly
  fs->Zero();

  if (theSection) {
    const Matrix &fSec = theSection->getSectionFlexibility();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*fs)(i,j) = fSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++) {
    double k = theAdditions[i-theSectionOrder]->getTangent();
    if (k == 0.0) {
      opserr << "SectionAggregator::getSectionFlexibility -- singular section stiffness\n";
			       
      (*fs)(i,i) = 1.e14;
    }
    else
      (*fs)(i,i) = 1/k;
  }	
  
  return *fs;
}

const Matrix &
SectionAggregator::getInitialFlexibility(void)
{
  int i = 0;
    
  int theSectionOrder = 0;

  // Zero before assembly
  fs->Zero();

  if (theSection) {
    const Matrix &fSec = theSection->getInitialFlexibility();
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*fs)(i,j) = fSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++) {
    double k = theAdditions[i-theSectionOrder]->getInitialTangent();
    (*fs)(i,i) = 1.0/k;
  }	
  
  return *fs;
}

const Vector &
SectionAggregator::getStressResultant(void)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &sSec = theSection->getStressResultant();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*s)(i) = sSec(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*s)(i) = theAdditions[i-theSectionOrder]->getStress();
  
  return *s;
}

SectionForceDeformation *
SectionAggregator::getCopy(void)
{
  SectionAggregator *theCopy = 0;
    
  if (theSection)
    theCopy = new SectionAggregator(this->getTag(), *theSection,
				    numMats, theAdditions, *matCodes);
  else
    theCopy = new SectionAggregator(this->getTag(), numMats,
				    theAdditions, *matCodes);
  
  if (theCopy == 0) {
    opserr << "SectionAggregator::getCopy -- failed to allocate copy\n";
    exit(-1);
  }
			  
  
  return theCopy;
}

const ID&
SectionAggregator::getType ()
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const ID &secType = theSection->getType();
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*theCode)(i) = secType(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*theCode)(i) = (*matCodes)(i-theSectionOrder);

  return *theCode;
}

int
SectionAggregator::getOrder () const
{
  int order = numMats;

  if (theSection != 0)
    order += theSection->getOrder();

  return order;
}

int
SectionAggregator::commitState(void)
{
  int err = 0;
    
  if (theSection)
    err += theSection->commitState();
  
  for (int i = 0; i < numMats; i++)
    err += theAdditions[i]->commitState();
  
  return err;
}

int
SectionAggregator::revertToLastCommit(void)
{
  int err = 0;
  
  int i = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToLastCommit();
  
  // Do the same for the uniaxial materials
  for (i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToLastCommit();
  
  return err;
}	

int
SectionAggregator::revertToStart(void)
{
  int err = 0;
  
  // Revert the section
  if (theSection)
    err += theSection->revertToStart();
  
  // Do the same for the uniaxial materials
  for (int i = 0; i < numMats; i++)
    err += theAdditions[i]->revertToStart();
  
  return err;
}


int
SectionAggregator::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  // Need otherDbTag since classTags ID and data ID may be the same size
  if (otherDbTag == 0) 
    otherDbTag = theChannel.getDbTag();
  
  // Create ID for tag and section order data
  static ID data(5);
  
  int order = this->getOrder();
  
  data(0) = this->getTag();
  data(1) = otherDbTag;
  data(2) = order;
  if (theSection == 0)
	  data(3) = 0;
  else 
	  data(3) = theSection->getOrder();
  data(3) = (theSection != 0) ? theSection->getOrder() : 0;
  data(4) = numMats;

  // Send the tag and section order data
  res += theChannel.sendID(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "SectionAggregator::sendSelf -- could not send data ID\n";
			    
    return res;
  }
  
  // Determine how many classTags there are and allocate ID vector
  // for the tags and section code
  int numTags = (theSection == 0) ? numMats : numMats + 1;
  ID classTags(2*numTags + numMats);
  
  // Loop over the UniaxialMaterials filling in class and db tags
  int i, dbTag;
  for (i = 0; i < numMats; i++) {
    classTags(i) = theAdditions[i]->getClassTag();
    
    dbTag = theAdditions[i]->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theAdditions[i]->setDbTag(dbTag);
    }
    
    classTags(i+numTags) = dbTag;
  }
  
  // Put the Section class and db tags into the ID vector
  if (theSection != 0) {
    classTags(numTags-1) = theSection->getClassTag();
    
    dbTag = theSection->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theSection->setDbTag(dbTag);
    }
    
    classTags(2*numTags-1) = dbTag;
  }
  
  // Put the UniaxialMaterial codes into the ID vector
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    classTags(j) = (*matCodes)(i);
  
  // Send the material class and db tags and section code
  res += theChannel.sendID(otherDbTag, cTag, classTags);
  if (res < 0) {
    opserr << "SectionAggregator::sendSelf -- could not send classTags ID\n";
    return res;
  }

  // Ask the UniaxialMaterials to send themselves
  for (i = 0; i < numMats; i++) {
    res += theAdditions[i]->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "SectionAggregator::sendSelf -- could not send UniaxialMaterial, i = " << i << endln;
      return res;
    }
  }
  
  // Ask the Section to send itself
  if (theSection != 0) {
    res += theSection->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "SectionAggregator::sendSelf -- could not send SectionForceDeformation\n";
      return res;
    }
  }
  
  return res;
}


int
SectionAggregator::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // Create an ID and receive tag and section order
  static ID data(5);
  res += theChannel.recvID(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "SectionAggregator::recvSelf -- could not receive data ID\n";
    return res;
  }
  
  this->setTag(data(0));
  otherDbTag = data(1);
  int order = data(2);
  int theSectionOrder = data(3);
  numMats = data(4);

  if (order > 0) {
    if (e == 0 || e->Size() != order) {
      if (e != 0) {
	delete e;
	delete s;
	delete ks;
	delete fs;
	delete theCode;
      }
      e = new Vector(workArea, order);
      s = new Vector(&workArea[maxOrder], order);
      ks = new Matrix(&workArea[2*maxOrder], order, order);
      fs = new Matrix(&workArea[maxOrder*(maxOrder+2)], order, order);
      theCode = new ID(codeArea, order);
    }
  }

  if (numMats > 0) {
    if (matCodes == 0 || matCodes->Size() != numMats) {
      if (matCodes != 0)
	delete matCodes;

      matCodes = new ID(numMats);
    }
  }

  // Determine how many classTags there are and allocate ID vector
  int numTags = (theSectionOrder == 0) ? numMats : numMats + 1;
  ID classTags(numTags*2 + numMats);
  
  // Receive the material class and db tags
  res += theChannel.recvID(otherDbTag, cTag, classTags);
  if (res < 0) {
    opserr << "SectionAggregator::recvSelf -- could not receive classTags ID\n";
    return res;
  }

  // Check if null pointer, allocate if so
  if (theAdditions == 0) {
    theAdditions = new UniaxialMaterial *[numMats];
    if (theAdditions == 0) {
      opserr << "SectionAggregator::recvSelf -- could not allocate UniaxialMaterial array\n";
      return -1;
    }
    // Set pointers to null ... will get allocated by theBroker
    for (int j = 0; j < numMats; j++)
      theAdditions[j] = 0;
  }
  
  // Loop over the UniaxialMaterials
  int i, classTag;
  for (i = 0; i < numMats; i++) {
    classTag = classTags(i);
    
    // Check if the UniaxialMaterial is null; if so, get a new one
    if (theAdditions[i] == 0)
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    
    // Check that the UniaxialMaterial is of the right type; if not, delete
    // the current one and get a new one of the right type
    else if (theAdditions[i]->getClassTag() != classTag) {
      delete theAdditions[i];
      theAdditions[i] = theBroker.getNewUniaxialMaterial(classTag);
    }
    
    // Check if either allocation failed
    if (theAdditions[i] == 0) {
      opserr << "SectionAggregator::recvSelf -- could not get UniaxialMaterial, i = " << i << endln;
      return -1;
    }
    
    // Now, receive the UniaxialMaterial
    theAdditions[i]->setDbTag(classTags(i+numTags));
    res += theAdditions[i]->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "SectionAggregator::recvSelf -- could not receive UniaxialMaterial, i = " << i << endln;
      return res;
    }
  }

  // If there is no Section to receive, return
  if (theSectionOrder != 0) {
	
 
	classTag = classTags(numTags-1);
  
	// Check if the Section is null; if so, get a new one
	if (theSection == 0)
		theSection = theBroker.getNewSection(classTag);
  
	// Check that the Section is of the right type; if not, delete
	// the current one and get a new one of the right type
	else if (theSection->getClassTag() != classTag) {
		delete theSection;
		theSection = theBroker.getNewSection(classTag);
	}
  
	// Check if either allocation failed
	if (theSection == 0) {
		opserr << "SectionAggregator::recvSelf -- could not get a SectionForceDeformation\n";
		return -1;
	}

	// Now, receive the Section
	theSection->setDbTag(classTags(2*numTags-1));
	res += theSection->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		opserr << "SectionAggregator::recvSelf -- could not receive SectionForceDeformation\n";
		return res;
	}
  }

  // Fill in the section code
  int j = 2*numTags;
  for (i = 0; i < numMats; i++, j++)
    (*matCodes)(i) = classTags(j);

  return res;
}

Response*
SectionAggregator::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  
  Response *theResponse =0;

  if ( (strcmp(argv[0],"deformations") == 0) || (strcmp(argv[0],"deformation") == 0) || 
	(strcmp(argv[0],"forces") == 0) || (strcmp(argv[0],"force") == 0) ||
       (strcmp(argv[0],"forceAndDeformation") == 0)) {
    
    return this->SectionForceDeformation::setResponse(argv, argc, output);
  } 
  // by SAJalali
  int num = numMats;
  if (theSection != 0)
	  num++;
  if ((strcmp(argv[0], "energy") == 0) || (strcmp(argv[0], "Energy") == 0)) {
	  return theResponse = new MaterialResponse(this, 8, Vector(num));
  }


  if (theSection != 0)
    return theSection->setResponse(argv, argc, output);
  else 
    return this->SectionForceDeformation::setResponse(argv, argc, output);

  return 0;
}

//by SAJalali
int
SectionAggregator::getResponse(int responseID, Information &sectInfo)
{
	FiberSection2d* sec = 0;
	int num = numMats;
	if (theSection != 0)
		num++;
	Vector res(num);
	if (responseID == 8) {
		for (int i = 0; i < numMats; i++)
			res(i) = theAdditions[i]->getEnergy();
		sec = (FiberSection2d*)theSection;
		if (sec != 0)
			res(numMats) = sec->getEnergy();
		sectInfo.setVector(res);
		//opserr << "energyVect=" << res << "\n";
	}
	else
		return SectionForceDeformation::getResponse(responseID, sectInfo);

	return 0;
}


void
SectionAggregator::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "\nSection Aggregator, tag: " << this->getTag() << endln;
        if (theSection) {
            s << "\tSection, tag: " << theSection->getTag() << endln;
            theSection->Print(s, flag);
        }
        s << "\tUniaxial Additions" << endln;
        for (int i = 0; i < numMats; i++)
            s << "\t\tUniaxial Material, tag: " << theAdditions[i]->getTag() << endln;
        s << "\tUniaxial codes " << *matCodes << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        theSection->Print(s, flag);
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"SectionAggregator\", ";
        if (theSection) {
            s << "\"section\": \"" << theSection->getTag() << "\", ";
        }
        s << "\"materials\": [";
        for (int i = 0; i < numMats - 1; i++)
            s << "\"" << theAdditions[i]->getTag() << "\", ";
        s << "\"" << theAdditions[numMats - 1]->getTag() << "\"], ";
        s << "\"dof\": [";
        for (int i = 0; i < numMats - 1; i++) {
            if ((*matCodes)(i) == 2)
                s << "\"P\", ";
            else if ((*matCodes)(i) == 3)
                s << "\"Vy\", ";
            else if ((*matCodes)(i) == 5)
                s << "\"Vz\", ";
            else if ((*matCodes)(i) == 6)
                s << "\"T\", ";
            else if ((*matCodes)(i) == 4)
                s << "\"My\", ";
            else if ((*matCodes)(i) == 1)
                s << "\"Mz\", ";
        }
        if ((*matCodes)(numMats - 1) == 2)
            s << "\"P\"]}";
        else if ((*matCodes)(numMats - 1) == 3)
            s << "\"Vy\"]}";
        else if ((*matCodes)(numMats - 1) == 5)
            s << "\"Vz\"]}";
        else if ((*matCodes)(numMats - 1) == 6)
            s << "\"T\"]}";
        else if ((*matCodes)(numMats - 1) == 4)
            s << "\"My\"]}";
        else if ((*matCodes)(numMats - 1) == 1)
            s << "\"Mz\"]}";
    }
}


int
SectionAggregator::getVariable(const char *argv, Information &info)
{

  info.theDouble = 0.0;
  int i;
  int order = numMats;
  if (theSection != 0)
    order += theSection->getOrder();

  const Vector &e = this->getSectionDeformation();
  const ID &code  = this->getType();

  if (strcmp(argv,"axialStrain") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_P)
	info.theDouble  += e(i);
  }  else if (strcmp(argv,"curvatureZ") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MZ)
	info.theDouble += e(i);
  } else if (strcmp(argv,"curvatureY") == 0) {
    for (i = 0; i < order; i++)
      if (code(i) == SECTION_RESPONSE_MY)
	info.theDouble += e(i);
  } else 
    return -1;

  return 0;
}

// AddingSensitivity:BEGIN ////////////////////////////////////
int
SectionAggregator::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // Check if the parameter belongs to only an aggregated material
  if (strstr(argv[0],"addition") != 0 || strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return -1;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop to find the right material
    int ok = 0;
    for (int i = 0; i < numMats; i++)
      if (materialTag == theAdditions[i]->getTag()) {
	ok = theAdditions[i]->setParameter(&argv[2], argc-2, param);
	if (ok != -1)
	  result = ok;
      }
    
    return result;
  } 
  
  // Check if the parameter belongs to only the section
  else if (strstr(argv[0],"section") != 0) {
    
    if (argc < 2) {
      opserr << "SectionAggregator::setParameter() - insufficient argc < 2 for section option. " << endln;
      return -1;
    }

    return theSection->setParameter(&argv[1], argc-1, param);
  } 

  else { // Default -- send to everything
    int ok = 0;
    for (int i = 0; i < numMats; i++) {
      ok = theAdditions[i]->setParameter(argv, argc, param);
      if (ok != -1)
	result = ok;
    }
    if (theSection != 0) {
      ok = theSection->setParameter(argv, argc, param);
      if (ok != -1)
	result = ok;
    }
  }

  return result;
}

const Vector &
SectionAggregator::getSectionDeformationSensitivity(int gradIndex)
{
  s->Zero();

  return *s;
}

const Vector &
SectionAggregator::getStressResultantSensitivity(int gradIndex,
						 bool conditional)
{
  int i = 0;

  int theSectionOrder = 0;
    
  if (theSection) {
    const Vector &dsdh = theSection->getStressResultantSensitivity(gradIndex,
								   conditional);
    theSectionOrder = theSection->getOrder();
    
    for (i = 0; i < theSectionOrder; i++)
      (*s)(i) = dsdh(i);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*s)(i) = theAdditions[i-theSectionOrder]->getStressSensitivity(gradIndex,
								    conditional);
  
  return *s;
}

const Matrix &
SectionAggregator::getSectionTangentSensitivity(int gradIndex)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getSectionTangentSensitivity(gradIndex);
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getTangentSensitivity(gradIndex);
  
  return *ks;
}

const Matrix&
SectionAggregator::getInitialTangentSensitivity(int gradIndex)
{
  int i = 0;

  int theSectionOrder = 0;

  // Zero before assembly
  ks->Zero();

  if (theSection) {
    const Matrix &kSec = theSection->getInitialTangentSensitivity(gradIndex);
    theSectionOrder = theSection->getOrder();

    for (i = 0; i < theSectionOrder; i++)
      for (int j = 0; j < theSectionOrder; j++)
	(*ks)(i,j) = kSec(i,j);
  }
  
  int order = theSectionOrder + numMats;

  for ( ; i < order; i++)
    (*ks)(i,i) = theAdditions[i-theSectionOrder]->getInitialTangentSensitivity(gradIndex);
  
  return *ks;
}

int
SectionAggregator::commitSensitivity(const Vector& defSens,
				     int gradIndex, int numGrads)
{
  int ret = 0;
  int i = 0;

  dedh = defSens;

  int theSectionOrder = 0;

  if (theSection) {
    theSectionOrder = theSection->getOrder();
    Vector dedh(workArea, theSectionOrder);
    
    for (i = 0; i < theSectionOrder; i++)
      dedh(i) = defSens(i);
    
    ret = theSection->commitSensitivity(dedh, gradIndex, numGrads);
  }

  int order = theSectionOrder + numMats;
  
  for ( ; i < order; i++)
    ret += theAdditions[i-theSectionOrder]->commitSensitivity(defSens(i),
							      gradIndex,
							      numGrads);
  
  return ret;
}

// AddingSensitivity:END ///////////////////////////////////

const Vector&
SectionAggregator::getdedh(void)
{
  return dedh;
}
