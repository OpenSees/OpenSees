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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
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
#include <ID.h>

#include <classTags.h>

// constructors:
SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      int numAdds, UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), numMats(numAdds),
  theSection(0), theAdditions(0), code(0), otherDbTag(0)
{
    theSection = theSec.getCopy();
    
    if (!theSection)
	g3ErrorHandler->fatal("SectionAggregator - failed to get copy of section");

    theSectionOrder = theSection->getOrder();

    order = numMats + theSectionOrder;

    ks = Matrix(order,order);
    fs = Matrix(order,order);
    s = Vector(order);
    e = Vector(order);
    
    if (!theAdds)
	g3ErrorHandler->fatal("SectionAggregator - null uniaxial material array passed");

    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions)
	g3ErrorHandler->fatal("SectionAggregator - failed to allocate pointers");
    
    int i;
    
    for (i = 0; i < numMats; i++) {
	if (!theAdds[i])
	    g3ErrorHandler->fatal("SectionAggregator - null uniaxial material pointer passed");
	
	theAdditions[i] = theAdds[i]->getCopy(this);
	
	if (!theAdditions[i])
	    g3ErrorHandler->fatal("SectionAggregator - failed to copy uniaxial material");
    }	
    
    code = new ID(order);
    
    if (!code)
    	g3ErrorHandler->fatal("SectionAggregator - failed to allocate ID");
    
    const ID &secCode = theSection->getType();
    
    for (i = 0; i < theSectionOrder; i++)
	(*code)(i) = secCode(i);
    
    for ( ; i < order; i++)
	(*code)(i) = addCodes(i-theSectionOrder);
}

SectionAggregator::SectionAggregator (int tag, int numAdds, UniaxialMaterial **theAdds,
				      const ID &addCodes): 
  SectionForceDeformation(tag, SEC_TAG_Aggregator), 
  theSection(0), theAdditions(0), code(0), 
  theSectionOrder(0), numMats(numAdds), order(numAdds), otherDbTag(0)
{
    ks = Matrix(order,order);
    fs = Matrix(order,order);
    s = Vector(order);
    e = Vector(order);
    
    if (!theAdds)
	g3ErrorHandler->fatal("SectionAggregator - null uniaxial material array passed");

    theAdditions = new UniaxialMaterial *[numMats];

    if (!theAdditions)
	g3ErrorHandler->fatal("SectionAggregator - failed to allocate pointers");
    
    int i;
    
    for (i = 0; i < numMats; i++) {
	if (!theAdds[i])
	    g3ErrorHandler->fatal("SectionAggregator - null uniaxial material pointer passed");
	
	theAdditions[i] = theAdds[i]->getCopy(this);
	
	if (!theAdditions[i])
	    g3ErrorHandler->fatal("SectionAggregator - failed to copy uniaxial material");
    }	
    
    code = new ID(addCodes);
    
    if (!code)
    	g3ErrorHandler->fatal("SectionAggregator - failed to allocate ID");
}

SectionAggregator::SectionAggregator (int tag, SectionForceDeformation &theSec,
				      UniaxialMaterial &theAddition, int c) :
  SectionForceDeformation(tag, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), code(0), numMats(1), otherDbTag(0)
{
    theSection = theSec.getCopy();
    
    if (!theSection)
	g3ErrorHandler->fatal("SectionAggregator - failed to get copy of section");

    theSectionOrder = theSection->getOrder();

    order = 1 + theSectionOrder;

    ks = Matrix(order,order);
    fs = Matrix(order,order);
    s = Vector(order);
    e = Vector(order);

    theAdditions = new UniaxialMaterial *[1];

    theAdditions[0] = theAddition.getCopy(this);
    
    if (!theAdditions[0])
	g3ErrorHandler->fatal("SectionAggregator - failed to copy uniaxial material");
    
    code = new ID(order);
    
    if (!code)
    	g3ErrorHandler->fatal("SectionAggregator - failed to allocate ID");
    
    const ID &secCode = theSection->getType();
    int i;
    
    for (i = 0; i < theSectionOrder; i++)
		(*code)(i) = secCode(i);
    
    (*code)(i) = c;
}



// constructor for blank object that recvSelf needs to be invoked upon
SectionAggregator::SectionAggregator():
  SectionForceDeformation(0, SEC_TAG_Aggregator),
  theSection(0), theAdditions(0), code(0), numMats(0),
  order(0), theSectionOrder(0), otherDbTag(0)
{

}

// destructor:
SectionAggregator::~SectionAggregator()
{
   int i;

   if (theSection)
       delete theSection;

   for (i = 0; i < numMats; i++)
       if (theAdditions[i])
	   delete theAdditions[i];

   if (theAdditions)
       delete [] theAdditions;
   
   if (code)
       delete code;
}

int SectionAggregator::setTrialSectionDeformation (const Vector &deforms)
{
	e = deforms;

	int ret = 0;
	int i = 0;

    if (theSection) {
		Vector v(theSectionOrder);
    
		for (i = 0; i < theSectionOrder; i++)
			v(i) = e(i);

		ret = theSection->setTrialSectionDeformation(v);
	}
	
	for ( ; i < order; i++)
		ret += theAdditions[i-theSectionOrder]->setTrialStrain(e(i));

	return ret;
}

const Vector &
SectionAggregator::getSectionDeformation(void)
{
   return e;
}

const Matrix &
SectionAggregator::getSectionTangent(void)
{
    int i = 0;

    if (theSection) {

		const Matrix &k = theSection->getSectionTangent();

		int j;
	
		for (i = 0; i < theSectionOrder; i++)
			for (j = 0; j < theSectionOrder; j++)
				ks(i,j) = k(i,j);
    }
    
    for ( ; i < order; i++) {
		ks(i,i) = theAdditions[i-theSectionOrder]->getTangent();
		if (ks(i,i) == 0.0) {
			g3ErrorHandler->warning("WARNING SectionAggregator::getSectionTangent - singular section stiffness");
			ks(i,i) = 1.e-12;
		}
	}

    return ks;
}

const Matrix &
SectionAggregator::getSectionFlexibility(void)
{
    int i = 0;
    
    if (theSection) {
    
		const Matrix &f = theSection->getSectionFlexibility();
    
		int j;
	
		for (i = 0; i < theSectionOrder; i++)
			for (j = 0; j < theSectionOrder; j++)
				fs(i,j) = f(i,j);
    }
    
    double k;
    
    for ( ; i < order; i++) {
		k = theAdditions[i-theSectionOrder]->getTangent();
		if (k == 0.0) {
			g3ErrorHandler->warning("WARNING SectionAggregator::getSectionFlexibility - singular section stiffness");
		    fs(i,i) = 1.e14;
		}
		else
			fs(i,i) = 1/k;
    }	
    
    return fs;
}

const Vector &
SectionAggregator::getStressResultant(void)
{
    int i = 0;
    
    if (theSection) {

		const Vector &p = theSection->getStressResultant();

		for (i = 0; i < theSectionOrder; i++)
			s(i) = p(i);
	}
	
    for ( ; i < order; i++)
		s(i) = theAdditions[i-theSectionOrder]->getStress();

    return s;
}

SectionForceDeformation *
SectionAggregator::getCopy(void)
{
    SectionAggregator *theCopy;
    
    if (theSection) {
	ID *c = new ID(numMats);

	for (int i = 0; i < numMats; i++)
	    (*c)(i) = (*code)(i+theSectionOrder);
    
	theCopy = new SectionAggregator (this->getTag(), *theSection, numMats, theAdditions, *c);
	
	delete c;
    }
    else
	theCopy = new SectionAggregator (this->getTag(), order, theAdditions, *code);
    
    if (!theCopy)
	g3ErrorHandler->fatal("SectionAggregator::getCopy - failed to allocate copy");
    
    theCopy->e = e;  // section deformations
    theCopy->s = s;  // section forces
    theCopy->ks = ks;  // section stiffness
    theCopy->fs = fs;  // section flexibility

    return theCopy;
}

const ID&
SectionAggregator::getType () const
{
    return *code;
}

int
SectionAggregator::getOrder () const
{
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

	// Revert the section, then reconstruct the section deformation vector
	// to its last committed value
    if (theSection) {
		err += theSection->revertToLastCommit();
		const Vector &esec = theSection->getSectionDeformation();
		
		for (i = 0; i < theSectionOrder; i++)
			e(i) = esec(i);
	}
	
	// Do the same for the uniaxial materials
    for ( ; i < order; i++) {
		int j = i-theSectionOrder;
		err += theAdditions[j]->revertToLastCommit();
		e(i) = theAdditions[j]->getStrain();
	}

	return err;
}	

int
SectionAggregator::revertToStart(void)
{
    int err = 0;
    
    if (theSection)
	err += theSection->revertToStart();
    
    for (int i = 0; i < numMats; i++)
	err += theAdditions[i]->revertToStart();
    
    // revert the section variables to start
    e.Zero();
    s.Zero();
    ks.Zero();
    fs.Zero(); 
    
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

	data(0) = this->getTag();
	data(1) = otherDbTag;
	data(2) = order;
	data(3) = theSectionOrder;
	data(4) = numMats;

	// Send the tag and section order data
	res += theChannel.sendID(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send data ID",
			"SectionAggregator::sendSelf");
		return res;
	}

	if (order > 0) {
		// Send the committed section deformations
		res += theChannel.sendVector(this->getDbTag(), cTag, e);
		if (res < 0) {
			g3ErrorHandler->warning("%s -- failed to send section deformations",
				"SectionAggregator::sendSelf");
			return res;
		}
	}

	// Determine how many classTags there are and allocate ID vector
	// for the tags and section code
	int numTags = (theSection == 0) ? numMats : numMats + 1;
	ID classTags(2*numTags + order);

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

	// Put the section code into the ID vector
	int j = 2*numTags;
	for (i = 0; i < order; i++, j++)
		classTags(j) = (*code)(i);

	// Send the material class and db tags and section code
	res += theChannel.sendID(otherDbTag, cTag, classTags);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send classTags ID",
			"SectionAggregator::sendSelf");
		return res;
	}

	// Ask the UniaxialMaterials to send themselves
	for (i = 0; i < numMats; i++) {
		res += theAdditions[i]->sendSelf(cTag, theChannel);
		if (res < 0) {
			g3ErrorHandler->warning("%s -- could not send UniaxialMaterial, i = %d",
				"SectionAggregator::sendSelf", i);
			return res;
		}
	}

	// Ask the Section to send itself
	if (theSection != 0) {
		res += theSection->sendSelf(cTag, theChannel);
		if (res < 0) {
			g3ErrorHandler->warning("%s -- could not send SectionForceDeformation",
				"SectionAggregator::sendSelf");
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
		g3ErrorHandler->warning("%s -- could not receive data ID",
			"SectionAggregator::recvSelf");
		return res;
	}
    
	this->setTag(data(0));
	otherDbTag = data(1);
	order = data(2);
	theSectionOrder = data(3);
	numMats = data(4);

	if (order > 0) {
		
		ks = Matrix(order,order);
	    fs = Matrix(order,order);
		s = Vector(order);
		e = Vector(order);

		// Receive the committed section deformations
		res += theChannel.recvVector(this->getDbTag(), cTag, e);
		if (res < 0) {
			g3ErrorHandler->warning("%s -- failed to receive section deformations",
				"SectionAggregator::recvSelf");
			return res;
		}
	}

	// Determine how many classTags there are and allocate ID vector
   	int numTags = (theSectionOrder == 0) ? numMats : numMats + 1;
	ID classTags(numTags*2 + order);

	// Receive the material class and db tags
    res += theChannel.recvID(otherDbTag, cTag, classTags);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive classTags ID",
			"SectionAggregator::recvSelf");
		return res;
	}

	// Check if null pointer, allocate if so
	if (theAdditions == 0) {
		theAdditions = new UniaxialMaterial *[numMats];
		if (theAdditions == 0) {
			g3ErrorHandler->warning("%s -- could not allocate UniaxialMaterial array",
					"SectionAggregator::recvSelf");
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
			g3ErrorHandler->warning("%s -- could not get UniaxialMaterial, i = %d",
				"SectionAggregator::recvSelf", i);
			return -1;
		}

		// Now, receive the UniaxialMaterial
		theAdditions[i]->setDbTag(classTags(i+numTags));
		res += theAdditions[i]->recvSelf(cTag, theChannel, theBroker);
		if (res < 0) {
			g3ErrorHandler->warning("%s -- could not receive UniaxialMaterial, i = %d",
				"SectionAggregator::recvSelf", i);
			return res;
		}
	}

	// If there is no Section to receive, return
	if (theSectionOrder == 0)
		return res;

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
		g3ErrorHandler->warning("%s -- could not get a SectionForceDeformation",
			"SectionAggregator::recvSelf");
		return -1;
	}

	// Now, receive the Section
	theSection->setDbTag(classTags(2*numTags-1));
	res += theSection->recvSelf(cTag, theChannel, theBroker);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive SectionForceDeformation",
			"SectionAggregator::recvSelf");
		return res;
	}

	// Check if section ID code is null or wrong size, reallocate if so
	if (code == 0)
		code = new ID(order);
	else if (code->Size() != order) {
		delete code;
		code = new ID(order);
	}
	if (code == 0) {
		g3ErrorHandler->warning("%s -- could not allocate new code ID",
			"SectionAggregator::recvSelf");
		return -1;
	}

	// Fill in the section code
	int j = 2*numTags;
	for (i = 0; i < order; i++, j++)
		(*code)(i) = classTags(j);

    return res;
}

void
SectionAggregator::Print(ostream &s, int flag)
{
    s << "\nSection Aggregator, tag: " << this->getTag() << endl;
	s << "\tSection code: " << *code;
    if (theSection)
		s << "\tSection, tag: " << theSection->getTag() << endl;
    s << "\tUniaxial Additions" << endl;
    for (int i = 0; i < numMats; i++)
		s << "\t\tUniaxial Material, tag: " << theAdditions[i]->getTag() << endl;
}

int
SectionAggregator::setVariable(const char *argv)
{
	// Axial strain
	if (strcmp(argv,"axialStrain") == 0)
		return 1;
	// Curvature about the section z-axis
	else if (strcmp(argv,"curvatureZ") == 0)
		return 2;
	// Curvature about the section y-axis
	else if (strcmp(argv,"curvatureY") == 0)
		return 3;
	else
		return -1;
}

int
SectionAggregator::getVariable(int variableID, double &info)
{
	int i;

	info = 0.0;

	switch (variableID) {
		case 1:	// Axial strain
			// Series model ... add all sources of deformation
			for (i = 0; i < order; i++)
				if ((*code)(i) == SECTION_RESPONSE_P)
					info += e(i);
			return 0;
		case 2:	// Curvature about the section z-axis
			for (i = 0; i < order; i++)
				if ((*code)(i) == SECTION_RESPONSE_MZ)
					info += e(i);
			return 0;
		case 3:	// Curvature about the section y-axis
			for (i = 0; i < order; i++)
				if ((*code)(i) == SECTION_RESPONSE_MY)
					info += e(i);
			return 0;
		default:
			return -1;
	}
}
