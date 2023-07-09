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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:01:47 $
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/PartitionedModelBuilder.cpp,v $
                                                                        
                                                                        
// File: ~/model/PartitionedModelBuilder.C
//
// Written: fmk 
// Created: 03/98
// Revision: A
//

#include <PartitionedModelBuilder.h>
#include <PartitionedDomain.h>
#include <Subdomain.h>
#include <SubdomainIter.h>

//  PartitionedModelBuilderModel(Domain &theDomain, int theClassTag);
//	constructor
PartitionedModelBuilder::PartitionedModelBuilder(PartitionedDomain &aPartitionedDomain,
						 int theClassTag)
:ModelBuilder(aPartitionedDomain), MovableObject(theClassTag),
thePartitionedDomain(&aPartitionedDomain)
{

}

PartitionedModelBuilder::PartitionedModelBuilder(Subdomain &aSubdomain,
						 int theClassTag)
:ModelBuilder(aSubdomain), MovableObject(theClassTag),
thePartitionedDomain(0)
{

}

PartitionedModelBuilder::~PartitionedModelBuilder()
{
    
}

int
PartitionedModelBuilder::buildFE_Model(void)
{
  int result;

  if (thePartitionedDomain == 0) {
    opserr << "PartitionedModelBuilder::buildFE_Model(void) -";
    opserr << "No PartitionedDomain associated with this object\n";
    return -1;
  }

  // we build the interface, i.e. nodes on boundaries and any constraints and loads
  int numSubdomains = thePartitionedDomain->getNumSubdomains();
  result = this->buildInterface(numSubdomains);
  if (result != 0) {
    opserr << "PartitionedModelBuilder::buildFE_Model(void) -";
    opserr << "buildInterface failed\n";
    return result;
  }

  // now build the subdomains, stopping if an error in building any subdomain
  SubdomainIter &theSubs = thePartitionedDomain->getSubdomains();
  Subdomain *theSubdomain;
  while ((theSubdomain = theSubs()) != 0) {
    result = theSubdomain->buildSubdomain(numSubdomains, *this);
    if (result != 0) {
	opserr << "PartitionedModelBuilder::buildFE_Model(void) -";
	opserr << "buildSubdomain failed for Subdomain " << theSubdomain->getTag();
	opserr << endln;
	return result;
    }

  }

  // if got here a PartitiondDomain has been populated
  return 0;
}

PartitionedDomain *
PartitionedModelBuilder::getPartitionedDomainPtr(void) const 
{
    return thePartitionedDomain;
}
    

