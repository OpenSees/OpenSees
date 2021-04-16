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

// $Revision: 1.21 $
// $Date: 2010-09-16 00:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/partitioned/PartitionedDomain.cpp,v $

// Written: fmk
// Revision: A
//
// Description: This file contains the class definition for PartitionedDomain.
// PartitionedDomain is an abstract class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints
// and MP_Constraints just like a normal domain. In addition the domain provides
// a method to partition the domain into Subdomains.
//
// ModelBuilder. There are no partitions in a PartitionedDomain.
//
// What: "@(#) PartitionedDomain.C, revA"

#include <PartitionedDomain.h>
#include <stdlib.h>

#include <Matrix.h>

#include <DomainPartitioner.h>
#include <Element.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <MP_Constraint.h>
#include <ArrayOfTaggedObjects.h>
#include <ArrayOfTaggedObjectsIter.h>
#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>
#include <Subdomain.h>
#include <PartitionedDomain.h>
#include <PartitionedDomainEleIter.h>
#include <PartitionedDomainSubIter.h>
#include <SingleDomEleIter.h>
#include <Vertex.h>
#include <Graph.h>
#include <LoadPattern.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <SP_Constraint.h>
#include <Recorder.h>
#include <Parameter.h>
#include <ParameterIter.h>

#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

#include <SubdomainIter.h>

#include <FileStream.h>

typedef map<int, int>         MAP_INT;
typedef MAP_INT::value_type   MAP_INT_TYPE;
typedef MAP_INT::iterator     MAP_INT_ITERATOR;

typedef map<int, ID *> MAP_ID;
typedef MAP_ID::value_type   MAP_ID_TYPE;
typedef MAP_ID::iterator     MAP_ID_ITERATOR;

typedef map<int, Vertex *> MAP_VERTEX;
typedef MAP_VERTEX::value_type   MAP_VERTEX_TYPE;
typedef MAP_VERTEX::iterator     MAP_VERTEX_ITERATOR;

PartitionedDomain::PartitionedDomain()
  : Domain(),
    theSubdomains(0), theDomainPartitioner(0),
    theSubdomainIter(0), mySubdomainGraph(0), has_sent_yet(false)
{
  elements = new MapOfTaggedObjects();//(1024);
  theSubdomains = new ArrayOfTaggedObjects(32);
  theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

  mainEleIter = new SingleDomEleIter(elements);
  theEleIter = new PartitionedDomainEleIter(this);

  if (theSubdomains == 0 || elements == 0 ||
      theSubdomainIter == 0 ||
      theEleIter == 0 || mainEleIter == 0) {

    opserr << "FATAL: PartitionedDomain::PartitionedDomain ";
    opserr << "  - ran out of memory\n";
    exit(-1);
  }
}


PartitionedDomain::PartitionedDomain(DomainPartitioner &thePartitioner)
  : Domain(),
    theSubdomains(0), theDomainPartitioner(&thePartitioner),
    theSubdomainIter(0), mySubdomainGraph(0), has_sent_yet(false)
{
  elements = new MapOfTaggedObjects();//(1024);
  theSubdomains = new ArrayOfTaggedObjects(32);
  theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

  mainEleIter = new SingleDomEleIter(elements);
  theEleIter = new PartitionedDomainEleIter(this);

  if (theSubdomains == 0 || elements == 0 ||
      theSubdomainIter == 0 || theDomainPartitioner == 0 ||
      theEleIter == 0 || mainEleIter == 0) {

    opserr << "FATAL: PartitionedDomain::PartitionedDomain ";
    opserr << "  - ran out of memory\n";
    exit(-1);
  }
}


PartitionedDomain::PartitionedDomain(int numNodes, int numElements,
                                     int numSPs, int numMPs, int numLoadPatterns,
                                     int numSubdomains,
                                     DomainPartitioner &thePartitioner)

  : Domain(numNodes, 0, numSPs, numMPs, numLoadPatterns),
    theSubdomains(0), theDomainPartitioner(&thePartitioner),
    theSubdomainIter(0), mySubdomainGraph(0), has_sent_yet(false)
{
  elements = new MapOfTaggedObjects();//(numElements);
  theSubdomains = new ArrayOfTaggedObjects(numSubdomains);
  theSubdomainIter = new PartitionedDomainSubIter(theSubdomains);

  mainEleIter = new SingleDomEleIter(elements);
  theEleIter = new PartitionedDomainEleIter(this);

  if (theSubdomains == 0 || elements == 0 ||
      theSubdomainIter == 0 ||
      theEleIter == 0 || mainEleIter == 0) {

    opserr << "FATAL: PartitionedDomain::PartitionedDomain(int ..) ";
    opserr << "  - ran out of memory\n";
    exit(-1);
  }
}




PartitionedDomain::~PartitionedDomain()
{
  this->PartitionedDomain::clearAll();

  if (elements != 0)
    delete elements;

  if (theSubdomains != 0)
    delete theSubdomains;

  if (theSubdomainIter != 0)
    delete theSubdomainIter;

  if (theEleIter != 0)
    delete theEleIter;
}

void
PartitionedDomain::clearAll(void)
{
  this->removeRecorders();

  SubdomainIter &mySubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = mySubdomains()) != 0)
    theSub->clearAll();

  theSubdomains->clearAll();
  this->Domain::clearAll();
  elements->clearAll();
}



bool
PartitionedDomain::addElement(Element *elePtr)
{
  if (elePtr->isSubdomain() == true)
    return this->addSubdomain((Subdomain *)elePtr);

  int eleTag = elePtr->getTag();
#ifdef _DEBUG

  // check ele Tag >= 0
  if (eleTag < 0) {
    opserr << "PartitionedDomain::addElement - Element " << eleTag;
    opserr << " tag must be >= 0\n";
    return false;
  }

  // check its not in this or any of the subdomains
  // MISSING CODE

  // check all the elements nodes exist in the domain
  const ID &nodes = elePtr->getExternalNodes();
  for (int i = 0; i < nodes.Size(); i++) {
    int nodeTag = nodes(i);
    Node *nodePtr = this->getNode(nodeTag);
    if (nodePtr == 0) {
      opserr << "PartitionedDomain::addElement - In element " << eleTag;
      opserr << " no node " << nodeTag << " exists in the domain\n";
      return false;
    }
  }

#endif

  TaggedObject *other = elements->getComponentPtr(eleTag);
  if (other != 0)
    return false;

  bool result = elements->addComponent(elePtr);
  if (result == true) {
    elePtr->setDomain(this);
    elePtr->update();
    this->domainChange();
  }

  return result;
}




bool
PartitionedDomain::addNode(Node *nodePtr)
{
#ifdef _DEBUG

#endif
  return (this->Domain::addNode(nodePtr));
}



bool
PartitionedDomain::addSP_Constraint(SP_Constraint *load)
{
  int nodeTag = load->getNodeTag();
  bool ok = false;

  if (!has_sent_yet)
  {
      Node *nodePtr = this->getNode(nodeTag);
      if (nodePtr != 0) {
        return this->Domain::addSP_Constraint(load);
      } else 
        return false;
  }

  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    ok = this->Domain::addSP_Constraint(load);
    if (ok == false) {
      return ok;
    }
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true) {
      ok = theSub->addSP_Constraint(load);
      if (ok == false) {
        opserr << "PartitiondDomain::addSP - failed to add to remote subdomain\n";
        load->Print(opserr);
        return ok;
      }
    }
  }

  // if no subdomain .. node not in model .. error message and return failure
  if (ok == false) {
    opserr << "PartitionedDomain::addSP_Constraint - cannot add as node with tag" <<
           nodeTag << "does not exist in model\n";
  }

  return ok;
}



int
PartitionedDomain::addSP_Constraint(int axisDirn, double axisValue,
                                    const ID &fixityCodes, double tol)
{
  int numAdded = 0;

  numAdded = this->Domain::addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);

  if (!has_sent_yet)
  {
    return numAdded;
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    numAdded += theSub->addSP_Constraint(axisDirn, axisValue, fixityCodes, tol);
  }

  return numAdded;
}


bool
PartitionedDomain::addSP_Constraint(SP_Constraint *load, int pattern)
{
  int nodeTag = load->getNodeTag();
  bool ok = false;

  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    ok = this->Domain::addSP_Constraint(load, pattern);
    if (ok == false)
      return false;
  }

  if (!has_sent_yet)
  {
    return ok;
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true)  {
      ok = theSub->addSP_Constraint(load, pattern);
      if (ok == false)
        return false;
    }
  }


  // if no subdomain .. node not in model .. error message and return failure
  if (ok == false) {
    opserr << "PartitionedDomain::addSP_Constraint - cannot add as node with tag" <<
           nodeTag << "does not exist in model\n";
  }

  return ok;
}


bool
PartitionedDomain::addMP_Constraint(MP_Constraint *load) {
  bool res = false;
  bool getRetained = false;
  bool addedMain = false;

  // to every domain with the constrained node we must
  // 1. add the retained node if not already there & any sp constraints on that node
  // 2. add the constraint.

  int retainedNodeTag = load->getNodeRetained();
  int constrainedNodeTag = load->getNodeConstrained();

  //
  // first we check the main domain
  // if has the constrained but not retained we mark as needing retained
  //

  Node *retainedNodePtr = this->Domain::getNode(retainedNodeTag);
  Node *constrainedNodePtr = this->Domain::getNode(constrainedNodeTag);
  if (constrainedNodePtr != 0) {
    if (retainedNodePtr != 0) {

      res = this->Domain::addMP_Constraint(load);
      if (res == false) {
        opserr << "PartitionedDomain::addMP_Constraint - problems adding to main domain\n";
        return res;
      } else {
        addedMain = true;
      }

    } else {
      getRetained = true;
    }
  }

  if (!has_sent_yet)
  {
    return addedMain;
  }

  //
  // now we check all subdomains
  // if a subdomain has both nodes we add the constraint, if only the
  // constrained node we mark as needing retained
  //

  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {

    bool hasConstrained = theSub->hasNode(constrainedNodeTag);

    if (hasConstrained == true) {

      bool hasRestrained = theSub->hasNode(retainedNodeTag);

      if (hasRestrained == true) {

        res = theSub->addMP_Constraint(load);

        if (res == false) {
          opserr << "PartitionedDomain::addMP_Constraint - problems adding to subdomain with retained\n";
          return res;
        } else {
          ;
        }
      } else {
        getRetained = true;
      }
    }
  }

  //
  // if getRetained is true, a subdomain or main domain has the constrained node
  // but no retained node .. we must go get it && SP constraints as well
  // 1. we first go get it
  // 2. then we add to main domain
  // 3. then we add to any subdomain
  //

  if (getRetained == true) {

    // we get a copy of the retained Node, set mass equal 0 (don't want to double it)
    if (retainedNodePtr == 0) {
      SubdomainIter &theSubdomains2 = this->getSubdomains();
      while ((theSub = theSubdomains2()) != 0 && retainedNodePtr == 0) {

        bool hasRestrained = theSub->hasNode(retainedNodeTag);

        if (hasRestrained == true) {
          retainedNodePtr = theSub->getNode(retainedNodeTag);

          Matrix mass(retainedNodePtr->getNumberDOF(), retainedNodePtr->getNumberDOF());
          mass.Zero();
          retainedNodePtr->setMass(mass);
        }
      }
    } else {
      // get a copy & zero the mass
      retainedNodePtr = new Node(*retainedNodePtr, false);
    }

    if (retainedNodePtr == 0) {
      opserr << "partitionedDomain::addMP_Constraint - can't find retained node anywhere!\n";
      return false;
    } else {
      ;
    }

    //
    // if main has it we add the retained to main & constraint
    //

    if (constrainedNodePtr != 0 && addedMain == false) {
      res = this->Domain::addNode(retainedNodePtr);

      if (res == false) {
        opserr << "PartitionedDomain::addMP_Constraint - problems adding retained to main domain\n";
        return res;
      }
      res = this->Domain::addMP_Constraint(load);

      if (res == false) {
        opserr << "PartitionedDomain::addMP_Constraint - problems adding constraint to main domain after adding node\n";
        return res;
      }
    }

    //
    // to subdmains that have the constrained but no retained
    // 1. we add a copy of retained
    // 2. we add the constraint
    //

    SubdomainIter &theSubdomains3 = this->getSubdomains();
    while ((theSub = theSubdomains3()) != 0) {
      bool hasConstrained = theSub->hasNode(constrainedNodeTag);
      if (hasConstrained == true) {

        bool hasRestrained = theSub->hasNode(retainedNodeTag);
        if (hasRestrained == false) {

          res = theSub->addExternalNode(retainedNodePtr);

          if (res == false) {
            opserr << "PartitionedDomain::addMP_Constraint - problems adding retained to subdomain\n";
            return res;
          }

          res = theSub->addMP_Constraint(load);

          if (res == false) {
            opserr << "PartitionedDomain::addMP_Constraint - problems adding constraint to subdomain after adding node\n";
            return res;
          }
        }
      }
    }

    // clean up memory
    if (addedMain == true && getRetained == true)
      delete retainedNodePtr;

  }

  return res;

}


bool
PartitionedDomain::addLoadPattern(LoadPattern *loadPattern)
{
  bool result = true;

  int tag = loadPattern->getTag();
  if (this->getLoadPattern(tag) != 0) {
    opserr << "PartitionedDomain::addLoadPattern - cannot add as LoadPattern with tag" <<
           tag << "already exists in model\n";
    return false;
  }

  if (has_sent_yet)
  {
    SubdomainIter &theSubdomains = this->getSubdomains();
    Subdomain *theSub;
    while ((theSub = theSubdomains()) != 0) {
      bool res = theSub->addLoadPattern(loadPattern);
      if (res != true) {
        opserr << "PartitionedDomain::addLoadPattern - cannot add as LoadPattern with tag: " <<
               tag << " to subdomain\n";
        result = res;
      }
    }
  }

  this->Domain::addLoadPattern(loadPattern);

  return result;
}


bool
PartitionedDomain::addNodalLoad(NodalLoad *load, int pattern)
{
  int nodeTag = load->getNodeTag();

  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Node *nodePtr = this->getNode(nodeTag);
  if (nodePtr != 0) {
    return (this->Domain::addNodalLoad(load, pattern));
  }

  if (!has_sent_yet)
  {
    return false;
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasNode(nodeTag);
    if (res == true) {
      // opserr << "PartitionedDomain::addLoadPattern(LoadPattern *loadPattern) SUB " << theSub->getTag() << *load;
      return theSub->addNodalLoad(load, pattern);
    }
  }

  // if no subdomain .. node not in model
  opserr << "PartitionedDomain::addNodalLoad - cannot add as node with tag" <<
         nodeTag << "does not exist in model\n";
  return false;
}


bool
PartitionedDomain::addElementalLoad(ElementalLoad *load, int pattern)
{
  int eleTag = load->getElementTag();

  // check the Node exists in the Domain or one of Subdomains

  // if in Domain add it as external .. ignore Subdomains
  Element *elePtr = this->getElement(eleTag);
  if (elePtr != 0) {
    return (this->Domain::addElementalLoad(load, pattern));
  }

  if (!has_sent_yet)
  {
    return false;
  }

  // find subdomain with node and add it .. break if find as internal node
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *theSub;
  while ((theSub = theSubdomains()) != 0) {
    bool res = theSub->hasElement(eleTag);
    if (res == true) {
      // opserr << "PartitionedDomain::addLoadPattern(LoadPattern *loadPattern) SUB " << theSub->getTag() << *load;
      return theSub->addElementalLoad(load, pattern);
    }
  }

  // if no subdomain .. node not in model
  opserr << "PartitionedDomain::addElementalLoad - cannot add as element with tag" <<
         eleTag << "does not exist in model\n";

  return false;
}


Element *
PartitionedDomain::removeElement(int tag)
{

  // opserr << " PartitionedDomain::removeElement() - Removing element with tag " << tag <<  endln;
  // we first see if its in the original domain
  TaggedObject *res = elements->removeComponent(tag);
  Element *result = 0;
  if (res != 0) {
    result = (Element *)res;
    this->domainChange();
    return result;
  }
  if (!has_sent_yet)
  {
    return 0;
  }
  // opserr << " PartitionedDomain::removeElement() - Searching subdomains for element #  " << tag <<  endln;

  // if not there we must check all the other subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result = theSub->removeElement(tag);
      if (result != 0) {
        return result;
      }
    }
  }
  // opserr << " PartitionedDomain::removeElement() - element not found, tag =  " << tag <<  endln;

  // its not there
  return 0;
}


Node *
PartitionedDomain::removeNode(int tag)
{
  // we first remove it form the original domain (in case on boundary)
  Node *result = this->Domain::removeNode(tag);
  if (!has_sent_yet)
  {
    return result;
  }
  // we must also try removing from the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      Node *res = theSub->removeNode(tag);
      if (res != 0)
        result = res;
    }
  }

  if (result != 0)
    this->domainChange();

  return result;
}

SP_Constraint *
PartitionedDomain::removeSP_Constraint(int tag)
{
  // we first see if its in the original domain
  SP_Constraint *result = this->Domain::removeSP_Constraint(tag);
  if (result != 0) {
    this->domainChange();
    return result;
  }
  if (!has_sent_yet)
  {
    return result;
  }

  // if not there we must check all the other subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result = theSub->removeSP_Constraint(tag);
      if (result != 0) {
        return result;
      }
    }
  }

  // its not there
  return 0;
}


int
PartitionedDomain::removeSP_Constraint(int nodeTag, int dof, int loadPatternTag)
{

  // we first see if its in the original domain
  int result = this->Domain::removeSP_Constraint(nodeTag, dof, loadPatternTag);
  if (!has_sent_yet)
  {
    return result;
  }

  // if not there we must check all the other subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result += theSub->removeSP_Constraint(nodeTag, dof, loadPatternTag);
    }
  }

  if (result != 0) {
    this->domainChange();
  }

  // its not there
  return result;
}



MP_Constraint *
PartitionedDomain::removeMP_Constraint(int tag)
{
  // we first see if its in the original domain
  MP_Constraint *result = this->Domain::removeMP_Constraint(tag);
  if (result != 0) {
    this->domainChange();
    return result;
  }
  if (!has_sent_yet)
  {
    return result;
  }

  // if not there we must check all the other subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result = theSub->removeMP_Constraint(tag);
      if (result != 0) {
        return result;
      }
    }
  }

  // its not there
  return 0;
}


int
PartitionedDomain::removeMP_Constraints(int tag)
{
  // we first see if its in the original domain
  int result = this->Domain::removeMP_Constraints(tag);
  if (!has_sent_yet)
  {
    return result;
  }
  // if not there we must check all the other subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result += theSub->removeMP_Constraints(tag);
    }
  }

  if (result != 0) {
    this->domainChange();
  }

  // its not there
  return result;
}


LoadPattern *
PartitionedDomain::removeLoadPattern(int tag)
{
  // we first see if its in the original domain
  LoadPattern *result = this->Domain::removeLoadPattern(tag);
  if (!has_sent_yet)
  {
    return result;
  }

  // we must also try removing from the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      LoadPattern *res = theSub->removeLoadPattern(tag);
      if (res != 0)
        result = res;
    }
  }

  if (result != 0)
    this->domainChange();

  return result;
}

// public member functions which have to be modified
ElementIter       &
PartitionedDomain::getElements()
{
  theEleIter->reset();
  return *theEleIter;
}


Element  *
PartitionedDomain::getElement(int tag)
{
  // we first see if its in the original domain
  TaggedObject *res = elements->getComponentPtr(tag);
  Element *result = 0;
  if (res != 0) {
    result = (Element *)res;
    return result;
  }

  /*
  // go through the other subdomains until we find it or we run out of subdomains
  if (theSubdomains != 0) {
  ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
  TaggedObject *theObject;
  while ((theObject = theSubsIter()) != 0) {
    Subdomain *theSub = (Subdomain *)theObject;
    result = theSub->getElement(tag);
    if (result != 0)
  return result;
  }
  }
  */

  // its not there
  return 0;
}


int
PartitionedDomain::getNumElements(void) const
{
  int result = elements->getNumComponents();

  // add the number of subdomains
  result +=  theSubdomains->getNumComponents();
  return result;
}

void
PartitionedDomain::applyLoad(double timeStep)
{
  this->Domain::applyLoad(timeStep);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->applyLoad(timeStep);
    }
  }
}


void
PartitionedDomain::setCommitTag(int newTag)
{
  this->Domain::setCommitTag(newTag);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->setCommitTag(newTag);
    }
  }
}



void
PartitionedDomain::setCurrentTime(double newTime)
{
  this->Domain::setCurrentTime(newTime);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->setCurrentTime(newTime);
    }
  }
}


void
PartitionedDomain::setCommittedTime(double newTime)
{
  this->Domain::setCommittedTime(newTime);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->setCommittedTime(newTime);
    }
  }
}


void
PartitionedDomain::setLoadConstant(void)
{
  this->Domain::setLoadConstant();

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->setLoadConstant();
    }
  }
}


int
PartitionedDomain::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
  this->Domain::setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    }
  }
  return 0;
}


int
PartitionedDomain::update(void)
{
  int res = this->Domain::update();

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->computeNodalResponse();
      res += theSub->update();
    }
  }

#ifdef _PARALLEL_PROCESSING
// opserr << "PartitionedDomain:: barrierCheck\n";
  return this->barrierCheck(res);
#endif
  return res;
}


#ifdef _PARALLEL_PROCESSING
int
PartitionedDomain::barrierCheck(int res)
{
  int result = res;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int subResult = theSub->barrierCheckIN();
      if (subResult != 0)
        result = subResult;
    }

    ArrayOfTaggedObjectsIter theSubsIter1(*theSubdomains);
    while ((theObject = theSubsIter1()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->barrierCheckOUT(result);
    }
  }

  return result;
}
#endif

int
PartitionedDomain::update(double newTime, double dT)
{
  this->applyLoad(newTime);
  int res = this->Domain::update();

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->computeNodalResponse();
      res += theSub->update(newTime, dT);
    }
  }

#ifdef _PARALLEL_PROCESSING
  return this->barrierCheck(res);
#endif
  return res;

  /*

  opserr << "PartitionedDomain::update(double newTime, double dT) -1\n";
  int result = 0;


  opserr << "PartitionedDomain::update(double newTime, double dT) -2\n";
  this->update();
  opserr << "PartitionedDomain::update(double newTime, double dT) -2a\n";

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->update(newTime, dT);
    }
    this->barrierCheck(result);
  }
  opserr << "PartitionedDomain::update(double newTime, double dT) -3\n";
  return result;

  */

}


int
PartitionedDomain::hasDomainChanged(void)
{
  return this->Domain::hasDomainChanged();
}

int
PartitionedDomain::analysisStep(double dT)
{
  // first we need to see if any subdomain has changed & mark the change in domain
  bool domainChangedAnySubdomain = this->getDomainChangeFlag();
  if (domainChangedAnySubdomain == false) {
    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
        Subdomain *theSub = (Subdomain *)theObject;
        domainChangedAnySubdomain = theSub->getDomainChangeFlag();
      }
    }
  }

  if (domainChangedAnySubdomain == true) {
    this->Domain::domainChange();
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
        Subdomain *theSub = (Subdomain *)theObject;
        theSub->domainChange();
      }
    }
  }

  this->Domain::analysisStep(dT);

  int res = 0;
  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      res += theSub->analysisStep(dT);
      if (res != 0)
        opserr << "PartitionedDomain::step - subdomain " << theSub->getTag() << " failed in step\n";
    }
  }
  return res;
}



int
PartitionedDomain::eigenAnalysis(int numModes, bool generalized, bool findSmallest)
{
  // first we need to see if any subdomain has changed & mark the change in domain
  bool domainChangedAnySubdomain = this->getDomainChangeFlag();
  if (domainChangedAnySubdomain == false) {
    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
        Subdomain *theSub = (Subdomain *)theObject;
        domainChangedAnySubdomain = theSub->getDomainChangeFlag();
      }
    }
  }

  if (domainChangedAnySubdomain == true) {
    this->Domain::domainChange();
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while (((theObject = theSubsIter()) != 0) && (domainChangedAnySubdomain == false)) {
        Subdomain *theSub = (Subdomain *)theObject;
        theSub->domainChange();
      }
    }
  }

  this->Domain::eigenAnalysis(numModes, generalized, findSmallest);

  int res = 0;
  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      res += theSub->eigenAnalysis(numModes, generalized, findSmallest);
      if (res != 0)
        opserr << "PartitionedDomain::step - subdomain " << theSub->getTag() << " failed in step\n";
    }
  }
  return res;
}


int
PartitionedDomain::record(bool fromAnalysis)
{
  int result = 0;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      result += theSub->record(fromAnalysis);
      if (result < 0) {
        opserr << "PartitionedDomain::record(void)";
        opserr << " - failed in Subdomain::record()\n";
      }
    }
  }

  result += this->Domain::record();

  return result;
}

int
PartitionedDomain::commit(void)
{
  int result = this->Domain::commit();

  // static int ctag=0;
  // char buffer[50];
  // sprintf(buffer, "domaininfo_%d.txt", ctag++);
  // FileStream fid(buffer);

  // fid << "MAIN DOMAIN --------------------------------------------------\n\n";

  // this->Print(fid);

  if (result < 0) {
    opserr << "PartitionedDomain::commit(void) - failed in Domain::commit()\n";
    return result;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->commit();
      // fid << "Sub-Domain # " << theSub->getTag() << " --------------------------------------------------\n\n";
      // theSub->Print(fid);

      if (res < 0) {
        opserr << "PartitionedDomain::commit(void)";
        opserr << " - failed in Subdomain::commit()\n";
        return res;
      }
    }
  }

  // opserr << "Subdomain # MASTER " << " update_time = " << this->Domain::update_time_committed << endln;


  // now we load balance if we have subdomains and a partitioner
  int numSubdomains = this->getNumSubdomains();
  if (numSubdomains != 0 && theDomainPartitioner != 0)  {
    // opserr << "Subdomain # MASTER " << " BALANCING! " << endln;
    Graph &theSubGraphs = this->getSubdomainGraph();
    theDomainPartitioner->balance(theSubGraphs);
  }

  return 0;
}


int
PartitionedDomain::revertToLastCommit(void)
{
  int result = this->Domain::revertToLastCommit();
  if (result < 0) {
    opserr << "PartitionedDomain::revertToLastCommit(void) - failed in Domain::revertToLastCommit()\n";
    return result;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->revertToLastCommit();
      if (res < 0) {
        opserr << "PartitionedDomain::revertToLastCommit(void)";
        opserr << " - failed in Subdomain::revertToLastCommit()\n";
        return res;
      }
    }
  }

  return 0;
}

int
PartitionedDomain::revertToStart(void)
{
  int result = this->Domain::revertToStart();
  if (result < 0) {
    opserr << "PartitionedDomain::revertToLastCommit(void) - failed in Domain::revertToLastCommit()\n";
    return result;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->revertToStart();
      if (res < 0) {
        opserr << "PartitionedDomain::revertToLastCommit(void)";
        opserr << " - failed in Subdomain::revertToLastCommit()\n";
        return res;
      }
    }
  }

  return 0;
}


int
PartitionedDomain::addRecorder(Recorder &theRecorder)
{
  int result = this->Domain::addRecorder(theRecorder);
  if (result < 0)
    return -1;

  if (!has_sent_yet)
  {
    return result;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->addRecorder(theRecorder);
      if (res < 0) {
        opserr << "PartitionedDomain::revertToLastCommit(void)";
        opserr << " - failed in Subdomain::revertToLastCommit()\n";
        return res;
      }
    }
  }
  return 0;
}

int
PartitionedDomain::removeRecorders(void)
{
  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->removeRecorders();
      if (res < 0) {
        opserr << "PartitionedDomain::removeRecorders(void)";
        opserr << " - failed in Subdomain::removeRecorders()\n";
        return res;
      }
    }
  }

  if (this->Domain::removeRecorders() < 0)
    return -1;

  this->barrierCheck(1.0);

  return 0;
}


int
PartitionedDomain::removeRecorder(int tag)
{
  if (this->Domain::removeRecorder(tag) < 0)
    return -1;

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      int res = theSub->removeRecorder(tag);
      if (res < 0) {
        opserr << "PartitionedDomain::removeRecorder(void)";
        opserr << " - failed in Subdomain::removeRecorder()\n";
        return res;
      }
    }
  }
  return 0;
}

void
PartitionedDomain::Print(OPS_Stream &s, int flag)
{
  this->Domain::Print(s, flag);

  s << "\nELEMENT DATA: NumEle: " << elements->getNumComponents() << "\n";
  elements->Print(s);

  // print all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      theObject->Print(s, flag);
    }
  }
}


void
PartitionedDomain::Print(OPS_Stream &s, ID *nodeTags, ID *eleTags, int flag)
{
  if (nodeTags != 0)
    this->Domain::Print(s, nodeTags, 0, flag);

  if (eleTags != 0) {
    int numEle = eleTags->Size();
    for (int i = 0; i < numEle; i++) {
      int eleTag = (*eleTags)(i);
      TaggedObject *theEle = elements->getComponentPtr(eleTag);
      if (theEle != 0)
        theEle->Print(s, flag);
    }
  }


  // print all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->Print(s, nodeTags, eleTags, flag);
    }
  }
}



int
PartitionedDomain::setPartitioner(DomainPartitioner *thePartitioner)
{
  theDomainPartitioner = thePartitioner;

  return 0;
}


int
PartitionedDomain::partition(int numPartitions, bool usingMain, int mainPartitionID, int specialElementTag)
{
  int result = 0;
  // need to create element graph before create new subdomains
  // DO NOT REMOVE THIS LINE __ EVEN IF COMPILER WARNING ABOUT UNUSED VARIABLE
  Graph &theEleGraph = this->getElementGraph();

  // now we call partition on the domainPartitioner which does the partitioning
  DomainPartitioner *thePartitioner = this->getPartitioner();
  if (thePartitioner != 0) {
    thePartitioner->setPartitionedDomain(*this);
    result =  thePartitioner->partition(numPartitions, usingMain, mainPartitionID, specialElementTag);
  } else {
    opserr << "PartitionedDomain::partition(int numPartitions) - no associated partitioner\n";
    return -1;
  }

  //
  // add recorder objects
  //

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      for (int i = 0; i < numRecorders; i++) {
        int res = theSub->addRecorder(*theRecorders[i]);
        if (res != 0) {
          opserr << "PartitionedDomain::partition(void)";
          opserr << " - failed to add Recorder to subdomain\n";
          return res;
        }
      }
    }
  }

  //
  // add parameters
  //
  ParameterIter     &theParameters = this->getParameters();
  Parameter *theParameter;
  while ((theParameter = theParameters()) != 0) {
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0) {
        Subdomain *theSub = (Subdomain *)theObject;
        int res = theSub->addParameter(theParameter);
        if (res != 0) {
          opserr << "PartitionedDomain::partition(void)";
          opserr << " - failed to add Parameter to subdomain\n";
          return res;
        }
      }
    }
    theParameter->setDomain(this);  // needed as some components will move
  }

  has_sent_yet = true;

  return result;
}


int PartitionedDomain::repartition(int numPartitions, bool usingMain, int mainPartitionID, int specialElementTag)
{
  opserr << "PartitionedDomain::repartition() - On Master Process\n";



  return 0;
}


bool
PartitionedDomain::addSubdomain(Subdomain *theSubdomain)
{
  int eleTag = theSubdomain->getTag();
  TaggedObject *other = theSubdomains->getComponentPtr(eleTag);
  if (other != 0)
    return false;

  bool result = theSubdomains->addComponent(theSubdomain);
  if (result == true) {
    theSubdomain->setDomain(this);
    this->domainChange();
  }

  return result;

}

int
PartitionedDomain::getNumSubdomains(void)
{
  return theSubdomains->getNumComponents();
}

Subdomain *
PartitionedDomain::getSubdomainPtr(int tag)
{
  TaggedObject *mc = theSubdomains->getComponentPtr(tag);
  if (mc == 0) return 0;
  Subdomain *result = (Subdomain *)mc;
  return result;
}

SubdomainIter &
PartitionedDomain::getSubdomains(void)
{
  theSubdomainIter->reset();
  return *theSubdomainIter;
}



DomainPartitioner *
PartitionedDomain::getPartitioner(void) const
{
  return theDomainPartitioner;
}



int
PartitionedDomain::buildEleGraph(Graph *theEleGraph)
{
  // see if quick return
  int numVertex = this->getNumElements();
  if (numVertex == 0)
    return 0;

  //
  // iterate over the lements of the domain
  //  create a vertex with a unique tag for each element
  //  also create a map to hold element tag - vertex tag mapping
  //

  MAP_INT theEleToVertexMap;
  MAP_INT_ITERATOR theEleToVertexMapEle;


  TaggedObject *theTagged;
  TaggedObjectIter &theElements = elements->getComponents();
  int count = START_VERTEX_NUM;
  while ((theTagged = theElements()) != 0) {
    int eleTag = theTagged->getTag();
    Vertex *vertexPtr = new Vertex(count, eleTag);

    if (vertexPtr == 0) {
      opserr << "WARNING Domain::buildEleGraph - Not Enough Memory to create the " << count << " vertex\n";
      return -1;
    }

    // Get the compute cost and communications cost.
    Element * theElement =  static_cast<Element *>(theTagged);
    double eleWeight = (double) theElement->getNumDOF(); //theElement->getTime();
    int eleCommCost = 0;//theElement->getMoveCost();
    vertexPtr->setWeight(eleWeight);
    // vertexPtr->setTmp(eleCommCost);

    theEleGraph->addVertex(vertexPtr);
    theEleToVertexMapEle = theEleToVertexMap.find(eleTag);
    if (theEleToVertexMapEle == theEleToVertexMap.end()) {
      theEleToVertexMap.insert(MAP_INT_TYPE(eleTag, count));

      // check if sucessfully added
      theEleToVertexMapEle = theEleToVertexMap.find(eleTag);
      if (theEleToVertexMapEle == theEleToVertexMap.end()) {
        opserr << "Domain::buildEleGraph - map STL failed to add object with tag : " << eleTag << endln;
        return false;
      }

      count++;
    }
  }

  //
  // We now need to determine which elements are asssociated with each node.
  // As this info is not in the Node interface we must build it;
  //
  // again we will use an stl map, index will be nodeTag, object will be Vertex
  // do using vertices for each node, when we addVertex at thes nodes we
  // will not be adding vertices but element tags.
  //

  MAP_ID theNodeToVertexMap;
  MAP_ID_ITERATOR theNodeEle;

  Node *nodPtr;

  // now create the vertices with a reference equal to the node number.
  // and a tag which ranges from 0 through numVertex-1 and placed in
  // theNodeTagVertices at a position equal to the node's tag.

  NodeIter &theNodes = this->getNodes();
  while ((nodPtr = theNodes()) != 0) {
    int nodeTag = nodPtr->getTag();
    ID *eleTags = new ID(0, 4);

    if (eleTags == 0) {
      opserr << "WARNING Domain::buildEleGraph - Not Enough Memory to create the " << count << " vertex\n";
      return -1;
    }

    theNodeEle = theNodeToVertexMap.find(nodeTag);
    if (theNodeEle == theNodeToVertexMap.end()) {
      theNodeToVertexMap.insert(MAP_ID_TYPE(nodeTag, eleTags));

      // check if sucessfully added
      theNodeEle = theNodeToVertexMap.find(nodeTag);
      if (theNodeEle == theNodeToVertexMap.end()) {
        opserr << "Domain::buildEleGraph - map STL failed to add object with tag : " << nodeTag << endln;
        return false;
      }
    }
  }

  // now add the the Elements to the node vertices
  Element *theEle;
  TaggedObjectIter &theEle3 = elements->getComponents();

  while ((theTagged = theEle3()) != 0) {
    theEle = (Element *)theTagged;
    int eleTag = theEle->getTag();
    const ID &id = theEle->getExternalNodes();

    int size = id.Size();
    for (int i = 0; i < size; i++) {
      int nodeTag = id(i);

      MAP_ID_ITERATOR theNodeEle;
      theNodeEle = theNodeToVertexMap.find(nodeTag);
      if (theNodeEle == theNodeToVertexMap.end()) {
        return -1;
      } else {
        ID *theNodeEleTags = (*theNodeEle).second;
        theNodeEleTags->insert(eleTag);
      }
    }
  }

  //
  // now add the edges to the vertices of our element graph;
  // this is done by looping over the Node vertices, getting their
  // Adjacenecy and adding edges between elements with common nodes
  //

  MAP_ID_ITERATOR currentComponent;
  currentComponent = theNodeToVertexMap.begin();
  while (currentComponent != theNodeToVertexMap.end()) {
    ID *id = (*currentComponent).second;

    int size = id->Size();
    for (int i = 0; i < size; i++) {
      int eleTag1 = (*id)(i);

      theEleToVertexMapEle = theEleToVertexMap.find(eleTag1);
      if (theEleToVertexMapEle != theEleToVertexMap.end()) {
        int vertexTag1 = (*theEleToVertexMapEle).second;

        for (int j = 0; j < size; j++)
          if (i != j) {
            int eleTag2 = (*id)(j);
            theEleToVertexMapEle = theEleToVertexMap.find(eleTag2);
            if (theEleToVertexMapEle != theEleToVertexMap.end()) {
              int vertexTag2 = (*theEleToVertexMapEle).second;

              // addEdge() adds for both vertices - do only once

              if (vertexTag1 > vertexTag2) {
                theEleGraph->addEdge(vertexTag1, vertexTag2);
                theEleGraph->addEdge(vertexTag2, vertexTag1);
              }
            }
          }
      }
    }
    currentComponent++;
  }

  // clean up - delete the ID's associated with the nodes
  currentComponent = theNodeToVertexMap.begin();
  while (currentComponent != theNodeToVertexMap.end()) {
    delete (*currentComponent).second;
    currentComponent++;
  }

  return 0;
}



// a method which will only remove a node from the partitioned domain
// it does not touch the subdomains .. can be dangerous to use.
Node *
PartitionedDomain::removeExternalNode(int tag)
{
  return (this->Domain::removeNode(tag));
}

Graph &
PartitionedDomain::getSubdomainGraph(void)
{
  //The subdomain graph has number of vertices equal to number of subdomains, including P0
  int numVertex = theSubdomains->getNumComponents() + 1;   // Use P0 as a vertex toooo

  // delete the old always - only object that will
  // use this is a DomainBalancer & it is always looking for latest
  if (mySubdomainGraph != 0) {
    delete mySubdomainGraph;
    mySubdomainGraph = 0;
  }

  // create a new graph
  if (mySubdomainGraph == 0)
    mySubdomainGraph = new Graph(numVertex + START_VERTEX_NUM);

  if (mySubdomainGraph == 0) // if still 0 try a smaller one
    mySubdomainGraph = new Graph();


  // see if quick return
  if (numVertex == 0)
    return *mySubdomainGraph;

  //Mapping element tags to vertex corresponding vertex
  MAP_INT subDTagToVtxTag;

  //Add P0 to the graph
  SubdomainIter &theSubdomains = this->getSubdomains();
  Subdomain *subDPtr = 0;
  int subDTag = 1;
  double myCostP0 = 0.0;//this->getUpdateTime();

  Vertex *selfvertexPtr = new Vertex(subDTag, subDTag, myCostP0);
  mySubdomainGraph->addVertex(selfvertexPtr);
  
  subDTagToVtxTag.insert(MAP_INT_TYPE(subDTag, subDTag));

  while ((subDPtr = theSubdomains()) != 0) {    
    int subDTag = subDPtr->getTag();
    
    //Vertex for current subD
    Vertex *vertexPtr = new Vertex(subDTag, subDTag, subDPtr->getCost());

    //Add to graph and map
    mySubdomainGraph->addVertex(vertexPtr);
    subDTagToVtxTag.insert(MAP_INT_TYPE(subDTag, subDTag));
  }

  // We now need to determine which theSubdomains are asssociated with each node.
  // As this info is not in the Node interface we must build it; which we
  // do using vertices for each node, when we addVertex at thes nodes we
  // will not be adding vertices but element tags.

  MAP_VERTEX   nodeTagToVtx;


  // now create the vertices with a reference equal to the node number.
  // and a tag which ranges from 0 through numVertex-1 and placed in
  // theNodeTagVertices at a position equal to the node's tag.

  NodeIter &niter = this->getNodes();
  Node *nodPtr = 0;
  int count = START_VERTEX_NUM;
  while ((nodPtr = niter()) != 0) {
    int nodeTag = nodPtr->getTag();
    Vertex *vertexPtr = new Vertex(count++, nodeTag);
    vertexPtr->addEdge(1);

    nodeTagToVtx.insert(MAP_VERTEX_TYPE(nodeTag, vertexPtr));
  }

  // now add the the TheSubdomains to the nodes
  SubdomainIter &theSubdomains2 = this->getSubdomains();
  while ((subDPtr = theSubdomains2()) != 0) {
    int subdTag = subDPtr->getTag();
    const ID &subdEdgeNodes = subDPtr->getExternalNodes();
    int numberExternalNodes = subdEdgeNodes.Size();
    for (int i = 0; i < numberExternalNodes; i++)
    {
      int currentNodeTag = subdEdgeNodes(i);

      MAP_VERTEX_ITERATOR mpvtx = nodeTagToVtx.find(currentNodeTag);
      if(mpvtx != nodeTagToVtx.end()) // If we find current node Tag....
      {
        mpvtx->second->addEdge(subdTag);
      }
      else //If not found
      {
        Vertex *newvtxptr = new Vertex(count++, currentNodeTag);
        newvtxptr->addEdge(subdTag);
        nodeTagToVtx.insert(MAP_VERTEX_TYPE(currentNodeTag, newvtxptr));
      }
    }
  }

  // now add the edges to the vertices of our element graph;
  // this is done by looping over the Node vertices, getting their
  // Adjacenecy and adding edges between theSubdomains with common nodes

  for (auto it=nodeTagToVtx.begin(); it!=nodeTagToVtx.end(); ++it)
  {
    int nodeTag = it->first;
    Vertex *vertexPtr = it->second;
    
    const ID &connectedSubdomains = vertexPtr->getAdjacency();

    int numberOfConnectedSubdomains = connectedSubdomains.Size();
    for (int i = 0; i < numberOfConnectedSubdomains; i++) {
      int subdomainTag1 = connectedSubdomains(i);
      int vertexTag1 = subDTagToVtxTag.find(subdomainTag1)->first;

      for (int j = 0; j < numberOfConnectedSubdomains; j++)
        if (i != j) {
          int subdomainTag2 = connectedSubdomains(j);

          int vertexTag2 = subDTagToVtxTag.find(subdomainTag2)->second;

          // addEdge() adds for both vertices - do only once
          if (vertexTag1 > vertexTag2)
          {
            mySubdomainGraph->addEdge(vertexTag1, vertexTag2);
            mySubdomainGraph->addEdge(vertexTag2, vertexTag1);
          }
        }
    }
  }

  for (auto it=nodeTagToVtx.begin(); it!=nodeTagToVtx.end(); ++it)
  {
    delete it->second;
  }


  return *mySubdomainGraph;
}


double
PartitionedDomain::getNodeDisp(int nodeTag, int dof, int &errorFlag)
{
  double result = this->Domain::getNodeDisp(nodeTag, dof, errorFlag);

  if (errorFlag != 0) {

    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0 && errorFlag != 0) {
        Subdomain *theSub = (Subdomain *)theObject;
        result = theSub->getNodeDisp(nodeTag, dof, errorFlag);
        if (errorFlag == 0)
          return result;
      }
    }
  }

  return result;
}


int
PartitionedDomain::setMass(const Matrix &mass, int nodeTag)
{
  int result = this->Domain::setMass(mass, nodeTag);

  if (result != 0) {

    // do the same for all the subdomains
    if (theSubdomains != 0) {
      ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
      TaggedObject *theObject;
      while ((theObject = theSubsIter()) != 0 && result != 0) {
        Subdomain *theSub = (Subdomain *)theObject;
        result = theSub->setMass(mass, nodeTag);
      }
    }
  }

  return result;
}

const Vector *
PartitionedDomain::getNodeResponse(int nodeTag, NodeResponseType response)
{
  const Vector *res = this->Domain::getNodeResponse(nodeTag, response);
  if (res != 0)
    return res;

  // opserr << "PartitionedDomain::getNodeResponse - Did not find response for node # " << nodeTag << " in main domain\n";

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      const Vector *result = theSub->getNodeResponse(nodeTag, response);
      if (result != 0)
        return result;
      // else
        // opserr << "PartitionedDomain::getNodeResponse - Did not find response for node # " << nodeTag << " in domain # " << theSub->getTag() << "\n";
    }
  }

  return NULL;
}

const Vector *
PartitionedDomain::getElementResponse(int eleTag, const char **argv, int argc) {

  const Vector *res = this->Domain::getElementResponse(eleTag, argv, argc);
  if (res != 0) {
    return res;
  }

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      const Vector *result = theSub->getElementResponse(eleTag, argv, argc);
      if (result != 0) {
        return result;
      }
    }
  }

  return NULL;
}

int
PartitionedDomain::calculateNodalReactions(bool inclInertia)
{
  int res = this->Domain::calculateNodalReactions(inclInertia);

  // do the same for all the subdomains
  /*
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      res += theSub->calculateNodalReactions(inclInertia);
    }
  }
  */
  return res;
}

bool
PartitionedDomain::addParameter(Parameter *param)
{
  bool res = this->Domain::addParameter(param);

  // do the same for all the subdomains
  if (theSubdomains != 0) {

    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->addParameter(param);

    }
  }

  return res;
}

Parameter *
PartitionedDomain::removeParameter(int tag)
{
  Parameter *res = this->Domain::removeParameter(tag);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      theSub->removeParameter(tag);
    }
  }

  return res;
}


int
PartitionedDomain::updateParameter(int tag, int value)
{
  int res = 0;

  // do the same for all the subdomains

  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {

      Subdomain *theSub = (Subdomain *)theObject;
      res += theSub->updateParameter(tag, value);

    }
  }

  res += this->Domain::updateParameter(tag, value);

  return res;
}


int
PartitionedDomain::updateParameter(int tag, double value)
{

  int res = 0;

  res += this->Domain::updateParameter(tag, value);

  // do the same for all the subdomains
  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {
      Subdomain *theSub = (Subdomain *)theObject;
      res += theSub->updateParameter(tag, value);
    }
  }

  res += this->Domain::updateParameter(tag, value);

  return res;
}

TaggedObjectStorage* PartitionedDomain::getElementsStorage()
{
  return elements;
}


GraphPartitioner* PartitionedDomain::getGraphPartitioner(void)
{
  return theDomainPartitioner->getGraphPartitioner();
}



int
PartitionedDomain::activateElements(const ID& elementList)
{
  int res = 0;

  // do the same for all the subdomains

  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {

      Subdomain *theSub = (Subdomain *)theObject;

      // opserr << "PartitionedDomain::activateElements elementList = " << elementList << endln;
      res += theSub->activateElements(elementList);

    }
  }

  res += this->Domain::activateElements(elementList);

  return res;
}

int
PartitionedDomain::deactivateElements(const ID& elementList)
{
  int res = 0;

  // do the same for all the subdomains

  if (theSubdomains != 0) {
    ArrayOfTaggedObjectsIter theSubsIter(*theSubdomains);
    TaggedObject *theObject;
    while ((theObject = theSubsIter()) != 0) {

      Subdomain *theSub = (Subdomain *)theObject;

      // opserr << "PartitionedDomain::activateElements elementList = " << elementList << endln;
      res += theSub->deactivateElements(elementList);

    }
  }

  res += this->Domain::deactivateElements(elementList);

  return res;
}
