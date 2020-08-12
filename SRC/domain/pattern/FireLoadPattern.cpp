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
                                                                        
// $Revision: 1.9 $
// $Date: 2006/09/05 20:51:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/EarthquakePattern.cpp,v $
                                                                        
// Based on load pattern 
// Purpose: This file contains the class definition for EarthquakePattern.
// EarthquakePattern is an abstract class.

//Modified by Panagistis Kotsoivnos,Liming Jiang[University of Edinburgh]

#include <FireLoadPattern.h>
#include <GroundMotion.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <ElementalLoad.h>
#include <LoadPattern.h>
#include <ID.h>
#include <TimeSeries.h>
#include <NodalLoad.h>
#include <SP_Constraint.h>
#include <ArrayOfTaggedObjects.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <SingleDomSP_Iter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <GroundMotion.h>

#include <OPS_Globals.h>

static int numFireLoadPattern = 0;

FireLoadPattern::FireLoadPattern(int tag, int classTag)
  :LoadPattern(tag, classTag), theSeries1(0), theSeries2(0), theSeries3(0), theSeries4(0), 
   theSeries5(0), theSeries6(0), theSeries7(0), theSeries8(0), theSeries9(0), loadFactors(9), currentTime(0.0) 
{

  loadFactors.Zero();

  if (numFireLoadPattern == 0) {
    numFireLoadPattern++;
    opserr << "Using OpenSees Thermal Extension \n\tfrom University of Edinburgh (UoE) OpenSees developers Group, Prof. Asif Usmani\n";
  }
}

FireLoadPattern::FireLoadPattern(int tag)
  :LoadPattern(tag, PATTERN_TAG_LoadPattern), theSeries1(0), theSeries2(0), theSeries3(0), theSeries4(0), 
   theSeries5(0), theSeries6(0), theSeries7(0), theSeries8(0), theSeries9(0), loadFactors(9), currentTime(0.0)
{
  loadFactors.Zero();

  if (numFireLoadPattern == 0) {
    numFireLoadPattern++;
    opserr << "Using OpenSees Thermal Extension \n\tfrom University of Edinburgh (UoE) OpenSees developers Group, Prof. Asif Usmani\n";
  }
}


FireLoadPattern::~FireLoadPattern()
{
  if (theSeries1 != 0)
    delete theSeries1;
  if (theSeries2 != 0)
    delete theSeries2;
  if (theSeries3 != 0)
    delete theSeries3;
  if (theSeries4 != 0)
    delete theSeries4;
  if (theSeries5 != 0)
    delete theSeries5;
  if (theSeries6 != 0)
    delete theSeries6;
  if (theSeries7 != 0)
    delete theSeries7;
  if (theSeries8 != 0)
    delete theSeries8;
  if (theSeries9 != 0)
    delete theSeries9;
}

//PK input from tcl for this version (2.1)
void
FireLoadPattern::setFireTimeSeries(TimeSeries *theTimeSeries1, TimeSeries *theTimeSeries2, 
				   TimeSeries *theTimeSeries3, TimeSeries *theTimeSeries4, 
				   TimeSeries *theTimeSeries5, TimeSeries *theTimeSeries6, 
				   TimeSeries *theTimeSeries7, TimeSeries *theTimeSeries8, 
				   TimeSeries *theTimeSeries9)
{
  // invoke the destructor on the old TimeSeries
  if (theSeries1 != 0)
    delete theSeries1;
  if (theSeries2 != 0)
    delete theSeries2;
  if (theSeries3 != 0)
    delete theSeries3;
  if (theSeries4 != 0)
    delete theSeries4;
  if (theSeries5 != 0)
    delete theSeries5;
  if (theSeries6 != 0)
    delete theSeries6;
  if (theSeries7 != 0)
    delete theSeries7;
  if (theSeries8 != 0)
    delete theSeries8;
  if (theSeries9 != 0)
    delete theSeries9;
  
  // set the pointer to the new series objects
  theSeries1 = theTimeSeries1;
  theSeries2 = theTimeSeries2;
  theSeries3 = theTimeSeries3;
  theSeries4 = theTimeSeries4;
  theSeries5 = theTimeSeries5;
  theSeries6 = theTimeSeries6;
  theSeries7 = theTimeSeries7;
  theSeries8 = theTimeSeries8;
  theSeries9 = theTimeSeries9;
  //opserr << "FireLoadPattern set firetimeseries called\n";
}

void 
FireLoadPattern::applyLoad(double time)
{
	// first determine the load factor
  if (theSeries1 != 0 && isConstant != 0) {
    loadFactors(0) = theSeries1->getFactor(time);
    loadFactors(1) = theSeries2->getFactor(time);
    loadFactors(2) = theSeries3->getFactor(time);
    loadFactors(3) = theSeries4->getFactor(time);
    loadFactors(4) = theSeries5->getFactor(time);
    loadFactors(5) = theSeries6->getFactor(time);
    loadFactors(6) = theSeries7->getFactor(time);
    loadFactors(7) = theSeries8->getFactor(time);
    loadFactors(8) = theSeries9->getFactor(time);
  }
  //this is a fire load pattern so we always need multiple timeseries
	NodalLoad *nodLoad;
	NodalLoadIter &theNodalIter = this->getNodalLoads();
	while ((nodLoad = theNodalIter()) != 0)
		 nodLoad->applyLoad(loadFactors);

  ElementalLoad *eleLoad;
  ElementalLoadIter &theElementalIter = this->getElementalLoads();
  while ((eleLoad = theElementalIter()) != 0)
    eleLoad->applyLoad(loadFactors);

}

bool
FireLoadPattern::addSP_Constraint(SP_Constraint *)
{
  opserr << "FireLoadPattern::addSP_Constraint() - cannot add SP_Constraint to FireLoadPattern\n";
  return false;
}
/*
bool
FireLoadPattern::addNodalLoad(NodalLoad *)
{
  opserr << "FireLoadPattern::addNodalLoad() - cannot add NodalLoad to FireLoadPattern\n";  
  return false;
}
*/

//***************************************************************************************** */

void 
FireLoadPattern::Print(OPS_Stream &s, int flag)
{
  opserr << "FireLoadPattern::Print() - not yet implemented\n";    
}

