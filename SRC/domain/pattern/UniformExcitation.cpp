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
                                                                        
// $Revision$
// $Date$
// $URL$
                                                                        
                                                                        
// File: ~/domain/load/UniformExcitation.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for UniformExcitation.
// UniformExcitation is an abstract class.

#include <UniformExcitation.h>
#include <GroundMotion.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SP_ConstraintIter.h>
#include <SP_Constraint.h>
#include <elementAPI.h>

void* OPS_TimeSeriesIntegrator();

void* OPS_UniformExcitationPattern()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING insufficient args : pattern UniformExcitation tag dir\n";
	return 0;
    }
    
    int patternID;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &patternID) < 0) {
	opserr << "WARNING invalid patternID\n";
	return 0;
    }
    
    int dir;
    if (OPS_GetIntInput(&numdata, &dir) < 0) {
	opserr << "WARNING invalid dir \n";
	return 0;
    }
    
    dir--; // subtract 1 for c indexing
    
    TimeSeries *accelSeries = 0;
    TimeSeries *velSeries = 0;
    TimeSeries *dispSeries = 0;
    TimeSeriesIntegrator *seriesIntegrator = 0;
    double vel0 = 0.0;
    double fact = 1.0;
    
    bool doneSeries = false;
    while (OPS_GetNumRemainingInputArgs()>1 && doneSeries == false) {

	const char* flag = OPS_GetString();
      
	if ((strcmp(flag,"-vel0") == 0) || (strcmp(flag,"-initialVel") == 0)) {
	
	    if (OPS_GetDoubleInput(&numdata, &vel0) < 0) {
		opserr << "WARNING invalid vel0: pattern type UniformExcitation\n";
		return 0;
	    }
	}
      
	else if ((strcmp(flag,"-fact") == 0) || (strcmp(flag,"-factor") == 0)) {
	
	    if (OPS_GetDoubleInput(&numdata, &fact) < 0) {
		opserr << "WARNING invalid fact: pattern type UniformExcitation\n";
		return 0;
	    }
	}
      
      
	else if ((strcmp(flag,"-accel") == 0) || (strcmp(flag,"-acceleration") == 0)) {

	    int tsTag;
	    if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		opserr << "WARNING invalid accel series tag\n";
		return 0;
	    }
	
	    accelSeries = OPS_getTimeSeries(tsTag);
	
	    if (accelSeries == 0) {
		opserr << "WARNING invalid accel series: " << tsTag;
		opserr << " pattern UniformExcitation -accel {series}\n";
		return 0;
	    }
	    
	} else if ((strcmp(flag,"-vel") == 0) || (strcmp(flag,"-velocity") == 0)) {

	    int tsTag;
	    if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		opserr << "WARNING invalid vel series tag\n";
		return 0;
	    }
	    velSeries = OPS_getTimeSeries(tsTag);
	
	    if (velSeries == 0) {
		opserr << "WARNING invalid vel series: " << tsTag;
		opserr << " pattern UniformExcitation -vel {series}\n";
		return 0;
	    }
	
	} else if ((strcmp(flag,"-disp") == 0) || (strcmp(flag,"-displacement") == 0)) {

	    int tsTag;
	    if (OPS_GetIntInput(&numdata, &tsTag) < 0) {
		opserr << "WARNING invalid disp series tag\n";
		return 0;
	    }
	
	    dispSeries = OPS_getTimeSeries(tsTag);
	
	    if (dispSeries == 0) {
		opserr << "WARNING invalid disp series: " << tsTag;
		opserr << " pattern UniformExcitation -disp {series}\n";
		return 0;
	    }
	
	} else if ((strcmp(flag,"-int") == 0) || (strcmp(flag,"-integrator") == 0)) {

	    seriesIntegrator = (TimeSeriesIntegrator*) OPS_TimeSeriesIntegrator();
	    if (seriesIntegrator == 0) return 0;
	}
    
	else 
	    doneSeries = true;
    }
    
    if (dispSeries == 0 && velSeries == 0 && accelSeries == 0) {
	opserr << "WARNING invalid series, want - pattern UniformExcitation";
	opserr << "-disp {dispSeries} -vel {velSeries} -accel {accelSeries} ";
	opserr << "-int {Series Integrator}\n";
	return 0;
    }
    
    GroundMotion *theMotion = new GroundMotion(dispSeries, velSeries,
					       accelSeries, seriesIntegrator);
    
    if (theMotion == 0) {
	opserr << "WARNING ran out of memory creating ground motion - pattern UniformExcitation ";
	opserr << patternID << endln;
      
	return 0;
    }
    
    // create the UniformExcitation Pattern
    UniformExcitation* thePattern = new UniformExcitation(*theMotion, dir, patternID, vel0, fact);
    
    if (thePattern == 0) {
	opserr << "WARNING ran out of memory creating load pattern - pattern UniformExcitation ";
	opserr << patternID << endln;
      
	// clean up memory allocated up to this point and return an error
	if (theMotion != 0)
	    delete theMotion;
      
	return 0;
    }

    return thePattern;
}

UniformExcitation::UniformExcitation()
:EarthquakePattern(0, PATTERN_TAG_UniformExcitation), 
 theMotion(0), theDof(0), vel0(0.0), fact(0.0)
{

}


UniformExcitation::UniformExcitation(GroundMotion &_theMotion, 
				     int dof, int tag, double velZero, double theFactor)
:EarthquakePattern(tag, PATTERN_TAG_UniformExcitation), 
 theMotion(&_theMotion), theDof(dof), vel0(velZero), fact(theFactor)
{
  // add the motion to the list of ground motions
  this->addMotion(*theMotion);
}


UniformExcitation::~UniformExcitation()
{

}


const GroundMotion *
UniformExcitation::getGroundMotion(void)
{
  return theMotion;
}

int
UniformExcitation::setParameter(const char **argv, int argc, Parameter &param)
{
  return theMotion->setParameter(argv, argc, param);
}

/*
int
UniformExcitation::updateParameter(int parameterID, Information &info)
{
  return theMotion->updateParameter(parameterID, info);
}

int
UniformExcitation::activateParameter(int pparameterID)
{
  return theMotion->activateParameter(pparameterID);
}
*/

void
UniformExcitation::setDomain(Domain *theDomain) 
{
  this->LoadPattern::setDomain(theDomain);

  // now we go through and set all the node velocities to be vel0 
  // for those nodes not fixed in the dirn!
  if (vel0 != 0.0) {

    SP_ConstraintIter &theSPs = theDomain->getSPs();
    SP_Constraint *theSP;
    ID constrainedNodes(0);
    int count = 0;
    while ((theSP=theSPs()) != 0) {
      if (theSP->getDOF_Number() == theDof) {
	constrainedNodes[count] = theSP->getNodeTag();
	count++;
      }
    }


    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    Vector newVel(1);
    int currentSize = 1;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      if (constrainedNodes.getLocation(tag) < 0) {
	int numDOF = theNode->getNumberDOF();
	if (numDOF != currentSize) 
	  newVel.resize(numDOF);
	
	newVel = theNode->getVel();
	newVel(theDof) = vel0;
	theNode->setTrialVel(newVel);
	theNode->commitState();
      }
    }
  }
}

void
UniformExcitation::applyLoad(double time)
{
    Domain *theDomain = this->getDomain();
    if (theDomain == 0)
        return;
    
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
        theNode->setNumColR(1);
        const Vector &crds=theNode->getCrds();
        int ndm = crds.Size();
        
        if (ndm == 1) {
            theNode->setR(theDof, 0, fact);
        }
        else if (ndm == 2) {
            if (theDof < 2) {
                theNode->setR(theDof, 0, fact);
            }
            else if (theDof == 2) {
                double xCrd = crds(0);
                double yCrd = crds(1);
                theNode->setR(0, 0, -fact*yCrd);
                theNode->setR(1, 0, fact*xCrd);
                theNode->setR(2, 0, fact);
            }
        }
        else if (ndm == 3) {
            if (theDof < 3) {
                theNode->setR(theDof, 0, fact);
            }
            else if (theDof == 3) {
                double yCrd = crds(1);
                double zCrd = crds(2);
                theNode->setR(1, 0, -fact*zCrd);
                theNode->setR(2, 0, fact*yCrd);
                theNode->setR(3, 0, fact);
            }
            else if (theDof == 4) {
                double xCrd = crds(0);
                double zCrd = crds(2);
                theNode->setR(0, 0, fact*zCrd);
                theNode->setR(2, 0, -fact*xCrd);
                theNode->setR(4, 0, fact);
            }
            else if (theDof == 5) {
                double xCrd = crds(0);
                double yCrd = crds(1);
                theNode->setR(0, 0, -fact*yCrd);
                theNode->setR(1, 0, fact*xCrd);
                theNode->setR(5, 0, fact);
            }
        }
    }
    
    this->EarthquakePattern::applyLoad(time);
    
    return;
}


void
UniformExcitation::applyLoadSensitivity(double time)
{
  Domain *theDomain = this->getDomain();
  if (theDomain == 0)
    return;

//  if (numNodes != theDomain->getNumNodes()) {
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      theNode->setNumColR(1);
      theNode->setR(theDof, 0, 1.0);
    }
//  }

  this->EarthquakePattern::applyLoadSensitivity(time);

  return;
}



int 
UniformExcitation::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector data(6);
  data(0) = this->getTag();
  data(1) = theDof;
  data(2) = vel0;
  data(5) = fact;
  data(3) = theMotion->getClassTag();
  
  int motionDbTag = theMotion->getDbTag();
  if (motionDbTag == 0) {
    motionDbTag = theChannel.getDbTag();
    theMotion->setDbTag(motionDbTag);
  }
  data(4) = motionDbTag;

  int res = theChannel.sendVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformExcitation::sendSelf() - channel failed to send data\n";
    return res;
  }
      
  res = theMotion->sendSelf(commitTag, theChannel);
  if (res < 0) {
    opserr << "UniformExcitation::sendSelf() - ground motion to send self\n";
    return res;
  }

  return 0;
}


int 
UniformExcitation::recvSelf(int commitTag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector data(6);
  int res = theChannel.recvVector(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "UniformExcitation::recvSelf() - channel failed to recv data\n";
    return res;
  }

  this->setTag(int(data(0)));
  theDof = int(data(1));
  vel0 = data(2);
  fact = data(5);
  int motionClassTag = int(data(3));
  int motionDbTag = int(data(4));

  if (theMotion == 0 || theMotion->getClassTag() != motionClassTag) {
    if (theMotion != 0)
      delete theMotion;
    theMotion = theBroker.getNewGroundMotion(motionClassTag);
    if (theMotion == 0) {
      opserr << "UniformExcitation::recvSelf() - could not create a grond motion\n";
      return -3;
    }

    // have to set the motion in EarthquakePattern base class
    if (numMotions == 0) 
      this->addMotion(*theMotion);
    else
      theMotions[0] = theMotion;
  }

  theMotion->setDbTag(motionDbTag);
  res = theMotion->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) {
      opserr << "UniformExcitation::recvSelf() - motion could not receive itself \n";
      return res;
  }

  return 0;
}


void 
UniformExcitation::Print(OPS_Stream &s, int flag)
{
  s << "UniformExcitation  " << this->getTag() << " - Not Printing the GroundMotion\n";
}

LoadPattern *
UniformExcitation::getCopy(void)
{
  LoadPattern *theCopy = new UniformExcitation(*theMotion, theDof, this->getTag());
   return theCopy;
}
//  LocalWords:  OpenSees
