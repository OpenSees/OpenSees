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
                                                                        
                                                                        
                                                                        
#ifndef FireLoadPattern_h
#define FireLoadPattern_h

// File: ~/domain/pattern/EarthquakePattern.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for EarthquakePattern.
// EarthquakePattern is an abstract class.

//Modified by Panagistis Kotsoivnos,[University of Edinburgh]

#include <LoadPattern.h>
#include <DomainComponent.h>
#include <Vector.h>

class GroundMotion;
class NodalLoad;
class TimeSeries;
class ElementalLoad;
class SP_Constraint;
class NodalLoadIter;
class ElementalLoadIter;
class SingleDomSP_Iter;
class SP_ConstraintIter;
class TaggedObjectStorage;
class GroundMotion;

class FireLoadPattern : public LoadPattern
{
  public:
    FireLoadPattern(int tag);
    ~FireLoadPattern();
    FireLoadPattern(int tag, int classTag);
    
    void applyLoad(double time);
	
    bool addSP_Constraint(SP_Constraint *);

    
    void setFireTimeSeries(TimeSeries *theSeries1, TimeSeries *theSeries2, 
			   TimeSeries *theSeries3, TimeSeries *theSeries4, TimeSeries *theSeries5, 
			   TimeSeries *theSeries6, TimeSeries *theSeries7, TimeSeries *theSeries8, TimeSeries *theSeries9);
    
    // methods for o/p
    // int sendSelf(int commitTag, Channel &theChannel) =0;
    // int recvSelf(int commitTag, Channel &theChannel, 
    //		 FEM_ObjectBroker &theBroker) =0;
    void Print(OPS_Stream &s, int flag =0);        
    
    // method to obtain a blank copy of the LoadPattern
    //FireLoadPattern *getCopy(void) =0;
    
 protected:

  private:
    //PK start
    TimeSeries *theSeries1;
    TimeSeries *theSeries2;
    TimeSeries *theSeries3;
    TimeSeries *theSeries4;
    TimeSeries *theSeries5;
    TimeSeries *theSeries6;
    TimeSeries *theSeries7;
    TimeSeries *theSeries8;
    TimeSeries *theSeries9;
    Vector loadFactors;

    //PK end
    double currentTime;
};

#endif
