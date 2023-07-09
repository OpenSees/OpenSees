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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-09-05 20:51:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/EarthquakePattern.h,v $
                                                                        
                                                                        
#ifndef EarthquakePattern_h
#define EarthquakePattern_h

// File: ~/domain/pattern/EarthquakePattern.h
//
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for EarthquakePattern.
// EarthquakePattern is an abstract class.

#include <LoadPattern.h>

class GroundMotion;
class Vector;

class EarthquakePattern : public LoadPattern
{
  public:
    EarthquakePattern(int tag, int classTag);
    virtual ~EarthquakePattern();

    virtual void applyLoad(double time);
    virtual bool addSP_Constraint(SP_Constraint *);
    virtual bool addNodalLoad(NodalLoad *);
    virtual bool addElementalLoad(ElementalLoad *);

    // methods for o/p
    virtual int sendSelf(int commitTag, Channel &theChannel) =0;
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker) =0;
    virtual void Print(OPS_Stream &s, int flag =0);        

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void) =0;
    
    // AddingSensitivity:BEGIN //////////////////////////////////////////
    virtual void applyLoadSensitivity(double pseudoTime = 0.0);
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int  updateParameter(int parameterID, Information &info);
    virtual int  activateParameter(int parameterID);
    // AddingSensitivity:END ///////////////////////////////////////////
    
 protected:
    int addMotion(GroundMotion &theMotion);
    GroundMotion **theMotions;
    int numMotions;

  private:
    Vector *uDotG, *uDotDotG;
    double currentTime;

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
