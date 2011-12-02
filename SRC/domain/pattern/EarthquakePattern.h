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
// $Date: 2000-09-15 08:23:19 $
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
    virtual void Print(ostream &s, int flag =0);        

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void) =0;
    
 protected:
    int addMotion(GroundMotion &theMotion);
    GroundMotion **theMotions;
    int numMotions;

  private:
    Vector *uDotG, *uDotDotG;
};

#endif
