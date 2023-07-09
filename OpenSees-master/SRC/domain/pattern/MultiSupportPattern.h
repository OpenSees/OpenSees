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
// $Date: 2005-11-22 19:44:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/MultiSupportPattern.h,v $
                                                                        
#ifndef MultiSupportPattern_h
#define MultiSupportPattern_h

// Written: fmk 11/00
// Revised:
//
// Purpose: This file contains the class definition for MultiSupportPattern.
// MultiSupportPattern is an abstract class.

#include <LoadPattern.h>
#include <ID.h>

class GroundMotion;
class Vector;

class MultiSupportPattern : public LoadPattern
{
  public:
    MultiSupportPattern(int tag, int classTag);
    MultiSupportPattern(int tag);    
    MultiSupportPattern();    
    virtual ~MultiSupportPattern();

    virtual void applyLoad(double time);
    virtual bool addNodalLoad(NodalLoad *);
    virtual bool addElementalLoad(ElementalLoad *);
    
    // methods for o/p
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    virtual void Print(OPS_Stream &s, int flag =0);        

    // method to obtain a blank copy of the LoadPattern
    virtual LoadPattern *getCopy(void);

    int addMotion(GroundMotion &theMotion, int tag);    
    GroundMotion *getMotion(int tag);        

 protected:

 private:
    GroundMotion **theMotions;
    ID theMotionTags;
    int numMotions;
    int dbMotions;
};

#endif
