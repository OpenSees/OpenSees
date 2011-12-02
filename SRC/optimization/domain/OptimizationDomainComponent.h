/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
**   Optimization module developed by:                                **
**   Quan Gu  (qgu@ucsd.edu)                                          **
**   Joel Conte (jpconte@ucsd.edu)                                    **
**   Philip Gill (pgill@ucsd.edu)                                     **
** ****************************************************************** */


//
// Written by Quan Gu (qgu@ucsd.edu)    March 2010
//

#ifndef OptimizationDomainComponent_h
#define OptimizationDomainComponent_h


#include <TaggedObject.h>
#include <MovableObject.h>

 
//#include <OptimizationDomain.h>
#include <Vector.h>


class OptimizationDomain;
      
class OptimizationDomainComponent: public TaggedObject, public MovableObject
{
//
//   #pragma message("OptimizationDomainComponent is included")
//
public:
    virtual ~OptimizationDomainComponent();
	virtual void setOptimizationDomain(OptimizationDomain *theOptimizationDomain);
    virtual OptimizationDomain *getOptimizationDomain(void) const;
    virtual void Print(OPS_Stream &s, int flag =0) =0;       

protected:
    OptimizationDomainComponent(int tag, int classTag);

    
//private:    // to save some uncessary functions
 
    OptimizationDomain *theOptimizationDomain;
};

#endif


