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
// $Date: 2005-02-28 20:39:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/CentralDifference.h,v $
                                                                        
                                                                        
#ifndef CentralDifference_h
#define CentralDifference_h

// Written: fmk 
// Created: 11/98
//
// Description: This file contains the class definition for CentralDifference.

// Written: fmk 
// Created: 11/98
//
// Description: This file contains the implementation of the CentralDifference class.
//
// equi at time n expressed by:
//   MAn + CVn + R(Un) = Rext n
//
// approximating Vn and An by:
//     Vn = 1/2*dT (Un+1 -Un)
//     An = 1/dT^2 (Un+1 - 2Un + Un-1)
//
// results in:
//   [(1/dT*dT) M + (1/2*dT) C] Un+1 = Rext n - R(Un) + [(1/dT*dT)M](2Un-Un-1) + [(2/dT)C] Un-1


// What: "@(#) CentralDifference.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class CentralDifference : public TransientIntegrator
{
  public:
    CentralDifference();
    ~CentralDifference();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        

    int domainChanged(void);    
    int newStep(double deltaT);    
    int update(const Vector &deltaU);

    int commit(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);        
    
  protected:
    
  private:
    int updateCount;            // method should only have one update per step
    double c1, c2, c3;          // some constants we need to keep
    Vector *Utm1;               // disp response quantity at time t - deltaT
    Vector *U, *Udot, *Udotdot; // response quantities at time t
    Vector *Y, *Z;              // garbage 
    double deltaT;
};

#endif

