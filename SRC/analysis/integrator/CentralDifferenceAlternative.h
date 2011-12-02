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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-01-27 04:32:07 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/CentralDifferenceAlternative.h,v $
                                                                        
#ifndef CentralDifferenceAlternative_h
#define CentralDifferenceAlternative_h

// Written: fmk 
// Created: 11/98
// Revision: A
//
// Description: This file contains the class definition for CentralDifferenceAlternative.
// CentralDifferenceAlternative is an algorithmic class for performing a transient 
// analysis using the alternative form of the Central Differenceintegration scheme, 
// which is an explicit direct integration scheme as outlined in the book 'Concepts
// and Applications of Finite Element Analysis' by Cook, Malkus & Plesha.
//
// What: "@(#) CentralDifferenceAlternative.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class CentralDifferenceAlternative : public TransientIntegrator
{
  public:
    CentralDifferenceAlternative();
    ~CentralDifferenceAlternative();

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
    int updateCount;    // method should only have one update per step
    Vector *Ut, *Utp1;  // disp response quantities at time t and t + deltaT
    Vector *Udot;      // vel response quantity at time t-1/2 delta t
    double deltaT;
};

#endif

