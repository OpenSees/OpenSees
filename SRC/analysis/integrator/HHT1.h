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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/HHT1.h,v $
                                                                        
                                                                        
#ifndef HHT1_h
#define HHT1_h

// Written: fmk 
//
// Description: This file contains the class definition for HHT1.
// HHT1 is an algorithmic class for performing a transient analysis
// using the HHT1 integration scheme.
//
// What: "@(#) HHT1.h, revA"

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class HHT1 : public TransientIntegrator
{
  public:
    HHT1();
    HHT1(double alpha);
    HHT1(double alpha, double alphaM, double betaK, double betaKi, double betaKc);        

    ~HHT1();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        

    int domainChanged(void);    
    int initialize(void);        
    int newStep(double deltaT);    
    int update(const Vector &deltaU);

    int commit(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);        
    
  protected:
    
  private:
    double alpha;
    double gamma;
    double beta;

    double alphaM, betaK, betaKi, betaKc;

    double c1, c2, c3;  // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot; // response quantities at time t
    Vector *U, *Udot, *Udotdot; // response quantities at time t + delta t
    Vector *Ualpha, *Udotalpha; // response quantities at time t+alpha delta t
};

#endif

