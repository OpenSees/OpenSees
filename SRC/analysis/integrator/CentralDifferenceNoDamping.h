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
                                                                        
#ifndef CentralDifferenceNoDamping_h
#define CentralDifferenceNoDamping_h

// Written: fmk 
// Created: 02/05
// Revision: A
//
// Description: This file contains the class definition for CentralDifferenceNoDamping.
// CentralDifferenceNoDamping is an algorithmic class for performing a transient 
// analysis using the Central Difference Scheme as implemented in Dyna
//       An = M(-1) (Pn - Fn)
//       Vn+1/2 = Vn-1/2 + dT * An
//       Dn+1   = Dn + deltaT * Vn+1/2
// which is an explicit direct integration scheme as outlined in the paper:
// Goudreau, G.L. and J.O. Hallquist, "Recent Developments in Large Scale Finite Element Lagrangian 
// Hydrocode Technology", Journal of Computer Methods in Applied Mechanics and Engineering, 30, 1982.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class CentralDifferenceNoDamping : public TransientIntegrator
{
  public:
    CentralDifferenceNoDamping();
    ~CentralDifferenceNoDamping();

    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);        
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);    

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
    Vector *U;          // disp response quantities at time t + deltaT
    Vector *Udot;       // vel response quantity at time t-1/2 delta t
    Vector *Udotdot;    // accel response at time t
    double deltaT;
};

#endif

