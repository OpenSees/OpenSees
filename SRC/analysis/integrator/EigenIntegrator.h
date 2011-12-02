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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-14 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/EigenIntegrator.h,v $
                                                                        
                                                                        
// File: ~/analysis/integrator/eigenIntegrator/EigenIntegrator.h
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenIntegrator.
// EigenIntegrator is an algorithmic class for setting up the finite element 
// equations for a eigen problem analysis. 
//
// This class is inheritanted from the base class of Integrator which was
// created by fmk (Frank).


#ifndef EigenIntegrator_h
#define EigenIntegrator_h

#include <Integrator.h>

class EigenSOE;
class AnalysisModel;
class FE_Element;
class DOF_Group;
class Vector;

class EigenIntegrator : public Integrator
{
  public:
     EigenIntegrator();
     virtual ~EigenIntegrator();
     
     virtual void setLinks(AnalysisModel &theModel,
			   EigenSOE &theSOE);
       
     // methods to form the M and K matrices.
     virtual int formK();
     virtual int formM();
     
     // methods to instruct the FE_Element and DOF_Group objects
     // how to determing their contribution to M and K
     virtual int formEleTangK(FE_Element *theEle);
     virtual int formEleTangM(FE_Element *theEle);
     virtual int formNodTangM(DOF_Group *theDof);
     virtual int update(const Vector &deltaU);

     virtual int formEleTangent(FE_Element *theEle);
     virtual int formNodTangent(DOF_Group *theDof);
     virtual int formEleResidual(FE_Element *theEle);
     virtual int formNodUnbalance(DOF_Group *theDof);

     virtual int newStep(void);
     
     virtual int getLastResponse(Vector &result, const ID &id);

     virtual int sendSelf(int commitTag, Channel &theChannel);
     virtual int recvSelf(int commitTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker);     
     virtual void Print(OPS_Stream &s, int flag = 0);
     
 protected:
     virtual EigenSOE *getEigenSOEPtr() const;
     virtual AnalysisModel *getAnalysisModelPtr() const;
  
 private:
     EigenSOE *theSOE;
     AnalysisModel *theAnalysisModel;
     int flagK;

};

#endif




