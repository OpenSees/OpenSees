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
                                                                        
// $Revision: 1.10 $
// $Date: 2009-08-25 22:08:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/IncrementalIntegrator.h,v $
                                                                        
#ifndef IncrementalIntegrator_h
#define IncrementalIntegrator_h

// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A
//
// Description: This file contains the interface for IncrementalIntegrator. 
// IncrementalIntegrator is an algorithmic class for setting up the finite 
// element equations in an incremental analysis and for updating the nodal
// response quantities based on the values in the soln vector.
//
// What: "@(#) IncrementalIntegrator.h, revA"

#include <Integrator.h>

class LinearSOE;
class EigenSOE;
class AnalysisModel;
class ConvergenceTest;
class FE_Element;
class DOF_Group;
class Vector;

#define CURRENT_TANGENT 0
#define INITIAL_TANGENT 1
#define CURRENT_SECANT  2
#define INITIAL_THEN_CURRENT_TANGENT  3
#define NO_TANGENT  4
#define SECOND_TANGENT 5

class IncrementalIntegrator : public Integrator
{
  public:
    IncrementalIntegrator(int classTag);
    virtual ~IncrementalIntegrator();

    virtual void setLinks(AnalysisModel &theModel,
			  LinearSOE &theSOE,
			  ConvergenceTest *theTest);

    virtual void setEigenSOE(EigenSOE *theSOE);


    // methods to set up the system of equations
    virtual int  formTangent(int statusFlag = CURRENT_TANGENT);    
    virtual int  formUnbalance(void);

    // pure virtual methods to define the FE_ELe and DOF_Group contributions
    virtual int formEleTangent(FE_Element *theEle) =0;
    virtual int formNodTangent(DOF_Group *theDof) =0;    
    virtual int formEleResidual(FE_Element *theEle) =0;
    virtual int formNodUnbalance(DOF_Group *theDof) =0;    

    // methods to update the domain
    virtual int newStep(double deltaT);
    virtual int update(const Vector &deltaU) =0;
    virtual int commit(void);
    virtual int revertToLastStep(void);
    virtual int initialize(void);

    //    int setModalDampingFactors(const Vector &);
    int setupModal(const Vector *modalDampingValues);
    int addModalDampingForce(const Vector *modalDampingValues);
    int addModalDampingMatrix(const Vector *modalDampingValues);
    virtual double getCFactor(void);
    
    virtual const Vector &getVel(void);
    int doMv(const Vector &v, Vector &res);

// AddingSensitivity:BEGIN //////////////////////////////////
    virtual int revertToStart();
    virtual int formIndependentSensitivityLHS(int statusFlag = CURRENT_TANGENT);
// AddingSensitivity:END ////////////////////////////////////
    
    // method introduced for domain decomposition
    virtual int getLastResponse(Vector &result, const ID &id);
    
  protected:
    LinearSOE *getLinearSOE(void) const;
    AnalysisModel *getAnalysisModel(void) const;
    ConvergenceTest *getConvergenceTest(void) const;

    virtual int  formNodalUnbalance(void);        
    virtual int  formElementResidual(void);            
    int statusFlag;

    //    Vector *modalDampingValues;
    EigenSOE *theEigenSOE;
    double *eigenVectors;
    Vector *eigenValues;
    Vector *dampingForces;
    bool isDiagonal;
    double *diagMass;
    Vector *mV;
    Vector *tmpV1;
    Vector *tmpV2;
    
  private:
    LinearSOE *theSOE;
    AnalysisModel *theAnalysisModel;
    ConvergenceTest *theTest;

};

#endif

