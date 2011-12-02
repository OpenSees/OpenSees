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
// $Date: 2003-02-22 01:02:06 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/transformation/TransformationFE.h,v $
                                                                        
                                                                        
#ifndef TransformationFE_h
#define TransformationFE_h

// File: ~/analysis/fe_ele/transformation/TransformationFE.h
// 
// Written: fmk 
// Created: 05/99
// Revision: A
//
// Description: This file contains the class definition for TransformationFE.
// TransformationFE objects handle MP_Constraints using the transformation
// method T^t K T. SP_Constraints are handled by the TransformationConstraintHandler.
//
// What: "@(#) TransformationFE.h, revA"

#include <FE_Element.h>
class SP_Constraint;
class DOF_Group;
class TransformationConstraintHandler;

class TransformationFE: public FE_Element
{
  public:
    TransformationFE(Element *theElement, 
		     TransformationConstraintHandler &theHandler);
    ~TransformationFE();    

    // public methods for setting/obtaining mapping information
    virtual const ID &getDOFtags(void) const;
    virtual const ID &getID(void) const;
    void setAnalysisModel(AnalysisModel &theModel);
    virtual int setID(void);
    
    // methods to form and obtain the tangent and residual
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual const Vector &getResidual(Integrator *theIntegrator);
    
    // methods for ele-by-ele strategies
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &accel, double fcat = 1.0);
    virtual const Vector &getD_Force(const Vector &vel, double fcat = 1.0);
    const Vector &getLastResponse(void);
    int addSP(SP_Constraint &theSP);
    
  protected:
    int transformResponse(const Vector &modResponse, Vector &unmodResponse);
    
  private:
    
    // private variables - a copy for each object of the class        
    DOF_Group **theDOFs;
    int numSPs;
    SP_Constraint **theSPs;
    ID *modID;
    Matrix *modTangent;
    Vector *modResidual;
    int numGroups;
    int numTransformedDOF;
    int numOriginalDOF;
    TransformationConstraintHandler *theHandler;
    
    // static variables - single copy for all objects of the class	
    static Matrix **modMatrices; // array of pointers to class wide matrices
    static Vector **modVectors;  // array of pointers to class widde vectors
    static Matrix **theTransformations; // for holding pointers to the T matrices
    static int numTransFE;     // number of objects    
    static int transCounter;   // a counter used to indicate when to do something
    static int sizeTransformations; // size of theTransformations array
    static double *dataBuffer;
    static double *localKbuffer;
    static int    *dofData;
    static int sizeBuffer;
};

#endif




