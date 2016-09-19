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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/Integrator.h,v $
                                                                        
                                                                        
#ifndef Integrator_h
#define Integrator_h

// File: ~/analysis/integrator/Integrator.h
// 
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class interface for Integrator.
// Integrator is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) Integrator.h, revA"

#include <MovableObject.h>
#include <OPS_Globals.h>
class FE_Element;
class DOF_Group;
class Vector;
class ID;
class FEM_ObjectBroker;
class Matrix;

class Integrator: public MovableObject
{
public:
     Integrator(int classTag);
    virtual ~Integrator();
    
    virtual int domainChanged(void);
    
    virtual int formEleTangent(FE_Element *theEle) =0;
    virtual int formNodTangent(DOF_Group *theDof) =0;    
    virtual int formEleResidual(FE_Element *theEle) =0;
    virtual int formNodUnbalance(DOF_Group *theDof) =0;    

    // Methods provided for Domain Decomposition
    virtual int getLastResponse(Vector &result, const ID &id) =0;

    // Method provided for Output
    virtual void Print(OPS_Stream &s, int flag =0) =0;

    // Sensitivity Integrator interface
    virtual int formSensitivityRHS(int gradNum);
    virtual int formIndependentSensitivityRHS();
    virtual int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    virtual int commitSensitivity (int gradNum, int numGrads);
    ////////////////////////////////Abbas//////////////////
    virtual int formEleTangentSensitivity(FE_Element *theEle, int gradNumber);  
    virtual double getLambdaSensitivity(int gradNumber);
    virtual int computeSensitivities();//Abbas
    int sensitivityDomainChanged();//Abbass
    bool shouldComputeAtEachStep(void);
    bool newAlgorithm(void) {return true;};
    virtual  bool computeSensitivityAtEachIteration();
    bool activateSensitivityKey();
    bool activateSensitivity( ){return SensitivityKey;}; 
     ///////////////////////////////Abbas//////////////////

 protected:
private:


int analysisTypeTag;
bool SensitivityKey; // toactivate the sensitivity ind 

};

#endif






