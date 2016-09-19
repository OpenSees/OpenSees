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
                                                                        
// $Revision: 1.14 $
// $Date: 2006-02-08 20:20:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/FE_Element.h,v $
                                                                        
                                                                        
#ifndef FE_Element_h
#define FE_Element_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for FE_Element.
//
// What: "@(#) FE_Element.h, revA"

#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <TaggedObject.h>

class TransientIntegrator;
class Element;
class Integrator;
class AnalysisModel;

class FE_Element: public TaggedObject
{
  public:
    FE_Element(int tag, Element *theElement);
    FE_Element(int tag, int numDOF_Group, int ndof);
    virtual ~FE_Element();    

    // public methods for setting/obtaining mapping information
    virtual const ID &getDOFtags(void) const;
    virtual const ID &getID(void) const;
    void setAnalysisModel(AnalysisModel &theModel);
    virtual int  setID(void);
    
    // methods to form and obtain the tangent and residual
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual const Vector &getResidual(Integrator *theIntegrator);

    // methods to allow integrator to build tangent
    virtual void  zeroTangent(void);
    virtual void  addKtToTang(double fact = 1.0);
    virtual void  addKiToTang(double fact = 1.0);
    virtual void  addKgToTang(double fact = 1.0);
    virtual void  addCtoTang(double fact = 1.0);    
    virtual void  addMtoTang(double fact = 1.0);    
    virtual void  addKpToTang(double fact = 1.0, int numP = 0);
    virtual int   storePreviousK(int numP);
    
    // methods to allow integrator to build residual    
    virtual void  zeroResidual(void);    
    virtual void  addRtoResidual(double fact = 1.0);
    virtual void  addRIncInertiaToResidual(double fact = 1.0);    

    // methods for ele-by-ele strategies
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
    virtual const Vector &getK_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);
    virtual void  addM_Force(const Vector &accel, double fact = 1.0);    
    virtual void  addD_Force(const Vector &vel, double fact = 1.0);    
    virtual void  addK_Force(const Vector &disp, double fact = 1.0);
    virtual void  addKg_Force(const Vector &disp, double fact = 1.0);    

    virtual int updateElement(void);

    virtual Integrator *getLastIntegrator(void);
    virtual const Vector &getLastResponse(void);
    Element *getElement(void);

    virtual void  Print(OPS_Stream&, int = 0) {return;};

    // AddingSensitivity:BEGIN ////////////////////////////////////
    virtual void addResistingForceSensitivity(int gradNumber, double fact = 1.0);
    virtual void addM_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual void addD_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual int  commitSensitivity           (int gradNum, int numGrads);
    // AddingSensitivity:END //////////////////////////////////////
    
  protected:
    void  addLocalM_Force(const Vector &accel, double fact = 1.0);    
    void  addLocalD_Force(const Vector &vel, double fact = 1.0);    
    void  addLocalM_ForceSensitivity(int gradNumber, const Vector &accel, double fact = 1.0);    
    void  addLocalD_ForceSensitivity(int gradNumber, const Vector &vel, double fact = 1.0);    


    // protected variables - a copy for each object of the class        
    ID myDOF_Groups;
    ID myID;

  private:
    // private variables - a copy for each object of the class    
    int numDOF;
    AnalysisModel *theModel;
    Element *myEle;
    Vector *theResidual;
    Matrix *theTangent;
    Integrator *theIntegrator; // need for Subdomain

    
    // static variables - single copy for all objects of the class	
    static Matrix errMatrix;
    static Vector errVector;
    static Matrix **theMatrices; // array of pointers to class wide matrices
    static Vector **theVectors;  // array of pointers to class widde vectors
    static int numFEs;           // number of objects
    

};

#endif


