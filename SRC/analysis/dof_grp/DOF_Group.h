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
                                                                        
// $Revision: 1.13 $
// $Date: 2007-02-14 18:44:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/dof_grp/DOF_Group.h,v $
                                                                        
                                                                        
#ifndef DOF_Group_h
#define DOF_Group_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the class definition for DOF_Group.
// A DOF_Group object is instantiated by the ConstraintHandler for 
// every unconstrained node in the domain. The constrained nodes require 
// specialised types of DOF_Group; which deal with the constraints. DOF_Group
// objects can handle 0 boundary constraints; if the eqn number of a DOF is 
// less than START_EQN_NUM a value of 0.0 is set for disp, vel and accel when
// a setNode*(Vector &) is invoked.
//
// What: "@(#) DOF_Group.h, revA"

#include <ID.h>
#include <TaggedObject.h>

class Node;
class Vector;
class Matrix;
class TransientIntegrator;
class Integrator;

class DOF_Group: public TaggedObject
{
  public:
    DOF_Group(int tag, Node *myNode);
    DOF_Group(int tag, int ndof);    
    virtual ~DOF_Group();    

    virtual void setID(int dof, int value);
    virtual void setID(const ID &values);
    virtual const ID &getID(void) const;
    virtual int doneID(void);    

    virtual int getNodeTag(void) const;
    virtual int getNumDOF(void) const;    
    virtual int getNumFreeDOF(void) const;
    virtual int getNumConstrainedDOF(void) const;

    // methods to form the tangent
    virtual const Matrix &getTangent(Integrator *theIntegrator);
    virtual void  zeroTangent(void);
    virtual void  addMtoTang(double fact = 1.0);    
    virtual void  addCtoTang(double fact = 1.0);    

    // methods to form the unbalance
    virtual const Vector &getUnbalance(Integrator *theIntegrator);
    virtual void  zeroUnbalance(void);
    virtual void  addPtoUnbalance(double fact = 1.0);
    virtual void  addPIncInertiaToUnbalance(double fact = 1.0);    
    virtual void  addM_Force(const Vector &Udotdot, double fact = 1.0);        

    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);

    // methods to obtain committed responses from the nodes
    virtual const Vector & getCommittedDisp(void);
    virtual const Vector & getCommittedVel(void);
    virtual const Vector & getCommittedAccel(void);
    
    // methods to update the trial response at the nodes
    virtual void setNodeDisp(const Vector &u);
    virtual void setNodeVel(const Vector &udot);
    virtual void setNodeAccel(const Vector &udotdot);

    virtual void incrNodeDisp(const Vector &u);
    virtual void incrNodeVel(const Vector &udot);
    virtual void incrNodeAccel(const Vector &udotdot);

    // methods to set the eigen vectors
    virtual void setEigenvector(int mode, const Vector &eigenvalue);
    virtual const Matrix &getEigenvectors(void);

    virtual double getDampingBetaFactor(int mode, double ratio, double wn);
    virtual const Vector &getDampingBetaForce(int mode, double beta);

	
    // method added for TransformationDOF_Groups
    virtual Matrix *getT(void);

// AddingSensitivity:BEGIN ////////////////////////////////////
    virtual void addM_ForceSensitivity(const Vector &Udotdot, double fact = 1.0);        
    virtual void addD_ForceSensitivity(const Vector &vel, double fact = 1.0);
    virtual void addD_Force(const Vector &vel, double fact = 1.0);

    virtual const Vector & getDispSensitivity(int gradNumber);
    virtual const Vector & getVelSensitivity(int gradNumber);
    virtual const Vector & getAccSensitivity(int gradNumber);
    virtual int saveDispSensitivity(const Vector &v,
				    int gradNum, int numGrads);
    virtual int saveVelSensitivity(const Vector &vdot,
				   int gradNum, int numGrads);
    virtual int saveAccSensitivity(const Vector &vdotdot,
				   int gradNum, int numGrads);
    virtual int saveSensitivity(const Vector &v, const Vector &vdot,
				const Vector &vdotdot, int gradNum, int numGrads);
// AddingSensitivity:END //////////////////////////////////////
    virtual void  Print(OPS_Stream&, int = 0) {return;};
    virtual void resetNodePtr(void);
  
   protected:
    void  addLocalM_Force(const Vector &Udotdot, double fact = 1.0);     

    // protected variables - a copy for each object of the class            
    Vector *unbalance;
    Matrix *tangent;
    Node *myNode;
    
  private:
    // private variables - a copy for each object of the class        
    ID 	myID;
    int numDOF;

    // static variables - single copy for all objects of the class	    
    static Matrix errMatrix;
    static Vector errVect;
    static Matrix **theMatrices; // array of pointers to class wide matrices
    static Vector **theVectors;  // array of pointers to class widde vectors
    static int numDOFs;           // number of objects    
};

#endif

