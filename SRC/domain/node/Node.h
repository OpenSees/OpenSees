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
// $Date: 2008-08-19 22:50:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/Node.h,v $
                                                                        
                                                                        
#ifndef Node_h
#define Node_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for Node.
// A Node provides the abstraction for a node in the FEM.
// Nodes have original position, trial displacement, velocity and 
// acceleration and committed displacement, velocity and acceleration 
// (the last committed() trial quantities).
//
// What: "@(#) Node.h, revA"

#include <DomainComponent.h>


class Element;
class Vector;
class Matrix;
class Channel;
class Renderer;

class DOF_Group;
class NodalThermalAction; //L.Jiang [ SIF ]

class Node : public DomainComponent
{
  public:
    // constructors
    Node(int classTag);
    Node(int tag, int classTag);
    Node(int tag, int ndof, double Crd1, Vector *displayLoc = 0);
    Node(int tag, int ndof, double Crd1, double Crd2, Vector *displayLoc = 0);
    Node(int tag, int ndof, double Crd1, double Crd2, double Crd3, Vector *displayLoc = 0);
    Node(const Node &theCopy, bool copyMass = true);        
    
    // destructor
    virtual ~Node();

    // public methods dealing with the DOF at the node
    virtual int  getNumberDOF(void) const;    
    virtual void setDOF_GroupPtr(DOF_Group *theDOF_Grp);
    virtual DOF_Group *getDOF_GroupPtr(void);

    // public methods for obtaining the nodal coordinates
    virtual const Vector &getCrds(void) const;
    virtual int getDisplayCrds(Vector &results, double fact, int displayMode=0);
    virtual int setDisplayCrds(const Vector &theCrds);

    // public methods for obtaining committed and trial 
    // response quantities of the node
    virtual const Vector &getDisp(void);
    virtual const Vector &getVel(void);
    virtual const Vector &getAccel(void);    
    virtual const Vector &getIncrDisp(void);
    virtual const Vector &getIncrDeltaDisp(void);

    virtual const Vector &getTrialDisp(void);    
    virtual const Vector &getTrialVel(void);    
    virtual const Vector &getTrialAccel(void);    

    // public methods for updating the trial response quantities
    virtual int setTrialDisp(double value, int dof);    
    virtual int setTrialDisp(const Vector &);    
    virtual int setTrialVel(const Vector &);    
    virtual int setTrialAccel(const Vector &);        

    virtual int incrTrialDisp(const Vector &);    
    virtual int incrTrialVel(const Vector &);    
    virtual int incrTrialAccel(const Vector &);        
    
    // public methods for adding and obtaining load information
    virtual void zeroUnbalancedLoad(void);
    virtual int addUnbalancedLoad(const Vector &load, double fact = 1.0); 
    virtual int addInertiaLoadToUnbalance(const Vector &accel, double fact = 1.0);    
    virtual const Vector &getUnbalancedLoad(void);     
    virtual const Vector &getUnbalancedLoadIncInertia(void);        

    // public methods dealing with the committed state of the node
    virtual int commitState();
    virtual int revertToLastCommit();    
    virtual int revertToStart();        

    // public methods for dynamic analysis
    virtual const Matrix &getMass(void);
    virtual int setMass(const Matrix &theMass);
    virtual int setNumColR(int numCol);
    virtual int setR(int row, int col, double Value);
    virtual const Vector &getRV(const Vector &V);        

    virtual int setRayleighDampingFactor(double alphaM);
    virtual const Matrix &getDamp(void);

    // public methods for eigen vector
    virtual int setNumEigenvectors(int numVectorsToStore);
    virtual int setEigenvector(int mode, const Vector &eigenVector);
    virtual const Matrix &getEigenvectors(void);
    
    // public methods for output
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    virtual void Print(OPS_Stream &s, int flag = 0);
    virtual int displaySelf(Renderer &theRenderer, int displayMode, float fact);

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    int addInertiaLoadSensitivityToUnbalance(const Vector &accel, 
					     double fact = 1.0, 
					     bool tag=false);    
    Matrix getMassSensitivity(void);
    virtual const Matrix &getDampSensitivity(void);
    int    getCrdsSensitivity(void);
    int	   saveDispSensitivity(const Vector &v, int gradNum, int numGrads);
    int	   saveVelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    int	   saveAccelSensitivity(const Vector &vdot, int gradNum, int numGrads);
    double getDispSensitivity(int dof, int gradNum);
    double getVelSensitivity(int dof, int gradNum);
    double getAccSensitivity(int dof, int gradNum);
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter(int parameterID, Information &info);
    int    activateParameter(int parameterID);

    // AddingSensitivity:END ///////////////////////////////////////////

	//Add Pointer to NodalThermalAction id applicable------begin-----L.Jiang [ SIF ]
	virtual NodalThermalAction* getNodalThermalActionPtr(void);
	virtual void setNodalThermalActionPtr(NodalThermalAction* theAction);
	//Add Pointer to NodalThermalAction id applicable-----end------L.Jiang [ SIF ]


    virtual const Vector &getReaction();
    virtual int   addReactionForce(const Vector &, double factor);
    virtual int   resetReactionForce(int flag);
    virtual const Vector *getResponse(NodeResponseType);

    virtual void setCrds(double Crd1);
    virtual void setCrds(double Crd1, double Crd2);
    virtual void setCrds(double Crd1, double Crd2, double Crd3);
    virtual void setCrds(const Vector &);

  protected:

  private:
    // priavte methods used to create the Vector objects 
    // for the committed and trial response quantaties.
    int createDisp(void);
    int createVel(void);
    int createAccel(void); 

    // private method to set up global matrices
    int setGlobalMatrices();

    // private data associated with each node object
    int numberDOF;                    // number of dof at Node
    DOF_Group *theDOF_GroupPtr;       // pointer to associated DOF_Group
    Vector *Crd;                      // original nodal coords
    Vector *commitDisp, *commitVel, *commitAccel; // commited quantities
    Vector *trialDisp, *trialVel, *trialAccel;     // trial quantities
    Vector *unbalLoad;                // unbalanced load
    Vector *incrDisp;
    Vector *incrDeltaDisp;
    
    double *disp, *vel, *accel; // double arrays holding the displ, 
                                // vel and accel values

    int dbTag1, dbTag2, dbTag3, dbTag4; // needed for database
    Matrix *R;                          // nodal participation matrix
    Matrix *mass;                       // pointer to mass matrix
    Vector *unbalLoadWithInertia;       
    double alphaM;                    // rayleigh damping factor
    
    Matrix *theEigenvectors;

    // AddingSensitivity:BEGIN /////////////////////////////////////////
    Matrix *dispSensitivity;
    Matrix *velSensitivity;
    Matrix *accSensitivity;
    int parameterID;
    // AddingSensitivity:END ///////////////////////////////////////////

    NodalThermalAction *theNodalThermalActionPtr; //Added by Liming Jiang for pointer to nodalThermalAction, [SIF]

    static Matrix **theMatrices;
    static int numMatrices;
    static Matrix **theVectors;
    static int numVectors;
    int index;

    Vector *reaction;
    Vector *displayLocation;
};

#endif

