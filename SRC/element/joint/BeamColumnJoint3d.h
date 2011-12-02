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
// $Date: 2003-06-11 18:20:12 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/BeamColumnJoint3d.h,v $
                                                                        
// Written: NM (nmitra@u.washington.edu)
// Created: Feb 2003
//
// Description: This file contains the class defination for beam-column joint.
// This element is a 4 noded 24 dof (6 dof at each node) finite area super-element, being a slight
// variation of the 2d one. The element takes in 13 different material types in order to simulate
// the inelastic action observed in a reinforced beam column joint. Though it has 6 dof per node 
// the out of the plane nodal dof are constrained or fixed and the inplane nodal dof are activated
                                                                        
#ifndef BeamColumnJoint3d_h
#define BeamColumnJoint3d_h


#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>


class Node;
class Channel;

class Response;
class Renderer;
class UniaxialMaterial;


class BeamColumnJoint3d : public Element
{
  public:
    // default constructor
	BeamColumnJoint3d(); 
	
	// defined constructor
    BeamColumnJoint3d(int tag,int Nd1, int Nd2, int Nd3, int Nd4,
		UniaxialMaterial& theMat1, UniaxialMaterial& theMat2,
		UniaxialMaterial& theMat3, UniaxialMaterial& theMat4,
		UniaxialMaterial& theMat5, UniaxialMaterial& theMat6,
		UniaxialMaterial& theMat7, UniaxialMaterial& theMat8,
		UniaxialMaterial& theMat9, UniaxialMaterial& theMat10,
		UniaxialMaterial& theMat11, UniaxialMaterial& theMat12,
		UniaxialMaterial& theMat13);
	
	// default destructor
	 ~BeamColumnJoint3d();

    ////////////// public methods to obtain information about dof & connectivity    
   bool	isSubdomain(void) { return false; } ;
    
	// get number of external nodes
	int getNumExternalNodes(void) const;

    // return connected external nodes
	const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    
	// return number of DOFs
	int getNumDOF(void);	
    
	// set domain performs check on dof and associativity with node
	void setDomain(Domain *theDomain);

    //////////////////////////// public methods to set the state of the element    
    
	// commit state
	int commitState(void);
    
	// revert to last commit
	int revertToLastCommit(void);        
    
	// revert to start
	int revertToStart(void);        
    
	// determine current strain and set strain in material
	int update(void);
    
    //////////////////////// public methods to obtain stiffness, mass, damping and 
    ////////////////////////////////////// residual information    
    
	// returns converged tangent stiffness matrix
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);           

	// not required for this element formulation
	const Matrix &getDamp(void);    
    const Matrix &getMass(void);    

    // not required for this element formulation
	void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);
    
	// get converged residual
	const Vector &getResistingForce(void);
    
	// get converged residual with inertia terms
	const Vector &getResistingForceIncInertia(void);            

    // public methods for element output for parallel and database processing
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
	// display element graphically
	int displaySelf(Renderer &theViewer, int displayMode, float fact);    
    
	// print out element data
	void Print(OPS_Stream &s, int flag =0);    

    // implemented to print into file
	Response *setResponse(const char **argv, int argc, Information &eleInfo);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter (char **argv, int argc, Information &info);
    int updateParameter (int parameterID, Information &info);


  protected:
    
  private:
    
	// private methods
    void getGlobalDispls(Vector &dg) ;
    void getBCJoint();
	void getdg_df();
	void getdDef_du();
	void matDiag(Vector, Matrix&);
	void getMatResponse(Vector, Vector&, Vector&);
	void formR(Vector);
	void formK(Vector);
	void formTransfMat();
	double getStepSize(double,double,Vector,Vector,Vector,Vector,double);

    // material info
	UniaxialMaterial* MaterialPtr[13];  // pointer to the 13 different materials

	// node info
	ID  connectedExternalNodes;   // contains the tags of the end nodes
	Node* nodePtr[4];             // pointers to four nodes


	// various other element parameters
	double nodeCrd[4][3]; 
	double elemWidth;
	double elemHeight;

	Vector Uecommit;             // vector of external commited displacements
	Vector UeIntcommit;          // vector of internal commited displacements   
	Vector UeprCommit;           // vector of previous external committed displacements
	Vector UeprIntCommit;        // vector of previous internal committed displacements  
	Matrix BCJoint;       // matrix describing relation between the component deformations and the external and internal deformations
	Matrix dg_df;         // matrix of derivative of internal equilibrium 
    Matrix dDef_du;       // matrix of a portion of BCJoint reqd. for static condensation

    Matrix K;               // element stiffness matrix
    Vector R;               // element residual matrix
	
	// static transformation matrices
	static Matrix Transf;
	static Matrix Tran;

};

#endif

