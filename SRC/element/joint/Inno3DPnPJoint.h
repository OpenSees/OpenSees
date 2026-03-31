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

// Code developed by: Cristian V. Miculaș  (github user name: cvmiculas)
// Element conceptualization: Cristian V. Miculaș (cristian.miculas@uc.pt), Ricardo J. Costa (rjcosta@dec.uc.pt) and Luís Simões da Silva (luisss@dec.uc.pt)
// Affiliation: Civil Engineering Department, ISISE, University of Coimbra, Portugal
// Acknowledgements: This work has been supported in part by national funds through FCT – Foundation for Science and Technology, Portugal, under grant agreement SFRH/BD/138151/2018 awarded to Cristian V. Miculaş

// Created: 04.02.2024
// Revised: 

// Description:
// This file contains the class implementation for the 3D beam-column-joint element object designed for the 3D innovative plug-and-play steel tubular joint configuration proposed within the INNO3DJOINTS project (https://ec.europa.eu/info/funding-tenders/opportunities/portal/screen/how-to-participate/org-details/960532413/project/749959/program/31061225/details).
// This element has 5 external nodes (6 DOFs/node) and 4 internal nodes (1 DOF/node), resulting in a total of 34 DOFs.
// This element can be viewed as a 2D plate in the 3D space.
// This element has 32 componenets (0D elements), each allowing for a different uniaxial material tag.

// References:
//
// 1. C.V. Miculaş, Innovative plug-and-play joints for hybrid tubular constructions (Ph.D. thesis), University of Coimbra, Portugal, 2023, https://estudogeral.uc.pt/handle/10316/110990
//
// 2. C. V. Miculaş, R. J. Costa, L. S. da Silva, R. Simões, H. Craveiro, T. Tankova, 3D macro-element for innovative plug-and-play joints, J. Constructional Steel Research 214 (2024), https://doi.org/10.1016/j.jcsr.2023.108436
//
// 3. C.V. Miculaş, R.J. Costa, L. Simões da Silva, R. Simões, H. Craveiro, T. Tankova, Macro-modelling of the three-dimensional interaction between the faces of a steel tubular column joint, in: F. Di Trapani, C. Demartino, G.C. Marano, G. Monti (Eds.), Proceedings of the 2022 Eurasian OpenSees Days, Springer Nature Switzerland, Cham, 2023, pp. 408–422, http://dx.doi.org/10.1007/978-3-031-30125-4_37
//

#ifndef Inno3DPnPJoint_h
#define Inno3DPnPJoint_h

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <FileStream.h>
#include <OPS_Stream.h>

class Node;
class Channel;
class FEM_ObjectBroker;
class Response;
class Renderer;
class UniaxialMaterial;

class Inno3DPnPJoint : public Element
{
  public:
    // default constructor
    Inno3DPnPJoint();
	
    // defined constructor
								

	Inno3DPnPJoint(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5,
				    UniaxialMaterial& theMat1,
				    UniaxialMaterial& theMat2,
				    UniaxialMaterial& theMat3,
				    UniaxialMaterial& theMat4,
				    UniaxialMaterial& theMat5,
				    UniaxialMaterial& theMat6,
				    UniaxialMaterial& theMat7,
				    UniaxialMaterial& theMat8,
				    UniaxialMaterial& theMat9,
				    UniaxialMaterial& theMat10,
				    UniaxialMaterial& theMat11,
					UniaxialMaterial& theMat12,
					UniaxialMaterial& theMat13,
					UniaxialMaterial& theMat14,
					UniaxialMaterial& theMat15,
					UniaxialMaterial& theMat16,
					UniaxialMaterial& theMat17,
					UniaxialMaterial& theMat18,
					UniaxialMaterial& theMat19,
					UniaxialMaterial& theMat20,
					UniaxialMaterial& theMat21,
					UniaxialMaterial& theMat22,
					UniaxialMaterial& theMat23,
					UniaxialMaterial& theMat24,
					UniaxialMaterial& theMat25,
					UniaxialMaterial& theMat26,
					UniaxialMaterial& theMat27,
					UniaxialMaterial& theMat28,
					UniaxialMaterial& theMat29,
					UniaxialMaterial& theMat30,
					UniaxialMaterial& theMat31,
					UniaxialMaterial& theMat32);

    // default destructor
    ~Inno3DPnPJoint();

	//####### public methods to obtain information about dof & connectivity
	
    bool isSubdomain(void) { return false; } ;
    
    // get number of external nodes
    int getNumExternalNodes(void) const;
    
    // return connected external nodes
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    
    // return number of DOFs
    int getNumDOF(void);    
    
    // set domain performs check on dof and associativity with node
    void setDomain(Domain *theDomain);
    
    //####### public methods to set the state of the element 
	
    // commit state
    int commitState(void);
    
    // revert to last commit
    int revertToLastCommit(void);                
    
    // revert to start
    int revertToStart(void);                
    
    // determine current strain and set strain in material
    int update(void);
	
	//####### public methods to obtain stiffness, mass, damping and
	//####### residual information  
	
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
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    // display element graphically
	int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    
    // print out element data
    void Print(OPS_Stream &s, int flag =0);        
    
    // implemented to print into file
    const char *getClassType(void) const {return "Inno3DPnPJoint";};
	
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter (char **argv, int argc, Information &info);
    int updateParameter (int parameterID, Information &info);
    
  
protected:

private:

	//private methods
    void getGlobalDispls(Vector&) ;
    void getmatA(void);
    void getdg_df(void);
    void getdDef_du(void);
    void matDiag(Vector, Matrix&);
    void getMatResponse(Vector, Vector&, Vector&);
    void formR(Vector);
    void formK(Vector);
	void formTransfMat();
    double getStepSize(double,double,Vector,Vector,Vector,Vector,double);

    // material info
    UniaxialMaterial **MaterialPtr;		// pointer to the 32 different materials

	// node info
	ID ExternalNodes;   				// contains the tags of the end nodes
	Node* nodePtr[5];         		    // pointers to four nodes
		
	int nodeDbTag, dofDbTag;
	
	// various other element parameters
	double dcX;
	double dcZ;
	
	Vector Uecommit;        // vector of external committed displacements
	Vector UeIntcommit;     // vector of internal committed displacements     
	Vector UeprCommit;      // vector of previous external committed displacements
	Vector UeprIntCommit;   // vector of previous internal committed displacements  
	Matrix matA;       		// matrix describing relation between the component deformations and the external and internal deformations
	Matrix dg_df;         	// matrix of derivative of internal equilibrium 
	Matrix dDef_du;       	// matrix of a portion of matA reqd. for static condensation
	
	Matrix K;             	// element stiffness matrix
	Vector R;             	// element residual matrix 

	// static transformation matrices
	static Matrix Transf;
	static Matrix Tran;	
};
#endif