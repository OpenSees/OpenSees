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
#ifndef PML3DGENERAL_H
#define PML3DGENERAL_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>
#include <Eigen/Dense>


#define PML3DGeneral_NUM_DOF 72
#define PML3DGeneral_NUM_PROPS 12
#define PML3DGeneral_NUM_NODES 8
class PML3DGeneral : public Element {

public:
	PML3DGeneral();                                                                  //null constructor
	PML3DGeneral(int    tag, int* nodeTags,
				 double _E, double _nu, double _rho,
                 double _gamma, double _beta, double _eta,
                 double _m_pml, double _L_pml, double _R_pml, 
                 double _x0_pml, double _y0_pml, double _z0_pml,
                 double _nx_pml, double _ny_pml, double _nz_pml);                    // full constructor
	virtual ~PML3DGeneral();                                                         //destructor
	const char* getClassType(void) const { return "PML3DGeneral"; };                 //return class type
	void setDomain(Domain* theDomain);                                               // set domain
	int getNumExternalNodes() const; 	   						                     // get number of external nodes
	const ID& getExternalNodes(); 								                     // get external nodes
	Node** getNodePtrs(void); 									                     // get external nodes
	int getNumDOF(); 											                     // get number of DOF
	int commitState(); 											                     // commit state
	int revertToLastCommit(); 									                     // revert to last commit
	int revertToStart(); 										                     // revert to start
	int update(void); 											                     // update
	void Print(OPS_Stream& s, int flag); 							                 // print out element data
	const Matrix& getTangentStiff(); 							                     // get stiffness matrix
	const Matrix& getInitialStiff(); 							                     // get initial stiffness matrix
	const Matrix& getMass(); 									                     // get mass matrix
	const Matrix& getDamp(); 									                     // get damping matrix
	void zeroLoad(); 											                     // set residual to 0
	int addLoad(ElementalLoad* theLoad, double loadFactor); 		                 // add element loads
	const Vector& getResistingForce(); 							                     // get residual
	const Vector& getResistingForceIncInertia(); 				                     // get residual including damping forces
	int sendSelf(int commitTag, Channel& theChannel); 			                     // send self
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker & theBroker);  // receive self
	Response* setResponse(const char** argv, int argc, OPS_Stream& s);               // set response
	int getResponse(int responseID, Information& eleInformation);                    // get response
	int setParameter(const char** argv, int argc, Parameter& param);                 // set parameter   
	int updateParameter(int parameterID, Information& info);                         // update parameter
	int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);
	Eigen::MatrixXd ComputeJacobianMatrix(const double ri, const double si, const double ti) const;
	Eigen::MatrixXd ComputeShapeFunctionMatrix(const double ri, const double si, const double ti) const;
	Eigen::MatrixXd ComputeStrainDisplacementMatrix(const double ri, const double si, const double ti, const Eigen::MatrixXd &Jij) const;
	Eigen::VectorXd ComputePMLStretchingFactors(const double ri, const double si, const double ti, const double rho, const double mu, const double lambda) const;
	void ComputePMLMatrix();
	void ComputeStiffnessMatrix();
	void ComputeMassMatrix();
	void ComputeDampingMatrix();

private:
	Domain* Domainptr;                                            // pointer to the domain
	double props[PML3DGeneral_NUM_PROPS];                         // material properties
	ID connectedExternalNodes;  					              // eight node numbers
	Node* nodePointers[PML3DGeneral_NUM_NODES];    	              // pointers to eight nodes
	static double eta;                                            // Newmark parameters: eta
	static double beta; 					  	                  // Newmark parameters: beta
	static double gamma; 					  	                  // Newmark parameters: gamma
	Matrix ImpedanceMatrix;								   // impedance matrix
	Matrix StiffnessMatrix;                                // tangent matrix
	Matrix MassMatrix;						               // mass matrix
	Matrix DampingMatrix;	 					           // damping matrix
	Matrix EffectiveStiffnessMatrix;					   // effective stiffness matrix
	static Vector resid; 						                   // residual vector
	Vector ubart; 				                                  // ubar at time t 
	Vector ubar; 				                                  // ubar at time t+dt
	static double dt; 								              // time step
	static int eleCount; 						                  // element count
    double m_pml;                                                 // PML parameter m 
    double L_pml;                                                 // PML parameter L
    double R_pml;                                                 // PML parameter R
    double x0_pml;                                                // PML parameter x0 (x coordinate of the reference point)
    double y0_pml;                                                // PML parameter y0 (y coordinate of the reference point)
    double z0_pml;                                                // PML parameter z0 (z coordinate of the reference point)
    double nx_pml;                                                // PML parameter nx (x component of the normal vector)
    double ny_pml;                                                // PML parameter ny (y component of the normal vector)
    double nz_pml;                                                // PML parameter nz (z component of the normal vector)
	double E;                                                     // Young's modulus
	double nu;                                                    // Poisson's ratio
	double rho;                                                   // density
};

#endif
