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


// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu)
//
// Eight node PML3D element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#ifndef PML3DVISCOUS_H
#define PML3DVISCOUS_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <Domain.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <ElasticIsotropicMaterial.h>

#define PML3DVISCOUS_NUM_DOF 72
#define PML3DVISCOUS_NUM_PROPS 12
#define PML3DVISCOUS_NUM_NODES 8

#ifdef _WIN32

#define pml3d_	      PML_3D

extern "C" void  pml3d_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	double* gMatrix,
	double* hMatrix,
	int* NDOFEL,
	double* PROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE,
	int* LFLAGS);

#else

#define pml3d_	      pml_3d_

extern "C" void  pml3d_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	double* gMatrix,
	double* hMatrix,
	int* NDOFEL,
	double* PROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE,
	int* LFLAGS);

#endif

class PML3DVISCOUS : public Element {

public:

	PML3DVISCOUS();                                                                         //null constructor
	// PML3DVISCOUS(int tag, int* nodeTags, double* newmarks, double* dData);                  // full constructor
	PML3DVISCOUS(int tag, int* nodeTags,
				 NDMaterial* theMat, double PMLThickness,
				 double* Xref, double* Normal,
				 double alpha_0, double beta_0, bool explicitAB,
				 double Cp, double m_coeff, double R,
				 double gammaN, double betaN, double etaN, double keisiN);                  // full constructor
	
	virtual ~PML3DVISCOUS();                                                                //destructor
	const char* getClassType(void) const { return "PML3DVISCOUS"; };                        //return class type
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
	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);                         // update parameter
	int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

private:
	void calculateMatrices(); 								// calculate matrices method
	Domain* Domainptr;                              // pointer to the domain
	// double props[PML3DVISCOUS_NUM_PROPS];                  // material properties
	ID connectedExternalNodes;  					//eight node numbers
	Node* nodePointers[PML3DVISCOUS_NUM_NODES];    	    //pointers to eight nodes
	static double K[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // stiffness matrix
	static double C[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // damping matrix
	static double M[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // mass matrix
    static double G[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // G matrix
    static double H[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];        // H matrix
	static double Keff[PML3DVISCOUS_NUM_DOF * PML3DVISCOUS_NUM_DOF];     // effective stiffness matrix
	static double gamma; 					  	    // Newmark parameters: gamma
	static double beta; 					  	    // Newmark parameters: beta
	static double eta;                              // Newmark parameters: eta
	static double keisi;                            // Newmark parameters: keisi
	static Matrix tangent;                          // tangent matrix
	static Vector resid; 						    // residual vector
	static Matrix mass;						        // mass matrix
	static Matrix damping;	 					    // damping matrix
	static Matrix Gmat;                             // G matrix
	static Matrix Hmat;                             // H matrix
	static Matrix keffmat;                          // effective stiffness matrix
	Vector ubart; 				                    // ubar at time t 
    Vector ubarbart;                                // ubarbar at time t
	Vector ubar; 				                    // ubar at time t+dt
    Vector ubarbar;                                 // ubarbar at time t+dt
	static double dt; 								// time step
	int updateflag; 								// update flag
	int update_dt;                                  // flag for updating dt
	static int eleCount; 						    // element count
	static int ComputedEleTag; 					    // computed element tag
	double coords[24];
	// int innertag; 								// inner tag
	// static int numberOfElements; 			    // number of elements
	double alpha0;								    // alpha0 pml parameter
	double beta0;								    // beta0 pml parameter
	double xref, yref, zref;					    // reference point
	double nx, ny, nz;							    // normal vector
	const char* meshtype;							// type of  mesh
	double PML_L;                                    // PML minimum thickness
	double m_coeff;                                 // m parameter for PML
	double R_coeff;									    // PML damping ratio
	double cp_ref;								    // reference wave speed
	bool ExplicitAlphaBeta;						 // flag for explicit alpha beta
	
	NDMaterial* theMaterial;         // pointer to the material
	double E, nu, rho;          // material properties: E, nu, rho

};

#endif

