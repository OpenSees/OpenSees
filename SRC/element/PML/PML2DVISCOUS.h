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
// Eight node PML2DVISCOUS element .. a c++ wrapper to fortran routine 
// provided by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2023


#ifndef PML2DVISCOUS_H
#define PML2DVISCOUS_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <Domain.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>

#define PML2DVISCOUS_NUM_DOF 20
#define PML2DVISCOUS_NUM_PROPS 11
#define PML2DVISCOUS_NUM_NODES 4

#ifdef _WIN32

#define pml2d_	      PML_2D

extern "C" void  pml2d_(double* kMatrix,
	double* cMatrix,
	double* mMatrix,
	double* gMatrix,
	int* NDOFEL,
	double* PROPS,
	int* NPROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE);

#else

#define pml2d_	      pml_2d_

extern "C" void  pml2d_(double* kMatrix,
	double* cMatrix,
	double* mMatrix,
	double* gMatrix,
	int* NDOFEL,
	double* PROPS,
	int* NPROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE);

#endif

class PML2DVISCOUS : public Element {

public:

	PML2DVISCOUS();                                                                         //null constructor
	PML2DVISCOUS(int tag, int* nodeTags, double* newmarks, double* dData);                  // full constructor
	virtual ~PML2DVISCOUS();                                                                //destructor
	const char* getClassType(void) const { return "PML2DVISCOUS"; };                        //return class type
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

	Domain* Domainptr;                              // pointer to the domain
	double props[PML2DVISCOUS_NUM_PROPS];                  // material properties
	ID connectedExternalNodes;  					//eight node numbers
	Node* nodePointers[PML2DVISCOUS_NUM_NODES];    	    //pointers to eight nodes
	double K[PML2DVISCOUS_NUM_DOF * PML2DVISCOUS_NUM_DOF];        // stiffness matrix
	double C[PML2DVISCOUS_NUM_DOF * PML2DVISCOUS_NUM_DOF];        // damping matrix
	double M[PML2DVISCOUS_NUM_DOF * PML2DVISCOUS_NUM_DOF];        // mass matrix
    double G[PML2DVISCOUS_NUM_DOF * PML2DVISCOUS_NUM_DOF];        // G matrix
	double Keff[PML2DVISCOUS_NUM_DOF * PML2DVISCOUS_NUM_DOF];     // effective stiffness matrix
	static double eta;                              // Newmark parameters: eta
	static double beta; 					  	    // Newmark parameters: beta
	static double gamma; 					  	// Newmark parameters: gamma
	static Matrix tangent;                          // tangent matrix
	static Vector resid; 						    // residual vector
	static Matrix mass;						        // mass matrix
	static Matrix damping;	 					    // damping matrix
	Vector ubart; 				                    // ubar at time t 
	Vector ubar; 				                    // ubar at time t+dt
	static double dt; 								// time step
	int updateflag; 								// update flag
	int update_dt;                                  // flag for updating dt
	double alpharayleigh; 						    // rayleigh damping parameter
	double betarayleigh; 							// rayleigh damping parameter
	static int eleCount; 						    // element count
	// int innertag; 								// inner tag
	// static int numberOfElements; 			    // number of elements
};

#endif

