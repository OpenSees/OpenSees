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


// Written by: Long Chen, Pedro Arduino (parduino@uw.edu), Wenyang Zhang and fmk
//
// Eight node PML3D element .. a c++ wrapper to fortran routine 
// providewd by Wenyang Zhang (zwyll@ucla.edu), University of California, Los Angeles
//
// University of Washington, UC. Los Angeles, U.C. Berkeley, 12, 2020


#ifndef PML3D_H
#define PML3D_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>

#define PML3D_NUM_DOF 144
#define PML3D_NUM_PROPS 12
#define PML3D_NUM_NODES 8

#ifdef _WIN32

#define pml_	      PML_3D

extern "C" void  pml_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	int* NDOFEL,
	double* PROPS,
	int* NPROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE);

#else

#define pml_	      pml_3d_

extern "C" void  pml_(double* mMatrix,
	double* cMatrix,
	double* kMatrix,
	int* NDOFEL,
	double* PROPS,
	int* NPROPS,
	double* COORDS,
	int* MCRD,
	int* NNODE);

#endif

class PML3D : public Element {

public:

	//null constructor
	PML3D();

	//full constructor
	PML3D(int tag,
		int* nodeTags,
		double* dData);

	//destructor 
	virtual ~PML3D();

	const char* getClassType(void) const { return "PML3D"; };

	//set domain
	void setDomain(Domain* theDomain);

	//get the number of external nodes
	int getNumExternalNodes() const;

	//return connected external nodes
	const ID& getExternalNodes();
	Node** getNodePtrs(void);

	//return number of dofs
	int getNumDOF();

	//commit state
	int commitState();

	//revert to last commit 
	int revertToLastCommit();

	//revert to start 
	int revertToStart();

	// update
	int update(void);

	//print out element data
	void Print(OPS_Stream& s, int flag);

	//return stiffness matrix 
	const Matrix& getTangentStiff();
	const Matrix& getInitialStiff();
	const Matrix& getMass();
	const Matrix& getDamp();

	void zeroLoad();
	int addLoad(ElementalLoad* theLoad, double loadFactor);
	// int addInertiaLoadToUnbalance(const Vector &accel);

	//get residual
	const Vector& getResistingForce();

	//get residual with inertia terms
	const Vector& getResistingForceIncInertia();

	// public methods for element output
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker
		& theBroker);

	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& eleInformation);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int parameterID, Information& info);

	//plotting 
	int displaySelf(Renderer&, int mode, float fact, const char** displayModes = 0, int numModes = 0);

private:

	//
	// private attributes
	//

	double props[PML3D_NUM_PROPS];
	ID connectedExternalNodes;  //eight node numbers
	Node* nodePointers[PML3D_NUM_NODES];      //pointers to eight nodes

	double K[PML3D_NUM_DOF * PML3D_NUM_DOF];
	double C[PML3D_NUM_DOF * PML3D_NUM_DOF];
	double M[PML3D_NUM_DOF * PML3D_NUM_DOF];

	//
	// static attributes
	//

	static Matrix tangent;
	static Vector resid;
	static Matrix mass;
	static Matrix damping;
};

#endif

