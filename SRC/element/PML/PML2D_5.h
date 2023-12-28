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


// Written by: Amin Pakzad, Pedro Arduino (parduino@uw.edu) and Adriano Trono
//
// Four node PML2D_5 element Derivation based on the Adriano Torino PhD thesis



#ifndef PML2D_5_H
#define PML2D_5_H

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <Domain.h>
#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>

#define PML2D_5_NUM_DOF   13
#define PML2D_5_NUM_NODES 5
#define PML2D_5_NUM_PROPS 8

class PML2D_5 : public Element {

public:
	PML2D_5();                                                                       //null constructor
	PML2D_5(int tag, int* nodeTags, double E, double nu,
	double rho, double pmlthicknessx, double pmlthicknessy, double Halfwidth, 
	double Depth, double r0, double R, double Vc);
	virtual ~PML2D_5();                                                              //destructor
	const char* getClassType(void) const { return "PML2D_5"; };                      //return class type
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
	void ComputeK(double* K,double* XYelement, double beta_0_x, double beta_0_y, double L_PML_x,
                                double L_PML_y, double xi, double yj, double rho, double E, double nu);
	void ComputeM(double* M,double* XYelement, double alpha_0_x, double alpha_0_y, double L_PML_x,
                                double L_PML_y, double xi, double yj, double rho, double E, double nu);
	void ComputeC(double* C,double* XYelement, double alpha_0_x, double alpha_0_y,
							  double beta_0_x, double beta_0_y, double L_PML_x,double L_PML_y, 
								double xi, double yj, double rho, double E, double nu);

private:
	double E;                
	double nu;				
	double rho;		     	
	double pmlthicknessx;    
	double pmlthicknessy;    
	double Halfwidth;	    
	double Depth;
	double r0;
	double R;
	double Vc;		    
	ID connectedExternalNodes;  					//five node numbers
	Node* nodePointers[PML2D_5_NUM_NODES];    	    //pointers to five nodes
	double K[PML2D_5_NUM_DOF * PML2D_5_NUM_DOF];        // stiffness matrix
	double C[PML2D_5_NUM_DOF * PML2D_5_NUM_DOF];        // damping matrix
	double M[PML2D_5_NUM_DOF * PML2D_5_NUM_DOF];        // mass matrix
	static Matrix tangent;                          // tangent matrix
	static Vector resid; 						    // residual vector
	static int eleCount;                            // element count

};

#endif

