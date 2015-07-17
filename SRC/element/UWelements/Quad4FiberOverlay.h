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

#ifndef Quad4FiberOverlay_h
#define Quad4FiberOverlay_h

// Created: M. Chiaramonte,  P. Arduino,  P.Mackenzie-Helnwein, UW, 03.29.2011
// Modified: Alborz Ghofrani, UW, 2011
//
// Description: This file contains the implementation of the Quad4FiberOverlay class


#include <Node.h>
#include <Element.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

// number of nodes per element
#define SL_NUM_NODE 4
// d.o.f. per node
#define SL_NUM_NDF  2
// degrees of freedom per element
#define SL_NUM_DOF  8
// displacement degrees of freedom per element
#define SL_NUM_DDOF  8

class Domain;
class Node;
class Channel;
class UniaxialMaterial;
class FEM_ObjectBroker;
class Response;

class Quad4FiberOverlay : public Element
{
	public:
		Quad4FiberOverlay(int tag, int nd1, int nd2, int nd3, int nd4,
				  UniaxialMaterial &m, double Af,double beta1, double beta2);
		Quad4FiberOverlay();
		~Quad4FiberOverlay();
		
		// initialization
		void setDomain(Domain *thedomain);
 
		// public methods to obtain information about dof & connectivity                 
		 int getNumExternalNodes(void) const;
		 const ID &getExternalNodes(void);
		 Node **getNodePtrs(void);
		 int getNumDOF(void);
 
		// public methods to set the state of the element                                 
		int commitState(void);
		int revertToLastCommit(void);
		int revertToStart(void);
		int update(void);
 
		// public methods to obtain stiffness  
		const Matrix &getTangentStiff(void);
		const Matrix &getInitialStiff(void);
 
		// public method to obtain resisting force
		const Vector &getResistingForce(void);
 
		// method for obtaining information specific to an element 
		void Print(OPS_Stream &s, int flag =0);
		Response *setResponse(const char **argv, int argc, OPS_Stream &s);
		int getResponse(int responseID, Information &eleInformation);
		int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

		// public methods for database/parallel processing                                                      
		int sendSelf(int commitTag, Channel &theChannel);
		int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
	
	protected:

	private:
		// private member functions - only available to objects of the class              
		
		double computeCurrentStrain(void) ;
		int Dual();
		int UpdateBase(double Xi, double Eta);
		int	getEltBb(double Xi, double Eta);

		// private attributes - a copy for each object of the class                       
		UniaxialMaterial *theMaterial;       // pointer to a material                     
		ID  externalNodes;                   // contains the id's of end nodes        
		Node *theNodes[4];

		// static data - single copy for all objects of the class  
		static Matrix FiberK;   // class wide matrix for returning stiffness                            
		static Vector P;		// class wide vector for returning residual 
		static double pts[2];	// Stores quadrature points
		static double wts;		// Stores quadrature weights

		Vector u;
		Vector g1;							// Basis
		Vector g2;							// Basis
		Vector dualg1;						// Dual Basis
		Vector dualg2;						// Dual Basis

		Vector Q1;							// Current crd 1
		Vector Q2;							// Current crd 2
		Vector Q3;							// Current crd 3
		Vector Q4;							// Current crd 4

		Vector Qfi;							// Fiber i Node
		Vector Qfj;							// Fiber j node
		Vector Vf;							// Fiber direction vector

		Vector nFi;
		Vector nFj;

		int iStartNode;						// Fiber starting interploation node
		int iEndNode;						// Fiber starting interploation node
		int jStartNode;						// Fiber ending interploation node
		int jEndNode;						// Fiber ending interploation node

		Matrix dNidxAlphai;					// J^-1 ^ dN/dxi	
		Vector A;							// Normalized fiber direction vector
		Vector AA;							// A tensor A
		Vector Bb;							// B mat 
		double beta1;						// Location of node for fiber
		double beta2;						// Location of node for fiber 
		double Af;							// Area of Fiber
		double Lf;							// Length of fiber
		double strain;						// strain
};

#endif


