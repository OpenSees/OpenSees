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
                                                                        
// Written: Alborz Ghofrani, Pedro Arduino, U.Washington
// Created: Oct 2014
// Description: This file contains the class definition for QuadBeamEmbedContact.

#ifndef QuadBeamEmbedContact_h
#define QuadBeamEmbedContact_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// number of nodes per element
#define QBEC_NUM_NODE 6
// number of dimensions
#define QBEC_NUM_DIM  2
// degrees of freedom per element
#define QBEC_NUM_DOF  14

class Node;
class NDMaterial;
class Response;

class QuadBeamEmbedContact : public Element
{
  public:
    QuadBeamEmbedContact(int tag, int Qnd1, int Qnd2, int Qnd3, int Qnd4,
		 int Bnd1, int Bnd2, double r, double mu, double ep, double et);
    QuadBeamEmbedContact();
    ~QuadBeamEmbedContact();

    const char *getClassType(void) const {return "QuadBeamEmbedContact";};

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    // public methods to set the state of the element    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    int update(void);

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    

    const Vector &getResistingForce(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
	
  protected:
    
  private:
    // private attributes - a copy for each object of the class

    ID externalNodes; // Tags of quad and beam nodes

    Node *theNodes[6];
	
	static Vector		m_ContactForces;	// force vector
	static Matrix		m_ContactStiffness;	// stiffness matrix
	static const double m_Pi;

	double	m_et;		// penalty parameter
	double	m_ep;		// penalty parameter
	double	m_radius;	// beam Radius
	double	m_friction;	// friction coefficient
	
	Vector	m_Q1;		// quad node 1 coordinates
	Vector	m_Q2;		// quad node 2 coordinates
	Vector	m_Q3;		// quad node 3 coordinates
	Vector	m_Q4;		// quad node 4 coordinates
	Vector	m_B1;		// beam node 1 coordinates
	Vector	m_B2;		// beam node 2 coordinates
	Vector	m_Ba;		// beam tangent at node 1
	Vector	m_Bb;		// beam tangent at node 2

	Vector	m_Q1_c;		// current quad node 1 coordinates
	Vector	m_Q2_c;		// current quad node 2 coordinates
	Vector	m_Q3_c;		// current quad node 3 coordinates
	Vector	m_Q4_c;		// current quad node 4 coordinates
	Vector	m_B1_c;		// current beam node 1 coordinates
	Vector	m_B2_c;		// current beam node 2 coordinates
	Vector	m_Ba_c;		// current beam tangent at node 1
	Vector	m_Bb_c;		// current beam tangent at node 2

	Vector	m_Q1_disp;	// quad node 1 displacements
	Vector	m_Q2_disp;	// quad node 2 displacements
	Vector	m_Q3_disp;	// quad node 3 displacements
	Vector	m_Q4_disp;	// quad node 4 displacements
	Vector	m_B1_disp;	// beam node 1 displacements
	Vector	m_B2_disp;	// beam node 2 displacements

	double	m_eta_n;	// beam contact point location parameter at time n
	Vector	m_y_n;		// beam contact point location vector at time n
	Vector	m_y_n1;		// beam contact point location vector at time n+1
	Vector	m_xi_n;		// solid contact point location parameter at time n
	Vector	m_x_n;		// solid contact point location vector at time n
	Vector	m_x_n1;		// solid contact point location vector at time n+1
	Vector	m_yc0;		// beam contact point at time 0
	
	double	m_Length;	// beam length
	double	m_Gap;		// gap and gap function
	double	m_Slip;		// slip
	double	m_Phi;		// slip function
	double	m_SigmaN;	// contact pressure
	double	m_SigmaT_n;	// friction pressure step n
	double	m_SigmaT;	// friction pressure step n+1
	double  m_SignSigmaT;	// friction pressure sign

	Vector	m_bShape;
	Vector	m_sShape;
	Vector	m_DbShape;

	Vector	m_t1;		// covariant tangent vector to beam centerline
	Vector	m_n;		// normal vector to beam centerline
	Vector	m_tau1;		// normalized t1
	Vector	m_nu;		// normalized n
	
	Vector	m_Bn, m_Bs;	// Bn and Bs
	Vector	m_beamIntersect;
	Vector	m_nFi, m_nFj;
	bool	createElement;
	bool	inContact;
	bool	normDepOnY;
	bool	isStuck, wasStuck;

	int		intersection(Vector& solidIntersect, Vector& beamIntersect);
	int		getContactPt(Vector& nFi, Vector& nFj, Vector& solidContactPt, double& beamContactPt);
	int		updateShapeFuncs(Vector xi, double eta);
	int		updateBase(double eta);
	int		project(Vector& xi, Vector& x_n1);
	void	computeB(void);
	double	getIntJacobian(void);

};

#endif

