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

/* Written by: Mohammad Salehi (mohammad.salehi@tamu.edu)
** Created: 07/19
** Description: The source code for the 2D gradient inelastic (GI) force-based beam-column element formulation
**
**
** References:
** 
** Mohammad Salehi and Petros Sideris (2017)
** “Refined Gradient Inelastic Flexibility-Based Formulation for Members Subjected to Arbitrary Loading”
** ASCE Journal of Engineering Mechanics, 143(9): 04017090
**
** Petros Sideris and Mohammad Salehi (2016)
** “A Gradient Inelastic Flexibility-Based Frame Element Formulation”
** ASCE Journal of Engineering Mechanics, 142(7): 04016039
*/

#ifndef GradientInelasticBeamColumn2d_H
#define GradientInelasticBeamColumn2d_H

// Directives
#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <BeamIntegration.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>

// Forward Declaration
class Response;

class GradientInelasticBeamColumn2d : public Element
{
public:
	// Constructors
	GradientInelasticBeamColumn2d();
	GradientInelasticBeamColumn2d(int tag, int nodeI, int nodeJ,
		int numSec, SectionForceDeformation &endSec1, SectionForceDeformation &sec, SectionForceDeformation &endSec2, double R1, double R2,
		BeamIntegration &BI, CrdTransf &CT, double LC,
		double minTolerance = 1.0e-10, double maxTolerance = 1.0e-8, int maxNumIters = 50,
		bool constH = false,
		bool corControl = false, double maxEps = 0.0, double maxPhi = 0.0);

	// Destructor
	~GradientInelasticBeamColumn2d();

	// Method to Get Class Type
	const char *getClassType() const { return "GradientInelasticBeamColumn2d"; };

	// Method to Initialize the Domain; base class: DomainComponent
	void setDomain(Domain *theDomain);

	// Methods to Obtain Information about DOFs and Connectivity; base class: Element
	int getNumExternalNodes(void) const;
	const ID &getExternalNodes(void);
	Node **getNodePtrs(void);
	int getNumDOF(void);

	// Methods to Set the State of Element; base class: Element
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);
	int update(void);

	// Methods to Obtain Stiffness Matrices; base class: Element
	const Matrix &getTangentStiff(void);
	const Matrix &getInitialStiff(void);
	const Matrix &getMass(void);

	// Methods to Obtain Resisting Forces; base class: Element
	const Vector &getResistingForce(void);
	const Vector &getResistingForceIncInertia(void);

	// Methods to Obtain Information Specific to Element; base class: Element
	void Print(OPS_Stream &s, int flag = 0);
	Response *setResponse(const char **argv, int argc, OPS_Stream &output);
	int getResponse(int responseID, Information &eleInfo);

	// Method to Display Element
	int displaySelf(Renderer& theViewer, int displayMode, float fact, const char** displayModes = 0, int numModes = 0);

	// Methods to Do Parallel Processing; base class: Channel
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

private:
	// Private Attributes
	ID connectedExternalNodes;              // contains tags of end nodes
	Node *theNodes[2];                      // pointer to nodes
	SectionForceDeformation **sections;     // pointers to sections
	BeamIntegration *beamIntegr;			// pointer to integration method
	CrdTransf *crdTransf;					// pointer to coordinate transformation method

	int numSections;		// number of sections
	int secOrder;			// section's order of behavior
	int maxIters;           // maximum number of element-level iterations

	double secLR1;			// length ratio for end section 1
	double secLR2;			// length ratio for end section 2

	bool correctionControl;	// indicates whether to update control size of iteration correction or not
	bool cnstH;				// indicates whether to keep [H] constant or not

	double lc;				// characteristic length
	double minTol;          // minimum acceptable tolerance for element-level iterations
	double maxTol;          // maximum acceptable tolerance for element-level iterations
	double F_tol_q;			// tolerance multiplier for q convergence test
	double F_tol_f_ms;		// tolerance multiplier for F_ms convergence test
	double maxEpsInc;
	double maxPhiInc;
	double L;				// element's initial length

	Matrix *B_q;			// total displacement integration matrix
	Matrix *B_Q;			// total force interpolation matrix
	Matrix *H;				// nonlocal averaging matrix
	Matrix *H_init;			// initial nonlocal averaging matrix
	Matrix *H_inv;			// inverse of nonlocal averaging matrix
	Matrix *B_q_H_inv_init;	// initial B_q/H
	Matrix *K0;				// pointer to initial stiffness matrix in the basic system

	// State Variables
	int initialFlag;        // flag which indicates if element has been initialized
	int iterNo;				// to record iteration numbers needed to converge
	int strIterNo;
	int totStrIterNo;
	int commitNo;			// total of times commitState has been called
	Vector iters;

	Matrix *J;				// Jacobian matrix
	Matrix *J_init;
	Matrix *J_commit;

	Vector Q;               // element resisting forces in the basic system
	Vector Q_commit;        // committed element resisting forces in the basic system
	Vector *d_tot;			// all section deformations
	Vector *d_tot_commit;	// committed section deformations
	Vector *d_nl_tot;		// all nonlocal strains
	Vector *d_nl_tot_commit;// committed nonlocal strains
	Vector *F_ms;			// section forces
	Vector *F_ms_commit;	// committed section forces
	Vector k_init;			// element's initial basic stiffness matrix's diagonal
	Vector *flex_ms_init;	// sections' flexibility matrix's diagonal
	Vector *trial_change;
	Vector *max_trial_change;

	Vector *hh;

	Vector *d_sec;          // array of section deformations
	Vector *d_sec_commit;   // array of committed section forces

	// complete the state variables

	// Private Methods
  void setSectionPointers(void);
  
	void assembleMatrix(Matrix &A, const Matrix &B, int rowStart, int rowEnd, int colStart, int colEnd, double fact);
	void assembleMatrix(Matrix &A, const Vector &B, int col, double fact);
	void assembleVector(Vector &A, const Vector &B, int rowStart, int rowEnd, double fact);

	void getSectionsTangentStiff(Matrix &tStiff);
	void getSectionsInitialStiff(Matrix &iStiff);

	const Matrix &getBasicStiff(void);
	const Matrix &getInitialBasicStiff(void);

	double weightedNorm(const Vector &W, const Vector &V, bool sqRt = true);

	bool qConvergence(const int &iter, const Vector &qt, const Vector &dnl_tot, Vector &Dq, double &dqNorm);
	bool fConvergence(const int &iter, const Vector &Qt, Vector &DF_ms, double &dfNorm);

	// Static Class Wide Variables
	static Matrix theMatrix;
	static Vector theVector;
};

#endif // GradientInelasticBeamColumn2d_H
