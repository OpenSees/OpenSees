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
// $Date: 2010-05-04 17:14:46 $
// $Source: /scratch/slocal/chroot/cvsroot/openseescomp/CompositePackages/MixedBeamColumnAsym3d/MixedBeamColumnAsym3d.h,v $

// Modified by: Xinlong Du and Jerome F. Hajjar, Northeastern University, USA; Year 2020
// Description: Adapted for analysis of asymmetric sections with introducing
// high-order axial terms for the basic element formulation
// References:
// Du, X., & Hajjar, J. F. (2021). Three-dimensional nonlinear mixed 6-DOF beam element 
// for thin-walled members. Thin-Walled Structures, 164, 107817. 

#ifndef MixedBeamColumnAsym3d_h
#define MixedBeamColumnAsym3d_h

// Written: Mark D. Denavit, University of Illinois at Urbana-Champaign
//
// Description: This file contains the interface for the MixedBeamColumnAsym3d class.
// It defines the class interface and the class attributes.
//
// What: "@(#) MixedBeamColumnAsym3d.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

#include <BeamIntegration.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>

class Node;
class Channel;
class Response;
class BeamIntegration;
class SectionForceDeformation;

class MixedBeamColumnAsym3d : public Element
{
  public:
    // constructors
    MixedBeamColumnAsym3d (int tag, int nodeI, int nodeJ,
            int numSections, SectionForceDeformation **sectionPtrs, BeamIntegration &bi,
            CrdTransf &coordTransf, double ys, double zs, double massDensPerUnitLength, int doRayleigh, bool geomLinear);
    MixedBeamColumnAsym3d ();

    // destructor
    ~MixedBeamColumnAsym3d();

    // public methods to obtain information about dof & connectivity
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
    const Matrix &getMass(void);
    const Matrix &getDamp(void);

    void zeroLoad(void);
    int addLoad(ElementalLoad *theLoad, double loadFactor);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);

    // public methods for output
    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag = 0);
    friend OPS_Stream &operator<<(OPS_Stream &s, MixedBeamColumnAsym3d &E);

    Response* setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInfo);

    const char *getClassType(void) const {return "MixedBeamColumnAsym3d";};

  protected:

  private:
    // Private Functions - Shape Functions
    Matrix getNld_hat(int sec, const Vector &v, double L, bool geomLinear);
    Vector getd_hat(int sec, const Vector &v, double L, bool geomLinear);
    Matrix getNd1(int sec, const Vector &v, double L, bool geomLinear);
    Matrix getNd2(int sec, double P, double L);
    Matrix getKg(int sec, Vector P, double L);
    Matrix getMd(int sec, Vector dShapeFcn, Vector dFibers, double L);

    // Private Functions - Interaction With The Sections
    //void getSectionTangent(int sec,int type,Matrix &kSection,double &GJ);
    //void getSectionStress(int sec,Vector &fSection,double &torsion);
    //void setSectionDeformation(int sec,Vector &defSection,double &twist);

    // Private Attributes - a copy for each object of the class
    ID connectedExternalNodes;              // tags of the end nodes
    Node *theNodes[2];                      // pointers to the nodes
    BeamIntegration *beamIntegr;            //
    int numSections;                        //
    SectionForceDeformation **sections;     // array of pointers to sections
    CrdTransf *crdTransf;                   // pointer to coordinate transformation object

    int doRayleigh;                         // flag for whether or not rayleigh damping is active for this element
    bool geomLinear;						            // flag for whether or not the internal geometric nonlinearity is active
    double rho;                             // mass density per unit length

    int itr;  // Counts the number of iterations from last committed state (not very sure)
    int initialFlag;

    // Attributes that do NOT change during the analysis
    double initialLength;
    Matrix *Ki;

    // Element Load Variables
    Matrix *sp;
    double p0[5]; // Reactions in the basic system due to element loads

    // Attributes that change during the analysis
    Vector V;
    Vector internalForce;          //Xinlong: internal resisting force (6)
    Vector naturalForce;           //Xinlong: generalized force degrees of freedom (7)
    Vector lastNaturalDisp;
    Matrix Hinv;
    Matrix GMH;
    Matrix kv;                     // stiffness matrix in the basic system
    Vector *sectionForceFibers;
    Vector *sectionDefFibers;
    Matrix *sectionFlexibility;
	Vector *sectionForceShapeFcn;

    // Committed versions
    Vector committedV;
    Vector committedInternalForce;
    Vector commitedNaturalForce;
    Vector commitedLastNaturalDisp;
    Matrix commitedHinv;
    Matrix commitedGMH;
    Matrix kvcommit;
    Vector *commitedSectionForceFibers;
    Vector *commitedSectionDefFibers;
    Matrix *commitedSectionFlexibility;

    // static data - single copy for all objects of the class
    //static int maxNumSections;
	enum { maxNumSections = 10 };
    static Matrix theMatrix;
    static Vector theVector;
    static double workArea[];
    //static Matrix transformNaturalCoords;
    //static Matrix transformNaturalCoordsT;
        // matrix to transform the natural coordinates from what the coordinate transformation uses and what the element uses

    // These variable are always recomputed, so there is no need to store them for each instance of the element
    static Vector *sectionDefShapeFcn;
    static Matrix *nldhat;
    static Matrix *nd1;
    static Matrix *nd2;
    static Matrix *nd1T;
    static Matrix *nd2T;

	double ys;     //Xinlong: y coord of shear center relative to centroid
	double zs;     //Xinlong: z coord of shear center relative to centroid
	
	enum { NDM = 3 };          // dimension of the problem (3d)
	enum { NND = 6 };          // number of nodal dof's
	enum { NGF = 7 };          //Xinlong: number of generalized force degrees of freedom
	enum { NSD = 5 };          //Xinlong: number of section dofs
	enum { NEGD = 12 };        // number of element global dof's
	enum { NEBD = 6 };         // number of element dof's in the basic system
};

#endif

