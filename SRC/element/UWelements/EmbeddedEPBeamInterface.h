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

// Written: Alborz Ghofrani, Diego Turello, Pedro Arduino, U.Washington 
// Created: May 2017
// Description: This file contains the class definition for EmbeddedEPBeamInterface.

#ifndef EmbedEPBeamInterface_h
#define EmbedEPBeamInterface_h

#include <Element.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

#include <vector>
#include <set>
#include <map>

// number of dimensions
#define EEPBIP_NUM_DIM  3

class Node;
class NDMaterial;
class Response;
class CrdTransf;
class Parameter;

class EmbeddedEPBeamInterface : public Element
{
public:
    EmbeddedEPBeamInterface(int tag);
    EmbeddedEPBeamInterface(int tag, std::vector <int> beamTag, std::vector <int> solidTag, int crdTransfTag, int matTag, 
    std::vector <double>  beamRho, std::vector <double>  beamTheta, std::vector <double>  solidXi, std::vector <double>  solidEta,
    std::vector <double>  solidZeta, double radius, std::vector <double> area, std::vector <double> length, 
    bool writeConnectivity = false, const char * connectivityFN = "", double width = 0.0);
    EmbeddedEPBeamInterface();
    ~EmbeddedEPBeamInterface();

    const char *getClassType(void) const { return "EmbeddedEPBeamInterface"; };

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
    void Print(OPS_Stream &s, int flag = 0);

    Response *setResponse(const char **argv, int argc,
        OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

protected:

private:

    int EEPBIP_numNodes, EEPBIP_numDOF, m_numSolidDOF;

    ID externalNodes; // Tags of beam and solid nodes

    int *theSolidTags;
    int *solidNodeTags;
    int *theBeamTags;
    int *beamNodeTags;
    int theMatTag;

    Node **theNodes;

    Vector		m_InterfaceForces;	// force vector
    Matrix		m_InterfaceStiffness;	// stiffness matrix

    std::map<int, int> m_solidNodeMap;
    std::map<int, int> m_beamNodeMap;

    // iso-parametric coordinates
    Vector m_solid_xi;
    Vector m_solid_eta;
    Vector m_solid_zeta;
    Vector m_beam_rho;
    Vector m_beam_theta;
    Vector m_area;
    Vector m_beamLength;

    // shape functions
    Vector  m_Ns;
    double  m_Hb1, m_Hb2, m_Hb3, m_Hb4, m_Nb1, m_Nb2;
    double  m_dH1, m_dH2, m_dH3, m_dH4;

    double	m_beam_radius;	 // beam Radius
    double  m_beam_length;   // beam length
    double  m_intWidth;      // interface width

    int     m_numBeamNodes, m_numSolidNodes, m_numEmbeddedPoints;

    CrdTransf*  crdTransf;  // pointer to coordinate transformation object
    NDMaterial** theMat;     // pointer to the interface material

    Vector  m_Lambda;

    // initial displacement
    Vector m_solidInitDisp;
    Vector m_beamInitDisp;

    Vector m_uRel;


    // copied from BeamContact3D
    Matrix mQa;                 // coordinate transform for node a
    Matrix mQb;                 // coordinate transform for node b
    Matrix mQc;
    Vector mc1;                 // tangent vector at project point c
    Matrix mBphi, mBu, mHf;
    Matrix mA, mB, mAt, mBt, mAAt, mBBt, mABt, mBinv, mBinvTAT;

    void ComputeBphiAndBu(Matrix &Bphi, Matrix &Bu, double L);            // method to compute Bphi and Bu, used in ComputeB and update
    void ComputeHf(Matrix &Hf, double theta);                   // method to compute Hf
    int	 updateShapeFuncs(double xi, double eta, double zeta, double rho, double L);             // method to update shape functions

    Matrix ComputeSkew(Vector theta);    // function returns skew matrix of given vector
    Matrix Transpose(int dim1, int dim2, const Matrix &M);   // functions returns the tranpose of Matrix M (does not exist in Matrix Class!)

   
    // define some functions for recording results purposes
    Vector GetInteractionPtDisp(int flag);
    Vector GetInteractionPtForce(int flag);
    Vector GetElementalForce(int flag);

};

extern NDMaterial *OPS_getNDMaterial(int tag);

#endif

