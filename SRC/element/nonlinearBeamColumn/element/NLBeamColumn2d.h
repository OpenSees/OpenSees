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
                                                                        
// $Revision: 1.15 $
// $Date: 2003-02-27 17:15:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/nonlinearBeamColumn/element/NLBeamColumn2d.h,v $
                                                                        
                                                                        
// File: ~/model/element/NLBeamColumn2d.h
//
// Written by Remo Magalhaes de Souza on  09/98 
// Revised: rms 10/98 (uses section class)
//          rms 01/99 (with distributed loads)
//          rms 07/99 (using setDomain)
//          rms 08/99 (included P-Delta effect)
//	    fmk 10/99 setResponse() & getResponse()
//          rms 11/99 (included rigid joint offsets)
//          rms 04/00 (using transformation class w/ linear or corotational transf)
//          rms 04/00 (generalized to iterative/non-iterative algorithm)
//          mhs 06/00 (using new section class w/ variable dimensions)
//          rms 06/00 (torsional stiffness considered at the section level)
//          rms 06/00 (making copy of the sections)
//          rms 06/00 (storing section history variables at the element level)
//
//
// Purpose: This file contains the class definition for NLBeamColumn2d.
// NLBeamColumn2d is a materially nonlinear flexibility based frame element.

#ifndef NLBeamColumn2d_h
#define NLBeamColumn2d_h

#include <Element.h>
#include <Node.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <SectionForceDeformation.h>
#include <CrdTransf2d.h>
#include <GaussLobattoQuadRule1d01.h>

class Response;

class NLBeamColumn2d: public Element
{
  public:
    NLBeamColumn2d ();
    NLBeamColumn2d (int tag, int nodeI, int nodeJ, 
		    int numSections, SectionForceDeformation *sectionPtrs[], 
		    CrdTransf2d &coordTransf, double massDensPerUnitLength = 0.0, 
		    int maxNumIters = 10, double tolerance = 1e-12, int maxSub = 10);
    
    ~NLBeamColumn2d();

    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);

    int getNumDOF(void);
    void setDomain(Domain *theDomain);

    int commitState(void);
    int revertToLastCommit(void);        
    int revertToStart(void);
    int update(void);    
    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getMass(void);

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            
    
    bool isSubdomain(void);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &theViewer, int displayMode, float fact);        

    friend OPS_Stream &operator<<(OPS_Stream &s, NLBeamColumn2d &E);        
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, Information &eleInformation);
    int getResponse(int responseID, Information &eleInformation);
    
    int setParameter(const char **argv, int argc, Information &info);
    int updateParameter(int parameterID, Information &info);

  private:
    void getGlobalDispls(Vector &dg) const;
    void getGlobalAccels(Vector &ag) const;      
    void getForceInterpolatMatrix(double xi, Matrix &b, const ID &code);
    void getDistrLoadInterpolatMatrix(double xi, Matrix &bp, const ID &code);
    void compSectionDisplacements(Vector sectionCoords[], Vector sectionDispls[]) const;
    void initializeSectionHistoryVariables (void);
 

    // internal data
    ID     connectedExternalNodes; // tags of the end nodes
    int    nSections;              // number of sections (integration 
                                   // points) along the element
    SectionForceDeformation **sections;          // array of pointers to sections
    CrdTransf2d *crdTransf;        // pointer to coordinate tranformation object 
	                           // (performs the transformation between the global and basic system)
    double rho;                    // mass density per unit length
    int    maxIters;               // maximum number of local iterations
    double tol;	                   // tolerance for relative energy norm for local iterations

    // THESE SHOULD BE REMOVED!!! -- MHS
    double cosTheta, sinTheta;     // cossine directors

    int    initialFlag;            // indicates if the element has been initialized
    
    Node *theNodes[2];

    Vector load;                   // equivalent nodal loads ????

    Matrix kv;                     // stiffness matrix in the basic system 
    Vector Se;                     // element resisting forces in the basic system

    Matrix kvcommit;               // commited stiffness matrix in the basic system
    Vector Secommit;               // commited element end forces in the basic system

    Matrix *fs;                    // array of section flexibility matrices
    Vector *vs;                    // array of section deformation vectors
    Vector *Ssr;                   // array of section resisting force vectors

    Vector *vscommit;              // array of commited section deformation vectors

    Matrix *sp;  // Applied section forces due to element loads, 3 x nSections
    double p0[3]; // Reactions in the basic system due to element loads

    Matrix *Ki;
    int maxSubdivisions;

    
    static Matrix theMatrix;
    static Vector theVector;
    static GaussLobattoQuadRule1d01 quadRule;
    static double workArea[];

    static Vector *vsSubdivide; 
    static Matrix *fsSubdivide; 
    static Vector *SsrSubdivide;
    static int NLBeamColumn2d::maxNumSections;
};

#endif


