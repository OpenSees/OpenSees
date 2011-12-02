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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/beamWithHinges/BeamWithHinges3d.h,v $
                                                                        
                                                                        
///////////////////////////////////////////////////////
// File:  ~/Src/element/beamWithHinges/BeamWithHinges3d.h
//
// Written by MHS
// 
// Written:  June 2000
//
// Purpose:  This header file contains the prototype
// for BeamWithHinges3d.  The element uses hinges (passed as
// object pointers) at the element ends.  The middle section
// is analyzed elastically, while the hinge ends return
// a section flexibility sampled at their middle points.  The element
// calculates a space frame (12x12) stiffness matrix from this info.
//

#ifndef BeamWithHinges3d_h
#define BeamWithHinges3d_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class Channel;
class FEM_ObjectBroker;

class SectionForceDeformation;

class CrdTransf3d;

class BeamWithHinges3d: public Element
{
  public:
    BeamWithHinges3d ();
    
    BeamWithHinges3d (int tag, int nodeI, int nodeJ,
		      double E, double Iz, double Iy, double A,
		      double G, double J,  double alpha,
		      SectionForceDeformation &sectionRefI, double hingeIlen, 
		      SectionForceDeformation &sectionRefJ, double hingeJlen,
		      CrdTransf3d &coordTrans, double shearL = 1.0,
		      double massDensPerUnitLength = 0.0, int max = 1, double tol = 1.0e-10);
    
    BeamWithHinges3d (int tag, int nodeI, int nodeJ,
		      double E, double Iz, double Iy, double A,
		      double G, double J, double alpha,
		      SectionForceDeformation &sectionRefI, double hingeIlen, 
		      SectionForceDeformation &sectionRefJ, double hingeJlen,
		      CrdTransf3d &coordTrans, const Vector &distLoad,
		      double shearL = 1.0, double massDensPerUnitLength = 0.0,
		      int max = 1, double tol = 1.0e-10);
    
    ~BeamWithHinges3d();
    
    int getNumExternalNodes (void) const;
    const ID &getExternalNodes (void);
    int getNumDOF (void);
    void setDomain (Domain *theDomain);

    int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);

    const Matrix &getTangentStiff (void);
    const Matrix &getSecantStiff (void);
    const Matrix &getDamp (void);
    const Matrix &getMass (void);

    void zeroLoad (void);
    int addLoad (const Vector &moreLoad);
	int addInertiaLoadToUnbalance(const Vector &accel);
    const Vector &getResistingForce (void);
    const Vector &getResistingForceIncInertia (void);

    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, 
		  FEM_ObjectBroker &theBroker);

    int setResponse (char **argv, int argc, Information &info);
    int getResponse (int responseID, Information &info);
    
    int setParameter (char **argv, int argc, Information &info);
    int updateParameter (int parameterID, Information &info);
    
    void Print (ostream &s, int flag = 0);
    
  protected:
    
  private:
    
    void setNodePtrs (Domain *theDomain);
    void getGlobalDispls (Vector &dg);
    void getGlobalAccels (Vector &ag);
    void setStiffMatrix (void);
    void setMass (void);
    
    CrdTransf3d *theCoordTransf;
    
    void getForceInterpMatrix (Matrix &b, double x, const ID &c, int &shearKeyVY, int &shearKeyVZ);
    void getDistrLoadInterpMatrix (Matrix &bp, double x, const ID &c);
    
    void setHinges (void);
    void setElasticFlex (void);
    
    ////////////////
    // Internal Data
	    
    double E, Iz, Iy, A, G, J, alpha, L;
    double massDens;
    double hingeIlen, hingeJlen;
		    
    ID connectedExternalNodes;    
    Node *node1Ptr, *node2Ptr; 
    SectionForceDeformation *sectionI, *sectionJ;
    Matrix K;
    Matrix m;
    Matrix d;

    Matrix fElastic;
    Matrix vElastic;
		
    Matrix b1, bp1;
    Matrix b3, bp3;

    Matrix fs1, fs3;
    
	Vector sr1, sr3;

    Vector e1, e3;
    
    Vector UePrev;
    Vector P;
    Vector Pinert;
    Matrix kb;
    Vector q;
    Vector load;
    Vector prevDistrLoad;
		
    Vector distrLoadCommit;
    Vector UeCommit;
    Matrix kbCommit;
    Vector qCommit;
    
    bool initialFlag;
    
    double shearLength;
    int shearIkeyVY, shearJkeyVY;
    int shearIkeyVZ, shearJkeyVZ;    

    double shearWeightIVY, shearWeightIVZ;
    double shearWeightJVY, shearWeightJVZ;
    
    int maxIter;
    double tolerance;
};

#endif
