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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        
#ifndef TFP_Bearing2d_h
#define TFP_Bearing2d_h

// Written: Tracey Becker 
// Conversion from matlab to c++: fmk
//
// Description: This file contains the interface for the TFP_Bearing2d class.
//
// What: "@(#) TFP_Bearing2d.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class UniaxialMaterial;

class TFP_Bearing2d : public Element
{
  public:
    // constructors
  TFP_Bearing2d(int tag, 
		int Nd1, int Nd2, 
		double *r, 
		double *dio,
		double *di,
		double *mu,
		double *h,
		double H0,
		double a,
		double K,
		double vYield);


    TFP_Bearing2d();    
    
    // destructor
    ~TFP_Bearing2d();

    
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

    // public methods to obtain stiffness
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);

    // public method to obtain resisting force
    const Vector &getResistingForce(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    
    int displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    int kt3Drma(double *v, double *vp, double *Fr, double A, double *P, double *vpi);

  private:
    // parameters
    double r[4];
    double dOut[4];
    double dIn[4];
    double mu[4];
    double h[4];
    double K;
    double vyield;

    // state variables
    double d[8];
    double vs[8];

    double vpCommit[8];
    double vpTrial[8];
    double vCommit[8];
    double vTrial[8];
    double FrCommit[8];
    double FrTrial[8];

    double PCommit[4];
    double PTrial[4];
    double UCommit[4];
    double UTrial[4];

    double HTrial;
    double HCommit;
    double H0;
    double Ac;
    double Ap;
    double dh;
    double N[4];

    ID  externalNodes;  // contains the id's of end nodes
    Node *theNodes[2];  // node pointers

    int numDOF;
    Matrix *theMatrix;
    Vector *theVector;

};
#endif

