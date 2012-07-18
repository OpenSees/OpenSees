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
// Header file for TripleFrictionPendulum element
// Written by Nhan Dao, nhan.unr@gmail.com
#ifndef CircularEPPGap_h
#define CircularEPPGap_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>

class Node;
class Channel;

class TripleFrictionPendulum : public Element
{
  public:
    // constructors
    TripleFrictionPendulum(int tag,
    	int Nd1, int Nd2,
    	double L1,
    	double L2,
    	double L3,
    	double Ubar1,
    	double Ubar2,
    	double Ubar3,
    	double Mu1slow,
    	double Mu1fast,
    	double N1slow,
    	double N1fast,
    	double Alpha10,
    	double Alpha11,
    	double Alpha12,
    	double Mu2slow,
    	double Mu2fast,
    	double N2slow,
    	double N2fast,
    	double Alpha20,
    	double Alpha21,
    	double Alpha22,
    	double Mu3slow,
    	double Mu3fast,
    	double N3slow,
    	double N3fast,
    	double Alpha30,
    	double Alpha31,
    	double Alpha32,
    	double W,
    	double Uy,
    	double Kvc,
    	double Kvt,
		double minFv,
		double maxMuFac,
		double tol);

    TripleFrictionPendulum();    
    
    // destructor
    ~TripleFrictionPendulum();

    
    // public methods to obtain inforrmation about dof & connectivity
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
	Element *getCopy(void);
    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);    
    const Matrix &getMass(void);    

    const Vector &getResistingForce(void);

    // public methods for output    
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

  protected:
    
  private:
    // private member functions - only available to objects of the class

    void CircularElasticGap(Matrix &kj, Vector &fj, double Ej,double Gapj,Vector di);
    void BidirectionalPlastic(Matrix &ki, Vector &fi, Vector &epitmp, Vector &qitmp, double Fyi, double Ei, double Hi, Vector epi, Vector qi, Vector di);
    void Segment(Vector &epitmp, Vector &qitmp, bool &conv, Matrix &kij, Vector &di, Vector epi, Vector qi, Vector f, Vector df, double Fyi, double Ei, double Hi, double Ej, double Gapj, double Tol, int Niter);
    void TFPElement(bool &Conv, Vector &ep1tmp, Vector &ep3tmp, Vector &ep5tmp, Vector &q1tmp, Vector &q3tmp, Vector &q5tmp, Matrix &K, Vector &f, Matrix &k12, Matrix &k34, Matrix &k56, Vector &d1, Vector &d3, Vector &d5, Vector ep1, Vector ep3, Vector ep5, Vector q1, Vector q3, Vector q5, Vector u, Vector dusub, double Fy1, double Fy3, double Fy5, double E1, double E3, double E5, double H1, double H3, double H5, double E2, double E4, double E6, double Gap2, double Gap4, double Gap6, double Tol, int Niter);
    void StiffnessForm(Matrix &K, Matrix k12, Matrix k34, Matrix k56);
	void MInverse(Matrix MInput,Matrix &Mout,int ncol);
    double L1;
	double L2;
	double L3;
	double Ubar1;
	double Ubar2;
	double Ubar3;
	double Mu1slow;
	double Mu1fast;
	double N1slow;
	double N1fast;
	double Alpha10;
	double Alpha11;
	double Alpha12;
	double Mu2slow;
	double Mu2fast;
	double N2slow;
	double N2fast;
	double Alpha20;
	double Alpha21;
	double Alpha22;
	double Mu3slow;
	double Mu3fast;
	double N3slow;
	double N3fast;
	double Alpha30;
	double Alpha31;
	double Alpha32;
	double W;
	double Uy;
	double Kvc;
	double Kvt;
	double MinFv;
	double MaxMuFac;
	double TOL;
	int Niter;
	double A1slow;
	double A1fast;
	double av1;
	double avbar1;
	double A2slow;
	double A2fast;
	double av2;
	double avbar2;
	double A3slow;
	double A3fast;
	double av3;
	double avbar3;
    Matrix K;
    Matrix Kpr;
    Vector f;
    Vector fpr;
    Matrix k12;
    Matrix k12pr;
    Matrix k34;
    Matrix k34pr;
    Matrix k56;
    Matrix k56pr;
    Vector d1;
    Vector d1pr;
    Vector d3;
    Vector d3pr;
    Vector d5;
    Vector d5pr;
    Vector ep1;
    Vector ep1pr;
    Vector ep3;
    Vector ep3pr;
    Vector ep5;
    Vector ep5pr;
    Vector q1;
    Vector q1pr;
    Vector q3;
    Vector q3pr;
    Vector q5;
    Vector q5pr;
    Vector ep1tmp;
    Vector ep3tmp;
    Vector ep5tmp;
	Vector q1tmp;
	Vector q3tmp;
	Vector q5tmp;
	double Vel1Avg, Vel3Avg,Vel5Avg;
	double Wcr, Wavg, Wpr;
    double Fy1, Fy3, Fy5, Fy1pr, Fy3pr, Fy5pr;
    double E1, E3, E5;
    double H1, H3, H5;
    double E2, E4, E6;
    double Gap2, Gap4, Gap6;    
    bool Conv;
    double Fvert;
    double Kvert;
    double Hisolator;
    double Dx, Dy;

    // private attributes - a copy for each object of the class    
    ID  externalNodes;          	 // contains the id's of end nodes
    Matrix trans;       // hold the transformation matrix, could use a Vector
	                // if ever bother to change the Vector interface for
			// x-product.

    Node *theNodes[2];  // node pointers

    // static data - single copy for all objects of the class
    static Matrix eleK;   // class wide matrix for returning stiffness
	static Matrix eleKinit;   // class wide matrix for returning initial stiffness
    static Matrix eleD;   // class wide matrix for returning damping
    static Matrix eleM;   // class wide matrix for returning mass 
    static Vector eleR;   // class wide vector for returning residual
};
#endif

