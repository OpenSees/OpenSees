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
                                                                        
// $Revision: 6049 $
// $Date: 2015-07-17 01:56:36 -0300 (Fri, 17 Jul 2015) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/CatenaryCable/CatenaryCable.h $
                                                                        
                                                                        
#ifndef CatenaryCable_h
#define CatenaryCable_h

// Written: jaabell  (Jose Abell)
// Created: May 2017
// Revision: A
//
// Description: This element is a catenary cable, suitable for static and dynamic analysis of 
//              cable structures including thermal effects. Based on:
//
//  Salehi Ahmad Abad, M., Shooshtari, A., Esmaeili, V., & Naghavi Riabi, A. (2013). 
//        Nonlinear analysis of cable structures under general loadings. Finite Elements in Analysis and Design, 
//        73, 11–19. https://doi.org/10.1016/j.finel.2013.05.002
//
//  With dynamical extensions (mass matrix).
// 
//  Verification suite can be found in www.joseabell.com
//
// What: "@(#) CatenaryCable.h, revA"

// #define USE_QUADMATH 0

#ifdef USE_QUADMATH
#define FLOATTYPE __float128
#define SQRT sqrtq
#define LOG logq
#else
#define FLOATTYPE double
#define SQRT sqrt
#define LOG log
#endif

#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class CatenaryCable : public Element
{
  public:
    CatenaryCable(int tag, int node1, int node2, double weight, double E, double A, double L0, double alpha, double temperature_change, double rho, double error_tol, int Nsubsteps, int massType);
    
    CatenaryCable();    
    ~CatenaryCable();

    const char *getClassType(void) const {return "CatenaryCable";};

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
    const Matrix &getKi(void);
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);
    const Matrix &getDamp(void);    
    const Matrix &getMass(void);    

    void zeroLoad(void);	
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);            

    // public methods for element output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
    void Print(OPS_Stream &s, int flag =0);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    Vector getEnergyVector();

  protected:
    
  private:
    void compute_lambda0(void) ;
    void compute_projected_lengths(void) ;
    void compute_flexibility_matrix(void) ;
    void computeMass();
    void computeMassLumped();
    void computeMassByIntegration();
    void computeMassCloughStyle();
    void computeMassEquivalentTruss();
    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;     // contains the tags of the end nodes

    Vector *theLoad;    // pointer to the load vector P
    Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
    Vector *theVector;  // pointer to objects vector (a class wide Vector)


    double weight;          // weight of the CatenaryCable
    double E;               // Young's modulus of the CatenaryCable material
    double A;               // area of CatenaryCable
    double L0;              // unstreched, unexpanded length of CatenaryCable 
    double alpha;           // Coefficient of thermal expansion for CatenaryCable
    double temperature_change; // Change in temperature of the CatenaryCable
    double rho;             // rho: mass density per unit length of the CatenaryCable
    double error_tol;
    int Nsubsteps;

    double lx0, ly0, lz0;
    double w1, w2, w3;
    double f1, f2, f3;
    double l1, l2, l3;

    double lambda0;
    double l[3];            //Projected lengths

    double KE, PE;
    double KE_n, PE_n;

    bool first_step;

    int massType;

    Node *theNodes[2];
    Vector *load;
    Vector *load_incl_inertia;
    Vector *load_lastcommit;
    
    // static data - single copy for all objects of the class   
    static Matrix Flexibility;       // class wide flexibility matrix for iterations
    static Matrix Stiffness;
    static Matrix Mass;
    static Matrix ZeroMatrix;
    static Vector Forces;
};

#endif
