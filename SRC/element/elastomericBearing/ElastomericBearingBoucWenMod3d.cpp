/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** ****************************************************************** */

// Written: Ping Tan ��Guangzhou University, China��
// Created: 2015-10-12 18:46:30
// Revision: A
//
// Description: This file contains the implementation of the
// ElastomericBearingBoucWenMod class.

// What: "@(#) ElastomericBearingBoucWenMod.cpp, revA"

#include "ElastomericBearingBoucWenMod3d.h"
#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <ElementResponse.h>
#include <UniaxialMaterial.h>

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <G3Globals.h>
#include <Message.h>
using namespace std;
#include <iostream>

#define PI 3.14159l

// initialize the class wide variables
Matrix ElastomericBearingBoucWenMod3d::theMatrix(12,12);
Vector ElastomericBearingBoucWenMod3d::theVector(12);


#include <elementAPI.h>
#define OPS_Export 

void *
OPS_ElastomericBearingBoucWenMod3d(void)
{
  Element *theElement = 0;
  
  int numRemainingArgs = OPS_GetNumRemainingInputArgs();
  
  if (numRemainingArgs < 18) {
    opserr << "Want: ElastomericBearingBoucWenMod eleTag iNode jNode kInit fy Gr Kbulk D1 D2 ts tr n alpha1 alpha2 mu eta beta gamma <-PMod a1 a2 > <-TMod T b1 b2 b3 b4> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol> <-orient <x1 x2 x3> y1 y2 y3>\n";
    return 0;	
  }
  int tag;
  int iNode, jNode;
  double kInit, fy, alpha1, Gr, Kbulk, D1, D2, ts, tr;
  int n;

  double alpha2 = 0.0;
  double mu = 2.0;
  double eta = 1.0;
  double beta = 0.5;
  double gamma = 0.5;
  double a1 = 0.0;
  double a2 = 1.0;
  double T = 23.0;
  double b1 = 1.0;
  double b2 = 0.0;
  double b3 = 0.0;
  double b4 = 0.0;
  double shearDistI = 0.5;
  
  int doRayleigh = 0;
  double mass = 0.0;
  int maxIter = 25;
  double tol = 1E-12;  
  
  Vector x(0);
  Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
  
  int iData[3];
  double dData[15];
  
  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d error reading integer data\n";
    return 0;
  }
  tag = iData[0];
  iNode = iData[1];
  jNode = iData[2];
  
  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading data\n";
    return 0;	
  }
  kInit = dData[0];
  fy = dData[1];
  Gr= dData[2];
  Kbulk = dData[3];
  D1 = dData[4];
  D2 = dData[5];
  ts = dData[6];
  tr = dData[7];

  numData = 1;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading data\n";
    return 0;	
  }
  n = iData[0];

  numData = 6;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading data\n";
    return 0;	
  }
  alpha1 = dData[0];
  alpha2 = dData[1];
  mu = dData[2];
  eta = dData[3];
  beta = dData[4];
  gamma = dData[5];

  // while not done, loop over options
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *option = OPS_GetString();
    if (strcmp(option,"-PMod") == 0)  {
      numData = 2;
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -PMod data\n";
	break;
      }
      a1 = dData[0];
      a2 = dData[1];
    } else if (strcmp(option,"-TMod") == 0)  {
      numData = 4;
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -TMod data\n";
	break;
      }
      b1 = dData[0];
      b1 = dData[1];
      b1 = dData[2];
      b1 = dData[3];
    } else if (strcmp(option,"-shearDist") == 0)  {
      numData = 1;
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -shearDist data\n";
	break;
      }
      shearDistI = dData[0];
    } else if (strcmp(option,"-doRayleigh") == 0)  {
      doRayleigh  = 1;
    } else if (strcmp(option,"-mass") == 0)  {
      numData = 1;
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -mass data\n";
	break;
      }
      mass = dData[0];
    } else if (strcmp(option,"-iter") == 0)  {
      numData = 1;
      if (OPS_GetInt(&numData, iData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -iter data\n";
	break;
      }
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -iter data\n";
	break;
      }
      maxIter = iData[0];
      tol = dData[0];
    } else if (strcmp(option,"-orient") == 0)  {
      int numOrient = OPS_GetNumRemainingInputArgs();
      numData = numOrient;
      if (numData != 3 || numData != 6) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -orient data, need 3 or 6 values\n";
      }
      if (OPS_GetDouble(&numData, dData) != 0) {
	opserr << "WARNING element ElastomericBeamringMod3d tag: " << tag << " error reading -orient data\n";
	break;
      }
      if (numOrient == 3) {
	y(0) = dData[0];
	y(1) = dData[1];
	y(2) = dData[2];
      } else if (numOrient == 6) {
	x.resize(3);
	x(0) = dData[0];
	x(1) = dData[1];
	x(2) = dData[2];
	y(0) = dData[3];
	y(1) = dData[4];
	y(2) = dData[5];
      }
    }
  }
  theElement = new ElastomericBearingBoucWenMod3d(tag, iNode, jNode, kInit, fy, Gr, 
						  Kbulk, D1, D2, ts, tr, n, alpha1, 
						  alpha2,  mu, eta, beta, gamma,
						  a1, a2, T, b1, b2, b3, b4, y, x, 
						  shearDistI, doRayleigh, mass, 
						  maxIter, tol);
  
  return theElement;
}

ElastomericBearingBoucWenMod3d::ElastomericBearingBoucWenMod3d(int tag,
							       int Nd1, int Nd2, double _kInit, double _fy,
							       double Gr, double Kbulk, double Di, double Do,
							       double ts, double tr, int n,double _alpha1,
							       double _alpha2, double _mu, double _eta, double _beta, double _gamma,
							       double _a1, double _a2, double _T,
							       double _b1, double _b2, double _b3, double _b4,
							       const Vector _y, const Vector _x,
							       double sdI, int addRay, double m, int maxiter, double _tol)
  : Element(tag, ELE_TAG_ElastomericBearingBoucWenMod3d),
    connectedExternalNodes(2), kInit(_kInit),fy(_fy) ,G(Gr), k0(0.0), qYield(0.0), k2(0.0), k3(0.0),
    D1(Di), D2(Do), alpha1(_alpha1), alpha2(_alpha2), mu(_mu), eta(_eta), beta(_beta), gamma(_gamma), 
    a1(_a1), a2(_a2), T(_T), b1(_b1), b2(_b2), b3(_b3), b4(_b4),A(1.0), x(_x), y(_y),
    shearDistI(sdI), addRayleigh(addRay), mass(m), maxIter(maxiter), tol(_tol),
    L(0.0), onP0(true), ub(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6), theLoad(12)
{
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "ElastomericBearingBoucWenMod3d::ElastomericBearingBoucWenMod3d() - element: "
	   << this->getTag() << " - failed to create an ID of size 2.\n";
    exit(-1);
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
  
  // initialize parameters
  k0 = (1.0-alpha1)*kInit;
  qYield = (1.0-alpha1-alpha2*pow(fy/kInit,mu-1.0))*fy;
  k2 = alpha1*kInit;
  k3 = alpha2*kInit;
  double Ac=(PI/4)*(D2*D2-D1*D1);
  S=(D2-D1)/(4*tr);
  Tr=n*tr;
  h=Tr + (n-1)*ts;
  double F;
  if(D1<DBL_EPSILON) {
    F=1.0;
  } else {
    double r=D2/D1;										// Outer to inner diameter ratio
    F=(r*r+1)/((r-1)*(r-1)) + (1+r)/((1-r)*log(r));		// Dimension modification factor
  }
  Ec=1.0/((1/(6*G*S*S*F))+(4.0/3.0)*(1/Kbulk));			// Compressive modulus of elastomeric bearing
  
  double I=(PI/64)*(pow(D2,4)-pow(D1,4));            // Moment of inertia of bearing 
  Kv0=Ac*Ec/Tr;											// Pre-cavitation stiffness at zero lateral displacement 
  double Er=Ec/3; // Rotation modulus of bearing
  double Is=I*h/Tr;										// Adjusted moment of intertia of bearing
  // Rotation
  Kr= Er*Is/h;
  
  // Torsion
  Kt=G*(2*Is)/h;
   
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = kbInit(2,2) = A*k0 + k2;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
}


ElastomericBearingBoucWenMod3d::ElastomericBearingBoucWenMod3d()
    : Element(0, ELE_TAG_ElastomericBearingBoucWenMod3d),
    connectedExternalNodes(2),  kInit(0.0), fy(0.0), k0(0.0), qYield(0.0), k2(0.0), k3(0.0),Kv0(0.0),
    Tr(0.0), D1(0.0), D2(0.0), alpha1(0.0), alpha2(0.0), mu(2.0), eta(1.0), beta(0.5), gamma(0.5), 
	a1(0.0), a2(1.0), T(23.0), b1(1.0), b2(0.0), b3(0.0), b4(0.0), A(1.0), x(0), y(0),
    shearDistI(0.5), addRayleigh(0), mass(0.0), maxIter(25), tol(1E-12),
    L(0.0), onP0(false), ub(6), z(2), dzdu(2,2), qb(6), kb(6,6), ul(12),
    Tgl(12,12), Tlb(6,12), ubC(6), zC(2), kbInit(6,6), theLoad(12)
{
  // ensure the connectedExternalNode ID is of correct size & set values
  if (connectedExternalNodes.Size() != 2)  {
    opserr << "ElastomericBearingBoucWenMod3d::ElastomericBearingBoucWenMod3d() - element: "
	   << this->getTag() << " - failed to create an ID of size 2.\n";
    exit(-1);
  }
  
  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;
}


ElastomericBearingBoucWenMod3d::~ElastomericBearingBoucWenMod3d()
{
  // invoke the destructor on any objects created by the object
  // that the object still holds a pointer to
}


int ElastomericBearingBoucWenMod3d::getNumExternalNodes() const
{
    return 2;
}


const ID& ElastomericBearingBoucWenMod3d::getExternalNodes()
{
    return connectedExternalNodes;
}


Node** ElastomericBearingBoucWenMod3d::getNodePtrs()
{
    return theNodes;
}


int ElastomericBearingBoucWenMod3d::getNumDOF()
{
    return 12;
}


void ElastomericBearingBoucWenMod3d::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (!theDomain)  {
        theNodes[0] = 0;
        theNodes[1] = 0;
        
        return;
    }
    
    // first set the node pointers
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));
    
    // if can't find both - send a warning message
    if (!theNodes[0] || !theNodes[1])  {
        if (!theNodes[0])  {
            opserr << "WARNING ElastomericBearingBoucWenMod3d::setDomain() - Nd1: "
                << connectedExternalNodes(0)
                << " does not exist in the model for";
        } else  {
            opserr << "WARNING ElastomericBearingBoucWenMod3d::setDomain() - Nd2: "
                << connectedExternalNodes(1)
                << " does not exist in the model for";
        }
        opserr << " element: " << this->getTag() << ".\n";
        
        return;
    }
    
    // now determine the number of dof and the dimension
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();
    
    // if differing dof at the ends - print a warning message
    if (dofNd1 != 6)  {
        opserr << "ElastomericBearingBoucWenMod3d::setDomain() - node 1: "
            << connectedExternalNodes(0)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    if (dofNd2 != 6)  {
        opserr << "ElastomericBearingBoucWenMod3d::setDomain() - node 2: "
            << connectedExternalNodes(1)
            << " has incorrect number of DOF (not 6).\n";
        return;
    }
    
    // call the base class method
    this->DomainComponent::setDomain(theDomain);
    
    // set up the transformation matrix for orientation
    this->setUp();
}


int ElastomericBearingBoucWenMod3d::commitState()
{
    int errCode = 0;
    
    // commit trial history variables
    ubC = ub;
    zC  = z;
    
    // commit the base class
    errCode += this->Element::commitState();
    
    return errCode;
}


int ElastomericBearingBoucWenMod3d::revertToLastCommit()
{
    int errCode = 0;
    
    return errCode;
}


int ElastomericBearingBoucWenMod3d::revertToStart()
{   
    int errCode=0;
    
    // reset trial history variables
    ub.Zero();
    z.Zero();
    qb.Zero();
    
    // reset committed history variables
    ubC.Zero();
    zC.Zero();
    
    // reset tangent of hysteretic evolution parameters
    dzdu(0,0) = dzdu(1,1) = A*k0/qYield;
    dzdu(1,0) = dzdu(0,1) = 0.0;
    
    // reset stiffness matrix in basic system
    kb = kbInit;
    
    return errCode;
}


int ElastomericBearingBoucWenMod3d::update()
{
    // get global trial displacements and velocities
    const Vector &dsp1 = theNodes[0]->getTrialDisp();
    const Vector &dsp2 = theNodes[1]->getTrialDisp();
    const Vector &vel1 = theNodes[0]->getTrialVel();
    const Vector &vel2 = theNodes[1]->getTrialVel();
    
    static Vector ug(12), ugdot(12), uldot(12), ubdot(6);
    for (int i=0; i<6; i++)  {
        ug(i)   = dsp1(i);  ugdot(i)   = vel1(i);
        ug(i+6) = dsp2(i);  ugdot(i+6) = vel2(i);
    }
    
    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);
    uldot.addMatrixVector(0.0, Tgl, ugdot, 1.0);
    
    // transform response from the local to the basic system
    ub.addMatrixVector(0.0, Tlb, ul, 1.0);
    ubdot.addMatrixVector(0.0, Tlb, uldot, 1.0);
    
    // 1) get axial force and stiffness in basic x-direction
            kb(0,0)=Kv0;
			qb(0)=kb(0,0)*ub(0);
    
    // 2) calculate shear forces and stiffnesses in basic y- and z-direction
    // get displacement increments (trial - committed)
    Vector delta_ub = ub - ubC;
	// calculate dependency parameters
	double alphaT = b1+b2*T+b3*pow(T,2)+b4*pow(T,3);
	double Ac=(PI/4)*(D2*D2-D1*D1);
    double PM=fabs(qb(0)/Ac);
    double alphaP = a1*log(PM)+a2;
	qYield = (1.0-alpha1-alpha2*pow(fy/(kInit*alphaT*alphaP),mu-1.0))*fy;
	k0 = (1.0-alpha1)*kInit*alphaT*alphaP;
	k2 = alpha1*kInit*alphaT*alphaP;
    k3 = alpha2*kInit*alphaT*alphaP;    
    
	if (sqrt(pow(delta_ub(1),2)+pow(delta_ub(2),2)) > 0.0)  {
    
        // get yield displacement
        double uy = qYield/k0;
        
        // calculate hysteretic evolution parameter z using Newton-Raphson
        int iter = 0;
        double zNrm, tmp1, tmp2, tmp3, tmp4;
        Vector f(2), delta_z(2);
        Matrix Df(2,2);
        do  {
            zNrm = z.Norm();
            if (zNrm == 0.0)  // check because of negative exponents
                zNrm = DBL_EPSILON;
            tmp1 = z(0)*delta_ub(1) + z(1)*delta_ub(2);
            tmp2 = gamma + beta*sgn(tmp1);
            tmp3 = pow(zNrm,eta-2.0)*tmp1*tmp2;
            tmp4 = pow(zNrm,eta-4.0)*tmp2;
            
            // function and derivative
            f(0) = z(0) - zC(0) - 1.0/uy*(A*delta_ub(1) - z(0)*tmp3);
            f(1) = z(1) - zC(1) - 1.0/uy*(A*delta_ub(2) - z(1)*tmp3);
            
            Df(0,0) = 1.0 + tmp4/uy*(pow(z(1),3)*delta_ub(2) + 2.0*pow(z(1),2)*z(0)*delta_ub(1)
                + z(1)*pow(z(0),2)*delta_ub(2)*(eta-1.0) + pow(z(0),3)*delta_ub(1)*eta);
            Df(1,0) = z(1)*tmp4/uy*(pow(z(1),2)*delta_ub(1) + z(0)*z(1)*delta_ub(2)*(eta-2.0)
                + pow(z(0),2)*delta_ub(1)*(eta-1.0));
            Df(0,1) = z(0)*tmp4/uy*(pow(z(0),2)*delta_ub(2) + z(0)*z(1)*delta_ub(1)*(eta-2.0)
                + pow(z(1),2)*delta_ub(2)*(eta-1.0));
            Df(1,1) = 1.0 + tmp4/uy*(pow(z(0),3)*delta_ub(1) + 2.0*pow(z(0),2)*z(1)*delta_ub(2)
                + z(0)*pow(z(1),2)*delta_ub(1)*(eta-1.0) + pow(z(1),3)*delta_ub(2)*eta);
            
            // issue warning if diagonal of derivative Df is zero
            if ((fabs(Df(0,0)) <= DBL_EPSILON) || (fabs(Df(1,1)) <= DBL_EPSILON))  {
                opserr << "WARNING: ElastomericBearingBoucWenMod3d::update() - "
                    << "zero Jacobian in Newton-Raphson scheme for hysteretic "
                    << "evolution parameter z.\n";
                return -1;
            }
            
            // advance one step
            delta_z = f/Df;
            z -= delta_z;
            iter++;
        } while ((delta_z.Norm() >= tol) && (iter < maxIter));
        
        // issue warning if Newton-Raphson scheme did not converge
        if (iter >= maxIter)   {
            opserr << "WARNING: ElastomericBearingBoucWenMod3d::update() - "
                << "did not find the hysteretic evolution parameters z after "
                << iter << " iterations and norm: " << delta_z.Norm() << endln;
            return -2;
        }
        
        // get derivative of hysteretic evolution parameter
        delta_z = z-zC;
        if (fabs(delta_ub(1)) > 0.0)  {
            dzdu(0,0) = delta_z(0)/delta_ub(1);
            dzdu(1,0) = delta_z(1)/delta_ub(1);
        }
        if (fabs(delta_ub(2)) > 0.0)  {
            dzdu(0,1) = delta_z(0)/delta_ub(2);
            dzdu(1,1) = delta_z(1)/delta_ub(2);
        }
		
		
       
		qb(1) = qYield*z(0) + k2*ub(1) + k3*sgn(ub(1))*pow(fabs(ub(1)),mu);
        qb(2) = qYield*z(1) + k2*ub(2) + k3*sgn(ub(2))*pow(fabs(ub(2)),mu);
        // set tangent stiffness
        kb(1,1) = qYield*dzdu(0,0) + k2 + k3*mu*pow(fabs(ub(1)),mu-1.0);
        kb(1,2) = qYield*dzdu(0,1);
        kb(2,1) = qYield*dzdu(1,0);
        kb(2,2) = qYield*dzdu(1,1) + k2 + k3*mu*pow(fabs(ub(2)),mu-1.0);
    }
    
    // 3) get moment and stiffness in basic x-direction
    qb(3) = Kt*ub(3);
    kb(3,3) = Kt;
    
    // 4) get moment and stiffness in basic y-direction
    qb(4) = Kr*ub(4);
    kb(4,4) = Kr;
    
    // 5) get moment and stiffness in basic z-direction
    qb(5) = Kr*ub(5);
    kb(5,5) = Kr;
    
    return 0;
}


const Matrix& ElastomericBearingBoucWenMod3d::getTangentStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kb, 1.0);
    
    // add geometric stiffness to local stiffness
    double kGeo1 = 0.5*qb(0);
    kl(5,1)   -= kGeo1;
    kl(5,7)   += kGeo1;
    kl(11,1)  -= kGeo1;
    kl(11,7)  += kGeo1;
    kl(4,2)   += kGeo1;
    kl(4,8)   -= kGeo1;
    kl(10,2)  += kGeo1;
    kl(10,8)  -= kGeo1;
    double kGeo2 = kGeo1*shearDistI*L;
    kl(5,5)   += kGeo2;
    kl(11,5)  -= kGeo2;
    kl(4,4)   += kGeo2;
    kl(10,4)  -= kGeo2;
    double kGeo3 = kGeo1*(1.0 - shearDistI)*L;
    kl(5,11)  -= kGeo3;
    kl(11,11) += kGeo3;
    kl(4,10)  -= kGeo3;
    kl(10,10) += kGeo3;
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWenMod3d::getInitialStiff()
{
    // zero the matrix
    theMatrix.Zero();
    
    // transform from basic to local system
    static Matrix kl(12,12);
    kl.addMatrixTripleProduct(0.0, Tlb, kbInit, 1.0);
    
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWenMod3d::getDamp()
{
    // zero the matrix
    theMatrix.Zero();
    
    // call base class to setup Rayleigh damping
    theMatrix = this->Element::getDamp();
    double factThis = 0.0;
    
    // now add damping tangent from materials
    static Matrix cb(6,6);
    cb.Zero();

    // transform from basic to local system
    static Matrix cl(12,12);
    cl.addMatrixTripleProduct(0.0, Tlb, cb, 1.0);
    
    // transform from local to global system and add to cg
    theMatrix.addMatrixTripleProduct(factThis, Tgl, cl, 1.0);
    
    return theMatrix;
}


const Matrix& ElastomericBearingBoucWenMod3d::getMass()
{
    // zero the matrix
    theMatrix.Zero();
    
    // check for quick return
    if (mass == 0.0)  {
        return theMatrix;
    }
    
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theMatrix(i,i)     = m;
        theMatrix(i+6,i+6) = m;
    }
    
    return theMatrix;
}


void ElastomericBearingBoucWenMod3d::zeroLoad()
{
    theLoad.Zero();
}


int ElastomericBearingBoucWenMod3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr <<"ElastomericBearingBoucWenMod3d::addLoad() - "
        << "load type unknown for element: "
        << this->getTag() << ".\n";
    
    return -1;
}


int ElastomericBearingBoucWenMod3d::addInertiaLoadToUnbalance(const Vector &accel)
{
    // check for quick return
    if (mass == 0.0)  {
        return 0;
    }
    
    // get R * accel from the nodes
    const Vector &Raccel1 = theNodes[0]->getRV(accel);
    const Vector &Raccel2 = theNodes[1]->getRV(accel);
    
    if (6 != Raccel1.Size() || 6 != Raccel2.Size())  {
        opserr << "ElastomericBearingBoucWenMod3d::addInertiaLoadToUnbalance() - "
            << "matrix and vector sizes are incompatible.\n";
        return -1;
    }
    
    // want to add ( - fact * M R * accel ) to unbalance
    // take advantage of lumped mass matrix
    double m = 0.5*mass;
    for (int i=0; i<3; i++)  {
        theLoad(i)   -= m * Raccel1(i);
        theLoad(i+6) -= m * Raccel2(i);
    }
    
    return 0;
}


const Vector& ElastomericBearingBoucWenMod3d::getResistingForce()
{
    // zero the residual
    theVector.Zero();
    
    // determine resisting forces in local system
    static Vector ql(12);
    ql.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
    
    // add P-Delta moments to local forces
    double kGeo1 = 0.5*qb(0);
    double MpDelta1 = kGeo1*(ul(7)-ul(1));
    ql(5)  += MpDelta1;
    ql(11) += MpDelta1;
    double MpDelta2 = kGeo1*shearDistI*L*ul(5);
    ql(5)  += MpDelta2;
    ql(11) -= MpDelta2;
    double MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(11);
    ql(5)  -= MpDelta3;
    ql(11) += MpDelta3;
    double MpDelta4 = kGeo1*(ul(8)-ul(2));
    ql(4)  -= MpDelta4;
    ql(10) -= MpDelta4;
    double MpDelta5 = kGeo1*shearDistI*L*ul(4);
    ql(4)  += MpDelta5;
    ql(10) -= MpDelta5;
    double MpDelta6 = kGeo1*(1.0 - shearDistI)*L*ul(10);
    ql(4)  -= MpDelta6;
    ql(10) += MpDelta6;
    
    // determine resisting forces in global system
    theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);
    
    // subtract external load
    theVector.addVector(1.0, theLoad, -1.0);
    
    return theVector;
}


const Vector& ElastomericBearingBoucWenMod3d::getResistingForceIncInertia()
{
    // this already includes damping forces from materials
    theVector = this->getResistingForce();
    
    // add the damping forces from rayleigh damping
    if (addRayleigh == 1)  {
        if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
            theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    }
    
    // add inertia forces from element mass
    if (mass != 0.0)  {
        const Vector &accel1 = theNodes[0]->getTrialAccel();
        const Vector &accel2 = theNodes[1]->getTrialAccel();
        
        double m = 0.5*mass;
        for (int i=0; i<3; i++)  {
            theVector(i)   += m * accel1(i);
            theVector(i+6) += m * accel2(i);
        }
    }
    
    return theVector;
}


int ElastomericBearingBoucWenMod3d::sendSelf(int commitTag, Channel &sChannel)
{
    // send element parameters
    static Vector data(17);
    data(0) = this->getTag();
    data(1) = k0;
    data(2) = qYield;
    data(3) = k2;
    data(4) = k3;
    data(5) = mu;
    data(6) = eta;
    data(7) = beta;
    data(8) = gamma;
    data(9) = A;
    data(10) = shearDistI;
    data(11) = addRayleigh;
    data(12) = mass;
    data(13) = maxIter;
    data(14) = tol;
    data(15) = x.Size();
    data(16) = y.Size();
    sChannel.sendVector(0, commitTag, data);
    
    // send the two end nodes
    sChannel.sendID(0, commitTag, connectedExternalNodes);
    
    
    
    
    // send remaining data
    if (x.Size() == 3)
        sChannel.sendVector(0, commitTag, x);
    if (y.Size() == 3)
        sChannel.sendVector(0, commitTag, y);
    
    return -1;
}


int ElastomericBearingBoucWenMod3d::recvSelf(int commitTag, Channel &rChannel,
    FEM_ObjectBroker &theBroker)
{
    
    
    // receive element parameters
    static Vector data(17);
    rChannel.recvVector(0, commitTag, data);
    this->setTag((int)data(0));
    k0 = data(1);
    qYield = data(2);
    k2 = data(3);
    k3 = data(4);
    mu = data(5);
    eta = data(6);
    beta = data(7);
    gamma = data(8);
    A = data(9);
    shearDistI = data(10);
    addRayleigh = (int)data(11);
    mass = data(12);
    maxIter = (int)data(13);
    tol = data(14);
    
    // receive the two end nodes
    rChannel.recvID(0, commitTag, connectedExternalNodes);
    
   
    
    // receive remaining data
    if ((int)data(15) == 3)  {
        x.resize(3);
        rChannel.recvVector(0, commitTag, x);
    }
    if ((int)data(16) == 3)  {
        y.resize(3);
        rChannel.recvVector(0, commitTag, y);
    }
    onP0 = false;
    
    // initialize initial stiffness matrix
    kbInit.Zero();
    kbInit.Zero();
    kbInit(0,0) = Kv0;
    kbInit(1,1) = kbInit(2,2) = A*k0 + k2;
    kbInit(3,3) = Kt;
    kbInit(4,4) = Kr;
    kbInit(5,5) = Kr;
    
    // initialize variables
    this->revertToStart();
    
    return -1;
}


int ElastomericBearingBoucWenMod3d::displaySelf(Renderer &theViewer,
    int displayMode, float fact)
{
    // first determine the end points of the element based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    
    static Vector v1(3);
    static Vector v2(3);
    
    if (displayMode >= 0)  {
        const Vector &end1Disp = theNodes[0]->getDisp();
        const Vector &end2Disp = theNodes[1]->getDisp();
        
        for (int i=0; i<3; i++)  {
            v1(i) = end1Crd(i) + end1Disp(i)*fact;
            v2(i) = end2Crd(i) + end2Disp(i)*fact;
        }
    } else  {
        int mode = displayMode * -1;
        const Matrix &eigen1 = theNodes[0]->getEigenvectors();
        const Matrix &eigen2 = theNodes[1]->getEigenvectors();
        
        if (eigen1.noCols() >= mode)  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
                v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;
            }
        } else  {
            for (int i=0; i<3; i++)  {
                v1(i) = end1Crd(i);
                v2(i) = end2Crd(i);
            }
        }
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0, this->getTag(), 0);
}


void ElastomericBearingBoucWenMod3d::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) {
        // print everything
        s << "Element: " << this->getTag() << endln; 
        s << "  type: ElastomericBearingBoucWenMod3d\n";
        s << "  iNode: " << connectedExternalNodes(0);
        s << "  jNode: " << connectedExternalNodes(1) << endln;
        s << "  k0: " << k0 << "  qYield: " << qYield << "  k2: " << k2 << endln;
        s << "  k3: " << k3 << "  mu: " << mu << endln;
        s << "  eta: " << eta << "  beta: " << beta << "  gamma: " << gamma << endln;
        s << "  shearDistI: " << shearDistI << "  addRayleigh: "
            << addRayleigh << "  mass: " << mass << endln;
        s << "  maxIter: " << maxIter << "  tol: " << tol << endln;
        // determine resisting forces in global system
        s << "  resisting force: " << this->getResistingForce() << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ElastomericBearingBoucWenMod3d\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"k0\": " << k0 << ", ";
        s << "\"qYield\": " << qYield << ", ";
        s << "\"k2\": " << k2 << ", ";
        s << "\"k3\": " << k3 << ", ";
        s << "\"mu\": " << mu << ", ";
        s << "\"eta\": " << eta << ", ";
        s << "\"beta\": " << beta << ", ";
        s << "\"gamma\": " << gamma << ", ";
        s << "\"shearDistI\": " << shearDistI << ", ";
        s << "\"addRayleigh\": " << addRayleigh << ", ";
        s << "\"mass\": " << mass << "}";
    }
}


Response* ElastomericBearingBoucWenMod3d::setResponse(const char **argv, int argc,
    OPS_Stream &output)
{
    Response *theResponse = 0;
    
    output.tag("ElementOutput");
    output.attr("eleType","ElastomericBearingBoucWenMod3d");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);
    
    // global forces
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
        strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0)
    {
        output.tag("ResponseType","Px_1");
        output.tag("ResponseType","Py_1");
        output.tag("ResponseType","Pz_1");
        output.tag("ResponseType","Mx_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Mz_1");
        output.tag("ResponseType","Px_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","Mx_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 1, theVector);
    }
    // local forces
    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0)
    {
        output.tag("ResponseType","N_ 1");
        output.tag("ResponseType","Vy_1");
        output.tag("ResponseType","Vz_1");
        output.tag("ResponseType","T_1");
        output.tag("ResponseType","My_1");
        output.tag("ResponseType","Tz_1");
        output.tag("ResponseType","N_2");
        output.tag("ResponseType","Py_2");
        output.tag("ResponseType","Pz_2");
        output.tag("ResponseType","T_2");
        output.tag("ResponseType","My_2");
        output.tag("ResponseType","Mz_2");
        
        theResponse = new ElementResponse(this, 2, theVector);
    }
    // basic forces
    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0)
    {
        output.tag("ResponseType","qb1");
        output.tag("ResponseType","qb2");
        output.tag("ResponseType","qb3");
        output.tag("ResponseType","qb4");
        output.tag("ResponseType","qb5");
        output.tag("ResponseType","qb6");
        
        theResponse = new ElementResponse(this, 3, Vector(6));
    }
    // local displacements
    else if (strcmp(argv[0],"localDisplacement") == 0 ||
        strcmp(argv[0],"localDisplacements") == 0)
    {
        output.tag("ResponseType","ux_1");
        output.tag("ResponseType","uy_1");
        output.tag("ResponseType","uz_1");
        output.tag("ResponseType","rx_1");
        output.tag("ResponseType","ry_1");
        output.tag("ResponseType","rz_1");
        output.tag("ResponseType","ux_2");
        output.tag("ResponseType","uy_2");
        output.tag("ResponseType","uz_2");
        output.tag("ResponseType","rx_2");
        output.tag("ResponseType","ry_2");
        output.tag("ResponseType","rz_2");
        
        theResponse = new ElementResponse(this, 4, theVector);
    }
    // basic displacements
    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0 || 
        strcmp(argv[0],"basicDeformation") == 0 || strcmp(argv[0],"basicDeformations") == 0 ||
        strcmp(argv[0],"basicDisplacement") == 0 || strcmp(argv[0],"basicDisplacements") == 0)
    {
        output.tag("ResponseType","ub1");
        output.tag("ResponseType","ub2");
        output.tag("ResponseType","ub3");
        output.tag("ResponseType","ub4");
        output.tag("ResponseType","ub5");
        output.tag("ResponseType","ub6");
        
        theResponse = new ElementResponse(this, 5, Vector(6));
    }
    // hysteretic evolution parameter
    else if (strcmp(argv[0],"hystereticParameter") == 0 || strcmp(argv[0],"hystParameter") == 0 || 
        strcmp(argv[0],"hystereticParam") == 0 || strcmp(argv[0],"hystParam") == 0 ||
        strcmp(argv[0],"z") == 0)
    {
        output.tag("ResponseType","z1");
        output.tag("ResponseType","z2");
        
        theResponse = new ElementResponse(this, 6, Vector(2));
    }
    
    
    output.endTag(); // ElementOutput
    
    return theResponse;
}


int ElastomericBearingBoucWenMod3d::getResponse(int responseID, Information &eleInfo)
{
    double kGeo1, MpDelta1, MpDelta2, MpDelta3, MpDelta4, MpDelta5, MpDelta6;
    
    switch (responseID)  {
    case 1:  // global forces
        return eleInfo.setVector(this->getResistingForce());
        
    case 2:  // local forces
        theVector.Zero();
        // determine resisting forces in local system
        theVector.addMatrixTransposeVector(0.0, Tlb, qb, 1.0);
        // add P-Delta moments
        kGeo1 = 0.5*qb(0);
        MpDelta1 = kGeo1*(ul(7)-ul(1));
        theVector(5)  += MpDelta1;
        theVector(11) += MpDelta1;
        MpDelta2 = kGeo1*shearDistI*L*ul(5);
        theVector(5)  += MpDelta2;
        theVector(11) -= MpDelta2;
        MpDelta3 = kGeo1*(1.0 - shearDistI)*L*ul(11);
        theVector(5)  -= MpDelta3;
        theVector(11) += MpDelta3;
        MpDelta4 = kGeo1*(ul(8)-ul(2));
        theVector(4)  -= MpDelta4;
        theVector(10) -= MpDelta4;
        MpDelta5 = kGeo1*shearDistI*L*ul(4);
        theVector(4)  += MpDelta5;
        theVector(10) -= MpDelta5;
        MpDelta6 = kGeo1*(1.0 - shearDistI)*L*ul(10);
        theVector(4)  -= MpDelta6;
        theVector(10) += MpDelta6;
        return eleInfo.setVector(theVector);
        
    case 3:  // basic forces
        return eleInfo.setVector(qb);
        
    case 4:  // local displacements
        return eleInfo.setVector(ul);
        
    case 5:  // basic displacements
        return eleInfo.setVector(ub);
        
    case 6:  // hysteretic evolution parameter
        return eleInfo.setVector(z);
        
    default:
        return -1;
    }
}


// set up the transformation matrix for orientation
void ElastomericBearingBoucWenMod3d::setUp()
{
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	
    Vector xp = end2Crd - end1Crd;
    L = xp.Norm();
    
    if (L > DBL_EPSILON)  {
        if (x.Size() == 0)  {
            x.resize(3);
            x = xp;
        } else if (onP0)  {
            opserr << "WARNING ElastomericBearingBoucWenMod3d::setUp() - "
                << "element: " << this->getTag()
                << " - ignoring nodes and using specified "
                << "local x vector to determine orientation.\n";
        }
    }
    // check that vectors for orientation are of correct size
    if (x.Size() != 3 || y.Size() != 3)  {
        opserr << "ElastomericBearingBoucWenMod3d::setUp() - "
            << "element: " << this->getTag() << endln
            << " - incorrect dimension of orientation vectors.\n";
        exit(-1);
    }
    
    // establish orientation of element for the transformation matrix
    // z = x cross y
    Vector z(3);
    z(0) = x(1)*y(2) - x(2)*y(1);
    z(1) = x(2)*y(0) - x(0)*y(2);
    z(2) = x(0)*y(1) - x(1)*y(0);
    
    // y = z cross x
    y(0) = z(1)*x(2) - z(2)*x(1);
    y(1) = z(2)*x(0) - z(0)*x(2);
    y(2) = z(0)*x(1) - z(1)*x(0);
    
    // compute length(norm) of vectors
    double xn = x.Norm();
    double yn = y.Norm();
    double zn = z.Norm();
    
    // check valid x and y vectors, i.e. not parallel and of zero length
    if (xn == 0 || yn == 0 || zn == 0)  {
        opserr << "ElastomericBearingBoucWenMod3d::setUp() - "
            << "element: " << this->getTag() << endln
            << " - invalid orientation vectors.\n";
        exit(-1);
    }
    
    // create transformation matrix from global to local system
    Tgl.Zero();
    Tgl(0,0) = Tgl(3,3) = Tgl(6,6) = Tgl(9,9)   = x(0)/xn;
    Tgl(0,1) = Tgl(3,4) = Tgl(6,7) = Tgl(9,10)  = x(1)/xn;
    Tgl(0,2) = Tgl(3,5) = Tgl(6,8) = Tgl(9,11)  = x(2)/xn;
    Tgl(1,0) = Tgl(4,3) = Tgl(7,6) = Tgl(10,9)  = y(0)/yn;
    Tgl(1,1) = Tgl(4,4) = Tgl(7,7) = Tgl(10,10) = y(1)/yn;
    Tgl(1,2) = Tgl(4,5) = Tgl(7,8) = Tgl(10,11) = y(2)/yn;
    Tgl(2,0) = Tgl(5,3) = Tgl(8,6) = Tgl(11,9)  = z(0)/zn;
    Tgl(2,1) = Tgl(5,4) = Tgl(8,7) = Tgl(11,10) = z(1)/zn;
    Tgl(2,2) = Tgl(5,5) = Tgl(8,8) = Tgl(11,11) = z(2)/zn;
    
    // create transformation matrix from local to basic system (linear)
    Tlb.Zero();
    Tlb(0,0) = Tlb(1,1) = Tlb(2,2) = Tlb(3,3) = Tlb(4,4) = Tlb(5,5) = -1.0;
    Tlb(0,6) = Tlb(1,7) = Tlb(2,8) = Tlb(3,9) = Tlb(4,10) = Tlb(5,11) = 1.0;
    Tlb(1,5) = -shearDistI*L;
    Tlb(1,11) = -(1.0 - shearDistI)*L;
    Tlb(2,4) = -Tlb(1,5);
    Tlb(2,10) = -Tlb(1,11);
}


double ElastomericBearingBoucWenMod3d::sgn(double x)
{
    if (x > 0)
        return 1.0;
    else if (x < 0)
        return -1.0;
    else
        return 0.0;
}
