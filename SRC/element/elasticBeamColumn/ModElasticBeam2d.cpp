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
// $URL$
                                                                        
// Written: fmk 11/95
// Revised: Dimitrios G. Lignos 01/31/2011
//
// Purpose: This file contains the class definition for ModElasticBeam2d.
// ModElasticBeam2d is a 2d beam element that considers modified elastic stiffness
// coefficients for nonprismatic sections
// Modification is done for "n" stiff springs adjustments 
// based on Zareian and Medina (2010), Computers and Structures, Vol. 88(1-2)

#include "ModElasticBeam2d.h"
#include <elementAPI.h>

#include <ElementalLoad.h>

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>

Matrix ModElasticBeam2d::K(6,6);
Vector ModElasticBeam2d::P(6);
Matrix ModElasticBeam2d::kb(3,3);


static int numModElasticBeam2d = 0;

void *
OPS_ModElasticBeam2d()
{
  // print out a message about who wrote this element & any copyright info wanted
  if (numModElasticBeam2d == 0) {
    opserr << "ModElasticBeam2d element -> for Stiffness Modification Factors by D.Lignos";
    numModElasticBeam2d++;
  }

  Element *theEle = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingArgs == 0) { // parallel processing
    theEle = new ModElasticBeam2d();
    return theEle;
  }

  if (numRemainingArgs < 10) {
    opserr << "ERROR not enough args provided, want: element ModElasticBeam2d tag? iNode? jNode? A? E? I? K11? K33? K44? transfType? <-alpha $alpha> <-d $d> <-rho $rho> <-cMass>\n";
    return 0;
  }

  int numData;
  int iData[5];     //tag, iNode, jNode, transfTag, cMass
  double dData[9];  // A, E, I, K11, K33, K44, alpha, d, rho

  iData[4] = 0;  // cMass
  dData[6] = 0.0; dData[7] = 0.0;  // alpha and d
  dData[8] = 0.0;  // rho

  numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid element data (tag, iNode, jNode) element ElasticBeamColumn2d\n";
    return 0;
  }

  int eleTag = iData[0];

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING error reading element data (A, E, I) element ElasticBeamColumn2d " << eleTag << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, &dData[3]) != 0) {
    opserr << "WARNING error reading element data (K11, K33, K44) element ElasticBeamColumn2d " << eleTag << endln;
    return 0;
  }
	
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
	opserr << "WARNING error reading data (transfTag) element ElasticBeamColumn2d " << eleTag << endln;
	return 0;
  }
	
  numRemainingArgs = OPS_GetNumRemainingInputArgs();
  while (numRemainingArgs > 1) {
    const char *argvLoc = OPS_GetString();;
    numData = 1;

    if ((strcmp(argvLoc, "-alpha") == 0) || (strcmp(argvLoc, "-Alpha") == 0) || (strcmp(argvLoc, "-ALPHA") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[6]) != 0) {
	opserr << "WARNING error reading element data (alpha) element ElasticBeamColumn2d " << eleTag << endln;
	return 0;
      }
    } else if ((strcmp(argvLoc, "-d") == 0) || (strcmp(argvLoc, "-D") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
	opserr << "WARNING error reading element data (D) element ElasticBeamColumn2d " << eleTag << endln;
	return 0;
      }
    } else if ((strcmp(argvLoc, "-rho") == 0) || (strcmp(argvLoc, "Rho") == 0) || (strcmp(argvLoc, "-RHO") == 0)) {
      if (OPS_GetDoubleInput(&numData, &dData[8]) != 0) {
	opserr << "WARNING error reading element data (rho) element ElasticBeamColumn2d " << eleTag << endln;
	return 0;
      }
    } else if ((strcmp(argvLoc, "-lMass") == 0) || (strcmp(argvLoc, "lMass") == 0)) {
      iData[4] = 0;  // lumped mass matrix (default)
    } else if ((strcmp(argvLoc, "-cMass") == 0) || (strcmp(argvLoc, "cMass") == 0)) {
      iData[4] = 1;  // consistent mass matrix
    }
    numRemainingArgs = OPS_GetNumRemainingInputArgs();      
  }

  CrdTransf *theTransf = OPS_getCrdTransf(iData[3]);
  if (theTransf == 0) {
    opserr << "WARNING error could not find a transformation with tag: " << iData[3] << "element ElasticBeamColumn2d " << eleTag << endln;
    return 0;
  } 

  theEle = new ModElasticBeam2d(iData[0], dData[0], dData[1], dData[2],
				iData[1], iData[2],
				dData[3], dData[4], dData[5],
				*theTransf,
				dData[6], dData[7], dData[8], iData[4]);

  //  delete theTransf;

  return theEle;
}


ModElasticBeam2d::ModElasticBeam2d()
  :Element(0,ELE_TAG_ModElasticBeam2d), 
   A(0.0), E(0.0), I(0.0), 
   K11(0.0), K33(0.0), K44(0.0), 
   alpha(0.0), d(0.0), rho(0.0), cMass(0),
   Q(6), q(3), 
   connectedExternalNodes(2), theCoordTransf(0)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}


// modified by Dimitrios Lignos for modified Stiffness Coefficients
ModElasticBeam2d::ModElasticBeam2d(int tag, double a, double e, double i, 
				   int Nd1, int Nd2, 
				   double k11, double k33, double k44, CrdTransf &coordTransf,
				   double Alpha, double depth, double r, int cm)
  :Element(tag,ELE_TAG_ModElasticBeam2d), 
   A(a), E(e), I(i), 
   K11(k11), K33(k33), K44(k44), 
   alpha(Alpha), d(depth), rho(r), cMass(cm),
   Q(6), q(3),
  connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  theCoordTransf = coordTransf.getCopy2d();

  if (!theCoordTransf) {
    opserr << "ModElasticBeam2d::ModElasticBeam2d -- failed to get copy of coordinate transformation\n";
    exit(01);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

ModElasticBeam2d::~ModElasticBeam2d()
{
    if (theCoordTransf)
	delete theCoordTransf;
}

int
ModElasticBeam2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ModElasticBeam2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ModElasticBeam2d::getNodePtrs(void) 
{
  return theNodes;
}

int
ModElasticBeam2d::getNumDOF(void)
{
    return 6;
}

void
ModElasticBeam2d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ModElasticBeam2d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "ModElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ModElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "ModElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "ModElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ModElasticBeam2d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ModElasticBeam2d::setDomain -- Element has zero length\n";
      exit(-1);
    }
}

int
ModElasticBeam2d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ModElasticBeam2d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
ModElasticBeam2d::revertToLastCommit()
{
    return theCoordTransf->revertToLastCommit();
}

int
ModElasticBeam2d::revertToStart()
{
    return theCoordTransf->revertToStart();
}

int
ModElasticBeam2d::update(void)
{
  return theCoordTransf->update();
}

const Matrix &
ModElasticBeam2d::getTangentStiff(void)
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = K44*I*EoverL;		// 2EI/L
  double EIoverL4 = K11*I*EoverL;       // 4EI/L
 // added by Dimitrios Lignos
  double EIoverL6 = K33*I*EoverL;		// 4EI/L
  
  // determine q = kv + q0
  q(0) = EAoverL*v(0);
  q(1) = EIoverL4*v(1) + EIoverL2*v(2);
  //q(2) = EIoverL2*v(1) + EIoverL4*v(2);
  q(2) = EIoverL2*v(1) + EIoverL6*v(2);

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  kb(0,0) = EAoverL;
  kb(1,1) = EIoverL4;
  //kb(2,2) = EIoverL4;
  kb(2,2) = EIoverL6;
  kb(2,1) = kb(1,2) = EIoverL2;
  
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
ModElasticBeam2d::getInitialStiff(void)
{
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = K44*I*EoverL;		// 2EI/L
  double EIoverL4 = K11*I*EoverL;	    // 4EI/L
  // Added By Dimitrios Lignos
  double EIoverL6 = K33*I*EoverL;       // 4EI/L
  
  kb(0,0) = EAoverL;
  kb(1,1) = EIoverL4; 
  //kb(2,2) = EIoverL4;
  kb(2,2) = EIoverL6;
  
  kb(2,1) = kb(1,2) = EIoverL2;
  
  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
ModElasticBeam2d::getMass(void)
{ 
    K.Zero();
    
    if (rho > 0.0)  {
        // get initial element length
        double L = theCoordTransf->getInitialLength();
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*L;
            K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;
        } else  {
            // consistent mass matrix
            static Matrix ml(6,6);
            double m = rho*L/420.0;
            ml(0,0) = ml(3,3) = m*140.0;
            ml(0,3) = ml(3,0) = m*70.0;

            ml(1,1) = ml(4,4) = m*156.0;
            ml(1,4) = ml(4,1) = m*54.0;
            ml(2,2) = ml(5,5) = m*4.0*L*L;
            ml(2,5) = ml(5,2) = -m*3.0*L*L;
            ml(1,2) = ml(2,1) = m*22.0*L;
            ml(4,5) = ml(5,4) = -ml(1,2);
            ml(1,5) = ml(5,1) = -m*13.0*L;
            ml(2,4) = ml(4,2) = -ml(1,5);
            
            // transform local mass matrix to global system
            K = theCoordTransf->getGlobalMatrixFromLocal(ml);
        }
    }
    
    return K;
}

void 
ModElasticBeam2d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  return;
}

int 
ModElasticBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wt*L;
    double M = V*L/6.0; // wt*L*L/12
    double P = wa*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    q0[1] -= M;
    q0[2] += M;
  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    double M1 = -a * b2 * P * L2;
    double M2 = a2 * b * P * L2;
    q0[1] += M1;
    q0[2] += M2;
  }
  
  else if (type == LOAD_TAG_Beam2dTempLoad) {
    double Ttop1 = data(0)* loadFactor;
    double Tbot1 = data(1)* loadFactor;
    double Ttop2 = data(2)* loadFactor;
    double Tbot2 = data(3)* loadFactor;
        
    // fixed end forces due to a linear thermal load
    double dT1 = Ttop1-Tbot1;
    double dT = (Ttop2-Tbot2)-(Ttop1-Tbot1);
    double a = alpha/d;  // constant based on temp difference at top and bottom, 
    // coefficient of thermal expansion and beam depth
    double M1 = a*E*I*(-dT1+(4.0/3.0)*dT); //Fixed End Moment end 1
    double M2 = a*E*I*(dT1+(5.0/3.0)*dT); //Fixed End Moment end 2
    double F = alpha*(((Ttop2+Ttop1)/2+(Tbot2+Tbot1)/2)/2)*E*A; // Fixed End Axial Force
    double M1M2divL =(M1+M2)/L; // Fixed End Shear
    
    // Reactions in basic system
    p0[0] += 0;
    p0[1] += M1M2divL;
    p0[2] -= M1M2divL;

    // Fixed end forces in basic system
    q0[0] -= F;
    q0[1] += M1;
    q0[2] += M2;
  }

  else {
    opserr << "ModElasticBeam2d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
ModElasticBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "ModElasticBeam2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);

    Q(3) -= m * Raccel2(0);
    Q(4) -= m * Raccel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(6);
    for (int i=0; i<3; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

const Vector &
ModElasticBeam2d::getResistingForceIncInertia()
{	
  P = this->getResistingForce();
  
  // subtract external load P = P - Q;
  P.addVector(1.0, Q, -1.0);
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
  
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5*rho*L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);

    P(3) += m * accel2(0);
    P(4) += m * accel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(6);
    for (int i=0; i<3; i++)  {
      accel(i)   = accel1(i);
      accel(i+3) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
  
  return P;
}


const Vector &
ModElasticBeam2d::getResistingForce()
{
  theCoordTransf->update();
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;			// EA/L
  double EIoverL2 = K44*I*EoverL;		// 2EI/L
  double EIoverL4 = K11*I*EoverL;	    // 4EI/L
  // Added by Dimitrios Lignos
  double EIoverL6 = K33*I*EoverL;        // 4EI/L
  
  // determine q = kv + q0
  q(0) = EAoverL*v(0);
  q(1) = EIoverL4*v(1) + EIoverL2*v(2);
//  q(2) = EIoverL2*v(1) + EIoverL4*v(2);
  q(2) = EIoverL2*v(1) + EIoverL6*v(2);
  
  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);
  
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);
  
  return P;
}

int
ModElasticBeam2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

    static Vector data(19);
    
    data(0) = A;
    data(1) = E; 
    data(2) = I; 
    data(3) = rho;
    data(4) = cMass;
    data(5) = this->getTag();
    data(6) = connectedExternalNodes(0);
    data(7) = connectedExternalNodes(1);
    data(8) = theCoordTransf->getClassTag();
    	
    int dbTag = theCoordTransf->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
	theCoordTransf->setDbTag(dbTag);
    }

    data(9) = dbTag;
    data(10) = alpha;
    data(11) = d;

    data(12) = alphaM;
    data(13) = betaK;
    data(14) = betaK0;
    data(15) = betaKc;
 //added by Dimitrios Lignos   
	data(16) = K11;
	data(17) = K33;
	data(18) = K44;
    // Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ModElasticBeam2d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "ModElasticBeam2d::sendSelf -- could not send CoordTransf\n";
      return res;
    }
    
    return res;
}

int
ModElasticBeam2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
	
    static Vector data(19);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ModElasticBeam2d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1); 
    I = data(2); 
    alpha = data(10);
    d = data(11);

    alphaM = data(12);
    betaK = data(13);
    betaK0 = data(14);
    betaKc = data(15);
//added by Dimitrios Lignos
	K11 = data(16);
	K33 = data(17);
	K44 = data(18);

    rho = data(3);
    cMass = (int)data(4);
    this->setTag((int)data(5));
    connectedExternalNodes(0) = (int)data(6);
    connectedExternalNodes(1) = (int)data(7);

    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(8);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ModElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
    
    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ModElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
	
    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(9));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ModElasticBeam2d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }
    
    // Revert the crdtrasf to its last committed state
    theCoordTransf->revertToLastCommit();
    
    return res;
}

void
ModElasticBeam2d::Print(OPS_Stream &s, int flag)
{
  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  }
  
  if (flag == OPS_PRINT_CURRENTSTATE) {
    this->getResistingForce();
    s << "\nModElasticBeam2d: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    double P  = q(0);
    double M1 = q(1);
    double M2 = q(2);
    double L = theCoordTransf->getInitialLength();
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << endln;
  }
  
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"ModElasticBeam2d\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
      s << "\"E\": " << E << ", ";
      s << "\"A\": " << A << ", ";
      s << "\"Iz\": " << I << ", ";
      s << "\"K11\": " << K11 << ", ";
      s << "\"K33\": " << K33 << ", ";
      s << "\"K44\": " << K44 << ", ";
      s << "\"massperlength\": " << rho << ", ";
      s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int
ModElasticBeam2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
    // first determine the end points of the beam based on
    // the display factor (a measure of the distorted image)
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();	

    static Vector v1(3);
    static Vector v2(3);

    if (displayMode >= 0) {
      const Vector &end1Disp = theNodes[0]->getDisp();
      const Vector &end2Disp = theNodes[1]->getDisp();
    
      for (int i = 0; i < 2; i++) {
	v1(i) = end1Crd(i) + end1Disp(i)*fact;
	v2(i) = end2Crd(i) + end2Disp(i)*fact;    
      }
    } else {
      int mode = displayMode  *  -1;
      const Matrix &eigen1 = theNodes[0]->getEigenvectors();
      const Matrix &eigen2 = theNodes[1]->getEigenvectors();
      if (eigen1.noCols() >= mode) {
	for (int i = 0; i < 2; i++) {
	  v1(i) = end1Crd(i) + eigen1(i,mode-1)*fact;
	  v2(i) = end2Crd(i) + eigen2(i,mode-1)*fact;    
	}    
      } else {
	for (int i = 0; i < 2; i++) {
	  v1(i) = end1Crd(i);
	  v2(i) = end2Crd(i);
	}    
      }
    }
    
    return theViewer.drawLine (v1, v2, 1.0, 1.0);
}

Response*
ModElasticBeam2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ModElasticBeam2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

    // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, P);
  
  // local forces
  }    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 3, P);

  // basic forces
  }    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 4, Vector(3));
  }  

  output.endTag(); // ElementOutput

  return theResponse;
}

int
ModElasticBeam2d::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // local forces
    // Axial
    N = q(0);
    P(3) =  N;
    P(0) = -N+p0[0];
    // Moment
    M1 = q(1);
    M2 = q(2);
    P(2) = M1;
    P(5) = M2;
    // Shear
    V = (M1+M2)/L;
    P(1) =  V+p0[1];
    P(4) = -V+p0[2];
    return eleInfo.setVector(P);
    
  case 4: // basic forces
    return eleInfo.setVector(q);

  default:
    return -1;
  }
}

int
ModElasticBeam2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  // A of the beam interior
  else if (strcmp(argv[0],"A") == 0)
    return param.addObject(2, this);
  
  // I of the beam interior
  else if (strcmp(argv[0],"I") == 0)
    return param.addObject(3, this);
   
  // K11 of the beam interior
    
  if (strcmp(argv[0],"K11") == 0) {      
    return param.addObject(4, this);
	  return 4;
    }
  // K33 of the beam interior  
  else if (strcmp(argv[0],"K33") == 0) {      
    return param.addObject(5, this);
	  return 5;
    }
  
  // K44 of the beam interior   
  else if (strcmp(argv[0],"K44") == 0) {      
    return param.addObject(6, this);
	  return 6;
    }

  else   
	  return -1;
}

int
ModElasticBeam2d::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		E = info.theDouble;
		return 0;
	case 2:
		A = info.theDouble;
		return 0;
	case 3:
		I = info.theDouble;
		return 0;
		// Added By Dimitrios Lignos
	case 4:
		this->K11 = info.theDouble;
		return 0;
	case 5:
		this->K33 = info.theDouble;
		return 0;
	case 6:
		this->K44 = info.theDouble;
		return 0;
	default:
		return -1;
	}
}

