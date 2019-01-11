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
                                                                        
// Written: fmk 09/15

#include "ComponentElement2d.h"
#include <ElementalLoad.h>
#include <UniaxialMaterial.h>

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <Renderer.h>

#include <UniaxialMaterial.h>

#include <math.h>
#include <stdlib.h>
#include <SolutionAlgorithm.h>
#include <IncrementalIntegrator.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

Vector ComponentElement2d::P(6);
Matrix ComponentElement2d::K(6,6);

#define ELE_TAG_ComponentElement2d 40

void *
OPS_ComponentElement2d(void)
{
  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr << "Invalid #args,  want: element CompositeElement tag iNode jNode A E I crdTag hinge1 hinge2 \n";
    return 0;
  }
  
  int iData[6];
  double dData[3];  
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids ints" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids double" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING ElasticComponent2d - invalids second set ints" << endln;
    return 0;
  }

  double mass = 0.0;
  int cMass = 0;
  while(OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if(type == "-rho") {
      int numData = 1;
      if(OPS_GetNumRemainingInputArgs() > 0) {
	if(OPS_GetDoubleInput(&numData,&mass) < 0) return 0;
      }
    } else if(type == "-cMass") {
      cMass = 1;
    }
  }

  CrdTransf *theTrans = OPS_getCrdTransf(iData[3]);

  UniaxialMaterial *end1 = OPS_getUniaxialMaterial(iData[4]);
  UniaxialMaterial *end2 = OPS_getUniaxialMaterial(iData[5]);

  // Parsing was successful, allocate the material
  theElement = new ComponentElement2d(iData[0], dData[0], dData[1], dData[2], 
				      iData[1], iData[2], 
				      *theTrans, end1, end2, 
				      mass,cMass);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type ComponentElement2d\n";
    return 0;
  }
  
  return theElement;
}


ComponentElement2d::ComponentElement2d()
  :Element(0,ELE_TAG_ComponentElement2d), 
   A(0.0), E(0.0), I(0.0), rho(0.0), cMass(0),
   Q(6), q(3), connectedExternalNodes(2), theCoordTransf(0)
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

ComponentElement2d::ComponentElement2d(int tag, double a, double e, double i, 
				       int Nd1, int Nd2, CrdTransf &coordTransf,
				       UniaxialMaterial *end1, UniaxialMaterial *end2,
				       double r, int cm)
  :Element(tag,ELE_TAG_ComponentElement2d), 
   A(a), E(e), I(i), rho(r), cMass(cm),
   Q(6), q(3), kb(3,3),
   connectedExternalNodes(2), theCoordTransf(0), end1Hinge(0), end2Hinge(0),
   kTrial(2,2), R(4), uTrial(4), uCommit(4), init(false)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy2d();
  if (!theCoordTransf) {
    opserr << "ComponentElement2d::ComponentElement2d -- failed to get copy of coordinate transformation\n";
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
  
  if (end1 != 0)
    end1Hinge = end1->getCopy();
  if (end2 != 0)
    end2Hinge = end2->getCopy();

  uTrial.Zero();
  uCommit.Zero();
}

ComponentElement2d::~ComponentElement2d()
{
  if (theCoordTransf)
    delete theCoordTransf;

  if (end1Hinge != 0)
    delete end1Hinge;
  
  if (end2Hinge != 0)
    delete end2Hinge;
}

int
ComponentElement2d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ComponentElement2d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ComponentElement2d::getNodePtrs(void) 
{
  return theNodes;
}

int
ComponentElement2d::getNumDOF(void)
{
    return 6;
}

void
ComponentElement2d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ComponentElement2d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "ComponentElement2d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ComponentElement2d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "ComponentElement2d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "ComponentElement2d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ComponentElement2d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ComponentElement2d::setDomain -- Element has zero length\n";
      exit(-1);
    }

    EAoverL  = A*E/L;		// EA/L
    EIoverL2 = 2.0*I*E/L;	// 2EI/L
    EIoverL4 = 4.0*E*I/L;	// 4EI/L
}

int
ComponentElement2d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ComponentElement2d::commitState () - failed in base class";
  }    
  uCommit = uTrial;

  retVal += theCoordTransf->commitState();

  end1Hinge->commitState();
  end2Hinge->commitState();

  return retVal;
}

int
ComponentElement2d::revertToLastCommit()
{
  uTrial = uCommit;

  end1Hinge->revertToLastCommit();
  end2Hinge->revertToLastCommit();

  return theCoordTransf->revertToLastCommit();
}

int
ComponentElement2d::revertToStart()
{
  uCommit.Zero();
  uTrial.Zero();
  init = false;
  end1Hinge->revertToStart();
  end2Hinge->revertToStart();
  return theCoordTransf->revertToStart();
}

int
ComponentElement2d::update(void)
{
  // get previous displacements and the new end delta displacements
  theCoordTransf->update();

  double u1 = uTrial(0);
  double u2 = uTrial(1);
  double u3 = uTrial(2);
  double u4 = uTrial(3);
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  const Vector &dv = theCoordTransf->getBasicIncrDeltaDisp();

  double du1 = dv(1);
  double du4 = dv(2);

  // get hinge forces and tangent
  // NOTE: need tangent used in solution algorithm
  double k1 = 0.;
  double F1 = 0.0;
  if (end1Hinge != 0) {
    F1 = end1Hinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1Hinge->getInitialTangent();      
    }
    else
      k1 = end1Hinge->getTangent();
  }
  //  k1 = end1Hinge->getTangent();

  double k2 = 0.;
  double F2 = 0.0;
  if (end2Hinge != 0) {
    F2 = end2Hinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2Hinge->getInitialTangent();      
    else
      k2 = end2Hinge->getTangent();
  }
  //  k2 = end2Hinge->getTangent();

  // calculate forces for our superelement structure
  double R1 = -F1;
  double R2 =  F1 + EIoverL2*(2*u2 + u3) + q0[1];
  double R3 = -F2 + EIoverL2*(u2 + 2*u3) + q0[2];
  double R4 =  F2;

  // determine change in internal dof, using last K
  // dUi = inv(Kii)*(Pi-Kie*dUe)
  double delta = 1.0/((k1+EIoverL4)*(k2+EIoverL4)-EIoverL2*EIoverL2);
  double du2 = delta*((k2+EIoverL4)*(k1*du1-R2) - EIoverL2*(k2*du4-R3));
  double du3 = delta*(-EIoverL2*(k1*du1-R2) + (k1+EIoverL4)*(k2*du4-R3));

  // update displacements at nodes
  u1 += du1;
  u2 += du2;
  u3 += du3;
  u4 += du4;

  bool converged = false;
  int count = 0;
  int maxCount = 10;
  double tol = 1.0e-10;

  // iterate, at least once, to remove internal node unbalance
  while (converged == false) {

    // set new strain in hinges
    end1Hinge->setTrialStrain(u2-u1);
    end2Hinge->setTrialStrain(u4-u3);

    // obtain new hinge forces and tangents
    k1 = 0.;
    F1 = 0.0;
    if (end1Hinge != 0) {
      F1 = end1Hinge->getStress();
      k1 = end1Hinge->getTangent();    
    }

    k2 = 0.;
    F2 = 0.0;
    if (end2Hinge != 0) {
      F2 = end2Hinge->getStress();
      k2 = end2Hinge->getTangent();
    }
    
    // determine nodal forces
    R1 = -F1;
    R2 =  F1 + EIoverL2 * (2*u2 + u3) + q0[1];
    R3 = -F2 + EIoverL2 * (u2 + 2*u3) + q0[2];
    R4 =  F2;

    // check if converged:
    //    norm resisting forces at internal dof or change in displacement
    //    at these internal dof is less than some tolerance

    if ((sqrt(R2*R2 + R3*R3) > tol) && 
	(sqrt(du2*du2+du3*du3) > tol) &&
	count < maxCount) {

      // if not converged we determine new internal dof displacements
      // note we have not changed du1 or du4 from previous step
      delta = 1.0/((k1+EIoverL4)*(k2+EIoverL4)-EIoverL2*EIoverL2);
      du2 = delta*((k2+EIoverL4)*R2 - EIoverL2*R3);
      du3 = delta*((k1+EIoverL4)*R3 - EIoverL2*R2);

      // unbalance was negative of P so subtract instead of add
      u2 -= du2;
      u3 -= du3;

      count++;

    } else
      converged = true;
  }

  delta = 1.0/((k1+EIoverL4)*(k2+EIoverL4)-EIoverL2*EIoverL2);
  
  // compute new condensed matrix
  kTrial(0,0) = k1 - (delta*k1*k1)*(k2+EIoverL4);
  kTrial(1,1) = k2 - (delta*k2*k2)*(k1+EIoverL4);
  kTrial(0,1) = delta*(k1*k2*EIoverL2);
  kTrial(1,0) = delta*(k1*k2*EIoverL2);

  // compute basic forces, leaving off q0's .. added in getResistingForce
  q(0) = EAoverL*v(0);
  q(1) = R1 + delta*k1*((k2+EIoverL4)*R2 - EIoverL2*R3);
  q(2) = R4 + delta*k2*((k1+EIoverL4)*R3 - EIoverL2*R2);

  // store new displacements
  uTrial(0) = u1;
  uTrial(1) = u2;
  uTrial(2) = u3;
  uTrial(3) = u4;

  return 0;
}

const Vector &
ComponentElement2d::getResistingForce()
{
  // get hinge forces and tangents
  // NOTE: condense out using same tangent as algorithm
  double k1 = 0.;
  double F1 = 0.0;
  if (end1Hinge != 0) {
    F1 = end1Hinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1Hinge->getInitialTangent();      
    }
    else
      k1 = end1Hinge->getTangent();
  }
  //  k1 = end1Hinge->getTangent();

  double k2 = 0.;
  double F2 = 0.0;
  if (end2Hinge != 0) {
    F2 = end2Hinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2Hinge->getInitialTangent();      
    else
      k2 = end2Hinge->getTangent();
  }
  //  k2 = end2Hinge->getTangent();

  double u2 = uTrial(1);
  double u3 = uTrial(2);

  // compute internal forces in our superelement structure
  double R1 = -F1;
  double R2 =  F1 + EIoverL2*(2*u2 + u3) + q0[1];
  double R3 = -F2 + EIoverL2*(u2 + 2*u3) + q0[2];
  double R4 =  F2;

  // condense out internal forces
  double delta = 1.0/((k1+EIoverL4)*(k2+EIoverL4)-EIoverL2*EIoverL2);

  q(0) += q0[0];
  q(1) = R1 + delta*k1*((k2+EIoverL4)*R2 - EIoverL2*R3);
  q(2) = R4 + delta*k2*((k1+EIoverL4)*R3 - EIoverL2*R2);

  // Vector for reactions in basic system
  Vector p0Vec(p0, 3);
  
  // Vector for reactions in basic system
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}


const Matrix &
ComponentElement2d::getTangentStiff(void)
{
  // determine q = kv + q0
  static Vector R(6);  

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  
  kb(0,0) = EAoverL;
  kb(1,1) = kTrial(0,0);
  kb(2,2) = kTrial(1,1);
  kb(1,2) = kTrial(0,1);
  kb(2,1) = kTrial(1,0);

  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
ComponentElement2d::getInitialStiff(void)
{
  double L = theCoordTransf->getInitialLength();

  double k1 = 0.;
  if (end1Hinge != 0) 
    k1 = end1Hinge->getInitialTangent();
  double k2 = 0.;
  if (end2Hinge != 0) 
    k2 = end2Hinge->getInitialTangent();

  double delta = 1.0/((k1+EIoverL4)*(k2+EIoverL4)-EIoverL2*EIoverL2);
  
  // compute new condensed matrix
  static Matrix kb0(3,3);
  kb0(0,0) = EAoverL;
  kb0(1,1) = k1 - delta*(k1*k1*(k2+EIoverL4));
  kb0(2,2) = k2 - delta*(k2*k2*(k1+EIoverL4));
  kb0(1,2) = delta*(k1*k2*EIoverL2);
  kb0(2,1) = delta*(k1*k2*EIoverL2);
  
  return theCoordTransf->getInitialGlobalStiffMatrix(kb0);
}

const Matrix &
ComponentElement2d::getMass(void)
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
ComponentElement2d::zeroLoad(void)
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
ComponentElement2d::addLoad(ElementalLoad *theLoad, double loadFactor)
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
  
  else {
    opserr << "ComponentElement2d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
ComponentElement2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "ComponentElement2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double L = theCoordTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  
  Q(3) -= m * Raccel2(0);
  Q(4) -= m * Raccel2(1);
  
  return 0;
}

const Vector &
ComponentElement2d::getResistingForceIncInertia()
{	
  P = this->getResistingForce();
  
  // subtract external load P = P - Q
  P.addVector(1.0, Q, -1.0);
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
  
  // take advantage of lumped mass matrix
  double L = theCoordTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  P(0) += m * accel1(0);
  P(1) += m * accel1(1);
  
  P(3) += m * accel2(0);
  P(4) += m * accel2(1);
  
  return P;
}



int
ComponentElement2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(16);
  
  data(0) = A;
  data(1) = E; 
  data(2) = I; 
  data(3) = rho;
  //  data(4) = cMass;
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

  // data(10) = alpha;
  //  data(11) = d;
  
  data(12) = alphaM;
  data(13) = betaK;
  data(14) = betaK0;
  data(15) = betaKc;
  
  // Send the data vector
  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
      opserr << "ComponentElement2d::sendSelf -- could not send data Vector\n";
      return res;
  }
  
  // Ask the CoordTransf to send itself
  res += theCoordTransf->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement2d::sendSelf -- could not send CoordTransf\n";
    return res;
  }
  
  return res;
}

int
ComponentElement2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
	
    static Vector data(16);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ComponentElement2d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1); 
    I = data(2); 

    //    alpha = data(10);
    //    d = data(11);

    alphaM = data(12);
    betaK  = data(13);
    betaK0 = data(14);
    betaKc = data(15);

    rho = data(3);
    //    cMass = (int)data(4);
    this->setTag((int)data(5));
    connectedExternalNodes(0) = (int)data(6);
    connectedExternalNodes(1) = (int)data(7);

    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(8);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ComponentElement2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
    
    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ComponentElement2d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
	
    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(9));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement2d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }
    
    return res;
}

void
ComponentElement2d::Print(OPS_Stream &s, int flag)
{
  // to update forces!
  this->getResistingForce();

  if (flag == -1) {
      int eleTag = this->getTag();
      s << "EL_BEAM\t" << eleTag << "\t";
      s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1);
      s << "0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
      this->getResistingForce();
      s << "\nComponentElement2d: " << this->getTag() << endln;
      s << "\tConnected Nodes: " << connectedExternalNodes;
      s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
      s << "\tmass density:  " << rho << endln;
      double P = q(0);
      double M1 = q(1);
      double M2 = q(2);
      double L = theCoordTransf->getInitialLength();
      double V = (M1 + M2) / L;
      s << "\tEnd 1 Forces (P V M): " << -P + p0[0]
          << " " << V + p0[1] << " " << M1 << endln;
      s << "\tEnd 2 Forces (P V M): " << P
          << " " << -V + p0[2] << " " << M2 << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{";
      s << "\"name\": " << this->getTag() << ", ";
      s << "\"type\": \"ComponentElement2d\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
      s << "\"E\": " << E << ", ";
      s << "\"A\": " << A << ", ";
      s << "\"Iz\": " << I << ", ";
      s << "\"massperlength\": " << rho << ", ";
      s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int
ComponentElement2d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  static Vector v1(3);
  static Vector v2(3);
  static Vector vp(3);

  theNodes[0]->getDisplayCrds(v1, fact);
  theNodes[1]->getDisplayCrds(v2, fact);

  float d1 = 0.0;
  float d2 = 0.0;
  float d3 = 0.0;

  int res = 0;

  if (displayMode > 0 && numMode == 0) {

    res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);
    
  } else if (displayMode < 0) {
    
    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);
    
    // add eigenvector values
    int mode = displayMode  *  -1;

    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 2; i++) {
	v1(i) += eigen1(i,mode-1)*fact;
	v2(i) += eigen2(i,mode-1)*fact;    
      }    
    }

    res = theViewer.drawLine (v1, v2, 0.0, 0.0, this->getTag(), 0);
  }

  if (numMode > 0) {
    // calculate q for potential need below
    this->getResistingForce();
    vp = theCoordTransf->getBasicTrialDisp();
  }
  
  for (int i=0; i<numMode; i++) {

    const char *theMode = modes[i];
    if (strcmp(theMode, "axialForce") == 0) {
      d1 = q(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endMoments") == 0) {

      d1 = q(1);
      d2 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localForces") == 0) {
      d1 = q(0);
      d2 = q(1);
      d3 = q(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "axialDeformation") == 0) {
      d1 = vp(0); 
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
      
    } else if (strcmp(theMode, "endRotations") == 0) {

      d1 = vp(1);
      d2 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20.;
      res += theViewer.drawPoint(v1+delta, d1, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d2, this->getTag(), i);

    } else if (strcmp(theMode, "localDeformations") == 0) {
      d1 = vp(0);
      d2 = vp(1);
      d3 = vp(2);
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);

    } else if (strcmp(theMode, "plasticDeformations") == 0) {
      d1 = 0.;
      d2 = 0.;
      d3 = 0.;
      static Vector delta(3); delta = v2-v1; delta/=20;
      res += theViewer.drawPoint(v1+delta, d2, this->getTag(), i);
      res += theViewer.drawPoint(v2-delta, d3, this->getTag(), i);
      res +=theViewer.drawLine(v1, v2, d1, d1, this->getTag(), i);
    }

  }    

  return res;
}

Response*
ComponentElement2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ComponentElement2d");
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
    
  
  } else if (strcmp(argv[0],"hingeDefoAndForce") == 0) {

    output.tag("ResponseType","end1_Defo");
    output.tag("ResponseType","end1_Force");
    output.tag("ResponseType","end2_Defo");
    output.tag("ResponseType","end2_Force");
    
    theResponse = new ElementResponse(this, 5, Vector(4));

  } else if (strcmp(argv[0],"hingeTangent") == 0) {

    output.tag("ResponseType","end1_Tangent");
    output.tag("ResponseType","end1_Tangent");
    
    theResponse = new ElementResponse(this, 6, Vector(2));
  }  


  output.endTag(); // ElementOutput

  return theResponse;
}

int
ComponentElement2d::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();
  this->getResistingForce();
  static Vector vect4(4);
  static Vector vect2(2);

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

  case 5: // basic forces
    vect4.Zero();
    if (end1Hinge != 0) {
      vect4(0) = end1Hinge->getStrain();
      vect4(1) = end1Hinge->getStress();
    }
    if (end1Hinge != 0) {
      vect4(2) = end2Hinge->getStrain();
      vect4(3) = end2Hinge->getStress();
    }
    return eleInfo.setVector(vect4);

  case 6: // basic forces
    if (end1Hinge != 0) {
      vect2(0) = end1Hinge->getTangent();
    }
    if (end1Hinge != 0) {
      vect2(1) = end2Hinge->getTangent();
    }
    return eleInfo.setVector(vect2);



  default:
    return -1;
  }
}

int
ComponentElement2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(1, this);

  // A of the beam interior
  if (strcmp(argv[0],"A") == 0)
    return param.addObject(2, this);
  
  // I of the beam interior
  if (strcmp(argv[0],"I") == 0)
    return param.addObject(3, this);
  
  return -1;
}

int
ComponentElement2d::updateParameter (int parameterID, Information &info)
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
	default:
		return -1;
	}
}

