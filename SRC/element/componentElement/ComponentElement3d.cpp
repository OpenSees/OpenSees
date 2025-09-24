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
                                                                        
// Written: fmk 09/15, MHS 11/22 for 3D

#include "ComponentElement3d.h"
#include <ElementalLoad.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>

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

Vector ComponentElement3d::P(12);
Matrix ComponentElement3d::K(12,12);

void *
OPS_ComponentElement3d(void)
{
  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 3) {
    opserr << "Invalid #args,  want: element componentElement tag iNode jNode A E G J Iy Iz crdTag hinge1z hinge2z hinge1y hinge2y \n";
    return 0;
  }
  
  int iData[8];
  double dData[6];  
  int numData = 3;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING componentElement - invalid ints" << endln;
    return 0;
  }

  numData = 6;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING componentElement - invalid double" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[3]) != 0) {
    opserr << "WARNING componentElement - invalid transformation tag" << endln;
    return 0;
  }

  bool useK = false;
  double k[4];

  std::string flag = OPS_GetString();
  if (flag == "-stiffness" || flag == "-k") {
    numData = 4;
    if (OPS_GetDoubleInput(&numData, k) != 0) {
      opserr << "WARNING componentElement - invalid stiffness values" << endln;
      return 0;
    }
    useK = true;
  }
  else {
    OPS_ResetCurrentInputArg(-1);
    numData = 4;
    if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
      opserr << "WARNING componentElement - invalid second material tag" << endln;
      return 0;
    }    
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

  if (useK) {
    theElement = new ComponentElement3d(iData[0], dData[0], dData[1], dData[5],
					dData[4], dData[2], dData[3],
					iData[1], iData[2], 
					*theTrans, k[0], k[1], k[2], k[3],
					mass,cMass);
  }
  else {
    UniaxialMaterial *end1z = OPS_getUniaxialMaterial(iData[4]);
    UniaxialMaterial *end2z = OPS_getUniaxialMaterial(iData[5]);
    UniaxialMaterial *end1y = OPS_getUniaxialMaterial(iData[6]);
    UniaxialMaterial *end2y = OPS_getUniaxialMaterial(iData[7]);  
    
    // Parsing was successful, allocate the material
    theElement = new ComponentElement3d(iData[0], dData[0], dData[1], dData[5],
					dData[4], dData[2], dData[3],
					iData[1], iData[2], 
					*theTrans, end1z, end2z, end1y, end2y, 
					mass,cMass);
  }
  
  if (theElement == 0) {
    opserr << "WARNING could not create element of type componentElement\n";
    return 0;
  }
  
  return theElement;
}


ComponentElement3d::ComponentElement3d()
  :Element(0,ELE_TAG_ComponentElement3d), 
   A(0.0), E(0.0), Iz(0.0), Iy(0.0), G(0.0), J(0.0), rho(0.0), cMass(0),
   Q(12), q(6), connectedExternalNodes(2), theCoordTransf(0),
   end1zHinge(0), end2zHinge(0), end1yHinge(0), end2yHinge(0),
   kzTrial(2,2), uzTrial(4), uzCommit(4),
   kyTrial(2,2), uyTrial(4), uyCommit(4), kb(6,6), init(false)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;  

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;        
}

ComponentElement3d::ComponentElement3d(int tag, double a, double e, double iz,
				       double iy, double g, double j,
				       int Nd1, int Nd2, CrdTransf &coordTransf,
				       UniaxialMaterial *end1z, UniaxialMaterial *end2z,
				       UniaxialMaterial *end1y, UniaxialMaterial *end2y,
				       double r, int cm)
  :Element(tag,ELE_TAG_ComponentElement3d), 
   A(a), E(e), Iz(iz), Iy(iy), G(g), J(j), rho(r), cMass(cm),
   Q(12), q(6), 
   connectedExternalNodes(2), theCoordTransf(0),
   end1zHinge(0), end2zHinge(0), end1yHinge(0), end2yHinge(0),
   kzTrial(2,2), uzTrial(4), uzCommit(4),
   kyTrial(2,2), uyTrial(4), uyCommit(4), kb(6,6), init(false)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy3d();
  if (!theCoordTransf) {
    opserr << "ComponentElement3d::ComponentElement3d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;    

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;    

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
  
  if (end1z != 0)
    end1zHinge = end1z->getCopy();
  if (end2z != 0)
    end2zHinge = end2z->getCopy();
  if (end1y != 0)
    end1yHinge = end1y->getCopy();
  if (end2y != 0)
    end2yHinge = end2y->getCopy();  

  uzTrial.Zero();
  uzCommit.Zero();
  uyTrial.Zero();
  uyCommit.Zero();  
}

ComponentElement3d::ComponentElement3d(int tag, double a, double e, double iz,
				       double iy, double g, double j,
				       int Nd1, int Nd2, CrdTransf &coordTransf,
				       double kzI, double kzJ, double kyI, double kyJ,
				       double r, int cm)
  :Element(tag,ELE_TAG_ComponentElement3d), 
   A(a), E(e), Iz(iz), Iy(iy), G(g), J(j), rho(r), cMass(cm),
   Q(12), q(6), 
   connectedExternalNodes(2), theCoordTransf(0),
   end1zHinge(0), end2zHinge(0), end1yHinge(0), end2yHinge(0),
   kzTrial(2,2), uzTrial(4), uzCommit(4),
   kyTrial(2,2), uyTrial(4), uyCommit(4), kb(6,6), init(false)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy3d();
  if (!theCoordTransf) {
    opserr << "ComponentElement3d::ComponentElement3d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;    

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;    

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
  
  if (kzI > 0.0)
    end1zHinge = new ElasticMaterial(0, kzI);
  if (kzJ > 0.0)
    end2zHinge = new ElasticMaterial(0, kzJ);
  if (kyI > 0.0)
    end1yHinge = new ElasticMaterial(0, kyI);
  if (kyJ > 0.0)
    end2yHinge = new ElasticMaterial(0, kyJ);  

  uzTrial.Zero();
  uzCommit.Zero();
  uyTrial.Zero();
  uyCommit.Zero();    
}

ComponentElement3d::~ComponentElement3d()
{
  if (theCoordTransf)
    delete theCoordTransf;

  if (end1zHinge != 0)
    delete end1zHinge;
  
  if (end2zHinge != 0)
    delete end2zHinge;

  if (end1yHinge != 0)
    delete end1yHinge;
  
  if (end2yHinge != 0)
    delete end2yHinge;  
}

int
ComponentElement3d::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
ComponentElement3d::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
ComponentElement3d::getNodePtrs(void) 
{
  return theNodes;
}

int
ComponentElement3d::getNumDOF(void)
{
    return 12;
}

void
ComponentElement3d::setDomain(Domain *theDomain)
{
  if (theDomain == 0) {
    opserr << "ComponentElement3d::setDomain -- Domain is null\n";
    exit(-1);
  }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == 0) {
      opserr << "ComponentElement3d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      exit(-1);
    }
			      
    if (theNodes[1] == 0) {
      opserr << "ComponentElement3d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      exit(-1);
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 6) {
      opserr << "ComponentElement3d::setDomain -- Node 1: " << connectedExternalNodes(0) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 6) {
      opserr << "ComponentElement3d::setDomain -- Node 2: " << connectedExternalNodes(1) 
	     << " has incorrect number of DOF\n";
      exit(-1);
    }
	
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
	opserr << "ComponentElement3d::setDomain -- Error initializing coordinate transformation\n";
	exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ComponentElement3d::setDomain -- Element has zero length\n";
      exit(-1);
    }

    EAoverL  = A*E/L;		// EA/L
    EIzoverL2 = 2.0*Iz*E/L;	// 2EI/L
    EIzoverL4 = 2.0*EIzoverL2;	// 4EI/L
    EIyoverL2 = 2.0*Iy*E/L;	// 2EI/L
    EIyoverL4 = 2.0*EIyoverL2;	// 4EI/L
    GJoverL = G*J/L;
}

int
ComponentElement3d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ComponentElement3d::commitState () - failed in base class";
  }    
  uzCommit = uzTrial;
  uyCommit = uyTrial;

  retVal += theCoordTransf->commitState();

  if (end1zHinge != 0)
    end1zHinge->commitState();
  if (end2zHinge != 0)
    end2zHinge->commitState();
  if (end1yHinge != 0)
    end1yHinge->commitState();
  if (end2yHinge != 0)
    end2yHinge->commitState();  

  return retVal;
}

int
ComponentElement3d::revertToLastCommit()
{
  uzTrial = uzCommit;
  uyTrial = uyCommit;

  if (end1zHinge != 0)  
    end1zHinge->revertToLastCommit();
  if (end2zHinge != 0)
    end2zHinge->revertToLastCommit();
  if (end1yHinge != 0)
    end1yHinge->revertToLastCommit();
  if (end2yHinge != 0)
    end2yHinge->revertToLastCommit();  

  return theCoordTransf->revertToLastCommit();
}

int
ComponentElement3d::revertToStart()
{
  uzCommit.Zero();
  uzTrial.Zero();
  uyCommit.Zero();
  uyTrial.Zero();  
  init = false;

  if (end1zHinge != 0)  
    end1zHinge->revertToStart();
  if (end2zHinge != 0)
    end2zHinge->revertToStart();
  if (end1yHinge != 0)
    end1yHinge->revertToStart();
  if (end2yHinge != 0)
    end2yHinge->revertToStart();
  
  return theCoordTransf->revertToStart();
}

int
ComponentElement3d::update(void)
{
  // get previous displacements and the new end delta displacements
  theCoordTransf->update();

  double u1 = uzTrial(0);
  double u2 = uzTrial(1);
  double u3 = uzTrial(2);
  double u4 = uzTrial(3);
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  const Vector &dv = theCoordTransf->getBasicIncrDeltaDisp();

  q(0) = EAoverL*v(0);
  q(5) = GJoverL*v(5);  

  
  double du1 = dv(1);
  double du4 = dv(2);

  // get hinge forces and tangent
  // NOTE: need tangent used in solution algorithm
  double k1 = 0.;
  double F1 = 0.0;
  if (end1zHinge != 0) {
    F1 = end1zHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1zHinge->getInitialTangent();      
    }
    else
      k1 = end1zHinge->getTangent();
  }
  //  k1 = end1zHinge->getTangent();

  double k2 = 0.;
  double F2 = 0.0;
  if (end2zHinge != 0) {
    F2 = end2zHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2zHinge->getInitialTangent();      
    else
      k2 = end2zHinge->getTangent();
  }
  //  k2 = end2zHinge->getTangent();

  // calculate forces for our superelement structure
  double R1 = -F1;
  double R2 =  F1 + EIzoverL2*(2*u2 + u3) + q0[1];
  double R3 = -F2 + EIzoverL2*(u2 + 2*u3) + q0[2];
  double R4 =  F2;

  // determine change in internal dof, using last K
  // dUi = inv(Kii)*(Pi-Kie*dUe)
  double delta = 1.0/((k1+EIzoverL4)*(k2+EIzoverL4)-EIzoverL2*EIzoverL2);
  double du2 = delta*((k2+EIzoverL4)*(k1*du1-R2) - EIzoverL2*(k2*du4-R3));
  double du3 = delta*(-EIzoverL2*(k1*du1-R2) + (k1+EIzoverL4)*(k2*du4-R3));

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
    if (end1zHinge != 0)
      end1zHinge->setTrialStrain(u2-u1);
    if (end2zHinge != 0)
      end2zHinge->setTrialStrain(u4-u3);

    // obtain new hinge forces and tangents
    k1 = 0.;
    F1 = 0.0;
    if (end1zHinge != 0) {
      F1 = end1zHinge->getStress();
      k1 = end1zHinge->getTangent();    
    }

    k2 = 0.;
    F2 = 0.0;
    if (end2zHinge != 0) {
      F2 = end2zHinge->getStress();
      k2 = end2zHinge->getTangent();
    }
    
    // determine nodal forces
    R1 = -F1;
    R2 =  F1 + EIzoverL2 * (2*u2 + u3) + q0[1];
    R3 = -F2 + EIzoverL2 * (u2 + 2*u3) + q0[2];
    R4 =  F2;

    // check if converged:
    //    norm resisting forces at internal dof or change in displacement
    //    at these internal dof is less than some tolerance

    if ((sqrt(R2*R2 + R3*R3) > tol) && 
	(sqrt(du2*du2+du3*du3) > tol) &&
	count < maxCount) {

      // if not converged we determine new internal dof displacements
      // note we have not changed du1 or du4 from previous step
      delta = 1.0/((k1+EIzoverL4)*(k2+EIzoverL4)-EIzoverL2*EIzoverL2);
      du2 = delta*((k2+EIzoverL4)*R2 - EIzoverL2*R3);
      du3 = delta*((k1+EIzoverL4)*R3 - EIzoverL2*R2);

      // unbalance was negative of P so subtract instead of add
      u2 -= du2;
      u3 -= du3;

      count++;

    } else
      converged = true;
  }

  delta = 1.0/((k1+EIzoverL4)*(k2+EIzoverL4)-EIzoverL2*EIzoverL2);
  
  // compute new condensed matrix
  kzTrial(0,0) = k1 - (delta*k1*k1)*(k2+EIzoverL4);
  kzTrial(1,1) = k2 - (delta*k2*k2)*(k1+EIzoverL4);
  kzTrial(0,1) = delta*(k1*k2*EIzoverL2);
  kzTrial(1,0) = delta*(k1*k2*EIzoverL2);

  // compute basic forces, leaving off q0's .. added in getResistingForce
  q(1) = R1 + delta*k1*((k2+EIzoverL4)*R2 - EIzoverL2*R3);
  q(2) = R4 + delta*k2*((k1+EIzoverL4)*R3 - EIzoverL2*R2);

  // store new displacements
  uzTrial(0) = u1;
  uzTrial(1) = u2;
  uzTrial(2) = u3;
  uzTrial(3) = u4;






  u1 = uyTrial(0);
  u2 = uyTrial(1);
  u3 = uyTrial(2);
  u4 = uyTrial(3);
  
  //const Vector &v = theCoordTransf->getBasicTrialDisp();
  //const Vector &dv = theCoordTransf->getBasicIncrDeltaDisp();

  du1 = dv(3);
  du4 = dv(4);

  // get hinge forces and tangent
  // NOTE: need tangent used in solution algorithm
  k1 = 0.;
  F1 = 0.0;
  if (end1yHinge != 0) {
    F1 = end1yHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1yHinge->getInitialTangent();      
    }
    else
      k1 = end1yHinge->getTangent();
  }
  //  k1 = end1yHinge->getTangent();

  k2 = 0.;
  F2 = 0.0;
  if (end2yHinge != 0) {
    F2 = end2yHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2yHinge->getInitialTangent();      
    else
      k2 = end2yHinge->getTangent();
  }
  //  k2 = end2yHinge->getTangent();

  // calculate forces for our superelement structure
  R1 = -F1;
  R2 =  F1 + EIyoverL2*(2*u2 + u3) + q0[3];
  R3 = -F2 + EIyoverL2*(u2 + 2*u3) + q0[4];
  R4 =  F2;

  // determine change in internal dof, using last K
  // dUi = inv(Kii)*(Pi-Kie*dUe)
  delta = 1.0/((k1+EIyoverL4)*(k2+EIyoverL4)-EIyoverL2*EIyoverL2);
  du2 = delta*((k2+EIyoverL4)*(k1*du1-R2) - EIyoverL2*(k2*du4-R3));
  du3 = delta*(-EIyoverL2*(k1*du1-R2) + (k1+EIyoverL4)*(k2*du4-R3));

  // update displacements at nodes
  u1 += du1;
  u2 += du2;
  u3 += du3;
  u4 += du4;

  converged = false;
  count = 0;
  maxCount = 10;
  tol = 1.0e-10;

  // iterate, at least once, to remove internal node unbalance
  while (converged == false) {

    // set new strain in hinges
    if (end1yHinge != 0)
      end1yHinge->setTrialStrain(u2-u1);
    if (end2yHinge != 0)
      end2yHinge->setTrialStrain(u4-u3);

    // obtain new hinge forces and tangents
    k1 = 0.;
    F1 = 0.0;
    if (end1yHinge != 0) {
      F1 = end1yHinge->getStress();
      k1 = end1yHinge->getTangent();    
    }

    k2 = 0.;
    F2 = 0.0;
    if (end2yHinge != 0) {
      F2 = end2yHinge->getStress();
      k2 = end2yHinge->getTangent();
    }
    
    // determine nodal forces
    R1 = -F1;
    R2 =  F1 + EIyoverL2 * (2*u2 + u3) + q0[3];
    R3 = -F2 + EIyoverL2 * (u2 + 2*u3) + q0[4];
    R4 =  F2;

    // check if converged:
    //    norm resisting forces at internal dof or change in displacement
    //    at these internal dof is less than some tolerance

    if ((sqrt(R2*R2 + R3*R3) > tol) && 
	(sqrt(du2*du2+du3*du3) > tol) &&
	count < maxCount) {

      // if not converged we determine new internal dof displacements
      // note we have not changed du1 or du4 from previous step
      delta = 1.0/((k1+EIyoverL4)*(k2+EIyoverL4)-EIyoverL2*EIyoverL2);
      du2 = delta*((k2+EIyoverL4)*R2 - EIyoverL2*R3);
      du3 = delta*((k1+EIyoverL4)*R3 - EIyoverL2*R2);

      // unbalance was negative of P so subtract instead of add
      u2 -= du2;
      u3 -= du3;

      count++;

    } else
      converged = true;
  }

  delta = 1.0/((k1+EIyoverL4)*(k2+EIyoverL4)-EIyoverL2*EIyoverL2);
  
  // compute new condensed matrix
  kyTrial(0,0) = k1 - (delta*k1*k1)*(k2+EIyoverL4);
  kyTrial(1,1) = k2 - (delta*k2*k2)*(k1+EIyoverL4);
  kyTrial(0,1) = delta*(k1*k2*EIyoverL2);
  kyTrial(1,0) = delta*(k1*k2*EIyoverL2);

  // compute basic forces, leaving off q0's .. added in getResistingForce
  q(3) = R1 + delta*k1*((k2+EIyoverL4)*R2 - EIyoverL2*R3);
  q(4) = R4 + delta*k2*((k1+EIyoverL4)*R3 - EIyoverL2*R2);

  // store new displacements
  uyTrial(0) = u1;
  uyTrial(1) = u2;
  uyTrial(2) = u3;
  uyTrial(3) = u4;



  
  return 0;
}

const Vector &
ComponentElement3d::getResistingForce()
{
  // get hinge forces and tangents
  // NOTE: condense out using same tangent as algorithm
  double k1 = 0.;
  double F1 = 0.0;
  if (end1zHinge != 0) {
    F1 = end1zHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1zHinge->getInitialTangent();      
    }
    else
      k1 = end1zHinge->getTangent();
  }
  //  k1 = end1zHinge->getTangent();

  double k2 = 0.;
  double F2 = 0.0;
  if (end2zHinge != 0) {
    F2 = end2zHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2zHinge->getInitialTangent();      
    else
      k2 = end2zHinge->getTangent();
  }
  //  k2 = end2zHinge->getTangent();

  double u2 = uzTrial(1);
  double u3 = uzTrial(2);

  // compute internal forces in our superelement structure
  double R1 = -F1;
  double R2 =  F1 + EIzoverL2*(2*u2 + u3) + q0[1];
  double R3 = -F2 + EIzoverL2*(u2 + 2*u3) + q0[2];
  double R4 =  F2;

  // condense out internal forces
  double delta = 1.0/((k1+EIzoverL4)*(k2+EIzoverL4)-EIzoverL2*EIzoverL2);

  q(0) += q0[0];
  q(1) = R1 + delta*k1*((k2+EIzoverL4)*R2 - EIzoverL2*R3);
  q(2) = R4 + delta*k2*((k1+EIzoverL4)*R3 - EIzoverL2*R2);


  // get hinge forces and tangents
  // NOTE: condense out using same tangent as algorithm
  k1 = 0.;
  F1 = 0.0;
  if (end1yHinge != 0) {
    F1 = end1yHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT) {
      k1 = end1yHinge->getInitialTangent();      
    }
    else
      k1 = end1yHinge->getTangent();
  }
  //  k1 = end1yHinge->getTangent();

  k2 = 0.;
  F2 = 0.0;
  if (end2yHinge != 0) {
    F2 = end2yHinge->getStress();
    if (SOLUTION_ALGORITHM_tangentFlag == INITIAL_TANGENT)
      k2 = end2yHinge->getInitialTangent();      
    else
      k2 = end2yHinge->getTangent();
  }
  //  k2 = end2yHinge->getTangent();

  u2 = uyTrial(1);
  u3 = uyTrial(2);

  // compute internal forces in our superelement structure
  R1 = -F1;
  R2 =  F1 + EIyoverL2*(2*u2 + u3) + q0[3];
  R3 = -F2 + EIyoverL2*(u2 + 2*u3) + q0[4];
  R4 =  F2;

  // condense out internal forces
  delta = 1.0/((k1+EIyoverL4)*(k2+EIyoverL4)-EIyoverL2*EIyoverL2);

  q(3) = R1 + delta*k1*((k2+EIyoverL4)*R2 - EIyoverL2*R3);
  q(4) = R4 + delta*k2*((k1+EIyoverL4)*R3 - EIyoverL2*R2);

  
  // Vector for reactions in basic system
  Vector p0Vec(p0, 5);
  
  // Vector for reactions in basic system
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}


const Matrix &
ComponentElement3d::getTangentStiff(void)
{
  // determine q = kv + q0
  static Vector R(12);  

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];
  
  kb(0,0) = EAoverL;
  kb(5,5) = GJoverL;
  
  kb(1,1) = kzTrial(0,0);
  kb(2,2) = kzTrial(1,1);
  kb(1,2) = kzTrial(0,1);
  kb(2,1) = kzTrial(1,0);

  kb(3,3) = kyTrial(0,0);
  kb(4,4) = kyTrial(1,1);
  kb(3,4) = kyTrial(0,1);
  kb(4,3) = kyTrial(1,0);  

  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
ComponentElement3d::getInitialStiff(void)
{
  double k1 = 0.;
  if (end1zHinge != 0) 
    k1 = end1zHinge->getInitialTangent();
  double k2 = 0.;
  if (end2zHinge != 0) 
    k2 = end2zHinge->getInitialTangent();

  double delta = 1.0/((k1+EIzoverL4)*(k2+EIzoverL4)-EIzoverL2*EIzoverL2);
  
  // compute new condensed matrix
  static Matrix kb0(6,6);
  kb0(0,0) = EAoverL;
  kb0(5,5) = GJoverL;
  
  kb0(1,1) = k1 - delta*(k1*k1*(k2+EIzoverL4));
  kb0(2,2) = k2 - delta*(k2*k2*(k1+EIzoverL4));
  kb0(1,2) = delta*(k1*k2*EIzoverL2);
  kb0(2,1) = delta*(k1*k2*EIzoverL2);


  k1 = 0.;
  if (end1yHinge != 0) 
    k1 = end1yHinge->getInitialTangent();
  k2 = 0.;
  if (end2yHinge != 0) 
    k2 = end2yHinge->getInitialTangent();

  delta = 1.0/((k1+EIyoverL4)*(k2+EIyoverL4)-EIyoverL2*EIyoverL2);
  
  kb0(3,3) = k1 - delta*(k1*k1*(k2+EIyoverL4));
  kb0(4,4) = k2 - delta*(k2*k2*(k1+EIyoverL4));
  kb0(3,4) = delta*(k1*k2*EIyoverL2);
  kb0(4,3) = delta*(k1*k2*EIyoverL2);

  
  return theCoordTransf->getInitialGlobalStiffMatrix(kb0);
}

const Matrix &
ComponentElement3d::getMass(void)
{ 
  K.Zero();
  
  if (rho > 0.0)  {
    // get initial element length
    double L = theCoordTransf->getInitialLength();
    
        if (cMass == 0)  {

            // lumped mass matrix
            double m = 0.5*rho*L;
            K(0,0) = K(1,1) = K(2,2) = K(6,6) = K(7,7) = K(8,8) = m;

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
ComponentElement3d::zeroLoad(void)
{
  Q.Zero();

  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;
  q0[3] = 0.0;
  q0[4] = 0.0;  

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;
  p0[3] = 0.0;
  p0[4] = 0.0;
  
  return;
}

int 
ComponentElement3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0)*loadFactor;  // Transverse
    double wz = data(1)*loadFactor;  // Transverse
    double wx = data(2)*loadFactor;  // Axial (+ve from node I to J)
    
    double Vy = 0.5*wy*L;
    double Mz = Vy*L/6.0; // wy*L*L/12
    double Vz = 0.5*wz*L;
    double My = Vz*L/6.0; // wz*L*L/12
    double P = wx*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    if (end1zHinge != 0 && end2zHinge != 0) {
      q0[1] -= Mz;
      q0[2] += Mz;
    }
    if (end1zHinge == 0 && end2zHinge != 0) {
      q0[2] += wy*L*L/8;
    }
    if (end1zHinge != 0 && end2zHinge == 0) {
      q0[1] -= wy*L*L/8;
    }
    
    if (end1yHinge != 0 && end2yHinge != 0) {
      q0[3] += My;
      q0[4] -= My;
    }
    if (end1yHinge == 0 && end2yHinge != 0) {
      q0[4] -= wz*L*L/8;
    }
    if (end1yHinge != 0 && end2yHinge == 0) {
      q0[3] += wz*L*L/8;
    }    
  }

  else {
    opserr << "ComponentElement3d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
ComponentElement3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
	
  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "ComponentElement3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double L = theCoordTransf->getInitialLength();
  double m = 0.5*rho*L;
  
  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  Q(2) -= m * Raccel1(2);  
  
  Q(6) -= m * Raccel2(0);
  Q(7) -= m * Raccel2(1);
  Q(8) -= m * Raccel2(2);  
  
  return 0;
}

const Vector &
ComponentElement3d::getResistingForceIncInertia()
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
  P(2) += m * accel1(2);  
  
  P(6) += m * accel2(0);
  P(7) += m * accel2(1);
  P(8) += m * accel2(2);  
  
  return P;
}



int
ComponentElement3d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(19+8+4+4); // 35
  
  data(0) = A;
  data(1) = E; 
  data(2) = Iz;
  data(16) = Iy;
  data(17) = G;
  data(18) = J;
  data(3) = rho;
  data(4) = cMass;
  data(5) = this->getTag();
  data(6) = connectedExternalNodes(0);
  data(7) = connectedExternalNodes(1);

  for (int i = 0; i < 4; i++) {
    data(27+i) = uzCommit(i);
    data(31+i) = uyCommit(i);
  }
  
  data(8) = theCoordTransf->getClassTag();
  int dbTag = theCoordTransf->getDbTag();
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      theCoordTransf->setDbTag(dbTag);
  }
  data(9) = dbTag;

  data(19) = end1zHinge->getClassTag();
  dbTag = end1zHinge->getDbTag();
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      end1zHinge->setDbTag(dbTag);
  }
  data(20) = dbTag;

  data(21) = end2zHinge->getClassTag();
  dbTag = end2zHinge->getDbTag();
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      end2zHinge->setDbTag(dbTag);
  }
  data(22) = dbTag;

  data(23) = end1yHinge->getClassTag();
  dbTag = end1yHinge->getDbTag();
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      end1yHinge->setDbTag(dbTag);
  }
  data(24) = dbTag;

  data(25) = end2yHinge->getClassTag();
  dbTag = end2yHinge->getDbTag();
  if (dbTag == 0) {
    dbTag = theChannel.getDbTag();
    if (dbTag != 0)
      end2yHinge->setDbTag(dbTag);
  }
  data(26) = dbTag;    

  
  // data(10) = alpha;
  //  data(11) = d;
  
  data(12) = alphaM;
  data(13) = betaK;
  data(14) = betaK0;
  data(15) = betaKc;
  
  // Send the data vector
  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
      opserr << "ComponentElement3d::sendSelf -- could not send data Vector\n";
      return res;
  }
  
  // Ask the CoordTransf to send itself
  res += theCoordTransf->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement3d::sendSelf -- could not send CoordTransf\n";
    return res;
  }

  // Ask hinge 1 to send itself
  res += end1zHinge->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement3d::sendSelf -- could not send hinge 1z\n";
    return res;
  }

  // Ask hinge 2 to send itself
  res += end2zHinge->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement3d::sendSelf -- could not send hinge 2z\n";
    return res;
  }

  // Ask hinge 1 to send itself
  res += end1yHinge->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement3d::sendSelf -- could not send hinge 1y\n";
    return res;
  }

  // Ask hinge 2 to send itself
  res += end2yHinge->sendSelf(cTag, theChannel);
  if (res < 0) {
    opserr << "ComponentElement3d::sendSelf -- could not send hinge 2y\n";
    return res;
  }  
  
  return res;
}

int
ComponentElement3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
	
    static Vector data(19+8+4+4); // 35

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1); 
    Iz = data(2);
    Iy = data(16);
    G = data(17);
    J = data(18);         

    //    alpha = data(10);
    //    d = data(11);

    alphaM = data(12);
    betaK  = data(13);
    betaK0 = data(14);
    betaKc = data(15);

    rho = data(3);
    cMass = (int)data(4);
    this->setTag((int)data(5));
    connectedExternalNodes(0) = (int)data(6);
    connectedExternalNodes(1) = (int)data(7);

    for (int i = 0; i < 4; i++) {
      uzCommit(i) = data(27+i);
      uyCommit(i) = data(31+i);
    }
    
    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(8);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
    
    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get a CrdTransf2d\n";
	exit(-1);
      }
    }
	
    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(9));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }



    // Check if hinge 1 is null; if so, get a new one
    int hinge1zTag = (int)data(19);
    if (end1zHinge == 0) {
      end1zHinge = theBroker.getNewUniaxialMaterial(hinge1zTag);
      if (end1zHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 1z UniaxialMaterial" << endln;
	exit(-1);
      }
    }
    
    // Check that hinge 1 is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (end1zHinge->getClassTag() != hinge1zTag) {
      delete end1zHinge;
      end1zHinge = theBroker.getNewUniaxialMaterial(hinge1zTag);
      if (end1zHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 1z UniaxialMaterial" << endln;
	exit(-1);
      }
    }
	
    // Now, receive hinge 1
    end1zHinge->setDbTag((int)data(20));
    res += end1zHinge->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive hinge 1z" << endln;
      return res;
    }


    // Check if hinge 2 is null; if so, get a new one
    int hinge2zTag = (int)data(21);
    if (end2zHinge == 0) {
      end2zHinge = theBroker.getNewUniaxialMaterial(hinge2zTag);
      if (end2zHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 2z UniaxialMaterial" << endln;
	exit(-1);
      }
    }
    
    // Check that hinge 2 is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (end2zHinge->getClassTag() != hinge2zTag) {
      delete end2zHinge;
      end2zHinge = theBroker.getNewUniaxialMaterial(hinge2zTag);
      if (end2zHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 2z UniaxialMaterial" << endln;
	exit(-1);
      }
    }
	
    // Now, receive hinge 2
    end2zHinge->setDbTag((int)data(22));
    res += end2zHinge->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive hinge 2z" << endln;
      return res;
    }




    // Check if hinge 1 is null; if so, get a new one
    int hinge1yTag = (int)data(23);
    if (end1yHinge == 0) {
      end1yHinge = theBroker.getNewUniaxialMaterial(hinge1yTag);
      if (end1yHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 1y UniaxialMaterial" << endln;
	exit(-1);
      }
    }
    
    // Check that hinge 1 is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (end1yHinge->getClassTag() != hinge1yTag) {
      delete end1yHinge;
      end1yHinge = theBroker.getNewUniaxialMaterial(hinge1yTag);
      if (end1yHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 1y UniaxialMaterial" << endln;
	exit(-1);
      }
    }
	
    // Now, receive hinge 1
    end1yHinge->setDbTag((int)data(24));
    res += end1yHinge->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive hinge y" << endln;
      return res;
    }


    // Check if hinge 2 is null; if so, get a new one
    int hinge2yTag = (int)data(25);
    if (end2yHinge == 0) {
      end2yHinge = theBroker.getNewUniaxialMaterial(hinge2yTag);
      if (end2yHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 2y UniaxialMaterial" << endln;
	exit(-1);
      }
    }
    
    // Check that hinge 2 is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (end2yHinge->getClassTag() != hinge2yTag) {
      delete end2yHinge;
      end2yHinge = theBroker.getNewUniaxialMaterial(hinge2yTag);
      if (end2yHinge == 0) {
	opserr << "ComponentElement3d::recvSelf -- could not get hinge 2y UniaxialMaterial" << endln;
	exit(-1);
      }
    }
	
    // Now, receive hinge 2
    end2yHinge->setDbTag((int)data(26));
    res += end2yHinge->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ComponentElement3d::recvSelf -- could not receive hinge 2y" << endln;
      return res;
    }

    this->revertToLastCommit();
    
    return res;
}

void
ComponentElement3d::Print(OPS_Stream &s, int flag)
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
      s << "\nComponentElement3d: " << this->getTag() << endln;
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
      s << "\"type\": \"ComponentElement3d\", ";
      s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
      s << "\"E\": " << E << ", ";
      s << "\"A\": " << A << ", ";
      s << "\"Iz\": " << Iz << ", ";
      s << "\"massperlength\": " << rho << ", ";
      s << "\"materials\": [" ;
      if (end1zHinge) 
        s << "\"" << end1zHinge->getTag() << "\", ";
      else
        s << "null, ";
      if (end2zHinge) 
        s << "\"" << end2zHinge->getTag() << "\"], ";
      else
        s << "null], ";
      s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}

int
ComponentElement3d::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{
  static Vector v1(3);
  static Vector v2(3);
  static Vector vp(3);

  theNodes[0]->getDisplayCrds(v1, fact, displayMode);
  theNodes[1]->getDisplayCrds(v2, fact, displayMode);

  float d1 = 0.0;
  float d2 = 0.0;
  float d3 = 0.0;

  int res = 0;

  if (displayMode > 0 && numMode == 0)
      res += theViewer.drawLine(v1, v2, d1, d1, this->getTag(), 0);
  else if (displayMode < 0)
      return theViewer.drawLine(v1, v2, 0.0, 0.0, this->getTag(), 0);

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
ComponentElement3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ComponentElement3d");
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
  // basic stiffness
  else if (strcmp(argv[0],"basicStiffness") == 0) {
    output.tag("ResponseType","N");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Mz_2");
    output.tag("ResponseType","My_1");
    output.tag("ResponseType","My_2");
    output.tag("ResponseType","T");
    
    theResponse = new ElementResponse(this, 19, Matrix(6,6));
  }
  else if (strcmp(argv[0],"hingeDefoAndForce") == 0) {

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

  if (theResponse == 0)
    theResponse = theCoordTransf->setResponse(argv, argc, output);
  
  return theResponse;
}

int
ComponentElement3d::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();
  this->getResistingForce();
  static Vector vect4(4);
  static Vector vect2(2);
  static Matrix kb(6,6);
  
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
    if (end1zHinge != 0) {
      vect4(0) = end1zHinge->getStrain();
      vect4(1) = end1zHinge->getStress();
    }
    if (end2zHinge != 0) {
      vect4(2) = end2zHinge->getStrain();
      vect4(3) = end2zHinge->getStress();
    }
    return eleInfo.setVector(vect4);

  case 6: // basic forces
    if (end1zHinge != 0) {
      vect2(0) = end1zHinge->getTangent();
    }
    if (end2zHinge != 0) {
      vect2(1) = end2zHinge->getTangent();
    }
    return eleInfo.setVector(vect2);

  case 19:
    kb.Zero();
    kb(0,0) = EAoverL;
    kb(5,5) = GJoverL;
    kb(1,1) = kzTrial(0,0);
    kb(2,2) = kzTrial(1,1);
    kb(1,2) = kzTrial(0,1);
    kb(2,1) = kzTrial(1,0);
    kb(3,3) = kyTrial(0,0);
    kb(4,4) = kyTrial(1,1);
    kb(3,4) = kyTrial(0,1);
    kb(4,3) = kyTrial(1,0);    
    return eleInfo.setMatrix(kb);

  default:
    return -1;
  }
}

int
ComponentElement3d::setParameter(const char **argv, int argc, Parameter &param)
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
  if (strcmp(argv[0],"Iz") == 0)
    return param.addObject(3, this);
  if (strcmp(argv[0],"Iy") == 0)
    return param.addObject(4, this);  

  if (strcmp(argv[0],"G") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0],"J") == 0)
    return param.addObject(6, this);  
  
  return -1;
}

int
ComponentElement3d::updateParameter (int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		E = info.theDouble;
		EAoverL = E*A/L;
		EIzoverL2 = 2*E*Iz/L;
		EIzoverL4 = 2*EIzoverL2;
		EIyoverL2 = 2*E*Iy/L;
		EIyoverL4 = 2*EIyoverL2;		
		return 0;
	case 2:
		A = info.theDouble;
		EAoverL = E*A/L;		
		return 0;
	case 3:
		Iz = info.theDouble;
		EIzoverL2 = 2*E*Iz/L;
		EIzoverL4 = 2*EIzoverL2;		
		return 0;
	case 4:
		Iy = info.theDouble;
		EIyoverL2 = 2*E*Iy/L;
		EIyoverL4 = 2*EIyoverL2;				
		return 0;
	case 5:
	  G = info.theDouble;
	  GJoverL = G*J/L;
	case 6:
	  J = info.theDouble;
	  GJoverL = G*J/L;	  
	default:
		return -1;
	}
}

