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
                                                                        
// $Revision: 6149 $
// $Date: 2015-11-18 16:04:51 -0300 (Wed, 18 Nov 2015) $
// $URL: svn://peera.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/CatenaryCable/CatenaryCable.cpp $
                                                                        
                                                                        
// Written: jaabell  (Jose Abell)
// Created: May 2017
// Revision: A
//
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
//
// What: "@(#) CatenaryCable.C, revA"

#include <CatenaryCable.h>
#include <Information.h>
#include <Parameter.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>
#include <Renderer.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <ElementResponse.h>

#include <ElementalLoad.h>

#ifdef _USEQUADMATH
//#include <quadmath.h>
#endif


//Types of mass matrices available 
#define CATENARY_CABLE_MASS_LUMPED 0
#define CATENARY_CABLE_MASS_INTEGRATION 1
#define CATENARY_CABLE_MASS_CLOUGHSTYLE 2
#define CATENARY_CABLE_MASS_EQUVALENTTRUSS 3

// initialise the class wide variables
Matrix CatenaryCable::Flexibility(3,3);
Matrix CatenaryCable::Stiffness(6,6);
Matrix CatenaryCable::Mass(6,6);
Matrix CatenaryCable::ZeroMatrix(6,6);
Vector CatenaryCable::Forces(6);

// constructor:
//  responsible for allocating the necessary space needed by each object
//  and storing the tags of the CatenaryCable end nodes.

#include <elementAPI.h>
#define OPS_Export 
static int num_CatenaryCableElement = 0;

OPS_Export void *
OPS_CatenaryCableElement()
{

  if (num_CatenaryCableElement == 0) {
    num_CatenaryCableElement++;
    opserr<<"CatenaryCableElement element - Written: P. Ibanez and J. A. Abell (UANDES). www.joseabell.com.\n";
  }

  Element *theElement = 0;

  int numRemainingArgs = OPS_GetNumRemainingInputArgs();


  //(int tag, int node1, int node2, double weight, double E, double A, double L0, double alpha, double temperature_change = 100., double rho)
  if (numRemainingArgs < 4) 
  {
    opserr << "Invalid Args want: element CatenaryCable $tag $iNode $jNode $weight $E $A $L0 $alpha $temperature_change $rho $errorTol $Nsubsteps\n";
    return 0; 
  }

  if (numRemainingArgs != 13 )
  {
    opserr << "Got " << numRemainingArgs << " args. Expected 13\n";
        opserr << "Invalid Args want: element CatenaryCable $tag $iNode $jNode $weight $E $A $L0 $alpha $temperature_change $rho $errorTol $Nsubsteps $massType\n";
        return 0; // it's a CatenaryCableSection
  }

  int iData[3];
  double dData[8];
  
  int numData = 3;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer (tag, iNode, jNode) in element CatenaryCable " << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING:  element CatenaryCable - invalid double data. Check $weight $E $A $L0 $alpha $temperature_change $rho $errorTol $Nsubsteps $massType\n";
    return 0; 
  }

  numData = 1;
  int Nsubsteps = 0;
  if (OPS_GetInt(&numData, &Nsubsteps) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer $Nsubsteps in element CatenaryCable " << endln;
    return 0;
  }

  int massType = 0;
  if (OPS_GetInt(&numData, &massType) != 0) {
    opserr << "WARNING element CatenaryCable - invalid integer $massType in element CatenaryCable " << endln;
    return 0;
  }


  // now create the CatenaryCable
  theElement = new CatenaryCable(iData[0], iData[1], iData[2], 
              dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], Nsubsteps, massType);

  
  if (theElement == 0) {
    opserr << "WARNING: out of memory: element CatenaryCable " << iData[0] << 
      " $iNode $jNode ...\n";
  }

  

  return theElement;
}

CatenaryCable::CatenaryCable(int tag, int node1, int node2, double weight_, double E_, double A_, double L0_, double alpha_, double temperature_change_, double rho_, double error_tol_, int Nsubsteps_, int massType_)
 :Element(tag,ELE_TAG_CatenaryCable),
  connectedExternalNodes(2),
  weight(weight_),
  E(E_),
  A(A_),
  L0(L0_),
  alpha(alpha_),
  temperature_change(temperature_change_),
  rho(rho_),
  error_tol(error_tol_),
  Nsubsteps(Nsubsteps_),
  first_step(true),
  massType(massType_)
{
    
    // ensure the connectedExternalNode ID is of correct size & set values
    if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL CatenaryCable::CatenaryCable - " <<  tag << "failed to create an ID of size 2\n";
      exit(-1);
    }

    connectedExternalNodes(0) = node1;
    connectedExternalNodes(1) = node2;        

    // set node pointers to NULL
    for (int i=0; i<2; i++)
      theNodes[i] = 0;

    load = 0;
    load_lastcommit = 0;
    load_incl_inertia = 0;
}

// constructor:
//   invoked by a FEM_ObjectBroker - blank object that recvSelf needs
//   to be invoked upon
CatenaryCable::CatenaryCable()
:Element(0,ELE_TAG_CatenaryCable),     
  connectedExternalNodes(2),
  weight(0),
  E(0),
  A(0),
  L0(0),
  alpha(0),
  temperature_change(0),
  rho(0),
  error_tol(0),
  Nsubsteps(0),
  first_step(true),
  massType(0)
{
    // ensure the connectedExternalNode ID is of correct size 
  if (connectedExternalNodes.Size() != 2) {
      opserr << "FATAL CatenaryCable::CatenaryCable - failed to create an ID of size 2\n";
      exit(-1);
  }

  for (int i=0; i<2; i++)
    theNodes[i] = 0;

  load = 0;
  load_lastcommit = 0;
  load_incl_inertia = 0;

}

//  destructor
//     delete must be invoked on any objects created by the object
//     and on the matertial object.
CatenaryCable::~CatenaryCable()
{

}


int
CatenaryCable::getNumExternalNodes(void) const
{
    return 2;
}

const ID &
CatenaryCable::getExternalNodes(void) 
{
    return connectedExternalNodes;
}

Node **
CatenaryCable::getNodePtrs(void) 
{
  return theNodes;
}

int
CatenaryCable::getNumDOF(void) 
{
    return 6;
}


// method: setDomain()
//    to set a link to the enclosing Domain and to set the node pointers.
//    also determines the number of dof associated
//    with the CatenaryCable element, we set matrix and vector pointers,
//    allocate space for t matrix, determine the length
//    and set the transformation matrix.
void
CatenaryCable::setDomain(Domain *theDomain)
{
    // check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0) {
    	theNodes[0] = 0;
    	theNodes[1] = 0;
    	return;
    }

    // first set the node pointers
    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    theNodes[0] = theDomain->getNode(Nd1);
    theNodes[1] = theDomain->getNode(Nd2);	
    
    // if can't find both - send a warning message
    if ((theNodes[0] == 0) || (theNodes[1] == 0)) {
      if (theNodes[0] == 0)
      	opserr <<"CatenaryCable::setDomain() - CatenaryCable" << this->getTag() << " node " << Nd1 << "does not exist in the model\n";
      else
	      opserr <<"CatenaryCable::setDomain() - CatenaryCable" << this->getTag() << " node " << Nd2 << "does not exist in the model\n";

      return;
    }

    // now determine the number of dof and the dimesnion    
    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();	

    if(L0 <= 0)
    {
      const Vector &end1Crd = theNodes[0]->getCrds();
      const Vector &end2Crd = theNodes[1]->getCrds();
      double dx = end2Crd(0) - end1Crd(0);
      double dy = end2Crd(1) - end1Crd(1);
      double dz = end2Crd(2) - end1Crd(2);
      L0 = sqrt(dx*dx + dy*dy + dz*dz);
    }

    // if differing dof at the ends - print a warning message
    if (dofNd1 != dofNd2) 
    {
      opserr <<"WARNING CatenaryCable::setDomain(): nodes " << Nd1 << " and " << Nd2 <<   "have differing dof at ends for CatenaryCable " << this->getTag() << endln;
      return;
    }	

    // call the base class method
    this->DomainComponent::setDomain(theDomain);



    // create the load vector
    if (load == 0)
    {
      load = new Vector(6);
      load_incl_inertia = new Vector(6);
    }
    if (load_lastcommit == 0)
    {
      load_lastcommit = new Vector(6);
    }

    Flexibility.Zero();
    Stiffness.Zero();
    Mass.Zero();
    ZeroMatrix.Zero();
    Forces.Zero();

    if (load == 0) {
      opserr << "CatenaryCable::setDomain - CatenaryCable " << this->getTag() <<  "out of memory creating vector of size" << 6 << endln;
      exit(-1);
      return;
    }          

    w1 = 0;
    w2 = 0;
    w3 = weight;
    
}   	 


int
CatenaryCable::commitState()
{
  int retVal = 0;

  *load_lastcommit = *load;
  KE_n = KE;
  PE_n = PE;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) 
  {
    opserr << "CatenaryCable::commitState () - failed in base class\n";
  }    

  return retVal;
}

int
CatenaryCable::revertToLastCommit()
{
  return 0;
}

int
CatenaryCable::revertToStart()
{
  return 0;
}

int
CatenaryCable::update(void)
{
  // opserr << "CatenaryCable::update\n";


  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();
  const Vector &end1Disp = theNodes[0]->getTrialDisp();
  const Vector &end2Disp = theNodes[1]->getTrialDisp();


  // double acc1 = 0;
  // double acc2 = 0;
  // double acc3 = 0;
  // int count = 0;
  // for (int i = 0; i < 2; i++) 
  // {
  //   const Vector &accel = theNodes[i]->getTrialAccel();
  //   acc1 += 0.5*accel(0);
  //   acc2 += 0.5*accel(1);
  //   acc3 += 0.5*accel(2);
  // }
  // w1 += -rho*acc1;
  // w2 += -rho*acc2;
  // w3 += -rho*acc3;


  // const Vector &end1Disp_commit = theNodes[0]->getTrialDisp();
  // const Vector &end2Disp_commit = theNodes[1]->getTrialDisp();

  lx0 = end2Crd(0) + end2Disp(0) - (end1Crd(0) + end1Disp(0));
  ly0 = end2Crd(1) + end2Disp(1) - (end1Crd(1) + end1Disp(1));
  lz0 = end2Crd(2) + end2Disp(2) - (end1Crd(2) + end1Disp(2));

  // opserr << "end1Crd = " << end1Crd << endln;
  // opserr << "end1Disp = " << end1Disp << endln;
  // opserr << "end2Crd = " << end2Crd << endln;
  // opserr << "end2Disp = " << end2Disp << endln;

  // opserr <<"   w1 = " <<  w1 << endln; 
  // opserr <<"   w2 = " <<  w2 << endln; 
  // opserr <<"   w3 = " <<  w3 << endln; 


  compute_lambda0();

  //inicializar fuerzas
  double l10 = SQRT(lx0*lx0 + ly0*ly0);
  
  double f10 = (w3*l10)/(2.*lambda0);
  double f20 = 0;
  double f30;

  if (lambda0 >10)
    f30 = -(w3/2.)*(-lz0*(1.) + L0);
  else
    f30 = -(w3/2.)*(-lz0*((cosh(lambda0))/(sinh(lambda0))) + L0);

  // #condicion de objetividad.
  double theta = atan2(ly0 , lx0);
  
  double f10n = cos(theta)*f10 - sin(theta)*f20;
  double f20n = sin(theta)*f10 + cos(theta)*f20;
  double f30n = f30;
  
  f1 = f10n;
  f2 = f20n;
  f3 = f30n;

  // opserr <<"   lambda0 = " <<  lambda0 << endln; 
  // opserr <<"   f10n = " <<  f10n << endln; 
  // opserr <<"   f20n = " <<  f20n << endln; 
  // opserr <<"   f30n = " <<  f30n << endln; 


  compute_projected_lengths();


  // opserr <<"   f1 = " <<  f1 << endln; 
  // opserr <<"   f2 = " <<  f2 << endln; 
  // opserr <<"   f3 = " <<  f3 << endln; 

  // opserr <<"   l1 = " <<  l1 << endln; 
  // opserr <<"   l2 = " <<  l2 << endln; 
  // opserr <<"   l3 = " <<  l3 << endln; 

  //Misclosure vector   dl = sp.matrix([[lx0-lxi],[ly0-lyi],[lz0-lzi]], dtype = sp.double())
  static Vector dl(3);
  dl.Zero();

  dl(0) = lx0 - l1;
  dl(1) = ly0 - l2;
  dl(2) = lz0 - l3;
  
  // const double error_tol = 1e-6;

  static Vector Fi0(3);
  Fi0(0) = f1;
  Fi0(1) = f2;
  Fi0(2) = f3;

  double dl_max = dl.pNorm(-1);   // Computes the infinity norm
  double relative_error = fabs(dl_max)/L0;
  double max_relative_error = 0;
  int iter_max = 0;
  double min_relative_error = 1/error_tol;
  int iter_min = 0;
  // opserr << " update: relative_error = " << relative_error << endln;

  int iter = 0;
  while( relative_error > error_tol)
  {
    if(relative_error < min_relative_error)
    {
      min_relative_error = relative_error;
      iter_min = iter;
    }
    if(relative_error > max_relative_error)
    {
      max_relative_error = relative_error;
      iter_max = iter;
    }
    // opserr << " update: relative_error = " << relative_error << endln;
    f1 = Fi0(0);
    f2 = Fi0(1);
    f3 = Fi0(2);
    // f4 = -f1 - w1*L0;
    // f5 = -f2 - w2*L0;
    // f6 = -f3 - w3*L0;

    compute_projected_lengths();

    //Update misclosure
    dl(0) = lx0 - l1;
    dl(1) = ly0 - l2;
    dl(2) = lz0 - l3;

    // const int Nsubstep = 1;
    for(int substep = 0; substep < Nsubsteps; substep++)
    {
        f1 = Fi0(0);
        f2 = Fi0(1);
        f3 = Fi0(2);
        // f4 = -f1
        // f5 = -f2
        // f6 = -f3 - w3*L0

        compute_flexibility_matrix();
    
        // static Matrix K(3,3);
        static Vector dF(3);
        
        Flexibility.Solve(dl, dF);
    
        dF  = dF / Nsubsteps;
          
        Fi0 = Fi0 + dF;
    }

    dl_max = dl.pNorm(-1);   // Computes the infinity norm
    relative_error = fabs(dl_max)/L0;

    iter+= 1;

    if(iter > 100)
    {
      opserr << "CatenaryCable::update() - Failed to converge.\n";
      opserr << "   tag = " << this->getTag() << endln;
      opserr << "   L0 = " << L0 << endln;
      opserr << "   relative_error = " << relative_error << endln;
      opserr << "   iteratations = " << iter << endln;
      opserr << "   min_relative_error = " << min_relative_error << " at iter = " << iter_min << endln;
      opserr << "   max_relative_error = " << max_relative_error << " at iter = " << iter_max << endln;
      opserr << "   Nsubsteps = " << Nsubsteps << endln;
      opserr << "   end1Crd = " << end1Crd << endln;
      opserr << "   end1Disp = " << end1Disp << endln;
      opserr << "   end2Crd = " << end2Crd << endln;
      opserr << "   end2Disp = " << end2Disp << endln;
      opserr << "    w1 = " <<  w1 << endln; 
      opserr << "    w2 = " <<  w2 << endln; 
      opserr << "    w3 = " <<  w3 << endln; 
      opserr << "    lambda0 = " <<  lambda0 << endln; 
      opserr << "    f10n = " <<  f10n << endln; 
      opserr << "    f20n = " <<  f20n << endln; 
      opserr << "    f30n = " <<  f30n << endln; 
      opserr << "    f1 = " <<  f1 << endln; 
      opserr << "    f2 = " <<  f2 << endln; 
      opserr << "    f3 = " <<  f3 << endln; 
      opserr << "    l1 = " <<  l1 << endln; 
      opserr << "    l2 = " <<  l2 << endln; 
      opserr << "    l3 = " <<  l3 << endln; 
      return -1;
    }
  }

  // opserr << "niter = " << iter << endln;

  return 0;
}


const Matrix &
CatenaryCable::getTangentStiff(void)
{
  static Matrix K(3,3);
  K.Zero();
  Stiffness.Zero();
  
  compute_flexibility_matrix();

  Flexibility.Invert(K);

  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
    {
      Stiffness(i,j) = -K(i,j);
      Stiffness(i+3,j+3) = -K(i,j);
      Stiffness(i,j+3) = K(i,j);
      Stiffness(i+3,j) = K(i,j);
    }
  }


  // opserr << "Stiffness = " << Stiffness << endln;

  return Stiffness;
}


const Matrix &
CatenaryCable::getInitialStiff(void)
{
    return Stiffness;
}

const Matrix &
CatenaryCable::getDamp(void)
{
  return ZeroMatrix;
}


const Matrix &
CatenaryCable::getMass(void)
{
  computeMass();
  return Mass;
}

void 
CatenaryCable::zeroLoad(void)
{
  load->Zero();
}

int 
CatenaryCable::addLoad(ElementalLoad *theLoad, double loadFactor)

{  
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  if (type == LOAD_TAG_Beam3dUniformLoad) 
  {
      // opserr <<"CatenaryCable::addLoad - Uniform Load - loadFactor = " << loadFactor << endln; 
      w1 = loadFactor*data(0);
      w2 = loadFactor*data(1);
      w3 = loadFactor*data(2);
      // opserr <<  "   w = ("  <<  w1 <<  ", " <<  w2 << ", " <<  w3 << ")" << endln; 
      
      return 0;
  }

  opserr <<"CatenaryCable::addLoad - load type (" << type <<") unknown for CatenaryCable with tag: " << this->getTag() << endln; 
  return -1;
}

int 
CatenaryCable::addInertiaLoadToUnbalance(const Vector &accel)
{
  opserr << "CatenaryCable::addInertiaLoadToUnbalance\n";
  /* 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  IS NOT GETTING CALLED

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  */
  // static Vector resid(6);
  // resid.Zero();

  getResistingForce();

  // // check for a quick return
  // if ( rho == 0.0) 
  //   return 0;

  // int count = 0;
  // for (int i = 0; i < 2; i++) 
  // {
  //   const Vector &Raccel = theNodes[i]->getRV(accel);
  //   for (int j = 0; j < 3; j++)
  //     resid(count++) = Raccel(j);
  // }
  


  // opserr << "accel = " << accel;
  // opserr << "load = " << accel;

  // // want to add ( - fact * M R * accel ) to unbalance
  // computeMass();
  // load->addMatrixVector(1.0, Mass, resid, -1.0);
  
  return 0;
}



const Vector &
CatenaryCable::getResistingForce()
{	


  double f4 = -f1 - w1*L0;
  double f5 = -f2 - w2*L0;
  double f6 = -f3 - w3*L0;
    
  (*load)(0) = f1;
  (*load)(1) = f2;
  (*load)(2) = f3;
  (*load)(3) = f4;
  (*load)(4) = f5;
  (*load)(5) = f6;

  static Vector disp(6);
  const Vector &end1Disp = theNodes[0]->getIncrDisp();
  const Vector &end2Disp = theNodes[1]->getIncrDisp();
  disp(0) = end1Disp(0);
  disp(1) = end1Disp(1);
  disp(2) = end1Disp(2);
  disp(3) = end2Disp(0);
  disp(4) = end2Disp(1);
  disp(5) = end2Disp(2);

  PE = PE_n + 0.5*((*load_lastcommit + *load)^(disp));

  return *load;
}


const Vector &
CatenaryCable::getResistingForceIncInertia()
{	
  
  // // subtract external load
  // (*theVector) -= *theLoad;
  
  // // now include the mass portion
  // if (L != 0.0 && rho != 0.0) {
    
  //   // add inertia forces from element mass
  //   const Vector &accel1 = theNodes[0]->getTrialAccel();
  //   const Vector &accel2 = theNodes[1]->getTrialAccel();	
    
  //   int numDOF2 = numDOF/2;
    
  //   if (cMass == 0)  {
  //     // lumped mass matrix
  //     double m = 0.5*rho*L;
  //     for (int i = 0; i < dimension; i++) {
  //       (*theVector)(i) += m*accel1(i);
  //       (*theVector)(i+numDOF2) += m*accel2(i);
  //     }
  //   } else  {
  //     // consistent mass matrix
  //     double m = rho*L/6.0;
  //     for (int i=0; i<dimension; i++) {
  //       (*theVector)(i) += 2.0*m*accel1(i) + m*accel2(i);
  //       (*theVector)(i+numDOF2) += m*accel1(i) + 2.0*m*accel2(i);
  //     }
  //   }
    
  //   // add the damping forces if rayleigh damping
  //   if (doRayleighDamping == 1 && (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
  //     theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
  // } else {
    
  //   // add the damping forces if rayleigh damping
  //   if (doRayleighDamping == 1 && (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0))
  //     theVector->addVector(1.0, this->getRayleighDampingForces(), 1.0);
  // }



  this->getResistingForce();



  computeMass();



  static Vector accel(6);
  static Vector veloc(6);
  accel.Zero();
  veloc.Zero();


  // check for a quick return
  if ( rho == 0.0) 
    return *load;

  int count = 0;
  for (int i = 0; i < 2; i++) 
  {
    const Vector &Raccel = theNodes[i]->getTrialAccel();
    const Vector &Rveloc = theNodes[i]->getTrialVel();
    for (int j = 0; j < 3; j++)
    {      
      accel(count) = Raccel(j);
      veloc(count) = Rveloc(j);
      count++;
    }
  }

  *load_incl_inertia = *load;
  load_incl_inertia->addMatrixVector(1.0, Mass, accel, 1.0);

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      *load_incl_inertia += this->getRayleighDampingForces();


  KE = 0;

  for(int i = 0; i < 6; i++)
  {    
    for(int j = 0; j < 6; j++)
    {
      KE += veloc(i)*Mass(i,j)*veloc(j)*0.5;
    }
  }
  return *load_incl_inertia;


  
}

int
CatenaryCable::sendSelf(int commitTag, Channel &theChannel)
{
  // int res;

  // // note: we don't check for dataTag == 0 for Element
  // // objects as that is taken care of in a commit by the Domain
  // // object - don't want to have to do the check if sending data
  // int dataTag = this->getDbTag();

  // // CatenaryCable packs it's data into a Vector and sends this to theChannel
  // // along with it's dbTag and the commitTag passed in the arguments

  // static Vector data(12);
  // data(0) = this->getTag();
  // data(1) = dimension;
  // data(2) = numDOF;
  // data(3) = A;
  // data(6) = rho;
  // data(7) = doRayleighDamping;
  // data(8) = cMass;
  
  // data(4) = theMaterial->getClassTag();
  // int matDbTag = theMaterial->getDbTag();
  
  // if (initialDisp != 0) {
  //   for (int i=0; i<dimension; i++) {
  //     data[9+i] = initialDisp[i];
  //   }
  // }

  // // NOTE: we do have to ensure that the material has a database
  // // tag if we are sending to a database channel.
  // if (matDbTag == 0) {
  //   matDbTag = theChannel.getDbTag();
  //   if (matDbTag != 0)
  //     theMaterial->setDbTag(matDbTag);
  // }
  // data(5) = matDbTag;

  // res = theChannel.sendVector(dataTag, commitTag, data);
  // if (res < 0) {
  //   opserr <<"WARNING CatenaryCable::sendSelf() - " << this->getTag() << " failed to send Vector\n";
  //   return -1;
  // }	      

  // // CatenaryCable then sends the tags of it's two end nodes
  // res = theChannel.sendID(dataTag, commitTag, connectedExternalNodes);
  // if (res < 0) {
  //   opserr <<"WARNING CatenaryCable::sendSelf() - " << this->getTag() << " failed to send Vector\n";
  //   return -2;
  // }

  // // finally CatenaryCable asks it's material object to send itself
  // res = theMaterial->sendSelf(commitTag, theChannel);
  // if (res < 0) {
  //   opserr <<"WARNING CatenaryCable::sendSelf() - " << this->getTag() << " failed to send its Material\n";
  //   return -3;
  // }

  return 0;
}

int
CatenaryCable::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
 //  int res;
 //  int dataTag = this->getDbTag();

 //  // CatenaryCable creates a Vector, receives the Vector and then sets the 
 //  // internal data with the data in the Vector

 //  static Vector data(12);
 //  res = theChannel.recvVector(dataTag, commitTag, data);
 //  if (res < 0) {
 //    opserr <<"WARNING CatenaryCable::recvSelf() - failed to receive Vector\n";
 //    return -1;
 //  }	      

 //  this->setTag((int)data(0));
 //  dimension = (int)data(1);
 //  numDOF = (int)data(2);
 //  A = data(3);
 //  rho = data(6);
 //  doRayleighDamping = (int)data(7);
 //  cMass = (int)data(8);

 //  initialDisp = new double[dimension];
 //  for (int i=0; i<dimension; i++)
 //    initialDisp[i] = 0.0;

 //  int initial = 0;
 //  for (int i=0; i<dimension; i++) {
 //    if (data(9+i) != 0.0) {
 //      initial = 1;
 //    }
 //  }
  
 //  if (initial != 0) {
 //    for (int i=0; i<dimension; i++) {
 //      initialDisp[i] = data(9+i);
 //    }    
 //  }
  
 //  // CatenaryCable now receives the tags of it's two external nodes
 //  res = theChannel.recvID(dataTag, commitTag, connectedExternalNodes);
 //  if (res < 0) {
 //    opserr <<"WARNING CatenaryCable::recvSelf() - " << this->getTag() << " failed to receive ID\n";
 //    return -2;
 //  }

 //  // finally CatenaryCable creates a material object of the correct type,
 //  // sets its database tag and asks this new object to recveive itself.

 //  int matClass = (int)data(4);
 //  int matDb = (int)data(5);

 //  // check if we have a material object already & if we do if of right type
 //  if ((theMaterial == 0) || (theMaterial->getClassTag() != matClass)) {

 //    // if old one .. delete it
 //    if (theMaterial != 0)
 //      delete theMaterial;

 //    // create a new material object
 //    theMaterial = theBroker.getNewUniaxialMaterial(matClass);
 //    if (theMaterial == 0) {
 //      opserr <<"WARNING CatenaryCable::recvSelf() - " << this->getTag() 
	// << " failed to get a blank Material of type " << matClass << endln;
 //      return -3;
 //    }
 //  }

 //  theMaterial->setDbTag(matDb); // note: we set the dbTag before we receive the material
 //  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
 //  if (res < 0) {
 //    opserr <<"WARNING CatenaryCable::recvSelf() - "<< this->getTag() << "failed to receive its Material\n";
 //    return -3;    
 //  }

  return 0;
}

int
CatenaryCable::displaySelf(Renderer &theViewer, int displayMode, float fact, 
		   const char **displayModes, int numModes)
{  
  return -1;
}



void
CatenaryCable::Print(OPS_Stream &s, int flag)
{
    /* compute the strain and axial force in the member
    double strain, force;
    strain = theMaterial->getStrain();
    force = A * theMaterial->getStress();

    if (flag == OPS_PRINT_CURRENTSTATE) { // print everything
        s << "Element: " << this->getTag();
        s << " type: CatenaryCable  iNode: " << connectedExternalNodes(0);
        s << " jNode: " << connectedExternalNodes(1);
        s << " Area: " << A << " Mass/Length: " << rho;
        s << " cMass: " << cMass;

        s << " \n\t strain: " << strain;
        if (initialDisp != 0) {
            s << " initialDisplacements: ";
            for (int i = 0; i < dimension; i++)
                s << initialDisp[i] << " ";
        }

        s << " axial load: " << force;

        if (L != 0.0) {
            int numDOF2 = numDOF / 2;
            double temp;
            for (int i = 0; i < dimension; i++) {
                temp = cosX[i] * force;
                (*theVector)(i) = -temp;
                (*theVector)(i + numDOF2) = temp;
            }
            s << " \n\t unbalanced load: " << *theVector;
        }

        s << " \t Material: " << *theMaterial;
        s << endln;
    }
    
    if (flag == 1) {
        s << this->getTag() << "  " << strain << "  ";
        s << force << endln;
    }*/
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"CatenaryCable\", ";
        s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"E\": " << E << ", ";
        s << "\"A\": " << A << ", ";
        s << "\"L0\": " << L0 << ", ";
        s << "\"alpha\": " << alpha << ", ";
        s << "\"deltaT\": " << temperature_change << ", ";
        s << "\"massperlength\": " << rho << "\"}";
    }
}


Response*
CatenaryCable::setResponse(const char **argv, int argc, OPS_Stream &output)
{

    Response *theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType","CatenaryCable");
    output.attr("eleTag",this->getTag());
    output.attr("node1",connectedExternalNodes[0]);
    output.attr("node2",connectedExternalNodes[1]);

    // //
    // // we compare argv[0] for known response types for the CatenaryCable
    // //


    if ((strcmp(argv[0],"force") == 0) || (strcmp(argv[0],"forces") == 0) 
        || (strcmp(argv[0],"globalForce") == 0) || (strcmp(argv[0],"globalForces") == 0))
    {
            output.tag("ResponseType", "f1");
            output.tag("ResponseType", "f2");
            output.tag("ResponseType", "f3");
            output.tag("ResponseType", "f4");
            output.tag("ResponseType", "f5");
            output.tag("ResponseType", "f6");
            theResponse =  new ElementResponse(this, 1, Vector(6));

    } 
    else if (strcmp(argv[0],"energy") == 0)
    {
            output.tag("ResponseType", "KineticEnergy");
            output.tag("ResponseType", "PotentialEnergy");
            theResponse =  new ElementResponse(this, 2, Vector(2));
    }

    return theResponse;
}

int 
CatenaryCable::getResponse(int responseID, Information &eleInfo)
{
    // double strain;

    switch (responseID) {
    case 1:
        return eleInfo.setVector(this->getResistingForce());

    case 2:
        return eleInfo.setVector(this->getEnergyVector());

    // case 3:
    //     if (L == 0.0) {
    //         strain = 0.0;
    //     } else {
    //         strain = theMaterial->getStrain();
    //     }
    //     return eleInfo.setDouble(L * strain);

    default:
        return 0;
    }
  return 0;
}


Vector CatenaryCable::getEnergyVector()
{
  Vector energy(2);
  energy(0) = KE_n;
  energy(1) = PE_n;

  return energy;
}



void CatenaryCable::compute_lambda0(void)
{

    lambda0 = 0;
    if (lx0*lx0 + ly0*ly0 == 0)
      lambda0 = 1e6;
    else if( L0*L0 <= (lx0*lx0 + ly0*ly0 + lz0*lz0))
      lambda0 = 0.2;
    else if( L0*L0 > lx0*lx0 + ly0*ly0 + lz0*lz0)
      lambda0 = SQRT(3*((L0*L0 - lz0*lz0) / (lx0*lx0 + ly0*ly0)) - 1);

    // opserr << "L0 = " << L0 << endln;
    // opserr << "lx0 = " << lx0 << endln;
    // opserr << "ly0 = " << ly0 << endln;
    // opserr << "lz0 = " << lz0 << endln;
}


void CatenaryCable::compute_projected_lengths(void)
{
    //Correct weights for acceleration
    // double acc1 = 0;
    // double acc2 = 0;
    // double acc3 = 0;
    // int count = 0;
    // for (int i = 0; i < 2; i++) 
    // {
    //   const Vector &accel = theNodes[i]->getTrialAccel();
    //   acc1 += 0.5*accel(0);
    //   acc2 += 0.5*accel(1);
    //   acc3 += 0.5*accel(2);
    // }
    FLOATTYPE w_1 = w1;// - rho*acc1;
    FLOATTYPE w_2 = w2;// - rho*acc2;
    FLOATTYPE w_3 = w3;// - rho*acc3;

    FLOATTYPE w = SQRT(w_1*w_1 + w_2*w_2 + w_3*w_3);
    FLOATTYPE a1 = (f1*w_1) + (f2*w_2) + (f3*w_3);
    FLOATTYPE t1 = SQRT((f1*f1 ) + (f2*f2) + (f3*f3));
    FLOATTYPE a = (-(w_1*L0) - f1);
    FLOATTYPE b = (-(w_2*L0) - f2);
    FLOATTYPE c = (-(w_3*L0) - f3);
    FLOATTYPE t2 = SQRT(a*a + b*b + c*c);
    // l1 = (-(L0*f1)/(E*A))-(((L0*L0)*w_1)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_1*(t1-t2))+(((w*w)*f1)-(a1*w_1))*(log((a1/w)+t1)-log((L0*w)+(a1/w)+t2)));
    // l2 = (-(L0*f2)/(E*A))-(((L0*L0)*w_2)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_2*(t1-t2))+(((w*w)*f2)-(a1*w_2))*(log((a1/w)+t1)-log((L0*w)+(a1/w)+t2)));
    // l3 = (-(L0*f3)/(E*A))-(((L0*L0)*w_3)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_3*(t1-t2))+(((w*w)*f3)-(a1*w_3))*(log((a1/w)+t1)-log((L0*w)+(a1/w)+t2)));

    l1 = (-(L0*f1)/(E*A))-(((L0*L0)*w_1)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_1*(t1-t2))+(((w*w)*f1)-(a1*w_1))*LOG(((a1/w)+t1)/((L0*w)+(a1/w)+t2)));
    l2 = (-(L0*f2)/(E*A))-(((L0*L0)*w_2)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_2*(t1-t2))+(((w*w)*f2)-(a1*w_2))*LOG(((a1/w)+t1)/((L0*w)+(a1/w)+t2)));
    l3 = (-(L0*f3)/(E*A))-(((L0*L0)*w_3)/(2.*E*A))+((1.+alpha*temperature_change)/(w*w*w))*((w*w_3*(t1-t2))+(((w*w)*f3)-(a1*w_3))*LOG(((a1/w)+t1)/((L0*w)+(a1/w)+t2)));

}


void CatenaryCable::compute_flexibility_matrix(void)
{
    //Correct weights for acceleration
    // double acc1 = 0;
    // double acc2 = 0;
    // double acc3 = 0;
    // int count = 0;
    // for (int i = 0; i < 2; i++) 
    // {
    //   const Vector &accel = theNodes[i]->getTrialAccel();
    //   acc1 += 0.5*accel(0);
    //   acc2 += 0.5*accel(1);
    //   acc3 += 0.5*accel(2);
    // }
    FLOATTYPE w_1 = w1;// - rho*acc1;
    FLOATTYPE w_2 = w2;// - rho*acc2;
    FLOATTYPE w_3 = w3;// - rho*acc3;

    FLOATTYPE w = sqrt(w_1*w_1 + w_2*w_2 + w_3*w_3);
    FLOATTYPE a1 = f1*w_1 + f2*w_2 + f3*w_3;
    FLOATTYPE t1 = sqrt(f1*f1+f2*f2+f3*f3);
    FLOATTYPE a = (-(w_1*L0)-f1);
    FLOATTYPE b = (-(w_2*L0)-f2);
    FLOATTYPE c = (-(w_3*L0)-f3);
    FLOATTYPE t2 = sqrt(a*a + b*b + c*c);

    FLOATTYPE W[3] = {w_1, w_2, w_3};
    FLOATTYPE f4 = -f1 - w_1*L0;
    FLOATTYPE f5 = -f2 - w_2*L0;
    FLOATTYPE f6 = -f3 - w_3*L0;
    
    FLOATTYPE F[6] = {f1, f2, f3, f4, f5, f6};

    FLOATTYPE b0 = 0.;
    FLOATTYPE b1 = 0.;
    FLOATTYPE b2 = 0.;

    
    for(int i = 0; i < 3; i++)
    {
      for(int j = 0; j < 3; j++)
      {
        b1 = (-w * W[i] * ( F[j+3]/t2 + F[j]/t1 ) + 
          ( w*w * F[i] - a1 * W[i]) * (
            (w*F[j] + W[j]*(L0*w + t2)) / (t2 * (L0*w*w + a1 + w*t2))
            - 
            (w*F[j] + W[j]*t1) / (t1 * (a1 + w*t1))
            )
          );
        if (i==j)
        {
                  b0 = -L0 / (E*A);
                  b2 = W[i]*W[i] - w*w;
        }
        else
        {
          b0 = 0.  ;
          b2 = W[i]*W[j];
        }  
        // Flexibility(i, j) = b0 - (1 + alpha*temperature_change)/(w*w*w)*(b1 + b2*(log( a1/w + t1) - log(a1/w + t2 + L0*w)));
        Flexibility(i, j) = b0 - (1 + alpha*temperature_change)/(w*w*w)*(b1 + b2*LOG(( a1/w + t1)/(a1/w + t2 + L0*w)));
      }
    }

}


void CatenaryCable::computeMass()
{
  switch (massType)
  {
    case CATENARY_CABLE_MASS_LUMPED:
      computeMassLumped();
      break;
    case CATENARY_CABLE_MASS_CLOUGHSTYLE:
      computeMassCloughStyle();
      break;
    case CATENARY_CABLE_MASS_INTEGRATION:
      computeMassByIntegration();
      break;    
    case CATENARY_CABLE_MASS_EQUVALENTTRUSS:
      computeMassEquivalentTruss();
      break;
    default:
      opserr << "CatenaryCable::computeMass() -- Unknown massType = " << massType << endln;
  }
}

void CatenaryCable::computeMassLumped()
{
  double nodal_mass = rho*L0/2;
  Mass(0,0) = nodal_mass;
  Mass(1,1) = nodal_mass;
  Mass(2,2) = nodal_mass;
  Mass(3,3) = nodal_mass;
  Mass(4,4) = nodal_mass;
  Mass(5,5) = nodal_mass;
}

void CatenaryCable::computeMassEquivalentTruss()
{

  double m = rho*L0/6.0;
  for (int i = 0; i < 3; i++) {
    Mass(i,i) = 2.0*m;
    Mass(i,i+3) = m;
    Mass(i+3,i) = m;
    Mass(i+3,i+3) = 2.0*m;
  }
}

void CatenaryCable::computeMassByIntegration()
{
    opserr << "CatenaryCable::computeMass() -- Mass by integration not yet available -- Defaulting to lumped " << endln;
    computeMassLumped();  // To be implemented
}

void CatenaryCable::computeMassCloughStyle()
{
  double total_mass = rho*L0;
  double f1x = fabs((*load)(0));
  double f2x = fabs((*load)(3));
  double f1y = fabs((*load)(1));
  double f2y = fabs((*load)(4));
  double f1z = fabs((*load)(2));
  double f2z = fabs((*load)(5));
  double f1 = sqrt(f1x*f1x + f1y*f1y + f1z*f1z);
  double f2 = sqrt(f2x*f2x + f2y*f2y + f2z*f2z);
  double m1 = total_mass*f1/(f1+f2);
  double m2 = total_mass*f1/(f1+f2);
  Mass(0,0) = m1;
  Mass(1,1) = m1;
  Mass(2,2) = m1;
  Mass(3,3) = m2;
  Mass(4,4) = m2;
  Mass(5,5) = m2;
}

